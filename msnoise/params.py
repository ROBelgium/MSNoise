"""
MSNoise Layered Parameters
==========================

:class:`LayeredParams` replaces the flat ``AttribDict`` previously returned by
:func:`~msnoise.api.get_merged_params_for_lineage`.  Each step's configuration
lives in its own namespace, so keys that share a name across categories (e.g.
``freqmin`` in ``filter`` *and* ``mwcs``) never collide.

Usage::

    # Category access
    params.cc.maxlag           # → float
    params.filter.freqmin      # → float  (different from params.mwcs.freqmin)
    params.mwcs_dtt.dtt_minlag # → float

    # Cross-step convenience properties (most compute steps need these)
    params.output_folder                         # from global
    params.mov_stack                             # from stack
    params.components_to_compute                 # from cc
    params.components_to_compute_single_station  # from cc

    # Dynamic innermost-layer access (used in aggregate_dvv_pairs)
    params.category          # → "mwcs_dtt_dvv"
    params.category_layer    # → AttribDict for that category

    # Serialisation
    params.to_yaml("run_params.yaml")
    p2 = LayeredParams.from_yaml("run_params.yaml")
"""

from __future__ import annotations

import datetime
from collections import OrderedDict
from obspy.core.util.attribdict import AttribDict


class LayeredParams:
    """Namespaced parameter container for an MSNoise lineage.

    Each step in the lineage contributes one layer keyed by its category
    (e.g. ``"global"``, ``"preprocess"``, ``"cc"``, ``"filter"`` …).
    Attributes on ``LayeredParams`` itself delegate to the appropriate layer;
    there is no flat merged view and no silent key collision.
    """

    # Keys whose category is fixed regardless of which step is running.
    # These are surfaced as top-level properties for convenience.
    _FIXED_PROPERTIES = {
        "output_folder": "global",
        "mov_stack": "stack",
        "components_to_compute": "cc",
        "components_to_compute_single_station": "cc",
    }

    # ------------------------------------------------------------------ #
    # Construction                                                         #
    # ------------------------------------------------------------------ #

    def __init__(self) -> None:
        # Ordered from global (outermost) to current step (innermost).
        object.__setattr__(self, "_layers", OrderedDict())
        object.__setattr__(self, "_lineage_names", [])

    def _add_layer(self, category: str, attrib_dict) -> None:
        """Add one step's config as a named layer.

        :param category: Config category name e.g. ``"mwcs_dtt"``.
        :param attrib_dict: ``AttribDict`` (or any mapping) of key→value for
            this category's config set.
        """
        self._layers[category] = attrib_dict

    def _set_lineage_names(self, names: list[str]) -> None:
        object.__setattr__(self, "_lineage_names", list(names))

    # ------------------------------------------------------------------ #
    # Category access                                                      #
    # ------------------------------------------------------------------ #

    def __getattr__(self, name: str):
        """``params.<category>`` returns the AttribDict for that layer.

        ``params.global_`` maps to the ``"global"`` layer (avoids the Python
        keyword).
        """
        layers = object.__getattribute__(self, "_layers")
        # global_ → "global"
        lookup = "global" if name == "global_" else name
        if lookup in layers:
            return layers[lookup]
        raise AttributeError(
            f"LayeredParams has no category {name!r}. "
            f"Available categories: {list(layers)}"
        )

    def __setattr__(self, name: str, value) -> None:
        raise AttributeError(
            "LayeredParams is immutable. Use _add_layer() during construction."
        )

    def __getitem__(self, category: str):
        """``params["mwcs_dtt"]`` — bracket-style category access."""
        layers = object.__getattribute__(self, "_layers")
        if category in layers:
            return layers[category]
        raise KeyError(
            f"LayeredParams has no category {category!r}. "
            f"Available: {list(layers)}"
        )

    # ------------------------------------------------------------------ #
    # Convenience properties                                               #
    # ------------------------------------------------------------------ #

    @property
    def lineage_names(self) -> list[str]:
        """Ordered step-name list e.g. ``['preprocess_1', 'cc_1', ...]``."""
        return object.__getattribute__(self, "_lineage_names")

    @property
    def category(self) -> str:
        """Category name of the innermost (current) step."""
        layers = object.__getattribute__(self, "_layers")
        if not layers:
            raise RuntimeError("LayeredParams has no layers.")
        return next(reversed(layers))

    @property
    def category_layer(self):
        """AttribDict for the innermost step (``params[params.category]``).

        Used in ``aggregate_dvv_pairs`` where the dvv category is determined
        at runtime.
        """
        layers = object.__getattribute__(self, "_layers")
        if not layers:
            raise RuntimeError("LayeredParams has no layers.")
        return layers[next(reversed(layers))]

    @property
    def step_name(self) -> str:
        """Step name of the innermost step e.g. ``'mwcs_dtt_1'``."""
        names = object.__getattribute__(self, "_lineage_names")
        return names[-1] if names else ""

    @property
    def categories(self) -> list[str]:
        """Ordered list of all category names in this params object."""
        return list(object.__getattribute__(self, "_layers"))

    @property
    def output_folder(self) -> str:
        """Output folder path — always read from the ``global`` layer."""
        return self._layers["global"].output_folder

    @property
    def mov_stack(self):
        """Moving-stack list — always read from the ``stack`` layer."""
        return self._layers["stack"].mov_stack

    @property
    def components_to_compute(self) -> list[str]:
        """Components between stations — always read from the ``cc`` layer."""
        return self._layers["cc"].components_to_compute

    @property
    def components_to_compute_single_station(self) -> list[str]:
        """Single-station components — always read from the ``cc`` layer."""
        return self._layers["cc"].components_to_compute_single_station

    # ------------------------------------------------------------------ #
    # Inspection / serialisation helpers                                   #
    # ------------------------------------------------------------------ #

    def as_flat_dict(self) -> dict:
        """Return a flat merged dict (innermost layer wins on collision).

        Use only for legacy compatibility or debugging.  Prefer explicit
        ``params.<category>.<key>`` access in all production code.
        """
        merged = {}
        for layer in object.__getattribute__(self, "_layers").values():
            merged.update(dict(layer))
        return merged

    def to_yaml_string(self) -> str:
        """Serialise the full lineage config chain to a YAML string."""
        import yaml

        layers = object.__getattribute__(self, "_layers")
        names = object.__getattribute__(self, "_lineage_names")

        doc = {
            "msnoise_params_version": 1,
            "generated": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds") + "Z",
            "lineage": "/".join(names),
        }
        for cat, layer in layers.items():
            doc[cat] = _attribdict_to_plain(layer)

        return yaml.dump(doc, default_flow_style=False, sort_keys=False,
                         allow_unicode=True)

    def to_yaml(self, path: str) -> None:
        """Write the full lineage config chain to a YAML file at *path*."""
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(self.to_yaml_string())

    @classmethod
    def from_yaml(cls, path: str) -> "LayeredParams":
        """Reconstruct a :class:`LayeredParams` from a YAML file.

        Does *not* require a database connection — useful for offline
        reproducibility checks.
        """
        import yaml
        with open(path, encoding="utf-8") as fh:
            doc = yaml.safe_load(fh)

        p = cls()
        lineage = doc.get("lineage", "")
        p._set_lineage_names(lineage.split("/") if lineage else [])

        skip = {"msnoise_params_version", "generated", "lineage"}
        for key, value in doc.items():
            if key in skip:
                continue
            if isinstance(value, dict):
                p._add_layer(key, AttribDict(value))

        return p

    def __repr__(self) -> str:
        cats = list(object.__getattribute__(self, "_layers"))
        names = object.__getattribute__(self, "_lineage_names")
        return (
            f"LayeredParams("
            f"lineage={'/'.join(names)!r}, "
            f"categories={cats})"
        )


# ------------------------------------------------------------------ #
# Internal helpers                                                     #
# ------------------------------------------------------------------ #

def _attribdict_to_plain(ad) -> dict:
    """Convert an AttribDict to a plain Python dict for YAML serialisation."""
    result = {}
    for k, v in ad.items():
        if isinstance(v, (list, tuple)):
            result[k] = [list(i) if isinstance(i, tuple) else i for i in v]
        elif hasattr(v, "items"):
            result[k] = _attribdict_to_plain(v)
        else:
            result[k] = v
    return result


def _build_layered_params(
    global_attrib,
    lineage_steps: list,
    lineage_names: list[str],
    step_configs: list,       # list of AttribDict, one per step in lineage_steps
) -> LayeredParams:
    """Build a :class:`LayeredParams` from a lineage.

    :param global_attrib: ``AttribDict`` from :func:`~msnoise.api.get_params`.
    :param lineage_steps: List of WorkflowStep ORM objects (upstream→downstream).
    :param lineage_names: Corresponding step-name strings.
    :param step_configs: Per-step ``AttribDict`` from
        :func:`~msnoise.api.get_config_set_details`, same order as
        *lineage_steps*.
    :returns: :class:`LayeredParams` with one layer per category.
    """
    p = LayeredParams()
    p._set_lineage_names(lineage_names)

    # Global layer first
    p._add_layer("global", global_attrib)

    # One layer per step in lineage order
    for step, cfg in zip(lineage_steps, step_configs):
        cat = step.category
        if cat == "global":
            continue  # already added

        # Post-process cc layer: split comma-separated component strings
        if cat == "cc":
            for key in ("components_to_compute",
                        "components_to_compute_single_station"):
                val = cfg.get(key, None)
                if val is not None and isinstance(val, str):
                    cfg[key] = [v for v in val.split(",") if v]

        # Post-process stack layer: normalise mov_stack to list-of-tuples
        if cat == "stack" and hasattr(cfg, "mov_stack"):
            ms = cfg.mov_stack
            if ms is not None:
                if not isinstance(ms, (list, tuple)):
                    ms = [ms]
                if ms and not isinstance(ms[0], (list, tuple)):
                    ms = [ms]
                cfg["mov_stack"] = [tuple(m) for m in ms]

        p._add_layer(cat, cfg)

    return p
