"""
MSNoiseResult — user-facing class for loading computed results.

Methods are dynamically gated: only methods whose required category is present
in the lineage are exposed.  Accessing a method whose category is absent raises
``AttributeError`` with a clear message, and the method is absent from ``dir()``
and tab-completion.

Typical notebook usage::

    from msnoise.results import MSNoiseResult
    from msnoise.core.db import connect

    db = connect()

    # Build from integer step IDs (most common)
    result = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                    stack=1, refstack=1)

    # Only methods valid for this lineage are visible
    result.get_ccf(...)   # works — stack is in lineage
    result.get_mwcs(...)  # AttributeError — mwcs not in lineage

    # Discover downstream branches (e.g. mwcs_1, stretching_1, wavelet_1)
    for branch in result.branches():
        print(branch.category, branch.lineage_names[-1])

    # Iterate all done results for a category
    for r in MSNoiseResult.list(db, "mwcs_dtt"):
        df = r.get_mwcs_dtt("YA.UV05.00:YA.UV06.00", "ZZ", ("6h", "6h"))
"""

from __future__ import annotations

import glob
import os
from typing import Optional

# ── Step prefix registry ───────────────────────────────────────────────────────

_STEP_PREFIXES = [
    "preprocess", "cc", "psd", "psd_rms", "filter", "stack", "refstack",
    "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
    "stretching", "stretching_dvv",
    "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
]

_PREFIX_TO_KWARG = {p: p for p in _STEP_PREFIXES}


def _step_prefix(step_name: str) -> str:
    """Return the category prefix of a step name, e.g. 'mwcs_dtt_1' -> 'mwcs_dtt'."""
    for prefix in sorted(_STEP_PREFIXES, key=len, reverse=True):
        if step_name.startswith(prefix + "_"):
            return prefix
    raise ValueError(f"Cannot determine category prefix for step name {step_name!r}")


# ── Dynamic method gating ──────────────────────────────────────────────────────

# Registry: method_name -> required category (str) or frozenset of alternatives
_METHOD_CATEGORIES: dict = {}


def _lineage_method(category):
    """Decorator: register a method as requiring *category* in the lineage.

    Pass a list to mean "any of these categories satisfies the requirement"
    (used for get_dvv which serves three DVV step types).
    """
    required = frozenset(category) if isinstance(category, list) else category

    def decorator(fn):
        _METHOD_CATEGORIES[fn.__name__] = required
        return fn

    return decorator


def _category_present(required, present: frozenset) -> bool:
    """Return True if *required* is satisfied by the *present* category set."""
    if isinstance(required, frozenset):
        return bool(required & present)
    return required in present


# ── MSNoiseResult ──────────────────────────────────────────────────────────────

class MSNoiseResult:
    """A resolved pipeline branch with dynamically gated data-loading methods.

    Only ``get_*`` methods whose required step category appears in
    ``lineage_names`` are accessible.  Accessing a gated method raises
    ``AttributeError`` with a helpful message and the method does not appear
    in ``dir()`` or tab-completion.

    Do not instantiate directly — use :meth:`from_ids`, :meth:`from_names`,
    or :meth:`list`.

    Attributes
    ----------
    lineage_names : list[str]
        Ordered step-name strings, e.g.
        ``['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']``.
    category : str
        Category of the terminal step, e.g. ``'refstack'``.
    params : LayeredParams
        Merged configuration parameters for this lineage.
    output_folder : str
        Root output folder.
    """

    def __init__(self, db, lineage_names: list):
        from .core.workflow import resolve_lineage_params
        self._db = db
        self.lineage_names: list = list(lineage_names)
        _, _, self.params = resolve_lineage_params(db, lineage_names)
        self.output_folder: str = self.params.output_folder
        self.category: str = _step_prefix(lineage_names[-1]) if lineage_names else ""
        self._present_categories: frozenset = frozenset(
            _step_prefix(n) for n in lineage_names
        )

    # ── constructors ──────────────────────────────────────────────────────────

    @classmethod
    def from_names(cls, db, names: list) -> "MSNoiseResult":
        """Build from an explicit list of step-name strings."""
        return cls(db, names)

    @classmethod
    def from_ids(cls, db, preprocess=None, cc=None, psd=None, psd_rms=None,
                 filter=None, stack=None, refstack=None, mwcs=None,
                 mwcs_dtt=None, mwcs_dtt_dvv=None, stretching=None,
                 stretching_dvv=None, wavelet=None, wavelet_dtt=None,
                 wavelet_dtt_dvv=None) -> "MSNoiseResult":
        """Build from integer configset IDs in canonical workflow order."""
        kwargs = dict(preprocess=preprocess, cc=cc, psd=psd, psd_rms=psd_rms,
                      filter=filter, stack=stack, refstack=refstack, mwcs=mwcs,
                      mwcs_dtt=mwcs_dtt, mwcs_dtt_dvv=mwcs_dtt_dvv,
                      stretching=stretching, stretching_dvv=stretching_dvv,
                      wavelet=wavelet, wavelet_dtt=wavelet_dtt,
                      wavelet_dtt_dvv=wavelet_dtt_dvv)
        parts = [f"{p}_{v}" for p in _STEP_PREFIXES if (v := kwargs.get(p)) is not None]
        if not parts:
            raise ValueError("At least one step ID must be provided.")
        return cls(db, parts)

    @classmethod
    def list(cls, db, category: str, include_empty: bool = False) -> list:
        """Return all done MSNoiseResult objects for a step category."""
        from .core.workflow import get_done_lineages_for_category, get_workflow_steps
        results = []
        if include_empty:
            for step in get_workflow_steps(db):
                if step.category == category:
                    for names in _upstream_lineage_for_step(db, step):
                        if names:
                            results.append(cls(db, names))
        else:
            for names in get_done_lineages_for_category(db, category):
                results.append(cls(db, names))
        return results

    # ── dynamic gating ────────────────────────────────────────────────────────

    def __dir__(self) -> list:
        base = super().__dir__()
        return [
            name for name in base
            if name not in _METHOD_CATEGORIES
            or _category_present(_METHOD_CATEGORIES[name], self._present_categories)
        ]

    def __getattr__(self, name: str):
        if name in _METHOD_CATEGORIES:
            required = _METHOD_CATEGORIES[name]
            req_str = (" or ".join(sorted(required))
                       if isinstance(required, frozenset) else repr(required))
            present_str = (", ".join(sorted(self._present_categories)) or "(none)")
            raise AttributeError(
                f"{name!r} requires {req_str} in the lineage, "
                f"but this result only covers: {present_str}. "
                f"Use .branches() to navigate to a downstream step."
            )
        raise AttributeError(f"{type(self).__name__!r} has no attribute {name!r}")

    # ── navigation ────────────────────────────────────────────────────────────

    def branches(self, include_empty: bool = False) -> list:
        """Return downstream MSNoiseResult objects one step below this lineage."""
        from .msnoise_table_def import declare_tables
        schema = declare_tables()
        Job = schema.Job
        WorkflowStep = schema.WorkflowStep
        WorkflowLink = schema.WorkflowLink

        terminal_step = (
            self._db.query(WorkflowStep)
            .filter(WorkflowStep.step_name == self.lineage_names[-1])
            .first()
        )
        if terminal_step is None:
            return []

        links = (
            self._db.query(WorkflowLink)
            .filter(WorkflowLink.from_step_id == terminal_step.step_id)
            .filter(WorkflowLink.is_active.is_(True))
            .all()
        )
        results = []
        for link in links:
            child_step = (
                self._db.query(WorkflowStep)
                .filter(WorkflowStep.step_id == link.to_step_id)
                .first()
            )
            if child_step is None:
                continue
            child_names = self.lineage_names + [child_step.step_name]
            if not include_empty:
                has_done = (
                    self._db.query(Job.ref)
                    .filter(Job.step_id == child_step.step_id)
                    .filter(Job.lineage == "/".join(child_names))
                    .filter(Job.flag == "D")
                    .first()
                )
                if has_done is None:
                    continue
            results.append(MSNoiseResult(self._db, child_names))
        return results

    def __repr__(self) -> str:
        valid = [m for m in sorted(_METHOD_CATEGORIES) if m in dir(self)]
        return (
            f"MSNoiseResult(category={self.category!r}, "
            f"lineage={'/'.join(self.lineage_names)!r}, "
            f"available_methods={valid})"
        )

    # ── internal helpers ──────────────────────────────────────────────────────

    def _lineage_upstream_of(self, category: str) -> list:
        result = []
        for name in self.lineage_names:
            if _step_prefix(name) == category:
                break
            result.append(name)
        return result

    def _lineage_through(self, category: str) -> list:
        result = []
        for name in self.lineage_names:
            result.append(name)
            if _step_prefix(name) == category:
                break
        return result

    def _step_name_for(self, category: str) -> str:
        for name in self.lineage_names:
            if _step_prefix(name) == category:
                return name
        raise ValueError(f"No '{category}' step in lineage.")

    def _discover(self, path_pattern: str) -> list:
        return sorted(set(glob.glob(path_pattern)))

    def _load_pair_comp_movstack(self, root, lineage, loader_fn,
                                  pair, components, mov_stack):
        """Generic discovery+load for (root, lineage, sta1, sta2, comp, ms) loaders."""
        base = os.path.join(root, *lineage, "_output")
        results = {}
        mov_stacks = (
            ["%s_%s" % (mov_stack[0], mov_stack[1])] if mov_stack is not None
            else [os.path.basename(p)
                  for p in self._discover(os.path.join(base, "*"))
                  if os.path.isdir(p)]
        )
        for ms_str in mov_stacks:
            ms_tuple = tuple(ms_str.split("_", 1))
            comp_list = (
                [components] if components is not None
                else [os.path.basename(p)
                      for p in self._discover(os.path.join(base, ms_str, "*"))
                      if os.path.isdir(p)]
            )
            for comp in comp_list:
                nc_dir = os.path.join(base, ms_str, comp)
                files = (
                    self._discover(os.path.join(nc_dir, f"{pair.split(':')[0]}_{pair.split(':')[1]}.nc"))
                    if pair is not None
                    else self._discover(os.path.join(nc_dir, "*.nc"))
                )
                for fpath in files:
                    fname = os.path.splitext(os.path.basename(fpath))[0]
                    sta1, sta2 = fname.split("_", 1)
                    pair_key = f"{sta1}:{sta2}"
                    try:
                        data = loader_fn(root, lineage, sta1, sta2, comp, ms_tuple)
                        results[(pair_key, comp, ms_tuple)] = data
                    except Exception:
                        pass
        if pair is not None and components is not None:
            return {k[2]: v for k, v in results.items()}
        if pair is not None and mov_stack is not None:
            return {k[1]: v for k, v in results.items()}
        if components is not None and mov_stack is not None:
            return {k[0]: v for k, v in results.items()}
        return results

    # ── load methods (gated by @_lineage_method) ──────────────────────────────

    @_lineage_method("stack")
    def get_ccf(self, pair=None, components=None, mov_stack=None, format="xarray"):
        """Load CCFs from the stack step. Requires 'stack' in lineage."""
        from .core.io import xr_get_ccf
        from .core.workflow import get_t_axis

        lineage = self._lineage_through("stack")
        taxis = get_t_axis(self.params)
        root = self.output_folder

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_ccf(root, lineage, sta1, sta2, components,
                               mov_stack, taxis, format=format)

        base = os.path.join(root, *lineage, "_output")
        results = {}
        mov_stacks = (
            ["%s_%s" % (mov_stack[0], mov_stack[1])] if mov_stack is not None
            else [os.path.basename(p)
                  for p in self._discover(os.path.join(base, "*")) if os.path.isdir(p)]
        )
        comps = [components] if components is not None else None
        for ms_str in mov_stacks:
            ms_tuple = tuple(ms_str.split("_", 1))
            comp_list = (comps if comps is not None
                         else [os.path.basename(p)
                               for p in self._discover(os.path.join(base, ms_str, "*"))
                               if os.path.isdir(p)])
            for comp in comp_list:
                nc_dir = os.path.join(base, ms_str, comp)
                files = (
                    self._discover(os.path.join(nc_dir, f"{pair.split(':')[0]}_{pair.split(':')[1]}.nc"))
                    if pair is not None
                    else self._discover(os.path.join(nc_dir, "*.nc"))
                )
                for fpath in files:
                    fname = os.path.splitext(os.path.basename(fpath))[0]
                    sta1, sta2 = fname.split("_", 1)
                    pair_key = f"{sta1}:{sta2}"
                    try:
                        df = xr_get_ccf(root, lineage, sta1, sta2, comp,
                                        ms_tuple, taxis, format=format)
                        results[(pair_key, comp, ms_tuple)] = df
                    except Exception:
                        pass
        if pair is not None and components is not None and mov_stack is None:
            return {k[2]: v for k, v in results.items()}
        return results

    @_lineage_method("refstack")
    def get_ref(self, pair=None, components=None, format="xarray"):
        """Load reference stacks. Requires 'refstack' in lineage."""
        from .core.io import xr_get_ref
        from .core.workflow import get_t_axis

        lineage = self._lineage_through("refstack")
        taxis = get_t_axis(self.params)
        root = self.output_folder

        def _fmt(ds):
            if format == "dataframe":
                import pandas as pd
                da = ds.REF
                return pd.Series(da.values, index=da.coords["taxis"].values, name="REF")
            return ds

        if pair is not None and components is not None:
            sta1, sta2 = pair.split(":")
            return _fmt(xr_get_ref(root, lineage, sta1, sta2, components, taxis))

        base = os.path.join(root, *lineage, "_output", "REF")
        results = {}
        comp_list = (
            [components] if components is not None
            else [os.path.basename(p)
                  for p in self._discover(os.path.join(base, "*")) if os.path.isdir(p)]
        )
        for comp in comp_list:
            nc_dir = os.path.join(base, comp)
            files = (
                self._discover(os.path.join(nc_dir, f"{pair.replace(':','_')}.nc"))
                if pair is not None
                else self._discover(os.path.join(nc_dir, "*.nc"))
            )
            for fpath in files:
                fname = os.path.splitext(os.path.basename(fpath))[0]
                sta1, sta2 = fname.split("_", 1)
                try:
                    ds = xr_get_ref(root, lineage, sta1, sta2, comp, taxis)
                    results[(f"{sta1}:{sta2}", comp)] = _fmt(ds)
                except Exception:
                    pass
        if pair is not None:
            return {k[1]: v for k, v in results.items()}
        return results

    @_lineage_method("mwcs")
    def get_mwcs(self, pair=None, components=None, mov_stack=None, format="xarray"):
        """Load MWCS results. Requires 'mwcs' in lineage."""
        from .core.io import xr_get_mwcs
        lineage = self._lineage_through("mwcs")
        root = self.output_folder
        fmt = "dataset" if format == "xarray" else "dataframe"
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_mwcs(root, lineage, sta1, sta2, components, mov_stack, format=fmt)
        return self._load_pair_comp_movstack(
            root, lineage, lambda *a, **kw: xr_get_mwcs(*a, format=fmt, **kw),
            pair, components, mov_stack)

    @_lineage_method("mwcs_dtt")
    def get_mwcs_dtt(self, pair=None, components=None, mov_stack=None, format="xarray"):
        """Load MWCS-DTT results. Requires 'mwcs_dtt' in lineage."""
        from .core.io import xr_get_dtt
        lineage = self._lineage_through("mwcs_dtt")
        root = self.output_folder
        fmt = "dataset" if format == "xarray" else "dataframe"
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_dtt(root, lineage, sta1, sta2, components, mov_stack, format=fmt)
        return self._load_pair_comp_movstack(
            root, lineage, lambda *a, **kw: xr_get_dtt(*a, format=fmt, **kw),
            pair, components, mov_stack)

    @_lineage_method("stretching")
    def get_stretching(self, pair=None, components=None, mov_stack=None, format="xarray"):
        """Load stretching results. Requires 'stretching' in lineage."""
        from .core.io import _xr_get_stretching
        lineage = self._lineage_through("stretching")
        root = self.output_folder
        fmt = "dataset" if format == "xarray" else "dataframe"
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return _xr_get_stretching(root, lineage, sta1, sta2, components, mov_stack, format=fmt)
        return self._load_pair_comp_movstack(
            root, lineage, lambda *a, **kw: _xr_get_stretching(*a, format=fmt, **kw),
            pair, components, mov_stack)

    @_lineage_method(["mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"])
    def get_dvv(self, pair_type="ALL", components=None, mov_stack=None, format="xarray"):
        """Load DVV aggregate results. Requires a DVV step in lineage."""
        from .core.io import xr_get_dvv_agg
        dvv_cat   = self.category
        lineage   = self._lineage_through(dvv_cat)
        step_name = self._step_name_for(dvv_cat)
        root      = self.output_folder
        fmt       = "dataset" if format == "xarray" else "dataframe"
        if components is not None and mov_stack is not None:
            return xr_get_dvv_agg(root, lineage, step_name, mov_stack,
                                   pair_type, components, format=fmt)
        base = os.path.join(root, *lineage, step_name, "_output")
        results = {}
        ms_dirs = (
            [os.path.join(base, "%s_%s" % (mov_stack[0], mov_stack[1]))]
            if mov_stack is not None
            else [p for p in self._discover(os.path.join(base, "*")) if os.path.isdir(p)]
        )
        for ms_dir in ms_dirs:
            ms_str   = os.path.basename(ms_dir)
            ms_tuple = tuple(ms_str.split("_", 1))
            pattern  = (f"dvv_{pair_type}_*.nc" if components is None
                        else f"dvv_{pair_type}_{components}.nc")
            for nc_file in self._discover(os.path.join(ms_dir, pattern)):
                fname = os.path.splitext(os.path.basename(nc_file))[0]
                parts = fname.split("_", 2)
                if len(parts) < 3:
                    continue
                _, pt, comp = parts
                try:
                    results[(pt, comp, ms_tuple)] = xr_get_dvv_agg(
                        root, lineage, step_name, ms_tuple, pt, comp, format=fmt)
                except Exception:
                    pass
        if components is not None and mov_stack is not None:
            return {k[0]: v for k, v in results.items()}
        if components is not None:
            return {(k[0], k[2]): v for k, v in results.items()}
        if mov_stack is not None:
            return {(k[0], k[1]): v for k, v in results.items()}
        return results

    @_lineage_method("wavelet")
    def get_wct(self, pair=None, components=None, mov_stack=None, format="xarray"):
        """Load WCT results. Requires 'wavelet' in lineage."""
        from .core.io import xr_load_wct
        lineage = self._lineage_through("wavelet")
        root = self.output_folder
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            ds = xr_load_wct(root, lineage, sta1, sta2, components, mov_stack)
            return ds.to_dataframe() if format == "dataframe" else ds
        return self._load_pair_comp_movstack(root, lineage, xr_load_wct,
                                              pair, components, mov_stack)

    @_lineage_method("wavelet_dtt")
    def get_wct_dtt(self, pair=None, components=None, mov_stack=None, format="xarray"):
        """Load WCT dt/t results. Requires 'wavelet_dtt' in lineage."""
        from .core.io import xr_get_wct_dtt
        lineage = self._lineage_through("wavelet_dtt")
        root = self.output_folder
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            ds = xr_get_wct_dtt(root, lineage, sta1, sta2, components, mov_stack)
            return ds.to_dataframe() if format == "dataframe" else ds
        return self._load_pair_comp_movstack(root, lineage, xr_get_wct_dtt,
                                              pair, components, mov_stack)

    @_lineage_method("psd")
    def get_psd(self, seed_id=None, day=None, format="xarray"):
        """Load PSD results. Requires 'psd' in lineage."""
        from .core.io import xr_load_psd
        step_name = self._step_name_for("psd")
        lineage   = self._lineage_upstream_of("psd")
        root      = self.output_folder
        fmt       = "dataset" if format == "xarray" else "dataframe"
        if seed_id is not None and day is not None:
            return xr_load_psd(root, lineage, step_name, seed_id, day, format=fmt)
        base = os.path.join(root, *lineage, step_name, "_output", "daily")
        results = {}
        seed_ids = (
            [seed_id] if seed_id is not None
            else [os.path.basename(p)
                  for p in self._discover(os.path.join(base, "*")) if os.path.isdir(p)]
        )
        for sid in seed_ids:
            day_files = (
                self._discover(os.path.join(base, sid, f"{day}.nc"))
                if day is not None
                else self._discover(os.path.join(base, sid, "*.nc"))
            )
            for fpath in day_files:
                day_key = os.path.splitext(os.path.basename(fpath))[0]
                r = xr_load_psd(root, lineage, step_name, sid, day_key, format=fmt)
                if r is not None:
                    results[(sid, day_key)] = r
        if seed_id is not None:
            return {k[1]: v for k, v in results.items()}
        return results

    @_lineage_method("psd_rms")
    def get_psd_rms(self, seed_id=None, format="xarray"):
        """Load PSD RMS results. Requires 'psd_rms' in lineage."""
        from .core.io import xr_load_rms
        step_name = self._step_name_for("psd_rms")
        lineage   = self._lineage_upstream_of("psd_rms")
        root      = self.output_folder
        fmt       = "dataset" if format == "xarray" else "dataframe"
        if seed_id is not None:
            return xr_load_rms(root, lineage, step_name, seed_id, format=fmt)
        base = os.path.join(root, *lineage, step_name, "_output")
        results = {}
        for sid_dir in self._discover(os.path.join(base, "*")):
            if not os.path.isdir(sid_dir):
                continue
            sid = os.path.basename(sid_dir)
            r = xr_load_rms(root, lineage, step_name, sid, format=fmt)
            if r is not None:
                results[sid] = r
        return results


# ── module-level helper ────────────────────────────────────────────────────────

def _upstream_lineage_for_step(db, step) -> list:
    """Walk WorkflowLinks upstream from *step*, return all root-to-step paths."""
    from .core.workflow import get_workflow_links
    from .msnoise_table_def import declare_tables
    schema = declare_tables()
    WorkflowStep = schema.WorkflowStep
    links = get_workflow_links(db)
    parent_map: dict = {}
    for lnk in links:
        parent_map.setdefault(lnk.to_step_id, []).append(lnk.from_step_id)

    def _walk(step_id, visited):
        step_obj = db.query(WorkflowStep).filter(
            WorkflowStep.step_id == step_id).first()
        if step_obj is None or step_obj.category == "global":
            return [[]]
        if step_id in visited:
            return [[step_obj.step_name]]
        visited = visited | {step_id}
        parent_ids = parent_map.get(step_id, [])
        if not parent_ids:
            return [[step_obj.step_name]]
        paths = []
        for parent_id in parent_ids:
            for upstream_path in _walk(parent_id, visited):
                paths.append(upstream_path + [step_obj.step_name])
        return paths

    return [p for p in _walk(step.step_id, set()) if p]
