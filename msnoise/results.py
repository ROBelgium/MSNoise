"""
MSNoiseResult — user-facing class for loading computed results.

``MSNoiseResult`` is the **recommended entry point** for reading any output
produced by the MSNoise pipeline from a notebook or script.  You do not need
to know the on-disk path layout or the low-level ``core.io`` functions: just
build a result object for your lineage and call the appropriate ``get_*``
method.

Two construction paths exist depending on whether a database is available:

* **DB-connected** (active project): use :meth:`~MSNoiseResult.from_ids`,
  :meth:`~MSNoiseResult.from_names`, or :meth:`~MSNoiseResult.list`.
* **DB-free** (portable bundle): use :meth:`~MSNoiseResult.from_bundle` —
  no database, no ``msnoise db init``, works anywhere the msnoise package is
  installed.

For the full reference see :class:`msnoise.results.MSNoiseResult`.

Dynamic method gating
---------------------

Methods are dynamically gated: only methods whose required step category is
present in the lineage are exposed.  Accessing an unavailable method raises
``AttributeError`` with a clear message, and the method is absent from
``dir()`` and tab-completion, so Jupyter's autocomplete only shows what is
actually readable.

Constructing a result object (DB-connected)
--------------------------------------------

Use :meth:`MSNoiseResult.from_ids` with the integer config-set numbers for
each step you want to include.  You only need to go as far down the pipeline
as the data you want to read::

    from msnoise.results import MSNoiseResult
    from msnoise.core.db import connect

    db = connect()

    # ----- read stacked CCFs and everything downstream -----
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                               stack=1, refstack=1)
    da = r.get_ccf("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", ("1D", "1D"))

    # ----- read raw (pre-stack) CC outputs only -----
    # Initialise only to filter level — get_ccf_raw is still available.
    r_cc = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)
    da = r_cc.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ",
                           date="2023-01-01", kind="all")

    # ----- read all the way to dv/v -----
    r_dvv = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                   stack=1, refstack=1,
                                   mwcs=1, mwcs_dtt=1, mwcs_dtt_dvv=1)
    ds = r_dvv.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D", "1D"))

    # Only methods valid for this lineage are visible
    r_cc.get_ccf(...)      # AttributeError — 'stack' not in lineage
    r_cc.get_ccf_raw(...)  # works — 'cc' is in lineage

Available get_* methods
------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Method
     - Requires in lineage
     - Returns
   * - ``get_ccf_raw``
     - ``cc``
     - Raw per-window or daily-stacked CCFs written by ``s03``
   * - ``get_ccf``
     - ``stack``
     - Moving-stacked CCFs written by ``s04``
   * - ``get_ref``
     - ``refstack``
     - Reference stacks written by ``s04_stack_refstack``
   * - ``get_mwcs``
     - ``mwcs``
     - MWCS results written by ``s05``
   * - ``get_mwcs_dtt``
     - ``mwcs_dtt``
     - MWCS dt/t results written by ``s06``
   * - ``get_stretching``
     - ``stretching``
     - Stretching results written by ``s10``
   * - ``get_wct``
     - ``wavelet``
     - Wavelet coherence results written by ``s08``
   * - ``get_wct_dtt``
     - ``wavelet_dtt``
     - WCT dt/t results written by ``s09``
   * - ``get_dvv``
     - any dvv step
     - Aggregated dv/v written by ``s07``
   * - ``get_psd``
     - ``psd``
     - Daily PSD written by ``psd_compute``
   * - ``get_psd_rms``
     - ``psd_rms``
     - PSD RMS written by ``psd_compute_rms``

Reading raw CC outputs (pre-stack)
------------------------------------

``get_ccf_raw`` reads the files written directly by ``s03_compute_no_rotation``
*before* any stacking.  Two storage layouts exist depending on the config:

* ``kind="all"`` — per-window CCFs under ``_output/all/`` with dims
  ``(times, taxis)`` where ``times`` is the window start time.
* ``kind="daily"`` — daily-stacked CCFs under ``_output/daily/`` with dim
  ``(taxis,)`` per file.

The method only requires ``cc`` in the lineage, so you can call it on a
filter-level result without needing ``stack`` or anything further::

    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)

    # Single day — per-window DataArray with dims (times, taxis)
    da = r.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ",
                        date="2023-01-01", kind="all")

    # All days for one pair, daily stacks — dict keyed by date string
    d = r.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", kind="daily")

    # Concatenate all days into a single DataArray (times, taxis)
    import xarray as xr
    da_all = xr.concat(d.values(), dim="times").sortby("times")

Navigating branches
--------------------

When a lineage has more than one downstream branch (e.g. ``mwcs_1`` *and*
``stretching_1`` both follow ``refstack_1``), use :meth:`branches` to
enumerate them::

    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                               stack=1, refstack=1)
    for branch in r.branches():
        print(branch)   # MSNoiseResult(category='mwcs', ...)

``branches()`` works identically on bundle-loaded results — it uses the DAG
topology from :func:`~msnoise.core.workflow.get_workflow_chains` and confirms
each branch by checking for an ``_output`` subdirectory on disk.

Iterating over all completed results
--------------------------------------

:meth:`list` returns every ``MSNoiseResult`` for which at least one Done job
exists in the given category::

    for r in MSNoiseResult.list(db, "mwcs_dtt"):
        ds = r.get_mwcs_dtt("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", ("1D", "1D"))

Exporting dv/v with provenance
--------------------------------

:meth:`export_dvv` writes dv/v NetCDF files that embed full provenance
(lineage string, all config parameters, station metadata)::

    r = MSNoiseResult.list(db, "mwcs_dtt_dvv")[0]
    written = r.export_dvv("exports/")
    for f in written:
        print(f)  # dvv_CC_ZZ__pre1-cc1-f1-stk1-ref1-mwcs1-dtt1-dvv1__m1D-1D.nc

    # Reload and inspect the embedded params
    import xarray as xr
    from msnoise.params import MSNoiseParams
    ds = xr.open_dataset(written[0])
    params = MSNoiseParams.from_yaml_string(ds.attrs["msnoise_params"])
    print(params.mwcs.mwcs_wlen)

Portable bundles — HPC to laptop, or journal supplementary data
----------------------------------------------------------------

:meth:`export_bundle` packages a selectable slice of the computed ``_output/``
tree into a self-describing directory (or ``.zip``), together with
``params.yaml`` (full lineage config) and ``MANIFEST.json`` (sha256 per file,
station list, original output path).  The bundle can be read on any machine
with MSNoise installed — **no database required**.

**Selecting what to include** with the ``from_step=`` argument:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - ``from_step=``
     - What is copied
   * - ``"stack"``
     - Moving-stacked CCFs + all downstream steps (large)
   * - ``"refstack"``
     - Reference stacks + all dv/v branches (moderate; **default for CC workflows**)
   * - ``"mwcs_dtt_dvv"``
     - Only the final MWCS dv/v aggregates (small; ideal for paper supplements)
   * - ``"psd"``
     - PSD daily files + RMS (moderate; **default for PSD-only workflows**)
   * - ``"cc"`` / ``"filter"``
     - Raw per-window CCFs + everything downstream (very large; opt-in)
   * - ``None``
     - Earliest step with ``_output`` in the lineage (auto-detected)

The full lineage directory nesting is preserved verbatim inside the bundle
root, so :meth:`from_bundle` can set ``output_folder = bundle_root`` and all
:func:`~msnoise.core.io.xr_get_*` calls resolve paths identically — nothing
in ``core/io.py`` needed to change.

**Typical HPC → laptop workflow**::

    # ── On the HPC (has DB) ───────────────────────────────────────────────
    from msnoise.results import MSNoiseResult
    from msnoise.core.db import connect

    db = connect()
    r = MSNoiseResult.list(db, "mwcs_dtt_dvv")[0]

    # Export reference stacks + all dv/v branches as a zip:
    bundle = r.export_bundle(
        "belgium_2023/",
        from_step="refstack",
        compress=True,          # → belgium_2023.zip  (~50 MB typical network)
    )

    # Export only the final dv/v for a journal data-availability statement:
    r.export_bundle("paper_SI/", from_step="mwcs_dtt_dvv", compress=True)

    # ── rsync to local machine ────────────────────────────────────────────
    # rsync -avz hpc:/scratch/user/msnoise/belgium_2023.zip ./

    # ── On the laptop / at the reviewer's machine (no DB needed) ─────────
    from msnoise.results import MSNoiseResult

    r2 = MSNoiseResult.from_bundle("belgium_2023.zip")

    # Optional integrity check (recommended for published data):
    r2.verify()
    # OK — 347 files verified.

    # Full MSNoiseResult API works immediately:
    ds  = r2.get_dvv("CC", "ZZ", ("1D", "1D"))
    da  = r2.get_ccf("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", ("1D", "1D"))
    ref = r2.get_ref("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ")

    # Navigate branches (folder scan — no DB):
    for branch in r2.branches():
        print(branch)
    # MSNoiseResult(category='mwcs_dtt_dvv', lineage='...', available_methods=[...])
    # MSNoiseResult(category='stretching_dvv', lineage='...', available_methods=[...])

    # Inspect the exact parameters that produced the result:
    print(r2.params.mwcs.mwcs_wlen)        # → 10.0
    print(r2.params.cc.cc_sampling_rate)   # → 20.0

    # params.yaml is also preserved in MANIFEST.json alongside the original
    # HPC output_folder path, msnoise version, and generation timestamp —
    # sufficient for a data-availability or reproducibility statement.
"""

from __future__ import annotations

import glob
import os

# ── Step prefix registry ───────────────────────────────────────────────────────
# Derived from the plugin-aware workflow order so plugin categories are
# automatically included without editing this file.

def _get_step_prefixes():
    from .core.workflow import get_workflow_order
    return get_workflow_order()

_PREFIX_TO_KWARG = {p: p for p in _get_step_prefixes()}


def _step_prefix(step_name: str) -> str:
    """Return the category prefix of a step name, e.g. 'mwcs_dtt_1' -> 'mwcs_dtt'."""
    for prefix in sorted(_get_step_prefixes(), key=len, reverse=True):
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
    params : MSNoiseParams
        Merged configuration parameters for this lineage.
    output_folder : str
        Root output folder.
    """

    def __init__(self, db, lineage_names: list, _params=None):
        self._db = db  # None in bundle / DB-free mode
        self.lineage_names: list = list(lineage_names)
        if _params is not None:
            self.params = _params
        else:
            from .core.workflow import resolve_lineage_params
            _, _, self.params = resolve_lineage_params(db, lineage_names)
        self.output_folder: str = self.params.global_.output_folder
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
        parts = [f"{p}_{v}" for p in _get_step_prefixes() if (v := kwargs.get(p)) is not None]
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

    @classmethod
    def from_bundle(cls, path: str) -> "MSNoiseResult":
        """Load a read-only :class:`MSNoiseResult` from a bundle directory or
        ``.zip`` file produced by :meth:`export_bundle`.

        No database connection is required.  All ``get_*`` methods work
        immediately after construction.  :meth:`branches` uses a folder scan
        rather than the DB.

        :param path: Path to a bundle directory **or** a ``.zip`` file.
        :returns: :class:`MSNoiseResult` with ``_db=None``.
        :raises FileNotFoundError: if ``params.yaml`` is absent from the bundle.

        Example — directory bundle::

            r = MSNoiseResult.from_bundle("/data/msnoise_bundle/")
            ds = r.get_dvv("CC", "ZZ", ("1D", "1D"))

        Example — zip bundle::

            r = MSNoiseResult.from_bundle("/data/bundle.zip")
            da = r.get_ccf("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", ("1D", "1D"))
            r.verify()
        """
        import os
        import zipfile as _zf
        import tempfile as _tmp

        # ── normalise: zip → extract to temp dir ──────────────────────────
        _tmpdir = None  # TemporaryDirectory — kept alive on the instance

        if isinstance(path, str) and path.endswith(".zip") and os.path.isfile(path):
            _tmpdir = _tmp.TemporaryDirectory(prefix="msnoise_bundle_")
            with _zf.ZipFile(path) as zf:
                zf.extractall(_tmpdir.name)
            # params.yaml may be at the zip root or one level inside a
            # top-level directory (the common "zip contains a single folder"
            # pattern produced by export_bundle compress=True).
            if os.path.isfile(os.path.join(_tmpdir.name, "params.yaml")):
                bundle_root = _tmpdir.name
            else:
                candidates = [
                    os.path.join(_tmpdir.name, d)
                    for d in os.listdir(_tmpdir.name)
                    if os.path.isfile(
                        os.path.join(_tmpdir.name, d, "params.yaml")
                    )
                ]
                if not candidates:
                    raise FileNotFoundError(
                        "params.yaml not found inside the zip bundle."
                    )
                bundle_root = candidates[0]
        else:
            bundle_root = os.path.abspath(path)

        # ── load params.yaml ──────────────────────────────────────────────
        yaml_path = os.path.join(bundle_root, "params.yaml")
        if not os.path.isfile(yaml_path):
            raise FileNotFoundError(
                f"params.yaml not found at {yaml_path!r}. "
                "Is this a valid MSNoise bundle directory?"
            )

        from .params import MSNoiseParams
        from obspy.core.util.attribdict import AttribDict

        params = MSNoiseParams.from_yaml(yaml_path)

        # Override output_folder with the bundle's actual location on this
        # machine.  The original HPC path is preserved in MANIFEST.json
        # (output_folder_original) for provenance.
        layers = object.__getattribute__(params, "_layers")
        if "global" in layers:
            global_dict = dict(layers["global"])
            global_dict["output_folder"] = bundle_root
            layers["global"] = AttribDict(global_dict)

        # ── construct without DB ──────────────────────────────────────────
        inst = cls.__new__(cls)
        inst._db = None
        inst._bundle_root = bundle_root
        inst._tmpdir = _tmpdir        # keeps TemporaryDirectory alive
        inst.lineage_names = list(params.lineage_names)
        inst.params = params
        inst.output_folder = bundle_root
        inst.category = (
            _step_prefix(inst.lineage_names[-1]) if inst.lineage_names else ""
        )
        inst._present_categories = frozenset(
            _step_prefix(n) for n in inst.lineage_names
        )
        return inst

    # ── dynamic gating ────────────────────────────────────────────────────────

    def __dir__(self) -> list:
        base = super().__dir__()
        return [
            name for name in base
            if name not in _METHOD_CATEGORIES
            or _category_present(_METHOD_CATEGORIES[name], self._present_categories)
        ]

    def __getattribute__(self, name: str):
        """Gate registered methods whose required category is absent from the lineage.

        Must use __getattribute__ (not __getattr__) because the get_* methods
        ARE defined on the class — __getattr__ is only the fallback for missing
        attributes and would never be called for them.
        """
        if name in _METHOD_CATEGORIES:
            present = object.__getattribute__(self, "_present_categories")
            required = _METHOD_CATEGORIES[name]
            if not _category_present(required, present):
                req_str = (" or ".join(sorted(required))
                           if isinstance(required, frozenset) else repr(required))
                present_str = ", ".join(sorted(present)) or "(none)"
                raise AttributeError(
                    f"{name!r} requires category {req_str} in the lineage, "
                    f"but this result only covers: {present_str}. "
                    f"Use .branches() to navigate to a downstream step."
                )
        return object.__getattribute__(self, name)

    # ── navigation ────────────────────────────────────────────────────────────

    def branches(self, include_empty: bool = False) -> list:
        """Return downstream MSNoiseResult objects one step below this lineage.

        In bundle / DB-free mode (``_db is None``) the DAG topology is taken
        from :func:`~msnoise.core.workflow.get_workflow_chains` and child
        directories are checked for the presence of an ``_output`` subdirectory
        to confirm they were actually computed.
        """
        if self._db is None:
            return self._branches_from_folders(include_empty=include_empty)

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
        # When branching from a refstack result to mwcs/stretching/wavelet, the
        # actual job lineage encodes both parents: …/stack_N/refstack_M/dvv_step.
        # (lineage convention A1 — refstack and stack are siblings, both children
        # of filter, but the mwcs lineage string contains both to distinguish
        # which stack CCFs were used alongside which reference.)
        # We look up the sibling stack step(s) via WorkflowLinks from the same
        # filter parent, then enumerate all (stack_N, dvv_step) combinations.
        terminal_cat = terminal_step.category
        _dvv_cats = {"mwcs", "stretching", "wavelet"}

        results = []
        for link in links:
            child_step = (
                self._db.query(WorkflowStep)
                .filter(WorkflowStep.step_id == link.to_step_id)
                .first()
            )
            if child_step is None:
                continue

            # Build candidate lineage(s) for this child step.
            if terminal_cat == "refstack" and child_step.category in _dvv_cats:
                # Find sibling stack steps: shared filter parent → stack children.
                # The filter step is the parent of refstack in the new topology.
                filter_parent_ids = (
                    self._db.query(WorkflowLink.from_step_id)
                    .filter(WorkflowLink.to_step_id == terminal_step.step_id)
                    .filter(WorkflowLink.is_active.is_(True))
                    .all()
                )
                sibling_stack_names = []
                for (fp_id,) in filter_parent_ids:
                    fp_step = self._db.query(WorkflowStep).filter(
                        WorkflowStep.step_id == fp_id).first()
                    if fp_step and fp_step.category == "filter":
                        # Stack steps linked from this filter
                        stk_links = (
                            self._db.query(WorkflowLink)
                            .filter(WorkflowLink.from_step_id == fp_id)
                            .filter(WorkflowLink.is_active.is_(True))
                            .all()
                        )
                        for sl in stk_links:
                            s = self._db.query(WorkflowStep).filter(
                                WorkflowStep.step_id == sl.to_step_id).first()
                            if s and s.category == "stack":
                                sibling_stack_names.append(s.step_name)

                if not sibling_stack_names:
                    # Fallback: no stack sibling found, use direct lineage
                    candidate_lineages = [self.lineage_names + [child_step.step_name]]
                else:
                    # One branch per stack_N sibling
                    candidate_lineages = []
                    filter_prefix = self.lineage_names[:-1]  # strip refstack_M
                    for stk_name in sibling_stack_names:
                        candidate_lineages.append(
                            filter_prefix + [stk_name,
                                             terminal_step.step_name,
                                             child_step.step_name]
                        )
            else:
                candidate_lineages = [self.lineage_names + [child_step.step_name]]

            for child_names in candidate_lineages:
                if not include_empty:
                    from .msnoise_table_def import Lineage as _Lineage
                    from sqlalchemy.orm import aliased as _aliased
                    _lin_alias = _aliased(_Lineage)
                    _lin_str = "/".join(child_names)
                    has_done = (
                        self._db.query(Job.ref)
                        .filter(Job.step_id == child_step.step_id)
                        .join(_lin_alias, Job.lineage_id == _lin_alias.lineage_id)
                        .filter(_lin_alias.lineage_str == _lin_str)
                        .filter(Job.flag == "D")
                        .first()
                    )
                    if has_done is None:
                        continue
                results.append(MSNoiseResult(self._db, child_names))
        return results

    def _branches_from_folders(self, include_empty: bool = False) -> list:
        """Discover child MSNoiseResult objects via folder scan (no DB).

        Uses :func:`~msnoise.core.workflow.get_workflow_chains` for DAG
        topology, then checks whether each expected child directory contains
        an ``_output`` subdirectory to confirm it was actually computed.
        """
        import os
        from .core.workflow import get_workflow_chains

        chains = get_workflow_chains()
        terminal_cat = self.category
        next_cats = chains.get(terminal_cat, {}).get("next_steps", [])
        if not next_cats:
            return []

        # Infer set number from the terminal step name, e.g. "refstack_1" → "1".
        # All steps in a single exported lineage share the same set number, so
        # this is reliable for bundles.  For mixed-set lineages the folder
        # existence check is the safety net.
        terminal_name = self.lineage_names[-1]
        set_num = terminal_name.rsplit("_", 1)[-1]

        results = []
        for child_cat in next_cats:
            child_step_name = f"{child_cat}_{set_num}"
            child_names = self.lineage_names + [child_step_name]
            child_dir = os.path.join(self.output_folder, *child_names)

            # Fallback: if the inferred name doesn't exist, glob for any
            # matching category directory (handles mixed set numbers).
            if not os.path.isdir(child_dir):
                import glob as _glob
                pattern = os.path.join(
                    self.output_folder, *self.lineage_names,
                    f"{child_cat}_*",
                )
                matches = sorted(
                    d for d in _glob.glob(pattern) if os.path.isdir(d)
                )
                if not matches:
                    if include_empty:
                        # Emit a placeholder with the inferred name
                        pass
                    else:
                        continue
                else:
                    child_dir = matches[0]
                    child_step_name = os.path.basename(child_dir)
                    child_names = self.lineage_names + [child_step_name]

            child_output = os.path.join(child_dir, "_output")
            if not include_empty and not os.path.isdir(child_output):
                continue

            child = MSNoiseResult.__new__(MSNoiseResult)
            child._db = None
            child._bundle_root = getattr(self, "_bundle_root", self.output_folder)
            child._tmpdir = getattr(self, "_tmpdir", None)
            child.lineage_names = child_names
            child.params = self.params          # shared MSNoiseParams (immutable)
            child.output_folder = self.output_folder
            child.category = child_cat
            child._present_categories = frozenset(
                _step_prefix(n) for n in child_names
            )
            results.append(child)
        return results

    def verify(self, verbose: bool = False) -> bool:
        """Verify bundle integrity against ``MANIFEST.json``.

        Re-hashes every file listed in the manifest and compares to the stored
        sha256 digest.  Prints a summary line; returns ``True`` if all pass.

        :param verbose: If ``True``, print every file path as it is checked.
        :returns: ``True`` if all checksums match, ``False`` otherwise.
        :raises FileNotFoundError: if ``MANIFEST.json`` is absent (not a
            bundle, or loaded via :meth:`from_ids` / :meth:`from_names`).

        Example::

            r = MSNoiseResult.from_bundle("/data/bundle/")
            assert r.verify(), "Bundle integrity check failed"
        """
        import hashlib
        import json

        manifest_path = os.path.join(self.output_folder, "MANIFEST.json")
        if not os.path.isfile(manifest_path):
            raise FileNotFoundError(
                f"MANIFEST.json not found at {manifest_path!r}. "
                "verify() is only available for bundles loaded via from_bundle()."
            )
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)

        failures = []
        entries = manifest.get("files", [])
        for entry in entries:
            fpath = os.path.join(self.output_folder, entry["path"])
            if not os.path.isfile(fpath):
                failures.append((entry["path"], "file missing"))
                continue
            h = hashlib.sha256()
            with open(fpath, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
            digest = h.hexdigest()
            if digest != entry["sha256"]:
                failures.append(
                    (entry["path"], f"sha256 mismatch (got {digest[:12]}…)")
                )
            elif verbose:
                print(f"  OK  {entry['path']}")

        total = len(entries)
        if failures:
            print(f"FAILED {len(failures)}/{total} files:")
            for path, reason in failures:
                print(f"  FAIL  {path}: {reason}")
            return False
        print(f"OK — {total} files verified.")
        return True

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
    def get_ccf(self, pair=None, components=None, mov_stack=None):
        """Load CCFs from the stack step. Requires 'stack' in lineage.

        :returns: :class:`xarray.DataArray` (single pair/comp/ms) or dict of
            DataArrays keyed by ``(pair, comp, mov_stack)``.
        """
        from .core.io import xr_get_ccf
        from .core.workflow import get_t_axis

        lineage = self._lineage_through("stack")
        taxis = get_t_axis(self.params)
        root = self.output_folder

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_ccf(root, lineage, sta1, sta2, components, mov_stack, taxis)

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
                        results[(pair_key, comp, ms_tuple)] = xr_get_ccf(
                            root, lineage, sta1, sta2, comp, ms_tuple, taxis)
                    except Exception:
                        pass
        if pair is not None and components is not None and mov_stack is None:
            return {k[2]: v for k, v in results.items()}
        return results

    @_lineage_method("cc")
    def get_ccf_raw(self, pair=None, components=None, date=None, kind="all"):
        """Load raw CC step outputs (per-window or daily-stacked CCFs).

        Reads the files written by ``s03_compute_no_rotation`` **before**
        stacking.  Requires only ``'cc'`` in the lineage, so a result
        initialised down to ``filter_N`` (e.g.
        ``MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)``) exposes
        this method.

        Path layout on disk::

            <o> / preprocess_1 / cc_1 / filter_1 / _output / all|daily
                         / <comp> / <sta1>_<sta2> / <YYYY-MM-DD>.nc

        Parameters
        ----------
        pair : str or None
            Station pair in ``"NET.STA.LOC:NET.STA.LOC"`` format.
            When *None* all available pairs are returned.
        components : str or None
            Component string, e.g. ``"ZZ"``.  *None* = all.
        date : str or None
            ISO date ``"YYYY-MM-DD"``.  *None* = all available dates.
        kind : {"all", "daily"}
            ``"all"``   — per-window CCFs (``_output/all/``), dims
            ``(times, taxis)`` per file.
            ``"daily"`` — daily-stacked CCFs (``_output/daily/``), dim
            ``(taxis,)`` per file, expanded with a ``times`` coordinate so
            results can be concatenated by the caller.

        Returns
        -------
        xarray.DataArray
            When *pair*, *components* **and** *date* are all specified.
        dict
            Otherwise: keys are ``(pair_key, comp, date_str)`` tuples,
            collapsed to 1-tuples for any dimension fixed by the caller.

        Raises
        ------
        ValueError
            If no ``cc_*`` or filter step is found in the lineage.

        Example::

            r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)

            # Single file — per-window DataArray (times, taxis)
            da = r.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ",
                                date="2023-01-01", kind="all")

            # All dates for one pair/comp — daily stacks
            d = r.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", kind="daily")
            # d is {date_str: DataArray(times=[T], taxis=...)}
        """
        import os
        import glob as _glob

        import pandas as _pd
        from .core.io import xr_get_ccf_all, xr_get_ccf_daily

        # ── derive cc_lineage and filter_step ────────────────────────────────
        cc_idx = next(
            (i for i, n in enumerate(self.lineage_names) if n.startswith("cc_")),
            None,
        )
        if cc_idx is None:
            raise ValueError("No 'cc_*' step found in lineage_names.")
        cc_lineage = self.lineage_names[: cc_idx + 1]   # e.g. ['preprocess_1', 'cc_1']
        if cc_idx + 1 >= len(self.lineage_names):
            raise ValueError(
                "No filter step found after cc_* in lineage_names. "
                "CC outputs are stored under the filter step name; "
                "initialise with at least filter=1."
            )
        filter_step = self.lineage_names[cc_idx + 1]    # e.g. 'filter_1'

        root = self.output_folder
        loader = xr_get_ccf_all if kind == "all" else xr_get_ccf_daily

        # ── fast path: fully specified ────────────────────────────────────────
        if pair is not None and components is not None and date is not None:
            sta1, sta2 = pair.split(":")
            return loader(root, cc_lineage, filter_step, sta1, sta2, components, date)

        # ── discovery ────────────────────────────────────────────────────────
        kind_dir = os.path.join(root, *cc_lineage, filter_step, "_output", kind)
        comp_list = (
            [components] if components is not None
            else [os.path.basename(p)
                  for p in sorted(_glob.glob(os.path.join(kind_dir, "*")))
                  if os.path.isdir(p)]
        )

        results = {}
        for comp in comp_list:
            if pair is not None:
                sta1_raw, sta2_raw = pair.split(":")
                pair_dirs = [os.path.join(kind_dir, comp, f"{sta1_raw}_{sta2_raw}")]
            else:
                pair_dirs = sorted(
                    p for p in _glob.glob(os.path.join(kind_dir, comp, "*"))
                    if os.path.isdir(p)
                )

            for pair_dir in pair_dirs:
                if not os.path.isdir(pair_dir):
                    continue
                dir_name = os.path.basename(pair_dir)
                sta1, sta2 = dir_name.split("_", 1)
                pair_key = f"{sta1}:{sta2}"

                date_files = (
                    [os.path.join(pair_dir, f"{date}.nc")]
                    if date is not None
                    else sorted(_glob.glob(os.path.join(pair_dir, "*.nc")))
                )

                for fpath in date_files:
                    if not os.path.isfile(fpath):
                        continue
                    date_str = os.path.splitext(os.path.basename(fpath))[0]
                    try:
                        da = loader(root, cc_lineage, filter_step,
                                    sta1, sta2, comp, date_str)
                        if kind == "daily":
                            # Expand so callers can xr.concat along times
                            da = da.expand_dims(
                                {"times": [_pd.Timestamp(date_str)]}
                            )
                        results[(pair_key, comp, date_str)] = da
                    except Exception:
                        pass

        # ── collapse fixed dimensions ─────────────────────────────────────────
        if pair is not None and components is not None:
            # Only date varies → {date_str: DataArray}
            return {k[2]: v for k, v in results.items()}
        if pair is not None and date is not None:
            # Only comp varies → {comp: DataArray}
            return {k[1]: v for k, v in results.items()}
        if components is not None and date is not None:
            # Only pair varies → {pair_key: DataArray}
            return {k[0]: v for k, v in results.items()}
        return results

    @_lineage_method("refstack")
    def get_ref(self, pair=None, components=None):
        """Load reference stacks. Requires 'refstack' in lineage.

        :returns: :class:`xarray.DataArray` (single pair/comp) or dict of
            DataArrays keyed by ``(pair, comp)``.
        """
        from .core.io import xr_get_ref
        from .core.workflow import get_t_axis

        lineage = self._lineage_through("refstack")
        taxis = get_t_axis(self.params)
        root = self.output_folder

        if pair is not None and components is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_ref(root, lineage, sta1, sta2, components, taxis)

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
                    results[(f"{sta1}:{sta2}", comp)] = xr_get_ref(
                        root, lineage, sta1, sta2, comp, taxis)
                except Exception:
                    pass
        if pair is not None:
            return {k[1]: v for k, v in results.items()}
        return results

    @_lineage_method("mwcs")
    def get_mwcs(self, pair=None, components=None, mov_stack=None):
        """Load MWCS results. Requires 'mwcs' in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets.
        """
        from .core.io import xr_get_mwcs
        lineage = self._lineage_through("mwcs")
        root = self.output_folder
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_mwcs(root, lineage, sta1, sta2, components, mov_stack)
        return self._load_pair_comp_movstack(root, lineage, xr_get_mwcs,
                                              pair, components, mov_stack)

    @_lineage_method("mwcs_dtt")
    def get_mwcs_dtt(self, pair=None, components=None, mov_stack=None):
        """Load MWCS-DTT results. Requires 'mwcs_dtt' in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets.
        """
        from .core.io import xr_get_dtt
        lineage = self._lineage_through("mwcs_dtt")
        root = self.output_folder
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_dtt(root, lineage, sta1, sta2, components, mov_stack)
        return self._load_pair_comp_movstack(root, lineage, xr_get_dtt,
                                              pair, components, mov_stack)

    @_lineage_method("stretching")
    def get_stretching(self, pair=None, components=None, mov_stack=None):
        """Load stretching results. Requires 'stretching' in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets.
        """
        from .core.io import _xr_get_stretching
        lineage = self._lineage_through("stretching")
        root = self.output_folder
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return _xr_get_stretching(root, lineage, sta1, sta2, components, mov_stack)
        return self._load_pair_comp_movstack(root, lineage, _xr_get_stretching,
                                              pair, components, mov_stack)

    @_lineage_method(["mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"])
    def get_dvv(self, pair_type="ALL", components=None, mov_stack=None):
        """Load DVV aggregate results. Requires a DVV step in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets keyed by
            ``(pair_type, comp, mov_stack)``.
        """
        from .core.io import xr_get_dvv_agg
        dvv_cat   = self.category
        lineage   = self._lineage_upstream_of(dvv_cat)
        step_name = self._step_name_for(dvv_cat)
        root      = self.output_folder
        if components is not None and mov_stack is not None:
            return xr_get_dvv_agg(root, lineage, step_name, mov_stack,
                                   pair_type, components)
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
                        root, lineage, step_name, ms_tuple, pt, comp)
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
    def get_wct(self, pair=None, components=None, mov_stack=None):
        """Load WCT results. Requires 'wavelet' in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets.
        """
        from .core.io import xr_load_wct
        lineage = self._lineage_through("wavelet")
        root = self.output_folder
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_load_wct(root, lineage, sta1, sta2, components, mov_stack)
        return self._load_pair_comp_movstack(root, lineage, xr_load_wct,
                                              pair, components, mov_stack)

    @_lineage_method("wavelet_dtt")
    def get_wct_dtt(self, pair=None, components=None, mov_stack=None):
        """Load WCT dt/t results. Requires 'wavelet_dtt' in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets.
        """
        from .core.io import xr_get_wct_dtt
        lineage = self._lineage_through("wavelet_dtt")
        root = self.output_folder
        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_wct_dtt(root, lineage, sta1, sta2, components, mov_stack)
        return self._load_pair_comp_movstack(root, lineage, xr_get_wct_dtt,
                                              pair, components, mov_stack)

    @_lineage_method("psd")
    def get_psd(self, seed_id=None, day=None):
        """Load PSD results. Requires 'psd' in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets keyed by
            ``(seed_id, day)``.
        """
        from .core.io import xr_load_psd
        step_name = self._step_name_for("psd")
        lineage   = self._lineage_upstream_of("psd")
        root      = self.output_folder
        if seed_id is not None and day is not None:
            return xr_load_psd(root, lineage, step_name, seed_id, day)
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
                r = xr_load_psd(root, lineage, step_name, sid, day_key)
                if r is not None:
                    results[(sid, day_key)] = r
        if seed_id is not None:
            return {k[1]: v for k, v in results.items()}
        return results

    @_lineage_method("psd_rms")
    def get_psd_rms(self, seed_id=None):
        """Load PSD RMS results. Requires 'psd_rms' in lineage.

        :returns: :class:`xarray.Dataset` or dict of Datasets keyed by
            seed_id.
        """
        from .core.io import xr_load_rms
        step_name = self._step_name_for("psd_rms")
        lineage   = self._lineage_upstream_of("psd_rms")
        root      = self.output_folder
        if seed_id is not None:
            return xr_load_rms(root, lineage, step_name, seed_id)
        base = os.path.join(root, *lineage, step_name, "_output")
        results = {}
        for sid_dir in self._discover(os.path.join(base, "*")):
            if not os.path.isdir(sid_dir):
                continue
            sid = os.path.basename(sid_dir)
            r = xr_load_rms(root, lineage, step_name, sid)
            if r is not None:
                results[sid] = r
        return results

    # ── helpers ───────────────────────────────────────────────────────────────

    def export_bundle(
        self,
        dest: str,
        from_step: str = None,
        compress: bool = False,
        overwrite: bool = False,
    ) -> str:
        """Export a portable, self-describing bundle of computed results.

        Copies the selected portion of the ``_output/`` tree together with
        ``params.yaml`` and ``MANIFEST.json`` (sha256 per file) into *dest*.
        The result can be read on any machine with MSNoise installed using
        :meth:`from_bundle` — no database connection required.

        The full lineage directory nesting is preserved verbatim, so that
        :meth:`from_bundle` can set ``output_folder = bundle_root`` and all
        existing :func:`~msnoise.core.io.xr_get_*` calls resolve paths
        identically.

        Parameters
        ----------
        dest:
            Output directory path.  Created if absent.  If *compress* is
            ``True``, a ``<dest>.zip`` file is written instead and the
            temporary directory is removed on success.
        from_step:
            Category name of the first step whose ``_output/`` to include.
            Must be present in the lineage.  All downstream steps are
            included automatically; ancestor directory names are created as
            empty path components to preserve the full nesting.

            Accepted values (examples): ``"stack"``, ``"refstack"``,
            ``"mwcs_dtt_dvv"``, ``"psd"``.  ``"cc"`` advances
            automatically to ``"filter"`` because raw CCFs are physically
            stored under the filter-step directory.

            ``None`` (default) picks the earliest step in the lineage that
            has an ``_output`` directory (``"stack"`` for CC workflows,
            ``"psd"`` for PSD-only workflows).
        compress:
            If ``True``, zip the bundle directory after writing and return
            the ``.zip`` path.  The directory is removed on success.
        overwrite:
            If ``True``, overwrite an existing *dest* directory or ``.zip``.

        Returns
        -------
        str
            Absolute path to the bundle directory or ``.zip`` file written.

        Raises
        ------
        ValueError
            If *from_step* is not found in the lineage.
        FileExistsError
            If *dest* (or ``<dest>.zip``) already exists and *overwrite*
            is ``False``.
        FileNotFoundError
            If the source directory for *from_step* does not exist on disk.

        Examples
        --------
        Export from refstack downwards, compressed::

            db = connect()
            r = MSNoiseResult.list(db, "mwcs_dtt_dvv")[0]
            path = r.export_bundle(
                "belgium_2023/",
                from_step="refstack",
                compress=True,
            )
            # → belgium_2023.zip

        Read on another machine::

            r2 = MSNoiseResult.from_bundle("belgium_2023.zip")
            r2.verify()
            ds = r2.get_dvv("CC", "ZZ", ("1D", "1D"))

        Journal supplementary data (dvv only)::

            r.export_bundle(
                "paper_SI/",
                from_step="mwcs_dtt_dvv",
                compress=True,
            )
        """
        import datetime
        import hashlib
        import json
        import shutil
        import zipfile

        try:
            import msnoise as _ms
            _version = _ms.__version__
        except Exception:
            _version = "unknown"

        dest = os.path.abspath(dest)
        zip_path = dest + ".zip" if compress else None

        # ── guard existing output ─────────────────────────────────────────
        if os.path.exists(dest) and not overwrite:
            raise FileExistsError(
                f"{dest!r} already exists.  Pass overwrite=True to replace it."
            )
        if zip_path and os.path.exists(zip_path) and not overwrite:
            raise FileExistsError(
                f"{zip_path!r} already exists.  Pass overwrite=True to replace it."
            )

        # ── resolve from_step ─────────────────────────────────────────────
        # Steps that write _output directly under their own step folder.
        # "cc" is special: its raw CCFs land under the filter-step folder,
        # so we advance from_step to "filter" to include that _output.
        _steps_with_output = frozenset({
            "filter", "stack", "refstack",
            "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
            "stretching", "stretching_dvv",
            "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
            "psd", "psd_rms",
        })

        if from_step is None:
            for name in self.lineage_names:
                cat = _step_prefix(name)
                if cat in _steps_with_output:
                    from_step = cat
                    break

        if from_step is None:
            raise ValueError(
                "Could not determine from_step — lineage has no exportable steps."
            )

        # "cc" advances to "filter" so that raw CCFs (stored under
        # filter_step/_output/all|daily/) are included.
        effective_from = "filter" if from_step == "cc" else from_step

        # Validate from_step exists in lineage
        found = any(
            _step_prefix(n) in (from_step, effective_from)
            for n in self.lineage_names
        )
        if not found:
            raise ValueError(
                f"from_step={from_step!r} not found in lineage "
                f"{'/'.join(self.lineage_names)!r}."
            )

        # ── determine source directory ────────────────────────────────────
        # lineage_through(effective_from) gives the path segment list up to
        # and including the from-step.  Everything at and below that directory
        # is the data we copy — child step folders live inside it.
        lineage_to_copy = self._lineage_through(effective_from)
        src_dir = os.path.join(self.output_folder, *lineage_to_copy)
        if not os.path.isdir(src_dir):
            raise FileNotFoundError(
                f"Source directory not found: {src_dir!r}.  "
                f"Has the pipeline been run up to '{from_step}'?"
            )

        # ── copy _output tree ─────────────────────────────────────────────
        if os.path.exists(dest) and overwrite:
            shutil.rmtree(dest)
        os.makedirs(dest)

        # Mirror full lineage nesting so output_folder = bundle_root works.
        dest_step_dir = os.path.join(dest, *lineage_to_copy)
        shutil.copytree(src_dir, dest_step_dir)

        # ── params.yaml ───────────────────────────────────────────────────
        self.params.to_yaml(os.path.join(dest, "params.yaml"))

        # ── provenance: stations + data_sources ───────────────────────────
        from .core.stations import get_station_pairs
        _station_ids: set = set()
        if self._db is not None:
            try:
                for _sta1, _sta2 in get_station_pairs(self._db):
                    _station_ids.add(f"{_sta1.net}.{_sta1.sta}")
                    _station_ids.add(f"{_sta2.net}.{_sta2.sta}")
            except Exception:
                pass
        data_sources_yaml, stations_yaml = _build_datasource_provenance(
            self._db, sorted(_station_ids)
        )

        import yaml as _yaml

        # ── MANIFEST.json ─────────────────────────────────────────────────
        generated = (
            datetime.datetime.now(datetime.timezone.utc)
            .isoformat(timespec="seconds") + "Z"
        )
        file_entries: list = []
        total_bytes = 0

        for dirpath, _dirs, filenames in os.walk(dest_step_dir):
            for fname in sorted(filenames):
                fpath = os.path.join(dirpath, fname)
                rel = os.path.relpath(fpath, dest).replace(os.sep, "/")
                size = os.path.getsize(fpath)
                total_bytes += size
                h = hashlib.sha256()
                with open(fpath, "rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
                file_entries.append({
                    "path":       rel,
                    "sha256":     h.hexdigest(),
                    "size_bytes": size,
                })

        manifest = {
            "msnoise_bundle_version":  1,
            "msnoise_version":         _version,
            "generated":               generated,
            "lineage":                 "/".join(self.lineage_names),
            "export_from_step":        from_step,
            "output_folder_original":  self.output_folder,
            "stations":    _yaml.safe_load(stations_yaml) if stations_yaml else [],
            "data_sources":_yaml.safe_load(data_sources_yaml) if data_sources_yaml else [],
            "files":           file_entries,
            "total_size_bytes": total_bytes,
            "total_files":      len(file_entries),
        }
        with open(os.path.join(dest, "MANIFEST.json"), "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2)

        # ── optional zip ──────────────────────────────────────────────────
        if compress:
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
                for dirpath, _dirs, filenames in os.walk(dest):
                    for fname in filenames:
                        fpath = os.path.join(dirpath, fname)
                        # arcname: bundle dir name + relative path inside
                        arcname = os.path.relpath(fpath, os.path.dirname(dest))
                        zf.write(fpath, arcname)
            shutil.rmtree(dest)
            return zip_path

        return dest

    @staticmethod
    def to_dataframe(ds):
        """Convert an xarray Dataset or DataArray returned by any ``get_*``
        method to a :class:`~pandas.DataFrame`.

        This is the recommended escape hatch for users who need pandas
        for custom analysis, plotting, or export.  All ``get_*`` methods
        return xarray objects; call this helper when you need a DataFrame::

            r = MSNoiseResult.from_ids(db, mwcs=1, mwcs_dtt=1)
            ds = r.get_mwcs_dtt("BE.UCC:BE.MEM", "ZZ", ("1D", "1D"))
            df = MSNoiseResult.to_dataframe(ds)

        For Datasets with multiple variables the result is a DataFrame with a
        :class:`~pandas.MultiIndex` on the columns.  For DataArrays the result
        is a flat DataFrame indexed by ``times``.

        :param ds: :class:`xarray.Dataset` or :class:`xarray.DataArray`.
        :returns: :class:`~pandas.DataFrame`.
        """
        import xarray as xr
        if isinstance(ds, xr.DataArray):
            return ds.to_dataframe(name=ds.name or "value")
        return ds.to_dataframe()

    # ── export ────────────────────────────────────────────────────────────────

    @_lineage_method(["mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"])
    def export_dvv(
        self,
        path: str,
        pair_type: str = "CC",
        components: str = None,
        mov_stack: tuple = None,
    ) -> list:
        """Export dv/v time series bundled with full parameter provenance.

        Writes one NetCDF file per ``(components, mov_stack)`` combination.
        Each file embeds the dv/v statistics (``mean``, ``std``, ``median``,
        ``n_pairs``, and any weighted / trimmed variants) plus global
        attributes for full reproducibility:

        ``lineage``
            Slash-separated step-name path.
        ``msnoise_params``
            Full YAML dump of :attr:`params` — one block per config category.
            Load offline with ``yaml.safe_load(ds.attrs["msnoise_params"])``.
        ``msnoise_version``
            Package version string.
        ``generated``
            ISO-8601 UTC timestamp.
        ``pair_type``, ``components``, ``mov_stack``
            Provenance of the specific slice exported.

        Parameters
        ----------
        path:
            Output directory (created if absent) **or** a full ``.nc`` path.
            When a directory is given, files are named
            ``dvv_<pair_type>_<comp>__<lineage_tag>__m<ms>.nc``.
            A full path requires exactly one ``(components, mov_stack)`` to be
            specified.
        pair_type:
            Pair-type filter (default ``"CC"``).
        components:
            Component string, e.g. ``"ZZ"``.  ``None`` exports all found.
        mov_stack:
            Tuple ``(window, step)``.  ``None`` exports all found.

        Returns
        -------
        list of str
            Paths of files written.

        Example::

            db = connect()
            r = MSNoiseResult.list(db, "mwcs_dtt_dvv")[0]
            written = r.export_dvv("exports/", pair_type="CC")
            for f in written:
                print(f)
                # dvv_CC_ZZ__pre1-cc1-f1-stk1-ref1-mwcs1-dtt1-dvv1__m1D-1D.nc

            # Reload and inspect provenance
            import xarray as xr, yaml
            ds = xr.open_dataset(written[0])
            params = yaml.safe_load(ds.attrs["msnoise_params"])
            print(params["mwcs"]["mwcs_wlen"])
        """
        import datetime
        try:
            import msnoise as _ms
            _version = _ms.__version__
        except Exception:
            _version = "unknown"

        from .core.config import lineage_to_plot_tag

        # Collect datasets via get_dvv (dict keyed by (pt, comp, ms))
        raw = self.get_dvv(pair_type=pair_type, components=components,
                           mov_stack=mov_stack)

        # Normalise to {(pt, comp, ms_tuple): ds}
        datasets = {}
        if isinstance(raw, dict):
            for key, ds in raw.items():
                if len(key) == 3:
                    datasets[key] = ds
                elif len(key) == 2:
                    # get_dvv collapsed one dimension
                    if mov_stack is not None:
                        pt, comp_k = key
                        datasets[(pt, comp_k, mov_stack)] = ds
                    else:
                        pt, ms_k = key
                        datasets[(pt, pair_type, ms_k)] = ds
        elif hasattr(raw, "times"):
            # Single xr.Dataset returned when both components and mov_stack given
            datasets[(pair_type,
                      components or "unknown",
                      mov_stack or ("?", "?"))] = raw

        if not datasets:
            raise FileNotFoundError(
                f"No dv/v data found for pair_type={pair_type!r}, "
                f"components={components!r}, mov_stack={mov_stack!r}."
            )

        lin_str     = "/".join(self.lineage_names)
        lin_tag     = lineage_to_plot_tag(self.lineage_names)
        params_yaml = self.params.to_yaml_string()
        generated   = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds") + "Z"

        # Collect station IDs from all pairs active in this project.
        # We use all station pairs rather than trying to parse pair names
        # from DVV output coordinates — this is more robust and captures
        # the full set of stations that contributed to the result.
        from .core.stations import get_station_pairs
        _station_ids = set()
        try:
            for _sta1, _sta2 in get_station_pairs(self._db):
                _station_ids.add(f"{_sta1.net}.{_sta1.sta}")
                _station_ids.add(f"{_sta2.net}.{_sta2.sta}")
        except Exception:
            pass  # non-fatal — provenance blocks will be empty if this fails
        data_sources_yaml, stations_yaml = _build_datasource_provenance(
            self._db, sorted(_station_ids)
        )

        is_dir = os.path.isdir(path) or not path.endswith(".nc")
        written = []

        for (pt, comp, ms), ds in sorted(datasets.items()):
            ms_str = (f"{ms[0]}_{ms[1]}" if isinstance(ms, (tuple, list))
                      else str(ms))

            if is_dir:
                os.makedirs(path, exist_ok=True)
                fname    = f"dvv_{pt}_{comp}__{lin_tag}__m{ms_str}.nc"
                out_path = os.path.join(path, fname)
            else:
                if len(datasets) > 1:
                    raise ValueError(
                        "A full file path can only be used when exactly one "
                        "(components, mov_stack) combination is selected."
                    )
                out_path = path

            ds = ds.copy()
            ds.attrs.update({
                "lineage":          lin_str,
                "msnoise_params":   params_yaml,
                "msnoise_version":  _version,
                "generated":        generated,
                "pair_type":        pt,
                "components":       comp,
                "mov_stack":        ms_str,
                "data_sources":     data_sources_yaml,
                "stations":         stations_yaml,
            })
            from .core.io import _f32_encoding
            ds.to_netcdf(out_path, encoding=_f32_encoding(ds))
            written.append(out_path)

        return written



# ── module-level helper ────────────────────────────────────────────────────────


def _build_datasource_provenance(db, station_ids: list) -> tuple:
    """Build data_sources and stations provenance dicts for export.

    :param db: SQLAlchemy session.
    :param station_ids: List of ``"NET.STA"`` strings involved in the result.
    :returns: Tuple ``(data_sources_yaml, stations_yaml)`` — YAML strings.
    """
    import yaml
    from .core.stations import resolve_data_source, get_station
    # Collect unique stations
    stations_data = []
    sources_seen  = {}

    for sid in sorted(set(station_ids)):
        parts = sid.split(".")
        net, sta = parts[0], parts[1]
        station = get_station(db, net, sta)
        if station is None:
            continue
        ds = resolve_data_source(db, station)

        if ds.ref not in sources_seen:
            sources_seen[ds.ref] = {
                "id":             ds.ref,
                "name":           ds.name,
                "uri":            ds.uri,
                "data_structure": ds.data_structure,
                "auth_env":       ds.auth_env,
                # credentials are never exported
            }

        stations_data.append({
            "net":            station.net,
            "sta":            station.sta,
            "X":              station.X,
            "Y":              station.Y,
            "altitude":       station.altitude,
            "coordinates":    station.coordinates,
            "data_source_id": ds.ref,
        })

    data_sources_yaml = yaml.dump(
        list(sources_seen.values()),
        default_flow_style=False, sort_keys=False, allow_unicode=True
    )
    stations_yaml = yaml.dump(
        stations_data,
        default_flow_style=False, sort_keys=False, allow_unicode=True
    )
    return data_sources_yaml, stations_yaml


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
