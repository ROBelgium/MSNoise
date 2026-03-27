"""
MSNoiseResult — user-facing class for loading computed results.

Typical notebook usage::

    from msnoise.results import MSNoiseResult
    from msnoise.api import connect

    db = connect()

    # Build from integer step IDs (most common)
    result = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                    stack=1, refstack=1)

    # Load a CCF for a specific pair / component / mov_stack
    df = result.get_ccf("YA.UV05.00:YA.UV06.00", "ZZ", ("6h", "6h"))

    # Omit any argument to get back a dict of all matches
    all_ccfs = result.get_ccf()   # {(pair, comp, mov_stack): DataFrame}

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

# WORKFLOW_ORDER defines the canonical step-name prefixes and their processing order.
_STEP_PREFIXES = [
    "preprocess",
    "cc",
    "psd",
    "psd_rms",
    "filter",
    "stack",
    "refstack",
    "mwcs",
    "mwcs_dtt",
    "mwcs_dtt_dvv",
    "stretching",
    "stretching_dvv",
    "wavelet",
    "wavelet_dtt",
    "wavelet_dtt_dvv",
]

# Map step-name prefix → keyword argument name used in from_ids()
_PREFIX_TO_KWARG = {p: p for p in _STEP_PREFIXES}


def _step_prefix(step_name: str) -> str:
    """Return the category prefix of a step name, e.g. 'mwcs_dtt_1' → 'mwcs_dtt'."""
    for prefix in sorted(_STEP_PREFIXES, key=len, reverse=True):
        if step_name.startswith(prefix + "_"):
            return prefix
    raise ValueError(f"Cannot determine category prefix for step name {step_name!r}")


class MSNoiseResult:
    """A resolved pipeline branch with convenience data-loading methods.

    Do not instantiate directly — use the class-method constructors
    :meth:`from_ids`, :meth:`from_names`, or :meth:`list`.

    Attributes
    ----------
    lineage_names : list[str]
        Ordered step-name strings from root to terminal step,
        e.g. ``['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']``.
    category : str
        Category of the terminal step, e.g. ``'refstack'``.
    params : AttribDict
        Merged configuration parameters for this lineage.
    output_folder : str
        Root output folder (``params.output_folder``).
    """

    # ------------------------------------------------------------------ #
    # Construction                                                         #
    # ------------------------------------------------------------------ #

    def __init__(self, db, lineage_names: list[str]):
        from .api import resolve_lineage_params
        self._db = db
        self.lineage_names: list[str] = list(lineage_names)
        _, _, self.params = resolve_lineage_params(db, lineage_names)
        self.output_folder: str = getattr(self.params, "output_folder", "OUTPUT")
        self.category: str = _step_prefix(lineage_names[-1]) if lineage_names else ""

    @classmethod
    def from_names(cls, db, names: list[str]) -> "MSNoiseResult":
        """Build from an explicit list of step-name strings.

        Example::

            r = MSNoiseResult.from_names(
                db, ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']
            )
        """
        return cls(db, names)

    @classmethod
    def from_ids(
        cls,
        db,
        preprocess: Optional[int] = None,
        cc: Optional[int] = None,
        psd: Optional[int] = None,
        psd_rms: Optional[int] = None,
        filter: Optional[int] = None,
        stack: Optional[int] = None,
        refstack: Optional[int] = None,
        mwcs: Optional[int] = None,
        mwcs_dtt: Optional[int] = None,
        mwcs_dtt_dvv: Optional[int] = None,
        stretching: Optional[int] = None,
        stretching_dvv: Optional[int] = None,
        wavelet: Optional[int] = None,
        wavelet_dtt: Optional[int] = None,
        wavelet_dtt_dvv: Optional[int] = None,
    ) -> "MSNoiseResult":
        """Build from integer configset IDs, stopping at the last non-None step.

        Only the steps you specify are included in the lineage.  Steps are
        added in canonical workflow order regardless of the keyword order you
        use.

        Example::

            # CC branch up to refstack
            r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                       stack=1, refstack=1)

            # PSD branch
            r = MSNoiseResult.from_ids(db, psd=1, psd_rms=1)

            # MWCS-DTT branch
            r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                       stack=1, refstack=1, mwcs=1, mwcs_dtt=1)

            # MWCS DVV aggregate
            r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                       stack=1, refstack=1, mwcs=1,
                                       mwcs_dtt=1, mwcs_dtt_dvv=1)
        """
        kwargs = {
            "preprocess":    preprocess,
            "cc":            cc,
            "psd":           psd,
            "psd_rms":       psd_rms,
            "filter":        filter,
            "stack":         stack,
            "refstack":      refstack,
            "mwcs":          mwcs,
            "mwcs_dtt":      mwcs_dtt,
            "mwcs_dtt_dvv":  mwcs_dtt_dvv,
            "stretching":    stretching,
            "stretching_dvv":stretching_dvv,
            "wavelet":       wavelet,
            "wavelet_dtt":   wavelet_dtt,
            "wavelet_dtt_dvv":   wavelet_dtt_dvv,
        }
        parts = []
        for prefix in _STEP_PREFIXES:
            val = kwargs.get(prefix)
            if val is not None:
                parts.append(f"{prefix}_{val}")
        if not parts:
            raise ValueError("At least one step ID must be provided.")
        return cls(db, parts)

    @classmethod
    def list(
        cls,
        db,
        category: str,
        include_empty: bool = False,
    ) -> list["MSNoiseResult"]:
        """Return all done :class:`MSNoiseResult` objects for a step category.

        Wraps :func:`~msnoise.api.get_done_lineages_for_category`.

        Parameters
        ----------
        category:
            Step category to query, e.g. ``'mwcs_dtt'``, ``'stretching'``.
        include_empty:
            If True, also include lineages that have no done jobs (useful for
            inspecting configured-but-unrun branches).  Default False.

        Example::

            for r in MSNoiseResult.list(db, 'mwcs_dtt'):
                print(r.lineage_names)
        """
        from .api import get_done_lineages_for_category, get_workflow_steps

        results = []

        if include_empty:
            # Return all active workflow steps of this category,
            # enumerating ALL distinct upstream paths (handles multi-configset DAGs).
            for step in get_workflow_steps(db):
                if step.category == category:
                    for names in _upstream_lineage_for_step(db, step):
                        if names:
                            results.append(cls(db, names))
        else:
            for names in get_done_lineages_for_category(db, category):
                results.append(cls(db, names))

        return results

    # ------------------------------------------------------------------ #
    # Navigation                                                           #
    # ------------------------------------------------------------------ #

    def branches(self, include_empty: bool = False) -> list["MSNoiseResult"]:
        """Return downstream :class:`MSNoiseResult` objects one step below.

        Follows active :class:`~msnoise.msnoise_table_def.WorkflowLink` rows
        from the terminal step of this lineage.  By default only returns
        branches that have at least one Done job.

        Parameters
        ----------
        include_empty:
            If True, return branches even if no jobs are done yet.

        Example::

            refstack_result = MSNoiseResult.from_ids(db, preprocess=1, cc=1,
                                                      filter=1, stack=1, refstack=1)
            for branch in refstack_result.branches():
                # branch.category is 'mwcs', 'stretching', or 'wavelet'
                print(branch.category, branch.lineage_names[-1])
        """
        from .msnoise_table_def import declare_tables

        schema = declare_tables()
        Job = schema.Job
        WorkflowStep = schema.WorkflowStep
        WorkflowLink = schema.WorkflowLink

        # Find the WorkflowStep for the terminal step of our lineage
        terminal_name = self.lineage_names[-1]
        terminal_step = (
            self._db.query(WorkflowStep)
            .filter(WorkflowStep.step_name == terminal_name)
            .first()
        )
        if terminal_step is None:
            return []

        # Follow outgoing active links
        links = (
            self._db.query(WorkflowLink)
            .filter(WorkflowLink.from_step_id == terminal_step.step_id)
            .filter(WorkflowLink.is_active == True)
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

            # Build child lineage = current lineage + child step name
            child_names = self.lineage_names + [child_step.step_name]
            child_lineage_str = "/".join(child_names)

            if not include_empty:
                # Check at least one Done job exists for this lineage+step
                has_done = (
                    self._db.query(Job.ref)
                    .filter(Job.step_id == child_step.step_id)
                    .filter(Job.lineage == child_lineage_str)
                    .filter(Job.flag == "D")
                    .first()
                )
                if has_done is None:
                    continue

            results.append(MSNoiseResult(self._db, child_names))

        return results

    def __repr__(self) -> str:
        return (
            f"MSNoiseResult(category={self.category!r}, "
            f"lineage={'/'.join(self.lineage_names)!r})"
        )

    # ------------------------------------------------------------------ #
    # Internal helpers                                                     #
    # ------------------------------------------------------------------ #

    def _require_category(self, *categories: str) -> None:
        """Raise ValueError if the lineage doesn't reach any of the given categories."""
        present = {_step_prefix(n) for n in self.lineage_names}
        for cat in categories:
            if cat not in present:
                tip = "call .branches() to reach downstream steps"
                raise ValueError(
                    f"Lineage {'/'.join(self.lineage_names)!r} does not include "
                    f"a '{cat}' step. {tip}."
                )

    def _lineage_upstream_of(self, category: str) -> list[str]:
        """Return lineage_names up to but not including the step of *category*."""
        result = []
        for name in self.lineage_names:
            if _step_prefix(name) == category:
                break
            result.append(name)
        return result

    def _lineage_through(self, category: str) -> list[str]:
        """Return lineage_names up to and including the step of *category*."""
        result = []
        for name in self.lineage_names:
            result.append(name)
            if _step_prefix(name) == category:
                break
        return result

    def _step_name_for(self, category: str) -> str:
        """Return the step name for a given category in this lineage."""
        for name in self.lineage_names:
            if _step_prefix(name) == category:
                return name
        raise ValueError(f"No '{category}' step in lineage.")

    def _discover(self, path_pattern: str) -> list[str]:
        """Glob *path_pattern* and return sorted unique matches."""
        return sorted(set(glob.glob(path_pattern)))

    # ------------------------------------------------------------------ #
    # Load methods                                                         #
    # ------------------------------------------------------------------ #

    def get_ccf(
        self,
        pair: Optional[str] = None,
        components: Optional[str] = None,
        mov_stack: Optional[tuple] = None,
        format: str = "xarray",
    ):
        """Load cross-correlation functions from the stack step output.

        Parameters
        ----------
        pair:
            Station pair string ``"NET.STA.LOC:NET.STA.LOC"``.
            If None, all pairs are returned.
        components:
            Component string e.g. ``"ZZ"``.  If None, all components returned.
        mov_stack:
            Moving-stack tuple e.g. ``("6h", "6h")``.  If None, all stacks returned.
        format:
            ``"xarray"`` (default) or ``"dataframe"``.

        Returns
        -------
        If all args are specified: a :class:`~pandas.DataFrame`.
        Otherwise: a dict keyed by ``(pair, components, mov_stack)``.
        """
        from .api import xr_get_ccf, get_t_axis

        self._require_category("stack")
        # CCFs live under the stack step folder
        lineage = self._lineage_through("stack")
        taxis = get_t_axis(self.params)
        root = self.output_folder

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_ccf(root, lineage, sta1, sta2, components,
                               mov_stack, taxis, format=format)

        # Discovery mode — glob the output tree
        base = os.path.join(root, *lineage, "_output")
        results = {}

        # Resolve wildcards
        mov_stacks = (
            ["%s_%s" % (mov_stack[0], mov_stack[1])] if mov_stack is not None
            else [os.path.basename(p)
                  for p in self._discover(os.path.join(base, "*"))
                  if os.path.isdir(p)]
        )
        comps = (
            [components] if components is not None
            else None  # resolved per mov_stack below
        )

        for ms_str in mov_stacks:
            ms_tuple = tuple(ms_str.split("_", 1))
            comp_list = (
                comps if comps is not None
                else [os.path.basename(p)
                      for p in self._discover(os.path.join(base, ms_str, "*"))
                      if os.path.isdir(p)]
            )
            for comp in comp_list:
                nc_dir = os.path.join(base, ms_str, comp)
                if pair is not None:
                    sta1, sta2 = pair.split(":")
                    files = self._discover(
                        os.path.join(nc_dir, f"{sta1}_{sta2}.nc"))
                else:
                    files = self._discover(os.path.join(nc_dir, "*.nc"))

                for fpath in files:
                    fname = os.path.splitext(os.path.basename(fpath))[0]
                    sta1, sta2 = fname.split("_", 1)
                    pair_key = f"{sta1}:{sta2}"
                    try:
                        df = xr_get_ccf(root, lineage, sta1, sta2, comp,
                                        ms_tuple, taxis, format=format)
                        results[(pair_key, comp, ms_tuple)] = df
                    except (FileNotFoundError, Exception):
                        pass

        if pair is not None and components is not None and mov_stack is None:
            # Simplify key to just mov_stack
            return {k[2]: v for k, v in results.items()}
        return results

    def get_ref(
        self,
        pair: Optional[str] = None,
        components: Optional[str] = None,
        format: str = "xarray",
    ):
        """Load reference stack from the refstack step output.

        Parameters
        ----------
        format:
            ``"xarray"`` (default) returns an :class:`xarray.Dataset`.
            ``"dataframe"`` calls ``.CCF.to_dataframe()`` on the Dataset.

        Returns
        -------
        If all args specified: result in requested format.
        Otherwise: dict keyed by ``(pair, components)``.
        """
        from .api import xr_get_ref, get_t_axis

        self._require_category("refstack")
        lineage = self._lineage_through("refstack")
        taxis = get_t_axis(self.params)
        root = self.output_folder

        def _apply_format(ds):
            if format == "dataframe":
                import pandas as pd_ref
                da = ds.REF
                return pd_ref.Series(
                    da.values,
                    index=da.coords["taxis"].values,
                    name="REF",
                )
            return ds

        if pair is not None and components is not None:
            sta1, sta2 = pair.split(":")
            return _apply_format(xr_get_ref(root, lineage, sta1, sta2, components, taxis))

        base = os.path.join(root, *lineage, "_output", "REF")
        results = {}
        comp_list = ([components] if components is not None
                     else [os.path.basename(p)
                            for p in self._discover(os.path.join(base, "*"))
                            if os.path.isdir(p)])

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
                pair_key = f"{sta1}:{sta2}"
                try:
                    ds = xr_get_ref(root, lineage, sta1, sta2, comp, taxis)
                    results[(pair_key, comp)] = _apply_format(ds)
                except (FileNotFoundError, Exception):
                    pass

        if pair is not None:
            return {k[1]: v for k, v in results.items()}
        return results

    def get_mwcs(
        self,
        pair: Optional[str] = None,
        components: Optional[str] = None,
        mov_stack: Optional[tuple] = None,
        format: str = "xarray",
    ):
        """Load MWCS results.

        Parameters
        ----------
        format:
            ``"xarray"`` (default) or ``"dataframe"``.

        Returns
        -------
        If all args specified: result in requested format.
        Otherwise: dict keyed by ``(pair, components, mov_stack)``.
        """
        from .api import xr_get_mwcs

        self._require_category("mwcs")
        lineage = self._lineage_through("mwcs")
        root = self.output_folder
        api_fmt = "dataset" if format == "xarray" else "dataframe"

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_mwcs(root, lineage, sta1, sta2, components, mov_stack,
                               format=api_fmt)

        return self._load_pair_comp_movstack(
            root, lineage,
            lambda *a, **kw: xr_get_mwcs(*a, format=api_fmt, **kw),
            pair, components, mov_stack)

    def get_mwcs_dtt(
        self,
        pair: Optional[str] = None,
        components: Optional[str] = None,
        mov_stack: Optional[tuple] = None,
        format: str = "xarray",
    ):
        """Load MWCS-DTT (dt/t) results.

        Parameters
        ----------
        format:
            ``"xarray"`` (default) or ``"dataframe"``.

        Returns
        -------
        If all args specified: result in requested format.
        Otherwise: dict keyed by ``(pair, components, mov_stack)``.
        """
        from .api import xr_get_dtt

        self._require_category("mwcs_dtt")
        lineage = self._lineage_through("mwcs_dtt")
        root = self.output_folder
        api_fmt = "dataset" if format == "xarray" else "dataframe"

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return xr_get_dtt(root, lineage, sta1, sta2, components, mov_stack,
                              format=api_fmt)

        return self._load_pair_comp_movstack(
            root, lineage,
            lambda *a, **kw: xr_get_dtt(*a, format=api_fmt, **kw),
            pair, components, mov_stack)

    def get_stretching(
        self,
        pair: Optional[str] = None,
        components: Optional[str] = None,
        mov_stack: Optional[tuple] = None,
        format: str = "xarray",
    ):
        """Load stretching results.

        Parameters
        ----------
        format:
            ``"xarray"`` (default) or ``"dataframe"``.

        Returns
        -------
        If all args specified: result in requested format.
        Otherwise: dict keyed by ``(pair, components, mov_stack)``.
        """
        from .api import _xr_get_stretching

        self._require_category("stretching")
        lineage = self._lineage_through("stretching")
        root = self.output_folder
        api_fmt = "dataset" if format == "xarray" else "dataframe"

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            return _xr_get_stretching(root, lineage, sta1, sta2,
                                      components, mov_stack, format=api_fmt)

        return self._load_pair_comp_movstack(
            root, lineage,
            lambda *a, **kw: _xr_get_stretching(*a, format=api_fmt, **kw),
            pair, components, mov_stack)

    def get_dvv(
        self,
        pair_type: str = "ALL",
        components: Optional[str] = None,
        mov_stack: Optional[tuple] = None,
        format: str = "xarray",
    ):
        """Load pre-aggregated network dv/v statistics.

        The :class:`MSNoiseResult` must be instantiated at a DVV aggregate step
        (``mwcs_dtt_dvv``, ``stretching_dvv``, or ``wavelet_dtt_dvv``).

        Parameters
        ----------
        pair_type:
            ``"ALL"`` (default), ``"CC"``, ``"SC"``, or ``"AC"``.
        components:
            Component string e.g. ``"ZZ"`` or ``"ALL"``.  If *None*, all
            available components are returned.
        mov_stack:
            Moving-stack tuple e.g. ``("1D", "1D")``.  If *None*, all
            available mov_stacks are returned.
        format:
            ``"xarray"`` (default) returns an :class:`xarray.Dataset`.
            ``"dataframe"`` returns a :class:`~pandas.DataFrame`.

        Returns
        -------
        If *components* and *mov_stack* are both specified: result in
        requested format.
        Otherwise: dict keyed by ``(pair_type, components, mov_stack)``.
        """
        from .api import xr_get_dvv_agg

        DVV_CATEGORIES = ("mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv")
        self._require_category(*DVV_CATEGORIES)

        # Resolve the DVV step name and its lineage
        dvv_cat = self.category
        lineage   = self._lineage_through(dvv_cat)
        step_name = self._step_name_for(dvv_cat)
        root      = self.output_folder
        api_fmt   = "dataset" if format == "xarray" else "dataframe"

        if components is not None and mov_stack is not None:
            return xr_get_dvv_agg(root, lineage, step_name,
                                   mov_stack, pair_type, components,
                                   format=api_fmt)

        # Discovery mode — find all available (pair_type, component, mov_stack)
        base = os.path.join(root, *lineage, step_name, "_output")
        results = {}

        ms_dirs = (
            [os.path.join(base, "%s_%s" % (mov_stack[0], mov_stack[1]))]
            if mov_stack is not None
            else [p for p in self._discover(os.path.join(base, "*"))
                  if os.path.isdir(p)]
        )
        for ms_dir in ms_dirs:
            ms_str   = os.path.basename(ms_dir)
            ms_tuple = tuple(ms_str.split("_", 1))
            pattern  = f"dvv_{pair_type}_*.nc" if components is None else \
                       f"dvv_{pair_type}_{components}.nc"
            for nc_file in self._discover(os.path.join(ms_dir, pattern)):
                fname = os.path.splitext(os.path.basename(nc_file))[0]
                # dvv_CC_ZZ  →  parts = ["dvv", "CC", "ZZ"]
                parts = fname.split("_", 2)
                if len(parts) < 3:
                    continue
                _, pt, comp = parts
                try:
                    results[(pt, comp, ms_tuple)] = xr_get_dvv_agg(
                        root, lineage, step_name,
                        ms_tuple, pt, comp, format=api_fmt)
                except (FileNotFoundError, Exception):
                    pass

        # Simplify keys when specific dimensions are constrained
        if components is not None and mov_stack is not None:
            return {k[0]: v for k, v in results.items()}   # keyed by pair_type
        if components is not None:
            return {(k[0], k[2]): v for k, v in results.items()}
        if mov_stack is not None:
            return {(k[0], k[1]): v for k, v in results.items()}
        return results

    def get_wct(
        self,
        pair: Optional[str] = None,
        components: Optional[str] = None,
        mov_stack: Optional[tuple] = None,
        format: str = "xarray",
    ):
        """Load Wavelet Coherence Transform results.

        Parameters
        ----------
        format:
            ``"xarray"`` (default) returns an :class:`xarray.Dataset`.
            ``"dataframe"`` flattens to a :class:`~pandas.DataFrame`.

        Returns
        -------
        If all args specified: result in requested format.
        Otherwise: dict keyed by ``(pair, components, mov_stack)``.
        """
        from .api import xr_load_wct

        self._require_category("wavelet")
        lineage = self._lineage_through("wavelet")
        root = self.output_folder

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            ds = xr_load_wct(root, lineage, sta1, sta2, components, mov_stack)
            if format == "dataframe":
                return ds.to_dataframe()
            return ds

        return self._load_pair_comp_movstack(
            root, lineage, xr_load_wct, pair, components, mov_stack)

    def get_wct_dtt(
        self,
        pair: Optional[str] = None,
        components: Optional[str] = None,
        mov_stack: Optional[tuple] = None,
        format: str = "xarray",
    ):
        """Load WCT dt/t results.

        Parameters
        ----------
        format:
            ``"xarray"`` (default) returns an :class:`xarray.Dataset`.
            ``"dataframe"`` flattens to a :class:`~pandas.DataFrame`.

        Returns
        -------
        If all args specified: result in requested format.
        Otherwise: dict keyed by ``(pair, components, mov_stack)``.
        """
        from .api import xr_get_wct_dtt

        self._require_category("wavelet_dtt")
        lineage = self._lineage_through("wavelet_dtt")
        root = self.output_folder

        if pair is not None and components is not None and mov_stack is not None:
            sta1, sta2 = pair.split(":")
            ds = xr_get_wct_dtt(root, lineage, sta1, sta2, components, mov_stack)
            if format == "dataframe":
                return ds.to_dataframe()
            return ds

        return self._load_pair_comp_movstack(
            root, lineage, xr_get_wct_dtt, pair, components, mov_stack)

    def get_psd(
        self,
        seed_id: Optional[str] = None,
        day: Optional[str] = None,
        format: str = "xarray",
    ):
        """Load daily PSD results.

        Parameters
        ----------
        seed_id:
            SEED identifier ``"NET.STA.LOC.CHAN"``.  If None, all channels.
        day:
            Date string ``"YYYY-MM-DD"``.  If None, all days.
        format:
            ``"xarray"`` (default) or ``"dataframe"``.

        Returns
        -------
        If all args specified: result in requested format.
        Otherwise: dict keyed by ``(seed_id, day)``.
        """
        from .api import xr_load_psd

        self._require_category("psd")
        step_name = self._step_name_for("psd")
        lineage = self._lineage_upstream_of("psd")
        root = self.output_folder
        api_fmt = "dataset" if format == "xarray" else "dataframe"

        if seed_id is not None and day is not None:
            return xr_load_psd(root, lineage, step_name, seed_id, day,
                               format=api_fmt)

        base = os.path.join(root, *lineage, step_name, "_output", "daily")
        results = {}
        seed_ids = (
            [seed_id] if seed_id is not None
            else [os.path.basename(p)
                  for p in self._discover(os.path.join(base, "*"))
                  if os.path.isdir(p)]
        )
        for sid in seed_ids:
            day_files = (
                self._discover(os.path.join(base, sid, f"{day}.nc"))
                if day is not None
                else self._discover(os.path.join(base, sid, "*.nc"))
            )
            for fpath in day_files:
                day_key = os.path.splitext(os.path.basename(fpath))[0]
                result = xr_load_psd(root, lineage, step_name, sid, day_key,
                                     format=api_fmt)
                if result is not None:
                    results[(sid, day_key)] = result

        if seed_id is not None:
            return {k[1]: v for k, v in results.items()}
        return results

    def get_psd_rms(
        self,
        seed_id: Optional[str] = None,
        format: str = "xarray",
    ):
        """Load PSD RMS results.

        Parameters
        ----------
        seed_id:
            SEED identifier ``"NET.STA.LOC.CHAN"``.  If None, all channels.
        format:
            ``"xarray"`` (default) or ``"dataframe"``.

        Returns
        -------
        If specified: result in requested format.
        Otherwise: dict keyed by ``seed_id``.
        """
        from .api import xr_load_rms

        self._require_category("psd_rms")
        step_name = self._step_name_for("psd_rms")
        lineage = self._lineage_upstream_of("psd_rms")
        root = self.output_folder
        api_fmt = "dataset" if format == "xarray" else "dataframe"

        if seed_id is not None:
            return xr_load_rms(root, lineage, step_name, seed_id, format=api_fmt)

        base = os.path.join(root, *lineage, step_name, "_output")
        results = {}
        for sid_dir in self._discover(os.path.join(base, "*")):
            if not os.path.isdir(sid_dir):
                continue
            sid = os.path.basename(sid_dir)
            result = xr_load_rms(root, lineage, step_name, sid, format=api_fmt)
            if result is not None:
                results[sid] = result
        return results

    # ------------------------------------------------------------------ #
    # Shared discovery helper                                              #
    # ------------------------------------------------------------------ #

    def _load_pair_comp_movstack(self, root, lineage, loader_fn,
                                  pair, components, mov_stack):
        """Generic discovery+load for functions with (root, lineage, sta1, sta2, comp, mov_stack) signature."""
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
                if pair is not None:
                    sta1, sta2 = pair.split(":")
                    files = self._discover(
                        os.path.join(nc_dir, f"{sta1}_{sta2}.nc"))
                else:
                    files = self._discover(os.path.join(nc_dir, "*.nc"))

                for fpath in files:
                    fname = os.path.splitext(os.path.basename(fpath))[0]
                    sta1, sta2 = fname.split("_", 1)
                    pair_key = f"{sta1}:{sta2}"
                    try:
                        data = loader_fn(root, lineage, sta1, sta2,
                                         comp, ms_tuple)
                        results[(pair_key, comp, ms_tuple)] = data
                    except (FileNotFoundError, Exception):
                        pass

        # Simplify keys when only one dimension is wildcarded
        if pair is not None and components is not None:
            return {k[2]: v for k, v in results.items()}
        if pair is not None and mov_stack is not None:
            return {k[1]: v for k, v in results.items()}
        if components is not None and mov_stack is not None:
            return {k[0]: v for k, v in results.items()}
        return results


# ------------------------------------------------------------------ #
# Module-level helper                                                  #
# ------------------------------------------------------------------ #

def _upstream_lineage_for_step(db, step) -> list[str]:
    """Walk WorkflowLinks upstream from *step* to build all full lineage name lists.

    Returns a list of lineage name lists (one per distinct upstream path).
    In a simple single-configset workflow this returns exactly one list.
    In a multi-configset DAG (multiple incoming links to one step), all
    paths are enumerated.
    """
    from .api import get_workflow_links
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    WorkflowStep = schema.WorkflowStep

    links = get_workflow_links(db)
    # Build multi-map: to_step_id → [from_step_id, ...]
    parent_map: dict[int, list[int]] = {}
    for lnk in links:
        parent_map.setdefault(lnk.to_step_id, []).append(lnk.from_step_id)

    def _walk(step_id, visited) -> list[list[str]]:
        """Recursively walk upstream, returning all root→current paths."""
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

    all_paths = _walk(step.step_id, set())
    # Filter out single-element paths (just the step itself, no real upstream)
    return [p for p in all_paths if p]
