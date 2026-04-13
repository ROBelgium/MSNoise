.. _msnoise_result:

**************
MSNoiseResult
**************

``MSNoiseResult`` is the **recommended way to read pipeline outputs** from a
notebook or script.  Rather than constructing file paths by hand or calling
low-level ``core.io`` functions directly, you build a result object that knows
which lineage branch it covers and exposes only the ``get_*`` methods that are
valid for that branch.

.. contents::
   :local:
   :depth: 2


Quickstart
==========

.. code-block:: python

    from msnoise.results import MSNoiseResult
    from msnoise.core.db import connect

    db = connect()

    # Build a result object for a full CCF → dv/v lineage
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                               stack=1, refstack=1,
                               mwcs=1, mwcs_dtt=1, mwcs_dtt_dvv=1)

    # Load a stacked CCF (xarray DataArray)
    da = r.get_ccf("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", ("1D", "1D"))

    # Load a dv/v curve (xarray Dataset)
    ds = r.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D", "1D"))

    # Convert to pandas when needed
    df = MSNoiseResult.to_dataframe(ds)

You do not need to know the on-disk path layout or the low-level I/O functions:
``MSNoiseResult`` handles both the path construction and the xarray loading.


Constructing a result object
=============================

Three constructors are available:

``from_ids`` (most common)
    Pass the integer config-set number for every step you need.  You only need
    to go as far down the pipeline as the data you want to read:

    .. code-block:: python

        # Read raw CC outputs only — no stacking needed
        r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)

        # Read stacked CCFs and the reference
        r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                   stack=1, refstack=1)

        # Full pipeline to dv/v via stretching
        r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                   stack=1, refstack=1,
                                   stretching=1, stretching_dvv=1)

        # PSD branch
        r = MSNoiseResult.from_ids(db, psd=1, psd_rms=1)

``from_names``
    Build from an explicit list of step-name strings (useful when scripting
    over results returned by :meth:`branches`):

    .. code-block:: python

        r = MSNoiseResult.from_names(
            db, ["preprocess_1", "cc_1", "filter_1", "stack_1"]
        )

``list``
    Iterate over every completed lineage for a given step category:

    .. code-block:: python

        for r in MSNoiseResult.list(db, "mwcs_dtt"):
            print(r)


Dynamic method gating
======================

Only ``get_*`` methods whose required step category is present in the lineage
are accessible.  Trying to call a gated method raises ``AttributeError`` with
a helpful message.  This also means Jupyter's tab-completion only shows
methods that can actually return data:

.. code-block:: python

    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)

    r.get_ccf_raw(...)  # ✓  'cc' is in lineage
    r.get_ccf(...)      # ✗  AttributeError: requires 'stack'
    r.get_mwcs(...)     # ✗  AttributeError: requires 'mwcs'

The table below maps each method to its minimum required lineage:

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Method
     - Requires
     - What it reads
   * - :meth:`get_ccf_raw <msnoise.results.MSNoiseResult.get_ccf_raw>`
     - ``cc``
     - Raw per-window or daily-stacked CCFs from ``s03``
   * - :meth:`get_ccf <msnoise.results.MSNoiseResult.get_ccf>`
     - ``stack``
     - Moving-stacked CCFs from ``s04``
   * - :meth:`get_ref <msnoise.results.MSNoiseResult.get_ref>`
     - ``refstack``
     - Reference stacks from ``s04_stack_refstack``
   * - :meth:`get_mwcs <msnoise.results.MSNoiseResult.get_mwcs>`
     - ``mwcs``
     - MWCS results from ``s05``
   * - :meth:`get_mwcs_dtt <msnoise.results.MSNoiseResult.get_mwcs_dtt>`
     - ``mwcs_dtt``
     - MWCS dt/t results from ``s06``
   * - :meth:`get_stretching <msnoise.results.MSNoiseResult.get_stretching>`
     - ``stretching``
     - Stretching results from ``s10``
   * - :meth:`get_wct <msnoise.results.MSNoiseResult.get_wct>`
     - ``wavelet``
     - Wavelet coherence results from ``s08``
   * - :meth:`get_wct_dtt <msnoise.results.MSNoiseResult.get_wct_dtt>`
     - ``wavelet_dtt``
     - WCT dt/t results from ``s09``
   * - :meth:`get_dvv <msnoise.results.MSNoiseResult.get_dvv>`
     - any ``*_dvv`` step
     - Aggregated dv/v from ``s07``
   * - :meth:`get_psd <msnoise.results.MSNoiseResult.get_psd>`
     - ``psd``
     - Daily PSD from ``psd_compute``
   * - :meth:`get_psd_rms <msnoise.results.MSNoiseResult.get_psd_rms>`
     - ``psd_rms``
     - PSD RMS from ``psd_compute_rms``


Reading raw CC outputs (pre-stack)
====================================

``get_ccf_raw`` reads the files written directly by ``s03_compute_no_rotation``
*before* any stacking.  It only requires ``cc`` in the lineage.

Two storage layouts exist depending on which config flags are active:

``kind="all"`` (default)
    Per-window CCFs stored under ``_output/all/``.  Each file covers one
    calendar day and has dimensions ``(times, taxis)`` where ``times``
    holds the window start times.

``kind="daily"``
    Daily-stacked CCFs stored under ``_output/daily/``.  Each file covers
    one day with dimension ``(taxis,)`` only.  When multiple days are
    returned the method expands each file with a ``times`` coordinate so
    they can be concatenated.

.. code-block:: python

    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)

    # Single day, per-window — DataArray with dims (times, taxis)
    da = r.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ",
                        date="2023-01-01", kind="all")

    # All days for one pair/comp, daily stacks — dict keyed by date string
    d = r.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", kind="daily")

    # Concatenate into a single (times, taxis) DataArray
    import xarray as xr
    da_all = xr.concat(d.values(), dim="times").sortby("times")

    # All pairs, all comps, single day — dict keyed by (pair, comp, date)
    everything = r.get_ccf_raw(date="2023-01-01", kind="all")

.. seealso::

    :doc:`workflow_ccf/005_compute_cc` — the step that writes these files.


Reading stacked CCFs
=====================

``get_ccf`` reads the moving-stacked CCFs written by ``s04_stack_mov`` and
requires ``stack`` in the lineage.  Results are keyed by
``(mov_stack_window, mov_stack_step)``:

.. code-block:: python

    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)

    # Single pair / component / moving-stack window
    da = r.get_ccf("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", ("1D", "1D"))

    # All pairs and components for a fixed moving-stack window
    d = r.get_ccf(mov_stack=("1D", "1D"))  # {(pair, comp, ms): DataArray}

.. seealso::

    :doc:`workflow_ccf/006_stack` — the step that writes these files.


Reading dv/v results
======================

``get_dvv`` is available when any of ``mwcs_dtt_dvv``, ``stretching_dvv``, or
``wavelet_dtt_dvv`` is in the lineage:

.. code-block:: python

    # Via MWCS
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                               stack=1, refstack=1,
                               mwcs=1, mwcs_dtt=1, mwcs_dtt_dvv=1)
    ds = r.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D", "1D"))

    # Via stretching
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                               stack=1, refstack=1,
                               stretching=1, stretching_dvv=1)
    ds = r.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D", "1D"))

.. seealso::

    :doc:`workflow_dvv/009_compute_dvv` — MWCS-based dv/v aggregation.
    :doc:`workflow_dvv/011c_compute_stretching_dvv` — stretching-based dv/v.
    :doc:`workflow_dvv/010d_compute_wct_dtt_dvv` — wavelet-based dv/v.


Reading PSD results
====================

.. code-block:: python

    r = MSNoiseResult.from_ids(db, psd=1, psd_rms=1)

    # Daily PSD for one seed ID (xarray Dataset)
    ds = r.get_psd("BE.UCC..HHZ", day="2023-01-01")

    # PSD RMS time series
    ds_rms = r.get_psd_rms("BE.UCC..HHZ")

.. seealso::

    :doc:`workflow_psd/009_compute_psd` — PSD computation.
    :doc:`workflow_psd/010_compute_psd_rms` — PSD RMS computation.


Navigating branches
====================

When a lineage has multiple downstream steps (e.g. both ``mwcs_1`` and
``stretching_1`` follow ``refstack_1``), use :meth:`branches` to enumerate
them:

.. code-block:: python

    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                               stack=1, refstack=1)
    for branch in r.branches():
        print(branch)
        # MSNoiseResult(category='mwcs', lineage='preprocess_1/.../mwcs_1', ...)


Exporting dv/v with full provenance
=====================================

:meth:`export_dvv` writes NetCDF files with embedded lineage, config
parameters, and station metadata — useful for archiving results or sharing
with collaborators who do not have access to the MSNoise project database:

.. code-block:: python

    r = MSNoiseResult.list(db, "mwcs_dtt_dvv")[0]
    written = r.export_dvv("exports/")
    for f in written:
        print(f)
    # dvv_CC_ZZ__pre1-cc1-f1-stk1-ref1-mwcs1-dtt1-dvv1__m1D-1D.nc

    # Reload and inspect provenance
    import xarray as xr, yaml
    ds = xr.open_dataset(written[0])
    params = yaml.safe_load(ds.attrs["msnoise_params"])
    print(params["mwcs"]["mwcs_wlen"])


API Reference
=============

.. autoclass:: msnoise.results.MSNoiseResult
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: msnoise.results
   :no-index: