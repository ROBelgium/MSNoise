"""
Reference Stack computation for MSNoise 2.x configsets workflow.

This step computes (or validates) the REF stack for a ``refstack`` configset.
It is positioned downstream of the MOV ``stack`` step and upstream of
``mwcs`` / ``stretching`` / ``wavelet``.

Two modes are supported, determined by ``ref_begin`` in the refstack configset:

**Mode A — Fixed REF** (``ref_begin`` is a date or ``"1970-01-01"``)
    Reads daily CCFs from the upstream ``stack_N`` output folder, stacks the
    windows falling within ``[ref_begin, ref_end]``, and writes a single REF
    NetCDF file under the ``refstack_M`` step folder::

        OUTPUT/.../stack_N/refstack_M/_output/REF/<components>/sta1_sta2.nc

**Mode B — Rolling REF** (``ref_begin`` is a negative integer string, e.g. ``"-5"``)
    No file is written.  The job validates that MOV data exists for the pair
    and marks itself Done immediately.  The actual rolling reference is
    computed on-the-fly inside the mwcs / stretching / wavelet workers via
    :func:`~msnoise.api.compute_rolling_ref`.

Configuration parameters (from the ``refstack`` configset):

* ``ref_begin``
* ``ref_end``
* ``stack_method``
* ``pws_timegate`` / ``pws_power``  (Mode A + pws only)
* ``export_format`` / ``sac_format``  (Mode A only)

To run this step:

.. code-block:: sh

    $ msnoise cc stack_refstack

This step also supports parallel processing:

.. code-block:: sh

    $ msnoise -t 4 cc stack_refstack
"""

import time

import numpy as np
import pandas as pd

from .api import (
    connect,
    get_logger,
    get_params,
    get_next_lineage_batch,
    is_next_job_for_step,
    get_results_all,
    get_results,
    massive_update_job,
    build_ref_datelist,
    refstack_is_rolling,
    xr_save_ref,
    get_t_axis,
    validate_stack_data,
)
from .wiener import wiener_filt


def main(loglevel="INFO"):
    """Compute REF stacks for all pending ``refstack`` jobs."""

    logger = get_logger("msnoise.refstack_child", loglevel, with_pid=True)
    logger.info("*** Starting: Compute REF Stack (refstack) ***")

    db = connect()
    orig_params = get_params(db)

    while is_next_job_for_step(db, step_category="refstack"):
        logger.info("Getting the next refstack job")

        batch = get_next_lineage_batch(
            db,
            step_category="refstack",
            group_by="pair_lineage",
            loglevel=loglevel,
            drop_current_step_name=True,
        )
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs        = batch["jobs"]
        pair        = batch["pair"]
        params      = batch["params"]
        # lineage_names has refstack_M dropped (drop_current_step_name=True),
        # so it ends with stack_N — exactly where the CCF data lives.
        lineage_names = batch["lineage_names"]
        step          = batch["step"]
        lineage_str   = batch["lineage_str"]

        logger.info(f"New REFSTACK Job: pair={pair} lineage={lineage_str}")

        # ── Mode B: rolling — no file to write, just validate data exists ──
        if refstack_is_rolling(params):
            logger.info(
                f"Mode B (rolling index): ref_begin={params.ref_begin}, "
                f"ref_end={params.ref_end}. "
                "No REF file will be written; reference computed on-the-fly "
                "at MWCS/stretching/WCT time."
            )
            massive_update_job(db, jobs, "D")
            # Downstream mwcs/stretching/wavelet jobs are triggered by
            # `msnoise new_jobs --after refstack` — see s02new_jobs.py
            continue

        # ── Mode A: fixed-date REF ──
        taxis = get_t_axis(params)

        sta1, sta2 = pair.split(":")
        filterid = 1  # kept for path compatibility; not used functionally

        # Wiener filter parameters (reuse stack logic if configured)
        wienerfilt  = getattr(params, "wienerfilt", False)
        wiener_M    = None
        wiener_N    = None
        gap_threshold = None
        if wienerfilt:
            wiener_M = int(
                pd.to_timedelta(params.wiener_mlen).total_seconds()
                / params.corr_duration
            )
            wiener_N = int(
                pd.to_timedelta(params.wiener_nlen).total_seconds()
                * params.cc_sampling_rate
            )
            gap_threshold = (
                wiener_M if params.keep_all
                else pd.to_timedelta(params.wiener_mlen).days
            )
            logger.info("Wiener filter enabled for REF stack")

        for components in params.components_to_compute:
            logger.info(
                f"Processing {pair}-{components} REF stack (Mode A, "
                f"ref_begin={params.ref_begin}, ref_end={params.ref_end})"
            )

            start, end, datelist = build_ref_datelist(params)

            if params.keep_all:
                c = get_results_all(
                    db, params.output_folder, lineage_names,
                    sta1, sta2, components, datelist,
                    format="xarray", params=params,
                )
            else:
                c = get_results(
                    db, sta1, sta2, filterid, components, datelist,
                    mov_stack=1, format="xarray", params=params,
                )

            is_valid, message = validate_stack_data(c, "reference")
            if not is_valid:
                logger.error(
                    f"Invalid reference data for {sta1}:{sta2}-{components}: "
                    f"{message}"
                )
                continue
            if "Warning" in message:
                logger.warning(f"{sta1}:{sta2}-{components}: {message}")

            if not c.data_vars:
                logger.debug(f"No data found for {sta1}:{sta2}-{components}")
                continue

            if wienerfilt:
                c = wiener_filt(c, wiener_M, wiener_N, gap_threshold)

            ref_stack = c.mean(dim="times")

            xr_save_ref(
                params.output_folder,
                lineage_names,    # ends with stack_N
                step.step_name,   # refstack_M  — becomes the step sub-folder
                sta1, sta2, components, filterid, taxis, ref_stack,
            )
            logger.info(
                f"REF stack written: "
                f"{params.output_folder}/"
                f"{'/'.join(lineage_names)}/{step.step_name}/_output/REF/"
                f"{components}/{sta1}_{sta2}.nc"
            )

        massive_update_job(db, jobs, "D")
        logger.info(
            f"Marked {len(jobs)} refstack job(s) as Done for pair={pair}"
        )
