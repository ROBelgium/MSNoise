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

from .core.db import connect, get_logger
from .core.workflow import (build_ref_datelist, get_next_lineage_batch, get_t_axis, is_next_job_for_step, massive_update_job, propagate_downstream, refstack_is_rolling, refstack_needs_recompute)
from .core.signal import validate_stack_data
from .core.io import xr_load_ccf_for_stack, xr_save_ref
from .wiener import wiener_filt


def main(loglevel="INFO"):
    """Compute REF stacks for all pending ``refstack`` jobs."""

    logger = get_logger("msnoise.refstack", loglevel, with_pid=True)
    logger.info("*** Starting: Compute REF Stack (refstack) ***")

    db = connect()
    while is_next_job_for_step(db, step_category="refstack"):
        logger.info("Getting the next refstack job")

        batch = get_next_lineage_batch(
            db,
            step_category="refstack",
            group_by="pair_lineage",
            loglevel=loglevel,
        )
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs        = batch["jobs"]
        pair        = batch["pair"]
        params      = batch["params"]
        # lineage_names_upstream ends with stack_N (current step refstack_M excluded)
        lineage_names = batch["lineage_names_upstream"]
        # lineage_names_cc strips stack_N too, ending at filter_N
        # where raw daily CC h5 files live under filter_N/_output/all/...
        # Note: cannot use lineage_names_mov here — that strips refstack_* entries
        # which don't exist at this level; we need to strip stack_N explicitly.
        lineage_names_cc = lineage_names[:-1]
        step          = batch["step"]
        lineage_str   = batch["lineage_str"]

        logger.info(f"New REFSTACK Job: pair={pair} lineage={lineage_str}")

        # ── Mode B: rolling — no file to write, just validate data exists ──
        if refstack_is_rolling(params):
            logger.info(
                f"Mode B (rolling index): ref_begin={params.refstack.ref_begin}, "
                f"ref_end={params.refstack.ref_end}. "
                "No REF file will be written; reference computed on-the-fly "
                "at MWCS/stretching/WCT time."
            )
            massive_update_job(db, jobs, "D")
            if not batch["params"].global_.hpc:
                propagate_downstream(db, batch)
            continue

        # ── Mode A: fixed-date REF ──

        # Check whether any Done stack days for this pair fall inside
        # the configured ref_begin..ref_end window.  If none do, the
        # reference waveform is unaffected by the new data — skip the
        # heavy computation and mark Done immediately.  Downstream
        # mwcs/stretching/wavelet T jobs were already created by
        # propagate_downstream (stack category) so they can start right away.
        if not refstack_is_rolling(params):
            if not refstack_needs_recompute(db, pair, batch["lineage_names_upstream"], params):
                logger.info(
                    f"REF skip for {pair}: no Done stack days fall inside "
                    f"ref_begin={params.refstack.ref_begin} .. "
                    f"ref_end={params.refstack.ref_end}. "
                    "Reference unchanged — marking Done without recomputation."
                )
                massive_update_job(db, jobs, "D")
                if not batch["params"].global_.hpc:
                    propagate_downstream(db, batch)
                continue

        taxis = get_t_axis(params)

        sta1, sta2 = pair.split(":")

        # Wiener filter parameters (reuse stack logic if configured)
        wienerfilt  = params.stack.wienerfilt
        wiener_M    = None
        wiener_N    = None
        gap_threshold = None
        if wienerfilt:
            wiener_M = int(
                pd.to_timedelta(params.stack.wiener_mlen).total_seconds()
                / params.cc.corr_duration
            )
            wiener_N = int(
                pd.to_timedelta(params.stack.wiener_nlen).total_seconds()
                * params.cc.cc_sampling_rate
            )
            gap_threshold = (
                wiener_M if params.cc.keep_all
                else pd.to_timedelta(params.stack.wiener_mlen).days
            )
            logger.info("Wiener filter enabled for REF stack")

        if sta1 == sta2:
            components_to_compute = params.cc.components_to_compute_single_station
        else:
            components_to_compute = params.cc.components_to_compute

        for components in components_to_compute:
            logger.info(
                f"Processing {pair}-{components} REF stack (Mode A, "
                f"ref_begin={params.refstack.ref_begin}, ref_end={params.refstack.ref_end})"
            )

            start, end, datelist = build_ref_datelist(params)

            if params.cc.keep_all:
                c = xr_load_ccf_for_stack(
                    params.global_.output_folder, lineage_names_cc,
                    sta1, sta2, components, datelist,
                )
            else:
                logger.warning(
                    "keep_all=N is unsupported in lineage workflow; "
                    "falling back to keep_days daily stacks"
                )
                c = xr_load_ccf_for_stack(
                    params.global_.output_folder, lineage_names_cc,
                    sta1, sta2, components, datelist,
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

            # c is a Dataset with variable "CCF" (dims: times, taxis)
            # mean over times → Dataset with variable "CCF" (dim: taxis)
            # xr_save_ref expects a "REF"-named DataArray/Dataset,
            # so extract and rename before saving.
            ref_stack = c["CCF"].mean(dim="times").rename("REF")

            xr_save_ref(
                params.global_.output_folder,
                lineage_names,    # ends with stack_N
                step.step_name,   # refstack_M  — becomes the step sub-folder
                sta1, sta2, components, taxis, ref_stack.to_dataset(),
            )
            logger.info(
                f"REF stack written: "
                f"{params.global_.output_folder}/"
                f"{'/'.join(lineage_names)}/{step.step_name}/_output/REF/"
                f"{components}/{sta1}_{sta2}.nc"
            )

        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)
        logger.info(
            f"Marked {len(jobs)} refstack job(s) as Done for pair={pair}"
        )
