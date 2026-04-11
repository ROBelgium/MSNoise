#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Wavelet DTT Computation

This module computes dv/v from the saved WCT (Wavelet Coherence Transform)
results using a lineage-based job approach. Each job processes one station
pair across all components and moving stacks.

The WCT data is loaded from the upstream ``wavelet`` step output and the
resulting dv/v, error, and coherence DataFrames are saved under the
``wavelet_dtt`` step path:

.. code-block:: text

    <output_folder>/<lineage>/<wavelet_dtt_step>/_output/<mov_stack>/<comp>/<sta1>_<sta2>.nc

To run this step:

.. code-block:: sh

    msnoise cc dtt compute_wct_dtt
"""

import time

import numpy as np

from .core.db import connect, get_logger
from .core.stations import get_interstation_distance, get_station_pairs
from .core.workflow import (extend_days, get_next_lineage_batch, is_next_job_for_step, massive_update_job, propagate_downstream)
from .core.compute import compute_wct_dtt_batch, build_wct_dtt_dataset
from .core.io import xr_load_wct, xr_save_wct_dtt


def main(loglevel="INFO"):
    """
    Main function to compute dv/v from WCT results using a lineage-based approach.

    Reads accumulated per-pair WCT data written by the ``wavelet`` step and
    computes dv/v for each (component, mov_stack) combination using the
    ``wavelet_dtt`` configuration parameters merged into the lineage.

    Output is stored at::

        <output_folder>/<full_lineage>/_output/<mov_stack>/<comp>/<sta1>_<sta2>.nc
    """
    logger = get_logger('msnoise.wavelet_dtt', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT DTT ***')

    # Pre-compute interstation distances once (needed for dynamic lag)
    db = connect()
    interstations = {}
    for sta1, sta2 in get_station_pairs(db):
        s1 = "%s_%s" % (sta1.net, sta1.sta)
        s2 = "%s_%s" % (sta2.net, sta2.sta)
        if s1 == s2:
            interstations["%s_%s" % (s1, s2)] = 0.0
        else:
            interstations["%s_%s" % (s1, s2)] = get_interstation_distance(
                sta1, sta2, sta1.coordinates
            )

    while is_next_job_for_step(db, step_category="wavelet_dtt"):
        batch = get_next_lineage_batch(
            db, step_category="wavelet_dtt", group_by="pair_lineage",
            loglevel=loglevel,
        )

        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        # lineage_names_upstream ends with the upstream "wavelet" step name,
        # used for xr_load_wct (read).  The full lineage (including wavelet_dtt)
        # is reconstructed via lineage_names_upstream + step.step_name inside
        # xr_save_wct_dtt, yielding: <lineage>/wavelet_dtt_N/_output/...
        lineage_names = batch["lineage_names_upstream"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]
        root = params.global_.output_folder

        logger.info(
            f"New WCT DTT Job: pair={pair} n_days={len(days)} lineage={lineage_str}"
        )

        station1, station2 = pair.split(":")
        components_to_compute = (
            params.cc.components_to_compute_single_station
            if station1 == station2
            else params.cc.components_to_compute
        )
        mov_stacks = params.stack.mov_stack

        # DVV computation parameters from the wavelet_dtt step config.
        # Fall back to wavelet-step freq bounds when dtt-specific bounds are unset.
        dtt_freqmin = params.wavelet_dtt.wct_dtt_freqmin
        dtt_freqmax = params.wavelet_dtt.wct_dtt_freqmax

        # Resolve interstation distance for dynamic lag
        n1, s1, _l1 = station1.split(".")
        n2, s2, _l2 = station2.split(".")
        dpair = "%s_%s_%s_%s" % (n1, s1, n2, s2)
        dist = interstations.get(dpair, 0.0)

        for component in components_to_compute:
            for mov_stack in mov_stacks:
                try:
                    ds = xr_load_wct(
                        root, lineage_names, station1, station2, component, mov_stack
                    )
                except FileNotFoundError:
                    logger.warning(
                        f"No WCT data found for {pair}/{component}/{mov_stack}, skipping"
                    )
                    continue
                except Exception as e:
                    logger.error(
                        f"Error loading WCT data for {pair}/{component}/{mov_stack}: {e}"
                    )
                    continue

                freqs = ds.coords["freqs"].values
                taxis = ds.coords["taxis"].values
                times = ds.coords["times"].values

                freq_inx = np.where((freqs >= dtt_freqmin) & (freqs <= dtt_freqmax))[0]
                if freq_inx.size == 0:
                    logger.warning(
                        f"No frequencies in [{dtt_freqmin}, {dtt_freqmax}] Hz "
                        f"for {pair}, skipping"
                    )
                    continue
                freqs_subset = freqs[freq_inx]

                # Use extend_days to handle label="right" resampling:
                # daily mov_stack stamps are at midnight of the *next* day,
                # so we extend the search window by one extra day.
                to_search = extend_days(days)  # pandas DatetimeIndex + 1 extra day
                times_np  = times.astype("datetime64[D]")
                to_search_d = np.array(to_search, dtype="datetime64[D]")
                valid_mask  = np.isin(times_np, to_search_d)

                # Materialise the full 3D WCT arrays once before the time loop
                # so that each ds[var].values[i] is a cheap numpy index, not a
                # repeated lazy read that caches the full array on first access.
                # Consistent with s05's data_arr = data.values pattern.
                wxamp_np = ds["WXamp"].values   # (times, freqs, taxis)
                wcoh_np  = ds["Wcoh"].values
                wxdt_np  = ds["WXdt"].values
                ds.close()
                del ds

                dates_out = []
                dtt_rows = []
                err_rows = []
                coh_rows = []

                for i, t in enumerate(times):
                    if not valid_mask[i]:
                        continue

                    WXamp_i = wxamp_np[i]
                    Wcoh_i  = wcoh_np[i]
                    WXdt_i  = wxdt_np[i]

                    try:
                        dtt, err, coh, freqs_subset = compute_wct_dtt_batch(
                            freqs, taxis, WXamp_i, Wcoh_i, WXdt_i,
                            params, dist=dist,
                        )
                    except Exception as e:
                        logger.error(
                            f"DVV computation error for {pair}/{component} at {t}: {e}"
                        )
                        continue

                    dates_out.append(t)
                    dtt_rows.append(dtt)
                    err_rows.append(err)
                    coh_rows.append(coh)

                del wxamp_np, wcoh_np, wxdt_np  # free 3D WCT arrays

                if not dates_out:
                    logger.warning(
                        f"No DTT results for {pair}/{component}/{mov_stack}"
                    )
                    continue

                ds_out = build_wct_dtt_dataset(
                    dates_out, dtt_rows, err_rows, coh_rows, freqs_subset,
                )

                try:
                    # Output path: root/<lineage_names>/<step.step_name>/_output/...
                    # e.g. root/preprocess_1/cc_1/filter_1/stack_1/wavelet_1/wavelet_dtt_1/_output/...
                    xr_save_wct_dtt(
                        root, lineage_names, step.step_name,
                        station1, station2, component, mov_stack,
                        taxis, ds_out,
                    )
                    logger.info(
                        f"Saved WCT DTT for {pair}/{component}/{mov_stack} "
                        f"({len(dates_out)} time steps)"
                    )
                except Exception as e:
                    logger.error(
                        f"Error saving WCT DTT for {pair}/{component}/{mov_stack}: {e}"
                    )

        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)

    db.close()
    logger.info('*** Finished: Compute WCT DTT ***')
