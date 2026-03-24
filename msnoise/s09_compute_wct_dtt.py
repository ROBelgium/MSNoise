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
import pandas as pd

from .api import (
    connect, get_logger, is_next_job_for_step, get_next_lineage_batch,
    get_station_pairs, get_interstation_distance,
    xr_load_wct, xr_save_wct_dtt2, massive_update_job,
    compute_wct_dvv, get_wct_avgcoh,
)


def main(loglevel="INFO", batch_size=None):
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
    db.close()

    while True:
        db = connect()
        if not is_next_job_for_step(db, step_category="wavelet_dtt"):
            db.close()
            break

        batch = get_next_lineage_batch(
            db, step_category="wavelet_dtt", group_by="pair_lineage",
            loglevel=loglevel, drop_current_step_name=False,
        )
        db.close()

        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        # lineage_names[:-1] ends with the upstream "wavelet" step name, used
        # for xr_load_wct (read).  The full lineage (including wavelet_dtt step)
        # is reconstructed via  lineage_names[:-1] + step.step_name  inside
        # xr_save_wct_dtt2, yielding:  <lineage>/wavelet_dtt_N/_output/...
        lineage_names = batch["lineage_names"][:-1]
        lineage_str = batch["lineage_str"]
        step = batch["step"]
        root = params.output_folder

        logger.info(
            f"New WCT DTT Job: pair={pair} n_days={len(days)} lineage={lineage_str}"
        )

        station1, station2 = pair.split(":")
        components_to_compute = (
            params.components_to_compute_single_station
            if station1 == station2
            else params.components_to_compute
        )
        mov_stacks = params.mov_stack

        # DVV computation parameters from the wavelet_dtt step config.
        # Fall back to wavelet-step freq bounds when dtt-specific bounds are unset.
        dtt_freqmin = float(getattr(params, 'wct_dtt_freqmin', None) or
                            getattr(params, 'wct_freqmin', 0.1))
        dtt_freqmax = float(getattr(params, 'wct_dtt_freqmax', None) or
                            getattr(params, 'wct_freqmax', 2.0))
        minlag = float(getattr(params, 'wct_minlag', 5.0))
        lag_type = str(getattr(params, 'wct_lag', None) or 'static')
        width = float(getattr(params, 'wct_width', 30.0))
        v = float(getattr(params, 'wct_v', 1.0))
        sides = str(getattr(params, 'wct_sides', None) or 'both')
        mincoh = float(getattr(params, 'wct_mincoh', 0.0))
        maxdt = float(getattr(params, 'wct_maxdt', 1.0))
        coda_cycles = int(getattr(params, 'wct_codacycles', 5))
        min_nonzero = float(getattr(params, 'wct_min_nonzero', 0.1))

        # Resolve interstation distance for dynamic lag
        n1, s1, _l1 = station1.split(".")
        n2, s2, _l2 = station2.split(".")
        dpair = "%s_%s_%s_%s" % (n1, s1, n2, s2)
        dist = interstations.get(dpair, 0.0)

        if lag_type == "dynamic" and v > 0:
            lag_min = dist / v
        else:
            lag_min = minlag

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

                times_dt = pd.DatetimeIndex(times)
                to_search = pd.to_datetime(days)

                start = to_search.min().normalize()
                end = (to_search.max() + pd.Timedelta(days=1)).normalize()

                valid_mask = (times_dt >= start) & (times_dt <= end)

                dates_out = []
                dvv_rows = []
                err_rows = []
                coh_rows = []

                for i, t in enumerate(times):
                    if not valid_mask[i]:
                        continue

                    WXamp_i = ds["WXamp"].values[i]
                    Wcoh_i = ds["Wcoh"].values[i]
                    WXdt_i = ds["WXdt"].values[i]

                    try:
                        dvv, err, _ = compute_wct_dvv(
                            freqs, taxis, WXamp_i, Wcoh_i, WXdt_i,
                            lag_min=lag_min, coda_cycles=coda_cycles,
                            mincoh=mincoh, maxdt=maxdt,
                            min_nonzero=min_nonzero,
                            freqmin=dtt_freqmin, freqmax=dtt_freqmax,
                        )
                        coh = get_wct_avgcoh(
                            freqs, taxis, Wcoh_i,
                            freqmin=dtt_freqmin, freqmax=dtt_freqmax,
                            lag_min=lag_min, coda_cycles=coda_cycles,
                        )
                    except Exception as e:
                        logger.error(
                            f"DVV computation error for {pair}/{component} at {t}: {e}"
                        )
                        continue

                    dates_out.append(pd.Timestamp(t))
                    dvv_rows.append(dvv)
                    err_rows.append(err)
                    coh_rows.append(coh)

                del ds

                if not dates_out:
                    logger.warning(
                        f"No DVV results for {pair}/{component}/{mov_stack}"
                    )
                    continue

                dvv_df = pd.DataFrame(dvv_rows, index=dates_out, columns=freqs_subset)
                err_df = pd.DataFrame(err_rows, index=dates_out, columns=freqs_subset)
                coh_df = pd.DataFrame(coh_rows, index=dates_out, columns=freqs_subset)
                dvv_df.sort_index(inplace=True)
                err_df.sort_index(inplace=True)
                coh_df.sort_index(inplace=True)

                try:
                    # Output path: root/<lineage_names>/<step.step_name>/_output/...
                    # e.g. root/preprocess_1/cc_1/filter_1/stack_1/wavelet_1/wavelet_dtt_1/_output/...
                    xr_save_wct_dtt2(
                        root, lineage_names, step.step_name,
                        station1, station2, component, mov_stack,
                        taxis, dvv_df, err_df, coh_df,
                    )
                    logger.info(
                        f"Saved WCT DTT for {pair}/{component}/{mov_stack} "
                        f"({len(dates_out)} time steps)"
                    )
                except Exception as e:
                    logger.error(
                        f"Error saving WCT DTT for {pair}/{component}/{mov_stack}: {e}"
                    )

        db = connect()
        massive_update_job(db, jobs, "D")
        db.close()

    logger.info('*** Finished: Compute WCT DTT ***')
