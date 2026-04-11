"""
Wavelet Coherence Transform (WCT) Computation
This script performs the computation of the Wavelet Coherence Transform (WCT), a tool used to analyze the correlation between two time series in the time-frequency domain. The script supports parallel processing and interacts with a database to manage job statuses.

Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* |dtt_minlag|
* |dtt_lag|
* |dtt_v|
* |dtt_width|
* |dtt_maxdt|
* |wct_mincoh|
* |wct_codacycles|
* |wct_ns|
* |wct_nt|
* |wct_vpo|
* |wct_nptsfreq|
* |wct_min_nonzero|
* |wct_norm|
* |hpc|

This process is job-based, so it is possible to run several instances in
parallel.

To run this step:

.. code-block:: sh

    $ msnoise cc dtt compute_wct

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dtt compute_wct

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

"""

import time
import numpy as np
import xarray as xr
from .core.db import connect, get_logger
from .core.workflow import (compute_rolling_ref, extend_days, get_next_lineage_batch, get_t_axis, is_next_job_for_step, massive_update_job, propagate_downstream, refstack_is_rolling)
from .core.signal import xwt, prepare_ref_wct, apply_wct, get_wavelet_type
from .core.io import xr_get_ccf, xr_get_ref, xr_save_wct

def main(loglevel="INFO"):
    """
    Main function to process WCT jobs using lineage-based approach
    """
    global logger
    logger = get_logger('msnoise.wavelet', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT ***')

    db = connect()

    while is_next_job_for_step(db, step_category="wavelet"):
        batch = get_next_lineage_batch(db, step_category="wavelet", group_by="pair_lineage",
                                       loglevel=loglevel)

        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        lineage_names = batch["lineage_names_upstream"]
        lineage_names_mov = batch["lineage_names_mov"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]
        root = params.global_.output_folder

        logger.info(f"New WCT Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        taxis = get_t_axis(params)

        station1, station2 = pair.split(":")

        if station1 == station2:
            components_to_compute = params.cc.components_to_compute_single_station
        else:
            components_to_compute = params.cc.components_to_compute

        mov_stacks = params.stack.mov_stack
        goal_sampling_rate = params.cc.cc_sampling_rate

        # Filter frequency range from the wavelet step config
        freqmin = params.wavelet.wct_freqmin
        freqmax = params.wavelet.wct_freqmax

        # Wavelet computation parameters from the wavelet step config
        ns = params.wavelet.wct_ns
        nt = params.wavelet.wct_nt
        vpo = params.wavelet.wct_vpo
        nptsfreq = params.wavelet.wct_nptsfreq
        wct_norm = params.wavelet.wct_norm
        _wt_raw = (params.wavelet.wavelet_type or "").strip()
        try:
            wavelet_type = tuple(eval(_wt_raw)) if _wt_raw else ("Morlet", 6.)
        except Exception:
            logger.warning(f"Could not parse wavelet_type {_wt_raw!r}, using Morlet(6)")
            wavelet_type = ("Morlet", 6.)

        logger.info(f"WCT params: freqmin={freqmin} freqmax={freqmax} ns={ns} nt={nt} vpo={vpo}")

        # Build the mother wavelet object once — same for all components/stacks/days
        mother = get_wavelet_type(wavelet_type)

        for component in components_to_compute:
            rolling_mode = refstack_is_rolling(params)
            ref_wct_data = None  # pre-computed ref CWT for Mode A (set below)

            if not rolling_mode:
                # Mode A: load fixed REF from disk
                try:
                    ref_da = xr_get_ref(root, lineage_names, station1, station2,
                                        component, taxis, ignore_network=True)
                    ref = ref_da.values
                    if wct_norm:
                        ori_waveform = ref / ref.max()
                    else:
                        ori_waveform = ref
                except FileNotFoundError as fp:
                    logger.error(f"FILE DOES NOT EXIST: {fp}, skipping")
                    continue
                except Exception as e:
                    logger.error(f"Error getting reference waveform: {str(e)}")
                    continue

                if not len(ref):
                    continue

                # Pre-compute CWT + smoothed power of the reference ONCE.
                # In Mode A the reference never changes, so we avoid repeating
                # this O(N·S) computation for every day in the time loop.
                try:
                    ref_wct_data = prepare_ref_wct(
                        ori_waveform, goal_sampling_rate,
                        int(ns), float(nt), int(vpo),
                        freqmin, freqmax, int(nptsfreq), mother
                    )
                    # Extract freqs/coi for Dataset construction later
                    _, _, _, freqs, coi, _, _, _ = ref_wct_data
                except Exception as e:
                    logger.error(f"Error preparing ref WCT: {str(e)}, falling back to xwt()")
                    ref_wct_data = None
            # Mode B: ori_waveform and ref_wct_data set per time-step below

            for mov_stack in mov_stacks:
                WXamp_list = []
                WXcoh_list = []
                WXdt_list = []
                dates_list = []

                try:
                    data = xr_get_ccf(root, lineage_names_mov, station1, station2,
                                      component, mov_stack, taxis)
                except FileNotFoundError as fp:
                    logger.error(f"FILE DOES NOT EXIST: {fp}, skipping")
                    continue

                # ── Time filtering ───────────────────────────────────────
                to_search = extend_days(days)
                times_floor = data.coords["times"].values.astype("datetime64[D]")
                mask = np.isin(times_floor, np.array(list(to_search), dtype="datetime64[D]"))
                data = data.isel(times=mask)
                data = data.dropna("times", how="all")

                if rolling_mode:
                    ref_rolling = compute_rolling_ref(
                        data, int(params.refstack.ref_begin), int(params.refstack.ref_end)
                    )

                # Materialise the full CCF array once before the time loop —
                # avoids a disk read per iteration now that data is lazy.
                data_np = data.values  # shape: (times, taxis)

                for _i_row, date in enumerate(data.coords["times"].values):
                    waveform = data_np[_i_row]
                    if wct_norm:
                        new_waveform = waveform / waveform.max()
                    else:
                        new_waveform = waveform

                    if rolling_mode:
                        _rr = ref_rolling[_i_row]
                        ori_waveform = _rr / _rr.max() if wct_norm else _rr
                        # Mode B: no pre-computed ref (changes every step)
                        cur_ref_data = None
                    else:
                        cur_ref_data = ref_wct_data

                    try:
                        if cur_ref_data is not None:
                            # Mode A fast path: ref CWT already computed
                            WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = apply_wct(
                                cur_ref_data, new_waveform, int(ns), float(nt)
                            )
                        else:
                            # Mode B or fallback: full xwt each call
                            WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                                ori_waveform, new_waveform, goal_sampling_rate,
                                int(ns), float(nt), int(vpo),
                                freqmin, freqmax, int(nptsfreq), wavelet_type
                            )
                        WXamp_list.append(WXamp)
                        WXcoh_list.append(Wcoh)
                        WXdt_list.append(WXdt)
                        dates_list.append(date)
                    except Exception as e:
                        logger.error(f"Error in WCT for {date}: {str(e)}")
                        continue

                if dates_list:
                    try:
                        wct_ds = xr.Dataset({
                            "WXamp": xr.DataArray(
                                np.array(WXamp_list).real,
                                dims=["times", "freqs", "taxis"],
                                coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
                            ),
                            "Wcoh": xr.DataArray(
                                np.array(WXcoh_list).real,
                                dims=["times", "freqs", "taxis"],
                                coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
                            ),
                            "WXdt": xr.DataArray(
                                np.array(WXdt_list).real,
                                dims=["times", "freqs", "taxis"],
                                coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
                            ),
                        })
                        xr_save_wct(root, lineage_names, step.step_name,
                                    station1, station2, component, mov_stack,
                                    wct_ds)
                    except Exception as e:
                        logger.error(f"Error saving WCT: {str(e)}")

        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)

    db.close()
    logger.info('*** Finished: Compute WCT ***')
