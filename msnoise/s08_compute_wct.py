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

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
from .api import *

def main(loglevel="INFO", batch_size=None):
    """
    Main function to process WCT jobs using lineage-based approach
    """
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel, with_pid=True)
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
        root = params.output_folder

        logger.info(f"New WCT Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        taxis = get_t_axis(params)

        station1, station2 = pair.split(":")

        if station1 == station2:
            components_to_compute = params.components_to_compute_single_station
        else:
            components_to_compute = params.components_to_compute

        mov_stacks = params.mov_stack
        goal_sampling_rate = params.cc_sampling_rate

        # Filter frequency range from the wavelet step config
        freqmin = params.wct_freqmin
        freqmax = params.wct_freqmax

        # Wavelet computation parameters from the wavelet step config
        ns = getattr(params, 'wct_ns', 5)
        nt = getattr(params, 'wct_nt', 5)
        vpo = getattr(params, 'wct_vpo', 20)
        nptsfreq = getattr(params, 'wct_nptsfreq', 300)
        wct_norm = getattr(params, 'wct_norm', "Y") == "Y"
        wavelet_type = eval(getattr(params, 'wavelet_type', '') or "('Morlet',6.)")

        logger.info(f"WCT params: freqmin={freqmin} freqmax={freqmax} ns={ns} nt={nt} vpo={vpo}")

        for component in components_to_compute:
            rolling_mode = refstack_is_rolling(params)
            if not rolling_mode:
                # Mode A: load fixed REF from disk
                try:
                    ref_data = xr_get_ref(root, lineage_names, station1, station2,
                                          component, taxis, ignore_network=True)
                    ref = ref_data.CCF.values
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
            # Mode B: ori_waveform set per time-step below after data is loaded

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

                to_search = _extend_days(days)
                data = data[data.index.floor('d').isin(to_search)]
                data = data.dropna()

                if rolling_mode:
                    ref_rolling = compute_rolling_ref(
                        data, int(params.ref_begin), int(params.ref_end)
                    )

                for _i_row, (date, row) in enumerate(data.iterrows()):
                    waveform = row.values
                    if wct_norm:
                        new_waveform = waveform / waveform.max()
                    else:
                        new_waveform = waveform

                    if rolling_mode:
                        _rr = ref_rolling[_i_row]
                        ori_waveform = _rr / _rr.max() if wct_norm else _rr

                    try:
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
                        xr_save_wct(root, lineage_names, step.step_name,
                                     station1, station2, component, mov_stack,
                                     taxis, freqs, WXamp_list, WXcoh_list,
                                     WXdt_list, dates_list)
                    except Exception as e:
                        logger.error(f"Error saving WCT: {str(e)}")

        massive_update_job(db, jobs, "D")

    db.close()
    logger.info('*** Finished: Compute WCT ***')
