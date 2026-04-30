"""Moving stack computation.

Reads per-day or per-window CCF NetCDF files written by the CC step and
produces **moving-window stacked CCFs** for every ``(pair, component,
mov_stack)`` combination.

The stacking window is defined by ``|stack.mov_stack|`` as a list of
``(window, step)`` tuples.  For example ``('7D', '1D')`` produces a
rolling 7-day mean, sampled daily.  Each output is written as a NetCDF file
under the lineage output path::

    OUTPUT/.../stack_N/_output/<mov_stack>/<component>/<sta1>_<sta2>.nc

When ``|stack.wienerfilt|`` is ``Y``, a Wiener filter is applied to the
CCF time series before the rolling mean, attenuating isolated spikes
and filling short gaps.

This step is **upstream of refstack**: once Done, ``propagate_downstream``
creates both a REF sentinel job (triggering ``s04_stack_refstack``)
and direct MWCS/stretching/wavelet T jobs for the new days.

To run this step:

.. code-block:: sh

    $ msnoise cc stack

Parallel processing:

.. code-block:: sh

    $ msnoise -t 4 cc stack

Configuration Parameters
------------------------

* |stack.mov_stack|
* |stack.wienerfilt|
* |stack.wiener_mlen|
* |stack.wiener_nlen|
* |cc.keep_all|
* |cc.corr_duration|
* |cc.cc_sampling_rate|
* |cc.components_to_compute|
* |cc.components_to_compute_single_station|
* |global.hpc|
"""
import datetime
import math
import time

import numpy as np
import pandas as pd
import xarray as xr
from .core.db import connect, get_logger
from .core.workflow import (get_next_lineage_batch, get_t_axis, is_next_job_for_step, massive_update_job, propagate_downstream)
from .core.signal import validate_stack_data, rolling_mean_nan, rolling_stack
from .core.io import xr_load_ccf_for_stack, xr_save_ccf
from .core.signal import wiener_filt



def main(stype, loglevel="INFO"):
    """Computes the REF/MOV stacks.

    Parameters
    ----------
    stype : {'mov', 'ref'}
        Defines which of the REF or Moving-window stacks must be exported
    interval : float, optional
        Number of days before now to search for modified CC jobs

    """
    logger = get_logger('msnoise.stack', loglevel, with_pid=True)
    logger.debug('Starting the %s stack' % stype)
    db = connect()


    while is_next_job_for_step(db, step_category="stack"):
        logger.debug("Getting the next batch")
        batch = get_next_lineage_batch(db, step_category="stack", group_by="pair_lineage", loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        lineage_names = batch["lineage_names_upstream"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]

        logger.info(f"New STACK Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        taxis = get_t_axis(params)

        mov_stacks   = params.stack.mov_stack
        stack_method = params.stack.stack_method
        pws_timegate = params.stack.pws_timegate
        pws_power    = params.stack.pws_power
        if stack_method != "linear":
            logger.info(f"Moving stack: stack_method={stack_method!r}")
        wiener_mlen = params.stack.wiener_mlen
        wiener_nlen = params.stack.wiener_nlen
        wienerfilt = params.stack.wienerfilt
        wiener_M = int(pd.to_timedelta(wiener_mlen).total_seconds() / params.cc.corr_duration)
        wiener_N = int(pd.to_timedelta(wiener_nlen).total_seconds() * params.cc.cc_sampling_rate)

        # is there a better alternative for threshold?
        if params.cc.keep_all:
            wiener_gap_threshold = wiener_M  # no. indices which will be considered adjacent by wiener
        else:
            wiener_gap_threshold = pd.to_timedelta(wiener_mlen).days

        if wienerfilt:
            logger.info('Wiener filter enabled, will apply to CCFs before stacking')

        logger.info(
            "There are STACKS jobs for some days to recompute for %s" % pair)

        sta1, sta2 = pair.split(':')
        if sta1 == sta2:
            components_to_compute = params.cc.components_to_compute_single_station
        else:
            components_to_compute = params.cc.components_to_compute

        for components in components_to_compute:
            logger.info('Processing %s-%s MOV stack' % (pair, components))

            # Calculate the maximum mov_rolling value (in days)
            max_mov_rolling = max(pd.to_timedelta(mov_stack[0]).total_seconds() for mov_stack in mov_stacks)

            if wienerfilt:
                wiener_mlen_days = math.ceil(pd.to_timedelta(wiener_mlen).total_seconds() / 86400)
                max_mov_rolling_days = max(1, math.ceil(max_mov_rolling / 86400), 2*wiener_mlen_days) #2*wiener to minimise edge effect
            else:
                max_mov_rolling_days = int(max(1, math.ceil(max_mov_rolling / 86400)))

            days = list(days)
            days.sort()
            days = [day if isinstance(day, datetime.datetime) else datetime.datetime.strptime(day, '%Y-%m-%d') for day in days]
            day_diffs = np.diff(days)
            gaps = [i+1 for i, diff in enumerate(day_diffs) if diff.days > 1] #get index of days with gaps
            del day_diffs  # no longer needed; free timedelta array
            gaps.insert(0,0) #zero index also 'gap' (need previous data for stacking)

            all_days = list(days)
            wiener_extra_days = []
            excess_days = [] #excess days added for padding, for later removal

            for gap_idx in gaps:
                #Add days before start beginning of new segment
                start = days[gap_idx]
                for j in range(1, max_mov_rolling_days+1):
                    preceding_day = start - datetime.timedelta(days=j)
                    if preceding_day not in days: #if not already present
                        all_days.append(preceding_day)
                        # Handle excess days
                        if not wienerfilt:  # If Wiener filter is not true, add all new days to excess
                            excess_days.append(preceding_day)
                        elif j <= wiener_mlen_days:
                            #keep first half of padding window
                            wiener_extra_days.append(preceding_day)
                            if preceding_day in excess_days:
                                #if day previously marked excess, remove
                                excess_days.remove(preceding_day)
                        elif (j > wiener_mlen_days) and (preceding_day not in wiener_extra_days):
                            #mark second half of padding window as excess data to be removed pre-save
                            excess_days.append(preceding_day)

                if wienerfilt:
                    #add days at end of previous segment (only necessary if wiener filt applied)
                    end = days[gap_idx-1]
                    for j in range(1, 2*wiener_mlen_days+1): #2*wiener to minimise edge effect
                        future_day = end + datetime.timedelta(days=j)
                        if future_day not in days: #if not already present
                            all_days.append(future_day)
                            if j <= wiener_mlen_days:
                                #keep first half of padding window
                                wiener_extra_days.append(future_day)
                                if future_day in excess_days:
                                    #if day previously marked excess, remove
                                    excess_days.remove(future_day)
                            elif (j > wiener_mlen_days) and (future_day not in wiener_extra_days):
                                #mark second half of padding window as excess data to be removed pre-save
                                excess_days.append(future_day)

            all_days = sorted(set(all_days))
            excess_days = sorted(set(excess_days))

            del gaps   # gap index only needed during build above
            del wiener_extra_days  # only needed during gap-padding build above

            if params.cc.keep_all:
                # Always use the eager loader (P14 path): opens each file once,
                # reads float32 data directly, single numpy concat.
                # The lazy/dask path (xr_open_ccf_mfdataset) opens every file
                # twice on GPFS (once for timestamp metadata, once for data) and
                # gains no I/O parallelism due to HDF5's global thread lock —
                # slower in practice despite lower peak RAM.
                _t0 = time.time()
                c = xr_load_ccf_for_stack(params.global_.output_folder, lineage_names,
                                          sta1, sta2, components, all_days)
                logger.debug(f"  [timing] load CCF ({len(all_days)} days): {time.time()-_t0:.2f}s")
                if not len(c):
                    logger.warning("No data found for %s-%s" % (sta1, sta2))
                    continue
                _t0 = time.time()
                dr = c.resample(times="%is" % params.cc.corr_duration).mean()
                c.close()
                del c  # free raw CCF data — dr is all we need from here
                dr = dr.compute()
                logger.debug(f"  [timing] resample+compute → dr {dr['CCF'].shape}: {time.time()-_t0:.2f}s")

            else:
                logger.warning("keep_all=N is unsupported in lineage workflow; "
                               "falling back to keep_days daily stacks")
                _t0 = time.time()
                c = xr_load_ccf_for_stack(params.global_.output_folder, lineage_names,
                                          sta1, sta2, components, all_days)
                if not len(c):
                    logger.warning("No data found for %s-%s" % (sta1, sta2))
                    continue
                logger.debug(f"  [timing] load CCF keep_days ({len(all_days)} days): {time.time()-_t0:.2f}s")
                _t0 = time.time()
                dr = c.resample(times="1D").mean()
                c.close()
                del c  # free raw CCF data — dr is all we need from here
                dr = dr.compute()
                logger.debug(f"  [timing] resample+compute → dr {dr['CCF'].shape}: {time.time()-_t0:.2f}s")

            if wienerfilt:
                _t0 = time.time()
                dr = wiener_filt(dr, wiener_M, wiener_N, wiener_gap_threshold)
                logger.debug(f"  [timing] wiener_filt: {time.time()-_t0:.2f}s")

            # Validate the resampled stack (dr), not the raw windowed data (c).
            # dr is what gets saved — it's smaller, already averaged, and the
            # NaN fraction here is what actually matters for downstream steps.
            is_valid, message = validate_stack_data(dr, "moving")
            if not is_valid:
                logger.error(f"Invalid moving stack data for {sta1}:{sta2}-{components}: {message}")
                del dr
                continue
            elif "Warning" in message:
                logger.warning(f"{sta1}:{sta2}-{components}: {message}")

            excess_dates = pd.to_datetime(excess_days).values
            del excess_days  # numpy excess_dates is all we need from here
            for mov_stack in mov_stacks:
                mov_rolling, mov_sample = mov_stack
                _t0 = time.time()

                if mov_rolling == mov_sample:
                    xx = dr.resample(times=mov_sample, label="right", skipna=True).mean().dropna("times", how="all")
                else:
                    mov_rolling = pd.to_timedelta(mov_rolling).total_seconds()
                    if params.cc.keep_all:
                        duration_to_windows = mov_rolling / params.cc.corr_duration
                    else:
                        duration_to_windows = mov_rolling / 86400.0
                    if not duration_to_windows.is_integer():
                        logger.print("Warning, rounding down the number of windows to roll over")
                    duration_to_windows = int(max(1, math.floor(duration_to_windows)))
                    # rolling_stack: O(N) for linear, O(N×W) for pws/tfpws.
                    # xarray rolling is O(N×W) when NaN handling disables
                    # its cumsum fast-path; rolling_mean_nan avoids that.
                    rolled = rolling_stack(
                        dr["CCF"].values, duration_to_windows,
                        stack_method=stack_method,
                        min_periods=1,
                        pws_timegate=pws_timegate,
                        pws_power=pws_power,
                        goal_sampling_rate=params.cc.cc_sampling_rate,
                        freqmin=params.filter.freqmin,
                        freqmax=params.filter.freqmax,
                    )
                    xx = xr.Dataset(
                        {"CCF": (dr["CCF"].dims, rolled)},
                        coords=dr["CCF"].coords,
                    )
                    xx = xx.resample(times=mov_sample, label="right", skipna=True).asfreq().dropna("times", how="all")

                logger.debug(f"  [timing] rolling/resample {mov_stack}: {time.time()-_t0:.2f}s  shape={xx['CCF'].shape}")
                mask = xx.times.dt.floor('D').isin(excess_dates)
                xx_cleaned = xx.where(~mask, drop=True)

                _t0 = time.time()
                xr_save_ccf(params.global_.output_folder, lineage_names, step.step_name,
                            sta1, sta2, components, mov_stack, taxis, xx_cleaned, overwrite=False)
                logger.debug(f"  [timing] save {mov_stack}: {time.time()-_t0:.2f}s")
                del xx, xx_cleaned

            del dr          # free resampled CCF dataset before next component
            del excess_dates, all_days

        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)

        # legacy hpc comment (pre-P48):
        # if stype != "step" and not params.global_.hpc: update MWCS/WCT/STR
