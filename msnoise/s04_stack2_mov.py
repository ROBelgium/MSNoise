"""
TODO
"""

import argparse

import scipy.signal

from .api import *
from .wiener import *
from obspy.core import AttribDict

import logbook
import matplotlib.pyplot as plt


def main(stype, loglevel="INFO"):
    """Computes the REF/MOV stacks.

    Parameters
    ----------
    stype : {'mov', 'ref'}
        Defines which of the REF or Moving-window stacks must be exported
    interval : float, optional
        Number of days before now to search for modified CC jobs

    """
    # logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.stack_child', loglevel, with_pid=True)
    logger.debug('Starting the %s stack' % stype)
    db = connect()

    orig_params = get_params(db)
    taxis = get_t_axis(db)

    while is_next_job_for_step(db, step_category="stack"):
        logger.info("Getting the next job")
        jobs, step = get_next_job_for_step(db, step_category="stack", group_by="pair_lineage")

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        lineage_str = jobs[0].lineage
        if not lineage_str:
            raise ValueError("STACK jobs must have a non-empty lineage (v2 assumption)")

        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        # 1) current step config (stack_*)
        step_params = get_config_set_details(
            db,
            jobs[0].config_category,
            jobs[0].config_set_number,
            format="AttribDict"
        )

        # 2) resolve lineage steps from the job, and merge params for THIS lineage only
        lineage_steps = lineage_str_to_steps(db, lineage_str, workflow_id=jobs[0].workflow_id, strict=True)

        lineage_steps, lineage_names, params = get_merged_params_for_lineage(
            db, orig_params, step_params, lineage_steps
        )

        # drop the "current step name" if you use lineage_names as folder segments upstream of stack output
        lineage_names = lineage_names[:-1]

        logger.info(f"New STACK Job: pair={pair} n_days={len(days)} lineage={lineage_str}")



        mov_stacks = params.mov_stack
        wiener_mlen = params.wiener_mlen
        wiener_nlen = params.wiener_nlen
        wienerfilt = params.wienerfilt
        wiener_M = int(pd.to_timedelta(wiener_mlen).total_seconds() / params.corr_duration)
        wiener_N = int(pd.to_timedelta(wiener_nlen).total_seconds() * params.cc_sampling_rate)

        # is there a better alternative for threshold?
        if params.keep_all:
            wiener_gap_threshold = wiener_M  # no. indices which will be considered adjacent by wiener
        else:
            wiener_gap_threshold = pd.to_timedelta(wiener_mlen).days

        if wienerfilt:
            logger.info('Wiener filter enabled, will apply to CCFs before stacking')

        logger.info(
            "There are STACKS jobs for some days to recompute for %s" % pair)

        sta1, sta2 = pair.split(':')
        filterid = 1
        for components in params.components_to_compute:
            logger.info('Processing %s-%s-%i MOV stack' %
                (pair, components, filterid))

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

            if params.keep_all:
                c = get_results_all(db, params.output_folder, lineage_names,
                                    sta1, sta2, components, all_days, format="xarray")
                if not len(c):
                    logger.warning("No data found for %s-%s-%i" % (sta1, sta2, filterid))
                    continue
                c = c.sortby('times')
                dr = c.resample(times="%is" % params.corr_duration).mean()

            else:
                # TODO this is not yet compatible with the lineage folder structure:
                logger.warning("keep_all=N used, sampling interval will be 1-day")
                c = get_results(db, lineage_names, sta1, sta2, filterid, components, all_days,  mov_stack=1, format="xarray",
                    params=params).sortby('times')
                dr = c.resample(times="1D").mean()

            is_valid, message = validate_stack_data(c, "moving")
            if not is_valid:
                logger.error(f"Invalid moving stack data for {sta1}:{sta2}-{components}-{filterid}: {message}")
                continue
            elif "Warning" in message:
                logger.warning(f"{sta1}:{sta2}-{components}-{filterid}: {message}")

            if wienerfilt:
                dr = wiener_filt(dr, wiener_M, wiener_N, wiener_gap_threshold)

            excess_dates = pd.to_datetime(excess_days).values
            for mov_stack in mov_stacks:
                # if mov_stack > len(dr.times):
                #     logger.error("not enough data for mov_stack=%i" % mov_stack)
                #     continue
                mov_rolling, mov_sample = mov_stack
                # print(mov_rolling, mov_sample)

                if mov_rolling == mov_sample:
                    # just resample & mean
                    xx = dr.resample(times=mov_sample, label="right", skipna=True).mean().dropna("times",
                                                                                                   how="all")
                else:
                    mov_rolling = pd.to_timedelta(mov_rolling).total_seconds()
                    # print("Will roll over %i seconds" % mov_rolling)
                    if params.keep_all:
                        duration_to_windows = mov_rolling / params.corr_duration
                    else:
                        duration_to_windows = mov_rolling / 86400.0
                    if not duration_to_windows.is_integer():
                        logger.print("Warning, rounding down the number of windows to roll over")
                    duration_to_windows = int(max(1, math.floor(duration_to_windows)))
                    # print("Which is %i windows of %i seconds duration" % (duration_to_windows, params.corr_duration))
                    xx = dr.rolling(times=duration_to_windows, min_periods=1).mean("win")
                    xx = xx.resample(times=mov_sample, label="right", skipna=True).asfreq().dropna("times", how="all")

                mask = xx.times.dt.floor('D').isin(excess_dates)
                xx_cleaned = xx.where(~mask, drop=True) #remove days not associated with current jobs

                xr_save_ccf(params.output_folder, lineage_names, step.step_name,
                            sta1, sta2, components, filterid, mov_stack, taxis, xx_cleaned, overwrite=False)
                del xx, xx_cleaned

        massive_update_job(db, jobs, "D")

        # if stype != "step" and not params.hpc:
        #     for job in jobs:
        #         update_job(db, job.day, job.pair, 'MWCS', 'T')
        #         update_job(db, job.day, job.pair, 'WCT', 'T')
        #         update_job(db, job.day, job.pair, 'STR', 'T')
