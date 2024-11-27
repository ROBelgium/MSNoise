"""MSNoise is capable of using a reference function defined by absolute or
relative dates span. For example, an absolute range could be "from 1 January
2010 to 31 December 2011" and a relative range could be "the last 200 days".
In the latter case, the REF will need to be exported at every run, meaning the
following steps (MWCS and DTT) will be executed on the whole configured period.
If the REF is defined between absolute dates, excluding "today", the MWCS and
DTT will only be calculated for new data (e.g. "yesterday" and "today").
The corresponding configuration bits are ``ref_begin`` and ``ref_end``. In the
future, we plan on allowing multiple references to be defined.

Only data for new/modified dates need to be exported. If any CC-job has been
marked "Done" within the last day and triggered the creation of STACK jobs,
the stacks will be calculated and a new MWCS job will be inserted in the
database. For dates in the period of interest, the moving-window stack will
only be exported if new/modified CCF is available.
The export directory are "REF/" and "DAY%03i/" where %03i will be replaced by
the number of days stacked together (DAYS_005 for a 5-days stack, e.g.).

Please note that within MSNoise, stacks are always *inclusive* of the time/day
mentioned. For example, a 5-days stack on January 10, will contain
cross-correlation functions computed for January 6, 7, 8, 9 AND 10!
The graphical representation centered on a "January 10" tick might then display
changes in the CCF that occurred *on* the 10th !

Moving-window stacks are configured using the ``mov_stack`` parameter in
``msnoise admin``.

If ``stack_method`` is 'linear', then a simple mean CFF of all daily is saved
as the mov or ref CCF. On the other hand, if ``stack_method`` is 'pws', then
all the Phase Weighted Stack (PWS) is computed and saved as the mov or ref CCF.
The PWS is done in two steps: first the mean coherence between the instantaneous
phases of all windows is calculated, and eventually serves a weighting factor
on the mean. The smoothness of this weighting array is defined using the
``pws_timegate`` parameter in the configuration. The weighting array is the
power of the mean coherence array. If ``pws_power`` is equal to 0, a linear
stack is done (then it's faster to do set ``stack_method`` = 'linear'). Usual
value is 2.

.. warning:: PWS is largely untested, not cross-validated. It looks good, but
    that doesn't mean a lot, does it? Use with Caution! And if you
    cross-validate it, please let us know!!

    Schimmel, M. and Paulssen H., "Noise reduction and detection
    of weak, coherent signals through phase-weighted stacks". Geophysical Journal
    International 130, 2 (1997): 497-505.

If ``wienerfilt`` set to 'Y' (default: 'N'), a Wiener filter is applied to CCFs
prior to stacking. This applies smoothing to CCFs in both datetime axis 
(window length corresponding to ``wiener_mlen``) and lag time axis (window 
length corresponding to ``wiener_nlen``). For applying the filter, CCFs are 
treated as continguous or separate segments depending on the temporal gap between 
them. For gaps shorter than ``wiener_mlen``, segments are treated as adjacent 
by the filter. For larger gaps, the segments are considered seperate and the 
filter is applied independently to each segment (thus avoiding smoothing over 
large temporal gaps). 

To reduce influence of edge effects in wiener filter, CCF segments are padded 
in datetime axis with an additional window of 2 * ``wiener_mlen`` CCFs (rounded 
up to nearest day, if they exist) at both ends. After filtering, the outer 
half of the padding window (``wiener_mlen``, farthest from original segment) is 
discarded. The inner half (cloeset to original segment) will (after stacking) 
overwrite previously stacked data. This ensures that data near the edges of a
segment, which may have initially been smoothed with limited neighbouring CCFs
(e.g. if processing data in real-time), are updated when additional data becomes 
available.

.. warning:: While padding "T"odo jobs with existing data (where possible) may 
    reduce the influence of edge effects of wiener filter, caution should be 
    exercised when processing data in separate steps/segments (e.g. real-time 
    processing of smaller data segments) as new data can produce changes further 
    away from the edge of the CCF 'image'. Applying the wiener filter may, 
    therefore, be better suited when possible to process all jobs at once 
    (i.e. set all STACK jobs as "T"odo).


Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |ref_begin|
* |ref_end|
* |mov_stack| | * changed in 2.0 !!
* |stack_method| | *new in 1.4*
* |pws_timegate| | *new in 1.4*
* |pws_power| | *new in 1.4*
* |wienerfilt| | *new in 2.0*
* |wiener_mlen| | *new in 2.0*
* |wiener_nlen| | *new in 2.0*
* |hpc| | *new in 1.6*


Once done, each job is marked "D"one in the database and, unless ``hpc`` is
``Y``, MWCS jobs are inserted/updated in the database.

Usage:
~~~~~~

.. include:: ../clickhelp/msnoise-cc-stack.rst


For most users, the REF stack will need to be computed only once for specific
dates and then, on routine basis, only compute the MOV stacks:


.. Computing the REF and MOV stacks is done as two separate actions. From MSNoise 2.0, 
    the REF stack is computed independant of the status (flag) of STACK jobs. 
    Computing the REF stack will therefore no longer set STACK jobs to "I"n
    progress, previously requiring them to be reset prior to running the MOV stacks.

.. code-block:: sh

    $ msnoise stack -r
    $ msnoise stack -m

as for all other steps, this procedure can be run in parallel:

.. code-block:: sh

    $ msnoise -t 4 stack -r
    $ msnoise -t 4 stack -m


.. versionadded:: 1.4
    The Phase Weighted Stack.

.. versionadded:: 1.6
    The ``hpc`` parameter that can prevent the automatic creation of MWCS jobs.
    The REF and MOV stacks have been separated and need to be run independently.

.. versionadded:: 2.0
    MOV stack defined as list of tuples (allowing greater flexibility in stack length/sample rate, including subdaily).
    The REF stack computation no longer depends on STACK jobs table.
    Wiener filter option added pre-stacking.
    
"""

import argparse

import scipy.signal

from .api import *
from .wiener import *


import logbook
import matplotlib.pyplot as plt

def main(stype, interval=1.0, loglevel="INFO"):
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
    logger = get_logger('msnoise.stack_child', loglevel,
                        with_pid=True)
    logger.debug('Starting the %s stack' % stype)
    db = connect()

    params = get_params(db)
    taxis = get_t_axis(db)

    mov_stacks = params.mov_stack
    # if 1 in mov_stacks:
    #     mov_stacks.remove(1)  # remove 1 day stack, it will be done automatically

    filters = get_filters(db, all=False)

    wiener_mlen = params.wiener_mlen
    wiener_nlen = params.wiener_nlen
    wienerfilt = params.wienerfilt
    wiener_M =  int(pd.to_timedelta(wiener_mlen).total_seconds() / params.corr_duration)
    wiener_N =  int(pd.to_timedelta(wiener_nlen).total_seconds() * params.cc_sampling_rate)
    
    # is there a better alternative for threshold?
    wiener_gap_threshold = wiener_M #no. indices which will be considered adjacent by wiener
    
    if wienerfilt:
        logger.info('Wiener filter enabled, will apply to CCFs before stacking')

    if stype == "ref":
        
        pairs = []
        for sta1, sta2 in get_station_pairs(db, used=False):
            for loc1 in sta1.locs():
                s1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
                for loc2 in sta2.locs():
                    s2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                    pairs.append((s1, s2))
   
        for sta1, sta2 in pairs:
            for f in filters:
                filterid = int(f.ref)

                if sta1 == sta2:
                    components_to_compute = params.components_to_compute_single_station
                else:
                    components_to_compute = params.components_to_compute

                for components in components_to_compute:
                    logger.info('Processing %s:%s-%s-%i REF stack' %
                        (sta1, sta2, components, filterid))

                    start, end, datelist = build_ref_datelist(db)
                    
                    if params.keep_all:
                        c = get_results_all(db, sta1, sta2, filterid, components, datelist, format="xarray")
                    else:
                        c = get_results(db, sta1, sta2, filterid, components, datelist,  mov_stack=1, format="xarray", params=params)

                    # dr = xr_save_ccf(sta1, sta2, components, filterid, 1, taxis, c)
                    dr = c
                    if not c.data_vars:
                        logger.debug('No data found for %s:%s-%s-%i' %
                            (sta1, sta2, components, filterid))
                        continue
                    
                    if wienerfilt:
                        dr = wiener_filt(dr, wiener_M, wiener_N, wiener_gap_threshold)

                    # TODO add other stack methods here! using apply?
                    _ = dr.mean(dim="times")                    
                    xr_save_ref(sta1, sta2, components, filterid, taxis, _)
    
    else:   #stype== 'mov'
        while is_dtt_next_job(db, flag='T', jobtype='STACK'):
            jobs = get_dtt_next_job(db, flag='T', jobtype='STACK')

            if not len(jobs):
                # edge case, should only occur when is_next returns true, but
                # get_next receives no jobs (heavily parallelised calls).
                time.sleep(np.random.random())
                continue
            pair = jobs[0].pair
            refs, days = zip(*[[job.ref, job.day] for job in jobs])

            logger.info(
                "There are STACKS jobs for some days to recompute for %s" % pair)

            sta1, sta2 = pair.split(':')
            for f in filters:
                filterid = int(f.ref)

                if sta1 == sta2:
                    components_to_compute = params.components_to_compute_single_station
                else:
                    components_to_compute = params.components_to_compute

                for components in components_to_compute:
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
                        c = get_results_all(db, sta1, sta2, filterid, components, all_days, format="xarray")
                    else:
                        logger.warning("keep_all=N used by default mov_stack=('1D','1D')")
                        c = get_results(db, sta1, sta2, filterid, components, all_days,  mov_stack=1, format="xarray", params=params)
                    
                    dr = c
                    dr = dr.sortby('times')
                    dr = dr.resample(times="%is" % params.corr_duration).mean()

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
                            duration_to_windows = mov_rolling / params.corr_duration
                            if not duration_to_windows.is_integer():
                                logger.print("Warning, rounding down the number of windows to roll over")
                            duration_to_windows = int(max(1, math.floor(duration_to_windows)))
                            # print("Which is %i windows of %i seconds duration" % (duration_to_windows, params.corr_duration))
                            # That construct thing shouldn't be necessary, but without it, tests fail when bottleneck is installed
                            # ref: https://github.com/pydata/xarray/issues/3165
                            xx = dr.rolling(times=duration_to_windows, min_periods=1).construct("win").mean("win")
                            xx = xx.resample(times=mov_sample, label="right", skipna=True).asfreq().dropna("times", how="all")
                       
                        mask = xx.times.dt.floor('D').isin(excess_dates)
                        xx_cleaned = xx.where(~mask, drop=True) #remove days not associated with current jobs

                        xr_save_ccf(sta1, sta2, components, filterid, mov_stack, taxis, xx_cleaned, overwrite=False)
                        del xx, xx_cleaned

            massive_update_job(db, jobs, "D")
            if stype != "step" and not params.hpc:
                for job in jobs:
                    update_job(db, job.day, job.pair, 'MWCS', 'T')
