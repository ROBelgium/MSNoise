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
mentionned. For example, a 5-days stack on January 10, will contain
cross-correlation functions computed for January 6, 7, 8, 9 AND 10!
The graphical representation centered on a "January 10" tick might then display
changes in the CCF that occurred *on* the 10th !

Moving-window stack length(s) are configured using the ``mov_stack`` bit.

If ``stack_method`` is 'linear', then a simple mean CFF of all daily is saved
as the mov or ref CCF. On the other hand, if ``stack_method`` is 'pws', then
all the Phase Weighted Stack (PWS) is computed and saved as the mov or ref CCF.
The PWS is done in two steps: first the mean coherence between the instataneous
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

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |ref_begin|
* |ref_end|
* |mov_stack|
* |stack_method| | *new in 1.4*
* |pws_timegate| | *new in 1.4*
* |pws_power| | *new in 1.4*
* |hpc| | *new in 1.6*


Once done, each job is marked "D"one in the database and, unless ``hpc`` is 
``Y``, MWCS jobs are inserted/updated in the database.

Usage:
~~~~~~
The best way to call this code is to start it from the console (-h shows the
help)

.. include:: ../clickhelp/msnoise-stack.rst


On a routine basis, one should thus run the following to compute REF *and* MOV
stacks:

.. todo:: finalize the difference between REF and MOV/STEP stacks.

.. code-block:: sh

    $ msnoise stack -r -m


.. versionadded:: 1.4
    The Phase Weighted Stack.

.. versionadded:: 1.6
    The ``hpc`` parameter that can prevent the automatic creation of MWCS jobs.
"""

import argparse

import scipy.signal

from .api import *


import logbook


def main(stype, interval=1.0, loglevel="INFO"):
    """Computes the REF/MOV stacks.
    
    Parameters
    ----------
    stype : {'mov', 'ref'}
        Defines which of the REF or Moving-window stacks must be exported
    interval : float, optional
        Number of days before now to search for modified CC jobs

    """
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.stack_child', loglevel,
                        with_pid=True)
    logger.debug('Starting the %s stack' % stype)
    db = connect()
    
    export_format = get_config(db, 'export_format')

    if export_format == "BOTH":
        mseed = True
        sac = True
    elif export_format == "SAC":
        mseed = False
        sac = True
    elif export_format == "MSEED":
        mseed = True
        sac = False
    
    maxlag = float(get_config(db, "maxlag"))
    cc_sampling_rate = float(get_config(db, "cc_sampling_rate"))

    stack_method = get_config(db, 'stack_method')
    pws_timegate = float(get_config(db, 'pws_timegate'))
    pws_power = float(get_config(db, 'pws_power'))
    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    # Get Configuration
    params = get_params(db)
    plugins = get_config(db, "plugins")
    extra_jobtypes = []
    if plugins:
        plugins = plugins.split(",")
        for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.jobtypes'):
            module_name = ep.module_name.split(".")[0]
            if module_name in plugins:
                jobtypes = ep.load()()
                for jobtype in jobtypes:
                    if jobtype["after"] == "refstack":
                        extra_jobtypes.append(jobtype["name"])
    
    if stype == "mov" or stype == "step":
        start, end, datelist = build_movstack_datelist(db)
        format = "matrix"
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]
        if 1 in mov_stacks:
            mov_stacks.remove(1)  # remove 1 day stack, it should exist already
    
    elif stype == "ref":
        start, end, datelist = build_ref_datelist(db)
        format = "stack"
    
    if stype == "step":
        datelists = {}
        for mov_stack in mov_stacks:
            if mov_stack == 7:
                rng = pd.date_range(start, end, freq="W")
            elif mov_stack == 31:
                rng = pd.date_range(start, end, freq="M")
            elif mov_stack == 91:
                rng = pd.date_range(start, end, freq="Q")
            else:
                rng = pd.date_range(start, end, freq="%iD"%mov_stack)
            datelists[mov_stack] = rng.map(lambda t: t.date())
        #~ print datelists
    biglist = []
    filters = get_filters(db, all=False)
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
            for components in params.all_components:
                pair = "%s:%s" % (sta1, sta2)
                sta1 = sta1.replace('.', '_')
                sta2 = sta2.replace('.', '_')
                logger.debug('Processing %s-%s-%i' %
                                  (pair, components, filterid))
                # updated_days = updated_days_for_dates(db, start, end, pair.replace('_', '.'), jobtype='CC', interval=datetime.timedelta(days=interval),returndays=True)
                updated_days = [UTCDateTime(d).datetime.date() for d in days]
                if len(updated_days) != 0:
                    logger.debug("New Data for %s-%s-%i" %
                                  (pair, components, filterid))
                    #~ print updated_days
                    # TODO: load only the updated dates +- max(mov_stack)+1
                    # Note: this would no longer be needed if the stack is h5
                    nstack, stack_total = get_results(
                        db, sta1, sta2, filterid, components, datelist, format=format, params=params)
                    if not nstack:
                        logger.debug("No new data found, hmmm")
                    logger.debug("Data loaded")
                    if nstack > 0:
                        if stype == "mov":
                            logger.debug("Mov Stack!")
                            for i, date in enumerate(datelist):
                                jobadded = False
                                for mov_stack in mov_stacks:
                                    if i < mov_stack:
                                        low = 0
                                        high = mov_stack
                                    else:
                                        low = i - mov_stack + 1
                                        high = i + 1
                                    newdata = False
                                    for uday in datelist[low:high]:
                                        if uday in updated_days:
                                            newdata = True
                                            break
                                    if newdata:
                                        corr = stack_total[low:high]
                                        if not np.all(np.isnan(corr)):
                                            day_name = "%s_%s" % (
                                                sta1, sta2)
                                            logger.debug("%s %s %s [%s - %s] (%i day stack)" % (
                                                day_name, components, date, datelist[low], datelist[i], mov_stack))
                                            corr = stack(corr, stack_method,
                                                         pws_timegate,
                                                         pws_power,
                                                         goal_sampling_rate)
                                            if not len(corr):
                                                continue

                                            corr = scipy.signal.detrend(
                                                corr).astype(np.float32)
                                            stack_path = os.path.join(
                                                "STACKS", "%02i" % filterid, "%03i_DAYS" % mov_stack, components, day_name)
                                            filename = os.path.join(
                                                stack_path, str(date))
                                            if mseed:
                                                export_mseed(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate, params=params)
                                            if sac:
                                                export_sac(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate, params=params)
                                            day_name = "%s:%s" % (
                                                sta1, sta2)
                                            if not jobadded and not params.hpc:
                                                update_job(
                                                    db, date, day_name.replace('_', '.'), 'MWCS', 'T')
                                                jobadded = True
                                        del corr
                        elif stype == "step":
                            jobs_done = []
                            for mov_stack in mov_stacks:
                                for i, date in enumerate(datelists[mov_stack]):
                                    if date not in datelist:
                                        continue
                                    if i < mov_stack:
                                        low = 0
                                        high = mov_stack
                                    else:
                                        low = datelist.index(date) - mov_stack + 1
                                        high = datelist.index(date) + 1
                                    newdata = False
                                    for uday in datelist[low:high]:
                                        if uday in updated_days:
                                            newdata = True
                                            break
                                    if newdata:
                                        corr = stack_total[low:high]
                                        if not np.all(np.isnan(corr)):
                                            day_name = "%s_%s" % (
                                                sta1, sta2)
                                            logger.debug("%s %s %s [%s - %s] (%i day stack)" % (
                                                day_name, components, date, datelist[low], datelist[high-1], mov_stack))
                                            corr = stack(corr, stack_method,
                                                         pws_timegate,
                                                         pws_power,
                                                         goal_sampling_rate)
                                            corr = scipy.signal.detrend(corr)
                                            stack_path = os.path.join(
                                                "STACKS", "%02i" % filterid, "%03i_DAYS" % mov_stack, components, day_name)
                                            filename = os.path.join(
                                                stack_path, str(date))
                                            if mseed:
                                                export_mseed(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate, params=params)
                                            if sac:
                                                export_sac(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate, params=params)
                                            day_name = "%s:%s" % (
                                                sta1, sta2)
                                            job = "%s %s" % (date, day_name)
                                            if job not in jobs_done and not params.hpc:
                                                update_job(
                                                    db, date, day_name.replace('_', '.'), 'MWCS', 'T')
                                                jobs_done.append(job)
                                        del corr

                        elif stype == "ref":
                            stack_path = os.path.join(
                                "STACKS", "%02i" % filterid, "REF", components)
                            ref_name = "%s_%s" % (sta1, sta2)
                            filename = os.path.join(stack_path, ref_name)
                            stack_total = scipy.signal.detrend(stack_total)

                            if mseed:
                                export_mseed(
                                    db, filename, pair, components, filterid, stack_total, params=params)
                            if sac:
                                export_sac(
                                    db, filename, pair, components, filterid, stack_total,params=params)
                            ref_name = "%s:%s" % (sta1, sta2)
                            update_job(
                                db, "REF", ref_name.replace('_', '.'), 'MWCS', 'T')
                            for jobtype in extra_jobtypes:
                                update_job(db, "REF", ref_name.replace('_', '.'), jobtype, 'T')
                            del stack_total

        # THIS SHOULD BE IN THE API
        # This doesn't set MWCS jobs for REF stacks
        if stype != "ref":
            massive_update_job(db, jobs, "D")
            if stype != "step" and not params.hpc:
                for job in jobs:
                    update_job(db, job.day, job.pair, 'MWCS', 'T')
        if stype == "ref":
            biglist += jobs

    if stype == "ref":
        logger.info("You just finished REF stacking, remember to reset the "
                    "STACK jobs if you need to compute a MOV stacks. "
                    "Run 'msnoise reset STACK' when all process have finished.")
        logger.info("The current STACK jobs have been intentionnaly left "
                    "'I'n progress so they can be reset.")
    #     massive_update_job(db, biglist, "T")

    logger.debug("Finished Stacking")


def refstack(interval):
    main('ref', interval)


def movstack(interval):
    main('mov', interval)
