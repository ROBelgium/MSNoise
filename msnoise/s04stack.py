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
marked "Done" within the last day, the stacks will be calculated and a new DTT
job will be inserted in the database. For dates in the period of interest, the
moving-window stack will only be exported if new/modified CCF is available.
The export directory are "REF/" and "DAY%03i/" where %03i will be replaced by
the number of days stacked together (DAYS_005 for a 5-days stack, e.g.).

Please note that within MSNoise, stacks are always *inclusive* of the time/day
mentionned. For example, a 5-days stack on January 10, will contain
cross-correlation functions computed for January 6, 7, 8, 9 AND 10!
The graphical representation centered on a "January 10" tick might then display
changes in the CCF that occurred *on* the 10th !

Moving-window stack length(s) are configured using the ``mov_stack`` bit.

Usage:
~~~~~~
The best way to call this code is to start it from the console (-h shows the
help)

.. code-block:: sh

    $ msnoise stack --help

    Usage: msnoise-script.py stack [OPTIONS]

      Stacks the [REF] and/or [MOV] windows

    Options:
      -r, --ref               Compute the REF Stack
      -m, --mov               Compute the MOV Stacks
      -s, --step              Compute the STEP Stacks
      -i, --interval INTEGER  Number of days before now to search for modified
                              Jobs
      --help                  Show this message and exit.

On a routine basis, one should thus run the following to compute REF *and* MOV
stacks:

.. code-block:: sh

    $ msnoise stack -r -m

While, when playing around with data, and surely on the first run, one should
define the *-i INTERVAL*, as jobs might have been marked "Done" more than 24
hours before running the stack. This, for example, will tell the code to search
for jobs marked in the last 10 days:

.. code-block:: sh

    $ msnoise stack -r -m -i 10
"""

import os
import argparse
import numpy as np
import pandas as pd

import scipy.signal
from numpy import nanmean
from api import *

import logging


def main(stype, interval=1):
    """Computes the REF/MOV stacks.
    
    Parameters
    ----------
    stype : {'mov', 'ref'}
        Defines which of the REF or Moving-window stacks must be exported
    interval : int, optional
        Number of days before now to search for modified CC jobs

    """
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

    
    logging.debug('Starting the %s stack' % stype)
    db = connect()
    components_to_compute = get_components_to_compute(db)
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

    for f in get_filters(db, all=False):
        filterid = int(f.ref)
        for components in components_to_compute:
            for station1, station2 in get_station_pairs(db, used=True):
                sta1 = "%s_%s" % (station1.net, station1.sta)
                sta2 = "%s_%s" % (station2.net, station2.sta)
                pair = "%s:%s" % (sta1, sta2)
                logging.debug('Processing %s-%s-%i' %
                                  (pair, components, filterid))
                updated_days = updated_days_for_dates(db, start, end, pair.replace('_', '.'), jobtype='CC', interval=datetime.timedelta(days=interval),returndays=True)
                if len(updated_days) != 0:
                    logging.debug("New Data for %s-%s-%i" %
                                  (pair, components, filterid))
                    #~ print updated_days
                    nstack, stack_total = get_results(
                        db, sta1, sta2, filterid, components, datelist, format=format)
                    if nstack > 0:
                        if stype == "mov":
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
                                            logging.debug("%s %s [%s - %s] (%i day stack)" % (
                                                day_name, date, datelist[low], datelist[i], mov_stack))
                                            corr = nanmean(corr, axis=0)
                                            corr = scipy.signal.detrend(
                                                corr)
                                            stack_path = os.path.join(
                                                "STACKS", "%02i" % filterid, "%03i_DAYS" % mov_stack, components, day_name)
                                            filename = os.path.join(
                                                stack_path, str(date))
                                            if mseed:
                                                export_mseed(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate)
                                            if sac:
                                                export_sac(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate)
                                            day_name = "%s:%s" % (
                                                sta1, sta2)
                                            if not jobadded:
                                                update_job(
                                                    db, date, day_name.replace('_', '.'), 'DTT', 'T')
                                                jobadded = True
                                        del corr
                        elif stype == "step":
                            jobs = []
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
                                            logging.debug("%s %s [%s - %s] (%i day stack)" % (
                                                day_name, date, datelist[low], datelist[high-1], mov_stack))
                                            corr = nanmean(corr, axis=0)
                                            corr = scipy.signal.detrend(
                                                corr)
                                            stack_path = os.path.join(
                                                "STACKS", "%02i" % filterid, "%03i_DAYS" % mov_stack, components, day_name)
                                            filename = os.path.join(
                                                stack_path, str(date))
                                            if mseed:
                                                export_mseed(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate)
                                            if sac:
                                                export_sac(
                                                    db, filename, pair, components, filterid, corr, maxlag=maxlag, cc_sampling_rate=cc_sampling_rate)
                                            day_name = "%s:%s" % (
                                                sta1, sta2)
                                            job = "%s %s" % (date, day_name)
                                            if job not in jobs:
                                                update_job(
                                                    db, date, day_name.replace('_', '.'), 'DTT', 'T')
                                                jobs.append(job)
                                        del corr
                            #~ for date in datelist:
                                #~ day_name = "%s:%s" % (
                                                    #~ sta1, sta2)
                                #~ job = "%s %s" % (date, day_name)
                                #~ if job not in jobs:
                                    #~ update_job(
                                        #~ db, date, day_name.replace('_', '.'), 'DTT', 'T')
                                    #~ jobs.append(job)
                        
                        elif stype == "ref":
                            stack_path = os.path.join(
                                "STACKS", "%02i" % filterid, "REF", components)
                            ref_name = "%s_%s" % (sta1, sta2)
                            filename = os.path.join(stack_path, ref_name)
                            stack_total = scipy.signal.detrend(stack_total)

                            if mseed:
                                export_mseed(
                                    db, filename, pair, components, filterid, stack_total)
                            if sac:
                                export_sac(
                                    db, filename, pair, components, filterid, stack_total)
                            ref_name = "%s:%s" % (sta1, sta2)
                            update_job(
                                db, "REF", ref_name.replace('_', '.'), 'DTT', 'T')
                            del stack_total


def refstack(interval):
    main('ref', interval)


def movstack(interval):
    main('mov', interval)


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

    parser = argparse.ArgumentParser(description='Compute [REF,MOV] stacks if\
                                jobs have been modified in the last i days.',
                                     epilog=__doc__)
    parser.add_argument('-r', '--ref', action="store_true",
                        help='Triggers the computation of REF stacks',
                        default=False)
    parser.add_argument('-m', '--mov', action="store_true",
                        help='Triggers the computation of MOV stacks',
                        default=False)
    parser.add_argument('-i', '--interval',
                        help='Number of days before now to search for\
                        modified CC jobs [default:1]', default=1, type=int)
    args = parser.parse_args()

    logging.info('Starting this program with: ref=%s, mov=%s, interval=%i'
                 % (args.ref, args.mov, args.interval))

    db = connect()
    if args.ref:
        logging.info("*** Starting: REF Stack ***")
        refstack(args.interval)
        logging.info("*** Finished: REF Stack ***")
    if args.mov:
        logging.info("*** Starting: MOV Stack ***")
        movstack(args.interval)
        logging.info("*** Finished: MOV Stack ***")
