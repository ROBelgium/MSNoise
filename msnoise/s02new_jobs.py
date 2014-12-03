""" This script searches the database for files flagged "N"ew or "M"odified.
For each date in the configured range, it checks if other stations are
available and defines the new jobs to be processed. Those are inserted in the
*jobs* table of the database.

To run it from the console:

.. code-block:: sh

    $ msnoise new_jobs
"""

from database_tools import *
import logging
import numpy as np
import os


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: New Jobs ***')

    db = connect()
    if get_config(db, name="autocorr") in ['Y', 'y', '1', 1]:
        AUTOCORR = True
    else:
        AUTOCORR = False

    stations_to_analyse = [sta.sta for sta in get_stations(db, all=False)]
    nfs = get_new_files(db)

    days = {}
    old_day = 0
    old_pair = ""
    day_pairs = []
    jobs = []
    all_jobs = []
    i = 0
    for nf in nfs:
        # logging.debug('%s.%s will be MASTER for %s-%s'% (nf.net, nf.sta, nf.starttime, nf.endtime))
        if nf.sta in stations_to_analyse:
            day = "%s" % (nf.starttime.date())
    
            available_stations = []
            for station in get_data_availability(db, starttime=nf.starttime, endtime=nf.endtime):
                if station.sta in stations_to_analyse:
                    if '%s.%s' % (station.net, station.sta) not in available_stations:
                        available_stations.append(
                            '%s.%s' % (station.net, station.sta))
                    # logging.debug('Will process %s.%s vs %s.%s'% (nf.net, nf.sta, station.net, station.sta))
    
            stations = np.array([])
            pairs = []
            nS = '%s.%s' % (nf.net, nf.sta)
            i = 0
            for aS in available_stations:
                if not AUTOCORR and nS == aS:
                    pass
                else:
                    if i == 0:
                        pairs = np.array(':'.join(sorted([nS, aS])))
                        i += 1
                    else:
                        pairs = np.vstack((pairs, ':'.join(sorted([nS, aS]))))
    
            pairs = np.unique(pairs)
            for pair in pairs:
                daypair = "%s=%s" % (day, pair)
                if daypair not in jobs:
                    all_jobs.append({"day":day,"pair":pair,"type":"CC","flag":"T"})
                    jobs.append(daypair)
    
    count = len(jobs)
    logging.debug("Found %i new jobs to do" % count)
    
    # for job in jobs:
        # day, pair = job.split("=")
        # alljobs.append({"day":day,"pair":pair,"type":"CC","flag":"T"})

    massive_insert_job(all_jobs)
    # update all _data_availability and mark files as "A"rchives
    for sta in get_stations(db, all=True):
        mark_data_availability(db, sta.net, sta.sta, flag='A')
    
    logging.info('*** Finished: New Jobs ***')
    return count


if __name__ == "__main__":
    main()