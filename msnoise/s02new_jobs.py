""" This script searches the database for files flagged "N"ew or "M"odified.
For each date in the configured range, it checks if other stations are
available and defines the new jobs to be processed. Those are inserted in the
*jobs* table of the database.

To run it from the console:

.. code-block:: sh

    $ msnoise new_jobs

.. warning:: Upon first run, if you expect the number of jobs to be large (many
    days, many stations), pass the ``--init`` parameter to optimize the insert.
    Only use this flag once, otherwise problems will arise from duplicate
    entries in the jobs table.
"""

from .api import *
import pandas as pd


def main(init=False, nocc=False):

    logging.info('*** Starting: New Jobs ***')

    db = connect()

    logging.debug("Checking plugins' entry points")
    plugins = get_config(db, "plugins")
    extra_jobtypes_scan_archive = []
    extra_jobtypes_new_files = []
    if plugins:
        plugins = plugins.split(",")
        for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.jobtypes'):
            module_name = ep.module_name.split(".")[0]
            if module_name in plugins:
                jobtypes = ep.load()()
                for jobtype in jobtypes:
                    if jobtype["after"] == "scan_archive":
                        extra_jobtypes_scan_archive.append(jobtype["name"])
                    elif jobtype["after"] == "new_files":
                        extra_jobtypes_new_files.append(jobtype["name"])

    autocorr = get_config(db, name="autocorr", isbool=True)
    logging.debug('Scanning New/Modified files')
    stations_to_analyse = ["%s.%s" % (sta.net, sta.sta) for sta in get_stations(db, all=False)]
    all_jobs = []
    updated_days = []
    nfs = get_new_files(db)
    now = datetime.datetime.utcnow()
    for nf in nfs:
        tmp = "%s.%s" % (nf.net, nf.sta)
        if tmp not in stations_to_analyse:
            continue

        start, end = nf.starttime.date(), nf.endtime.date()
        for date in pd.date_range(start, end, freq="D"):
            updated_days.append(date.date())
            for jobtype in extra_jobtypes_new_files:
                all_jobs.append({"day": date.date(),
                                 "pair": "%s.%s" % (nf.net, nf.sta),
                                 "jobtype": jobtype,
                                 "flag": "T", "lastmod": now})

    all_jobs = list(np.unique(all_jobs))
    updated_days = np.asarray(updated_days)
    updated_days = np.unique(updated_days)
    logging.debug('Determining available data for each "updated date"')
    count = 0
    if len(extra_jobtypes_scan_archive) != 0 or not nocc:
        for day in updated_days:
            jobs = []
            modified = []
            available = []
            for data in get_data_availability(db, starttime=day, endtime=day+datetime.timedelta(days=1)):
                sta = "%s.%s" % (data.net, data.sta)
                if sta in stations_to_analyse:
                    available.append(sta)
                    if data.flag in ["N", "M"]:
                        modified.append(sta)
            modified = np.unique(modified)
            available = np.unique(available)
            for m in modified:
                for a in available:
                    if m != a or autocorr:
                        pair = ':'.join(sorted([m, a]))
                        if pair not in jobs:
                            if not nocc:
                                all_jobs.append({"day": day, "pair": pair,
                                                 "jobtype": "CC", "flag": "T",
                                                 "lastmod": now})
                            for jobtype in extra_jobtypes_scan_archive:
                                all_jobs.append({"day": day, "pair": pair,
                                             "jobtype": jobtype, "flag": "T",
                                             "lastmod": now})
                            jobs.append(pair)

            if init and len(all_jobs) > 1e5:
                logging.debug('Already 100.000 jobs, inserting/updating')
                massive_insert_job(all_jobs)
                all_jobs = []
                count += 1e5
    else:
        logging.debug("skipping the CC jobs creation & the extrajobtype creation")
    if len(all_jobs) != 0:
        logging.debug('Inserting/Updating %i jobs' % len(all_jobs))
        if init:
            massive_insert_job(all_jobs)
        else:
            for job in all_jobs:
                update_job(db, job['day'], job['pair'],
                           job['jobtype'], job['flag'],
                           commit=False)
    db.commit()
    count += len(all_jobs)

    for sta in get_stations(db, all=True):
        mark_data_availability(db, sta.net, sta.sta, flag='A')

    db.commit()
    logging.info('*** Finished: New Jobs ***')

    return count


if __name__ == "__main__":
    main()