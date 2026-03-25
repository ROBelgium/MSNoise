# TODO add documentation :-)

import os
import traceback
import matplotlib
matplotlib.use("Agg")
import numpy as np
import datetime
from obspy.core import UTCDateTime, read, Stream
from obspy.signal import PPSD

import warnings
warnings.filterwarnings("ignore")

from .api import *

import logbook

def main(loglevel="INFO", njobs_per_worker=9999):
    logger = logbook.Logger("msnoise")
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_psd_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: HDF to RMS ***')
    db = connect()
    logger.debug('Preloading all instrument response')

    params = get_params(db)
    qc_params = get_config_set_details(db, 'qc', 1, format='AttribDict')
    if qc_params:
        params.update(qc_params)
    ppsd_components = params.qc_components
    if not os.path.isdir(os.path.join("PSD", "RMS", params.qc_rms_type)):
        os.makedirs(os.path.join("PSD", "RMS", params.qc_rms_type))

    while is_dtt_next_job(db, jobtype='HDF2RMS'):
        logger.info("Getting the next job")
        jobs = get_dtt_next_job(db, jobtype='HDF2RMS')
        logger.debug("I will process %i jobs" % len(jobs))
        if len(jobs) == 0:
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised code)
            continue
        datelists = {}
        station = None
        for job in jobs:
            net, sta, loc = job.pair.split('.')
            if station is None:
                logger.debug("Processing %s" % (job.pair))
                station = get_station(db, net, sta)
            for chan in station.chans():
                if chan not in datelists:
                    datelists[chan] = []
                if chan[-1] in ppsd_components:
                    datelists[chan].append(datetime.datetime.strptime(job.day, "%Y-%m-%d"))
        for chan in datelists:
            if not len(datelists[chan]):
                continue
            seed_id = "%s.%s.%s.%s" % (net, sta, loc, chan)
            logger.debug("Will open HDFstore: %s" % seed_id)
            store = hdf_open_store(seed_id, mode="r")

            s = datelists[chan][0].strftime("%Y-%m-%d %H:%M:%S")
            e = (datelists[chan][-1] + datetime.timedelta(days=1)).strftime(
                "%Y-%m-%d %H:%M:%S")
            logger.debug("Selecting data between %s and %s" % (s, e))
            data = store.select("PSD", "(index >= '%s') & (index <= '%s')" % (s, e))
            data = store.PSD
            # only need to compute RMS for new/updated PSD data
            data = data.sort_index(axis=1)

            RMS = psd_df_rms(data, freqs=params.qc_rms_frequency_ranges,
                         output=params.qc_rms_type)
            RMSstore = hdf_open_store(seed_id,
                                      location=os.path.join("PSD", "RMS",
                                                            params.qc_rms_type))
            hdf_insert_or_update(RMSstore, "RMS", RMS)
            hdf_close_store(RMSstore)
            del RMS, RMSstore
            del data
        massive_update_job(db, jobs, "D")