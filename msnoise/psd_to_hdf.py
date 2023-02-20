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
    logger.info('*** Starting: PSD to HDF ***')
    db = connect()
    logger.debug('Preloading all instrument response')
    responses = preload_instrument_responses(db, return_format="inventory")

    params = get_params(db)
    ppsd_components = params.qc_components

    while is_dtt_next_job(db, jobtype='PSD2HDF'):
        logger.info("Getting the next job")
        jobs = get_dtt_next_job(db, jobtype='PSD2HDF')
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
                print("Processing %s" % (job.pair))
                station = get_station(db, net, sta)
            for chan in station.chans():
                if chan not in datelists:
                    datelists[chan] = []
                if chan[-1] in ppsd_components:
                    datelists[chan].append(datetime.datetime.strptime(job.day, "%Y-%m-%d"))
        for chan in datelists:
            if not len(datelists[chan]):
                continue
            # TODO ADD (optional) INTERPOLATION OF THE PERIOD BINS, TO AVOID CRASHES WHEN CHANGING THE SPS OF THE STATION -> HDF WILL NOT LIKE IT !
            ppsd = psd_read_results(net, sta, loc, chan, datelists[chan], use_cache=False)
            if ppsd is None or not len(ppsd.times_processed):
                continue
            new = psd_ppsd_to_dataframe(ppsd)
            store = hdf_open_store("%s.%s.%s.%s" % (net, sta, loc, chan))
            hdf_insert_or_update(store, "PSD", new)
            hdf_close_store(store)
            del ppsd, new, store
        massive_update_job(db, jobs, "D")
        if not params.hpc:
            for job in jobs:
                update_job(db, job.day, job.pair, 'HDF2RMS', 'T')