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


def dfrms(a):
    return np.sqrt(np.trapz(a.values, a.index))

def rms(s, f):
    return np.sqrt(np.trapz(s, f))

def df_rms(d, freqs, output="VEL"):
    d = d.dropna(axis=1, how='all')
    RMS = {}
    for fmin, fmax in freqs:
        pmin = 1. / fmax
        pmax = 1. / fmin
        ix = np.where((d.columns >= pmin) & (d.columns <= pmax))[0]
        spec = d.iloc[:, ix]
        f = d.columns[ix]

        w2f = (2.0 * np.pi * f)

        # The acceleration power spectrum (dB to Power! = divide by 10 and not 20!)
        amp = 10.0 ** (spec / 10.)
        if output == "ACC":
            RMS["%.1f-%.1f" % (fmin, fmax)] = amp.apply(dfrms, axis=1)
            continue

        # velocity power spectrum (divide by omega**2)
        vamp = amp / w2f ** 2
        if output == "VEL":
            RMS["%.1f-%.1f" % (fmin, fmax)] = vamp.apply(dfrms, axis=1)
            continue

        # displacement power spectrum (divide by omega**2)
        damp = vamp / w2f ** 2

        RMS["%.1f-%.1f" % (fmin, fmax)] = damp.apply(dfrms, axis=1)

    return pd.DataFrame(RMS, index=d.index)  # .tz_localize("UTC")#.dropna()




def main(loglevel="INFO", njobs_per_worker=9999):
    logger = logbook.Logger("msnoise")
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_psd_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: HDF to RMS ***')
    db = connect()
    logger.debug('Preloading all instrument response')

    params = get_params(db)
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
            seed_id = "%s.%s.%s.%s" % (net, sta, loc, chan)
            print("Will open HDFstore: %s" % seed_id)
            store = hdf_open_store(seed_id, mode="r")

            s = datelists[chan][0].strftime("%Y-%m-%d %H:%M:%S")
            e = (datelists[chan][-1] + datetime.timedelta(days=1)).strftime(
                "%Y-%m-%d %H:%M:%S")
            print("Selecting data between %s and %s" % (s, e))
            data = store.select("PSD", "(index >= '%s') & (index <= '%s')" % (s, e))
            data = store.PSD
            # only need to compute RMS for new/updated PSD data
            data = data.sort_index(axis=1)

            RMS = df_rms(data, freqs=params.qc_rms_frequency_ranges,
                         output=params.qc_rms_type)
            RMSstore = hdf_open_store(seed_id,
                                      location=os.path.join("PSD", "RMS",
                                                            params.qc_rms_type))
            hdf_insert_or_update(RMSstore, "RMS", RMS)
            hdf_close_store(RMSstore)
            del RMS, RMSstore
            del data
        massive_update_job(db, jobs, "D")