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

        ix = np.where((d.columns >= fmin) & (d.columns <= fmax))[0]
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
    responses = preload_instrument_responses(db, return_format="inventory")

    params = get_params(db)
    ppsd_components = params.qc_components
    ppsd_length = params.qc_ppsd_length
    ppsd_overlap = params.qc_ppsd_overlap
    ppsd_period_smoothing_width_octaves = params.qc_ppsd_period_smoothing_width_octaves
    ppsd_period_step_octaves = params.qc_ppsd_period_step_octaves
    ppsd_period_limits = params.qc_ppsd_period_limits
    ppsd_db_bins = params.qc_ppsd_db_bins
    for out in ["DISP", "VEL", "ACC"]:
        if not os.path.isdir(os.path.join("PSD", "RMS", out)):
            os.makedirs(os.path.join("PSD", "RMS", out))

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
            store = hdf_open_store(seed_id)
            print(store.PSD.head())

            freqs = [(1.0, 20.0), (4.0, 14.0), (4.0, 40.0), (4.0, 9.0)]

            data = store.PSD
            data = data.sort_index(axis=1)

            displacement_spectrum = df_rms(data, freqs, output="DISP")
            displacement_spectrum.to_csv(os.path.join("PSD","RMS","DISP","%s.csv" % seed_id))

            #     velocity_spectrum = df_rms(data, freqs, output="VEL")
            #     velocity_spectrum.reset_index().to_feather("VEL_%s.feather" % netsta)

            # acc_spectrum = df_rms(data, freqs, output="ACC")
            # acc_spectrum.reset_index().to_feather(
            #     "PSD/RMS/ACC/%s.feather" % seed_id)

            hdf_close_store(store)
            del store, data, displacement_spectrum
        massive_update_job(db, jobs, "D")