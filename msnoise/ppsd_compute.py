# TODO add documentation :-)

import os
import traceback
import numpy as np
import datetime
from obspy.core import UTCDateTime, read, Stream
from obspy.signal import PPSD


from .api import to_sds, get_logger, preload_instrument_responses

import logbook


from msnoise.api import connect, is_next_job, get_next_job, \
    get_data_availability, get_config, update_job, to_sds, get_params


def main(loglevel="INFO", njobs_per_worker=9999):
    logger = logbook.Logger("msnoise")
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_psd_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute PPSD ***')
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

    while is_next_job(db, jobtype='QC'):
        logger.info("Getting the next job")
        jobs = get_next_job(db, jobtype='QC', limit=njobs_per_worker)
        logger.debug("I will process %i jobs" % len(jobs))
        if len(jobs) == 0:
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised code)
            continue
        for job in jobs:
            net, sta = job.pair.split('.')
            gd = UTCDateTime(job.day).datetime
            files = get_data_availability(
                db, net=net, sta=sta,
                starttime=(UTCDateTime(job.day) - 1.5 * 3600).datetime,
                endtime=gd)
            if len(files) == 0:
                continue
            for comp in ppsd_components:
                toprocess = []
                for file in files:
                    if file.comp[-1] != comp:
                        continue
                    tmp = os.path.join(file.path, file.file)
                    toprocess.append(tmp)
                if len(toprocess) == 0:
                    continue
                st = Stream()
                for tmp in np.unique(toprocess):
                    logger.debug("Reading %s" % tmp)
                    try:
                        st += read(tmp, starttime=UTCDateTime(gd) - 1.5 * 3600,
                                   endtime=UTCDateTime(
                                       gd + datetime.timedelta(days=1)) - 0.01)
                    except:
                        logger.debug("Problem loading %s" % tmp)
                if not len(st):
                    continue

                st.merge()
                st = st.split()
                for tr in st:
                    tr.stats.network = tr.stats.network.upper()
                    tr.stats.station = tr.stats.station.upper()
                    tr.stats.channel = tr.stats.channel.upper()

                tr = st.select(component=comp)[0]
                out = to_sds(tr.stats, gd.year, int(gd.strftime('%j')))
                npzdout = os.path.join("PSD", "NPZ", out)
                logger.debug("ppsd will be output to:", npzdout)
                ppsd = PPSD(tr.stats, metadata=responses,
                            ppsd_length=ppsd_length, overlap=ppsd_overlap,
                            period_smoothing_width_octaves=ppsd_period_smoothing_width_octaves,
                            period_step_octaves=ppsd_period_step_octaves,
                            period_limits=ppsd_period_limits,
                            db_bins=ppsd_db_bins
                            )
                # TODO handle when the response for this trace is not in the inv
                ppsd.add(st)
                out = to_sds(tr.stats, gd.year, int(gd.strftime('%j')))

                pngout = os.path.join("PSD", "PNG", out)
                if not os.path.isdir(os.path.split(npzdout)[0]):
                    os.makedirs(os.path.split(npzdout)[0])
                    os.makedirs(os.path.split(pngout)[0])

                ppsd.save_npz(npzdout + ".npz")
                update_job(db, job.day, job.pair, 'QC', 'D', ref=job.ref)
                try:
                    ppsd.plot(pngout + ".png")
                except:
                    logger.debug("Error saving PNG image")
                    traceback.print_exc()

                del ppsd

        logger.debug('Day (job) "D"one')

        # TODO Add massive update / HPC flag
        # for job in jobs:
        #     update_job(db, job.day, job.pair, 'QC', 'D', ref=job.ref)
