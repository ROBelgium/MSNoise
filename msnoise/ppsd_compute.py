"""
This code is responsible for the computation of the Power Spectral Densities.

This script will search for QC *jobs* marked "T"odo in the database
by day and process them. As soon as one day is selected, the
corresponding jobs are marked "I"n Progress in the database. This allows
running several instances of this script in parallel.

The PSD are calculated using the implementation in ObsPy.
Parameters can be defined in the Config to control the windowing and the
smoothness of the PSD results (see below).

The PSD computation produces two SDS-like structures that contain the results,
both under the ``PSD/`` folder:

* ``PSD/NPZ/``: contains the NPZ, or "compressed numpy arrays", with all the
  individual psd spectra
* ``PSD/PNG/``: contains the images representing the daily PPSD (probabilistic
  power spectral densities), which are useful for a routine, or rapid check.

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |qc_components|
* |qc_ppsd_length|
* |qc_ppsd_overlap|
* |qc_ppsd_period_smoothing_width_octaves|
* |qc_ppsd_period_step_octaves|
* |qc_ppsd_period_limits|
* |qc_ppsd_db_bins|


Instrument Response
~~~~~~~~~~~~~~~~~~~

To be able to compute the PSD, the instrument responses need to be provided.
TODO: add link to "how to check if my responses are OK"


To run this script:

.. code-block:: sh

    $ msnoise qc compute_psd


This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 qc compute_psd

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.


.. versionadded:: 2.0
    New in 2.0


"""
# TODO add documentation :-)

import os
import traceback
import matplotlib
matplotlib.use("Agg")
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

    while is_next_job(db, jobtype='PSD'):
        logger.info("Getting the next job")
        jobs = get_next_job(db, jobtype='PSD', limit=njobs_per_worker)
        logger.debug("I will process %i jobs" % len(jobs))
        if len(jobs) == 0:
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised code)
            continue
        for job in jobs:
            net, sta, loc = job.pair.split('.')
            print("Processing %s"% job.pair)
            gd = UTCDateTime(job.day).datetime
            files = get_data_availability(
                db, net=net, sta=sta, loc=loc,
                starttime=(UTCDateTime(job.day) - 1.5 * ppsd_length).datetime,
                endtime=gd)
            if len(files) == 0:
                print("No files found for %s" % job.day)
                continue

            for comp in ppsd_components:
                toprocess = []
                for file in files:
                    if file.chan[-1] != comp:
                        continue
                    tmp = os.path.join(file.path, file.file)
                    toprocess.append(tmp)
                if len(toprocess) == 0:
                    continue
                st = Stream()
                for tmp in np.unique(toprocess):
                    logger.debug("Reading %s" % tmp)
                    try:
                        st += read(tmp, starttime=UTCDateTime(gd) - 1.5 * ppsd_length,
                                   endtime=UTCDateTime(
                                       gd + datetime.timedelta(days=1)) - 0.001)
                    except:
                        logger.debug("Problem loading %s" % tmp)
                if not len(st):
                    continue

                try:
                    st.merge()
                except:
                    logger.info("Failed merging streams:")
                    traceback.print_exc()
                    continue

                st = st.split()
                for tr in st:
                    tr.stats.network = tr.stats.network.upper()
                    tr.stats.station = tr.stats.station.upper()
                    tr.stats.channel = tr.stats.channel.upper()

                tr = st.select(component=comp)[0]
                out = to_sds(tr.stats, gd.year, int(gd.strftime('%j')))
                npzdout = os.path.join("PSD", "NPZ", out)
                logger.debug("ppsd will be output to: %s" % npzdout)
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
                update_job(db, job.day, job.pair, 'PSD', 'D', ref=job.ref)
                if not params.hpc:
                    for job in jobs:
                        update_job(db, job.day, job.pair, 'PSD2HDF', 'T')
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
