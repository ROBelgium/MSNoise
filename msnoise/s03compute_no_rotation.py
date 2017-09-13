""" This code is responsible for the computation of the cross-correlation
functions.

This script will group *jobs* marked "T"odo in the database by day and process
them using the following scheme. As soon as one day is selected, the
corresponding jobs are marked "I"n Progress in the database. This allows
running several instances of this script in parallel.

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |cc_sampling_rate|
* |analysis_duration|
* |overlap|
* |maxlag|
* |corr_duration|
* |windsorizing|
* |resampling_method|
* |remove_response|
* |response_format|
* |response_path|
* |response_prefilt|
* |preprocess_lowpass|
* |preprocess_highpass|
* |keep_all|
* |keep_days|
* |stack_method|
* |pws_timegate|
* |pws_power|
* |whitening|  | *new in 1.5*

Waveform Pre-processing
~~~~~~~~~~~~~~~~~~~~~~~
Pairs are first split and a station list is created. The database is then
queried to get file paths. For each station, all files potentially containing
data for the day are opened. The traces are then merged and splitted, to obtain
the most continuous chunks possible. The different chunks are then demeaned,
tapered and merged again to a 1-day long trace. If a chunk is not aligned
on the sampling grid (that is, start at a integer times the sample spacing in s)
, the chunk is phase-shifted in the frequency domain. This requires tapering and
fft/ifft. If the gap between two chunks is small, compared to a currently
hard-coded value (10 samples), the gap is filled with interpolated values.
Larger gaps will not be filled with interpolated values.

.. warning::
    As from MSNoise 1.5, traces are no longer padded by or merged with 0s.

Each 1-day long trace is then low-passed (at ``preprocess_lowpass`` Hz),
high-passed (at ``preprocess_highpass`` Hz), then if needed,
decimated/downsampled. Decimation/Downsampling are configurable
(``resampling_method``) and users are advised testing Decimate. One advantage of
Downsampling over Decimation is that it is able to downsample the data by any
factor, not only integer factors. Downsampling can be achieved with the new
ObsPy Lanczos resampler, giving results similar to those by scikits.samplerate.

.. note:: Python 3 users will most probably struggle installing
    scikits.samplerate, and therefore will have to use either Decimate or
    Lanczos instead of Resample. This is not a problem because the Lanczos
    resampling gives results similar to those by scikits.samplerate.


If configured, each 1-day long trace is corrected for its instrument response.
Currently, only dataless seed and inventory XML are supported.

As from MSNoise 1.5, the preprocessing routine is separated from the compute_cc
and can be used by plugins with their own parameters. The routine returns a
Stream object containing all the traces for all the stations/components.

Processing
~~~~~~~~~~

Once all traces are preprocessed, station pairs are processed sequentially.
If a component different from *ZZ* is to be computed, the traces are first
rotated. This supposes the user has provided the station coordinates in the
*station* table. The rotation is computed for Radial and Transverse components.

Then, for each ``corr_duration`` window in the signal, and for each filter
configured in the database, the traces are clipped to ``windsorizing`` times
the RMS (or 1-bit converted) and then whitened in the frequency domain
(see :ref:`whiten`) between the frequency bounds. The whitening procedure can be
skipped by setting the ``whitening`` configuration to `None`. The two other
``whitening`` modes are "[A]ll except for auto-correlation" or "Only if
[C]omponents are different". This allows skipping the whitening when, for
example, computing ZZ components for very close by stations (much closer than
the wavelength sampled), leading to spatial autocorrelation issues.

When both traces are ready, the cross-correlation function is computed
(see :ref:`mycorr`). The function returned contains data for time lags
corresponding to ``maxlag`` in the acausal (negative lags) and causal
(positive lags) parts.

Stacking and Saving Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If configured (setting ``keep_all`` to 'Y'), each ``corr_duration`` CCF is
saved to the hard disk. By default, the ``keep_days`` setting is set to True
and so "N = 1 day / corr_duration" CCF are stacked and saved to the hard disk
in the STACKS/001_DAYS folder.

.. note:: Currently, the keep-all data (every CCF) are not used by next steps.

If ``stack_method`` is 'linear', then a simple mean CCF of all windows is saved
as the daily CCF. On the other hand, if ``stack_method`` is 'pws', then
all the Phase Weighted Stack (PWS) is computed and saved as the daily CCF. The
PWS is done in two steps: first the mean coherence between the instataneous
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
    of weak, coherent signals through phase-weighted stacks". Geophysical
    Journal International 130, 2 (1997): 497-505.

Once done, each job is marked "D"one in the database.

To run this script:

.. code-block:: sh

    $ msnoise compute_cc


This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 compute_cc

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.


.. versionadded:: 1.4
    The Instrument Response removal & The Phase Weighted Stack &
    Parallel Processing

.. versionadded:: 1.5
    The Obspy Lanczos resampling method, gives similar results as the
    scikits.samplerate package, thus removing the requirement for it.
    This method is defined by default.

.. versionadded:: 1.5
    The preprocessing routine is separated from the compute_cc and can be called
    by external plugins.

"""
import sys
import time

import matplotlib.mlab as mlab

try:
    from scikits.samplerate import resample
except:
    pass

from .api import *
from .move2obspy import myCorr, myCorr2
from .move2obspy import whiten, whiten2

from .preprocessing import preprocess


class Params():
    pass


def main():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: Compute CC ***')

    # Connection to the DB
    db = connect()

    if len(get_filters(db, all=False)) == 0:
        logging.info("NO FILTERS DEFINED, exiting")
        sys.exit()

    # Get Configuration
    params = Params()
    params.goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    params.goal_duration = float(get_config(db, "analysis_duration"))
    params.overlap = float(get_config(db, "overlap"))
    params.maxlag = float(get_config(db, "maxlag"))
    params.corr_duration = float(get_config(db, "corr_duration"))
    params.min30 = float(get_config(db, "corr_duration")) *\
                   params.goal_sampling_rate
    params.windsorizing = float(get_config(db, "windsorizing"))
    params.whitening = get_config(db, 'whitening')
    params.resampling_method = get_config(db, "resampling_method")
    params.preprocess_lowpass = float(get_config(db, "preprocess_lowpass"))
    params.preprocess_highpass = float(get_config(db, "preprocess_highpass"))
    params.keep_all = get_config(db, 'keep_all', isbool=True)
    params.keep_days = get_config(db, 'keep_days', isbool=True)
    params.components_to_compute = get_components_to_compute(db)

    params.stack_method = get_config(db, 'stack_method')
    params.pws_timegate = float(get_config(db, 'pws_timegate'))
    params.pws_power = float(get_config(db, 'pws_power'))

    logging.info("Will compute %s" % " ".join(params.components_to_compute))

    if get_config(db, 'remove_response', isbool=True):
        logging.debug('Pre-loading all instrument response')
        responses = preload_instrument_responses(db)
    else:
        responses = None

    while is_next_job(db, jobtype='CC'):
        jobs = get_next_job(db, jobtype='CC')

        if len(jobs) == 0:
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            continue

        stations = []
        pairs = []
        refs = []

        for job in jobs:
            refs.append(job.ref)
            pairs.append(job.pair)
            netsta1, netsta2 = job.pair.split(':')
            stations.append(netsta1)
            stations.append(netsta2)
            goal_day = job.day


        stations = np.unique(stations)

        logging.info("New CC Job: %s (%i pairs with %i stations)" %
                     (goal_day, len(pairs), len(stations)))
        jt = time.time()

        comps = []
        for comp in params.components_to_compute:
            if comp[0] in ["R", "T"] or comp[1] in ["R", "T"]:
                comps.append("E")
                comps.append("N")
            else:
                comps.append(comp[0])
                comps.append(comp[1])
        comps = np.unique(comps)
        basetime, stream = preprocess(db, stations, comps, goal_day, params,
                                      responses)
        if not len(stream):
            logging.info("Not enough data for this day !")
            logging.info("Marking job Done and continuing with next !")
            for job in jobs:
                update_job(db, job.day, job.pair, 'CC', 'D')
            continue
        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.goal_sampling_rate
        logging.info("Starting slides")
        start_processing = time.time()
        allcorr = {}
        for tmp in stream.slide(params.corr_duration,
                                params.corr_duration * (1 - params.overlap)):
            logging.info("Processing %s - %s" % (tmp[0].stats.starttime,
                                                 tmp[0].stats.endtime))
            tmp = tmp.copy().sort()

            channels_to_remove = []
            for gap in tmp.get_gaps(min_gap=0):
                if gap[-2] > 0:
                    channels_to_remove.append(
                        ".".join([gap[0], gap[1], gap[2], gap[3]]))

            for chan in np.unique(channels_to_remove):
                logging.debug("%s contains gap(s), removing it" % chan)
                net, sta, loc, chan = chan.split(".")
                for tr in tmp.select(network=net,
                                     station=sta,
                                     location=loc,
                                     channel=chan):
                    tmp.remove(tr)
            if len(tmp) == 0:
                logging.debug("No traces without gaps")
                continue

            base = np.amax([tr.stats.npts for tr in tmp])
            for tr in tmp:
                if tr.stats.npts != base:
                    tmp.remove(tr)
                    logging.debug("One trace is too short, removing it")

            if len(tmp) == 0:
                logging.debug("No traces left in slice")
                continue

            nfft = next_fast_len(tmp[0].stats.npts)
            tmp.detrend("demean")

            for tr in tmp:
                if params.windsorizing == -1:
                    np.sign(tr.data, tr.data)  # inplace
                elif params.windsorizing != 0:
                    rms = tr.data.std() * params.windsorizing
                    np.clip(tr.data, -rms, rms, tr.data)  # inplace
            tmp.taper(0.04)
            Napod = 100

            data = np.asarray([tr.data for tr in tmp])
            names = [tr.id.split(".") for tr in tmp]

            # index net.sta comps for energy later
            channel_index = {}
            psds = []
            for i, name in enumerate(names):
                n1, s1, l1, c1 = name
                netsta = "%s.%s" % (n1, s1)
                if netsta not in channel_index:
                    channel_index[netsta] = {}
                channel_index[netsta][c1[-1]] = i

                Pxx, freqs = mlab.psd(tmp[i].data,
                                      Fs=tmp[i].stats.sampling_rate,
                                      NFFT=nfft)
                psds.append(np.sqrt(Pxx))
            psds = np.asarray(psds)

            for chan in channel_index:
                comps = channel_index[chan].keys()
                if "E" in comps and "N" in comps:
                    iE = channel_index[chan]["E"]
                    iN = channel_index[chan]["N"]
                    # iZ = channel_index[chan]["Z"]
                    mm = psds[[iE,iN]].mean(axis=0)
                    psds[iE] = mm
                    psds[iN] = mm
                    # psds[iZ] = mm

            # define pairwise CCs
            tmptime = tmp[0].stats.starttime.datetime
            thisdate = tmptime.strftime("%Y-%m-%d")
            thistime = tmptime.strftime("%Y-%m-%d %H:%M:%S")
            pair_index = []
            for sta1, sta2 in itertools.combinations(names, 2):
                n1, s1, l1, c1 = sta1
                n2, s2, l2, c2 = sta2
                comp = "%s%s" % (c1[-1], c2[-1])
                if n1 == n2 and s1 == s2:
                    continue
                if comp in params.components_to_compute:
                    pair_index.append(
                        ["%s.%s_%s.%s_%s" % (n1, s1, n2, s2, comp),
                         names.index(sta1), names.index(sta2)])

            for filterdb in get_filters(db, all=False):
                filterid = filterdb.ref
                low = float(filterdb.low)
                high = float(filterdb.high)
                rms_threshold = filterdb.rms_threshold
                fffff = scipy.fftpack.fftfreq(nfft, d=dt)
                freqVec = scipy.fftpack.fftfreq(nfft, d=dt)[:nfft // 2]
                J = np.where((freqVec >= low) & (freqVec <= high))[0]
                low = J[0] - Napod
                if low <= 0:
                    low = 1
                porte1 = J[0]
                porte2 = J[-1]
                high = J[-1] + Napod
                if high > nfft / 2:
                    high = int(nfft // 2)
                ffts = scipy.fftpack.fftn(data, shape=[nfft, ], axes=[1, ])
                whiten2(ffts, nfft, low, high, porte1, porte2, psds)  # inplace

                energy = np.ones(ffts.shape[0])

                # logging.info("Pre-whitened %i traces"%(i+1))
                corr = myCorr2(ffts,
                               np.ceil(params.maxlag / dt),
                               energy,
                               pair_index,
                               plot=False,
                               nfft=nfft)

                for key in corr:
                    ccfid = key + "_%02i" % filterid + "_" + thisdate
                    if ccfid not in allcorr:
                        allcorr[ccfid] = {}
                    allcorr[ccfid][thistime] = corr[key]
                del corr

        for ccfid in allcorr.keys():
            station1, station2, components, filterid, date = ccfid.split('_')

            corrs = np.asarray(list(allcorr[ccfid].values()))
            corr = stack(db, corrs)
            thisdate = goal_day
            thistime = "0_0"
            add_corr(
                db, station1.replace('.', '_'),
                station2.replace('.', '_'), int(filterid),
                thisdate, thistime, params.min30 /
                                    params.goal_sampling_rate,
                components, corr,
                params.goal_sampling_rate, day=True,
                ncorr=corrs.shape[0])

        for job in jobs:
            update_job(db, job.day, job.pair, 'CC', 'D')

        logging.info(
            "Job Finished. It took %.2f seconds (preprocess: %.2f s & process %.2f s)" % (
            (time.time() - jt), start_processing - jt,
            time.time() - start_processing))
        del stream
    logging.info('*** Finished: Compute CC ***')

if __name__ == "__main__":
    main()    

