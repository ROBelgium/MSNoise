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
* |preprocess_max_gap|  | *new in 1.6*
* |keep_all|
* |keep_days|
* |stack_method|
* |pws_timegate|
* |pws_power|
* |whitening|  | *new in 1.5*
* |hpc| | *new in 1.6*

Waveform Pre-processing
~~~~~~~~~~~~~~~~~~~~~~~
Pairs are first split and a station list is created. The database is then
queried to get file paths. For each station, all files potentially containing
data for the day are opened. The traces are then merged and splitted, to obtain
the most continuous chunks possible. The different chunks are then demeaned,
tapered and merged again to a 1-day long trace. If a chunk is not aligned
on the sampling grid (that is, start at a integer times the sample spacing in s)
, the chunk is phase-shifted in the frequency domain. This requires tapering and
fft/ifft. If the gap between two chunks is small, compared to a configurable
value (``preprocess_max_gap``), the gap is filled with interpolated values.
Larger gaps will not be filled with interpolated values.

.. warning::
    As from MSNoise 1.5, traces are no longer padded by or merged with 0s.

Each 1-day long trace is then high-passed (at ``preprocess_highpass`` Hz),
then if needed low-passed (at ``preprocess_lowpass`` Hz) and
decimated/downsampled. Decimation/Downsampling are configurable
(``resampling_method``) and users are advised testing Decimate. One advantage of
Downsampling over Decimation is that it is able to downsample the data by any
factor, not only integer factors. Downsampling can be achieved with the new
ObsPy Lanczos resampler.

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

Once done, each job is marked "D"one in the database and, unless ``hpc`` is
``Y``, STACK jobs are inserted/updated in the database.

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

.. versionadded:: 1.6
    The ``preprocess_max_gap`` parameter.
    The ``hpc`` parameter that can prevent the automatic creation of STACK jobs.
    


"""
import sys
import time

from .api import *
from .move2obspy import myCorr
from .move2obspy import whiten

from .preprocessing import preprocess


import logbook


def main(loglevel="INFO"):
    logger = logbook.Logger("msnoise")
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_cc_child', loglevel,
                        with_pid=True)

    logger.info('*** Starting: Compute CC ***')

    # Connection to the DB
    db = connect()

    if len(get_filters(db, all=False)) == 0:
        logger.info("NO FILTERS DEFINED, exiting")
        sys.exit()

    # Get Configuration
    params = get_params(db)
    filters = get_filters(db, all=False)

    logger.info("Will compute %s" % " ".join(params.components_to_compute))
    
    if "R" not in ''.join(params.components_to_compute) and "T" not in ''.join(params.components_to_compute):
        logger.info("You seem to have configured no R nor T components, thus no rotation are needed. You should therefore use the 'msnoise compute_cc' instead, which is much faster")

    if params.remove_response:
        logger.debug('Pre-loading all instrument response')
        responses = preload_instrument_responses(db)
    else:
        responses = None

    while is_next_job(db, jobtype='CC'):
        jobs = get_next_job(db, jobtype='CC')

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

        logger.info("New CC Job: %s (%i pairs with %i stations)" %
                     (goal_day, len(pairs), len(stations)))
        jt = time.time()

        comps = []
        for comp in params.all_components:
            for c in comp:
                if c in ["R", "T"]:
                    comps.append("E")
                    comps.append("N")
                else:
                    comps.append(c)

        comps = np.unique(comps)
        stream = preprocess(db, stations, comps, goal_day, params, responses)

        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.goal_sampling_rate
        start_processing = time.time()
        # ITERATING OVER PAIRS #####
        logger.info("Will do %i pairs" % len(pairs))
        for job in jobs:
            pair = job.pair
            orig_pair = pair

            logger.info('Processing pair: %s' % pair.replace(':', ' vs '))
            tt = time.time()
            station1, station2 = pair.split(':')
            pair = (np.where(stations == station1)
                    [0][0], np.where(stations == station2)[0][0])

            s1 = get_station(db, station1.split('.')[0], station1.split('.')[1])
            s2 = get_station(db, station2.split('.')[0], station2.split('.')[1])

            if s1.X and params.components_to_compute != ["ZZ", ]:
                x0, y0, c0 = (s1.X, s1.Y, s1.coordinates)
                x1, y1, c1 = (s2.X, s2.Y, s1.coordinates)

                if c0 == c1:
                    coordinates = c0
                else:
                    coordinates = 'MIX'

                cpl_az = azimuth(coordinates, x0, y0, x1, y1)
                logger.info("Azimuth=%.1f"%cpl_az)
            else:
                cpl_az = 0.

            for components in params.components_to_compute:
                logger.debug("Processing {!s}".format(components))
                t1 = stream.select(network=s1.net, station=s1.sta)
                t2 = stream.select(network=s2.net, station=s2.sta)
                if (components == "ZZ") \
                        or ("E" in components)\
                        or ("N" in components)\
                        or ("1" in components)\
                        or ("2" in components):
                    t1 = t1.select(component=components[0])
                    t2 = t2.select(component=components[1])
                else:
                    logger.debug('Rotating streams, making sure they are'
                                  ' aligned')

                    if components[0] == "Z":
                        t1 = t1.select(component=components[0])
                    else:
                        t1_novert = t1.copy()
                        for tr in t1_novert.select(component="Z"):
                            t1_novert.remove(tr)
                        if len(t1_novert):
                            # Make these streams contain the same gaps
                            t1_novert = make_same_length(t1_novert)
                            t1 = t1_novert.rotate("NE->RT", cpl_az).\
                                select(component=components[0])
                        else:
                            t1 = t1_novert

                    if components[1] == "Z":
                        t2 = t2.select(component=components[1])
                    else:
                        t2_novert = t2.copy()
                        for tr in t2_novert.select(component="Z"):
                            t2_novert.remove(tr)
                        if len(t2_novert):
                            # Make these streams contain the same gaps
                            t2_novert = make_same_length(t2_novert)
                            t2 = t2_novert.rotate("NE->RT", cpl_az).\
                                select(component=components[1])
                        else:
                            t2 = t2_novert

                if not len(t1):
                    logger.info("No Data for %s.%s..%s" % (
                        s1.net, s1.sta, components[0]))
                    continue
                if not len(t2):
                    logger.info("No Data for %s.%s..%s" % (
                        s2.net, s2.sta, components[1]))
                    continue

                current = t1+t2

                allcorr = {}
                for tmp in current.slide(params.corr_duration,
                                         params.corr_duration *
                                         (1-params.overlap)):
                    gaps = []
                    for gap in tmp.get_gaps(min_gap=0):
                        if gap[-2] > 0:
                            gaps.append(gap)

                    if len(gaps) > 0:
                        logger.debug("Sliding Windows %s contains gaps,"
                                      " skipping..." % tmp[0].stats.starttime)
                        continue
                    if tmp[0].stats.npts < 2*(params.maxlag *
                                              params.goal_sampling_rate) + 1:
                        continue
                    if len(tmp) < 2:
                        continue
                    tmp = tmp.copy()
                    tmp.detrend("demean")

                    for tr in tmp:
                        if params.windsorizing == -1:
                            np.sign(tr.data, tr.data)  # inplace
                        elif params.windsorizing != 0:
                            rms = tr.data.std() * params.windsorizing
                            np.clip(tr.data, -rms, rms, tr.data)  # inplace

                    # TODO should not hardcode 4 percent!
                    tmp.taper(0.04)
                    tmp1 = tmp.select(network=s1.net, station=s1.sta,
                                      component=components[0])
                    if len(tmp1) == 0:
                        continue

                    tmp2 = tmp.select(network=s2.net, station=s2.sta,
                                      component=components[1])
                    if len(tmp2) == 0:
                        continue

                    tmp1 = tmp1[0]
                    tmp2 = tmp2[0]
                    nfft = next_fast_len(tmp1.stats.npts)

                    whitening = True
                    if params.whitening == "A":
                        if (s1.net == s2.net) and (s1.sta == s2.sta) and (
                           components[0] == components[1]):
                            whitening = False
                    elif params.whitening == "C":
                        if components[0] == components[1]:
                            whitening = False
                    elif params.whitening == "N":
                        whitening = False

                    for filterdb in filters:
                        filterid = filterdb.ref
                        low = float(filterdb.low)
                        high = float(filterdb.high)
                        rms_threshold = filterdb.rms_threshold

                        trames2hWb = np.zeros((2, int(nfft)), dtype=np.complex)
                        skip = False
                        for i, station in enumerate(pair):
                            if tmp[i].data.std() > rms_threshold:
                                if whitening:
                                    #logger.debug("Whitening %s" % components)
                                    trames2hWb[i] = whiten(tmp[i].data, nfft,
                                                           dt, low, high,
                                                           plot=False)
                                else:
                                    #logger.debug("Autocorr %s"%components)
                                    X = tmp[i].copy()
                                    X.filter("bandpass", freqmin=low,
                                             freqmax=high, zerophase=True)
                                    trames2hWb[i] = scipy.fftpack.fft(X.data, nfft)
                            else:
                                skip = True
                                logger.debug('Slice RMS is smaller (%e) than rms_threshold (%e)!'
                                              % (tmp[i].data.std(), rms_threshold))
                        if not skip:
                            corr = myCorr(trames2hWb, np.ceil(params.maxlag / dt), plot=False, nfft=nfft)
                            if not np.all(np.isfinite(corr)):
                                logger.debug("corr object contains NaNs, skipping")
                                continue
                            if len(corr) < 2 * (params.maxlag * params.goal_sampling_rate) + 1:
                                logger.debug(
                                    "corr object is too small, skipping")
                                continue
                            tmptime = tmp[0].stats.starttime.datetime
                            thisdate = tmptime.strftime("%Y-%m-%d")
                            thistime = tmptime.strftime("%Y-%m-%d %H:%M:%S")
                            if params.keep_all or params.keep_days:
                                ccfid = "%s_%s_%s_%s_%s" % (station1, station2,
                                                         filterid, components,
                                                         thisdate)
                                if ccfid not in allcorr:
                                    allcorr[ccfid] = {}
                                allcorr[ccfid][thistime] = corr

                            del corr, thistime, trames2hWb, tmptime
                        del low, high, rms_threshold
                    del tmp, tmp1, tmp2

                if params.keep_all:
                    for ccfid in allcorr.keys():
                        export_allcorr(db, ccfid, allcorr[ccfid])

                if params.keep_days:
                    for ccfid in allcorr.keys():
                        station1, station2, filterid, components, date = ccfid.split('_')

                        corrs = np.asarray(list(allcorr[ccfid].values()))
                        corr = stack(corrs, params.stack_method,
                                     params.pws_timegate, params.pws_power,
                                     params.goal_sampling_rate)

                        thisdate = goal_day
                        thistime = "0_0"
                        add_corr(
                                db, station1.replace('.', '_'),
                                station2.replace('.', '_'), int(filterid),
                                thisdate, thistime,  params.min30 /
                                params.goal_sampling_rate,
                                components, corr,
                                params.goal_sampling_rate, day=True,
                                ncorr=corrs.shape[0], params=params)
                        del corrs, corr, thisdate, thistime
                del current, allcorr, t1, t2
            logger.debug("Updating Job")
            # Would be better after the massive update at the end of the day job
            # but here it allows to only insert the stack job if the CC was 
            # successful.
            if not params.hpc:
                update_job(db, goal_day, orig_pair, 'STACK', 'T')

            logger.info("Finished processing this pair. It took %.2f seconds" % (time.time() - tt))
        massive_update_job(db, jobs, "D")
        clean_scipy_cache()

        logger.info("Job Finished. It took %.2f seconds (preprocess: %.2f s & process %.2f s)" % ((time.time() - jt), start_processing-jt, time.time()-start_processing))
        del stream
    logger.info('*** Finished: Compute CC ***')
