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
* |decimation_factor|
* |remove_response|
* |response_format|
* |response_path|
* |response_prefilt|
* |preprocess_lowpass|
* |preprocess_highpass|
* |keep_all|
* |keep_days|
* |stack_method| | *new in 1.4*
* |pws_timegate| | *new in 1.4*
* |pws_power| | *new in 1.4*


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
Larger gaps will not be filled with interpolated values and remaining chunks
will be tapered and then merged with 0 values in the gaps.

If shorter than 1-day, the trace final is padded with zeros. If longer, it is
cut to match the start/end of the day.

If configured, each 1-day long trace is corrected for its instrument response.
Currently, only dataless seed and inventory XML are supported.

.. note:: Removing the instrument response is a computationally very expensive
   task and *not* useful for dv/v iff your instruments didn't change during the
   analysed period. It is also not needed for tomography iff all instruments are
   the same, or at least have an identical phase response in the frequency band
   of interest.

Each 1-day long trace is then low-passed (at ``preprocess_lowpass`` Hz),
high-passed (at ``preprocess_highpass`` Hz), then if needed,
decimated/downsampled. Decimation/Downsampling are configurable
(``resampling_method``) and users are advised testing both. One advantage of
Downsampling over Decimation is that it is able to downsample the data by any
factor, not only integer factors.

.. note:: Python 3 users will most probably struggle installing
    scikits.samplerate, and therefore will have to use either Decimate or
    Lanczos instead of Resample. This is not a problem because the Lanczos
    resampling give very similar results as the scikits.samplerate package.

Processing
~~~~~~~~~~

Once all traces are preprocessed, station pairs are processed sequentially.
If a component different from *ZZ* is to be computed, the traces are first
rotated. This supposes the user has provided the station coordinates in the
*station* table. The rotation is computed for Radial and Transverse components:

.. code-block:: python

    R = N * np.cos(Az * np.pi / 180.) + E * np.sin(Az * np.pi / 180.)
    T = N * np.sin(Az * np.pi / 180.) - E * np.cos(Az * np.pi / 180.)

Then, for each ``corr_duration`` window in the signal, and for each filter
configured in the database, the traces are clipped to ``windsorizing`` times
the RMS (or 1-bit converted) and then whitened (see :ref:`whiten`) between the
frequency bounds.
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

If ``stack_method`` is 'linear', then a simple mean CFF of all windows is saved
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

.. seealso:: Schimmel, M. and Paulssen H., "Noise reduction and detection
    of weak, coherent signals through phase-weighted stacks". Geophysical Journal
    International 130, 2 (1997): 497-505.

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

"""
import sys
import time

try:
    from scikits.samplerate import resample
except:
    pass

from .api import *
from .move2obspy import myCorr
from .move2obspy import whiten

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
    params.min30 = float(get_config(db, "corr_duration")) * params.goal_sampling_rate
    params.windsorizing = float(get_config(db, "windsorizing"))
    params.resampling_method = get_config(db, "resampling_method")
    params.decimation_factor = int(get_config(db, "decimation_factor"))
    params.preprocess_lowpass = float(get_config(db, "preprocess_lowpass"))
    params.preprocess_highpass = float(get_config(db, "preprocess_highpass"))
    params.keep_all = get_config(db, 'keep_all', isbool=True)
    params.keep_days = get_config(db, 'keep_days', isbool=True)
    params.components_to_compute = get_components_to_compute(db)

    params.stack_method = get_config(db, 'stack_method')
    params.pws_timegate = float(get_config(db, 'pws_timegate'))
    params.pws_power = float(get_config(db, 'pws_power'))

    logging.info("Will compute %s" % " ".join(params.components_to_compute))

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
        del jobs

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
        basetime, stream = preprocess(db, stations, comps, goal_day, params)

        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.goal_sampling_rate

        # ITERATING OVER PAIRS #####
        for pair in pairs:
            orig_pair = pair

            logging.info('Processing pair: %s' % pair.replace(':', ' vs '))
            tt = time.time()
            station1, station2 = pair.split(':')
            pair = (np.where(stations == station1)
                    [0][0], np.where(stations == station2)[0][0])

            s1 = get_station(db, station1.split('.')[0], station1.split('.')[1])
            s2 = get_station(db, station2.split('.')[0], station2.split('.')[1])

            if s1.X and params.components_to_compute != ["ZZ",]:
                X0, Y0, c0 = (s1.X, s1.Y, s1.coordinates)
                X1, Y1, c1 = (s2.X, s2.Y, s1.coordinates)

                if c0 == c1:
                    coordinates = c0
                else:
                    coordinates = 'MIX'

                cplAz = azimuth(coordinates, X0, Y0, X1, Y1)
                logging.info("Azimuth=%.1f"%cplAz)
            else:
                # logging.debug('No Coordinates found! Skipping azimuth calculation!')
                cplAz = 0.

            for components in params.components_to_compute:
                t1 = stream.select(station=s1.sta)
                t2 = stream.select(station=s2.sta)
                if (components == "ZZ") \
                        or ("E" in components)\
                        or ("N" in components)\
                        or ("1" in components)\
                        or ("2" in components):
                    t1 = t1.select(component=components[0])
                    t2 = t2.select(component=components[1])
                else:
                    t1 = t1.rotate("NE->RT", cplAz).select(component=components[0])
                    t2 = t2.rotate("NE->RT", cplAz).select(component=components[1])

                if not len(t1):
                    logging.info("No Data for %s.%s..%s" % (
                        s1.net, s1.sta, components[0]))
                    continue
                if not len(t2):
                    logging.info("No Data for %s.%s..%s" % (
                        s2.net, s2.sta, components[1]))
                    continue

                current = t1+t2
                # print(current)

                allcorr = {}
                for tmp in current.slide(params.corr_duration, params.corr_duration*(1-params.overlap)):
                    gaps = []
                    for gap in tmp.get_gaps(min_gap=0):
                        if gap[-2] > 0:
                            gaps.append(gap)

                    if len(gaps) > 0:
                        logging.debug("Sliding Windows %s contains gaps, skipping..." % (tmp[0].stats.starttime))
                        continue
                    if tmp[0].stats.npts < 2*(params.maxlag * params.goal_sampling_rate) + 1:
                        continue
                    if len(tmp) < 2:
                        continue
                    tmp = tmp.copy()
                    tmp.detrend("demean")

                    for tr in tmp:
                        if params.windsorizing == -1:
                            tr.data = np.sign(tr.data)
                        elif params.windsorizing != 0:
                            rms = tr.data.std()
                            indexes = np.where(np.abs(tr.data) > (params.windsorizing * rms))[0]
                            # clipping at windsorizing*rms
                            tr.data[indexes] = (tr.data[indexes] / np.abs(
                                tr.data[indexes])) * params.windsorizing * rms
                    tmp.taper(0.04)
                    tmp1 = tmp.select(station=s1.sta, component=components[0])
                    if len(tmp1) == 0:
                        continue
                    else:
                        tmp1 = tmp1[0]
                    tmp2 = tmp.select(station=s2.sta, component=components[1])
                    if len(tmp2) == 0:
                        continue
                    else:
                        tmp2 = tmp2[0]
                    nfft = next_fast_len(tmp1.stats.npts)
                    autocorr = False
                    if (s1.net == s2.net) and (s1.sta == s2.sta) and (
                        components[0] == components[1]):
                        autocorr = True

                    for filterdb in get_filters(db, all=False):
                        filterid = filterdb.ref
                        low = float(filterdb.low)
                        high = float(filterdb.high)
                        rms_threshold = filterdb.rms_threshold

                        trames2hWb = np.zeros((2, int(nfft)), dtype=np.complex)
                        skip = False
                        for i, station in enumerate(pair):
                            if tmp[i].data.std() > rms_threshold:
                                if autocorr:
                                    # logging.debug("Autocorr %s"%components)
                                    tmp[i].filter("bandpass", freqmin=low, freqmax=high, zerophase=True)
                                    trames2hWb[i] = scipy.fftpack.fft(tmp[i].data, nfft)
                                else:
                                    # logging.debug("Whitening %s" % components)
                                    trames2hWb[i] = whiten(tmp[i].data, nfft,
                                                           dt, low, high,
                                                           plot=False)
                            else:
                                skip = True
                                logging.debug('Slice RMS is smaller (%e) than rms_threshold (%e)!'
                                              % (tmp[i].data.std(), rms_threshold))
                        if not skip:
                            corr = myCorr(trames2hWb, np.ceil(params.maxlag / dt), plot=False, nfft=nfft)
                            if not np.all(np.isfinite(corr)):
                                logging.debug("corr object contains NaNs, skipping")
                                continue
                            if len(corr) < 2* (params.maxlag * params.goal_sampling_rate) + 1:
                                logging.debug(
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
                        corr = stack(db, corrs)

                        thisdate = time.strftime(
                                    "%Y-%m-%d", time.gmtime(basetime))
                        thistime = time.strftime(
                                    "%H_%M", time.gmtime(basetime))
                        add_corr(
                                db, station1.replace('.', '_'),
                                station2.replace('.', '_'), int(filterid),
                                thisdate, thistime,  params.min30 /
                                params.goal_sampling_rate,
                                components, corr,
                                params.goal_sampling_rate, day=True,
                                ncorr=corrs.shape[0])
                        del corrs, corr, thisdate, thistime
                del current, allcorr, t1, t2
            logging.debug("Updating Job")
            update_job(db, goal_day, orig_pair, 'CC', 'D')

            logging.info("Finished processing this pair. It took %.2f seconds" % (time.time() - tt))
        clean_scipy_cache()
        logging.info("Job Finished. It took %.2f seconds" % (time.time() - jt))
    logging.info('*** Finished: Compute CC ***')

if __name__ == "__main__":
    main()    

