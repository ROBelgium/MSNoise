""" This code is responsible for the computation of the cross-correlation
functions.

This script will group *jobs* marked "T"odo in the database by day and process
them using the following scheme. As soon as one day is selected, the
corresponding jobs are marked "I"n Progress in the database. This allows
running several instances of this script in parallel.

As of MSNoise 1.6, the ``compute`` step has been completely rewritten:

* The ``compute_cc`` step has been completely rewritten to make use of 2D arrays
  holding the data, processing them "in place" for the different steps (FFT,
  whitening, etc). This results in much more efficient computation. The process
  slides on time windows and computes the correlations using indexes in a 2D
  array, therefore avoiding an exponential number of identical operations on
  data windows.

* This new code is the default ``compute_cc``, and it doesn't allow computing
  rotated components. For users needing ``R`` or ``T`` components, there are two
  options: either use the old code, now named ``compute_cc_rot``, or compute the
  full (6 components actually are enough) tensor using the new code, and rotate
  the components afterwards. From initial tests, this latter solution is a lot
  faster than the first, thanks to the new processing in 2D.

* It is now possible to do the Cross-Correlation (classic "CC"), the Auto-
  Correlation ("AC") or the Cross-Components within the same station ("SC").
  To achieve this, we removed the `ZZ`, `ZT`, etc parameters from the
  configuration and replaced it with ``components_to_compute`` which takes a
  list: e.g. `ZZ,ZE,ZN,EZ,EE,EN,NZ,NE,NN` for the full non-rotated tensor
  between stations. Adding components to the new
  ``components_to_compute_single_station`` will allow computing the
  cross-components (SC) or auto-correlation (AC) of each station.

* The cross-correlation is done on sliding windows on the available data. For
  each window, if one trace contains a gap, it is eliminated from the
  computation. This corrects previous errors linked with gaps synchronised in
  time that lead to perfect sinc autocorrelation functions. The windows should
  have a duration of at least "2 times the `maxlag`+1" to be computable.


Configuration Parameters
------------------------

The following parameters (modifiable via ```msnoise admin```) are used for
this step:

* |components_to_compute|
* |components_to_compute_single_station|  | *new in 1.6*
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
* |preprocess_taper_length|  | *new in 1.6*
* |keep_all|
* |keep_days|
* |stack_method|
* |pws_timegate|
* |pws_power|
* |whitening|  | *new in 1.5*
* |whitening_type|  | *new in 1.6*
* |hpc| | *new in 1.6*

.. automodule:: msnoise.preprocessing


Computing the Cross-Correlations
--------------------------------

Processing using ``msnoise compute_cc``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: We still need to describe the workflow in plain text, but the
    following graph should help you understand how the code is structured


.. image:: ../.static/compute_cc.png
    :align: center


.. image:: ../.static/MyCorr2.png
    :align: center




Processing using ``msnoise compute_cc_rot``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Saving Results (stacking the daily correlations)
------------------------------------------------

If configured (setting ``keep_all`` to 'Y'), each ``corr_duration`` CCF is
saved to the hard disk in the ``output_folder``. By default, the ``keep_days``
setting is set to True and so "N = 1 day / corr_duration" CCF are stacked to
produce a daily cross-correlation function, which is saved to the hard disk in
the ``STACKS/001_DAYS`` folder.

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


Usage
-----

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
    The compute_cc has been completely rewritten to be much faster, taking
    advantage from 2D FFT computation and in-place array modifications.
    The standard compute_cc does process CC, AC and SC in the same code. Only
    if users need to compute R and/or T components, they will have to use the
    slower previous code, now called ``compute_cc_rot``.

"""
#TODO docstring
import sys
import time

import matplotlib.mlab as mlab

from .api import *
from .move2obspy import myCorr2
from .move2obspy import whiten2
from .move2obspy import pcc_xcorr

from .preprocessing import preprocess

from scipy.stats import scoreatpercentile

import scipy.signal
import scipy.fftpack as sf
from scipy.fftpack import next_fast_len
from obspy.signal.filter import bandpass


import logbook


def winsorizing(data, params, input="timeseries", nfft=0):
    input1D = False
    if len(data.shape) == 1:
        data = data.reshape(-1, data.shape[0])
        input1D = True
    if input == "fft":
        # data = np.real(sf.ifft(data, n=nfft, axis=1))
        data = sf.ifftn(data, [nfft, ], axes=[1, ]).astype(np.float64)

    for i in range(data.shape[0]):
        if params.windsorizing == -1:
            np.sign(data[i], data[i])  # inplace
        elif params.windsorizing != 0:

            # imin, imax = scoreatpercentile(data[i], [0.1, 99.9])
            # not_outliers = np.where((data[i] >= imin) &
            #                         (data[i] <= imax))[0]
            # rms = data[i][not_outliers].std() * params.windsorizing
            rms = data[i].std() * params.windsorizing
            np.clip(data[i], -rms, rms, data[i])  # inplace
    if input == "fft":
        data = sf.fftn(data, [nfft, ], axes=[1, ])
    if input1D:
        data = data[0]
    return data


def main(loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_cc_norot_child', loglevel,
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
    logger.info("Will compute [%s] for different stations" % " ".join(params.components_to_compute))
    logger.info("Will compute [%s] for single stations" % " ".join(params.components_to_compute_single_station))

    if "R" in ''.join(params.components_to_compute) or "T" in ''.join(params.components_to_compute):
        logger.info("You seem to have configured R and/or T components, thus rotations ARE needed. You should therefore use the 'msnoise compute_cc_rot' instead.")
        return()
    
    if params.whitening not in ["A", "N"]:
        logger.info("The 'whitening' parameter is set to '%s', which is not supported by this process. Set it to 'A' or 'N', or use the 'msnoise compute_cc_rot' instead." % params.whitening)
        return ()

    if params.remove_response:
        logger.debug('Pre-loading all instrument response')
        responses = preload_instrument_responses(db, return_format="inventory")
    else:
        responses = None
    logger.info("Checking if there are jobs to do")
    while is_next_job(db, jobtype='CC'):
        logger.info("Getting the next job")
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
            comps.append(comp[0])
            comps.append(comp[1])

        comps = np.unique(comps)
        stream = preprocess(db, stations, comps, goal_day, params, responses)
        if not len(stream):
            logger.debug("Not enough data for this day !")
            logger.debug("Marking job Done and continuing with next !")
            for job in jobs:
                update_job(db, job.day, job.pair, 'CC', 'D', ref=job.ref)
            continue
        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.goal_sampling_rate
        logger.debug("Starting slides")
        start_processing = time.time()
        allcorr = {}
        for tmp in stream.slide(params.corr_duration,
                                params.corr_duration * (1 - params.overlap)):
            logger.debug("Processing %s - %s" % (tmp[0].stats.starttime,
                                                 tmp[0].stats.endtime))
            tmp = tmp.copy().sort()

            channels_to_remove = []
            for gap in tmp.get_gaps(min_gap=0):
                if gap[-2] > 0:
                    channels_to_remove.append(
                        ".".join([gap[0], gap[1], gap[2], gap[3]]))

            for chan in np.unique(channels_to_remove):
                logger.debug("%s contains gap(s), removing it" % chan)
                net, sta, loc, chan = chan.split(".")
                for tr in tmp.select(network=net,
                                     station=sta,
                                     location=loc,
                                     channel=chan):
                    tmp.remove(tr)
            if len(tmp) == 0:
                logger.debug("No traces without gaps")
                continue

            base = np.amax([tr.stats.npts for tr in tmp])
            if base <= (params.maxlag*params.goal_sampling_rate*2+1):
                logger.debug("All traces shorter are too short to export"
                              " +-maxlag")
                continue

            for tr in tmp:
                if tr.stats.npts != base:
                    tmp.remove(tr)
                    logger.debug("One trace is too short, removing it")

            if len(tmp) == 0:
                logger.debug("No traces left in slice")
                continue

            nfft = next_fast_len(tmp[0].stats.npts)
            tmp.detrend("demean")

            # TODO should not hardcode 100 taper points in spectrum
            napod = 100

            data = np.asarray([tr.data for tr in tmp])
            names = [tr.id.split(".") for tr in tmp]

            if not params.clip_after_whiten:
                logger.debug("Winsorizing (clipping) data before whiten")
                data = winsorizing(data, params) #inplace

            # TODO should not hardcode 4 percent!
            wlen = int(0.04 * data.shape[1])
            taper_sides = scipy.signal.hann(2 * wlen + 1)
            taper = np.hstack(
                (taper_sides[:wlen], np.ones(data.shape[1] - 2 * wlen),
                 taper_sides[len(taper_sides) - wlen:]))
            for i in range(data.shape[0]):
                data[i] *= taper
            # index net.sta comps for energy later
            channel_index = {}
            if params.whitening != "N" and params.whitening_type == "PSD":
                psds = []
                for i, name in enumerate(names):
                    n1, s1, l1, c1 = name
                    netsta = "%s.%s.%s" % (n1, s1, l1)
                    if netsta not in channel_index:
                        channel_index[netsta] = {}
                    channel_index[netsta][c1[-1]] = i
    
                    pxx, freqs = mlab.psd(tmp[i].data,
                                          Fs=tmp[i].stats.sampling_rate,
                                          NFFT=nfft,
                                          detrend='mean')
                    psds.append(np.sqrt(pxx))
                psds = np.asarray(psds)
            else:
                psds = np.zeros(1)

            for chan in channel_index:
                comps = channel_index[chan].keys()
                if "E" in comps and "N" in comps:
                    i_e = channel_index[chan]["E"]
                    i_n = channel_index[chan]["N"]
                    # iZ = channel_index[chan]["Z"]
                    mm = psds[[i_e,i_n]].mean(axis=0)
                    psds[i_e] = mm
                    psds[i_n] = mm
                    # psds[iZ] = mm

            # define pairwise CCs
            tmptime = tmp[0].stats.starttime.datetime
            thisdate = tmptime.strftime("%Y-%m-%d")
            thistime = tmptime.strftime("%Y-%m-%d %H:%M:%S")

            # Standard operator for CC
            cc_index = []
            if len(params.components_to_compute):
                for sta1, sta2 in itertools.combinations(names, 2):
                    n1, s1, l1, c1 = sta1
                    n2, s2, l2, c2 = sta2
                    pair = "%s.%s.%s:%s.%s.%s" % (n1, s1, l1, n2, s2, l2)
                    if pair not in pairs:
                        continue
                    comp = "%s%s" % (c1[-1], c2[-1])
                    if comp in params.components_to_compute:
                        cc_index.append(
                            ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                             names.index(sta1), names.index(sta2)])

            # Different iterator func for single station AC or SC:
            single_station_pair_index_sc = []
            single_station_pair_index_ac = []

            if len(params.components_to_compute_single_station):
                for sta1, sta2 in itertools.combinations_with_replacement(names, 2):
                    n1, s1, l1, c1 = sta1
                    n2, s2, l2, c2 = sta2
                    if n1 != n2 or s1 != s2:
                        continue
                    pair = "%s.%s.%s:%s.%s.%s" % (n1, s1, l1, n2, s2, l2)
                    if pair not in pairs:
                        continue
                    comp = "%s%s" % (c1[-1], c2[-1])
                    if comp in params.components_to_compute_single_station:
                        if c1[-1] == c2[-1]:
                            single_station_pair_index_ac.append(
                                ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                 names.index(sta1), names.index(sta2)])
                        else:
                        # If the components are different, we can just
                        # process them using the default CC code (should warn)
                            single_station_pair_index_sc.append(
                                ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                 names.index(sta1), names.index(sta2)])
                    if comp[::-1] in params.components_to_compute_single_station:
                        if c1[-1] != c2[-1]:
                            # If the components are different, we can just
                            # process them using the default CC code (should warn)
                            single_station_pair_index_sc.append(
                                ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp[::-1]),
                                 names.index(sta2), names.index(sta1)])

            for filterdb in filters:
                filterid = filterdb.ref
                filterlow = float(filterdb.low)
                filterhigh = float(filterdb.high)

                freq_vec = sf.fftfreq(nfft, d=dt)[:nfft // 2]
                freq_sel = np.where((freq_vec >= filterlow) & (freq_vec <= filterhigh))[0]
                low = freq_sel[0] - napod
                if low <= 0:
                    low = 0
                p1 = freq_sel[0]
                p2 = freq_sel[-1]
                high = freq_sel[-1] + napod
                if high > nfft / 2:
                    high = int(nfft // 2)

                # Make a copy of the original data to prevent modifying it
                _data = data.copy()
                if params.whitening == "N":
                    # if the data doesn't need to be whitened, we can simply
                    # band-pass filter the traces now:
                    for i, _ in enumerate(_data):
                        _data[i] = bandpass(_, freqmin=filterlow,
                                            freqmax=filterhigh,
                                            df=params.goal_sampling_rate,
                                            corners=8)

                # First let's compute the AC and SC
                if len(single_station_pair_index_ac):
                    tmp = _data.copy()
                    if params.whitening == "A":
                        # if the data isn't already filtered, we still need to
                        # do it for the AutoCorrelation:
                        for i, _ in enumerate(tmp):
                            tmp[i] = bandpass(_, freqmin=filterlow,
                                              freqmax=filterhigh,
                                              df=params.goal_sampling_rate,
                                              corners=8)
                            if params.clip_after_whiten:
                                logger.debug("Winsorizing (clipping) data after bandpass (AC)")
                                tmp[i] = winsorizing(tmp[i], params, input="timeseries")


                    if params.cc_type_single_station_AC == "CC":
                        ffts = sf.fftn(tmp, [nfft, ], axes=[1, ])
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2,
                            axis=1)))

                        # Computing standard CC
                        corr = myCorr2(ffts,
                                       np.ceil(params.maxlag / dt),
                                       energy,
                                       single_station_pair_index_ac,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc_normalisation)

                    elif params.cc_type_single_station_AC == "PCC":
                        corr = pcc_xcorr(tmp, np.ceil(params.maxlag / dt),
                                         None, single_station_pair_index_ac)
                    else:
                        print("cc_type_single_station_AC = %s not implemented, "
                              "exiting")
                        exit(1)

                    for key in corr:
                        ccfid = key + "_%02i" % filterid + "_" + thisdate
                        if ccfid not in allcorr:
                            allcorr[ccfid] = {}
                        allcorr[ccfid][thistime] = corr[key]
                    del corr, energy

                if len(cc_index):
                    if params.cc_type == "CC":
                        ffts = sf.fftn(_data, [nfft, ], axes=[1, ])
                        if params.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.whitening_type)  # inplace
                        if params.clip_after_whiten:
                            logger.debug(
                                "Winsorizing (clipping) data after whiten")
                            ffts = winsorizing(ffts, params, input="fft", nfft=nfft)

                        # energy = np.sqrt(np.sum(np.abs(ffts)**2, axis=1)/nfft)
                        energy = np.real(np.sqrt( np.mean(sf.ifft(ffts, n=nfft, axis=1) ** 2, axis=1)))
        
                        # logger.info("Pre-whitened %i traces"%(i+1))
                        # Computing standard CC
                        corr = myCorr2(ffts,
                                       np.ceil(params.maxlag / dt),
                                       energy,
                                       cc_index,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc_normalisation)
                        
                        for key in corr:
                            ccfid = key + "_%02i" % filterid + "_" + thisdate
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, energy, ffts
                    else:
                        print("cc_type = %s not implemented, "
                              "exiting")
                        exit(1)

                if len(single_station_pair_index_sc):
                    if params.cc_type_single_station_SC == "CC":
                        # logger.debug("Compute SC using %s" % params.cc_type)
                        ffts = sf.fftn(_data, [nfft, ], axes=[1, ])
                        if params.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.whitening_type)  # inplace
                        if params.clip_after_whiten:
                            ffts = winsorizing(ffts, params, input="fft",
                                               nfft=nfft)
                        # energy = np.sqrt(np.sum(np.abs(ffts)**2, axis=1)/nfft)
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2,
                            axis=1)))

                        # logger.info("Pre-whitened %i traces"%(i+1))
                        # Computing standard CC
                        corr = myCorr2(ffts,
                                       np.ceil(params.maxlag / dt),
                                       energy,
                                       single_station_pair_index_sc,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc_normalisation)

                        for key in corr:
                            ccfid = key + "_%02i" % filterid + "_" + thisdate
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, energy, ffts
                    else:
                        print("cc_type_single_station_SC = %s not implemented, "
                              "exiting")
                        exit(1)
            del psds

        if params.keep_all:
            for ccfid in allcorr.keys():
                export_allcorr2(db, ccfid, allcorr[ccfid])

        if params.keep_days:
            for ccfid in allcorr.keys():
                print("Exporting %s" % ccfid)
                station1, station2, components, filterid, date = \
                    ccfid.split('_')

                corrs = np.asarray(list(allcorr[ccfid].values()))
                if not len(corrs):
                    logger.debug("No data to stack.")
                    continue
                corr = stack(corrs, params.stack_method, params.pws_timegate,
                             params.pws_power, params.goal_sampling_rate)
                if not len(corr):
                    logger.debug("No data to save.")
                    continue
                thisdate = goal_day
                thistime = "0_0"
                add_corr(
                    db, station1, station2, int(filterid),
                    thisdate, thistime, params.min30 /
                                        params.goal_sampling_rate,
                    components, corr,
                    params.goal_sampling_rate, day=True,
                    ncorr=corrs.shape[0],
                    params=params)

        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")
        if not params.hpc:
            for job in jobs:
                update_job(db, job.day, job.pair, 'STACK', 'T')

        logger.info("Job Finished. It took %.2f seconds (preprocess: %.2f s & "
                     "process %.2f s)" % ((time.time() - jt),
                                          start_processing - jt,
                                          time.time() - start_processing))
        del stream, allcorr
    logger.info('*** Finished: Compute CC ***')
