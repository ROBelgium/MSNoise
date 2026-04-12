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

The following parameters (modifiable via ``msnoise admin``) are used for
this step:

* |cc.components_to_compute|
* |cc.components_to_compute_single_station|
* |cc.cc_sampling_rate|
* |cc.cc_normalisation|
* |cc.cc_type|
* |cc.cc_type_single_station_AC|
* |cc.cc_type_single_station_SC|
* |cc.cc_taper_fraction|
* |cc.clip_after_whiten|
* |cc.overlap|
* |cc.maxlag|
* |cc.corr_duration|
* |cc.winsorizing|
* |cc.whitening|
* |cc.whitening_type|
* |cc.keep_all|
* |cc.keep_days|
* |cc.stack_method|
* |cc.pws_timegate|
* |cc.pws_power|
* |preprocess.resampling_method|
* |preprocess.remove_response|
* |global.response_path|
* |preprocess.response_prefilt|
* |preprocess.preprocess_lowpass|
* |preprocess.preprocess_highpass|
* |preprocess.preprocess_max_gap|
* |preprocess.preprocess_taper_length|
* |global.hpc|

.. automodule:: msnoise.core.preprocessing

Computing the Cross-Correlations
--------------------------------

Processing using ``msnoise cc compute_cc``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: We still need to describe the workflow in plain text, but the
    following graph should help you understand how the code is structured

.. image:: ../.static/compute_cc.png
    :align: center

.. image:: ../.static/MyCorr2.png
    :align: center

Processing using ``msnoise cc compute_cc_rot``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once all traces are preprocessed, station pairs are processed sequentially.
If a component different from *ZZ* is to be computed, the traces are first
rotated. This supposes the user has provided the station coordinates in the
*station* table. The rotation is computed for Radial and Transverse components.

Then, for each ``corr_duration`` window in the signal, and for each filter
configured in the database, the traces are clipped to ``winsorizing`` times
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

    $ msnoise cc compute_cc

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc compute_cc

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

    The Instrument Response removal & The Phase Weighted Stack &
    Parallel Processing

    The Obspy Lanczos resampling method, gives similar results as the
    scikits.samplerate package, thus removing the requirement for it.
    This method is defined by default.

    The preprocessing routine is separated from the compute_cc and can be called
    by external plugins.

    The compute_cc has been completely rewritten to be much faster, taking
    advantage from 2D FFT computation and in-place array modifications.
    The standard compute_cc does process CC, AC and SC in the same code. Only
    if users need to compute R and/or T components, they will have to use the
    slower previous code, now called ``compute_cc_rot``.

"""
#TODO docstring
import itertools
import logging
import os
import time

import numpy as np

from .core.db import connect, get_logger
from .core.config import get_config_set_details
from .core.workflow import (get_filter_steps_for_cc_step, get_next_lineage_batch, get_t_axis, is_next_job_for_step, massive_update_job, propagate_downstream)
from .core.signal import stack, winsorizing, get_preprocessed_stream
from .core.io import save_daily_ccf, xr_save_ccf_all
from .core.compute import myCorr2
from .core.compute import whiten2
from .core.compute import pcc_xcorr



import scipy.signal
import scipy.fft as sf
from scipy.fft import next_fast_len
from obspy.signal.filter import bandpass


def main(loglevel="INFO", chunk_size=0):
    global logger
    logger = get_logger('msnoise.cc', loglevel, with_pid=True)
    logger.info('*** Starting: Compute CC ***')

    # Connection to the DB
    db = connect()

    # if len(get_filters(db, all=False)) == 0:
    #     logger.info("NO FILTERS DEFINED, exiting")
    #     sys.exit()

    # Get Configuration
    # filters = get_filters(db, all=False)
    # logger.info("Will compute [%s] for different stations" % " ".join(params.cc.components_to_compute))
    # logger.info("Will compute [%s] for single stations" % " ".join(params.cc.components_to_compute_single_station))

    # if "R" in ''.join(params.cc.components_to_compute) or "T" in ''.join(params.cc.components_to_compute):
    #     logger.info("You seem to have configured R and/or T components, thus rotations ARE needed. You should therefore use the 'msnoise compute_cc_rot' instead.")
    #     return()
    #
    # if params.cc.whitening not in ["A", "N"]:
    #     logger.info("The 'whitening' parameter is set to '%s', which is not supported by this process. Set it to 'A' or 'N', or use the 'msnoise compute_cc_rot' instead." % params.cc.whitening)
    #     return ()
    #
    # if params.preprocess.remove_response:
    #     logger.debug('Pre-loading all instrument response')
    #     responses = preload_instrument_responses(db, return_format="inventory")
    # else:
    #     responses = None
    logger.debug("Checking if there are jobs to do")
    if chunk_size > 0:
        logger.info(f"CC chunk_size={chunk_size}: each worker claims up to {chunk_size} pairs per day")

    while is_next_job_for_step(db, step_category="cc"):
        batch = get_next_lineage_batch(db, step_category="cc", group_by="day_lineage",
                                       chunk_size=chunk_size, loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        params = batch["params"]
        lineage_names_upstream = batch["lineage_names_upstream"]
        lineage_names = batch["lineage_names"]  # full: e.g. [preprocess_1, cc_1]
        step = batch["step"]

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
        t_axis = get_t_axis(params)

        # Filters:
        filter_steps = get_filter_steps_for_cc_step(db, step.step_id)
        filters = []
        for filter in filter_steps:
            filters.append([filter.step_name, get_config_set_details(db, filter.category, filter.set_number, format='AttribDict')])
        logger.info("New CC Job: %s (%i pairs with %i stations)" %
                     (goal_day, len(pairs), len(stations)))
        jt = time.time()

        # Read per-station preprocessed files — one file per NET.STA.LOC
        # under _output/<goal_day>/<NET.STA.LOC>.mseed (v2 layout).
        _preprocess_step = lineage_names_upstream[-1] if lineage_names_upstream else ""
        _preprocess_out  = os.path.join(params.global_.output_folder,
                                        *lineage_names_upstream[:-1])
        stream = get_preprocessed_stream(
            _preprocess_out, _preprocess_step, goal_day, stations
        )
        # TODO PREPROCESS IF THE "PREPROCESS_ON_THE_FLY" config?
        # comps = np.unique(comps)
        # with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        #     stream = executor.submit(preprocess, stations, comps, goal_day, params, responses, loglevel).result()
        # logger.info("Received preprocessed traces")


        # stream = preprocess(db, stations, comps, goal_day, params, responses)
        if not len(stream):
            logger.debug("Not enough data for this day !")
            logger.debug("Marking jobs 'Fail' and continuing with next !")
            massive_update_job(db, jobs, "F")
            continue
        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.cc.cc_sampling_rate
        logger.debug("Starting slides")
        start_processing = time.time()
        allcorr = {}
        for tmp in stream.slide(params.cc.corr_duration,
                                params.cc.corr_duration * (1 - params.cc.overlap)):
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
            if base <= (params.cc.maxlag*params.cc.cc_sampling_rate*2+1):
                logger.debug("All traces shorter are too short to export"
                              " +-maxlag")
                continue

            for tr in tmp:
                if tr.stats.npts != base:
                    tmp.remove(tr)
                    logger.debug("%s trace is too short, removing it" % tr.id)

            if len(tmp) == 0:
                logger.debug("No traces left in slice")
                continue

            nfft = next_fast_len(tmp[0].stats.npts)
            tmp.detrend("demean")

            # Spectral end-taper width: ~0.5% of nfft, minimum 50 bins.
            # Tapers the edges of the FFT output before whitening to
            # suppress spectral leakage at DC and Nyquist.
            napod = max(50, nfft // 200)

            data = np.asarray([tr.data for tr in tmp])
            names = [tr.id.split(".") for tr in tmp]

            if not params.cc.clip_after_whiten:
                # logger.debug("Winsorizing (clipping) data before whiten")
                data = winsorizing(data, params) #inplace

            # Time-domain cosine taper fraction (config: cc_taper_fraction, default 4%).
            # Applied symmetrically to both ends of each trace before CC.
            wlen = int(params.cc.cc_taper_fraction * data.shape[1])
            taper_sides = scipy.signal.windows.hann(2 * wlen + 1)
            taper = np.hstack(
                (taper_sides[:wlen], np.ones(data.shape[1] - 2 * wlen),
                 taper_sides[len(taper_sides) - wlen:]))
            for i in range(data.shape[0]):
                data[i] *= taper
            # index net.sta comps for energy later
            channel_index = {}
            if params.cc.whitening != "N" and params.cc.whitening_type == "PSD":
                psds = []
                for i, name in enumerate(names):
                    n1, s1, l1, c1 = name
                    netsta = "%s.%s.%s" % (n1, s1, l1)
                    if netsta not in channel_index:
                        channel_index[netsta] = {}
                    channel_index[netsta][c1[-1]] = i

                    freqs, pxx = scipy.signal.welch(tmp[i].data,
                                               fs=tmp[i].stats.sampling_rate,
                                               nperseg=nfft,
                                               detrend='constant')
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
            # Pairwise CC index: for each (sta1, sta2) combination we build
            # a list of (i, j, component_pair) tuples so the inner FFT/IFFT
            # loop can be fully vectorised across all pairs in one batch
            # (see myCorr2 in core/compute.py — replaces the old per-pair loop).

            tmptime = tmp[0].stats.starttime.datetime
            thisdate = tmptime.strftime("%Y-%m-%d")
            tmptime = tmp[0].stats.endtime.datetime
            thistime = tmptime.strftime("%Y-%m-%d %H:%M:%S")

            for filter_name, filter in filters:
                # logger.debug("Processing filter %s" % filter_name)
                # print(filter)

                # Standard operator for CC
                cc_index = []
                if filter.CC:
                    if len(params.cc.components_to_compute):
                        for sta1, sta2 in itertools.combinations(names, 2):
                            n1, s1, l1, c1 = sta1
                            n2, s2, l2, c2 = sta2
                            if n1 == n2 and s1 == s2 and l1 == l2:
                                continue
                            pair = "%s.%s.%s:%s.%s.%s" % (n1, s1, l1, n2, s2, l2)
                            if pair not in pairs:
                                continue
                            comp = "%s%s" % (c1[-1], c2[-1])
                            if comp in params.cc.components_to_compute:
                                cc_index.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                    names.index(sta1), names.index(sta2)])

                # Different iterator func for single station AC or SC:
                single_station_pair_index_sc = []
                single_station_pair_index_ac = []

                if len(params.cc.components_to_compute_single_station):
                    for sta1, sta2 in itertools.combinations_with_replacement(names, 2):
                        n1, s1, l1, c1 = sta1
                        n2, s2, l2, c2 = sta2
                        if n1 != n2 or s1 != s2:
                            continue
                        pair = "%s.%s.%s:%s.%s.%s" % (n1, s1, l1, n2, s2, l2)
                        if pair not in pairs:
                            continue
                        comp = "%s%s" % (c1[-1], c2[-1])
                        if comp in params.cc.components_to_compute_single_station:
                            if filter.AC and c1[-1] == c2[-1]:
                                single_station_pair_index_ac.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                    names.index(sta1), names.index(sta2)])
                            elif filter.SC:
                            # If the components are different, we can just
                            # process them using the default CC code (should warn)
                                single_station_pair_index_sc.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                    names.index(sta1), names.index(sta2)])
                        if comp[::-1] in params.cc.components_to_compute_single_station:
                            if filter.SC and c1[-1] != c2[-1]:
                                # If the components are different, we can just
                                # process them using the default CC code (should warn)
                                single_station_pair_index_sc.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp[::-1]),
                                    names.index(sta2), names.index(sta1)])

                # logger.debug("CC index: %s" % cc_index)
                # logger.debug("single station AC: ", single_station_pair_index_ac)
                # logger.debug("single station SC: ", single_station_pair_index_sc)

                filterlow = float(filter.freqmin)
                filterhigh = float(filter.freqmax)
                # logger.debug("Filter: %s - %s Hz" % (filterlow, filterhigh))
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
                if params.cc.whitening == "N":
                    # if the data doesn't need to be whitened, we can simply
                    # band-pass filter the traces now:
                    for i, _ in enumerate(_data):
                        _data[i] = bandpass(_, freqmin=filterlow,
                                            freqmax=filterhigh,
                                            df=params.cc.cc_sampling_rate,
                                            corners=8)

                # First let's compute the AC and SC
                if len(single_station_pair_index_ac):
                    tmp = _data.copy()
                    if params.cc.whitening == "A":
                        # if the data isn't already filtered, we still need to
                        # do it for the AutoCorrelation:
                        for i, _ in enumerate(tmp):
                            tmp[i] = bandpass(_, freqmin=filterlow,
                                              freqmax=filterhigh,
                                              df=params.cc.cc_sampling_rate,
                                              corners=8)
                    if params.cc.clip_after_whiten:
                        # logger.debug("Winsorizing (clipping) data after bandpass (AC)")
                        tmp = winsorizing(tmp, params, input="timeseries")

                    if params.cc.cc_type_single_station_AC == "CC":
                        ffts = sf.fftn(tmp, [nfft, ], axes=[1, ])
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2,
                            axis=1)))
                        # Computing standard CC
                        corr = myCorr2(ffts,
                                       np.ceil(params.cc.maxlag / dt),
                                       energy,
                                       single_station_pair_index_ac,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc.cc_normalisation)
                        del energy, ffts

                    elif params.cc.cc_type_single_station_AC == "PCC":
                        corr = pcc_xcorr(tmp,
                                         np.ceil(params.cc.maxlag / dt),
                                         None,
                                         single_station_pair_index_ac,
                                         normalized=params.cc.cc_normalisation)
                    else:
                        logging.error("cc_type_single_station_AC = %s not implemented, "
                              "exiting")
                        exit(1)

                    for key in corr:
                        ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                        logger.debug("CCF ID: %s" % ccfid)
                        if ccfid not in allcorr:
                            allcorr[ccfid] = {}
                        allcorr[ccfid][thistime] = corr[key]
                    del corr

                if len(cc_index):
                    if params.cc.cc_type == "CC":
                        ffts = sf.fftn(_data, [nfft, ], axes=[1, ])
                        if params.cc.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                        if params.cc.clip_after_whiten:
                            # logger.debug(
                            #     "Winsorizing (clipping) data after whiten")
                            ffts = winsorizing(ffts, params, input="fft", nfft=nfft)

                        # energy = np.sqrt(np.sum(np.abs(ffts)**2, axis=1)/nfft)
                        energy = np.real(np.sqrt( np.mean(sf.ifft(ffts, n=nfft, axis=1) ** 2, axis=1)))

                        # logger.info("Pre-whitened %i traces"%(i+1))
                        # Computing standard CC
                        corr = myCorr2(ffts,
                                       np.ceil(params.cc.maxlag / dt),
                                       energy,
                                       cc_index,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc.cc_normalisation)

                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            # logger.debug("CCF ID - CC: %s" % ccfid)
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, energy, ffts

                    elif params.cc.cc_type == "PCC":
                        # Phase Cross-Correlation (v=2, FFT-accelerated).
                        # Operates on time-domain data; amplitude is discarded
                        # per-sample → insensitive to transients without
                        # explicit temporal normalisation.
                        corr = pcc_xcorr(_data,
                                         np.ceil(params.cc.maxlag / dt),
                                         None,
                                         cc_index,
                                         normalized=params.cc.cc_normalisation)
                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr

                    else:
                        logging.error("cc_type = %s not implemented, "
                              "exiting")
                        exit(1)

                if len(single_station_pair_index_sc):
                    if params.cc.cc_type_single_station_SC == "CC":
                        # logger.debug("Compute SC using %s" % params.cc.cc_type)
                        ffts = sf.fftn(_data, [nfft, ], axes=[1, ])
                        if params.cc.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                        if params.cc.clip_after_whiten:
                            ffts = winsorizing(ffts, params, input="fft",
                                               nfft=nfft)
                        # energy = np.sqrt(np.sum(np.abs(ffts)**2, axis=1)/nfft)
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2,
                            axis=1)))

                        # logger.info("Pre-whitened %i traces"%(i+1))
                        # Computing standard CC
                        corr = myCorr2(ffts,
                                       np.ceil(params.cc.maxlag / dt),
                                       energy,
                                       single_station_pair_index_sc,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc.cc_normalisation)

                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            logger.debug("CCF ID - SC: %s" % ccfid)
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, energy, ffts

                    elif params.cc.cc_type_single_station_SC == "PCC":
                        # Phase Cross-Correlation (v=2, FFT-accelerated).
                        corr = pcc_xcorr(_data,
                                         np.ceil(params.cc.maxlag / dt),
                                         None,
                                         single_station_pair_index_sc,
                                         normalized=params.cc.cc_normalisation)
                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            logger.debug("CCF ID - SC PCC: %s" % ccfid)
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr

                    else:
                        logging.error("cc_type_single_station_SC = %s not implemented, "
                              "exiting")
                        exit(1)
            del psds

        if params.cc.keep_all:
            for ccfid, windows in allcorr.items():
                station1, station2, components, filterid, date = ccfid.split('+')
                window_times = list(windows.keys())
                corrs = np.asarray(list(windows.values()))
                xr_save_ccf_all(
                    root=params.global_.output_folder,
                    lineage=lineage_names,
                    step_name=filterid,  # filterid from ccfid = filter step name e.g. 'filter_1'
                    station1=station1,
                    station2=station2,
                    components=components,
                    date=date,
                    window_times=window_times,
                    taxis=t_axis,
                    corrs=corrs,
                )

        if params.cc.keep_days:
            for ccfid in allcorr.keys():
                logging.debug("Exporting %s" % ccfid)
                station1, station2, components, filter_name, date = \
                    ccfid.split('+')

                corrs = np.asarray(list(allcorr[ccfid].values()))
                if not len(corrs):
                    logger.debug("No data to stack.")
                    continue
                corr = stack(corrs, params.cc.stack_method, params.cc.pws_timegate,
                             params.cc.pws_power, params.cc.cc_sampling_rate)
                if not len(corr):
                    logger.debug("No data to save.")
                    continue
                thisdate = goal_day
                thistime = "0_0"
                save_daily_ccf(
                    root=params.global_.output_folder,
                    lineage=lineage_names,
                    step_name=filter_name,  # filter_name from ccfid = filter step name e.g. 'filter_1'
                    station1=station1,
                    station2=station2,
                    components=components,
                    date=thisdate,
                    corr=corr,
                    taxis=t_axis,
                )

        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)

        logger.info("Job Finished. It took %.2f seconds (preprocess: %.2f s & "
                     "process %.2f s)" % ((time.time() - jt),
                                          start_processing - jt,
                                          time.time() - start_processing))
        del stream, allcorr
    logger.info('*** Finished: Compute CC ***')
