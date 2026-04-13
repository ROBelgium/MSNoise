r"""Computation of inter-station and single-station cross-correlation functions.

This module implements the ``msnoise cc compute`` worker.  It groups *jobs*
marked **T**odo in the database by day, marks them **I**n-progress, processes
all sliding windows and filters for that day, then marks them **D**one.
Multiple instances can run in parallel because each worker atomically claims
its own batch of jobs.


Overview
--------

The worker operates on 2-D arrays of shape ``(n_stations, N)`` — one row per
pre-processed trace, one column per sample — and processes all station pairs
simultaneously for each time window.  Two independent correlation algorithms
are available:

* **CC** — classic Generalised Normalised Cross-Correlation (GNCC), computed
  in the frequency domain via FFT (see `Cross-Correlation (CC)`_ below).
* **PCC** — Phase Cross-Correlation v2 (PCC2), a transient-robust
  alternative that operates in the time domain before taking any FFT (see
  `Phase Cross-Correlation (PCC2)`_ below).

The algorithm is selected per correlation mode through three independent
configuration parameters: ``cc_type`` (inter-station CC),
``cc_type_single_station_AC`` (auto-correlation) and
``cc_type_single_station_SC`` (same-station cross-component).


Cross-Correlation (CC)
----------------------

The classic ambient-noise cross-correlation (Shapiro & Campillo 2004;
Bensen et al. 2007) is computed using the tiled-batch implementation
:func:`~msnoise.core.compute.myCorr2`.

**Processing chain for each time window and filter band:**

1. **Temporal normalisation** — optional Winsorising (clipping to
   ``winsorizing`` × RMS) applied either *before* whitening (default) or
   *after* (``clip_after_whiten = Y``).

2. **Time-domain taper** — a cosine (Hann) taper of width
   ``cc_taper_fraction × N`` is applied symmetrically to both ends to
   suppress spectral leakage.

3. **FFT** — each trace is zero-padded to ``nfft = next_fast_len(N)`` and
   transformed to the frequency domain:

   .. math::

       X_i[k] = \mathcal{F}\{x_i[n]\}, \quad k = 0, \ldots, N_\text{fft}-1

4. **Spectral whitening** — controlled by ``whitening`` and ``whitening_type``
   (see :func:`~msnoise.core.compute.whiten2`).  Three modes are available:

   - ``whitening = N`` — no whitening; a bandpass filter is applied in the
     time domain instead to restrict each window to the current filter band.
   - ``whitening = A`` — whiten all traces (including auto-correlations).
   - ``whitening = C`` — whiten only when the two components differ.

   Three spectral shapes are supported via ``whitening_type``:

   - *Brutal* (default) — sets ``|X[k]| = 1`` inside the passband with a
     cosine taper at the edges; equivalent to retaining only the spectral
     phase.
   - *HANN* — Hann-weighted one-bit normalisation inside the passband.
   - *PSD* — divides by the smoothed power spectral density, then clips
     outlier bins at the 5th–95th percentile.

5. **Cross-spectrum and IFFT** — for each requested pair :math:`(i, j)`:

   .. math::

       C_{ij}[k] = X_i^*[k] \cdot X_j[k]

   .. math::

       c_{ij}[\tau] = \mathcal{F}^{-1}\{C_{ij}[k]\} \,/\, N

   The output is folded so that negative lags :math:`[-\tau_\text{max}, 0)`
   precede positive lags :math:`[0, \tau_\text{max}]`, yielding a CCF of
   length :math:`2\,\tau_\text{max} + 1` samples.

6. **Normalisation** — optional post-processing controlled by
   ``cc_normalisation``:

   - ``POW`` — divide by the product of the per-station RMS energies
     :math:`(e_i \cdot e_j)`.
   - ``MAX`` — divide by the maximum value of the CCF.
   - ``ABSMAX`` — divide by the absolute maximum.

**References**

Bensen, G. D., Ritzwoller, M. H., Barmin, M. P., Levshin, A. L., Lin, F.,
Moschetti, M. P., Shapiro, N. M., & Yang, Y. (2007).
Processing seismic ambient noise data to obtain reliable broad-band surface
wave dispersion measurements.
*Geophysical Journal International*, 169(3), 1239–1260.
https://doi.org/10.1111/j.1365-246X.2007.03374.x

Shapiro, N. M., & Campillo, M. (2004).
Emergence of broadband Rayleigh waves from correlations of the ambient
seismic noise.
*Geophysical Research Letters*, 31(7), L07614.
https://doi.org/10.1029/2004GL019491


Phase Cross-Correlation (PCC2)
------------------------------

Phase Cross-Correlation (Schimmel 1999; Schimmel et al. 2011) is an
alternative to the classic GNCC that discards amplitude information entirely
at every sample, making it intrinsically robust against transient noise
(earthquakes, instrumental glitches) *without* requiring explicit temporal
normalisation such as one-bit or clipping.  MSNoise implements PCC version 2
(PCC2), the FFT-accelerated formulation introduced by Ventosa et al. (2019,
2023) as part of the FastPCC package.  The implementation is self-contained
in :func:`~msnoise.core.compute.pcc_xcorr` and does **not** depend on the
``phasecorr`` or ``fastpcc`` external packages.

**Algorithm**

For each trace :math:`x_i[n]` of length :math:`N`:

1. **Band-pass filter** — the trace is band-pass filtered (8-pole Butterworth,
   zero-phase) to the current filter band :math:`[f_\text{low}, f_\text{high}]`
   *before* any further processing.  This is essential: because PCC discards
   spectral amplitudes, the filter band is the *only* mechanism that makes each
   ``filter_N`` configuration produce a distinct cross-correlation.  (For CC,
   the same role is played by the frequency-domain passband of ``whiten2``.)

2. **Optional spectral whitening** — if ``whitening ≠ N``, the bandpass output
   is FFT-whitened within the filter band and transformed back to the time
   domain.  This distributes energy evenly across all frequencies inside the
   passband *before* amplitude normalisation, so that the PCC phase signal
   reflects broadband phase coherence rather than being dominated by the most
   energetic frequency in the band.

3. **Analytic signal and amplitude normalisation** — the complex analytic
   signal :math:`x_i^{(a)}[n]` is computed via the Hilbert transform (FFT-based,
   :func:`scipy.signal.hilbert`).  Each sample is then divided by its own
   amplitude to produce the *phase signal* :math:`\varphi_i`:

   .. math::

       \varphi_i[n] = \frac{x_i^{(a)}[n]}{|x_i^{(a)}[n]| + \varepsilon}

   where :math:`\varepsilon = 10^{-6} \max_n |x_i^{(a)}[n]|` is a numerical
   stability floor (matching FastPCC's ``AmpNormf`` convention).
   By construction, :math:`|\varphi_i[n]| \leq 1` for all :math:`n`, and
   :math:`|\varphi_i[n]| \approx 1` wherever the signal is not near zero.
   This per-sample normalisation is what makes PCC insensitive to amplitude
   transients: a spike 1000 × larger than the ambient noise is reduced to
   exactly the same weight as any other sample.

4. **Zero-padding** — each phase signal is zero-padded to
   :math:`N_z = \text{next\_fast\_len}(N + \tau_\text{max})` to compute a
   linear (non-circular) cross-correlation and avoid wrap-around artefacts.

5. **FFT cross-spectrum and IFFT** — all phase signals are pre-transformed
   once, then for each pair :math:`(i, j)`:

   .. math::

       \text{PCC2}_{ij}[\tau] = \mathcal{F}^{-1}\!\left\{
           \Phi_i^*[k] \cdot \Phi_j[k]
       \right\} \Big/ N

   where :math:`\Phi_i[k] = \mathcal{F}\{\varphi_i[n]\}` (zero-padded to
   :math:`N_z`).  Division by :math:`N` (not :math:`N_z`) gives a peak
   amplitude of approximately 1 for identical signals, consistent with
   PCC2 as a cross-coherence measure.

6. **Lag unwrapping** — the IFFT output of length :math:`N_z` stores positive
   lags :math:`\tau = 0, \ldots, \tau_\text{max}` at indices
   :math:`0, \ldots, \tau_\text{max}`, and negative lags
   :math:`\tau = -\tau_\text{max}, \ldots, -1` at indices
   :math:`N_z - \tau_\text{max}, \ldots, N_z - 1` (Ventosa's convention).
   The two slices are concatenated to give the final CCF of length
   :math:`2\,\tau_\text{max} + 1`:

   .. math::

       \text{PCC2}_{ij}[-\tau_\text{max}{:}\tau_\text{max}]
       = \bigl[\text{full}[N_z - \tau_\text{max}{:}N_z],\;
                \text{full}[0{:}\tau_\text{max}+1]\bigr]

   .. note::

       Because :math:`N_z > N`, the extracted CCF is **not** symmetric around
       lag 0 even for a signal correlated with itself.  This is a mathematical
       property of zero-padded linear correlation and does not indicate an
       implementation error.  The self-correlation *does* have its global
       maximum at lag 0.

7. **Normalisation** — same ``cc_normalisation`` options as for CC
   (``MAX``, ``ABSMAX``, or none).  The ``POW`` option is silently ignored
   for PCC2 because amplitude information has already been discarded in step 3.

**Computational cost**

PCC2 requires one additional Hilbert transform (one FFT + IFFT per trace)
relative to CC, but the dominant cost — :math:`O(P \cdot N_z \log N_z)` for
:math:`P` pairs — is identical.  In practice the overhead is negligible for
large networks.

**References**

Schimmel, M. (1999).
Phase cross-correlations: design, comparisons, and applications.
*Bulletin of the Seismological Society of America*, 89(5), 1366–1378.
https://doi.org/10.1785/BSSA0890051366

Schimmel, M., Stutzmann, E., & Gallart, J. (2011).
Using instantaneous phase coherence for signal extraction from ambient noise
data at a local to a global scale.
*Geophysical Journal International*, 184(1), 494–506.
https://doi.org/10.1111/j.1365-246X.2010.04861.x

Ventosa, S., Schimmel, M., & Stutzmann, E. (2019).
Towards the processing of large data volumes with phase cross-correlation.
*Seismological Research Letters*, 90(4), 1663–1669.
https://doi.org/10.1785/0220190022

Ventosa, S., & Schimmel, M. (2023).
FastPCC: Fast phase cross-correlation algorithm for large seismic datasets.
*IEEE Transactions on Geoscience and Remote Sensing*, 61, 1–17.
https://doi.org/10.1109/TGRS.2023.3294302


Stacking daily windows
-----------------------

For each ``corr_duration``-long window, and for each configured filter, the
CCF (CC or PCC2) is computed and accumulated.  If ``keep_all = Y`` the
individual window CCFs are written to disk.  By default (``keep_days = Y``),
all windows for the day are stacked into a single daily CCF using either a
linear mean or Phase-Weighted Stack (PWS; Schimmel & Paulssen 1997):

.. note::

    PWS is provided as an experimental option.  It has not been
    systematically cross-validated.  Use with caution.

Schimmel, M., & Paulssen, H. (1997).
Noise reduction and detection of weak, coherent signals through
phase-weighted stacks.
*Geophysical Journal International*, 130(2), 497–505.
https://doi.org/10.1111/j.1365-246X.1997.tb05664.x


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
* |global.hpc|


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
                    # Always bandpass to the current filter band — required for
                    # both CC and PCC so that each filter_N produces a distinct
                    # result.  For CC this also replaces the old whitening==A
                    # guard which silently skipped the bandpass when whitening=N.
                    tmp = np.array([
                        bandpass(tr, freqmin=filterlow, freqmax=filterhigh,
                                 df=params.cc.cc_sampling_rate, corners=8)
                        for tr in _data
                    ])
                    if params.cc.clip_after_whiten:
                        tmp = winsorizing(tmp, params, input="timeseries")

                    if params.cc.cc_type_single_station_AC == "CC":
                        ffts = sf.fftn(tmp, [nfft, ], axes=[1, ])
                        if params.cc.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2,
                            axis=1)))
                        corr = myCorr2(ffts,
                                       np.ceil(params.cc.maxlag / dt),
                                       energy,
                                       single_station_pair_index_ac,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc.cc_normalisation)
                        del energy, ffts

                    elif params.cc.cc_type_single_station_AC == "PCC":
                        tmp_pcc = tmp
                        if params.cc.whitening != "N":
                            ffts_pcc = sf.fftn(tmp, [nfft, ], axes=[1, ])
                            whiten2(ffts_pcc, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                            tmp_pcc = np.real(sf.ifft(ffts_pcc, n=nfft, axis=1)[:, :tmp.shape[1]])
                        corr = pcc_xcorr(tmp_pcc,
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
                        # Band-pass filter to the current filter's frequency band
                        # before computing the analytic signal — this is what
                        # makes each filter_N produce a distinct PCC result.
                        # (The CC path achieves the same via whiten2's frequency
                        # window; PCC bypasses whitening so we apply it explicitly.)
                        _data_pcc = np.array([
                            bandpass(tr, freqmin=filterlow, freqmax=filterhigh,
                                     df=params.cc.cc_sampling_rate, corners=8)
                            for tr in _data
                        ])
                        if params.cc.whitening != "N":
                            # Optional: spectrally whiten within the band before
                            # AmpNorm so the phase signal is broadband across the
                            # full filter passband rather than dominated by the
                            # most energetic frequency in the band.
                            ffts_pcc = sf.fftn(_data_pcc, [nfft, ], axes=[1, ])
                            whiten2(ffts_pcc, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                            _data_pcc = np.real(sf.ifft(ffts_pcc, n=nfft, axis=1)[:, :len(_data_pcc[0])])
                        corr = pcc_xcorr(_data_pcc,
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
