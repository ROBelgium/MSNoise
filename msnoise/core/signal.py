"""MSNoise signal processing, preprocessing helpers, and stacking utilities."""

__all__ = [
    "check_and_phase_shift",
    "compute_wct_dtt",
    "find_segments",
    "getCoherence",
    "getGaps",
    "get_preprocessed_stream",
    "get_wavelet_type",
    "get_wct_avgcoh",
    "get_window",
    "make_same_length",
    "nextpow2",
    "preload_instrument_responses",
    "prepare_abs_positive_fft",
    "psd_df_rms",
    "psd_rms",
    "save_preprocessed_streams",
    "smoothCFS",
    "stack",
    "validate_stack_data",
    "wiener_filt",
    "winsorizing",
    "xwt",
]

import copy
import glob
import logging
import math
import os
import warnings

import numpy as np
import pandas as pd

from .config import get_config

def validate_stack_data(dataset, stack_type="reference"):
    """Validates stack data before processing

    Parameters:
        dataset: xarray Dataset to validate
        stack_type: Type of stack ("reference" or "moving") for error messages
    Returns:
        (is_valid, message) tuple
    """
    if dataset is None or not dataset.data_vars:
        return False, f"No data found for {stack_type} stack"

    if not hasattr(dataset, 'CCF'):
        return False, f"Missing CCF data in {stack_type} stack"

    data = dataset.CCF
    if data.size == 0:
        return False, f"Empty dataset in {stack_type} stack"

    nan_count = np.isnan(data.values).sum()
    total_points = data.values.size

    if nan_count == total_points:
        return False, f"{stack_type.capitalize()} stack contains only NaN values"

    if nan_count > 0:
        percent_nan = (nan_count / total_points) * 100
        return True, f"Warning: {stack_type.capitalize()} stack contains {percent_nan:.1f}% NaN values"

    return True, "OK"

# ============================================================


def nextpow2(x):
    """
    Returns the next power of 2 of `x`.

    :type x: int
    :param x: any value

    :rtype: int
    :returns: the next power of 2 of `x`
    """

    return np.ceil(np.log2(np.abs(x)))



def check_and_phase_shift(trace, taper_length=20.0):
    # TODO replace this hard coded taper length

    import scipy.fft as sf
    from scipy.fft import next_fast_len
    if trace.stats.npts < 4 * taper_length*trace.stats.sampling_rate:
        trace.data = np.zeros(trace.stats.npts)
        return trace

    dt = np.mod(trace.stats.starttime.datetime.microsecond*1.0e-6,
                trace.stats.delta)
    if (trace.stats.delta - dt) <= np.finfo(float).eps:
        dt = 0.
    if dt != 0.:
        if dt <= (trace.stats.delta / 2.):
            dt = -dt
#            direction = "left"
        else:
            dt = (trace.stats.delta - dt)
#            direction = "right"
        logging.debug("correcting time by %.6fs"%dt)
        trace.detrend(type="demean")
        trace.detrend(type="simple")
        trace.taper(max_percentage=None, max_length=1.0)

        n = next_fast_len(int(trace.stats.npts))
        FFTdata = sf.fft(trace.data, n=n)
        fftfreq = sf.fftfreq(n, d=trace.stats.delta)
        FFTdata = FFTdata * np.exp(1j * 2. * np.pi * fftfreq * dt)
        FFTdata = FFTdata.astype(np.complex64)
        sf.ifft(FFTdata, n=n, overwrite_x=True)
        trace.data = np.real(FFTdata[:len(trace.data)]).astype(float)
        trace.stats.starttime += dt
        del FFTdata, fftfreq
        # clean_scipy_cache()
        return trace
    else:
        return trace



def getGaps(stream, min_gap=None, max_gap=None):
    # Create shallow copy of the traces to be able to sort them later on.
    copied_traces = copy.copy(stream.traces)
    stream.sort()
    gap_list = []
    for _i in range(len(stream.traces) - 1):
        # skip traces with different network, station, location or channel
        if stream.traces[_i].id != stream.traces[_i + 1].id:
            continue
        # different sampling rates should always result in a gap or overlap
        if stream.traces[_i].stats.delta == stream.traces[_i + 1].stats.delta:
            flag = True
        else:
            flag = False
        stats = stream.traces[_i].stats
        stime = stats['endtime']
        etime = stream.traces[_i + 1].stats['starttime']
        delta = etime.timestamp - stime.timestamp
        # Check that any overlap is not larger than the trace coverage
        if delta < 0:
            temp = stream.traces[_i + 1].stats['endtime'].timestamp - \
                etime.timestamp
            if (delta * -1) > temp:
                delta = -1 * temp
        # Check gap/overlap criteria
        if min_gap and delta < min_gap:
            continue
        if max_gap and delta > max_gap:
            continue
        # Number of missing samples
        nsamples = int(round(math.fabs(delta) * stats['sampling_rate']))
        # skip if is equal to delta (1 / sampling rate)
        if flag and nsamples == 1:
            continue
        elif delta > 0:
            nsamples -= 1
        else:
            nsamples += 1
        gap_list.append([_i, _i+1,
                        stats['network'], stats['station'],
                        stats['location'], stats['channel'],
                        stime, etime, delta, nsamples])
    # Set the original traces to not alter the stream object.
    stream.traces = copied_traces
    del copied_traces
    return gap_list



def winsorizing(data, params, input="timeseries", nfft=0):
    """Clip (Winsorise) a 2-D data array in the time or frequency domain.

    Supports both one-shot sign-clipping (``winsorizing == -1``) and
    RMS-based clipping (``winsorizing > 0``).  When *input* is ``"fft"``
    the array is temporarily transformed back to the time domain, clipped,
    then re-transformed.

    :param data: 1-D or 2-D array of shape ``(n_traces, n_samples)``.
    :param params: MSNoise params object; must expose ``params.cc.winsorizing``.
    :param input: ``"timeseries"`` (default) or ``"fft"``.
    :param nfft: FFT length used when *input* is ``"fft"``; ignored otherwise.
    :returns: Clipped array (same shape as input).
    """
    import scipy.fft as sf
    input1D = False
    if len(data.shape) == 1:
        data = data.reshape(-1, data.shape[0])
        input1D = True
    if input == "fft":
        data = sf.ifftn(data, [nfft, ], axes=[1, ]).astype(float)
    for i in range(data.shape[0]):
        if params.cc.winsorizing == -1:
            np.sign(data[i], data[i])  # inplace
        elif params.cc.winsorizing != 0:
            rms = data[i].std() * params.cc.winsorizing
            np.clip(data[i], -rms, rms, data[i])  # inplace
    if input == "fft":
        data = sf.fftn(data, [nfft, ], axes=[1, ])
    if input1D:
        data = data[0]
    return data



def get_window(window="boxcar", half_win=3):
    """Return a normalised complex smoothing window for MWCS processing.

    :param window: ``"boxcar"`` (default) or ``"hanning"``.
    :param half_win: Half-width in samples (full window = ``2*half_win+1``).
    :returns: Complex numpy array of length ``2*half_win+1``, sum-normalised.
    """
    import scipy.signal
    window_len = 2 * half_win + 1
    if window == "boxcar":
        w = scipy.signal.windows.boxcar(window_len).astype("complex")
    else:
        w = scipy.signal.windows.hann(window_len).astype("complex")
    return w / window_len



def getCoherence(dcs, ds1, ds2):
    """Compute cross-coherence between two spectra.

    :param dcs: Cross-spectrum magnitudes (1-D array, length *n*).
    :param ds1: Auto-spectrum of signal 1 (1-D array, length *n*).
    :param ds2: Auto-spectrum of signal 2 (1-D array, length *n*).
    :returns: Complex coherence array of length *n*, clipped to ``|coh| <= 1``.
    """
    n = len(dcs)
    coh = np.zeros(n).astype("complex")
    valids = np.argwhere(np.logical_and(np.abs(ds1) > 0, np.abs(ds2) > 0))
    coh[valids] = dcs[valids] / (ds1[valids] * ds2[valids])
    coh[coh > (1.0 + 0j)] = 1.0 + 0j
    return coh



def prepare_abs_positive_fft(line, sampling_rate):
    """Return the positive-frequency part of the absolute FFT of *line*.

    :param line: 1-D signal array.
    :param sampling_rate: Sampling rate in Hz.
    :returns: ``(freq, val)`` - positive-frequency vector and absolute FFT values.
    """
    val = np.fft.fft(line)
    val = np.abs(val)
    freq = np.fft.fftfreq(len(line), 1.0 / sampling_rate)
    idx = np.where(freq >= 0)
    return freq[idx], val[idx]

# ============================================================


def _conv2(x, y, mode="same"):
    """2-D convolution using :func:`scipy.signal.convolve2d` with 180-degree rotations."""
    from scipy.signal import convolve2d
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)



def smoothCFS(cfs, scales, dt, ns, nt):
    """Smooth CWT coefficients along both time and scale axes.

    :param cfs: CWT coefficient array, shape ``(n_scales, n_times)``.
    :param scales: 1-D array of wavelet scales.
    :param dt: Sampling interval in seconds.
    :param ns: Length of the moving-average filter across scales.
    :param nt: Gaussian width parameter along time.
    :returns: Smoothed coefficient array, same shape as *cfs*.
    """
    import scipy.fft as sf
    N = cfs.shape[1]
    npad = sf.next_fast_len(N, real=True)
    omega = np.arange(1, np.fix(npad / 2) + 1, 1).tolist()
    omega = np.array(omega) * ((2 * np.pi) / npad)
    omega_save = -omega[int(np.fix((npad - 1) / 2)) - 1:0:-1]
    omega_2 = np.concatenate((0., omega), axis=None)
    omega_2 = np.concatenate((omega_2, omega_save), axis=None)
    omega = np.concatenate((omega_2, -omega[0]), axis=None)
    normscales = scales / dt
    for kk in range(0, cfs.shape[0]):
        F = np.exp(-nt * (normscales[kk] ** 2) * omega ** 2)
        smooth = np.fft.ifft(F * np.fft.fft(cfs[kk - 1], npad))
        cfs[kk - 1] = smooth[0:N]
    H = 1 / ns * np.ones((ns, 1))
    cfs = _conv2(cfs, H)
    return cfs



def get_wavelet_type(wavelet_type):
    """Return a :mod:`pycwt` wavelet object for the given type/parameter pair.

    :param wavelet_type: Tuple ``(name, param)`` or ``(name,)``.
        Supported names: ``"Morlet"``, ``"Paul"``, ``"DOG"``, ``"MexicanHat"``.
    :returns: Corresponding :mod:`pycwt` wavelet instance.
    """
    import pycwt as wavelet
    defaults = {"Morlet": 6, "Paul": 4, "DOG": 2, "MexicanHat": 2}
    name = wavelet_type[0]
    param = float(wavelet_type[1]) if len(wavelet_type) == 2 else defaults[name]
    if name == "Morlet":
        return wavelet.Morlet(param)
    elif name == "Paul":
        return wavelet.Paul(param)
    elif name == "DOG":
        return wavelet.DOG(param)
    elif name == "MexicanHat":
        return wavelet.MexicanHat()
    else:
        raise ValueError(f"Unknown wavelet type: {name!r}")



def prepare_ref_wct(trace_ref, fs, ns=3, nt=0.25, vpo=12,
                    freqmin=0.1, freqmax=8.0, nptsfreq=100, mother=None):
    """Pre-compute the CWT and smoothed power of the **reference** trace.

    Call this once for a fixed reference waveform (Mode A), then pass the
    returned tuple to :func:`apply_wct` for each current-day trace.  This
    avoids recomputing ``cwt_reference`` and ``smoothCFS(power_ref)`` — the
    two most expensive operations — on every iteration of the time loop.

    :param trace_ref: Reference signal (1-D numpy array).
    :param fs: Sampling frequency in Hz.
    :param ns: Scale-axis smoothing parameter.
    :param nt: Time-axis smoothing parameter.
    :param vpo: Voices-per-octave.
    :param freqmin: Lowest frequency of interest (Hz).
    :param freqmax: Highest frequency of interest (Hz).
    :param nptsfreq: Number of frequency points.
    :param mother: pycwt wavelet object (from :func:`get_wavelet_type`).
        If ``None``, defaults to ``Morlet(6)``.
    :returns: ``(cwt_ref, cfs1, scales, freqs, coi, invscales, dt, freqlim)``
        — an opaque tuple to pass directly to :func:`apply_wct`.
    """
    import pycwt as wavelet_lib
    if mother is None:
        mother = wavelet_lib.Morlet(6)
    nx = np.size(trace_ref)
    x_ref = np.transpose(trace_ref)
    dt = 1.0 / fs
    dj = 1.0 / vpo
    s0 = 2 * dt
    freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True)
    cwt_ref, scales, freqs, coi, _, _ = wavelet_lib.cwt(
        x_ref, dt, dj, s0, -1, mother, freqs=freqlim)
    scales_col = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1.0 / scales_col)
    power_ref = (invscales * abs(cwt_ref) ** 2).astype(complex)
    cfs1 = smoothCFS(power_ref, scales_col, dt, ns, nt)
    return cwt_ref, cfs1, scales_col, freqs, coi, invscales, dt, freqlim


def apply_wct(ref_wct_data, trace_current, ns=3, nt=0.25):
    """Compute WCT for one current-day trace given pre-computed reference data.

    :param ref_wct_data: Tuple returned by :func:`prepare_ref_wct`.
    :param trace_current: Current-day signal (1-D numpy array, same length as ref).
    :param ns: Scale-axis smoothing parameter (must match the value used in
        :func:`prepare_ref_wct`).
    :param nt: Time-axis smoothing parameter (idem).
    :returns: ``(WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi)``
    """
    import pycwt as wavelet_lib
    cwt_ref, cfs1, scales_col, freqs, coi, invscales, dt, freqlim = ref_wct_data
    x_cur = np.transpose(trace_current)
    # Re-derive wavelet params consistent with the reference CWT
    dj = 1.0 / (scales_col.shape[0] - 1) if scales_col.shape[0] > 1 else 1.0
    s0 = 2 * dt
    nx = np.size(trace_current)
    # Use the same freqlim from the reference call so scales are identical
    cwt_cur, _, _, _, _, _ = wavelet_lib.cwt(
        x_cur, dt, dj, s0, -1, wavelet_lib.Morlet(6), freqs=freqlim)
    # Recompute invscales for current length (same as ref if length unchanged)
    inv_cur = np.kron(np.ones((1, nx)), 1.0 / scales_col)
    power_cur = (inv_cur * abs(cwt_cur) ** 2).astype(complex)
    crossCFS_raw = cwt_ref * np.conj(cwt_cur)
    WXamp = abs(crossCFS_raw)
    cross_spectrum = (inv_cur * crossCFS_raw).astype(complex)
    cfs2 = smoothCFS(power_cur, scales_col, dt, ns, nt)
    crossCFS = smoothCFS(cross_spectrum, scales_col, dt, ns, nt)
    valid_mask = (cfs1 > 0) & (cfs2 > 0)
    WXspec = np.full_like(crossCFS, np.nan, dtype=complex)
    Wcoh   = np.full_like(crossCFS, np.nan)
    WXspec[valid_mask] = crossCFS[valid_mask] / (
        np.sqrt(cfs1[valid_mask]) * np.sqrt(cfs2[valid_mask]))
    Wcoh[valid_mask] = (abs(crossCFS[valid_mask]) ** 2
                        / (cfs1[valid_mask] * cfs2[valid_mask]))
    WXangle = np.angle(WXspec)
    Wcoh = np.clip(Wcoh, 0.0, 1.0)
    pp2 = np.array([[2 * np.pi * f] for f in freqs])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)
    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi


def xwt(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12,
        freqmin=0.1, freqmax=8.0, nptsfreq=100, wavelet_type=("Morlet", 6.)):
    """Wavelet Coherence Transform (WCT) between two time series.

    Convenience wrapper around :func:`prepare_ref_wct` + :func:`apply_wct`.
    Use those two functions directly in hot loops (Mode A fixed-REF) to avoid
    recomputing the reference CWT on every call.

    :param trace_ref: Reference signal (1-D array).
    :param trace_current: Current signal (1-D array, same length).
    :param fs: Sampling frequency in Hz.
    :param ns: Scale-axis smoothing parameter.
    :param nt: Time-axis smoothing parameter.
    :param vpo: Voices-per-octave; higher = finer scale resolution.
    :param freqmin: Lowest frequency of interest (Hz).
    :param freqmax: Highest frequency of interest (Hz).
    :param nptsfreq: Number of frequency points between *freqmin* and *freqmax*.
    :param wavelet_type: ``(name, param)`` tuple passed to :func:`get_wavelet_type`.
    :returns: ``(WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi)``
    """
    mother = get_wavelet_type(wavelet_type)
    ref_data = prepare_ref_wct(trace_ref, fs, ns, nt, vpo,
                               freqmin, freqmax, nptsfreq, mother)
    return apply_wct(ref_data, trace_current, ns, nt)



def compute_wct_dtt(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20,
                    mincoh=0.5, maxdt=0.2, min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Compute dv/v and associated errors from wavelet coherence transform results.

    :param freqs: Frequency values from the WCT.
    :param tvec: Time axis.
    :param WXamp: Cross-wavelet amplitude array (freqs × taxis).
    :param Wcoh: Wavelet coherence array (freqs × taxis).
    :param delta_t: Time delay array (freqs × taxis).
    :param lag_min: Minimum coda lag in seconds.
    :param coda_cycles: Number of periods to use as coda window width.
    :param mincoh: Minimum coherence threshold.
    :param maxdt: Maximum allowed time delay.
    :param min_nonzero: Minimum fraction of valid (non-zero weight) samples required.
    :param freqmin: Lower frequency bound for regression.
    :param freqmax: Upper frequency bound for regression.
    :returns: Tuple of (dt/t, err, weighting_function).
    """
    import warnings
    from scipy.optimize import OptimizeWarning
    from obspy.signal.regression import linear_regression

    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
    dvv = np.zeros(len(inx[0]))
    err = np.zeros(len(inx[0]))

    weight_func = np.log(np.abs(WXamp)) / np.log(np.abs(WXamp)).max()
    zero_idx = np.where((Wcoh < mincoh) | (delta_t > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / weight_func.max()
    wf[zero_idx] = 0

    problematic_freqs = []
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)
        tindex = np.where(
            ((tvec >= -lag_max) & (tvec <= -lag_min)) |
            ((tvec >= lag_min) & (tvec <= lag_max))
        )[0]
        if len(tvec) > 2:
            if not np.any(delta_t[ifreq]):
                continue
            delta_t[ifreq][tindex] = np.nan_to_num(delta_t[ifreq][tindex])
            w = wf[ifreq]
            nzc_perc = np.count_nonzero(w[tindex] > 0) / len(tindex)
            if nzc_perc >= min_nonzero:
                with warnings.catch_warnings(record=True) as w_catcher:
                    warnings.simplefilter("always", OptimizeWarning)
                    m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex],
                                              w[tindex], intercept_origin=True)
                    if any(issubclass(warning.category, OptimizeWarning)
                           for warning in w_catcher):
                        problematic_freqs.append(freqs[ifreq])
                dvv[ii], err[ii] = m, em
            else:
                dvv[ii], err[ii] = np.nan, np.nan
    if problematic_freqs:
        logging.warning(
            f"Covariance issues at {min(problematic_freqs):.2f}-{max(problematic_freqs):.2f} Hz: "
            f"consider adjusting min_nonzero={min_nonzero}, mincoh={mincoh}, "
            f"maxdt={maxdt}, coda_cycles={coda_cycles}"
        )
    return dvv, err, wf



def get_wct_avgcoh(freqs, tvec, wcoh, freqmin, freqmax, lag_min=5, coda_cycles=20):
    """
    Calculate average wavelet coherence over a frequency range and coda window.

    :param freqs: Frequency array.
    :param tvec: Time axis.
    :param wcoh: Wavelet coherence array (freqs × taxis).
    :param freqmin: Lower frequency bound.
    :param freqmax: Upper frequency bound.
    :param lag_min: Minimum coda lag in seconds.
    :param coda_cycles: Number of periods to use as coda window width.
    :returns: Average coherence per frequency bin within [freqmin, freqmax].
    """
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
    coh = np.zeros(inx[0].shape)
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)
        tindex = np.where(
            ((tvec >= -lag_max) & (tvec <= -lag_min)) |
            ((tvec >= lag_min) & (tvec <= lag_max))
        )[0]
        if len(tvec) > 2:
            if not np.any(wcoh[ifreq]) or wcoh[ifreq][tindex].size == 0:
                coh[ii] = np.nan
                continue
            coh[ii] = np.nanmean(np.abs(wcoh[ifreq][tindex]))
        else:
            coh[ii] = np.nan
    return coh

# ============================================================


def preload_instrument_responses(session, return_format="dataframe"):
    """
    This function preloads all instrument responses from ``response_path``
    and stores the seed ids, start and end dates, and paz for every channel
    in a DataFrame. Any file readable by obspy's read_inventory will be processed.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type return_format: str
    :param return_format: The format of the returned object, either
        ``dataframe`` or ``inventory``.

    :rtype: :class:`~pandas.DataFrame` or :class:`~obspy.core.inventory.inventory.Inventory`
    :returns: A table containing all channels with the time of operation and
        poles and zeros (DataFrame), or an obspy Inventory object.

    """
    from obspy.core.inventory import Inventory
    from obspy import read_inventory, UTCDateTime
    logging.debug('Preloading instrument response')
    files = glob.glob(os.path.join(get_config(session, 'response_path'), "*"))
    channels = []
    all_inv = Inventory()
    for file in files:
        logging.debug("Processing %s" % file)
        try:
            inv = read_inventory(file)

            if return_format == "inventory":
                all_inv += inv
                continue

            for net in inv.networks:
                for sta in net.stations:
                    for cha in sta.channels:
                        seed_id = "%s.%s.%s.%s" % (net.code, sta.code,
                                                   cha.location_code,
                                                   cha.code)
                        pzdict = {}
                        try:
                            resp = inv.get_response(seed_id, cha.start_date + 10)
                            polezerostage = resp.get_paz()
                        except Exception as e:
                            logging.warning(
                                'Failed to get PAZ for SEED ID "%s", this '
                                'SEED ID will have an empty dictionary '
                                'for Poles and Zeros '
                                'information (Error message: %s).' % (
                                    seed_id, str(e)))
                        else:
                            totalsensitivity = resp.instrument_sensitivity
                            pzdict['poles'] = polezerostage.poles
                            pzdict['zeros'] = polezerostage.zeros
                            pzdict['gain'] = polezerostage.normalization_factor
                            pzdict['sensitivity'] = totalsensitivity.value
                        lat = cha.latitude
                        lon = cha.longitude
                        elevation = cha.elevation
                        if lat is None or lon is None or elevation is None:
                            lat = sta.latitude
                            lon = sta.longitude
                            elevation = sta.elevation
                        if lat is None or lon is None or elevation is None:
                            logging.error(
                                'Failed to look up coordinates for SEED '
                                'ID: %s' % seed_id)
                        channels.append([seed_id, cha.start_date,
                                         cha.end_date or UTCDateTime(),
                                         pzdict, lat, lon, elevation])

        except Exception as e:
            logging.error('Failed to process file %s: %s' % (file, str(e)))


    logging.debug('Finished Loading instrument responses')
    if return_format == "inventory":
        return all_inv

    if return_format == "dataframe":
        channels = pd.DataFrame(channels, columns=["channel_id", "start_date",
                                                   "end_date", "paz",
                                                   "latitude", "longitude", "elevation"],)
        return channels







def save_preprocessed_streams(stream, output_dir, step_name, goal_day):
    """Write preprocessed traces to per-station files.

    Output layout::

        <output_dir>/<step_name>/_output/<goal_day>/<NET.STA.LOC>.mseed

    One file per station (all channels for that station in the same file).
    This is concurrency-safe: each station writes its own file with no
    shared state between parallel workers.

    :param stream: :class:`~obspy.core.stream.Stream` to write.
    :param output_dir: Base output directory (``params.global_.output_folder``).
    :param step_name: Workflow step name (e.g. ``"preprocess_1"``).
    :param goal_day: Processing date string (``YYYY-MM-DD``).
    :returns: List of written file paths (one per station).
    """
    from obspy import Stream as _Stream
    day_dir = os.path.join(output_dir, step_name, "_output", goal_day)
    os.makedirs(day_dir, exist_ok=True)

    # Group traces by NET.STA.LOC
    by_station = {}
    for tr in stream:
        sid = f"{tr.stats.network}.{tr.stats.station}.{tr.stats.location}"
        by_station.setdefault(sid, _Stream())
        tr.data = tr.data.astype(np.float32)
        by_station[sid].append(tr)

    saved = []
    for sid, st in by_station.items():
        output_path = os.path.join(day_dir, f"{sid}.mseed")
        st.write(output_path, format="MSEED")
        saved.append(output_path)

    return saved


def get_preprocessed_stream(output_dir, step_name, goal_day, stations):
    """Read per-station preprocessed files and return a merged Stream.

    Counterpart to :func:`save_preprocessed_streams`.  Reads only the
    station files needed for *stations* (a list of ``"NET.STA.LOC"``
    strings) and returns them merged into a single
    :class:`~obspy.core.stream.Stream`.

    Missing station files are silently skipped (the cc worker checks
    for empty streams downstream).

    :param output_dir: Base output directory.
    :param step_name: Workflow step name (e.g. ``"preprocess_1"``).
    :param goal_day: Processing date string (``YYYY-MM-DD``).
    :param stations: List of ``"NET.STA.LOC"`` strings to load.
    :returns: :class:`~obspy.core.stream.Stream`.
    """
    from obspy import Stream as _Stream, read as _read
    day_dir = os.path.join(output_dir, step_name, "_output", goal_day)
    merged = _Stream()
    for sid in stations:
        fpath = os.path.join(day_dir, f"{sid}.mseed")
        if os.path.isfile(fpath):
            merged += _read(fpath)
        else:
            import logging as _log
            _log.getLogger("msnoise.signal").debug(
                f"Preprocessed file not found: {fpath}"
            )
    return merged


# ── PSD computation helpers ──────────────────────────────────────────────────

def psd_rms(s, f):
    """Compute RMS from a power spectrum array and frequency array.

    :param s: Power spectral density values (1-D array).
    :param f: Frequency values (1-D array, same length as *s*).
    :returns: Float - square-root of the integrated power.
    """
    return np.sqrt(np.trapezoid(s, f))


def psd_df_rms(d, freqs, output="VEL"):
    """Compute per-frequency-band RMS from PSD data.

    :param d: :class:`xarray.Dataset` with a ``PSD`` variable and dims
        ``(times, periods)``, as returned by :func:`~msnoise.core.io.xr_load_psd`.
    :param freqs: List of ``(fmin, fmax)`` tuples defining frequency bands.
    :param output: Physical unit - ``"VEL"`` (default), ``"ACC"``, or ``"DISP"``.
    :returns: :class:`xarray.Dataset` with one ``RMS`` variable and dims
        ``(times, bands)``, ready for :func:`~msnoise.core.io.xr_save_rms`.
    """
    import xarray as xr_mod
    da = d.PSD.load()  # materialise — we need numpy ops
    periods = da.coords["periods"].values.astype(float)
    times   = da.coords["times"].values

    rms_bands = {}
    for fmin, fmax in freqs:
        pmin = 1.0 / fmax
        pmax = 1.0 / fmin
        ix = np.where((periods >= pmin) & (periods <= pmax))[0]
        if ix.size == 0:
            continue
        f    = periods[ix]
        w2f  = 2.0 * np.pi * f
        amp  = 10.0 ** (da.values[:, ix] / 10.0)   # (times, n_periods)
        if output == "ACC":
            vals = np.sqrt(np.trapezoid(amp, f, axis=1))
        elif output == "VEL":
            vamp = amp / w2f ** 2
            vals = np.sqrt(np.trapezoid(vamp, f, axis=1))
        else:
            vamp = amp / w2f ** 2
            damp = vamp / w2f ** 2
            vals = np.sqrt(np.trapezoid(damp, f, axis=1))
        rms_bands[f"{fmin:.1f}-{fmax:.1f}"] = vals

    bands = list(rms_bands.keys())
    data  = np.column_stack(list(rms_bands.values())) if bands else np.empty((len(times), 0))
    da_rms = xr_mod.DataArray(
        data, coords=[times, bands], dims=["times", "bands"], name="RMS"
    )
    return da_rms.to_dataset()

def make_same_length(st):
    """
    This function takes a stream of equal sampling rate and makes sure that all
    channels have the same length and the same gaps.
    """
    warnings.warn("make_same_length() is deprecated and will be removed in a future MSNoise release.", DeprecationWarning, stacklevel=2)
    from obspy import Stream
    # Merge traces
    st.merge()

    # Initialize arrays to be filled with start+endtimes of all traces
    starttimes = []
    endtimes = []

    # Loop over all traces of the stream
    for tr in st:
        # Force conversion to masked arrays
        if not np.ma.count_masked(tr.data):
            tr.data = np.ma.array(tr.data, mask=False)
        # Read out start+endtimes of traces to trim
        starttimes.append(tr.stats.starttime)
        endtimes.append(tr.stats.endtime)

    # trim stream to common starttimes
    if max(starttimes) >= min(endtimes):
        return Stream()
    st.trim(max(starttimes), min(endtimes))

    # get the mask of all traces, i.e. the parts where at least one trace has
    # a gap

    if len(st) < 2:
        return st

    masks=[]
    for tr in st:
        masks.append(tr.data.mask)
    mask =  np.any(masks,axis=0)

    # apply the mask to all traces
    for tr in st:
        tr.data.mask = mask

    st = st.split()
    return st

# ============================================================


def stack(data, stack_method="linear", pws_timegate=10.0, pws_power=2,
          goal_sampling_rate=20.0):
    """
    :type data: :class:`numpy.ndarray`
    :param data: the data to stack, each row being one CCF
    :type stack_method: str
    :param stack_method: either ``linear``: average of all CCF or ``pws`` to
        compute the phase weigthed stack. If ``pws`` is selected,
        the function expects the ``pws_timegate`` and ``pws_power``.
    :type pws_timegate: float
    :param pws_timegate: PWS time gate in seconds. Width of the smoothing
         window to convolve with the PWS spectrum.
    :type pws_power: float
    :param pws_power: Power of the PWS weights to be applied to the CCF stack.
    :type goal_sampling_rate: float
    :param goal_sampling_rate: Sampling rate of the CCF array submitted
    :rtype: :class:`numpy.array`
    :return: the stacked CCF.
    """
    if len(data) == 0:
        logging.debug("No data to stack.")
        return []
    data = data[~np.isnan(data).any(axis=1)]
    sanitize = False
    # TODO clean sanitize function, add param to config and make sure not to
    # kill the data[i] if all data are corrcoeff >0.9 (either very stable
    # corr or autocorr, then this sanitize should not occur.
    if len(data) != 1 and sanitize:
        threshold = 0.99
        corr = data.mean(axis=0)
        corrcoefs = np.array([np.corrcoef(di, corr)[1][0] for di in data])
        toolarge = np.where(corrcoefs >= threshold)[0]
        if len(toolarge):
            data = data[np.where(corrcoefs <= threshold)[0]]

    if len(data) == 0:
        return []
    if stack_method == "linear":
        # logging.debug("Doing a linear stack")
        corr = data.mean(axis=0)

    elif stack_method == "pws":
        import scipy.signal as ss
        # logging.debug("Doing a PWS stack")
        corr = np.zeros(data.shape[1], dtype='f8')
        phasestack = np.zeros(data.shape[1], dtype='c8')
        for i in range(data.shape[0]):
            data[i] -= data[i].mean()
        for c in data:
            phase = np.angle(ss.hilbert(c))
            phasestack.real += np.cos(phase)
            phasestack.imag += np.sin(phase)
        coh = 1. / data.shape[0] * np.abs(phasestack)

        timegate_samples = int(pws_timegate * goal_sampling_rate)
        coh = np.convolve(ss.windows.boxcar(timegate_samples) /
                          timegate_samples, coh, 'same')
        coh = np.power(coh, pws_power)
        for c in data:
            corr += c * coh
        corr /= data.shape[0]

    return corr

# ── Wiener filter helpers (moved from msnoise/wiener.py) ──────────────────

def find_segments(data, gap_threshold):
    """Identify continuous non-NaN segments in an xarray DataArray.

    :param data: 2-D xarray DataArray (times × lags).
    :param gap_threshold: Maximum index gap before treating as a new segment.
    :returns: List of lists of row indices forming each continuous segment.
    """
    current_segment = []
    continuous_segments = []
    prev_idx = None

    for i in range(data.shape[0]):
        if not data[i, :].isnull().all():
            if prev_idx is not None and (i - prev_idx > gap_threshold):
                continuous_segments.append(current_segment)
                current_segment = []
            current_segment.append(i)
            prev_idx = i
        else:
            prev_idx = None

    if current_segment:
        continuous_segments.append(current_segment)

    return continuous_segments


def wiener_filt(data, M, N, gap_threshold):
    """Apply a 2-D Wiener filter to the CCF dataset, segment by segment.

    Operates only on continuous (non-NaN) segments of the time axis to avoid
    smearing across data gaps.

    :param data: xarray Dataset containing a ``CCF`` variable (times × taxis).
    :param M: Wiener filter window size along the time axis.
    :param N: Wiener filter window size along the lag axis.
    :param gap_threshold: Passed to :func:`find_segments`.
    :returns: Copy of *data* with the filtered ``CCF`` variable.
    """
    from scipy.signal import wiener as _wiener
    import xarray as xr

    ccfs = data["CCF"]
    segments = find_segments(ccfs, gap_threshold)

    # Work on a single numpy array (one allocation) rather than two deep
    # DataArray/Dataset copies; peak memory is now ~2× CCF instead of ~3×.
    filtered_values = ccfs.values.copy()
    for segment in segments:
        filtered_values[segment, :] = _wiener(ccfs.values[segment, :], (M, N))

    filtered_data = data.copy(deep=False)   # shallow: copies Dataset shell only
    filtered_data["CCF"] = xr.DataArray(
        filtered_values, coords=ccfs.coords, dims=ccfs.dims, attrs=ccfs.attrs
    )
    del filtered_values
    return filtered_data
