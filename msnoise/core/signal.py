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

    # Use xarray .isnull() so this works on both numpy and dask-backed arrays.
    # data.values on a dask array would load the entire dataset into RAM.
    nan_count = int(data.isnull().sum().values)
    total_points = data.size

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
        # scipy's overwrite_x is only a hint, never a guarantee: relying on
        # it and dropping the return value would silently leave FFTdata in
        # the frequency domain on any build that declines to overwrite.
        FFTdata = sf.ifft(FFTdata, n=n)
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
        # np.real(), not .astype(float): casting a complex array to float
        # discards the imaginary part behind a ComplexWarning.  The input is
        # Hermitian (whiten2 enforces it) so the imaginary part is numerical
        # noise, but the intent should be explicit.
        # NOTE: the array is nfft long while the trace is only npts long, so
        # the trailing zero-padding slightly dilutes the RMS below and makes
        # the clip threshold marginally tighter than the clip_after_whiten=N
        # path.  next_fast_len keeps nfft/npts under ~1.05, so the effect is
        # under ~2.5 % on the threshold.
        data = np.real(sf.ifftn(data, [nfft, ], axes=[1, ]))
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

    .. note::

        The normalisation used to be ``w / window_len``, which sum-normalises
        a boxcar but leaves a Hann window summing to about 0.45.  The
        constant cancels exactly in the MWCS coherence (numerator and
        denominator are smoothed with the same window) and in the
        origin-forced WLS slope and its error, so this is a consistency fix
        with :func:`~msnoise.core.compute.smooth`, not a numerical change.
    """
    import scipy.signal
    window_len = 2 * half_win + 1
    if window == "boxcar":
        w = scipy.signal.windows.boxcar(window_len).astype("complex")
    else:
        w = scipy.signal.windows.hann(window_len).astype("complex")
    return w / w.sum()



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

    Time-axis smoothing is a Gaussian applied in the Fourier domain, with a
    width proportional to the (dt-normalised) scale of *that* row.  Scale-axis
    smoothing is a length-``ns`` moving average across rows.  Port of MATLAB
    ``wcoherence``'s ``smoothCFS``.

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
    for kk in range(cfs.shape[0]):
        F = np.exp(-nt * (normscales[kk] ** 2) * omega ** 2)
        smooth = np.fft.ifft(F * np.fft.fft(cfs[kk], npad))
        cfs[kk] = smooth[0:N]
    H = 1 / ns * np.ones((ns, 1))
    cfs = _conv2(cfs, H)
    return cfs




# ── Wavelet classes and CWT ──────────────────────────────────────────────────
# Self-contained implementation of the mother wavelets and CWT used by the
# WCT pipeline (s08_compute_wct).  Algorithm: :footcite:t:`TorrenceCompo1998`,
# Tables 1–2.

class _Morlet:
    """Morlet wavelet (:footcite:t:`TorrenceCompo1998`, Table 1, row 1).

    :param f0: Central angular frequency. Default 6.
    """
    name = "Morlet"

    def __init__(self, f0=6):
        self.f0 = float(f0)
        self.dofmin = 2
        if f0 == 6:
            self.cdelta, self.gamma, self.deltaj0 = 0.776, 2.32, 0.60
        else:
            self.cdelta = self.gamma = self.deltaj0 = -1

    def psi_ft(self, f):
        """Fourier transform of the Morlet wavelet."""
        return np.pi ** -0.25 * np.exp(-0.5 * (f - self.f0) ** 2)

    def flambda(self):
        """Fourier wavelength (T&C eq. 1)."""
        return (4 * np.pi) / (self.f0 + np.sqrt(2 + self.f0 ** 2))

    def coi(self):
        """e-Folding time (T&C Table 1)."""
        return 1.0 / np.sqrt(2)

    def smooth(self, W, dt, dj, scales):
        """Smooth CWT coefficients (time + scale axes) for coherence."""
        from scipy.signal import convolve2d
        m, n = W.shape
        N = int(2 ** np.ceil(np.log2(n)))
        k = 2 * np.pi * np.fft.fftfreq(N)
        snorm = scales / dt
        F = np.exp(-0.5 * (snorm[:, np.newaxis] ** 2) * k ** 2)
        T = np.fft.ifft(F * np.fft.fft(W, n=N, axis=1), axis=1)[:, :n]
        if np.isreal(W).all():
            T = T.real
        wsize = self.deltaj0 / dj * 2
        win = np.ones(max(1, int(np.round(wsize))))
        win /= win.sum()
        T = convolve2d(T, win[:, np.newaxis], "same")
        return T


class _Paul:
    """Paul wavelet of order *m* (:footcite:t:`TorrenceCompo1998`, Table 1, row 2)."""
    name = "Paul"

    def __init__(self, m=4):
        self.m = int(m)
        self.dofmin = 2
        if m == 4:
            self.cdelta, self.gamma, self.deltaj0 = 1.132, 1.17, 1.50
        else:
            self.cdelta = self.gamma = self.deltaj0 = -1

    def psi_ft(self, f):
        m = self.m
        norm = 2 ** m / np.sqrt(m * np.prod(np.arange(2, 2 * m, dtype=float)))
        return norm * f ** m * np.exp(-f) * (f > 0)

    def flambda(self):
        return 4 * np.pi / (2 * self.m + 1)

    def coi(self):
        return np.sqrt(2)


class _DOG:
    """Derivative-of-Gaussian wavelet of order *m* (:footcite:t:`TorrenceCompo1998`, Table 1, row 3)."""
    name = "DOG"

    def __init__(self, m=2):
        from scipy.special import gamma as _gamma
        self.m = int(m)
        self.dofmin = 1
        self._norm = 1.0 / np.sqrt(_gamma(m + 0.5))
        if m == 2:
            self.cdelta, self.gamma, self.deltaj0 = 3.541, 1.43, 1.40
        elif m == 6:
            self.cdelta, self.gamma, self.deltaj0 = 1.966, 1.37, 0.97
        else:
            self.cdelta = self.gamma = self.deltaj0 = -1

    def psi_ft(self, f):
        return -(1j ** self.m) * self._norm * f ** self.m * np.exp(-0.5 * f ** 2)

    def flambda(self):
        return 2 * np.pi / np.sqrt(self.m + 0.5)

    def coi(self):
        return 1.0 / np.sqrt(2)


class _MexicanHat(_DOG):
    """Mexican hat wavelet (DOG with m=2)."""
    name = "MexicanHat"

    def __init__(self):
        super().__init__(m=2)


def _cwt(signal, dt, dj=1 / 12, s0=-1, J=-1, wavelet=None, freqs=None):
    """Continuous wavelet transform (:footcite:t:`TorrenceCompo1998`).

    Uses FFT-domain convolution for efficiency.

    :param signal: 1-D input array.
    :param dt: Sampling interval (seconds).
    :param dj: Scale spacing. Default 1/12.
    :param s0: Smallest scale. Default ``2*dt / wavelet.flambda()``.
    :param J: Number of scales minus one. Default derived from signal length.
    :param wavelet: Mother wavelet instance. Defaults to ``_Morlet(6)``.
    :param freqs: Optional custom frequency array (Hz); overrides dj/s0/J.
    :returns: ``(W, sj, freqs, coi, signal_ft, ftfreqs)``
    """
    if wavelet is None:
        wavelet = _Morlet()
    n0 = len(signal)
    if freqs is not None:
        sj = 1.0 / (wavelet.flambda() * np.asarray(freqs))
    else:
        if s0 == -1:
            s0 = 2 * dt / wavelet.flambda()
        if J == -1:
            J = int(np.round(np.log2(n0 * dt / s0) / dj))
        sj = s0 * 2 ** (np.arange(0, J + 1) * dj)
        freqs = 1.0 / (wavelet.flambda() * sj)
    N = int(2 ** np.ceil(np.log2(n0)))
    signal_ft = np.fft.fft(signal, n=N)
    ftfreqs = 2 * np.pi * np.fft.fftfreq(N, dt)
    sj_col = sj[:, np.newaxis]
    psi_ft_bar = ((sj_col * ftfreqs[1] * N) ** 0.5
                  * np.conj(wavelet.psi_ft(sj_col * ftfreqs)))
    W = np.fft.ifft(signal_ft * psi_ft_bar, axis=1)
    sel = ~np.isnan(W).all(axis=1)
    sj, freqs, W = sj[sel], freqs[sel], W[sel, :n0]
    coi = (wavelet.flambda() * wavelet.coi() * dt
           * (n0 / 2 - np.abs(np.arange(n0) - (n0 - 1) / 2)))
    ft_out = signal_ft[1: N // 2] / N ** 0.5
    ftfreqs_out = ftfreqs[1: N // 2] / (2 * np.pi)
    return W, sj, freqs, coi, ft_out, ftfreqs_out

def get_wavelet_type(wavelet_type):
    """Return an internal wavelet object for the given type/parameter pair.

    :param wavelet_type: Tuple ``(name, param)`` or ``(name,)``.
        Supported names: ``"Morlet"``, ``"Paul"``, ``"DOG"``, ``"MexicanHat"``.
    :returns: Corresponding wavelet instance (_Morlet, _Paul, _DOG, _MexicanHat).
    """
    defaults = {"Morlet": 6, "Paul": 4, "DOG": 2, "MexicanHat": 2}
    name = wavelet_type[0]
    if name not in defaults:
        raise ValueError(f"Unknown wavelet type: {name!r}")
    param = float(wavelet_type[1]) if len(wavelet_type) == 2 else defaults[name]
    if name == "Morlet":
        return _Morlet(param)
    elif name == "Paul":
        return _Paul(int(param))
    elif name == "DOG":
        return _DOG(int(param))
    elif name == "MexicanHat":
        return _MexicanHat()



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
    :param mother: Wavelet instance from :func:`get_wavelet_type`.
        If ``None``, defaults to ``_Morlet(6)``.
    :returns: ``(cwt_ref, cfs1, scales, freqs, coi, invscales, dt, freqlim,
        mother)`` — an opaque tuple to pass directly to :func:`apply_wct`.
        Treat it as opaque: unpack by index, not by fixed-width tuple
        assignment, so future additions do not break callers.
    """
    if mother is None:
        mother = _Morlet(6)
    nx = np.size(trace_ref)
    x_ref = np.transpose(trace_ref)
    dt = 1.0 / fs
    dj = 1.0 / vpo
    s0 = 2 * dt
    freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True)
    cwt_ref, scales, freqs, coi, _, _ = _cwt(
        x_ref, dt, dj, s0, -1, mother, freqs=freqlim)
    scales_col = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1.0 / scales_col)
    power_ref = (invscales * abs(cwt_ref) ** 2).astype(complex)
    cfs1 = smoothCFS(power_ref, scales_col, dt, ns, nt)
    return (cwt_ref, cfs1, scales_col, freqs, coi, invscales, dt, freqlim,
            mother)


def apply_wct(ref_wct_data, trace_current, ns=3, nt=0.25):
    """Compute WCT for one current-day trace given pre-computed reference data.

    Implements the wavelet cross-spectrum following :footcite:t:`Grinsted2004`.
    The time-delay ``WXdt = phase(WXspec) / (2πf)`` follows
    :footcite:t:`Mao2020` (their eq. 1).

    :param ref_wct_data: Tuple returned by :func:`prepare_ref_wct`.
    :param trace_current: Current-day signal (1-D numpy array, same length as ref).
    :param ns: Scale-axis smoothing parameter (must match the value used in
        :func:`prepare_ref_wct`).
    :param nt: Time-axis smoothing parameter (idem).
    :returns: ``(WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi)``
    """
    (cwt_ref, cfs1, scales_col, freqs, coi, invscales, dt, freqlim,
     mother) = ref_wct_data
    x_cur = np.transpose(trace_current)
    nx = np.size(trace_current)
    # ``freqs=freqlim`` overrides dj/s0/J inside _cwt, so they are unused here.
    # The mother wavelet MUST be the same one used for the reference: the
    # scales are derived as 1/(mother.flambda()*f), so a mismatch would make
    # cwt_ref and cwt_cur live on different scale grids.
    cwt_cur, _, _, _, _, _ = _cwt(
        x_cur, dt, 1.0 / 12, -1, -1, mother, freqs=freqlim)
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
    """Wavelet Coherence Transform (WCT) between two time series
    (:footcite:t:`Grinsted2004`; traveltime estimation: :footcite:t:`Mao2020`).

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
    """Compute dv/v and associated errors from WCT results following
    :footcite:t:`Mao2020`.

    For each frequency band, fits a weighted linear regression
    ``delta_t(t) = -(dv/v) * t`` using log-amplitude weights from the
    cross-wavelet spectrum (their eq. 3–4).


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


def _band_label(v):
    """Format a band edge, preferring the historical one-decimal form.

    Returns ``"%.1f" % v`` when that is exact (so existing band names such as
    ``"1.0-20.0"`` are unchanged) and a lossless ``"%g"`` otherwise.
    """
    one_dp = f"{v:.1f}"
    return one_dp if float(one_dp) == float(v) else f"{v:g}"


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
        # ObsPy's PPSD reports an acceleration PSD **per unit frequency**
        # against a PERIOD axis.  Both the omega factor and the integration
        # variable must therefore be frequencies:
        #     rms^2 = int S(f) df,   S_vel = S_acc / (2*pi*f)^2
        # The previous code passed ``periods`` in as ``f``, so it used
        # omega = 2*pi*T instead of 2*pi/T and integrated dT instead of df.
        fq    = 1.0 / periods[ix]                  # Hz
        order = np.argsort(fq)                     # trapezoid needs ascending x
        fq    = fq[order]
        w2f   = 2.0 * np.pi * fq
        amp   = 10.0 ** (da.values[:, ix][:, order] / 10.0)   # (times, n_freq)
        if output == "ACC":
            vals = np.sqrt(np.trapezoid(amp, fq, axis=1))
        elif output == "VEL":
            vamp = amp / w2f ** 2
            vals = np.sqrt(np.trapezoid(vamp, fq, axis=1))
        else:
            vamp = amp / w2f ** 2
            damp = vamp / w2f ** 2
            vals = np.sqrt(np.trapezoid(damp, fq, axis=1))
        # "%.1f" is lossy below 0.05 Hz: (0.05, 0.1) and (0.06, 0.1) both
        # label as "0.1-0.1" and the second silently overwrites the first.
        # Keep the historical one-decimal label whenever it round-trips
        # exactly (all round-numbered bands, e.g. "1.0-20.0") and fall back
        # to a lossless one otherwise.
        rms_bands[f"{_band_label(fmin)}-{_band_label(fmax)}"] = vals

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



def _morlet_wavelet(M, s, w=5.0):
    """Complex Morlet wavelet (L2-normalised) for use in :func:`tfpws_stack`.

    The abscissa is the **sample offset** from the centre of the window, so
    the wavelet has a centre frequency of ``w / (2*pi*s)`` cycles per sample
    and a Gaussian envelope of standard deviation ``s`` samples — which is
    what :func:`tfpws_stack`'s ``scales = w / (2*pi*f/fs)`` mapping assumes.

    .. note::

        This used to be ``x = np.linspace(-10, 10, M)``, i.e. a fixed
        abscissa independent of *M* and *s*.  The sample spacing was then
        ``20/(M-1) ~ 2/s``, making the actual centre frequency scale as
        ``1/s**2`` and the envelope width as ``s**2/2``.  The scale set
        requested by :func:`tfpws_stack` therefore analysed a completely
        different frequency band from the one asked for.

    :param M: Number of samples (odd).
    :param s: Scale parameter in samples (controls dilation).
    :param w: Central angular frequency (default 5.0).
    :returns: Complex 1-D array of length *M*.
    """
    x = np.arange(M) - (M - 1) / 2.0        # sample offsets from the centre
    wav = np.exp(1j * w * x / s) * np.exp(-0.5 * (x / s) ** 2)
    return wav / (np.pi ** 0.25 * np.sqrt(s))


def tfpws_stack(data, fs, freqmin, freqmax, power=2, nscales=20):
    """Time-frequency phase-weighted stack (Schimmel & Gallart 2007).

    Computes the CWT of every input trace with a complex Morlet wavelet at
    *nscales* log-spaced scales, derives the instantaneous-phase coherence
    :math:`c(a, t)` across traces at each (scale *a*, lag *t*) point,
    averages over scales to produce a per-lag weight :math:`w(t)`, then
    returns the linear mean multiplied by :math:`w(t)^{\\textit{power}}`.

    The Morlet wavelet used here is the standard complex Morlet in the
    convolution (not the FFT-based :class:`_Morlet` used by the WCT
    pipeline):

    .. math::

        \\psi_{M,s}(t) =
            \\frac{1}{\\pi^{1/4} \\sqrt{s}}\\,
            e^{i\\omega_0 t/s}\\,
            e^{-t^2/(2s^2)}, \\qquad \\omega_0 = 5

    The scale–frequency relationship is :math:`f = \\omega_0 f_s / (2\\pi s)`,
    so scales are chosen as
    :math:`s_k = \\omega_0 / (2\\pi f_k / f_s)` for *nscales* frequencies
    :math:`f_k` log-spaced in [*freqmin*, *freqmax*].

    The weight is:

    .. math::

        w(t) = \\left[
            \\frac{1}{A} \\sum_{k=1}^{A}
            \\frac{1}{N} \\left|
                \\sum_{j=1}^{N} e^{i\\,\\arg\\mathcal{W}_j(s_k, t)}
            \\right|
        \\right]^{\\textit{power}}

    where :math:`\\mathcal{W}_j(s_k, t)` is the CWT coefficient of trace
    *j* at scale *k* and lag *t*, and *A* = *nscales*.

    .. note::

        This function is called by :func:`stack` when
        ``stack_method="tfpws"``.  Prefer that entry-point in production
        code; call this function directly only when you need to tune
        *nscales* or *power* without going through the full stack
        dispatcher.

    :param data: 2-D array of shape ``(N_traces, N_lags)``.  NaN rows
        must be removed before calling (handled by :func:`stack`).
    :param fs: Sampling rate of the CCF traces (Hz).
    :param freqmin: Lower frequency bound for the CWT scale range (Hz).
        Should match the parent filter's *freqmin*.
    :param freqmax: Upper frequency bound for the CWT scale range (Hz).
        Should match the parent filter's *freqmax*.
    :param power: Exponent applied to the coherence weight (default 2).
        Equivalent to the ``pws_power`` parameter in :func:`stack`.
    :param nscales: Number of log-spaced CWT scales (default 20).
    :returns: 1-D array of length ``N_lags``.

    .. footcite:p:`Schimmel2007`
    """
    from scipy.signal import fftconvolve

    N, T = data.shape
    W = 5.0
    freqs  = np.logspace(np.log10(freqmin), np.log10(freqmax), nscales)
    scales = W / (2.0 * np.pi * freqs / fs)

    phase_sum = np.zeros((nscales, T), dtype=complex)
    for trace in data:
        for k, s in enumerate(scales):
            M = max(int(10 * s), 3) | 1
            wav = _morlet_wavelet(M, s, w=W)
            coeffs = fftconvolve(trace, wav[::-1].conj(), mode="same")
            phase_sum[k] += np.exp(1j * np.angle(coeffs))

    coh    = np.abs(phase_sum) / N
    weight = coh.mean(axis=0) ** power
    return data.mean(axis=0) * weight

def stack(data, stack_method="linear", pws_timegate=10.0, pws_power=2,
          goal_sampling_rate=20.0, freqmin=1.0, freqmax=10.0,
          tfpws_nscales=20):
    """Stack an array of CCF traces into a single representative trace.

    Three methods are available, selected via *stack_method*:

    **Linear stack** (``"linear"``)
        The arithmetic mean across all input traces:

        .. math::

            s(t) = \\frac{1}{N} \\sum_{j=1}^{N} d_j(t)

        Incoherent noise cancels as :math:`1/\\sqrt{N}`.  Fastest and most
        transparent, but offers no protection against high-amplitude
        transients (earthquakes, instrumental glitches) that survive
        pre-processing.

    **Phase-weighted stack** (``"pws"``)
        Introduced by Schimmel & Paulssen (1997) :footcite:p:`Schimmel1997`.
        Each sample is weighted by the instantaneous phase coherence
        :math:`c(t)` of the analytic signal across all traces:

        .. math::

            c(t) = \\frac{1}{N} \\left|
                \\sum_{j=1}^{N} e^{i\\,\\phi_j(t)}
            \\right|, \\qquad c(t) \\in [0, 1]

        where :math:`\\phi_j(t) = \\arg\\bigl(d_j(t) + i\\,\\mathcal{H}\\{d_j\\}(t)\\bigr)`
        is the instantaneous phase of trace *j* obtained via the Hilbert
        transform :math:`\\mathcal{H}`.  The coherence is smoothed with a
        boxcar window of *pws_timegate* seconds before being raised to the
        power *v* = *pws_power*:

        .. math::

            s(t) = \\frac{1}{N} \\sum_{j=1}^{N} d_j(t) \\cdot c(t)^v

        High-amplitude incoherent transients have random phases across
        traces, so :math:`c(t) \\approx 0` at those times; coherent
        arrivals have :math:`c(t) \\approx 1` and are preserved.

    **Time-frequency phase-weighted stack** (``"tfpws"``)
        The TF extension of PWS by Schimmel & Gallart (2007)
        :footcite:p:`Schimmel2007`.  Phase coherence is computed in the
        time-frequency domain via a continuous wavelet transform (CWT)
        with a complex Morlet wavelet, giving a coherence map
        :math:`c(a, t)` that is both scale- and time-dependent.
        Averaging over the *nscales* log-spaced scales spanning
        [*freqmin*, *freqmax*] Hz yields a single per-lag weight:

        .. math::

            W_j(a, t) = \\mathcal{W}\\{d_j\\}(a, t)

        .. math::

            c(a, t) = \\frac{1}{N} \\left|
                \\sum_{j=1}^{N} e^{i\\,\\arg W_j(a,t)}
            \\right|

        .. math::

            w(t) = \\left[
                \\frac{1}{A} \\sum_{a} c(a, t)
            \\right]^v, \\qquad
            s(t) = \\frac{1}{N} \\sum_{j=1}^{N} d_j(t) \\cdot w(t)

        where :math:`A` is the number of scales.  Because coherence is
        evaluated independently at each scale, tf-PWS is more sensitive
        to narrow-band coherent arrivals than time-domain PWS.  It is
        particularly effective for noise autocorrelations where body-wave
        reflections occupy a limited frequency band
        (Romero & Schimmel 2018 :footcite:p:`Romero2018`).

        .. note::

            Memory scales as :math:`O(N \\times A \\times T)` complex128.
            For long CCFs or large archives consider reducing *nscales*
            (default 20) or chunking pairs outside this function.

    

    :type data: :class:`numpy.ndarray`
    :param data: 2-D array of shape ``(N_traces, N_lags)``, each row one CCF.
    :type stack_method: str
    :param stack_method: ``"linear"``, ``"pws"``, or ``"tfpws"``.
    :type pws_timegate: float
    :param pws_timegate: Boxcar smoothing window for the PWS coherence
        estimate, in seconds (``"pws"`` only).  Default 10 s.
    :type pws_power: float
    :param pws_power: Exponent *v* applied to the coherence weight.
        Larger values increase selectivity.  Shared by ``"pws"`` and
        ``"tfpws"``.  Default 2.
    :type goal_sampling_rate: float
    :param goal_sampling_rate: Sampling rate of the CCF array (Hz).
    :type freqmin: float
    :param freqmin: Lower frequency bound (Hz) for the CWT scale range
        — ``"tfpws"`` only.  Should match the parent filter's *freqmin*.
    :type freqmax: float
    :param freqmax: Upper frequency bound (Hz) for the CWT scale range
        — ``"tfpws"`` only.  Should match the parent filter's *freqmax*.
    :type tfpws_nscales: int
    :param tfpws_nscales: Number of log-spaced CWT scales between
        *freqmin* and *freqmax* — ``"tfpws"`` only.  Default 20.
    :rtype: :class:`numpy.ndarray`
    :return: 1-D stacked CCF of length ``N_lags``, or ``[]`` if no valid
        (non-NaN) traces are present.
    """
    if len(data) == 0:
        logging.debug("No data to stack.")
        return []
    data = data[~np.isnan(data).any(axis=1)]
    # NOTE: a correlation-coefficient "sanitize" block used to sit here behind
    # a hardcoded ``sanitize = False``.  It was never reachable and its own
    # TODO flagged that it would wrongly discard stable pairs and
    # autocorrelations, so it has been removed rather than wired in.

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

        timegate_samples = max(1, int(pws_timegate * goal_sampling_rate))
        coh = np.convolve(ss.windows.boxcar(timegate_samples) /
                          timegate_samples, coh, 'same')
        coh = np.power(coh, pws_power)
        for c in data:
            corr += c * coh
        corr /= data.shape[0]

    elif stack_method == "tfpws":
        corr = tfpws_stack(
            data, fs=goal_sampling_rate,
            freqmin=freqmin, freqmax=freqmax,
            power=pws_power, nscales=tfpws_nscales,
        )

    else:
        logging.warning(f"Unknown stack_method {stack_method!r}; falling back to linear.")
        corr = data.mean(axis=0)

    return corr

# ── Wiener filter helpers (moved from msnoise/wiener.py) ──────────────────

def find_segments(data, gap_threshold):
    """Identify continuous non-NaN segments in an xarray DataArray.

    :param data: 2-D xarray DataArray (times × lags).
    :param gap_threshold: Maximum index gap before treating as a new segment.
        A run of ``> gap_threshold`` consecutive all-NaN rows starts a new
        segment; shorter runs are bridged.
    :returns: List of lists of row indices forming each continuous segment.
    """
    current_segment = []
    continuous_segments = []
    prev_idx = None

    for i in range(data.shape[0]):
        if not data[i, :].isnull().all():
            # NOTE: an ``else: prev_idx = None`` branch used to reset the
            # tracker on every all-NaN row.  Because ``i`` is a contiguous
            # range, ``i - prev_idx`` was then always 1 on the next valid row
            # and the gap test could never fire — find_segments always
            # returned a single segment and wiener_filt smeared straight
            # across data gaps, which is exactly what it exists to prevent.
            if prev_idx is not None and (i - prev_idx > gap_threshold):
                continuous_segments.append(current_segment)
                current_segment = []
            current_segment.append(i)
            prev_idx = i

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


def rolling_mean_nan(arr: "np.ndarray", window: int, min_periods: int = 1) -> "np.ndarray":
    """O(N) NaN-aware rolling mean along axis 0 using cumulative sums.

    Drop-in replacement for xarray's ``rolling(times=W).mean()`` on a 2-D
    numpy array.  Complexity is O(N × F) regardless of *window* size —
    unlike xarray's default implementation which is O(N × W × F) when NaN
    handling forces a fallback from the cumsum fast-path.

    :param arr:         2-D array of shape ``(T, F)`` (time × taxis).
    :param window:      Rolling window size (number of time steps).
    :param min_periods: Minimum number of non-NaN values required to produce
                        a result; positions with fewer valid values yield NaN.
    :returns:           Array of same shape and dtype as *arr*.
    """
    import numpy as np

    T, F = arr.shape

    # ── validity mask and value fill ──────────────────────────────────────
    valid  = ~np.isnan(arr)                          # bool  (T, F)
    filled = np.where(valid, arr, np.float32(0.0))   # float32 — no promotion

    # ── padded cumulative sums ─────────────────────────────────────────────
    # float32 cumsum is accurate for CCF values in [-1, 1]: worst-case
    # accumulated rounding error over 30 000 steps < 0.004, negligible.
    cv = np.empty((T + 1, F), dtype=np.float32)
    cv[0] = 0.0
    np.cumsum(filled, axis=0, out=cv[1:])
    del filled                                        # free early

    # int32 is sufficient: max count = window, well within int32.
    cc = np.empty((T + 1, F), dtype=np.int32)
    cc[0] = 0
    np.cumsum(valid, axis=0, out=cc[1:])
    del valid                                         # free early

    # ── window sums and counts ─────────────────────────────────────────────
    i1 = np.arange(1, T + 1)
    i0 = np.maximum(0, i1 - window)

    win_val = cv[i1] - cv[i0]                        # float32 (T, F)
    del cv
    win_cnt = cc[i1] - cc[i0]                        # int32  (T, F)
    del cc

    # ── masked mean ────────────────────────────────────────────────────────
    enough   = win_cnt >= min_periods
    safe_cnt = np.where(win_cnt > 0, win_cnt, np.int32(1))
    del win_cnt
    result   = np.where(enough, win_val / safe_cnt, np.float32(np.nan))
    del win_val, enough, safe_cnt
    return result.astype(arr.dtype)


def rolling_stack(arr: "np.ndarray", window: int,
                  stack_method: str = "linear",
                  min_periods: int = 1,
                  **stack_kwargs) -> "np.ndarray":
    """Rolling stack along axis 0 using the requested stacking method.

    Dispatches to the appropriate implementation:

    * ``"linear"``:  :func:`rolling_mean_nan` — O(N), constant in window size.
    * ``"pws"`` / ``"tfpws"``:  sliding-window loop calling :func:`stack`
      — inherently O(N × W) because phase coherence is non-linear, but
      avoids xarray overhead and skips all-NaN windows early.

    :param arr:          2-D float32 array ``(T, F)`` — time × taxis.
    :param window:       Rolling window size in time steps.
    :param stack_method: ``"linear"``, ``"pws"``, or ``"tfpws"``.
    :param min_periods:  Minimum valid (non-NaN) traces required; positions
                         with fewer yield NaN.
    :param stack_kwargs: Extra keyword arguments forwarded to :func:`stack`
                         (``pws_timegate``, ``pws_power``, ``freqmin``, …).
    :returns:            Array of same shape and dtype as *arr*.
    """
    import numpy as np

    if stack_method == "linear":
        return rolling_mean_nan(arr, window, min_periods)

    T, F = arr.shape
    result = np.full_like(arr, np.nan)

    for i in range(T):
        start = max(0, i - window + 1)
        window_data = arr[start:i + 1]
        # Drop whole-trace NaN rows (days with no data)
        valid_rows = ~np.all(np.isnan(window_data), axis=1)
        window_valid = window_data[valid_rows]
        if len(window_valid) < min_periods:
            continue
        result[i] = stack(window_valid, stack_method=stack_method, **stack_kwargs)

    return result
