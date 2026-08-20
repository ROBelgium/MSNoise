"""MSNoise core computation functions — cross-correlation, whitening, MWCS.

Moved from ``msnoise/move2obspy.py``.
"""

__all__ = [
    "mwcs",
    "myCorr2",
    "pcc_xcorr",
    "smooth",
    "whiten",
    "whiten2",
    "compute_wct_dtt_batch",
    "resolve_wct_lag_min",
    "build_wct_dtt_dataset",
]

import logging

import numpy as np
import scipy
import scipy.fft as sf

import scipy.optimize
import scipy.signal
from scipy.stats import scoreatpercentile
from obspy.signal.invsim import cosine_taper
from obspy.signal.regression import linear_regression
from .signal import getCoherence


logger = logging.getLogger(__name__)

_MYCORR2_CHUNK = 64  # pairs per tile — tuned for L2 cache locality


def myCorr2(data, maxlag, energy, index, plot=False, nfft=None,
            normalized=False):
    """Compute cross-correlations for all requested station pairs.

    Tiled-batch implementation: processes pairs in chunks of ``_MYCORR2_CHUNK``
    to amortise Python loop overhead while keeping peak memory bounded.
    Faster than the original per-pair loop for N > ~50 stations (~2-3x for
    N=200-500); falls back gracefully to near-loop speed for small N.

    :type data: :class:`numpy.ndarray`
    :param data: 2-D array ``(n_stations, nfft)`` containing the FFT of each
        pre-whitened time series.
    :type maxlag: int
    :param maxlag: Output half-length in samples; CCF returned over
        ``[-maxlag : maxlag]`` (length ``2*maxlag + 1``).
    :type energy: :class:`numpy.ndarray`
    :param energy: Per-station RMS energy ``(n_stations,)`` used for POW
        normalisation.
    :type index: list
    :param index: List of ``(ccf_id, sta1_idx, sta2_idx)`` tuples.
    :type normalized: str or bool
    :param normalized: ``"POW"``, ``"MAX"``, ``"ABSMAX"``, or falsy for none.
    :rtype: dict
    :returns: ``{ccf_id: ccf_array}`` for every pair in *index*.

    .. note::

        **Amplitude scale.** ``sf.ifft`` already applies its own ``1/nfft``,
        and the result is divided by ``Nt == nfft`` a second time, so the CCF
        amplitude scales as ``1/nfft**2``.  This is the historical MSNoise
        convention and is kept for backward compatibility, but it means the
        absolute CCF amplitude depends on the FFT length.  Since ``nfft =
        next_fast_len(npts)`` is recomputed per window in
        :mod:`~msnoise.s03_compute_no_rotation`, two windows of a day with
        different ``npts`` enter the daily stack at slightly different
        amplitudes.  This is irrelevant for every dv/v method (MWCS, WCT and
        stretching are all scale-invariant) and for any
        ``cc_normalisation`` other than ``NO``, but it does affect raw CCF
        amplitudes.
    """
    if not index:
        return {}

    maxlag  = int(np.round(maxlag))
    Nt      = data.shape[1]
    ids     = [item[0] for item in index]
    sta1s   = np.array([item[1] for item in index], dtype=int)
    sta2s   = np.array([item[2] for item in index], dtype=int)
    n_pairs = len(ids)

    # Lag-window index (same for all pairs).  ``folded`` below is
    # ``[raw[1:Nt], raw[0:Nt]]`` → 2*Nt - 1 samples spanning lags
    # -(Nt-1) … +(Nt-1), so maxlag cannot exceed Nt-1.  The old code had a
    # separate ``maxlag == Nt`` branch that sized the output buffer at 2*Nt
    # and therefore raised a broadcast error; clamping instead keeps a single
    # code path and always returns a usable CCF.
    if maxlag > Nt - 1:
        logger.warning(
            "myCorr2: maxlag=%i exceeds the %i available lags for nfft=%i; "
            "clamping to %i", maxlag, Nt - 1, Nt, Nt - 1)
        maxlag = Nt - 1
    tcorr   = np.arange(-(Nt - 1), Nt)
    dN      = np.where(np.abs(tcorr) <= maxlag)[0]
    out_len = len(dN)
    min_len = 2 * maxlag + 1

    folded_all = np.empty((n_pairs, out_len), dtype=np.float64)

    # ── Tiled batch: chunk pairs to stay cache-friendly ─────────────────
    for start in range(0, n_pairs, _MYCORR2_CHUNK):
        end = min(start + _MYCORR2_CHUNK, n_pairs)
        s1  = sta1s[start:end]
        s2  = sta2s[start:end]

        # Cross-spectrum for this chunk
        cross  = np.conj(data[s1]) * data[s2]                # (chunk, nfft)
        raw    = np.real(sf.ifft(cross, n=nfft, axis=1)) / Nt
        folded = np.concatenate(
            [raw[:, -(Nt - 1):], raw[:, :Nt + 1]], axis=1
        )                                                     # (chunk, 2*Nt)

        # Normalisation applied before lag trim (preserves legacy behaviour:
        # the max is computed over the full folded CCF, not just the lag window)
        if normalized == "POW":
            norms = energy[s1] * energy[s2]
            norms = np.where(norms != 0, norms, 1.0)
            folded /= norms[:, None]
        elif normalized == "MAX":
            mx = folded.max(axis=1, keepdims=True)
            mx[mx == 0] = 1.0
            folded /= mx
        elif normalized == "ABSMAX":
            mx = np.abs(folded).max(axis=1, keepdims=True)
            mx[mx == 0] = 1.0
            folded /= mx

        # Lag trim
        folded_all[start:end] = folded[:, dN]

    return {
        ids[k]: folded_all[k]
        for k in range(n_pairs)
        if folded_all.shape[1] >= min_len
    }

def _analytic_phase(x: np.ndarray, eps_rel: float = 1e-6) -> np.ndarray:
    """Return the amplitude-normalised analytic signal (phase signal) of *x*.

    This is a direct Python translation of ``AnalyticSignal`` + ``AmpNorm``
    from FastPCC (:footcite:t:`Ventosa2019` / :footcite:t:`Ventosa2023`):

    1. FFT of real *x*
    2. Zero negative-frequency bins, double positive-frequency bins → analytic
    3. IFFT → complex envelope  ``X_a[n] = x[n] + i·H{x[n]}``
    4. Divide each sample by its magnitude (+ eps for numerical safety)
       so every sample lies on the complex unit circle.

    :param x: 1-D real time series.
    :param eps_rel: Stability floor = ``eps_rel * max|X_a|``. Matches
        FastPCC's ``AmpNormf`` (which uses ``1e-6 * sqrt(max_power)``).
    :returns: Complex array of the same length as *x*, ``|y[n]| ≤ 1``.
    """
    xa = scipy.signal.hilbert(x)          # scipy Hilbert = FFT-based analytic
    amp = np.abs(xa)
    eps = eps_rel * max(amp.max(), 1e-30)  # guard: all-zero input
    return xa / (amp + eps)


def _analytic_phase_batch(X: np.ndarray, eps_rel: float = 1e-6) -> np.ndarray:
    """Vectorised :func:`_analytic_phase` over rows of *X* (n_traces, N).

    :param X: 2-D real array ``(n_traces, N)``.
    :returns: Complex array ``(n_traces, N)`` of phase signals.
    """
    Xa = scipy.signal.hilbert(X, axis=1)
    amp = np.abs(Xa)
    row_max = amp.max(axis=1, keepdims=True)
    eps = eps_rel * np.where(row_max > 0, row_max, 1e-30)
    return Xa / (amp + eps)


def pcc_xcorr(data, maxlag, energy, index, plot=False, nfft=None,
              normalized=False):
    """Phase Cross-Correlation v=2 (PCC2) — pure NumPy/SciPy implementation.

    Replaces the former dependency on the unmaintained ``phasecorr`` package
    with a self-contained translation of the FFT-accelerated ``pcc2_set``
    routine from FastPCC (:footcite:t:`Ventosa2019`; :footcite:t:`Ventosa2023`).

    **Algorithm** (matches FastPCC ``pcc2_set``):

    1. Compute the amplitude-normalised analytic signal (phase signal)
       ``φ[n] = X_a[n] / |X_a[n]|`` for each trace — amplitude information
       is discarded entirely, so the result is insensitive to amplitude
       transients (earthquakes, glitches) without explicit temporal
       normalisation.
    2. Zero-pad to ``Nz = next_fast_len(N + maxlag)`` to avoid circular
       wrap-around (linear cross-correlation).
    3. Compute ``PCC2(lag) = IFFT(conj(FFT(φ1)) · FFT(φ2)) / (Nz · N)``
       for every pair in *index* — O(N log N), same cost as GNCC.

    :type data: :class:`numpy.ndarray`
    :param data: 2-D **time-domain** array ``(n_stations, N)``; real-valued.
        Unlike :func:`myCorr2`, PCC2 requires the time-domain input because
        the Hilbert transform must be computed before any FFT.
    :type maxlag: int or float
    :param maxlag: Half-length of output CCF in **samples**.
    :param energy: Unused (kept for API compatibility with :func:`myCorr2`).
    :param index: List of ``(ccf_id, sta1_idx, sta2_idx)`` tuples.
    :param normalized: ``"MAX"`` or ``"ABSMAX"`` to normalise output; falsy
        for none. (``"POW"`` is meaningless for PCC2 since amplitudes are
        discarded; it is silently ignored.)
    :rtype: dict
    :returns: ``{ccf_id: ccf_array}`` of length ``2*maxlag + 1`` per pair.

    .. note::

        **Biased estimator.** The sum is divided by the fixed ``N`` rather
        than by the ``N - |lag|`` samples that actually overlap, so PCC2
        carries a triangular taper of ``(N - |lag|) / N`` across the lag
        window (≈3 % at ``maxlag = 0.03 N``).  This matches the fixed
        ``1/nfft`` used by :func:`myCorr2`, so CC and PCC remain directly
        comparable.

    .. rubric:: References

    .. footcite:p:`Ventosa2019,Ventosa2023`
    """
    if not index:
        return {}

    ml  = int(np.round(maxlag))
    N   = data.shape[1]
    Nz  = sf.next_fast_len(N + ml)          # zero-pad length (avoids aliasing)
    norm_factor = float(N)   # IFFT already divides by Nz; divide by N for unit peak

    # Compute all phase signals at once (vectorised over traces) ──────────
    phase = _analytic_phase_batch(data)      # (n_stations, N), complex

    # Pre-FFT all phase signals into the padded length ───────────────────
    PHASE = sf.fft(phase, n=Nz, axis=1)     # (n_stations, Nz), complex

    corr = {}
    for ccf_id, sta1, sta2 in index:
        # Cross-spectrum → IFFT → real part
        xcorr_full = np.real(sf.ifft(np.conj(PHASE[sta1]) * PHASE[sta2]))
        xcorr_full /= norm_factor

        # Unwrap lags: [0 .. N-1] is positive lags, [Nz-ml .. Nz-1] is negative
        # Ventosa's convention: for lag<0 read from tail, for lag>=0 from head
        pos = xcorr_full[:ml + 1]       # lags  0 … +ml
        neg = xcorr_full[Nz - ml:Nz]   # lags -ml … -1
        ccf = np.concatenate([neg, pos])  # full: -ml … 0 … +ml

        if normalized == "MAX":
            mx = ccf.max()
            if mx != 0:
                ccf /= mx
        elif normalized == "ABSMAX":
            mx = np.abs(ccf).max()
            if mx != 0:
                ccf /= mx

        corr[ccf_id] = ccf

    return corr

def whiten(data, Nfft, delta, freqmin, freqmax, plot=False, returntime=False):
    """Spectral whitening (1-bit / brutal mode) for a single real trace.

    Computes the FFT of *data*, normalises the amplitude to unity in the
    passband ``[freqmin, freqmax]``, and returns the result.  A 100-sample
    cosine taper transitions smoothly to zero outside the passband, ensuring
    no sharp spectral edges.  This is the "brutal" whitening described by
    :footcite:t:`Bensen2007` (their Section 2.2).

    .. math::

        \\tilde{X}[k] = \\exp\\!\\left(i\\,\\arg X[k]\\right)
        \\qquad \\text{for } \\nu_k \\in [f_\\text{low},\\, f_\\text{high}]

    Bins outside the passband are zeroed (with a cosine-tapered transition
    region); Hermitian symmetry is enforced so that the inverse FFT yields a
    real signal.

    :type data: :class:`numpy.ndarray`
    :param data: 1-D real time series.
    :type Nfft: int
    :param Nfft: FFT length (zero-padding if larger than ``len(data)``).
    :type delta: float
    :param delta: Sampling interval in seconds (``1 / sampling_rate``).
    :type freqmin: float
    :param freqmin: Lower passband frequency (Hz).
    :type freqmax: float
    :param freqmax: Upper passband frequency (Hz).
    :type plot: bool
    :param plot: Show a diagnostic plot of the whitening stages (default False).
    :type returntime: bool
    :param returntime: If True, return the whitened time-domain signal instead
        of the frequency-domain array.
    :rtype: :class:`numpy.ndarray`
    :returns: Whitened one-sided FFT array (complex) of length ``Nfft``, or the
        corresponding real time-domain signal if *returntime* is True.
    """
    if plot:
        import matplotlib.pyplot as plt
        plt.subplot(411)
        plt.plot(np.arange(len(data)) * delta, data)
        plt.xlim(0, len(data) * delta)
        plt.title('Input trace')

    Napod = 100
    Nfft = int(Nfft)
    # Positive-frequency bins are 0 … (Nfft-1)//2 — which is Nfft//2 - 1 for
    # even Nfft but Nfft//2 for odd Nfft (next_fast_len can return 9, 15, 25,
    # 27, 45, 75, 81, 125 …), hence (Nfft - 1) // 2 + 1 rather than Nfft // 2.
    freqVec = sf.fftfreq(Nfft, d=delta)[:(Nfft - 1) // 2 + 1]
    J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
    if J.size == 0:
        raise ValueError(
            f"whiten: no FFT bin falls inside the [{freqmin:g}, {freqmax:g}] Hz "
            f"passband (Nfft={Nfft}, delta={delta:g} s → bin width "
            f"{1.0 / (Nfft * delta):g} Hz, Nyquist {0.5 / delta:g} Hz)")
    low = J[0] - Napod
    if low <= 0:
        low = 1

    porte1 = J[0]
    porte2 = J[-1]
    high = J[-1] + Napod
    if high > Nfft / 2:
        high = int(Nfft // 2)

    FFTRawSign = sf.fft(data, Nfft)

    if plot:
        plt.subplot(412)
        axis = np.arange(len(FFTRawSign))
        plt.plot(axis[1:], np.abs(FFTRawSign[1:]))
        plt.xlim(0, max(axis))
        plt.title('FFTRawSign')

    # Left tapering:
    FFTRawSign[0:low] *= 0
    FFTRawSign[low:porte1] = np.cos(
        np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[low:porte1]))
    # Pass band:
    FFTRawSign[porte1:porte2] = np.exp(1j * np.angle(FFTRawSign[porte1:porte2]))
    # Right tapering:
    FFTRawSign[porte2:high] = np.cos(
        np.linspace(0., np.pi / 2., high - porte2)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[porte2:high]))
    FFTRawSign[high:] *= 0

    # Hermitian symmetry (because the input is real): X[Nfft-k] = conj(X[k])
    # for k = 1 … (Nfft-1)//2.  Writing it in terms of (Nfft-1)//2 rather than
    # Nfft//2 keeps it correct for odd Nfft, where the old slice left the
    # first negative-frequency bin untouched.
    nneg = (Nfft - 1) // 2
    FFTRawSign[Nfft - nneg:] = FFTRawSign[1:nneg + 1].conjugate()[::-1]

    if plot:
        plt.subplot(413)
        axis = np.arange(len(FFTRawSign))
        plt.axvline(low, c='g')
        plt.axvline(porte1, c='g')
        plt.axvline(porte2, c='r')
        plt.axvline(high, c='r')

        plt.axvline(Nfft - high, c='r')
        plt.axvline(Nfft - porte2, c='r')
        plt.axvline(Nfft - porte1, c='g')
        plt.axvline(Nfft - low, c='g')

        plt.plot(axis, np.abs(FFTRawSign))
        plt.xlim(0, max(axis))

        wdata = np.real(sf.ifft(FFTRawSign, Nfft))
        plt.subplot(414)
        plt.plot(np.arange(len(wdata)) * delta, wdata)
        plt.xlim(0, len(wdata) * delta)
        plt.show()
    if returntime:
        return np.real(sf.ifft(FFTRawSign, Nfft))[:len(data)]
    return FFTRawSign

def whiten2(fft, Nfft, low, high, porte1, porte2, psds, whiten_type):
    """Vectorised in-place spectral whitening for a batch of pre-computed FFTs.

    Operates on the one-sided positive-frequency half of each FFT row and
    enforces Hermitian symmetry afterward. Three modes are available via
    *whiten_type*, corresponding to the normalisation strategies described
    by :footcite:t:`Bensen2007`:

    - **brutal** (default) — one-bit normalisation: amplitude set to unity
      inside the passband with a cosine taper at both edges:

      .. math::

          \\tilde{X}[k] = \\exp\\!\\left(i\\,\\arg X[k]\\right)

    - **HANN** — one-bit normalisation weighted by a Hann window across the
      passband, smoothly tapering the spectral amplitude:

      .. math::

          \\tilde{X}[k] = \\frac{X[k]}{|X[k]|}\\cdot w_\\text{Hann}[k]

    - **PSD** — divide by a pre-computed smoothed PSD, then clip outlier bins
      at the 5th–95th percentile to suppress spectral spikes:

      .. math::

          \\tilde{X}[k] = \\operatorname{clip}\\!\\left(
              \\frac{X[k]}{S[k]},\\,-A,\\,A
          \\right)

      where :math:`S[k]` is the smoothed PSD and :math:`A` is the RMS of
      the non-outlier bins.

    :type fft: :class:`numpy.ndarray`
    :param fft: 2-D complex array ``(n_traces, Nfft)`` of pre-computed FFTs.
        **Modified in-place.**
    :type Nfft: int
    :param Nfft: FFT length (must match ``fft.shape[1]``).
    :type low: int
    :param low: Bin index where the left cosine taper begins (below passband).
    :type high: int
    :param high: Bin index where the right cosine taper ends (above passband).
    :type porte1: int
    :param porte1: First bin of the flat passband.
    :type porte2: int
    :param porte2: Last bin of the flat passband.
    :type psds: :class:`numpy.ndarray` or None
    :param psds: Smoothed PSD array ``(n_traces, Nfft//2+1)``; used only
        when *whiten_type* is ``"PSD"``.
    :type whiten_type: str
    :param whiten_type: One of ``"brutal"`` (default), ``"HANN"``, or
        ``"PSD"``.
    :returns: None (modifies *fft* in-place).
    """
    taper = np.ones(Nfft // 2 + 1)
    taper[0:low] *= 0
    taper[low:porte1] *= np.cos(np.linspace(np.pi / 2., 0, porte1 - low)) ** 2
    taper[porte2:high] *= np.cos(
        np.linspace(0., np.pi / 2., high - porte2)) ** 2
    taper[high:] *= 0
    # NOTE: a stray ``taper *= taper`` used to sit here, making the PSD-branch
    # taper cos^4 while whiten() and the brutal branch use cos^2.  Removed so
    # all three whitening shapes share the same spectral window.

    hann = scipy.signal.windows.hann(porte2 - porte1 + 1)  # / float(porte2-porte1)
    nneg = (Nfft - 1) // 2   # number of negative-frequency bins

    for i in range(fft.shape[0]):
        if whiten_type == "PSD":
            fft[i][:Nfft // 2 + 1] /= psds[i]
            fft[i][:Nfft // 2 + 1] *= taper
            tmp = fft[i, porte1:porte2]
            # Percentiles, the outlier mask and the clip must all be taken on
            # the MAGNITUDE.  Applied to the complex array (as they were),
            # numpy resolves them lexicographically on the real part, and
            # np.clip() replaces every out-of-range bin by the real scalar
            # +/-rms — i.e. it zeroes the imaginary part and destroys the phase
            # of exactly the bins the spike suppression is aimed at.
            amp = np.abs(tmp)
            imin = scoreatpercentile(amp, 5)
            imax = scoreatpercentile(amp, 95)
            not_outliers = np.where((amp >= imin) & (amp <= imax))[0]
            # std() of a (zero-mean) complex spectrum == RMS of its magnitudes,
            # which is the "A" of the docstring formula.
            rms = float(tmp[not_outliers].std()) if not_outliers.size else 0.0
            if rms > 0:
                scale = np.ones_like(amp)
                np.divide(rms, amp, out=scale, where=amp > rms)
                fft[i, porte1:porte2] = tmp * scale
            fft[i, 0:low] *= 0
            fft[i, high:] *= 0
        elif whiten_type == "HANN":
            np.divide(fft[i], np.abs(fft[i]), out=fft[i], where=fft[i]!=0)
            fft[i][:porte1] *= 0.0
            fft[i][porte1:porte2 + 1] *= hann
            fft[i][porte2 + 1:] *= 0.0
        else:
            # print("Doing the classic Brutal Whiten")
            # Left tapering:
            fft[i, 0:low] *= 0
            fft[i, low:porte1] = np.cos(
                np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * np.exp(
                1j * np.angle(fft[i, low:porte1]))
            # Pass band:
            fft[i, porte1:porte2] = np.exp(1j * np.angle(fft[i, porte1:porte2]))
            # Right tapering:
            fft[i, porte2:high] = np.cos(
                np.linspace(0., np.pi / 2., high - porte2)) ** 2 * np.exp(
                1j * np.angle(fft[i, porte2:high]))
            fft[i, high:] *= 0

        # Hermitian symmetry (because the input is real): X[Nfft-k] =
        # conj(X[k]) for k = 1 … (Nfft-1)//2 — correct for odd Nfft too.
        fft[i, Nfft - nneg:] = np.conjugate(fft[i, 1:nneg + 1])[::-1]

def smooth(x, window='boxcar', half_win=3):
    """Smooth a 1-D array with a symmetric window.

    Pads the signal by reflection at both ends (length ``2*half_win + 1``)
    to reduce boundary effects, then convolves with a normalised window.

    :type x: :class:`numpy.ndarray`
    :param x: 1-D input array.
    :type window: str
    :param window: ``"boxcar"`` (uniform, default) or ``"hanning"``.
    :type half_win: int
    :param half_win: Half-width of the window; full width = ``2*half_win+1``.
    :rtype: :class:`numpy.ndarray`
    :returns: Smoothed array, same length as *x*.
    """
    window_len = 2 * half_win + 1
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    if window == "boxcar":
        w = scipy.signal.windows.boxcar(window_len).astype('complex')
    else:
        w = scipy.signal.windows.hann(window_len).astype('complex')
    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[half_win:len(y) - half_win]

def mwcs(current, reference, freqmin, freqmax, df, tmin, window_length, step,
         smoothing_half_win=5):
    """Moving Window Cross-Spectrum (MWCS) time-delay measurement.

    Both time series are sliced into overlapping windows. Each window is
    mean-adjusted and cosine-tapered (85%) before being Fourier-transformed.

    The cross-spectrum between the reference and current windows is:

    .. math::

        X(\\nu) = F_\\text{ref}(\\nu)\\, F_\\text{cur}^*(\\nu)

    where :math:`{}^*` denotes complex conjugation. :math:`X(\\nu)` is smoothed
    by convolution with a Hanning window. Cross-coherency is then:

    .. math::

        C(\\nu) = \\frac{
            \\left| \\overline{X(\\nu)} \\right|
        }{
            \\sqrt{
                \\overline{\\left|F_\\text{ref}(\\nu)\\right|^2}\\;
                \\overline{\\left|F_\\text{cur}(\\nu)\\right|^2}
            }
        }

    where the over-bar denotes smoothing. The mean coherence per window is the
    mean of :math:`C(\\nu)` over the frequency band.

    The unwrapped phase of :math:`X(\\nu)` is linearly proportional to frequency:

    .. math::

        \\phi_j = 2\\pi\\,\\delta t\\,\\nu_j

    so the time delay :math:`\\delta t` is the slope of a weighted linear
    regression over the frequency band (:footcite:t:`Clarke2011`, extending
    :footcite:t:`Poupinet1984`). Weights incorporate both cross-spectral
    amplitude and coherence. The slope error is:

    .. math::

        e_m = \\sqrt{
            \\sum_j \\left(
                \\frac{w_j\\,\\nu_j}{\\sum_i w_i\\,\\nu_i^2}
            \\right)^2 \\sigma_\\phi^2
        }, \\qquad
        \\sigma_\\phi^2 = \\frac{\\sum_j (\\phi_j - m\\,\\nu_j)^2}{N-1}

    where :math:`w_j` are the per-sample weights and :math:`\\nu_j` are the
    cross-coherences.

    Returns one row per moving window: central lag time, :math:`\\delta t`,
    error, and mean coherence.

.. warning::

    The time series will not be filtered before computing the cross-spectrum!
    They should be band-pass filtered around the `freqmin`-`freqmax` band of
    interest beforehand.

.. warning::

    The phase is unwrapped **inside** ``[freqmin, freqmax]`` only. This is
    valid as long as :math:`|\\delta t| < 1 / (2\\,f_\\text{min})`, which holds
    by many orders of magnitude for dv/v monitoring. Before MSNoise 2.0 the
    unwrap started at DC, which let cycle slips in the incoherent sub-``freqmin``
    bins add a constant :math:`2\\pi k` to the whole band and produce isolated
    spurious ``dt`` values.

:type current: :class:`numpy.ndarray`
:param current: The "Current" timeseries
:type reference: :class:`numpy.ndarray`
:param reference: The "Reference" timeseries
:type freqmin: float
:param freqmin: The lower frequency bound to compute the dephasing (in Hz)
:type freqmax: float
:param freqmax: The higher frequency bound to compute the dephasing (in Hz)
:type df: float
:param df: The sampling rate of the input timeseries (in Hz)
:type tmin: float
:param tmin: The leftmost time lag (used to compute the "time lags array")
:type window_length: float
:param window_length: The moving window length (in seconds)
:type step: float
:param step: The step to jump for the moving window (in seconds)
:type smoothing_half_win: int
:param smoothing_half_win: If different from 0, defines the half length of
    the smoothing hanning window.

:rtype: :class:`numpy.ndarray`
:returns: [time_axis,delta_t,delta_err,delta_mcoh]. time_axis contains the
    central times of the windows. The three other columns contain dt, error and
    mean coherence for each window.

    .. footbibliography::
    """
    delta_t = []
    delta_err = []
    delta_mcoh = []
    time_axis = []

    window_length_samples = int(window_length * df)
    step_samples = int(step * df)
    # try:
    #     from sf.helper import next_fast_len
    # except ImportError:
    #     from obspy.signal.util import next_pow_2 as next_fast_len
    from .signal import nextpow2
    padd = int(2 ** (nextpow2(window_length_samples) + 2))
    # padd = next_fast_len(window_length_samples)

    # Frequency axis and measurement band are loop-invariant: build them once
    # and fail early with a readable message rather than deep inside the loop.
    freq_vec = sf.fftfreq(padd, 1. / df)[:padd // 2]
    index_range = np.argwhere(np.logical_and(freq_vec >= freqmin,
                                             freq_vec <= freqmax)).flatten()
    if index_range.size < 2:
        raise ValueError(
            f"mwcs: the [{freqmin:g}, {freqmax:g}] Hz band contains "
            f"{index_range.size} FFT bin(s) for a {window_length:g} s window "
            f"at {df:g} Hz (bin width {df / padd:g} Hz). Widen the band or "
            f"increase mwcs_wlen.")

    count = 0
    tp = cosine_taper(window_length_samples, 0.85)
    minind = 0
    maxind = window_length_samples
    while maxind <= len(current):
        cci = current[minind:(minind + window_length_samples)]
        cci = scipy.signal.detrend(cci, type='linear')
        cci *= tp

        cri = reference[minind:(minind + window_length_samples)]
        cri = scipy.signal.detrend(cri, type='linear')
        cri *= tp

        minind += step_samples
        maxind += step_samples

        fcur = sf.fft(cci, n=padd)[:padd // 2]
        fref = sf.fft(cri, n=padd)[:padd // 2]

        fcur2 = np.real(fcur) ** 2 + np.imag(fcur) ** 2
        fref2 = np.real(fref) ** 2 + np.imag(fref) ** 2

        # Calculate the cross-spectrum
        X = fref * (fcur.conj())
        if smoothing_half_win != 0:
            dcur = np.sqrt(smooth(fcur2, window='hanning',
                                  half_win=smoothing_half_win))
            dref = np.sqrt(smooth(fref2, window='hanning',
                                  half_win=smoothing_half_win))
            X = smooth(X, window='hanning',
                       half_win=smoothing_half_win)
        else:
            dcur = np.sqrt(fcur2)
            dref = np.sqrt(fref2)

        dcs = np.abs(X)

        # Get Coherence and its mean value
        coh = getCoherence(dcs, dref, dcur)
        mcoh = np.mean(coh[index_range])

        # Get Weights
        w = 1.0 / (1.0 / (coh[index_range] ** 2) - 1.0)
        w[coh[index_range] >= 0.99] = 1.0 / (1.0 / 0.9801 - 1.0)
        w = np.sqrt(w * np.sqrt(dcs[index_range]))
        w = np.real(w)

        # Frequency array:
        v = np.real(freq_vec[index_range]) * 2 * np.pi

        # Phase — unwrapped INSIDE the measurement band only.
        #
        # Unwrapping from DC (as this used to do) walks through the bins below
        # ``freqmin``, where a band-passed CCF carries only filter leakage and
        # the phase is essentially random.  np.unwrap then performs a random
        # walk there, and any cycle slip it picks up shifts EVERY in-band phase
        # by a constant 2*pi*k.  Because the regression below is forced through
        # the origin, that constant offset maps directly into a spurious
        # dt ≈ k / f — the classic isolated MWCS outlier.
        #
        # Anchoring the unwrap on the first in-band bin instead assumes only
        # |dt| < 1 / (2 * freqmin) (e.g. 0.5 s for freqmin = 1 Hz), which is
        # orders of magnitude above any dv/v measurement.  Within the band the
        # bin-to-bin phase step is 2*pi*dt*df/padd, so continuity there is
        # never in question.
        phi = np.unwrap(np.angle(X)[index_range])

        # Calculate the slope with a weighted least square linear regression
        # forced through the origin
        # weights for the WLS must be the variance !
        m, em = linear_regression(v, phi, w)

        delta_t.append(m)

        # print phi.shape, v.shape, w.shape
        e = np.sum((phi - m * v) ** 2) / (np.size(v) - 1)
        s2x2 = np.sum(v ** 2 * w ** 2)
        sx2 = np.sum(w * v ** 2)
        e = np.sqrt(e * s2x2 / sx2 ** 2)

        delta_err.append(e)
        delta_mcoh.append(np.real(mcoh))
        time_axis.append(tmin + window_length / 2. + count * (step_samples/df))
        count += 1

        del fcur, fref
        del X
        del w, v, e, s2x2, sx2, m, em

    # The loop only emits complete windows; report whatever tail was dropped.
    # (The previous check could never fire: the loop exits as soon as
    # maxind > len(current), so maxind <= len(current) + step_samples always.)
    if count:
        used = (count - 1) * step_samples + window_length_samples
    else:
        used = 0
    tail = len(current) - used
    if tail > 0:
        logger.debug(
            "mwcs: %i trailing sample(s) did not fill a complete %i-sample "
            "window and were not used", tail, window_length_samples)

    return np.array([time_axis, delta_t, delta_err, delta_mcoh]).T



# ── WCT-DTT helpers (shared between s08 fused mode and s09) ─────────────────


def resolve_wct_lag_min(dtt_params, dist: float) -> float:
    """Resolve the coda lag minimum for WCT dt/t computation.

    :param dtt_params: Params object with ``wavelet_dtt`` attributes.
    :param dist: Interstation distance in km (used for dynamic lag).
    :returns: Lag minimum in seconds.
    """
    lag_type = str(getattr(dtt_params.wavelet_dtt, "wct_lag", None) or "static")
    v = float(getattr(dtt_params.wavelet_dtt, "wct_v", 1.0) or 1.0)
    minlag = float(getattr(dtt_params.wavelet_dtt, "wct_minlag", 5.0))
    if lag_type == "dynamic" and v > 0:
        return dist / v
    return minlag


def compute_wct_dtt_batch(freqs, taxis, WXamp, Wcoh, WXdt, dtt_params, dist: float = 0.0):
    """Compute WCT dt/t and average coherence for one time-step.

    Shared implementation used by both the fused path in s08 (where WXamp/Wcoh/WXdt
    are freshly computed in memory) and the standalone s09 path (where they are
    loaded from a WCT NetCDF file).

    :param freqs: 1-D frequency array from the WCT (Hz).
    :param taxis: 1-D lag-time axis array (s).
    :param WXamp: 2-D cross-wavelet amplitude array ``(freqs, taxis)``.
    :param Wcoh: 2-D wavelet coherence array ``(freqs, taxis)``.
    :param WXdt: 2-D time-delay array ``(freqs, taxis)``.
    :param dtt_params: Merged params object containing ``wavelet_dtt.*`` attributes.
    :param dist: Interstation distance in km (for dynamic lag). Default 0.
    :returns: Tuple ``(dtt_row, err_row, coh_row, freqs_subset)`` where
        *dtt_row* and *err_row* are 1-D arrays over the DTT frequency subset
        and *coh_row* is the average coherence per frequency bin.
    """
    from .signal import compute_wct_dtt, get_wct_avgcoh

    fp = dtt_params.wavelet_dtt
    lag_min = resolve_wct_lag_min(dtt_params, dist)
    freqmin = fp.wct_dtt_freqmin
    freqmax = fp.wct_dtt_freqmax
    coda_cycles = int(fp.wct_codacycles)

    dtt_row, err_row, _ = compute_wct_dtt(
        freqs, taxis, WXamp, Wcoh, WXdt,
        lag_min=lag_min,
        coda_cycles=coda_cycles,
        mincoh=fp.wct_mincoh,
        maxdt=fp.wct_maxdt,
        min_nonzero=fp.wct_min_nonzero,
        freqmin=freqmin,
        freqmax=freqmax,
    )
    coh_row = get_wct_avgcoh(
        freqs, taxis, Wcoh,
        freqmin=freqmin,
        freqmax=freqmax,
        lag_min=lag_min,
        coda_cycles=coda_cycles,
    )
    mask = (freqs >= freqmin) & (freqs <= freqmax)
    freqs_subset = freqs[mask]
    return dtt_row, err_row, coh_row, freqs_subset


def build_wct_dtt_dataset(dates_list, dtt_rows, err_rows, coh_rows, freqs_subset):
    """Build a sorted xarray Dataset of WCT dt/t results.

    Shared between the fused save path in s08 and the save path in s09.

    :param dates_list: List of datetime64 timestamps.
    :param dtt_rows: List of 1-D dt/t arrays (one per date).
    :param err_rows: List of 1-D error arrays.
    :param coh_rows: List of 1-D average-coherence arrays.
    :param freqs_subset: 1-D frequency array for the DTT band.
    :returns: :class:`xarray.Dataset` with variables ``DTT``, ``ERR``, ``COH``
        and dims ``(times, frequency)``, sorted by ``times``.
    """
    import numpy as np
    import xarray as xr

    dates_arr = np.array(dates_list, dtype="datetime64[ns]")
    sort_idx  = np.argsort(dates_arr)
    return xr.Dataset({
        "DTT": xr.DataArray(
            np.array(dtt_rows)[sort_idx],
            dims=["times", "frequency"],
            coords={"times": dates_arr[sort_idx], "frequency": freqs_subset},
        ),
        "ERR": xr.DataArray(
            np.array(err_rows)[sort_idx],
            dims=["times", "frequency"],
            coords={"times": dates_arr[sort_idx], "frequency": freqs_subset},
        ),
        "COH": xr.DataArray(
            np.array(coh_rows)[sort_idx],
            dims=["times", "frequency"],
            coords={"times": dates_arr[sort_idx], "frequency": freqs_subset},
        ),
    })
