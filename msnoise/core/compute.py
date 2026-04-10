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
    """
    if not index:
        return {}

    maxlag  = int(np.round(maxlag))
    Nt      = data.shape[1]
    ids     = [item[0] for item in index]
    sta1s   = np.array([item[1] for item in index], dtype=int)
    sta2s   = np.array([item[2] for item in index], dtype=int)
    n_pairs = len(ids)

    # Lag-window index (same for all pairs)
    if maxlag != Nt:
        tcorr = np.arange(-(Nt - 1), Nt)
        dN    = np.where(np.abs(tcorr) <= maxlag)[0]
    else:
        dN = None
    out_len = len(dN) if dN is not None else 2 * Nt
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
        folded_all[start:end] = folded[:, dN] if dN is not None else folded

    return {
        ids[k]: folded_all[k]
        for k in range(n_pairs)
        if folded_all.shape[1] >= min_len
    }

def pcc_xcorr(data, maxlag, energy, index, plot=False, nfft=None,
              normalized=False):
    """

    :param data:
    :param maxlag:
    :param energy:
    :param index:
    :param plot:
    :param nfft:
    :param normalized:
    :return:
    """
    from phasecorr.phasecorr import xcorr
    corr = {}
    ml = int(maxlag)
    for id, sta1, sta2 in index:
        corr[id] = xcorr(data[sta1], data[sta2],
                         lags=range(-ml, ml + 1),
                         parallel=True)
    return corr

def whiten(data, Nfft, delta, freqmin, freqmax, plot=False, returntime=False):
    """This function takes 1-dimensional *data* timeseries array,
    goes to frequency domain using fft, whitens the amplitude of the spectrum
    in frequency domain between *freqmin* and *freqmax*
    and returns the whitened fft.

    :type data: :class:`numpy.ndarray`
    :param data: Contains the 1D time series to whiten
    :type Nfft: int
    :param Nfft: The number of points to compute the FFT
    :type delta: float
    :param delta: The sampling frequency of the `data`
    :type freqmin: float
    :param freqmin: The lower frequency bound
    :type freqmax: float
    :param freqmax: The upper frequency bound
    :type plot: bool
    :param plot: Whether to show a raw plot of the action (default: False)

    :rtype: :class:`numpy.ndarray`
    :returns: The FFT of the input trace, whitened between the frequency bounds
    """

    # TODO: docsting
    if plot:
        import matplotlib.pyplot as plt
        plt.subplot(411)
        plt.plot(np.arange(len(data)) * delta, data)
        plt.xlim(0, len(data) * delta)
        plt.title('Input trace')

    Napod = 100
    Nfft = int(Nfft)
    freqVec = sf.fftfreq(Nfft, d=delta)[:Nfft // 2]
    J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
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
    FFTRawSign[high:Nfft + 1] *= 0

    # Hermitian symmetry (because the input is real)
    FFTRawSign[-(Nfft // 2) + 1:] = FFTRawSign[1:(Nfft // 2)].conjugate()[::-1]

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
    """This function takes 1-dimensional *data* timeseries array,
    goes to frequency domain using fft, whitens the amplitude of the spectrum
    in frequency domain between *freqmin* and *freqmax*
    and returns the whitened fft.

    :type data: :class:`numpy.ndarray`
    :param data: Contains the 1D time series to whiten
    :type Nfft: int
    :param Nfft: The number of points to compute the FFT
    :type delta: float
    :param delta: The sampling frequency of the `data`
    :type freqmin: float
    :param freqmin: The lower frequency bound
    :type freqmax: float
    :param freqmax: The upper frequency bound
    :type plot: bool
    :param plot: Whether to show a raw plot of the action (default: False)

    :rtype: :class:`numpy.ndarray`
    :returns: The FFT of the input trace, whitened between the frequency bounds
    """
    # TODO: docsting
    taper = np.ones(Nfft // 2 + 1)
    taper[0:low] *= 0
    taper[low:porte1] *= np.cos(np.linspace(np.pi / 2., 0, porte1 - low)) ** 2
    taper[porte2:high] *= np.cos(
        np.linspace(0., np.pi / 2., high - porte2)) ** 2
    taper[high:] *= 0
    taper *= taper

    hann = scipy.signal.windows.hann(porte2 - porte1 + 1)  # / float(porte2-porte1)

    for i in range(fft.shape[0]):
        if whiten_type == "PSD":
            fft[i][:Nfft // 2 + 1] /= psds[i]
            fft[i][:Nfft // 2 + 1] *= taper
            tmp = fft[i, porte1:porte2]
            imin = scoreatpercentile(tmp, 5)
            imax = scoreatpercentile(tmp, 95)
            not_outliers = np.where((tmp >= imin) & (tmp <= imax))[0]
            rms = tmp[not_outliers].std() * 1.0
            np.clip(fft[i, porte1:porte2], -rms, rms,
                    fft[i, porte1:porte2])  # inplace
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

        # Hermitian symmetry (because the input is real)
        fft[i, -(Nfft // 2) + 1:] = np.conjugate(fft[i, 1:(Nfft // 2)])[::-1]

def smooth(x, window='boxcar', half_win=3):
    """ some window smoothing """
    # TODO: docsting
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
    """The `current` time series is compared to the `reference`.
Both time series are sliced in several overlapping windows.
Each slice is mean-adjusted and cosine-tapered (85% taper) before being Fourier-
transformed to the frequency domain.
:math:`F_{cur}(\\nu)` and :math:`F_{ref}(\\nu)` are the first halves of the
Hermitian symmetric Fourier-transformed segments. The cross-spectrum
:math:`X(\\nu)` is defined as
:math:`X(\\nu) = F_{ref}(\\nu) F_{cur}^*(\\nu)`

in which :math:`{}^*` denotes the complex conjugation.
:math:`X(\\nu)` is then smoothed by convolution with a Hanning window.
The similarity of the two time-series is assessed using the cross-coherency
between energy densities in the frequency domain:

:math:`C(\\nu) = \\frac{|\\overline{X(\\nu))}|}{\\sqrt{|\\overline{F_{ref}(\\nu)|^2} |\\overline{F_{cur}(\\nu)|^2}}}`

in which the over-line here represents the smoothing of the energy spectra for
:math:`F_{ref}` and :math:`F_{cur}` and of the spectrum of :math:`X`. The mean
coherence for the segment is defined as the mean of :math:`C(\\nu)` in the
frequency range of interest. The time-delay between the two cross correlations
is found in the unwrapped phase, :math:`\\phi(\\nu)`, of the cross spectrum and is
linearly proportional to frequency:

:math:`\\phi_j = m. \\nu_j, m = 2 \\pi \\delta t`

The time shift for each window between two signals is the slope :math:`m` of a
weighted linear regression of the samples within the frequency band of interest.
The weights are those introduced by [Clarke2011]_,
which incorporate both the cross-spectral amplitude and cross-coherence, unlike
[Poupinet1984]_. The errors are estimated using the weights (thus the
coherence) and the squared misfit to the modelled slope:

:math:`e_m = \\sqrt{\\sum_j{(\\frac{w_j \\nu_j}{\\sum_i{w_i \\nu_i^2}})^2}\\sigma_{\\phi}^2}`

where :math:`w` are weights, :math:`\\nu` are cross-coherences and
:math:`\\sigma_{\\phi}^2` is the squared misfit of the data to the modelled slope
and is calculated as :math:`\\sigma_{\\phi}^2 = \\frac{\\sum_j(\\phi_j - m \\nu_j)^2}{N-1}`

The output of this process is a table containing, for each moving window: the
central time lag, the measured delay, its error and the mean coherence of the
segment.

.. warning::

    The time series will not be filtered before computing the cross-spectrum!
    They should be band-pass filtered around the `freqmin`-`freqmax` band of
    interest beforehand.

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
    from .core.signal import nextpow2
    padd = int(2 ** (nextpow2(window_length_samples) + 2))
    # padd = next_fast_len(window_length_samples)
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

        # Find the values the frequency range of interest
        freq_vec = sf.fftfreq(len(X) * 2, 1. / df)[:padd // 2]
        index_range = np.argwhere(np.logical_and(freq_vec >= freqmin,
                                                 freq_vec <= freqmax))

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

        # Phase:
        phi = np.angle(X)
        phi[0] = 0.
        phi = np.unwrap(phi)
        phi = phi[index_range]

        # Calculate the slope with a weighted least square linear regression
        # forced through the origin
        # weights for the WLS must be the variance !
        m, em = linear_regression(v.flatten(), phi.flatten(), w.flatten())

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
        del freq_vec
        del index_range
        del w, v, e, s2x2, sx2, m, em

    if maxind > len(current) + step * df:
        logging.warning("The last window was too small, but was computed")

    return np.array([time_axis, delta_t, delta_err, delta_mcoh]).T

