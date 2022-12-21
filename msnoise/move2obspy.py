import logging


import numpy as np
import scipy
if scipy.__version__ < "1.4.0":
    import scipy.fftpack as sf
    from scipy.fftpack.helper import next_fast_len
    import scipy.fftpack._fftpack as sff
else:
    import scipy.fft as sf
    from scipy.fft import next_fast_len

import scipy.optimize
import scipy.signal
from scipy.stats import scoreatpercentile
from obspy.signal.invsim import cosine_taper
from obspy.signal.regression import linear_regression



def myCorr(data, maxlag, plot=False, nfft=None):
    """This function takes ndimensional *data* array, computes the cross-correlation in the frequency domain
    and returns the cross-correlation function between [-*maxlag*:*maxlag*].

    :type data: :class:`numpy.ndarray`
    :param data: This array contains the fft of each timeseries to be cross-correlated.
    :type maxlag: int
    :param maxlag: This number defines the number of samples (N=2*maxlag + 1) of the CCF that will be returned.

    :rtype: :class:`numpy.ndarray`
    :returns: The cross-correlation function between [-maxlag:maxlag]
    """
    # TODO: docsting
    if nfft is None:
        s1 = np.array(data[0].shape)
        shape = s1 - 1
        # Speed up FFT by padding to optimal size for FFTPACK
        fshape = [next_fast_len(int(d)) for d in shape]
        nfft = fshape[0]

    normalized = True
    allCpl = False

    maxlag = np.round(maxlag)

    Nt = data.shape[1]

    # data = scipy.fftpack.fft(data,int(Nfft),axis=1)

    if plot:
        import matplotlib.pyplot as plt
        plt.subplot(211)
        plt.plot(np.arange(len(data[0])) * 0.05, np.abs(data[0]))
        plt.subplot(212)
        plt.plot(np.arange(len(data[1])) * 0.05, np.abs(data[1]))

    corr = np.conj(data[0]) * data[1]
    corr = np.real(sf.ifft(corr, nfft)) / Nt
    corr = np.concatenate((corr[-Nt + 1:], corr[:Nt + 1]))

    if plot:
        plt.figure()
        plt.plot(corr)

    if normalized:
        E = np.prod(np.real(np.sqrt(
            np.mean(sf.ifft(data, n=nfft, axis=1) ** 2, axis=1))))
        corr /= np.real(E)

    if maxlag != Nt:
        tcorr = np.arange(-Nt + 1, Nt)
        dN = np.where(np.abs(tcorr) <= maxlag)[0]
        corr = corr[dN]

    del data
    return corr


def myCorr2(data, maxlag, energy, index, plot=False, nfft=None, 
            normalized=False):
    """This function takes ndimensional *data* array, computes the
    cross-correlation in the frequency domain and returns the cross-correlation
    function between [-*maxlag*:*maxlag*].

    :type data: :class:`numpy.ndarray`
    :param data: This array contains the fft of each timeseries to be cross-correlated.
    :type maxlag: int
    :param maxlag: This number defines the number of samples (N=2*maxlag + 1) of the CCF that will be returned.

    :rtype: :class:`numpy.ndarray`
    :returns: The cross-correlation function between [-maxlag:maxlag]
    """
    # TODO: docsting

    maxlag = np.round(maxlag)
    Nt = data.shape[1]

    if maxlag != Nt:
        tcorr = np.arange(-Nt + 1, Nt)
        dN = np.where(np.abs(tcorr) <= maxlag)[0]
    corrs = {}
    for id, sta1, sta2 in index:
        corr = np.conj(data[sta1]) * data[sta2]
        corr = np.real(sf.ifft(corr, nfft)) / Nt
        corr = np.concatenate((corr[-Nt + 1:], corr[:Nt + 1]))

        if normalized:
            corr /= (energy[sta1] * energy[sta2])

        if maxlag != Nt:
            corr = corr[dN]
        if len(corr) < (2*maxlag)+1:
            continue
        corrs[id] = corr

    return corrs


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

def whiten(data, Nfft, delta, freqmin, freqmax, plot=False):
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
    taper = np.ones(Nfft//2+1)
    taper[0:low] *= 0
    taper[low:porte1] *= np.cos(np.linspace(np.pi / 2., 0, porte1 - low))**2
    taper[porte2:high] *= np.cos(np.linspace(0., np.pi / 2., high - porte2))**2
    taper[high:] *= 0
    taper *= taper
    for i in range(fft.shape[0]):
        if whiten_type == "PSD":
            fft[i][:Nfft//2+1] /= psds[i]
            fft[i][:Nfft//2+1] *= taper
            tmp = fft[i, porte1:porte2]
            imin = scoreatpercentile(tmp, 5)
            imax = scoreatpercentile(tmp, 95)
            not_outliers = np.where((tmp >= imin) & (tmp <= imax))[0]
            rms = tmp[not_outliers].std() * 1.0
            np.clip(fft[i, porte1:porte2], -rms, rms, fft[i, porte1:porte2])  # inplace
            fft[i, 0:low] *= 0
            fft[i, high:] *= 0
        else:
            # print("Doing the classic Brutal Whiten")
            # Left tapering:
            fft[i,0:low] *= 0
            fft[i,low:porte1] = np.cos(np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * np.exp(1j * np.angle(fft[i,low:porte1]))
            # Pass band:
            fft[i,porte1:porte2] = np.exp(1j * np.angle(fft[i,porte1:porte2]))
            # Right tapering:
            fft[i,porte2:high] = np.cos(np.linspace(0., np.pi / 2., high - porte2)) ** 2 * np.exp(1j * np.angle(fft[i,porte2:high]))
            fft[i,high:] *= 0

        # Hermitian symmetry (because the input is real)
        fft[i,-(Nfft // 2) + 1:] = np.conjugate(fft[i,1:(Nfft // 2)])[::-1]



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


def getCoherence(dcs, ds1, ds2):
    # TODO: docsting
    n = len(dcs)
    coh = np.zeros(n).astype('complex')
    valids = np.argwhere(np.logical_and(np.abs(ds1) > 0, np.abs(ds2 > 0)))
    coh[valids] = dcs[valids] / (ds1[valids] * ds2[valids])
    coh[coh > (1.0 + 0j)] = 1.0 + 0j
    return coh


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

:math:`C(\\nu) = \\frac{|\overline{X(\\nu))}|}{\sqrt{|\overline{F_{ref}(\\nu)|^2} |\overline{F_{cur}(\\nu)|^2}}}`


in which the over-line here represents the smoothing of the energy spectra for
:math:`F_{ref}` and :math:`F_{cur}` and of the spectrum of :math:`X`. The mean
coherence for the segment is defined as the mean of :math:`C(\\nu)` in the
frequency range of interest. The time-delay between the two cross correlations
is found in the unwrapped phase, :math:`\phi(\nu)`, of the cross spectrum and is
linearly proportional to frequency:

:math:`\phi_j = m. \nu_j, m = 2 \pi \delta t`

The time shift for each window between two signals is the slope :math:`m` of a
weighted linear regression of the samples within the frequency band of interest.
The weights are those introduced by [Clarke2011]_,
which incorporate both the cross-spectral amplitude and cross-coherence, unlike
[Poupinet1984]_. The errors are estimated using the weights (thus the
coherence) and the squared misfit to the modelled slope:

:math:`e_m = \sqrt{\sum_j{(\\frac{w_j \\nu_j}{\sum_i{w_i \\nu_i^2}})^2}\sigma_{\phi}^2}`

where :math:`w` are weights, :math:`\\nu` are cross-coherences and
:math:`\sigma_{\phi}^2` is the squared misfit of the data to the modelled slope
and is calculated as :math:`\sigma_{\phi}^2 = \\frac{\sum_j(\phi_j - m \\nu_j)^2}{N-1}`

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

    window_length_samples = np.int(window_length * df)
    # try:
    #     from sf.helper import next_fast_len
    # except ImportError:
    #     from obspy.signal.util import next_pow_2 as next_fast_len
    from msnoise.api import nextpow2
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

        minind += int(step*df)
        maxind += int(step*df)

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
        time_axis.append(tmin+window_length/2.+count*step)
        count += 1

        del fcur, fref
        del X
        del freq_vec
        del index_range
        del w, v, e, s2x2, sx2, m, em

    if maxind > len(current) + step*df:
        logging.warning("The last window was too small, but was computed")

    return np.array([time_axis, delta_t, delta_err, delta_mcoh]).T

