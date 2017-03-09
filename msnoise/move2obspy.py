import logging

import matplotlib.pyplot as plt
import numpy as np
import scipy.fftpack
import scipy.optimize
import scipy.signal
import statsmodels.api as sm
from obspy.signal.invsim import cosine_taper
from scipy.fftpack.helper import next_fast_len

from .api import nextpow2


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
    if nfft is None:
        s1 = np.array(data[0].shape)
        shape = s1 - 1
        # Speed up FFT by padding to optimal size for FFTPACK
        fshape = [next_fast_len(int(d)) for d in shape]
        nfft = fshape[0]

    normalized = True
    allCpl = False

    maxlag = np.round(maxlag)
    # ~ print "np.shape(data)",np.shape(data)
    if data.shape[0] == 2:
        # ~ print "2 matrix to correlate"
        if allCpl:
            # Skipped this unused part
            pass
        else:
            K = data.shape[0]
            # couples de stations
            couples = np.concatenate((np.arange(0, K), K + np.arange(0, K)))

    Nt = data.shape[1]
    Nc = 2 * Nt - 1

    # corr = scipy.fftpack.fft(data,int(Nfft),axis=1)
    corr = data

    if plot:
        plt.subplot(211)
        plt.plot(np.arange(len(corr[0])) * 0.05, np.abs(corr[0]))
        plt.subplot(212)
        plt.plot(np.arange(len(corr[1])) * 0.05, np.abs(corr[1]))

    corr = np.conj(corr[couples[0]]) * corr[couples[1]]
    corr = np.real(scipy.fftpack.ifft(corr, nfft)) / (Nt)
    corr = np.concatenate((corr[-Nt + 1:], corr[:Nt + 1]))

    if plot:
        plt.figure()
        plt.plot(corr)

    if normalized:
        E = np.real(np.sqrt(
            np.mean(scipy.fftpack.ifft(data, n=nfft, axis=1) ** 2, axis=1)))
        normFact = E[0] * E[1]
        corr /= np.real(normFact)

    if maxlag != Nt:
        tcorr = np.arange(-Nt + 1, Nt)
        dN = np.where(np.abs(tcorr) <= maxlag)[0]
        corr = corr[dN]

    del data
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

    if plot:
        plt.subplot(411)
        plt.plot(np.arange(len(data)) * delta, data)
        plt.xlim(0, len(data) * delta)
        plt.title('Input trace')

    Napod = 100
    Nfft = int(Nfft)
    freqVec = scipy.fftpack.fftfreq(Nfft, d=delta)[:Nfft // 2]

    J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
    low = J[0] - Napod
    if low <= 0:
        low = 1

    porte1 = J[0]
    porte2 = J[-1]
    high = J[-1] + Napod
    if high > Nfft / 2:
        high = int(Nfft // 2)

    FFTRawSign = scipy.fftpack.fft(data, Nfft)

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

        wdata = np.real(scipy.fftpack.ifft(FFTRawSign, Nfft))
        plt.subplot(414)
        plt.plot(np.arange(len(wdata)) * delta, wdata)
        plt.xlim(0, len(wdata) * delta)
        plt.show()
    return FFTRawSign


def smooth(x, window='boxcar', half_win=3):
    """ some window smoothing """
    window_len = 2 * half_win + 1
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = np.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    if window == "boxcar":
        w = scipy.signal.boxcar(window_len).astype('complex')
    else:
        w = scipy.signal.hanning(window_len).astype('complex')
    y = np.convolve(w / w.sum(), s, mode='valid')
    return y[half_win:len(y) - half_win]


def getCoherence(dcs, ds1, ds2):
    n = len(dcs)
    coh = np.zeros(n).astype('complex')
    valids = np.argwhere(np.logical_and(np.abs(ds1) > 0, np.abs(ds2 > 0)))
    coh[valids] = dcs[valids] / (ds1[valids] * ds2[valids])
    coh[coh > (1.0 + 0j)] = 1.0 + 0j
    return coh


def mwcs_old(ccCurrent, ccReference, fmin, fmax, sampRate, tmin, windL, step,
         plot=False):
    """...

    :type ccCurrent: :class:`numpy.ndarray`
    :param ccCurrent: The "Current" timeseries
    :type ccReference: :class:`numpy.ndarray`
    :param ccReference: The "Reference" timeseries
    :type fmin: float
    :param fmin: The lower frequency bound to compute the dephasing
    :type fmax: float
    :param fmax: The higher frequency bound to compute the dephasing
    :type sampRate: float
    :param sampRate: The sample rate of the input timeseries
    :type tmin: float
    :param tmin: The leftmost time lag (used to compute the "time lags array")
    :type windL: float
    :param windL: The moving window length
    :type step: float
    :param step: The step to jump for the moving window
    :type plot: bool
    :param plot: If True, plots the MWCS result for each window. Defaults to
        False

    :rtype: :class:`numpy.ndarray`
    :returns: [Taxis,deltaT,deltaErr,deltaMcoh]. Taxis contains the central
        times of the windows. The three other columns contain dt, error and
        mean coherence for each window.
    """

    windL = np.int(windL * sampRate)
    step = np.int(step * sampRate)
    count = 0
    deltaT = []
    deltaErr = []
    deltaMcoh = []
    Taxis = []
    # padd = 2 ** (nextpow2(windL) + 2)
    padd = next_fast_len(windL)

    # Tentative checking if enough point are used to compute the FFT
    freqVec = scipy.fftpack.fftfreq(int(padd), 1. / sampRate)[:int(padd) // 2]
    indRange = np.argwhere(np.logical_and(freqVec >= fmin,
                                          freqVec <= fmax))
    if len(indRange) < 2:
        padd = 2 ** (nextpow2(windL) + 3)

    tp = cosine_taper(windL, .85)

    timeaxis = (np.arange(len(ccCurrent)) / float(sampRate)) + tmin
    minind = 0
    maxind = windL
    while maxind <= len(ccCurrent):
        ind = minind

        cci = ccCurrent[ind:(ind + windL)].copy()
        cci = scipy.signal.detrend(cci, type='linear')
        # cci -= cci.min()
        # cci /= cci.max()
        # cci -= np.mean(cci)
        cci *= tp

        cri = ccReference[ind:(ind + windL)].copy()
        cri = scipy.signal.detrend(cri, type='linear')
        # cri -= cri.min()
        # cri /= cri.max()
        # cri -= np.mean(cri)
        cri *= tp

        Fcur = scipy.fftpack.fft(cci, n=int(padd))[:int(padd) // 2]
        Fref = scipy.fftpack.fft(cri, n=int(padd))[:int(padd) // 2]

        Fcur2 = np.real(Fcur) ** 2 + np.imag(Fcur) ** 2
        Fref2 = np.real(Fref) ** 2 + np.imag(Fref) ** 2

        smoother = 5

        dcur = np.sqrt(smooth(Fcur2, window='hanning', half_win=smoother))
        dref = np.sqrt(smooth(Fref2, window='hanning', half_win=smoother))

        # Calculate the cross-spectrum
        X = Fref * (Fcur.conj())
        X = smooth(X, window='hanning', half_win=smoother)
        dcs = np.abs(X)

        # Find the values the frequency range of interest
        freqVec = scipy.fftpack.fftfreq(len(X) * 2, 1. / sampRate)[
                  :int(padd) // 2]
        indRange = np.argwhere(np.logical_and(freqVec >= fmin,
                                              freqVec <= fmax))

        # Get Coherence and its mean value
        coh = getCoherence(dcs, dref, dcur)
        mcoh = np.mean(coh[indRange])

        # Get Weights
        w = 1.0 / (1.0 / (coh[indRange] ** 2) - 1.0)
        w[coh[indRange] >= 0.99] = 1.0 / (1.0 / 0.9801 - 1.0)
        w = np.sqrt(w * np.sqrt(dcs[indRange]))
        # w /= (np.sum(w)/len(w)) #normalize
        w = np.real(w)

        # Frequency array:
        v = np.real(freqVec[indRange]) * 2 * np.pi
        vo = np.real(freqVec) * 2 * np.pi

        # Phase:
        phi = np.angle(X)
        phi[0] = 0.
        phi = np.unwrap(phi)
        # phio = phi.copy()
        phi = phi[indRange]

        # Calculate the slope with a weighted least square linear regression
        # forced through the origin
        # weights for the WLS must be the variance !
        res = sm.regression.linear_model.WLS(phi, v, w ** 2).fit()

        # print "forced", np.real(res.params[0])
        # print "!forced", np.real(res2.params[0])

        m = np.real(res.params[0])
        deltaT.append(m)

        # print phi.shape, v.shape, w.shape
        e = np.sum((phi - m * v) ** 2) / (np.size(v) - 1)
        s2x2 = np.sum(v ** 2 * w ** 2)
        sx2 = np.sum(w * v ** 2)
        e = np.sqrt(e * s2x2 / sx2 ** 2)
        # print w.shape
        if plot:
            plt.figure()
            plt.suptitle('%.1fs' % (timeaxis[ind + windL // 2]))
            plt.subplot(311)
            plt.plot(cci)
            plt.plot(cri)
            ax = plt.subplot(312)
            plt.plot(vo / (2 * np.pi), phio)
            plt.scatter(v / (2 * np.pi), phi, c=w, edgecolor='none',
                        vmin=0.6, vmax=1)
            plt.subplot(313, sharex=ax)
            plt.plot(v / (2 * np.pi), coh[indRange])
            plt.axhline(mcoh, c='r')
            plt.axhline(1.0, c='k', ls='--')
            plt.xlim(-0.1, 1.5)
            plt.ylim(0, 1.5)
            plt.show()

        deltaErr.append(e)
        deltaMcoh.append(np.real(mcoh))
        Taxis.append(timeaxis[ind + windL // 2])
        count += 1

        minind += step
        maxind += step
        del Fcur, Fref
        del X
        del freqVec
        del indRange
        del w, v, e, s2x2
        del res

    if maxind > len(ccCurrent) + step:
        logging.warning("The last window was too small, but was computed")

    return np.array([Taxis, deltaT, deltaErr, deltaMcoh]).T


def linear_regression(xdata, ydata, weights=None, p0 = None, intercept=False):
    """ Use non-linear least squares to fit a function, f, to data. This method
    is a generalized version of :meth:`scipy.optimize.minpack.curve_fit`;
    allowing for:

    * OLS without intercept : ``linear_regression(xdata, ydata)``
    * OLS with intercept : ``linear_regression(xdata, ydata, intercept=True)``
    * WLS without intercept : ``linear_regression(xdata, ydata, weights)``
    * WLS with intercept : ``linear_regression(xdata, ydata, weights, intercept=True)``

    If the expected values of slope (and intercept) are different from 0.0,
    provide the p0 value(s).

    :param xdata: The independent variable where the data is measured.
    :param ydata: The dependent data - nominally f(xdata, ...)
    :param weights: If not None, the uncertainties in the ydata array. These are
     used as weights in the least-squares problem. If None, the uncertainties
     are assumed to be 1. In SciPy vocabulary, our weights are 1/sigma.
    :param p0: Initial guess for the parameters. If None, then the initial
     values will all be 0 (Different from SciPy where all are 1)
    :param intercept: If False: solves y=a*x ; if True: solves y=a*x+b.
    :return:
    """
    if weights is not None:
        sigma = 1./weights
    else:
        sigma = None
    if intercept:
        p, cov = scipy.optimize.curve_fit(lambda x, a, b: a * x + b,
                                          xdata, ydata,
                                          [0, 0], sigma=sigma,
                                          absolute_sigma=False,
                                          xtol=1e-20)
        slope, intercept = p
        std_slope = np.sqrt(cov[0, 0])
        std_intercept = np.sqrt(cov[1, 1])
        return slope, intercept, std_slope, std_intercept

    else:
        p, cov = scipy.optimize.curve_fit(lambda x, a: a * x,
                                          xdata, ydata,
                                          0, sigma=sigma,
                                          absolute_sigma=False,
                                          xtol=1e-20)
        slope = p[0]
        std_slope = np.sqrt(cov[0, 0])
        return slope, std_slope


def mwcs(current, reference, freqmin, freqmax, df, tmin, window_length, step, smoothing_half_win=5):
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

:math:`C(\\nu) = \\frac{|\overline{X(\\nu))}|}{\sqrt{|\overline{F_{ref}(\\nu)}|^2.|\overline{F_{cur}(\\nu)}|^2}}`


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

:math:`e_m = \sqrt{\sum_j{\\frac{w_j \\nu_j}{\sum_i{w_i \\nu_i^2}}^2}\sigma_{\phi}^2}`

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
:param smoothing_half_win: If different from 0, defines the half length of the
    smoothing hanning window.


:rtype: :class:`numpy.ndarray`
:returns: [Taxis,deltaT,deltaErr,deltaMcoh]. Taxis contains the central
    times of the windows. The three other columns contain dt, error and
    mean coherence for each window.
    """
    deltaT = []
    deltaErr = []
    deltaMcoh = []
    Taxis = []

    windL2 = np.int(window_length * df)
    padd = next_fast_len(windL2)

    count = 0
    tp = cosine_taper(windL2, 0.85)
    minind = 0
    maxind = windL2
    while maxind <= len(current):
        cci = current[minind:(minind + windL2)]
        cci = scipy.signal.detrend(cci, type='linear')
        cci *= tp

        cri = reference[minind:(minind + windL2)]
        cri = scipy.signal.detrend(cri, type='linear')
        cri *= tp

        minind += int(step*df)
        maxind += int(step*df)

        Fcur = scipy.fftpack.fft(cci, n=padd)[:padd // 2]
        Fref = scipy.fftpack.fft(cri, n=padd)[:padd // 2]

        Fcur2 = np.real(Fcur) ** 2 + np.imag(Fcur) ** 2
        Fref2 = np.real(Fref) ** 2 + np.imag(Fref) ** 2

        # Calculate the cross-spectrum
        X = Fref * (Fcur.conj())
        if smoothing_half_win != 0:
            dcur = np.sqrt(smooth(Fcur2, window='hanning', half_win=smoothing_half_win))
            dref = np.sqrt(smooth(Fref2, window='hanning', half_win=smoothing_half_win))
            X = smooth(X, window='hanning', half_win=smoothing_half_win)
        else:
            dcur = np.sqrt(Fcur2)
            dref = np.sqrt(Fref2)

        dcs = np.abs(X)

        # Find the values the frequency range of interest
        freq_vec = scipy.fftpack.fftfreq(len(X) * 2, 1. / df)[:padd // 2]
        index_range = np.argwhere(np.logical_and(freq_vec >= freqmin,
                                                 freq_vec <= freqmax))

        # Get Coherence and its mean value
        coh = getCoherence(dcs, dref, dcur)
        mcoh = np.mean(coh[index_range])

        # Get Weights
        w = 1.0 / (1.0 / (coh[index_range] ** 2) - 1.0)
        w[coh[index_range] >= 0.99] = 1.0 / (1.0 / 0.9801 - 1.0)
        w = np.sqrt(w * np.sqrt(dcs[index_range]))
        # w /= (np.sum(w)/len(w)) #normalize
        w = np.real(w)

        # Frequency array:
        v = np.real(freq_vec[index_range]) * 2 * np.pi
        vo = np.real(freq_vec) * 2 * np.pi

        # Phase:
        phi = np.angle(X)
        phi[0] = 0.
        phi = np.unwrap(phi)
        # phio = phi.copy()
        phi = phi[index_range]

        # Calculate the slope with a weighted least square linear regression
        # forced through the origin
        # weights for the WLS must be the variance !
        m, em = linear_regression(v.flatten(), phi.flatten(), w.flatten())

        deltaT.append(m)

        # print phi.shape, v.shape, w.shape
        e = np.sum((phi - m * v) ** 2) / (np.size(v) - 1)
        s2x2 = np.sum(v ** 2 * w ** 2)
        sx2 = np.sum(w * v ** 2)
        e = np.sqrt(e * s2x2 / sx2 ** 2)

        deltaErr.append(e)
        deltaMcoh.append(np.real(mcoh))
        Taxis.append(tmin+window_length/2.+count*step)
        count += 1

        del Fcur, Fref
        del X
        del freq_vec
        del index_range
        del w, v, e, s2x2, sx2, m, em

    if maxind > len(current) + step*df:
        logging.warning("The last window was too small, but was computed")

    return np.array([Taxis, deltaT, deltaErr, deltaMcoh]).T

