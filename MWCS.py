"""
This function gives the relative travel time variations (dt/t) between cross-
correlations along with the associated error in the measurement. you can set
the frequency range in which you want to determine dt/t, i.e. you can measure
dt/t in a narrow frequency band even for a broadband signal.
The function follows the procedure set out in Clarke et al. (2011)."""

import statsmodels.api as sm
from obspy.signal import cosTaper
import scipy.signal
import numpy as np
# import matplotlib.pyplot as plt
import logging
import scipy.fftpack


def nextpow2(x):
    return np.ceil(np.log2(np.abs(x)))


def smooth(x, window='boxcar'):
    """ some window smoothing """
    half_win = 11
    window_len = 2*half_win+1
    # extending the data at beginning and at the end
    # to apply the window at the borders
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    if window == "boxcar":
        w = scipy.signal.boxcar(window_len).astype('complex')
    else:
        w = scipy.signal.hanning(window_len).astype('complex')

    y = np.convolve(w/w.sum(), s, mode='valid')
    return y[half_win:len(y)-half_win]


def getCoherence(dcs, ds1, ds2):
    # ! number of data
    n = len(dcs)
    coh = np.zeros(n).astype('complex')
    # ! calculate coherence
    valids = np.argwhere(np.logical_and(np.abs(ds1) > 0, np.abs(ds2 > 0)))
    coh[valids] = dcs[valids]/(ds1[valids]*ds2[valids])
    coh[coh > (1.0+0j)] = 1.0+0j
    return coh


def mwcs(ccCurrent, ccReference, fmin, fmax, sampRate, tmin, windL, step):
    """...

    Parameters
    ----------
    ccCurrent : numpy.ndarray
        The "Current" timeseries
    ccReference : numpy.ndarray
        The "Reference" timeseries
    fmin : int
        The lower frequency bound to compute the dephasing
    fmax : int
        The higher frequency bound to compute the dephasing
    sampRate : int
        The sample rate of the input timeseries
    tmin : int
        The leftmost time lag (used to compute the "time lags array")
    windL : int
        The moving window length
    step : int
        The step to jump for the moving window


    Returns
    -------
    data : numpy.ndarray
        Taxis,deltaT,deltaErr,deltaMcoh"""

    windL = np.int(windL*sampRate)
    step = np.int(step*sampRate)
    count = 0
    deltaT = []
    deltaErr = []
    deltaMcoh = []
    Taxis = []
    padd = 2**(nextpow2(windL)+1)
    tp = cosTaper(windL, 0.02)
    timeaxis = (np.arange(len(ccCurrent)) / float(sampRate))+tmin

    minind = 0
    maxind = windL
    while maxind <= len(ccCurrent):
        ind = minind
        cci = ccCurrent[ind:(ind+windL)].copy()
        cci -= np.mean(cci)
        cci *= tp

        cri = ccReference[ind:(ind+windL)].copy()
        cri -= np.mean(cri)
        cri *= tp

        Fcur = scipy.fftpack.fft(cci, n=int(padd))[:padd/2]
        Fref = scipy.fftpack.fft(cri, n=int(padd))[:padd/2]

        Fcur2 = np.real(Fcur)**2 + np.imag(Fcur)**2
        Fref2 = np.real(Fref)**2 + np.imag(Fref)**2

        dcur = np.sqrt(smooth(Fcur2, window='hanning'))
        dref = np.sqrt(smooth(Fref2, window='hanning'))

        #Calculate the cross-spectrum
        X = Fref*(Fcur.conj())
        X = smooth(X, window='hanning')
        dcs = np.abs(X)

        #Find the values the frequency range of interest
        freqVec = scipy.fftpack.fftfreq(len(X)*2, 1./sampRate)[:padd/2]
        indRange = np.argwhere(np.logical_and(freqVec >= fmin,
                                              freqVec <= fmax))

        # Get Coherence and its mean value
        coh = getCoherence(dcs, dref, dcur)
        mcoh = np.mean(coh[indRange])

        #Get Weights
        w = 1.0 / (1.0 / (coh[indRange]**2) - 1.0)
        w[coh[indRange] >= 0.99] = 1.0 / (1.0 / 0.9801 - 1.0)
        w = np.sqrt(w * np.sqrt(dcs[indRange]))
        # w /= (np.sum(w)/len(w)) #normalize
        w = np.real(w)

        #Frequency array:
        v = np.real(freqVec[indRange]) * 2 * np.pi
        # vo = np.real(freqVec) * 2 * np.pi

        #Phase:
        phi = np.angle(X)
        phi[0] = 0
        phi = np.unwrap(phi)
        # print phi[0]
        # phio = phi
        phi = phi[indRange]

        #Calculate the slope with a weighted least square linear regression
        #forced through the origin
        #weights for the WLS must be the variance !
        res = sm.regression.linear_model.WLS(phi, v, w**2).fit()

        m = np.real(res.params[0])
        deltaT.append(m)

        # print phi.shape, v.shape, w.shape
        e = np.sum((phi-m*v)**2) / (np.size(v)-1)
        s2x2 = np.sum(v**2 * w**2)
        sx2 = np.sum(w * v**2)
        e = np.sqrt(e*s2x2 / sx2**2)
        # print w.shape
        # plt.plot(vo, phio)
        # plt.scatter(v,phi,c=w)
        # plt.plot(vo,vo*m)
        # plt.xlim(-1,10)
        # plt.ylim(-5,5)
        # plt.show()

        deltaErr.append(e)
        # print m, e, res.bse[0]
        deltaMcoh.append(np.real(mcoh))
        Taxis.append(timeaxis[ind + windL/2])
        count = count + 1

        minind += step
        maxind += step
        del Fcur, Fref
        del X
        del freqVec
        del indRange
        del w, v, e, s2x2
        del res

    if maxind > len(ccCurrent)+step:
        logging.warning("The last window was too small, but was computed")

    return np.array([Taxis, deltaT, deltaErr, deltaMcoh]).T
