import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

from .api import nextpow2

def myCorr(data, maxlag, plot=False):
    """This function takes ndimensional *data* array, computes the cross-correlation in the frequency domain
    and returns the cross-correlation function between [-*maxlag*:*maxlag*].

    :type data: :class:`numpy.ndarray`
    :param data: This array contains the fft of each timeseries to be cross-correlated.
    :type maxlag: int
    :param maxlag: This number defines the number of samples (N=2*maxlag + 1) of the CCF that will be returned.

    :rtype: :class:`numpy.ndarray`
    :returns: The cross-correlation function between [-maxlag:maxlag]
    """

    normalized = True
    allCpl = False

    maxlag = np.round(maxlag)
    #~ print "np.shape(data)",np.shape(data)
    if data.shape[0] == 2:
        #~ print "2 matrix to correlate"
        if allCpl:
            # Skipped this unused part
            pass
        else:
            K = data.shape[0]
            #couples de stations
            couples = np.concatenate((np.arange(0, K), K + np.arange(0, K)))

    Nt = data.shape[1]
    Nc = 2 * Nt - 1
    Nfft = 2 ** nextpow2(Nc)

    # corr = scipy.fftpack.fft(data,int(Nfft),axis=1)
    corr = data

    if plot:
            plt.subplot(211)
            plt.plot(np.arange(len(corr[0])) * 0.05, np.abs(corr[0]))
            plt.subplot(212)
            plt.plot(np.arange(len(corr[1])) * 0.05, np.abs(corr[1]))

    corr = np.conj(corr[couples[0]]) * corr[couples[1]]
    corr = np.real(scipy.fftpack.ifft(corr)) / Nt
    corr = np.concatenate((corr[-Nt + 1:], corr[:Nt + 1]))

    if plot:
        plt.figure()
        plt.plot(corr)

    E = np.sqrt(np.mean(scipy.fftpack.ifft(data, axis=1) ** 2, axis=1))
    normFact = E[0] * E[1]

    if normalized:
        corr /= np.real(normFact)

    if maxlag != Nt:
        tcorr = np.arange(-Nt + 1, Nt)
        dN = np.where(np.abs(tcorr) <= maxlag)[0]
        corr = corr[dN]

    del data
    return corr

if __name__ == "__main__":

    import time

    data = np.random.random((2, 86400 * 20))
    print(data.shape)
    t = time.time()
    corr = myCorr(data, maxlag=25, plot=False)
    print( "Time:", time.time() - t)
    print(np.mean(corr))

    plt.figure()
    plt.plot(corr)
    plt.axhline(1)
    plt.show()
