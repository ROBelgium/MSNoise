import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack


def nextpow2(x):
    return np.ceil(np.log2(np.abs(x)))


def myCorr(data, maxlag, plot=False):
    normalized = True
    allCpl = False

    maxlag = np.round(maxlag)
    #~ print "np.shape(data)", np.shape(data)
    if np.shape(data)[0] == 2:
        #~ print "2 matrices to correlate"
        if allCpl:
            # Skipped this unused part
            pass
        else:
            K = np.shape(data)[0]
            # couples de stations
            couples = np.concatenate((np.arange(0, K), K + np.arange(0, K)))

    Nt = np.shape(data)[1]
    Nc = 2 * Nt - 1
    Nfft = 2 ** nextpow2(Nc)

    corr = scipy.fftpack.fft(data, int(Nfft), axis=1)

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

    E = np.sqrt(np.mean(data ** 2, axis=1))
    normFact = E[0] * E[1]

    if normalized:
        corr /= normFact

    if maxlag != Nt:
        tcorr = np.arange(-Nt + 1, Nt)
        dN = np.where(np.abs(tcorr) <= maxlag)[0]
        corr = corr[dN]

    del data
    return corr

if __name__ == "__main__":

    import time

    data = np.random.random((2, 86400 * 20))
    print data.shape
    t = time.time()
    corr = myCorr(data, maxlag=25, plot=False)
    print("Time: %.3f" % (time.time() - t))
    print(np.mean(corr))

    plt.figure()
    plt.plot(corr)
    plt.xlabel("time lag")
    plt.ylabel("covariance")
    plt.axhline(1)
    plt.show()
