import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack

def nextpow2(x):
    return np.ceil(np.log2(np.abs(x)))


def whiten(data, Nfft, delta, freqmin, freqmax, plot=False):
    """This function takes 1-dimensional *data* timeseries array,
    goes to frequency domain using fft, whitens the amplitude of the spectrum
    in frequency domain between *freqmin* and *freqmax*
    and returns the whitened fft.

    Parameters
    ----------
    data : numpy.ndarray
        Contains the 1D time series to whiten
    Nfft : int
        The number of points to compute the FFT
    delta : int
        The sampling frequency of the `data`
    freqmin : int
        The lower frequency bound
    freqmax : int
        The upper frequency bound
    plot : bool
        Whether to show a raw plot of the action (default: False)

    Returns
    -------
    data : numpy.ndarray
        The FFT of the input trace, whitened between the frequency bounds
"""

    if plot:
        plt.subplot(411)
        plt.plot(np.arange(len(data)) * delta, data)
        plt.xlim(0, len(data) * delta)
        plt.title('Input trace')

    Napod = 100
    Nfft = int(Nfft)
    freqVec = scipy.fftpack.fftfreq(Nfft,d=delta)[:Nfft/2]
    
    J = np.where((freqVec >= freqmin) & (freqVec <= freqmax))[0]
    low = J[0] - Napod
    if low <= 0:
        low = 1

    porte1 = J[0]
    porte2 = J[-1]
    high = J[-1] + Napod
    if high > Nfft / 2:
        high = Nfft // 2

    FFTRawSign = scipy.fftpack.fft(data, Nfft)
    
    if plot:
        plt.subplot(412)
        axis = np.arange(len(FFTRawSign))
        plt.plot(axis[1:], np.abs(FFTRawSign[1:]))
        plt.xlim(0, max(axis))
        plt.title('FFTRawSign')

    # Left tapering:
    FFTRawSign[0:low] *= 0
    FFTRawSign[low:porte1] = np.cos(np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * np.exp(1j * np.angle(FFTRawSign[low:porte1]))
    # Pass band:
    FFTRawSign[porte1:porte2] = np.exp(1j * np.angle(FFTRawSign[porte1:porte2]))
    # Right tapering:
    FFTRawSign[porte2:high] = np.cos(np.linspace(0., np.pi / 2., high - porte2)) ** 2 * np.exp(1j * np.angle(FFTRawSign[porte2:high]))
    FFTRawSign[high:Nfft+1] *= 0
    
    # Hermitian symmetry (because the input is real)
    FFTRawSign[-Nfft/2+1:] = FFTRawSign[1:Nfft/2].conjugate()[::-1]

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

        wdata = np.real(scipy.fftpack.ifft(FFTRawSign))
        plt.subplot(414)
        plt.plot(np.arange(len(wdata)) * delta, wdata)
        plt.xlim(0, len(wdata) * delta)
        plt.show()
    
    return FFTRawSign

    
if __name__ == '__main__':
    import time
    N = 2048
    np.random.seed(1234)
    a = np.random.random(N)
    a = np.sin(a) + np.sin(a / 4.) + np.sin(a / 16.)
    a -= a.mean()
    t = time.clock()
    for i in range(1000):
        whiten(a.copy(), N, 0.05, 1.0, 5.9, plot=False)
    print "1000 loops:", (time.clock()-t) * 1000, "ms"
