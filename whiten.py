import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import warnings
warnings.simplefilter("ignore")
# from tempfile import mkdtemp
# cachedir = mkdtemp()
# from joblib import Memory
# memory = Memory(cachedir=cachedir, verbose=0)


def nextpow2(x):
    return np.ceil(np.log2(np.abs(x)))


# @memory.cache
def whiten(matsign, Nfft, tau, frebas, frehaut, plot=False, tmp=False):
    """This function takes 1-dimensional *matsign* timeseries array,
    goes to frequency domain using fft, whitens the amplitude of the spectrum
    in frequency domain between *frebas* and *frehaut*
    and returns the whitened fft.

    Parameters
    ----------
    matsign : numpy.ndarray
        Contains the 1D time series to whiten
    Nfft : int
        The number of points to compute the FFT
    tau : int
        The sampling frequency of the `matsign`
    frebas : int
        The lower frequency bound
    frehaut : int
        The upper frequency bound
    plot : bool
        Whether to show a raw plot of the action (default: False)
    tmp : int
        unused.


    Returns
    -------
    data : numpy.ndarray
        The FFT of the input trace, whitened between the frequency bounds
"""
    # if len(matsign)/2 %2 != 0:
        # matsign = np.append(matsign,[0,0])

    if plot:
        plt.subplot(411)
        plt.plot(np.arange(len(matsign)) * tau, matsign)
        plt.xlim(0, len(matsign) * tau)
        plt.title('Input trace')

    Napod = 300

    freqVec = np.arange(0., Nfft / 2.0) / (tau * (Nfft - 1))
    J = np.where((freqVec >= frebas) & (freqVec <= frehaut))[0]
    low = J[0] - Napod
    if low < 0:
        low = 0

    porte1 = J[0]
    porte2 = J[-1]
    high = J[-1] + Napod
    if high > Nfft / 2:
        high = Nfft / 2

    FFTRawSign = scipy.fftpack.fft(matsign, Nfft)

    if plot:
        plt.subplot(412)
        axis = np.arange(len(FFTRawSign))
        plt.plot(axis[1:], np.abs(FFTRawSign[1:]))
        plt.xlim(0, max(axis))
        plt.title('FFTRawSign')

    # Apodisation a gauche en cos2
    FFTRawSign[0:low] *= 0
    FFTRawSign[low:porte1] = np.cos(np.linspace(np.pi / 2., np.pi, porte1 - low)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[low:porte1]))
    # Porte
    FFTRawSign[porte1:porte2] = np.exp(
        1j * np.angle(FFTRawSign[porte1:porte2]))
    # Apodisation a droite en cos2
    FFTRawSign[porte2:high] = np.cos(np.linspace(0., np.pi / 2., high - porte2)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[porte2:high]))

    if low == 0:
        low = 1

    FFTRawSign[-low:] *= 0
    # Apodisation a gauche en cos2
    FFTRawSign[-porte1:-low] = np.cos(np.linspace(0., np.pi / 2., porte1 - low)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[-porte1:-low]))
    # Porte
    FFTRawSign[-porte2:-porte1] = np.exp(
        1j * np.angle(FFTRawSign[-porte2:-porte1]))
    # ~ # Apodisation a droite en cos2
    FFTRawSign[-high:-porte2] = np.cos(np.linspace(np.pi / 2., np.pi, high - porte2)) ** 2 * np.exp(
        1j * np.angle(FFTRawSign[-high:-porte2]))

    FFTRawSign[high:-high] *= 0

    FFTRawSign[-1] *= 0.
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

    wmatsign = np.real(scipy.fftpack.ifft(FFTRawSign))
    del matsign
    if plot:
        plt.subplot(414)
        plt.plot(np.arange(len(wmatsign)) * tau, wmatsign)
        plt.xlim(0, len(wmatsign) * tau)
        plt.show()
    # return wmatsign
    return FFTRawSign


if __name__ == '__main__':
    import time
    N = 86438

    a = np.random.random(N)
    a = np.sin(a) + np.sin(a / 4.) + np.sin(a / 16.)

    t = time.time()
    whiten(a, N, 0.05, 1.0, 5.9, plot=False)
    print time.time() - t

    t = time.time()
    whiten(a, N, 0.05, 1.0, 5.9, plot=False)
    print time.time() - t
