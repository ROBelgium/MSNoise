import matplotlib.pyplot as plt
import numpy as np
import scipy.fftpack
from msnoise.myCorr import myCorr


a = np.random.random(2048)-0.5
b = np.random.random(2048)-0.5

a = scipy.fftpack.fft(a)[:1024]
b = scipy.fftpack.fft(b)[:1024]
data = np.array([a,b])
maxlag = 100
corr = myCorr(data,maxlag)
timelags = np.arange(-maxlag,maxlag+1)
plt.plot(timelags, corr)
plt.show()