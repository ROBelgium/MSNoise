import matplotlib.pyplot as plt
import numpy as np
from MSNoise.myCorr import myCorr

data = np.random.random((2,2048))+0.5
maxlag = 100
corr = myCorr(data,maxlag)
timelags = np.arange(-maxlag,maxlag+1)
plt.plot(timelags, corr)
plt.show()