# -*- coding: utf-8 -*-
"""
Compute the H/V spectral ratio from the imaginary part of the CCFs
==================================================================

"""

###########################################################################
# Import & getting the data from the computed autocorrelations (ZZ, EE, NN)

import os
if "SPHINX_DOC_BUILD" in os.environ or 1:
    os.chdir(r"C:\tmp\MSNOISE_DOC")

import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import datetime
register_matplotlib_converters()

plt.style.use("ggplot")

from msnoise.api import connect, get_results, build_movstack_datelist, get_params, get_t_axis, xr_get_ccf


# connect to the database
db = connect()

# Obtain a list of dates between ``start_date`` and ``enddate``
start, end, datelist = build_movstack_datelist(db)

# Get the list of parameters from the DB:
params = get_params(db)

# Get the time axis for plotting the CCF:
taxis = get_t_axis(db)

# Get the first mov_stack configured:
mov_stack = params.mov_stack[0]

# Get the results for two station, filter id=1, ZZ component, mov_stack=1 and the results as a 2D array:
ZZ = xr_get_ccf("PF.FJS.00", "PF.FJS.00", "ZZ", 2, mov_stack, taxis, format="dataframe")
EE = xr_get_ccf("PF.FJS.00", "PF.FJS.00", "EE", 2, mov_stack, taxis, format="dataframe")
NN = xr_get_ccf("PF.FJS.00", "PF.FJS.00", "NN", 2, mov_stack, taxis, format="dataframe")

##############################################################
# Checking the data has the same time base (index) and joining


r = "1D"
rZZ = ZZ.resample(r).mean()
rEE = EE.resample(r).mean()
rNN = NN.resample(r).mean()

merged = pd.concat({'ZZ': rZZ, 'EE': rEE, 'NN': rNN}, axis=0)
merged.index.names = ['channel', 'date']

# Swap the levels of the MultiIndex
result = merged.swaplevel('channel', 'date').sort_index(level='date') #.iloc[:500]

##############################################################
# Helper functions


def get_imag(d, fs):
    NFFT = 1024*32
    iX = np.imag(np.fft.fft(d, n=NFFT)[:NFFT//2])
    freqs = np.fft.fftfreq(NFFT, d=1/fs)[:NFFT//2]

    return freqs, iX

##############################################################
# Doing the job !


HVSRs = {}
for date, group in result.groupby(level='date'):
    print(f"Date: {date}")
    group= group.droplevel(0)
    try:
        f, iZZ = get_imag(group.loc["ZZ"].values, 20)
        f, iEE = get_imag(group.loc["EE"].values, 20)
        f, iNN = get_imag(group.loc["NN"].values, 20)
    except:
        continue
    hvsr = np.sqrt((iEE+iNN)/iZZ)
    HVSRs[date] = hvsr

HVSR = pd.DataFrame(HVSRs, index=f)
X = HVSR.copy().fillna(0.0)
X = X.T.resample(r).mean().T


##############################################################
# Plotting

plt.subplots(figsize=(8,10))
plt.pcolormesh(X.index, X.columns, X.T, rasterized=True, vmax=5)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Date")
plt.grid()
plt.xlim(1,10)
plt.show()


##############################################################
# Plotting & zooming around 1-4 Hz

plt.subplots(figsize=(8,10))
plt.pcolormesh(X.index, X.columns, X.T, rasterized=True, vmax=3)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Date")
plt.grid()
plt.xlim(1,4)
plt.show()

##############################################################
# Plotting the HVSR curve

plt.subplots(figsize=(8,10))
hvsr = X.mean(axis=1)
plt.plot(hvsr.index, hvsr)
plt.show()


#EOF