# -*- coding: utf-8 -*-
"""
Compute the H/V spectral ratio from the imaginary part of the CCFs
==================================================================

"""

###########################################################################
# Import & getting the data from the computed autocorrelations (ZZ, EE, NN)

import os
if "SPHINX_DOC_BUILD" in os.environ:
    if "MSNOISE_DOC" in os.environ:
        os.chdir(os.environ["MSNOISE_DOC"])

import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("ggplot")

from msnoise.api import connect
from msnoise.results import MSNoiseResult


# connect to the database
db = connect()

# Build a result object at the stack step (filter 2 for HVSR)
result = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=2, stack=1)
params = result.params

# Get the first mov_stack configured:
mov_stack = params.mov_stack[0]

# Get the autocorrelation CCFs (ZZ, EE, NN) for station PF.FJS.00
ZZ = result.get_ccf(pair="PF.FJS.00:PF.FJS.00", components="ZZ",
                    mov_stack=mov_stack, format="dataframe")
EE = result.get_ccf(pair="PF.FJS.00:PF.FJS.00", components="EE",
                    mov_stack=mov_stack, format="dataframe")
NN = result.get_ccf(pair="PF.FJS.00:PF.FJS.00", components="NN",
                    mov_stack=mov_stack, format="dataframe")

##############################################################
# Checking the data has the same time base (index) and joining

r = "1D"
rZZ = ZZ.resample(r).mean()
rEE = EE.resample(r).mean()
rNN = NN.resample(r).mean()

merged = pd.concat({"ZZ": rZZ, "EE": rEE, "NN": rNN}, axis=0)
merged.index.names = ["channel", "date"]

# Swap the levels of the MultiIndex
result_df = merged.swaplevel("channel", "date").sort_index(level="date")

##############################################################
# Helper functions


def get_imag(d, fs):
    NFFT = 1024 * 32
    iX = np.imag(np.fft.fft(d, n=NFFT)[:NFFT // 2])
    freqs = np.fft.fftfreq(NFFT, d=1 / fs)[:NFFT // 2]
    return freqs, iX

##############################################################
# Doing the job!


HVSRs = {}
for date, group in result_df.groupby(level="date"):
    print(f"Date: {date}")
    group = group.droplevel(0)
    try:
        f, iZZ = get_imag(group.loc["ZZ"].values, 20)
        f, iEE = get_imag(group.loc["EE"].values, 20)
        f, iNN = get_imag(group.loc["NN"].values, 20)
    except Exception:
        continue
    hvsr = np.sqrt((iEE + iNN) / iZZ)
    HVSRs[date] = hvsr

HVSR = pd.DataFrame(HVSRs, index=f)
X = HVSR.copy().fillna(0.0)
X = X.T.resample(r).mean().T


##############################################################
# Plotting

plt.subplots(figsize=(8, 10))
plt.pcolormesh(X.index, X.columns, X.T, rasterized=True, vmax=5)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Date")
plt.grid()
plt.xlim(1, 10)
plt.show()


##############################################################
# Plotting & zooming around 1-4 Hz

plt.subplots(figsize=(8, 10))
plt.pcolormesh(X.index, X.columns, X.T, rasterized=True, vmax=3)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Date")
plt.grid()
plt.xlim(1, 4)
plt.show()

##############################################################
# Plotting the HVSR curve

plt.subplots(figsize=(8, 10))
hvsr = X.mean(axis=1)
plt.plot(hvsr.index, hvsr)
plt.show()


#EOF
