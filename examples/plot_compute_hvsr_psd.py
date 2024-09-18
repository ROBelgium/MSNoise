#!/usr/bin/env python
# coding: utf-8

# In[4]:


# -*- coding: utf-8 -*-
"""
Compute the H/V spectral ratio from the ratios of PSDs
======================================================

"""

import os

if "SPHINX_DOC_BUILD" in os.environ or 1:
    os.chdir(r"C:\tmp\MSNOISE_DOC")

import matplotlib
# matplotlib.use("agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import datetime

register_matplotlib_converters()

plt.style.use("ggplot")

from msnoise.api import *


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



def convert_to_velocity(df):
    df = df.resample("30Min").mean()
    df.columns = 1. / df.columns
    df = df.sort_index(axis=1)
    df = df.dropna(axis=1, how='all')

    w2f = (2.0 * np.pi * df.columns)
    # The acceleration amplitude spectrum and velocity spectral amplitude (not power)
    vamp = np.sqrt(10.0 ** (df / 10.) / w2f ** 2)
    return vamp



Z = hdf_open_store("PF.FJS.00.HHZ", mode="r").PSD
Z = convert_to_velocity(Z)
E = hdf_open_store("PF.FJS.00.HHE", mode="r").PSD
E = convert_to_velocity(E)
N = hdf_open_store("PF.FJS.00.HHN", mode="r").PSD
N = convert_to_velocity(N)



r = "6h"
rZ = Z.resample(r).mean()
rE = E.resample(r).mean()
rN = N.resample(r).mean()

merged = pd.concat({'Z': rZ, 'E': rE, 'N': rN}, axis=0)
merged.index.names = ['channel', 'date']

# Swap the levels of the MultiIndex
result = merged.swaplevel('channel', 'date').sort_index(level='date')  # .iloc[:500]
result



HVSRs = {}
for date, group in result.groupby(level='date'):
    print(f"Date: {date}")
    group = group.droplevel(0)
    try:
        iZ = group.loc["Z"].values
        iE = group.loc["E"].values
        iN = group.loc["N"].values
    except:
        continue
    hvsr = np.sqrt((iE + iN) / iZ)
    HVSRs[date] = hvsr



HVSR = pd.DataFrame(HVSRs, index=result.columns)



X = HVSR.copy().fillna(0.0)
X = X.T.resample(r).mean().T

##############################################################
# Plotting
plt.subplots(figsize=(18, 18))
plt.pcolormesh(X.index, X.columns, X.T, rasterized=True, vmax=4)
plt.colorbar()
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.show()

##############################################################
# Plotting & zooming around 0.5-4 Hz
plt.subplots(figsize=(18, 18))
plt.pcolormesh(X.index, X.columns, X.T, rasterized=True, vmax=4)
plt.colorbar()
plt.xlim(0.5, 4)
plt.grid()
plt.xlabel("Frequency (Hz)")
plt.show()

##############################################################
# Plotting the HVSR curve (truncating the low frequency here):
X.loc[0.2:20].mean(axis=1).plot()
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude")
plt.show()



