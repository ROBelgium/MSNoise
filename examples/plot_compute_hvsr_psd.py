# -*- coding: utf-8 -*-
"""
Compute the H/V spectral ratio from the ratios of PSDs
======================================================

"""

import os

if "SPHINX_DOC_BUILD" in os.environ:
    if "MSNOISE_DOC" in os.environ:
        os.chdir(os.environ["MSNOISE_DOC"])

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.style.use("ggplot")

from msnoise.core.db import connect
from msnoise.results import MSNoiseResult

import matplotlib
matplotlib.use("agg")

# connect to the database
db = connect()

# Build a result object at the psd step (psd=1 only; no preprocess needed for PSD)
result = MSNoiseResult.from_ids(db, psd=1)

def convert_to_velocity(df):
    df = df.resample("30Min").mean()
    df.columns = 1.0 / df.columns
    df = df.sort_index(axis=1)
    df = df.dropna(axis=1, how="all")
    w2f = 2.0 * np.pi * df.columns
    # Velocity spectral amplitude (not power)
    vamp = np.sqrt(10.0 ** (df / 10.0) / w2f ** 2)
    return vamp

# Load PSDs via MSNoiseResult.get_psd — returns xarray Dataset directly
Z_ds = result.get_psd(seed_id="PF.FJS.00.HHZ")
E_ds = result.get_psd(seed_id="PF.FJS.00.HHE")
N_ds = result.get_psd(seed_id="PF.FJS.00.HHN")

Z = convert_to_velocity(Z_ds.PSD.to_dataframe().unstack())
E = convert_to_velocity(E_ds.PSD.to_dataframe().unstack())
N = convert_to_velocity(N_ds.PSD.to_dataframe().unstack())

r = "6h"
rZ = Z.resample(r).mean()
rE = E.resample(r).mean()
rN = N.resample(r).mean()

merged = pd.concat({"Z": rZ, "E": rE, "N": rN}, axis=0)
merged.index.names = ["channel", "date"]

result_df = merged.swaplevel("channel", "date").sort_index(level="date")

HVSRs = {}
for date, group in result_df.groupby(level="date"):
    print(f"Date: {date}")
    group = group.droplevel(0)
    try:
        iZ = group.loc["Z"].values
        iE = group.loc["E"].values
        iN = group.loc["N"].values
    except Exception:
        continue
    hvsr = np.sqrt((iE + iN) / iZ)
    HVSRs[date] = hvsr

HVSR = pd.DataFrame(HVSRs, index=result_df.columns)
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
