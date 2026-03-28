# -*- coding: utf-8 -*-
"""
Plot an interferogram
=====================

"""

import os
if "SPHINX_DOC_BUILD" in os.environ:
    if "MSNOISE_DOC" in os.environ:
        os.chdir(os.environ["MSNOISE_DOC"])

import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt

plt.style.use("ggplot")

from msnoise.core.db import connect
from msnoise.results import MSNoiseResult

# connect to the database
db = connect()

# Build a result object at the default stack step
result = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
params = result.params

# Get the first mov_stack configured:
mov_stack = params.mov_stack[0]

# Get the CCF for two stations, filter 1, ZZ component, first mov_stack:
ccf = result.get_ccf(pair="PF.FJS.00:PF.FOR.00", components="ZZ",
                     mov_stack=mov_stack, format="xarray")

# Plot the interferogram
ccf.plot(robust=True, figsize=(10, 8))


##############################################################################
# Running a simple moving window average can be done with xarray's rolling:

smooth = ccf.rolling(times=5).mean()

smooth.plot(robust=True, figsize=(10, 8))


#############################################################
# Plotting the sub-daily stacks and zooming on +- 20 seconds:

if "SPHINX_DOC_BUILD" in os.environ:
    if "MSNOISE_DOC" in os.environ:
        os.chdir(os.environ["MSNOISE_DOC"])

# Get the last mov_stack configured:
mov_stack = params.mov_stack[-1]

ccf = result.get_ccf(pair="PF.FJS.00:PF.FOR.00", components="ZZ",
                     mov_stack=mov_stack, format="xarray")

# Plot and zoom to +-20 seconds lag
ccf.loc[:, -20:20].plot(robust=True, figsize=(10, 8))


#############################################################
# Finally, stacking to a Reference function

ccf.mean(axis=0).plot(figsize=(10, 8))

#EOF
