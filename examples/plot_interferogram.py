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
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
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
ccf = xr_get_ccf("PF.FJS.00", "PF.FOR.00", "ZZ", 1, mov_stack, taxis, format="xarray")

# Plot the interferogram
ccf.plot(robust=True, figsize=(10,8))


##############################################################################
# Running a simple moving window average can be done with pandas's functions:

smooth = ccf.rolling(times=5).mean()

smooth.plot(robust=True, figsize=(10,8))


#############################################################
# Plotting the sub-daily stacks and zooming on +- 20 seconds:

if "SPHINX_DOC_BUILD" in os.environ:
    if "MSNOISE_DOC" in os.environ:
        os.chdir(os.environ["MSNOISE_DOC"])

# Get the last mov_stack configured:
mov_stack = params.mov_stack[-1]

# Get the results for two station, filter id=1, ZZ component, mov_stack=1 and the results as a 2D array:
ccf = xr_get_ccf("PF.FJS.00", "PF.FOR.00", "ZZ", 1, mov_stack, taxis, format="xarray")

# Plot the interferogram and zoom in around +- 10 seconds lag
ccf.loc[:,-20:20].plot(robust=True, figsize=(10,8))


#############################################################
# Finally, stacking to a Reference function

ccf.mean(axis=0).plot(figsize=(10,8))

#EOF