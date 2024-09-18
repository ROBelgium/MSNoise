# -*- coding: utf-8 -*-
"""
Plot a Reference CCF
====================

"""

# The following two lines are only needed for building this documentation
# Delete them and run the code in your project folder.

import os
if "SPHINX_DOC_BUILD" in os.environ or 1:
    os.chdir(r"C:\tmp\MSNOISE_DOC")


import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()


plt.style.use("ggplot")

from msnoise.api import connect, get_results, build_movstack_datelist, get_params, get_t_axis, xr_get_ref

# connect to the database
db = connect()

# Obtain a list of dates between ``start_date`` and ``enddate``
start, end, datelist = build_movstack_datelist(db)

# Get the list of parameters from the DB:
params = get_params(db)

# Get the time axis for plotting the CCF:
taxis = get_t_axis(db)

# Get the results for two station, filter id=1, ZZ component, mov_stack=("1d","1d") and stack the results:
ccf = xr_get_ref("PF.FJS.00", "PF.FOR.00", "ZZ", 1, taxis)


ccf.plot(figsize=(8,8))


#EOF