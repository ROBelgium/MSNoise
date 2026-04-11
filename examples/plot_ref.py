# -*- coding: utf-8 -*-
"""
Plot a Reference CCF
====================

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

# Build a result object at the refstack step
result = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1, refstack=1)

# Get the reference CCF for two stations, filter 1, ZZ component:
ref = result.get_ref(pair="PF.FJS.00:PF.FOR.00", components="ZZ",
                     )

ref.plot(figsize=(8, 8))

#EOF
