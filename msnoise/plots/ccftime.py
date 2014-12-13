"""MSNoise is ...

Usage:
~~~~~~

.. code-block:: sh

    $ msnoise plot ccftime

"""
# plot interferogram
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter, DayLocator, MonthLocator, YearLocator
from scipy.stats import scoreatpercentile
from scipy.stats.stats import nanmean

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
import numpy as np
import sys
import scipy.signal
from obspy.core import read, Stream, Trace

from ..api import *


def main(sta1, sta2, filterid, components, mov_stack=1, ampli=5, seismic=False, show=False):
    db = connect()
    components_to_compute = get_components_to_compute(db)
    maxlag = float(get_config(db,'maxlag'))
    cc_sampling_rate = float(get_config(db,'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db)
   
    plt.figure(figsize=(16,16))
    sta1 = sta1.replace('.','_')
    sta2 = sta2.replace('.','_')
    if sta2 > sta1: # alphabetical order filtering!
        pair = "%s:%s"%(sta1,sta2)
        
        print "New Data for %s-%s-%i-%i"%(pair,components,filterid, mov_stack)
        format = "matrix"
        nstack, stack_total = get_results(db,sta1,sta2,filterid,components,datelist,mov_stack, format=format)
        ax = plt.subplot(111)
        for i, line in enumerate(stack_total):
            line /= line.max()
            plt.plot(line * ampli + i, c='k')
            if seismic:
                y1 = np.ones(len(line)) * i
                y2 = line*ampli + i
                plt.fill_between(np.arange(len(line)), y1, y2, where=y2>=y1, facecolor='k', interpolate=True)

        plt.ylabel("Lag Time (s)")
        plt.axhline(0,lw=0.5,c='k')
        plt.grid()
        plt.title('%s : %s'%(sta1,sta2))
        name = '%i.%s_%s.png'%(filterid,sta1,sta2)

        #~ plt.savefig('interfero_publi.png',dpi=300)
        plt.show()
                            

if __name__ == "__main__":
    main()