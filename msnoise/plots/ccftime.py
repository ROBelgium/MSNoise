"""
This plot shows the cross-correlation functions (CCF) vs time. The parameters
allow to plot the daily or the mov-stacked CCF. Filters and components are
selectable too. The ``--ampli`` argument allows to increase the vertical scale
of the CCFs. The ``--seismic`` shows the up-going wiggles with a black-filled
background (very heavy !).


.. code-block:: sh

    msnoise plot ccftime --hel
    Usage: msnoise-script.py plot ccftime [OPTIONS] STA1 STA2

      Plots the dv/v (parses the dt/t results)

    Options:
      -f, --filterid INTEGER   Filter ID
      -c, --comp TEXT          Components (ZZ, ZR,...)
      -m, --mov_stack INTEGER  Mov Stack to read from disk
      -a, --ampli FLOAT        Amplification
      -S, --seismic            Seismic style
      -s, --show BOOLEAN       Show interactively?
      --help                   Show this message and exit.

Example:

``msnoise plot ccftime ID.KWUI ID.POSI`` will plot all defaults:

.. image:: .static/ccftime.png
"""
# plot interferogram
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter, DayLocator, MonthLocator, YearLocator
from scipy.stats import scoreatpercentile
from scipy.stats.stats import nanmean

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from matplotlib.widgets import Cursor

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
    samples = get_maxlag_samples(db)
    cc_sampling_rate = float(get_config(db,'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db)
    base = mdates.date2num(start) 
    fig = plt.figure(figsize=(16,16))
    sta1 = sta1.replace('.','_')
    sta2 = sta2.replace('.','_')
    t = np.arange(samples)/cc_sampling_rate - maxlag

    if sta2 > sta1: # alphabetical order filtering!
        pair = "%s:%s"%(sta1,sta2)
        
        print "New Data for %s-%s-%i-%i"%(pair,components,filterid, mov_stack)
        format = "matrix"
        nstack, stack_total = get_results(db,sta1,sta2,filterid,components,datelist,mov_stack, format=format)
        ax = plt.subplot(111)
        for i, line in enumerate(stack_total):
            line /= line.max()
            plt.plot(t, line * ampli + i + base , c='k')
            if seismic:
                y1 = np.ones(len(line)) * i
                y2 = line*ampli + i + base
                plt.fill_between(t, y1, y2, where=y2>=y1, facecolor='k', interpolate=True)

        plt.ylabel("Lag Time (s)")
        plt.axhline(0,lw=0.5,c='k')
        plt.grid()
        plt.title('%s : %s'%(sta1,sta2))
        name = '%i.%s_%s.png'%(filterid,sta1,sta2)
        plt.scatter(0,[start,])
        plt.ylim(start, end)
        plt.xlim(-maxlag, maxlag)
        ax.fmt_ydata = mdates.DateFormatter('%Y-%m-%d')
        cursor = Cursor(ax, useblit=True, color='red', linewidth=1.2)
        # ax.yaxis.auto_
        # ax.yaxis.set_major_locator( YearLocator() )
        # ax.yaxis.set_major_formatter(  DateFormatter('%Y-%m') )
        # ax.yaxis.set_minor_locator( MonthLocator(interval=2) )
        # ax.yaxis.set_minor_formatter(  DateFormatter('%Y-%m-%d') )
        
        
        # fig.canvas.draw()
        # labels = [item.get_text() for item in ax.get_yticklabels()]
        # for i, label in enumerate(labels):
            # if label != '':
                # labels[i] = datelist[int(label)].strftime('%Y-%m-%d')

        # ax.set_yticklabels(labels)
        
        #~ plt.savefig('interfero_publi.png',dpi=300)
        plt.show()
                            

if __name__ == "__main__":
    main()