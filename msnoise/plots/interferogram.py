"""MSNoise is ...

Usage:
~~~~~~

.. code-block:: sh

    $ msnoise plot interferogram

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


def main(sta1, sta2, filterid, components, mov_stack=1):
    db = connect()
    components_to_compute = get_components_to_compute(db)
    maxlag = float(get_config(db,'maxlag'))
    cc_sampling_rate = float(get_config(db,'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db)
    # mov_stack = get_config(db,"mov_stack")
 
    ## if end date later than today, adjust end and datelist
    end = min(end, datetime.date.today())
    after_today = [days for days in range(len(datelist)) if datelist[days] > end]
    if after_today:
        del datelist[after_today[0]:]
    ##

   
    plt.figure(figsize=(16,16))
    sta1 = sta1.replace('.','_')
    sta2 = sta2.replace('.','_')
    if sta2 > sta1: # alphabetical order filtering!
        pair = "%s:%s"%(sta1,sta2)
        
        print "New Data for %s-%s-%i-%i"%(pair,components,filterid, mov_stack)
        format = "matrix"
        nstack, stack_total = get_results(db,sta1,sta2,filterid,components,datelist,mov_stack, format=format)
        # vmax = scoreatpercentile(np.abs(stack_total[np.isnan(stack_total)==False]) , 98)
        # for i in range(stack_total.shape[0]):
            # if not np.all( np.isnan(stack_total[i,:])):
                # print np.max(stack_total[i,:])
                # stack_total[i,:] /= np.max(stack_total[i,:])
        # stack_total /= np.max(stack_total, axis=0)
        xextent = (date2num(start), date2num(end),-120,120)
        ax = plt.subplot(111)
        plt.imshow(stack_total.T, extent=xextent, aspect="auto",interpolation='none',origin='lower',cmap='seismic',
                vmin=-1e-1,vmax=1e-1)
        plt.ylabel("Lag Time (s)")
        plt.axhline(0,lw=0.5,c='k')
        plt.grid()

        ax.xaxis.set_major_locator( YearLocator() )
        ax.xaxis.set_major_formatter(  DateFormatter('%Y-%m') )

        # ax.xaxis.set_minor_locator( MonthLocator(interval=2) )
        # ax.xaxis.set_minor_formatter(  DateFormatter('%Y-%m-%d') )
        
        lag = 120
        plt.ylim(-lag,lag)
        plt.title('%s : %s'%(sta1,sta2))
        name = '%i.%s_%s.png'%(filterid,sta1,sta2)

        #~ plt.savefig('interfero_publi.png',dpi=300)
        # plt.figure()
        # maxx = np.argmax(stack_total, axis=0)
        # plt.plot(maxx)
        
        
        plt.show()
        
        
                            

if __name__ == "__main__":
    main()
