"""
This plot shows the result of the MWCS calculations in two superposed images.
One is the dt calculated vs time lag and the other one is the coherence. The
image is constructed by horizontally stacking the MWCS of different days. The
two right panels show the mean and standard deviation per time lag of the whole
image. The selected time lags for the dt/t calculation are presented with green
horizontal lines, and the minimum coherence or the maximum dt are in red.

The ``filterid``, ``comp`` and ``mov_stack`` allow filtering the data used.

.. include:: clickhelp/msnoise-plot-mwcs.rst

Example:

``msnoise plot mwcs ID.KWUI ID.POSI -m 3`` will plot all defaults with the
mov_stack = 3:

.. image:: .static/mwcs.png

"""

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

from matplotlib import pyplot, colors, cm

def cmap_powerlaw_adjust(cmap, a):
    '''
    returns a new colormap based on the one given
    but adjusted via power-law:

    newcmap = oldcmap**a
    '''
    if a < 0.:
        return cmap
    cdict = copy.copy(cmap._segmentdata)
    fn = lambda x : (x[0]**a, x[1], x[2])
    for key in ('red','green','blue'):
        cdict[key] = map(fn, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1), \
            "Resulting indices extend out of the [0, 1] segment."
    return colors.LinearSegmentedColormap('colormap',cdict,1024)

def cmap_center_adjust(cmap, center_ratio):
    '''
    returns a new colormap based on the one given
    but adjusted so that the old center point higher
    (>0.5) or lower (<0.5)
    '''
    if not (0. < center_ratio) & (center_ratio < 1.):
        return cmap
    a = math.log(center_ratio) / math.log(0.5)
    return cmap_powerlaw_adjust(cmap, a)

def cmap_center_point_adjust(cmap, range, center):
    '''
    converts center to a ratio between 0 and 1 of the
    range given and calls cmap_center_adjust(). returns
    a new adjusted colormap accordingly
    '''
    if not ((range[0] < center) and (center < range[1])):
        print("beu")
        return cmap
    return cmap_center_adjust(cmap,
        abs(center - range[0]) / abs(range[1] - range[0]))

def main(sta1, sta2, filterid, components, mov_stack=1, show=True, outfile=None):
    db = connect()
    maxlag = float(get_config(db,'maxlag'))
    start, end, datelist = build_movstack_datelist(db)
    
    dtt_lag = get_config(db, "dtt_lag")
    dtt_v = float(get_config(db, "dtt_v"))
    dtt_minlag = float(get_config(db, "dtt_minlag"))
    dtt_width = float(get_config(db, "dtt_width"))
    dtt_sides = get_config(db, "dtt_sides")
    minCoh = float(get_config(db, "dtt_mincoh"))
    maxErr = float(get_config(db, "dtt_maxerr"))
    maxDt = float(get_config(db, "dtt_maxdt"))

    def plot_lags(minlag,maxlag):
        plt.axhline(minlag, c='g')
        plt.axhline(-minlag, c='g')
        plt.axhline(maxlag, c='g')
        plt.axhline(-maxlag, c='g')
    
    sta1 = sta1.replace('.','_')
    sta2 = sta2.replace('.','_')
    if sta2 > sta1: # alphabetical order filtering!
        pair = "%s_%s"%(sta1,sta2)
        station1 = sta1.split("_")
        station2 = sta2.split("_")
        
        station1 = get_station(db, station1[0], station1[1])
        station2 = get_station(db, station2[0], station2[1])
        
        
        if dtt_lag == "static":
            minlag = dtt_minlag
        else:
            minlag = get_interstation_distance(station1, station2, station1.coordinates) / dtt_v
            print( minlag)
        
        maxlag2 = minlag + dtt_width
        
       
        print("New Data for %s-%s-%i-%i"%(pair,components,filterid, mov_stack))
        format = "matrix"
        
        alldf=[]
        id = []
        alldt = []
        allcoh = []
        for day in datelist:
            fname = os.path.join('MWCS', "%02i" % filterid, "%03i_DAYS" %
                                           mov_stack, components, pair, '%s.txt' % day)
            #~ print fname
            if os.path.isfile(fname):
                df = pd.read_csv(fname, delimiter=' ', header=None, index_col=0, names=['t', 'dt', 'err', 'coh']) 
                alldt.append( df["dt"])
                allcoh.append( df["coh"])
                id.append(day)
                del df
        print( len(alldt[0]))
                
        alldt = pd.DataFrame(alldt,index=pd.DatetimeIndex(id))
        allcoh = pd.DataFrame(allcoh,index=pd.DatetimeIndex(id))
        
        crange = [np.amin(alldt.values),np.amax(alldt.values)]
        
        alldt = alldt.resample('D', how='mean')
        allcoh = allcoh.resample('D', how='mean')
        
        xextent = (date2num(id[0]), date2num(id[-1]), -maxlag, maxlag)
        
        
        gs = gridspec.GridSpec(2, 2,
                       width_ratios=[3,1],
                       height_ratios=[1,1]
                       )
        
        plt.figure()
        ax = plt.subplot(gs[0])
        plt.imshow(alldt.T, extent=xextent,aspect="auto",interpolation='none',origin='lower',cmap=cmap_center_point_adjust(cm.seismic,crange,0))
        cb = plt.colorbar()
        cb.set_label('dt')
        plt.ylabel("Lag Time (s)")
        plt.axhline(0,lw=0.5,c='k')
        plt.grid()
        ax.xaxis.set_major_locator( YearLocator() )
        ax.xaxis.set_major_formatter(  DateFormatter('%Y-%m-%d') )
        plt.title('%s : %s : dt'%(sta1,sta2))
        plot_lags(minlag, maxlag2)
        
        plt.subplot(gs[1])
        plt.plot(alldt.mean( axis=0), alldt.columns, c='k')
        plt.grid()
        plot_lags(minlag, maxlag2)
        plt.axvline(-maxDt, c='r', ls='--')
        plt.axvline(maxDt, c='r', ls='--')
        plt.xlabel('dt')
        plt.ylabel("Lag Time (s)")
        
        ax = plt.subplot(gs[2], sharex=ax, sharey=ax)
        plt.imshow(allcoh.T, extent=xextent,aspect="auto",interpolation='none',origin='lower',cmap='hot',vmin=minCoh, vmax=1)
        #~ plt.imshow(allcoh.T, extent=xextent,aspect="auto",interpolation='none',origin='lower',cmap='hot')
        cb = plt.colorbar()
        cb.set_label('mean coherence')
        plt.ylabel("Lag Time (s)")
        plt.axhline(0,lw=0.5,c='k')
        plt.grid()
        ax.xaxis.set_major_locator( YearLocator() )
        ax.xaxis.set_major_formatter(  DateFormatter('%Y-%m-%d') )
        plt.title('%s : %s : mean coherence'%(sta1,sta2))
        plot_lags(minlag, maxlag2)
        
        plt.subplot(gs[3])
        m = allcoh.mean( axis=0)
        s = allcoh.std(axis=0)
        plt.plot(m, allcoh.columns, c='k')
        plt.fill_betweenx( allcoh.columns, m-s, m+s, color='silver',)

        plt.grid()
        plot_lags(minlag, maxlag2)
        plt.axvline(minCoh, c='r', ls='--')
        plt.xlabel('Coherence')
        plt.ylabel("Lag Time (s)")
        
        name = '%s-%s f%i m%i' % (sta1,sta2, filterid, mov_stack)
        name = name.replace('_', '.')
        
        plt.suptitle(name)

        if outfile:
            if outfile.startswith("?"):
                pair = pair.replace(':','-')
                outfile = outfile.replace('?', '%s-%s-f%i-m%i' % (pair,
                                                                  components,
                                                                  filterid,
                                                                  mov_stack))
            outfile = "mwcs " + outfile
            print("output to:", outfile)
            plt.savefig(outfile)
        if show:
            plt.show()
        
        
                            

if __name__ == "__main__":
    main()