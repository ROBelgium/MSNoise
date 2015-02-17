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


def main(filterid, components, mov_stack=1, ampli=1, show=True):
    db = connect()

    pairs = get_station_pairs(db, used=1)
    maxlag = float(get_config(db, 'maxlag'))
    maxlagsamples = get_maxlag_samples(db)
    t = np.linspace(-maxlag,maxlag, maxlagsamples)

    plt.figure()
    dists=[]
    for pair in pairs:
        station1, station2 = pair

        dist = get_interstation_distance(station1, station2, station1.coordinates)
        #~ print station1.sta, station2.sta, dist
        dists.append(dist)

        sta1 = "%s.%s" % (station1.net, station1.sta)
        sta2 = "%s.%s" % (station2.net, station2.sta)
        pair = "%s:%s" % (sta1, sta2)
        ref_name = pair.replace('.', '_').replace(':', '_')
        rf = os.path.join("STACKS", "%02i" %
                          filterid, "REF", components, ref_name + ".MSEED")
        ref = read(rf)[0]
        ref.normalize()
        ref = ref.data * ampli

        plt.plot(t, ref+dist,c='k')
        
    plt.ylabel("Interstation Distance in km")
    plt.xlabel("Lag Time")
    plt.title("Filter = %02i" % filterid)
    
    for velocity in [3.0, 2.0, 1.0]:
        plt.plot([0,-max(dists)/velocity,], [0, max(dists)],c='r',label='%.1f $km s^{-1}$'%velocity)
        plt.plot([0,max(dists)/velocity,], [0, max(dists)],c='r',)
    
    

    plt.xlim(-maxlag, maxlag)
    plt.legend(loc=4)
    plt.show()
            
        
                            

if __name__ == "__main__":
    main()