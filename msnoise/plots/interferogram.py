"""
This plot shows the cross-correlation functions (CCF) vs time in a very similar
manner as on the *ccftime* plot above, but shows an image instead of wiggles.
The parameters allow to plot the daily or the mov-stacked CCF. Filters and
components are selectable too.

.. include:: clickhelp/msnoise-plot-interferogram.rst

Example:

``msnoise plot interferogram YA.UV06 YA.UV10 -m5`` will plot the ZZ component
(default), filter 1 (default) and mov_stack 5:

.. image:: .static/interferogram.png

"""
# plot interferogram
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter, YearLocator
from matplotlib.widgets import Cursor

from ..api import *

def main(sta1, sta2, filterid, components, mov_stack=1, show=True, outfile=None):
    db = connect()
    components_to_compute = get_components_to_compute(db)
    maxlag = float(get_config(db,'maxlag'))
    cc_sampling_rate = float(get_config(db,'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db)
    # mov_stack = get_config(db,"mov_stack")
 
   
    fig = plt.figure(figsize=(16,16))
    sta1 = sta1.replace('.','_')
    sta2 = sta2.replace('.','_')
    if sta2 >= sta1: # alphabetical order filtering!
        pair = "%s:%s"%(sta1,sta2)
        
        print("New Data for %s-%s-%i-%i"%(pair,components,filterid, mov_stack))
        format = "matrix"
        nstack, stack_total = get_results(db,sta1,sta2,filterid,components,
                                          datelist,mov_stack, format=format)

        xextent = (date2num(start), date2num(end),-maxlag,maxlag)
        ax = plt.subplot(111)
        plt.imshow(stack_total.T, extent=xextent, aspect="auto",
                   interpolation='none',origin='lower',cmap='seismic',
                   vmin=-1e-2,vmax=1e-2)
        plt.ylabel("Lag Time (s)")
        plt.axhline(0,lw=0.5,c='k')
        plt.grid()

        ax.xaxis.set_major_locator( YearLocator() )
        ax.xaxis.set_major_formatter(  DateFormatter('%Y-%m') )

        for filterdb in get_filters(db, all=True):
            if filterid == filterdb.ref:
                low = float(filterdb.low)
                high = float(filterdb.high)
                break
        
        lag = 120
        plt.ylim(-lag,lag)
        plt.title('%s : %s, %s, Filter %d (%.2f - %.2f Hz), Stack %d' %
                  (sta1.replace('_', '.'), sta2.replace('_', '.'), components,
                   filterid, low, high, mov_stack))
        cursor = Cursor(ax, useblit=True, color='black', linewidth=1.2)
        if outfile:
            if outfile.startswith("?"):
                pair = pair.replace(':','-')
                outfile = outfile.replace('?', '%s-%s-f%i-m%i' % (pair,
                                                                  components,
                                                                  filterid,
                                                                  mov_stack))
            outfile = "interferogram " + outfile
            print("output to:", outfile)
            plt.savefig(outfile)
        if show:
            plt.show()


if __name__ == "__main__":
    main()
