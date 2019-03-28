"""
This plot shows the cross-correlation functions (CCF) vs time in a very similar
manner as on the *ccftime* plot above, but shows an image instead of wiggles.
The parameters allow to plot the daily or the mov-stacked CCF. Filters and
components are selectable too. Passing ``--refilter`` allows to bandpass filter
CCFs before plotting (new in 1.5).

.. include:: clickhelp/msnoise-plot-interferogram.rst

Example:

``msnoise plot interferogram YA.UV06 YA.UV10 -m5`` will plot the ZZ component
(default), filter 1 (default) and mov_stack 5:

.. image:: .static/interferogram.png

"""
# plot interferogram
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter, YearLocator, DayLocator,\
    HourLocator
from matplotlib.widgets import Cursor

from obspy.signal.filter import bandpass

from ..api import *


def main(sta1, sta2, filterid, components, mov_stack=1, show=True,
         outfile=None, refilter=None, **kwargs):
    db = connect()
    maxlag = float(get_config(db, 'maxlag'))
    cc_sampling_rate = float(get_config(db, 'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db)
    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)
    fig = plt.figure(figsize=(12, 9))
    sta1 = sta1.replace('.', '_')
    sta2 = sta2.replace('.', '_')
    if sta2 >= sta1:
        pair = "%s:%s" % (sta1, sta2)
        
        print("New Data for %s-%s-%i-%i" % (pair, components, filterid,
                                            mov_stack))

        nstack, stack_total = get_results(db, sta1, sta2, filterid, components,
                                          datelist, mov_stack, format="matrix")

        xextent = (date2num(start), date2num(end), -maxlag, maxlag)
        ax = plt.subplot(111)
        data = stack_total
        if refilter:
            for i, d in enumerate(data):
                data[i] = bandpass(data[i], freqmin, freqmax, cc_sampling_rate,
                                   zerophase=True)
        vmax = np.nanmax(data) * 0.9
        plt.imshow(data.T, extent=xextent, aspect="auto",
                   interpolation='none', origin='lower', cmap='seismic',
                   vmin=-vmax, vmax=vmax)
        plt.ylabel("Lag Time (s)")
        plt.axhline(0, lw=0.5, c='k')
        plt.grid()

        # ax.xaxis.set_major_locator(DayLocator())
        # ax.xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d"))
        # ax.xaxis.set_minor_locator(DayLocator())
        # ax.xaxis.set_minor_formatter(DateFormatter('%Y-%m-%d %H:%M'))

        for filterdb in get_filters(db, all=True):
            if filterid == filterdb.ref:
                low = float(filterdb.low)
                high = float(filterdb.high)
                break

        if "ylim" in kwargs:
            plt.ylim(kwargs["ylim"][0],kwargs["ylim"][1])
        else:
            plt.ylim(-maxlag, maxlag)

        title = '%s : %s, %s, Filter %d (%.2f - %.2f Hz), Stack %d' % \
                (sta1.replace('_', '.'), sta2.replace('_', '.'), components,
                 filterid, low, high, mov_stack)
        if refilter:
            title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
        plt.title(title)
        fig.autofmt_xdate()
        cursor = Cursor(ax, useblit=True, color='black', linewidth=1.2)
        if outfile:
            if outfile.startswith("?"):
                pair = pair.replace(':', '-')
                outfile = outfile.replace('?', '%s-%s-f%i-m%i' % (pair,
                                                                  components,
                                                                  filterid,
                                                                  mov_stack))
            outfile = "interferogram " + outfile
            print("output to:", outfile)
            plt.savefig(outfile)
        if show:
            plt.show()
