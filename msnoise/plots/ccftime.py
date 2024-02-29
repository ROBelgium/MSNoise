"""
This plot shows the cross-correlation functions (CCF) vs time. The parameters
allow to plot the daily or the mov-stacked CCF. Filters and components are
selectable too. The ``--ampli`` argument allows to increase the vertical scale
of the CCFs. The ``--seismic`` shows the up-going wiggles with a black-filled
background (very heavy !). Passing ``--refilter`` allows to bandpass filter
CCFs before plotting (new in 1.5).

.. include:: ../clickhelp/msnoise-cc-plot-ccftime.rst


Example:

``msnoise cc plot ccftime YA.UV06 YA.UV11`` will plot all defaults:

.. image:: ../.static/ccftime.png

For zooming in the CCFs:

``msnoise cc plot ccftime YA.UV05 YA.UV11 --xlim=-10,10 --ampli=15``:

.. image:: ../.static/ccftime_zoom.png


It is sometimes useful to refilter the CCFs on the fly:

``msnoise cc plot ccftime YA.UV05 YA.UV11 -r 0.5:1.0``:

.. image:: ../.static/ccftime_refilter.png



"""
# plot interferogram

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

from obspy.signal.filter import envelope as obspy_envelope
from obspy.signal.filter import bandpass
from ..api import *


def main(sta1, sta2, filterid, components, mov_stackid=1, ampli=5, seismic=False,
         show=False, outfile=None, envelope=False, refilter=None,
         normalize=None, loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.cc_plot_ccftime', loglevel,
                        with_pid=True)
    db = connect()
    params = get_params(db)
    mov_stack = params.mov_stack[mov_stackid-1]
    maxlag = float(get_config(db, 'maxlag'))
    samples = get_maxlag_samples(db)
    cc_sampling_rate = float(get_config(db, 'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db)
    base = mdates.date2num(start) 
    plt.figure(figsize=(12, 9))
    sta1 = sta1 #.replace('.', '_')
    sta2 = sta2 #.replace('.', '_')
    t = np.arange(samples)/cc_sampling_rate - maxlag
    taxis = get_t_axis(db)

    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    pair = "%s:%s" % (sta1, sta2)

    logger.info("Fetching CCF data for %s-%s-%i-%s" % (pair, components, filterid,
                                                 str(mov_stack)))
    try:
        stack_total = xr_get_ccf(sta1, sta2, components, filterid, mov_stack, taxis)
    except FileNotFoundError as fullpath:
        logger.error("FILE DOES NOT EXIST: %s, exiting" % fullpath)
        sys.exit(1)

    # convert index to mdates
    stack_total.index = mdates.date2num(stack_total.index.to_pydatetime())

    if len(stack_total) == 0:
        logger.error("No CCF found for this request")
        return

    if normalize == "common":
        stack_total /= np.nanmax(stack_total)
    ax = plt.subplot(111)
    for i, line in stack_total.iterrows():
        if np.all(np.isnan(line)):
            continue
        if refilter:
            line = bandpass(line, freqmin, freqmax, cc_sampling_rate,
                            zerophase=True)
        if envelope:
            line = obspy_envelope(line)
        if normalize == "individual":
            line /= line.max()
        plt.plot(t, line * ampli + i, c='k', lw=0.5)
        if seismic:
            y1 = np.ones(len(line)) * i
            y2 = line*ampli + i
            plt.fill_between(t, y1, y2, where=y2 >= y1, facecolor='k',
                             interpolate=True)

    filter = get_filters(db, ref=filterid)
    low = float(filter.low)
    high = float(filter.high)

    plt.xlabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c='k')
    plt.grid()
    title = '%s : %s, %s, Filter %d (%.2f - %.2f Hz), Stack %i (%s_%s)' %\
            (sta1, sta2, components,
             filterid, low, high, mov_stackid, mov_stack[0], mov_stack[1])
    if refilter:
        title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
    plt.title(title)
    plt.scatter(0, [start, ], alpha=0)
    plt.xlabel("Time Lag (s)")
    plt.ylim(start-datetime.timedelta(days=10),
             end+datetime.timedelta(days=10))
    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"][0],kwargs["xlim"][1])
    else:
        plt.xlim(-maxlag, maxlag)
    ax.fmt_ydata = mdates.DateFormatter('%Y-%m-%d')
    cursor = Cursor(ax, useblit=True, color='red', linewidth=1.2)
    plt.tight_layout()
    if outfile:
        if outfile.startswith("?"):
            pair = pair.replace(':', '-')
            outfile = outfile.replace('?', '%s-%s-f%i-m%s_%s' % (pair,
                                                              components,
                                                              filterid,
                                                              mov_stack[0],
                                                              mov_stack[1]))
        outfile = "ccftime " + outfile
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
