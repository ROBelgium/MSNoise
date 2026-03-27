"""
This plot shows the cross-correlation functions (CCF) vs time in a very similar
manner as on the *ccftime* plot above, but shows an image instead of wiggles.
The parameters allow to plot the daily or the mov-stacked CCF. Filters and
components are selectable too. Passing ``--refilter`` allows to bandpass filter
CCFs before plotting (new in 1.5).

.. include:: ../clickhelp/msnoise-cc-plot-interferogram.rst

Example:

``msnoise cc plot interferogram YA.UV06 YA.UV11 -m5`` will plot the ZZ component
(default), filter 1 (default) and mov_stack 5:

.. image:: ../.static/interferogram.png

"""
# plot interferogram
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter
from matplotlib.widgets import Cursor

from obspy.signal.filter import bandpass

from ..api import (
    build_movstack_datelist,
    check_stations_uniqueness,
    connect,
    get_config_set_details,
    get_logger,
)
from ..results import MSNoiseResult


def main(sta1, sta2, preprocessid=1, ccid=1, filterid=1, stackid=1, stackid_item=1,
         components="ZZ", show=True,
         outfile=None, refilter=None, loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.cc_plot_interferogram', loglevel,
                        with_pid=True)
    db = connect()
    result = MSNoiseResult.from_ids(db, preprocess=preprocessid, cc=ccid,
                                    filter=filterid, stack=stackid)
    params = result.params
    mov_stack = params.mov_stack[stackid_item - 1]
    start, end, datelist = build_movstack_datelist(db)

    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)
    fig = plt.figure(figsize=(12, 9))

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    pair = "%s:%s" % (sta1, sta2)

    logger.info("Fetching CCF data for %s-%s-%i-%s" % (pair, components, filterid,
                                        mov_stack))


    try:
        data = result.get_ccf(f"{sta1}:{sta2}", components, mov_stack)
    except FileNotFoundError as fullpath:
        logger.error("FILE DOES NOT EXIST: %s, exiting" % fullpath)
        return
    xextent = (date2num(data.index[0]), date2num(data.index[-1]), -params.maxlag, params.maxlag)
    ax = plt.subplot(111)
    # data = stack_total
    if refilter:
        for i, d in enumerate(data):
            data.iloc[i] = bandpass(data.iloc[i], freqmin, freqmax, params.cc_sampling_rate,
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

    filter_params = get_config_set_details(db, 'filter', filterid, format='AttribDict')
    if filter_params:
        low = float(filter_params.freqmin)
        high = float(filter_params.freqmax)
    else:
        low = high = 0.0

    if "ylim" in kwargs:
        plt.ylim(kwargs["ylim"][0],kwargs["ylim"][1])
    else:
        plt.ylim(-params.maxlag, params.maxlag)

    title = '%s : %s, %s, Filter %d (%.2f - %.2f Hz), Stack %i (%s_%s)' % \
            (sta1, sta2, components,
             filterid, low, high, stackid, mov_stack[0], mov_stack[1])
    if refilter:
        title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
    plt.title(title)
    fig.autofmt_xdate()
    plt.tight_layout()
    cursor = Cursor(ax, useblit=True, color='black', linewidth=1.2)
    if outfile:
        if outfile.startswith("?"):
            pair = pair.replace(':', '-')
            outfile = outfile.replace('?', '%s-%s-f%i-m%s_%s' % (pair,
                                                              components,
                                                              filterid,
                                                              mov_stack[0],
                                                              mov_stack[1]))
        outfile = "interferogram " + outfile
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()
