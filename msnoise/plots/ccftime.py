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
import datetime
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

from pandas.plotting import register_matplotlib_converters
from matplotlib.widgets import Cursor
register_matplotlib_converters()

from obspy.signal.filter import envelope as obspy_envelope
from obspy.signal.filter import bandpass
from ....db import connect, get_logger
from ....config import build_plot_outfile
from ....stations import check_stations_uniqueness
from ....workflow import build_movstack_datelist
from ..results import MSNoiseResult


def main(sta1, sta2, preprocessid=1, ccid=1, filterid=1, stackid=1, stackid_item=1,
         components="ZZ", ampli=5, seismic=False,
         show=False, outfile=None, envelope=False, refilter=None,
         normalize=None, loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.cc_plot_ccftime', loglevel,
                        with_pid=True)
    db = connect()
    result = MSNoiseResult.from_ids(db, preprocess=preprocessid, cc=ccid,
                                    filter=filterid, stack=stackid)
    params = result.params
    mov_stack = params.mov_stack[stackid_item-1]
    start, end, datelist = build_movstack_datelist(db)

    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    pair = "_".join([sta1, sta2])
    try:
        stack_total = result.get_ccf(f"{sta1}:{sta2}", components, mov_stack, format="dataframe")
        t = stack_total.columns.values
    except FileNotFoundError as fullpath:
        logger.error("FILE DOES NOT EXIST: %s, exiting" % fullpath)
        return

    # convert index to mdates
    stack_total.index = mdates.date2num(stack_total.index.to_pydatetime())

    if len(stack_total) == 0:
        logger.error("No CCF found for this request")
        return

    if normalize == "common":
        stack_total /= np.nanmax(stack_total)

    fig, ax = plt.subplots(1, 1,figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    for i, line in stack_total.iterrows():
        if np.all(np.isnan(line)):
            continue
        if refilter:
            line = bandpass(line, freqmin, freqmax, params.cc.cc_sampling_rate,
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

    low = float(params.filter.freqmin)
    high = float(params.filter.freqmax)

    plt.xlabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c='k')
    plt.grid()
    title = '%s : %s, %s\n Preprocess %i - CC %i - Filter %d (%.2f - %.2f Hz) - Stack %i (%s_%s)' %\
            (sta1, sta2, components, preprocessid, ccid,
             filterid, low, high, stackid, mov_stack[0], mov_stack[1])
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
        plt.xlim(-params.cc.maxlag, params.cc.maxlag)
    ax.fmt_ydata = mdates.DateFormatter('%Y-%m-%d')
    plt.tight_layout()
    if outfile:
        outfile = build_plot_outfile(
            outfile, "ccftime", result.lineage_names,
            pair=pair, components=components, mov_stack=mov_stack)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        cursor = Cursor(ax, useblit=True, color='red', linewidth=1)  # noqa: F841
        plt.show()
    else:
        plt.close()
