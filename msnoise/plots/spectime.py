"""
This plot shows the cross-correlation functions' spectrum vs time. The
parameters allow to plot the daily or the mov-stacked CCF. Filters and
components are selectable too. The ``--ampli`` argument allows to increase the
vertical scale of the CCFs. Passing ``--refilter`` allows to bandpass filter
CCFs before computing the FFT and plotting. Passing ``--startdate`` and
``--enddate`` parameters allows to specify which period of data should be
plotted. By default the plot uses dates determined in database.

.. include:: ../clickhelp/msnoise-cc-plot-spectime.rst

Example:

``msnoise cc plot spectime YA.UV05 YA.UV11`` will plot all defaults:

.. image:: ../.static/spectime.png

Zooming in the X-axis and playing with the amplitude:

``msnoise cc plot spectime YA.UV05 YA.UV11 --xlim=0.08,1.1 --ampli=10``:

.. image:: ../.static/spectime_zoom.png

And refiltering to enhance high frequency content:

``msnoise cc plot spectime YA.UV05 YA.UV11 --xlim=0.5,1.1 --ampli=10 -r0.7:1.0``:

.. image:: ../.static/spectime_refilter.png

"""

import datetime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor
import numpy as np
from obspy.signal.filter import bandpass

from ..core.db import connect, get_logger
from ..core.config import build_plot_outfile
from ..core.stations import check_stations_uniqueness
from ..core.workflow import build_movstack_datelist
from ..core.signal import prepare_abs_positive_fft
from ..results import MSNoiseResult

def main(sta1, sta2, preprocessid=1, ccid=1, filterid=1, stackid=1, stackid_item=1,
         components="ZZ", ampli=5, show=True, outfile=None,  refilter=None,
         loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.cc_plot_spectime', loglevel,
                        with_pid=True)

    db = connect()
    result = MSNoiseResult.from_ids(db, preprocess=preprocessid, cc=ccid,
                                    filter=filterid, stack=stackid)
    params = result.params
    mov_stack = params.stack.mov_stack[stackid_item - 1]
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
        stack_total = result.get_ccf(f"{sta1}:{sta2}", components, mov_stack)
    except FileNotFoundError as fullpath:
        logger.error("FILE DOES NOT EXIST: %s, exiting" % fullpath)
        return

    # Convert times coord to matplotlib date numbers for y-axis
    _times_num = mdates.date2num(stack_total.coords["times"].values.astype("datetime64[ms]").astype(object))

    if stack_total.sizes["times"] == 0:
        logger.error("No CCF found for this request")
        return

    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)

    for i, row in zip(_times_num, stack_total.values):
        line = row.copy()
        if np.all(np.isnan(line)):
            continue

        if refilter:
            line = bandpass(line, freqmin, freqmax, params.cc.cc_sampling_rate,
                            zerophase=True)

        freq, line = prepare_abs_positive_fft(line, params.cc.cc_sampling_rate)
        line /= line.max()

        ax.plot(freq, line * ampli + i, c='k', lw=1)

    low = float(params.filter.freqmin)
    high = float(params.filter.freqmax)

    ax.set_ylim(start-datetime.timedelta(days=ampli),
                end+datetime.timedelta(days=ampli))
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"][0],kwargs["xlim"][1])

    ax.set_xlabel("Frequency [Hz]")
    ax.set_xscale('log')
    ax.grid()

    title = '%s : %s, %s\n Preprocess %i - CC %i - Filter %d (%.2f - %.2f Hz) - Stack %i (%s_%s)' %\
            (sta1, sta2, components, preprocessid, ccid,
             filterid, low, high, stackid, mov_stack[0], mov_stack[1])
    if refilter:
        title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
    ax.set_title(title)

    if outfile:
        outfile = build_plot_outfile(
            outfile, "spectime", result.lineage_names,
            pair=pair, components=components, mov_stack=mov_stack)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        cursor = Cursor(ax, useblit=True, color='red', linewidth=1)  # noqa: F841
        plt.show()
    else:
        plt.close(fig)
