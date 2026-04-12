"""
This plot shows the cross-correlation functions (CCF) vs time in a very similar
.. include:: /clickhelp/msnoise-cc-plot-interferogram.rst

manner as on the *ccftime* plot above, but shows an image instead of wiggles.
The parameters allow to plot the daily or the mov-stacked CCF. Filters and
components are selectable too. Passing ``--refilter`` allows to bandpass filter
CCFs before plotting .


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

from ..core.db import connect, get_logger
from ..core.config import build_plot_outfile, get_config_set_details
from ..core.stations import check_stations_uniqueness
from ..core.workflow import build_movstack_datelist
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
    mov_stack = params.stack.mov_stack[stackid_item - 1]
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
    _times = data.coords["times"].values.astype("datetime64[ms]").astype(object)
    xextent = (date2num(_times[0]), date2num(_times[-1]), -params.cc.maxlag, params.cc.maxlag)
    ax = plt.subplot(111)
    # data = stack_total
    if refilter:
        _arr = data.values.copy()
        for i in range(_arr.shape[0]):
            _arr[i] = bandpass(_arr[i], freqmin, freqmax, params.cc.cc_sampling_rate,
                               zerophase=True)
        data = data.copy(data=_arr)
    vmax = np.nanmax(data.values) * 0.9
    plt.imshow(data.values.T, extent=xextent, aspect="auto",
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
        plt.ylim(-params.cc.maxlag, params.cc.maxlag)

    title = '%s : %s, %s, Filter %d (%.2f - %.2f Hz), Stack %i (%s_%s)' % \
            (sta1, sta2, components,
             filterid, low, high, stackid, mov_stack[0], mov_stack[1])
    if refilter:
        title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
    plt.title(title)
    fig.autofmt_xdate()
    plt.tight_layout()
    if outfile:
        outfile = build_plot_outfile(
            outfile, "interferogram", result.lineage_names,
            pair=pair, components=components, mov_stack=mov_stack)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        cursor = Cursor(ax, useblit=True, color='red', linewidth=1)  # noqa: F841
        plt.show()
    else:
        plt.close()
