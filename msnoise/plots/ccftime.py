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
from matplotlib.widgets import Cursor

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

from obspy.signal.filter import envelope as obspy_envelope
from obspy.signal.filter import bandpass
from ..api import (
    build_movstack_datelist,
    check_stations_uniqueness,
    connect,
    get_logger,
    get_merged_params_for_lineage,
    get_params,
    lineage_str_to_steps,
    xr_get_ccf,
)


def main(sta1, sta2, preprocess_id=1, cc_id=1, filter_id=1, stack_id=1, stack_item=1,
         components="ZZ", ampli=5, seismic=False,
         show=False, outfile=None, envelope=False, refilter=None,
         normalize=None, loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.cc_plot_ccftime', loglevel,
                        with_pid=True)
    db = connect()
    params = get_params(db)
    lineage_names = [f"preprocess_{preprocess_id}",f"cc_{cc_id}",f"filter_{filter_id}", f"stack_{stack_id}"]
    lineage_str = "/".join(lineage_names)
    steps = lineage_str_to_steps(db, lineage_str, "/")
    paralineage, lineage_names, params = get_merged_params_for_lineage(db, params, {}, steps)
    mov_stack = params.mov_stack[stack_item-1]
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
        stack_total = xr_get_ccf(params.output_folder, lineage_names,
               sta1, sta2, components, mov_stack, None)
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
            line = bandpass(line, freqmin, freqmax, params.cc_sampling_rate,
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

    low = float(params.freqmin)
    high = float(params.freqmax)

    plt.xlabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c='k')
    plt.grid()
    title = '%s : %s, %s\n Preprocess %i - CC %i - Filter %d (%.2f - %.2f Hz) - Stack %i (%s_%s)' %\
            (sta1, sta2, components, preprocess_id, cc_id,
             filter_id, low, high, stack_id, mov_stack[0], mov_stack[1])
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
        plt.xlim(-params.maxlag, params.maxlag)
    ax.fmt_ydata = mdates.DateFormatter('%Y-%m-%d')
    cursor = Cursor(ax, useblit=True, color='red', linewidth=1.2)
    plt.tight_layout()
    if outfile:
        if outfile.startswith("?"):
            pair = pair.replace(':', '-')
            # TODO outfile naming -> make it a helper based on lineage??
            outfile = outfile.replace('?', '%s-%s-f%i-m%s_%s' % (pair,
                                                              components,
                                                              filter_id,
                                                              mov_stack[0],
                                                              mov_stack[1]))
        outfile = "ccftime " + outfile
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()
