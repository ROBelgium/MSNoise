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
import sys
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Cursor
from obspy.signal.filter import bandpass

from ..api import (
    build_movstack_datelist,
    check_stations_uniqueness,
    connect,
    get_logger,
    get_merged_params_for_lineage,
    get_params,
    lineage_str_to_steps,
    prepare_abs_positive_fft,
    xr_get_ccf,
)

def main(sta1, sta2, preprocess_id=1, cc_id=1, filter_id=1, stack_id=1, stack_item=1,
         components="ZZ", ampli=5, show=True, outfile=None,  refilter=None,
         loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.cc_plot_spectime', loglevel,
                        with_pid=True)

    db = connect()
    params = get_params(db)
    lineage_names = [f"preprocess_{preprocess_id}", f"cc_{cc_id}", f"filter_{filter_id}", f"stack_{stack_id}"]
    lineage_str = "/".join(lineage_names)
    steps = lineage_str_to_steps(db, lineage_str, "/")
    paralineage, lineage_names, params = get_merged_params_for_lineage(db, params, {}, steps)
    mov_stack = params.mov_stack[stack_item - 1]
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

    fig, ax = plt.subplots(1, 1, figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)

    for i, line in stack_total.iterrows():
        if np.all(np.isnan(line)):
            continue

        if refilter:
            line = bandpass(line, freqmin, freqmax, params.cc_sampling_rate,
                            zerophase=True)

        freq, line = prepare_abs_positive_fft(line, params.cc_sampling_rate)
        line /= line.max()

        ax.plot(freq, line * ampli + i, c='k', lw=1)

    low = float(params.freqmin)
    high = float(params.freqmax)

    ax.set_ylim(start-datetime.timedelta(days=ampli),
                end+datetime.timedelta(days=ampli))
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"][0],kwargs["xlim"][1])

    ax.set_xlabel("Frequency [Hz]")
    ax.set_xscale('log')
    ax.grid()

    title = '%s : %s, %s\n Preprocess %i - CC %i - Filter %d (%.2f - %.2f Hz) - Stack %i (%s_%s)' %\
            (sta1, sta2, components, preprocess_id, cc_id,
             filter_id, low, high, stack_id, mov_stack[0], mov_stack[1])
    if refilter:
        title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
    ax.set_title(title)

    cursor = Cursor(ax, useblit=True, color='red', linewidth=1.2)
    if outfile:
        if outfile.startswith("?"):
            pair = pair.replace(':', '-')
            # TODO outfile naming -> make it a helper based on lineage??
            outfile = outfile.replace('?', '%s-%s-f%i-m%s_%s' % (pair,
                                                              components,
                                                              filter_id,
                                                              mov_stack[0],
                                                              mov_stack[1]))
        outfile = "spectime " + outfile
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close(fig)
