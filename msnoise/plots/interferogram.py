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
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, DateFormatter, YearLocator, DayLocator,\
    HourLocator
from matplotlib.widgets import Cursor

from obspy.signal.filter import bandpass

from ..api import (
    build_movstack_datelist,
    check_stations_uniqueness,
    connect,
    get_config,
    get_config_set_details,
    get_logger,
    get_merged_params_for_lineage,
    get_params,
    lineage_str_to_steps,
    xr_get_ccf,
)


def main(sta1, sta2, preprocess_id=1, cc_id=1, filter_id=1, stack_id=1, stack_item=1,
         components="ZZ", show=True,
         outfile=None, refilter=None, loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.cc_plot_interferogram', loglevel,
                        with_pid=True)
    db = connect()
    params = get_params(db)
    lineage_names = [f"preprocess_{preprocess_id}", f"cc_{cc_id}", f"filter_{filter_id}", f"stack_{stack_id}"]
    lineage_str = "/".join(lineage_names)
    steps = lineage_str_to_steps(db, lineage_str, "/")
    paralineage, lineage_names, params = get_merged_params_for_lineage(db, params, {}, steps)
    mov_stack = params.mov_stack[stack_item - 1]
    start, end, datelist = build_movstack_datelist(db)


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

    # print(mov_stack)
    output_folder = get_config(db, 'output_folder') or 'OUTPUT'

    logger.info("Fetching CCF data for %s-%s-%i-%s" % (pair, components, filter_id,
                                        mov_stack))


    try:
        data = xr_get_ccf(params.output_folder, lineage_names,
               sta1, sta2, components, mov_stack, None)
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

    filter_params = get_config_set_details(db, 'filter', filter_id, format='AttribDict')
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
             filter_id, low, high, stack_id, mov_stack[0], mov_stack[1])
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
                                                              filter_id,
                                                              mov_stack[0],
                                                              mov_stack[1]))
        outfile = "interferogram " + outfile
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()
