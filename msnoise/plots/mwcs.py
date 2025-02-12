"""
This plot shows the result of the MWCS calculations in two superposed images.
One is the dt calculated vs time lag and the other one is the coherence. The
image is constructed by horizontally stacking the MWCS of different days. The
two right panels show the mean and standard deviation per time lag of the whole
image. The selected time lags for the dt/t calculation are presented with green
horizontal lines, and the minimum coherence or the maximum dt are in red.

The ``filterid``, ``comp`` and ``mov_stack`` allow filtering the data used.

.. include:: ../clickhelp/msnoise-cc-dvv-plot-mwcs.rst

Example:

``msnoise cc dvv plot mwcs ID.KWUI ID.POSI -m 3`` will plot all defaults with the
mov_stack = 3:

.. image:: ../.static/mwcs.png

"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from matplotlib.dates import date2num, AutoDateFormatter, AutoDateLocator

from ..api import *


def main(sta1, sta2, filterid, components, mov_stack=1, show=True,
         outfile=None, loglevel='INFO'):
    logger = get_logger('msnoise.cc_dvv_plot_mwcs', loglevel,
                        with_pid=True)
    db = connect()
    maxlag = float(get_config(db, 'maxlag'))
    start, end, datelist = build_movstack_datelist(db)
    
    dtt_lag = get_config(db, "mwcs_dtt_lag")
    dtt_v = float(get_config(db, "dtt_v"))
    dtt_minlag = float(get_config(db, "dtt_minlag"))
    dtt_width = float(get_config(db, "dtt_width"))
    dtt_sides = get_config(db, "mwcs_dtt_sides")
    minCoh = float(get_config(db, "mwcs_dtt_mincoh"))
    maxErr = float(get_config(db, "mwcs_dtt_maxerr"))
    maxDt = float(get_config(db, "mwcs_dtt_maxdt"))

    def plot_lags(minlag, maxlag):
        plt.axhline(minlag, c='g')
        plt.axhline(-minlag, c='g')
        plt.axhline(maxlag, c='g')
        plt.axhline(-maxlag, c='g')

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)
    pair = "%s_%s" % (sta1, sta2)

    station1 = sta1.split(".")
    station2 = sta2.split(".")

    station1 = get_station(db, station1[0], station1[1])
    station2 = get_station(db, station2[0], station2[1])

    if dtt_lag == "static":
        minlag = dtt_minlag
    else:
        minlag = get_interstation_distance(station1, station2,
                                           station1.coordinates) / dtt_v


    maxlag2 = minlag + dtt_width
    params = get_params(db)
    mov_stack = params.mov_stack[mov_stack-1]
    logger.info("Fetching CCF data for %s-%s-%i-%s" % (pair, components, filterid,
                                                       mov_stack))

    id = []
    alldt = []
    allcoh = []
    # for day in datelist:
    #     fname = os.path.join('MWCS', "%02i" % filterid, "%03i_DAYS" %
    #                          mov_stack, components, pair, '%s.txt' % day)
    #     if os.path.isfile(fname):
    #         df = pd.read_csv(fname, delimiter=' ', header=None, index_col=0,
    #                          names=['t', 'dt', 'err', 'coh'])
    #         alldt.append(df["dt"])
    #         allcoh.append(df["coh"])
    #         id.append(day)
    #         del df

    try:
        mwcs = xr_get_mwcs(sta1, sta2, components, filterid, mov_stack)
    except FileNotFoundError as fullpath:
        logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
        return

    alldt = mwcs.M
    allcoh = mwcs.MCOH
    id = mwcs.index
    # print(mwcs.M)

    alldt = alldt.resample('D').mean()
    allcoh = allcoh.resample('D').mean()

    xextent = (date2num(id[0]), date2num(id[-1]), -maxlag, maxlag)

    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 1])

    plt.figure()
    ax1 = plt.subplot(gs[0])
    im = plt.imshow(alldt.T, extent=xextent, aspect="auto",
                    interpolation='none', origin='lower', cmap=cm.seismic)
    cscale = np.nanpercentile(alldt, q=99)
    im.set_clim((cscale * -1, cscale))
    cb = plt.colorbar()
    cb.set_label('dt')
    plt.ylabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c='k')
    plt.grid()
    plt.title('%s : %s : dt' % (sta1, sta2))
    plot_lags(minlag, maxlag2)
    plt.setp(ax1.get_xticklabels(), visible=False)

    plt.subplot(gs[1], sharey=ax1)
    plt.plot(alldt.mean(axis=0), alldt.columns, c='k')
    plt.grid()
    plot_lags(minlag, maxlag2)
    plt.axvline(-maxDt, c='r', ls='--')
    plt.axvline(maxDt, c='r', ls='--')
    plt.xlabel('dt')
    plt.ylabel("Lag Time (s)")

    ax2 = plt.subplot(gs[2], sharex=ax1, sharey=ax1)
    plt.imshow(allcoh.T, extent=xextent, aspect="auto",
               interpolation='none', origin='lower', cmap='hot',
               vmin=minCoh, vmax=1)

    cb = plt.colorbar()
    cb.set_label('mean coherence')
    plt.ylabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c='k')
    plt.grid()
    locator = AutoDateLocator()
    ax2.xaxis.set_major_locator(locator)
    ax2.xaxis.set_major_formatter(AutoDateFormatter(locator))
    plt.setp(plt.xticks()[1], rotation=30, ha='right')
    plt.title('%s : %s : mean coherence' % (sta1, sta2))
    plot_lags(minlag, maxlag2)

    plt.subplot(gs[3], sharey=ax1)
    m = allcoh.mean(axis=0)
    s = allcoh.std(axis=0)
    plt.plot(m, allcoh.columns, c='k')
    plt.fill_betweenx(allcoh.columns, m-s, m+s, color='silver',)

    plt.grid()
    plot_lags(minlag, maxlag2)
    plt.axvline(minCoh, c='r', ls='--')
    plt.xlabel('Coherence')
    plt.ylabel("Lag Time (s)")

    name = '%s-%s f%i m%s' % (sta1, sta2, filterid, mov_stack)
    name = name.replace('_', '.')

    plt.suptitle(name)

    if outfile:
        if outfile.startswith("?"):
            pair = pair.replace(':', '-')
            outfile = outfile.replace('?', '%s-%s-f%i-m%s_%s' % (pair,
                                                              components,
                                                              filterid,
                                                              mov_stack[0],
                                                              mov_stack[1]))
        outfile = "mwcs " + outfile
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
