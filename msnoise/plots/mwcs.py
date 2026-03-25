"""
This plot shows the result of the MWCS calculations in two superposed images.
One is the dt calculated vs time lag and the other one is the coherence. The
image is constructed by horizontally stacking the MWCS of different days. The
two right panels show the mean and standard deviation per time lag of the whole
image. The selected time lags for the dt/t calculation are presented with green
horizontal lines, and the minimum coherence or the maximum dt are in red.

.. include:: ../clickhelp/msnoise-cc-dtt-plot-mwcs.rst

Example:

``msnoise cc dtt plot mwcs BE.UCC.-- BE.MEM.-- -f 1 -w 1 -mi 1`` will plot
MWCS for filter 1, mwcs set 1, first mov_stack:

.. image:: ../.static/mwcs.png
"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.dates import date2num, AutoDateFormatter, AutoDateLocator
import numpy as np

from ..api import (
    connect, get_params, get_logger, build_movstack_datelist,
    get_station, get_interstation_distance, check_stations_uniqueness,
    get_config_set_details, get_done_lineages_for_category,
    xr_get_mwcs, get_merged_params_for_lineage, resolve_lineage_from_ids,
)

def main(sta1, sta2, preprocess_id=1, cc_id=1, filter_id=1, stack_id=1,
         stack_item=1, mwcs_id=1,
         components="ZZ", show=True, outfile=None, loglevel="INFO"):
    """Plot MWCS dt and coherence images for a station pair.

    :param sta1: Station 1 in NET.STA.LOC format.
    :param sta2: Station 2 in NET.STA.LOC format (must be ≥ sta1 alphabetically).
    :param preprocess_id: Preprocessing step set number.
    :param cc_id: Cross-correlation step set number.
    :param filter_id: Filter step set number.
    :param stack_id: Stack step set number.
    :param stack_item: 1-based index into ``params.mov_stack``.
    :param mwcs_id: MWCS step set number.
    :param components: Component pair string (e.g. ``'ZZ'``).
    :param show: Display the figure interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_mwcs", loglevel, with_pid=True)
    db = connect()
    params = get_params(db)

    lineage_names, params = resolve_lineage_from_ids(db, params, preprocess_id=preprocess_id, cc_id=cc_id, filter_id=filter_id, stack_id=stack_id, mwcs_id=mwcs_id)
    mov_stack = params.mov_stack[stack_item - 1]

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    # DTT lag window parameters (from merged params or defaults)
    dtt_lag   = getattr(params, "dtt_lag",    "static")
    dtt_v     = float(getattr(params, "dtt_v",     1.0))
    dtt_minlag= float(getattr(params, "dtt_minlag", 5.0))
    dtt_width = float(getattr(params, "dtt_width",  30.0))
    minCoh    = float(getattr(params, "dtt_mincoh", 0.0))
    maxErr    = float(getattr(params, "dtt_maxerr", 0.0))
    maxDt     = float(getattr(params, "dtt_maxdtt", 1.0))
    maxlag    = float(getattr(params, "maxlag",     120.0))

    if dtt_lag == "dynamic":
        st1 = get_station(db, *sta1.split(".")[:2])
        st2 = get_station(db, *sta2.split(".")[:2])
        minlag = get_interstation_distance(st1, st2, st1.coordinates) / dtt_v
    else:
        minlag = dtt_minlag
    maxlag2 = minlag + dtt_width

    def plot_lags(mn, mx):
        for v in (mn, -mn, mx, -mx):
            plt.axhline(v, c="g", lw=0.8)

    logger.info("Loading MWCS for %s-%s comp=%s mov_stack=%s lineage=%s",
                sta1, sta2, components, mov_stack, "/".join(lineage_names))

    try:
        mwcs = xr_get_mwcs(params.output_folder, lineage_names,
                             sta1, sta2, components, mov_stack)
    except FileNotFoundError as fp:
        logger.error("FILE DOES NOT EXIST: %s", fp)
        return

    alldt  = mwcs["M"].resample("D").mean()
    allcoh = mwcs["MCOH"].resample("D").mean()

    xextent = (date2num(alldt.index[0]), date2num(alldt.index[-1]),
               -maxlag, maxlag)

    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 1])
    plt.figure(figsize=(12, 9))

    ax1 = plt.subplot(gs[0])
    im = plt.imshow(alldt.T, extent=xextent, aspect="auto",
                    interpolation="none", origin="lower", cmap=cm.seismic)
    cscale = np.nanpercentile(alldt.values, 99)
    im.set_clim(-cscale, cscale)
    cb = plt.colorbar()
    cb.set_label("dt")
    plt.ylabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c="k")
    plt.grid()
    title = ("%s : %s  |  comp=%s  |  "
             "Preprocess %i – CC %i – Filter %i – Stack %i – MWCS %i  |  "
             "mov_stack=%s_%s") % (
        sta1, sta2, components,
        preprocess_id, cc_id, filter_id, stack_id, mwcs_id,
        mov_stack[0], mov_stack[1],
    )
    plt.title(title)
    plot_lags(minlag, maxlag2)
    plt.setp(ax1.get_xticklabels(), visible=False)

    plt.subplot(gs[1], sharey=ax1)
    plt.plot(alldt.mean(axis=0), alldt.columns, c="k")
    plt.grid()
    plot_lags(minlag, maxlag2)
    plt.axvline(-maxDt, c="r", ls="--")
    plt.axvline( maxDt, c="r", ls="--")
    plt.xlabel("dt")
    plt.ylabel("Lag Time (s)")

    ax2 = plt.subplot(gs[2], sharex=ax1, sharey=ax1)
    plt.imshow(allcoh.T, extent=xextent, aspect="auto",
               interpolation="none", origin="lower", cmap="hot",
               vmin=minCoh, vmax=1)
    cb = plt.colorbar()
    cb.set_label("mean coherence")
    plt.ylabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c="k")
    plt.grid()
    locator = AutoDateLocator()
    ax2.xaxis.set_major_locator(locator)
    ax2.xaxis.set_major_formatter(AutoDateFormatter(locator))
    plt.setp(plt.xticks()[1], rotation=30, ha="right")
    plt.title("%s : %s : mean coherence" % (sta1, sta2))
    plot_lags(minlag, maxlag2)

    plt.subplot(gs[3], sharey=ax1)
    m = allcoh.mean(axis=0)
    s = allcoh.std(axis=0)
    plt.plot(m, allcoh.columns, c="k")
    plt.fill_betweenx(allcoh.columns, m - s, m + s, color="silver")
    plt.grid()
    plot_lags(minlag, maxlag2)
    plt.axvline(minCoh, c="r", ls="--")
    plt.xlabel("Coherence")
    plt.ylabel("Lag Time (s)")

    plt.tight_layout()

    if outfile:
        if outfile.startswith("?"):
            outfile = outfile.replace(
                "?", "%s-%s-%s-f%i-w%i-m%s_%s" % (
                    sta1.replace(".", "-"), sta2.replace(".", "-"),
                    components, filter_id, mwcs_id,
                    mov_stack[0], mov_stack[1],
                )
            )
        outfile = "mwcs_" + outfile
        logger.info("Saving to: %s", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()

if __name__ == "__main__":
    main("", "")
