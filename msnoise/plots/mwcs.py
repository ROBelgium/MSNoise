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

from ..core.db import connect, get_logger
from ..core.config import build_plot_outfile
from ..core.stations import check_stations_uniqueness, get_interstation_distance, get_station
from ..results import MSNoiseResult

def main(sta1, sta2, preprocessid=1, ccid=1, filterid=1, stackid=1,
         stackid_item=1, refstackid=1, mwcsid=1, mwcsdttid=1,
         components="ZZ", show=True, outfile=None, loglevel="INFO"):
    """Plot MWCS dt and coherence images for a station pair.

    :param sta1: Station 1 in NET.STA.LOC format.
    :param sta2: Station 2 in NET.STA.LOC format (must be ≥ sta1 alphabetically).
    :param preprocessid: Preprocessing step set number.
    :param ccid: Cross-correlation step set number.
    :param filterid: Filter step set number.
    :param stackid: Stack step set number.
    :param stackid_item: 1-based index into ``params.stack.mov_stack``.
    :param mwcsid: MWCS step set number.
    :param mwcsdttid: MWCS-DTT step set number (used for DTT lag window params).
    :param components: Component pair string (e.g. ``'ZZ'``).
    :param show: Display the figure interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_mwcs", loglevel, with_pid=True)
    db = connect()

    # Full lineage including mwcs_dtt — for DTT lag window params
    dtt_result = MSNoiseResult.from_ids(db, preprocess=preprocessid, cc=ccid,
                                        filter=filterid, stack=stackid,
                                        refstack=refstackid, mwcs=mwcsid,
                                        mwcs_dtt=mwcsdttid)
    dtt_params = dtt_result.params

    # MWCS-only lineage for data loading
    result = MSNoiseResult.from_ids(db, preprocess=preprocessid, cc=ccid,
                                    filter=filterid, stack=stackid,
                                    refstack=refstackid, mwcs=mwcsid)
    params = result.params

    mov_stack = params.stack.mov_stack[stackid_item-1]

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    # DTT lag window parameters (from merged params or defaults)
    dtt_lag   = dtt_params.mwcs_dtt.dtt_lag
    dtt_v     = dtt_params.mwcs_dtt.dtt_v
    dtt_minlag= dtt_params.mwcs_dtt.dtt_minlag
    dtt_width = dtt_params.mwcs_dtt.dtt_width
    minCoh    = dtt_params.mwcs_dtt.dtt_mincoh
    maxDt     = dtt_params.mwcs_dtt.dtt_maxdtt
    maxlag    = dtt_params.cc.maxlag

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

    logger.info(f"Loading MWCS for {sta1}-{sta2} comp={components} mov_stack={mov_stack} lineage={'/'.join(result.lineage_names)}")

    try:
        mwcs = result.get_mwcs(f"{sta1}:{sta2}", components, mov_stack)
    except FileNotFoundError as fp:
        logger.error(f"FILE DOES NOT EXIST: {fp}")
        return

    # mwcs is a Dataset with MWCS variable (times, taxis, keys)
    _alldt  = mwcs.MWCS.sel(keys="M").resample(times="1D").mean()
    _allcoh = mwcs.MWCS.sel(keys="MCOH").resample(times="1D").mean()
    _t_num  = date2num(_alldt.coords["times"].values.astype("datetime64[ms]").astype(object))
    _taxis  = _alldt.coords["taxis"].values

    xextent = (_t_num[0], _t_num[-1], -maxlag, maxlag)

    gs = gridspec.GridSpec(2, 2, width_ratios=[3, 1], height_ratios=[1, 1])
    plt.figure(figsize=(12, 9))

    ax1 = plt.subplot(gs[0])
    im = plt.imshow(_alldt.values.T, extent=xextent, aspect="auto",
                    interpolation="none", origin="lower", cmap=cm.seismic)
    cscale = np.nanpercentile(_alldt.values, 99)
    im.set_clim(-cscale, cscale)
    cb = plt.colorbar()
    cb.set_label("dt")
    plt.ylabel("Lag Time (s)")
    plt.axhline(0, lw=0.5, c="k")
    plt.grid()
    title = ("%s : %s  |  %s\n"
             "Preprocess %i – CC %i – Filter %i – Stack %i – REF Stack %i - MWCS %i  |  "
             "mov_stack=%s_%s\n Reference lines from (future) MWCS-DTT %i") % (
        sta1, sta2, components,
        preprocessid, ccid, filterid, stackid, refstackid, mwcsid,
        mov_stack[0], mov_stack[1], mwcsdttid
    )
    plt.suptitle(title)
    plot_lags(minlag, maxlag2)
    plt.setp(ax1.get_xticklabels(), visible=False)

    plt.subplot(gs[1], sharey=ax1)
    plt.plot(_alldt.mean("times").values, _taxis, c="k")
    plt.grid()
    plot_lags(minlag, maxlag2)
    plt.axvline(-maxDt, c="r", ls="--")
    plt.axvline( maxDt, c="r", ls="--")
    plt.xlabel("dt")
    plt.ylabel("Lag Time (s)")

    ax2 = plt.subplot(gs[2], sharex=ax1, sharey=ax1)
    plt.imshow(_allcoh.values.T, extent=xextent, aspect="auto",
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
    _m = _allcoh.mean("times").values
    _s = _allcoh.std("times").values
    plt.plot(_m, _taxis, c="k")
    plt.fill_betweenx(_taxis, _m - _s, _m + _s, color="silver")
    plt.grid()
    plot_lags(minlag, maxlag2)
    plt.axvline(minCoh, c="r", ls="--")
    plt.xlabel("Coherence")
    plt.ylabel("Lag Time (s)")

    plt.tight_layout()

    if outfile:
        outfile = build_plot_outfile(
            outfile, "mwcs", result.lineage_names,
            pair=f"{sta1}:{sta2}", components=components, mov_stack=mov_stack)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()

if __name__ == "__main__":
    main("", "")
