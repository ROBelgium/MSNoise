"""
This plots dt (delay time) against t (time lag) for a single day. It shows
the raw MWCS measurements as a scatter plot plus the M and M0 regression lines
computed by the MWCS-DTT step.

.. include:: clickhelp/msnoise-cc-dtt-plot-dtt.rst

Example:

``msnoise cc dtt plot dtt BE.UCC.-- BE.MEM.-- 2023-06-15 -f 1 -w 1 -d 1``
will plot MWCS scatter + DTT regression for that day:

.. image:: .static/dtt.png

"""

import matplotlib.pyplot as plt
import numpy as np

from ..core.db import connect, get_logger
from ..core.config import build_plot_outfile
from ..core.stations import check_stations_uniqueness, get_interstation_distance, get_station
from ..results import MSNoiseResult

def main(sta1, sta2, filterid=1, components="ZZ", day=None,
         preprocessid=1, ccid=1, stackid=1, stackid_item=1,
         mwcsid=1, mwcsdttid=1,
         show=True, outfile=None, loglevel="INFO"):
    """Plot dt vs t scatter and regression lines for a single day.

    :param sta1: Station 1 (NET.STA.LOC).
    :param sta2: Station 2 (NET.STA.LOC, ≥ sta1 alphabetically).
    :param filterid: Filter set number.
    :param components: Component pair string.
    :param day: Date string in ``YYYY-MM-DD`` format.
    :param preprocessid: Preprocessing step set number.
    :param ccid: CC step set number.
    :param stackid: Stack step set number.
    :param stackid_item: 1-based index into ``params.stack.mov_stack``.
    :param mwcsid: MWCS step set number.
    :param mwcsdttid: MWCS-DTT step set number.
    :param show: Display interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_dtt", loglevel, with_pid=True)
    db = connect()

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    # MWCS result — for raw scatter points
    mwcs_result = MSNoiseResult.from_ids(db, preprocess=preprocessid, cc=ccid,
                                         filter=filterid, stack=stackid,
                                         refstack=1, mwcs=mwcsid)
    params = mwcs_result.params
    mov_stack = params.stack.mov_stack[stackid_item - 1]

    # DTT result — for regression lines
    dtt_result = MSNoiseResult.from_ids(db, preprocess=preprocessid, cc=ccid,
                                        filter=filterid, stack=stackid,
                                        refstack=1, mwcs=mwcsid,
                                        mwcs_dtt=mwcsdttid)

    # DTT lag window params
    dtt_lag    = params.mwcs_dtt.dtt_lag
    dtt_v      = params.mwcs_dtt.dtt_v
    dtt_minlag = params.mwcs_dtt.dtt_minlag
    dtt_width  = params.mwcs_dtt.dtt_width
    maxlag     = params.cc.maxlag

    if dtt_lag == "dynamic":
        st1 = get_station(db, *sta1.split(".")[:2])
        st2 = get_station(db, *sta2.split(".")[:2])
        minlag = get_interstation_distance(st1, st2, st1.coordinates) / dtt_v
    else:
        minlag = dtt_minlag
    maxlag2 = minlag + dtt_width

    # --- Load MWCS for the requested day ---
    logger.info(f"Loading MWCS for {sta1}-{sta2} comp={components} day={day} mov={mov_stack}")
    try:
        mwcs_all = mwcs_result.get_mwcs(f"{sta1}:{sta2}", components, mov_stack)
    except FileNotFoundError as fp:
        logger.error(f"MWCS FILE DOES NOT EXIST: {fp}")
        return

    import numpy as _np
    _day_mask = _np.array(
        [str(t)[:10] == day for t in mwcs_all.MWCS.coords["times"].values]
    )
    if not _day_mask.any():
        logger.error(f"No MWCS data found for day {day}")
        return

    _mwcs_day = mwcs_all.MWCS.isel(times=_day_mask)
    t   = _mwcs_day.coords["taxis"].values.astype(float)
    dt  = _mwcs_day.sel(keys="M").values.squeeze()
    err = _mwcs_day.sel(keys="EM").values.squeeze()

    # --- Load DTT regression for the requested day ---
    try:
        dtt_all = dtt_result.get_mwcs_dtt(f"{sta1}:{sta2}", components, mov_stack)
    except FileNotFoundError as fp:
        logger.error(f"DTT FILE DOES NOT EXIST: {fp}")
        dtt_all = None

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(t, dt, zorder=3, label="MWCS dt")
    ax.errorbar(t, dt, yerr=err, linestyle="None", alpha=0.5)

    ax.axvspan(-maxlag2, -minlag, alpha=0.15, color="b", label="regression window")
    ax.axvspan( minlag,  maxlag2, alpha=0.15, color="b")

    xline = np.linspace(-maxlag, maxlag, 200)
    if dtt_all is not None:
        _dtt_mask = _np.array(
            [str(t)[:10] == day for t in dtt_all.DTT.coords["times"].values]
        )
        if _dtt_mask.any():
            _dtt_day = dtt_all.DTT.isel(times=_dtt_mask)
            M   = float(_dtt_day.sel(keys="m").values[0])
            M0  = float(_dtt_day.sel(keys="m0").values[0])
            A   = float(_dtt_day.sel(keys="a").values[0])
            EM  = float(_dtt_day.sel(keys="em").values[0])
            EM0 = float(_dtt_day.sel(keys="em0").values[0])

            ax.plot(xline, M0 * xline,             "r",   label="M0=%.4f" % M0)
            ax.plot(xline, (M0-EM0) * xline,       "r",   alpha=0.3)
            ax.plot(xline, (M0+EM0) * xline,       "r",   alpha=0.3)
            ax.plot(xline, M * xline + A,           "k",   label="M=%.4f" % M)
            ax.plot(xline, (M-EM) * xline + A,     "k",   alpha=0.3)
            ax.plot(xline, (M+EM) * xline + A,     "k",   alpha=0.3)

    ax.set_xlim(-maxlag, maxlag)
    ax.set_xlabel("Time lag (s)")
    ax.set_ylabel("Delay time dt (s)")
    ax.legend()
    ax.grid(True, ls="-", lw=0.2)

    title = ("%s – %s  |  comp=%s  |  day=%s\n"
             "Preprocess %i – CC %i – Filter %i – Stack %i item %i"
             " – MWCS %i – DTT %i") % (
        sta1, sta2, components, day,
        preprocessid, ccid, filterid, stackid, stackid_item,
        mwcsid, mwcsdttid,
    )
    ax.set_title(title)
    plt.tight_layout()

    if outfile:
        outfile = build_plot_outfile(
            outfile, "dtt", dtt_result.lineage_names,
            pair=f"{sta1}:{sta2}", components=components, extra=day)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()

if __name__ == "__main__":
    main("", "")
