"""
This plots dt (delay time) against t (time lag) for a single day. It shows
the raw MWCS measurements as a scatter plot plus the M and M0 regression lines
computed by the MWCS-DTT step.

.. include:: ../clickhelp/msnoise-cc-dtt-plot-dtt.rst

Example:

``msnoise cc dtt plot dtt BE.UCC.-- BE.MEM.-- 2023-06-15 -f 1 -w 1 -d 1``
will plot MWCS scatter + DTT regression for that day:

.. image:: .static/dtt.png

.. versionadded:: 1.4 (Thanks to C.G. Donaldson)
"""

import matplotlib.pyplot as plt
import numpy as np

from ..api import (
    connect, get_params, get_logger,
    get_station, get_interstation_distance, check_stations_uniqueness,
    get_config_set_details, xr_get_mwcs, xr_get_dtt,
    lineage_str_to_steps, get_merged_params_for_lineage, resolve_lineage_from_ids,
)

def main(sta1, sta2, filter_id=1, components="ZZ", day=None,
         preprocess_id=1, cc_id=1, stack_id=1, stack_item=1,
         mwcs_id=1, dtt_id=1,
         show=True, outfile=None, loglevel="INFO"):
    """Plot dt vs t scatter and regression lines for a single day.

    :param sta1: Station 1 (NET.STA.LOC).
    :param sta2: Station 2 (NET.STA.LOC, ≥ sta1 alphabetically).
    :param filter_id: Filter set number.
    :param components: Component pair string.
    :param day: Date string in ``YYYY-MM-DD`` format.
    :param preprocess_id: Preprocessing step set number.
    :param cc_id: CC step set number.
    :param stack_id: Stack step set number.
    :param stack_item: 1-based index into ``params.mov_stack``.
    :param mwcs_id: MWCS step set number.
    :param dtt_id: MWCS-DTT step set number.
    :param show: Display interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_dtt", loglevel, with_pid=True)
    db = connect()
    params = get_params(db)

    if sta2 < sta1:
        logger.error("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    # Resolve MWCS lineage (for raw scatter points)
    mwcs_lineage, params = resolve_lineage_from_ids(db, params, preprocess_id=preprocess_id, cc_id=cc_id, filter_id=filter_id, stack_id=stack_id, mwcs_id=mwcs_id)
    mov_stack = params.mov_stack[stack_item - 1]

    # Resolve DTT lineage (for regression lines)
    dtt_lineage, _ = resolve_lineage_from_ids(db, params, preprocess_id=preprocess_id, cc_id=cc_id, filter_id=filter_id, stack_id=stack_id, mwcs_id=mwcs_id, mwcs_dtt_id=dtt_id)

    # DTT lag window params
    dtt_lag    = getattr(params, "dtt_lag",    "static")
    dtt_v      = float(getattr(params, "dtt_v",     1.0))
    dtt_minlag = float(getattr(params, "dtt_minlag", 5.0))
    dtt_width  = float(getattr(params, "dtt_width",  30.0))
    maxlag     = float(getattr(params, "maxlag",     120.0))

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
        mwcs_all = xr_get_mwcs(params.output_folder, mwcs_lineage,
                                  sta1, sta2, components, mov_stack)
    except FileNotFoundError as fp:
        logger.error(f"MWCS FILE DOES NOT EXIST: {fp}")
        return

    import pandas as pd
    day_ts = pd.Timestamp(day)
    day_data = mwcs_all[mwcs_all.index.floor("D") == day_ts]
    if day_data.empty:
        logger.error(f"No MWCS data found for day {day}")
        return

    t   = day_data.columns.get_level_values("taxis").astype(float)
    dt  = day_data["M"].values.squeeze()
    err = day_data["EM"].values.squeeze()

    # --- Load DTT regression for the requested day ---
    try:
        dtt_all = xr_get_dtt(params.output_folder, dtt_lineage,
                               sta1, sta2, components, mov_stack)
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
        day_dtt = dtt_all[dtt_all.index.floor("D") == day_ts]
        if not day_dtt.empty:
            M   = float(day_dtt["m"].iloc[0])
            M0  = float(day_dtt["m0"].iloc[0])
            A   = float(day_dtt["a"].iloc[0])
            EM  = float(day_dtt["em"].iloc[0])
            EM0 = float(day_dtt["em0"].iloc[0])

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
        preprocess_id, cc_id, filter_id, stack_id, stack_item,
        mwcs_id, dtt_id,
    )
    ax.set_title(title)
    plt.tight_layout()

    if outfile:
        if outfile.startswith("?"):
            outfile = outfile.replace(
                "?", "%s-%s-%s-f%i-w%i-d%i-%s" % (
                    sta1.replace(".", "-"), sta2.replace(".", "-"),
                    components, filter_id, mwcs_id, dtt_id, day,
                )
            )
        outfile = "dtt_" + outfile
        logger.info(f"Saving to: {outfile}")
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()

if __name__ == "__main__":
    main("", "")
