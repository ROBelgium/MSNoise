"""
This plot shows the per-pair and network-mean dt/t timeseries from the
MWCS-DTT step, equivalent to the old ``timing`` plot but reading from the
lineage-aware NetCDF store.

.. include:: ../clickhelp/msnoise-cc-dtt-plot-timing.rst

Example:

``msnoise cc dtt plot timing`` will plot all defaults:

.. image:: .static/dvv.png
"""

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from ..api import (
    connect, get_params, get_logger, build_movstack_datelist,
    get_config, get_config_set_details,
    get_done_lineages_for_category, get_station_pairs, xr_get_dtt2,
    resolve_lineage_params,
)


def main(mov_stackid=None, dttname="m", components="ZZ",
         filterid=1, mwcsid=1, dttid=1,
         pairs=None, showALL=False,
         show=False, outfile=None, loglevel="INFO"):
    """Plot per-pair and network-mean dt/t timeseries.

    :param mov_stackid: 1-based index into ``params.mov_stack`` (0/None = all).
    :param dttname: DTT column to display: ``'m'`` (slope = dt/t) or
        ``'m0'`` (zero-intercept slope).
    :param components: Component pair string.
    :param filterid: Filter set number.
    :param mwcsid: MWCS set number.
    :param dttid: MWCS-DTT set number.
    :param pairs: List of ``'NET.STA.LOC:NET.STA.LOC'`` strings to highlight.
    :param showALL: Unused (kept for signature compatibility).
    :param show: Display interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_timing", loglevel, with_pid=True)
    db = connect()
    params = get_params(db)
    root = get_config(db, "output_folder") or "OUTPUT"

    # ------------------------------------------------------------------ #
    # Resolve lineage                                                      #
    # ------------------------------------------------------------------ #
    all_lineages = get_done_lineages_for_category(db, "mwcs_dtt")
    if not all_lineages:
        logger.error("No completed mwcs_dtt jobs found in the database.")
        return

    filter_step   = "filter_%i"   % filterid
    mwcs_step     = "mwcs_%i"     % mwcsid
    dtt_step      = "mwcs_dtt_%i" % dttid
    lineage = None
    for lin in all_lineages:
        if filter_step in lin and mwcs_step in lin and lin[-1] == dtt_step:
            lineage = lin
            break

    if lineage is None:
        logger.error(
            "No completed mwcs_dtt lineage for filter_%i/mwcs_%i/mwcs_dtt_%i. "
            "Available: %s", filterid, mwcsid, dttid,
            ["/".join(l) for l in all_lineages],
        )
        return

    logger.info("Using lineage: %s", "/".join(lineage))

    # ------------------------------------------------------------------ #
    # Moving stacks                                                        #
    # ------------------------------------------------------------------ #
    build_movstack_datelist(db)
    _, _, params = resolve_lineage_params(db, lineage)

    if mov_stackid and mov_stackid != 0:
        mov_stacks = [params.mov_stack[mov_stackid - 1]]
    else:
        mov_stacks = params.mov_stack

    # ------------------------------------------------------------------ #
    # Filter frequency bounds for title                                    #
    # ------------------------------------------------------------------ #
    low = high = 0.0
    filter_params = get_config_set_details(db, "filter", filterid,
                                           format="AttribDict")
    if filter_params:
        low  = float(filter_params.freqmin)
        high = float(filter_params.freqmax)

    # ------------------------------------------------------------------ #
    # Build pair list                                                      #
    # ------------------------------------------------------------------ #
    all_pairs = []
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            s1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
            for loc2 in sta2.locs():
                s2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                all_pairs.append((s1, s2))

    highlight = set(pairs) if pairs else set()

    # ------------------------------------------------------------------ #
    # Plot                                                                 #
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(len(mov_stacks), 1, sharex=True, figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    left = right = None

    col       = dttname.lower()   # 'm' or 'm0'
    err_col   = "e" + col         # 'em' or 'em0'

    for i, mov_stack in enumerate(mov_stacks):
        ax = axes[i] if len(mov_stacks) > 1 else axes
        plt.sca(ax)

        pair_series = []
        err_series  = []

        for (s1, s2) in all_pairs:
            try:
                df = xr_get_dtt2(root, lineage, s1, s2, components, mov_stack)
            except FileNotFoundError:
                continue

            if col not in df.columns:
                continue

            ts = df[col]
            pair_series.append(ts)

            # Highlighted pairs
            pair_key = "%s:%s" % (s1, s2)
            if pair_key in highlight:
                es = df.get(err_col)
                ax.plot(ts.index, ts.values * -100,
                        label=pair_key, alpha=0.8)
                if es is not None:
                    ax.fill_between(ts.index,
                                    (ts - es) * -100,
                                    (ts + es) * -100,
                                    alpha=0.2)

        if not pair_series:
            logger.warning("No data for mov_stack=%s comp=%s", mov_stack, components)
            continue

        import pandas as pd
        all_df = pd.concat(pair_series, axis=1)
        t = all_df.index
        mean_ts   = all_df.mean(axis=1) * -100
        median_ts = all_df.median(axis=1) * -100

        ax.plot(t, mean_ts,   label="mean",   lw=1.5)
        ax.plot(t, median_ts, label="median", lw=1.5, linestyle="--")

        ax.set_ylabel("dv/v (%)")
        ax.grid(True)
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d"))
        stack_label = "%s_%s" % (mov_stack[0], mov_stack[1])

        if i == 0:
            ax.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=4,
                      ncol=3, borderaxespad=0.0)
            ax.set_title("Stack 1 (%s)" % stack_label)
            if len(t):
                left, right = t[0], t[-1]
        else:
            ax.set_title("Stack %i (%s)" % (i + 1, stack_label))
            if left is not None:
                ax.set_xlim(left, right)

    fig.autofmt_xdate()
    plt.suptitle(
        "%s  |  Filter %d (%.2f – %.2f Hz)  |  MWCS %i – DTT %i  |  col=%s"
        % (components, filterid, low, high, mwcsid, dttid, col)
    )

    if outfile:
        if outfile.startswith("?"):
            tag = "%s-f%i-w%i-d%i-%s" % (components, filterid, mwcsid, dttid, col)
            if len(mov_stacks) == 1:
                tag += "-m%s_%s" % (mov_stacks[0][0], mov_stacks[0][1])
            outfile = outfile.replace("?", tag)
        outfile = "timing_" + outfile
        logger.info("Saving to: %s", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()


if __name__ == "__main__":
    main()
