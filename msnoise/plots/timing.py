"""
This plot shows the per-pair and network-mean dt/t timeseries from the
MWCS-DTT step, equivalent to the old ``timing`` plot but reading from the
lineage-aware NetCDF store.

.. include:: clickhelp/msnoise-cc-dtt-plot-timing.rst

Example:

``msnoise cc dtt plot timing`` will plot all defaults:

.. image:: .static/dvv.png
"""

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from ..core.db import connect, get_logger
from ..core.config import build_plot_outfile, get_config_set_details
from ..core.stations import get_station_pairs
from ..core.workflow import build_movstack_datelist
from ..results import MSNoiseResult


def main(mov_stackid=None, dttname="m", components="ZZ",
         filterid=1, mwcsid=1, dttid=1,
         pairs=None, showALL=False,
         show=False, outfile=None, loglevel="INFO"):
    """Plot per-pair and network-mean dt/t timeseries.

    :param mov_stackid: 1-based index into ``params.stack.mov_stack`` (0/None = all).
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

    # ------------------------------------------------------------------ #
    # Resolve lineage via MSNoiseResult                                    #
    # ------------------------------------------------------------------ #
    filter_step = f"filter_{filterid}"
    mwcs_step   = f"mwcs_{mwcsid}"
    dtt_step    = f"mwcs_dtt_{dttid}"

    all_results = MSNoiseResult.list(db, "mwcs_dtt")
    if not all_results:
        logger.error("No completed mwcs_dtt jobs found in the database.")
        return

    result = next(
        (r for r in all_results
         if filter_step in r.lineage_names
         and mwcs_step in r.lineage_names
         and r.lineage_names[-1] == dtt_step),
        None
    )
    if result is None:
        logger.error(
            f"No completed mwcs_dtt lineage for filter_{filterid}/mwcs_{mwcsid}/mwcs_dtt_{dttid}. "
            f"Available: {['/'.join(r.lineage_names) for r in all_results]}"
        )
        return

    logger.info(f"Using lineage: {'/'.join(result.lineage_names)}")
    params = result.params

    # ------------------------------------------------------------------ #
    # Moving stacks                                                        #
    # ------------------------------------------------------------------ #
    build_movstack_datelist(db)

    if mov_stackid and mov_stackid != 0:
        mov_stacks = [params.stack.mov_stack[mov_stackid - 1]]
    else:
        mov_stacks = params.stack.mov_stack

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

        for (s1, s2) in all_pairs:
            try:
                df = result.get_mwcs_dtt(f"{s1}:{s2}", components, mov_stack,
                                          )
            except FileNotFoundError:
                continue

            _keys = list(df.DTT.coords["keys"].values)
            if col not in _keys:
                continue

            ts = df.DTT.sel(keys=col)   # DataArray (times,)
            pair_series.append(ts)

            # Highlighted pairs
            pair_key = "%s:%s" % (s1, s2)
            if pair_key in highlight:
                es = df.DTT.sel(keys=err_col) if err_col in _keys else None
                ax.plot(ts.coords["times"].values, ts.values * -100,
                        label=pair_key, alpha=0.8)
                if es is not None:
                    ax.fill_between(ts.coords["times"].values,
                                    (ts.values - es.values) * -100,
                                    (ts.values + es.values) * -100,
                                    alpha=0.2)

        if not pair_series:
            logger.warning(f"No data for mov_stack={mov_stack} comp={components}")
            continue

        import xarray as _xr
        import numpy as _np
        all_da = _xr.concat(pair_series, dim="pair")
        t = all_da.coords["times"].values
        mean_ts   = all_da.mean("pair").values * -100
        median_ts = _np.nanmedian(all_da.values, axis=0) * -100

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
        outfile = build_plot_outfile(
            outfile, "timing", result.lineage_names,
            components=components, extra=col)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()


if __name__ == "__main__":
    main()
