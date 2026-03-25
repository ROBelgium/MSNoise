"""
Plot dt/t (and dv/v) from the Wavelet Coherence Transform method.

Reads per-pair WCT dt/t results (stored as NetCDF by the ``wavelet_dtt``
step) and plots network-aggregated statistics using
:func:`msnoise.api.compute_dtt_wct`.

Two plot styles are available via the ``--visualize`` option:

* ``"timeseries"`` (default): weighted-mean dv/v vs time, one subplot per
  moving-stack window.
* ``"heatmap"``: 2-D time × frequency heatmap of the first available pair's
  ``dvv`` DataArray (useful for inspecting frequency dependence).

Example:

``msnoise cc dtt plot wct`` will plot all defaults (timeseries).

``msnoise cc dtt plot wct -v heatmap -f 1 -w 1 -d 1`` plots a heatmap for
filter 1, WCT set 1, DTT set 1.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

from ..api import (
    connect, get_params, get_logger, build_movstack_datelist,
    get_done_lineages_for_category, compute_dtt_wct, xr_get_wct_dtt,
    resolve_lineage_params,
    get_config, get_config_set_details, get_station_pairs,
)


def _plot_timeseries(db, root, lineage, mov_stacks, comp_list,
                     filterid, wctid, dttid, freqmin, freqmax,
                     params, logger):
    """Weighted-mean dv/v timeseries, one subplot per moving stack."""
    low = params.freqmin
    high = params.freqmax

    fig, axes = plt.subplots(len(mov_stacks), 1, sharex=True, figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    left = right = None

    for i, mov_stack in enumerate(mov_stacks):
        ax = axes[i] if len(mov_stacks) > 1 else axes
        plt.sca(ax)

        for comp in comp_list:
            try:
                stats = compute_dtt_wct(
                    db, root, lineage, mov_stack,
                    components=comp, params=params,
                    freqmin=freqmin, freqmax=freqmax,
                )
            except ValueError:
                logger.warning(
                    "No WCT-DTT data for mov_stack=%s comp=%s", mov_stack, comp
                )
                continue

            t = stats.index
            # WCT stores dt/t; dv/v = -dt/t (same sign convention as MWCS)
            ax.plot(t, -stats[("dvv", "weighted_mean")],
                    label="%s: weighted mean" % comp)
            ax.plot(t, -stats[("dvv", "mean")],
                    label="%s: mean" % comp, alpha=0.6)
            ax.plot(t, -stats[("dvv", "50%")],
                    label="%s: median" % comp, alpha=0.6)
            ax.plot(t, -stats[("dvv", "trimmed_mean")],
                    label="%s: trimmed mean" % comp, alpha=0.6)

            if left is None and len(t):
                left, right = t[0], t[-1]

        ax.set_ylabel("dv/v (%)")
        ax.grid(True)
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
        stack_label = "%s_%s" % (mov_stack[0], mov_stack[1])
        if i == 0:
            ax.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=4,
                      ncol=2, borderaxespad=0.0)
            ax.set_title("Stack 1 (%s)" % stack_label)
        else:
            ax.set_title("Stack %i (%s)" % (i + 1, stack_label))
            if left is not None:
                ax.set_xlim(left, right)

    fig.autofmt_xdate()

    freq_label = " [%.2f–%.2f Hz band]" % (params.wct_dtt_freqmin, params.wct_dtt_freqmax)
    plt.suptitle(
        "%s  |  Filter %d (%.2f–%.2f Hz)  |  WCT%s"
        % (",".join(comp_list), filterid, low, high, freq_label)
    )
    return fig


def _plot_heatmap(db, root, lineage, mov_stack, comp, logger):
    """2-D time × frequency heatmap for the first available pair."""
    # Find first pair that has data
    ds = None
    pair_label = ""
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            s1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
            for loc2 in sta2.locs():
                s2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                try:
                    ds = xr_get_wct_dtt(root, lineage, s1, s2, comp, mov_stack)
                    pair_label = "%s – %s" % (s1, s2)
                    break
                except FileNotFoundError:
                    continue
            if ds is not None:
                break
        if ds is not None:
            break

    if ds is None:
        logger.error("No WCT-DTT data found for heatmap.")
        return None

    dvv_df = (ds["dvv"].to_dataframe()
              .unstack(level="frequency")
              .droplevel(0, axis=1))
    full_idx = dvv_df.index
    freqs = dvv_df.columns.astype(float)
    data = np.ma.masked_invalid(dvv_df.values)

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    variables = [("dvv", "dt/t (%)", "RdBu_r", None, None),
                 ("err", "Error (%)", "viridis", 0, None),
                 ("coh", "Coherence", "RdYlGn", 0, 1)]

    for ax, (var, label, cmap, vmin, vmax) in zip(axes, variables):
        vdf = (ds[var].to_dataframe()
               .unstack(level="frequency")
               .droplevel(0, axis=1))
        vdata = np.ma.masked_invalid(vdf.values)
        if vmin is None:
            vmin = np.nanpercentile(vdf.values, 2)
        if vmax is None:
            vmax = np.nanpercentile(vdf.values, 98)
        mesh = ax.pcolormesh(full_idx, freqs, vdata.T,
                             cmap=cmap, vmin=vmin, vmax=vmax)
        cb = fig.colorbar(mesh, ax=ax, shrink=0.8)
        cb.set_label(label, rotation=270, labelpad=15)
        ax.set_ylabel("Frequency (Hz)")
        ax.set_title(var.upper())

    axes[-1].set_xlabel("Date")
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.setp(axes[-1].get_xticklabels(), rotation=45)

    stack_label = "%s_%s" % (mov_stack[0], mov_stack[1])
    fig.suptitle("WCT dt/t – %s  %s  stack=%s" % (pair_label, comp, stack_label),
                 fontsize=13)
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    return fig


def main(mov_stackid=0, components="ZZ", filterid=1, wctid=1, dttid=1,
         pairs=None, showALL=False, start="1970-01-01", end="2100-01-01",
         visualize="timeseries", ranges="[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]",
         show=True, outfile=None, loglevel="INFO"):
    """Plot dt/t / dv/v from WCT results.

    :param mov_stackid: 1-based moving-stack index (0 = all).
    :param components: Component pair string.
    :param filterid: Filter set number.
    :param wctid: WCT config set number.
    :param dttid: WCT-DTT config set number.
    :param visualize: ``'timeseries'`` or ``'heatmap'``.
    :param show: Display the figure interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_wct", loglevel, with_pid=True)
    db = connect()
    params = get_params(db)
    root = get_config(db, "output_folder") or "OUTPUT"

    # ------------------------------------------------------------------ #
    # Resolve lineage                                                      #
    # ------------------------------------------------------------------ #
    all_lineages = get_done_lineages_for_category(db, "wavelet_dtt")
    if not all_lineages:
        logger.error("No completed wavelet_dtt jobs found in the database.")
        return

    filter_step = "filter_%i" % filterid
    wct_step = "wavelet_%i" % wctid
    wct_dtt_step = "wavelet_dtt_%i" % dttid
    lineage = None
    for lin in all_lineages:
        if (filter_step in lin and wct_step in lin
                and lin[-1] == wct_dtt_step):
            lineage = lin
            break

    if lineage is None:
        logger.error(
            "No completed wavelet_dtt lineage found for "
            "filter_%i / wavelet_%i / wavelet_dtt_%i. Available: %s",
            filterid, wctid, dttid,
            ["/".join(l) for l in all_lineages],
        )
        return

    logger.info("Using lineage: %s" % "/".join(lineage))

    # ------------------------------------------------------------------ #
    # Moving stacks                                                        #
    # ------------------------------------------------------------------ #
    build_movstack_datelist(db)
    _, _, params = resolve_lineage_params(db, lineage)

    if mov_stackid and mov_stackid != 0:
        mov_stacks = [params.mov_stack[mov_stackid - 1]]
    else:
        mov_stacks = params.mov_stack

    comp_list = [c.strip() for c in components.split(",")]

    # Parse optional frequency band from ranges string (first range used)
    freqmin = freqmax = None
    try:
        import ast
        parsed = ast.literal_eval("[%s]" % ranges.strip("[]"))
        if parsed and isinstance(parsed[0], (list, tuple)):
            freqmin, freqmax = float(parsed[0][0]), float(parsed[0][1])
    except Exception:
        pass

    # ------------------------------------------------------------------ #
    # Dispatch to plot style                                               #
    # ------------------------------------------------------------------ #
    if visualize == "heatmap":
        fig = _plot_heatmap(db, root, lineage, mov_stacks[0],
                            comp_list[0], logger)
    else:
        fig = _plot_timeseries(db, root, lineage, mov_stacks, comp_list,
                               filterid, wctid, dttid,
                               freqmin, freqmax, params, logger)

    if fig is None:
        return

    if outfile:
        if outfile.startswith("?"):
            tag = "%s-f%i-w%i-d%i" % (",".join(comp_list), filterid, wctid, dttid)
            if len(mov_stacks) == 1:
                tag += "-m%s_%s" % (mov_stacks[0][0], mov_stacks[0][1])
            outfile = outfile.replace("?", tag)
        outfile = "wct_%s_%s" % (visualize, outfile)
        logger.info("Saving to: %s", outfile)
        fig.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close(fig)
