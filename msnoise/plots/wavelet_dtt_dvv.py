"""
Plot dt/t (and dv/v) from the Wavelet Coherence Transform method.

The ``timeseries`` style reads pre-aggregated network dv/v from the
``wavelet_dtt_dvv`` step output written by :mod:`msnoise.s07_compute_dvv`.
The ``heatmap`` style reads per-pair WCT dt/t NetCDF files directly.

Two plot styles are available via the ``--visualize`` option:

* ``"timeseries"`` (default): weighted-mean dv/v vs time, one subplot per
  moving-stack window.
* ``"heatmap"``: 2-D time × frequency heatmap of the first available pair's
  ``dtt`` DataArray (useful for inspecting frequency dependence).

Example:

``msnoise cc dtt dvv plot wavelet_dtt_dvv`` will plot all defaults (timeseries).

``msnoise cc dtt dvv plot wavelet_dtt_dvv -v heatmap -f 1 -w 1 -d 1`` plots a
heatmap for filter 1, WCT set 1, DTT set 1.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.dates import DateFormatter

from ..core.db import connect, get_logger
from ..core.config import build_plot_outfile
from ..core.stations import get_station_pairs
from ..core.workflow import build_movstack_datelist
from ..results import MSNoiseResult


def _plot_timeseries(result, mov_stacks, comp_list,
                     filterid, wctid, dttid, freqmin, freqmax, pair_type,
                     logger):
    """Weighted-mean dv/v timeseries, one subplot per moving stack.

    Reads from the pre-aggregated ``wavelet_dtt_dvv`` step output.
    The result must be instantiated at the ``wavelet_dtt_dvv`` step.
    """
    params = result.params
    low  = params.filter.freqmin
    high = params.filter.freqmax

    fig, axes = plt.subplots(len(mov_stacks), 1, sharex=True, figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    left = right = None

    for i, mov_stack in enumerate(mov_stacks):
        ax = axes[i] if len(mov_stacks) > 1 else axes
        plt.sca(ax)

        for comp in comp_list:
            try:
                ds = result.get_dvv(pair_type=pair_type, components=comp,
                                     mov_stack=mov_stack, format="xarray")
            except (FileNotFoundError, ValueError):
                logger.warning(
                    f"No wavelet_dtt_dvv data for mov_stack={mov_stack} comp={comp}. "
                    "Run 'msnoise dtt compute_wavelet_dtt_dvv' first."
                )
                continue

            t = ds.coords["times"].values
            for stat_name, alpha in [("weighted_mean", 1.0), ("mean", 0.6),
                                      ("median", 0.6), ("trimmed_mean", 0.6)]:
                if stat_name in ds:
                    ax.plot(t, ds[stat_name].values,
                            label=f"{comp}: {stat_name}", alpha=alpha)

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

    freq_label = " [%.2f–%.2f Hz band]" % (params.wavelet_dtt.wct_dtt_freqmin, params.wavelet_dtt.wct_dtt_freqmax)
    plt.suptitle(
        "%s  |  Filter %d (%.2f–%.2f Hz)  |  WCT%s"
        % (",".join(comp_list), filterid, low, high, freq_label)
    )
    return fig


def _plot_heatmap(result, mov_stack, comp, logger):
    """2-D time × frequency heatmap for the first available pair."""
    db = result._db
    ds = None
    pair_label = ""
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            s1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
            for loc2 in sta2.locs():
                s2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                try:
                    ds = result.get_wct_dtt(f"{s1}:{s2}", comp, mov_stack)
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

    dtt_df = (ds["dtt"].to_dataframe()
              .unstack(level="frequency")
              .droplevel(0, axis=1))
    full_idx = dtt_df.index
    freqs = dtt_df.columns.astype(float)

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    variables = [("dtt", "dt/t", "RdBu_r", None, None),
                 ("err", "Error", "viridis", 0, None),
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


def main(mov_stackid=0, components="ZZ", pair_type="CC", filterid=1, wctid=1, dttid=1,
         dvvid=1, pairs=None, showALL=False, start="1970-01-01", end="2100-01-01",
         visualize="timeseries", ranges="[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]",
         show=True, outfile=None, loglevel="INFO"):
    """Plot dt/t from WCT DTT results.

    :param mov_stackid: 1-based moving-stack index (0 = all).
    :param components: Component pair string.
    :param filterid: Filter set number.
    :param wctid: WCT config set number.
    :param dttid: WCT-DTT config set number.
    :param dvvid: ``wavelet_dtt_dvv`` config set number.
    :param visualize: ``'timeseries'`` (uses pre-aggregated dvv) or ``'heatmap'``
        (uses per-pair wct_dtt data directly).
    :param show: Display the figure interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_wct", loglevel, with_pid=True)
    db = connect()

    # ------------------------------------------------------------------ #
    # Resolve lineage via MSNoiseResult                                    #
    # ------------------------------------------------------------------ #
    filter_step    = f"filter_{filterid}"
    wct_step       = f"wavelet_{wctid}"
    wct_dtt_step   = f"wavelet_dtt_{dttid}"
    dvv_step       = f"wavelet_dtt_dvv_{dvvid}"

    if visualize == "heatmap":
        # Heatmap reads per-pair WCT-DTT data — resolve at wavelet_dtt level
        all_results = MSNoiseResult.list(db, "wavelet_dtt")
        target_cat  = "wavelet_dtt"
        target_step = wct_dtt_step
    else:
        # Timeseries reads pre-aggregated dvv — resolve at wavelet_dtt_dvv level
        all_results = MSNoiseResult.list(db, "wavelet_dtt_dvv")
        target_cat  = "wavelet_dtt_dvv"
        target_step = dvv_step

    if not all_results:
        logger.error(
            f"No completed {target_cat} results found. "
            + ("Run 'msnoise dtt compute_wavelet_dtt_dvv' first."
               if target_cat == "wavelet_dtt_dvv" else "")
        )
        return

    result = next(
        (r for r in all_results
         if filter_step  in r.lineage_names
         and wct_step    in r.lineage_names
         and wct_dtt_step in r.lineage_names
         and r.lineage_names[-1] == target_step),
        None
    )
    if result is None:
        logger.error(
            f"No {target_cat} lineage found for "
            f"filter_{filterid} / wavelet_{wctid} / wavelet_dtt_{dttid}"
            + (f" / dvv_{dvvid}" if target_cat == "wavelet_dtt_dvv" else "")
            + f". Available: {['/'.join(r.lineage_names) for r in all_results]}"
        )
        return

    logger.info("Using lineage: %s" % "/".join(result.lineage_names))
    params = result.params

    # ------------------------------------------------------------------ #
    # Moving stacks                                                        #
    # ------------------------------------------------------------------ #
    build_movstack_datelist(db)

    if mov_stackid and mov_stackid != 0:
        mov_stacks = [params.stack.mov_stack[mov_stackid - 1]]
    else:
        mov_stacks = params.stack.mov_stack

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
        fig = _plot_heatmap(result, mov_stacks[0], comp_list[0], logger)
    else:
        fig = _plot_timeseries(result, mov_stacks, comp_list,
                               filterid, wctid, dttid,
                               freqmin, freqmax, pair_type, logger)

    if fig is None:
        return

    if outfile:
        outfile = build_plot_outfile(
            outfile, "dvv_wavelet", result.lineage_names,
            components=",".join(comp_list), extra=visualize)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            fig.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close(fig)
