"""
Plot dv/v from the Stretching method.

Reads per-pair stretching results (stored as NetCDF by the stretching step)
and plots network-aggregated statistics using :func:`msnoise.api.compute_dtt_stretching`.

Example:

``msnoise cc dtt plot dvvs`` will plot all defaults.

``msnoise cc dtt plot dvvs -f 2 -m 1 -c ZZ`` will plot filter 2, mov_stack 1,
component ZZ.
"""

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from ..api import (
    connect, get_params, get_logger, build_movstack_datelist,
    get_done_lineages_for_category, compute_dtt_stretching,
    resolve_lineage_params,
    get_config, get_config_set_details,
)


def main(mov_stackid=None, components="ZZ", filterid=1, stretchingid=1,
         pairs=None, show=False, outfile=None, loglevel="INFO"):
    """Plot network-level dv/v from the stretching method.

    :param mov_stackid: 1-based index into ``params.mov_stack`` to plot a
        single moving-stack window.  ``None`` or ``0`` plots all.
    :param components: Component pair string, comma-separated for multiple.
    :param filterid: Filter set number (used to select the lineage).
    :param stretchingid: Stretching config set number.
    :param pairs: Unused (kept for CLI signature compatibility).
    :param show: Display the figure interactively.
    :param outfile: Save figure to this path (``?`` triggers auto-naming).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_dvvs", loglevel, with_pid=True)
    db = connect()
    params = get_params(db)
    root = get_config(db, "output_folder") or "OUTPUT"

    # ------------------------------------------------------------------ #
    # Resolve lineage                                                      #
    # ------------------------------------------------------------------ #
    all_lineages = get_done_lineages_for_category(db, "stretching")
    if not all_lineages:
        logger.error("No completed stretching jobs found in the database.")
        return

    filter_step = "filter_%i" % filterid
    stretching_step = "stretching_%i" % stretchingid
    lineage = None
    for lin in all_lineages:
        if filter_step in lin and lin[-1] == stretching_step:
            lineage = lin
            break

    if lineage is None:
        logger.error(
            f"No completed stretching lineage found for filter_{filterid} / stretching_{stretchingid}. "
            f"Available: {['/'.join(l) for l in all_lineages]}"
        )
        return

    logger.info(f"Using lineage: {'/'.join(lineage)}")

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
    # Component list and filter bounds                                     #
    # ------------------------------------------------------------------ #
    comp_list = [c.strip() for c in components.split(",")]

    low = high = 0.0
    filter_params = get_config_set_details(db, "filter", filterid, format="AttribDict")
    if filter_params:
        low = float(filter_params.freqmin)
        high = float(filter_params.freqmax)

    # ------------------------------------------------------------------ #
    # Plot                                                                 #
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(len(mov_stacks), 1, sharex=True, figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    left = right = None

    for i, mov_stack in enumerate(mov_stacks):
        ax = axes[i] if len(mov_stacks) > 1 else axes
        plt.sca(ax)

        for comp in comp_list:
            try:
                stats = compute_dtt_stretching(
                    db, root, lineage, mov_stack,
                    components=comp, params=params,
                )
            except ValueError:
                logger.warning(f"No stretching data for mov_stack={mov_stack} comp={comp}")
                continue

            # dv/v (%) = (Delta - 1) * 100
            t = stats.index
            ax.plot(t, (stats[("Delta", "weighted_mean")] - 1) * 100,
                    label="%s: weighted mean" % comp)
            ax.plot(t, (stats[("Delta", "mean")] - 1) * 100,
                    label="%s: mean" % comp, alpha=0.6)
            ax.plot(t, (stats[("Delta", "50%")] - 1) * 100,
                    label="%s: median" % comp, alpha=0.6)
            ax.plot(t, (stats[("Delta", "trimmed_mean")] - 1) * 100,
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
            ax.set_title("Stack %i (%s)" % (mov_stackid or 1, stack_label))
        else:
            ax.set_title("Stack %i (%s)" % (i + 1, stack_label))
            if left is not None:
                ax.set_xlim(left, right)

    fig.autofmt_xdate()
    plt.suptitle(
        "%s  |  Filter %d (%.2f – %.2f Hz)  |  Stretching"
        % (",".join(comp_list), filterid, low, high)
    )

    if outfile:
        if outfile.startswith("?"):
            tag = "%s-f%i-s%i" % (",".join(comp_list), filterid, stretchingid)
            if len(mov_stacks) == 1:
                tag += "-m%s_%s" % (mov_stacks[0][0], mov_stacks[0][1])
            outfile = outfile.replace("?", tag)
        outfile = "dvvs_" + outfile
        logger.info(f"Saving to: {outfile}")
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()


if __name__ == "__main__":
    main()
