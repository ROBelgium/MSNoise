"""
Plot dv/v from the Stretching method.

Reads pre-aggregated network dv/v from the ``stretching_dvv`` step output
written by :mod:`msnoise.s07_compute_dvv`.

Example:

``msnoise cc dtt plot dvvs`` will plot all defaults.

``msnoise cc dtt plot dvvs -f 2 -m 1 -c ZZ`` will plot filter 2, mov_stack 1,
component ZZ.
"""

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from ..db import connect, get_logger
from ..config import build_plot_outfile, get_config_set_details
from ..results import MSNoiseResult


def main(mov_stackid=None, components="ZZ", pair_type="CC", filterid=1, stretchingid=1,
         dvvid=1, pairs=None, show=False, outfile=None, loglevel="INFO"):
    """Plot network-level dv/v from the stretching aggregate.

    Requires the ``stretching_dvv`` step to have been run first
    (``msnoise dtt compute_stretching_dvv``).

    :param mov_stackid: 1-based index into ``params.mov_stack``.
        ``None`` or ``0`` plots all moving-stack windows.
    :param components: Component pair string, comma-separated for multiple.
    :param filterid: Filter set number.
    :param stretchingid: Stretching config set number.
    :param dvvid: ``stretching_dvv`` config set number.
    :param show: Display the figure interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger("msnoise.cc_dtt_plot_dvvs", loglevel, with_pid=True)
    db = connect()

    # Resolve lineage: find a stretching_dvv result that matches
    filter_step     = f"filter_{filterid}"
    stretching_step = f"stretching_{stretchingid}"
    dvv_step        = f"stretching_dvv_{dvvid}"

    all_results = MSNoiseResult.list(db, "stretching_dvv")
    if not all_results:
        logger.error(
            "No completed stretching_dvv results found. "
            "Run 'msnoise dtt compute_stretching_dvv' first."
        )
        return

    result = next(
        (r for r in all_results
         if filter_step     in r.lineage_names
         and stretching_step in r.lineage_names
         and r.lineage_names[-1] == dvv_step),
        None
    )
    if result is None:
        logger.error(
            f"No stretching_dvv result for filter_{filterid} / "
            f"stretching_{stretchingid} / dvv_{dvvid}. "
            f"Available: {['/'.join(r.lineage_names) for r in all_results]}"
        )
        return

    logger.info(f"Using lineage: {'/'.join(result.lineage_names)}")
    params = result.params

    if mov_stackid and mov_stackid != 0:
        mov_stacks = [params.mov_stack[mov_stackid - 1]]
    else:
        mov_stacks = params.mov_stack

    comp_list = [c.strip() for c in components.split(",")]

    low = high = 0.0
    filter_params = get_config_set_details(db, "filter", filterid,
                                           format="AttribDict")
    if filter_params:
        low  = float(filter_params.freqmin)
        high = float(filter_params.freqmax)

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
                    f"No stretching_dvv data for mov_stack={mov_stack} comp={comp}. "
                    "Run 'msnoise dtt compute_stretching_dvv' first."
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
            ax.set_title("Stack %i (%s)" % (mov_stackid or 1, stack_label))
        else:
            ax.set_title("Stack %i (%s)" % (i + 1, stack_label))
            if left is not None:
                ax.set_xlim(left, right)

    fig.autofmt_xdate()
    plt.suptitle(
        "%s  |  Filter %d (%.2f – %.2f Hz)  |  Stretching dv/v"
        % (",".join(comp_list), filterid, low, high)
    )

    if outfile:
        outfile = build_plot_outfile(
            outfile, "dvv_stretching", result.lineage_names,
            components=",".join(comp_list))
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()


if __name__ == "__main__":
    main()
