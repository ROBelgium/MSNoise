"""
Plot dv/v from the MWCS → dt/t method.

Reads per-pair MWCS-DTT results and plots network-aggregated statistics
using :func:`msnoise.api.compute_dvv`.

Example:

``msnoise cc dtt plot dvv`` will plot all defaults.

``msnoise cc dtt plot dvv -f 2 -m 1 -c ZZ`` will plot filter 2, mov_stack 1,
component ZZ.

.. image:: ../.static/dvv.png
"""
import traceback

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from ..api import (
    connect, get_params, get_logger, build_movstack_datelist,
    get_done_lineages_for_category, compute_dvv,
    resolve_lineage_params,
    get_config, get_config_set_details, lineage_str_to_steps, get_merged_params_for_lineage
)


def main(preprocess_id=1, cc_id=1, filter_id=1, stack_id=1, stack_item=None, refstack_id=1,
         mwcs_id=1, mwcs_dtt_id=1,
         dttname="M", components='ZZ', show=False, outfile=None, loglevel="INFO"):
    """Plot network-level dv/v from the MWCS-DTT method.

    :param mov_stackid: 1-based index into ``params.mov_stack``.
        ``None`` or ``0`` plots all moving-stack windows.
    :param dttname: Column to use for dv/v (``'m'`` = slope, ``'m0'`` =
        zero-intercept slope).  Default ``'M'`` is remapped to ``'m'``.
    :param components: Component pair string, comma-separated for multiple.
    :param filterid: Filter set number (selects the lineage).
    :param mwcsid: MWCS config set number.
    :param dttid: MWCS-DTT config set number.
    :param pairs: Unused (kept for CLI compatibility).
    :param showALL: Unused (kept for CLI compatibility).
    :param show: Display the figure interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger('msnoise.cc_dtt_plot_dvv', loglevel, with_pid=True)

    db = connect()
    params = get_params(db)
    lineage_names = [f"preprocess_{preprocess_id}", f"cc_{cc_id}",
                     f"filter_{filter_id}", f"stack_{stack_id}", f"refstack_{refstack_id}",
                     f"mwcs_{mwcs_id}", f"mwcs_dtt_{mwcs_dtt_id}"]
    lineage_str = "/".join(lineage_names)
    steps = lineage_str_to_steps(db, lineage_str, "/")
    paralineage, lineage_names, params = get_merged_params_for_lineage(db, params, {}, steps)

    # Normalise dttname: legacy 'M' → 'm'
    col = dttname.lower() if dttname else 'm'

    logger.info("Using lineage: %s" % "/".join(lineage_names))

    # ------------------------------------------------------------------ #
    # Moving stacks                                                        #
    # ------------------------------------------------------------------ #


    if stack_item and stack_item != 0:
        mov_stacks = [params.mov_stack[stack_item - 1]]
    else:
        mov_stacks = params.mov_stack

    # ------------------------------------------------------------------ #
    # Component list and filter bounds                                     #
    # ------------------------------------------------------------------ #
    comp_list = [c.strip() for c in components.split(",")]

    low = float(params.freqmin)
    high = float(params.freqmax)

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
                stats = compute_dvv(
                    db, params.output_folder, lineage_names, mov_stack,
                    components=comp, params=params,
                )
            except ValueError:
                logger.warning(
                    "No mwcs_dtt data for mov_stack=%s comp=%s" % (mov_stack, comp)
                )
                continue
            except Exception:
                logger.error(traceback.format_exc())
                continue

            t = stats.index
            for stat_name in ("mean", "50%", "trimmed_mean", "weighted_mean"):
                try:
                    ax.plot(t, stats.loc[:, (col, stat_name)] * -100,
                            label="%s: %s" % (comp, stat_name))
                except KeyError:
                    pass
            for stat_name in ("5%", "95%"):
                try:
                    ax.plot(t, stats.loc[:, (col, stat_name)] * -100,
                            label="%s: %s" % (comp, stat_name), alpha=0.4)
                except KeyError:
                    pass

            if left is None and len(t):
                left, right = t[0], t[-1]

        ax.set_ylabel("dv/v (%)")
        ax.grid(True)
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
        stack_label = "%s_%s" % (mov_stack[0], mov_stack[1])
        if i == 0:
            ax.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=4,
                      ncol=2, borderaxespad=0.0)
            ax.set_title("Stack %i (%s)" % (stack_item or 1, stack_label))
        else:
            ax.set_title("Stack %i (%s)" % (i + 1, stack_label))
            if left is not None:
                ax.set_xlim(left, right)

    fig.autofmt_xdate()
    plt.suptitle(
        "%s  |  Filter %d (%.2f – %.2f Hz)  |  MWCS-DTT"
        % (",".join(comp_list), filter_id, low, high)
    )

    if outfile:
        if outfile.startswith("?"):
            # Use str(comp_list) to match original naming convention
            # e.g. ['ZZ']-f1-MM.png  so tests can assert the filename
            tag = "%s-f%i-M%s" % (str(comp_list), filter_id, dttname)
            if len(mov_stacks) == 1:
                tag += "-m%s_%s" % (mov_stacks[0][0], mov_stacks[0][1])
            outfile = outfile.replace("?", tag)
        outfile = "dvv " + outfile
        logger.info(f"Saving to: {outfile}")
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()


if __name__ == "__main__":
    main()
