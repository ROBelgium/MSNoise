"""
Plot dv/v from the MWCS → dt/t method.

Reads pre-aggregated network dv/v from the ``mwcs_dtt_dvv`` step output
written by :mod:`msnoise.s07_compute_dvv`.

Example:

``msnoise cc dtt plot mwcs_dtt`` will plot all defaults.

``msnoise cc dtt plot mwcs_dtt -f 2 -m 1 -c ZZ`` will plot filter 2, mov_stack 1,
component ZZ.

.. image:: ../.static/dvv.png
"""
import traceback

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from ..core.db import connect, get_logger
from ..core.config import build_plot_outfile
from ..results import MSNoiseResult


def main(preprocessid=1, ccid=1, filterid=1, stackid=1, stackid_item=None, refstackid=1,
         mwcsid=1, mwcsdttid=1, dvvid=1,
         dttname="m", components='ZZ', pair_type="CC", show=False, outfile=None, loglevel="INFO"):
    """Plot network-level dv/v from the MWCS-DTT aggregate.

    Requires the ``mwcs_dtt_dvv`` step to have been run first
    (``msnoise dtt dvv compute_mwcs_dtt_dvv``).

    :param dvvid: ``mwcs_dtt_dvv`` config set number.
    :param stackid_item: 1-based index into ``params.stack.mov_stack``.
        ``None`` or ``0`` plots all moving-stack windows.
    :param dttname: Ignored (kept for CLI compat). The aggregate always shows
        all available statistics (mean, weighted_mean, trimmed_mean, median).
    :param components: Component pair string, comma-separated for multiple.
    :param show: Display the figure interactively.
    :param outfile: Save path (``?`` = auto-name).
    :param loglevel: Logging verbosity.
    """
    logger = get_logger('msnoise.cc_dtt_plot_dvv', loglevel, with_pid=True)

    db = connect()
    result = MSNoiseResult.from_ids(
        db, preprocess=preprocessid, cc=ccid,
        filter=filterid, stack=stackid,
        refstack=refstackid, mwcs=mwcsid,
        mwcs_dtt=mwcsdttid, mwcs_dtt_dvv=dvvid,
    )
    params = result.params

    logger.info("Using lineage: %s" % "/".join(result.lineage_names))

    if stackid_item and stackid_item != 0:
        mov_stacks = [params.stack.mov_stack[stackid_item - 1]]
    else:
        mov_stacks = params.stack.mov_stack

    comp_list = [c.strip() for c in components.split(",")]

    low = float(params.filter.freqmin)
    high = float(params.filter.freqmax)

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
                    "No mwcs_dtt_dvv data for mov_stack=%s comp=%s pair_type=%s. "
                    "Run 'msnoise dtt dvv compute_mwcs_dtt_dvv' first." % (mov_stack, comp, pair_type)
                )
                continue
            except Exception:
                logger.error(traceback.format_exc())
                continue

            t = ds.coords["times"].values
            for stat_name in ("mean", "median", "trimmed_mean", "weighted_mean"):
                if stat_name in ds:
                    ax.plot(t, ds[stat_name].values,
                            label="%s: %s" % (comp, stat_name))
            for stat_name in ("q05", "q95"):
                if stat_name in ds:
                    ax.plot(t, ds[stat_name].values,
                            label="%s: %s" % (comp, stat_name), alpha=0.4)

            if left is None and len(t):
                left, right = t[0], t[-1]

        ax.set_ylabel("dv/v (%)")
        ax.grid(True)
        ax.xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
        stack_label = "%s_%s" % (mov_stack[0], mov_stack[1])
        if i == 0:
            ax.legend(bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=4,
                      ncol=2, borderaxespad=0.0)
            ax.set_title("Stack %i (%s)" % (stackid_item or 1, stack_label))
        else:
            ax.set_title("Stack %i (%s)" % (i + 1, stack_label))
            if left is not None:
                ax.set_xlim(left, right)

    fig.autofmt_xdate()
    plt.suptitle(
        "dv/v — MWCS-DTT | filter %i (%.3f-%.3f Hz) | preprocess %i / "
        "cc %i / stack %i / refstack %i / mwcs %i / mwcs_dtt %i / dvv %i" % (
            filterid, low, high,
            preprocessid, ccid, stackid, refstackid, mwcsid, mwcsdttid, dvvid
        )
    )

    if outfile:
        outfile = build_plot_outfile(
            outfile, "dvv_mwcs", result.lineage_names,
            components=components)
        if outfile:
            logger.info(f"Saving to: {outfile}")
            plt.savefig(outfile)
    if show:
        plt.show()
    plt.close()
    db.close()

