"""
This plot shows the final output of MSNoise.


.. include:: ../clickhelp/msnoise-cc-dvv-plot-dvv.rst


Example:

``msnoise cc dvv plot dvv`` will plot all defaults:

.. image:: ../.static/dvv.png

"""
import traceback

import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from ..api import *


def main(mov_stackid=None, dttname="M", components='ZZ', filterid=1,
         pairs=[], showALL=False, show=False, outfile=None, loglevel="INFO"):
    logger = get_logger('msnoise.cc_dvv_plot_dvv', loglevel,
                        with_pid=True)
    db = connect()
    params = get_params(db)
    start, end, datelist = build_movstack_datelist(db)

    if mov_stackid and mov_stackid != "":
        mov_stack = params.mov_stack[mov_stackid - 1]
        mov_stacks = [mov_stack, ]
        print(mov_stack)
    else:
        mov_stacks = params.mov_stack

    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]
    print(mov_stacks)
    low = high = 0.0
    filter = get_filters(db, ref=filterid)
    low = float(filter.low)
    high = float(filter.high)

    fig, axes = plt.subplots(len(mov_stacks), 1, sharex=True, figsize=(12, 9))

    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    for i, mov_stack in enumerate(mov_stacks):
        current = start
        try:
            plt.sca(axes[i])
        except:
            plt.sca(axes)
        # plt.title('%s Moving Window' % mov_stack)
        for comps in components:
            try:
                dvv = xr_get_dvv(comps, filterid, mov_stack)
            except FileNotFoundError as fullpath:
                logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                continue
            for _ in ["mean", "50%", "trimmed_mean", "weighted_mean"]:
                plt.plot(dvv.index, dvv.loc[:, ("m", _)] * -100, label="%s: %s" % (comps,_ ))
            for _ in ["5%","95%"]:
                plt.plot(dvv.index, dvv.loc[:, ("m", _)] * -100, label="%s: %s" % (comps,_ ))

        plt.ylabel('dv/v (%)')
        if i == 0:
            plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                       ncol=2, borderaxespad=0.)
            left, right = dvv.index[0], dvv.index[-1]
            plt.title("Stack %i (%s_%s)"% (mov_stackid or i+1, mov_stack[0], mov_stack[1]))
        else:
            plt.title("Stack %i (%s_%s)"% (i+1, mov_stack[0], mov_stack[1]))
            plt.xlim(left, right)

        plt.grid(True)
        plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
    fig.autofmt_xdate()
    title = '%s, Filter %d (%.2f - %.2f Hz)' % \
            (",".join(components), filterid, low, high)
    plt.suptitle(title)

    if outfile:
        if outfile.startswith("?"):
            if len(mov_stacks) == 1:
                outfile = outfile.replace('?', '%s-f%i-m%s_%s-M%s' % (components,
                                                                   filterid,
                                                                   mov_stack[0],
                                                                   mov_stack[1],
                                                                   dttname))
            else:
                outfile = outfile.replace('?', '%s-f%i-M%s' % (components,
                                                               filterid,
                                                               dttname))
        outfile = "dvv " + outfile
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()


if __name__ == "__main__":
    main()
