"""
This plot shows the final output of MSNoise.


.. include:: clickhelp/msnoise-plot-dvv.rst


Example:

``msnoise plot dvv`` will plot all defaults:

.. image:: .static/dvv.png

"""

import pandas as pd

import numpy as np
import datetime
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from statsmodels.tsa.tsatools import detrend

from ..api import *



def main(mov_stack=None, dttname="A", components='ZZ', filterid=1,
         pairs=[], showALL=False, show=False, outfile=None):
    db = connect()

    if get_config(db, name="autocorr", isbool=True):
        condition = "sta2 >= sta1"
    else:
        condition = "sta2 > sta1"

    start, end, datelist = build_movstack_datelist(db)

    if mov_stack != 0:
        mov_stacks = [mov_stack,]
    else:
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    gs = gridspec.GridSpec(len(mov_stacks), 1)
    plt.figure(figsize=(15, 10))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    for i, mov_stack in enumerate(mov_stacks):
        current = start
        first = True
        alldf = []
        while current <= end:
            day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                               mov_stack, components, '%s.txt' % current)
            if os.path.isfile(day):
                df = pd.read_csv(day, header=0, index_col=0, parse_dates=True)
                alldf.append(df)
            current += datetime.timedelta(days=1)
        if len(alldf) == 0:
            print("No Data for %s m%i f%i" % (components, mov_stack, filterid))
            continue

        alldf = pd.concat(alldf)
        if 'alldf' in locals():
            errname = "E" + dttname

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()

            plt.subplot(gs[i])

            for pair in pairs:
                print(pair)
                pair1 = alldf[alldf['Pairs'] == pair].copy()
                print(pair1.head())
                plt.plot(pair1.index, pair1[dttname], label=pair)
                plt.fill_between(pair1.index, pair1[dttname]-pair1[errname],
                                 pair1[dttname]+pair1[errname], zorder=-1,
                                 alpha=0.5)
                pair1.to_csv('%s-m%i-f%i.csv'%(pair, mov_stack, filterid))

            if showALL:
                plt.plot(ALL.index, ALL[dttname], c='r',
                         label='ALL: $\delta v/v$ of the mean network')

            tmp2 = allbut[dttname].resample('D', how='mean')
            tmp2.plot(label='mean')

            tmp3 = allbut[dttname].resample('D', how='median')
            tmp3.plot(label='median')

            if first_plot == 1:
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                           ncol=2, borderaxespad=0.)
                left, right = plt.xlim()
                if mov_stack == 1:
                    plt.title('1 Day')
                else:
                    plt.title('%i Days Moving Window' % mov_stack)
                first_plot = False
            else:
                plt.xlim(left, right)
                plt.title('%i Days Moving Window' % mov_stack)

            plt.grid(True)
            del alldf
    if outfile:
        if outfile.startswith("?"):
            if len(mov_stacks) == 1:
                outfile = outfile.replace('?', '%s-f%i-m%i-M%s' % (components,
                                                                   filterid,
                                                                   mov_stack,
                                                                   dttname))
            else:
                outfile = outfile.replace('?', '%s-f%i-M%s' % (components,
                                                               filterid,
                                                               dttname))
        outfile = "timing " + outfile
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
if __name__ == "__main__":
    main()