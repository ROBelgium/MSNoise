# PARAMETERS FOR MATPLOTLIB :
import pandas as pd

import numpy as np
import datetime
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from statsmodels.tsa.tsatools import detrend

from ..database_tools import *

def wavg(group):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    return wavg

def wstd(group):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wstd


def get_wavgwstd(data):
    grouped = data.groupby(level=0)
    g = grouped.apply(wavg)
    h = grouped.apply(wstd)
    return g, h

def main():
    db = connect()

    if get_config(db, name="autocorr") in ['Y', 'y', '1', 1]:
        condition = "sta2 >= sta1"
    else:
        condition = "sta2 > sta1"

    start, end, datelist = build_movstack_datelist(db)


    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    filterid = 1
    components = "ZZ"
    datatype = 'msnoise'

    gs = gridspec.GridSpec(len(mov_stacks), 1)
    plt.figure(figsize=(15, 20))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    for i, mov_stack in enumerate(mov_stacks):
        current = start
        first = True
        while current <= end:
            day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                               mov_stack, components, '%s.txt' % current)
            if os.path.isfile(day):
                df = pd.read_csv(day, header=0, index_col=0, parse_dates=True)
                if first:
                    alldf = df
                    first = False
                else:
                    alldf = alldf.append(df)
                del df
            current += datetime.timedelta(days=1)
        if 'alldf' in locals():
            global dttname, errname
            dttname = "M0"
            errname = "EM0"

            alldf[dttname] *= -100
            alldf[errname] *= -100

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()

            py1_wmean, py1_wstd = get_wavgwstd(allbut)
            py1_wmean = py1_wmean.resample('D', how='mean')
            py1_wstd = py1_wstd.resample('D', how='mean').fillna(0.0)

            data = detrend(py1_wmean)

            plt.subplot(gs[i])
            plt.plot(ALL.index, ALL[dttname],c='r',label='ALL: $\delta v/v$ of the mean network')
            #plt.fill_between(ALL.index,ALL[dttname]-ALL[errname],ALL[dttname]+ALL[errname],lw=1,color='red',zorder=-1,alpha=0.3)
            plt.plot(py1_wmean.index, data,c='g',lw=1,zorder=11,label='Weighted mean of $\delta v/v$ of individual pairs')
            plt.fill_between(py1_wmean.index, data+py1_wstd,data-py1_wstd,color='g',lw=1,zorder=-1,alpha=0.3)
            plt.ylabel('$\delta v/v$ in %')

            plt.ylim(0.5, -0.5)

            if first_plot == 1:
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                           ncol=1, borderaxespad=0.)
                left, right = plt.xlim()
                if mov_stack == 1:
                    plt.title('1 Day')
                else:
                    plt.title('%i Days Moving Window' % mov_stack)
                first_plot = False
            else:
                plt.xlim(left, right)
                plt.title('%i Days Moving Window' % mov_stack)

            plt.grid()
            del alldf

    plt.savefig('dtt_allmovstacksNEW_%s.png' % dttname, dpi=300)
    plt.show()

if __name__ == "__main__":
    main()