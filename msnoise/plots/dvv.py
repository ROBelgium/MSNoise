"""MSNoise is ...

Usage:
~~~~~~

.. code-block:: sh

    $ msnoise plot dvv

"""
# PARAMETERS FOR MATPLOTLIB :
import pandas as pd

import numpy as np
import datetime
import os

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from statsmodels.tsa.tsatools import detrend

from ..api import *

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

def main(all=True, mov_stack=None, savefig=False, show=False):
    db = connect()

    if get_config(db, name="autocorr") in ['Y', 'y', '1', 1]:
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
        alldf = []
        while current <= end:
            day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                               mov_stack, components, '%s.txt' % current)
            if os.path.isfile(day):
                df = pd.read_csv(day, header=0, index_col=0, parse_dates=True)
                alldf.append(df)
            current += datetime.timedelta(days=1)
        alldf = pd.concat(alldf)
        if 'alldf' in locals():
            global dttname, errname
            dttname = "M0"
            errname = "EM0"

            alldf[dttname] *= -100
            alldf[errname] *= -100

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()

            groups = {}
            groups['CRATER'] = ["UV11","UV15","FJS","FLR","SNE","UV12","FOR","RVL","UV06"]
            groups['GPENTES'] = ["UV03","UV08","UV04","UV02","HDL"]
            groups['VOLCAN'] = groups['CRATER'] + groups['GPENTES'] + ['HIM','VIL']
            
            
            plt.subplot(gs[i])
            x = {}
            for group in groups.keys():
                pairindex = []
                for j, pair in enumerate(allbut['Pairs']):
                    net1, sta1, net2, sta2 = pair.split('_')
                    if sta1 in groups[group] and sta2 in groups[group]:
                        pairindex.append(j)
                tmp = allbut.iloc[np.array(pairindex)]
                tmp = tmp.resample('D', how='mean')
                #~ plt.plot(tmp.index, tmp[dttname], label=group)
                x[group] = tmp
            
            tmp = x["CRATER"] - x["VOLCAN"]
            #~ plt.plot(tmp.index, tmp[dttname], label="Crater - Volcan")
            
            py1_wmean, py1_wstd = get_wavgwstd(allbut)
            py1_wmean = py1_wmean.resample('D', how='median')
            py1_wstd = py1_wstd.resample('D', how='mean').fillna(0.0)

            data = detrend(py1_wmean)

            pair1 = alldf[alldf['Pairs'] == 'YA_FJS_YA_SNE'].copy()
            
            #~ plt.plot(ALL.index, ALL[dttname],c='r',label='ALL: $\delta v/v$ of the mean network')
            plt.plot(pair1.index, pair1[dttname], c='b',label='pair')
            #~ plt.plot(pair2.index, pair2[dttname], c='magenta',label='pair')
            #~ r = pd.rolling_mean(pair1[dttname], 30)
            #~ plt.plot(r.index, r, c='k')
            #plt.fill_between(ALL.index,ALL[dttname]-ALL[errname],ALL[dttname]+ALL[errname],lw=1,color='red',zorder=-1,alpha=0.3)
            #~ plt.plot(py1_wmean.index, data,c='g',lw=1,zorder=11,label='Weighted mean of $\delta v/v$ of individual pairs')
            #~ plt.fill_between(py1_wmean.index, data+py1_wstd,data-py1_wstd,color='g',lw=1,zorder=-1,alpha=0.3)
            #~ plt.fill_between(py1_wmean.index, data+3*py1_wstd,data-3*py1_wstd,color='g',lw=1,zorder=-1,alpha=0.1)
            plt.ylabel('$\delta v/v$ in %')
            #~ for pair in allbut['Pairs']:
                #~ tmp = allbut[allbut['Pairs']==pair]
                #~ plt.plot(tmp.index, tmp[dttname],lw=0.5)
            # plt.ylim(0.5, -0.5)

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
    if savefig:
        plt.savefig('dtt_allmovstacksNEW_%s.png' % dttname, dpi=300)
    plt.show()

if __name__ == "__main__":
    main()