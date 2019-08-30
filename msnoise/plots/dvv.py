"""
This plot shows the final output of MSNoise.


.. include:: clickhelp/msnoise-plot-dvv.rst


Example:

``msnoise plot dvv`` will plot all defaults:

.. image:: .static/dvv.png

"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from matplotlib.dates import DateFormatter

from ..api import *


def wavg(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    return wavg


def wstd(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wstd


def get_wavgwstd(data, dttname, errname):
    grouped = data.groupby(level=0)
    g = grouped.apply(wavg, dttname=dttname, errname=errname)
    h = grouped.apply(wstd, dttname=dttname, errname=errname)
    return g, h


def main(mov_stack=None, dttname="M", components='ZZ', filterid=1,
         pairs=[], showALL=False, show=False, outfile=None):
    db = connect()

    start, end, datelist = build_movstack_datelist(db)

    if mov_stack != 0:
        mov_stacks = [mov_stack, ]
    else:
        mov_stack = get_config(db, "mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack), ]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    if components.count(","):
        components = components.split(",")
    else:
        components = [components, ]

    low = high = 0.0
    for filterdb in get_filters(db, all=True):
        if filterid == filterdb.ref:
            low = float(filterdb.low)
            high = float(filterdb.high)
            break

    gs = gridspec.GridSpec(len(mov_stacks), 1)
    fig = plt.figure(figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    for i, mov_stack in enumerate(mov_stacks):
        current = start
        first = True
        alldf = []
        while current <= end:
            for comp in components:
                day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                                   mov_stack, comp, '%s.txt' % current)
                if os.path.isfile(day):
                    df = pd.read_csv(day, header=0, index_col=0,
                                     parse_dates=True)
                    alldf.append(df)
            current += datetime.timedelta(days=1)
        if len(alldf) == 0:
            print("No Data for %s m%i f%i" % (components, mov_stack, filterid))
            continue

        alldf = pd.concat(alldf)
        print(mov_stack, alldf.head())
        if 'alldf' in locals():
            errname = "E" + dttname

            alldf[dttname] *= -100
            alldf[errname] *= -100

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()

            if first_plot == 1:
                ax = plt.subplot(gs[i])
            else:
                plt.subplot(gs[i], sharex=ax)
            # x = {}
            # for group in groups.keys():
            #     pairindex = []
            #     for j, pair in enumerate(allbut['Pairs']):
            #         net1, sta1, net2, sta2 = pair.split('_')
            #         if sta1 in groups[group] and sta2 in groups[group]:
            #             pairindex.append(j)
            #     tmp = allbut.iloc[np.array(pairindex)]
            #     tmp = tmp.resample('D', how='mean')
            #     #~ plt.plot(tmp.index, tmp[dttname], label=group)
            #     x[group] = tmp
            #
            # tmp = x["CRATER"] - x["VOLCAN"]
            # plt.plot(tmp.index, tmp[dttname], label="Crater - Volcan")

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

            tmp2 = allbut[dttname].resample('D').mean()
            plt.plot(tmp2.index, tmp2, label="mean")
            # tmp2.plot(label='mean',)

            tmp3 = allbut[dttname].resample('D').median()
            # tmp3.plot(label='median')
            plt.plot(tmp3.index, tmp3, label="median")
            plt.ylabel('dv/v (%)')

            if first_plot == 1:
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                           ncol=2, borderaxespad=0.)
                left, right = tmp2.index[0], tmp2.index[-1]
                if mov_stack == 1:
                    plt.title('1 Day')
                else:
                    plt.title('%i Days Moving Window' % mov_stack)
                first_plot = False
            else:
                plt.xlim(left, right)
                plt.title('%i Days Moving Window' % mov_stack)

            plt.grid(True)
            plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
            fig.autofmt_xdate()
            title = '%s, Filter %d (%.2f - %.2f Hz)' % \
                    (",".join(components), filterid, low, high)
            plt.suptitle(title)
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
        outfile = "dvv " + outfile
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()


if __name__ == "__main__":
    main()
