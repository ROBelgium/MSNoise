
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from msnoise.api import *

eruptions = []
eruptions.append((datetime.date(2014, 6, 21), datetime.date(2014, 6, 22), "1"))
eruptions.append((datetime.date(2015, 2, 4), datetime.date(2015, 2, 16), "2"))
eruptions.append((datetime.date(2015, 5, 17), datetime.date(2015, 5, 31), "3"))
eruptions.append((datetime.date(2015, 7, 31), datetime.date(2015, 8, 2), "4"))
eruptions.append(
    (datetime.date(2015, 8, 24), datetime.date(2015, 10, 18), "5.1"))
eruptions.append(
    (datetime.date(2015, 10, 22), datetime.date(2015, 10, 2), "5.2"))
eruptions.append(
    (datetime.date(2015, 10, 29), datetime.date(2015, 11, 1), "5.3"))
eruptions.append((datetime.date(2016, 5, 26), datetime.date(2016, 5, 27), "6"))
eruptions.append((datetime.date(2016, 9, 11), datetime.date(2016, 9, 19), "7"))
eruptions.append((datetime.date(2017, 1, 31), datetime.date(2017, 2, 28), "8"))
eruptions.append((datetime.date(2017, 5, 17), datetime.date(2017, 5, 19), "9"))
eruptions.append((datetime.date(2017, 7, 14), datetime.date(2017, 8, 30), "10"))


def plot_eruptions():
    span_a = 0.5
    span_c = "pink"
    for (start, end, label) in eruptions:
        plt.axvspan(start, end, facecolor=span_c, alpha=span_a, zorder=-1, )
        plt.axvline(start, c='k', lw=1)
        plt.axvline(end, c='k', lw=1)

    plt.axvline(datetime.date.today() - datetime.timedelta(days=1), color='k',
                alpha=0.5, zorder=-1)

    plt.axvline(datetime.date.today(), color='k', alpha=0.5, zorder=-1)

# 
def wavg(group, dttname, errname):
    d = group[dttname]
    group[errname][group[errname] == 0] = 1e-6
    w = 1. / group[errname]
    try:
        wavg = (d * w).sum() / w.sum()
    except:
        wavg = d.mean()
    return wavg

# def wavg(group):
#     print("group:", group.shape)
#     d = group["M"]
#     group["EM"][group["EM"] == 0] = 1e-6
#     w = 1. / group["EM"]
#     wavg = (d * w).sum() / w.sum()
#     return wavg

from scipy.stats import trim_mean
def trim(group):
    d = group["M"]
    return trim_mean(d, 0.1)

def wstd(group, dttname, errname):
    print("group", group)
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

from statsmodels.tsa.filters.hp_filter import hpfilter
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

    gs = gridspec.GridSpec(len(mov_stacks), 1)
    plt.figure(figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    groups = []
    try:
        from groups_config import groups
    except:
        pass
    
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
        # print(mov_stack, alldf.head())
        if 'alldf' in locals():
            errname = "E" + dttname

            alldf[dttname] *= -100
            alldf[errname] *= -100

            ALL = alldf[alldf['Pairs'] == 'ALL'].copy()
            allbut = alldf[alldf['Pairs'] != 'ALL'].copy()
            allbut.to_csv('all-m%i-f%i.csv' % (mov_stack, filterid))
            if first_plot == 1:
                ax = plt.subplot(gs[i])
            else:
                plt.subplot(gs[i], sharex=ax)
            x = {}
            for group in groups.keys():
                pairindex = []
                for j, pair in enumerate(allbut['Pairs']):
                    net1, sta1, net2, sta2 = pair.split('_')
                    if sta1 in groups[group] and sta2 in groups[group]:
                        pairindex.append(j)
                tmp = allbut.iloc[np.array(pairindex)]
                G = tmp.groupby(pd.TimeGrouper('D'))
                tmp = G.mean()
                plt.plot(tmp.index, tmp["M"], label=group)

                # tmp = G.apply(trim)
                # plt.plot(tmp.index, tmp, label=group+"-Tmean")
                # cycle, trend = hpfilter(tmp, 2)
                # plt.plot(tmp.index, trend, ls="-", lw=1)
                # x[group] = tmp

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

            # tmp2 = allbut[dttname].resample('D').mean()
            # tmp2.plot(label='mean',)
            # 
            # tmp3 = allbut[dttname].resample('D').median()
            # tmp3.plot(label='median')
            plt.ylabel('dv/v (%)')
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
            plot_eruptions()
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
        outfile = "dvv " + outfile
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
if __name__ == "__main__":
    main()
