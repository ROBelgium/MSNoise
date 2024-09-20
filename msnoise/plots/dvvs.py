"""
This plot shows the final output of MSNoise for the stretching method.

Example:

``msnoise plot dvvs`` will plot all defaults.

"""

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from matplotlib.dates import DateFormatter
from matplotlib.dates import MonthLocator

from msnoise.api import *



def main(mov_stackid=None, components='ZZ', filterid=1, pairs=[],
         show=False, outfile=None):

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

    low = high = 0.0
    low = high = 0.0
    filter = get_filters(db, ref=filterid)
    low = float(filter.low)
    high = float(filter.high)

    gs = gridspec.GridSpec(len(mov_stacks), 1)
    fig = plt.figure(figsize=(12, 9))
    plt.subplots_adjust(bottom=0.06, hspace=0.3)
    first_plot = True
    for i, mov_stack in enumerate(mov_stacks):
        alldf = []
        for comp in components:
            filedir = os.path.join("STR2","%02i" % filterid,
                               "%s_%s" % (mov_stack[0], mov_stack[1]), comp)

            listfiles = os.listdir(path=filedir)
            for file in listfiles:
                rf = os.path.join("STR2","%02i" % filterid,
                                  "%s_%s" % (mov_stack[0], mov_stack[1]), comp, file)

                # Append all series and give them the pair names
                s = pd.read_csv(rf, index_col=0, parse_dates=True).iloc[:,0]
                s = pd.Series(s, name=file[:-4])

                alldf.append(s)

        if len(alldf) == 0:
            print("No Data for %s m%s f%i" % (components, mov_stack, filterid))
            continue

        alldf = pd.concat(alldf, axis=1)

        if 'alldf' in locals():
            # Plot dvv mean and median or multiple dvv
            if first_plot == 1:
                ax = plt.subplot(gs[i])
            else:
                plt.subplot(gs[i], sharex=ax)

            alldf_mean = (alldf.groupby(axis=1, level=0).mean()-1)*100
            for pair in pairs:
                print(pair)
                pair1 = alldf_mean[pair].copy()
                print(pair1.head())
                plt.plot(pair1.index, pair1, label=pair)

            tmp2 = (alldf.mean(axis=1)-1)*100
            print(mov_stack, tmp2.head())
            plt.plot(tmp2.index, tmp2, label="mean")
            tmp3 = (alldf.median(axis=1)-1)*100
            plt.plot(tmp3.index, tmp3, label="median")
            plt.ylabel('dv/v (%)')

            if first_plot == 1:
                plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=4,
                           ncol=2, borderaxespad=0.)
                left, right = tmp2.index[0], tmp2.index[-1]
                if mov_stack == 1:
                    plt.title('1 Day')
                else:
                    plt.title('%s Moving Window' % str(mov_stack))
                first_plot = False
            else:
                plt.xlim(left, right)
                plt.title('%s Moving Window' % str(mov_stack))

            plt.grid(True)
            plt.gca().xaxis.set_major_formatter(DateFormatter("%Y-%m-%d %H:%M"))
            fig.autofmt_xdate()
            title = '%s, Filter %d (%.2f - %.2f Hz)' % \
                    (",".join(components), filterid, low, high)
            plt.suptitle(title)
            del alldf, alldf_mean

    if outfile:
        if outfile.startswith("?"):
            if len(mov_stacks) == 1:
                outfile = outfile.replace('?', '%s-f%i-m%i-M%s' % (components,
                                                                   filterid,
                                                                   mov_stack,
                                                                   "STR"))
            else:
                outfile = outfile.replace('?', '%s-f%i-M%s' % (components,
                                                               filterid,
                                                               "STR"))
        outfile = "dvv_" + outfile
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()


if __name__ == "__main__":
    main()
