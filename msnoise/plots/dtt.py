"""
This plots dt (delay time) against t (time lag). It shows the results
from the MWCS step, plus the calculated regression lines M0 and M.
The errors in the regression lines are also plotted as fainter lines.
The time lags used to calculate the regression are shown in blue.


.. include:: clickhelp/msnoise-plot-dtt.rst

Example

``msnoise plot dtt Z7.HRIM Z7.LIND 2014-08-10 -f 14 -m 20`` will plot:

.. image:: .static/dtt.png

.. versionadded:: 1.4 (Thanks to C.G. Donaldson)
"""

import matplotlib.pyplot as plt

from ..api import *


def main(sta1, sta2, filterid, components, day, mov_stack=1, show=True, outfile=None):
    db = connect()
    dtt_lag = get_config(db, "dtt_lag")
    dtt_v = float(get_config(db, "dtt_v"))
    dtt_minlag = float(get_config(db, "dtt_minlag"))
    dtt_width = float(get_config(db, "dtt_width"))
    dtt_sides = get_config(db, "dtt_sides")
    dbmaxlag = int(float(get_config(db, "maxlag")))
    sta1 = sta1.replace('.','_')
    sta2 = sta2.replace('.','_')
    if sta2 > sta1: # alphabetical order filtering!
        pair = "%s_%s"%(sta1,sta2)
        station1 = sta1.split("_")
        station2 = sta2.split("_")

        station1 = get_station(db, station1[0], station1[1])
        station2 = get_station(db, station2[0], station2[1])

        if dtt_lag == "static":
            minlag = dtt_minlag
            maxlag = minlag + dtt_width
        else:
            minlag = get_interstation_distance(station1, station2, station1.coordinates) / dtt_v
            maxlag = minlag + dtt_width

        fname = os.path.join('MWCS', "%02i" % filterid, "%03i_DAYS" %mov_stack, components, pair, '%s.txt' % day)
        print( fname)
        t = []
        dt = []
        err = []
        if os.path.isfile(fname):
            #df = pd.read_csv(fname, delimiter=' ', header=None, index_col=0, names=['t', 'dt', 'err', 'coh'])
            df = pd.read_csv(fname, delimiter=' ', header=None, names=['t', 'dt', 'err', 'coh'])
            t = df["t"].tolist()
            dt = df["dt"].tolist()
            err = df["err"].tolist()
            del df

        alldf = []
        fname = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %mov_stack, components, '%s.txt' % day)
        if not os.path.isfile(fname):
            return
        df = pd.read_csv(fname, delimiter=',')
        alldf.append(df)
        alldf = pd.concat(alldf)
        line = alldf[alldf['Pairs'] == pair].copy()
        print(line)
        M = (line["M"])
        M0 = (line["M0"])
        A = (line["A"])
        EA = (line["EA"])
        EM = (line["EM"])
        EM0 = (line["EM0"])

        plt.scatter(t, dt)
        plt.errorbar(t, dt, yerr=err, linestyle="None")
        plt.xlabel("Time (s)")
        plt.ylabel("Delay time (s)")
        plt.axvspan(-maxlag, -minlag, 0,1, color='b', alpha=0.5)
        plt.axvspan(minlag, maxlag, 0,1, color='b', alpha=0.5)
        xlineM0 = range(-dbmaxlag, dbmaxlag, 5)
        ylineM0 = []
        ylineEM0min = []
        ylineEM0max = []
        for i in range(-dbmaxlag, dbmaxlag, 5):
            ylineM0.append(M0 * i)
            ylineEM0min.append((M0-EM0) * i)
            ylineEM0max.append((M0+EM0) * i)
        plt.plot(xlineM0, ylineM0, 'r', label='M0')
        plt.plot(xlineM0, ylineEM0min, 'r', alpha=0.3)
        plt.plot(xlineM0, ylineEM0max, 'r', alpha=0.3)
        xlineM = range(-dbmaxlag, dbmaxlag, 5)
        ylineM = []
        ylineEMmax = []
        ylineEMmin = []
        for i in range(-dbmaxlag, dbmaxlag, 5):
            ylineM.append((M * i) + A)
            ylineEMmin.append(((M-EM) * i) + A)
            ylineEMmax.append(((M+EM) * i) + A)
        plt.plot(xlineM, ylineM, 'k', label='M')
        plt.plot(xlineM, ylineEMmin, 'k', alpha=0.3)
        plt.plot(xlineM, ylineEMmax, 'k', alpha=0.3)
        name = '%s-%s f%i m%i %s' % (sta1, sta2, filterid, mov_stack, day)
        name = name.replace('_', '.')
        plt.suptitle(name)
        plt.legend()
        plt.grid(True, ls="-",lw=0.2)
        if outfile:
            if outfile.startswith("?"):
                basename = '%s-%s-f%i-m%i-%s' % (sta1, sta2, filterid,
                                                 mov_stack,day)
                outfile = outfile.replace('?', basename)
            outfile = "dtt_" + outfile
            print("output to: %s" % outfile)
            plt.savefig(outfile)

        if show:
            plt.show()


if __name__ == "__main__":
    main()