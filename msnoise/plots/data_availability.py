"""
Plots the data availability, as contained in the database. Every day which
has a least some data will be coloured in red. Days with no data remain blank.


.. include:: ../clickhelp/msnoise-plot-data_availability.rst


Example:

``msnoise plot data_availability`` :

.. image:: ../.static/data_availability.png

"""

import matplotlib.dates
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.dates import date2num

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

import re

from ..api import *


def main(chan, show=False, outfile=None):
    db = connect()
    start, end, datelist = build_movstack_datelist(db)
    dates = []
    stations = []
    used_stations = get_stations(db, format="seed_id")
    if not len(used_stations):
        logging.info("There seems to be no configured stations\n"
                     "Did you run 'msnoise db update_loc_chan'\n"
                     "to populate the channels/location codes?")
        sys.exit(1)
    # this is needed in case users use classic FDSN requests with ? as wildcard
    chan = chan.replace("?",".")
    matcher = re.compile(chan)
    for day in datelist:
        daystart = datetime.datetime.combine(day, datetime.time(0, 0, 0))
        dayend = datetime.datetime.combine(day, datetime.time(23, 59, 59))
        data = get_data_availability(db, starttime=daystart, endtime=dayend)
        for di in data:
            if chan.count('.') and matcher.search(di.chan):
                #print("%s matches %s" % (di.chan, chan))
                pass
            elif di.chan != chan:
                continue
            _ = "%s.%s.%s.%s" % (di.net, di.sta, di.loc, di.chan)
            if _ in used_stations:
                stations.append(_)
                dates.append(di.starttime)

    data = pd.DataFrame({"stations": stations}, index=dates)
    data = data.groupby('stations')

    llen = (end-start).days + 1
    ngroups = len(data.groups.keys())
    matrix = np.zeros((ngroups, llen))
    start = datetime.datetime.combine(start, datetime.time(0, 0, 0))

    for i, group in enumerate(sorted(data.groups.keys())):
        new = True
        for di in data.groups[group]:
            if new:
                print(group, di)
                new = False
            dt = (di-start).days
            matrix[i, dt] = 1
        print(di)

    gs = gridspec.GridSpec(2, 1, height_ratios=[4, 1])

    plt.figure(figsize=(12, 9))
    ax = plt.subplot(gs[0])
    plt.imshow(matrix, interpolation="none", aspect='auto', cmap='bwr',
               vmin=-1, vmax=1, extent=(date2num(start), date2num(end),
                                        0, ngroups),
               origin='lower')

    plt.yticks(np.arange(ngroups)+0.5, sorted(data.groups.keys()))
    ax.xaxis.set_major_locator(
        matplotlib.dates.MonthLocator())

    ax.xaxis.set_major_formatter(
        matplotlib.dates.DateFormatter('%Y-%m-%d')
    )
    plt.gcf().autofmt_xdate()
    plt.grid()

    ax = plt.subplot(gs[1], sharex=ax)
    plt.plot(datelist, np.sum(matrix, axis=0))
    ax.set_ylim((-0.1, np.amax(np.sum(matrix, axis=0))+0.1))
    plt.ylabel('N stations')
    plt.gcf().autofmt_xdate()
    plt.grid()
    plt.tight_layout()
    if outfile:
        if outfile.startswith("?"):
            now = datetime.datetime.now()
            now = now.strftime('data availability on %Y-%m-%d %H.%M.%S')
            outfile = outfile.replace('?', now)
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
