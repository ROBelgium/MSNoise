"""
Plots the data availability, as contained in the database. Every day which
has a least some data will be coloured in red. Days with no data remain blank.


.. include:: clickhelp/msnoise-plot-data_availability.rst


Example:

``msnoise plot data_availability`` :

.. image:: .static/data_availability.png

"""

import glob
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.gridspec import GridSpec
from obspy import read_inventory, UTCDateTime
from obspy.signal import PPSD
from msnoise.api import *


def main(net, sta, comp, time_of_weekday=None, period_lim=None, show=False,
         outfile=None, cmap="viridis", color_lim=None):

    db = connect()
    logging.debug('Preloading all instrument response')
    response_format = get_config(db, 'response_format')
    response_files = glob.glob(os.path.join(get_config(db, 'response_path'),"%s.%s.*.*"%(net, sta)))

    start, end, datelist = build_movstack_datelist(db)
    first = True
    ppsd = None
    for day in datelist:
        jday = int(day.strftime("%j"))
        toglob = os.path.join('PSD', 'NPZ', "%s" % day.year, net, sta, comp+".D", "%s.%s.*.%s.D.%s.%03i.npz"%(net, sta, comp, day.year, jday))
        files = glob.glob(toglob)
        if not len(files):
            print("No files found for %s.%s %s: %s" % (net,sta,comp, day))
            continue
        file = files[0]
        if os.path.isfile(file):
            if first:
                ppsd = PPSD.load_npz(file)
                first = False
            else:
                try:
                    ppsd.add_npz(file)
                except:
                    pass
    if not ppsd:
        # f = open(os.path.join(os.path.split(__file__)[0], "nodata.png"), 'rb').read()
        # outfile.write(f)
        return

    print("Calculate histogram...")
    ppsd.calculate_histogram( time_of_weekday=time_of_weekday)
    if len(ppsd._times_processed) < 2:
        print("Not enough data for %s.%s - %s"%(net, sta, comp))
        return
    print("Plotting PPSD...")
    fig = ppsd.plot(show_mean=True, show=False, period_lim=period_lim,show_coverage=False)
    print(fig)

    print("Plotting spectrogram...")
    gs = GridSpec(2, 2, width_ratios=[10,1], hspace=0.3, left=0.2, right=0.85, top=0.80, bottom=0.2)
    fig.set_size_inches(8.2, 11.6, forward=True)
    fig.axes[0].set_subplotspec(gs[0])
    fig.axes[0].set_position(gs[0].get_position(fig))
    fig.axes[1].set_subplotspec(gs[1])
    fig.axes[1].set_position(gs[1].get_position(fig))

    fig.axes[0].set_title("")
    ax2 = plt.subplot(gs[2])
    print("Creating Time Index...")
    print(UTCDateTime(ppsd.times_processed[1]).datetime-UTCDateTime(ppsd.times_processed[0]).datetime)
    new_times = pd.date_range(UTCDateTime(ppsd.times_processed[0]).datetime,
                              UTCDateTime(ppsd.times_processed[-1]).datetime,
                              freq="1H")

    ind_times = np.array([UTCDateTime(t).datetime for t in ppsd.current_times_used])
    data = np.asarray(ppsd._binned_psds)
    print("Creating DataFrame...")
    data = pd.DataFrame(data, index=ind_times)
    # TODO : why resample here ?
    # data = data.resample("H").fillna(method="ffill")

    if color_lim:
        vmin, vmax = color_lim
        plt.pcolormesh(data.index, ppsd.period_bin_centers, data.T, cmap=cmap,
                       vmin=vmin, vmax=vmax, rasterized=True)
    else:
        plt.imshow(data.T, cmap=cmap, aspect='auto', origin='lower',
                   extent=(date2num(new_times[0].to_pydatetime()),
                           date2num(new_times[-1].to_pydatetime()),
                           ppsd.period_bin_centers[0],
                           ppsd.period_bin_centers[-1]),
                   )

    cb = plt.colorbar(cax=plt.subplot(gs[3]), use_gridspec=True)
    plt.sca(ax2)
    cb.set_label("Amplitude [dB]")
    plt.ylabel("Period [s]")
    plt.semilogy((UTCDateTime(ppsd.times_processed[0]).datetime), 0, alpha=0)

    plt.xlim(UTCDateTime(ppsd.times_processed[0]).datetime, UTCDateTime(ppsd.times_processed[-1]).datetime)
    plt.ylim(period_lim[0], period_lim[1])
    plt.suptitle(ppsd._get_plot_title())
    fig.autofmt_xdate()
    xax = plt.gca().get_xaxis()  # get the x-axis
    adf = xax.get_major_formatter()  # the the auto-formatter

    adf.scaled[1. / 24] = '%Y-%m-%d %H:%M'  # set the < 1d scale to H:M
    adf.scaled[1.0] = '%Y-%m-%d'  # set the > 1d < 1m scale to Y-m-d
    adf.scaled[30.] = '%Y-%m'  # set the > 1m < 1Y scale to Y-m
    adf.scaled[365.] = '%Y'  # set the > 1y scale to Y

    if outfile:
        print("output to:", outfile)
        plt.savefig(outfile)

    if show:
        plt.show()
    plt.gcf()
    plt.gca()
    plt.close(fig)
    del fig
