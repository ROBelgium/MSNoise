"""
This plot shows the cross-correlation functions' spectrum vs time. The
parameters allow to plot the daily or the mov-stacked CCF. Filters and
components are selectable too. The ``--ampli`` argument allows to increase the
vertical scale of the CCFs. Passing ``--refilter`` allows to bandpass filter
CCFs before computing the FFT and plotting. Passing ``--startdate`` and
``--enddate`` parameters allows to specify which period of data should be
plotted. By default the plot uses dates determined in database.

.. include:: ../clickhelp/msnoise-cc-plot-spectime.rst


Example:

``msnoise cc plot spectime YA.UV05 YA.UV11`` will plot all defaults:

.. image:: ../.static/spectime.png


Zooming in the X-axis and playing with the amplitude:

``msnoise cc plot spectime YA.UV05 YA.UV11 --xlim=0.08,1.1 --ampli=10``:

.. image:: ../.static/spectime_zoom.png


And refiltering to enhance high frequency content:

``msnoise cc plot spectime YA.UV05 YA.UV11 --xlim=0.5,1.1 --ampli=10 -r0.7:1.0``:

.. image:: ../.static/spectime_refilter.png



"""


import datetime

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Cursor
from obspy.signal.filter import bandpass

from msnoise.api import build_movstack_datelist, connect, get_config, \
    get_filters, get_results, check_stations_uniqueness, xr_get_ccf,\
    get_t_axis


def main(sta1, sta2, filterid, components, mov_stack=1, ampli=5, show=False,
         outfile=False, refilter=None, startdate=None, enddate=None, **kwargs):

    db = connect()
    cc_sampling_rate = float(get_config(db, 'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db)
    base = mdates.date2num(start)
    taxis = get_t_axis(db)
    # TODO: Height adjustment of the plot for large number of stacks.
    # Preferably interactive
    fig = plt.figure(figsize=(12, 9))

    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)

    if sta2 < sta1:
        print("Stations STA1 STA2 should be sorted alphabetically")
        return

    sta1 = check_stations_uniqueness(db, sta1)
    sta2 = check_stations_uniqueness(db, sta2)

    pair = "%s:%s" % (sta1, sta2)

    print("New Data for %s-%s-%i-%i" % (pair, components, filterid,
                                        mov_stack))
    stack_total = xr_get_ccf(sta1, sta2, components, filterid, mov_stack, taxis)

    # convert index to mdates
    stack_total.index = mdates.date2num(stack_total.index.to_pydatetime())

    if len(stack_total) == 0:
        print("No CCF found for this request")
        return
    ax = plt.subplot(111)
    for i, line in stack_total.iterrows():
        if np.all(np.isnan(line)):
            continue

        if refilter:
            line = bandpass(line, freqmin, freqmax, cc_sampling_rate,
                            zerophase=True)

        freq, line = prepare_abs_postitive_fft(line, cc_sampling_rate)
        line /= line.max()

        ax.plot(freq, line * ampli + i, c='k', lw=1)

    filter = get_filters(db, ref=filterid)
    low = float(filter.low)
    high = float(filter.high)


    ax.set_ylim(start-datetime.timedelta(days=ampli),
                end+datetime.timedelta(days=ampli))
    ax.yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"][0],kwargs["xlim"][1])

    ax.set_xlabel("Frequency [Hz]")
    ax.set_xscale('log')
    ax.grid()

    title = '%s : %s, %s, Filter %d (%.2f - %.2f Hz), Stack %d' %\
            (sta1.replace('_', '.'), sta2.replace('_', '.'), components,
             filterid, low, high, mov_stack)
    if refilter:
        title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
    ax.set_title(title)

    cursor = Cursor(ax, useblit=True, color='red', linewidth=1.2)
    print(outfile)
    if outfile:
        if outfile.startswith("?"):
            pair = pair.replace(':', '-')
            outfile = outfile.replace('?', '%s-%s-f%i-m%i' % (pair,
                                                              components,
                                                              filterid,
                                                              mov_stack))
        outfile = "spectime " + outfile
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close(fig)


def prepare_abs_postitive_fft(line, sampling_rate):
    """
    Method that returns a positive part of FFT of provided signal along with
    a corresponding frequency vector.

    :type line: numpy.ndarray
    :param line: Signal to calculate fft.
    :type sampling_rate: float
    :param sampling_rate: Sampling rate of provided signal

    :rtype: tuple(numpy.ndarray, numpy.ndarray)
    :return: Tuple of two arrays. One contains frequency vector for positive
    part of FFT, second contains positive and absolute FFT of input array.
    """
    val = np.fft.fft(line)
    val = np.abs(val)

    freq = np.fft.fftfreq(len(line), (1/sampling_rate))
    idx = np.where(freq >= 0)
    freq = freq[idx]
    val = val[idx]

    return freq, val