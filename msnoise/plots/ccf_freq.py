"""

TODO FILL THIS

"""
# plot ccf_freq

import datetime

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Cursor
from obspy.signal.filter import bandpass

from msnoise.api import build_movstack_datelist, connect, get_config, \
    get_filters, get_results


def main(sta1, sta2, filterid, components, mov_stack=1, ampli=5, show=False,
         outfile=None, refilter=None, startdate=None, enddate=None):

    db = connect()
    cc_sampling_rate = float(get_config(db, 'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db, startdate, enddate)
    base = mdates.date2num(start)
    sta1 = sta1.replace('.', '_')
    sta2 = sta2.replace('.', '_')

    # TODO: Height adjustment of the plot for large number of stacks.
    # Preferably interactive
    fig = plt.figure(figsize=(12, 9))

    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)

    if sta2 >= sta1:
        pair = "%s:%s" % (sta1, sta2)

        print("New Data for %s-%s-%i-%i" % (pair, components, filterid,
                                            mov_stack))
        nstack, stack_total = get_results(db, sta1, sta2, filterid, components,
                                          datelist, mov_stack, format="matrix")
        ax = fig.add_subplot(111)
        for i, line in enumerate(stack_total):
            if np.all(np.isnan(line)):
                continue

            freq, line = prepare_abs_postitive_fft(line, cc_sampling_rate)

            if refilter:
                line = bandpass(line, freqmin, freqmax, cc_sampling_rate,
                                zerophase=True)

            ax.plot(freq, line * ampli + i + base, c='k')

        for filterdb in get_filters(db, all=True):
            if filterid == filterdb.ref:
                low = float(filterdb.low)
                high = float(filterdb.high)
                break

        ax.set_ylim(start-datetime.timedelta(days=ampli),
                    end+datetime.timedelta(days=ampli))
        ax.yaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))

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

        if outfile:
            if outfile.startswith("?"):
                pair = pair.replace(':', '-')
                outfile = outfile.replace('?', '%s-%s-f%i-m%i' % (pair,
                                                                  components,
                                                                  filterid,
                                                                  mov_stack))
            outfile = "ccftime " + outfile
            print("output to:", outfile)
            fig.savefig(outfile)
        if show:
            fig.show()
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
    freq = [x for x in freq if x >= 0]

    val = val[:len(freq)]

    return freq, val
