"""
This plot shows the cross-correlation functions (CCF) vs time. The parameters
allow to plot the daily or the mov-stacked CCF. Filters and components are
selectable too. The ``--ampli`` argument allows to increase the vertical scale
of the CCFs. The ``--seismic`` shows the up-going wiggles with a black-filled
background (very heavy !). Passing ``--refilter`` allows to bandpass filter
CCFs before plotting (new in 1.5).

.. include:: clickhelp/msnoise-plot-ccftime.rst


Example:

``msnoise plot ccftime ID.KWUI ID.POSI`` will plot all defaults:

.. image:: .static/ccftime.png
"""
# plot interferogram

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.widgets import Cursor

from obspy.signal.filter import envelope as obspy_envelope
from obspy.signal.filter import bandpass
from ..api import *


def main(sta1, sta2, filterid, components, mov_stack=1, ampli=5, seismic=False,
         show=False, outfile=None, envelope=False, refilter=None,
         startdate=None, enddate=None):

    db = connect()
    cc_sampling_rate = float(get_config(db, 'cc_sampling_rate'))
    start, end, datelist = build_movstack_datelist(db, startdate, enddate)
    base = mdates.date2num(start)
    sta1 = sta1.replace('.', '_')
    sta2 = sta2.replace('.', '_')

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
            else:
                freq, line = prepare_fft(line, cc_sampling_rate)

            if refilter:
                line = bandpass(line, freqmin, freqmax, cc_sampling_rate,
                                zerophase=True)

            if envelope:
                line = obspy_envelope(line)

            # line /= line.max()
            ax.plot(freq, line * ampli + i + base, c='k')

            if seismic:
                y1 = np.ones(len(line)) * i
                y2 = line*ampli + i + base
                plt.fill_between(t, y1, y2, where=y2 >= y1, facecolor='k',
                                 interpolate=True)

        for filterdb in get_filters(db, all=True):
            if filterid == filterdb.ref:
                low = float(filterdb.low)
                high = float(filterdb.high)
                break

        ax.set_xlabel("Frequency [Hz]")
        ax.set_xscale('log')
        ax.grid()

        title = '%s : %s, %s, Filter %d (%.2f - %.2f Hz), Stack %d' %\
                (sta1.replace('_', '.'), sta2.replace('_', '.'), components,
                 filterid, low, high, mov_stack)
        if refilter:
            title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
        ax.set_title(title)

        ax.set_ylim(start-datetime.timedelta(days=ampli),
                    end+datetime.timedelta(days=ampli))
        ax.fmt_ydata = mdates.DateFormatter('%Y-%m-%d')


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
            plt.show()

def prepare_fft(line, sampling_rate):
    """
    Method that returns a positive part of FFT of provided signal along with
    a corresponding frequency vector.

    :type line: todo
    :param line: Signal to calculate fft.
    :type sampling_rate: float
    :param sampling_rate: Sampling rate of provided signal

    :rtype: tuple #TODO
    :return: TODO
    """
    val = np.fft.fft(line)
    val = np.abs(val)

    freq = np.fft.fftfreq(len(line),(1/sampling_rate))
    freq = [x for x in freq if x >=0]

    val = val[:len(freq)]

    return freq, val