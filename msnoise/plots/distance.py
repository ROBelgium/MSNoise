"""
Plots the REF stacks vs interstation distance. This could help deciding which
parameters to use in the dt/t calculation step. Passing ``--refilter`` allows
to bandpass filter CCFs before plotting (new in 1.5). It is also possible to
only draw CCFs for pairs including one station by passing ``--virtual-pair``
followed by the desired ``NET.STA`` (new in 1.5).

.. include:: clickhelp/msnoise-plot-distance.rst

Example:

``msnoise plot distance`` will plot all defaults:

.. image:: .static/distance.png

"""
# plot interferogram

import matplotlib.pyplot as plt

from ..api import *


def main(filterid, components, ampli=1, show=True, outfile=None,
         refilter=None, virtual_source=None, **kwargs):
    db = connect()

    pairs = get_station_pairs(db, used=1)
    cc_sampling_rate = float(get_config(db, 'cc_sampling_rate'))
    export_format = get_config(db, 'export_format')
    if export_format == "BOTH":
        extension = ".MSEED"
    else:
        extension = "."+export_format
    maxlag = float(get_config(db, 'maxlag'))
    maxlagsamples = get_maxlag_samples(db)
    t = np.linspace(-maxlag, maxlag, maxlagsamples)

    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)

    plt.figure()
    dists = []
    for pair in pairs:
        station1, station2 = pair

        dist = get_interstation_distance(station1, station2,
                                         station1.coordinates)
        dists.append(dist)

        sta1 = "%s.%s" % (station1.net, station1.sta)
        sta2 = "%s.%s" % (station2.net, station2.sta)

        if virtual_source is not None:
            if virtual_source not in [sta1, sta2]:
                continue

        pair = "%s:%s" % (sta1, sta2)
        print(pair, dist)
        ref_name = pair.replace('.', '_').replace(':', '_')
        rf = os.path.join("STACKS", "%02i" %
                          filterid, "REF", components, ref_name + extension)
        if os.path.isfile(rf):
            ref = read(rf)[0]
            if refilter:
                ref.detrend("simple")
                ref.taper(0.02)
                ref.filter("bandpass", freqmin=freqmin, freqmax=freqmax,
                           zerophase=True)
            ref.normalize()
            ref = ref.data * ampli
            plt.plot(t, ref+dist, c='k', lw=0.4)
        
    plt.ylabel("Interstation Distance in km")
    plt.xlabel("Lag Time")
    low = high = 0.0
    for filterdb in get_filters(db, all=True):
        if filterid == filterdb.ref:
            low = float(filterdb.low)
            high = float(filterdb.high)
            break
    title = '%s, Filter %d (%.2f - %.2f Hz)' % \
            (components, filterid, low, high,)
    if refilter:
        title += ", Re-filtered (%.2f - %.2f Hz)" % (freqmin, freqmax)
    plt.title(title)


    colors = ['r', 'g', 'b']
    for i, velocity in enumerate([3.0, 2.0, 1.0]):
        plt.plot([0, -1*np.max(dists)/velocity], [0, np.max(dists)],
                 c=colors[i], label='%.1f $km s^{-1}$' % velocity)
        plt.plot([0, np.max(dists)/velocity], [0, np.max(dists)], c=colors[i])

    if "xlim" in kwargs:
        plt.xlim(kwargs["xlim"][0], kwargs["xlim"][1])
    else:
        plt.xlim(-maxlag, maxlag)
    plt.xlabel("Time Lag (s)")
    plt.grid(True)
    plt.legend(loc=4)
    if outfile:
        if outfile.startswith("?"):
            newname = 'distance %s-f%i' % (components,
                                           filterid)
            outfile = outfile.replace('?', newname)
        print("output to:", outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
