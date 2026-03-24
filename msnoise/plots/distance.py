"""
Plots the REF stacks vs interstation distance. This could help deciding which
parameters to use in the dt/t calculation step. Passing ``--refilter`` allows
to bandpass filter CCFs before plotting (new in 1.5). It is also possible to
only draw CCFs for pairs including one station by passing ``--virtual-pair``
followed by the desired ``NET.STA`` (new in 1.5).

.. include:: ../clickhelp/msnoise-cc-plot-distance.rst

Example:

``msnoise cc plot distance`` will plot all defaults:

.. image:: ../.static/distance.png

"""

import matplotlib.pyplot as plt
from obspy import Trace
from ..api import *


def main(filterid, components, ampli=1, show=True, outfile=None,
         refilter=None, virtual_source=None, loglevel="INFO", **kwargs):
    logger = get_logger('msnoise.plotdistance_child', loglevel,
                        with_pid=True)
    db = connect()

    # Build merged params from all configsets in the lineage — same approach
    # as batch["params"] in the processing steps.
    # global_1 is excluded from the lineage by get_stack_lineage_for_filter;
    # global params are already in get_params(db).
    # Use get_refstack_lineage_for_filter so xr_get_ref finds files under
    # the refstack step folder (.../stack_1/refstack_1/_output/REF/...).
    lineage = get_refstack_lineage_for_filter(db, filterid)
    stack_lineage = get_stack_lineage_for_filter(db, filterid)
    params = get_params(db)
    for step_name in lineage:
        category, set_num = step_name.rsplit('_', 1)
        cfg = get_config_set_details(db, category, int(set_num), format='AttribDict')
        if cfg:
            params.update(cfg)

    cc_sampling_rate = float(params.cc_sampling_rate)
    maxlag = float(params.maxlag)
    output_folder = params.output_folder or 'OUTPUT'
    taxis = np.linspace(-maxlag, maxlag, int(2 * maxlag * cc_sampling_rate) + 1)

    pairs = get_station_pairs(db, used=1)

    if refilter:
        freqmin, freqmax = refilter.split(':')
        freqmin = float(freqmin)
        freqmax = float(freqmax)

    plt.figure()
    dists = []
    for pair in pairs:
        station1, station2 = pair
        # TODO get distance for LOCids!!
        dist = get_interstation_distance(station1, station2,
                                         station1.coordinates)
        if dist == 0 and station1 != station2:
            logger.warning("Distance is 0.0 km for %s.%s:%s.%s" % (station1.net, station1.sta, station2.net, station2.sta))
        dists.append(dist)
        for loc1 in station1.locs():
            for loc2 in station2.locs():
                sta1 = "%s.%s.%s" % (station1.net, station1.sta, loc1)
                sta2 = "%s.%s.%s" % (station2.net, station2.sta, loc2)

                if virtual_source is not None:
                    if virtual_source not in [sta1, sta2]:
                        continue

                try:
                    ref = xr_get_ref(output_folder, lineage, sta1, sta2, components, filterid, taxis)
                    ref = Trace(data=ref.CCF.values.copy())
                    ref.stats.sampling_rate = cc_sampling_rate
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue

                if refilter:
                    ref.detrend("simple")
                    ref.taper(0.02)
                    ref.filter("bandpass", freqmin=freqmin, freqmax=freqmax,
                               zerophase=True)
                ref.normalize()
                ref = ref.data * ampli
                plt.plot(taxis, ref + dist, c='k', lw=0.4)

    plt.ylabel("Interstation Distance in km")
    plt.xlabel("Lag Time")
    low = float(params.freqmin)
    high = float(params.freqmax)
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
        logger.info("output to: %s" % outfile)
        plt.savefig(outfile)
    if show:
        plt.show()
    else:
        plt.close()
