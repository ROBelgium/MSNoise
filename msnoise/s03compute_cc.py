""" This code is responsible for the computation of the cross-correlation
functions.

This script will group *jobs* marked "T"odo in the database by day and process
them using the following scheme. As soon as one day is selected, the
corresponding jobs are marked "I"n Progress in the database. This allows
running several instances of this script in parallel.


Waveform Preprocessing
~~~~~~~~~~~~~~~~~~~~~~
Pairs are first split and a station list is created. The database is then
queried to get file paths. For each station, all files potentially containing
data for the day are opened. The traces are then merged and splitted, to obtain
the most continuous chunks possible. The different chunks are then demeaned,
tapered and merged again to a 1-day long trace. If shorter than 1-day, the
trace is padded with zeros. If longer, it is cut to match the start/end of the
day.

Each 1-day long trace is then low-passed (at ``preprocess_lowpass`` Hz),
high-passed (at ``preprocess_highpass`` Hz), then decimated/downsampled.
Decimation/Downsampling are configurable (``resampling_method``) and users are
advised testing both. One advantage of Downsampling over Decimation is that
it is able to downsample the data by any factor, not only integer factors.

.. warning::
    For an unknown reason, the PAZ-correction has disappeard from the
    current sqlvolution on GitHub: CHECK!


Processing
~~~~~~~~~~

Once all traces are preprocessed, station pairs are processed sequentially.
If a component different from *ZZ* is to be computed, the traces are first
rotated. This supposes the user has provided the station coordinates in the
*station* table. The rotation is computed for Radial and Transverse components:

.. code-block:: python

    R = tramef_N * np.cos(cplAz * np.pi / 180.) + tramef_E * np.sin(cplAz * np.pi / 180.)
    T = tramef_N * np.sin(cplAz * np.pi / 180.) - tramef_E * np.cos(cplAz * np.pi / 180.)

Then, for each ``corr_duration`` window in the signal, and for each filter
configured in the database, the traces are clipped to ``windsorizing`` times
the RMS and then whitened (see :ref:`whiten`) between the frequency bounds.
When both traces are ready, the cross-correlation function is computed
(see :ref:`mycorr`). The function returned contains data for time lags 
corresponding to ``maxlag`` in the acausal (negative lags) and causal
(positive lags) parts.

If configured (setting ``keep_all`` to 'Y'), each ``corr_duration`` CCF is
saved to the hard disk. By default, the ``keep_days`` setting is set to True
and so "N = 1 day / corr_duration" CCF are stacked and saved to the hard disk
in the STACKS/001_DAYS folder.

Once done, each job is marked "D"one in the database.

To run this script:

.. code-block:: sh

    $ msnoise compute_cc
"""

import time
import calendar
import sys

from obspy.core import utcdatetime
from scikits.samplerate import resample

from api import *
from myCorr import myCorr
from whiten import whiten


def preprocess(db, stations, comps, goal_day, params, tramef_Z, tramef_E = np.array([]), tramef_N = np.array([])):

    datafilesZ = {}
    datafilesE = {}
    datafilesN = {}

    for station in stations:
        datafilesZ[station] = []
        datafilesE[station] = []
        datafilesN[station] = []
        net, sta = station.split('.')
        gd = datetime.datetime.strptime(goal_day, '%Y-%m-%d')
        files = get_data_availability(
            db, net=net, sta=sta, starttime=gd, endtime=gd)
        for file in files:
            comp = file.comp
            fullpath = os.path.join(file.path, file.file)
            if comp[-1] == 'Z':
                datafilesZ[station].append(fullpath)
            elif comp[-1] == 'E':
                datafilesE[station].append(fullpath)
            elif comp[-1] == 'N':
                datafilesN[station].append(fullpath)

    j = 0
    for istation, station in enumerate(stations):
        for comp in comps:
            files = eval("datafiles%s['%s']" % (comp, station))
            if len(files) != 0:
                logging.debug("%s.%s Reading %i Files" %
                              (station, comp, len(files)))
                stream = Stream()
                for file in sorted(files):
                    st = read(file, dytpe=np.float)
                    for tr in st:
                        tr.data = tr.data.astype(np.float)
                    stream += st
                    del st

                logging.debug("Checking sample alignment")
                for i, trace in enumerate(stream):
                    stream[i] = check_and_phase_shift(trace)

                stream.sort()
                logging.debug("Checking Gaps")
                if len(getGaps(stream)) > 0:
                    max_gap = 10
                    only_too_long=False
                    while getGaps(stream) and not only_too_long:
                        too_long = 0
                        gaps = getGaps(stream)
                        for gap in gaps:
                            if int(gap[-1]) <= max_gap:
                                stream[gap[0]] = stream[gap[0]].__add__(stream[gap[1]], method=0, fill_value="interpolate")
                                stream.remove(stream[gap[1]])
                                break
                            else:
                                too_long += 1
                        if too_long == len(gaps):
                            only_too_long = True

                taper_length = 20.0 #seconds
                for trace in stream:
                    if trace.stats.npts < 4 * taper_length*trace.stats.sampling_rate:
                        trace.data = np.zeros(trace.stats.npts)
                    else:
                        trace.detrend(type="demean")
                        trace.detrend(type="linear")
                        taper_1s = taper_length * float(trace.stats.sampling_rate) / trace.stats.npts
                        cp = cosTaper(trace.stats.npts, taper_1s)
                        trace.data *= cp
                stream.merge(method=0, fill_value=0.0)

                logging.debug("%s.%s Slicing Stream to %s:%s" % (station, comp, utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')), utcdatetime.UTCDateTime(goal_day.replace('-', '')) + params.goal_duration - stream[0].stats.delta))
                stream[0].trim(utcdatetime.UTCDateTime(goal_day.replace('-', '')), utcdatetime.UTCDateTime(
                    goal_day.replace('-', '')) + params.goal_duration - stream[0].stats.delta, pad=True, fill_value=0.0,
                    nearest_sample=False)
                trace = stream[0]
                
                logging.debug(
                    "%s.%s Highpass at %.2f Hz" % (station, comp, params.preprocess_highpass))
                trace.filter("highpass", freq=params.preprocess_highpass, zerophase=True)
                
                if trace.stats.sampling_rate != params.goal_sampling_rate:
                    logging.debug(
                        "%s.%s Lowpass at %.2f Hz" % (station, comp, params.preprocess_lowpass))
                    trace.filter("lowpass", freq=params.preprocess_lowpass, zerophase=True)

                    

                    if params.resampling_method == "Resample":
                        logging.debug("%s.%s Downsample to %.1f Hz" %
                                      (station, comp, params.goal_sampling_rate))
                        trace.data = resample(
                            trace.data, params.goal_sampling_rate / trace.stats.sampling_rate, 'sinc_fastest')

                    elif params.resampling_method == "Decimate":
                        logging.debug("%s.%s Decimate by a factor of %i" %
                                      (station, comp, params.decimation_factor))
                        trace.data = trace.data[::params.decimation_factor]
                    trace.stats.sampling_rate = params.goal_sampling_rate

                year, month, day, hourf, minf, secf, wday, yday, isdst = trace.stats.starttime.utctimetuple()

                if j == 0:
                    t = time.strptime("%04i:%02i:%02i:%02i:%02i:%02i" %
                                      (year, month, day, hourf, minf, secf), "%Y:%m:%d:%H:%M:%S")
                    basetime = calendar.timegm(t)

                if len(trace.data) % 2 != 0:
                    trace.data = np.append(trace.data, 0.)

                if comp == "Z":
                    tramef_Z[istation] = trace.data
                elif comp == "E":
                    tramef_E[istation] = trace.data
                elif comp == "N":
                    tramef_N[istation] = trace.data

                del trace, stream
    if len(tramef_E) != 0:
        return basetime, tramef_Z, tramef_E, tramef_N
    else:
        return basetime, tramef_Z

class Params():
    pass


def main():
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: Compute CC ***')

    # Connection to the DB
    db = connect()

    if len(get_filters(db, all=False)) == 0:
        logging.info("NO FILTERS DEFINED, exiting")
        sys.exit()

    # Get Configuration
    params = Params()
    params.goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    params.goal_duration = float(get_config(db, "analysis_duration"))
    params.overlap = float(get_config(db, "overlap"))
    params.maxlag = float(get_config(db, "maxlag"))
    params.min30 = float(get_config(db, "corr_duration")) * params.goal_sampling_rate
    params.windsorizing = float(get_config(db, "windsorizing"))
    params.resampling_method = get_config(db, "resampling_method")
    params.decimation_factor = int(get_config(db, "decimation_factor"))
    params.preprocess_lowpass = float(get_config(db, "preprocess_lowpass"))
    params.preprocess_highpass = float(get_config(db, "preprocess_highpass"))
    params.keep_all = get_config(db, 'keep_all', isbool=True)
    params.keep_days = get_config(db, 'keep_days', isbool=True)
    params.components_to_compute = []
    for comp in ['ZZ', 'RR', 'TT', 'TR', 'RT', 'ZR', 'RZ', 'TZ', 'ZT']:
        if get_config(db, comp, isbool=True):
            params.components_to_compute.append(comp)
    logging.info("Will compute %s" % " ".join(params.components_to_compute))

    while is_next_job(db, type='CC'):
        jobs = get_next_job(db, type='CC')
        stations = []
        pairs = []
        refs = []

        for job in jobs:
            refs.append(job.ref)
            pairs.append(job.pair)
            netsta1, netsta2 = job.pair.split(':')
            stations.append(netsta1)
            stations.append(netsta2)
            goal_day = job.day

        stations = np.unique(stations)

        logging.info("New CC Job: %s (%i pairs with %i stations)" %
                     (goal_day, len(pairs), len(stations)))
        jt = time.time()

        xlen = int(params.goal_duration * params.goal_sampling_rate)

        if ''.join(params.components_to_compute).count('R') > 0 or ''.join(params.components_to_compute).count('T') > 0:
            comps = ['Z', 'E', 'N']
            tramef_Z = np.zeros((len(stations), xlen))
            tramef_E = np.zeros((len(stations), xlen))
            tramef_N = np.zeros((len(stations), xlen))
            basetime, tramef_Z, tramef_E, tramef_N = preprocess(db, stations, comps, goal_day, params, tramef_Z, tramef_E, tramef_N)

        else:
            comps = ['Z']
            tramef_Z = np.zeros((len(stations), xlen))
            basetime, tramef_Z = preprocess(db, stations, comps, goal_day, params, tramef_Z)


        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.goal_sampling_rate
        # Calculate the number of slices

        slices = int(params.goal_duration * params.goal_sampling_rate / params.min30)
        begins = []
        ends = []
        i = 0
        while i <=  (params.goal_duration - params.min30/params.goal_sampling_rate):
            begins.append(int(i * params.goal_sampling_rate))
            ends.append(int(i * params.goal_sampling_rate + params.min30))
            i += int(params.min30/params.goal_sampling_rate * (1.0-params.overlap))
        slices = len(begins)

        #
        # Computing only ZZ components ? Then we can be much faster:
        #

        #if False:
        if len(params.components_to_compute) == 1 and params.components_to_compute[0] == "ZZ":
            Nfft = params.min30
            if params.min30 / 2 % 2 != 0:
                Nfft = params.min30 + 2
            cp = cosTaper(int(params.min30), 0.04)

            logging.info("Pre-Whitening Traces")
            whitened_slices = np.zeros((len(stations), len(get_filters(db, all=False)), slices, int(Nfft)), dtype=np.complex)
            for istation, station in enumerate(stations):
                for islice, (begin, end) in enumerate(zip(begins,ends)):
                    tmp = tramef_Z[istation, begin:end]
                    rmsmat = np.std(np.abs(tmp))
                    if params.windsorizing == -1:
                        tmp = np.sign(tmp)
                    elif params.windsorizing != 0:
                        indexes = np.where(
                            np.abs(tmp) > (params.windsorizing * rmsmat))[0]
                        tmp[indexes] = (tmp[indexes] / np.abs(
                            tmp[indexes])) * params.windsorizing * rmsmat
                    tmp *= cp
                    for ifilter, filter in enumerate(get_filters(db, all=False)):
                        whitened_slices[istation, ifilter, islice,:] = whiten(tmp, Nfft, dt, float(filter.low), float(filter.high), plot=False)
                    del tmp
            del tramef_Z
            logging.info("Processing CC")
            for ifilter, filter in enumerate(get_filters(db, all=False)):
                for pair in pairs:
                    orig_pair = pair
                    if params.keep_all:
                        allcorr = {}
                    if params.keep_days:
                        daycorr = np.zeros(get_maxlag_samples(db,))
                        ndaycorr = 0
                    station1, station2 = pair.split(':')
                    pair = (np.where(stations == station1)
                            [0][0], np.where(stations == station2)[0][0])
                    for islice in range(slices):
                        tmp = np.vstack((whitened_slices[pair[0], ifilter, islice],
                                         whitened_slices[pair[1], ifilter, islice]))
                        corr = myCorr(tmp, np.ceil(params.maxlag / dt), plot=False)
                        tmptime = time.gmtime(basetime + begins[islice] /
                                                  params.goal_sampling_rate)
                        thisdate = time.strftime("%Y-%m-%d", tmptime)
                        thistime = time.strftime("%Y-%m-%d %H:%M:%S",
                                                 tmptime)
                        if not np.any(np.isnan(corr)) and not np.any(np.isinf(corr)):
                            if params.keep_all:
                                ccfid = "%s_%s_%s_%s_%s" % (station1, station2,
                                                            filter.ref, 'ZZ',
                                                            thisdate)
                                if ccfid not in allcorr:
                                    allcorr[ccfid] = {}
                                allcorr[ccfid][thistime] = corr

                            if params.keep_days:
                                    daycorr += corr
                                    ndaycorr += 1

                    if params.keep_all:
                        for ccfid in allcorr.keys():
                            export_allcorr(db, ccfid, allcorr[ccfid])

                    if params.keep_days:
                        thisdate = time.strftime(
                            "%Y-%m-%d", time.gmtime(basetime))
                        thistime = time.strftime(
                            "%H_%M", time.gmtime(basetime))
                        add_corr(
                            db, station1.replace(
                                '.', '_'), station2.replace('.', '_'), filter.ref,
                            thisdate, thistime, params.min30 / params.goal_sampling_rate, 'ZZ', daycorr, params.goal_sampling_rate, day=True, ncorr=ndaycorr)
                    update_job(db, goal_day, orig_pair, 'CC', 'D')
            logging.info("Job Finished. It took %.2f seconds" % (time.time() - jt))

        else:
        # if 1:
        # print '##### ITERATING OVER PAIRS #####'
            for pair in pairs:
                orig_pair = pair

                logging.info('Processing pair: %s' % pair.replace(':', ' vs '))
                tt = time.time()
                station1, station2 = pair.split(':')
                pair = (np.where(stations == station1)
                        [0][0], np.where(stations == station2)[0][0])

                s1 = get_station(db, station1.split('.')[0], station1.split('.')[1])
                s2 = get_station(db, station2.split('.')[0], station2.split('.')[1])

                if s1.X:
                    X0 = s1.X
                    Y0 = s1.Y
                    c0 = s1.coordinates

                    X1 = s2.X
                    Y1 = s2.Y
                    c1 = s2.coordinates

                    if c0 == c1:
                        coordinates = c0
                    else:
                        coordinates = 'MIX'

                    cplAz = np.deg2rad(azimuth(coordinates, X0, Y0, X1, Y1))
                    logging.debug("Azimuth=%.1f"%np.rad2deg(cplAz))
                else:
                    # logging.debug('No Coordinates found! Skipping azimuth calculation!')
                    cplAz = 0.

                for components in params.components_to_compute:
                    
                    if components == "ZZ":
                        t1 = tramef_Z[pair[0]]
                        t2 = tramef_Z[pair[1]]
                    elif components[0] == "Z":
                        t1 = tramef_Z[pair[0]]
                        t2 = tramef_E[pair[1]]
                    elif components[1] == "Z":
                        t1 = tramef_E[pair[0]]
                        t2 = tramef_Z[pair[1]]
                    else:
                        t1 = tramef_E[pair[0]]
                        t2 = tramef_E[pair[1]]
                    if np.all(t1 == 0) or np.all(t2 == 0):
                        logging.debug("%s contains empty trace(s), skipping"%components)
                        continue
                    del t1, t2
                    
                    if components[0] == "Z":
                        t1 = tramef_Z[pair[0]]
                    elif components[0] == "R":
                        if cplAz != 0:
                            t1 = tramef_N[pair[0]] * np.cos(cplAz) +\
                                 tramef_E[pair[0]] * np.sin(cplAz)
                        else:
                            t1 = tramef_E[pair[0]]

                    elif components[0] == "T":
                        if cplAz != 0:
                            t1 = tramef_N[pair[0]] * np.sin(cplAz) -\
                                 tramef_E[pair[0]] * np.cos(cplAz)
                        else:
                            t1 = tramef_N[pair[0]]

                    if components[1] == "Z":
                        t2 = tramef_Z[pair[1]]
                    elif components[1] == "R":
                        if cplAz != 0:
                            t2 = tramef_N[pair[1]] * np.cos(cplAz) +\
                                 tramef_E[pair[1]] * np.sin(cplAz)
                        else:
                            t2 = tramef_E[pair[1]]
                    elif components[1] == "T":
                        if cplAz != 0:
                            t2 = tramef_N[pair[1]] * np.sin(cplAz) -\
                                 tramef_E[pair[1]] * np.cos(cplAz)
                        else:
                            t2 = tramef_N[pair[1]]

                    trames = np.vstack((t1, t2))
                    del t1, t2

                    daycorr = {}
                    ndaycorr = {}
                    allcorr = {}
                    for filterdb in get_filters(db, all=False):
                        filterid = filterdb.ref
                        daycorr[filterid] = np.zeros(get_maxlag_samples(db,))
                        ndaycorr[filterid] = 0

                    for islice, (begin, end) in enumerate(zip(begins, ends)):
                        # print "Progress: %#2d/%2d"% (islice+1,slices)
                        trame2h = trames[:, begin:end]

                        rmsmat = np.std(np.abs(trame2h), axis=1)
                        for filterdb in get_filters(db, all=False):
                            filterid = filterdb.ref
                            low = float(filterdb.low)
                            high = float(filterdb.high)
                            rms_threshold = filterdb.rms_threshold

                            Nfft = params.min30
                            if params.min30 / 2 % 2 != 0:
                                Nfft = params.min30 + 2

                            trames2hWb = np.zeros((2, int(Nfft)), dtype=np.complex)
                            skip = False
                            for i, station in enumerate(pair):
                                if rmsmat[i] > rms_threshold:
                                    cp = cosTaper(len(trame2h[i]),0.04)
                                    trame2h[i] -= trame2h[i].mean()
                                    
                                    if params.windsorizing == -1:
                                        trame2h[i] = np.sign(trame2h[i])
                                    elif params.windsorizing != 0:
                                        indexes = np.where(
                                            np.abs(trame2h[i]) > (params.windsorizing * rmsmat[i]))[0]
                                        # clipping at windsorizing*rms
                                        trame2h[i][indexes] = (trame2h[i][indexes] / np.abs(
                                            trame2h[i][indexes])) * params.windsorizing * rmsmat[i]

                                    trames2hWb[i] = whiten(
                                        trame2h[i]*cp, Nfft, dt, low, high, plot=False)
                                else:
                                    trames2hWb[i] = np.zeros(int(Nfft))
                                    skip = True
                                    logging.debug('Slice is Zeros!')
                            if not skip:
                                corr = myCorr(trames2hWb, np.ceil(params.maxlag / dt), plot=False)
                                tmptime = time.gmtime(basetime + begin /
                                                      params.goal_sampling_rate)
                                thisdate = time.strftime("%Y-%m-%d", tmptime)
                                thistime = time.strftime("%Y-%m-%d %H:%M:%S",
                                                         tmptime)
                                if params.keep_all:
                                    ccfid = "%s_%s_%s_%s_%s" % (station1, station2,
                                                             filterid, components,
                                                             thisdate)
                                    if ccfid not in allcorr:
                                        allcorr[ccfid] = {}
                                    allcorr[ccfid][thistime] = corr

                                if params.keep_days:
                                    if not np.any(np.isnan(corr)) and \
                                            not np.any(np.isinf(corr)):
                                        daycorr[filterid] += corr
                                        ndaycorr[filterid] += 1

                                del corr, thistime, trames2hWb

                    if params.keep_all:
                        for ccfid in allcorr.keys():
                            export_allcorr(db, ccfid, allcorr[ccfid])

                    if params.keep_days:
                        try:
                            for filterdb in get_filters(db, all=False):
                                filterid = filterdb.ref
                                corr = daycorr[filterid]
                                ncorr = ndaycorr[filterid]
                                if ncorr > 0:
                                    logging.debug(
                                        "Saving daily CCF for filter %02i, comp %s (stack of %02i CCF)" % (filterid, components, ncorr))

                                    thisdate = time.strftime(
                                        "%Y-%m-%d", time.gmtime(basetime))
                                    thistime = time.strftime(
                                        "%H_%M", time.gmtime(basetime))
                                    add_corr(
                                        db, station1.replace('.', '_'),
                                        station2.replace('.', '_'), filterid,
                                        thisdate, thistime,  params.min30 /
                                        params.goal_sampling_rate,
                                        components, corr,
                                        params.goal_sampling_rate, day=True,
                                        ncorr=ncorr)
                                del corr, ncorr
                        except Exception as e:
                            logging.debug(str(e))
                    del trames, daycorr, ndaycorr
                logging.debug("Updating Job")
                update_job(db, goal_day, orig_pair, 'CC', 'D')

                logging.info("Finished processing this pair. It took %.2f seconds" %
                              (time.time() - tt))
            logging.info("Job Finished. It took %.2f seconds" % (time.time() - jt))
    logging.info('*** Finished: Compute CC ***')

if __name__ == "__main__":
    main()    

