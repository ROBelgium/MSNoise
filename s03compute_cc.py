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
saved to the hard disk. By default, the ``keep_day`` setting is set to True
and so "N = 1 day / corr_duration" CCF are stacked and saved to the hard disk
in the STACKS/001_DAYS folder.

Once done, each job is marked "D"one in the database.

To run this script:

.. code-block:: sh

    python s03compute_cc.py
"""

import numpy as np
from obspy.core import read, utcdatetime, Stream
from obspy.signal import cosTaper
from obspy.signal.filter import lowpass, highpass
from scikits.samplerate import resample
import time
import calendar
import datetime
import sys
import os
import scipy.fftpack
from database_tools import *
from myCorr import myCorr
from whiten import whiten

import logging

if __name__ == "__main__":
    
    logging.basicConfig(level=logging.DEBUG,
                        filename="./compute_cc.log",
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        filemode='w')
    
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)
    
    logging.info('*** Starting: Compute CC ***')
    
    # Connection to the DB
    db = connect()
    
    # Get Configuration
    components_to_compute = []
    for comp in ['ZZ', 'RR', 'TT', 'TR', 'RT', 'ZR', 'RZ', 'TZ', 'ZT']:
        if get_config(db, comp) in ['Y', 'y', '1', 1]:
            components_to_compute.append(comp)
    
    logging.info("Will compute %s" % " ".join(components_to_compute))
    
    # allow_large_concats(db)
    
    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))  # was 20.0
    goal_duration = float(get_config(db, "analysis_duration"))  # was 86400
    maxlag = float(get_config(db, "maxlag"))
    min30 = float(get_config(db, "corr_duration")) * goal_sampling_rate
    windsorizing = float(get_config(db, "windsorizing"))
    resampling_method = get_config(db, "resampling_method")
    decimation_factor = int(get_config(db, "decimation_factor"))
    preprocess_lowpass = float(get_config(db, "preprocess_lowpass"))
    preprocess_highpass = float(get_config(db, "preprocess_highpass"))
    
    keep_all = False
    if get_config(db, 'keep_all') in ['Y', 'y', '1', 1]:
        keep_all = True
    
    keep_days = False
    if get_config(db, 'keep_days') in ['Y', 'y', '1', 1]:
        keep_days = True
    
    # Process !
    while is_next_job(db, type='CC'):
        jobs = get_next_job(db, type='CC')
        stations = []
    
        pairs = []
        refs = []
        for job in jobs:
            refs.append(job.ref)
            pairs.append(job.pair)
            goal_day = job.day
    
        for pair in pairs:
            netsta1, netsta2 = pair.split(':')
            stations.append(netsta1)
            stations.append(netsta2)
            update_job(db, goal_day, pair, 'CC', 'I')
    
        fi = len(get_filters(db, all=False))
        if fi == 0:
            logging.info("NO FILTERS DEFINED, exiting")
            sys.exit()
    
        stations = np.unique(stations)
    
        logging.info("New CC Job: %s (%i pairs with %i stations)" %
                     (goal_day, len(pairs), len(stations)))
        jt = time.time()
    
        datafilesZ = {}
        datafilesE = {}
        datafilesN = {}
    
        durations = []
    
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
    
        TimeVec = np.arange(0., goal_duration, 1. / goal_sampling_rate)
    
        if ''.join(components_to_compute).count('R') > 0 or ''.join(components_to_compute).count('T') > 0:
            comps = ['Z', 'E', 'N']
            tramef_Z = np.zeros((len(stations), len(TimeVec)))
            tramef_E = np.zeros((len(stations), len(TimeVec)))
            tramef_N = np.zeros((len(stations), len(TimeVec)))
        else:
            comps = ['Z']
            tramef_Z = np.zeros((len(stations), len(TimeVec)))
    
        j = 0
        for istation, station in enumerate(stations):
            for comp in comps:
                files = eval("datafiles%s['%s']" % (comp, station))
                if len(files) != 0:
                    logging.debug("%s.%s Reading %i Files" %
                                  (station, comp, len(files)))
                    stream = Stream()
                    for file in sorted(files):
                        st = read(file, format="MSEED")
                        stream += st
                        del st
                    stream.merge()
                    stream = stream.split()
                    for trace in stream:
                        data = trace.data
                        if len(data) > 2:
                            tp = cosTaper(len(data), 0.01)
                            data -= np.mean(data)
                            data *= tp
                            trace.data = data
                        else:
                            trace.data *= 0
                        del data
                    logging.debug("%s.%s Merging Stream" % (station, comp))
                    #fills gaps with 0s and gives only one 'Trace'
                    stream.merge(fill_value=0)
                    logging.debug("%s.%s Slicing Stream to %s:%s" % (station, comp, utcdatetime.UTCDateTime(
                        goal_day.replace('-', '')), utcdatetime.UTCDateTime(goal_day.replace('-', '')) + goal_duration - stream[0].stats.delta))
    
                    stream[0].trim(utcdatetime.UTCDateTime(goal_day.replace('-', '')), utcdatetime.UTCDateTime(
                        goal_day.replace('-', '')) + goal_duration - stream[0].stats.delta, pad=True, fill_value=0.0)
                    trace = stream[0]
    
                    data = trace.data
                    freq = preprocess_lowpass
                    logging.debug(
                        "%s.%s Lowpass at %.2f Hz" % (station, comp, freq))
                    data = lowpass(
                        trace.data, freq, trace.stats.sampling_rate, zerophase=True)
    
                    freq = preprocess_highpass
                    logging.debug(
                        "%s.%s Highpass at %.2f Hz" % (station, comp, freq))
                    data = highpass(
                        data, freq, trace.stats.sampling_rate, zerophase=True)
    
                    samplerate = trace.stats['sampling_rate']
                    if samplerate != goal_sampling_rate:
                        if resampling_method == "Resample":
                            logging.debug("%s.%s Downsample to %.1f Hz" %
                                          (station, comp, goal_sampling_rate))
                            data = resample(
                                data, goal_sampling_rate / trace.stats.sampling_rate, 'sinc_fastest')
                        elif resampling_method == "Decimate":
                            logging.debug("%s.%s Decimate by a factor of %i" %
                                          (station, comp, decimation_factor))
                            data = data[::decimation_factor]
    
                    # logging.debug('Data for %s: %s - %s' % (station, trace.stats.starttime , trace.stats.endtime))
                    # print 'Data for %s: %s - %s' % (station,
                    # trace.stats.starttime , trace.stats.endtime)
                    year, month, day, hourf, minf, secf, wday, yday, isdst = trace.stats.starttime.utctimetuple(
                    )
    
                    TimeVec = np.arange(0., goal_duration, 1. / goal_sampling_rate)
                    # trame = np.zeros(len(TimeVec))
    
                    trame = data
    
                    if j == 0:
                        FirstFile_TimeSecDebFic = hourf * \
                            60 * 60 + minf * 60 + secf
                        t = time.strptime("%04i:%02i:%02i:%02i:%02i:%02i" %
                                          (year, month, day, hourf, minf, secf), "%Y:%m:%d:%H:%M:%S")
                        basetime = calendar.timegm(t)
    
                    # TimeSecDebFic = hourf*60*60+minf*60+secf
                    # Relative_TimeSecDebFic = TimeSecDebFic - FirstFile_TimeSecDebFic
    
                    # VecDiff=Relative_TimeSecDebFic-TimeVec
                    # Valdmin = np.amin(abs(VecDiff))
                    # Indmin = np.where(VecDiff==Valdmin)[0][0]
                    # if np.round(Valdmin*1e5)/1e5 != 0:
                        # print "Correction decalage en temps"
                        # FFTdata = np.fft.fft(data)
                        # FFTdata[np.ceil(len(data)/2):] *= 0.
                        # VecFre = np.arange(0,len(data)-1) / (samplerate/ (len(data)-1))
                        # FFTcorr = FFTdata * np.exp(1j * 2. * np.pi * VecFre * Valdmin).T
                        # datac=2. * np.real(np.fft.ifft(FFTcorr))
    
                        # trame[Indmin:len(datac)+Indmin-1]=datac
                    # else:
                        # trame[Indmin:Indmin+len(data)]=data
    
                    # del VecDiff
    
                    if len(trame) % 2 != 0:
                        trame = np.append(trame, 0.)
                    if comp == "Z":
                        tramef_Z[istation] = trame
                    elif comp == "E":
                        tramef_E[istation] = trame
                    elif comp == "N":
                        tramef_N[istation] = trame
    
                    del data, trace, stream, trame
    
        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / goal_sampling_rate
        fe = goal_sampling_rate
        # Calculate the number of slices
        tranches = int(goal_duration * fe / min30)
        # print
        # print '##### ITERATING OVER PAIRS #####'
    
        for pair in pairs:
            orig_pair = pair
            logging.debug('Processing pair: %s' % pair.replace(':', ' vs '))
            tt = time.time()
            # print ">PROCESSING PAIR %s"%pair.replace(':',' vs ')
            station1, station2 = pair.split(':')
            pair = (np.where(stations == station1)
                    [0][0], np.where(stations == station2)[0][0])
    
            s1 = get_station(db, station1.split('.')[0], station1.split('.')[1])
            s2 = get_station(db, station2.split('.')[0], station2.split('.')[1])
    
            X0 = s1.X
            Y0 = s1.Y
            c0 = s1.coordinates
    
            X1 = s2.X
            Y1 = s2.Y
            c1 = s2.coordinates
    
            if c0 == c1:
                if c0 == 'DEG':
                    # print "> I will compute the azimut based on degrees"
                    coordinates = 'DEG'
                else:
                    # print "> I will compute the azimut based on meters"
                    coordinates = 'UTM'
            else:
                # print "> Coordinates type don't match, I will need to compute
                # more stuff !!"
                coordinates = 'MIX'
            # print "X0,Y0 ; X1,Y1:", X0, Y0, X1, Y1
            cplAz = azimuth(coordinates, X0, Y0, X1, Y1)
    
            for components in components_to_compute:
                # we create the two parts of the correlation array checking for the
                # right components :
                if components[0] == "Z":
                    t1 = tramef_Z[pair[0]]
                elif components[0] == "R":
                    t1 = tramef_N[pair[0]] * np.cos(cplAz * np.pi / 180.) + tramef_E[
                        pair[0]] * np.sin(cplAz * np.pi / 180.)
                elif components[0] == "T":
                    t1 = tramef_N[pair[0]] * np.sin(cplAz * np.pi / 180.) - tramef_E[
                        pair[0]] * np.cos(cplAz * np.pi / 180.)
    
                if components[1] == "Z":
                    t2 = tramef_Z[pair[1]]
                elif components[1] == "R":
                    t2 = tramef_N[pair[1]] * np.cos(cplAz * np.pi / 180.) + tramef_E[
                        pair[1]] * np.sin(cplAz * np.pi / 180.)
                elif components[1] == "T":
                    t2 = tramef_N[pair[1]] * np.sin(cplAz * np.pi / 180.) - tramef_E[
                        pair[1]] * np.cos(cplAz * np.pi / 180.)
    
                trames = np.vstack((t1, t2))
                del t1, t2
                ncorr = 0
    
                daycorr = {}
                ndaycorr = {}
                for filterdb in get_filters(db, all=False):
                    filterid = filterdb.ref
                    daycorr[filterid] = np.zeros(get_maxlag_samples(db,))
                    ndaycorr[filterid] = 0
    
                for itranche in range(0, tranches):
                    # print "Avancement: %#2d/%2d"% (itranche+1,tranches)
                    trame2h = trames[:, itranche * min30:(itranche + 1) * min30]
                    rmsmat = np.std(np.abs(trame2h), axis=1)
                    for filterdb in get_filters(db, all=False):
                        filterid = filterdb.ref
                        low = float(filterdb.low)
                        high = float(filterdb.high)
                        rms_threshold = filterdb.rms_threshold
                        # print "Filter Bounds used:", filterid, low, high
                        # Npts = min30
                        # Nc = 2* Npts - 1
                        # Nfft = 2**nextpow2(Nc)
    
                        Nfft = min30
                        if min30 / 2 % 2 != 0:
                            Nfft = min30 + 2
    
                        trames2hWb = np.zeros((2, Nfft), dtype=np.complex)
                        for i, station in enumerate(pair):
                            # print "USING rms threshold = %f" % rms_threshold
                            # logging.debug("rmsmat[i] = %f" % rmsmat[i])
                            if rmsmat[i] > rms_threshold:
                                if windsorizing != 0:
                                    indexes = np.where(
                                        np.abs(trame2h[i]) > (windsorizing * rmsmat[i]))[0]
                                    # clipping at windsorizing*rms
                                    trame2h[i][indexes] = (trame2h[i][indexes] / np.abs(
                                        trame2h[i][indexes])) * windsorizing * rmsmat[i]
    
                                # logging.debug('whiten')
    
                                trames2hWb[i] = whiten(
                                    trame2h[i], Nfft, dt, low, high, plot=False)
                            else:
                                # logging.debug("Station no %d, pas de pretraitement car rms < %f ou NaN"% (i, rms_threshold))
                                trames2hWb[i] = np.zeros(Nfft)
    
                        corr = myCorr(trames2hWb, np.ceil(maxlag / dt), plot=False)
    
                        thisdate = time.strftime(
                            "%Y-%m-%d", time.gmtime(basetime + itranche * min30 / fe))
                        thistime = time.strftime(
                            "%H_%M", time.gmtime(basetime + itranche * min30 / fe))
                        if keep_all:
                            add_corr(db, station1.replace('.', '_'), station2.replace(
                                '.', '_'), filterid, thisdate, thistime, min30 / fe, components, corr, fe)
    
                        if keep_days:
                            if not np.any(np.isnan(corr)) and not np.any(np.isinf(corr)):
                                daycorr[filterid] += corr
                                ndaycorr[filterid] += 1
    
                        del corr, thistime, trames2hWb
    
                if keep_days:
                    try:
                        for filterdb in get_filters(db, all=False):
                            filterid = filterdb.ref
                            corr = daycorr[filterid]
                            ncorr = ndaycorr[filterid]
                            if ncorr > 0:
                                logging.debug(
                                    "Saving daily CCF for filter %02i (stack of %02i CCF)" % (filterid, ncorr))
    
                                # corr /= ncorr
                                thisdate = time.strftime(
                                    "%Y-%m-%d", time.gmtime(basetime))
                                thistime = time.strftime(
                                    "%H_%M", time.gmtime(basetime))
                                add_corr(
                                    db, station1.replace(
                                        '.', '_'), station2.replace('.', '_'), filterid,
                                    thisdate, thistime, min30 / fe, components, corr, fe, day=True, ncorr=ncorr)
                            del corr, ncorr
                    except Exception as e:
                        logging.debug(str(e))
                del trames, daycorr, ndaycorr
            logging.debug("Updating Job")
            update_job(db, goal_day, orig_pair, 'CC', 'D')
    
            logging.debug("Finished processing this pair. It took %.2f seconds" %
                          (time.time() - tt))
        logging.info("Job Finished. It took %.2f seconds" % (time.time() - jt))
    logging.info('*** Finished: Compute CC ***')
    
    
    # EOF
