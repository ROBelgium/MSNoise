# database_tools.py
import os
import logging
import copy
import datetime
import itertools
import cPickle

from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool
import numpy as np
import scipy.fftpack
import math

from obspy.core import Stream, Trace, read
from obspy.signal import cosTaper
from obspy.sac import SacIO
from obspy.core.util import gps2DistAzimuth

from msnoise_table_def import *

def get_tech():
    tech, hostname, database, username, password = read_database_inifile()
    return tech


def connect(inifile=os.path.join(os.getcwd(), 'db.ini')):
    tech, hostname, database, username, password = read_database_inifile(inifile)
    if tech == 1:
        engine = create_engine('sqlite:///%s' % hostname, echo=False)
    else:
        engine = create_engine('mysql://%s:%s@%s/%s' % (username, password,
                                                        hostname, database),
                               echo=False, poolclass=NullPool)
    Session = sessionmaker(bind=engine)
    return Session()


def create_database_inifile(tech, hostname, database, username, password):
    f = open(os.path.join(os.getcwd(), 'db.ini'), 'w')
    cPickle.dump([tech, hostname, database, username, password], f)
    f.close()


def read_database_inifile(inifile=os.path.join(os.getcwd(), 'db.ini')):
    f = open(inifile, 'r')
    tech, hostname, database, username, password = cPickle.load(f)
    f.close()
    return [tech, hostname, database, username, password]
############ CONFIG ############


def get_config(session, name=None):
    if name:
        config = session.query(Config).filter(Config.name == name).first()
        if config is not None:
            config = config.value
        else:
            config = ''
    else:
        config = {}
        configs = session.query(Config).all()
        for c in configs:
            config[c.name] = c.value
    return config


def update_config(session, name, value):
    config = session.query(Config).filter(Config.name == name).first()
    config.value = value
    session.commit()
    return

############ FILTERS ############


def get_filters(session, all=False):
    if all:
        filters = session.query(Filter).all()
    else:
        filters = session.query(Filter).filter(Filter.used == True).all()
    return filters


def update_filter(session, ref, low, high, mwcs_low, mwcs_high,
                  rms_threshold, mwcs_wlen, mwcs_step, used):
    filter = session.query(Filter).filter(Filter.ref == ref).first()
    if filter is None:
        filter = Filter(low, mwcs_low, high, mwcs_high, rms_threshold,
                        mwcs_wlen, mwcs_step, used)
        session.add(filter)
    else:
        filter.low = low
        filter.high = high
        filter.mwcs_low = mwcs_low
        filter.mwcs_high = mwcs_high
        filter.rms_threshold = rms_threshold
        filter.mwcs_wlen = mwcs_wlen
        filter.mwcs_step = mwcs_step
        filter.used = used
    session.commit()
    return

############ NETWORK AND STATION ############


def get_networks(session, all=False):
    if all:
        networks = session.query(Station).group_by(Station.net).all()
    else:
        networks = session.query(Station).filter(Station.used == True).group_by(Station.net)
    return [net.net for net in networks]


def get_stations(session, all=False, net=None):
    if all:
        if net is not None:
            stations = session.query(Station).filter(Station.net == net).order_by(Station.net).order_by(Station.sta)
        else:
            stations = session.query(Station).order_by(Station.net).order_by(Station.sta).all()
    else:
        stations = session.query(Station).filter(Station.used == True).order_by(Station.net).order_by(Station.sta)
        if net is not None:
            stations = stations.filter(Station.net == net).order_by(Station.net).order_by(Station.sta)
    return stations


def get_station(session, net, sta):
    station = session.query(Station).filter(Station.net == net).filter(Station.sta == sta).first()
    return station


def update_station(session, net, sta, X, Y, altitude, coordinates='UTM', instrument='N/A', used=1):
    station = session.query(Station).filter(Station.net == net).filter(Station.sta == sta).first()
    if station is None:
        station = Station(net, sta, X, Y, altitude, coordinates, instrument, used)
        session.add(station)
    else:
        station.X = X
        station.Y = Y
        station.altitude = altitude
        station.coordinates = coordinates
        station.instrument = instrument
        station.used = used
    session.commit()
    return True


def get_station_pairs(session, used=None, net=None):
    stations = get_stations(session, all=False, net=net)
    if get_config(session, name="autocorr") in ['Y', 'y', '1', 1]:
        return itertools.combinations_with_replacement(stations, 2)
    else:
        return itertools.combinations(stations, 2)

############ DATA AVAILABILITY ############


def update_data_availability(session, net, sta, comp, path, file, starttime, endtime, data_duration, gaps_duration, samplerate):
    data = session.query(DataAvailability).filter(DataAvailability.file == file).first()
    if data is None:
        flag = "N"
        data = DataAvailability(net, sta, comp, path, file, starttime, endtime, data_duration, gaps_duration, samplerate, flag)
        session.add(data)
        toreturn = True
    else:
        modified = False
        for item in ['net', 'sta', 'comp', 'path', 'starttime', 'endtime',
                     'data_duration', 'gaps_duration', 'samplerate']:
            if eval("data.%s != %s" % (item, item)):
                modified = True
                break
        if modified:
            data.net = net
            data.sta = sta
            data.comp = comp
            data.path = path
            data.starttime = starttime
            data.endtime = endtime
            data.data_duration = data_duration
            data.gaps_duration = gaps_duration
            data.samplerate = samplerate
            data.flag = "M"
        toreturn = False
    session.commit()
    return toreturn


def get_new_files(session):
    files = session.query(DataAvailability).filter(DataAvailability.flag != 'A').order_by(DataAvailability.starttime).all()
    return files

def get_data_availability(session, net=None, sta=None, comp=None, starttime=None, endtime=None):
    if not starttime:
        data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta).filter(DataAvailability.comp == comp).all()
    elif not net:
        data = session.query(DataAvailability).filter(DataAvailability.starttime <= endtime).filter(DataAvailability.endtime >= starttime).all()
    else:
        data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta).filter(func.DATE(DataAvailability.starttime) <= starttime.date()).filter(func.DATE(DataAvailability.endtime) >= endtime.date()).all()
    return data


def mark_data_availability(session,net,sta,flag):
    data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta).all()
    for d in data:
        d.flag = flag
    session.commit()

def count_data_availability_flags(session):
    return session.query(func.count(DataAvailability.flag),DataAvailability.flag).group_by(DataAvailability.flag).all()
    

############ JOBS ############


def update_job(session, day, pair, type, flag, commit=True, returnjob=True):
    job = session.query(Job).filter(Job.day == day).filter(Job.pair == pair).filter(Job.type == type).first()
    if job is None:
        job = Job(day, pair, type, 'T')
        if commit:
            session.add(job)
    else:
        job.flag = flag
        job.lastmod = datetime.datetime.utcnow()
    if commit:
        session.commit()
    if returnjob:
        return job


def is_next_job(session, flag='T', type='CC'):
    job = session.query(Job).filter(Job.type == type).filter(Job.flag == flag).first()
    if job is None:
        return False
    else:
        return True


def get_next_job(session, flag='T', type='CC'):
    day = session.query(Job).filter(Job.type == type).filter(Job.flag == flag).order_by(Job.day).first().day
    # print day
    jobs = session.query(Job).filter(Job.type == type).filter(Job.flag == flag).filter(Job.day == day).all()
    return jobs


def is_dtt_next_job(session, flag='T', type='DTT', ref=False):
    if ref:
        job = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.pair == ref).filter(Job.day == 'REF').first()
    else:
        job = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.day != 'REF').first()
    if job is None:
        return False
    else:
        return True


def get_dtt_next_job(session, flag='T', type='DTT'):
    pair = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.day != 'REF').first().pair
    jobs = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.day != 'REF').filter(Job.pair == pair).all()
    refs = [job.ref for job in jobs]
    days = [job.day for job in jobs]
    for job in jobs:
        job.flag = 'I'
    session.commit()
    return pair, days, refs


def reset_dtt_jobs(session, pair):
    jobs = session.query(Job).filter(Job.pair == pair).filter(Job.type == "DTT").all()
    for job in jobs:
        job.flag = "T"
    session.commit()


def get_job_types(session, type='CC'):
    return session.query(func.count(Job.flag),Job.flag).filter(Job.type == type).group_by(Job.flag).all()


############ CORRELATIONS ############


def add_corr(db,station1, station2, filterid, date, time, duration, components, CF, sampling_rate, day=False, ncorr=0):
    output_folder = get_config(db, 'output_folder')
    export_format = get_config(db, 'export_format')
    if export_format == "BOTH":
        mseed = True
        sac = True
    elif export_format == "SAC":
        mseed = False
        sac = True
    elif export_format == "MSEED":
        mseed = True
        sac = False

    if day:
        path = os.path.join("STACKS", "%02i" % filterid, "001_DAYS", components, "%s_%s" % (station1, station2), str(date))
        pair = "%s:%s" % (station1, station2)
        if mseed:
            export_mseed(db, path, pair, components, filterid, CF/ncorr, ncorr)
        if sac:
            export_sac(db, path, pair, components, filterid, CF/ncorr, ncorr)

    else:
        file = '%s.cc' % time
        path = os.path.join(output_folder, "%02i" % filterid, station1, station2, components, date)
        if not os.path.isdir(path):
            os.makedirs(path)

        t = Trace()
        t.data = CF
        t.stats.sampling_rate = sampling_rate
        t.stats.starttime = -float(get_config(db, 'maxlag'))
        t.stats.components = components
        # if ncorr != 0:
            # t.stats.location = "%02i"%ncorr
        st = Stream(traces=[t, ])
        st.write(os.path.join(path, file), format='mseed')
        del t, st


def export_sac(db, filename, pair, components, filterid, corr, ncorr=0, sac_format=None, maxlag=None, cc_sampling_rate=None):
    if sac_format is None:
        sac_format = get_config(db, "sac_format")
    if maxlag is None:
        maxlag = float(get_config(db, "maxlag"))
    if cc_sampling_rate is None:
        cc_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    try:
        os.makedirs(os.path.split(filename)[0])
    except:
        pass
    filename += ".SAC"
    mytrace = Trace(data=corr)
    mytrace.stats['station'] = pair
    mytrace.stats['sampling_rate'] = cc_sampling_rate

    st = Stream(traces=[mytrace, ])
    st.write(filename, format='SAC')
    tr = SacIO(filename)
    if sac_format == "doublets":
        tr.SetHvalue('A', 120)
    else:
        tr.SetHvalue('B', -maxlag)
        tr.SetHvalue('DEPMIN', np.min(corr))
        tr.SetHvalue('DEPMAX', np.max(corr))
        tr.SetHvalue('DEPMEN', np.mean(corr))
        tr.SetHvalue('SCALE', 1)
        tr.SetHvalue('NPTS', len(corr))
    tr.WriteSacBinary(filename)
    del st, tr
    return


def export_mseed(db, filename, pair, components, filterid, corr, ncorr=0, maxlag=None, cc_sampling_rate=None):
    try:
        os.makedirs(os.path.split(filename)[0])
    except:
        pass
    filename += ".MSEED"
    
    if maxlag is None:
        maxlag = float(get_config(db, "maxlag"))
    if cc_sampling_rate is None:
        cc_sampling_rate = float(get_config(db, "cc_sampling_rate"))

    mytrace = Trace(data=corr)
    mytrace.stats['station'] = pair[:11]
    mytrace.stats['sampling_rate'] = cc_sampling_rate
    mytrace.stats['start_time'] = -maxlag
    mytrace.stats['location'] = "%02i" % ncorr

    st = Stream(traces=[mytrace, ])
    st.write(filename, format='MSEED')
    del st
    return


def get_results(session, station1, station2, filterid, components, dates, format="stack"):
    export_format = get_config(session, 'export_format')
    if format == "stack":
        stack = np.zeros(get_maxlag_samples(session))
        i = 0
        for date in dates:
            daystack = os.path.join("STACKS", "%02i" % filterid, "001_DAYS", components, "%s_%s"%(station1, station2), str(date))
            # logging.debug('reading: %s' % daystack)
            if export_format == "BOTH":
                daystack += ".MSEED"
            elif export_format == "SAC":
                daystack += ".SAC"
            elif export_format == "MSEED":
                daystack += ".MSEED"
            if os.path.isfile(daystack):
                st = read(daystack)
                if not np.any(np.isnan(st[0].data)) and not np.any(np.isinf(st[0].data)):
                    stack += st[0].data
                    i += 1
                else:
                    # print "NaN ! or Inf"
                    pass
        if i > 0:
            return i, stack / i
        else:
            return 0, None

    elif format=="matrix":
        stack = np.zeros((len(dates), get_maxlag_samples(session)))
        i = 0
        for j, date in enumerate(dates):
            daystack = os.path.join("STACKS", "%02i"%filterid, "001_DAYS", components, "%s_%s"%(station1, station2), str(date))
            # logging.debug('reading: %s' % daystack)
            if export_format == "BOTH":
                daystack += ".MSEED"
            elif export_format == "SAC":
                daystack += ".SAC"
            elif export_format == "MSEED":
                daystack += ".MSEED"
            if os.path.isfile(daystack):
                st = read(daystack)
                stack[j] = st[0].data
                i += 1
            else:
                stack[j] *= np.nan
        return i, stack


############ session MISC ############


def get_maxlag_samples(session):
    maxlag = float(get_config(session, 'maxlag'))
    return int(2*maxlag*float(get_config(session, 'cc_sampling_rate'))+1)


def get_components_to_compute(session):
    components_to_compute = []
    for comp in ['ZZ', 'RR', 'TT', 'TR', 'RT', 'ZR', 'RZ', 'TZ', 'ZT']:
        if get_config(session, comp) in ['Y', 'y', '1', 1]:
            components_to_compute.append(comp)
    return components_to_compute


def build_ref_datelist(session):
    begin = get_config(session, "ref_begin")
    end = get_config(session, "ref_end")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d').date()
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()

    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]


def build_movstack_datelist(session):
    begin = get_config(session, "startdate")
    end = get_config(session, "enddate")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d').date()
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()

    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]


def build_daystack_datelist(session):
    refstart, refend, refdates = build_ref_datelist(session)
    movstart, movend, movdates = build_movstack_datelist(session)
    start = min(refstart, movstart)
    end = max(refend, movend)
    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]


def updated_days_for_dates(session, date1, date2, pair, type='CC', interval=datetime.timedelta(days=1), returndays=False):
    lastmod = datetime.datetime.now() - interval
    if pair == '%':
        days = session.query(Job).filter(Job.day >= date1).filter(Job.day <= date2).filter(Job.type == type).filter(Job.lastmod >= lastmod).group_by(Job.day).order_by(Job.day).all()
    else:
        days = session.query(Job).filter(Job.pair == pair).filter(Job.day >= date1).filter(Job.day <= date2).filter(Job.type == type).filter(Job.lastmod >= lastmod).group_by(Job.day).order_by(Job.day).all()
    logging.debug('Found %03i updated days' % len(days))
    if returndays and len(days) != 0:
        return [datetime.datetime.strptime(day.day,'%Y-%m-%d').date() for day in days] ## RETURN DATE LIST !!!
    elif returndays and len(days) == 0:
        return []
    else:
        return True

############ MISCS ############


def azimuth(coordinates, x0, y0, x1, y1):
    if coordinates == "DEG":
        dist, azim, bazim = gps2DistAzimuth(y0, x0, y1, x1)
        # print dist, azim, bazi
        return azim
    elif coordinates == 'UTM':
        azim = 90. - np.arctan2((y1 - y0), (x1 - x0)) * 180. / np.pi
        # print azim
        return azim
    else:
        print "woooooow, please consider having a single coordinate system for all stations"
        return 0


def nextpow2(x):
    return np.ceil(np.log2(np.abs(x)))


def check_and_phase_shift(trace):
    print trace
    taper_length = 20.0
    if trace.stats.npts < 4 * taper_length*trace.stats.sampling_rate:
        trace.data = np.zeros(trace.stats.npts)
        return trace
    
    dt = np.mod(trace.stats.starttime.datetime.microsecond*1.0e-6, trace.stats.delta)
    if (trace.stats.delta -dt) <= np.finfo(float).eps:
        dt = 0
    if dt != 0:
        if dt <= (trace.stats.delta / 2.):
            dt = -dt
            direction = "left"
        else:
            dt = (trace.stats.delta - dt)
            direction = "right"
        trace.detrend(type="demean")
        trace.detrend(type="simple")
        taper_1s = taper_length * float(trace.stats.sampling_rate) / trace.stats.npts
        cp = cosTaper(trace.stats.npts, taper_1s)
        trace.data *= cp
        print "Trace is offset by %.6f s from closest delta (%s)" % (dt, direction)
        n = int(2**nextpow2(len(trace.data)))
        FFTdata = scipy.fftpack.fft(trace.data, n=n)
        fftfreq = scipy.fftpack.fftfreq(n,d=trace.stats.delta)
        FFTdata = FFTdata * np.exp(1j * 2. * np.pi * fftfreq * dt)
        trace.data = np.real(scipy.fftpack.ifft(FFTdata, n=n)[:len(trace.data)])
        trace.stats.starttime += dt
        return trace
    else: 
        print "No Offset"
        return trace
    

def getGaps(stream, min_gap=None, max_gap=None):
    # Create shallow copy of the traces to be able to sort them later on.
    copied_traces = copy.copy(stream.traces)
    stream.sort()
    gap_list = []
    for _i in xrange(len(stream.traces) - 1):
        # skip traces with different network, station, location or channel
        if stream.traces[_i].id != stream.traces[_i + 1].id:
            continue
        # different sampling rates should always result in a gap or overlap
        if stream.traces[_i].stats.delta == stream.traces[_i + 1].stats.delta:
            flag = True
        else:
            flag = False
        stats = stream.traces[_i].stats
        stime = stats['endtime']
        etime = stream.traces[_i + 1].stats['starttime']
        delta = etime.timestamp - stime.timestamp
        # Check that any overlap is not larger than the trace coverage
        if delta < 0:
            temp = stream.traces[_i + 1].stats['endtime'].timestamp - \
                etime.timestamp
            if (delta * -1) > temp:
                delta = -1 * temp
        # Check gap/overlap criteria
        if min_gap and delta < min_gap:
            continue
        if max_gap and delta > max_gap:
            continue
        # Number of missing samples
        nsamples = int(round(math.fabs(delta) * stats['sampling_rate']))
        # skip if is equal to delta (1 / sampling rate)
        if flag and nsamples == 1:
            continue
        elif delta > 0:
            nsamples -= 1
        else:
            nsamples += 1
        gap_list.append([_i,_i+1, 
                        stats['network'], stats['station'],
                        stats['location'], stats['channel'],
                        stime, etime, delta, nsamples])
    # Set the original traces to not alter the stream object.
    stream.traces = copied_traces
    return gap_list


################## TEST
if __name__ == "__main__":
    s = connect()
    for filter in get_filters(s, False):
        print filter.ref, filter.low, filter.high

    print get_networks(s)

    for station in get_stations(s, False, net='BE'):
        print station.net, station.sta

    print get_config(s)
    print get_config(s, 'data_folder')
