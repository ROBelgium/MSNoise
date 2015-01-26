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

from msnoise_table_def import Filter, Job, Station, Config, DataAvailability

def get_tech():
    """Returns the current DB technology used (reads from the db.ini file)
    
    Returns
    -------
    tech : int
        The database technology used: 1=sqlite 2=mysql
    """
    tech, hostname, database, username, password = read_database_inifile()
    return tech


def get_engine(inifile=os.path.join(os.getcwd(), 'db.ini')):
    """Returns the a SQLAlchemy Engine
    
    Parameters
    ----------
    inifile : string
        The path to the db.ini file to use. Defaults to os.cwd()
    
    Returns
    -------
    engine : object
    """
    tech, hostname, database, username, password = read_database_inifile(inifile)
    if tech == 1:
        engine = create_engine('sqlite:///%s' % hostname, echo=False)
    else:
        engine = create_engine('mysql://%s:%s@%s/%s' % (username, password,
                                                        hostname, database),
                               echo=False, poolclass=NullPool)
    return engine


def connect(inifile=os.path.join(os.getcwd(), 'db.ini')):
    """Returns the a SQLAlchemy Session
    
    Parameters
    ----------
    inifile : string
        The path to the db.ini file to use. Defaults to os.cwd()
    
    Returns
    -------
    session : object
    """
    engine = get_engine(inifile)
    Session = sessionmaker(bind=engine)
    return Session()


def create_database_inifile(tech, hostname, database, username, password):
    """Creates the db.ini file based on supplied parameters.
    
    Parameters
    ----------
    tech : int
        The database technology used: 1=sqlite 2=mysql
    hostname : string
        The hostname of the server (if tech=2) or the name of the sqlite file
        if tech=1)
    database : string
        The database name
    username : string
    password : string
    
    """
    f = open(os.path.join(os.getcwd(), 'db.ini'), 'w')
    cPickle.dump([tech, hostname, database, username, password], f)
    f.close()


def read_database_inifile(inifile=os.path.join(os.getcwd(), 'db.ini')):
    """Reads the parameters from the db.ini file.
    
    Parameters
    ----------
    inifile : string
        The path to the db.ini file to use. Defaults to os.cwd()
    
    
    Returns
    ----------
    tech : int
        The database technology used: 1=sqlite 2=mysql
    hostname : string
        The hostname of the server (if tech=2) or the name of the sqlite file
        if tech=1)
    database : string
        The database name
    username : string
    password : string
    
    """
    f = open(inifile, 'r')
    tech, hostname, database, username, password = cPickle.load(f)
    f.close()
    return [tech, hostname, database, username, password]


############ CONFIG ############


def get_config(session, name=None, bool=False):
    """Get the value of one or all config bits from the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    name : string, optional
        The name of the config bit to get. If omitted, a dictionnary with all
        config items will be returned
    
    
    Returns
    ----------
    config : string or dict
    
    """
    
    if name:
        config = session.query(Config).filter(Config.name == name).first()
        if config is not None:
            if bool:
                if config.value in [True,'true','Y','y','1',1]:
                    config = True
                else:
                    config = False
            else:
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
    """Update one config bit in the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    name : string
        The name of the config bit to set.
    value : string
        The value of the config bit to set.
    
    """
    
    config = session.query(Config).filter(Config.name == name).first()
    config.value = value
    session.commit()
    return

############ FILTERS ############


def get_filters(session, all=False):
    """Get Filters from the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    all : bool, optional
        Returns all filters from the database if True, or only filters where
        `used` = 1 if False (default)
    
    
    Returns
    ----------
    filters : list of Filter
    
    """
    if all:
        filters = session.query(Filter).all()
    else:
        filters = session.query(Filter).filter(Filter.used == True).all()
    return filters


def update_filter(session, ref, low, high, mwcs_low, mwcs_high,
                  rms_threshold, mwcs_wlen, mwcs_step, used):
    """Updates or Insert a new Filter in the database.
    
    See also
    --------
    :class:`msnoise.msnoise_table_def.Filter`

    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    ref : int
        The id of the Filter in the database
    low : float
        The lower frequency bound of the Whiten function (in Hz)
    high : float
        The upper frequency bound of the Whiten function (in Hz)
    mwcs_low : float
        The lower frequency bound of the linear regression done in MWCS (in Hz)
    mwcs_high : float
        The upper frequency bound of the linear regression done in MWCS (in Hz)
    rms_threshold : float
        Not used anymore
    mwcs_wlen : float
        Window length (in seconds) to perform MWCS
    mwcs_step : float
        Step (in seconds) of the windowing procedure in MWCS
    used : bool
        Is the filter activated for the processing
    
    
    """
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
    """Get Networks from the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    all : bool, optional
        Returns all networks from the database if True, or only networks where
        at least one station has `used` = 1 if False (default)
    
    
    Returns
    ----------
    networks : list of string
    
    """
    if all:
        networks = session.query(Station).group_by(Station.net).all()
    else:
        networks = session.query(Station).filter(Station.used == True).group_by(Station.net)
    return [net.net for net in networks]


def get_stations(session, all=False, net=None):
    """Get Stations from the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    all : bool, optional
        Returns all stations from the database if True, or only stations where
        `used` = 1 if False (default)
    net : string, optional
        if set, limits the stations returned to this network
    
    
    Returns
    ----------
    stations : list of Station
    
    """
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
    """Get one Station from the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    net : string
        The network code
    sta : string
        The station code
    
    
    Returns
    ----------
    station : Station
    
    """
    station = session.query(Station).filter(Station.net == net).filter(Station.sta == sta).first()
    return station


def update_station(session, net, sta, X, Y, altitude, coordinates='UTM', instrument='N/A', used=1):
    """Updates or Insert a new Station in the database.
    
    See also
    --------
    :class:`msnoise.msnoise_table_def.Station`
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
    ref : int
        The Station ID in the database
    net : string
        The network code of the Station
    sta : string
        The station code
    X : float
        The X coordinate of the station
    Y : float
        The Y coordinate of the station
    altitude : float
        The altitude of the station
    coordinates : {'DEG', 'UTM'}
        The coordinates system. DEG is WGS84 latitude/longitude in degrees. 
        UTM is expressed in meters.
    instrument : string
        The instrument code, useful with PAZ correction
    used : bool
        Whether this station must be used in the computations.
    """
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
    """Returns an iterator over all possible station pairs.
    If auto-correlation is configured in the database, returns N*N pairs,
    otherwise returns N*(N-1)/2 pairs.
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
    used : bool
        Select only stations marked used if False (default) or all stations
        present in the database if True
    net : string
        Network code to filter for the pairs.
        
    Returns
    -------
    iterable
        Iterable of Station object pairs.
    """
    stations = get_stations(session, all=False, net=net)
    if get_config(session, name="autocorr") in ['Y', 'y', '1', 1]:
        return itertools.combinations_with_replacement(stations, 2)
    else:
        return itertools.combinations(stations, 2)


def get_interstation_distance(station1, station2, coordinates="DEG"):
    """Returns the distance in km between `station1` and `station2`.
    
    Parameters
    ----------
    station1 : Station object
    station2 : Station object
    coordinates : str
        The coordinate system.
    
    Returns
    -------
    float
        the distance in km.
    """
    
    if coordinates == "DEG":
        dist, azim, bazim = gps2DistAzimuth(station1.Y, station1.X, station2.Y, station2.X)
        return dist /1.e3
    else:
        print "woooooow, UTM system distance not computed"
        return 0

############ DATA AVAILABILITY ############


def update_data_availability(session, net, sta, comp, path, file, starttime, endtime, data_duration, gaps_duration, samplerate):
    """
    Updates a DataAvailability object in the database
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    ref : int
        The Station ID in the database
    net : string
        The network code of the Station
    sta : string
        The station code
    comp : string
        The component (channel)
    path : string
    file : string
    starttime : datetime
        Start time of the file
    endtime : datetime
        End time of the file
    data_duration : float
        Cumulative duration of available data in the file
    gaps_duration : float
        Cumulative duration of gaps in the file
    samplerate : float
        Sample rate of the data in the file (in Hz)

    """
    
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
    """
    Returns the files marked "N"ew or "M"odified in the database
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
        
    Returns
    -------
    files : list of DataAvailability
    """
    
    files = session.query(DataAvailability).filter(DataAvailability.flag != 'A').order_by(DataAvailability.starttime).all()
    return files

def get_data_availability(session, net=None, sta=None, comp=None, starttime=None, endtime=None):
    """
    Returns the DataAvailability objects for specific `net`, `sta`, `starttime` or `endtime`
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
    net : string
        Network code
    sta : string
        Station code
    starttime : datetime
        Start time of the search
    endtime : datetime
        End time of the search

    Returns
    -------
    data : list of DataAvailability
    """
    
    if not starttime:
        data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta).filter(DataAvailability.comp == comp).all()
    elif not net:
        data = session.query(DataAvailability).filter(DataAvailability.starttime <= endtime).filter(DataAvailability.endtime >= starttime).all()
    else:
        data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta).filter(func.DATE(DataAvailability.starttime) <= starttime.date()).filter(func.DATE(DataAvailability.endtime) >= endtime.date()).all()
    return data


def mark_data_availability(session, net, sta, flag):
    """
    Updates the flag of all DataAvailability objects matching `net`.`sta` in the database
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    net : string
        The network code of the Station
    sta : string
        The station code
    flag : {'N', 'M', 'A'}
        Status of the DataAvailability object: New, Modified or Archive
    """
    
    data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta)
    data.update({DataAvailability.flag: flag})
    session.commit()

def count_data_availability_flags(session):
    """
    Count the number of DataAvailability, grouped by `flag`
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    
    Returns
    -------
    data : list of [count, flag] pairs
    """
    
    return session.query(func.count(DataAvailability.flag),DataAvailability.flag).group_by(DataAvailability.flag).all()
    

############ JOBS ############


def update_job(session, day, pair, type, flag, commit=True, returnjob=True):
    """
    Updates or Inserts a new job in the database.
    
    Parameters
    ----------
    
    ref : int
        The Job ID in the database
    day : string
        The day in YYYY-MM-DD format
    pair : string
    type : {'CC', 'DTT'}
    flag : {'T', 'I', 'D'}
        Status of the Job: "T"odo, "I"n Progress, "D"one.
    commit : bool
        Whether to directly commit (True, default) or not (False)
    returnjob : bool
        Return the modified/inserted Job (True, default) or not (False)
        
    Return
    ------
    
    job : Job
        If returnjob is True, returns the modified/inserted Job.
    
    """
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


def massive_insert_job(jobs):
    """
    Routine to use a low level function to insert much faster a list of jobs.
    
    Parameters
    ----------
    
    jobs : list of Job
        A list of Job to insert
    
    """
    engine = get_engine()
    engine.execute(
        Job.__table__.insert(),
        jobs)


def is_next_job(session, flag='T', type='CC'):
    """
    Are there any Job in the database, with flag=`flag` and type=`type`
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
    flag : {'T', 'I', 'D'}
        Status of the Job: "T"odo, "I"n Progress, "D"one.
    type : {'CC', 'DTT'}
    
    Returns
    -------
    bool
        True if at least one Job matches, False otherwise.
    """
    job = session.query(Job).filter(Job.type == type).filter(Job.flag == flag).first()
    if job is None:
        return False
    else:
        return True


def get_next_job(session, flag='T', type='CC'):
    """
    Get the next Job in the database, with flag=`flag` and type=`type`. Jobs
    of the same `type` are grouped per day. This function also sets the flag of
    all selected Jobs to "I"n progress.
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
    flag : {'T', 'I', 'D'}
        Status of the Job: "T"odo, "I"n Progress, "D"one.
    type : {'CC', 'DTT'}
    
    Returns
    -------
    jobs : list of Job
    """
    day = session.query(Job).filter(Job.type == type).filter(Job.flag == flag).order_by(Job.day).first().day
    jobs = session.query(Job).filter(Job.type == type).filter(Job.flag == flag).filter(Job.day == day)
    tmp = jobs.all()
    jobs.update({Job.flag: 'I'})
    session.commit()
    return tmp


def is_dtt_next_job(session, flag='T', type='DTT', ref=False):
    """
    Are there any DTT Job in the database, with flag=`flag` and type=`type`. If
    `ref` is provided, checks if a DTT "REF" job is present.
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
    flag : {'T', 'I', 'D'}
        Status of the Job: "T"odo, "I"n Progress, "D"one.
    type : {'CC', 'DTT'}
    ref : bool
        Whether to check for a REF job (True) or not (False, default)
    
    Returns
    -------
    bool
        True if at least one Job matches, False otherwise.
    """
    if ref:
        job = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.pair == ref).filter(Job.day == 'REF').first()
    else:
        job = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.day != 'REF').first()
    if job is None:
        return False
    else:
        return True


def get_dtt_next_job(session, flag='T', type='DTT'):
    """
    Get the next DTT Job in the database, with flag=`flag` and type=`type`. Jobs
    are then grouped per pair. This function also sets the flag of
    all selected Jobs to "I"n progress.
    
    Parameters
    -----------
    session : object
        A Session object, as obtained using `connect()`
    flag : {'T', 'I', 'D'}
        Status of the Job: "T"odo, "I"n Progress, "D"one.
    type : {'CC', 'DTT'}
    
    Returns
    -------
    pairs : list
        List of station pair names
    days : list of string
        Days of the next DTT jobs
    refs : list of int
        Job IDs (for later being able to update their flag)
    """
    pair = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.day != 'REF').first().pair
    jobs = session.query(Job).filter(Job.flag == flag).filter(Job.type == type).filter(Job.day != 'REF').filter(Job.pair == pair)
    refs, days = zip(*[[job.ref,job.day] for job in jobs.all()])
    jobs.update({Job.flag: 'I'})
    session.commit()
    return pair, days, refs


def reset_jobs(session, jobtype):
    """
    Sets the flag of all `jobtype` Jobs to "T"odo.
    
    Parameters
    ----------
    jobtype : string
        The type to reset: CC or DTT
    """
    jobs = session.query(Job).filter(Job.type == jobtype)
    jobs.update({Job.flag: 'T'})
    session.commit()

def reset_dtt_jobs(session, pair):
    """
    Sets the flag of all DTT Jobs of one `pair` to "T"odo.
    
    Parameters
    ----------
    pair : string
        The pair to update
    """
    
    jobs = session.query(Job).filter(Job.pair == pair).filter(Job.type == "DTT")
    jobs.update({Job.flag: 'T'})
    session.commit()


def get_job_types(session, type='CC'):
    """
    Count the number of DataAvailability of a specific `type`,
    grouped by `flag`.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    type : {'CC', 'DTT'}
        
    
    Returns
    -------
    data : list of [count, flag] pairs
    """
    
    return session.query(func.count(Job.flag),Job.flag).filter(Job.type == type).group_by(Job.flag).all()


############ CORRELATIONS ############


def add_corr(session,station1, station2, filterid, date, time, duration, components, CF, sampling_rate, day=False, ncorr=0):
    """
    Adds a CCF to the data archive on disk.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    station1 : str
        The name of station 1 (formatted NET.STA)
    station2 : str
        The name of station 2 (formatted NET.STA)
    filterid : int
        The ID (ref) of the filter
    date : datetime.date or str
        The date of the CCF
    time : datetime.time or str
        The time of the CCF
    duration : float
        The total duration of the exported CCF
    components : str
        The name of the components used (ZZ, ZR, ...)
    sampling_rate : float
        The sampling rate of the exported CCF
    day : bool
        Whether this function is called to export a daily stack (True) or each
        CCF (when keep_all parameter is set to True in the configuration).
        Defaults to True.
    ncorr : int
        Number of CCF that have been stacked for this CCF.
        
    
    Returns
    -------
    None
    """
    
    output_folder = get_config(session, 'output_folder')
    export_format = get_config(session, 'export_format')
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
            export_mseed(session, path, pair, components, filterid, CF/ncorr, ncorr)
        if sac:
            export_sac(session, path, pair, components, filterid, CF/ncorr, ncorr)

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


def get_results(session, station1, station2, filterid, components, dates, mov_stack = 1, format="stack"):
    export_format = get_config(session, 'export_format')
    if format == "stack":
        stack = np.zeros(get_maxlag_samples(session))
        i = 0
        for date in dates:
            daystack = os.path.join("STACKS", "%02i" % filterid, "%03i_DAYS"%mov_stack, components, "%s_%s"%(station1, station2), str(date))
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
            daystack = os.path.join("STACKS", "%02i"%filterid, "%03i_DAYS"%mov_stack, components, "%s_%s"%(station1, station2), str(date))
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
    """
    Returns the length of the CC functions. Gets the maxlag and sampling rate
    from the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    
    Returns
    -------
    maxlag_samples : int
    """
    
    maxlag = float(get_config(session, 'maxlag'))
    return int(2*maxlag*float(get_config(session, 'cc_sampling_rate'))+1)


def get_components_to_compute(session):
    """
    Returns the components configured in the database.
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    
    Returns
    -------
    components_to_compute : list of string
    """
    
    components_to_compute = []
    for comp in ['ZZ', 'RR', 'TT', 'TR', 'RT', 'ZR', 'RZ', 'TZ', 'ZT']:
        if get_config(session, comp) in ['Y', 'y', '1', 1]:
            components_to_compute.append(comp)
    return components_to_compute


def build_ref_datelist(session):
    """
    Creates a date array for the REF
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    
    Returns
    -------
    start : datetime
        Start of the REF
    end : datetime
        End of the REF
    datelist : list of datetime
        All dates between start and end
    """
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
    """
    Creates a date array for the analyse period
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    
    Returns
    -------
    start : datetime
        Start of the analyse
    end : datetime
        End of the analyse
    datelist : list of datetime
        All dates between start and end
    """
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


def updated_days_for_dates(session, date1, date2, pair, type='CC', interval=datetime.timedelta(days=1), returndays=False):
    """
    Determines if any Job of type=`type` and for pair=`pair`, concerning a date
    between `date1` and `date2` has been modified in the last interval=`interval`.
    
    
    Parameters
    ----------
    session : object
        A Session object, as obtained using `connect()`
    date1 : datetime
        Beginning of the period of interest
    date2 : datetime
        End of the period of interest
    pair : string
        Pair of interest
    type : {'CC', 'DTT'}
    interval : datetime.timedelta
        Interval of time before now to search for updated days
    returndays : bool
        Whether to return a list of days (True) or not (False, default)
    
    Returns
    -------
    return : list or bool
        List of days if returndays is True, True if not. (not clear!)
    """
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
    """
    Returns the azimuth between two coordinate sets.
    
    Parameters
    ----------
    
    coordinates : {'DEG', 'UTM', 'MIX'}
    x0 : float
        X coordinate of station 1
    y0 : float
        Y coordinate of station 1
    x1 : float
        X coordinate of station 2
    y1 : float
        Y coordinate of station 2
    
    
    Returns
    -------
    azimuth : float
        Azimuth in degrees
    """
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
    """
    Returns the next power of 2 of `x`.
    
    Parameters
    ----------
    x : int
    
    Returns
    -------
    nextpow2 : int
    """
    
    return np.ceil(np.log2(np.abs(x)))


def check_and_phase_shift(trace):
    # print trace
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
        # print "Trace is offset by %.6f s from closest delta (%s)" % (dt, direction)
        n = int(2**nextpow2(len(trace.data)))
        FFTdata = scipy.fftpack.fft(trace.data, n=n)
        fftfreq = scipy.fftpack.fftfreq(n,d=trace.stats.delta)
        FFTdata = FFTdata * np.exp(1j * 2. * np.pi * fftfreq * dt)
        trace.data = np.real(scipy.fftpack.ifft(FFTdata, n=n)[:len(trace.data)])
        trace.stats.starttime += dt
        return trace
    else: 
        # print "No Offset"
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
