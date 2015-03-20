import os
import logging
import copy
import datetime
import itertools
import cPickle
import math

from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool
import numpy as np
import pandas as pd
import scipy.fftpack
from obspy.core import Stream, Trace, read, AttribDict
from obspy.signal import cosTaper

from obspy.core.util import gps2DistAzimuth

from msnoise_table_def import Filter, Job, Station, Config, DataAvailability


def get_tech():
    """Returns the current DB technology used (reads from the db.ini file)

    :rtype: int
    :returns: The database technology used: 1=sqlite 2=mysql
    """
    tech, hostname, database, username, password = read_database_inifile()
    return tech


def get_engine(inifile=None):
    """Returns the a SQLAlchemy Engine

    :type inifile: str
    :param inifile: The path to the db.ini file to use. Defaults to os.cwd() +
        db.ini

    :rtype: :class:`sqlalchemy.engine.Engine`
    :returns: An :class:`~sqlalchemy.engine.Engine` Object
    """
    if not inifile:
        inifile = os.path.join(os.getcwd(), 'db.ini')
    tech, hostname, database, user, passwd = read_database_inifile(inifile)
    if tech == 1:
        engine = create_engine('sqlite:///%s' % hostname, echo=False)
    else:
        engine = create_engine('mysql://%s:%s@%s/%s' % (user, passwd, hostname,
                                                        database),
                               echo=False, poolclass=NullPool)
    return engine


def connect(inifile=None):
    """Establishes a connection to the database and returns a Session object.

    :type inifile: string
    :param inifile: The path to the db.ini file to use. Defaults to os.cwd() +
        db.ini

    :rtype: :class:`sqlalchemy.orm.session.Session`
    :returns: A :class:`~sqlalchemy.orm.session.Session` object, needed for
        many of the other API methods.
    """
    if not inifile:
        inifile = os.path.join(os.getcwd(), 'db.ini')
    engine = get_engine(inifile)
    Session = sessionmaker(bind=engine)
    return Session()


def create_database_inifile(tech, hostname, database, username, password):
    """Creates the db.ini file based on supplied parameters.

    :type tech: int
    :param tech: The database technology used: 1=sqlite 2=mysql
    :type hostname: string
    :param hostname: The hostname of the server (if tech=2) or the name of the
        sqlite file if tech=1)
    :type database: string
    :param database: The database name
    :type username: string
    :param username: The user name
    :type password: string
    :param password: The password of `user`

    :return: None
    """
    f = open(os.path.join(os.getcwd(), 'db.ini'), 'w')
    cPickle.dump([tech, hostname, database, username, password], f)
    f.close()


def read_database_inifile(inifile=None):
    """Reads the parameters from the db.ini file.

    :type inifile: string
    :param inifile: The path to the db.ini file to use. Defaults to os.cwd() +
        db.ini

    :rtype: tuple
    :returns: tech, hostname, database, username, password
    """
    if not inifile:
        inifile = os.path.join(os.getcwd(), 'db.ini')

    f = open(inifile, 'r')
    tech, hostname, database, username, password = cPickle.load(f)
    f.close()
    return [tech, hostname, database, username, password]


# CONFIG


def get_config(session, name=None, isbool=False):
    """Get the value of one or all config bits from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type name: str
    :param name: The name of the config bit to get. If omitted, a dictionnary
        with all config items will be returned
    :type isbool: bool
    :param isbool: if True, returns True/False for config `name`. Defaults to
        False

    :rtype: str, bool or dict
    :returns: the value for `name` or a dict of all config values
    """

    if name:
        config = session.query(Config).filter(Config.name == name).first()
        if config is not None:
            if isbool:
                if config.value in [True, 'True', 'true', 'Y', 'y', '1', 1]:
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

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type name: str
    :param name: The name of the config bit to set.

    :type value: str
    :param value: The value of parameter `name`

    """

    config = session.query(Config).filter(Config.name == name).first()
    config.value = value
    session.commit()
    return

# FILTERS PART


def get_filters(session, all=False):
    """Get Filters from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type all: bool
    :param all: Returns all filters from the database if True, or only filters
        where `used` = 1 if False (default)

    :rtype: list of :class:`~msnoise.msnoise_table_def.Filter`
    :returns: a list of Filter
    """

    if all:
        filters = session.query(Filter).all()
    else:
        filters = session.query(Filter).filter(Filter.used == True).all()
    return filters


def update_filter(session, ref, low, mwcs_low, high, mwcs_high,
                  rms_threshold, mwcs_wlen, mwcs_step, used):
    """Updates or Insert a new Filter in the database.

    .. seealso:: :class:`msnoise.msnoise_table_def.Filter`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type ref: int
    :param ref: The id of the Filter in the database
    :type low: float
    :param low: The lower frequency bound of the Whiten function (in Hz)
    :type high: float
    :param high: The upper frequency bound of the Whiten function (in Hz)
    :type mwcs_low: float
    :type mwcs_low: The lower frequency bound of the linear regression done in
        MWCS (in Hz)
    :type mwcs_high: float
    :type mwcs_high: The upper frequency bound of the linear regression done in
        WCS (in Hz)
    :type rms_threshold: float
    :param rms_threshold: Not used anymore
    :type mwcs_wlen: float
    :param mwcs_wlen: Window length (in seconds) to perform MWCS
    :type mwcs_step: float
    :param mwcs_step: Step (in seconds) of the windowing procedure in MWCS
    :type used: bool
    :param used: Is the filter activated for the processing
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

# NETWORK AND STATION


def get_networks(session, all=False):
    """Get Networks from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type all: bool
    :param all: Returns all networks from the database if True, or only
        networks at least one station has `used` = 1 if False (default)

    :rtype: list of str
    :returns: a list of network codes
    """
    if all:
        networks = session.query(Station).group_by(Station.net).all()
    else:
        networks = session.query(Station).filter(Station.used == True).group_by(Station.net)
    return [net.net for net in networks]


def get_stations(session, all=False, net=None):
    """Get Stations from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type all: bool
    :param all: Returns all stations from the database if True, or only
        stations where `used` = 1 if False (default)
    :type net: str
    :param net: if set, limits the stations returned to this network

    :rtype: list of :class:`msnoise.msnoise_table_def.Station`
    :returns: list of :class:`~msnoise.msnoise_table_def.Station`
    """
    q = session.query(Station)
    if all:
        if net is not None:
            stations = q.filter(Station.net == net).order_by(Station.net).order_by(Station.sta)
        else:
            stations = q.order_by(Station.net).order_by(Station.sta).all()
    else:
        stations = q.filter(Station.used == True).order_by(Station.net).order_by(Station.sta)
        if net is not None:
            stations = stations.filter(Station.net == net).order_by(Station.net).order_by(Station.sta)
    return stations


def get_station(session, net, sta):
    """Get one Station from the database.
    
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: the network code
    :type sta: str
    :param sta: the station code

    :rtype: :class:`msnoise.msnoise_table_def.Station`
    :returns: a :class:`~msnoise.msnoise_table_def.Station` Object

    """
    station = session.query(Station).filter(Station.net == net).filter(Station.sta == sta).first()
    return station


def update_station(session, net, sta, X, Y, altitude, coordinates='UTM',
                   instrument='N/A', used=1):
    """Updates or Insert a new Station in the database.

    .. seealso :: :class:`msnoise.msnoise_table_def.Station`
    
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: The network code of the Station
    :type sta: str
    :param sta: The station code
    :type X: float
    :param X: The X coordinate of the station
    :type Y: float
    :param Y: The Y coordinate of the station
    :type altitude: float
    :param altitude: The altitude of the station
    :type coordinates: str
    :param coordinates: The coordinates system. "DEG" is WGS84 latitude/
        longitude in degrees. "UTM" is expressed in meters.
    :type instrument: str
    :param instrument: The instrument code, useful with PAZ correction
    :type used: bool
    :param used: Whether this station must be used in the computations.
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

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type used: bool, int
    :param used: Select only stations marked used if False (default) or all
        stations present in the database if True
    :type net: str
    :param net: Network code to filter for the pairs.

    :rtype: iterable
    :returns: An iterable of :class:`~msnoise.msnoise_table_def.Station` object
        pairs
    """
    stations = get_stations(session, all=False, net=net)
    if get_config(session, name="autocorr", isbool=True):
        return itertools.combinations_with_replacement(stations, 2)
    else:
        return itertools.combinations(stations, 2)


def get_interstation_distance(station1, station2, coordinates="DEG"):
    """Returns the distance in km between `station1` and `station2`.

    .. warning:: Currently the stations coordinates system have to be the same!

    :type station1: :class:`~msnoise.msnoise_table_def.Station`
    :param station1: A Station object
    :type station2: :class:`~msnoise.msnoise_table_def.Station`
    :param station2: A Station object
    :type coordinates: str
    :param coordinates: The coordinates system. "DEG" is WGS84 latitude/
        longitude in degrees. "UTM" is expressed in meters.



    :rtype: float
    :returns: The interstation distance in km
    """

    if coordinates == "DEG":
        dist, azim, bazim = gps2DistAzimuth(station1.Y, station1.X,
                                            station2.Y, station2.X)
        return dist / 1.e3
    else:
        dist = np.hypot(float(station1.X - station2.X),
                        float(station1.Y - station2.Y)) / 1.e3
        return dist


# DATA AVAILABILITY


def update_data_availability(session, net, sta, comp, path, file, starttime,
                             endtime, data_duration, gaps_duration,
                             samplerate):
    """
    Updates a DataAvailability object in the database

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: The network code of the Station
    :type sta: str
    :param sta: The station code
    :type comp: str
    :param comp: The component (channel)
    :type path: str
    :param path: The full path to the folder containing the file
    :type file: str
    :param file: The name of the file
    :type starttime: datetime.datetime
    :param starttime: Start time of the file
    :type endtime: datetime.datetime
    :param endtime: End time of the file
    :type data_duration: float
    :param data_duration: Cumulative duration of available data in the file
    :type gaps_duration: float
    :param gaps_duration: Cumulative duration of gaps in the file
    :type samplerate: float
    :param samplerate: Sample rate of the data in the file (in Hz)
    """

    data = session.query(DataAvailability).filter(DataAvailability.file == file).first()
    if data is None:
        flag = "N"
        data = DataAvailability(net, sta, comp, path, file, starttime, endtime,
                                data_duration, gaps_duration, samplerate, flag)
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

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list
    :returns: list of :class:`~msnoise.msnoise_table_def.DataAvailability`
    """

    files = session.query(DataAvailability).filter(DataAvailability.flag != 'A').order_by(DataAvailability.starttime).all()
    return files


def get_data_availability(session, net=None, sta=None, comp=None,
                          starttime=None, endtime=None):
    """
    Returns the :class:`~msnoise.msnoise_table_def.DataAvailability` objects
    for specific `net`, `sta`, `starttime` or `endtime`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: Network code
    :type sta: str
    :param sta: Station code
    :type starttime: datetime.datetime, datetime.date
    :param starttime: Start time of the search
    :type endtime: datetime.datetime, datetime.date
    :param endtime: End time of the search

    :rtype: list
    :returns: list of :class:`~msnoise.msnoise_table_def.DataAvailability`
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
    Updates the flag of all
    :class:`~msnoise.msnoise_table_def.DataAvailability` objects matching
    `net.sta` in the database

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type net: str
    :param net: Network code
    :type sta: str
    :param sta: Station code
    :type flag: str
    :param flag: Status of the DataAvailability object: New, Modified or
        Archive. Values accepted are {'N', 'M', 'A'}
    """

    data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta)
    data.update({DataAvailability.flag: flag})
    session.commit()


def count_data_availability_flags(session):
    """
    Count the number of :class:`~msnoise.msnoise_table_def.DataAvailability`,
    grouped by `flag`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list
    :returns: list of [count, flag] pairs
    """

    return session.query(func.count(DataAvailability.flag),DataAvailability.flag).group_by(DataAvailability.flag).all()


# Jobs


def update_job(session, day, pair, jobtype, flag, commit=True, returnjob=True):
    """
    Updates or Inserts a new :class:`~msnoise.msnoise_table_def.Job` in the
    database.

    :type day: str
    :param day: The day in YYYY-MM-DD format
    :type pair: str
    :param pair: the name of the pair (EXAMPLE?)
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.
    :type commit: bool
    :param commit: Whether to directly commit (True, default) or not (False)
    :type returnjob: bool
    :param returnjob: Return the modified/inserted Job (True, default) or not
        (False)


    :rtype: :class:`~msnoise.msnoise_table_def.Job` or None
    :returns: If returnjob is True, returns the modified/inserted Job.
    """
    
    job = session.query(Job).filter(Job.day == day).filter(Job.pair == pair).filter(Job.jobtype == jobtype).first()
    if job is None:
        job = Job(day, pair, jobtype, 'T')
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
    Routine to use a low level function to insert much faster a list of
    :class:`~msnoise.msnoise_table_def.Job`. This method uses the Engine
    directly, no need to pass a Session object.

    :type jobs: list
    :param jobs: a list of :class:`~msnoise.msnoise_table_def.Job` to insert.
    """
    engine = get_engine()
    engine.execute(
        Job.__table__.insert(),
        jobs)


def is_next_job(session, flag='T', jobtype='CC'):
    """
    Are there any :class:`~msnoise.msnoise_table_def.Job` in the database,
    with flag=`flag` and jobtype=`type`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.

    :rtype: bool
    :returns: True if at least one :class:`~msnoise.msnoise_table_def.Job`
        matches, False otherwise.
    """
    job = session.query(Job).filter(Job.jobtype == jobtype).filter(Job.flag == flag).first()
    if job is None:
        return False
    else:
        return True


def get_next_job(session, flag='T', jobtype='CC'):
    """
    Get the next :class:`~msnoise.msnoise_table_def.Job` in the database,
    with flag=`flag` and jobtype=`jobtype`. Jobs of the same `type` are grouped per
    day. This function also sets the flag of all selected Jobs to "I"n progress.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.

    :rtype: list
    :returns: list of :class:`~msnoise.msnoise_table_def.Job`
    """
    day = session.query(Job).filter(Job.jobtype == jobtype).filter(Job.flag == flag).order_by(Job.day).first().day
    jobs = session.query(Job).filter(Job.jobtype == jobtype).filter(Job.flag == flag).filter(Job.day == day)
    tmp = jobs.all()
    jobs.update({Job.flag: 'I'})
    session.commit()
    return tmp


def is_dtt_next_job(session, flag='T', jobtype='DTT', ref=False):
    """
    Are there any DTT :class:`~msnoise.msnoise_table_def.Job` in the database,
    with flag=`flag` and jobtype=`jobtype`. If `ref` is provided, checks if a DTT
    "REF" job is present.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.
    :type ref: bool
    :param ref: Whether to check for a REF job (True) or not (False, default)

    :rtype: bool
    :returns: True if at least one Job matches, False otherwise.
    """
    q = session.query(Job.ref).filter(Job.flag == flag).filter(Job.jobtype == jobtype)
    if ref:
        job = q.filter(Job.pair == ref).filter(Job.day == 'REF').count()
    else:
        job = q.filter(Job.day != 'REF').count()
    if job == 0:
        return False
    else:
        return True


def get_dtt_next_job(session, flag='T', jobtype='DTT'):
    """
    Get the next DTT :class:`~msnoise.msnoise_table_def.Job` in the database,
    with flag=`flag` and jobtype=`jobtype`. Jobs are then grouped per station pair.
    This function also sets the flag of all selected Jobs to "I"n progress.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.

    :rtype: tuple
    :returns: (pairs, days, refs):
        List of station pair names -
        Days of the next DTT jobs -
        Job IDs (for later being able to update their flag).
    """
    pair = session.query(Job).filter(Job.flag == flag).filter(Job.jobtype == jobtype).filter(Job.day != 'REF').first().pair
    jobs = session.query(Job.ref, Job.day).filter(Job.flag == flag).filter(Job.jobtype == jobtype).filter(Job.day != 'REF').filter(Job.pair == pair)
    tmp = list(jobs)
    jobs.update({Job.flag: 'I'}, synchronize_session=False)
    session.commit()
    refs, days = zip(*[[job.ref,job.day] for job in tmp])
    return pair, days, refs


def reset_jobs(session, jobtype, alljobs=False):
    """
    Sets the flag of all `jobtype` Jobs to "T"odo.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type alljobs: bool
    :param alljobs: If True, resets all jobs. If False (default), only resets
        jobs "I"n progress.
    """
    jobs = session.query(Job).filter(Job.jobtype == jobtype)
    if not alljobs:
        jobs = jobs.filter(Job.flag == "I")
    jobs.update({Job.flag: 'T'})
    session.commit()


def reset_dtt_jobs(session, pair):
    """
    Sets the flag of all DTT Jobs of one `pair` to "T"odo.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type pair: str
    :param pair: The pair to update
    """

    jobs = session.query(Job).filter(Job.pair == pair).filter(Job.jobtype == "DTT")
    jobs.update({Job.flag: 'T'})
    session.commit()


def get_job_types(session, jobtype='CC'):
    """
    Count the number of Jobs of a specific `type`,
    grouped by `flag`.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?

    :rtype: list
    :returns: list of [count, flag] pairs
    """

    return session.query(func.count(Job.flag), Job.flag).filter(Job.jobtype == jobtype).group_by(Job.flag).all()


# CORRELATIONS


def export_allcorr(session, ccfid, data):
    output_folder = get_config(session, 'output_folder')
    station1, station2, filterid, components, date = ccfid.split('_')

    path = os.path.join(output_folder, "%02i" % int(filterid),
                        station1, station2, components)
    if not os.path.isdir(path):
        os.makedirs(path)

    df = pd.DataFrame().from_dict(data).T
    df.columns = get_t_axis(session)
    df.to_hdf(os.path.join(path, date+'.h5'), 'data')
    del df
    return


def add_corr(session, station1, station2, filterid, date, time, duration, components, CF, sampling_rate, day=False, ncorr=0):
    """
    Adds a CCF to the data archive on disk.
    
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET.STA)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA)
    :type filterid: int
    :param filterid: The ID (ref) of the filter
    :type date: datetime.date or str
    :param date: The date of the CCF
    :type time: datetime.time or str
    :param time: The time of the CCF
    :type duration: float
    :param duration: The total duration of the exported CCF
    :type components: str
    :param components: The name of the components used (ZZ, ZR, ...)
    :type sampling_rate: float
    :param sampling_rate: The sampling rate of the exported CCF
    :type day: bool
    :param day: Whether this function is called to export a daily stack (True)
        or each CCF (when keep_all parameter is set to True in the
        configuration). Defaults to True.
    :type ncorr: int
    :param ncorr: Number of CCF that have been stacked for this CCF.
    """

    output_folder = get_config(session, 'output_folder')
    export_format = get_config(session, 'export_format')
    sac, mseed = False, False
    if export_format == "BOTH":
        mseed = True
        sac = True
    elif export_format == "SAC":
        sac = True
    elif export_format == "MSEED":
        mseed = True

    if day:
        path = os.path.join("STACKS", "%02i" % filterid, "001_DAYS", components,
                            "%s_%s" % (station1, station2), str(date))
        pair = "%s:%s" % (station1, station2)
        if mseed:
            export_mseed(session, path, pair, components, filterid, CF/ncorr,
                         ncorr)
        if sac:
            export_sac(session, path, pair, components, filterid, CF/ncorr,
                       ncorr)

    else:
        file = '%s.cc' % time
        path = os.path.join(output_folder, "%02i" % filterid, station1,
                            station2, components, date)
        if not os.path.isdir(path):
            os.makedirs(path)

        t = Trace()
        t.data = CF
        t.stats.sampling_rate = sampling_rate
        t.stats.starttime = -float(get_config(session, 'maxlag'))
        t.stats.components = components
        # if ncorr != 0:
            # t.stats.location = "%02i"%ncorr
        st = Stream(traces=[t, ])
        st.write(os.path.join(path, file), format='mseed')
        del t, st


def export_sac(db, filename, pair, components, filterid, corr, ncorr=0,
               sac_format=None, maxlag=None, cc_sampling_rate=None):
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
    mytrace.stats.sac = AttribDict()

    mytrace.stats.sac.b = -maxlag
    mytrace.stats.sac.depmin = np.min(corr)
    mytrace.stats.sac.depmax = np.max(corr)
    mytrace.stats.sac.depmen = np.mean(corr)
    mytrace.stats.sac.scale = 1
    mytrace.stats.sac.npts = len(corr)

    st = Stream(traces=[mytrace, ])
    st.write(filename, format='SAC')
    del st
    return


def export_mseed(db, filename, pair, components, filterid, corr, ncorr=0,
                 maxlag=None, cc_sampling_rate=None):
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


def get_results(session, station1, station2, filterid, components, dates,
                mov_stack = 1, format="stack"):
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
            try:
                st = read(daystack)
                if not np.any(np.isnan(st[0].data)) and not np.any(np.isinf(st[0].data)):
                    stack += st[0].data
                    i += 1
            except:
                pass
        if i > 0:
            return i, stack / i
        else:
            return 0, None

    elif format == "matrix":
        stack = np.zeros((len(dates), get_maxlag_samples(session))) * np.nan
        i = 0
        base = os.path.join("STACKS", "%02i"%filterid,
                                    "%03i_DAYS"%mov_stack, components,
                                    "%s_%s"%(station1, station2), "%s")
        if export_format == "BOTH":
            base += ".MSEED"
        elif export_format == "SAC":
            base += ".SAC"
        elif export_format == "MSEED":
            base += ".MSEED"
        for j, date in enumerate(dates):
            daystack = base % str(date)
            try:
                stack[j][:] = read(daystack)[0].data
                i += 1
            except:
                pass
        return i, stack


def get_results_all(session, station1, station2, filterid, components, dates,
                mov_stack = 1, format="stack"):

    output_folder = get_config(session, 'output_folder')
    path = os.path.join(output_folder, "%02i" % int(filterid),
                        station1, station2, components)
    results = []
    for date in dates:
        fname = os.path.join(path, date.strftime('%Y-%m-%d.h5'))
        if os.path.isfile(fname):
            df = pd.read_hdf(fname, 'data', parse_dates=True)
            df.index = pd.to_datetime(df.index)
            results.append(df)
    result = pd.concat(results)
    del results
    return result

# Some helper functions


def get_maxlag_samples(session):
    """
    Returns the length of the CC functions. Gets the maxlag and sampling rate
    from the database.


    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: int
    :returns: the length of the CCF
    """

    maxlag = float(get_config(session, 'maxlag'))
    cc_sampling_rate = float(get_config(session, 'cc_sampling_rate'))
    return int(2*maxlag*cc_sampling_rate)+1


def get_t_axis(session):
    """
    Returns the time axis (in seconds) of the CC functions.
    Gets the maxlag from the database and uses `get_maxlag_samples` function.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: :class:`numpy.array`
    :returns: the time axis
    """

    maxlag = float(get_config(session, 'maxlag'))
    samples = get_maxlag_samples(session)
    return np.linspace(-maxlag, maxlag, samples)


def get_components_to_compute(session):
    """
    Returns the components configured in the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list of str
    :returns: a list of components to compute
    """

    components_to_compute = []
    for comp in ['ZZ', 'RR', 'TT', 'TR', 'RT', 'ZR', 'RZ', 'TZ', 'ZT']:
        if get_config(session, comp, isbool=True):
            components_to_compute.append(comp)
    return components_to_compute


def build_ref_datelist(session):
    """
    Creates a date array for the REF.
    The returned tuple contains a start and an end date, and a list of
    individual dates between the two.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: tuple
    :returns: (start, end, datelist)
    """
    begin = get_config(session, "ref_begin")
    end = get_config(session, "ref_end")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d').date()
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()
    end = min(end, datetime.date.today())
    datelist = pd.date_range(start, end).map(lambda x: x.date())
    return start, end, datelist.tolist()


def build_movstack_datelist(session):
    """
    Creates a date array for the analyse period.
    The returned tuple contains a start and an end date, and a list of
    individual dates between the two.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: tuple
    :returns: (start, end, datelist)
    """
    begin = get_config(session, "startdate")
    end = get_config(session, "enddate")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d').date()
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()
    end = min(end, datetime.date.today())
    datelist = pd.date_range(start, end).map(lambda x: x.date())
    return start, end, datelist.tolist()


def updated_days_for_dates(session, date1, date2, pair, jobtype='CC',
                           interval=datetime.timedelta(days=1),
                           returndays=False):
    """
    Determines if any Job of jobtype=`jobtype` and for pair=`pair`,
    concerning a date between `date1` and `date2` has been modified in the last
    interval=`interval`.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type date1: datetime.datetime
    :param date1: Beginning of the period of interest
    :type date2: datetime.datetime
    :param date2: End of the period of interest
    :type pair: str
    :param pair: Pair of interest
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type interval: datetime.timedelta
    :param interval: Interval of time before now to search for updated days
    :type returndays: bool
    :param returndays: Whether to return a list of days (True) or not (False,
        default)

    :rtype: list or bool
    :returns: List of days if returndays is True, only "True" if not.
        (not clear!)
    """
    lastmod = datetime.datetime.now() - interval
    if pair == '%':
        days = session.query(Job).filter(Job.day >= date1).filter(Job.day <= date2).filter(Job.jobtype == jobtype).filter(Job.lastmod >= lastmod).group_by(Job.day).order_by(Job.day).all()
    else:
        days = session.query(Job).filter(Job.pair == pair).filter(Job.day >= date1).filter(Job.day <= date2).filter(Job.jobtype == jobtype).filter(Job.lastmod >= lastmod).group_by(Job.day).order_by(Job.day).all()
    logging.debug('Found %03i updated days' % len(days))
    if returndays and len(days) != 0:
        return [datetime.datetime.strptime(day.day,'%Y-%m-%d').date() for day in days] ## RETURN DATE LIST !!!
    elif returndays and len(days) == 0:
        return []
    else:
        return True

# MISC


def azimuth(coordinates, x0, y0, x1, y1):
    """
    Returns the azimuth between two coordinate sets.

    :type coordinates: str
    :param coordinates: {'DEG', 'UTM', 'MIX'}
    :type x0: float
    :param x0: X coordinate of station 1
    :type y0: float
    :param y0: Y coordinate of station 1
    :type x1: float
    :param x1: X coordinate of station 2
    :type y1: float
    :param y1: Y coordinate of station 2

    :rtype: float
    :returns: The azimuth in degrees
    """
    if coordinates == "DEG":
        dist, azim, bazim = gps2DistAzimuth(y0, x0, y1, x1)
        return azim
    elif coordinates == 'UTM':
        azim = 90. - np.arctan2((y1 - y0), (x1 - x0)) * 180. / np.pi
        return azim
    else:
        print "Please consider having a single coordinate system for\
            all stations"
        return 0


def nextpow2(x):
    """
    Returns the next power of 2 of `x`.

    :type x: int
    :param x: any value

    :rtype: int
    :returns: the next power of 2 of `x`
    """

    return np.ceil(np.log2(np.abs(x)))


def check_and_phase_shift(trace):
    # print trace
    taper_length = 20.0
    if trace.stats.npts < 4 * taper_length*trace.stats.sampling_rate:
        trace.data = np.zeros(trace.stats.npts)
        return trace

    dt = np.mod(trace.stats.starttime.datetime.microsecond*1.0e-6,
                trace.stats.delta)
    if (trace.stats.delta - dt) <= np.finfo(float).eps:
        dt = 0
    if dt != 0:
        if dt <= (trace.stats.delta / 2.):
            dt = -dt
#            direction = "left"
        else:
            dt = (trace.stats.delta - dt)
#            direction = "right"
        trace.detrend(type="demean")
        trace.detrend(type="simple")
        taper_1s = taper_length * float(trace.stats.sampling_rate) / trace.stats.npts
        cp = cosTaper(trace.stats.npts, taper_1s)
        trace.data *= cp

        n = int(2**nextpow2(len(trace.data)))
        FFTdata = scipy.fftpack.fft(trace.data, n=n)
        fftfreq = scipy.fftpack.fftfreq(n, d=trace.stats.delta)
        FFTdata = FFTdata * np.exp(1j * 2. * np.pi * fftfreq * dt)
        trace.data = np.real(scipy.fftpack.ifft(FFTdata, n=n)[:len(trace.data)])
        trace.stats.starttime += dt
        return trace
    else:
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
        gap_list.append([_i, _i+1,
                        stats['network'], stats['station'],
                        stats['location'], stats['channel'],
                        stime, etime, delta, nsamples])
    # Set the original traces to not alter the stream object.
    stream.traces = copied_traces
    return gap_list


if __name__ == "__main__":
    s = connect()
    for filter in get_filters(s, False):
        print filter.ref, filter.low, filter.high

    print get_networks(s)

    for station in get_stations(s, False, net='BE'):
        print station.net, station.sta

    print get_config(s)
    print get_config(s, 'data_folder')
