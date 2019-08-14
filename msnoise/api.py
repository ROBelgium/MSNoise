import collections
import copy
import datetime
import itertools
import logging
import os
import glob
import traceback
try:
    import cPickle
except:
    import pickle as cPickle
import math
import pkg_resources
import sys

from logbook import Logger, StreamHandler
import sys

from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool
from sqlalchemy.sql.expression import func
import numpy as np
import pandas as pd
import scipy as sp
import scipy.fftpack
from scipy.fftpack.helper import next_fast_len
import scipy.fftpack._fftpack as sff
import scipy.optimize

from obspy.core import Stream, Trace, read, AttribDict, UTCDateTime
from obspy import read_inventory
from obspy.io.xseed import Parser
from obspy.geodetics import gps2dist_azimuth

from . import DBConfigNotFoundError
from .msnoise_table_def import Filter, Job, Station, Config, DataAvailability


def get_logger(name, loglevel=None, with_pid=False):
    """
    Returns the current configured logger or configure a new one.
    """
    # if with_pid:
    #     log_fmt='%(asctime)s msnoise [pid %(process)d] '\
    #             '[%(levelname)s] %(message)s'
    # else:
    #     log_fmt='%(asctime)s msnoise [%(levelname)s] %(message)s'
    # logger = logging.getLogger(name)
    # # Remove any inherited StreamHandler to avoid duplicate lines
    # for h in logger.handlers:
    #     if isinstance(h, logging.StreamHandler):
    #         logger.removeHandler(h)
    # handler = logging.StreamHandler(sys.stderr)
    # handler.setFormatter(
    #         logging.Formatter(fmt=log_fmt, datefmt='%Y-%m-%d %H:%M:%S'))
    # logger.addHandler(handler)
    # logger.setLevel(loglevel)
    # logger.propagate = False

    if with_pid:
        log_fmt="{record.time} msnoise [pid {record.process}]" \
                "[{record.level_name}]: {record.message}"
    else:
        log_fmt="{record.time} msnoise [{record.level_name}]: {record.message}"

    StreamHandler(sys.stdout, format_string=log_fmt, 
                  level=loglevel).push_application()
    logger = Logger(name)

    return logger


def get_engine(inifile=None):
    """Returns the a SQLAlchemy Engine

    :type inifile: str
    :param inifile: The path to the db.ini file to use. Defaults to os.cwd() +
        db.ini

    :rtype: :class:`sqlalchemy.engine.Engine`
    :returns: An :class:`~sqlalchemy.engine.Engine` Object
    """
    dbini = read_db_inifile(inifile)
    if dbini.tech == 1:
        engine = create_engine('sqlite:///%s' % dbini.hostname, echo=False,
                               connect_args={'check_same_thread': False})
    else:
        engine = create_engine('mysql+pymysql://%s:%s@%s/%s'
                               % (dbini.username, dbini.password,
                                  dbini.hostname, dbini.database),
                               echo=False, poolclass=NullPool,
                               connect_args={'connect_timeout': 15})
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


def create_database_inifile(tech, hostname, database, username, password,
                            prefix=""):
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
    :param prefix: The prefix to use for all tables
    :type prefix: string
    :param password: The password of `user`

    :return: None
    """
    f = open(os.path.join(os.getcwd(), 'db.ini'), 'wb')
    cPickle.dump([tech, hostname, database, username, password, prefix], f,
                 protocol=2)
    f.close()


def read_db_inifile(inifile=None):
    """Reads the parameters from the db.ini file.

    :type inifile: string
    :param inifile: The path to the db.ini file to use. Defaults to os.cwd() +
        db.ini

    :rtype: tuple
    :returns: tech, hostname, database, username, password
    """
    IniFile = collections.namedtuple('IniFile', ['tech', 'hostname',
        'database', 'username', 'password', 'prefix'])
    if not inifile:
        inifile = os.path.join(os.getcwd(), 'db.ini')

    try:
        f = open(inifile, 'rb')
    #except FileNotFoundError:  # This is better but only for python3
    except IOError:
        raise DBConfigNotFoundError(
                "No db.ini file in this directory, please run "
                "'msnoise db init' in this folder to initialize it as "
                "an MSNoise project folder.")
    try:
        # New ini file with prefix support
        tech, hostname, database, username, password, prefix = cPickle.load(f)
    except:
        # Old ini file without prefix
        tech, hostname, database, username, password = cPickle.load(f)
        prefix = ""
    f.close()
    return IniFile(tech, hostname, database, username, password, prefix)


# CONFIG


def get_config(session, name=None, isbool=False, plugin=None):
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
    :type plugin: str
    :param plugin: if provided, gives the name of the Plugin config to use. E.g.
        if "Amazing" is provided, MSNoise will try to load the "AmazingConfig"
        entry point. See :doc:`plugins` for details.

    :rtype: str, bool or dict
    :returns: the value for `name` or a dict of all config values
    """
    if plugin:
        for ep in pkg_resources.iter_entry_points(
                group='msnoise.plugins.table_def'):
            if ep.name.replace("Config", "") == plugin:
                table = ep.load()
    else:
        table = Config
    if name:
        config = session.query(table).filter(table.name == name).first()
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


def update_config(session, name, value, plugin=None):
    """Update one config bit in the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type name: str
    :param name: The name of the config bit to set.

    :type value: str
    :param value: The value of parameter `name`. Can also be NULL if you don't
        want to use this particular parameter.

    :type plugin: str
    :param plugin: if provided, gives the name of the Plugin config to use. E.g.
        if "Amazing" is provided, MSNoise will try to load the "AmazingConfig"
        entry point. See :doc:`plugins` for details.

    """
    if plugin:
        for ep in pkg_resources.iter_entry_points(
                group='msnoise.plugins.table_def'):
            if ep.name.replace("Config", "") == plugin:
                table = ep.load()
    else:
        table = Config
    config = session.query(table).filter(table.name == name).first()
    if "NULL" in value:
        config.value = None
    else:
        config.value = value
    session.commit()
    return


def get_params(session):
    """Get config parameters from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :returns: a Param class containing the parameters
    """
    # TODO: this could be populated automatically from defauts iff defaults
    # would mention types
    from obspy.core.util.attribdict import AttribDict
    from .default import default
    s = session
    params = AttribDict()
    for name in default.keys():
        if len(default[name]) == 2:
            params[name] = get_config(s, name)
        else:
            txt, _, itemtype = default[name]
            if itemtype is bool:
                params[name] = get_config(s, name, isbool=True)
            else:
                params[name] = itemtype(get_config(s, name))

    # TODO remove reference to goal_sampling_rate
    params.goal_sampling_rate = params.cc_sampling_rate
    params.min30 = params.corr_duration * params.goal_sampling_rate
    params.components_to_compute = get_components_to_compute(s)
    params.components_to_compute_single_station = get_components_to_compute_single_station(s)
    params.all_components = np.unique(params.components_to_compute_single_station + \
                            params.components_to_compute)
    
    return params

# FILTERS PART


def get_filters(session, all=False):
    """Get Filters from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type all: bool
    :param all: Returns all filters from the database if True, or only filters
        where `used` = 1 if False (default)

    :rtype: list of :class:`~msnoise.msnoise_table_def.declare_tables.Filter`
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

    .. seealso:: :class:`msnoise.msnoise_table_def.declare_tables.Filter`

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
        filter = Filter()
        filter.low = low
        filter.high = high
        filter.mwcs_low = mwcs_low
        filter.mwcs_high = mwcs_high
        filter.rms_threshold = rms_threshold
        filter.mwcs_wlen = mwcs_wlen
        filter.mwcs_step = mwcs_step
        filter.used = used
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
        networks = session.query(Station).filter(Station.used == True).\
            group_by(Station.net)
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

    :rtype: list of :class:`msnoise.msnoise_table_def.declare_tables.Station`
    :returns: list of :class:`~msnoise.msnoise_table_def.declare_tables.Station`
    """
    q = session.query(Station)
    if all:
        if net is not None:
            stations = q.filter(Station.net == net).order_by(Station.net).\
                order_by(Station.sta)
        else:
            stations = q.order_by(Station.net).order_by(Station.sta).all()
    else:
        stations = q.filter(Station.used == True).order_by(Station.net).\
            order_by(Station.sta)
        if net is not None:
            stations = stations.filter(Station.net == net).\
                order_by(Station.net).order_by(Station.sta)
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

    :rtype: :class:`msnoise.msnoise_table_def.declare_tables.Station`
    :returns: a :class:`~msnoise.msnoise_table_def.declare_tables.Station` Object

    """
    station = session.query(Station).filter(Station.net == net).\
        filter(Station.sta == sta).first()
    return station


def update_station(session, net, sta, X, Y, altitude, coordinates='UTM',
                   instrument='N/A', used=1):
    """Updates or Insert a new Station in the database.

    .. seealso :: :class:`msnoise.msnoise_table_def.declare_tables.Station`
    
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
    station = session.query(Station).filter(Station.net == net).\
        filter(Station.sta == sta).first()
    if station is None:
        station = Station(net, sta, X, Y, altitude, coordinates, instrument,
                          used)
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
    :returns: An iterable of :class:`~msnoise.msnoise_table_def.declare_tables.Station` object
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

    :type station1: :class:`~msnoise.msnoise_table_def.declare_tables.Station`
    :param station1: A Station object
    :type station2: :class:`~msnoise.msnoise_table_def.declare_tables.Station`
    :param station2: A Station object
    :type coordinates: str
    :param coordinates: The coordinates system. "DEG" is WGS84 latitude/
        longitude in degrees. "UTM" is expressed in meters.

    :rtype: float
    :returns: The interstation distance in km
    """

    if coordinates == "DEG":
        dist, azim, bazim = gps2dist_azimuth(station1.Y, station1.X,
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

    data = session.query(DataAvailability).\
        with_hint(DataAvailability, 'USE INDEX (da_index)'). \
        filter(DataAvailability.path == path). \
        filter(DataAvailability.file == file).\
        filter(DataAvailability.net == net).\
        filter(DataAvailability.sta == sta).\
        filter(DataAvailability.comp == comp).first()
    if data is None:
        flag = "N"
        data = DataAvailability(net, sta, comp, path, file, starttime, endtime,
                                data_duration, gaps_duration, samplerate, flag)
        session.add(data)
        toreturn = 1
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
            toreturn = -1
        else:
            toreturn = 0
    session.commit()
    return toreturn


def get_new_files(session):
    """
    Returns the files marked "N"ew or "M"odified in the database

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list
    :returns: list of :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability`
    """

    files = session.query(DataAvailability).\
        filter(DataAvailability.flag != 'A').\
        order_by(DataAvailability.starttime).all()
    return files


def get_data_availability(session, net=None, sta=None, comp=None,
                          starttime=None, endtime=None):
    """
    Returns the :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability` objects
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
    :returns: list of :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability`
    """

    if not starttime:
        data = session.query(DataAvailability).\
            filter(DataAvailability.net == net).\
            filter(DataAvailability.sta == sta).\
            filter(DataAvailability.comp == comp).all()
    elif not net:
        data = session.query(DataAvailability).\
            filter(DataAvailability.starttime <= endtime).\
            filter(DataAvailability.endtime >= starttime).all()
    else:
        data = session.query(DataAvailability).\
            filter(DataAvailability.net == net).\
            filter(DataAvailability.sta == sta).\
            filter(func.DATE(DataAvailability.starttime) <= endtime.date()).\
            filter(func.DATE(DataAvailability.endtime) >= starttime.date()).all()
        if not len(data):
            data = session.query(DataAvailability). \
                filter(DataAvailability.sta == "MULTIPLEX"). \
                filter(func.DATE(DataAvailability.starttime) <= endtime.date()).\
                filter(
                func.DATE(DataAvailability.endtime) >= starttime.date()).all()
    return data


def mark_data_availability(session, net, sta, flag):
    """
    Updates the flag of all
    :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability` objects matching
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
    logging.debug("Updating: %s %s to flag=%s" %(net, sta, flag))
    da = DataAvailability.__table__
    stmt = da.update().where(da.c.sta==sta).where(da.c.net==net).\
        values(flag=flag)
    session.execute(stmt)
    session.commit()


def count_data_availability_flags(session):
    """
    Count the number of :class:`~msnoise.msnoise_table_def.declare_tables.DataAvailability`,
    grouped by `flag`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list
    :returns: list of [count, flag] pairs
    """

    return session.query(func.count(DataAvailability.flag),
                         DataAvailability.flag).\
        group_by(DataAvailability.flag).all()


# Jobs
# TODO bad doing this here!
import time
def update_job(session, day, pair, jobtype, flag, commit=True, returnjob=True,
               ref=None):
    """
    Updates or Inserts a new :class:`~msnoise.msnoise_table_def.declare_tables.Job` in the
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


    :rtype: :class:`~msnoise.msnoise_table_def.declare_tables.Job` or None
    :returns: If returnjob is True, returns the modified/inserted Job.
    """
    if ref:
        job = session.query(Job).filter(Job.ref == ref).first()
    else:
        job = session.query(Job)\
            .with_hint(Job, 'USE INDEX (job_index)')\
            .filter(Job.day == day)\
            .filter(Job.pair == pair)\
            .filter(Job.jobtype == jobtype).first()
    if job is None:
        job = Job(day, pair, jobtype, 'T')
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
    :class:`~msnoise.msnoise_table_def.declare_tables.Job`. This method uses the Engine
    directly, no need to pass a Session object.

    :type jobs: list
    :param jobs: a list of :class:`~msnoise.msnoise_table_def.declare_tables.Job` to insert.
    """
    engine = get_engine()
    engine.execute(
        Job.__table__.insert(),
        jobs)


def massive_update_job(session, jobs, flag="D"):
    """
    Routine to use a low level function to update much faster a list of
    :class:`~msnoise.msnoise_table_def.declare_tables.Job`. This method uses the Job.ref
    which is unique.

    :type jobs: list
    :param jobs: a list of :class:`~msnoise.msnoise_table_def.declare_tables.Job` to update.
    :type flag: str
    :param flag: The destination flag.
    """
    updated = False
    mappings = [{'ref': job.ref, 'flag': flag} for job in jobs]
    while not updated:
        try:
            session.bulk_update_mappings(Job, mappings)
            session.commit()
            updated = True
        except:
            time.sleep(np.random.random())
            pass
    return

def is_next_job(session, flag='T', jobtype='CC'):
    """
    Are there any :class:`~msnoise.msnoise_table_def.declare_tables.Job` in the database,
    with flag=`flag` and jobtype=`type`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.

    :rtype: bool
    :returns: True if at least one :class:`~msnoise.msnoise_table_def.declare_tables.Job`
        matches, False otherwise.
    """
    job = session.query(Job).with_hint(Job, 'USE INDEX (job_index2)').\
        filter(Job.jobtype == jobtype).\
        filter(Job.flag == flag).first()
    if job is None:
        return False
    else:
        return True


def get_next_job(session, flag='T', jobtype='CC'):
    """
    Get the next :class:`~msnoise.msnoise_table_def.declare_tables.Job` in the database,
    with flag=`flag` and jobtype=`jobtype`. Jobs of the same `type` are grouped
    per day. This function also sets the flag of all selected Jobs to "I"n
    progress.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.

    :rtype: list
    :returns: list of :class:`~msnoise.msnoise_table_def.declare_tables.Job`
    """
    tmp = []
    while not len(tmp):
        jobs = session.query(Job).filter(Job.jobtype == jobtype).\
            filter(Job.flag == flag).\
            filter(Job.day == session.query(Job).
                   with_hint(Job, 'USE INDEX (job_index)').
                   filter(Job.jobtype == jobtype).
                   filter(Job.flag == flag).first().day).\
            with_for_update()
        # print(jobs.statement.compile(compile_kwargs={"literal_binds": True}))
        tmp = jobs.all()
        jobs.update({Job.flag: 'I'})
        session.commit()
    return tmp


def is_dtt_next_job(session, flag='T', jobtype='DTT', ref=False):
    """
    Are there any DTT :class:`~msnoise.msnoise_table_def.declare_tables.Job` in the database,
    with flag=`flag` and jobtype=`jobtype`. If `ref` is provided, checks if a
    DTT "REF" job is present.

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
    q = session.query(Job.ref).with_hint(Job, 'USE INDEX (job_index2)').\
        filter(Job.flag == flag).\
        filter(Job.jobtype == jobtype)
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
    Get the next DTT :class:`~msnoise.msnoise_table_def.declare_tables.Job` in the database,
    with flag=`flag` and jobtype=`jobtype`. Jobs are then grouped per station
    pair. This function also sets the flag of all selected Jobs to "I"n
    progress.

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
    if read_db_inifile().tech == 1:
        rand = func.random
    else:
        rand = func.rand
    try:
        jobs = session.query(Job.ref, Job.day, Job.pair, Job.flag).filter(Job.flag == flag).\
            filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
            filter(Job.pair == session.query(Job).
                   with_hint(Job, 'USE INDEX (job_index)').
                   filter(Job.flag == flag).
                   filter(Job.jobtype == jobtype).
                   filter(Job.day != 'REF').
                   order_by(rand()).first().pair).\
            with_for_update()
        tmp = jobs.all()
        mappings = [{'ref': job.ref, 'flag': "I"} for job in tmp]
        updated = False
        while not updated:
            try:
                session.bulk_update_mappings(Job, mappings)
                session.commit()
                updated=True
            except:
                traceback.print_exc()
                time.sleep(np.random.random())
                pass
        return tmp
    except:
        traceback.print_exc()
        return []


def reset_jobs(session, jobtype, alljobs=False, rule=None):
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
    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    jobs = session.query(Job).filter(Job.jobtype == jobtype)
    if rule:
        session.execute("UPDATE %sjobs set flag='T' where jobtype='%s' and  %s"
                        % (prefix, jobtype, rule))
        session.commit()
        return
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

    jobs = session.query(Job).filter(Job.pair == pair).\
        filter(Job.jobtype == "DTT")
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

    return session.query(func.count(Job.flag), Job.flag).\
        filter(Job.jobtype == jobtype).group_by(Job.flag).all()


def get_jobs_by_lastmod(session, jobtype='CC', lastmod=datetime.datetime.now()):
    """

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type jobtype: str
    :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
    :type lastmod: datetime.datetime
    :param lastmod: Jobs' modification time

    :rtype: list
    :return: list of Job objects.
    """
    jobs = session.query(Job).filter(Job.jobtype == jobtype).\
        filter(Job.lastmod >= lastmod).all()
    return jobs


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


def export_allcorr2(session, ccfid, data):
    output_folder = get_config(session, 'output_folder')
    station1, station2, components, filterid, date = ccfid.split('_')

    path = os.path.join(output_folder, "%02i" % int(filterid),
                        station1, station2, components)
    if not os.path.isdir(path):
        os.makedirs(path)

    df = pd.DataFrame().from_dict(data).T
    df.columns = get_t_axis(session)
    df.to_hdf(os.path.join(path, date+'.h5'), 'data')
    del df
    return


def add_corr(session, station1, station2, filterid, date, time, duration,
             components, CF, sampling_rate, day=False, ncorr=0, params=None):
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
    :type params: dict
    :param params: A dictionnary of MSNoise config parameters as returned by
        :func:`get_params`.
    """

    output_folder = params.output_folder
    export_format = params.export_format
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
            export_mseed(session, path, pair, components, filterid, CF,
                         ncorr, params=params)
        if sac:
            export_sac(session, path, pair, components, filterid, CF,
                       ncorr, params=params)

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
               sac_format=None, maxlag=None, cc_sampling_rate=None,
               params=None):
    maxlag = params.maxlag
    cc_sampling_rate = params.goal_sampling_rate
    sac_format = params.sac_format
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
    if maxlag:
        mytrace.stats.starttime = -maxlag
    mytrace.stats.sac = AttribDict()
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
                 maxlag=None, cc_sampling_rate=None, params=None):
    try:
        os.makedirs(os.path.split(filename)[0])
    except:
        pass
    filename += ".MSEED"
    maxlag = params.maxlag
    cc_sampling_rate = params.goal_sampling_rate
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


def stack(data, stack_method="linear", pws_timegate=10.0, pws_power=2,
          goal_sampling_rate=20.0):
    """
    :type data: :class:`numpy.ndarray`
    :param data: the data to stack, each row being one CCF
    :type stack_method: str
    :param stack_method: either ``linear``: average of all CCF or ``pws`` to
        compute the phase weigthed stack. If ``pws`` is selected,
        the function expects the ``pws_timegate`` and ``pws_power``.
    :type pws_timegate: float
    :param pws_timegate: PWS time gate in seconds. Width of the smoothing
         window to convolve with the PWS spectrum.
    :type pws_power: float
    :param pws_power: Power of the PWS weights to be applied to the CCF stack.
    :type goal_sampling_rate: float
    :param goal_sampling_rate: Sampling rate of the CCF array submitted
    :rtype: :class:`numpy.array`
    :return: the stacked CCF.
    """
    if len(data) == 0:
        logging.debug("No data to stack.")
        return []
    data = data[~np.isnan(data).any(axis=1)]
    sanitize = False
    # TODO clean sanitize function, add param to config and make sure not to
    # kill the data[i] if all data are corrcoeff >0.9 (either very stable
    # corr or autocorr, then this sanitize should not occur.
    if len(data) != 1 and sanitize:
        threshold = 0.99
        npts = data.shape[1]
        corr = data.mean(axis=0)
        corrcoefs = np.array([np.corrcoef(di, corr)[1][0] for di in data])
        toolarge = np.where(corrcoefs >= threshold)[0]
        if len(toolarge):
            data = data[np.where(corrcoefs <= threshold)[0]]

    if len(data) == 0:
        return []
    if stack_method == "linear":
        # logging.debug("Doing a linear stack")
        corr = data.mean(axis=0)

    elif stack_method == "pws":
        # logging.debug("Doing a PWS stack")
        corr = np.zeros(data.shape[1], dtype='f8')
        phasestack = np.zeros(data.shape[1], dtype='c8')
        for i in range(data.shape[0]):
            data[i] -= data[i].mean()
        for c in data:
            phase = np.angle(sp.signal.hilbert(c))
            phasestack.real += np.cos(phase)
            phasestack.imag += np.sin(phase)
        coh = 1. / data.shape[0] * np.abs(phasestack)

        timegate_samples = int(pws_timegate * goal_sampling_rate)
        coh = np.convolve(sp.signal.boxcar(timegate_samples) /
                          timegate_samples, coh, 'same')
        coh = np.power(coh, pws_power)
        for c in data:
            corr += c * coh
        corr /= data.shape[0]

    return corr


def get_results(session, station1, station2, filterid, components, dates,
                mov_stack = 1, format="stack", params=None):
    """

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET_STA)
    :type station2: str
    :param station2: The name of station 2 (formatted NET_STA)
    :type filterid: int
    :param filterid: The ID (ref) of the filter
    :type components: str
    :param components: The name of the components used (ZZ, ZR, ...)
    :type dates: list
    :param dates: List of TODO datetime.datetime
    :type mov_stack: int
    :param mov_stack: Moving window stack.
    :type format: str
    :param format: Either ``stack``: the data will be stacked according to
        the parameters passed with ``params`` or ``matrix``: to get a 2D
        array of CCF.
    :type params: dict
    :param params: A dictionnary of MSNoise config parameters as returned by
        :func:`get_params`.
    :rtype: :class:`numpy.ndarray`
    :return: Either a 1D CCF (if format is ``stack`` or a 2D array (if format=
        ``matrix``).
    """
    if not params:
        export_format = get_config(session, 'export_format')
    else:
        export_format = params.export_format

    stack_data = np.zeros((len(dates), get_maxlag_samples(session))) * np.nan
    i = 0
    base = os.path.join("STACKS", "%02i" % filterid,
                        "%03i_DAYS" % mov_stack, components,
                        "%s_%s" % (station1, station2), "%s")
    if export_format == "BOTH":
        base += ".MSEED"
        export_format = "MSEED"
    elif export_format == "SAC":
        base += ".SAC"
    elif export_format == "MSEED":
        base += ".MSEED"
    logging.debug("Reading files...")
    for j, date in enumerate(dates):
        daystack = base % str(date)
        try:
            stack_data[j, :] = read(daystack, format=export_format)[0].data[:]
            i += 1
        except:
            pass

    if format == "matrix":
        return i, stack_data

    elif format == "stack":
        logging.debug("Stacking...")

        corr = stack(stack_data, params.stack_method, params.pws_timegate,
                     params.pws_power, params.goal_sampling_rate)

        if i > 0:
            return i, corr
        else:
            return 0, None


def get_results_all(session, station1, station2, filterid, components, dates):
    """
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET_STA)
    :type station2: str
    :param station2: The name of station 2 (formatted NET_STA)
    :type filterid: int
    :param filterid: The ID (ref) of the filter
    :type components: str
    :param components: The name of the components used (ZZ, ZR, ...)
    :type dates: list
    :param dates: List of TODO datetime.datetime
    :rtype: :class:`pandas.DataFrame`
    :return: All CCF results in a :class:`pandas.DataFrame`, where the index
        is the time of the CCF and the columns are the times in the coda.
    """

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
    :returns: the length of the CCF in samples
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
    :returns: the time axis in seconds
    """

    maxlag = float(get_config(session, 'maxlag'))
    samples = get_maxlag_samples(session)
    return np.linspace(-maxlag, maxlag, samples)


def get_components_to_compute(session, plugin=None):
    """
    Returns the components configured in the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list of str
    :returns: a list of components to compute
    """

    components_to_compute = get_config(session, "components_to_compute",
                                       plugin=plugin)
    if len(components_to_compute) == 0:
        return []
    elif components_to_compute.count(",") == 0:
        components_to_compute = [components_to_compute,]
    else:
        components_to_compute = components_to_compute.split(",")
    return components_to_compute


def get_components_to_compute_single_station(session, plugin=None):
    """
    Returns the components configured in the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: list of str
    :returns: a list of components to compute
    """

    components_to_compute = get_config(session, "components_to_compute_single_station",
                                       plugin=plugin)
    if len(components_to_compute) == 0:
        return []
    elif components_to_compute.count(",") == 0:
        components_to_compute = [components_to_compute,]
    else:
        components_to_compute = components_to_compute.split(",")
    return list(np.unique([''.join(sorted(a)) for a in components_to_compute]))


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
    elif begin == "1970-01-01":
        start = session.query(DataAvailability).order_by(
            DataAvailability.starttime).first().starttime.date()
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()
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
    elif begin == "1970-01-01": # TODO this fails when the DA is empty
        start = session.query(DataAvailability).order_by(
            DataAvailability.starttime).first().starttime.date()
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()
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
        days = session.query(Job).filter(Job.day >= date1).\
            filter(Job.day <= date2).filter(Job.jobtype == jobtype).\
            filter(Job.lastmod >= lastmod).group_by(Job.day).\
            order_by(Job.day).all()
    else:
        days = session.query(Job).filter(Job.pair == pair).\
            filter(Job.day >= date1).filter(Job.day <= date2).\
            filter(Job.jobtype == jobtype).filter(Job.lastmod >= lastmod).\
            group_by(Job.day).order_by(Job.day).all()
    logging.debug('Found %03i updated days' % len(days))
    if returndays and len(days) != 0:
        return [datetime.datetime.strptime(day.day,'%Y-%m-%d').date() for day
                in days] ## RETURN DATE LIST !!!
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
        dist, azim, bazim = gps2dist_azimuth(y0, x0, y1, x1)
        return azim
    elif coordinates == 'UTM':
        if (np.isclose(y0, y1) & np.isclose(x0, x1)):
            return 0
        else:
            azim = 90. - np.arctan2((y1 - y0), (x1 - x0)) * 180. / np.pi
            return azim % 360
    else:
        logging.warning("Please consider having a single coordinate system for"
                        " all stations")
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
    # TODO replace this hard coded taper length
    taper_length = 20.0
    if trace.stats.npts < 4 * taper_length*trace.stats.sampling_rate:
        trace.data = np.zeros(trace.stats.npts)
        return trace

    dt = np.mod(trace.stats.starttime.datetime.microsecond*1.0e-6,
                trace.stats.delta)
    if (trace.stats.delta - dt) <= np.finfo(float).eps:
        dt = 0.
    if dt != 0.:
        if dt <= (trace.stats.delta / 2.):
            dt = -dt
#            direction = "left"
        else:
            dt = (trace.stats.delta - dt)
#            direction = "right"
        logging.debug("correcting time by %.6fs"%dt)
        trace.detrend(type="demean")
        trace.detrend(type="simple")
        trace.taper(max_percentage=None, max_length=1.0)

        n = next_fast_len(int(trace.stats.npts))
        FFTdata = scipy.fftpack.fft(trace.data, n=n)
        fftfreq = scipy.fftpack.fftfreq(n, d=trace.stats.delta)
        FFTdata = FFTdata * np.exp(1j * 2. * np.pi * fftfreq * dt)
        FFTdata = FFTdata.astype(np.complex64)
        scipy.fftpack.ifft(FFTdata, n=n, overwrite_x=True)
        trace.data = np.real(FFTdata[:len(trace.data)]).astype(np.float)
        trace.stats.starttime += dt
        del FFTdata, fftfreq
        clean_scipy_cache()
        return trace
    else:
        return trace


def getGaps(stream, min_gap=None, max_gap=None):
    # Create shallow copy of the traces to be able to sort them later on.
    copied_traces = copy.copy(stream.traces)
    stream.sort()
    gap_list = []
    for _i in range(len(stream.traces) - 1):
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
    del copied_traces
    return gap_list


def make_same_length(st):
    """
    This function takes a stream of equal sampling rate and makes sure that all
    channels have the same length and the same gaps.
    """

    # Merge traces
    st.merge()

    # Initialize arrays to be filled with start+endtimes of all traces
    starttimes = []
    endtimes = []

    # Loop over all traces of the stream
    for tr in st:
        # Force conversion to masked arrays
        if not np.ma.count_masked(tr.data):
            tr.data = np.ma.array(tr.data, mask=False)
        # Read out start+endtimes of traces to trim
        starttimes.append(tr.stats.starttime)
        endtimes.append(tr.stats.endtime)

    # trim stream to common starttimes
    if max(starttimes) >= min(endtimes):
        return Stream()
    st.trim(max(starttimes), min(endtimes))

    # get the mask of all traces, i.e. the parts where at least one trace has
    # a gap
    # TODO add cases with more than 2 or 3 traces (could append?)
    # TODO is there a better way to AND masks ?
    if len(st) < 2:
        return st
    elif len(st) == 2:
        mask = np.logical_or(st[0].data.mask, st[1].data.mask)
    elif len(st) == 3:
        mask = np.logical_or(st[0].data.mask, st[1].data.mask, st[2].data.mask)

    # apply the mask to all traces
    for tr in st:
        tr.data.mask = mask

    st = st.split()
    return st


def clean_scipy_cache():
    """This functions wraps all destroy scipy cache at once. It is a workaround
    to the memory leak induced by the "caching" functions in scipy fft."""
    sff.destroy_zfft_cache()
    sff.destroy_zfftnd_cache()
    sff.destroy_drfft_cache()
    sff.destroy_cfft_cache()
    sff.destroy_cfftnd_cache()
    sff.destroy_rfft_cache()
    sff.destroy_ddct2_cache()
    sff.destroy_ddct1_cache()
    # sff.destroy_ddct4_cache()
    sff.destroy_dct2_cache()
    sff.destroy_dct1_cache()
    # sff.destroy_dct4_cache()
    sff.destroy_ddst2_cache()
    sff.destroy_ddst1_cache()
    sff.destroy_dst2_cache()
    sff.destroy_dst1_cache()
    scipy.fftpack.convolve.destroy_convolve_cache()


def preload_instrument_responses(session):
    """
    This function preloads all instrument responses from ``response_format``
    and stores the seed ids, start and end dates, and paz for every channel
    in a DataFrame.
    
    .. warning::
        This function only works for ``response_format`` being "inventory"
        or "dataless".
    
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect` 
    
    :rtype: pandas.DataFrame
    :returns: A table containing all channels with the time of operation and
        poles and zeros.
    
    """
    logging.debug('Preloading instrument response')
    response_format = get_config(session, 'response_format')
    files = glob.glob(os.path.join(get_config(session, 'response_path'), "*"))
    channels = []
    if response_format == "inventory":
        for file in files:
            logging.debug("Processing %s" % file)
            try:
                inv = read_inventory(file, format='STATIONXML')
                for net in inv.networks:
                    for sta in net.stations:
                        for cha in sta.channels:
                            seed_id = "%s.%s.%s.%s" % (net.code, sta.code,
                                                       cha.location_code,
                                                       cha.code)
                            resp = inv.get_response(seed_id, cha.start_date+10)
                            polezerostage = resp.get_paz()
                            totalsensitivity = resp.instrument_sensitivity
                            pzdict = {}
                            pzdict['poles'] = polezerostage.poles
                            pzdict['zeros'] = polezerostage.zeros
                            pzdict['gain'] = polezerostage.normalization_factor
                            pzdict['sensitivity'] = totalsensitivity.value
                            channels.append([seed_id, cha.start_date,
                                             cha.end_date or UTCDateTime(),
                                             pzdict, cha.latitude,
                                             cha.longitude])
            except:
                pass

    elif response_format == "dataless":
        for file in files:
            try:
                p = Parser(file)
                for channel in p.get_inventory()["channels"]:
                    resp = p.get_paz(channel["channel_id"],
                                     channel["start_date"]+10)
                    channels.append([channel["channel_id"],
                                     channel["start_date"],
                                     channel["end_date"] or UTCDateTime(),
                                     resp,
                                     channel["latitude"],
                                     channel["longitude"]])
            except:
                pass
    channels = pd.DataFrame(channels, columns=["channel_id", "start_date",
                                               "end_date", "paz", "latitude",
                                               "longitude"],)
    logging.debug('Finished Loading instrument responses')
    return channels
