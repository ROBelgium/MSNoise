"""
SQLAlchemy table definition.
"""

import datetime
import os
from collections import namedtuple
from sqlalchemy import Column, Integer, String, Float, Boolean, DateTime,\
    text, TIMESTAMP, Enum, REAL, UniqueConstraint, Index
from sqlalchemy.ext.declarative import declarative_base, declared_attr

try:
    import cPickle
except:
    import pickle as cPickle


# Name of the file in the current directory that keeps
# information pointing to the database to use.
DB_INI = 'db.ini'


def read_prefix(db_ini_filename=DB_INI):
    """
    Returns the table prefix set in the db_ini_filename file.
    """
    try:
        with open(db_ini_filename, 'rb') as ini_file:
            data = cPickle.load(ini_file)
            # Note: data should contain the following table
            # [tech hostname database username password prefix]
            # with prefix being optional for backward compatibility
            return data[5]
    except (IndexError, IOError):
        return ''


def declare_tables(prefix=None):
    """
    Define classes mapped to relational database tables and return the declared
    classes in a numedtuple.
    """
    # Implementation note: these definitions must be done in a function to
    # allow for the prefix to be evaluated at run time. If this was done at
    # module level, the code would be evaluated at compile time (when the
    # interpreter starts) and the prefix could not be specified at run time.

    # If no prefix is given, read it from the db.ini file
    if prefix is None:
        prefix = read_prefix()

    # Define the namedtuple to return
    sqlschema = namedtuple('SQLSchema', ['Base', 'PrefixerBase', 'Filter',
        'Job', 'Station', 'Config', 'DataAvailability'])

    # Create the SQLAlchemy base and subclass it to prefix the table names
    Base = declarative_base()

    class PrefixerBase(Base):
        """
        A Base abstract subclass that will add our prefix to table names.
        """
        __abstract__ = True

        @declared_attr
        def __tablename__(cls):
            table_prefix = prefix
            if prefix:
                table_prefix += '_'
            return table_prefix + cls.__incomplete_tablename__

    ########################################################################

    class Filter(PrefixerBase):
        """
        Filter base class.

        :type ref: int
        :param ref: The id of the Filter in the database
        :type low: float
        :param low: The lower frequency bound of the Whiten function (in Hz)
        :type high: float
        :param high: The upper frequency bound of the Whiten function (in Hz)
        :type mwcs_low: float
        :param mwcs_low: The lower frequency bound of the linear regression done in
            MWCS (in Hz)
        :type mwcs_high: float
        :param mwcs_high: The upper frequency bound of the linear regression done in
            MWCS (in Hz)
        :type mwcs_wlen: float
        :param mwcs_wlen: Window length (in seconds) to perform MWCS
        :type mwcs_step: float
        :param mwcs_step: Step (in seconds) of the windowing procedure in MWCS
        :type dtt_minlag: float
        :param dtt_minlag: If ``dtt_lag`` =static (in config table): min lag time (in seconds)
        :type dtt_width: float
        :param dtt_width: Width of the time lag window (in seconds)
        :type dtt_v: float
        :param dtt_v: If ``dttlag`` =dynamic (in config table): what velocity to use to avoid ballistic waves [1.0] km/s (default=1.0)
        :type used: bool
        :param used: Is the filter activated for the processing
        """

        __incomplete_tablename__ = "filters"

        ref = Column(Integer, primary_key=True)
        low = Column(Float())
        mwcs_low = Column(Float())
        high = Column(Float())
        mwcs_high = Column(Float())
        mwcs_wlen = Column(Float())
        mwcs_step = Column(Float())
        dtt_minlag = Column(Float())
        dtt_width = Column(Float())
        dtt_v = Column(Float())
        used = Column(Boolean(True))

        def __init__(self, **kwargs):
            """"""
            # self.low = low
            # self.mwcs_low = mwcs_low
            # self.high = high
            # self.mwcs_high = mwcs_high
            # self.mwcs_wlen = mwcs_wlen
            # self.mwcs_step = mwcs_step
            # self.dtt_minlag = dtt_minlag
            # self.dtt_width = dtt_width
            # self.dtt_v = dtt_v
            # self.used = used

    ########################################################################

    class Job(PrefixerBase):
        """
        Job Object

        :type ref: int
        :param ref: The Job ID in the database
        :type day: str
        :param day: The day in YYYY-MM-DD format
        :type pair: str
        :param pair: the name of the pair (EXAMPLE?)
        :type jobtype: str
        :param jobtype: CrossCorrelation (CC) or dt/t (DTT) Job?
        :type flag: str
        :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.
        """
        __incomplete_tablename__ = "jobs"

        ref = Column(Integer, primary_key=True)
        day = Column(String(10))
        pair = Column(String(23))
        jobtype = Column(String(10))
        flag = Column(String(1))
        lastmod = Column(TIMESTAMP, server_onupdate=text('CURRENT_TIMESTAMP'),
                         server_default=text("CURRENT_TIMESTAMP"))

        table_args__ = (Index('job_index', "day", "pair", "jobtype", unique=True),
                        Index('job_index2', "jobtype", "flag", unique=False))

        def __init__(self, day, pair, jobtype, flag,
                     lastmod=datetime.datetime.utcnow()):
            """"""
            self.day = day
            self.pair = pair
            self.jobtype = jobtype
            self.flag = flag
            self.lastmod = lastmod

    class Station(PrefixerBase):
        """
        Station Object

        :type ref: int
        :param ref: The Station ID in the database
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
        __incomplete_tablename__ = "stations"
        ref = Column(Integer, primary_key=True)
        net = Column(String(10))
        sta = Column(String(10))
        used_location_codes = Column(String(200))
        used_channel_names = Column(String(200))

        X = Column(REAL())
        Y = Column(REAL())
        altitude = Column(Float())
        coordinates = Column(Enum('DEG', 'UTM', name="coordinate_type"))
        used = Column(Boolean)

        def __init__(self, *args):
            """"""
            if len(args):
                self.net = args[0]
                self.sta = args[1]
                self.X = args[2]
                self.Y = args[3]
                self.altitude = args[4]
                self.coordinates = args[5]
                self.used = args[7]

        def locs(self):
            if self.used_location_codes is None:
                location_codes = []
            else:
                location_codes = sorted(self.used_location_codes.split(","))
            return location_codes

        def chans(self):
            if self.used_channel_names is None:
                channels = []
            else:
                channels = sorted(self.used_channel_names.split(","))
            return channels

    ########################################################################

    class Config(PrefixerBase):
        """
        Config Object

        :type name: str
        :param name: The name of the config bit to set.

        :type value: str
        :param value: The value of parameter `name`
        """
        __incomplete_tablename__ = "config"
        name = Column(String(255), primary_key=True)
        value = Column(String(255))      
        used_in = Column(String(255))

        def __init__(self, name, value, used_in):
            """"""
            self.name = name
            self.value = value
            self.used_in = used_in

    ########################################################################

    class DataAvailability(PrefixerBase):
        """
        DataAvailability Object

        :type ref: int
        :param ref: The Station ID in the database
        :type net: str
        :param net: The network code of the Station
        :type sta: str
        :param sta: The station code
        :type chan: str
        :param chan: The component (channel)
        :type path: str
        :param path: The full path to the folder containing the file
        :type file: str
        :param file: The name of the file
        :type starttime: datetime
        :param starttime: Start time of the file
        :type endtime: datetime
        :param endtime: End time of the file
        :type data_duration: float
        :param data_duation: Cumulative duration of available data in the file
        :type gaps_duration: float
        :param gaps_duration: Cumulative duration of gaps in the file
        :type samplerate: float
        :param samplerate: Sample rate of the data in the file (in Hz)
        :type flag: str
        :param flag: The status of the entry: "N"ew, "M"odified or "A"rchive
        """
        __incomplete_tablename__ = "data_availability"
        ref = Column(Integer, primary_key=True, autoincrement=True)
        net = Column(String(10))
        sta = Column(String(10))
        loc = Column(String(10))
        chan = Column(String(20))
        path = Column(String(255))
        file = Column(String(255))
        starttime = Column(DateTime)
        endtime = Column(DateTime)
        data_duration = Column(Float)
        gaps_duration = Column(Float)
        samplerate = Column(Float)
        flag = Column(String(1))
        # UniqueConstraint('net', 'sta', 'comp', 'filename', name='uix_1')
        table_args__ = (Index('da_index',
                              "path",
                              "file",
                              "net",
                              "sta",
                              "loc",
                              "chan", unique=True),)

        def __init__(self, net, sta, loc, chan, path, file, starttime, endtime,
                     data_duration, gaps_duration, samplerate, flag):
            """"""
            self.net = net
            self.sta = sta
            self.loc = loc
            self.chan = chan
            self.path = path
            self.file = file
            self.starttime = starttime
            self.endtime = endtime
            self.data_duration = data_duration
            self.gaps_duration = gaps_duration
            self.samplerate = samplerate
            self.flag = flag

    ########################################################################

    return sqlschema(Base, PrefixerBase,
                     Filter, Job, Station, Config, DataAvailability)
    # end of declare_tables()


# These module objects only use the prefix defined in db.ini.
# They should be re-defined if the prefix is to be changed.
Base, PrefixerBase, Filter, Job, Station, Config, DataAvailability = declare_tables()
