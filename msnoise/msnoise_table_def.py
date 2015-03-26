from sqlalchemy import Column, Integer, String, Float, Boolean, DateTime,\
    text, TIMESTAMP, Enum
from sqlalchemy.ext.declarative import declarative_base
import datetime

Base = declarative_base()


class Filter(Base):
    """
    Filter base class.

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
        MWCS (in Hz)
    :type rms_threshold: float
    :param rms_threshold: Not used anymore
    :type mwcs_wlen: float
    :param mwcs_wlen: Window length (in seconds) to perform MWCS
    :type mwcs_step: float
    :param mwcs_step: Step (in seconds) of the windowing procedure in MWCS
    :type used: bool
    :param used: Is the filter activated for the processing
    """

    __tablename__ = "filters"

    ref = Column(Integer, primary_key=True)
    low = Column(Float)
    mwcs_low = Column(Float)
    high = Column(Float)
    mwcs_high = Column(Float)
    rms_threshold = Column(Float)
    mwcs_wlen = Column(Float)
    mwcs_step = Column(Float)
    used = Column(Boolean)

    def __init__(self, low, mwcs_low, high, mwcs_high, rms_threshold,
                 mwcs_wlen, mwcs_step, used):
        """"""
        self.low = low
        self.mwcs_low = mwcs_low
        self.high = high
        self.mwcs_high = mwcs_high
        self.rms_threshold = rms_threshold
        self.mwcs_wlen = mwcs_wlen
        self.mwcs_step = mwcs_step
        self.used = used


class Job(Base):
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
    __tablename__ = "jobs"

    ref = Column(Integer, primary_key=True)
    day = Column(String(10))
    pair = Column(String(20))
    jobtype = Column(String(10))
    flag = Column(String(1))
    lastmod = Column(TIMESTAMP, server_onupdate=text('CURRENT_TIMESTAMP'))

    def __init__(self, day, pair, jobtype, flag,
                 lastmod=datetime.datetime.utcnow()):
        """"""
        self.day = day
        self.pair = pair
        self.jobtype = jobtype
        self.flag = flag
        self.lastmod = lastmod


class Station(Base):
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
    __tablename__ = "stations"
    ref = Column(Integer, primary_key=True)
    net = Column(String(10))
    sta = Column(String(10))
    X = Column(Float)
    Y = Column(Float)
    altitude = Column(Float)
    coordinates = Column(Enum('DEG', 'UTM'))
    instrument = Column(String(20))
    used = Column(Boolean)

    def __init__(self, net, sta, X, Y, altitude, coordinates, instrument,
                 used):
        """"""
        self.net = net
        self.sta = sta
        self.X = X
        self.Y = Y
        self.altitude = altitude
        self.coordinates = coordinates
        self.instrument = instrument
        self.used = used

########################################################################


class Config(Base):
    """
    Config Object

    :type name: str
    :param name: The name of the config bit to set.

    :type value: str
    :param value: The value of parameter `name`
    """
    __tablename__ = "config"
    name = Column(String(255), primary_key=True)
    value = Column(String(255))

    def __init__(self, name, value):
        """"""
        self.name = name
        self.value = value

########################################################################


class DataAvailability(Base):
    """
    DataAvailability Object

    :type ref: int
    :param ref: The Station ID in the database
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
    __tablename__ = "data_availability"
    ref = Column(Integer, primary_key=True, autoincrement=True)
    net = Column(String(10))
    sta = Column(String(10))
    comp = Column(String(20))
    path = Column(String(255))
    file = Column(String(255))
    starttime = Column(DateTime)
    endtime = Column(DateTime)
    data_duration = Column(Float)
    gaps_duration = Column(Float)
    samplerate = Column(Float)
    flag = Column(String(1))

    def __init__(self, net, sta, comp, path, file, starttime, endtime,
                 data_duration, gaps_duration, samplerate, flag):
        """"""
        self.net = net
        self.sta = sta
        self.comp = comp
        self.path = path
        self.file = file
        self.starttime = starttime
        self.endtime = endtime
        self.data_duration = data_duration
        self.gaps_duration = gaps_duration
        self.samplerate = samplerate
        self.flag = flag
