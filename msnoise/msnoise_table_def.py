"""
SQLAlchemy table definition.
"""

import datetime
import os
from collections import namedtuple
from sqlalchemy import Column, Integer, String, Float, Boolean, DateTime,\
    text, TIMESTAMP, Enum, REAL, UniqueConstraint, Index, ForeignKey, Table

from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm import relationship

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
        'Job', 'Station', 'Config', 'DataAvailability', 'DvvMwcs', 'DvvMwcsDtt', 
        'DvvStretching', 'DvvWct', 'DvvWctDtt', 'filter_mwcs_assoc', 'mwcs_dtt_assoc',
        'filter_stretching_assoc', 'filter_wct_assoc', 'wct_dtt_assoc'])

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
    # **Association Tables**
    
    filter_mwcs_assoc = Table(
        "filter_mwcs_assoc", PrefixerBase.metadata,
        Column("dvv_mwcs_ref", Integer, ForeignKey("dvv_mwcs.ref"), primary_key=True),
        Column("filt_ref", Integer, ForeignKey("filters.ref"), primary_key=True)
    )

    mwcs_dtt_assoc = Table(
        "mwcs_dtt_assoc", PrefixerBase.metadata,
        Column("dvv_mwcs_dtt_ref", Integer, ForeignKey("dvv_mwcs_dtt.ref"), primary_key=True),
        Column("dvv_mwcs_ref", Integer, ForeignKey("dvv_mwcs.ref"), primary_key=True)
    )

    filter_stretching_assoc = Table(
        "filter_stretching_assoc", PrefixerBase.metadata,
        Column("dvv_stretching_ref", Integer, ForeignKey("dvv_stretching.ref"), primary_key=True),
        Column("filt_ref", Integer, ForeignKey("filters.ref"), primary_key=True)
    )

    filter_wct_assoc = Table(
        "filter_wct_assoc", PrefixerBase.metadata,
        Column("dvv_wct_ref", Integer, ForeignKey("dvv_wct.ref"), primary_key=True),
        Column("filt_ref", Integer, ForeignKey("filters.ref"), primary_key=True)
    )

    wct_dtt_assoc = Table(
        "wct_dtt_assoc", PrefixerBase.metadata,
        Column("dvv_wct_dtt_ref", Integer, ForeignKey("dvv_wct_dtt.ref"), primary_key=True),
        Column("dvv_wct_ref", Integer, ForeignKey("dvv_wct.ref"), primary_key=True)
    )

    ########################################################################

    class DvvMwcs(PrefixerBase):
        """
        Dvv_mwcs base class.

        :type ref: int
        :param ref: The id of the MWCS_params in the database
        :type filt_ref: int
        :param filt_ref: The id of the a filter in the filter table
        :type freqmin: float
        :param freqmin: The lower frequency bound to apply MWCS (in Hz)
        :type freqmax: float
        :param : The upper frequency bound to apply MWCS (in Hz)        
        :type mwcs_wlen: float
        :param mwcs_wlen: Window length (in seconds) to perform MWCS
        :type mwcs_step: float
        :param mwcs_step: Step (in seconds) of the windowing procedure in MWCS
        :type used: bool
        :param used: Is the parameter set activated for the processing
        """

        __incomplete_tablename__ = "dvv_mwcs"

        ref = Column(Integer, primary_key=True)
        freqmin = Column(Float())
        freqmax = Column(Float())
        mwcs_wlen = Column(Float())
        mwcs_step = Column(Float())
        used = Column(Boolean(), default=True)

        # Many-to-Many relationship with Filters
        filters = relationship("Filter", secondary=filter_mwcs_assoc, back_populates="mwcs_params")

        # Many-to-Many relationship with MWCS DTT settings
        dtt_params = relationship("DvvMwcsDtt", secondary=mwcs_dtt_assoc, back_populates="mwcs_params")

        def __str__(self):
            return f"MWCS Params {self.ref} ({self.freqmin}-{self.freqmax} Hz, mwcs_wlen:{self.mwcs_wlen}, mwcs_step:{self.mwcs_step})"

    class DvvMwcsDtt(PrefixerBase):
        """
        Dvv_mwcs_dtt base class.

        :type ref: int
        :param ref: The id of the dtt params for mwcs in the database
        :type dvv_mwcs_ref: int
        :param dvv_mwcs_ref: The id of mwcs parameters from dvv_mwcs table
        :type dtt_minlag: float
        :param dtt_minlag: If ``dtt_lag`` =static (in config table): min lag time (in seconds)
        :type dtt_width: float
        :param dtt_width: Width of the time lag window (in seconds)
        :type dtt_lag: string
        :params dtt_lag: How is the lag window defined for MWCS [static]/dynamic.
        :type dtt_v: float
        :param dtt_v: If ``dttlag`` =dynamic (in config table): what velocity to use to avoid ballistic waves [1.0] km/s (default=1.0)
        :type dtt_sides: string
        :params dtt_sides: Which sides to use,str,both/left/right
        :type dtt_mincoh: float
        :params dtt_mincoh: Minimum coherence on dt measurement, MWCS points with values lower than that will not be used in the WLS, [0:1] 
        :type dtt_maxerr: float
        :params dtt_maxerr: Maximum error on dt measurement, MWCS points with values larger than that will not be used in the WLS [0:1]
        :type dtt_maxdt: float
        :params dtt_maxdt: Maximum dt values, MWCS points with values larger than that will not be used in the WLS (in seconds)
        :type used: bool
        :param used: Is the parameter set activated for the processing
        """

        __incomplete_tablename__ = "dvv_mwcs_dtt"

        ref = Column(Integer, primary_key=True)
        dtt_minlag = Column(Float())
        dtt_width = Column(Float())
        dtt_lag = Column(String(255))
        dtt_v = Column(Float())
        dtt_sides = Column(String(255))
        dtt_mincoh = Column(Float())
        dtt_maxerr = Column(Float())
        dtt_maxdt = Column(Float())
        used = Column(Boolean(), default=True)

        # Many-to-Many relationship with MWCS parameter sets
        mwcs_params = relationship("DvvMwcs", secondary=mwcs_dtt_assoc, back_populates="dtt_params")

    class DvvStretching(PrefixerBase):
        """
        Dvv_stretching base class.

        :type ref: int
        :param ref: The ID of the stretching parameters in the database
        :type filt_ref: int
        :param filt_ref: The ID of a filter in the filters table
        :type stretching_minlag: float
        :param stretching_minlag: Minimum lag time for stretching analysis (in seconds)
        :type stretching_width: float
        :param stretching_width: Width of the time lag window (in seconds)
        :type stretching_lag: string
        :param stretching_lag: How the lag window is defined for stretching [static]/dynamic
        :type stretching_v: float
        :param stretching_v: If ``stretching_lag`` = dynamic, what velocity to use to avoid ballistic waves [1.0] km/s
        :type stretching_sides: string
        :param stretching_sides: Which sides to use, options: both/left/right
        :type stretching_max: float
        :param stretching_max: Maximum stretching coefficient, e.g. 0.5 = 50%, 0.01 = 1%
        :type stretching_nsteps: int
        :param stretching_nsteps: Number of stretching steps between 1-``stretching_max`` and 1+``stretching_max``
        :type used: bool
        :param used: Is the parameter set activated for the processing
        """

        __incomplete_tablename__ = "dvv_stretching"

        ref = Column(Integer, primary_key=True)
        filt_ref = Column(Integer, ForeignKey('filters.ref'))
        stretching_minlag = Column(Float())
        stretching_width = Column(Float())
        stretching_lag = Column(String(255))
        stretching_v = Column(Float())
        stretching_sides = Column(String(255))
        stretching_max = Column(Float())
        stretching_nsteps = Column(Integer())
        used = Column(Boolean(), default=True)

        # Many-to-Many relationship with Filters
        filters = relationship("Filter", secondary=filter_stretching_assoc, back_populates="stretching_params")

    class DvvWct(PrefixerBase):
        """
        Dvv_wct base class.

        :type ref: int
        :param ref: The id of the WCT parameters in the database
        :type filt_ref: int
        :param filt_ref: The id of a filter in the filter table
        :type wct_ns: float
        :param wct_ns: Smoothing parameter in frequency
        :type wct_nt: float
        :param wct_nt: Smoothing parameter in time
        :type wct_vpo: float
        :param wct_vpo: Spacing parameter between discrete scales
        :type wct_nptsfreq: int
        :param wct_nptsfreq: Number of frequency points between min and max
        :type wct_norm: bool
        :param wct_norm: If the REF and CCF are normalized before computing wavelet? [Y]/N
        :type wavelet_type: string
        :param wavelet_type: Type of wavelet function used (e.g., Morlet, Paul, DOG, MexicanHat)
        :type used: bool
        :param used: Is the parameter set activated for the processing?
        """

        __incomplete_tablename__ = "dvv_wct"

        ref = Column(Integer, primary_key=True)
        filt_ref = Column(Integer, ForeignKey('filters.ref'))
        wct_ns = Column(Float())
        wct_nt = Column(Float())
        wct_vpo = Column(Float())
        wct_nptsfreq = Column(Integer())
        wct_norm = Column(Boolean(), default=True)
        wavelet_type = Column(String(255))
        used = Column(Boolean(), default=True)

        # Many-to-Many relationship with Filters
        filters = relationship("Filter", secondary=filter_wct_assoc, back_populates="wct_params")

        # Many-to-Many relationship with WCT DTT settings
        dtt_params = relationship("DvvWctDtt", secondary=wct_dtt_assoc, back_populates="wct_params")

    class DvvWctDtt(PrefixerBase):
        """
        Dvv_wct_dtt base class.

        :type ref: int
        :param ref: The id of the dtt params for wct in the database
        :type dvv_wct_ref: int
        :param dvv_wct_ref: The id of wct parameters from dvv_wct table
        :type wct_minlag: float
        :param wct_minlag: Minimum lag time (in seconds)
        :type wct_width: float
        :param wct_width: Width of the time lag window (in seconds)
        :type wct_lag: string
        :param wct_lag: How is the lag window defined for WCT [static]/dynamic.
        :type wct_v: float
        :param wct_v: Velocity parameter used for dynamic lag calculation
        :type wct_sides: string
        :param wct_sides: Which sides to use (both/left/right)
        :type wct_mincoh: float
        :param wct_mincoh: Minimum coherence on dt measurement
        :type wct_maxdt: float
        :param wct_maxdt: Maximum dt values (in seconds)
        :type wct_codacycles: int
        :param wct_codacycles: Number of cycles of period (1/freq) between lag_min and lag_max
        :type wct_min_nonzero: float
        :param wct_min_nonzero: Percentage of data points with non-zero weighting required for regression
        :type used: bool
        :param used: Is the parameter set activated for the processing?
        """

        __incomplete_tablename__ = "dvv_wct_dtt"

        ref = Column(Integer, primary_key=True)
        dvv_wct_ref = Column(Integer, ForeignKey('dvv_wct.ref'))
        wct_minlag = Column(Float())
        wct_width = Column(Float())
        wct_lag = Column(String(255))
        wct_v = Column(Float())
        wct_sides = Column(String(255))
        wct_mincoh = Column(Float())
        wct_maxdt = Column(Float())
        wct_codacycles = Column(Integer())
        wct_min_nonzero = Column(Float())
        used = Column(Boolean(), default=True)

        # Many-to-Many relationship with MWCS parameter sets
        wct_params = relationship("DvvWct", secondary=wct_dtt_assoc, back_populates="dtt_params")

    class Filter(PrefixerBase):
        """
        Filter base class.

        :type ref: int
        :param ref: The id of the Filter in the database
        :type freqmin: float
        :param freqmin: The lower frequency bound of the Whiten function (in Hz)
        :type freqmax: float
        :param : The upper frequency bound of the Whiten function (in Hz)
        :type CC: bool
        :param CC: Compute cross-correlation functions between different pairs?
        :type SC: bool
        :param SC: Compute cross-correlation functions between different components of single station?
        :type AC: bool
        :param AC: Compute auto-correlation functions from single station components?        
        :type used: bool
        :param used: Is the filter activated for the processing?
        """

        __incomplete_tablename__ = "filters"

        ref = Column(Integer, primary_key=True)
        freqmin = Column(Float())
        freqmax = Column(Float())
        CC = Column(Boolean(True))
        SC = Column(Boolean(True))
        AC = Column(Boolean(True))
        used = Column(Boolean(True))

        # Many-to-Many relationship with MWCS parameter sets
        mwcs_params = relationship("DvvMwcs", secondary=filter_mwcs_assoc, back_populates="filters")
        stretching_params = relationship("DvvStretching", secondary=filter_stretching_assoc, back_populates="filters")
        wct_params = relationship("DvvWct", secondary=filter_wct_assoc, back_populates="filters")

        def __str__(self):
            return f"Filter {self.ref} ({self.freqmin}-{self.freqmax} Hz)"

        def __init__(self, **kwargs):
            """"""
            # self.freqmin = freqmin
            # self.freqmax = freqmax
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
                     Filter, Job, Station, Config, DataAvailability, DvvMwcs, DvvMwcsDtt, DvvStretching, DvvWct, DvvWctDtt,
                     filter_mwcs_assoc, mwcs_dtt_assoc, filter_stretching_assoc, filter_wct_assoc, wct_dtt_assoc)
    # end of declare_tables()


# These module objects only use the prefix defined in db.ini.
# They should be re-defined if the prefix is to be changed.
Base, PrefixerBase, Filter, Job, Station, Config, DataAvailability, DvvMwcs, DvvMwcsDtt, DvvStretching, DvvWct, DvvWctDtt, \
     filter_mwcs_assoc, mwcs_dtt_assoc, filter_stretching_assoc, filter_wct_assoc, wct_dtt_assoc = declare_tables()
