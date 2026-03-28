"""
Refactored MSNoise table definitions with unified configuration system
and improved workflow chain management using PrefixerBase and declare_tables.
"""

import datetime
from collections import namedtuple
from sqlalchemy import Column, Integer, String, Float, Boolean, \
    DateTime, Text, ForeignKey, UniqueConstraint, Index, Enum, REAL, TIMESTAMP, CheckConstraint
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm import relationship, backref  # noqa: F401
from sqlalchemy.sql import text
# from .api import read_prefix

WORKFLOW_CHAINS = {
    'global': {
        'next_steps': ['preprocess', 'psd'],
        'is_entry_point': True,
        'is_terminal': False
    },
    'preprocess': {
        'next_steps': ['cc'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'psd': {
        'next_steps': ['psd_rms'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'psd_rms': {
        'next_steps': [],
        'is_entry_point': False,
        'is_terminal': True
    },
    'cc': {
        'next_steps': ['filter'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'filter': {
        'next_steps': ['stack'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'stack': {
        'next_steps': ['refstack'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'refstack': {
        'next_steps': ['mwcs', 'stretching', 'wavelet'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'mwcs': {
        'next_steps': ['mwcs_dtt'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'mwcs_dtt': {
        'next_steps': ['mwcs_dtt_dvv'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'stretching': {
        'next_steps': ['stretching_dvv'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'wavelet': {
        'next_steps': ['wavelet_dtt'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'wavelet_dtt': {
        'next_steps': ['wavelet_dtt_dvv'],
        'is_entry_point': False,
        'is_terminal': False
    },
    'mwcs_dtt_dvv': {
        'next_steps': [],
        'is_entry_point': False,
        'is_terminal': True
    },
    'stretching_dvv': {
        'next_steps': [],
        'is_entry_point': False,
        'is_terminal': True
    },
    'wavelet_dtt_dvv': {
        'next_steps': [],
        'is_entry_point': False,
        'is_terminal': True
    }
}

# Canonical processing order used for sorting workflow steps in the UI and
# during step creation.  Keep in sync with WORKFLOW_CHAINS above.
WORKFLOW_ORDER = [
    'global',
    'preprocess',
    'cc',
    'psd',
    'psd_rms',
    'filter',
    'stack',
    'refstack',
    'mwcs',
    'mwcs_dtt',
    'mwcs_dtt_dvv',
    'stretching',
    'stretching_dvv',
    'wavelet',
    'wavelet_dtt',
    'wavelet_dtt_dvv',
]


def declare_tables(prefix=None):
    """
    Define classes mapped to relational database tables and return the declared
    classes in a namedtuple.
    """
    # Implementation note: these definitions must be done in a function to
    # allow for the prefix to be evaluated at run time. If this was done at
    # module level, the code would be evaluated at compile time (when the
    # interpreter starts) and the prefix could not be specified at run time.

    # If no prefix is given, read it from the db.ini file
    if prefix is None:
        prefix = ""

    # Define the namedtuple to return
    sqlschema = namedtuple('SQLSchema', ['Base', 'PrefixerBase',
        'Job', 'Station', 'Config', 'DataAvailability', 'WorkflowStep', 'WorkflowLink',
        'Lineage', 'DataSource'])

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

    class Config(PrefixerBase):
        """
        Unified configuration parameter storage for both global and workflow-specific settings.
        
        This class replaces the old separate Config and ConfigSets tables, providing
        a unified approach to configuration management.
        
        :type ref: int
        :param ref: The configuration parameter ID in the database
        :type name: str
        :param name: The parameter name (e.g., 'maxlag', 'dtt_minlag')
        :type category: str
        :param category: The parameter category ('global', 'mwcs', 'mwcs_dtt', 'stretching', etc.)
        :type set_number: int
        :param set_number: Configuration set number (NULL for global, number for workflow sets)
        :type value: str
        :param value: The parameter value as string
        :type param_type: str
        :param param_type: The parameter type ('str', 'int', 'float', 'bool')
        :type default_value: str
        :param default_value: The default value for this parameter
        :type description: str
        :param description: Description of the parameter
        :type units: str
        :param units: Units of measurement for the parameter
        :type possible_values: str
        :param possible_values: slash-separated list of possible values
        """
        
        __incomplete_tablename__ = "config"
        
        ref = Column(Integer, primary_key=True)
        
        # Parameter identification
        name = Column(String(255), nullable=False)
        category = Column(String(255), nullable=False, default='global')
        set_number = Column(Integer, nullable=True)  # NULL for global, number for workflow sets
        
        # Parameter value and type
        value = Column(String(255), nullable=False)
        param_type = Column(String(50), nullable=False, default='str')
        default_value = Column(String(255), nullable=True)
        
        # Documentation and validation
        description = Column(Text, nullable=True)
        units = Column(String(50), nullable=True)
        possible_values = Column(Text, nullable=True)  # slash-separated or JSON

        # Unique constraint ensuring one parameter per name/category/set combination
        __table_args__ = (
            UniqueConstraint('name', 'category', 'set_number', name='_config_param_unique'),
            Index('idx_config_category_set', 'category', 'set_number'),
            Index('idx_config_name', 'name'),
        )
        
        def __init__(self, name=None, value=None, category='global', set_number=None, 
                     param_type='str', description=None, **kwargs):
            """Initialize configuration parameter"""
            self.name = name
            self.value = str(value) if value is not None else ""
            self.category = category
            self.set_number = set_number
            self.param_type = param_type
            self.description = description
            for key, val in kwargs.items():
                setattr(self, key, val)
        
        def get_typed_value(self):
            """Convert string value to the appropriate type"""
            if self.param_type == 'int':
                return int(self.value)
            elif self.param_type == 'float':
                return float(self.value)
            elif self.param_type == 'bool':
                return self.value.lower() in ('true', 'yes', 'y', '1')
            else:
                return self.value
        
        def is_global_config(self):
            """Check if this is a global configuration parameter"""
            return self.category == 'global'
        
        def is_workflow_config(self):
            """Check if this is a workflow-specific configuration parameter"""
            return self.category != 'global'
        
        def validate_value(self):
            """Validate the value against possible_values if defined"""
            if self.possible_values:
                valid_values = [v.strip() for v in self.possible_values.split('/')]
                if self.value not in valid_values:
                    raise ValueError(f"Invalid value '{self.value}' for parameter '{self.name}'. "
                                   f"Valid values are: {valid_values}")
        
        @classmethod
        def get_global_config(cls, session, name):
            """Get a global configuration parameter"""
            return session.query(cls).filter(
                cls.name == name,
                cls.category == 'global'
            ).first()
        
        @classmethod
        def get_workflow_config(cls, session, category, set_number):
            """Get all parameters for a workflow configuration set"""
            return session.query(cls).filter(
                cls.category == category,
                cls.set_number == set_number,
            ).all()
        
        @classmethod
        def get_config_dict(cls, session, category, set_number=None):
            """Get configuration as a dictionary"""
            if set_number is None:
                # Global config
                params = session.query(cls).filter(
                    cls.category == 'global',
                ).all()
            else:
                # Workflow config
                params = cls.get_workflow_config(session, category, set_number)
            
            return {param.name: param.get_typed_value() for param in params}
        
        def __str__(self):
            return f"Config({self.name}={self.value}, category={self.category}, set={self.set_number})"
        
        def __repr__(self):
            return f"<Config(name='{self.name}', category='{self.category}', " \
                   f"set_number={self.set_number}, value='{self.value}')>"

    ########################################################################

    class Lineage(PrefixerBase):
        """Normalised lineage-string table.

        Each distinct lineage path (e.g.
        ``"preprocess_1/cc_1/filter_1/stack_1/refstack_1/mwcs_1/mwcs_dtt_1"``)
        is stored exactly once.  :class:`Job` rows reference it via a small
        integer foreign key instead of repeating the full string across
        potentially millions of rows (~20× storage saving for the column).

        :attr lineage_id: Auto-incremented primary key.
        :attr lineage_str: Slash-separated step-name path, unique.
        """
        __incomplete_tablename__ = "lineages"

        lineage_id  = Column(Integer, primary_key=True, autoincrement=True)
        lineage_str = Column(String(512), nullable=False, unique=True)

        def __init__(self, lineage_str: str):
            self.lineage_str = lineage_str

        def __repr__(self):
            return f"<Lineage({self.lineage_id}: {self.lineage_str!r})>"


    class Job(PrefixerBase):
        """
        Enhanced Job Object with workflow support

        This class maintains backward compatibility while adding workflow awareness.
        Jobs are now linked to specific workflow steps and their associated config sets.

        :type ref: int
        :param ref: The Job ID in the database
        :type day: str
        :param day: The day in YYYY-MM-DD format
        :type pair: str
        :param pair: the name of the pair (STATION1:STATION2)
        :type flag: str
        :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.
        :type step_id: int
        :param step_id: Foreign key to WorkflowStep table
        :type jobtype: str
        :param jobtype: Legacy job type, now derived from step info (for backward compatibility)
        :type priority: int
        :param priority: Job priority (higher number = higher priority)
        """
        __incomplete_tablename__ = "jobs"

        ref = Column(Integer, primary_key=True)
        day = Column(String(10))
        pair = Column(String(23))
        flag = Column(String(1))
        lastmod = Column(TIMESTAMP, server_onupdate=text('CURRENT_TIMESTAMP'),
                         server_default=text("CURRENT_TIMESTAMP"))

        # New workflow-aware fields
        step_id = Column(Integer, ForeignKey(f"{prefix}workflow_steps.step_id"), nullable=False)

        # Normalised lineage FK (String column replaced by integer reference)
        lineage_id  = Column(Integer,
                             ForeignKey(f"{prefix}lineages.lineage_id"),
                             nullable=True, index=True)
        jobtype = Column(String(50))  # Now derived from step info, but kept for compatibility
        priority = Column(Integer, default=0)  # Job priority

        # Relationships
        workflow_step = relationship("WorkflowStep", backref="jobs")
        lineage_ref   = relationship("Lineage",
                                     foreign_keys=[lineage_id],
                                     backref="jobs")

        @property
        def lineage(self):
            """The lineage string, resolved through the Lineage FK."""
            return self.lineage_ref.lineage_str if self.lineage_ref else None

        @lineage.setter
        def lineage(self, value):
            """Set lineage by string.

            If the object is already attached to a session, resolves
            immediately via ``_get_or_create_lineage_id``.  Otherwise
            stores the string in ``_pending_lineage_str`` for resolution
            by the ``before_insert`` event registered at module level.
            """
            if value is None:
                self.lineage_id = None
                self.lineage_ref = None
                self._pending_lineage_str = None
                return
            from sqlalchemy.orm import object_session
            session = object_session(self)
            if session is not None:
                from msnoise.core.workflow import _get_or_create_lineage_id
                with session.no_autoflush:
                    lid = _get_or_create_lineage_id(session, value)
                if lid is not None:
                    self.lineage_id = lid
                else:
                    self._pending_lineage_str = value
            else:
                # Not yet in a session — defer resolution to before_insert
                self._pending_lineage_str = value


        __table_args__ = (
            # Updated unique constraint to include workflow context
            Index('job_index', "day", "pair", "step_id", "lineage_id", unique=True),
            Index('job_index2', "flag", "step_id", "priority", unique=False),
            # Legacy index for backward compatibility
            Index('job_legacy_index', "day", "pair", "jobtype", unique=False),
        )

        def __init__(self, day=None, pair=None, flag=None,
                     step_id=None, jobtype=None, priority=0, lastmod=None, **kwargs):
            """Initialize job with workflow support"""
            self.day = day
            self.pair = pair
            self.flag = flag
            self.step_id = step_id
            self.priority = priority
            self.lastmod = lastmod or datetime.datetime.now(datetime.timezone.utc)

            # Handle jobtype - derive from step if not provided
            if jobtype is not None:
                self.jobtype = jobtype
            elif step_id is not None:
                # Will be set via property or method once step relationship is loaded
                self.jobtype = None
            else:
                self.jobtype = None

            for key, val in kwargs.items():
                setattr(self, key, val)

        @property
        def derived_jobtype(self):
            """Derive jobtype from workflow step information"""
            if self.workflow_step:
                return f"{self.workflow_step.step_name}_{self.workflow_step.set_number}"
            return self.jobtype

        @property
        def config_category(self):
            """Get the config category for this job"""
            if self.workflow_step:
                return self.workflow_step.category
            return None

        @property
        def config_set_number(self):
            """Get the config set number for this job"""
            if self.workflow_step:
                return self.workflow_step.set_number
            return None

        @property
        def step_name(self):
            """Get the step name for this job"""
            if self.workflow_step:
                return self.workflow_step.step_name
            return None

        def __str__(self):
            return f"Job({self.day}, {self.pair}, {self.derived_jobtype or self.jobtype}, {self.flag})"

        def __repr__(self):
            return f"<Job(day='{self.day}', pair='{self.pair}', " \
                   f"step_id={self.step_id}, flag='{self.flag}')>"

    ########################################################################

    class DataSource(PrefixerBase):
        """
        DataSource Object — defines where raw waveform data comes from.

        :type ref: int
        :param ref: Primary key
        :type name: str
        :param name: Human label e.g. "local", "IRIS", "GEOFON", "EIDA"
        :type uri: str
        :param uri: Data location. Schemes:
            - bare path or ``sds:///path`` → local SDS archive
            - ``fdsn://http://...``         → FDSN web service
            - ``eida://http://...``         → EIDA routing client
        :type data_structure: str
        :param data_structure: SDS sub-path format for local/SDS sources
            (e.g. "SDS", "BUD"). Ignored for FDSN/EIDA.
        :type auth_env: str
        :param auth_env: Environment variable prefix for credentials.
            Worker looks up ``{auth_env}_FDSN_USER``,
            ``{auth_env}_FDSN_PASSWORD``, ``{auth_env}_FDSN_TOKEN``.
            Default "MSNOISE".
        """
        __incomplete_tablename__ = "data_sources"

        ref = Column(Integer, primary_key=True, autoincrement=True)
        name = Column(String(64), nullable=False, unique=True)
        uri = Column(String(512), nullable=False, default="")
        data_structure = Column(String(64), nullable=False, default="SDS")
        auth_env = Column(String(64), nullable=False, default="MSNOISE")
        archive_format = Column(String(32), nullable=True, default=None)
        # Force ObsPy read format for local/SDS sources (None = auto-detect)

        __table_args__ = (
            Index('idx_datasource_name', 'name'),
        )

        def __init__(self, name=None, uri="", data_structure="SDS",
                     auth_env="MSNOISE", archive_format=None, **kwargs):
            self.name = name
            self.uri = uri
            self.data_structure = data_structure
            self.auth_env = auth_env
            self.archive_format = archive_format
            for key, val in kwargs.items():
                setattr(self, key, val)

        def __str__(self):
            return f"DataSource({self.name!r}, uri={self.uri!r})"

        def __repr__(self):
            return (f"DataSource(name={self.name!r}, uri={self.uri!r}, "
                    f"data_structure={self.data_structure!r}, "
                    f"auth_env={self.auth_env!r}, "
                    f"archive_format={self.archive_format!r})")

    ########################################################################

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
        :type used: bool
        :param used: Whether this station must be used in the computations.
        """
        __incomplete_tablename__ = "stations"
        
        ref = Column(Integer, primary_key=True)
        net = Column(String(10))
        sta = Column(String(10))
        used_location_codes = Column(String(200))
        used_channel_names = Column(String(200))
        X = Column(REAL)
        Y = Column(REAL)
        altitude = Column(Float)
        coordinates = Column(Enum('DEG', 'UTM', name="coordinate_type"))
        used = Column(Boolean, default=True)
        data_source_id = Column(Integer, ForeignKey(f"{prefix}data_sources.ref"), nullable=True)
        # NULL → use project default DataSource (id=1)

        __table_args__ = (
            UniqueConstraint('net', 'sta', name='_station_unique'),
            Index('idx_station_datasource', 'data_source_id'),
            Index('idx_station_net', 'net'),
            Index('idx_station_used', 'used'),
        )

        def __init__(self, net=None, sta=None, X=None, Y=None, altitude=None,
                     coordinates=None, used=True, data_source_id=None, **kwargs):
            """Initialize station"""
            self.net = net
            self.sta = sta
            self.X = X
            self.Y = Y
            self.altitude = altitude
            self.coordinates = coordinates
            self.used = used
            self.data_source_id = data_source_id
            for key, val in kwargs.items():
                setattr(self, key, val)

        def __str__(self):
            return f"Station({self.net}.{self.sta})"

        def locs(self):
            """Get list of location codes"""
            if self.used_location_codes is None:
                return []
            return sorted(self.used_location_codes.split(","))

        def chans(self):
            """Get list of channel names"""
            if self.used_channel_names is None:
                return []
            return sorted(self.used_channel_names.split(","))

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
        path = Column(String(512))
        # Path RELATIVE to DataSource.uri — reconstruct full path at runtime
        file = Column(String(255))
        data_source_id = Column(Integer, ForeignKey(f"{prefix}data_sources.ref"), nullable=True)
        starttime = Column(DateTime)
        endtime = Column(DateTime)
        data_duration = Column(Float)
        gaps_duration = Column(Float)
        samplerate = Column(Float)
        flag = Column(String(1))

        __table_args__ = (
            Index('da_index', "path", "file", "net", "sta", "loc", "chan", unique=True),
            Index('idx_da_net_sta', 'net', 'sta'),
            Index('idx_da_time', 'starttime', 'endtime'),
            Index('idx_da_flag', 'flag'),
            Index('idx_da_datasource', 'data_source_id'),
        )

        def __init__(self, net=None, sta=None, loc=None, chan=None, path=None, file=None,
                     starttime=None, endtime=None, data_duration=None, gaps_duration=None,
                     samplerate=None, flag=None, data_source_id=None, **kwargs):
            """Initialize data availability record"""
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
            self.data_source_id = data_source_id
            for key, val in kwargs.items():
                setattr(self, key, val)

        def __str__(self):
            return f"DataAvailability({self.net}.{self.sta}.{self.loc}.{self.chan})"

    ########################################################################

    class WorkflowStep(Base):
        """
        Simple workflow step container - just wraps a config category+set_number
        """
        __tablename__ = f"{prefix}workflow_steps"

        # Primary key
        step_id = Column(Integer, primary_key=True, autoincrement=True)

        # Step identification
        step_name = Column(String(50), nullable=False)

        # Reference to config set
        category = Column(String(20), nullable=False)
        set_number = Column(Integer, nullable=False)

        # Metadata
        description = Column(String(200), nullable=True)
        is_active = Column(Boolean, default=True)
        created_at = Column(DateTime, default=lambda: datetime.datetime.now(datetime.timezone.utc))
        updated_at = Column(DateTime, default=lambda: datetime.datetime.now(datetime.timezone.utc), onupdate=lambda: datetime.datetime.now(datetime.timezone.utc))

        # Constraints
        __table_args__ = (
            UniqueConstraint('step_name', name='unique_step_name'),
            UniqueConstraint('category', 'set_number', name='unique_config_per_category'),
        )

        def __repr__(self):
            return f"<WorkflowStep('{self.step_name}', {self.category}:{self.set_number})>"

        def get_config_params(self, session):
            """Get configuration parameters for this step"""
            return session.query(Config).filter(
                Config.category == self.category,
                Config.set_number == self.set_number
            ).all()

    class WorkflowLink(Base):
        """
        Links between workflow steps - supports all relationship types
        """
        __tablename__ = f"{prefix}workflow_links"

        # Primary key
        link_id = Column(Integer, primary_key=True, autoincrement=True)

        # Relationship definition
        from_step_id = Column(Integer, ForeignKey(f"{prefix}workflow_steps.step_id"), nullable=False)
        to_step_id = Column(Integer, ForeignKey(f"{prefix}workflow_steps.step_id"), nullable=False)

        # Metadata
        link_type = Column(String(20), default="default")  # For different types of connections
        is_active = Column(Boolean, default=True)
        created_at = Column(DateTime, default=lambda: datetime.datetime.now(datetime.timezone.utc))

        # Relationships
        from_step = relationship("WorkflowStep", foreign_keys=[from_step_id], backref="outgoing_links")
        to_step = relationship("WorkflowStep", foreign_keys=[to_step_id], backref="incoming_links")

        # Constraints
        __table_args__ = (
            UniqueConstraint('from_step_id', 'to_step_id', name='unique_link'),
            CheckConstraint('from_step_id != to_step_id', name='no_self_links'),
        )

        def __repr__(self):
            return f"<WorkflowLink({self.from_step_id} -> {self.to_step_id})>"

    # Register before_insert on THIS call's Job class so the event
    # always targets the class actually used by active sessions.
    from sqlalchemy import event as _event, text as _text

    @_event.listens_for(Job, "before_insert", propagate=True)
    def _resolve_pending_lineage(mapper, connection, target):
        """Resolve _pending_lineage_str → lineage_id via raw SQL before INSERT."""
        pending = getattr(target, "_pending_lineage_str", None)
        if not pending:
            return
        # Determine lineages table name (handles prefix)
        lin_table = mapper.persist_selectable.metadata.tables.get(
            f"{prefix}lineages" if prefix else "lineages")
        if lin_table is None:
            return
        tname = lin_table.name
        row = connection.execute(
            _text(f"SELECT lineage_id FROM {tname} WHERE lineage_str = :s"),
            {"s": pending}
        ).fetchone()
        if row is not None:
            target.lineage_id = row[0]
        else:
            result = connection.execute(
                _text(f"INSERT INTO {tname} (lineage_str) VALUES (:s)"),
                {"s": pending}
            )
            target.lineage_id = result.lastrowid
        target._pending_lineage_str = None

    # Return the schema namedtuple
    return sqlschema(Base, PrefixerBase,
                     Job, Station, Config, DataAvailability, WorkflowStep, WorkflowLink,
                     Lineage, DataSource)
    # end of declare_tables()

Base, PrefixerBase, Job, Station, Config, DataAvailability, WorkflowStep, WorkflowLink, Lineage, DataSource = declare_tables()







# Helper functions for configuration management
def get_config(session, name, category='global', set_number=None):
    """Get a configuration parameter value"""
    # This function needs to be updated to work with the declare_tables pattern
    schema = declare_tables()
    config_dict = schema.Config.get_config_dict(session, category, set_number)
    return config_dict.get(name)


def set_config(session, name, value, category='global', set_number=None, 
               param_type='str', description=None):
    """Set a configuration parameter"""
    # This function needs to be updated to work with the declare_tables pattern
    schema = declare_tables()
    Config = schema.Config
    
    param = session.query(Config).filter(
        Config.name == name,
        Config.category == category,
        Config.set_number == set_number
    ).first()
    
    if not param:
        param = Config(
            name=name,
            category=category,
            set_number=set_number,
            param_type=param_type,
            description=description
        )
        session.add(param)
    
    param.value = str(value)
    param.updated_at = datetime.datetime.now(datetime.timezone.utc)
    param.validate_value()
    session.commit()
    return param

