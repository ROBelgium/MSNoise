"""
Refactored MSNoise table definitions with unified configuration system
and improved workflow chain management using PrefixerBase and declare_tables.
"""

import datetime
from collections import namedtuple
from sqlalchemy import create_engine, Column, Integer, String, Float, Boolean, \
    DateTime, Text, ForeignKey, Table, UniqueConstraint, Index, Enum, REAL, TIMESTAMP
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm import relationship, sessionmaker, backref
from sqlalchemy.sql import func, text
# from .api import read_prefix


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
    sqlschema = namedtuple('SQLSchema', ['Base', 'PrefixerBase', 'Filter',
        'Job', 'Station', 'Config', 'DataAvailability', 'WorkflowSteps',
        'filter_workflow_assoc'])

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
    
    filter_workflow_assoc = Table(
        "filter_workflow_assoc", PrefixerBase.metadata,
        Column("workflow_step_ref", Integer, ForeignKey("workflow_steps.ref"), primary_key=True),
        Column("filt_ref", Integer, ForeignKey("filters.ref"), primary_key=True)
    )

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
        :type used: bool
        :param used: Whether this parameter is active
        :type used_in: str
        :param used_in: Which processing steps use this parameter
        :type created_at: datetime
        :param created_at: When this parameter was created
        :type updated_at: datetime
        :param updated_at: When this parameter was last updated
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

        # Usage metadata
        used = Column(Boolean, nullable=False, default=True)
        used_in = Column(String(255), nullable=True)  # Which processing steps use this

        # Unique constraint ensuring one parameter per name/category/set combination
        __table_args__ = (
            UniqueConstraint('name', 'category', 'set_number', name='_config_param_unique'),
            Index('idx_config_category_set', 'category', 'set_number'),
            Index('idx_config_name', 'name'),
            Index('idx_config_used', 'used'),
        )
        
        def __init__(self, name=None, value=None, category='global', set_number=None, 
                     param_type='str', description=None, used_in=None, **kwargs):
            """Initialize configuration parameter"""
            self.name = name
            self.value = str(value) if value is not None else ""
            self.category = category
            self.set_number = set_number
            self.param_type = param_type
            self.description = description
            self.used_in = used_in
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
                cls.used == True
            ).all()
        
        @classmethod
        def get_config_dict(cls, session, category, set_number=None):
            """Get configuration as a dictionary"""
            if set_number is None:
                # Global config
                params = session.query(cls).filter(
                    cls.category == 'global',
                    cls.used == True
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

    class WorkflowSteps(PrefixerBase):
        """
        WorkflowSteps base class - Defines workflow processing steps that can be chained.

        :type ref: int
        :param ref: The id of the workflow step in the database
        :type step_type: str
        :param step_type: The type of workflow step (mwcs, mwcs_dtt, stretching, wavelet, etc.)
        :type config_set_number: int
        :param config_set_number: Reference to the configuration set number
        :type parent_step_ref: int
        :param parent_step_ref: Reference to the parent workflow step (for chaining)
        :type name: str
        :param name: User-friendly name for this workflow step
        :type description: str
        :param description: Description of this workflow step
        :type chain_order: int
        :param chain_order: Order in the processing chain (0 for root steps)
        :type used: bool
        :param used: Is this workflow step activated for processing
        :type created_at: datetime
        :param created_at: When this workflow step was created
        :type updated_at: datetime
        :param updated_at: When this workflow step was last updated
        """

        __incomplete_tablename__ = "workflow_steps"

        ref = Column(Integer, primary_key=True)
        step_type = Column(String(50), nullable=False)
        config_set_number = Column(Integer, nullable=False)
        parent_step_ref = Column(Integer, ForeignKey('workflow_steps.ref'), nullable=True)
        name = Column(String(255))
        description = Column(String(500))
        chain_order = Column(Integer, default=0)
        used = Column(Boolean, default=True)
        created_at = Column(DateTime, default=datetime.datetime.utcnow)
        updated_at = Column(DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow)

        # Self-referential relationship for chaining
        parent_step = relationship("WorkflowSteps", remote_side=[ref], backref="child_steps")

        # Many-to-Many relationship with Filters (only for root steps that connect to filters)
        filters = relationship("Filter", secondary=filter_workflow_assoc, back_populates="workflow_steps",
                               info={'description': 'The filters used for this workflow step',
                                     'note': "Select one or more filters to apply this step to"})

        __table_args__ = (
            Index('idx_workflow_steps_type', 'step_type'),
            Index('idx_workflow_steps_config', 'config_set_number'),
            Index('idx_workflow_steps_parent', 'parent_step_ref'),
            Index('idx_workflow_steps_used', 'used'),
            Index('idx_workflow_steps_order', 'chain_order'),
        )

        # Valid workflow chains - defines what step types can follow what
        VALID_CHAINS = {
            'preprocess': ['qc', 'cc'],
            'qc': [], # Terminal step
            'cc': ['filter'],
            'filter': ['stack'],
            'stack': ['mwcs', 'stretching', 'wavelet'],
            'mwcs': ['mwcs_dtt'],
            'mwcs_dtt': [],  # Terminal step
            'stretching': [],  # Terminal step
            'wavelet': ['wavelet_dtt'],
            'wavelet_dtt': [],  # Terminal step
        }

        # Configuration type mapping - maps step types to their config categories
        CONFIG_TYPE_MAP = {
            'preprocess': 'preprocess',
            'cc': 'cc',
            'qc': 'qc',
            'filter': 'filter',
            'stack': 'stack',
            'mwcs': 'mwcs',
            'mwcs_dtt': 'mwcs_dtt',
            'stretching': 'stretching',
            'wavelet': 'wavelet',
            'wavelet_dtt': 'wavelet_dtt',
            'wct': 'wct',
            'wct_dtt': 'wct_dtt',
        }

        def __init__(self, step_type=None, config_set_number=None, parent_step_ref=None,
                     name=None, description=None, used=True, **kwargs):
            """Initialize workflow step"""
            self.step_type = step_type
            self.config_set_number = config_set_number
            self.parent_step_ref = parent_step_ref
            self.name = name
            self.description = description
            self.used = used

            # Set chain order based on parent
            if parent_step_ref is None:
                self.chain_order = 0
            else:
                # This will be set properly when parent relationship is established
                self.chain_order = 1
            
            for key, val in kwargs.items():
                setattr(self, key, val)

        def __str__(self):
            return f"WorkflowStep {self.ref} ({self.step_type}): {self.name or 'Unnamed'}"

        def is_root_step(self):
            """Check if this is a root step (no parent)"""
            return self.parent_step_ref is None

        def is_terminal_step(self):
            """Check if this is a terminal step (cannot have children)"""
            return len(self.VALID_CHAINS.get(self.step_type, [])) == 0

        def can_have_child(self, child_step_type):
            """Check if this step can have a child of given type"""
            return child_step_type in self.VALID_CHAINS.get(self.step_type, [])

        def get_valid_child_types(self):
            """Get list of step types that can be added as children"""
            return self.VALID_CHAINS.get(self.step_type, [])

        def validate_chain(self, session):
            """Validate that this step can be chained to its parent"""
            if self.parent_step_ref is None:
                # Root step - must be connected to filters
                if not self.filters:
                    raise ValueError("Root workflow steps must be connected to at least one filter")
                return True
            
            parent = session.query(WorkflowSteps).filter(WorkflowSteps.ref == self.parent_step_ref).first()
            if parent is None:
                raise ValueError(f"Parent step {self.parent_step_ref} not found")
            
            # Check if this step type is valid for the parent
            if not parent.can_have_child(self.step_type):
                raise ValueError(f"Invalid chain: {parent.step_type} cannot be followed by {self.step_type}")
            
            # Child steps should not have direct filter connections
            if self.filters:
                raise ValueError("Child workflow steps should not have direct filter connections")
            
            return True

        def get_root_step(self, session):
            """Get the root step of this chain"""
            if self.parent_step_ref is None:
                return self
            
            parent = session.query(WorkflowSteps).filter(WorkflowSteps.ref == self.parent_step_ref).first()
            if parent:
                return parent.get_root_step(session)
            return self

        def get_chain_filters(self, session):
            """Get all filters associated with this chain (from root step)"""
            root_step = self.get_root_step(session)
            return root_step.filters

        def get_full_chain(self, session):
            """Get the complete chain from root to this step"""
            chain = []
            current = self
            
            # Build chain backwards to root
            while current:
                chain.insert(0, current)
                if current.parent_step_ref:
                    current = session.query(WorkflowSteps).filter(WorkflowSteps.ref == current.parent_step_ref).first()
                else:
                    break
            
            return chain

        def get_chain_path(self, session):
            """Get the chain path as a string (e.g., 'mwcs -> mwcs_dtt')"""
            chain = self.get_full_chain(session)
            return ' -> '.join([step.step_type for step in chain])

        def get_expected_config_type(self):
            """Get the expected configuration type for this step"""
            return self.CONFIG_TYPE_MAP.get(self.step_type, self.step_type)

        def get_config_params(self, session):
            """Get all configuration parameters for this workflow step"""
            config_type = self.get_expected_config_type()
            return session.query(Config).filter(
                Config.category == config_type,
                Config.set_number == self.config_set_number,
                Config.used == True
            ).all()

        def get_config_dict(self, session):
            """Get configuration parameters as a dictionary with typed values"""
            config_type = self.get_expected_config_type()
            return Config.get_config_dict(session, config_type, self.config_set_number)

        def get_processing_steps(self, session):
            """Get all steps in order for processing this workflow chain"""
            chain = self.get_full_chain(session)
            return sorted(chain, key=lambda x: x.chain_order)

        def get_all_descendant_steps(self, session):
            """Get all descendant steps recursively"""
            descendants = []
            for child in self.child_steps:
                descendants.append(child)
                descendants.extend(child.get_all_descendant_steps(session))
            return descendants

        def get_terminal_steps(self, session):
            """Get all terminal steps reachable from this step"""
            if self.is_terminal_step():
                return [self]
            
            terminals = []
            for child in self.child_steps:
                if child.used:
                    terminals.extend(child.get_terminal_steps(session))
            return terminals

    ########################################################################

    class Filter(PrefixerBase):
        """
        Filter base class.

        :type ref: int
        :param ref: The id of the Filter in the database
        :type freqmin: float
        :param freqmin: The lower frequency bound of the Whiten function (in Hz)
        :type freqmax: float
        :param freqmax: The upper frequency bound of the Whiten function (in Hz)
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
        freqmin = Column(Float, nullable=False)
        freqmax = Column(Float, nullable=False)
        CC = Column(Boolean, default=True)
        SC = Column(Boolean, default=True)
        AC = Column(Boolean, default=True)
        used = Column(Boolean, default=True)

        # Many-to-Many relationship with WorkflowSteps (only root steps)
        workflow_steps = relationship("WorkflowSteps", secondary=filter_workflow_assoc, back_populates="filters")

        __table_args__ = (
            Index('idx_filter_freqs', 'freqmin', 'freqmax'),
            Index('idx_filter_used', 'used'),
        )

        def __init__(self, freqmin=None, freqmax=None, CC=True, SC=True, AC=True, used=True, **kwargs):
            """Initialize filter"""
            self.freqmin = freqmin
            self.freqmax = freqmax
            self.CC = CC
            self.SC = SC
            self.AC = AC
            self.used = used
            for key, val in kwargs.items():
                setattr(self, key, val)

        def __str__(self):
            return f"Filter {self.ref} ({self.freqmin}-{self.freqmax} Hz)"

        def validate_frequencies(self):
            """Validate frequency bounds"""
            if self.freqmin >= self.freqmax:
                raise ValueError(f"Minimum frequency ({self.freqmin}) must be less than "
                               f"maximum frequency ({self.freqmax})")
            if self.freqmin < 0:
                raise ValueError("Minimum frequency cannot be negative")

        def get_workflow_chains(self, session):
            """Get all complete workflow chains that use this filter"""
            chains = []
            for step in self.workflow_steps:
                if step.parent_step_ref is None:  # Root step
                    # Get all terminal steps in chains starting from this root
                    terminal_steps = self._get_terminal_steps(step, session)
                    chains.extend(terminal_steps)
            return chains

        def _get_terminal_steps(self, step, session):
            """Recursively find all terminal steps in the workflow chain"""
            if step.is_terminal_step():
                return [step]
            
            terminal_steps = []
            for child in step.child_steps:
                if child.used:
                    terminal_steps.extend(self._get_terminal_steps(child, session))
            
            return terminal_steps

        def get_processing_workflows(self, session):
            """Get all processing workflows (complete chains) for this filter"""
            workflows = []
            for root_step in self.workflow_steps:
                if root_step.parent_step_ref is None and root_step.used:
                    workflows.extend(self._build_workflows(root_step, session))
            return workflows

        def _build_workflows(self, step, session, current_chain=None):
            """Build complete workflow chains from a root step"""
            if current_chain is None:
                current_chain = []
            
            current_chain = current_chain + [step]
            
            if step.is_terminal_step():
                return [current_chain]
            
            workflows = []
            for child in step.child_steps:
                if child.used:
                    workflows.extend(self._build_workflows(child, session, current_chain))
            
            return workflows

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

        __table_args__ = (
            Index('job_index', "day", "pair", "jobtype", unique=True),
            Index('job_index2', "jobtype", "flag", unique=False),
        )

        def __init__(self, day=None, pair=None, jobtype=None, flag=None,
                     lastmod=None, **kwargs):
            """Initialize job"""
            self.day = day
            self.pair = pair
            self.jobtype = jobtype
            self.flag = flag
            self.lastmod = lastmod or datetime.datetime.utcnow()
            for key, val in kwargs.items():
                setattr(self, key, val)

        def __str__(self):
            return f"Job({self.day}, {self.pair}, {self.jobtype}, {self.flag})"

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
        X = Column(REAL)
        Y = Column(REAL)
        altitude = Column(Float)
        coordinates = Column(Enum('DEG', 'UTM', name="coordinate_type"))
        used = Column(Boolean, default=True)

        __table_args__ = (
            UniqueConstraint('net', 'sta', name='_station_unique'),
            Index('idx_station_net', 'net'),
            Index('idx_station_used', 'used'),
        )

        def __init__(self, net=None, sta=None, X=None, Y=None, altitude=None, 
                     coordinates=None, used=True, **kwargs):
            """Initialize station"""
            self.net = net
            self.sta = sta
            self.X = X
            self.Y = Y
            self.altitude = altitude
            self.coordinates = coordinates
            self.used = used
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
        path = Column(String(255))
        file = Column(String(255))
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
        )

        def __init__(self, net=None, sta=None, loc=None, chan=None, path=None, file=None, 
                     starttime=None, endtime=None, data_duration=None, gaps_duration=None, 
                     samplerate=None, flag=None, **kwargs):
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
            for key, val in kwargs.items():
                setattr(self, key, val)

        def __str__(self):
            return f"DataAvailability({self.net}.{self.sta}.{self.loc}.{self.chan})"

    ########################################################################

    # Return the schema namedtuple
    return sqlschema(Base, PrefixerBase,
                     Filter, Job, Station, Config, DataAvailability, WorkflowSteps,
                     filter_workflow_assoc)
    # end of declare_tables()

Base, PrefixerBase, Filter, Job, Station, Config, DataAvailability, WorkflowSteps, filter_workflow_assoc = declare_tables()


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
    param.updated_at = datetime.datetime.utcnow()
    param.validate_value()
    session.commit()
    return param


def get_workflow_config(session, step_type, config_set_number):
    """Get workflow configuration dictionary"""
    schema = declare_tables()
    return schema.Config.get_config_dict(session, step_type, config_set_number)


def validate_workflow_chain(session, steps):
    """Validate a complete workflow chain"""
    for i, step in enumerate(steps):
        step.validate_chain(session)
        if i > 0:
            if not steps[i-1].can_have_child(step.step_type):
                raise ValueError(f"Invalid chain: {steps[i-1].step_type} cannot be "
                               f"followed by {step.step_type}")