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

from logbook import Logger, StreamHandler
import sys

import numpy as np
import pandas as pd
import xarray as xr


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
    from sqlalchemy import create_engine
    from sqlalchemy.pool import NullPool
    dbini = read_db_inifile(inifile)
    if dbini.tech == 1:
        engine = create_engine('sqlite:///%s' % dbini.hostname, echo=False,
                               connect_args={'check_same_thread': False})
    elif dbini.tech == 2:
        engine = create_engine('mysql+pymysql://%s:%s@%s/%s'
                               % (dbini.username, dbini.password,
                                  dbini.hostname, dbini.database),
                               echo=False, poolclass=NullPool,
                               connect_args={'connect_timeout': 15})
    elif dbini.tech == 3:
        engine = create_engine('postgresql+psycopg2://%s:%s@%s/%s'
                               % (dbini.username, dbini.password,
                                  dbini.hostname, dbini.database),
                               echo=False, poolclass=NullPool,
                               connect_args={'connect_timeout': 15})
    else:
        raise ValueError("tech value must be 1, 2 or 3")
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
    from sqlalchemy.orm import sessionmaker
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

    :rtype: collections.namedtuple
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
        import pkg_resources
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
        import pkg_resources
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

def create_config_set(session, set_name):
    """
    Create a configuration set for a workflow step.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type set_name: str
    :param set_name: The name of the workflow step (e.g., 'mwcs', 'mwcs_dtt', etc.)

    :rtype: int or None
    :returns: The set_number if set was created successfully, None otherwise
    """
    import os
    import csv
    from sqlalchemy import func

    # Define the config file path
    config_file = os.path.join(os.path.dirname(__file__), 'config', f'config_{set_name}.csv')

    if not os.path.exists(config_file):
        return None

    # Find the next available set_number for this category
    max_set_number = session.query(func.max(Config.set_number)).filter(
        Config.category == set_name
    ).scalar()
    next_set_number = (max_set_number + 1) if max_set_number is not None else 1

    # Read the CSV file and add config entries
    with open(config_file, 'r', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            name = row['name']
            default_value = row['default']
            definition = row.get('definition', '')
            param_type = row.get('type', 'str')
            possible_values = row.get('possible_values', '')

            # Create new config entry with set_number
            config = Config(
                name=name,
                category=set_name,
                set_number=next_set_number,
                value=default_value,
                param_type=param_type,
                default_value=default_value,
                description=definition,
                possible_values=possible_values,
                used_in=f"[{set_name}]",
                used=True
            )
            session.add(config)

    session.commit()
    return next_set_number

def delete_config_set(session, set_name, set_number):
    """
    Delete a configuration set for a workflow step.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type set_name: str
    :param set_name: The name of the workflow step (e.g., 'mwcs', 'mwcs_dtt', etc.)
    :type set_number: int
    :param set_number: The set number to delete
    :rtype: bool
    :returns: True if set was deleted successfully, False otherwise
    """
    try:
        from sqlalchemy import func
        from .msnoise_table_def import Config

        # Validate input parameters
        if not isinstance(set_number, int) or set_number < 1:
            raise ValueError("set_number must be a positive integer")

        if not set_name or not isinstance(set_name, str):
            raise ValueError("set_name must be a non-empty string")

        # Check if the config set exists
        config_count = session.query(func.count(Config.ref)).filter(
            Config.category == set_name,
            Config.set_number == set_number
        ).scalar()

        if config_count == 0:
            return False  # Set doesn't exist

        # Delete all config entries for this set
        deleted_count = session.query(Config).filter(
            Config.category == set_name,
            Config.set_number == set_number
        ).delete()

        session.commit()

        get_logger("msnoise", "INFO").info(f"Deleted config set '{set_name}' set_number {set_number} "
                          f"({deleted_count} entries)")

        return True

    except Exception as e:
        session.rollback()
        traceback.print_exc()
        get_logger("msnoise").error(f"Failed to delete config set '{set_name}' set_number {set_number}: {str(e)}")
        return False


def list_config_sets(session, set_name=None):
    """
    List all configuration sets, optionally filtered by category.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type set_name: str or None
    :param set_name: Optional category filter (e.g., 'mwcs', 'mwcs_dtt', etc.)
    :rtype: list
    :returns: List of tuples (category, set_number, entry_count)
    """
    try:
        from sqlalchemy import func
        from .msnoise_table_def import Config

        query = session.query(
            Config.category,
            Config.set_number,
            func.count(Config.ref).label('entry_count')
        ).filter(
            Config.set_number.isnot(None)  # Exclude global configs
        )

        if set_name:
            query = query.filter(Config.category == set_name)

        results = query.group_by(
            Config.category,
            Config.set_number
        ).order_by(
            Config.category,
            Config.set_number
        ).all()

        return [(row.category, row.set_number, row.entry_count) for row in results]

    except Exception as e:
        traceback.print_exc()
        get_logger("msnoise").error(f"Failed to list config sets: {str(e)}")

        return []


def get_config_set_details(session, set_name, set_number):
    """
    Get details of a specific configuration set.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type set_name: str
    :param set_name: The category name
    :type set_number: int
    :param set_number: The set number
    :rtype: list
    :returns: List of config entries in the set
    """
    try:
        from .msnoise_table_def import Config

        configs = session.query(Config).filter(
            Config.category == set_name,
            Config.set_number == set_number
        ).order_by(Config.name).all()

        return [{'name': c.name, 'value': c.value, 'description': c.description}
                for c in configs]

    except Exception as e:
        get_logger("msnoise").error(f"Failed to get config set details: {str(e)}")
        return []

def get_config_categories_definition():
    """Get the standard configuration categories with their display names and order"""
    return [
        ('global', 'Global Parameters'),
        ('preprocess', 'Preprocessing'),
        ('cc', 'Cross-Correlation'),
        ('filter', 'Filters'),
        ('stack', 'Stacks'),
        ('mwcs', 'MWCS'),
        ('mwcs_dtt', 'MWCS dt/t'),
        ('stretching', 'Stretching'),
        ('wavelet', 'Wavelet'),
        ('wavelet_dtt', 'Wavelet dt/t'),
        ('qc', 'QC')
    ]

def get_config_sets_organized(session):
    """Get configuration sets organized by category in the standard order"""
    from sqlalchemy import func
    from . import msnoise_table_def as schema

    # Get category definitions
    category_order = get_config_categories_definition()

    # Get all config sets grouped by category and set_number
    all_sets = session.query(
        schema.Config.category,
        schema.Config.set_number,
        func.count(schema.Config.ref).label('param_count')
    ).group_by(
        schema.Config.category,
        schema.Config.set_number
    ).all()

    # Group sets by category
    sets_by_category = {}
    for category, set_number, param_count in all_sets:
        if category not in sets_by_category:
            sets_by_category[category] = []
        sets_by_category[category].append({
            'category': category,
            'set_number': set_number,
            'param_count': param_count
        })

    # Sort sets within each category by set_number
    for category in sets_by_category:
        sets_by_category[category].sort(key=lambda x: x['set_number'])

    # Create ordered list of categories with their sets
    ordered_categories = []
    for category_key, display_name in category_order:
        if category_key in sets_by_category:
            ordered_categories.append({
                'category': category_key,
                'display_name': display_name,
                'sets': sets_by_category[category_key]
            })

    # Add any additional categories not in the predefined order
    for category in sets_by_category:
        if category not in [cat[0] for cat in category_order]:
            ordered_categories.append({
                'category': category,
                'display_name': category.title(),
                'sets': sets_by_category[category]
            })

    return ordered_categories

def get_params(session):
    """Get config parameters from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: :class:`obspy.core.util.attribdict.AttribDict`
    :returns: a Param object containing the parameters
    """
    # TODO: this could be populated automatically from defauts iff defaults
    # would mention types
    from obspy.core.util.attribdict import AttribDict
    from .default import default
    s = session
    params = AttribDict()
    for name in default.keys():
        itemtype = default[name].type
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

    if not isinstance(params.mov_stack[0], tuple):
        params.mov_stack = [params.mov_stack, ]
    else:
        params.mov_stack = params.mov_stack
    return params

# FILTERS PART


def get_filters(session, all=False, ref=None):
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
    if ref:
        filter = session.query(Filter).filter(Filter.ref == ref).first()
        return filter
    if all:
        filters = session.query(Filter).all()
    else:
        filters = session.query(Filter).filter(Filter.used == True).all()
    return filters


def update_filter(session, ref, freqmin, freqmax,
                  CC, SC, AC, used):
    """Updates or Insert a new Filter in the database.

    .. seealso:: :class:`msnoise.msnoise_table_def.declare_tables.Filter`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

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
    :param used: Is the filter activated for the processing
    """
    filter = session.query(Filter).filter(Filter.ref == ref).first()
    if filter is None:
        filter = Filter()
        filter.freqmin = freqmin
        filter.freqmax = freqmax
        filter.CC = CC
        filter.SC = SC
        filter.AC = AC
        filter.used = used
        session.add(filter)
    else:
        filter.freqmin = freqmin
        filter.freqmax = freqmax
        filter.CC = CC
        filter.SC = SC
        filter.AC = AC
        filter.used = used
    session.commit()
    return

# DVV_MWCS

def get_dvv_mwcs_jobs(session, all=False):
    """Get all MWCS parameter sets linked to used filters."""

    query = (
        session.query(Filter.ref, DvvMwcs)
        .select_from(Filter)  # Explicitly select Filter as the starting table
        .join(filter_mwcs_assoc, Filter.ref == filter_mwcs_assoc.c.filt_ref)
        .join(DvvMwcs, filter_mwcs_assoc.c.dvv_mwcs_ref == DvvMwcs.ref)
    )

    if not all:
        query = query.filter(Filter.used == True)

    try:
        results = query.all()
    except Exception as e:
        print(f"ERROR: Query failed with exception: {e}")
        return {}

    mwcs_mapping = {}

    for filter_ref, mwcs in results:
        if filter_ref not in mwcs_mapping:
            mwcs_mapping[filter_ref] = []
        
        mwcs_mapping[filter_ref].append(mwcs)

    return mwcs_mapping

def get_dvv_mwcs(session, all=False):
    """Retrieve all MWCS parameter sets as ORM objects."""
    query = session.query(DvvMwcs)

    if not all:
        query = query.filter(DvvMwcs.used == True)

    results = query.all()  # Returns ORM objects

    return results  # No conversion to dictionaries

def update_dvv_mwcs(session, ref, freqmin, freqmax, mwcs_wlen, mwcs_step, used, filter_refs):
    """Updates or Inserts a new MWCS parameter set in the database and links it to filters.

    .. seealso:: :class:`msnoise.msnoise_table_def.declare_tables.DvvMwcs`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as obtained by :func:`connect`

    :type ref: int
    :param ref: The id of the MWCS parameter set in the database

    :type freqmin: float
    :param freqmin: The lower frequency bound of the MWCS step (in Hz)

    :type freqmax: float
    :param freqmax: The upper frequency bound of the MWCS step (in Hz)

    :type mwcs_wlen: float
    :param mwcs_wlen: Window length (in seconds) to perform MWCS

    :type mwcs_step: float
    :param mwcs_step: Step (in seconds) of the windowing procedure in MWCS

    :type used: bool
    :param used: Is the MWCS parameter set activated for processing

    :type filter_refs: list[int]
    :param filter_refs: A list of filter IDs to associate with this MWCS parameter set
    """

    # Debugging: Check the ref value before querying
    print(f"DEBUG: Received ref = {ref}")

    # Check if MWCS parameters already exist
    mwcs = session.query(DvvMwcs).filter(DvvMwcs.ref == ref).first()

    if mwcs is None:
        # Create a new MWCS parameter set
        print("DEBUG: Creating new MWCS parameter set")
        mwcs = DvvMwcs()
        mwcs.freqmin = freqmin
        mwcs.freqmax = freqmax
        mwcs.mwcs_wlen = mwcs_wlen
        mwcs.mwcs_step = mwcs_step
        mwcs.used = used
        session.add(mwcs)
        session.commit()  # Ensure the ref is assigned

        print(f"DEBUG: Assigned ref (after commit) = {mwcs.ref}")
    else:
        # Update existing MWCS parameter set
        mwcs.freqmin = freqmin
        mwcs.freqmax = freqmax
        mwcs.mwcs_wlen = mwcs_wlen
        mwcs.mwcs_step = mwcs_step
        mwcs.used = used

    # Print the MWCS object before linking filters
    print(f"DEBUG: MWCS Object: {mwcs}")

    # Update MWCS-Filter associations
    mwcs.filters = session.query(Filter).filter(Filter.ref.in_(filter_refs)).all()

    session.commit()
    print(f"DEBUG: Successfully updated MWCS with ref {mwcs.ref}")

    return

# DVV_MWCS_DTT

def get_dvv_mwcs_dtt_jobs(session, all=False):
    """Get all MWCS DTT parameter sets linked to used filters and MWCS sets.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A SQLAlchemy session object

    :type all: bool
    :param all: If True, fetch all filters; otherwise, only where `used` = True

    :rtype: dict
    :returns: A nested dictionary where:
              - `filter_ref`: The filter ID
              - `mwcs_ref`: The associated MWCS parameter ID
              - `dtt_params`: A list of associated MWCS DTT parameter sets
    """

    query = (
        session.query(Filter.ref, DvvMwcs.ref, DvvMwcsDtt)
        .join(filter_mwcs_assoc, Filter.ref == filter_mwcs_assoc.c.filt_ref)  # Explicit join
        .join(DvvMwcs, filter_mwcs_assoc.c.dvv_mwcs_ref == DvvMwcs.ref)
        .join(mwcs_dtt_assoc, DvvMwcs.ref == mwcs_dtt_assoc.c.dvv_mwcs_ref)
        .join(DvvMwcsDtt, mwcs_dtt_assoc.c.dvv_mwcs_dtt_ref == DvvMwcsDtt.ref)
    )

    if not all:
        query = query.filter(Filter.used == True)

    results = query.all()

    mwcs_dtt_mapping = {}

    for filter_ref, mwcs_ref, dtt in results:
        if filter_ref not in mwcs_dtt_mapping:
            mwcs_dtt_mapping[filter_ref] = {}

        if mwcs_ref not in mwcs_dtt_mapping[filter_ref]:
            mwcs_dtt_mapping[filter_ref][mwcs_ref] = []

        mwcs_dtt_mapping[filter_ref][mwcs_ref].append(dtt)  # Directly append the ORM object)

    return mwcs_dtt_mapping

def get_dvv_mwcs_dtt(session, all=False):
    """Retrieve all DTT parameter sets as ORM objects."""
    query = session.query(DvvMwcsDtt)

    if not all:
        query = query.filter(DvvMwcsDtt.used == True)

    results = query.all()  # Returns ORM objects

    return results  # No conversion to dictionaries

def update_dvv_mwcs_dtt(session, ref, dtt_minlag, dtt_width, dtt_lag, dtt_v,
                         dtt_sides, dtt_mincoh, dtt_maxerr, dtt_maxdtt, used,
                         mwcs_refs):
    """Updates or Inserts a new MWCS DTT parameter set in the database and links it to MWCS sets.

    .. seealso:: :class:`msnoise.msnoise_table_def.declare_tables.DvvMwcsDtt`

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as obtained by :func:`connect`

    :type ref: int
    :param ref: The id of the MWCS DTT parameter set in the database

    :type dtt_minlag: float
    :param dtt_minlag: Minimum lag time for MWCS DTT (in seconds)

    :type dtt_width: float
    :param dtt_width: Width of the time lag window (in seconds)

    :type dtt_lag: str
    :param dtt_lag: Method to define lag window [static/dynamic]

    :type dtt_v: float
    :param dtt_v: Velocity threshold for avoiding ballistic waves (if `dtt_lag=dynamic`)

    :type dtt_sides: str
    :param dtt_sides: Sides to use for DTT analysis [both/left/right]

    :type dtt_mincoh: float
    :param dtt_mincoh: Minimum coherence on dt measurement [0:1]

    :type dtt_maxerr: float
    :param dtt_maxerr: Maximum allowable error on dt measurement [0:1]

    :type dtt_maxdtt: float
    :param dtt_maxdtt: Maximum dt/t values for MWCS DTT (in seconds)

    :type used: bool
    :param used: Is the MWCS DTT parameter set activated for processing

    :type mwcs_refs: list[int]
    :param mwcs_refs: A list of MWCS parameter set IDs to associate with this DTT parameter set
    """

    # Check if MWCS DTT parameters already exist
    dtt = session.query(DvvMwcsDtt).filter(DvvMwcsDtt.ref == ref).first()

    if dtt is None:
        # Create a new MWCS DTT parameter set
        dtt = DvvMwcsDtt()
        dtt.dtt_minlag = dtt_minlag
        dtt.dtt_width = dtt_width
        dtt.dtt_lag = dtt_lag
        dtt.dtt_v = dtt_v
        dtt.dtt_sides = dtt_sides
        dtt.dtt_mincoh = dtt_mincoh
        dtt.dtt_maxerr = dtt_maxerr
        dtt.dtt_maxdtt = dtt_maxdtt
        dtt.used = used
        session.add(dtt)
    else:
        # Update existing MWCS DTT parameter set
        dtt.dtt_minlag = dtt_minlag
        dtt.dtt_width = dtt_width
        dtt.dtt_lag = dtt_lag
        dtt.dtt_v = dtt_v
        dtt.dtt_sides = dtt_sides
        dtt.dtt_mincoh = dtt_mincoh
        dtt.dtt_maxerr = dtt_maxerr
        dtt.dtt_maxdtt = dtt_maxdtt
        dtt.used = used

    # Update MWCS-DTT associations
    dtt.mwcs_params = session.query(DvvMwcs).filter(DvvMwcs.ref.in_(mwcs_refs)).all()

    session.commit()
    return

# DVV STRETCHING

def get_dvv_stretching_jobs(session, all=False):
    """Get all stretching parameter sets linked to used filters."""

    query = (
        session.query(Filter.ref, DvvStretching)
        .select_from(Filter)
        .join(filter_stretching_assoc, Filter.ref == filter_stretching_assoc.c.filt_ref)
        .join(DvvStretching, filter_stretching_assoc.c.dvv_stretching_ref == DvvStretching.ref)
    )

    if not all:
        query = query.filter(Filter.used == True)
    try:
        results = query.all()
    except Exception as e:
        print(f"ERROR: Query failed with exception: {e}")
        return {}

    stretching_mapping = {}

    for filter_ref, stretching in results:
        if filter_ref not in stretching_mapping:
            stretching_mapping[filter_ref] = []
        
        stretching_mapping[filter_ref].append(stretching)
    return stretching_mapping


def get_dvv_stretching(session, all=False):
    """Retrieve all Stretching parameter sets as ORM objects."""
    query = session.query(DvvStretching)

    if not all:
        query = query.filter(DvvStretching.used == True)

    results = query.all()  # Returns ORM objects

    return results 


def update_dvv_stretching(session, ref, stretching_minlag, stretching_width, stretching_lag,
                           stretching_v, stretching_sides, stretching_max, stretching_nsteps, used, filter_refs):
    """Updates or Inserts a new Stretching parameter set in the database and links it to filters."""

    stretching = session.query(DvvStretching).filter(DvvStretching.ref == ref).first()

    if stretching is None:
        stretching = DvvStretching()
        stretching.stretching_minlag = stretching_minlag
        stretching.stretching_width = stretching_width
        stretching.stretching_lag = stretching_lag
        stretching.stretching_v = stretching_v
        stretching.stretching_sides = stretching_sides
        stretching.stretching_max = stretching_max
        stretching.stretching_nsteps = stretching_nsteps
        stretching.used = used
        session.add(stretching)
        session.commit()
    else:
        stretching.stretching_minlag = stretching_minlag
        stretching.stretching_width = stretching_width
        stretching.stretching_lag = stretching_lag
        stretching.stretching_v = stretching_v
        stretching.stretching_sides = stretching_sides
        stretching.stretching_max = stretching_max
        stretching.stretching_nsteps = stretching_nsteps
        stretching.used = used

    stretching.filters = session.query(Filter).filter(Filter.ref.in_(filter_refs)).all()
    session.commit()

    return

# DVV WCT

def get_dvv_wct_jobs(session, all=False):
    """Get all WCT parameter sets linked to used filters."""

    query = (
        session.query(Filter.ref, DvvWct)
        .select_from(Filter)
        .join(filter_wct_assoc, Filter.ref == filter_wct_assoc.c.filt_ref)
        .join(DvvWct, filter_wct_assoc.c.dvv_wct_ref == DvvWct.ref)
    )

    if not all:
        query = query.filter(Filter.used == True)

    try:
        results = query.all()
    except Exception as e:
        print(f"ERROR: Query failed with exception: {e}")
        return {}

    wct_mapping = {}

    for filter_ref, wct in results:
        if filter_ref not in wct_mapping:
            wct_mapping[filter_ref] = []
        
        wct_mapping[filter_ref].append(wct)
    return wct_mapping


def get_dvv_wct(session, all=False):
    """Retrieve all WCT parameter sets as ORM objects."""
    query = session.query(DvvWct)

    if not all:
        query = query.filter(DvvWct.used == True)

    results = query.all()

    return results


def update_dvv_wct(session, ref, wct_freqmin, wct_freqmax, wct_ns, wct_nt, wct_vpo, wct_nptsfreq,
                    wct_norm, wavelet_type, used, filter_refs):
    """Updates or Inserts a new WCT parameter set in the database and links it to filters."""

    wct = session.query(DvvWct).filter(DvvWct.ref == ref).first()

    if wct is None:
        wct = DvvWct()
        wct.wct_freqmin = wct_freqmin
        wct.wct_freqmax = wct_freqmax
        wct.wct_ns = wct_ns
        wct.wct_nt = wct_nt
        wct.wct_vpo = wct_vpo
        wct.wct_nptsfreq = wct_nptsfreq
        wct.wct_norm = wct_norm
        wct.wavelet_type = wavelet_type
        wct.used = used
        session.add(wct)
        session.commit()
    else:
        wct.wct_freqmin = wct_freqmin
        wct.wct_freqmax = wct_freqmax
        wct.wct_ns = wct_ns
        wct.wct_nt = wct_nt
        wct.wct_vpo = wct_vpo
        wct.wct_nptsfreq = wct_nptsfreq
        wct.wct_norm = wct_norm
        wct.wavelet_type = wavelet_type
        wct.used = used

    wct.filters = session.query(Filter).filter(Filter.ref.in_(filter_refs)).all()
    session.commit()

    return

def get_dvv_wct_dtt_jobs(session, all=False):
    """Get all WCT DTT parameter sets linked to used WCT parameter sets."""
    
    query = (
        session.query(Filter.ref, DvvWct.ref, DvvWctDtt)
        .join(filter_wct_assoc, Filter.ref == filter_wct_assoc.c.filt_ref)  # Explicitly join filters
        .join(DvvWct, filter_wct_assoc.c.dvv_wct_ref == DvvWct.ref)  # Join WCT parameters
        .join(wct_dtt_assoc, DvvWct.ref == wct_dtt_assoc.c.dvv_wct_ref)  # Join WCT-DTT association
        .join(DvvWctDtt, wct_dtt_assoc.c.dvv_wct_dtt_ref == DvvWctDtt.ref)  # Join WCT DTT parameters
    )

    if not all:
        query = query.filter(DvvWct.used == True)

    results = query.all()

    wct_dtt_mapping = {}

    for filter_ref, wct_ref, dtt in results:
        if filter_ref not in wct_dtt_mapping:
            wct_dtt_mapping[filter_ref] = {}

        if wct_ref not in wct_dtt_mapping[filter_ref]:
            wct_dtt_mapping[filter_ref][wct_ref] = []

        wct_dtt_mapping[filter_ref][wct_ref].append(dtt) 

    return wct_dtt_mapping

def update_dvv_wct_dtt(session, ref, wct_dtt_freqmin, wct_dtt_freqmax, wct_minlag, wct_width,
                       wct_lag, wct_v, wct_sides, wct_mincoh, wct_maxdt, wct_codacycles,
                       wct_min_nonzero, used, wct_refs):
    """
    Updates or Inserts a new WCT-DTT parameter set in the database and links it to WCT parameter sets.

    Parameters:
    - session: SQLAlchemy session object.
    - ref: int, Reference ID for the WCT-DTT parameter set.
    - wct_dtt_freqmin, wct_dtt_freqmax: float, Frequency bounds.
    - wct_minlag: float, Minimum lag time (seconds).
    - wct_width: float, Width of the time lag window (seconds).
    - wct_lag: str, Lag window definition ('static' or 'dynamic').
    - wct_v: float, Velocity parameter for dynamic lag calculation.
    - wct_sides: str, Which sides to use ('both', 'left', 'right').
    - wct_mincoh: float, Minimum coherence threshold.
    - wct_maxdt: float, Maximum time delay (seconds).
    - wct_codacycles: int, Number of cycles between lag_min and lag_max.
    - wct_min_nonzero: float, Minimum nonzero weighting percentage for regression.
    - used: bool, Whether the parameter set is activated.
    - wct_refs: list, References to `DvvWct` parameter sets.

    Returns:
    - None
    """

    # Check if the WCT-DTT parameter set already exists
    wct_dtt = session.query(DvvWctDtt).filter(DvvWctDtt.ref == ref).first()

    if wct_dtt is None:
        # If not found, create a new one
        wct_dtt = DvvWctDtt()
        wct_dtt.wct_dtt_freqmin = wct_dtt_freqmin
        wct_dtt.wct_dtt_freqmax = wct_dtt_freqmax
        wct_dtt.wct_minlag = wct_minlag
        wct_dtt.wct_width = wct_width
        wct_dtt.wct_lag = wct_lag
        wct_dtt.wct_v = wct_v
        wct_dtt.wct_sides = wct_sides
        wct_dtt.wct_mincoh = wct_mincoh
        wct_dtt.wct_maxdt = wct_maxdt
        wct_dtt.wct_codacycles = wct_codacycles
        wct_dtt.wct_min_nonzero = wct_min_nonzero
        wct_dtt.used = used
        session.add(wct_dtt)
        session.commit()
    else:
        # Update existing entry
        wct_dtt.wct_dtt_freqmin = wct_dtt_freqmin
        wct_dtt.wct_dtt_freqmax = wct_dtt_freqmax
        wct_dtt.wct_minlag = wct_minlag
        wct_dtt.wct_width = wct_width
        wct_dtt.wct_lag = wct_lag
        wct_dtt.wct_v = wct_v
        wct_dtt.wct_sides = wct_sides
        wct_dtt.wct_mincoh = wct_mincoh
        wct_dtt.wct_maxdt = wct_maxdt
        wct_dtt.wct_codacycles = wct_codacycles
        wct_dtt.wct_min_nonzero = wct_min_nonzero
        wct_dtt.used = used

    # Link WCT-DTT parameters to WCT parameter sets
    wct_dtt.wct_params = session.query(DvvWct).filter(DvvWct.ref.in_(wct_refs)).all()

    session.commit()

def get_dvv_wct_dtt(session, all=False):
    """Retrieve all DTT parameter sets as ORM objects."""
    query = session.query(DvvWctDtt)

    if not all:
        query = query.filter(DvvWctDtt.used == True)

    results = query.all() 

    return results  


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


def get_stations(session, all=False, net=None, format="raw"):
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
    if format == "raw":
        return stations
    if format == "seed_id":
        output = []
        for sta in stations:
            if sta.used_location_codes is None:
                location_codes = []
            else:
                location_codes = sta.used_location_codes.split(",")
            if sta.used_channel_names is None:
                channels = []
            else:
                channels = []
                for i, chan in enumerate(sta.used_channel_names.split(",")):
                    if chan.count("?"):
                        for comp in ["Z", "N", "E", "1", "2"]:
                            channels.append(chan.replace("?", comp))
                    else:
                        channels.append(chan)

            for loc in location_codes:
                for chan in channels:
                    output.append("%s.%s.%s.%s" % (sta.net, sta.sta, loc, chan))
        return output


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
    :param X: The X coordinate of the station (Easting or Longitude)
    :type Y: float
    :param Y: The Y coordinate of the station (Northing or Latitude)
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

    if coordinates == "DEG" and (not -90 <= Y <= 90 or not -180 <= X <= 180):
        raise ValueError("Coordinates must be valid WGS84 latitude (%.4f) and longitude (%.4f). " % (Y, X))

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
    if len(get_config(session, name="components_to_compute_single_station")):
        return itertools.combinations_with_replacement(stations, 2)
    else:
        return itertools.combinations(stations, 2)


def check_stations_uniqueness(session, station):
    """

    :param session:
    :param station:
    :return:
    """
    # if the station is net.sta.loc, nothing to do
    if station.count(".") == 2:
        return station

    logging.info("It seems you're voluntarily missing the location code for"
                 " \"%s\". We'll handle this automatically, if there are no "
                 "conflicts." % station)
    net, sta = station.split(".")
    locs = get_station(session, net, sta).locs()
    if len(locs) != 1:
        logging.info("There are more than 1 location codes for this station: "
                     "%s" % locs)
        return station
    station += ".%s" % locs[0]
    logging.info("Found %s to be the unique solution for this station" % station)
    return station



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
    from obspy.geodetics import gps2dist_azimuth
    if coordinates == "DEG":
        dist, azim, bazim = gps2dist_azimuth(station1.Y, station1.X,
                                             station2.Y, station2.X)
        return dist / 1.e3
    else:
        dist = np.hypot(float(station1.X - station2.X),
                        float(station1.Y - station2.Y)) / 1.e3
        return dist


# DATA AVAILABILITY


def update_data_availability(session, net, sta, loc, chan, path, file, starttime,
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
    :type chan: str
    :param chan: The component (channel)
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
        filter(DataAvailability.path == path). \
        filter(DataAvailability.file == file).\
        filter(DataAvailability.net == net).\
        filter(DataAvailability.sta == sta). \
        filter(DataAvailability.loc == loc). \
        filter(DataAvailability.chan == chan).first()
    if data is None:
        flag = "N"
        data = DataAvailability(net, sta, loc, chan, path, file, starttime, endtime,
                                data_duration, gaps_duration, samplerate, flag)
        session.add(data)
        toreturn = 1
    else:
        modified = False
        for item in ['net', 'sta', 'loc', 'chan', 'path', 'starttime',
                     'endtime', 'data_duration', 'gaps_duration', 'samplerate']:
            if eval("data.%s != %s" % (item, item)):
                modified = True
                break
        if modified:
            data.net = net
            data.sta = sta
            data.loc = loc
            data.chan = chan
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


def get_data_availability(session, net=None, sta=None, loc=None, chan=None,
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
    from sqlalchemy.sql.expression import func
    if not starttime:
        data = session.query(DataAvailability).\
            filter(DataAvailability.net == net).\
            filter(DataAvailability.sta == sta). \
            filter(DataAvailability.loc == loc). \
            filter(DataAvailability.chan == chan).all()
    elif not net:
        data = session.query(DataAvailability).\
            filter(DataAvailability.starttime <= endtime).\
            filter(DataAvailability.endtime >= starttime).all()
    else:
        data = session.query(DataAvailability).\
            filter(DataAvailability.net == net).\
            filter(DataAvailability.sta == sta). \
            filter(DataAvailability.loc == loc). \
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
    from sqlalchemy.sql.expression import func
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
    from sqlalchemy import text
    if ref:
        job = session.query(Job).filter(text("ref=:ref")).params(ref=ref).first()
    else:
        job = session.query(Job)\
            .filter(text("day=:day"))\
            .filter(text("pair=:pair"))\
            .filter(text("jobtype=:jobtype"))\
            .params(day=day, pair=pair, jobtype=jobtype).first()
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
    with engine.connect() as conn:
        conn.execute(
            Job.__table__.insert(),
            jobs)
        try:
            conn.commit()
        except:
            pass


def massive_update_job(session, jobs, flag="D"):
    """
    Routine to use a low level function to update much faster a list of
    :class:`~msnoise.msnoise_table_def.declare_tables.Job`. This method uses the Job.ref
    which is unique.
    :type session: Session
    :param session: the database connection object
    :type jobs: list or tuple
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
    job = session.query(Job).\
        filter(Job.jobtype == jobtype).\
        filter(Job.flag == flag).first()
    if job is None:
        return False
    else:
        return True


def get_next_job(session, flag='T', jobtype='CC', limit=99999):
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
    from sqlalchemy import update
    tmp = []
    while not len(tmp):
        jobs = session.query(Job).filter(Job.jobtype == jobtype).\
            filter(Job.flag == flag).\
            filter(Job.day == session.query(Job).
                   filter(Job.jobtype == jobtype).
                   filter(Job.flag == flag).first().day).\
            limit(limit).with_for_update()

        tmp = jobs.all()
        refs = [_.ref for _ in tmp]
        q = update(Job).values({"flag":"I"}).where(Job.ref.in_(refs))
        session.execute(q)
        # jobs.update({Job.flag: 'I'})
        session.commit()
    return tmp


def get_dvv_jobs(session, flag='T', jobtype='DVV', limit=99999):
    from sqlalchemy import update
    tmp = []
    while not len(tmp):
        jobs = session.query(Job).filter(Job.jobtype == jobtype). \
            filter(Job.flag == flag). \
            limit(limit).with_for_update()

        tmp = jobs.all()
        refs = [_.ref for _ in tmp]
        q = update(Job).values({"flag": "I"}).where(Job.ref.in_(refs))
        session.execute(q)
        # jobs.update({Job.flag: 'I'})
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
    q = session.query(Job.ref).\
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
    from sqlalchemy.sql.expression import func
    if read_db_inifile().tech == 2:
        rand = func.rand
    else:
        rand = func.random
    try:
        jobs = session.query(Job.ref, Job.day, Job.pair, Job.flag).filter(Job.flag == flag).\
            filter(Job.jobtype == jobtype).filter(Job.day != 'REF').\
            filter(Job.pair == session.query(Job).
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
    from sqlalchemy.sql.expression import func
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
        os.makedirs(path, exist_ok=True)

    df = pd.DataFrame().from_dict(data).T
    df.columns = get_t_axis(session)
    df.to_hdf(os.path.join(path, date+'.h5'), key='data')
    del df
    return


def export_allcorr2(session, ccfid, data):
    output_folder = get_config(session, 'output_folder')
    station1, station2, components, filterid, date = ccfid.split('_')

    path = os.path.join(output_folder, "%02i" % int(filterid),
                        station1, station2, components)
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

    df = pd.DataFrame().from_dict(data).T
    df.columns = get_t_axis(session)
    df.to_hdf(os.path.join(path, date+'.h5'), key='data')
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
    :param station1: The name of station 1 (formatted NET.STA.LOC)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA.LOC)
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
    :type params: dict, :class:`obspy.core.util.attribdict.AttribDict`
    :param params: A dictionnary of MSNoise config parameters as returned by
        :func:`get_params`.
    """
    from obspy import Stream, Trace
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
            os.makedirs(path, exist_ok=True)

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
    from obspy.core.util.attribdict import AttribDict
    from obspy import Stream, Trace
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
        os.makedirs(os.path.split(filename)[0], exist_ok=True)
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
    from obspy import Trace, Stream
    try:
        os.makedirs(os.path.split(filename)[0], exist_ok=True)
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
        import scipy.signal as ss
        # logging.debug("Doing a PWS stack")
        corr = np.zeros(data.shape[1], dtype='f8')
        phasestack = np.zeros(data.shape[1], dtype='c8')
        for i in range(data.shape[0]):
            data[i] -= data[i].mean()
        for c in data:
            phase = np.angle(ss.hilbert(c))
            phasestack.real += np.cos(phase)
            phasestack.imag += np.sin(phase)
        coh = 1. / data.shape[0] * np.abs(phasestack)

        timegate_samples = int(pws_timegate * goal_sampling_rate)
        coh = np.convolve(ss.windows.boxcar(timegate_samples) /
                          timegate_samples, coh, 'same')
        coh = np.power(coh, pws_power)
        for c in data:
            corr += c * coh
        corr /= data.shape[0]

    return corr


def get_extension(export_format):
    if export_format == "BOTH":
        return ".MSEED"
    elif export_format == "SAC":
        return ".SAC"
    elif export_format == "MSEED":
        return ".MSEED"
    else:
        return ".MSEED"


def get_ref(session, station1, station2, filterid, components, params=None):
    """
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET.STA.LOC)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA.LOC)
    :type filterid: int
    :param filterid: The ID (ref) of the filter
    :type components: str
    :param components: The name of the components used (ZZ, ZR, ...)
    :type params: dict, :class:`obspy.core.util.attribdict.AttribDict`
    :param params: A dictionnary of MSNoise config parameters as returned by
        :func:`get_params`.
    :rtype: :class:`obspy.trace`
    :return: A Trace object containing the ref
    """
    from obspy import Trace, read
    if not params:
        export_format = get_config(session, 'export_format')
        extension = get_extension(export_format)
    else:
        extension = get_extension(params.export_format)

    ref_name = "%s_%s" % (station1, station2)
    ref_name = ref_name
    rf = os.path.join("STACKS", "%02i" %
                      filterid, "REF", components,
                      ref_name + extension)
    if not os.path.isfile(rf):
        logging.debug("No REF file named %s, skipping." % rf)
        return Trace()

    return read(rf)[0]


def get_results(session, station1, station2, filterid, components, dates,
                mov_stack=1, format="stack", params=None):
    """
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET.STA.LOC)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA.LOC)
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
    :type params: dict, :class:`obspy.core.util.attribdict.AttribDict`
    :param params: A dictionnary of MSNoise config parameters as returned by
        :func:`get_params`.
    :rtype: :class:`numpy.ndarray`
    :return: Either a 1D CCF (if format is ``stack`` or a 2D array (if format=
        ``matrix``).
    """
    from obspy import read
    if not params:
        export_format = get_config(session, 'export_format')
        extension = get_extension(export_format)
    else:
        export_format = params.export_format
        extension = get_extension(params.export_format)
    if export_format == "BOTH":
        export_format = "MSEED"

    stack_data = np.zeros((len(dates), get_maxlag_samples(session))) * np.nan
    i = 0
    base = os.path.join("STACKS", "%02i" % filterid,
                        "%03i_DAYS" % mov_stack, components,
                        "%s_%s" % (station1, station2), "%s") + extension
    logging.debug("Reading files... in %s" % base)
    lastday = dates[0]
    for j, date in enumerate(dates):

        if isinstance(date, str):
            daystack = base % str(date)
        else:
            daystack = base % date.strftime('%Y-%m-%d')

        try:
            stack_data[j, :] = read(daystack, format=export_format)[0].data[:]
            lastday = str(date)
            i += 1
        except:
            # traceback.print_exc()
            pass

    if format == "matrix":
        return i, stack_data

    elif format == "dataframe":
        taxis = get_t_axis(session)
        return pd.DataFrame(stack_data, index=pd.DatetimeIndex(dates),
                            columns=taxis).loc[:lastday]
    elif format == "xarray":
        taxis = get_t_axis(session)
        times = pd.DatetimeIndex(dates)
        dr = xr.DataArray(stack_data, coords=[times, taxis],
                          dims=["times", "taxis"]).dropna("times", how="all")
        dr.name = "CCF"
        return dr.to_dataset()

    elif format == "stack":
        logging.debug("Stacking...")

        corr = stack(stack_data, params.stack_method, params.pws_timegate,
                     params.pws_power, params.goal_sampling_rate)

        if i > 0:
            return i, corr
        else:
            return 0, None



def get_mwcs(session, station1, station2, filterid, components, date,
                mov_stack=1):
    """
    TODO
    """
    file = os.path.join('MWCS', "%02i" % filterid, "%s_%s" % (mov_stack[0], mov_stack[1]),
                        components, "%s_%s" % (station1, station2),
                        '%s.txt' % date)
    if os.path.isfile(file):
        df = pd.read_csv(
            file, delimiter=' ', header=None, index_col=0,
            names=['t', 'dt', 'err', 'coh'])
        return df
    else:
        return pd.DataFrame()


def get_results_all(session, station1, station2, filterid, components, dates,
                    format="dataframe"):
    """
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET.STA.LOC)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA.LOC)
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
        if isinstance(date, str):
            fname = os.path.join(path, date+".h5")
        else:
            fname = os.path.join(path, date.strftime('%Y-%m-%d.h5'))
        if os.path.isfile(fname):
            df = pd.read_hdf(fname, 'data', parse_dates=True)
            df.index = pd.to_datetime(df.index)
            results.append(df)
    if len(results):
        result = pd.concat(results)
        del results
        if format == "dataframe":
            return result
        elif format == "xarray":
            taxis = get_t_axis(session)
            times = result.index
            dr = xr.DataArray(result, coords=[times, taxis],
                              dims=["times", "taxis"]).dropna("times", how="all")
            dr.name = "CCF"
            dr = dr.sortby('times')
            return dr.to_dataset()
    else:
        if format == "xarray":
            return xr.Dataset()
        else:
            return pd.DataFrame()

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

def get_filter_components_to_compute(session, filterid, params):
    """
    Returns the components to compute for a specific filter based on whether CC, SC, or AC are enabled.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as obtained by :func:`connect`

    :type filterid: int
    :param filterid: The ID of the filter to get components for.

    :type params: obspy.core.util.attribdict.AttribDict
    :param params: A parameters object containing components_to_compute and components_to_compute_single_station.

    :rtype: tuple of lists
    :returns: (components_to_compute, components_to_compute_single_station)
    """

    # Get the filter object
    filter_obj = session.query(Filter).filter(Filter.ref == filterid).first()
    if not filter_obj:
        return [], []

    # Get the global component lists
    global_components = params.components_to_compute
    global_components_single_station = params.components_to_compute_single_station

    # Initialize empty lists for selected components
    filt_components = []
    filt_components_single_station = []

    # **Filter components based on filter settings**
    if filter_obj.CC:
        filt_components.extend(global_components)

    if filter_obj.SC or filter_obj.AC:
        for comp in global_components_single_station:
            if comp[0] == comp[1]:  # Auto-correlation (AC) → same component (e.g., ZZ, EE, NN)
                if filter_obj.AC:
                    filt_components_single_station.append(comp)
            else:  # Cross-component (SC) → different components (e.g., ZN, EZ, NE)
                if filter_obj.SC:
                    filt_components_single_station.append(comp)

    # Remove duplicates and sort for consistency
    filt_components = list(np.unique(filt_components))
    filt_components_single_station = list(np.unique(filt_components_single_station))

    return filt_components, filt_components_single_station

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

def validate_stack_data(dataset, stack_type="reference"):
    """Validates stack data before processing
    
    Parameters:
        dataset: xarray Dataset to validate 
        stack_type: Type of stack ("reference" or "moving") for error messages
    Returns:
        (is_valid, message) tuple
    """
    if dataset is None or not dataset.data_vars:
        return False, f"No data found for {stack_type} stack"
        
    if not hasattr(dataset, 'CCF'):
        return False, f"Missing CCF data in {stack_type} stack"
        
    data = dataset.CCF
    if data.size == 0:
        return False, f"Empty dataset in {stack_type} stack"
        
    nan_count = np.isnan(data.values).sum()
    total_points = data.values.size
    
    if nan_count == total_points:
        return False, f"{stack_type.capitalize()} stack contains only NaN values"
    
    if nan_count > 0:
        percent_nan = (nan_count / total_points) * 100
        return True, f"Warning: {stack_type.capitalize()} stack contains {percent_nan:.1f}% NaN values"
        
    return True, "OK"
    
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
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d').date()

    if end == "2100-01-01":
        end = session.query(DataAvailability).order_by(
            DataAvailability.endtime.desc()).first().endtime.date()
    else:
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
    from sqlalchemy import text
    lastmod = datetime.datetime.now() - interval
    if pair == '%':
        days = session.query(Job).\
            filter(text("day>=:date1")).\
            filter(text("day<=:date2")).\
            filter(Job.jobtype == jobtype).\
            filter(Job.lastmod >= lastmod).group_by(Job.day).\
            order_by(Job.day).params(date1=date1.strftime("%Y-%m-%d"),
                                     date2=date2.strftime("%Y-%m-%d")).all()
    else:
        days = session.query(Job).filter(Job.pair == pair).\
            filter(Job.day >= date1.strftime("%Y-%m-%d")). \
            filter(Job.day <= date2.strftime("%Y-%m-%d")). \
            filter(Job.jobtype == jobtype).filter(Job.lastmod >= lastmod).\
            group_by(Job.day).order_by(Job.day).with_entities("day").all()
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
    from obspy.geodetics import gps2dist_azimuth
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


def check_and_phase_shift(trace, taper_length=20.0):
    # TODO replace this hard coded taper length

    import scipy.fft as sf
    from scipy.fft import next_fast_len
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
        FFTdata = sf.fft(trace.data, n=n)
        fftfreq = sf.fftfreq(n, d=trace.stats.delta)
        FFTdata = FFTdata * np.exp(1j * 2. * np.pi * fftfreq * dt)
        FFTdata = FFTdata.astype(np.complex64)
        sf.ifft(FFTdata, n=n, overwrite_x=True)
        trace.data = np.real(FFTdata[:len(trace.data)]).astype(float)
        trace.stats.starttime += dt
        del FFTdata, fftfreq
        # clean_scipy_cache()
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
    from obspy import Stream
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

    if len(st) < 2:
        return st

    masks=[]
    for tr in st:
        masks.append(tr.data.mask)
    mask =  np.any(masks,axis=0)

    # apply the mask to all traces
    for tr in st:
        tr.data.mask = mask

    st = st.split()
    return st

def preload_instrument_responses(session, return_format="dataframe"):
    """
    This function preloads all instrument responses from ``response_path``
    and stores the seed ids, start and end dates, and paz for every channel
    in a DataFrame. Any file readable by obspy's read_inventory will be processed.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :type return_format: str
    :param return_format: The format of the returned object, either
        ``dataframe`` or ``inventory``.

    :rtype: :class:`~pandas.DataFrame` or :class:`~obspy.core.inventory.inventory.Inventory`
    :returns: A table containing all channels with the time of operation and
        poles and zeros (DataFrame), or an obspy Inventory object.

    """
    from obspy.core.inventory import Inventory
    from obspy import read_inventory, UTCDateTime
    logging.debug('Preloading instrument response')
    response_format = get_config(session, 'response_format')
    files = glob.glob(os.path.join(get_config(session, 'response_path'), "*"))
    channels = []
    all_inv = Inventory()
    for file in files:
        logging.debug("Processing %s" % file)
        try:
            inv = read_inventory(file)

            if return_format == "inventory":
                all_inv += inv
                continue

            for net in inv.networks:
                for sta in net.stations:
                    for cha in sta.channels:
                        seed_id = "%s.%s.%s.%s" % (net.code, sta.code,
                                                   cha.location_code,
                                                   cha.code)
                        pzdict = {}
                        try:
                            resp = inv.get_response(seed_id, cha.start_date + 10)
                            polezerostage = resp.get_paz()
                        except Exception as e:
                            logging.warning(
                                'Failed to get PAZ for SEED ID "%s", this '
                                'SEED ID will have an empty dictionary '
                                'for Poles and Zeros '
                                'information (Error message: %s).' % (
                                    seed_id, str(e)))
                        else:
                            totalsensitivity = resp.instrument_sensitivity
                            pzdict['poles'] = polezerostage.poles
                            pzdict['zeros'] = polezerostage.zeros
                            pzdict['gain'] = polezerostage.normalization_factor
                            pzdict['sensitivity'] = totalsensitivity.value
                        lat = cha.latitude
                        lon = cha.longitude
                        elevation = cha.elevation
                        if lat is None or lon is None or elevation is None:
                            lat = sta.latitude
                            lon = sta.longitude
                            elevation = sta.elevation
                        if lat is None or lon is None or elevation is None:
                            logging.error(
                                'Failed to look up coordinates for SEED '
                                'ID: %s' % seed_id)
                        channels.append([seed_id, cha.start_date,
                                         cha.end_date or UTCDateTime(),
                                         pzdict, lat, lon, elevation])

        except Exception as e:
            logging.error('Failed to process file %s: %s' % (file, str(e)))


    logging.debug('Finished Loading instrument responses')
    if return_format == "inventory":
        return all_inv

    if return_format == "dataframe":
        channels = pd.DataFrame(channels, columns=["channel_id", "start_date",
                                                   "end_date", "paz",
                                                   "latitude", "longitude", "elevation"],)
        return channels



# Needed for QC/PPSD, should re-order the defs


def to_sds(stats,year, jday):
    SDS="YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.JDAY"
    file=SDS.replace('YEAR', "%04i"%year)
    file=file.replace('NET', stats.network)
    file=file.replace('STA', stats.station)
    file=file.replace('LOC', stats.location)
    file=file.replace('CHAN', stats.channel)
    file=file.replace('JDAY', "%03i"%jday)
    file=file.replace('TYPE', "D")
    return file

## PSD part (not sure it'll end up here but easier to handle for now)
from functools import lru_cache


def psd_read_results(net, sta, loc, chan, datelist, format='PPSD', use_cache=True):
    from obspy.signal import PPSD
    if loc == "--":
        loc = ""
    fn = "%s.%s.%s.%s-%s_%s.npz" % (net, sta, loc, chan, datelist[0], datelist[-1])
    import tempfile
    fn = os.path.join(tempfile.gettempdir(), "MSNOISE-PSD", fn)
    if use_cache and os.path.isfile(fn):
        logging.debug("I found this cool file: %s" % fn)
        ppsd = PPSD.load_npz(fn)
    else:
        first = True
        ppsd = None
        for day in datelist:
            jday = int(day.strftime("%j"))
            toglob = os.path.join('PSD', 'NPZ', "%s" % day.year, net, sta,
                                  chan + ".D", "%s.%s.%s.%s.D.%s.%03i.npz" % (
                                  net, sta, loc, chan, day.year, jday))
            files = glob.glob(toglob)
            if not len(files):
                logging.error("No files found for %s.%s.%s.%s: %s" % (
                net, sta, loc, chan, day))
                continue
            file = files[0]
            if os.path.isfile(file):
                if first:
                    ppsd = PPSD.load_npz(file)
                    first = False
                else:
                    try:
                        ppsd.add_npz(file)
                    except:
                        pass
    if not ppsd:
        return None
    if use_cache:
        if not os.path.isdir(os.path.split(fn)[0]):
            os.makedirs(os.path.split(fn)[0], exist_ok=True)
        ppsd.save_npz(fn[:-4])
    return ppsd


def psd_ppsd_to_dataframe(ppsd):
    from obspy import UTCDateTime
    ind_times = np.array(
        [UTCDateTime(t).datetime for t in ppsd.current_times_used])
    data = np.asarray(ppsd._binned_psds)
    return pd.DataFrame(data, index=ind_times, columns=ppsd.period_bin_centers)


def hdf_open_store_from_fn(fn, mode="a"):
    store = pd.HDFStore(fn, complevel=9, complib="blosc:blosclz", mode=mode)
    return store


def hdf_open_store(filename, location=os.path.join("PSD", "HDF"), mode="a",
                   format='table'):
    if ".h5" in filename:
        filename = filename.replace(".h5", "")
    pd.set_option('io.hdf.default_format', format)
    if not os.path.isdir(location):
        os.makedirs(location, exist_ok=True)
    fn = os.path.join(location, filename + ".h5")
    store = pd.HDFStore(fn, complevel=9, complib="blosc:blosclz", mode=mode)
    return store


def hdf_insert_or_update(store, key, new):
    if key in store:
        filter = store[key].index.intersection(new.index)
        if len(filter):
            coordinates = store.select_as_coordinates(key, "index=filter")
            store.remove(key, where=coordinates)
            store.append(key, new, format='t', data_columns=True, append=True)
    else:
        store.append(key, new)


def hdf_close_store(store):
    store.close()
    del store


def xr_create_or_open(fn, taxis=[], name="CCF"):
    if os.path.isfile(fn):
        # load_dataset works (it loads content in mem and closes, open_dataset
        # failed, the file handle was still open and it failed later.
        ds = xr.load_dataset(fn)
        return ds
    times = pd.date_range("2000-01-01", freq="h", periods=0)
    if name == "CCF":
        data = np.random.random((len(times), len(taxis)))
        dr = xr.DataArray(data, coords=[times, taxis], dims=["times", "taxis"])
    elif name == "REF":
        data = np.random.random(len(taxis))
        dr = xr.DataArray(data, coords=[taxis], dims=["taxis"])
    elif name == "MWCS":
        keys = ["M", "EM", "MCOH"]
        data = np.random.random((len(times), len(taxis), len(keys)))
        dr = xr.DataArray(data, coords=[times, taxis, keys],
                          dims=["times", "taxis", "keys"])
    elif name == "DTT":
        keys = ["m", "em", "a", "ea", "m0", "em0"]
        data = np.random.random((len(times), len(keys)))
        dr = xr.DataArray(data, coords=[times, keys],
                          dims=["times", "keys"])
    elif name == "DVV":
        level0 = ["m", "em", "a", "ea", "m0", "em0"]
        level1 = ['10%', '25%', '5%', '50%', '75%', '90%', '95%', 'count', 'max', 'mean',
                  'min', 'std', 'trimmed_mean', 'trimmed_std', 'weighted_mean', 'weighted_std']
        data = np.random.random((len(times), len(level0), len(level1)))
        dr = xr.DataArray(data, coords=[times, level0, level1],
                          dims=["times", "level0", "level1"])
    elif name == "WCT":
        dvv_data = np.random.random((len(times), len(taxis)))
        err_data = np.random.random((len(times), len(taxis)))
        coh_data = np.random.random((len(times), len(taxis)))

        dvv_da = xr.DataArray(dvv_data, coords=[times, taxis], dims=['times', 'frequency'])
        err_da = xr.DataArray(err_data, coords=[times, taxis], dims=['times', 'frequency'])
        coh_da = xr.DataArray(coh_data, coords=[times, taxis], dims=['times', 'frequency'])

        # Combine into a Dataset
        ds = xr.Dataset({'dvv': dvv_da,'err': err_da,'coh': coh_da})
        return ds
        
    else:
        logging.error("Not implemented, name=%s invalid." % name)
        sys.exit(1)
    dr.name = name
    return dr.to_dataset()


def xr_insert_or_update(dataset, new):
    tt = new.merge(dataset, compat='override', combine_attrs="drop_conflicts")
    return tt.combine_first(dataset)


def xr_save_and_close(dataset, fn):
    if not os.path.isdir(os.path.split(fn)[0]):
        os.makedirs(os.path.split(fn)[0], exist_ok=True)
    dataset.to_netcdf(fn, mode="w")
    dataset.close()
    del dataset


def get_dvv(session, filterid, components, dates,
            mov_stack=1, aggregation="median", dttname="M"):
    pairs = []

    current = dates[0]
    end = dates[-1]

    alldf = []
    while current <= end:
        for comp in components:
            day = os.path.join('DTT', "%02i" % filterid, "%03i_DAYS" %
                               mov_stack, components, '%s.txt' % current)
            if os.path.isfile(day):
                df = pd.read_csv(day, header=0, index_col=0,
                                 parse_dates=True)
                alldf.append(df)
        current += datetime.timedelta(days=1)
    if len(alldf) == 0:
        print("No Data for %s m%i f%i" % (components, mov_stack, filterid))

    alldf = pd.concat(alldf)
    # print(mov_stack, alldf.head())
    if 'alldf' in locals():
        errname = "E" + dttname
        alldf.to_csv("tt.csv")
        alldf[dttname] *= -100
        alldf[errname] *= -100

        allbut = alldf[alldf['Pairs'] != 'ALL'].copy()

        for pair in pairs:
            print(pair)
            pair1 = alldf[alldf['Pairs'] == pair].copy()
            print(pair1.head())

            pair1.to_csv('%s-m%i-f%i.csv' % (pair, mov_stack, filterid))

        if aggregation == "median":
            tmp3 = allbut[dttname].resample('D').median()
            etmp3 = allbut[errname].resample('D').median()

        elif aggregation == "mean":
            tmp3 = allbut[dttname].resample('D').mean()
            etmp3 = allbut[errname].resample('D').mean()
        else:
            print('Choose median or mean')
    return tmp3, etmp3

def xr_save_ccf(station1, station2, components, filterid, mov_stack, taxis, new, overwrite=False):
    path = os.path.join("STACKS2", "%02i" % filterid,
                        "%s_%s" % (mov_stack[0], mov_stack[1]), "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)
    fullpath = os.path.join(path, fn)
    if overwrite:
        xr_save_and_close(new, fullpath)
    else:
        dr = xr_create_or_open(fullpath, taxis, name="CCF")
        dr = xr_insert_or_update(dr, new)
        dr = dr.sortby("times")
        xr_save_and_close(dr, fullpath)
        return dr


def xr_get_ccf(station1, station2, components, filterid, mov_stack, taxis, format="dataframe"):
    path = os.path.join("STACKS2", "%02i" % filterid,
                        "%s_%s" % (mov_stack[0], mov_stack[1]), "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)

    fullpath = os.path.join(path, fn)
    if not os.path.isfile(fullpath):
        # logging.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
        raise FileNotFoundError(fullpath)
    data = xr_create_or_open(fullpath, taxis, name="CCF")
    if format == "xarray":
        return data.CCF
    else:
        return data.CCF.to_dataframe().unstack().droplevel(0, axis=1)


def xr_save_ref(station1, station2, components, filterid, taxis, new, overwrite=False):
    path = os.path.join("STACKS2", "%02i" % filterid,
                        "REF", "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)
    fullpath = os.path.join(path, fn)
    if overwrite:
        xr_save_and_close(new, fullpath)
    else:
        dr = xr_create_or_open(fullpath, taxis, name="REF")
        dr = xr_insert_or_update(dr, new)
        xr_save_and_close(dr, fullpath)
        return dr


def xr_get_ref(station1, station2, components, filterid, taxis):
    path = os.path.join("STACKS2", "%02i" % filterid,
                        "REF", "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)

    fullpath = os.path.join(path, fn)
    if not os.path.isfile(fullpath):
        # logging.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
        raise FileNotFoundError(fullpath)
    data = xr_create_or_open(fullpath, taxis, name="REF")
    return data.CCF.to_dataframe()


def xr_save_mwcs2(station1, station2, components, filterid, mwcsid, mov_stack, taxis, dataframe):
    fn = os.path.join("DVV/MWCS/MWCS", "f%02i" % filterid, "mwcs%02i" % mwcsid,
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))
    if not os.path.isdir(os.path.split(fn)[0]):
        os.makedirs(os.path.split(fn)[0], exist_ok=True)
    d = dataframe.stack(future_stack=True).stack(future_stack=True)
    d.index = d.index.set_names(["times", "keys", "taxis"])
    d = d.reorder_levels(["times", "taxis", "keys"])
    d.columns = ["MWCS"]
    taxis = np.unique(d.index.get_level_values('taxis'))
    dr = xr_create_or_open(fn, taxis=taxis, name="MWCS")
    rr = d.to_xarray().to_dataset(name="MWCS")
    rr = xr_insert_or_update(dr, rr)
    xr_save_and_close(rr, fn)
    del dr, rr, d

def xr_get_mwcs2(station1, station2, components, filterid, mwcsid, mov_stack):
    fn = os.path.join("DVV/MWCS/MWCS", "f%02i" % filterid, "mwcs%02i" % mwcsid,
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))
    if not os.path.isfile(fn):
        # logging.error("FILE DOES NOT EXIST: %s, skipping" % fn)
        raise FileNotFoundError(fn)
    data = xr_create_or_open(fn, name="MWCS")
    data = data.MWCS.to_dataframe().reorder_levels(['times', 'taxis', 'keys']).unstack().droplevel(0, axis=1).unstack()
    return data


def xr_save_dtt2(station1, station2, components, filterid, mwcsid, dttid, mov_stack, dataframe):
    """
    :param station1: string, name of station 1
    :param station2: string, name of station 2
    :param components: string, name of the components
    :param filterid: int, filter id
    :param mov_stack: int, number of days in the moving stack
    :param dataframe: pandas DataFrame containing the data
    :return: None

    This method saves the given data in a NetCDF file using the specified parameters. The file path is constructed based on the station names, components, filter id, and moving stack number
    *. The data in the DataFrame is stacked, and the index is set to include "times" and "keys" as names. The column in the DataFrame is renamed to "DTT". A new or existing NetCDF file is
    * opened using the given file path, and the stacked data is inserted or updated in the file. The resulting dataset is then saved and the file is closed.
    """
    fn = os.path.join("DVV/MWCS/DTT", "f%02i" % filterid,
                      "mwcs%02i" % mwcsid, "dtt%02i" % dttid,
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))
    if not os.path.isdir(os.path.split(fn)[0]):
        os.makedirs(os.path.split(fn)[0], exist_ok=True)
    d = dataframe.stack(future_stack=True)
    d.index = d.index.set_names(["times", "keys"])
    d.columns = ["DTT"]
    dr = xr_create_or_open(fn, taxis=[], name="DTT")
    rr = d.to_xarray().to_dataset(name="DTT")
    rr = xr_insert_or_update(dr, rr)
    xr_save_and_close(rr, fn)

def xr_get_dtt2(station1, station2, components, filterid, mwcsid, dttid, mov_stack):
    """
    :param station1: The first station name
    :param station2: The second station name
    :param components: The components to be used
    :param filterid: The filter ID
    :param mov_stack: The movement stack
    :return: The extracted data

    This method retrieves the DTT data from a NetCDF file based on the given inputs. It constructs the file path using the provided parameters and checks if the file exists. If the file
    * does not exist, it raises a FileNotFoundError. Otherwise, it opens the NetCDF file and extracts the DTT variable as a dataframe. The dataframe is then rearranged and returned as the
    * result.
    """
    fn = os.path.join("DVV/MWCS/DTT", "f%02i" % filterid,
                      "mwcs%02i" % mwcsid, "dtt%02i" % dttid,                      
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))
    if not os.path.isfile(fn):
        # logging.error("FILE DOES NOT EXIST: %s, skipping" % fn)
        raise FileNotFoundError(fn)
    dr = xr_create_or_open(fn, taxis=[], name="DTT")
    data = dr.DTT.to_dataframe().reorder_levels(['times', 'keys']).unstack().droplevel(0, axis=1)
    return data

def xr_save_dvv2(components, filterid, mwcsid, dttid, mov_stack, dataframe):
    fn = os.path.join("DVV/MWCS/DVV", "f%02i" % filterid,
                      "mwcs%02i" % mwcsid, "dtt%02i" % dttid,
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s.nc" % components)
    if not os.path.isdir(os.path.split(fn)[0]):
        os.makedirs(os.path.split(fn)[0], exist_ok=True)

    if dataframe.columns.nlevels > 1:
        d = dataframe.stack(future_stack=True).stack(future_stack=True)
    else:
        d = dataframe.stack(future_stack=True)
    
    level_names = ["times", "level1", "level0"]
    d.index = d.index.set_names(level_names[:d.index.nlevels])

    if d.index.nlevels == 3:
        d = d.reorder_levels(["times", "level0", "level1"])

    d.columns = ["DVV"]
    # taxis = np.unique(d.index.get_level_values('taxis'))
    dr = xr_create_or_open(fn, taxis=[], name="DVV")
    rr = d.to_xarray().to_dataset(name="DVV")
    rr = xr_insert_or_update(dr, rr)
    xr_save_and_close(rr, fn)
    del dr, rr, d


def xr_get_dvv2(components, filterid, mwcsid, dttid, mov_stack):
    fn = os.path.join("DVV/MWCS/DVV", "f%02i" % filterid,
                      "mwcs%02i" % mwcsid, "dtt%02i" % dttid,
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s.nc" % components)
    if not os.path.isfile(fn):
        # logging.error("FILE DOES NOT EXIST: %s, skipping" % fn)
        raise FileNotFoundError(fn)
    data = xr_create_or_open(fn, name="DVV")
    data = data.DVV.to_dataframe().reorder_levels(['times', 'level1', 'level0']).unstack().droplevel(0, axis=1).unstack()
    return data

def wavg(group, dttname, errname):
    """
    Calculate the weighted average of a given group using the provided parameters.

    :param group: A pandas DataFrame or Series representing the group of data.
    :param dttname: The name of the column containing the data to be averaged.
    :param errname: The name of the column containing the error for each data point.
    :return: The weighted average of the data.

    """
    d = group[dttname]
    group.loc[group[errname] == 0,errname] = 1e-6
    w = 1. / group[errname]
    try:
        wavg = (d * w).sum() / w.sum()
    except:
        wavg = d.mean()
    return wavg


def wstd(group, dttname, errname):
    """
    :param group: A dictionary containing data for different groups.
    :param dttname: The key in the `group` dictionary that corresponds to the data array.
    :param errname: The key in the `group` dictionary that corresponds to the error array.
    :return: The weighted standard deviation of the data array.

    This method calculates the weighted standard deviation of the data array specified by `dttname` in the `group` dictionary.
    The weights are derived from the error array specified by `errname` in the `group` dictionary.

    The weighted standard deviation is computed using the following formula:
        wstd = sqrt(sum(w * (d - wavg) ** 2) / ((N - 1) * sum(w) / N))
    where:
        - d is the data array specified by `dttname` in the `group` dictionary.
        - w is the weight array derived from the error array specified by `errname` in the `group` dictionary.
        - wavg is the weighted average of the data array.
        - N is the number of non-zero weights.

    Note: This method uses the `np` module from NumPy.
    """
    d = group[dttname]
    group.loc[group[errname] == 0,errname] = 1e-6
    w = 1. / group[errname]
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w.values)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wstd


def get_wavgwstd(data, dttname, errname):
    """
    Calculate the weighted average and weighted standard deviation for a given data.

    :param data: The data to calculate the weighted average and weighted standard deviation.
    :type data: pandas.DataFrame

    :param dttname: The name of the column in the data frame containing the weights for the weighted average and weighted standard deviation calculation.
    :type dttname: str

    :param errname: The name of the column in the data frame containing the errors on the data.
    :type errname: str

    :return: A tuple containing the calculated weighted average and weighted standard deviation.
    :rtype: tuple
    """
    grouped = data.groupby(level=0)
    g = grouped.apply(wavg, dttname=dttname, errname=errname)
    h = grouped.apply(wstd, dttname=dttname, errname=errname)
    return g, h


def trim(data, dttname, limits=0.1):
    """
    Trimmed mean and standard deviation calculation.

    :param data: DataFrame containing the data.
    :param dttname: Name of the column used for grouping.
    :param limits: Trimming limits (default is 0.1).
    :return: Tuple containing the trimmed mean and trimmed standard deviation.
    """
    from scipy.stats.mstats import trimmed_mean, trimmed_std
    grouped = data[dttname].groupby(level=0)
    if limits == 0:
        g = grouped.mean()
        h = grouped.std()
    else:
        g = grouped.apply(trimmed_mean, limits=limits)
        h = grouped.apply(trimmed_std, limits=limits)
    return g, h


def compute_dvv(session, filterid, mov_stack, pairs=None, components=None, params=None, method=None, **kwargs):
    """
    Compute DVV statistics for a given seismic session, filter ID, and movement stack. The method calculates various
    statistical measures for DVV values based on the provided parameters.

    :param session: Seismic session object.
    :type session: object
    :param filterid: ID of the filter.
    :type filterid: str
    :param mov_stack: Movement stack object.
    :type mov_stack: object
    :param pairs: List of station pairs. Defaults to None.
    :type pairs: list, optional
    :param components: Component(s) to compute DVV for. Defaults to None.
    :type components: str, optional
    :param params: Parameters object containing configuration settings. Defaults to None.
    :type params: object, optional
    :param method: Method to use for computation. Defaults to None.
    :type method: str, optional
    :param **kwargs: Additional keyword arguments
    :return: DataFrame containing the computed statistical measures for DVV values.
    :rtype: DataFrame
    :raises ValueError: If no DVV values are computed.
    """
   
    if pairs == None:
        pairs = []
        for sta1, sta2 in get_station_pairs(session):
            for loc1 in sta1.locs():
                s1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
                for loc2 in sta2.locs():
                    s2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                    pairs.append((s1, s2))
    all = []
    for (s1, s2) in pairs:
        if components == None:
            if s1 == s2:
                # if not provided, we'll load all:
                comps = params.components_to_compute_single_station
            else:
                comps = params.components_to_compute
        else:
            if components.count(',') == 0:
                comps = [components, ]
            else:
                comps = components.split(',')

        for comp in comps:

            if (s1 == s2) and (comp not in params.components_to_compute_single_station):
                continue
            if (s1 != s2) and (comp not in params.components_to_compute):
                continue

            try:
                dtt = xr_get_dtt(s1, s2, comp, filterid, mov_stack)
                all.append(dtt)
            except FileNotFoundError:
                traceback.print_exc()
                continue
    if not len(all):
        raise ValueError
    if len(all) == 1:
        return all[0]
    all = pd.concat(all)
    percentiles = kwargs.get("percentiles", [.05, .10, .25, .5, .75, .90, .95])
    stats = all.groupby(level=0).describe(percentiles=percentiles)
    for c in ["m", "m0", "a"]:
        stats[(c, "weighted_mean")], stats[(c, "weighted_std")] = get_wavgwstd(all, c, 'e'+c)
        stats[(c, "trimmed_mean")], stats[(c, "trimmed_std")] = trim(all, c, kwargs.get("limits", None))

    return stats.sort_index(axis=1)

def compute_dvv2(session, filterid, mwcsid, dttid, mov_stack, pairs=None, components=None, params=None, method=None, **kwargs):
    """
    Compute DVV statistics for a given seismic session, filter ID, and movement stack. The method calculates various
    statistical measures for DVV values based on the provided parameters.

    :param session: Seismic session object.
    :type session: object
    :param filterid: ID of the filter.
    :type filterid: str
    :param mov_stack: Movement stack object.
    :type mov_stack: object
    :param pairs: List of station pairs. Defaults to None.
    :type pairs: list, optional
    :param components: Component(s) to compute DVV for. Defaults to None.
    :type components: str, optional
    :param params: Parameters object containing configuration settings. Defaults to None.
    :type params: object, optional
    :param method: Method to use for computation. Defaults to None.
    :type method: str, optional
    :param **kwargs: Additional keyword arguments
    :return: DataFrame containing the computed statistical measures for DVV values.
    :rtype: DataFrame
    :raises ValueError: If no DVV values are computed.
    """
    if pairs == None:
        pairs = []
        for sta1, sta2 in get_station_pairs(session):
            for loc1 in sta1.locs():
                s1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
                for loc2 in sta2.locs():
                    s2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                    pairs.append((s1, s2))
    all = []
    for (s1, s2) in pairs:
        if components == None:
            if s1 == s2:
                # if not provided, we'll load all:
                comps = params.components_to_compute_single_station
            else:
                comps = params.components_to_compute
        else:
            if components.count(',') == 0:
                comps = [components, ]
            else:
                comps = components.split(',')

        for comp in comps:

            if (s1 == s2) and (comp not in params.components_to_compute_single_station):
                continue
            if (s1 != s2) and (comp not in params.components_to_compute):
                continue

            try:
                dtt = xr_get_dtt2(s1, s2, comp, filterid, mwcsid, dttid, mov_stack)
                all.append(dtt)
            except FileNotFoundError:
                traceback.print_exc()
                continue
    if not len(all):
        raise ValueError
    if len(all) == 1:
        return all[0]
    all = pd.concat(all)
    percentiles = kwargs.get("percentiles", [.05, .10, .25, .5, .75, .90, .95])
    stats = all.groupby(level=0).describe(percentiles=percentiles)
    for c in ["m", "m0", "a"]:
        stats[(c, "weighted_mean")], stats[(c, "weighted_std")] = get_wavgwstd(all, c, 'e'+c)
        stats[(c, "trimmed_mean")], stats[(c, "trimmed_std")] = trim(all, c, kwargs.get("limits", None))

    return stats.sort_index(axis=1)

def xr_save_wct2(station1, station2, components, filterid, wctid, mov_stack, taxis, freqs, WXamp_list, WXcoh_list, WXdt_list, dates_list):
    """
    Save WCT results into an xarray Dataset and store it as a NetCDF file.

    Parameters:
    - station1, station2: str, Station pair
    - components: str, Seismic component (e.g., ZZ)
    - filterid: int, Filter reference ID
    - wctid: int, WCT parameter ID
    - mov_stack: tuple, Moving stack window (e.g., ('1d', '1d'))
    - taxis, freqs: np.array, Time axis and frequency axis
    - WXamp_list, Wcoh_list, WXdt_list: list of np.array, WCT outputs
    - dates_list: list of datetime, Timestamps for each WCT calculation
    """

    # Convert lists to xarray DataArrays
    WXamp_da = xr.DataArray(
        data=np.array(WXamp_list),
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="WXamp"
    )

    Wcoh_da = xr.DataArray(
        data=np.array(WXcoh_list),
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="Wcoh"
    )

    WXdt_da = xr.DataArray(
        data=np.array(WXdt_list),
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="WXdt"
    )

    # Combine into an xarray Dataset
    ds = xr.Dataset({"WXamp": WXamp_da, "Wcoh": Wcoh_da, "WXdt": WXdt_da})

    # Define output directory
    fn = os.path.join("DVV/WCT/WCT", "f%02i" % filterid, "wct%02i" % wctid, 
                            "%s_%s" % (mov_stack[0], mov_stack[1]), "%s" % components, 
                            f"{station1}_{station2}.nc")
    
    os.makedirs(os.path.dirname(fn), exist_ok=True)

    # Save to NetCDF
    xr_save_and_close(ds, fn)
    
    print(f"Saved WCT data for {station1}-{station2} to {fn}")

    # Cleanup memory
    del ds, WXamp_da, Wcoh_da, WXdt_da

def xr_load_wct(station1, station2, components, filterid, wctid, mov_stack):
    """
    Load WCT results from an xarray Dataset stored in a NetCDF file.

    Parameters:
    - station1, station2: str, Station pair
    - components: str, Seismic component (e.g., ZZ)
    - filterid: int, Filter reference ID
    - wctid: int, WCT parameter ID
    - mov_stack: tuple, Moving stack window (e.g., ('1d', '1d'))

    Returns:
    - ds: xarray.Dataset containing the WCT data (WXamp, Wcoh, WXdt)
    """

    # Construct the file path
    fn = os.path.join(
        "DVV/WCT/WCT", f"f{filterid:02d}", f"wct{wctid:02d}",
        f"{mov_stack[0]}_{mov_stack[1]}", components,
        f"{station1}_{station2}.nc"
    )

    # Check if the file exists
    if not os.path.exists(fn):
        raise FileNotFoundError(f"File not found: {fn}")

    # Load and return the dataset
    ds = xr.load_dataset(fn)
    return ds

def xr_save_wct_dtt2(station1, station2, components, filterid, wctid, dttid, mov_stack, taxis, dvv_df, err_df, coh_df):
    """
    Save the Wavelet Coherence Transform (WCT) results as a NetCDF file.

    :param station1: The first station in the pair.
    :type station1: str
    :param station2: The second station in the pair.
    :type station2: str
    :param components: The components (e.g., Z, N, E) being analyzed.
    :type components: str
    :param filterid: Filter ID used in the analysis.
    :type filterid: int
    :param mov_stack: Tuple of (start, end) representing the moving stack window.
    :type mov_stack: tuple
    :param taxis: Time axis corresponding to the WCT data.
    :type taxis: array-like
    :param dvv_df: DataFrame containing dvv data (2D).
    :type dvv_df: pandas.DataFrame
    :param err_df: DataFrame containing err data (2D).
    :type err_df: pandas.DataFrame
    :param coh_df: DataFrame containing coh data (2D).
    :type coh_df: pandas.DataFrame
    :returns: None
    """

    # Construct the file path
    fn = os.path.join(
        "DVV/WCT/DTT", f"f{filterid:02d}", f"wct{wctid:02d}", f"dtt{dttid:02d}",
        f"{mov_stack[0]}_{mov_stack[1]}", components,
        f"{station1}_{station2}.nc"
    )

    # Ensure the directory exists
    os.makedirs(os.path.dirname(fn), exist_ok=True)

    # Convert DataFrames to xarray.DataArrays
    dvv_da = xr.DataArray(dvv_df.values, coords=[dvv_df.index, dvv_df.columns], dims=['times', 'frequency'])
    err_da = xr.DataArray(err_df.values, coords=[err_df.index, err_df.columns], dims=['times', 'frequency'])
    coh_da = xr.DataArray(coh_df.values, coords=[coh_df.index, coh_df.columns], dims=['times', 'frequency'])

    # Combine into a single xarray.Dataset
    ds = xr.Dataset({
        'dvv': dvv_da,
        'err': err_da,
        'coh': coh_da
    })

    existing_ds = xr_create_or_open(fn, name="WCT")
    updated_ds = xr_insert_or_update(existing_ds, ds)
    xr_save_and_close(updated_ds, fn)


    logging.debug(f"Saved WCT data to {fn}")
    # Clean up
    del dvv_da, err_da, coh_da, ds

def filter_within_daterange(date, start_date, end_date):
    """Check if a date falls within the configured range"""
    return start_date <= date <= end_date
