import collections
import copy
import datetime
import itertools
import logging
import os
import glob
import traceback
import pickle
import math
import time

from logbook import Logger, StreamHandler
import sys

import numpy as np
import pandas as pd
import xarray as xr

from . import DBConfigNotFoundError
from .msnoise_table_def import (Job, Station, Config, DataAvailability,
                                WorkflowStep, WorkflowLink,
                                WORKFLOW_CHAINS, WORKFLOW_ORDER)


# ============================================================
# Section 1 — Infrastructure: Logging, DB, Configuration
# ============================================================


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
    :returns: :class:`~sqlalchemy.orm.session.Session` object, needed for
        many of the other API methods.
    """
    from sqlalchemy.orm import sessionmaker
    if not inifile:
        inifile = os.path.join(os.getcwd(), 'db.ini')

    engine = get_engine(inifile)
    session = sessionmaker(bind=engine)
    return session()


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
    pickle.dump([tech, hostname, database, username, password, prefix], f,
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
        tech, hostname, database, username, password, prefix = pickle.load(f)
    except:
        # Old ini file without prefix
        tech, hostname, database, username, password = pickle.load(f)
        prefix = ""
    f.close()
    return IniFile(tech, hostname, database, username, password, prefix)

# ── Configuration ──────────────────────────────────────────


def get_config(session, name=None, isbool=False, plugin=None, category='global', set_number=None):
    """Get the value of one or all config bits from the database.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type name: str
    :param name: The name of the config bit to get. If omitted, a dictionary
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
        from importlib.metadata import entry_points
        for ep in entry_points(group='msnoise.plugins.table_def'):
            if ep.name.replace("Config", "") == plugin:
                table = ep.load()
    else:
        table = Config
    if name:
        query = session.query(table).filter(table.name == name)
        if hasattr(table, 'category') and category is not None:
            query = query.filter(table.category == category)
            if category == 'global':
                set_number = 1
        if hasattr(table, 'set_number') and set_number is not None:
            query = query.filter(table.set_number == set_number)
        config = query.first()
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


def update_config(session, name, value, plugin=None, category='global', set_number=None):
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

    :type category: str
    :param category: The config category (default 'global'). Use e.g. 'stack',
        'filter', 'mwcs' to target a specific config set.

    :type set_number: int or None
    :param set_number: The config set number within the category (default None
        for global config). Use e.g. 1 for stack_1.

    """
    if plugin:
        from importlib.metadata import entry_points
        for ep in entry_points(group='msnoise.plugins.table_def'):
            if ep.name.replace("Config", "") == plugin:
                table = ep.load()
    else:
        table = Config
    query = session.query(table).filter(table.name == name)
    if hasattr(table, 'category') and category is not None:
        query = query.filter(table.category == category)
    if hasattr(table, 'set_number') and set_number is not None:
        query = query.filter(table.set_number == set_number)
    config = query.first()
    if config is None:
        return
    if "NULL" in str(value):
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


def get_config_set_details(session, set_name, set_number, format="list"):
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

        if format == "list":
            return [{'name': c.name, 'value': c.value, 'description': c.description}
             for c in configs]
        elif format == "AttribDict":
            from obspy.core import AttribDict
            import pydoc
            attrib_dict = AttribDict()
            for c in configs:
                if c.value is None or c.value == '':
                    if c.param_type == 'str':
                        attrib_dict[c.name] = ''
                    # else: skip — can't convert empty string to int/float/eval
                    continue
                itemtype = pydoc.locate(c.param_type)
                if itemtype is bool:
                    if c.value in [True, 'True', 'true', 'Y', 'y', '1', 1]:
                        attrib_dict[c.name] = True
                    else:
                        attrib_dict[c.name] = False
                else:
                    attrib_dict[c.name] = itemtype(c.value)
            return attrib_dict


    except Exception as e:
        from obspy.core import AttribDict
        get_logger("msnoise", loglevel="ERROR").error(f"Failed to get config set details: {str(e)}")
        return AttribDict() if format == "AttribDict" else []


def get_config_categories_definition():
    """Get the standard configuration categories with display names, order, and indent level.

    Each entry is ``(category_key, display_name, level)`` where *level* is the depth
    relative to ``global`` (0).  Used by the config-sets admin page for visual
    indentation that mirrors the workflow graph.
    """
    return [
        # Tree order: children immediately follow their parent
        ('global',      'Global Parameters',  0),
        ('preprocess',  'Preprocessing',       1),
        ('cc',          'Cross-Correlation',   2),
        ('filter',      'Filters',             3),
        ('stack',       'Moving Stacks',       4),
        ('refstack',    'Reference Stacks',    5),
        ('mwcs',          'MWCS',                    6),
        ('mwcs_dtt',      'MWCS dt/t',               7),
        ('mwcs_dtt_dvv',  'MWCS dv/v Aggregate',     8),
        ('stretching',    'Stretching',               6),
        ('stretching_dvv','Stretching dv/v Aggregate',7),
        ('wavelet',       'Wavelet',                  6),
        ('wavelet_dtt',   'Wavelet dt/t',             7),
        ('wct_dtt_dvv',   'WCT dv/v Aggregate',       8),
        ('psd',         'PSD',                 1),
        ('psd_rms',     'PSD RMS',             2),
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

    # Migrate legacy NULL set_number rows (written by pre-refactor code
    # where global config had no set_number).  Fold them into set_number=1
    # so they are accessible via the normal workflow path.  Skip if a
    # set_number=1 row already exists for that category (no overwrite).
    null_categories = {cat for cat, sn, _ in all_sets if sn is None}
    existing_1 = {cat for cat, sn, _ in all_sets if sn == 1}
    for null_cat in null_categories - existing_1:
        session.query(schema.Config).filter(
            schema.Config.category == null_cat,
            schema.Config.set_number == None,  # noqa: E711
        ).update({"set_number": 1}, synchronize_session=False)
    if null_categories - existing_1:
        session.commit()
        # Re-query after migration
        all_sets = session.query(
            schema.Config.category,
            schema.Config.set_number,
            func.count(schema.Config.ref).label('param_count')
        ).group_by(
            schema.Config.category,
            schema.Config.set_number
        ).all()

    # Group sets by category, skipping any remaining NULL rows (e.g.
    # where set_number=1 already existed so migration was skipped).
    sets_by_category = {}
    for category, set_number, param_count in all_sets:
        if set_number is None:
            continue
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
    for category_key, display_name, level in category_order:
        if category_key in sets_by_category:
            ordered_categories.append({
                'category': category_key,
                'display_name': display_name,
                'level': level,
                'sets': sets_by_category[category_key]
            })

    # Add any additional categories not in the predefined order
    for category in sets_by_category:
        if category not in [cat[0] for cat in category_order]:
            ordered_categories.append({
                'category': category,
                'display_name': category.title(),
                'level': 0,
                'sets': sets_by_category[category]
            })

    return ordered_categories

# ── Parameters ─────────────────────────────────────────────


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
    # params.goal_sampling_rate = params.cc_sampling_rate
    # params.min30 = params.corr_duration * params.goal_sampling_rate
    # params.components_to_compute = get_components_to_compute(s)
    # params.components_to_compute_single_station = get_components_to_compute_single_station(s)
    # params.all_components = np.unique(params.components_to_compute_single_station + \
    #                         params.components_to_compute)

    # if not isinstance(params.mov_stack[0], tuple):
    #     params.mov_stack = [params.mov_stack, ]
    # else:
    #     params.mov_stack = params.mov_stack
    return params


# ============================================================
# Section 2 — Network / Station
# ============================================================


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
        station = Station(net, sta, X, Y, altitude, coordinates, used)
        session.add(station)
    else:
        station.X = X
        station.Y = Y
        station.altitude = altitude
        station.coordinates = coordinates
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


# ============================================================
# Section 3 — Data Availability
# ============================================================


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


# ============================================================
# Section 4 — Workflow Structure
# ============================================================


def get_workflow_steps(session):
    """Get all steps in a workflow"""
    from .msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowStep) \
        .filter(schema.WorkflowStep.is_active == True) \
        .order_by(schema.WorkflowStep.step_name).all()


def get_workflow_links(session):
    """Get all links in a workflow"""
    from .msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowLink) \
        .filter(schema.WorkflowLink.is_active == True).all()


def get_workflow_graph(session):
    """Return workflow as nodes and edges for visualization, sorted by workflow order"""
    steps = get_workflow_steps(session)
    links = get_workflow_links(session)

    # Sort steps by workflow order (WORKFLOW_ORDER imported from msnoise_table_def)
    def get_workflow_order_key(step):
        try:
            category_order = WORKFLOW_ORDER.index(step.category)
        except ValueError:
            # If category not in predefined order, put it at the end
            category_order = len(WORKFLOW_ORDER)

        # Sort by category first, then by set_number
        return (category_order, step.set_number or 0)

    sorted_steps = sorted(steps, key=get_workflow_order_key)

    nodes = []
    for step in sorted_steps:
        nodes.append({
            "id": step.step_id,
            "name": step.step_name,
            "category": step.category,
            "set_number": step.set_number,
            "description": step.description
        })

    edges = []
    for link in links:
        edges.append({
            "from": link.from_step_id,
            "to": link.to_step_id,
            "type": link.link_type
        })

    return {"nodes": nodes, "edges": edges}


def create_workflow_step(session, step_name, category, set_number, description=None):
    """Create a new workflow step"""
    from .msnoise_table_def import declare_tables
    schema = declare_tables()

    step = schema.WorkflowStep(
        step_name=step_name,
        category=category,
        set_number=set_number,
        description=description
    )

    session.add(step)
    session.commit()
    return step


def create_workflow_link(session, from_step_id, to_step_id, link_type="default"):
    """Create a link between two workflow steps"""
    from .msnoise_table_def import declare_tables
    schema = declare_tables()

    # Check if link already exists
    existing = session.query(schema.WorkflowLink).filter(
        schema.WorkflowLink.from_step_id == from_step_id,
        schema.WorkflowLink.to_step_id == to_step_id
    ).first()

    if existing:
        return existing

    link = schema.WorkflowLink(
        from_step_id=from_step_id,
        to_step_id=to_step_id,
        link_type=link_type
    )

    session.add(link)
    session.commit()
    return link


def get_step_successors(session, step_id):
    """Get all steps that this step feeds into"""
    from .msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowStep) \
        .join(schema.WorkflowLink, schema.WorkflowStep.step_id == schema.WorkflowLink.to_step_id) \
        .filter(schema.WorkflowLink.from_step_id == step_id) \
        .filter(schema.WorkflowLink.is_active == True).all()


def get_first_runnable_steps_per_branch(session, source_step_id, skip_categories=None):
    """
    Return the first runnable step on each outgoing branch from source_step_id,
    skipping steps in skip_categories (default: {"filter"}).

    - "Branch" roots are the immediate successors of source_step_id.
    - If skipped nodes fan out, each fan-out is treated as a separate branch.
    - Results are deduped (set behavior) and returned in stable order.
    """
    if skip_categories is None:
        skip_categories = {"filter"}

    # We dedupe by step_id, but preserve stable ordering at the end
    result_by_step_id = {}

    def immediate_successors(step_id):
        # Reuse existing helper (one-hop)
        return [
            s for s in get_step_successors(session, step_id)
            if getattr(s, "is_active", True)
        ]

    def descend_until_runnable(start_step):
        # Explore forward until the first non-skipped step is found for each path
        stack = [start_step]
        visited = set()

        while stack:
            step = stack.pop()
            if step.step_id in visited:
                continue
            visited.add(step.step_id)

            if step.category not in skip_categories:
                result_by_step_id[step.step_id] = step
                continue  # stop this path at first runnable

            for nxt in immediate_successors(step.step_id):
                stack.append(nxt)

    for succ in immediate_successors(source_step_id):
        descend_until_runnable(succ)

    return [result_by_step_id[k] for k in sorted(result_by_step_id)]


def _get_step_predecessors(session, step_id):
    """Get all steps that feed into this step"""
    from .msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowStep) \
        .join(schema.WorkflowLink, schema.WorkflowStep.step_id == schema.WorkflowLink.from_step_id) \
        .filter(schema.WorkflowLink.to_step_id == step_id) \
        .filter(schema.WorkflowLink.is_active == True).all()


def get_upstream_steps_for_step_id(session, step_id, topo_order=True, include_self=False):
    """
    Returns all upstream WorkflowStep nodes (recursive predecessors) for `step_id`.

    Parameters
    ----------
    session : sqlalchemy.orm.session.Session
    step_id : int
        The step_id of the node whose upstream lineage you want.
    topo_order : bool
        If True (default), returns nodes in dependency order (upstream first),
        suitable for "merge config" application.
        If False, returns in discovery order (still de-duplicated).
    include_self : bool
        If True, includes the node identified by `step_id` in the returned list.

    Returns
    -------
    list[WorkflowStep]
        De-duplicated list of upstream steps.
    """
    visited = set()
    ordered = []

    def dfs(current_step_id):
        # Direct predecessors of the current node:
        preds = _get_step_predecessors(session, current_step_id) or []
        for p in preds:
            if p.step_id in visited:
                continue
            visited.add(p.step_id)
            dfs(p.step_id)
            if topo_order:
                # post-order ensures: preprocess -> cc -> filter -> ...
                ordered.append(p)
            else:
                ordered.append(p)

    dfs(step_id)

    if include_self:
        self_step = session.query(WorkflowStep).filter(WorkflowStep.step_id == step_id).first()
        if self_step is not None and self_step.step_id not in visited:
            ordered.append(self_step)

    return ordered


def create_workflow_steps_from_config_sets(session):
    """
    Create workflow steps automatically from all existing config sets,
    sorted by natural workflow order.

    Returns:
        tuple: (created_count, existing_count, error_message)
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()

    try:
        # Get all unique category+set_number combinations
        config_sets = session.query(
            schema.Config.category,
            schema.Config.set_number
        ).filter(
            schema.Config.set_number.isnot(None)  # Exclude global configs
        ).distinct().all()

        # Sort by workflow order
        def get_workflow_order(config_set):
            category, set_number = config_set
            try:
                return WORKFLOW_ORDER.index(category)
            except ValueError:
                # If category not in predefined order, put it at the end
                return len(WORKFLOW_ORDER)

        config_sets = sorted(config_sets, key=get_workflow_order)
        created_count = 0
        existing_count = 0

        for category, set_number in config_sets:
            # Check if step already exists
            existing_step = session.query(schema.WorkflowStep).filter(
                schema.WorkflowStep.category == category,
                schema.WorkflowStep.set_number == set_number,
            ).first()

            if not existing_step:
                step_name = f"{category}_{set_number}"
                description = f"Auto-generated step for {category} configuration set {set_number}"

                create_workflow_step(
                    session,
                    step_name,
                    category,
                    set_number,
                    description
                )
                created_count += 1
            else:
                existing_count += 1

        return created_count, existing_count, None

    except Exception as e:
        return 0, 0, str(e)


def create_workflow_links_from_steps(session):
    """
    Create workflow links automatically between existing workflow steps,
    following natural workflow progression.

    Returns:
        tuple: (created_count, existing_count, error_message)
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()

    try:
        # Get all workflow steps
        steps = session.query(schema.WorkflowStep).all()

        # Group steps by category and set_number
        steps_by_category = {}
        for step in steps:
            category = step.category
            if category not in steps_by_category:
                steps_by_category[category] = {}
            steps_by_category[category][step.set_number] = step

        created_count = 0
        existing_count = 0

        # Create links based on workflow chains
        # WORKFLOW_CHAINS imported from msnoise_table_def; use simple next_steps list
        for source_category, target_categories in WORKFLOW_CHAINS.items():
            # msnoise_table_def WORKFLOW_CHAINS values are dicts with 'next_steps'
            if isinstance(target_categories, dict):
                target_categories = target_categories.get('next_steps', [])
            if source_category not in steps_by_category:
                continue

            # For each source step in this category
            for source_set_number, source_step in steps_by_category[source_category].items():

                # Link to each target category
                for target_category in target_categories:
                    if target_category not in steps_by_category:
                        continue

                    # Strategy: Create links based on workflow logic
                    target_steps_to_link = []

                    if source_category == 'global':
                        # Global steps link to all steps in target categories
                        target_steps_to_link = list(steps_by_category[target_category].values())

                    elif target_category in [
                        'filter', 'refstack',
                        'mwcs', 'stretching', 'wavelet',
                        'mwcs_dtt', 'wavelet_dtt',
                        'mwcs_dtt_dvv', 'stretching_dvv', 'wct_dtt_dvv',
                    ]:
                        # For all processing steps that can have multiple instances,
                        # link to all target steps in the category
                        # This includes DTT steps which can process results from multiple MWCS/wavelet sets
                        target_steps_to_link = list(steps_by_category[target_category].values())

                    else:
                        # Default: link to matching set numbers if available, otherwise to all
                        if source_set_number in steps_by_category[target_category]:
                            target_steps_to_link = [steps_by_category[target_category][source_set_number]]
                        else:
                            target_steps_to_link = list(steps_by_category[target_category].values())

                    # Create the links
                    for target_step in target_steps_to_link:
                        # Check if link already exists
                        existing_link = session.query(schema.WorkflowLink).filter(
                            schema.WorkflowLink.from_step_id == source_step.step_id,
                            schema.WorkflowLink.to_step_id == target_step.step_id
                        ).first()

                        if not existing_link:
                            create_workflow_link(
                                session,
                                source_step.step_id,
                                target_step.step_id,
                                'default'
                            )
                            created_count += 1
                        else:
                            existing_count += 1

        return created_count, existing_count, None

    except Exception as e:
        return 0, 0, str(e)

# Job management functions for workflow-aware jobs


def get_workflow_step_config(session, step_name):
    """Return the :class:`WorkflowStep` for *step_name*, or ``None``.

    :param session: Active SQLAlchemy session.
    :param step_name: Exact step name string, e.g. ``"preprocess_1"``.
    """
    for step in get_workflow_steps(session):
        if step.step_name == step_name:
            return step
    return None


# ============================================================
# Section 5 — Jobs
# ============================================================


def update_job(session, day, pair, jobtype, flag,
               step_id=None, priority=0, lineage=None,
               commit=True, returnjob=True, ref=None):
    """
    Updates or Inserts a :class:`~msnoise.msnoise_table_def.declare_tables.Job`
    in the database.  Workflow-aware: handles ``step_id`` and ``lineage`` fields.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type day: str
    :param day: The day in YYYY-MM-DD format
    :type pair: str
    :param pair: the name of the station pair
    :type jobtype: str
    :param jobtype: Job type string, e.g. ``"preprocess_1"``
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one, "F"ailed.
    :type step_id: int or None
    :param step_id: WorkflowStep primary key, or None
    :type priority: int
    :param priority: Job priority (default 0)
    :type lineage: str or None
    :param lineage: Lineage string encoding upstream configset chain
    :type commit: bool
    :param commit: Whether to directly commit (True, default) or not (False)
    :type returnjob: bool
    :param returnjob: Return the modified/inserted Job (True, default) or not (False)
    :type ref: int or None
    :param ref: If provided, look up the job by its primary key instead of
        (day, pair, jobtype, lineage).

    :rtype: :class:`~msnoise.msnoise_table_def.declare_tables.Job` or None
    :returns: If returnjob is True, returns the modified/inserted Job.
    """
    from sqlalchemy import text

    if ref:
        job = session.query(Job).filter(text("ref=:ref")).params(ref=ref).first()
    else:
        job = (session.query(Job)
               .filter(text("day=:day"))
               .filter(text("pair=:pair"))
               .filter(text("jobtype=:jobtype"))
               .filter(text("lineage=:lineage"))
               .params(day=day, pair=pair, jobtype=jobtype, lineage=lineage)
               .first())

    if job is None:
        job = Job()
        job.day = day
        job.pair = pair
        job.jobtype = jobtype
        job.step_id = step_id
        job.priority = priority
        job.flag = flag
        job.lastmod = datetime.datetime.utcnow()
        job.lineage = lineage
        session.add(job)
    else:
        # Never demote a Done job back to Todo
        if not (job.flag == "D" and flag == "T"):
            job.flag = flag
        job.step_id = step_id
        job.priority = priority
        job.lastmod = datetime.datetime.utcnow()
        job.lineage = lineage

    if commit:
        session.commit()
    if returnjob:
        return job


def massive_insert_job(session, jobs):
    """
    Bulk-insert a list of job dicts into the jobs table.

    Each dict must contain at minimum ``day``, ``pair``, ``jobtype``,
    ``flag`` and ``lastmod`` keys.  Optional keys: ``step_id``, ``priority``,
    ``lineage``.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobs: list[dict]
    :param jobs: Job records to insert.
    """
    job_records = [
        {
            'day':      j['day'],
            'pair':     j['pair'],
            'jobtype':  j['jobtype'],
            'step_id':  j.get('step_id'),
            'priority': j.get('priority', 0),
            'flag':     j['flag'],
            'lastmod':  j['lastmod'],
            'lineage':  j.get('lineage'),
        }
        for j in jobs
    ]
    session.bulk_insert_mappings(Job, job_records)
    session.commit()


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


def reset_jobs(session, jobtype, alljobs=False, reset_i=True, reset_e=True):
    """Reset jobs with the given ``jobtype`` string back to "T"odo.

    Works with the v2 workflow model where ``jobtype`` is a step name such
    as ``"cc_1"`` or ``"refstack_1"``.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobtype: str
    :param jobtype: Step name to reset (e.g. ``"cc_1"``)
    :type alljobs: bool
    :param alljobs: If True reset all jobs regardless of current flag;
        otherwise only resets "I" and/or "E" flagged jobs.
    :type reset_i: bool
    :param reset_i: Reset "I"n-progress jobs (default True)
    :type reset_e: bool
    :param reset_e: Reset "E"rror/failed jobs (default True)
    """
    from sqlalchemy import update as sa_update
    q = sa_update(Job).where(Job.jobtype == jobtype)
    if alljobs:
        session.execute(q.values(flag='T'))
    else:
        flags = []
        if reset_i:
            flags.append('I')
        if reset_e:
            flags.append('F')
        if flags:
            session.execute(q.where(Job.flag.in_(flags)).values(flag='T'))
    session.commit()


def get_job_types(session, jobtype):
    """Return job counts grouped by flag for a given ``jobtype`` string.

    Works with the v2 workflow model where ``jobtype`` is a step name such
    as ``"cc_1"``.  Returns a list of ``(count, flag)`` tuples, matching the
    interface expected by the test suite.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobtype: str
    :param jobtype: Step name to query (e.g. ``"cc_1"``)
    :rtype: list of (int, str)
    :returns: List of (count, flag) pairs
    """
    from sqlalchemy import func
    rows = (session.query(func.count(Job.flag), Job.flag)
            .filter(Job.jobtype == jobtype)
            .group_by(Job.flag)
            .all())
    return rows


# ── Workflow-aware job API ──────────────────────────────────


def get_next_job_for_step(
        session,
        step_category="preprocess",
        flag="T",
        group_by="day",
        limit_days=None,
):
    """
    Return a claimed batch of jobs for a workflow step category.

    group_by:
      - "day": claim all jobs for the selected (step_id, jobtype, day)
      - "pair": claim all jobs for the selected (step_id, jobtype, pair)
      - "pair_lineage": claim all jobs for the selected (step_id, jobtype, pair, lineage)
      - "day_lineage": claim all jobs for the selected (step_id, jobtype, day, lineage)
    """
    from .msnoise_table_def import declare_tables
    from sqlalchemy import update

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    # refstack jobs always have day="REF"; all other categories never do.
    if step_category == "refstack":
        day_filter = (Job.day == "REF")
    else:
        day_filter = (Job.day != "REF")

    next_job = (
        session.query(Job, WorkflowStep)
        .join(WorkflowStep, Job.jobtype == WorkflowStep.step_name)
        .filter(WorkflowStep.category == step_category)
        .filter(Job.flag == flag)
        .filter(day_filter)
        .order_by(Job.priority.desc(), Job.lastmod)
        .first()
    )
    if not next_job:
        return [], None

    seed_job, step = next_job
    batch_q = (
        session.query(Job)
        .filter(Job.step_id == step.step_id)
        .filter(Job.jobtype == seed_job.jobtype)
        .filter(day_filter)
        .filter(Job.flag == flag)
    )

    if group_by == "day":
        batch_q = batch_q.filter(Job.day == seed_job.day)
        batch_q = batch_q.order_by(Job.pair.asc(), Job.lineage.asc())
    elif group_by == "pair":
        batch_q = batch_q.filter(Job.pair == seed_job.pair)
        batch_q = batch_q.order_by(Job.day.asc(), Job.lineage.asc())
    elif group_by == "pair_lineage":
        batch_q = batch_q.filter(Job.pair == seed_job.pair).filter(Job.lineage == seed_job.lineage)
        batch_q = batch_q.order_by(Job.day.asc())
        if limit_days is not None:
            batch_q = batch_q.limit(int(limit_days))
    elif group_by == "day_lineage":
        batch_q = batch_q.filter(Job.day == seed_job.day).filter(Job.lineage == seed_job.lineage)
        batch_q = batch_q.order_by(Job.pair.asc())
    else:
        raise ValueError(f"Unsupported group_by={group_by!r}")

    jobs = batch_q.with_for_update().all()
    if not jobs:
        return [], step

    upd = (
        update(Job)
        .where(Job.step_id == step.step_id)
        .where(Job.jobtype == seed_job.jobtype)
        .where(day_filter)
        .where(Job.flag == flag)
    )

    if group_by == "day":
        upd = upd.where(Job.day == seed_job.day)
    elif group_by == "pair":
        upd = upd.where(Job.pair == seed_job.pair)
    elif group_by == "pair_lineage":
        upd = upd.where(Job.pair == seed_job.pair).where(Job.lineage == seed_job.lineage)
        if limit_days is not None:
            # Note: SQL UPDATE with LIMIT is not portable. Keep limit_days for SELECT batching only.
            # The UPDATE still claims the whole (pair, lineage) batch.
            pass
    else:  # "day_lineage"
        upd = upd.where(Job.day == seed_job.day).where(Job.lineage == seed_job.lineage)

    session.execute(upd.values(flag="I"))
    session.commit()

    return jobs, step


def is_next_job_for_step(session, step_category="preprocess", flag='T'):
    """
    Are there any workflow-aware Jobs in the database with the specified
    flag and step category?

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type step_category: str
    :param step_category: Workflow step category (e.g., "preprocess", "qc")
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.

    :rtype: bool
    :returns: True if at least one Job matches the criteria, False otherwise.
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    if step_category == "refstack":
        day_filter = (Job.day == "REF")
    else:
        day_filter = (Job.day != "REF")

    job = session.query(Job) \
        .join(WorkflowStep, Job.jobtype == WorkflowStep.step_name) \
        .filter(WorkflowStep.category == step_category) \
        .filter(Job.flag == flag) \
        .filter(day_filter) \
        .first()

    if job is None:
        return False
    else:
        return True


def get_workflow_job_counts(db):
    """
    Get job counts by status for the workflow

    :param db: Database connection
    :return: Dictionary with job counts by status
    """
    from sqlalchemy import func

    counts = db.query(
        Job.flag,
        func.count(Job.ref).label('count')
    ).group_by(Job.flag).all()

    result = {'T': 0, 'I': 0, 'D': 0}
    for flag, count in counts:
        result[flag] = count

    return result


# ============================================================
# Section 6 — Lineage Resolution
# ============================================================


def get_filter_steps_for_cc_step(session, cc_step_id):
    """
    Get all filter steps that are children of a specific CC step.

    :param session: Database session
    :param cc_step_id: The step_id of the CC step
    :return: List of filter workflow steps that are successors of the CC step
    """
    from .api import get_step_successors

    # Get all steps that are successors of the CC step
    successor_steps = get_step_successors(session, cc_step_id)

    # Filter to only include steps with category "filter"
    filter_steps = [step for step in successor_steps if step.category == "filter"]

    return filter_steps


def lineage_to_string(lineage_steps):
    return "/".join(s.step_name for s in lineage_steps if s.step_name)


def lineage_str_to_step_names(lineage_str, sep="/"):
    """
    Convert a lineage string like "preprocess_1/cc_1/filter_1" into a list of
    step_name strings, preserving order.
    """
    if lineage_str is None:
        raise ValueError("lineage_str is None")
    lineage_str = str(lineage_str).strip()
    if not lineage_str:
        return []
    return [p.strip() for p in lineage_str.split(sep) if p.strip()]


def lineage_str_to_steps(session, lineage_str, sep="/", strict=True):
    """
    Resolve a lineage string to a list of WorkflowStep ORM objects (ordered).

    Parameters
    ----------
    session : sqlalchemy.orm.session.Session
    lineage_str : str
        e.g. "preprocess_1/cc_1/filter_1"
    sep : str
        Separator used in lineage strings.
    strict : bool
        If True, raises if any step_name can't be resolved.
        If False, silently skips missing steps.

    Returns
    -------
    list[WorkflowStep]
    """
    names = lineage_str_to_step_names(lineage_str, sep=sep)
    if not names:
        return []

    rows = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.step_name.in_(names))
        .all()
    )
    by_name = {s.step_name: s for s in rows}

    steps = []
    missing = []
    for name in names:
        s = by_name.get(name)
        if s is None:
            missing.append(name)
            if not strict:
                continue
        else:
            steps.append(s)

    if missing and strict:
        raise ValueError(
            "Could not resolve lineage steps: %s"
            % ", ".join(missing)
        )

    return steps


def get_lineages_to_step_id(
    session,
    step_id,
    include_self=True,
    max_depth=50,
    max_paths=5000,
):
    """
    Enumerate upstream lineages (all distinct paths) ending at `step_id`.

    Returns a list of paths, each path being a list[WorkflowStep] ordered
    upstream -> downstream. This preserves branch structure (unlike
    get_upstream_steps_for_step_id which de-duplicates nodes).

    Safety:
      - max_depth prevents infinite loops in case of bad/cyclic graphs
      - max_paths prevents combinatorial explosion

    Parameters
    ----------
    session : sqlalchemy.orm.session.Session
    step_id : int
    include_self : bool
        If True, each path ends with the step itself.
    max_depth : int
    max_paths : int

    Returns
    -------
    list[list[WorkflowStep]]
    """
    # Resolve node objects on demand
    def get_step(sid):
        return session.query(WorkflowStep).filter(WorkflowStep.step_id == sid).first()

    # We assume the workflow graph is a DAG. If it isn't, we guard with max_depth
    # and a per-path "seen" set.
    paths = []

    def dfs(current_id, suffix_path, seen_ids, depth):
        if depth > max_depth:
            raise RuntimeError(f"Exceeded max_depth={max_depth} while expanding lineage for step_id={step_id}")

        preds = _get_step_predecessors(session, current_id) or []

        # Leaf: no predecessors => we have a complete path
        if not preds:
            # suffix_path currently holds downstream part (from current -> ... -> target)
            final_path = list(reversed(suffix_path))
            paths.append(final_path)
            if len(paths) > max_paths:
                raise RuntimeError(f"Exceeded max_paths={max_paths} while expanding lineage for step_id={step_id}")
            return

        for p in preds:
            if p.step_id in seen_ids:
                # cycle protection (shouldn't happen in a DAG)
                continue
            dfs(
                p.step_id,
                suffix_path + [p],
                seen_ids | {p.step_id},
                depth + 1,
            )

    start = get_step(step_id)
    if start is None:
        return []

    if include_self:
        dfs(step_id, [start], {step_id}, 0)
    else:
        dfs(step_id, [], set(), 0)

    return paths


def get_done_lineages_for_category(session, category):
    """Return all distinct computed lineages for a given workflow step category.

    Queries ``Job`` rows whose associated ``WorkflowStep.category`` matches
    *category* and whose flag is ``'D'`` (done), then de-duplicates and
    resolves each ``lineage`` string into an ordered list of step-name
    strings (upstream → downstream, including the step itself).

    This is the correct way to enumerate output folders for a step that may
    be reached via multiple upstream paths (e.g. multiple filters, multiple
    MWCS configs), because it reflects what was *actually computed* rather
    than what the DAG topology suggests.

    Example::

        lineages = get_done_lineages_for_category(db, 'stretching')
        # → [
        #     ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'stretching_1'],
        #     ['preprocess_1', 'cc_1', 'filter_2', 'stack_1', 'stretching_1'],
        # ]

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param category: Workflow step category, e.g. ``'stretching'``,
        ``'mwcs_dtt'``, ``'wavelet_dtt'``.
    :type category: str
    :rtype: list[list[str]]
    :returns: Sorted list of unique lineage-name lists.
    """
    rows = (
        session.query(Job.lineage)
        .join(Job.workflow_step)
        .filter(WorkflowStep.category == category)
        .filter(Job.flag == "D")
        .filter(Job.lineage.isnot(None))
        .distinct()
        .all()
    )
    seen = set()
    result = []
    for (lineage_str,) in rows:
        if not lineage_str or lineage_str in seen:
            continue
        seen.add(lineage_str)
        names = lineage_str_to_step_names(lineage_str)
        # Guard: skip empty or single-entry lineages (global-only, no real steps)
        if len(names) >= 2:
            result.append(names)
    result.sort()
    return result


def _load_step_config(db, step):
    # step must contain config_category + config_set_number (or equivalent)
    return get_config_set_details(
        db,
        step.category,
        step.set_number,
        format="AttribDict",
    )


def _merge_params(orig_params, configs_in_order):
    from obspy.core import AttribDict
    merged = dict(orig_params)
    for cfg in configs_in_order:
        merged.update(cfg)  # later overrides earlier
    return AttribDict(merged)


def get_merged_params_for_lineage(db, orig_params, step_params, lineage):
    # lineage is upstream -> downstream
    lineage = [s for s in lineage if s.category not in {"global"}]

    lineage_names = [s.step_name for s in lineage]
    # print("Lineage Names:", lineage_names)

    lineage_cfgs = [_load_step_config(db, s) for s in lineage]
    # print("Lineage Configs:", lineage_cfgs)

    # TODO: gather all those into a "sanitize params" def?
    params = _merge_params(orig_params, lineage_cfgs + [step_params])

    # Split comma-separated component strings only when the param is present
    # and not already a list (preprocess/psd steps don't carry these keys).
    for key in ('components_to_compute', 'components_to_compute_single_station'):
        val = getattr(params, key, None)
        if val is not None and isinstance(val, str):
            setattr(params, key, val.split(','))

    if hasattr(params, "mov_stack"):
        if not isinstance(params.mov_stack[0], tuple):
            params.mov_stack = [params.mov_stack]

    return lineage, lineage_names, params


def resolve_lineage_params(session, lineage_names):
    """Resolve a lineage name-list into a fully merged params object.

    Given a list of step-name strings (as returned by
    :func:`get_done_lineages_for_category`), resolves them to
    :class:`~msnoise.msnoise_table_def.WorkflowStep` ORM objects and merges
    every step's configuration into the global params, exactly as the
    processing steps themselves do via :func:`get_next_lineage_batch`.

    Returns ``(lineage_steps, lineage_names, params)`` — the same tuple as
    :func:`get_merged_params_for_lineage` — so callers can use ``params``
    directly (it will have ``components_to_compute``, ``mov_stack``, etc.).

    Example::

        lineage_names = get_done_lineages_for_category(db, 'mwcs_dtt')[0]
        _, _, params = resolve_lineage_params(db, lineage_names)
        mov_stack = params.mov_stack[0]

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param lineage_names: Ordered list of step-name strings, e.g.
        ``['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'mwcs_1', 'mwcs_dtt_1']``.
    :rtype: tuple(list, list[str], AttribDict)
    """
    orig_params = get_params(session)
    lineage_str = "/".join(lineage_names)
    steps = lineage_str_to_steps(session, lineage_str, sep="/", strict=False)
    return get_merged_params_for_lineage(session, orig_params, {}, steps)


# ============================================================
# Stretching xarray I/O
# ============================================================


def resolve_lineage_from_ids(session, params, preprocess_id=1, cc_id=None,
                              filter_id=None, stack_id=None, refstack_id=None,
                              mwcs_id=None, mwcs_dtt_id=None,
                              stretching_id=None, wavelet_id=None,
                              wavelet_dtt_id=None):
    """Build a lineage name list from integer step-set IDs and resolve params.

    Constructs step-name strings (e.g. ``"preprocess_1"``, ``"cc_2"``) from
    the supplied integer IDs in workflow order, then delegates to
    :func:`resolve_lineage_params` to merge all step configurations.
    Only IDs that are not ``None`` are included.

    :returns: ``(lineage_steps, lineage_names, merged_params)`` - same tuple
        as :func:`get_merged_params_for_lineage`.
    """
    parts = [f"preprocess_{preprocess_id}"]
    if cc_id is not None:
        parts.append(f"cc_{cc_id}")
    if filter_id is not None:
        parts.append(f"filter_{filter_id}")
    if stack_id is not None:
        parts.append(f"stack_{stack_id}")
    if refstack_id is not None:
        parts.append(f"refstack_{refstack_id}")
    if mwcs_id is not None:
        parts.append(f"mwcs_{mwcs_id}")
    if mwcs_dtt_id is not None:
        parts.append(f"mwcs_dtt_{mwcs_dtt_id}")
    if stretching_id is not None:
        parts.append(f"stretching_{stretching_id}")
    if wavelet_id is not None:
        parts.append(f"wavelet_{wavelet_id}")
    if wavelet_dtt_id is not None:
        parts.append(f"wavelet_dtt_{wavelet_dtt_id}")
    return resolve_lineage_params(session, parts)


# ============================================================
# Lineage enumeration helpers
# ============================================================


def get_next_lineage_batch(
        db,
        step_category,
        group_by="pair_lineage",
        loglevel="INFO",
        day_value=None,
):
    """
    Standard worker prolog for lineage-aware steps.

    - Claims the next batch of jobs for `step_category` using `get_next_job_for_step`.
    - Extracts (pair, lineage_str, refs, days).
    - Loads current step config (from the job row).
    - Resolves lineage_str -> WorkflowStep objects.
    - Merges params for that lineage + current step config.

    Returns
    -------
    dict with keys:
      - jobs, step
      - pair, lineage_str, lineage_steps
      - lineage_names          — full list including current step name
      - lineage_names_upstream — full minus current step (replaces manual [:-1])
      - lineage_names_mov      — upstream with any refstack_* entries stripped
                                 (used by mwcs/wct/stretching to find MOV CCFs)
      - refs, days
      - step_params, params
    or None if no jobs were claimed (caller should continue/sleep).
    """
    # Important: keep logging policy consistent with your scripts
    logger = get_logger(f"msnoise.worker.{step_category}", loglevel)

    # Ensure there is work (fast check)
    if not is_next_job_for_step(db, step_category=step_category, flag="T"):
        return None

    jobs, step = get_next_job_for_step(
        db,
        step_category=step_category,
        group_by=group_by,
        flag="T",
        # day_value=day_value,  # enable once you add day filtering to scheduler
    )

    if not jobs:
        return None

    pair = jobs[0].pair
    lineage_str = getattr(jobs[0], "lineage", None)
    if not lineage_str:
        raise ValueError(f"{step_category.upper()} jobs must have a non-empty lineage (v2 assumption)")

    refs = [job.ref for job in jobs]
    days = [job.day for job in jobs]

    # Load current step config from job row (job carries config_category & config_set_number)
    step_params = get_config_set_details(
        db,
        jobs[0].config_category,
        jobs[0].config_set_number,
        format="AttribDict",
    )

    # Merge params for THIS lineage only
    orig_params = get_params(db)
    lineage_steps = lineage_str_to_steps(db, lineage_str, strict=True)
    lineage_steps, lineage_names, params = get_merged_params_for_lineage(db, orig_params, step_params, lineage_steps)

    # Derived lineage lists — no manual slicing needed in callers
    lineage_names_upstream = lineage_names[:-1] if lineage_names else []
    lineage_names_mov = _strip_refstack_from_lineage(lineage_names_upstream)

    logger.info(f"New {step_category.upper()} batch: pair={pair} n={len(jobs)} group_by={group_by} lineage={lineage_str}")

    return {
        "jobs": jobs,
        "step": step,
        "pair": pair,
        "lineage_str": lineage_str,
        "lineage_steps": lineage_steps,
        "lineage_names": lineage_names,
        "lineage_names_upstream": lineage_names_upstream,
        "lineage_names_mov": lineage_names_mov,
        "refs": refs,
        "days": days,
        "step_params": step_params,
        "params": params,
    }


def get_stack_lineage_for_filter(session, filterid):
    """Get the full lineage path through a specific filter step to its downstream stack.

    Traverses the workflow graph upstream from ``filter_{filterid}`` to the root
    preprocess step, then appends the stack step that is immediately downstream of
    the filter.  Returns a list of step names suitable for use as the *lineage*
    argument to :func:`xr_get_ref`, :func:`xr_get_ccf`, etc.

    Example for the default single-pipeline::

        get_stack_lineage_for_filter(db, 1)
        # → ['preprocess_1', 'cc_1', 'filter_1', 'stack_1']

    :type session: :class:`sqlalchemy.orm.session.Session`
    :type filterid: int
    :param filterid: The filter set_number (e.g. 1 for filter_1).
    :rtype: list of str
    """
    steps = get_workflow_steps(session)
    links = get_workflow_links(session)
    step_map = {s.step_id: s for s in steps}

    filter_step = next(
        (s for s in steps if s.category == 'filter' and s.set_number == filterid),
        None,
    )
    if filter_step is None:
        return []

    # Walk upstream from filter to the root step, skipping global config steps
    parent_map = {link.to_step_id: link.from_step_id for link in links}
    path = [filter_step.step_name]
    current_id = filter_step.step_id
    visited = set()
    while current_id in parent_map:
        if current_id in visited:
            break  # guard against cycles
        visited.add(current_id)
        parent_id = parent_map[current_id]
        parent_step = step_map[parent_id]
        if parent_step.category != 'global':
            path.insert(0, parent_step.step_name)
        current_id = parent_id

    # Append the stack step immediately downstream of this filter (if any)
    for link in links:
        if link.from_step_id == filter_step.step_id:
            child = step_map.get(link.to_step_id)
            if child and child.category == 'stack':
                path.append(child.step_name)
                break

    return path


def get_refstack_lineage_for_filter(session, filterid, refstack_set_number=1):
    """Get the full lineage path through a filter step down to its refstack.

    Extends :func:`get_stack_lineage_for_filter` by also appending the
    ``refstack_M`` step that is immediately downstream of the stack step.
    This is the correct lineage to pass to :func:`xr_get_ref`, since REF
    files now live under the refstack step folder.

    Example for the default single-pipeline::

        get_refstack_lineage_for_filter(db, 1)
        # → ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']

    :type session: :class:`sqlalchemy.orm.session.Session`
    :type filterid: int
    :param filterid: The filter set_number (e.g. 1 for filter_1).
    :type refstack_set_number: int
    :param refstack_set_number: Which refstack set to use (default 1).
    :rtype: list of str
    """
    path = get_stack_lineage_for_filter(session, filterid)
    if not path:
        return path

    steps = get_workflow_steps(session)
    links = get_workflow_links(session)

    # Find the stack step at the end of the path
    stack_step_name = path[-1]
    stack_step = next(
        (s for s in steps if s.step_name == stack_step_name), None
    )
    if stack_step is None:
        return path

    # Find the refstack step downstream of this stack step
    # Prefer the requested refstack_set_number; fall back to the first available
    refstack_steps = [
        s for s in steps
        if s.category == 'refstack'
        and any(lk.from_step_id == stack_step.step_id and lk.to_step_id == s.step_id
                for lk in links)
    ]
    if not refstack_steps:
        return path  # no refstack in this workflow, return stack-level lineage

    # Pick requested set_number if available, else first
    refstack_step = next(
        (s for s in refstack_steps if s.set_number == refstack_set_number),
        refstack_steps[0]
    )
    return path + [refstack_step.step_name]


# ============================================================
# Section 7 — Time and Date Utilities
# ============================================================


def extend_days(days):
    """Return a :class:`~pandas.DatetimeIndex` from *days* extended by one
    extra day at the end.

    Replaces the pandas 1.x pattern::

        idx = pd.to_datetime(days)
        idx = idx.append(pd.DatetimeIndex([idx[-1] + pd.Timedelta("1d")]))

    which was removed in pandas 2.0.

    :param days: sequence of date-like values (strings, dates, datetimes…)
    :rtype: :class:`pandas.DatetimeIndex`
    """
    idx = pd.to_datetime(days)
    return pd.DatetimeIndex(list(idx) + [idx[-1] + pd.Timedelta("1d")])


def get_t_axis(params):
    """
    Returns the time axis (in seconds) of the CC functions.

    :rtype: :class:`numpy.array`
    :returns: the time axis in seconds
    """
    samples = int(2 * params.maxlag * params.cc_sampling_rate) + 1
    return np.linspace(-params.maxlag, params.maxlag, samples)


def get_maxlag_samples(maxlag, cc_sampling_rate):
    """
    Returns the length of the CC functions. Gets the maxlag and sampling rate
    from the database.


    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: int
    :returns: the length of the CCF in samples
    """

    return int(2*maxlag*cc_sampling_rate)+1


def build_ref_datelist(params):
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
    begin = params.ref_begin
    end = params.ref_end
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
    begin = get_config(session, "startdate", category='global', set_number=1)
    end = get_config(session, "enddate", category='global', set_number=1)
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


def refstack_is_rolling(params):
    """
    Return True if the refstack configset uses rolling-index mode.

    Rolling mode is indicated by ``ref_begin`` being a negative integer string
    (e.g. ``"-5"``). In this mode no REF file is written to disk; the reference
    is computed on-the-fly at MWCS/stretching/WCT time via
    :func:`compute_rolling_ref`.

    :type params: :class:`obspy.core.AttribDict`
    :param params: Merged parameter set containing ``ref_begin``.
    :rtype: bool
    """
    val = str(getattr(params, "ref_begin", "1970-01-01")).strip()
    return val.startswith("-")


def compute_rolling_ref(data, ref_begin, ref_end):
    """Compute a per-index rolling reference from CCF data.

    For each time index ``i``, the reference is::

        mean(data[i + ref_begin : i + ref_end])

    Both ``ref_begin`` and ``ref_end`` must be negative integers with
    ``ref_begin < ref_end <= -1``.  Use ``ref_end=-1`` to exclude the current
    window (compare to the previous N windows).

    Uses ``min_periods=1`` semantics: the first few steps receive whatever
    partial window is available rather than NaN.

    :type data: :class:`pandas.DataFrame` or :class:`xarray.DataArray`
        Shape ``(n_times, n_lag_samples)``.
    :type ref_begin: int
        Negative offset for the start of the rolling window (e.g. ``-5``).
    :type ref_end: int
        Negative offset for the end of the rolling window (e.g. ``-1``).
        Must satisfy ``ref_begin < ref_end <= 0``.
    :rtype: :class:`numpy.ndarray`
        Shape ``(n_times, n_lag_samples)``.  Row ``i`` is the reference for
        time step ``i``.
    """
    import xarray as xr_mod
    if isinstance(data, xr_mod.DataArray):
        arr = data.values
    else:
        arr = data.values  # DataFrame.values
    n = arr.shape[0]
    refs = np.full_like(arr, np.nan, dtype=float)

    for i in range(n):
        start_idx = i + ref_begin   # e.g. i - 5
        end_idx   = i + ref_end     # e.g. i - 1  (exclusive in slice)
        if end_idx <= 0:
            continue                # not enough history yet
        start_idx = max(0, start_idx)
        if start_idx >= end_idx:
            continue
        refs[i] = np.nanmean(arr[start_idx:end_idx], axis=0)

    return refs


def _strip_refstack_from_lineage(lineage_names):
    """
    Return a copy of ``lineage_names`` with any ``refstack_*`` entries removed.

    Used internally by :func:`get_next_lineage_batch` to build
    ``lineage_names_mov`` — the path for :func:`xr_get_ccf`, which lives
    under the ``stack_N`` step folder rather than under ``refstack_M``.

    Example::

        _strip_refstack_from_lineage(
            ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']
        )
        # → ['preprocess_1', 'cc_1', 'filter_1', 'stack_1']

    :type lineage_names: list of str
    :rtype: list of str
    """
    return [n for n in lineage_names if not n.startswith("refstack_")]


# Public alias kept for external code that imported the old name
strip_refstack_from_lineage = _strip_refstack_from_lineage


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


def filter_within_daterange(date, start_date, end_date):
    """Check if a date falls within the configured range"""
    return start_date <= date <= end_date


# ============================================================================
# Workflow Management Functions
# ============================================================================


# ============================================================
# Section 8 — Signal Processing Utilities
# ============================================================


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


def winsorizing(data, params, input="timeseries", nfft=0):
    """Clip (Winsorise) a 2-D data array in the time or frequency domain.

    Supports both one-shot sign-clipping (``winsorizing == -1``) and
    RMS-based clipping (``winsorizing > 0``).  When *input* is ``"fft"``
    the array is temporarily transformed back to the time domain, clipped,
    then re-transformed.

    :param data: 1-D or 2-D array of shape ``(n_traces, n_samples)``.
    :param params: MSNoise params object; must expose ``params.winsorizing``.
    :param input: ``"timeseries"`` (default) or ``"fft"``.
    :param nfft: FFT length used when *input* is ``"fft"``; ignored otherwise.
    :returns: Clipped array (same shape as input).
    """
    import scipy.fft as sf
    input1D = False
    if len(data.shape) == 1:
        data = data.reshape(-1, data.shape[0])
        input1D = True
    if input == "fft":
        data = sf.ifftn(data, [nfft, ], axes=[1, ]).astype(float)
    for i in range(data.shape[0]):
        if params.winsorizing == -1:
            np.sign(data[i], data[i])  # inplace
        elif params.winsorizing != 0:
            rms = data[i].std() * params.winsorizing
            np.clip(data[i], -rms, rms, data[i])  # inplace
    if input == "fft":
        data = sf.fftn(data, [nfft, ], axes=[1, ])
    if input1D:
        data = data[0]
    return data


def get_window(window="boxcar", half_win=3):
    """Return a normalised complex smoothing window for MWCS processing.

    :param window: ``"boxcar"`` (default) or ``"hanning"``.
    :param half_win: Half-width in samples (full window = ``2*half_win+1``).
    :returns: Complex numpy array of length ``2*half_win+1``, sum-normalised.
    """
    import scipy.signal
    window_len = 2 * half_win + 1
    if window == "boxcar":
        w = scipy.signal.windows.boxcar(window_len).astype("complex")
    else:
        w = scipy.signal.windows.hann(window_len).astype("complex")
    return w / window_len


def getCoherence(dcs, ds1, ds2):
    """Compute cross-coherence between two spectra.

    :param dcs: Cross-spectrum magnitudes (1-D array, length *n*).
    :param ds1: Auto-spectrum of signal 1 (1-D array, length *n*).
    :param ds2: Auto-spectrum of signal 2 (1-D array, length *n*).
    :returns: Complex coherence array of length *n*, clipped to ``|coh| <= 1``.
    """
    n = len(dcs)
    coh = np.zeros(n).astype("complex")
    valids = np.argwhere(np.logical_and(np.abs(ds1) > 0, np.abs(ds2) > 0))
    coh[valids] = dcs[valids] / (ds1[valids] * ds2[valids])
    coh[coh > (1.0 + 0j)] = 1.0 + 0j
    return coh


def prepare_abs_positive_fft(line, sampling_rate):
    """Return the positive-frequency part of the absolute FFT of *line*.

    :param line: 1-D signal array.
    :param sampling_rate: Sampling rate in Hz.
    :returns: ``(freq, val)`` - positive-frequency vector and absolute FFT values.
    """
    val = np.fft.fft(line)
    val = np.abs(val)
    freq = np.fft.fftfreq(len(line), 1.0 / sampling_rate)
    idx = np.where(freq >= 0)
    return freq[idx], val[idx]


def wavg_wstd(data, errors):
    """Weighted average and weighted standard deviation of equal-length arrays.

    :param data: 1-D array of measurement values.
    :param errors: 1-D array of measurement errors (zeros replaced by ``1e-6``).
    :returns: ``(weighted_mean, weighted_std)`` as floats.
    """
    errors = np.where(errors == 0, 1e-6, errors)
    w = 1.0 / errors
    wavg = (data * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (data - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wavg, wstd


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


def _get_wavgwstd(data, dttname, errname):
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


def _trim(data, dttname, limits=0.1):
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


# ============================================================
# Section 9 — Wavelet Coherence Transform (WCT)
# ============================================================


def _conv2(x, y, mode="same"):
    """2-D convolution using :func:`scipy.signal.convolve2d` with 180-degree rotations."""
    from scipy.signal import convolve2d
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)


def smoothCFS(cfs, scales, dt, ns, nt):
    """Smooth CWT coefficients along both time and scale axes.

    :param cfs: CWT coefficient array, shape ``(n_scales, n_times)``.
    :param scales: 1-D array of wavelet scales.
    :param dt: Sampling interval in seconds.
    :param ns: Length of the moving-average filter across scales.
    :param nt: Gaussian width parameter along time.
    :returns: Smoothed coefficient array, same shape as *cfs*.
    """
    import scipy.fft as sf
    N = cfs.shape[1]
    npad = sf.next_fast_len(N, real=True)
    omega = np.arange(1, np.fix(npad / 2) + 1, 1).tolist()
    omega = np.array(omega) * ((2 * np.pi) / npad)
    omega_save = -omega[int(np.fix((npad - 1) / 2)) - 1:0:-1]
    omega_2 = np.concatenate((0., omega), axis=None)
    omega_2 = np.concatenate((omega_2, omega_save), axis=None)
    omega = np.concatenate((omega_2, -omega[0]), axis=None)
    normscales = scales / dt
    for kk in range(0, cfs.shape[0]):
        F = np.exp(-nt * (normscales[kk] ** 2) * omega ** 2)
        smooth = np.fft.ifft(F * np.fft.fft(cfs[kk - 1], npad))
        cfs[kk - 1] = smooth[0:N]
    H = 1 / ns * np.ones((ns, 1))
    cfs = _conv2(cfs, H)
    return cfs


def get_wavelet_type(wavelet_type):
    """Return a :mod:`pycwt` wavelet object for the given type/parameter pair.

    :param wavelet_type: Tuple ``(name, param)`` or ``(name,)``.
        Supported names: ``"Morlet"``, ``"Paul"``, ``"DOG"``, ``"MexicanHat"``.
    :returns: Corresponding :mod:`pycwt` wavelet instance.
    """
    import pycwt as wavelet
    defaults = {"Morlet": 6, "Paul": 4, "DOG": 2, "MexicanHat": 2}
    name = wavelet_type[0]
    param = float(wavelet_type[1]) if len(wavelet_type) == 2 else defaults[name]
    if name == "Morlet":
        return wavelet.Morlet(param)
    elif name == "Paul":
        return wavelet.Paul(param)
    elif name == "DOG":
        return wavelet.DOG(param)
    elif name == "MexicanHat":
        return wavelet.MexicanHat()
    else:
        raise ValueError(f"Unknown wavelet type: {name!r}")


def xwt(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12,
         freqmin=0.1, freqmax=8.0, nptsfreq=100, wavelet_type=("Morlet", 6.)):
    """Wavelet Coherence Transform (WCT) between two time series.

    :param trace_ref: Reference signal (1-D array).
    :param trace_current: Current signal (1-D array, same length).
    :param fs: Sampling frequency in Hz.
    :param ns: Scale-axis smoothing parameter.
    :param nt: Time-axis smoothing parameter.
    :param vpo: Voices-per-octave; higher = finer scale resolution.
    :param freqmin: Lowest frequency of interest (Hz).
    :param freqmax: Highest frequency of interest (Hz).
    :param nptsfreq: Number of frequency points between *freqmin* and *freqmax*.
    :param wavelet_type: ``(name, param)`` tuple passed to :func:`get_wavelet_type`.
    :returns: ``(WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi)``
    """
    import pycwt as wavelet
    mother = get_wavelet_type(wavelet_type)
    nx = np.size(trace_current)
    x_reference = np.transpose(trace_ref)
    x_current = np.transpose(trace_current)
    dt = 1 / fs
    dj = 1 / vpo
    J = -1
    s0 = 2 * dt
    freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True)
    cwt_reference, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(
        x_reference, dt, dj, s0, J, mother, freqs=freqlim)
    cwt_current, _, _, _, _, _ = wavelet.cwt(
        x_current, dt, dj, s0, J, mother, freqs=freqlim)
    scales = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1 / scales)
    power_ref = (invscales * abs(cwt_reference) ** 2).astype(complex)
    power_cur = (invscales * abs(cwt_current) ** 2).astype(complex)
    crossCFS = cwt_reference * np.conj(cwt_current)
    WXamp = abs(crossCFS)
    cross_spectrum = (invscales * crossCFS).astype(complex)
    cfs1 = smoothCFS(power_ref, scales, dt, ns, nt)
    cfs2 = smoothCFS(power_cur, scales, dt, ns, nt)
    crossCFS = smoothCFS(cross_spectrum, scales, dt, ns, nt)
    mask1 = cfs1 > 0
    mask2 = cfs2 > 0
    valid_mask = mask1 & mask2
    WXspec = np.full_like(crossCFS, np.nan, dtype=complex)
    Wcoh = np.full_like(crossCFS, np.nan)
    WXspec[valid_mask] = crossCFS[valid_mask] / (
        np.sqrt(cfs1[valid_mask]) * np.sqrt(cfs2[valid_mask]))
    Wcoh[valid_mask] = (abs(crossCFS[valid_mask]) ** 2
                        / (cfs1[valid_mask] * cfs2[valid_mask]))
    WXangle = np.angle(WXspec)
    Wcoh = np.clip(Wcoh, 0.0, 1.0)
    pp = 2 * np.pi * freqs
    pp2 = np.array([[kk] for kk in pp])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)
    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi


def compute_wct_dtt(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20,
                    mincoh=0.5, maxdt=0.2, min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Compute dv/v and associated errors from wavelet coherence transform results.

    :param freqs: Frequency values from the WCT.
    :param tvec: Time axis.
    :param WXamp: Cross-wavelet amplitude array (freqs × taxis).
    :param Wcoh: Wavelet coherence array (freqs × taxis).
    :param delta_t: Time delay array (freqs × taxis).
    :param lag_min: Minimum coda lag in seconds.
    :param coda_cycles: Number of periods to use as coda window width.
    :param mincoh: Minimum coherence threshold.
    :param maxdt: Maximum allowed time delay.
    :param min_nonzero: Minimum fraction of valid (non-zero weight) samples required.
    :param freqmin: Lower frequency bound for regression.
    :param freqmax: Upper frequency bound for regression.
    :returns: Tuple of (dt/t, err, weighting_function).
    """
    import warnings
    from scipy.optimize import OptimizeWarning
    from obspy.signal.regression import linear_regression

    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
    dvv = np.zeros(len(inx[0]))
    err = np.zeros(len(inx[0]))

    weight_func = np.log(np.abs(WXamp)) / np.log(np.abs(WXamp)).max()
    zero_idx = np.where((Wcoh < mincoh) | (delta_t > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / weight_func.max()
    wf[zero_idx] = 0

    problematic_freqs = []
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)
        tindex = np.where(
            ((tvec >= -lag_max) & (tvec <= -lag_min)) |
            ((tvec >= lag_min) & (tvec <= lag_max))
        )[0]
        if len(tvec) > 2:
            if not np.any(delta_t[ifreq]):
                continue
            delta_t[ifreq][tindex] = np.nan_to_num(delta_t[ifreq][tindex])
            w = wf[ifreq]
            nzc_perc = np.count_nonzero(w[tindex] > 0) / len(tindex)
            if nzc_perc >= min_nonzero:
                with warnings.catch_warnings(record=True) as w_catcher:
                    warnings.simplefilter("always", OptimizeWarning)
                    m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex],
                                              w[tindex], intercept_origin=True)
                    if any(issubclass(warning.category, OptimizeWarning)
                           for warning in w_catcher):
                        problematic_freqs.append(freqs[ifreq])
                dvv[ii], err[ii] = m, em
            else:
                dvv[ii], err[ii] = np.nan, np.nan
    if problematic_freqs:
        logging.warning(
            f"Covariance issues at {min(problematic_freqs):.2f}-{max(problematic_freqs):.2f} Hz: "
            f"consider adjusting min_nonzero={min_nonzero}, mincoh={mincoh}, "
            f"maxdt={maxdt}, coda_cycles={coda_cycles}"
        )
    return dvv, err, wf


def get_wct_avgcoh(freqs, tvec, wcoh, freqmin, freqmax, lag_min=5, coda_cycles=20):
    """
    Calculate average wavelet coherence over a frequency range and coda window.

    :param freqs: Frequency array.
    :param tvec: Time axis.
    :param wcoh: Wavelet coherence array (freqs × taxis).
    :param freqmin: Lower frequency bound.
    :param freqmax: Upper frequency bound.
    :param lag_min: Minimum coda lag in seconds.
    :param coda_cycles: Number of periods to use as coda window width.
    :returns: Average coherence per frequency bin within [freqmin, freqmax].
    """
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
    coh = np.zeros(inx[0].shape)
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)
        tindex = np.where(
            ((tvec >= -lag_max) & (tvec <= -lag_min)) |
            ((tvec >= lag_min) & (tvec <= lag_max))
        )[0]
        if len(tvec) > 2:
            if not np.any(wcoh[ifreq]) or wcoh[ifreq][tindex].size == 0:
                coh[ii] = np.nan
                continue
            coh[ii] = np.nanmean(np.abs(wcoh[ifreq][tindex]))
        else:
            coh[ii] = np.nan
    return coh


# ============================================================
# Section 10 — Preprocessing I/O
# ============================================================


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


def save_preprocessed_streams(stream, output_dir, step_name, goal_day):
    """Write a preprocessed Stream to disk in workflow-aware path structure.

    Creates ``<output_dir>/<step_name>/_output/`` and writes
    ``<goal_day>.mseed``.

    :param stream: :class:`~obspy.core.stream.Stream` to write.
    :param output_dir: Base output directory.
    :param step_name: Workflow step name (e.g. ``"preprocess_1"``).
    :param goal_day: Processing date string (``YYYY-MM-DD``).
    :returns: List containing the single saved file path.
    """
    workflow_dir = os.path.join(output_dir, step_name, "_output")
    os.makedirs(workflow_dir, exist_ok=True)
    output_path = os.path.join(workflow_dir, f"{goal_day}.mseed")
    for tr in stream:
        tr.data = tr.data.astype(np.float32)
    stream.write(output_path, format="MSEED")
    return [output_path]


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


# ============================================================
# Section 11 — xarray I/O
# ============================================================
# ── Core helpers ───────────────────────────────────────────


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
    elif name == "STR":
        keys = ["Delta", "Coeff", "Error"]
        data = np.random.random((len(times), len(keys)))
        dr = xr.DataArray(data, coords=[times, keys],
                          dims=["times", "keys"])
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

# ── CCF ─────────────────────────────────────────────────────


def xr_save_ccf(root, lineage, step_name, station1, station2, components, mov_stack, taxis, new, overwrite=False):
    path = os.path.join(root, *lineage, step_name, "_output",
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


def xr_get_ccf(root, lineage, station1, station2, components, mov_stack, taxis, format="dataset"):
    """Load CCF results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataset"`` (default) returns a :class:`xarray.DataArray` (CCF variable).
        ``"dataframe"`` returns a :class:`~pandas.DataFrame` (legacy).
        ``"xarray"`` is an alias for ``"dataset"`` (deprecated, kept for compat).
    """
    path = os.path.join(root, *lineage, "_output",
                        "%s_%s" % (mov_stack[0], mov_stack[1]), "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)

    fullpath = os.path.join(path, fn)
    if not os.path.isfile(fullpath):
        raise FileNotFoundError(fullpath)
    data = xr_create_or_open(fullpath, taxis, name="CCF")
    if format in ("dataset", "xarray"):
        return data.CCF
    # ── DataFrame (legacy) ──────────────────────────────────────────────
    return data.CCF.to_dataframe().unstack().droplevel(0, axis=1)


def xr_save_ref(root, lineage, step_name, station1, station2, components, taxis, new, overwrite=False):
    path = os.path.join(root, *lineage, step_name, "_output",
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


def xr_get_ref(root, lineage, station1, station2, components, taxis, ignore_network=False):
    path = os.path.join(root, *lineage, "_output",
                        "REF", "%s" % components)
    # If ignore_network is True, strip the network code from the station names
    if ignore_network:
        s1_parts = station1.split('.')
        s2_parts = station2.split('.')

        available_files = glob.glob(os.path.join(path, "*.%s.%s_*.%s.%s.nc" % (s1_parts[1],s1_parts[2], s2_parts[1], s1_parts[2])))

        if available_files:
            # Use the first available reference file
            fullpath = available_files[0]
        else:
            raise FileNotFoundError(f"No reference file found for station {s1_parts[1]} and {s2_parts[1]}")
    else:
        fn = "%s_%s.nc" % (station1, station2)

        fullpath = os.path.join(path, fn)
        if not os.path.isfile(fullpath):
            # logging.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
            raise FileNotFoundError(fullpath)
    data = xr_create_or_open(fullpath, taxis, name="REF")
    return data


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


def get_results_all(session, root, lineage_names, station1, station2, components, dates,
                    format="dataframe", params=None):
    """
    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type station1: str
    :param station1: The name of station 1 (formatted NET.STA.LOC)
    :type station2: str
    :param station2: The name of station 2 (formatted NET.STA.LOC)
    :type components: str
    :param components: The name of the components used (ZZ, ZR, ...)
    :type dates: list
    :param dates: List of TODO datetime.datetime
    :rtype: :class:`pandas.DataFrame`
    :return: All CCF results in a :class:`pandas.DataFrame`, where the index
        is the time of the CCF and the columns are the times in the coda.
    """

    path = os.path.join(root, *lineage_names, "_output", "all", components,  station1, station2)
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

    if len(results) > 0:
        result = pd.concat(results)
        del results
        if format == "dataframe":
            return result
        elif format == "xarray":
            taxis = get_t_axis(params)
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

# ── MWCS ────────────────────────────────────────────────────


def xr_save_mwcs(root, lineage, step_name, station1, station2, components, mov_stack, taxis, dataset):
    """Save MWCS results to a NetCDF file.

    :param dataset: :class:`xarray.Dataset` with a ``MWCS`` variable and
        dims ``(times, taxis, keys)``, as built by
        :mod:`~msnoise.s05_compute_mwcs`.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                        "%s_%s" % (mov_stack[0], mov_stack[1]),
                       "%s" % components,
                       "%s_%s.nc" % (station1, station2))
    os.makedirs(os.path.split(fn)[0], exist_ok=True)
    dr = xr_create_or_open(fn, taxis=dataset.coords.get("taxis", taxis), name="MWCS")
    rr = xr_insert_or_update(dr, dataset)
    xr_save_and_close(rr, fn)


def xr_get_mwcs(root, lineage, station1, station2, components, mov_stack, format="dataset"):
    """Load MWCS results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        with MultiIndex columns ``(keys, taxis)``.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset`.
    """
    fn = os.path.join(root, *lineage, "_output",
                        "%s_%s" % (mov_stack[0], mov_stack[1]),
                       "%s" % components,
                       "%s_%s.nc" % (station1, station2))
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    data = xr_create_or_open(fn, name="MWCS")

    if format == "dataset":
        return data

    # ── DataFrame (legacy) ──────────────────────────────────────
    da = data.MWCS  # DataArray with dims (times, taxis, keys)
    # Build DataFrame directly — pandas-version-safe
    times_vals = da.coords["times"].values
    taxis_vals = da.coords["taxis"].values
    keys_vals  = da.coords["keys"].values
    n_t, n_tx, n_k = da.values.shape
    # Transpose to (times, keys, taxis) so MultiIndex is (keys, taxis)
    midx = pd.MultiIndex.from_product([keys_vals, taxis_vals], names=["keys", "taxis"])
    arr = da.values.transpose(0, 2, 1).reshape(n_t, n_k * n_tx)
    return pd.DataFrame(arr, index=pd.DatetimeIndex(times_vals), columns=midx)

# ── DTT ─────────────────────────────────────────────────────


def xr_save_dtt(root, lineage, step_name, station1, station2, components, mov_stack, dataset):
    """Save DTT results to a NetCDF file.

    :param dataset: :class:`xarray.Dataset` with a ``DTT`` variable and
        dims ``(times, keys)``, as built by :mod:`~msnoise.s06_compute_mwcs_dtt`.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))
    os.makedirs(os.path.split(fn)[0], exist_ok=True)
    dr = xr_create_or_open(fn, taxis=[], name="DTT")
    rr = xr_insert_or_update(dr, dataset)
    xr_save_and_close(rr, fn)


def xr_get_dtt(root, lineage, station1, station2, components, mov_stack, format="dataset"):
    """Load DTT results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset`.
    """
    fn = os.path.join(root, *lineage, "_output",
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))

    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    dr = xr_create_or_open(fn, taxis=[], name="DTT")

    if format == "dataset":
        return dr

    # ── DataFrame (legacy) ──────────────────────────────────────────────
    da = dr.DTT  # DataArray with dims (times, keys)
    return pd.DataFrame(
        da.values,
        index=pd.DatetimeIndex(da.coords["times"].values),
        columns=list(da.coords["keys"].values),
    )

# ── Stretching ──────────────────────────────────────────────


def xr_save_stretching(root, lineage, step_name, station1, station2,
                        components, mov_stack, dataset):
    """Save per-pair stretching results to a NetCDF file.

    Path layout::

        <root>/<lineage>/<step_name>/_output/<mov_stack[0]>_<mov_stack[1]>/<components>/<sta1>_<sta2>.nc

    :param dataset: :class:`xarray.Dataset` with a ``STR`` variable and
        dims ``(times, keys)`` where keys = ``['Delta', 'Coeff', 'Error']``,
        as built by :mod:`~msnoise.s10_stretching`.
    """
    fn = os.path.join(
        root, *lineage, step_name, "_output",
        "%s_%s" % (mov_stack[0], mov_stack[1]),
        components,
        "%s_%s.nc" % (station1, station2),
    )
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    dr = xr_create_or_open(fn, taxis=[], name="STR")
    rr = xr_insert_or_update(dr, dataset)
    xr_save_and_close(rr, fn)


def _xr_get_stretching(root, lineage, station1, station2, components, mov_stack, format="dataset"):
    """Load per-pair stretching results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        with columns ``Delta``, ``Coeff``, ``Error``.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset`.

    :raises FileNotFoundError: if the NetCDF file does not exist.
    """
    fn = os.path.join(
        root, *lineage, "_output",
        "%s_%s" % (mov_stack[0], mov_stack[1]),
        components,
        "%s_%s.nc" % (station1, station2),
    )
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    dr = xr_create_or_open(fn, taxis=[], name="STR")

    if format == "dataset":
        return dr

    # ── DataFrame (legacy) ──────────────────────────────────────────────
    da = dr.STR  # DataArray with dims (times, keys)
    return pd.DataFrame(
        da.values,
        index=pd.DatetimeIndex(da.coords["times"].values),
        columns=list(da.coords["keys"].values),
    )


# ============================================================
# Stretching aggregate (all pairs → stats DataFrame)
# ============================================================

# ── DVV ─────────────────────────────────────────────────────





# ── WCT ─────────────────────────────────────────────────────


def xr_save_wct(root, lineage, step_name, station1, station2, components, mov_stack, taxis, freqs, WXamp_list, WXcoh_list, WXdt_list, dates_list):
    """
    Save WCT results into an xarray Dataset and store it as a NetCDF file.

    Parameters:
    - root: str, Root output folder path
    - lineage: list, Lineage path components
    - step_name: str, Step name for output path
    - station1, station2: str, Station pair
    - components: str, Seismic component (e.g., ZZ)
    - mov_stack: tuple, Moving stack window (e.g., ('1d', '1d'))
    - taxis, freqs: np.array, Time axis and frequency axis
    - WXamp_list, WXcoh_list, WXdt_list: list of np.array, WCT outputs
    - dates_list: list of datetime, Timestamps for each WCT calculation
    """

    # Convert lists to xarray DataArrays (all WCT outputs are real-valued;
    # cast explicitly in case intermediate complex dtype was inherited)
    WXamp_da = xr.DataArray(
        data=np.array(WXamp_list).real,
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="WXamp"
    )

    Wcoh_da = xr.DataArray(
        data=np.array(WXcoh_list).real,
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="Wcoh"
    )

    WXdt_da = xr.DataArray(
        data=np.array(WXdt_list).real,
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="WXdt"
    )

    # Combine into an xarray Dataset
    ds = xr.Dataset({"WXamp": WXamp_da, "Wcoh": Wcoh_da, "WXdt": WXdt_da})

    # Define output directory
    fn = os.path.join(root, *lineage, step_name, "_output",
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      components, f"{station1}_{station2}.nc")

    os.makedirs(os.path.dirname(fn), exist_ok=True)

    # Save to NetCDF
    xr_save_and_close(ds, fn)


    # Cleanup memory
    del ds, WXamp_da, Wcoh_da, WXdt_da


def xr_load_wct(root, lineage, station1, station2, components, mov_stack):
    """
    Load WCT results from an xarray Dataset stored in a NetCDF file.

    Parameters:
    - root: str, Root output folder path
    - lineage: list, Lineage path components
    - station1, station2: str, Station pair
    - components: str, Seismic component (e.g., ZZ)
    - mov_stack: tuple, Moving stack window (e.g., ('1d', '1d'))

    Returns:
    - ds: xarray.Dataset containing the WCT data (WXamp, Wcoh, WXdt)
    """

    # Construct the file path
    fn = os.path.join(root, *lineage, "_output",
                      f"{mov_stack[0]}_{mov_stack[1]}", components,
                      f"{station1}_{station2}.nc")

    # Check if the file exists
    if not os.path.exists(fn):
        raise FileNotFoundError(f"File not found: {fn}")

    # Load and return the dataset
    ds = xr.load_dataset(fn)
    return ds


def xr_save_wct_dtt(root, lineage, step_name, station1, station2, components, mov_stack, taxis, dataset):
    """Save WCT-DTT results to a NetCDF file.

    :param dataset: :class:`xarray.Dataset` with variables ``dtt``, ``err``,
        ``coh`` and dims ``(times, frequency)``, as built by
        :mod:`~msnoise.s09_compute_wct_dtt`.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                      f"{mov_stack[0]}_{mov_stack[1]}", components,
                      f"{station1}_{station2}.nc")
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    existing_ds = xr_create_or_open(fn, name="WCT")
    updated_ds = xr_insert_or_update(existing_ds, dataset)
    xr_save_and_close(updated_ds, fn)
    logging.debug(f"Saved WCT DTT data to {fn}")


def xr_get_wct_dtt(root, lineage, station1, station2, components, mov_stack):
    """Load per-pair WCT dt/t results from a NetCDF file.

    Returns an :class:`xarray.Dataset` with variables ``dtt``, ``err``,
    ``coh`` each of shape ``(times, frequency)``, as written by
    :func:`xr_save_wct_dtt`.

    :raises FileNotFoundError: if the NetCDF file does not exist.
    :rtype: :class:`xarray.Dataset`
    """
    fn = os.path.join(
        root, *lineage, "_output",
        "%s_%s" % (mov_stack[0], mov_stack[1]),
        components,
        "%s_%s.nc" % (station1, station2),
    )
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    return xr.load_dataset(fn)


# ============================================================
# WCT-DTT aggregate (all pairs → frequency-band stats DataFrame)
# ============================================================


# ============================================================
# Section 12 — CCF and Stack Computation
# ============================================================


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


def export_allcorr(session, ccfid, data, base_folder=None, params=None, t_axis=None):
    if base_folder is None:
        # Legacy behaviour:
        output_folder = get_config(session, 'output_folder')
        base_folder = output_folder
    else:
        # If params.output_folder is an absolute project root you want to respect,
        # you can anchor relative base_folder under it:
        if params is not None and hasattr(params, "output_folder") and not os.path.isabs(base_folder):
            base_folder = os.path.join(params.output_folder, base_folder)

    station1, station2, components, filterid, date = ccfid.split('+')

    path = os.path.join(
        base_folder,
        filterid,
        '_output',
        'all',
        components,
        station1, station2
    )
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

    df = pd.DataFrame().from_dict(data).T
    df.columns = t_axis
    df.to_hdf(os.path.join(path, date+'.h5'), key='data')
    del df
    return


def add_corr(session, station1, station2, filterid, date, time, duration,
             components, CF, sampling_rate, day=False, ncorr=0, params=None, base_folder=None):
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
        if base_folder is None:
            # Legacy behaviour:
            path = os.path.join(
                "STACKS", "%02i" % filterid, "001_DAYS", components,
                          "%s_%s" % (station1, station2), str(date)
            )
        else:
            # Workflow-aware behaviour:
            if params is not None and hasattr(params, "output_folder") and not os.path.isabs(base_folder):
                base_folder = os.path.join(params.output_folder, base_folder)

            path = os.path.join(
                base_folder,
                filterid,
                '_output',
                'daily',
                components,
                station1, station2,
                str(date),
            )
        pair = "%s:%s" % (station1, station2)
        if mseed:
            _export_mseed(session, path, pair, components, filterid, CF,
                         ncorr, params=params)
        if sac:
            _export_sac(session, path, pair, components, filterid, CF,
                       ncorr, params=params)

    else:
        file = '%s.cc' % time
        path = os.path.join(output_folder, filterid, station1,
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


def _export_sac(db, filename, pair, components, filterid, corr, ncorr=0,
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


def _export_mseed(db, filename, pair, components, filterid, corr, ncorr=0,
                 maxlag=None, cc_sampling_rate=None, params=None):
    from obspy import Trace, Stream
    try:
        os.makedirs(os.path.split(filename)[0], exist_ok=True)
    except:
        pass
    filename += ".MSEED"
    maxlag = params.maxlag
    cc_sampling_rate = params.cc_sampling_rate
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






# ============================================================
# Section 13b — DVV Aggregate I/O + Aggregation Helpers
# ============================================================


def _classify_pair_type(sta1: str, sta2: str, component: str = "") -> str:
    """Classify a station pair as CC, SC, or AC.

    Classification is based on the SEED ``NET.STA.LOC`` ids and the component
    pair string (e.g. ``"ZZ"``, ``"ZN"``):

    - **AC** (autocorrelation): ``sta1 == sta2`` *and* both component letters
      are identical (e.g. ``"ZZ"``, ``"EE"``).  Same instrument correlating
      with itself.
    - **SC** (single-channel cross-component): ``sta1 == sta2`` *and*
      component letters differ (e.g. ``"ZN"``).  Different channels at the
      same location.
    - **CC** (cross-correlation): ``sta1 != sta2``.  Different physical
      locations regardless of component.

    :param sta1: First station SEED id ``NET.STA.LOC``.
    :param sta2: Second station SEED id ``NET.STA.LOC``.
    :param component: Component-pair string (e.g. ``"ZZ"``, ``"ZN"``).
        If empty or length < 2, falls back to comparing only the SEED ids.
    :returns: ``"AC"``, ``"SC"``, or ``"CC"``.
    """
    if sta1 != sta2:
        return "CC"
    # Same NET.STA.LOC — distinguish by component
    if len(component) >= 2 and component[0] == component[-1]:
        return "AC"
    return "SC"


def _dvv_column_spec(parent_category: str, pair_type: str, params) -> tuple:
    """Return ``(dv_col, err_col, quality_col)`` for the given parent category.

    For ``mwcs_dtt``, the columns depend on the pair type and user config
    (smart defaults: CC→free-intercept m/em, SC/AC→forced-zero m0/em0).
    For ``stretching`` and ``wct_dtt`` the columns are fixed.

    :param parent_category: ``"mwcs_dtt"``, ``"stretching"``, or ``"wavelet_dtt"``.
    :param pair_type: ``"CC"``, ``"SC"``, ``"AC"``, or ``"ALL"``.
    :param params: Params object from :func:`get_params`.
    :returns: Tuple ``(dv_col, err_col, quality_col)`` where any of the
        quality columns may be ``None`` if not available.
    """
    if parent_category == "mwcs_dtt":
        pt = pair_type if pair_type in ("CC", "SC", "AC") else "CC"
        pt_lower = pt.lower()
        dv_col  = getattr(params, f"dvv_{pt_lower}_value",  "m" if pt == "CC" else "m0")
        err_col = getattr(params, f"dvv_{pt_lower}_error", "em" if pt == "CC" else "em0")
        return dv_col, err_col, "mcoh"
    elif parent_category == "stretching":
        return "Delta", "Error", "Coeff"
    elif parent_category == "wavelet_dtt":
        # WCT stores dtt/err/coh per frequency; scalar extracted separately
        return "dtt", "err", "coh"
    else:
        raise ValueError(f"Unknown parent_category: {parent_category!r}")


def _freq_average_wct(ds, freqmin: float, freqmax: float,
                       quality_min: float, freq_agg: str = "mean"):
    """Collapse the ``(times, frequency)`` WCT-DTT Dataset to scalar time series.

    Selects the frequency band ``[freqmin, freqmax]``, masks cells where
    ``coh < quality_min``, then collapses the frequency axis using *freq_agg*.

    :param ds: :class:`xarray.Dataset` with variables ``dtt``, ``err``, ``coh``
        and dims ``(times, frequency)``.
    :param freqmin: Lower frequency bound (Hz).
    :param freqmax: Upper frequency bound (Hz).
    :param quality_min: Minimum coherence; cells below are set to NaN.
    :param freq_agg: ``"mean"`` or ``"median"`` over the frequency axis.
    :returns: Tuple ``(dv_da, err_da)`` — two 1-D :class:`xarray.DataArray`
        with dim ``times``.
    :raises ValueError: if no frequencies fall within [freqmin, freqmax].
    """
    freqs = ds.coords["frequency"].values
    mask_freq = (freqs >= freqmin) & (freqs <= freqmax)
    if not mask_freq.any():
        raise ValueError(
            f"No frequencies in [{freqmin}, {freqmax}] Hz in WCT-DTT output. "
            f"Available: {freqs.min():.3f}–{freqs.max():.3f} Hz"
        )
    dtt_sel = ds["dtt"].isel(frequency=mask_freq)
    err_sel = ds["err"].isel(frequency=mask_freq)
    coh_sel = ds["coh"].isel(frequency=mask_freq)

    # Mask low-coherence cells
    bad = coh_sel < quality_min
    dtt_sel = dtt_sel.where(~bad)
    err_sel = err_sel.where(~bad)

    # Collapse frequency axis
    if freq_agg == "median":
        dv_da  = dtt_sel.median(dim="frequency")
        err_da = err_sel.median(dim="frequency")
    else:
        dv_da  = dtt_sel.mean(dim="frequency")
        err_da = err_sel.mean(dim="frequency")

    return dv_da, err_da


def xr_save_dvv_agg(root, lineage, step_name, mov_stack,
                    pair_type: str, component: str, dataset):
    """Save a DVV aggregate result to a NetCDF file.

    Path layout::

        <root>/<lineage>/<step_name>/_output/<mov_stack>/dvv_<pair_type>_<component>.nc

    :param dataset: :class:`xarray.Dataset` with dim ``times`` and stat
        variables as built by :func:`aggregate_dvv_pairs`.
    """
    ms_str = "%s_%s" % (mov_stack[0], mov_stack[1])
    fn = os.path.join(root, *lineage, step_name, "_output",
                      ms_str, f"dvv_{pair_type}_{component}.nc")
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    dataset.to_netcdf(fn, mode="w")


def xr_get_dvv_agg(root, lineage, step_name, mov_stack,
                   pair_type: str, component: str, format: str = "dataset"):
    """Load a DVV aggregate result from a NetCDF file.

    :param format: ``"dataset"`` (default) or ``"dataframe"``.
    :raises FileNotFoundError: if the file does not exist.
    """
    ms_str = "%s_%s" % (mov_stack[0], mov_stack[1])
    fn = os.path.join(root, *lineage, step_name, "_output",
                      ms_str, f"dvv_{pair_type}_{component}.nc")
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    ds = xr.load_dataset(fn)
    if format == "dataframe":
        return ds.to_dataframe()
    return ds


def aggregate_dvv_pairs(root, parent_lineage, parent_step_name,
                        parent_category: str, mov_stack, component: str,
                        pair_type: str, pairs, params) -> xr.Dataset:
    """Aggregate per-pair DTT/STR/WCT-DTT results into network-level dv/v statistics.

    Reads all per-pair output files for the given ``(mov_stack, component)``
    combination, extracts the appropriate dv/v and error columns per pair type,
    applies quality filtering, then computes across-pair statistics at each
    time step.

    :param root: Output folder root.
    :param parent_lineage: Lineage name list up to and including the parent
        DTT step (e.g. ``["preprocess_1", ..., "mwcs_dtt_1"]``).
    :param parent_step_name: Step name of the parent DTT step.
    :param parent_category: ``"mwcs_dtt"``, ``"stretching"``, or ``"wavelet_dtt"``.
    :param mov_stack: Tuple ``(window, step)`` e.g. ``("1D", "1D")``.
    :param component: Component string e.g. ``"ZZ"``.
    :param pair_type: ``"CC"``, ``"SC"``, ``"AC"``, or ``"ALL"``.
    :param pairs: Iterable of ``(sta1, sta2)`` SEED-id string tuples.
    :param params: Params object from :func:`get_params`.
    :returns: :class:`xarray.Dataset` with dim ``times`` and stat variables.
    :raises ValueError: if no data files are found for the given combination.
    """
    dv_col, err_col, quality_col = _dvv_column_spec(
        parent_category, pair_type, params)

    quality_min    = float(getattr(params, "dvv_quality_min", 0.0))
    do_weighted    = str(getattr(params, "dvv_weighted_mean", "Y")).upper() == "Y"
    do_trimmed     = str(getattr(params, "dvv_trimmed_mean", "Y")).upper() == "Y"
    trim_sigma     = float(getattr(params, "dvv_trim_limit", 3.0))
    out_percent    = str(getattr(params, "dvv_output_percent", "Y")).upper() == "Y"
    percentile_str = getattr(params, "dvv_percentiles", "5,25,75,95")
    percentiles    = [float(p) for p in str(percentile_str).split(",") if p.strip()]

    # WCT-specific params
    if parent_category == "wavelet_dtt":
        wct_freqmin  = float(getattr(params, "dvv_freqmin", 0.1))
        wct_freqmax  = float(getattr(params, "dvv_freqmax", 2.0))
        wct_freq_agg = str(getattr(params, "dvv_freq_agg", "mean"))

    # ── 1. Collect per-pair 1-D time series ─────────────────────────────
    pair_dvv  = []   # list of (times_array, dv_array, err_array)

    for sta1, sta2 in pairs:
        # Filter by pair_type
        pt = _classify_pair_type(sta1, sta2, component)
        if pair_type != "ALL" and pt != pair_type:
            continue

        try:
            if parent_category == "mwcs_dtt":
                ds = xr_get_dtt(root, parent_lineage, sta1, sta2,
                                 component, mov_stack, format="dataset")
                da_dv  = ds["DTT"].sel(keys=dv_col)
                da_err = ds["DTT"].sel(keys=err_col)
                if quality_col in ds["DTT"].coords["keys"].values:
                    da_q = ds["DTT"].sel(keys=quality_col)
                    bad = da_q < quality_min
                    da_dv  = da_dv.where(~bad)
                    da_err = da_err.where(~bad)

            elif parent_category == "stretching":
                ds = _xr_get_stretching(root, parent_lineage, sta1, sta2,
                                         component, mov_stack, format="dataset")
                da_dv  = ds["STR"].sel(keys=dv_col)
                da_err = ds["STR"].sel(keys=err_col)
                if quality_col:
                    da_q = ds["STR"].sel(keys=quality_col)
                    bad = da_q < quality_min
                    da_dv  = da_dv.where(~bad)
                    da_err = da_err.where(~bad)
                # Convert Delta → dv/v: dv/v = Delta - 1
                da_dv = da_dv - 1.0

            elif parent_category == "wavelet_dtt":
                ds = xr_get_wct_dtt(root, parent_lineage, sta1, sta2,
                                     component, mov_stack)
                da_dv, da_err = _freq_average_wct(
                    ds, wct_freqmin, wct_freqmax, quality_min, wct_freq_agg)

        except FileNotFoundError:
            continue
        except Exception:
            continue

        if out_percent:
            da_dv = da_dv * 100.0

        times = da_dv.coords["times"].values
        dv    = da_dv.values.astype(float)
        err   = da_err.values.astype(float)
        pair_dvv.append((times, dv, err))

    if not pair_dvv:
        raise ValueError(
            f"No data found for parent={parent_category} "
            f"step={parent_step_name} mov_stack={mov_stack} "
            f"component={component} pair_type={pair_type}"
        )

    # ── 2. Build (pairs × times) arrays on a common time axis ───────────
    all_times_sorted = np.array(sorted({t for times, _, _ in pair_dvv
                                         for t in times}),
                                dtype="datetime64[ns]")
    n_t = len(all_times_sorted)
    n_p = len(pair_dvv)

    dv_mat  = np.full((n_p, n_t), np.nan)
    err_mat = np.full((n_p, n_t), np.nan)

    time_idx = {t: i for i, t in enumerate(all_times_sorted)}
    for p, (times, dv, err) in enumerate(pair_dvv):
        for ti, t in enumerate(times):
            key = np.datetime64(t, "ns")
            if key in time_idx:
                j = time_idx[key]
                dv_mat[p, j]  = dv[ti]
                err_mat[p, j] = err[ti]

    # ── 3. Compute statistics across pairs at each time step ─────────────
    # n_pairs: count of non-NaN pairs per time step
    n_pairs = np.sum(~np.isnan(dv_mat), axis=0).astype(float)

    mean_dv   = np.nanmean(dv_mat, axis=0)
    std_dv    = np.nanstd(dv_mat, axis=0, ddof=1)
    median_dv = np.nanmedian(dv_mat, axis=0)

    data_vars = {
        "mean":     ("times", mean_dv),
        "std":      ("times", std_dv),
        "median":   ("times", median_dv),
        "n_pairs":  ("times", n_pairs),
    }

    # Percentiles
    for pct in percentiles:
        key = f"q{int(pct):02d}"
        pct_vals = np.nanpercentile(dv_mat, pct, axis=0)
        data_vars[key] = ("times", pct_vals)

    # Weighted mean/std
    if do_weighted:
        w_mat = np.where(err_mat > 0, 1.0 / err_mat**2, 0.0)
        w_sum = np.nansum(w_mat, axis=0)
        w_sum_safe = np.where(w_sum > 0, w_sum, np.nan)
        wmean = np.nansum(w_mat * np.nan_to_num(dv_mat), axis=0) / w_sum_safe
        # weighted std
        wvar  = np.nansum(w_mat * (np.nan_to_num(dv_mat) - wmean[None, :])**2,
                           axis=0) / w_sum_safe
        wstd  = np.sqrt(wvar)
        data_vars["weighted_mean"] = ("times", wmean)
        data_vars["weighted_std"]  = ("times", wstd)

    # Trimmed mean/std (sigma-based: remove pairs > trim_sigma*std from mean)
    if do_trimmed:
        tmean = np.full(n_t, np.nan)
        tstd  = np.full(n_t, np.nan)
        for j in range(n_t):
            col = dv_mat[:, j]
            valid = col[~np.isnan(col)]
            if len(valid) < 2:
                tmean[j] = np.nan if len(valid) == 0 else valid[0]
                tstd[j]  = np.nan
                continue
            mu = np.mean(valid)
            sigma = np.std(valid, ddof=1)
            kept = valid[np.abs(valid - mu) <= trim_sigma * sigma]
            tmean[j] = np.mean(kept) if len(kept) > 0 else np.nan
            tstd[j]  = np.std(kept, ddof=1) if len(kept) > 1 else np.nan
        data_vars["trimmed_mean"] = ("times", tmean)
        data_vars["trimmed_std"]  = ("times", tstd)

    # ── 4. Build output Dataset ──────────────────────────────────────────
    ds_out = xr.Dataset(
        {k: xr.DataArray(v[1], dims=["times"],
                          coords={"times": all_times_sorted})
         for k, v in data_vars.items()},
        attrs={
            "parent_category": parent_category,
            "parent_step":     parent_step_name,
            "pair_type":       pair_type,
            "component":       component,
            "mov_stack":       "%s_%s" % (mov_stack[0], mov_stack[1]),
            "dvv_unit":        "percent" if out_percent else "fraction",
            "created":         str(np.datetime64("now")),
        }
    )
    return ds_out


# ============================================================
# Section 14 — PSD and HDF
# ============================================================


def psd_dfrms(a):
    """Integrate a PSD Series over its period axis using the trapezoid rule.

    :param a: :class:`pandas.Series` whose index is period values.
    :returns: Square-root of the integrated power (RMS-equivalent).
    """
    return np.sqrt(np.trapezoid(a.values, a.index))


def psd_rms(s, f):
    """Compute RMS from a power spectrum array and frequency array.

    :param s: Power spectral density values (1-D array).
    :param f: Frequency values (1-D array, same length as *s*).
    :returns: Float - square-root of the integrated power.
    """
    return np.sqrt(np.trapezoid(s, f))


def psd_df_rms(d, freqs, output="VEL"):
    """Compute per-frequency-band RMS from a PPSD DataFrame.

    :param d: :class:`pandas.DataFrame` with period columns and time index.
    :param freqs: List of ``(fmin, fmax)`` tuples defining frequency bands.
    :param output: Physical unit - ``"VEL"`` (default), ``"ACC"``, or ``"DISP"``.
    :returns: :class:`pandas.DataFrame` with one column per frequency band.
    """
    d = d.dropna(axis=1, how="all")
    RMS = {}
    for fmin, fmax in freqs:
        pmin = 1.0 / fmax
        pmax = 1.0 / fmin
        ix = np.where((d.columns >= pmin) & (d.columns <= pmax))[0]
        spec = d.iloc[:, ix]
        f = d.columns[ix]
        w2f = 2.0 * np.pi * f
        amp = 10.0 ** (spec / 10.0)
        if output == "ACC":
            RMS[f"{fmin:.1f}-{fmax:.1f}"] = amp.apply(psd_dfrms, axis=1)
        elif output == "VEL":
            vamp = amp / w2f ** 2
            RMS[f"{fmin:.1f}-{fmax:.1f}"] = vamp.apply(psd_dfrms, axis=1)
        else:
            vamp = amp / w2f ** 2
            damp = vamp / w2f ** 2
            RMS[f"{fmin:.1f}-{fmax:.1f}"] = damp.apply(psd_dfrms, axis=1)
    return pd.DataFrame(RMS, index=d.index)


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


def psd_ppsd_to_dataset(ppsd):
    """Convert an ObsPy PPSD object to an :class:`xarray.Dataset`.

    Builds the same data as :func:`psd_ppsd_to_dataframe` but returns an
    ``xr.Dataset`` with a ``PSD`` variable of dims ``(times, periods)`` —
    ready to pass directly to :func:`xr_save_psd` without a DataFrame
    round-trip.

    :param ppsd: :class:`~obspy.signal.spectral_estimation.PPSD` object.
    :returns: :class:`xarray.Dataset`.
    """
    from obspy import UTCDateTime
    times = np.array([UTCDateTime(t).datetime for t in ppsd.current_times_used])
    periods = np.asarray(ppsd.period_bin_centers, dtype=float)
    data = np.asarray(ppsd._binned_psds, dtype=float)
    da = xr.DataArray(
        data,
        coords=[times, periods],
        dims=["times", "periods"],
        name="PSD",
    )
    return da.to_dataset()


def xr_save_psd(root, lineage, step_name, seed_id, day, dataset):
    """Save a daily PSD result to a NetCDF file.

    Path layout::

        <root>/<lineage>/<step_name>/_output/daily/<seed_id>/<YYYY-MM-DD>.nc

    :param dataset: :class:`xarray.Dataset` with a ``PSD`` variable and
        dims ``(times, periods)``, as built by :func:`psd_ppsd_to_dataset`.
    """
    day_str = day if isinstance(day, str) else day.strftime("%Y-%m-%d")
    fn = os.path.join(
        root, *lineage, step_name, "_output", "daily",
        seed_id, f"{day_str}.nc",
    )
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    dataset.to_netcdf(fn, mode="w")


def xr_load_psd(root, lineage, step_name, seed_id, day, format="dataset"):
    """Load a daily PSD NetCDF written by :func:`xr_save_psd`.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        or ``None`` if file not found.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset` or ``None``.
    """
    day_str = day if isinstance(day, str) else day.strftime("%Y-%m-%d")
    fn = os.path.join(
        root, *lineage, step_name, "_output", "daily",
        seed_id, f"{day_str}.nc",
    )
    if not os.path.isfile(fn):
        return None
    ds = xr.load_dataset(fn)

    if format == "dataset":
        return ds

    # ── DataFrame (legacy) ──────────────────────────────────────
    da = ds.PSD
    df = pd.DataFrame(
        da.values,
        index=pd.DatetimeIndex(da.coords["times"].values),
        columns=da.coords["periods"].values.astype(float),
    )
    ds.close()
    return df

def xr_save_rms(root, lineage, step_name, seed_id, dataframe):
    """Save per-station PSD RMS results to a NetCDF file.

    Accepts either a :class:`~pandas.DataFrame` (legacy, index = DatetimeIndex,
    columns = band labels) or an :class:`xarray.Dataset` with a ``RMS``
    variable and dims ``(times, bands)``.
    """
    fn = os.path.join(root, *lineage, step_name, "_output", seed_id, "RMS.nc")
    os.makedirs(os.path.dirname(fn), exist_ok=True)

    # ── Dataset path (new, xarray-native) ──────────────────────────────
    if isinstance(dataframe, xr.Dataset):
        if os.path.isfile(fn):
            existing_ds = xr.load_dataset(fn)
            merged = xr_insert_or_update(existing_ds, dataframe)
            merged.to_netcdf(fn, mode="w")
        else:
            dataframe.to_netcdf(fn, mode="w")
        return

    # ── DataFrame path (legacy compat) ─────────────────────────────────
    bands = list(dataframe.columns.astype(str))
    times = pd.DatetimeIndex(dataframe.index)

    if os.path.isfile(fn):
        existing = xr_load_rms(root, lineage, step_name, seed_id,
                               format="dataframe")
        if existing is not None:
            existing = existing[~existing.index.isin(times)]
            dataframe = pd.concat([existing, dataframe]).sort_index()
            bands = list(dataframe.columns.astype(str))
            times = pd.DatetimeIndex(dataframe.index)

    da = xr.DataArray(
        dataframe.values.astype(float),
        coords=[times, bands],
        dims=["times", "bands"],
        name="RMS",
    )
    da.to_dataset().to_netcdf(fn, mode="w")


def xr_load_rms(root, lineage, step_name, seed_id, format="dataset"):
    """Load per-station PSD RMS results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        or ``None`` if file not found.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset` or ``None``.
    """
    fn = os.path.join(root, *lineage, step_name, "_output", seed_id, "RMS.nc")
    if not os.path.isfile(fn):
        return None
    ds = xr.load_dataset(fn)

    if format == "dataset":
        return ds

    # ── DataFrame (legacy) ──────────────────────────────────────
    df = pd.DataFrame(
        ds.RMS.values,
        index=pd.DatetimeIndex(ds.RMS.coords["times"].values),
        columns=list(ds.RMS.coords["bands"].values),
    )
    ds.close()
    return df


# ============================================================
# Section 15 — Filters (legacy)
# ============================================================



