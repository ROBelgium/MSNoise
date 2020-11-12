import collections
import os
try:
    import cPickle
except:
    import pickle as cPickle

from logbook import Logger, StreamHandler
import sys

from sqlalchemy import create_engine, func
from sqlalchemy.orm import sessionmaker
from sqlalchemy.pool import NullPool
from sqlalchemy.sql.expression import func


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
    station = session.query(Station).filter(Station.net == net). \
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
