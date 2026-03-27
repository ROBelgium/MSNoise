import traceback
import logging
import os
import sys

# import importlib_metadata
import sqlalchemy
from sqlalchemy import text
import time

import click
# from click_plugins import with_plugins
# import importlib.metadata

try:
    from .._version import version as __version__
except ImportError:
    __version__ = "0.0.0.dev0"


from .. import MSNoiseError, DBConfigNotFoundError
from ....db import connect, get_logger
from ....config import (delete_config_set, get_config, get_config_set_details, list_config_sets)
from ....stations import update_station

from ..msnoise_table_def import DataAvailability


class OrderedGroup(click.Group):
    def list_commands(self, ctx):
        return self.commands.keys()


def run_threaded(main_func, ctx, *args, **kwargs):
    """Run *main_func* once (threads=1) or in parallel worker processes.

    Reads ``MSNOISE_threads``, ``MSNOISE_threadsdelay`` and
    ``MSNOISE_verbosity`` from *ctx.obj* so callers don't have to repeat
    that boilerplate.  Any positional *args* and extra keyword arguments
    are forwarded to *main_func* alongside ``loglevel``.
    """
    threads  = ctx.obj['MSNOISE_threads']
    delay    = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main_func(*args, loglevel=loglevel, **kwargs)
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main_func, args=args,
                        kwargs={"loglevel": loglevel, **kwargs})
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


def parse_extra_args(ctx, param, extra_args):
    # extra_args = extra_args.split(" ")
    kwargs = {}
    for e in extra_args:
        kw, value = e.split("=")
        kwargs[kw.replace("--", "")] = eval(value)
    return kwargs


def validate_verbosity(ctx, param, value):
    """
    Validate the --quiet and --verbose options to have it conflict with one
    another.
    """
    if param.name == 'quiet':
        excluded = 'verbose'
    elif param.name == 'verbose':
        excluded = 'quiet'
    if value and excluded in ctx.params:
        raise click.BadParameter('Cannot use both --quiet and --verbose option.')
    return value


def show_config_values(db, names):
    """
    Show configuration value of parameters provided in the 'names' list.
    """
    from ..config import get_config
    from ..default import default
    for key in names:
        display_value = value = get_config(db, key)
        if value == '':
            # Use a more explicit representation of the empty string
            display_value = "''"
        try:
            default_value = default[key].default
        except KeyError:
            click.secho("Error: unknown parameter '%s'" % key)
            continue
        if value == default_value:
            click.secho("   %s: %s" % (key, display_value))
        else:
            click.secho(" M %s: %s" % (key, display_value), fg='green')


def info_db_ini():
    """
    Show information stored in the db.ini file.
    """
    from ..db import read_db_inifile
    dbini = read_db_inifile()
    click.echo('Database information stored in the db.ini file:')
    if dbini.tech == 1:
        click.echo(' - database type: SQLite')
        click.echo(' - filename: {}'.format(dbini.hostname))
    elif dbini.tech == 2:
        click.echo(' - database type: MySQL')
        click.echo(' - hostname: {}'.format(dbini.hostname))
        click.echo(' - database: {}'.format(dbini.database))
        click.echo(' - username: {}'.format(dbini.username))
        click.echo(' - password: {}'.format('*' * len(dbini.password)))
        click.echo(' - table prefix: {}'.format(dbini.prefix if dbini.prefix else '(none)'))


def info_folders(db):
    """
    Show information about folders used by MSNoise.
    """
    from ..config import get_config
    click.echo('')
    click.echo('General:')

    click.echo(' - the MSNoise module is installed in: %s'
               % os.path.abspath(os.path.join(
                   os.path.dirname(__file__), '..', '..')))

    if os.path.isfile('db.ini'):
        click.echo(' - db.ini is present in the current directory')
    else:
        click.secho(' - db.ini not found in the present directory, is MSNoise installed here ?',
                    fg='red')
        return
    info_db_ini()

    click.echo('')
    click.echo('Configuration:')

    data_folder = get_config(db, "data_folder")
    if os.path.isdir(data_folder):
        click.echo(" - %s exists" % data_folder)
    else:
        click.secho(" - %s does not exists !" % data_folder, fg='red')

    output_folder = get_config(db, "output_folder")
    if os.path.isdir(output_folder):
        click.echo(" - %s exists" % output_folder)
    else:
        if get_config(db, 'keep_all') in ['Y', 'y']:
            from ..workflow import get_workflow_job_counts
            counts = get_workflow_job_counts(db)
            done_count = counts.get('D', 0)
            if done_count > 0:
                click.secho(
                    " - %s does not exists and that is not normal"
                    " (%i jobs done)" % (output_folder, done_count),
                    fg='red')
            else:
                click.secho(
                    " - %s does not exists and that is normal"
                    " (%i jobs done)" % (output_folder, done_count))
        else:
            click.secho(
                " - %s does not exists (and that is normal because"
                " keep_all=False)" % output_folder)


def info_parameters(db):
    """
    Show values of each configuration parameters.
    """
    from ..default import default
    from ..config import get_config_set_details
    from ..workflow import get_workflow_steps
    click.echo('')
    click.echo('Configuration values:'
            '   | Normal colour indicates that the default value is used'
            '   | Green indicates "M"odified values')
    # TODO: add plugins params
    show_config_values(db, default.keys())

    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    if filter_steps:
        click.echo('')
        click.echo('Filter configsets:')
        for step in filter_steps:
            details = get_config_set_details(db, 'filter', step.set_number)
            click.echo(f'  filter_{step.set_number} ({step.step_name}):')
            for row in details:
                click.echo(f'    {row["name"]} = {row["value"]}')


def info_stations(db):
    """
    Show information about configured stations.
    """
    from ..stations import get_stations
    click.echo('')
    click.echo('Stations:')
    click.echo('  NET.STA    Long.     Lat.    Alt.   Coord. Used?')
    s = None
    na_sign = '-'
    for s in get_stations(db, all=True):
        click.echo(' {:>8s}  {:>9s} {:>8s}  {:>6s}  {:3s}    {:1s}'.format(
                '{}.{}'.format(s.net, s.sta),
                '{:.4f}'.format(s.X) if s.X is not None else na_sign,
                '{:.4f}'.format(s.Y) if s.Y is not None else na_sign,
                '{:.1f}'.format(s.altitude) if s.altitude is not None else na_sign,
                s.coordinates or na_sign,
                'Y' if s.used else 'N'))
        click.echo("  | Location code(s): %s" % s.used_location_codes)
        click.echo("  | Channel names(s): %s" % s.used_channel_names)
    if s is None:
        click.echo(' ')


def info_jobs(db):
    """
    Show information about jobs registered in database.
    """
    from ..workflow import get_workflow_job_counts
    counts = get_workflow_job_counts(db)
    click.echo("Jobs:")
    for flag, label in [('T', 'Todo'), ('I', 'In Progress'), ('D', 'Done'), ('F', 'Failed')]:
        click.echo('  %s: %i' % (label, counts.get(flag, 0)))


def info_plugins(db):
    """
    Show information about configured plugins.
    """
    plugins = get_config(db, "plugins")
    if not plugins:
        return
    plugins = plugins.split(",")
    from importlib.metadata import entry_points
    for ep in list(entry_points(group='msnoise.plugins.jobtypes')):
        module_name = ep.value.split(".")[0]
        if module_name in plugins:
            click.echo('')
            click.echo('Plugin: %s' % module_name)
            for row in ep.load()():
                click.echo(' %s' % row["name"])


# if sys.version_info < (3, 11):
#     click_command_tree_entry_points = importlib.metadata.entry_points().get('click_command_tree', [])
# else:
#     click_command_tree_entry_points = importlib.metadata.entry_points.select(group='click_command_tree')

# @with_plugins(click_command_tree_entry_points)
@click.group(context_settings=dict(max_content_width=120), cls=OrderedGroup)
@click.option('-t', '--threads', default=1, help='Number of threads to use \
(only affects modules that are designed to do parallel processing)')
@click.option('-d', '--delay', default=1,  help='In the case of multi-threading'
                    ', defines the number of seconds to wait before lauching '
                    'the next thread. Defaults to [1] second ')
@click.option('-c', '--custom', default=False, is_flag=True, help='Use custom \
 file for plots. To use this, copy the plot script here and edit it.')
@click.option('-v', '--verbose', is_flag=True, callback=validate_verbosity)
@click.option('-q', '--quiet', is_flag=True, default=False,
              callback=validate_verbosity)
@click.version_option(__version__)
@click.pass_context
def cli(ctx, threads, delay, custom, verbose, quiet):
    ctx.obj['MSNOISE_threads'] = threads
    ctx.obj['MSNOISE_threadsdelay'] = delay
    ctx.obj['MSNOISE_custom'] = custom
    ctx.obj['MSNOISE_verbosity'] = "INFO"
    if quiet:
        ctx.obj['MSNOISE_verbosity'] = "WARNING"
    elif verbose:
        ctx.obj['MSNOISE_verbosity'] = "DEBUG"
    global logger
    if threads > 1:
        logger = get_logger('msnoise', ctx.obj['MSNOISE_verbosity'], with_pid=True)
    else:
        logger = get_logger('msnoise', ctx.obj['MSNOISE_verbosity'])
    # Is this really needed?
    if custom:
        sys.path.append(os.getcwd())



@cli.command()
@click.option('-p', '--port', default=5000, help='Port to open')
def admin(port):
    """Starts the Web Admin on http://localhost:5000 by default"""
    from ..msnoise_admin import main
    main(port)

@cli.group(cls=OrderedGroup)
def db():
    """Commands to interact with the database"""
    pass


@db.command(name="init")
@click.option('--tech', help='Database technology: 1=SQLite 2=MySQL/MariaDB 3=PostgreSQL',
              default=None)
@click.option('--auto-workflow', is_flag=True, default=False,
              help='Automatically create all default config sets, workflow steps and links '
                   'without prompting.')
def db_init(tech, auto_workflow):
    """This command initializes the current folder to be a MSNoise Project
    by creating a database and a db.ini file."""
    click.echo('Launching the init')
    from ..s00_installer import main
    result = main(tech)
    if result != 0:
        return

    if auto_workflow or click.confirm(
        '\nWould you like to automatically create all default config sets, '
        'workflow steps and links?',
        default=True
    ):
        _create_default_workflow()


def _create_default_workflow():
    """Create default config sets, workflow steps and links after db init."""
    from ..db import connect
    from ..config import create_config_set
    from ..workflow import (create_workflow_steps_from_config_sets,
                            create_workflow_links_from_steps)

    ALL_CATEGORIES = [
        'preprocess', 'cc', 'filter', 'stack',
        'refstack',
        'mwcs', 'mwcs_dtt', 'mwcs_dtt_dvv',
        'stretching', 'stretching_dvv',
        'wavelet', 'wavelet_dtt', 'wavelet_dtt_dvv',
        'psd', 'psd_rms',
    ]

    db = connect()

    # 1. Config sets
    click.echo('\n[1/3] Creating default config sets...')
    created_sets = []
    for cat in ALL_CATEGORIES:
        set_number = create_config_set(db, cat)
        if set_number is not None:
            created_sets.append(f'{cat}_{set_number}')
    db.commit()
    if created_sets:
        click.echo(f'      Created: {", ".join(created_sets)}')
    else:
        click.echo('      All config sets already exist.')

    # 2. Workflow steps
    click.echo('[2/3] Creating workflow steps from config sets...')
    created, existing, err = create_workflow_steps_from_config_sets(db)
    if err:
        click.echo(f'      Error: {err}', err=True)
    else:
        click.echo(f'      Created {created} step(s), {existing} already existed.')

    # 3. Workflow links
    click.echo('[3/3] Creating workflow links between steps...')
    created, existing, err = create_workflow_links_from_steps(db)
    if err:
        click.echo(f'      Error: {err}', err=True)
    else:
        click.echo(f'      Created {created} link(s), {existing} already existed.')

    db.close()
    click.echo('\nSetup complete! Your workflow is ready to use.')


@db.command(name="update_loc_chan")
@click.pass_context
def db_da_stations_update_loc_chan(ctx):
    """Populates the Location & Channel from the Data Availability
    table. Warning: rewrites automatically, no confirmation."""
    from msnoise.db import connect
    from msnoise.stations import get_stations
    session = connect()
    stations = get_stations(session)
    for sta in stations:
        data = session.query(DataAvailability). \
            filter(text("net=:net")). \
            filter(text("sta=:sta")). \
            group_by(DataAvailability.net, DataAvailability.sta,
                     DataAvailability.loc, DataAvailability.chan). \
            params(net=sta.net, sta=sta.sta).all()
        locids = list(set(sorted([d.loc for d in data])))
        chans = list(set(sorted([d.chan for d in data])))
        # logger.info("%s.%s has locids:%s and chans:%s" % (sta.net, sta.sta,
        #                                                   locids, chans))
        sta.used_location_codes = ",".join(locids)
        sta.used_channel_names = ",".join(chans)
        try:
            session.commit()
        except Exception:
            traceback.print_exc()

@db.command(name="execute")
@click.argument('sql_command')
@click.option('-o', '--outfile', help='Output filename (?="request.csv")',
              default=None, type=str)
@click.option('-s', '--show', help='Show output (in case of SELECT statement)?',
              default=True, type=bool)
@click.pass_context
def db_execute(ctx, sql_command, outfile=None, show=True):
    """EXPERT MODE: Executes 'sql_command' on the database. Use this command
    at your own risk!!"""
    from msnoise.db import connect
    db = connect()
    for cmd in sql_command.split(";"):
        if not len(cmd):
            continue
        logger.info("Executing '%s'" % cmd)
        r = db.execute(text(cmd))
        if cmd.count("select") or cmd.count("SELECT"):
            result = r.fetchall()
            if not len(result):
                logger.info("The query returned no results, sorry.")
            else:
                import pandas as pd
                df = pd.DataFrame(result, columns=r.keys())
                if show:
                    pd.set_option('display.max_rows', None)
                    pd.set_option('display.max_columns', None)
                    pd.set_option('display.width', None)
                    pd.set_option('display.max_colwidth', None)
                    print(df)
                if outfile:
                    if outfile == "?":
                        df.to_csv("request.csv")
                    else:
                        df.to_csv("%s" % outfile)
    db.commit()
    db.close()



@db.command(name="upgrade")
def db_upgrade():
    """Upgrade the database from previous to a new version.
    This procedure adds new parameters with their default value
    in the config database.
    """
    from ..db import connect, read_db_inifile
    from ..default import default
    db = connect()
    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    for name in default.keys():
        try:
            db.add(Config(name=name, value=default[name].default))
            db.commit()
        except Exception:
            db.rollback()
            # print("Passing %s: already in DB" % name)
            continue
    try:
        db.execute(text("CREATE UNIQUE INDEX job_index ON %sjobs (day, pair, "
                   "jobtype)" %
                   prefix))
        db.commit()
    except Exception:
        logging.info("It looks like the v1.5 'job_index' is already in the DB")
        db.rollback()

    try:
        db.execute(text("CREATE INDEX job_index2 ON %sjobs (jobtype, flag)" %
                   prefix))
        db.commit()
    except Exception:
        logging.info("It looks like the v1.6 'job_index2' is already in the DB")
        db.rollback()

    try:
        db.execute(text("CREATE UNIQUE INDEX da_index ON %sdata_availability (path, "
                   "file, net, sta, loc, chan)" %
                   prefix))
        db.commit()
    except Exception:
        logging.info("It looks like the v1.5 'da_index' is already in the DB")
        db.rollback()

    db.close()


@db.command(name='clean_duplicates')
def db_clean_duplicates():
    """Checks the Jobs table and deletes duplicate entries"""
    from msnoise.db import connect, read_db_inifile

    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    db = connect()
    if dbini.tech == 1:
        # SQlite case
        query = "WITH cte AS (SELECT *, ROW_NUMBER() OVER (PARTITION BY day, pair, jobtype ORDER BY ref) AS rn FROM {0}jobs) DELETE FROM {0}jobs WHERE ref IN (SELECT ref FROM cte WHERE rn > 1);".format(prefix)
    elif dbini.tech == 2:
        # Mysql/mariadb
        query = "DELETE j1 FROM {0}jobs j1 INNER JOIN {0}jobs j2 ON j1.ref > j2.ref AND j1.day = j2.day AND j1.pair = j2.pair AND j1.jobtype = j2.jobtype;".format(prefix)
    else:
        # postgresql
        query = "DELETE FROM {0}jobs WHERE ref NOT IN (SELECT MIN(ref) FROM {0}jobs GROUP BY day, pair, jobtype);".format(prefix)
    db.execute(text(query))
    db.commit()
    db.close()


@db.command(name="dump")
@click.option("--format", default="csv")
def db_dump(format):
    """Dumps the complete database in formatted files, defaults to CSV.
    """
    from ..db import get_engine
    from sqlalchemy import MetaData
    import pandas as pd

    if format == "csv":
        engine = get_engine(inifile=os.path.join(os.getcwd(), 'db.ini'))

        meta = MetaData()
        meta.reflect(bind=engine)

        for table in meta.sorted_tables:
            with engine.connect() as conn:
                rows = conn.execute(table.select()).all()
            df = pd.DataFrame(rows)
            click.echo("Dumping table %s to %s.csv" % (table.name, table.name))
            df.to_csv("%s.csv" % table.name, index=False)
    else:
        click.echo("Currently only the csv format is supported, sorry.", err=True)


@db.command(name="import")
@click.argument("table")
@click.option("--format", default="csv")
@click.option("--force", is_flag=True, default=False)
def db_import(table, format, force):
    """
    Imports msnoise tables from formatted files (CSV).
    """
    from ..db import get_engine, read_db_inifile
    import pandas as pd
    dbini = read_db_inifile(inifile=os.path.join(os.getcwd(), 'db.ini'))

    if format == "csv":
        engine = get_engine(inifile=os.path.join(os.getcwd(), 'db.ini'))
        logger.info("Loading table %s from %s.csv" % (table, table))
        df = pd.read_csv("%s.csv" % table)
        if force:
            if_exists = "replace"
        else:
            if_exists = "fail"
        try:
            if dbini.tech == 1:
                with engine.connect() as conn:
                    df.to_sql(table, conn.connection, if_exists=if_exists)
            else:
                df.to_sql(table, engine, if_exists=if_exists)
        except ValueError:
            traceback.print_exc()
            logger.info("You're probably getting the error above because the "
                  "table already exists, if you want to replace the table "
                  "with the imported data, then pass the --force option")
    else:
        logger.error("Currently only the csv format is supported, sorry.")


@cli.command()
@click.option('-j', '--jobs', is_flag=True, help='Jobs Info only')
def info(jobs):
    """Outputs general information about the current install and config, plus
    information about jobs and their status."""
    from ..db import connect

    if not os.path.isfile('db.ini'):
        click.secho(' - db.ini is not present, is MSNoise installed here ?',
                    fg='red')
        return

    db = connect()
    if not jobs:
        info_folders(db)
        info_parameters(db)
        info_stations(db)
    info_jobs(db)
    info_plugins(db)
    db.close()


@cli.group()
def config():
    """
    This command allows to set a parameter value in the database, show
    parameter values, or synchronise station metadata, depending on the
    invoked subcommands.
    """
    pass


@config.command(name='sync')
def config_sync():
    """
    Synchronise station metadata from inventory/dataless.
    """
    from ..db import connect
    from ..stations import get_stations, update_station
    from ..signal import preload_instrument_responses

    db = connect()
    responses = preload_instrument_responses(db)
    netsta = []
    for id, row in responses.iterrows():
        net, sta, loc, chan = row["channel_id"].split(".")
        netsta.append("%s.%s"%(net,sta))
    responses["netsta"] = netsta

    for station in get_stations(db):
        id = "%s.%s" % (station.net, station.sta)
        coords = responses[responses["netsta"] == id]
        if not len(coords):
            logger.error("No coords for %s, skipping,..." % id)
            continue
        try:
            lon = float(coords["longitude"].values[0])
            lat = float(coords["latitude"].values[0])
            elevation = float(coords["elevation"].values[0])
        except Exception:
            logging.warning(
                'Problem getting coordinates for '
                '"%s": %s' % (id, str(coords)))
            continue
        update_station(db, net=station.net, sta=station.sta, X=lon, Y=lat,
                       altitude=elevation, coordinates="DEG")
        logging.info("Added coordinates (%.5f %.5f %.1f) for station %s.%s" %
                     (lon, lat, elevation, station.net, station.sta))
    db.close()


@config.command(name='set')
@click.argument('name_value')
def config_set(name_value):
    """
    Set a configuration value. The argument should be of the form
    'variable=value'.
    """
    from ..default import default
    if not name_value.count("="):
        click.echo("!! format of the set command is name=value !!")
        return
    name, value = name_value.split("=")
    if name not in default:
        click.echo("!! unknown parameter %s !!" % name)
        return
    from ..config import update_config
    from ..db import connect
    db = connect()
    update_config(db, name, value)
    db.commit()
    db.close()
    click.echo("Successfully updated parameter %s = %s" % (name, value))


@config.command(name='get')
@click.argument('names', nargs=-1)
def config_get(names):
    """
    Display the value of the given configuration variable(s).
    """
    from ..db import connect
    db = connect()
    show_config_values(db, names)
    db.close()


@config.command(name='reset')
@click.argument('names', nargs=-1)
def config_reset(names):
    """
    Reset the value of the given configuration variable(s) to their default.
    """
    from ..default import default
    from ..config import update_config
    from ..db import connect
    for key in names:
        default_value = default[key].default
        db = connect()
        update_config(db, key, default_value)
        # db.commit()
        db.close()
        click.echo("Successfully reset parameter %s = %s" % (key, default_value))

@config.command()
@click.argument('set_name')
@click.pass_context
def create_set(ctx, set_name):
    """Create a configuration set for a workflow step.

    SET_NAME: Name of the workflow step (e.g., mwcs, mwcs_dtt, etc.)
    """
    from ..config import create_config_set
    from ..db import connect

    db = connect()

    set_number = create_config_set(db, set_name)
    if set_number is not None:
        click.echo(f"Configuration set '{set_name}' created successfully with set_number: {set_number}")
    else:
        click.echo(f"Error: Configuration file for '{set_name}' not found.")

    db.close()


@config.command()
@click.argument('set_name')
@click.argument('set_number', type=int)
@click.option('--confirm', is_flag=True,
              help='Skip confirmation prompt')
def delete_set(set_name, set_number, confirm):
    """Delete a configuration set.

    SET_NAME: The category name (e.g., mwcs, mwcs_dtt)
    SET_NUMBER: The set number to delete
    """
    session = connect()

    try:
        # Show details of what will be deleted
        details = get_config_set_details(session, set_name, set_number)

        if not details:
            click.echo(f"Config set '{set_name}' set_number {set_number} not found.")
            return

        click.echo(f"Config set '{set_name}' set_number {set_number} contains {len(details)} entries:")
        for detail in details:
            click.echo(f"  - {detail['name']}: {detail['value']}")

        if not confirm:
            if not click.confirm("\nAre you sure you want to delete this config set?"):
                click.echo("Deletion cancelled.")
                return

        # Perform deletion
        if delete_config_set(session, set_name, set_number):
            click.echo(f"✓ Successfully deleted config set '{set_name}' set_number {set_number}")
        else:
            click.echo(f"✗ Failed to delete config set '{set_name}' set_number {set_number}")

    except Exception as e:
        click.echo(f"Error: {str(e)}")
    finally:
        session.close()


@config.command()
@click.option('--category', '-c', help='Filter by category name')
def list_sets(category):
    """List all configuration sets."""
    session = connect()

    try:
        sets = list_config_sets(session, category)

        if not sets:
            if category:
                click.echo(f"No config sets found for category '{category}'")
            else:
                click.echo("No config sets found")
            return

        click.echo("Configuration Sets:")
        click.echo("=" * 60)
        click.echo(f"{'Category':<20} {'Set #':<8} {'Entries':<8}")
        click.echo("-" * 60)

        for cat, set_num, count in sets:
            click.echo(f"{cat:<20} {set_num:<8} {count:<8}")

    except Exception as e:
        click.echo(f"Error: {str(e)}")
    finally:
        session.close()


@config.command()
@click.argument('set_name')
@click.argument('set_number', type=int)
def show_set(set_name, set_number):
    """Show details of a configuration set.

    SET_NAME: The category name (e.g., mwcs, mwcs_dtt)
    SET_NUMBER: The set number to show
    """
    session = connect()

    try:
        details = get_config_set_details(session, set_name, set_number)

        if not details:
            click.echo(f"Config set '{set_name}' set_number {set_number} not found.")
            return

        click.echo(f"Config Set: {set_name} #{set_number}")
        click.echo("=" * 60)

        for detail in details:
            click.echo(f"Name: {detail['name']}")
            click.echo(f"Value: {detail['value']}")
            if detail['description']:
                click.echo(f"Description: {detail['description']}")
            click.echo("-" * 40)

    except Exception as e:
        click.echo(f"Error: {str(e)}")
    finally:
        session.close()


@config.command()
@click.argument('old_set_name')
@click.argument('old_set_number', type=int)
@click.argument('new_set_name')
@click.argument('new_set_number', type=int)
def copy_set(old_set_name, old_set_number, new_set_name, new_set_number):
    """Copy a configuration set to a new set.

    Useful for creating variations of existing configurations.
    """
    session = connect()

    try:
        from .msnoise_table_def import Config

        # Get source config set
        source_configs = session.query(Config).filter(
            Config.category == old_set_name,
            Config.set_number == old_set_number
        ).all()

        if not source_configs:
            click.echo(f"Source config set '{old_set_name}' set_number {old_set_number} not found.")
            return

        # Check if destination already exists
        existing = session.query(Config).filter(
            Config.category == new_set_name,
            Config.set_number == new_set_number
        ).first()

        if existing:
            if not click.confirm(f"Config set '{new_set_name}' set_number {new_set_number} already exists. Overwrite?"):
                click.echo("Copy cancelled.")
                return

            # Delete existing entries
            delete_config_set(session, new_set_name, new_set_number)

        # Copy configs
        for config in source_configs:
            new_config = Config(
                name=config.name,
                category=new_set_name,
                set_number=new_set_number,
                value=config.value,
                param_type=config.param_type,
                default_value=config.default_value,
                description=config.description,
                units=config.units,
                possible_values=config.possible_values,
                used_in=config.used_in,
                used=config.used
            )
            session.add(new_config)

        session.commit()
        click.echo(f"✓ Successfully copied {len(source_configs)} config entries from "
                   f"'{old_set_name}' #{old_set_number} to '{new_set_name}' #{new_set_number}")

    except Exception as e:
        session.rollback()
        click.echo(f"Error: {str(e)}")
    finally:
        session.close()

@config.command()
@click.option('--force', is_flag=True,
              help='Force creation even if config sets already exist')
@click.option('--dry-run', is_flag=True,
              help='Show what would be created without actually creating it')
def create_all_sets(force, dry_run):
    """Create one configuration set for each workflow category.

    This command creates a complete set of workflow configurations for each category.
    """

    categories = [
        # 'global',
        'psd',
        'psd_rms',
        'preprocess',
        'cc',
        'filter',
        'stack',
        'mwcs',
        'mwcs_dtt',
        'stretching',
        'wavelet',
        'wavelet_dtt',
    ]

    from ..config import create_config_set
    from ..db import connect

    db = connect()
    for set_name in categories:
        set_number = create_config_set(db, set_name)
        if set_number is not None:
            click.echo(f"Configuration set '{set_name}' created successfully with set_number: {set_number}")
        else:
            click.echo(f"Error: Configuration file for '{set_name}' not found.")

    db.close()

@config.command()
@click.pass_context
def create_workflow_step(ctx):
    """Create a new workflow step interactively"""
    from ..workflow import create_workflow_step

    session = ctx.obj['session']

    step_name = click.prompt('Step name')
    category = click.prompt('Category (preprocess, cc, filter, stack, etc.)')
    set_number = click.prompt('Set number', type=int)
    description = click.prompt('Description (optional)', default='')

    try:
        step = create_workflow_step(session, step_name, category, set_number, description)
        click.echo(f'Created workflow step: {step.step_name} ({step.category}:{step.set_number})')
    except Exception as e:
        click.echo(f'Error: {str(e)}', err=True)

@config.command()
@click.option('--verbose', '-v', is_flag=True, help='Show detailed output')
def create_workflow_steps_from_configs(verbose):
    """Create workflow steps automatically from all existing config sets.

    This command scans all configuration sets in the database and creates
    corresponding workflow steps, sorted by natural workflow order.
    """
    from ..db import connect
    from ..workflow import create_workflow_steps_from_config_sets

    if verbose:
        click.echo("Creating workflow steps")

    try:
        db = connect()
        created_count, existing_count, error_message = create_workflow_steps_from_config_sets(
            db
        )

        if error_message:
            click.echo(f"Error creating workflow steps: {error_message}", err=True)
            return 1

        if created_count > 0:
            click.echo(f"✓ Created {created_count} workflow steps from config sets", color='green')

        if existing_count > 0:
            click.echo(f"ℹ {existing_count} workflow steps already existed", color='yellow')

        if created_count == 0 and existing_count == 0:
            click.echo("⚠ No config sets found to create workflow steps from", color='yellow')

        if verbose and (created_count > 0 or existing_count > 0):
            click.echo(f"Total workflow steps processed: {created_count + existing_count}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        return 1

    return 0

@config.command()
@click.pass_context
def list_workflow_steps(ctx):
    """List all workflow steps"""
    from ..db import connect
    from ..workflow import get_workflow_steps

    session = connect()
    steps = get_workflow_steps(session)
    session.close()

    if not steps:
        click.echo('No workflow steps found.')
        return

    click.echo('Workflow steps:')
    for step in steps:
        click.echo(f'  {step.step_id}: {step.step_name} ({step.category}:{step.set_number})')

@config.command()
@click.pass_context
def show_workflow_graph(ctx):
    """Show workflow graph"""
    from ..db import connect
    from ..workflow import get_workflow_graph

    session = connect()
    graph = get_workflow_graph(session)

    click.echo('Workflow graph:')
    click.echo('\nNodes:')
    for node in graph['nodes']:
        click.echo(f'  {node["id"]}: {node["name"]} ({node["category"]}:{node["set_number"]})')

    click.echo('\nEdges:')
    for edge in graph['edges']:
        click.echo(f'  {edge["from"]} -> {edge["to"]} ({edge["type"]})')

@config.command()
@click.option('--verbose', '-v', is_flag=True, help='Show detailed output')
def create_workflow_links(verbose):
    """Create workflow links automatically between existing workflow steps.

    This command creates links between workflow steps following the natural
    workflow progression: preprocess -> cc -> filter -> stack -> mwcs/stretching/wavelet
    -> mwcs_dtt/wavelet_dtt.

    Links are created based on matching set numbers and workflow logic.
    """
    from ..db import connect
    from ..workflow import create_workflow_links_from_steps

    if verbose:
        click.echo("Creating workflow links")

    try:
        db = connect()
        created_count, existing_count, error_message = create_workflow_links_from_steps(
            db
        )

        if error_message:
            click.echo(f"Error creating workflow links: {error_message}", err=True)
            return 1

        if created_count > 0:
            click.echo(f"✓ Created {created_count} workflow links", color='green')

        if existing_count > 0:
            click.echo(f"ℹ {existing_count} workflow links already existed", color='yellow')

        if created_count == 0 and existing_count == 0:
            click.echo("⚠ No workflow steps found to create links between", color='yellow')

        if verbose and (created_count > 0 or existing_count > 0):
            click.echo(f"Total workflow links processed: {created_count + existing_count}")

    except Exception as e:
        click.echo(f"Error: {str(e)}", err=True)
        return 1

    return 0

@cli.command()
@click.argument('jobtype')
@click.option('-a', '--all', is_flag=True, help='Reset all jobs')
@click.option('-r', '--rule', help='Reset job that match this SQL rule')
def reset(jobtype, all, rule):
    """Resets the jobs to "T"odo. JOBTYPE is the acronym of the job type.
    By default only resets jobs "I"n progress. --all resets all jobs, whatever
    the flag value. Standard Job Types are CC, STACK, MWCS and DTT, but
    plugins can define their own."""
    from ..db import connect, read_db_inifile
    from ..workflow import reset_jobs
    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    session = connect()
    if jobtype == "DA":
        session.execute(text("UPDATE {0}data_availability SET flag='M'"
                        .format(prefix)))
    elif jobtype != jobtype.upper():
        logging.info("The jobtype %s is not uppercase (usually jobtypes"
                     " are uppercase...)"%jobtype)
    reset_jobs(session, jobtype, all, rule)
    session.close()




@cli.command()
@click.option('--fromDA',  help='Populates the station table '
                                'using network and station codes'
                                ' found in the data_availability'
                                ' table, overrides the default'
                                ' workflow step.',
              is_flag=True)
@click.pass_context
def populate(ctx, fromda):
    """Rapidly scan the archive filenames and find Network/Stations, only works
    with known archive structures, or with a custom code provided by the user.
    """
    loglevel = ctx.obj['MSNOISE_verbosity']
    if fromda:
        logger.info("Populating the Station table")
        logger.info("Overriding workflow...")
        db = connect()
        stations = db.query(DataAvailability.net, DataAvailability.sta). \
            group_by(DataAvailability.net, DataAvailability.sta)

        for net, sta in stations:
            logger.info('Adding: %s.%s' % (net, sta))
            X = 0.0
            Y = 0.0
            altitude = 0.0
            coordinates = 'UTM'
            instrument = 'N/A'
            update_station(db, net=net, sta=sta, X=X, Y=Y,
                           altitude=altitude, coordinates=coordinates,
                           instrument=instrument)
        logger.info("Checking the available loc ids and chans...")
        ctx.invoke(db_da_stations_update_loc_chan)
    else:
        from ..s00_populate_station_table import main
        main(loglevel=loglevel)


@cli.command(name='scan_archive')
@click.option('-i', '--init', is_flag=True, help='First run ?')
@click.option('--path',  help='Scan all files in specific folder, overrides the'
                              ' default workflow step.')
@click.option('-r', '--recursively', is_flag=True,
              help='When scanning a path, walk subfolders automatically ?')
@click.option('--crondays',
              help='Number of past days to monitor, typically used in cron jobs '
                   "(overrides the 'crondays' configuration value). "
                   'Must be a float representing a number of days, or designate '
                   " weeks, days, and/or hours using the format 'Xw Xd Xh'.")
@click.pass_context
def scan_archive(ctx, init, crondays, path, recursively):
    """Scan the archive and insert into the Data Availability table."""
    from .. import s01_scan_archive as s01scan_archive
    nthreads = ctx.obj['MSNOISE_threads']
    if path:
        if not os.path.isdir(path):
            logging.critical('Cannot scan from %s: not such directory.' % path)
            sys.exit(1)
        logging.info('Overriding workflow: only scanning path'
                     ' %s (%srecursively)'
                     % (path, '' if recursively else 'non-'))
        s01scan_archive.main(init, nthreads, crondays, path, recursively)
    else:
        s01scan_archive.main(init, nthreads, crondays)


@cli.group(name="plot")
def plot():
    """Commands to trigger plots (data availability, station map)"""
    pass


# Shared plot options — reused across all plot commands
show_option    = click.option('-s', '--show', default=True, type=bool,
                               help='Show figure interactively?')
outfile_option = click.option('-o', '--outfile', default=None, type=str,
                               help='Output filename (?=auto). Supports any matplotlib format, '
                                    'e.g. ?.pdf for PDF with automatic naming.')
pair_type_option = click.option('-p', '--pair_type', default="CC",
                                 type=click.Choice(["CC", "SC", "AC"]),
                                 help='Pair type to plot (CC/SC/AC). Default: CC')
dvvid_option   = click.option('-D', '--dvvid', default=1,
                               help='DVV aggregate config set number')


@plot.command(name='data_availability')
@click.option('-c', '--chan', default="?HZ", help="Channel, you can use the ? wildcard, e.g. '?HZ' (default) or "
                                                  "'HH?', etc.")
@show_option
@outfile_option
@click.pass_context
def plot_data_availability(ctx, chan, show, outfile):
    """Plots the Data Availability vs time"""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from data_availability import main # NOQA
    else:
        from ..plots.data_availability import main
    main(chan, show, outfile, loglevel=loglevel)



@plot.command(name='station_map')
@show_option
@outfile_option
@click.pass_context
def plot_station_map(ctx, show, outfile):
    """Plots the station map (very very basic)"""
    if ctx.obj['MSNOISE_custom']:
        from station_map import main # NOQA
    else:
        from ..plots.station_map import main
    main(show, outfile)


@cli.command(name='new_jobs')
@click.option('-i', '--init', is_flag=True, help='First run ? This disables '
                                                 'the check for existing jobs.')
@click.option('--nocc', is_flag=True, default=False, help='Disable the creation'
                                                          ' of CC jobs.')
@click.option(
    '--after',
    default="",
    help=(
        "Create the next runnable jobs in the workflow based on DONE jobs of the given "
        "config-set type/category (e.g. 'cc', 'stack', 'mwcs'), skipping filter steps in between. "
        "Example: 'msnoise new_jobs --after cc' will create STACK jobs from CC jobs marked D "
        "(via CC -> filter -> stack)."
    ),
)

def new_jobs(init, nocc, after=""):
    """Determines if new jobs are to be defined"""

    from ..s02_new_jobs import main
    main(init, nocc, after)


@cli.group(cls=OrderedGroup)
def cc():
    """Commands for the "Cross-Correlations" Workflow"""
    pass


@cc.command()
@click.option('--threads', default=1, help='Number of threads to use for processing')
@click.option('-v', '--verbose', count=True, help='Increase verbosity')
@click.pass_context
def preprocess(ctx, threads, verbose):
    """Run preprocessing computations on workflow jobs"""
    from ..s02_preprocessing import main
    run_threaded(main, ctx)


@cc.command(name='compute_cc')
@click.pass_context
def cc_compute_cc(ctx):
    """Computes the CC jobs (based on the "New Jobs" identified)"""
    from ..s03_compute_no_rotation import main
    run_threaded(main, ctx)


@cc.command(name='compute_cc_rot')
@click.pass_context
def cc_compute_cc_rot(ctx):
    """Computes the CC jobs too (allows for R or T components)"""
    from msnoise.unused.s03compute_cc import main
    run_threaded(main, ctx)


@cc.command(name="stack")
@click.pass_context
@click.option('-r', '--ref', is_flag=True, help='Compute the REF Stack')
@click.option('-m', '--mov', is_flag=True, help='Compute the MOV Stacks')
@click.option('-s', '--step', is_flag=True, help='Compute the STEP Stacks')
def cc_stack(ctx, ref, mov, step):
    """Stacks the [REF] or [MOV] windows.
    Computes the STACK jobs.
    """
    click.secho('Lets STACK !', fg='green')

    if ref:
        click.secho(
            "ERROR: 'msnoise cc stack -r' is no longer supported.\n"
            "Use 'msnoise cc stack_refstack' instead.",
            fg='red', err=True)
        sys.exit(1)

    if ref and mov:
        click.secho("With MSNoise 1.6, you can't run REF & MOV stacks"
                    "simultaneously, please run them one after the other.")
        sys.exit()

    from ..s04_stack_mov import main
    if mov:
        run_threaded(main, ctx, 'mov')
    if step:
        run_threaded(main, ctx, 'step')


@cc.command(name="stack_refstack")
@click.pass_context
def cc_stack_refstack(ctx):
    """Compute REF stacks for all pending refstack configset jobs.

    Reads ref_begin/ref_end from the refstack configset:

    - Absolute date / 1970-01-01  → Mode A: writes a mean REF file to disk.

    - Negative integer (e.g. -5)  → Mode B: no file written; rolling reference
      computed on-the-fly inside the MWCS/stretching/WCT workers.
    """
    from ..s04_stack_refstack import main
    run_threaded(main, ctx)


@cc.group(name="plot")
def cc_plot():
    """Commands to trigger different plots"""
    pass

# Define reusable option groups as decorator lists
def common_options(*decorators):
    """Combine multiple decorators into one."""
    def decorator(f):
        for dec in reversed(decorators):
            f = dec(f)
        return f
    return decorator


# Define individual options as reusable decorators
preprocessid_option = click.option('-p', '--preprocessid', default=1, help='Preprocessing step ID')
ccid_option = click.option('-cc', '--ccid', default=1, help='CC step ID')
filterid_option = click.option('-f', '--filterid', default=1, help='Filter ID')
stackid_option = click.option('-m', '--stackid', default=1, help='Stack step ID')
stackid_item_option = click.option('-mi', '--stackid_item', default=None, type=int, help='Mov Stack item within that Stack step ID')
refstackid_option = click.option('-rs', '--refstackid', default=1, help='REF Stack step ID')
comp_option = click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ')

mwcsid_option = click.option('-w', '--mwcsid', default=1, help='MWCS config set number')
mwcsdttid_option   = click.option('-wi', '--mwcsdttid', default=1, help='MWCS-DTT config set number')
stretchingid_option = click.option('-S', '--stretchingid', default=1, help='Stretching config set number')
wctid_option        = click.option('-w', '--wctid', default=1, help='WCT config set number')
wctdttid_option     = click.option('-d', '--wctdttid', default=1, help='WCT-DTT config set number')

# Build incremental option bundles (CC side)
base_options = common_options(preprocessid_option, ccid_option, filterid_option)
stack_options = common_options(base_options, stackid_option, stackid_item_option, refstackid_option)
full_options  = common_options(stack_options, comp_option)

mwcs_options = common_options(stack_options, mwcsid_option, comp_option)
mwcsdtt_options = common_options(mwcs_options, mwcsdttid_option)



# DTT-side lineage options (preprocess → cc → filter → stack → mwcs → mwcs_dtt)
mwcsid_option  = click.option('-w', '--mwcsid',  default=1, help='MWCS step set number')
dttid_option   = click.option('-d', '--dttid',   default=1, help='MWCS-DTT step set number')
# Full DTT lineage bundle: preprocess / cc / filter / stack / stack_item / mwcs / dtt / comp
dtt_lineage_options = common_options(
    preprocessid_option, ccid_option, filterid_option,
    stackid_option, stackid_item_option, refstackid_option,
    mwcsid_option, dttid_option, comp_option,
)

@cc_plot.command(name="distance",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZE, NZ, 1E,...). Defaults to ZZ')
@click.option('-a', '--ampli', default=1.0, help='Amplification of the individual lines on the vertical axis ('
                                                 'default=1)')
@show_option
@outfile_option
@click.option('-r', '--refilter', default=None,
              help='Refilter CCFs before plotting (e.g. 4:8 for filtering CCFs '
                   'between 4.0 and 8.0 Hz. This will update the plot title.')
@click.option('--virtual-source', default=None,
              help='Use only pairs including this station. Format must be '
                   'NET.STA')
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED,
                callback=parse_extra_args)
@click.pass_context
def cc_plot_distance(ctx, filterid, comp, ampli, show, outfile, refilter,
                     virtual_source, extra_args):
    """Plots the REFs of all pairs vs distance"""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from distance import main # NOQA
    else:
        from ..plots.distance import main
    main(filterid, comp, ampli, show, outfile, refilter, virtual_source,
         loglevel=loglevel, **extra_args)


@cc_plot.command(name="interferogram",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.argument('sta1')
@click.argument('sta2')
@full_options
@show_option
@outfile_option
@click.option('-r', '--refilter', default=None,
              help='Refilter CCFs before plotting (e.g. 4:8 for filtering CCFs '
                   'between 4.0 and 8.0 Hz. This will update the plot title.')
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED,
                callback=parse_extra_args)
@click.pass_context
def cc_plot_interferogram(ctx, sta1, sta2,  preprocessid, ccid, filterid, stackid, stackid_item,
                    comp, show,
                          outfile,
                          refilter, extra_args):
    """Plots the interferogram between sta1 and sta2 (parses the CCFs)
    STA1 and STA2 must be provided with this format: NET.STA.LOC !"""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from interferogram import main # NOQA
    else:
        from ..plots.interferogram import main
    main(sta1, sta2, preprocessid, ccid, filterid, stackid, stackid_item,
         comp, show, outfile, refilter,
         loglevel=loglevel, **extra_args)


@cc_plot.command(name="ccftime",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.argument('sta1')
@click.argument('sta2')
@full_options
@click.option('-a', '--ampli', default=5.0, help='Amplification of the individual lines on the vertical axis ('
                                                 'default=1)')
@click.option('-S', '--seismic', is_flag=True, help='Seismic style: fill the space between the zero and the positive '
                                                    'wiggles')
@show_option
@outfile_option
@click.option('-e', '--envelope', is_flag=True, help='Plot envelope instead of '
                                                     'time series')
@click.option('-r', '--refilter', default=None,
              help='Refilter CCFs before plotting (e.g. 4:8 for filtering CCFs '
                   'between 4.0 and 8.0 Hz. This will update the plot title.')
@click.option("--normalize", default="individual")
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED,
                callback=parse_extra_args)
@click.pass_context
def cc_plot_ccftime(ctx, sta1, sta2, preprocessid, ccid, filterid, stackid, stackid_item,
                    comp,
                    ampli, seismic, show, outfile, envelope, refilter,
                    normalize, extra_args):
    """Plots the ccf vs time between sta1 and sta2
    STA1 and STA2 must be provided with this format: NET.STA.LOC !"""
    loglevel = ctx.obj['MSNOISE_verbosity']
    # if sta1 > sta2:
    #     click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
    #     return
    if ctx.obj['MSNOISE_custom']:
        from ccftime import main # NOQA
    else:
        from ..plots.ccftime import main
    main(sta1, sta2, preprocessid, ccid, filterid, stackid, stackid_item, comp,
         ampli, seismic, show, outfile,
         envelope, refilter, normalize, loglevel=loglevel, **extra_args)


@cc_plot.command(name="spectime",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.argument('sta1')
@click.argument('sta2')
@full_options
@click.option('-a', '--ampli', default=5.0, help='Amplification of the individual lines on the vertical axis ('
                                                 'default=1)')
@show_option
@outfile_option
@click.option('-r', '--refilter', default=None,
              help='Refilter CCFs before plotting (e.g. 4:8 for filtering CCFs '
                   'between 4.0 and 8.0 Hz. This will update the plot title.')
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED,
                callback=parse_extra_args)
@click.pass_context
def cc_plot_spectime(ctx, sta1, sta2, preprocessid, ccid, filterid, stackid, stackid_item,
                    comp, ampli, show, outfile, refilter, extra_args):
    """Plots the ccf's spectrum vs time between sta1 and sta2
    STA1 and STA2 must be provided with this format: NET.STA.LOC !"""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from spectime import main # NOQA
    else:
        from ..plots.spectime import main
    main(sta1, sta2, preprocessid, ccid, filterid, stackid, stackid_item, comp,
         ampli,  show=show, outfile=outfile,
         refilter=refilter, loglevel=loglevel, **extra_args)


@cc.group(cls=OrderedGroup)
def dtt():
    """Commands for the "Delay Time / dt/t" Workflow"""
    pass


@dtt.command(name='compute_mwcs')
@click.pass_context
def dtt_compute_mwcs(ctx):
    """Computes the MWCS jobs"""
    from ..s05_compute_mwcs import main
    run_threaded(main, ctx)


@dtt.command(name='compute_mwcs_dtt')
@click.pass_context
def dtt_compute_dtt(ctx):
    """Computes the dt/t jobs based on the new MWCS data"""
    from ..s06_compute_mwcs_dtt import main
    run_threaded(main, ctx)


@dtt.command(name='compute_stretching')
@click.pass_context
def dtt_compute_stretching2(ctx):
    """Computes the stretching based on the new stacked data"""
    from ..s10_stretching import main
    run_threaded(main, ctx)


@dtt.command(name='compute_wct')
@click.pass_context
def dtt_compute_wct(ctx):
    """Computes the wavelet jobs based on the new STACK data"""
    from ..s08_compute_wct import main
    run_threaded(main, ctx)

@dtt.command(name='compute_wct_dtt')
@click.pass_context
def dtt_compute_wct_dtt(ctx):
    """Computes dv/v from WCT results (wavelet_dtt step, lineage-based)"""
    from ..s09_compute_wct_dtt import main
    run_threaded(main, ctx)

@dtt.group(name="plot")
def dtt_plot():
    """Commands to trigger different plots"""
    pass


@dtt.group(name="dvv")
def dtt_dvv():
    """Commands for dv/v aggregation across station pairs."""
    pass


@dtt_dvv.command(name='compute_mwcs_dtt_dvv')
@click.pass_context
def dtt_dvv_compute_mwcs(ctx):
    """Aggregate MWCS dv/v across station pairs (mwcs_dtt_dvv step)."""
    from ..s07_compute_dvv import main
    run_threaded(main, ctx, "mwcs_dtt_dvv")


@dtt_dvv.command(name='compute_stretching_dvv')
@click.pass_context
def dtt_dvv_compute_stretching(ctx):
    """Aggregate Stretching dv/v across station pairs (stretching_dvv step)."""
    from ..s07_compute_dvv import main
    run_threaded(main, ctx, "stretching_dvv")


@dtt_dvv.command(name='compute_wavelet_dtt_dvv')
@click.pass_context
def dtt_dvv_compute_wct(ctx):
    """Aggregate WCT dv/v across station pairs (wavelet_dtt_dvv step). Supports multi-band extraction."""
    from ..s07_compute_dvv import main
    run_threaded(main, ctx, "wavelet_dtt_dvv")


@dtt_dvv.group(name="plot")
def dtt_dvv_plot():
    """Commands to plot aggregated dv/v results."""
    pass


@dtt_dvv_plot.command(name="mwcs_dvv")
@mwcsdtt_options
@dvvid_option
@click.option('-M', '--dttname', default="m", help='DTT column: m (slope) or m0 (zero-intercept slope)')
@pair_type_option
@show_option
@outfile_option
@click.pass_context
def dtt_dvv_plot_mwcs(ctx, preprocessid, ccid, filterid, stackid, stackid_item,
                      refstackid, mwcsid, mwcsdttid, dvvid, comp, dttname, pair_type, show, outfile):
    """Plot dv/v from MWCS-DTT aggregate. Requires mwcs_dtt_dvv step."""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from dvv import main  # NOQA
    else:
        from ..plots.mwcs_dtt_dvv import main
    main(preprocessid, ccid, filterid, stackid, stackid_item, refstackid,
         mwcsid, mwcsdttid, dvvid=dvvid, components=comp,
         dttname=dttname, pair_type=pair_type, show=show, outfile=outfile, loglevel=loglevel)


@dtt_dvv_plot.command(name="stretching_dvv")
@filterid_option
@stretchingid_option
@dvvid_option
@comp_option
@click.option('-m', '--mov_stack', default=0, help='Plot specific mov stack (1-based index, 0=all)')
@pair_type_option
@show_option
@outfile_option
@click.pass_context
def dtt_dvv_plot_stretching(ctx, mov_stack, comp, filterid, stretchingid,
                             dvvid, pair_type, show, outfile):
    """Plot dv/v from Stretching aggregate. Requires stretching_dvv step."""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from dvvs import main  # NOQA
    else:
        from ..plots.stretching_dvv import main
    main(mov_stackid=mov_stack, components=comp, filterid=filterid,
         stretchingid=stretchingid, dvvid=dvvid, pair_type=pair_type,
         show=show, outfile=outfile, loglevel=loglevel)


@dtt_dvv_plot.command(name="wavelet_dvv")
@filterid_option
@comp_option
@click.option('-m', '--mov_stack', default=0, help='Plot specific mov stack (1-based index, 0=all)')
@wctid_option
@wctdttid_option
@dvvid_option
@pair_type_option
@click.option('-v', '--visualize', default="timeseries",
              type=click.Choice(["timeseries", "heatmap"]),
              help='Plot style: timeseries (uses dvv aggregate) or heatmap (uses per-pair data)')
@click.option('-r', '--ranges', default="[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]",
              help='Frequency ranges for band averaging (first range used for timeseries)')
@show_option
@outfile_option
@click.pass_context
def dtt_dvv_plot_wct(ctx, filterid, comp, mov_stack, wctid, wctdttid, dvvid,
                     pair_type, visualize, ranges, show, outfile):
    """Plot dv/v from WCT-DTT aggregate. Use -v heatmap for per-pair frequency view."""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from wavelet_dtt_dvv import main  # NOQA
    else:
        from ..plots.wavelet_dtt_dvv import main
    main(mov_stackid=mov_stack, components=comp, filterid=filterid,
         wctid=wctid, dttid=wctdttid, dvvid=dvvid, pair_type=pair_type,
         visualize=visualize, ranges=ranges, show=show, outfile=outfile, loglevel=loglevel)




@dtt_plot.command(name='mwcs')
@click.argument('sta1')
@click.argument('sta2')
@dtt_lineage_options
@show_option
@outfile_option
@click.pass_context
def dtt_plot_mwcs(ctx, sta1, sta2, preprocessid, ccid, filterid, stackid,
                  stackid_item, refstackid, mwcsid, dttid, comp, show, outfile):
    """Plots the MWCS dt and coherence images for a station pair.
    STA1 and STA2 must be provided as NET.STA.LOC.
    Lineage: -p preprocess, -cc cc, -f filter, -m stack, -mi stack_item, -w mwcs."""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from mwcs import main  # NOQA
    else:
        from ..plots.mwcs import main
    main(sta1, sta2,
         preprocessid=preprocessid, ccid=ccid, filterid=filterid,
         stackid=stackid, stackid_item=stackid_item or 1, refstackid=refstackid,
         mwcsid=mwcsid, mwcsdttid=dttid,
         components=comp, show=show, outfile=outfile, loglevel=loglevel)


@dtt_plot.command(name="mwcs_dtt_day")
@click.argument('sta1')
@click.argument('sta2')
@click.argument('day')
@dtt_lineage_options
@show_option
@outfile_option
@click.pass_context
def dtt_plot_dtt(ctx, sta1, sta2, day, preprocessid, ccid, filterid,
                 stackid, stackid_item, mwcsid, dttid, comp, show, outfile):
    """Plots dt against t (scatter + regression) for a single day.
    STA1, STA2: NET.STA.LOC. DAY: YYYY-MM-DD.
    Lineage: -p preprocess, -cc cc, -f filter, -m stack, -mi stack_item, -w mwcs, -d dtt."""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from dtt import main  # NOQA
    else:
        from ..plots.dtt import main
    main(sta1, sta2, filterid=filterid, components=comp, day=day,
         preprocessid=preprocessid, ccid=ccid, stackid=stackid,
         stackid_item=stackid_item, mwcsid=mwcsid, mwcsdttid=dttid,
         show=show, outfile=outfile, loglevel=loglevel)


@dtt_plot.command(name='mwcs_dtt_timing')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-w', '--mwcsid', default=1, help='MWCS step set number')
@click.option('-d', '--dttid', default=1, help='MWCS-DTT step set number')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...). Defaults to ZZ')
@click.option('-m', '--mov_stack', default=0, help='Plot specific mov stack (1-based index, 0=all)')
@click.option('-p', '--pair', default=None,
              help='Highlight a specific pair (NET.STA.LOC:NET.STA.LOC)', multiple=True)
@click.option('-M', '--dttname', default="m",
              help='DTT column: m (slope=dt/t) or m0 (zero-intercept)')
@show_option
@outfile_option
@click.pass_context
def dtt_plot_mwcs_timing(ctx, mov_stack, comp, dttname, filterid, mwcsid, dttid,
                         pair, show, outfile):
    """Plots network-mean dt/t timeseries from MWCS-DTT results.
    Optionally highlight specific pairs with -p NET.STA.LOC:NET.STA.LOC."""
    loglevel = ctx.obj['MSNOISE_verbosity']
    if ctx.obj['MSNOISE_custom']:
        from timing import main  # NOQA
    else:
        from ..plots.timing import main
    main(mov_stackid=mov_stack, dttname=dttname, components=comp,
         filterid=filterid, mwcsid=mwcsid, dttid=dttid,
         pairs=list(pair), show=show, outfile=outfile, loglevel=loglevel)


# PLOT GROUP
#


@cli.group(cls=OrderedGroup)
def qc():
    """Commands for computing PSD, RMS, etc..."""
    pass


@qc.command(name='compute_psd')
@click.option('-n', '--njobs_per_worker', default=9999,
              help='Reduce this number when processing a small number of days '
                   'but a large number of stations')
@click.pass_context
def qc_compute_psd(ctx, njobs_per_worker):
    """Computes the PSD jobs, saves results as NetCDF files.
       Based on New or Modified files identified by the new_jobs step."""
    from ..psd_compute import main
    run_threaded(main, ctx, njobs_per_worker=njobs_per_worker)




@qc.command(name='plot_psd')
@click.argument('seed_id')
@click.pass_context
def qc_plot_psd(ctx, seed_id):
    """Plots the PSD and spectrogram based on NPZ files"""
    from ..plots.ppsd import main
    net,sta,loc,chan = seed_id.split(".")
    main(net, sta, loc, chan, time_of_weekday=None,
         period_lim=(0.02, 50.0), cmap="viridis",
         color_lim=None, show=True)


@qc.command(name='compute_psd_rms')
@click.pass_context
def qc_compute_rms(ctx):
    """Computes the RMS from PSD NetCDF files."""
    from ..psd_compute_rms import main
    run_threaded(main, ctx)






# @with_plugins(iter_entry_points('msnoise.plugins'))
@cli.group()
def plugin():
    """Runs a command in a named plugin"""
    pass

try:
    db = connect()
    plugins = get_config(db, "plugins")
    db.close()
except DBConfigNotFoundError:
    plugins = None
except sqlalchemy.exc.OperationalError as e:
    logging.critical('Unable to read project configuration: error connecting to the database:{}'.format(str(e)))
    sys.exit(1)

if plugins:
    from importlib.metadata import entry_points
    plugins = plugins.split(",")
    for ep in list(entry_points(group='msnoise.plugins.commands')):
        module_name = ep.value.split(".")[0]
        if module_name in plugins:
            plugin.add_command(ep.load())


@cli.group(cls=OrderedGroup)
def utils():
    """Command group for smaller tools"""
    pass

@utils.command(name="bugreport")
@click.option('-s', '--sys', is_flag=True, help='System Info')
@click.option('-m', '--modules', is_flag=True, help='Modules Info')
@click.option('-e', '--env', is_flag=True, help='Environment Info')
@click.option('-a', '--all', is_flag=True, help='All Info')
@click.pass_context
def utils_bugreport(ctx, sys, modules, env, all):
    """This command launches the Bug Report script."""
    click.echo('Let\'s Bug Report MSNoise !')
    # click.echo('Working on %i threads' % ctx.obj['MSNOISE_threads'])
    from ..bugreport import main
    main(sys, modules, env, all)


@utils.command(name="test")
@click.option('-p', '--prefix', default="", help='Prefix for tables')
@click.option('--tech', default=1, help='Test using (1) SQLite or (2) MariaDB (you need to start that server before!)')
@click.option('-c', '--content', default=False, is_flag=True)
def utils_test(prefix, tech, content):
    """Runs the test suite in a temporary folder"""
    import matplotlib.pyplot as plt
    import pytest
    plt.switch_backend("agg")

    # Prepare environment variables for the test session
    os.environ["PREFIX"] = prefix
    os.environ["TECH"] = str(tech)

    # Determine which test suite to run
    test_module = 'content_tests' if content else 'tests'

    # Construct the path to the test module
    test_path = os.path.join(os.path.dirname(__file__), '..', 'test', f'{test_module}.py')

    # Run pytest on the selected test module
    exit_code = pytest.main(['-s', test_path])

    # Handle the exit code as needed
    if exit_code != 0:
        print("Tests failed.")

@utils.command(name="jupyter")
def utils_jupyter():
    """Launches an jupyter notebook in the current folder"""
    os.system("jupyter notebook --ip 0.0.0.0 --no-browser")


@utils.command(name="export-params")
@click.option('-l', '--lineage', default=None, type=str,
              help='Lineage string e.g. "preprocess_1/cc_1/filter_1/stack_1/mwcs_1/mwcs_dtt_1"')
@preprocessid_option
@ccid_option
@filterid_option
@stackid_option
@click.option('-rs', '--refstackid', default=None, type=int,
              help='REF Stack step ID (optional)')
@click.option('-w', '--mwcsid', default=None, type=int,
              help='MWCS step ID (optional)')
@click.option('-wi', '--mwcsdttid', default=None, type=int,
              help='MWCS-DTT step ID (optional)')
@click.option('-S', '--stretchingid', default=None, type=int,
              help='Stretching step ID (optional)')
@click.option('-sd', '--stretchingdvvid', default=None, type=int,
              help='Stretching-DVV step ID (optional)')
@click.option('-W', '--wctid', default=None, type=int,
              help='WCT step ID (optional)')
@click.option('-wd', '--wctdttid', default=None, type=int,
              help='WCT-DTT step ID (optional)')
@click.option('-wdv', '--waveletdvvid', default=None, type=int,
              help='Wavelet-DVV step ID (optional)')
@click.option('-o', '--output', default=None, type=str,
              help='Output YAML path (default: params_<lineage>.yaml in current dir)')
@click.pass_context
def utils_export_params(ctx, lineage, preprocessid, ccid, filterid, stackid,
                        refstackid, mwcsid, mwcsdttid,
                        stretchingid, stretchingdvvid,
                        wctid, wctdttid, waveletdvvid, output):
    """Export the full layered parameter chain for a lineage to YAML.

    Either supply --lineage as a slash-separated string, or build it from
    individual step IDs.  The resulting YAML contains one block per
    config category in lineage order — no key collisions, fully
    self-describing for reproducibility.

    Examples::

        msnoise utils export-params --lineage preprocess_1/cc_1/filter_1/stack_1
        msnoise utils export-params -p 1 -cc 1 -f 1 -m 1 -w 1 -wi 1
        msnoise utils export-params -p 1 -cc 1 -f 1 -m 1 -S 1 -sd 1
        msnoise utils export-params -p 1 -cc 1 -f 1 -m 1 -W 1 -wd 1 -wdv 1
    """
    from ..config import get_params
    from ..db import connect
    from ..workflow import lineage_str_to_steps
    from ..config import get_merged_params_for_lineage

    db = connect()
    orig_params = get_params(db)

    if lineage:
        lin_str = lineage
    else:
        # Build from integer IDs
        parts = [f"preprocess_{preprocessid}", f"cc_{ccid}",
                 f"filter_{filterid}", f"stack_{stackid}"]
        if refstackid:     parts.append(f"refstack_{refstackid}")
        if mwcsid:         parts.append(f"mwcs_{mwcsid}")
        if mwcsdttid:      parts.append(f"mwcs_dtt_{mwcsdttid}")
        if stretchingid:   parts.append(f"stretching_{stretchingid}")
        if stretchingdvvid: parts.append(f"stretching_dvv_{stretchingdvvid}")
        if wctid:          parts.append(f"wavelet_{wctid}")
        if wctdttid:       parts.append(f"wavelet_dtt_{wctdttid}")
        if waveletdvvid:   parts.append(f"wavelet_dtt_dvv_{waveletdvvid}")
        lin_str = "/".join(parts)

    steps = lineage_str_to_steps(db, lin_str, sep="/", strict=False)
    _, lineage_names, params = get_merged_params_for_lineage(db, orig_params, {}, steps)

    out_path = output or f"params_{lin_str.replace('/', '_')}.yaml"
    params.to_yaml(out_path)
    click.echo(f"Exported params for lineage '{lin_str}' → {out_path}")
    db.close()




def run():
    try:
        cli(obj={})
    except MSNoiseError as e:
        logging.critical(str(e))
