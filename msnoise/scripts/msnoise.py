import traceback
import logging
import os
import sys
import sqlalchemy
from sqlalchemy import text
import time

import click


from .._version import get_git_version
__version__ = get_git_version(dirty=True)

from .. import MSNoiseError, DBConfigNotFoundError
from ..api import connect, get_config, update_station, get_logger, get_job_types
from ..msnoise_table_def import DataAvailability

class OrderedGroup(click.Group):
    def list_commands(self, ctx):
        return self.commands.keys()


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
    from ..api import get_config
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
    from ..api import read_db_inifile
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
    from ..api import get_config
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
            for job in get_job_types(db):
                if job[1] == 'D':
                    if job[0] > 0:
                        click.secho(
                            " - %s does not exists and that is not normal"
                            " (%i CC jobs done)" % (output_folder, job[0]),
                            fg='red')
                    else:
                        click.secho(
                            " - %s does not exists and that is normal"
                            " (%i CC jobs done)" % (output_folder, job[0]))
        else:
            click.secho(
                " - %s does not exists (and that is normal because"
                " keep_all=False)" % output_folder)


def info_parameters(db):
    """
    Show values of each configuration parameters.
    """
    from ..api import get_filters
    from ..default import default
    click.echo('')
    click.echo('Configuration values:'
            '   | Normal colour indicates that the default value is used'
            '   | Green indicates "M"odified values')
    # TODO: add plugins params
    show_config_values(db, default.keys())

    click.echo('')
    click.echo('Filters:')
    click.echo(' ID:   [low:high]   [mwcs_low:mwcs_high] mwcs_wlen mwcs_step Used?')

    for f in get_filters(db, all=True):
        click.echo(' {:2d}: {:^15s} {:^20s} {:^9s} {:^9s}  {:1s}'
            .format(f.ref,
                '[{:.3f}:{:.3f}]'.format(f.low, f.high),
                '[{:.3f}:{:.3f}]'.format(f.mwcs_low, f.mwcs_high),
                '{:.0f}'.format(f.mwcs_wlen),
                '{:.0f}'.format(f.mwcs_step),
                'Y' if f.used else 'N'))


def info_stations(db):
    """
    Show information about configured stations.
    """
    from ..api import get_stations
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
    if s is None:
        click.echo(' ')


def info_jobs(db):
    """
    Show information about jobs registered in database.
    """
    from ..api import get_job_types

    jobtypes = {}
    jobtypes["QC"] = ["PSD", "PSD2HDF", "HDF2RMS"]
    jobtypes["CC"] = ["CC", "STACK", "MWCS", "DTT", "DVV"]

    click.echo("Jobs:")
    for category in ["QC", "CC"]:
        click.echo(' %s:' % category)
        for jobtype in jobtypes[category]:
            click.echo('  %s:' % jobtype)
            n = None
            for (n, jobtype) in get_job_types(db, jobtype):
                click.echo("   %s : %i" % (jobtype, n))
            if n is None:
                click.echo('   none defined')


def info_plugins(db):
    """
    Show information about configured plugins.
    """
    from ..api import get_config, get_job_types
    plugins = get_config(db, "plugins")
    if not plugins:
        return
    plugins = plugins.split(",")
    import pkg_resources
    for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.jobtypes'):
        module_name = ep.module_name.split(".")[0]
        if module_name in plugins:
            click.echo('')
            click.echo('Plugin: %s' % module_name)
            for row in ep.load()():
                click.echo(' %s:' % row["name"])
                for (n, jobtype) in get_job_types(db, row["name"]):
                    click.echo("  %s : %i" % (jobtype, n))

from pkg_resources import iter_entry_points

# import click
# from click_plugins import with_plugins
# @with_plugins(iter_entry_points('click_command_tree'))
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
    # sys.path.append(os.getcwd())



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
@click.option('--tech', help='Database technology: 1=SQLite 2=MySQL 3=PostgreSQL',
              default=None)
def db_init(tech):
    """This command initializes the current folder to be a MSNoise Project
    by creating a database and a db.ini file."""
    click.echo('Launching the init')
    from ..s000installer import main
    main(tech)


@db.command(name="update_loc_chan")
def db_da_stations_update_loc_chan():
    """Populates the Location & Channel from the Data Availability
    table. Warning: rewrites automatically, no confirmation."""
    from msnoise.api import connect, get_stations

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
        print("%s.%s has locids:%s and chans:%s" % (sta.net, sta.sta,
                                                    locids, chans))
        sta.used_location_codes = ",".join(locids)
        sta.used_channel_names = ",".join(chans)
        session.commit()

@db.command(name="execute")
@click.argument('sql_command')
@click.option('-o', '--outfile', help='Output filename (?="request.csv")',
              default=None, type=str)
@click.option('-s', '--show', help='Show output (in case of SELECT statement)?',
              default=True, type=bool)
def db_execute(sql_command, outfile=None, show=True):
    """EXPERT MODE: Executes 'sql_command' on the database. Use this command
    at your own risk!!"""
    from msnoise.api import connect

    db = connect()
    for cmd in sql_command.split(";"):
        if not len(cmd):
            continue
        print("Executing '%s'" % cmd)
        r = db.execute(cmd)
        if cmd.count("select") or cmd.count("SELECT"):
            result = r.fetchall()
            if not len(result):
                print("The query returned no results, sorry.")
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
    from ..api import connect, Config, read_db_inifile
    from ..default import default
    db = connect()
    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    for name in default.keys():
        try:
            db.add(Config(name=name, value=default[name].default))
            db.commit()
        except:
            db.rollback()
            # print("Passing %s: already in DB" % name)
            continue
    try:
        db.execute("CREATE UNIQUE INDEX job_index ON %sjobs (day, pair, "
                   "jobtype)" %
                   prefix)
        db.commit()
    except:
        logging.info("It looks like the v1.5 'job_index' is already in the DB")
        db.rollback()

    try:
        db.execute("CREATE INDEX job_index2 ON %sjobs (jobtype, flag)" %
                   prefix)
        db.commit()
    except:
        logging.info("It looks like the v1.6 'job_index2' is already in the DB")
        db.rollback()

    try:
        db.execute("CREATE UNIQUE INDEX da_index ON %sdata_availability (path, "
                   "file, net, sta, loc, chan)" %
                   prefix)
        db.commit()
    except:
        logging.info("It looks like the v1.5 'da_index' is already in the DB")
        db.rollback()

    db.close()


@db.command(name='clean_duplicates')
def db_clean_duplicates():
    """Checks the Jobs table and deletes duplicate entries"""
    from msnoise.api import connect, read_db_inifile

    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    db = connect()
    if dbini.tech == 1:
        query = 'DELETE FROM {0}jobs WHERE rowid NOT IN '\
                '(SELECT MIN(rowid) FROM {0}jobs GROUP BY day,pair,jobtype)'\
                .format(prefix)
    else:
        query = 'DELETE from {0}jobs USING {0}jobs as j1, {0}jobs as j2 '\
                'WHERE (j1.ref > j2.ref) AND (j1.day=j2.day) '\
                'AND (j1.pair=j2.pair) AND (j1.jobtype=j2.jobtype)'\
                .format(prefix)
    db.execute(query)
    db.commit()
    db.close()


@db.command(name="dump")
@click.option("--format", default="csv")
def db_dump(format):
    """Dumps the complete database in formatted files, defaults to CSV.
    """
    from ..api import connect, get_engine
    from sqlalchemy import MetaData
    import pandas as pd

    if format == "csv":
        engine = get_engine(inifile=os.path.join(os.getcwd(), 'db.ini'))

        meta = MetaData()
        meta.reflect(bind=engine)

        for table in meta.sorted_tables:
            r = [dict(row) for row in engine.execute(table.select())]
            df = pd.DataFrame(r)
            print("Dumping table %s to %s.csv" % (table.name, table.name))
            df.to_csv("%s.csv" % table.name, index=False)
    else:
        print("Currently only the csv format is supported, sorry.")


@db.command(name="import")
@click.argument("table")
@click.option("--format", default="csv")
@click.option("--force", is_flag=True, default=False)
def db_import(table, format, force):
    """
    Imports msnoise tables from formatted files (CSV).
    """
    from ..api import connect, get_engine
    from sqlalchemy import MetaData
    import pandas as pd

    if format == "csv":
        engine = get_engine(inifile=os.path.join(os.getcwd(), 'db.ini'))
        print("Loading table %s from %s.csv" % (table, table))
        df = pd.read_csv("%s.csv" % table)
        if force:
            df.to_sql(table, engine, if_exists="replace")
        else:
            try:
                df.to_sql(table, engine)
            except ValueError:
                traceback.print_exc()
                print("!"*80)
                print("You're probably getting the error above because the "
                      "table already exists, if you want to replace the table "
                      "with the imported data, then pass the --force option")
    else:
        print("Currently only the csv format is supported, sorry.")



@cli.command()
@click.option('-j', '--jobs', is_flag=True, help='Jobs Info only')
def info(jobs):
    """Outputs general information about the current install and config, plus
    information about jobs and their status."""
    from ..api import connect

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
    import glob
    from ..api import connect, get_config, get_stations, update_station,\
        preload_instrument_responses

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
            print("No coords for %s, skipping,..." % id)
            continue
        lon = float(coords["longitude"].values[0])
        lat = float(coords["latitude"].values[0])
        update_station(db, station.net, station.sta, lon, lat, 0, "DEG", )
        logging.info("Added coordinates (%.5f %.5f) for station %s.%s" %
                    (lon, lat, station.net, station.sta))
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
    from ..api import connect, update_config
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
    from ..api import connect, get_config
    db = connect()
    show_config_values(db, names)
    db.close()



@cli.command()
@click.argument('jobtype')
@click.option('-a', '--all', is_flag=True, help='Reset all jobs')
@click.option('-r', '--rule', help='Reset job that match this SQL rule')
def reset(jobtype, all, rule):
    """Resets the jobs to "T"odo. JOBTYPE is the acronym of the job type.
    By default only resets jobs "I"n progress. --all resets all jobs, whatever
    the flag value. Standard Job Types are CC, STACK, MWCS and DTT, but
    plugins can define their own."""
    from ..api import connect, reset_jobs, read_db_inifile
    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    session = connect()
    if jobtype == "DA":
        session.execute("UPDATE {0}data_availability SET flag='M'"
                        .format(prefix))
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
def populate(fromda):
    """Rapidly scan the archive filenames and find Network/Stations, only works
    with known archive structures, or with a custom code provided by the user.
    """
    if fromda:
        logging.info("Overriding workflow...")
        db = connect()
        stations = db.query(DataAvailability.net, DataAvailability.sta). \
            group_by(DataAvailability.net, DataAvailability.sta)

        for net, sta in stations:
            print('Adding:', net, sta)
            X = 0.0
            Y = 0.0
            altitude = 0.0
            coordinates = 'UTM'
            instrument = 'N/A'
            update_station(db, net, sta, X, Y, altitude,
                           coordinates=coordinates, instrument=instrument)
    else:
        from ..s002populate_station_table import main
        main()


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
    from .. import s01scan_archive
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



@plot.command(name='data_availability')
@click.option('-c', '--chan', default="?HZ", help="Channel, you can use the ? wildcard, e.g. '?HZ' (default) or 'HH?', etc.")
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def plot_data_availability(ctx, chan, show, outfile):
    """Plots the Data Availability vs time"""
    if ctx.obj['MSNOISE_custom']:
        from data_availability import main
    else:
        from ..plots.data_availability import main
    main(chan, show, outfile)



@plot.command(name='station_map')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def plot_station_map(ctx, show, outfile):
    """Plots the station map (very very basic)"""
    if ctx.obj['MSNOISE_custom']:
        from station_map import main
    else:
        from ..plots.station_map import main
    main(show, outfile)


@cli.command(name='new_jobs')
@click.option('-i', '--init', is_flag=True, help='First run ? This disables '
                                                 'the check for existing jobs.')
@click.option('--nocc', is_flag=True, default=False, help='Disable the creation'
                                                          ' of CC jobs.')
@click.option('--hpc', help='Format PREVIOUS:NEXT. When running on HPC, '
                            'create the next jobs in the workflow based on the'
                            'previous step mentioned here. Example:'
                            '"msnoise new_jobs --hpc CC:STACK" will create '
                            'STACK jobs based on CC jobs marked "D"one.')
def new_jobs(init, nocc, hpc=""):
    """Determines if new CC/QC jobs are to be defined"""
    if not hpc:
        from ..s02new_jobs import main
        main(init, nocc)
    if hpc:
        from ..api import connect, read_db_inifile
        dbini = read_db_inifile()
        prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
        left, right = hpc.split(':')
        db = connect()
        db.execute("INSERT INTO {prefix}jobs (pair, day, jobtype, flag) "
                   "SELECT pair, day, '{right_type}', 'T' FROM {prefix}jobs "
                   "WHERE jobtype='{left_type}' AND flag='D';"
                   .format(prefix=prefix, right_type=right, left_type=left))
        db.commit()
        db.close()



@cli.group(cls=OrderedGroup)
def cc():
    """Commands for the "Cross-Correlations" Workflow"""
    pass


@cc.command(name='compute_cc')
@click.pass_context
def cc_compute_cc(ctx):
    """Computes the CC jobs (based on the "New Jobs" identified)"""
    from ..s03compute_no_rotation import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main, kwargs={"loglevel": loglevel})
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


@cc.command(name='compute_cc_rot')
@click.pass_context
def cc_compute_cc_rot(ctx):
    """Computes the CC jobs (based on the "New Jobs" identified)"""
    from ..s03compute_cc import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main, kwargs={"loglevel": loglevel})
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


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
    from ..s04stack import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']

    if ref and mov:
        click.secho("With MSNoise 1.6, you can't run REF & MOV stacks"
                    "simultaneously, please run them one after the other.")
        sys.exit()

    if threads == 1:
        if ref:
            main('ref', loglevel=loglevel)
        if mov:
            main('mov', loglevel=loglevel)
        if step:
            main('step', loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        if ref:
            for i in range(threads):
                p = Process(target=main, args=["ref",],
                            kwargs={"loglevel": loglevel})
                p.start()
                processes.append(p)
                time.sleep(delay)
        for p in processes:
            p.join()
        if mov:
            for i in range(threads):
                p = Process(target=main, args=["mov",],
                            kwargs={"loglevel": loglevel})
                p.start()
                processes.append(p)
                time.sleep(delay)
        for p in processes:
            p.join()
        if step:
            for i in range(threads):
                p = Process(target=main, args=["step",],
                            kwargs={"loglevel": loglevel})
                p.start()
                processes.append(p)
                time.sleep(delay)
        for p in processes:
            p.join()


@cc.command(name="stack2")
@click.pass_context
@click.option('-r', '--ref', is_flag=True, help='Compute the REF Stack')
@click.option('-m', '--mov', is_flag=True, help='Compute the MOV Stacks')
@click.option('-s', '--step', is_flag=True, help='Compute the STEP Stacks')
def cc_stack2(ctx, ref, mov, step):
    """Stacks the [REF] or [MOV] windows.
    Computes the STACK jobs.
    """
    click.secho('Lets STACK !', fg='green')
    from ..s04_stack2 import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']

    if ref and mov:
        click.secho("With MSNoise 1.6, you can't run REF & MOV stacks"
                    "simultaneously, please run them one after the other.")
        sys.exit()

    if threads == 1:
        if ref:
            main('ref', loglevel=loglevel)
        if mov:
            main('mov', loglevel=loglevel)
        if step:
            main('step', loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        if ref:
            for i in range(threads):
                p = Process(target=main, args=["ref",],
                            kwargs={"loglevel": loglevel})
                p.start()
                processes.append(p)
                time.sleep(delay)
        for p in processes:
            p.join()
        if mov:
            for i in range(threads):
                p = Process(target=main, args=["mov",],
                            kwargs={"loglevel": loglevel})
                p.start()
                processes.append(p)
                time.sleep(delay)
        for p in processes:
            p.join()
        if step:
            for i in range(threads):
                p = Process(target=main, args=["step",],
                            kwargs={"loglevel": loglevel})
                p.start()
                processes.append(p)
                time.sleep(delay)
        for p in processes:
            p.join()


@cc.group(name="plot")
def cc_plot():
    """Commands to trigger different plots"""
    pass


@cc_plot.command(name="distance",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-a', '--ampli', default=1.0, help='Amplification')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
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

    if ctx.obj['MSNOISE_custom']:
        from distance import main
    else:
        from ..plots.distance import main
    main(filterid, comp, ampli, show, outfile, refilter, virtual_source,
         **extra_args)


@cc_plot.command(name="interferogram",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.argument('sta1')
@click.argument('sta2')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1,
              help='Mov Stack to read from disk')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.option('-r', '--refilter', default=None,
              help='Refilter CCFs before plotting (e.g. 4:8 for filtering CCFs '
                   'between 4.0 and 8.0 Hz. This will update the plot title.')
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED,
                callback=parse_extra_args)
@click.pass_context
def cc_plot_interferogram(ctx, sta1, sta2, filterid, comp, mov_stack, show,
                          outfile,
                          refilter, extra_args):
    """Plots the interferogram between sta1 and sta2 (parses the CCFs)
    STA1 and STA2 must be provided with this format: NET.STA !"""

    if sta1 > sta2:
        click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
        return
    if ctx.obj['MSNOISE_custom']:
        from interferogram import main
    else:
        from ..plots.interferogram import main
    main(sta1, sta2, filterid, comp, mov_stack, show, outfile, refilter,
         **extra_args)


@cc_plot.command(name="ccftime",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.argument('sta1')
@click.argument('sta2')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1,
              help='Mov Stack to read from disk')
@click.option('-a', '--ampli', default=5.0, help='Amplification')
@click.option('-S', '--seismic', is_flag=True, help='Seismic style')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.option('-e', '--envelope', is_flag=True, help='Plot envelope instead of '
                                                     'time series')
@click.option('-r', '--refilter', default=None,
              help='Refilter CCFs before plotting (e.g. 4:8 for filtering CCFs '
                   'between 4.0 and 8.0 Hz. This will update the plot title.')
@click.option("--normalize", default="individual")
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED,
                callback=parse_extra_args)
@click.pass_context
def cc_plot_ccftime(ctx, sta1, sta2, filterid, comp, mov_stack,
                    ampli, seismic, show, outfile, envelope, refilter,
                    normalize, extra_args):
    """Plots the ccf vs time between sta1 and sta2
    STA1 and STA2 must be provided with this format: NET.STA !"""

    if sta1 > sta2:
        click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
        return
    if ctx.obj['MSNOISE_custom']:
        from ccftime import main
    else:
        from ..plots.ccftime import main
    main(sta1, sta2, filterid, comp, mov_stack, ampli, seismic, show, outfile,
         envelope, refilter, normalize, **extra_args)


@cc_plot.command(name="spectime",
                 context_settings=dict(ignore_unknown_options=True, ))
@click.argument('sta1')
@click.argument('sta2')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1,
              help='Mov Stack to read from disk')
@click.option('-a', '--ampli', default=5.0, help='Amplification')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.option('-r', '--refilter', default=None,
              help='Refilter CCFs before plotting (e.g. 4:8 for filtering CCFs '
                   'between 4.0 and 8.0 Hz. This will update the plot title.')
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED,
                callback=parse_extra_args)
@click.pass_context
def cc_plot_spectime(ctx, sta1, sta2, filterid, comp, mov_stack,
                     ampli, show, outfile, refilter, extra_args):
    """Plots the ccf's spectrum vs time between sta1 and sta2
    STA1 and STA2 must be provided with this format: NET.STA !"""

    if sta1 > sta2:
        click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
        return
    if ctx.obj['MSNOISE_custom']:
        from spectime import main
    else:
        from ..plots.spectime import main
    main(sta1, sta2, filterid, comp, mov_stack, ampli, show, outfile,
         refilter, **extra_args)


@cc.group(cls=OrderedGroup)
def dvv():
    """Commands for the "Relative Velocity Variations" Workflow"""
    pass


@dvv.command(name='compute_mwcs')
@click.pass_context
def dvv_compute_mwcs(ctx):
    """Computes the MWCS jobs"""
    from ..s05compute_mwcs import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main, kwargs={"loglevel": loglevel})
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()

@dvv.command(name='compute_mwcs2')
@click.pass_context
def dvv_compute_mwcs2(ctx):
    """Computes the MWCS jobs"""
    from ..s05compute_mwcs2 import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main, kwargs={"loglevel": loglevel})
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


@dvv.command(name='compute_stretching')
@click.pass_context
def dvv_compute_stretching(ctx):
    """[experimental] Computes the stretching based on the new stacked data"""
    from ..stretch import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main()
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main)
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


@dvv.command(name='compute_dtt')
@click.pass_context
@click.option('-i', '--interval', default=1.0, help='Number of days before now to\
 search for modified Jobs')
def dvv_compute_dtt(ctx, interval):
    """Computes the dt/t jobs based on the new MWCS data"""
    from ..s06compute_dtt import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main, kwargs={"loglevel": loglevel})
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


@dvv.command(name='compute_dtt2')
@click.pass_context
@click.option('-i', '--interval', default=1.0, help='Number of days before now to\
 search for modified Jobs')
def dvv_compute_dtt2(ctx, interval):
    """Computes the dt/t jobs based on the new MWCS data"""
    from ..s06compute_dtt2 import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        for i in range(threads):
            p = Process(target=main, kwargs={"loglevel": loglevel})
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()

@dvv.group(name="plot")
def dvv_plot():
    """Commands to trigger different plots"""
    pass


@dvv_plot.command(name='mwcs')
@click.argument('sta1')
@click.argument('sta2')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1,
              help='Mov Stack to read from disk')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def dvv_plot_mwcs(ctx, sta1, sta2, filterid, comp, mov_stack, show, outfile):
    """Plots the mwcs results between sta1 and sta2 (parses the CCFs)
    STA1 and STA2 must be provided with this format: NET.STA !"""
    if sta1 > sta2:
        click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
        return
    if ctx.obj['MSNOISE_custom']:
        from mwcs import main
    else:
        from ..plots.mwcs import main
    main(sta1, sta2, filterid, comp, mov_stack, show, outfile)


@dvv_plot.command(name="dvv")
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=0, help='Plot specific mov stacks')
@click.option('-p', '--pair', default=None, help='Plot a specific pair',
              multiple=True)
@click.option('-A', '--all', help='Show the ALL line?', is_flag=True)
@click.option('-M', '--dttname', default="M", help='Plot M or M0?')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def dvv_plot_dvv(ctx, mov_stack, comp, dttname, filterid, pair, all, show, outfile):
    """Plots the dv/v (parses the dt/t results)
    Individual pairs can be plotted extra using the -p flag one or more times.
    Example: msnoise plot dvv -p ID_KWUI_ID_POSI
    Example: msnoise plot dvv -p ID_KWUI_ID_POSI -p ID_KWUI_ID_TRWI
    Remember to order stations alphabetically !
    """
    if ctx.obj['MSNOISE_custom']:
        from dvv import main
    else:
        from ..plots.dvv import main
    main(mov_stack, dttname, comp, filterid, pair, all, show, outfile)




@dvv_plot.command(name="dtt")
@click.argument('sta1')
@click.argument('sta2')
@click.argument('day')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1,
              help='Mov Stack to read from disk')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def dvv_plot_dtt(ctx, sta1, sta2, filterid, day, comp, mov_stack, show, outfile):
    """Plots a graph of dt against t
    STA1 and STA2 must be provided with this format: NET.STA !
    DAY must be provided in the ISO format: YYYY-MM-DD"""
    if sta1 > sta2:
        click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
        return
    if ctx.obj['MSNOISE_custom']:
        from dtt import main
    else:
        from ..plots.dtt import main
    main(sta1, sta2, filterid, comp, day, mov_stack, show, outfile)



@dvv_plot.command(name="timing")
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=0, help='Plot specific mov stacks')
@click.option('-p', '--pair', default=None, help='Plot a specific pair',
              multiple=True)
@click.option('-A', '--all', help='Show the ALL line?', is_flag=True)
@click.option('-M', '--dttname', default="A", help='Plot M or M0?')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def dvv_plot_timing(ctx, mov_stack, comp, dttname, filterid, pair, all, show, outfile):
    """Plots the timing (parses the dt/t results)
    Individual pairs can be plotted extra using the -p flag one or more times.
    Example: msnoise plot timing -p ID_KWUI_ID_POSI
    Example: msnoise plot timing -p ID_KWUI_ID_POSI -p ID_KWUI_ID_TRWI
    Remember to order stations alphabetically !
    """
    if ctx.obj['MSNOISE_custom']:
        from timing import main
    else:
        from ..plots.timing import main
    main(mov_stack, dttname, comp, filterid, pair, all, show, outfile)


# @cli.command(name='compute_dvv')
# @click.option('-f', '--filterid', default=1, help='Filter ID')
# @click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
# @click.option('-m', '--mov_stack', default=0, help='Plot specific mov stacks')
# @click.option('-p', '--pair', default=None, help='Plot a specific pair',
#               multiple=True)
# @click.option('-A', '--all', help='Show the ALL line?', is_flag=True)
# @click.option('-M', '--dttname', default="M", help='Plot M or M0?')
# @click.option('-s', '--show', help='Show interactively?',
#               default=True, type=bool)
# @click.option('-o', '--outfile', help='Output filename (?=auto)',
#               default=None, type=str)
# @click.pass_context
# def compute_dvv(ctx, mov_stack, comp, dttname, filterid, pair, all, show, outfile):
#     """Plots the dv/v (parses the dt/t results)
#     Individual pairs can be plotted extra using the -p flag one or more times.
#     Example: msnoise plot dvv -p ID_KWUI_ID_POSI
#     Example: msnoise plot dvv -p ID_KWUI_ID_POSI -p ID_KWUI_ID_TRWI
#     Remember to order stations alphabetically !
#     """
#     if ctx.obj['MSNOISE_custom']:
#         from s07_compute_dvv import main
#     else:com
#         from ..s07_compute_dvv import main
#     main(mov_stack, dttname, comp, filterid, pair, all, show, outfile)



#
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
    """Computes the PSD jobs, based on New or Modified files identified by
       the new_jobs step"""
    from ..ppsd_compute import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    print(loglevel)
    if threads == 1:
        main(loglevel=loglevel, njobs_per_worker=njobs_per_worker)
    else:
        from multiprocessing import Process
        processes = []
        kwargs = {"loglevel": loglevel, "njobs_per_worker": njobs_per_worker}
        for i in range(threads):
            p = Process(target=main, kwargs=kwargs)
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


@qc.command(name='psd_to_hdf')
@click.option('-n', '--njobs_per_worker', default=9999,
              help='Reduce this number when processing a small number of days '
                   'but a large number of stations')
@click.pass_context
def qc_psd_to_hdf(ctx, njobs_per_worker):
    """Groups the PSD calculated as NPZ to HDF"""
    from ..psd_to_hdf import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    print(loglevel)
    if threads == 1:
        main(loglevel=loglevel, njobs_per_worker=njobs_per_worker)
    else:
        from multiprocessing import Process
        processes = []
        kwargs = {"loglevel": loglevel, "njobs_per_worker": njobs_per_worker}
        for i in range(threads):
            p = Process(target=main, kwargs=kwargs)
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


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


@qc.command(name='hdf_to_rms')
@click.pass_context
def qc_compute_rms(ctx):
    """Computes the RMS based on HDFs"""
    from ..psd_compute_rms import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    print(loglevel)
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        kwargs = {"loglevel": loglevel}
        for i in range(threads):
            p = Process(target=main, kwargs=kwargs)
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


@qc.command(name='export_rms')
@click.pass_context
def qc_export_rms(ctx):
    """Exports the RMS dataframes as CSV files"""
    from ..psd_export_rms import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
    print(loglevel)
    if threads == 1:
        main(loglevel=loglevel)
    else:
        from multiprocessing import Process
        processes = []
        kwargs = {"loglevel": loglevel}
        for i in range(threads):
            p = Process(target=main, kwargs=kwargs)
            p.start()
            processes.append(p)
            time.sleep(delay)
        for p in processes:
            p.join()


@qc.command(name='optimize')
@click.pass_context
def qc_psd_optimize(ctx):
    """Optimizes the HDFs using ptrepack (should be used periodically)"""
    import os, glob
    for file in sorted(glob.glob("PSD/HDF/*")):
        print("Optimizing %s" % file)
        os.system("ptrepack --chunkshape=auto --propindexes --complevel=9 --complib=blosc %s %s" % (file, file.replace(".h5", '_r.h5')))
        os.system("mv %s %s " % ( file.replace(".h5", '_r.h5'), file,) )


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
    import pkg_resources
    plugins = plugins.split(",")
    for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.commands'):
        module_name = ep.module_name.split(".")[0]
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
@click.option('-c', '--content', default=False, is_flag=True)
def utils_test(prefix, content):
    """Runs the test suite, should be executed in an empty folder!"""
    import matplotlib.pyplot as plt
    plt.switch_backend("agg")
    if not content:
        from ..test.tests import main
    else:
        from ..test.content_tests import main
    main(prefix=prefix)


@utils.command(name="jupyter")
def utils_jupyter():
    """Launches an jupyter notebook in the current folder"""
    os.system("jupyter notebook --ip 0.0.0.0 --no-browser")



def run():
    try:
        cli(obj={})
    except MSNoiseError as e:
        logging.critical(str(e))
