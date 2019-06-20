import traceback
import logging
import os
import sys
import sqlalchemy
import time

import click
import pkg_resources

from .. import MSNoiseError, DBConfigNotFoundError
from ..api import connect, get_config, update_station, get_logger, get_job_types
from ..msnoise_table_def import DataAvailability

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
            default_value = default[key][1]
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
    click.echo('\nDatabase information stored in the db.ini file:')
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
    click.echo('Configuration values:\n'
            '   | Normal colour indicates that the default value is used\n'
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
    click.echo("\nJobs:")
    for jobtype in ["CC", "STACK", "MWCS", "DTT"]:
        click.echo(' %s:' % jobtype)
        n = None
        for (n, jobtype) in get_job_types(db, jobtype):
            click.echo("  %s : %i" % (jobtype, n))
        if n is None:
            click.echo('  none defined')


def info_plugins(db):
    """
    Show information about configured plugins.
    """
    from ..api import get_config, get_job_types
    plugins = get_config(db, "plugins")
    if not plugins:
        return
    plugins = plugins.split(",")
    for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.jobtypes'):
        module_name = ep.module_name.split(".")[0]
        if module_name in plugins:
            click.echo('')
            click.echo('Plugin: %s' % module_name)
            for row in ep.load()():
                click.echo(' %s:' % row["name"])
                for (n, jobtype) in get_job_types(db, row["name"]):
                    click.echo("  %s : %i" % (jobtype, n))


@click.group()
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


# @with_plugins(iter_entry_points('msnoise.plugins'))
@cli.group()
def plugin():
    """Runs a command in a named plugin"""
    pass


# @with_plugins(iter_entry_points('msnoise.plugins'))
@cli.group()
def p():
    """Short cut for plugins"""
    pass


@cli.command()
@click.option('-p', '--prefix', default="", help='Prefix for tables')
def test(prefix):
    """Runs the test suite, should be executed in an empty folder!"""
    import matplotlib.pyplot as plt
    plt.switch_backend("agg")
    from ..test.tests import main
    main(prefix=prefix)


@cli.command()
@click.option('-p', '--port', default=5000, help='Port to open')
def admin(port):
    """Starts the Web Admin on http://localhost:5000 by default"""
    from ..msnoise_admin import main
    main(port)


@cli.command(name='upgrade-db')
def upgrade_db():
    """DEPRECATED: since MSNoise 1.6, please use "msnoise db upgrade" instead"""
    click.echo(upgrade_db.__doc__)
    pass


@cli.group()
def db():
    """Top level command to interact with the database"""
    pass


@db.command(name='clean_duplicates')
def clean_duplicates():
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


@db.command()
def upgrade():
    """Upgrade the database from previous to a new version.\n
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
            db.add(Config(name=name, value=default[name][1]))
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
                   "file, net, sta, comp)" %
                   prefix)
        db.commit()
    except:
        logging.info("It looks like the v1.5 'da_index' is already in the DB")
        db.rollback()

    db.close()


@db.command()
@click.option('--tech', help='Database technology: 1=SQLite 2=MySQL',
              default=None)
def init(tech):
    """This command initializes the current folder to be a MSNoise Project
    by creating a database and a db.ini file."""
    click.echo('Launching the init')
    from ..s000installer import main
    main(tech)


@db.command()
@click.argument('sql_command')
def execute(sql_command):
    """EXPERT MODE: Executes 'sql_command' on the database. Use this command
    at your own risk!!"""
    from msnoise.api import connect

    db = connect()
    r = db.execute(sql_command)

    if sql_command.count("select") or sql_command.count("SELECT"):
        print(r.fetchall())
    db.commit()
    db.close()


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


@cli.command()
def install():
    """DEPRECATED: since MSNoise 1.6, please use "msnoise db init" instead"""
    click.echo(install.__doc__)
    pass


@cli.group()
def config():
    """
    This command allows to set a parameter value in the database, show
    parameter values, or synchronise station metadata, depending on the
    invoked subcommands.

    Called without argument, it used to launch the Configurator (now
    accessible using 'msnoise config gui') but the recommended way to
    configure MSNoise is now to use the web interface through the command
    'msnoise admin'.
    """
    pass


@config.command(name='gui')
def config_gui():
    """
    Run the deprecated configuration GUI tool.  Please use the configuration
    web interface using 'msnoise admin' instead.
    """
    from ..s001configurator import main
    click.echo("Let's Configure MSNoise !")
    main()


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
@click.option('-s', '--sys', is_flag=True, help='System Info')
@click.option('-m', '--modules', is_flag=True, help='Modules Info')
@click.option('-e', '--env', is_flag=True, help='Environment Info')
@click.option('-a', '--all', is_flag=True, help='All Info')
@click.pass_context
def bugreport(ctx, sys, modules, env, all):
    """This command launches the Bug Report script."""
    click.echo('Let\'s Bug Report MSNoise !')
    # click.echo('Working on %i threads' % ctx.obj['MSNOISE_threads'])
    from ..bugreport import main
    main(sys, modules, env, all)


@cli.command()
@click.option('--fromDA',  help='Populates the station table '
                                'using network and station codes'
                                ' found in the data_availability'
                                ' table, overrides the default'
                                ' workflow step.',
              is_flag=True)
def populate(fromda):
    """Rapidly scan the archive filenames and find Network/Stations"""
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
    """Determines if new CC jobs are to be defined"""
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


@cli.command(name='compute_cc')
@click.pass_context
def compute_cc(ctx):
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


@cli.command(name='compute_cc2')
@click.pass_context
def compute_cc2(ctx):
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


@cli.command()
@click.pass_context
@click.option('-r', '--ref', is_flag=True, help='Compute the REF Stack')
@click.option('-m', '--mov', is_flag=True, help='Compute the MOV Stacks')
@click.option('-s', '--step', is_flag=True, help='Compute the STEP Stacks')
def stack(ctx, ref, mov, step):
    """Stacks the [REF] and/or [MOV] windows"""
    click.secho('Lets STACK !', fg='green')
    from ..s04stack import main
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    loglevel = ctx.obj['MSNOISE_verbosity']
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


@cli.command(name='compute_mwcs')
@click.pass_context
def compute_mwcs(ctx):
    """Computes the MWCS based on the new stacked data"""
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


@cli.command(name='compute_stretching')
def compute_stretching():
    """[experimental] Computes the stretching based on the new stacked data"""
    from ..stretch import main
    main()


@cli.command(name='compute_dtt')
@click.pass_context
@click.option('-i', '--interval', default=1.0, help='Number of days before now to\
 search for modified Jobs')
def compute_dtt(ctx, interval):
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



@cli.command(name='compute_dvv')
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
def compute_dvv(ctx, mov_stack, comp, dttname, filterid, pair, all, show, outfile):
    """Plots the dv/v (parses the dt/t results)\n
    Individual pairs can be plotted extra using the -p flag one or more times.\n
    Example: msnoise plot dvv -p ID_KWUI_ID_POSI\n
    Example: msnoise plot dvv -p ID_KWUI_ID_POSI -p ID_KWUI_ID_TRWI\n
    Remember to order stations alphabetically !
    """
    if ctx.obj['MSNOISE_custom']:
        from s07_compute_dvv import main
    else:
        from ..s07_compute_dvv import main
    main(mov_stack, dttname, comp, filterid, pair, all, show, outfile)


@cli.command()
@click.argument('jobtype')
@click.option('-a', '--all', is_flag=True, help='Reset all jobs')
@click.option('-r', '--rule', help='Reset job that match this SQL rule')
def reset(jobtype, all, rule):
    """Resets the job to "T"odo. ARG is [CC] or [DTT]. By default
    only resets jobs "I"n progress. --all resets all jobs, whatever
    the flag value"""
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
def jupyter():
    """Launches an jupyter notebook in the current folder"""
    os.system("jupyter notebook --ip 0.0.0.0 --no-browser")


#
# PLOT GROUP
#

@cli.group()
def plot():
    """Top level command to trigger different plots"""
    pass


@plot.command(name='data_availability')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def data_availability(ctx, show, outfile):
    """Plots the Data Availability vs time"""
    if ctx.obj['MSNOISE_custom']:
        from data_availability import main
    else:
        from ..plots.data_availability import main
    main(show, outfile)


@plot.command()
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
def dvv(ctx, mov_stack, comp, dttname, filterid, pair, all, show, outfile):
    """Plots the dv/v (parses the dt/t results)\n
    Individual pairs can be plotted extra using the -p flag one or more times.\n
    Example: msnoise plot dvv -p ID_KWUI_ID_POSI\n
    Example: msnoise plot dvv -p ID_KWUI_ID_POSI -p ID_KWUI_ID_TRWI\n
    Remember to order stations alphabetically !
    """
    if ctx.obj['MSNOISE_custom']:
        from dvv import main
    else:
        from ..plots.dvv import main
    main(mov_stack, dttname, comp, filterid, pair, all, show, outfile)


@plot.command()
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
def timing(ctx, mov_stack, comp, dttname, filterid, pair, all, show, outfile):
    """Plots the timing (parses the dt/t results)\n
    Individual pairs can be plotted extra using the -p flag one or more times.\n
    Example: msnoise plot timing -p ID_KWUI_ID_POSI\n
    Example: msnoise plot timing -p ID_KWUI_ID_POSI -p ID_KWUI_ID_TRWI\n
    Remember to order stations alphabetically !
    """
    if ctx.obj['MSNOISE_custom']:
        from timing import main
    else:
        from ..plots.timing import main
    main(mov_stack, dttname, comp, filterid, pair, all, show, outfile)


@plot.command(context_settings=dict(ignore_unknown_options=True,))
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
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED, callback=parse_extra_args)
@click.pass_context
def interferogram(ctx, sta1, sta2, filterid, comp, mov_stack, show, outfile,
                  refilter, extra_args):
    """Plots the interferogram between sta1 and sta2 (parses the CCFs)\n
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



@plot.command(context_settings=dict(ignore_unknown_options=True,))
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
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED, callback=parse_extra_args)
@click.pass_context
def ccftime(ctx, sta1, sta2, filterid, comp, mov_stack,
            ampli, seismic, show, outfile, envelope, refilter, normalize, extra_args):
    """Plots the ccf vs time between sta1 and sta2\n
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


@plot.command(context_settings=dict(ignore_unknown_options=True,))
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
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED, callback=parse_extra_args)
@click.pass_context
def spectime(ctx, sta1, sta2, filterid, comp, mov_stack,
            ampli, show, outfile, refilter, extra_args):
    """Plots the ccf's spectrum vs time between sta1 and sta2\n
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


@plot.command()
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
def mwcs(ctx, sta1, sta2, filterid, comp, mov_stack, show, outfile):
    """Plots the mwcs results between sta1 and sta2 (parses the CCFs)\n
    STA1 and STA2 must be provided with this format: NET.STA !"""
    if sta1 > sta2:
        click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
        return
    if ctx.obj['MSNOISE_custom']:
        from mwcs import main
    else:
        from ..plots.mwcs import main
    main(sta1, sta2, filterid, comp, mov_stack, show, outfile)


@plot.command(context_settings=dict(ignore_unknown_options=True,))
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
@click.argument('extra_args', nargs=-1, type=click.UNPROCESSED, callback=parse_extra_args)
@click.pass_context
def distance(ctx, filterid, comp, ampli, show, outfile, refilter,
             virtual_source, extra_args):
    """Plots the REFs of all pairs vs distance"""
    
        
        
    if ctx.obj['MSNOISE_custom']:
        from distance import main
    else:
        from ..plots.distance import main
    main(filterid, comp, ampli, show, outfile, refilter, virtual_source,
         **extra_args)


@plot.command(name='station_map')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def station_map(ctx, show, outfile):
    """Plots the station map (very very basic)"""
    if ctx.obj['MSNOISE_custom']:
        from station_map import main
    else:
        from ..plots.station_map import main
    main(show, outfile)


@plot.command()
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
def dtt(ctx, sta1, sta2, filterid, day, comp, mov_stack, show, outfile):
    """Plots a graph of dt against t\n
    STA1 and STA2 must be provided with this format: NET.STA !\n
    DAY must be provided in the ISO format: YYYY-MM-DD"""
    if sta1 > sta2:
        click.echo("Stations STA1 and STA2 must be sorted alphabetically.")
        return
    if ctx.obj['MSNOISE_custom']:
        from dtt import main
    else:
        from ..plots.dtt import main
    main(sta1, sta2, filterid, comp, day, mov_stack, show, outfile)


## Main script

try:
    db = connect()
    plugins = get_config(db, "plugins")
    db.close()
except DBConfigNotFoundError:
    plugins = None
except sqlalchemy.exc.OperationalError as e:
    logging.critical('Unable to read project configuration: error connecting to the database:\n{}'.format(str(e)))
    sys.exit(1)

if plugins:
    plugins = plugins.split(",")
    for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.commands'):
        module_name = ep.module_name.split(".")[0]
        if module_name in plugins:
            plugin.add_command(ep.load())
            p.add_command(ep.load())


def run():
    try:
        cli(obj={})
    except MSNoiseError as e:
        logging.critical(str(e))
