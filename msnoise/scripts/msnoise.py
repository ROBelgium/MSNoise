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
from ..core.db import connect, get_logger
from ..core.config import (delete_config_set, get_config, get_config_set_details, list_config_sets)
from ..core.stations import update_station

from ..msnoise_table_def import DataAvailability

from ..msnoise_table_def import Config


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


def show_config_values(db):
    """Print global config values; green = modified from default."""
    from ..msnoise_table_def import Config as _Config
    rows = (db.query(_Config)
              .filter(_Config.category == "global")
              .filter(_Config.set_number == 1)
              .order_by(_Config.name)
              .all())
    for row in rows:
        value = row.value if row.value is not None else ""
        default_value = row.default_value if row.default_value is not None else ""
        display_value = value if value != "" else "''"
        modified = value != default_value
        line = (" M " if modified else "   ") + f"{row.name}: {display_value}"
        click.secho(line, fg="green" if modified else None)


def info_db_ini():
    """
    Show information stored in the db.ini file.
    """
    from ..core.db import read_db_inifile
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
    from ..core.config import get_config
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

    from ..msnoise_table_def import DataSource
    all_ds = db.query(DataSource).order_by(DataSource.ref).all()
    for _ds in all_ds:
        data_folder = _ds.uri or ""
        if data_folder and data_folder.startswith("fdsn://"):
            click.echo(" - DataSource '%s': FDSN → %s" % (_ds.name, data_folder[7:]))
        elif data_folder and data_folder.startswith("eida://"):
            click.echo(" - DataSource '%s': EIDA → %s" % (_ds.name, data_folder[7:]))
        elif data_folder and os.path.isdir(data_folder):
            click.echo(" - DataSource '%s': %s exists" % (_ds.name, data_folder))
        elif data_folder:
            click.secho(" - DataSource '%s': %s does not exist !" % (_ds.name, data_folder), fg='red')
        else:
            click.secho(" - DataSource '%s': no URI configured" % _ds.name, fg='yellow')

    output_folder = get_config(db, "output_folder")
    if os.path.isdir(output_folder):
        click.echo(" - %s exists" % output_folder)
    else:
        if get_config(db, 'keep_all') in ['Y', 'y']:
            from ..core.workflow import get_workflow_job_counts
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
    """Show global config values. Use 'msnoise config list' for full per-category detail."""
    click.echo('')
    click.echo('Global config  (green = modified | run "msnoise config list" for full detail)')
    show_config_values(db)


def info_stations(db):
    """
    Show information about configured stations.
    """
    from ..core.stations import get_stations
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
    from ..core.workflow import get_workflow_job_counts
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
@click.option('--from-yaml', 'from_yaml', default=None, metavar='PATH',
              help='Seed config sets and links from a project YAML (category_N + after keys).')
@click.option('--yes', '-y', is_flag=True, default=False,
              help='Skip interactive prompts, accepting all defaults (e.g. auto-create jobs).')
def db_init(tech, auto_workflow, from_yaml, yes):
    """This command initializes the current folder to be a MSNoise Project
    by creating a database and a db.ini file."""
    click.echo('Launching the init')
    from ..s00_installer import main
    result = main(tech)
    if result != 0:
        return

    if from_yaml:
        import os
        import datetime
        if not os.path.isfile(from_yaml):
            click.echo(f"Error: {from_yaml!r} not found.", err=True)
            return
        click.echo(f'\nSeeding project config from {from_yaml!r}...')
        from ..core.db import connect
        from ..core.config import create_project_from_yaml, get_config
        from ..core.fdsn import is_remote_source
        from ..msnoise_table_def import DataSource, Station
        from ..core.stations import resolve_data_source, get_stations
        from ..core.workflow import get_workflow_steps
        from ..s02_new_jobs import update_job
        db = connect()
        try:
            created, warnings = create_project_from_yaml(db, from_yaml)
        except ValueError as e:
            click.echo(f"Error: {e}", err=True)
            db.close()
            return
        click.echo(f'Created: {", ".join(created)}')
        for w in warnings:
            click.echo(f'Warning: {w}', err=True)
        click.echo('Setup complete!' if not warnings
                   else 'Setup complete (with warnings — check links in admin UI).')

        _DEFAULT_START, _DEFAULT_END = "1970-01-01", "2100-01-01"
        startdate = get_config(db, "startdate")
        enddate   = get_config(db, "enddate")
        all_ds = db.query(DataSource).all()
        remote_ds = [d for d in all_ds if is_remote_source(d.uri or "")]
        if remote_ds:
            ds_names = ", ".join(d.name for d in remote_ds)
            click.echo(f"\nDataSource(s) '{ds_names}' use remote FDSN/EIDA.")
            if startdate == _DEFAULT_START or enddate == _DEFAULT_END:
                click.echo(
                    "  ⚠  startdate/enddate are still at default values "
                    f"({startdate} → {enddate}).\n"
                    "  Set them in your project.yaml or via `msnoise config set` "
                    "before downloading or creating jobs.",
                    err=True,
                )
            else:
                click.echo("  How would you like to proceed?")
                click.echo("  [1] Bulk download first  — fetch all raw waveforms into SDS, then run pipeline")
                click.echo("  [2] Stream from FDSN     — preprocess fetches on-the-fly (creates jobs now)")
                click.echo("  [3] Skip                 — I'll handle it manually")
                choice = click.prompt("Choice", type=click.Choice(["1", "2", "3"]),
                                      default="1", show_default=True)

                if choice == "1":
                    click.echo("\n  Next steps:")
                    click.echo(f"    msnoise utils download          # fetch raw waveforms into SDS")
                    click.echo(f"    msnoise scan_archive            # populate data availability")
                    click.echo(f"    msnoise utils run_workflow      # run full pipeline")

                elif choice == "2":
                    pre_steps = [s for s in get_workflow_steps(db)
                                 if s.category == "preprocess"]
                    if not pre_steps:
                        click.echo("  No preprocess step found — skipping job creation.", err=True)
                    else:
                        pre_step = pre_steps[0]
                        start = datetime.date.fromisoformat(startdate)
                        end   = datetime.date.fromisoformat(enddate)
                        dates, cur = [], start
                        while cur <= end:
                            dates.append(cur.isoformat())
                            cur += datetime.timedelta(days=1)
                        from ..core.workflow import massive_insert_job
                        now = datetime.datetime.utcnow()
                        jobs = [
                            {"day": day, "pair": f"{sta.net}.{sta.sta}.{loc}",
                             "jobtype": pre_step.step_name, "flag": "T",
                             "step_id": pre_step.step_id,
                             "lineage": pre_step.step_name, "lastmod": now}
                            for day in dates
                            for sta in get_stations(db, all=False)
                            for loc in (sta.locs() or [""])
                        ]
                        massive_insert_job(db, jobs)
                        n = len(jobs)
                        click.echo(f"  Created {n} preprocess job(s) over {len(dates)} day(s).")
                        click.echo(f"\n  Next step:")
                        click.echo(f"    msnoise utils run_workflow")

        db.close()
        return

    if auto_workflow or click.confirm(
        '\nWould you like to automatically create all default config sets, '
        'workflow steps and links?',
        default=True
    ):
        _create_default_workflow()


def _create_default_workflow():
    """Create default config sets, workflow steps and links after db init."""
    from ..core.db import connect
    from ..core.config import create_config_set
    from ..core.workflow import (create_workflow_steps_from_config_sets,
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
    from msnoise.core.db import connect
    from msnoise.core.stations import get_stations
    session = connect()
    stations = get_stations(session)
    for sta in stations:
        data = session.query(DataAvailability.loc, DataAvailability.chan). \
            filter(text("net=:net")). \
            filter(text("sta=:sta")). \
            group_by(DataAvailability.loc, DataAvailability.chan). \
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
    from msnoise.core.db import connect
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
    """Upgrade the database from a previous version.

    Ensures every parameter defined in the config CSV files is present in the
    database with its default value.  Covers all categories (global, cc, mwcs,
    psd, …) for every config set already in the DB — not just global params.

    Safe to run on an already up-to-date project: existing values are never
    overwritten.
    """
    import csv
    import glob
    import os
    from ..core.db import connect, read_db_inifile

    db = connect()
    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''

    config_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "config"
    )
    csv_files = sorted(glob.glob(os.path.join(config_dir, "config_*.csv")))

    added = 0
    for csv_path in csv_files:
        # Derive category name from filename: config_cc.csv → "cc"
        category = os.path.basename(csv_path)[len("config_"):-len(".csv")]

        # Find every set_number that already exists for this category in the DB
        existing_sets = [
            row[0] for row in
            db.query(Config.set_number)
              .filter(Config.category == category)
              .distinct()
              .all()
        ]
        # Global always has set 1; if no sets exist yet for other categories,
        # there's nothing to upgrade (create_set handles first-time creation).
        if not existing_sets:
            if category == "global":
                existing_sets = [1]
            else:
                continue

        with open(csv_path, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))

        for set_number in existing_sets:
            for row in rows:
                name = row["name"]
                default_value = row.get("default", "")
                definition = row.get("definition", "")
                param_type = row.get("type", "str")
                possible_values = row.get("possible_values", "")

                # Only insert if the (name, category, set_number) row is absent
                exists = (
                    db.query(Config)
                      .filter(Config.name == name)
                      .filter(Config.category == category)
                      .filter(Config.set_number == set_number)
                      .first()
                )
                if exists is None:
                    db.add(Config(
                        name=name,
                        category=category,
                        set_number=set_number,
                        value=default_value,
                        default_value=default_value,
                        param_type=param_type,
                        description=definition,
                        possible_values=possible_values,
                        used=True,
                    ))
                    added += 1

        db.commit()

    click.echo(f"✓ DB upgrade complete: {added} new parameter(s) added across "
               f"{len(csv_files)} categories.")

    # Legacy index creation (idempotent — errors are expected on up-to-date DBs)
    for stmt, label in [
        (f"CREATE UNIQUE INDEX job_index ON {prefix}jobs "
         "(day, pair, jobtype)", "v1.5 job_index"),
        (f"CREATE INDEX job_index2 ON {prefix}jobs "
         "(jobtype, flag)", "v1.6 job_index2"),
        (f"CREATE UNIQUE INDEX da_index ON {prefix}data_availability "
         "(path, file, net, sta, loc, chan)", "v1.5 da_index"),
    ]:
        try:
            db.execute(text(stmt))
            db.commit()
        except Exception:
            db.rollback()

    db.close()


@db.command(name='clean_duplicates')
def db_clean_duplicates():
    """Checks the Jobs table and deletes duplicate entries"""
    from msnoise.core.db import connect, read_db_inifile

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


@db.command(name="export-yaml")
@click.argument("output", default="project.yaml", metavar="PATH")
@click.option("--all-values", is_flag=True, default=False,
              help="Write every config key, not just non-default values.")
def db_export_yaml(output, all_values):
    """Export the current project config to a project YAML file.

    The output can be re-applied to a fresh DB with:
    msnoise db init --from-yaml PATH
    """
    from ..core.db import connect
    from ..core.config import export_project_to_yaml
    db = connect()
    path = export_project_to_yaml(db, output, only_non_defaults=not all_values)
    db.close()
    click.echo(f"Project exported to {path!r}")


@db.command(name="dump")
@click.option("--format", default="csv")
def db_dump(format):
    """Dumps the complete database in formatted files, defaults to CSV.
    """
    from ..core.db import get_engine
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
    from ..core.db import get_engine, read_db_inifile
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
    from ..core.db import connect

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
    from ..core.db import connect
    from ..core.stations import get_stations, update_station
    from ..core.signal import preload_instrument_responses

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
@click.argument('key')
@click.argument('value')
@click.pass_context
def config_set(ctx, key, value):
    """Set a configuration parameter using dot notation.

    \b
    KEY format:
      name                   →  global parameter (e.g. output_folder)
      category.name          →  category set 1  (e.g. cc.cc_sampling_rate)
      category.N.name        →  explicit set N  (e.g. mwcs.2.mwcs_wlen)

    \b
    Examples:
      msnoise config set output_folder /data/output
      msnoise config set cc.cc_sampling_rate 25
      msnoise config set cc.2.cc_sampling_rate 25
      msnoise config set mwcs.2.mwcs_wlen 10
      msnoise config set global.hpc Y
    """
    from ..core.config import parse_config_key, _cast_config_value, update_config
    from ..core.db import connect

    try:
        category, set_number, name = parse_config_key(key)
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
        return

    db = connect()
    try:
        # Look up the param_type from the DB for validation
        row = (db.query(Config)
               .filter(Config.name == name)
               .filter(Config.category == category)
               .filter(Config.set_number == set_number)
               .first())
        if row is None:
            click.echo(
                f"Error: parameter '{name}' not found in {category} set {set_number}.\n"
                f"Use 'msnoise config list {category}' to see available parameters.",
                err=True
            )
            ctx.exit(1)
            return

        try:
            validated = _cast_config_value(name, value, row.param_type or "str")
        except ValueError as e:
            click.echo(f"Error: {e}", err=True)
            ctx.exit(1)
            return

        update_config(db, name, validated, category=category, set_number=set_number)
        db.commit()
        click.echo(f"✓  {category}.{set_number}.{name}  =  {validated}")
    finally:
        db.close()


@config.command(name='get')
@click.argument('key')
@click.pass_context
def config_get(ctx, key):
    """Get a configuration parameter value using dot notation.

    \b
    KEY format:
      name                   →  global parameter
      category.name          →  category set 1
      category.N.name        →  explicit set N

    \b
    Examples:
      msnoise config get output_folder
      msnoise config get cc.cc_sampling_rate
      msnoise config get mwcs.2.mwcs_wlen
    """
    from ..core.config import parse_config_key
    from ..core.db import connect

    try:
        category, set_number, name = parse_config_key(key)
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
        return

    db = connect()
    try:
        row = (db.query(Config)
               .filter(Config.name == name)
               .filter(Config.category == category)
               .filter(Config.set_number == set_number)
               .first())
        if row is None:
            click.echo(
                f"Error: parameter '{name}' not found in {category} set {set_number}.",
                err=True
            )
            ctx.exit(1)
            return
        default_marker = " *" if row.value != row.default_value else ""
        click.echo(f"{category}.{set_number}.{name}  =  {row.value}{default_marker}  "
                   f"[{row.param_type or 'str'}]  (default: {row.default_value})")
    finally:
        db.close()


@config.command(name='list')
@click.argument('filter_key', required=False, default=None,
                metavar='[CATEGORY[.N]]')
@click.pass_context
def config_list(ctx, filter_key):
    """List configuration parameters, optionally filtered.

    \b
    Examples:
      msnoise config list              # all categories, all sets
      msnoise config list cc           # all cc sets
      msnoise config list mwcs.2       # mwcs set 2 only

    Non-default values are marked with *.
    """
    from ..core.db import connect

    # Parse the optional filter
    filter_category = None
    filter_set = None
    if filter_key:
        parts = filter_key.split(".")
        filter_category = parts[0]
        if len(parts) >= 2:
            try:
                filter_set = int(parts[1])
            except ValueError:
                click.echo(f"Error: invalid set number in {filter_key!r}", err=True)
                ctx.exit(1)
                return

    db = connect()
    try:
        q = db.query(Config).order_by(
            Config.category, Config.set_number, Config.name)
        if filter_category:
            q = q.filter(Config.category == filter_category)
        if filter_set is not None:
            q = q.filter(Config.set_number == filter_set)

        rows = q.all()
        if not rows:
            click.echo("No configuration parameters found.")
            return

        current_header = None
        for row in rows:
            header = f"{row.category}  (set {row.set_number})"
            if header != current_header:
                if current_header is not None:
                    click.echo("")
                click.echo(f"\n{header}")
                click.echo("─" * len(header))
                current_header = header

            marker = " *" if row.value != row.default_value else "  "
            name_col  = f"{row.name:<40}"
            val_col   = f"{str(row.value):<20}"
            type_col  = f"[{row.param_type or 'str'}]"
            click.echo(f"  {name_col}{val_col}{marker} {type_col}")
    finally:
        db.close()


@config.command(name='reset')
@click.argument('key')
@click.pass_context
def config_reset(ctx, key):
    """Reset a configuration parameter to its default value.

    \b
    KEY format: same dot notation as 'config set'.

    \b
    Examples:
      msnoise config reset cc.cc_sampling_rate
      msnoise config reset mwcs.2.mwcs_wlen
    """
    from ..core.config import parse_config_key, update_config
    from ..core.db import connect

    try:
        category, set_number, name = parse_config_key(key)
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
        ctx.exit(1)
        return

    db = connect()
    try:
        row = (db.query(Config)
               .filter(Config.name == name)
               .filter(Config.category == category)
               .filter(Config.set_number == set_number)
               .first())
        if row is None:
            click.echo(
                f"Error: parameter '{name}' not found in {category} set {set_number}.",
                err=True
            )
            ctx.exit(1)
            return
        update_config(db, name, row.default_value,
                      category=category, set_number=set_number)
        db.commit()
        click.echo(f"✓  {category}.{set_number}.{name}  reset to  {row.default_value}")
    finally:
        db.close()



@config.command()
@click.argument('set_name')
@click.pass_context
def create_set(ctx, set_name):
    """Create a configuration set for a workflow step.

    SET_NAME: Name of the workflow step (e.g., mwcs, mwcs_dtt, etc.)
    """
    from ..core.config import create_config_set
    from ..core.db import connect

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

    from ..core.config import create_config_set
    from ..core.db import connect

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
    from ..core.workflow import create_workflow_step

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
    from ..core.db import connect
    from ..core.workflow import create_workflow_steps_from_config_sets

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
    from ..core.db import connect
    from ..core.workflow import get_workflow_steps

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
    from ..core.db import connect
    from ..core.workflow import get_workflow_graph

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
    from ..core.db import connect
    from ..core.workflow import create_workflow_links_from_steps

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
@click.option('-a', '--all', 'alljobs', is_flag=True,
              help='Reset all jobs regardless of flag (default: only I and F).')
@click.option('-d', '--downstream', is_flag=True,
              help='Also reset all steps reachable downstream from JOBTYPE.')
@click.option('-r', '--rule', help='Reset jobs matching this SQL rule (legacy).')
def reset(jobtype, alljobs, downstream, rule):
    """Reset jobs to "T"odo.

    JOBTYPE can be a step name (e.g. cc_1) or a category (e.g. cc, stack).
    A category resets all active steps of that category.
    --downstream additionally resets every step reachable via WorkflowLinks.

    Examples:

    \b
      msnoise reset cc_1          # reset in-progress cc_1 jobs
      msnoise reset cc            # reset all cc_N steps
      msnoise reset cc_1 -a       # reset all cc_1 jobs (any flag)
      msnoise reset cc_1 -a -d    # reset cc_1 and everything downstream
      msnoise reset cc -a -d      # reset all cc_N and everything downstream
    """
    from ..core.db import connect, read_db_inifile
    from ..core.workflow import reset_jobs
    dbini = read_db_inifile()
    prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
    session = connect()
    if jobtype == "DA":
        session.execute(text("UPDATE {0}data_availability SET flag='M'"
                        .format(prefix)))
        session.commit()
        session.close()
        return
    try:
        reset_steps = reset_jobs(session, jobtype, alljobs=alljobs,
                                 downstream=downstream)
        if downstream:
            click.echo(f"Reset {len(reset_steps)} step(s): {', '.join(reset_steps)}")
        else:
            click.echo(f"Reset step(s): {', '.join(reset_steps)}")
    except ValueError as e:
        click.echo(f"Error: {e}", err=True)
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
    db = connect()
    if fromda:
        logger.info("Populating the Station table")
        logger.info("Overriding workflow...")
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
        from ..core.stations import populate_stations
        populate_stations(db, loglevel=loglevel)


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
    """Plots the station map"""
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
@click.option('--chunk-size', default=0, type=int,
              help='Max pairs to process per worker per day (0 = all, default). '
                   'Set >0 to share a day across multiple parallel workers '
                   'without write conflicts (recommended for >50 stations).')
@click.pass_context
def cc_compute_cc(ctx, chunk_size):
    """Computes the CC jobs (based on the "New Jobs" identified)"""
    from ..s03_compute_no_rotation import main
    run_threaded(main, ctx, chunk_size=chunk_size)


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

    - Absolute date / 1970-01-01  -> Mode A: writes a mean REF file to disk.

    - Negative integer (e.g. -5)  -> Mode B: no file written; rolling reference
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



# DTT-side lineage options (preprocess -> cc -> filter -> stack -> mwcs -> mwcs_dtt)
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
                    refstackid, comp, show,
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
                    refstackid, comp,
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
                    refstackid, comp, ampli, show, outfile, refilter, extra_args):
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
@click.option('--chunk-size', default=0, type=int,
              help='Max stations to process per worker per day (0 = all, default). '
                   'Set >0 to share a day across multiple parallel workers '
                   'without write conflicts (recommended for >50 stations).')
@click.pass_context
def qc_compute_psd(ctx, chunk_size):
    """Computes the PSD jobs, saves results as NetCDF files.
       Based on New or Modified files identified by the new_jobs step."""
    from ..s20_psd_compute import main
    run_threaded(main, ctx, chunk_size=chunk_size)




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
    from ..s21_psd_compute_rms import main
    run_threaded(main, ctx)


@qc.command(name='plot_psd_rms')
@click.argument('seed_ids', nargs=-1, required=False, metavar='[SEED_ID...]')
@click.option('--psd-id', default=1, show_default=True, type=int,
              help='Config-set number for the psd step (matches psd_N in workflow)')
@click.option('--psd-rms-id', default=1, show_default=True, type=int,
              help='Config-set number for the psd_rms step (matches psd_rms_N in workflow)')
@click.option('--type', '-t', 'plot_type', default='timeseries', show_default=True,
              type=click.Choice(['timeseries', 'clockplot', 'hourmap',
                                  'gridmap', 'dailyplot', 'all'],
                                case_sensitive=False),
              help='Plot type to produce')
@click.option('--band', '-b', default=None,
              help='Frequency band label, e.g. "1.0-10.0" (default: first available)')
@click.option('--scale', default=1e9, show_default=True, type=float,
              help='Amplitude scale factor applied before display')
@click.option('--unit', '-u', default='nm', show_default=True,
              help='Amplitude unit label for axis/colorbar')
@click.option('--timezone', '-z', default='UTC', show_default=True,
              help='Time zone for local-time plots (e.g. Europe/Brussels)')
@click.option('--day-start', default=7.0, show_default=True, type=float,
              help='Start of daytime in decimal hours (local time), e.g. 7.0 = 07:00')
@click.option('--day-end', default=19.0, show_default=True, type=float,
              help='End of daytime in decimal hours (local time), e.g. 19.0 = 19:00')
@click.option('--split-date', default=None,
              help='ISO date (YYYY-MM-DD) for before/after split (clockplot, dailyplot)')
@click.option('--annotate', '-a', default=None,
              help='Event annotations as DATE=LABEL,DATE=LABEL')
@click.option('--resample-freq', '-R', default='30min', show_default=True,
              help='Pandas offset string for time-bin resampling '
                   '(e.g. "15min", "1h", "2h").  Default: 30min')
@click.option('--agg-func', '-A', default='mean', show_default=True,
              help='Aggregation method for each resampled bin: '
                   'mean, median, max, min, std, sum, or any pandas agg string.  '
                   'Default: mean')
@click.option('--outfile', '-o', default=None,
              help='Base filename for saved figure(s). Type suffix is appended.')
@click.option('--no-show', is_flag=True, default=False,
              help='Do not call plt.show() (useful in non-interactive environments)')
@click.pass_context
def qc_plot_psd_rms(ctx, seed_ids, psd_id, psd_rms_id, plot_type, band, scale, unit,
                    timezone, day_start, day_end, split_date, annotate,
                    resample_freq, agg_func, outfile, no_show):
    """Plot PSD-RMS results: time-series, clock plots, hour-maps, grid-maps.

    SEED_ID can be zero or more SEED identifiers (NET.STA.LOC.CHAN).  When
    omitted, all stations found on disk for the given psd/psd_rms lineage are
    plotted.  Uses MSNoiseResult for lineage-correct path resolution.

    \b
        msnoise qc plot_psd_rms BE.UCC..HHZ
        msnoise qc plot_psd_rms                           # all stations
        msnoise qc plot_psd_rms BE.UCC..HHZ BE.MEM..HHZ --type timeseries
        msnoise qc plot_psd_rms BE.UCC..HHZ --type clockplot --timezone Europe/Brussels
        msnoise qc plot_psd_rms BE.UCC..HHZ --type all --outfile noise.pdf
        msnoise qc plot_psd_rms --psd-id 2 --psd-rms-id 2 BE.UCC..HHZ
        msnoise qc plot_psd_rms BE.UCC..HHZ --day-start 8 --day-end 20
        msnoise qc plot_psd_rms BE.UCC..HHZ --day-start 0 --day-end 24  # disable
        msnoise qc plot_psd_rms BE.UCC..HHZ --resample-freq 1h --agg-func median
        msnoise qc plot_psd_rms BE.UCC..HHZ -R 15min -A max
    """
    import matplotlib
    if no_show:
        matplotlib.use('Agg')
    from ..plots.psd_rms import main as plot_main

    annotations = {}
    if annotate:
        for pair in annotate.split(','):
            date_str, _, label = pair.partition('=')
            annotations[date_str.strip()] = label.strip()

    loglevel = ctx.obj.get('LOGLEVEL', 'INFO') if ctx.obj else 'INFO'
    plot_main(
        seed_ids=list(seed_ids),
        psd_id=psd_id,
        psd_rms_id=psd_rms_id,
        plot_type=plot_type,
        band=band,
        scale=scale,
        unit=unit,
        time_zone=timezone,
        day_start=day_start,
        day_end=day_end,
        resample_freq=resample_freq,
        agg_func=agg_func,
        split_date=split_date,
        annotations=annotations or None,
        outfile=outfile,
        show=not no_show,
        loglevel=loglevel,
    )




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
@click.option('--tech', default=1, help='Test using (1) SQLite or (2) MariaDB')
@click.option('-c', '--content', default=False, is_flag=True,
              help='Run content tests instead of standard tests')
@click.option('--fast', default=False, is_flag=True,
              help='Run the fast smoke test suite (no real data needed)')
def utils_test(prefix, tech, content, fast):
    """Runs the test suite in a temporary folder.

    Use --fast for a quick smoke test that exercises the full workflow
    lifecycle using stub compute functions (no seismic data required,
    completes in < 30 seconds).
    """
    import matplotlib.pyplot as plt
    import pytest
    plt.switch_backend("agg")

    os.environ["PREFIX"] = prefix
    os.environ["TECH"] = str(tech)

    test_dir = os.path.join(os.path.dirname(__file__), '..', 'test')

    if fast:
        test_path = os.path.join(test_dir, 'test_smoke.py')
        args = ['-s', '-v', '--log-cli-level=WARNING', test_path]
    elif content:
        test_path = os.path.join(test_dir, 'content_tests.py')
        args = ['-s', test_path]
    else:
        test_path = os.path.join(test_dir, 'tests.py')
        args = ['-s', '--log-cli-level=WARNING', test_path]

    exit_code = pytest.main(args)
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
    from ..core.config import get_params
    from ..core.db import connect
    from ..core.workflow import lineage_str_to_steps
    from ..core.config import get_merged_params_for_lineage

    db = connect()
    orig_params = get_params(db)

    if lineage:
        lin_str = lineage
    else:
        # Build from integer IDs.
        # Lineage convention (A1): mwcs jobs encode both parents as
        # …/stack_N/refstack_M/mwcs_1.  A standalone refstack export
        # (no downstream dvv step) uses …/filter_N/refstack_M instead.
        _has_dvv = any([mwcsid, mwcsdttid, stretchingid, stretchingdvvid, wctid, wctdttid, waveletdvvid])
        parts = [f"preprocess_{preprocessid}", f"cc_{ccid}", f"filter_{filterid}"]
        if _has_dvv:
            # Full dvv lineage: stack and refstack both encoded
            parts.append(f"stack_{stackid}")
            if refstackid:
                parts.append(f"refstack_{refstackid}")
        else:
            # Stack-only or refstack-only export
            parts.append(f"stack_{stackid}")
            if refstackid:
                # Standalone refstack: hangs directly off filter, not stack
                # Replace stack with filter-rooted refstack path
                parts = parts[:-1] + [f"refstack_{refstackid}"]
        if mwcsid:
            parts.append(f"mwcs_{mwcsid}")
        if mwcsdttid:
            parts.append(f"mwcs_dtt_{mwcsdttid}")
        if stretchingid:
            parts.append(f"stretching_{stretchingid}")
        if stretchingdvvid:
            parts.append(f"stretching_dvv_{stretchingdvvid}")
        if wctid:
            parts.append(f"wavelet_{wctid}")
        if wctdttid:
            parts.append(f"wavelet_dtt_{wctdttid}")
        if waveletdvvid:
            parts.append(f"wavelet_dtt_dvv_{waveletdvvid}")
        lin_str = "/".join(parts)

    steps = lineage_str_to_steps(db, lin_str, sep="/", strict=False)
    _, lineage_names, params = get_merged_params_for_lineage(db, orig_params, {}, steps)

    out_path = output or f"params_{lin_str.replace('/', '_')}.yaml"
    params.to_yaml(out_path)
    click.echo(f"Exported params for lineage '{lin_str}' -> {out_path}")
    db.close()



@utils.command(name="export-dvv")
@click.option("-l", "--lineage", default=None, type=str,
              help="Slash-separated lineage string ending at a DVV step")
@click.option("-p",   "--preprocessid",    default=1,    type=int)
@click.option("-cc",  "--ccid",            default=1,    type=int)
@click.option("-f",   "--filterid",        default=1,    type=int)
@click.option("-s",   "--stackid",         default=1,    type=int)
@click.option("-r",   "--refstackid",      default=1,    type=int)
@click.option("-m",   "--mwcsid",          default=None, type=int)
@click.option("-md",  "--mwcsdttid",       default=None, type=int)
@click.option("-mdd", "--mwcsdvvid",       default=None, type=int)
@click.option("-S",   "--stretchingid",    default=None, type=int)
@click.option("-sd",  "--stretchingdvvid", default=None, type=int)
@click.option("-W",   "--wctid",           default=None, type=int)
@click.option("-wd",  "--wctdttid",        default=None, type=int)
@click.option("-wdv", "--waveletdvvid",    default=None, type=int)
@click.option("--pair-type",  default="CC", show_default=True)
@click.option("--components", default=None, type=str,
              help="Component to export (e.g. ZZ). Default: all.")
@click.option("--mov-stack",  default=None, type=str,
              help="Moving stack as window_step (e.g. 1D_1D). Default: all.")
@click.option("-o", "--output", default=".", type=str, show_default=True,
              help="Output directory (or full .nc path for a single export)")
@click.pass_context
def utils_export_dvv(ctx, lineage, preprocessid, ccid, filterid, stackid,
                     refstackid, mwcsid, mwcsdttid, mwcsdvvid,
                     stretchingid, stretchingdvvid, wctid, wctdttid,
                     waveletdvvid, pair_type, components, mov_stack, output):
    """Export dv/v time series as NetCDF with embedded parameter provenance.

    Each output file contains dv/v statistics (mean, std, median, n_pairs,
    weighted/trimmed variants) on a ``times`` dimension, with the complete
    processing parameter chain embedded as a global YAML attribute for full
    reproducibility.

    Either supply --lineage as a slash-separated string ending at a DVV step,
    or build it from individual step IDs.

    Examples::

        # MWCS dv/v, all components and mov_stacks, output to current directory
        msnoise utils export-dvv -p 1 -cc 1 -f 1 -s 1 -r 1 -m 1 -md 1 -mdd 1

        # Stretching dv/v, ZZ only, custom output directory
        msnoise utils export-dvv -p 1 -cc 1 -f 1 -s 1 -r 1 -S 1 -sd 1 \
            --components ZZ --output /data/release/

        # Using a lineage string
        msnoise utils export-dvv \
            --lineage preprocess_1/cc_1/filter_1/stack_1/refstack_1/mwcs_1/mwcs_dtt_1/mwcs_dtt_dvv_1

        # Reload the result in Python
        # >>> import xarray as xr, yaml
        # >>> ds = xr.open_dataset("dvv_CC_ZZ__pre1-...nc")
        # >>> params = yaml.safe_load(ds.attrs["msnoise_params"])
        # >>> print(params["mwcs"]["mwcs_wlen"])
    """
    from ..core.db import connect
    from ..core.workflow import lineage_str_to_steps
    from ..results import MSNoiseResult

    db = connect()

    if lineage:
        lin_str = lineage
    else:
        parts = [f"preprocess_{preprocessid}", f"cc_{ccid}",
                 f"filter_{filterid}", f"stack_{stackid}",
                 f"refstack_{refstackid}"]
        if mwcsid:
            parts.append(f"mwcs_{mwcsid}")
        if mwcsdttid:
            parts.append(f"mwcs_dtt_{mwcsdttid}")
        if mwcsdvvid:
            parts.append(f"mwcs_dtt_dvv_{mwcsdvvid}")
        if stretchingid:
            parts.append(f"stretching_{stretchingid}")
        if stretchingdvvid:
            parts.append(f"stretching_dvv_{stretchingdvvid}")
        if wctid:
            parts.append(f"wavelet_{wctid}")
        if wctdttid:
            parts.append(f"wavelet_dtt_{wctdttid}")
        if waveletdvvid:
            parts.append(f"wavelet_dtt_dvv_{waveletdvvid}")
        lin_str = "/".join(parts)

    steps      = lineage_str_to_steps(db, lin_str, sep="/", strict=False)
    step_names = [s.step_name for s in steps]
    result     = MSNoiseResult(db, step_names)

    ms_tuple = None
    if mov_stack:
        parts_ms = mov_stack.split("_", 1)
        if len(parts_ms) == 2:
            ms_tuple = tuple(parts_ms)
        else:
            raise click.BadParameter(
                f"--mov-stack must be window_step (e.g. 1D_1D), got {mov_stack!r}")

    written = result.export_dvv(
        output,
        pair_type=pair_type,
        components=components,
        mov_stack=ms_tuple,
    )
    for fpath in written:
        click.echo(f"Written: {fpath}")
    db.close()





@utils.command(name="create_preprocess_jobs")
@click.option("--date", default=None, type=str,
              help="Single date to create jobs for (YYYY-MM-DD).")
@click.option("--date_range", nargs=2, default=(None, None), type=str,
              metavar="START END",
              help="Date range: START END (inclusive, YYYY-MM-DD).")
@click.option("--set-number", default=1, type=int,
              help="Preprocess config set number (default 1).")
def create_preprocess_jobs_cmd(date, date_range, set_number):
    """Create preprocess T jobs for FDSN/EIDA sources (bypasses scan_archive).

    For local/SDS sources, only creates jobs where DataAvailability records
    already exist for the requested date(s).  For FDSN/EIDA sources, creates
    jobs for all active stations regardless (the preprocess worker fetches on
    demand).

    Examples::

        msnoise utils create_preprocess_jobs --date 2026-03-28

        msnoise utils create_preprocess_jobs --date_range 2026-01-01 2026-03-28
    """
    import datetime
    from ..core.db import connect
    from ..core.stations import (get_stations, resolve_data_source,
                                 get_data_availability)
    from ..core.fdsn import is_remote_source
    from ..core.workflow import get_workflow_steps
    from ..s02_new_jobs import update_job

    if not date and not any(date_range):
        raise click.UsageError("Provide --date or --date_range START END.")
    if date and any(date_range):
        raise click.UsageError("Provide --date OR --date_range, not both.")

    # Build date list
    if date:
        dates = [date]
    else:
        start = datetime.date.fromisoformat(date_range[0])
        end   = datetime.date.fromisoformat(date_range[1])
        dates = []
        cur = start
        while cur <= end:
            dates.append(cur.isoformat())
            cur += datetime.timedelta(days=1)

    db = connect()

    # Find the preprocess step
    steps = get_workflow_steps(db)
    pre_step = next(
        (s for s in steps
         if s.category == "preprocess" and s.set_number == set_number),
        None
    )
    if pre_step is None:
        click.echo(f"No preprocess step with set_number={set_number} found.", err=True)
        db.close()
        return

    stations = list(get_stations(db, all=False))
    if not stations:
        click.echo("No active stations found.", err=True)
        db.close()
        return

    created = 0
    skipped = 0

    for goal_day in dates:
        gd = datetime.datetime.strptime(goal_day, "%Y-%m-%d")
        gd_end = gd + datetime.timedelta(seconds=86401)

        for station in stations:
            ds  = resolve_data_source(db, station)

            if is_remote_source(ds.uri):
                # FDSN/EIDA — create job unconditionally
                locs = station.locs() or [""]
                for loc in locs:
                    pair = f"{station.net}.{station.sta}.{loc}"
                    update_job(db, goal_day, pair, pre_step.step_name, "T",
                               step_id=pre_step.step_id,
                               lineage=pre_step.step_name,
                               commit=False)
                    created += 1
            else:
                # Local/SDS — only create if DA records exist for this day
                for loc in (station.locs() or [""]):
                    da_records = get_data_availability(
                        db, net=station.net, sta=station.sta,
                        loc=loc, starttime=gd, endtime=gd_end
                    )
                    if da_records:
                        pair = f"{station.net}.{station.sta}.{loc}"
                        update_job(db, goal_day, pair, pre_step.step_name, "T",
                                   step_id=pre_step.step_id,
                                   lineage=pre_step.step_name,
                                   commit=False)
                        created += 1
                    else:
                        skipped += 1

        db.commit()

    click.echo(f"Created {created} preprocess T job(s) across {len(dates)} day(s) "
               f"({skipped} skipped — no DA records for local sources).")
    db.close()


@utils.command(name="create_psd_jobs")
@click.option("--date", default=None, type=str,
              help="Single date to create jobs for (YYYY-MM-DD).")
@click.option("--date_range", nargs=2, default=(None, None), type=str,
              metavar="START END",
              help="Date range: START END (inclusive, YYYY-MM-DD).")
@click.option("--set-number", default=1, type=int,
              help="PSD config set number (default 1).")
def create_psd_jobs_cmd(date, date_range, set_number):
    """Create psd T jobs for FDSN/EIDA sources (bypasses scan_archive).

    For FDSN/EIDA sources, creates jobs for all active stations for the
    requested date(s) regardless of DataAvailability (the psd worker fetches
    on demand).  For local/SDS sources, only creates jobs where
    DataAvailability records already exist for the requested date(s).

    Examples::

        msnoise utils create_psd_jobs --date 2026-03-28

        msnoise utils create_psd_jobs --date_range 2026-01-01 2026-03-28
    """
    import datetime
    from ..core.db import connect
    from ..core.stations import (get_stations, resolve_data_source,
                                 get_data_availability)
    from ..core.fdsn import is_remote_source
    from ..core.workflow import get_workflow_steps
    from ..s02_new_jobs import update_job

    if not date and not any(date_range):
        raise click.UsageError("Provide --date or --date_range START END.")
    if date and any(date_range):
        raise click.UsageError("Provide --date OR --date_range, not both.")

    # Build date list
    if date:
        dates = [date]
    else:
        start = datetime.date.fromisoformat(date_range[0])
        end   = datetime.date.fromisoformat(date_range[1])
        dates = []
        cur = start
        while cur <= end:
            dates.append(cur.isoformat())
            cur += datetime.timedelta(days=1)

    db = connect()

    # Find the psd step
    steps = get_workflow_steps(db)
    psd_step = next(
        (s for s in steps
         if s.category == "psd" and s.set_number == set_number),
        None
    )
    if psd_step is None:
        click.echo(f"No psd step with set_number={set_number} found.", err=True)
        db.close()
        return

    stations = list(get_stations(db, all=False))
    if not stations:
        click.echo("No active stations found.", err=True)
        db.close()
        return

    created = 0
    skipped = 0

    for goal_day in dates:
        gd     = datetime.datetime.strptime(goal_day, "%Y-%m-%d")
        gd_end = gd + datetime.timedelta(seconds=86401)

        for station in stations:
            ds     = resolve_data_source(db, station)
            locs   = station.locs() or [""]

            if is_remote_source(ds.uri):
                # FDSN/EIDA — create job unconditionally
                for loc in locs:
                    pair = f"{station.net}.{station.sta}.{loc}"
                    update_job(db, goal_day, pair, psd_step.step_name, "T",
                               step_id=psd_step.step_id,
                               lineage=psd_step.step_name,
                               commit=False)
                    created += 1
            else:
                # Local/SDS — only create if DA records exist for this day
                for loc in locs:
                    da_records = get_data_availability(
                        db, net=station.net, sta=station.sta,
                        loc=loc, starttime=gd, endtime=gd_end
                    )
                    if da_records:
                        pair = f"{station.net}.{station.sta}.{loc}"
                        update_job(db, goal_day, pair, psd_step.step_name, "T",
                                   step_id=psd_step.step_id,
                                   lineage=psd_step.step_name,
                                   commit=False)
                        created += 1
                    else:
                        skipped += 1

        db.commit()

    click.echo(f"Created {created} psd T job(s) across {len(dates)} day(s) "
               f"({skipped} skipped — no DA records for local sources).")
    db.close()




@utils.command(name="import-stationxml")
@click.argument("source")
@click.option("--data-source-id", "-d", default=None, type=int,
              help="DataSource ID to assign to imported stations (default: project default).")
@click.option("--no-save", is_flag=True, default=False,
              help="Do NOT save the inventory to response_path after import.")
@click.pass_context
def utils_import_stationxml(ctx, source, data_source_id, no_save):
    """Import stations from a StationXML file or FDSN URL.

    SOURCE can be a local file path or a URL (e.g. an FDSN station web
    service query with level=channel).

    The parsed inventory is written to the project's ``response_path``
    directory by default, making instrument responses immediately available
    to the preprocessing step.  Use --no-save to skip this.

    \b
    Examples:
      msnoise utils import-stationxml inventory.xml
      msnoise utils import-stationxml https://eida.ethz.ch/fdsnws/station/1/query?...
      msnoise utils import-stationxml inventory.xml --data-source-id 2 --no-save
    """
    from ..core.db import connect
    from ..core.stations import import_stationxml

    db = connect()
    try:
        created, updated, saved_path = import_stationxml(
            db, source,
            data_source_id=data_source_id,
            save_to_response_path=not no_save,
        )
        click.echo(f"✓ Import complete: {created} station(s) created, {updated} updated.")
        if saved_path:
            click.echo(f"  Inventory saved to: {saved_path}")
        elif not no_save:
            click.echo("  (Inventory not saved — response_path not configured.)")
    except Exception as exc:
        click.echo(f"✗ Import failed: {exc}", err=True)
        ctx.exit(1)
    finally:
        db.close()





@utils.command(name="run_workflow")
@click.option("-t", "--threads", default=1, show_default=True,
              help="Number of parallel workers per step.")
@click.option("--hpc", is_flag=True, default=None,
              help="Force HPC mode: insert 'new_jobs --after' between steps "
                   "(overrides db config).")
@click.option("--from", "from_category", default=None, metavar="CATEGORY",
              help="Start from this category (skip earlier steps).")
@click.option("--until", "until_category", default=None, metavar="CATEGORY",
              help="Stop after this category.")
@click.option("--on-failure", "on_failure",
              type=click.Choice(["stop", "continue"]), default="stop",
              show_default=True,
              help="What to do when a step exits non-zero.")
@click.option("--chunk-size", "chunk_size", default=0, show_default=True,
              help="--chunk-size passed to 'cc compute' and 'qc compute_psd'.")
@click.option("--dry-run", is_flag=True,
              help="Print the execution plan without running anything.")
@click.option("--export-script", "export_script", default=None,
              metavar="PATH",
              help="Write a shell script to PATH instead of executing.")
@click.pass_context
def utils_run_workflow(ctx, threads, hpc, from_category, until_category,
                       on_failure, chunk_size, dry_run, export_script):
    """Run all workflow steps in topological order.

    Discovers active steps from the database, sorts them by dependency order,
    and executes the appropriate ``msnoise`` sub-command for each category.
    Stops (or continues, see --on-failure) if any step fails.

    Examples:

    \b
      msnoise utils run_workflow                    # run everything, 1 worker
      msnoise utils run_workflow -t 8               # 8 parallel workers
      msnoise utils run_workflow --from stack        # resume from stack
      msnoise utils run_workflow --until mwcs_dtt    # stop after mwcs_dtt
      msnoise utils run_workflow --dry-run           # preview plan
      msnoise utils run_workflow --export-script run.sh  # write shell script
    """
    import subprocess
    import datetime as _dt
    from ..core.db import connect
    from ..core.workflow import run_workflow_plan, get_params

    session = connect()

    # Resolve HPC: CLI flag wins; otherwise read from config
    hpc_mode = hpc
    if hpc_mode is None:
        try:
            params = get_params(session)
            hpc_mode = bool(getattr(getattr(params, "global_", None), "hpc", False))
        except Exception:
            hpc_mode = False

    plan = run_workflow_plan(
        session,
        threads=threads,
        hpc=hpc_mode,
        from_category=from_category,
        until_category=until_category,
        chunk_size=chunk_size,
    )
    session.close()

    if not plan:
        click.echo("No runnable steps found in the active workflow.")
        return

    # ── Build the command list for each step ──────────────────────────────
    # Each entry: (display_label, [full argv list], hpc_after bool)
    import shutil
    _msnoise_bin = shutil.which("msnoise")
    if _msnoise_bin is None:
        import pathlib
        _msnoise_bin = str(pathlib.Path(sys.executable).parent / "msnoise")
    msnoise_exe = [_msnoise_bin]

    def _step_argv(step):
        """Full argv for a RunStep.  -t N and verbosity flags are top-level msnoise options."""
        base = msnoise_exe[:]
        if threads != 1:
            base += ["-t", str(threads)]
        verbosity = ctx.obj.get("MSNOISE_verbosity", "INFO")
        if verbosity == "DEBUG":
            base += ["-v"]
        elif verbosity == "WARNING":
            base += ["-q"]
        return base + step.cmd_tokens

    def _hpc_argv(category):
        return msnoise_exe + ["new_jobs", "--after", category]

    # ── Dry-run / export-script mode ─────────────────────────────────────
    def _fmt_cmd(argv):
        """Human-readable command: replace the full msnoise path with 'msnoise'."""
        tokens = argv[1:]   # strip the msnoise binary path
        return "msnoise " + " ".join(tokens)
        return "msnoise " + " ".join(tokens)

    if dry_run or export_script:
        lines = []
        if export_script:
            ts = _dt.datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
            cwd = os.getcwd()
            lines += [
                "#!/bin/bash",
                f"# MSNoise workflow script — generated {ts}",
                f"# Project: {cwd}",
                f"# Steps: {' '.join(s.category for s in plan)}",
                "set -e" if on_failure == "stop" else "# set -e  (--on-failure=continue)",
                "",
            ]

        for i, step in enumerate(plan, 1):
            tag = f"[{i:02d}/{len(plan):02d}] {step.category}"
            jobs_hint = f"  # {step.job_count} T jobs" if step.job_count else ""
            cmd_str = _fmt_cmd(_step_argv(step))
            if export_script:
                lines.append(f"echo '{tag}'")
                lines.append(cmd_str + jobs_hint)
                if step.hpc_after:
                    hpc_str = _fmt_cmd(_hpc_argv(step.category))
                    lines.append(f"echo '    [HPC] new_jobs --after {step.category}'")
                    lines.append(hpc_str)
                lines.append("")
            else:
                click.echo(f"{tag:<30} {cmd_str}{jobs_hint}")
                if step.hpc_after:
                    hpc_str = _fmt_cmd(_hpc_argv(step.category))
                    click.echo(f"{'':30} {hpc_str}  # HPC propagation")

        if export_script:
            script_text = "\n".join(lines) + "\n"
            with open(export_script, "w") as fh:
                fh.write(script_text)
            os.chmod(export_script, 0o755)
            click.echo(f"✓ Script written to: {export_script}")
        return

    # ── Execute ──────────────────────────────────────────────────────────
    results = []   # (category, status, elapsed, detail)
    t_total_start = time.time()
    failed_cats = []

    for i, step in enumerate(plan, 1):
        n = len(plan)
        jobs_hint = f" ({step.job_count} T jobs)" if step.job_count else " (0 T jobs — may be a no-op)"
        click.echo(f"\n[{i:02d}/{n:02d}] {step.category}{jobs_hint}")

        argv = _step_argv(step)
        click.echo(f"      → {_fmt_cmd(argv)}")

        t0 = time.time()
        result = subprocess.run(argv, cwd=os.getcwd())
        elapsed = time.time() - t0

        if result.returncode == 0:
            click.echo(f"      ✓ done in {elapsed:.1f}s")
            results.append((step.category, "ok", elapsed, None))
        else:
            detail = f"exit {result.returncode}"
            click.echo(f"      ✗ FAILED ({detail}) in {elapsed:.1f}s", err=True)
            results.append((step.category, "fail", elapsed, detail))
            failed_cats.append(step.category)
            if on_failure == "stop":
                click.echo("      Stopping (--on-failure=stop).", err=True)
                break

        # HPC propagation step
        if step.hpc_after and (result.returncode == 0 or on_failure == "continue"):
            hpc_argv = _hpc_argv(step.category)
            click.echo(f"      → {_fmt_cmd(hpc_argv)}  [HPC propagation]")
            hpc_result = subprocess.run(hpc_argv, cwd=os.getcwd())
            if hpc_result.returncode != 0:
                click.echo(f"      ✗ new_jobs --after {step.category} failed", err=True)
                failed_cats.append(f"new_jobs --after {step.category}")

    # ── Summary ──────────────────────────────────────────────────────────
    total_elapsed = time.time() - t_total_start
    n_ok   = sum(1 for _, s, _, _ in results if s == "ok")
    n_fail = sum(1 for _, s, _, _ in results if s == "fail")
    n_skip = len(plan) - len(results)

    click.echo(f"\n{'═' * 54}")
    click.echo(f" Run complete in {total_elapsed:.1f}s: "
               f"{n_ok} ✓  {n_fail} ✗  {n_skip} skipped")
    if failed_cats:
        click.echo(f" Failed: {', '.join(failed_cats)}", err=True)
        ctx.exit(1)


@utils.command(name="download")
@click.option("--sds-path", "sds_path", default=None,
              type=click.Path(file_okay=False, writable=True),
              metavar="PATH",
              help="SDS archive root to write into.  Auto-detected from the "
                   "DataSource table when unambiguous; falls back to ./SDS.")
@click.option("--startdate", default=None, metavar="YYYY-MM-DD",
              help="Override project startdate.")
@click.option("--enddate", default=None, metavar="YYYY-MM-DD",
              help="Override project enddate.")
@click.pass_context
def utils_download(ctx, sds_path, startdate, enddate):
    """Download waveforms from FDSN/EIDA sources into an SDS archive.

    Reads remote DataSource entries and the station table from the database,
    builds an ObsPy Inventory for server-side station gating, and runs
    ObsPy MassDownloader for the project date range.  Waveforms are written
    as day-aligned MiniSEED files in SDS layout.

    The SDS write root is resolved in order:
    --sds-path → single unambiguous local SDS DataSource → ./SDS (warning).

    StationXML is stored alongside and never overwritten.  Traces are never
    discarded for missing instrument response.

    Example::

        msnoise utils download
        msnoise utils download --sds-path /data/SDS
        msnoise utils download --startdate 2013-06-01 --enddate 2013-06-30
    """
    from ..core.db import connect
    from ..msnoise_table_def import declare_tables
    from ..core.fdsn import mass_download
    from pathlib import Path

    db      = connect()
    _schema = declare_tables()
    mass_download(db, _schema,
                  sds_root=Path(sds_path) if sds_path else None,
                  startdate_override=startdate,
                  enddate_override=enddate)


# ---------------------------------------------------------------------------
# msnoise project  — project archive commands (no DB required)
# ---------------------------------------------------------------------------

@cli.group(cls=OrderedGroup)
def project():
    """Commands for exporting and importing MSNoise project archives."""
    pass


@project.command(name="export")
@click.option(
    "--level",
    "levels",
    multiple=True,
    metavar="LEVEL",
    help="Level(s) to export.  Repeat for multiple (e.g. --level stack --level dvv).  "
         "Pass 'all' to export every level that has content on disk.",
)
@click.option(
    "--output", "-o",
    default=None,
    metavar="PATH",
    help="Destination .tar.zst file (single level only).",
)
@click.option(
    "--output-dir",
    default=None,
    metavar="DIR",
    help="Destination directory for multi-level export.  One .tar.zst per level "
         "plus bundle_pointer.yaml are written here.",
)
@click.option(
    "--url-base",
    default="",
    metavar="URL",
    help="URL prefix for bundle_pointer.yaml (e.g. https://ftp.example.org/msnoise/study).  "
         "Placeholder strings are written when omitted.",
)
@click.option(
    "--project-dir",
    default=".",
    show_default=True,
    metavar="DIR",
    help="MSNoise project root (default: current directory).",
)
def project_export(levels, output, output_dir, url_base, project_dir):
    """Export project archive(s) at one or more entry levels.

    Single level — writes one .tar.zst and prints its SHA-256::

        msnoise project export --level stack -o /data/level_stack.tar.zst

    Multiple levels — writes one archive per level plus bundle_pointer.yaml::

        msnoise project export --level stack --level dvv \\
            --output-dir /data/bundles/

    All levels with content on disk::

        msnoise project export --level all --output-dir /data/bundles/ \\
            --url-base https://ftp.seismology.be/msnoise/study
    """
    from ..core.project_io import export_project, export_project_levels

    if not levels:
        raise click.UsageError("Provide at least one --level (or --level all).")

    # Single-level shorthand: --level X -o file.tar.zst
    if list(levels) != ["all"] and len(levels) == 1 and output and not output_dir:
        level = levels[0]
        click.echo(f"Exporting level={level!r} from {project_dir!r} → {output!r}")
        sha = export_project(project_dir, level, output)
        click.echo(f"Done.  SHA-256: {sha}")
        click.echo(f"Paste into bundle_pointer.yaml:\n  sha256: \"{sha}\"")
        return

    # Multi-level path — requires --output-dir
    if not output_dir:
        raise click.UsageError(
            "Multi-level export requires --output-dir.  "
            "Use -o FILE only for a single --level."
        )

    level_arg = "all" if list(levels) == ["all"] else list(levels)
    click.echo(f"Exporting level(s)={level_arg!r} → {output_dir!r}")
    exported = export_project_levels(project_dir, level_arg, output_dir, url_base=url_base)
    if not exported:
        click.echo("No levels had content to export.")
    else:
        click.echo(f"\nExported {len(exported)} archive(s): {list(exported)}")
        click.echo("Review and update the URLs in bundle_pointer.yaml before publishing.")


@project.command(name="import")
@click.option(
    "--from", "pointer",
    required=True,
    metavar="PATH",
    help="Path to bundle_pointer.yaml.",
)
@click.option(
    "--level",
    "levels",
    multiple=True,
    metavar="LEVEL",
    help="Entry level(s) to download and import.  Repeat for multiple levels "
         "(e.g. --level stack --level dvv).  Pass 'all' to download every "
         "level listed in bundle_pointer.yaml.",
)
@click.option(
    "--project-dir",
    default=".",
    show_default=True,
    metavar="DIR",
    help="Destination directory (created if absent).",
)
@click.option(
    "--with-jobs",
    is_flag=True,
    default=False,
    help="Reconstruct flag=D jobs from the extracted _output/ tree "
         "so the pipeline can be resumed immediately.",
)
def project_import(pointer, levels, project_dir, with_jobs):
    """Download and import project archive(s) from a bundle_pointer.yaml.

    Downloads .tar.zst archive(s) for the requested level(s), verifies
    SHA-256, extracts all into the same project directory, initialises
    the database, and optionally reconstructs flag=D jobs.

    Examples::

        # single level
        msnoise project import --from bundle_pointer.yaml --level stack

        # multiple levels — all extracted into the same project dir
        msnoise project import --from bundle_pointer.yaml \\
            --level stack --level dvv --project-dir ./my_project --with-jobs

        # everything in bundle_pointer.yaml
        msnoise project import --from bundle_pointer.yaml \\
            --level all --project-dir ./my_project --with-jobs
    """
    from ..core.project_io import import_project_archive
    from ..project import MSNoiseProject

    if not levels:
        raise click.UsageError("Provide at least one --level (or --level all).")

    # "all" as a single value means download everything
    level_arg = "all" if list(levels) == ["all"] else list(levels)

    root = import_project_archive(pointer, level_arg, project_dir)
    click.echo(f"Extracted to {root}")

    proj = MSNoiseProject.from_project_dir(root)
    # Record imported levels for init_db with_jobs reconstruction
    proj._imported_levels = (
        None if level_arg == "all" else level_arg  # None → reads meta.yaml per archive
    )
    proj.init_db(with_jobs=with_jobs)
    click.echo("Database initialised.")

    if with_jobs:
        effective = level_arg if isinstance(level_arg, list) else list(levels)
        click.echo(
            f"\nReady.  Resume the pipeline with:\n"
            + "\n".join(f"  msnoise new_jobs --after {lv}" for lv in effective)
        )


def run():
    try:
        cli(obj={})
    except MSNoiseError as e:
        logging.critical(str(e))
