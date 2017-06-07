import traceback
import logging
import os
import sys
import time

import click
import pkg_resources


# from click_plugins import with_plugins


@click.group()
@click.option('-t', '--threads', default=1, help='Number of threads to use \
(only affects modules that are designed to do parallel processing)')
@click.option('-d', '--delay', default=1,  help='In the case of multi-threading'
                    ', defines the number of seconds to wait before lauching '
                    'the next thread. Defaults to [1] second ')
@click.option('-c', '--custom', default=False, is_flag=True, help='Use custom \
 file for plots. To use this, copy the plot script here and edit it.')
@click.option('-v', '--verbose', default=2, count=True)
@click.pass_context
def cli(ctx, threads, delay, custom, verbose):
    ctx.obj['MSNOISE_threads'] = threads
    ctx.obj['MSNOISE_threadsdelay'] = delay
    if verbose == 0:
        ctx.obj['MSNOISE_verbosity'] = "WARNING"
    elif verbose == 1:
        ctx.obj['MSNOISE_verbosity'] = "INFO"
    elif verbose > 1:
        ctx.obj['MSNOISE_verbosity'] = "DEBUG"
    logging.basicConfig(level=ctx.obj['MSNOISE_verbosity'],
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    sys.path.append(os.getcwd())
    ctx.obj['MSNOISE_custom'] = custom

    pass


# @with_plugins(iter_entry_points('msnoise.plugins'))
@click.group()
def plugin():
    """Runs a command in a named plugin"""
    pass


# @with_plugins(iter_entry_points('msnoise.plugins'))
@click.group()
def p():
    """Short cut for plugins"""
    pass


@click.command()
def test():
    """Runs the test suite, should be executed in an empty folder!"""
    import matplotlib.pyplot as plt
    plt.switch_backend("agg")
    from ..test.tests import main
    main()


@click.command()
@click.option('-p', '--port', default=5000, help='Port to open')
def admin(port):
    from ..msnoise_admin import main
    main(port)


@click.command()
def upgrade_db():
    """Upgrade the database from previous to a new version.\n
    This procedure adds new parameters with their default value
    in the config database.
    """
    from ..api import connect, Config
    from ..default import default
    db = connect()
    for name in default.keys():
        try:
            db.add(Config(name=name, value=default[name][-1]))
            db.commit()
        except:
            db.rollback()
            # print("Passing %s: already in DB" % name)
            continue
    try:
        db.execute("CREATE INDEX job_index ON jobs (day, pair, jobtype)")
        db.commit()
    except:
        logging.info("It looks like the v1.5 'job_index' is already in the DB")
        db.rollback()

    try:
        db.execute("CREATE INDEX da_index ON data_availability (path, file, net, sta, comp)")
        db.commit()
    except:
        logging.info("It looks like the v1.5 'da_index' is already in the DB")
        db.rollback()

    db.close()


@click.command()
@click.option('-j', '--jobs', is_flag=True, help='Jobs Info only')
def info(jobs):
    """Outputs general information about the current install and config, plus
    information about jobs and their status."""
    from ..api import connect, get_config, get_job_types, get_filters, \
        get_stations
    from ..default import default

    click.echo('')
    click.echo('General:')

    def d(path):
        return os.path.split(path)[0]

    click.echo('MSNoise is installed in: %s'
               % d(d(d(os.path.abspath(__file__)))))

    if os.path.isfile('db.ini'):
        click.echo(' - db.ini is present')
    else:
        click.secho(' - db.ini is not present, is MSNoise installed here ?',
                    fg='red')
        return
    click.echo('')
    db = connect()

    if not jobs:
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

        click.echo('')
        click.echo('Raw config bits: "D"efault or "M"odified (green)')
        for key in default.keys():
            tmp = get_config(db, key)
            if tmp == default[key][1]:
                click.secho(" D %s: %s" % (key, tmp))
            else:
                click.secho(" M %s: %s" % (key, tmp), fg='green')

        click.echo('')
        click.echo('Filters:')
        print('ID: [low:high]  [mwcs_low:mwcs_high]    mwcs_wlen    mwcs_step'
              '   used')
        for f in get_filters(db, all=True):
            data = (f.ref,
                    f.low,
                    f.high,
                    f.mwcs_low,
                    f.mwcs_high,
                    f.mwcs_wlen,
                    f.mwcs_step,
                    ['N', 'Y'][f.used])
            print('%02i: [%.03f:%.03f] [%.03f:%.03f] %.03i %.03i %s' % data)

        click.echo('')
        click.echo('Stations:')
        for s in get_stations(db, all=True):
            data = (s.net, s.sta, s.X, s.Y, s.altitude, s.coordinates,
                    ['N', 'Y'][s.used])
            print('%s.%s %.4f %.4f %.1f %s %s' % data)

    click.echo('')
    click.echo('CC Jobs:')
    for (n, jobtype) in get_job_types(db, 'CC'):
        click.echo(" %s : %i" % (jobtype, n))

    click.echo('')
    click.echo('DTT Jobs:')
    for (n, jobtype) in get_job_types(db, 'DTT'):
        click.echo(" %s : %i" % (jobtype, n))


@click.command()
def install():
    """This command launches the installer."""
    click.echo('Launching the installer')
    from ..s000installer import main
    main()


@click.command()
@click.option('-s', '--set', help='Modify config value: usage --set name=value')
@click.option('-S', '--sync', is_flag=True, help='Sync station metadata from'
                                                 ' inventory/dataless')
def config(set, sync):
    """This command should now only be used to use the command line to set
    a parameter value in the data base. It used to launch the Configurator but
    the recommended way to configure MSNoise is to use the "msnoise admin" web
    interface."""
    if set:
        from ..default import default
        if not set.count("="):
            click.echo("!! format of the set command is name=value !!")
            return
        name, value = set.split("=")
        if name not in default:
            click.echo("!! unknown parameter %s !!" % name)
            return
        from ..api import connect, update_config
        db = connect()
        update_config(db, name, value)
        db.commit()
        db.close()
        click.echo("Successfully updated parameter %s = %s" % (name, value))
    elif sync:
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

    else:
        from ..s001configurator import main
        click.echo('Let\'s Configure MSNoise !')
        main()


@click.command()
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


@click.command()
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
        from ..msnoise_table_def import DataAvailability
        from ..api import update_station
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


@click.command()
@click.option('-i', '--init', is_flag=True, help='First run ?')
@click.option('--path',  help='Scan all files in specific folder, overrides the'
                              ' default workflow step.')
@click.pass_context
def scan_archive(ctx, init, path):
    """Scan the archive and insert into the Data Availability table."""
    if path:
        logging.info("Overriding workflow...")
        from obspy import UTCDateTime
        from ..s01scan_archive import worker
        db = connect()
        startdate = UTCDateTime(get_config(db, "startdate"))
        enddate = UTCDateTime(get_config(db, "enddate"))
        cc_sampling_rate = float(get_config(db, "cc_sampling_rate"))
        db.close()

        worker(sorted(os.listdir(path)), path, startdate, enddate,
               cc_sampling_rate, init=True)
    else:
        from ..s01scan_archive import main
        main(init, threads=ctx.obj['MSNOISE_threads'])


@click.command()
@click.option('-i', '--init', is_flag=True, help='First run ?')
@click.option('--nocc', is_flag=True, default=False, help='Disable the creation'
                                                          ' of CC jobs')
def new_jobs(init, nocc):
    """Determines if new CC jobs are to be defined"""
    from ..s02new_jobs import main
    main(init, nocc)


@click.command()
@click.pass_context
def compute_cc(ctx):
    """Computes the CC jobs (based on the "New Jobs" identified)"""
    from ..s03compute_cc import main
    from multiprocessing import Process
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']
    processes = []
    for i in range(threads):
        p = Process(target=main)
        p.start()
        processes.append(p)
        time.sleep(delay)
    for p in processes:
        p.join()


@click.command()
@click.option('-r', '--ref', is_flag=True, help='Compute the REF Stack')
@click.option('-m', '--mov', is_flag=True, help='Compute the MOV Stacks')
@click.option('-s', '--step', is_flag=True, help='Compute the STEP Stacks')
@click.option('-i', '--interval', default=1.0, help='Number of days before now to'
                                                  ' search for modified Jobs')
def stack(ref, mov, step, interval):
    """Stacks the [REF] and/or [MOV] windows"""
    click.secho('Lets STACK !', fg='green')
    from ..s04stack import main
    if ref:
        main('ref', interval)
    if mov:
        main('mov', interval)
    if step:
        main('step', interval)


@click.command()
@click.pass_context
def compute_mwcs(ctx):
    """Computes the MWCS based on the new stacked data"""
    from ..s05compute_mwcs import main
    from multiprocessing import Process
    threads = ctx.obj['MSNOISE_threads']
    delay = ctx.obj['MSNOISE_threadsdelay']

    processes = []
    for i in range(threads):
        p = Process(target=main)
        p.start()
        processes.append(p)
        time.sleep(delay)
    for p in processes:
        p.join()


@click.command()
def compute_stretching():
    """[experimental] Computes the stretching based on the new stacked data"""
    from ..stretch import main
    main()


@click.command()
@click.option('-i', '--interval', default=1.0, help='Number of days before now to\
 search for modified Jobs')
def compute_dtt(interval):
    """Computes the dt/t jobs based on the new MWCS data"""
    from ..s06compute_dtt import main
    main(interval)


@click.command()
@click.argument('jobtype')
@click.option('-a', '--all', is_flag=True, help='Reset all jobs')
@click.option('-r', '--rule', help='Reset job that match this SQL rule')
def reset(jobtype, all, rule):
    """Resets the job to "T"odo. ARG is [CC] or [DTT]. By default
    only resets jobs "I"n progress. --all resets all jobs, whatever
    the flag value"""
    from ..api import connect, reset_jobs
    session = connect()
    if jobtype == "DA":
        session.execute("update data_availability set flag='M' where 1")
    elif jobtype != jobtype.upper():
        logging.info("The jobtype %s is not uppercase (usually jobtypes"
                     " are uppercase...)"%jobtype)
    reset_jobs(session, jobtype, all, rule)
    session.close()


@click.command()
def ipython():
    """Launches an ipython notebook in the current folder"""
    os.system("ipython notebook --pylab inline --ip 0.0.0.0")


#
# PLOT GROUP
#

@click.group()
def plot():
    """Top level command to trigger different plots"""
    pass


@click.command()
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


@click.command()
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


@click.command()
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


@click.command()
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
@click.pass_context
def interferogram(ctx, sta1, sta2, filterid, comp, mov_stack, show, outfile,
                  refilter):
    """Plots the interferogram between sta1 and sta2 (parses the CCFs)\n
    STA1 and STA2 must be provided with this format: NET.STA !"""
    if ctx.obj['MSNOISE_custom']:
        from interferogram import main
    else:
        from ..plots.interferogram import main
    main(sta1, sta2, filterid, comp, mov_stack, show, outfile, refilter)


@click.command()
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
@click.pass_context
def ccftime(ctx, sta1, sta2, filterid, comp, mov_stack,
            ampli, seismic, show, outfile, envelope, refilter):
    """Plots the ccf vs time between sta1 and sta2 (parses the dt/t results)\n
    STA1 and STA2 must be provided with this format: NET.STA !"""
    if ctx.obj['MSNOISE_custom']:
        from ccftime import main
    else:
        from ..plots.ccftime import main
    main(sta1, sta2, filterid, comp, mov_stack, ampli, seismic, show, outfile,
         envelope, refilter)


@click.command()
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
    if ctx.obj['MSNOISE_custom']:
        from mwcs import main
    else:
        from ..plots.mwcs import main
    main(sta1, sta2, filterid, comp, mov_stack, show, outfile)


@click.command()
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
@click.pass_context
def distance(ctx, filterid, comp, ampli, show, outfile, refilter,
             virtual_source):
    """Plots the REFs of all pairs vs distance"""
    if ctx.obj['MSNOISE_custom']:
        from distance import main
    else:
        from ..plots.distance import main
    main(filterid, comp, ampli, show, outfile, refilter, virtual_source)


@click.command()
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


@click.command()
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
    if ctx.obj['MSNOISE_custom']:
        from dtt import main
    else:
        from ..plots.dtt import main
    main(sta1, sta2, filterid, comp, day, mov_stack, show, outfile)


# Add plot commands to the plot group:
plot.add_command(data_availability)
plot.add_command(dvv)
plot.add_command(interferogram)
plot.add_command(ccftime)
plot.add_command(mwcs)
plot.add_command(distance)
plot.add_command(station_map)
plot.add_command(timing)
plot.add_command(dtt)

# Add all commands to the cli group:
cli.add_command(info)
cli.add_command(admin)
cli.add_command(upgrade_db)
cli.add_command(install)
cli.add_command(config)
cli.add_command(populate)
cli.add_command(bugreport)
cli.add_command(scan_archive)
cli.add_command(new_jobs)
cli.add_command(compute_cc)
cli.add_command(stack)
cli.add_command(compute_mwcs)
cli.add_command(compute_stretching)
cli.add_command(compute_dtt)
cli.add_command(reset)
cli.add_command(ipython)
cli.add_command(test)
# Finally add the plot group too:
cli.add_command(plot)

try:
    from ..api import connect, get_config

    db = connect()
    plugins = get_config(db, "plugins")
except:
    plugins = None

if plugins:
    plugins = plugins.split(",")
    for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.commands'):
        module_name = ep.module_name.split(".")[0]
        if module_name in plugins:
            plugin.add_command(ep.load())
            p.add_command(ep.load())

cli.add_command(plugin)
cli.add_command(p)


def run():
    cli(obj={})
