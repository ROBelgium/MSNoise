import os
import sys
import click
import pkg_resources
import logging
import time
import traceback
from pkg_resources import iter_entry_points
# from click_plugins import with_plugins


@click.group()
@click.option('-t', '--threads', default=1, help='Number of threads to use \
(only affects modules that are designed to do parallel processing)')
@click.option('-c', '--custom', default=False, is_flag=True, help='Use custom \
 file for plots. To use this, copy the plot script here and edit it.')
@click.option('-v', '--verbose', default=2, count=True)
@click.pass_context
def cli(ctx, threads, custom, verbose):
    ctx.obj['MSNOISE_threads'] = threads
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


@click.command()
def test():
    """Runs the test suite, should be executed in an empty folder!"""
    from ..test.tests import main
    main()


@click.command()
@click.option('-p', '--port', default=5000, help='Port to open')
def admin(port):
    from ..msnoise_admin import main
    main(port)


@click.command()
def upgrade_db():
    """Upgrade the database from pre-1.3 to MSNoise 1.3. This should only
    be ran once."""
    from sqlalchemy.exc import IntegrityError, OperationalError, InvalidRequestError
    from ..api import connect, Config, get_tech, get_engine
    from ..default import default
    db = connect()
    for name in default.keys():
        try:
            db.add(Config(name=name, value=default[name][-1]))
            db.commit()
        except:
            db.rollback()
            print("Passing %s: already in DB" % name)
            continue

    db.close()

    if get_tech() == 2:
        try:
            e = get_engine()
            e.execute('ALTER TABLE `jobs` CHANGE `type` `jobtype` VARCHAR( 10 )')
        except OperationalError:
            print( "The jobs table seems already up-to-date, exiting.")
    else:
        print("OK, the new config parameters have been inserted, but you need" \
              "to edit the `jobs` table manually in oder to match the new" \
              "columns naming.")
        print ("Please read http://msnoise.org/doc/releasenotes/msnoise-1.3.html")


@click.command()
@click.option('-j', '--jobs', is_flag=True, help='Jobs Info only')
def info(jobs):
    """Outputs general information about the current install and config, plus
    information about jobs and their status."""
    from ..api import connect, get_config, get_job_types
    from ..default import default

    click.echo('')
    click.echo('General:')

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
            if get_config(db, 'keep_all') in ['Y','y']:
                for job in get_job_types(db):
                    if job[1] == 'D':
                        if job[0] > 0:
                            click.secho(" - %s does not exists and that is not normal (%i CC jobs done)" % (output_folder, job[0]), fg='red')
                        else:
                            click.secho(" - %s does not exists and that is normal (%i CC jobs done)" % (output_folder, job[0]))
            else:
                click.secho(" - %s does not exists (and that is normal because keep_all=False)" % output_folder)


        click.echo('')
        click.echo('Raw config bits: "D"efault or "M"odified (green)')
        for key in default.keys():
            tmp = get_config(db, key)
            if tmp == default[key][1]:
                click.secho(" D %s: %s" %(key, tmp ))
            else:
                click.secho(" M %s: %s" %(key, tmp ), fg='green')

    click.echo('CC Jobs:')
    for (n,jobtype) in get_job_types(db,'CC'):
        click.echo(" %s : %i" % (jobtype, n))

    click.echo('')
    click.echo('DTT Jobs:')
    for (n,jobtype) in get_job_types(db,'DTT'):
        click.echo(" %s : %i" % (jobtype, n))


@click.command()
def install():
    """This command launches the installer."""
    click.echo('Launching the installer')
    from ..s000installer import main
    main()


@click.command()
@click.option('-s', '--set', help='Modify config value: usage --set name=value')
def config(set):
    """This command launches the Configurator or... TODO"""
    if set:
        from ..default import default
        if not set.count("="):
            click.echo("!! format of the set command is name=value !!")
            return
        name, value = set.split("=")
        if not name in default:
            click.echo("!! unknown parameter %s !!"%name)
            return
        from ..api import connect, update_config
        db = connect()
        update_config(db, name, value)
        db.commit()
        db.close()
        click.echo("Successfully updated parameter %s = %s"%(name,value))
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
    #click.echo('Working on %i threads' % ctx.obj['MSNOISE_threads'])
    from ..bugreport import main
    main(sys, modules, env, all)


@click.command()
def populate():
    """Rapidly scan the archive filenames and find Network/Stations"""
    from ..s002populate_station_table import main
    main()


@click.command()
@click.option('-i', '--init', is_flag=True, help='First run ?')
@click.pass_context
def scan_archive(ctx, init):
    """Scan the archive and insert into the Data Availability table."""
    from ..s01scan_archive import main
    main(init, threads=ctx.obj['MSNOISE_threads'])


@click.command()
@click.option('-i', '--init', is_flag=True, help='First run ?')
def new_jobs(init):
    """Determines if new CC jobs are to be defined"""
    from ..s02new_jobs import main
    main(init)


@click.command()
@click.pass_context
def compute_cc(ctx):
    """Computes the CC jobs (based on the "New Jobs" identified)"""
    from ..s03compute_cc import main
    from multiprocessing import Process
    threads = ctx.obj['MSNOISE_threads']
    processes = []
    for i in range(threads):
        p = Process(target=main)
        p.start()
        processes.append(p)
        time.sleep(1)
    for p in processes:
        p.join()


@click.command()
@click.option('-r', '--ref', is_flag=True, help='Compute the REF Stack')
@click.option('-m', '--mov', is_flag=True, help='Compute the MOV Stacks')
@click.option('-s', '--step', is_flag=True, help='Compute the STEP Stacks')
@click.option('-i', '--interval', default=1, help='Number of days before now to\
 search for modified Jobs')
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
    processes = []
    for i in range(threads):
        p = Process(target=main)
        p.start()
        processes.append(p)

    for p in processes:
        p.join()


@click.command()
def compute_stretching():
    """[experimental] Computes the stretching based on the new stacked data"""
    from ..stretch import main
    main()


@click.command()
@click.option('-i', '--interval', default=1, help='Number of days before now to\
 search for modified Jobs')
def compute_dtt(interval):
    """Computes the dt/t jobs based on the new MWCS data"""
    from ..s06compute_dtt import main
    main(interval)


@click.command()
@click.argument('jobtype')
@click.option('-a', '--all', is_flag=True, help='Reset all jobs')
def reset(jobtype, all):
    """Resets the job to "T"odo. ARG is [CC] or [DTT]. By default
    only resets jobs "I"n progress. --all resets all jobs, whatever
    the flag value"""
    from ..api import connect, reset_jobs
    session = connect()
    reset_jobs(session, jobtype, all)
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
@click.option('-m', '--mov_stack', default=0,  help='Plot specific mov stacks')
@click.option('-p', '--pair', default=None,  help='Plot a specific pair',
              multiple=True)
@click.option('-A', '--all', help='Show the ALL line?', is_flag=True)
@click.option('-M', '--dttname', default="M",  help='Plot M or M0?')
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
@click.option('-m', '--mov_stack', default=0,  help='Plot specific mov stacks')
@click.option('-p', '--pair', default=None,  help='Plot a specific pair',
              multiple=True)
@click.option('-A', '--all', help='Show the ALL line?', is_flag=True)
@click.option('-M', '--dttname', default="A",  help='Plot M or M0?')
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
@click.option('-m', '--mov_stack', default=1, help='Mov Stack to read from disk')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def interferogram(ctx, sta1, sta2, filterid, comp, mov_stack, show, outfile):
    """Plots the interferogram between sta1 and sta2 (parses the CCFs)\n
    STA1 and STA2 must be provided with this format: NET.STA !"""
    if ctx.obj['MSNOISE_custom']:
        from interferogram import main
    else:
        from ..plots.interferogram import main
    main(sta1, sta2, filterid, comp, mov_stack, show, outfile)

@click.command()
@click.argument('sta1')
@click.argument('sta2')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1, help='Mov Stack to read from disk')
@click.option('-a', '--ampli', default=5.0, help='Amplification')
@click.option('-S', '--seismic', is_flag=True, help='Seismic style')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def ccftime(ctx, sta1, sta2, filterid, comp, mov_stack,
            ampli, seismic, show, outfile):
    """Plots the ccf vs time between sta1 and sta2 (parses the dt/t results)\n
    STA1 and STA2 must be provided with this format: NET.STA !"""
    if ctx.obj['MSNOISE_custom']:
        from ccftime import main
    else:
        from ..plots.ccftime import main
    main(sta1, sta2, filterid, comp, mov_stack, ampli, seismic, show, outfile)


@click.command()
@click.argument('sta1')
@click.argument('sta2')
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1, help='Mov Stack to read from disk')
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
@click.pass_context
def distance(ctx, filterid, comp,  ampli, show, outfile):
    """Plots the REFs of all pairs vs distance"""
    if ctx.obj['MSNOISE_custom']:
        from distance import main
    else:
        from ..plots.distance import main
    main(filterid, comp, ampli, show, outfile)


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
@click.option('-f', '--filterid', default=1, help='Filter ID')
@click.option('-d', '--day', default="2015-01-01", help='Day')
@click.option('-c', '--comp', default="ZZ", help='Components (ZZ, ZR,...)')
@click.option('-m', '--mov_stack', default=1, help='Mov Stack to read from disk')
@click.option('-s', '--show', help='Show interactively?',
              default=True, type=bool)
@click.option('-o', '--outfile', help='Output filename (?=auto)',
              default=None, type=str)
@click.pass_context
def dtt(ctx, sta1, sta2, filterid, day, comp, mov_stack, show, outfile):
    """Plots a graph of dt against t\n
    STA1 and STA2 must be provided with this format: NET.STA !"""
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

cli.add_command(plugin)


def run():
    cli(obj={})
