"""
One advantage of MSNoise is its ability to be used as an automated monitoring
tool. In order to run every night on the data acquired during the previous day,
MSNoise needs to check the data archive for new or modified files.

Those files could have been acquired during the last day, but be data of a
previously offline station and contain useful information for, say, a month ago.
The time to search for is defined in the config from the 'crondays' value. For
convenience, this parameter can be temporarily redefined on the command line
using the `--crondays` option of the `scan_archive` subcommand. In both cases,
it can be a float designating a number of days in the past, or a string
designating a number of weeks, days, and/or hours in the format 'Xw Xd Xh'
(each group being optional, as well as the separating blanks).

The scan_archive script inspects the modified time attribute ('mtime') of files
in the archives to locate new or modified files. Once located, they are inserted
(if new) or updated (if modified) in the data availability table.

To run the code on two Process, execute the following in console:

.. code-block:: sh

    $ msnoise -t 2 scan_archive

Special case: first run
~~~~~~~~~~~~~~~~~~~~~~~~

This script is the same as for the routine, but one has to pass the --init
option:

.. code-block:: sh

    $ msnoise -t 2 scan_archive --init

This will scan the data_archive folder the configured stations and will insert
all files found in the data_availability table in the database. As usual,
calling the script with a --help argument will show its usage.

.. _scan-archive-expert:

Expert (lazy) mode:
~~~~~~~~~~~~~~~~~~~

Sometimes, you only want to scan a few files and run MSNoise on them. To do this
simply run:

.. code-block:: sh

    $ msnoise scan_archive --path /path/to/where/files/are

and MSNoise will read anything ObsPy can (provided the files have a proper
header (network code, station code and channel code). Then, once done, simply
run the :ref:`"populate from DataAvailability"<populate-expert>` procedure.

This command can also scan folders recursively:

.. code-block:: sh

    $ msnoise scan_archive --path /path/to/archive --recursively
"""


import argparse
import datetime
import glob
import logging
import multiprocessing
import obspy
import os
import re
import sys
import time

import traceback

# Use the built-in version of scandir/walk if possible,
# otherwise use the scandir module version
try:
    from os import scandir, walk
except ImportError:
    from scandir import scandir, walk

from . import api
from . import FatalError
from . import data_structures


import logbook
logger = logbook.Logger(__name__)


def update_availability(db, folder, basename, data):
    """
    Updates the availability of data in the database from a miniseed data
    stream. Returns the return code for the update.

    :param db: the sqlalchemy session object
    :param folder: the directory where lies the miniseed file
    :param basename: the basename of the miniseed file to read data from
    :param data: obspy.core.stream.Stream object for the channel
    """
    gaps = data.get_gaps()
    gaps_duration = 0
    for gap in gaps:
        gaps_duration += gap[6]

    data_duration = 0
    start = datetime.datetime(year=2100, month=1, day=1)
    stop = datetime.datetime(year=1900, month=1, day=1)
    for trace in data:
        data_duration += trace.stats.delta * trace.stats.npts
        if trace.stats.starttime.datetime < start:
            starttime = trace.stats.starttime
            start = trace.stats.starttime.datetime
        if trace.stats.endtime.datetime > stop:
            endtime = trace.stats.endtime
            stop = trace.stats.endtime.datetime

    net = trace.stats.network.upper()
    sta = trace.stats.station.upper()
    comp = trace.stats.channel.upper()
    path = folder.replace('\\', '/')
    starttime = starttime.datetime.replace(microsecond=0)
    endtime = endtime.datetime.replace(microsecond=0)

    return api.update_data_availability(db, net, sta, comp, path, basename,
            starttime, endtime, data_duration, gaps_duration,
            data[0].stats.sampling_rate)


def process_stream(db, folder, basename, stream, id_, startdate, enddate,
                   goal_sampling_rate):
    """
    Reads a data stream from an archive file and update the availability in the
    database. Returns the return code of the update_availability function.

    :param db: the sqlalchemy session object
    :param folder: the directory where lies the miniseed file
    :param basename: the basename of the miniseed file to read data from
    :param stream: unfiltered obspy.core.stream.Stream object
    :param id_: the id of the stream to read
    :param startdate: the startdate configuration value
    :param enddate: the enddate configuration value
    :param goal_sampling_rate: the sampling rate
    """
    data = stream.select(id=id_)
    if data[-1].stats.endtime.date < startdate:
        logger.debug('ignoring %s: before start date!' % id_)
        return 0
    if data[0].stats.starttime.date > enddate:
        logger.debug('ignoring %s: after end date!' % id_)
        return 0
    if data[0].stats.sampling_rate < goal_sampling_rate:
        logger.debug('ignoring %s: sampling rate smaller than '
                     'CC sampling rate' % id_)
        return 0
    return update_availability(db, folder, basename, data)


def scan_data_files(db, folder, files, startdate, enddate, goal_sampling_rate,
        archive_format, logger):
    """
    Processes a list of files from a folder, and update the data availability
    table in the database whenever their data matches our dates and sampling
    rate parameters.

    :param db: the sqlalchemy session object.
    :param folder: the directory where lies the miniseed files.
    :param files: the list of files of the folder to read.
    :param stream: unfiltered obspy.core.stream.Stream object.
    :param id_: the id of the stream to read.
    :param startdate: the startdate configuration value.
    :param enddate: the enddate configuration value.
    :param goal_sampling_rate: the sampling rate.
    :paran logger: the logger instance to use for logging.
    """
    added = 0
    modified = 0
    unchanged = 0
    for basename in files:
        pathname = os.path.join(folder, basename)
        # logger.debug('reading file %s' % pathname)
        try:
            # Note: if format is None or unknown, obspy will use auto-detection.
            # See https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html
            stream = obspy.core.read(pathname, headonly=True, format=archive_format or None)
            for id in set([t.id for t in stream]):
                update_rv = process_stream(db, folder, basename, stream, id,
                                           startdate, enddate,
                                           goal_sampling_rate)
                if update_rv == 1:
                    added += 1
                elif update_rv == -1:
                    modified += 1
                else:
                    unchanged += 1
        except obspy.io.mseed.ObsPyMSEEDFilesizeTooSmallError as e:
            logger.warning("Ignoring possible empty file '%s'."
                           ' Got error %s: %s' %
                           (pathname, e.__class__.__name__, str(e)))
        except OSError as e:
            # This should catch errors about file accesses
            logger.error("Error while processing file '%s': %s" %
                         (pathname, str(e)))
            # db.close()
            # Re-raise the exception to end the scan
            # raise
        except Exception as e:
            logger.error("Error while processing file '%s': %s" %
                         (pathname, str(e)))
    logger.info('%s: Added %i | Modified %i | Unchanged %i' %
                (folder, added, modified, unchanged))


def list_directory(folder, mintime):
    """
    Builds and returns a list of files from a given directory that were
    modified after the mintime epoch (or all file in this directory if mintime
    is None).

    :param folder: the directory where lies the miniseed files
    :param mintime: the datetime.datetime object representing the minimal
        modification time of miniseed files that should be analysed
    """
    files = []
    try:
        for entry in scandir(folder):
            if not entry.is_file() or entry.name.startswith('.'):
                # silently ignore directories and files starting with '.'
                continue
            if mintime is None or entry.stat().st_mtime >= mintime:
                files.append(entry.name)
    except OSError as e:
        logger.error('Error while reading folder %s (%s)'
                     % (folder, str(e)))
        # Re-raise the exception to end the scan
        # raise
    except Exception as e:
        logger.error('Error while reading folder %s (%s)'
                     % (folder, str(e)))
        # Re-raise the exception to end the scan
    return files


def scan_folders(folders, mintime, startdate, enddate, goal_sampling_rate,
        archive_format, loglevel=None):
    """
    Reads files in a list of folders and updates their data availability in
    database, silently ignoring non-matching files and empty folders.
    If mintime is not None, only files modified since the mintime epoch will be
    considered.

    :param db: the sqlalchemy session object.
    :param folder: the directory where lies the miniseed files.
    :param files: the list of files of the folder to read.
    :param stream: unfiltered obspy.core.stream.Stream object.
    :param id_: the id of the stream to read.
    :param startdate: the startdate configuration value.
    :param enddate: the enddate configuration value.
    :param goal_sampling_rate: the sampling rate.
    :param loglevel: the loglevel to use in the children processes (if None,
        use the global logger instance).
    """
    global logger
    # Reconfigure logger to show the pid number in log records
    logger = api.get_logger('msnoise.scan_archive_child', logger.level,
                            with_pid=True)
    db = api.connect()
    for folder in folders:
        logger.debug('scanning dir %s' % folder)
        if not os.path.isdir(folder):
            logger.warning('Ignoring untidy file %s' % folder)
            continue

        # Get the files of the directory matching our requirements
        files = list_directory(folder, mintime)
        if mintime is not None:
            debug_msg = 'Found %d files with mtime >= %d in %s' \
                            % (len(files), mintime, folder)
        else:
            debug_msg = 'Found %d files in %s' % (len(files), folder)
        logger.debug(debug_msg)
        if files:
            scan_data_files(db, folder, files, startdate, enddate,
                            goal_sampling_rate, archive_format, logger)
        # else: no matching files found in this directory, nothing to do
    db.close()


def get_archives_folders(data_folder, data_structure,
                         years, stations, channels):
    """
    Builds and returns the list of directories of the archive to scan for data.

    Uses the data_structure description string (in the form
    'NET/STA/YEAR/NET.STA.YEAR.DAY.MSEED') and replaces each element to build
    the directory list.

    :param data_folder: the top directory of the data files
    :param data_structure: the data_structure format string
    :param years: the years to process
    :param stations: the stations to process
    :param channels: the channel to process
    """
    folders = []
    for year in years:
        for channel in channels:
            stafolder = os.path.dirname(data_structure) \
                        .replace('YEAR', '%04i' % year) \
                        .replace('DAY', '*') \
                        .replace('HOUR', '*') \
                        .replace('CHAN', channel) \
                        .replace('TYPE', 'D') \
                        .replace('LOC', '*')
            for sta in stations:
                folders.append(os.path.join(
                    data_folder,
                    stafolder.replace('NET', sta.net) \
                             .replace('STA', sta.sta)))
    # returns a set to make sure we only get unique folder globs
    # (to be verified: could there really be duplicates?)
    return set(folders)


def get_data_structure(config_data_st):
    """
    Returns the data structure description string according to the
    'data_structure' configuration variable.

    This variable holds either a pre-defined template name (defined in the
    data_structures.py module), or the description string itself (in a form
    similar to 'YEAR/NET/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY').

    If none of these formats are found, and attempt is made to read the
    'rawpath' variable from a 'custom.py' module in the current directory.

    :param config_sta_st: the value of the data_structure setting in the
        configuration

    Returns None if the description string cannot be found.
    """
    if config_data_st in data_structures.data_structure:
        return data_structures.data_structure[config_data_st]

    if config_data_st.count('/') != 0:
        return config_data_st

    logger.info("Can't parse the archive for format %s !" % config_data_st)
    logger.info('trying to import local parser (should return a station list)')
    # Try to load custom.py
    # warning: this could lead to load any 'custom' module on sys.path...
    if os.path.isfile(os.path.join(os.getcwd(), 'custom.py')):
        try:
            sys.path.insert(0, os.getcwd())
            from custom import data_structure as rawpath
            return rawpath
        except ImportError:
            pass
    return None


def spawn_processes(pool, nproc, dir_list, scan_func, scan_args):
    """
    Spawn nproc processes that will execute scan_func(scan_args).

    Return a list of multiprocessing.pool.AsyncResult objects that
    represent the list of children to wait for.
    """
    # Get the ceiling of (len(dir_list) / nproc) (w/o importing math.ceil)
    n = int((len(dir_list) + nproc - 1) / nproc)

    children = []
    for i in range(nproc):
        folder_slice = dir_list[n*i:n*(i+1)]
        logger.debug('Spawning a child to process {} of the {} directories.'
                     .format(len(folder_slice), len(dir_list)))
        children.append(pool.apply_async(scan_func,
            [folder_slice] + list(scan_args)))
    return children


def await_children(pool, children):
    """
    Wait for children to return.

    Watch multiprocessing.pool.AsyncResult object in 'children' until all
    children are finished.  If any of the children raise an exception or
    exits anormally, terminate the pool.  Re-raise any exception raised by a
    child.
    """
    try:
        while children:
            # Wait a few seconds for each child to return but don't block
            # so we catch any child error early and terminate all children.
            # Note: once python 2 support in msnoise will be dropped,
            # consider replacing this loop by error_callback arguments to
            # apply_async (unavailable in python 2).
            for result in children:
                try:
                    rc = result.get(timeout=1)
                    # The child has terminated without raising any exception:
                    # remove it from the list of children to wait for
                    children.remove(result)
                except multiprocessing.TimeoutError:
                    # AsyncResult.get() reached its wait timeout: do nothing
                    pass
    except Exception:
        # AsyncResult.get() re-raised an exception raised in the child
        logger.debug('A child process has raised an exception:'
                     ' terminating and re-raising it.')
        pool.terminate()
        pool.join()
        raise  # this re-raises the last active exception


def scan_archive(folder_globs, nproc, mintime, startdate, enddate,
                 goal_sampling_rate, archive_format):
    """
    For each files in the archive folders, fork a process to read its data, and
    updates the availibility in the database. If mintime is not None, only
    consider files whose 'mtime' timestamp is equal of older than mintime (in
    seconds since 01/01/1970).

    :param folder_globs: the glob string designating the folders to scan
    :param nproc: the maximum number of concurrent process to spawn to scan the
        files
    :param mintime: the datetime.datetime object representing the minimal
        modification time of miniseed files that should be analysed
    :param startdate: the startdate configuration value
    :param enddate: the enddate configuration value
    :param goal_sampling_rate: the sampling rate
    """
    # develop the list of directories matching the globs
    logger.debug('Building directory list from archive folders {}...'
                 .format(', '.join(folder_globs)))
    dir_list = []
    for f in sorted(folder_globs):
        dir_list.extend(glob.glob(f))

    logger.info('Scanning {} directories...'.format(len(dir_list)))
    if nproc == 1:
        # In single process mode, we simply call scan_folder()
        scan_folders(dir_list, mintime, startdate, enddate, goal_sampling_rate,
                     archive_format)
    else:
        # In multiprocessing mode, we split the folders into nproc lists of
        # similar size and have them processed by as many child processes.
        # (This reduces the number of process spawning, as it is a slow
        # operation.)

        # Note: consider using a context manager instead of try/except
        # once python 2 support will be dropped in msnoise:
        # with multiprocessing.Pool(processes=nproc) as pool:
        pool = multiprocessing.Pool(processes=nproc)

        # Launch nproc children working on a mostly equal number of folders
        children = spawn_processes(pool, nproc, dir_list, scan_folders,
                (mintime, startdate, enddate, goal_sampling_rate,
                 archive_format, logger.level))
        # Wait for children to finish, or terminate them if one crashes
        await_children(pool, children)


def parse_crondays(crondays):
    """
    Parse the 'crondays' option and return it as a datetime.timedelta object.
    """
    try:
        crondays = float(crondays)
    except ValueError:
        # The provided value is not a float: it must match '[Xw][Xd][Xh]'
        # (with optional blank characters between groups)
        match = re.search('^(?:(\d+)w\s*)?(?:(\d+)d\s*)?(?:(\d+)h\s*)?$', crondays)
        if not match:
            raise FatalError("Unrecognized format for "
                    "configuration parameter 'crondays'")
        delta_days = 7 * int(match.group(1) or 0) + int(match.group(2) or 0)
        delta_seconds = 3600 * int(match.group(3) or 0)
    else:
        # Use abs() to accept negative values for
        # backward compatibility with MSNoise < 1.6.
        delta_days = abs(crondays)
        delta_seconds = 0
    return datetime.timedelta(days=delta_days, seconds=delta_seconds)


def main(init=False, threads=1, crondays=None, forced_path=None,
         forced_path_recursive=True):
    """
    Update data availibility information from modified miniseed files.

    Scans the archives for miniseed files recently modified (or all if init is
    true) and updates data availaibility information in the database.

    :param init: scan all files in the archive, disregarding their modification
        time (if not None, crondays is ignored).
    :param threads: the maximum number of threads/processes to use while
        scanning the files.
    :para crondays: override the 'crondays' configuration value for this run to
        limit the scan to files modified since the designated relative date in
        the past (according to the file's modification time); if None, use the
        value from the configuration.
    :para forced_path: scan files from this path instead of the configured
        archive path.
    :para forced_path_recursive: whether to also scan files from subdirectories
        of forced_path.
    """
    scanning_starttime = time.time()
    logger.info('*** Starting: Scan Archive ***')
    db = api.connect()

    # if api.read_db_inifile().tech == 1 and threads != 1:
    #     logger.warning('You can not work on %i threads because SQLite only'
    #                    ' supports 1 connection at a time. Continuing with '
    #                    ' 1 process only.' % threads)
    #     threads = 1
    logger.info('Will work on %i thread(s)' % threads)

    if init:
        logger.info('Initializing: updating availability using the whole'
                    ' archive (should be run only once)')
        modification_delta = 0
    else:
        if crondays is None:
            crondays = api.get_config(db, 'crondays')
        modification_delta = parse_crondays(crondays)
        logger.info('Updating availability: scanning the archive for files'
                    ' modified {} ago or less.'.format(modification_delta))

    startdate = datetime.datetime.strptime(
            api.get_config(db, 'startdate'),'%Y-%m-%d').date()
    enddate = datetime.datetime.strptime(
            api.get_config(db, 'enddate'), '%Y-%m-%d').date()
    archive_format = api.get_config(db, 'archive_format')
    goal_sampling_rate = float(api.get_config(db, 'cc_sampling_rate'))
    search_info_log = 'Will search for files between {} and {}'\
                      .format(startdate, enddate)
    if archive_format:
        search_info_log += ", forcing format '{}'".format(archive_format)
    logger.info(search_info_log + '.')

    if init:
        mintime = None
    else:
        # Note: avoid datetime.timestamp() below as it is python3 only *and*
        # does not work correctly with naive datetime representing UTC.
        # (See the official doc for datetime.timestamp())
        mintime = (datetime.datetime.utcnow() - modification_delta
                   - datetime.datetime(1970, 1, 1, 0, 0)).total_seconds()

    if forced_path is None:
        data_folder = os.path.realpath(api.get_config(db, 'data_folder'))
        channels = api.get_config(db, 'channels').split(',')
        logger.debug('Will search for channels: %s' % channels)


        rawpath = get_data_structure(api.get_config(db, 'data_structure'))
        if rawpath is None:
            raise FatalError("Cannot read configured data_structure '%s' anywhere "
                    "(tried file custom.py in folder '%s')."
                    % (api.get_config(db, 'data_structure'), os.getcwd()))
            return

        if not os.path.isdir(data_folder):
            raise FatalError("Cannot find directory '{}'. Aborting."
                             .format(data_folder))
        folders_to_glob = get_archives_folders(
                data_folder,
                rawpath,
                range(startdate.year,
                      min(datetime.datetime.utcnow().year, enddate.year) + 1),
                api.get_stations(db, all=False),
                channels)
        #logger.debug('Folders to glob: %s' % ','.join(folders_to_glob))
    elif forced_path_recursive:
        # Scan directory and all subdirectories in the tree
        folders_to_glob = [d[0] for d in walk(forced_path)]
    else:
        # Only scan the forced directory
        folders_to_glob = [forced_path]

    # Close the db connection as we don't need it in this process any more.
    db.close()

    # Run the main scan
    try:
        scan_archive(folders_to_glob, threads, mintime, startdate, enddate,
                     goal_sampling_rate, archive_format)
    except Exception as e:
        logger.critical('Scan aborted because the following error occured '
                        'while scanning the archive:\n{}'.format(e))
        traceback.print_exc()
    else:
        logger.info('*** Finished: Scan Archive ***')
        logger.info('It took %.2f seconds' % (time.time() - scanning_starttime))

