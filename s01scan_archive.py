"""One advantage of MSNoise is its ability to be used as an automated
monitoring tool. To run every night on the data acquired during the
previous day, MSNoise needs to check the data archive for new or modified
files. Those files could have been acquired during the last day, but be data of
a previously offline station and contain useful information for, say, a month
ago. The time to search for is defined in the config. The scan_archive script
uses the find command (gnufind on Windows) with the -mtime argument to locate
new or modified files. Once located, they are inserted (if new) or updated (if
modified) in the data availability table.

To run the code on two Process, execute the following in console:

.. code-block:: sh

    python s01scan_archive.py -t 2

Special case: first run
~~~~~~~~~~~~~~~~~~~~~~~~

This script is the same as for the routine, but one has to pass the init
argument:

.. code-block:: sh

    python s01scan_archive.py --init -t 2

This will scan the data_archive folder the configured stations and will insert
all files found in the data_availability table in the database. As usual,
calling the script with a --help argument will show its usage.
"""
from obspy.core import read
import glob
import os
import datetime
import time
import logging
import logging.handlers
from subprocess import Popen, PIPE
from multiprocessing import Process
import multiprocessing
import argparse

from database_tools import *
from data_structures import data_structure


def worker(files, folder, startdate, enddate):
    # logger = logging.getLogger('worker')

    db = connect()
    for file in files:
        file = os.path.join(folder, file)
        try:
            r0 = time.time()
            name = os.path.split(file)[1]
            data = read(file, headonly=True)
            # print data
            if data[0].stats.starttime.date < startdate:
                r2 = time.time()
                logging.info(
                    '%s: Before Start-Date! (%.2f)' % (name, r2 - r0))
            elif data[-1].stats.endtime.date > enddate:
                r2 = time.time()
                logging.info('%s: After End-Date! (%.2f)' % (name, r2 - r0))
            else:
                gaps = data.getGaps()
                gaps_duration = 0
                for gap in gaps:
                    gaps_duration += gap[6]
                data_duration = 0
                start = datetime.datetime.strptime('2100-01-01', '%Y-%m-%d')
                stop = datetime.datetime.strptime('1900-01-01', '%Y-%m-%d')
                for trace in data:
                    data_duration += trace.stats.delta * trace.stats.npts
                    if trace.stats.starttime.datetime < start:
                        starttime = trace.stats.starttime
                        start = trace.stats.starttime.datetime
                    if trace.stats.endtime.datetime > stop:
                        endtime = trace.stats.endtime
                        stop = trace.stats.endtime.datetime

                net = trace.stats.network
                sta = trace.stats.station
                comp = trace.stats.channel
                path = folder.replace('\\', '/')
                r1 = time.time()

                result = update_data_availability(
                    db, net, sta, comp, path, name, starttime.datetime,
                    endtime.datetime, data_duration, gaps_duration,
                    data[0].stats.sampling_rate)

                r2 = time.time()
                if result:
                    logging.info(
                        'Added: "%s" (read:%.2f (%.2f) seconds | save:%.4f seconds)' %
                        (name, r1 - r0, r2 - r0, r2 - r1))
                else:
                    logging.info(
                        'Already Exists: "%s" (read:%.2f (%.2f) seconds | save:%.4f seconds)' %
                        (name, r1 - r0, r2 - r0, r2 - r1))
        except Exception as e:
            print "Problem", e
        db.close()

if __name__ == "__main__":
    t = time.time()
    parser = argparse.ArgumentParser(description='Scan the data archive and insert the\
    metadata in the database')
    parser.add_argument('-i', '--init', action="store_true",
                        help='Initialize the archive: should only be done upon first run.\
                        Will read all files in the archive that match the station/component\
                        (check that)',
                        default=False)
    parser.add_argument('-t', '--threads',
                        help='Number of parellel threads to use [default:1]', default=1, type=int)
    args = parser.parse_args()

    # rootLogger = logging.getLogger('')
    # rootLogger.setLevel(logging.DEBUG)
    # socketHandler = logging.handlers.SocketHandler('localhost',
                        # logging.handlers.DEFAULT_TCP_LOGGING_PORT)
    # rootLogger.addHandler(socketHandler)
    # global logger
    # logger = logging.getLogger('main')

    multiprocessing.log_to_stderr()
    global logger
    logger = multiprocessing.get_logger()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    for handler in logger.handlers:
        handler.setFormatter(formatter)

    logger.info('*** Starting: Scan Archive ***')
    db = connect()

    init = False
    mtime = -2

    if args.init:
        logger.info("Initializing (should be run only once)")
        mtime = "-20000"
        init = True
    else:
        mtime = "%s" % mtime

    nthreads = 1
    if args.threads:
        nthreads = args.threads
    if get_tech() == 1:
        logger.info("You can not work on %i threads because SQLite only\
 supports 1 connection at a time" % nthreads)
        nthreads = 1

    logger.info("Will work on %i threads" % nthreads)

    if os.name == "nt":
        find = "gnufind"
    else:
        find = "find"
    startdate = get_config(db, 'startdate')
    startdate = datetime.datetime.strptime(startdate, '%Y-%m-%d').date()
    enddate = get_config(db, 'enddate')
    enddate = datetime.datetime.strptime(enddate, '%Y-%m-%d').date()

    data_folder = get_config(db, 'data_folder')
    data_struc = get_config(db, 'data_structure')
    channels = [c for c in get_config(db, 'channels').split(',')]

    folders_to_glob = []
    rawpath = data_structure[data_struc]
    for year in range(startdate.year, min(datetime.datetime.now().year, enddate.year) + 1):
        for channel in channels:
            stafol = os.path.split(rawpath)[0].replace('YEAR', "%04i" % year).replace('DAY', '*').replace(
                'HOUR', '*').replace('CHAN', channel).replace('TYPE', '*').replace('LOC', '*')
            for sta in get_stations(db, all=False):
                tmp = os.path.join(data_folder, stafol.replace(
                    'NET', sta.net).replace('STA', sta.sta))
                folders_to_glob.append(os.path.join(data_folder, tmp))

    clients = []
    for fi in sorted(folders_to_glob):
        folders = glob.glob(fi)
        for folder in sorted(folders):
            if init:
                files = os.listdir(folder)
            else:
                proc = Popen(
                    [find, folder, "-type", "f", "-mtime", mtime, "-print"], stdout=PIPE, stderr=PIPE)
                stdout, stderr = proc.communicate()

                if len(stdout) != 0:
                    files = sorted(stdout.split('\n'))
                else:
                    files = []

            if '' in files:
                files.remove('')

            if len(files) != 0:
                logger.info('Started: %s'%folder)
                client = Process(target=worker, args=([files,folder,startdate,enddate]))
                client.start()
                clients.append(client)
            while len(clients) >= nthreads:
                for client in clients:
                    client.join(0.01)
                    if not client.is_alive():
                        client.join(0.01)
                        clients.remove(client)
                
    while len(clients) != 0:
        for client in clients:
            client.join(0.01)
            if not client.is_alive():
                client.join(0.01)
                clients.remove(client)

    logger.info('*** Finished: Scan Archive ***')
    logger.info('It took %.2f seconds' % (time.time() - t))
