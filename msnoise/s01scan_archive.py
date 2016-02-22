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

    $ msnoise -t 2 scan_archive

Special case: first run
~~~~~~~~~~~~~~~~~~~~~~~~

This script is the same as for the routine, but one has to pass the init
argument:

.. code-block:: sh

    $ msnoise -t 2 scan_archive --init

This will scan the data_archive folder the configured stations and will insert
all files found in the data_availability table in the database. As usual,
calling the script with a --help argument will show its usage.

"""
from obspy.core import read
import glob
import os
import sys
import datetime
import time
import logging
import logging.handlers
from subprocess import Popen, PIPE
from multiprocessing import Process
import multiprocessing
import multiprocessing_logging
import numpy as np

import argparse
import traceback



from .api import *
from .data_structures import data_structure

def worker(files, folder, startdate, enddate, goal_sampling_rate):
    import logging
    logging = logging.getLogger("worker-logger")
    db = connect()
    added = 0
    modified = 0
    for file in files:
        file = os.path.join(folder, file)
        try:
            name = os.path.split(file)[1]
            data = read(file, headonly=True)
            if data[0].stats.starttime.date < startdate:
                logging.info(
                    '%s: Before Start-Date!' % (name))
            elif data[-1].stats.endtime.date > enddate:
                logging.info('%s: After End-Date!' % (name))
            elif data[0].stats.sampling_rate < goal_sampling_rate:
                logging.info("%s: Sampling rate smaller than CC sampling rate"
                             % name)
            else:
                gaps = data.get_gaps()
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

                net = trace.stats.network.upper()
                sta = trace.stats.station.upper()
                comp = trace.stats.channel.upper()
                path = folder.replace('\\', '/')

                starttime = starttime.datetime.replace(microsecond=0)
                endtime = endtime.datetime.replace(microsecond=0)
                result = update_data_availability(
                    db, net, sta, comp, path, name, starttime,
                    endtime, data_duration, gaps_duration,
                    data[0].stats.sampling_rate)

                if result:
                    added += 1
                else:
                    modified += 1

        except Exception as e:
            logging.debug("ERROR: %s", e)
    db.close()
    logging.debug("%s: Added %i | Modified %i", str(folder), added, modified)
    return

def main(init=False, threads=1):
    logger = logging.getLogger('')
    # logger.addHandler(
    # multiprocessing_logging.MultiProcessingHandler('worker-logger'))

    t = time.time()

    logging.info('*** Starting: Scan Archive ***')
    db = connect()

    mtime = -2

    if init:
        logging.info("Initializing (should be run only once)")
        mtime = "-20000"
        init = True
    else:
        mtime = "%s" % mtime

    nthreads = 1
    if threads:
        nthreads = threads
    if get_tech() == 1 and nthreads > 1:
        logging.info("You can not work on %i threads because SQLite only\
 supports 1 connection at a time" % nthreads)
        nthreads = 1

    logging.info("Will work on %i threads" % nthreads)

    if os.name == "nt":
        find = "find"
    else:
        find = "find"
    startdate = get_config(db, 'startdate')
    startdate = datetime.datetime.strptime(startdate, '%Y-%m-%d').date()
    enddate = get_config(db, 'enddate')
    enddate = datetime.datetime.strptime(enddate, '%Y-%m-%d').date()
    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))

    data_folder = get_config(db, 'data_folder')
    data_struc = get_config(db, 'data_structure')
    channels = [c for c in get_config(db, 'channels').split(',')]
    folders_to_glob = []
    if data_struc in data_structure.keys():
        rawpath = data_structure[data_struc]
    else:
        print("Can't parse the archive for format %s !" % data_struc)
        print("trying to import local parser (should return a station list)")
        print()
        if os.path.isfile(os.path.join(os.getcwd(), 'custom.py')):
            sys.path.append(os.getcwd())
            from custom import data_structure as rawpath
            rawpath=data_structure[data_struc]
        else:
            print("No file named custom.py in the %s folder" % os.getcwd())
            return
    for year in range(startdate.year, min(datetime.datetime.now().year, enddate.year) + 1):
        for channel in channels:
            stafol = os.path.split(rawpath)[0].replace('YEAR', "%04i" % year).replace('DAY', '*').replace(
                'HOUR', '*').replace('CHAN', channel).replace('TYPE', '*').replace('LOC', '*')
            for sta in get_stations(db, all=False):
                tmp = os.path.join(data_folder, stafol.replace(
                'NET', sta.net).replace('STA', sta.sta))
                folders_to_glob.append(tmp)
    folders_to_glob = np.unique(folders_to_glob)

    clients = []
    for fi in sorted(folders_to_glob):
        folders = glob.glob(fi)
        for folder in sorted(folders):
            if init:
                files = os.listdir(folder)
            else:
                proc = Popen(
                    [find, folder, "-type", "f", "-mtime", mtime, "-print"],
                    stdout=PIPE, stderr=PIPE)
                stdout, stderr = proc.communicate()

                if len(stdout) != 0:
                    files = sorted(stdout.replace('\r', '').split('\n'))
                else:
                    files = []

            if '' in files:
                files.remove('')

            if len(files) != 0:
                logging.info('%s: Started' % folder)
                client = Process(target=worker, args=([files, folder,
                                                       startdate, enddate,
                                                       goal_sampling_rate]))
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

    logging.info('*** Finished: Scan Archive ***')
    logging.info('It took %.2f seconds' % (time.time() - t))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Scan the data archive and insert the\
    metadata in the database')
    parser.add_argument('-i', '--init', action="store_true",
                        help='Initialize the archive: should only be done upon first run.\
                        Will read all files in the archive that match the station/component\
                        (check that)',
                        default=False)
    parser.add_argument('-t', '--threads',
                        help='Number of parellel threads to use [default:1]',
                        default=1, type=int)
    args = parser.parse_args()
    main(args.init, args.threads)
