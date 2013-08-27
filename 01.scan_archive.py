from obspy.core import read
import glob,sys,os,datetime
import time
import logging
import threading
from subprocess import Popen, PIPE

from database_tools import *
from data_structures import data_structure

class ActivePool(object):
    def __init__(self):
        super(ActivePool, self).__init__()
        self.active = []
        self.lock = threading.Lock()
    def makeActive(self, name):
        with self.lock:
            self.active.append(name)
            # logging.debug('Running: %s', self.active)
    def makeInactive(self, name):
        with self.lock:
            self.active.remove(name)
            # logging.debug('Finished: %s', self.active)

def worker(s, pool):
    # logging.debug('Waiting to join the pool')
    with s:
        Fname = threading.currentThread().getName()
        folder = threading.currentThread().folder
        pool.makeActive(Fname)
        db = connect()
        try:
            r0 = time.time()
            source = Fname
            name = os.path.split(source)[1]
            data = read(source,quality=False)
            # print data
            if data[0].stats.starttime.date < startdate:
                r2 = time.time()
                logging.debug('%s: Before Start-Date! (%.2f)'%(name,r2-r0))
            elif data[-1].stats.endtime.date > enddate:
                r2 = time.time()
                logging.debug('%s: After End-Date! (%.2f)'%(name,r2-r0))
            else:
                gaps = data.getGaps()
                gaps_duration = 0
                for gap in gaps:
                    gaps_duration += gap[6]
                data_duration = 0
                start  = 1e20
                stop = 0

                for trace in data:
                    data_duration += trace.stats.delta * trace.stats.npts
                    if trace.stats.starttime.timestamp < start:
                        starttime = trace.stats.starttime
                        start = trace.stats.starttime.datetime
                    if trace.stats.endtime.timestamp > stop:
                        endtime = trace.stats.endtime
                        stop = trace.stats.endtime.datetime
                
                net = trace.stats.network
                sta = trace.stats.station
                comp = trace.stats.channel
                path=folder.replace('\\','/')
                r1 = time.time()
                
                result = update_data_availability(db, net, sta, comp, path, name, starttime.datetime, endtime.datetime,data_duration, gaps_duration, data[0].stats.sampling_rate)
                
                r2 = time.time()
                if result :
                    logging.debug('Added: "%s" (read:%.2f (%.2f) seconds | save:%.4f seconds)'% (name,r1-r0,r2-r0,r2-r1))
                else:
                    logging.debug('Already Exists: "%s" (read:%.2f (%.2f) seconds | save:%.4f seconds)'% (name,r1-r0,r2-r0,r2-r1))
        except Exception as e:
            print "Problem", e
        pool.makeInactive(Fname)
        db.close()

if __name__ == "__main__":
    db = connect()

    logging.basicConfig(level=logging.DEBUG,
                        filename="./scan_archive_threaded.log",
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        filemode='w')

    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    logging.info('*** Starting: Scan Archive ***')
    
    init = False
    mtime = get_config(db,"crondays")
    print "mtime:",mtime
    if len(sys.argv) > 1:
        if sys.argv[1] == 'init':
            print "init man!"
            mtime = "-20000"
            init = True
    else:
        mtime = "%s"%mtime

    if os.name == "nt":
        find = "gnufind"
    else:
        find = "find"
    startdate = get_config(db,'startdate')
    startdate = datetime.datetime.strptime(startdate,'%Y-%m-%d').date()
    enddate = get_config(db,'enddate')
    enddate = datetime.datetime.strptime(enddate,'%Y-%m-%d').date()

    data_folder = get_config(db,'data_folder')
    data_struc = get_config(db,'data_structure')
    channels = [c for c in get_config(db,'channels').split(',')]
    
    folders_to_glob = []
    rawpath = data_structure[data_struc]
    for year in range(startdate.year, min(datetime.datetime.now().year , enddate.year) + 1):
        for channel in channels:
            stafol = os.path.split(rawpath)[0].replace('YEAR',"%04i"%year).replace('DAY','*').replace('HOUR','*').replace('CHAN',channel).replace('TYPE','*').replace('LOC','*')
            for sta in get_stations(db, all=False):
                tmp = os.path.join( data_folder, stafol.replace('NET',sta.net).replace('STA',sta.sta))
                folders_to_glob.append(os.path.join( data_folder, tmp))
    
    pool = ActivePool()
    s = threading.Semaphore(5)
    for fi in sorted(folders_to_glob):
        folders = glob.glob(fi)
        for folder in sorted(folders):
            if init:
                proc = Popen(["ls", "-1",folder], stdout=PIPE, stderr=PIPE)
            else:
                proc = Popen([find,folder,"-type","f", "-mtime", mtime, "-print"], stdout=PIPE, stderr=PIPE)
            
            stdout,stderr = proc.communicate()

            if len(stdout) != 0:
                sources = sorted(stdout.split('\n'))
                for source in sources:
                    if len(source) != 0:
                        source = os.path.join(folder,source) #.replace('\r','').replace('\n','').replace('\\','/')
                        t = threading.Thread(target=worker, name=source, args=(s, pool))
                        t.folder = folder
                        t.start()
    while threading.activeCount() != 1:
        time.sleep(0.1)
    logging.info('*** Finished: Scan Archive ***')