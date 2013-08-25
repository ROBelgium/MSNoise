from database_tools import *
import logging
logging.basicConfig(level=logging.DEBUG,
                    filename="./new_jobs.log",
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)


logging.info('*** Starting: New Jobs ***')

db = connect()

if get_config(db, name="autocorr") in ['Y','y','1',1]:
    AUTOCORR = True
else:
    AUTOCORR = False

stations_to_analyse = [sta.sta for sta in get_stations(db,used=True)]
new_files = get_new_files(db)

days = {}
old_day = 0
old_pair = ""
day_pairs = []
for new_file in new_files:
    Mnet, Msta, comp, starttime, endtime, data_duration = new_file
    # logging.debug('%s.%s will be MASTER for %s-%s'% (Mnet, Msta, starttime, endtime))
    if Msta in stations_to_analyse:
        day = "%s" % (starttime.date())
        # logging.debug("%s : %s" %(day, Msta))
        if day != old_day:
            day_pairs = np.unique(day_pairs)
            for pair in day_pairs:
                logging.debug('New Job for: %s - %s'%(day,pair))
                add_job(db, day, pair,type='CC')
            day_pairs = []
            old_day = day
        
        # print day,
        # print '%s.%s' % (net,sta),
        
        available_stations = []
        for station in get_data_availability(db, starttime=starttime.date(), endtime=endtime.date()):
            net, sta, comp, datetime, endtime, duration, gaps, samplerate, flag = station
            if sta in stations_to_analyse:
                if '%s.%s' % (net,sta) not in available_stations:
                    available_stations.append('%s.%s'%(net,sta))
        # print ":", ','.join(available_stations)
        stations = np.array([])
        pairs = []
        nS = '%s.%s' % (Mnet,Msta)
        i = 0
        for aS in available_stations:
            if not AUTOCORR and nS == aS:
                pass
            else:                    
                if i == 0:
                    pairs = np.array(':'.join(sorted([nS,aS])))
                    # stations = np.array([nS, aS])
                    i+=1
                else:
                    pairs = np.vstack((pairs,':'.join(sorted([nS,aS]))))
                    # stations = np.append(stations, np.array([nS,aS]))

        pairs = np.unique(pairs)
        for pair in pairs:
            day_pairs.append(pair)
# update all _data_availability and mark files as "A"rchives
for sta in get_stations(db,used=True):
    mark_data_availability(db,condition='net="%s" and sta="%s"'%(sta.net, sta.sta), flag='A')
    
logging.info('*** Finished: New Jobs ***')