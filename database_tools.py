## Database Tools
import MySQLdb as SQL
import cPickle
import numpy as np
import os
from stat import *
from collections import OrderedDict

from obspy.core import Stream, Trace, read
from obspy.sac import SacIO
from obspy.core.util import gps2DistAzimuth

import datetime
import itertools

#Data Structure Definitions:
data_structure = {}
data_structure['SDS'] = "YEAR/NET/STA/CHAN.D/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"
data_structure['BUD'] = "NET/STA/STA.NET.LOC.CHAN.YEAR.DAY"
data_structure['IDDS'] = "YEAR/NET/STA/CHAN.TYPE/DAY/NET.STA.LOC.CHAN.TYPE.YEAR.DAY.HOUR"
data_structure['PDF'] = "YEAR/STA/CHAN.TYPE/NET.STA.LOC.CHAN.TYPE.YEAR.DAY"

default = OrderedDict()
default['data_folder'] = ["Data Folder",'']
default['output_folder'] = ["CC Output Folder",'CROSS_CORRELATIONS']
default['data_structure'] = ["Data Structure [SDS]/BUD/IDDS",'SDS']
default['network'] = ["Network to analyse [*]",'*']
default['channels'] = ["Channels need to match the value (ex: [*], *Z, BH*, HHZ,...)",'*']

default['startdate'] = ["Start Date to process: [1970-01-01]='since beginning of the archive'","1970-01-01"]
default['enddate'] = ["End Date to process: [2100-01-01]='No end'","2100-01-01"]

default['analysis_duration'] = ["Duration of the Analysis (total in seconds : 3600, [86400])",'86400']
default['cc_sampling_rate'] = ["Sampling Rate for the CrossCorrelation [20.0]",'20.0']
default['resampling_method'] = ["Resampling method [Resample]/Decimate",'Resample']
default['decimation_factor'] = ["If Resampling mether=Decimate, decimation factor [5]",'5']
default['preprocess_lowpass'] = ["Preprocessing Low-pass value in Hz [8.0]",'8.0']
default['preprocess_highpass'] = ["Preprocessing High-pass value in Hz [0.01]",'0.01']


default['maxlag'] = ["Maximum lag (in seconds) [120.0]",'120.']
default['corr_duration'] = ["Data windows to correlate (in seconds) [1800.]",'1800.']
default['windsorizing'] = ["Windsorizing at N time RMS (in unit), 0 disables windsorizing [3]",'3']

default['crondays'] = ["Number of days to monitors with cron [-1]",'-1']


default['ZZ'] = ["Compute ZZ correlation [Y]/N",'Y']
default['ZR'] = ["Compute ZR correlation [Y]/N",'Y']
default['ZT'] = ["Compute ZT correlation [Y]/N",'Y']
default['RZ'] = ["Compute RZ correlation [Y]/N",'Y']
default['RR'] = ["Compute RR correlation [Y]/N",'Y']
default['RT'] = ["Compute RT correlation [Y]/N",'Y']
default['TZ'] = ["Compute TZ correlation [Y]/N",'Y']
default['TR'] = ["Compute TR correlation [Y]/N",'Y']
default['TT'] = ["Compute TT correlation [Y]/N",'Y']

default['autocorr'] = ["Compute Auto correlation [Y]/N",'N']
default['PAZ'] = ["Correct instrumental responce from paz [Y]/N",'N']
default['keep_all'] = ["Keep all 30 seconds cross-corr [Y]/N",'N']
default['keep_days'] = ["Keep all daily cross-corr [Y]/N",'Y']

default['ref_begin'] = ["Beginning or REF stacks. Can be absolute (2012-01-01) or relative (-100) days",'-100']
default['ref_end'] = ["End or REF stacks. Same as ref_begin",'0']

default['mov_stack'] = ["Number of days to stack for the Moving-window stacks ([5]= [day-4:day]), can be a comma-separated list 2,5,10","5"]

default['export_format'] = ["Export stacks in which format(s) ? SAC/MSEED/[BOTH]","MSEED"]
default['sac_format'] = ["Format for SAC stacks ? [doublets]/clarke","doublets"]



# MsNoise classes:
# Raw Stations & Filters

class Station:
    ref = None
    net = None
    sta = None
    name = None
    X = None
    Y = None
    altitude = None
    coordinates = None
    data_folder = None
    file_format = None
    params = None
    format = None
    
    def __init__(self,*args,**kwargs):
        for kwarg in kwargs:
            if kwarg == 'ref':
                self.ref = kwargs[kwarg]
            elif kwarg == 'net':
                self.net = kwargs[kwarg]
            elif kwarg == 'sta':
                self.sta = kwargs[kwarg]
            elif kwarg == 'X':
                self.X = kwargs[kwarg]
            elif kwarg == 'Y':
                self.Y = kwargs[kwarg]
            elif kwarg == 'coordinates':
                self.coordinates = kwargs[kwarg]
            elif kwarg == 'altitude':
                self.altitude = kwargs[kwarg]
            elif kwarg == 'data_folder':
                self.data_folder = kwargs[kwarg]
            elif kwarg == 'file_format':
                self.file_format = kwargs[kwarg]
            elif kwarg == 'params':
                self.params = kwargs[kwarg]
            elif kwarg == 'format':
                self.format = kwargs[kwarg]
            

class Filter:
    ref = None
    low = None
    high = None
    mwcs_low = None
    mwcs_high = None
    rms_threshold = None
    mwcs_wlen = None
    mwcs_step = None
    used = None
    
    def __init__(self,*args,**kwargs):
        for kwarg in kwargs:
            if kwarg == 'ref':
                self.ref = kwargs[kwarg]
            elif kwarg == 'low':
                self.low = float(kwargs[kwarg])
            elif kwarg == 'mwcs_low':
                self.mwcs_low = float(kwargs[kwarg])
            elif kwarg == 'high':
                self.high = float(kwargs[kwarg])
            elif kwarg == 'mwcs_high':
                self.mwcs_high = float(kwargs[kwarg])
            elif kwarg == 'mwcs_wlen':
                self.mwcs_wlen = float(kwargs[kwarg])
            elif kwarg == 'mwcs_step':
                self.mwcs_step = float(kwargs[kwarg])
            elif kwarg == 'rms_threshold':
                self.rms_threshold = kwargs[kwarg]
            elif kwarg == 'used':
                self.used = kwargs[kwarg]

def connect():
    f = open('db.ini','r')
    hostname,database,username,password = cPickle.load(f)
    f.close()
    db = SQL.connect(host=hostname,db=database,user=username,passwd=password)
    return db

def disconnect(db):
    if 1:
        db.close()
        del db
    # except:
        # print sys.exc_info()
    
def update_station(db, net, sta, X, Y, altitude,coordinates='UTM',instrument='N/A',used=1,):
    c = db.cursor()
    c.execute('SELECT ref from _stations where net="%s" and sta="%s"'%(net,sta))
    r = c.fetchall()
    if len(r) != 0:
        ref = r[0][0]
        request = 'update _stations set net ="%s" , sta="%s", X=%f, Y=%f, altitude=%f, coordinates="%s", instrument="%s", used="%i" where ref=%i'%(net, sta, X, Y, altitude, coordinates, instrument,used,ref)
    else:
        request = 'insert into _stations (net, sta, X, Y, altitude, coordinates, instrument) values ("%s","%s",%f,%f,%f,"%s","%s")'%(net, sta, X, Y, altitude, coordinates, instrument)
    c.execute(request)
    db.commit()
    
def get_stations(db,used=None,net=None):
    c = db.cursor()
    stations = []
    if not used:
        if net != None:
            c.execute("select ref, net, sta, X, Y, altitude, coordinates, instrument, used from _stations where net = '%s' order by net,sta"%net)
        else:
            c.execute("select ref, net, sta, X, Y, altitude, coordinates, instrument, used from _stations order by net,sta")
        for istation in c.fetchall():
            s = Station()
            s.ref, s.net, s.sta, s.X, s.Y, s.altitude, s.coordinates, s.instrument, s.used = istation
            stations.append(s)
    else:
        c.execute("select ref, net, sta, X, Y, altitude, coordinates, instrument, used from _stations where used=1 order by net,sta")
        for istation in c.fetchall():
            s = Station()
            s.ref, s.net, s.sta, s.X, s.Y, s.altitude, s.coordinates, s.instrument,s.used = istation
            stations.append(s)

    return stations

def get_station_pairs(db,used=None,net=None):
    if get_config(db, name="autocorr") in ['Y','y','1',1]:
        return itertools.combinations_with_replacement(get_stations(db,used=used,net=net),2)
    else:
        return itertools.combinations(get_stations(db,used=used,net=net),2)


def get_networks(db,used=None):
    c = db.cursor()
    if not used:
        c.execute("SELECT net FROM _stations GROUP BY net")
    else:
        c.execute("SELECT net FROM _stations where used=1 GROUP BY net")
    return [net[0] for net in c.fetchall()]

def get_station(db,net,sta):
    c = db.cursor()
    c.execute("select ref, net, sta, X, Y, altitude, coordinates, instrument from _stations where net='%s' and sta='%s'"%(net,sta))
    stations = []
    for istation in c.fetchall():
        s = Station()
        s.ref, s.net, s.sta, s.X, s.Y, s.altitude, s.coordinates, s.instrument = istation
        stations.append(s)
    return stations[0]

def get_filters(db,all=True):
    c = db.cursor()
    if all:
        c.execute("select ref, low, high, mwcs_low, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used from _filters")
    else:
        c.execute("select ref, low, high, mwcs_low, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used from _filters where used=1")
    filters = []
    for ifilter in c.fetchall():
        ref, low, high, mwcs_low, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used = ifilter
        f = Filter(ref=ref, low=low, high=high,mwcs_low=mwcs_low, mwcs_high=mwcs_high, rms_threshold=rms_threshold, mwcs_wlen=mwcs_wlen, mwcs_step=mwcs_step,used=used)
        filters.append(f)
    return filters

def update_filter(db,ref, low, high, mwcs_low, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used):
    c = db.cursor()
    c.execute('SELECT ref from _filters where ref = %i'%(ref))
    r = c.fetchall()
    if len(r) != 0:
        ref = r[0][0]
        request = 'update _filters set low = "%f", high="%f", mwcs_low="%f", mwcs_high="%f", rms_threshold="%f",  mwcs_wlen="%f", mwcs_step="%f", used="%i"  where ref=%i'%(low,high,mwcs_low, mwcs_high,rms_threshold, mwcs_wlen, mwcs_step,used,ref)
    else:
        request = 'insert into _filters (low,high,mwcs_low, mwcs_high,rms_threshold, mwcs_wlen, mwcs_step,used) values (%f,%f,%f,%f,%f,%f,%f,%i)'%(low,high,mwcs_low, mwcs_high,rms_threshold, mwcs_wlen, mwcs_step,used)
    c.execute(request)
    db.commit()

def get_filter_bounds(db,filterid):
    c = db.cursor()
    c.execute("select ref, low, high from _filters where ref = %i"%filterid)
    filter = c.fetchall()[0]
    ref, low, high = filter

    return (float(low), float(high))

def create_database_inifile(hostname,database,username,password):
    f = open('db.ini','w')
    cPickle.dump([hostname,database,username,password],f)
    f.close()

def create_database(root,passwd):
    f = open('db.ini','r')
    hostname,database,username,password = cPickle.load(f)
    f.close()
    db = SQL.connect(host=hostname,user=root,passwd=passwd)
    c = db.cursor()
    c.execute("CREATE DATABASE `%s` DEFAULT CHARACTER SET latin1 COLLATE latin1_swedish_ci;"%database)
    c.execute("GRANT ALL PRIVILEGES ON `%s` . * TO '%s'@'%%'"%(database,username))
    db.close()

def create_stations_table(db):
    c = db.cursor()
    c.execute(""" CREATE TABLE IF NOT EXISTS `_stations` (
    `ref` int(11) NOT NULL AUTO_INCREMENT,
    `net` varchar(10) NOT NULL,
    `sta` varchar(10) NOT NULL,
    `X` float NOT NULL,
    `Y` float NOT NULL,
    `altitude` float NOT NULL,
    `coordinates` enum('DEG','UTM') NOT NULL,
    `instrument` varchar(10) NOT NULL,
    `used` BOOL NOT NULL DEFAULT 1,
    PRIMARY KEY (`ref`)
    ) ENGINE=InnoDB  DEFAULT CHARSET=latin1;""")
    
def create_filters_table(db):
    c = db.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS `_filters` (
    `ref` int(11) NOT NULL AUTO_INCREMENT,
    `low` varchar(10) NOT NULL,
    `mwcs_low` varchar(10) NOT NULL,
    `high` varchar(10) NOT NULL,
    `mwcs_high` varchar(10) NOT NULL,
    `rms_threshold` float NOT NULL,
    `mwcs_wlen` float NOT NULL,
    `mwcs_step` float NOT NULL,
    `used` BOOL NOT NULL DEFAULT 1,
    PRIMARY KEY (`ref`)
    ) ENGINE=InnoDB  DEFAULT CHARSET=latin1;""")

def create_config_table(db):
    c = db.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS `_config` (
    `name` varchar(255) NOT NULL,
    `value` varchar(255) NOT NULL,
    PRIMARY KEY (`name`)
    ) ENGINE=InnoDB  DEFAULT CHARSET=latin1;""")

    for name in default.keys():
        c.execute('insert into _config (name,value) values ("%s","%s")'%(name,default[name][1]))
    db.commit()
    
def create_data_availability_table(db):
    c = db.cursor()
    c.execute("""CREATE TABLE IF NOT EXISTS `_data_availability` (
    `ref` int(11) NOT NULL AUTO_INCREMENT,
    `net` varchar(20) NOT NULL,
    `sta` varchar(20) NOT NULL,
    `comp` varchar(20) NOT NULL,
    `path` varchar(255) NOT NULL,
    `file` varchar(255) NOT NULL,
    `starttime` datetime NOT NULL,
    `endtime` datetime NOT NULL,
    `data_duration` float(9,3) NOT NULL,
    `gaps_duration` float(9,3) NOT NULL,
    `samplerate` float(9,3) NOT NULL,
    `flag` varchar(1) NOT NULL DEFAULT 'N',
    PRIMARY KEY (`ref`)
    ) ENGINE=InnoDB  DEFAULT CHARSET=latin1;""")

def create_data_availability_views(db):
    c = db.cursor()
    c.execute("""
    CREATE view _data_availability_A as SELECT starttime, endtime, data_duration, gaps_duration, GROUP_CONCAT(net,'.',sta) as A 
        FROM `_data_availability` where flag='A' group by starttime,endtime, data_duration , gaps_duration
        order by starttime, net, sta asc;
        
    CREATE view _data_availability_N as SELECT starttime, endtime, data_duration, gaps_duration, GROUP_CONCAT(net,'.',sta) as N 
        FROM `_data_availability` where (flag='N' or flag='M') group by starttime,endtime, data_duration , gaps_duration
        order by starttime, net, sta asc;
    """)
  
def add_data_availability(db,net, sta,comp, path,file,starttime, endtime, data_duration, gaps_duration, samplerate):
    starttime = starttime.isoformat().replace('T',' ')
    endtime = endtime.isoformat().replace('T',' ')
    c = db.cursor()
    request = 'select ref,net, sta,comp, path,file ,starttime, endtime, data_duration, gaps_duration, samplerate from _data_availability where net = "%s" and sta = "%s" and comp ="%s" and path = "%s" and file ="%s" and starttime ="%s" and endtime = "%s" and data_duration = %.3f and gaps_duration = %.3f and samplerate = %.3f '%(net, sta, comp, path,file,starttime, endtime,data_duration,gaps_duration,samplerate)
    c.execute(request)
    r = c.fetchall()
    if len(r) == 1:
        # the very same file is in the DB, not doing anything
        return 1
    if len(r) == 0:
        # the record is either absent or partially the same
        request = 'select ref,net, sta,comp, path,file ,starttime, endtime, data_duration, gaps_duration, samplerate from _data_availability where net = "%s" and sta = "%s" and comp ="%s" and path = "%s" and file ="%s" '%(net, sta, comp, path,file)
        c.execute(request)
        r = c.fetchall()
        if len(r) == 1:
            # data is partially there, updating the content of the row
            ref = r[0][0]
            request = "update _data_availability set starttime='%s', endtime='%s', data_duration='%s', gaps_duration='%s', samplerate='%s',flag='M' where ref=%i"%(starttime,endtime,data_duration,gaps_duration,samplerate,ref)
        else:
            # data is absent, inserting
            request = 'insert into _data_availability (net, sta,comp, path,file,starttime, endtime,data_duration,gaps_duration,samplerate) values ("%s","%s","%s","%s","%s","%s","%s",%.3f,%.3f,%.3f)'%(net, sta, comp, path,file,starttime, endtime,data_duration,gaps_duration,samplerate)
        c.execute(request)
        db.commit()
        return 1
    else:
        return 0

def update_data_availability(db,net, sta,comp, path,file,starttime, endtime, data_duration, gaps_duration, samplerate):
    c = db.cursor()
    c.execute('select ref from _data_availability where file="%s" order by starttime asc'%(file))
    result = c.fetchall()
    if len(result) !=0:
        ref = result[0][0]
        for item in ["net", "sta","comp", "path","file","starttime","endtime", "data_duration","gaps_duration", "samplerate"]:
            c.execute('update _data_availability set %s="%s" where ref=%i'%(item, eval(item),ref))
        c.execute('update _data_availability set flag="M" where ref=%i'%ref)
        db.commit()
        return 0
    else:
        request = 'insert into _data_availability (net, sta,comp, path,file,starttime,endtime,data_duration,gaps_duration,samplerate) values ("%s","%s","%s","%s","%s","%s","%s",%f,%f,%f)'%(net, sta, comp, path,file,starttime,endtime,data_duration,gaps_duration,samplerate)
        c.execute(request)
        db.commit()
        return 1
        

def mark_data_availability(db,condition='1', flag='A'):
    c = db.cursor()
    request = 'update _data_availability set flag="%s" where %s'%(flag, condition)
    c.execute(request)
    db.commit()

def get_data_availability(db, net=None, sta=None, comp=None, starttime=None,endtime=None):
    c = db.cursor()
    if not starttime:
        request = 'select starttime, endtime, data_duration, gaps_duration, samplerate, flag from _data_availability where net = "%s" and sta = "%s" and comp = "%s" order by starttime asc'%(net,sta,comp)
    else:
        request = 'select net, sta, comp, starttime, endtime, data_duration, gaps_duration, samplerate, flag from _data_availability where date(starttime) <="%s" and date(endtime) >= "%s"'%(starttime,endtime)
    c.execute(request)
    data = c.fetchall()
    return data

def pop_data_availability(db):
    c = db.cursor()
    request = 'select net, sta, comp from _data_availability GROUP BY net, sta, comp'
    c.execute(request)
    data = c.fetchall()
    return data

def get_new_files(db):
    c = db.cursor()
    c.execute("select net, sta, comp, starttime, endtime, data_duration from _data_availability where flag='N' or flag='M' order by date(starttime), net, sta asc")
    results = c.fetchall()
    return results

def get_config(db,name=None):
    c = db.cursor()
    if name:
        c.execute('select value from _config where name="%s"'%name)
        value = c.fetchall()
        if len(value) != 0:
            return value[0][0]
        else:
            return ''
    else:
        c.execute('select name, value from _config')
        pairs = c.fetchall()
        config = {}
        for pair in pairs:
            name,value=pair
            config[name]=value
        return config

def update_config(db,name, value):
    c = db.cursor()
    c.execute('select value from _config where name = "%s"'%name)
    if len(c.fetchall()) !=0:
        c.execute('update _config set value="%s" where name="%s"'%(value,name))
    else:
        c.execute('insert into _config (name,value) values ("%s","%s")'%(name,value))
    db.commit()

def get_maxlag_samples(db):
    maxlag = float(get_config(db,'maxlag'))
    return int(2*maxlag*float(get_config(db,'cc_sampling_rate'))+1)

def get_results(db,station1, station2, filterid, components, dates,format="stack"):
    export_format = get_config(db,'export_format')
    if format=="stack":
        stack = np.zeros(get_maxlag_samples(db))
        i = 0
        for date in dates:
            daystack = os.path.join("STACKS","%02i"%filterid,"001_DAYS",components,"%s_%s"%(station1, station2),str(date))
            if export_format == "BOTH":
                daystack += ".MSEED"
            elif export_format == "SAC":
                daystack += ".SAC"
            elif export_format == "MSEED":
                daystack += ".MSEED"
            if os.path.isfile(daystack):
                st = read(daystack)
                if not np.any(np.isnan(st[0].data)) and not np.any(np.isinf(st[0].data)):
                    stack += st[0].data
                    i += 1
                else:
                    # print "NaN ! or Inf"
                    pass
        if i > 0:
            return i, stack / i
        else: 
            return 0, None
    
    elif format=="matrix":
        stack = np.zeros((len(dates),get_maxlag_samples(db)))
        i = 0
        for j, date in enumerate(dates):
            daystack = os.path.join("STACKS","%02i"%filterid,"001_DAYS",components,"%s_%s"%(station1, station2),str(date))
            if export_format == "BOTH":
                daystack += ".MSEED"
            elif export_format == "SAC":
                daystack += ".SAC"
            elif export_format == "MSEED":
                daystack += ".MSEED"
            if os.path.isfile(daystack):
                st = read(daystack)
                stack[j] = st[0].data
                i += 1
            else:
                stack[j] *= np.nan
        return i, stack



def get_cf_availability(db, station1, station2, filterid, components, year, month, day):
    #SHOULD BE REPLACED WITH A FOLDER-SCANNING METHOD
    c = db.cursor()
    request = "select count from corr_%s_%s_%i_status where components=\"%s\" and YEAR(date) = %i and MONTH(date)=%i and DAYOFMONTH(date)=%i"%(station1, station2, filterid,components,year, month, day)
    c.execute(request)
    data = c.fetchall()
    if len(data) != 0:
        count = data[0][0]
    else:
        count = 0
    return count

def get_filenames(db, day, net, sta):
    c = db.cursor()
    request = 'select net,sta,comp,path,file,starttime, endtime,data_duration,samplerate from _data_availability where date(starttime)<="%s" and date(endtime) >= "%s" and net="%s" and sta="%s"'%(day,day,net,sta)
    c.execute(request)
    r = c.fetchall()
    return np.array(r)

def create_dtt_table(db):
    c = db.cursor()
    request = """CREATE TABLE IF NOT EXISTS `_dtts` (
    `ref` int(11) NOT NULL AUTO_INCREMENT,
    `day` varchar(10) NOT NULL,
    `pair` varchar(20) NOT NULL,
    `filter` int(2) NOT NULL,
    `components` varchar(2) NOT NULL,
    `mov_stack` int(3) NOT NULL,
    `method` varchar(20) NOT NULL,
    `M` double NOT NULL,
    `EM` double NOT NULL,
    `A` double NOT NULL,
    `EA` double NOT NULL,
    `M0` double NOT NULL,
    `EM0` double NOT NULL,
    `lastmod` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    PRIMARY KEY (`ref`)
    ) ENGINE=InnoDB  DEFAULT CHARSET=latin1;"""
    c.execute(request)


def add_dtt(db,day,pair,filterid,components,mov_stack,method,M,EM,A,EA,M0,EM0):
    c = db.cursor()
    request = "select 1 from _dtts where day='%s' and pair='%s' and filter='%i' and components='%s' and mov_stack=%i and method='%s' "%(day,pair,filterid,components,mov_stack,method)
    c.execute(request)
    if len(c.fetchall()) == 0:
        request = "insert into _dtts (day, pair,filter,components,mov_stack,method,M,EM,A,EA,M0,EM0) values ('%s','%s',%i,'%s',%i,'%s',%f,%f,%f,%f,%f,%f)"%(day,pair,filterid,components,mov_stack,method,M,EM,A,EA,M0,EM0)
        c.execute(request)
        db.commit()
    else:
        update_dtt(db,day,pair,filterid,components,mov_stack,method,M,EM,A,EA,M0,EM0)

def update_dtt(db,day,pair,filterid,components,mov_stack,method,M,EM,A,EA,M0,EM0):
    c = db.cursor()
    request = "update _dtts set M=%f, EM=%f,A=%f, EA=%f,M0=%f, EM0=%f where day = '%s' and pair ='%s' and filter =%i and components='%s' and mov_stack=%i and method='%s'"%(M,EM,A,EA,M0,EM0,day,pair,filterid,components,mov_stack,method)
    c.execute(request)
    db.commit()

def create_jobs_table(db):
    c = db.cursor()
    request = """CREATE TABLE IF NOT EXISTS `_jobs` (
    `ref` int(11) NOT NULL AUTO_INCREMENT,
    `day` varchar(10) NOT NULL,
    `pair` varchar(20) NOT NULL,
    `type` varchar(10) NOT NULL,
    `flag` varchar(1) NOT NULL DEFAULT 'T',
    `lastmod` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    PRIMARY KEY (`ref`)
    ) ENGINE=InnoDB  DEFAULT CHARSET=latin1;"""
    c.execute(request)

def add_job(db,day,pair,type):
    c = db.cursor()
    request = "select 1 from _jobs where day='%s' and pair='%s' and type='%s'"%(day,pair,type)
    c.execute(request)
    if len(c.fetchall()) == 0:
        request = "insert into _jobs (day, pair,type) values ('%s','%s','%s')"%(day,pair,type)
        c.execute(request)
        db.commit()
    else:
        update_job(db,day,pair,type,flag='T')

def is_next_job(db, flag = 'T',type='CC'):
    c = db.cursor()
    request = 'select day, GROUP_CONCAT(pair), GROUP_CONCAT(ref) from _jobs where flag = "%s" and type="%s" group by day,flag limit 1' % (flag,type)
    c.execute(request)
    r = c.fetchall()
    if len(r) != 0:
        return True
    else:
        return False

def get_next_job(db, flag = 'T',type='CC'):
    c = db.cursor()
    request = 'select day, GROUP_CONCAT(pair), GROUP_CONCAT(ref) from _jobs where flag = "%s" and type="%s" group by day,flag limit 1' % (flag,type)
    c.execute(request)
    r = c.fetchall()[0]
    if len(r) != 0:
        return np.array(r)
    else:
        return False

def update_job(db,day,pair,type='CC',flag='I'):
    c = db.cursor()
    request = "update _jobs set flag='%s' where day = '%s' and pair ='%s' and type ='%s' "%(flag,day,pair,type)
    c.execute(request)
    c.close()
    db.commit()

def get_job_status(db):
    c = db.cursor()
    request = "SELECT flag, COUNT(*) FROM _jobs GROUP BY flag"
    c.execute(request)
    r = c.fetchall()
    return r

def get_all_jobs_times(db):
    c = db.cursor()
    request = "SELECT flag, lastmod FROM _jobs"
    c.execute(request)
    r = c.fetchall()
    return r

def is_dtt_ref_job(db, pair, flag = 'T',type='DTT'):
    c = db.cursor()
    request = 'select 1 from _jobs where pair="%s" and flag = "%s" and type="%s" and day="REF" ' % (pair, flag,type)
    c.execute(request)
    r = c.fetchall()
    c.close()
    if len(r) != 0:
        return True
    else:
        return False

def is_dtt_mov_job(db, flag = 'T',type='DTT'):
    c = db.cursor()
    request = 'select pair, GROUP_CONCAT(day) from _jobs where flag = "%s" and type="%s" and day!="REF" group by pair limit 1' % (flag,type)
    c.execute(request)
    r = c.fetchall()
    c.close()
    if len(r) != 0:
        return True
    else:
        return False

def update_dtt_job(db,refs,type='CC',flag='I'):
    c = db.cursor()
    request = "update _jobs set flag='%s' where ref in (%s) "%(flag,refs)
    c.execute(request)
    db.commit()
    c.close()
    

def get_dtt_next_job(db, flag = 'T',type='DTT'):
    c = db.cursor()
    request = 'select pair, GROUP_CONCAT(day), GROUP_CONCAT(ref) from _jobs where flag = "%s" and type="%s" and day!="REF" group by pair order by day asc limit 1 for update' % (flag,type)
    c.execute(request)
    r = c.fetchall()
    
    if len(r) != 0:
        pair,days,refs=r[0]
        update_dtt_job(db,refs,type='DTT',flag='I')
        c.close()
        return np.array(r[0])
    else:
        return False, False

def reset_dtt_jobs(db, pair):
    c = db.cursor()
    request = 'update _jobs set flag = "T" where type="DTT" and pair="%s"' % (pair)
    c.execute(request)
    db.commit()
    c.close()

def add_corr(db,station1, station2, filterid, date, time, duration, components, CF, sampling_rate,day=False,ncorr=0):
    output_folder = get_config(db, 'output_folder')
    export_format = get_config(db,'export_format')
    if export_format == "BOTH":
        mseed = True
        sac = True
    elif export_format == "SAC":
        mseed = False
        sac = True
    elif export_format == "MSEED":
        mseed = True
        sac = False
    
    if day:
        path = os.path.join("STACKS","%02i"%filterid,"001_DAYS",components,"%s_%s"%(station1,station2),str(date))
        pair = "%s:%s"%(station1,station2)
        if mseed:
            export_mseed(db, path, pair, components, filterid,CF/ncorr,ncorr)
        if sac:
            export_sac(db, path, pair, components, filterid,CF/ncorr,ncorr)
    
    else:
        file = '%s.cc' % time
        path = os.path.join(output_folder, "%02i"% filterid, station1, station2,  components,  date)
        if not os.path.isdir(path):
            os.makedirs(path)
        
        t = Trace()
        t.data = CF
        t.stats.sampling_rate = sampling_rate
        t.stats.starttime=-float(get_config(db,'maxlag'))
        t.stats.components = components
        # if ncorr != 0:
            # t.stats.location = "%02i"%ncorr
        st = Stream(traces= [t,])
        st.write(os.path.join(path, file),format='mseed')
        del t, st

def allow_large_concats(db):
    c=db.cursor()
    request = "SET SESSION group_concat_max_len = 1000000000;"
    c.execute(request)
    

def build_ref_datelist(db):
    begin = get_config(db, "ref_begin")
    end = get_config(db, "ref_end")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin,'%Y-%m-%d').date()
        end = datetime.datetime.strptime(end,'%Y-%m-%d').date()
    
    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]

def build_movstack_datelist(db):
    begin = get_config(db, "startdate")
    end = get_config(db, "enddate")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin,'%Y-%m-%d').date()
        end = datetime.datetime.strptime(end,'%Y-%m-%d').date()
    
    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]

def build_daystack_datelist(db):
    refstart, refend, refdates = build_ref_datelist(db)
    movstart, movend, movdates = build_movstack_datelist(db)
    start = min(refstart,movstart)
    end = max(refend, movend)
    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]
    

def get_components_to_compute(db):
    components_to_compute = []
    for comp in ['ZZ','RR','TT','TR','RT','ZR','RZ','TZ','ZT']:
        if get_config(db, comp) in ['Y','y','1',1]:
            components_to_compute.append(comp)
    return components_to_compute

def updated_days_for_dates(db, date1, date2, pair,type='CC',interval='1 DAY',returndays=False):
    c = db.cursor()
    if pair != "%":
        pair = "AND pair='%s'" % pair
    else:
        pair = ""

    request = "SELECT ref,pair,day,pair FROM _jobs where day BETWEEN '%s' and '%s' %s AND type='%s' AND lastmod >= DATE_SUB(NOW(),INTERVAL %s)" %(date1, date2, pair, type, interval)
    c.execute(request)
    result = c.fetchall()
    if len(result) > 0:
        if not returndays:
            return True
        else:
            request = "SELECT DATE(day) FROM _jobs where day BETWEEN '%s' and '%s' %s AND type='%s'  AND lastmod >= DATE_SUB(NOW(),INTERVAL %s) group by day" %(date1, date2, pair, type, interval)
            c.execute(request)
            return [d[0] for d in c.fetchall()]
    else:
        if not returndays:
            return False
        else:
            return []


def export_sac(db, filename, pair, components, filterid,corr,ncorr=0,sac_format=None,maxlag=None,cc_sampling_rate=None):
    if sac_format is None:
        sac_format=get_config(db,"sac_format")
    if maxlag is None:
        maxlag = float(get_config(db,"maxlag"))
    if cc_sampling_rate is None:
        cc_sampling_rate = float(get_config(db,"cc_sampling_rate"))
    try:
        os.makedirs(os.path.split(filename)[0])
    except:
        pass
    filename += ".SAC"
    mytrace = Trace(data=corr)
    mytrace.stats['station'] = pair
    mytrace.stats['sampling_rate'] = cc_sampling_rate


    st = Stream(traces = [mytrace,])            
    st.write(filename,format='SAC')
    tr = SacIO(filename)
    if sac_format == "doublets":
        tr.SetHvalue('A',120)
    else:
        tr.SetHvalue('B',-maxlag)
        tr.SetHvalue('DEPMIN',np.min(corr))
        tr.SetHvalue('DEPMAX',np.max(corr))
        tr.SetHvalue('DEPMEN',np.mean(corr))
        tr.SetHvalue('SCALE',1)
        tr.SetHvalue('NPTS',len(corr))
    tr.WriteSacBinary(filename)
    del st, tr
    return

def export_mseed(db, filename, pair, components, filterid,corr,ncorr=0):
    try:
        os.makedirs(os.path.split(filename)[0])
    except:
        pass
    filename += ".MSEED"
    maxlag = float(get_config(db,"maxlag"))
    cc_sampling_rate = float(get_config(db,"cc_sampling_rate"))
    
    mytrace = Trace(data=corr)
    mytrace.stats['station'] = pair[:11]
    mytrace.stats['sampling_rate'] = cc_sampling_rate
    mytrace.stats['start_time'] = -maxlag
    mytrace.stats['location'] = "%02i" % ncorr

    st = Stream(traces = [mytrace,])            
    st.write(filename,format='MSEED')
    del st
    return



def oldazimuth(coordinates, x0, y0, x1, y1):
    if coordinates == 'DEG':
        lon1 = x0
        lat1 = y0
        
        lon2 = x1
        lat2 = y1

        R = 6371.0
        dLat = (lat2-lat1) * np.pi / 180.
        dLon = (lon2-lon1)* np.pi / 180.
        a = np.sin(dLat/2) * np.sin(dLat/2) + np.cos(lat1 * np.pi/180.) * np.cos(lat2 * np.pi/180.) * np.sin(dLon/2) * np.sin(dLon/2); 
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        d = R * c

        y = np.sin(dLon) * np.cos(lat2);
        x = np.cos(lat1* np.pi/180.)*np.sin(lat2* np.pi/180.) - np.sin(lat1* np.pi/180.)*np.cos(lat2* np.pi/180.)*np.cos(dLon);
        brng = np.arctan2(y, x) * 180. / np.pi
        return brng % 360.
    
    elif coordinates == 'UTM':
        azim = 90. - np.arctan2((y1-y0),(x1-x0)) *180./np.pi
        # print ">>> AZIMUTH:", azim
        
        return azim
    
    else:
        print 'mixed'
        
    return 0

def azimuth(coordinates, x0, y0, x1, y1):
    if coordinates == "DEG":
        dist, azim, bazim = gps2DistAzimuth(y0,x0,y1,x1)
        # print dist, azim, bazi
        return azim
    elif coordinates == 'UTM':
        azim = 90. - np.arctan2((y1-y0),(x1-x0)) *180./np.pi
        # print azim
        return azim
    else:
        print "woooooow, please consider having a single coordinate system for all stations"
        return 0



#EOF
