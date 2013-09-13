# queries.py
from sqlalchemy import create_engine, Date, func
from sqlalchemy.orm import sessionmaker
import numpy as np
import datetime
import itertools

from obspy.core import Stream, Trace, read
from obspy.sac import SacIO
import os

from msnoise_table_def import *

def connect():
    engine = create_engine('sqlite:///msnoise.sqlite', echo=False)
    Session = sessionmaker(bind=engine)
    return Session()


############ CONFIG ############
def get_config(session,name=None):
    if name:
        config = session.query(Config).filter(Config.name==name).first()
        if config is not None:
            config = config.value
        else:
            config = ''
    else:
        config = {}
        configs = session.query(Config).all()
        for c in configs:
            config[c.name] = c.value
    return config

def update_config(session,name, value):
    config = session.query(Config).filter(Config.name==name).first()
    config.value = value
    session.commit()

############ FILTERS ############
def get_filters(session, all=False):
    if all:
        filters = session.query(Filter).all()
    else:
        filters = session.query(Filter).filter(Filter.used==True).all()
    return filters

def update_filter(session,ref, low, high, mwcs_low, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used):
    filter = session.query(Filter).filter(Filter.ref==ref).first()
    if filter is None:
        filter = Filter(low, mwcs_low, high, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used)
        session.add(filter)
    else:
        filter.low = low
        filter.high = high
        filter.mwcs_low = mwcs_low
        filter.mwcs_high = mwcs_high
        filter.rms_threshold = rms_threshold
        filter.mwcs_wlen = mwcs_wlen
        filter.mwcs_step = mwcs_step
        filter.used = used
    session.commit()

############ NETWORK AND STATION ############
def get_networks(session, all=False):
    if all:
        networks = session.query(Station).group_by(Station.net).all()
    else:
        networks = session.query(Station).filter(Station.used ==True).group_by(Station.net)
    return [net.net for net in networks]

def get_stations(session, all=False,net=None):
    if all:
        if net != None:
            stations = session.query(Station).filter(Station.net==net)
        else:
            stations = session.query(Station).all()
    else:
        stations = session.query(Station).filter(Station.used==True)
        if net != None:
            stations = stations.filter(Station.net==net)
    return stations

def get_station(session, net, sta):
    station = session.query(Station).filter(Station.net==net).filter(Station.sta==sta).first()
    return station

def update_station(session, net, sta, X, Y, altitude,coordinates='UTM',instrument='N/A',used=1):
    station = session.query(Station).filter(Station.net==net).filter(Station.sta==sta).first()
    if station is None:
        station=Station(net, sta, X, Y, altitude,coordinates,instrument,used)
        session.add(station)
    else:
        station.X = X
        station.Y = Y
        station.altitude = altitude
        station.instrument = instrument
        station.used = used
    session.commit()

def get_station_pairs(session,used=None,net=None):
    if get_config(session, name="autocorr") in ['Y','y','1',1]:
        return itertools.combinations_with_replacement(get_stations(session,all=False,net=net),2)
    else:
        return itertools.combinations(get_stations(session,all=False,net=net),2)

############ DATA AVAILABILITY ############

def update_data_availability(session, net, sta, comp, path, file, starttime, endtime, data_duration, gaps_duration, samplerate):
    data = session.query(DataAvailability).filter(DataAvailability.file==file).first()
    if data is None:
        flag = "N"
        data = DataAvailability(net,sta,comp,path,file,starttime,endtime,data_duration,gaps_duration,samplerate,flag)
        session.add(data)
        toreturn = True
    else:
        modified = False
        for item in ['net','sta','comp','path','starttime','data_duration','gaps_duration','samplerate']:
            if eval("data.%s != %s"%(item,item)):
                modified = True
                break
        if modified:
            data.net = net
            data.sta = sta
            data.comp = comp
            data.path = path
            data.starttime = starttime
            data.data_duration = data_duration
            data.gaps_duration = gaps_duration
            data.samplerate = samplerate
            data.flag="M"
        
        toreturn = False
    session.commit()
    return toreturn

def get_new_files(session):
    files = session.query(DataAvailability).filter(DataAvailability.flag != 'A')
    return files

def get_data_availability(session, net=None, sta=None, comp=None, starttime=None,endtime=None):
    if not starttime:
        data = session.query(DataAvailability).filter(DataAvailability.net==net).filter(DataAvailability.sta==sta).filter(DataAvailability.comp==comp).all()
    if not net:
        data = session.query(DataAvailability).filter(DataAvailability.starttime <= starttime).filter(DataAvailability.endtime>=endtime).all()
    else:
        data = session.query(DataAvailability).filter(DataAvailability.net==net).filter(DataAvailability.sta==sta).filter(func.DATE(DataAvailability.starttime) <= starttime.date()).filter(func.DATE(DataAvailability.endtime)>=endtime.date()).all()
    return data

def mark_data_availability(session,net,sta,flag):
    data = session.query(DataAvailability).filter(DataAvailability.net == net).filter(DataAvailability.sta == sta).all()
    for d in data:
        d.flag=flag
    session.commit()

############ JOBS ############

def update_job(session,day,pair,type,flag,commit=True,returnjob=True):
    job = session.query(Job).filter(Job.day==day).filter(Job.pair==pair).filter(Job.type==type).first()
    if job is None:
        job = Job(day,pair,type,'T')
        if commit:
            session.add(job)
    else:
        job.flag = flag
        job.lastmod = datetime.datetime.utcnow()
    if commit:
        session.commit()
    if returnjob:
        return job

def is_next_job(session, flag = 'T',type='CC'):
    job = session.query(Job).filter(Job.type==type).filter(Job.flag==flag).first()
    if job is None:
        return False
    else:
        return True

def get_next_job(session, flag = 'T',type='CC'):
    day = session.query(Job).filter(Job.type==type).filter(Job.flag==flag).first().day
    # print day
    jobs = session.query(Job).filter(Job.type==type).filter(Job.flag==flag).filter(Job.day==day).all()
    return jobs


def is_dtt_next_job(session, flag = 'T',type='DTT',ref=False):
    if ref:
        job = session.query(Job).filter(Job.flag==flag).filter(Job.type==type).filter(Job.day=='REF').first()
    else:
        job = session.query(Job).filter(Job.flag==flag).filter(Job.type==type).filter(Job.day!='REF').first()
    if job is None:
        return False
    else:
        return True

def get_dtt_next_job(session, flag = 'T',type='DTT'):
    pair = session.query(Job).filter(Job.flag==flag).filter(Job.type==type).filter(Job.day!='REF').first().pair
    jobs = session.query(Job).filter(Job.flag==flag).filter(Job.type==type).filter(Job.day!='REF').filter(Job.pair==pair).all()
    refs = [job.ref for job in jobs]
    days = [job.day for job in jobs]
    for job in jobs:
        job.flag='I'
    session.commit()
    return pair, days, refs
    


############ CORRELATIONS ############

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

def get_results(session,station1, station2, filterid, components, dates,format="stack"):
    export_format = get_config(session,'export_format')
    if format=="stack":
        stack = np.zeros(get_maxlag_samples(session))
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
        stack = np.zeros((len(dates),get_maxlag_samples(session)))
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



############ session MISC ############
def get_maxlag_samples(session):
    maxlag = float(get_config(session,'maxlag'))
    return int(2*maxlag*float(get_config(session,'cc_sampling_rate'))+1)

def get_components_to_compute(session):
    components_to_compute = []
    for comp in ['ZZ','RR','TT','TR','RT','ZR','RZ','TZ','ZT']:
        if get_config(session, comp) in ['Y','y','1',1]:
            components_to_compute.append(comp)
    return components_to_compute

def build_ref_datelist(session):
    begin = get_config(session, "ref_begin")
    end = get_config(session, "ref_end")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin,'%Y-%m-%d').date()
        end = datetime.datetime.strptime(end,'%Y-%m-%d').date()
    
    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]

def build_movstack_datelist(session):
    begin = get_config(session, "startdate")
    end = get_config(session, "enddate")
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    else:
        start = datetime.datetime.strptime(begin,'%Y-%m-%d').date()
        end = datetime.datetime.strptime(end,'%Y-%m-%d').date()
    
    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]

def build_daystack_datelist(session):
    refstart, refend, refdates = build_ref_datelist(session)
    movstart, movend, movdates = build_movstack_datelist(session)
    start = min(refstart,movstart)
    end = max(refend, movend)
    r = (end+datetime.timedelta(days=1)-start).days
    return start, end, [start+datetime.timedelta(days=i) for i in range(r)]

def updated_days_for_dates(session, date1, date2, pair,type='CC',interval=datetime.timedelta(days=1),returndays=False):
    lastmod = datetime.datetime.now() - interval
    days = session.query(Job).filter(Job.day >= date1).filter(Job.day <=date2).filter(Job.type==type).filter(Job.lastmod >= lastmod).all()
    if returndays:
        return [datetime.datetime.strptime(day.day,'%Y-%m-%d').date() for day in days] ## RETURN DATE LIST !!!
    else:
        return True

############ MISCS ############

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


def nextpow2(x):
    return np.ceil(np.log2(np.abs(x)))


################## TEST
if __name__ == "__main__":
    s = connect()
    for filter in get_filters(s,False):
        print filter.ref, filter.low, filter.high

    print get_networks(s)

    for station in get_stations(s,False,net='BE'):
        print station.net, station.sta

    print get_config(s)
    print get_config(s,'data_folder')

    # low = 0.1
    # mwcs_low = 0.12
    # mwcs_high = 0.98
    # high=1.0
    # rms_threshold=0
    # mwcs_wlen=10.0
    # mwcs_step=5.0
    # used=True
    # f = Filter(low, mwcs_low, high, mwcs_high, rms_threshold, mwcs_wlen, mwcs_step, used)
    # session.add(f)
    # session.commit()
    # f.add()

