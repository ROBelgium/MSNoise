import os
import numpy as np
import scipy.signal
from scipy.stats.stats import nanmean

from database_tools import *

import logging
logging.basicConfig(level=logging.DEBUG,
                    filename="./stack.log",
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)



def stack(stype='day'):
    db = connect()
    components_to_compute = get_components_to_compute(db)
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
    
    
    if stype == "day":
        start, end, datelist = build_daystack_datelist(db)
        format = "stack"
    elif stype == "mov":
        start, end, datelist = build_movstack_datelist(db)
        format = "matrix"
        mov_stack = get_config(db,"mov_stack")
        if mov_stack.count(',') == 0:
            mov_stacks = [int(mov_stack),]
        else:
            mov_stacks = [int(mi) for mi in mov_stack.split(',')]
        if 1 in mov_stacks:
            mov_stacks.remove(1) #remove 1 day stack, it should exist already
    elif stype == "ref":
        start, end, datelist = build_ref_datelist(db)
        format = "stack"
    
    for f in get_filters(db,all=False):
        filterid = int(f.ref)
        for components in components_to_compute:
            for station1, station2 in get_station_pairs(db,used=True):
                sta1 = "%s_%s"%(station1.net, station1.sta)
                sta2 = "%s_%s"%(station2.net, station2.sta)
                pair = "%s:%s"%(sta1,sta2)
                if updated_days_for_dates(db,start,end,pair.replace('_','.'),type='CC',interval=datetime.timedelta(days=1)):
                    logging.debug("New Data for %s-%s-%i"%(pair,components,filterid))
                    if stype in ['mov','ref']:
                        nstack, stack_total = get_results(db,sta1,sta2,filterid,components,datelist,format=format)
                        updated_days = updated_days_for_dates(db,start,end,pair.replace('_','.'),type='CC',interval=datetime.timedelta(days=1),returndays=True)
                        if nstack > 0:
                            if stype == "mov":
                                for i, date in enumerate(datelist):
                                    jobadded= False
                                    for mov_stack in mov_stacks:
                                        if i < mov_stack:
                                            low = 0
                                            high = mov_stack
                                        else:
                                            low = i-mov_stack+1
                                            high = i+1
                                        newdata = False
                                        for uday in updated_days:
                                            if uday in datelist[low:high]:
                                                newdata = True
                                                break
                                        if newdata:
                                            corr = stack_total[low:high]
                                            if not np.all(np.isnan(corr)):
                                                day_name = "%s_%s"%(sta1,sta2)
                                                logging.debug("%s %s [%s - %s] (%i day stack)" % (day_name, date, datelist[low],datelist[i],mov_stack))
                                                corr = nanmean(corr,axis=0)
                                                corr = scipy.signal.detrend(corr)
                                                stack_path = os.path.join("STACKS","%02i"%filterid,"%03i_DAYS"%mov_stack,components,day_name)
                                                filename = os.path.join(stack_path,str(date))
                                                if mseed:
                                                    export_mseed(db,filename,pair,components,filterid,corr)
                                                if sac:
                                                    export_sac(db,filename,pair,components,filterid,corr)
                                                day_name = "%s:%s"%(sta1,sta2)
                                                if not jobadded:
                                                    update_job(db,date,day_name.replace('_','.'),'DTT','T')
                                                    jobadded= True
                                            del corr
                            
                            elif stype == "ref":
                                stack_path = os.path.join("STACKS","%02i"%filterid,"REF",components)
                                ref_name = "%s_%s"%(sta1,sta2)
                                filename = os.path.join(stack_path,ref_name)
                                stack_total = scipy.signal.detrend(stack_total)

                                if mseed:
                                    export_mseed(db,filename,pair,components,filterid,stack_total)
                                if sac:
                                    export_sac(db,filename,pair,components,filterid,stack_total)
                                ref_name = "%s:%s"%(sta1,sta2)
                                update_job(db,"REF",ref_name.replace('_','.'),'DTT','T')
                                del stack_total
                    elif stype == 'day':
                        updated_days = updated_days_for_dates(db,start,end,pair.replace('_','.'),type='CC',interval='1 DAY',returndays=True)
                        for date in updated_days:
                            date = date
                            print "Stacking %s day=%s"% (pair, date)
                            daystack = os.path.join("STACKS","%02i"%filterid,"001_DAYS",components,"%s_%s"%(sta1,sta2),str(date))
                            stack = np.zeros(get_maxlag_samples(db))
                            ncorr = 0
                            try:
                                os.makedirs(os.path.split(daystack)[0])
                            except:
                                pass
                            path = os.path.join(output_folder, "%02i"% filterid, sta1,sta2,components, str(date))
                            if os.path.isdir(path):
                                for file in os.listdir(path):
                                    if len(file) == 8:
                                        st = read(os.path.join(path,file))
                                        if not np.any(np.isnan(st[0].data)) and not np.any(np.isinf(st[0].data)):
                                            stack += st[0].data
                                            ncorr += 1
                                if ncorr > 0:
                                    if mseed:
                                        export_mseed(db, daystack, pair, components, filterid,stack/ncorr,ncorr)
                                    if sac:
                                        export_sac(db, daystack, pair, components, filterid,stack/ncorr,ncorr)
                            del stack
def daystack():
    stack(stype='day')

def refstack():
    stack(stype='ref')

def movstack():
    stack(stype='mov')

logging.info('*** Starting: Stack ***')
db = connect()
# daystack()
refstack()
movstack()
logging.info('*** Finished: Stack ***')


# EOF
