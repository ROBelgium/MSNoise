from obspy.core import read

from database_tools import *
from MWCS import mwcs
import logging
logging.basicConfig(level=logging.DEBUG,
                    filename="./compute_mwcs.log",
                    format='%(asctime)s [%(levelname)s] %(message)s',
                    filemode='w')

console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

logging.info('*** Starting: Compute MWCS ***')

db = connect()
components_to_compute = get_components_to_compute(db)

mov_stack = get_config(db,"mov_stack")
if mov_stack.count(',') == 0:
    mov_stacks = [int(mov_stack),]
else:
    mov_stacks = [int(mi) for mi in mov_stack.split(',')]

goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
maxlag = float(get_config(db, "maxlag"))
start, end, datelist = build_movstack_datelist(db)
allow_large_concats(db)

# First we reset all DTT jobs to "T"odo if the REF is new for a given pair

for station1, station2 in get_station_pairs(db,used=True):
    sta1 = "%s.%s"%(station1.net, station1.sta)
    sta2 = "%s.%s"%(station2.net, station2.sta)
    pair = "%s:%s"%(sta1,sta2)
    if is_dtt_ref_job(db,pair,type='DTT'):
        logging.info("We will recompute all MWCS based on the new REF for %s" % pair)
        reset_dtt_jobs(db,pair)
        update_job(db,"REF",pair,type='DTT',flag='D')

# Then we compute the jobs
while is_dtt_mov_job(db,type='DTT'):
    pair,days, refs = get_dtt_next_job(db, flag='T', type='DTT')
    logging.info("There are MWCS jobs for some days to recompute for %s"%pair)
    for day in days.split(','):
        for f in get_filters(db,all=False):
            filterid = int(f.ref)
            for components in components_to_compute:
                ref_name = pair.replace('.','_').replace(':','_')
                rf = os.path.join("STACKS","%02i"%filterid,"REF",components,ref_name+".MSEED")
                ref = read(rf)[0].data
                for mov_stack in mov_stacks:
                    df = os.path.join("STACKS","%02i"%filterid,"%03i_DAYS"%mov_stack,components,ref_name,str(day)+".MSEED")
                    if os.path.isfile(df):
                        cur = read(df)[0].data
                        logging.debug('Processing MWCS for: %s.%s.%02i - %s - %02i days'%(ref_name, components, filterid, day, mov_stack))
                        output = mwcs(cur,ref,f.mwcs_low,f.mwcs_high,goal_sampling_rate,-maxlag,f.mwcs_wlen,f.mwcs_step)
                        outfolder = os.path.join('MWCS',"%02i"%filterid,"%03i_DAYS"%mov_stack,components,ref_name)
                        if not os.path.isdir(outfolder):
                            os.makedirs(outfolder)
                        np.savetxt(os.path.join(outfolder,"%s.txt"%str(day)),output)
                        del output
        update_job(db,day,pair,type='DTT',flag='D')

logging.info('*** Finished: Compute MWCS ***')


