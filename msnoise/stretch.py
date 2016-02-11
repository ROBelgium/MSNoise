"""
Strechting...
"""

from obspy.core import read
import numpy as np
# import pandas as pd
from .api import *
# from MWCS import mwcs
import logging
import matplotlib.pyplot as plt


def stretch_mat_creation(refcc, str_range=0.01, nstr=1000):
    """ Matrix of stretched instance of a reference trace.

    The reference trace is stretched using a cubic spline interpolation
    algorithm form ``-str_range`` to ``str_range`` (in %) for totally
    ``nstr`` steps.
    The output of this function is a matrix containing the stretched version
    of the reference trace (one each row) (``strrefmat``) and the corresponding
    stretching amount (`strvec```).

    :type refcc: :class:`~numpy.ndarray`
    :param refcc: 1d ndarray. The reference trace that will be stretched
    :type str_range: float
    :param str_range: Amount, in percent, of the desired stretching (one side)
    :type nstr: int
    :param nstr: Number of stretching steps (one side)

    :rtype: :class:`~numpy.ndarray` and float
    :return: **strrefmat**:
        - 2d ndarray of stretched version of the reference trace.
        Its size is ``(nstr,len(refcc)/2)`` if ``signle_side==True``
        otherwise it is ``(nstr,len(refcc))``
    :rtype: float
    :return: **strvec**: List of float, stretch amount for each row
        of ``strrefmat``
    """

    n = len(refcc)
    samples_idx = np.arange(n) - n // 2
    strvec = 1 + np.linspace(-str_range, str_range, nstr)
    str_timemat = np.zeros((nstr, n))
    for ii in np.arange(nstr):
        str_timemat[ii, :] = samples_idx / strvec[nstr - 1 - ii]
    strrefmat = np.zeros_like(str_timemat)
    coord = np.zeros((2, n))
    for (i, row) in enumerate(str_timemat):
        coord[0, :] = row + n // 2
        strrefmat[i, :] = scipy.ndimage.map_coordinates(\
                          refcc.reshape((len(refcc), 1)), coord)
    return strrefmat, strvec

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: Compute MWCS ***')
    
    db = connect()
    components_to_compute = get_components_to_compute(db)
    
    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]
    
    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    maxlag = float(get_config(db, "maxlag"))
    start, end, datelist = build_movstack_datelist(db)
    # allow_large_concats(db)
    
    # First we reset all DTT jobs to "T"odo if the REF is new for a given pair
    
    #~ for station1, station2 in get_station_pairs(db, used=True):
        #~ sta1 = "%s.%s" % (station1.net, station1.sta)
        #~ sta2 = "%s.%s" % (station2.net, station2.sta)
        #~ pair = "%s:%s" % (sta1, sta2)
        #~ if is_dtt_next_job(db, type='DTT', ref=pair):
            #~ logging.info(
                #~ "We will recompute all MWCS based on the new REF for %s" % pair)
            #~ reset_dtt_jobs(db, pair)
            #~ update_job(db, "REF", pair, type='DTT', flag='D')
    
    # Then we compute the jobs
    
    while is_dtt_next_job(db, flag='T', jobtype='DTT'):
        pair, days, refs = get_dtt_next_job(db, flag='T', jobtype='DTT')
        logging.info(
            "There are MWCS jobs for some days to recompute for %s" % pair)
        alldays = []
        alldeltas = []
        allcoefs = []
        ref_name = pair.replace('.', '_').replace(':', '_')
        sta1, sta2 = pair.split(':')
        station1 = sta1.split(".")
        station2 = sta2.split(".")
        
        station1 = get_station(db, station1[0], station1[1])
        station2 = get_station(db, station2[0], station2[1])
        
        
        dtt_lag = get_config(db, "dtt_lag")
        dtt_v = float(get_config(db, "dtt_v"))
        dtt_minlag = float(get_config(db, "dtt_minlag"))
        dtt_width = float(get_config(db, "dtt_width"))
        dtt_sides = get_config(db, "dtt_sides")
            
        if dtt_lag == "static":
            minlag = dtt_minlag
        else:
            minlag = get_interstation_distance(station1, station2, station1.coordinates) / dtt_v
            print(minlag)
        
        maxlag2 = minlag + dtt_width
        print("betweeen", minlag, "and", maxlag2)
        
        rf = os.path.join("STACKS", "%02i" %
                          3, "REF", "ZZ", ref_name + ".MSEED")
        if os.path.isfile(rf):
            ref = read(rf)[0].data
            mid = int(goal_sampling_rate*maxlag)
            ref[mid-int(minlag*goal_sampling_rate):mid+int(minlag*goal_sampling_rate)] *= 0.
            ref[:mid-int(maxlag2*goal_sampling_rate)] *= 0.
            ref[mid+int(maxlag2*goal_sampling_rate):] *= 0.
            #~ plt.plot(ref)
            #~ plt.show()
            str_range = 0.01
            nstr = 1001
            ref_stretched, deltas = stretch_mat_creation(ref,str_range=str_range, nstr=nstr)
            for day in days:
                # logging.debug('Day=%s'%day)
                filters = get_filters(db, all=False)
                if 1:
                    f = filters[2]
                    filterid = int(f.ref)
                    for components in components_to_compute:
                        #removed ref from here
                        
                        for mov_stack in [20,]:
                            df = os.path.join(
                                "STACKS", "%02i" % filterid, "%03i_DAYS" %
                                mov_stack, components, ref_name, str(day) + ".MSEED")
                            if os.path.isfile(df):
                                cur = read(df)[0].data
                                
                                cur[mid-int(minlag*goal_sampling_rate):mid+int(minlag*goal_sampling_rate)] *= 0.
                                cur[:mid-int(maxlag2*goal_sampling_rate)] *= 0.
                                cur[mid+int(maxlag2*goal_sampling_rate):] *= 0.
                                logging.debug(
                                    'Processing Stretching for: %s.%s.%02i - %s - %02i days' %
                                    (ref_name, components, filterid, day, mov_stack))
                                coeffs =[]
                                for i in range(ref_stretched.shape[0]):
                                    ci = np.corrcoef(cur,ref_stretched[i])[0,1]
                                    coeffs.append(ci)
                                # output = mwcs(
                                    # cur, ref, f.mwcs_low, f.mwcs_high, goal_sampling_rate, -maxlag, f.mwcs_wlen, f.mwcs_step)
                                # outfolder = os.path.join(
                                    # 'STRETCH', "%02i" % filterid, "%03i_DAYS" % mov_stack, components, ref_name)
                                # if not os.path.isdir(outfolder):
                                    # os.makedirs(outfolder)
                                # np.savetxt(
                                    # os.path.join(outfolder, "%s.txt" % str(day)), output)
                                # del output

                                if mov_stack == 20:
                                    tday = datetime.datetime.strptime(day, "%Y-%m-%d")
                                    # print type(day), day
                                    alldays.append(tday)
                                    alldeltas.append(deltas[np.argmax(coeffs)])
                                    allcoefs.append(np.max(coeffs))
                                    #print day, deltas[np.argmax(coeffs)]

                                # plt.subplot(211)
                                # plt.figure()
                                # plt.plot(ref)
                                # plt.plot(cur)
                                # plt.plot(ref_stretched[np.argmax(coeffs)])
                                # plt.show()
                                # plt.subplot(223)
                                # plt.imshow(ref_stretched,aspect='auto',origin='upper',cmap='bwr')
                                # plt.grid()
                                # plt.axhline(nstr/2)
                                # plt.subplot(224)
                                # plt.plot(coeffs,np.arange(len(coeffs)))
                                # plt.axhline(np.argmax(coeffs))
                                # plt.title(deltas[np.argmax(coeffs)])
                                # plt.ylim(nstr-1,0)
                                # plt.grid()
                                # plt.show()
                            
                update_job(db, day, pair, jobtype='DTT', flag='D')
            #~ print np.array(alldeltas).shape
            #~ print allcoefs
            df = pd.DataFrame(np.array([alldeltas,allcoefs]).T, index=alldays)
            print(df.head())
            df.to_csv(os.path.join("STR", "%s.csv"%ref_name))
            #~ p = plt.plot(alldays, np.array(alldeltas).flatten(), c='k')
            #~ plt.scatter(alldays, np.array(alldeltas).flatten(), c=allcoefs,s=200)
            #~ plt.colorbar()
            #~ plt.title(pair)
            #~ plt.show()

    logging.info('*** Finished: Compute MWCS ***')

if __name__ == "__main__":
    main()