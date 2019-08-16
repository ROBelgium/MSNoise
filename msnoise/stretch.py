"""
Stretching...
"""

from .api import *
from scipy import asarray as ar
from scipy.optimize import curve_fit
from scipy.ndimage import map_coordinates

def stretch_mat_creation(refcc, str_range=0.01, nstr=1001):
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
        strrefmat[i, :] = map_coordinates(refcc.reshape((len(refcc), 1)), coord)
    return strrefmat, strvec


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('*** Starting: Compute STR ***')

    db = connect()
    components_to_compute = get_components_to_compute(db)

    mov_stack = get_config(db, "mov_stack")
    if mov_stack.count(',') == 0:
        mov_stacks = [int(mov_stack), ]
    else:
        mov_stacks = [int(mi) for mi in mov_stack.split(',')]

    goal_sampling_rate = float(get_config(db, "cc_sampling_rate"))
    maxlag = float(get_config(db, "maxlag"))
    export_format = get_config(db, 'export_format')
    if export_format == "BOTH":
        extension = ".MSEED"
    else:
        extension = "."+export_format

    # First we reset all DTT jobs to "T"odo if the REF is new for a given pair
    # for station1, station2 in get_station_pairs(db, used=True):
    #     sta1 = "%s.%s" % (station1.net, station1.sta)
    #     sta2 = "%s.%s" % (station2.net, station2.sta)
    #     pair = "%s:%s" % (sta1, sta2)
    #     if is_dtt_next_job(db, jobtype='DTT', ref=pair):
    #         logging.info(
    #           "We will recompute all STR based on the new REF for %s" % pair)
    #         reset_dtt_jobs(db, pair)
    #         update_job(db, "REF", pair, jobtype='DTT', flag='D')

    filters = get_filters(db, all=False)
    params = get_params(db)
    # Then we compute the jobs
    while is_dtt_next_job(db, flag='T', jobtype='MWCS'):
        jobs = get_dtt_next_job(db, flag='T', jobtype='MWCS')

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        logging.info(
            "There are MWCS jobs for some days to recompute for %s" % pair)
        
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

        maxlag2 = minlag + dtt_width
        print("betweeen", minlag, "and", maxlag2)

        for f in get_filters(db, all=False):
            filterid = int(f.ref)

            for mov_stack in mov_stacks:
                for components in params.all_components:
                        rf = os.path.join("STACKS", "%02i" %
                                          filterid, "REF", components, ref_name + extension)
                        if os.path.isfile(rf):
                            ref = read(rf)[0].data
                            mid = int(goal_sampling_rate*maxlag)
                            ref[mid-int(minlag*goal_sampling_rate):mid+int(minlag*goal_sampling_rate)] *= 0.
                            ref[:mid-int(maxlag2*goal_sampling_rate)] *= 0.
                            ref[mid+int(maxlag2*goal_sampling_rate):] *= 0.
                        else:
                            logging.debug(
                                "No REF file named %s, skipping." % rf)
                            continue
                        alldays = []
                        alldeltas = []
                        allcoefs = []
                        allerrs = []
                        str_range = 0.5  ### HARD CODE!!! ###
                        nstr = 1001  ### HARD CODE!!! ###
                        ref_stretched, deltas = stretch_mat_creation(ref,
                                                                     str_range=str_range,
                                                                     nstr=nstr)
                        for day in days:
                            df = os.path.join(
                                "STACKS", "%02i" % filterid, "%03i_DAYS" %
                                mov_stack, components, ref_name, str(day) + extension)

                            if os.path.isfile(df):
                                cur = read(df)[0].data   ### read the current mseed file ###
                                cur[mid-int(minlag*goal_sampling_rate):mid+int(minlag*goal_sampling_rate)] *= 0.
                                cur[:mid-int(maxlag2*goal_sampling_rate)] *= 0.
                                cur[mid+int(maxlag2*goal_sampling_rate):] *= 0.  ### replace with zeroes at all
                                                                                 ### times outside minlag to maxlag
                                logging.debug(
                                    'Processing Stretching for: %s.%s.%02i - %s - %02i days' %
                                    (ref_name, components, filterid, day, mov_stack))

                                coeffs =[]
                                for i in range(ref_stretched.shape[0]):
                                    ci = np.corrcoef(cur,ref_stretched[i])[0,1]
                                    coeffs.append(ci)

                                tday = datetime.datetime.strptime(day, "%Y-%m-%d")
                                alldays.append(tday)
                                alldeltas.append(deltas[np.argmax(coeffs)])
                                allcoefs.append(np.max(coeffs))

                                ###### gaussian fit ######
                                def gauss_function(x, a, x0, sigma):
                                    return a*np.exp(-(x-x0)**2/(2*sigma**2))
                                x = ar(range(len(coeffs)))
                                ymax_index = coeffs.index(np.max(coeffs))
                                ymin = np.min(coeffs)
                                coeffs_shift = []
                                for i in coeffs:
                                    i += np.absolute(ymin) # make all points above zero
                                    coeffs_shift.append(i)
                                n = len(coeffs)
                                x0 = sum(x)/n
                                sigma = (sum((x-x0)**2)/n)**0.5
                                try:
                                    popt, pcov = curve_fit(gauss_function, x, coeffs_shift, [ymax_index, x0, sigma])
                                    FWHM = 2 * ((2*np.log(2))**0.5)*popt[2] # convert sigma (popt[2]) to FWHM
                                    error = FWHM / 2  ### error is half width at full maximum
                                except RuntimeError:
                                    error = np.nan # gaussian fit failed

                                allerrs.append(error)

                        df = pd.DataFrame(np.array([alldeltas,allcoefs,allerrs]).T, index=alldays, columns=["Delta", "Coeff", "Error"],)
                        output = os.path.join('STR', "%02i" % filterid, "%03i_DAYS" % mov_stack, components)
                        if not os.path.isdir(output):
                            os.makedirs(output)
                        df.to_csv(os.path.join(output, "%s.csv" % ref_name), index_label="Date")

        # THIS SHOULD BE IN THE API
        updated = False
        mappings = [{'ref': job.ref, 'flag': "D"} for job in jobs]
        while not updated:
            try:
                db.bulk_update_mappings(Job, mappings)
                db.commit()
                updated = True
            except:
                time.sleep(np.random.random())
                pass

    logging.info('*** Finished: Compute STR ***')

if __name__ == "__main__":
    main()
