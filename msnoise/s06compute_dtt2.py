from obspy.signal.regression import linear_regression
from .api import *

import logbook


def wavg_wstd(data, errors):
    d = data
    errors[errors == 0] = 1e-6
    w = 1. / errors
    wavg = (d * w).sum() / w.sum()
    N = len(np.nonzero(w)[0])
    wstd = np.sqrt(np.sum(w * (d - wavg) ** 2) / ((N - 1) * np.sum(w) / N))
    return wavg, wstd


def main(interval=1, loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_dtt_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute DT/T ***')
    db = connect()
    params = get_params(db)

    start, end, datelist = build_movstack_datelist(db)

    mov_stacks = params.mov_stack

    components_to_compute = get_components_to_compute(db)
    updated_dtt = updated_days_for_dates(
        db, start, end, '%', jobtype='DTT', returndays=True,
        interval=datetime.timedelta(days=interval))
    interstations = {}
    for sta1, sta2 in get_station_pairs(db):
        s1 = "%s_%s" % (sta1.net, sta1.sta)
        s2 = "%s_%s" % (sta2.net, sta2.sta)
        if s1 == s2:
            interstations["%s_%s" % (s1, s2)] = 0.0
        else:
            interstations["%s_%s" % (s1, s2)] = get_interstation_distance(sta1,
                                                                          sta2,
                                                                          sta1.coordinates)

    filters = get_filters(db, all=False)
    while is_dtt_next_job(db, flag='T', jobtype='DTT'):
        # TODO would it be possible to make the next 8 lines in the API ?
        jobs = get_dtt_next_job(db, flag='T', jobtype='DTT')

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])
        netsta1, netsta2 = pair.split(':')
        station1, station2 = pair.split(":")
        n1, s1, l1 = netsta1.split(".")
        n2, s2, l2 = netsta2.split(".")
        # todo, include location code for computing distances?
        dpair = "%s_%s_%s_%s" % (n1, s1, n2, s2)
        dist = interstations[dpair] if dpair in interstations else 0.0
        logger.info(
            "There are DTT jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)
            freqmin = f.mwcs_low
            freqmax = f.mwcs_high
            if station1 == station2:
                components_to_compute = params.components_to_compute_single_station
            else:
                components_to_compute = params.components_to_compute

            for components in components_to_compute:
                for mov_stack in mov_stacks:
                    output = []
                    try:
                        mwcs = xr_get_mwcs(station1, station2, components, filterid, mov_stack)
                    except FileNotFoundError as fullpath:
                        logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                        continue

                    # todo = mwcs.index.intersection(pd.to_datetime(days))
                    # mwcs = mwcs.loc[todo]
                    # mwcs = mwcs.dropna()

                    to_search = pd.to_datetime(days)
                    to_search = to_search.append(pd.DatetimeIndex([to_search[-1] + pd.Timedelta("1d"), ]))
                    mwcs = mwcs[mwcs.index.floor('d').isin(to_search)]
                    mwcs = mwcs.dropna()


                    M = mwcs.xs("M", level=0, axis=1).copy()
                    EM = mwcs.xs("EM", level=0, axis=1).copy()
                    MCOH = mwcs.xs("MCOH", level=0, axis=1).copy()
                    tArray = M.columns.values

                    if params.dtt_lag == "static":
                        lmlag = -f.dtt_minlag
                        rmlag = f.dtt_minlag
                    else:
                        lmlag = -dist / f.dtt_v
                        rmlag = dist / f.dtt_v
                    lMlag = lmlag - f.dtt_width
                    rMlag = rmlag + f.dtt_width

                    if params.dtt_sides == "both":
                        tindex = np.where(
                            ((tArray >= lMlag) & (tArray <= lmlag)) | (
                                        (tArray >= rmlag) & (tArray <= rMlag)))[
                            0]
                    elif params.dtt_sides == "left":
                        tindex = \
                        np.where((tArray >= lMlag) & (tArray <= lmlag))[0]
                    else:
                        tindex = \
                        np.where((tArray >= rmlag) & (tArray <= rMlag))[0]
                    tmp = np.setdiff1d(np.arange(len(tArray)), tindex)
                    EM.iloc[:, tmp] = 1.0
                    MCOH.iloc[:, tmp] *= 0.0

                    MCOH[MCOH < params.dtt_mincoh] = 0.0
                    EM[EM > params.dtt_maxerr] = 1.0

                    # TODO missing check on max_dt !!


                    values = []
                    dates = []
                    for i in range(len(M.index)):
                        errArray = EM.iloc[i]
                        dtArray = M.iloc[i]
                        cohArray = MCOH.iloc[i]

                        index = np.where((errArray != 1.0) & (cohArray != 0.0))[0]
                        errArray = errArray.iloc[index]
                        dtArray = dtArray.iloc[index]

                        w = 1.0 / errArray
                        w[~np.isfinite(w)] = 1.0
                        VecXfilt = tArray[index]
                        VecYfilt = dtArray
                        if len(VecYfilt) >= 2:
                            m, a, em, ea = linear_regression(
                                VecXfilt, VecYfilt, w,
                                intercept_origin=False)
                            m0, em0 = linear_regression(
                                VecXfilt, VecYfilt, w,
                                intercept_origin=True)
                            values.append([m, em, a, ea, m0, em0])
                            dates.append(M.index[i])
                    output = pd.DataFrame(values, index=dates,
                                          columns=["m", "em", "a", "ea", "m0", "em0"])

                    xr_save_dtt(station1, station2, components, filterid, mov_stack, output)

        massive_update_job(db, jobs, "D")

    logger.info('*** Finished: Compute DTT ***')
