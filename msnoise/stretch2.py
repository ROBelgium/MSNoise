"""
.. warning:: if using only ``mov_stack`` = 1, no STR jobs is inserted in the
    database and consequently, no STR calculation will be done! FIX!


Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``mwcs_low``: The lower frequency bound of the linear regression done in
  MWCS (in Hz)
* ``mwcs_high``: The upper frequency bound of the linear regression done in
  MWCS (in Hz)
* ``mwcs_wlen``: Window length (in seconds) to perform MWCS
* ``mwcs_step``: Step (in seconds) of the windowing procedure in MWCS

* |hpc| | *new in 1.6*

In short, both time series are sliced in several overlapping windows and
preprocessed. The similarity of the two time-series is assessed using the
cross-coherence between energy densities in the frequency domain. The time
delay between the two cross correlations is found in the unwrapped phase of
the cross spectrum and is linearly proportional to frequency. This "Delay" for
each window between two signals is the slope of a weighted linear regression
(WLS) of the samples within the frequency band of interest.

For each filter, the frequency band can be configured using ``mwcs_low``
and ``mwcs_high``, and the window and overlap lengths using ``mwcs_wlen`` and
``mwcs_step``.

The output of this process is a table of delays measured at each window in the
functions. The following is an example for lag times between -115 and -90.
In this case, the window length was 10 seconds with an overlap of 5 seconds.

.. code-block:: python

          LAG_TIME          DELAY           ERROR         MEAN COHERENCE
    -1.1500000000e+02 -1.4781146383e-01 5.3727119135e-02 2.7585243911e-01
    -1.1000000000e+02 -6.8207526992e-02 2.0546644311e-02 3.1620999352e-01
    -1.0500000000e+02 -1.0337029577e-01 8.6645155402e-03 4.2439269880e-01
    -1.0000000000e+02 -2.8668775696e-02 6.2522215988e-03 5.7159849528e-01
    -9.5000000000e+01  4.1803941008e-02 1.5102285789e-02 4.1238557789e-01
    -9.0000000000e+01  4.8139400233e-02 3.2700657018e-02 3.0586187792e-01

This process is job-based, so it is possible to run several instances in
parallel.

Once done, each job is marked "D"one in the database and, unless ``hpc`` is 
``Y``, DTT jobs are inserted/updated in the database.

To run this step:

.. code-block:: sh

    $ msnoise cc dvv compute_stretching

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dvv compute_stretching

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

.. versionadded:: 1.4
    Parallel Processing
"""

from .api import *

import logbook
from scipy.fft import next_fast_len
from numpy import asarray as ar
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
    :param str_range: Amount of the desired stretching (one side)
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


def main(loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_str_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute Stretching ***')

    db = connect()
    params = get_params(db)
    export_format = params.export_format
    if export_format == "BOTH":
        extension = ".MSEED"
    else:
        extension = "." + export_format
    mov_stacks = params.mov_stack

    goal_sampling_rate = params.cc_sampling_rate
    maxlag = params.maxlag

    # First we reset all DTT jobs to "T"odo if the REF is new for a given pair
    # for station1, station2 in get_station_pairs(db, used=True):
    #     sta1 = "%s.%s" % (station1.net, station1.sta)
    #     sta2 = "%s.%s" % (station2.net, station2.sta)
    #     pair = "%s:%s" % (sta1, sta2)
    #     if is_dtt_next_job(db, jobtype='DTT', ref=pair):
    #         logger.info(
    #             "We will recompute all MWCS based on the new REF for %s" % pair)
    #         reset_dtt_jobs(db, pair)
    #         update_job(db, "REF", pair, jobtype='DTT', flag='D')
    # 
    logger.debug('Ready to compute')
    # Then we compute the jobs
    outfolders = []
    filters = get_filters(db, all=False)
    time.sleep(np.random.random() * 5)
    taxis = get_t_axis(db)
    smoothing_half_win= 5
    # hanningwindow = get_window("hanning", smoothing_half_win)
    while is_dtt_next_job(db, flag='T', jobtype='STR'):
        # TODO would it be possible to make the next 8 lines in the API ?
        jobs = get_dtt_next_job(db, flag='T', jobtype='STR')

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        logger.info(
            "There are STR (stretching) jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)
            freqmin = f.mwcs_low
            freqmax = f.mwcs_high
            low = f.low
            high = f.high

            def ww(a):
                from .move2obspy import whiten
                n = next_fast_len(len(a))
                return whiten(a, n, 1./params.cc_sampling_rate,
                              low, high, returntime=True)
            for components in params.all_components:
                ref_name = pair.replace(':', '_')
                station1, station2 = pair.split(":")
                try:
                    ref = xr_get_ref(station1, station2, components, filterid,
                                     taxis)
                    ref = ref.CCF.values
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
                # print("Whitening ref")
                # ref = ww(ref)

                # zero the data outside of the minlag-maxlag timing
                if params.dtt_lag == "static":
                    minlag = params.dtt_minlag
                else:
                    SS1 = station1.split(".")
                    SS2 = station2.split(".")

                    SS1 = get_station(db, SS1[0], SS1[1])
                    SS2 = get_station(db, SS2[0], SS2[1])
                    minlag = get_interstation_distance(SS1, SS2,
                                                       SS1.coordinates) / params.dtt_v
                maxlag2 = minlag + params.dtt_width
                mid = int(params.goal_sampling_rate * params.maxlag)
                print("betweeen", minlag, "and", maxlag2    )
                ref[mid - int(minlag * goal_sampling_rate):mid + int(
                    minlag * goal_sampling_rate)] *= 0.
                ref[:mid - int(maxlag2 * goal_sampling_rate)] *= 0.
                ref[mid + int(maxlag2 * goal_sampling_rate):] *= 0.
                # TODO ADD the def here or in the API
                str_range = params.stretching_max
                nstr = params.stretching_nsteps
                ref_stretched, deltas = stretch_mat_creation(ref,str_range=str_range, nstr=nstr)

                for mov_stack in mov_stacks:
                    output = []
                    #alldays = []
                    #alldeltas = []
                    #allcoefs = []
                    #allerrs = []

                    fn = r"STACKS2/%02i/%s_%s/%s/%s_%s.nc" % (
                    filterid, mov_stack[0], mov_stack[1], components, station1, station2)
                    print("Reading %s" % fn)
                    if not os.path.isfile(fn):
                        print("FILE DOES NOT EXIST: %s, skipping" % fn)
                        continue
                    data = xr_create_or_open(fn)
                    data = data.CCF.to_dataframe().unstack().droplevel(0, axis=1)
                    to_search = pd.to_datetime(days)
                    data = data[data.index.floor('d').isin(to_search)]
                    data = data.dropna()

                    # print("Whitening %s" % fn)
                    # data = pd.DataFrame(data)
                    # data = data.apply(ww, axis=1, result_type="broadcast")

                    data.iloc[:,mid - int(minlag * params.goal_sampling_rate):mid + int(
                        minlag * params.goal_sampling_rate)] *= 0.
                    data.iloc[:,mid - int(maxlag2 * params.goal_sampling_rate)] *= 0.
                    data.iloc[:,mid + int(maxlag2 * params.goal_sampling_rate):] *= 0.

                    data_values = data.values
                    num_days = data_values.shape[0]
                    num_stretch = ref_stretched.shape[0]

                    # Normalizing the data for correlation
                    data_norm = (data_values - data_values.mean(axis=1, keepdims=True)) / data_values.std(axis=1, keepdims=True)
                    ref_stretched_norm = (ref_stretched - ref_stretched.mean(axis=1, keepdims=True)) / ref_stretched.std(axis=1, keepdims=True)

                    # Compute the correlation coefficients
                    corr_coeffs = np.dot(ref_stretched_norm, data_norm.T) / ref_stretched.shape[1]

                    max_corr_indices = np.argmax(corr_coeffs, axis=0)
                    max_corr_values = corr_coeffs[max_corr_indices, np.arange(num_days)]

                    alldays = data.index
                    alldeltas = deltas[max_corr_indices]
                    allcoefs = max_corr_values

                    allerrs = []

                    for day_idx in range(data_values.shape[0]):

                        ###### gaussian fit ######
                        def gauss_function(x, a, x0, sigma):
                            return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

                        coeffs = corr_coeffs[:, day_idx]
                        x = ar(range(len(coeffs)))
                        ymax_index = np.argmax(coeffs)
                        ymin = np.min(coeffs)
                        coeffs_shift = []
                        for j in coeffs:
                            j += np.absolute(ymin)  # make all points above zero
                            coeffs_shift.append(j)
                        n = len(coeffs)
                        x0 = sum(x) / n
                        sigma = (sum((x - x0) ** 2) / n) ** 0.5
                        try:
                            popt, pcov = curve_fit(gauss_function, x,
                                                   coeffs_shift,
                                                   [ymax_index, x0, sigma])
                            FWHM = 2 * ((2 * np.log(2)) ** 0.5) * popt[
                                2]  # convert sigma (popt[2]) to FWHM
                            error = FWHM / 2  ### error is half width at full maximum
                        except RuntimeError:
                            error = np.nan  # gaussian fit failed

                        allerrs.append(error)
                    df = pd.DataFrame(
                        np.array([alldeltas, allcoefs, allerrs]).T,
                        index=alldays, columns=["Delta", "Coeff", "Error"], )
                    output = os.path.join('STR2', "%02i" % filterid,
                                          "%s_%s" % (mov_stack[0],mov_stack[1]), components)
                    # print(df.head())
                    if not os.path.isdir(output):
                        os.makedirs(output)
                    fn = os.path.join(output, "%s.csv" % ref_name)
                    if not os.path.isfile(fn):
                        df.to_csv(fn, index_label="Date")
                    else:
                        dest = pd.read_csv(fn, index_col=0,
                                           parse_dates=True)
                        final = pd.concat([dest, df])
                        final = final[~final.index.duplicated(keep='last')]
                        final = final.sort_index()
                        final.to_csv(fn, index_label="Date")
                    ### TODO END

        massive_update_job(db, jobs, "D")
        if not params.hpc:
            for job in jobs:
                update_job(db, job.day, job.pair, 'DTT', 'T')
    logger.info('*** Finished: Compute Stretching ***')
