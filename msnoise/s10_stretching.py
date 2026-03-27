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

    $ msnoise cc dtt compute_stretching

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dtt compute_stretching

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

.. versionadded:: 1.4
    Parallel Processing
"""

from .core.db import connect, get_logger
from .core.config import get_params
from .core.stations import get_interstation_distance, get_station
from .core.workflow import (compute_rolling_ref, get_next_lineage_batch, get_t_axis, is_next_job_for_step, massive_update_job, refstack_is_rolling)
from .core.io import xr_get_ccf, xr_get_ref, xr_save_stretching

import time

import numpy as np
import pandas as pd
import xarray as xr
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
    logger = get_logger('msnoise.stretching', loglevel, with_pid=True)
    logger.info('*** Starting: Compute Stretching ***')

    db = connect()
    params = get_params(db)
    time.sleep(np.random.random() * 5)

    while is_next_job_for_step(db, step_category="stretching"):
        logger.debug("Getting the next batch")
        batch = get_next_lineage_batch(db, step_category="stretching", group_by="pair_lineage", loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        lineage_names = batch["lineage_names_upstream"]
        lineage_names_mov = batch["lineage_names_mov"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]

        taxis = get_t_axis(params)

        logger.info(f"New Stretching Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        root = params.output_folder
        mov_stacks = params.mov_stack
        goal_sampling_rate = params.cc.cc_sampling_rate

        netsta1, netsta2 = pair.split(':')
        station1, station2 = pair.split(":")
        if station1 == station2:
            components_to_compute = params.components_to_compute_single_station
        else:
            components_to_compute = params.components_to_compute
        
        for components in components_to_compute:
            station1, station2 = pair.split(":")
            rolling_mode = refstack_is_rolling(params)
            if not rolling_mode:
                # Mode A: load fixed REF from disk
                try:
                    ref = xr_get_ref(root, lineage_names,
                                     station1, station2, components, taxis)
                    ref = ref.REF.values.copy()
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue

            # zero the data outside of the minlag-maxlag timing
            if params.stretching.stretching_lag == "static":
                minlag = params.stretching.stretching_minlag
            else:
                SS1 = station1.split(".")
                SS2 = station2.split(".")

                SS1 = get_station(db, SS1[0], SS1[1])
                SS2 = get_station(db, SS2[0], SS2[1])
                minlag = get_interstation_distance(SS1, SS2,
                                                SS1.coordinates) / params.stretching.stretching_v
            maxlag2 = minlag + params.stretching.stretching_width
            mid = int(params.cc.cc_sampling_rate * params.cc.maxlag)
            ref[mid - int(minlag * goal_sampling_rate):mid + int(minlag * goal_sampling_rate)] = 0.
            ref[:mid - int(maxlag2 * goal_sampling_rate)] = 0.
            ref[mid + int(maxlag2 * goal_sampling_rate):] = 0.

            # TODO ADD the def here or in the API
            str_range = params.stretching.stretching_max
            nstr = params.stretching.stretching_nsteps
            ref_stretched, deltas = stretch_mat_creation(ref,str_range=str_range, nstr=nstr)

            for mov_stack in mov_stacks:
                #alldays = []
                #alldeltas = []
                #allcoefs = []
                #allerrs = []

                try:
                    data = xr_get_ccf(root, lineage_names_mov,
                                      station1, station2, components, mov_stack, taxis,
                                      format="dataframe")
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
                logger.debug("Processing %s:%s m%s %s" % (station1, station2, mov_stack, components))

                to_search = pd.to_datetime(days)
                data = data[data.index.floor('d').isin(to_search)]
                data = data.dropna()

                if rolling_mode:
                    ref_rolling = compute_rolling_ref(  # noqa: F841 — computed but not yet used in loop
                        data, int(params.refstack.ref_begin), int(params.refstack.ref_end)
                    )

                # print("Whitening %s" % fn)
                # data = pd.DataFrame(data)
                # data = data.apply(ww, axis=1, result_type="broadcast")

                data.iloc[:,mid - int(minlag * params.cc.cc_sampling_rate):mid + int(
                    minlag * params.cc.cc_sampling_rate)] *= 0.
                data.iloc[:,mid - int(maxlag2 * params.cc.cc_sampling_rate)] *= 0.
                data.iloc[:,mid + int(maxlag2 * params.cc.cc_sampling_rate):] *= 0.

                data_values = data.values
                num_days = data_values.shape[0]

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
                # Build xarray Dataset directly — no DataFrame round-trip
                ds_out = xr.Dataset(
                    {
                        "STR": xr.DataArray(
                            np.column_stack([alldeltas, allcoefs, allerrs]),
                            dims=["times", "keys"],
                            coords={
                                "times": alldays.values,
                                "keys":  ["Delta", "Coeff", "Error"],
                            },
                        )
                    }
                )
                xr_save_stretching(
                    root, lineage_names, step.step_name,
                    station1, station2, components, mov_stack, ds_out,
                )

        massive_update_job(db, jobs, "D")
    #    if not params.global_.hpc:
   #        for job in jobs:
     #           update_job(db, job.day, job.pair, 'DTT', 'T')
    logger.info('*** Finished: Compute Stretching ***')
