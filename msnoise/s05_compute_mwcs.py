"""
.. warning:: if using only ``mov_stack`` = 1, no MWCS jobs is inserted in the
    database and consequently, no MWCS calculation will be done! FIX!

Following Clarke et al (2011), we apply the :ref:`mwcs`
to study the relative dephasing between Moving-Window stacks ("Current") and a
Reference using Moving-Window Cross-Spectral analysis. The *jobs* "T"o do have
been inserted in the datavase during the stack procedure.

Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``freqmin``: The lower frequency bound of the linear regression done in
  MWCS (in Hz)
* ``freqmax``: The upper frequency bound of the linear regression done in
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

For each filter, the linear regression done in MWCS is performed between
``freqmin`` and ``freqmax`` and the window and overlap lengths configured
using ``mwcs_wlen`` and ``mwcs_step``.

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

    $ msnoise cc dtt compute_mwcs

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dtt compute_mwcs

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

.. versionadded:: 1.4
    Parallel Processing
"""

import time

import numpy as np
import xarray as xr
from obspy.signal.invsim import cosine_taper
from obspy.signal.regression import linear_regression
import scipy
import scipy.fft as sf
from scipy.fft import next_fast_len
import scipy.optimize
import scipy.signal

from .api import (
    compute_rolling_ref,
    connect,
    extend_days,
    getCoherence,
    get_logger,
    get_next_lineage_batch,
    get_t_axis,
    get_window,
    is_next_job_for_step,
    massive_update_job,
    nextpow2,
    refstack_is_rolling,
    xr_get_ccf,
    xr_get_ref,
    xr_save_mwcs,
)

def main(loglevel="INFO"):
    logger = get_logger('msnoise.mwcs', loglevel, with_pid=True)
    logger.info('*** Starting: Compute MWCS ***')

    db = connect()
    logger.debug('Ready to compute')
    # Then we compute the jobs

    smoothing_half_win= 5
    hanningwindow = get_window("hanning", smoothing_half_win)

    while is_next_job_for_step(db, step_category="mwcs"):
        logger.debug("Getting the next batch")
        batch = get_next_lineage_batch(db, step_category="mwcs", group_by="pair_lineage", loglevel=loglevel)
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

        logger.info(f"New MWCS Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        mov_stacks = params.mov_stack

        goal_sampling_rate = params.cc.cc_sampling_rate

        logger.info(
            "There are MWCS jobs for some days to recompute for %s" % pair)

        # mwcsid = int(params.ref)
        freqmin = params.mwcs.freqmin
        freqmax = params.mwcs.freqmax

        def ww(a):
            from .move2obspy import whiten
            n = next_fast_len(len(a))
            return whiten(a, n, 1./params.cc.cc_sampling_rate,
                          freqmin, freqmax, returntime=True)
        station1, station2 = pair.split(":")

        if station1 == station2:
            components_to_compute = params.components_to_compute_single_station
        else:
            components_to_compute = params.components_to_compute

        for components in components_to_compute:
            rolling_mode = refstack_is_rolling(params)
            if not rolling_mode:
                # Mode A: load fixed REF from disk (under refstack_M folder)
                try:
                    ref = xr_get_ref(params.output_folder, lineage_names,
                                     station1, station2, components, taxis)
                    ref = ref.REF.values
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
                if not len(ref):
                    continue
            # Mode B: ref computed per time-step below after data is loaded

            for mov_stack in mov_stacks:
                output = []
                try:
                    data = xr_get_ccf(params.output_folder, lineage_names_mov,
                                      station1, station2, components, mov_stack, taxis,
                                      format="dataframe")
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
                logger.debug("Processing %s:%s m%s %s" % (station1, station2, mov_stack, components))
                # todo = data.index.intersection()
                # data = data.loc[todo]

                to_search = extend_days(days)
                # data = data[(data.index.floor('d').isin(to_search) or data.index.ceil('d').isin(to_search))]
                data = data[data.index.floor('d').isin(to_search)]
                data = data.dropna()

                if rolling_mode:
                    # Mode B: compute per-index rolling reference from MOV data
                    ref_rolling = compute_rolling_ref(
                        data, int(params.refstack.ref_begin), int(params.refstack.ref_end)
                    )
                    # ref_rolling shape: (n_times, n_lag_samples)

                # print("Whitening %s" % fn)
                # data = pd.DataFrame(data)
                # data = data.apply(ww, axis=1, result_type="broadcast")

                # work on 2D mwcs:
                window_length_samples = int(
                    params.mwcs.mwcs_wlen * goal_sampling_rate)

                padd = int(2 ** (nextpow2(window_length_samples) + 2))

                count = 0
                tp = cosine_taper(window_length_samples, 0.85)
                minind = 0
                maxind = window_length_samples
                step_samples = int(params.mwcs.mwcs_step * goal_sampling_rate)
                if step_samples != (params.mwcs.mwcs_step * goal_sampling_rate):
                    logger.warning('mwcs_step of %g s incompatible with %i Hz sampling rate. Step size of %g s used instead' %
                                    (params.mwcs.mwcs_step, goal_sampling_rate, step_samples/goal_sampling_rate))

                freq_vec = sf.fftfreq(padd, 1. / goal_sampling_rate)[
                        :padd // 2]
                # Find the values the frequency range of interest
                index_range = np.argwhere(
                    np.logical_and(freq_vec >= freqmin,
                                freq_vec <= freqmax)).flatten()
                cci = np.empty((data.shape[0], window_length_samples))
                while maxind <= data.shape[1]:
                    cci[:] = data.iloc[:,
                            minind:(minind + window_length_samples)].values
                    scipy.signal.detrend(cci, type="linear", axis=1,
                                        overwrite_data=True)
                    for i in range(cci.shape[0]):
                        cci[i] *= tp

                    if rolling_mode:
                        # Mode B: per-row reference (n_times, window_length_samples)
                        cri = ref_rolling[:, minind:(minind + window_length_samples)].copy()
                        scipy.signal.detrend(cri, type='linear', axis=1, overwrite_data=True)
                        for _i in range(cri.shape[0]):
                            cri[_i] *= tp
                    else:
                        cri = ref[
                            minind:(minind + window_length_samples)].copy()
                        scipy.signal.detrend(cri, type='linear',
                                            overwrite_data=True)
                        cri *= tp

                    minind += step_samples
                    maxind += step_samples

                    fcur = sf.fft(cci, axis=1, n=padd)[:, :padd // 2]
                    if rolling_mode:
                        # cri is 2D (n_times, window); fref shape: (n_times, padd//2)
                        fref = sf.fft(cri, axis=1, n=padd)[:, :padd // 2]
                    else:
                        fref = sf.fft(cri, n=padd)[:padd // 2]

                    fcur2 = np.real(fcur) ** 2 + np.imag(fcur) ** 2
                    fcur2 = fcur2.astype(float)
                    fref2 = np.real(fref) ** 2 + np.imag(fref) ** 2
                    fref2 = fref2.astype(float)

                    X = fref * fcur.conj()
                    if smoothing_half_win != 0:
                        for i in range(fcur2.shape[0]):
                            fcur2[i] = np.sqrt(
                                scipy.signal.convolve(fcur2[i],
                                                    hanningwindow.real,
                                                    "same"))

                        fref2 = np.sqrt(scipy.signal.convolve(fref2,
                                                            hanningwindow.real,
                                                            "same"))

                        for i in range(X.shape[0]):
                            X[i] = scipy.signal.convolve(X[i],
                                                        hanningwindow,
                                                        "same")

                    else:
                        fcur2 = fcur2.apply(np.sqrt)
                        fref2 = fref2.apply(np.sqrt)

                    dcs = np.abs(X)

                    # Get Coherence and its mean value
                    W = []
                    MCOH = []
                    for i in range(dcs.shape[0]):
                        coh = getCoherence(dcs[i, index_range],
                                        fref2[index_range],
                                        fcur2[i, index_range])
                        mcoh = np.mean(coh)
                        MCOH.append(np.real(mcoh))
                        # Get Weights, avoid zero division here below:
                        coh[coh==1] = 1.0-1e-9
                        w = 1.0 / (1.0 / (coh ** 2) - 1.0)
                        w[coh >= 0.99] = 1.0 / (1.0 / 0.9801 - 1.0)
                        w = np.sqrt(w * np.sqrt(dcs[i][index_range]))
                        w = np.real(w)
                        W.append(w)

                    W = np.asarray(W)
                    #                     # Frequency array:
                    v = np.real(freq_vec[index_range]) * 2 * np.pi

                    # Phase:

                    phi = np.angle(X)
                    phi = phi.astype(float)
                    phi[:, 0] = 0.0
                    phi = np.unwrap(phi, axis=1)
                    phi = phi[:, index_range]

                    # Calculate the slope with a weighted least square linear regression
                    # forced through the origin
                    # weights for the WLS must be the variance !

                    result = np.array([linear_regression(v.flatten(),
                                                        phi[i].flatten(),
                                                        W[i].flatten(),
                                                        ) for i in
                                    range(phi.shape[0])])

                    M = result[:, 0]
                    e = np.sum((phi - np.outer(M, v)) ** 2, axis=1) / (
                                len(v) - 1)
                    s2x2 = np.sum(v ** 2 * W ** 2, axis=1)
                    sx2 = np.sum(W * v ** 2, axis=1)
                    E = np.sqrt(e * s2x2 / sx2 ** 2)

                    ti = -params.cc.maxlag + params.mwcs.mwcs_wlen / 2. + count * (step_samples/goal_sampling_rate)
                    # print("Finished processing t_center=", ti, "s")
                    output.append((ti, M, E, MCOH))
                    count += 1
                    del fcur, fref, fcur2, fref2, result, cri
                    del X
                    del M, E, MCOH

                # Build xarray Dataset: dims (times, taxis, keys)
                t_centers = np.array([o[0] for o in output])
                M_arr    = np.stack([o[1] for o in output], axis=1)  # (n_times, n_taxis)
                E_arr    = np.stack([o[2] for o in output], axis=1)
                MCOH_arr = np.stack([o[3] for o in output], axis=1)

                times_coord = data.index.values
                ds_out = xr.Dataset(
                    {
                        "MWCS": xr.DataArray(
                            np.stack([M_arr, E_arr, MCOH_arr], axis=2),
                            dims=["times", "taxis", "keys"],
                            coords={
                                "times": times_coord,
                                "taxis": t_centers,
                                "keys":  ["M", "EM", "MCOH"],
                            },
                        )
                    }
                )
                xr_save_mwcs(params.output_folder, lineage_names, step.step_name,
                              station1, station2, components, mov_stack, taxis, ds_out)
                del data, output, ds_out

        massive_update_job(db, jobs, "D")
        # if not params.global_.hpc:
        #     for job in jobs:
        #         update_job(db, job.day, job.pair, 'DTT', 'T')
    logger.info('*** Finished: Compute MWCS ***')
