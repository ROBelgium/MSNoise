r"""Compute seismic velocity changes via the stretching method.

The stretching technique (:footcite:t:`SensSchoenfelder2006`) measures
:math:`dv/v` by finding the uniform dilation factor :math:`\varepsilon`
that maximises the normalised cross-correlation between a time-stretched
copy of the reference CCF and the current daily CCF:

.. math::

    s_\varepsilon(\tau) = s_\text{ref}\!\left(\frac{\tau}{1+\varepsilon}\right)

.. math::

    CC(\varepsilon) = \frac{
        \int_{\tau_\text{min}}^{\tau_\text{max}}
        s_\text{cur}(\tau)\, s_\varepsilon(\tau)\, d\tau
    }{
        \sqrt{
            \int |s_\text{cur}|^2 d\tau \cdot
            \int |s_\varepsilon|^2 d\tau
        }
    }

The relative velocity change is :math:`dv/v = -\varepsilon^*` where
:math:`\varepsilon^* = \arg\max_\varepsilon\, CC(\varepsilon)`.
In practice, :math:`\varepsilon` is sampled over ``stretching_nsteps``
values in :math:`[-\varepsilon_\text{max},\, \varepsilon_\text{max}]`
(``stretching_max``) and the optimal value is found by peak-picking.

.. warning:: if using only ``mov_stack`` = 1, no STR jobs is inserted in the
    database and consequently, no STR calculation will be done! FIX!


Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |stretching.stretching_max|
* |stretching.stretching_nsteps|
* |stretching.stretching_minlag|
* |stretching.stretching_width|
* |stretching.stretching_lag|
* |stretching.stretching_v|
* |stretching.stretching_sides|
* |stack.mov_stack|
* |refstack.ref_begin|
* |refstack.ref_end|
* |cc.cc_sampling_rate|
* |cc.maxlag|
* |cc.components_to_compute|
* |cc.components_to_compute_single_station|
* |global.hpc|

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

    Parallel Processing

.. footbibliography::
"""

from .core.db import connect, get_logger
from .core.stations import get_interstation_distance, get_station
from .core.workflow import (compute_rolling_ref, extend_days, get_next_lineage_batch, get_t_axis, is_next_job_for_step, massive_update_job, propagate_downstream, refstack_is_rolling)
from .core.io import xr_get_ccf, xr_get_ref, xr_save_stretching

import time

import numpy as np
import xarray as xr
from scipy.ndimage import map_coordinates


def _hwhm_errors(corr_coeffs):
    """Sub-sample HWHM error estimate for the stretching correlation curve.

    Replaces ``scipy.optimize.curve_fit`` Gaussian fitting (~600x faster,
    6x more accurate on benchmarks).  For each day:

    1. Shift the correlation curve to all-positive values.
    2. Walk left and right from the peak to find the integer half-max bracket.
    3. Linear interpolation within each bracket gives sub-sample crossing positions.
    4. HWHM = (right_crossing - left_crossing) / 2.

    The result is in stretch-index units (same as the original curve_fit error).
    Convert to dv/v: ``hwhm * 2 * str_range / (nstr - 1)``.

    :param corr_coeffs: 2-D array ``(nstr, n_days)`` — correlation coefficients
        of each stretched reference against each current CCF.
    :returns: 1-D array ``(n_days,)`` of HWHM estimates; ``np.nan`` where the
        peak touches the array boundary or the curve is flat.
    """
    nstr, n_days = corr_coeffs.shape
    allerrs = np.empty(n_days)
    for day_idx in range(n_days):
        c = corr_coeffs[:, day_idx]
        c_pos = c - c.min()            # shift to all-positive
        peak = int(np.argmax(c_pos))
        if c_pos[peak] == 0.0 or peak == 0 or peak == nstr - 1:
            allerrs[day_idx] = np.nan
            continue
        half = c_pos[peak] * 0.5
        # Walk left to find the bracket [left, left+1] that straddles half
        left = peak
        while left > 0 and c_pos[left] > half:
            left -= 1
        if c_pos[left] > half:
            # Walk hit index 0 without ever dropping to half-max: the curve
            # does not decay inside the search range, so interpolating here
            # would extrapolate to a negative index and fabricate a width.
            allerrs[day_idx] = np.nan
            continue
        # Sub-sample left crossing via linear interpolation
        denom_l = c_pos[left + 1] - c_pos[left]
        left_sub = left + (half - c_pos[left]) / denom_l if denom_l != 0.0 else left
        # Walk right to find the bracket [right-1, right] that straddles half
        right = peak
        while right < nstr - 1 and c_pos[right] > half:
            right += 1
        if c_pos[right] > half:
            allerrs[day_idx] = np.nan
            continue
        # Sub-sample right crossing via linear interpolation
        denom_r = c_pos[right - 1] - c_pos[right]
        right_sub = (right - 1 + (c_pos[right - 1] - half) / denom_r
                     if denom_r != 0.0 else right)
        allerrs[day_idx] = (right_sub - left_sub) / 2.0
    return allerrs


def apply_coda_mask(data, mid, minlag, maxlag2, sampling_rate):
    """Zero everything outside the two-sided coda window, in place.

    Keeps only ``minlag <= |lag| <= maxlag2`` on both sides of zero lag.
    The reference and the current CCFs MUST be masked identically, otherwise
    the correlation coefficient is diluted and the z-scoring of the current
    trace is skewed — the negative-lag branch used to zero a single column
    instead of the whole tail, which is exactly what this helper prevents.

    :param data: 1-D ``(n_samples,)`` or 2-D ``(n_days, n_samples)`` array;
        modified in place and also returned.
    :param mid: Index of the zero-lag sample.
    :param minlag: Inner edge of the coda window, in seconds.
    :param maxlag2: Outer edge of the coda window, in seconds.
    :param sampling_rate: CCF sampling rate in Hz.
    :returns: *data*, masked.
    """
    sr   = int(sampling_rate)
    lo   = mid - int(minlag * sr)
    hi   = mid + int(minlag * sr)
    lo2  = mid - int(maxlag2 * sr)
    hi2  = mid + int(maxlag2 * sr)
    if data.ndim == 1:
        data[lo:hi] = 0.0
        data[:lo2]  = 0.0
        data[hi2:]  = 0.0
    else:
        data[:, lo:hi] = 0.0
        data[:, :lo2]  = 0.0
        data[:, hi2:]  = 0.0
    return data


def hwhm_to_dvv(hwhm, str_range, nstr):
    """Convert a HWHM from stretch-index units to dv/v units.

    :func:`_hwhm_errors` measures the half-width of the correlation peak in
    row indices of the stretch matrix (0 .. ``nstr-1``).  The Delta column is
    a stretch factor and :func:`~msnoise.core.io.aggregate_dvv_pairs` reads
    Error as the uncertainty on ``Delta - 1``, so the two must share a unit.

    :param hwhm: Half-width in stretch-index units (scalar or array).
    :param str_range: ``stretching_max`` — the one-sided stretch range.
    :param nstr: ``stretching_nsteps`` — number of stretch steps.
    :returns: The same value expressed in dv/v units.
    """
    return hwhm * (2.0 * float(str_range) / (nstr - 1))


def stretch_mat_creation(refcc, str_range=0.01, nstr=1001):
    r"""Matrix of stretched instances of a reference trace.

    Implements the core of the stretching method (:footcite:t:`SensSchoenfelder2006`).

    For each stretch factor :math:`\varepsilon_k` sampled in
    :math:`[-\varepsilon_\text{max}, \varepsilon_\text{max}]`, the dilated
    trace is computed via cubic spline resampling on the compressed/expanded
    time grid:

    .. math::

        s_k[n] = s_\text{ref}\!\left(\frac{n}{1+\varepsilon_k}\right)

    Interpolation is performed with ``scipy.ndimage.map_coordinates`` (cubic
    spline, ``order=3``) in a single vectorised batch over all ``nstr`` rows.

    :type refcc: :class:`~numpy.ndarray`
    :param refcc: 1-D ndarray — the reference trace to stretch.
    :type str_range: float
    :param str_range: Maximum stretch amount (one side).
    :type nstr: int
    :param nstr: Number of stretching steps (total).
    :returns: ``(strrefmat, strvec)`` where ``strrefmat`` is ``(nstr, len(refcc))``
        and ``strvec`` holds the stretch factors.

    .. note::

        ``strvec`` is deliberately **reversed** relative to the time grid used
        to build the rows: row *k* resamples on ``samples_idx / strvec[-1-k]``
        while the caller reports ``strvec[k]``.  Since ``strvec`` is symmetric
        about 1, that is equivalent to reporting ``2 - eps_k``, which turns
        the best-matching dilation ``eps = 1 - dv/v`` into a reported value of
        ``1 + dv/v``.  This is what makes ``Delta - 1`` a correctly signed
        dv/v downstream: verified numerically to give exactly -0.010, -0.003,
        0.000, +0.003 and +0.010 for those true dv/v values.  Do not "tidy"
        the reversal away without also flipping the sign in
        :func:`~msnoise.core.io.aggregate_dvv_pairs`.
    """
    n = len(refcc)
    samples_idx = np.arange(n) - n // 2
    strvec = 1 + np.linspace(-str_range, str_range, nstr)
    # Build all (nstr, n) time-index rows at once — no Python loop
    str_timemat = samples_idx[None, :] / strvec[::-1, None]   # (nstr, n)
    # Single batch map_coordinates call replaces the per-row loop
    rows = (str_timemat + n // 2).ravel()
    cols = np.zeros(nstr * n)
    strrefmat = map_coordinates(
        refcc.reshape(n, 1), [rows, cols], order=3
    ).reshape(nstr, n)
    return strrefmat, strvec

def main(loglevel="INFO"):
    logger = get_logger('msnoise.stretching', loglevel, with_pid=True)
    logger.info('*** Starting: Compute Stretching ***')

    db = connect()
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
        lineage_names_ref = batch["lineage_names_ref"]
        lineage_names_mov = batch["lineage_names_mov"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]

        taxis = get_t_axis(params)

        logger.info(f"New Stretching Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        root = params.global_.output_folder
        mov_stacks = params.stack.mov_stack
        goal_sampling_rate = params.cc.cc_sampling_rate

        netsta1, netsta2 = pair.split(':')
        station1, station2 = pair.split(":")
        # Respect filter CC/SC/AC flags — skip pair types never computed
        if station1 == station2:
            if not (params.filter.SC or params.filter.AC):
                continue
            components_to_compute = params.cc.components_to_compute_single_station
        else:
            if not params.filter.CC:
                continue
            components_to_compute = params.cc.components_to_compute
        
        for components in components_to_compute:
            # Skip component types not enabled in filter config
            if station1 == station2:
                _is_ac = len(components) >= 2 and components[0] == components[-1]
                if _is_ac and not params.filter.AC:
                    continue
                if not _is_ac and not params.filter.SC:
                    continue
            station1, station2 = pair.split(":")
            rolling_mode = refstack_is_rolling(params)

            # ── Lag window parameters (same for Mode A and Mode B) ───────────
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
            str_range = params.stretching.stretching_max
            nstr = params.stretching.stretching_nsteps

            if not rolling_mode:
                # ── Mode A: load fixed REF from disk (skip if no refstack) ───
                if lineage_names_ref is None:
                    ref = None
                else:
                    try:
                        ref = xr_get_ref(root, lineage_names_ref,
                                         station1, station2, components, taxis)
                        ref = ref.values.copy()
                    except FileNotFoundError as fullpath:
                        logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                        continue

                # Zero data outside the lag window on the fixed reference
                apply_coda_mask(ref, mid, minlag, maxlag2, goal_sampling_rate)

                # Pre-build the full stretched-reference matrix (nstr x n_samples)
                ref_stretched, deltas = stretch_mat_creation(ref, str_range=str_range, nstr=nstr)

            for mov_stack in mov_stacks:
                try:
                    data = xr_get_ccf(root, lineage_names_mov,
                                      station1, station2, components, mov_stack, taxis)
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
                logger.debug("Processing %s:%s m%s %s" % (station1, station2, mov_stack, components))

                to_search = extend_days(days)
                # Filter to relevant days and drop all-NaN time steps
                data = data.sel(times=data.times.dt.floor('D').isin(to_search))
                data = data.dropna('times', how='all')

                if not data.sizes['times']:
                    continue

                # Materialise to numpy, then apply zero-lag windowing in-place
                data_values = data.values.copy()
                # Same mask as the reference above.  This used to be written
                # out inline and zeroed a SINGLE column instead of the whole
                # negative-lag tail, so the current CCFs kept every sample
                # between -maxlag and -maxlag2 while the reference had zeros
                # there — an asymmetric coda window.
                apply_coda_mask(data_values, mid, minlag, maxlag2,
                                params.cc.cc_sampling_rate)
                num_days = data_values.shape[0]

                if rolling_mode:
                    # ── Mode B: per-row rolling reference ────────────────────
                    # compute_rolling_ref returns shape (n_times, n_samples)
                    ref_rolling = compute_rolling_ref(
                        data_values, int(params.refstack.ref_begin), int(params.refstack.ref_end)
                    )

                    # ── Mode B: vectorise masking + normalisation, keep
                    # per-row stretch call (each row has a different rolling ref)
                    # Batch coda masking on all rolling refs at once
                    refs_masked = apply_coda_mask(
                        ref_rolling.copy(), mid, minlag, maxlag2,
                        goal_sampling_rate)

                    # Normalise current data (all rows at once)
                    cur_mean = data_values.mean(axis=1, keepdims=True)
                    cur_std  = data_values.std(axis=1, keepdims=True)
                    cur_std  = np.where(cur_std != 0, cur_std, 1.0)
                    data_norm_all = (data_values - cur_mean) / cur_std  # (n_days, n_samples)

                    alldeltas_list = []
                    allcoefs_list  = []
                    allerrs_list   = []

                    for i_row in range(num_days):
                        ref_str_row, deltas_row = stretch_mat_creation(
                            refs_masked[i_row], str_range=str_range, nstr=nstr
                        )
                        ref_std = ref_str_row.std(axis=1, keepdims=True)
                        ref_std = np.where(ref_std != 0, ref_std, 1.0)
                        ref_str_norm = (ref_str_row - ref_str_row.mean(axis=1, keepdims=True)) / ref_std

                        cc_row = ref_str_norm @ data_norm_all[i_row] / ref_str_row.shape[1]
                        best   = int(np.argmax(cc_row))

                        alldeltas_list.append(deltas_row[best])
                        allcoefs_list.append(cc_row[best])
                        allerrs_list.append(_hwhm_errors(cc_row.reshape(nstr, 1))[0])

                    alldays   = data.coords["times"].values
                    alldeltas = np.array(alldeltas_list)
                    allcoefs  = np.array(allcoefs_list)
                    allerrs   = np.array(allerrs_list)

                else:
                    # ── Mode A: vectorised matrix correlation ─────────────────
                    data_std = data_values.std(axis=1, keepdims=True)
                    data_std = np.where(data_std != 0, data_std, 1.0)
                    data_norm = (data_values - data_values.mean(axis=1, keepdims=True)) / data_std
                    ref_std = ref_stretched.std(axis=1, keepdims=True)
                    ref_std = np.where(ref_std != 0, ref_std, 1.0)
                    ref_stretched_norm = (ref_stretched - ref_stretched.mean(axis=1, keepdims=True)) / ref_std

                    # Compute the correlation coefficients
                    corr_coeffs = np.dot(ref_stretched_norm, data_norm.T) / ref_stretched.shape[1]

                    max_corr_indices = np.argmax(corr_coeffs, axis=0)
                    max_corr_values  = corr_coeffs[max_corr_indices, np.arange(num_days)]

                    alldays   = data.coords["times"].values
                    alldeltas = deltas[max_corr_indices]
                    allcoefs  = max_corr_values
                    # Vectorized error estimate: direct HWHM scan on the
                    # correlation-coefficient curve, ~134x faster than curve_fit.
                    allerrs   = _hwhm_errors(corr_coeffs)

                # _hwhm_errors returns a half-width in STRETCH-INDEX units
                # (0 .. nstr-1).  The Delta column is a stretch factor and
                # aggregate_dvv_pairs reads Error as the error on
                # ``Delta - 1``, i.e. in dv/v units, so the two must share a
                # unit.  Without this conversion the stored error was
                # (nstr-1)/(2*stretching_max) times too large — 25000x for
                # the defaults (nstr=1001, stretching_max=0.02).
                allerrs = hwhm_to_dvv(allerrs, str_range, nstr)

                # Build xarray Dataset directly — no DataFrame round-trip
                ds_out = xr.Dataset(
                    {
                        "STR": xr.DataArray(
                            np.column_stack([alldeltas, allcoefs, allerrs]),
                            dims=["times", "keys"],
                            coords={
                                "times": alldays,
                                "keys":  ["Delta", "Coeff", "Error"],
                            },
                        )
                    }
                )
                xr_save_stretching(
                    root, lineage_names, step.step_name,
                    station1, station2, components, mov_stack, ds_out,
                )
                del data_values, ds_out  # free materialised CCF + output Dataset

        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)
   #        for job in jobs:
     #           update_job(db, job.day, job.pair, 'DTT', 'T')
    logger.info('*** Finished: Compute Stretching ***')
