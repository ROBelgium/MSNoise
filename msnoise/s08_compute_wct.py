"""
Wavelet Coherence Transform (WCT) Computation
This script performs the computation of the Wavelet Coherence Transform (WCT), a tool used to analyze the correlation between two time series in the time-frequency domain. The script supports parallel processing and interacts with a database to manage job statuses.

Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* |dtt_minlag|
* |dtt_lag|
* |dtt_v|
* |dtt_width|
* |dtt_maxdt|
* |wct_mincoh|
* |wct_codacycles|
* |wct_ns|
* |wct_nt|
* |wct_vpo|
* |wct_nptsfreq|
* |wct_min_nonzero|
* |wct_norm|
* |hpc|

This process is job-based, so it is possible to run several instances in
parallel.


To run this step:

.. code-block:: sh

    $ msnoise cc dtt compute_wct

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dtt compute_wct

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

"""

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import scipy.signal
import pycwt as wavelet
from scipy.signal import convolve2d
from .api import *
import scipy
import scipy.fft as sf

def smoothCFS(cfs, scales, dt, ns, nt):
    """
    Smooth the continuous wavelet transform coefficients using a Fourier domain approach.
    Parameters
    ----------
    cfs : numpy.ndarray
        Continuous wavelet transform coefficients.
    scales : numpy.ndarray
        Scales used in the wavelet transform.
    dt : float
        Sampling interval.
    ns : int
        Smoothing parameter for the moving average filter.
    nt : float
        Smoothing parameter for the Gaussian filter.
    Returns
    -------
    numpy.ndarray
        Smoothed continuous wavelet transform coefficients.
    """
    N = cfs.shape[1]
    npad = sf.next_fast_len(N, real=True)
    omega = np.arange(1, np.fix(npad / 2) + 1, 1).tolist()
    omega = np.array(omega) * ((2 * np.pi) / npad)
    omega_save = -omega[int(np.fix((npad - 1) / 2)) - 1:0:-1]
    omega_2 = np.concatenate((0., omega), axis=None)
    omega_2 = np.concatenate((omega_2, omega_save), axis=None)
    omega = np.concatenate((omega_2, -omega[0]), axis=None)
    # Normalize scales by DT because we are not including DT in the angular frequencies here.
    # The smoothing is done by multiplication in the Fourier domain.
    normscales = scales / dt

    for kk in range(0, cfs.shape[0]):
        F = np.exp(-nt * (normscales[kk] ** 2) * omega ** 2)
        smooth = np.fft.ifft(F * np.fft.fft(cfs[kk - 1], npad))
        cfs[kk - 1] = smooth[0:N]
    # Convolve the coefficients with a moving average smoothing filter across scales.
    H = 1 / ns * np.ones((ns, 1))

    cfs = conv2(cfs, H)
    return cfs

def conv2(x, y, mode='same'):
    """
    Perform 2D convolution of matrices x and y
    """
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)

def get_wavelet_type(wavelet_type):
    """
    return a wavelet object based on the specified wavelet type and associated parameter
    """
    # Default parameters for each wavelet type
    default_params = {
        'Morlet': 6,
        'Paul': 4,
        'DOG': 2,
        'MexicanHat': 2  # MexicanHat inherits from DOG with m=2
    }

    wavelet_name = wavelet_type[0]

    # If a second argument is provided, use it; otherwise, use the default value
    if len(wavelet_type) == 2:
        param = float(wavelet_type[1])
    else:
        param = default_params[wavelet_name]

    # Get the corresponding wavelet object
    if wavelet_name == 'Morlet':
        return wavelet.Morlet(param)
    elif wavelet_name == 'Paul':
        return wavelet.Paul(param)
    elif wavelet_name == 'DOG':
        return wavelet.DOG(param)
    elif wavelet_name == 'MexicanHat':
        return wavelet.MexicanHat()  # Uses m=2, so no need for param
    else:
        raise logger.error(f"Unknown wavelet type: {wavelet_name}")



def xwt(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12, freqmin=0.1, freqmax=8.0, nptsfreq=100, wavelet_type=('Morlet',6.)):
    """
    Wavelet coherence transform (WCT) on two time series..
    The WCT finds regions in time frequency space where the two time
    series co-vary, but do not necessarily have high power.
    
    Modified from https://github.com/Qhig/cross-wavelet-transform
    Parameters
    ----------
    trace_ref, trace_current : numpy.ndarray, list
        Input signals.
    fs : float
        Sampling frequency.
    ns : smoothing parameter. 
        Default value is 3
    nt : smoothing parameter. 
        Default value is 0.25
    vpo : float,
        Spacing parameter between discrete scales. Default value is 12.
        Higher values will result in better scale resolution, but
        slower calculation and plot.
        
    freqmin : float,
        Smallest frequency
        Default value is 0.1 Hz
    freqmax : float,
        Highest frequency
        Default value is 8.0 Hz
    nptsfreq : int,
        Number of frequency points between freqmin and freqmax.
        Default value is 100 points
    wavelet_type: list,
        Wavelet type and associated parameter.
        Default Morlet wavelet with a central frequency w0 = 6
       
    Returns
        ----------
    WXamp : numpy.ndarray
        Amplitude of the cross-wavelet transform.
    WXspec : numpy.ndarray
        Complex cross-wavelet transform, representing both magnitude and phase information.
    WXangle : numpy.ndarray
        Phase angles of the cross-wavelet transform, indicating the phase relationship between the input signals.
    Wcoh : numpy.ndarray
        Wavelet coherence, representing the degree of correlation between the two signals in time-frequency space.
    WXdt : numpy.ndarray
        Time delay between the signals, estimated from the phase angles.
    freqs : numpy.ndarray
        Frequencies corresponding to the scales of the wavelet transform.
    coi : numpy.ndarray
        Cone of influence, representing the region of the wavelet spectrum where edge effects become significant.
    
    """
    
    mother = get_wavelet_type(wavelet_type) # mother wavelet class: Morlet, Paul, DOG, MexicanHat param 
    # nx represent the number of element in the trace_current array
    nx = np.size(trace_current)
    x_reference = np.transpose(trace_ref)
    x_current = np.transpose(trace_current)
    # Sampling interval
    dt = 1 / fs
    # Spacing between discrete scales, the default value is 1/12
    dj = 1 / vpo 
    # Number of scales less one, -1 refers to the default value which is J = (log2(N * dt / so)) / dj.
    J = -1
    # Smallest scale of the wavelet, default value is 2*dt
    s0 = 2 * dt  # Smallest scale of the wavelet, default value is 2*dt

    # Creation of the frequency vector that we will use in the continuous wavelet transform 
    freqlim = np.linspace(freqmax, freqmin, num=nptsfreq, endpoint=True, retstep=False, dtype=None, axis=0)

    # Calculation of the two wavelet transform independently
    # scales are calculated using the wavelet Fourier wavelength
    # fft : Normalized fast Fourier transform of the input trace
    # fftfreqs : Fourier frequencies for the calculated FFT spectrum.
    ###############################################################################################################
    try:
        cwt_reference, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(x_reference, dt, dj, s0, J, mother, freqs=freqlim)
        cwt_current, _, _, _, _, _ = wavelet.cwt(x_current, dt, dj, s0, J, mother, freqs=freqlim)
    except Exception as e:
        logger.error(f"Error in wavelet transform: {str(e)}")
        raise
    ###############################################################################################################

    scales = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1 / scales)

    # Apply smoothing
    power_ref = (invscales * abs(cwt_reference) ** 2).astype(complex)
    power_cur = (invscales * abs(cwt_current) ** 2).astype(complex)

    crossCFS = cwt_reference * np.conj(cwt_current)
    WXamp = abs(crossCFS)
    cross_spectrum = (invscales * crossCFS).astype(complex)

    # smoothCFS with complex arrays:
    cfs1 = smoothCFS(power_ref, scales, dt, ns, nt)
    cfs2 = smoothCFS(power_cur, scales, dt, ns, nt)
    crossCFS = smoothCFS(cross_spectrum, scales, dt, ns, nt)
    
    # Handle zeros to prevent divide-by-zero
    mask1 = cfs1 > 0
    mask2 = cfs2 > 0
    valid_mask = mask1 & mask2  # Areas where both have power
            
    # Initialize with NaNs
    WXspec = np.full_like(crossCFS, np.nan, dtype=complex)
    Wcoh = np.full_like(crossCFS, np.nan)
    
    # Division operations
    WXspec[valid_mask] = crossCFS[valid_mask] / (np.sqrt(cfs1[valid_mask]) * np.sqrt(cfs2[valid_mask]))
    Wcoh[valid_mask] = abs(crossCFS[valid_mask]) ** 2 / (cfs1[valid_mask] * cfs2[valid_mask])
    WXangle = np.angle(WXspec)

    # Clip coherence values to [0,1]
    Wcoh = np.clip(Wcoh, 0.0, 1.0)
       
    # Calculate time delay
    pp = 2 * np.pi * freqs
    pp2 = np.array([[kk] for kk in pp])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)


    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi



def main(loglevel="INFO", batch_size=None):
    """
    Main function to process WCT jobs using lineage-based approach
    """
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT ***')

    db = connect()


    db.close()

    while True:
        db = connect()
        if not is_next_job_for_step(db, step_category="wavelet"):
            db.close()
            break

        batch = get_next_lineage_batch(db, step_category="wavelet", group_by="pair_lineage",
                                       loglevel=loglevel, drop_current_step_name=False)
        db.close()

        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        lineage_names = batch["lineage_names"][:-1]
        lineage_names_mov = strip_refstack_from_lineage(lineage_names)
        lineage_str = batch["lineage_str"]
        step = batch["step"]
        root = params.output_folder

        logger.info(f"New WCT Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        taxis = get_t_axis(params)

        station1, station2 = pair.split(":")

        if station1 == station2:
            components_to_compute = params.components_to_compute_single_station
        else:
            components_to_compute = params.components_to_compute

        mov_stacks = params.mov_stack
        goal_sampling_rate = params.cc_sampling_rate

        # Filter frequency range from the wavelet step config
        freqmin = params.wct_freqmin
        freqmax = params.wct_freqmax

        # Wavelet computation parameters from the wavelet step config
        ns = getattr(params, 'wct_ns', 5)
        nt = getattr(params, 'wct_nt', 5)
        vpo = getattr(params, 'wct_vpo', 20)
        nptsfreq = getattr(params, 'wct_nptsfreq', 300)
        wct_norm = getattr(params, 'wct_norm', "Y")
        wavelet_type = eval(getattr(params, 'wavelet_type', '') or "('Morlet',6.)")

        logger.info(f"WCT params: freqmin={freqmin} freqmax={freqmax} ns={ns} nt={nt} vpo={vpo}")

        for component in components_to_compute:
            rolling_mode = refstack_is_rolling(params)
            if not rolling_mode:
                # Mode A: load fixed REF from disk
                try:
                    ref_data = xr_get_ref(root, lineage_names, station1, station2,
                                          component, None, taxis, ignore_network=True)
                    ref = ref_data.CCF.values
                    if wct_norm:
                        ori_waveform = ref / ref.max()
                    else:
                        ori_waveform = ref
                except FileNotFoundError as fp:
                    logger.error(f"FILE DOES NOT EXIST: {fp}, skipping")
                    continue
                except Exception as e:
                    logger.error(f"Error getting reference waveform: {str(e)}")
                    continue

                if not len(ref):
                    continue
            # Mode B: ori_waveform set per time-step below after data is loaded

            for mov_stack in mov_stacks:
                WXamp_list = []
                WXcoh_list = []
                WXdt_list = []
                dates_list = []

                try:
                    data = xr_get_ccf(root, lineage_names_mov, station1, station2,
                                      component, None, mov_stack, taxis)
                except FileNotFoundError as fp:
                    logger.error(f"FILE DOES NOT EXIST: {fp}, skipping")
                    continue

                to_search = pd.to_datetime(days)
                to_search = to_search.append(
                    pd.DatetimeIndex([to_search[-1] + pd.Timedelta("1d")])
                )
                data = data[data.index.floor('d').isin(to_search)]
                data = data.dropna()

                if rolling_mode:
                    ref_rolling = compute_rolling_ref(
                        data, int(params.ref_begin), int(params.ref_end)
                    )

                for _i_row, (date, row) in enumerate(data.iterrows()):
                    waveform = row.values
                    if wct_norm:
                        new_waveform = waveform / waveform.max()
                    else:
                        new_waveform = waveform

                    if rolling_mode:
                        _rr = ref_rolling[_i_row]
                        ori_waveform = _rr / _rr.max() if wct_norm else _rr

                    try:
                        WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                            ori_waveform, new_waveform, goal_sampling_rate,
                            int(ns), float(nt), int(vpo),
                            freqmin, freqmax, int(nptsfreq), wavelet_type
                        )
                        WXamp_list.append(WXamp)
                        WXcoh_list.append(Wcoh)
                        WXdt_list.append(WXdt)
                        dates_list.append(date)
                    except Exception as e:
                        logger.error(f"Error in WCT for {date}: {str(e)}")
                        continue

                if dates_list:
                    try:
                        xr_save_wct2(root, lineage_names, step.step_name,
                                     station1, station2, component, mov_stack,
                                     taxis, freqs, WXamp_list, WXcoh_list,
                                     WXdt_list, dates_list)
                    except Exception as e:
                        logger.error(f"Error saving WCT: {str(e)}")

        db = connect()
        massive_update_job(db, jobs, "D")
        db.close()

    logger.info('*** Finished: Compute WCT ***')
