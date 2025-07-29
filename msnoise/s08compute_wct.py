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

    $ msnoise cc dvv compute_wct

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dvv compute_wct

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

"""

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import scipy.optimize
import scipy.signal
import pycwt as wavelet
from scipy.signal import convolve2d
from obspy.signal.regression import linear_regression
from .api import *
import logbook
import scipy
import scipy.fft as sf
import warnings
from scipy.optimize import OptimizeWarning
import re

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

def compute_wct_dvv(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20, mincoh=0.5, maxdt=0.2, 
            min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
    """
    Compute the dv/v values and associated errors from the wavelet transform results.
    Parameters
    ----------
    freqs : numpy.ndarray
        Frequency values corresponding to the wavelet transform.
    tvec : numpy.ndarray
        Time vector.
    WXamp : numpy.ndarray
        Amplitude of the cross-wavelet transform.
    Wcoh : numpy.ndarray
        Wavelet coherence.
    delta_t : numpy.ndarray
        Time delays between signals.
    lag_min : int, optional
        Minimum lag in seconds. Default is 5.
    coda_cycles : int, optional
        Number of coda cycles to consider. Default is 20.
    mincoh : float, optional
        Minimum coherence value for weighting. Default is 0.5.
    maxdt : float, optional
        Maximum time delay for weighting. Default is 0.2.
    min_nonzero : float, optional
        Minimum percentage of non-zero weights required for valid estimation. Default is 0.25.
    freqmin : float, optional
        Minimum frequency for calculation. Default is 0.1 Hz.
    freqmax : float, optional
        Maximum frequency for calculation. Default is 2.0 Hz.
    Returns
    -------
    tuple
        dvv values (percentage), errors (percentage), and weighting function used.
    """   
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))  # Filter frequencies within the specified range
    dvv, err = np.zeros(len(inx[0])), np.zeros(len(inx[0])) # Initialize dvv and err arrays

    # Weighting function based on WXamp
    weight_func = np.log(np.abs(WXamp)) / np.log(np.abs(WXamp)).max()
    zero_idx = np.where((Wcoh < mincoh) | (delta_t > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / weight_func.max()
    wf[zero_idx] = 0

    # Track frequencies with estimation issues
    problematic_freqs = []
                
    # Loop through frequency indices for linear regression
    for ii, ifreq in enumerate(inx[0]):
        period = 1.0 / freqs[ifreq]
        lag_max = lag_min + (period * coda_cycles)

        # Coda selection
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0]

        if len(tvec) > 2:
            if not np.any(delta_t[ifreq]):
                continue

            delta_t[ifreq][tindex] = np.nan_to_num(delta_t[ifreq][tindex])
            w = wf[ifreq]  # Weighting function for the specific frequency
            non_zero_mask = w > 0
            nzc_perc = np.count_nonzero(non_zero_mask[tindex]) / len(tindex)
            if nzc_perc >= min_nonzero:
                with warnings.catch_warnings(record=True) as w_catcher:
                    warnings.simplefilter("always", OptimizeWarning)
                    m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex], w[tindex], intercept_origin=True)
                    if any(issubclass(warning.category, OptimizeWarning) for warning in w_catcher):
                        problematic_freqs.append(freqs[ifreq])
                        
                dvv[ii], err[ii] = -m, em
            else:
                dvv[ii], err[ii] = np.nan, np.nan
        else:
            logger.debug('Not enough points to estimate dv/v for WCT')
            
    if problematic_freqs:
        logger.warning(f"Covariance estimation issues between {min(problematic_freqs):.2f}-{max(problematic_freqs):.2f} Hz: Modify min_nonzero (current: {min_nonzero}), mincoh (current: {mincoh}), maxdt (current: {maxdt}) and/or coda_cycles (current: {coda_cycles})")
    
    return dvv * 100, err * 100, wf

def get_avgcoh(freqs, tvec, wcoh, freqmin, freqmax, lag_min=5, coda_cycles=20):
    """
    Calculate the average wavelet coherence over a specified frequency range and time lags.

    :param freqs: A numpy array that represents frequency values.
    :type freqs: numpy.ndarray
    :param tvec: A time vector represented as a numpy array.
    :type tvec: numpy.ndarray
    :param wcoh: The wavelet coherence array, represented as a numpy array.
    :type wcoh: numpy.ndarray
    :param freqmin: The minimum frequency for coherence calculation, represented as a floating-point number.
    :type freqmin: float
    :param freqmax: The maximum frequency for coherence calculation, represented as a floating-point number.
    :type freqmax: float
    :param lag_min: The minimum lag in seconds for coherence calculation. This is optional and it defaults to 5.
    :type lag_min: int, optional
    :param coda_cycles: The number of coda cycles to consider. This is optional and it defaults to 20.
    :type coda_cycles: int, optional
    :returns: A numpy array of average coherence values computed over the specified frequency range and time lags.
    :rtype: numpy.ndarray
    """
    inx = np.where((freqs>=freqmin) & (freqs<=freqmax)) 
    coh = np.zeros(inx[0].shape) # Create empty vector for coherence

    for ii, ifreq in enumerate(inx[0]): # Loop through frequencies index     
        period = 1.0/freqs[ifreq]
        lag_max = lag_min + (period*coda_cycles) 
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0] # Index of the coda

        if len(tvec)>2: # check time vector size
            if not np.any(wcoh[ifreq]) or wcoh[ifreq][tindex].size == 0:
                coh[ii] = np.nan
                continue
            c = np.nanmean(np.abs(wcoh[ifreq][tindex]))
            coh[ii] = c

        else:
            logger.debug('Not enough points to compute average coherence') #not sure why it would ever get here, but just in case.
            coh[ii] = np.nan

    return coh

def save_day_wct_results(station1, station2, date, component, filterid, wctid, dttid, mov_stack, taxis, dvv, err, coh, freqs, inx):
    """
    Save the WCT results for a single day to a file in a hierarchical folder structure:
    DVV/WCT/WCT/f{filterid}/wct{wctid}/{mov_stack}/{component}/{station_pair}/{date}.npz
    
    Parameters:
    ----------
    station1, station2 : str
        Names of the stations
    date : datetime
        The date of the measurement
    component : str
        Component code
    filterid : int
        Filter ID
    mov_stack : int or str
        Moving stack value
    taxis : numpy.ndarray
        Time axis
    dvv, err, coh : numpy.ndarray
        The dv/v values, errors, and coherence values
    freqs : numpy.ndarray
        Frequency values
    inx : tuple
        Indices of valid frequencies
        
    Returns:
    -------
    str
        Path to the saved file
    """
    # Format filterid with leading zeros (e.g., 1 -> "01")
    filter_str = f"{int(filterid):02d}"
    wct_str = f"wct{int(wctid):02d}"
    dtt_str = f"dtt{int(dttid):02d}"

    # Format mov_stack (could be string like '10d_1d' or int)
    if isinstance(mov_stack, (list, tuple)):
        mov_stack_str = '_'.join(str(m) for m in mov_stack)
    else:
        mov_stack_str = str(mov_stack)
    
    # Station pair string
    pair_str = f"{station1}_{station2}"
    
    # Create directory structure
    #DVV/WCT/WCT/f{filterid}/wct{wctid}/{mov_stack}/{component}/{station_pair}/
    dir_path = os.path.join('DVV', 'WCT','WCT', filter_str, wct_str, dtt_str, mov_stack_str, component, pair_str)
    os.makedirs(dir_path, exist_ok=True)
    
    # Filename is just the date: YYYY-MM-DD.npz
    date_str = date.strftime('%Y-%m-%d')
    filename = f"{date_str}.npz"
    filepath = os.path.join(dir_path, filename)
    
    # Create DataFrames
    freqs_subset = freqs[inx]
    
    # Save to numpy compressed format
    np.savez_compressed(
        filepath,
        date=date,
        freqs=freqs_subset,
        dvv=dvv,
        err=err,
        coh=coh,
        taxis=taxis,
        station1=station1,
        station2=station2,
        component=component,
        filterid=filterid,
        mov_stack=mov_stack
    )
    
    logger.debug(f"Saved WCT results to {filepath}")
    return filepath

def get_all_station_distances(db):
    """
    Pre-calculate distances for all station pairs using MSNoise's get_interstation_distance API.
    
    Parameters:
    -----------
    db : database session
        MSNoise database session
        
    Returns:
    --------
    dict
        Dictionary with station pair keys and distance values in km
    """
    from itertools import combinations
    
    logger.info("Pre-calculating all station distances...")
    
    # Get all stations from database
    stations = get_stations(db, all=True)
    
    # Filter stations that have coordinates
    valid_stations = []
    for station in stations:
        if station.X is not None and station.Y is not None:
            valid_stations.append(station)
        else:
            logger.debug(f"Station {station.net}.{station.sta} has no coordinates, skipping")
    
    logger.info(f"Found {len(valid_stations)} stations with coordinates")
    
    # Calculate distances for all station pairs
    distances = {}
    
    for station1, station2 in combinations(valid_stations, 2):
        try:
            # Calculate distance using MSNoise's built-in function
            distance = get_interstation_distance(station1, station2, coordinates="DEG")
            
            # Create comprehensive key formats for this pair
            # Basic formats
            basic_keys = [
                f"{station1.net}.{station1.sta}:{station2.net}.{station2.sta}",
                f"{station2.net}.{station2.sta}:{station1.net}.{station1.sta}",
                f"{station1.net}.{station1.sta}_{station2.net}.{station2.sta}",
                f"{station2.net}.{station2.sta}_{station1.net}.{station1.sta}",
                # Station codes only
                f"{station1.sta}:{station2.sta}",
                f"{station2.sta}:{station1.sta}",
                f"{station1.sta}_{station2.sta}",
                f"{station2.sta}_{station1.sta}",
            ]
            
            # Add formats with location codes (common in MSNoise jobs)
            location_formats = [
                f"{station1.net}.{station1.sta}.--:{station2.net}.{station2.sta}.--",
                f"{station2.net}.{station2.sta}.--:{station1.net}.{station1.sta}.--",
                f"{station1.net}.{station1.sta}.---{station2.net}.{station2.sta}.--",
                f"{station2.net}.{station2.sta}.---{station1.net}.{station1.sta}.--",
                f"{station1.net}.{station1.sta}.--_{station2.net}.{station2.sta}.--",
                f"{station2.net}.{station2.sta}.--_{station1.net}.{station1.sta}.--",
            ]
            
            # Combine all key formats
            all_keys = basic_keys + location_formats
            
            # Store distance with all possible key formats
            for key in all_keys:
                distances[key] = distance
                
            logger.debug(f"Distance {station1.net}.{station1.sta} - {station2.net}.{station2.sta}: {distance:.2f} km")
            
        except Exception as e:
            logger.error(f"Error calculating distance between {station1.net}.{station1.sta} and {station2.net}.{station2.sta}: {e}")
            continue
    
    logger.info(f"Pre-calculated {len(distances)//(len(basic_keys) + len(location_formats))} unique station pair distances")
    logger.info(f"Total cache entries: {len(distances)}")
    
    # Debug: show a few cache keys
    sample_keys = list(distances.keys())[:10]
    logger.debug(f"Sample cache keys: {sample_keys}")
    
    return distances



def get_distance_from_cache(distances_cache, station1, station2):
    """
    Get distance from pre-calculated cache.
    
    Parameters:
    -----------
    distances_cache : dict
        Pre-calculated distances dictionary
    station1, station2 : str
        Station names (can include location codes like 8K.ASK.-- or just 8K.ASK)
        
    Returns:
    --------
    float or None
        Distance in kilometers, or None if not found or if it's an autocorrelation
    """
    def clean_station_name(station_name):
        """Clean station name by removing location/channel codes"""
        # Remove location codes (anything after the second dot)
        parts = station_name.split('.')
        if len(parts) >= 2:
            # Keep only NET.STA, remove location/channel codes
            return f"{parts[0]}.{parts[1]}"
        return station_name
    
    # Clean the station names
    clean_sta1 = clean_station_name(station1)
    clean_sta2 = clean_station_name(station2)
    
    logger.debug(f"Looking for distance between {station1} ({clean_sta1}) and {station2} ({clean_sta2})")
    
    # Handle autocorrelations (same station)
    if clean_sta1 == clean_sta2:
        logger.debug(f"Autocorrelation detected for {clean_sta1}, distance = 0.0 km")
        return 0.0
    
    # For cross-correlations, try different key formats with both original and cleaned names
    possible_keys = [
        # Try with cleaned names first
        f"{clean_sta1}:{clean_sta2}",
        f"{clean_sta2}:{clean_sta1}",
        f"{clean_sta1}_{clean_sta2}",
        f"{clean_sta2}_{clean_sta1}",
        # Try with original names
        f"{station1}:{station2}",
        f"{station2}:{station1}",
        f"{station1}_{station2}",
        f"{station2}_{station1}",
        # Try extracting just station codes (without network)
        f"{clean_sta1.split('.')[1]}:{clean_sta2.split('.')[1]}" if '.' in clean_sta1 and '.' in clean_sta2 else None,
        f"{clean_sta2.split('.')[1]}:{clean_sta1.split('.')[1]}" if '.' in clean_sta1 and '.' in clean_sta2 else None,
    ]
    
    # Filter out None values
    possible_keys = [k for k in possible_keys if k is not None]
    
    logger.debug(f"Trying keys: {possible_keys[:4]}...")  # Log first few keys
    
    for key in possible_keys:
        if key in distances_cache:
            logger.debug(f"Found cached distance for {station1}-{station2} using key '{key}': {distances_cache[key]:.2f} km")
            return distances_cache[key]
    
    logger.debug(f"No cached distance found for {station1}-{station2} after trying {len(possible_keys)} key formats")
    return None
      
def process_wct_job(pair, day, params, taxis, filters, wct_params_dict, distances_cache= None):  
    """
    Process a single WCT job without database access.
    
    Parameters:
    ----------
    job_data : dict
        Dictionary with job information (ref, pair, day, flag)
    params : dict-like
        Configuration parameters
    taxis : numpy.ndarray
        Time axis
    filters : list
        List of filters to use
    distances_cache : dict, optional
        Pre-calculated station distances dictionary
    Returns:
    -------
    bool
        True if job was successfully processed, False otherwise
    """
    # Parse pair and day
    station1, station2 = pair.split(":")
    
    dist = 0.0  # Default for autocorrelations
    if distances_cache:
        cached_dist = get_distance_from_cache(distances_cache, station1, station2)
        if cached_dist is not None:
            dist = cached_dist
        else:
            logger.warning(f"Distance not found in cache for {pair}, using default 0.0")
    else:
        # Fallback to original calculation if no cache available
        try:
            n1, s1, l1 = station1.split(".")
            n2, s2, l2 = station2.split(".")
            s1 = "%s_%s" % (n1, s1)
            s2 = "%s_%s" % (n2, s2)
            if s1 != s2:
                # This would require database access - not ideal in this approach
                logger.warning(f"No distance cache available and cannot calculate distance for {pair}")
        except:
            logger.warning(f"Could not parse station names for distance calculation: {pair}")
    
    logger.debug(f"Using distance {dist:.2f} km for pair {pair}")

    # Determine components to compute
    if station1 == station2:
        components_to_compute = params.components_to_compute_single_station
    else:
        components_to_compute = params.components_to_compute
        
    logger.info(f"Processing Pair: {pair}, Day: {day} for {len(filters)} filters, {len(components_to_compute)} components and {len(params.mov_stack)} moving windows")
        
    # Extract parameters
    #ns = params.wct_ns
    #nt = params.wct_nt 
    #vpo = params.wct_vpo 
    #nptsfreq = params.wct_nptsfreq
    #coda_cycles = params.dtt_codacycles 
    #min_nonzero = params.dvv_min_nonzero
    #wct_norm = params.wct_norm
    #wavelet_type = params.wavelet_type
    
    mov_stacks = params.mov_stack
    goal_sampling_rate = params.cc_sampling_rate
    #maxdt = params.dtt_maxdt
    #mincoh = params.dtt_mincoh

    # Convert date string to datetime
    try:
        date = pd.to_datetime(day)
    except:
        logger.error(f"Failed to parse date: {day}")
        return False
    
    operations_succeeded = 0  # Track if all operations for this job succeed

    # Process each filter
    for f in filters:
        filterid = int(f['ref'])
        freqmin = f['low']
        freqmax = f['high']

        #### due to strcucure change####
        if filterid not in wct_params_dict:
            logger.debug(f"No WCT parameters found for filter {filterid}, skipping")
            continue
        
        # Process each WCT parameter set for this filter
        for wctid, dtt_param_list in wct_params_dict[filterid].items():
            
            # Process each DTT parameter set
            for dtt_params in dtt_param_list:
                if not dtt_params['used']:
                    continue
                
                dttid = dtt_params['ref']
                
                # Extract WCT parameters (with defaults if not in DTT params)
                ns = getattr(params, 'wct_ns', 5)
                nt = getattr(params, 'wct_nt', 5) 
                vpo = getattr(params, 'wct_vpo', 20)
                nptsfreq = getattr(params, 'wct_nptsfreq', 300)
                wct_norm = getattr(params, 'wct_norm', "Y")
                wavelet_type = getattr(params, 'wavelet_type', "('Morlet',6.)")
                
                # Use DTT parameters for processing
                coda_cycles = dtt_params['wct_codacycles']
                min_nonzero = dtt_params['wct_min_nonzero'] / 100.0  # Convert percentage to fraction
                maxdt = dtt_params['wct_maxdt']
                mincoh = dtt_params['wct_mincoh']
                
                # Calculate lag parameters based on DTT settings
                if dtt_params['wct_lag'] == "static":
                    lag_min = dtt_params['wct_minlag']
                else:
                    if dtt_params['wct_v'] > 0:
                        lag_min = dist / dtt_params['wct_v']
                    else:
                        logger.warning(f"Invalid velocity {dtt_params['wct_v']} for dynamic lag, using static")
                        lag_min = dtt_params['wct_minlag']
                
                lag_max = lag_min + dtt_params['wct_width']
                
                # Use frequency bounds from DTT parameters
                processing_freqmin = dtt_params['wct_dtt_freqmin']
                processing_freqmax = dtt_params['wct_dtt_freqmax']
                
                logger.debug(f"Using DTT params - lag_min: {lag_min}, lag_max: {lag_max}, "
                           f"freq_range: {processing_freqmin}-{processing_freqmax} Hz, "
                           f"sides: {dtt_params['wct_sides']}")
        #### end: due to strcucure change####
                
                # Process each component
                for component in components_to_compute:
                    try:
                        # Get reference waveform
                        ref = xr_get_ref(station1, station2, component, filterid, taxis, ignore_network=True)
                        ref = ref.CCF.values
                        if wct_norm:
                            ori_waveform = (ref/ref.max()) 
                        else:
                            ori_waveform = ref
                    except FileNotFoundError as fullpath:
                        logger.error(f"FILE DOES NOT EXIST: {fullpath}, skipping")
                        continue
                    except Exception as e:
                        logger.error(f"Error getting reference waveform: {str(e)}")
                        continue
                    
                    if not len(ref):
                        logger.warning(f"Empty reference waveform for {station1}:{station2} {component} f{filterid}")
                        continue
                    
                    # Process each moving stack
                    for mov_stack in mov_stacks:
                        try:
                            # Get the data for this day and moving stack
                            data = xr_get_ccf(station1, station2, component, filterid, mov_stack, taxis)
                            
                            # Filter data to get only this specific day
                            day_data = data[data.index == date]
                            
                            if day_data.empty:
                                logger.warning(f"No data for {date} in {station1}:{station2} {component} f{filterid} m{mov_stack}")
                                continue
                            
                            logger.debug(f"Processing {station1}:{station2} {component} f{filterid} m{mov_stack} - {date}")
                            
                            # Get the waveform for this day
                            waveform = day_data.iloc[0].values
                            if wct_norm:
                                new_waveform = (waveform/waveform.max()) 
                            else:
                                new_waveform = waveform
                            
                            # Calculate wavelet coherence transform
                            WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                                ori_waveform, new_waveform, goal_sampling_rate, 
                                int(ns), float(nt), int(vpo), 
                                freqmin, freqmax, int(nptsfreq), wavelet_type, maxdt
                            )
                            
                            # Compute dv/v values
                            dvv, err, wf = compute_wct_dvv(
                                freqs, taxis, WXamp, Wcoh, WXdt, 
                                lag_min=int(lag_min), coda_cycles=coda_cycles, 
                                mincoh=mincoh, maxdt=maxdt, min_nonzero=min_nonzero, 
                                freqmin=freqmin, freqmax=freqmax
                            )
                            
                            # Calculate average coherence
                            coh = get_avgcoh(
                                freqs, taxis, Wcoh, freqmin, freqmax, 
                                lag_min=int(lag_min), coda_cycles=coda_cycles
                            )
                            
                            # Get frequency indices
                            inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
        
                            # Save results for this day
                            save_day_wct_results(
                                station1, station2, date, component, 
                                filterid, wctid, dttid, mov_stack, taxis, 
                                dvv, err, coh, freqs, inx
                            )
                            
                            operations_succeeded += 1
                        
                        except Exception as e:
                            logger.error(f"Error processing {component} f{filterid} m{mov_stack}: {str(e)}")
                            continue
    
    # Determine overall success status
    if operations_succeeded == 0:
        logger.error(f"Job failed completely: {pair} - {day}.")
        return False
    else:
        logger.info(f"Job processing complete: {pair} - {day}. {operations_succeeded} operations succeeded.")
        return True

def claim_available_jobs(db, batch_size=5, max_attempts=3):
    """
    Claims a batch of available jobs - SIMPLE VERSION
    """
    import time
    from datetime import datetime
        
    for attempt in range(max_attempts):
        try:
            # Find available jobs
            available_jobs = db.query(Job).filter(
                Job.flag == 'T',
                Job.jobtype == 'WCT',
                Job.day != 'REF'
            ).limit(batch_size).all()
            
            if not available_jobs:
                if attempt < max_attempts - 1:
                    time.sleep(1 * (attempt + 1))
                    continue
                else:
                    return []
            
            # Get job IDs
            job_ids = [job.ref for job in available_jobs]
            
            # Use ORM update - much simpler and more reliable
            updated_count = db.query(Job).filter(
                Job.ref.in_(job_ids),
                Job.flag == 'T'
            ).update(
                {Job.flag: 'I', Job.lastmod: datetime.utcnow()},
                synchronize_session=False
            )
            
            db.commit()
            
            if updated_count > 0:
                # Return the updated jobs
                return available_jobs[:updated_count]
            else:
                if attempt < max_attempts - 1:
                    time.sleep(1 * (attempt + 1))
                    continue
                else:
                    return []
                
        except Exception as e:
            db.rollback()
            logger.error(f"Error claiming jobs: {str(e)}")
            if attempt < max_attempts - 1:
                time.sleep(1)
            else:
                raise
    
    return []

def main(loglevel="INFO", batch_size=5):
    """
    Main function to process WCT jobs using day job distribution approach
    """
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT ***')
    
    # Add small random delay to prevent all processes starting at exactly the same time
    time.sleep(np.random.random() * 2)
    
    try:
        # Initial database connection to get configuration
        db = connect()
        
        # Ensure no lingering transactions
        if hasattr(db, 'in_transaction') and db.in_transaction():
            db.rollback()
            
        params = get_params(db)
        taxis = get_t_axis(db)
        filters = get_filters(db, all=False)
        wct_dtt_params = get_dvv_wct_dtt_jobs(db, all=False)
        
        # Convert filters to dictionaries before closing the session
        filter_dicts = []
        for f in filters:
            filter_str = str(f)
            match = re.match(r'Filter (\d+) \(([\d.]+)-([\d.]+) Hz\)', filter_str)
    
            filter_dicts.append({
                'ref': match.group(1),
                'low': float(match.group(2)),
                'high': float(match.group(3))
            })
        
        wct_params = {}
        for filter_ref, wct_list in wct_dtt_params.items():
            filterid = int(filter_ref)
            wct_params[filterid] = {}
            
            for wct_ref, dtt_list in wct_list.items():
                wctid = int(wct_ref)
                wct_params[filterid][wctid] = []
                
                for dtt_params in dtt_list:
                    # Convert ORM object to dictionary for easier handling
                    dtt_dict = {
                        'ref': dtt_params.ref,
                        'wct_dtt_freqmin': dtt_params.wct_dtt_freqmin,
                        'wct_dtt_freqmax': dtt_params.wct_dtt_freqmax,
                        'wct_minlag': dtt_params.wct_minlag,
                        'wct_width': dtt_params.wct_width,
                        'wct_lag': dtt_params.wct_lag,
                        'wct_v': dtt_params.wct_v,
                        'wct_sides': dtt_params.wct_sides,
                        'wct_mincoh': dtt_params.wct_mincoh,
                        'wct_maxdt': dtt_params.wct_maxdt,
                        'wct_codacycles': dtt_params.wct_codacycles,
                        'wct_min_nonzero': dtt_params.wct_min_nonzero,
                        'used': dtt_params.used
                    }
                    wct_params[filterid][wctid].append(dtt_dict)
                    
        #params.coda_velocity = 2.0
        #params.coda_safety_factor = 1.3
        
        # Pre-calculate all station distances
        try:
            distances_cache = get_all_station_distances(db)
            logger.info(f"Successfully cached distances for station pairs")
        except Exception as e:
            logger.error(f"Error pre-calculating distances: {e}")
            distances_cache = {}
            
        db.close()
        
        # Process jobs in batches to avoid too many database connections
        jobs_processed = 0
        max_jobs_per_process = 1000  # Safety limit
        consecutive_empty_batches = 0
        max_consecutive_empty = 3  # Exit after this many empty batches
        
        while jobs_processed < max_jobs_per_process:
            # Connect to the database and get a batch of jobs
            db = connect()
            
            # Ensure no lingering transactions
            if hasattr(db, 'in_transaction') and db.in_transaction():
                db.rollback()
            
            # Check if there are any jobs left
            job_exists = is_dtt_next_job(db, flag='T', jobtype='WCT')
            if not job_exists:
                logger.info("No more WCT jobs to process")
                db.close()
                break
                
            # Try to claim a batch of jobs
            jobs = claim_available_jobs(db, batch_size)
            
            if not jobs:
                consecutive_empty_batches += 1
                if consecutive_empty_batches >= max_consecutive_empty:
                    logger.info(f"Exiting after {consecutive_empty_batches} empty batches")
                    break
                time.sleep(2)  # Wait before next attempt
                continue
            else:
                consecutive_empty_batches = 0  # Reset counter when we get jobs
                
            # Process each job in the batch
            for job in jobs:
                pair = job.pair
                day = job.day

                # Process this job
                success = process_wct_job(pair, day, params, taxis, filter_dicts, wct_params, distances_cache)
                                
                # Update job status
                db = connect()
                if success:
                    update_job(db, day, pair, 'WCT', 'D')
                else:
                    update_job(db, day, pair, 'WCT', 'E')  # Mark as error
                db.close()
                
                jobs_processed += 1 
                
            # Small delay between batches
            time.sleep(0.5)
            
        logger.info(f"Process completed {jobs_processed} jobs")
                        
    except Exception as e:
        logger.error(f"Error in main processing: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
    finally:
        # Final cleanup
        try:
            if 'db' in locals() and db is not None:
                # Ensure we don't leave transactions open
                if hasattr(db, 'in_transaction') and db.in_transaction():
                    db.rollback()
                db.close()
        except:
            pass

    logger.info(f'*** Finished: Compute WCT ***')

if __name__ == "__main__":
    main()
