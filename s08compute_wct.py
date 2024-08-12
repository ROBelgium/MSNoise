"""
WCT


Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* |hpc| | 

This process is job-based, so it is possible to run several instances in
parallel.

Once done, each job is marked "D"one in the database and, unless ``hpc`` is 
``Y``, DTT jobs are inserted/updated in the database.

To run this step:

.. code-block:: sh

    $ msnoise cc dvv compute_wct

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dvv compute_wct

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.
    Parallel Processing
"""

import os
import time
import numpy as np
import pandas as pd
import xarray as xr
import scipy.optimize
import scipy.signal
import pycwt as wavelet
from scipy.fft import next_fast_len
from scipy.signal import convolve2d
from scipy.stats import scoreatpercentile
from obspy.signal.regression import linear_regression
from .api import *

import logbook
import scipy
import scipy.fft as sf

def get_avgcoh(freqs, tvec, wcoh, freqmin, freqmax, lag_min=5, coda_cycles=20):

    inx = np.where((freqs>=freqmin) & (freqs<=freqmax)) 
    coh = np.zeros(inx[0].shape) # Create empty vector for coherence

    for ii, ifreq in enumerate(inx[0]): # Loop through frequencies index     
        period = 1.0/freqs[ifreq]
        lag_max = lag_min + (period*coda_cycles) 
        tindex = np.where(((tvec >= -lag_max) & (tvec <= -lag_min)) | ((tvec >= lag_min) & (tvec <= lag_max)))[0] # Index of the coda

        if len(tvec)>2: # check time vector size
            if not np.any(wcoh[ifreq]): # check non-empty dt array
                continue
            c = np.nanmean(wcoh[ifreq][tindex])
            coh[ii] = c

        else:
            print('not enough points to compute average coherence') #not sure why it would ever get here, but just in case.
            coh[ii] = np.nan

    return coh

def smoothCFS(cfs, scales, dt, ns, nt):
    """
    Smoothing function
    """
    N = cfs.shape[1]
    npad = int(2 ** nextpow2(N))
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


## nextpow2 function
# Returns the exponents p for the smallest powers of two that satisfy the relation  : 2**p >= abs(x)
def nextpow2(x):
    res = np.ceil(np.log2(x))
    return res.astype('int')

## conv2 function
# Returns the two-dimensional convolution of matrices x and y
def conv2(x, y, mode='same'):
    return np.rot90(convolve2d(np.rot90(x, 2), np.rot90(y, 2), mode=mode), 2)


def get_dvv(freqs, tvec, WXamp, Wcoh, delta_t, lag_min=5, coda_cycles=20, mincoh=0.5, maxdt=0.2, 
            min_nonzero=0.25, freqmin=0.1, freqmax=2.0):
   
    inx = np.where((freqs >= freqmin) & (freqs <= freqmax))  # Filter frequencies within the specified range
    dvv, err = np.zeros(len(inx[0])), np.zeros(len(inx[0])) # Initialize dvv and err arrays

    # Weighting function based on WXamp
    weight_func = np.log(np.abs(WXamp)) / np.log(np.abs(WXamp)).max()
    zero_idx = np.where((Wcoh < mincoh) | (delta_t > maxdt))
    wf = (weight_func + abs(np.nanmin(weight_func))) / weight_func.max()
    wf[zero_idx] = 0

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
            w[~np.isfinite(w)] = 1.0

            # Percentage of non-zero weights
            nzc_perc = np.count_nonzero(w[tindex]) / len(tindex)

            if nzc_perc >= min_nonzero:
                m, em = linear_regression(tvec[tindex], delta_t[ifreq][tindex], w[tindex], intercept_origin=True)
                dvv[ii], err[ii] = -m, em
            else:
                dvv[ii], err[ii] = np.nan, np.nan
        else:
            print('Not enough points to estimate dv/v for WCT')

    return dvv * 100, err * 100, wf

def xwt(trace_ref, trace_current, fs, ns=3, nt=0.25, vpo=12, freqmin=0.1, freqmax=8.0, nptsfreq=100):
    """
    Wavelet coherence transform (WCT).

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
       
    Returns
    """
    # Choosing a Morlet wavelet with a central frequency w0 = 6
    mother = wavelet.Morlet(6.)
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
    cwt_reference, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(x_reference, dt, dj, s0, J, mother, freqs=freqlim)
    cwt_current, _, _, _, _, _ = wavelet.cwt(x_current, dt, dj, s0, J, mother, freqs=freqlim)
    ###############################################################################################################

    scales = np.array([[kk] for kk in scales])
    invscales = np.kron(np.ones((1, nx)), 1 / scales)
    
    cfs2 = smoothCFS(invscales * abs(cwt_current) ** 2, scales, dt, ns, nt)
    cfs1 = smoothCFS(invscales * abs(cwt_reference) ** 2, scales, dt, ns, nt)
    
    crossCFS = cwt_reference * np.conj(cwt_current)
    WXamp = abs(crossCFS)
    # cross-wavelet transform operation with smoothing
    crossCFS = smoothCFS(invscales * crossCFS, scales, dt, ns, nt)
    WXspec = crossCFS / (np.sqrt(cfs1) * np.sqrt(cfs2))
    WXangle = np.angle(WXspec)
    Wcoh = abs(crossCFS) ** 2 / (cfs1 * cfs2)
    pp = 2 * np.pi * freqs
    pp2 = np.array([[kk] for kk in pp])
    WXdt = WXangle / np.kron(np.ones((1, nx)), pp2)


    return WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi

def xr_save_wct(station1, station2, components, filterid, mov_stack, taxis, dvv_df, err_df, coh_df):
    """
    Save the Wavelet Coherence Transform (WCT) results as a NetCDF file.
    
    Parameters
    ----------
    station1 : str
        The first station in the pair.
    station2 : str
        The second station in the pair.
    components : str
        The components (e.g., Z, N, E) being analyzed.
    filterid : int
        Filter ID used in the analysis.
    mov_stack : tuple
        Tuple of (start, end) representing the moving stack window.
    taxis : array-like
        Time axis corresponding to the WCT data.
    dvv_df : pandas.DataFrame
        DataFrame containing dvv data (2D).
    err_df : pandas.DataFrame
        DataFrame containing err data (2D).
    coh_df : pandas.DataFrame
        DataFrame containing coh data (2D).
    
    Returns
    -------
    None
    """
    # Construct the file path
    fn = os.path.join("WCT", f"{filterid:02d}", f"{mov_stack[0]}_{mov_stack[1]}",
                      components, f"{station1}_{station2}.nc")
    
    # Ensure the directory exists
    os.makedirs(os.path.dirname(fn), exist_ok=True)

    # Convert DataFrames to xarray.DataArrays
    dvv_da = xr.DataArray(dvv_df.values, coords=[dvv_df.index, dvv_df.columns], dims=['times', 'frequency'])
    err_da = xr.DataArray(err_df.values, coords=[err_df.index, err_df.columns], dims=['times', 'frequency'])
    coh_da = xr.DataArray(coh_df.values, coords=[coh_df.index, coh_df.columns], dims=['times', 'frequency'])

    # Combine into a single xarray.Dataset
    ds = xr.Dataset({
        'dvv': dvv_da,
        'err': err_da,
        'coh': coh_da
    })

    # Save the dataset to a NetCDF file
    ds.to_netcdf(fn)
    
    print(f"Saved WCT data to {fn}")

    # Clean up
    del dvv_da, err_da, coh_da, ds




def main(loglevel="INFO"):

    # Reconfigure logger to show the pid number in log records
    global logger
    logger = get_logger('msnoise.compute_wct_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute WCT ***')

    db = connect()
    params = get_params(db)
    taxis = get_t_axis(db)

    ns = params.wct_ns
    nt = params.wct_nt 
    vpo = params.wct_vpo 
    nptsfreq = params.wct_nptsfreq
    coda_cycles = params.dtt_codacycles 
    min_nonzero = params.dvv_min_nonzero
    wct_norm = params.wct_norm

    mov_stacks = params.mov_stack
    goal_sampling_rate = params.cc_sampling_rate
    lag_min = params.dtt_minlag
    maxdt = params.dtt_maxdt
    mincoh = params.dtt_mincoh

    logger.debug('Ready to compute')
    # Then we compute the jobs
    filters = get_filters(db, all=False)
    time.sleep(np.random.random() * 5)

    while is_dtt_next_job(db, flag='T', jobtype='WCT'):
        # TODO would it be possible to make the next 8 lines in the API ?
        jobs = get_dtt_next_job(db, flag='T', jobtype='WCT')

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        logger.info(
            "There are WCT jobs for some days to recompute for %s" % pair)
        for f in filters:
            filterid = int(f.ref)
            freqmin = f.low
            freqmax = f.high

            ref_name = pair.replace(':', '_')
            station1, station2 = pair.split(":")
            
            if station1 == station2:
                components_to_compute = params.components_to_compute_single_station
            else:
                components_to_compute = params.components_to_compute


            for components in components_to_compute:
                try:
                    ref = xr_get_ref(station1, station2, components, filterid, taxis)
                    ref = ref.CCF.values

                    if wct_norm:
                        ori_waveform = (ref/ref.max()) 
                    else:
                        ori_waveform = ref
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue
                if not len(ref):
                    continue
                

                for mov_stack in mov_stacks:
                    dvv_list = []
                    err_list = []
                    coh_list = []
                    data_dates=[]
                    try:
                        data = xr_get_ccf(station1, station2, components, filterid, mov_stack, taxis)
                    except FileNotFoundError as fullpath:
                        logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                        continue
                    logger.debug("Processing %s:%s f%i m%s %s" % (station1, station2, filterid, mov_stack, components))

                    to_search = pd.to_datetime(days)
                    to_search = to_search.append(pd.DatetimeIndex([to_search[-1]+pd.Timedelta("1d"),]))
                    # data = data[(data.index.floor('d').isin(to_search) or data.index.ceil('d').isin(to_search))]
                    data = data[data.index.floor('d').isin(to_search)]
                    data = data.dropna()
                    
                    cur = data#.CCF.values
                    if wct_norm:
                        new_waveform = (cur/cur.max()) 
                    else:
                        new_waveform = cur

                    for date, row in new_waveform.iterrows():
                        WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(ori_waveform, row.values, goal_sampling_rate, int(ns), int(nt), int(vpo), freqmin, freqmax, int(nptsfreq))
                        dvv, err, wf = get_dvv(freqs, taxis, WXamp, Wcoh, WXdt, lag_min=int(lag_min), coda_cycles=coda_cycles, mincoh=mincoh, maxdt=maxdt, min_nonzero=min_nonzero, freqmin=freqmin, freqmax=freqmax)
                        coh = get_avgcoh(freqs, taxis, Wcoh, freqmin, freqmax, lag_min=int(lag_min), coda_cycles=coda_cycles)

                        dvv_list.append(dvv)
                        err_list.append(err)
                        coh_list.append(coh)
                        data_dates.append(date)


                    if len(dvv_list) > 1:
                        inx = np.where((freqs >= freqmin) & (freqs <= freqmax))
                        
                        dvv_df = pd.DataFrame(dvv_list, columns=freqs[inx], index=data_dates)
                        err_df = pd.DataFrame(err_list, columns=freqs[inx], index=data_dates)
                        coh_df = pd.DataFrame(coh_list, columns=freqs[inx], index=data_dates)

                        # Save using the xr_save_wct function
                        xr_save_wct(station1, station2, components, filterid, mov_stack, taxis, dvv_df, err_df, coh_df)

                        del dvv_df, err_df, coh_df
                    del cur

        massive_update_job(db, jobs, "D")

    logger.info('*** Finished: Compute WCT ***')



