import numpy as np
from scipy.signal import wiener, butter, filtfilt
import glob
import os
import matplotlib.pyplot as plt
from obspy import Stream, read, Trace
from concurrent.futures import ProcessPoolExecutor

def find_segments(data):

    """ Identify continuous non-Nan segments in xarray """

    current_segment = []
    continuous_segments = []

    for i in range(data.shape[0]):
        if not data[i, :].isnull().all():
            current_segment.append(i) #add index of array 
        else:
            if len(current_segment) > 0:
                continuous_segments.append(current_segment)
                current_segment = []

    if len(current_segment) > 0:
        continuous_segments.append(current_segment)

    return continuous_segments

def wiener_filt(data, M, N):

    ccfs = data['CCF']    
    segments = find_segments(ccfs) #get segments of continuous data to apply wiener filter

    filtered_ccfs = ccfs.copy(deep=True)

    for segment in segments:
        segment_data = ccfs[segment, :].values
        filtered_segment = wiener(segment_data, (M,N))
        filtered_ccfs[segment,:] = filtered_segment

    filtered_data = data.copy(deep=True)
    filtered_data['CCF'] = filtered_ccfs

    return filtered_data

