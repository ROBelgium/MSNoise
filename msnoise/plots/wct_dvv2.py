#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
plot_wct_dvv.py

Script to plot DVV results from merged WCT files created by s09merge_wct.py.
Creates heatmaps showing DVV variations over time and frequency.

Usage:
    python plot_wct_dvv.py --filter-id 1 --wct-id 1 --dtt-id 1 --station-pair UV01_UV02 --component ZZ
    python plot_wct_dvv.py --merged-dir DVV/WCT/WCT_MERGED --output-dir plots
"""

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import Normalize
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings('ignore')

# Import functions from the API if available
try:
    from api import get_logger, connect, get_params
    HAS_API = True
except ImportError:
    HAS_API = False
    print("Warning: MSNoise API not available. Using standalone mode.")

def setup_logging(loglevel="INFO"):
    """Set up basic logging if MSNoise API is not available"""
    if HAS_API:
        return get_logger('wct_dvv_plotter', loglevel)
    else:
        import logging
        logging.basicConfig(
            level=getattr(logging, loglevel.upper()),
            format='%(asctime)s [%(levelname)s] %(message)s'
        )
        return logging.getLogger('wct_dvv_plotter')

def find_merged_files(merged_dir="DVV/WCT/WCT_MERGED", filter_id=None, 
                     wct_id=None, dtt_id=None, station_pair=None, component=None):
    """
    Find merged WCT files based on search criteria
    
    Parameters:
    -----------
    merged_dir : str
        Directory containing merged WCT results
    filter_id : int, optional
        Filter ID to search for
    wct_id : int, optional  
        WCT parameter ID to search for
    dtt_id : int, optional
        DTT parameter ID to search for
    station_pair : str, optional
        Station pair (e.g., 'UV01_UV02')
    component : str, optional
        Component (e.g., 'ZZ')
        
    Returns:
    --------
    list
        List of file paths matching the criteria
    """
    
    # Build search pattern
    pattern_parts = []
    
    if filter_id is not None:
        pattern_parts.append(f"f{filter_id:02d}")
    else:
        pattern_parts.append("f*")
        
    if wct_id is not None:
        pattern_parts.append(f"wct{wct_id:02d}")
    else:
        pattern_parts.append("wct*")
        
    if dtt_id is not None:
        pattern_parts.append(f"dtt{dtt_id:02d}")
    else:
        pattern_parts.append("dtt*")
        
    pattern_parts.append("*")  # moving stack
    
    if component is not None:
        pattern_parts.append(component)
    else:
        pattern_parts.append("*")
        
    if station_pair is not None:
        pattern_parts.append(f"{station_pair}.nc")
    else:
        pattern_parts.append("*.nc")
    
    search_pattern = os.path.join(merged_dir, *pattern_parts)
    files = glob.glob(search_pattern)
    
    return sorted(files)

def inspect_data_structure(data, filepath):
    """
    Inspect and print the structure of loaded data for debugging
    
    Parameters:
    -----------
    data : xr.Dataset
        Dataset to inspect
    filepath : str
        Path to the file for context
    """
    print(f"\n=== Data Structure for {os.path.basename(filepath)} ===")
    print(f"Variables: {list(data.data_vars.keys())}")
    print(f"Coordinates: {list(data.coords.keys())}")
    
    for var in data.data_vars:
        da = data[var]
        print(f"\n{var}:")
        print(f"  Dimensions: {da.dims}")
        print(f"  Shape: {da.shape}")
        print(f"  Coords: {list(da.coords.keys())}")
        
        # Print coordinate info
        for coord in da.coords:
            coord_data = da.coords[coord]
            if coord_data.size <= 10:
                print(f"  {coord}: {coord_data.values}")
            else:
                print(f"  {coord}: {coord_data.size} values, range [{coord_data.values[0]} ... {coord_data.values[-1]}]")
    print("=" * 50)

def load_wct_data(filepath):
    """
    Load WCT data from NetCDF file
    
    Parameters:
    -----------
    filepath : str
        Path to the NetCDF file
        
    Returns:
    --------
    xr.Dataset
        Dataset containing DVV, error, and coherence data
    """
    try:
        ds = xr.load_dataset(filepath)
        
        # Debug: inspect the data structure
        inspect_data_structure(ds, filepath)
        
        return ds
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None

def parse_filename(filepath):
    """
    Parse information from the file path
    
    Parameters:
    -----------
    filepath : str
        Path to the NetCDF file
        
    Returns:
    --------
    dict
        Dictionary with parsed information
    """
    parts = filepath.split(os.sep)
    filename = os.path.basename(filepath)
    
    info = {
        'filepath': filepath,
        'filename': filename,
        'station_pair': filename.replace('.nc', ''),
    }
    
    # Extract information from path
    for part in parts:
        if part.startswith('f') and part[1:].isdigit():
            info['filter_id'] = int(part[1:])
        elif part.startswith('wct') and part[3:].isdigit():
            info['wct_id'] = int(part[3:])
        elif part.startswith('dtt') and part[3:].isdigit():
            info['dtt_id'] = int(part[3:])
        elif '_' in part and any(char.isdigit() for char in part):
            # This might be the moving stack
            info['moving_stack'] = part
        elif part in ['ZZ', 'ZN', 'ZE', 'NZ', 'NN', 'NE', 'EZ', 'EN', 'EE']:
            info['component'] = part
    
    return info

def create_dvv_heatmap(data, info, variable='dvv', output_dir='plots', 
                      figsize=(16, 10), cmap='RdBu_r', vmin=None, vmax=None):
    """
    Create DVV heatmap plot similar to the reference implementation
    
    Parameters:
    -----------
    data : xr.Dataset
        Dataset containing the WCT results
    info : dict
        Information about the data (from parse_filename)
    variable : str
        Variable to plot ('dvv', 'err', or 'coh')
    output_dir : str
        Directory to save plots
    figsize : tuple
        Figure size (width, height)
    cmap : str
        Colormap name
    vmin, vmax : float, optional
        Color scale limits
    """
    
    if variable not in data.data_vars:
        print(f"Variable '{variable}' not found in dataset. Available: {list(data.data_vars.keys())}")
        return
    
    # Get the data array
    da = data[variable]
    
    print(f"Data dimensions: {da.dims}")
    print(f"Data shape: {da.shape}")
    
    # Convert to pandas DataFrame for easier manipulation
    if da.shape[0] == 1:
        # Single time point - need to handle specially
        print("Single time point detected - creating heatmap with single time slice")
        
        # Extract the data
        freq_vals = da.coords['frequency'].values
        time_vals = da.coords['times'].values
        data_vals = da.values
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set colormap and limits based on variable type
        if variable == 'dvv':
            cmap = 'RdBu_r'
            if vmin is None or vmax is None:
                data_clean = data_vals[~np.isnan(data_vals)]
                if len(data_clean) > 0:
                    if vmin is None:
                        vmin = np.percentile(data_clean, 2)
                    if vmax is None:
                        vmax = np.percentile(data_clean, 98)
                else:
                    vmin, vmax = -1, 1
            cbar_label = 'dv/v (%)'
        elif variable == 'coh':
            cmap = 'RdYlGn'
            vmin, vmax = 0, 1
            cbar_label = 'Coherence'
        else:  # err
            cmap = 'viridis'
            if vmin is None or vmax is None:
                data_clean = data_vals[~np.isnan(data_vals)]
                if len(data_clean) > 0:
                    if vmin is None:
                        vmin = np.percentile(data_clean, 2)
                    if vmax is None:
                        vmax = np.percentile(data_clean, 98)
                else:
                    vmin, vmax = 0, 1
            cbar_label = 'Error'
        
        # For single time point, create a thin time window around the point
        time_center = pd.to_datetime(time_vals[0])
        time_start = time_center - pd.Timedelta(hours=12)
        time_end = time_center + pd.Timedelta(hours=12)
        time_array = pd.date_range(time_start, time_end, periods=3)
        
        # Replicate the data across the small time window
        data_2d = np.tile(data_vals, (len(time_array), 1))
        
        # Create meshgrid for plotting
        Time, Freq = np.meshgrid(time_array, freq_vals)
        
        # Create the heatmap
        mesh = ax.pcolormesh(Time, Freq, data_2d.T, cmap=cmap, vmin=vmin, vmax=vmax)
        
    else:
        # Multiple time points - normal heatmap
        fig, ax = plt.subplots(figsize=figsize)
        
        # Convert to DataFrame
        df = da.to_dataframe().unstack(level='frequency')
        df = df.droplevel(0, axis=1)  # Remove variable name level
        
        # Set colormap and limits
        if variable == 'dvv':
            cmap = 'RdBu_r'
            if vmin is None or vmax is None:
                data_clean = df.values[~np.isnan(df.values)]
                if len(data_clean) > 0:
                    if vmin is None:
                        vmin = np.percentile(data_clean, 2)
                    if vmax is None:
                        vmax = np.percentile(data_clean, 98)
                else:
                    vmin, vmax = -1, 1
            cbar_label = 'dv/v (%)'
        elif variable == 'coh':
            cmap = 'RdYlGn'
            vmin, vmax = 0, 1
            cbar_label = 'Coherence'
        else:  # err
            cmap = 'viridis'
            if vmin is None or vmax is None:
                data_clean = df.values[~np.isnan(df.values)]
                if len(data_clean) > 0:
                    if vmin is None:
                        vmin = np.percentile(data_clean, 2)
                    if vmax is None:
                        vmax = np.percentile(data_clean, 98)
                else:
                    vmin, vmax = 0, 1
            cbar_label = 'Error'
        
        # Create full date range (fill gaps with NaN)
        full_index = pd.date_range(start=df.index.min(), end=df.index.max(), freq='D')
        df_reindexed = df.reindex(full_index)
        
        time = pd.to_datetime(full_index)
        freq = df.columns.astype(float)
        data_array = np.ma.masked_invalid(df_reindexed.values)
        
        # Create the heatmap
        mesh = ax.pcolormesh(time, freq, data_array.T, cmap=cmap, vmin=vmin, vmax=vmax)
    
    # Add colorbar
    cbar = plt.colorbar(mesh, ax=ax)
    cbar.set_label(cbar_label, rotation=270, labelpad=15, fontsize=18)
    
    # Customize the plot (similar to your reference)
    station_pair = info.get('station_pair', 'Unknown')
    ax.set_title(f"{station_pair} - {variable.upper()}", fontsize=28)
    ax.set_ylabel('Frequency (Hz)', fontsize=24)
    ax.set_xlabel('Date', fontsize=24)
    
    # Tick parameters
    ax.tick_params(axis='x', which='major', length=10, width=2, labelsize=20)
    ax.tick_params(axis='x', which='minor', length=6, width=1)
    ax.tick_params(axis='y', which='major', length=8, width=1.5, labelsize=20)
    
    # Format dates on x-axis
    if da.shape[0] > 1:  # Only for multi-time data
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        plt.setp(ax.get_xticklabels(), rotation=45)
    else:
        # For single time point, format differently
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    
    plt.tight_layout()
    
    # Save plot
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_filename = f"{variable}_{info['filename'].replace('.nc', '')}"
    if 'filter_id' in info:
        output_filename += f"_f{info['filter_id']:02d}"
    if 'wct_id' in info and 'dtt_id' in info:
        output_filename += f"_wct{info['wct_id']:02d}_dtt{info['dtt_id']:02d}"
    output_filename += '.png'
    
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved heatmap plot: {output_path}")
    
    plt.close()

def create_summary_plot(data, info, output_dir='plots', figsize=(16, 12)):
    """
    Create a summary plot with DVV, error, and coherence heatmap subplots
    
    Parameters:
    -----------
    data : xr.Dataset
        Dataset containing the WCT results
    info : dict
        Information about the data
    output_dir : str
        Directory to save plots
    figsize : tuple
        Figure size (width, height)
    """
    
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
    
    variables = ['dvv', 'err', 'coh']
    cmaps = ['RdBu_r', 'viridis', 'RdYlGn']
    labels = ['dv/v (%)', 'Error', 'Coherence']
    
    for i, (var, cmap, label) in enumerate(zip(variables, cmaps, labels)):
        if var not in data.data_vars:
            continue
            
        da = data[var]
        
        # Handle single time point by creating small time window
        if da.shape[0] == 1:
            # Single time point
            freq_vals = da.coords['frequency'].values
            time_vals = da.coords['times'].values
            data_vals = da.values
            
            # Create small time window
            time_center = pd.to_datetime(time_vals[0])
            time_start = time_center - pd.Timedelta(hours=12)
            time_end = time_center + pd.Timedelta(hours=12)
            time_array = pd.date_range(time_start, time_end, periods=3)
            
            # Replicate data
            data_2d = np.tile(data_vals, (len(time_array), 1))
            
            # Set color limits
            if var == 'dvv':
                data_clean = data_vals[~np.isnan(data_vals)]
                if len(data_clean) > 0:
                    vmin, vmax = np.percentile(data_clean, [2, 98])
                else:
                    vmin, vmax = -1, 1
            elif var == 'coh':
                vmin, vmax = 0, 1
            else:  # err
                data_clean = data_vals[~np.isnan(data_vals)]
                if len(data_clean) > 0:
                    vmin, vmax = np.percentile(data_clean, [2, 98])
                else:
                    vmin, vmax = 0, 1
            
            # Create meshgrid and plot
            Time, Freq = np.meshgrid(time_array, freq_vals)
            mesh = axes[i].pcolormesh(Time, Freq, data_2d.T, cmap=cmap, vmin=vmin, vmax=vmax)
            
        else:
            # Multiple time points - normal heatmap
            # Convert to DataFrame
            df = da.to_dataframe().unstack(level='frequency')
            df = df.droplevel(0, axis=1)
            
            # Create full date range
            full_index = pd.date_range(start=df.index.min(), end=df.index.max(), freq='D')
            df_reindexed = df.reindex(full_index)
            
            time = pd.to_datetime(full_index)
            freq = df.columns.astype(float)
            data_array = np.ma.masked_invalid(df_reindexed.values)
            
            # Set color limits
            if var == 'dvv':
                data_clean = df.values[~np.isnan(df.values)]
                if len(data_clean) > 0:
                    vmin, vmax = np.percentile(data_clean, [2, 98])
                else:
                    vmin, vmax = -1, 1
            elif var == 'coh':
                vmin, vmax = 0, 1
            else:  # err
                data_clean = df.values[~np.isnan(df.values)]
                if len(data_clean) > 0:
                    vmin, vmax = np.percentile(data_clean, [2, 98])
                else:
                    vmin, vmax = 0, 1
            
            mesh = axes[i].pcolormesh(time, freq, data_array.T, cmap=cmap, vmin=vmin, vmax=vmax)
        
        # Customize subplot
        axes[i].set_ylabel('Frequency (Hz)', fontsize=18)
        axes[i].set_title(f'{var.upper()}', fontsize=20)
        axes[i].tick_params(axis='both', which='major', labelsize=16)
        
        # Add colorbar
        cbar = fig.colorbar(mesh, ax=axes[i], shrink=0.8)
        cbar.set_label(label, rotation=270, labelpad=15, fontsize=16)
    
    # Format x-axis only on bottom subplot
    axes[-1].set_xlabel('Date', fontsize=18)
    if data[variables[0]].shape[0] > 1:
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        axes[-1].xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        plt.setp(axes[-1].get_xticklabels(), rotation=45)
    else:
        axes[-1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    
    # Main title
    station_pair = info.get('station_pair', 'Unknown')
    component = info.get('component', 'Unknown')
    main_title = f"WCT Results: {station_pair} - {component}"
    fig.suptitle(main_title, fontsize=24, y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    
    # Save plot
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_filename = f"summary_{info['filename'].replace('.nc', '')}"
    if 'filter_id' in info:
        output_filename += f"_f{info['filter_id']:02d}"
    if 'wct_id' in info and 'dtt_id' in info:
        output_filename += f"_wct{info['wct_id']:02d}_dtt{info['dtt_id']:02d}"
    output_filename += '.png'
    
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved summary heatmap plot: {output_path}")
    
    plt.close()

def plot_frequency_timeseries(data, info, frequencies=None, output_dir='plots', figsize=(12, 6)):
    """
    Plot DVV time series for specific frequencies
    
    Parameters:
    -----------
    data : xr.Dataset
        Dataset containing the WCT results
    info : dict
        Information about the data
    frequencies : list, optional
        List of frequencies to plot. If None, uses a few representative frequencies
    output_dir : str
        Directory to save plots
    figsize : tuple
        Figure size (width, height)
    """
    
    if 'dvv' not in data.data_vars:
        print("DVV data not found in dataset")
        return
    
    da = data['dvv']
    
    # Select frequencies to plot
    if frequencies is None:
        all_freqs = da.coords['frequency'].values
        if len(all_freqs) >= 5:
            # Select 5 representative frequencies
            indices = np.linspace(0, len(all_freqs)-1, 5, dtype=int)
            frequencies = all_freqs[indices]
        else:
            frequencies = all_freqs
    
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = plt.cm.viridis(np.linspace(0, 1, len(frequencies)))
    
    for freq, color in zip(frequencies, colors):
        # Find closest frequency
        freq_idx = np.argmin(np.abs(da.coords['frequency'].values - freq))
        actual_freq = da.coords['frequency'].values[freq_idx]
        
        # Extract time series for this frequency
        ts = da.isel(frequency=freq_idx)
        
        # Plot with error bars if available
        if 'err' in data.data_vars:
            err_ts = data['err'].isel(frequency=freq_idx)
            ax.errorbar(ts.coords['times'], ts.values, yerr=err_ts.values, 
                       color=color, alpha=0.7, label=f'{actual_freq:.2f} Hz')
        else:
            ax.plot(ts.coords['times'], ts.values, color=color, 
                   label=f'{actual_freq:.2f} Hz')
    
    ax.set_xlabel('Date')
    ax.set_ylabel('dv/v (%)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Format dates
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.setp(ax.get_xticklabels(), rotation=45)
    
    # Title
    title_parts = []
    if 'station_pair' in info:
        title_parts.append(f"Stations: {info['station_pair']}")
    if 'component' in info:
        title_parts.append(f"Component: {info['component']}")
    
    title = f"DVV Time Series: {' | '.join(title_parts)}"
    ax.set_title(title)
    
    plt.tight_layout()
    
    # Save plot
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_filename = f"timeseries_{info['filename'].replace('.nc', '')}"
    if 'filter_id' in info:
        output_filename += f"_f{info['filter_id']:02d}"
    if 'wct_id' in info and 'dtt_id' in info:
        output_filename += f"_wct{info['wct_id']:02d}_dtt{info['dtt_id']:02d}"
    output_filename += '.png'
    
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved time series plot: {output_path}")
    
    plt.close()

def parse_frequency_ranges(ranges_str):
    """
    Parse frequency ranges string into list of tuples
    
    Parameters:
    -----------
    ranges_str : str
        String like "[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]"
        
    Returns:
    --------
    list of tuples
        List of (min_freq, max_freq) tuples
    """
    import ast
    try:
        # Remove outer brackets if present and split by '], ['
        ranges_str = ranges_str.strip()
        if ranges_str.startswith('[') and ranges_str.endswith(']'):
            # Handle nested list format
            ranges = ast.literal_eval(ranges_str)
            if isinstance(ranges[0], list):
                return [tuple(r) for r in ranges]
            else:
                # Single range
                return [tuple(ranges)]
        else:
            # Parse manually
            ranges = []
            parts = ranges_str.split('], [')
            for part in parts:
                part = part.strip('[]')
                values = [float(x.strip()) for x in part.split(',')]
                if len(values) == 2:
                    ranges.append(tuple(values))
            return ranges
    except Exception as e:
        print(f"Warning: Could not parse frequency ranges '{ranges_str}': {e}")
        return [(0.5, 1.0), (1.0, 2.0), (2.0, 4.0)]  # Default ranges

def filter_data_by_date_range(data, start_date, end_date):
    """
    Filter data by date range
    
    Parameters:
    -----------
    data : xr.Dataset
        Dataset to filter
    start_date : str
        Start date in YYYY-MM-DD format
    end_date : str  
        End date in YYYY-MM-DD format
        
    Returns:
    --------
    xr.Dataset
        Filtered dataset
    """
    try:
        if start_date != "1970-01-01":
            start_pd = pd.to_datetime(start_date)
            data = data.sel(times=data.times >= start_pd)
        
        if end_date != "2100-01-01":
            end_pd = pd.to_datetime(end_date)
            data = data.sel(times=data.times <= end_pd)
            
        return data
    except Exception as e:
        print(f"Warning: Could not filter by date range: {e}")
        return data

def create_spectrum_plot(data, info, output_dir='plots', figsize=(12, 6)):
    """
    Create a spectrum plot for single time point data
    
    Parameters:
    -----------
    data : xr.Dataset
        Dataset containing single time point WCT results
    info : dict
        Information about the data
    output_dir : str
        Directory to save plots
    figsize : tuple
        Figure size (width, height)
    """
    
    fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
    
    variables = ['dvv', 'err', 'coh'] 
    colors = ['red', 'blue', 'green']
    labels = ['dv/v (%)', 'Error', 'Coherence']
    
    for i, (var, color, label) in enumerate(zip(variables, colors, labels)):
        if var not in data.data_vars:
            continue
            
        da = data[var]
        
        if da.shape[0] == 1:  # Single time point
            freqs = da.coords['frequency'].values
            values = da.values.squeeze()
            
            axes[i].plot(freqs, values, color=color, linewidth=2, label=label)
            axes[i].set_ylabel(label)
            axes[i].grid(True, alpha=0.3)
            axes[i].set_title(f'{var.upper()} Spectrum', fontsize=11)
            
            # Set appropriate y-limits
            if var == 'dvv':
                # Center around zero for DVV
                abs_max = np.nanmax(np.abs(values))
                if abs_max > 0:
                    axes[i].set_ylim(-abs_max*1.1, abs_max*1.1)
            elif var == 'coh':
                # Coherence should be between 0 and 1
                axes[i].set_ylim(0, 1.1)
    
    # Format x-axis only on bottom subplot
    axes[-1].set_xlabel('Frequency (Hz)')
    axes[-1].set_xlim(freqs.min(), freqs.max())
    
    # Main title
    title_parts = []
    if 'station_pair' in info:
        title_parts.append(f"Stations: {info['station_pair']}")
    if 'component' in info:
        title_parts.append(f"Component: {info['component']}")
    
    # Get the single time point for the title
    time_str = pd.to_datetime(data.coords['times'].values[0]).strftime('%Y-%m-%d')
    title_parts.append(f"Date: {time_str}")
    
    main_title = f"WCT Spectrum: {' | '.join(title_parts)}"
    fig.suptitle(main_title, fontsize=14, y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    
    # Save plot
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_filename = f"spectrum_{info['filename'].replace('.nc', '')}"
    if 'filter_id' in info:
        output_filename += f"_f{info['filter_id']:02d}"
    if 'wct_id' in info and 'dtt_id' in info:
        output_filename += f"_wct{info['wct_id']:02d}_dtt{info['dtt_id']:02d}"
    output_filename += '.png'
    
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved spectrum plot: {output_path}")
    
    plt.close()

def create_frequency_range_plot(data, info, freq_ranges, output_dir='plots', 
                               start_date=None, end_date=None, figsize=(12, 8)):
    """
    Create DVV plot averaged over frequency ranges
    
    Parameters:
    -----------
    data : xr.Dataset
        Dataset containing the WCT results
    info : dict
        Information about the data
    freq_ranges : list of tuples
        List of (min_freq, max_freq) frequency ranges
    output_dir : str
        Directory to save plots
    start_date, end_date : str, optional
        Date range for filtering
    figsize : tuple
        Figure size (width, height)
    """
    
    if 'dvv' not in data.data_vars:
        print("DVV data not found in dataset")
        return
    
    # Filter by date range
    if start_date or end_date:
        data = filter_data_by_date_range(data, start_date or "1970-01-01", 
                                       end_date or "2100-01-01")
    
    da = data['dvv']
    
    fig, ax = plt.subplots(figsize=figsize)
    
    colors = plt.cm.Set1(np.linspace(0, 1, len(freq_ranges)))
    
    for (fmin, fmax), color in zip(freq_ranges, colors):
        # Select frequency range
        freq_mask = (da.coords['frequency'] >= fmin) & (da.coords['frequency'] <= fmax)
        da_range = da.sel(frequency=freq_mask)
        
        if da_range.sizes['frequency'] == 0:
            print(f"Warning: No frequencies found in range {fmin}-{fmax} Hz")
            continue
        
        # Average over frequency range
        dvv_avg = da_range.mean(dim='frequency')
        
        # Plot with error bars if available
        if 'err' in data.data_vars:
            err_range = data['err'].sel(frequency=freq_mask)
            err_avg = err_range.mean(dim='frequency')
            ax.errorbar(dvv_avg.coords['times'], dvv_avg.values, yerr=err_avg.values,
                       color=color, alpha=0.7, label=f'{fmin}-{fmax} Hz')
        else:
            ax.plot(dvv_avg.coords['times'], dvv_avg.values, color=color,
                   label=f'{fmin}-{fmax} Hz', linewidth=2)
    
    ax.set_xlabel('Date')
    ax.set_ylabel('dv/v (%)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Format dates
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    plt.setp(ax.get_xticklabels(), rotation=45)
    
    # Title
    title_parts = []
    if 'station_pair' in info:
        title_parts.append(f"Stations: {info['station_pair']}")
    if 'component' in info:
        title_parts.append(f"Component: {info['component']}")
    
    title = f"DVV Frequency Ranges: {' | '.join(title_parts)}"
    ax.set_title(title)
    
    plt.tight_layout()
    
    # Save plot
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    output_filename = f"freq_ranges_{info['filename'].replace('.nc', '')}"
    if 'filter_id' in info:
        output_filename += f"_f{info['filter_id']:02d}"
    if 'wct_id' in info and 'dtt_id' in info:
        output_filename += f"_wct{info['wct_id']:02d}_dtt{info['dtt_id']:02d}"
    output_filename += '.png'
    
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved frequency ranges plot: {output_path}")
    
    if output_filename:  # show parameter equivalent
        plt.show()
    else:
        plt.close()

def main(mov_stackid=0, components="ZZ", filterid=1, wctid=1, dttid=1, 
         pairs=None, showALL=False, start="1970-01-01", end="2100-01-01", 
         visualize="dvv", ranges="[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]", 
         show=True, outfile=None, loglevel="INFO"):
    """
    Main function adapted for MSNoise click command structure
    
    Parameters:
    -----------
    mov_stackid : int
        Moving stack ID (0 for all)
    components : str
        Components to plot (e.g., "ZZ")
    filterid : int
        Filter ID
    wctid : int
        WCT parameter ID
    dttid : int
        DTT parameter ID
    pairs : list, optional
        List of station pairs to plot
    showALL : bool
        Show ALL line (not used in this implementation)
    start : str
        Start date (YYYY-MM-DD)
    end : str
        End date (YYYY-MM-DD)
    visualize : str
        Visualization type ('dvv', 'heatmap', 'timeseries', 'ranges')
    ranges : str
        Frequency ranges string
    show : bool
        Show plots interactively
    outfile : str, optional
        Output filename
    loglevel : str
        Logging level
    """
    
    # Set up logging
    logger = setup_logging(loglevel)
    logger.info("Starting WCT DVV plotting")
    
    # Determine merged directory - try both merged and direct paths
    merged_dirs = [
        "DVV/WCT/WCT_MERGED",
        "DVV/WCT/DTT"
    ]
    
    merged_dir = None
    for d in merged_dirs:
        if os.path.exists(d):
            merged_dir = d
            break
    
    if merged_dir is None:
        merged_dir = "DVV/WCT/WCT_MERGED"  # Default
        logger.warning(f"Using default merged directory: {merged_dir}")
    
    logger.info(f"Searching for files in {merged_dir}")
    
    # Convert mov_stackid to moving stack pattern
    mov_stack_pattern = None
    if mov_stackid != 0:
        mov_stack_pattern = f"*{mov_stackid}*"
    
    # Find files for each pair or all pairs
    all_files = []
    
    if pairs:
        # Plot specific pairs
        for pair in pairs:
            files = find_merged_files(
                merged_dir=merged_dir,
                filter_id=filterid,
                wct_id=wctid,
                dtt_id=dttid,
                station_pair=pair,
                component=components
            )
            all_files.extend(files)
    else:
        # Find all files matching criteria
        all_files = find_merged_files(
            merged_dir=merged_dir,
            filter_id=filterid,
            wct_id=wctid,
            dtt_id=dttid,
            station_pair=None,
            component=components
        )
    
    if not all_files:
        logger.error("No files found matching the criteria")
        return 1
    
    logger.info(f"Found {len(all_files)} files to process")
    
    # Parse frequency ranges for range plots
    freq_ranges = parse_frequency_ranges(ranges)
    
    # Set output directory
    output_dir = "plots"
    if outfile:
        if os.path.dirname(outfile):
            output_dir = os.path.dirname(outfile)
    
    # Process each file
    for filepath in all_files:
        logger.info(f"Processing {filepath}")
        
        # Load data
        data = load_wct_data(filepath)
        if data is None:
            continue
        
        # Parse file information
        info = parse_filename(filepath)
        logger.info(f"File info: {info}")
        
        # Create plots based on visualization type
        try:
            if visualize == "dvv" or visualize == "heatmap":
                create_dvv_heatmap(data, info, variable='dvv', 
                                 output_dir=output_dir)
                create_summary_plot(data, info, output_dir=output_dir)
                
            elif visualize == "timeseries":
                plot_frequency_timeseries(data, info, output_dir=output_dir)
                
            elif visualize == "ranges":
                create_frequency_range_plot(data, info, freq_ranges, 
                                          output_dir=output_dir, 
                                          start_date=start, end_date=end)
            else:
                # Default: create all plot types
                create_dvv_heatmap(data, info, variable='dvv', 
                                 output_dir=output_dir)
                create_summary_plot(data, info, output_dir=output_dir)
                plot_frequency_timeseries(data, info, output_dir=output_dir)
                create_frequency_range_plot(data, info, freq_ranges, 
                                          output_dir=output_dir,
                                          start_date=start, end_date=end)
                
        except Exception as e:
            logger.error(f"Error creating plots for {filepath}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    logger.info(f"Plotting completed. Check {output_dir} for results.")
    return 0

# Standalone execution for testing
def standalone_main():
    """Standalone main function for testing with argparse"""
    parser = argparse.ArgumentParser(description='Plot DVV results from merged WCT files')
    parser.add_argument('--merged-dir', default='DVV/WCT/WCT_MERGED',
                       help='Directory containing merged WCT results')
    parser.add_argument('--output-dir', default='plots',
                       help='Directory to save plots')
    parser.add_argument('--filter-id', type=int, default=1,
                       help='Filter ID to plot')
    parser.add_argument('--wct-id', type=int, default=1,
                       help='WCT parameter ID to plot')
    parser.add_argument('--dtt-id', type=int, default=1,
                       help='DTT parameter ID to plot')
    parser.add_argument('--station-pair',
                       help='Station pair to plot (e.g., UV01_UV02)')
    parser.add_argument('--component', default='ZZ',
                       help='Component to plot (e.g., ZZ)')
    parser.add_argument('--visualize', default='dvv', 
                       choices=['dvv', 'heatmap', 'timeseries', 'ranges'],
                       help='Type of visualization to create')
    parser.add_argument('--ranges', default="[0.5, 1.0], [1.0, 2.0], [2.0, 4.0]",
                       help='Frequency ranges for range plots')
    parser.add_argument('--start', default="1970-01-01",
                       help='Start date (YYYY-MM-DD)')
    parser.add_argument('--end', default="2100-01-01", 
                       help='End date (YYYY-MM-DD)')
    parser.add_argument('--loglevel', default='INFO',
                       help='Logging level')
    
    args = parser.parse_args()
    
    # Convert to main function parameters
    pairs = [args.station_pair] if args.station_pair else None
    
    return main(
        mov_stackid=0,
        components=args.component,
        filterid=args.filter_id,
        wctid=args.wct_id,
        dttid=args.dtt_id,
        pairs=pairs,
        showALL=False,
        start=args.start,
        end=args.end,
        visualize=args.visualize,
        ranges=args.ranges,
        show=True,
        outfile=None,
        loglevel=args.loglevel
    )

if __name__ == "__main__":
    sys.exit(standalone_main())
