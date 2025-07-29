#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
s09merge_wct.py

This module merges daily WCT files for each station pair and parameter set
into a consolidated file with hierarchical output directories:
WCT_MERGED/[filter_id]/[mov_stack]/[component]/[station_pair].nc

This script is designed to be called by MSNoise through the command line
interface using the 'msnoise cc dvv merge_wct' command.

Usage:
msnoise cc dvv merge_wct
msnoise cc dvv merge_wct --wct-dir /path/to/wct --output-dir /path/to/output
msnoise -t 8 cc dvv merge_wct
"""

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import datetime
import traceback
from .api import (get_logger, xr_create_or_open, xr_save_and_close, 
                 xr_insert_or_update, validate_stack_data, xr_save_wct2)

def setup_progress_tracking(output_dir):
    """Set up progress tracking with log file"""
    log_file = os.path.join(output_dir, 'merge_progress.log')
    
    # Initialize log file if it doesn't exist
    if not os.path.exists(log_file):
        try:
            os.makedirs(os.path.dirname(log_file), exist_ok=True)
            with open(log_file, 'w') as f:
                f.write(f"# WCT Merge Progress Log - Created {datetime.datetime.now()}\n")
                f.write("# Format: filter|movstack|component|pair\n")
        except Exception as e:
            print(f"Warning: Could not create log file {log_file}: {e}")
    
    return log_file

def load_completed_merges(log_file):
    """Load the set of completed merges from log file"""
    completed = set()
    if os.path.exists(log_file):
        try:
            with open(log_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        completed.add(line)
        except Exception as e:
            print(f"Warning: Could not read log file {log_file}: {e}")
    return completed

def save_completed_merge(log_file, merge_key):
    """Save a completed merge to the log file"""
    try:
        with open(log_file, 'a') as f:
            f.write(f"{merge_key}\n")
    except Exception as e:
        print(f"Warning: Could not write to log file {log_file}: {e}")

def create_merge_key(filter_dir, movstack_dir, component_dir, pair_dir):
    """Create a unique key for a merge operation"""
    return f"{filter_dir}|{movstack_dir}|{component_dir}|{pair_dir}"

def parse_pair_name(pair_folder):
    """Parse station pair folder name into station1 and station2"""
    return pair_folder.split('_', 1)

def get_filter_directories(wct_dir):
    """Get list of filter directories in the WCT directory"""
    if not os.path.exists(wct_dir):
        return []
    
    # Look for filter directories (numbered folders like "01", "02", etc.)
    filters = []
    for item in os.listdir(wct_dir):
        if os.path.isdir(os.path.join(wct_dir, item)) and item.isdigit():
            filters.append(item)
    
    return sorted(filters)

def get_movstack_directories(wct_dir, filter_dir):
    """Get list of moving stack directories for a filter"""
    filter_path = os.path.join(wct_dir, filter_dir)
    movstacks = []
    for item in os.listdir(filter_path):
        if os.path.isdir(os.path.join(filter_path, item)):
            movstacks.append(item)
    
    return sorted(movstacks)

def get_component_directories(wct_dir, filter_dir, movstack_dir):
    """Get list of component directories for a filter and moving stack"""
    path = os.path.join(wct_dir, filter_dir, movstack_dir)
    components = []
    for item in os.listdir(path):
        if os.path.isdir(os.path.join(path, item)):
            components.append(item)
    
    return sorted(components)

def get_station_pair_directories(wct_dir, filter_dir, movstack_dir, component_dir):
    """Get list of station pair directories"""
    path = os.path.join(wct_dir, filter_dir, movstack_dir, component_dir)
    pairs = []
    for item in os.listdir(path):
        if os.path.isdir(os.path.join(path, item)) and '_' in item:
            pairs.append(item)
    
    return sorted(pairs)

def merge_daily_files(wct_dir, output_dir, filter_dir, movstack_dir, component_dir, pair_dir, logger):
    """
    Merge daily files for a specific station pair and parameter set
    
    Parameters:
    -----------
    wct_dir : str
        Path to the WCT directory
    output_dir : str
        Path to the output directory
    filter_dir : str
        Name of the filter directory
    movstack_dir : str
        Name of the moving stack directory
    component_dir : str
        Name of the component directory
    pair_dir : str
        Name of the station pair directory
    logger : logging.Logger
        Logger object
    
    Returns:
    --------
    bool
        True if successful, False otherwise
    """
    full_path = os.path.join(wct_dir, filter_dir, movstack_dir, component_dir, pair_dir)
    station1, station2 = parse_pair_name(pair_dir)
    
    # Find all date files (*.npz)
    files = glob.glob(os.path.join(full_path, "*.npz"))
    
    if not files:
        logger.warning(f"No .npz files found in {full_path}")
        return False
    
    logger.info(f"Found {len(files)} files for filter {filter_dir}, movstack {movstack_dir}, component {component_dir}, pair {pair_dir}")
    
    # Initialize lists to store data
    dates = []
    dvv_values = []
    err_values = []
    coh_values = []
    
    # Pre-define variables that will be used later
    filterid = int(filter_dir)
    component = component_dir
    taxis = None
    
    # Load each file and extract data
    for file in sorted(files):
        try:
            data = np.load(file, allow_pickle=True)
            
            # Try to extract date from filename first (more reliable)
            filename = os.path.basename(file)
            date_str = filename.split('.')[0]  # Remove .npz extension
            
            try:
                date = pd.to_datetime(date_str)
            except ValueError:
                # If filename doesn't contain a valid date, try to get it from the data
                if 'date' in data:
                    date_val = data['date']
                    # Safely convert various date formats
                    if isinstance(date_val, np.ndarray):
                        # If it's an array of size 1, extract the item
                        if date_val.size == 1:
                            try:
                                date = pd.to_datetime(date_val.item())
                            except:
                                date = pd.to_datetime(str(date_val.item()))
                        else:
                            # If it's a larger array, use the first element
                            date = pd.to_datetime(date_val[0])
                    else:
                        # If it's already a scalar
                        date = pd.to_datetime(date_val)
                else:
                    logger.warning(f"Could not determine date for {file}, skipping")
                    continue
            
            # Extract data values, ensuring they're in the correct format
            if 'dvv' in data:
                dvv = data['dvv']
                if isinstance(dvv, np.ndarray) and dvv.ndim > 0:
                    dvv_values.append(dvv)
                else:
                    logger.warning(f"Invalid dvv data in {file}, skipping")
                    continue
            else:
                logger.warning(f"No dvv data in {file}, skipping")
                continue
                
            if 'err' in data:
                err = data['err']
                if isinstance(err, np.ndarray) and err.ndim > 0:
                    err_values.append(err)
                else:
                    # Use zeros array of same shape as dvv if err is invalid
                    err_values.append(np.zeros_like(dvv))
            else:
                # Use zeros array of same shape as dvv if err is missing
                err_values.append(np.zeros_like(dvv))
                
            if 'coh' in data:
                coh = data['coh']
                if isinstance(coh, np.ndarray) and coh.ndim > 0:
                    coh_values.append(coh)
                else:
                    # Use ones array of same shape as dvv if coh is invalid
                    coh_values.append(np.ones_like(dvv))
            else:
                # Use ones array of same shape as dvv if coh is missing
                coh_values.append(np.ones_like(dvv))
            
            dates.append(date)
            
            # Get additional metadata from first file
            if len(dates) == 1:
                # Get taxis from first file
                if 'taxis' in data:
                    taxis = data['taxis']
                else:
                    logger.warning(f"No taxis found in {file}, will use default")
                    taxis = np.linspace(-120, 120, 241)  # Default taxis
                    
                # Override metadata if available in the file
                if 'filterid' in data:
                    try:
                        filterid_val = data['filterid']
                        if isinstance(filterid_val, np.ndarray) and filterid_val.size == 1:
                            filterid = int(filterid_val.item())
                        else:
                            filterid = int(filterid_val)
                    except:
                        filterid = int(filter_dir)
                        
                if 'component' in data:
                    try:
                        component_val = data['component']
                        if isinstance(component_val, np.ndarray) and component_val.size == 1:
                            component = str(component_val.item())
                        else:
                            component = str(component_val)
                    except:
                        component = component_dir
        
        except Exception as e:
            logger.error(f"Error loading file {file}: {str(e)}")
            continue
    
    if not dates or len(dates) != len(dvv_values):
        logger.warning(f"No valid data found or data mismatch for {pair_dir} - {component_dir} f{filter_dir} m{movstack_dir}")
        return False
    
    # Get frequencies from the first file
    try:
        first_file = np.load(files[0], allow_pickle=True)
        if 'freqs' in first_file:
            freqs_data = first_file['freqs']
            if isinstance(freqs_data, np.ndarray) and freqs_data.size > 0:
                freqs = freqs_data
            else:
                # Generate default frequencies if invalid
                logger.warning(f"Invalid freqs data in {files[0]}, using defaults")
                freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
        else:
            # Generate default frequencies if missing
            logger.warning(f"No freqs data in {files[0]}, using defaults")
            freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
    except Exception as e:
        logger.error(f"Error extracting frequencies: {str(e)}")
        # Generate default frequencies based on data shape
        freqs = np.linspace(0.1, 2.0, dvv_values[0].shape[0])
    
    # Create DataFrames with robust shape validation
    try:
        # Ensure all data arrays have the same shape
        common_shape = dvv_values[0].shape[0]
        valid_indices = []
        
        for i, (date, dvv, err, coh) in enumerate(zip(dates, dvv_values, err_values, coh_values)):
            if (dvv.shape[0] == common_shape and 
                err.shape[0] == common_shape and 
                coh.shape[0] == common_shape):
                valid_indices.append(i)
            else:
                logger.warning(f"Skipping data point {date} due to shape mismatch")
        
        if not valid_indices:
            logger.error(f"No valid data points with consistent shapes")
            return False
            
        # Filter to only use valid data points
        filtered_dates = [dates[i] for i in valid_indices]
        filtered_dvv = [dvv_values[i] for i in valid_indices]
        filtered_err = [err_values[i] for i in valid_indices]
        filtered_coh = [coh_values[i] for i in valid_indices]
        
        dvv_df = pd.DataFrame(filtered_dvv, index=filtered_dates, columns=freqs)
        err_df = pd.DataFrame(filtered_err, index=filtered_dates, columns=freqs)
        coh_df = pd.DataFrame(filtered_coh, index=filtered_dates, columns=freqs)
        
        # Sort by date
        dvv_df = dvv_df.sort_index()
        err_df = err_df.sort_index()
        coh_df = coh_df.sort_index()
        
        logger.info(f"Successfully created DataFrames with {len(filtered_dates)} time points and {len(freqs)} frequencies")
        
    except Exception as e:
        logger.error(f"Error creating DataFrames: {str(e)}")
        return False
    
    # Parse mov_stack format for MSNoise
    if '_' in movstack_dir:
        parts = movstack_dir.split('_')
        mov_stack = (parts[0], parts[1])
    else:
        mov_stack = (movstack_dir, movstack_dir)
    
    # Use api xr_save_wct function
    try:
        logger.info(f"Attempting to save using xr_save_wct function")
        xr_save_wct2(station1, station2, component, filterid, mov_stack, taxis, dvv_df, err_df, coh_df, output_dir=output_dir)
        logger.info(f"Successfully saved WCT data using xr_save_wct for {pair_dir}")
        return True
        
    except Exception as e:
        logger.error(f"Error with xr_save_wct: {str(e)}")
        logger.info("Attempting fallback to direct netCDF save")
        
        # Fallback to direct xarray save if MSNoise function fails
        try:
            # Create output directory matching the input hierarchy
            output_subdir = os.path.join(output_dir, filter_dir, movstack_dir, component_dir)
            os.makedirs(output_subdir, exist_ok=True)
            
            # Create xarray dataset
            coords = {
                'time': dvv_df.index,
                'freq': freqs
            }
            
            attrs = {
                'station1': station1,
                'station2': station2,
                'component': component,
                'filterid': filterid,
                'mov_stack': movstack_dir,
                'created': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                'merged_from': f"{len(filtered_dates)} daily files"
            }
            
            dvv_da = xr.DataArray(
                dvv_df.values,
                dims=['time', 'freq'],
                coords=coords,
                name='dvv'
            )
            
            err_da = xr.DataArray(
                err_df.values,
                dims=['time', 'freq'],
                coords=coords,
                name='err'
            )
            
            coh_da = xr.DataArray(
                coh_df.values,
                dims=['time', 'freq'],
                coords=coords,
                name='coh'
            )
            
            ds = xr.Dataset(
                {'dvv': dvv_da, 'err': err_da, 'coh': coh_da},
                attrs=attrs
            )
            
            # Save to file - use the same station pair name for the output file
            filename = f"{pair_dir}.nc"
            output_path = os.path.join(output_subdir, filename)
            
            ds.to_netcdf(output_path)
            logger.info(f"Saved merged file using fallback method: {output_path}")
            return True
            
        except Exception as e2:
            logger.error(f"Error with fallback save method: {str(e2)}")
            logger.debug(traceback.format_exc())
            return False

def main(wct_dir='WCT', output_dir='WCT_MERGED', loglevel="INFO", num_processes=1, use_tqdm=False):
    """
    Main function to merge all daily WCT files
    
    Parameters:
    -----------
    wct_dir : str
        Directory containing WCT results [default: WCT]
    output_dir : str
        Directory to save merged results [default: WCT_MERGED]
    loglevel : str
        Logging level [default: INFO]
    num_processes : int
        Number of parallel processes to use [default: 1 for MSNoise integration]
    use_tqdm : bool
        Use tqdm progress bar (False for MSNoise integration) [default: False]
    
    Returns:
    --------
    int
        0 if successful, 1 if error
    """
    # Configure logger using MSNoise API
    logger = get_logger('msnoise.merge_wct', loglevel)
    logger.info(f"Starting merge process, reading from {wct_dir}, writing to {output_dir}")
    
    # Set up progress tracking
    log_file = setup_progress_tracking(output_dir)
    completed_merges = load_completed_merges(log_file)
    logger.info(f"Found {len(completed_merges)} previously completed merges")
    
    try:
        # Get list of filter directories
        filters = get_filter_directories(wct_dir)
        if not filters:
            logger.error(f"No filter directories found in {wct_dir}")
            return 1
        
        logger.info(f"Found {len(filters)} filter directories to process")
        
        # Count total combinations
        total_combinations = 0
        for filter_dir in filters:
            try:
                movstacks = get_movstack_directories(wct_dir, filter_dir)
                for movstack_dir in movstacks:
                    try:
                        components = get_component_directories(wct_dir, filter_dir, movstack_dir)
                        for component_dir in components:
                            try:
                                pairs = get_station_pair_directories(wct_dir, filter_dir, movstack_dir, component_dir)
                                total_combinations += len(pairs)
                            except Exception as e:
                                logger.warning(f"Error accessing pairs in {filter_dir}/{movstack_dir}/{component_dir}: {e}")
                    except Exception as e:
                        logger.warning(f"Error accessing components in {filter_dir}/{movstack_dir}: {e}")
            except Exception as e:
                logger.warning(f"Error accessing movstacks in {filter_dir}: {e}")
        
        logger.info(f"Total combinations to process: {total_combinations}")
        
        if total_combinations == 0:
            logger.warning("No data combinations found to merge")
            return 0
        
        # Start timing
        start_time = datetime.datetime.now()
        
        # Process each parameter combination sequentially (better for MSNoise integration)
        merged_count = 0
        skipped_count = 0
        failed_count = 0
        
        processed = 0
        for filter_dir in filters:
            logger.info(f"Processing filter: {filter_dir}")
            
            try:
                # Get moving stack directories
                movstacks = get_movstack_directories(wct_dir, filter_dir)
                for movstack_dir in movstacks:
                    logger.info(f"Processing movstack: {movstack_dir}")
                    
                    try:
                        # Get component directories
                        components = get_component_directories(wct_dir, filter_dir, movstack_dir)
                        for component_dir in components:
                            logger.info(f"Processing component: {component_dir}")
                            
                            try:
                                # Get station pair directories
                                pairs = get_station_pair_directories(wct_dir, filter_dir, movstack_dir, component_dir)
                                for pair_dir in pairs:
                                    processed += 1
                                    
                                    # Check if this merge was already completed
                                    merge_key = create_merge_key(filter_dir, movstack_dir, component_dir, pair_dir)
                                    
                                    if merge_key in completed_merges:
                                        logger.debug(f"Skipping already completed merge: {merge_key}")
                                        skipped_count += 1
                                        continue
                                    
                                    logger.info(f"Merging files for filter {filter_dir}, movstack {movstack_dir}, component {component_dir}, pair {pair_dir} ({processed}/{total_combinations})")
                                    
                                    result = merge_daily_files(
                                        wct_dir, output_dir, filter_dir, 
                                        movstack_dir, component_dir, pair_dir, logger
                                    )
                                    
                                    if result:
                                        merged_count += 1
                                        # Log the completed merge
                                        save_completed_merge(log_file, merge_key)
                                        logger.info(f"Completed merge {merged_count}: {merge_key}")
                                    else:
                                        failed_count += 1
                                        logger.warning(f"Failed to merge: {merge_key}")
                            except Exception as e:
                                logger.error(f"Error processing pairs in {filter_dir}/{movstack_dir}/{component_dir}: {e}")
                                continue
                    except Exception as e:
                        logger.error(f"Error processing components in {filter_dir}/{movstack_dir}: {e}")
                        continue
            except Exception as e:
                logger.error(f"Error processing movstacks in {filter_dir}: {e}")
                continue
        
        # Calculate summary statistics
        elapsed = datetime.datetime.now() - start_time
        avg_time = elapsed.total_seconds() / total_combinations if total_combinations > 0 else 0
        
        logger.info(f"Merge process completed in {str(elapsed).split('.')[0]}")
        logger.info(f"Processed {total_combinations} combinations: {merged_count} merged, {skipped_count} skipped, {failed_count} failed")
        logger.info(f"Average processing time: {avg_time:.2f} seconds per combination")
        logger.info(f"Progress log: {log_file}")
        
        return 0
    except Exception as e:
        logger.error(f"Error in merge process: {str(e)}")
        logger.debug(traceback.format_exc())
        return 1

def dvv_merge_wct(ctx, wct_dir, output_dir):
    """
    MSNoise command function to merge WCT files
    
    Parameters:
    -----------
    ctx : click.Context
        Click context containing MSNoise configuration
    wct_dir : str
        Directory containing WCT results
    output_dir : str
        Directory to save merged results
    """
    loglevel = ctx.obj['MSNOISE_verbosity']
    threads = ctx.obj.get('MSNOISE_threads', 1)
    
    # For MSNoise integration, we use sequential processing and MSNoise logging
    # Multiprocessing can cause issues with MSNoise's database connections
    return main(wct_dir=wct_dir, output_dir=output_dir, loglevel=loglevel, 
                num_processes=1, use_tqdm=False)

if __name__ == "__main__":
    import sys
    sys.exit(main())
