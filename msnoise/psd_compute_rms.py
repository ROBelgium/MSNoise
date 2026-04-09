"""
Compute RMS from PSD NetCDF files and save results as NetCDF.

This step reads the per-day PSD NetCDF files written by compute_psd,
aggregates them per station, computes per-frequency-band RMS values, and
writes the results to a NetCDF file hierarchically below the upstream PSD
step:

    <output_folder>/<psd_step_name>/<psd_rms_step_name>/_output/<NET.STA.LOC.CHAN>/RMS.nc

The NetCDF file has two dimensions: times (one row per PPSD window)
and bands (one column per frequency band configured in
psd_rms_frequency_ranges).

.. versionchanged:: 2.1
    NC output replaces CSV.  Output path is now hierarchically below psd.
.. versionchanged:: 2.3
    Migrated from get_next_job_for_step to the canonical get_next_lineage_batch
    worker loop, consistent with all other worker scripts.  output_folder and
    step config are now sourced from the LayeredParams object returned by
    get_next_lineage_batch instead of raw get_config / get_config_set_details
    calls.
"""

import time
import traceback

import numpy as np
import pandas as pd

from .core.db import connect, get_logger
from .core.stations import get_station
from .core.workflow import get_next_lineage_batch, is_next_job_for_step, massive_update_job
from .core.io import psd_df_rms, xr_load_psd, xr_save_rms

CATEGORY = "psd_rms"


def main(loglevel="INFO", njobs_per_worker=9999):
    logger = get_logger(f"msnoise.{CATEGORY}", loglevel, with_pid=True)
    logger.info("*** Starting: Compute PSD RMS ***")

    db = connect()

    while is_next_job_for_step(db, step_category=CATEGORY):
        logger.debug("Getting the next batch")
        batch = get_next_lineage_batch(
            db, step_category=CATEGORY, group_by="pair_lineage",
            loglevel=loglevel,
        )
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs   = batch["jobs"]
        step   = batch["step"]
        params = batch["params"]

        output_folder = params.global_.output_folder
        step_name     = step.step_name

        rms_freq_ranges = params.psd_rms.psd_rms_frequency_ranges
        rms_type        = params.psd_rms.psd_rms_type

        first_job = jobs[0]
        net, sta, loc = first_job.pair.split(".")
        if loc == "--":
            loc = ""

        # With the corrected lineage "psd_1/psd_rms_1", lineage_names_upstream
        # is ["psd_1"] — i.e. it correctly ends with the upstream psd step name.
        # PSD NC files are written by psd_compute.py with lineage=[], so
        # we pass psd_lineage=[] to xr_load_psd.  The psd step name is
        # the last element of lineage_names_upstream (e.g. "psd_1").
        lineage_upstream = batch["lineage_names_upstream"]
        psd_step_name    = lineage_upstream[-1] if lineage_upstream else "psd_1"
        psd_lineage      = []   # PSD files live directly under output_folder/<step_name>/

        # RMS files live one level higher: output_folder/psd_1/psd_rms_1/_output/
        rms_lineage = [psd_step_name]

        station = get_station(db, net, sta)
        psd_components = params.psd.psd_components
        if psd_components:
            channels = [ch for ch in station.chans()
                        if ch[-1] in psd_components]
        else:
            channels = list(station.chans())

        days = sorted({job.day for job in jobs})

        logger.info(
            f"New PSD_RMS batch: {net}.{sta}.{loc} "
            f"n_days={len(days)} lineage={batch['lineage_str']}"
        )

        for chan in channels:
            seed_id = f"{net}.{sta}.{loc}.{chan}"
            logger.debug(f"Processing {seed_id}")

            frames = []
            for day in days:
                df = xr_load_psd(
                    output_folder, psd_lineage, psd_step_name, seed_id, day,
                    format="dataframe",
                )
                if df is not None and not df.empty:
                    frames.append(df)
                else:
                    logger.debug(f"No PSD NC for {seed_id} {day}")

            if not frames:
                logger.warning(f"No PSD data found for {seed_id}")
                continue

            data = pd.concat(frames).sort_index().sort_index(axis=1)

            try:
                rms = psd_df_rms(data, freqs=rms_freq_ranges, output=rms_type)
            except Exception:
                logger.error(f"psd_df_rms failed for {seed_id}")
                traceback.print_exc()
                continue

            try:
                xr_save_rms(output_folder, rms_lineage, step_name, seed_id, rms)
                logger.info(
                    f"Saved RMS NC for {seed_id} under "
                    f"{psd_step_name}/{step_name}/"
                )
            except Exception:
                logger.error(f"Failed saving RMS NC for {seed_id}")
                traceback.print_exc()

            del data, rms

        massive_update_job(db, jobs, "D")
        logger.debug("Batch done")

    logger.info("*** Finished: Compute PSD RMS ***")