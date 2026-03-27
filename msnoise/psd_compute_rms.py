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
"""

import traceback
import pandas as pd

from .db import connect, get_logger
from .config import get_config, get_config_set_details
from .stations import get_station
from .workflow import get_next_job_for_step, is_next_job_for_step, massive_update_job
from .io import psd_df_rms, xr_load_psd, xr_save_rms


def main(loglevel="INFO", njobs_per_worker=9999):
    logger = get_logger("msnoise.psd_compute_rms", loglevel, with_pid=True)
    logger.info("*** Starting: Compute PSD RMS ***")

    db = connect()
    output_folder = get_config(db, "output_folder") or "OUTPUT"

    while is_next_job_for_step(db, step_category="psd_rms"):
        logger.debug("Getting the next batch")
        result = get_next_job_for_step(
            db, step_category="psd_rms", group_by="pair_lineage"
        )

        if result is None:
            break
        jobs, step = result
        if not jobs:
            continue

        step_config = get_config_set_details(
            db,
            jobs[0].config_category,
            jobs[0].config_set_number,
            format="AttribDict",
        )

        rms_freq_ranges = step_config.psd_rms_frequency_ranges
        rms_type        = step_config.psd_rms_type

        first_job = jobs[0]
        net, sta, loc = first_job.pair.split(".")

        # job.lineage = upstream psd step name e.g. "psd_1"
        # Output hierarchy: OUTPUT/psd_1/psd_rms_1/_output/<seed_id>/RMS.nc
        psd_step_name = first_job.lineage or "psd_1"
        lineage       = [psd_step_name]
        step_name     = step.step_name

        # PSD NC files have no lineage above psd
        psd_lineage = []

        station = get_station(db, net, sta)
        if hasattr(step_config, "psd_components"):
            channels = [ch for ch in station.chans()
                        if ch[-1] in step_config.psd_components]
        else:
            channels = list(station.chans())

        days = sorted({job.day for job in jobs})

        for chan in channels:
            seed_id = f"{net}.{sta}.{loc}.{chan}"
            logger.debug(f"Processing {seed_id}")

            frames = []
            for day in days:
                df = xr_load_psd(
                    output_folder, psd_lineage, psd_step_name, seed_id, day,
                    format="dataframe"
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
                xr_save_rms(output_folder, lineage, step_name, seed_id, rms)
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