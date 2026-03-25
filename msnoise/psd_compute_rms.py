"""
Compute RMS from PSD NetCDF files and export as CSV.

This step reads the per-day PSD NetCDF files written by ``compute_psd``,
aggregates them per station, computes per-frequency-band RMS values, and
writes the results to CSV files in the output folder.

The output path is::

    <output_folder>/<psd_rms_step_name>/_output/<NET.STA.LOC.CHAN>/RMS_<type>.csv

Each CSV has a DatetimeIndex (one row per PPSD window) and one column per
frequency band configured in ``psd_rms_frequency_ranges``.

Configuration Parameters (from ``psd_rms`` config set)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``psd_rms_frequency_ranges``
* ``psd_rms_type``

To run this step:

.. code-block:: sh

    $ msnoise qc compute_rms

.. versionadded:: 2.0
.. versionchanged:: 2.1
    Reads NC files instead of HDF stores.  The intermediate
    ``psd_to_hdf`` step is no longer needed.
"""

import datetime
import os
import traceback

import numpy as np
import pandas as pd

from .api import (
    connect,
    get_config_set_details,
    get_logger,
    get_next_job_for_step,
    get_params,
    get_station,
    is_next_job_for_step,
    massive_update_job,
    psd_df_rms,
    xr_load_psd,
)


def main(loglevel="INFO", njobs_per_worker=9999):
    logger = get_logger("msnoise.psd_compute_rms", loglevel, with_pid=True)
    logger.info("*** Starting: Compute PSD RMS ***")

    db = connect()
    orig_params = get_params(db)
    output_folder = getattr(orig_params, "output_folder", "OUTPUT")

    while is_next_job_for_step(db, step_category="psd_rms"):
        logger.info("Getting the next job")
        result = get_next_job_for_step(
            db, step_category="psd_rms", group_by="pair"
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

        # All jobs in the batch share the same station (group_by="pair")
        first_job = jobs[0]
        net, sta, loc = first_job.pair.split(".")

        # Derive the upstream psd step name from the job lineage.
        # The lineage field is set by the propagation function (step 6)
        # as e.g. "psd_1".  If not set, fall back to the step's
        # predecessor name from the workflow.
        psd_step_name = first_job.lineage or "psd_1"
        # Lineage list for xr_load_psd is empty (psd is a root step).
        psd_lineage   = []

        step_name = step.step_name  # e.g. "psd_rms_1"

        # Collect all days and channels in this batch
        station  = get_station(db, net, sta)
        if hasattr(step_config, "psd_components"):
            channels = [ch for ch in station.chans()
                        if ch[-1] in step_config.psd_components]
        else:
            channels = list(station.chans())

        days = sorted({job.day for job in jobs})

        for chan in channels:
            seed_id = f"{net}.{sta}.{loc}.{chan}"
            logger.debug(f"Processing {seed_id}")

            # Load all daily NC files for this channel
            frames = []
            for day in days:
                df = xr_load_psd(
                    output_folder, psd_lineage, psd_step_name, seed_id, day
                )
                if df is not None and not df.empty:
                    frames.append(df)
                else:
                    logger.debug(f"No PSD NC for {seed_id} {day}")

            if not frames:
                logger.warning(f"No PSD data found for {seed_id}")
                continue

            data = pd.concat(frames).sort_index()
            data = data.sort_index(axis=1)

            try:
                rms = psd_df_rms(data, freqs=rms_freq_ranges, output=rms_type)
            except Exception:
                logger.error(f"psd_df_rms failed for {seed_id}")
                traceback.print_exc()
                continue

            # Write CSV (append / update existing rows)
            out_dir = os.path.join(
                output_folder, step_name, "_output", seed_id
            )
            os.makedirs(out_dir, exist_ok=True)
            csv_path = os.path.join(out_dir, f"RMS_{rms_type}.csv")

            if os.path.isfile(csv_path):
                existing = pd.read_csv(csv_path, index_col=0, parse_dates=True)
                # Drop rows whose index is in the new batch then concat
                existing = existing[~existing.index.isin(rms.index)]
                rms = pd.concat([existing, rms]).sort_index()

            rms.to_csv(csv_path)
            logger.debug(f"Saved RMS CSV: {csv_path}")
            del data, rms

        massive_update_job(db, jobs, "D")
        logger.debug('Batch "D"one')

    logger.info("*** Finished: Compute PSD RMS ***")
