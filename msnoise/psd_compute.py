"""
Compute Power Spectral Densities and save results as NetCDF files.

This step processes **psd** workflow jobs — one job per station per day —
and writes a single NetCDF file per station-channel-day containing all
individual PSD window estimates.  The file layout is::

    <output_folder>/<psd_step_name>/_output/daily/<NET.STA.LOC.CHAN>/<YYYY-MM-DD>.nc

The NetCDF file has two dimensions: ``times`` (UTC datetimes of each PPSD
window used) and ``periods`` (period bin centres in seconds).  The single
data variable is ``PSD`` with units dB.

Concurrency is safe: each worker processes a different (station, day) job
selected atomically from the database, and each worker writes to a
distinct file path.

Configuration Parameters (from ``psd`` config set)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``psd_components``
* ``psd_ppsd_length``
* ``psd_ppsd_overlap``
* ``psd_ppsd_period_smoothing_width_octaves``
* ``psd_ppsd_period_step_octaves``
* ``psd_ppsd_period_limits``
* ``psd_ppsd_db_bins``

To run this step:

.. code-block:: sh

    $ msnoise qc compute_psd

Parallel execution:

.. code-block:: sh

    $ msnoise -t 4 qc compute_psd

.. versionadded:: 2.0
"""

import os
import time
import traceback

import matplotlib
matplotlib.use("Agg")
import numpy as np
from obspy.core import UTCDateTime
from obspy.signal import PPSD

from .core.db import connect, get_logger
from .core.stations import get_data_availability, get_station, resolve_data_source, read_waveforms_from_availability
from .core.fdsn import is_remote_source, fetch_raw_waveforms
from .core.workflow import get_next_lineage_batch, is_next_job_for_step, massive_update_job, propagate_downstream
from .core.signal import preload_instrument_responses
from .core.stations import to_sds
from .core.io import psd_ppsd_to_dataset, xr_save_psd

CATEGORY = "psd"


def main(loglevel="INFO", chunk_size=0):
    logger = get_logger(f"msnoise.{CATEGORY}", loglevel, with_pid=True)
    logger.info("*** Starting: Compute PSD ***")

    db = connect()
    responses = preload_instrument_responses(db, return_format="inventory")

    if chunk_size > 0:
        logger.info(f"PSD chunk_size={chunk_size}: each worker claims up to {chunk_size} stations per day")

    while is_next_job_for_step(db, step_category=CATEGORY):
        batch = get_next_lineage_batch(db, step_category=CATEGORY,
                                       group_by="day_lineage",
                                       chunk_size=chunk_size,
                                       loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs          = batch["jobs"]
        step          = batch["step"]
        params        = batch["params"]

        step_name     = step.step_name
        output_folder = params.global_.output_folder

        psd_components = params.psd.psd_components
        ppsd_length    = params.psd.psd_ppsd_length
        ppsd_overlap   = params.psd.psd_ppsd_overlap
        period_smooth  = params.psd.psd_ppsd_period_smoothing_width_octaves
        period_step    = params.psd.psd_ppsd_period_step_octaves
        period_limits  = params.psd.psd_ppsd_period_limits
        db_bins        = params.psd.psd_ppsd_db_bins

        # PSD is a root step (global → psd); there is no upstream lineage.
        lineage = []

        # ── Detect DataSource from the first job's station (all jobs in a
        # day_lineage batch share the same station/datasource) ────────────────
        first_net, first_sta, _ = jobs[0].pair.split(".")
        first_station = get_station(db, first_net, first_sta)
        ds = resolve_data_source(db, first_station)
        remote = is_remote_source(ds.uri)

        done_jobs   = []
        failed_jobs = []

        for job in jobs:
            net, sta, loc = job.pair.split(".")
            loc_clean = "" if loc == "--" else loc
            logger.debug(f"Processing {job.pair} {job.day}")

            t_start = UTCDateTime(job.day) - 1.5 * ppsd_length
            t_end   = UTCDateTime(job.day) + 86400
            gd      = UTCDateTime(job.day).datetime

            # ── Acquire raw waveforms ────────────────────────────────────────
            if remote:
                st = fetch_raw_waveforms(
                    db, [job], job.day, params,
                    t_start=t_start, t_end=t_end,
                )
                if not st:
                    logger.warning(
                        f"No FDSN data for {job.pair} {job.day} — marking Failed"
                    )
                    failed_jobs.append(job)
                    continue
                # Normalise "" location → "--"
                for tr in st:
                    if tr.stats.location == "":
                        tr.stats.location = "--"
            else:
                da_records = get_data_availability(
                    db,
                    net=net, sta=sta, loc=loc_clean,
                    starttime=t_start.datetime,
                    endtime=gd,
                )
                if not da_records:
                    logger.warning(
                        f"No files in DataAvailability for {job.pair} {job.day} "
                        f"— marking Failed"
                    )
                    failed_jobs.append(job)
                    continue
                st = read_waveforms_from_availability(
                    db, da_records, t_start, t_end, logger=logger
                )
                if not st:
                    logger.warning(
                        f"Could not read any waveforms for {job.pair} {job.day} "
                        f"— marking Failed"
                    )
                    failed_jobs.append(job)
                    continue

            # ── Normalise and merge ──────────────────────────────────────────
            try:
                st.merge()
            except Exception:
                logger.warning(f"Stream merge failed for {job.pair} {job.day}:")
                traceback.print_exc()
                failed_jobs.append(job)
                continue

            st = st.split()
            for tr in st:
                tr.stats.network = tr.stats.network.upper()
                tr.stats.station = tr.stats.station.upper()
                tr.stats.channel = tr.stats.channel.upper()

            # ── Per-component PPSD ───────────────────────────────────────────
            job_failed = False
            for comp in psd_components:
                sel = st.select(component=comp)
                if not sel:
                    continue
                tr = sel[0]

                ppsd = PPSD(
                    tr.stats,
                    metadata=responses,
                    ppsd_length=ppsd_length,
                    overlap=ppsd_overlap,
                    period_smoothing_width_octaves=period_smooth,
                    period_step_octaves=period_step,
                    period_limits=period_limits,
                    db_bins=db_bins,
                )
                try:
                    ppsd.add(st)
                except Exception:
                    logger.warning(f"PPSD.add failed for {job.pair} {comp}")
                    traceback.print_exc()
                    continue

                if not ppsd.times_processed:
                    logger.debug(f"No PPSD windows for {job.pair} {comp} {job.day}")
                    del ppsd
                    continue

                ds_nc = psd_ppsd_to_dataset(ppsd)
                seed_id = (f"{tr.stats.network}.{tr.stats.station}"
                           f".{tr.stats.location}.{tr.stats.channel}")

                try:
                    xr_save_psd(output_folder, lineage, step_name,
                                seed_id, job.day, ds_nc)
                    logger.debug(f"Saved PSD NC for {seed_id} {job.day}")
                except Exception:
                    logger.error(f"Failed saving PSD NC for {seed_id} {job.day}")
                    traceback.print_exc()
                    job_failed = True

                try:
                    out = to_sds(tr.stats, gd.year, int(gd.strftime("%j")))
                    pngout = os.path.join(output_folder, "PSD", "PNG", out)
                    os.makedirs(os.path.dirname(pngout), exist_ok=True)
                    ppsd.plot(pngout + ".png")
                except Exception:
                    logger.debug("Error saving PNG image")

                del ppsd

            del st  # release stream for this station/day before moving to next job
            if job_failed:
                failed_jobs.append(job)
            else:
                done_jobs.append(job)

        if done_jobs:
            massive_update_job(db, done_jobs, "D")
            logger.debug(f"Marked {len(done_jobs)} PSD job(s) Done")
        if failed_jobs:
            massive_update_job(db, failed_jobs, "F")
            logger.warning(f"Marked {len(failed_jobs)} PSD job(s) Failed")
        if done_jobs and not params.global_.hpc:
            propagate_downstream(db, batch)

    logger.info("*** Finished: Compute PSD ***")
