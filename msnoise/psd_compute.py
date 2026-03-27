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
.. versionchanged:: 2.1
    Output format changed from NPZ+PNG to NetCDF.  The intermediate
    ``psd_to_hdf`` step is no longer needed.
.. versionchanged:: 2.2
    Migrated to canonical get_next_lineage_batch worker loop.
    Instrument responses are only preloaded when remove_response is set.
    PNG output path is now relative to output_folder.
"""

import datetime
import os
import time
import traceback

import matplotlib
matplotlib.use("Agg")
import numpy as np
from obspy import Stream
from obspy.core import UTCDateTime
from obspy.signal import PPSD

from ...db import connect, get_logger
from ...stations import get_data_availability
from ...workflow import get_next_lineage_batch, is_next_job_for_step, massive_update_job
from ...signal import preload_instrument_responses, to_sds
from ...io import psd_ppsd_to_dataset, xr_save_psd

CATEGORY = "psd"


def main(loglevel="INFO", njobs_per_worker=9999):
    logger = get_logger(f"msnoise.{CATEGORY}", loglevel, with_pid=True)
    logger.info("*** Starting: Compute PSD ***")

    db = connect()
    # Always load the response, they are needed for the PSD calculations
    responses = preload_instrument_responses(db, return_format="inventory")

    while is_next_job_for_step(db, step_category=CATEGORY):
        batch = get_next_lineage_batch(db, step_category=CATEGORY,
                                       group_by="day", loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs          = batch["jobs"]
        step          = batch["step"]
        params        = batch["params"]

        step_name     = step.step_name
        output_folder = params.output_folder

        psd_components = params.psd.psd_components
        ppsd_length    = params.psd.psd_ppsd_length
        ppsd_overlap   = params.psd.psd_ppsd_overlap
        period_smooth  = params.psd.psd_ppsd_period_smoothing_width_octaves
        period_step    = params.psd.psd_ppsd_period_step_octaves
        period_limits  = params.psd.psd_ppsd_period_limits
        db_bins        = params.psd.psd_ppsd_db_bins

        # PSD is a root step (global → psd); there is no upstream lineage.
        # xr_save_psd / xr_load_psd use `lineage=[]` so files land at:
        #   <output_folder>/<step_name>/_output/daily/<seed_id>/<day>.nc
        # The batch lineage_names (e.g. ['global_1', 'psd_1']) includes
        # the global step which is not a real output folder — using [] is
        # consistent with how psd_compute_rms.py reads these files.
        lineage = []

        done_jobs   = []
        failed_jobs = []

        for job in jobs:
            net, sta, loc = job.pair.split(".")
            logger.debug(f"Processing {job.pair} {job.day}")
            gd = UTCDateTime(job.day).datetime

            files = get_data_availability(
                db,
                net=net, sta=sta, loc=loc,
                starttime=(UTCDateTime(job.day) - 1.5 * ppsd_length).datetime,
                endtime=gd,
            )
            if not files:
                logger.warning(f"No files found for {job.pair} {job.day} — marking Failed")
                failed_jobs.append(job)
                continue

            job_failed = False
            for comp in psd_components:
                toprocess = [
                    os.path.join(f.path, f.file)
                    for f in files
                    if f.chan[-1] == comp
                ]
                if not toprocess:
                    continue

                st = Stream()
                for fpath in np.unique(toprocess):
                    logger.debug(f"Reading {fpath}")
                    try:
                        st += __import__("obspy").read(
                            fpath,
                            starttime=UTCDateTime(gd) - 1.5 * ppsd_length,
                            endtime=UTCDateTime(
                                gd + datetime.timedelta(days=1)
                            ) - 0.001,
                        )
                    except Exception:
                        logger.debug(f"Problem loading {fpath}")

                if not st:
                    continue

                try:
                    st.merge()
                except Exception:
                    logger.warning("Failed merging streams:")
                    traceback.print_exc()
                    continue

                st = st.split()
                for tr in st:
                    tr.stats.network = tr.stats.network.upper()
                    tr.stats.station = tr.stats.station.upper()
                    tr.stats.channel = tr.stats.channel.upper()

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
                    logger.debug(
                        f"No PPSD windows for {job.pair} {comp} {job.day}"
                    )
                    del ppsd
                    continue

                ds = psd_ppsd_to_dataset(ppsd)
                seed_id = (f"{tr.stats.network}.{tr.stats.station}"
                           f".{tr.stats.location}.{tr.stats.channel}")

                try:
                    xr_save_psd(output_folder, lineage, step_name,
                                seed_id, job.day, ds)
                    logger.debug(f"Saved PSD NC for {seed_id} {job.day}")
                except Exception:
                    logger.error(f"Failed saving PSD NC for {seed_id} {job.day}")
                    traceback.print_exc()
                    job_failed = True

                # Optional PNG — path anchored under output_folder
                try:
                    out = to_sds(tr.stats, gd.year, int(gd.strftime("%j")))
                    pngout = os.path.join(output_folder, "PSD", "PNG", out)
                    os.makedirs(os.path.dirname(pngout), exist_ok=True)
                    ppsd.plot(pngout + ".png")
                except Exception:
                    logger.debug("Error saving PNG image")

                del ppsd

            if job_failed:
                failed_jobs.append(job)
            else:
                done_jobs.append(job)

        # Batch-update all job flags in two calls instead of N commits
        if done_jobs:
            massive_update_job(db, done_jobs, "D")
            logger.debug(f"Marked {len(done_jobs)} PSD job(s) Done")
        if failed_jobs:
            massive_update_job(db, failed_jobs, "F")
            logger.warning(f"Marked {len(failed_jobs)} PSD job(s) Failed")

    logger.info("*** Finished: Compute PSD ***")
