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

from .api import (
    connect,
    get_data_availability,
    get_logger,
    is_next_job_for_step,
    get_next_lineage_batch,
    massive_update_job,
    preload_instrument_responses,
    psd_ppsd_to_dataframe,
    to_sds,
    xr_save_psd,
)

CATEGORY = "psd"


def main(loglevel="INFO"):
    logger = get_logger(f"msnoise.{CATEGORY}", loglevel, with_pid=True)
    logger.info("*** Starting: Compute PSD ***")

    db = connect()
    # Always load the response, they are needed for the PSD calculations
    responses = preload_instrument_responses(db, return_format="inventory")

    while is_next_job_for_step(db, step_category=CATEGORY):
        batch = get_next_lineage_batch(db, step_category=CATEGORY,
                                       group_by="day", loglevel=loglevel,
                                       drop_current_step_name=False)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs          = batch["jobs"]
        step          = batch["step"]
        params        = batch["params"]
        lineage_names = batch["lineage_names"]
        days          = batch["days"]

        step_name     = step.step_name
        output_folder = getattr(params, "output_folder", "OUTPUT")

        psd_components = params.psd_components
        ppsd_length    = params.psd_ppsd_length
        ppsd_overlap   = params.psd_ppsd_overlap
        period_smooth  = params.psd_ppsd_period_smoothing_width_octaves
        period_step    = params.psd_ppsd_period_step_octaves
        period_limits  = params.psd_ppsd_period_limits
        db_bins        = params.psd_ppsd_db_bins

        # lineage_names already contains the current step name
        # (drop_current_step_name=False), so use it directly for path
        # construction once xr_save_psd accepts it; for now keep empty list
        # as psd is a root step (global -> psd).
        lineage = []

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
                logger.error(f"No files found for {job.pair} {job.day}")
                job.flag = "F"
                db.commit()
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
                    logger.info("Failed merging streams:")
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
                    logger.debug(f"PPSD.add failed for {job.pair} {comp}")
                    traceback.print_exc()
                    continue

                if not ppsd.times_processed:
                    logger.debug(
                        f"No PPSD windows for {job.pair} {comp} {job.day}"
                    )
                    del ppsd
                    continue

                df = psd_ppsd_to_dataframe(ppsd)
                seed_id = (f"{tr.stats.network}.{tr.stats.station}"
                           f".{tr.stats.location}.{tr.stats.channel}")

                try:
                    xr_save_psd(output_folder, lineage, step_name,
                                seed_id, job.day, df)
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

            job.flag = "F" if job_failed else "D"
            db.commit()
            logger.debug('Job done')

    logger.info("*** Finished: Compute PSD ***")
