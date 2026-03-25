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
"""

import datetime
import os
import traceback

import matplotlib
matplotlib.use("Agg")
import numpy as np
from obspy import Stream
from obspy.core import AttribDict, UTCDateTime
from obspy.signal import PPSD

from .api import (
    connect,
    get_config,
    get_config_set_details,
    get_data_availability,
    get_logger,
    get_params,
    is_next_job_for_step,
    get_next_job_for_step,
    preload_instrument_responses,
    psd_ppsd_to_dataframe,
    to_sds,
    update_job,
    xr_save_psd,
)


def main(loglevel="INFO", njobs_per_worker=9999):
    logger = get_logger("msnoise.psd_compute", loglevel, with_pid=True)
    logger.info("*** Starting: Compute PSD ***")

    db = connect()
    logger.debug("Preloading all instrument responses")
    responses = preload_instrument_responses(db, return_format="inventory")

    orig_params = get_params(db)
    output_folder = getattr(orig_params, "output_folder", "OUTPUT")

    while is_next_job_for_step(db, step_category="psd"):
        logger.info("Getting the next job")
        result = get_next_job_for_step(db, step_category="psd")
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
        params = AttribDict(**orig_params, **step_config)

        psd_components   = params.psd_components
        ppsd_length      = params.psd_ppsd_length
        ppsd_overlap     = params.psd_ppsd_overlap
        period_smooth    = params.psd_ppsd_period_smoothing_width_octaves
        period_step      = params.psd_ppsd_period_step_octaves
        period_limits    = params.psd_ppsd_period_limits
        db_bins          = params.psd_ppsd_db_bins

        # Lineage for output path: upstream steps leading to this one.
        # For a psd_N step the lineage is just [step.step_name] since psd
        # is a root-level workflow (global → psd).
        lineage = []  # no upstream steps
        step_name = step.step_name  # e.g. "psd_1"

        for job in jobs:
            net, sta, loc = job.pair.split(".")
            logger.debug(f"Processing {job.pair}")
            gd = UTCDateTime(job.day).datetime

            files = get_data_availability(
                db,
                net=net, sta=sta, loc=loc,
                starttime=(UTCDateTime(job.day) - 1.5 * ppsd_length).datetime,
                endtime=gd,
            )
            if not files:
                logger.error(f"No files found for {job.pair} {job.day}")
                update_job(db, job.day, job.pair, job.jobtype, "F", ref=job.ref)
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

                # PPSD computation
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
                        f"No PPSD windows processed for {job.pair} {comp} {job.day}"
                    )
                    del ppsd
                    continue

                # Convert to DataFrame and save as NetCDF
                df = psd_ppsd_to_dataframe(ppsd)
                seed_id = f"{tr.stats.network}.{tr.stats.station}" \
                          f".{tr.stats.location}.{tr.stats.channel}"

                try:
                    xr_save_psd(
                        output_folder,
                        lineage,
                        step_name,
                        seed_id,
                        job.day,
                        df,
                    )
                    logger.debug(
                        f"Saved PSD NC for {seed_id} {job.day}"
                    )
                except Exception:
                    logger.error(
                        f"Failed saving PSD NC for {seed_id} {job.day}"
                    )
                    traceback.print_exc()
                    job_failed = True

                # Optional PNG for quick visual check
                try:
                    out = to_sds(tr.stats, gd.year, int(gd.strftime("%j")))
                    pngout = os.path.join("PSD", "PNG", out)
                    os.makedirs(os.path.dirname(pngout), exist_ok=True)
                    ppsd.plot(pngout + ".png")
                except Exception:
                    logger.debug("Error saving PNG image")

                del ppsd

            flag = "F" if job_failed else "D"
            update_job(db, job.day, job.pair, job.jobtype, flag, ref=job.ref)

        logger.debug('Day (job) "D"one')

    logger.info("*** Finished: Compute PSD ***")
