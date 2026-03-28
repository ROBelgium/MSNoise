"""
This step preprocesses waveforms using the preprocessing.py module and saves
the resulting Stream objects to disk in a workflow-aware folder structure.

For local/SDS sources, waveforms are read from the archive via
DataAvailability records.  For FDSN/EIDA sources, waveforms are fetched from
the remote service using ``get_waveforms_bulk`` and optionally cached as raw
files (``fdsn_keep_raw=Y``) before preprocessing.
"""

import time
import traceback
import numpy as np
from .core.db import connect, get_logger
from .core.workflow import get_next_lineage_batch, is_next_job_for_step, massive_update_job
from .core.signal import preload_instrument_responses, save_preprocessed_streams
from .core.stations import resolve_data_source, get_station
from .core.fdsn import is_remote_source, fetch_and_preprocess
from .preprocessing import preprocess

CATEGORY = "preprocess"


def main(loglevel="INFO"):
    """
    Main preprocessing workflow function.

    Dispatches to the local-archive or FDSN/EIDA fetch path depending on
    the station's DataSource URI scheme, then preprocesses and writes
    per-station output files.
    """
    logger = get_logger(f"msnoise.{CATEGORY}", loglevel, with_pid=True)
    logger.info('*** Starting: Preprocessing Step ***')

    db = connect()

    while is_next_job_for_step(db, step_category=CATEGORY):
        batch = get_next_lineage_batch(db, step_category=CATEGORY,
                                       group_by="day", loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs       = batch["jobs"]
        step       = batch["step"]
        params     = batch["params"]
        days       = batch["days"]

        goal_day   = days[0]
        step_name  = step.step_name
        output_dir = params.global_.output_folder

        logger.info(f"Processing {len(jobs)} jobs for step '{step_name}' on {goal_day}")

        # Mark all in-progress atomically before processing
        massive_update_job(db, jobs, "I")

        try:
            raw = params.preprocess.preprocess_components
            components = raw.split(',') if isinstance(raw, str) else (list(raw) if raw else ['Z'])

            if params.preprocess.remove_response in ('Y', 'y', True):
                logger.debug('Pre-loading all instrument responses')
                responses = preload_instrument_responses(db, return_format="inventory")
            else:
                responses = None

            # ── Detect DataSource scheme from the first job's station ────────
            first_net, first_sta, _ = jobs[0].pair.split(".")
            first_station = get_station(db, first_net, first_sta)
            ds = resolve_data_source(db, first_station)
            remote = is_remote_source(ds.uri)

            if remote:
                # ── FDSN / EIDA path ─────────────────────────────────────────
                logger.info(
                    f"DataSource {ds.name!r} is remote ({ds.uri}) — "
                    f"fetching via FDSN/EIDA"
                )
                stream, done_jobs, failed_jobs = fetch_and_preprocess(
                    db, jobs, goal_day, params,
                    responses=responses, loglevel=loglevel
                )
            else:
                # ── Local / SDS path ─────────────────────────────────────────
                stations = [job.pair for job in jobs]
                logger.debug(f"Processing stations (local): {stations}")
                logger.debug(f"Components: {components}")

                stream = preprocess(stations, components, goal_day, params,
                                    responses=responses, loglevel=loglevel)

                ids        = {f"{tr.stats.network}.{tr.stats.station}.{tr.stats.location}"
                              for tr in stream}
                done_jobs  = [j for j in jobs if j.pair in ids]
                failed_jobs = [j for j in jobs if j.pair not in ids]

            # ── Write per-station output files ───────────────────────────────
            if stream:
                saved_files = save_preprocessed_streams(
                    stream, output_dir, step_name, goal_day)
                logger.info(f"Saved {len(saved_files)} preprocessed file(s)")

            if done_jobs:
                massive_update_job(db, done_jobs, "D")
            if failed_jobs:
                logger.warning(
                    f"{len(failed_jobs)} job(s) for {step_name} on {goal_day} "
                    f"marked Failed (station not in stream)"
                )
                massive_update_job(db, failed_jobs, "F")

            if not batch["params"].global_.hpc:
                from .core.workflow import propagate_downstream
                propagate_downstream(db, batch)

        except Exception:
            logger.error(f"Error processing step {step_name} on {goal_day}:")
            logger.error(traceback.format_exc())
            massive_update_job(db, jobs, "F")

    logger.info('*** Finished: Preprocessing Step ***')


if __name__ == "__main__":
    main()
