"""
This step preprocesses waveforms using the preprocessing.py module and saves
the resulting Stream objects to disk in a workflow-aware folder structure.
"""

import time
import traceback
import numpy as np
from .api import (connect, get_logger, is_next_job_for_step,
                  get_next_lineage_batch, massive_update_job,
                  preload_instrument_responses,
                  save_preprocessed_streams)
from .preprocessing import preprocess

CATEGORY = "preprocess"


def main(loglevel="INFO"):
    """
    Main preprocessing workflow function.

    Processes waveforms using the preprocessing.py module and saves
    the resulting Stream objects to disk.
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

        jobs          = batch["jobs"]
        step          = batch["step"]
        params        = batch["params"]
        days          = batch["days"]

        goal_day  = days[0]
        step_name = step.step_name
        output_dir = params.output_folder

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

            stations = [job.pair for job in jobs]
            logger.debug(f"Processing stations: {stations}")
            logger.debug(f"Components: {components}")

            stream = preprocess(stations, components, goal_day, params,
                                responses=responses, loglevel=loglevel)

            ids = {f"{tr.stats.network}.{tr.stats.station}.{tr.stats.location}"
                   for tr in stream}

            saved_files = save_preprocessed_streams(
                stream, output_dir, step_name, goal_day)
            logger.info(f"Saved {len(saved_files)} preprocessed files")

            # Partition jobs into done/failed and batch-update both
            done_jobs   = [j for j in jobs if j.pair in ids]
            failed_jobs = [j for j in jobs if j.pair not in ids]
            if done_jobs:
                massive_update_job(db, done_jobs, "D")
            if failed_jobs:
                logger.warning(f"{len(failed_jobs)} job(s) for {step_name} on {goal_day} marked Failed (station not in stream)")
                massive_update_job(db, failed_jobs, "F")

        except Exception:
            logger.error(f"Error processing step {step_name} on {goal_day}:")
            logger.error(traceback.format_exc())
            massive_update_job(db, jobs, "F")

    logger.info('*** Finished: Preprocessing Step ***')


if __name__ == "__main__":
    main()
