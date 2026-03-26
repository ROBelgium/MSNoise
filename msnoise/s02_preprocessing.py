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

        goal_day  = days[0]
        step_name = step.step_name
        output_dir = getattr(params, 'output_folder', "OUTPUT")

        logger.info(f"Processing {len(jobs)} jobs for step '{step_name}' on {goal_day}")

        # Mark all in-progress BEFORE the try block
        for job in jobs:
            job.flag = 'I'
        db.commit()

        try:
            components = (params.preprocess_components.split(',')
                          if getattr(params, 'preprocess_components', None)
                          else ['Z'])

            if getattr(params, 'remove_response', 'N') in ('Y', 'y'):
                logger.debug('Pre-loading all instrument responses')
                responses = preload_instrument_responses(db, return_format="inventory")
            else:
                responses = None

            stations = [job.pair for job in jobs]
            logger.info(f"Processing stations: {stations}")
            logger.info(f"Components: {components}")

            stream = preprocess(stations, components, goal_day, params,
                                responses=responses, loglevel=loglevel)

            ids = {f"{tr.stats.network}.{tr.stats.station}.{tr.stats.location}"
                   for tr in stream}

            saved_files = save_preprocessed_streams(
                stream, output_dir, step_name, goal_day)
            logger.info(f"Saved {len(saved_files)} preprocessed files")

            # Mark per-job D or F depending on whether station was processed
            for job in jobs:
                if job.pair in ids:
                    job.flag = 'D'
                    logger.debug(f"Job {step_name} for {job.pair} marked Done")
                else:
                    job.flag = 'F'
                    logger.warning(f"Job {step_name} for {job.pair} marked Failed")
            db.commit()

        except Exception:
            logger.error(f"Error processing step {step_name} on {goal_day}:")
            logger.error(traceback.format_exc())
            for job in jobs:
                job.flag = 'F'
            db.commit()

    logger.info('*** Finished: Preprocessing Step ***')


if __name__ == "__main__":
    main()
