"""
This step preprocesses waveforms using the preprocessing.py module and saves
the resulting Stream objects to disk in a workflow-aware folder structure.
"""

import os
import datetime
import pickle
import traceback
import numpy as np
from .api import (connect, get_logger, get_params, update_job,
                  preload_instrument_responses, get_workflow_steps, get_config_set_details, get_next_job_for_step)
from .preprocessing import preprocess
from obspy.core import AttribDict

def get_workflow_step_config(session, step_name):
    """Get workflow step configuration by step name."""
    steps = get_workflow_steps(session)
    for step in steps:
        if step.step_name == step_name:
            return step
    return None


def save_preprocessed_streams(stream, output_dir, step_name, goal_day):
    """
    Save preprocessed streams to disk in workflow-aware structure.

    :param stream: Dictionary {station: stream}
    :param output_dir: Base output directory
    :param step_name: Workflow step name
    :param goal_day: Processing date
    :return: List of saved file paths
    """
    # Create workflow-aware directory structure
    workflow_dir = os.path.join(output_dir, step_name)
    os.makedirs(workflow_dir, exist_ok=True)
    filename = f"{goal_day}.mseed"
    output_path = os.path.join(workflow_dir, filename)
    stream.write(output_path, format="MSEED")


    return [output_path]


def main(init=False, threads=1, loglevel="INFO"):
    """
    Main preprocessing workflow function.

    Processes waveforms using the preprocessing.py module and saves
    the resulting Stream objects to disk.
    """
    logger = get_logger('msnoise.step_preprocessing', loglevel, with_pid=True)
    logger.info('*** Starting: Preprocessing Step ***')

    # Connect to database and get params
    db = connect()
    params = get_params(db)

    # # Load instrument responses if needed
    # responses = None
    # if params.remove_response in ["Y", "y"]:
    #     logger.info("Loading instrument responses...")
    #     responses = preload_instrument_responses(db)

    # Get output directory
    output_dir = getattr(params, 'output_folder', 'PREPROCESSED')

    job_count = 0

    while True:
        # Get next set of preprocessing jobs (same step + day)
        jobs, step = get_next_job_for_step(db, step_category="preprocess")

        if not jobs:
            logger.info("No more preprocessing jobs to process")
            break


        # All jobs in the set have the same step and day
        first_job = jobs[0]



        step_name = first_job.jobtype
        goal_day = first_job.day

        logger.info(f"Processing {len(jobs)} jobs for step '{step_name}' on {goal_day}")

        try:
            # Get workflow step configuration
            # step_config = get_workflow_step_config(db, step_name)
            step_config = get_config_set_details(db, first_job.config_category, first_job.config_set_number, format='AttribDict')
            print(step_config)
            if not step_config:
                logger.error(f"No workflow step configuration found for: {step_name}")
                for job in jobs:
                    update_job(db, job.day, job.pair, job.jobtype, 'F')
                continue

            # Parse components from step configuration
            components = step_config.preprocess_components.split(',') if step_config.preprocess_components else ['Z']

            if step_config.remove_response:
                logger.debug('Pre-loading all instrument response')
                responses = preload_instrument_responses(db, return_format="inventory")
            else:
                responses = None
            # Extract stations from jobs
            stations = [job.pair for job in jobs]

            logger.info(f"Processing stations: {stations}")
            logger.info(f"Components: {components}")

            # Mark all jobs as in progress
            for job in jobs:
                update_job(db, job.day, job.pair, job.jobtype, 'I')


            stream = preprocess(stations, components, goal_day, AttribDict(**params, **step_config),
                                responses=responses, loglevel=loglevel)

            print(stream)
            ids = [f"{tr.stats.network}.{tr.stats.station}.{tr.stats.location}" for tr in stream]

            # Save all preprocessed streams
            saved_files = save_preprocessed_streams(stream, output_dir,
                                                    step_name, goal_day)

            logger.info(f"Saved {len(saved_files)} preprocessed files")

            # Update job statuses
            for job in jobs:
                station = job.pair
                if station in ids:
                    update_job(db, job.day, job.pair, job.jobtype, 'D')
                    logger.info(f"Job {job.jobtype} for {job.pair} marked as done")
                else:
                    update_job(db, job.day, job.pair, job.jobtype, 'F')
                    logger.warning(f"Job {job.jobtype} for {job.pair} marked as failed")

            job_count += len(jobs)

        except Exception as e:
            logger.error(f"Error processing job set for step {step_name} on {goal_day}: {str(e)}")
            logger.error(traceback.format_exc())

            # Mark all jobs as failed
            for job in jobs:
                update_job(db, job.day, job.pair, job.jobtype, 'F')

            continue
            break

    logger.info(f"*** Finished: Preprocessing Step - Processed {job_count} jobs ***")


if __name__ == "__main__":
    main()