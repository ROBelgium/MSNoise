""" This script searches the database for files flagged "N"ew or "M"odified.
For each date in the configured range, it checks if other stations are
available and defines the new jobs to be processed.  Only jobs within the
configured ``startdate`` and ``enddate`` are considered avoiding unnecessary
job creation. Those are inserted in the
*jobs* table of the database.

To run it from the console:

.. code-block:: sh

    $ msnoise new_jobs

Upon first run, if you expect the number of jobs to be large (many days,
many stations), pass the ``--init`` parameter to optimize the insert. Only use
this flag once, otherwise problems will arise from duplicate entries in
the jobs table.

.. code-block:: sh

    $ msnoise new_jobs --init

Performance / running on HPC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By setting the ``hpc`` configuration parameter to ``Y``, you will disable the
automatic creation of jobs during the workflow, to avoid numerous
interactions with the database (select & update or insert). The jobs have
then to be inserted manually:

.. code-block:: sh

    $ msnoise new_jobs --hpc CC:STACK

should be run after the ``msnoise compute_cc`` step in order to create the
``STACK`` jobs.
"""

from .api import *
import pandas as pd


import logbook
logger = logbook.Logger(__name__)

def create_cc_jobs_from_preprocess(session, workflow_id="default"):
    """
    Create CC jobs based on completed preprocess jobs.

    This function looks for preprocess jobs marked as 'D'one and creates
    corresponding CC jobs based on workflow links and CC step configurations.
    """
    from .msnoise_table_def import declare_tables
    from .api import get_workflow_steps, get_step_successors, get_config_set_details

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    logger.info('Creating CC jobs from completed preprocess jobs')

    # Get all preprocess steps
    preprocess_steps = session.query(WorkflowStep) \
        .filter(WorkflowStep.workflow_id == workflow_id) \
        .filter(WorkflowStep.category == "preprocess") \
        .filter(WorkflowStep.is_active == True) \
        .all()

    all_jobs = []
    crap_all_jobs_text = []
    count = 0
    now = datetime.datetime.utcnow()

    for preprocess_step in preprocess_steps:
        logger.debug(f'Processing preprocess step: {preprocess_step.step_name}')

        # Get CC steps that follow this preprocess step
        cc_successor_steps = get_step_successors(session, preprocess_step.step_id)
        cc_steps = [step for step in cc_successor_steps if step.category == "cc"]

        if not cc_steps:
            logger.debug(f'No CC successors found for {preprocess_step.step_name}')
            continue

        # Get all completed preprocess jobs for this step
        completed_jobs = session.query(Job) \
            .filter(Job.jobtype == preprocess_step.step_name) \
            .filter(Job.workflow_id == workflow_id) \
            .filter(Job.step_id == preprocess_step.step_id) \
            .filter(Job.flag == 'D') \
            .all()

        if not completed_jobs:
            logger.debug(f'No completed jobs found for {preprocess_step.step_name}')
            continue

        # Group completed jobs by day
        jobs_by_day = {}
        for job in completed_jobs:
            day = job.day
            if day not in jobs_by_day:
                jobs_by_day[day] = []
            jobs_by_day[day].append(job)

        # Create CC jobs for each CC step
        for cc_step in cc_steps:
            logger.debug(f'Creating jobs for CC step: {cc_step.step_name}')

            # Get CC step configuration
            cc_config = get_config_set_details(session, cc_step.category,
                                               cc_step.set_number, format='AttribDict')

            if not cc_config:
                logger.error(f'No configuration found for CC step: {cc_step.step_name}')
                continue

            # Determine job type based on configuration
            has_cross_station = bool(getattr(cc_config, 'components_to_compute', None))
            has_single_station = bool(getattr(cc_config, 'components_to_compute_single_station', None))

            for day, day_jobs in jobs_by_day.items():
                if has_cross_station:
                    # Create cross-station CC jobs (pair-based)
                    stations = [job.pair for job in day_jobs]
                    pairs = create_station_pairs(stations)

                    for pair in pairs:
                        job_data = {
                            "day": day,
                            "pair": pair,
                            "jobtype": cc_step.step_name,
                            "workflow_id": cc_step.workflow_id,
                            "step_id": cc_step.step_id,
                            "priority": getattr(cc_step, 'priority', 0),
                            "flag": "T",
                            "lastmod": now,
                            "category": cc_step.category,
                            "set_number": cc_step.set_number
                        }

                        jobtxt = ''.join(str(x) for x in job_data.values())
                        if jobtxt not in crap_all_jobs_text:
                            all_jobs.append(job_data)
                            crap_all_jobs_text.append(jobtxt)
                            count += 1

                if has_single_station:
                    # Create single-station CC jobs
                    for job in day_jobs:
                        job_data = {
                            "day": day,
                            "pair": f"{job.pair}:{job.pair}",  # Single station
                            "jobtype": cc_step.step_name,
                            "workflow_id": cc_step.workflow_id,
                            "step_id": cc_step.step_id,
                            "priority": getattr(cc_step, 'priority', 0),
                            "flag": "T",
                            "lastmod": now,
                            "category": cc_step.category,
                            "set_number": cc_step.set_number
                        }

                        jobtxt = ''.join(str(x) for x in job_data.values())
                        if jobtxt not in crap_all_jobs_text:
                            all_jobs.append(job_data)
                            crap_all_jobs_text.append(jobtxt)
                            count += 1

    return all_jobs, count

def create_station_pairs(stations):
    """
    Create station pairs for cross-station correlation.

    :param stations: List of station identifiers (NET.STA.LOC format)
    :return: List of station pairs
    """
    pairs = []
    for i in range(len(stations)):
        for j in range(i + 1, len(stations)):
            pair = f"{stations[i]}:{stations[j]}"
            pairs.append(pair)
    return pairs


def main(init=False, nocc=False):
    logger.info('*** Starting: New Jobs (Workflow-aware) ***')

    db = connect()
    # params = get_params(db)
    logger.debug("Checking plugins' entry points")
    plugins = get_config(db, "plugins")
    extra_jobtypes_scan_archive = []
    extra_jobtypes_new_files = ["PSD"]
    if plugins:
        import pkg_resources
        plugins = plugins.split(",")
        for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.jobtypes'):
            module_name = ep.module_name.split(".")[0]
            if module_name in plugins:
                jobtypes = ep.load()()
                for jobtype in jobtypes:
                    if jobtype["after"] == "scan_archive":
                        extra_jobtypes_scan_archive.append(jobtype["name"])
                    elif jobtype["after"] == "new_files":
                        extra_jobtypes_new_files.append(jobtype["name"])

    # Get all workflow steps with category "preprocess"
    workflow_steps = get_workflow_steps(db, workflow_id="default")
    preprocess_steps = [step for step in workflow_steps if step.category in ["preprocess", "qc"]]

    logger.info(f'Found {len(preprocess_steps)} preprocess workflow steps')
    for step in preprocess_steps:
        logger.debug(f'  - {step.step_name} (ID: {step.step_id})')

    logger.info('Scanning New/Modified files')
    stations_to_analyse = []
    error = False
    for sta in get_stations(db, all=False):

        if not len(sta.locs()):
            logger.error("You haven't defined location codes to use for %s.%s, "
                         "you should run 'msnoise db update_loc_chan'; exiting." %
                         (sta.net, sta.sta))
            error = True
        for loc in sta.locs():
            stations_to_analyse.append("%s.%s.%s" % (sta.net, sta.sta, loc))
    if error:
        return

    all_jobs = []
    crap_all_jobs_text = []
    updated_days = []
    nfs = get_new_files(db)
    now = datetime.datetime.utcnow()
    start_date, end_date, datelist = build_movstack_datelist(db)
    count = 0
    # Create jobs for single-station workflow steps (PSD, etc.)
    for nf in nfs:
        tmp = "%s.%s.%s" % (nf.net, nf.sta, nf.loc)
        if tmp not in stations_to_analyse:
            continue

        start, end = nf.starttime.date(), nf.endtime.date()
        for date in pd.date_range(start, end, freq="D"):
            if filter_within_daterange(date.date(), start_date, end_date):
                updated_days.append(date.date())

                # Create jobs for single-station preprocess steps
                for step in preprocess_steps:
                    #todo add filter on component if nf.chan[-1] in step.
                    if 1:
                        job = {
                            "day": date.date().strftime("%Y-%m-%d"),
                            "pair": "%s.%s.%s" % (nf.net, nf.sta, nf.loc),
                            "jobtype": step.step_name,  # Use step name as jobtype
                            "workflow_id": step.workflow_id,
                            "step_id": step.step_id,
                            "priority": 0,
                            "flag": "T",
                            "lastmod": now
                        }
                        jobtxt = ''.join(str(x) for x in job.values())
                        if jobtxt not in crap_all_jobs_text:
                            all_jobs.append(job)
                            crap_all_jobs_text.append(jobtxt)
                            count += 1

            if init and len(all_jobs) > 1e5:
                logger.debug('Already 100.000 jobs, inserting/updating')
                massive_insert_job_workflow(db, all_jobs)
                all_jobs = []
                count += 1e5


    if len(all_jobs) != 0:
        logger.debug('Inserting/Updating %i jobs' % len(all_jobs))
        if init:
            massive_insert_job_workflow(db, all_jobs)
        else:
            for job in all_jobs:
                update_job_workflow(db, job['day'], job['pair'],
                                    job['jobtype'], job['flag'],
                                    workflow_id=job.get('workflow_id', 'default'),
                                    step_id=job.get('step_id'),
                                    priority=job.get('priority', 0),
                                    commit=False)
    db.commit()
    count += len(all_jobs)

    for sta in get_stations(db, all=True):
        mark_data_availability(db, sta.net, sta.sta, flag='A')

    db.commit()
    logger.info("Inserted %i jobs" % count)
    logger.info('*** Finished: New Jobs (Workflow-aware) ***')

    if not nocc:
        logger.info('Creating CC jobs from completed preprocess jobs')
        cc_jobs, cc_count = create_cc_jobs_from_preprocess(db)

        if cc_jobs:
            logger.debug(f'Inserting/Updating {len(cc_jobs)} CC jobs')
            if init:
                massive_insert_job_workflow(db, cc_jobs)
            else:
                for job in cc_jobs:
                    update_job_workflow(db, job['day'], job['pair'],
                                        job['jobtype'], job['flag'],
                                        workflow_id=job.get('workflow_id', 'default'),
                                        step_id=job.get('step_id'),
                                        priority=job.get('priority', 0),
                                        commit=False)
            db.commit()
            logger.info(f"Created {cc_count} CC jobs")

    return count


def update_job_workflow(session, day, pair, jobtype, flag, workflow_id="default",
                        step_id=None, priority=0, commit=True, returnjob=True, ref=None):
    """
    Updates or Inserts a new workflow-aware Job in the database.

    Extended version of update_job that handles workflow fields.
    """
    from sqlalchemy import text
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job

    if ref:
        job = session.query(Job).filter(text("ref=:ref")).params(ref=ref).first()
    else:
        job = session.query(Job) \
            .filter(text("day=:day")) \
            .filter(text("pair=:pair")) \
            .filter(text("jobtype=:jobtype")) \
            .filter(text("workflow_id=:workflow_id")) \
            .params(day=day, pair=pair, jobtype=jobtype, workflow_id=workflow_id).first()

    if job is None:
        # Create new job with workflow fields
        job = Job()
        job.day = day
        job.pair = pair
        job.jobtype = jobtype
        job.workflow_id = workflow_id
        job.step_id = step_id
        job.priority = priority
        job.flag = flag
        job.lastmod = datetime.datetime.utcnow()
        session.add(job)
    else:
        # Update existing job
        job.flag = flag
        job.workflow_id = workflow_id
        job.step_id = step_id
        job.priority = priority
        job.lastmod = datetime.datetime.utcnow()

    if commit:
        session.commit()
    if returnjob:
        return job


def massive_insert_job_workflow(session, jobs):
    """
    Massive insert of workflow-aware jobs using bulk operations.

    Extended version of massive_insert_job that handles workflow fields.
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job

    # Convert jobs to format expected by SQLAlchemy bulk operations
    job_records = []
    for job in jobs:
        job_records.append({
            'day': job['day'],
            'pair': job['pair'],
            'jobtype': job['jobtype'],
            'workflow_id': job.get('workflow_id', 'default'),
            'step_id': job.get('step_id'),
            'priority': job.get('priority', 0),
            'flag': job['flag'],
            'lastmod': job['lastmod']
        })

    # Use bulk insert for better performance
    session.bulk_insert_mappings(Job, job_records)
    session.commit()

