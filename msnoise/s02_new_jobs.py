from .msnoise_table_def import Lineage
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

import datetime

from .core.db import connect, get_logger
from .core.config import get_config, get_config_set_details
from .core.stations import get_new_files, get_stations, mark_data_availability
from .core.workflow import (build_movstack_datelist, _lineage_id_for, filter_within_daterange, get_lineages_to_step_id, get_workflow_steps, lineage_str_to_step_names, massive_insert_job, update_job)
import pandas as pd


# Module-level logger used by propagation functions (propagate_stack_jobs_from_cc_done, etc.)
# that are called before main() has a chance to configure a logger.
logger = get_logger('msnoise.new_jobs', loglevel="INFO")

def propagate_stack_jobs_from_cc_done(session):
    """
    From DONE CC jobs, create/update TODO MOV STACK jobs only:
      - day=<cc job day>  (normal per-day MOV stack recompute)

    REF jobs are now handled by ``refstack`` steps, triggered separately via
    :func:`propagate_refstack_jobs_from_stack_done` (``--after stack``).

    Jobs created are STACK jobs (jobtype = stack_step.step_name,
    step_id = stack_step.step_id) with lineage ending at that stack step
    (e.g. ``preprocess_1/cc_1/filter_1/stack_1``).
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.utcnow()
    created = 0
    bumped = 0

    # Be tolerant to naming (your UI shows "STACKING")
    cc_category_names = {"cc", "crosscorrelation", "cross-correlation"}
    stack_category_names = {"stack", "stacking"}

    cc_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category.in_(sorted(cc_category_names)))
        .all()
    )
    stack_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category.in_(sorted(stack_category_names)))
        .all()
    )

    logger.info(f"[--after cc] active cc_steps={len(cc_steps)} active stack_steps={len(stack_steps)}")

    if not cc_steps:
        logger.warning("[--after cc] No active CC steps found. Check WorkflowStep.category values.")
        return 0
    if not stack_steps:
        logger.warning("[--after cc] No active STACK steps found. Check WorkflowStep.category values (stack vs stacking).")
        return 0

    def _pair_is_auto(pair_str):
        a, b = pair_str.split(":")
        return a == b

    def _lineage_allowed_for_pair(lineage_step_names, pair_str):
        """
        Enforce filter step applicability (CC/AC/SC) based on pair type.
        If no filter step exists in the lineage, allow.
        """
        is_auto = _pair_is_auto(pair_str)
        is_cross = not is_auto

        filter_name = None
        for name in reversed(lineage_step_names):
            if name.startswith("filter_"):
                filter_name = name
                break
        if not filter_name:
            return True

        filt_step = (
            session.query(WorkflowStep)
            .filter(WorkflowStep.step_name == filter_name)
            .first()
        )
        if filt_step is None:
            raise ValueError(f"Lineage references filter step '{filter_name}' which cannot be resolved")

        filt_cfg = get_config_set_details(session, filt_step.category, filt_step.set_number, format="AttribDict")
        allow_cc = bool(getattr(filt_cfg, "CC", False))
        allow_ac = bool(getattr(filt_cfg, "AC", False))
        allow_sc = bool(getattr(filt_cfg, "SC", False))

        if is_cross:
            return allow_cc
        return allow_ac or allow_sc

    def _upsert_job(day_value, pair_value, lineage_to_stack, stack_step):
        nonlocal created, bumped

        existing = (
            session.query(Job)
            .filter(Job.step_id == stack_step.step_id)
            .filter(Job.jobtype == stack_step.step_name)
            .filter(Job.day == day_value)
            .filter(Job.pair == pair_value)
            .filter(Job.lineage_id == _lineage_id_for(session, lineage_to_stack))
            .first()
        )
        if existing is not None:
            if existing.flag != "T":
                existing.flag = "T"
                existing.lastmod = now
                bumped += 1
            return

        session.add(
            Job(
                day=day_value,
                pair=pair_value,
                flag="T",
                step_id=stack_step.step_id,
                jobtype=stack_step.step_name,
                priority=getattr(stack_step, "priority", 0) or 0,
                lastmod=now,
                lineage=lineage_to_stack,
            )
        )
        created += 1

    # For each CC step, map DONE CC jobs to all matching stack lineages
    for cc_step in cc_steps:
        done_cc_jobs = (
            session.query(Job)
            .filter(Job.step_id == cc_step.step_id)
            .filter(Job.flag == "D")
            .all()
        )
        logger.info(f"[--after cc] cc_step={cc_step.step_name} done_cc_jobs={len(done_cc_jobs)}")
        if not done_cc_jobs:
            continue

        for pj in done_cc_jobs:
            if not pj.lineage:
                raise ValueError("DONE CC job has empty lineage (v2 assumption)")

            parent_names = [n for n in lineage_str_to_step_names(pj.lineage, sep="/") if "global" not in n]
            if not parent_names:
                continue

            for stack_step in stack_steps:
                # enumerate full paths to this stack step (graph truth)
                stack_paths = get_lineages_to_step_id(session, step_id=stack_step.step_id, include_self=True)

                for path in stack_paths:
                    lin_names = [s.step_name for s in path if s.step_name and "global" not in s.step_name]

                    # must match the CC job lineage prefix
                    if len(lin_names) < len(parent_names):
                        continue
                    if lin_names[:len(parent_names)] != parent_names:
                        continue

                    # apply filter gating (CC/AC/SC) based on pair type
                    if not _lineage_allowed_for_pair(lin_names, pj.pair):
                        continue

                    lineage_to_stack = "/".join(lin_names)  # ends with stack_*

                    # MOV stack day job (same day as CC job)
                    # Note: REF jobs are now created by propagate_refstack_jobs_from_stack_done
                    _upsert_job(day_value=pj.day, pair_value=pj.pair, lineage_to_stack=lineage_to_stack, stack_step=stack_step)

    session.commit()
    logger.info(f"[--after cc] MOV STACK jobs: created={created}, bumped_to_T={bumped}")
    return created + bumped


def propagate_refstack_jobs_from_stack_done(session):
    """
    From DONE MOV stack jobs, create TODO refstack jobs with day="REF".

    For each DONE stack_N job (any day, any pair), find all refstack_M
    steps that are direct successors of stack_N in the workflow graph and
    create a single day="REF" job per (pair, refstack_M, lineage) if one
    does not already exist.

    Called via msnoise new_jobs --after stack.
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.utcnow()
    created = 0
    bumped = 0

    stack_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category == "stack")
        .all()
    )

    if not stack_steps:
        logger.warning("[--after stack] No active stack steps found.")
        return 0

    for stack_step in stack_steps:
        refstack_successors = (
            session.query(WorkflowStep)
            .join(schema.WorkflowLink, schema.WorkflowLink.to_step_id == WorkflowStep.step_id)
            .filter(schema.WorkflowLink.from_step_id == stack_step.step_id)
            .filter(schema.WorkflowLink.is_active.is_(True))
            .filter(WorkflowStep.category == "refstack")
            .all()
        )
        if not refstack_successors:
            continue

        done_pairs = (
            session.query(Job.pair)
            .filter(Job.step_id == stack_step.step_id)
            .filter(Job.flag == "D")
            .filter(Job.day != "REF")
            .distinct()
            .all()
        )
        if not done_pairs:
            continue

        for (pair,) in done_pairs:
            sample_stack_job = (
                session.query(Job)
                .filter(Job.step_id == stack_step.step_id)
                .filter(Job.pair == pair)
                .filter(Job.flag == "D")
                .filter(Job.day != "REF")
                .first()
            )
            if not sample_stack_job or not sample_stack_job.lineage:
                continue

            stack_lineage = sample_stack_job.lineage

            for ref_step in refstack_successors:
                refstack_lineage = stack_lineage + "/" + ref_step.step_name

                existing = (
                    session.query(Job)
                    .filter(Job.step_id == ref_step.step_id)
                    .filter(Job.day == "REF")
                    .filter(Job.pair == pair)
                    .filter(Job.lineage_id == _lineage_id_for(session, refstack_lineage))
                    .first()
                )
                if existing is not None:
                    if existing.flag != "T":
                        existing.flag = "T"
                        existing.lastmod = now
                        bumped += 1
                    continue

                session.add(Job(
                    day="REF",
                    pair=pair,
                    flag="T",
                    step_id=ref_step.step_id,
                    jobtype=ref_step.step_name,
                    priority=0,
                    lastmod=now,
                    lineage=refstack_lineage,
                ))
                created += 1

    session.commit()
    logger.info(f"[--after stack] REFSTACK jobs: created={created}, bumped_to_T={bumped}")
    return created + bumped


def propagate_mwcs_jobs_from_refstack_done(session):
    """
    From DONE refstack jobs (day="REF"), create TODO mwcs/stretching/wavelet jobs.

    For each DONE refstack_M REF job, find all downstream dv/v steps
    (mwcs, stretching, wavelet) and create one job per (pair, day) found among
    the completed MOV stack jobs that are upstream of this refstack.

    Called via msnoise new_jobs --after refstack.
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.utcnow()
    created = 0

    refstack_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category == "refstack")
        .all()
    )

    if not refstack_steps:
        logger.warning("[--after refstack] No active refstack steps found.")
        return 0

    for ref_step in refstack_steps:
        dvv_successors = (
            session.query(WorkflowStep)
            .join(schema.WorkflowLink, schema.WorkflowLink.to_step_id == WorkflowStep.step_id)
            .filter(schema.WorkflowLink.from_step_id == ref_step.step_id)
            .filter(schema.WorkflowLink.is_active.is_(True))
            .filter(WorkflowStep.category.in_(["mwcs", "stretching", "wavelet"]))
            .all()
        )
        if not dvv_successors:
            continue

        stack_predecessors = (
            session.query(WorkflowStep)
            .join(schema.WorkflowLink, schema.WorkflowLink.from_step_id == WorkflowStep.step_id)
            .filter(schema.WorkflowLink.to_step_id == ref_step.step_id)
            .filter(schema.WorkflowLink.is_active.is_(True))
            .filter(WorkflowStep.category == "stack")
            .all()
        )
        if not stack_predecessors:
            continue

        for stack_step in stack_predecessors:
            done_ref_jobs = (
                session.query(Job)
                .filter(Job.step_id == ref_step.step_id)
                .filter(Job.flag == "D")
                .filter(Job.day == "REF")
                .all()
            )
            if not done_ref_jobs:
                continue

            for ref_job in done_ref_jobs:
                pair = ref_job.pair
                refstack_lineage = ref_job.lineage

                done_stack_days = (
                    session.query(Job.day)
                    .filter(Job.step_id == stack_step.step_id)
                    .filter(Job.pair == pair)
                    .filter(Job.flag == "D")
                    .filter(Job.day != "REF")
                    .distinct()
                    .all()
                )

                for (day,) in done_stack_days:
                    for dvv_step in dvv_successors:
                        dvv_lineage = refstack_lineage + "/" + dvv_step.step_name

                        existing = (
                            session.query(Job.ref)
                            .filter(Job.step_id == dvv_step.step_id)
                            .filter(Job.day == day)
                            .filter(Job.pair == pair)
                            .filter(Job.lineage_id == _lineage_id_for(session, dvv_lineage))
                            .first()
                        )
                        if existing:
                            continue

                        session.add(Job(
                            day=day,
                            pair=pair,
                            flag="T",
                            step_id=dvv_step.step_id,
                            jobtype=dvv_step.step_name,
                            priority=0,
                            lastmod=now,
                            lineage=dvv_lineage,
                        ))
                        created += 1

    session.commit()
    logger.info(f"[--after refstack] DVV jobs created: {created}")
    return created


def propagate_dvv_jobs_from_dtt_done(session, source_category: str) -> int:
    """Insert one ``day="DVV"`` job per unique lineage into each downstream
    DVV step when the parent DTT step has at least one DONE job.

    Unlike the per-day passthrough used for MWCS/stretching/WCT jobs,
    DVV aggregation operates over the **full time series** of all pairs at
    once, so a single sentinel job per lineage is sufficient.

    :param source_category: ``"mwcs_dtt"``, ``"stretching"``, or
        ``"wct_dtt"`` — the parent DTT step category.
    :returns: Number of new DVV jobs created or reset to ``"T"``.
    """
    DVV_TARGET = {
        "mwcs_dtt":    "mwcs_dtt_dvv",
        "stretching":  "stretching_dvv",
        "wavelet_dtt": "wavelet_dtt_dvv",
    }
    target_category = DVV_TARGET.get(source_category)
    if target_category is None:
        raise ValueError(f"Unknown source_category {source_category!r}")

    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.utcnow()
    created = 0

    # Find all active parent DTT steps
    parent_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category == source_category)
        .all()
    )

    for parent_step in parent_steps:
        # Find downstream DVV steps linked from this parent
        dvv_steps = (
            session.query(WorkflowStep)
            .join(schema.WorkflowLink,
                  schema.WorkflowLink.to_step_id == WorkflowStep.step_id)
            .filter(schema.WorkflowLink.from_step_id == parent_step.step_id)
            .filter(schema.WorkflowLink.is_active.is_(True))
            .filter(WorkflowStep.category == target_category)
            .all()
        )
        if not dvv_steps:
            continue

        # Collect unique lineages from DONE parent jobs (any pair, any day)
        done_lineages = (
            session.query(Lineage.lineage_str)
            .join(Job, Job.lineage_id == Lineage.lineage_id)
            .filter(Job.step_id == parent_step.step_id)
            .filter(Job.flag == "D")
            .distinct()
            .all()
        )
        done_lineages = [row[0] for row in done_lineages if row[0]]
        if not done_lineages:
            continue

        for dvv_step in dvv_steps:
            for parent_lineage_str in done_lineages:
                # DVV lineage = parent lineage + "/" + dvv step name
                dvv_lineage_str = (
                    parent_lineage_str.rstrip("/") + "/" + dvv_step.step_name
                )

                # If a job already exists, reset to TODO if not already TODO
                existing = (
                    session.query(Job)
                    .filter(Job.step_id == dvv_step.step_id)
                    .filter(Job.day == "DVV")
                    .filter(Job.lineage == dvv_lineage_str)
                    .first()
                )
                if existing:
                    if existing.flag != "T":
                        existing.flag = "T"
                        existing.lastmod = now
                    continue

                session.add(Job(
                    day="DVV",
                    pair="ALL",
                    flag="T",
                    step_id=dvv_step.step_id,
                    jobtype=dvv_step.step_name,
                    priority=0,
                    lastmod=now,
                    lineage=dvv_lineage_str,
                ))
                created += 1

    session.commit()
    logger.info(
        f"[--after {source_category}] DVV aggregate jobs created: {created}"
    )
    return created


def create_passthrough_jobs_from_done_parent(session, parent_step, child_step):
    """
    DONE parent_step jobs -> TODO child_step jobs, preserving (day, pair).

    This is the "HPC passthrough" rule used from CC downstream in your setup.
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job

    now = datetime.datetime.utcnow()
    created = 0

    done_parent_jobs = (
        session.query(Job)
        .filter(Job.step_id == parent_step.step_id)
        .filter(Job.flag == "D")
        .filter(Job.day != "REF")
        .filter(Job.day != "DVV")
        .all()
    )

    def _pair_is_auto(pair_str):
        a, b = pair_str.split(":")
        return a == b

    def _lineage_allowed_for_pair(lineage_step_names, pair_str):
        """
        Enforce filter step applicability (CC/AC/SC) based on pair type.
        If no filter step exists in the lineage, allow.
        """
        is_auto = _pair_is_auto(pair_str)
        is_cross = not is_auto

        # Find last filter_* in the lineage
        filter_name = None
        for name in reversed(lineage_step_names):
            if name.startswith("filter_"):
                filter_name = name
                break
        if not filter_name:
            return True

        filt_step = (
            session.query(schema.WorkflowStep)
            .filter(schema.WorkflowStep.step_name == filter_name)
            .first()
        )
        if filt_step is None:
            raise ValueError(f"Lineage references filter step '{filter_name}' which cannot be resolved")

        filt_cfg = get_config_set_details(session, filt_step.category, filt_step.set_number, format="AttribDict")

        allow_cc = bool(getattr(filt_cfg, "CC", False))
        allow_ac = bool(getattr(filt_cfg, "AC", False))
        allow_sc = bool(getattr(filt_cfg, "SC", False))

        if is_cross:
            return allow_cc
        return allow_ac or allow_sc

    def build_child_lineages(parent_lineage_str, pair_str):
        """
        Return *all* valid child lineage strings for this parent lineage and pair,
        including intermediate steps (e.g., filter_1) but WITHOUT duplicating nodes.

        Key rule: select full graph paths to `child_step` whose prefix matches
        the parent job lineage prefix.
        """
        parent_names = [
            n for n in lineage_str_to_step_names(parent_lineage_str, sep="/")
            if "global" not in n
        ]
        if not parent_names:
            raise ValueError("Parent job has empty lineage (v2 assumption)")

        # Enumerate all possible full paths to the child step (graph truth)
        child_paths = get_lineages_to_step_id(session, step_id=child_step.step_id, include_self=True)

        candidates = []
        for path in child_paths:
            # Convert to step_name list, drop global_* nodes
            lin_names = [s.step_name for s in path if s.step_name and "global" not in s.step_name]

            # Must start with the parent lineage (prefix match), otherwise it's a different branch
            if len(lin_names) < len(parent_names):
                continue
            if lin_names[:len(parent_names)] != parent_names:
                continue

            # This is already the FULL lineage to child, so do NOT do parent + suffix
            # (that’s how you got preprocess_1/cc_1/cc_1/...)
            if not _lineage_allowed_for_pair(lin_names, pair_str):
                continue

            candidates.append("/".join(lin_names))

        # de-dup while preserving order
        unique = []
        seen = set()
        for x in candidates:
            if x in seen:
                continue
            seen.add(x)
            unique.append(x)

        if not unique:
            raise ValueError(
                f"No valid lineage for passthrough {parent_step.step_name} -> {child_step.step_name} "
                f"from parent_lineage='{parent_lineage_str}' pair='{pair_str}'"
            )

        return unique

    for pj in done_parent_jobs:
        for child_lineage in build_child_lineages(pj.lineage, pj.pair):
            exists = (
                session.query(Job.ref)
                .filter(Job.step_id == child_step.step_id)
                .filter(Job.day == pj.day)
                .filter(Job.pair == pj.pair)
                .filter(Job.lineage_id == _lineage_id_for(session, child_lineage))
                .first()
            )
            if exists:
                continue

            session.add(
                Job(
                    day=pj.day,
                    pair=pj.pair,
                    flag="T",
                    step_id=child_step.step_id,
                    jobtype=child_step.step_name,
                    priority=getattr(child_step, "priority", 0) or 0,
                    lastmod=now,
                    lineage=child_lineage,
                )
            )
            created += 1

    return created


def propagate_psd_rms_jobs_from_psd_done(session):
    """Create TODO psd_rms jobs from DONE psd jobs.

    For each DONE ``psd_N`` job (station, day), find all active ``psd_rms``
    successor steps linked to that psd step and create a matching
    ``psd_rms_M`` job for the same (station, day) if one does not already
    exist.

    The job ``lineage`` field is set to the upstream psd step name
    (e.g. ``"psd_1"``), so that :func:`~msnoise.psd_compute_rms.main`
    knows where to look for NC files.

    Called via ``msnoise new_jobs --after psd``.

    :returns: Number of jobs created.
    :rtype: int
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.utcnow()
    created = 0

    psd_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category == "psd")
        .all()
    )

    if not psd_steps:
        logger.warning("[--after psd] No active psd steps found.")
        return 0

    for psd_step in psd_steps:
        # Find all psd_rms steps directly downstream of this psd step
        psd_rms_steps = (
            session.query(WorkflowStep)
            .join(schema.WorkflowLink,
                  schema.WorkflowLink.to_step_id == WorkflowStep.step_id)
            .filter(schema.WorkflowLink.from_step_id == psd_step.step_id)
            .filter(schema.WorkflowLink.is_active.is_(True))
            .filter(WorkflowStep.category == "psd_rms")
            .all()
        )
        if not psd_rms_steps:
            continue

        # All DONE psd jobs for this psd step
        done_psd_jobs = (
            session.query(Job)
            .filter(Job.step_id == psd_step.step_id)
            .filter(Job.flag == "D")
            .all()
        )
        if not done_psd_jobs:
            continue

        for psd_job in done_psd_jobs:
            for psd_rms_step in psd_rms_steps:
                # Lineage encodes the upstream psd step name so that
                # psd_compute_rms knows which NC output folder to read.
                lineage = psd_step.step_name

                existing = (
                    session.query(Job.ref)
                    .filter(Job.step_id == psd_rms_step.step_id)
                    .filter(Job.day == psd_job.day)
                    .filter(Job.pair == psd_job.pair)
                    .join(Job.lineage_ref)
                    .filter(Lineage.lineage_str == lineage)
                    .first()
                )
                if existing:
                    continue

                session.add(Job(
                    day=psd_job.day,
                    pair=psd_job.pair,
                    flag="T",
                    step_id=psd_rms_step.step_id,
                    jobtype=psd_rms_step.step_name,
                    priority=0,
                    lastmod=now,
                    lineage=lineage,
                ))
                created += 1

    session.commit()
    logger.info(f"[--after psd] PSD_RMS jobs created: {created}")
    return created


def propagate_first_runnable_from_category(session, source_category, skip_categories=None):
    """
    For each WorkflowStep in `source_category`:
      - find first runnable steps per branch (skipping categories like "filter")
      - create passthrough jobs for those runnable steps from DONE parent jobs

    Returns total number of created jobs.
    """
    if skip_categories is None:
        skip_categories = {"filter"}

    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    from .core.workflow import get_first_runnable_steps_per_branch
    WorkflowStep = schema.WorkflowStep

    parent_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.category == source_category)
        .filter(WorkflowStep.is_active.is_(True))
        .all()
    )

    total_created = 0
    for parent_step in parent_steps:
        runnable_children = get_first_runnable_steps_per_branch(
            session,
            source_step_id=parent_step.step_id,
            skip_categories=skip_categories,
        )
        for child_step in runnable_children:
            total_created += create_passthrough_jobs_from_done_parent(
                session=session,
                parent_step=parent_step,
                child_step=child_step,
            )

    session.commit()
    return total_created


def create_cc_jobs_from_preprocess(session):
    """
    Create CC jobs based on completed preprocess jobs.

    This function looks for preprocess jobs marked as 'D'one and creates
    corresponding CC jobs based on workflow links and CC step configurations.
    """
    from .msnoise_table_def import declare_tables
    from .core.config import get_config_set_details
    from .core.workflow import get_step_successors

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    logger.info('Creating CC jobs from completed preprocess jobs')

    # Get all preprocess steps
    preprocess_steps = session.query(WorkflowStep) \
        .filter(WorkflowStep.category == "preprocess") \
        .filter(WorkflowStep.is_active.is_(True)) \
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
                        # NEW: for CC jobs created from preprocess outputs, lineage begins at preprocess step
                        # and should include cc_step too.
                        # We keep it aligned with your folder tree: preprocess_X/cc_Y
                        lineage = f"{preprocess_step.step_name}/{cc_step.step_name}"
                        job_data = {
                            "day": day,
                            "pair": pair,
                            "jobtype": cc_step.step_name,
                            "step_id": cc_step.step_id,
                            "priority": getattr(cc_step, 'priority', 0),
                            "flag": "T",
                            "lastmod": now,
                            "category": cc_step.category,
                            "set_number": cc_step.set_number,
                            "lineage": lineage,  # <-- NEW
                        }

                        jobtxt = ''.join(str(x) for x in job_data.values())
                        if jobtxt not in crap_all_jobs_text:
                            all_jobs.append(job_data)
                            crap_all_jobs_text.append(jobtxt)
                            count += 1

                if has_single_station:
                    # Create single-station CC jobs
                    for job in day_jobs:
                        lineage = f"{preprocess_step.step_name}/{cc_step.step_name}"  # <-- NEW
                        job_data = {
                            "day": day,
                            "pair": f"{job.pair}:{job.pair}",  # Single station
                            "jobtype": cc_step.step_name,
                            "step_id": cc_step.step_id,
                            "priority": getattr(cc_step, 'priority', 0),
                            "flag": "T",
                            "lastmod": now,
                            "category": cc_step.category,
                            "set_number": cc_step.set_number,
                            "lineage": lineage,  # <-- NEW
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


def main(init=False, nocc=False, after=False):
    logger.info('*** Starting: New Jobs (Workflow-aware) ***')

    db = connect()
    # params = get_params(db)

    if after:
        source_category = str(after).strip().lower()

        # Optional validation: ensure it's a known config-set type/category
        allowed_categories = {
            "global", "preprocess", "psd", "psd_rms", "cc", "filter",
            "stack", "refstack", "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
            "stretching", "stretching_dvv",
            "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
        }
        if source_category not in allowed_categories:
            raise ValueError(
                f"Invalid --after value '{after}'. Expected a config set type/category, "
                f"one of: {sorted(allowed_categories)}"
            )

        if source_category == "cc":
            created = propagate_stack_jobs_from_cc_done(session=db)
            logger.info(f'Propagation from category "cc" created/updated {created} MOV STACK job(s)')
            return

        # preprocess→cc is a fan-out: one single-station job → many station-pair jobs.
        # propagate_first_runnable_from_category can't handle that (it assumes pair format
        # is preserved), so we use the dedicated function instead.
        if source_category == "preprocess":
            cc_jobs, created = create_cc_jobs_from_preprocess(session=db)
            if cc_jobs:
                for job in cc_jobs:
                    update_job(db, job['day'], job['pair'],
                                        job['jobtype'], job['flag'],
                                        step_id=job.get('step_id'),
                                        priority=job.get('priority', 0),
                                        lineage=job.get('lineage'),
                                        commit=False)
                db.commit()
            logger.info(f'Propagation from category "preprocess" created {created} CC job(s)')
            return

        if source_category == "psd":
            created = propagate_psd_rms_jobs_from_psd_done(session=db)
            logger.info(f'Propagation from category "psd" created {created} PSD_RMS job(s)')
            return

        if source_category == "stack":
            created = propagate_refstack_jobs_from_stack_done(session=db)
            logger.info(f'Propagation from category "stack" created/updated {created} REFSTACK job(s)')
            return

        if source_category == "refstack":
            created = propagate_mwcs_jobs_from_refstack_done(session=db)
            logger.info(f'Propagation from category "refstack" created {created} DVV job(s)')
            return

        if source_category in ("mwcs_dtt", "stretching", "wavelet_dtt"):
            created = propagate_dvv_jobs_from_dtt_done(session=db,
                                                        source_category=source_category)
            logger.info(
                f'Propagation from category "{source_category}" '
                f'created {created} DVV aggregate job(s)'
            )
            return

        created = propagate_first_runnable_from_category(
            session=db,
            source_category=source_category,
            skip_categories={"filter"},
        )
        logger.info(f"HPC propagation from category \"{source_category}\" created {created} job(s)")
        return

    logger.debug("Checking plugins' entry points")
    plugins = get_config(db, "plugins", category='global', set_number=1)
    extra_jobtypes_scan_archive = []
    extra_jobtypes_new_files = ["PSD"]
    if plugins:
        from importlib.metadata import entry_points
        plugins = plugins.split(",")
        for ep in entry_points(group='msnoise.plugins.jobtypes'):
            module_name = ep.value.split(":")[0].split(".")[0]
            if module_name in plugins:
                jobtypes = ep.load()()
                for jobtype in jobtypes:
                    if jobtype["after"] == "scan_archive":
                        extra_jobtypes_scan_archive.append(jobtype["name"])
                    elif jobtype["after"] == "new_files":
                        extra_jobtypes_new_files.append(jobtype["name"])

    # Get all workflow steps with category "preprocess"
    workflow_steps = get_workflow_steps(db)
    preprocess_steps = [step for step in workflow_steps if step.category in ["preprocess", "psd"]]

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
                        lineage = step.step_name  # <-- NEW: preprocess/qc jobs start lineage at themselves
                        job = {
                            "day": date.date().strftime("%Y-%m-%d"),
                            "pair": "%s.%s.%s" % (nf.net, nf.sta, nf.loc),
                            "jobtype": step.step_name,  # Use step name as jobtype
                            "step_id": step.step_id,
                            "priority": 0,
                            "flag": "T",
                            "lastmod": now,
                            "lineage": lineage,        # <-- NEW
                        }
                        jobtxt = ''.join(str(x) for x in job.values())
                        if jobtxt not in crap_all_jobs_text:
                            all_jobs.append(job)
                            crap_all_jobs_text.append(jobtxt)
                            count += 1

            if init and len(all_jobs) > 1e5:
                logger.debug('Already 100.000 jobs, inserting/updating')
                massive_insert_job(db, all_jobs)
                all_jobs = []
                count += 1e5


    if len(all_jobs) != 0:
        logger.debug('Inserting/Updating %i jobs' % len(all_jobs))
        if init:
            massive_insert_job(db, all_jobs)
        else:
            for job in all_jobs:
                update_job(db, job['day'], job['pair'],
                                    job['jobtype'], job['flag'],
                                    step_id=job.get('step_id'),
                                    priority=job.get('priority', 0),
                                    lineage=job.get('lineage'),  # <-- NEW
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
                massive_insert_job(db, cc_jobs)
            else:
                for job in cc_jobs:
                    update_job(db, job['day'], job['pair'],
                                        job['jobtype'], job['flag'],
                                        step_id=job.get('step_id'),
                                        priority=job.get('priority', 0),
                                        lineage=job.get('lineage'),  # <-- NEW
                                        commit=False)
            db.commit()
            logger.info(f"Created {cc_count} CC jobs")

    return count
