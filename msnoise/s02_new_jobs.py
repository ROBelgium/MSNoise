from .msnoise_table_def import Lineage
from sqlalchemy.orm import aliased
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
from .core.config import get_config, get_config_set_details, get_params
from .core.stations import get_new_files, get_stations, mark_data_availability
from .core.workflow import (build_movstack_datelist, _lineage_id_for, _get_or_create_lineage_id, filter_within_daterange, get_lineages_to_step_id, get_workflow_steps, lineage_str_to_step_names, massive_insert_job, update_job)
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

    Performance: fully batched — one SELECT to load all DONE CC jobs, one
    bulk SELECT to find existing STACK jobs, then bulk_insert_mappings +
    a single UPDATE … WHERE ref IN (…) for bumps.  No per-job DB round-trips.
    """
    from sqlalchemy import tuple_
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.now(datetime.timezone.utc)
    created = 0
    bumped = 0

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

    # --- Hoist invariant per-filter-step queries out of the pair loop ----------
    # Cache: filter_step_name -> (allow_cross, allow_auto)
    _filter_allowed_cache: dict = {}

    def _lineage_allowed_for_pair(lineage_step_names, is_auto):
        """
        Enforce filter step applicability (CC/AC/SC) based on pair type.
        If no filter step exists in the lineage, allow.
        """
        filter_name = None
        for name in reversed(lineage_step_names):
            if name.startswith("filter_"):
                filter_name = name
                break
        if not filter_name:
            return True

        if filter_name not in _filter_allowed_cache:
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
            _filter_allowed_cache[filter_name] = (allow_cc, allow_ac or allow_sc)

        allow_cross, allow_auto = _filter_allowed_cache[filter_name]
        return allow_auto if is_auto else allow_cross

    # --- Pre-resolve all stack lineage strings and lineage_ids ----------------
    # Build: cc_lineage_prefix -> list of (lineage_to_stack, stack_step, lineage_id)
    # This is O(stack_steps * stack_paths) — tiny, done once.
    _cc_prefix_to_targets: dict = {}  # cc_lineage_str -> list[(lin_to_stack, stack_step, lid, priority)]

    for stack_step in stack_steps:
        stack_paths = get_lineages_to_step_id(session, step_id=stack_step.step_id, include_self=True)
        priority = getattr(stack_step, "priority", 0) or 0

        for path in stack_paths:
            lin_names = [s.step_name for s in path if s.step_name and "global" not in s.step_name]
            if not lin_names:
                continue
            lineage_to_stack = "/".join(lin_names)
            lineage_id = _get_or_create_lineage_id(session, lineage_to_stack)

            # The CC prefix is everything up to (but not including) the stack step.
            # lin_names ends with stack_N; everything before it is the CC lineage.
            cc_prefix_names = lin_names[:-1]
            # Also strip any trailing filter_ pass-through nodes so we can match
            # actual CC job lineage strings (which end at cc_N, not filter_N).
            while cc_prefix_names and cc_prefix_names[-1].startswith("filter_"):
                cc_prefix_names = cc_prefix_names[:-1]
            cc_prefix = "/".join(cc_prefix_names)

            if cc_prefix not in _cc_prefix_to_targets:
                _cc_prefix_to_targets[cc_prefix] = []
            _cc_prefix_to_targets[cc_prefix].append(
                (lineage_to_stack, lin_names, stack_step, lineage_id, priority)
            )

    logger.debug(f"[--after cc] resolved {len(_cc_prefix_to_targets)} distinct CC→STACK lineage mappings")

    # --- Build the full desired set of (step_id, day, pair, lineage_id) -------
    # desired: dict keyed by (step_id, day, pair, lineage_id) ->
    #          (lineage_to_stack, stack_step, priority)
    desired: dict = {}

    for cc_step in cc_steps:
        # Single query: fetch only the columns we need (avoid loading full ORM objects)
        done_cc_rows = (
            session.query(Job.day, Job.pair, Job.lineage_id)
            .filter(Job.step_id == cc_step.step_id)
            .filter(Job.flag == "D")
            .all()
        )
        logger.info(f"[--after cc] cc_step={cc_step.step_name} done_cc_jobs={len(done_cc_rows)}")
        if not done_cc_rows:
            continue

        # Build a map lineage_id -> lineage_str for all CC jobs at once
        cc_lineage_ids = {row.lineage_id for row in done_cc_rows if row.lineage_id is not None}
        from .msnoise_table_def import Lineage as _Lineage
        lid_to_str: dict = {}
        if cc_lineage_ids:
            rows = (
                session.query(_Lineage.lineage_id, _Lineage.lineage_str)
                .filter(_Lineage.lineage_id.in_(cc_lineage_ids))
                .all()
            )
            lid_to_str = {r.lineage_id: r.lineage_str for r in rows}

        for row in done_cc_rows:
            lineage_str = lid_to_str.get(row.lineage_id, "")
            if not lineage_str:
                raise ValueError("DONE CC job has empty lineage (v2 assumption)")

            # Normalise: strip global_ nodes, keep the rest
            parent_names = [n for n in lineage_str_to_step_names(lineage_str, sep="/") if "global" not in n]
            cc_prefix = "/".join(parent_names)

            targets = _cc_prefix_to_targets.get(cc_prefix)
            if not targets:
                continue

            is_auto = _pair_is_auto(row.pair)

            for (lineage_to_stack, lin_names, stack_step, lineage_id, priority) in targets:
                if not _lineage_allowed_for_pair(lin_names, is_auto):
                    continue
                key = (stack_step.step_id, row.day, row.pair, lineage_id)
                if key not in desired:
                    desired[key] = (lineage_to_stack, stack_step, priority, lineage_id)

    logger.info(f"[--after cc] desired STACK jobs: {len(desired)}")
    if not desired:
        session.commit()
        return 0

    # --- Find already-existing STACK jobs for these step_ids ------------------
    stack_step_ids = list({k[0] for k in desired})

    # Fetch all existing jobs for those step_ids in one query
    existing_rows = (
        session.query(Job.ref, Job.step_id, Job.day, Job.pair, Job.lineage_id, Job.flag)
        .filter(Job.step_id.in_(stack_step_ids))
        .all()
    )
    existing_map: dict = {}   # (step_id, day, pair, lineage_id) -> (ref, flag)
    for r in existing_rows:
        existing_map[(r.step_id, r.day, r.pair, r.lineage_id)] = (r.ref, r.flag)

    # --- Split desired into insert / bump / skip ------------------------------
    to_insert = []
    to_bump_refs = []

    for key, (lineage_to_stack, stack_step, priority, lineage_id) in desired.items():
        step_id, day, pair, lid = key
        existing = existing_map.get(key)
        if existing is None:
            to_insert.append({
                "day":        day,
                "pair":       pair,
                "jobtype":    stack_step.step_name,
                "step_id":    step_id,
                "priority":   priority,
                "flag":       "T",
                "lastmod":    now,
                "lineage_id": lineage_id,
            })
        else:
            ref, flag = existing
            if flag != "T":
                to_bump_refs.append(ref)

    logger.info(f"[--after cc] to_insert={len(to_insert)}, to_bump={len(to_bump_refs)}")

    # --- Bulk insert new jobs -------------------------------------------------
    if to_insert:
        session.bulk_insert_mappings(Job, to_insert)
        created = len(to_insert)

    # --- Bulk bump existing non-T jobs back to T ------------------------------
    if to_bump_refs:
        from sqlalchemy import update as sa_update
        # Chunk to avoid hitting SQLite/MySQL parameter limits
        _CHUNK = 900
        for i in range(0, len(to_bump_refs), _CHUNK):
            chunk = to_bump_refs[i:i + _CHUNK]
            session.execute(
                sa_update(Job)
                .where(Job.ref.in_(chunk))
                .values(flag="T", lastmod=now)
                .execution_options(synchronize_session=False)
            )
        bumped = len(to_bump_refs)

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

    now = datetime.datetime.now(datetime.timezone.utc)
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

    Performance: fully batched — no per-job DB round-trips.
    """
    from sqlalchemy import update as sa_update
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.now(datetime.timezone.utc)
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

    # Accumulate all desired (step_id, day, pair, lineage_id) across all refstack steps
    # desired: key -> (jobtype, priority)
    desired: dict = {}

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

        # Fetch all DONE REF jobs for this refstack step (one query)
        done_ref_rows = (
            session.query(Job.pair, Job.lineage_id)
            .filter(Job.step_id == ref_step.step_id)
            .filter(Job.flag == "D")
            .filter(Job.day == "REF")
            .all()
        )
        if not done_ref_rows:
            continue

        # Resolve lineage_id -> lineage_str for refstack jobs (one query)
        ref_lineage_ids = {r.lineage_id for r in done_ref_rows if r.lineage_id is not None}
        from .msnoise_table_def import Lineage as _Lineage
        ref_lid_to_str: dict = {}
        if ref_lineage_ids:
            rows = (
                session.query(_Lineage.lineage_id, _Lineage.lineage_str)
                .filter(_Lineage.lineage_id.in_(ref_lineage_ids))
                .all()
            )
            ref_lid_to_str = {r.lineage_id: r.lineage_str for r in rows}

        for stack_step in stack_predecessors:
            # Fetch all DONE stack (pair, day) tuples in one query — no per-pair loop
            done_stack_rows = (
                session.query(Job.pair, Job.day)
                .filter(Job.step_id == stack_step.step_id)
                .filter(Job.flag == "D")
                .filter(Job.day != "REF")
                .all()
            )
            # Build set of done days per pair
            stack_days_by_pair: dict = {}
            for r in done_stack_rows:
                stack_days_by_pair.setdefault(r.pair, set()).add(r.day)

            for ref_row in done_ref_rows:
                pair = ref_row.pair
                refstack_lineage = ref_lid_to_str.get(ref_row.lineage_id, "")
                if not refstack_lineage:
                    continue

                days = stack_days_by_pair.get(pair, set())
                for day in days:
                    for dvv_step in dvv_successors:
                        dvv_lineage = refstack_lineage + "/" + dvv_step.step_name
                        lid = _get_or_create_lineage_id(session, dvv_lineage)
                        key = (dvv_step.step_id, day, pair, lid)
                        if key not in desired:
                            desired[key] = (dvv_step.step_name, 0)

    logger.info(f"[--after refstack] desired MWCS/stretching/wavelet jobs: {len(desired)}")
    if not desired:
        session.commit()
        return 0

    # One bulk SELECT for all existing jobs matching these step_ids
    dvv_step_ids = list({k[0] for k in desired})
    existing_keys: set = set()
    existing_rows = (
        session.query(Job.step_id, Job.day, Job.pair, Job.lineage_id)
        .filter(Job.step_id.in_(dvv_step_ids))
        .filter(Job.day != "REF")
        .filter(Job.day != "DVV")
        .all()
    )
    for r in existing_rows:
        existing_keys.add((r.step_id, r.day, r.pair, r.lineage_id))

    to_insert = []
    for key, (jobtype, priority) in desired.items():
        if key in existing_keys:
            continue
        step_id, day, pair, lid = key
        to_insert.append({
            "day":        day,
            "pair":       pair,
            "jobtype":    jobtype,
            "step_id":    step_id,
            "priority":   priority,
            "flag":       "T",
            "lastmod":    now,
            "lineage_id": lid,
        })

    if to_insert:
        session.bulk_insert_mappings(Job, to_insert)
        created = len(to_insert)

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

    Performance: fully batched — no per-job DB round-trips.
    """
    DVV_TARGET = {
        "mwcs_dtt":    "mwcs_dtt_dvv",
        "stretching":  "stretching_dvv",
        "wavelet_dtt": "wavelet_dtt_dvv",
    }
    target_category = DVV_TARGET.get(source_category)
    if target_category is None:
        raise ValueError(f"Unknown source_category {source_category!r}")

    from sqlalchemy import update as sa_update
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.now(datetime.timezone.utc)
    created = 0

    parent_steps = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category == source_category)
        .all()
    )

    # desired: dvv_lineage_str -> (dvv_step_id, dvv_step_name, dvv_lineage_id)
    desired: dict = {}

    for parent_step in parent_steps:
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
                dvv_lineage_str = parent_lineage_str.rstrip("/") + "/" + dvv_step.step_name
                if dvv_lineage_str not in desired:
                    lid = _get_or_create_lineage_id(session, dvv_lineage_str)
                    desired[dvv_lineage_str] = (dvv_step.step_id, dvv_step.step_name, lid)

    if not desired:
        session.commit()
        logger.info(f"[--after {source_category}] DVV aggregate jobs created: 0")
        return 0

    # One bulk SELECT for all existing DVV sentinel jobs
    dvv_step_ids = list({v[0] for v in desired.values()})
    existing_rows = (
        session.query(Job.ref, Job.lineage_id, Job.flag)
        .filter(Job.step_id.in_(dvv_step_ids))
        .filter(Job.day == "DVV")
        .all()
    )
    existing_by_lid: dict = {}  # lineage_id -> (ref, flag)
    for r in existing_rows:
        existing_by_lid[r.lineage_id] = (r.ref, r.flag)

    to_insert = []
    to_bump_refs = []

    for dvv_lineage_str, (step_id, step_name, lid) in desired.items():
        existing = existing_by_lid.get(lid)
        if existing is None:
            to_insert.append({
                "day":        "DVV",
                "pair":       "ALL",
                "jobtype":    step_name,
                "step_id":    step_id,
                "priority":   0,
                "flag":       "T",
                "lastmod":    now,
                "lineage_id": lid,
            })
        else:
            ref, flag = existing
            if flag != "T":
                to_bump_refs.append(ref)

    if to_insert:
        session.bulk_insert_mappings(Job, to_insert)
        created = len(to_insert)

    if to_bump_refs:
        _CHUNK = 900
        for i in range(0, len(to_bump_refs), _CHUNK):
            chunk = to_bump_refs[i:i + _CHUNK]
            session.execute(
                sa_update(Job)
                .where(Job.ref.in_(chunk))
                .values(flag="T", lastmod=now)
                .execution_options(synchronize_session=False)
            )

    session.commit()
    logger.info(
        f"[--after {source_category}] DVV aggregate jobs created: {created}, bumped: {len(to_bump_refs)}"
    )
    return created + len(to_bump_refs)


def create_passthrough_jobs_from_done_parent(session, parent_step, child_step):
    """
    DONE parent_step jobs -> TODO child_step jobs, preserving (day, pair).

    This is the "HPC passthrough" rule used from CC downstream in your setup.

    Performance: fully batched — lineage paths and filter gating computed once,
    existing jobs fetched in one query, inserts via bulk_insert_mappings.
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job

    now = datetime.datetime.now(datetime.timezone.utc)
    created = 0

    # Fetch only the columns we need
    done_parent_rows = (
        session.query(Job.day, Job.pair, Job.lineage_id)
        .filter(Job.step_id == parent_step.step_id)
        .filter(Job.flag == "D")
        .filter(Job.day != "REF")
        .filter(Job.day != "DVV")
        .all()
    )
    if not done_parent_rows:
        return 0

    # Resolve all parent lineage_ids to strings in one query
    parent_lids = {r.lineage_id for r in done_parent_rows if r.lineage_id is not None}
    from .msnoise_table_def import Lineage as _Lineage
    lid_to_str: dict = {}
    if parent_lids:
        rows = (
            session.query(_Lineage.lineage_id, _Lineage.lineage_str)
            .filter(_Lineage.lineage_id.in_(parent_lids))
            .all()
        )
        lid_to_str = {r.lineage_id: r.lineage_str for r in rows}

    # --- Hoist invariant work: graph paths and filter config ------------------
    child_paths = get_lineages_to_step_id(session, step_id=child_step.step_id, include_self=True)
    child_priority = getattr(child_step, "priority", 0) or 0

    # Cache filter gating: filter_step_name -> (allow_cross, allow_auto)
    _filter_cache: dict = {}

    def _filter_allowed(lin_names, is_auto):
        filter_name = None
        for name in reversed(lin_names):
            if name.startswith("filter_"):
                filter_name = name
                break
        if not filter_name:
            return True
        if filter_name not in _filter_cache:
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
            _filter_cache[filter_name] = (allow_cc, allow_ac or allow_sc)
        allow_cross, allow_auto = _filter_cache[filter_name]
        return allow_auto if is_auto else allow_cross

    # Pre-build: (parent_lineage_str, is_auto) -> list of (child_lin, child_lid)
    # Constant across all jobs with the same lineage+pair-type — computed once.
    _prefix_to_children: dict = {}

    def _build_children_for_prefix(parent_lineage_str, is_auto):
        parent_names = [n for n in lineage_str_to_step_names(parent_lineage_str, sep="/") if "global" not in n]
        if not parent_names:
            raise ValueError("Parent job has empty lineage (v2 assumption)")
        results = []
        seen = set()
        for path in child_paths:
            lin_names = [s.step_name for s in path if s.step_name and "global" not in s.step_name]
            if len(lin_names) < len(parent_names):
                continue
            if lin_names[:len(parent_names)] != parent_names:
                continue
            if not _filter_allowed(lin_names, is_auto):
                continue
            child_lin = "/".join(lin_names)
            if child_lin not in seen:
                seen.add(child_lin)
                lid = _get_or_create_lineage_id(session, child_lin)
                results.append((child_lin, lid))
        if not results:
            raise ValueError(
                f"No valid lineage for passthrough {parent_step.step_name} -> {child_step.step_name} "
                f"from parent_lineage='{parent_lineage_str}' is_auto={is_auto}"
            )
        return results

    # Build desired set in pure Python
    desired: dict = {}  # (step_id, day, pair, lineage_id) -> child_lid

    for row in done_parent_rows:
        parent_lineage_str = lid_to_str.get(row.lineage_id, "")
        if not parent_lineage_str:
            raise ValueError("DONE parent job has empty lineage (v2 assumption)")

        is_auto = ":" in row.pair and row.pair.split(":")[0] == row.pair.split(":")[1]

        cache_key = (parent_lineage_str, is_auto)
        if cache_key not in _prefix_to_children:
            _prefix_to_children[cache_key] = _build_children_for_prefix(parent_lineage_str, is_auto)

        for (child_lin, child_lid) in _prefix_to_children[cache_key]:
            key = (child_step.step_id, row.day, row.pair, child_lid)
            if key not in desired:
                desired[key] = child_lid

    if not desired:
        return 0

    # One bulk SELECT for existing child jobs
    existing_keys: set = set()
    existing_rows = (
        session.query(Job.step_id, Job.day, Job.pair, Job.lineage_id)
        .filter(Job.step_id == child_step.step_id)
        .all()
    )
    for r in existing_rows:
        existing_keys.add((r.step_id, r.day, r.pair, r.lineage_id))

    to_insert = []
    for key, child_lid in desired.items():
        if key in existing_keys:
            continue
        step_id, day, pair, lid = key
        to_insert.append({
            "day":        day,
            "pair":       pair,
            "jobtype":    child_step.step_name,
            "step_id":    step_id,
            "priority":   child_priority,
            "flag":       "T",
            "lastmod":    now,
            "lineage_id": child_lid,
        })

    if to_insert:
        session.bulk_insert_mappings(Job, to_insert)
        created = len(to_insert)

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

    Performance: fully batched — no per-job DB round-trips.

    :returns: Number of jobs created.
    :rtype: int
    """
    from .msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    now = datetime.datetime.now(datetime.timezone.utc)
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

    # desired: (psd_rms_step_id, day, pair, lineage_id) -> (jobtype, priority)
    desired: dict = {}

    for psd_step in psd_steps:
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

        # Fetch only the columns we need (one query per psd_step)
        done_psd_rows = (
            session.query(Job.day, Job.pair)
            .filter(Job.step_id == psd_step.step_id)
            .filter(Job.flag == "D")
            .all()
        )
        if not done_psd_rows:
            continue

        # Lineage is the psd step name — constant for this psd_step
        lineage = psd_step.step_name
        lineage_id = _get_or_create_lineage_id(session, lineage)

        for psd_rms_step in psd_rms_steps:
            for row in done_psd_rows:
                key = (psd_rms_step.step_id, row.day, row.pair, lineage_id)
                if key not in desired:
                    desired[key] = (psd_rms_step.step_name, 0)

    if not desired:
        session.commit()
        logger.info("[--after psd] PSD_RMS jobs created: 0")
        return 0

    # One bulk SELECT for all existing psd_rms jobs
    psd_rms_step_ids = list({k[0] for k in desired})
    existing_keys: set = set()
    existing_rows = (
        session.query(Job.step_id, Job.day, Job.pair, Job.lineage_id)
        .filter(Job.step_id.in_(psd_rms_step_ids))
        .all()
    )
    for r in existing_rows:
        existing_keys.add((r.step_id, r.day, r.pair, r.lineage_id))

    to_insert = []
    for key, (jobtype, priority) in desired.items():
        if key in existing_keys:
            continue
        step_id, day, pair, lid = key
        to_insert.append({
            "day":        day,
            "pair":       pair,
            "jobtype":    jobtype,
            "step_id":    step_id,
            "priority":   priority,
            "flag":       "T",
            "lastmod":    now,
            "lineage_id": lid,
        })

    if to_insert:
        session.bulk_insert_mappings(Job, to_insert)
        created = len(to_insert)

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
    crap_all_jobs_text = set()  # set: O(1) lookup vs O(N) list — critical at large N
    count = 0
    now = datetime.datetime.now(datetime.timezone.utc)

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
                            crap_all_jobs_text.add(jobtxt)
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
                            crap_all_jobs_text.add(jobtxt)
                            count += 1

    return all_jobs, count

def create_station_pairs(stations):
    """
    Create station pairs for cross-station correlation.

    Station pairs are always sorted alphabetically (sta1 < sta2) to ensure
    consistent pair naming regardless of processing order — this is MSNoise's
    convention throughout the codebase.

    :param stations: List of station identifiers (NET.STA.LOC format)
    :return: List of station pairs
    """
    pairs = []
    sorted_stations = sorted(stations)
    for i in range(len(sorted_stations)):
        for j in range(i + 1, len(sorted_stations)):
            sta1, sta2 = sorted_stations[i], sorted_stations[j]
            pair = f"{sta1}:{sta2}"
            pairs.append(pair)
    return pairs


def main(init=False, nocc=False, after=False):
    logger.info('*** Starting: New Jobs (Workflow-aware) ***')

    db = connect()
    # params = get_params(db)

    if after:
        source_category = str(after).strip().lower()

        params = get_params(db)
        if not params.global_.hpc:
            logger.debug(
                f"new_jobs --after {source_category}: hpc=False, propagation is "
                "normally handled inline by worker scripts via propagate_downstream(). "
                "Running anyway for compatibility (e.g. manual re-runs, tests)."
            )
            # Fall through — run the propagation. In hpc=False mode this is
            # a reconciliation pass; inline workers already handle the hot path.

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
    crap_all_jobs_text = set()  # set: O(1) lookup vs O(N) list — critical at large N
    updated_days = []
    nfs = get_new_files(db)
    now = datetime.datetime.now(datetime.timezone.utc)
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
                            crap_all_jobs_text.add(jobtxt)
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

    # In hpc=True mode: propagate preprocess→cc explicitly.
    # In hpc=False mode: preprocess workers call propagate_downstream() inline.
    if not nocc:
        _params = get_params(db)
        if _params.global_.hpc:
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
                                            lineage=job.get('lineage'),
                                            commit=False)
                db.commit()
                logger.info(f"Created {cc_count} CC jobs")
        else:
            logger.debug('hpc=False: preprocess→cc propagation handled by propagate_downstream()')

    return count
