"""MSNoise workflow topology, job management, lineage resolution and scheduling."""
import datetime
import time

import numpy as np
import pandas as pd

from .db import connect, get_logger
from .config import (get_config, get_config_set_details,
                     get_merged_params_for_lineage, get_params)
from ..msnoise_table_def import (Job, WorkflowStep, DataAvailability,
                                WORKFLOW_CHAINS, WORKFLOW_ORDER, Lineage)

def get_workflow_steps(session):
    """Get all steps in a workflow"""
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowStep) \
        .filter(schema.WorkflowStep.is_active.is_(True)) \
        .order_by(schema.WorkflowStep.step_name).all()



def get_workflow_links(session):
    """Get all links in a workflow"""
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowLink) \
        .filter(schema.WorkflowLink.is_active.is_(True)).all()



def get_workflow_graph(session):
    """Return workflow as nodes and edges for visualization, sorted by workflow order"""
    steps = get_workflow_steps(session)
    links = get_workflow_links(session)

    # Sort steps by workflow order (WORKFLOW_ORDER imported from msnoise_table_def)
    def get_workflow_order_key(step):
        try:
            category_order = WORKFLOW_ORDER.index(step.category)
        except ValueError:
            # If category not in predefined order, put it at the end
            category_order = len(WORKFLOW_ORDER)

        # Sort by category first, then by set_number
        return (category_order, step.set_number or 0)

    sorted_steps = sorted(steps, key=get_workflow_order_key)

    nodes = []
    for step in sorted_steps:
        nodes.append({
            "id": step.step_id,
            "name": step.step_name,
            "category": step.category,
            "set_number": step.set_number,
            "description": step.description
        })

    edges = []
    for link in links:
        edges.append({
            "from": link.from_step_id,
            "to": link.to_step_id,
            "type": link.link_type
        })

    return {"nodes": nodes, "edges": edges}



def create_workflow_step(session, step_name, category, set_number, description=None):
    """Create a new workflow step"""
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()

    step = schema.WorkflowStep(
        step_name=step_name,
        category=category,
        set_number=set_number,
        description=description
    )

    session.add(step)
    session.commit()
    return step



def create_workflow_link(session, from_step_id, to_step_id, link_type="default"):
    """Create a link between two workflow steps"""
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()

    # Check if link already exists
    existing = session.query(schema.WorkflowLink).filter(
        schema.WorkflowLink.from_step_id == from_step_id,
        schema.WorkflowLink.to_step_id == to_step_id
    ).first()

    if existing:
        return existing

    link = schema.WorkflowLink(
        from_step_id=from_step_id,
        to_step_id=to_step_id,
        link_type=link_type
    )

    session.add(link)
    session.commit()
    return link



def get_step_successors(session, step_id):
    """Get all steps that this step feeds into"""
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowStep) \
        .join(schema.WorkflowLink, schema.WorkflowStep.step_id == schema.WorkflowLink.to_step_id) \
        .filter(schema.WorkflowLink.from_step_id == step_id) \
        .filter(schema.WorkflowLink.is_active.is_(True)).all()



def get_first_runnable_steps_per_branch(session, source_step_id, skip_categories=None):
    """
    Return the first runnable step on each outgoing branch from source_step_id,
    skipping steps in skip_categories (default: {"filter"}).

    - "Branch" roots are the immediate successors of source_step_id.
    - If skipped nodes fan out, each fan-out is treated as a separate branch.
    - Results are deduped (set behavior) and returned in stable order.
    """
    if skip_categories is None:
        skip_categories = {"filter"}

    # We dedupe by step_id, but preserve stable ordering at the end
    result_by_step_id = {}

    def immediate_successors(step_id):
        # Reuse existing helper (one-hop)
        return [
            s for s in get_step_successors(session, step_id)
            if getattr(s, "is_active", True)
        ]

    def descend_until_runnable(start_step):
        # Explore forward until the first non-skipped step is found for each path
        stack = [start_step]
        visited = set()

        while stack:
            step = stack.pop()
            if step.step_id in visited:
                continue
            visited.add(step.step_id)

            if step.category not in skip_categories:
                result_by_step_id[step.step_id] = step
                continue  # stop this path at first runnable

            for nxt in immediate_successors(step.step_id):
                stack.append(nxt)

    for succ in immediate_successors(source_step_id):
        descend_until_runnable(succ)

    return [result_by_step_id[k] for k in sorted(result_by_step_id)]



def _get_step_predecessors(session, step_id):
    """Get all steps that feed into this step"""
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()

    return session.query(schema.WorkflowStep) \
        .join(schema.WorkflowLink, schema.WorkflowStep.step_id == schema.WorkflowLink.from_step_id) \
        .filter(schema.WorkflowLink.to_step_id == step_id) \
        .filter(schema.WorkflowLink.is_active.is_(True)).all()



def get_upstream_steps_for_step_id(session, step_id, topo_order=True, include_self=False):
    """
    Returns all upstream WorkflowStep nodes (recursive predecessors) for `step_id`.

    Parameters
    ----------
    session : sqlalchemy.orm.session.Session
    step_id : int
        The step_id of the node whose upstream lineage you want.
    topo_order : bool
        If True (default), returns nodes in dependency order (upstream first),
        suitable for "merge config" application.
        If False, returns in discovery order (still de-duplicated).
    include_self : bool
        If True, includes the node identified by `step_id` in the returned list.

    Returns
    -------
    list[WorkflowStep]
        De-duplicated list of upstream steps.
    """
    visited = set()
    ordered = []

    def dfs(current_step_id):
        # Direct predecessors of the current node:
        preds = _get_step_predecessors(session, current_step_id) or []
        for p in preds:
            if p.step_id in visited:
                continue
            visited.add(p.step_id)
            dfs(p.step_id)
            if topo_order:
                # post-order ensures: preprocess -> cc -> filter -> ...
                ordered.append(p)
            else:
                ordered.append(p)

    dfs(step_id)

    if include_self:
        self_step = session.query(WorkflowStep).filter(WorkflowStep.step_id == step_id).first()
        if self_step is not None and self_step.step_id not in visited:
            ordered.append(self_step)

    return ordered



def create_workflow_steps_from_config_sets(session):
    """
    Create workflow steps automatically from all existing config sets,
    sorted by natural workflow order.

    Returns:
        tuple: (created_count, existing_count, error_message)
    """
    from ..msnoise_table_def import declare_tables

    schema = declare_tables()

    try:
        # Get all unique category+set_number combinations
        config_sets = session.query(
            schema.Config.category,
            schema.Config.set_number
        ).filter(
            schema.Config.set_number.isnot(None)  # Exclude global configs
        ).distinct().all()

        # Sort by workflow order
        def get_workflow_order(config_set):
            category, set_number = config_set
            try:
                return WORKFLOW_ORDER.index(category)
            except ValueError:
                # If category not in predefined order, put it at the end
                return len(WORKFLOW_ORDER)

        config_sets = sorted(config_sets, key=get_workflow_order)
        created_count = 0
        existing_count = 0

        for category, set_number in config_sets:
            # Check if step already exists
            existing_step = session.query(schema.WorkflowStep).filter(
                schema.WorkflowStep.category == category,
                schema.WorkflowStep.set_number == set_number,
            ).first()

            if not existing_step:
                step_name = f"{category}_{set_number}"
                description = f"Auto-generated step for {category} configuration set {set_number}"

                create_workflow_step(
                    session,
                    step_name,
                    category,
                    set_number,
                    description
                )
                created_count += 1
            else:
                existing_count += 1

        return created_count, existing_count, None

    except Exception as e:
        return 0, 0, str(e)



def create_workflow_links_from_steps(session):
    """
    Create workflow links automatically between existing workflow steps,
    following natural workflow progression.

    Returns:
        tuple: (created_count, existing_count, error_message)
    """
    from ..msnoise_table_def import declare_tables

    schema = declare_tables()

    try:
        # Get all workflow steps
        steps = session.query(schema.WorkflowStep).all()

        # Group steps by category and set_number
        steps_by_category = {}
        for step in steps:
            category = step.category
            if category not in steps_by_category:
                steps_by_category[category] = {}
            steps_by_category[category][step.set_number] = step

        created_count = 0
        existing_count = 0

        # Create links based on workflow chains
        # WORKFLOW_CHAINS imported from msnoise_table_def; use simple next_steps list
        for source_category, target_categories in WORKFLOW_CHAINS.items():
            # msnoise_table_def WORKFLOW_CHAINS values are dicts with 'next_steps'
            if isinstance(target_categories, dict):
                target_categories = target_categories.get('next_steps', [])
            if source_category not in steps_by_category:
                continue

            # For each source step in this category
            for source_set_number, source_step in steps_by_category[source_category].items():

                # Link to each target category
                for target_category in target_categories:
                    if target_category not in steps_by_category:
                        continue

                    # Strategy: Create links based on workflow logic
                    target_steps_to_link = []

                    if source_category == 'global':
                        # Global steps link to all steps in target categories
                        target_steps_to_link = list(steps_by_category[target_category].values())

                    elif target_category in [
                        'filter', 'refstack',
                        'mwcs', 'stretching', 'wavelet',
                        'mwcs_dtt', 'wavelet_dtt',
                        'mwcs_dtt_dvv', 'stretching_dvv', 'wavelet_dtt_dvv',
                    ]:
                        # For all processing steps that can have multiple instances,
                        # link to all target steps in the category
                        # This includes DTT steps which can process results from multiple MWCS/wavelet sets
                        target_steps_to_link = list(steps_by_category[target_category].values())

                    else:
                        # Default: link to matching set numbers if available, otherwise to all
                        if source_set_number in steps_by_category[target_category]:
                            target_steps_to_link = [steps_by_category[target_category][source_set_number]]
                        else:
                            target_steps_to_link = list(steps_by_category[target_category].values())

                    # Create the links
                    for target_step in target_steps_to_link:
                        # Check if link already exists
                        existing_link = session.query(schema.WorkflowLink).filter(
                            schema.WorkflowLink.from_step_id == source_step.step_id,
                            schema.WorkflowLink.to_step_id == target_step.step_id
                        ).first()

                        if not existing_link:
                            create_workflow_link(
                                session,
                                source_step.step_id,
                                target_step.step_id,
                                'default'
                            )
                            created_count += 1
                        else:
                            existing_count += 1

        return created_count, existing_count, None

    except Exception as e:
        return 0, 0, str(e)

def _get_or_create_lineage_id(session, lineage_str):
    """Return the lineage_id for *lineage_str*, creating a Lineage row if needed.

    Checks (in order):
    1. Already-loaded persistent objects in the session identity map.
    2. Pending new objects in session.new (not yet flushed).
    3. DB query wrapped in no_autoflush.
    This avoids UNIQUE constraint violations from duplicate inserts.
    """
    if lineage_str is None:
        return None

    # 1. Check identity map (persistent rows already loaded this session)
    for obj in session.identity_map.values():
        if isinstance(obj, Lineage) and obj.lineage_str == lineage_str:
            return obj.lineage_id

    # 2. Check pending new objects (added but not yet flushed)
    for obj in list(session.new):
        if isinstance(obj, Lineage) and obj.lineage_str == lineage_str:
            # Already pending — flush to get its ID, then return it
            session.flush()
            return obj.lineage_id

    # 3. Query the DB (no_autoflush prevents re-entrant flush)
    with session.no_autoflush:
        row = session.query(Lineage).filter(
            Lineage.lineage_str == lineage_str).first()
    if row is not None:
        return row.lineage_id

    # 4. Truly new — insert and flush to get the generated ID
    row = Lineage(lineage_str=lineage_str)
    session.add(row)
    session.flush()
    return row.lineage_id


def _lineage_id_for(session, lineage_str):
    """Return lineage_id for *lineage_str* without creating, or None."""
    if lineage_str is None:
        return None
    for obj in session.identity_map.values():
        if isinstance(obj, Lineage) and obj.lineage_str == lineage_str:
            return obj.lineage_id
    for obj in list(session.new):
        if isinstance(obj, Lineage) and obj.lineage_str == lineage_str:
            return obj.lineage_id  # may be None if not yet flushed
    with session.no_autoflush:
        row = session.query(Lineage).filter(
            Lineage.lineage_str == lineage_str).first()
    return row.lineage_id if row else None


def _lineage_str_from_id(session, lineage_id):
    """Resolve *lineage_id* back to the lineage string.  Returns None if not found.

    Checks identity map first (no SQL), falls back to a no_autoflush query.
    """
    if lineage_id is None:
        return None
    for obj in session.identity_map.values():
        if isinstance(obj, Lineage) and obj.lineage_id == lineage_id:
            return obj.lineage_str
    with session.no_autoflush:
        row = session.query(Lineage).filter(Lineage.lineage_id == lineage_id).first()
    return row.lineage_str if row else None


def propagate_downstream(session, batch: dict) -> int:
    """Propagate a just-completed batch to all immediate downstream worker steps.

    Called immediately after :func:`massive_update_job` marks a batch Done,
    **only when** ``params.global_.hpc`` is ``False`` (the default).  In HPC
    mode the operator runs ``msnoise new_jobs --after <category>`` manually.

    Delegates to the specialised ``propagate_*`` functions from
    ``s02_new_jobs`` for transitions that require non-trivial logic (fan-outs,
    sentinel jobs, etc.).  For simple pair-preserving transitions it performs a
    targeted bulk upsert keyed on the batch's exact (pair, day) tuples.

    Parameters
    ----------
    session : sqlalchemy.orm.Session
    batch : dict
        Return value of :func:`get_next_lineage_batch`.

    Returns
    -------
    int
        Total jobs created or bumped to T.
    """
    completed_step  = batch["step"]
    category        = completed_step.category
    total           = 0

    # ── Special transitions: delegate to existing propagate_* functions ──────
    # These functions are already correct, tested, and handle all the edge cases
    # (sentinel jobs, fan-outs, day="REF" semantics, etc.).

    if category == "preprocess":
        # Fan-out: one station → many station-pairs
        from ..s02_new_jobs import create_cc_jobs_from_preprocess
        cc_jobs, n_cc = create_cc_jobs_from_preprocess(session)
        for _j in cc_jobs:
            update_job(session, _j["day"], _j["pair"],
                       _j["jobtype"], _j["flag"],
                       step_id=_j.get("step_id"),
                       lineage=_j.get("lineage"),
                       commit=False)
        if cc_jobs:
            session.commit()
        return n_cc

    if category == "stack":
        # Two things happen when stack completes:
        #
        # 1. Create a day="REF" sentinel per pair so refstack can check
        #    whether the reference window needs recomputing.
        #
        # 2. ALSO directly create mwcs/stretching/wavelet T jobs for the
        #    new (pair, day) tuples in this batch.  These jobs can run
        #    immediately if the reference is already up-to-date (the
        #    common daily-operations case where new days fall outside the
        #    fixed ref_begin..ref_end window).
        #
        # Refstack will check whether any Done stack days fall inside
        # ref_begin..ref_end.  If none do, it marks the REF job Done
        # without recomputing — the mwcs jobs created here are already
        # waiting.  If the window IS affected, refstack recomputes the
        # reference and propagate_downstream fires again from refstack,
        # bumping existing mwcs T jobs (idempotent) or creating new ones.
        from ..s02_new_jobs import (propagate_refstack_jobs_from_stack_done,
                                     propagate_mwcs_jobs_from_refstack_done)
        total  = propagate_refstack_jobs_from_stack_done(session)
        # Direct propagation for the new days — bypasses the refstack gate
        # when the reference is already up-to-date.
        total += propagate_mwcs_jobs_from_refstack_done(session)
        return total

    if category == "refstack":
        # Looks up MOV stack days and creates per-(pair,day) mwcs/str/wct jobs
        from ..s02_new_jobs import propagate_mwcs_jobs_from_refstack_done
        return propagate_mwcs_jobs_from_refstack_done(session)

    if category in ("mwcs_dtt", "stretching", "wavelet_dtt"):
        # DVV sentinel: day="DVV", pair="ALL"
        from ..s02_new_jobs import propagate_dvv_jobs_from_dtt_done
        return propagate_dvv_jobs_from_dtt_done(session, category)

    if category == "psd":
        # Single-station → single-station, simple pass-through
        from ..s02_new_jobs import propagate_psd_rms_jobs_from_psd_done
        return propagate_psd_rms_jobs_from_psd_done(session)

    # ── Generic pair-preserving transition ────────────────────────────────────
    # For cc→filter(→stack), mwcs→dtt, wavelet→wct_dtt, etc.
    # The completed batch's (pair, day) tuples are used directly.

    from ..msnoise_table_def import declare_tables as _dt
    schema   = _dt()
    Job      = schema.Job
    WFStep   = schema.WorkflowStep
    WFLink   = schema.WorkflowLink

    completed_names = batch["lineage_names"]
    days  = list(dict.fromkeys(batch["days"]))           # unique, order-preserving
    pairs = list(dict.fromkeys(j.pair for j in batch["jobs"]))

    PASSTHROUGH = frozenset({"filter", "global"})

    successors = (
        session.query(WFStep)
        .join(WFLink, WFStep.step_id == WFLink.to_step_id)
        .filter(WFLink.from_step_id == completed_step.step_id)
        .filter(WFLink.is_active.is_(True))
        .all()
    )

    now = datetime.datetime.now(datetime.timezone.utc)

    for succ in successors:
        # Collect worker steps reachable through pass-throughs
        worker_targets: list = []

        def _collect(step, intermediates):
            if step.category in PASSTHROUGH:
                nxt = (
                    session.query(WFStep)
                    .join(WFLink, WFStep.step_id == WFLink.to_step_id)
                    .filter(WFLink.from_step_id == step.step_id)
                    .filter(WFLink.is_active.is_(True))
                    .all()
                )
                for ns in nxt:
                    _collect(ns, intermediates + [step.step_name])
            else:
                worker_targets.append((step, intermediates))

        _collect(succ, [])

        for worker_step, intermediates in worker_targets:
            downstream_names = completed_names + intermediates + [worker_step.step_name]
            downstream_lin   = "/".join(downstream_names)
            downstream_lid   = _get_or_create_lineage_id(session, downstream_lin)

            existing = {
                (j.pair, j.day): j
                for j in session.query(Job)
                .filter(Job.step_id == worker_step.step_id)
                .filter(Job.lineage_id == downstream_lid)
                .filter(Job.pair.in_(pairs))
                .filter(Job.day.in_(days))
                .all()
            }

            to_insert = []
            for pair in pairs:
                for day in days:
                    ej = existing.get((pair, day))
                    if ej is None:
                        to_insert.append({
                            "day":        day,
                            "pair":       pair,
                            "flag":       "T",
                            "jobtype":    worker_step.step_name,
                            "step_id":    worker_step.step_id,
                            "priority":   getattr(worker_step, "priority", 0) or 0,
                            "lastmod":    now,
                            "lineage_id": downstream_lid,
                        })
                    elif ej.flag not in ("T", "I"):
                        ej.flag    = "T"
                        ej.lastmod = now
                        total += 1

            if to_insert:
                session.bulk_insert_mappings(Job, to_insert)
                total += len(to_insert)

    if total:
        session.commit()

    return total


def update_job(session, day, pair, jobtype, flag,
               step_id=None, priority=0, lineage=None,
               commit=True, returnjob=True, ref=None):
    """
    Updates or Inserts a :class:`~msnoise.msnoise_table_def.declare_tables.Job`
    in the database.  Workflow-aware: handles ``step_id`` and ``lineage`` fields.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type day: str
    :param day: The day in YYYY-MM-DD format
    :type pair: str
    :param pair: the name of the station pair
    :type jobtype: str
    :param jobtype: Job type string, e.g. ``"preprocess_1"``
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one, "F"ailed.
    :type step_id: int or None
    :param step_id: WorkflowStep primary key, or None
    :type priority: int
    :param priority: Job priority (default 0)
    :type lineage: str or None
    :param lineage: Lineage string encoding upstream configset chain
    :type commit: bool
    :param commit: Whether to directly commit (True, default) or not (False)
    :type returnjob: bool
    :param returnjob: Return the modified/inserted Job (True, default) or not (False)
    :type ref: int or None
    :param ref: If provided, look up the job by its primary key instead of
        (day, pair, jobtype, lineage).

    :rtype: :class:`~msnoise.msnoise_table_def.declare_tables.Job` or None
    :returns: If returnjob is True, returns the modified/inserted Job.
    """
    from sqlalchemy import text

    if ref:
        job = session.query(Job).filter(text("ref=:ref")).params(ref=ref).first()
    else:
        _lineage_id = _get_or_create_lineage_id(session, lineage)
        job = (session.query(Job)
               .filter(Job.day == day)
               .filter(Job.pair == pair)
               .filter(Job.jobtype == jobtype)
               .filter(
                   Job.lineage_id == _lineage_id if _lineage_id is not None
                   else Job.lineage_id.is_(None)
               )
               .first())

    if job is None:
        job = Job()
        job.day = day
        job.pair = pair
        job.jobtype = jobtype
        job.step_id = step_id
        job.priority = priority
        job.flag = flag
        job.lastmod = datetime.datetime.now(datetime.timezone.utc)
        job.lineage_id = _lineage_id
        session.add(job)
    else:
        # Never demote a Done job back to Todo
        if not (job.flag == "D" and flag == "T"):
            job.flag = flag
        job.step_id = step_id
        job.priority = priority
        job.lastmod = datetime.datetime.now(datetime.timezone.utc)
        job.lineage_id = _lineage_id

    if commit:
        session.commit()
    if returnjob:
        return job



def massive_insert_job(session, jobs):
    """
    Bulk-insert a list of job dicts into the jobs table.

    Each dict must contain at minimum ``day``, ``pair``, ``jobtype``,
    ``flag`` and ``lastmod`` keys.  Optional keys: ``step_id``, ``priority``,
    ``lineage``.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobs: list[dict]
    :param jobs: Job records to insert.
    """
    # Resolve lineage strings to IDs before bulk insert
    # (bulk_insert_mappings bypasses ORM events so we must do this explicitly)
    _lin_cache: dict = {}
    def _resolve(lin_str):
        if lin_str is None:
            return None
        if lin_str not in _lin_cache:
            _lin_cache[lin_str] = _get_or_create_lineage_id(session, lin_str)
        return _lin_cache[lin_str]

    job_records = [
        {
            'day':        j['day'],
            'pair':       j['pair'],
            'jobtype':    j['jobtype'],
            'step_id':    j.get('step_id'),
            'priority':   j.get('priority', 0),
            'flag':       j['flag'],
            'lastmod':    j['lastmod'],
            'lineage_id': _resolve(j.get('lineage')),
        }
        for j in jobs
    ]
    session.bulk_insert_mappings(Job, job_records)
    session.commit()



def massive_update_job(session, jobs, flag="D"):
    """
    Routine to use a low level function to update much faster a list of
    :class:`~msnoise.msnoise_table_def.declare_tables.Job`. This method uses the Job.ref
    which is unique.
    :type session: Session
    :param session: the database connection object
    :type jobs: list or tuple
    :param jobs: a list of :class:`~msnoise.msnoise_table_def.declare_tables.Job` to update.
    :type flag: str
    :param flag: The destination flag.
    """
    updated = False
    mappings = [{'ref': job.ref, 'flag': flag} for job in jobs]
    while not updated:
        try:
            session.bulk_update_mappings(Job, mappings)
            session.commit()
            updated = True
        except Exception:
            time.sleep(np.random.random())
            pass
    return



def reset_jobs(session, jobtype, alljobs=False, reset_i=True, reset_e=True):
    """Reset jobs with the given ``jobtype`` string back to "T"odo.

    Works with the v2 workflow model where ``jobtype`` is a step name such
    as ``"cc_1"`` or ``"refstack_1"``.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobtype: str
    :param jobtype: Step name to reset (e.g. ``"cc_1"``)
    :type alljobs: bool
    :param alljobs: If True reset all jobs regardless of current flag;
        otherwise only resets "I" and/or "E" flagged jobs.
    :type reset_i: bool
    :param reset_i: Reset "I"n-progress jobs (default True)
    :type reset_e: bool
    :param reset_e: Reset "E"rror/failed jobs (default True)
    """
    from sqlalchemy import update as sa_update
    q = sa_update(Job).where(Job.jobtype == jobtype)
    if alljobs:
        session.execute(q.values(flag='T'))
    else:
        flags = []
        if reset_i:
            flags.append('I')
        if reset_e:
            flags.append('F')
        if flags:
            session.execute(q.where(Job.flag.in_(flags)).values(flag='T'))
    session.commit()



def get_job_types(session, jobtype):
    """Return job counts grouped by flag for a given ``jobtype`` string.

    Works with the v2 workflow model where ``jobtype`` is a step name such
    as ``"cc_1"``.  Returns a list of ``(count, flag)`` tuples.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobtype: str
    :param jobtype: Step name to query (e.g. ``"cc_1"``)
    :rtype: list of (int, str)
    :returns: List of (count, flag) pairs
    """
    from sqlalchemy import func
    rows = (session.query(func.count(Job.flag), Job.flag)
            .filter(Job.jobtype == jobtype)
            .group_by(Job.flag)
            .all())
    return rows


# ── Workflow-aware job API ──────────────────────────────────


def get_next_job_for_step(
        session,
        step_category="preprocess",
        flag="T",
        group_by="day",
        limit_days=None,
):
    """
    Return a claimed batch of jobs for a workflow step category.

    group_by:
      - "day": claim all jobs for the selected (step_id, jobtype, day)
      - "pair": claim all jobs for the selected (step_id, jobtype, pair)
      - "pair_lineage": claim all jobs for the selected (step_id, jobtype, pair, lineage)
      - "day_lineage": claim all jobs for the selected (step_id, jobtype, day, lineage)
    """
    from ..msnoise_table_def import declare_tables
    from sqlalchemy import update

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    # refstack jobs always have day="REF"; all other categories never do.
    if step_category == "refstack":
        day_filter = (Job.day == "REF")
    else:
        day_filter = (Job.day != "REF")

    next_job = (
        session.query(Job, WorkflowStep)
        .join(WorkflowStep, Job.jobtype == WorkflowStep.step_name)
        .filter(WorkflowStep.category == step_category)
        .filter(Job.flag == flag)
        .filter(day_filter)
        .order_by(Job.priority.desc(), Job.lastmod)
        .first()
    )
    if not next_job:
        return [], None

    seed_job, step = next_job
    batch_q = (
        session.query(Job)
        .filter(Job.step_id == step.step_id)
        .filter(Job.jobtype == seed_job.jobtype)
        .filter(day_filter)
        .filter(Job.flag == flag)
    )

    if group_by == "day":
        batch_q = batch_q.filter(Job.day == seed_job.day)
        batch_q = batch_q.order_by(Job.pair.asc(), Job.lineage_id.asc())
    elif group_by == "pair":
        batch_q = batch_q.filter(Job.pair == seed_job.pair)
        batch_q = batch_q.order_by(Job.day.asc(), Job.lineage_id.asc())
    elif group_by == "pair_lineage":
        batch_q = batch_q.filter(Job.pair == seed_job.pair).filter(Job.lineage_id == seed_job.lineage_id)
        batch_q = batch_q.order_by(Job.day.asc())
        if limit_days is not None:
            batch_q = batch_q.limit(int(limit_days))
    elif group_by == "day_lineage":
        batch_q = batch_q.filter(Job.day == seed_job.day).filter(Job.lineage_id == seed_job.lineage_id)
        batch_q = batch_q.order_by(Job.pair.asc())
    else:
        raise ValueError(f"Unsupported group_by={group_by!r}")

    jobs = batch_q.with_for_update().all()
    if not jobs:
        return [], step

    upd = (
        update(Job)
        .where(Job.step_id == step.step_id)
        .where(Job.jobtype == seed_job.jobtype)
        .where(day_filter)
        .where(Job.flag == flag)
    )

    if group_by == "day":
        upd = upd.where(Job.day == seed_job.day)
    elif group_by == "pair":
        upd = upd.where(Job.pair == seed_job.pair)
    elif group_by == "pair_lineage":
        upd = upd.where(Job.pair == seed_job.pair).where(Job.lineage_id == seed_job.lineage_id)
        if limit_days is not None:
            # Note: SQL UPDATE with LIMIT is not portable. Keep limit_days for SELECT batching only.
            # The UPDATE still claims the whole (pair, lineage) batch.
            pass
    else:  # "day_lineage"
        upd = upd.where(Job.day == seed_job.day).where(Job.lineage_id == seed_job.lineage_id)

    session.execute(upd.values(flag="I"))
    session.commit()

    return jobs, step



def is_next_job_for_step(session, step_category="preprocess", flag='T'):
    """
    Are there any workflow-aware Jobs in the database with the specified
    flag and step category?

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`
    :type step_category: str
    :param step_category: Workflow step category (e.g., "preprocess", "qc")
    :type flag: str
    :param flag: Status of the Job: "T"odo, "I"n Progress, "D"one.

    :rtype: bool
    :returns: True if at least one Job matches the criteria, False otherwise.
    """
    from ..msnoise_table_def import declare_tables

    schema = declare_tables()
    Job = schema.Job
    WorkflowStep = schema.WorkflowStep

    if step_category == "refstack":
        day_filter = (Job.day == "REF")
    else:
        day_filter = (Job.day != "REF")

    job = session.query(Job) \
        .join(WorkflowStep, Job.jobtype == WorkflowStep.step_name) \
        .filter(WorkflowStep.category == step_category) \
        .filter(Job.flag == flag) \
        .filter(day_filter) \
        .first()

    if job is None:
        return False
    else:
        return True



def get_workflow_job_counts(db):
    """
    Get job counts by status for the workflow

    :param db: Database connection
    :return: Dictionary with job counts by status
    """
    from sqlalchemy import func

    counts = db.query(
        Job.flag,
        func.count(Job.ref).label('count')
    ).group_by(Job.flag).all()

    result = {'T': 0, 'I': 0, 'D': 0}
    for flag, count in counts:
        result[flag] = count

    return result

# ============================================================


def get_filter_steps_for_cc_step(session, cc_step_id):
    """
    Get all filter steps that are children of a specific CC step.

    :param session: Database session
    :param cc_step_id: The step_id of the CC step
    :return: List of filter workflow steps that are successors of the CC step
    """
    from ..api import get_step_successors

    # Get all steps that are successors of the CC step
    successor_steps = get_step_successors(session, cc_step_id)

    # Filter to only include steps with category "filter"
    filter_steps = [step for step in successor_steps if step.category == "filter"]

    return filter_steps



def lineage_str_to_step_names(lineage_str, sep="/"):
    """
    Convert a lineage string like "preprocess_1/cc_1/filter_1" into a list of
    step_name strings, preserving order.
    """
    if lineage_str is None:
        raise ValueError("lineage_str is None")
    lineage_str = str(lineage_str).strip()
    if not lineage_str:
        return []
    return [p.strip() for p in lineage_str.split(sep) if p.strip()]



def lineage_str_to_steps(session, lineage_str, sep="/", strict=True):
    """
    Resolve a lineage string to a list of WorkflowStep ORM objects (ordered).

    Parameters
    ----------
    session : sqlalchemy.orm.session.Session
    lineage_str : str
        e.g. "preprocess_1/cc_1/filter_1"
    sep : str
        Separator used in lineage strings.
    strict : bool
        If True, raises if any step_name can't be resolved.
        If False, silently skips missing steps.

    Returns
    -------
    list[WorkflowStep]
    """
    names = lineage_str_to_step_names(lineage_str, sep=sep)
    if not names:
        return []

    rows = (
        session.query(WorkflowStep)
        .filter(WorkflowStep.step_name.in_(names))
        .all()
    )
    by_name = {s.step_name: s for s in rows}

    steps = []
    missing = []
    for name in names:
        s = by_name.get(name)
        if s is None:
            missing.append(name)
            if not strict:
                continue
        else:
            steps.append(s)

    if missing and strict:
        raise ValueError(
            "Could not resolve lineage steps: %s"
            % ", ".join(missing)
        )

    return steps



def get_lineages_to_step_id(
    session,
    step_id,
    include_self=True,
    max_depth=50,
    max_paths=5000,
):
    """
    Enumerate upstream lineages (all distinct paths) ending at `step_id`.

    Returns a list of paths, each path being a list[WorkflowStep] ordered
    upstream -> downstream. This preserves branch structure (unlike
    get_upstream_steps_for_step_id which de-duplicates nodes).

    Safety:
      - max_depth prevents infinite loops in case of bad/cyclic graphs
      - max_paths prevents combinatorial explosion

    Parameters
    ----------
    session : sqlalchemy.orm.session.Session
    step_id : int
    include_self : bool
        If True, each path ends with the step itself.
    max_depth : int
    max_paths : int

    Returns
    -------
    list[list[WorkflowStep]]
    """
    # Resolve node objects on demand
    def get_step(sid):
        return session.query(WorkflowStep).filter(WorkflowStep.step_id == sid).first()

    # We assume the workflow graph is a DAG. If it isn't, we guard with max_depth
    # and a per-path "seen" set.
    paths = []

    def dfs(current_id, suffix_path, seen_ids, depth):
        if depth > max_depth:
            raise RuntimeError(f"Exceeded max_depth={max_depth} while expanding lineage for step_id={step_id}")

        preds = _get_step_predecessors(session, current_id) or []

        # Leaf: no predecessors => we have a complete path
        if not preds:
            # suffix_path currently holds downstream part (from current -> ... -> target)
            final_path = list(reversed(suffix_path))
            paths.append(final_path)
            if len(paths) > max_paths:
                raise RuntimeError(f"Exceeded max_paths={max_paths} while expanding lineage for step_id={step_id}")
            return

        for p in preds:
            if p.step_id in seen_ids:
                # cycle protection (shouldn't happen in a DAG)
                continue
            dfs(
                p.step_id,
                suffix_path + [p],
                seen_ids | {p.step_id},
                depth + 1,
            )

    start = get_step(step_id)
    if start is None:
        return []

    if include_self:
        dfs(step_id, [start], {step_id}, 0)
    else:
        dfs(step_id, [], set(), 0)

    return paths



def get_done_lineages_for_category(session, category):
    """Return all distinct computed lineages for a given workflow step category.

    Queries ``Job`` rows whose associated ``WorkflowStep.category`` matches
    *category* and whose flag is ``'D'`` (done), then de-duplicates and
    resolves each ``lineage`` string into an ordered list of step-name
    strings (upstream → downstream, including the step itself).

    This is the correct way to enumerate output folders for a step that may
    be reached via multiple upstream paths (e.g. multiple filters, multiple
    MWCS configs), because it reflects what was *actually computed* rather
    than what the DAG topology suggests.

    Example::

        lineages = get_done_lineages_for_category(db, 'stretching')
        # → [
        #     ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'stretching_1'],
        #     ['preprocess_1', 'cc_1', 'filter_2', 'stack_1', 'stretching_1'],
        # ]

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param category: Workflow step category, e.g. ``'stretching'``,
        ``'mwcs_dtt'``, ``'wavelet_dtt'``.
    :type category: str
    :rtype: list[list[str]]
    :returns: Sorted list of unique lineage-name lists.
    """
    from sqlalchemy.orm import aliased
    Lin = aliased(Lineage)
    rows = (
        session.query(Lin.lineage_str)
        .join(Job, Job.lineage_id == Lin.lineage_id)
        .join(Job.workflow_step)
        .filter(WorkflowStep.category == category)
        .filter(Job.flag == "D")
        .distinct()
        .all()
    )
    seen = set()
    result = []
    for (lineage_str,) in rows:
        if not lineage_str or lineage_str in seen:
            continue
        seen.add(lineage_str)
        names = lineage_str_to_step_names(lineage_str)
        # Guard: skip empty or single-entry lineages (global-only, no real steps)
        if len(names) >= 2:
            result.append(names)
    result.sort()
    return result



def resolve_lineage_params(session, lineage_names):
    """Resolve a lineage name-list into a fully merged params object.

    Given a list of step-name strings (as returned by
    :func:`get_done_lineages_for_category`), resolves them to
    :class:`~msnoise.msnoise_table_def.WorkflowStep` ORM objects and merges
    every step's configuration into the global params, exactly as the
    processing steps themselves do via :func:`get_next_lineage_batch`.

    Returns ``(lineage_steps, lineage_names, params)`` — the same tuple as
    :func:`get_merged_params_for_lineage` — so callers can use ``params``
    directly (it will have ``components_to_compute``, ``mov_stack``, etc.).

    Example::

        lineage_names = get_done_lineages_for_category(db, 'mwcs_dtt')[0]
        _, _, params = resolve_lineage_params(db, lineage_names)
        mov_stack = params.stack.mov_stack[0]

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param lineage_names: Ordered list of step-name strings, e.g.
        ``['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'mwcs_1', 'mwcs_dtt_1']``.
    :rtype: tuple(list, list[str], LayeredParams)
    """
    orig_params = get_params(session)
    lineage_str = "/".join(lineage_names)
    steps = lineage_str_to_steps(session, lineage_str, sep="/", strict=False)
    return get_merged_params_for_lineage(session, orig_params, {}, steps)

# ============================================================


def get_next_lineage_batch(
        db,
        step_category,
        group_by="pair_lineage",
        loglevel="INFO",
        day_value=None,
):
    """
    Standard worker prolog for lineage-aware steps.

    - Claims the next batch of jobs for `step_category` using `get_next_job_for_step`.
    - Extracts (pair, lineage_str, refs, days).
    - Loads current step config (from the job row).
    - Resolves lineage_str -> WorkflowStep objects.
    - Builds a LayeredParams for that lineage (one layer per category).

    Returns
    -------
    dict with keys:
      - jobs, step
      - pair, lineage_str, lineage_steps
      - lineage_names          — full list including current step name
      - lineage_names_upstream — full minus current step (replaces manual [:-1])
      - lineage_names_mov      — upstream with any refstack_* entries stripped
                                 (used by mwcs/wct/stretching to find MOV CCFs)
      - refs, days
      - step_params, params    — params is a LayeredParams instance
    or None if no jobs were claimed (caller should continue/sleep).
    """
    # Important: keep logging policy consistent with your scripts
    logger = get_logger(f"msnoise.worker.{step_category}", loglevel)

    # Ensure there is work (fast check)
    if not is_next_job_for_step(db, step_category=step_category, flag="T"):
        return None

    jobs, step = get_next_job_for_step(
        db,
        step_category=step_category,
        group_by=group_by,
        flag="T",
        # day_value=day_value,  # enable once you add day filtering to scheduler
    )

    if not jobs:
        return None

    pair = jobs[0].pair
    # Use _lineage_str_from_id to safely resolve after session.commit() expiry
    lineage_str = _lineage_str_from_id(db, jobs[0].lineage_id)
    if not lineage_str:
        raise ValueError(
            f"{step_category.upper()} jobs must have a non-empty lineage "
            f"(lineage_id={jobs[0].lineage_id!r})"
        )

    refs = [job.ref for job in jobs]
    days = [job.day for job in jobs]

    # Load current step config from job row (job carries config_category & config_set_number)
    step_params = get_config_set_details(
        db,
        jobs[0].config_category,
        jobs[0].config_set_number,
        format="AttribDict",
    )

    # Merge params for THIS lineage only
    orig_params = get_params(db)
    lineage_steps = lineage_str_to_steps(db, lineage_str, strict=True)
    lineage_steps, lineage_names, params = get_merged_params_for_lineage(db, orig_params, step_params, lineage_steps)

    # Derived lineage lists — no manual slicing needed in callers
    lineage_names_upstream = lineage_names[:-1] if lineage_names else []
    lineage_names_mov = _strip_refstack_from_lineage(lineage_names_upstream)

    logger.info(f"New {step_category.upper()} batch: pair={pair} n={len(jobs)} group_by={group_by} lineage={lineage_str}")

    return {
        "jobs": jobs,
        "step": step,
        "pair": pair,
        "lineage_str": lineage_str,
        "lineage_steps": lineage_steps,
        "lineage_names": lineage_names,
        "lineage_names_upstream": lineage_names_upstream,
        "lineage_names_mov": lineage_names_mov,
        "refs": refs,
        "days": days,
        "step_params": step_params,
        "params": params,
    }



def get_refstack_lineage_for_filter(session, filterid, refstack_set_number=1):
    """Get the full lineage path through a filter step down to its refstack.

    Extends :func:`get_stack_lineage_for_filter` by also appending the
    ``refstack_M`` step that is immediately downstream of the stack step.
    This is the correct lineage to pass to :func:`xr_get_ref`, since REF
    files now live under the refstack step folder.

    Example for the default single-pipeline::

        get_refstack_lineage_for_filter(db, 1)
        # → ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']

    :type session: :class:`sqlalchemy.orm.session.Session`
    :type filterid: int
    :param filterid: The filter set_number (e.g. 1 for filter_1).
    :type refstack_set_number: int
    :param refstack_set_number: Which refstack set to use (default 1).
    :rtype: list of str
    """
    path = _get_stack_lineage_for_filter(session, filterid)
    if not path:
        return path

    steps = get_workflow_steps(session)
    links = get_workflow_links(session)

    # Find the stack step at the end of the path
    stack_step_name = path[-1]
    stack_step = next(
        (s for s in steps if s.step_name == stack_step_name), None
    )
    if stack_step is None:
        return path

    # Find the refstack step downstream of this stack step
    # Prefer the requested refstack_set_number; fall back to the first available
    refstack_steps = [
        s for s in steps
        if s.category == 'refstack'
        and any(lk.from_step_id == stack_step.step_id and lk.to_step_id == s.step_id
                for lk in links)
    ]
    if not refstack_steps:
        return path  # no refstack in this workflow, return stack-level lineage

    # Pick requested set_number if available, else first
    refstack_step = next(
        (s for s in refstack_steps if s.set_number == refstack_set_number),
        refstack_steps[0]
    )
    return path + [refstack_step.step_name]

# ============================================================


def extend_days(days):
    """Return a :class:`~pandas.DatetimeIndex` from *days* extended by one
    extra day at the end.

    Replaces the pandas 1.x pattern::

        idx = pd.to_datetime(days)
        idx = idx.append(pd.DatetimeIndex([idx[-1] + pd.Timedelta("1d")]))

    which was removed in pandas 2.0.

    :param days: sequence of date-like values (strings, dates, datetimes…)
    :rtype: :class:`pandas.DatetimeIndex`
    """
    idx = pd.to_datetime(days)
    return pd.DatetimeIndex(list(idx) + [idx[-1] + pd.Timedelta("1d")])



def get_t_axis(params):
    """
    Returns the time axis (in seconds) of the CC functions.

    :rtype: :class:`numpy.array`
    :returns: the time axis in seconds
    """
    samples = int(2 * params.cc.maxlag * params.cc.cc_sampling_rate) + 1
    return np.linspace(-params.cc.maxlag, params.cc.maxlag, samples)



def build_ref_datelist(params, session=None):
    """
    Creates a date array for the REF.
    The returned tuple contains a start and an end date, and a list of
    individual dates between the two.

    :type params: :class:`~msnoise.api.AttribDict`
    :param params: A params object as obtained by :func:`get_params`.
    :type session: :class:`sqlalchemy.orm.session.Session`, optional
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`. Required only when ``ref_begin`` is
        ``"1970-01-01"`` (auto-detect from data availability).

    :rtype: tuple
    :returns: (start, end, datelist)
    """
    begin = params.refstack.ref_begin
    end = params.refstack.ref_end
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    elif begin == "1970-01-01":
        if session is None:
            session = connect()
        start = session.query(DataAvailability).order_by(
            DataAvailability.starttime).first()
        if start:
            start = start.starttime.date()
        else:
            start = datetime.date(1970, 1, 1)
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d').date()
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()
    end = min(end, datetime.date.today())
    datelist = pd.date_range(start, end).map(lambda x: x.date())
    return start, end, datelist.tolist()



def build_movstack_datelist(session):
    """
    Creates a date array for the analyse period.
    The returned tuple contains a start and an end date, and a list of
    individual dates between the two.

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object, as
        obtained by :func:`connect`

    :rtype: tuple
    :returns: (start, end, datelist)
    """
    begin = get_config(session, "startdate", category='global', set_number=1)
    end = get_config(session, "enddate", category='global', set_number=1)
    if begin[0] == '-':
        start = datetime.date.today() + datetime.timedelta(days=int(begin))
        end = datetime.date.today() + datetime.timedelta(days=int(end))
    elif begin == "1970-01-01": # TODO this fails when the DA is empty
        try:
            start = session.query(DataAvailability).order_by(
                DataAvailability.starttime).first().starttime.date()
        except Exception:
            start = datetime.date(1970, 1, 1)
    else:
        start = datetime.datetime.strptime(begin, '%Y-%m-%d').date()

    if end == "2100-01-01":
        try:
            end = session.query(DataAvailability).order_by(
                DataAvailability.endtime.desc()).first().endtime.date()
        except Exception:
            end = datetime.date.today()
    else:
        end = datetime.datetime.strptime(end, '%Y-%m-%d').date()
    end = min(end, datetime.date.today())
    datelist = pd.date_range(start, end).map(lambda x: x.date())
    return start, end, datelist.tolist()



def refstack_is_rolling(params):
    """
    Return True if the refstack configset uses rolling-index mode.

    Rolling mode is indicated by ``ref_begin`` being a negative integer string
    (e.g. ``"-5"``). In this mode no REF file is written to disk; the reference
    is computed on-the-fly at MWCS/stretching/WCT time via
    :func:`compute_rolling_ref`.

    :type params: :class:`obspy.core.AttribDict`
    :param params: Merged parameter set containing ``ref_begin``.
    :rtype: bool
    """
    val = str(params.refstack.ref_begin).strip()
    return val.startswith("-")



def compute_rolling_ref(data, ref_begin, ref_end):
    """Compute a per-index rolling reference from CCF data.

    For each time index ``i``, the reference is::

        mean(data[i + ref_begin : i + ref_end])

    Both ``ref_begin`` and ``ref_end`` must be negative integers with
    ``ref_begin < ref_end <= -1``.  Use ``ref_end=-1`` to exclude the current
    window (compare to the previous N windows).

    Uses ``min_periods=1`` semantics: the first few steps receive whatever
    partial window is available rather than NaN.

    :type data: :class:`pandas.DataFrame` or :class:`xarray.DataArray`
        Shape ``(n_times, n_lag_samples)``.
    :type ref_begin: int
        Negative offset for the start of the rolling window (e.g. ``-5``).
    :type ref_end: int
        Negative offset for the end of the rolling window (e.g. ``-1``).
        Must satisfy ``ref_begin < ref_end <= 0``.
    :rtype: :class:`numpy.ndarray`
        Shape ``(n_times, n_lag_samples)``.  Row ``i`` is the reference for
        time step ``i``.
    """
    import xarray as xr_mod
    if isinstance(data, xr_mod.DataArray):
        arr = data.values
    else:
        arr = data.values  # DataFrame.values
    n = arr.shape[0]
    refs = np.full_like(arr, np.nan, dtype=float)

    for i in range(n):
        start_idx = i + ref_begin   # e.g. i - 5
        end_idx   = i + ref_end     # e.g. i - 1  (exclusive in slice)
        if end_idx <= 0:
            continue                # not enough history yet
        start_idx = max(0, start_idx)
        if start_idx >= end_idx:
            continue
        refs[i] = np.nanmean(arr[start_idx:end_idx], axis=0)

    return refs



def _strip_refstack_from_lineage(lineage_names):
    """
    Return a copy of ``lineage_names`` with any ``refstack_*`` entries removed.

    Used internally by :func:`get_next_lineage_batch` to build
    ``lineage_names_mov`` — the path for :func:`xr_get_ccf`, which lives
    under the ``stack_N`` step folder rather than under ``refstack_M``.

    Example::

        _strip_refstack_from_lineage(
            ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']
        )
        # → ['preprocess_1', 'cc_1', 'filter_1', 'stack_1']

    :type lineage_names: list of str
    :rtype: list of str
    """
    return [n for n in lineage_names if not n.startswith("refstack_")]



def filter_within_daterange(date, start_date, end_date):
    """Check if a date falls within the configured range"""
    return start_date <= date <= end_date


# ── CCF ─────────────────────────────────────────────────────


def _get_stack_lineage_for_filter(session, filterid):
    """Get the full lineage path through a specific filter step to its downstream stack.

    Traverses the workflow graph upstream from ``filter_{filterid}`` to the root
    preprocess step, then appends the stack step that is immediately downstream of
    the filter.  Returns a list of step names suitable for use as the *lineage*
    argument to :func:`xr_get_ref`, :func:`xr_get_ccf`, etc.

    Example for the default single-pipeline::

        get_stack_lineage_for_filter(db, 1)
        # → ['preprocess_1', 'cc_1', 'filter_1', 'stack_1']

    :type session: :class:`sqlalchemy.orm.session.Session`
    :type filterid: int
    :param filterid: The filter set_number (e.g. 1 for filter_1).
    :rtype: list of str
    """
    steps = get_workflow_steps(session)
    links = get_workflow_links(session)
    step_map = {s.step_id: s for s in steps}

    filter_step = next(
        (s for s in steps if s.category == 'filter' and s.set_number == filterid),
        None,
    )
    if filter_step is None:
        return []

    # Walk upstream from filter to the root step, skipping global config steps
    parent_map = {link.to_step_id: link.from_step_id for link in links}
    path = [filter_step.step_name]
    current_id = filter_step.step_id
    visited = set()
    while current_id in parent_map:
        if current_id in visited:
            break  # guard against cycles
        visited.add(current_id)
        parent_id = parent_map[current_id]
        parent_step = step_map[parent_id]
        if parent_step.category != 'global':
            path.insert(0, parent_step.step_name)
        current_id = parent_id

    # Append the stack step immediately downstream of this filter (if any)
    for link in links:
        if link.from_step_id == filter_step.step_id:
            child = step_map.get(link.to_step_id)
            if child and child.category == 'stack':
                path.append(child.step_name)
                break

    return path
