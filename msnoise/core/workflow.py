"""MSNoise workflow topology, job management, lineage resolution and scheduling."""

__all__ = [
    "build_movstack_datelist",
    "build_ref_datelist",
    "compute_rolling_ref",
    "create_workflow_link",
    "create_workflow_links_from_steps",
    "create_workflow_step",
    "create_workflow_steps_from_config_sets",
    "extend_days",
    "get_done_lineages_for_category",
    "get_filter_steps_for_cc_step",
    "get_job_types",
    "get_lineages_to_step_id",
    "get_next_job_for_step",
    "get_next_lineage_batch",
    "get_refstack_lineage_for_filter",
    "get_t_axis",
    "get_workflow_job_counts",
    "get_workflow_links",
    "get_workflow_steps",
    "is_next_job_for_step",
    "lineage_str_to_step_names",
    "lineage_str_to_steps",
    "massive_insert_job",
    "massive_update_job",
    "propagate_downstream",
    "refstack_is_rolling",
    "_strip_stack_from_lineage",
    "refstack_needs_recompute",
    "reset_jobs",
    "resolve_lineage_params",
    "update_job",
    "filter_within_daterange",
    "get_first_runnable_steps_per_branch",
    "get_step_successors",
    "get_workflow_graph",
    "get_workflow_chains",
    "get_workflow_order",
    "get_step_abbrevs",
    "get_category_display_info",
    "is_terminal_category",
    "is_entry_category",
    # Job-propagation helpers (used by s02_new_jobs propagators)
    "query_active_steps",
    "resolve_lineage_ids",
    "bulk_upsert_jobs",
    # Workflow runner
    "get_category_cli_command",
    "run_workflow_plan",
]

import collections
import datetime
import time

import numpy as np
import pandas as pd

from .db import connect, get_logger
from .config import (get_config, get_config_set_details,
                     get_merged_params_for_lineage, get_params)
from ..msnoise_table_def import Job, WorkflowStep, DataAvailability, Lineage

# Module-level schema — declare_tables() is idempotent; this gives us WorkflowLink
# and Config without scattering local declare_tables() calls everywhere.
# §10: class identity is safe here because msnoise_table_def itself calls
# declare_tables() at module level, so all calls return the same cached classes.
from ..msnoise_table_def import declare_tables as _declare_tables
_schema = _declare_tables()
WorkflowLink = _schema.WorkflowLink
Config       = _schema.Config

# ── Built-in workflow topology ────────────────────────────────────────────────
# Defines the directed adjacency of workflow categories. Plugin packages can
# extend this via the ``msnoise.plugins.workflow_chains`` entry point group.
# Each callable registered there must return a dict with the same schema:
#   { category: { "next_steps": [...], "is_entry_point": bool, "is_terminal": bool } }

_BUILTIN_WORKFLOW_CHAINS = {
    # Keys: next_steps, is_entry_point, is_terminal
    #       abbrev       → short tag used in plot filenames (lineage_to_plot_tag)
    #       display_name → human label used in admin UI config-sets page
    #       level        → depth in the DAG tree (global=0, its children=1, …)
    #                      psd and preprocess are both level=1 (both branch off global)
    #       cli          → msnoise sub-command token list to run this category's worker
    #                      None for pass-through categories (filter, global) — no jobs
    #
    # Plugin packages can add entries via msnoise.plugins.workflow_chains:
    #   each callable returns a dict with the same schema.
    "global": {
        "next_steps": ["preprocess", "psd"],
        "is_entry_point": True, "is_terminal": False,
        "abbrev": "g", "display_name": "Global Parameters", "level": 0,
        "cli": None,
    },
    "preprocess": {
        "next_steps": ["cc"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "pre", "display_name": "Preprocessing", "level": 1,
        "cli": ["cc", "preprocess"],
    },
    "cc": {
        "next_steps": ["filter"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "cc", "display_name": "Cross-Correlation", "level": 2,
        "cli": ["cc", "compute_cc"],
    },
    "filter": {
        "next_steps": ["stack", "refstack"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "f", "display_name": "Filters", "level": 3,
        "cli": None,
    },
    "stack": {
        # Stack feeds refstack by default. Direct stack→mwcs links (no refstack)
        # are supported but must be added manually via the admin UI — they are NOT
        # auto-created by db upgrade to avoid spurious direct jobs when a refstack
        # step is present in the same chain.
        "next_steps": ["refstack"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "stk", "display_name": "Moving Stacks", "level": 4,
        "cli": ["cc", "stack", "-m"],
    },
    "refstack": {
        # Sibling of stack AND optional intermediate between stack and mwcs.
        # Lineage convention A1: mwcs jobs use  …/stack_N/refstack_M/mwcs_1
        # so that lineage_names_mov (strip refstack) resolves the CCF path.
        # If stack links directly to mwcs (no refstack), lineage is …/stack_N/mwcs_1.
        "next_steps": ["mwcs", "stretching", "wavelet"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "ref", "display_name": "Reference Stacks", "level": 4,
        "cli": ["cc", "stack_refstack"],
    },
    "mwcs": {
        "next_steps": ["mwcs_dtt"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "mwcs", "display_name": "MWCS", "level": 5,
        "cli": ["cc", "dtt", "compute_mwcs"],
    },
    "mwcs_dtt": {
        "next_steps": ["mwcs_dtt_dvv"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "dtt", "display_name": "MWCS dt/t", "level": 6,
        "cli": ["cc", "dtt", "compute_mwcs_dtt"],
    },
    "mwcs_dtt_dvv": {
        "next_steps": [],
        "is_entry_point": False, "is_terminal": True,
        "abbrev": "dvv", "display_name": "MWCS dv/v Aggregate", "level": 7,
        "cli": ["cc", "dtt", "dvv", "compute_mwcs_dtt_dvv"],
    },
    "stretching": {
        "next_steps": ["stretching_dvv"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "str", "display_name": "Stretching", "level": 5,
        "cli": ["cc", "dtt", "compute_stretching"],
    },
    "stretching_dvv": {
        "next_steps": [],
        "is_entry_point": False, "is_terminal": True,
        "abbrev": "sdvv", "display_name": "Stretching dv/v Aggregate", "level": 6,
        "cli": ["cc", "dtt", "dvv", "compute_stretching_dvv"],
    },
    "wavelet": {
        "next_steps": ["wavelet_dtt"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "wct", "display_name": "Wavelet", "level": 5,
        "cli": ["cc", "dtt", "compute_wct"],
    },
    "wavelet_dtt": {
        "next_steps": ["wavelet_dtt_dvv"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "wdtt", "display_name": "Wavelet dt/t", "level": 6,
        "cli": ["cc", "dtt", "compute_wct_dtt"],
    },
    "wavelet_dtt_dvv": {
        "next_steps": [],
        "is_entry_point": False, "is_terminal": True,
        "abbrev": "wdvv", "display_name": "WCT dv/v Aggregate", "level": 7,
        "cli": ["cc", "dtt", "dvv", "compute_wavelet_dtt_dvv"],
    },
    # ── PSD branch (also level=1: child of global, sibling of preprocess) ──
    "psd": {
        "next_steps": ["psd_rms"],
        "is_entry_point": False, "is_terminal": False,
        "abbrev": "psd", "display_name": "PSD", "level": 1,
        "cli": ["qc", "compute_psd"],
    },
    "psd_rms": {
        "next_steps": [],
        "is_entry_point": False, "is_terminal": True,
        "abbrev": "rms", "display_name": "PSD RMS", "level": 2,
        "cli": ["qc", "compute_psd_rms"],
    },
}

# Canonical display order for UI sorting and config-set creation.
# Plugin packages can append extra categories via the
# ``msnoise.plugins.workflow_order`` entry point group.  Each callable must
# return a list of category name strings in the desired insertion order.
_BUILTIN_WORKFLOW_ORDER = [
    "global", "preprocess", "cc", "psd", "psd_rms",
    "filter", "stack", "refstack",
    "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
    "stretching", "stretching_dvv",
    "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
]

def get_workflow_chains():
    """Return the full workflow adjacency map, including any plugin extensions.

    Merges the built-in topology with addenda registered via the
    ``msnoise.plugins.workflow_chains`` entry point group.  Each registered
    callable must return a dict of the same schema as the built-in chains::

        {
            "my_step": {
                "next_steps": ["other_step"],
                "is_entry_point": False,
                "is_terminal": False,
            }
        }

    Plugin entries are merged in load order; later plugins can override
    earlier ones for the same category key.

    :rtype: dict[str, dict]
    """
    from importlib.metadata import entry_points
    chains = dict(_BUILTIN_WORKFLOW_CHAINS)
    for ep in entry_points(group="msnoise.plugins.workflow_chains"):
        try:
            addendum = ep.load()()
            chains.update(addendum)
        except Exception as exc:
            import logging
            logging.getLogger("msnoise.workflow").warning(
                f"Failed to load workflow_chains plugin {ep.name!r}: {exc}"
            )
    return chains


def get_workflow_order():
    """Return the canonical category display order, including plugin extensions.

    Merges the built-in order with any extra categories registered via the
    ``msnoise.plugins.workflow_order`` entry point group.  Each registered
    callable must return a list of category name strings in the desired
    insertion order (appended after the built-in list).

    :rtype: list[str]
    """
    from importlib.metadata import entry_points
    order = list(_BUILTIN_WORKFLOW_ORDER)
    for ep in entry_points(group="msnoise.plugins.workflow_order"):
        try:
            extra = ep.load()()
            for cat in extra:
                if cat not in order:
                    order.append(cat)
        except Exception as exc:
            import logging
            logging.getLogger("msnoise.workflow").warning(
                f"Failed to load workflow_order plugin {ep.name!r}: {exc}"
            )
    return order


def get_category_cli_command(category: str) -> list | None:
    """Return the ``msnoise`` sub-command token list for *category*, or None.

    Reads the ``"cli"`` key from :func:`get_workflow_chains` so that the
    mapping stays in the single source of truth.  ``None`` means the category
    is a pass-through (no worker, no jobs).

    Plugin packages extend the mapping by adding a ``"cli"`` key to the dict
    returned by their ``msnoise.plugins.workflow_chains`` entry point — no
    separate entry point is needed.

    :param category: Workflow category string (e.g. ``"cc"``, ``"stack"``).
    :rtype: list[str] or None
    """
    chains = get_workflow_chains()
    entry = chains.get(category, {})
    if isinstance(entry, dict):
        return entry.get("cli", None)
    return None


# Named tuple for a single step in the execution plan.
RunStep = collections.namedtuple(
    "RunStep",
    ["category", "step_names", "cmd_tokens", "job_count", "hpc_after"],
)
"""
A single entry in the workflow execution plan.

:param category:   Workflow category string (e.g. ``"cc"``).
:param step_names: List of active step names in this category (e.g. ``["cc_1"]``).
:param cmd_tokens: ``msnoise`` sub-command tokens (e.g. ``["cc", "compute"]``),
                   or ``None`` if pass-through.
:param job_count:  Number of ``T`` jobs waiting at plan-generation time.
:param hpc_after:  ``True`` if a ``msnoise new_jobs --after <category>`` call
                   should be inserted after this step (HPC mode).
"""


def run_workflow_plan(
    session,
    *,
    threads: int = 1,
    hpc: bool | None = None,
    from_category: str | None = None,
    until_category: str | None = None,
    chunk_size: int = 0,
) -> list:
    """Build an ordered execution plan for all active workflow steps.

    Performs a topological sort of active ``WorkflowStep`` categories via
    ``WorkflowLink`` rows and returns an ordered list of :class:`RunStep`
    namedtuples.  Pass-through categories (``filter``, ``global``) are
    omitted.  Categories with no runnable command (see
    :func:`get_category_cli_command`) are also omitted.

    The returned plan is pure data — no subprocesses are started.  The CLI
    command ``msnoise utils run_workflow`` drives the actual execution.

    :param session: SQLAlchemy session.
    :param threads: Number of parallel worker threads/processes (``-t N``).
    :param hpc: Override ``hpc`` flag from config.  If ``None``, reads from
        ``params.global_.hpc``.
    :param from_category: Skip all categories before this one (inclusive start).
    :param until_category: Stop after this category (inclusive end).
    :param chunk_size: ``--chunk-size`` value passed to ``cc compute`` and
        ``qc compute_psd``.
    :returns: Ordered list of :class:`RunStep` namedtuples.
    :rtype: list[RunStep]
    """
    from sqlalchemy import func as _func

    # -- Resolve HPC flag -----------------------------------------------------
    if hpc is None:
        params = get_params(session)
        hpc = bool(getattr(getattr(params, "global_", None), "hpc", False))

    # -- Load active steps and links ------------------------------------------
    all_steps = get_workflow_steps(session)
    all_links = get_workflow_links(session)

    if not all_steps:
        return []

    # -- Topological sort at CATEGORY level -----------------------------------
    # Build category adjacency from active links.
    step_to_cat = {s.step_id: s.category for s in all_steps}
    cat_deps: dict = {}   # category -> set of predecessor categories
    active_cats = {s.category for s in all_steps}

    for cat in active_cats:
        cat_deps.setdefault(cat, set())

    for lk in all_links:
        src_cat = step_to_cat.get(lk.from_step_id)
        tgt_cat = step_to_cat.get(lk.to_step_id)
        if src_cat and tgt_cat and src_cat != tgt_cat:
            cat_deps.setdefault(tgt_cat, set()).add(src_cat)

    # Kahn's algorithm — tie-break by _BUILTIN_WORKFLOW_ORDER so the CC chain
    # runs before the independent PSD branch when both are ready simultaneously.
    _wf_order = get_workflow_order()
    def _order_key(cat):
        try:
            return _wf_order.index(cat)
        except ValueError:
            return len(_wf_order)

    in_degree = {cat: len(deps) for cat, deps in cat_deps.items()}
    queue = sorted([c for c, d in in_degree.items() if d == 0], key=_order_key)
    topo_order = []
    while queue:
        cat = queue.pop(0)
        topo_order.append(cat)
        for other, deps in cat_deps.items():
            if cat in deps:
                in_degree[other] -= 1
                if in_degree[other] == 0:
                    queue.append(other)
                    queue.sort(key=_order_key)

    # -- Filter: skip pass-throughs and no-command categories -----------------
    PASSTHROUGH = frozenset({"filter", "global"})
    runnable_order = [
        cat for cat in topo_order
        if cat not in PASSTHROUGH
        and get_category_cli_command(cat) is not None
    ]

    # -- Apply --from / --until ------------------------------------------------
    if from_category and from_category in runnable_order:
        idx = runnable_order.index(from_category)
        runnable_order = runnable_order[idx:]
    if until_category and until_category in runnable_order:
        idx = runnable_order.index(until_category)
        runnable_order = runnable_order[:idx + 1]

    # -- Build per-category step lists and job counts -------------------------
    schema = _declare_tables()
    Job = schema.Job

    steps_by_cat: dict = {}
    for s in all_steps:
        steps_by_cat.setdefault(s.category, []).append(s.step_name)

    # Categories that feed into the next step in HPC mode (need new_jobs --after)
    # = all non-terminal runnable categories
    chains = get_workflow_chains()

    plan: list = []
    for cat in runnable_order:
        cmd_tokens = list(get_category_cli_command(cat))

        # Inject --chunk-size for cc and psd (day_lineage steps)
        if chunk_size > 0 and cat in ("cc", "psd"):
            cmd_tokens += ["--chunk-size", str(chunk_size)]

        # Inject --threads (msnoise -t N is a top-level flag, handled by caller)

        step_names = sorted(steps_by_cat.get(cat, []))

        # Count T jobs for this category
        step_ids = [
            s.step_id for s in all_steps if s.category == cat
        ]
        if step_ids:
            row = (
                session.query(_func.count(Job.ref))
                .filter(Job.step_id.in_(step_ids))
                .filter(Job.flag == "T")
                .scalar()
            )
            job_count = row or 0
        else:
            job_count = 0

        # HPC: insert new_jobs --after unless this is the last runnable step.
        # Walk transitively through PASSTHROUGH categories (e.g. filter → stack)
        # so that cc correctly gets hpc_after=True even though its immediate
        # successor is the passthrough "filter" step.
        def _has_runnable_descendant(cat, _seen=None):
            if _seen is None:
                _seen = set()
            if cat in _seen:
                return False
            _seen.add(cat)
            info = chains.get(cat, {})
            for nc in (info.get("next_steps", []) if isinstance(info, dict) else []):
                if nc in PASSTHROUGH:
                    if _has_runnable_descendant(nc, _seen):
                        return True
                elif get_category_cli_command(nc) is not None:
                    return True
            return False

        hpc_after = hpc and _has_runnable_descendant(cat)

        plan.append(RunStep(
            category=cat,
            step_names=step_names,
            cmd_tokens=cmd_tokens,
            job_count=job_count,
            hpc_after=hpc_after,
        ))

    return plan




def get_step_abbrevs():
    """Return a mapping of workflow category → short plot-filename abbreviation.

    Merges the built-in abbreviations (from :data:`_BUILTIN_WORKFLOW_CHAINS`)
    with any ``abbrev`` keys provided by plugin entry points via
    ``msnoise.plugins.workflow_chains``.  Falls back to the raw category name
    for any category without an ``abbrev`` entry.

    :rtype: dict[str, str]
    """
    chains = get_workflow_chains()
    return {
        cat: entry.get("abbrev", cat) if isinstance(entry, dict) else cat
        for cat, entry in chains.items()
    }


def get_category_display_info():
    """Return ordered display metadata for all workflow categories.

    Each entry is a dict with keys ``category``, ``display_name``, ``level``
    (DAG depth: global=0, preprocess/psd=1, cc=2, …).  Ordered for UI
    rendering: the main CC branch first (preprocess→cc→…→dvv), then the
    PSD branch (psd→psd_rms).  Plugin categories are appended at the end.

    :rtype: list[dict]
    """
    chains = get_workflow_chains()

    # Canonical UI display order: CC branch first, PSD branch at end.
    # This differs from WORKFLOW_ORDER (which is pipeline/topological order)
    # but matches the tree structure users see in the admin config-sets page.
    _UI_ORDER = [
        "global",
        "preprocess", "cc", "filter", "stack", "refstack",
        "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
        "stretching", "stretching_dvv",
        "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
        "psd", "psd_rms",  # ← PSD branch: level=1 like preprocess
    ]

    result = []
    seen = set()
    for cat in _UI_ORDER:
        if cat not in chains:
            continue
        entry = chains[cat]
        if isinstance(entry, dict):
            result.append({
                "category":     cat,
                "display_name": entry.get("display_name", cat.replace("_", " ").title()),
                "level":        entry.get("level", 0),
            })
        seen.add(cat)

    # Append plugin-added categories not in _UI_ORDER
    for cat, entry in chains.items():
        if cat in seen:
            continue
        if isinstance(entry, dict):
            result.append({
                "category":     cat,
                "display_name": entry.get("display_name", cat.replace("_", " ").title()),
                "level":        entry.get("level", 0),
            })

    return result


def is_terminal_category(category):
    """Return ``True`` if *category* is a terminal workflow step (no successors).

    Terminal categories are those with an empty ``next_steps`` list in
    :func:`get_workflow_chains`.  Examples: ``psd_rms``, ``mwcs_dtt_dvv``.

    :param category: Workflow category string.
    :rtype: bool
    """
    chains = get_workflow_chains()
    entry = chains.get(category, {})
    if isinstance(entry, dict):
        return entry.get("is_terminal", not entry.get("next_steps", True))
    return not entry  # flat list format: empty = terminal


def is_entry_category(category):
    """Return ``True`` if *category* is a DAG entry point (no predecessors).

    Entry categories are those with ``is_entry_point: True`` in
    :func:`get_workflow_chains`.  Currently only ``global``.

    :param category: Workflow category string.
    :rtype: bool
    """
    chains = get_workflow_chains()
    entry = chains.get(category, {})
    if isinstance(entry, dict):
        return entry.get("is_entry_point", False)
    return False


# Backward-compat aliases — code that imported these from msnoise_table_def
# (or from workflow.py directly) still works.
WORKFLOW_CHAINS = _BUILTIN_WORKFLOW_CHAINS
WORKFLOW_ORDER  = _BUILTIN_WORKFLOW_ORDER



def get_workflow_steps(session):
    """Get all steps in a workflow"""
    return (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .order_by(WorkflowStep.step_name)
        .all()
    )



def get_workflow_links(session):
    """Get all links in a workflow"""
    return (
        session.query(WorkflowLink)
        .filter(WorkflowLink.is_active.is_(True))
        .all()
    )



def get_workflow_graph(session):
    """Return workflow as nodes and edges for visualization, sorted by workflow order"""
    steps = get_workflow_steps(session)
    links = get_workflow_links(session)

    # Sort steps by workflow order (WORKFLOW_ORDER imported from msnoise_table_def)
    _wf_order = get_workflow_order()
    def get_workflow_order_key(step):
        try:
            category_order = _wf_order.index(step.category)
        except ValueError:
            category_order = len(_wf_order)
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
            "type": link.link_type,
        })

    # No virtual edges needed — stack→mwcs links are now real WorkflowLinks
    # created by create_workflow_links_from_steps (set-number matched).

    return {"nodes": nodes, "edges": edges}



def create_workflow_step(session, step_name, category, set_number, description=None):
    """Create a new workflow step"""
    step = WorkflowStep(
        step_name=step_name, category=category,
        set_number=set_number, description=description,
    )
    session.add(step)
    session.commit()
    return step



def create_workflow_link(session, from_step_id, to_step_id, link_type="default"):
    """Create a link between two workflow steps"""
    existing = session.query(WorkflowLink).filter(
        WorkflowLink.from_step_id == from_step_id,
        WorkflowLink.to_step_id == to_step_id,
    ).first()
    if existing:
        return existing
    link = WorkflowLink(from_step_id=from_step_id, to_step_id=to_step_id, link_type=link_type)
    session.add(link)
    session.commit()
    return link



def get_step_successors(session, step_id):
    """Get all steps that this step feeds into"""
    return (
        session.query(WorkflowStep)
        .join(WorkflowLink, WorkflowStep.step_id == WorkflowLink.to_step_id)
        .filter(WorkflowLink.from_step_id == step_id)
        .filter(WorkflowLink.is_active.is_(True))
        .all()
    )



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
    return (
        session.query(WorkflowStep)
        .join(WorkflowLink, WorkflowStep.step_id == WorkflowLink.from_step_id)
        .filter(WorkflowLink.to_step_id == step_id)
        .filter(WorkflowLink.is_active.is_(True))
        .all()
    )



def create_workflow_steps_from_config_sets(session):
    """
    Create workflow steps automatically from all existing config sets,
    sorted by natural workflow order.

    Returns:
        tuple: (created_count, existing_count, error_message)
    """
    try:
        config_sets = session.query(
            Config.category, Config.set_number,
        ).filter(Config.set_number.isnot(None)).distinct().all()

        _order = get_workflow_order()
        def _key(cs):
            try:
                return _order.index(cs[0])
            except ValueError:
                return len(_order)

        config_sets = sorted(config_sets, key=_key)
        created_count = 0
        existing_count = 0

        for category, set_number in config_sets:
            existing_step = session.query(WorkflowStep).filter(
                WorkflowStep.category == category,
                WorkflowStep.set_number == set_number,
            ).first()

            if not existing_step:
                create_workflow_step(
                    session, f"{category}_{set_number}", category, set_number,
                    f"Auto-generated step for {category} configuration set {set_number}",
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
    following the DAG defined by :func:`get_workflow_chains`.

    **Linking rule**: every step in a source category is linked to every step
    in each of its successor categories (ALL×ALL).  Set numbers carry no
    topological meaning — all combinations are created and the user removes
    unwanted links via the admin UI.  The DAG is determined entirely by
    ``WorkflowLink`` rows.

    Returns:
        tuple: (created_count, existing_count, error_message)
    """
    try:
        steps = session.query(WorkflowStep).all()
        steps_by_category: dict[str, list] = {}
        for step in steps:
            steps_by_category.setdefault(step.category, []).append(step)

        created_count = 0
        existing_count = 0

        for source_category, chain_info in get_workflow_chains().items():
            target_categories = chain_info.get('next_steps', []) if isinstance(chain_info, dict) else chain_info
            if source_category not in steps_by_category:
                continue
            for source_step in steps_by_category[source_category]:
                for target_category in target_categories:
                    if target_category not in steps_by_category:
                        continue
                    for target_step in steps_by_category[target_category]:
                        existing_link = session.query(WorkflowLink).filter(
                            WorkflowLink.from_step_id == source_step.step_id,
                            WorkflowLink.to_step_id == target_step.step_id,
                        ).first()
                        if not existing_link:
                            create_workflow_link(session, source_step.step_id, target_step.step_id, 'default')
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
    try:
        row = Lineage(lineage_str=lineage_str)
        session.add(row)
        session.flush()
        return row.lineage_id
    except Exception:
        # Another worker inserted it concurrently — roll back and re-query
        session.rollback()
        with session.no_autoflush:
            row = session.query(Lineage).filter(
                Lineage.lineage_str == lineage_str).first()
        return row.lineage_id if row else None


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


def query_active_steps(session, categories):
    """Return all active WorkflowStep rows whose category is in *categories*.

    :param categories: A single category string or an iterable of strings.
    :rtype: list[WorkflowStep]
    """
    if isinstance(categories, str):
        categories = [categories]
    categories = list(categories)
    if len(categories) == 1:
        return (
            session.query(WorkflowStep)
            .filter(WorkflowStep.is_active.is_(True))
            .filter(WorkflowStep.category == categories[0])
            .all()
        )
    return (
        session.query(WorkflowStep)
        .filter(WorkflowStep.is_active.is_(True))
        .filter(WorkflowStep.category.in_(categories))
        .all()
    )


def resolve_lineage_ids(session, lids):
    """Batch-resolve a set of lineage_ids to their lineage strings.

    :param lids: Iterable of integer lineage_ids (``None`` values are skipped).
    :returns: ``{lineage_id: lineage_str}`` for every id that exists in the DB.
    :rtype: dict[int, str]
    """
    lids = {lid for lid in lids if lid is not None}
    if not lids:
        return {}
    rows = (
        session.query(Lineage.lineage_id, Lineage.lineage_str)
        .filter(Lineage.lineage_id.in_(lids))
        .all()
    )
    return {r.lineage_id: r.lineage_str for r in rows}


def bulk_upsert_jobs(session, to_insert, to_bump_refs, now, *, bump_flag_filter=None):
    """Bulk-insert new jobs and bump existing ones back to ``"T"``.

    :param to_insert: List of job dicts (keys: day, pair, jobtype, step_id,
        priority, flag, lastmod, lineage_id).
    :param to_bump_refs: List of ``Job.ref`` primary-key values to reset to T.
    :param now: ``datetime`` used for ``lastmod`` on bumped rows.
    :param bump_flag_filter: If given, only rows whose current flag equals this
        value are bumped.  Pass ``None`` to bump unconditionally.
    :returns: ``(n_inserted, n_bumped)``
    :rtype: tuple[int, int]
    """
    from sqlalchemy import update as _sa_update

    created = 0
    bumped = 0
    _CHUNK = 900

    if to_insert:
        for chunk_start in range(0, len(to_insert), _CHUNK):
            chunk = to_insert[chunk_start:chunk_start + _CHUNK]
            try:
                session.bulk_insert_mappings(Job, chunk)
                session.flush()
            except Exception:
                session.rollback()
        created = len(to_insert)

    if to_bump_refs:
        for i in range(0, len(to_bump_refs), _CHUNK):
            chunk = to_bump_refs[i : i + _CHUNK]
            stmt = (
                _sa_update(Job)
                .where(Job.ref.in_(chunk))
                .values(flag="T", lastmod=now)
                .execution_options(synchronize_session=False)
            )
            if bump_flag_filter is not None:
                stmt = stmt.where(Job.flag == bump_flag_filter)
            session.execute(stmt)
        bumped = len(to_bump_refs)

    return created, bumped


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

    if category == "cc":
        # cc Done → create both stack AND refstack T jobs (siblings under filter).
        # Stack: one per (pair, day). Refstack: one REF sentinel per pair.
        from ..s02_new_jobs import (propagate_stack_jobs_from_cc_done,
                                     propagate_refstack_jobs_from_cc_done)
        total  = propagate_stack_jobs_from_cc_done(session)
        total += propagate_refstack_jobs_from_cc_done(session)
        return total

    if category == "stack":
        # Stack can feed mwcs either:
        #   a) via refstack (join gate: both Done stack days + Done refstack REF)
        #   b) directly (no refstack): stack Done → mwcs T
        # Both are triggered here; each function checks WorkflowLinks to decide
        # which paths apply.
        from ..s02_new_jobs import (propagate_mwcs_jobs_from_refstack_done,
                                     propagate_mwcs_jobs_from_stack_done)
        total  = propagate_mwcs_jobs_from_refstack_done(session)
        total += propagate_mwcs_jobs_from_stack_done(session)
        return total

    if category == "refstack":
        # REF Done — create mwcs/stretching/wavelet for all (pair, day) combos
        # where Done stack jobs also exist (join logic inside the function).
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

    completed_names = batch["lineage_names"]
    days  = list(dict.fromkeys(batch["days"]))
    pairs = list(dict.fromkeys(j.pair for j in batch["jobs"]))

    PASSTHROUGH = frozenset({"filter", "global"})

    successors = (
        session.query(WorkflowStep)
        .join(WorkflowLink, WorkflowStep.step_id == WorkflowLink.to_step_id)
        .filter(WorkflowLink.from_step_id == completed_step.step_id)
        .filter(WorkflowLink.is_active.is_(True))
        .all()
    )

    now = datetime.datetime.now(datetime.timezone.utc)

    for succ in successors:
        worker_targets: list = []

        def _collect(step, intermediates):
            if step.category in PASSTHROUGH:
                nxt = (
                    session.query(WorkflowStep)
                    .join(WorkflowLink, WorkflowStep.step_id == WorkflowLink.to_step_id)
                    .filter(WorkflowLink.from_step_id == step.step_id)
                    .filter(WorkflowLink.is_active.is_(True))
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
            to_bump   = []
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
                        to_bump.append(ej.ref)

            n_ins, n_bump = bulk_upsert_jobs(session, to_insert, to_bump, now)
            total += n_ins + n_bump

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

    Note: If session becomes invalid, will attempt to get a fresh session from the engine.

    :type session: Session
    :param session: the database connection object
    :type jobs: list or tuple
    :param jobs: a list of :class:`~msnoise.msnoise_table_def.declare_tables.Job` to update.
    :type flag: str
    :param flag: The destination flag.
    """
    from sqlalchemy.exc import OperationalError, InvalidRequestError

    updated = False
    mappings = [{'ref': job.ref, 'flag': flag} for job in jobs]
    max_retries = 10
    retry_count = 0
    session_reconnected = False
    active_session = session

    while not updated and retry_count < max_retries:
        try:
            active_session.bulk_update_mappings(Job, mappings)
            active_session.commit()
            updated = True

        except OperationalError as e:
            # Database locked or temporary connection issue - retry with backoff
            retry_count += 1
            if retry_count >= max_retries:
                raise RuntimeError(
                    f"Failed to update jobs after {max_retries} retries due to database lock/connection issues"
                ) from e
            time.sleep(np.random.random() * (2 ** retry_count))  # exponential backoff

            # After a few retries, try to rollback the session
            if retry_count >= 3:
                try:
                    active_session.rollback()
                except Exception:
                    pass  # rollback failed, continue with retries

        except (InvalidRequestError, Exception) as e:
            # Session is invalid/closed or other error - try to get a fresh session from engine
            if not session_reconnected:
                try:
                    # Get the engine from the current session
                    engine = active_session.get_bind()

                    # Close the old session
                    try:
                        active_session.close()
                    except Exception:
                        pass

                    # Create a new session from the engine
                    from sqlalchemy.orm import sessionmaker
                    Session = sessionmaker(bind=engine)
                    active_session = Session()
                    session_reconnected = True
                    retry_count += 1

                    if retry_count >= max_retries:
                        raise RuntimeError(
                            f"Failed to update jobs after session reconnection and {max_retries} retries"
                        ) from e

                    time.sleep(1.0)
                    continue

                except Exception as recovery_error:
                    raise RuntimeError(
                        f"Cannot reconnect to database: {recovery_error}"
                    ) from e
            else:
                # Already reconnected session, this is a permanent failure
                raise RuntimeError(
                    f"Error persists after session reconnection: {type(e).__name__}: {e}"
                ) from e

    if not updated:
        raise RuntimeError(
            f"Failed to update jobs after {max_retries} retries"
        )

    # Clean up: if we created a new session, close it
    if session_reconnected:
        try:
            active_session.close()
        except Exception:
            pass

    return



def reset_jobs(session, jobtype, alljobs=False, reset_i=True, reset_e=True,
               downstream=False):
    """Reset jobs back to ``"T"``odo.

    *jobtype* can be:

    - A specific step name, e.g. ``"cc_1"`` — resets only that step.
    - A category name, e.g. ``"cc"`` — resets **all** active steps in that
      category (``cc_1``, ``cc_2``, …).

    If *downstream* is ``True``, all steps reachable via ``WorkflowLink``
    from the resolved step(s) are also reset (breadth-first, respects
    ``is_active``).

    :type session: :class:`sqlalchemy.orm.session.Session`
    :param session: A :class:`~sqlalchemy.orm.session.Session` object
    :type jobtype: str
    :param jobtype: Step name **or** category name to reset.
    :type alljobs: bool
    :param alljobs: If True reset all jobs regardless of current flag;
        otherwise only resets ``"I"`` and/or ``"F"`` flagged jobs.
    :type reset_i: bool
    :param reset_i: Reset ``"I"``n-progress jobs (default True).
    :type reset_e: bool
    :param reset_e: Reset ``"E"``rror/failed jobs (default True).
    :type downstream: bool
    :param downstream: Also reset every step reachable downstream via active
        ``WorkflowLink`` rows (default False).
    """
    from sqlalchemy import update as sa_update

    # -- Resolve jobtype to a list of step objects ----------------------------
    all_steps = get_workflow_steps(session)
    step_by_name = {s.step_name: s for s in all_steps}
    steps_by_cat = {}
    for s in all_steps:
        steps_by_cat.setdefault(s.category, []).append(s)

    if jobtype in step_by_name:
        seed_steps = [step_by_name[jobtype]]
    elif jobtype in steps_by_cat:
        seed_steps = steps_by_cat[jobtype]
    else:
        raise ValueError(
            f"Unknown step name or category: {jobtype!r}. "
            f"Expected e.g. 'cc_1' (step) or 'cc' (category)."
        )

    # -- Optionally expand to all downstream steps ----------------------------
    if downstream:
        all_links = get_workflow_links(session)
        link_map: dict = {}
        for lk in all_links:
            link_map.setdefault(lk.from_step_id, []).append(lk.to_step_id)
        step_by_id = {s.step_id: s for s in all_steps}

        visited: set = set()
        queue = list(seed_steps)
        target_steps = []
        while queue:
            step = queue.pop(0)
            if step.step_id in visited:
                continue
            visited.add(step.step_id)
            target_steps.append(step)
            for child_id in link_map.get(step.step_id, []):
                child = step_by_id.get(child_id)
                if child and child.step_id not in visited:
                    queue.append(child)
    else:
        target_steps = seed_steps

    # -- Reset flags ----------------------------------------------------------
    flags = []
    if not alljobs:
        if reset_i:
            flags.append('I')
        if reset_e:
            flags.append('F')

    for step in target_steps:
        q = sa_update(Job).where(Job.jobtype == step.step_name)
        if alljobs:
            session.execute(q.values(flag='T'))
        elif flags:
            session.execute(q.where(Job.flag.in_(flags)).values(flag='T'))

    session.commit()
    return [s.step_name for s in target_steps]



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
        chunk_size=0,
):
    """
    Return a claimed batch of jobs for a workflow step category.

    .. note::
        Most callers should use :func:`get_next_lineage_batch` instead, which
        wraps this function and also resolves lineage strings, params, and
        station-pair metadata into a ready-to-use batch dict.

    group_by:
      - "day": claim all jobs for the selected (step_id, jobtype, day)
      - "pair": claim all jobs for the selected (step_id, jobtype, pair)
      - "pair_lineage": claim all jobs for the selected (step_id, jobtype, pair, lineage)
      - "day_lineage": claim all jobs for the selected (step_id, jobtype, day, lineage)

    chunk_size:
      Only applies to ``group_by="day_lineage"``.  When > 0, at most *chunk_size*
      jobs are claimed per batch instead of the full day.  Useful for steps with
      O(N²) work per day (CC) or large per-day station counts (PSD) so that
      multiple workers can share the same day in parallel without write conflicts.

      ``chunk_size=0`` (default) claims everything — identical to the original
      behaviour.
    """
    from sqlalchemy import update

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
        if chunk_size > 0:
            # Claim at most chunk_size jobs for this day.  The seed job is
            # included in the count, so we limit to chunk_size total.
            # Deterministic pair ordering ensures different workers claim
            # non-overlapping subsets (each grabs the next unclaimed block).
            batch_q = batch_q.limit(int(chunk_size))
    else:
        raise ValueError(f"Unsupported group_by={group_by!r}")

    jobs = batch_q.with_for_update().all()
    if not jobs:
        return [], step

    # Collect the ref IDs that were actually selected so the UPDATE is
    # scoped to exactly these rows — critical for chunk_size correctness.
    claimed_refs = [j.ref for j in jobs]

    upd = (
        update(Job)
        .where(Job.ref.in_(claimed_refs))
        .values(flag="I")
        .execution_options(synchronize_session=False)
    )
    session.execute(upd)
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
    if step_category == "refstack":
        day_filter = (Job.day == "REF")
    else:
        day_filter = (Job.day != "REF")

    return session.query(Job) \
        .join(WorkflowStep, Job.jobtype == WorkflowStep.step_name) \
        .filter(WorkflowStep.category == step_category) \
        .filter(Job.flag == flag) \
        .filter(day_filter) \
        .first() is not None



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
    # get_step_successors is defined in this module — no import needed

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
    :rtype: tuple(list, list[str], MSNoiseParams)
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
        chunk_size=0,
):
    """
    Standard worker prolog for lineage-aware steps.

    - Claims the next batch of jobs for `step_category` using `get_next_job_for_step`.
    - Extracts (pair, lineage_str, refs, days).
    - Loads current step config (from the job row).
    - Resolves lineage_str -> WorkflowStep objects.
    - Builds a MSNoiseParams for that lineage (one layer per category).

    Parameters
    ----------
    chunk_size : int, optional
        Passed to :func:`get_next_job_for_step`.  Only effective for
        ``group_by="day_lineage"``.  When > 0, at most *chunk_size* jobs are
        claimed per batch, enabling multiple workers to share the same day in
        parallel.  Default ``0`` = claim everything (original behaviour).

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
      - step_params, params    — params is a MSNoiseParams instance
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
        chunk_size=chunk_size,
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
    # lineage_names_ref: path to the REF file folder.
    # - If refstack is in the lineage (convention A1: …/stack_N/refstack_M/mwcs_1):
    #   strip stack_* → resolves to …/filter_N/refstack_M/_output/REF/
    # - If no refstack in lineage (direct stack→mwcs: …/stack_N/mwcs_1):
    #   set to None — workers must skip REF loading for this topology.
    _has_refstack = any(n.startswith("refstack_") for n in lineage_names_upstream)
    lineage_names_ref = (
        _strip_stack_from_lineage(lineage_names_upstream) if _has_refstack else None
    )

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
        "lineage_names_ref": lineage_names_ref,
        "refs": refs,
        "days": days,
        "step_params": step_params,
        "params": params,
    }



def get_refstack_lineage_for_filter(session, filterid, refstack_set_number=1):
    """Get the full lineage path from root through a filter step to its refstack.

    Refstack is now a sibling of stack (both children of the filter pass-through).
    Returns the path ending at ``refstack_M`` directly downstream of
    ``filter_{filterid}`` — suitable for :func:`xr_get_ref`.

    Example for the default single-pipeline::

        get_refstack_lineage_for_filter(db, 1)
        # → ['preprocess_1', 'cc_1', 'filter_1', 'refstack_1']

    :type session: :class:`sqlalchemy.orm.session.Session`
    :type filterid: int
    :param filterid: The filter set_number (e.g. 1 for filter_1).
    :type refstack_set_number: int
    :param refstack_set_number: Which refstack set to use (default 1).
    :rtype: list of str
    """
    filter_step, step_map, links, path = _walk_upstream_from_filter(session, filterid)
    if filter_step is None or not path:
        return path or []

    # Find refstack steps directly linked from this filter step.
    refstack_steps = [
        step_map[lk.to_step_id]
        for lk in links
        if lk.from_step_id == filter_step.step_id
        and lk.to_step_id in step_map
        and step_map[lk.to_step_id].category == 'refstack'
    ]
    if not refstack_steps:
        return path  # no refstack in this workflow

    refstack_step = next(
        (s for s in refstack_steps if s.set_number == refstack_set_number),
        refstack_steps[0],
    )
    return path + [refstack_step.step_name]


def _walk_upstream_from_filter(session, filterid):
    """Walk upstream from ``filter_{filterid}`` to the DAG root.

    Returns ``(filter_step, step_map, links, path)`` where *path* is the
    ordered list of step-name strings from root → filter (global nodes
    excluded).  Returns ``(None, {}, [], [])`` if the filter step is not found.

    This is the shared kernel for :func:`_get_filter_lineage` and
    :func:`_get_stack_lineage_for_filter`.
    """
    steps = get_workflow_steps(session)
    links = get_workflow_links(session)
    step_map = {s.step_id: s for s in steps}

    filter_step = next(
        (s for s in steps if s.category == 'filter' and s.set_number == filterid),
        None,
    )
    if filter_step is None:
        return None, {}, [], []

    parent_map = {link.to_step_id: link.from_step_id for link in links}
    path = [filter_step.step_name]
    current_id = filter_step.step_id
    visited = set()
    while current_id in parent_map:
        if current_id in visited:
            break
        visited.add(current_id)
        parent_id = parent_map[current_id]
        parent_step = step_map[parent_id]
        if parent_step.category != 'global':
            path.insert(0, parent_step.step_name)
        current_id = parent_id

    return filter_step, step_map, links, path


def _get_filter_lineage(session, filterid):
    """Walk upstream from filter_{filterid} to root; return path including filter.

    Unlike :func:`_get_stack_lineage_for_filter`, this does NOT append a
    downstream step — it stops at the filter node.  Used by
    :func:`get_refstack_lineage_for_filter` now that refstack hangs directly
    off filter rather than off stack.
    """
    _, _, _, path = _walk_upstream_from_filter(session, filterid)
    return path

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



def refstack_needs_recompute(session, pair, cc_lineage_prefix, params):
    """Return True if the REF stack needs to be recomputed for this pair.

    Queries Done ``stack`` jobs for *pair* across all active stack_N steps
    whose lineage starts with *cc_lineage_prefix* (e.g.
    ``["preprocess_1", "cc_1", "filter_1"]``), and checks whether any date
    falls inside the configured ``[ref_begin, ref_end]`` window.

    Stack and refstack are now siblings (both children of filter/cc), so the
    refstack worker no longer has a stack step in its upstream lineage.
    Instead, the cc/filter prefix is shared by both branches and is used here
    to locate the sibling stack jobs.

    Should only be called for Mode A (fixed-date) refstack jobs.

    :param session: SQLAlchemy session.
    :param pair: Station pair string ``"NET.STA.LOC:NET.STA.LOC"``.
    :param cc_lineage_prefix: Lineage name list up to (but not including) the
        stack/refstack level, e.g. ``["preprocess_1", "cc_1", "filter_1"]``.
        Pass ``batch["lineage_names_upstream"]`` from the refstack worker —
        that list ends at filter_N (the pass-through above refstack).
    :param params: :class:`~msnoise.params.MSNoiseParams` for this lineage.
    :rtype: bool
    """
    import datetime as _dt

    ref_start, ref_end, _ = build_ref_datelist(params, session)

    stack_steps = query_active_steps(session, "stack")
    if not stack_steps:
        return True  # no stack steps — assume recompute needed

    # Build expected lineage prefix for each stack_N:
    # e.g. "preprocess_1/cc_1/filter_1/stack_1"
    prefix_str = "/".join(cc_lineage_prefix) if cc_lineage_prefix else ""

    for stack_step in stack_steps:
        expected_lineage = (
            f"{prefix_str}/{stack_step.step_name}" if prefix_str
            else stack_step.step_name
        )
        stack_lin_id = _lineage_id_for(session, expected_lineage)
        if stack_lin_id is None:
            continue  # this stack_N has no jobs with this lineage

        done_days = (
            session.query(Job.day)
            .filter(Job.step_id == stack_step.step_id)
            .filter(Job.pair == pair)
            .filter(Job.flag == "D")
            .filter(Job.day != "REF")
            .filter(Job.lineage_id == stack_lin_id)
            .all()
        )

        for (day_str,) in done_days:
            try:
                day = _dt.datetime.strptime(day_str, "%Y-%m-%d").date()
            except ValueError:
                continue
            if ref_start <= day <= ref_end:
                return True

    return False



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



def _strip_stack_from_lineage(lineage_names):
    """Return a copy of ``lineage_names`` with any ``stack_*`` entries removed.

    Used to build ``lineage_names_ref`` — the path for :func:`xr_get_ref`,
    which lives under the ``refstack_M`` step folder.  In convention A1 the
    mwcs lineage string encodes both parents as ``…/stack_N/refstack_M/mwcs_1``,
    but the REF file is written by the refstack worker under
    ``…/filter_N/refstack_M/_output/REF/`` — i.e. without the ``stack_N``
    prefix that appears in the mwcs lineage.

    Example::

        _strip_stack_from_lineage(
            ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']
        )
        # → ['preprocess_1', 'cc_1', 'filter_1', 'refstack_1']

    :type lineage_names: list of str
    :rtype: list of str
    """
    return [n for n in lineage_names if not n.startswith("stack_")]


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
    filter_step, step_map, links, path = _walk_upstream_from_filter(session, filterid)
    if filter_step is None:
        return []

    # Append the stack step immediately downstream of this filter (if any)
    for link in links:
        if link.from_step_id == filter_step.step_id:
            child = step_map.get(link.to_step_id)
            if child and child.category == 'stack':
                path.append(child.step_name)
                break

    return path
