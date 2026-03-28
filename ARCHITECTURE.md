# MSNoise Architecture Reference

*Self-reference for future Claude sessions. Read this before touching workflow, jobs, or lineage.*

---

## 1. The Workflow Graph

MSNoise uses a directed acyclic graph (DAG) of **WorkflowSteps** connected by **WorkflowLinks**.

```
global_1 ──► preprocess_1 ──► cc_1 ──► filter_1 ──► stack_1 ──► refstack_1
                                                                      │
         ──► psd_1 ──► psd_rms_1                           ┌─────────┼─────────┐
                                                            ▼         ▼         ▼
                                                          mwcs_1  stretch_1  wavelet_1
                                                            │         │         │
                                                        dtt_1     str_dvv  wct_dtt_1
                                                            │         │         │
                                                        dvv_1     (done)   wct_dvv_1
```

**Pass-through nodes**: `filter_N` and `global_N` appear in lineage strings but have **no worker scripts** and **no jobs**. They are purely parameter namespaces. `propagate_downstream` must recurse through them transparently.

---

## 2. Jobs — Two Distinct Concepts

### `Job` (ORM object, `msnoise_table_def.declare_tables().Job`)
SQLAlchemy mapped class. **Always call `declare_tables()` to get a fresh class** — do NOT use module-level imports like `from msnoise_table_def import Job` for DB queries, as multiple `declare_tables()` calls create different class objects. Use `schema = declare_tables(); Job = schema.Job`.

Key columns: `ref` (PK), `day`, `pair`, `flag` (T/I/D/F), `step_id` (FK→WorkflowStep), `lineage_id` (FK→Lineage), `jobtype` (= step_name), `priority`, `lastmod`.

### Job flags lifecycle
```
T (Todo) → I (In Progress, claimed by worker) → D (Done)
                                              → F (Failed)
D → T  only allowed explicitly (reset) or by propagate_downstream (upstream changed)
```

### Special job types
- **Preprocess/PSD jobs**: `pair` = single station `"YA.UV05.00"`, `day` = date string
- **CC/stack/mwcs jobs**: `pair` = station pair `"YA.UV05.00:YA.UV06.00"`, `day` = date string  
- **Refstack REF jobs**: `pair` = station pair, `day = "REF"` — sentinel, one per pair
- **DVV sentinel jobs**: `pair = "ALL"`, `day = "DVV"` — aggregation trigger, one per lineage

`is_next_job_for_step(db, "refstack")` specifically filters `day == "REF"`. All others filter `day != "REF"`.

---

## 3. Lineage — The Path Through the Graph

Every job carries a `lineage_id` FK → `Lineage(lineage_str)`. The string is a `/`-separated path of step names from root to the current step:

```
preprocess_1
preprocess_1/cc_1
preprocess_1/cc_1/filter_1/stack_1      ← filter_1 is pass-through but appears in path!
preprocess_1/cc_1/filter_1/stack_1/refstack_1
preprocess_1/cc_1/filter_1/stack_1/refstack_1/mwcs_1/mwcs_dtt_1/mwcs_dtt_dvv_1
psd_1
psd_1/psd_rms_1
```

**Deduplication**: multiple jobs share the same `Lineage` row. The `before_insert` event on `Job` resolves `lineage_str` → `lineage_id` via raw SQL (not ORM session) to avoid re-entrant flush.

**Key helpers in `core/workflow.py`**:
- `_get_or_create_lineage_id(session, str)` — 4-step lookup: identity_map → session.new → DB → INSERT
- `_lineage_id_for(session, str)` — read-only, returns None if not found
- `_lineage_str_from_id(session, id)` — reverse lookup, safe after session.commit() expiry
- `get_lineages_to_step_id(session, step_id)` — DFS returning all upstream paths to a step

---

## 4. Worker Script Pattern

Every worker script follows this canonical loop:

```python
while is_next_job_for_step(db, step_category="cc"):
    batch = get_next_lineage_batch(db, "cc", group_by="day_lineage")
    if batch is None:
        break
    # ... compute and write output files using batch["params"] ...
    massive_update_job(db, batch["jobs"], flag="D")
    if not batch["params"].global_.hpc:        # hpc=False: propagate inline
        propagate_downstream(db, batch)
```

`get_next_lineage_batch` returns a dict with:
- `jobs`, `step`, `pair`, `lineage_str`, `lineage_names`
- `lineage_names_upstream` = lineage_names[:-1] (for output paths)
- `lineage_names_mov` = upstream with refstack_* stripped (for reading MOV CCFs)
- `days`, `refs`
- `params` — LayeredParams (access global config as `params.global_.hpc`)

**`group_by` parameter**:
- `"day_lineage"` — batch all pairs for a given day+lineage (preprocess, cc, psd)
- `"pair_lineage"` — batch all days for a given pair+lineage (stack, refstack, mwcs, stretching, wct)

---

## 5. `propagate_downstream` — The Dispatch Table

After marking jobs Done in hpc=False mode, `propagate_downstream(session, batch)` fires:

| Completed category | Delegates to | Reason |
|---|---|---|
| `preprocess` | `create_cc_jobs_from_preprocess` | Fan-out: 1 station → N pairs |
| `stack` | `propagate_refstack_jobs_from_stack_done` | Creates `day="REF"` sentinel per pair |
| `refstack` | `propagate_mwcs_jobs_from_refstack_done` | Looks up MOV stack days → creates per-(pair,day) mwcs/str/wct jobs |
| `mwcs_dtt`/`stretching`/`wavelet_dtt` | `propagate_dvv_jobs_from_dtt_done` | Creates `day="DVV", pair="ALL"` sentinel |
| `psd` | `propagate_psd_rms_jobs_from_psd_done` | Simple single-station passthrough |
| `cc` (and others) | Generic pair×day bulk upsert | Crosses filter_1 passthrough via `_collect_workers` recursion |

**HPC mode** (`hpc=Y`): `propagate_downstream` is NOT called by workers. Operator runs `msnoise new_jobs --after X` manually. In hpc=False mode, `new_jobs --after X` still works as a reconciliation/safety pass (logs debug message).

---

## 6. Output File Paths

All output paths follow: `root / lineage_names[:-1] / step_name / _output / ...`

Note: `lineage_names[:-1]` = **upstream** of the current step (not including it). `step_name` = last element. Together they form the full path:

```
OUTPUT/preprocess_1/cc_1/filter_1/_output/daily/ZZ/YA.UV05.00_YA.UV06.00/2010-09-01.nc
       ─────────────────────────── ────────────
       lineage_names_upstream       step_name
```

DVV aggregates:
```
OUTPUT / upstream_lineage / dvv_step_name / _output / 1D_1D / dvv_CC_ZZ.nc
```

---

## 7. The REF Sentinel Design

`refstack` uses `day="REF"` (one job per pair, not per day). This is intentional:

1. Stack_N completes → `propagate_refstack_jobs_from_stack_done` creates one REF job per pair
2. Refstack worker processes REF job — computes new reference from ALL available stack days
3. After REF Done → `propagate_mwcs_jobs_from_refstack_done` re-queries stack table for all Done MOV days → creates one mwcs/stretching/wct job per (pair, day)

**Why**: any new stack day potentially changes the reference → all pair comparisons must be recomputed. The REF job is a synchronization point, not a per-day operation.

---

## 8. Config & Params

Config lives in DB tables, one `ConfigSet` per category (global, preprocess, cc, filter, stack, refstack, mwcs, etc.).

- `get_params(db)` → flat `AttribDict` (all global config). Access `params.hpc`, `params.cc_sampling_rate`, etc.
- `get_merged_params_for_lineage(db, ...)` → `LayeredParams`. Access per-category via `params.global_.hpc`, `params.mwcs.mwcs_wlen`, `params["stack"].mov_stack`

**`LayeredParams` access**: `params.global_` (note underscore — avoids Python keyword), `params.cc`, `params.mwcs` etc. Raises `AttributeError` for unknown categories.

---

## 9. Common Pitfalls

1. **`declare_tables()` class identity**: `Job` from one call ≠ `Job` from another. Always call `declare_tables()` locally, or use the session's mapped class via `schema = declare_tables()`.

2. **Session cache staleness**: After commits from other sessions (e.g., `massive_update_job`), call `db.expire_all()` OR use a fresh `connect()` for queries in tests.

3. **`propagate_downstream` is idempotent**: calling it twice for the same batch is safe — existing T/I jobs are left unchanged, existing D/F jobs are bumped to T.

4. **`new_jobs --after X` in hpc=False**: logs a debug message and runs anyway (reconciliation). Does NOT block. This preserves backward compatibility with scripts and tests that call it explicitly.

5. **`_lineage_through(cat)` vs `_lineage_upstream_of(cat)`**: DVV files are saved with `upstream` lineage + step_name. Use `_lineage_upstream_of(dvv_cat)` in `get_dvv()` to match the save path.
