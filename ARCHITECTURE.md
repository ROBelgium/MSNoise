# MSNoise Architecture Reference

*Self-reference for future Claude sessions. Read this before touching workflow, jobs, or lineage.*

---

## 1. The Workflow Graph

MSNoise uses a directed acyclic graph (DAG) of **WorkflowSteps** connected by **WorkflowLinks**.

```
global_1 ──► preprocess_1 ──► cc_1 ──► filter_1 ──► stack_1
                                                         │
         ──► psd_1 ──► psd_rms_1          ┌─────────────┼─────────────┐
                                           │             │             │
                                      REF sentinel   (direct)     (direct)
                                           │             │             │
                                      refstack_1     mwcs_1      stretching_1
                                           │         wavelet_1       │
                                    (re-propagates)     │        str_dvv_1
                                                     dtt_1
                                                        │
                                                    dvv_1
```

**Pass-through nodes**: `filter_N` and `global_N` appear in lineage strings but have **no worker scripts** and **no jobs**. They are purely parameter namespaces. `propagate_downstream` recurses through them transparently via `_collect()`.

**Canonical category order** (from `WORKFLOW_ORDER` in `msnoise_table_def.py`):
`global → preprocess → cc → psd → psd_rms → filter → stack → refstack → mwcs → mwcs_dtt → mwcs_dtt_dvv → stretching → stretching_dvv → wavelet → wavelet_dtt → wavelet_dtt_dvv`

---

## 2. Package Layout & Import Paths

```
msnoise/
  __init__.py            # MSNoiseError, DBConfigNotFoundError, FatalError; re-exports connect
  api.py                 # Backward-compat re-export of msnoise.core.*
  core/
    __init__.py          # Re-exports all 6 submodules via *
    db.py                # connect, get_engine, read_db_inifile, get_logger
    config.py            # get_config, update_config, create_config_set, get_params,
                         # get_merged_params_for_lineage, build_plot_outfile,
                         # lineage_to_plot_tag, get_config_categories_definition
    stations.py          # get_stations, get_station_pairs, update_data_availability,
                         # get_data_availability, add_data_source, get_data_source,
                         # resolve_data_source, get_waveform_path, import_stationxml,
                         # set_station_source, set_network_source, set_all_stations_source
    workflow.py          # get_next_lineage_batch, propagate_downstream,
                         # is_next_job_for_step, massive_insert_job, massive_update_job,
                         # reset_jobs, build_ref_datelist, get_lineages_to_step_id,
                         # _get_or_create_lineage_id, create_workflow_steps_from_config_sets,
                         # create_workflow_links_from_steps, get_t_axis,
                         # get_filter_steps_for_cc_step, get_refstack_lineage_for_filter
    signal.py            # winsorizing, stack, xwt, compute_wct_dtt, get_wct_avgcoh,
                         # preload_instrument_responses, save_preprocessed_streams,
                         # get_preprocessed_stream, make_same_length, validate_stack_data
    io.py                # xr_save_ccf, xr_get_ccf, xr_save_ccf_daily, xr_get_ccf_daily,
                         # xr_save_ref, xr_get_ref, xr_load_ccf_for_stack,
                         # xr_save_mwcs, xr_get_mwcs, xr_save_dtt, xr_get_dtt,
                         # xr_save_stretching, xr_save_wct, xr_load_wct,
                         # xr_save_wct_dtt, xr_get_wct_dtt,
                         # xr_save_dvv_agg, xr_get_dvv_agg, aggregate_dvv_pairs,
                         # xr_save_psd, xr_load_psd, xr_save_rms, xr_load_rms,
                         # psd_rms, psd_df_rms, psd_read_results
    fdsn.py              # parse_datasource_scheme, is_remote_source, build_client,
                         # fetch_waveforms_bulk, fetch_and_preprocess, _write_raw_cache
  msnoise_table_def.py   # declare_tables() → Config, Lineage, Job, DataSource,
                         # Station, DataAvailability, WorkflowStep, WorkflowLink
                         # WORKFLOW_CHAINS dict, WORKFLOW_ORDER list
  params.py              # LayeredParams, _build_layered_params
  results.py             # MSNoiseResult (high-level result reader)
  data_structures.py     # data_structure dict (SDS/BUD/IDDS/PDF path templates)
  move2obspy.py          # myCorr2, pcc_xcorr, whiten, whiten2, mwcs, smooth
  preprocessing.py       # apply_preprocessing_to_stream, preprocess
  wiener.py              # wiener_filt, find_segments
  s01_scan_archive.py    # scan_archive main, get_archives_folders, parse_crondays
  s02_new_jobs.py        # new_jobs main; all propagate_* functions
  s02_preprocessing.py   # preprocess worker main
  s03_compute_no_rotation.py  # CC computation main
  s04_stack_mov.py       # moving stack main (stype='mov'|'ref')
  s04_stack_refstack.py  # reference stack main
  s05_compute_mwcs.py    # MWCS main
  s06_compute_mwcs_dtt.py # MWCS DTT main
  s07_compute_dvv.py     # DVV aggregate main (mwcs_dtt_dvv / stretching_dvv / wavelet_dtt_dvv)
  s08_compute_wct.py     # Wavelet coherence main
  s09_compute_wct_dtt.py # Wavelet DTT main
  s10_stretching.py      # Stretching main; stretch_mat_creation
  psd_compute.py         # PSD computation main
  psd_compute_rms.py     # PSD RMS computation main
  config/                # One config_<category>.csv per workflow category
  plots/                 # ccftime, data_availability, distance, dtt, interferogram,
                         # mwcs, mwcs_dtt_dvv, ppsd, psd_rms, spectime,
                         # station_map, stretching_dvv, timing, wavelet_dtt_dvv
  scripts/msnoise.py     # Full Click CLI; entry point msnoise = scripts.msnoise:run
  msnoise_admin.py       # Flask-Admin web UI (msnoise admin)
```

**The one universal entry point**: `from msnoise.api import connect` (or `from msnoise import connect`).
All other symbols live in `msnoise.core.*` and are re-exported through `api.py` for backward compat.

---

## 3. Jobs — Two Distinct Concepts

### `Job` (ORM object, `msnoise_table_def.declare_tables().Job`)
SQLAlchemy mapped class. **Always call `declare_tables()` to get a fresh class** — do NOT use module-level imports like `from msnoise_table_def import Job` for DB queries, as multiple `declare_tables()` calls create different class objects. Use `schema = declare_tables(); Job = schema.Job`.

Key columns: `ref` (PK), `day`, `pair`, `flag` (T/I/D/F), `step_id` (FK→WorkflowStep), `lineage_id` (FK→Lineage), `jobtype` (= step_name), `priority`, `lastmod`.

Key `@property` fields (derived from `workflow_step` relationship — require the relationship to be loaded):
- `job.config_category` → `workflow_step.category`
- `job.config_set_number` → `workflow_step.set_number`
- `job.step_name` → `workflow_step.step_name`
- `job.lineage` → resolves `lineage_id` → `lineage_str` via the `Lineage` FK

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

## 4. Lineage — The Path Through the Graph

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

## 5. Worker Script Pattern

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
- `jobs`, `step`, `pair`, `lineage_str`, `lineage_steps`
- `lineage_names` — full list including current step name
- `lineage_names_upstream` — `lineage_names[:-1]` (for output paths — excludes current step)
- `lineage_names_mov` — upstream with any `refstack_*` entries stripped (for reading MOV CCFs in mwcs/wct/stretching)
- `days`, `refs`
- `step_params` — raw `AttribDict` for the current step's config set
- `params` — `LayeredParams` (access global config as `params.global_.hpc`)

**`group_by` parameter**:
- `"day_lineage"` — batch all pairs for a given day+lineage (preprocess, cc, psd)
- `"pair_lineage"` — batch all days for a given pair+lineage (stack, refstack, mwcs, stretching, wct)

---

## 6. `propagate_downstream` — The Dispatch Table

After marking jobs Done in hpc=False mode, `propagate_downstream(session, batch)` fires.

### The stack → downstream split (KEY DESIGN DECISION)

When stack completes, **two things happen simultaneously**:

1. **REF sentinel** (`day="REF"`) created per pair → triggers refstack
2. **Direct mwcs/stretching/wavelet T jobs** created for the new (pair, day) tuples

This dual propagation enables the common observatory case where new days fall **outside** the fixed `ref_begin..ref_end` window. In that case:
- Refstack picks up the REF job, checks the window, finds no days inside → **skips recomputation**, marks REF Done immediately
- MWCS/stretching/wavelet T jobs were already waiting → start immediately on the new data

If a new day DOES fall inside the ref window, refstack recomputes and `propagate_downstream` fires again from refstack — the existing mwcs T jobs are bumped (idempotent upsert).

This mirrors the old `stack -m` (daily) + `stack -r` (occasional) workflow without requiring manual intervention.

### Full dispatch table

| Completed category | Delegates to | Reason |
|---|---|---|
| `preprocess` | `create_cc_jobs_from_preprocess` | Fan-out: 1 station → N pairs |
| `stack` | `propagate_refstack_jobs_from_stack_done` + `propagate_mwcs_jobs_from_refstack_done` | Dual: REF sentinel + direct mwcs jobs |
| `refstack` | `propagate_mwcs_jobs_from_refstack_done` | REF just changed: re-propagate all days |
| `mwcs_dtt`/`stretching`/`wavelet_dtt` | `propagate_dvv_jobs_from_dtt_done` | Creates `day="DVV", pair="ALL"` sentinel |
| `psd` | `propagate_psd_rms_jobs_from_psd_done` | Simple single-station passthrough |
| `cc` (and others) | Generic pair×day bulk upsert via `_collect()` recursion | Crosses filter_1 passthrough transparently |

### HPC mode (`hpc=Y`)

`propagate_downstream` is NOT called by workers. Operator runs `msnoise new_jobs --after X` manually. In hpc=False mode, `new_jobs --after X` still works as a reconciliation/safety pass (logs debug message + runs anyway). This preserves backward compatibility.

---

## 7. The REF Sentinel — When It Computes vs Skips

`s04_stack_refstack` checks on every REF job:

```python
ref_start, ref_end, _ = build_ref_datelist(params, db)
done_stack_days = query all Done stack days for this pair
any_in_window = any(ref_start <= day <= ref_end for day in done_stack_days)

if not any_in_window:
    # New data outside ref window — reference unchanged
    # Mark REF Done without computing. MWCS jobs already waiting.
    mark_done(); propagate_downstream(); continue

# New data inside window — recompute reference
compute_ref_stack(); write_ref_file(); mark_done(); propagate_downstream()
```

**Mode B (rolling ref)**: `ref_begin` is a negative integer (e.g. `"-5"` = last 5 days). No REF file written; reference computed on-the-fly in MWCS/stretching/WCT. REF job still exists as synchronization point but is always fast (no file I/O).

---

## 8. Output File Paths

All output paths follow: `root / lineage_names_upstream / step_name / _output / ...`

`lineage_names_upstream` = `lineage_names[:-1]` (excludes current step). `step_name` = last element.

```
OUTPUT/preprocess_1/cc_1/filter_1/_output/daily/ZZ/YA.UV05.00_YA.UV06.00/2010-09-01.nc
       ─────────────────────────── ────────────
       lineage_names_upstream       step_name (= filter_1 for cc)

OUTPUT/...stack_1/refstack_1/_output/REF/ZZ/YA.UV05.00_YA.UV06.00.nc
OUTPUT/...refstack_1/mwcs_1/_output/1D_1D/ZZ/YA.UV05.00_YA.UV06.00.nc

# DVV aggregates:
OUTPUT / upstream_lineage / dvv_step_name / _output / 1D_1D / dvv_CC_ZZ.nc
```

**Critical**: `get_dvv()` uses `_lineage_upstream_of(dvv_cat)` (not `_lineage_through`) to match the save path.

---

## 9. Config & Params

Config lives in the unified `Config` table. One row per `(name, category, set_number)`. Global config has `category='global'`, `set_number=1` (legacy NULL rows are migrated to 1 by `get_config_sets_organized`).

Key columns on `Config`: `name`, `category`, `set_number`, `value` (always string), `param_type` ('str'/'int'/'float'/'bool'), `default_value`, `description`, `possible_values` (slash-separated).

- `get_params(db)` → single-layer `LayeredParams`. Access `params.global_.hpc`, `params.global_.output_folder`, etc.
- `get_merged_params_for_lineage(db, orig_params, step_params, lineage)` → `(lineage, lineage_names, LayeredParams)`. Access: `params.global_` (underscore avoids Python keyword), `params.cc`, `params.mwcs` etc.
- `get_config_set_details(db, category, set_number, format="AttribDict")` → plain `AttribDict` for a single config set.

**`LayeredParams` access**:
- `params.global_.hpc` — global layer
- `params["stack"].mov_stack` — dict-style also works
- `params.category_layer` — the *current step's* layer (useful in DVV aggregation functions that receive params from the worker)
- `params.step_name` — name of the innermost step
- `params.lineage_names` — full list of step names in this lineage
- `params.categories` — list of all loaded category names
- `params.as_flat_dict()` — all params flattened (later layers win on name collision)
- `params.to_yaml(path)` / `LayeredParams.from_yaml(path)` — serialization

**NOT** `params.hpc` — that's `params.global_.hpc`.

Config CSV files in `msnoise/config/config_<category>.csv` define defaults. Columns: `name`, `default`, `definition`, `type`, `possible_values`.

---

## 10. ORM Tables

All tables go through `declare_tables(prefix=None)` which returns a namespace object. **Never import ORM classes at module level for DB queries** — always:

```python
from msnoise.msnoise_table_def import declare_tables
schema = declare_tables()
Job = schema.Job
```

| Table | Key columns |
|---|---|
| `Config` | `ref`, `name`, `category`, `set_number`, `value`, `param_type`, `default_value`, `description`, `possible_values` |
| `Lineage` | `lineage_id` (PK), `lineage_str` (unique) |
| `Job` | `ref` (PK), `day`, `pair`, `flag`, `step_id` (FK→WorkflowStep), `lineage_id` (FK→Lineage), `jobtype`, `priority`, `lastmod` |
| `DataSource` | `ref` (PK), `name` (unique), `uri`, `data_structure`, `auth_env`, `archive_format`, `network_code`, `channels` |
| `Station` | `ref` (PK), `net`, `sta`, `X`, `Y`, `altitude`, `coordinates`, `used`, `data_source_id` (FK→DataSource, nullable) |
| `DataAvailability` | `ref` (PK), `net`, `sta`, `loc`, `chan`, `path`, `file`, `starttime`, `endtime`, `data_source_id` |
| `WorkflowStep` | `step_id` (PK), `step_name` (unique), `category`, `set_number`, `description`, `is_active` |
| `WorkflowLink` | `link_id` (PK), `from_step_id`, `to_step_id`, `link_type`, `is_active` |

`WorkflowStep` has a `UniqueConstraint('category', 'set_number')` — one step per config set.
`Job` has `Index('job_index', 'day', 'pair', 'step_id', 'lineage_id', unique=True)`.

---

## 11. DataSource Abstraction

Every `Station` has an optional `data_source_id` FK. NULL → use `DataSource.ref=1` (project default, created by installer).

**URI schemes** (parsed by `fdsn.parse_datasource_scheme`):
- bare path or `sds:///path` → local SDS/BUD/IDDS/PDF archive
- `fdsn://http://...` → FDSN web service
- `eida://http://...` → EIDA routing client

**`fdsn.py` key functions**:
- `build_client(ds)` → ObsPy `Client` or EIDA routing client from a `DataSource`
- `fetch_waveforms_bulk(client, bulk_request, retries=3)` → `Stream`
- `fetch_and_preprocess(...)` → fetches + preprocesses in one shot for remote sources
- `_write_raw_cache(stream, ...)` — caches raw waveforms when `fdsn_keep_raw=Y` (global config)

**Station-to-source assignment**:
- `set_station_source(db, net, sta, data_source_id)` — per-station override
- `set_network_source(db, net, data_source_id)` — all stations in a network
- `set_all_stations_source(db, data_source_id)` — project-wide
- `resolve_data_source(db, station)` — returns effective DS (station override or default)
- `get_waveform_path(db, da)` — joins `DataSource.uri + da.path + da.file`

---

## 12. MSNoiseResult — High-Level Result Reader

`MSNoiseResult` (`results.py`) is the user-facing class for reading computed results without needing to know file paths or lineage strings.

```python
# Construct from integer set IDs
r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                            stack=1, refstack=1, mwcs=1, mwcs_dtt=1,
                            mwcs_dtt_dvv=1)

# List all available result sets for a category
for r in MSNoiseResult.list(db, "mwcs_dtt"):
    print(r)

# Access results — all return xarray Dataset/DataArray by default
ds = r.get_ccf(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_ref(pair="BE.UCC--BE.MEM", components="ZZ")
ds = r.get_mwcs(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_mwcs_dtt(...)
ds = r.get_stretching(...)
ds = r.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_wct(...)
ds = r.get_wct_dtt(...)
ds = r.get_psd(seed_id="BE.UCC..HHZ", day="2023-01-01")
ds = r.get_psd_rms(seed_id="BE.UCC..HHZ")
df = r.export_dvv(...)   # returns DataFrame

# Navigate branches
for branch in r.branches():   # other lineages reachable from same root
    print(branch)
```

**Internal helpers**:
- `r._lineage_upstream_of(category)` — lineage list up to but NOT including `category`
- `r._lineage_through(category)` — lineage list including `category`
- `r._step_name_for(category)` — step_name string for a category in this lineage

**DVV path rule**: `get_dvv()` uses `_lineage_upstream_of(dvv_cat)` + step_name — NOT `_lineage_through`. Using `_lineage_through` produces a doubled path `...dvv_1/dvv_1/_output/...`.

---

## 13. DVV Aggregation (`aggregate_dvv_pairs`)

Called from `s07_compute_dvv.py`. Reads per-pair DTT/STR/WCT-DTT files and aggregates across pairs.

**Output variables** in the returned `xr.Dataset` (dim: `times`):
- Always: `mean`, `std`, `median`, `n_pairs`
- If `dvv_output_percent=Y`: values multiplied by 100 before aggregation
- If `dvv_weighted_mean=Y`: `weighted_mean`, `weighted_std` (weights = 1/err²)
- If `dvv_trimmed_mean=Y`: `trimmed_mean`, `trimmed_std` (sigma-clip at `dvv_trim_limit`)
- If `dvv_percentiles` set (comma-separated floats): `q05`, `q25`, `q50`, `q75`, `q95` etc.

**Quality filtering**: pairs with `quality_col < dvv_quality_min` are masked to NaN before aggregation.

**Sign convention**:
- MWCS: `dv/v = -1 * dt/t` (sign flip applied in aggregator)
- Stretching: `dv/v = stretching_factor - 1`
- Wavelet: `dv/v = -1 * wct_dt/t` (sign flip applied)

**WCT-specific**: frequency averaging over `[dvv_freqmin, dvv_freqmax]` via `_freq_average_wct`, method controlled by `dvv_freq_agg`.

---

## 14. Common Pitfalls

1. **`declare_tables()` class identity**: `Job` from one call ≠ `Job` from another Python object. Always `schema = declare_tables(); Job = schema.Job`. Module-level imports work IF the module is already cached (same process), but can fail across test sessions.

2. **Session cache staleness**: After commits from other sessions or `massive_update_job`, call `db.expire_all()` OR use a fresh `connect()` for subsequent queries. `bulk_insert_mappings` bypasses ORM events entirely — lineage strings must be resolved to IDs BEFORE the call.

3. **`propagate_downstream` is idempotent**: calling twice for the same batch is safe — existing T/I jobs are left unchanged; existing D/F jobs are bumped to T (upstream changed semantics).

4. **`new_jobs --after X` in hpc=False**: logs a debug message and runs anyway (reconciliation pass). Does NOT block. This preserves backward compatibility.

5. **`_lineage_through(cat)` vs `_lineage_upstream_of(cat)`**: DVV files are saved with `upstream` lineage + step_name. Use `_lineage_upstream_of(dvv_cat)` in `get_dvv()` to match the save path. Using `_lineage_through` produces a doubled path like `.../dvv_1/dvv_1/_output/...`.

6. **Preprocess→CC is a fan-out**: one preprocess job (`pair="YA.UV05.00"`) → multiple CC jobs (one per station pair). `propagate_downstream` handles this via `create_cc_jobs_from_preprocess`, NOT the generic pair×day loop. The generic loop would create CC jobs with single-station pairs (wrong).

7. **`days` and `pairs` in `propagate_downstream`**: always deduplicate with `dict.fromkeys()` before iterating. The `batch["days"]` list has one entry per job, so with `group_by="day_lineage"` (3 pairs × 1 day = 3 entries), all three have the same day.

8. **`params.category_layer`**: gives the config layer of the *innermost* (current) step in the LayeredParams. Used in `aggregate_dvv_pairs` to access `dvv_quality_min`, `dvv_weighted_mean`, etc. without knowing the exact category name.

9. **`Config.set_number` migration**: old projects may have global config rows with `set_number=NULL`. `get_config_sets_organized()` migrates these to `set_number=1` on first call. Subsequent code always uses `set_number=1` for global config.

10. **`job.config_category` / `job.config_set_number`**: these are `@property` fields that go through the `workflow_step` relationship. The relationship must be loaded (it is in normal worker flow, since `get_next_job_for_step` does a join). Do not access them on bare Job objects constructed without a session.
