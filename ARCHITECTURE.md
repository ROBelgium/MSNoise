# MSNoise Architecture Reference

*Self-reference for future Claude sessions. Read this before touching workflow, jobs, or lineage.*

---

## 0. Claude Session Protocol

**Always do this at the start of a session:**

1. `git pull` to get current master
2. Read this file
3. Run `msnoise utils test --fast` from a test project dir to confirm baseline (see §15)
4. Never modify the DB schema or `declare_tables()` without understanding §10
5. Patch workflow: Claude generates `git diff`-style patch → Thomas applies locally, tests, pushes → Claude pulls before next patch

**Key rule**: Always `git pull` before generating a patch. Never stack patches without pulling in between.

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
    __init__.py          # Re-exports all 6 submodules via *; signal=DSP/PSD, io=file I/O
    db.py                # connect, get_engine, read_db_inifile, get_logger
    config.py            # get_config, update_config, create_config_set, get_params,
                         # get_merged_params_for_lineage, build_plot_outfile,
                         # lineage_to_plot_tag, get_config_categories_definition
    stations.py          # get_stations, get_station_pairs, update_data_availability,
                         # get_data_availability, add_data_source, get_data_source,
                         # resolve_data_source, get_waveform_path,
                         # read_waveforms_from_availability,   ← DA records → Stream
                         # import_stationxml,
                         # set_station_source, set_network_source, set_all_stations_source
    workflow.py          # get_next_lineage_batch, propagate_downstream,
                         # is_next_job_for_step, massive_insert_job, massive_update_job,
                         # update_job, reset_jobs, build_ref_datelist,
                         # get_lineages_to_step_id, get_done_lineages_for_category,
                         # _get_or_create_lineage_id, _lineage_id_for, _lineage_str_from_id,
                         # create_workflow_steps_from_config_sets,
                         # create_workflow_links_from_steps,
                         # get_workflow_steps, get_workflow_links, get_workflow_job_counts,
                         # get_t_axis, build_movstack_datelist,
                         # get_filter_steps_for_cc_step, get_refstack_lineage_for_filter,
                         # refstack_is_rolling, refstack_needs_recompute, extend_days,
                         # compute_rolling_ref, lineage_str_to_step_names, lineage_str_to_steps
    signal.py            # winsorizing, stack, xwt, compute_wct_dtt, get_wct_avgcoh,
                         # preload_instrument_responses, save_preprocessed_streams,
                         # get_preprocessed_stream, validate_stack_data,
                         # psd_rms, psd_df_rms,       ← pure DSP; psd_df_rms takes xarray Dataset, returns Dataset
                         # make_same_length            ← kept as public API for plugins
    io.py                # xr_save_ccf, xr_get_ccf, xr_save_ccf_daily, xr_get_ccf_daily,
                         # xr_save_ccf_all, xr_get_ccf_all,
                         # xr_save_ref, xr_get_ref, xr_load_ccf_for_stack,
                         # xr_save_mwcs, xr_get_mwcs, xr_save_dtt, xr_get_dtt,
                         # xr_save_stretching, xr_save_wct, xr_load_wct,
                         # xr_save_wct_dtt, xr_get_wct_dtt,
                         # xr_save_dvv_agg, xr_get_dvv_agg, aggregate_dvv_pairs,
                         # xr_save_psd, xr_load_psd, xr_save_rms, xr_load_rms,
                         # psd_ppsd_to_dataframe, psd_ppsd_to_dataset, psd_read_results
    fdsn.py              # parse_datasource_scheme, is_remote_source, build_client,
                         # fetch_waveforms_bulk,
                         # fetch_raw_waveforms,        ← raw only, no preprocessing (for PSD)
                         # fetch_and_preprocess,       ← fetch + full preprocessing pipeline
                         # _write_raw_cache
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
  s02_new_jobs.py        # new_jobs main; all propagate_* functions (see §6)
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
  psd_compute.py         # PSD computation main — DataSource-aware (local SDS or FDSN)
  psd_compute_rms.py     # PSD RMS computation main
  config/                # One config_<category>.csv per workflow category
  plots/                 # ccftime, data_availability, distance, dtt, interferogram,
                         # mwcs, mwcs_dtt_dvv, ppsd, psd_rms, spectime,
                         # station_map, stretching_dvv, timing, wavelet_dtt_dvv
  scripts/msnoise.py     # Full Click CLI; entry point msnoise = scripts.msnoise:run
  msnoise_admin.py       # Flask-Admin web UI (msnoise admin)
  test/
    test_smoke.py        # Fast lifecycle smoke tests (31 tests, ~20s) — run with --fast
    tests.py             # Full integration tests (101 tests, ~2 min) — uses real data
```

**The one universal entry point**: `from msnoise.api import connect` (or `from msnoise import connect`).
All other symbols live in `msnoise.core.*` and are re-exported through `api.py` for backward compat.

---

## 3. Jobs — Two Distinct Concepts

### `Job` (ORM object, `msnoise_table_def.declare_tables().Job`)
SQLAlchemy mapped class. **Always call `declare_tables()` to get a fresh class** — do NOT use module-level imports like `from msnoise_table_def import Job` for DB queries, as multiple `declare_tables()` calls create different class objects. Use `schema = declare_tables(); Job = schema.Job`.

Key columns: `ref` (PK), `day`, `pair`, `flag` (T/I/D/F), `step_id` (FK→WorkflowStep), `lineage_id` (FK→Lineage), `jobtype` (= step_name, used as join key to WorkflowStep), `priority`, `lastmod`.

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

Every job carries a `lineage_id` FK → `Lineage(lineage_str)`. The string is a `/`-separated path of step names from root to the current step (inclusive):

```
preprocess_1
preprocess_1/cc_1
preprocess_1/cc_1/filter_1/stack_1      ← filter_1 is pass-through but appears in path!
preprocess_1/cc_1/filter_1/stack_1/refstack_1
preprocess_1/cc_1/filter_1/stack_1/refstack_1/mwcs_1/mwcs_dtt_1/mwcs_dtt_dvv_1
psd_1
psd_1/psd_rms_1                         ← psd_rms must include its own step name
```

**Deduplication**: multiple jobs share the same `Lineage` row. The `before_insert` event on `Job` resolves `lineage_str` → `lineage_id` via raw SQL (not ORM session) to avoid re-entrant flush.

**Key helpers in `core/workflow.py`**:
- `_get_or_create_lineage_id(session, str)` — 4-step lookup: identity_map → session.new → DB → INSERT
- `_lineage_id_for(session, str)` — read-only, returns None if not found
- `_lineage_str_from_id(session, id)` — reverse lookup, safe after session.commit() expiry
- `get_lineages_to_step_id(session, step_id)` — DFS returning all upstream paths to a step

**Origin steps** (no upstream): `preprocess_N` and `psd_N` use `lineage = step.step_name` (e.g. `"psd_1"`). All downstream steps use the full path. `psd_rms_N` jobs must use `"psd_1/psd_rms_1"` — if only `"psd_1"` is stored, `get_next_lineage_batch` builds `LayeredParams(['global','psd'])` and `params.psd_rms.*` raises `AttributeError`.

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

### Full dispatch table

| Completed category | `s02_new_jobs.py` function | Notes |
|---|---|---|
| `preprocess` | `propagate_stack_jobs_from_cc_done` | Fan-out: 1 station → N pairs |
| `cc` (and pass-throughs) | `propagate_first_runnable_from_category` | Crosses filter_1 transparently via `_collect()` |
| `stack` | `propagate_refstack_jobs_from_stack_done` + `propagate_mwcs_jobs_from_refstack_done` | Dual: REF sentinel + direct mwcs/str/wavelet jobs |
| `refstack` | `propagate_mwcs_jobs_from_refstack_done` | REF just changed: re-propagate all days |
| `mwcs_dtt`/`stretching`/`wavelet_dtt` | `propagate_dvv_jobs_from_dtt_done` | Creates `day="DVV", pair="ALL"` sentinel |
| `psd` | `propagate_psd_rms_jobs_from_psd_done` | Single-station passthrough |

### HPC mode (`hpc=Y`)

`propagate_downstream` is NOT called by workers. Operator runs `msnoise new_jobs --after X` manually. In hpc=False mode, `new_jobs --after X` still works as a reconciliation/safety pass (logs debug message + runs anyway).

---

## 7. The REF Sentinel — When It Computes vs Skips

`s04_stack_refstack` checks on every REF job using `refstack_needs_recompute()` from `core/workflow.py`:

```python
if not refstack_is_rolling(params):
    if not refstack_needs_recompute(db, pair, batch["lineage_names_upstream"], params):
        # No Done stack days fall inside ref_begin..ref_end window.
        # Reference unchanged — mark Done without recomputation.
        mark_done(); propagate_downstream(); continue

# New data inside window (or rolling mode) — recompute reference
compute_ref_stack(); write_ref_file(); mark_done(); propagate_downstream()
```

`refstack_needs_recompute(session, pair, lineage_names_upstream, params)` lives in `core/workflow.py`. It queries Done stack jobs for the pair, resolves their lineage, and returns `True` if any date falls inside `[ref_begin, ref_end]`. It uses `declare_tables()` internally (not the module-level `WorkflowStep`) so it is session-safe.

**Mode B (rolling ref)**: `ref_begin` is a negative integer (e.g. `"-5"` = last 5 days). No REF file written; reference computed on-the-fly in MWCS/stretching/WCT. The Mode B path marks Done and calls `propagate_downstream` immediately (in `hpc=False`) — downstream MWCS/stretching/wavelet jobs are created inline, not via `new_jobs --after refstack`.

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

# PSD:
OUTPUT/psd_1/_output/daily/YA.UV05.00.HHZ/2010-09-01.nc
OUTPUT/psd_1/psd_rms_1/_output/YA.UV05.00.HHZ.nc

# DVV aggregates:
OUTPUT / upstream_lineage / dvv_step_name / _output / 1D_1D / dvv_CC_ZZ.nc
```

**Critical**: `get_dvv()` uses `_lineage_upstream_of(dvv_cat)` (not `_lineage_through`) to match the save path. PSD uses `lineage=[]` because it is a root step with no upstream folder.

---

## 9. Config & Params

Config lives in the unified `Config` table. One row per `(name, category, set_number)`. Global config has `category='global'`, `set_number=1`.

Key columns on `Config`: `name`, `category`, `set_number`, `value` (always string), `param_type` ('str'/'int'/'float'/'bool'), `default_value`, `description`, `possible_values` (slash-separated).

- `get_params(db)` → single-layer `LayeredParams`. Access `params.global_.hpc`, `params.global_.output_folder`, etc.
- `get_merged_params_for_lineage(db, orig_params, step_params, lineage)` → `(lineage, lineage_names, LayeredParams)`. Access: `params.global_` (underscore avoids Python keyword), `params.cc`, `params.mwcs` etc.
- `get_config_set_details(db, category, set_number, format="AttribDict")` → plain `AttribDict` for a single config set.

**`LayeredParams` access**:
- `params.global_.hpc` — global layer
- `params["stack"].mov_stack` — dict-style also works
- `params.category_layer` — the *current step's* layer (useful when the step name isn't known statically)
- `params.step_name` — name of the innermost step
- `params.lineage_names` — full list of step names in this lineage
- `params.categories` — list of all loaded category names
- `params.as_flat_dict()` — all params flattened (later layers win on name collision)
- `params.to_yaml(path)` / `LayeredParams.from_yaml(path)` — serialization

**NOT** `params.hpc` — that's `params.global_.hpc`.

Config CSV files in `msnoise/config/config_<category>.csv` define defaults. Columns: `name`, `default`, `definition`, `type`, `possible_values`.

**CLI config helpers** (`core/config.py`):
- `parse_config_key(key)` → `(category, set_number, name)`. Accepts `name` (global), `category.name` (set 1), `category.N.name` (explicit set).
- `_cast_config_value(name, value_str, param_type)` → validated string. Validates against `param_type` (`str/int/float/bool/eval`); normalises bools to `Y`/`N`. Raises `ValueError` with a clear message on bad input. Always call before `update_config()` in CLI code.

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

**`core/fdsn.py` key functions**:
- `is_remote_source(uri)` → True if FDSN/EIDA scheme
- `build_client(ds)` → ObsPy `Client` or EIDA routing client from a `DataSource`
- `fetch_waveforms_bulk(client, bulk_request, retries=3)` → `Stream` (low-level, retry logic)
- `fetch_raw_waveforms(db, jobs, goal_day, params, t_start, t_end)` → raw `Stream` with no preprocessing. Used by `psd_compute.py` because PPSD needs raw waveforms + handles its own response correction. Accepts optional time window overrides (PSD needs padding before midnight).
- `fetch_and_preprocess(db, jobs, goal_day, params, responses, loglevel)` → `(stream, done_jobs, failed_jobs)` — delegates fetch to `fetch_raw_waveforms`, then applies the full preprocessing pipeline
- `_write_raw_cache(stream, ...)` — caches raw waveforms when `fdsn_keep_raw=Y` (global config)

**`core/stations.py` key functions**:
- `resolve_data_source(db, station)` — returns effective DS (station override or default)
- `get_waveform_path(db, da)` — joins `DataSource.uri + da.path + da.file` (never use `os.path.join(da.path, da.file)` directly — it ignores the DataSource root)
- `read_waveforms_from_availability(db, da_records, t_start, t_end, logger)` → `Stream`. The canonical "DA records → waveforms" helper. Uses `get_waveform_path` internally.

**Station-to-source assignment**:
- `set_station_source(db, net, sta, data_source_id)` — per-station override
- `set_network_source(db, net, data_source_id)` — all stations in a network
- `set_all_stations_source(db, data_source_id)` — project-wide

**PSD DataSource pattern** (mirrors preprocess step):
```python
ds = resolve_data_source(db, get_station(db, net, sta))
if is_remote_source(ds.uri):
    st = fetch_raw_waveforms(db, [job], goal_day, params, t_start=..., t_end=...)
else:
    da_records = get_data_availability(db, net=net, sta=sta, ...)
    st = read_waveforms_from_availability(db, da_records, t_start, t_end)
```

---

## 12. MSNoiseResult — xarray-Only Result Reader (v2.x)

`MSNoiseResult` (`results.py`) is the user-facing class for reading computed results. **All methods return xarray Dataset or DataArray objects** (no `format=` parameter). For pandas conversion, use the static helper `MSNoiseResult.to_dataframe(ds)`.

```python
# Construct from integer set IDs
r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                            stack=1, refstack=1, mwcs=1, mwcs_dtt=1,
                            mwcs_dtt_dvv=1)

# List all available result sets for a category
for r in MSNoiseResult.list(db, "mwcs_dtt"):
    print(r)

# Access results — all return xarray Dataset or DataArray
da = r.get_ccf(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_ref(pair="BE.UCC--BE.MEM", components="ZZ")
ds = r.get_mwcs(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_mwcs_dtt(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_stretching(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_wct(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_wct_dtt(pair="BE.UCC--BE.MEM", components="ZZ", mov_stack=("1D","1D"))
ds = r.get_psd(seed_id="BE.UCC..HHZ", day="2023-01-01")
ds = r.get_psd_rms(seed_id="BE.UCC..HHZ")

# Convert to pandas when needed (escape hatch for all steps)
df = MSNoiseResult.to_dataframe(ds)

# Export DVV with full provenance (NetCDF + YAML metadata)
written = r.export_dvv("output/")   # returns list of file paths

# Navigate branches
for branch in r.branches():   # other lineages reachable from same root
    print(branch)
```

**xarray-only API** (since v2.x): The `format=` parameter has been removed from all workflow step methods (CC, stack, ref, MWCS, stretching, WCT, DVV aggregates). All `get_*` methods return xarray Dataset or DataArray objects. For pandas users:

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

**WCT-DTT variable names**: `DTT`, `ERR`, `COH` (uppercase) — written by `s09_compute_wct_dtt.py` and read by `_freq_average_wct` in `core/io.py`. Do not use lowercase `dtt/err/coh`.

**Sign convention**:
- MWCS: `dv/v = -1 * dt/t` (sign flip applied in aggregator)
- Stretching: `dv/v = stretching_factor - 1`
- Wavelet: `dv/v = -1 * wct_dt/t` (sign flip applied)

---

## 14. CLI — Key Commands

Entry point: `msnoise` (maps to `scripts/msnoise.py:run`).

**Job management**:
```sh
msnoise new_jobs                         # scan DA + create preprocess/cc/stack jobs
msnoise new_jobs --after psd             # propagate psd→psd_rms
msnoise reset cc_1 --all                 # reset all cc_1 jobs to T
msnoise utils create_preprocess_jobs --date 2026-03-28        # FDSN bypass
msnoise utils create_preprocess_jobs --date_range START END   # FDSN bypass, range
msnoise utils create_psd_jobs --date 2026-03-28               # PSD jobs for FDSN station
msnoise utils create_psd_jobs --date_range START END --set-number 1
```

`create_preprocess_jobs` and `create_psd_jobs` create T jobs directly without scanning DataAvailability — essential when adding a new FDSN-sourced station. For local/SDS sources they only create jobs where DA records already exist.

**Compute workers**:
```sh
msnoise cc compute                       # runs s03_compute_no_rotation
msnoise cc compute --chunk-size 50       # claim 50 pairs/day (large networks, ≥50 stations)
msnoise -t 8 cc compute --chunk-size 50  # 8 parallel workers × 50 pairs each
msnoise qc compute_psd --chunk-size 20   # claim 20 stations/day for PSD
msnoise cc stack -m                      # moving stack
msnoise cc stack_refstack                # reference stack
msnoise cc dtt compute_mwcs             # MWCS
msnoise cc dtt compute_mwcs_dtt         # MWCS DTT
msnoise cc dtt dvv compute_mwcs_dtt_dvv # DVV aggregate
msnoise qc compute_psd                  # PSD
msnoise qc compute_psd_rms              # PSD RMS
msnoise -t 4 qc compute_psd             # parallel (4 workers)
```

**Configuration** (dot-notation, see §9):
```sh
# Read
msnoise config get output_folder             # global shorthand
msnoise config get cc.cc_sampling_rate       # category.name → set 1
msnoise config get mwcs.2.mwcs_wlen          # explicit set number

# Write (validates type before writing)
msnoise config set output_folder /data/output
msnoise config set cc.cc_sampling_rate 25
msnoise config set mwcs.2.mwcs_wlen 10

# Reset to default
msnoise config reset cc.cc_sampling_rate
msnoise config reset mwcs.2.mwcs_wlen

# List / inspect
msnoise config list                          # all categories, all sets
msnoise config list cc                       # all cc sets
msnoise config list mwcs.2                   # mwcs set 2 only
                                             # * marks non-default values

# Config set management
msnoise config create_set mwcs               # add a second mwcs set
msnoise config delete_set mwcs 2             # remove mwcs set 2
msnoise config list_sets                     # show all category × set_number combos
msnoise config show_set mwcs 1               # detailed view of one set
msnoise config copy_set mwcs 1 mwcs 2        # copy set 1 → set 2
```

**DB utilities**:
```sh
msnoise db upgrade          # add new config params to existing DB
msnoise db clean_duplicates # remove duplicate job rows
msnoise db dump             # export DB to CSV
```

---

## 15. Testing

```sh
# From a project directory (with msnoise.ini):
msnoise utils test --fast    # 32 smoke tests, ~20s — run after every patch
msnoise utils test           # full integration tests, ~2 min — requires test data

# Direct pytest:
cd /path/to/msnoise_project
python -m pytest /path/to/msnoise/msnoise/test/test_smoke.py -v
python -m pytest /path/to/msnoise/msnoise/test/tests.py -v

# Specific test:
python -m pytest /path/to/msnoise/msnoise/test/test_smoke.py::test_smoke_172_psd_rms_params_layer -v
```

**Smoke test coverage** (`test_smoke.py`, `@pytest.mark.order(N)`):
- 01-03: schema, workflow creation, direct job seeding (preprocess + psd)
- 04-17: full pipeline stubs (preprocess→cc→stack→refstack→mwcs/str/wct→dvv→psd→psd_rms)
- 171: psd_rms job lineage must be `"psd_1/psd_rms_1"` (not just `"psd_1"`)
- 172: `params.psd_rms.*` accessible via `get_next_lineage_batch` (LayeredParams layer test)
- 173: `get_waveform_path` prepends DataSource.uri correctly
- 174: `create_psd_jobs` CLI creates T jobs
- 175: `get_next_lineage_batch(chunk_size=2)` claims ≤2 CC jobs per day, leaves rest as T
- 18-19: lineage normalisation, DVV discoverability
- 185: DataSource station resolution
- 20-25: second-day processing, filter pass-through, idempotency

---

## 16. Common Pitfalls

1. **`declare_tables()` class identity**: `Job` from one call ≠ `Job` from another Python object. Always `schema = declare_tables(); Job = schema.Job`. Module-level imports work IF the module is already cached (same process), but can fail across test sessions.

2. **Session cache staleness**: After commits from other sessions or `massive_update_job`, call `db.expire_all()` OR use a fresh `connect()` for subsequent queries. `bulk_insert_mappings` bypasses ORM events entirely — lineage strings must be resolved to IDs BEFORE the call.

3. **`propagate_downstream` is idempotent**: calling twice for the same batch is safe — existing T/I jobs are left unchanged; existing D/F jobs are bumped to T (upstream changed semantics).

4. **`new_jobs --after X` in hpc=False**: logs a debug message and runs anyway (reconciliation pass). Does NOT block. This preserves backward compatibility.

5. **`_lineage_through(cat)` vs `_lineage_upstream_of(cat)`**: DVV files are saved with `upstream` lineage + step_name. Use `_lineage_upstream_of(dvv_cat)` in `get_dvv()` to match the save path. Using `_lineage_through` produces a doubled path like `.../dvv_1/dvv_1/_output/...`.

6. **Preprocess→CC is a fan-out**: one preprocess job (`pair="YA.UV05.00"`) → multiple CC jobs (one per station pair). `propagate_downstream` handles this via `create_cc_jobs_from_preprocess`, NOT the generic pair×day loop. The generic loop would create CC jobs with single-station pairs (wrong).

7. **`days` and `pairs` in `propagate_downstream`**: always deduplicate with `dict.fromkeys()` before iterating. The `batch["days"]` list has one entry per job, so with `group_by="day_lineage"` (3 pairs × 1 day = 3 entries), all three have the same day.

8. **`params.category_layer`**: gives the config layer of the *innermost* (current) step in the LayeredParams. Useful in `aggregate_dvv_pairs` to access `dvv_quality_min`, `dvv_weighted_mean`, etc. without knowing the exact category name. Do NOT use it as a substitute for a missing layer (e.g. if `params.psd_rms` raises AttributeError, the lineage is wrong — fix the lineage, not the accessor).

9. **`job.config_category` / `job.config_set_number`**: these are `@property` fields that go through the `workflow_step` relationship. The relationship must be loaded (it is in normal worker flow, since `get_next_job_for_step` does a join). Do not access them on bare Job objects constructed without a session.

10. **DataAvailability paths**: never do `os.path.join(da.path, da.file)` directly — this ignores the DataSource root. Always use `get_waveform_path(db, da)` or `read_waveforms_from_availability(db, da_records, ...)`.

11. **`make_same_length` is public API**: it has a deprecation warning in its docstring but is kept in `core/signal.py` because external plugins depend on it. Do not remove it.

12. **WCT-DTT variable names are uppercase**: `DTT`, `ERR`, `COH` (not `dtt`, `err`, `coh`). Test stubs must use uppercase or `xr_get_wct_dtt` will fail silently.

13. **Lazy dataset handle lifetime**: all `get_*` readers in `core/io.py` use `xr.open_dataset()` (lazy). The file handle stays open until `.close()` is called or the object is garbage-collected. **Two specific readers materialise + close immediately** because their file is also a write target in the same pipeline:
   - `xr_get_ccf` — the stacked CCF file is merged and rewritten by `xr_save_ccf` (stack worker)
   - `xr_get_ref` — the REF file is rewritten by `xr_save_ref` (refstack worker)
   All other lazy readers (MWCS, DTT, STR, WCT) read from files that are **never written back to by the same caller** — those handles are safe to keep lazy, but callers must `ds.close()` explicitly once done to avoid leaking handles (especially in loops over many pairs). Pattern: extract all needed `.values`, then call `ds.close()`.

14. **`xr_get_ccf` returns an in-memory DataArray**: despite using `open_dataset` internally, `xr_get_ccf` calls `.load()` + `.close()` before returning — the result is fully in-memory (not lazy-backed). Do not assume it needs `.load()` again.

15. **`chunk_size` is a CLI-only parameter** (not a config CSV key): pass `--chunk-size N` to `msnoise cc compute` or `msnoise qc compute_psd`. It controls how many pairs (CC) or stations (PSD) a single worker claims per day. Default 0 = claim all (original behaviour). Only effective for `group_by="day_lineage"` steps. Do NOT add it to downstream steps (stack, MWCS, stretching) — those use `pair_lineage` and write to per-pair accumulated files where concurrent writes would corrupt data.
