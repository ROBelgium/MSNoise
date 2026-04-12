# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # CC vs PCC — Side-by-side Comparison
#
# This notebook builds a minimal MSNoise project from scratch and runs two
# parallel correlation lineages on the same dataset:
#
# | Step name       | Method | Notes                              |
# |-----------------|--------|------------------------------------|
# | `preprocess_1`  | —      | shared preprocessing               |
# | `cc_1`          | CC     | standard geometrically-normalised  |
# | `cc_2`          | PCC    | phase cross-correlation (v=2)      |
# | `filter_1`      | —      | shared band-pass filter            |
# | `stack_1`       | linear | shared stacking config             |
#
# Both `cc_1` and `cc_2` feed the same `filter_1` → `stack_1` steps.
# The final section retrieves stacked CCFs from both lineages and plots
# them side-by-side for a direct comparison.
#
# ## Prerequisites
#
# * MSNoise installed (with the PCC patch applied).
# * A one-day MiniSEED dataset for **at least two stations**.
#   The classic MSNoise tutorial dataset (network `YA`, stations `UV05`,
#   `UV06`, `UV10`, day 2010-09-01 / Julian day 244, stored in PDF layout)
#   is used as the reference, but **any SDS or PDF archive works** — just
#   adjust the USER SETTINGS cell below.

# %% [markdown]
# ## 0 · Imports

# %%
import os
import shutil
import logging
import tempfile

import numpy as np
import matplotlib
matplotlib.use("Agg")          # headless — switch to "Qt5Agg" for interactive
import matplotlib.pyplot as plt
import xarray as xr

# MSNoise core
from msnoise.core.db import connect, create_database_inifile, get_engine
from msnoise.core.config import (
    create_config_set, update_config, get_config, get_params,
    get_config_set_details,
)
from msnoise.core.stations import update_station
from msnoise.core.workflow import (
    create_workflow_steps_from_config_sets,
    create_workflow_links_from_steps,
    create_workflow_step,
    create_workflow_link,
    get_workflow_steps,
)
from msnoise.msnoise_table_def import declare_tables, DataAvailability
from msnoise.results import MSNoiseResult

# MSNoise compute steps
from msnoise.s01_scan_archive import main as scan_archive
from msnoise.s02_new_jobs import main as new_jobs
from msnoise.s02_preprocessing import main as preprocess
from msnoise.s03_compute_no_rotation import main as compute_cc
from msnoise.s04_stack_mov import main as stack_mov

logging.basicConfig(level=logging.WARNING)   # suppress INFO noise in notebook

# %% [markdown]
# ## 1 · User Settings
#
# **Edit this cell to match your dataset.**
#
# The defaults match the classic MSNoise tutorial dataset (YA network,
# stations UV05/UV06/UV10, PDF archive layout).

# %%
# ── Path to your waveform archive ────────────────────────────────────────────
# Can be an absolute path or relative.
# The installer will store it verbatim in the DataSource table.
DATA_PATH = "/path/to/your/data"          # ← EDIT THIS

# ── Archive layout ────────────────────────────────────────────────────────────
# "PDF"  →  YEAR/STA/CHAN.D/NET.STA.LOC.CHAN.D.YEAR.JDAY
# "SDS"  →  YEAR/NET/STA/CHAN.D/NET.STA.LOC.CHAN.D.YEAR.JDAY
DATA_STRUCTURE = "PDF"                    # ← "SDS" or "PDF"

# ── Network / channel filter used when scanning the archive ──────────────────
NETWORK_CODE = "YA"
CHANNELS     = "HHZ"

# ── Date range to process (ISO format, inclusive) ────────────────────────────
STARTDATE = "2010-09-01"
ENDDATE   = "2010-09-01"

# ── Station coordinates (net, sta, lon°E, lat°N, elev_m) ─────────────────────
# Fill in your actual stations.  These are the YA tutorial stations.
STATIONS = [
    ("YA", "UV05",  29.735,  -17.817, 1174.0),
    ("YA", "UV06",  29.785,  -17.827, 1162.0),
    ("YA", "UV10",  29.790,  -17.847, 1180.0),
]

# ── CC config (shared by both cc_1 and cc_2) ─────────────────────────────────
CC_SAMPLING_RATE = 20.0     # Hz — resample target
MAXLAG           = 120.0    # seconds
CORR_DURATION    = 1800.0   # seconds per window
COMPONENTS       = "ZZ"

# ── Filter passband ───────────────────────────────────────────────────────────
FREQMIN = 0.1   # Hz
FREQMAX = 1.0   # Hz

# ── Stack config ──────────────────────────────────────────────────────────────
MOV_STACK = "(('1D','1D'),)"   # one 1-day stack

# ── Working directory (project folder) ───────────────────────────────────────
# A fresh temporary directory is created automatically.
# Set to a fixed path if you want a persistent project between kernel restarts.
WORK_DIR = None   # None → auto temp dir

# %% [markdown]
# ## 2 · Create Project Directory and MSNoise Database

# %%
if WORK_DIR is None:
    WORK_DIR = tempfile.mkdtemp(prefix="msnoise_ccvspcc_")
    print(f"Created temporary project directory: {WORK_DIR}")
else:
    os.makedirs(WORK_DIR, exist_ok=True)
    print(f"Using project directory: {WORK_DIR}")

os.chdir(WORK_DIR)

# Initialise SQLite database
create_database_inifile(
    tech=1,
    hostname=os.path.join(WORK_DIR, "msnoise.sqlite"),
    database="", username="", password="", prefix="",
)

# Create all tables
db = connect()
declare_tables().Base.metadata.create_all(db.get_bind())

print("Database created:", os.path.join(WORK_DIR, "msnoise.sqlite"))

# %% [markdown]
# ## 3 · Create Default Config Sets and Workflow Skeleton
#
# We create **one** config set for every category except `cc`, for which we
# will create **two** (cc_1 = CC, cc_2 = PCC).

# %%
ALL_CATEGORIES = [
    "global", "preprocess", "filter", "stack", "refstack",
    "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
    "stretching", "stretching_dvv",
    "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
    "psd", "psd_rms",
]

for cat in ALL_CATEGORIES:
    sn = create_config_set(db, cat)
    print(f"  {cat}_1 created (set_number={sn})")

# CC set 1 — standard cross-correlation
cc1_sn = create_config_set(db, "cc")
print(f"  cc_{cc1_sn} created  ← standard CC")

# CC set 2 — phase cross-correlation
cc2_sn = create_config_set(db, "cc")
print(f"  cc_{cc2_sn} created  ← PCC")

db.commit()

# %% [markdown]
# ## 4 · Configure Global Parameters

# %%
OUTPUT_FOLDER = os.path.join(WORK_DIR, "OUTPUT")

global_params = {
    "output_folder": OUTPUT_FOLDER,
    "startdate":     STARTDATE,
    "enddate":       ENDDATE,
    "hpc":           "N",     # propagate_downstream fires automatically
}
for k, v in global_params.items():
    update_config(db, k, v, category="global", set_number=1)
    print(f"  global.{k} = {v!r}")

# %% [markdown]
# ## 5 · Configure the Two CC Sets

# %%
# ── Shared CC parameters ──────────────────────────────────────────────────────
cc_shared = {
    "cc_sampling_rate":  str(CC_SAMPLING_RATE),
    "maxlag":            str(MAXLAG),
    "corr_duration":     str(CORR_DURATION),
    "components_to_compute": COMPONENTS,
    "keep_all":          "Y",
    "keep_days":         "Y",
    "whitening":         "A",        # whiten all inter-station pairs
    "winsorizing":       "3.0",      # 3× RMS clipping
    "clip_after_whiten": "N",
    "stack_method":      "linear",
}

for k, v in cc_shared.items():
    update_config(db, k, v, category="cc", set_number=cc1_sn)
    update_config(db, k, v, category="cc", set_number=cc2_sn)

# ── cc_1: standard CC ─────────────────────────────────────────────────────────
update_config(db, "cc_type", "CC", category="cc", set_number=cc1_sn)
print(f"cc_{cc1_sn}: cc_type = CC  (standard cross-correlation)")

# ── cc_2: PCC — amplitude normalisation implicit, clipping/whitening optional ─
# PCC2 discards amplitude per-sample so temporal normalisation is less critical.
# We keep whitening for broad spectral content but disable winsorising clipping.
update_config(db, "cc_type",    "PCC", category="cc", set_number=cc2_sn)
update_config(db, "winsorizing", "0",   category="cc", set_number=cc2_sn)  # off
print(f"cc_{cc2_sn}: cc_type = PCC (phase cross-correlation v=2)")
print(f"          winsorizing disabled (amplitude discarded by PCC)")

db.commit()

# %% [markdown]
# ## 6 · Configure Filter and Stack

# %%
update_config(db, "freqmin", str(FREQMIN), category="filter", set_number=1)
update_config(db, "freqmax", str(FREQMAX), category="filter", set_number=1)
update_config(db, "CC",      "Y",          category="filter", set_number=1)

update_config(db, "mov_stack", MOV_STACK, category="stack", set_number=1)
update_config(db, "ref_begin", STARTDATE,  category="refstack", set_number=1)
update_config(db, "ref_end",   ENDDATE,    category="refstack", set_number=1)

db.commit()
print(f"filter_1: {FREQMIN}–{FREQMAX} Hz")
print(f"stack_1:  {MOV_STACK}")

# %% [markdown]
# ## 7 · Configure DataSource and Station Table

# %%
# Create the default DataSource (the installer normally does this; when
# initialising the DB manually we must insert it ourselves).
DataSource = declare_tables().DataSource
ds = DataSource(
    name="local",
    uri=os.path.realpath(DATA_PATH),
    data_structure=DATA_STRUCTURE,
    auth_env="MSNOISE",
    network_code=NETWORK_CODE,
    channels=CHANNELS,
)
db.add(ds)
db.commit()
print(f"DataSource created: uri={os.path.realpath(DATA_PATH)!r}")
print(f"                    data_structure={DATA_STRUCTURE!r}")

# Add stations
for net, sta, lon, lat, elev in STATIONS:
    update_station(db, net=net, sta=sta, X=lon, Y=lat, altitude=elev,
                   coordinates="DEG", used=1)
    print(f"  Added station {net}.{sta}  ({lon:.3f}°E, {lat:.3f}°N)")

db.commit()

# %% [markdown]
# ## 8 · Build the Workflow Graph
#
# The topology we want:
#
# ```
# preprocess_1 ──► cc_1 ──► filter_1 ──► stack_1
#              └──► cc_2 ──►    ↑
# ```
#
# Both cc steps feed the **same** filter_1 and stack_1.

# %%
# Create WorkflowSteps from all config sets (auto-discovers preprocess_1,
# cc_1, cc_2, filter_1, stack_1, …)
created, existing, err = create_workflow_steps_from_config_sets(db)
assert err is None, f"Error creating workflow steps: {err}"
print(f"Workflow steps: {created} created, {existing} already existed")

# Auto-link following MSNoise's canonical dependency rules
created_links, existing_links, err = create_workflow_links_from_steps(db)
assert err is None, f"Error creating workflow links: {err}"
print(f"Workflow links: {created_links} created, {existing_links} already existed")

# Verify the topology
steps = {s.step_name: s for s in get_workflow_steps(db)}
print("\nWorkflow steps present:", sorted(steps.keys()))

# %% [markdown]
# ## 9 · Scan Archive and Seed Jobs

# %%
# Scan the waveform archive → populate DataAvailability table
scan_archive(init=True, threads=1)

# Update loc/chan on Station rows from DataAvailability
from sqlalchemy import text as _text
_db2 = connect()
for sta in _db2.query(declare_tables().Station):
    data = _db2.query(DataAvailability). \
        filter(_text("net=:net")).filter(_text("sta=:sta")). \
        group_by(DataAvailability.net, DataAvailability.sta,
                 DataAvailability.loc, DataAvailability.chan). \
        params(net=sta.net, sta=sta.sta).all()
    sta.used_location_codes = ",".join(sorted({d.loc for d in data}))
    sta.used_channel_names  = ",".join(sorted({d.chan for d in data}))
_db2.commit()
_db2.close()

# Verify availability
_db3 = connect()
n_da = _db3.query(DataAvailability).count()
print(f"DataAvailability rows: {n_da}")
_db3.close()

# Seed initial jobs (creates preprocess_1 T-jobs)
new_jobs(init=True)

# %% [markdown]
# ## 10 · Run the Pipeline
#
# Steps run sequentially in-process.  For large datasets use `msnoise` CLI
# or HPC submission instead.

# %%
print("── preprocess ─────────────────────────────────────────────────────")
preprocess()

print("── new_jobs --after preprocess ────────────────────────────────────")
new_jobs(after="preprocess")

print("── compute_cc (cc_1 = CC  and  cc_2 = PCC) ───────────────────────")
compute_cc()

print("── new_jobs --after cc ────────────────────────────────────────────")
new_jobs(after="cc")

print("── stack_mov ──────────────────────────────────────────────────────")
stack_mov(stype="mov")

print("── new_jobs --after stack ─────────────────────────────────────────")
new_jobs(after="stack")

print("Done.")

# %% [markdown]
# ## 11 · Gather MSNoiseResult Objects

# %%
db = connect()

# Lineage for CC branch:  preprocess_1 / cc_1 / filter_1 / stack_1
result_cc = MSNoiseResult.from_ids(
    db,
    preprocess=1,
    cc=cc1_sn,        # cc_1
    filter=1,
    stack=1,
)
print("CC  lineage:", result_cc.lineage_names)

# Lineage for PCC branch: preprocess_1 / cc_2 / filter_1 / stack_1
result_pcc = MSNoiseResult.from_ids(
    db,
    preprocess=1,
    cc=cc2_sn,        # cc_2
    filter=1,
    stack=1,
)
print("PCC lineage:", result_pcc.lineage_names)

# %% [markdown]
# ## 12 · Compare Stacked CCFs

# %%
# Retrieve all stacked CCFs for both lineages
ccfs_cc  = result_cc.get_ccf()
ccfs_pcc = result_pcc.get_ccf()

print(f"CC  result: {len(ccfs_cc)}  CCF(s) found")
print(f"PCC result: {len(ccfs_pcc)} CCF(s) found")

for key in sorted(ccfs_cc):
    pair, comp, ms = key
    print(f"  {pair}  {comp}  {ms[0]}–{ms[1]}")

# %% [markdown]
# ## 13 · Plot: CC vs PCC per Pair

# %%
fig_dir = os.path.join(WORK_DIR, "figures")
os.makedirs(fig_dir, exist_ok=True)

common_keys = sorted(set(ccfs_cc) & set(ccfs_pcc))
if not common_keys:
    print("No common keys found between CC and PCC results — check that "
          "the pipeline completed successfully.")
else:
    n = len(common_keys)
    fig, axes = plt.subplots(n, 2, figsize=(14, 3.5 * n), squeeze=False)
    fig.suptitle(
        f"Stacked CCF comparison: CC (left) vs PCC (right)\n"
        f"Filter {FREQMIN}–{FREQMAX} Hz  |  {STARTDATE}",
        fontsize=13, y=1.01,
    )

    for row, key in enumerate(common_keys):
        pair, comp, ms = key
        da_cc  = ccfs_cc[key]
        da_pcc = ccfs_pcc[key]

        # Time axis (seconds)
        taxis = da_cc.coords["taxis"].values if "taxis" in da_cc.coords \
            else np.linspace(-MAXLAG, MAXLAG, da_cc.shape[-1])

        # Stack over the time dimension if multiple days present
        ccf_cc  = float(da_cc.mean())  if da_cc.ndim == 0 else da_cc.values
        ccf_pcc = float(da_pcc.mean()) if da_pcc.ndim == 0 else da_pcc.values

        if ccf_cc.ndim > 1:
            ccf_cc  = ccf_cc.mean(axis=0)
        if ccf_pcc.ndim > 1:
            ccf_pcc = ccf_pcc.mean(axis=0)

        label = f"{pair}  {comp}"

        # ── CC panel ──────────────────────────────────────────────────────
        ax = axes[row, 0]
        ax.plot(taxis, ccf_cc, color="#1f77b4", lw=0.9)
        ax.axvline(0, color="k", lw=0.6, ls="--", alpha=0.4)
        ax.set_title(f"{label}  —  CC")
        ax.set_xlabel("Lag (s)")
        ax.set_ylabel("Amplitude")
        ax.set_xlim(-MAXLAG, MAXLAG)

        # ── PCC panel ─────────────────────────────────────────────────────
        ax = axes[row, 1]
        ax.plot(taxis, ccf_pcc, color="#d62728", lw=0.9)
        ax.axvline(0, color="k", lw=0.6, ls="--", alpha=0.4)
        ax.set_title(f"{label}  —  PCC")
        ax.set_xlabel("Lag (s)")
        ax.set_ylabel("Amplitude")
        ax.set_xlim(-MAXLAG, MAXLAG)

    plt.tight_layout()
    outfig = os.path.join(fig_dir, "cc_vs_pcc.png")
    fig.savefig(outfig, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {outfig}")

# %% [markdown]
# ## 14 · Overlay: CC vs PCC on the same axes (one panel per pair)

# %%
if common_keys:
    fig2, axes2 = plt.subplots(1, len(common_keys),
                               figsize=(6 * len(common_keys), 4),
                               squeeze=False)
    fig2.suptitle(
        f"CC (blue) vs PCC (red)  |  {FREQMIN}–{FREQMAX} Hz  |  {STARTDATE}",
        fontsize=12,
    )

    for col, key in enumerate(common_keys):
        pair, comp, ms = key
        da_cc  = ccfs_cc[key]
        da_pcc = ccfs_pcc[key]

        taxis = da_cc.coords["taxis"].values if "taxis" in da_cc.coords \
            else np.linspace(-MAXLAG, MAXLAG, da_cc.values.shape[-1])

        ccf_cc  = da_cc.values
        ccf_pcc = da_pcc.values
        if ccf_cc.ndim > 1:  ccf_cc  = ccf_cc.mean(axis=0)
        if ccf_pcc.ndim > 1: ccf_pcc = ccf_pcc.mean(axis=0)

        # Normalise to ABSMAX for visual overlay
        def _norm(x):
            m = np.abs(x).max()
            return x / m if m > 0 else x

        ax = axes2[0, col]
        ax.plot(taxis, _norm(ccf_cc),  color="#1f77b4", lw=1.0,
                label="CC",  alpha=0.85)
        ax.plot(taxis, _norm(ccf_pcc), color="#d62728", lw=1.0,
                label="PCC", alpha=0.85)
        ax.axvline(0, color="k", lw=0.6, ls="--", alpha=0.3)
        ax.set_title(f"{pair}  {comp}")
        ax.set_xlabel("Lag (s)")
        ax.set_ylabel("Normalised amplitude")
        ax.set_xlim(-MAXLAG, MAXLAG)
        ax.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    outfig2 = os.path.join(fig_dir, "cc_vs_pcc_overlay.png")
    fig2.savefig(outfig2, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {outfig2}")

# %% [markdown]
# ## 15 · Summary Table

# %%
import pandas as pd

rows = []
for key in common_keys:
    pair, comp, ms = key
    cc_arr  = ccfs_cc[key].values
    pcc_arr = ccfs_pcc[key].values
    if cc_arr.ndim  > 1: cc_arr  = cc_arr.mean(axis=0)
    if pcc_arr.ndim > 1: pcc_arr = pcc_arr.mean(axis=0)
    rows.append({
        "Pair":        pair,
        "Component":   comp,
        "Mov. stack":  f"{ms[0]}–{ms[1]}",
        "CC peak":     f"{np.abs(cc_arr).max():.4f}",
        "PCC peak":    f"{np.abs(pcc_arr).max():.4f}",
        "CC  lag@peak (s)":  f"{np.where(np.abs(cc_arr)  == np.abs(cc_arr).max())[0][0] - len(cc_arr)//2:.1f}",
        "PCC lag@peak (s)":  f"{np.where(np.abs(pcc_arr) == np.abs(pcc_arr).max())[0][0] - len(pcc_arr)//2:.1f}",
    })

df = pd.DataFrame(rows)
print(df.to_string(index=False))

# %% [markdown]
# ## Notes
#
# ### Why PCC output looks different from CC
#
# * **CC** peak amplitude depends on the product of the two trace energies.
#   It is sensitive to amplitude transients (earthquakes, instrumental
#   glitches) even after one-bit normalisation or spectral whitening.
#
# * **PCC** (v=2) projects every sample of the analytic signal onto the unit
#   circle before correlating, so amplitude information is discarded entirely.
#   The output is bounded to [−1, 1] by construction.  Temporal normalisation
#   (winsorising) is therefore largely redundant and was disabled for `cc_2`.
#
# ### Adjusting the comparison
#
# * Change `FREQMIN` / `FREQMAX` in the User Settings cell and re-run from
#   §6 onwards (filter_1 is the same config set for both branches).
# * To compare **PCC v=1** instead of v=2, add
#   `update_config(db, "cc_type", "PCC1", category="cc", set_number=cc2_sn)`
#   — note that PCC1 uses the slower time-domain path.
# * The `result_cc` and `result_pcc` objects expose the full MSNoiseResult
#   API: call `.get_mwcs()`, `.get_dvv()`, etc. after running the downstream
#   steps on each lineage.
