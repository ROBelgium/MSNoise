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
# # PSD-RMS Visualisation — Multi-Station Analysis
#
# This notebook shows how to load RMS results for **multiple stations** and
# compare them using overlaid time-series, side-by-side clock plots, and
# combined grid-maps.
#
# It also demonstrates:
# - Loading RMS data directly from NetCDF (no DB required for post-processing)
# - Before/after comparison via `split_date`
# - Annotating figures with events (e.g. maintenance windows, earthquakes)
# - Batch export of all figures
#
# **Pre-requisites**: `msnoise qc compute_psd` and `msnoise qc compute_psd_rms`
# must have completed for all stations listed below.

# %% [markdown]
# ## 1. Imports and configuration

# %%
import os
import matplotlib.pyplot as plt
import pandas as pd

from msnoise.core.db import connect
from msnoise.results import MSNoiseResult
from msnoise.plots.psd_rms import (
    DayNightWindow,
    load_rms,
    plot_timeseries,
    plot_clockplot,
    plot_hourmap,
    plot_gridmap,
    plot_dailyplot,
    localize_tz_and_reindex,
    stack_wday_time,
    radial_hours,
    clock24_plot_commons,
)

# ── Project settings — edit these ────────────────────────────────────────────
SEED_IDS = [
    "BE.UCC..HHZ",
    "BE.MEM..HHZ",
    "BE.RCHB..HHZ",
]
PSD_ID     = 1   # config-set number of psd step     (psd_1)
PSD_RMS_ID = 1   # config-set number of psd_rms step (psd_rms_1)
BAND         = None          # None = first available band; or e.g. "1.0-10.0"
SCALE        = 1e9           # to nanometres (DISP output)
UNIT         = "nm"
TIMEZONE     = "Europe/Brussels"
RESAMPLE     = "30min"

# Optional: mark events on figures
ANNOTATIONS  = {
    # "2024-01-15": "Maintenance",
    # "2024-03-01": "Storm event",
}

# Optional: split date for before/after comparison
SPLIT_DATE   = None   # e.g. "2024-02-01"

SAVE_DIR     = "./psd_rms_figures"
os.makedirs(SAVE_DIR, exist_ok=True)

# Day / Night window -- controls shading and daytime median overlay.
# Adjust to your station's local anthropogenic noise schedule.
DAY_NIGHT    = DayNightWindow(day_start=7, day_end=19)
# ─────────────────────────────────────────────────────────────────────────────

# %% [markdown]
# ## 2. Load RMS data

# %%
db = connect()
# MSNoiseResult resolves the correct psd_N/psd_rms_N/_output/ path automatically.
result = MSNoiseResult.from_ids(db, psd=PSD_ID, psd_rms=PSD_RMS_ID)
data = load_rms(result, SEED_IDS)
db.close()

available = list(data.keys())
print(f"Loaded data for: {available}")

if not data:
    raise RuntimeError(
        "No RMS data found. Run `msnoise qc compute_psd_rms` first."
    )

# Determine band to use
first_df = next(iter(data.values()))
BAND = BAND or first_df.columns[0]
print(f"Using band: {BAND}")
print(f"Time range: {first_df.index[0]}  →  {first_df.index[-1]}")

# %% [markdown]
# ## 3. Overview: RMS statistics per station

# %%
rows = []
for sid, df in data.items():
    col = df[BAND] if BAND in df.columns else df.iloc[:, 0]
    rows.append({
        "SEED ID": sid,
        "Band": BAND,
        "N samples": len(col),
        f"Median ({UNIT})": f"{col.median() * SCALE:.3f}",
        f"P5 ({UNIT})": f"{col.quantile(0.05) * SCALE:.3f}",
        f"P95 ({UNIT})": f"{col.quantile(0.95) * SCALE:.3f}",
    })

pd.DataFrame(rows).set_index("SEED ID")

# %% [markdown]
# ## 4. Multi-station time-series overlay
#
# All stations on the same axes — useful for spotting correlated noise sources
# (anthropogenic, meteorological) vs station-specific artefacts.

# %%
fig = plot_timeseries(
    data,
    band=BAND,
    scale=SCALE,
    unit=UNIT,
    annotations=ANNOTATIONS,
    time_zone=TIMEZONE,
    resample_freq=RESAMPLE,
    day_night=DAY_NIGHT,
    show=False,
    outfile=os.path.join(SAVE_DIR, "multi_timeseries.png"),
)
plt.show()

# %% [markdown]
# ## 5. Side-by-side clock plots
#
# One subplot per station, all on the same radial scale.  Reveals which
# stations are dominated by diurnal anthropogenic noise (clear rush-hour
# peaks on weekdays) vs ambient/oceanic noise (no clear daily pattern).

# %%
n = len(data)
fig, axes = plt.subplots(1, n, figsize=(6 * n, 6),
                          subplot_kw={"projection": "polar"})
if n == 1:
    axes = [axes]

rmax_global = 0.0
for sid, df in data.items():
    s = localize_tz_and_reindex(df[[BAND]], RESAMPLE, TIMEZONE) * SCALE
    p95 = float(s.iloc[:, 0].quantile(0.95))
    rmax_global = max(rmax_global, p95 * 1.5)

for ax, (sid, df) in zip(axes, data.items()):
    s = localize_tz_and_reindex(df[[BAND]], RESAMPLE, TIMEZONE) * SCALE
    s.columns = ["rms"]
    stacked = stack_wday_time(s, scale=1.0)
    closed = stacked.copy()
    closed.loc[len(closed) + 1] = closed.iloc[0]
    closed.index = radial_hours(len(closed))
    closed.plot(ax=ax)
    ax.set_title(sid.split(".")[1], fontsize=11)
    clock24_plot_commons(ax, unit=UNIT)
    ax.set_rmax(rmax_global)
    ax.set_rmin(0)

plt.suptitle(f"24-h Clock Plots — Band [{BAND}] Hz", fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(os.path.join(SAVE_DIR, "multi_clockplots.png"),
            bbox_inches="tight", facecolor="w")
plt.show()

# %% [markdown]
# ## 6. Before/after clock plots (if SPLIT_DATE is set)
#
# Set `SPLIT_DATE` in the configuration cell to compare noise levels across an
# event or intervention.  Each station gets two polar panels.

# %%
if SPLIT_DATE:
    for sid, df in data.items():
        fig = plot_clockplot(
            {sid: df},
            band=BAND,
            scale=SCALE,
            unit=UNIT,
            split_date=SPLIT_DATE,
            time_zone=TIMEZONE,
            resample_freq=RESAMPLE,
            day_night=DAY_NIGHT,
            show=False,
            outfile=os.path.join(SAVE_DIR, f"{sid}_clockplot_split.png"),
        )
        plt.suptitle(f"{sid}  [{BAND}] Hz — split at {SPLIT_DATE}", y=1.02)
        plt.tight_layout()
        plt.show()
else:
    print("SPLIT_DATE not set — skipping before/after comparison.")

# %% [markdown]
# ## 7. Grid-maps for each station
#
# Cartesian view: x = calendar date, y = hour of day.  Best for exploring
# long time series where the polar hour-map becomes crowded.

# %%
for sid, df in data.items():
    fig = plot_gridmap(
        {sid: df},
        band=BAND,
        scale=SCALE,
        unit=UNIT,
        annotations=ANNOTATIONS,
        time_zone=TIMEZONE,
        resample_freq=RESAMPLE,
        day_night=DAY_NIGHT,
        show=False,
        outfile=os.path.join(SAVE_DIR, f"{sid}_gridmap.png"),
    )
    plt.show()

# %% [markdown]
# ## 8. Hour-maps for each station
#
# Polar view: radial = elapsed days, angular = hour of day.  Good for
# visually judging how seasonal or episodic changes affect the diurnal cycle.

# %%
for sid, df in data.items():
    fig = plot_hourmap(
        {sid: df},
        band=BAND,
        scale=SCALE,
        unit=UNIT,
        annotations=ANNOTATIONS,
        time_zone=TIMEZONE,
        resample_freq=RESAMPLE,
        day_night=DAY_NIGHT,
        show=False,
        outfile=os.path.join(SAVE_DIR, f"{sid}_hourmap.png"),
    )
    plt.show()

# %% [markdown]
# ## 9. Daily noise profiles (weekday stacks)
#
# For each station: median RMS per hour-of-day, stratified by weekday.
# A classic way to characterise the anthropogenic noise signature.

# %%
for sid, df in data.items():
    fig = plot_dailyplot(
        {sid: df},
        band=BAND,
        scale=SCALE,
        unit=UNIT,
        split_date=SPLIT_DATE,
        time_zone=TIMEZONE,
        resample_freq=RESAMPLE,
        day_night=DAY_NIGHT,
        show=False,
        outfile=os.path.join(SAVE_DIR, f"{sid}_dailyplot.png"),
    )
    plt.show()

# %% [markdown]
# ## 10. Direct NetCDF access (no DB required)
#
# For post-processing or sharing results outside a full MSNoise installation,
# you can load the NetCDF files directly using ``result.get_psd_rms()``.

# %%
# Example: load without a DB connection
# Reconstruct path components from MSNoiseResult for direct access
_db2 = connect()
_r = MSNoiseResult.from_ids(_db2, psd=PSD_ID, psd_rms=PSD_RMS_ID)
OUTPUT_FOLDER = _r.output_folder
_db2.close()


for sid in SEED_IDS:
    ds = _r.get_psd_rms(seed_id=sid)
    if ds is not None:
        print(f"{sid}:\n{ds}\n")
    else:
        print(f"{sid}: no data found")

# %% [markdown]
# ## 11. Batch CLI export
#
# The full workflow can be scripted from the command line.  After running
# `msnoise qc compute_psd_rms`, produce all figures with:
#
# ```bash
# for SID in BE.UCC..HHZ BE.MEM..HHZ BE.RCHB..HHZ; do
#     msnoise qc plot_psd_rms $SID \
#         --type all \
#         --band 1.0-10.0 \
#         --timezone Europe/Brussels \
#         --outfile figures/${SID}.png \
#         --no-show
# done
# ```
#
# Or for a single station with split-date comparison and event annotations:
#
# ```bash
# msnoise qc plot_psd_rms BE.UCC..HHZ \
#     --type clockplot \
#     --band 1.0-10.0 \
#     --split-date 2024-02-01 \
#     --annotate "2024-01-15=Maintenance,2024-03-01=Storm" \
#     --timezone Europe/Brussels \
#     --outfile BE.UCC..HHZ_clockplot.pdf \
#     --no-show
# ```

# %% [markdown]
# ---
# *Generated with `msnoise.plots.psd_rms` — ported and adapted from*
# *[SeismoRMS](https://github.com/ThomasLecocq/SeismoRMS).*
