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
# # PSD-RMS Visualisation — Quickstart
#
# This notebook demonstrates how to load the PSD-RMS output produced by
# `msnoise qc compute_psd_rms` and reproduce all plot types provided by the
# `msnoise.plots.psd_rms` module.
#
# ## Background
#
# The `psd_rms` workflow step reads the per-day PSD NetCDF files produced by
# `compute_psd`, integrates each PSD over user-defined frequency bands, and
# stores the resulting root-mean-square (RMS) amplitude time series in a
# NetCDF file:
#
# ```
# OUTPUT/
#   psd_1/
#     psd_rms_1/
#       _output/
#         BE.UCC..HHZ/
#           RMS.nc        ← times × bands DataArray
# ```
#
# Each column of `RMS.nc` is labelled `"fmin-fmax"` (e.g. `"1.0-10.0"`).  The
# amplitude unit depends on `psd_rms_type` in the step config (VEL, ACC, or
# DISP); in the examples below we assume displacement and scale to nanometres.
#
# ### Plot types
#
# | Command / function | Description |
# |--------------------|-------------|
# | `plot_timeseries`  | Calendar time-series with daytime median overlay |
# | `plot_clockplot`   | Polar 24-h clock, median per weekday/hour |
# | `plot_hourmap`     | Polar time×hour colour-mesh (radial = elapsed days) |
# | `plot_gridmap`     | Cartesian date×hour colour-mesh |
# | `plot_dailyplot`   | Hour-of-day stacked curves, one per weekday |

# %% [markdown]
# ## 1. Setup

# %%
# Standard imports
import matplotlib.pyplot as plt

# MSNoise database connection and result API
from msnoise.core.db import connect
from msnoise.results import MSNoiseResult

# RMS plotting helpers
from msnoise.plots.psd_rms import (
    DayNightWindow,
    load_rms, plot_timeseries, plot_clockplot,
    plot_hourmap, plot_gridmap, plot_dailyplot,
)

# %% [markdown]
# ## 2. Connect to the MSNoise database and load RMS data
#
# `load_rms` reads the NetCDF file(s) produced by `compute_psd_rms`.
# It returns a dict `{seed_id: DataFrame}` where the DataFrame columns are
# the frequency band labels and the index is a DatetimeIndex (UTC).

# %%
db = connect()

# ── Edit these to match your project ─────────────────────────────────────────
SEED_IDS  = ["BE.UCC..HHZ"]   # one or more SEED IDs
PSD_ID    = 1                  # config-set number of your psd step     (psd_1)
PSD_RMS_ID = 1                 # config-set number of your psd_rms step (psd_rms_1)
# ─────────────────────────────────────────────────────────────────────────────

# MSNoiseResult resolves the correct output path (psd_N/psd_rms_N/_output/)
# without requiring hard-coded step-name strings.
result = MSNoiseResult.from_ids(db, psd=PSD_ID, psd_rms=PSD_RMS_ID)
data = load_rms(result, SEED_IDS)
db.close()

for sid, df in data.items():
    print(f"{sid}: {len(df)} rows, bands = {list(df.columns)}")
    print(df.head(3))

# %% [markdown]
# ## 3. Choose a frequency band
#
# The `band` parameter selects which column of the RMS DataFrame to plot.
# Use any of the labels printed above, or leave it `None` to auto-select the
# first available band.

# %%
# Pick a band — adjust to match your psd_rms_frequency_ranges config
BAND = None   # e.g. "1.0-10.0"

# Amplitude scaling: 1e9 converts SI (m) to nanometres
SCALE = 1e9
UNIT  = "nm"

# Local time zone — used for clockplot/dailyplot hour-of-day grouping
TIMEZONE = "Europe/Brussels"   # or "UTC"

# Day / Night window — controls shading and the daytime median overlay.
# Adjust to match your station's local anthropogenic noise schedule.
#   DayNightWindow(day_start=0, day_end=24)  — disables shading entirely
#   DayNightWindow(day_start=22, day_end=6)  — nighttime-city station
DAY_NIGHT = DayNightWindow(day_start=7, day_end=19)

# %% [markdown]
# ## 4. Time-series plot
#
# Shows the RMS amplitude over the full time range.  Green shading marks
# business days; dashed lines show the 6–16 h daytime median per day.

# %%
fig = plot_timeseries(
    data,
    band=BAND,
    scale=SCALE,
    unit=UNIT,
    day_night=DAY_NIGHT,
    time_zone=TIMEZONE,
    resample_freq="30min",
    agg_func="mean",
    # Optional event annotations:
    # annotations={"2024-01-15": "Maintenance", "2024-02-01": "New sensor"},
    show=False,
)
plt.show()

# %% [markdown]
# ## 5. Clock plot
#
# Polar plot with 24 hours around the circumference and one curve per weekday.
# Reveals diurnal patterns (e.g. rush-hour noise peaks).
#
# Supply `split_date` to show before/after comparison in two panels.

# %%
fig = plot_clockplot(
    data,
    band=BAND,
    scale=SCALE,
    unit=UNIT,
    day_night=DAY_NIGHT,
    time_zone=TIMEZONE,
    resample_freq="30min",
    agg_func="mean",
    # split_date="2024-02-01",   # uncomment for before/after panels
    show=False,
)
plt.show()

# %% [markdown]
# ## 6. Hour-map (polar colour-mesh)
#
# Radial axis = elapsed days since the start of the dataset.
# Angular axis = hour of day.  Colour encodes RMS amplitude.
# Reveals how diurnal patterns evolve over time.

# %%
fig = plot_hourmap(
    data,
    band=BAND,
    scale=SCALE,
    unit=UNIT,
    day_night=DAY_NIGHT,
    time_zone=TIMEZONE,
    resample_freq="30min",
    agg_func="mean",
    show=False,
)
plt.show()

# %% [markdown]
# ## 7. Grid-map (Cartesian colour-mesh)
#
# Calendar date on x-axis, hour of day on y-axis.
# Easier to read for long time series than the polar hour-map.

# %%
fig = plot_gridmap(
    data,
    band=BAND,
    scale=SCALE,
    unit=UNIT,
    day_night=DAY_NIGHT,
    time_zone=TIMEZONE,
    resample_freq="30min",
    agg_func="mean",
    show=False,
)
plt.show()

# %% [markdown]
# ## 8. Daily plot (hour-of-day stacked by weekday)
#
# Each curve = median RMS for one weekday at each hour.
# Solid = full dataset; dashed = after `split_date` (if provided).

# %%
fig = plot_dailyplot(
    data,
    band=BAND,
    scale=SCALE,
    unit=UNIT,
    day_night=DAY_NIGHT,
    time_zone=TIMEZONE,
    resample_freq="30min",
    agg_func="mean",
    # split_date="2024-02-01",
    show=False,
)
plt.show()

# %% [markdown]
# ## 9. Saving figures
#
# Pass `outfile` to any of the plot functions to save to disk.
# `plot_type="all"` in the CLI will append a suffix for each type.

# %%
# Example: save all plot types to PNG
for pt, fn in [
    ("timeseries", "psd_rms_timeseries.png"),
    ("clockplot",  "psd_rms_clockplot.png"),
    ("hourmap",    "psd_rms_hourmap.png"),
    ("gridmap",    "psd_rms_gridmap.png"),
    ("dailyplot",  "psd_rms_dailyplot.png"),
]:
    funcs = {
        "timeseries": plot_timeseries,
        "clockplot":  plot_clockplot,
        "hourmap":    plot_hourmap,
        "gridmap":    plot_gridmap,
        "dailyplot":  plot_dailyplot,
    }
    try:
        funcs[pt](data, band=BAND, scale=SCALE, unit=UNIT,
                  time_zone=TIMEZONE, show=False, outfile=fn)
        print(f"Saved {fn}")
# %% [markdown]
# ## 11. CLI equivalent — with resample/agg options
#
# ```bash
# msnoise qc plot_psd_rms BE.UCC..HHZ --type timeseries --resample-freq 1h --agg-func median
# msnoise qc plot_psd_rms BE.UCC..HHZ --type all -R 30min -A mean --no-show
# ```
    except Exception as e:
        print(f"Skipped {fn}: {e}")

# %% [markdown]
# ## 10. CLI equivalent
#
# All of the above can be reproduced directly from the command line once
# `compute_psd_rms` has been run:
#
# ```bash
# # Time-series for one station
# msnoise qc plot_psd_rms BE.UCC..HHZ --type timeseries --band 1.0-10.0
#
# # Clock plot in local time, save to PDF
# msnoise qc plot_psd_rms BE.UCC..HHZ \
#     --type clockplot --timezone Europe/Brussels \
#     --outfile noise.pdf --no-show
#
# # All plot types at once
# msnoise qc plot_psd_rms BE.UCC..HHZ --type all --outfile station.png --no-show
#
# # Before/after comparison with event annotation
# msnoise qc plot_psd_rms BE.UCC..HHZ \
#     --type clockplot --split-date 2024-01-15 \
#     --annotate "2024-01-15=Maintenance,2024-02-01=New sensor"
# ```
