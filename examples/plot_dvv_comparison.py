# -*- coding: utf-8 -*-
"""
Compare dv/v from MWCS, Stretching, and WCT
============================================

Plots the network-level dv/v time series produced by the three available
methods side by side — MWCS dt/t, Stretching, and Wavelet Cross-Transform
(WCT) — so their agreement and systematic differences can be assessed at a
glance.

This example is designed to run after the full MSNoise test workflow (all
compute steps completed).  Each ``(components, mov_stack)`` combination for
which at least one method has results produces its own figure.  Methods that
were not computed are silently skipped.

Each figure has two panels:

* **Top** — dv/v time series for each available method, with a ±1σ shaded
  band.
* **Bottom** — residuals relative to the first available method (e.g.
  MWCS − Stretching, MWCS − WCT) to highlight systematic offsets.
"""

import os

if "SPHINX_DOC_BUILD" in os.environ:
    if "MSNOISE_DOC" in os.environ:
        os.chdir(os.environ["MSNOISE_DOC"])

import matplotlib
matplotlib.use("agg")

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

from msnoise.api import connect
from msnoise.results import MSNoiseResult
from msnoise.core.config import lineage_to_plot_tag, build_plot_outfile

# ── palette & style ──────────────────────────────────────────────────────────
plt.style.use("ggplot")

METHOD_STYLE = {
    "mwcs_dtt_dvv":    dict(color="#1f77b4", label="MWCS",       lw=1.8, zorder=3),
    "stretching_dvv":  dict(color="#d62728", label="Stretching", lw=1.8, zorder=2),
    "wavelet_dtt_dvv": dict(color="#2ca02c", label="WCT",        lw=1.8, zorder=1),
}

# ── connect ───────────────────────────────────────────────────────────────────
db = connect()

# ── discover available DVV results ───────────────────────────────────────────
results_by_method = {}
for category in ("mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"):
    found = MSNoiseResult.list(db, category)
    if found:
        results_by_method[category] = found[0]
        tag = lineage_to_plot_tag(found[0].lineage_names)
        print(f"  {category}: {tag}")
    else:
        print(f"  {category}: not computed — skipping")

if not results_by_method:
    print("No DVV results found — run the full pipeline first.")
    db.close()
    raise SystemExit(0)

PAIR_TYPE = "CC"

# ── collect all (components, mov_stack) combinations across methods ───────────
all_keys = set()
raw = {}
for category, result in results_by_method.items():
    data = result.get_dvv(pair_type=PAIR_TYPE)
    filtered = {(comp, ms): ds
                for (pt, comp, ms), ds in data.items()
                if pt == PAIR_TYPE}
    raw[category] = filtered
    all_keys |= set(filtered.keys())

if not all_keys:
    print("No CC DVV data found.")
    db.close()
    raise SystemExit(0)

# ── one figure per (components, mov_stack) ────────────────────────────────────
for comp, ms in sorted(all_keys):
    ms_label = f"{ms[0]}-{ms[1]}"
    title    = f"dv/v comparison  |  comp={comp}  |  mov_stack={ms_label}"

    # Collect time series (converted to %)
    series = {}
    for category in ("mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"):
        if category not in raw or (comp, ms) not in raw[category]:
            continue
        ds = raw[category][(comp, ms)]
        if "mean" not in ds:
            continue
        mean_da = ds["mean"]
        std_da  = ds.get("std", ds.get("weighted_std", None))
        times   = pd.DatetimeIndex(mean_da.coords["times"].values)
        mean_s  = pd.Series(mean_da.values * 100.0, index=times, name=category)
        std_s   = (pd.Series(std_da.values * 100.0, index=times)
                   if std_da is not None else None)
        series[category] = (mean_s, std_s)

    if not series:
        print(f"  No data for comp={comp}, mov_stack={ms_label} — skipping")
        continue

    # ── figure: 2 rows ────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(12, 7))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.08)
    ax_dvv = fig.add_subplot(gs[0])
    ax_res = fig.add_subplot(gs[1], sharex=ax_dvv)

    # Top panel
    reference_s = None
    for category, (mean_s, std_s) in series.items():
        style = METHOD_STYLE[category]
        ax_dvv.plot(mean_s.index, mean_s.values,
                    color=style["color"], lw=style["lw"],
                    label=style["label"], zorder=style["zorder"])
        if std_s is not None:
            ax_dvv.fill_between(mean_s.index,
                                mean_s.values - std_s.values,
                                mean_s.values + std_s.values,
                                color=style["color"], alpha=0.15,
                                zorder=style["zorder"] - 1)
        if reference_s is None:
            reference_s = mean_s

    ax_dvv.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax_dvv.set_ylabel("dv/v (%)")
    ax_dvv.set_title(title)
    ax_dvv.legend(loc="upper right", framealpha=0.9)
    ax_dvv.tick_params(labelbottom=False)

    # Bottom panel — residuals vs first method
    ref_label = METHOD_STYLE[next(iter(series))]["label"]
    has_residual = False
    for category, (mean_s, _) in series.items():
        if mean_s is reference_s:
            continue
        style = METHOD_STYLE[category]
        common = reference_s.index.intersection(mean_s.index)
        if len(common) == 0:
            continue
        resid = reference_s.loc[common] - mean_s.loc[common]
        ax_res.plot(resid.index, resid.values,
                    color=style["color"], lw=1.4,
                    label=f"{ref_label} \u2212 {style['label']}",
                    zorder=style["zorder"])
        has_residual = True

    ax_res.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax_res.set_ylabel("Residual (%)")
    ax_res.set_xlabel("Date")
    if has_residual:
        ax_res.legend(loc="upper right", framealpha=0.9, fontsize=8)

    fig.autofmt_xdate()

    # ── save ──────────────────────────────────────────────────────────────────
    first_result = next(iter(results_by_method.values()))
    outfile = build_plot_outfile(
        "?.png", "dvv_comparison",
        first_result.lineage_names,
        components=comp,
        mov_stack=ms,
    )
    print(f"  Saving: {outfile}")
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    plt.show()

db.close()

# EOF
