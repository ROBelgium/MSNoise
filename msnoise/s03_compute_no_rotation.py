r"""Computation of inter-station and single-station cross-correlation functions.

This module implements the ``msnoise cc compute`` worker.  It groups *jobs*
marked **T**odo in the database by day, marks them **I**n-progress, processes
all sliding windows and filters for that day, then marks them **D**one.
Multiple instances can run in parallel because each worker atomically claims
its own batch of jobs.


Overview
--------

The worker operates on 2-D arrays of shape ``(n_stations, N)`` — one row per
pre-processed trace, one column per sample — and processes all station pairs
simultaneously for each time window.  Two independent correlation algorithms
are available:

* **CC** — classic Generalised Normalised Cross-Correlation (GNCC), computed
  in the frequency domain via FFT (see `Cross-Correlation (CC)`_ below).
* **PCC** — Phase Cross-Correlation v2 (PCC2), a transient-robust
  alternative that operates in the time domain before taking any FFT (see
  `Phase Cross-Correlation (PCC2)`_ below).

The algorithm is selected per correlation mode through three independent
configuration parameters: ``cc_type`` (inter-station CC),
``cc_type_single_station_AC`` (auto-correlation) and
``cc_type_single_station_SC`` (same-station cross-component).


Cross-Correlation (CC)
----------------------

The classic ambient-noise cross-correlation (Shapiro & Campillo 2004;
Bensen et al. 2007) is computed using the tiled-batch implementation
:func:`~msnoise.core.compute.myCorr2`.

**Shared pre-processing per time window** (all modes — CC, AC, SC):

1. **Temporal normalisation** — optional Winsorising (clipping to
   ``winsorizing`` × RMS) applied either *before* whitening (default) or
   *after* (``clip_after_whiten = Y``).

2. **Time-domain taper** — a cosine (Hann) taper of width
   ``cc_taper_fraction × N`` is applied symmetrically to both ends to
   suppress spectral leakage.

3. **Single bandpass per filter band** — one bandpassed copy ``_data_bp``
   is computed from the raw tapered data at the start of each filter
   iteration and reused by all three modes (AC, CC, SC) and both
   algorithms (CC, PCC).  No mode bandpasses the data a second time.

**Processing chain — CC algorithm:**

For each requested pair :math:`(i, j)` in each filter band:

- ``whitening = N`` — ``_data_bp`` fed directly into the cross-spectrum:

  .. math::

      C_{ij}[k] = X_i^*[k] \cdot X_j[k], \quad
      X_i = \mathcal{F}\{x_i^\text{bp}[n]\}

- ``whitening ≠ N`` — raw (non-bandpassed) data fed into FFT, then
  spectral whitening applied in-place by
  :func:`~msnoise.core.compute.whiten2` (which itself bandlimits to
  ``[f_\text{low}, f_\text{high}]``, making a prior bandpass redundant):

  .. math::

      \tilde{X}_i[k] = \operatorname{whiten2}(\mathcal{F}\{x_i[n]\})

  Three whitening shapes are available via ``whitening_type``:

  - *Brutal* (default) — ``|\tilde{X}[k]| = 1`` inside the passband
    with a cosine-tapered transition.
  - *HANN* — Hann-weighted one-bit normalisation inside the passband.
  - *PSD* — divide by smoothed PSD, then clip outlier bins at the
    5th–95th percentile.

The IFFT, folding and optional normalisation follow in all cases:

.. math::

    c_{ij}[\tau] = \mathcal{F}^{-1}\{C_{ij}[k]\} / N

**Normalisation** (``cc_normalisation``): ``POW`` — divide by the
product of RMS energies :math:`(e_i \cdot e_j)`; ``MAX`` — divide by
the CCF maximum; ``ABSMAX`` — divide by the absolute maximum.

**References**

Bensen, G. D., Ritzwoller, M. H., Barmin, M. P., Levshin, A. L., Lin, F.,
Moschetti, M. P., Shapiro, N. M., & Yang, Y. (2007).
Processing seismic ambient noise data to obtain reliable broad-band surface
wave dispersion measurements.
*Geophysical Journal International*, 169(3), 1239–1260.
https://doi.org/10.1111/j.1365-246X.2007.03374.x

Shapiro, N. M., & Campillo, M. (2004).
Emergence of broadband Rayleigh waves from correlations of the ambient
seismic noise.
*Geophysical Research Letters*, 31(7), L07614.
https://doi.org/10.1029/2004GL019491


Phase Cross-Correlation (PCC2)
------------------------------

Phase Cross-Correlation (Schimmel 1999; Schimmel et al. 2011) is an
alternative to the classic GNCC that discards amplitude information entirely
at every sample, making it intrinsically robust against transient noise
(earthquakes, instrumental glitches) *without* requiring explicit temporal
normalisation such as one-bit or clipping.  MSNoise implements PCC version 2
(PCC2), the FFT-accelerated formulation introduced by Ventosa et al. (2019,
2023) as part of the FastPCC package.  The implementation is self-contained
in :func:`~msnoise.core.compute.pcc_xcorr` and does **not** depend on the
``phasecorr`` or ``fastpcc`` external packages.

**Processing chain — PCC2 algorithm:**

For each trace the shared ``_data_bp`` (already bandpassed to
:math:`[f_\text{low}, f_\text{high}]`) is used directly:

1. **Optional spectral whitening** — if ``whitening ≠ N``, ``_data_bp``
   is FFT-whitened within the band and transformed back to the time domain.
   This distributes energy evenly across the passband *before* the per-sample
   amplitude normalisation, so the resulting phase signal reflects broadband
   coherence rather than being dominated by the most energetic frequency.

2. **Analytic signal and amplitude normalisation** — the complex analytic
   signal :math:`x_i^{(a)}[n]` is computed via the Hilbert transform
   (:func:`scipy.signal.hilbert`).  Each sample is divided by its own
   amplitude to produce the *phase signal* :math:`\varphi_i`:

   .. math::

       \varphi_i[n] = \frac{x_i^{(a)}[n]}{|x_i^{(a)}[n]| + \varepsilon}

   where :math:`\varepsilon = 10^{-6} \max_n |x_i^{(a)}[n]|` is a numerical
   stability floor (matching FastPCC's ``AmpNormf`` convention).
   By construction :math:`|\varphi_i[n]| \leq 1`, and amplitude transients
   (earthquakes, glitches) are reduced to exactly the same per-sample weight
   as the ambient noise — no explicit temporal normalisation is needed.

3. **Zero-padding** — each phase signal is zero-padded to
   :math:`N_z = \text{next\_fast\_len}(N + \tau_\text{max})` to compute a
   linear (non-circular) cross-correlation.

4. **FFT cross-spectrum and IFFT** — all phase signals are pre-transformed
   once, then for each pair :math:`(i, j)`:

   .. math::

       \text{PCC2}_{ij}[\tau] = \mathcal{F}^{-1}\!\left\{
           \Phi_i^*[k] \cdot \Phi_j[k]
       \right\} \Big/ N

   Division by :math:`N` (not :math:`N_z`) gives a peak amplitude of
   approximately 1 for identical signals.

5. **Lag unwrapping** — positive lags :math:`0 \ldots \tau_\text{max}`
   come from ``full[0:\tau_\text{max}+1]``, negative lags
   :math:`-\tau_\text{max} \ldots -1` from ``full[N_z-\tau_\text{max}:N_z]``.

6. **Normalisation** — ``MAX`` or ``ABSMAX`` as for CC; ``POW`` is silently
   ignored because amplitudes are discarded in step 2.

**References**

Schimmel, M. (1999).
Phase cross-correlations: design, comparisons, and applications.
*Bulletin of the Seismological Society of America*, 89(5), 1366–1378.
https://doi.org/10.1785/BSSA0890051366

Schimmel, M., Stutzmann, E., & Gallart, J. (2011).
Using instantaneous phase coherence for signal extraction from ambient noise
data at a local to a global scale.
*Geophysical Journal International*, 184(1), 494–506.
https://doi.org/10.1111/j.1365-246X.2010.04861.x

Ventosa, S., Schimmel, M., & Stutzmann, E. (2019).
Towards the processing of large data volumes with phase cross-correlation.
*Seismological Research Letters*, 90(4), 1663–1669.
https://doi.org/10.1785/0220190022

Ventosa, S., & Schimmel, M. (2023).
FastPCC: Fast phase cross-correlation algorithm for large seismic datasets.
*IEEE Transactions on Geoscience and Remote Sensing*, 61, 1–17.
https://doi.org/10.1109/TGRS.2023.3294302


Stacking daily windows
-----------------------

For each ``corr_duration``-long window, and for each configured filter, the
CCF (CC or PCC2) is computed and accumulated.  If ``keep_all = Y`` the
individual window CCFs are written to disk.  By default (``keep_days = Y``),
all windows for the day are stacked into a single daily CCF using either a
linear mean or Phase-Weighted Stack (PWS; Schimmel & Paulssen 1997):

.. note::

    PWS is provided as an experimental option.  It has not been
    systematically cross-validated.  Use with caution.

Schimmel, M., & Paulssen, H. (1997).
Noise reduction and detection of weak, coherent signals through
phase-weighted stacks.
*Geophysical Journal International*, 130(2), 497–505.
https://doi.org/10.1111/j.1365-246X.1997.tb05664.x


Configuration Parameters
------------------------

The following parameters (modifiable via ``msnoise admin``) are used for
this step:

* |cc.components_to_compute|
* |cc.components_to_compute_single_station|
* |cc.cc_sampling_rate|
* |cc.cc_normalisation|
* |cc.cc_type|
* |cc.cc_type_single_station_AC|
* |cc.cc_type_single_station_SC|
* |cc.cc_taper_fraction|
* |cc.clip_after_whiten|
* |cc.overlap|
* |cc.maxlag|
* |cc.corr_duration|
* |cc.winsorizing|
* |cc.whitening|
* |cc.whitening_type|
* |cc.keep_all|
* |cc.keep_days|
* |cc.stack_method|
* |cc.pws_timegate|
* |cc.pws_power|
* |global.hpc|


Computing the Cross-Correlations
---------------------------------

The following two diagrams describe the complete execution flow of
``msnoise cc compute_cc``.  Both are rendered from the DOT source embedded
below, so they stay in sync with the code automatically.

**Figure 1 — Outer job loop** shows how days are consumed, time windows
slid, and per-filter results accumulated and saved.

.. graphviz::

   digraph cc_outer {
       graph [
           rankdir=LR
           fontname="Helvetica,Arial,sans-serif"
           fontsize=12
           nodesep=0.4
           ranksep=0.6
           bgcolor="white"
           pad=0.4
           label="Figure 1 – Outer job loop (s03_compute_no_rotation)"
           labelloc=t
           labeljust=l
       ]
       node [fontname="Helvetica,Arial,sans-serif" fontsize=11 style=filled penwidth=1.5]
       edge [fontname="Helvetica,Arial,sans-serif" fontsize=10 color="#444444" penwidth=1.2]

       node [shape=oval fillcolor="#263238" fontcolor=white color="#263238"]
       START [label="begin"]
       END   [label="end"]

       node [shape=diamond fillcolor="#E3F2FD" fontcolor="#0D1B2A" color="#1565C0"]
       D_JOB   [label="New CC\njob?"]
       D_SLIDE [label="Next time\nwindow?"]
       D_GAPS  [label="Any trace\nhas gaps?"]
       D_SHORT [label="All traces\ntoo short?"]
       D_WTYPE [label="whitening_type\n== PSD?"]
       D_WMODE [label="whitening\n!= N?"]
       D_FILTER[label="Next filter\nband?"]

       node [shape=diamond fillcolor="#E8EAF6" fontcolor="#0D1B2A" color="#283593"]
       D_KEEPALL [label="keep_all?"]
       D_KEEPDAY [label="keep_days?"]

       node [shape=box fillcolor="#FAFAFA" fontcolor="#0D1B2A" color="#607D8B" style="filled,rounded"]
       P_STREAM    [label="Load preprocessed\nstreams for day"]
       P_SLIDE     [label="Slide window\n(corr_duration, overlap)"]
       P_RMGAP     [label="Remove gapped\ntraces"]
       P_PREP      [label="Detrend / demean\nWinsorise (if !clip_after_whiten)\nCosine taper (cc_taper_fraction %)"]
       P_PSD       [label="Precompute Welch\nPSD per station"]
       P_BP_N      [label="Bandpass _data\n[f_low, f_high]"]
       P_FIDX      [label="Build pair indices\n(AC · CC · SC)"]
       P_ACC       [label="Accumulate window CCF\ninto allcorr dict"]

       node [shape=box fillcolor="#FFF8E1" fontcolor="#0D1B2A" color="#F9A825" style="filled,rounded,bold"]
       P_PERFILTER [label="▶  Per-filter processing\n(see Figure 2)"]

       node [shape=box fillcolor="#E8EAF6" fontcolor="#0D1B2A" color="#283593" style="filled,rounded"]
       P_KEEPALL_S [label="Save all-window CCFs\n(xr_save_ccf_all)"]
       P_STACK     [label="Stack windows\n(linear or PWS)"]
       P_SAVEDAY   [label="Save daily CCF\n(xr_save_ccf_daily)"]
       P_JOBS      [label="massive_update_job → Done\npropagate_downstream"]

       // Force two rows via an invisible break node
       node [shape=point width=0 height=0 style=invis]
       BREAK

       // Row 1: start → window loop → prep/whitening
       { rank=same; START; D_JOB; END }
       { rank=same; P_STREAM; D_SLIDE }
       { rank=same; P_SLIDE; P_RMGAP; D_GAPS; D_SHORT }
       { rank=same; P_PREP; D_WTYPE; P_PSD }

       // Row break sentinel — pin BREAK at same rank as P_FIDX
       // Row 2: filter loop → per-filter → save
       { rank=same; P_FIDX; D_FILTER }
       { rank=same; BREAK; D_WMODE; P_BP_N; P_PERFILTER; P_ACC }
       { rank=same; D_KEEPALL; P_KEEPALL_S; D_KEEPDAY; P_STACK; P_SAVEDAY; P_JOBS }

       // Invisible structural edge to force rank break after row 1
       D_WTYPE -> BREAK [style=invis]
       BREAK -> D_WMODE [style=invis]

       START    -> D_JOB
       D_JOB    -> END          [label="no"]
       D_JOB    -> P_STREAM     [label="yes"]
       P_STREAM -> D_SLIDE

       D_SLIDE  -> D_KEEPALL    [label="no more\nwindows" constraint=false]
       D_SLIDE  -> P_SLIDE      [label="yes"]
       P_SLIDE  -> P_RMGAP
       P_RMGAP  -> D_GAPS
       D_GAPS   -> D_SLIDE      [label="yes" constraint=false]
       D_GAPS   -> D_SHORT      [label="no"]
       D_SHORT  -> D_SLIDE      [label="yes" constraint=false]
       D_SHORT  -> P_PREP       [label="no"]

       P_PREP   -> D_WTYPE
       D_WTYPE  -> P_PSD        [label="yes (PSD)"]
       D_WTYPE  -> P_FIDX       [label="no"]
       P_PSD    -> P_FIDX

       P_FIDX   -> D_FILTER
       D_FILTER -> D_KEEPALL    [label="no more filters" constraint=false]
       D_FILTER -> D_WMODE      [label="yes"]
       D_WMODE  -> P_BP_N       [label="no\n(whitening=N)"]
       D_WMODE  -> P_PERFILTER  [label="yes"]
       P_BP_N   -> P_PERFILTER
       P_PERFILTER -> P_ACC
       P_ACC    -> D_FILTER     [constraint=false]

       D_KEEPALL   -> P_KEEPALL_S [label="yes"]
       D_KEEPALL   -> D_KEEPDAY   [label="no"]
       P_KEEPALL_S -> D_KEEPDAY
       D_KEEPDAY   -> P_STACK     [label="yes"]
       D_KEEPDAY   -> P_JOBS      [label="no"]
       P_STACK     -> P_SAVEDAY
       P_SAVEDAY   -> P_JOBS
       P_JOBS      -> D_JOB       [constraint=false]
   }

**Figure 2 — Per-filter processing** shows the three parallel correlation
modes (AC, CC, SC) and the CC/PCC2 algorithm branch within each.  The entry
node carries ``_data_bp`` (bandpassed once per filter iteration, shared by
all modes) and ``_data_raw`` (used only by CC/SC with ``whitening != N``,
where ``whiten2`` applies the spectral window itself).

.. graphviz::

   digraph cc_perfilter {
       graph [
           rankdir=TB
           fontname="Helvetica,Arial,sans-serif"
           fontsize=12
           nodesep=0.4
           ranksep=0.5
           bgcolor="white"
           pad=0.4
           label="Figure 2 – Per-filter processing (each time window × filter band)"
           labelloc=t
           labeljust=l
       ]
       node [fontname="Helvetica,Arial,sans-serif" fontsize=11 style=filled penwidth=1.5]
       edge [fontname="Helvetica,Arial,sans-serif" fontsize=10 color="#444444" penwidth=1.2]

       node [shape=oval fillcolor="#607D8B" fontcolor=white color="#607D8B"]
       ENTRY [label="enter\n_data_bp  (bandpassed once)\n_data_raw (for CC+whiten path)"]
       EXIT  [label="return\n(to outer loop)"]

       node [shape=box fillcolor="#F0F4F8" fontcolor="#0D1B2A" color="#607D8B" style="filled,rounded"]
       P_ACC [label="Accumulate CCF\ninto allcorr"]

       subgraph cluster_ac {
           label="Auto-Correlation (AC)" labelloc=t
           color="#F9A825" style=rounded penwidth=2.0

           node [shape=diamond fillcolor="#FFF9C4" fontcolor="#0D1B2A" color="#F9A825"]
           D_AC     [label="AC pairs?"]
           D_ACTYPE [label="cc_type\nsingle_AC?"]
           D_CAW_AC [label="clip_after\nwhiten?"]

           node [shape=box fillcolor="#FFFDE7" fontcolor="#0D1B2A" color="#F9A825" style="filled,rounded"]
           P_AC_CAW [label="Winsorise (time-series)"]

           node [shape=box fillcolor="#FFF9C4" fontcolor="#0D1B2A" color="#F9A825" style="filled,rounded"]
           P_AC_CC_F  [label="FFT → nfft"]
           P_AC_CC_W  [label="whiten2  (inplace)\nif whitening != N"]
           P_AC_CC_E  [label="Compute RMS energy"]
           P_AC_CC_C  [label="myCorr2\n(IFFT · fold · normalise)"]
           P_AC_PCC_W [label="FFT → whiten2 → IFFT\n(optional, whitening != N)"]
           P_AC_PCC_C [label="pcc_xcorr\nHilbert → φ(t) → FFT\ncross-spec → IFFT / N"]
       }

       subgraph cluster_cc {
           label="Cross-Correlation (CC)" labelloc=t
           color="#2E7D32" style=rounded penwidth=2.0

           node [shape=diamond fillcolor="#E8F5E9" fontcolor="#0D1B2A" color="#2E7D32"]
           D_CC     [label="CC pairs?"]
           D_CCTYPE [label="cc_type?"]
           D_CAW_CC [label="clip_after\nwhiten?"]
           D_CC_WN  [label="whitening\n!= N?"]

           node [shape=box fillcolor="#E8F5E9" fontcolor="#0D1B2A" color="#2E7D32" style="filled,rounded"]
           P_CC_RAW   [label="FFT(_data_raw)\n→ whiten2 (inplace)"]
           P_CC_BP_F  [label="FFT(_data_bp)"]
           P_CC_CAW   [label="Winsorise (FFT domain)"]
           P_CC_E     [label="Compute RMS energy"]
           P_CC_CORR  [label="myCorr2\n(IFFT · fold · normalise)"]
           P_CC_PW    [label="FFT → whiten2 → IFFT\n(optional, whitening != N)"]
           P_CC_PCORR [label="pcc_xcorr\nHilbert → φ(t) → FFT\ncross-spec → IFFT / N"]
       }

       subgraph cluster_sc {
           label="Same-station Cross-Component (SC)" labelloc=t
           color="#6A1B9A" style=rounded penwidth=2.0

           node [shape=diamond fillcolor="#F3E5F5" fontcolor="#0D1B2A" color="#6A1B9A"]
           D_SC     [label="SC pairs?"]
           D_SCTYPE [label="cc_type\nsingle_SC?"]
           D_CAW_SC [label="clip_after\nwhiten?"]
           D_SC_WN  [label="whitening\n!= N?"]

           node [shape=box fillcolor="#F3E5F5" fontcolor="#0D1B2A" color="#6A1B9A" style="filled,rounded"]
           P_SC_RAW   [label="FFT(_data_raw)\n→ whiten2 (inplace)"]
           P_SC_BP_F  [label="FFT(_data_bp)"]
           P_SC_CAW   [label="Winsorise (FFT domain)"]
           P_SC_E     [label="Compute RMS energy"]
           P_SC_CORR  [label="myCorr2\n(IFFT · fold · normalise)"]
           P_SC_PW    [label="FFT → whiten2 → IFFT\n(optional, whitening != N)"]
           P_SC_PCORR [label="pcc_xcorr\nHilbert → φ(t) → FFT\ncross-spec → IFFT / N"]
       }

       ENTRY -> D_AC

       // AC
       D_AC     -> D_CC      [label="no"]
       D_AC     -> D_CAW_AC  [label="yes"]
       D_CAW_AC -> P_AC_CAW  [label="yes"]
       D_CAW_AC -> D_ACTYPE  [label="no"]
       P_AC_CAW -> D_ACTYPE
       D_ACTYPE -> P_AC_CC_F  [label="CC"]
       D_ACTYPE -> P_AC_PCC_W [label="PCC"]
       P_AC_CC_F  -> P_AC_CC_W
       P_AC_CC_W  -> P_AC_CC_E
       P_AC_CC_E  -> P_AC_CC_C
       P_AC_PCC_W -> P_AC_PCC_C
       P_AC_CC_C  -> P_ACC
       P_AC_PCC_C -> P_ACC

       P_ACC -> D_CC

       // CC
       D_CC     -> D_SC      [label="no"]
       D_CC     -> D_CCTYPE  [label="yes"]
       D_CCTYPE -> D_CC_WN   [label="CC"]
       D_CCTYPE -> P_CC_PW   [label="PCC"]
       D_CC_WN  -> P_CC_RAW  [label="yes\n(_data_raw)"]
       D_CC_WN  -> P_CC_BP_F [label="no\n(_data_bp)"]
       P_CC_RAW -> D_CAW_CC
       P_CC_BP_F -> D_CAW_CC
       D_CAW_CC -> P_CC_CAW  [label="yes"]
       D_CAW_CC -> P_CC_E    [label="no"]
       P_CC_CAW -> P_CC_E
       P_CC_E   -> P_CC_CORR
       P_CC_PW  -> P_CC_PCORR
       P_CC_CORR  -> P_ACC
       P_CC_PCORR -> P_ACC

       P_ACC -> D_SC

       // SC
       D_SC     -> EXIT      [label="no"]
       D_SC     -> D_SCTYPE  [label="yes"]
       D_SCTYPE -> D_SC_WN   [label="CC"]
       D_SCTYPE -> P_SC_PW   [label="PCC"]
       D_SC_WN  -> P_SC_RAW  [label="yes\n(_data_raw)"]
       D_SC_WN  -> P_SC_BP_F [label="no\n(_data_bp)"]
       P_SC_RAW -> D_CAW_SC
       P_SC_BP_F -> D_CAW_SC
       D_CAW_SC -> P_SC_CAW  [label="yes"]
       D_CAW_SC -> P_SC_E    [label="no"]
       P_SC_CAW -> P_SC_E
       P_SC_E   -> P_SC_CORR
       P_SC_PW  -> P_SC_PCORR
       P_SC_CORR  -> P_ACC
       P_SC_PCORR -> P_ACC
       P_ACC -> EXIT [constraint=false]
   }


Usage
-----

To run this script:

.. code-block:: sh

    $ msnoise cc compute_cc

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc compute_cc

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

"""

import itertools
import logging
import os
import time

import numpy as np

from .core.db import connect, get_logger
from .core.config import get_config_set_details
from .core.workflow import (get_filter_steps_for_cc_step, get_next_lineage_batch, get_t_axis, is_next_job_for_step, massive_update_job, propagate_downstream)
from .core.signal import stack, winsorizing, get_preprocessed_stream
from .core.io import save_daily_ccf, xr_save_ccf_all
from .core.compute import myCorr2
from .core.compute import whiten2
from .core.compute import pcc_xcorr



import scipy.signal
import scipy.fft as sf
from scipy.fft import next_fast_len
from obspy.signal.filter import bandpass


def main(loglevel="INFO", chunk_size=0):
    global logger
    logger = get_logger('msnoise.cc', loglevel, with_pid=True)
    logger.info('*** Starting: Compute CC ***')

    # Connection to the DB
    db = connect()

    logger.debug("Checking if there are jobs to do")
    if chunk_size > 0:
        logger.info(f"CC chunk_size={chunk_size}: each worker claims up to {chunk_size} pairs per day")

    while is_next_job_for_step(db, step_category="cc"):
        batch = get_next_lineage_batch(db, step_category="cc", group_by="day_lineage",
                                       chunk_size=chunk_size, loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        params = batch["params"]
        lineage_names_upstream = batch["lineage_names_upstream"]
        lineage_names = batch["lineage_names"]  # full: e.g. [preprocess_1, cc_1]
        step = batch["step"]

        stations = []
        pairs = []
        refs = []

        for job in jobs:
            refs.append(job.ref)
            pairs.append(job.pair)
            netsta1, netsta2 = job.pair.split(':')
            stations.append(netsta1)
            stations.append(netsta2)
            goal_day = job.day

        stations = np.unique(stations)
        t_axis = get_t_axis(params)

        # Filters:
        filter_steps = get_filter_steps_for_cc_step(db, step.step_id)
        filters = []
        for filter in filter_steps:
            filters.append([filter.step_name, get_config_set_details(db, filter.category, filter.set_number, format='AttribDict')])
        logger.info("New CC Job: %s (%i pairs with %i stations)" %
                     (goal_day, len(pairs), len(stations)))
        jt = time.time()

        # Read per-station preprocessed files — one file per NET.STA.LOC
        # under _output/<goal_day>/<NET.STA.LOC>.mseed (v2 layout).
        _preprocess_step = lineage_names_upstream[-1] if lineage_names_upstream else ""
        _preprocess_out  = os.path.join(params.global_.output_folder,
                                        *lineage_names_upstream[:-1])
        stream = get_preprocessed_stream(
            _preprocess_out, _preprocess_step, goal_day, stations
        )
        # TODO PREPROCESS IF THE "PREPROCESS_ON_THE_FLY" config?
        # comps = np.unique(comps)
        # with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        #     stream = executor.submit(preprocess, stations, comps, goal_day, params, responses, loglevel).result()
        # logger.info("Received preprocessed traces")


        # stream = preprocess(db, stations, comps, goal_day, params, responses)
        if not len(stream):
            logger.debug("Not enough data for this day !")
            logger.debug("Marking jobs 'Fail' and continuing with next !")
            massive_update_job(db, jobs, "F")
            continue
        # print '##### STREAMS ARE ALL PREPARED AT goal Hz #####'
        dt = 1. / params.cc.cc_sampling_rate
        logger.debug("Starting slides")
        start_processing = time.time()
        allcorr = {}
        for tmp in stream.slide(params.cc.corr_duration,
                                params.cc.corr_duration * (1 - params.cc.overlap)):
            logger.debug("Processing %s - %s" % (tmp[0].stats.starttime,
                                                 tmp[0].stats.endtime))
            tmp = tmp.copy().sort()

            channels_to_remove = []
            for gap in tmp.get_gaps(min_gap=0):
                if gap[-2] > 0:
                    channels_to_remove.append(
                        ".".join([gap[0], gap[1], gap[2], gap[3]]))

            for chan in np.unique(channels_to_remove):
                logger.debug("%s contains gap(s), removing it" % chan)
                net, sta, loc, chan = chan.split(".")
                for tr in tmp.select(network=net,
                                     station=sta,
                                     location=loc,
                                     channel=chan):
                    tmp.remove(tr)
            if len(tmp) == 0:
                logger.debug("No traces without gaps")
                continue

            base = np.amax([tr.stats.npts for tr in tmp])
            if base <= (params.cc.maxlag*params.cc.cc_sampling_rate*2+1):
                logger.debug("All traces shorter are too short to export"
                              " +-maxlag")
                continue

            for tr in tmp:
                if tr.stats.npts != base:
                    tmp.remove(tr)
                    logger.debug("%s trace is too short, removing it" % tr.id)

            if len(tmp) == 0:
                logger.debug("No traces left in slice")
                continue

            nfft = next_fast_len(tmp[0].stats.npts)
            tmp.detrend("demean")

            # Spectral end-taper width: ~0.5% of nfft, minimum 50 bins.
            # Tapers the edges of the FFT output before whitening to
            # suppress spectral leakage at DC and Nyquist.
            napod = max(50, nfft // 200)

            data = np.asarray([tr.data for tr in tmp])
            names = [tr.id.split(".") for tr in tmp]

            if not params.cc.clip_after_whiten:
                # logger.debug("Winsorizing (clipping) data before whiten")
                data = winsorizing(data, params) #inplace

            # Time-domain cosine taper fraction (config: cc_taper_fraction, default 4%).
            # Applied symmetrically to both ends of each trace before CC.
            wlen = int(params.cc.cc_taper_fraction * data.shape[1])
            taper_sides = scipy.signal.windows.hann(2 * wlen + 1)
            taper = np.hstack(
                (taper_sides[:wlen], np.ones(data.shape[1] - 2 * wlen),
                 taper_sides[len(taper_sides) - wlen:]))
            for i in range(data.shape[0]):
                data[i] *= taper
            # index net.sta comps for energy later
            channel_index = {}
            if params.cc.whitening != "N" and params.cc.whitening_type == "PSD":
                psds = []
                for i, name in enumerate(names):
                    n1, s1, l1, c1 = name
                    netsta = "%s.%s.%s" % (n1, s1, l1)
                    if netsta not in channel_index:
                        channel_index[netsta] = {}
                    channel_index[netsta][c1[-1]] = i

                    freqs, pxx = scipy.signal.welch(tmp[i].data,
                                               fs=tmp[i].stats.sampling_rate,
                                               nperseg=nfft,
                                               detrend='constant')
                    psds.append(np.sqrt(pxx))
                psds = np.asarray(psds)
            else:
                psds = np.zeros(1)

            for chan in channel_index:
                comps = channel_index[chan].keys()
                if "E" in comps and "N" in comps:
                    i_e = channel_index[chan]["E"]
                    i_n = channel_index[chan]["N"]
                    # iZ = channel_index[chan]["Z"]
                    mm = psds[[i_e,i_n]].mean(axis=0)
                    psds[i_e] = mm
                    psds[i_n] = mm
                    # psds[iZ] = mm

            # define pairwise CCs
            # Pairwise CC index: for each (sta1, sta2) combination we build
            # a list of (i, j, component_pair) tuples so the inner FFT/IFFT
            # loop can be fully vectorised across all pairs in one batch
            # (see myCorr2 in core/compute.py — replaces the old per-pair loop).

            tmptime = tmp[0].stats.starttime.datetime
            thisdate = tmptime.strftime("%Y-%m-%d")
            tmptime = tmp[0].stats.endtime.datetime
            thistime = tmptime.strftime("%Y-%m-%d %H:%M:%S")

            for filter_name, filter in filters:
                # logger.debug("Processing filter %s" % filter_name)
                # print(filter)

                # Standard operator for CC
                cc_index = []
                if filter.CC:
                    if len(params.cc.components_to_compute):
                        for sta1, sta2 in itertools.combinations(names, 2):
                            n1, s1, l1, c1 = sta1
                            n2, s2, l2, c2 = sta2
                            if n1 == n2 and s1 == s2 and l1 == l2:
                                continue
                            pair = "%s.%s.%s:%s.%s.%s" % (n1, s1, l1, n2, s2, l2)
                            if pair not in pairs:
                                continue
                            comp = "%s%s" % (c1[-1], c2[-1])
                            if comp in params.cc.components_to_compute:
                                cc_index.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                    names.index(sta1), names.index(sta2)])

                # Different iterator func for single station AC or SC:
                single_station_pair_index_sc = []
                single_station_pair_index_ac = []

                if len(params.cc.components_to_compute_single_station):
                    for sta1, sta2 in itertools.combinations_with_replacement(names, 2):
                        n1, s1, l1, c1 = sta1
                        n2, s2, l2, c2 = sta2
                        if n1 != n2 or s1 != s2:
                            continue
                        pair = "%s.%s.%s:%s.%s.%s" % (n1, s1, l1, n2, s2, l2)
                        if pair not in pairs:
                            continue
                        comp = "%s%s" % (c1[-1], c2[-1])
                        if comp in params.cc.components_to_compute_single_station:
                            if filter.AC and c1[-1] == c2[-1]:
                                single_station_pair_index_ac.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                    names.index(sta1), names.index(sta2)])
                            elif filter.SC:
                            # If the components are different, we can just
                            # process them using the default CC code (should warn)
                                single_station_pair_index_sc.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp),
                                    names.index(sta1), names.index(sta2)])
                        if comp[::-1] in params.cc.components_to_compute_single_station:
                            if filter.SC and c1[-1] != c2[-1]:
                                # If the components are different, we can just
                                # process them using the default CC code (should warn)
                                single_station_pair_index_sc.append(
                                    ["%s.%s.%s_%s.%s.%s_%s" % (n1, s1, l1, n2, s2, l2, comp[::-1]),
                                    names.index(sta2), names.index(sta1)])

                # logger.debug("CC index: %s" % cc_index)
                # logger.debug("single station AC: ", single_station_pair_index_ac)
                # logger.debug("single station SC: ", single_station_pair_index_sc)

                filterlow = float(filter.freqmin)
                filterhigh = float(filter.freqmax)
                # logger.debug("Filter: %s - %s Hz" % (filterlow, filterhigh))
                freq_vec = sf.fftfreq(nfft, d=dt)[:nfft // 2]
                freq_sel = np.where((freq_vec >= filterlow) & (freq_vec <= filterhigh))[0]
                low = freq_sel[0] - napod
                if low <= 0:
                    low = 0
                p1 = freq_sel[0]
                p2 = freq_sel[-1]
                high = freq_sel[-1] + napod
                if high > nfft / 2:
                    high = int(nfft // 2)

                # Bandpassed copy of data for this filter band.
                # Computed once here and reused by all modes (AC, CC, SC) and
                # algorithms (CC, PCC) to avoid redundant bandpass calls.
                # For CC+whiten: _data_raw feeds FFT→whiten2 (whiten2 applies
                # the spectral window itself, so no prior bandpass is needed).
                # For PCC and AC: _data_bp is the starting point.
                # For CC+no-whiten: _data_bp is used directly.
                _data_raw = data.copy()
                _data_bp = np.array([
                    bandpass(tr, freqmin=filterlow, freqmax=filterhigh,
                             df=params.cc.cc_sampling_rate, corners=8)
                    for tr in _data_raw
                ])

                # ── Auto-Correlation (AC) ─────────────────────────────────
                if len(single_station_pair_index_ac):
                    tmp = _data_bp.copy()
                    if params.cc.clip_after_whiten:
                        tmp = winsorizing(tmp, params, input="timeseries")

                    if params.cc.cc_type_single_station_AC == "CC":
                        ffts = sf.fftn(tmp, [nfft, ], axes=[1, ])
                        if params.cc.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2,
                            axis=1)))
                        corr = myCorr2(ffts,
                                       np.ceil(params.cc.maxlag / dt),
                                       energy,
                                       single_station_pair_index_ac,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc.cc_normalisation)
                        del energy, ffts

                    elif params.cc.cc_type_single_station_AC == "PCC":
                        tmp_pcc = tmp
                        if params.cc.whitening != "N":
                            ffts_pcc = sf.fftn(tmp, [nfft, ], axes=[1, ])
                            whiten2(ffts_pcc, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                            tmp_pcc = np.real(sf.ifft(ffts_pcc, n=nfft, axis=1)[:, :tmp.shape[1]])
                            del ffts_pcc
                        corr = pcc_xcorr(tmp_pcc,
                                         np.ceil(params.cc.maxlag / dt),
                                         None,
                                         single_station_pair_index_ac,
                                         normalized=params.cc.cc_normalisation)
                    else:
                        logging.error("cc_type_single_station_AC = %s not implemented, "
                              "exiting")
                        exit(1)

                    for key in corr:
                        ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                        logger.debug("CCF ID: %s" % ccfid)
                        if ccfid not in allcorr:
                            allcorr[ccfid] = {}
                        allcorr[ccfid][thistime] = corr[key]
                    del corr, tmp
                    # tmp_pcc is either an alias of tmp or a whitened copy;
                    # both are released above (tmp) or go out of scope here.

                # ── Cross-Correlation (CC) ────────────────────────────────
                if len(cc_index):
                    if params.cc.cc_type == "CC":
                        # CC+whiten: whiten2 applies the spectral window, so
                        # feed raw data into FFT (no prior bandpass needed).
                        # CC+no-whiten: use the pre-bandpassed copy.
                        _cc_src = _data_raw if params.cc.whitening != "N" else _data_bp
                        ffts = sf.fftn(_cc_src, [nfft, ], axes=[1, ])
                        if params.cc.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                        if params.cc.clip_after_whiten:
                            ffts = winsorizing(ffts, params, input="fft", nfft=nfft)
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2, axis=1)))
                        corr = myCorr2(ffts,
                                       np.ceil(params.cc.maxlag / dt),
                                       energy,
                                       cc_index,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc.cc_normalisation)
                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, energy, ffts

                    elif params.cc.cc_type == "PCC":
                        # PCC: bandpass is already done (_data_bp). If whitening
                        # is requested, whiten2 adds spectral shaping on top.
                        _data_pcc = _data_bp.copy()
                        if params.cc.whitening != "N":
                            ffts_pcc = sf.fftn(_data_pcc, [nfft, ], axes=[1, ])
                            whiten2(ffts_pcc, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                            _data_pcc = np.real(sf.ifft(ffts_pcc, n=nfft, axis=1)[:, :_data_pcc.shape[1]])
                            del ffts_pcc
                        corr = pcc_xcorr(_data_pcc,
                                         np.ceil(params.cc.maxlag / dt),
                                         None,
                                         cc_index,
                                         normalized=params.cc.cc_normalisation)
                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, _data_pcc

                    else:
                        logging.error("cc_type = %s not implemented, "
                              "exiting")
                        exit(1)

                # ── Same-station Cross-Component (SC) ─────────────────────
                if len(single_station_pair_index_sc):
                    if params.cc.cc_type_single_station_SC == "CC":
                        # Same logic as CC/CC above: raw → whiten2, or bp → no-whiten.
                        _sc_src = _data_raw if params.cc.whitening != "N" else _data_bp
                        ffts = sf.fftn(_sc_src, [nfft, ], axes=[1, ])
                        if params.cc.whitening != "N":
                            whiten2(ffts, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                        if params.cc.clip_after_whiten:
                            ffts = winsorizing(ffts, params, input="fft", nfft=nfft)
                        energy = np.real(np.sqrt(np.mean(
                            sf.ifft(ffts, n=nfft, axis=1) ** 2, axis=1)))
                        corr = myCorr2(ffts,
                                       np.ceil(params.cc.maxlag / dt),
                                       energy,
                                       single_station_pair_index_sc,
                                       plot=False,
                                       nfft=nfft,
                                       normalized=params.cc.cc_normalisation)
                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            logger.debug("CCF ID - SC: %s" % ccfid)
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, energy, ffts

                    elif params.cc.cc_type_single_station_SC == "PCC":
                        # Same as CC/PCC: start from _data_bp, optionally whiten.
                        _data_sc_pcc = _data_bp.copy()
                        if params.cc.whitening != "N":
                            ffts_sc_pcc = sf.fftn(_data_sc_pcc, [nfft, ], axes=[1, ])
                            whiten2(ffts_sc_pcc, nfft, low, high, p1, p2, psds,
                                    params.cc.whitening_type)  # inplace
                            _data_sc_pcc = np.real(sf.ifft(ffts_sc_pcc, n=nfft, axis=1)[:, :_data_sc_pcc.shape[1]])
                            del ffts_sc_pcc
                        corr = pcc_xcorr(_data_sc_pcc,
                                         np.ceil(params.cc.maxlag / dt),
                                         None,
                                         single_station_pair_index_sc,
                                         normalized=params.cc.cc_normalisation)
                        for key in corr:
                            ccfid = key.replace("_","+") + "+" + filter_name + "+" + thisdate
                            logger.debug("CCF ID - SC PCC: %s" % ccfid)
                            if ccfid not in allcorr:
                                allcorr[ccfid] = {}
                            allcorr[ccfid][thistime] = corr[key]
                        del corr, _data_sc_pcc

                    else:
                        logging.error("cc_type_single_station_SC = %s not implemented, "
                              "exiting")
                        exit(1)
                del _data_raw, _data_bp
            del psds

        if params.cc.keep_all:
            for ccfid, windows in allcorr.items():
                station1, station2, components, filterid, date = ccfid.split('+')
                window_times = list(windows.keys())
                corrs = np.asarray(list(windows.values()))
                xr_save_ccf_all(
                    root=params.global_.output_folder,
                    lineage=lineage_names,
                    step_name=filterid,  # filterid from ccfid = filter step name e.g. 'filter_1'
                    station1=station1,
                    station2=station2,
                    components=components,
                    date=date,
                    window_times=window_times,
                    taxis=t_axis,
                    corrs=corrs,
                )

        if params.cc.keep_days:
            for ccfid in allcorr.keys():
                logging.debug("Exporting %s" % ccfid)
                station1, station2, components, filter_name, date = \
                    ccfid.split('+')

                corrs = np.asarray(list(allcorr[ccfid].values()))
                if not len(corrs):
                    logger.debug("No data to stack.")
                    continue
                corr = stack(corrs, params.cc.stack_method, params.cc.pws_timegate,
                             params.cc.pws_power, params.cc.cc_sampling_rate)
                if not len(corr):
                    logger.debug("No data to save.")
                    continue
                thisdate = goal_day
                thistime = "0_0"
                save_daily_ccf(
                    root=params.global_.output_folder,
                    lineage=lineage_names,
                    step_name=filter_name,  # filter_name from ccfid = filter step name e.g. 'filter_1'
                    station1=station1,
                    station2=station2,
                    components=components,
                    date=thisdate,
                    corr=corr,
                    taxis=t_axis,
                )

        # THIS SHOULD BE IN THE API
        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)

        logger.info("Job Finished. It took %.2f seconds (preprocess: %.2f s & "
                     "process %.2f s)" % ((time.time() - jt),
                                          start_processing - jt,
                                          time.time() - start_processing))
        del stream, allcorr
    logger.info('*** Finished: Compute CC ***')
