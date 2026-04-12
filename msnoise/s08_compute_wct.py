"""
Wavelet Coherence Transform (WCT) Computation
===============================================

This script performs the computation of the Wavelet Coherence Transform (WCT), a tool used to analyze the correlation between two time series in the time-frequency domain. The script supports parallel processing and interacts with a database to manage job statuses.

Filter Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* |wavelet.wct_freqmin|
* |wavelet.wct_freqmax|
* |wavelet.wct_ns|
* |wavelet.wct_nt|
* |wavelet.wct_vpo|
* |wavelet.wct_nptsfreq|
* |wavelet.wct_norm|
* |wavelet.wavelet_type|
* |wavelet.wct_compute_dtt|
* |stack.mov_stack|
* |refstack.ref_begin|
* |refstack.ref_end|
* |cc.cc_sampling_rate|
* |cc.components_to_compute|
* |cc.components_to_compute_single_station|
* |global.hpc|

This process is job-based, so it is possible to run several instances in
parallel.

To run this step:

.. code-block:: sh

    $ msnoise cc dtt compute_wct

This step also supports parallel processing/threading:

.. code-block:: sh

    $ msnoise -t 4 cc dtt compute_wct

will start 4 instances of the code (after 1 second delay to avoid database
conflicts). This works both with SQLite and MySQL but be aware problems
could occur with SQLite.

"""

import time
import numpy as np
import xarray as xr
from .core.db import connect, get_logger
from .core.workflow import (compute_rolling_ref, extend_days, get_next_lineage_batch,
                             get_step_successors, get_t_axis, is_next_job_for_step,
                             massive_update_job, propagate_downstream,
                             refstack_is_rolling, resolve_lineage_params)
from .core.signal import xwt, prepare_ref_wct, apply_wct, get_wavelet_type
from .core.compute import compute_wct_dtt_batch, build_wct_dtt_dataset
from .core.io import xr_get_ccf, xr_get_ref, xr_save_wct, xr_save_wct_dtt

def main(loglevel="INFO"):
    """Compute Wavelet Coherence Transform (WCT) and optionally dt/t inline.

    **Frequency sampling** (``wct_nptsfreq`` and ``wct_vpo`` config parameters)

    The CWT samples ``nptsfreq`` scales linearly between ``wct_freqmin`` and
    ``wct_freqmax``.  Because the Morlet wavelet only resolves ~1/vpo octaves
    per scale, adjacent bins are correlated — the number of *independent* scales
    is approximately ``vpo × log2(freqmax / freqmin)``.  For the default 0.1–1.0
    Hz band at ``vpo=12`` that is ~40 independent scales; the previous default
    of ``nptsfreq=300`` at ``vpo=20`` computed >7× more scales than the wavelet
    can independently resolve, tripling CWT cost with no scientific benefit.
    Default values have been revised: ``nptsfreq=100``, ``vpo=12``.

    **Fused mode** (``wct_compute_dtt=True``, the default)

    In fused mode all downstream ``wavelet_dtt_N`` config sets are resolved at
    the start of each batch.  For each day, immediately after calling
    :func:`~msnoise.core.signal.apply_wct`, the WCT arrays (``WXamp``,
    ``Wcoh``, ``WXdt``) are passed directly to
    :func:`~msnoise.core.signal.compute_wct_dtt` and
    :func:`~msnoise.core.signal.get_wct_avgcoh` for every downstream config.
    Only the compact ``DTT/ERR/COH`` dataset (per frequency, per day) is written
    to disk — the full 3-D WCT arrays are **never stored**.

    Storage comparison (1 pair × 1 component × 1 year, default params):

    * WCT file (``WXamp + Wcoh + WXdt``, 365 × 100 × 4801 samples):
      ~700 MB raw, ~280 MB compressed (float32 + zlib-4).
    * DTT file (``DTT + ERR + COH``, 365 × freq_subset):
      ~150 KB compressed — roughly **2000× smaller**.

    For a network with 50 pairs × 3 components × 3 mov_stacks the saving is
    ~125 GB of intermediate WCT storage per year.

    Multiple downstream ``wavelet_dtt`` config sets (e.g. ``wavelet_dtt_1``
    with 0.1–0.5 Hz and ``wavelet_dtt_2`` with 0.5–1.0 Hz) are all computed
    in a single pass through the data — each WCT array slice is used by all
    configs before being discarded.

    **Standard mode** (``wct_compute_dtt=False``)

    The full WCT arrays are stored in NetCDF files under the ``wavelet`` step
    output folder.  Use this when you want to re-run ``wavelet_dtt`` with
    different parameters (freqmin, freqmax, mincoh, coda_cycles, etc.) without
    recomputing the CWT, or when you need the 2-D WCT images for inspection.
    """
    global logger
    logger = get_logger('msnoise.wavelet', loglevel, with_pid=True)
    logger.info('*** Starting: Compute WCT ***')

    db = connect()

    while is_next_job_for_step(db, step_category="wavelet"):
        batch = get_next_lineage_batch(db, step_category="wavelet", group_by="pair_lineage",
                                       loglevel=loglevel)

        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        lineage_names = batch["lineage_names_upstream"]
        lineage_names_mov = batch["lineage_names_mov"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]
        root = params.global_.output_folder

        logger.info(f"New WCT Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        taxis = get_t_axis(params)

        station1, station2 = pair.split(":")

        # Respect filter CC/SC/AC flags — skip pair types never computed
        if station1 == station2:
            if not (params.filter.SC or params.filter.AC):
                continue
            components_to_compute = params.cc.components_to_compute_single_station
        else:
            if not params.filter.CC:
                continue
            components_to_compute = params.cc.components_to_compute

        mov_stacks = params.stack.mov_stack
        goal_sampling_rate = params.cc.cc_sampling_rate

        # Filter frequency range from the wavelet step config
        freqmin = params.wavelet.wct_freqmin
        freqmax = params.wavelet.wct_freqmax

        # Wavelet computation parameters from the wavelet step config
        ns = params.wavelet.wct_ns
        nt = params.wavelet.wct_nt
        vpo = params.wavelet.wct_vpo
        nptsfreq = params.wavelet.wct_nptsfreq
        wct_norm = params.wavelet.wct_norm
        _wt_raw = (params.wavelet.wavelet_type or "").strip()
        try:
            wavelet_type = tuple(eval(_wt_raw)) if _wt_raw else ("Morlet", 6.)
        except Exception:
            logger.warning(f"Could not parse wavelet_type {_wt_raw!r}, using Morlet(6)")
            wavelet_type = ("Morlet", 6.)

        logger.info(f"WCT params: freqmin={freqmin} freqmax={freqmax} ns={ns} nt={nt} vpo={vpo}")

        # ── Fused mode: resolve all downstream wavelet_dtt configs ───────────
        compute_dtt_inline = bool(getattr(params.wavelet, "wct_compute_dtt", True))
        dtt_configs = []  # list of (dtt_step, dtt_params) for all wavelet_dtt_N
        if compute_dtt_inline:
            for dtt_step in get_step_successors(db, step.step_id):
                if dtt_step.category != "wavelet_dtt":
                    continue
                dtt_lineage = batch["lineage_names"] + [dtt_step.step_name]
                try:
                    _, _, dtt_params = resolve_lineage_params(db, dtt_lineage)
                    dtt_configs.append((dtt_step, dtt_params, dtt_lineage))
                except Exception as e:
                    logger.warning(
                        f"Could not resolve params for {dtt_step.step_name}: {e}; "
                        "will fall back to storing WCT files."
                    )
            if dtt_configs:
                logger.info(
                    f"Fused mode: computing DTT inline for "
                    f"{[d[0].step_name for d in dtt_configs]}"
                )
            else:
                logger.warning(
                    "wct_compute_dtt=True but no wavelet_dtt steps found downstream; "
                    "falling back to storing WCT files."
                )
                compute_dtt_inline = False

        # Build the mother wavelet object once — same for all components/stacks/days
        mother = get_wavelet_type(wavelet_type)

        # Pre-compute interstation distance for fused dynamic-lag DTT
        _dist = 0.0
        if compute_dtt_inline and dtt_configs:
            _needs_dist = any(
                str(d[1].wavelet_dtt.wct_lag or "static") == "dynamic"
                for d in dtt_configs
            )
            if _needs_dist:
                try:
                    from .core.stations import get_station, get_interstation_distance
                    n1, s1_, l1 = station1.split(".")
                    n2, s2_, l2 = station2.split(".")
                    _sta1 = get_station(db, n1, s1_)
                    _sta2 = get_station(db, n2, s2_)
                    if _sta1 and _sta2 and station1 != station2:
                        _dist = get_interstation_distance(
                            _sta1, _sta2, _sta1.coordinates)
                except Exception as _e:
                    logger.debug(f"Could not compute interstation distance: {_e}")

        for component in components_to_compute:
            # Skip component types not enabled in filter config
            if station1 == station2:
                _is_ac = len(component) >= 2 and component[0] == component[-1]
                if _is_ac and not params.filter.AC:
                    continue
                if not _is_ac and not params.filter.SC:
                    continue
            rolling_mode = refstack_is_rolling(params)
            ref_wct_data = None  # pre-computed ref CWT for Mode A (set below)

            if not rolling_mode:
                # Mode A: load fixed REF from disk
                try:
                    ref_da = xr_get_ref(root, lineage_names, station1, station2,
                                        component, taxis, ignore_network=True)
                    ref = ref_da.values
                    if wct_norm:
                        ori_waveform = ref / ref.max()
                    else:
                        ori_waveform = ref
                except FileNotFoundError as fp:
                    logger.error(f"FILE DOES NOT EXIST: {fp}, skipping")
                    continue
                except Exception as e:
                    logger.error(f"Error getting reference waveform: {str(e)}")
                    continue

                if not len(ref):
                    continue

                # Pre-compute CWT + smoothed power of the reference ONCE.
                # In Mode A the reference never changes, so we avoid repeating
                # this O(N·S) computation for every day in the time loop.
                try:
                    ref_wct_data = prepare_ref_wct(
                        ori_waveform, goal_sampling_rate,
                        int(ns), float(nt), int(vpo),
                        freqmin, freqmax, int(nptsfreq), mother
                    )
                    # Extract freqs/coi for Dataset construction later
                    _, _, _, freqs, coi, _, _, _ = ref_wct_data
                except Exception as e:
                    logger.error(f"Error preparing ref WCT: {str(e)}, falling back to xwt()")
                    ref_wct_data = None
            # Mode B: ori_waveform and ref_wct_data set per time-step below

            for mov_stack in mov_stacks:
                WXamp_list = []
                WXcoh_list = []
                WXdt_list = []
                dates_list = []
                # Fused mode: per-wavelet_dtt accumulation — {dtt_step_name: {rows}}
                dtt_accum = {
                    d[0].step_name: {"dtt": [], "err": [], "coh": [], "freqs_subset": None}
                    for d in dtt_configs
                } if compute_dtt_inline else {}

                try:
                    data = xr_get_ccf(root, lineage_names_mov, station1, station2,
                                      component, mov_stack, taxis)
                except FileNotFoundError as fp:
                    logger.error(f"FILE DOES NOT EXIST: {fp}, skipping")
                    continue

                # ── Time filtering ───────────────────────────────────────
                to_search = extend_days(days)
                times_floor = data.coords["times"].values.astype("datetime64[D]")
                mask = np.isin(times_floor, np.array(list(to_search), dtype="datetime64[D]"))
                data = data.isel(times=mask)
                data = data.dropna("times", how="all")

                if rolling_mode:
                    ref_rolling = compute_rolling_ref(
                        data, int(params.refstack.ref_begin), int(params.refstack.ref_end)
                    )

                # Materialise the full CCF array once before the time loop —
                # avoids a disk read per iteration now that data is lazy.
                data_np = data.values  # shape: (times, taxis)

                for _i_row, date in enumerate(data.coords["times"].values):
                    waveform = data_np[_i_row]
                    if wct_norm:
                        new_waveform = waveform / waveform.max()
                    else:
                        new_waveform = waveform

                    if rolling_mode:
                        _rr = ref_rolling[_i_row]
                        ori_waveform = _rr / _rr.max() if wct_norm else _rr
                        # Mode B: no pre-computed ref (changes every step)
                        cur_ref_data = None
                    else:
                        cur_ref_data = ref_wct_data

                    try:
                        if cur_ref_data is not None:
                            # Mode A fast path: ref CWT already computed
                            WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = apply_wct(
                                cur_ref_data, new_waveform, int(ns), float(nt)
                            )
                        else:
                            # Mode B or fallback: full xwt each call
                            WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = xwt(
                                ori_waveform, new_waveform, goal_sampling_rate,
                                int(ns), float(nt), int(vpo),
                                freqmin, freqmax, int(nptsfreq), wavelet_type
                            )
                        if compute_dtt_inline and dtt_configs:
                            # Fused: compute DTT for each downstream config now,
                            # while WXamp/Wcoh/WXdt are hot in memory
                            for dtt_step, dtt_params, _dtt_lin in dtt_configs:
                                sname = dtt_step.step_name
                                acc   = dtt_accum[sname]
                                try:
                                    dtt_row, err_row, coh_row, freqs_sub = compute_wct_dtt_batch(
                                        freqs, taxis, WXamp, Wcoh, WXdt,
                                        dtt_params, dist=_dist,
                                    )
                                    acc["dtt"].append(dtt_row)
                                    acc["err"].append(err_row)
                                    acc["coh"].append(coh_row)
                                    if acc["freqs_subset"] is None:
                                        acc["freqs_subset"] = freqs_sub
                                except Exception as e_dtt:
                                    logger.error(
                                        f"Fused DTT error for {sname}/{date}: {e_dtt}"
                                    )
                            dates_list.append(date)
                        else:
                            # Standard mode: accumulate WCT arrays for file storage
                            WXamp_list.append(WXamp)
                            WXcoh_list.append(Wcoh)
                            WXdt_list.append(WXdt)
                            dates_list.append(date)
                    except Exception as e:
                        logger.error(f"Error in WCT for {date}: {str(e)}")
                        continue

                if not dates_list:
                    continue

                if compute_dtt_inline and dtt_configs:
                    # ── Fused mode: save one WCT-DTT file per wavelet_dtt config ──
                    for dtt_step, _dtt_params, dtt_lineage in dtt_configs:
                        sname = dtt_step.step_name
                        acc   = dtt_accum[sname]
                        if not acc["dtt"]:
                            continue
                        freqs_sub = acc["freqs_subset"]
                        if freqs_sub is None or len(freqs_sub) == 0:
                            logger.warning(
                                f"No frequencies in DTT band for {sname}, skipping"
                            )
                            continue
                        ds_dtt = build_wct_dtt_dataset(
                            dates_list, acc["dtt"], acc["err"], acc["coh"],
                            acc["freqs_subset"],
                        )
                        try:
                            # Write to the wavelet_dtt step's own output path,
                            # using the upstream lineage (without wavelet_dtt_N)
                            # so the path matches what s09 / wavelet_dtt_dvv expects.
                            xr_save_wct_dtt(
                                root, batch["lineage_names"], sname,
                                station1, station2, component, mov_stack,
                                taxis, ds_dtt,
                            )
                            logger.debug(
                                f"Fused DTT saved: {sname}/{station1}:{station2}"
                                f"/{component}/{mov_stack}"
                            )
                        except Exception as e:
                            logger.error(
                                f"Error saving fused DTT for {sname}: {e}"
                            )
                else:
                    # ── Standard mode: store full WCT arrays ─────────────────────
                    try:
                        wct_ds = xr.Dataset({
                            "WXamp": xr.DataArray(
                                np.array(WXamp_list).real,
                                dims=["times", "freqs", "taxis"],
                                coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
                            ),
                            "Wcoh": xr.DataArray(
                                np.array(WXcoh_list).real,
                                dims=["times", "freqs", "taxis"],
                                coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
                            ),
                            "WXdt": xr.DataArray(
                                np.array(WXdt_list).real,
                                dims=["times", "freqs", "taxis"],
                                coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
                            ),
                        })
                        xr_save_wct(root, lineage_names, step.step_name,
                                    station1, station2, component, mov_stack,
                                    wct_ds)
                    except Exception as e:
                        logger.error(f"Error saving WCT: {str(e)}")

        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            # In fused mode, wavelet_dtt output files are already written above.
            # propagate_downstream still creates wavelet_dtt jobs — s09 will run
            # and find its output already present (harmless, idempotent via
            # _xr_insert_or_update), then propagate to wavelet_dtt_dvv.
            propagate_downstream(db, batch)

    db.close()
    logger.info('*** Finished: Compute WCT ***')
