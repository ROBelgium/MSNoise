"""Compute dt/t from MWCS measurements using weighted linear regression.

Reads the MWCS output (delay, error, coherence per lag window per time step),
applies lag-window masking and quality thresholds, then fits a slope
(dv/v = -dt/t) using vectorized WLS — both origin-forced and with intercept.

Configuration Parameters
------------------------

* |mwcs_dtt.dtt_minlag|
* |mwcs_dtt.dtt_width|
* |mwcs_dtt.dtt_lag|
* |mwcs_dtt.dtt_v|
* |mwcs_dtt.dtt_sides|
* |mwcs_dtt.dtt_mincoh|
* |mwcs_dtt.dtt_maxerr|
* |mwcs_dtt.dtt_maxdtt|
* |stack.mov_stack|
* |cc.components_to_compute|
* |cc.components_to_compute_single_station|
* |global.hpc|
"""

import time

import numpy as np
import xarray as xr
from .core.db import connect, get_logger
from .core.stations import get_interstation_distance, get_station_pairs
from .core.workflow import (extend_days, get_next_lineage_batch, is_next_job_for_step, massive_update_job, propagate_downstream)
from .core.io import xr_get_mwcs, xr_save_dtt


def main(loglevel="INFO"):
    logger = get_logger('msnoise.mwcs_dtt', loglevel, with_pid=True)
    logger.info('*** Starting: Compute DT/T ***')
    db = connect()

    # start, end, datelist = build_movstack_datelist(db)

    interstations = {}
    for sta1, sta2 in get_station_pairs(db):
        s1 = "%s_%s" % (sta1.net, sta1.sta)
        s2 = "%s_%s" % (sta2.net, sta2.sta)
        if s1 == s2:
            interstations["%s_%s" % (s1, s2)] = 0.0
        else:
            interstations["%s_%s" % (s1, s2)] = get_interstation_distance(sta1,
                                                                          sta2,
                                                                          sta1.coordinates)

    while is_next_job_for_step(db, step_category="mwcs_dtt"):
        logger.debug("Getting the next batch")
        batch = get_next_lineage_batch(db, step_category="mwcs_dtt", group_by="pair_lineage", loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        pair = batch["pair"]
        days = batch["days"]
        params = batch["params"]
        lineage_names = batch["lineage_names_upstream"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]

        logger.info(f"New MWCS-DTT Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        mov_stacks = params.stack.mov_stack


        netsta1, netsta2 = pair.split(':')
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
            
        n1, s1, l1 = netsta1.split(".")
        n2, s2, l2 = netsta2.split(".")
        # todo, include location code for computing distances?
        dpair = "%s_%s_%s_%s" % (n1, s1, n2, s2)
        dist = interstations[dpair] if dpair in interstations else 0.0
        logger.info(
            "There are DTT jobs for some days to recompute for %s" % pair)
        for components in components_to_compute:
            # Skip component types not enabled in filter config
            if station1 == station2:
                _is_ac = len(components) >= 2 and components[0] == components[-1]
                if _is_ac and not params.filter.AC:
                    continue
                if not _is_ac and not params.filter.SC:
                    continue
            for mov_stack in mov_stacks:
                try:
                    ds = xr_get_mwcs(params.global_.output_folder, lineage_names,
                                     station1, station2, components, mov_stack)
                except FileNotFoundError as fullpath:
                    logger.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
                    continue

                # ── Time filtering ───────────────────────────────────────
                to_search = extend_days(days)
                times_floor = ds.coords["times"].values.astype("datetime64[D]")
                mask = np.isin(times_floor, np.array(list(to_search), dtype="datetime64[D]"))
                ds = ds.isel(times=mask)
                # Drop times where all MWCS values are NaN
                ds = ds.dropna("times", how="all")
                if ds.sizes["times"] == 0:
                    continue

                # ── Extract per-key arrays: shape (times, taxis) ─────────
                da = ds["MWCS"]  # dims: (times, taxis, keys)
                M    = da.sel(keys="M").values.copy()     # (times, taxis)
                EM   = da.sel(keys="EM").values.copy()
                MCOH = da.sel(keys="MCOH").values.copy()
                tArray = da.coords["taxis"].values        # (taxis,)
                times  = ds.coords["times"].values
                ds.close()  # release lazy file handle — all needed data is now in numpy

                # ── Lag window masking ───────────────────────────────────
                if params.mwcs_dtt.dtt_lag == "static":
                    lmlag = -params.mwcs_dtt.dtt_minlag
                    rmlag = params.mwcs_dtt.dtt_minlag
                else:
                    lmlag = -dist / params.mwcs_dtt.dtt_v
                    rmlag = dist / params.mwcs_dtt.dtt_v
                lMlag = lmlag - params.mwcs_dtt.dtt_width
                rMlag = rmlag + params.mwcs_dtt.dtt_width

                if params.mwcs_dtt.dtt_sides == "both":
                    tindex = np.where(
                        ((tArray >= lMlag) & (tArray <= lmlag)) |
                        ((tArray >= rmlag) & (tArray <= rMlag)))[0]
                elif params.mwcs_dtt.dtt_sides == "left":
                    tindex = np.where((tArray >= lMlag) & (tArray <= lmlag))[0]
                else:
                    tindex = np.where((tArray >= rmlag) & (tArray <= rMlag))[0]

                tmp = np.setdiff1d(np.arange(len(tArray)), tindex)
                EM[:, tmp]   = 1.0
                MCOH[:, tmp] = 0.0

                MCOH[MCOH < params.mwcs_dtt.dtt_mincoh] = 0.0
                if params.mwcs_dtt.dtt_maxerr > 0:
                    EM[EM > params.mwcs_dtt.dtt_maxerr] = 1.0

                # Exclude values exceeding dtt_maxdtt
                dtt_values = np.abs(M / tArray)
                row_indices, col_indices = np.where(dtt_values > params.mwcs_dtt.dtt_maxdtt)
                EM[row_indices, col_indices]   = 1.0
                MCOH[row_indices, col_indices] = 0.0

                # ── Vectorized WLS regression (all times at once) ────────
                # ObsPy linear_regression minimises sum(w² * residuals²) where
                # w = weights = 1/sigma.  Closed-form matching that exactly:
                #   origin-forced:  m0 = sum(w²·t·y) / sum(w²·t²)
                #   with intercept: solved via 2×2 normal equations
                # Valid mask: EM != 1.0 (not flagged) AND MCOH != 0.0 (coherent)
                valid_mask = (EM != 1.0) & (MCOH != 0.0)   # (n_times, n_lag)

                # Weight matrix: w = 1/EM where valid, 0 elsewhere
                w_raw = np.where(valid_mask, 1.0 / np.where(EM != 0, EM, 1.0), 0.0)
                w_raw = np.where(np.isfinite(w_raw), w_raw, 1.0)
                # Zero out invalid entries so they do not contribute to sums
                w_raw = w_raw * valid_mask

                # Require at least 2 valid lag samples per time step
                n_valid = valid_mask.sum(axis=1)             # (n_times,)
                good = n_valid >= 2

                t = tArray                                    # (n_lag,)
                y = M                                         # (n_times, n_lag)
                w2 = w_raw ** 2                              # (n_times, n_lag)

                # ── Origin-forced: m0 = sum(w²·t·y) / sum(w²·t²) ──────
                Sw2t2 = (w2 * t[None, :] ** 2).sum(axis=1)  # (n_times,)
                Sw2ty = (w2 * t[None, :] * y).sum(axis=1)
                m0_all = np.where(Sw2t2 != 0, Sw2ty / Sw2t2, np.nan)
                # Error: sqrt(chi2 / ((n-1) * Sw2t2))
                res0   = y - m0_all[:, None] * t[None, :]
                chi20  = (w2 * res0 ** 2).sum(axis=1)
                em0_all = np.where(
                    (Sw2t2 != 0) & (n_valid > 1),
                    np.sqrt(np.abs(chi20) / np.maximum(n_valid - 1, 1) / np.maximum(Sw2t2, 1e-30)),
                    np.nan,
                )

                # ── With intercept: 2×2 normal equations ─────────────
                Sw2   = w2.sum(axis=1)
                Sw2t  = (w2 * t[None, :]).sum(axis=1)
                Sw2y  = (w2 * y).sum(axis=1)
                denom = Sw2 * Sw2t2 - Sw2t ** 2
                m_all  = np.where(denom != 0, (Sw2 * Sw2ty - Sw2t * Sw2y) / denom, np.nan)
                a_all  = np.where(denom != 0, (Sw2t2 * Sw2y - Sw2t * Sw2ty) / denom, np.nan)
                res1   = y - m_all[:, None] * t[None, :] - a_all[:, None]
                chi21  = (w2 * res1 ** 2).sum(axis=1)
                scale  = np.where(
                    (denom != 0) & (n_valid > 2),
                    np.abs(chi21) / np.maximum(n_valid - 2, 1),
                    np.nan,
                )
                em_all  = np.where(denom != 0, np.sqrt(np.abs(scale * Sw2t2 / np.maximum(denom, 1e-30))), np.nan)
                ea_all  = np.where(denom != 0, np.sqrt(np.abs(scale * Sw2   / np.maximum(denom, 1e-30))), np.nan)

                # ── Mean coherence over valid lag window ──────────────
                coh_sum = MCOH[:, tindex].sum(axis=1) if len(tindex) else np.zeros(len(times))
                coh_cnt = np.maximum((MCOH[:, tindex] != 0).sum(axis=1), 1) if len(tindex) else np.ones(len(times))
                mcoh_all = coh_sum / coh_cnt

                # ── Collect results for valid time steps ──────────────
                out_mask   = good
                out_times  = times[out_mask]
                m_vals     = m_all[out_mask]
                em_vals    = em_all[out_mask]
                a_vals     = a_all[out_mask]
                ea_vals    = ea_all[out_mask]
                m0_vals    = m0_all[out_mask]
                em0_vals   = em0_all[out_mask]
                mcoh_vals  = mcoh_all[out_mask]

                if len(out_times) == 0:
                    # Free all intermediates before continuing
                    del M, EM, MCOH, tArray, times, dtt_values, tmp, tindex
                    del valid_mask, w_raw, w2, n_valid, good, t, y
                    del Sw2t2, Sw2ty, m0_all, res0, chi20, em0_all
                    del Sw2, Sw2t, Sw2y, denom, m_all, a_all
                    del res1, chi21, scale, em_all, ea_all
                    del coh_sum, coh_cnt, mcoh_all
                    del out_mask, out_times
                    continue

                # Free large intermediates — only the _vals arrays are still needed
                del M, EM, MCOH, tArray, times, dtt_values, tmp, tindex
                del valid_mask, w_raw, w2, n_valid, good, t, y
                del Sw2t2, Sw2ty, m0_all, res0, chi20, em0_all
                del Sw2, Sw2t, Sw2y, denom, m_all, a_all
                del res1, chi21, scale, em_all, ea_all
                del coh_sum, coh_cnt, mcoh_all, out_mask

                # ── Build output Dataset ─────────────────────────────────
                out_times = np.array(out_times)
                ds_out = xr.Dataset(
                    {
                        "DTT": xr.DataArray(
                            np.column_stack([m_vals, em_vals, a_vals,
                                            ea_vals, m0_vals, em0_vals,
                                            mcoh_vals]),
                            dims=["times", "keys"],
                            coords={
                                "times": out_times,
                                "keys":  ["m", "em", "a", "ea", "m0", "em0",
                                          "mcoh"],
                            },
                        )
                    }
                )
                xr_save_dtt(params.global_.output_folder, lineage_names, step.step_name,
                            station1, station2, components, mov_stack, ds_out)
                del m_vals, em_vals, a_vals, ea_vals, m0_vals, em0_vals, mcoh_vals
                del out_times, ds_out

        massive_update_job(db, jobs, "D")
        if not batch["params"].global_.hpc:
            propagate_downstream(db, batch)

    logger.info('*** Finished: Compute DTT ***')
