import time

import numpy as np
import xarray as xr
from obspy.signal.regression import linear_regression
from ...db import connect, get_logger
from ...stations import get_interstation_distance, get_station_pairs
from ...workflow import (extend_days, get_next_lineage_batch, is_next_job_for_step, massive_update_job)
from ...io import xr_get_mwcs, xr_save_dtt


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

        mov_stacks = params.mov_stack


        netsta1, netsta2 = pair.split(':')
        station1, station2 = pair.split(":")
        if station1 == station2:
            components_to_compute = params.components_to_compute_single_station
        else:
            components_to_compute = params.components_to_compute
            
        n1, s1, l1 = netsta1.split(".")
        n2, s2, l2 = netsta2.split(".")
        # todo, include location code for computing distances?
        dpair = "%s_%s_%s_%s" % (n1, s1, n2, s2)
        dist = interstations[dpair] if dpair in interstations else 0.0
        logger.info(
            "There are DTT jobs for some days to recompute for %s" % pair)
        for components in components_to_compute:
            for mov_stack in mov_stacks:
                try:
                    ds = xr_get_mwcs(params.output_folder, lineage_names,
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

                # ── Row-by-row WLS regression ────────────────────────────
                m_vals   = []
                em_vals  = []
                a_vals   = []
                ea_vals  = []
                m0_vals  = []
                em0_vals = []
                mcoh_vals = []   # mean coherence over valid lag window
                out_times = []

                for i in range(len(times)):
                    errArray = EM[i]
                    dtArray  = M[i]
                    cohArray = MCOH[i]

                    index = np.where((errArray != 1.0) & (cohArray != 0.0))[0]
                    if len(index) < 2:
                        continue

                    w = 1.0 / errArray[index]
                    w[~np.isfinite(w)] = 1.0
                    VecXfilt = tArray[index]
                    VecYfilt = dtArray[index]

                    m, a, em, ea = linear_regression(
                        VecXfilt, VecYfilt, w, intercept_origin=False)
                    m0, em0 = linear_regression(
                        VecXfilt, VecYfilt, w, intercept_origin=True)

                    # Mean coherence over the valid lag window (tindex)
                    mcoh = float(np.nanmean(cohArray[tindex])) if len(tindex) else 0.0

                    m_vals.append(m)
                    em_vals.append(em)
                    a_vals.append(a)
                    ea_vals.append(ea)
                    m0_vals.append(m0)
                    em0_vals.append(em0)
                    mcoh_vals.append(mcoh)
                    out_times.append(times[i])

                if not out_times:
                    continue

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
                xr_save_dtt(params.output_folder, lineage_names, step.step_name,
                            station1, station2, components, mov_stack, ds_out)

        massive_update_job(db, jobs, "D")

    logger.info('*** Finished: Compute DTT ***')
