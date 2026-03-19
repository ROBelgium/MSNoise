"""
TODO
"""

import argparse

import scipy.signal

from .api import *
from .wiener import *
from obspy.core import AttribDict

import logbook
import matplotlib.pyplot as plt


def main(stype, loglevel="INFO"):
    """Computes the REF/MOV stacks.

    Parameters
    ----------
    stype : {'mov', 'ref'}
        Defines which of the REF or Moving-window stacks must be exported
    interval : float, optional
        Number of days before now to search for modified CC jobs

    """
    # logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.stack_child', loglevel, with_pid=True)
    logger.debug('Starting the %s stack' % stype)
    db = connect()

    orig_params = get_params(db)
    taxis = get_t_axis(db)

    while is_next_job_for_step(db, step_category="stack_ref"):
        logger.info("Getting the next job")
        jobs, step = get_next_job_for_step(db, step_category="stack_ref", group_by="pair_lineage")

        if not len(jobs):
            # edge case, should only occur when is_next returns true, but
            # get_next receives no jobs (heavily parallelised calls).
            time.sleep(np.random.random())
            continue
        pair = jobs[0].pair
        lineage_str = jobs[0].lineage
        if not lineage_str:
            raise ValueError("STACK jobs must have a non-empty lineage (v2 assumption)")

        refs, days = zip(*[[job.ref, job.day] for job in jobs])

        # 1) current step config (stack_*)
        step_params = get_config_set_details(
            db,
            jobs[0].config_category,
            jobs[0].config_set_number,
            format="AttribDict"
        )

        # 2) resolve lineage steps from the job, and merge params for THIS lineage only
        lineage_steps = lineage_str_to_steps(db, lineage_str, strict=True)

        lineage_steps, lineage_names, params = get_merged_params_for_lineage(
            db, orig_params, step_params, lineage_steps
        )

        # drop the "current step name" if you use lineage_names as folder segments upstream of stack output
        lineage_names = lineage_names[:-1]

        logger.info(f"New STACK Job: pair={pair} n_days={len(days)} lineage={lineage_str}")

        mov_stacks = params.mov_stack
        wiener_mlen = params.wiener_mlen
        wiener_nlen = params.wiener_nlen
        wienerfilt = params.wienerfilt
        wiener_M = int(pd.to_timedelta(wiener_mlen).total_seconds() / params.corr_duration)
        wiener_N = int(pd.to_timedelta(wiener_nlen).total_seconds() * params.cc_sampling_rate)

        # is there a better alternative for threshold?
        if params.keep_all:
            wiener_gap_threshold = wiener_M  # no. indices which will be considered adjacent by wiener
        else:
            wiener_gap_threshold = pd.to_timedelta(wiener_mlen).days

        if wienerfilt:
            logger.info('Wiener filter enabled, will apply to CCFs before stacking')

        logger.info(
            "There are STACKS jobs for some days to recompute for %s" % pair)

        sta1, sta2 = pair.split(':')
        filterid = 1
        for components in params.components_to_compute:
            logger.info('Processing %s-%s-%i REF stack' %
                (pair, components, filterid))

            start, end, datelist = build_ref_datelist(db)

            if params.keep_all:
                c = get_results_all(db, params.output_folder, lineage_names,
                                    sta1, sta2, components, datelist, format="xarray")
                if not len(c):
                    logger.warning("No data found for %s-%s-%i" % (sta1, sta2, filterid))
                    continue

            else:
                c = get_results(db, sta1, sta2, filterid, components, datelist, mov_stack=1, format="xarray",
                                params=params)

            is_valid, message = validate_stack_data(c, "reference")
            if not is_valid:
                logger.error(f"Invalid reference data for {sta1}:{sta2}-{components}-{filterid}: {message}")
                continue
            elif "Warning" in message:
                logger.warning(f"{sta1}:{sta2}-{components}-{filterid}: {message}")

            # dr = xr_save_ccf(sta1, sta2, components, filterid, 1, taxis, c)
            dr = c
            if not c.data_vars:
                logger.debug('No data found for %s:%s-%s-%i' %
                             (sta1, sta2, components, filterid))
                continue

            if wienerfilt:
                dr = wiener_filt(dr, wiener_M, wiener_N, wiener_gap_threshold)

            # TODO add other stack methods here! using apply?
            _ = dr.mean(dim="times")
            xr_save_ref(params.output_folder, lineage_names, step.step_name,
                        sta1, sta2, components, filterid, taxis, _)


        massive_update_job(db, jobs, "D")

        # if stype != "step" and not params.hpc:
        #     for job in jobs:
        #         update_job(db, job.day, job.pair, 'MWCS', 'T')
        #         update_job(db, job.day, job.pair, 'WCT', 'T')
        #         update_job(db, job.day, job.pair, 'STR', 'T')
