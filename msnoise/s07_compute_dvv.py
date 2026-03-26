import traceback
import time
import numpy as np

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.signal.regression import linear_regression
from .api import (
    compute_dvv,
    connect,
    get_logger,
    get_next_lineage_batch,
    is_next_job_for_step,
    massive_update_job,
    xr_save_dvv,
)


def main(loglevel="INFO"):
    logger = get_logger('msnoise.stretching', loglevel, with_pid=True)
    logger.info('*** Starting: Compute DV/V ***')
    db = connect()

    while is_next_job_for_step(db, step_category="stretching"):
        logger.info("Getting the next job")
        batch = get_next_lineage_batch(db, step_category="stretching", group_by="pair_lineage",
                                       loglevel=loglevel)
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs = batch["jobs"]
        params = batch["params"]
        lineage_names = batch["lineage_names_upstream"]
        lineage_str = batch["lineage_str"]
        step = batch["step"]

        root = params.output_folder
        mov_stacks = params.mov_stack

        filt_all_components = np.unique(
            params.components_to_compute + params.components_to_compute_single_station
        )

        logger.info(f"New DVV Job: lineage={lineage_str}")

        for mov_stack in mov_stacks:
            for components in filt_all_components:
                logger.debug("Processing m%s %s" % (mov_stack, components))
                try:
                    dvv = compute_dvv(db, root, lineage_names, mov_stack,
                                       pairs=None, components=components, params=params)
                except ValueError:
                    traceback.print_exc()
                    logger.error("No data for m%s: %s" % (mov_stack, components))
                    continue
                xr_save_dvv(root, lineage_names, step.step_name, components, mov_stack, dvv)
                del dvv
            try:
                dvv = compute_dvv(db, root, lineage_names, mov_stack,
                                   pairs=None, components=None, params=params)
            except ValueError:
                logger.error("No data for any component: m%s" % str(mov_stack))
                continue
            xr_save_dvv(root, lineage_names, step.step_name, "ALL", mov_stack, dvv)
            del dvv

        massive_update_job(db, jobs, "D")

    logger.info('*** Finished: Compute DV/V ***')


if __name__ == "__main__":
    main()
