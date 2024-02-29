import traceback

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from obspy.signal.regression import linear_regression
from .api import *

import logbook


def main(interval=1, loglevel="INFO"):
    logger = logbook.Logger(__name__)
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_dtt_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: Compute DV/V ***')
    db = connect()
    params = get_params(db)
    mov_stacks = params.mov_stack
    filters = get_filters(db, all=False)

    for f in filters:
        filterid = int(f.ref)

        for mov_stack in mov_stacks:
            for components in params.all_components:
                logger.debug("Processing f%i m%s %s" % (filterid,   mov_stack, components))
                try:
                    dvv = compute_dvv(db, filterid, mov_stack, pairs=None,
                                  components=components, params=params)
                except ValueError:
                    traceback.print_exc()
                    logger.error("No data for f%i m%s: %s" % (filterid, mov_stack, components))
                    continue
                xr_save_dvv(components, filterid, mov_stack, dvv)
                del dvv
            try:
                dvv = compute_dvv(db, filterid, mov_stack, pairs=None,
                              components=None, params=params)
            except ValueError:
                logger.error("No data for any component: f%i m%s" % (filterid, mov_stack))
                continue
            xr_save_dvv("ALL", filterid, mov_stack, dvv)
            del dvv
    
    logger.info('*** Finished: Compute DV/V ***')


if __name__ == "__main__":
    main()
