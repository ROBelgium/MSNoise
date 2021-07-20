# TODO add documentation :-)

import os
import traceback
import matplotlib
matplotlib.use("Agg")
import numpy as np
import datetime
from obspy.core import UTCDateTime, read, Stream
from obspy.signal import PPSD

import warnings
warnings.filterwarnings("ignore")


from .api import *

import logbook




def main(loglevel="INFO", njobs_per_worker=9999):
    logger = logbook.Logger("msnoise")
    # Reconfigure logger to show the pid number in log records
    logger = get_logger('msnoise.compute_psd_child', loglevel,
                        with_pid=True)
    logger.info('*** Starting: HDF to RMS ***')
    db = connect()
    params = get_params(db)
    files = glob.glob(os.path.join("PSD", "RMS", params.qc_rms_type, '*.h5'))
    for file in sorted(files):
        print(file)
        store = hdf_open_store_from_fn(file, "r")\

        store.RMS.to_csv(file.replace(".h5", ".csv"))
        store.close()
        del store


    # hdf_insert_or_update(RMSstore, "PSD", RMS)
    # hdf_close_store(RMSstore)
    # del RMS, RMSstore
    # del data
    #     massive_update_job(db, jobs, "D")