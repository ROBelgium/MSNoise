"""MSNoise core infrastructure — re-exported for convenient access.

Submodules
----------
db             database connection, engine, logging
config         configuration CRUD, parameters, plot filename helpers
stations       station and data-availability management
workflow       workflow topology, job management, lineage, scheduling
signal         DSP, preprocessing helpers, stacking, PSD RMS, Wiener filter
io             xarray I/O for all result types (CCF/MWCS/DTT/STR/WCT/DVV/PSD)
fdsn           FDSN/EIDA waveform fetching and DataSource helpers
compute        core CC/whitening/MWCS computation (moved from move2obspy)
preprocessing  waveform preprocessing pipeline

Usage
-----
    from msnoise.core.db           import connect, get_logger
    from msnoise.core.config       import get_config, get_params, build_plot_outfile
    from msnoise.core.stations     import get_stations, get_station_pairs
    from msnoise.core.workflow     import get_next_lineage_batch, reset_jobs
    from msnoise.core.signal       import stack, winsorizing, xwt, psd_rms
    from msnoise.core.io           import xr_save_ccf, xr_get_ccf, aggregate_dvv_pairs
    from msnoise.core.compute      import myCorr2, whiten2, mwcs
    from msnoise.core.preprocessing import preprocess
"""

from .db import *            # noqa: F401,F403
from .config import *        # noqa: F401,F403
from .stations import *      # noqa: F401,F403
from .workflow import *      # noqa: F401,F403
from .signal import *        # noqa: F401,F403
from .io import *            # noqa: F401,F403
from .fdsn import *          # noqa: F401,F403
from .compute import *       # noqa: F401,F403
from .preprocessing import * # noqa: F401,F403
