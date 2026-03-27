"""MSNoise API — compatibility shim.

.. deprecated::
    This module is retained only as a convenience re-export for third-party
    code written against earlier MSNoise versions.  New code should import
    directly from the focused submodules::

        from msnoise.db      import connect, get_logger
        from msnoise.config  import get_config, get_params, build_plot_outfile
        from msnoise.stations import get_stations, get_station_pairs
        from msnoise.workflow import get_next_lineage_batch, reset_jobs
        from msnoise.signal  import stack, winsorizing, xwt
        from msnoise.io      import xr_save_ccf, xr_get_ccf, aggregate_dvv_pairs

Internal MSNoise code imports directly from the submodules above.
"""
import warnings as _warnings

_warnings.warn(
    "Importing from msnoise.api is deprecated. "
    "Import directly from msnoise.db / .config / .stations / .workflow / .signal / .io",
    DeprecationWarning,
    stacklevel=2,
)

from .db import *        # noqa: F401,F403
from .config import *    # noqa: F401,F403
from .stations import *  # noqa: F401,F403
from .workflow import *  # noqa: F401,F403
from .signal import *    # noqa: F401,F403
from .io import *        # noqa: F401,F403

from .msnoise_table_def import (  # noqa: F401
    Job, Station, Config, DataAvailability, WorkflowStep,
    WORKFLOW_CHAINS, WORKFLOW_ORDER,
)
