"""MSNoise public API — compatibility shim.

.. deprecated::
    This module is retained for third-party code written against earlier
    MSNoise versions.  New code should import from ``msnoise.core``::

        from msnoise.core.db      import connect, get_logger
        from msnoise.core.config  import get_config, get_params, build_plot_outfile
        from msnoise.core.stations import get_stations, get_station_pairs
        from msnoise.core.workflow import get_next_lineage_batch, reset_jobs
        from msnoise.core.signal  import stack, winsorizing, xwt
        from msnoise.core.io      import xr_save_ccf, xr_get_ccf, aggregate_dvv_pairs

Internal MSNoise code imports directly from ``msnoise.core.*``.
"""
import warnings as _warnings

_warnings.warn(
    "Importing from msnoise.api is deprecated. "
    "Import directly from msnoise.core.db / .config / .stations / .workflow / .signal / .io",
    DeprecationWarning,
    stacklevel=2,
)

from .core import *  # noqa: F401,F403

from .msnoise_table_def import (  # noqa: F401
    Job, Station, Config, DataAvailability, WorkflowStep,
    WORKFLOW_CHAINS, WORKFLOW_ORDER,
)
