"""MSNoise plugin stable re-export surface.

Plugin packages should import from here rather than from ``msnoise.core.*``
directly, so that internal restructuring in MSNoise does not silently break
plugins.

**Entry point groups** (declare in your ``pyproject.toml``):

.. code-block:: toml

    [project.entry-points."msnoise.plugins.table_def"]
    myplugin = "myplugin.tabledef:get_table_def"

    [project.entry-points."msnoise.plugins.workflow_chains"]
    myplugin = "myplugin.workflow:get_chains"

    [project.entry-points."msnoise.plugins.workflow_order"]
    myplugin = "myplugin.workflow:get_order"

    [project.entry-points."msnoise.plugins.jobtypes"]
    myplugin = "myplugin.jobs:get_jobtypes"

    [project.entry-points."msnoise.plugins.commands"]
    myplugin = "myplugin.cli:cli"
"""

# ── Database & session ────────────────────────────────────────────────────────
from ..core.db import connect, get_logger                          # noqa: F401

# ── Config & params ───────────────────────────────────────────────────────────
from ..core.config import (                                        # noqa: F401
    get_config,
    update_config,
    create_config_set,
    get_params,
    get_merged_params_for_lineage,
    get_config_set_details,
)

# ── Workflow topology ─────────────────────────────────────────────────────────
from ..core.workflow import (                                       # noqa: F401
    get_workflow_chains,
    get_workflow_order,
    get_workflow_steps,
    get_workflow_links,
    get_step_abbrevs,
    is_terminal_category,
    is_entry_category,
)

# ── Jobs ──────────────────────────────────────────────────────────────────────
from ..core.workflow import (                                       # noqa: F401
    get_next_lineage_batch,
    is_next_job_for_step,
    massive_update_job,
    massive_insert_job,
    propagate_downstream,
    update_job,
    reset_jobs,
    extend_days,
    get_t_axis,
)

# ── I/O ───────────────────────────────────────────────────────────────────────
from ..core.io import (                                            # noqa: F401
    xr_get_ccf,
    xr_save_ccf,
    xr_get_ref,
    xr_save_ref,
    xr_get_mwcs,
    xr_save_mwcs,
    xr_get_dtt,
    xr_save_dtt,
    xr_save_stretching,
    xr_load_wct,
    xr_save_wct,
    xr_get_wct_dtt,
    xr_save_wct_dtt,
    xr_get_dvv_agg,
    xr_save_dvv_agg,
    xr_save_psd,
    xr_load_psd,
    xr_save_rms,
    xr_load_rms,
    aggregate_dvv_pairs,
)

# ── Stations & data ───────────────────────────────────────────────────────────
from ..core.stations import (                                      # noqa: F401
    get_stations,
    get_station_pairs,
    get_data_availability,
    get_interstation_distance,
    get_station,
)

# ── Signal / compute ──────────────────────────────────────────────────────────
from ..core.signal import validate_stack_data, psd_df_rms          # noqa: F401
from ..core.compute import (                                       # noqa: F401
    compute_wct_dtt_batch,
    build_wct_dtt_dataset,
    resolve_wct_lag_min,
)
