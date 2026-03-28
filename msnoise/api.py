"""MSNoise API compatibility module.

The only symbol most code needs is ``connect`` — the universal entry point::

    from msnoise.api import connect
    # or equivalently:
    from msnoise import connect
    # or directly:
    from msnoise.core.db import connect

All other symbols are re-exported here for backward compatibility with code
written against earlier MSNoise versions.  New code should import directly
from the appropriate ``msnoise.core.*`` submodule.
"""

# The one symbol worth importing from here
from .core.db import connect  # noqa: F401

# Full backward-compat re-export (no deprecation — just convenience)
from .core import *  # noqa: F401,F403

from .msnoise_table_def import (  # noqa: F401
    Job, Station, Config, DataAvailability, WorkflowStep, DataSource,
    WORKFLOW_CHAINS, WORKFLOW_ORDER,
)
