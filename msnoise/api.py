"""MSNoise public API.

The universal entry point::

    from msnoise.api import connect
    # or equivalently:
    from msnoise import connect
    # or directly:
    from msnoise.core.db import connect

All other public symbols live in ``msnoise.core.*`` and are importable from
there directly.  This module re-exports the full ``core`` namespace for
convenience::

    from msnoise.api import get_params, get_next_lineage_batch, xr_get_ccf
"""

# The one symbol most code needs
from .core.db import connect  # noqa: F401

# Full core namespace re-export (respects __all__ in each submodule)
from .core import *  # noqa: F401,F403
