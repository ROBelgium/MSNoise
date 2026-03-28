__author__ = "Thomas LECOCQ, Corentin CAUDRON and Florent BRENGUIER"
__copyright__ = "Copyright 2015-2021, The Authors"
__credits__ = []
__license__ = "GPL"
try:
    from ._version import version as __version__
except ImportError:
    __version__ = "0.0.0.dev0"
__maintainer__ = "Thomas LECOCQ"
__email__ = "Thomas.Lecocq at seismology.be"
__status__ = "Production"


class MSNoiseError(Exception):
    pass

class DBConfigNotFoundError(MSNoiseError):
    pass

class FatalError(MSNoiseError):
    pass

# Convenience: connect is the universal entry point
from .core.db import connect  # noqa: F401

# Convenience: connect is the universal entry point
from .core.db import connect  # noqa: F401
