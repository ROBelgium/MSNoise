__author__ = "Thomas LECOCQ, Corentin CAUDRON and Florent BRENGUIER"
__copyright__ = "Copyright 2015, The Authors"
__credits__ = []
__license__ = "GPL"
__version__ = "1.3"
__maintainer__ = "Thomas LECOCQ"
__email__ = "Thomas.Lecocq at seismology.be"
__status__ = "Production"

import os

MSNoisePATH = os.path.realpath(os.path.dirname(__file__))
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
