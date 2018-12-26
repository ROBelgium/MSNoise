__author__ = "Thomas LECOCQ, Corentin CAUDRON and Florent BRENGUIER"
__copyright__ = "Copyright 2015, The Authors"
__credits__ = []
__license__ = "GPL"
__version__ = "1.3"
__maintainer__ = "Thomas LECOCQ"
__email__ = "Thomas.Lecocq at seismology.be"
__status__ = "Production"


class MSNoiseError(Exception):
    pass

class DBConfigNotFoundError(MSNoiseError):
    pass

class FatalError(MSNoiseError):
    pass
