from .base import *
from . import fluids, grids, models, scalers, utils

VERSION = (0, 0, 1)
__version__ = ".".join([str(x) for x in VERSION])
