from .base import *
from . import grids, fluids, wells, models, utils, profme, plots, pinns, scalers

VERSION = (0, 0, 1)
__version__ = ".".join([str(x) for x in VERSION])
