from os.path import join, dirname, abspath

from .experimental import *
from .io import *
from .model import *
from .utils import *

SRC_PATH = abspath(join(dirname(__file__), '../'))
DATA_PATH = abspath(join(dirname(__file__), '../../data'))

