from os.path import join, dirname, abspath

from .experimental import *
from .io import *
from .model import *

SRC_PATH = abspath(join(dirname(__file__), '../'))
DATA_PATH = abspath(join(dirname(__file__), '../../data'))
DFBALAB_PATH = abspath(join(dirname(__file__), '../matlab/DFBAlab'))
CONFIG_PATH = abspath(join(dirname(__file__), '../../config'))