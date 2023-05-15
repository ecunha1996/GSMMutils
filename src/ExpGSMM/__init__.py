
from os.path import join, dirname, abspath
from .model import *
from .io import *
from .experimental import *

SRC_PATH = abspath(join(dirname(__file__), '../'))
DATA_PATH = abspath(join(dirname(__file__), '../../data'))
DFBALAB_PATH = abspath(join(dirname(__file__), '../matlab/DFBAlab'))

