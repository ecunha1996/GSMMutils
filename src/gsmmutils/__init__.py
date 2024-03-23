from os.path import join, dirname, abspath

SRC_PATH = abspath(join(dirname(__file__), '../'))
DATA_PATH = abspath(join(dirname(__file__), '../../data'))
CONFIG_PATH = abspath(join(dirname(__file__), '../../../config'))

from .experimental import *
from .io import *
from .model import *
from .utils import *
from .stats import *
from .annotation import *


def welcome():
    from NicePrinter import cyan, title
    import pyjokes
    print(cyan(title("Welcome To GSMMutils", 90, '=')))
    print(pyjokes.get_joke(language="en"))
