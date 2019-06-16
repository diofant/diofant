"""Diofant is a Python library for symbolic mathematics. """

import os
DIOFANT_DEBUG = os.getenv('DIOFANT_DEBUG', 'False') != 'False'
del os

import pkg_resources
__version__ = pkg_resources.get_distribution(__name__).version
del pkg_resources

from .core import *
from .logic import *
from .polys import *
from .domains import *
from .series import *
from .functions import *
from .ntheory import *
from .concrete import *
from .simplify import *
from .sets import *
from .solvers import *
from .matrices import *
from .geometry import *
from .utilities import *
from .integrals import *
from .tensor import *
from .parsing import *
from .calculus import *
from .combinatorics import *
from .plotting import *
from .printing import *
from .interactive import *
