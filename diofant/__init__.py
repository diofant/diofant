"""Diofant is a Python library for symbolic mathematics. """

import os
DIOFANT_DEBUG = os.getenv('DIOFANT_DEBUG', 'False') != 'False'
del os

import pkg_resources
__version__ = pkg_resources.get_distribution(__name__).version
del pkg_resources

from .core import *  # noqa: F401,F403
from .logic import *  # noqa: F401,F403
from .polys import *  # noqa: F401,F403
from .domains import *  # noqa: F401,F403
from .series import *  # noqa: F401,F403
from .functions import *  # noqa: F401,F403
from .ntheory import *  # noqa: F401,F403
from .concrete import *  # noqa: F401,F403
from .simplify import *  # noqa: F401,F403
from .sets import *  # noqa: F401,F403
from .solvers import *  # noqa: F401,F403
from .matrices import *  # noqa: F401,F403
from .geometry import *  # noqa: F401,F403
from .utilities import *  # noqa: F401,F403
from .integrals import *  # noqa: F401,F403
from .tensor import *  # noqa: F401,F403
from .parsing import *  # noqa: F401,F403
from .calculus import *  # noqa: F401,F403
from .combinatorics import *  # noqa: F401,F403
from .plotting import *  # noqa: F401,F403
from .printing import *  # noqa: F401,F403
from .interactive import *  # noqa: F401,F403
