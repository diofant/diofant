"""Diofant is a Python library for symbolic mathematics. It aims to become a
full-featured computer algebra system (CAS) while keeping the code as
simple as possible in order to be comprehensible and easily extensible.
"""

__version__ = "0.8.0"

import os
DIOFANT_DEBUG = os.getenv('DIOFANT_DEBUG', 'False') != 'False'
del os

from .core import *  # noqa: F403
from .logic import *  # noqa: F403
from .polys import *  # noqa: F403
from .series import *  # noqa: F403
from .functions import *  # noqa: F403
from .ntheory import *  # noqa: F403
from .concrete import *  # noqa: F403
from .simplify import *  # noqa: F403
from .sets import *  # noqa: F403
from .solvers import *  # noqa: F403
from .matrices import *  # noqa: F403
from .geometry import *  # noqa: F403
from .utilities import *  # noqa: F403
from .integrals import *  # noqa: F403
from .tensor import *  # noqa: F403
from .parsing import *  # noqa: F403
from .calculus import *  # noqa: F403
from .combinatorics import *  # noqa: F403
from .plotting import *  # noqa: F403
from .printing import *  # noqa: F403
from .interactive import *  # noqa: F403
