"""Diofant is a Python library for symbolic mathematics. It aims to become a
full-featured computer algebra system (CAS) while keeping the code as
simple as possible in order to be comprehensible and easily extensible.
Diofant is written entirely in Python and does not require any external
libraries, except optionally for plotting support.
"""

__version__ = "0.8.0a1"


def __diofant_debug():
    # helper function so we don't import os globally
    import os
    debug_str = os.getenv('DIOFANT_DEBUG', 'False')
    if debug_str in ('True', 'False'):
        return eval(debug_str)
    else:
        raise RuntimeError("unrecognized value for DIOFANT_DEBUG: %s" %
                           debug_str)
DIOFANT_DEBUG = __diofant_debug()

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

evalf._create_evalf_table()
