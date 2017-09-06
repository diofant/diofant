"""Some calculus-related methods waiting to find a better place in the
Diofant modules tree.
"""

from .euler import euler_equations  # noqa: F401
from .singularities import singularities  # noqa: F401
from .finite_diff import finite_diff_weights, apply_finite_diff, as_finite_diff  # noqa: F401
from .optimization import minimize, maximize  # noqa: F401
