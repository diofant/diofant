"""Some calculus-related methods waiting to find a better place in the
Diofant modules tree.
"""

from .euler import euler_equations
from .finite_diff import apply_finite_diff, as_finite_diff, finite_diff_weights
from .optimization import maximize, minimize
from .singularities import singularities
