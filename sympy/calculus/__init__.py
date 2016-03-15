"""Some calculus-related methods waiting to find a better place in the
SymPy modules tree.
"""

from .euler import euler_equations
from .singularities import singularities
from .finite_diff import finite_diff_weights, apply_finite_diff, as_finite_diff
from .optimization import minimize, maximize
