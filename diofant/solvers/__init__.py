"""A module for solving all kinds of equations.

Examples
========

>>> solve(x**5 + 5*x**4 + 10*x**3 + 10*x**2 + 5*x + 1, x)
[{x: -1}]
"""

from .solvers import solve, solve_linear, checksol  # noqa: F401
from .diophantine import diophantine  # noqa: F401
from .recurr import rsolve, rsolve_poly, rsolve_ratio, rsolve_hyper  # noqa: F401
from .ode import checkodesol, classify_ode, dsolve, homogeneous_order  # noqa: F401
from .polysys import solve_poly_system, solve_linear_system  # noqa: F401
from .pde import (pde_separate, pde_separate_add, pde_separate_mul,  # noqa: F401
                  pdsolve, classify_pde, checkpdesol)
from .deutils import ode_order  # noqa: F401
from .inequalities import reduce_inequalities  # noqa: F401
