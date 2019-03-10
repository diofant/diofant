"""
This module implements methods for solving equations and inequalities.
"""

from .deutils import ode_order
from .diophantine import diophantine
from .inequalities import reduce_inequalities
from .ode import checkodesol, classify_ode, dsolve, homogeneous_order
from .pde import (checkpdesol, classify_pde, pde_separate, pde_separate_add,
                  pde_separate_mul, pdsolve)
from .polysys import solve_linear_system, solve_poly_system
from .recurr import rsolve, rsolve_hyper, rsolve_poly, rsolve_ratio
from .solvers import checksol, solve, solve_linear
