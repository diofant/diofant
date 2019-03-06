"""
This module implements methods for solving equations and inequalities.
"""

from .solvers import solve, solve_linear, checksol
from .diophantine import diophantine
from .recurr import rsolve, rsolve_poly, rsolve_ratio, rsolve_hyper
from .ode import checkodesol, classify_ode, dsolve, homogeneous_order
from .polysys import solve_poly_system, solve_linear_system
from .pde import (pde_separate, pde_separate_add, pde_separate_mul,
                  pdsolve, classify_pde, checkpdesol)
from .deutils import ode_order
from .inequalities import reduce_inequalities
