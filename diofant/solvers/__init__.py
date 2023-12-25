"""This module implements methods for solving equations and inequalities."""

from .inequalities import reduce_inequalities
from .ode import dsolve
from .recurr import rsolve
from .solvers import solve


__all__ = 'reduce_inequalities', 'dsolve', 'rsolve', 'solve'
