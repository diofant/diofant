"""
Package for symbolic vector algebra in 3D.
"""

from .coordsysrect import CoordSysCartesian
from .deloperator import Del
from .dyadic import BaseDyadic, Dyadic, DyadicAdd, DyadicMul, DyadicZero
from .functions import (curl, divergence, express, gradient, is_conservative,
                        is_solenoidal, matrix_to_vector, scalar_potential,
                        scalar_potential_difference)
from .orienters import (AxisOrienter, BodyOrienter, QuaternionOrienter,
                        SpaceOrienter)
from .point import Point
from .scalar import BaseScalar
from .vector import BaseVector, Vector, VectorAdd, VectorMul, VectorZero


__all__ = ('CoordSysCartesian', 'Del', 'BaseDyadic', 'Dyadic', 'DyadicAdd',
           'DyadicMul', 'DyadicZero', 'curl', 'divergence', 'express',
           'gradient', 'is_conservative', 'is_solenoidal', 'matrix_to_vector',
           'scalar_potential', 'scalar_potential_difference', 'AxisOrienter',
           'BodyOrienter', 'QuaternionOrienter', 'SpaceOrienter', 'Point',
           'BaseScalar', 'BaseVector', 'Vector', 'VectorAdd',
           'VectorMul', 'VectorZero')
