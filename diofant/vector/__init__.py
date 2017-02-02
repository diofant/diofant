"""
Package for symbolic vector algebra in 3D.
"""

from .vector import Vector, VectorAdd, VectorMul, BaseVector, VectorZero
from .dyadic import Dyadic, DyadicAdd, DyadicMul, BaseDyadic, DyadicZero
from .scalar import BaseScalar
from .deloperator import Del
from .coordsysrect import CoordSysCartesian
from .functions import (express, matrix_to_vector, curl, divergence,
                        gradient, is_conservative, is_solenoidal,
                        scalar_potential, scalar_potential_difference)
from .point import Point
from .orienters import (AxisOrienter, BodyOrienter, SpaceOrienter,
                        QuaternionOrienter)
