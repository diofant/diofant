"""
Package for symbolic vector algebra in 3D.
"""

from .vector import Vector, VectorAdd, VectorMul, BaseVector, VectorZero  # noqa: F401
from .dyadic import Dyadic, DyadicAdd, DyadicMul, BaseDyadic, DyadicZero  # noqa: F401
from .scalar import BaseScalar  # noqa: F401
from .deloperator import Del  # noqa: F401
from .coordsysrect import CoordSysCartesian  # noqa: F401
from .functions import (express, matrix_to_vector, curl, divergence,  # noqa: F401
                        gradient, is_conservative, is_solenoidal,
                        scalar_potential, scalar_potential_difference)
from .point import Point  # noqa: F401
from .orienters import (AxisOrienter, BodyOrienter, SpaceOrienter,  # noqa: F401
                        QuaternionOrienter)
