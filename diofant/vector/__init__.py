from diofant.vector.vector import (Vector, VectorAdd, VectorMul,
                                   BaseVector, VectorZero)
from diofant.vector.dyadic import (Dyadic, DyadicAdd, DyadicMul,
                                   BaseDyadic, DyadicZero)
from diofant.vector.scalar import BaseScalar
from diofant.vector.deloperator import Del
from diofant.vector.coordsysrect import CoordSysCartesian
from diofant.vector.functions import (express, matrix_to_vector,
                                      curl, divergence, gradient,
                                      is_conservative, is_solenoidal,
                                      scalar_potential,
                                      scalar_potential_difference)
from diofant.vector.point import Point
from diofant.vector.orienters import (AxisOrienter, BodyOrienter,
                                      SpaceOrienter, QuaternionOrienter)
