"""Implementation of mathematical domains. """

from ..core.compatibility import GROUND_TYPES
from . import (algebraicfield, complexfield, domain, expressiondomain,
               finitefield, integerring, rationalfield, realfield)
from .algebraicfield import (AlgebraicField, ComplexAlgebraicField,
                             RealAlgebraicField)
from .complexfield import CC, ComplexField
from .domain import Domain
from .expressiondomain import EX, ExpressionDomain
from .finitefield import FiniteField
from .finitefield import GMPYFiniteField as FF_gmpy
from .finitefield import PythonFiniteField as FF_python
from .groundtypes import PythonRational
from .integerring import IntegerRing, ZZ_gmpy, ZZ_python
from .rationalfield import QQ_gmpy, QQ_python, RationalField
from .realfield import RR, RealField


_GROUND_TYPES_MAP = {'gmpy': (FF_gmpy, ZZ_gmpy, QQ_gmpy),
                     'python': (FF_python, ZZ_python, QQ_python)}
FF, ZZ, QQ = _GROUND_TYPES_MAP[GROUND_TYPES]
GF = FF

# Remove clash with functions of polys module:
del ring
del field
