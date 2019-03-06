"""Implementation of mathematical domains. """

from . import domain
from . import finitefield
from . import integerring
from . import rationalfield
from . import realfield
from . import complexfield
from . import algebraicfield
from . import expressiondomain

from .domain import Domain
from .finitefield import (FiniteField, GMPYFiniteField as FF_gmpy,
                          PythonFiniteField as FF_python)
from .integerring import IntegerRing,  ZZ_gmpy, ZZ_python
from .rationalfield import RationalField, QQ_gmpy, QQ_python
from .realfield import RR, RealField
from .complexfield import CC, ComplexField
from .algebraicfield import (AlgebraicField, ComplexAlgebraicField,
                             RealAlgebraicField)
from .expressiondomain import EX, ExpressionDomain
from .groundtypes import PythonRational
from ..core.compatibility import GROUND_TYPES

_GROUND_TYPES_MAP = {'gmpy': (FF_gmpy, ZZ_gmpy, QQ_gmpy),
                     'python': (FF_python, ZZ_python, QQ_python)}
FF, ZZ, QQ = _GROUND_TYPES_MAP[GROUND_TYPES]
GF = FF

# Remove clash with functions of polys module:
del ring
del field
