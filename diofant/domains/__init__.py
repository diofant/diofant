"""Implementation of mathematical domains. """

from . import domain  # noqa: F401
from . import finitefield  # noqa: F401
from . import integerring  # noqa: F401
from . import rationalfield  # noqa: F401
from . import realfield  # noqa: F401
from . import complexfield  # noqa: F401
from . import algebraicfield  # noqa: F401
from . import expressiondomain  # noqa: F401

from .domain import Domain  # noqa: F401
from .finitefield import (FiniteField, GMPYFiniteField as FF_gmpy,  # noqa: F401
                          PythonFiniteField as FF_python)
from .integerring import IntegerRing,  ZZ_gmpy, ZZ_python  # noqa: F401
from .rationalfield import RationalField, QQ_gmpy, QQ_python  # noqa: F401
from .realfield import RR, RealField  # noqa: F401
from .complexfield import CC, ComplexField  # noqa: F401
from .algebraicfield import (AlgebraicField, ComplexAlgebraicField,  # noqa: F401
                             RealAlgebraicField)
from .expressiondomain import EX, ExpressionDomain  # noqa: F401
from .groundtypes import PythonRational  # noqa: F401
from ..core.compatibility import GROUND_TYPES

_GROUND_TYPES_MAP = {'gmpy': (FF_gmpy, ZZ_gmpy, QQ_gmpy),
                     'python': (FF_python, ZZ_python, QQ_python)}
FF, ZZ, QQ = _GROUND_TYPES_MAP[GROUND_TYPES]
GF = FF

# Remove clash with functions of polys module:
del ring
del field
