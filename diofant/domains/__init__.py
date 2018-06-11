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
from .finitefield import FiniteField, GMPYFiniteField, PythonFiniteField  # noqa: F401
from .integerring import GMPYIntegerRing, IntegerRing, PythonIntegerRing  # noqa: F401
from .rationalfield import GMPYRationalField, PythonRationalField, RationalField  # noqa: F401
from .realfield import RealField
from .complexfield import ComplexField
from .algebraicfield import AlgebraicField  # noqa: F401
from .expressiondomain import ExpressionDomain
from .groundtypes import PythonRational  # noqa: F401

FF_python = PythonFiniteField
FF_gmpy = GMPYFiniteField

ZZ_python = PythonIntegerRing()
ZZ_gmpy = GMPYIntegerRing()

QQ_python = PythonRationalField()
QQ_gmpy = GMPYRationalField()

RR = RealField()
CC = ComplexField()

from ..core.compatibility import GROUND_TYPES

_GROUND_TYPES_MAP = {'gmpy': (FF_gmpy, ZZ_gmpy, QQ_gmpy),
                     'python': (FF_python, ZZ_python, QQ_python)}

FF, ZZ, QQ = _GROUND_TYPES_MAP[GROUND_TYPES]
GF = FF
EX = ExpressionDomain()

# Remove clash with functions of polys module:
del ring
del field
