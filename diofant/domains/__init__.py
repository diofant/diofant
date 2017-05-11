"""Implementation of mathematical domains. """

from . import domain
from . import finitefield
from . import integerring
from . import rationalfield
from . import realfield
from . import complexfield
from . import pythonfinitefield
from . import gmpyfinitefield
from . import pythonintegerring
from . import gmpyintegerring
from . import pythonrationalfield
from . import gmpyrationalfield
from . import algebraicfield
from . import polynomialring
from . import fractionfield
from . import expressiondomain

from .domain import Domain
from .finitefield import FiniteField
from .integerring import IntegerRing
from .rationalfield import RationalField
from .realfield import RealField
from .complexfield import ComplexField
from .pythonfinitefield import PythonFiniteField
from .gmpyfinitefield import GMPYFiniteField
from .pythonintegerring import PythonIntegerRing
from .gmpyintegerring import GMPYIntegerRing
from .pythonrationalfield import PythonRationalField
from .gmpyrationalfield import GMPYRationalField
from .algebraicfield import AlgebraicField
from .polynomialring import PolynomialRing
from .fractionfield import FractionField
from .expressiondomain import ExpressionDomain
from .groundtypes import PythonRational

FF_python = PythonFiniteField
FF_gmpy = GMPYFiniteField

ZZ_python = PythonIntegerRing
ZZ_gmpy = GMPYIntegerRing

QQ_python = PythonRationalField
QQ_gmpy = GMPYRationalField

RR = RealField()
CC = ComplexField()

from ..core.compatibility import GROUND_TYPES

_GROUND_TYPES_MAP = {'gmpy': (FF_gmpy, ZZ_gmpy(), QQ_gmpy()),
                     'python': (FF_python, ZZ_python(), QQ_python())}

FF, ZZ, QQ = _GROUND_TYPES_MAP[GROUND_TYPES]
GF = FF
EX = ExpressionDomain()

# Remove clash with functions of polys module:
del ring
del field
