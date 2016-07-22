"""Implementation of mathematical domains. """

__all__ = []

from . import domain
__all__.extend(domain.__all__)
from .domain import *  # noqa: F403

from . import finitefield
__all__.extend(finitefield.__all__)
from .finitefield import *  # noqa: F403

from . import integerring
__all__.extend(integerring.__all__)
from .integerring import *  # noqa: F403

from . import rationalfield
__all__.extend(rationalfield.__all__)
from .rationalfield import *  # noqa: F403

from . import realfield
__all__.extend(realfield.__all__)
from .realfield import *  # noqa: F403

from . import complexfield
__all__.extend(complexfield.__all__)
from .complexfield import *  # noqa: F403

from . import pythonfinitefield
__all__.extend(pythonfinitefield.__all__)
from .pythonfinitefield import *  # noqa: F403

from . import gmpyfinitefield
__all__.extend(gmpyfinitefield.__all__)
from .gmpyfinitefield import *  # noqa: F403

from . import pythonintegerring
__all__.extend(pythonintegerring.__all__)
from .pythonintegerring import *  # noqa: F403

from . import gmpyintegerring
__all__.extend(gmpyintegerring.__all__)
from .gmpyintegerring import *  # noqa: F403

from . import pythonrationalfield
__all__.extend(pythonrationalfield.__all__)
from .pythonrationalfield import *  # noqa: F403

from . import gmpyrationalfield
__all__.extend(gmpyrationalfield.__all__)
from .gmpyrationalfield import *  # noqa: F403

from . import algebraicfield
__all__.extend(algebraicfield.__all__)
from .algebraicfield import *  # noqa: F403

from . import polynomialring
__all__.extend(polynomialring.__all__)
from .polynomialring import *  # noqa: F403

from . import fractionfield
__all__.extend(fractionfield.__all__)
from .fractionfield import *  # noqa: F403

from . import expressiondomain
__all__.extend(expressiondomain.__all__)
from .expressiondomain import *  # noqa: F403

FF_python = PythonFiniteField
FF_gmpy = GMPYFiniteField

ZZ_python = PythonIntegerRing
ZZ_gmpy = GMPYIntegerRing

QQ_python = PythonRationalField
QQ_gmpy = GMPYRationalField

RR = RealField()
CC = ComplexField()

from .groundtypes import PythonRational

from diofant.core.compatibility import GROUND_TYPES

_GROUND_TYPES_MAP = {
    'gmpy': (FF_gmpy, ZZ_gmpy(), QQ_gmpy()),
    'python': (FF_python, ZZ_python(), QQ_python()),
}

try:
    FF, ZZ, QQ = _GROUND_TYPES_MAP[GROUND_TYPES]
except KeyError:
    raise ValueError("invalid ground types: %s" % GROUND_TYPES)

GF = FF

EX = ExpressionDomain()

__all__.extend([
    "FF_python", "FF_gmpy",
    "ZZ_python", "ZZ_gmpy",
    "QQ_python", "QQ_gmpy",
    "GF", "FF", "ZZ", "QQ", "RR", "CC", "EX",
])
