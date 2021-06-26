"""Ground types for various mathematical domains in Diofant."""

import builtins
import fractions
from math import factorial as python_factorial  # noqa: F401
from math import gcd as python_gcd  # noqa: F401
from math import isqrt as python_sqrt  # noqa: F401

from ..core.compatibility import HAS_GMPY
from ..core.numbers import Float as DiofantReal  # noqa: F401
from ..core.numbers import Integer as DiofantInteger  # noqa: F401
from ..core.numbers import Rational as DiofantRational  # noqa: F401
from ..core.numbers import igcdex as python_gcdex  # noqa: F401


PythonInteger = builtins.int
PythonReal = builtins.float
PythonComplex = builtins.complex
PythonRational = fractions.Fraction


if HAS_GMPY:
    from gmpy2 import denom as gmpy_denom
    from gmpy2 import fac as gmpy_factorial
    from gmpy2 import gcd as gmpy_gcd
    from gmpy2 import gcdext as gmpy_gcdex
    from gmpy2 import isqrt as gmpy_sqrt
    from gmpy2 import mpq as GMPYRational  # noqa: N812
    from gmpy2 import mpz as GMPYInteger  # noqa: N812
    from gmpy2 import numer as gmpy_numer
    from gmpy2 import qdiv as gmpy_qdiv
else:
    class _GMPYInteger:
        def __init__(self, obj):
            pass

    class _GMPYRational:
        def __init__(self, obj):
            pass

    GMPYInteger = _GMPYInteger
    GMPYRational = _GMPYRational
    gmpy_factorial = None
    gmpy_numer = None
    gmpy_denom = None
    gmpy_gcdex = None
    gmpy_gcd = None
    gmpy_sqrt = None
    gmpy_qdiv = None
