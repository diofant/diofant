"""Ground types for various mathematical domains in Diofant."""

import builtins
import fractions
from math import factorial as python_factorial
from math import gcd as python_gcd
from math import isqrt as python_sqrt

from ..core.compatibility import HAS_GMPY
from ..core.numbers import Float as DiofantReal
from ..core.numbers import Integer as DiofantInteger
from ..core.numbers import Rational as DiofantRational
from ..core.numbers import igcdex as python_gcdex


__all__ = ('python_factorial', 'python_gcd', 'python_sqrt', 'DiofantReal',
           'DiofantInteger', 'DiofantRational', 'python_gcdex', 'PythonInteger',
           'PythonInteger', 'PythonReal', 'PythonComplex', 'PythonRational',
           'gmpy_denom', 'gmpy_factorial', 'gmpy_gcd', 'gmpy_gcdex',
           'gmpy_sqrt', 'GMPYRational', 'GMPYInteger',
           'gmpy_numer', 'gmpy_qdiv')


PythonInteger = builtins.int
PythonReal = builtins.float
PythonComplex = builtins.complex
PythonRational = fractions.Fraction


if HAS_GMPY:
    # pylint: disable=no-name-in-module
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
