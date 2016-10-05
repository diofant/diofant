"""Ground types for various mathematical domains in Diofant. """

__all__ = ()

import builtins
import fractions

import mpmath.libmp as mlib

from diofant.core.compatibility import HAS_GMPY

PythonInteger = builtins.int
PythonReal = builtins.float
PythonComplex = builtins.complex
PythonRational = fractions.Fraction

from diofant.core.numbers import (
    igcdex as python_gcdex,
    igcd as python_gcd,
    ilcm as python_lcm)
from diofant import (
    Float as DiofantReal,
    Integer as DiofantInteger,
    Rational as DiofantRational)

if HAS_GMPY:
    from gmpy2 import (
        mpz as GMPYInteger,
        mpq as GMPYRational,
        fac as gmpy_factorial,
        numer as gmpy_numer,
        denom as gmpy_denom,
        gcdext as gmpy_gcdex,
        gcd as gmpy_gcd,
        lcm as gmpy_lcm,
        isqrt as gmpy_sqrt,
        qdiv as gmpy_qdiv)
else:
    class GMPYInteger:
        def __init__(self, obj):
            pass

    class GMPYRational:
        def __init__(self, obj):
            pass

    gmpy_factorial = None
    gmpy_numer = None
    gmpy_denom = None
    gmpy_gcdex = None
    gmpy_gcd = None
    gmpy_lcm = None
    gmpy_sqrt = None
    gmpy_qdiv = None


def python_sqrt(n):
    return int(mlib.isqrt(n))


def python_factorial(n):
    return int(mlib.ifac(n))
