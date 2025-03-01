"""Ground types for various mathematical domains in Diofant."""

import builtins
import fractions
from math import factorial as python_factorial
from math import gcd as python_gcd
from math import isqrt as python_sqrt

from ..core.numbers import Float as DiofantReal
from ..core.numbers import Integer as DiofantInteger
from ..core.numbers import Rational as DiofantRational
from ..core.numbers import igcdex as python_gcdex
from ..external import gmpy


__all__ = ('python_factorial', 'python_gcd', 'python_sqrt', 'DiofantReal',
           'DiofantInteger', 'DiofantRational', 'python_gcdex', 'PythonInteger',
           'PythonInteger', 'PythonReal', 'PythonComplex', 'PythonRational',
           'gmpy_factorial', 'gmpy_gcd', 'gmpy_gcdex',
           'gmpy_sqrt', 'GMPYRational', 'GMPYInteger')


PythonInteger = builtins.int
PythonReal = builtins.float
PythonComplex = builtins.complex
PythonRational = fractions.Fraction


if gmpy:
    gmpy_gcd = gmpy.gcd
    gmpy_gcdex = gmpy.gcdext
    gmpy_factorial = gmpy.fac
    gmpy_sqrt = gmpy.isqrt
    GMPYRational = gmpy.mpq
    GMPYInteger = gmpy.mpz
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
    gmpy_gcdex = None
    gmpy_gcd = None
    gmpy_sqrt = None
