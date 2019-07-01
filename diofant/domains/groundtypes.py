"""Ground types for various mathematical domains in Diofant. """

__all__ = ()

import builtins
import fractions
from math import factorial as python_factorial  # noqa: F401
from math import gcd as python_gcd  # noqa: F401

from ..core.compatibility import HAS_GMPY
from ..core.numbers import Float as DiofantReal  # noqa: F401
from ..core.numbers import Integer as DiofantInteger  # noqa: F401
from ..core.numbers import Rational as DiofantRational  # noqa: F401
from ..core.numbers import igcdex as python_gcdex  # noqa: F401
from ..core.numbers import ilcm as python_lcm  # noqa: F401
from ..core.power import isqrt as python_sqrt  # noqa: F401


PythonInteger = builtins.int
PythonReal = builtins.float
PythonComplex = builtins.complex
PythonRational = fractions.Fraction


if HAS_GMPY:
    from gmpy2 import (  # noqa: N812
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
