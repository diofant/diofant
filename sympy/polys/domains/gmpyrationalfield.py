"""Implementaton of :class:`GMPYRationalField` class. """

from __future__ import print_function, division

from sympy.polys.domains.rationalfield import RationalField
from sympy.polys.domains.groundtypes import (
    GMPYRational, SymPyRational,
    gmpy_numer, gmpy_denom, gmpy_factorial, gmpy_qdiv,
)

from sympy.polys.polyerrors import CoercionFailed
from sympy.utilities import public


@public
class GMPYRationalField(RationalField):
    """Rational field based on GMPY mpq class. """

    dtype = GMPYRational
    zero = dtype(0)
    one = dtype(1)
    tp = type(one)
    alias = 'QQ_gmpy'

    def __init__(self):
        pass

    def get_ring(self):
        """Returns ring associated with ``self``. """
        from sympy.polys.domains import GMPYIntegerRing
        return GMPYIntegerRing()

    def to_sympy(self, a):
        """Convert `a` to a SymPy object. """
        return SymPyRational(int(gmpy_numer(a)),
                             int(gmpy_denom(a)))

    def from_sympy(self, a):
        """Convert SymPy's Integer to `dtype`. """
        if a.is_Rational:
            return GMPYRational(a.p, a.q)
        elif a.is_Float:
            from sympy.polys.domains import RR
            return GMPYRational(*RR.to_rational(a))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def from_ZZ_python(K1, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return GMPYRational(a)

    def from_QQ_python(K1, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return GMPYRational(a.numerator, a.denominator)

    def from_ZZ_gmpy(K1, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return GMPYRational(a)

    def from_QQ_gmpy(K1, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return a

    def from_RealField(K1, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return GMPYRational(*K0.to_rational(a))

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies `__div__`.  """
        return GMPYRational(gmpy_qdiv(a, b))

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies `__div__`. """
        return GMPYRational(gmpy_qdiv(a, b))

    def rem(self, a, b):
        """Remainder of `a` and `b`, implies nothing.  """
        return self.zero

    def div(self, a, b):
        """Division of `a` and `b`, implies `__div__`. """
        return GMPYRational(gmpy_qdiv(a, b)), self.zero

    def numer(self, a):
        """Returns numerator of `a`. """
        return a.numerator

    def denom(self, a):
        """Returns denominator of `a`. """
        return a.denominator

    def factorial(self, a):
        """Returns factorial of `a`. """
        return GMPYRational(gmpy_factorial(int(a)))
