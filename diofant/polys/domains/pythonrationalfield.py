"""Implementation of :class:`PythonRationalField` class. """

from diofant.polys.domains.rationalfield import RationalField
from diofant.polys.domains.groundtypes import (PythonInteger, PythonRational,
                                               DiofantRational)
from diofant.polys.polyerrors import CoercionFailed
from diofant.utilities import public


@public
class PythonRationalField(RationalField):
    """Rational field based on Python rational number type. """

    dtype = PythonRational
    zero = dtype(0)
    one = dtype(1)
    alias = 'QQ_python'

    def __init__(self):
        pass

    def get_ring(self):
        """Returns ring associated with ``self``. """
        from diofant.polys.domains import PythonIntegerRing
        return PythonIntegerRing()

    def to_diofant(self, a):
        """Convert `a` to a Diofant object. """
        return DiofantRational(a.numerator, a.denominator)

    def from_diofant(self, a):
        """Convert Diofant's Rational to `dtype`. """
        if a.is_Rational:
            return PythonRational(a.p, a.q)
        elif a.is_Float:
            from diofant.polys.domains import RR
            p, q = RR.to_rational(a)
            return PythonRational(int(p), int(q))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def from_ZZ_python(self, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return PythonRational(a)

    def from_QQ_python(self, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return a

    def from_ZZ_gmpy(self, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return PythonRational(PythonInteger(a))

    def from_QQ_gmpy(self, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return PythonRational(PythonInteger(a.numer()),
                              PythonInteger(a.denom()))

    def from_RealField(self, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        p, q = K0.to_rational(a)
        return PythonRational(int(p), int(q))

    def numer(self, a):
        """Returns numerator of `a`. """
        return a.numerator

    def denom(self, a):
        """Returns denominator of `a`. """
        return a.denominator
