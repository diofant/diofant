"""Implementaton of :class:`GMPYIntegerRing` class. """

from diofant.polys.domains.integerring import IntegerRing
from diofant.polys.domains.groundtypes import (
    GMPYInteger, DiofantInteger,
    gmpy_factorial,
    gmpy_gcdex, gmpy_gcd, gmpy_lcm, gmpy_sqrt)
from diofant.polys.polyerrors import CoercionFailed
from diofant.utilities import public


@public
class GMPYIntegerRing(IntegerRing):
    """Integer ring based on GMPY's ``mpz`` type. """

    dtype = GMPYInteger
    zero = dtype(0)
    one = dtype(1)
    tp = type(one)
    alias = 'ZZ_gmpy'

    def __init__(self):
        """Allow instantiation of this domain. """

    def to_diofant(self, a):
        """Convert ``a`` to a Diofant object. """
        return DiofantInteger(int(a))

    def from_diofant(self, a):
        """Convert Diofant's Integer to ``dtype``. """
        if a.is_Integer:
            return GMPYInteger(a.p)
        elif a.is_Float and int(a) == a:
            return GMPYInteger(int(a))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def from_FF_python(self, a, K0):
        """Convert ``ModularInteger(int)`` to GMPY's ``mpz``. """
        return GMPYInteger(a.to_int())

    def from_ZZ_python(self, a, K0):
        """Convert Python's ``int`` to GMPY's ``mpz``. """
        return GMPYInteger(a)

    def from_QQ_python(self, a, K0):
        """Convert Python's ``Fraction`` to GMPY's ``mpz``. """
        if a.denominator == 1:
            return GMPYInteger(a.numerator)

    def from_FF_gmpy(self, a, K0):
        """Convert ``ModularInteger(mpz)`` to GMPY's ``mpz``. """
        return a.to_int()

    def from_ZZ_gmpy(self, a, K0):
        """Convert GMPY's ``mpz`` to GMPY's ``mpz``. """
        return a

    def from_QQ_gmpy(self, a, K0):
        """Convert GMPY ``mpq`` to GMPY's ``mpz``. """
        if a.denominator == 1:
            return a.numerator

    def from_RealField(self, a, K0):
        """Convert mpmath's ``mpf`` to GMPY's ``mpz``. """
        p, q = K0.to_rational(a)

        if q == 1:
            return GMPYInteger(p)

    def gcdex(self, a, b):
        """Compute extended GCD of ``a`` and ``b``. """
        h, s, t = gmpy_gcdex(a, b)
        return s, t, h

    def gcd(self, a, b):
        """Compute GCD of ``a`` and ``b``. """
        return gmpy_gcd(a, b)

    def lcm(self, a, b):
        """Compute LCM of ``a`` and ``b``. """
        return gmpy_lcm(a, b)

    def sqrt(self, a):
        """Compute square root of ``a``. """
        return gmpy_sqrt(a)

    def factorial(self, a):
        """Compute factorial of ``a``. """
        return gmpy_factorial(a)
