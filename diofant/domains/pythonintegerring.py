"""Implementaton of :class:`PythonIntegerRing` class. """

from ..polys.polyerrors import CoercionFailed
from .groundtypes import (DiofantInteger, PythonInteger, python_factorial,
                          python_gcd, python_gcdex, python_lcm, python_sqrt)
from .integerring import IntegerRing


__all__ = ('PythonIntegerRing',)


class PythonIntegerRing(IntegerRing):
    """Integer ring based on Python's integers. """

    dtype = PythonInteger
    zero = dtype(0)
    one = dtype(1)

    def __init__(self):
        """Allow instantiation of this domain. """

    def to_diofant(self, a):
        """Convert ``a`` to a Diofant object. """
        return DiofantInteger(a)

    def from_diofant(self, a):
        """Convert Diofant's Integer to ``dtype``. """
        if a.is_Integer:
            return PythonInteger(a.p)
        elif a.is_Float and int(a) == a:
            return PythonInteger(int(a))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def from_PythonFiniteField(self, a, K0):
        """Convert ``ModularInteger(int)`` to Python's ``int``. """
        return a.to_int()

    def from_PythonIntegerRing(self, a, K0):
        """Convert Python's ``int`` to Python's ``int``. """
        return a

    def from_PythonRationalField(self, a, K0):
        """Convert Python's ``Fraction`` to Python's ``int``. """
        if a.denominator == 1:
            return a.numerator

    def from_GMPYFiniteField(self, a, K0):
        """Convert ``ModularInteger(mpz)`` to Python's ``int``. """
        return PythonInteger(a.to_int())

    def from_GMPYIntegerRing(self, a, K0):
        """Convert GMPY's ``mpz`` to Python's ``int``. """
        return PythonInteger(a)

    def from_GMPYRationalField(self, a, K0):
        """Convert GMPY's ``mpq`` to Python's ``int``. """
        if a.denominator == 1:
            return PythonInteger(a.numerator)

    def from_RealField(self, a, K0):
        """Convert mpmath's ``mpf`` to Python's ``int``. """
        p, q = K0.to_rational(a)

        if q == 1:
            return PythonInteger(p)

    def gcdex(self, a, b):
        """Compute extended GCD of ``a`` and ``b``. """
        return python_gcdex(a, b)

    def gcd(self, a, b):
        """Compute GCD of ``a`` and ``b``. """
        return python_gcd(a, b)

    def lcm(self, a, b):
        """Compute LCM of ``a`` and ``b``. """
        return python_lcm(a, b)

    def sqrt(self, a):
        """Compute square root of ``a``. """
        return python_sqrt(a)

    def factorial(self, a):
        """Compute factorial of ``a``. """
        return python_factorial(a)
