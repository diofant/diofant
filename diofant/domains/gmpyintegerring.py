"""Implementaton of :class:`GMPYIntegerRing` class. """

from .groundtypes import (GMPYInteger, gmpy_factorial, gmpy_gcd, gmpy_gcdex,
                          gmpy_lcm, gmpy_sqrt)
from .integerring import IntegerRing


__all__ = ('GMPYIntegerRing',)


class GMPYIntegerRing(IntegerRing):
    """Integer ring based on GMPY's integers. """

    dtype = GMPYInteger
    zero = dtype(0)
    one = dtype(1)

    def __init__(self):
        """Allow instantiation of this domain. """

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
