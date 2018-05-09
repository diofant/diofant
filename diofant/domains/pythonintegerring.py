"""Implementaton of :class:`PythonIntegerRing` class. """

from .groundtypes import (PythonInteger, python_factorial, python_gcd,
                          python_gcdex, python_lcm, python_sqrt)
from .integerring import IntegerRing


__all__ = ('PythonIntegerRing',)


class PythonIntegerRing(IntegerRing):
    """Integer ring based on Python's integers. """

    dtype = PythonInteger
    zero = dtype(0)
    one = dtype(1)

    def __init__(self):
        """Allow instantiation of this domain. """

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
