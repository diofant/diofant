"""Implementation of :class:`IntegerRing` class."""

import abc

from ..polys.polyerrors import CoercionFailedError
from .characteristiczero import CharacteristicZero
from .groundtypes import (DiofantInteger, GMPYInteger, PythonInteger,
                          gmpy_factorial, gmpy_gcd, gmpy_gcdex, gmpy_sqrt,
                          python_factorial, python_gcd, python_gcdex,
                          python_sqrt)
from .ring import CommutativeRing
from .simpledomain import SimpleDomain


class IntegerRing(CharacteristicZero, SimpleDomain, CommutativeRing):
    """General class for integer rings."""

    rep = 'ZZ'

    is_IntegerRing = True
    is_Numerical = True

    has_assoc_Ring = True

    @property
    def field(self):
        """Return a field associated with ``self``."""
        from . import QQ
        return QQ

    def to_expr(self, element):
        return DiofantInteger(element)

    def from_expr(self, expr):
        if expr.is_Integer:
            return self.dtype(expr.numerator)
        if expr.is_Float and int(expr) == expr:
            return self.dtype(expr)
        raise CoercionFailedError(f'expected an integer, got {expr}')

    def _from_PythonIntegerRing(self, a, K0):
        return self.dtype(a)
    _from_GMPYIntegerRing = _from_PythonIntegerRing
    _from_PythonFiniteField = _from_PythonIntegerRing
    _from_GMPYFiniteField = _from_PythonIntegerRing
    _from_GMPYIntegerModRing = _from_PythonIntegerRing
    _from_PythonIntegerModRing = _from_PythonIntegerRing

    def _from_PythonRationalField(self, a, K0):
        if a.denominator == 1:
            return self.dtype(a.numerator)
    _from_GMPYRationalField = _from_PythonRationalField

    def _from_RealField(self, a, K0):
        p, q = K0.to_rational(a)

        if q == 1:
            return self.dtype(p)

    def _from_AlgebraicField(self, a, K0):
        if a.is_ground:
            return self.convert(a.rep.LC, K0.domain)

    @abc.abstractmethod
    def finite_field(self, p):
        """Return a finite field."""

    @abc.abstractmethod
    def finite_ring(self, n):
        """Return a finite ring."""


class PythonIntegerRing(IntegerRing):
    """Integer ring based on Python's integers."""

    dtype = PythonInteger
    zero = dtype(0)
    one = dtype(1)

    def gcdex(self, a, b):
        """Compute extended GCD of ``a`` and ``b``."""
        return python_gcdex(a, b)

    def gcd(self, a, b):
        """Compute GCD of ``a`` and ``b``."""
        return python_gcd(a, b)

    def sqrt(self, a):
        """Compute square root of ``a``."""
        return python_sqrt(a)

    def factorial(self, a):
        """Compute factorial of ``a``."""
        return python_factorial(a)

    def finite_field(self, p):
        from .finitefield import PythonFiniteField
        return PythonFiniteField(p)

    def finite_ring(self, n):
        from .finitefield import PythonIntegerModRing
        return PythonIntegerModRing(n)


class GMPYIntegerRing(IntegerRing):
    """Integer ring based on GMPY's integers."""

    dtype = GMPYInteger
    zero = dtype(0)
    one = dtype(1)

    def gcdex(self, a, b):
        """Compute extended GCD of ``a`` and ``b``."""
        h, s, t = gmpy_gcdex(a, b)
        return s, t, h

    def gcd(self, a, b):
        """Compute GCD of ``a`` and ``b``."""
        return gmpy_gcd(a, b)

    def sqrt(self, a):
        """Compute square root of ``a``."""
        return gmpy_sqrt(a)

    def factorial(self, a):
        """Compute factorial of ``a``."""
        return gmpy_factorial(a)

    def finite_field(self, p):
        from .finitefield import GMPYFiniteField
        return GMPYFiniteField(p)

    def finite_ring(self, n):
        from .finitefield import GMPYIntegerModRing
        return GMPYIntegerModRing(n)


ZZ_python = PythonIntegerRing()
ZZ_gmpy = GMPYIntegerRing()
