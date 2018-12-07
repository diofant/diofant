"""Implementation of :class:`IntegerRing` class. """

import abc

from ..polys.polyerrors import CoercionFailed
from .characteristiczero import CharacteristicZero
from .groundtypes import (DiofantInteger, GMPYInteger, PythonInteger,
                          gmpy_factorial, gmpy_gcd, gmpy_gcdex, gmpy_lcm,
                          gmpy_sqrt, python_factorial, python_gcd,
                          python_gcdex, python_lcm, python_sqrt)
from .ring import Ring
from .simpledomain import SimpleDomain


__all__ = ('GMPYIntegerRing', 'IntegerRing', 'PythonIntegerRing')


class IntegerRing(Ring, CharacteristicZero, SimpleDomain):
    """General class for integer rings. """

    rep = 'ZZ'

    is_IntegerRing = is_ZZ = True
    is_Numerical = True

    has_assoc_Ring = True
    has_assoc_Field = True

    @property
    def field(self):
        """Returns a field associated with ``self``. """
        from . import QQ
        return QQ

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return self.field.algebraic_field(*extension)

    def to_expr(self, a):
        """Convert ``a`` to a Diofant object. """
        return DiofantInteger(a)

    def from_expr(self, a):
        """Convert Diofant's Integer to ``dtype``. """
        if a.is_Integer:
            return self.dtype(a.numerator)
        elif a.is_Float and int(a) == a:
            return self.dtype(int(a))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def _from_PythonIntegerRing(self, a, K0):
        return self.dtype(a)

    def _from_PythonFiniteField(self, a, K0):
        return self.dtype(a.to_int())

    def _from_PythonRationalField(self, a, K0):
        if a.denominator == 1:
            return self.dtype(a.numerator)

    def _from_GMPYFiniteField(self, a, K0):
        return self.dtype(a.to_int())

    def _from_GMPYIntegerRing(self, a, K0):
        return self.dtype(a)

    def _from_GMPYRationalField(self, a, K0):
        if a.denominator == 1:
            return self.dtype(a.numerator)

    def _from_RealField(self, a, K0):
        p, q = K0.to_rational(a)

        if q == 1:
            return self.dtype(p)

    def _from_AlgebraicField(self, a, K0):
        if a.is_ground:
            return self.convert(a.LC(), K0.domain)

    @abc.abstractmethod
    def finite_field(self, p):
        """Returns a finite field. """
        raise NotImplementedError


class PythonIntegerRing(IntegerRing):
    """Integer ring based on Python's integers. """

    dtype = PythonInteger
    zero = dtype(0)
    one = dtype(1)

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

    def finite_field(self, p):
        """Returns a finite field. """
        from .finitefield import PythonFiniteField
        return PythonFiniteField(p, self)


class GMPYIntegerRing(IntegerRing):
    """Integer ring based on GMPY's integers. """

    dtype = GMPYInteger
    zero = dtype(0)
    one = dtype(1)

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

    def finite_field(self, p):
        """Returns a finite field. """
        from .finitefield import GMPYFiniteField
        return GMPYFiniteField(p, self)


ZZ_python = PythonIntegerRing()
ZZ_gmpy = GMPYIntegerRing()
