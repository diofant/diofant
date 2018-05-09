"""Implementation of :class:`IntegerRing` class. """

import math

from ..polys.polyerrors import CoercionFailed
from .characteristiczero import CharacteristicZero
from .groundtypes import DiofantInteger
from .ring import Ring
from .simpledomain import SimpleDomain


__all__ = ('IntegerRing',)


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

    def to_diofant(self, a):
        """Convert ``a`` to a Diofant object. """
        return DiofantInteger(a)

    def from_diofant(self, a):
        """Convert Diofant's Integer to ``dtype``. """
        if a.is_Integer:
            return self.dtype(a.p)
        elif a.is_Float and int(a) == a:
            return self.dtype(int(a))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def from_PythonIntegerRing(self, a, K0):
        """Convert Python's ``int`` to ``dtype``. """
        return self.dtype(a)

    def from_PythonFiniteField(self, a, K0):
        """Convert ``ModularInteger(int)`` to ``dtype``. """
        return self.dtype(a.to_int())

    def from_PythonRationalField(self, a, K0):
        """Convert Python's ``Fraction`` to ``dtype``. """
        if a.denominator == 1:
            return self.dtype(a.numerator)

    def from_GMPYFiniteField(self, a, K0):
        """Convert ``ModularInteger(mpz)`` to ``dtype``. """
        return self.dtype(a.to_int())

    def from_GMPYIntegerRing(self, a, K0):
        """Convert GMPY's ``mpz`` to ``dtype``. """
        return self.dtype(a)

    def from_GMPYRationalField(self, a, K0):
        """Convert GMPY's ``mpq`` to ``dtype``. """
        if a.denominator == 1:
            return self.dtype(a.numerator)

    def from_RealField(self, a, K0):
        """Convert mpmath's ``mpf`` to ``dtype``. """
        p, q = K0.to_rational(a)

        if q == 1:
            return self.dtype(p)

    def from_AlgebraicField(self, a, K0):
        """Convert an algebraic number to ``dtype``. """
        if a.is_ground:
            return self.convert(a.LC(), K0.domain)

    def log(self, a, b):
        """Returns b-base logarithm of ``a``. """
        return self.dtype(math.log(int(a), b))
