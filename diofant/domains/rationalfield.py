"""Implementation of :class:`RationalField` class."""

from ..polys.polyerrors import CoercionFailed
from .characteristiczero import CharacteristicZero
from .field import Field
from .groundtypes import (DiofantRational, GMPYRational, PythonRational,
                          gmpy_qdiv)
from .simpledomain import SimpleDomain


__all__ = 'GMPYRationalField', 'PythonRationalField', 'RationalField'


class RationalField(CharacteristicZero, SimpleDomain, Field):
    """General class for rational fields."""

    rep = 'QQ'

    is_RationalField = is_QQ = True
    is_Numerical = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`."""
        from . import AlgebraicField
        return AlgebraicField(self, *extension)

    def to_expr(self, a):
        return DiofantRational(a.numerator, a.denominator)

    def from_expr(self, a):
        if a.is_Rational:
            return self.dtype(a.numerator, a.denominator)
        elif a.is_Float:
            from . import RR
            return self.dtype(*RR.to_rational(a))
        else:
            raise CoercionFailed("expected `Rational` object, got %s" % a)

    def _from_PythonIntegerRing(self, a, K0):
        return self.dtype(a)

    def _from_PythonRationalField(self, a, K0):
        return self.dtype(a.numerator, a.denominator)

    def _from_GMPYIntegerRing(self, a, K0):
        return self.dtype(a)

    def _from_GMPYRationalField(self, a, K0):
        return self.dtype(a.numerator, a.denominator)

    def _from_RealField(self, a, K0):
        return self.dtype(*K0.to_rational(a))

    def _from_AlgebraicField(self, a, K0):
        if a.is_ground:
            return self.convert(a.LC(), K0.domain)


class PythonRationalField(RationalField):
    """Rational field based on Python's rationals."""

    dtype = PythonRational
    zero = dtype(0)
    one = dtype(1)

    @property
    def ring(self):
        """Returns ring associated with ``self``."""
        from .integerring import PythonIntegerRing
        return PythonIntegerRing()


class GMPYRationalField(RationalField):
    """Rational field based on GMPY's rationals."""

    dtype = GMPYRational
    zero = dtype(0)
    one = dtype(1)

    @property
    def ring(self):
        """Returns ring associated with ``self``."""
        from .integerring import GMPYIntegerRing
        return GMPYIntegerRing()

    def exquo(self, a, b):
        """Exact quotient of `a` and `b`, implies `__truediv__`."""
        return self.dtype(gmpy_qdiv(a, b))

    def quo(self, a, b):
        """Quotient of `a` and `b`, implies `__truediv__`."""
        return self.dtype(gmpy_qdiv(a, b))


QQ_python = PythonRationalField()
QQ_gmpy = GMPYRationalField()
