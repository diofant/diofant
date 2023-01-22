"""Implementation of :class:`RationalField` class."""

from ..polys.polyerrors import CoercionFailedError
from .characteristiczero import CharacteristicZero
from .field import Field
from .groundtypes import DiofantRational, GMPYRational, PythonRational
from .simpledomain import SimpleDomain


class RationalField(CharacteristicZero, SimpleDomain, Field):
    """General class for rational fields."""

    rep = 'QQ'

    is_RationalField = True
    is_Numerical = True

    has_assoc_Ring = True

    def algebraic_field(self, *extension):
        r"""Return an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`."""
        from . import AlgebraicField
        return AlgebraicField(self, *extension)

    def to_expr(self, element):
        return DiofantRational(element.numerator, element.denominator)

    def from_expr(self, expr):
        if expr.is_Rational:
            return self.dtype(expr.numerator, expr.denominator)
        if expr.is_Float:
            from . import RR
            return self.dtype(*RR.to_rational(expr))
        raise CoercionFailedError(f'expected `Rational` object, got {expr}')

    def _from_PythonIntegerRing(self, a, K0):
        return self.dtype(a)
    _from_GMPYIntegerRing = _from_PythonIntegerRing

    def _from_PythonRationalField(self, a, K0):
        return self.dtype(a.numerator, a.denominator)
    _from_GMPYRationalField = _from_PythonRationalField

    def _from_RealField(self, a, K0):
        return self.dtype(*K0.to_rational(a))

    def _from_AlgebraicField(self, a, K0):
        if a.is_ground:
            return self.convert(a.rep.LC, K0.domain)


class PythonRationalField(RationalField):
    """Rational field based on Python's rationals."""

    dtype = PythonRational
    zero = dtype(0)
    one = dtype(1)

    @property
    def ring(self):
        from .integerring import PythonIntegerRing
        return PythonIntegerRing()


class GMPYRationalField(RationalField):
    """Rational field based on GMPY's rationals."""

    dtype = GMPYRational
    zero = dtype(0)
    one = dtype(1)

    @property
    def ring(self):
        from .integerring import GMPYIntegerRing
        return GMPYIntegerRing()


QQ_python = PythonRationalField()
QQ_gmpy = GMPYRationalField()
