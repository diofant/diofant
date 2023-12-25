"""Implementation of :class:`RealField` class."""

from __future__ import annotations

import mpmath

from ..core import Float
from ..polys.polyerrors import CoercionFailedError
from .characteristiczero import CharacteristicZero
from .field import Field
from .mpelements import MPContext
from .simpledomain import SimpleDomain


class RealField(CharacteristicZero, SimpleDomain, Field):
    """Real numbers up to the given precision."""

    rep = 'RR'

    is_RealField = True

    is_Exact = False
    is_Numerical = True

    @property
    def precision(self):
        return self._context.prec

    @property
    def dps(self):
        return self._context.dps

    @property
    def tolerance(self):
        return self._context.tolerance

    def __new__(cls, prec=53, dps=None, tol=None):
        context = MPContext(prec, dps, tol)

        obj = super().__new__(cls)

        try:
            obj.dtype = _reals_cache[(context.prec, context.tolerance)]
        except KeyError:
            _reals_cache[(context.prec, context.tolerance)] = obj.dtype = context.mpf

        context._parent = obj
        obj._context = context
        obj._hash = hash((cls.__name__, obj.dtype, context.prec, context.tolerance))

        obj.zero = obj.dtype(0)
        obj.one = obj.dtype(1)

        return obj

    def __getnewargs_ex__(self):
        return (), {'prec': self.precision,
                    'tol': mpmath.mpf(self.tolerance._mpf_)}

    def __eq__(self, other):
        return (isinstance(other, RealField)
                and self.precision == other.precision
                and self.tolerance == other.tolerance)

    def __hash__(self):
        return self._hash

    def to_expr(self, element):
        return Float(element, self.dps)

    def from_expr(self, expr):
        number = expr.evalf(self.dps)

        if number.is_Number:
            return self.dtype(number)
        raise CoercionFailedError(f'expected real number, got {expr}')

    def _from_PythonIntegerRing(self, element, base):
        return self.dtype(element)
    _from_GMPYIntegerRing = _from_PythonIntegerRing

    def _from_PythonRationalField(self, element, base):
        return self.dtype(element.numerator) / element.denominator
    _from_GMPYRationalField = _from_PythonRationalField

    def _from_AlgebraicField(self, element, base):
        return self.from_expr(base.to_expr(element))

    def _from_RealField(self, element, base):
        if self == base:
            return element
        return self.dtype(element)

    def _from_ComplexField(self, element, base):
        if not element.imag:
            return self.dtype(element.real)

    def to_rational(self, element, limit=True):
        """Convert a real number to rational number."""
        return self._context.to_rational(element, limit)

    def get_exact(self):
        from . import QQ
        return QQ

    def gcd(self, a, b):
        return self.one

    def almosteq(self, a, b, tolerance=None):
        """Check if ``a`` and ``b`` are almost equal."""
        return self._context.almosteq(a, b, tolerance)


_reals_cache: dict[tuple, RealField] = {}


RR = RealField()
