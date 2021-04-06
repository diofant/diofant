"""Implementation of :class:`ComplexField` class."""

from __future__ import annotations

import mpmath

from ..core import Float, I
from ..polys.polyerrors import CoercionFailed
from .characteristiczero import CharacteristicZero
from .field import Field
from .mpelements import MPContext
from .simpledomain import SimpleDomain


class ComplexField(CharacteristicZero, SimpleDomain, Field):
    """Complex numbers up to the given precision."""

    rep = 'CC'

    is_ComplexField = True

    is_Exact = False
    is_Numerical = True

    _default_precision = 53

    @property
    def has_default_precision(self):
        return self.precision == self._default_precision

    @property
    def precision(self):
        return self._context.prec

    @property
    def dps(self):
        return self._context.dps

    @property
    def tolerance(self):
        return self._context.tolerance

    def __new__(cls, prec=_default_precision, dps=None, tol=None):
        context = MPContext(prec, dps, tol)

        obj = super().__new__(cls)

        try:
            obj.dtype = _complexes_cache[(context.prec, context.tolerance)]
        except KeyError:
            _complexes_cache[(context.prec, context.tolerance)] = obj.dtype = context.mpc

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
        return (isinstance(other, ComplexField)
                and self.precision == other.precision
                and self.tolerance == other.tolerance)

    def __hash__(self):
        return self._hash

    def to_expr(self, element):
        return Float(element.real, self.dps) + I*Float(element.imag, self.dps)

    def from_expr(self, expr):
        number = expr.evalf(self.dps)
        real, imag = number.as_real_imag()

        if real.is_Number and imag.is_Number:
            return self.dtype(real, imag)
        else:
            raise CoercionFailed(f'expected complex number, got {expr}')

    def _from_PythonIntegerRing(self, element, base):
        return self.dtype(element)
    _from_GMPYIntegerRing = _from_PythonIntegerRing
    _from_RealField = _from_PythonIntegerRing
    _from_ComplexField = _from_PythonIntegerRing

    def _from_PythonRationalField(self, element, base):
        return self.dtype(element.numerator) / element.denominator
    _from_GMPYRationalField = _from_PythonRationalField

    def _from_AlgebraicField(self, element, base):
        return self.from_expr(base.to_expr(element))

    def get_exact(self):
        from . import QQ
        return QQ.algebraic_field(I)

    def gcd(self, a, b):
        return self.one

    def almosteq(self, a, b, tolerance=None):
        """Check if ``a`` and ``b`` are almost equal."""
        return self._context.almosteq(a, b, tolerance)

    def is_normal(self, element):
        return True


_complexes_cache: dict[tuple, ComplexField] = {}


CC = ComplexField()
