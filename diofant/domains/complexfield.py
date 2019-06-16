"""Implementation of :class:`ComplexField` class. """

from ..core import Float, I
from ..polys.polyerrors import CoercionFailed, DomainError
from .characteristiczero import CharacteristicZero
from .field import Field
from .mpelements import MPContext
from .simpledomain import SimpleDomain


__all__ = 'ComplexField',


_complexes_cache = {}


class ComplexField(Field, CharacteristicZero, SimpleDomain):
    """Complex numbers up to the given precision. """

    rep = 'CC'

    is_ComplexField = is_CC = True

    is_Exact = False
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

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

    def __eq__(self, other):
        return (isinstance(other, ComplexField)
                and self.precision == other.precision
                and self.tolerance == other.tolerance)

    def __hash__(self):
        return self._hash

    def to_expr(self, element):
        """Convert ``element`` to Diofant number. """
        return Float(element.real, self.dps) + I*Float(element.imag, self.dps)

    def from_expr(self, expr):
        """Convert Diofant's number to ``dtype``. """
        number = expr.evalf(self.dps)
        real, imag = number.as_real_imag()

        if real.is_Number and imag.is_Number:
            return self.dtype(real, imag)
        else:
            raise CoercionFailed("expected complex number, got %s" % expr)

    def _from_PythonIntegerRing(self, element, base):
        return self.dtype(element)

    def _from_PythonRationalField(self, element, base):
        return self.dtype(element.numerator) / element.denominator

    def _from_GMPYIntegerRing(self, element, base):
        return self.dtype(int(element))

    def _from_GMPYRationalField(self, element, base):
        return self.dtype(int(element.numerator)) / int(element.denominator)

    def _from_AlgebraicField(self, element, base):
        return self.from_expr(base.to_expr(element))

    def _from_RealField(self, element, base):
        return self.dtype(element)

    def _from_ComplexField(self, element, base):
        return self.dtype(element)

    def get_exact(self):
        """Returns an exact domain associated with ``self``. """
        raise DomainError("there is no exact domain associated with %s" % self)

    def gcd(self, a, b):
        """Returns GCD of ``a`` and ``b``. """
        return self.one

    def lcm(self, a, b):
        """Returns LCM of ``a`` and ``b``. """
        return a*b

    def almosteq(self, a, b, tolerance=None):
        """Check if ``a`` and ``b`` are almost equal. """
        return self._context.almosteq(a, b, tolerance)


CC = ComplexField()
