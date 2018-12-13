"""Implementation of :class:`ExpressionDomain` class. """

from ..core import SympifyError, sympify
from .characteristiczero import CharacteristicZero
from .field import Field
from .simpledomain import SimpleDomain


__all__ = ('ExpressionDomain',)


class ExpressionDomain(Field, CharacteristicZero, SimpleDomain):
    """A class for arbitrary expressions. """

    is_SymbolicDomain = is_EX = True

    class Expression:
        """An arbitrary expression. """

        def __init__(self, ex):
            if not isinstance(ex, self.__class__):
                self.ex = sympify(ex)
            else:
                self.ex = ex.ex

        def __str__(self):
            return 'EX(%s)' % str(self.ex)

        def __hash__(self):
            return hash((self.__class__.__name__, self.ex))

        def as_expr(self):
            return self.ex

        @property
        def numerator(self):
            return self.__class__(self.ex.as_numer_denom()[0])

        @property
        def denominator(self):
            return self.__class__(self.ex.as_numer_denom()[1])

        def simplify(self, ex):
            return self.__class__(ex.cancel())

        def __abs__(self):
            return self.__class__(abs(self.ex))

        def __neg__(self):
            return self.__class__(-self.ex)

        def _to_ex(self, other):
            try:
                return self.__class__(other)
            except SympifyError:
                return

        def __add__(self, other):
            other = self._to_ex(other)

            if other is not None:
                return self.simplify(self.ex + other.ex)
            else:
                return NotImplemented

        def __radd__(self, other):
            return self.simplify(self.__class__(other).ex + self.ex)

        def __sub__(self, other):
            other = self._to_ex(other)

            if other is not None:
                return self.simplify(self.ex - other.ex)
            else:
                return NotImplemented

        def __rsub__(self, other):
            return self.simplify(self.__class__(other).ex - self.ex)

        def __mul__(self, other):
            other = self._to_ex(other)

            if other is not None:
                return self.simplify(self.ex*other.ex)
            else:
                return NotImplemented

        def __rmul__(self, other):
            return self.simplify(self.__class__(other).ex*self.ex)

        def __pow__(self, n):
            n = self._to_ex(n)

            if n is not None:
                return self.simplify(self.ex**n.ex)
            else:
                return NotImplemented

        def __truediv__(self, other):
            other = self._to_ex(other)

            if other is not None:
                return self.simplify(self.ex/other.ex)
            else:
                return NotImplemented

        def __rtruediv__(self, other):
            return self.simplify(self.__class__(other).ex/self.ex)

        def __eq__(self, other):
            return self.ex == self.__class__(other).ex

        def __bool__(self):
            return self.ex != 0

        def gcd(self, other):
            from ..polys import gcd
            return self.__class__(gcd(self.ex, self.__class__(other).ex))

        def lcm(self, other):
            from ..polys import lcm
            return self.__class__(lcm(self.ex, self.__class__(other).ex))

    dtype = Expression

    zero = Expression(0)
    one = Expression(1)

    rep = 'EX'

    has_assoc_Ring = False
    has_assoc_Field = True

    def to_expr(self, a):
        """Convert ``a`` to a Diofant object. """
        return a.as_expr()

    def from_expr(self, a):
        """Convert Diofant's expression to ``dtype``. """
        return self.dtype(a)

    def _from_PythonIntegerRing(self, a, K0):
        return self(K0.to_expr(a))

    def _from_PythonRationalField(self, a, K0):
        return self(K0.to_expr(a))

    def _from_GMPYIntegerRing(self, a, K0):
        return self(K0.to_expr(a))

    def _from_GMPYRationalField(self, a, K0):
        return self(K0.to_expr(a))

    def _from_RealField(self, a, K0):
        return self(K0.to_expr(a))

    def _from_PolynomialRing(self, a, K0):
        return self(K0.to_expr(a))

    def _from_FractionField(self, a, K0):
        return self(K0.to_expr(a))

    def _from_AlgebraicField(self, a, K0):
        return self(K0.to_expr(a))

    @property
    def ring(self):
        """Returns a ring associated with ``self``. """
        return self  # XXX: EX is not a ring but we don't have much choice here.

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return a.ex.as_coeff_mul()[0].is_positive

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return a.ex.as_coeff_mul()[0].is_negative

    def is_nonpositive(self, a):
        """Returns True if ``a`` is non-positive. """
        return a.ex.as_coeff_mul()[0].is_nonpositive

    def is_nonnegative(self, a):
        """Returns True if ``a`` is non-negative. """
        return a.ex.as_coeff_mul()[0].is_nonnegative

    def gcd(self, a, b):
        return a.gcd(b)

    def lcm(self, a, b):
        return a.lcm(b)


EX = ExpressionDomain()
