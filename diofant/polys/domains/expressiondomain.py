"""Implementation of :class:`ExpressionDomain` class. """

from diofant.polys.domains.field import Field
from diofant.polys.domains.simpledomain import SimpleDomain
from diofant.polys.domains.characteristiczero import CharacteristicZero
from diofant.core import sympify, SympifyError
from diofant.utilities import public


@public
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

        def __repr__(self):
            return 'EX(%s)' % repr(self.ex)

        def __str__(self):
            return 'EX(%s)' % str(self.ex)

        def __hash__(self):
            return hash((self.__class__.__name__, self.ex))

        def as_expr(self):
            return self.ex

        def numer(self):
            return self.__class__(self.ex.as_numer_denom()[0])

        def denom(self):
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

        __div__ = __truediv__
        __rdiv__ = __rtruediv__

        def __eq__(self, other):
            return self.ex == self.__class__(other).ex

        def __ne__(self, other):
            return not self.__eq__(other)

        def __bool__(self):
            return self.ex != 0

        def gcd(self, other):
            from diofant.polys import gcd
            return self.__class__(gcd(self.ex, self.__class__(other).ex))

        def lcm(self, other):
            from diofant.polys import lcm
            return self.__class__(lcm(self.ex, self.__class__(other).ex))

    dtype = Expression

    zero = Expression(0)
    one = Expression(1)

    rep = 'EX'

    has_assoc_Ring = False
    has_assoc_Field = True

    def __init__(self):
        pass

    def to_diofant(self, a):
        """Convert ``a`` to a Diofant object. """
        return a.as_expr()

    def from_diofant(self, a):
        """Convert Diofant's expression to ``dtype``. """
        return self.dtype(a)

    def from_ZZ_python(self, a, K0):
        """Convert a Python ``int`` object to ``dtype``. """
        return self(K0.to_diofant(a))

    def from_QQ_python(self, a, K0):
        """Convert a Python ``Fraction`` object to ``dtype``. """
        return self(K0.to_diofant(a))

    def from_ZZ_gmpy(self, a, K0):
        """Convert a GMPY ``mpz`` object to ``dtype``. """
        return self(K0.to_diofant(a))

    def from_QQ_gmpy(self, a, K0):
        """Convert a GMPY ``mpq`` object to ``dtype``. """
        return self(K0.to_diofant(a))

    def from_RealField(self, a, K0):
        """Convert a mpmath ``mpf`` object to ``dtype``. """
        return self(K0.to_diofant(a))

    def from_PolynomialRing(self, a, K0):
        """Convert a ``DMP`` object to ``dtype``. """
        return self(K0.to_diofant(a))

    def from_FractionField(self, a, K0):
        """Convert a ``DMF`` object to ``dtype``. """
        return self(K0.to_diofant(a))

    def from_ExpressionDomain(self, a, K0):
        """Convert a ``EX`` object to ``dtype``. """
        return a

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        return self  # XXX: EX is not a ring but we don't have much choice here.

    def get_field(self):
        """Returns a field associated with ``self``. """
        return self

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

    def numer(self, a):
        """Returns numerator of ``a``. """
        return a.numer()

    def denom(self, a):
        """Returns denominator of ``a``. """
        return a.denom()

    def gcd(self, a, b):
        return a.gcd(b)

    def lcm(self, a, b):
        return a.lcm(b)
