"""Implementation of :class:`FractionField` class. """

from ..polys.polyerrors import CoercionFailed, GeneratorsError
from .compositedomain import CompositeDomain
from .field import Field


__all__ = ('FractionField',)


class FractionField(Field, CompositeDomain):
    """A class for representing multivariate rational function fields. """

    is_FractionField = is_Frac = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def __init__(self, domain_or_field, symbols=None, order=None):
        from ..polys.fields import FracField

        if isinstance(domain_or_field, FracField) and symbols is None and order is None:
            field = domain_or_field
        else:
            field = FracField(symbols, domain_or_field, order)

        self.field = field
        self.dtype = field.dtype

        self.gens = field.gens
        self.ngens = field.ngens
        self.symbols = field.symbols
        self.domain = field.domain

    def new(self, element):
        return self.field.field_new(element)

    @property
    def zero(self):
        return self.field.zero

    @property
    def one(self):
        return self.field.one

    @property
    def order(self):
        return self.field.order

    def __str__(self):
        return str(self.domain) + '(' + ','.join(map(str, self.symbols)) + ')'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.domain, self.symbols))

    def __eq__(self, other):
        """Returns `True` if two domains are equivalent. """
        return isinstance(other, FractionField) and \
            self.dtype == other.dtype and self.field == other.field

    def to_diofant(self, a):
        """Convert `a` to a Diofant object. """
        return a.as_expr()

    def from_diofant(self, a):
        """Convert Diofant's expression to `dtype`. """
        return self.field.from_expr(a)

    def from_ZZ_python(self, a, K0):
        """Convert a Python `int` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_QQ_python(self, a, K0):
        """Convert a Python `Fraction` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_ZZ_gmpy(self, a, K0):
        """Convert a GMPY `mpz` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_QQ_gmpy(self, a, K0):
        """Convert a GMPY `mpq` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_RealField(self, a, K0):
        """Convert a mpmath `mpf` object to `dtype`. """
        return self(self.domain.convert(a, K0))

    def from_PolynomialRing(self, a, K0):
        """Convert a polynomial to ``dtype``. """
        try:
            return self.new(a)
        except (CoercionFailed, GeneratorsError):
            return

    def from_FractionField(self, a, K0):
        """Convert a rational function to ``dtype``. """
        try:
            return a.set_field(self.field)
        except (CoercionFailed, GeneratorsError):
            return

    def get_ring(self):
        """Returns a field associated with `self`. """
        return self.field.to_ring().to_domain()

    def is_positive(self, a):
        """Returns True if `LC(a)` is positive. """
        return self.domain.is_positive(a.numer.LC)

    def is_negative(self, a):
        """Returns True if `LC(a)` is negative. """
        return self.domain.is_negative(a.numer.LC)

    def is_nonpositive(self, a):
        """Returns True if `LC(a)` is non-positive. """
        return self.domain.is_nonpositive(a.numer.LC)

    def is_nonnegative(self, a):
        """Returns True if `LC(a)` is non-negative. """
        return self.domain.is_nonnegative(a.numer.LC)

    def numer(self, a):
        """Returns numerator of ``a``. """
        return a.numer

    def denom(self, a):
        """Returns denominator of ``a``. """
        return a.denom

    def factorial(self, a):
        """Returns factorial of `a`. """
        return self.dtype(self.domain.factorial(a))
