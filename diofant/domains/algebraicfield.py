"""Implementation of :class:`AlgebraicField` class. """

from ..polys.polyclasses import ANP
from ..polys.polyerrors import (CoercionFailed, DomainError, IsomorphismFailed,
                                NotAlgebraic)
from .characteristiczero import CharacteristicZero
from .field import Field
from .simpledomain import SimpleDomain


__all__ = ('AlgebraicField',)


class AlgebraicField(Field, CharacteristicZero, SimpleDomain):
    """A class for representing algebraic number fields. """

    dtype = ANP

    is_AlgebraicField = is_Algebraic = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    def __init__(self, dom, *ext):
        if not dom.is_QQ:
            raise DomainError("ground domain must be a rational field")

        from ..polys.numberfields import to_number_field

        self.ext = to_number_field(ext)
        self.mod = self.ext.minpoly.rep
        self.domain = dom

        self.ngens = 1
        self.symbols = self.gens = (self.ext.as_expr(),)
        self.unit = self([dom(1), dom(0)])

        self.zero = self.dtype.zero(self.mod.rep, dom)
        self.one = self.dtype.one(self.mod.rep, dom)

    def new(self, element):
        if isinstance(element, list):
            return self.dtype(element, self.mod.rep, self.domain)
        else:
            return self.convert(element)

    def __str__(self):
        return str(self.domain) + '<' + str(self.ext) + '>'

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.domain, self.ext))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, AlgebraicField) and \
            self.dtype == other.dtype and self.ext == other.ext

    def algebraic_field(self, *extension):
        r"""Returns an algebraic field, i.e. `\mathbb{Q}(\alpha, \ldots)`. """
        return AlgebraicField(self.domain, *((self.ext,) + extension))

    def to_diofant(self, a):
        """Convert ``a`` to a Diofant object. """
        from ..polys.numberfields import AlgebraicNumber
        return AlgebraicNumber(self.ext, a).as_expr()

    def from_diofant(self, a):
        """Convert Diofant's expression to ``dtype``. """
        try:
            return self([self.domain.from_diofant(a)])
        except CoercionFailed:
            pass

        from ..polys.numberfields import to_number_field

        try:
            return self(to_number_field(a, self.ext).native_coeffs())
        except (NotAlgebraic, IsomorphismFailed):
            raise CoercionFailed(
                "%s is not a valid algebraic number in %s" % (a, self))

    def from_ZZ_python(self, a, K0):
        """Convert a Python ``int`` object to ``dtype``. """
        return self([self.domain.convert(a, K0)])

    def from_QQ_python(self, a, K0):
        """Convert a Python ``Fraction`` object to ``dtype``. """
        return self([self.domain.convert(a, K0)])

    def from_ZZ_gmpy(self, a, K0):
        """Convert a GMPY ``mpz`` object to ``dtype``. """
        return self([self.domain.convert(a, K0)])

    def from_QQ_gmpy(self, a, K0):
        """Convert a GMPY ``mpq`` object to ``dtype``. """
        return self([self.domain.convert(a, K0)])

    def from_RealField(self, a, K0):
        """Convert a mpmath ``mpf`` object to ``dtype``. """
        return self([self.domain.convert(a, K0)])

    def from_AlgebraicField(self, a, K0):
        return self.from_diofant(K0.to_diofant(a))

    def get_ring(self):
        """Returns a ring associated with ``self``. """
        raise DomainError('there is no ring associated with %s' % self)

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return self.domain.is_positive(a.LC())

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return self.domain.is_negative(a.LC())

    def is_nonpositive(self, a):
        """Returns True if ``a`` is non-positive. """
        return self.domain.is_nonpositive(a.LC())

    def is_nonnegative(self, a):
        """Returns True if ``a`` is non-negative. """
        return self.domain.is_nonnegative(a.LC())

    def numer(self, a):
        """Returns numerator of ``a``. """
        return a

    def denom(self, a):
        """Returns denominator of ``a``. """
        return self.one
