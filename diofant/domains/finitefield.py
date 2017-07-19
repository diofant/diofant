"""Implementation of :class:`FiniteField` class. """

from ..polys.polyerrors import CoercionFailed
from .field import Field
from .groundtypes import DiofantInteger
from .modularinteger import ModularIntegerFactory
from .simpledomain import SimpleDomain


__all__ = ('FiniteField',)


class FiniteField(Field, SimpleDomain):
    """General class for finite fields. """

    rep = 'FF'

    is_FiniteField = is_FF = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    mod = None

    def __init__(self, mod, dom, symmetric=True):
        if mod <= 0:
            raise ValueError('modulus must be a positive integer, got %s' % mod)

        self.dtype = ModularIntegerFactory(mod, dom, symmetric, self)
        self.zero = self.dtype(0)
        self.one = self.dtype(1)
        self.domain = dom
        self.mod = mod

    def __str__(self):
        return 'GF(%s)' % self.mod

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.mod, self.domain))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, FiniteField) and \
            self.mod == other.mod and self.domain == other.domain

    def characteristic(self):
        """Return the characteristic of this domain. """
        return self.mod

    def get_field(self):
        """Returns a field associated with ``self``. """
        return self

    def to_diofant(self, a):
        """Convert ``a`` to a Diofant object. """
        return DiofantInteger(int(a))

    def from_diofant(self, a):
        """Convert Diofant's Integer to Diofant's ``Integer``. """
        if a.is_Integer:
            return self.dtype(self.domain.dtype(int(a)))
        elif a.is_Float and int(a) == a:
            return self.dtype(self.domain.dtype(int(a)))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def from_FF_python(self, a, K0=None):
        """Convert ``ModularInteger(int)`` to ``dtype``. """
        return self.dtype(self.domain.from_ZZ_python(a.val, K0.domain))

    def from_ZZ_python(self, a, K0=None):
        """Convert Python's ``int`` to ``dtype``. """
        return self.dtype(self.domain.from_ZZ_python(a, K0))

    def from_QQ_python(self, a, K0=None):
        """Convert Python's ``Fraction`` to ``dtype``. """
        if a.denominator == 1:
            return self.from_ZZ_python(a.numerator)

    def from_FF_gmpy(self, a, K0=None):
        """Convert ``ModularInteger(mpz)`` to ``dtype``. """
        return self.dtype(self.domain.from_ZZ_gmpy(a.val, K0.domain))

    def from_ZZ_gmpy(self, a, K0=None):
        """Convert GMPY's ``mpz`` to ``dtype``. """
        return self.dtype(self.domain.from_ZZ_gmpy(a, K0))

    def from_QQ_gmpy(self, a, K0=None):
        """Convert GMPY's ``mpq`` to ``dtype``. """
        if a.denominator == 1:
            return self.from_ZZ_gmpy(a.numerator)

    def from_RealField(self, a, K0):
        """Convert mpmath's ``mpf`` to ``dtype``. """
        p, q = K0.to_rational(a)

        if q == 1:
            return self.dtype(self.domain.dtype(p))
