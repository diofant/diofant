"""Implementation of :class:`FiniteField` class. """

from ..ntheory import isprime, perfect_power
from ..polys.polyerrors import CoercionFailed
from .field import Field
from .groundtypes import DiofantInteger
from .integerring import GMPYIntegerRing, PythonIntegerRing
from .modularinteger import ModularIntegerFactory
from .simpledomain import SimpleDomain


__all__ = ('FiniteField', 'GMPYFiniteField', 'PythonFiniteField')


class FiniteField(Field, SimpleDomain):
    """General class for finite fields. """

    is_FiniteField = is_FF = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    mod = None

    def __init__(self, mod, dom, symmetric=True):
        if not isprime(mod):
            if perfect_power(mod):  # pragma: no cover
                raise NotImplementedError
            raise ValueError('modulus must be a positive prime number, got %s' % mod)

        self.dtype = ModularIntegerFactory(mod, dom, symmetric, self)
        self.zero = self.dtype(0)
        self.one = self.dtype(1)
        self.domain = dom
        self.mod = mod

        self.rep = 'GF(%s)' % self.mod

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.mod, self.domain))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, FiniteField) and \
            self.mod == other.mod and self.domain == other.domain

    @property
    def characteristic(self):
        """Return the characteristic of this domain. """
        return self.mod

    def to_expr(self, a):
        """Convert ``a`` to a Diofant object. """
        return DiofantInteger(int(a))

    def from_expr(self, a):
        """Convert Diofant's Integer to ``dtype``. """
        if a.is_Integer:
            return self.dtype(self.domain.dtype(int(a)))
        elif a.is_Float and int(a) == a:
            return self.dtype(self.domain.dtype(int(a)))
        else:
            raise CoercionFailed("expected an integer, got %s" % a)

    def _from_PythonFiniteField(self, a, K0=None):
        return self.dtype(self.domain.convert(a.val, K0.domain))

    def _from_PythonIntegerRing(self, a, K0=None):
        return self.dtype(self.domain.convert(a, K0))

    def _from_PythonRationalField(self, a, K0=None):
        if a.denominator == 1:
            return self.convert(a.numerator)

    def _from_GMPYFiniteField(self, a, K0=None):
        return self.dtype(self.domain.convert(a.val, K0.domain))

    def _from_GMPYIntegerRing(self, a, K0=None):
        return self.dtype(self.domain.convert(a, K0))

    def _from_GMPYRationalField(self, a, K0=None):
        if a.denominator == 1:
            return self.convert(a.numerator)

    def _from_RealField(self, a, K0):
        p, q = K0.to_rational(a)

        if q == 1:
            return self.dtype(self.domain.dtype(p))


class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    def __init__(self, mod, symmetric=True):
        super().__init__(mod, PythonIntegerRing(), symmetric)


class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY's integers. """

    def __init__(self, mod, symmetric=True):
        super().__init__(mod, GMPYIntegerRing(), symmetric)
