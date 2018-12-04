"""Implementation of :class:`FiniteField` class. """

from ..core.cache import cacheit
from ..core.compatibility import DIOFANT_INTS
from ..ntheory import isprime, perfect_power
from ..polys.polyerrors import CoercionFailed
from .domainelement import DomainElement
from .field import Field
from .groundtypes import DiofantInteger
from .integerring import GMPYIntegerRing, PythonIntegerRing
from .simpledomain import SimpleDomain


__all__ = 'FiniteField', 'GMPYFiniteField', 'PythonFiniteField'


class FiniteField(Field, SimpleDomain):
    """General class for finite fields. """

    is_FiniteField = is_FF = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    mod = None

    @cacheit
    def __new__(cls, mod, dom):
        if not (isinstance(mod, DIOFANT_INTS) and isprime(mod)):
            pp = perfect_power(mod)
            if not pp:
                raise ValueError('modulus must be a positive prime number, got %s' % mod)
            mod, deg = pp
            raise NotImplementedError
        else:
            deg = 1

        mod = dom.convert(mod)
        order = mod**deg

        obj = super().__new__(cls)

        obj.dtype = type("ModularIntegerMod%s" % mod, (ModularInteger,),
                         {"_parent": obj, "mod": mod, "domain": dom})

        obj.domain = dom
        obj.mod = mod
        obj.order = order

        obj.rep = 'GF(%s)' % obj.order

        obj.zero = obj.dtype(0)
        obj.one = obj.dtype(1)

        return obj

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.order, self.domain))

    def __eq__(self, other):
        """Returns ``True`` if two domains are equivalent. """
        return isinstance(other, FiniteField) and \
            self.order == other.order and self.domain == other.domain

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
        return self.dtype(self.domain.convert(a.rep, K0.domain))

    def _from_PythonIntegerRing(self, a, K0=None):
        return self.dtype(self.domain.convert(a, K0))

    def _from_PythonRationalField(self, a, K0=None):
        if a.denominator == 1:
            return self.convert(a.numerator)

    def _from_GMPYFiniteField(self, a, K0=None):
        return self.dtype(self.domain.convert(a.rep, K0.domain))

    def _from_GMPYIntegerRing(self, a, K0=None):
        return self.dtype(self.domain.convert(a, K0))

    def _from_GMPYRationalField(self, a, K0=None):
        if a.denominator == 1:
            return self.convert(a.numerator)

    def _from_RealField(self, a, K0):
        p, q = K0.to_rational(a)

        if q == 1:
            return self.dtype(self.domain.dtype(p))

    def is_positive(self, a):
        """Returns True if ``a`` is positive. """
        return a.rep != 0

    def is_negative(self, a):
        """Returns True if ``a`` is negative. """
        return False


class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    def __new__(cls, mod):
        return super().__new__(cls, mod, PythonIntegerRing())


class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY's integers. """

    def __new__(cls, mod):
        return super().__new__(cls, mod, GMPYIntegerRing())


class ModularInteger(DomainElement):
    """A class representing a modular integer. """

    mod, domain, _parent = None, None, None

    @property
    def parent(self):
        return self._parent

    def __init__(self, rep):
        if isinstance(rep, self.__class__):
            self.rep = rep.rep % self.mod
        else:
            self.rep = self.domain.convert(rep) % self.mod

    def __hash__(self):
        return hash((self.rep, self.mod))

    def __int__(self):
        return int(self.rep)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.rep)

    def __add__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.__class__(self.rep + other.rep)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.__class__(self.rep - other.rep)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.__class__(self.rep * other.rep)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self * other**-1

    def __rtruediv__(self, other):
        return (self**-1).__mul__(other)

    def __mod__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.__class__(self.rep % other.rep)

    def __rmod__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return other.__mod__(self)

    def __pow__(self, exp):
        if not isinstance(exp, DIOFANT_INTS):
            raise TypeError("Integer exponent expected, got %s" % type(exp))
        if exp < 0:
            rep, exp = self.domain.invert(self.rep, self.mod), -exp
        else:
            rep = self.rep
        return self.__class__(rep**exp)

    def __eq__(self, other):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return self.rep == other.rep

    def __bool__(self):
        return bool(self.rep)
