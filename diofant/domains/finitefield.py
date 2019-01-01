"""Implementation of :class:`FiniteField` class. """

import functools
import numbers
import operator

from ..ntheory import isprime, perfect_power
from ..polys.polyerrors import CoercionFailed
from .domainelement import DomainElement
from .field import Field
from .groundtypes import DiofantInteger
from .integerring import GMPYIntegerRing, PythonIntegerRing
from .simpledomain import SimpleDomain


__all__ = 'FiniteField', 'GMPYFiniteField', 'PythonFiniteField'


_modular_integer_cache = {}


class FiniteField(Field, SimpleDomain):
    """General class for finite fields. """

    is_FiniteField = is_FF = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    mod = None

    def __new__(cls, mod, dom):
        if not (isinstance(mod, numbers.Integral) and isprime(mod)):
            if perfect_power(mod):  # pragma: no cover
                raise NotImplementedError
            raise ValueError('modulus must be a positive prime number, got %s' % mod)

        mod = dom.convert(mod)
        key = mod, dom

        obj = super().__new__(cls)

        try:
            obj.dtype = _modular_integer_cache[key]
        except KeyError:
            obj.dtype = type("ModularIntegerMod%s" % mod, (ModularInteger,),
                             {"mod": mod, "domain": dom, "_parent": obj})
            _modular_integer_cache[key] = obj.dtype

        obj.zero = obj.dtype(0)
        obj.one = obj.dtype(1)
        obj.domain = dom
        obj.mod = mod

        obj.rep = 'GF(%s)' % obj.mod

        return obj

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


class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers. """

    def __new__(cls, mod):
        return super().__new__(cls, mod, PythonIntegerRing())


class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY's integers. """

    def __new__(cls, mod):
        return super().__new__(cls, mod, GMPYIntegerRing())


@functools.total_ordering
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

    def __str__(self):
        return "%s mod %s" % (self.rep, self.mod)

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
        return self.__class__(self.rep * self._invert(other.rep))

    def __rtruediv__(self, other):
        return self.invert().__mul__(other)

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
        if not exp:
            return self.parent.one

        if exp < 0:
            rep, exp = self.invert(), -exp
        else:
            rep = self.rep

        return self.__class__(rep**exp)

    def _compare(self, other, op):
        try:
            other = self.parent.convert(other)
        except CoercionFailed:
            return NotImplemented
        return op(self.rep, other.rep)

    def __eq__(self, other):
        return self._compare(other, operator.eq)

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __bool__(self):
        return bool(self.rep)

    @classmethod
    def _invert(cls, value):
        return cls.domain.invert(value, cls.mod)

    def invert(self):
        return self.__class__(self._invert(self.rep))
