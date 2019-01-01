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

    def __init__(self, val):
        if isinstance(val, self.__class__):
            self.val = val.val % self.mod
        else:
            self.val = self.domain.convert(val) % self.mod

    def __hash__(self):
        return hash((self.val, self.mod))

    def __str__(self):
        return "%s mod %s" % (self.val, self.mod)

    def __int__(self):
        return int(self.val)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.val)

    @classmethod
    def _get_val(cls, other):
        if isinstance(other, cls):
            return other.val
        else:
            try:
                return cls.domain.convert(other)
            except CoercionFailed:
                return

    def __add__(self, other):
        val = self._get_val(other)

        if val is not None:
            return self.__class__(self.val + val)
        else:
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        val = self._get_val(other)

        if val is not None:
            return self.__class__(self.val - val)
        else:
            return NotImplemented

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        val = self._get_val(other)

        if val is not None:
            return self.__class__(self.val * val)
        else:
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        val = self._get_val(other)

        if val is not None:
            return self.__class__(self.val * self._invert(val))
        else:
            return NotImplemented

    def __rtruediv__(self, other):
        return self.invert().__mul__(other)

    def __mod__(self, other):
        val = self._get_val(other)

        if val is not None:
            return self.__class__(self.val % val)
        else:
            return NotImplemented

    def __rmod__(self, other):
        val = self._get_val(other)

        if val is not None:
            return self.__class__(val % self.val)
        else:
            return NotImplemented

    def __pow__(self, exp):
        if not exp:
            return self.__class__(self.domain.one)

        if exp < 0:
            val, exp = self.invert(), -exp
        else:
            val = self.val

        return self.__class__(val**exp)

    def _compare(self, other, op):
        val = self._get_val(other)

        if val is not None:
            return op(self.val, val % self.mod)
        else:
            return NotImplemented

    def __eq__(self, other):
        return self._compare(other, operator.eq)

    def __lt__(self, other):
        return self._compare(other, operator.lt)

    def __bool__(self):
        return bool(self.val)

    @classmethod
    def _invert(cls, value):
        return cls.domain.invert(value, cls.mod)

    def invert(self):
        return self.__class__(self._invert(self.val))
