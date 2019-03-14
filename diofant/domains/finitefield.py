"""Implementation of :class:`FiniteField` class. """

import numbers
import random

from ..core import Dummy, integer_digits
from ..ntheory import isprime, perfect_power
from ..polys.galoistools import gf_irreducible
from ..polys.polyerrors import CoercionFailed
from .field import Field
from .groundtypes import DiofantInteger
from .integerring import GMPYIntegerRing, PythonIntegerRing
from .quotientring import QuotientRingElement
from .simpledomain import SimpleDomain


__all__ = 'FiniteField', 'GMPYFiniteField', 'PythonFiniteField'


_modular_integer_cache = {}


class FiniteField(Field, SimpleDomain):
    """General class for finite fields. """

    is_FiniteField = is_FF = True
    is_Numerical = True

    has_assoc_Ring = False
    has_assoc_Field = True

    def __new__(cls, order, dom, modulus=None):
        if not (isinstance(order, numbers.Integral) and isprime(order)):
            pp = perfect_power(order)
            if not pp:
                raise ValueError('order must be a prime power, '
                                 'got %s' % order)
            mod, deg = pp
        else:
            mod, deg = order, 1
            if modulus is None:
                modulus = [1, 0]
            else:
                deg = len(modulus) - 1
                order = mod**deg

        if modulus is None:
            random.seed(0)
            modulus = gf_irreducible(deg, mod, dom)
        elif deg != len(modulus) - 1:
            raise ValueError('degree of a defining polynomial for the field'
                             ' does not match extension degree')
        modulus = tuple(map(dom.dtype, modulus))

        mod = dom.convert(mod)

        key = order, dom, mod, modulus

        obj = super().__new__(cls)

        obj.domain = dom
        obj.mod = mod
        obj.order = order

        obj.rep = 'GF(%s)' % obj.order

        try:
            obj.dtype = _modular_integer_cache[key]
        except KeyError:
            if deg == 1:
                obj.dtype = type("ModularInteger", (ModularInteger,),
                                 {"mod": mod, "domain": dom, "_parent": obj})
            else:
                ff = dom.finite_field(mod).poly_ring(Dummy('x'))
                mod = ff.from_dense(modulus)
                if not mod.is_irreducible:
                    raise ValueError('defining polynomial must be irreducible')
                obj.dtype = type("GaloisFieldElement", (GaloisFieldElement,),
                                 {"mod": mod, "domain": ff, "_parent": obj})
            _modular_integer_cache[key] = obj.dtype

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

    def __new__(cls, order, modulus=None):
        return super().__new__(cls, order, PythonIntegerRing(), modulus)


class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY's integers. """

    def __new__(cls, order, modulus=None):
        return super().__new__(cls, order, GMPYIntegerRing(), modulus)


class ModularInteger(QuotientRingElement):
    """A class representing a modular integer. """


class GaloisFieldElement(ModularInteger):
    def __init__(self, rep):
        if isinstance(rep, numbers.Integral):
            rep = rep % self.parent.order
            rep = integer_digits(rep, self.parent.mod)

        if isinstance(rep, (list, tuple)):
            rep = self.domain.from_dense(rep)
            self.rep = rep = rep % self.mod
        else:
            super().__init__(rep)

        self._int_rep = self.parent.domain.poly_ring(*self.rep.parent.symbols)(dict(self.rep))

    def __int__(self):
        return int(self._int_rep.eval(0, self.parent.mod))
