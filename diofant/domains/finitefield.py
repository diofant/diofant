"""Implementation of :class:`FiniteField` class."""

from __future__ import annotations

import numbers
import random

from ..core import Dummy, integer_digits
from ..ntheory import factorint, is_primitive_root, isprime
from ..polys.polyerrors import CoercionFailedError
from .field import Field
from .groundtypes import DiofantInteger
from .integerring import GMPYIntegerRing, PythonIntegerRing, ZZ_python
from .quotientring import QuotientRingElement
from .ring import CommutativeRing
from .simpledomain import SimpleDomain


class IntegerModRing(CommutativeRing, SimpleDomain):
    """General class for quotient rings over integers."""

    is_Numerical = True

    def __new__(cls, order, dom):
        if isprime(order):
            return dom.finite_field(order)

        mod = dom.convert(order)

        key = cls, order, dom

        obj = super().__new__(cls)

        obj.domain = dom
        obj.mod = mod
        obj.order = order

        obj.rep = f'IntegerModRing({obj.order})'

        try:
            obj.dtype = _modular_integer_cache[key]
        except KeyError:
            obj.dtype = type('ModularInteger', (ModularInteger,),
                             {'mod': mod, 'domain': dom, '_parent': obj})
            _modular_integer_cache[key] = obj.dtype

        obj.zero = obj.dtype(0)
        obj.one = obj.dtype(1)

        return obj

    def __hash__(self):
        return hash((self.__class__.__name__, self.dtype, self.order, self.domain))

    def __eq__(self, other):
        return isinstance(other, self.__class__) and \
            self.order == other.order and self.domain == other.domain

    def __getnewargs_ex__(self):
        return (self.order,), {}

    @property
    def characteristic(self):
        return self.order

    def to_expr(self, element):
        return DiofantInteger(int(element))

    def from_expr(self, expr):
        if expr.is_Integer:
            return self.dtype(self.domain.dtype(int(expr)))
        if expr.is_Float and int(expr) == expr:
            return self.dtype(self.domain.dtype(int(expr)))
        raise CoercionFailedError(f'expected an integer, got {expr}')

    def _from_PythonFiniteField(self, a, K0=None):
        return self.dtype(self.domain.convert(a.rep, K0.domain))
    _from_GMPYFiniteField = _from_PythonFiniteField

    def _from_PythonIntegerRing(self, a, K0=None):
        return self.dtype(self.domain.convert(a, K0) % self.characteristic)
    _from_GMPYIntegerRing = _from_PythonIntegerRing

    def _from_PythonRationalField(self, a, K0=None):
        if a.denominator == 1:
            return self.convert(a.numerator)
    _from_GMPYRationalField = _from_PythonRationalField

    def _from_RealField(self, a, K0):
        p, q = K0.to_rational(a)

        if q == 1:
            return self.dtype(self.domain.dtype(p))

    def is_normal(self, a):
        return True


class FiniteField(Field, IntegerModRing):
    """General class for finite fields."""

    is_FiniteField = True

    def __new__(cls, order, dom, modulus=None):
        try:
            pp = factorint(order)
            if not order or len(pp) != 1:
                raise ValueError
            mod, deg = pp.popitem()
        except ValueError as exc:
            raise ValueError('order must be a prime power, '
                             f'got {order}') from exc

        if deg == 1:
            if modulus:
                deg = len(modulus) - 1
            else:
                modulus = [0, 1]

        order = mod**deg

        if modulus is None:
            random.seed(0)
            ring = ZZ_python.finite_field(mod).inject(Dummy('x'))
            modulus = ring._gf_random(deg, irreducible=True).all_coeffs()
        elif deg != len(modulus) - 1:
            raise ValueError('degree of a defining polynomial for the field'
                             ' does not match extension degree')

        modulus = tuple(map(dom.dtype, modulus))

        mod = dom.convert(mod)

        key = cls, order, dom, mod, modulus

        obj = super(IntegerModRing, cls).__new__(cls)  # pylint: disable=bad-super-call

        obj.domain = dom
        obj.mod = mod
        obj.order = order

        if order > mod:
            obj.rep = f'GF({obj.mod}, {list(map(ZZ_python, modulus))})'
        else:
            obj.rep = f'GF({obj.mod})'

        try:
            obj.dtype = _modular_integer_cache[key]
        except KeyError as exc:
            if deg == 1:
                obj.dtype = type('ModularInteger', (ModularInteger,),
                                 {'mod': mod, 'domain': dom, '_parent': obj})
            else:
                ff = dom.finite_field(mod).inject(Dummy('x'))
                mod = ff.from_list(modulus)
                if not mod.is_irreducible:
                    raise ValueError('defining polynomial must be '
                                     'irreducible') from exc
                obj.dtype = type('GaloisFieldElement', (GaloisFieldElement,),
                                 {'mod': mod, 'domain': ff, '_parent': obj})
            _modular_integer_cache[key] = obj.dtype

        obj.zero = obj.dtype(0)
        obj.one = obj.dtype(1)

        return obj

    @property
    def characteristic(self):
        return self.mod


_modular_integer_cache: dict[tuple, IntegerModRing] = {}


class PythonIntegerModRing(IntegerModRing):
    """Quotient ring based on Python's integers."""

    def __new__(cls, order):
        return super().__new__(cls, order, PythonIntegerRing())


class GMPYIntegerModRing(IntegerModRing):
    """Quotient ring based on GMPY's integers."""

    def __new__(cls, order):
        return super().__new__(cls, order, GMPYIntegerRing())


class PythonFiniteField(FiniteField):
    """Finite field based on Python's integers."""

    def __new__(cls, order, modulus=None):
        return super().__new__(cls, order, PythonIntegerRing(), modulus)


class GMPYFiniteField(FiniteField):
    """Finite field based on GMPY's integers."""

    def __new__(cls, order, modulus=None):
        return super().__new__(cls, order, GMPYIntegerRing(), modulus)


class ModularInteger(QuotientRingElement):
    """A class representing a modular integer."""

    @property
    def numerator(self):
        return self

    @property
    def denominator(self):
        return self.parent.one

    @property
    def is_primitive(self):
        """Test if this is a primitive element."""
        parent = self.parent
        return is_primitive_root(int(self), parent.order)


class GaloisFieldElement(ModularInteger):
    """A class representing a Galois field element."""

    def __init__(self, rep):
        """Initialize self."""
        if isinstance(rep, numbers.Integral):
            rep = list(reversed(integer_digits(rep % self.parent.order, self.parent.mod)))

        if isinstance(rep, (list, tuple)):
            rep = self.domain.from_list(rep)

        super().__init__(rep)

    def __int__(self):
        rep = self.rep.set_domain(self.parent.domain)
        return int(rep(self.parent.mod))

    @property
    def is_primitive(self):
        parent = self.parent
        p = parent.characteristic
        f = self.rep
        domain = self.domain
        x = domain.gens[0]
        n = f.degree()

        if not (f.is_irreducible and n):
            return False

        t = x**n

        for _ in range(n, p**n - 1):
            r = t % f
            if r == 1:
                return False
            t = r*x
        return True
