"""Sparse rational function fields."""

from __future__ import annotations

import functools
import operator

from ..core import Expr, Symbol
from ..core.sympify import CantSympify, sympify
from ..domains.compositedomain import CompositeDomain
from ..domains.domainelement import DomainElement
from ..domains.field import Field
from .orderings import lex
from .polyerrors import CoercionFailed, GeneratorsError
from .rings import PolyElement, PolynomialRing


def field(symbols, domain, order=lex):
    """Construct new rational function field returning (field, x1, ..., xn)."""
    _field = FractionField(domain, symbols, order)
    return (_field,) + _field.gens


class FractionField(Field, CompositeDomain):
    """A class for representing multivariate rational function fields."""

    is_FractionField = True

    has_assoc_Ring = True

    def __new__(cls, domain, symbols, order=lex):
        ring = PolynomialRing(domain, symbols, order)
        symbols = ring.symbols
        ngens = ring.ngens
        domain = ring.domain
        order = ring.order

        key = cls.__name__, symbols, ngens, domain, order
        obj = _field_cache.get(key)

        if obj is None:
            obj = object.__new__(cls)
            obj._hash = hash(key)
            obj.dtype = type('FracElement', (FracElement,), {'field': obj})
            obj.symbols = symbols
            obj.ngens = ngens
            obj.domain = domain
            obj.order = order

            obj.zero = obj.dtype(ring.zero)
            obj.one = obj.dtype(ring.one)

            obj.gens = obj._gens()

            obj.rep = str(domain) + '(' + ','.join(map(str, symbols)) + ')'

            for symbol, generator in zip(obj.symbols, obj.gens):
                if isinstance(symbol, Symbol):
                    name = symbol.name

                    if not hasattr(obj, name):
                        setattr(obj, name, generator)

            _field_cache[key] = obj

        return obj

    def __getnewargs_ex__(self):
        return (self.domain, self.symbols), {'order': self.order}

    @property
    def characteristic(self):
        return self.domain.characteristic

    def _gens(self):
        """Return a list of polynomial generators."""
        return tuple(self.dtype(gen) for gen in self.ring.gens)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return self is other

    def clone(self, symbols=None, domain=None, order=None):
        return self.__class__(domain or self.domain, symbols or self.symbols, order or self.order)

    def __ne__(self, other):
        return self is not other

    def raw_new(self, numer, denom=None):
        return self.dtype(numer, denom)

    def domain_new(self, element):
        return self.domain.convert(element)

    def ground_new(self, element):
        try:
            return self(self.ring.ground_new(element))
        except CoercionFailed:
            domain = self.domain

            if not domain.is_Field and hasattr(domain, 'field'):
                ring = self.ring
                ground_field = domain.field
                element = ground_field.convert(element)
                numer = ring.ground_new(element.numerator)
                denom = ring.ground_new(element.denominator)
                return self.raw_new(numer, denom)
            else:
                raise NotImplementedError

    def __call__(self, element):
        if isinstance(element, FracElement):
            if self == element.field:
                return element
            else:
                raise NotImplementedError('conversion')
        elif isinstance(element, PolyElement):
            denom, numer = element.clear_denoms()
            numer = numer.set_ring(self.ring)
            denom = self.ring.ground_new(denom)
            return self.raw_new(numer, denom)
        elif isinstance(element, tuple) and len(element) == 2:
            numer, denom = list(map(self.ring.__call__, element))
            numer, denom = numer.cancel(denom)
            return self.raw_new(numer, denom)
        elif isinstance(element, str):
            raise NotImplementedError('parsing')
        elif isinstance(element, Expr):
            return self.convert(element)
        else:
            return self.ground_new(element)

    def from_expr(self, expr):
        expr = sympify(expr)
        domain = self.domain
        mapping = dict(zip(self.symbols, self.gens))

        def _rebuild(expr):
            if (generator := mapping.get(expr)) is not None:
                return generator
            elif expr.is_Add:
                return functools.reduce(operator.add, list(map(_rebuild, expr.args)))
            elif expr.is_Mul:
                return functools.reduce(operator.mul, list(map(_rebuild, expr.args)))
            elif expr.is_Pow:
                c, a = expr.exp.as_coeff_Mul(rational=True)
                if c.is_Integer and c != 1:
                    return _rebuild(expr.base**a)**int(c)

            if not domain.is_Field and hasattr(domain, 'field'):
                frac = domain.field.convert(expr)
            else:
                frac = domain.convert(expr)

            return self(frac)

        try:
            return _rebuild(expr)
        except CoercionFailed:
            raise ValueError('expected an expression convertible to a '
                             f'rational function in {self}, got {expr}')

    def to_ring(self):
        return self.domain.poly_ring(*self.symbols, order=self.order)

    def to_expr(self, element):
        ring = self.ring
        return ring.to_expr(element.numerator)/ring.to_expr(element.denominator)

    def _from_PythonIntegerRing(self, a, K0):
        return self(self.domain.convert(a, K0))
    _from_GMPYIntegerRing = _from_PythonIntegerRing
    _from_PythonRationalField = _from_PythonIntegerRing
    _from_GMPYRationalField = _from_PythonIntegerRing
    _from_RealField = _from_PythonIntegerRing
    _from_ComplexField = _from_PythonIntegerRing

    def _from_PolynomialRing(self, a, K0):
        try:
            return self(a)
        except (CoercionFailed, GeneratorsError):
            return

    def _from_FractionField(self, a, K0):
        try:
            return a.set_field(self)
        except (CoercionFailed, GeneratorsError):
            return

    @property
    def ring(self):
        return self.to_ring()

    def is_normal(self, a):
        return self.domain.is_normal(a.numerator.LC)


_field_cache: dict[tuple, FractionField] = {}


class FracElement(DomainElement, CantSympify):
    """Element of multivariate distributed rational function field.

    See Also
    ========

    FractionField

    """

    def __init__(self, numer, denom=None):
        if denom is None:
            denom = self.field.ring.one
        elif not denom:
            raise ZeroDivisionError('zero denominator')

        self._numerator = numer
        self._denominator = denom

    def __reduce__(self):
        return self.parent.__call__, ((self.numerator, self.denominator),)

    def raw_new(self, numer, denom):
        return self.__class__(numer, denom)

    def new(self, numer, denom):
        return self.raw_new(*numer.cancel(denom))

    def to_poly(self):
        if self.denominator != 1:
            raise ValueError('self.denominator should be 1')
        return self.numerator

    @property
    def numerator(self):
        return self._numerator

    @property
    def denominator(self):
        return self._denominator

    @property
    def parent(self):
        return self.field

    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.field, self.numerator, self.denominator))
        return _hash

    def copy(self):
        return self.raw_new(self.numerator.copy(), self.denominator.copy())

    def set_field(self, new_field):
        if self.field == new_field:
            return self
        else:
            new_ring = new_field.ring
            numer = self.numerator.set_ring(new_ring)
            denom = self.denominator.set_ring(new_ring)
            return new_field((numer, denom))

    def __eq__(self, other):
        if isinstance(other, self.field.dtype):
            return self.numerator == other.numerator and self.denominator == other.denominator
        else:
            return self.numerator == other and self.denominator == 1

    def __bool__(self):
        return bool(self.numerator)

    def __pos__(self):
        return self.raw_new(self.numerator, self.denominator)

    def __neg__(self):
        """Negate all coefficients in ``self``."""
        return self.raw_new(-self.numerator, self.denominator)

    def _extract_ground(self, element):
        domain = self.field.domain

        try:
            element = domain.convert(element)
        except CoercionFailed:
            ground_field = domain.field

            try:
                element = ground_field.convert(element)
            except CoercionFailed:
                return 0, None, None
            else:
                return -1, element.numerator, element.denominator
        else:
            return 1, element, None

    def __add__(self, other):
        """Add rational functions ``self`` and ``other``."""
        field = self.field

        if not other:
            return self
        elif not self:
            return other
        elif isinstance(other, field.dtype):
            if self.denominator == other.denominator:
                return self.new(self.numerator + other.numerator, self.denominator)
            else:
                return self.new(self.numerator*other.denominator + self.denominator*other.numerator,
                                self.denominator*other.denominator)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numerator + self.denominator*other, self.denominator)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__radd__(self)
                else:
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__radd__(self)

        return self.__radd__(other)

    def __radd__(self, other):
        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numerator + self.denominator*other_numer, self.denominator)
        elif not op:
            return NotImplemented
        else:
            return self.new(self.numerator*other_denom + self.denominator*other_numer,
                            self.denominator*other_denom)

    def __sub__(self, other):
        """Subtract rational functions ``self`` and ``other``."""
        field = self.field

        if not other:
            return self
        elif not self:
            return -other
        elif isinstance(other, field.dtype):
            if self.denominator == other.denominator:
                return self.new(self.numerator - other.numerator, self.denominator)
            else:
                return self.new(self.numerator*other.denominator - self.denominator*other.numerator,
                                self.denominator*other.denominator)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numerator - self.denominator*other, self.denominator)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__rsub__(self)
                else:
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__rsub__(self)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numerator - self.denominator*other_numer, self.denominator)
        elif not op:
            return NotImplemented
        else:
            return self.new(self.numerator*other_denom - self.denominator*other_numer,
                            self.denominator*other_denom)

    def __rsub__(self, other):
        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(-self.numerator + self.denominator*other_numer, self.denominator)
        elif not op:
            return NotImplemented
        else:
            return self.new(-self.numerator*other_denom + self.denominator*other_numer,
                            self.denominator*other_denom)

    def __mul__(self, other):
        """Multiply rational functions ``self`` and ``other``."""
        field = self.field

        if not self or not other:
            return field.zero
        elif isinstance(other, field.dtype):
            return self.new(self.numerator*other.numerator, self.denominator*other.denominator)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numerator*other, self.denominator)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__rmul__(self)
                else:
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__rmul__(self)

        return self.__rmul__(other)

    def __rmul__(self, other):
        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numerator*other_numer, self.denominator)
        elif not op:
            return NotImplemented
        else:
            return self.new(self.numerator*other_numer, self.denominator*other_denom)

    def __truediv__(self, other):
        """Computes quotient of fractions ``self`` and ``other``."""
        field = self.field

        if not other:
            raise ZeroDivisionError
        elif isinstance(other, field.dtype):
            return self.new(self.numerator*other.denominator, self.denominator*other.numerator)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numerator, self.denominator*other)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__rtruediv__(self)
                else:
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return NotImplemented

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numerator, self.denominator*other_numer)
        elif not op:
            return NotImplemented
        else:
            return self.new(self.numerator*other_denom, self.denominator*other_numer)

    def __rtruediv__(self, other):
        if not self:
            raise ZeroDivisionError
        elif isinstance(other, self.field.ring.dtype):
            return self.new(self.denominator*other, self.numerator)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.denominator*other_numer, self.numerator)
        elif not op:
            return NotImplemented
        else:
            return self.new(self.denominator*other_numer, self.numerator*other_denom)

    def __pow__(self, n):
        """Raise ``self`` to a non-negative power ``n``."""
        if n >= 0:
            return self.raw_new(self.numerator**n, self.denominator**n)
        elif not self:
            raise ZeroDivisionError
        else:
            return self.raw_new(self.denominator**-n, self.numerator**-n)

    def diff(self, x):
        """Computes partial derivative in ``x``.

        Examples
        ========

        >>> _, x, y, z = field('x y z', ZZ)
        >>> ((x**2 + y)/(z + 1)).diff(x)
        2*x/(z + 1)

        """
        x = x.to_poly()
        return self.new(self.numerator.diff(x)*self.denominator -
                        self.numerator*self.denominator.diff(x), self.denominator**2)

    def __call__(self, *values):
        if 0 < len(values) <= self.field.ngens:
            return self.eval(list(zip(self.field.gens, values)))
        else:
            raise ValueError(f'expected at least 1 and at most {self.field.ngens} values, got {len(values)}')

    def eval(self, x, a=None):
        if isinstance(x, list) and a is None:
            x = [(X.to_poly(), a) for X, a in x]
            numer, denom = self.numerator.eval(x), self.denominator.eval(x)
        else:
            x = x.to_poly()
            numer, denom = self.numerator.eval(x, a), self.denominator.eval(x, a)

        if self._extract_ground(denom) == (1, 1, None):
            return numer
        if isinstance(numer, PolyElement):
            field = numer.ring.field
        else:
            field = self.field
        return field((field.ring(numer), field.ring(denom)))

    def compose(self, x, a=None):
        """Computes the functional composition."""
        field = self.field

        if isinstance(x, list) and a is None:
            x = [(X.to_poly(), a) for X, a in x]
            numer = (self.numerator.compose([(X, a.numerator) for X, a in x]) *
                     self.denominator.compose([(X, a.denominator) for X, a in x]))
            denom = (self.numerator.compose([(X, a.denominator) for X, a in x]) *
                     self.denominator.compose([(X, a.numerator) for X, a in x]))
        else:
            x = x.to_poly()
            numer = (self.numerator.compose(x, a.numerator) *
                     self.denominator.compose(x, a.denominator))
            denom = (self.numerator.compose(x, a.denominator) *
                     self.denominator.compose(x, a.numerator))

        return field((field.ring(numer), field.ring(denom)))
