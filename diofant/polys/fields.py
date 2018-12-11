"""Sparse rational function fields. """

import functools
import operator

from ..core import Expr, Symbol, sympify
from ..core.sympify import CantSympify
from ..domains.compositedomain import CompositeDomain
from ..domains.domainelement import DomainElement
from ..domains.field import Field
from ..printing.defaults import DefaultPrinting
from ..utilities.magic import pollute
from .orderings import lex
from .polyerrors import CoercionFailed, GeneratorsError
from .rings import PolyElement, PolynomialRing


__all__ = ('FractionField', 'field', 'vfield')


def field(symbols, domain, order=lex):
    """Construct new rational function field returning (field, x1, ..., xn). """
    _field = FractionField(domain, symbols, order)
    return (_field,) + _field.gens


def vfield(symbols, domain, order=lex):
    """Construct new rational function field and inject generators into global namespace. """
    _field = FractionField(domain, symbols, order)
    pollute([sym.name for sym in _field.symbols], _field.gens)
    return _field


_field_cache = {}


class FractionField(Field, CompositeDomain):
    """A class for representing multivariate rational function fields. """

    is_FractionField = is_Frac = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def __new__(cls, domain, symbols, order=lex):
        ring = PolynomialRing(domain, symbols, order)
        symbols = ring.symbols
        ngens = ring.ngens
        domain = ring.domain
        order = ring.order

        _hash = hash((cls.__name__, symbols, ngens, domain, order))
        obj = _field_cache.get(_hash)

        if obj is None:
            obj = object.__new__(cls)
            obj._hash = _hash
            obj.dtype = type("FracElement", (FracElement,), {"field": obj})
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

            _field_cache[_hash] = obj

        return obj

    def _gens(self):
        """Return a list of polynomial generators. """
        return tuple(self.dtype(gen) for gen in self.ring.gens)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return self is other

    def __ne__(self, other):
        return self is not other

    def raw_new(self, numer, denom=None):
        return self.dtype(numer, denom)

    def new(self, numer, denom=None):
        if denom is None:
            denom = self.ring.one
        numer, denom = numer.cancel(denom)
        return self.raw_new(numer, denom)

    def domain_new(self, element):
        return self.domain.convert(element)

    def ground_new(self, element):
        try:
            return self.new(self.ring.ground_new(element))
        except CoercionFailed:
            domain = self.domain

            if not domain.is_Field and domain.has_assoc_Field:
                ring = self.ring
                ground_field = domain.field
                element = ground_field.convert(element)
                numer = ring.ground_new(element.numerator)
                denom = ring.ground_new(element.denominator)
                return self.raw_new(numer, denom)
            else:
                raise NotImplementedError

    def field_new(self, element):
        if isinstance(element, FracElement):
            if self == element.field:
                return element
            else:
                raise NotImplementedError("conversion")
        elif isinstance(element, PolyElement):
            denom, numer = element.clear_denoms()
            numer = numer.set_ring(self.ring)
            denom = self.ring.ground_new(denom)
            return self.raw_new(numer, denom)
        elif isinstance(element, tuple) and len(element) == 2:
            numer, denom = list(map(self.ring.ring_new, element))
            return self.new(numer, denom)
        elif isinstance(element, str):
            raise NotImplementedError("parsing")
        elif isinstance(element, Expr):
            return self.convert(element)
        else:
            return self.ground_new(element)

    __call__ = field_new

    def _rebuild_expr(self, expr, mapping):
        domain = self.domain

        def _rebuild(expr):
            generator = mapping.get(expr)

            if generator is not None:
                return generator
            elif expr.is_Add:
                return functools.reduce(operator.add, list(map(_rebuild, expr.args)))
            elif expr.is_Mul:
                return functools.reduce(operator.mul, list(map(_rebuild, expr.args)))
            elif expr.is_Pow:
                c, a = expr.exp.as_coeff_Mul(rational=True)
                if c.is_Integer and c != 1:
                    return _rebuild(expr.base**a)**int(c)

            if not domain.is_Field and domain.has_assoc_Field:
                return domain.field.convert(expr)
            else:
                return domain.convert(expr)

        return _rebuild(sympify(expr))

    def from_expr(self, expr):
        """Convert Diofant's expression to ``dtype``. """
        mapping = dict(zip(self.symbols, self.gens))

        try:
            frac = self._rebuild_expr(expr, mapping)
        except CoercionFailed:
            raise ValueError("expected an expression convertible to a rational function in %s, got %s" % (self, expr))
        else:
            return self.field_new(frac)

    def to_ring(self):
        return self.domain.poly_ring(*self.symbols, order=self.order)

    def to_expr(self, a):
        """Convert ``a`` to a Diofant object. """
        return a.as_expr()

    def _from_PythonIntegerRing(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_PythonRationalField(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_GMPYIntegerRing(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_GMPYRationalField(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_RealField(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_PolynomialRing(self, a, K0):
        try:
            return self.field_new(a)
        except (CoercionFailed, GeneratorsError):
            return

    def _from_FractionField(self, a, K0):
        try:
            return a.set_field(self)
        except (CoercionFailed, GeneratorsError):
            return

    @property
    def ring(self):
        """Returns a field associated with ``self``. """
        return self.to_ring()

    def is_positive(self, a):
        """Returns True if ``LC(a)`` is positive. """
        return self.domain.is_positive(a.numer.LC)

    def is_negative(self, a):
        """Returns True if ``LC(a)`` is negative. """
        return self.domain.is_negative(a.numer.LC)

    def is_nonpositive(self, a):
        """Returns True if ``LC(a)`` is non-positive. """
        return self.domain.is_nonpositive(a.numer.LC)

    def is_nonnegative(self, a):
        """Returns True if ``LC(a)`` is non-negative. """
        return self.domain.is_nonnegative(a.numer.LC)

    def factorial(self, a):
        """Returns factorial of ``a``. """
        return self.dtype(self.domain.factorial(a))


class FracElement(DomainElement, DefaultPrinting, CantSympify):
    """Element of multivariate distributed rational function field.

    See Also
    ========

    FractionField
    """

    def __init__(self, numer, denom=None):
        if denom is None:
            denom = self.field.ring.one
        elif not denom:
            raise ZeroDivisionError("zero denominator")

        self.numer = numer
        self.denom = denom

    def raw_new(self, numer, denom):
        return self.__class__(numer, denom)

    def new(self, numer, denom):
        return self.raw_new(*numer.cancel(denom))

    def to_poly(self):
        if self.denom != self.field.ring.one:
            raise ValueError("self.denom should be 1")
        return self.numer

    @property
    def numerator(self):
        return self.numer

    @property
    def denominator(self):
        return self.denom

    @property
    def parent(self):
        return self.field

    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.field, self.numer, self.denom))
        return _hash

    def copy(self):
        return self.raw_new(self.numer.copy(), self.denom.copy())

    def set_field(self, new_field):
        if self.field == new_field:
            return self
        else:
            new_ring = new_field.ring
            numer = self.numer.set_ring(new_ring)
            denom = self.denom.set_ring(new_ring)
            return new_field.new(numer, denom)

    def as_expr(self, *symbols):
        return self.numer.as_expr(*symbols)/self.denom.as_expr(*symbols)

    def __eq__(self, other):
        if isinstance(other, self.field.dtype):
            return self.numer == other.numer and self.denom == other.denom
        else:
            return self.numer == other and self.denom == self.field.ring.one

    def __bool__(self):
        return bool(self.numer)

    def sort_key(self):
        return self.denom.sort_key(), self.numer.sort_key()

    def _cmp(self, other, op):
        if isinstance(other, self.field.dtype):
            return op(self.sort_key(), other.sort_key())
        else:  # pragma: no cover
            return NotImplemented

    def __lt__(self, other):
        return self._cmp(other, operator.lt)

    def __le__(self, other):
        return self._cmp(other, operator.le)

    def __gt__(self, other):
        return self._cmp(other, operator.gt)

    def __ge__(self, other):
        return self._cmp(other, operator.ge)

    def __pos__(self):
        return self.raw_new(self.numer, self.denom)

    def __neg__(self):
        """Negate all coefficients in ``self``. """
        return self.raw_new(-self.numer, self.denom)

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
        """Add rational functions ``self`` and ``other``. """
        field = self.field

        if not other:
            return self
        elif not self:
            return other
        elif isinstance(other, field.dtype):
            if self.denom == other.denom:
                return self.new(self.numer + other.numer, self.denom)
            else:
                return self.new(self.numer*other.denom + self.denom*other.numer,
                                self.denom*other.denom)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numer + self.denom*other, self.denom)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__radd__(self)
                else:  # pragma: no cover
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__radd__(self)

        return self.__radd__(other)

    def __radd__(self, other):
        if isinstance(other, self.field.ring.dtype):
            return self.new(self.numer + self.denom*other, self.denom)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numer + self.denom*other_numer, self.denom)
        elif not op:  # pragma: no cover
            return NotImplemented
        else:
            return self.new(self.numer*other_denom + self.denom*other_numer,
                            self.denom*other_denom)

    def __sub__(self, other):
        """Subtract rational functions ``self`` and ``other``. """
        field = self.field

        if not other:
            return self
        elif not self:
            return -other
        elif isinstance(other, field.dtype):
            if self.denom == other.denom:
                return self.new(self.numer - other.numer, self.denom)
            else:
                return self.new(self.numer*other.denom - self.denom*other.numer,
                                self.denom*other.denom)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numer - self.denom*other, self.denom)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__rsub__(self)
                else:  # pragma: no cover
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__rsub__(self)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numer - self.denom*other_numer, self.denom)
        elif not op:  # pragma: no cover
            return NotImplemented
        else:
            return self.new(self.numer*other_denom - self.denom*other_numer,
                            self.denom*other_denom)

    def __rsub__(self, other):
        if isinstance(other, self.field.ring.dtype):
            return self.new(-self.numer + self.denom*other, self.denom)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(-self.numer + self.denom*other_numer, self.denom)
        elif not op:  # pragma: no cover
            return NotImplemented
        else:
            return self.new(-self.numer*other_denom + self.denom*other_numer,
                            self.denom*other_denom)

    def __mul__(self, other):
        """Multiply rational functions ``self`` and ``other``. """
        field = self.field

        if not self or not other:
            return field.zero
        elif isinstance(other, field.dtype):
            return self.new(self.numer*other.numer, self.denom*other.denom)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numer*other, self.denom)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__rmul__(self)
                else:  # pragma: no cover
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__rmul__(self)

        return self.__rmul__(other)

    def __rmul__(self, other):
        if isinstance(other, self.field.ring.dtype):
            return self.new(self.numer*other, self.denom)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numer*other_numer, self.denom)
        elif not op:  # pragma: no cover
            return NotImplemented
        else:
            return self.new(self.numer*other_numer, self.denom*other_denom)

    def __truediv__(self, other):
        """Computes quotient of fractions ``self`` and ``other``. """
        field = self.field

        if not other:
            raise ZeroDivisionError
        elif isinstance(other, field.dtype):
            return self.new(self.numer*other.denom, self.denom*other.numer)
        elif isinstance(other, field.ring.dtype):
            return self.new(self.numer, self.denom*other)
        else:
            if isinstance(other, FracElement):
                if isinstance(field.domain, FractionField) and field.domain.field == other.field:
                    pass
                elif isinstance(other.field.domain, FractionField) and other.field.domain.field == field:
                    return other.__rtruediv__(self)
                else:  # pragma: no cover
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__rtruediv__(self)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numer, self.denom*other_numer)
        elif not op:  # pragma: no cover
            return NotImplemented
        else:
            return self.new(self.numer*other_denom, self.denom*other_numer)

    def __rtruediv__(self, other):
        if not self:
            raise ZeroDivisionError
        elif isinstance(other, self.field.ring.dtype):
            return self.new(self.denom*other, self.numer)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.denom*other_numer, self.numer)
        elif not op:  # pragma: no cover
            return NotImplemented
        else:
            return self.new(self.denom*other_numer, self.numer*other_denom)

    def __pow__(self, n):
        """Raise ``self`` to a non-negative power ``n``. """
        if n >= 0:
            return self.raw_new(self.numer**n, self.denom**n)
        elif not self:
            raise ZeroDivisionError
        else:
            return self.raw_new(self.denom**-n, self.numer**-n)

    def diff(self, x):
        """Computes partial derivative in ``x``.

        Examples
        ========

        >>> _, x, y, z = field("x y z", ZZ)
        >>> ((x**2 + y)/(z + 1)).diff(x)
        2*x/(z + 1)
        """
        x = x.to_poly()
        return self.new(self.numer.diff(x)*self.denom -
                        self.numer*self.denom.diff(x), self.denom**2)

    def __call__(self, *values):
        if 0 < len(values) <= self.field.ngens:
            return self.evaluate(list(zip(self.field.gens, values)))
        else:
            raise ValueError("expected at least 1 and at most %s values, got %s" % (self.field.ngens, len(values)))

    def evaluate(self, x, a=None):
        if isinstance(x, list) and a is None:
            x = [(X.to_poly(), a) for X, a in x]
            numer, denom = self.numer.evaluate(x), self.denom.evaluate(x)
        else:
            x = x.to_poly()
            numer, denom = self.numer.evaluate(x, a), self.denom.evaluate(x, a)

        if self._extract_ground(denom) == (1, 1, None):
            return numer
        if isinstance(numer, PolyElement):
            field = numer.ring.to_field()
        else:
            field = self.field
        return field.new(field.ring(numer), field.ring(denom))

    def subs(self, x):
        if isinstance(x, (list, tuple)):
            x = [(X.to_poly(), a) for X, a in x]
            numer, denom = self.numer.subs(x), self.denom.subs(x)
        elif isinstance(x, (set, frozenset)):
            x = sorted((X.to_poly(), a) for X, a in x)
            numer, denom = self.numer.subs(x), self.denom.subs(x)
        elif isinstance(x, dict):
            x = sorted((k.to_poly(), x[k]) for k in x)
            numer, denom = self.numer.subs(x), self.denom.subs(x)
        else:
            raise ValueError("subs argument should be an iterable of pairs")

        return self.new(numer, denom)

    def compose(self, x, a=None):
        raise NotImplementedError
