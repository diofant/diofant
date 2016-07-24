"""Sparse rational function fields. """

from functools import reduce
from operator import add, mul, lt, le, gt, ge

from diofant.core.expr import Expr
from diofant.core.symbol import Symbol
from diofant.core.sympify import CantSympify, sympify
from diofant.polys.rings import PolyElement
from diofant.polys.orderings import lex
from diofant.polys.polyerrors import CoercionFailed
from diofant.polys.domains.domainelement import DomainElement
from diofant.polys.domains.polynomialring import PolynomialRing
from diofant.polys.domains.fractionfield import FractionField
from diofant.printing.defaults import DefaultPrinting
from diofant.utilities import public
from diofant.utilities.magic import pollute


@public
def field(symbols, domain, order=lex):
    """Construct new rational function field returning (field, x1, ..., xn). """
    _field = FracField(symbols, domain, order)
    return (_field,) + _field.gens


@public
def xfield(symbols, domain, order=lex):
    """Construct new rational function field returning (field, (x1, ..., xn)). """
    _field = FracField(symbols, domain, order)
    return _field, _field.gens


@public
def vfield(symbols, domain, order=lex):
    """Construct new rational function field and inject generators into global namespace. """
    _field = FracField(symbols, domain, order)
    pollute([ sym.name for sym in _field.symbols ], _field.gens)
    return _field


@public
def sfield(exprs, *symbols, **options):
    """Construct a field deriving generators and domain from options and input expressions. """
    raise NotImplementedError

_field_cache = {}


class FracField(DefaultPrinting):
    """Multivariate distributed rational function field. """

    def __new__(cls, symbols, domain, order=lex):
        from diofant.polys.rings import PolyRing
        ring = PolyRing(symbols, domain, order)
        symbols = ring.symbols
        ngens = ring.ngens
        domain = ring.domain
        order = ring.order

        _hash = hash((cls.__name__, symbols, ngens, domain, order))
        obj = _field_cache.get(_hash)

        if obj is None:
            obj = object.__new__(cls)
            obj._hash = _hash
            obj.ring = ring
            obj.dtype = type("FracElement", (FracElement,), {"field": obj})
            obj.symbols = symbols
            obj.ngens = ngens
            obj.domain = domain
            obj.order = order

            obj.zero = obj.dtype(ring.zero)
            obj.one = obj.dtype(ring.one)

            obj.gens = obj._gens()

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

    def __getnewargs__(self):
        return self.symbols, self.domain, self.order

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

            if not domain.has_Field and domain.has_assoc_Field:
                ring = self.ring
                ground_field = domain.get_field()
                element = ground_field.convert(element)
                numer = ring.ground_new(ground_field.numer(element))
                denom = ring.ground_new(ground_field.denom(element))
                return self.raw_new(numer, denom)
            else:
                raise

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
            return self.from_expr(element)
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
                return reduce(add, list(map(_rebuild, expr.args)))
            elif expr.is_Mul:
                return reduce(mul, list(map(_rebuild, expr.args)))
            elif expr.is_Pow:
                c, a = expr.exp.as_coeff_Mul(rational=True)
                if c.is_Integer and c != 1:
                    return _rebuild(expr.base**a)**int(c)

            try:
                return domain.convert(expr)
            except CoercionFailed:
                if not domain.has_Field and domain.has_assoc_Field:
                    return domain.get_field().convert(expr)
                else:
                    raise

        return _rebuild(sympify(expr))

    def from_expr(self, expr):
        mapping = dict(zip(self.symbols, self.gens))

        try:
            frac = self._rebuild_expr(expr, mapping)
        except CoercionFailed:
            raise ValueError("expected an expression convertible to a rational function in %s, got %s" % (self, expr))
        else:
            return self.field_new(frac)

    def to_domain(self):
        return FractionField(self)

    def to_ring(self):
        from diofant.polys.rings import PolyRing
        return PolyRing(self.symbols, self.domain, self.order)


class FracElement(DomainElement, DefaultPrinting, CantSympify):
    """Element of multivariate distributed rational function field. """

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
        if self.denom != 1:
            raise ValueError("self.denom should be 1")
        return self.numer

    def parent(self):
        return self.field.to_domain()

    def __getnewargs__(self):
        return self.field, self.numer, self.denom

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

    def __ne__(self, other):
        return not self.__eq__(other)

    def __bool__(self):
        return bool(self.numer)

    def sort_key(self):
        return self.denom.sort_key(), self.numer.sort_key()

    def _cmp(self, other, op):
        if isinstance(other, self.field.dtype):
            return op(self.sort_key(), other.sort_key())
        else:
            return NotImplemented

    def __lt__(self, other):
        return self._cmp(other, lt)

    def __le__(self, other):
        return self._cmp(other, le)

    def __gt__(self, other):
        return self._cmp(other, gt)

    def __ge__(self, other):
        return self._cmp(other, ge)

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
            if not domain.has_Field and domain.has_assoc_Field:
                ground_field = domain.get_field()

                try:
                    element = ground_field.convert(element)
                except CoercionFailed:
                    pass
                else:
                    return -1, ground_field.numer(element), ground_field.denom(element)

            return 0, None, None
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
                else:
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
        elif not op:
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
                else:
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__rsub__(self)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numer - self.denom*other_numer, self.denom)
        elif not op:
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
        elif not op:
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
                else:
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
        elif not op:
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
                else:
                    return NotImplemented
            elif isinstance(other, PolyElement):
                if isinstance(field.domain, PolynomialRing) and field.domain.ring == other.ring:
                    pass
                else:
                    return other.__rtruediv__(self)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.numer, self.denom*other_numer)
        elif not op:
            return NotImplemented
        else:
            return self.new(self.numer*other_denom, self.denom*other_numer)

    __div__ = __truediv__

    def __rtruediv__(self, other):
        if not self:
            raise ZeroDivisionError
        elif isinstance(other, self.field.ring.dtype):
            return self.new(self.denom*other, self.numer)

        op, other_numer, other_denom = self._extract_ground(other)

        if op == 1:
            return self.new(self.denom*other_numer, self.numer)
        elif not op:
            return NotImplemented
        else:
            return self.new(self.denom*other_numer, self.numer*other_denom)

    __rdiv__ = __rtruediv__

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

        >>> from diofant.polys.domains import ZZ

        >>> _, x, y, z = field("x,y,z", ZZ)
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

        field = numer.ring.to_field()
        return field.new(numer, denom)

    def subs(self, x, a=None):
        if isinstance(x, list) and a is None:
            x = [(X.to_poly(), a) for X, a in x]
            numer, denom = self.numer.subs(x), self.denom.subs(x)
        else:
            x = x.to_poly()
            numer, denom = self.numer.subs(x, a), self.denom.subs(x, a)

        return self.new(numer, denom)

    def compose(self, x, a=None):
        raise NotImplementedError
