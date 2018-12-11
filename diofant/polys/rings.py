"""Sparse polynomial rings. """

import functools
import math
import operator

from ..core import Expr, Symbol, oo
from ..core import symbols as _symbols
from ..core import sympify
from ..core.compatibility import is_sequence
from ..core.sympify import CantSympify
from ..domains.compositedomain import CompositeDomain
from ..domains.domainelement import DomainElement
from ..domains.ring import Ring
from ..ntheory import multinomial_coefficients
from ..printing.defaults import DefaultPrinting
from ..utilities.magic import pollute
from .compatibility import IPolys
from .constructor import construct_domain
from .densebasic import dmp_from_dict, dmp_to_dict
from .heuristicgcd import heugcd
from .modulargcd import func_field_modgcd, modgcd
from .monomials import (monomial_div, monomial_gcd, monomial_ldiv,
                        monomial_mul, monomial_pow)
from .orderings import lex
from .polyconfig import query
from .polyerrors import (CoercionFailed, ExactQuotientFailed, GeneratorsError,
                         GeneratorsNeeded, HeuristicGCDFailed,
                         MultivariatePolynomialError)
from .polyoptions import Domain as DomainOpt
from .polyoptions import Order as OrderOpt
from .polyoptions import build_options
from .polyutils import _dict_reorder, _parallel_dict_from_expr, expr_from_dict


__all__ = ('PolynomialRing', 'ring', 'sring', 'vring')


def ring(symbols, domain, order=lex):
    """Construct a polynomial ring returning ``(ring, x_1, ..., x_n)``.

    Parameters
    ==========

    symbols : str, Symbol/Expr or sequence of str, Symbol/Expr (non-empty)
    domain : :class:`~diofant.domains.domain.Domain` or coercible
    order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional, defaults to ``lex``

    Examples
    ========

    >>> R, x, y, z = ring("x y z", ZZ)
    >>> R
    ZZ[x,y,z]
    >>> x + y + z
    x + y + z
    """
    _ring = PolynomialRing(domain, symbols, order)
    return (_ring,) + _ring.gens


def vring(symbols, domain, order=lex):
    """Construct a polynomial ring and inject ``x_1, ..., x_n`` into the global namespace.

    Parameters
    ==========

    symbols : str, Symbol/Expr or sequence of str, Symbol/Expr (non-empty)
    domain : :class:`~diofant.domains.domain.Domain` or coercible
    order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional, defaults to ``lex``

    Examples
    ========

    >>> vring("x y z", ZZ)
    ZZ[x,y,z]
    >>> x + y + z
    x + y + z
    """
    _ring = PolynomialRing(domain, symbols, order)
    pollute([sym.name for sym in _ring.symbols], _ring.gens)
    return _ring


def sring(exprs, *symbols, **options):
    """Construct a ring deriving generators and domain from options and input expressions.

    Parameters
    ==========

    exprs : :class:`~diofant.core.expr.Expr` or sequence of :class:`~diofant.core.expr.Expr` (sympifiable)
    symbols : sequence of :class:`~diofant.core.symbol.Symbol`/:class:`~diofant.core.expr.Expr`
    options : keyword arguments understood by :class:`~diofant.polys.polyoptions.Options`

    Examples
    ========

    >>> R, f = sring(x + 2*y + 3*z)
    >>> R
    ZZ[x,y,z]
    >>> f
    x + 2*y + 3*z
    """
    single = False

    if not is_sequence(exprs):
        exprs, single = [exprs], True

    exprs = list(map(sympify, exprs))
    opt = build_options(symbols, options)

    # TODO: rewrite this so that it doesn't use expand() (see poly()).
    reps, opt = _parallel_dict_from_expr(exprs, opt)

    if opt.domain is None:
        # NOTE: this is inefficient because construct_domain() automatically
        # performs conversion to the target domain. It shouldn't do this.
        coeffs = sum((list(rep.values()) for rep in reps), [])
        opt.domain, _ = construct_domain(coeffs, opt=opt)

    _ring = PolynomialRing(opt.domain, opt.gens, opt.order)
    polys = list(map(_ring.from_dict, reps))

    if single:
        return _ring, polys[0]
    else:
        return _ring, polys


def _parse_symbols(symbols):
    if not symbols:
        raise GeneratorsNeeded("generators weren't specified")

    if isinstance(symbols, str):
        return _symbols(symbols, seq=True)
    elif isinstance(symbols, Expr):
        return symbols,
    elif is_sequence(symbols):
        if all(isinstance(s, str) for s in symbols):
            return _symbols(symbols)
        elif all(isinstance(s, Expr) for s in symbols):
            return symbols

    raise GeneratorsError("expected a string, Symbol or expression "
                          "or a non-empty sequence of strings, "
                          "Symbols or expressions")


_ring_cache = {}


class PolynomialRing(Ring, CompositeDomain, IPolys):
    """A class for representing multivariate polynomial rings. """

    is_PolynomialRing = is_Poly = True

    has_assoc_Ring = True
    has_assoc_Field = True

    def __new__(cls, domain, symbols, order=lex):
        symbols = tuple(_parse_symbols(symbols))
        ngens = len(symbols)
        domain = DomainOpt.preprocess(domain)
        order = OrderOpt.preprocess(order)

        _hash = hash((cls.__name__, symbols, ngens, domain, order))
        obj = _ring_cache.get(_hash)

        if obj is None:
            if domain.is_Composite and set(symbols) & set(domain.symbols):
                raise GeneratorsError("polynomial ring and it's ground domain share generators")

            obj = object.__new__(cls)
            obj._hash = _hash
            obj.dtype = type("PolyElement", (PolyElement,), {"ring": obj})
            obj.symbols = symbols
            obj.ngens = ngens
            obj.domain = domain
            obj.order = order

            obj.zero_monom = (0,)*ngens
            obj.gens = obj._gens()
            obj._gens_set = set(obj.gens)

            obj._one = [(obj.zero_monom, domain.one)]

            obj.leading_expv = lambda f: max(f, key=order)

            obj.rep = str(domain) + '[' + ','.join(map(str, symbols)) + ']'

            for symbol, generator in zip(obj.symbols, obj.gens):
                if isinstance(symbol, Symbol):
                    name = symbol.name

                    if not hasattr(obj, name):
                        setattr(obj, name, generator)

            _ring_cache[_hash] = obj

        return obj

    def _gens(self):
        """Return a list of polynomial generators. """
        one = self.domain.one
        _gens = []
        for i in range(self.ngens):
            expv = self.monomial_basis(i)
            poly = self.zero
            poly[expv] = one
            _gens.append(poly)
        return tuple(_gens)

    def __hash__(self):
        return self._hash

    def __eq__(self, other):
        return self is other

    def __ne__(self, other):
        return self is not other

    def clone(self, symbols=None, domain=None, order=None):
        return self.__class__(domain or self.domain, symbols or self.symbols, order or self.order)

    def monomial_basis(self, i):
        """Return the ith-basis element. """
        basis = [0]*self.ngens
        basis[i] = 1
        return tuple(basis)

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        return self.dtype(self._one)

    def domain_new(self, element, orig_domain=None):
        return self.domain.convert(element, orig_domain)

    def ground_new(self, coeff):
        return self.term_new(self.zero_monom, coeff)

    def term_new(self, monom, coeff):
        coeff = self.domain_new(coeff)
        poly = self.zero
        if coeff:
            poly[monom] = coeff
        return poly

    def ring_new(self, element):
        if isinstance(element, PolyElement):
            if self == element.ring:
                return element
            elif isinstance(self.domain, PolynomialRing) and self.domain.ring == element.ring:
                return self.ground_new(element)
            else:
                raise NotImplementedError
        elif isinstance(element, str):
            raise NotImplementedError
        elif isinstance(element, dict):
            return self.from_dict(element)
        elif isinstance(element, list):
            try:
                return self.from_terms(element)
            except ValueError:
                return self.from_list(element)
        elif isinstance(element, Expr):
            return self.convert(element)
        else:
            return self.ground_new(element)

    __call__ = new = ring_new

    def from_dict(self, element):
        domain_new = self.domain_new
        poly = self.zero

        for monom, coeff in element.items():
            coeff = domain_new(coeff)
            if coeff:
                poly[monom] = coeff

        return poly

    def from_terms(self, element):
        return self.from_dict(dict(element))

    def from_list(self, element):
        return self.from_dict(dmp_to_dict(element, self.ngens-1, self.domain))

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
                if c.is_Integer and c > 1:
                    return _rebuild(expr.base**a)**int(c)

            return domain.convert(expr)

        return _rebuild(sympify(expr))

    def from_expr(self, expr):
        """Convert Diofant's expression to ``dtype``. """
        mapping = dict(zip(self.symbols, self.gens))

        try:
            poly = self._rebuild_expr(expr, mapping)
        except CoercionFailed:
            raise ValueError("expected an expression convertible to a polynomial in %s, got %s" % (self, expr))
        else:
            return self.ring_new(poly)

    def index(self, gen):
        """Compute index of ``gen`` in ``self.gens``. """
        if isinstance(gen, int):
            i = gen

            if 0 <= i and i < self.ngens:
                pass
            elif -self.ngens <= i and i <= -1:
                i = self.ngens + i
            else:
                raise ValueError("invalid generator index: %s" % gen)
        elif isinstance(gen, self.dtype):
            try:
                i = self.gens.index(gen)
            except ValueError:
                raise ValueError("invalid generator: %s" % gen)
        else:
            raise ValueError("expected a polynomial generator, an integer, a string or None, got %s" % gen)

        return i

    def drop(self, *gens):
        """Remove specified generators from this ring. """
        indices = set(map(self.index, gens))
        symbols = [ s for i, s in enumerate(self.symbols) if i not in indices ]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    def to_ground(self):
        if self.domain.is_Composite or self.domain.is_AlgebraicField:
            return self.clone(domain=self.domain.domain)
        else:
            raise ValueError("%s is not a composite or algebraic "
                             "domain" % self.domain)

    def to_field(self):
        return self.domain.frac_field(*self.symbols, order=self.order)

    @property
    def is_univariate(self):
        return len(self.gens) == 1

    @property
    def is_multivariate(self):
        return len(self.gens) > 1

    def drop_to_ground(self, *gens):
        r"""
        Remove specified generators from the ring and inject them into
        its domain.
        """
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]
        gens = [gen for i, gen in enumerate(self.gens) if i not in indices]

        if not symbols:
            return self
        else:
            return self.clone(symbols=symbols, domain=self.drop(*gens))

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

    def _from_AlgebraicField(self, a, K0):
        if self.domain == K0:
            return self.new(a)

    def _from_PolynomialRing(self, a, K0):
        try:
            return a.set_ring(self)
        except (CoercionFailed, GeneratorsError):
            return

    def _from_FractionField(self, a, K0):
        denom = a.denominator

        if denom.is_ground:
            return self.convert(a.numerator/denom, K0.ring)

    @property
    def field(self):
        """Returns a field associated with ``self``. """
        return self.ring.to_field()

    def is_positive(self, a):
        """Returns True if ``LC(a)`` is positive. """
        return self.domain.is_positive(a.LC)

    def is_negative(self, a):
        """Returns True if ``LC(a)`` is negative. """
        return self.domain.is_negative(a.LC)

    def is_nonpositive(self, a):
        """Returns True if ``LC(a)`` is non-positive. """
        return self.domain.is_nonpositive(a.LC)

    def is_nonnegative(self, a):
        """Returns True if ``LC(a)`` is non-negative. """
        return self.domain.is_nonnegative(a.LC)

    def gcdex(self, a, b):
        """Extended GCD of ``a`` and ``b``. """
        return a.gcdex(b)

    def gcd(self, a, b):
        """Returns GCD of ``a`` and ``b``. """
        return a.gcd(b)

    def lcm(self, a, b):
        """Returns LCM of ``a`` and ``b``. """
        return a.lcm(b)

    def factorial(self, a):
        """Returns factorial of ``a``. """
        return self.new(self.domain.factorial(a))


class PolyElement(DomainElement, DefaultPrinting, CantSympify, dict):
    """Element of multivariate distributed polynomial ring.

    See Also
    ========

    PolynomialRing
    """

    def new(self, init):
        return self.__class__(init)

    @property
    def parent(self):
        return self.ring

    _hash = None

    def __hash__(self):
        # XXX: This computes a hash of a dictionary, but currently we don't
        # protect dictionary from being changed so any use site modifications
        # will make hashing go wrong. Use this feature with caution until we
        # figure out how to make a safe API without compromising speed of this
        # low-level class.
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.ring, frozenset(self.items())))
        return _hash

    def copy(self):
        """Return a copy of polynomial self.

        Polynomials are mutable; if one is interested in preserving
        a polynomial, and one plans to use inplace operations, one
        can copy the polynomial. This method makes a shallow copy.

        Examples
        ========

        >>> R, x, y = ring('x, y', ZZ)
        >>> p = (x + y)**2
        >>> p1 = p.copy()
        >>> p2 = p
        >>> p[R.zero_monom] = 3
        >>> p
        x**2 + 2*x*y + y**2 + 3
        >>> p1
        x**2 + 2*x*y + y**2
        >>> p2
        x**2 + 2*x*y + y**2 + 3
        """
        return self.new(self)

    def set_ring(self, new_ring):
        if self.ring == new_ring:
            return self
        elif self.ring.symbols != new_ring.symbols:
            terms = list(zip(*_dict_reorder(self, self.ring.symbols, new_ring.symbols)))
            return new_ring.from_terms(terms)
        else:
            return new_ring.from_dict(self)

    def as_expr(self, *symbols):
        if symbols and len(symbols) != self.ring.ngens:
            raise ValueError("not enough symbols, expected %s got %s" % (self.ring.ngens, len(symbols)))
        else:
            symbols = self.ring.symbols

        return expr_from_dict(self.as_expr_dict(), *symbols)

    def as_expr_dict(self):
        to_expr = self.ring.domain.to_expr
        return {monom: to_expr(coeff) for monom, coeff in self.items()}

    def clear_denoms(self):
        domain = self.ring.domain

        if not domain.is_Field or not domain.has_assoc_Ring:
            return domain.one, self

        ground_ring = domain.ring
        common = ground_ring.one
        lcm = ground_ring.lcm

        for coeff in self.values():
            common = lcm(common, coeff.denominator)

        poly = self.new([ (k, v*common) for k, v in self.items() ])
        return common, poly

    def strip_zero(self):
        """Eliminate monomials with zero coefficient. """
        for k, v in list(self.items()):
            if not v:
                del self[k]

    def __eq__(self, other):
        """Equality test for polynomials.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p1 = (x + y)**2 + (x - y)**2
        >>> p1 == 4*x*y
        False
        >>> p1 == 2*(x**2 + y**2)
        True
        """
        if not other:
            return not self
        elif isinstance(other, self.ring.dtype):
            return dict.__eq__(self, other)
        elif len(self) > 1:
            return False
        else:
            return self.get(self.ring.zero_monom) == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def almosteq(self, other, tolerance=None):
        """Approximate equality test for polynomials. """
        ring = self.ring

        if isinstance(other, ring.dtype):
            if set(self) != set(other):
                return False

            almosteq = ring.domain.almosteq

            for k in self:
                if not almosteq(self[k], other[k], tolerance):
                    return False
            else:
                return True
        elif len(self) > 1:
            return False
        else:
            try:
                other = ring.domain.convert(other)
            except CoercionFailed:
                return False
            else:
                return ring.domain.almosteq(self.coeff(1), other, tolerance)

    def sort_key(self):
        return len(self), self.terms()

    def _cmp(self, other, op):
        if isinstance(other, self.ring.dtype):
            return op(self.sort_key(), other.sort_key())
        else:
            return NotImplemented

    def __lt__(self, other):
        return self._cmp(other, operator.lt)

    def __le__(self, other):
        return self._cmp(other, operator.le)

    def __gt__(self, other):
        return self._cmp(other, operator.gt)

    def __ge__(self, other):
        return self._cmp(other, operator.ge)

    def _drop(self, gen):
        ring = self.ring
        i = ring.index(gen)

        if ring.is_univariate:
            return i, ring.domain
        else:
            symbols = list(ring.symbols)
            del symbols[i]
            return i, ring.clone(symbols=symbols)

    def drop(self, gen):
        i, ring = self._drop(gen)

        if self.ring.is_univariate:
            if self.is_ground:
                return self.coeff(1)
            else:
                raise ValueError("can't drop %s" % gen)
        else:
            poly = ring.zero

            for k, v in self.items():
                if k[i] == 0:
                    K = list(k)
                    del K[i]
                    poly[tuple(K)] = v
                else:
                    raise ValueError("can't drop %s" % gen)

            return poly

    def _drop_to_ground(self, gen):
        ring = self.ring
        i = ring.index(gen)

        symbols = list(ring.symbols)
        symbol = symbols[i]
        del symbols[i]
        return i, ring.clone(symbols=symbols, domain=ring.clone([symbol]))

    def drop_to_ground(self, gen):
        if self.ring.is_univariate:
            raise ValueError("can't drop only generator to ground")

        i, ring = self._drop_to_ground(gen)
        poly = ring.zero
        gen = ring.domain.gens[0]

        for monom, coeff in self.items():
            mon = monom[:i] + monom[i+1:]
            if mon not in poly:
                poly[mon] = (gen**monom[i]).mul_ground(coeff)
            else:
                poly[mon] += (gen**monom[i]).mul_ground(coeff)

        return poly

    def to_dense(self):
        return dmp_from_dict(self, self.ring.ngens-1, self.ring.domain)

    def to_dict(self):
        return dict(self)

    def str(self, printer, precedence, exp_pattern, mul_symbol):
        if not self:
            return printer._print(self.ring.domain.zero)
        prec_add = precedence["Add"]
        prec_atom = precedence["Atom"]
        ring = self.ring
        symbols = ring.symbols
        ngens = ring.ngens
        zm = ring.zero_monom
        sexpvs = []
        for expv, coeff in self.terms():
            positive = ring.domain.is_positive(coeff)
            sign = " + " if positive else " - "
            sexpvs.append(sign)
            if expv == zm:
                scoeff = printer._print(coeff)
                if scoeff.startswith("-"):
                    scoeff = scoeff[1:]
            else:
                if not positive:
                    coeff = -coeff
                if coeff != 1:
                    scoeff = printer.parenthesize(coeff, prec_add)
                else:
                    scoeff = ''
            sexpv = []
            for i in range(ngens):
                exp = expv[i]
                if not exp:
                    continue
                symbol = printer.parenthesize(symbols[i], prec_atom-1)
                if exp != 1:
                    sexpv.append(exp_pattern % (symbol, exp))
                else:
                    sexpv.append('%s' % symbol)
            if scoeff:
                sexpv = [scoeff] + sexpv
            sexpvs.append(mul_symbol.join(sexpv))
        head = sexpvs.pop(0)
        if head == " - ":
            sexpvs.insert(0, "-")
        return "".join(sexpvs)

    @property
    def is_generator(self):
        return self in self.ring._gens_set

    @property
    def is_ground(self):
        return not self or (len(self) == 1 and self.ring.zero_monom in self)

    @property
    def is_monomial(self):
        return not self or (len(self) == 1 and self.LC == 1)

    @property
    def is_term(self):
        return len(self) <= 1

    @property
    def is_negative(self):
        return self.ring.domain.is_negative(self.LC)

    @property
    def is_positive(self):
        return self.ring.domain.is_positive(self.LC)

    @property
    def is_nonnegative(self):
        return self.ring.domain.is_nonnegative(self.LC)

    @property
    def is_nonpositive(self):
        return self.ring.domain.is_nonpositive(self.LC)

    @property
    def is_zero(self):
        return not self

    @property
    def is_one(self):
        return self == self.ring.one

    @property
    def is_monic(self):
        return self.LC == self.ring.domain.one

    @property
    def is_primitive(self):
        return self.content() == self.ring.domain.one

    @property
    def is_linear(self):
        return all(sum(monom) <= 1 for monom in self)

    @property
    def is_quadratic(self):
        return all(sum(monom) <= 2 for monom in self)

    @property
    def is_squarefree(self):
        return self.ring.dmp_sqf_p(self)

    @property
    def is_irreducible(self):
        _, factors = self.factor_list()
        if not factors:
            return True
        elif len(factors) > 1:
            return False
        else:
            _, k = factors[0]
            return k == 1

    @property
    def is_cyclotomic(self):
        if self.ring.is_univariate:
            return self.ring.dup_cyclotomic_p(self)
        else:
            raise MultivariatePolynomialError("cyclotomic polynomial")

    def __neg__(self):
        return self.new([(monom, -coeff) for monom, coeff in self.items()])

    def __pos__(self):
        return self

    def __add__(self, other):
        """Add two polynomials.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> (x + y)**2 + (x - y)**2
        2*x**2 + 2*y**2
        """
        if not other:
            return self.copy()
        ring = self.ring
        if isinstance(other, ring.dtype):
            p = self.copy()
            get = p.get
            zero = ring.domain.zero
            for k, v in other.items():
                v = get(k, zero) + v
                if v:
                    p[k] = v
                else:
                    del p[k]
            return p
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif isinstance(other.ring.domain, PolynomialRing) and other.ring.domain.ring == ring:
                return other.__radd__(self)
            else:
                return NotImplemented

        try:
            cp2 = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            p = self.copy()
            zm = ring.zero_monom
            if zm not in self:
                p[zm] = cp2
            else:
                if other == -p[zm]:
                    del p[zm]
                else:
                    p[zm] += cp2
            return p

    def __radd__(self, n):
        p = self.copy()
        if not n:
            return p
        ring = self.ring
        try:
            n = ring.domain_new(n)
        except CoercionFailed:
            return NotImplemented
        else:
            zm = ring.zero_monom
            if zm not in self:
                p[zm] = n
            else:
                if n == -p[zm]:
                    del p[zm]
                else:
                    p[zm] += n
            return p

    def __sub__(self, other):
        """Subtract polynomial other from self.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p1 = x + y**2
        >>> p2 = x*y + y**2
        >>> p1 - p2
        -x*y + x
        """
        if not other:
            return self.copy()
        ring = self.ring
        if isinstance(other, ring.dtype):
            p = self.copy()
            get = p.get
            zero = ring.domain.zero
            for k, v in other.items():
                v = get(k, zero) - v
                if v:
                    p[k] = v
                else:
                    del p[k]
            return p
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif isinstance(other.ring.domain, PolynomialRing) and other.ring.domain.ring == ring:
                return other.__rsub__(self)
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            p = self.copy()
            zm = ring.zero_monom
            if zm not in self:
                p[zm] = -other
            else:
                if other == p[zm]:
                    del p[zm]
                else:
                    p[zm] -= other
            return p

    def __rsub__(self, n):
        """n - self with n convertible to the coefficient domain.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y
        >>> 4 - p
        -x - y + 4
        """
        ring = self.ring
        try:
            n = ring.domain_new(n)
        except CoercionFailed:
            return NotImplemented
        else:
            p = ring.zero
            for expv in self:
                p[expv] = -self[expv]
            p += n
            return p

    def __mul__(self, other):
        """Multiply two polynomials.

        Examples
        ========

        >>> _, x, y = ring('x, y', QQ)
        >>> p1 = x + y
        >>> p2 = x - y
        >>> p1*p2
        x**2 - y**2
        """
        ring = self.ring
        p = ring.zero
        if not self or not other:
            return p
        elif isinstance(other, ring.dtype):
            get = p.get
            zero = ring.domain.zero
            for exp1, v1 in self.items():
                for exp2, v2 in other.items():
                    exp = monomial_mul(exp1, exp2)
                    p[exp] = get(exp, zero) + v1*v2
            p.strip_zero()
            return p
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif isinstance(other.ring.domain, PolynomialRing) and other.ring.domain.ring == ring:
                return other.__rmul__(self)
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            for exp1, v1 in self.items():
                v = v1*other
                p[exp1] = v
            return p

    def __rmul__(self, other):
        """other * self with other in the coefficient domain of self.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y
        >>> 4 * p
        4*x + 4*y
        """
        p = self.ring.zero
        if not other:
            return p
        try:
            other = p.ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            for exp1, v1 in self.items():
                v = other*v1
                if v:
                    p[exp1] = v
            return p

    def __pow__(self, n):
        """raise polynomial to power `n`

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y**2
        >>> p**3
        x**3 + 3*x**2*y**2 + 3*x*y**4 + y**6
        """
        ring = self.ring
        n = int(n)

        if n < 0:
            raise ValueError("negative exponent")
        elif not n:
            if self:
                return ring.one
            else:
                raise ValueError("0**0")
        elif len(self) == 1:
            monom, coeff = list(self.items())[0]
            p = ring.zero
            p[monomial_pow(monom, n)] = coeff**n
            return p
        elif n == 1:
            return self.copy()
        elif n == 2:
            return self._square()
        elif n == 3:
            return self*self._square()
        elif len(self) <= 5:  # TODO: use an actuall density measure
            return self._pow_multinomial(n)
        else:
            return self._pow_generic(n)

    def _pow_generic(self, n):
        p = self.ring.one
        c = self

        while True:
            if n & 1:
                p = p*c
                n -= 1
                if not n:
                    break

            c = c._square()
            n = n // 2

        return p

    def _pow_multinomial(self, n):
        multinomials = list(multinomial_coefficients(len(self), n).items())
        zero_monom = self.ring.zero_monom
        terms = list(self.items())
        zero = self.ring.domain.zero
        poly = self.ring.zero

        for multinomial, multinomial_coeff in multinomials:
            product_monom = zero_monom
            product_coeff = multinomial_coeff

            for exp, (monom, coeff) in zip(multinomial, terms):
                if exp:
                    product_monom = monomial_mul(product_monom,
                                                 monomial_pow(monom, exp))
                    product_coeff *= coeff**exp

            monom = tuple(product_monom)
            coeff = product_coeff

            coeff = poly.get(monom, zero) + coeff

            if coeff:
                poly[monom] = coeff
            else:
                del poly[monom]

        return poly

    def _square(self):
        """square of a polynomial

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y**2
        >>> p._square()
        x**2 + 2*x*y**2 + y**4
        """
        ring = self.ring
        p = ring.zero
        get = p.get
        keys = list(self)
        zero = ring.domain.zero
        for i in range(len(keys)):
            k1 = keys[i]
            pk = self[k1]
            for j in range(i):
                k2 = keys[j]
                exp = monomial_mul(k1, k2)
                p[exp] = get(exp, zero) + pk*self[k2]
        p = p._imul_num(2)
        get = p.get
        for k, v in self.items():
            k2 = monomial_mul(k, k)
            p[k2] = get(k2, zero) + v**2
        p.strip_zero()
        return p

    def __divmod__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif isinstance(other, ring.dtype):
            return self.div(other)
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif isinstance(other.ring.domain, PolynomialRing) and other.ring.domain.ring == ring:
                return other.__rdivmod__(self)
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.quo_ground(other), self.rem_ground(other)

    def __rdivmod__(self, other):
        return NotImplemented

    def __mod__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif isinstance(other, ring.dtype):
            return divmod(self, other)[1]
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif isinstance(other.ring.domain, PolynomialRing) and other.ring.domain.ring == ring:
                return other.__rmod__(self)
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.rem_ground(other)

    def __rmod__(self, other):
        return NotImplemented

    def __floordiv__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif isinstance(other, ring.dtype):
            return divmod(self, other)[0]
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            elif isinstance(other.ring.domain, PolynomialRing) and other.ring.domain.ring == ring:
                return other.__rfloordiv__(self)
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.quo_ground(other)

    def __rfloordiv__(self, other):
        return NotImplemented

    def __truediv__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError("polynomial division")
        elif isinstance(other, ring.domain.dtype):
            return self.quo_ground(other)

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.quo_ground(other)

    def __rtruediv__(self, other):
        return NotImplemented

    def _term_div(self):
        zm = self.ring.zero_monom
        domain = self.ring.domain
        domain_quo = domain.quo

        if domain.is_Field:
            def term_div(a_lm_a_lc, b_lm_b_lc):
                a_lm, a_lc = a_lm_a_lc
                b_lm, b_lc = b_lm_b_lc
                if b_lm == zm:  # apparently this is a very common case
                    monom = a_lm
                else:
                    monom = monomial_div(a_lm, b_lm)
                if monom is not None:
                    return monom, domain_quo(a_lc, b_lc)
                else:
                    return
        else:
            def term_div(a_lm_a_lc, b_lm_b_lc):
                a_lm, a_lc = a_lm_a_lc
                b_lm, b_lc = b_lm_b_lc
                if b_lm == zm:  # apparently this is a very common case
                    monom = a_lm
                else:
                    monom = monomial_div(a_lm, b_lm)
                if not (monom is None or a_lc % b_lc):
                    return monom, domain_quo(a_lc, b_lc)
                else:
                    return

        return term_div

    def div(self, fv):
        """Division algorithm, see [CLO] p64.

        fv array of polynomials
           return qv, r such that
           self = sum(fv[i]*qv[i]) + r

        All polynomials are required not to be Laurent polynomials.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> f = x**3
        >>> f0 = x - y**2
        >>> f1 = x - y
        >>> qv, r = f.div((f0, f1))
        >>> qv[0]
        x**2 + x*y**2 + y**4
        >>> qv[1]
        0
        >>> r
        y**6
        """
        ring = self.ring
        ret_single = False
        if isinstance(fv, PolyElement):
            ret_single = True
            fv = [fv]
        if any(not f for f in fv):
            raise ZeroDivisionError("polynomial division")
        if not self:
            if ret_single:
                return ring.zero, ring.zero
            else:
                return [], ring.zero
        for f in fv:
            if f.ring != ring:
                raise ValueError('self and f must have the same ring')
        s = len(fv)
        qv = [ring.zero for i in range(s)]
        p = self.copy()
        r = ring.zero
        term_div = self._term_div()
        expvs = [fx.leading_expv() for fx in fv]
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                expv = p.leading_expv()
                term = term_div((expv, p[expv]), (expvs[i], fv[i][expvs[i]]))
                if term is not None:
                    expv1, c = term
                    qv[i] = qv[i]._iadd_monom((expv1, c))
                    p = p._iadd_poly_monom(fv[i], (expv1, -c))
                    divoccurred = 1
                else:
                    i += 1
            if not divoccurred:
                expv = p.leading_expv()
                r = r._iadd_monom((expv, p[expv]))
                del p[expv]
        if expv == ring.zero_monom:
            r += p
        if ret_single:
            return qv[0], r
        else:
            return qv, r

    def exquo(self, other):
        q, r = divmod(self, other)

        if not r:
            return q
        else:
            raise ExactQuotientFailed(self, other)

    def _iadd_monom(self, mc):
        """add to self the monomial coeff*x0**i0*x1**i1*...
        unless self is a generator -- then just return the sum of the two.

        mc is a tuple, (monom, coeff), where monomial is (i0, i1, ...)

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x**4 + 2*y
        >>> m = (1, 2)
        >>> p1 = p._iadd_monom((m, 5))
        >>> p1
        x**4 + 5*x*y**2 + 2*y
        >>> p1 is p
        True
        >>> p = x
        >>> p1 = p._iadd_monom((m, 5))
        >>> p1
        5*x*y**2 + x
        >>> p1 is p
        False
        """
        if self in self.ring._gens_set:
            cpself = self.copy()
        else:
            cpself = self
        expv, coeff = mc
        c = cpself.get(expv)
        if c is None:
            cpself[expv] = coeff
        else:
            c += coeff
            if c:
                cpself[expv] = c
            else:
                del cpself[expv]
        return cpself

    def _iadd_poly_monom(self, p2, mc):
        """add to self the product of (p)*(coeff*x0**i0*x1**i1*...)
        unless self is a generator -- then just return the sum of the two.

        mc is a tuple, (monom, coeff), where monomial is (i0, i1, ...)

        Examples
        ========

        >>> _, x, y, z = ring('x, y, z', ZZ)
        >>> p1 = x**4 + 2*y
        >>> p2 = y + z
        >>> m = (1, 2, 3)
        >>> p1 = p1._iadd_poly_monom(p2, (m, 3))
        >>> p1
        x**4 + 3*x*y**3*z**3 + 3*x*y**2*z**4 + 2*y
        """
        p1 = self
        if p1 in p1.ring._gens_set:
            p1 = p1.copy()
        m, c = mc
        get = p1.get
        zero = p1.ring.domain.zero
        for k, v in p2.items():
            ka = monomial_mul(k, m)
            coeff = get(ka, zero) + v*c
            if coeff:
                p1[ka] = coeff
            else:
                del p1[ka]
        return p1

    def degree(self, x=0):
        """
        The leading degree in ``x`` or the main variable.

        Note that the degree of 0 is negative infinity (the Diofant object -oo).
        """
        i = self.ring.index(x)

        if not self:
            return -oo
        else:
            return max(monom[i] for monom in self)

    def degrees(self):
        """
        A tuple containing leading degrees in all variables.

        Note that the degree of 0 is negative infinity (the Diofant object -oo)
        """
        if not self:
            return (-oo,)*self.ring.ngens
        else:
            return tuple(map(max, zip(*self)))

    def tail_degree(self, x=0):
        """
        The tail degree in ``x`` or the main variable.

        Note that the degree of 0 is negative infinity (the Diofant object -oo)
        """
        i = self.ring.index(x)

        if not self:
            return -oo
        else:
            return min(monom[i] for monom in self)

    def tail_degrees(self):
        """
        A tuple containing tail degrees in all variables.

        Note that the degree of 0 is negative infinity (the Diofant object -oo)
        """
        if not self:
            return (-oo,)*self.ring.ngens
        else:
            return tuple(map(min, zip(*self)))

    def leading_expv(self):
        """Leading monomial tuple according to the monomial ordering.

        Examples
        ========

        >>> _, x, y, z = ring('x, y, z', ZZ)
        >>> p = x**4 + x**3*y + x**2*z**2 + z**7
        >>> p.leading_expv()
        (4, 0, 0)
        """
        if self:
            return self.ring.leading_expv(self)
        else:
            return

    def _get_coeff(self, expv):
        return self.get(expv, self.ring.domain.zero)

    def coeff(self, element):
        """
        Returns the coefficient that stands next to the given monomial.

        Parameters
        ==========

        element : PolyElement (with ``is_monomial = True``) or 1

        Examples
        ========

        >>> _, x, y, z = ring("x y z", ZZ)
        >>> f = 3*x**2*y - x*y*z + 7*z**3 + 23

        >>> f.coeff(x**2*y)
        3
        >>> f.coeff(x*y)
        0
        >>> f.coeff(1)
        23
        """
        if element == 1:
            return self._get_coeff(self.ring.zero_monom)
        elif isinstance(element, self.ring.dtype):
            terms = list(element.items())
            if len(terms) == 1:
                monom, coeff = terms[0]
                if coeff == self.ring.domain.one:
                    return self._get_coeff(monom)

        raise ValueError("expected a monomial, got %s" % element)

    @property
    def LC(self):
        return self._get_coeff(self.leading_expv())

    @property
    def LM(self):
        expv = self.leading_expv()
        if expv is None:
            return self.ring.zero_monom
        else:
            return expv

    def leading_monom(self):
        """
        Leading monomial as a polynomial element.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> (3*x*y + y**2).leading_monom()
        x*y
        """
        p = self.ring.zero
        expv = self.leading_expv()
        if expv:
            p[expv] = self.ring.domain.one
        return p

    @property
    def LT(self):
        expv = self.leading_expv()
        if expv is None:
            return self.ring.zero_monom, self.ring.domain.zero
        else:
            return expv, self._get_coeff(expv)

    def leading_term(self):
        """Leading term as a polynomial element.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> (3*x*y + y**2).leading_term()
        3*x*y
        """
        p = self.ring.zero
        expv = self.leading_expv()
        if expv:
            p[expv] = self[expv]
        return p

    def _sorted(self, seq, order):
        if order is None:
            order = self.ring.order
        else:
            order = OrderOpt.preprocess(order)

        return sorted(seq, key=lambda monom: order(monom[0]), reverse=True)

    def coeffs(self, order=None):
        """Ordered list of polynomial coefficients.

        Parameters
        ==========

        order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional

        Examples
        ========

        >>> _, x, y = ring("x, y", ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.coeffs()
        [2, 1]
        >>> f.coeffs(grlex)
        [1, 2]
        """
        return [ coeff for _, coeff in self.terms(order) ]

    def monoms(self, order=None):
        """Ordered list of polynomial monomials.

        Parameters
        ==========

        order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional

        Examples
        ========

        >>> _, x, y = ring("x, y", ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.monoms()
        [(2, 3), (1, 7)]
        >>> f.monoms(grlex)
        [(1, 7), (2, 3)]
        """
        return [ monom for monom, _ in self.terms(order) ]

    def terms(self, order=None):
        """Ordered list of polynomial terms.

        Parameters
        ==========

        order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional

        Examples
        ========

        >>> _, x, y = ring("x, y", ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.terms()
        [((2, 3), 2), ((1, 7), 1)]
        >>> f.terms(grlex)
        [((1, 7), 1), ((2, 3), 2)]
        """
        return self._sorted(self.items(), order)

    def _imul_num(self, c):
        """multiply inplace the polynomial self by an element in the
        coefficient ring, provided self is not one of the generators;
        else multiply not inplace

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y**2
        >>> p1 = p._imul_num(3)
        >>> p1
        3*x + 3*y**2
        >>> p1 is p
        True
        >>> p = x
        >>> p1 = p._imul_num(3)
        >>> p1
        3*x
        >>> p1 is p
        False
        """
        if self in self.ring._gens_set:
            return self*c
        if not c:
            self.clear()
            return self
        for exp in self:
            self[exp] *= c
        return self

    def content(self):
        """Returns GCD of polynomial's coefficients. """
        domain = self.ring.domain
        cont = domain.zero
        gcd = domain.gcd

        for coeff in self.values():
            cont = gcd(cont, coeff)

        return cont

    def primitive(self):
        """Returns content and a primitive polynomial. """
        cont = self.content()
        return cont, self.quo_ground(cont)

    def monic(self):
        """Divides all coefficients by the leading coefficient. """
        if not self:
            return self
        else:
            return self.exquo_ground(self.LC)

    def mul_ground(self, x):
        if not x:
            return self.ring.zero

        terms = [(monom, coeff*x) for monom, coeff in self.items()]
        return self.new(terms)

    def mul_monom(self, monom):
        terms = [(monomial_mul(f_monom, monom), f_coeff) for f_monom, f_coeff in self.items()]
        return self.new(terms)

    def mul_term(self, term):
        monom, coeff = term

        if not self or not coeff:
            return self.ring.zero
        elif monom == self.ring.zero_monom:
            return self.mul_ground(coeff)

        terms = [(monomial_mul(f_monom, monom), f_coeff*coeff) for f_monom, f_coeff in self.items()]
        return self.new(terms)

    def quo_ground(self, x):
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not self or x == domain.one:
            return self

        if domain.is_Field:
            quo = domain.quo
            terms = [(monom, quo(coeff, x)) for monom, coeff in self.items()]
        else:
            terms = [(monom, coeff//x) for monom, coeff in self.items()]

        p = self.new(terms)
        p.strip_zero()
        return p

    def exquo_ground(self, x):
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not self or x == domain.one:
            return self

        terms = [(monom, domain.exquo(coeff, x)) for monom, coeff in self.items()]

        p = self.new(terms)
        p.strip_zero()
        return p

    def quo_term(self, term):
        monom, coeff = term

        if not coeff:
            raise ZeroDivisionError("polynomial division")
        elif not self:
            return self.ring.zero
        elif monom == self.ring.zero_monom:
            return self.quo_ground(coeff)

        term_div = self._term_div()

        terms = [term_div(t, term) for t in self.items()]
        return self.new([t for t in terms if t is not None])

    def trunc_ground(self, p):
        if self.ring.domain.is_IntegerRing:
            terms = []

            for monom, coeff in self.items():
                coeff = coeff % p

                if coeff > p // 2:
                    coeff = coeff - p

                terms.append((monom, coeff))
        else:
            terms = [(monom, coeff % p) for monom, coeff in self.items()]

        poly = self.new(terms)
        poly.strip_zero()
        return poly

    rem_ground = trunc_ground

    def extract_ground(self, g):
        f = self
        fc = f.content()
        gc = g.content()

        gcd = f.ring.domain.gcd(fc, gc)

        f = f.quo_ground(gcd)
        g = g.quo_ground(gcd)

        return gcd, f, g

    def _norm(self, norm_func):
        if not self:
            return self.ring.domain.zero
        else:
            return norm_func([abs(coeff) for coeff in self.values()])

    def max_norm(self):
        return self._norm(max)

    def l1_norm(self):
        return self._norm(sum)

    def deflate(self, *G):
        ring = self.ring
        polys = [self] + list(G)

        J = [0]*ring.ngens

        for p in polys:
            for monom in p:
                for i, m in enumerate(monom):
                    J[i] = math.gcd(J[i], m)

        for i, b in enumerate(J):
            if not b:
                J[i] = 1

        J = tuple(J)

        if all(b == 1 for b in J):
            return J, polys

        H = []

        for p in polys:
            h = ring.zero

            for I, coeff in p.items():
                N = [i//j for i, j in zip(I, J)]
                h[tuple(N)] = coeff

            H.append(h)

        return J, H

    def inflate(self, J):
        poly = self.ring.zero

        for I, coeff in self.items():
            N = [i*j for i, j in zip(I, J)]
            poly[tuple(N)] = coeff

        return poly

    def lcm(self, g):
        f = self
        domain = f.ring.domain

        if not domain.is_Field:
            fc, f = f.primitive()
            gc, g = g.primitive()
            c = domain.lcm(fc, gc)

        h = (f*g)//f.gcd(g)

        if not domain.is_Field:
            return h.mul_ground(c)
        else:
            return h.monic()

    def gcd(self, other):
        return self.cofactors(other)[0]

    def cofactors(self, other):
        if not self and not other:
            zero = self.ring.zero
            return zero, zero, zero
        elif not self:
            h, cff, cfg = self._gcd_zero(other)
            return h, cff, cfg
        elif not other:
            h, cfg, cff = other._gcd_zero(self)
            return h, cff, cfg
        elif len(self) == 1:
            h, cff, cfg = self._gcd_monom(other)
            return h, cff, cfg
        elif len(other) == 1:
            h, cfg, cff = other._gcd_monom(self)
            return h, cff, cfg

        J, (f, g) = self.deflate(other)
        h, cff, cfg = f._gcd(g)

        return h.inflate(J), cff.inflate(J), cfg.inflate(J)

    def _gcd_zero(self, other):
        one, zero = self.ring.one, self.ring.zero
        if other.is_nonnegative:
            return other, zero, one
        else:
            return -other, zero, -one

    def _gcd_monom(self, other):
        ring = self.ring
        ground_gcd = ring.domain.gcd
        ground_quo = ring.domain.quo
        mf, cf = list(self.items())[0]
        _mgcd, _cgcd = mf, cf
        for mg, cg in other.items():
            _mgcd = monomial_gcd(_mgcd, mg)
            _cgcd = ground_gcd(_cgcd, cg)
        h = self.new([(_mgcd, _cgcd)])
        cff = self.new([(monomial_ldiv(mf, _mgcd), ground_quo(cf, _cgcd))])
        cfg = self.new([(monomial_ldiv(mg, _mgcd), ground_quo(cg, _cgcd)) for mg, cg in other.items()])
        return h, cff, cfg

    def _gcd(self, other):
        ring = self.ring

        if ring.domain.is_RationalField:
            return self._gcd_QQ(other)
        elif ring.domain.is_IntegerRing:
            return self._gcd_ZZ(other)
        elif ring.domain.is_AlgebraicField:
            return self._gcd_AA(other)
        else:  # TODO: don't use dense representation (port PRS algorithms)
            return ring.dmp_inner_gcd(self, other)

    def _gcd_ZZ(self, other):
        if query('USE_HEU_GCD'):
            try:
                return heugcd(self, other)
            except HeuristicGCDFailed:  # pragma: no cover
                pass

        _gcd_zz_methods = {'modgcd': modgcd,
                           'prs': lambda f, g: self.ring.dmp_rr_prs_gcd(f, g)}

        method = _gcd_zz_methods[query('FALLBACK_GCD_ZZ_METHOD')]
        return method(self, other)

    def _gcd_QQ(self, g):
        f = self
        ring = f.ring
        new_ring = ring.clone(domain=ring.domain.ring)

        cf, f = f.clear_denoms()
        cg, g = g.clear_denoms()

        f = f.set_ring(new_ring)
        g = g.set_ring(new_ring)

        h, cff, cfg = f._gcd_ZZ(g)

        h = h.set_ring(ring)
        c, h = h.LC, h.monic()

        cff = cff.set_ring(ring).mul_ground(ring.domain.quo(c, cf))
        cfg = cfg.set_ring(ring).mul_ground(ring.domain.quo(c, cg))

        return h, cff, cfg

    def _gcd_AA(self, g):
        _gcd_aa_methods = {'modgcd': func_field_modgcd,
                           'prs': lambda f, g: self.ring.dmp_ff_prs_gcd(f, g)}

        method = _gcd_aa_methods[query('GCD_AA_METHOD')]
        return method(self, g)

    def cancel(self, g):
        """
        Cancel common factors in a rational function ``f/g``.

        Examples
        ========

        >>> R, x, y = ring("x y", ZZ)

        >>> (2*x**2 - 2).cancel(x**2 - 2*x + 1)
        (2*x + 2, x - 1)
        """
        f = self
        ring = f.ring

        if not f:
            return f, ring.one

        domain = ring.domain

        if not (domain.is_Field and domain.has_assoc_Ring):
            _, p, q = f.cofactors(g)

            if q.is_negative:
                p, q = -p, -q
        else:
            new_ring = ring.clone(domain=domain.ring)

            cq, f = f.clear_denoms()
            cp, g = g.clear_denoms()

            f = f.set_ring(new_ring)
            g = g.set_ring(new_ring)

            _, p, q = f.cofactors(g)
            _, cp, cq = new_ring.domain.cofactors(cp, cq)

            p = p.set_ring(ring)
            q = q.set_ring(ring)

            p_neg = p.is_negative
            q_neg = q.is_negative

            if p_neg and q_neg:
                p, q = -p, -q
            elif p_neg:
                cp, p = -cp, -p
            elif q_neg:
                cp, q = -cp, -q

            p = p.mul_ground(cp)
            q = q.mul_ground(cq)

        return p, q

    def diff(self, x):
        """Computes partial derivative in ``x``.

        Examples
        ========

        >>> _, x, y = ring("x y", ZZ)
        >>> p = x + x**2*y**3
        >>> p.diff(x)
        2*x*y**3 + 1
        """
        ring = self.ring
        i = ring.index(x)
        m = ring.monomial_basis(i)
        g = ring.zero
        for expv, coeff in self.items():
            if expv[i]:
                e = monomial_ldiv(expv, m)
                g[e] = coeff*expv[i]
        g.strip_zero()
        return g

    def __call__(self, *values):
        if 0 < len(values) <= self.ring.ngens:
            return self.evaluate(list(zip(self.ring.gens, values)))
        else:
            raise ValueError("expected at least 1 and at most %s values, got %s" % (self.ring.ngens, len(values)))

    def evaluate(self, x, a=None):
        f = self

        if isinstance(x, list) and a is None:
            (X, a), x = x[0], x[1:]
            f = f.evaluate(X, a)

            if not x:
                return f
            else:
                x = [ (Y.drop(X), a) for (Y, a) in x ]
                return f.evaluate(x)

        ring = f.ring
        i = ring.index(x)
        a = ring.domain.convert(a)

        if ring.is_univariate:
            result = ring.domain.zero

            for (n,), coeff in f.items():
                result += coeff*a**n

            return result
        else:
            poly = ring.drop(x).zero

            for monom, coeff in f.items():
                n, monom = monom[i], monom[:i] + monom[i+1:]
                coeff = coeff*a**n

                if monom in poly:
                    coeff = coeff + poly[monom]

                    if coeff:
                        poly[monom] = coeff
                    else:
                        del poly[monom]
                else:
                    if coeff:
                        poly[monom] = coeff

            return poly

    def subs(self, x):
        f = self
        a = None

        if isinstance(x, dict):
            return f.subs(sorted(x.items()))
        elif isinstance(x, (set, frozenset)):
            return f.subs(sorted(x))
        elif isinstance(x, (list, tuple)):
            if len(x) == 0:
                return f
            x = list(x)
            while len(x) > 1:
                f = f.subs([x[0]])
                del x[0]
            else:
                x, a = x[0]
        else:
            raise ValueError("subs argument should be an iterable of pairs")

        ring = f.ring
        i = ring.index(x)
        a = ring.domain.convert(a)

        if ring.is_univariate:
            result = ring.domain.zero

            for (n,), coeff in f.items():
                result += coeff*a**n

            return ring.ground_new(result)
        else:
            poly = ring.zero

            for monom, coeff in f.items():
                n, monom = monom[i], monom[:i] + (0,) + monom[i+1:]
                coeff = coeff*a**n

                if monom in poly:
                    coeff = coeff + poly[monom]

                    if coeff:
                        poly[monom] = coeff
                    else:
                        del poly[monom]
                else:
                    if coeff:
                        poly[monom] = coeff

            return poly

    def compose(self, x, a=None):
        """Computes the functional composition."""
        ring = self.ring
        poly = ring.zero
        gens_map = dict(zip(ring.gens, range(ring.ngens)))

        if a is not None:
            replacements = [(x, a)]
        else:
            if isinstance(x, list):
                replacements = list(x)
            elif isinstance(x, dict):
                replacements = sorted(x.items(), key=lambda k: gens_map[k[0]])
            else:
                raise ValueError("expected a generator, value pair a sequence of such pairs")

        for k, (x, g) in enumerate(replacements):
            replacements[k] = (gens_map[x], ring.ring_new(g))

        for monom, coeff in self.items():
            monom = list(monom)
            subpoly = ring.one

            for i, g in replacements:
                n, monom[i] = monom[i], 0
                if n:
                    subpoly *= g**n

            subpoly = subpoly.mul_term((tuple(monom), coeff))
            poly += subpoly

        return poly

    # TODO: following methods should point to polynomial
    # representation independent algorithm implementations.

    def pdiv(self, other):
        return self.ring.dmp_pdiv(self, other)

    def prem(self, other):
        return self.ring.dmp_prem(self, other)

    def pquo(self, other):
        return self.ring.dmp_quo(self, other)

    def pexquo(self, other):
        return self.ring.dmp_exquo(self, other)

    def half_gcdex(self, other):
        if self.ring.is_univariate:
            return self.ring.dup_half_gcdex(self, other)
        else:
            raise MultivariatePolynomialError("half extended Euclidean algorithm")

    def gcdex(self, other):
        if self.ring.is_univariate:
            return self.ring.dup_gcdex(self, other)
        else:
            raise MultivariatePolynomialError("extended Euclidean algorithm")

    def subresultants(self, other):
        return self.ring.dmp_subresultants(self, other)

    def resultant(self, other):
        return self.ring.dmp_resultant(self, other)

    def discriminant(self):
        return self.ring.dmp_discriminant(self)

    def decompose(self):
        if self.ring.is_univariate:
            return self.ring.dup_decompose(self)
        else:
            raise MultivariatePolynomialError("polynomial decomposition")

    def shift(self, a):
        if self.ring.is_univariate:
            return self.ring.dup_shift(self, a)
        else:
            raise MultivariatePolynomialError("polynomial shift")

    def sturm(self):
        if self.ring.is_univariate:
            return self.ring.dup_sturm(self)
        else:
            raise MultivariatePolynomialError("sturm sequence")

    def sqf_norm(self):
        return self.ring.dmp_sqf_norm(self)

    def sqf_part(self):
        return self.ring.dmp_sqf_part(self)

    def sqf_list(self):
        return self.ring.dmp_sqf_list(self)

    def factor_list(self):
        return self.ring.dmp_factor_list(self)
