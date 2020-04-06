"""Sparse polynomial rings."""

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
from ..ntheory.modular import symmetric_residue
from ..utilities.magic import pollute
from .compatibility import IPolys
from .constructor import construct_domain
from .densebasic import dmp_from_dict, dmp_to_dict
from .heuristicgcd import heugcd
from .modulargcd import func_field_modgcd, modgcd
from .monomials import Monomial
from .orderings import lex
from .polyconfig import query
from .polyerrors import (CoercionFailed, DomainError, ExactQuotientFailed,
                         GeneratorsError, GeneratorsNeeded, HeuristicGCDFailed,
                         MultivariatePolynomialError, PolynomialDivisionFailed,
                         PolynomialError)
from .polyoptions import Domain as DomainOpt
from .polyoptions import Order as OrderOpt
from .polyoptions import build_options
from .polyutils import _dict_reorder, _parallel_dict_from_expr, expr_from_dict


__all__ = 'PolynomialRing', 'ring', 'sring', 'vring'


def ring(symbols, domain, order=lex):
    """Construct a polynomial ring returning ``(ring, x_1, ..., x_n)``.

    Parameters
    ==========

    symbols : str, Symbol/Expr or sequence of str, Symbol/Expr (non-empty)
    domain : :class:`~diofant.domains.domain.Domain` or coercible
    order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional, defaults to ``lex``

    Examples
    ========

    >>> R, x, y, z = ring('x y z', ZZ)
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

    >>> vring('x y z', ZZ)
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

    raise GeneratorsError('expected a string, Symbol or expression '
                          'or a non-empty sequence of strings, '
                          'Symbols or expressions')


_ring_cache = {}


class PolynomialRing(Ring, CompositeDomain, IPolys):
    """A class for representing multivariate polynomial rings."""

    is_PolynomialRing = True

    has_assoc_Ring = True

    def __new__(cls, domain, symbols, order=lex):
        symbols = tuple(_parse_symbols(symbols))
        ngens = len(symbols)
        domain = DomainOpt.preprocess(domain)
        order = OrderOpt.preprocess(order)

        key = (cls.__name__, symbols, ngens, domain, order)
        obj = _ring_cache.get(key)

        if obj is None:
            if domain.is_Composite and set(symbols) & set(domain.symbols):
                raise GeneratorsError("polynomial ring and it's ground domain share generators")

            obj = object.__new__(cls)
            obj._hash = hash(key)
            obj.dtype = type('PolyElement', (PolyElement,), {'ring': obj})
            obj.symbols = symbols
            obj.ngens = ngens
            obj.domain = domain
            obj.order = order

            obj.zero_monom = Monomial((0,)*ngens)
            obj.gens = obj._gens()
            obj._gens_set = set(obj.gens)

            obj._one = [(obj.zero_monom, domain.one)]

            obj.leading_expv = lambda f: Monomial(max(f, key=order))

            obj.rep = str(domain) + '[' + ','.join(map(str, symbols)) + ']'

            for symbol, generator in zip(obj.symbols, obj.gens):
                if isinstance(symbol, Symbol):
                    name = symbol.name

                    if not hasattr(obj, name):
                        setattr(obj, name, generator)

            _ring_cache[key] = obj

        return obj

    def __getnewargs_ex__(self):
        return (self.domain, self.symbols), {'order': self.order}

    @property
    def characteristic(self):
        return self.domain.characteristic

    def _gens(self):
        """Return a list of polynomial generators."""
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

    def clone(self, symbols=None, domain=None, order=None):
        return self.__class__(domain or self.domain, symbols or self.symbols, order or self.order)

    def monomial_basis(self, i):
        """Return the ith-basis element."""
        basis = [0]*self.ngens
        basis[i] = 1
        return Monomial(basis)

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

    def __call__(self, element):
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

    def from_dict(self, element):
        domain_new = self.domain_new
        poly = self.zero

        for monom, coeff in element.items():
            if isinstance(monom, int):
                monom = monom,
            coeff = domain_new(coeff)
            if coeff:
                poly[monom] = coeff

        return poly

    def from_terms(self, element):
        return self.from_dict(dict(element))

    def from_list(self, element):
        return self.from_dict(dmp_to_dict(element, self.ngens-1))

    def from_expr(self, expr):
        expr = sympify(expr)

        domain = self.domain
        mapping = dict(zip(self.symbols, self.gens))

        def _rebuild(expr):
            generator = mapping.get(expr)

            if generator is not None:
                return generator
            elif expr.is_Add:
                return functools.reduce(operator.add, map(_rebuild, expr.args))
            elif expr.is_Mul:
                return functools.reduce(operator.mul, map(_rebuild, expr.args))
            elif expr.is_Pow:
                c, a = expr.exp.as_coeff_Mul(rational=True)
                if c.is_Integer and c > 1:
                    return _rebuild(expr.base**a)**int(c)

            return domain.convert(expr)

        try:
            poly = _rebuild(expr)
        except CoercionFailed:
            raise ValueError('expected an expression convertible to a '
                             f'polynomial in {self}, got {expr}')
        else:
            return self(poly)

    def index(self, gen):
        """Compute index of ``gen`` in ``self.gens``."""
        if isinstance(gen, int):
            i = gen

            if 0 <= i and i < self.ngens:
                pass
            elif -self.ngens <= i and i <= -1:
                i = self.ngens + i
            else:
                raise ValueError(f'invalid generator index: {gen}')
        elif isinstance(gen, self.dtype):
            try:
                i = self.gens.index(gen)
            except ValueError:
                raise ValueError(f'invalid generator: {gen}')
        elif isinstance(gen, str):
            try:
                i = self.symbols.index(Symbol(gen))
            except ValueError:
                raise ValueError(f'invalid generator: {gen}')
        elif isinstance(gen, Expr):
            try:
                i = self.symbols.index(gen)
            except ValueError:
                raise ValueError(f'invalid generator: {gen}')
        else:
            raise ValueError(f'expected a polynomial generator, an integer, a string, an expression or None, got {gen}')

        return i

    def drop(self, *gens):
        """Remove specified generators from this ring."""
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]

        if not symbols:
            return self.domain
        else:
            return self.clone(symbols=symbols)

    def to_ground(self):
        if self.domain.is_Composite or self.domain.is_AlgebraicField:
            return self.clone(domain=self.domain.domain)
        else:
            raise ValueError(f'{self.domain} is not a composite or algebraic domain')

    @property
    def is_univariate(self):
        return len(self.gens) == 1

    @property
    def is_multivariate(self):
        return len(self.gens) > 1

    def eject(self, *gens):
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

    def to_expr(self, element):
        return element.as_expr()

    def _from_PythonIntegerRing(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_PythonFiniteField(self, a, K0):
        if self.domain == K0:
            return self(a)
    _from_GMPYFiniteField = _from_PythonFiniteField

    def _from_PythonRationalField(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_GMPYIntegerRing(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_GMPYRationalField(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_RealField(self, a, K0):
        return self(self.domain.convert(a, K0))

    def _from_ExpressionDomain(self, a, K0):
        if self.domain == K0:
            return self(a)

    def _from_AlgebraicField(self, a, K0):
        if self.domain == K0:
            return self(a)

    def _from_PolynomialRing(self, a, K0):
        try:
            return a.set_ring(self)
        except (CoercionFailed, GeneratorsError):
            return

    def _from_FractionField(self, a, K0):
        if self.domain == K0:
            return self.ground_new(a)

        (q,), r = a.numerator.div([a.denominator])

        if r.is_zero:
            return self.convert(q, K0.field.ring)

    @property
    def field(self):
        """Returns a field associated with ``self``."""
        return self.domain.frac_field(*self.symbols, order=self.order)

    def is_normal(self, a):
        return self.domain.is_normal(a.LC)

    def gcdex(self, a, b):
        """Extended GCD of ``a`` and ``b``."""
        return a.gcdex(b)

    def half_gcdex(self, a, b):
        """Half extended GCD of ``a`` and ``b``."""
        return a.half_gcdex(b)

    def gcd(self, a, b):
        """Returns GCD of ``a`` and ``b``."""
        return a.gcd(b)

    def lcm(self, a, b):
        """Returns LCM of ``a`` and ``b``."""
        return a.lcm(b)


class PolyElement(DomainElement, CantSympify, dict):
    """Element of multivariate distributed polynomial ring.

    See Also
    ========

    PolynomialRing

    """

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

    def __reduce__(self):
        return self.parent.__call__, (dict(self),)

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
        return self.__class__(self)

    def set_ring(self, new_ring):
        if self.ring == new_ring:
            return self
        elif self.ring == new_ring.domain:
            return new_ring.ground_new(self)
        elif set(new_ring.symbols).issuperset(self.ring.symbols):
            terms = zip(*_dict_reorder(self, self.ring.symbols, new_ring.symbols))
            return new_ring.from_terms(terms)
        else:
            raise CoercionFailed(f"Can't set element ring to {new_ring}")

    def set_domain(self, new_domain):
        if self.ring.domain == new_domain:
            return self
        else:
            new_ring = self.ring.clone(domain=new_domain)
            return self.set_ring(new_ring)

    def as_expr(self, *symbols):
        if not symbols:
            symbols = self.ring.symbols
        elif len(symbols) != self.ring.ngens:
            raise ValueError(f'not enough symbols, expected {self.ring.ngens} got {len(symbols)}')

        to_expr = self.ring.domain.to_expr
        return expr_from_dict({monom: to_expr(self[monom]) for monom in self}, *symbols)

    def clear_denoms(self):
        domain = self.ring.domain

        if not domain.is_Field:
            return domain.one, self

        if domain.has_assoc_Ring:
            ground_ring = domain.ring
        else:
            ground_ring = domain

        common = ground_ring.one
        lcm = ground_ring.lcm

        for coeff in self.values():
            common = lcm(common, coeff.denominator)

        return common, self.__class__({k: self[k]*common for k in self})

    def _strip_zero(self):
        """Eliminate monomials with zero coefficient."""
        for k, v in list(self.items()):
            if not v:
                del self[k]

    def __setitem__(self, key, item):
        if not isinstance(key, Monomial):
            key = Monomial(key)
        super().__setitem__(key, item)

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
        elif isinstance(other, self.ring.field.dtype):
            return other.__eq__(self)
        elif len(self) > 1:
            return False
        else:
            return self.get(self.ring.zero_monom) == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def almosteq(self, other, tolerance=None):
        """Approximate equality test for polynomials."""
        ring = self.ring

        if isinstance(other, ring.dtype):
            if set(self) != set(other):
                return False

            almosteq = ring.domain.almosteq

            for k in self:
                if not almosteq(self[k], other[k], tolerance):
                    return False
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

    def drop(self, gen):
        ring = self.ring
        i = ring.index(gen)

        if ring.is_univariate:
            if self.is_ground:
                return self.coeff(1)
            else:
                raise ValueError(f"can't drop {gen}")
        else:
            symbols = list(ring.symbols)
            del symbols[i]
            ring = ring.clone(symbols=symbols)

            poly = ring.zero

            for k, v in self.items():
                if k[i] == 0:
                    K = list(k)
                    del K[i]
                    poly[K] = v
                else:
                    raise ValueError(f"can't drop {gen}")

            return poly

    def eject(self, *gens):
        ring = self.ring

        if not gens:
            return self

        if ring.is_univariate:
            raise ValueError("can't drop only generator to ground")

        indexes = [ring.index(gen) for gen in gens]
        ring = ring.eject(*indexes)

        poly = ring.zero
        gens = ring.domain.gens[0:len(indexes)]

        for monom, coeff in self.items():
            mon = tuple(monom[i] for i in range(self.ring.ngens) if i not in indexes)
            gc = functools.reduce(operator.mul, [x**n for x, n in zip(gens, (monom[i] for i in indexes))])
            if mon in poly:
                poly[mon] += gc.mul_ground(coeff)
            else:
                poly[mon] = gc.mul_ground(coeff)

        return poly

    def inject(self, front=False):
        ring = self.ring
        domain = ring.domain

        if not (domain.is_Composite or domain.is_AlgebraicField):
            return self

        new_ring = ring.to_ground()
        new_ring = new_ring.inject(*domain.symbols, front=front)

        poly = new_ring.zero

        for monom, coeff in self.items():
            coeff = coeff.to_dict()
            for cmonom, ccoeff in coeff.items():
                if front:
                    cmonom += monom
                else:
                    cmonom = monom + cmonom

                poly[cmonom] = ccoeff

        return poly

    def to_dense(self):
        return dmp_from_dict(self, self.ring.ngens-1, self.ring.domain)

    def to_dict(self):
        return dict(self)

    def str(self, printer, precedence, exp_pattern, mul_symbol):
        if not self:
            return printer._print(self.ring.domain.zero)
        prec_add = precedence['Add']
        prec_atom = precedence['Atom']
        ring = self.ring
        symbols = ring.symbols
        ngens = ring.ngens
        zm = ring.zero_monom
        sexpvs = []
        for expv, coeff in self.terms():
            normal = ring.domain.is_normal(coeff)
            sign = ' + ' if normal else ' - '
            sexpvs.append(sign)
            if expv == zm:
                scoeff = printer._print(coeff)
                if scoeff.startswith('-'):
                    scoeff = scoeff[1:]
            else:
                if not normal:
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
                    sexpv.append(f'{symbol}')
            if scoeff:
                sexpv = [scoeff] + sexpv
            sexpvs.append(mul_symbol.join(sexpv))
        head = sexpvs.pop(0)
        if head == ' - ':
            sexpvs.insert(0, '-')
        return ''.join(sexpvs)

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
    def is_irreducible(self):
        _, factors = self.factor_list()
        if not factors:
            return True
        elif len(factors) > 1:
            return False
        else:
            return factors[0][1] == 1

    @property
    def is_homogeneous(self):
        if self.is_zero:
            return True

        tdeg = sum(self.LM)

        for monom in self:
            _tdeg = sum(monom)

            if _tdeg != tdeg:
                return False

        return True

    def __neg__(self):
        return self.__class__({monom: -self[monom] for monom in self})

    def __pos__(self):
        return self

    def __abs__(self):
        return self.__class__({monom: abs(self[monom]) for monom in self})

    def __add__(self, other):
        """Add two polynomials.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> (x + y)**2 + (x - y)**2
        2*x**2 + 2*y**2

        """
        ring = self.ring
        if not isinstance(other, ring.dtype):
            try:
                other = ring.convert(other)
            except CoercionFailed:
                return NotImplemented
        p = self.copy()
        if not other:
            return p
        get = p.get
        zero = ring.domain.zero
        for k, v in other.items():
            p[k] = get(k, zero) + v
        p._strip_zero()
        return p

    def __radd__(self, other):
        return self.__add__(other)

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
        ring = self.ring
        if not isinstance(other, ring.dtype):
            try:
                other = ring.convert(other)
            except CoercionFailed:
                return NotImplemented
        p = self.copy()
        if not other:
            return p
        get = p.get
        zero = ring.domain.zero
        for k, v in other.items():
            p[k] = get(k, zero) - v
        p._strip_zero()
        return p

    def __rsub__(self, other):
        """Substract self from other, with other convertible to the coefficient domain.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y
        >>> 4 - p
        -x - y + 4

        """
        return (-self).__add__(other)

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
        if not isinstance(other, ring.dtype):
            try:
                other = ring.convert(other)
            except CoercionFailed:
                return NotImplemented
        p = ring.zero
        if not self or not other:
            return p
        get = p.get
        zero = ring.domain.zero
        for exp1, v1 in self.items():
            for exp2, v2 in other.items():
                exp = exp1*exp2
                p[exp] = get(exp, zero) + v1*v2
        p._strip_zero()
        return p

    def __rmul__(self, other):
        """Multiply other to self with other in the coefficient domain of self.

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ)
        >>> p = x + y
        >>> 4 * p
        4*x + 4*y

        """
        return self.__mul__(other)

    def __pow__(self, n):
        """Raise polynomial to power `n`.

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
            raise ValueError('negative exponent')
        elif not n:
            return ring.one
        elif len(self) == 1:
            monom, coeff = list(self.items())[0]
            p = ring.zero
            p[monom**n] = coeff**n
            return p
        elif n == 1:
            return self.copy()
        elif n == 2:
            return self._square()
        elif n == 3:
            return self*self._square()
        elif len(self) <= 5:
            return self._pow_multinomial(n)
        else:
            return self._pow_generic(n)

    def _pow_generic(self, n):
        p = self.ring.one
        c = self

        while n:
            if n & 1:
                p *= c
                n -= 1

            c = c._square()
            n //= 2

        return p

    def _pow_multinomial(self, n):
        multinomials = multinomial_coefficients(len(self), n).items()
        zero_monom = self.ring.zero_monom
        terms = self.items()
        zero = self.ring.domain.zero
        poly = self.ring.zero

        for multinomial, multinomial_coeff in multinomials:
            product_monom = zero_monom
            product_coeff = multinomial_coeff

            for exp, (monom, coeff) in zip(multinomial, terms):
                if exp:
                    product_monom *= monom**exp
                    product_coeff *= coeff**exp

            monom = product_monom
            coeff = poly.get(monom, zero) + product_coeff

            if coeff:
                poly[monom] = coeff
            elif monom in poly:
                del poly[monom]

        return poly

    def _square(self):
        """Square of a polynomial.

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
                exp = k1*k2
                p[exp] = get(exp, zero) + pk*self[k2]
        p += p
        get = p.get
        for k, v in self.items():
            k2 = k**2
            p[k2] = get(k2, zero) + v**2
        p._strip_zero()
        return p

    def __divmod__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError('polynomial division')
        elif isinstance(other, ring.dtype):
            (q,), r = self.div([other])
            return q, r
        elif isinstance(other, PolyElement):
            if isinstance(ring.domain, PolynomialRing) and ring.domain.ring == other.ring:
                pass
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.quo_ground(other), self.trunc_ground(other)

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __truediv__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError('polynomial division')
        elif isinstance(other, ring.domain.dtype):
            return self.quo_ground(other)

        try:
            other = ring.domain_new(other)
        except CoercionFailed:
            return NotImplemented
        else:
            return self.quo_ground(other)

    def div(self, fv):
        """Division algorithm, see :cite:`Cox2015ideals`, p. 64.

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
        order = ring.order
        if any(not f for f in fv):
            raise ZeroDivisionError('polynomial division')
        if any(f.ring != ring for f in fv):
            raise ValueError('self and f must have the same ring')
        if not self:
            return [ring.zero], ring.zero
        s = len(fv)
        qv = [ring.zero for i in range(s)]
        p = self.copy()
        r = ring.zero
        expvs = [fx.leading_expv() for fx in fv]
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                lt = p.leading_term()
                term = lt.quo_term((expvs[i], fv[i][expvs[i]]))
                if term:
                    expv1, c = term.LT
                    qv[i] = qv[i]._iadd_monom((expv1, c))
                    p = p._iadd_poly_monom(fv[i], (expv1, -c))
                    divoccurred = 1
                    if p and order(p.LM) >= order(lt.LM):
                        raise PolynomialDivisionFailed(self, fv[i], self.ring)
                else:
                    i += 1
            if not divoccurred:
                expv = p.leading_expv()
                r = r._iadd_monom((expv, p[expv]))
                del p[expv]
        r._hash = None
        return qv, r

    def exquo(self, other):
        q, r = divmod(self, other)

        if not r:
            return q
        else:
            raise ExactQuotientFailed(self, other)

    def _iadd_monom(self, mc):
        """Add to self the monomial coeff*x0**i0*x1**i1*...
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
        """Add to self the product of (p)*(coeff*x0**i0*x1**i1*...)
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
            ka = k*m
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
        return max((monom[i] for monom in self), default=-oo)

    def degree_list(self):
        """
        A tuple containing leading degrees in all variables.

        Note that the degree of 0 is negative infinity (the Diofant object -oo)

        """
        return tuple(self.degree(i) for i in range(self.ring.ngens))

    def tail_degree(self, x=0):
        """
        The tail degree in ``x`` or the main variable.

        Note that the degree of 0 is negative infinity (the Diofant object -oo)

        """
        i = self.ring.index(x)
        return min((monom[i] for monom in self), default=-oo)

    def tail_degrees(self):
        """
        A tuple containing tail degrees in all variables.

        Note that the degree of 0 is negative infinity (the Diofant object -oo)

        """
        return tuple(self.tail_degree(i) for i in range(self.ring.ngens))

    def total_degree(self):
        """Returns the total degree."""
        return max(sum(m) for m in self.monoms())

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

        >>> _, x, y, z = ring('x y z', ZZ)
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
        elif isinstance(element, self.ring.dtype) and element.is_monomial:
            monom = element.monoms().pop()
            return self._get_coeff(monom)
        elif is_sequence(element) and all(isinstance(n, int) for n in element):
            return self._get_coeff(element)

        raise ValueError(f'expected a monomial, got {element}')

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

        >>> _, x, y = ring('x, y', ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.coeffs()
        [2, 1]
        >>> f.coeffs(grlex)
        [1, 2]

        """
        return [coeff for _, coeff in self.terms(order)]

    def all_coeffs(self):
        if self.ring.is_univariate:
            if self.is_zero:
                return [self.parent.domain.zero]
            else:
                return self.to_dense()
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def monoms(self, order=None):
        """Ordered list of polynomial monomials.

        Parameters
        ==========

        order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.monoms()
        [(2, 3), (1, 7)]
        >>> f.monoms(grlex)
        [(1, 7), (2, 3)]

        """
        return [monom for monom, _ in self.terms(order)]

    def terms(self, order=None):
        """Ordered list of polynomial terms.

        Parameters
        ==========

        order : :class:`~diofant.polys.polyoptions.Order` or coercible, optional

        Examples
        ========

        >>> _, x, y = ring('x, y', ZZ, lex)
        >>> f = x*y**7 + 2*x**2*y**3

        >>> f.terms()
        [((2, 3), 2), ((1, 7), 1)]
        >>> f.terms(grlex)
        [((1, 7), 1), ((2, 3), 2)]

        """
        return self._sorted(self.items(), order)

    def content(self):
        """Returns GCD of polynomial's coefficients."""
        ring = self.ring
        domain = ring.domain
        cont = domain.zero
        gcd = domain.gcd

        for coeff in self.values():
            cont = gcd(cont, coeff)

            if cont == domain.one:
                break

        if not ring.is_normal(self):
            cont = -cont

        return cont

    def primitive(self):
        """Returns content and a primitive polynomial."""
        cont = self.content()
        prim = self.copy()
        if not prim.is_zero:
            prim = prim.quo_ground(cont)
        return cont, prim

    def monic(self):
        """Divides all coefficients by the leading coefficient."""
        if not self:
            return self
        else:
            return self.exquo_ground(self.LC)

    def mul_ground(self, x):
        if not x:
            return self.ring.zero

        return self.__class__({monom: self[monom]*x for monom in self})

    def mul_monom(self, monom):
        terms = {f_monom*monom: self[f_monom] for f_monom in self}
        return self.__class__(terms)

    def mul_term(self, term):
        monom, coeff = term

        if not self or not coeff:
            return self.ring.zero
        elif monom == self.ring.zero_monom:
            return self.mul_ground(coeff)

        terms = {f_monom*monom: self[f_monom]*coeff for f_monom in self}
        return self.__class__(terms)

    def quo_ground(self, x):
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not self or x == domain.one:
            return self

        if domain.is_Field:
            quo = domain.quo
            terms = {monom: quo(self[monom], x) for monom in self}
        else:
            terms = {monom: self[monom]//x for monom in self}

        p = self.__class__(terms)
        p._strip_zero()
        return p

    def exquo_ground(self, x):
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not self or x == domain.one:
            return self

        terms = {monom: domain.exquo(self[monom], x) for monom in self}

        p = self.__class__(terms)
        p._strip_zero()
        return p

    def quo_term(self, term):
        monom, coeff = term

        if not coeff:
            raise ZeroDivisionError('polynomial division')
        elif not self:
            return self.ring.zero

        ring = self.ring
        domain = ring.domain
        p = ring.zero

        for tm, tc in self.items():
            if monom != self.ring.zero_monom:
                tm /= monom
            if any(_ < 0 for _ in tm):
                continue
            if domain.is_Field or not tc % coeff:
                p[tm] = domain.quo(tc, coeff)

        return p

    def trunc_ground(self, p):
        if self.ring.domain.is_IntegerRing:
            terms = {}

            for monom, coeff in self.items():
                coeff = coeff % p
                terms[monom] = symmetric_residue(coeff, p)
        else:
            terms = {monom: self[monom] % p for monom in self}

        poly = self.__class__(terms)
        poly._strip_zero()
        return poly

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
                h[N] = coeff

            H.append(h)

        return J, H

    def inflate(self, J):
        poly = self.ring.zero

        for I, coeff in self.items():
            N = [i*j for i, j in zip(I, J)]
            poly[N] = coeff

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

        J, (f, g) = self.deflate(other)
        h, cff, cfg = f._gcd(g)

        return h.inflate(J), cff.inflate(J), cfg.inflate(J)

    def _gcd_zero(self, other):
        one, zero = self.ring.one, self.ring.zero
        if self.ring.domain.is_Field:
            return other.monic(), zero, self.ring.ground_new(other.LC)
        else:
            if not self.ring.is_normal(other):
                return -other, zero, -one
            else:
                return other, zero, one

    def _gcd(self, other):
        ring = self.ring

        if ring.domain.is_RationalField:
            return self._gcd_QQ(other)
        elif ring.domain.is_IntegerRing:
            return self._gcd_ZZ(other)
        elif ring.domain.is_AlgebraicField:
            return self._gcd_AA(other)
        elif not ring.domain.is_Exact:
            try:
                exact = ring.domain.get_exact()
            except DomainError:
                return ring.one, self, other

            f, g = map(lambda x: x.set_domain(exact), (self, other))

            return tuple(map(lambda x: x.set_domain(ring.domain), f.cofactors(g)))
        elif ring.domain.is_Field:
            return self.ring.dmp_ff_prs_gcd(self, other)
        else:
            return self.ring.dmp_rr_prs_gcd(self, other)

    def _gcd_ZZ(self, other):
        if query('USE_HEU_GCD'):
            try:
                return heugcd(self, other)
            except HeuristicGCDFailed:  # pragma: no cover
                pass

        _gcd_zz_methods = {'modgcd': modgcd,
                           'prs': self.ring.dmp_rr_prs_gcd}

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
                           'prs': self.ring.dmp_ff_prs_gcd}

        method = _gcd_aa_methods[query('GCD_AA_METHOD')]
        return method(self, g)

    def terms_gcd(self):
        if self.is_zero:
            return (0,)*self.ring.ngens, self

        G = functools.reduce(Monomial.gcd, self)

        if all(g == 0 for g in G):
            return G, self

        f = self.ring.zero

        for monom, coeff in self.items():
            f[monom/G] = coeff

        return G, f

    def cancel(self, g, include=True):
        """
        Cancel common factors in a rational function ``f/g``.

        Examples
        ========

        >>> R, x, y = ring('x y', ZZ)

        >>> (2*x**2 - 2).cancel(x**2 - 2*x + 1)
        (2*x + 2, x - 1)

        """
        f = self
        ring = f.ring
        domain = ring.domain

        if not (domain.is_Field and domain.has_assoc_Ring):
            _, p, q = f.cofactors(g)
            cp, cq = domain.one, domain.one
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

        p_neg = not ring.is_normal(p)
        q_neg = not ring.is_normal(q)

        if p_neg and q_neg:
            p, q = -p, -q
        elif p_neg:
            cp, p = -cp, -p
        elif q_neg:
            cp, q = -cp, -q

        if not include:
            return cp, cq, p, q

        p = p.mul_ground(cp)
        q = q.mul_ground(cq)

        return p, q

    def diff(self, x=0, m=1):
        """Computes partial derivative in ``x``.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)
        >>> p = x + x**2*y**3
        >>> p.diff(x)
        2*x*y**3 + 1

        """
        ring = self.ring
        i = ring.index(x)
        x = ring.monomial_basis(i)
        x = x**m
        g = ring.zero if m else self.compose(ring.gens[i], ring.zero)
        for expv, coeff in self.items():
            if expv[i]:
                e = expv/x
                for j in range(expv[i], expv[i] - m, -1):
                    coeff *= j
                g[e] = coeff
        g._strip_zero()
        return g

    def integrate(self, x=0, m=1):
        """Computes indefinite integral in ``x``."""
        ring = self.ring
        i = ring.index(x)
        x = ring.monomial_basis(i)
        x = x**m
        g = ring.zero
        for expv, coeff in self.items():
            e = expv*x
            for j in range(expv[i] + m, expv[i], -1):
                coeff /= j
            g[e] = coeff
        g._strip_zero()
        return g

    def __call__(self, *values):
        if 0 < len(values) <= self.ring.ngens:
            return self.eval(list(zip(self.ring.gens, values)))
        else:
            raise ValueError(f'expected at least 1 and at most {self.ring.ngens} values, got {len(values)}')

    def eval(self, x=0, a=0):
        f = self

        if isinstance(x, list) and not a:
            (X, a), x = x[0], x[1:]
            f = f.eval(X, a)

            if not x:
                return f
            else:
                x = [(Y.drop(X), a) for (Y, a) in x]
                return f.eval(x)

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
                    coeff += poly[monom]

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

        if a is not None:
            replacements = [(x, a)]
        else:
            if isinstance(x, list):
                replacements = list(x)
            elif isinstance(x, dict):
                replacements = sorted(x.items(), key=lambda k: ring.index(k[0]))
            else:
                raise ValueError('expected a generator, value pair a sequence of such pairs')

        for monom, coeff in self.items():
            monom = list(monom)
            subpoly = ring.one

            for x, g in replacements:
                i, g = ring.index(x), ring(g)
                n, monom[i] = monom[i], 0
                if n:
                    subpoly *= g**n

            subpoly = subpoly.mul_term((monom, coeff))
            poly += subpoly

        return poly

    def discriminant(self):
        """Computes discriminant of a polynomial."""
        ring = self.ring

        d = self.degree()

        if d <= 0:
            return ring.zero.drop(0)
        else:
            s = (-1)**((d*(d - 1)) // 2)
            c = self.eject(*ring.gens[1:]).LC

            return self.resultant(self.diff()) // (c*s)

    def shift(self, a):
        if self.ring.is_univariate:
            return self.compose(0, self.ring.gens[0] + a)
        else:
            raise MultivariatePolynomialError('polynomial shift')

    def slice(self, m, n, x=0):
        ring = self.ring
        poly = ring.zero
        j = ring.index(x)

        for monom, coeff in self.items():
            if not n > monom[j] >= m:
                if ring.is_univariate:
                    continue
                else:
                    monom = monom[:j] + (0,) + monom[j + 1:]

            if monom in poly:
                poly[monom] += coeff
            else:
                poly[monom] = coeff

        return poly

    def prem(self, other):
        """Polynomial pseudo-remainder.

        Examples
        ========

        >>> R, x, y = ring('x y', ZZ)

        >>> (x**2 + x*y).prem(2*x + 2)
        -4*y + 4

        References
        ==========

        * :cite:`Knuth1985seminumerical`, p. 407.

        """
        ring = self.ring

        if not isinstance(other, ring.dtype):
            other = ring.convert(other)

        f, g = self, other

        if ring.is_multivariate:
            f, g = map(lambda _: _.eject(*ring.gens[1:]), (f, g))
            r = f.prem(g)
            return r.inject()

        ring = f.ring

        df = f.degree()
        dg = g.degree()

        if dg < 0:
            raise ZeroDivisionError('polynomial division')

        r, dr = f, df

        if df < dg:
            return r

        x = ring.gens[0]
        n = df - dg + 1
        lc_g = g.LC

        while True:
            lc_r = r.LC
            n -= 1

            r *= lc_g
            r -= g*x**(dr - dg)*lc_r
            dr = r.degree()

            if dr < dg:
                break

        r *= lc_g**n

        return r

    def subresultants(self, other):
        """
        Computes subresultant PRS of two polynomials in `K[X]`.

        Examples
        ========

        >>> R, x, y = ring('x y', ZZ)

        >>> f = 3*x**2*y - y**3 - 4
        >>> g = x**2 + x*y**3 - 9

        >>> a = 3*x*y**4 + y**3 - 27*y + 4
        >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

        >>> f.subresultants(g) == [f, g, a, b]
        True

        """
        return self.resultant(other, includePRS=True)[1]

    def half_gcdex(self, other):
        """
        Half extended Euclidean algorithm in `F[x]`.

        Returns ``(s, h)`` such that ``h = gcd(self, other)``
        and ``s*self = h (mod other)``.

        Examples
        ========

        >>> R, x = ring('x', QQ)

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> f.half_gcdex(g)
        (-1/5*x + 3/5, x + 1)

        """
        ring = self.ring
        if ring.is_multivariate:
            raise MultivariatePolynomialError('half extended Euclidean algorithm')

        domain = ring.domain

        if not domain.is_Field:
            raise DomainError(f"can't compute half extended GCD over {domain}")

        a, b = ring.one, ring.zero
        f, g = self, other

        while g:
            q, r = divmod(f, g)
            f, g = g, r
            a, b = b, a - q*b

        a = a.quo_ground(f.LC)
        f = f.monic()

        return a, f

    def gcdex(self, other):
        """
        Extended Euclidean algorithm in `F[x]`.

        Returns ``(s, t, h)`` such that ``h = gcd(self, other)`` and
        ``s*self + t*other = h``.

        Examples
        ========

        >>> R, x = ring('x', QQ)

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> f.gcdex(g)
        (-1/5*x + 3/5, 1/5*x**2 - 6/5*x + 2, x + 1)

        """
        s, h = self.half_gcdex(other)
        t = h - self*s
        t //= other
        return s, t, h

    # The following methods aren't ported (yet) to polynomial
    # representation independent algorithm implementations.

    def resultant(self, other, includePRS=False):
        return self.ring.dmp_resultant(self, other, includePRS=includePRS)

    def decompose(self):
        if self.ring.is_univariate:
            return self.ring.dup_decompose(self)
        else:
            raise MultivariatePolynomialError('polynomial decomposition')

    def sturm(self):
        if self.ring.is_univariate:
            return self.ring.dup_sturm(self)
        else:
            raise MultivariatePolynomialError('sturm sequence')

    @property
    def is_cyclotomic(self):
        if self.ring.is_univariate:
            return self.ring.dup_cyclotomic_p(self)
        else:
            raise AttributeError('cyclotomic polynomial')

    @property
    def is_squarefree(self):
        return self.ring.dmp_sqf_p(self)

    def sqf_norm(self):
        return self.ring.dmp_sqf_norm(self)

    def sqf_part(self):
        return self.ring.dmp_sqf_part(self)

    def sqf_list(self):
        return self.ring.dmp_sqf_list(self)

    def factor_list(self):
        return self.ring.dmp_factor_list(self)
