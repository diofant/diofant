"""Sparse polynomial rings."""

from __future__ import annotations

import functools
import math
import operator

from ..config import query
from ..core import Add, Expr, Symbol, cacheit
from ..core import symbols as _symbols
from ..core.sympify import CantSympify, sympify
from ..domains.compositedomain import CompositeDomain
from ..domains.domainelement import DomainElement
from ..domains.ring import CommutativeRing
from ..ntheory import multinomial_coefficients
from ..ntheory.modular import symmetric_residue
from ..utilities.iterables import is_sequence
from .euclidtools import _GCD
from .factortools import _Factor
from .monomials import Monomial
from .orderings import lex
from .polyerrors import (CoercionFailedError, DomainError,
                         ExactQuotientFailedError, GeneratorsError,
                         GeneratorsNeededError, PolynomialDivisionFailedError)
from .polyoptions import Domain as DomainOpt
from .polyoptions import Order as OrderOpt
from .specialpolys import _TestPolys
from .sqfreetools import _SQF


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


def _parse_symbols(symbols):
    if not symbols:
        raise GeneratorsNeededError("generators weren't specified")

    if isinstance(symbols, str):
        return _symbols(symbols, seq=True)
    if isinstance(symbols, Expr):
        return symbols,
    if is_sequence(symbols):
        if all(isinstance(s, str) for s in symbols):
            return _symbols(symbols)
        if all(isinstance(s, Expr) for s in symbols):
            return tuple(symbols)

    raise GeneratorsError('expected a string, Symbol or expression '
                          'or a non-empty sequence of strings, '
                          'Symbols or expressions')


class PolynomialRing(_GCD, CommutativeRing, CompositeDomain, _SQF, _Factor, _TestPolys):
    """Return a multivariate polynomial ring."""

    is_PolynomialRing = True

    has_assoc_Ring = True

    def __new__(cls, domain, symbols, order=lex):
        from .univar import UnivarPolyElement, UnivarPolynomialRing

        symbols = _parse_symbols(symbols)
        ngens = len(symbols)
        domain = DomainOpt.preprocess(domain)
        order = OrderOpt.preprocess(order)

        new_cls = PolynomialRing if ngens > 1 else UnivarPolynomialRing

        key = new_cls.__name__, symbols, ngens, domain, order
        obj = _ring_cache.get(key)

        if obj is None:
            if isinstance(domain, CompositeDomain) and set(symbols) & set(domain.symbols):
                raise GeneratorsError("polynomial ring and it's ground domain share generators")

            obj = object.__new__(new_cls)
            obj._hash = hash(key)

            if new_cls == UnivarPolynomialRing:
                dtype = UnivarPolyElement
            else:
                dtype = PolyElement
            obj.dtype = type(dtype.__name__, (dtype,), {'ring': obj})

            obj.symbols = symbols
            obj.ngens = ngens
            obj.domain = domain
            obj.order = order

            obj.zero_monom = Monomial((0,)*ngens)

            gens = []
            one = domain.one
            expv = [0]*ngens
            for i in range(ngens):
                expv[i] = 1
                poly = obj.zero
                poly[expv] = one
                gens.append(poly)
                expv[i] = 0
            obj.gens = tuple(gens)

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

    def leading_expv(self, f, order=None):
        order = self.order if order is None else OrderOpt.preprocess(order)
        return Monomial(max(f, key=order, default=self.zero_monom))

    @property
    def characteristic(self):
        return self.domain.characteristic

    def __hash__(self):
        return self._hash

    def clone(self, symbols=None, domain=None, order=None):
        return self.__class__(domain or self.domain, symbols or self.symbols, order or self.order)

    @property
    def zero(self):
        return self.dtype()

    @property
    def one(self):
        return self.ground_new(self.domain.one)

    def domain_new(self, element, orig_domain=None):
        return self.domain.convert(element, orig_domain)

    def ground_new(self, coeff):
        return self.term_new(self.zero_monom, coeff)

    def term_new(self, monom, coeff):
        poly = self.zero
        if coeff := self.domain.convert(coeff):
            poly[monom] = coeff
        return poly

    def __call__(self, element):
        if isinstance(element, PolyElement):
            if self == element.ring:
                return element
            if isinstance(self.domain, PolynomialRing) and self.domain.ring == element.ring:
                return self.ground_new(element)
            raise NotImplementedError
        if isinstance(element, str):
            raise NotImplementedError
        if isinstance(element, dict):
            return self.from_dict(element)
        if isinstance(element, list):
            try:
                return self.from_terms(element)
            except (TypeError, ValueError):
                return self.from_list(element)
        if isinstance(element, Expr):
            return self.convert(element)
        return self.ground_new(element)

    def from_dict(self, element):
        domain_new = self.domain.convert
        poly = self.zero

        for monom, coeff in element.items():
            if isinstance(monom, int):
                monom = monom,
            if coeff := domain_new(coeff):
                poly[monom] = coeff

        return poly

    def from_terms(self, element):
        return self.from_dict(dict(element))

    def from_list(self, element):
        if self.is_univariate:
            domain = self.domain
            if any(isinstance(c, list) for c in element):
                return self.from_dict({(i,): domain.from_list(c)
                                       for i, c in enumerate(element)})
            return self.from_dict({(i,): domain.convert(c)
                                   for i, c in enumerate(element)})
        new_ring = self.eject(*self.gens[1:])
        poly = new_ring.from_list(element)
        return poly.inject()

    def from_expr(self, expr):
        expr = sympify(expr)
        mapping = dict(zip(self.symbols, self.gens))

        def _rebuild(expr):
            if generator := mapping.get(expr):
                return generator
            if expr.is_Add:
                return sum(map(_rebuild, expr.args))
            if expr.is_Mul:
                return math.prod(map(_rebuild, expr.args))
            if expr.is_Pow:
                c, a = expr.exp.as_coeff_Mul(rational=True)
                if c.is_Integer and c > 1:
                    return _rebuild(expr.base**a)**int(c)
            return self.ground_new(self.domain.convert(expr))

        try:
            return _rebuild(expr)
        except CoercionFailedError as exc:
            raise ValueError('expected an expression convertible to a '
                             f'polynomial in {self}, got {expr}') from exc

    def index(self, gen):
        """Compute index of ``gen`` in ``self.gens``."""
        if isinstance(gen, int) and -self.ngens <= gen < self.ngens:
            return gen % self.ngens
        if isinstance(gen, self.dtype):
            return self.gens.index(gen)
        if isinstance(gen, str):
            gen = Symbol(gen)
        if isinstance(gen, Expr):
            return self.symbols.index(gen)
        raise ValueError('expected a polynomial generator, an integer, '
                         f'a string or an expression, got {gen}')

    def drop(self, *gens):
        """Remove specified generators from this ring."""
        indices = set(map(self.index, gens))
        symbols = [s for i, s in enumerate(self.symbols) if i not in indices]

        if not symbols:
            return self.domain
        return self.clone(symbols=symbols)

    def to_ground(self):
        domain = self.domain
        if isinstance(domain, CompositeDomain) or domain.is_AlgebraicField:
            return self.clone(domain=domain.domain)
        raise ValueError(f'{domain} is not a composite or algebraic domain')

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
        return self.clone(symbols=symbols, domain=self.drop(*gens))

    def to_expr(self, element):
        symbols = self.symbols
        to_expr = self.domain.to_expr
        return Add(*(to_expr(element[k])*k.as_expr(*symbols) for k in element))

    def _from_PythonFiniteField(self, a, K0):
        if self.domain == K0:
            return self(a)
    _from_GMPYFiniteField = _from_PythonFiniteField
    _from_ExpressionDomain = _from_PythonFiniteField

    def _from_AlgebraicField(self, a, K0):
        e = self.domain._from_AlgebraicField(a, K0)
        if e is not None:
            return self(a)

    def _from_PythonIntegerRing(self, a, K0):
        return self(self.domain.convert(a, K0))
    _from_GMPYIntegerRing = _from_PythonIntegerRing
    _from_PythonRationalField = _from_PythonIntegerRing
    _from_GMPYRationalField = _from_PythonIntegerRing
    _from_RealField = _from_PythonIntegerRing
    _from_ComplexField = _from_PythonIntegerRing

    def _from_PolynomialRing(self, a, K0):
        try:
            return a.set_ring(self)
        except (CoercionFailedError, GeneratorsError):
            return

    def _from_FractionField(self, a, K0):
        if self.domain == K0:
            return self.ground_new(a)

        (q,), r = a.numerator.div([a.denominator])

        if not r:
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


_ring_cache: dict[tuple, PolynomialRing] = {}


class PolyElement(DomainElement, CantSympify, dict):
    """Represent a polynomial in a multivariate polynomial ring.

    A polynomial is mutable, until its hash is computed, e.g. for
    using an instance as a dictionary key.

    See Also
    ========

    PolynomialRing

    """

    @property
    def parent(self):
        return self.ring

    _hash = None

    def __hash__(self):
        _hash = self._hash
        if _hash is None:
            self._hash = _hash = hash((self.ring, frozenset(self.items())))
        return _hash

    def __reduce__(self):
        return self.parent.__call__, (dict(self),)

    def copy(self):
        """Return a shallow copy of self."""
        return self.__class__(self)

    def set_ring(self, new_ring):
        ring = self.ring
        symbols = ring.symbols
        new_symbols = new_ring.symbols

        if ring == new_ring:
            return self
        if ring == new_ring.domain:
            return new_ring.ground_new(self)
        if set(new_symbols).issuperset(symbols):
            coeffs = self.values()

            new_monoms = [[] for _ in range(len(self))]

            for gen in new_symbols:
                try:
                    j = symbols.index(gen)

                    for M, new_M in zip(self, new_monoms):
                        new_M.append(M[j])
                except ValueError:
                    for new_M in new_monoms:
                        new_M.append(0)

            terms = zip(map(Monomial, new_monoms), coeffs)

            return new_ring.from_terms(terms)
        raise CoercionFailedError(f"Can't set element ring to {new_ring}")

    def set_domain(self, new_domain):
        if self.ring.domain == new_domain:
            return self
        new_ring = self.ring.clone(domain=new_domain)
        return self.set_ring(new_ring)

    def clear_denoms(self, convert=False):
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

        f = self*domain.convert(common)

        if convert:
            f = f.set_domain(ground_ring)

        return common, f

    def _strip_zero(self):
        """Eliminate monomials with zero coefficient."""
        for k, v in list(self.items()):
            if not v:
                del self[k]

    def __setitem__(self, monom, coeff, /):
        """Set the coefficient for the given monomial.

        Parameters
        ==========

        monom : Monomial or PolyElement (with ``is_monomial = True``) or 1
        coeff : DomainElement

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)

        >>> p = (x + y)**2
        >>> p1 = p.copy()
        >>> p2 = p

        >>> p[x*y] = 0
        >>> p1
        x**2 + 2*x*y + y**2
        >>> p2
        x**2 + y**2

        >>> _ = hash(p)
        >>> p[x*y] = 1
        Traceback (most recent call last):
        ...
        RuntimeError: Polynomial x**2 + y**2 can't be modified

        """
        if self._hash is not None:
            raise RuntimeError(f"polynomial {self} can't be modified")

        ring = self.ring

        if isinstance(monom, Monomial):
            pass
        elif isinstance(monom, (tuple, list)):
            monom = Monomial(monom)
        elif isinstance(monom, ring.dtype) and monom.is_monomial:
            monom, = monom
        elif monom == 1:
            monom = ring.zero_monom
        else:
            raise TypeError(f'monomial expected, got {monom}')

        if coeff:
            super().__setitem__(monom, coeff)
        elif monom in self:
            del self[monom]

    def __eq__(self, other):
        """Equality test for polynomials.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)
        >>> p1 = (x + y)**2 + (x - y)**2
        >>> p1 == 4*x*y
        False
        >>> p1 == 2*(x**2 + y**2)
        True

        """
        ring = self.ring
        if not other:
            return not self
        if isinstance(other, ring.dtype):
            return dict.__eq__(self, other)
        if isinstance(other, ring.field.dtype):
            return other.__eq__(self)
        if len(self) > 1:
            return False
        return self[ring.zero_monom] == other

    def __ne__(self, other):
        return not self.__eq__(other)

    def drop(self, *gens):
        ring = self.ring
        indexes = sorted((ring.index(gen) for gen in gens), reverse=True)
        new_ring = ring.drop(*indexes)

        if new_ring == ring.domain:
            if self.is_ground:
                return self[1]
            raise ValueError(f"can't drop {gens}")

        poly = new_ring.zero

        for k, v in self.items():
            K = list(k)
            for i in indexes:
                if k[i] == 0:
                    del K[i]
                else:
                    raise ValueError(f"can't drop {gens}")
            poly[K] = v

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
            gc = math.prod((x**n for x, n in zip(gens, (monom[i] for i in indexes))))
            if mon in poly:
                poly[mon] += gc*coeff
            else:
                poly[mon] = gc*coeff

        return poly

    def inject(self, front=False):
        ring = self.ring
        domain = ring.domain

        if not (isinstance(domain, CompositeDomain) or domain.is_AlgebraicField):
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

    def to_dict(self):
        return dict(self)

    def all_coeffs(self):
        ring = self.ring
        if not (ground_gens := ring.gens[1:]):
            if self:
                return [self[(i,)] for i in range(self.degree() + 1)]
            return [self[(0,)]]
        poly = self.eject(*ground_gens)
        if poly:
            return [poly[(i,)].all_coeffs() for i in range(poly.degree() + 1)]
        return [poly[(0,)].all_coeffs()]

    def _str(self, printer, precedence, exp_pattern, mul_symbol):
        if not self:
            return printer._print(self.ring.domain.zero)
        prec_add = precedence['Add']
        prec_atom = precedence['Atom']
        ring = self.ring
        symbols = ring.symbols
        ngens = ring.ngens
        zm = ring.zero_monom
        order = ring.order
        sexpvs = []
        for expv, coeff in sorted(self.items(),
                                  key=lambda m: order(m[0]),
                                  reverse=True):
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
        return self.is_monomial and self.total_degree() == 1

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
    def is_linear(self):
        return all(sum(monom) <= 1 for monom in self)

    @property
    def is_quadratic(self):
        return all(sum(monom) <= 2 for monom in self)

    @property
    def is_irreducible(self):
        ring = self.ring
        domain = ring.domain

        if ring.is_univariate:
            if domain.is_FiniteField:
                method = query('GF_IRRED_METHOD')
                _irred_methods = {'ben-or': ring._gf_irreducible_p_ben_or,
                                  'rabin': ring._gf_irreducible_p_rabin}
                return _irred_methods[method](self)
            if domain.is_IntegerRing:
                res = ring._zz_irreducible_p(self)
                if res is not None:
                    return res

        _, factors = self.factor_list()

        if not factors:
            return True
        if len(factors) > 1:
            return False
        return factors[0][1] == 1

    @property
    def is_homogeneous(self):
        if not self:
            return True

        lm = self.LM
        tdeg = sum(lm)
        return all(sum(monom) == tdeg for monom in self if monom != lm)

    def __neg__(self):
        return self.__class__({monom: -self[monom] for monom in self})

    def __pos__(self):
        return self

    def __abs__(self):
        return self.__class__({monom: abs(self[monom]) for monom in self})

    def __add__(self, other):
        """Add two polynomials."""
        ring = self.ring
        try:
            other = ring.convert(other)
        except CoercionFailedError:
            return NotImplemented
        result = self.copy()
        for t in other.items():
            result = result._iadd_term(t)
        return result
    __radd__ = __add__

    def __sub__(self, other):
        """Subtract polynomial other from self."""
        ring = self.ring
        try:
            other = ring.convert(other)
        except CoercionFailedError:
            return NotImplemented
        result = self.copy()
        for k, v in other.items():
            result = result._iadd_term((k, -v))
        return result

    def __rsub__(self, other):
        """Substract self from other, with other convertible to the coefficient domain."""
        return (-self).__add__(other)

    def __mul__(self, other):
        """Multiply two polynomials."""
        ring = self.ring
        domain = ring.domain
        zero = ring.zero

        if not other:
            return zero
        if isinstance(other, domain.dtype):
            result = ring.dtype({monom: self[monom]*other for monom in self})
            if not domain.is_Field and not domain.is_IntegerRing:
                result._strip_zero()
            return result

        try:
            other = ring.convert(other)
        except CoercionFailedError:
            return NotImplemented

        if len(other) == 1:
            [(m, c)] = other.items()
            return self.__class__({monom*m: self[monom]*c for monom in self})

        result = zero
        for t in self.items():
            result = result._iadd_poly_term(other, t)
        return result
    __rmul__ = __mul__

    def __pow__(self, n, mod=None):
        """Raise polynomial to power `n`."""
        ring = self.ring

        if not isinstance(n, int) or n < 0:
            raise ValueError('exponent must be a nonnegative integer')

        if not n:
            return ring.one
        if len(self) > 5 or mod:
            return self._pow_generic(n, mod)
        if len(self) == 1:
            [(monom, coeff)] = self.items()
            p = ring.zero
            p[monom**n] = coeff**n
            return p
        if n == 1:
            return self.copy()
        if n == 2:
            return self._square()
        if n == 3:
            return self*self._square()
        return self._pow_multinomial(n)

    def _pow_generic(self, n, mod=None):
        p = self.ring.one
        c = self

        while n:
            if n & 1:
                p *= c
                if mod:
                    p %= mod
                n -= 1

            c = c._square()
            if mod:
                c %= mod
            n //= 2

        return p

    def _pow_multinomial(self, n):
        multinomials = multinomial_coefficients(len(self), n).items()
        ring = self.ring
        zero_monom = ring.zero_monom
        terms = self.items()
        poly = ring.zero

        for multinomial, multinomial_coeff in multinomials:
            product_monom = zero_monom
            product_coeff = multinomial_coeff

            for exp, (monom, coeff) in zip(multinomial, terms):
                if exp:
                    product_monom *= monom**exp
                    product_coeff *= coeff**exp

            poly = poly._iadd_term((product_monom, product_coeff))

        return poly

    def _square(self):
        """Square of a polynomial."""
        ring = self.ring
        p = ring.zero
        keys = list(self)
        for i, k1 in enumerate(keys):
            pk = self[k1]
            for j in range(i):
                k2 = keys[j]
                exp = k1*k2
                p[exp] += pk*self[k2]
        p += p
        for k, v in self.items():
            k2 = k**2
            p[k2] += v**2
        p._strip_zero()
        return p

    def __divmod__(self, other):
        ring = self.ring
        domain = ring.domain

        if not other:
            raise ZeroDivisionError('polynomial division')
        if isinstance(other, ring.dtype):
            (q,), r = self.div([other])
            return q, r
        if isinstance(other, PolyElement):
            if isinstance(domain, PolynomialRing) and domain.ring == other.ring:
                pass
            else:
                return NotImplemented

        try:
            other = ring.domain_new(other)
        except CoercionFailedError:
            return NotImplemented
        return self.quo_ground(other), self.trunc_ground(other)

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __truediv__(self, other):
        ring = self.ring

        if not other:
            raise ZeroDivisionError('polynomial division')
        if isinstance(other, ring.domain.dtype):
            return self.quo_ground(other)

        try:
            other = ring.domain_new(other)
        except CoercionFailedError:
            return NotImplemented
        return self.quo_ground(other)

    def div(self, fv):
        """Division algorithm for multivariate polynomials.

        Parameters
        ==========

        fv : sequence of PolyElement's
            List of divsors.

        Returns
        =======

        (qv, r) : tuple
            Where qv is the sequence of quotients and r is the remainder.

        Notes
        =====

        For multivariate polynomials the remainder is not uniquely
        determined, unless divisors form a Gröbner basis.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)
        >>> f = x**2*y
        >>> f1, f2 = x**2 - y, x*y - 1
        >>> f.div([f1, f2])
        ([y, 0], y**2)
        >>> f.div([f2, f1])
        ([x, 0], x)

        >>> g1, g2 = x - y**2, y**3 - 1
        >>> f.div([g1, g2])[1] == f.div([g2, g1])[1]
        True

        References
        ==========

        * :cite:`Cox2015ideals`, p. 64.

        """
        ring = self.ring
        domain = ring.domain
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
        expvs = [fx.LM for fx in fv]
        while p:
            i = 0
            divoccurred = 0
            while i < s and divoccurred == 0:
                lt = p.leading_term()
                [(expv1, c)] = lt.items()
                expv1 /= expvs[i]
                ci = fv[i][expvs[i]]
                if all(_ >= 0 for _ in expv1) and (domain.is_Field or not c % ci):
                    c = domain.quo(c, ci)
                    qv[i] = qv[i]._iadd_term((expv1, c))
                    p = p._iadd_poly_term(fv[i], (expv1, -c))
                    divoccurred = 1
                    if p and order(p.LM) >= order(lt.LM):
                        raise PolynomialDivisionFailedError(self, fv[i], self.ring)
                else:
                    i += 1
            if not divoccurred:
                expv = p.LM
                r = r._iadd_term((expv, p[expv]))
                del p[expv]
        r._hash = None
        return qv, r

    def exquo(self, other):
        q, r = divmod(self, other)

        if not r:
            return q
        raise ExactQuotientFailedError(self, other)

    def _iadd_term(self, term):
        """Add to self the term inplace.

        If self is a generator -- then just return the sum of the two.

        """
        p1 = self
        if p1.is_generator:
            p1 = p1.copy()
        monom, coeff = term
        coeff += p1[monom]
        if coeff:
            p1[monom] = coeff
        elif monom in p1:
            del p1[monom]
        return p1

    def _iadd_poly_term(self, p2, term):
        """Add inplace to self the product of p2 and term.

        If self is a generator -- then just return the sum of the two.

        """
        p1 = self
        if p1.is_generator:
            p1 = p1.copy()
        monom, coeff = term
        for m, c in p2.items():
            m *= monom
            c *= coeff
            p1 = p1._iadd_term((m, c))
        return p1

    def degree(self, x=0):
        """
        The leading degree in ``x`` or the main variable.

        Note that the degree of 0 is negative floating-point infinity.

        """
        i = self.ring.index(x)
        return max((monom[i] for monom in self), default=-math.inf)

    def tail_degree(self, x=0):
        """
        The tail degree in ``x`` or the main variable.

        Note that the degree of 0 is negative floating-point infinity.

        """
        i = self.ring.index(x)
        return min((monom[i] for monom in self), default=-math.inf)

    def total_degree(self):
        """Returns the total degree."""
        return max((sum(m) for m in self), default=-math.inf)

    def leading_expv(self, order=None):
        """Leading monomial tuple according to the monomial ordering.

        Examples
        ========

        >>> _, x, y, z = ring('x y z', ZZ)
        >>> p = x**4 + x**3*y + x**2*z**2 + z**7
        >>> p.leading_expv()
        (4, 0, 0)

        """
        return self.ring.leading_expv(self, order=order)

    def __getitem__(self, monom, /):
        """Return the coefficient for the given monomial.

        Parameters
        ==========

        monom : Monomial or PolyElement (with ``is_monomial = True``) or 1

        Examples
        ========

        >>> _, x, y, z = ring('x y z', ZZ)
        >>> f = 3*x**2*y - x*y*z + 7*z**3 + 23

        >>> f[x**2*y]
        3
        >>> f[x*y]
        0
        >>> f[1]
        23

        """
        ring = self.ring

        if isinstance(monom, tuple):
            return self.get(monom, ring.domain.zero)
        if monom == 1:
            return self.get(ring.zero_monom, ring.domain.zero)
        if isinstance(monom, ring.dtype) and monom.is_monomial:
            monom, = monom
            return self.get(monom, ring.domain.zero)

        raise TypeError(f'expected a monomial, got {monom}')

    @property
    def LC(self):
        return self[self.LM]

    @property
    def LM(self):
        return self.ring.leading_expv(self)

    @property
    def LT(self):
        expv = self.LM
        return expv, self[expv]

    def leading_term(self, order=None):
        """Leading term as a polynomial element.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)
        >>> (3*x*y + y**2).leading_term()
        3*x*y

        """
        ring = self.ring
        expv = ring.leading_expv(self, order=order)
        return ring.term_new(expv, self[expv])

    def content(self):
        """Returns GCD of polynomial's coefficients."""
        ring = self.ring
        domain = ring.domain
        gcd = domain.gcd

        cont = functools.reduce(gcd, self.values(), domain.zero)

        if not ring.is_normal(self):
            cont = -cont

        return cont

    def primitive(self):
        """Returns content and a primitive polynomial."""
        cont = self.content()
        prim = self.copy()
        if prim:
            prim = prim.quo_ground(cont)
        return cont, prim

    def monic(self):
        """Divides all coefficients by the leading coefficient."""
        if not self:
            return self
        return self.exquo_ground(self.LC)

    def quo_ground(self, x):
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not self or x == 1:
            return self

        if domain.is_Field:
            quo = domain.quo
            p = self.__class__({monom: quo(self[monom], x) for monom in self})
        else:
            p = self.__class__({monom: self[monom]//x for monom in self})

        p._strip_zero()
        return p

    def exquo_ground(self, x):
        domain = self.ring.domain

        if not x:
            raise ZeroDivisionError('polynomial division')
        if not self or x == 1:
            return self

        p = self.__class__({monom: domain.exquo(self[monom], x) for monom in self})
        p._strip_zero()
        return p

    def trunc_ground(self, p):
        if self.ring.domain.is_IntegerRing:
            terms = {}

            for monom, coeff in self.items():
                coeff %= p
                terms[monom] = symmetric_residue(coeff, p)
        else:
            terms = {monom: self[monom] % p for monom in self}

        poly = self.__class__(terms)
        poly._strip_zero()
        return poly

    def _norm(self, norm_func):
        if not self:
            return self.ring.domain.zero
        return norm_func([abs(coeff) for coeff in self.values()])

    def max_norm(self):
        return self._norm(max)

    def l1_norm(self):
        return self._norm(sum)

    def gcd(self, other):
        return self.ring.gcd(self, other)

    def lcm(self, other):
        return self.ring.lcm(self, other)

    def cofactors(self, other):
        return self.ring.cofactors(self, other)

    def terms_gcd(self):
        ring = self.ring

        if not self:
            return (0,)*ring.ngens, self

        G = functools.reduce(Monomial.gcd, self)

        if all(g == 0 for g in G):
            return G, self

        f = ring.zero

        for monom, coeff in self.items():
            f[monom/G] = coeff

        return G, f

    def cancel(self, g, include=True):
        """
        Cancel common factors in a rational function ``f/g``.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)

        >>> (2*x**2 - 2).cancel(x**2 - 2*x + 1)
        (2*x + 2, x - 1)

        """
        f = self
        ring = f.ring
        domain = ring.domain

        if not domain.is_Field or not domain.has_assoc_Ring:
            _, p, q = f.cofactors(g)
            cp, cq = domain.one, domain.one
        else:
            new_ring = ring.clone(domain=domain.ring)

            cq, f = f.clear_denoms(convert=True)
            cp, g = g.clear_denoms(convert=True)

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

        p *= domain(cp)
        q *= domain(cq)

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
        x, = ring.gens[i]
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
        x, = ring.gens[i]
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
        ring = self.ring
        ngens = ring.ngens
        if 0 < (nval := len(values)) <= ngens:
            return self.eval(list(zip(ring.gens, values)))
        raise ValueError(f'expected at least 1 and at most {ngens} values, got {nval}')

    def eval(self, x=0, a=0):
        if isinstance(x, list) and not a:
            (X, a), x = x[0], x[1:]
            f = self.eval(X, a)

            if x:
                return f.eval([(Y.drop(X), a) for (Y, a) in x])
            return f

        return self.compose(x, a).drop(x)

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
                replacements = list(x.items())
            else:
                raise ValueError('expected a generator, value pair a '
                                 'sequence of such pairs')

        replacements = [(ring.index(x), ring(g)) for x, g in replacements]
        replacements.sort(key=lambda k: k[0])

        if ring.is_univariate:
            [(i, g)] = replacements
            acc, d = ring.one, 0
            for monom, coeff in sorted(self.items(), key=lambda x: x[0]):
                n = monom[i]
                acc *= g**(n - d)
                d = n
                poly += acc*coeff
            return poly

        for monom, coeff in self.items():
            monom = list(monom)
            subpoly = ring.one

            for i, g in replacements:
                n, monom[i] = monom[i], 0
                subpoly *= g**n

            monom = Monomial(monom)
            subpoly *= ring.from_terms([(monom, coeff)])
            poly += subpoly

        return poly

    def discriminant(self):
        """Computes discriminant of a polynomial."""
        ring = self.ring

        if (d := self.degree()) <= 0:
            return ring.zero.drop(0)

        s = (-1)**((d*(d - 1)) // 2)
        c = self.eject(*ring.gens[1:]).LC

        return self.resultant(self.diff()) // (c*s)

    def slice(self, m, n, x=0):
        ring = self.ring
        poly = ring.zero
        j = ring.index(x)

        for monom, coeff in self.items():
            if not n > monom[j] >= m:
                if ring.is_univariate:
                    continue
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

        >>> _, x, y = ring('x y', ZZ)

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
            f, g = map(operator.methodcaller('eject', *ring.gens[1:]), (f, g))
            r = f.prem(g)
            return r.inject()

        df = f.degree()
        dg = g.degree()

        if dg < 0:
            raise ZeroDivisionError('polynomial division')

        r, dr = f, df

        if dr < dg:
            return r

        x = ring.gens[0]
        n = df - dg + 1
        lc_g = g.LC

        while dr >= dg:
            lc_r = r.LC
            n -= 1

            r *= lc_g
            r -= g*x**(dr - dg)*lc_r
            dr = r.degree()

        r *= lc_g**n

        return r

    @cacheit
    def resultant(self, other, includePRS=False):
        """
        Computes resultant of two polynomials in `K[X]`.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)

        >>> f = 3*x**2*y - y**3 - 4
        >>> g = x**2 + x*y**3 - 9

        >>> f.resultant(g)
        -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

        """
        ring = self.ring
        domain = ring.domain

        if (not includePRS and query('USE_COLLINS_RESULTANT') and
                (domain.is_IntegerRing or domain.is_RationalField)):
            return ring._collins_resultant(self, other)

        res = ring._primitive_prs(self, other)

        return res if includePRS else res[0]

    def subresultants(self, other):
        """
        Computes subresultant PRS of two polynomials in `K[X]`.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)

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

        >>> _, x = ring('x', QQ)

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> f.half_gcdex(g)
        (-1/5*x + 3/5, x + 1)

        """
        ring = self.ring
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

        >>> _, x = ring('x', QQ)

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> f.gcdex(g)
        (-1/5*x + 3/5, 1/5*x**2 - 6/5*x + 2, x + 1)

        """
        s, h = self.half_gcdex(other)
        t = h - self*s
        t //= other
        return s, t, h

    def sqf_list(self):
        return self.ring.sqf_list(self)

    def sqf_part(self):
        return self.ring.sqf_part(self)

    @property
    def is_squarefree(self):
        return self.ring.is_squarefree(self)

    def sqf_norm(self):
        return self.ring.sqf_norm(self)

    def factor_list(self):
        return self.ring.factor_list(self)
