"""User-friendly public interface to polynomial functions. """

import mpmath

from ..core import (Add, Basic, Derivative, Dummy, E, Expr, I, Integer, Mul,
                    Tuple, oo, preorder_traversal, sympify)
from ..core.compatibility import default_sort_key, iterable
from ..core.decorators import _sympifyit
from ..core.mul import _keep_coeff
from ..core.relational import Relational
from ..domains import FF, QQ, ZZ
from ..logic.boolalg import BooleanAtom
from ..utilities import group, sift
from . import polyoptions as options
from .constructor import construct_domain
from .fglmtools import matrix_fglm
from .groebnertools import groebner as _groebner
from .monomials import Monomial
from .orderings import monomial_key
from .polyerrors import (CoercionFailed, ComputationFailed, DomainError,
                         ExactQuotientFailed, GeneratorsError,
                         GeneratorsNeeded, MultivariatePolynomialError,
                         PolificationFailed, PolynomialError,
                         UnificationFailed)
from .polyutils import (_dict_from_expr, _dict_reorder,
                        _parallel_dict_from_expr, _sort_gens)
from .rationaltools import together
from .rings import PolyElement


__all__ = ('Poly', 'PurePoly', 'poly_from_expr', 'parallel_poly_from_expr',
           'degree', 'degree_list', 'LC', 'LM', 'LT', 'prem',
           'div', 'rem', 'quo', 'exquo', 'half_gcdex', 'gcdex',
           'invert', 'subresultants', 'resultant', 'discriminant', 'cofactors',
           'gcd_list', 'gcd', 'lcm_list', 'lcm', 'terms_gcd', 'trunc',
           'monic', 'content', 'primitive', 'compose', 'decompose', 'sturm',
           'sqf_norm', 'sqf_part', 'sqf_list', 'sqf',
           'factor_list', 'factor', 'intervals', 'refine_root', 'count_roots',
           'real_roots', 'nroots', 'ground_roots', 'nth_power_roots_poly',
           'cancel', 'reduced', 'groebner', 'GroebnerBasis', 'poly')


class Poly(Expr):
    """Generic class for representing polynomial expressions."""

    is_commutative = True
    is_Poly = True

    _op_priority = 10.1

    def __new__(cls, rep, *gens, **args):
        """Create a new polynomial instance out of something useful."""
        opt = options.build_options(gens, args)

        if 'order' in opt:
            raise NotImplementedError("'order' keyword is not implemented yet")

        if iterable(rep, exclude=str):
            if isinstance(rep, dict):
                return cls._from_dict(rep, opt)
            else:
                return cls._from_list(list(rep), opt)
        else:
            rep = sympify(rep)

            if rep.is_Poly:
                return cls._from_poly(rep, opt)
            else:
                return cls._from_expr(rep, opt)

    @classmethod
    def new(cls, rep, *gens):
        """Construct :class:`Poly` instance from raw representation."""
        if not isinstance(rep, PolyElement):
            raise PolynomialError(
                "invalid polynomial representation: %s" % rep)
        elif rep.ring.ngens != len(gens):
            raise PolynomialError("invalid arguments: %s, %s" % (rep, gens))

        obj = Expr.__new__(cls)

        obj.rep = rep
        obj.gens = gens

        return obj

    @classmethod
    def from_dict(cls, rep, *gens, **args):
        """Construct a polynomial from a :class:`dict`."""
        opt = options.build_options(gens, args)
        return cls._from_dict(rep, opt)

    @classmethod
    def from_list(cls, rep, *gens, **args):
        """Construct a polynomial from a :class:`list`."""
        opt = options.build_options(gens, args)
        return cls._from_list(rep, opt)

    @classmethod
    def from_poly(cls, rep, *gens, **args):
        """Construct a polynomial from a polynomial."""
        opt = options.build_options(gens, args)
        return cls._from_poly(rep, opt)

    @classmethod
    def from_expr(cls, rep, *gens, **args):
        """Construct a polynomial from an expression."""
        opt = options.build_options(gens, args)
        return cls._from_expr(rep, opt)

    @classmethod
    def _from_dict(cls, rep, opt):
        """Construct a polynomial from a :class:`dict`."""
        gens = opt.gens

        if not gens:
            raise GeneratorsNeeded(
                "can't initialize from 'dict' without generators")

        domain = opt.domain

        if domain is None:
            domain, rep = construct_domain(rep, opt=opt)
        else:
            for monom, coeff in rep.items():
                rep[monom] = domain.convert(coeff)

        ring = domain.poly_ring(*gens)

        return cls.new(ring.from_dict(rep), *gens)

    @classmethod
    def _from_list(cls, rep, opt):
        """Construct a polynomial from a :class:`list`."""
        gens = opt.gens

        if not gens:
            raise GeneratorsNeeded(
                "can't initialize from 'list' without generators")
        elif len(gens) != 1:
            raise MultivariatePolynomialError(
                "'list' representation not supported")

        domain = opt.domain

        if domain is None:
            domain, rep = construct_domain(rep, opt=opt)
        else:
            rep = list(map(domain.convert, rep))

        ring = domain.poly_ring(*gens)

        return cls.new(ring.from_list(rep), *gens)

    @classmethod
    def _from_poly(cls, rep, opt):
        """Construct a polynomial from a polynomial."""
        if cls != rep.__class__:
            rep = cls.new(rep.rep, *rep.gens)

        gens = opt.gens
        field = opt.field
        domain = opt.domain

        if gens and rep.gens != gens:
            if set(rep.gens) != set(gens):
                return cls._from_expr(rep.as_expr(), opt)
            else:
                rep = rep.reorder(*gens)

        if 'domain' in opt and domain:
            rep = rep.set_domain(domain)
        elif field is True:
            rep = rep.to_field()

        return rep

    @classmethod
    def _from_expr(cls, rep, opt):
        """Construct a polynomial from an expression."""
        rep, opt = _dict_from_expr(rep, opt)
        return cls._from_dict(rep, opt)

    def _hashable_content(self):
        """Allow Diofant to hash Poly instances."""
        return self.rep, self.gens

    def __hash__(self):
        return super().__hash__()

    @property
    def free_symbols(self):
        """
        Free symbols of a polynomial expression.

        Examples
        ========

        >>> Poly(x**2 + 1).free_symbols
        {x}
        >>> Poly(x**2 + y).free_symbols
        {x, y}
        >>> Poly(x**2 + y, x).free_symbols
        {x, y}

        """
        symbols = set()

        for gen in self.gens:
            symbols |= gen.free_symbols

        return symbols | self.free_symbols_in_domain

    @property
    def free_symbols_in_domain(self):
        """
        Free symbols of the domain of ``self``.

        Examples
        ========

        >>> Poly(x**2 + 1).free_symbols_in_domain
        set()
        >>> Poly(x**2 + y).free_symbols_in_domain
        set()
        >>> Poly(x**2 + y, x).free_symbols_in_domain
        {y}

        """
        domain, symbols = self.domain, set()

        if domain.is_Composite:
            for gen in domain.symbols:
                symbols |= gen.free_symbols
        elif domain.is_SymbolicDomain:
            for coeff in self.coeffs():
                symbols |= coeff.free_symbols

        return symbols

    @property
    def args(self):
        """
        Don't mess up with the core.

        Examples
        ========

        >>> Poly(x**2 + 1, x).args
        (x**2 + 1, x)

        """
        return (self.as_expr(),) + self.gens

    @property
    def is_number(self):
        return self.as_expr().is_number

    @property
    def gen(self):
        """
        Return the principal generator.

        Examples
        ========

        >>> Poly(x**2 + 1, x).gen
        x

        """
        return self.gens[0]

    @property
    def domain(self):
        """Get the ground domain of ``self``."""
        return self.rep.ring.domain

    @property
    def zero(self):
        """Return zero polynomial with ``self``'s properties."""
        return self.new(self.rep.ring.zero, *self.gens)

    @property
    def one(self):
        """Return one polynomial with ``self``'s properties."""
        return self.new(self.rep.ring.one, *self.gens)

    @property
    def unit(self):
        """Return unit polynomial with ``self``'s properties."""
        return self.new(self.rep.unit(self.rep.lev, self.domain), *self.gens)

    def unify(self, other):
        """
        Make ``self`` and ``other`` belong to the same domain.

        Examples
        ========

        >>> f, g = Poly(x/2 + 1), Poly(2*x + 1)

        >>> f
        Poly(1/2*x + 1, x, domain='QQ')
        >>> g
        Poly(2*x + 1, x, domain='ZZ')

        >>> F, G = f.unify(g)

        >>> F
        Poly(1/2*x + 1, x, domain='QQ')
        >>> G
        Poly(2*x + 1, x, domain='QQ')

        """
        _, per, F, G = self._unify(other)
        return per(F), per(G)

    def _unify(self, other):
        other = sympify(other)

        if not other.is_Poly:
            try:
                return (self.domain, self.per, self.rep,
                        self.rep.ring(self.domain.convert(other)))
            except CoercionFailed:
                raise UnificationFailed("can't unify %s with %s" % (self, other))

        newring = self.rep.ring.unify(other.rep.ring)
        gens = newring.symbols
        F, G = self.rep.set_ring(newring), other.rep.set_ring(newring)

        cls = self.__class__
        dom = newring.domain

        def per(rep, dom=dom, gens=gens, remove=None):
            if remove is not None:
                gens = gens[:remove] + gens[remove + 1:]

                if not gens:
                    return dom.to_expr(rep)

            return cls.new(rep, *gens)

        return dom, per, F, G

    def per(self, rep, gens=None, remove=None):
        """
        Create a Poly out of the given representation.

        Examples
        ========

        >>> a = Poly(x**2 + 1)
        >>> R = ZZ.poly_ring(x)

        >>> a.per(R.from_dense([ZZ(1), ZZ(1)]), gens=[y])
        Poly(y + 1, y, domain='ZZ')

        """
        if gens is None:
            gens = self.gens

        if remove is not None:
            gens = gens[:remove] + gens[remove + 1:]

            if not gens:
                return self.domain.to_expr(rep)

        return self.__class__.new(rep, *gens)

    def set_domain(self, domain):
        """Set the ground domain of ``self``."""
        opt = options.build_options(self.gens, {'domain': domain})
        newrep = self.rep.set_domain(opt.domain)
        return self.per(newrep)

    def set_modulus(self, modulus):
        """
        Set the modulus of ``self``.

        Examples
        ========

        >>> Poly(5*x**2 + 2*x - 1, x).set_modulus(2)
        Poly(x**2 + 1, x, modulus=2)

        """
        modulus = options.Modulus.preprocess(modulus)
        return self.set_domain(FF(modulus))

    def get_modulus(self):
        """
        Get the modulus of ``self``.

        Examples
        ========

        >>> Poly(x**2 + 1, modulus=2).get_modulus()
        2

        """
        domain = self.domain

        if domain.is_FiniteField:
            return Integer(domain.characteristic)
        else:
            raise PolynomialError("not a polynomial over a Galois field")

    def _eval_subs(self, old, new):
        """Internal implementation of :func:`~diofant.core.basic.Basic.subs`."""
        if old in self.gens:
            if new.is_number:
                return self.eval(old, new)
            else:
                try:
                    return self.replace(old, new)
                except PolynomialError:
                    pass

        return self.as_expr().subs({old: new})

    def exclude(self):
        """
        Remove unnecessary generators from ``self``.

        Examples
        ========

        >>> Poly(a + x, a, b, c, d, x).exclude()
        Poly(a + x, a, x, domain='ZZ')

        """
        rep = self.rep
        if rep.is_ground:
            return self
        for x in rep.ring.symbols:
            try:
                rep = rep.drop(x)
            except ValueError:
                pass

        return self.per(rep, gens=rep.ring.symbols)

    def replace(self, x, y=None):
        """
        Replace ``x`` with ``y`` in generators list.

        Examples
        ========

        >>> Poly(x**2 + 1, x).replace(x, y)
        Poly(y**2 + 1, y, domain='ZZ')

        """
        if y is None:
            if self.is_univariate:
                x, y = self.gen, x
            else:
                raise PolynomialError(
                    "syntax supported only in univariate case")

        if x == y:
            return self

        if x in self.gens and y not in self.gens:
            dom = self.domain

            if not dom.is_Composite or y not in dom.symbols:
                gens = list(self.gens)
                gens[gens.index(x)] = y
                rep = dom.poly_ring(*gens).from_dict(dict(self.rep))
                return self.per(rep, gens=gens)

        raise PolynomialError("can't replace %s with %s in %s" % (x, y, self))

    def reorder(self, *gens, **args):
        """
        Efficiently apply new order of generators.

        Examples
        ========

        >>> Poly(x**2 + x*y**2, x, y).reorder(y, x)
        Poly(y**2*x + x**2, y, x, domain='ZZ')

        """
        opt = options.Options((), args)

        if not gens:
            gens = _sort_gens(self.gens, opt=opt)
        elif set(self.gens) != set(gens):
            raise PolynomialError(
                "generators list can differ only up to order of elements")

        rep = dict(zip(*_dict_reorder(dict(self.rep), self.gens, gens)))

        newring = self.domain.poly_ring(*gens)
        rep = newring.from_dict(rep)

        return self.per(rep, gens=gens)

    def ltrim(self, gen):
        """
        Remove dummy generators from the "left" of ``self``.

        Examples
        ========

        >>> Poly(y**2 + y*z**2, x, y, z).ltrim(y)
        Poly(y**2 + y*z**2, y, z, domain='ZZ')

        """
        rep = self.as_dict(native=True)
        j = self._gen_to_level(gen)
        terms = {}

        for monom, coeff in rep.items():
            monom = monom[j:]

            if monom not in terms:
                terms[monom] = coeff
            else:
                raise PolynomialError("can't left trim %s" % self)

        gens = self.gens[j:]

        newring = self.domain.poly_ring(*gens)
        rep = newring.from_dict(terms)

        return self.new(rep, *gens)

    def has_only_gens(self, *gens):
        """
        Return ``True`` if ``Poly(f, *gens)`` retains ground domain.

        Examples
        ========

        >>> Poly(x*y + 1, x, y, z).has_only_gens(x, y)
        True
        >>> Poly(x*y + z, x, y, z).has_only_gens(x, y)
        False

        """
        indices = set()

        for gen in gens:
            try:
                index = self.gens.index(gen)
            except ValueError:
                raise GeneratorsError(
                    "%s doesn't have %s as generator" % (self, gen))
            else:
                indices.add(index)

        for monom in self.monoms():
            for i, elt in enumerate(monom):
                if i not in indices and elt:
                    return False

        return True

    def to_ring(self):
        """
        Make the ground domain a ring.

        Examples
        ========

        >>> Poly(x**2 + 1, field=True).to_ring()
        Poly(x**2 + 1, x, domain='ZZ')

        """
        return self.set_domain(self.domain.ring)

    def to_field(self):
        """
        Make the ground domain a field.

        Examples
        ========

        >>> Poly(x**2 + 1).to_field()
        Poly(x**2 + 1, x, domain='QQ')

        """
        return self.set_domain(self.domain.field)

    def to_exact(self):
        """
        Make the ground domain exact.

        Examples
        ========

        >>> Poly(x**2 + 1.0).to_exact()
        Poly(x**2 + 1, x, domain='QQ')

        """
        return self.set_domain(self.domain.get_exact())

    def retract(self, field=None):
        """
        Recalculate the ground domain of a polynomial.

        Examples
        ========

        >>> f = Poly(x**2 + 1, domain=QQ.poly_ring(y))
        >>> f
        Poly(x**2 + 1, x, domain='QQ[y]')

        >>> f.retract()
        Poly(x**2 + 1, x, domain='ZZ')
        >>> f.retract(field=True)
        Poly(x**2 + 1, x, domain='QQ')

        """
        dom, rep = construct_domain(self.as_dict(),
                                    field=field,
                                    composite=self.domain.is_Composite or None,
                                    extension=False if self.domain.is_SymbolicDomain else True)
        return self.from_dict(rep, self.gens, domain=dom)

    def slice(self, x, m, n=None):
        """Take a continuous subsequence of terms of ``self``."""
        if n is None:
            j, m, n = 0, x, m
        else:
            j = self._gen_to_level(x)

        m, n = int(m), int(n)

        result = self.rep.slice(m, n, j)
        return self.per(result)

    def coeffs(self, order=None):
        """
        Returns all non-zero coefficients from ``self`` in lex order.

        Examples
        ========

        >>> Poly(x**3 + 2*x + 3, x).coeffs()
        [1, 2, 3]

        See Also
        ========

        all_coeffs
        coeff_monomial

        """
        return [self.domain.to_expr(c) for c in self.rep.coeffs(order=order)]

    def monoms(self, order=None):
        """
        Returns all non-zero monomials from ``self`` in lex order.

        Examples
        ========

        >>> Poly(x**2 + 2*x*y**2 + x*y + 3*y, x, y).monoms()
        [(2, 0), (1, 2), (1, 1), (0, 1)]

        """
        return self.rep.monoms(order=order)

    def terms(self, order=None):
        """
        Returns all non-zero terms from ``self`` in lex order.

        Examples
        ========

        >>> Poly(x**2 + 2*x*y**2 + x*y + 3*y, x, y).terms()
        [((2, 0), 1), ((1, 2), 2), ((1, 1), 1), ((0, 1), 3)]

        """
        return [(m, self.domain.to_expr(c)) for m, c in self.rep.terms(order=order)]

    def all_coeffs(self):
        """
        Returns all coefficients from a univariate polynomial ``self``.

        Examples
        ========

        >>> Poly(x**3 + 2*x - 1, x).all_coeffs()
        [1, 0, 2, -1]

        """
        return [self.domain.to_expr(c) for c in self.rep.all_coeffs()]

    def termwise(self, func, *gens, **args):
        """
        Apply a function to all terms of ``self``.

        Examples
        ========

        >>> def func(k, coeff):
        ...     k = k[0]
        ...     return coeff//10**(2-k)

        >>> Poly(x**2 + 20*x + 400).termwise(func)
        Poly(x**2 + 2*x + 4, x, domain='ZZ')

        """
        terms = {}

        for monom, coeff in self.terms():
            result = func(monom, coeff)

            if isinstance(result, tuple):
                monom, coeff = result
            else:
                coeff = result

            if coeff:
                if monom not in terms:
                    terms[monom] = coeff
                else:
                    raise PolynomialError("%s monomial was generated"
                                          " twice" % monom)

        return self.from_dict(terms, *(gens or self.gens), **args)

    def length(self):
        """
        Returns the number of non-zero terms in ``self``.

        Examples
        ========

        >>> Poly(x**2 + 2*x - 1).length()
        3

        """
        return len(self.as_dict())

    def as_dict(self, native=False):
        """
        Switch to a :class:`dict` representation.

        Examples
        ========

        >>> Poly(x**2 + 2*x*y**2 - y, x, y).as_dict()
        {(0, 1): -1, (1, 2): 2, (2, 0): 1}

        """
        if native:
            return dict(self.rep)
        else:
            return {k: self.domain.to_expr(v) for k, v in self.rep.items()}

    def as_expr(self, *gens):
        """
        Convert a Poly instance to an Expr instance.

        Examples
        ========

        >>> f = Poly(x**2 + 2*x*y**2 - y, x, y)

        >>> f.as_expr()
        x**2 + 2*x*y**2 - y
        >>> f.as_expr({x: 5})
        10*y**2 - y + 25
        >>> f.as_expr(5, 6)
        379

        """
        if not gens:
            gens = self.gens
        elif len(gens) == 1 and isinstance(gens[0], dict):
            mapping = gens[0]
            gens = list(self.gens)

            for gen, value in mapping.items():
                try:
                    index = gens.index(gen)
                except ValueError:
                    raise GeneratorsError(
                        "%s doesn't have %s as generator" % (self, gen))
                else:
                    gens[index] = value

        return self.rep.as_expr(*gens)

    def deflate(self):
        """
        Reduce degree of ``self`` by mapping ``x_i**m`` to ``y_i``.

        Examples
        ========

        >>> Poly(x**6*y**2 + x**3 + 1, x, y).deflate()
        ((3, 2), Poly(x**2*y + x + 1, x, y, domain='ZZ'))

        """
        J, result = self.rep.deflate()
        return J, self.per(result[0])

    def inject(self, front=False):
        """
        Inject ground domain generators into ``self``.

        Examples
        ========

        >>> f = Poly(x**2*y + x*y**3 + x*y + 1, x)

        >>> f.inject()
        Poly(x**2*y + x*y**3 + x*y + 1, x, y, domain='ZZ')
        >>> f.inject(front=True)
        Poly(y**3*x + y*x**2 + y*x + 1, y, x, domain='ZZ')

        """
        dom = self.domain

        if dom.is_Numerical:
            return self

        if front:
            gens = dom.symbols + self.gens
        else:
            gens = self.gens + dom.symbols

        newring = dom.domain.poly_ring(*gens)
        result = newring.from_expr(self.rep.as_expr())

        return self.new(result, *gens)

    def eject(self, *gens):
        """
        Eject selected generators into the ground domain.

        Examples
        ========

        >>> f = Poly(x**2*y + x*y**3 + x*y + 1, x, y)

        >>> f.eject(x)
        Poly(x*y**3 + (x**2 + x)*y + 1, y, domain='ZZ[x]')
        >>> f.eject(y)
        Poly(y*x**2 + (y**3 + y)*x + 1, x, domain='ZZ[y]')

        """
        dom = self.domain

        if not dom.is_Numerical:
            raise DomainError("can't eject generators over %s" % dom)

        result = self.rep.copy()
        result = result.drop_to_ground(*gens)

        return self.new(result, *result.ring.symbols)

    def terms_gcd(self):
        """
        Remove GCD of terms from the polynomial ``self``.

        Examples
        ========

        >>> Poly(x**6*y**2 + x**3*y, x, y).terms_gcd()
        ((3, 1), Poly(x**3*y + 1, x, y, domain='ZZ'))

        """
        J, result = self.rep.terms_gcd()
        return J, self.per(result)

    def quo_ground(self, coeff):
        """
        Quotient of ``self`` by a an element of the ground domain.

        Examples
        ========

        >>> Poly(2*x + 4).quo_ground(2)
        Poly(x + 2, x, domain='ZZ')

        >>> Poly(2*x + 3).quo_ground(2)
        Poly(x + 1, x, domain='ZZ')

        """
        result = self.rep.quo_ground(coeff)
        return self.per(result)

    def exquo_ground(self, coeff):
        """
        Exact quotient of ``self`` by a an element of the ground domain.

        Examples
        ========

        >>> Poly(2*x + 4).exquo_ground(2)
        Poly(x + 2, x, domain='ZZ')

        >>> Poly(2*x + 3).exquo_ground(2)
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2 does not divide 3 in ZZ

        """
        result = self.rep.exquo_ground(coeff)
        return self.per(result)

    def prem(self, other):
        """
        Polynomial pseudo-remainder of ``self`` by ``other``.

        Examples
        ========

        >>> Poly(x**2 + 1, x).prem(Poly(2*x - 4, x))
        Poly(20, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        result = F.prem(G)
        return per(result)

    def div(self, other, auto=True):
        """
        Polynomial division with remainder of ``self`` by ``other``.

        Examples
        ========

        >>> Poly(x**2 + 1, x).div(Poly(2*x - 4, x))
        (Poly(1/2*x + 1, x, domain='QQ'), Poly(5, x, domain='QQ'))

        >>> Poly(x**2 + 1, x).div(Poly(2*x - 4, x), auto=False)
        (Poly(0, x, domain='ZZ'), Poly(x**2 + 1, x, domain='ZZ'))

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.is_Ring and not dom.is_Field:
            F, G = F.set_domain(F.ring.domain.field), G.set_domain(G.ring.domain.field)
            retract = True

        q, r = divmod(F, G)

        if retract:
            try:
                Q, R = q.set_domain(q.ring.domain.ring), r.set_domain(r.ring.domain.ring)
            except CoercionFailed:
                pass
            else:
                q, r = Q, R

        return per(q), per(r)

    def rem(self, other, auto=True):
        """
        Computes the polynomial remainder of ``self`` by ``other``.

        Examples
        ========

        >>> Poly(x**2 + 1, x).rem(Poly(2*x - 4, x))
        Poly(5, x, domain='ZZ')

        >>> Poly(x**2 + 1, x).rem(Poly(2*x - 4, x), auto=False)
        Poly(x**2 + 1, x, domain='ZZ')

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.is_Ring and not dom.is_Field:
            F, G = F.set_domain(F.ring.domain.field), G.set_domain(G.ring.domain.field)
            retract = True

        r = F % G

        if retract:
            try:
                r = r.set_domain(r.ring.domain.ring)
            except CoercionFailed:
                pass

        return per(r)

    def quo(self, other, auto=True):
        """
        Computes polynomial quotient of ``self`` by ``other``.

        Examples
        ========

        >>> Poly(x**2 + 1, x).quo(Poly(2*x - 4, x))
        Poly(1/2*x + 1, x, domain='QQ')

        >>> Poly(x**2 - 1, x).quo(Poly(x - 1, x))
        Poly(x + 1, x, domain='ZZ')

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.is_Ring and not dom.is_Field:
            F, G = F.set_domain(F.ring.domain.field), G.set_domain(G.ring.domain.field)
            retract = True

        q = F // G

        if retract:
            try:
                q = q.set_domain(q.ring.domain.ring)
            except CoercionFailed:
                pass

        return per(q)

    def exquo(self, other, auto=True):
        """
        Computes polynomial exact quotient of ``self`` by ``other``.

        Examples
        ========

        >>> Poly(x**2 - 1, x).exquo(Poly(x - 1, x))
        Poly(x + 1, x, domain='ZZ')

        >>> Poly(x**2 + 1, x).exquo(Poly(2*x - 4, x))
        Traceback (most recent call last):
        ...
        ExactQuotientFailed: 2*x - 4 does not divide x**2 + 1

        """
        dom, per, F, G = self._unify(other)
        retract = False

        if auto and dom.is_Ring and not dom.is_Field:
            F, G = F.set_domain(F.ring.domain.field), G.set_domain(G.ring.domain.field)
            retract = True

        try:
            q = F.exquo(G)
        except ExactQuotientFailed as exc:
            raise exc.new(self.as_expr(), other.as_expr())

        if retract:
            try:
                q = q.set_domain(q.ring.domain.ring)
            except CoercionFailed:
                pass

        return per(q)

    def _gen_to_level(self, gen):
        """Returns level associated with the given generator."""
        if isinstance(gen, int):
            length = len(self.gens)

            if -length <= gen < length:
                if gen < 0:
                    return length + gen
                else:
                    return gen
            else:
                raise PolynomialError("-%s <= gen < %s expected, got %s" %
                                      (length, length, gen))
        else:
            try:
                return self.gens.index(sympify(gen))
            except ValueError:
                raise PolynomialError(
                    "a valid generator expected, got %s" % gen)

    def degree(self, gen=0):
        """
        Returns degree of ``self`` in ``x_j``.

        The degree of 0 is negative infinity.

        Examples
        ========

        >>> Poly(x**2 + y*x + 1, x, y).degree()
        2
        >>> Poly(x**2 + y*x + y, x, y).degree(y)
        1
        >>> Poly(0, x).degree()
        -oo

        """
        j = self._gen_to_level(gen)

        return self.rep.degree(j)

    def degree_list(self):
        """
        Returns a list of degrees of ``self``.

        Examples
        ========

        >>> Poly(x**2 + y*x + 1, x, y).degree_list()
        (2, 1)

        """
        return self.rep.degree_list()

    def total_degree(self):
        """
        Returns the total degree of ``self``.

        Examples
        ========

        >>> Poly(x**2 + y*x + 1, x, y).total_degree()
        2
        >>> Poly(x + y**5, x, y).total_degree()
        5

        """
        return self.rep.total_degree()

    def LC(self, order=None):
        """
        Returns the leading coefficient of ``self``.

        Examples
        ========

        >>> Poly(4*x**3 + 2*x**2 + 3*x, x).LC()
        4

        """
        if order is not None:
            return self.coeffs(order)[0]

        result = self.rep.LC
        return self.domain.to_expr(result)

    def TC(self):
        """
        Returns the trailing coefficient of ``self``.

        Examples
        ========

        >>> Poly(x**3 + 2*x**2 + 3*x, x).TC()
        0

        """
        result = self.rep.ring.dmp_ground_TC(self.rep)
        return self.domain.to_expr(result)

    def EC(self, order=None):
        """
        Returns the last non-zero coefficient of ``self``.

        Examples
        ========

        >>> Poly(x**3 + 2*x**2 + 3*x, x).EC()
        3

        """
        EM = self.EM(order)
        return self.coeff_monomial(tuple(EM))

    def coeff_monomial(self, monom):
        """
        Returns the coefficient of ``monom`` in ``self`` if there, else None.

        Examples
        ========

        >>> p = Poly(24*x*y*exp(8) + 23*x, x, y)

        >>> p.coeff_monomial(x)
        23
        >>> p.coeff_monomial(y)
        0
        >>> p.coeff_monomial(x*y)
        24*E**8
        >>> p.coeff_monomial((1, 1))
        24*E**8

        Note that ``Expr.coeff()`` behaves differently, collecting terms
        if possible; the Poly must be converted to an Expr to use that
        method, however:

        >>> p.as_expr().coeff(x)
        24*E**8*y + 23
        >>> p.as_expr().coeff(y)
        24*E**8*x
        >>> p.as_expr().coeff(x*y)
        24*E**8

        """
        N = Monomial(monom, self.gens)
        if len(N) != len(self.gens):
            raise ValueError('exponent of each generator must be specified')

        result = self.rep.coeff(N)
        return self.domain.to_expr(result)

    def coeff(self, x, n=1, right=False):
        # the semantics of coeff_monomial and Expr.coeff are different;
        # if someone is working with a Poly, they should be aware of the
        # differences and chose the method best suited for the query.
        # Alternatively, a pure-polys method could be written here but
        # at this time the ``right`` keyword would be ignored because Poly
        # doesn't work with non-commutatives.
        raise NotImplementedError(
            'Either convert to Expr with `as_expr` method '
            'to use Expr\'s coeff method or else use the '
            '`coeff_monomial` method of Polys.')

    def LM(self, order=None):
        """
        Returns the leading monomial of ``self``.

        The leading monomial signifies the the monomial having the highest
        power of the principal generator in the polynomial expression.

        Examples
        ========

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).LM()
        x**2*y**0

        """
        LM = (0,)*len(self.gens) if self.is_zero else self.monoms(order)[0]
        return Monomial(LM, self.gens)

    def EM(self, order=None):
        """
        Returns the last non-zero monomial of ``self``.

        Examples
        ========

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).EM()
        x**0*y**1

        """
        EM = (0,)*len(self.gens) if self.is_zero else self.monoms(order)[-1]
        return Monomial(EM, self.gens)

    def LT(self, order=None):
        """
        Returns the leading term of ``self``.

        The leading term signifies the term having the highest power
        of the principal generator in the polynomial expression.

        Examples
        ========

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).LT()
        (x**2*y**0, 4)

        """
        LM = self.LM(order)
        return LM, self.coeff_monomial(tuple(LM))

    def ET(self, order=None):
        """
        Returns the last non-zero term of ``self``.

        Examples
        ========

        >>> Poly(4*x**2 + 2*x*y**2 + x*y + 3*y, x, y).ET()
        (x**0*y**1, 3)

        """
        EM = self.EM(order)
        return EM, self.coeff_monomial(tuple(EM))

    def max_norm(self):
        """
        Returns maximum norm of ``self``.

        Examples
        ========

        >>> Poly(-x**2 + 2*x - 3, x).max_norm()
        3

        """
        result = self.rep.max_norm()
        return self.domain.to_expr(result)

    def l1_norm(self):
        """
        Returns l1 norm of ``self``.

        Examples
        ========

        >>> Poly(-x**2 + 2*x - 3, x).l1_norm()
        6

        """
        result = self.rep.l1_norm()
        return self.domain.to_expr(result)

    def clear_denoms(self, convert=False):
        """
        Clear denominators, but keep the ground domain.

        Examples
        ========

        >>> f = Poly(x/2 + Rational(1, 3))

        >>> f.clear_denoms()
        (6, Poly(3*x + 2, x, domain='QQ'))
        >>> f.clear_denoms(convert=True)
        (6, Poly(3*x + 2, x, domain='ZZ'))

        """
        f = self

        if not f.domain.is_Field:
            return Integer(1), f

        dom = f.domain
        if dom.has_assoc_Ring:
            dom = f.domain.ring

        coeff, result = f.rep.clear_denoms()
        coeff, f = dom.to_expr(coeff), f.per(result)

        if not convert or not dom.has_assoc_Ring:
            return coeff, f
        else:
            return coeff, f.to_ring()

    def rat_clear_denoms(self, other):
        """
        Clear denominators in a rational function ``self/other``.

        Examples
        ========

        >>> f = Poly(x**2/y + 1, x)
        >>> g = Poly(x**3 + y, x)

        >>> p, q = f.rat_clear_denoms(g)

        >>> p
        Poly(x**2 + y, x, domain='ZZ[y]')
        >>> q
        Poly(y*x**3 + y**2, x, domain='ZZ[y]')

        """
        f, g = self, other

        dom, per, f, g = f._unify(g)

        f = per(f)
        g = per(g)

        if not (dom.is_Field and dom.has_assoc_Ring):
            return f, g

        a, f = f.clear_denoms(convert=True)
        b, g = g.clear_denoms(convert=True)

        f *= b
        g *= a

        return f, g

    def integrate(self, *specs, **args):
        """
        Computes indefinite integral of ``self``.

        Examples
        ========

        >>> Poly(x**2 + 2*x + 1, x).integrate()
        Poly(1/3*x**3 + x**2 + x, x, domain='QQ')

        >>> Poly(x*y**2 + x, x, y).integrate((0, 1), (1, 0))
        Poly(1/2*x**2*y**2 + 1/2*x**2, x, y, domain='QQ')

        """
        f = self

        if args.get('auto', True) and f.domain.is_Ring:
            f = f.to_field()

        if not specs:
            return f.per(f.rep.integrate(m=1))

        rep = f.rep

        for spec in specs:
            if type(spec) is tuple:
                gen, m = spec
            else:
                gen, m = spec, 1

            rep = rep.integrate(f._gen_to_level(gen), int(m))

        return f.per(rep)

    def diff(self, *specs, **kwargs):
        """
        Computes partial derivative of ``self``.

        Examples
        ========

        >>> Poly(x**2 + 2*x + 1, x).diff()
        Poly(2*x + 2, x, domain='ZZ')

        >>> Poly(x*y**2 + x, x, y).diff((0, 0), (1, 1))
        Poly(2*x*y, x, y, domain='ZZ')

        """
        if not kwargs.get('evaluate', True):
            return Derivative(self, *specs, **kwargs)

        if not specs:
            return self.per(self.rep.diff())

        rep = self.rep

        for spec in specs:
            if type(spec) is tuple:
                gen, m = spec
            else:
                gen, m = spec, 1

            rep = rep.diff(self._gen_to_level(gen), int(m))

        return self.per(rep)

    _eval_derivative = diff

    def eval(self, x, a=None, auto=True):
        """
        Evaluate ``self`` at ``a`` in the given variable.

        Examples
        ========

        >>> Poly(x**2 + 2*x + 3, x).eval(2)
        11

        >>> Poly(2*x*y + 3*x + y + 2, x, y).eval(x, 2)
        Poly(5*y + 8, y, domain='ZZ')

        >>> f = Poly(2*x*y + 3*x + y + 2*z, x, y, z)

        >>> f.eval({x: 2})
        Poly(5*y + 2*z + 6, y, z, domain='ZZ')
        >>> f.eval({x: 2, y: 5})
        Poly(2*z + 31, z, domain='ZZ')
        >>> f.eval({x: 2, y: 5, z: 7})
        45

        >>> f.eval((2, 5))
        Poly(2*z + 31, z, domain='ZZ')
        >>> f(2, 5)
        Poly(2*z + 31, z, domain='ZZ')

        """
        f = self

        if a is None:
            if isinstance(x, dict):
                mapping = x

                for gen, value in mapping.items():
                    f = f.eval(gen, value)

                return f
            elif isinstance(x, (tuple, list)):
                values = x

                if len(values) > len(f.gens):
                    raise ValueError("too many values provided")

                for gen, value in zip(f.gens, values):
                    f = f.eval(gen, value)

                return f
            else:
                j, a = 0, x
        else:
            j = f._gen_to_level(x)

        try:
            result = f.rep.eval(j, a)
        except CoercionFailed:
            if not auto:
                raise DomainError("can't evaluate at %s in %s" % (a, f.domain))
            else:
                a_domain, [a] = construct_domain([a])
                new_domain = f.domain.unify(a_domain, f.gens)

                f = f.set_domain(new_domain)
                a = new_domain.convert(a, a_domain)

                result = f.rep.eval(j, a)

        return f.per(result, remove=j)

    def __call__(self, *values):
        """
        Evaluate ``self`` at the give values.

        Examples
        ========

        >>> f = Poly(2*x*y + 3*x + y + 2*z, x, y, z)

        >>> f(2)
        Poly(5*y + 2*z + 6, y, z, domain='ZZ')
        >>> f(2, 5)
        Poly(2*z + 31, z, domain='ZZ')
        >>> f(2, 5, 7)
        45

        """
        return self.eval(values)

    def half_gcdex(self, other, auto=True):
        """
        Half extended Euclidean algorithm of ``self`` and ``other``.

        Returns ``(s, h)`` such that ``h = gcd(f, g)`` and ``s*f = h (mod g)``.

        Examples
        ========

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> Poly(f).half_gcdex(Poly(g))
        (Poly(-1/5*x + 3/5, x, domain='QQ'), Poly(x + 1, x, domain='QQ'))

        """
        dom, per, F, G = self._unify(other)

        if auto and dom.is_Ring:
            F, G = F.set_domain(F.ring.domain.field), G.set_domain(G.ring.domain.field)

        s, h = F.half_gcdex(G)
        return per(s), per(h)

    def gcdex(self, other, auto=True):
        """
        Extended Euclidean algorithm of ``self`` and ``other``.

        Returns ``(s, t, h)`` such that ``h = gcd(f, g)`` and ``s*f + t*g = h``.

        Examples
        ========

        >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
        >>> g = x**3 + x**2 - 4*x - 4

        >>> Poly(f).gcdex(Poly(g))
        (Poly(-1/5*x + 3/5, x, domain='QQ'),
         Poly(1/5*x**2 - 6/5*x + 2, x, domain='QQ'),
         Poly(x + 1, x, domain='QQ'))

        """
        dom, per, F, G = self._unify(other)

        if auto and dom.is_Ring:
            F, G = F.set_domain(F.ring.domain.field), G.set_domain(G.ring.domain.field)

        s, t, h = F.gcdex(G)
        return per(s), per(t), per(h)

    def invert(self, other, auto=True):
        """
        Invert ``self`` modulo ``other`` when possible.

        Examples
        ========

        >>> Poly(x**2 - 1, x).invert(Poly(2*x - 1, x))
        Poly(-4/3, x, domain='QQ')

        >>> Poly(x**2 - 1, x).invert(Poly(x - 1, x))
        Traceback (most recent call last):
        ...
        NotInvertible: zero divisor

        """
        dom, per, F, G = self._unify(other)

        if auto and dom.is_Ring:
            F, G = F.set_domain(F.ring.domain.field), G.set_domain(G.ring.domain.field)

        result = F.ring.dup_invert(F, G)
        return per(result)

    def subresultants(self, other):
        """
        Computes the subresultant PRS of ``self`` and ``other``.

        Examples
        ========

        >>> Poly(x**2 + 1, x).subresultants(Poly(x**2 - 1, x))
        [Poly(x**2 + 1, x, domain='ZZ'),
         Poly(x**2 - 1, x, domain='ZZ'),
         Poly(-2, x, domain='ZZ')]

        """
        _, per, F, G = self._unify(other)

        result = F.subresultants(G)
        return list(map(per, result))

    def resultant(self, other, includePRS=False):
        """
        Computes the resultant of ``self`` and ``other`` via PRS.

        If includePRS=True, it includes the subresultant PRS in the result.
        Because the PRS is used to calculate the resultant, this is more
        efficient than calling :func:`subresultants` separately.

        Examples
        ========

        >>> f = Poly(x**2 + 1, x)

        >>> f.resultant(Poly(x**2 - 1, x))
        4
        >>> f.resultant(Poly(x**2 - 1, x), includePRS=True)
        (4, [Poly(x**2 + 1, x, domain='ZZ'), Poly(x**2 - 1, x, domain='ZZ'),
             Poly(-2, x, domain='ZZ')])

        """
        _, per, F, G = self._unify(other)

        if includePRS:
            result, R = F.resultant(G, includePRS=includePRS)
            return per(result, remove=0), list(map(per, R))
        else:
            result = F.resultant(G)
            return per(result, remove=0)

    def discriminant(self):
        """
        Computes the discriminant of ``self``.

        Examples
        ========

        >>> Poly(x**2 + 2*x + 3, x).discriminant()
        -8

        """
        result = self.rep.discriminant()
        return self.per(result, remove=0)

    def dispersionset(self, other=None):
        r"""Compute the *dispersion set* of two polynomials.

        Examples
        ========

        >>> sorted(Poly((x - 3)*(x + 3)).dispersionset())
        [0, 6]

        See Also
        ========

        dispersion
        diofant.polys.dispersion.dispersionset

        """
        from .dispersion import dispersionset
        return dispersionset(self, other)

    def dispersion(self, other=None):
        r"""Compute the *dispersion* of polynomials.

        Examples
        ========

        >>> Poly((x - 3)*(x + 3)).dispersion()
        6

        See Also
        ========

        dispersionset
        diofant.polys.dispersion.dispersion

        References
        ==========

        * :cite:`Man1994disp`
        * :cite:`Koepf98`
        * :cite:`Abramov71rat`
        * :cite:`Man1993indefsum`

        """
        from .dispersion import dispersion
        return dispersion(self, other)

    def cofactors(self, other):
        """
        Returns the GCD of ``self`` and ``other`` and their cofactors.

        For two polynomials ``f`` and ``g`` it returns polynomials
        ``(h, cff, cfg)`` such that ``h = gcd(f, g)``, and ``cff = quo(f, h)``
        and ``cfg = quo(g, h)`` are, so called, cofactors of ``f`` and ``g``.

        Examples
        ========

        >>> Poly(x**2 - 1, x).cofactors(Poly(x**2 - 3*x + 2, x))
        (Poly(x - 1, x, domain='ZZ'),
         Poly(x + 1, x, domain='ZZ'),
         Poly(x - 2, x, domain='ZZ'))

        """
        _, per, F, G = self._unify(other)

        h, cff, cfg = F.cofactors(G)
        return per(h), per(cff), per(cfg)

    def gcd(self, other):
        """
        Returns the polynomial GCD of ``self`` and ``other``.

        Examples
        ========

        >>> Poly(x**2 - 1, x).gcd(Poly(x**2 - 3*x + 2, x))
        Poly(x - 1, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        result = F.gcd(G)
        return per(result)

    def lcm(self, other):
        """
        Returns polynomial LCM of ``self`` and ``other``.

        Examples
        ========

        >>> Poly(x**2 - 1, x).lcm(Poly(x**2 - 3*x + 2, x))
        Poly(x**3 - 2*x**2 - x + 2, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        result = F.lcm(G)
        return per(result)

    def trunc(self, p):
        """
        Reduce ``self`` modulo a constant ``p``.

        Examples
        ========

        >>> Poly(2*x**3 + 3*x**2 + 5*x + 7, x).trunc(3)
        Poly(-x**3 - x + 1, x, domain='ZZ')

        """
        p = self.domain.convert(p)

        result = self.rep.trunc_ground(p)
        return self.per(result)

    def monic(self, auto=True):
        """
        Divides all coefficients by ``LC(f)``.

        Examples
        ========

        >>> Poly(3*x**2 + 6*x + 9).monic()
        Poly(x**2 + 2*x + 3, x, domain='QQ')

        >>> Poly(3*x**2 + 4*x + 2).monic()
        Poly(x**2 + 4/3*x + 2/3, x, domain='QQ')

        """
        f = self

        if auto and f.domain.is_Ring:
            f = f.to_field()

        result = f.rep.monic()
        return f.per(result)

    def content(self):
        """
        Returns the GCD of polynomial coefficients.

        Examples
        ========

        >>> Poly(6*x**2 + 8*x + 12, x).content()
        2

        """
        result = self.rep.content()
        return self.domain.to_expr(result)

    def primitive(self):
        """
        Returns the content and a primitive form of ``self``.

        Examples
        ========

        >>> Poly(2*x**2 + 8*x + 12, x).primitive()
        (2, Poly(x**2 + 4*x + 6, x, domain='ZZ'))

        """
        cont, result = self.rep.primitive()
        return self.domain.to_expr(cont), self.per(result)

    def compose(self, other):
        """
        Computes the functional composition of ``self`` and ``other``.

        Examples
        ========

        >>> Poly(x**2 + x, x).compose(Poly(x - 1, x))
        Poly(x**2 - x, x, domain='ZZ')

        """
        _, per, F, G = self._unify(other)

        result = F.compose(G.ring.gens[0], G)
        return per(result)

    def decompose(self):
        """
        Computes a functional decomposition of ``self``.

        Examples
        ========

        >>> Poly(x**4 + 2*x**3 - x - 1).decompose()
        [Poly(x**2 - x - 1, x, domain='ZZ'), Poly(x**2 + x, x, domain='ZZ')]

        """
        result = self.rep.decompose()
        return list(map(self.per, result))

    def shift(self, a):
        """
        Efficiently compute Taylor shift ``f(x + a)``.

        Examples
        ========

        >>> Poly(x**2 - 2*x + 1, x).shift(2)
        Poly(x**2 + 2*x + 1, x, domain='ZZ')

        """
        result = self.rep.shift(a)
        return self.per(result)

    def sturm(self, auto=True):
        """
        Computes the Sturm sequence of ``self``.

        Examples
        ========

        >>> Poly(x**3 - 2*x**2 + x - 3, x).sturm()
        [Poly(x**3 - 2*x**2 + x - 3, x, domain='QQ'),
         Poly(3*x**2 - 4*x + 1, x, domain='QQ'),
         Poly(2/9*x + 25/9, x, domain='QQ'),
         Poly(-2079/4, x, domain='QQ')]

        """
        f = self

        if auto and f.domain.is_Ring:
            f = f.to_field()

        result = f.rep.sturm()
        return list(map(f.per, result))

    def sqf_norm(self):
        """
        Computes square-free norm of ``self``.

        Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and
        ``r(x) = Norm(g(x))`` is a square-free polynomial over ``K``,
        where ``a`` is the algebraic extension of the ground domain.

        Examples
        ========

        >>> s, f, r = Poly(x**2 + 1, x, extension=[sqrt(3)]).sqf_norm()

        >>> s
        1
        >>> f
        Poly(x**2 - 2*sqrt(3)*x + 4, x, domain='QQ<sqrt(3)>')
        >>> r
        Poly(x**4 - 4*x**2 + 16, x, domain='QQ')

        """
        s, g, r = self.rep.sqf_norm()
        return s, self.per(g), self.per(r)

    def sqf_part(self):
        """
        Computes square-free part of ``self``.

        Examples
        ========

        >>> Poly(x**3 - 3*x - 2, x).sqf_part()
        Poly(x**2 - x - 2, x, domain='ZZ')

        """
        result = self.rep.sqf_part()
        return self.per(result)

    def sqf_list(self):
        """
        Returns a list of square-free factors of ``self``.

        Examples
        ========

        >>> f = 2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16

        >>> Poly(f).sqf_list()
        (2, [(Poly(x + 1, x, domain='ZZ'), 2),
             (Poly(x + 2, x, domain='ZZ'), 3)])

        """
        coeff, factors = self.rep.sqf_list()
        return (self.domain.to_expr(coeff),
                [(self.per(g), k) for g, k in factors])

    def factor_list(self):
        """
        Returns a list of irreducible factors of ``self``.

        Examples
        ========

        >>> f = 2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y

        >>> Poly(f).factor_list()
        (2, [(Poly(x + y, x, y, domain='ZZ'), 1),
             (Poly(x**2 + 1, x, y, domain='ZZ'), 2)])

        """
        try:
            coeff, factors = self.rep.factor_list()
        except DomainError:
            return Integer(1), [(self, 1)]

        return (self.domain.to_expr(coeff),
                [(self.per(g), k) for g, k in factors])

    def intervals(self, all=False, eps=None, inf=None, sup=None, sqf=False):
        """
        Compute isolating intervals for roots of ``self``.

        For real roots the Vincent-Akritas-Strzebonski (VAS) continued fractions method is used.

        References
        ==========

        * :cite:`Alkiviadis2005comp`
        * :cite:`Alkiviadis2008cf`

        Examples
        ========

        >>> Poly(x**2 - 3, x).intervals()
        [((-2, -1), 1), ((1, 2), 1)]
        >>> Poly(x**2 - 3, x).intervals(eps=1e-2)
        [((-26/15, -19/11), 1), ((19/11, 26/15), 1)]

        """
        if eps is not None:
            eps = QQ.convert(eps)

            if eps <= 0:
                raise ValueError("'eps' must be a positive rational")

        if inf is not None:
            inf = QQ.convert(inf)
        if sup is not None:
            sup = QQ.convert(sup)

        R = self.rep.ring

        if self.is_univariate:
            if not all:
                if not sqf:
                    result = R.dup_isolate_real_roots(self.rep, eps=eps,
                                                      inf=inf, sup=sup)
                else:
                    result = R.dup_isolate_real_roots_sqf(self.rep,
                                                          eps=eps, inf=inf,
                                                          sup=sup)
            else:
                if not sqf:
                    result = R.dup_isolate_all_roots(self.rep, eps=eps,
                                                     inf=inf, sup=sup)
                else:
                    result = R.dup_isolate_all_roots_sqf(self.rep,
                                                         eps=eps, inf=inf, sup=sup)
        else:
            raise MultivariatePolynomialError("can't isolate roots of a multivariate polynomial")

        if sqf:
            def _real(interval):
                s, t = interval
                return QQ.to_expr(s), QQ.to_expr(t)

            if not all:
                return list(map(_real, result))

            def _complex(rectangle):
                (u, v), (s, t) = rectangle
                return (QQ.to_expr(u) + I*QQ.to_expr(v),
                        QQ.to_expr(s) + I*QQ.to_expr(t))

            real_part, complex_part = result

            return list(map(_real, real_part)), list(map(_complex, complex_part))
        else:
            def _real(interval):
                (s, t), k = interval
                return (QQ.to_expr(s), QQ.to_expr(t)), k

            if not all:
                return list(map(_real, result))

            def _complex(rectangle):
                ((u, v), (s, t)), k = rectangle
                return ((QQ.to_expr(u) + I*QQ.to_expr(v),
                         QQ.to_expr(s) + I*QQ.to_expr(t)), k)

            real_part, complex_part = result

            return list(map(_real, real_part)), list(map(_complex, complex_part))

    def refine_root(self, s, t, eps=None, steps=None, check_sqf=False):
        """
        Refine an isolating interval of a root to the given precision.

        Examples
        ========

        >>> Poly(x**2 - 3, x).refine_root(1, 2, eps=1e-2)
        (19/11, 26/15)

        """
        if check_sqf and not self.is_squarefree:
            raise PolynomialError("only square-free polynomials supported")

        s, t = QQ.convert(s), QQ.convert(t)

        if eps is not None:
            eps = QQ.convert(eps)

            if eps <= 0:
                raise ValueError("'eps' must be a positive rational")

        if steps is not None:
            steps = int(steps)
        elif eps is None:
            steps = 1

        S, T = self.rep.ring.dup_refine_real_root(self.rep, s, t, eps=eps, steps=steps)
        return QQ.to_expr(S), QQ.to_expr(T)

    def count_roots(self, inf=None, sup=None):
        """
        Return the number of roots of ``self`` in ``[inf, sup]`` interval.

        Examples
        ========

        >>> Poly(x**4 - 4, x).count_roots(-3, 3)
        2
        >>> Poly(x**4 - 4, x).count_roots(0, 1 + 3*I)
        1

        """
        inf_real, sup_real = True, True

        if inf is not None:
            inf = sympify(inf)

            if inf == -oo:
                inf = None
            else:
                re, im = inf.as_real_imag()

                if not im:
                    inf = QQ.convert(inf)
                else:
                    inf, inf_real = tuple(map(QQ.convert, (re, im))), False

        if sup is not None:
            sup = sympify(sup)

            if sup is oo:
                sup = None
            else:
                re, im = sup.as_real_imag()

                if not im:
                    sup = QQ.convert(sup)
                else:
                    sup, sup_real = tuple(map(QQ.convert, (re, im))), False

        if inf_real and sup_real:
            count = self.rep.ring.dup_count_real_roots(self.rep, inf=inf, sup=sup)
        else:
            if inf_real and inf is not None:
                inf = (inf, QQ.zero)

            if sup_real and sup is not None:
                sup = (sup, QQ.zero)

            count = self.rep.ring.dup_count_complex_roots(self.rep, inf=inf, sup=sup)

        return Integer(count)

    def root(self, index, radicals=True):
        """
        Get an indexed root of a polynomial.

        Examples
        ========

        >>> f = Poly(2*x**3 - 7*x**2 + 4*x + 4)

        >>> f.root(0)
        -1/2
        >>> f.root(1)
        2
        >>> f.root(2)
        2
        >>> f.root(3)
        Traceback (most recent call last):
        ...
        IndexError: root index out of [-3, 2] range, got 3

        >>> Poly(x**5 + x + 1).root(0)
        RootOf(x**3 - x**2 + 1, 0)

        """
        from .rootoftools import RootOf
        return RootOf(self, index, radicals=radicals)

    def real_roots(self, multiple=True, radicals=True):
        """
        Return a list of real roots with multiplicities.

        Examples
        ========

        >>> Poly(2*x**3 - 7*x**2 + 4*x + 4).real_roots()
        [-1/2, 2, 2]
        >>> Poly(x**3 + x + 1).real_roots()
        [RootOf(x**3 + x + 1, 0)]

        """
        from .rootoftools import RootOf
        reals = RootOf.real_roots(self, radicals=radicals)

        if multiple:
            return reals
        else:
            return group(reals, multiple=False)

    def all_roots(self, multiple=True, radicals=True):
        """
        Return a list of real and complex roots with multiplicities.

        Examples
        ========

        >>> Poly(2*x**3 - 7*x**2 + 4*x + 4).all_roots()
        [-1/2, 2, 2]
        >>> Poly(x**3 + x + 1).all_roots()
        [RootOf(x**3 + x + 1, 0), RootOf(x**3 + x + 1, 1),
         RootOf(x**3 + x + 1, 2)]

        """
        from .rootoftools import RootOf
        roots = RootOf.all_roots(self, radicals=radicals)

        if multiple:
            return roots
        else:
            return group(roots, multiple=False)

    def nroots(self, n=15, maxsteps=50, cleanup=True):
        """
        Compute numerical approximations of roots of ``self``.

        Parameters
        ==========

        n ... the number of digits to calculate
        maxsteps ... the maximum number of iterations to do

        If the accuracy `n` cannot be reached in `maxsteps`, it will raise an
        exception. You need to rerun with higher maxsteps.

        Examples
        ========

        >>> Poly(x**2 - 3).nroots(n=15)
        [-1.73205080756888, 1.73205080756888]
        >>> Poly(x**2 - 3).nroots(n=30)
        [-1.73205080756887729352744634151, 1.73205080756887729352744634151]

        """
        if self.is_multivariate:
            raise MultivariatePolynomialError(
                "can't compute numerical roots of %s" % self)

        if self.degree() <= 0:
            return []

        # For integer and rational coefficients, convert them to integers only
        # (for accuracy). Otherwise just try to convert the coefficients to
        # mpmath.mpc and raise an exception if the conversion fails.
        if self.domain is ZZ:
            coeffs = [int(coeff) for coeff in self.all_coeffs()]
        elif self.domain is QQ:
            denoms = [coeff.denominator for coeff in self.all_coeffs()]
            from ..core import ilcm
            fac = ilcm(*denoms)
            coeffs = [int(coeff*fac) for coeff in self.all_coeffs()]
        else:
            coeffs = [coeff.evalf(n, strict=False).as_real_imag()
                      for coeff in self.all_coeffs()]
            try:
                coeffs = [mpmath.mpc(*coeff) for coeff in coeffs]
            except TypeError:
                raise DomainError("Numerical domain expected, got %s" %
                                  self.domain)

        dps = mpmath.mp.dps
        mpmath.mp.dps = n

        try:
            # We need to add extra precision to guard against losing accuracy.
            # 10 times the degree of the polynomial seems to work well.
            roots = mpmath.polyroots(coeffs, maxsteps=maxsteps,
                                     cleanup=cleanup, error=False,
                                     extraprec=self.degree()*10)

            # Mpmath puts real roots first, then complex ones (as does all_roots)
            # so we make sure this convention holds here, too.
            roots = list(map(sympify,
                             sorted(roots, key=lambda r: (1 if r.imag else 0, r.real, r.imag))))
        except mpmath.libmp.NoConvergence:
            raise mpmath.libmp.NoConvergence(
                'convergence to root failed; try n < %s or maxsteps > %s' % (
                    n, maxsteps))
        finally:
            mpmath.mp.dps = dps

        return roots

    def ground_roots(self):
        """
        Compute roots of ``self`` by factorization in the ground domain.

        Examples
        ========

        >>> Poly(x**6 - 4*x**4 + 4*x**3 - x**2).ground_roots()
        {0: 2, 1: 2}

        """
        if self.is_multivariate:
            raise MultivariatePolynomialError(
                "can't compute ground roots of %s" % self)

        roots = {}

        for factor, k in self.factor_list()[1]:
            if factor.is_linear:
                a, b = factor.all_coeffs()
                roots[-b/a] = k

        return roots

    def nth_power_roots_poly(self, n):
        """
        Construct a polynomial with n-th powers of roots of ``self``.

        Examples
        ========

        >>> f = Poly(x**4 - x**2 + 1)

        >>> f.nth_power_roots_poly(2)
        Poly(x**4 - 2*x**3 + 3*x**2 - 2*x + 1, x, domain='ZZ')
        >>> f.nth_power_roots_poly(3)
        Poly(x**4 + 2*x**2 + 1, x, domain='ZZ')
        >>> f.nth_power_roots_poly(4)
        Poly(x**4 + 2*x**3 + 3*x**2 + 2*x + 1, x, domain='ZZ')
        >>> f.nth_power_roots_poly(12)
        Poly(x**4 - 4*x**3 + 6*x**2 - 4*x + 1, x, domain='ZZ')

        """
        if self.is_multivariate:
            raise MultivariatePolynomialError("must be a univariate polynomial")

        N = sympify(n)

        if N.is_Integer and N >= 1:
            n = int(N)
        else:
            raise ValueError("'n' must an integer and n >= 1, got %s" % n)

        x = self.gen
        t = Dummy('t')

        r = self.resultant(self.__class__.from_expr(x**n - t, x, t))

        return r.replace(t, x)

    def cancel(self, other, include=False):
        """
        Cancel common factors in a rational function ``self/other``.

        Examples
        ========

        >>> Poly(2*x**2 - 2, x).cancel(Poly(x**2 - 2*x + 1, x))
        (1, Poly(2*x + 2, x, domain='ZZ'), Poly(x - 1, x, domain='ZZ'))

        >>> Poly(2*x**2 - 2, x).cancel(Poly(x**2 - 2*x + 1, x), include=True)
        (Poly(2*x + 2, x, domain='ZZ'), Poly(x - 1, x, domain='ZZ'))

        """
        dom, per, F, G = self._unify(other)

        result = F.cancel(G, include=include)

        if not include:
            if dom.has_assoc_Ring:
                dom = dom.ring

            cp, cq, p, q = result

            cp = dom.to_expr(cp)
            cq = dom.to_expr(cq)

            return cp/cq, per(p), per(q)
        else:
            return tuple(map(per, result))

    @property
    def is_zero(self):
        """
        Returns ``True`` if ``self`` is a zero polynomial.

        Examples
        ========

        >>> Poly(0, x).is_zero
        True
        >>> Poly(1, x).is_zero
        False

        """
        return self.rep.is_zero

    @property
    def is_one(self):
        """
        Returns ``True`` if ``self`` is a unit polynomial.

        Examples
        ========

        >>> Poly(0, x).is_one
        False
        >>> Poly(1, x).is_one
        True

        """
        return self.rep.is_one

    @property
    def is_squarefree(self):
        """
        Returns ``True`` if ``self`` is a square-free polynomial.

        Examples
        ========

        >>> Poly(x**2 - 2*x + 1, x).is_squarefree
        False
        >>> Poly(x**2 - 1, x).is_squarefree
        True

        """
        return self.rep.is_squarefree

    @property
    def is_monic(self):
        """
        Returns ``True`` if the leading coefficient of ``self`` is one.

        Examples
        ========

        >>> Poly(x + 2, x).is_monic
        True
        >>> Poly(2*x + 2, x).is_monic
        False

        """
        return self.rep.is_monic

    @property
    def is_primitive(self):
        """
        Returns ``True`` if GCD of the coefficients of ``self`` is one.

        Examples
        ========

        >>> Poly(2*x**2 + 6*x + 12, x).is_primitive
        False
        >>> Poly(x**2 + 3*x + 6, x).is_primitive
        True

        """
        return self.rep.is_primitive

    @property
    def is_ground(self):
        """
        Returns ``True`` if ``self`` is an element of the ground domain.

        Examples
        ========

        >>> Poly(x, x).is_ground
        False
        >>> Poly(2, x).is_ground
        True
        >>> Poly(y, x).is_ground
        True

        """
        return self.rep.is_ground

    @property
    def is_linear(self):
        """
        Returns ``True`` if ``self`` is linear in all its variables.

        Examples
        ========

        >>> Poly(x + y + 2, x, y).is_linear
        True
        >>> Poly(x*y + 2, x, y).is_linear
        False

        """
        return self.rep.is_linear

    @property
    def is_quadratic(self):
        """
        Returns ``True`` if ``self`` is quadratic in all its variables.

        Examples
        ========

        >>> Poly(x*y + 2, x, y).is_quadratic
        True
        >>> Poly(x*y**2 + 2, x, y).is_quadratic
        False

        """
        return self.rep.is_quadratic

    @property
    def is_term(self):
        """
        Returns ``True`` if ``self`` is zero or has only one term.

        Examples
        ========

        >>> Poly(3*x**2, x).is_term
        True
        >>> Poly(3*x**2 + 1, x).is_term
        False

        """
        return self.rep.is_term

    @property
    def is_homogeneous(self):
        """
        Returns ``True`` if ``self`` is a homogeneous polynomial.

        A homogeneous polynomial is a polynomial whose all monomials with
        non-zero coefficients have the same total degree.

        Examples
        ========

        >>> Poly(x**2 + x*y, x, y).is_homogeneous
        True
        >>> Poly(x**3 + x*y, x, y).is_homogeneous
        False

        """
        return self.rep.is_homogeneous

    @property
    def is_irreducible(self):
        """
        Returns ``True`` if ``self`` has no factors over its domain.

        Examples
        ========

        >>> Poly(x**2 + x + 1, x, modulus=2).is_irreducible
        True
        >>> Poly(x**2 + 1, x, modulus=2).is_irreducible
        False

        """
        return self.rep.is_irreducible

    @property
    def is_univariate(self):
        """
        Returns ``True`` if ``self`` is a univariate polynomial.

        Examples
        ========

        >>> Poly(x**2 + x + 1, x).is_univariate
        True
        >>> Poly(x*y**2 + x*y + 1, x, y).is_univariate
        False
        >>> Poly(x*y**2 + x*y + 1, x).is_univariate
        True
        >>> Poly(x**2 + x + 1, x, y).is_univariate
        False

        """
        return len(self.gens) == 1

    @property
    def is_multivariate(self):
        """
        Returns ``True`` if ``self`` is a multivariate polynomial.

        Examples
        ========

        >>> Poly(x**2 + x + 1, x).is_multivariate
        False
        >>> Poly(x*y**2 + x*y + 1, x, y).is_multivariate
        True
        >>> Poly(x*y**2 + x*y + 1, x).is_multivariate
        False
        >>> Poly(x**2 + x + 1, x, y).is_multivariate
        True

        """
        return len(self.gens) != 1

    @property
    def is_cyclotomic(self):
        """
        Returns ``True`` if ``self`` is a cyclotomic polynomial.

        Examples
        ========

        >>> f = x**16 + x**14 - x**10 + x**8 - x**6 + x**2 + 1
        >>> Poly(f).is_cyclotomic
        False

        >>> g = x**16 + x**14 - x**10 - x**8 - x**6 + x**2 + 1
        >>> Poly(g).is_cyclotomic
        True

        """
        return self.rep.is_cyclotomic

    def __abs__(self):
        """
        Make all coefficients in ``self`` positive.

        Examples
        ========

        >>> abs(Poly(x**2 - 1, x))
        Poly(x**2 + 1, x, domain='ZZ')

        """
        result = abs(self.rep)
        return self.per(result)

    def __neg__(self):
        """
        Negate all coefficients in ``self``.

        Examples
        ========

        >>> -Poly(x**2 - 1, x)
        Poly(-x**2 + 1, x, domain='ZZ')

        """
        result = -self.rep
        return self.per(result)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if not other.is_Poly:
            try:
                other = self.__class__(other, *self.gens)
            except PolynomialError:
                return self.as_expr() + other

        _, per, F, G = self._unify(other)
        result = F + G
        return per(result)

    @_sympifyit('other', NotImplemented)
    def __radd__(self, other):
        try:
            other = self.__class__(other, *self.gens)
        except PolynomialError:
            return other + self.as_expr()

        return other + self

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if not other.is_Poly:
            try:
                other = self.__class__(other, *self.gens)
            except PolynomialError:
                return self.as_expr() - other

        _, per, F, G = self._unify(other)
        result = F - G
        return per(result)

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        try:
            other = self.__class__(other, *self.gens)
        except PolynomialError:
            return other - self.as_expr()

        return other - self

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if not other.is_Poly:
            try:
                other = self.__class__(other, *self.gens)
            except PolynomialError:
                return self.as_expr()*other

        _, per, F, G = self._unify(other)
        result = F * G
        return per(result)

    @_sympifyit('other', NotImplemented)
    def __rmul__(self, other):
        try:
            other = self.__class__(other, *self.gens)
        except PolynomialError:
            return other*self.as_expr()

        return other*self

    @_sympifyit('n', NotImplemented)
    def __pow__(self, n):
        if n.is_Integer and n >= 0:
            n = int(n)
            result = self.rep**n
            return self.per(result)
        else:
            return self.as_expr()**n

    @_sympifyit('other', NotImplemented)
    def __divmod__(self, other):
        if not other.is_Poly:
            other = self.__class__(other, *self.gens)

        return self.div(other)

    @_sympifyit('other', NotImplemented)
    def __rdivmod__(self, other):
        other = self.__class__(other, *self.gens)

        return other.div(self)

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if not other.is_Poly:
            other = self.__class__(other, *self.gens)

        return self.rem(other)

    @_sympifyit('other', NotImplemented)
    def __rmod__(self, other):
        other = self.__class__(other, *self.gens)

        return other.rem(self)

    @_sympifyit('other', NotImplemented)
    def __floordiv__(self, other):
        if not other.is_Poly:
            other = self.__class__(other, *self.gens)

        return self.quo(other)

    @_sympifyit('other', NotImplemented)
    def __rfloordiv__(self, other):
        other = self.__class__(other, *self.gens)

        return other.quo(self)

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        return self.as_expr()/other.as_expr()

    @_sympifyit('other', NotImplemented)
    def __eq__(self, other):
        f, g = self, other

        if not g.is_Poly:
            try:
                g = f.__class__(g, f.gens, domain=f.domain)
            except (PolynomialError, DomainError, CoercionFailed):
                return False

        if f.gens != g.gens:
            return False

        if f.domain != g.domain:
            try:
                dom = f.domain.unify(g.domain, f.gens)
            except UnificationFailed:  # pragma: no cover
                return NotImplemented

            f = f.set_domain(dom)
            g = g.set_domain(dom)

        return f.rep == g.rep

    def __bool__(self):
        return not self.is_zero


class PurePoly(Poly):
    """Class for representing pure polynomials."""

    def _hashable_content(self):
        """Allow Diofant to hash Poly instances."""
        return self.domain, frozenset(self.rep.items())

    def __hash__(self):
        return super().__hash__()

    @property
    def free_symbols(self):
        """
        Free symbols of a polynomial.

        Examples
        ========

        >>> PurePoly(x**2 + 1).free_symbols
        set()
        >>> PurePoly(x**2 + y).free_symbols
        set()
        >>> PurePoly(x**2 + y, x).free_symbols
        {y}

        """
        return self.free_symbols_in_domain

    @_sympifyit('other', NotImplemented)
    def __eq__(self, other):
        f, g = self, other

        if not g.is_Poly:
            try:
                g = f.__class__(g, f.gens, domain=f.domain)
            except (PolynomialError, DomainError, CoercionFailed):
                return False

        if len(f.gens) != len(g.gens):
            return False

        if f.domain != g.domain:
            try:
                dom = f.domain.unify(g.domain, f.gens)
            except UnificationFailed:  # pragma: no cover
                return NotImplemented

            f = f.set_domain(dom)
            g = g.set_domain(dom)

        return f.rep.items() == g.rep.items()

    def _unify(self, other):
        other = sympify(other)

        if not other.is_Poly:
            try:
                return (self.domain, self.per, self.rep,
                        self.rep.ring(self.domain.convert(other)))
            except CoercionFailed:
                raise UnificationFailed("can't unify %s with %s" % (self, other))

        if len(self.gens) != len(other.gens):
            raise UnificationFailed("can't unify %s with %s" % (self, other))

        newring = self.rep.ring.unify(other.rep.ring)
        gens = newring.symbols
        F, G = self.rep.set_ring(newring), other.rep.set_ring(newring)

        cls = self.__class__
        dom = newring.domain

        def per(rep, dom=dom, gens=gens, remove=None):
            if remove is not None:
                gens = gens[:remove] + gens[remove + 1:]

                if not gens:
                    return dom.to_expr(rep)

            return cls.new(rep, *gens)

        return dom, per, F, G


def poly_from_expr(expr, *gens, **args):
    """Construct a polynomial from an expression."""
    opt = options.build_options(gens, args)
    return _poly_from_expr(expr, opt)


def _poly_from_expr(expr, opt):
    """Construct a polynomial from an expression."""
    orig, expr = expr, sympify(expr)

    if not isinstance(expr, Basic):
        raise PolificationFailed(opt, orig, expr)
    elif expr.is_Poly:
        poly = expr.__class__._from_poly(expr, opt)

        opt.gens = poly.gens
        opt.domain = poly.domain

        if opt.polys is None:
            opt.polys = True

        return poly, opt
    elif opt.expand:
        expr = expr.expand()

    try:
        rep, opt = _dict_from_expr(expr, opt)
    except GeneratorsNeeded:
        raise PolificationFailed(opt, orig, expr)

    monoms, coeffs = list(zip(*list(rep.items())))
    domain = opt.domain

    if domain is None:
        opt.domain, coeffs = construct_domain(coeffs, opt=opt)
    else:
        coeffs = list(map(domain.convert, coeffs))

    rep = dict(zip(monoms, coeffs))
    poly = Poly._from_dict(rep, opt)

    if opt.polys is None:
        opt.polys = False

    return poly, opt


def parallel_poly_from_expr(exprs, *gens, **args):
    """Construct polynomials from expressions."""
    opt = options.build_options(gens, args)
    return _parallel_poly_from_expr(exprs, opt)


def _parallel_poly_from_expr(exprs, opt):
    """Construct polynomials from expressions."""
    from ..functions import Piecewise

    if len(exprs) == 2:
        f, g = exprs

        if isinstance(f, Poly) and isinstance(g, Poly):
            f = f.__class__._from_poly(f, opt)
            g = g.__class__._from_poly(g, opt)

            f, g = f.unify(g)

            opt.gens = f.gens
            opt.domain = f.domain

            if opt.polys is None:
                opt.polys = True

            return [f, g], opt

    origs, exprs = list(exprs), []
    _exprs, _polys = [], []

    failed = False

    for i, expr in enumerate(origs):
        expr = sympify(expr)

        if isinstance(expr, Basic):
            if expr.is_Poly:
                _polys.append(i)
            else:
                _exprs.append(i)

                if opt.expand:
                    expr = expr.expand()
        else:
            failed = True

        exprs.append(expr)

    if failed:
        raise PolificationFailed(opt, origs, exprs, True)

    if _polys:
        # XXX: this is a temporary solution
        for i in _polys:
            exprs[i] = exprs[i].as_expr()

    try:
        reps, opt = _parallel_dict_from_expr(exprs, opt)
    except GeneratorsNeeded:
        raise PolificationFailed(opt, origs, exprs, True)

    for k in opt.gens:
        if isinstance(k, Piecewise):
            raise PolynomialError("Piecewise generators do not make sense")

    coeffs_list, lengths = [], []

    all_monoms = []
    all_coeffs = []

    for rep in reps:
        monoms, coeffs = list(zip(*list(rep.items())))

        coeffs_list.extend(coeffs)
        all_monoms.append(monoms)

        lengths.append(len(coeffs))

    domain = opt.domain

    if domain is None:
        opt.domain, coeffs_list = construct_domain(coeffs_list, opt=opt)
    else:
        coeffs_list = list(map(domain.convert, coeffs_list))

    for k in lengths:
        all_coeffs.append(coeffs_list[:k])
        coeffs_list = coeffs_list[k:]

    polys = []

    for monoms, coeffs in zip(all_monoms, all_coeffs):
        rep = dict(zip(monoms, coeffs))
        poly = Poly._from_dict(rep, opt)
        polys.append(poly)

    if opt.polys is None:
        opt.polys = bool(_polys)

    return polys, opt


def degree(f, *gens, **args):
    """
    Return the degree of ``f`` in the given variable.

    The degree of 0 is negative infinity.

    Examples
    ========

    >>> degree(x**2 + y*x + 1, gen=x)
    2
    >>> degree(x**2 + y*x + 1, gen=y)
    1
    >>> degree(0, x)
    -oo

    """
    options.allowed_flags(args, ['gen', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('degree', 1, exc)

    return sympify(F.degree(opt.gen))


def degree_list(f, *gens, **args):
    """
    Return a list of degrees of ``f`` in all variables.

    Examples
    ========

    >>> degree_list(x**2 + y*x + 1)
    (2, 1)

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('degree_list', 1, exc)

    degrees = F.degree_list()

    return tuple(map(Integer, degrees))


def LC(f, *gens, **args):
    """
    Return the leading coefficient of ``f``.

    Examples
    ========

    >>> LC(4*x**2 + 2*x*y**2 + x*y + 3*y)
    4

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('LC', 1, exc)

    return F.LC(order=opt.order)


def LM(f, *gens, **args):
    """
    Return the leading monomial of ``f``.

    Examples
    ========

    >>> LM(4*x**2 + 2*x*y**2 + x*y + 3*y)
    x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('LM', 1, exc)

    monom = F.LM(order=opt.order)
    return monom.as_expr()


def LT(f, *gens, **args):
    """
    Return the leading term of ``f``.

    Examples
    ========

    >>> LT(4*x**2 + 2*x*y**2 + x*y + 3*y)
    4*x**2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('LT', 1, exc)

    monom, coeff = F.LT(order=opt.order)
    return coeff*monom.as_expr()


def prem(f, g, *gens, **args):
    """
    Compute polynomial pseudo-remainder of ``f`` and ``g``.

    Examples
    ========

    >>> prem(x**2 + 1, 2*x - 4)
    20

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('prem', 2, exc)

    r = F.prem(G)

    if not opt.polys:
        return r.as_expr()
    else:
        return r


def div(f, g, *gens, **args):
    """
    Compute polynomial division of ``f`` and ``g``.

    Examples
    ========

    >>> div(x**2 + 1, 2*x - 4, field=False)
    (0, x**2 + 1)
    >>> div(x**2 + 1, 2*x - 4)
    (x/2 + 1, 5)

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('div', 2, exc)

    q, r = F.div(G, auto=opt.auto)

    if not opt.polys:
        return q.as_expr(), r.as_expr()
    else:
        return q, r


def rem(f, g, *gens, **args):
    """
    Compute polynomial remainder of ``f`` and ``g``.

    Examples
    ========

    >>> rem(x**2 + 1, 2*x - 4, field=False)
    x**2 + 1
    >>> rem(x**2 + 1, 2*x - 4)
    5

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('rem', 2, exc)

    r = F.rem(G, auto=opt.auto)

    if not opt.polys:
        return r.as_expr()
    else:
        return r


def quo(f, g, *gens, **args):
    """
    Compute polynomial quotient of ``f`` and ``g``.

    Examples
    ========

    >>> quo(x**2 + 1, 2*x - 4)
    x/2 + 1
    >>> quo(x**2 - 1, x - 1)
    x + 1

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('quo', 2, exc)

    q = F.quo(G, auto=opt.auto)

    if not opt.polys:
        return q.as_expr()
    else:
        return q


def exquo(f, g, *gens, **args):
    """
    Compute polynomial exact quotient of ``f`` and ``g``.

    Examples
    ========

    >>> exquo(x**2 - 1, x - 1)
    x + 1

    >>> exquo(x**2 + 1, 2*x - 4)
    Traceback (most recent call last):
    ...
    ExactQuotientFailed: 2*x - 4 does not divide x**2 + 1

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('exquo', 2, exc)

    q = F.exquo(G, auto=opt.auto)

    if not opt.polys:
        return q.as_expr()
    else:
        return q


def half_gcdex(f, g, *gens, **args):
    """
    Half extended Euclidean algorithm of ``f`` and ``g``.

    Returns ``(s, h)`` such that ``h = gcd(f, g)`` and ``s*f = h (mod g)``.

    Examples
    ========

    >>> half_gcdex(x**4 - 2*x**3 - 6*x**2 + 12*x + 15, x**3 + x**2 - 4*x - 4)
    (-x/5 + 3/5, x + 1)

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        domain, (a, b) = construct_domain(exc.exprs)

        s, h = domain.half_gcdex(a, b)
        return domain.to_expr(s), domain.to_expr(h)

    s, h = F.half_gcdex(G, auto=opt.auto)

    if not opt.polys:
        return s.as_expr(), h.as_expr()
    else:
        return s, h


def gcdex(f, g, *gens, **args):
    """
    Extended Euclidean algorithm of ``f`` and ``g``.

    Returns ``(s, t, h)`` such that ``h = gcd(f, g)`` and ``s*f + t*g = h``.

    Examples
    ========

    >>> gcdex(x**4 - 2*x**3 - 6*x**2 + 12*x + 15, x**3 + x**2 - 4*x - 4)
    (-x/5 + 3/5, x**2/5 - 6*x/5 + 2, x + 1)

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        domain, (a, b) = construct_domain(exc.exprs)

        s, t, h = domain.gcdex(a, b)
        return domain.to_expr(s), domain.to_expr(t), domain.to_expr(h)

    s, t, h = F.gcdex(G, auto=opt.auto)

    if not opt.polys:
        return s.as_expr(), t.as_expr(), h.as_expr()
    else:
        return s, t, h


def invert(f, g, *gens, **args):
    """
    Invert ``f`` modulo ``g`` when possible.

    Examples
    ========

    >>> invert(x**2 - 1, 2*x - 1)
    -4/3

    >>> invert(x**2 - 1, x - 1)
    Traceback (most recent call last):
    ...
    NotInvertible: zero divisor

    For more efficient inversion of Rationals,
    use the ``mod_inverse`` function:

    >>> mod_inverse(3, 5)
    2
    >>> (Integer(2)/5).invert(Integer(7)/3)
    5/2

    See Also
    ========

    diofant.core.numbers.mod_inverse

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        domain, (a, b) = construct_domain(exc.exprs)

        return domain.to_expr(domain.invert(a, b))

    h = F.invert(G, auto=opt.auto)

    if not opt.polys:
        return h.as_expr()
    else:
        return h


def subresultants(f, g, *gens, **args):
    """
    Compute subresultant PRS of ``f`` and ``g``.

    Examples
    ========

    >>> subresultants(x**2 + 1, x**2 - 1)
    [x**2 + 1, x**2 - 1, -2]

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('subresultants', 2, exc)

    result = F.subresultants(G)

    if not opt.polys:
        return [r.as_expr() for r in result]
    else:
        return result


def resultant(f, g, *gens, **args):
    """
    Compute resultant of ``f`` and ``g``.

    Examples
    ========

    >>> resultant(x**2 + 1, x**2 - 1)
    4

    """
    includePRS = args.pop('includePRS', False)
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('resultant', 2, exc)

    if includePRS:
        result, R = F.resultant(G, includePRS=includePRS)
    else:
        result = F.resultant(G)

    if not opt.polys:
        if includePRS:
            return result.as_expr(), [r.as_expr() for r in R]
        return result.as_expr()
    else:
        if includePRS:
            return result, R
        return result


def discriminant(f, *gens, **args):
    """
    Compute discriminant of ``f``.

    Examples
    ========

    >>> discriminant(x**2 + 2*x + 3)
    -8

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('discriminant', 1, exc)

    result = F.discriminant()

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def cofactors(f, g, *gens, **args):
    """
    Compute GCD and cofactors of ``f`` and ``g``.

    Returns polynomials ``(h, cff, cfg)`` such that ``h = gcd(f, g)``, and
    ``cff = quo(f, h)`` and ``cfg = quo(g, h)`` are, so called, cofactors
    of ``f`` and ``g``.

    Examples
    ========

    >>> cofactors(x**2 - 1, x**2 - 3*x + 2)
    (x - 1, x + 1, x - 2)

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        domain, (a, b) = construct_domain(exc.exprs)

        h, cff, cfg = domain.cofactors(a, b)
        return tuple(map(domain.to_expr, (h, cff, cfg)))

    h, cff, cfg = F.cofactors(G)

    if not opt.polys:
        return h.as_expr(), cff.as_expr(), cfg.as_expr()
    else:
        return h, cff, cfg


def gcd_list(seq, *gens, **args):
    """
    Compute GCD of a list of polynomials.

    Examples
    ========

    >>> gcd_list([x**3 - 1, x**2 - 1, x**2 - 3*x + 2])
    x - 1

    """
    seq = sympify(seq)

    def try_non_polynomial_gcd(seq):
        if not gens and not args:
            domain, numbers = construct_domain(seq)

            if not numbers:
                return domain.zero
            elif domain.is_Numerical:
                result, numbers = numbers[0], numbers[1:]

                for number in numbers:
                    result = domain.gcd(result, number)

                    if result == domain.one:
                        break

                return domain.to_expr(result)

    result = try_non_polynomial_gcd(seq)

    if result is not None:
        return result

    options.allowed_flags(args, ['polys'])

    try:
        polys, opt = parallel_poly_from_expr(seq, *gens, **args)
    except PolificationFailed as exc:
        result = try_non_polynomial_gcd(exc.exprs)

        if result is not None:
            return result
        else:
            raise ComputationFailed('gcd_list', len(seq), exc)

    if not polys:
        if not opt.polys:
            return Integer(0)
        else:
            return Poly(0, opt=opt)

    result, polys = polys[0], polys[1:]

    for poly in polys:
        result = result.gcd(poly)

        if result.is_one:
            break

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def gcd(f, g, *gens, **args):
    """
    Compute GCD of ``f`` and ``g``.

    Examples
    ========

    >>> gcd(x**2 - 1, x**2 - 3*x + 2)
    x - 1

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        domain, (a, b) = construct_domain(exc.exprs)
        return domain.to_expr(domain.gcd(a, b))

    result = F.gcd(G)

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def lcm_list(seq, *gens, **args):
    """
    Compute LCM of a list of polynomials.

    Examples
    ========

    >>> lcm_list([x**3 - 1, x**2 - 1, x**2 - 3*x + 2])
    x**5 - x**4 - 2*x**3 - x**2 + x + 2

    """
    seq = sympify(seq)

    def try_non_polynomial_lcm(seq):
        if not gens and not args:
            domain, numbers = construct_domain(seq)

            if not numbers:
                return domain.one
            elif domain.is_Numerical:
                result, numbers = numbers[0], numbers[1:]

                for number in numbers:
                    result = domain.lcm(result, number)

                return domain.to_expr(result)

    result = try_non_polynomial_lcm(seq)

    if result is not None:
        return result

    options.allowed_flags(args, ['polys'])

    try:
        polys, opt = parallel_poly_from_expr(seq, *gens, **args)
    except PolificationFailed as exc:
        result = try_non_polynomial_lcm(exc.exprs)

        if result is not None:
            return result
        else:
            raise ComputationFailed('lcm_list', len(seq), exc)

    if not polys:
        if not opt.polys:
            return Integer(1)
        else:
            return Poly(1, opt=opt)

    result, polys = polys[0], polys[1:]

    for poly in polys:
        result = result.lcm(poly)

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def lcm(f, g, *gens, **args):
    """
    Compute LCM of ``f`` and ``g``.

    Examples
    ========

    >>> lcm(x**2 - 1, x**2 - 3*x + 2)
    x**3 - 2*x**2 - x + 2

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        domain, (a, b) = construct_domain(exc.exprs)
        return domain.to_expr(domain.lcm(a, b))

    result = F.lcm(G)

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def terms_gcd(f, *gens, **args):
    """
    Remove GCD of terms from ``f``.

    If the ``deep`` flag is True, then the arguments of ``f`` will have
    terms_gcd applied to them.

    If a fraction is factored out of ``f`` and ``f`` is an Add, then
    an unevaluated Mul will be returned so that automatic simplification
    does not redistribute it. The hint ``clear``, when set to False, can be
    used to prevent such factoring when all coefficients are not fractions.

    Examples
    ========

    >>> terms_gcd(x**6*y**2 + x**3*y, x, y)
    x**3*y*(x**3*y + 1)

    The default action of polys routines is to expand the expression
    given to them. terms_gcd follows this behavior:

    >>> terms_gcd((3+3*x)*(x+x*y))
    3*x*(x*y + x + y + 1)

    If this is not desired then the hint ``expand`` can be set to False.
    In this case the expression will be treated as though it were comprised
    of one or more terms:

    >>> terms_gcd((3+3*x)*(x+x*y), expand=False)
    (3*x + 3)*(x*y + x)

    In order to traverse factors of a Mul or the arguments of other
    functions, the ``deep`` hint can be used:

    >>> terms_gcd((3 + 3*x)*(x + x*y), expand=False, deep=True)
    3*x*(x + 1)*(y + 1)
    >>> terms_gcd(cos(x + x*y), deep=True)
    cos(x*(y + 1))

    Rationals are factored out by default:

    >>> terms_gcd(x + y/2)
    (2*x + y)/2

    Only the y-term had a coefficient that was a fraction; if one
    does not want to factor out the 1/2 in cases like this, the
    flag ``clear`` can be set to False:

    >>> terms_gcd(x + y/2, clear=False)
    x + y/2
    >>> terms_gcd(x*y/2 + y**2, clear=False)
    y*(x/2 + y)

    The ``clear`` flag is ignored if all coefficients are fractions:

    >>> terms_gcd(x/3 + y/2, clear=False)
    (2*x + 3*y)/6

    See Also
    ========
    diofant.core.exprtools.gcd_terms, diofant.core.exprtools.factor_terms

    """
    from ..core import Equality

    orig = sympify(f)
    if not isinstance(f, Expr) or f.is_Atom:
        return orig

    if args.get('deep', False):
        new = f.func(*[terms_gcd(a, *gens, **args) for a in f.args])
        args.pop('deep')
        args['expand'] = False
        return terms_gcd(new, *gens, **args)

    if isinstance(f, Equality):
        return f

    clear = args.pop('clear', True)
    options.allowed_flags(args, ['polys'])

    F, opt = poly_from_expr(f, *gens, **args)

    J, f = F.terms_gcd()

    if opt.domain.is_Field:
        denom, f = f.clear_denoms(convert=True)

    coeff, f = f.primitive()

    if opt.domain.is_Field:
        coeff /= denom

    term = Mul(*[x**j for x, j in zip(f.gens, J)])
    if coeff == 1:
        coeff = Integer(1)
        if term == 1:
            return orig

    if clear:
        return _keep_coeff(coeff, term*f.as_expr())
    # base the clearing on the form of the original expression, not
    # the (perhaps) Mul that we have now
    coeff, f = _keep_coeff(coeff, f.as_expr(), clear=False).as_coeff_Mul()
    return _keep_coeff(coeff, term*f, clear=False)


def trunc(f, p, *gens, **args):
    """
    Reduce ``f`` modulo a constant ``p``.

    Examples
    ========

    >>> trunc(2*x**3 + 3*x**2 + 5*x + 7, 3)
    -x**3 - x + 1

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('trunc', 1, exc)

    result = F.trunc(sympify(p))

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def monic(f, *gens, **args):
    """
    Divide all coefficients of ``f`` by ``LC(f)``.

    Examples
    ========

    >>> monic(3*x**2 + 4*x + 2)
    x**2 + 4*x/3 + 2/3

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('monic', 1, exc)

    result = F.monic(auto=opt.auto)

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def content(f, *gens, **args):
    """
    Compute GCD of coefficients of ``f``.

    Examples
    ========

    >>> content(6*x**2 + 8*x + 12)
    2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('content', 1, exc)

    return F.content()


def primitive(f, *gens, **args):
    """
    Compute content and the primitive form of ``f``.

    Examples
    ========

    >>> primitive(6*x**2 + 8*x + 12)
    (2, 3*x**2 + 4*x + 6)

    >>> eq = (2 + 2*x)*x + 2

    Expansion is performed by default:

    >>> primitive(eq)
    (2, x**2 + x + 1)

    Set ``expand`` to False to shut this off. Note that the
    extraction will not be recursive; use the as_content_primitive method
    for recursive, non-destructive Rational extraction.

    >>> primitive(eq, expand=False)
    (1, x*(2*x + 2) + 2)

    >>> eq.as_content_primitive()
    (2, x*(x + 1) + 1)

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('primitive', 1, exc)

    cont, result = F.primitive()
    if not opt.polys:
        return cont, result.as_expr()
    else:
        return cont, result


def compose(f, g, *gens, **args):
    """
    Compute functional composition ``f(g)``.

    Examples
    ========

    >>> compose(x**2 + x, x - 1)
    x**2 - x

    """
    options.allowed_flags(args, ['polys'])

    try:
        (F, G), opt = parallel_poly_from_expr((f, g), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('compose', 2, exc)

    result = F.compose(G)

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def decompose(f, *gens, **args):
    """
    Compute functional decomposition of ``f``.

    Examples
    ========

    >>> decompose(x**4 + 2*x**3 - x - 1)
    [x**2 - x - 1, x**2 + x]

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('decompose', 1, exc)

    result = F.decompose()

    if not opt.polys:
        return [r.as_expr() for r in result]
    else:
        return result


def sturm(f, *gens, **args):
    """
    Compute Sturm sequence of ``f``.

    Examples
    ========

    >>> sturm(x**3 - 2*x**2 + x - 3)
    [x**3 - 2*x**2 + x - 3, 3*x**2 - 4*x + 1, 2*x/9 + 25/9, -2079/4]

    """
    options.allowed_flags(args, ['auto', 'polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('sturm', 1, exc)

    result = F.sturm(auto=opt.auto)

    if not opt.polys:
        return [r.as_expr() for r in result]
    else:
        return result


def sqf_norm(f, *gens, **args):
    """
    Compute square-free norm of ``f``.

    Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and
    ``r(x) = Norm(g(x))`` is a square-free polynomial over ``K``,
    where ``a`` is the algebraic extension of the ground domain.

    Examples
    ========

    >>> sqf_norm(x**2 + 1, extension=[sqrt(3)])
    (1, x**2 - 2*sqrt(3)*x + 4, x**4 - 4*x**2 + 16)

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('sqf_norm', 1, exc)

    s, g, r = F.sqf_norm()

    if not opt.polys:
        return Integer(s), g.as_expr(), r.as_expr()
    else:
        return Integer(s), g, r


def sqf_part(f, *gens, **args):
    """
    Compute square-free part of ``f``.

    Examples
    ========

    >>> sqf_part(x**3 - 3*x - 2)
    x**2 - x - 2

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('sqf_part', 1, exc)

    result = F.sqf_part()

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def _sorted_factors(factors, method):
    """Sort a list of ``(expr, exp)`` pairs."""
    if method == 'sqf':
        def key(obj):
            poly, exp = obj
            rep = poly.rep.to_dense()
            return exp, len(rep), len(poly.gens), default_sort_key(rep)
    else:
        def key(obj):
            poly, exp = obj
            rep = poly.rep.to_dense()
            return len(rep), len(poly.gens), exp, default_sort_key(rep)

    return sorted(factors, key=key)


def _factors_product(factors):
    """Multiply a list of ``(expr, exp)`` pairs."""
    return Mul(*[f.as_expr()**k for f, k in factors])


def _symbolic_factor_list(expr, opt, method):
    """Helper function for :func:`_symbolic_factor`."""
    coeff, factors = Integer(1), []

    args = [i._eval_factor() if hasattr(i, '_eval_factor') else i
            for i in Mul.make_args(expr)]
    for arg in args:
        if arg.is_Number:
            coeff *= arg
            continue
        if arg.is_Mul:
            args.extend(arg.args)
            continue
        if arg.is_Pow and arg.base is not E:
            base, exp = arg.args
            if base.is_Number:
                factors.append((base, exp))
                continue
        else:
            base, exp = arg, Integer(1)

        try:
            poly, _ = _poly_from_expr(base, opt)
        except PolificationFailed as exc:
            factors.append((exc.expr, exp))
        else:
            func = getattr(poly, method + '_list')

            _coeff, _factors = func()
            if _coeff != 1:
                if exp.is_Integer:
                    coeff *= _coeff**exp
                elif _coeff.is_positive:
                    factors.append((_coeff, exp))
                else:
                    _factors.append((_coeff, Integer(1)))

            if exp == 1:
                factors.extend(_factors)
            elif exp.is_integer:
                factors.extend([(f, k*exp) for f, k in _factors])
            else:
                other = []

                for f, k in _factors:
                    if f.as_expr().is_positive:
                        factors.append((f, k*exp))
                    else:
                        other.append((f, k))

                factors.append((_factors_product(other), exp))

    return coeff, factors


def _symbolic_factor(expr, opt, method):
    """Helper function for :func:`_factor`."""
    if isinstance(expr, Expr) and not expr.is_Relational:
        if hasattr(expr, '_eval_factor'):
            return expr._eval_factor()
        coeff, factors = _symbolic_factor_list(together(expr), opt, method)
        return _keep_coeff(coeff, _factors_product(factors))
    elif hasattr(expr, 'args'):
        return expr.func(*[_symbolic_factor(arg, opt, method) for arg in expr.args])
    elif hasattr(expr, '__iter__'):
        return expr.__class__([_symbolic_factor(arg, opt, method) for arg in expr])
    else:
        raise NotImplementedError


def _generic_factor_list(expr, gens, args, method):
    """Helper function for :func:`sqf_list` and :func:`factor_list`."""
    options.allowed_flags(args, ['frac', 'polys'])
    opt = options.build_options(gens, args)

    expr = sympify(expr)

    if isinstance(expr, Expr) and not expr.is_Relational:
        numer, denom = together(expr).as_numer_denom()

        cp, fp = _symbolic_factor_list(numer, opt, method)
        cq, fq = _symbolic_factor_list(denom, opt, method)

        if fq and not opt.frac:
            raise PolynomialError("a polynomial expected, got %s" % expr)

        _opt = opt.clone({'expand': True})
        if not _opt.get('gens'):
            _opt['gens'] = set().union(*[set(f.gens)
                                         for f, _ in fp + fq if f.is_Poly])

        for factors in (fp, fq):
            for i, (f, k) in enumerate(factors):
                if not f.is_Poly:
                    f, _ = _poly_from_expr(f, _opt)
                    factors[i] = (f, k)

        fp = _sorted_factors(fp, method)
        fq = _sorted_factors(fq, method)

        if not opt.polys:
            fp = [(f.as_expr(), k) for f, k in fp]
            fq = [(f.as_expr(), k) for f, k in fq]

        coeff = cp/cq

        if not opt.frac:
            return coeff, fp
        else:
            return coeff, fp, fq
    else:
        raise PolynomialError("a polynomial expected, got %s" % expr)


def _generic_factor(expr, gens, args, method):
    """Helper function for :func:`sqf` and :func:`factor`."""
    options.allowed_flags(args, [])
    opt = options.build_options(gens, args)
    return _symbolic_factor(sympify(expr), opt, method)


def to_rational_coeffs(f):
    """
    try to transform a polynomial to have rational coefficients

    try to find a transformation ``x = alpha*y``

    ``f(x) = lc*alpha**n * g(y)`` where ``g`` is a polynomial with
    rational coefficients, ``lc`` the leading coefficient.

    If this fails, try ``x = y + beta``
    ``f(x) = g(y)``

    Returns ``None`` if ``g`` not found;
    ``(lc, alpha, None, g)`` in case of rescaling
    ``(None, None, beta, g)`` in case of translation

    Notes
    =====

    Currently it transforms only polynomials without roots larger than 2.

    Examples
    ========

    >>> p = Poly(((x**2-1)*(x-2)).subs({x: x*(1 + sqrt(2))}), x, domain=EX)
    >>> lc, r, _, g = to_rational_coeffs(p)
    >>> lc, r
    (7 + 5*sqrt(2), -2*sqrt(2) + 2)
    >>> g
    Poly(x**3 + x**2 - 1/4*x - 1/4, x, domain='QQ')
    >>> r1 = simplify(1/r)
    >>> Poly(lc*r**3*(g.as_expr()).subs({x: x*r1}), x, domain=EX) == p
    True

    """
    from ..simplify import simplify

    def _try_rescale(f, f1=None):
        """
        try rescaling ``x -> alpha*x`` to convert f to a polynomial
        with rational coefficients.
        Returns ``alpha, f``; if the rescaling is successful,
        ``alpha`` is the rescaling factor, and ``f`` is the rescaled
        polynomial; else ``alpha`` is ``None``.

        """
        from ..core import Add
        if f.is_multivariate or not (f.gens[0]).is_Atom:
            return
        n = f.degree()
        lc = f.LC()
        f1 = f1 or f1.monic()
        coeffs = f1.all_coeffs()[1:]
        coeffs = [simplify(coeffx) for coeffx in coeffs]
        if coeffs[-2]:
            rescale1_x = simplify(coeffs[-2]/coeffs[-1])
            coeffs1 = []
            for i in range(len(coeffs)):
                coeffx = simplify(coeffs[i]*rescale1_x**(i + 1))
                if not coeffx.is_rational:
                    break
                coeffs1.append(coeffx)
            else:
                rescale_x = simplify(1/rescale1_x)
                x = f.gens[0]
                v = [x**n]
                for i in range(1, n + 1):
                    v.append(coeffs1[i - 1]*x**(n - i))
                f = Add(*v)
                f = Poly(f)
                return lc, rescale_x, f

    def _try_translate(f, f1=None):
        """
        try translating ``x -> x + alpha`` to convert f to a polynomial
        with rational coefficients.
        Returns ``alpha, f``; if the translating is successful,
        ``alpha`` is the translating factor, and ``f`` is the shifted
        polynomial; else ``alpha`` is ``None``.

        """
        from ..core import Add
        if f.is_multivariate or not (f.gens[0]).is_Atom:
            return
        n = f.degree()
        f1 = f1 or f1.monic()
        coeffs = f1.all_coeffs()[1:]
        c = simplify(coeffs[0])
        if c and not c.is_rational:
            func = Add
            if c.is_Add:
                args = c.args
                func = c.func
            else:
                args = [c]
            sifted = sift(args, lambda z: z.is_rational)
            c2 = sifted[False]
            alpha = -func(*c2)/n
            f2 = f1.shift(alpha)
            return alpha, f2

    def _has_square_roots(p):
        """Return True if ``f`` is a sum with square roots but no other root."""
        from ..core.exprtools import Factors
        coeffs = p.coeffs()
        has_sq = False
        for y in coeffs:
            for x in Add.make_args(y):
                f = Factors(x).factors
                r = [wx.denominator for b, wx in f.items() if
                     b.is_number and wx.is_Rational and wx.denominator >= 2]
                if not r:
                    continue
                if min(r) == 2:
                    has_sq = True
                if max(r) > 2:
                    return False
        return has_sq

    if f.domain.is_SymbolicDomain and _has_square_roots(f):
        f1 = f.monic()
        r = _try_rescale(f, f1)
        if r:
            return r[0], r[1], None, r[2]
        else:
            r = _try_translate(f, f1)
            if r:
                return None, None, r[0], r[1]


def _torational_factor_list(p, x):
    """
    helper function to factor polynomial using to_rational_coeffs

    Examples
    ========

    >>> p = expand(((x**2-1)*(x-2)).subs({x: x*(1 + sqrt(2))}))
    >>> factors = _torational_factor_list(p, x); factors
    (-2, [(-x*(1 + sqrt(2))/2 + 1, 1), (-x*(1 + sqrt(2)) - 1, 1), (-x*(1 + sqrt(2)) + 1, 1)])
    >>> expand(factors[0]*Mul(*[z[0] for z in factors[1]])) == p
    True
    >>> p = expand(((x**2-1)*(x-2)).subs({x: x + sqrt(2)}))
    >>> factors = _torational_factor_list(p, x); factors
    (1, [(x - 2 + sqrt(2), 1), (x - 1 + sqrt(2), 1), (x + 1 + sqrt(2), 1)])
    >>> expand(factors[0]*Mul(*[z[0] for z in factors[1]])) == p
    True

    """
    from ..simplify import simplify
    p1 = Poly(p, x, domain='EX')
    n = p1.degree()
    res = to_rational_coeffs(p1)
    if not res:
        return
    lc, r, t, g = res
    factors = factor_list(g.as_expr())
    if lc:
        c = simplify(factors[0]*lc*r**n)
        r1 = simplify(1/r)
        a = []
        for z in factors[1:][0]:
            a.append((simplify(z[0].subs({x: x*r1})), z[1]))
    else:
        c = factors[0]
        a = []
        for z in factors[1:][0]:
            a.append((z[0].subs({x: x - t}), z[1]))
    return c, a


def sqf_list(f, *gens, **args):
    """
    Compute a list of square-free factors of ``f``.

    Examples
    ========

    >>> sqf_list(2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16)
    (2, [(x + 1, 2), (x + 2, 3)])

    """
    return _generic_factor_list(f, gens, args, method='sqf')


def sqf(f, *gens, **args):
    """
    Compute square-free factorization of ``f``.

    Examples
    ========

    >>> sqf(2*x**5 + 16*x**4 + 50*x**3 + 76*x**2 + 56*x + 16)
    2*(x + 1)**2*(x + 2)**3

    """
    return _generic_factor(f, gens, args, method='sqf')


def factor_list(f, *gens, **args):
    """
    Compute a list of irreducible factors of ``f``.

    Examples
    ========

    >>> factor_list(2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y)
    (2, [(x + y, 1), (x**2 + 1, 2)])

    """
    return _generic_factor_list(f, gens, args, method='factor')


def factor(f, *gens, **args):
    """
    Compute the factorization of expression, ``f``, into irreducibles. (To
    factor an integer into primes, use ``factorint``.)

    There two modes implemented: symbolic and formal. If ``f`` is not an
    instance of :class:`Poly` and generators are not specified, then the
    former mode is used. Otherwise, the formal mode is used.

    In symbolic mode, :func:`factor` will traverse the expression tree and
    factor its components without any prior expansion, unless an instance
    of :class:`~diofant.core.add.Add` is encountered (in this case formal factorization is
    used). This way :func:`factor` can handle large or symbolic exponents.

    By default, the factorization is computed over the rationals. To factor
    over other domain, e.g. an algebraic or finite field, use appropriate
    options: ``extension``, ``modulus`` or ``domain``.

    Examples
    ========

    >>> factor(2*x**5 + 2*x**4*y + 4*x**3 + 4*x**2*y + 2*x + 2*y)
    2*(x + y)*(x**2 + 1)**2

    >>> factor(x**2 + 1)
    x**2 + 1
    >>> factor(x**2 + 1, modulus=2)
    (x + 1)**2
    >>> factor(x**2 + 1, gaussian=True)
    (x - I)*(x + I)

    >>> factor(x**2 - 2, extension=sqrt(2))
    (x - sqrt(2))*(x + sqrt(2))

    >>> factor((x**2 - 1)/(x**2 + 4*x + 4))
    (x - 1)*(x + 1)/(x + 2)**2
    >>> factor((x**2 + 4*x + 4)**10000000*(x**2 + 1))
    (x + 2)**20000000*(x**2 + 1)

    By default, factor deals with an expression as a whole:

    >>> eq = 2**(x**2 + 2*x + 1)
    >>> factor(eq)
    2**(x**2 + 2*x + 1)

    If the ``deep`` flag is True then subexpressions will
    be factored:

    >>> factor(eq, deep=True)
    2**((x + 1)**2)

    See Also
    ========
    diofant.ntheory.factor_.factorint

    """
    f = sympify(f)
    if args.pop('deep', False):
        partials = {}
        muladd = f.atoms(Mul, Add)
        for p in muladd:
            fac = factor(p, *gens, **args)
            if (fac.is_Mul or fac.is_Pow) and fac != p:
                partials[p] = fac
        return f.xreplace(partials)

    try:
        return _generic_factor(f, gens, args, method='factor')
    except PolynomialError as msg:
        if not f.is_commutative:
            from ..core.exprtools import factor_nc
            return factor_nc(f)
        else:
            raise PolynomialError(msg)


def intervals(F, all=False, eps=None, inf=None, sup=None, strict=False, sqf=False):
    """
    Compute isolating intervals for roots of ``f``.

    Examples
    ========

    >>> intervals(x**2 - 3)
    [((-2, -1), 1), ((1, 2), 1)]
    >>> intervals(x**2 - 3, eps=1e-2)
    [((-26/15, -19/11), 1), ((19/11, 26/15), 1)]

    """
    if not hasattr(F, '__iter__'):
        try:
            F = Poly(F)
        except GeneratorsNeeded:
            return []

        return F.intervals(all=all, eps=eps, inf=inf, sup=sup, sqf=sqf)
    else:
        polys, opt = parallel_poly_from_expr(F, domain='QQ')

        if len(opt.gens) > 1:
            raise MultivariatePolynomialError

        R = polys[0].rep.ring
        polys = [p.rep for p in polys]

        if eps is not None:
            eps = opt.domain.convert(eps)

            if eps <= 0:
                raise ValueError("'eps' must be a positive rational")

        if inf is not None:
            inf = opt.domain.convert(inf)
        if sup is not None:
            sup = opt.domain.convert(sup)

        intervals = R.dup_isolate_real_roots_list(polys, eps=eps, inf=inf, sup=sup,
                                                  strict=strict)

        result = []

        for (s, t), indices in intervals:
            s, t = opt.domain.to_expr(s), opt.domain.to_expr(t)
            result.append(((s, t), indices))

        return result


def refine_root(f, s, t, eps=None, steps=None, check_sqf=False):
    """
    Refine an isolating interval of a root to the given precision.

    Examples
    ========

    >>> refine_root(x**2 - 3, 1, 2, eps=1e-2)
    (19/11, 26/15)

    """
    try:
        F = Poly(f)
    except GeneratorsNeeded:
        raise PolynomialError(
            "can't refine a root of %s, not a polynomial" % f)

    return F.refine_root(s, t, eps=eps, steps=steps, check_sqf=check_sqf)


def count_roots(f, inf=None, sup=None):
    """
    Return the number of roots of ``f`` in ``[inf, sup]`` interval.

    If one of ``inf`` or ``sup`` is complex, it will return the number of roots
    in the complex rectangle with corners at ``inf`` and ``sup``.

    Examples
    ========

    >>> count_roots(x**4 - 4, -3, 3)
    2
    >>> count_roots(x**4 - 4, 0, 1 + 3*I)
    1

    """
    try:
        F = Poly(f, greedy=False)
    except GeneratorsNeeded:
        raise PolynomialError("can't count roots of %s, not a polynomial" % f)

    return F.count_roots(inf=inf, sup=sup)


def real_roots(f, multiple=True):
    """
    Return a list of real roots with multiplicities of ``f``.

    Examples
    ========

    >>> real_roots(2*x**3 - 7*x**2 + 4*x + 4)
    [-1/2, 2, 2]

    """
    try:
        F = Poly(f, greedy=False)
    except GeneratorsNeeded:
        raise PolynomialError("can't compute real roots of %s, "
                              "not a polynomial" % f)

    return F.real_roots(multiple=multiple)


def nroots(f, n=15, maxsteps=50, cleanup=True):
    """
    Compute numerical approximations of roots of ``f``.

    Examples
    ========

    >>> nroots(x**2 - 3, n=15)
    [-1.73205080756888, 1.73205080756888]
    >>> nroots(x**2 - 3, n=30)
    [-1.73205080756887729352744634151, 1.73205080756887729352744634151]

    """
    try:
        F = Poly(f, greedy=False)
    except GeneratorsNeeded:
        raise PolynomialError(
            "can't compute numerical roots of %s, not a polynomial" % f)

    return F.nroots(n=n, maxsteps=maxsteps, cleanup=cleanup)


def ground_roots(f, *gens, **args):
    """
    Compute roots of ``f`` by factorization in the ground domain.

    Examples
    ========

    >>> ground_roots(x**6 - 4*x**4 + 4*x**3 - x**2)
    {0: 2, 1: 2}

    """
    options.allowed_flags(args, [])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('ground_roots', 1, exc)

    return F.ground_roots()


def nth_power_roots_poly(f, n, *gens, **args):
    """
    Construct a polynomial with n-th powers of roots of ``f``.

    Examples
    ========

    >>> f = x**4 - x**2 + 1
    >>> g = factor(nth_power_roots_poly(f, 2))

    >>> g
    (x**2 - x + 1)**2

    >>> R_f = [(r**2).expand() for r in roots(f)]
    >>> R_g = roots(g)

    >>> set(R_f) == set(R_g)
    True

    """
    options.allowed_flags(args, ['polys'])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('nth_power_roots_poly', 1, exc)

    result = F.nth_power_roots_poly(n)

    if not opt.polys:
        return result.as_expr()
    else:
        return result


def cancel(f, *gens, **args):
    """
    Cancel common factors in a rational function ``f``.

    Examples
    ========

    >>> A = Symbol('A', commutative=False)

    >>> cancel((2*x**2 - 2)/(x**2 - 2*x + 1))
    (2*x + 2)/(x - 1)
    >>> cancel((sqrt(3) + sqrt(15)*A)/(sqrt(2) + sqrt(10)*A))
    sqrt(6)/2

    """
    from ..core.exprtools import factor_terms
    from ..functions import Piecewise
    options.allowed_flags(args, ['polys'])

    f = sympify(f)

    if not isinstance(f, Tuple):
        if f.is_Atom or isinstance(f, Relational) or not isinstance(f, Expr):
            return f
        f = factor_terms(f, radical=True)
        p, q = f.as_numer_denom()

    elif len(f) == 2:
        p, q = f
    else:
        return factor_terms(f)

    try:
        (F, G), opt = parallel_poly_from_expr((p, q), *gens, **args)
    except PolificationFailed:
        if not isinstance(f, (tuple, Tuple)):
            return f
        else:
            return Integer(1), p, q
    except PolynomialError:
        assert not f.is_commutative or f.has(Piecewise)
        # Handling of noncommutative and/or piecewise expressions
        if f.is_Add or f.is_Mul:
            sifted = sift(f.args, lambda x: x.is_commutative is True and not x.has(Piecewise))
            c, nc = sifted[True], sifted[False]
            nc = [cancel(i) for i in nc]
            return f.func(cancel(f.func._from_args(c)), *nc)
        else:
            reps = []
            pot = preorder_traversal(f)
            next(pot)
            for e in pot:
                # XXX: This should really skip anything that's not Expr.
                if isinstance(e, (tuple, Tuple, BooleanAtom)):
                    continue
                reps.append((e, cancel(e)))
                pot.skip()  # this was handled successfully
            return f.xreplace(dict(reps))

    c, P, Q = F.cancel(G)

    if not isinstance(f, (tuple, Tuple)):
        return c*(P.as_expr()/Q.as_expr())
    else:
        if not opt.polys:
            return c, P.as_expr(), Q.as_expr()
        else:
            return c, P, Q


def reduced(f, G, *gens, **args):
    """
    Reduces a polynomial ``f`` modulo a set of polynomials ``G``.

    Given a polynomial ``f`` and a set of polynomials ``G = (g_1, ..., g_n)``,
    computes a set of quotients ``q = (q_1, ..., q_n)`` and the remainder ``r``
    such that ``f = q_1*g_1 + ... + q_n*g_n + r``, where ``r`` vanishes or ``r``
    is a completely reduced polynomial with respect to ``G``.

    Examples
    ========

    >>> reduced(2*x**4 + y**2 - x**2 + y**3, [x**3 - x, y**3 - y])
    ([2*x, 1], x**2 + y**2 + y)

    """
    options.allowed_flags(args, ['polys', 'auto'])

    try:
        polys, opt = parallel_poly_from_expr([f] + list(G), *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('reduced', 0, exc)

    domain = opt.domain
    retract = False

    if opt.auto and domain.is_Ring and not domain.is_Field:
        opt = opt.clone({'domain': domain.field})
        retract = True

    _ring = opt.domain.poly_ring(*opt.gens, order=opt.order)

    for i, poly in enumerate(polys):
        poly = dict(poly.set_domain(opt.domain).rep)
        polys[i] = _ring.from_dict(poly)

    Q, r = polys[0].div(polys[1:])

    Q = [Poly._from_dict(dict(q), opt) for q in Q]
    r = Poly._from_dict(dict(r), opt)

    if retract:
        try:
            _Q, _r = [q.to_ring() for q in Q], r.to_ring()
        except CoercionFailed:
            pass
        else:
            Q, r = _Q, _r

    if not opt.polys:
        return [q.as_expr() for q in Q], r.as_expr()
    else:
        return Q, r


def groebner(F, *gens, **args):
    r"""
    Computes the reduced Gröbner basis for a set of polynomials.

    Parameters
    ==========

    F : list
        a set of polynomials
    \*gens : tuple
        polynomial generators
    \**args : dict
        a dictionary of parameters, namely

        order : str, optional
            Monomial order, defaults to ``lex``.
        method : {'buchberger', 'f5b'}, optional
            Set algorithm to compute Gröbner basis.  By default, an improved
            implementation of the Buchberger algorithm is used.
        field : bool, optional
            Force coefficients domain to be a field.  Defaults to False.

    Examples
    ========

    >>> F = [x*y - 2*x, 2*x**2 - y**2]

    >>> groebner(F)
    GroebnerBasis([2*x**2 - y**2, x*y - 2*x, y**3 - 2*y**2],
                  x, y, domain='ZZ', order='lex')

    >>> groebner(F, order=grevlex)
    GroebnerBasis([y**3 - 2*y**2, 2*x**2 - y**2, x*y - 2*x],
                  x, y, domain='ZZ', order='grevlex')

    >>> groebner(F, field=True)
    GroebnerBasis([x**2 - y**2/2, x*y - 2*x, y**3 - 2*y**2],
                  x, y, domain='QQ', order='lex')

    References
    ==========

    * :cite:`Buchberger2001systems`
    * :cite:`Cox2015ideals`

    See Also
    ========

    diofant.solvers.polysys.solve_poly_system

    """
    return GroebnerBasis(F, *gens, **args)


class GroebnerBasis(Basic):
    """Represents a reduced Gröbner basis."""

    def __new__(cls, F, *gens, **args):
        """Compute a reduced Gröbner basis for a system of polynomials."""
        options.allowed_flags(args, ['polys', 'method'])

        try:
            polys, opt = parallel_poly_from_expr(F, *gens, **args)
        except PolificationFailed as exc:
            raise ComputationFailed('groebner', len(F), exc)

        ring = opt.domain.poly_ring(*opt.gens, order=opt.order)

        if not ring.domain.is_Exact:
            raise ValueError('Domain must be exact, got %s' % ring.domain)

        polys = [ring.from_dict(dict(_.rep))
                 for _ in polys if not _.is_zero]

        G = _groebner(polys, ring, method=opt.method)
        G = [Poly._from_dict(g, opt) for g in G]

        return cls._new(G, opt)

    @classmethod
    def _new(cls, basis, options):
        obj = Basic.__new__(cls)

        obj._basis = tuple(basis)
        obj._options = options

        return obj

    @property
    def args(self):
        return Tuple(*self.exprs), Tuple(*self.gens)

    @property
    def exprs(self):
        return [poly.as_expr() for poly in self._basis]

    @property
    def polys(self):
        return list(self._basis)

    @property
    def gens(self):
        return self._options.gens

    @property
    def domain(self):
        return self._options.domain

    @property
    def order(self):
        return self._options.order

    def __len__(self):
        return len(self._basis)

    def __iter__(self):
        if self._options.polys:
            return iter(self.polys)
        else:
            return iter(self.exprs)

    def __getitem__(self, item):
        if self._options.polys:
            basis = self.polys
        else:
            basis = self.exprs

        return basis[item]

    def __hash__(self):
        return hash((self._basis, tuple(sorted(self._options.items()))))

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._basis == other._basis and self._options == other._options
        elif iterable(other):
            return self.polys == list(other) or self.exprs == list(other)
        else:
            return False

    @property
    def dimension(self):
        """Dimension of the ideal, generated by a Gröbner basis."""
        sets = self.independent_sets
        if sets is not None:
            return max(len(s) for s in sets)

    @property
    def independent_sets(self):
        """Compute independent sets for ideal, generated by a Gröbner basis.

        References
        ==========

        * :cite:`Kredel1988indep`

        """
        if self.contains(Integer(1)):
            return

        HTG = [_.LM(order=self.order) for _ in self.polys]

        def dimrec(S, U, M):
            U1 = U.copy()
            while U1:
                x = U1.pop(0)
                S1 = S + [x]
                t = Monomial(Mul(*S1), self.gens)
                for ht in HTG:
                    if all(a and b or not a for a, b in zip(ht, t)):
                        break
                else:
                    M = dimrec(S1, U1, M)
            if any(all(_ in m for _ in S) for m in M):
                return M
            else:
                return [S] + M

        return dimrec([], list(self.gens), [])

    def set_order(self, order):
        """
        Convert a Gröbner basis from one ordering to another.

        Notes
        =====

        The FGLM algorithm :cite:`Faugere1993groebner` used to convert reduced Gröbner bases
        of zero-dimensional ideals from one ordering to another.  Sometimes it
        is infeasible to compute a Gröbner basis with respect to a particular
        ordering directly.

        Examples
        ========

        >>> F = [x**2 - 3*y - x + 1, y**2 - 2*x + y - 1]
        >>> G = groebner(F, x, y, order='grlex')

        >>> G.set_order('lex') == groebner(F, x, y, order='lex')
        True

        """
        src_order = self.order
        dst_order = monomial_key(order)

        if src_order == dst_order:
            return self

        if self.dimension != 0:
            raise NotImplementedError("can't convert Gröbner bases of "
                                      "ideals with positive dimension")

        polys = self.polys
        domain = self.domain

        opt = self._options.clone({'domain': domain.field,
                                   'order': dst_order})

        _ring = opt.domain.poly_ring(*opt.gens, order=src_order)

        for i, poly in enumerate(polys):
            poly = dict(poly.set_domain(opt.domain).rep)
            polys[i] = _ring.from_dict(poly)

        G = matrix_fglm(polys, _ring, dst_order)
        G = [Poly._from_dict(dict(g), opt) for g in G]

        if not domain.is_Field:
            G = [g.clear_denoms(convert=True)[1] for g in G]
            opt.domain = domain

        return self._new(G, opt)

    def reduce(self, expr, auto=True):
        """
        Reduces a polynomial modulo a Gröbner basis.

        Given a polynomial ``f`` and a set of polynomials ``G = (g_1, ..., g_n)``,
        computes a set of quotients ``q = (q_1, ..., q_n)`` and the remainder ``r``
        such that ``f = q_1*f_1 + ... + q_n*f_n + r``, where ``r`` vanishes or ``r``
        is a completely reduced polynomial with respect to ``G``.

        Examples
        ========

        >>> f = 2*x**4 - x**2 + y**3 + y**2
        >>> G = groebner([x**3 - x, y**3 - y])

        >>> G.reduce(f)
        ([2*x, 1], x**2 + y**2 + y)
        >>> Q, r = _

        >>> expand(sum(q*g for q, g in zip(Q, G)) + r)
        2*x**4 - x**2 + y**3 + y**2
        >>> _ == f
        True

        """
        poly = Poly._from_expr(expr, self._options)
        polys = [poly] + list(self._basis)

        opt = self._options
        domain = self.domain

        retract = False

        if auto and domain.is_Ring and not domain.is_Field:
            opt = self._options.clone({'domain': domain.field})
            retract = True

        _ring = opt.domain.poly_ring(*opt.gens, order=opt.order)

        for i, poly in enumerate(polys):
            poly = dict(poly.set_domain(opt.domain).rep)
            polys[i] = _ring.from_dict(poly)

        Q, r = polys[0].div(polys[1:])

        Q = [Poly._from_dict(dict(q), opt) for q in Q]
        r = Poly._from_dict(dict(r), opt)

        if retract:
            try:
                _Q, _r = [q.to_ring() for q in Q], r.to_ring()
            except CoercionFailed:
                pass
            else:
                Q, r = _Q, _r

        if not opt.polys:
            return [q.as_expr() for q in Q], r.as_expr()
        else:
            return Q, r

    def contains(self, poly):
        """
        Check if ``poly`` belongs the ideal generated by ``self``.

        Examples
        ========

        >>> f = 2*x**3 + y**3 + 3*y
        >>> G = groebner([x**2 + y**2 - 1, x*y - 2])

        >>> G.contains(f)
        True
        >>> G.contains(f + 1)
        False

        """
        return self.reduce(poly)[1] == 0


def poly(expr, *gens, **args):
    """
    Efficiently transform an expression into a polynomial.

    Examples
    ========

    >>> poly(x*(x**2 + x - 1)**2)
    Poly(x**5 + 2*x**4 - x**3 - 2*x**2 + x, x, domain='ZZ')

    """
    options.allowed_flags(args, [])

    def _poly(expr, opt):
        terms, poly_terms = [], []

        for term in Add.make_args(expr):
            factors, poly_factors = [], []

            for factor in Mul.make_args(term):
                if factor.is_Add:
                    poly_factors.append(_poly(factor, opt))
                elif (factor.is_Pow and factor.base.is_Add and
                      factor.exp.is_Integer and factor.exp >= 0):
                    poly_factors.append(_poly(factor.base,
                                              opt)**factor.exp)
                else:
                    factors.append(factor)

            if not poly_factors:
                terms.append(term)
            else:
                product = poly_factors[0]

                for factor in poly_factors[1:]:
                    product *= factor

                if factors:
                    factor = Mul(*factors)

                    if factor.is_Number:
                        product *= factor
                    else:
                        product *= Poly._from_expr(factor, opt)

                poly_terms.append(product)

        if not poly_terms:
            result = Poly._from_expr(expr, opt)
        else:
            result = poly_terms[0]

            for term in poly_terms[1:]:
                result += term

            if terms:
                term = Add(*terms)

                if term.is_Number:
                    result += term
                else:
                    result += Poly._from_expr(term, opt)

        return result.reorder(*opt.get('gens', ()), **args)

    expr = sympify(expr)

    if expr.is_Poly:
        return Poly(expr, *gens, **args)

    if 'expand' not in args:
        args['expand'] = False

    opt = options.build_options(gens, args)

    return _poly(expr, opt)
