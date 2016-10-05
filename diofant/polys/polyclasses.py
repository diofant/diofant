"""OO layer for several polynomial representations. """

from diofant.core.sympify import CantSympify
from diofant.polys.polyerrors import CoercionFailed, NotReversible
from diofant import oo


class GenericPoly:
    """Base class for low-level polynomial representations. """

    def ground_to_ring(self):
        """Make the ground domain a ring. """
        return self.set_domain(self.domain.get_ring())

    def ground_to_field(self):
        """Make the ground domain a field. """
        return self.set_domain(self.domain.get_field())

    def ground_to_exact(self):
        """Make the ground domain exact. """
        return self.set_domain(self.domain.get_exact())

    @classmethod
    def _perify_factors(cls, result, include):
        if include:
            coeff, factors = result
        else:
            coeff = result

        factors = [(cls(g), k) for g, k in factors]

        if include:
            return coeff, factors
        else:
            return factors

from diofant.polys.densebasic import (
    dmp_validate,
    dup_normal, dmp_normal,
    dup_convert, dmp_convert,
    dmp_from_diofant,
    dup_strip,
    dup_degree, dmp_degree_in,
    dmp_degree_list,
    dmp_negative_p,
    dup_LC, dmp_ground_LC,
    dup_TC, dmp_ground_TC,
    dmp_ground_nth,
    dmp_one, dmp_ground,
    dmp_zero_p, dmp_one_p, dmp_ground_p,
    dup_from_dict, dmp_from_dict,
    dmp_to_dict,
    dmp_deflate,
    dmp_inject, dmp_eject,
    dmp_terms_gcd,
    dmp_list_terms, dmp_exclude,
    dmp_slice_in, dmp_permute,
    dmp_to_tuple,)

from diofant.polys.densearith import (
    dmp_add_ground,
    dmp_sub_ground,
    dmp_mul_ground,
    dmp_quo_ground,
    dmp_exquo_ground,
    dmp_abs,
    dup_neg, dmp_neg,
    dup_add, dmp_add,
    dup_sub, dmp_sub,
    dup_mul, dmp_mul,
    dmp_sqr,
    dup_pow, dmp_pow,
    dmp_pdiv,
    dmp_prem,
    dmp_pquo,
    dmp_pexquo,
    dmp_div,
    dup_rem, dmp_rem,
    dmp_quo,
    dmp_exquo,
    dmp_add_mul, dmp_sub_mul,
    dmp_max_norm,
    dmp_l1_norm)

from diofant.polys.densetools import (
    dmp_clear_denoms,
    dmp_integrate_in,
    dmp_diff_in,
    dmp_eval_in,
    dup_revert,
    dmp_ground_trunc,
    dmp_ground_content,
    dmp_ground_primitive,
    dmp_ground_monic,
    dmp_compose,
    dup_decompose,
    dup_shift,
    dmp_lift)

from diofant.polys.euclidtools import (
    dup_half_gcdex, dup_gcdex, dup_invert,
    dmp_subresultants,
    dmp_resultant,
    dmp_discriminant,
    dmp_inner_gcd,
    dmp_gcd,
    dmp_lcm,
    dmp_cancel)

from diofant.polys.sqfreetools import (
    dup_gff_list,
    dmp_sqf_p,
    dmp_sqf_norm,
    dmp_sqf_part,
    dmp_sqf_list, dmp_sqf_list_include)

from diofant.polys.factortools import (
    dup_cyclotomic_p, dmp_irreducible_p,
    dmp_factor_list, dmp_factor_list_include)

from diofant.polys.rootisolation import (
    dup_isolate_real_roots_sqf,
    dup_isolate_real_roots,
    dup_isolate_all_roots_sqf,
    dup_isolate_all_roots,
    dup_refine_real_root,
    dup_count_real_roots,
    dup_count_complex_roots,
    dup_sturm)

from diofant.polys.polyerrors import (
    UnificationFailed,
    PolynomialError)


def init_normal_DMP(rep, lev, dom):
    return DMP(dmp_normal(rep, lev, dom), dom, lev)


class DMP(CantSympify):
    """Dense Multivariate Polynomials over `K`. """

    def __init__(self, rep, dom, lev=None, ring=None):
        if lev is not None:
            if type(rep) is dict:
                rep = dmp_from_dict(rep, lev, dom)
            elif type(rep) is not list:
                rep = dmp_ground(dom.convert(rep), lev)
        else:
            rep, lev = dmp_validate(rep)

        self.rep = rep
        self.lev = lev
        self.domain = dom
        self.ring = ring

    def __repr__(self):
        return "%s(%s, %s, %s)" % (self.__class__.__name__, self.rep,
                                   self.domain, self.ring)

    def __hash__(self):
        return hash((self.__class__.__name__, self.to_tuple(),
                     self.lev, self.domain, self.ring))

    def unify(self, other):
        """Unify representations of two multivariate polynomials. """
        if not isinstance(other, DMP) or self.lev != other.lev:
            raise UnificationFailed("can't unify %s with %s" % (self, other))

        if self.domain == other.domain and self.ring == other.ring:
            return self.lev, self.domain, self.per, self.rep, other.rep
        else:
            lev, dom = self.lev, self.domain.unify(other.domain)
            ring = self.ring
            if other.ring is not None:
                if ring is not None:
                    ring = ring.unify(other.ring)
                else:
                    ring = other.ring

            F = dmp_convert(self.rep, lev, self.domain, dom)
            G = dmp_convert(other.rep, lev, other.domain, dom)

            def per(rep, dom=dom, lev=lev, kill=False):
                if kill:
                    if not lev:
                        return rep
                    else:
                        lev -= 1

                return DMP(rep, dom, lev, ring)

            return lev, dom, per, F, G

    def per(self, rep, dom=None, kill=False, ring=None):
        """Create a DMP out of the given representation. """
        lev = self.lev

        if kill:
            if not lev:
                return rep
            else:
                lev -= 1

        if dom is None:
            dom = self.domain

        if ring is None:
            ring = self.ring

        return DMP(rep, dom, lev, ring)

    @classmethod
    def zero(cls, lev, dom, ring=None):
        return DMP(0, dom, lev, ring)

    @classmethod
    def one(cls, lev, dom, ring=None):
        return DMP(1, dom, lev, ring)

    @classmethod
    def from_list(cls, rep, lev, dom):
        """Create an instance of ``cls`` given a list of native coefficients. """
        return cls(dmp_convert(rep, lev, None, dom), dom, lev)

    @classmethod
    def from_diofant_list(cls, rep, lev, dom):
        """Create an instance of ``cls`` given a list of Diofant coefficients. """
        return cls(dmp_from_diofant(rep, lev, dom), dom, lev)

    def to_dict(self, zero=False):
        """Convert ``self`` to a dict representation with native coefficients. """
        return dmp_to_dict(self.rep, self.lev, self.domain, zero=zero)

    def to_diofant_dict(self, zero=False):
        """Convert ``self`` to a dict representation with Diofant coefficients. """
        rep = dmp_to_dict(self.rep, self.lev, self.domain, zero=zero)

        for k, v in rep.items():
            rep[k] = self.domain.to_diofant(v)

        return rep

    def to_tuple(self):
        """
        Convert ``self`` to a tuple representation with native coefficients.

        This is needed for hashing.
        """
        return dmp_to_tuple(self.rep, self.lev)

    @classmethod
    def from_dict(cls, rep, lev, dom):
        """Construct and instance of ``cls`` from a ``dict`` representation. """
        return cls(dmp_from_dict(rep, lev, dom), dom, lev)

    @classmethod
    def from_monoms_coeffs(cls, monoms, coeffs, lev, dom, ring=None):
        return DMP(dict(zip(monoms, coeffs)), dom, lev, ring)

    def to_ring(self):
        """Make the ground domain a ring. """
        return self.convert(self.domain.get_ring())

    def to_field(self):
        """Make the ground domain a field. """
        return self.convert(self.domain.get_field())

    def to_exact(self):
        """Make the ground domain exact. """
        return self.convert(self.domain.get_exact())

    def convert(self, dom):
        """Convert the ground domain of ``self``. """
        if self.domain == dom:
            return self
        else:
            return DMP(dmp_convert(self.rep, self.lev, self.domain, dom),
                       dom, self.lev)

    def slice(self, m, n, j=0):
        """Take a continuous subsequence of terms of ``self``. """
        return self.per(dmp_slice_in(self.rep, m, n, j, self.lev, self.domain))

    def coeffs(self, order=None):
        """Returns all non-zero coefficients from ``self`` in lex order. """
        return [c for _, c in dmp_list_terms(self.rep, self.lev,
                                             self.domain, order=order)]

    def monoms(self, order=None):
        """Returns all non-zero monomials from ``self`` in lex order. """
        return [m for m, _ in dmp_list_terms(self.rep, self.lev,
                                             self.domain, order=order)]

    def terms(self, order=None):
        """Returns all non-zero terms from ``self`` in lex order. """
        return dmp_list_terms(self.rep, self.lev, self.domain, order=order)

    def all_coeffs(self):
        """Returns all coefficients from ``self``. """
        if not self.lev:
            if not self:
                return [self.domain.zero]
            else:
                return [c for c in self.rep]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def all_monoms(self):
        """Returns all monomials from ``self``. """
        if not self.lev:
            n = dup_degree(self.rep)

            if n < 0:
                return [(0,)]
            else:
                return [(n - i,) for i, c in enumerate(self.rep)]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def all_terms(self):
        """Returns all terms from a ``self``. """
        if not self.lev:
            n = dup_degree(self.rep)

            if n < 0:
                return [((0,), self.domain.zero)]
            else:
                return [((n - i,), c) for i, c in enumerate(self.rep)]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def lift(self):
        """Convert algebraic coefficients to rationals. """
        return self.per(dmp_lift(self.rep, self.lev, self.domain),
                        dom=self.domain.domain)

    def deflate(self):
        """Reduce degree of `self` by mapping `x_i^m` to `y_i`. """
        J, F = dmp_deflate(self.rep, self.lev, self.domain)
        return J, self.per(F)

    def inject(self, front=False):
        """Inject ground domain generators into ``self``. """
        F, lev = dmp_inject(self.rep, self.lev, self.domain, front=front)
        return self.__class__(F, self.domain.domain, lev)

    def eject(self, dom, front=False):
        """Eject selected generators into the ground domain. """
        F = dmp_eject(self.rep, self.lev, dom, front=front)
        return self.__class__(F, dom, self.lev - len(dom.symbols))

    def exclude(self):
        r"""
        Remove useless generators from ``self``.

        Returns the removed generators and the new excluded ``self``.

        Examples
        ========

        >>> from diofant.polys.domains import ZZ

        >>> DMP([[[ZZ(1)]], [[ZZ(1)], [ZZ(2)]]], ZZ).exclude()
        ([2], DMP([[1], [1, 2]], ZZ, None))
        """
        J, F, u = dmp_exclude(self.rep, self.lev, self.domain)
        return J, self.__class__(F, self.domain, u)

    def permute(self, P):
        r"""
        Returns a polynomial in `K[x_{P(1)}, ..., x_{P(n)}]`.

        Examples
        ========

        >>> from diofant.polys.domains import ZZ

        >>> DMP([[[ZZ(2)], [ZZ(1), ZZ(0)]], [[]]], ZZ).permute([1, 0, 2])
        DMP([[[2], []], [[1, 0], []]], ZZ, None)

        >>> DMP([[[ZZ(2)], [ZZ(1), ZZ(0)]], [[]]], ZZ).permute([1, 2, 0])
        DMP([[[1], []], [[2, 0], []]], ZZ, None)
        """
        return self.per(dmp_permute(self.rep, P, self.lev, self.domain))

    def terms_gcd(self):
        """Remove GCD of terms from the polynomial ``self``. """
        J, F = dmp_terms_gcd(self.rep, self.lev, self.domain)
        return J, self.per(F)

    def add_ground(self, c):
        """Add an element of the ground domain to ``self``. """
        return self.per(dmp_add_ground(self.rep, self.domain.convert(c),
                                       self.lev, self.domain))

    def sub_ground(self, c):
        """Subtract an element of the ground domain from ``self``. """
        return self.per(dmp_sub_ground(self.rep, self.domain.convert(c),
                                       self.lev, self.domain))

    def mul_ground(self, c):
        """Multiply ``self`` by a an element of the ground domain. """
        return self.per(dmp_mul_ground(self.rep, self.domain.convert(c),
                                       self.lev, self.domain))

    def quo_ground(self, c):
        """Quotient of ``self`` by a an element of the ground domain. """
        return self.per(dmp_quo_ground(self.rep, self.domain.convert(c),
                                       self.lev, self.domain))

    def exquo_ground(self, c):
        """Exact quotient of ``self`` by a an element of the ground domain. """
        return self.per(dmp_exquo_ground(self.rep, self.domain.convert(c),
                                         self.lev, self.domain))

    def abs(self):
        """Make all coefficients in ``self`` positive. """
        return self.per(dmp_abs(self.rep, self.lev, self.domain))

    def neg(self):
        """Negate all coefficients in ``self``. """
        return self.per(dmp_neg(self.rep, self.lev, self.domain))

    def add(self, other):
        """Add two multivariate polynomials ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_add(F, G, lev, dom))

    def sub(self, other):
        """Subtract two multivariate polynomials ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_sub(F, G, lev, dom))

    def mul(self, other):
        """Multiply two multivariate polynomials ``f`` and ``g``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_mul(F, G, lev, dom))

    def sqr(self):
        """Square a multivariate polynomial ``self``. """
        return self.per(dmp_sqr(self.rep, self.lev, self.domain))

    def pow(self, n):
        """Raise ``self`` to a non-negative power ``n``. """
        if isinstance(n, int):
            return self.per(dmp_pow(self.rep, n, self.lev, self.domain))
        else:
            raise TypeError("``int`` expected, got %s" % type(n))

    def pdiv(self, other):
        """Polynomial pseudo-division of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        q, r = dmp_pdiv(F, G, lev, dom)
        return per(q), per(r)

    def prem(self, other):
        """Polynomial pseudo-remainder of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_prem(F, G, lev, dom))

    def pquo(self, other):
        """Polynomial pseudo-quotient of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_pquo(F, G, lev, dom))

    def pexquo(self, other):
        """Polynomial exact pseudo-quotient of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_pexquo(F, G, lev, dom))

    def div(self, other):
        """Polynomial division with remainder of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        q, r = dmp_div(F, G, lev, dom)
        return per(q), per(r)

    def rem(self, other):
        """Computes polynomial remainder of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_rem(F, G, lev, dom))

    def quo(self, other):
        """Computes polynomial quotient of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_quo(F, G, lev, dom))

    def exquo(self, other):
        """Computes polynomial exact quotient of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        res = per(dmp_exquo(F, G, lev, dom))
        if self.ring and res not in self.ring:
            from diofant.polys.polyerrors import ExactQuotientFailed
            raise ExactQuotientFailed(self, other, self.ring)
        return res

    def degree(self, j=0):
        """Returns the leading degree of ``self`` in ``x_j``. """
        if isinstance(j, int):
            return dmp_degree_in(self.rep, j, self.lev)
        else:
            raise TypeError("``int`` expected, got %s" % type(j))

    def degree_list(self):
        """Returns a list of degrees of ``self``. """
        return dmp_degree_list(self.rep, self.lev)

    def total_degree(self):
        """Returns the total degree of ``self``. """
        return max(sum(m) for m in self.monoms())

    def homogenize(self, s):
        """Return homogeneous polynomial of ``self``"""
        td = self.total_degree()
        result = {}
        new_symbol = (s == len(self.terms()[0][0]))
        for term in self.terms():
            d = sum(term[0])
            if d < td:
                i = td - d
            else:
                i = 0
            if new_symbol:
                result[term[0] + (i,)] = term[1]
            else:
                l = list(term[0])
                l[s] += i
                result[tuple(l)] = term[1]
        return DMP(result, self.domain, self.lev + int(new_symbol), self.ring)

    def homogeneous_order(self):
        """Returns the homogeneous order of ``self``. """
        if self.is_zero:
            return -oo

        monoms = self.monoms()
        tdeg = sum(monoms[0])

        for monom in monoms:
            _tdeg = sum(monom)

            if _tdeg != tdeg:
                return

        return tdeg

    def LC(self):
        """Returns the leading coefficient of ``self``. """
        return dmp_ground_LC(self.rep, self.lev, self.domain)

    def TC(self):
        """Returns the trailing coefficient of ``self``. """
        return dmp_ground_TC(self.rep, self.lev, self.domain)

    def nth(self, *N):
        """Returns the ``n``-th coefficient of ``self``. """
        if all(isinstance(n, int) for n in N):
            return dmp_ground_nth(self.rep, N, self.lev, self.domain)
        else:
            raise TypeError("a sequence of integers expected")

    def max_norm(self):
        """Returns maximum norm of ``self``. """
        return dmp_max_norm(self.rep, self.lev, self.domain)

    def l1_norm(self):
        """Returns l1 norm of ``self``. """
        return dmp_l1_norm(self.rep, self.lev, self.domain)

    def clear_denoms(self):
        """Clear denominators, but keep the ground domain. """
        coeff, F = dmp_clear_denoms(self.rep, self.lev, self.domain)
        return coeff, self.per(F)

    def integrate(self, m=1, j=0):
        """Computes the ``m``-th order indefinite integral of ``self`` in ``x_j``. """
        if not isinstance(m, int):
            raise TypeError("``int`` expected, got %s" % type(m))

        if not isinstance(j, int):
            raise TypeError("``int`` expected, got %s" % type(j))

        return self.per(dmp_integrate_in(self.rep, m, j, self.lev, self.domain))

    def diff(self, m=1, j=0):
        """Computes the ``m``-th order derivative of ``self`` in ``x_j``. """
        if not isinstance(m, int):
            raise TypeError("``int`` expected, got %s" % type(m))

        if not isinstance(j, int):
            raise TypeError("``int`` expected, got %s" % type(j))

        return self.per(dmp_diff_in(self.rep, m, j, self.lev, self.domain))

    def eval(self, a, j=0):
        """Evaluates ``self`` at the given point ``a`` in ``x_j``. """
        if not isinstance(j, int):
            raise TypeError("``int`` expected, got %s" % type(j))

        return self.per(dmp_eval_in(self.rep, self.domain.convert(a),
                                    j, self.lev, self.domain), kill=True)

    def half_gcdex(self, other):
        """Half extended Euclidean algorithm, if univariate. """
        lev, dom, per, F, G = self.unify(other)

        if not lev:
            s, h = dup_half_gcdex(F, G, dom)
            return per(s), per(h)
        else:
            raise ValueError('univariate polynomial expected')

    def gcdex(self, other):
        """Extended Euclidean algorithm, if univariate. """
        lev, dom, per, F, G = self.unify(other)

        if not lev:
            s, t, h = dup_gcdex(F, G, dom)
            return per(s), per(t), per(h)
        else:
            raise ValueError('univariate polynomial expected')

    def invert(self, other):
        """Invert ``self`` modulo ``other``, if possible. """
        lev, dom, per, F, G = self.unify(other)

        if not lev:
            return per(dup_invert(F, G, dom))
        else:
            raise ValueError('univariate polynomial expected')

    def revert(self, n):
        """Compute ``self**(-1)`` mod ``x**n``. """
        if not self.lev:
            return self.per(dup_revert(self.rep, n, self.domain))
        else:
            raise ValueError('univariate polynomial expected')

    def subresultants(self, other):
        """Computes subresultant PRS sequence of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        R = dmp_subresultants(F, G, lev, dom)
        return list(map(per, R))

    def resultant(self, other, includePRS=False):
        """Computes resultant of ``self`` and ``other`` via PRS. """
        lev, dom, per, F, G = self.unify(other)
        if includePRS:
            res, R = dmp_resultant(F, G, lev, dom, includePRS=includePRS)
            return per(res, kill=True), list(map(per, R))
        return per(dmp_resultant(F, G, lev, dom), kill=True)

    def discriminant(self):
        """Computes discriminant of ``self``. """
        return self.per(dmp_discriminant(self.rep, self.lev,
                                         self.domain), kill=True)

    def cofactors(self, other):
        """Returns GCD of ``self`` and ``other`` and their cofactors. """
        lev, dom, per, F, G = self.unify(other)
        h, cff, cfg = dmp_inner_gcd(F, G, lev, dom)
        return per(h), per(cff), per(cfg)

    def gcd(self, other):
        """Returns polynomial GCD of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_gcd(F, G, lev, dom))

    def lcm(self, other):
        """Returns polynomial LCM of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_lcm(F, G, lev, dom))

    def cancel(self, other, include=True):
        """Cancel common factors in a rational function ``self/other``. """
        lev, dom, per, F, G = self.unify(other)

        if include:
            F, G = dmp_cancel(F, G, lev, dom, include=True)
        else:
            cF, cG, F, G = dmp_cancel(F, G, lev, dom, include=False)

        F, G = per(F), per(G)

        if include:
            return F, G
        else:
            return cF, cG, F, G

    def trunc(self, p):
        """Reduce ``self`` modulo a constant ``p``. """
        return self.per(dmp_ground_trunc(self.rep, self.domain.convert(p),
                                         self.lev, self.domain))

    def monic(self):
        """Divides all coefficients by ``LC(self)``. """
        return self.per(dmp_ground_monic(self.rep, self.lev, self.domain))

    def content(self):
        """Returns GCD of polynomial coefficients. """
        return dmp_ground_content(self.rep, self.lev, self.domain)

    def primitive(self):
        """Returns content and a primitive form of ``self``. """
        cont, F = dmp_ground_primitive(self.rep, self.lev, self.domain)
        return cont, self.per(F)

    def compose(self, other):
        """Computes functional composition of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_compose(F, G, lev, dom))

    def decompose(self):
        """Computes functional decomposition of ``self``. """
        if not self.lev:
            return list(map(self.per, dup_decompose(self.rep, self.domain)))
        else:
            raise ValueError('univariate polynomial expected')

    def shift(self, a):
        """Efficiently compute Taylor shift ``self(x + a)``. """
        if not self.lev:
            return self.per(dup_shift(self.rep, self.domain.convert(a), self.domain))
        else:
            raise ValueError('univariate polynomial expected')

    def sturm(self):
        """Computes the Sturm sequence of ``self``. """
        if not self.lev:
            return list(map(self.per, dup_sturm(self.rep, self.domain)))
        else:
            raise ValueError('univariate polynomial expected')

    def gff_list(self):
        """Computes greatest factorial factorization of ``self``. """
        if not self.lev:
            return [(self.per(g), k) for g, k in dup_gff_list(self.rep,
                                                              self.domain)]
        else:
            raise ValueError('univariate polynomial expected')

    def sqf_norm(self):
        """Computes square-free norm of ``self``. """
        s, g, r = dmp_sqf_norm(self.rep, self.lev, self.domain)
        return s, self.per(g), self.per(r, dom=self.domain.domain)

    def sqf_part(self):
        """Computes square-free part of ``self``. """
        return self.per(dmp_sqf_part(self.rep, self.lev, self.domain))

    def sqf_list(self, all=False):
        """Returns a list of square-free factors of ``self``. """
        coeff, factors = dmp_sqf_list(self.rep, self.lev, self.domain, all)
        return coeff, [(self.per(g), k) for g, k in factors]

    def sqf_list_include(self, all=False):
        """Returns a list of square-free factors of ``self``. """
        factors = dmp_sqf_list_include(self.rep, self.lev, self.domain, all)
        return [(self.per(g), k) for g, k in factors]

    def factor_list(self):
        """Returns a list of irreducible factors of ``self``. """
        coeff, factors = dmp_factor_list(self.rep, self.lev, self.domain)
        return coeff, [(self.per(g), k) for g, k in factors]

    def factor_list_include(self):
        """Returns a list of irreducible factors of ``self``. """
        factors = dmp_factor_list_include(self.rep, self.lev, self.domain)
        return [(self.per(g), k) for g, k in factors]

    def intervals(self, all=False, eps=None, inf=None, sup=None, fast=False, sqf=False):
        """Compute isolating intervals for roots of ``self``. """
        if not self.lev:
            if not all:
                if not sqf:
                    return dup_isolate_real_roots(self.rep, self.domain, eps=eps,
                                                  inf=inf, sup=sup, fast=fast)
                else:
                    return dup_isolate_real_roots_sqf(self.rep, self.domain,
                                                      eps=eps, inf=inf,
                                                      sup=sup, fast=fast)
            else:
                if not sqf:
                    return dup_isolate_all_roots(self.rep, self.domain, eps=eps,
                                                 inf=inf, sup=sup, fast=fast)
                else:
                    return dup_isolate_all_roots_sqf(self.rep, self.domain,
                                                     eps=eps, inf=inf, sup=sup,
                                                     fast=fast)
        else:
            raise PolynomialError(
                "can't isolate roots of a multivariate polynomial")

    def refine_root(self, s, t, eps=None, steps=None, fast=False):
        """
        Refine an isolating interval to the given precision.

        ``eps`` should be a rational number.
        """
        if not self.lev:
            return dup_refine_real_root(self.rep, s, t, self.domain, eps=eps,
                                        steps=steps, fast=fast)
        else:
            raise PolynomialError(
                "can't refine a root of a multivariate polynomial")

    def count_real_roots(self, inf=None, sup=None):
        """Return the number of real roots of ``self`` in ``[inf, sup]``. """
        return dup_count_real_roots(self.rep, self.domain, inf=inf, sup=sup)

    def count_complex_roots(self, inf=None, sup=None):
        """Return the number of complex roots of ``self`` in ``[inf, sup]``. """
        return dup_count_complex_roots(self.rep, self.domain, inf=inf, sup=sup)

    @property
    def is_zero(self):
        """Returns ``True`` if ``self`` is a zero polynomial. """
        return dmp_zero_p(self.rep, self.lev)

    @property
    def is_one(self):
        """Returns ``True`` if ``self`` is a unit polynomial. """
        return dmp_one_p(self.rep, self.lev, self.domain)

    @property
    def is_ground(self):
        """Returns ``True`` if ``self`` is an element of the ground domain. """
        return dmp_ground_p(self.rep, None, self.lev)

    @property
    def is_sqf(self):
        """Returns ``True`` if ``self`` is a square-free polynomial. """
        return dmp_sqf_p(self.rep, self.lev, self.domain)

    @property
    def is_monic(self):
        """Returns ``True`` if the leading coefficient of ``self`` is one. """
        return self.domain.is_one(dmp_ground_LC(self.rep, self.lev, self.domain))

    @property
    def is_primitive(self):
        """Returns ``True`` if the GCD of the coefficients of ``self`` is one. """
        return self.domain.is_one(dmp_ground_content(self.rep, self.lev, self.domain))

    @property
    def is_linear(self):
        """Returns ``True`` if ``self`` is linear in all its variables. """
        return all(sum(monom) <= 1 for monom in dmp_to_dict(self.rep, self.lev,
                                                            self.domain).keys())

    @property
    def is_quadratic(self):
        """Returns ``True`` if ``self`` is quadratic in all its variables. """
        return all(sum(monom) <= 2 for monom in dmp_to_dict(self.rep, self.lev,
                                                            self.domain).keys())

    @property
    def is_monomial(self):
        """Returns ``True`` if ``self`` is zero or has only one term. """
        return len(self.to_dict()) <= 1

    @property
    def is_homogeneous(self):
        """Returns ``True`` if ``self`` is a homogeneous polynomial. """
        return self.homogeneous_order() is not None

    @property
    def is_irreducible(self):
        """Returns ``True`` if ``self`` has no factors over its domain. """
        return dmp_irreducible_p(self.rep, self.lev, self.domain)

    @property
    def is_cyclotomic(self):
        """Returns ``True`` if ``self`` is a cyclotomic polynomial. """
        if not self.lev:
            return dup_cyclotomic_p(self.rep, self.domain)
        else:
            return False

    def __abs__(self):
        return self.abs()

    def __neg__(self):
        return self.neg()

    def __add__(self, other):
        if not isinstance(other, DMP):
            try:
                other = self.per(dmp_ground(self.domain.convert(other), self.lev))
            except TypeError:
                return NotImplemented
            except (CoercionFailed, NotImplementedError):
                if self.ring is not None:
                    try:
                        other = self.ring.convert(other)
                    except (CoercionFailed, NotImplementedError):
                        return NotImplemented

        return self.add(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if not isinstance(other, DMP):
            try:
                other = self.per(dmp_ground(self.domain.convert(other), self.lev))
            except TypeError:
                return NotImplemented
            except (CoercionFailed, NotImplementedError):
                if self.ring is not None:
                    try:
                        other = self.ring.convert(other)
                    except (CoercionFailed, NotImplementedError):
                        return NotImplemented

        return self.sub(other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, DMP):
            return self.mul(other)
        else:
            try:
                return self.mul_ground(other)
            except TypeError:
                return NotImplemented
            except (CoercionFailed, NotImplementedError):
                if self.ring is not None:
                    try:
                        return self.mul(self.ring.convert(other))
                    except (CoercionFailed, NotImplementedError):
                        pass
                return NotImplemented

    def __div__(self, other):
        if isinstance(other, DMP):
            return self.exquo(other)
        else:
            try:
                return self.mul_ground(other)
            except TypeError:
                return NotImplemented
            except (CoercionFailed, NotImplementedError):
                if self.ring is not None:
                    try:
                        return self.exquo(self.ring.convert(other))
                    except (CoercionFailed, NotImplementedError):
                        pass
                return NotImplemented

    def __rdiv__(self, other):
        if isinstance(other, DMP):
            return other.exquo(self)
        elif self.ring is not None:
            try:
                return self.ring.convert(other).exquo(self)
            except (CoercionFailed, NotImplementedError):
                pass
        return NotImplemented

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        return self.pow(n)

    def __divmod__(self, other):
        return self.div(other)

    def __mod__(self, other):
        return self.rem(other)

    def __floordiv__(self, other):
        if isinstance(other, DMP):
            return self.quo(other)
        else:
            try:
                return self.quo_ground(other)
            except TypeError:
                return NotImplemented

    def __eq__(self, other):
        try:
            _, _, _, F, G = self.unify(other)

            if self.lev == other.lev:
                return F == G
        except UnificationFailed:
            pass

        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def eq(self, other, strict=False):
        if not strict:
            return self.__eq__(other)
        else:
            return self._strict_eq(other)

    def ne(self, other, strict=False):
        return not self.eq(other, strict=strict)

    def _strict_eq(self, other):
        return (isinstance(other, self.__class__) and self.lev == other.lev
                           and self.domain == other.domain and self.rep == other.rep)

    def __lt__(self, other):
        _, _, _, F, G = self.unify(other)
        return F.__lt__(G)

    def __le__(self, other):
        _, _, _, F, G = self.unify(other)
        return F.__le__(G)

    def __gt__(self, other):
        _, _, _, F, G = self.unify(other)
        return F.__gt__(G)

    def __ge__(self, other):
        _, _, _, F, G = self.unify(other)
        return F.__ge__(G)

    def __bool__(self):
        return not dmp_zero_p(self.rep, self.lev)


def init_normal_DMF(num, den, lev, dom):
    return DMF(dmp_normal(num, lev, dom),
               dmp_normal(den, lev, dom), dom, lev)


class DMF(CantSympify):
    """Dense Multivariate Fractions over `K`. """

    def __init__(self, rep, dom, lev=None, ring=None):
        num, den, lev = self._parse(rep, dom, lev)
        num, den = dmp_cancel(num, den, lev, dom)

        self.num = num
        self.den = den
        self.lev = lev
        self.domain = dom
        self.ring = ring

    @classmethod
    def new(cls, rep, dom, lev=None, ring=None):
        num, den, lev = cls._parse(rep, dom, lev)

        obj = object.__new__(cls)

        obj.num = num
        obj.den = den
        obj.lev = lev
        obj.domain = dom
        obj.ring = ring

        return obj

    @classmethod
    def _parse(cls, rep, dom, lev=None):
        if type(rep) is tuple:
            num, den = rep

            if lev is not None:
                if type(num) is dict:
                    num = dmp_from_dict(num, lev, dom)

                if type(den) is dict:
                    den = dmp_from_dict(den, lev, dom)
            else:
                num, num_lev = dmp_validate(num)
                den, den_lev = dmp_validate(den)

                if num_lev == den_lev:
                    lev = num_lev
                else:
                    raise ValueError('inconsistent number of levels')

            if dmp_zero_p(den, lev):
                raise ZeroDivisionError('fraction denominator')

            if dmp_zero_p(num, lev):
                den = dmp_one(lev, dom)
            else:
                if dmp_negative_p(den, lev, dom):
                    num = dmp_neg(num, lev, dom)
                    den = dmp_neg(den, lev, dom)
        else:
            num = rep

            if lev is not None:
                if type(num) is dict:
                    num = dmp_from_dict(num, lev, dom)
                elif type(num) is not list:
                    num = dmp_ground(dom.convert(num), lev)
            else:
                num, lev = dmp_validate(num)

            den = dmp_one(lev, dom)

        return num, den, lev

    def __repr__(self):
        return "%s((%s, %s), %s, %s)" % (self.__class__.__name__, self.num,
                                         self.den, self.domain, self.ring)

    def __hash__(self):
        return hash((self.__class__.__name__, dmp_to_tuple(self.num, self.lev),
                     dmp_to_tuple(self.den, self.lev), self.lev, self.domain,
                     self.ring))

    def poly_unify(self, other):
        """Unify a multivariate fraction and a polynomial. """
        if not isinstance(other, DMP) or self.lev != other.lev:
            raise UnificationFailed("can't unify %s with %s" % (self, other))

        if self.domain == other.domain and self.ring == other.ring:
            return (self.lev, self.domain, self.per,
                    (self.num, self.den), other.rep)
        else:
            lev, dom = self.lev, self.domain.unify(other.domain)
            ring = self.ring
            if other.ring is not None:
                if ring is not None:
                    ring = ring.unify(other.ring)
                else:
                    ring = other.ring

            F = (dmp_convert(self.num, lev, self.domain, dom),
                 dmp_convert(self.den, lev, self.domain, dom))

            G = dmp_convert(other.rep, lev, other.domain, dom)

            def per(num, den, cancel=True, kill=False, lev=lev):
                if kill:
                    if not lev:
                        return num/den
                    else:
                        lev = lev - 1

                if cancel:
                    num, den = dmp_cancel(num, den, lev, dom)

                return self.__class__.new((num, den), dom, lev, ring=ring)

            return lev, dom, per, F, G

    def frac_unify(self, other):
        """Unify representations of two multivariate fractions. """
        if not isinstance(other, DMF) or self.lev != other.lev:
            raise UnificationFailed("can't unify %s with %s" % (self, other))

        if self.domain == other.domain and self.ring == other.ring:
            return (self.lev, self.domain, self.per, (self.num, self.den),
                    (other.num, other.den))
        else:
            lev, dom = self.lev, self.domain.unify(other.domain)
            ring = other.ring
            if other.ring is not None:
                if ring is not None:
                    ring = ring.unify(other.ring)
                else:
                    ring = other.ring

            F = (dmp_convert(self.num, lev, self.domain, dom),
                 dmp_convert(self.den, lev, self.domain, dom))

            G = (dmp_convert(other.num, lev, other.domain, dom),
                 dmp_convert(other.den, lev, other.domain, dom))

            def per(num, den, cancel=True, kill=False, lev=lev):
                if kill:
                    if not lev:
                        return num/den
                    else:
                        lev = lev - 1

                if cancel:
                    num, den = dmp_cancel(num, den, lev, dom)

                return self.__class__.new((num, den), dom, lev, ring=ring)

            return lev, dom, per, F, G

    def per(self, num, den, cancel=True, kill=False, ring=None):
        """Create a DMF out of the given representation. """
        lev, dom = self.lev, self.domain

        if kill:
            if not lev:
                return num/den
            else:
                lev -= 1

        if cancel:
            num, den = dmp_cancel(num, den, lev, dom)

        if ring is None:
            ring = self.ring

        return self.__class__.new((num, den), dom, lev, ring=ring)

    def half_per(self, rep, kill=False):
        """Create a DMP out of the given representation. """
        lev = self.lev

        if kill:
            if not lev:
                return rep
            else:
                lev -= 1

        return DMP(rep, self.domain, lev)

    @classmethod
    def zero(cls, lev, dom, ring=None):
        return cls.new(0, dom, lev, ring=ring)

    @classmethod
    def one(cls, lev, dom, ring=None):
        return cls.new(1, dom, lev, ring=ring)

    def numer(self):
        """Returns the numerator of ``self``. """
        return self.half_per(self.num)

    def denom(self):
        """Returns the denominator of ``self``. """
        return self.half_per(self.den)

    def cancel(self):
        """Remove common factors from ``self.num`` and ``self.den``. """
        return self.per(self.num, self.den)

    def neg(self):
        """Negate all coefficients in ``self``. """
        return self.per(dmp_neg(self.num, self.lev, self.domain),
                        self.den, cancel=False)

    def add(self, other):
        """Add two multivariate fractions ``self`` and ``other``. """
        if isinstance(other, DMP):
            lev, dom, per, (F_num, F_den), G = self.poly_unify(other)
            num, den = dmp_add_mul(F_num, F_den, G, lev, dom), F_den
        else:
            lev, dom, per, F, G = self.frac_unify(other)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_add(dmp_mul(F_num, G_den, lev, dom),
                          dmp_mul(F_den, G_num, lev, dom), lev, dom)
            den = dmp_mul(F_den, G_den, lev, dom)

        return per(num, den)

    def sub(self, other):
        """Subtract two multivariate fractions ``self`` and ``other``. """
        if isinstance(other, DMP):
            lev, dom, per, (F_num, F_den), G = self.poly_unify(other)
            num, den = dmp_sub_mul(F_num, F_den, G, lev, dom), F_den
        else:
            lev, dom, per, F, G = self.frac_unify(other)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_sub(dmp_mul(F_num, G_den, lev, dom),
                          dmp_mul(F_den, G_num, lev, dom), lev, dom)
            den = dmp_mul(F_den, G_den, lev, dom)

        return per(num, den)

    def mul(self, other):
        """Multiply two multivariate fractions ``self`` and ``other``. """
        if isinstance(other, DMP):
            lev, dom, per, (F_num, F_den), G = self.poly_unify(other)
            num, den = dmp_mul(F_num, G, lev, dom), F_den
        else:
            lev, dom, per, F, G = self.frac_unify(other)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_mul(F_num, G_num, lev, dom)
            den = dmp_mul(F_den, G_den, lev, dom)

        return per(num, den)

    def pow(self, n):
        """Raise ``self`` to a non-negative power ``n``. """
        if isinstance(n, int):
            return self.per(dmp_pow(self.num, n, self.lev, self.domain),
                            dmp_pow(self.den, n, self.lev, self.domain),
                            cancel=False)
        else:
            raise TypeError("``int`` expected, got %s" % type(n))

    def quo(self, other):
        """Computes quotient of fractions ``self`` and ``other``. """
        if isinstance(other, DMP):
            lev, dom, per, (F_num, F_den), G = self.poly_unify(other)
            num, den = F_num, dmp_mul(F_den, G, lev, dom)
        else:
            lev, dom, per, F, G = self.frac_unify(other)
            (F_num, F_den), (G_num, G_den) = F, G

            num = dmp_mul(F_num, G_den, lev, dom)
            den = dmp_mul(F_den, G_num, lev, dom)

        res = per(num, den)
        if self.ring is not None and res not in self.ring:
            from diofant.polys.polyerrors import ExactQuotientFailed
            raise ExactQuotientFailed(self, other, self.ring)
        return res

    exquo = quo

    def invert(self, check=True):
        """Computes inverse of a fraction ``self``. """
        if check and self.ring is not None and not self.ring.is_unit(self):
            raise NotReversible(self, self.ring)
        return self.per(self.den, self.num, cancel=False)

    @property
    def is_zero(self):
        """Returns ``True`` if ``self`` is a zero fraction. """
        return dmp_zero_p(self.num, self.lev)

    @property
    def is_one(self):
        """Returns ``True`` if ``self`` is a unit fraction. """
        return (dmp_one_p(self.num, self.lev, self.domain) and
                dmp_one_p(self.den, self.lev, self.domain))

    def __neg__(self):
        return self.neg()

    def __add__(self, other):
        if isinstance(other, (DMP, DMF)):
            return self.add(other)

        try:
            return self.add(self.half_per(other))
        except TypeError:
            return NotImplemented
        except (CoercionFailed, NotImplementedError):
            if self.ring is not None:
                try:
                    return self.add(self.ring.convert(other))
                except (CoercionFailed, NotImplementedError):
                    pass
            return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, (DMP, DMF)):
            return self.sub(other)

        try:
            return self.sub(self.half_per(other))
        except TypeError:
            return NotImplemented
        except (CoercionFailed, NotImplementedError):
            if self.ring is not None:
                try:
                    return self.sub(self.ring.convert(other))
                except (CoercionFailed, NotImplementedError):
                    pass
            return NotImplemented

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, (DMP, DMF)):
            return self.mul(other)

        try:
            return self.mul(self.half_per(other))
        except TypeError:
            return NotImplemented
        except (CoercionFailed, NotImplementedError):
            if self.ring is not None:
                try:
                    return self.mul(self.ring.convert(other))
                except (CoercionFailed, NotImplementedError):
                    pass
            return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        return self.pow(n)

    def __div__(self, other):
        if isinstance(other, (DMP, DMF)):
            return self.quo(other)

        try:
            return self.quo(self.half_per(other))
        except TypeError:
            return NotImplemented
        except (CoercionFailed, NotImplementedError):
            if self.ring is not None:
                try:
                    return self.quo(self.ring.convert(other))
                except (CoercionFailed, NotImplementedError):
                    pass
            return NotImplemented

    def __rdiv__(self, other):
        r = self.invert(check=False)*other
        if self.ring and r not in self.ring:
            from diofant.polys.polyerrors import ExactQuotientFailed
            raise ExactQuotientFailed(other, self, self.ring)
        return r

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __eq__(self, other):
        try:
            if isinstance(other, DMP):
                _, _, _, (F_num, F_den), G = self.poly_unify(other)

                if self.lev == other.lev:
                    return dmp_one_p(F_den, self.lev, self.domain) and F_num == G
            else:
                _, _, _, F, G = self.frac_unify(other)

                if self.lev == other.lev:
                    return F == G
        except UnificationFailed:
            pass

        return False

    def __ne__(self, other):
        try:
            if isinstance(other, DMP):
                _, _, _, (F_num, F_den), G = self.poly_unify(other)

                if self.lev == other.lev:
                    return not (dmp_one_p(F_den, self.lev, self.domain) and F_num == G)
            else:
                _, _, _, F, G = self.frac_unify(other)

                if self.lev == other.lev:
                    return F != G
        except UnificationFailed:
            pass

        return True

    def __lt__(self, other):
        _, _, _, F, G = self.frac_unify(other)
        return F.__lt__(G)

    def __le__(self, other):
        _, _, _, F, G = self.frac_unify(other)
        return F.__le__(G)

    def __gt__(self, other):
        _, _, _, F, G = self.frac_unify(other)
        return F.__gt__(G)

    def __ge__(self, other):
        _, _, _, F, G = self.frac_unify(other)
        return F.__ge__(G)

    def __bool__(self):
        return not dmp_zero_p(self.num, self.lev)


def init_normal_ANP(rep, mod, dom):
    return ANP(dup_normal(rep, dom),
               dup_normal(mod, dom), dom)


class ANP(CantSympify):
    """Dense Algebraic Number Polynomials over a field. """

    def __init__(self, rep, mod, dom):
        if type(rep) is dict:
            self.rep = dup_from_dict(rep, dom)
        else:
            if type(rep) is not list:
                rep = [dom.convert(rep)]

            self.rep = dup_strip(rep)

        if isinstance(mod, DMP):
            self.mod = mod.rep
        else:
            if type(mod) is dict:
                self.mod = dup_from_dict(mod, dom)
            else:
                self.mod = dup_strip(mod)

        self.domain = dom

    def __repr__(self):
        return "%s(%s, %s, %s)" % (self.__class__.__name__, self.rep,
                                   self.mod, self.domain)

    def __hash__(self):
        return hash((self.__class__.__name__, self.to_tuple(),
                     dmp_to_tuple(self.mod, 0), self.domain))

    def unify(self, other):
        """Unify representations of two algebraic numbers. """
        if not isinstance(other, ANP) or self.mod != other.mod:
            raise UnificationFailed("can't unify %s with %s" % (self, other))

        if self.domain == other.domain:
            return self.domain, self.per, self.rep, other.rep, self.mod
        else:
            dom = self.domain.unify(other.domain)

            F = dup_convert(self.rep, self.domain, dom)
            G = dup_convert(other.rep, other.domain, dom)

            if dom != self.domain and dom != other.domain:
                mod = dup_convert(self.mod, self.domain, dom)
            else:
                if dom == self.domain:
                    mod = self.mod
                else:
                    mod = other.mod

            def per(rep):
                return ANP(rep, mod, dom)

        return dom, per, F, G, mod

    def per(self, rep, mod=None, dom=None):
        return ANP(rep, mod or self.mod, dom or self.domain)

    @classmethod
    def zero(cls, mod, dom):
        return ANP(0, mod, dom)

    @classmethod
    def one(cls, mod, dom):
        return ANP(1, mod, dom)

    def to_dict(self):
        """Convert ``self`` to a dict representation with native coefficients. """
        return dmp_to_dict(self.rep, 0, self.domain)

    def to_diofant_dict(self):
        """Convert ``self`` to a dict representation with Diofant coefficients. """
        rep = dmp_to_dict(self.rep, 0, self.domain)

        for k, v in rep.items():
            rep[k] = self.domain.to_diofant(v)

        return rep

    def to_list(self):
        """Convert ``self`` to a list representation with native coefficients. """
        return self.rep

    def to_diofant_list(self):
        """Convert ``self`` to a list representation with Diofant coefficients. """
        return [self.domain.to_diofant(c) for c in self.rep]

    def to_tuple(self):
        """
        Convert ``self`` to a tuple representation with native coefficients.

        This is needed for hashing.
        """
        return dmp_to_tuple(self.rep, 0)

    @classmethod
    def from_list(cls, rep, mod, dom):
        return ANP(dup_strip(list(map(dom.convert, rep))), mod, dom)

    def neg(self):
        return self.per(dup_neg(self.rep, self.domain))

    def add(self, other):
        dom, per, F, G, mod = self.unify(other)
        return per(dup_add(F, G, dom))

    def sub(self, other):
        dom, per, F, G, mod = self.unify(other)
        return per(dup_sub(F, G, dom))

    def mul(self, other):
        dom, per, F, G, mod = self.unify(other)
        return per(dup_rem(dup_mul(F, G, dom), mod, dom))

    def pow(self, n):
        """Raise ``self`` to a non-negative power ``n``. """
        if isinstance(n, int):
            if n < 0:
                F, n = dup_invert(self.rep, self.mod, self.domain), -n
            else:
                F = self.rep

            return self.per(dup_rem(dup_pow(F, n, self.domain),
                                    self.mod, self.domain))
        else:
            raise TypeError("``int`` expected, got %s" % type(n))

    def div(self, other):
        dom, per, F, G, mod = self.unify(other)
        return (per(dup_rem(dup_mul(F, dup_invert(G, mod, dom),
                                    dom), mod, dom)), self.zero(mod, dom))

    def rem(self, other):
        dom, _, _, _, mod = self.unify(other)
        return self.zero(mod, dom)

    def quo(self, other):
        dom, per, F, G, mod = self.unify(other)
        return per(dup_rem(dup_mul(F, dup_invert(G, mod, dom), dom), mod, dom))

    exquo = quo

    def LC(self):
        """Returns the leading coefficient of ``self``. """
        return dup_LC(self.rep, self.domain)

    def TC(self):
        """Returns the trailing coefficient of ``self``. """
        return dup_TC(self.rep, self.domain)

    @property
    def is_zero(self):
        """Returns ``True`` if ``self`` is a zero algebraic number. """
        return not self

    @property
    def is_one(self):
        """Returns ``True`` if ``self`` is a unit algebraic number. """
        return self.rep == [self.domain.one]

    @property
    def is_ground(self):
        """Returns ``True`` if ``self`` is an element of the ground domain. """
        return not self.rep or len(self.rep) == 1

    def __neg__(self):
        return self.neg()

    def __add__(self, other):
        if isinstance(other, ANP):
            return self.add(other)
        else:
            try:
                return self.add(self.per(other))
            except (CoercionFailed, TypeError):
                return NotImplemented

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, ANP):
            return self.sub(other)
        else:
            try:
                return self.sub(self.per(other))
            except (CoercionFailed, TypeError):
                return NotImplemented

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, ANP):
            return self.mul(other)
        else:
            try:
                return self.mul(self.per(other))
            except (CoercionFailed, TypeError):
                return NotImplemented

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        return self.pow(n)

    def __divmod__(self, other):
        return self.div(other)

    def __mod__(self, other):
        return self.rem(other)

    def __div__(self, other):
        if isinstance(other, ANP):
            return self.quo(other)
        else:
            try:
                return self.quo(self.per(other))
            except (CoercionFailed, TypeError):
                return NotImplemented

    __truediv__ = __div__

    def __eq__(self, other):
        try:
            _, _, F, G, _ = self.unify(other)

            return F == G
        except UnificationFailed:
            return False

    def __ne__(self, other):
        try:
            _, _, F, G, _ = self.unify(other)

            return F != G
        except UnificationFailed:
            return True

    def __lt__(self, other):
        _, _, F, G, _ = self.unify(other)
        return F.__lt__(G)

    def __le__(self, other):
        _, _, F, G, _ = self.unify(other)
        return F.__le__(G)

    def __gt__(self, other):
        _, _, F, G, _ = self.unify(other)
        return F.__gt__(G)

    def __ge__(self, other):
        _, _, F, G, _ = self.unify(other)
        return F.__ge__(G)

    def __bool__(self):
        return bool(self.rep)
