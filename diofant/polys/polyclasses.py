"""OO layer for several polynomial representations. """

from ..core import oo
from ..core.sympify import CantSympify
from .densearith import (dmp_abs, dmp_add, dmp_div, dmp_exquo,
                         dmp_exquo_ground, dmp_l1_norm, dmp_max_norm, dmp_mul,
                         dmp_mul_ground, dmp_neg, dmp_pdiv, dmp_pexquo,
                         dmp_pow, dmp_pquo, dmp_prem, dmp_quo, dmp_quo_ground,
                         dmp_rem, dmp_sub)
from .densebasic import (dmp_convert, dmp_deflate, dmp_degree_in,
                         dmp_degree_list, dmp_eject, dmp_exclude,
                         dmp_from_dict, dmp_ground, dmp_ground_LC,
                         dmp_ground_nth, dmp_ground_p, dmp_ground_TC,
                         dmp_inject, dmp_list_terms, dmp_one_p, dmp_permute,
                         dmp_slice_in, dmp_terms_gcd, dmp_to_dict,
                         dmp_to_tuple, dmp_validate, dmp_zero_p)
from .densetools import (dmp_clear_denoms, dmp_compose, dmp_diff_in,
                         dmp_eval_in, dmp_ground_content, dmp_ground_monic,
                         dmp_ground_primitive, dmp_ground_trunc,
                         dmp_integrate_in, dmp_lift, dup_decompose, dup_shift)
from .euclidtools import (dmp_cancel, dmp_discriminant, dmp_gcd, dmp_inner_gcd,
                          dmp_lcm, dmp_resultant, dmp_subresultants, dup_gcdex,
                          dup_half_gcdex, dup_invert)
from .factortools import dmp_factor_list, dmp_irreducible_p, dup_cyclotomic_p
from .polyerrors import PolynomialError, UnificationFailed
from .rootisolation import (dup_count_complex_roots, dup_count_real_roots,
                            dup_isolate_all_roots, dup_isolate_all_roots_sqf,
                            dup_isolate_real_roots, dup_isolate_real_roots_sqf,
                            dup_refine_real_root, dup_sturm)
from .sqfreetools import (dmp_sqf_list, dmp_sqf_list_include, dmp_sqf_norm,
                          dmp_sqf_p, dmp_sqf_part)


class DMP(CantSympify):
    """Dense Multivariate Polynomials over `K`. """

    def __init__(self, rep, dom, lev=None):
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

    def __hash__(self):
        return hash((self.__class__.__name__, self.to_tuple(),
                     self.lev, self.domain))

    def unify(self, other):
        """Unify representations of two multivariate polynomials. """
        if not isinstance(other, DMP) or self.lev != other.lev:
            raise UnificationFailed("can't unify %s with %s" % (self, other))

        if self.domain == other.domain:
            return self.lev, self.domain, self.per, self.rep, other.rep
        else:
            lev, dom = self.lev, self.domain.unify(other.domain)

            F = dmp_convert(self.rep, lev, self.domain, dom)
            G = dmp_convert(other.rep, lev, other.domain, dom)

            def per(rep, dom=dom, lev=lev, kill=False):
                if kill:
                    if not lev:
                        return rep
                    else:
                        lev -= 1
                return DMP(rep, dom, lev)

            return lev, dom, per, F, G

    def per(self, rep, dom=None, kill=False):
        """Create a DMP out of the given representation. """
        lev = self.lev

        if kill:
            if not lev:
                return rep
            else:
                lev -= 1

        if dom is None:
            dom = self.domain

        return DMP(rep, dom, lev)

    @classmethod
    def zero(cls, lev, dom):
        return DMP(0, dom, lev)

    @classmethod
    def one(cls, lev, dom):
        return DMP(1, dom, lev)

    @classmethod
    def from_list(cls, rep, lev, dom):
        """Create an instance of ``cls`` given a list of native coefficients. """
        return cls(dmp_convert(rep, lev, None, dom), dom, lev)

    def to_dict(self, zero=False):
        """Convert ``self`` to a dict representation with native coefficients. """
        return dmp_to_dict(self.rep, self.lev, self.domain, zero=zero)

    def to_diofant_dict(self, zero=False):
        """Convert ``self`` to a dict representation with Diofant coefficients. """
        rep = dmp_to_dict(self.rep, self.lev, self.domain, zero=zero)

        for k, v in rep.items():
            rep[k] = self.domain.to_expr(v)

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
    def from_monoms_coeffs(cls, monoms, coeffs, lev, dom):
        return DMP(dict(zip(monoms, coeffs)), dom, lev)

    def to_ring(self):
        """Make the ground domain a ring. """
        return self.convert(self.domain.ring)

    def to_field(self):
        """Make the ground domain a field. """
        return self.convert(self.domain.field)

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
            n = dmp_degree_in(self.rep, 0, 0)

            if n < 0:
                return [(0,)]
            else:
                return [(n - i,) for i, c in enumerate(self.rep)]
        else:
            raise PolynomialError('multivariate polynomials not supported')

    def all_terms(self):
        """Returns all terms from a ``self``. """
        if not self.lev:
            n = dmp_degree_in(self.rep, 0, 0)

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

        >>> DMP([[[ZZ(1)]], [[ZZ(1)], [ZZ(2)]]], ZZ).exclude()
        ([2], DMP([[1], [1, 2]], ZZ))
        """
        J, F, u = dmp_exclude(self.rep, self.lev, self.domain)
        return J, self.__class__(F, self.domain, u)

    def permute(self, P):
        r"""
        Returns a polynomial in `K[x_{P(1)}, ..., x_{P(n)}]`.

        Examples
        ========

        >>> DMP([[[ZZ(2)], [ZZ(1), ZZ(0)]], [[]]], ZZ).permute([1, 0, 2])
        DMP([[[2], []], [[1, 0], []]], ZZ)

        >>> DMP([[[ZZ(2)], [ZZ(1), ZZ(0)]], [[]]], ZZ).permute([1, 2, 0])
        DMP([[[1], []], [[2, 0], []]], ZZ)
        """
        return self.per(dmp_permute(self.rep, P, self.lev, self.domain))

    def terms_gcd(self):
        """Remove GCD of terms from the polynomial ``self``. """
        J, F = dmp_terms_gcd(self.rep, self.lev, self.domain)
        return J, self.per(F)

    def quo_ground(self, c):
        """Quotient of ``self`` by a an element of the ground domain. """
        return self.per(dmp_quo_ground(self.rep, self.domain.convert(c),
                                       self.lev, self.domain))

    def exquo_ground(self, c):
        """Exact quotient of ``self`` by a an element of the ground domain. """
        return self.per(dmp_exquo_ground(self.rep, self.domain.convert(c),
                                         self.lev, self.domain))

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

    def quo(self, other):
        """Computes polynomial quotient of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_quo(F, G, lev, dom))

    def exquo(self, other):
        """Computes polynomial exact quotient of ``self`` and ``other``. """
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_exquo(F, G, lev, dom))

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
        return DMP(result, self.domain, self.lev + int(new_symbol))

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

    def sqf_norm(self):
        """Computes square-free norm of ``self``. """
        s, g, r = dmp_sqf_norm(self.rep, self.lev, self.domain)
        return s, self.per(g), self.per(r, dom=self.domain.domain)

    def sqf_part(self):
        """Computes square-free part of ``self``. """
        return self.per(dmp_sqf_part(self.rep, self.lev, self.domain))

    def sqf_list(self):
        """Returns a list of square-free factors of ``self``. """
        coeff, factors = dmp_sqf_list(self.rep, self.lev, self.domain)
        return coeff, [(self.per(g), k) for g, k in factors]

    def sqf_list_include(self):
        """Returns a list of square-free factors of ``self``. """
        factors = dmp_sqf_list_include(self.rep, self.lev, self.domain)
        return [(self.per(g), k) for g, k in factors]

    def factor_list(self):
        """Returns a list of irreducible factors of ``self``. """
        coeff, factors = dmp_factor_list(self.rep, self.lev, self.domain)
        return coeff, [(self.per(g), k) for g, k in factors]

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
    def is_squarefree(self):
        """Returns ``True`` if ``self`` is a square-free polynomial. """
        return dmp_sqf_p(self.rep, self.lev, self.domain)

    @property
    def is_monic(self):
        """Returns ``True`` if the leading coefficient of ``self`` is one. """
        return dmp_ground_LC(self.rep, self.lev, self.domain) == self.domain.one

    @property
    def is_primitive(self):
        """Returns ``True`` if the GCD of the coefficients of ``self`` is one. """
        return dmp_ground_content(self.rep, self.lev, self.domain) == self.domain.one

    @property
    def is_linear(self):
        """Returns ``True`` if ``self`` is linear in all its variables. """
        return all(sum(monom) <= 1 for monom in dmp_to_dict(self.rep, self.lev,
                                                            self.domain))

    @property
    def is_quadratic(self):
        """Returns ``True`` if ``self`` is quadratic in all its variables. """
        return all(sum(monom) <= 2 for monom in dmp_to_dict(self.rep, self.lev,
                                                            self.domain))

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
        return self.per(dmp_abs(self.rep, self.lev, self.domain))

    def __neg__(self):
        return self.per(dmp_neg(self.rep, self.lev, self.domain))

    def __add__(self, other):
        if not isinstance(other, DMP):
            other = self.per(dmp_ground(self.domain.convert(other), self.lev))

        lev, dom, per, F, G = self.unify(other)
        return per(dmp_add(F, G, lev, dom))

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if not isinstance(other, DMP):
            other = self.per(dmp_ground(self.domain.convert(other), self.lev))

        lev, dom, per, F, G = self.unify(other)
        return per(dmp_sub(F, G, lev, dom))

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, DMP):
            lev, dom, per, F, G = self.unify(other)
            return per(dmp_mul(F, G, lev, dom))
        else:
            return self.per(dmp_mul_ground(self.rep, self.domain.convert(other),
                                           self.lev, self.domain))

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, n):
        if isinstance(n, int):
            return self.per(dmp_pow(self.rep, n, self.lev, self.domain))
        else:
            raise TypeError("``int`` expected, got %s" % type(n))

    def __divmod__(self, other):
        lev, dom, per, F, G = self.unify(other)
        q, r = dmp_div(F, G, lev, dom)
        return per(q), per(r)

    def __mod__(self, other):
        lev, dom, per, F, G = self.unify(other)
        return per(dmp_rem(F, G, lev, dom))

    def __floordiv__(self, other):
        if isinstance(other, DMP):
            return self.quo(other)
        else:
            return self.quo_ground(other)

    def __eq__(self, other):
        try:
            _, _, _, F, G = self.unify(other)

            return F == G
        except UnificationFailed:
            pass

        return False

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

    def __bool__(self):
        return not dmp_zero_p(self.rep, self.lev)
