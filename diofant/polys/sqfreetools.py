"""Square-free decomposition algorithms and related tools. """

from .densearith import dmp_mul, dmp_mul_ground, dmp_neg, dmp_quo, dmp_sub
from .densebasic import (dmp_convert, dmp_ground, dmp_ground_LC, dmp_ground_p,
                         dmp_inject, dmp_one_p, dmp_raise, dmp_swap,
                         dmp_zero_p)
from .densetools import (dmp_compose, dmp_diff_in, dmp_ground_monic,
                         dmp_ground_primitive)
from .euclidtools import dmp_gcd, dmp_resultant
from .galoistools import gf_sqf_list
from .polyerrors import DomainError


def dmp_sqf_p(f, u, K):
    """
    Return ``True`` if ``f`` is a square-free polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_sqf_p(R(0))
    True
    >>> R.dmp_sqf_p((x + y)**2)
    False
    >>> R.dmp_sqf_p(x**2 + y**2)
    True
    """
    if dmp_ground_p(f, None, u):
        return True
    else:
        g = f
        for i in range(u + 1):
            g = dmp_gcd(g, dmp_diff_in(f, 1, i, u, K), u, K)
            if dmp_ground_p(g, None, u):
                return True
        else:
            return False


def dmp_sqf_norm(f, u, K):
    """
    Square-free norm of ``f`` in ``K[X]``, useful over algebraic domains.

    Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and ``r(x) = Norm(g(x))``
    is a square-free polynomial over K, where ``a`` is the algebraic extension of ``K``.

    Examples
    ========

    >>> K = QQ.algebraic_field(I)
    >>> R, x, y = ring("x y", K)

    >>> R.dmp_sqf_norm(x*y + y**2)
    (1, x*y - I*x + y**2 - 3*I*y - 2,
     x**2*y**2 + x**2 + 2*x*y**3 + 2*x*y + y**4 + 5*y**2 + 4)
    """
    if not K.is_AlgebraicField:
        raise DomainError("ground domain must be algebraic")

    g = dmp_raise(K.mod.rep, u + 1, 0, K.domain)
    F = dmp_raise([K.one, -K.unit], u, 0, K)

    s = 0

    while True:
        h, _ = dmp_inject(f, u, K, front=True)
        r = dmp_resultant(g, h, u + 1, K.domain)

        if dmp_sqf_p(r, u, K.domain):
            return s, f, r
        else:
            for j in range(u + 1):
                f = dmp_swap(f, 0, j, u, K)
                f = dmp_compose(f, F, u, K)
                f = dmp_swap(f, 0, j, u, K)
            s += 1


def dmp_sqf_part(f, u, K):
    """
    Returns square-free part of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_sqf_part(x**3 + 2*x**2*y + x*y**2)
    x**2 + x*y

    """
    if K.is_FiniteField:
        _, sqf = dmp_sqf_list(f, u, K)

        g = [K.one]
        for f, _ in sqf:
            g = dmp_mul(g, f, u, K)

        return g

    if dmp_zero_p(f, u):
        return f

    if K.is_negative(dmp_ground_LC(f, u, K)):
        f = dmp_neg(f, u, K)

    gcd = f
    for i in range(u + 1):
        gcd = dmp_gcd(gcd, dmp_diff_in(f, 1, i, u, K), u, K)
    sqf = dmp_quo(f, gcd, u, K)

    if K.is_Field:
        return dmp_ground_monic(sqf, u, K)
    else:
        return dmp_ground_primitive(sqf, u, K)[1]


def dmp_gf_sqf_list(f, u, K):
    """Compute square-free decomposition of ``f`` in ``GF(p)[X]``. """
    if not u:
        f = dmp_convert(f, u, K, K.domain)

        coeff, factors = gf_sqf_list(f, K.mod, K.domain)

        for i, (f, k) in enumerate(factors):
            factors[i] = (dmp_convert(f, u, K.domain, K), k)

        return K.convert(coeff, K.domain), factors

    else:
        raise NotImplementedError('multivariate polynomials over finite fields')


def dmp_rr_yun0_sqf_list(f, u, K):
    """Compute square-free decomposition of ``f`` in zero-characteristic ring ``K``.

    References
    ==========

    * [LeeM13]_, page 8
    """
    if dmp_ground_p(f, None, u):
        return []

    result, count = [], 1
    qs = [dmp_diff_in(f, 1, i, u, K) for i in range(u + 1)]

    g = f
    for q in qs:
        g = dmp_gcd(g, q, u, K)

    while not dmp_one_p(f, u, K):
        for i in range(u + 1):
            qs[i] = dmp_quo(qs[i], g, u, K)
        f = dmp_quo(f, g, u, K)
        for i in range(u + 1):
            qs[i] = dmp_sub(qs[i], dmp_diff_in(f, 1, i, u, K), u, K)

        g = f
        for q in qs:
            g = dmp_gcd(g, q, u, K)
        if not dmp_one_p(g, u, K):
            result.append((g, count))

        count += 1

    return result


def dmp_sqf_list(f, u, K):
    """
    Return square-free decomposition of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_sqf_list(x**5 + 2*x**4*y + x**3*y**2)
    (1, [(x + y, 2), (x, 3)])
    """
    if K.is_FiniteField:
        return dmp_gf_sqf_list(f, u, K)

    if K.is_Field:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        coeff, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
            coeff = -coeff

    return coeff, dmp_rr_yun0_sqf_list(f, u, K)


def dmp_sqf_list_include(f, u, K):
    """
    Return square-free decomposition of a polynomial in ``K[x]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_sqf_list_include(x**5 + 2*x**4*y + x**3*y**2)
    [(1, 1), (x + y, 2), (x, 3)]
    """
    coeff, factors = dmp_sqf_list(f, u, K)

    if factors and factors[0][1] == 1:
        g = dmp_mul_ground(factors[0][0], coeff, u, K)
        return [(g, 1)] + factors[1:]
    else:
        g = dmp_ground(coeff, u)
        return [(g, 1)] + factors
