"""Square-free decomposition algorithms and related tools. """

from .densearith import dmp_mul_ground, dmp_neg, dmp_quo, dmp_sub, dup_mul
from .densebasic import (dmp_convert, dmp_degree, dmp_ground, dmp_ground_LC,
                         dmp_inject, dmp_raise, dmp_zero_p)
from .densetools import (dmp_compose, dmp_diff, dmp_ground_monic,
                         dmp_ground_primitive, dup_monic, dup_shift)
from .euclidtools import dmp_gcd, dmp_inner_gcd, dmp_resultant, dup_gcd
from .galoistools import gf_sqf_list, gf_sqf_part
from .polyerrors import DomainError


def dmp_sqf_p(f, u, K):
    """
    Return ``True`` if ``f`` is a square-free polynomial in ``K[X]``.

    Examples
    ========

    >>> dmp_sqf_p([[]], 1, ZZ)
    True
    >>> dmp_sqf_p([[1], [2, 0], [1, 0, 0]], 1, ZZ)
    False
    >>> dmp_sqf_p([[1], [], [1, 0, 0]], 1, ZZ)
    True
    """
    if dmp_zero_p(f, u):
        return True
    else:
        return not dmp_degree(dmp_gcd(f, dmp_diff(f, 1, u, K), u, K), u)


def dmp_sqf_norm(f, u, K):
    """
    Square-free norm of ``f`` in ``K[X]``, useful over algebraic domains.

    Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and ``r(x) = Norm(g(x))``
    is a square-free polynomial over K, where ``a`` is the algebraic extension of ``K``.

    Examples
    ========

    >>> K = QQ.algebraic_field(I)
    >>> R, x, y = ring("x y", K)
    >>> _, X, Y = ring("x y", QQ)

    >>> s, f, r = R.dmp_sqf_norm(x*y + y**2)

    >>> s == 1
    True
    >>> f == x*y + y**2 + K([QQ(-1), QQ(0)])*y
    True
    >>> r == X**2*Y**2 + 2*X*Y**3 + Y**4 + Y**2
    True

    """
    if not K.is_Algebraic:
        raise DomainError("ground domain must be algebraic")

    g = dmp_raise(K.mod.rep, u + 1, 0, K.domain)
    F = dmp_raise([K.one, -K.unit], u, 0, K)

    s = 0

    while True:
        h, _ = dmp_inject(f, u, K, front=True)
        r = dmp_resultant(g, h, u + 1, K.domain)

        if dmp_sqf_p(r, u, K.domain):
            break
        else:
            f, s = dmp_compose(f, F, u, K), s + 1

    return s, f, r


def dmp_gf_sqf_part(f, u, K):
    """Compute square-free part of ``f`` in ``GF(p)[X]``. """
    if not u:
        f = dmp_convert(f, u, K, K.domain)
        g = gf_sqf_part(f, K.mod, K.domain)
        return dmp_convert(g, u, K.domain, K)
    else:  # pragma: no cover
        raise NotImplementedError('multivariate polynomials over finite fields')


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
        return dmp_gf_sqf_part(f, u, K)

    if dmp_zero_p(f, u):
        return f

    if K.is_negative(dmp_ground_LC(f, u, K)):
        f = dmp_neg(f, u, K)

    gcd = dmp_gcd(f, dmp_diff(f, 1, u, K), u, K)
    sqf = dmp_quo(f, gcd, u, K)

    if K.has_Field:
        return dmp_ground_monic(sqf, u, K)
    else:
        return dmp_ground_primitive(sqf, u, K)[1]


def dmp_gf_sqf_list(f, u, K, all=False):
    """Compute square-free decomposition of ``f`` in ``GF(p)[X]``. """
    if not u:
        f = dmp_convert(f, u, K, K.domain)

        coeff, factors = gf_sqf_list(f, K.mod, K.domain, all=all)

        for i, (f, k) in enumerate(factors):
            factors[i] = (dmp_convert(f, u, K.domain, K), k)

        return K.convert(coeff, K.domain), factors

    else:
        raise NotImplementedError('multivariate polynomials over finite fields')


def dmp_sqf_list(f, u, K, all=False):
    """
    Return square-free decomposition of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**5 + 2*x**4*y + x**3*y**2

    >>> R.dmp_sqf_list(f)
    (1, [(x + y, 2), (x, 3)])
    >>> R.dmp_sqf_list(f, all=True)
    (1, [(1, 1), (x + y, 2), (x, 3)])
    """
    if K.is_FiniteField:
        return dmp_gf_sqf_list(f, u, K, all=all)

    if K.has_Field:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        coeff, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
            coeff = -coeff

    if dmp_degree(f, u) <= 0:
        return coeff, []

    result, i = [], 1

    h = dmp_diff(f, 1, u, K)
    g, p, q = dmp_inner_gcd(f, h, u, K)

    while True:
        d = dmp_diff(p, 1, u, K)
        h = dmp_sub(q, d, u, K)

        if dmp_zero_p(h, u):
            result.append((p, i))
            break

        g, p, q = dmp_inner_gcd(p, h, u, K)

        if all or dmp_degree(g, u) > 0:
            result.append((g, i))

        i += 1

    return coeff, result


def dmp_sqf_list_include(f, u, K, all=False):
    """
    Return square-free decomposition of a polynomial in ``K[x]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**5 + 2*x**4*y + x**3*y**2

    >>> R.dmp_sqf_list_include(f)
    [(1, 1), (x + y, 2), (x, 3)]
    >>> R.dmp_sqf_list_include(f, all=True)
    [(1, 1), (x + y, 2), (x, 3)]

    """
    coeff, factors = dmp_sqf_list(f, u, K, all=all)

    if factors and factors[0][1] == 1:
        g = dmp_mul_ground(factors[0][0], coeff, u, K)
        return [(g, 1)] + factors[1:]
    else:
        g = dmp_ground(coeff, u)
        return [(g, 1)] + factors


def dup_gff_list(f, K):
    """
    Compute greatest factorial factorization of ``f`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_gff_list(x**5 + 2*x**4 - x**3 - 2*x**2)
    [(x, 1), (x + 2, 4)]

    """
    if not f:
        raise ValueError("greatest factorial factorization doesn't exist for a zero polynomial")

    f = dup_monic(f, K)

    if not dmp_degree(f, 0):
        return []
    else:
        g = dup_gcd(f, dup_shift(f, K.one, K), K)
        H = dup_gff_list(g, K)

        for i, (h, k) in enumerate(H):
            g = dup_mul(g, dup_shift(h, -K(k), K), K)
            H[i] = (h, k + 1)

        f = dmp_quo(f, g, 0, K)

        if not dmp_degree(f, 0):
            return H
        else:
            return [(f, 1)] + H
