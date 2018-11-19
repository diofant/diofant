"""Square-free decomposition algorithms and related tools. """

from .densearith import (dmp_mul, dmp_mul_ground, dmp_neg, dmp_pow, dmp_quo,
                         dmp_sub, dup_mul)
from .densebasic import (dmp_convert, dmp_degree, dmp_degree_in, dmp_ground_LC,
                         dmp_inject, dmp_one_p, dmp_raise, dmp_zero_p)
from .densetools import (dmp_compose, dmp_diff_in, dmp_ground_monic,
                         dmp_ground_primitive, dup_shift)
from .euclidtools import dmp_gcd, dmp_resultant
from .galoistools import gf_sqf_list, gf_sqf_part
from .polyerrors import DomainError


def dmp_sqf_p_in(f, N, u, K):
    """
    Return ``True`` if ``f`` is a square-free polynomial in ``K[X]``.

    Examples
    ========

    >>> dmp_sqf_p_in([[]], [0], 1, ZZ)
    True
    >>> dmp_sqf_p_in([[1], [2, 0], [1, 0, 0]], [0], 1, ZZ)
    False
    >>> dmp_sqf_p_in([[1], [], [1, 0, 0]], [0], 1, ZZ)
    True
    """
    if dmp_zero_p(f, u):
        return True
    else:
        g = f
        for i in N:
            g = dmp_gcd(g, dmp_diff_in(f, 1, i, u, K), u, K)
            if max(dmp_degree_in(g, j, u) for j in N) <= 0:
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
    >>> _, X, Y = ring("x y", QQ)

    >>> s, f, r = R.dmp_sqf_norm(x*y + y**2)

    >>> s == 1
    True
    >>> f == x*y + y**2 - I*y
    True
    >>> r == X**2*Y**2 + 2*X*Y**3 + Y**4 + Y**2
    True

    """
    if not K.is_AlgebraicField:
        raise DomainError("ground domain must be algebraic")

    g = dmp_raise(K.mod.rep, u + 1, 0, K.domain)
    F = dmp_raise([K.one, -K.unit], u, 0, K)

    s = 0

    while True:
        h, _ = dmp_inject(f, u, K, front=True)
        r = dmp_resultant(g, h, u + 1, K.domain)

        if dmp_sqf_p_in(r, [0], u, K.domain):
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


def dmp_sqf_part_in(f, N, u, K):
    """
    Returns square-free part of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_sqf_part_in(x**3 + 2*x**2*y + x*y**2, [0, 1])
    x**2 + x*y

    """
    if K.is_FiniteField:
        return dmp_gf_sqf_part(f, u, K)

    if dmp_zero_p(f, u):
        return f

    if K.is_negative(dmp_ground_LC(f, u, K)):
        f = dmp_neg(f, u, K)

    gcd = f
    for i in N:
        gcd = dmp_gcd(gcd, dmp_diff_in(f, 1, i, u, K), u, K)
    sqf = dmp_quo(f, gcd, u, K)

    if K.has_Field:
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

        return [K.convert(coeff, K.domain)], factors

    else:
        raise NotImplementedError('multivariate polynomials over finite fields')


def dmp_sqf_list_in(f, N, u, K):
    """
    Return square-free decomposition of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**5 + 2*x**4*y + x**3*y**2

    >>> R.dmp_sqf_list_in(f, [0, 1])
    (1, [(x + y, 2), (x, 3)])
    """
    if K.is_FiniteField:
        return dmp_gf_sqf_list(f, u, K)

    if K.has_Field:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        coeff, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)
            coeff = -coeff

    coeff = dmp_mul_ground(f, coeff, u, K)

    result, count = [], 1
    qs = [dmp_diff_in(f, 1, i, u, K) for i in N]

    g = f
    for q in qs:
        g = dmp_gcd(g, q, u, K)

    while max(dmp_degree_in(f, j, u) for j in N) > 0:
        for i in range(len(N)):
            qs[i] = dmp_quo(qs[i], g, u, K)
        f = dmp_quo(f, g, u, K)
        for i, j in enumerate(N):
            qs[i] = dmp_sub(qs[i], dmp_diff_in(f, 1, j, u, K), u, K)

        g = f
        for q in qs:
            g = dmp_gcd(g, q, u, K)
        if not dmp_one_p(g, u, K):
            result.append((g, count))
            coeff = dmp_quo(coeff, dmp_pow(g, count, u, K), u, K)

        count += 1

    return coeff, result


def dmp_sqf_list_include_in(f, N, u, K):
    """
    Return square-free decomposition of a polynomial in ``K[x]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**5 + 2*x**4*y + x**3*y**2

    >>> R.dmp_sqf_list_include_in(f, [0, 1])
    [(1, 1), (x + y, 2), (x, 3)]
    """
    coeff, factors = dmp_sqf_list_in(f, N, u, K)

    if factors and factors[0][1] == 1:
        g = dmp_mul(factors[0][0], coeff, u, K)
        return [(g, 1)] + factors[1:]
    else:
        return [(coeff, 1)] + factors


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

    f = dmp_ground_monic(f, 0, K)

    if not dmp_degree(f, 0):
        return []
    else:
        g = dmp_gcd(f, dup_shift(f, K.one, K), 0, K)
        H = dup_gff_list(g, K)

        for i, (h, k) in enumerate(H):
            g = dup_mul(g, dup_shift(h, -K(k), K), K)
            H[i] = (h, k + 1)

        f = dmp_quo(f, g, 0, K)

        if not dmp_degree(f, 0):
            return H
        else:
            return [(f, 1)] + H
