"""Square-free decomposition algorithms and related tools."""

from .densearith import dmp_mul, dmp_quo, dmp_sub
from .densebasic import (dmp_degree_in, dmp_ground_LC, dmp_ground_p,
                         dmp_inject, dmp_one_p, dmp_raise, dmp_swap,
                         dmp_zero_p)
from .densetools import (dmp_compose, dmp_diff_in, dmp_ground_monic,
                         dmp_ground_primitive)
from .euclidtools import dmp_gcd, dmp_resultant
from .polyerrors import DomainError


def dmp_sqf_p(f, u, K):
    """
    Return ``True`` if ``f`` is a square-free polynomial in ``K[X]``.

    Examples
    ========

    >>> _, x, y = ring('x y', ZZ)

    >>> ((x + y)**2).is_squarefree
    False
    >>> (x**2 + y**2).is_squarefree
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
        return False


def dmp_sqf_norm(f, u, K):
    """
    Square-free norm of ``f`` in ``K[X]``, useful over algebraic domains.

    Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and ``r(x) = Norm(g(x))``
    is a square-free polynomial over K, where ``a`` is the algebraic extension of ``K``.

    Examples
    ========

    >>> _, x, y = ring('x y', QQ.algebraic_field(I))

    >>> (x*y + y**2).sqf_norm()
    (1, x*y - I*x + y**2 - 3*I*y - 2,
     x**2*y**2 + x**2 + 2*x*y**3 + 2*x*y + y**4 + 5*y**2 + 4)

    """
    if not K.is_AlgebraicField:
        raise DomainError('ground domain must be algebraic')

    g = dmp_raise(K.mod.to_dense(), u + 1, 0, K.domain)
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

    >>> _, x, y = ring('x y', ZZ)

    >>> (x**3 + 2*x**2*y + x*y**2).sqf_part()
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

    gcd = f
    for i in range(u + 1):
        gcd = dmp_gcd(gcd, dmp_diff_in(f, 1, i, u, K), u, K)
    sqf = dmp_quo(f, gcd, u, K)

    if K.is_Field:
        return dmp_ground_monic(sqf, u, K)
    else:
        return dmp_ground_primitive(sqf, u, K)[1]


def dup_gf_musser_sqf_list(f, K):
    """Compute square-free decomposition of the monic ``f`` in ``GF(q)[x]``.

    References
    ==========

    * :cite:`Geddes1992algorithms`, algorithm 8.3

    """
    n, factors, p = 1, [], K.characteristic
    m = K.order // p

    while not dmp_ground_p(f, None, 0):
        df = dmp_diff_in(f, 1, 0, 0, K)

        if not dmp_zero_p(df, 0):
            g = dmp_gcd(f, df, 0, K)
            h = dmp_quo(f, g, 0, K)
            f = g
            i = 1

            while not dmp_one_p(h, 0, K):
                g = dmp_gcd(f, h, 0, K)
                h = dmp_quo(h, g, 0, K)

                if dmp_degree_in(h, 0, 0) > 0:
                    factors.append((h, i*n))

                f, h, i = dmp_quo(f, g, 0, K), g, i + 1

        d = dmp_degree_in(f, 0, 0) // p
        n *= p

        for i in range(d + 1):
            f[i] = f[i*p]**m
        del f[d + 1:]

    return factors


def dmp_gf_sqf_list(f, u, K):
    """Compute square-free decomposition of the monic ``f`` in ``GF(q)[X]``.

    Examples
    ========

    >>> _, x = ring('x', FF(11))
    >>> f = x**11 + 1

    Note that:

    >>> f.diff()
    0 mod 11

    This phenomenon doesn't happen in characteristic zero. However we can
    still compute square-free decomposition of ``f``:

    >>> f.sqf_list()
    (1 mod 11, [(x + 1 mod 11, 11)])

    """
    if not u:
        return dup_gf_musser_sqf_list(f, K)

    raise NotImplementedError('multivariate polynomials over finite fields')


def dmp_rr_yun0_sqf_list(f, u, K):
    """Compute square-free decomposition of ``f`` in zero-characteristic ring ``K[X]``.

    References
    ==========

    * :cite:`Geddes1992algorithms`, algorithm 8.2
    * :cite:`LeeM2013factor`, algorithm 3.1

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

    >>> _, x, y = ring('x y', ZZ)

    >>> (x**5 + 2*x**4*y + x**3*y**2).sqf_list()
    (1, [(x + y, 2), (x, 3)])

    """
    if K.is_Field:
        coeff = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)
    else:
        coeff, f = dmp_ground_primitive(f, u, K)

    if K.is_FiniteField:
        return coeff, dmp_gf_sqf_list(f, u, K)

    return coeff, dmp_rr_yun0_sqf_list(f, u, K)
