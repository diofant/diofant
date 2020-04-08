"""Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences."""

from ..core import cacheit
from ..ntheory import nextprime
from ..ntheory.modular import crt, symmetric_residue
from .densearith import (dmp_add, dmp_max_norm, dmp_mul, dmp_mul_ground,
                         dmp_mul_term, dmp_neg, dmp_pow, dmp_quo,
                         dmp_quo_ground, dmp_sub, dup_mul)
from .densebasic import (dmp_apply_pairs, dmp_convert, dmp_degree_in,
                         dmp_ground, dmp_ground_LC, dmp_LC, dmp_raise,
                         dmp_strip, dmp_zero, dmp_zero_p)
from .densetools import dmp_clear_denoms, dmp_eval_in, dmp_ground_trunc
from .polyconfig import query
from .polyerrors import HomomorphismFailed


def dup_gcdex(f, g, K):
    """Extended Euclidean algorithm in `F[x]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(1)])
    f, g = map(ring.from_dense, (f, g))
    return tuple(map(ring.to_dense, f.gcdex(g)))


def dmp_prem(f, g, u, K):
    """Polynomial pseudo-remainder in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_dense, (f, g))
    return ring.to_dense(f.prem(g))


def dup_inner_subresultants(f, g, K):
    """
    Subresultant PRS algorithm in `K[x]`.

    Computes the subresultant polynomial remainder sequence (PRS)
    and the last non-zero scalar subresultant of `f` and `g`.

    By [1] Thm. 3, these are the constants '-c' (- to optimize
    computation of sign).
    The first subdeterminant is set to 1 by convention to match
    the polynomial and the scalar subdeterminants.
    If 'deg(f) < deg(g)', the subresultants of '(g,f)' are computed.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> (x**2 + 1).resultant(x**2 - 1, includePRS=True)
    (4, [x**2 + 1, x**2 - 1, -2])

    References
    ==========

    * :cite:`Brown1978prs`

    """
    n = dmp_degree_in(f, 0, 0)
    m = dmp_degree_in(g, 0, 0)

    if n < m:
        f, g = g, f
        n, m = m, n

    if not f:
        return [], 0

    if not g:
        return [f], K.one

    R = [f, g]
    d = n - m

    b = (-K.one)**(d + 1)

    h = dmp_prem(f, g, 0, K)
    h = dmp_mul_ground(h, b, 0, K)

    lc = dmp_LC(g, K)
    c = lc**d

    c = -c

    while h:
        k = dmp_degree_in(h, 0, 0)
        R.append(h)

        f, g, m, d = g, h, k, m - k

        b = -lc * c**d

        h = dmp_prem(f, g, 0, K)
        h = dmp_quo_ground(h, b, 0, K)

        lc = dmp_LC(g, K)

        q = c**(d - 1)
        c = K.quo((-lc)**d, q)

    return R, -c


def dup_prs_resultant(f, g, K):
    """
    Resultant algorithm in `K[x]` using subresultant PRS.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> (x**2 + 1).resultant(x**2 - 1, includePRS=True)
    (4, [x**2 + 1, x**2 - 1, -2])

    """
    R, S = dmp_inner_subresultants(f, g, 0, K)

    if len(R) < 2 or dmp_degree_in(R[-1], 0, 0) > 0:
        return K.zero, R

    return S, R


def dmp_inner_subresultants(f, g, u, K):
    """
    Subresultant PRS algorithm in `K[X]`.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> f.resultant(g, includePRS=True) == (b.drop(0), [f, g, a, b])
    True

    """
    if not u:
        return dup_inner_subresultants(f, g, K)

    n = dmp_degree_in(f, 0, u)
    m = dmp_degree_in(g, 0, u)

    if n < m:
        f, g = g, f
        n, m = m, n

    if dmp_zero_p(f, u):
        return [], dmp_zero(u)

    v = u - 1
    if dmp_zero_p(g, u):
        return [f], dmp_ground(K.one, v)

    R = [f, g]
    d = n - m

    b = dmp_pow(dmp_ground(-K.one, v), d + 1, v, K)

    h = dmp_prem(f, g, u, K)
    h = dmp_mul_term(h, b, 0, u, K)

    lc = dmp_LC(g, K)
    c = dmp_pow(lc, d, v, K)

    c = dmp_neg(c, v, K)

    while not dmp_zero_p(h, u):
        k = dmp_degree_in(h, 0, u)
        R.append(h)

        f, g, m, d = g, h, k, m - k

        b = dmp_mul(dmp_neg(lc, v, K),
                    dmp_pow(c, d, v, K), v, K)

        h = dmp_prem(f, g, u, K)
        h = [dmp_quo(ch, b, v, K) for ch in h]

        lc = dmp_LC(g, K)

        p = dmp_pow(dmp_neg(lc, v, K), d, v, K)
        q = dmp_pow(c, d - 1, v, K)
        c = dmp_quo(p, q, v, K)

    return R, dmp_neg(c, v, K)


def dmp_prs_resultant(f, g, u, K):
    """
    Resultant algorithm in `K[X]` using subresultant PRS.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> res, prs = f.resultant(g, includePRS=True)

    >>> res == b  # resultant has n-1 variables
    False
    >>> res == b.drop(x)
    True
    >>> prs == [f, g, a, b]
    True

    """
    if not u:
        return dup_prs_resultant(f, g, K)

    R, S = dmp_inner_subresultants(f, g, u, K)

    if len(R) < 2 or dmp_degree_in(R[-1], 0, u) > 0:
        return dmp_zero(u - 1), R

    return S, R


def dmp_zz_modular_resultant(f, g, p, u, K):
    """
    Compute resultant of `f` and `g` modulo a prime `p`.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> f = x + y + 2
    >>> g = 2*x*y + x + 3

    >>> R.dmp_zz_modular_resultant(f, g, 5)
    -2*y**2 + 1

    """
    if not u:
        return symmetric_residue(dup_prs_resultant(f, g, K)[0] % p, p)

    v = u - 1

    n = dmp_degree_in(f, 0, u)
    m = dmp_degree_in(g, 0, u)

    N = dmp_degree_in(f, 1, u)
    M = dmp_degree_in(g, 1, u)

    B = n*M + m*N

    D, a = [K.one], -K.one
    r = dmp_zero(v)

    while dmp_degree_in(D, 0, 0) <= B:
        while True:
            a += K.one

            if a == p:
                raise HomomorphismFailed('no luck')

            F = dmp_eval_in(f, symmetric_residue(a, p), 1, u, K)

            if dmp_degree_in(F, 0, v) == n:
                G = dmp_eval_in(g, symmetric_residue(a, p), 1, u, K)

                if dmp_degree_in(G, 0, v) == m:
                    break

        R = dmp_zz_modular_resultant(F, G, p, v, K)
        e = dmp_eval_in(r, a, 0, v, K)

        if not v:
            R = dmp_strip([R], 0)
            e = dmp_strip([e], 0)
        else:
            R = [R]
            e = [e]

        d = K.invert(dmp_eval_in(D, a, 0, 0, K), p)
        d = dmp_mul_ground(D, d, 0, K)
        d = dmp_raise(d, v, 0, K)

        c = dmp_mul(d, dmp_sub(R, e, v, K), v, K)
        r = dmp_add(r, c, v, K)

        r = dmp_ground_trunc(r, p, v, K)

        D = dup_mul(D, [K.one, -a], K)
        D = dmp_ground_trunc(D, p, 0, K)

    return r


def _collins_crt(r, R, P, p, K):
    """Wrapper of CRT for Collins's resultant algorithm."""
    return K(crt([P, p], [r, R], check=False, symmetric=True)[0])


def dmp_zz_collins_resultant(f, g, u, K):
    """
    Collins's modular resultant algorithm in `Z[X]`.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> f = x + y + 2
    >>> g = 2*x*y + x + 3

    >>> f.resultant(g)
    -2*y**2 - 5*y + 1

    """
    n = dmp_degree_in(f, 0, u)
    m = dmp_degree_in(g, 0, u)

    if n < 0 or m < 0:
        return dmp_zero(u - 1)

    A = dmp_max_norm(f, u, K)
    B = dmp_max_norm(g, u, K)

    a = dmp_ground_LC(f, u, K)
    b = dmp_ground_LC(g, u, K)

    v = u - 1

    B = K(2)*K.factorial(K(n + m))*A**m*B**n
    r, p, P = dmp_zero(v), K.one, K.one

    while P <= B:
        p = K(nextprime(p))

        while not (a % p) or not (b % p):
            p = K(nextprime(p))

        F = dmp_ground_trunc(f, p, u, K)
        G = dmp_ground_trunc(g, p, u, K)

        try:
            R = dmp_zz_modular_resultant(F, G, p, u, K)
        except HomomorphismFailed:
            continue

        if P == K.one:
            r = R
        else:
            r = dmp_apply_pairs(r, R, _collins_crt, (P, p, K), v, K)

        P *= p

    return r


def dmp_qq_collins_resultant(f, g, u, K0):
    """
    Collins's modular resultant algorithm in `Q[X]`.

    Examples
    ========

    >>> R, x, y = ring('x y', QQ)

    >>> f = x/2 + y + QQ(2, 3)
    >>> g = 2*x*y + x + 3

    >>> f.resultant(g)
    -2*y**2 - 7/3*y + 5/6

    """
    n = dmp_degree_in(f, 0, u)
    m = dmp_degree_in(g, 0, u)

    if n < 0 or m < 0:
        return dmp_zero(u - 1)

    K1 = K0.ring

    cf, f = dmp_clear_denoms(f, u, K0, K1)
    cg, g = dmp_clear_denoms(g, u, K0, K1)

    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)

    r = dmp_zz_collins_resultant(f, g, u, K1)
    r = dmp_convert(r, u - 1, K1, K0)

    c = K0.convert(cf**m * cg**n, K1)

    return dmp_quo_ground(r, c, u - 1, K0)


@cacheit
def dmp_resultant(f, g, u, K, includePRS=False):
    """
    Computes resultant of two polynomials in `K[X]`.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> f.resultant(g)
    -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    """
    if includePRS:
        return dmp_prs_resultant(f, g, u, K)

    if K.is_Field:
        if K.is_RationalField and query('USE_COLLINS_RESULTANT'):
            return dmp_qq_collins_resultant(f, g, u, K)
    else:
        if K.is_IntegerRing and query('USE_COLLINS_RESULTANT'):
            return dmp_zz_collins_resultant(f, g, u, K)

    return dmp_prs_resultant(f, g, u, K)[0]


def dmp_inner_gcd(f, g, u, K):
    """Computes polynomial GCD and cofactors of `f` and `g` in `K[X]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_dense, (f, g))
    return tuple(map(ring.to_dense, f.cofactors(g)))


def dmp_gcd(f, g, u, K):
    """Computes polynomial GCD of `f` and `g` in `K[X]`."""
    return dmp_inner_gcd(f, g, u, K)[0]


def dmp_primitive(f, u, K):
    """Returns multivariate content and a primitive polynomial."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_dense(f)
    new_ring, f = map(lambda _: _.eject(*ring.gens[1:]), (ring, f))
    c, f = f.primitive()
    f = f.inject()
    return new_ring.domain.to_dense(c), ring.to_dense(f)


class _GCD:
    """Mixin class for computing gcd."""

    def _rr_prs_gcd(self, f, g):
        """Computes polynomial GCD using subresultants over a ring."""
        ring = self

        if self.is_multivariate:
            ring, f, g = map(lambda _: _.eject(*self.gens[1:]), (ring, f, g))
            return tuple(map(lambda _: _.inject(), ring._rr_prs_gcd(f, g)))

        domain = ring.domain

        fc, ff = f.primitive()
        gc, fg = g.primitive()

        h = ff.subresultants(fg)[-1]
        _, h = h.primitive()

        c = domain.gcd(fc, gc)
        h *= c

        return h, f // h, g // h

    def _ff_prs_gcd(self, f, g):
        """Computes polynomial GCD using subresultants over a field."""
        ring = self

        if ring.is_multivariate:
            ring, f, g = map(lambda _: _.eject(*self.gens[1:]), (ring, f, g))
            h, f, g = map(lambda _: _.inject(), ring._rr_prs_gcd(f, g))
            c, h = h.LC, h.monic()
            return h, f.quo_ground(c), g.quo_ground(c)

        h = f.subresultants(g)[-1]
        h = h.monic()

        return h, f // h, g // h
