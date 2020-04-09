"""Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences."""

from ..core import cacheit
from ..ntheory import nextprime
from ..ntheory.modular import crt, symmetric_residue
from .densearith import (dmp_add, dmp_max_norm, dmp_mul, dmp_mul_ground,
                         dmp_quo_ground, dmp_sub, dup_mul)
from .densebasic import (dmp_apply_pairs, dmp_convert, dmp_degree_in,
                         dmp_ground_LC, dmp_raise, dmp_strip, dmp_zero)
from .densetools import dmp_clear_denoms, dmp_eval_in, dmp_ground_trunc
from .polyconfig import query
from .polyerrors import HomomorphismFailed


def dup_gcdex(f, g, K):
    """Extended Euclidean algorithm in `F[x]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(1)])
    f, g = map(ring.from_dense, (f, g))
    return tuple(map(ring.to_dense, f.gcdex(g)))


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
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_dense, (f, g))
    res = ring._primitive_prs(f, g)
    res0 = res[0]
    if ring.is_multivariate:
        res0 = ring.drop(0).to_dense(res0)
    return res0, list(map(ring.to_dense, res[1]))


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
        return symmetric_residue(dmp_prs_resultant(f, g, 0, K)[0] % p, p)

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

    def _primitive_prs(self, f, g):
        """
        Subresultant PRS algorithm in `K[X]`.

        Computes the last non-zero scalar subresultant of `f` and `g`
        and subresultant polynomial remainder sequence (PRS).

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
        * :cite:`Geddes1992algorithms`, example 7.6

        """
        ring = self
        domain = ring.domain

        if ring.is_multivariate:
            ring, f, g = map(lambda _: _.eject(*ring.gens[1:]), (ring, f, g))
            res = ring._primitive_prs(f, g)
            return res[0], [_.inject() for _ in res[1]]

        n = f.degree()
        m = g.degree()

        if n < m:
            f, g = g, f
            n, m = m, n

        c, r = domain.zero, []

        if f.is_zero:
            return c, r

        r.append(f)

        if g.is_zero:
            return c, r

        r.append(g)
        d = n - m

        h = f.prem(g)
        h *= (-domain.one)**(d + 1)

        lc = g.LC
        c = -lc**d

        while not h.is_zero:
            k = h.degree()
            r.append(h)

            f, g, m, d = g, h, k, m - k

            h = f.prem(g)
            h = h.quo_ground(-lc*c**d)

            lc = g.LC
            c = domain.quo((-lc)**d, c**(d - 1))

        if r[-1].degree() > 0:
            c = domain.zero
        else:
            c = -c

        return c, r
