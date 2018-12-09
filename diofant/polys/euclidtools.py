"""Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences. """

from ..core import cacheit
from ..ntheory import nextprime
from .densearith import (dmp_add, dmp_div, dmp_max_norm, dmp_mul,
                         dmp_mul_ground, dmp_mul_term, dmp_neg, dmp_pow,
                         dmp_prem, dmp_quo, dmp_quo_ground, dmp_rem, dmp_sub,
                         dmp_sub_mul, dup_mul)
from .densebasic import (dmp_apply_pairs, dmp_convert, dmp_degree_in,
                         dmp_ground, dmp_ground_LC, dmp_inflate, dmp_LC,
                         dmp_multi_deflate, dmp_one, dmp_one_p, dmp_raise,
                         dmp_strip, dmp_zero, dmp_zero_p, dmp_zeros)
from .densetools import (dmp_clear_denoms, dmp_diff_in, dmp_eval_in,
                         dmp_ground_monic, dmp_ground_primitive,
                         dmp_ground_trunc, dup_trunc)
from .galoistools import gf_crt, gf_int
from .heuristicgcd import heugcd
from .polyconfig import query
from .polyerrors import (DomainError, HeuristicGCDFailed, HomomorphismFailed,
                         NotInvertible)


def dup_half_gcdex(f, g, K):
    """
    Half extended Euclidean algorithm in `F[x]`.

    Returns ``(s, h)`` such that ``h = gcd(f, g)`` and ``s*f = h (mod g)``.

    Examples
    ========

    >>> R, x = ring("x", QQ)

    >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
    >>> g = x**3 + x**2 - 4*x - 4

    >>> R.dup_half_gcdex(f, g)
    (-1/5*x + 3/5, x + 1)
    """
    if not K.is_Field:
        raise DomainError("can't compute half extended GCD over %s" % K)

    a, b = [K.one], []

    while g:
        q, r = dmp_div(f, g, 0, K)
        f, g = g, r
        a, b = b, dmp_sub_mul(a, q, b, 0, K)

    a = dmp_quo_ground(a, dmp_LC(f, K), 0, K)
    f = dmp_ground_monic(f, 0, K)

    return a, f


def dup_gcdex(f, g, K):
    """
    Extended Euclidean algorithm in `F[x]`.

    Returns ``(s, t, h)`` such that ``h = gcd(f, g)`` and ``s*f + t*g = h``.

    Examples
    ========

    >>> R, x = ring("x", QQ)

    >>> f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
    >>> g = x**3 + x**2 - 4*x - 4

    >>> R.dup_gcdex(f, g)
    (-1/5*x + 3/5, 1/5*x**2 - 6/5*x + 2, x + 1)

    """
    s, h = dup_half_gcdex(f, g, K)

    F = dmp_sub_mul(h, s, f, 0, K)
    t = dmp_quo(F, g, 0, K)

    return s, t, h


def dup_invert(f, g, K):
    """
    Compute multiplicative inverse of `f` modulo `g` in `F[x]`.

    Examples
    ========

    >>> R, x = ring("x", QQ)

    >>> f = x**2 - 1
    >>> g = 2*x - 1
    >>> h = x - 1

    >>> R.dup_invert(f, g)
    -4/3

    >>> R.dup_invert(f, h)
    Traceback (most recent call last):
    ...
    NotInvertible: zero divisor

    """
    s, h = dup_half_gcdex(f, g, K)

    if h == [K.one]:
        return dmp_rem(s, g, 0, K)
    else:
        raise NotInvertible("zero divisor")


def dup_euclidean_prs(f, g, K):
    """
    Euclidean polynomial remainder sequence (PRS) in `K[x]`.

    Examples
    ========

    >>> R, x = ring("x", QQ)

    >>> f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    >>> prs = R.dup_euclidean_prs(f, g)

    >>> prs[0]
    x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> prs[1]
    3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21
    >>> prs[2]
    -5/9*x**4 + 1/9*x**2 - 1/3
    >>> prs[3]
    -117/25*x**2 - 9*x + 441/25
    >>> prs[4]
    233150/19773*x - 102500/6591
    >>> prs[5]
    -1288744821/543589225

    """
    prs = [f, g]
    h = dmp_rem(f, g, 0, K)

    while h:
        prs.append(h)
        f, g = g, h
        h = dmp_rem(f, g, 0, K)

    return prs


def dup_primitive_prs(f, g, K):
    """
    Primitive polynomial remainder sequence (PRS) in `K[x]`.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    >>> prs = R.dup_primitive_prs(f, g)

    >>> prs[0]
    x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    >>> prs[1]
    3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21
    >>> prs[2]
    -5*x**4 + x**2 - 3
    >>> prs[3]
    13*x**2 + 25*x - 49
    >>> prs[4]
    4663*x - 6150
    >>> prs[5]
    1

    """
    prs = [f, g]
    _, h = dmp_ground_primitive(dmp_prem(f, g, 0, K), 0, K)

    while h:
        prs.append(h)
        f, g = g, h
        _, h = dmp_ground_primitive(dmp_prem(f, g, 0, K), 0, K)

    return prs


def dup_inner_subresultants(f, g, K):
    """
    Subresultant PRS algorithm in `K[x]`.

    Computes the subresultant polynomial remainder sequence (PRS)
    and the non-zero scalar subresultants of `f` and `g`.
    By [1] Thm. 3, these are the constants '-c' (- to optimize
    computation of sign).
    The first subdeterminant is set to 1 by convention to match
    the polynomial and the scalar subdeterminants.
    If 'deg(f) < deg(g)', the subresultants of '(g,f)' are computed.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_inner_subresultants(x**2 + 1, x**2 - 1)
    ([x**2 + 1, x**2 - 1, -2], [1, 1, 4])

    References
    ==========

    * [Brown78]_
    """
    n = dmp_degree_in(f, 0, 0)
    m = dmp_degree_in(g, 0, 0)

    if n < m:
        f, g = g, f
        n, m = m, n

    if not f:
        return [], []

    if not g:
        return [f], [K.one]

    R = [f, g]
    d = n - m

    b = (-K.one)**(d + 1)

    h = dmp_prem(f, g, 0, K)
    h = dmp_mul_ground(h, b, 0, K)

    lc = dmp_LC(g, K)
    c = lc**d

    # Conventional first scalar subdeterminant is 1
    S = [K.one, c]
    c = -c

    while h:
        k = dmp_degree_in(h, 0, 0)
        R.append(h)

        f, g, m, d = g, h, k, m - k

        b = -lc * c**d

        h = dmp_prem(f, g, 0, K)
        h = dmp_quo_ground(h, b, 0, K)

        lc = dmp_LC(g, K)

        if d > 1:        # abnormal case
            q = c**(d - 1)
            c = K.quo((-lc)**d, q)
        else:
            c = -lc

        S.append(-c)

    return R, S


def dup_prs_resultant(f, g, K):
    """
    Resultant algorithm in `K[x]` using subresultant PRS.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_prs_resultant(x**2 + 1, x**2 - 1)
    (4, [x**2 + 1, x**2 - 1, -2])
    """
    if not f or not g:
        return K.zero, []

    R, S = dup_inner_subresultants(f, g, K)

    if dmp_degree_in(R[-1], 0, 0) > 0:
        return K.zero, R

    return S[-1], R


def dup_resultant(f, g, K, includePRS=False):
    """
    Computes resultant of two polynomials in `K[x]`.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_resultant(x**2 + 1, x**2 - 1)
    4
    """
    if includePRS:
        return dup_prs_resultant(f, g, K)
    return dup_prs_resultant(f, g, K)[0]


def dmp_inner_subresultants(f, g, u, K):
    """
    Subresultant PRS algorithm in `K[X]`.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> prs = [f, g, a, b]
    >>> sres = [[1], [1], [3, 0, 0, 0, 0], [-3, 0, 0, -12, 1, 0, -54, 8, 729, -216, 16]]

    >>> R.dmp_inner_subresultants(f, g) == (prs, sres)
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
        return [], []

    v = u - 1
    if dmp_zero_p(g, u):
        return [f], [dmp_ground(K.one, v)]

    R = [f, g]
    d = n - m

    b = dmp_pow(dmp_ground(-K.one, v), d + 1, v, K)

    h = dmp_prem(f, g, u, K)
    h = dmp_mul_term(h, b, 0, u, K)

    lc = dmp_LC(g, K)
    c = dmp_pow(lc, d, v, K)

    S = [dmp_ground(K.one, v), c]
    c = dmp_neg(c, v, K)

    while not dmp_zero_p(h, u):
        k = dmp_degree_in(h, 0, u)
        R.append(h)

        f, g, m, d = g, h, k, m - k

        b = dmp_mul(dmp_neg(lc, v, K),
                    dmp_pow(c, d, v, K), v, K)

        h = dmp_prem(f, g, u, K)
        h = [ dmp_quo(ch, b, v, K) for ch in h ]

        lc = dmp_LC(g, K)

        if d > 1:
            p = dmp_pow(dmp_neg(lc, v, K), d, v, K)
            q = dmp_pow(c, d - 1, v, K)
            c = dmp_quo(p, q, v, K)
        else:
            c = dmp_neg(lc, v, K)

        S.append(dmp_neg(c, v, K))

    return R, S


def dmp_subresultants(f, g, u, K):
    """
    Computes subresultant PRS of two polynomials in `K[X]`.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> R.dmp_subresultants(f, g) == [f, g, a, b]
    True

    """
    return dmp_inner_subresultants(f, g, u, K)[0]


def dmp_prs_resultant(f, g, u, K):
    """
    Resultant algorithm in `K[X]` using subresultant PRS.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> a = 3*x*y**4 + y**3 - 27*y + 4
    >>> b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    >>> res, prs = R.dmp_prs_resultant(f, g)

    >>> res == b             # resultant has n-1 variables
    False
    >>> res == b.drop(x)
    True
    >>> prs == [f, g, a, b]
    True

    """
    if not u:
        return dup_prs_resultant(f, g, K)

    if dmp_zero_p(f, u) or dmp_zero_p(g, u):
        return dmp_zero(u - 1), []

    R, S = dmp_inner_subresultants(f, g, u, K)

    if dmp_degree_in(R[-1], 0, u) > 0:
        return dmp_zero(u - 1), R

    return S[-1], R


def dmp_zz_modular_resultant(f, g, p, u, K):
    """
    Compute resultant of `f` and `g` modulo a prime `p`.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x + y + 2
    >>> g = 2*x*y + x + 3

    >>> R.dmp_zz_modular_resultant(f, g, 5)
    -2*y**2 + 1

    """
    if not u:
        return gf_int(dup_prs_resultant(f, g, K)[0] % p, p)

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

            F = dmp_eval_in(f, gf_int(a, p), 1, u, K)

            if dmp_degree_in(F, 0, v) == n:
                G = dmp_eval_in(g, gf_int(a, p), 1, u, K)

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
        D = dup_trunc(D, p, K)

    return r


def _collins_crt(r, R, P, p, K):
    """Wrapper of CRT for Collins's resultant algorithm. """
    return gf_int(gf_crt([r, R], [P, p], K), P*p)


def dmp_zz_collins_resultant(f, g, u, K):
    """
    Collins's modular resultant algorithm in `Z[X]`.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x + y + 2
    >>> g = 2*x*y + x + 3

    >>> R.dmp_zz_collins_resultant(f, g)
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

    >>> R, x, y = ring("x y", QQ)

    >>> f = x/2 + y + QQ(2, 3)
    >>> g = 2*x*y + x + 3

    >>> R.dmp_qq_collins_resultant(f, g)
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

    >>> R, x, y = ring("x y", ZZ)

    >>> f = 3*x**2*y - y**3 - 4
    >>> g = x**2 + x*y**3 - 9

    >>> R.dmp_resultant(f, g)
    -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    """
    if not u:
        return dup_resultant(f, g, K, includePRS=includePRS)

    if includePRS:
        return dmp_prs_resultant(f, g, u, K)

    if K.is_Field:
        if K.is_RationalField and query('USE_COLLINS_RESULTANT'):
            return dmp_qq_collins_resultant(f, g, u, K)
    else:
        if K.is_IntegerRing and query('USE_COLLINS_RESULTANT'):
            return dmp_zz_collins_resultant(f, g, u, K)

    return dmp_prs_resultant(f, g, u, K)[0]


def dmp_discriminant(f, u, K):
    """
    Computes discriminant of a polynomial in `K[X]`.

    Examples
    ========

    >>> R, x, y, z, t = ring("x y z t", ZZ)

    >>> R.dmp_discriminant(x**2*y + x*z + t)
    -4*y*t + z**2
    """
    d, v = dmp_degree_in(f, 0, u), u - 1

    if d <= 0:
        return dmp_zero(v) if u else K.zero
    else:
        s = (-1)**((d*(d - 1)) // 2)
        c = dmp_LC(f, K)

        r = dmp_resultant(f, dmp_diff_in(f, 1, 0, u, K), u, K)

        if u:
            c = dmp_mul_ground(c, K(s), v, K)
            return dmp_quo(r, c, v, K)
        else:
            return K.quo(r, c*K(s))


def _dmp_rr_trivial_gcd(f, g, u, K):
    """Handle trivial cases in GCD algorithm over a ring. """
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)

    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        if K.is_nonnegative(dmp_ground_LC(g, u, K)):
            return g, dmp_zero(u), dmp_one(u, K)
        else:
            return dmp_neg(g, u, K), dmp_zero(u), dmp_ground(-K.one, u)
    elif zero_g:
        if K.is_nonnegative(dmp_ground_LC(f, u, K)):
            return f, dmp_one(u, K), dmp_zero(u)
        else:
            return dmp_neg(f, u, K), dmp_ground(-K.one, u), dmp_zero(u)
    elif dmp_one_p(f, u, K) or dmp_one_p(g, u, K):
        return dmp_one(u, K), f, g
    elif u and query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)


def _dmp_ff_trivial_gcd(f, g, u, K):
    """Handle trivial cases in GCD algorithm over a field. """
    zero_f = dmp_zero_p(f, u)
    zero_g = dmp_zero_p(g, u)

    if zero_f and zero_g:
        return tuple(dmp_zeros(3, u, K))
    elif zero_f:
        return (dmp_ground_monic(g, u, K),
                dmp_zero(u),
                dmp_ground(dmp_ground_LC(g, u, K), u))
    elif zero_g:
        return (dmp_ground_monic(f, u, K),
                dmp_ground(dmp_ground_LC(f, u, K), u),
                dmp_zero(u))
    elif u and query('USE_SIMPLIFY_GCD'):
        return _dmp_simplify_gcd(f, g, u, K)


def _dmp_simplify_gcd(f, g, u, K):
    """Try to eliminate `x_0` from GCD computation in `K[X]`. """
    df = dmp_degree_in(f, 0, u)
    dg = dmp_degree_in(g, 0, u)

    if df > 0 and dg > 0:
        return

    if not (df or dg):
        F = dmp_LC(f, K)
        G = dmp_LC(g, K)
    else:
        if not df:
            F = dmp_LC(f, K)
            G = dmp_content(g, u, K)
        else:
            F = dmp_content(f, u, K)
            G = dmp_LC(g, K)

    v = u - 1
    h = dmp_gcd(F, G, v, K)

    cff = [ dmp_quo(cf, h, v, K) for cf in f ]
    cfg = [ dmp_quo(cg, h, v, K) for cg in g ]

    return [h], cff, cfg


def dup_rr_prs_gcd(f, g, K):
    """
    Computes polynomial GCD using subresultants over a ring.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_rr_prs_gcd(x**2 - 1, x**2 - 3*x + 2)
    (x - 1, x + 1, x - 2)

    """
    result = _dmp_rr_trivial_gcd(f, g, 0, K)

    if result is not None:
        return result

    fc, F = dmp_ground_primitive(f, 0, K)
    gc, G = dmp_ground_primitive(g, 0, K)

    c = K.gcd(fc, gc)

    h = dmp_subresultants(F, G, 0, K)[-1]
    _, h = dmp_ground_primitive(h, 0, K)

    if K.is_negative(dmp_LC(h, K)):
        c = -c

    h = dmp_mul_ground(h, c, 0, K)

    cff = dmp_quo(f, h, 0, K)
    cfg = dmp_quo(g, h, 0, K)

    return h, cff, cfg


def dup_ff_prs_gcd(f, g, K):
    """
    Computes polynomial GCD using subresultants over a field.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> R, x = ring("x", QQ)

    >>> R.dup_ff_prs_gcd(x**2 - 1, x**2 - 3*x + 2)
    (x - 1, x + 1, x - 2)

    """
    result = _dmp_ff_trivial_gcd(f, g, 0, K)

    if result is not None:
        return result

    h = dmp_subresultants(f, g, 0, K)[-1]
    h = dmp_ground_monic(h, 0, K)

    cff = dmp_quo(f, h, 0, K)
    cfg = dmp_quo(g, h, 0, K)

    return h, cff, cfg


def dmp_rr_prs_gcd(f, g, u, K):
    """
    Computes polynomial GCD using subresultants over a ring.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_rr_prs_gcd(f, g)
    (x + y, x + y, x)

    """
    if not u:
        return dup_rr_prs_gcd(f, g, K)

    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, F = dmp_primitive(f, u, K)
    gc, G = dmp_primitive(g, u, K)

    h = dmp_subresultants(F, G, u, K)[-1]
    c, _, _ = dmp_rr_prs_gcd(fc, gc, u - 1, K)

    if K.is_negative(dmp_ground_LC(h, u, K)):
        h = dmp_neg(h, u, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)

    cff = dmp_quo(f, h, u, K)
    cfg = dmp_quo(g, h, u, K)

    return h, cff, cfg


def dmp_ff_prs_gcd(f, g, u, K):
    """
    Computes polynomial GCD using subresultants over a field.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``, ``cff = quo(f, h)``,
    and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> R, x, y = ring("x y", QQ)

    >>> f = x**2/2 + x*y + y**2/2
    >>> g = x**2 + x*y

    >>> R.dmp_ff_prs_gcd(f, g)
    (x + y, 1/2*x + 1/2*y, x)

    """
    if not u:
        return dup_ff_prs_gcd(f, g, K)

    result = _dmp_ff_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    fc, F = dmp_primitive(f, u, K)
    gc, G = dmp_primitive(g, u, K)

    h = dmp_subresultants(F, G, u, K)[-1]
    c, _, _ = dmp_ff_prs_gcd(fc, gc, u - 1, K)

    _, h = dmp_primitive(h, u, K)
    h = dmp_mul_term(h, c, 0, u, K)
    h = dmp_ground_monic(h, u, K)

    cff = dmp_quo(f, h, u, K)
    cfg = dmp_quo(g, h, u, K)

    return h, cff, cfg


def dmp_zz_heu_gcd(f, g, u, K):
    """
    Heuristic polynomial GCD in `Z[X]`.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_zz_heu_gcd(f, g)
    (x + y, x + y, x)

    See Also
    ========

    diofant.polys.heuristicgcd.heugcd

    References
    ==========

    * [Liao95]_
    """
    result = _dmp_rr_trivial_gcd(f, g, u, K)

    if result is not None:
        return result

    ring = K.poly_ring(*["_%d" % i for i in range(u + 1)])
    f, g = map(ring.from_dense, (f, g))
    h, cff, cfg = heugcd(f, g)
    return tuple(map(ring.to_dense, f.cofactors(g)))


def dmp_qq_heu_gcd(f, g, u, K0):
    """
    Heuristic polynomial GCD in `Q[X]`.

    Returns ``(h, cff, cfg)`` such that ``a = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> R, x, y = ring("x y", QQ)

    >>> f = x**2/4 + x*y + y**2
    >>> g = x**2/2 + x*y

    >>> R.dmp_qq_heu_gcd(f, g)
    (x + 2*y, 1/4*x + 1/2*y, 1/2*x)

    """
    result = _dmp_ff_trivial_gcd(f, g, u, K0)

    if result is not None:
        return result

    K1 = K0.ring

    cf, f = dmp_clear_denoms(f, u, K0, K1)
    cg, g = dmp_clear_denoms(g, u, K0, K1)

    f = dmp_convert(f, u, K0, K1)
    g = dmp_convert(g, u, K0, K1)

    h, cff, cfg = dmp_zz_heu_gcd(f, g, u, K1)

    h = dmp_convert(h, u, K1, K0)

    c = dmp_ground_LC(h, u, K0)
    h = dmp_ground_monic(h, u, K0)

    cff = dmp_convert(cff, u, K1, K0)
    cfg = dmp_convert(cfg, u, K1, K0)

    cff = dmp_mul_ground(cff, K0.quo(c, cf), u, K0)
    cfg = dmp_mul_ground(cfg, K0.quo(c, cg), u, K0)

    return h, cff, cfg


def _dmp_zz_modgcd(f, g, u, K):
    from .modulargcd import modgcd
    ring = K.poly_ring(*["_%d" % i for i in range(u + 1)])
    f, g = map(ring.from_dense, (f, g))
    h, cff, cfg = modgcd(f, g)
    return tuple(map(ring.to_dense, f.cofactors(g)))


_gcd_zz_methods = {'modgcd': _dmp_zz_modgcd,
                   'prs': dmp_rr_prs_gcd}


def _dmp_aa_modgcd(f, g, u, K):
    from .modulargcd import func_field_modgcd
    ring = K.poly_ring(*["_%d" % i for i in range(u + 1)])
    f, g = map(ring.from_dense, (f, g))
    h, cff, cfg = func_field_modgcd(f, g)
    return tuple(map(ring.to_dense, f.cofactors(g)))


_gcd_aa_methods = {'modgcd': _dmp_aa_modgcd,
                   'prs': dmp_ff_prs_gcd}


def _dmp_inner_gcd(f, g, u, K):
    """Helper function for `dmp_inner_gcd()`. """
    if not K.is_Exact:
        try:
            exact = K.get_exact()
        except DomainError:
            return dmp_one(u, K), f, g

        f = dmp_convert(f, u, K, exact)
        g = dmp_convert(g, u, K, exact)

        h, cff, cfg = _dmp_inner_gcd(f, g, u, exact)

        h = dmp_convert(h, u, exact, K)
        cff = dmp_convert(cff, u, exact, K)
        cfg = dmp_convert(cfg, u, exact, K)

        return h, cff, cfg
    elif K.is_Field:
        if K.is_RationalField:
            if query('USE_HEU_GCD'):
                try:
                    return dmp_qq_heu_gcd(f, g, u, K)
                except HeuristicGCDFailed:  # pragma: no cover
                    pass
        elif K.is_AlgebraicField:
            method = _gcd_aa_methods[query('GCD_AA_METHOD')]
            return method(f, g, u, K)

        return dmp_ff_prs_gcd(f, g, u, K)
    else:
        if K.is_IntegerRing:
            if query('USE_HEU_GCD'):
                try:
                    return dmp_zz_heu_gcd(f, g, u, K)
                except HeuristicGCDFailed:  # pragma: no cover
                    pass

            method = _gcd_zz_methods[query('FALLBACK_GCD_ZZ_METHOD')]
            return method(f, g, u, K)

        return dmp_rr_prs_gcd(f, g, u, K)


def dmp_inner_gcd(f, g, u, K):
    """
    Computes polynomial GCD and cofactors of `f` and `g` in `K[X]`.

    Returns ``(h, cff, cfg)`` such that ``h = gcd(f, g)``,
    ``cff = quo(f, h)``, and ``cfg = quo(g, h)``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_inner_gcd(f, g)
    (x + y, x + y, x)

    """
    J, (f, g) = dmp_multi_deflate((f, g), u, K)
    h, cff, cfg = _dmp_inner_gcd(f, g, u, K)

    return (dmp_inflate(h, J, u, K),
            dmp_inflate(cff, J, u, K),
            dmp_inflate(cfg, J, u, K))


def dmp_gcd(f, g, u, K):
    """
    Computes polynomial GCD of `f` and `g` in `K[X]`.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_gcd(f, g)
    x + y

    """
    return dmp_inner_gcd(f, g, u, K)[0]


def dmp_rr_lcm(f, g, u, K):
    """
    Computes polynomial LCM over a ring in `K[X]`.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dmp_rr_lcm(x**2 - 1, x**2 - 3*x + 2)
    x**3 - 2*x**2 - x + 2

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_rr_lcm(f, g)
    x**3 + 2*x**2*y + x*y**2

    """
    fc, f = dmp_ground_primitive(f, u, K)
    gc, g = dmp_ground_primitive(g, u, K)

    c = K.lcm(fc, gc)

    h = dmp_quo(dmp_mul(f, g, u, K),
                dmp_gcd(f, g, u, K), u, K)

    return dmp_mul_ground(h, c, u, K)


def dmp_ff_lcm(f, g, u, K):
    """
    Computes polynomial LCM over a field in `K[X]`.

    Examples
    ========

    >>> R, x = ring("x", QQ)

    >>> f = (x**2 + 7*x/2 + 3)/2
    >>> g = x**2/2 + x

    >>> R.dmp_ff_lcm(f, g)
    x**3 + 7/2*x**2 + 3*x

    >>> R, x, y = ring("x y", QQ)

    >>> f = x**2/4 + x*y + y**2
    >>> g = x**2/2 + x*y

    >>> R.dmp_ff_lcm(f, g)
    x**3 + 4*x**2*y + 4*x*y**2
    """
    h = dmp_quo(dmp_mul(f, g, u, K),
                dmp_gcd(f, g, u, K), u, K)

    return dmp_ground_monic(h, u, K)


def dmp_lcm(f, g, u, K):
    """
    Computes polynomial LCM of `f` and `g` in `K[X]`.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x**2 + 2*x*y + y**2
    >>> g = x**2 + x*y

    >>> R.dmp_lcm(f, g)
    x**3 + 2*x**2*y + x*y**2
    """
    if K.is_Field:
        return dmp_ff_lcm(f, g, u, K)
    else:
        return dmp_rr_lcm(f, g, u, K)


def dmp_content(f, u, K):
    """
    Returns GCD of multivariate coefficients.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_content(2*x*y + 6*x + 4*y + 12)
    2*y + 6

    """
    cont, v = dmp_LC(f, K), u - 1

    if dmp_zero_p(f, u):
        return cont

    for c in f[1:]:
        cont = dmp_gcd(cont, c, v, K)

        if dmp_one_p(cont, v, K):
            break

    if K.is_negative(dmp_ground_LC(cont, v, K)):
        return dmp_neg(cont, v, K)
    else:
        return cont


def dmp_primitive(f, u, K):
    """
    Returns multivariate content and a primitive polynomial.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_primitive(2*x*y + 6*x + 4*y + 12)
    (2*y + 6, x + 2)

    """
    cont, v = dmp_content(f, u, K), u - 1

    if dmp_zero_p(f, u) or dmp_one_p(cont, v, K):
        return cont, f
    else:
        return cont, [ dmp_quo(c, cont, v, K) for c in f ]


def dmp_cancel(f, g, u, K, include=True):
    """
    Cancel common factors in a rational function `f/g`.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_cancel(2*x**2 - 2, x**2 - 2*x + 1)
    (2*x + 2, x - 1)

    """
    K0 = None

    if K.is_Field and K.has_assoc_Ring:
        K0, K = K, K.ring

        cq, f = dmp_clear_denoms(f, u, K0, K, convert=True)
        cp, g = dmp_clear_denoms(g, u, K0, K, convert=True)
    else:
        cp, cq = K.one, K.one

    _, p, q = dmp_inner_gcd(f, g, u, K)

    if K0 is not None:
        _, cp, cq = K.cofactors(cp, cq)

        p = dmp_convert(p, u, K, K0)
        q = dmp_convert(q, u, K, K0)

        K = K0

    p_neg = K.is_negative(dmp_ground_LC(p, u, K))
    q_neg = K.is_negative(dmp_ground_LC(q, u, K))

    if p_neg and q_neg:
        p, q = dmp_neg(p, u, K), dmp_neg(q, u, K)
    elif p_neg:
        cp, p = -cp, dmp_neg(p, u, K)
    elif q_neg:
        cp, q = -cp, dmp_neg(q, u, K)

    if not include:
        return cp, cq, p, q

    p = dmp_mul_ground(p, cp, u, K)
    q = dmp_mul_ground(q, cq, u, K)

    return p, q
