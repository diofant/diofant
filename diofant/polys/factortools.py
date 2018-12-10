"""Polynomial factorization routines in characteristic zero. """

import math

from ..ntheory import factorint, isprime, nextprime
from ..utilities import subsets
from .densearith import (dmp_add, dmp_add_mul, dmp_div, dmp_expand,
                         dmp_l1_norm, dmp_max_norm, dmp_mul, dmp_mul_ground,
                         dmp_neg, dmp_pow, dmp_quo, dmp_quo_ground, dmp_sub,
                         dmp_sub_mul, dup_add, dup_lshift, dup_mul, dup_sqr,
                         dup_sub)
from .densebasic import (dmp_convert, dmp_degree_in, dmp_degree_list,
                         dmp_eject, dmp_exclude, dmp_from_dict, dmp_ground_LC,
                         dmp_ground_p, dmp_include, dmp_inject, dmp_LC,
                         dmp_nest, dmp_normal, dmp_one, dmp_raise, dmp_strip,
                         dmp_swap, dmp_TC, dmp_terms_gcd, dmp_zero_p,
                         dup_inflate)
from .densetools import (dmp_clear_denoms, dmp_compose, dmp_diff_eval_in,
                         dmp_eval_in, dmp_eval_tail, dmp_ground_content,
                         dmp_ground_monic, dmp_ground_primitive,
                         dmp_ground_trunc, dup_mirror, dup_trunc)
from .euclidtools import dmp_inner_gcd, dmp_primitive
from .galoistools import (gf_add_mul, gf_div, gf_factor, gf_factor_sqf,
                          gf_from_int_poly, gf_gcdex, gf_mul, gf_rem,
                          gf_to_int_poly)
from .polyconfig import query
from .polyerrors import (CoercionFailed, DomainError, EvaluationFailed,
                         ExtraneousFactors)
from .polyutils import _sort_factors
from .sqfreetools import dmp_sqf_norm, dmp_sqf_p, dmp_sqf_part


def dmp_trial_division(f, factors, u, K):
    """Determine multiplicities of factors using trial division. """
    result = []

    for factor in factors:
        k = 0

        while True:
            q, r = dmp_div(f, factor, u, K)

            if dmp_zero_p(r, u):
                f, k = q, k + 1
            else:
                break

        result.append((factor, k))

    return _sort_factors(result)


def dmp_zz_mignotte_bound(f, u, K):
    """Mignotte bound for multivariate polynomials in `K[X]`. """
    a = dmp_max_norm(f, u, K)
    b = abs(dmp_ground_LC(f, u, K))
    n = sum(dmp_degree_list(f, u))

    return K.sqrt(K(n + 1))*2**n*a*b


def dup_zz_hensel_step(m, f, g, h, s, t, K):
    """
    One step in Hensel lifting in `Z[x]`.

    Given positive integer `m` and `Z[x]` polynomials `f`, `g`, `h`, `s`
    and `t` such that::

        f == g*h (mod m)
        s*g + t*h == 1 (mod m)

        lc(f) is not a zero divisor (mod m)
        lc(h) == 1

        deg(f) == deg(g) + deg(h)
        deg(s) < deg(h)
        deg(t) < deg(g)

    returns polynomials `G`, `H`, `S` and `T`, such that::

        f == G*H (mod m**2)
        S*G + T**H == 1 (mod m**2)

    References
    ==========

    * [Gathen99]_
    """
    M = m**2

    e = dmp_sub_mul(f, g, h, 0, K)
    e = dup_trunc(e, M, K)

    q, r = dmp_div(dup_mul(s, e, K), h, 0, K)

    q = dup_trunc(q, M, K)
    r = dup_trunc(r, M, K)

    u = dup_add(dup_mul(t, e, K), dup_mul(q, g, K), K)
    G = dup_trunc(dup_add(g, u, K), M, K)
    H = dup_trunc(dup_add(h, r, K), M, K)

    u = dup_add(dup_mul(s, G, K), dup_mul(t, H, K), K)
    b = dup_trunc(dup_sub(u, [K.one], K), M, K)

    c, d = dmp_div(dup_mul(s, b, K), H, 0, K)

    c = dup_trunc(c, M, K)
    d = dup_trunc(d, M, K)

    u = dup_add(dup_mul(t, b, K), dup_mul(c, G, K), K)
    S = dup_trunc(dup_sub(s, d, K), M, K)
    T = dup_trunc(dup_sub(t, u, K), M, K)

    return G, H, S, T


def dup_zz_hensel_lift(p, f, f_list, l, K):
    """
    Multifactor Hensel lifting in `Z[x]`.

    Given a prime `p`, polynomial `f` over `Z[x]` such that `lc(f)`
    is a unit modulo `p`, monic pair-wise coprime polynomials `f_i`
    over `Z[x]` satisfying::

        f = lc(f) f_1 ... f_r (mod p)

    and a positive integer `l`, returns a list of monic polynomials
    `F_1`, `F_2`, ..., `F_r` satisfying::

       f = lc(f) F_1 ... F_r (mod p**l)

       F_i = f_i (mod p), i = 1..r

    References
    ==========

    * [Gathen99]_
    """
    r = len(f_list)
    lc = dmp_LC(f, K)

    if r == 1:
        F = dmp_mul_ground(f, K.gcdex(lc, p**l)[0], 0, K)
        return [ dup_trunc(F, p**l, K) ]

    m = p
    k = r // 2
    d = int(math.ceil(math.log(l, 2)))

    g = gf_from_int_poly([lc], p)

    for f_i in f_list[:k]:
        g = gf_mul(g, gf_from_int_poly(f_i, p), p, K)

    h = gf_from_int_poly(f_list[k], p)

    for f_i in f_list[k + 1:]:
        h = gf_mul(h, gf_from_int_poly(f_i, p), p, K)

    s, t, _ = gf_gcdex(g, h, p, K)

    g = gf_to_int_poly(g, p)
    h = gf_to_int_poly(h, p)
    s = gf_to_int_poly(s, p)
    t = gf_to_int_poly(t, p)

    for _ in range(1, d + 1):
        (g, h, s, t), m = dup_zz_hensel_step(m, f, g, h, s, t, K), m**2

    return dup_zz_hensel_lift(p, g, f_list[:k], l, K) \
        + dup_zz_hensel_lift(p, h, f_list[k:], l, K)


def _test_pl(fc, q, pl):
    if q > pl // 2:
        q = q - pl
    if not q:
        return True
    return fc % q == 0


def dup_zz_zassenhaus(f, K):
    """Factor primitive square-free polynomials in `Z[x]`. """
    n = dmp_degree_in(f, 0, 0)

    if n == 1:
        return [f]

    fc = f[-1]
    A = dmp_max_norm(f, 0, K)
    b = dmp_LC(f, K)
    B = int(abs(K.sqrt(K(n + 1))*2**n*A*b))
    C = int((n + 1)**(2*n)*A**(2*n - 1))
    gamma = int(math.ceil(2*math.log(C, 2)))
    bound = int(2*gamma*math.log(gamma))
    a = []
    # choose a prime number `p` such that `f` be square free in Z_p
    # if there are many factors in Z_p, choose among a few different `p`
    # the one with fewer factors
    for px in range(3, bound + 1):  # pragma: no branch
        if not isprime(px) or b % px == 0:
            continue

        px = K.convert(px)
        F = dmp_normal(f, 0, K.finite_field(px))

        if not dmp_sqf_p(F, 0, K.finite_field(px)):
            continue

        F = gf_from_int_poly(f, px)
        fsqfx = gf_factor_sqf(F, px, K)[1]
        a.append((px, fsqfx))
        if len(fsqfx) < 15 or len(a) > 4:
            break
    p, fsqf = min(a, key=lambda x: len(x[1]))

    l = int(math.ceil(math.log(2*B + 1, p)))

    modular = [gf_to_int_poly(ff, p) for ff in fsqf]

    g = dup_zz_hensel_lift(p, f, modular, l, K)

    sorted_T = range(len(g))
    T = set(sorted_T)
    factors, s = [], 1
    pl = p**l

    while 2*s <= len(T):
        for S in subsets(sorted_T, s):
            # lift the constant coefficient of the product `G` of the factors
            # in the subset `S`; if it is does not divide `fc`, `G` does
            # not divide the input polynomial

            if b == 1:
                q = 1
                for i in S:
                    q = q*g[i][-1]
                q = q % pl
                if not _test_pl(fc, q, pl):
                    continue
            else:
                G = [b]
                for i in S:
                    G = dup_mul(G, g[i], K)
                G = dup_trunc(G, pl, K)
                G = dmp_ground_primitive(G, 0, K)[1]
                q = G[-1]
                if q and fc % q != 0:
                    continue

            H = [b]
            S = set(S)
            T_S = T - S

            if b == 1:
                G = [b]
                for i in S:
                    G = dup_mul(G, g[i], K)
                G = dup_trunc(G, pl, K)

            for i in T_S:
                H = dup_mul(H, g[i], K)

            H = dup_trunc(H, pl, K)

            G_norm = dmp_l1_norm(G, 0, K)
            H_norm = dmp_l1_norm(H, 0, K)

            if G_norm*H_norm <= B:
                T = T_S
                sorted_T = [i for i in sorted_T if i not in S]

                G = dmp_ground_primitive(G, 0, K)[1]
                f = dmp_ground_primitive(H, 0, K)[1]

                factors.append(G)
                b = dmp_LC(f, K)

                break
        else:
            s += 1

    return factors + [f]


def dup_zz_irreducible_p(f, K):
    """Test irreducibility using Eisenstein's criterion. """
    lc = dmp_LC(f, K)
    tc = dmp_TC(f, K)

    e_fc = dmp_ground_content(f[1:], 0, K)

    if e_fc:
        e_ff = factorint(int(e_fc))

        for p in e_ff:
            if (lc % p) and (tc % p**2):
                return True


def dup_cyclotomic_p(f, K, irreducible=False):
    """
    Efficiently test if ``f`` is a cyclotomic polynomial.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> f = x**16 + x**14 - x**10 + x**8 - x**6 + x**2 + 1
    >>> R.dup_cyclotomic_p(f)
    False

    >>> g = x**16 + x**14 - x**10 - x**8 - x**6 + x**2 + 1
    >>> R.dup_cyclotomic_p(g)
    True
    """
    if K.is_RationalField:
        try:
            K0, K = K, K.ring
            f = dmp_convert(f, 0, K0, K)
        except CoercionFailed:
            return False
    elif not K.is_IntegerRing:
        return False

    lc = dmp_LC(f, K)
    tc = dmp_TC(f, K)

    if lc != 1 or (tc != -1 and tc != 1):
        return False

    if not irreducible:
        coeff, factors = dmp_factor_list(f, 0, K)

        if coeff != K.one or factors != [(f, 1)]:
            return False

    n = dmp_degree_in(f, 0, 0)
    g, h = [], []

    for i in range(n, -1, -2):
        g.insert(0, f[i])

    for i in range(n - 1, -1, -2):
        h.insert(0, f[i])

    g = dup_sqr(dmp_strip(g, 0), K)
    h = dup_sqr(dmp_strip(h, 0), K)

    F = dup_sub(g, dup_lshift(h, 1, K), K)

    if K.is_negative(dmp_LC(F, K)):
        F = dmp_neg(F, 0, K)

    if F == f:
        return True

    g = dup_mirror(f, K)

    if K.is_negative(dmp_LC(g, K)):
        g = dmp_neg(g, 0, K)

    if F == g and dup_cyclotomic_p(g, K):
        return True

    G = dmp_sqf_part(F, 0, K)

    if dup_sqr(G, K) == F and dup_cyclotomic_p(G, K):
        return True

    return False


def dup_zz_cyclotomic_poly(n, K):
    """Efficiently generate n-th cyclotomic polynomial. """
    h = [K.one, -K.one]

    for p, k in factorint(n).items():
        h = dmp_quo(dup_inflate(h, p, K), h, 0, K)
        h = dup_inflate(h, p**(k - 1), K)

    return h


def _dup_cyclotomic_decompose(n, K):
    H = [[K.one, -K.one]]

    for p, k in factorint(n).items():
        Q = [dmp_quo(dup_inflate(h, p, K), h, 0, K) for h in H]
        H.extend(Q)

        for i in range(1, k):
            Q = [ dup_inflate(q, p, K) for q in Q ]
            H.extend(Q)

    return H


def dup_zz_cyclotomic_factor(f, K):
    """
    Efficiently factor polynomials `x**n - 1` and `x**n + 1` in `Z[x]`.

    Given a univariate polynomial `f` in `Z[x]` returns a list of factors
    of `f`, provided that `f` is in the form `x**n - 1` or `x**n + 1` for
    `n >= 1`. Otherwise returns None.

    Factorization is performed using using cyclotomic decomposition of `f`,
    which makes this method much faster that any other direct factorization
    approach (e.g. Zassenhaus's).

    References
    ==========

    * [Weisstein09]_
    """
    lc_f, tc_f = dmp_LC(f, K), dmp_TC(f, K)

    if dmp_ground_p(f, None, 0):
        return

    if lc_f != 1 or tc_f not in [-1, 1]:
        return

    if any(bool(cf) for cf in f[1:-1]):
        return

    n = dmp_degree_in(f, 0, 0)
    F = _dup_cyclotomic_decompose(n, K)

    if tc_f != K.one:
        return F
    else:
        H = []

        for h in _dup_cyclotomic_decompose(2*n, K):
            if h not in F:
                H.append(h)

        return H


def dup_zz_factor_sqf(f, K):
    """Factor square-free (non-primitive) polynomials in `Z[x]`. """
    cont, g = dmp_ground_primitive(f, 0, K)

    n = dmp_degree_in(g, 0, 0)

    if dmp_LC(g, K) < 0:
        cont, g = -cont, dmp_neg(g, 0, K)

    if n <= 0:
        return cont, []
    elif n == 1:
        return cont, [g]

    if query('USE_IRREDUCIBLE_IN_FACTOR'):
        if dup_zz_irreducible_p(g, K):
            return cont, [g]

    factors = None

    if query('USE_CYCLOTOMIC_FACTOR'):
        factors = dup_zz_cyclotomic_factor(g, K)

    if factors is None:
        factors = dup_zz_zassenhaus(g, K)

    return cont, _sort_factors(factors, multiple=False)


def dup_zz_factor(f, K):
    """
    Factor (non square-free) polynomials in `Z[x]`.

    Given a univariate polynomial `f` in `Z[x]` computes its complete
    factorization `f_1, ..., f_n` into irreducibles over integers::

                f = content(f) f_1**k_1 ... f_n**k_n

    The factorization is computed by reducing the input polynomial
    into a primitive square-free polynomial and factoring it using
    Zassenhaus algorithm. Trial division is used to recover the
    multiplicities of factors.

    The result is returned as a tuple consisting of::

              (content(f), [(f_1, k_1), ..., (f_n, k_n))

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_zz_factor(2*x**4 - 2)
    (2, [(x - 1, 1), (x + 1, 1), (x**2 + 1, 1)])

    Note that this is a complete factorization over integers,
    however over Gaussian integers we can factor the last term.

    By default, polynomials `x**n - 1` and `x**n + 1` are factored
    using cyclotomic decomposition to speedup computations. To
    disable this behaviour set cyclotomic=False.

    References
    ==========

    * [Gathen99]_
    """
    cont, g = dmp_ground_primitive(f, 0, K)

    n = dmp_degree_in(g, 0, 0)

    if dmp_LC(g, K) < 0:
        cont, g = -cont, dmp_neg(g, 0, K)

    if n <= 0:
        return cont, []
    elif n == 1:
        return cont, [(g, 1)]

    if query('USE_IRREDUCIBLE_IN_FACTOR'):
        if dup_zz_irreducible_p(g, K):
            return cont, [(g, 1)]

    g = dmp_sqf_part(g, 0, K)
    H = None

    if query('USE_CYCLOTOMIC_FACTOR'):
        H = dup_zz_cyclotomic_factor(g, K)

    if H is None:
        H = dup_zz_zassenhaus(g, K)

    factors = dmp_trial_division(f, H, 0, K)
    return cont, factors


def dmp_zz_wang_non_divisors(E, cs, ct, K):
    """Wang/EEZ: Compute a set of valid divisors.  """
    result = [ cs*ct ]

    for q in E:
        q = abs(q)

        for r in reversed(result):
            while r != 1:
                r = K.gcd(r, q)
                q = q // r

            if q == K.one:
                return

        result.append(q)

    return result[1:]


def dmp_zz_wang_test_points(f, T, ct, A, u, K):
    """Wang/EEZ: Test evaluation points for suitability. """
    if not dmp_eval_tail(dmp_LC(f, K), A, u - 1, K):
        raise EvaluationFailed('no luck')

    g = dmp_eval_tail(f, A, u, K)

    if not dmp_sqf_p(g, 0, K):
        raise EvaluationFailed('no luck')

    c, h = dmp_ground_primitive(g, 0, K)

    if K.is_negative(dmp_LC(h, K)):
        c, h = -c, dmp_neg(h, 0, K)

    v = u - 1

    E = [ dmp_eval_tail(t, A, v, K) for t, _ in T ]
    D = dmp_zz_wang_non_divisors(E, c, ct, K)

    if D is not None:
        return c, h, E
    else:
        raise EvaluationFailed('no luck')


def dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K):
    """Wang/EEZ: Compute correct leading coefficients. """
    C, J, v = [], [0]*len(E), u - 1

    for h in H:
        c = dmp_one(v, K)
        d = dmp_LC(h, K)*cs

        for i in reversed(range(len(E))):
            k, e, t = 0, E[i], T[i][0]

            while not (d % e):
                d, k = d//e, k + 1

            if k != 0:
                c, J[i] = dmp_mul(c, dmp_pow(t, k, v, K), v, K), 1

        C.append(c)

    if any(not j for j in J):  # pragma: no cover
        raise ExtraneousFactors

    CC, HH = [], []

    for c, h in zip(C, H):
        d = dmp_eval_tail(c, A, v, K)
        lc = dmp_LC(h, K)

        if cs == K.one:
            cc = lc//d
        else:
            g = K.gcd(lc, d)
            d, cc = d//g, lc//g
            h, cs = dmp_mul_ground(h, d, 0, K), cs//d

        c = dmp_mul_ground(c, cc, v, K)

        CC.append(c)
        HH.append(h)

    if cs == K.one:
        return f, HH, CC

    CCC, HHH = [], []

    for c, h in zip(CC, HH):
        CCC.append(dmp_mul_ground(c, cs, v, K))
        HHH.append(dmp_mul_ground(h, cs, 0, K))

    f = dmp_mul_ground(f, cs**(len(H) - 1), u, K)

    return f, HHH, CCC


def dup_zz_diophantine(F, m, p, K):
    """Wang/EEZ: Solve univariate Diophantine equations. """
    if len(F) == 2:
        a, b = F

        f = gf_from_int_poly(a, p)
        g = gf_from_int_poly(b, p)

        s, t, G = gf_gcdex(g, f, p, K)

        s = dup_lshift(s, m, K)
        t = dup_lshift(t, m, K)

        q, s = gf_div(s, f, p, K)

        t = gf_add_mul(t, q, g, p, K)

        s = gf_to_int_poly(s, p)
        t = gf_to_int_poly(t, p)

        result = [s, t]
    else:
        G = [F[-1]]

        for f in reversed(F[1:-1]):
            G.insert(0, dup_mul(f, G[0], K))

        S, T = [], [[1]]

        for f, g in zip(F, G):
            t, s = dmp_zz_diophantine([g, f], T[-1], [], 0, p, 1, K)
            T.append(t)
            S.append(s)

        result, S = [], S + [T[-1]]

        for s, f in zip(S, F):
            s = gf_from_int_poly(s, p)
            f = gf_from_int_poly(f, p)

            r = gf_rem(dup_lshift(s, m, K), f, p, K)
            s = gf_to_int_poly(r, p)

            result.append(s)

    return result


def dmp_zz_diophantine(F, c, A, d, p, u, K):
    """Wang/EEZ: Solve multivariate Diophantine equations. """
    if not A:
        S = [ [] for _ in F ]
        n = dmp_degree_in(c, 0, 0)

        for i, coeff in enumerate(c):
            if not coeff:
                continue

            T = dup_zz_diophantine(F, n - i, p, K)

            for j, (s, t) in enumerate(zip(S, T)):
                t = dmp_mul_ground(t, coeff, 0, K)
                S[j] = dup_trunc(dup_add(s, t, K), p, K)
    else:
        n = len(A)
        e = dmp_expand(F, u, K)

        a, A = A[-1], A[:-1]
        B, G = [], []

        for f in F:
            B.append(dmp_quo(e, f, u, K))
            G.append(dmp_eval_in(f, a, n, u, K))

        C = dmp_eval_in(c, a, n, u, K)

        v = u - 1

        S = dmp_zz_diophantine(G, C, A, d, p, v, K)
        S = [ dmp_raise(s, 1, v, K) for s in S ]

        for s, b in zip(S, B):
            c = dmp_sub_mul(c, s, b, u, K)

        c = dmp_ground_trunc(c, p, u, K)

        m = dmp_nest([K.one, -a], n, K)
        M = dmp_one(n, K)

        for k in range(d):
            k = K(k)
            if dmp_zero_p(c, u):
                break

            M = dmp_mul(M, m, u, K)
            C = dmp_diff_eval_in(c, k + 1, a, n, u, K)

            if not dmp_zero_p(C, v):
                C = dmp_quo_ground(C, K.factorial(k + 1), v, K)
                T = dmp_zz_diophantine(G, C, A, d, p, v, K)

                for i, t in enumerate(T):
                    T[i] = dmp_mul(dmp_raise(t, 1, v, K), M, u, K)

                for i, (s, t) in enumerate(zip(S, T)):
                    S[i] = dmp_add(s, t, u, K)

                for t, b in zip(T, B):
                    c = dmp_sub_mul(c, t, b, u, K)

                c = dmp_ground_trunc(c, p, u, K)

        S = [ dmp_ground_trunc(s, p, u, K) for s in S ]

    return S


def dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K):
    """Wang/EEZ: Parallel Hensel lifting algorithm. """
    S, n, v = [f], len(A), u - 1

    H = list(H)

    for i, a in enumerate(reversed(A[1:])):
        s = dmp_eval_in(S[0], a, n - i, u - i, K)
        S.insert(0, dmp_ground_trunc(s, p, v - i, K))

    d = max(dmp_degree_list(f, u)[1:])

    for j, s, a in zip(range(2, n + 2), S, A):
        G, w = list(H), j - 1

        I, J = A[:j - 2], A[j - 1:]

        for i, (h, lc) in enumerate(zip(H, LC)):
            lc = dmp_ground_trunc(dmp_eval_tail(lc, J, v, K), p, w - 1, K)
            H[i] = [lc] + dmp_raise(h[1:], 1, w - 1, K)

        m = dmp_nest([K.one, -a], w, K)
        M = dmp_one(w, K)

        c = dmp_sub(s, dmp_expand(H, w, K), w, K)

        dj = dmp_degree_in(s, w, w)

        for k in range(dj):
            k = K(k)
            if dmp_zero_p(c, w):
                break

            M = dmp_mul(M, m, w, K)
            C = dmp_diff_eval_in(c, k + 1, a, w, w, K)

            if not dmp_zero_p(C, w - 1):
                C = dmp_quo_ground(C, K.factorial(k + 1), w - 1, K)
                T = dmp_zz_diophantine(G, C, I, d, p, w - 1, K)

                for i, (h, t) in enumerate(zip(H, T)):
                    h = dmp_add_mul(h, dmp_raise(t, 1, w - 1, K), M, w, K)
                    H[i] = dmp_ground_trunc(h, p, w, K)

                h = dmp_sub(s, dmp_expand(H, w, K), w, K)
                c = dmp_ground_trunc(h, p, w, K)

    if dmp_expand(H, u, K) != f:
        raise ExtraneousFactors  # pragma: no cover
    else:
        return H


def dmp_zz_wang(f, u, K, mod=None, seed=None):
    """
    Factor primitive square-free polynomials in `Z[X]`.

    Given a multivariate polynomial `f` in `Z[x_1,...,x_n]`, which is
    primitive and square-free in `x_1`, computes factorization of `f` into
    irreducibles over integers.

    The procedure is based on Wang's Enhanced Extended Zassenhaus
    algorithm. The algorithm works by viewing `f` as a univariate polynomial
    in `Z[x_2,...,x_n][x_1]`, for which an evaluation mapping is computed::

                      x_2 -> a_2, ..., x_n -> a_n

    where `a_i`, for `i = 2, ..., n`, are carefully chosen integers.  The
    mapping is used to transform `f` into a univariate polynomial in `Z[x_1]`,
    which can be factored efficiently using Zassenhaus algorithm. The last
    step is to lift univariate factors to obtain true multivariate
    factors. For this purpose a parallel Hensel lifting procedure is used.

    The parameter ``seed`` is passed to _randint and can be used to seed randint
    (when an integer) or (for testing purposes) can be a sequence of numbers.

    References
    ==========

    * [Wang78]_
    * [Geddes92]_
    """
    from ..utilities.randtest import _randint

    randint = _randint(seed)

    ct, T = dmp_zz_factor(dmp_LC(f, K), u - 1, K)

    b = dmp_zz_mignotte_bound(f, u, K)
    p = K(nextprime(b))

    if mod is None:
        if u == 1:
            mod = 2
        else:
            mod = 1

    history, configs, A, r = set(), [], [K.zero]*u, None

    try:
        cs, s, E = dmp_zz_wang_test_points(f, T, ct, A, u, K)

        _, H = dup_zz_factor_sqf(s, K)

        r = len(H)

        if r == 1:
            return [f]

        configs = [(s, cs, E, H, A)]
    except EvaluationFailed:
        pass

    eez_num_configs = query('EEZ_NUMBER_OF_CONFIGS')
    eez_num_tries = query('EEZ_NUMBER_OF_TRIES')
    eez_mod_step = query('EEZ_MODULUS_STEP')

    while len(configs) < eez_num_configs:
        for _ in range(eez_num_tries):
            A = [ K(randint(-mod, mod)) for _ in range(u) ]

            if tuple(A) in history:
                continue
            else:
                history.add(tuple(A))

            try:
                cs, s, E = dmp_zz_wang_test_points(f, T, ct, A, u, K)
            except EvaluationFailed:
                continue

            _, H = dup_zz_factor_sqf(s, K)

            rr = len(H)

            if r is not None:
                if rr != r:  # pragma: no cover
                    if rr < r:
                        configs, r = [], rr
                    else:
                        continue
            else:
                r = rr

            if r == 1:
                return [f]

            configs.append((s, cs, E, H, A))

            if len(configs) == eez_num_configs:
                break
        else:
            mod += eez_mod_step

    s_norm, s_arg, i = None, 0, 0

    for s, _, _, _, _ in configs:
        _s_norm = dmp_max_norm(s, 0, K)

        if s_norm is not None:
            if _s_norm < s_norm:
                s_norm = _s_norm
                s_arg = i
        else:
            s_norm = _s_norm

        i += 1

    _, cs, E, H, A = configs[s_arg]
    orig_f = f

    try:
        f, H, LC = dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K)
        factors = dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K)
    except ExtraneousFactors:  # pragma: no cover
        if query('EEZ_RESTART_IF_NEEDED'):
            return dmp_zz_wang(orig_f, u, K, mod + 1)
        else:
            raise ExtraneousFactors(
                "we need to restart algorithm with better parameters")

    result = []

    for f in factors:
        _, f = dmp_ground_primitive(f, u, K)

        if K.is_negative(dmp_ground_LC(f, u, K)):
            f = dmp_neg(f, u, K)

        result.append(f)

    return result


def dmp_zz_factor(f, u, K):
    """
    Factor (non square-free) polynomials in `Z[X]`.

    Given a multivariate polynomial `f` in `Z[x]` computes its complete
    factorization `f_1, ..., f_n` into irreducibles over integers::

                 f = content(f) f_1**k_1 ... f_n**k_n

    The factorization is computed by reducing the input polynomial
    into a primitive square-free polynomial and factoring it using
    Enhanced Extended Zassenhaus (EEZ) algorithm. Trial division
    is used to recover the multiplicities of factors.

    The result is returned as a tuple consisting of::

             (content(f), [(f_1, k_1), ..., (f_n, k_n))

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_zz_factor(2*x**2 - 2*y**2)
    (2, [(x - y, 1), (x + y, 1)])

    References
    ==========

    * [Gathen99]_
    """
    if not u:
        return dup_zz_factor(f, K)

    if dmp_zero_p(f, u):
        return K.zero, []

    cont, g = dmp_ground_primitive(f, u, K)

    if dmp_ground_LC(g, u, K) < 0:
        cont, g = -cont, dmp_neg(g, u, K)

    if dmp_ground_p(g, None, u):
        return cont, []

    G, g = dmp_primitive(g, u, K)

    factors = []

    if dmp_degree_in(g, 0, u) > 0:
        g = dmp_sqf_part(g, u, K)
        H = dmp_zz_wang(g, u, K)
        factors = dmp_trial_division(f, H, u, K)

    for g, k in dmp_zz_factor(G, u - 1, K)[1]:
        factors.insert(0, ([g], k))

    return cont, _sort_factors(factors)


def dmp_ext_factor(f, u, K):
    """Factor multivariate polynomials over algebraic number fields. """
    lc = dmp_ground_LC(f, u, K)
    f = dmp_ground_monic(f, u, K)

    if dmp_ground_p(f, None, u):
        return lc, []

    f, F = dmp_sqf_part(f, u, K), f
    s, g, r = dmp_sqf_norm(f, u, K)

    _, factors = dmp_factor_list(r, u, K.domain)

    if len(factors) == 1:
        factors = [f]
    else:
        H = dmp_raise([K.one, s*K.unit], u, 0, K)

        for i, (factor, _) in enumerate(factors):
            h = dmp_convert(factor, u, K.domain, K)
            h, _, g = dmp_inner_gcd(h, g, u, K)
            for j in range(u + 1):
                h = dmp_swap(h, 0, j, u, K)
                h = dmp_compose(h, H, u, K)
                h = dmp_swap(h, 0, j, u, K)
            factors[i] = h

    return lc, dmp_trial_division(F, factors, u, K)


def dmp_gf_factor(f, u, K):
    """Factor multivariate polynomials over finite fields. """
    if u == 0:
        f = dmp_convert(f, 0, K, K.domain)

        coeff, factors = gf_factor(f, K.mod, K.domain)

        for i, (f, k) in enumerate(factors):
            factors[i] = (dmp_convert(f, 0, K.domain, K), k)

        return K.convert(coeff, K.domain), factors
    else:
        raise NotImplementedError('multivariate polynomials over finite fields')


def dmp_factor_list(f, u, K0):
    """Factor polynomials into irreducibles in `K[X]`. """
    J, f = dmp_terms_gcd(f, u, K0)
    cont, f = dmp_ground_primitive(f, u, K0)

    if K0.is_FiniteField:
        coeff, factors = dmp_gf_factor(f, u, K0)
    elif K0.is_AlgebraicField:
        coeff, factors = dmp_ext_factor(f, u, K0)
    else:
        if not K0.is_Exact:
            K0_inexact, K0 = K0, K0.get_exact()
            f = dmp_convert(f, u, K0_inexact, K0)
        else:
            K0_inexact = None

        if K0.is_Field:
            K = K0.ring

            denom, f = dmp_clear_denoms(f, u, K0, K)
            f = dmp_convert(f, u, K0, K)
        else:
            K = K0

        if K.is_IntegerRing:
            levels, f, v = dmp_exclude(f, u, K)
            coeff, factors = dmp_zz_factor(f, v, K)

            for i, (f, k) in enumerate(factors):
                factors[i] = (dmp_include(f, levels, v, K), k)
        elif K.is_Poly:
            f, v = dmp_inject(f, u, K)

            coeff, factors = dmp_factor_list(f, v, K.domain)

            for i, (f, k) in enumerate(factors):
                factors[i] = (dmp_eject(f, v, K), k)

            coeff = K.convert(coeff, K.domain)
        else:  # pragma: no cover
            raise DomainError('factorization not supported over %s' % K0)

        if K0.is_Field:
            for i, (f, k) in enumerate(factors):
                factors[i] = (dmp_convert(f, u, K, K0), k)

            coeff = K0.convert(coeff, K)
            coeff = K0.quo(coeff, denom)

            if K0_inexact:
                for i, (f, k) in enumerate(factors):
                    f = dmp_quo_ground(f, denom, u, K0)
                    f = dmp_convert(f, u, K0, K0_inexact)
                    factors[i] = (f, k)
                    coeff *= denom**k

                coeff = K0_inexact.convert(coeff, K0)
                K0 = K0_inexact

    for i, j in enumerate(reversed(J)):
        if not j:
            continue

        term = {(0,)*(u - i) + (1,) + (0,)*i: K0.one}
        factors.insert(0, (dmp_from_dict(term, u, K0), j))

    return coeff*cont, _sort_factors(factors)


def dmp_irreducible_p(f, u, K):
    """Returns ``True`` if ``f`` has no factors over its domain. """
    _, factors = dmp_factor_list(f, u, K)

    if not factors:
        return True
    elif len(factors) > 1:
        return False
    else:
        _, k = factors[0]
        return k == 1
