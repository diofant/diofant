"""Dense univariate polynomials with coefficients in Galois fields. """

import math
import random

from ..ntheory import factorint
from .densearith import (dmp_add, dmp_add_term, dmp_mul, dmp_mul_ground,
                         dmp_sqr, dmp_sub, dup_lshift)
from .densebasic import (dmp_convert, dmp_degree_in, dmp_from_dict, dmp_normal,
                         dmp_strip)
from .polyconfig import query
from .polyutils import _sort_factors


def gf_from_dict(f, p, K):
    """
    Create a ``GF(p)[x]`` polynomial from a dict.

    Examples
    ========

    >>> gf_from_dict({10: 4, 4: 33, 0: -1}, 5, ZZ)
    [4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 4]

    """
    f = dmp_from_dict(f, 0, K.finite_field(p))
    f = dmp_normal(f, 0, K.finite_field(p))

    return dmp_convert(f, 0, K.finite_field(p), K)


def gf_trunc(f, p):
    """
    Reduce all coefficients modulo ``p``.

    Examples
    ========

    >>> gf_trunc([7, -2, 3], 5)
    [2, 3, 3]

    """
    return dmp_strip([a % p for a in f], 0)


def gf_add_ground(f, a, p, K):
    """
    Compute ``f + a`` where ``f`` in ``GF(p)[x]`` and ``a`` in ``GF(p)``.

    Examples
    ========

    >>> gf_add_ground([3, 2, 4], 2, 5, ZZ)
    [3, 2, 1]

    """
    return gf_trunc(dmp_add_term(f, a, 0, 0, K), p)


def gf_sub_ground(f, a, p, K):
    """
    Compute ``f - a`` where ``f`` in ``GF(p)[x]`` and ``a`` in ``GF(p)``.

    Examples
    ========

    >>> gf_sub_ground([3, 2, 4], 2, 5, ZZ)
    [3, 2, 2]

    """
    return gf_trunc(dmp_add_term(f, -a, 0, 0, K), p)


def gf_mul_ground(f, a, p, K):
    """
    Compute ``f * a`` where ``f`` in ``GF(p)[x]`` and ``a`` in ``GF(p)``.

    Examples
    ========

    >>> gf_mul_ground([3, 2, 4], 2, 5, ZZ)
    [1, 4, 3]

    """
    return gf_trunc(dmp_mul_ground(f, a, 0, K), p)


def gf_quo_ground(f, a, p, K):
    """
    Compute ``f/a`` where ``f`` in ``GF(p)[x]`` and ``a`` in ``GF(p)``.

    Examples
    ========

    >>> gf_quo_ground([3, 2, 4], 2, 5, ZZ)
    [4, 1, 2]

    """
    return gf_mul_ground(f, K.invert(a, p), p, K)


def gf_add(f, g, p, K):
    """
    Add polynomials in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_add([3, 2, 4], [2, 2, 2], 5, ZZ)
    [4, 1]

    """
    return gf_trunc(dmp_add(f, g, 0, K), p)


def gf_sub(f, g, p, K):
    """
    Subtract polynomials in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_sub([3, 2, 4], [2, 2, 2], 5, ZZ)
    [1, 0, 2]

    """
    return gf_trunc(dmp_sub(f, g, 0, K), p)


def gf_mul(f, g, p, K):
    """
    Multiply polynomials in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_mul([3, 2, 4], [2, 2, 2], 5, ZZ)
    [1, 0, 3, 2, 3]

    """
    return gf_trunc(dmp_mul(f, g, 0, K), p)


def gf_sqr(f, p, K):
    """
    Square polynomials in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_sqr([3, 2, 4], 5, ZZ)
    [4, 2, 3, 1, 1]

    """
    return gf_trunc(dmp_sqr(f, 0, K), p)


def gf_div(f, g, p, K):
    """
    Division with remainder in ``GF(p)[x]``.

    Given univariate polynomials ``f`` and ``g`` with coefficients in a
    finite field with ``p`` elements, returns polynomials ``q`` and ``r``
    (quotient and remainder) such that ``f = q*g + r``.

    Examples
    ========

    >>> gf_div([1, 0, 1, 1], [1, 1, 0], 2, ZZ)
    ([1, 1], [1])

    References
    ==========

    * :cite:`Monagan1993inplace`
    * :cite:`Gathen1999modern`

    """
    df = dmp_degree_in(f, 0, 0)
    dg = dmp_degree_in(g, 0, 0)

    if not g:
        raise ZeroDivisionError("polynomial division")
    elif df < dg:
        return [], f

    inv = K.invert(g[0], p)

    h, dq, dr = list(f), df - dg, dg - 1

    for i in range(df + 1):
        coeff = h[i]

        for j in range(max(0, dg - i), min(df - i, dr) + 1):
            coeff -= h[i + j - dg] * g[dg - j]

        if i <= dq:
            coeff *= inv

        h[i] = coeff % p

    return h[:dq + 1], dmp_strip(h[dq + 1:], 0)


def gf_rem(f, g, p, K):
    """
    Compute polynomial remainder in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_rem([1, 0, 1, 1], [1, 1, 0], 2, ZZ)
    [1]

    """
    return gf_div(f, g, p, K)[1]


def gf_quo(f, g, p, K):
    """
    Compute exact quotient in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_quo([1, 0, 1, 1], [1, 1, 0], 2, ZZ)
    [1, 1]
    >>> gf_quo([1, 0, 3, 2, 3], [2, 2, 2], 5, ZZ)
    [3, 2, 4]

    """
    return gf_div(f, g, p, K)[0]


def gf_frobenius_monomial_base(g, p, K):
    """
    return the list of ``x**(i*p) mod g in Z_p`` for ``i = 0, .., n - 1``
    where ``n = dmp_degree_in(g, 0, 0)``

    Examples
    ========

    >>> gf_frobenius_monomial_base([1, 0, 2, 1], 5, ZZ)
    [[1], [4, 4, 2], [1, 2]]

    """
    n = dmp_degree_in(g, 0, 0)
    if n == 0:
        return []
    b = [0]*n
    b[0] = [1]
    if p < n:
        for i in range(1, n):
            mon = dup_lshift(b[i - 1], p, K)
            b[i] = gf_rem(mon, g, p, K)
    elif n > 1:
        b[1] = gf_pow_mod([K.one, K.zero], p, g, p, K)
        for i in range(2, n):
            b[i] = gf_mul(b[i - 1], b[1], p, K)
            b[i] = gf_rem(b[i], g, p, K)

    return b


def gf_frobenius_map(f, g, b, p, K):
    """
    compute gf_pow_mod(f, p, g, p, K) using the Frobenius map

    Parameters
    ==========

    f, g : polynomials in ``GF(p)[x]``
    b : frobenius monomial base
    p : prime number
    K : domain

    Examples
    ========

    >>> f = [2, 1, 0, 1]
    >>> g = [1, 0, 2, 1]
    >>> p = 5
    >>> b = gf_frobenius_monomial_base(g, p, ZZ)
    >>> r = gf_frobenius_map(f, g, b, p, ZZ)
    >>> gf_frobenius_map(f, g, b, p, ZZ)
    [4, 0, 3]

    """
    m = dmp_degree_in(g, 0, 0)
    if dmp_degree_in(f, 0, 0) >= m:
        f = gf_rem(f, g, p, K)
    if not f:
        return []
    n = dmp_degree_in(f, 0, 0)
    sf = [f[-1]]
    for i in range(1, n + 1):
        v = gf_mul_ground(b[i], f[n - i], p, K)
        sf = gf_add(sf, v, p, K)
    return sf


def _gf_pow_pnm1d2(f, n, g, b, p, K):
    """
    utility function for ``gf_edf_zassenhaus``
    Compute ``f**((p**n - 1) // 2)`` in ``GF(p)[x]/(g)``
    ``f**((p**n - 1) // 2) = (f*f**p*...*f**(p**n - 1))**((p - 1) // 2)``

    """
    f = gf_rem(f, g, p, K)
    h = f
    r = f
    for i in range(1, n):
        h = gf_frobenius_map(h, g, b, p, K)
        r = gf_mul(r, h, p, K)
        r = gf_rem(r, g, p, K)

    res = gf_pow_mod(r, (p - 1)//2, g, p, K)
    return res


def gf_pow_mod(f, n, g, p, K):
    """
    Compute ``f**n`` in ``GF(p)[x]/(g)`` using repeated squaring.

    Given polynomials ``f`` and ``g`` in ``GF(p)[x]`` and a non-negative
    integer ``n``, efficiently computes ``f**n (mod g)`` i.e. the remainder
    of ``f**n`` from division by ``g``, using the repeated squaring algorithm.

    Examples
    ========

    >>> gf_pow_mod([3, 2, 4], 3, [1, 1], 5, ZZ)
    []

    References
    ==========

    * :cite:`Gathen1999modern`

    """
    if not n:
        return [K.one]
    elif n == 1:
        return gf_rem(f, g, p, K)
    elif n == 2:
        return gf_rem(gf_sqr(f, p, K), g, p, K)

    h = [K.one]

    while True:
        if n & 1:
            h = gf_mul(h, f, p, K)
            h = gf_rem(h, g, p, K)
            n -= 1

        n >>= 1

        if not n:
            break

        f = gf_sqr(f, p, K)
        f = gf_rem(f, g, p, K)

    return h


def gf_gcd(f, g, p, K):
    """
    Euclidean Algorithm in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_gcd([3, 2, 4], [2, 2, 3], 5, ZZ)
    [1, 3]

    """
    while g:
        f, g = g, gf_rem(f, g, p, K)

    return gf_monic(f, p, K)[1]


def gf_monic(f, p, K):
    """
    Compute LC and a monic polynomial in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_monic([3, 2, 4], 5, ZZ)
    (3, [1, 4, 3])

    """
    if not f:
        return K.zero, []
    else:
        lc = f[0]

        if lc == K.one:
            return lc, list(f)
        else:
            return lc, gf_quo_ground(f, lc, p, K)


def gf_compose_mod(g, h, f, p, K):
    """
    Compute polynomial composition ``g(h)`` in ``GF(p)[x]/(f)``.

    Examples
    ========

    >>> gf_compose_mod([3, 2, 4], [2, 2, 2], [4, 3], 5, ZZ)
    [4]

    """
    if not g:
        return []

    comp = [g[0]]

    for a in g[1:]:
        comp = gf_mul(comp, h, p, K)
        comp = gf_add_ground(comp, a, p, K)
        comp = gf_rem(comp, f, p, K)

    return comp


def gf_trace_map(a, b, c, n, f, p, K):
    """
    Compute polynomial trace map in ``GF(p)[x]/(f)``.

    Given a polynomial ``f`` in ``GF(p)[x]``, polynomials ``a``, ``b``,
    ``c`` in the quotient ring ``GF(p)[x]/(f)`` such that ``b = c**t
    (mod f)`` for some positive power ``t`` of ``p``, and a positive
    integer ``n``, returns a mapping::

       a -> a**t**n, a + a**t + a**t**2 + ... + a**t**n (mod f)

    In factorization context, ``b = x**p mod f`` and ``c = x mod f``.
    This way we can efficiently compute trace polynomials in equal
    degree factorization routine, much faster than with other methods,
    like iterated Frobenius algorithm, for large degrees.

    Examples
    ========

    >>> gf_trace_map([1, 2], [4, 4], [1, 1], 4, [3, 2, 4], 5, ZZ)
    ([1, 3], [1, 3])

    References
    ==========

    * :cite:`Gathen1992frobenious`

    """
    u = gf_compose_mod(a, b, f, p, K)
    v = b

    if n & 1:
        U = gf_add(a, u, p, K)
        V = b
    else:
        U = a
        V = c

    n >>= 1

    while n:
        u = gf_add(u, gf_compose_mod(u, v, f, p, K), p, K)
        v = gf_compose_mod(v, v, f, p, K)

        if n & 1:
            U = gf_add(U, gf_compose_mod(u, V, f, p, K), p, K)
            V = gf_compose_mod(v, V, f, p, K)

        n >>= 1

    return gf_compose_mod(a, V, f, p, K), U


def _gf_trace_map(f, n, g, b, p, K):
    """
    utility for ``gf_edf_shoup``

    """
    f = gf_rem(f, g, p, K)
    h = f
    r = f
    for i in range(1, n):
        h = gf_frobenius_map(h, g, b, p, K)
        r = gf_add(r, h, p, K)
        r = gf_rem(r, g, p, K)
    return r


def gf_random(n, p, K):
    """
    Generate a random polynomial in ``GF(p)[x]`` of degree ``n``.

    Examples
    ========

    >>> gf_random(10, 5, ZZ) #doctest: +SKIP
    [1, 2, 3, 2, 1, 1, 1, 2, 0, 4, 2]

    """
    return [K.one] + [K(int(random.uniform(0, p))) for i in range(n)]


def gf_irreducible(n, p, K):
    """
    Generate random irreducible polynomial of degree ``n`` in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_irreducible(10, 5, ZZ) #doctest: +SKIP
    [1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4]

    """
    while True:
        f = gf_random(n, p, K)
        if gf_irreducible_p(f, p, K):
            return f


def gf_irred_p_ben_or(f, p, K):
    """
    Ben-Or's polynomial irreducibility test over finite fields.

    Examples
    ========

    >>> gf_irred_p_ben_or([1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4], 5, ZZ)
    True
    >>> gf_irred_p_ben_or([3, 2, 4], 5, ZZ)
    False

    References
    ==========

    * :cite:`Ben-Or1981ff`

    """
    n = dmp_degree_in(f, 0, 0)

    if n <= 1:
        return True

    _, f = gf_monic(f, p, K)
    if n < 5:
        H = h = gf_pow_mod([K.one, K.zero], p, f, p, K)

        for i in range(n//2):
            g = gf_sub(h, [K.one, K.zero], p, K)

            if gf_gcd(f, g, p, K) == [K.one]:
                h = gf_compose_mod(h, H, f, p, K)
            else:
                return False
    else:
        b = gf_frobenius_monomial_base(f, p, K)
        H = h = gf_frobenius_map([K.one, K.zero], f, b, p, K)
        for i in range(n//2):
            g = gf_sub(h, [K.one, K.zero], p, K)
            if gf_gcd(f, g, p, K) == [K.one]:
                h = gf_frobenius_map(h, f, b, p, K)
            else:
                return False

    return True


def gf_irred_p_rabin(f, p, K):
    """
    Rabin's polynomial irreducibility test over finite fields.

    Examples
    ========

    >>> gf_irred_p_rabin([1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4], 5, ZZ)
    True
    >>> gf_irred_p_rabin([3, 2, 4], 5, ZZ)
    False

    """
    n = dmp_degree_in(f, 0, 0)

    if n <= 1:
        return True

    _, f = gf_monic(f, p, K)

    x = [K.one, K.zero]

    indices = {n//d for d in factorint(n)}

    b = gf_frobenius_monomial_base(f, p, K)
    h = b[1]

    for i in range(1, n):
        if i in indices:
            g = gf_sub(h, x, p, K)

            if gf_gcd(f, g, p, K) != [K.one]:
                return False

        h = gf_frobenius_map(h, f, b, p, K)

    return h == x


_irred_methods = {
    'ben-or': gf_irred_p_ben_or,
    'rabin': gf_irred_p_rabin,
}


def gf_irreducible_p(f, p, K):
    """
    Test irreducibility of a polynomial ``f`` in ``GF(p)[x]``.

    Examples
    ========

    >>> gf_irreducible_p([1, 4, 2, 2, 3, 2, 4, 1, 4, 0, 4], 5, ZZ)
    True
    >>> gf_irreducible_p([3, 2, 4], 5, ZZ)
    False

    """
    method = query('GF_IRRED_METHOD')

    return _irred_methods[method](f, p, K)


def gf_Qmatrix(f, p, K):
    """
    Calculate Berlekamp's ``Q`` matrix.

    Examples
    ========

    >>> gf_Qmatrix([3, 2, 4], 5, ZZ)
    [[1, 0],
     [3, 4]]

    >>> gf_Qmatrix([1, 0, 0, 0, 1], 5, ZZ)
    [[1, 0, 0, 0],
     [0, 4, 0, 0],
     [0, 0, 1, 0],
     [0, 0, 0, 4]]

    """
    n, r = dmp_degree_in(f, 0, 0), int(p)

    q = [K.one] + [K.zero]*(n - 1)
    Q = [list(q)] + [[]]*(n - 1)

    for i in range(1, (n - 1)*r + 1):
        qq, c = [(-q[-1]*f[-1]) % p], q[-1]

        for j in range(1, n):
            qq.append((q[j - 1] - c*f[-j - 1]) % p)

        if not (i % r):
            Q[i//r] = list(qq)

        q = qq

    return Q


def gf_Qbasis(Q, p, K):
    """
    Compute a basis of the kernel of ``Q``.

    Examples
    ========

    >>> gf_Qbasis(gf_Qmatrix([1, 0, 0, 0, 1], 5, ZZ), 5, ZZ)
    [[1, 0, 0, 0], [0, 0, 1, 0]]

    >>> gf_Qbasis(gf_Qmatrix([3, 2, 4], 5, ZZ), 5, ZZ)
    [[1, 0]]

    """
    Q, n = [list(q) for q in Q], len(Q)

    for k in range(n):
        Q[k][k] = (Q[k][k] - K.one) % p

    for k in range(n):
        for i in range(k, n):
            if Q[k][i]:
                break
        else:
            continue

        inv = K.invert(Q[k][i], p)

        for j in range(n):
            Q[j][i] = (Q[j][i]*inv) % p

        for j in range(n):
            t = Q[j][k]
            Q[j][k] = Q[j][i]
            Q[j][i] = t

        for i in range(n):
            if i != k:
                q = Q[k][i]

                for j in range(n):
                    Q[j][i] = (Q[j][i] - Q[j][k]*q) % p

    for i in range(n):
        for j in range(n):
            if i == j:
                Q[i][j] = (K.one - Q[i][j]) % p
            else:
                Q[i][j] = (-Q[i][j]) % p

    basis = []

    for q in Q:
        if any(q):
            basis.append(q)

    return basis


def gf_berlekamp(f, p, K):
    """
    Factor a square-free ``f`` in ``GF(p)[x]`` for small ``p``.

    Examples
    ========

    >>> gf_berlekamp([1, 0, 0, 0, 1], 5, ZZ)
    [[1, 0, 2], [1, 0, 3]]

    """
    Q = gf_Qmatrix(f, p, K)
    V = gf_Qbasis(Q, p, K)

    for i, v in enumerate(V):
        V[i] = dmp_strip(list(reversed(v)), 0)

    factors = [f]

    for k in range(1, len(V)):
        for f in list(factors):
            s = K.zero

            while s < p:
                g = gf_sub_ground(V[k], s, p, K)
                h = gf_gcd(f, g, p, K)

                if h != [K.one] and h != f:
                    factors.remove(f)

                    f = gf_quo(f, h, p, K)
                    factors.extend([f, h])

                if len(factors) == len(V):
                    return _sort_factors(factors, multiple=False)

                s += K.one

    return _sort_factors(factors, multiple=False)


def gf_ddf_zassenhaus(f, p, K):
    """
    Cantor-Zassenhaus: Deterministic Distinct Degree Factorization

    Given a monic square-free polynomial ``f`` in ``GF(p)[x]``, computes
    partial distinct degree factorization ``f_1 ... f_d`` of ``f`` where
    ``deg(f_i) != deg(f_j)`` for ``i != j``. The result is returned as a
    list of pairs ``(f_i, e_i)`` where ``deg(f_i) > 0`` and ``e_i > 0``
    is an argument to the equal degree factorization routine.

    Examples
    ========

    >>> f = gf_from_dict({15: ZZ(1), 0: ZZ(-1)}, 11, ZZ)

    >>> gf_ddf_zassenhaus(f, 11, ZZ)
    [([1, 0, 0, 0, 0, 10], 1), ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], 2)]

    To obtain factorization into irreducibles, use equal degree factorization
    procedure (EDF) with each of the factors.

    References
    ==========

    * :cite:`Gathen1999modern`
    * :cite:`Geddes1992algorithms`

    """
    i, g, factors = 1, [K.one, K.zero], []

    b = gf_frobenius_monomial_base(f, p, K)
    while 2*i <= dmp_degree_in(f, 0, 0):
        g = gf_frobenius_map(g, f, b, p, K)
        h = gf_gcd(f, gf_sub(g, [K.one, K.zero], p, K), p, K)

        if h != [K.one]:
            factors.append((h, i))

            f = gf_quo(f, h, p, K)
            g = gf_rem(g, f, p, K)
            b = gf_frobenius_monomial_base(f, p, K)

        i += 1

    if f != [K.one]:
        return factors + [(f, dmp_degree_in(f, 0, 0))]
    else:
        return factors


def gf_edf_zassenhaus(f, n, p, K):
    """
    Cantor-Zassenhaus: Probabilistic Equal Degree Factorization

    Given a monic square-free polynomial ``f`` in ``GF(p)[x]`` and
    an integer ``n``, such that ``n`` divides ``deg(f)``, returns all
    irreducible factors ``f_1,...,f_d`` of ``f``, each of degree ``n``.
    EDF procedure gives complete factorization over Galois fields.

    Examples
    ========

    >>> gf_edf_zassenhaus([1, 1, 1, 1], 1, 5, ZZ)
    [[1, 1], [1, 2], [1, 3]]

    References
    ==========

    * :cite:`Gathen1999modern`
    * :cite:`Geddes1992algorithms`

    """
    factors = [f]

    if dmp_degree_in(f, 0, 0) <= n:
        return factors

    N = dmp_degree_in(f, 0, 0) // n
    if p != 2:
        b = gf_frobenius_monomial_base(f, p, K)

    while len(factors) < N:
        r = gf_random(2*n - 1, p, K)

        if p == 2:
            h = r

            for i in range(2**(n*N - 1)):
                r = gf_pow_mod(r, 2, f, p, K)
                h = gf_add(h, r, p, K)

            g = gf_gcd(f, h, p, K)
        else:
            h = _gf_pow_pnm1d2(r, n, f, b, p, K)
            g = gf_gcd(f, gf_sub_ground(h, K.one, p, K), p, K)

        if g != [K.one] and g != f:
            factors = gf_edf_zassenhaus(g, n, p, K) \
                + gf_edf_zassenhaus(gf_quo(f, g, p, K), n, p, K)

    return _sort_factors(factors, multiple=False)


def gf_ddf_shoup(f, p, K):
    """
    Kaltofen-Shoup: Deterministic Distinct Degree Factorization

    Given a monic square-free polynomial ``f`` in ``GF(p)[x]``, computes
    partial distinct degree factorization ``f_1,...,f_d`` of ``f`` where
    ``deg(f_i) != deg(f_j)`` for ``i != j``. The result is returned as a
    list of pairs ``(f_i, e_i)`` where ``deg(f_i) > 0`` and ``e_i > 0``
    is an argument to the equal degree factorization routine.

    This algorithm is an improved version of Zassenhaus algorithm for
    large ``deg(f)`` and modulus ``p`` (especially for ``deg(f) ~ lg(p)``).

    Examples
    ========

    >>> f = gf_from_dict({6: ZZ(1), 5: ZZ(-1), 4: ZZ(1), 3: ZZ(1), 1: ZZ(-1)}, 3, ZZ)

    >>> gf_ddf_shoup(f, 3, ZZ)
    [([1, 1, 0], 1), ([1, 1, 0, 1, 2], 2)]

    References
    ==========

    * :cite:`Kaltofen1998subquadratic`
    * :cite:`Shoup1995factor`
    * :cite:`Gathen1992frobenious`

    """
    n = dmp_degree_in(f, 0, 0)
    k = math.ceil(math.sqrt(n//2))
    b = gf_frobenius_monomial_base(f, p, K)
    h = gf_frobenius_map([K.one, K.zero], f, b, p, K)
    # U[i] = x**(p**i)
    U = [[K.one, K.zero], h] + [K.zero]*(k - 1)

    for i in range(2, k + 1):
        U[i] = gf_frobenius_map(U[i-1], f, b, p, K)

    h, U = U[k], U[:k]
    # V[i] = x**(p**(k*(i+1)))
    V = [h] + [K.zero]*(k - 1)

    for i in range(1, k):
        V[i] = gf_compose_mod(V[i - 1], h, f, p, K)

    factors = []

    for i, v in enumerate(V):
        h, j = [K.one], k - 1

        for u in U:
            g = gf_sub(v, u, p, K)
            h = gf_mul(h, g, p, K)
            h = gf_rem(h, f, p, K)

        g = gf_gcd(f, h, p, K)
        f = gf_quo(f, g, p, K)

        for u in reversed(U):
            h = gf_sub(v, u, p, K)
            F = gf_gcd(g, h, p, K)

            if F != [K.one]:
                factors.append((F, k*(i + 1) - j))

            g, j = gf_quo(g, F, p, K), j - 1

    if f != [K.one]:
        factors.append((f, dmp_degree_in(f, 0, 0)))

    return factors


def gf_edf_shoup(f, n, p, K):
    """
    Gathen-Shoup: Probabilistic Equal Degree Factorization

    Given a monic square-free polynomial ``f`` in ``GF(p)[x]`` and integer
    ``n`` such that ``n`` divides ``deg(f)``, returns all irreducible factors
    ``f_1,...,f_d`` of ``f``, each of degree ``n``. This is a complete
    factorization over Galois fields.

    This algorithm is an improved version of Zassenhaus algorithm for
    large ``deg(f)`` and modulus ``p`` (especially for ``deg(f) ~ lg(p)``).

    Examples
    ========

    >>> gf_edf_shoup([1, 2837, 2277], 1, 2917, ZZ)
    [[1, 852], [1, 1985]]

    References
    ==========

    * :cite:`Shoup1991ffactor`
    * :cite:`Gathen1992frobenious`

    """
    N, q = dmp_degree_in(f, 0, 0), int(p)

    if not N:
        return []
    if N <= n:
        return [f]

    factors, x = [f], [K.one, K.zero]

    r = gf_random(N - 1, p, K)

    if p == 2:
        h = gf_pow_mod(x, q, f, p, K)
        H = gf_trace_map(r, h, x, n - 1, f, p, K)[1]
        h1 = gf_gcd(f, H, p, K)
        h2 = gf_quo(f, h1, p, K)

        factors = gf_edf_shoup(h1, n, p, K) \
            + gf_edf_shoup(h2, n, p, K)
    else:
        b = gf_frobenius_monomial_base(f, p, K)
        H = _gf_trace_map(r, n, f, b, p, K)
        h = gf_pow_mod(H, (q - 1)//2, f, p, K)

        h1 = gf_gcd(f, h, p, K)
        h2 = gf_gcd(f, gf_sub_ground(h, K.one, p, K), p, K)
        h3 = gf_quo(f, gf_mul(h1, h2, p, K), p, K)

        factors = gf_edf_shoup(h1, n, p, K) \
            + gf_edf_shoup(h2, n, p, K) \
            + gf_edf_shoup(h3, n, p, K)

    return _sort_factors(factors, multiple=False)


def gf_zassenhaus(f, p, K):
    """
    Factor a square-free ``f`` in ``GF(p)[x]`` for medium ``p``.

    Examples
    ========

    >>> gf_zassenhaus([1, 4, 3], 5, ZZ)
    [[1, 1], [1, 3]]

    """
    factors = []

    for factor, n in gf_ddf_zassenhaus(f, p, K):
        factors += gf_edf_zassenhaus(factor, n, p, K)

    return _sort_factors(factors, multiple=False)


def gf_shoup(f, p, K):
    """
    Factor a square-free ``f`` in ``GF(p)[x]`` for large ``p``.

    Examples
    ========

    >>> gf_shoup([1, 4, 3], 5, ZZ)
    [[1, 1], [1, 3]]

    """
    factors = []

    for factor, n in gf_ddf_shoup(f, p, K):
        factors += gf_edf_shoup(factor, n, p, K)

    return _sort_factors(factors, multiple=False)


_factor_methods = {
    'berlekamp': gf_berlekamp,  # ``p`` : small
    'zassenhaus': gf_zassenhaus,  # ``p`` : medium
    'shoup': gf_shoup,      # ``p`` : large
}


def gf_factor_sqf(f, p, K):
    """
    Factor a square-free polynomial ``f`` in ``GF(p)[x]``.

    Returns its complete factorization into irreducibles::

                 f_1(x) f_2(x) ... f_d(x)

    where each ``f_i`` is a monic polynomial and ``gcd(f_i, f_j) == 1``,
    for ``i != j``.  The result is given as a tuple consisting of the
    leading coefficient of ``f`` and a list of factors of ``f``.

    Square-free factors of ``f`` can be factored into irreducibles over
    ``GF(p)`` using three very different methods:

    Berlekamp
        efficient for very small values of ``p`` (usually ``p < 25``)
    Cantor-Zassenhaus
        efficient on average input and with "typical" ``p``
    Shoup-Kaltofen-Gathen
        efficient with very large inputs and modulus

    If you want to use a specific factorization method, instead of the default
    one, set ``GF_FACTOR_METHOD`` with one of ``berlekamp``, ``zassenhaus`` or
    ``shoup`` values.

    Examples
    ========

    >>> gf_factor_sqf([3, 2, 4], 5, ZZ)
    (3, [[1, 1], [1, 3]])

    References
    ==========

    * :cite:`Gathen1999modern`

    """
    lc, f = gf_monic(f, p, K)

    if dmp_degree_in(f, 0, 0) < 1:
        return lc, []

    method = query('GF_FACTOR_METHOD')

    return lc, _factor_methods[method](f, p, K)
