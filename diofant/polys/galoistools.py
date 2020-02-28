"""Dense univariate polynomials with coefficients in Galois fields."""

import math
import random

from ..ntheory import factorint
from .densearith import (dmp_add, dmp_add_term, dmp_mul, dmp_mul_ground,
                         dmp_quo, dmp_rem, dmp_sqr, dmp_sub, dup_lshift)
from .densebasic import dmp_degree_in, dmp_ground_LC, dmp_strip
from .densetools import dmp_ground_monic
from .euclidtools import dmp_gcd
from .polyconfig import query
from .polyutils import _sort_factors


def dup_gf_frobenius_monomial_base(g, K):
    """
    Return the list of ``x**(i*p) mod g in Z_p`` for ``i = 0, .., n - 1``
    where ``n = dmp_degree_in(g, 0, 0)``.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**3 + 2*x + 1)
    >>> dup_gf_frobenius_monomial_base(f, R.domain)
    [[1 mod 5], [4 mod 5, 4 mod 5, 2 mod 5], [1 mod 5, 2 mod 5]]

    """
    n = dmp_degree_in(g, 0, 0)
    if n == 0:
        return []
    b = [0]*n
    b[0] = [K.one]
    p = K.mod
    if p < n:
        for i in range(1, n):
            mon = dup_lshift(b[i - 1], p, K)
            b[i] = dmp_rem(mon, g, 0, K)
    elif n > 1:
        b[1] = dup_gf_pow_mod([K.one, K.zero], p, g, K)
        for i in range(2, n):
            b[i] = dmp_mul(b[i - 1], b[1], 0, K)
            b[i] = dmp_rem(b[i], g, 0, K)

    return b


def dup_gf_frobenius_map(f, g, b, K):
    """
    Compute dup_gf_pow_mod(f, p, g, K) using the Frobenius map.

    Parameters
    ==========

    f, g : polynomials in ``GF(p)[x]``
    b : frobenius monomial base
    K : domain

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(2*x**3 + x**2 + 1)
    >>> g = R.to_dense(x**3 + 2*x + 1)
    >>> b = dup_gf_frobenius_monomial_base(g, R.domain)
    >>> dup_gf_frobenius_map(f, g, b, R.domain)
    [4 mod 5, 0 mod 5, 3 mod 5]

    """
    m = dmp_degree_in(g, 0, 0)
    if dmp_degree_in(f, 0, 0) >= m:
        f = dmp_rem(f, g, 0, K)
    if not f:
        return []
    n = dmp_degree_in(f, 0, 0)
    sf = [f[-1]]
    for i in range(1, n + 1):
        v = dmp_mul_ground(b[i], f[n - i], 0, K)
        sf = dmp_add(sf, v, 0, K)
    return sf


def _dup_gf_pow_pnm1d2(f, n, g, b, K):
    """
    Utility function for ``gf_edf_zassenhaus``.

    Compute ``f**((p**n - 1) // 2)`` in ``GF(p)[x]/(g)``
    ``f**((p**n - 1) // 2) = (f*f**p*...*f**(p**n - 1))**((p - 1) // 2)``

    """
    f = dmp_rem(f, g, 0, K)
    h = f
    r = f
    for i in range(1, n):
        h = dup_gf_frobenius_map(h, g, b, K)
        r = dmp_mul(r, h, 0, K)
        r = dmp_rem(r, g, 0, K)

    res = dup_gf_pow_mod(r, (K.mod - 1)//2, g, K)
    return res


def dup_gf_pow_mod(f, n, g, K):
    """
    Compute ``f**n`` in ``GF(p)[x]/(g)`` using repeated squaring.

    Given polynomials ``f`` and ``g`` in ``GF(p)[x]`` and a non-negative
    integer ``n``, efficiently computes ``f**n (mod g)`` i.e. the remainder
    of ``f**n`` from division by ``g``, using the repeated squaring algorithm.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> g = R.to_dense(x + 1)
    >>> dup_gf_pow_mod(f, 3, g, R.domain)
    []

    References
    ==========

    * :cite:`Gathen1999modern`

    """
    if not n:
        return [K.one]
    elif n == 1:
        return dmp_rem(f, g, 0, K)
    elif n == 2:
        return dmp_rem(dmp_sqr(f, 0, K), g, 0, K)

    h = [K.one]

    while True:
        if n & 1:
            h = dmp_mul(h, f, 0, K)
            h = dmp_rem(h, g, 0, K)
            n -= 1

        n >>= 1

        if not n:
            break

        f = dmp_sqr(f, 0, K)
        f = dmp_rem(f, g, 0, K)

    return h


def dup_gf_compose_mod(g, h, f, K):
    """
    Compute polynomial composition ``g(h)`` in ``GF(p)[x]/(f)``.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> g = R.to_dense(3*x**2 + 2*x + 4)
    >>> h = R.to_dense(2*x**2 + 2*x + 2)
    >>> f = R.to_dense(4*x + 3)
    >>> dup_gf_compose_mod(g, h, f, R.domain)
    [4 mod 5]

    """
    if not g:
        return []

    comp = [g[0]]

    for a in g[1:]:
        comp = dmp_mul(comp, h, 0, K)
        comp = dmp_add_term(comp, a, 0, 0, K)
        comp = dmp_rem(comp, f, 0, K)

    return comp


def dup_gf_trace_map(a, b, c, n, f, K):
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

    >>> R, x = ring('x', FF(5))
    >>> a = R.to_dense(x + 2)
    >>> b = R.to_dense(4*x + 4)
    >>> c = R.to_dense(x + 1)
    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> dup_gf_trace_map(a, b, c, 4, f, R.domain)
    ([1 mod 5, 3 mod 5], [1 mod 5, 3 mod 5])

    References
    ==========

    * :cite:`Gathen1992frobenious`

    """
    u = dup_gf_compose_mod(a, b, f, K)
    v = b

    if n & 1:
        U = dmp_add(a, u, 0, K)
        V = b
    else:
        U = a
        V = c

    n >>= 1

    while n:
        u = dmp_add(u, dup_gf_compose_mod(u, v, f, K), 0, K)
        v = dup_gf_compose_mod(v, v, f, K)

        if n & 1:
            U = dmp_add(U, dup_gf_compose_mod(u, V, f, K), 0, K)
            V = dup_gf_compose_mod(v, V, f, K)

        n >>= 1

    return dup_gf_compose_mod(a, V, f, K), U


def _gf_trace_map(f, n, g, b, K):
    """
    Utility for ``gf_edf_shoup``.

    """
    f = dmp_rem(f, g, 0, K)
    h = f
    r = f
    for i in range(1, n):
        h = dup_gf_frobenius_map(h, g, b, K)
        r = dmp_add(r, h, 0, K)
        r = dmp_rem(r, g, 0, K)
    return r


def dup_gf_random(n, K):
    """
    Generate a random polynomial in ``GF(p)[x]`` of degree ``n``.

    Examples
    ========

    >>> dup_gf_random(4, FF(5))  # doctest: +SKIP
    [1 mod 5, 4 mod 5, 4 mod 5, 2 mod 5, 1 mod 5]

    """
    return [K.one] + [K(int(random.uniform(0, K.order))) for i in range(n)]


def dup_gf_irreducible(n, K):
    """
    Generate random irreducible polynomial of degree ``n`` in ``GF(p)[x]``.

    Examples
    ========

    >>> dup_gf_irreducible(4, FF(5))  # doctest: +SKIP
    [1 mod 5, 2 mod 5, 4 mod 5, 4 mod 5, 3 mod 5]
    >>> dup_gf_irreducible_p(_, FF(5))
    True

    """
    while True:
        f = dup_gf_random(n, K)
        if dup_gf_irreducible_p(f, K):
            return f


def dup_gf_irred_p_ben_or(f, K):
    """
    Ben-Or's polynomial irreducibility test over finite fields.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**10 + 4*x**9 + 2*x**8 + 2*x**7 + 3*x**6 +
    ...                2*x**5 + 4*x**4 + x**3 + 4*x**2 + 4)
    >>> dup_gf_irred_p_ben_or(f, R.domain)
    True
    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> dup_gf_irred_p_ben_or(f, R.domain)
    False

    References
    ==========

    * :cite:`Ben-Or1981ff`

    """
    n = dmp_degree_in(f, 0, 0)

    if n <= 1:
        return True

    f = dmp_ground_monic(f, 0, K)
    p = K.mod
    if n < 5:
        H = h = dup_gf_pow_mod([K.one, K.zero], p, f, K)

        for i in range(n//2):
            g = dmp_sub(h, [K.one, K.zero], 0, K)

            if dmp_gcd(f, g, 0, K) == [K.one]:
                h = dup_gf_compose_mod(h, H, f, K)
            else:
                return False
    else:
        b = dup_gf_frobenius_monomial_base(f, K)
        H = h = dup_gf_frobenius_map([K.one, K.zero], f, b, K)
        for i in range(n//2):
            g = dmp_sub(h, [K.one, K.zero], 0, K)
            if dmp_gcd(f, g, 0, K) == [K.one]:
                h = dup_gf_frobenius_map(h, f, b, K)
            else:
                return False

    return True


def dup_gf_irred_p_rabin(f, K):
    """
    Rabin's polynomial irreducibility test over finite fields.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**10 + 4*x**9 + 2*x**8 + 2*x**7 + 3*x**6 +
    ...                2*x**5 + 4*x**4 + x**3 + 4*x**2 + 4)
    >>> dup_gf_irred_p_rabin(f, R.domain)
    True
    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> dup_gf_irred_p_rabin(f, R.domain)
    False

    """
    n = dmp_degree_in(f, 0, 0)

    if n <= 1:
        return True

    f = dmp_ground_monic(f, 0, K)

    x = [K.one, K.zero]

    indices = {n//d for d in factorint(n)}

    b = dup_gf_frobenius_monomial_base(f, K)
    h = b[1]

    for i in range(1, n):
        if i in indices:
            g = dmp_sub(h, x, 0, K)

            if dmp_gcd(f, g, 0, K) != [K.one]:
                return False

        h = dup_gf_frobenius_map(h, f, b, K)

    return h == x


_irred_methods = {
    'ben-or': dup_gf_irred_p_ben_or,
    'rabin': dup_gf_irred_p_rabin,
}


def dup_gf_irreducible_p(f, K):
    """
    Test irreducibility of a polynomial ``f`` in ``GF(p)[x]``.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**10 + 4*x**9 + 2*x**8 + 2*x**7 + 3*x**6 +
    ...                2*x**5 + 4*x**4 + x**3 + 4*x**2 + 4)
    >>> dup_gf_irreducible_p(f, R.domain)
    True
    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> dup_gf_irreducible_p(f, R.domain)
    False

    """
    method = query('GF_IRRED_METHOD')

    return _irred_methods[method](f, K)


def dup_gf_Qmatrix(f, K):
    """
    Calculate Berlekamp's ``Q`` matrix.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> dup_gf_Qmatrix(f, R.domain)
    [[1 mod 5, 0 mod 5],
     [3 mod 5, 4 mod 5]]

    >>> f = R.to_dense(x**4 + 1)
    >>> dup_gf_Qmatrix(f, R.domain)
    [[1 mod 5, 0 mod 5, 0 mod 5, 0 mod 5],
     [0 mod 5, 4 mod 5, 0 mod 5, 0 mod 5],
     [0 mod 5, 0 mod 5, 1 mod 5, 0 mod 5],
     [0 mod 5, 0 mod 5, 0 mod 5, 4 mod 5]]

    References
    ==========

    * :cite:`Geddes1992algorithms`, algorithm 8.5.

    """
    n, q = dmp_degree_in(f, 0, 0), K.order

    r = [K.one] + [K.zero]*(n - 1)
    Q = [r.copy()] + [[]]*(n - 1)

    for i in range(1, (n - 1)*q + 1):
        c, r[1:], r[0] = r[-1], r[:-1], K.zero
        for j in range(n):
            r[j] -= c*f[-j - 1]

        if not (i % q):
            Q[i//q] = r.copy()

    return Q


def dup_gf_Qbasis(Q, K):
    """
    Compute a basis of the kernel of ``Q - I``.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**4 + 1)
    >>> dup_gf_Qbasis(dup_gf_Qmatrix(f, R.domain), R.domain)
    [[1 mod 5, 0 mod 5, 0 mod 5, 0 mod 5],
     [0 mod 5, 0 mod 5, 1 mod 5, 0 mod 5]]

    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> dup_gf_Qbasis(dup_gf_Qmatrix(f, R.domain), R.domain)
    [[1 mod 5, 0 mod 5]]

    References
    ==========

    * :cite:`Knuth1985seminumerical`, section 4.6.2, algorithm N.

    """
    Q, n = [list(q) for q in Q], len(Q)
    for k in range(n):
        Q[k][k] -= K.one

    c, basis = [-1]*n, []

    for k in range(n):
        for j in range(n):
            if Q[k][j] and c[j] < 0:
                c[j] = k

                q = -K.one/Q[k][j]
                for i in range(k, n):
                    Q[i][j] *= q

                for i in range(n):
                    q = Q[k][i]
                    if i != j:
                        for s in range(k, n):
                            Q[s][i] += q*Q[s][j]

                break
        else:
            v = [K.zero]*n
            v[k] = K.one
            for j in range(n):
                for s in range(n):
                    if c[s] == j:
                        v[j] = Q[k][s]
            basis.append(v)

    return basis


def dup_gf_berlekamp(f, K):
    """
    Factor a square-free ``f`` in ``GF(p^m)[x]`` for small order.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**4 + 1)
    >>> dup_gf_berlekamp([1, 0, 0, 0, 1], R.domain)
    [[1 mod 5, 0 mod 5, 2 mod 5], [1 mod 5, 0 mod 5, 3 mod 5]]

    References
    ==========

    * :cite:`Geddes1992algorithms`, algorithm 8.4.
    * :cite:`Knuth1985seminumerical`, section 4.6.2.

    """
    Q = dup_gf_Qmatrix(f, K)
    V = dup_gf_Qbasis(Q, K)

    for i, v in enumerate(V):
        V[i] = dmp_strip(list(reversed(v)), 0)

    factors = [f]

    for v in V[1:]:
        for f in list(factors):
            for s in range(K.order):
                h = dmp_add_term(v, -K(s), 0, 0, K)
                g = dmp_gcd(f, h, 0, K)

                if g != [K.one] and g != f:
                    factors.remove(f)

                    f = dmp_quo(f, g, 0, K)
                    factors.extend([f, g])

                if len(factors) == len(V):
                    return _sort_factors(factors, multiple=False)

    return _sort_factors(factors, multiple=False)


def dup_gf_ddf_zassenhaus(f, K):
    """
    Cantor-Zassenhaus: Deterministic Distinct Degree Factorization

    Given a monic square-free polynomial ``f`` in ``GF(p)[x]``, computes
    partial distinct degree factorization ``f_1 ... f_d`` of ``f`` where
    ``deg(f_i) != deg(f_j)`` for ``i != j``. The result is returned as a
    list of pairs ``(f_i, e_i)`` where ``deg(f_i) > 0`` and ``e_i > 0``
    is an argument to the equal degree factorization routine.

    Examples
    ========

    >>> R, x = ring('x', FF(11))
    >>> f = R.to_dense(x**15 - 1)
    >>> dup_gf_ddf_zassenhaus(f, R.domain)
    [([1 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 10 mod 11], 1),
     ([1 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 1 mod 11, 0 mod 11,
       0 mod 11, 0 mod 11, 0 mod 11, 1 mod 11], 2)]

    To obtain factorization into irreducibles, use equal degree factorization
    procedure (EDF) with each of the factors.

    References
    ==========

    * :cite:`Gathen1999modern`
    * :cite:`Geddes1992algorithms`

    """
    i, g, factors = 1, [K.one, K.zero], []

    b = dup_gf_frobenius_monomial_base(f, K)
    while 2*i <= dmp_degree_in(f, 0, 0):
        g = dup_gf_frobenius_map(g, f, b, K)
        h = dmp_gcd(f, dmp_sub(g, [K.one, K.zero], 0, K), 0, K)

        if h != [K.one]:
            factors.append((h, i))

            f = dmp_quo(f, h, 0, K)
            g = dmp_rem(g, f, 0, K)
            b = dup_gf_frobenius_monomial_base(f, K)

        i += 1

    if f != [K.one]:
        return factors + [(f, dmp_degree_in(f, 0, 0))]
    else:
        return factors


def dup_gf_edf_zassenhaus(f, n, K):
    """
    Cantor-Zassenhaus: Probabilistic Equal Degree Factorization

    Given a monic square-free polynomial ``f`` in ``GF(p)[x]`` and
    an integer ``n``, such that ``n`` divides ``deg(f)``, returns all
    irreducible factors ``f_1,...,f_d`` of ``f``, each of degree ``n``.
    EDF procedure gives complete factorization over Galois fields.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**3 + x**2 + x + 1)
    >>> dup_gf_edf_zassenhaus(f, 1, R.domain)
    [[1 mod 5, 1 mod 5], [1 mod 5, 2 mod 5], [1 mod 5, 3 mod 5]]

    References
    ==========

    * :cite:`Gathen1999modern`
    * :cite:`Geddes1992algorithms`

    """
    factors = [f]

    if dmp_degree_in(f, 0, 0) <= n:
        return factors

    N = dmp_degree_in(f, 0, 0) // n
    p = K.mod
    if p != 2:
        b = dup_gf_frobenius_monomial_base(f, K)

    while len(factors) < N:
        r = dup_gf_random(2*n - 1, K)

        if p == 2:
            h = r

            for i in range(2**(n*N - 1)):
                r = dup_gf_pow_mod(r, 2, f, K)
                h = dmp_add(h, r, 0, K)

            g = dmp_gcd(f, h, 0, K)
        else:
            h = _dup_gf_pow_pnm1d2(r, n, f, b, K)
            g = dmp_gcd(f, dmp_add_term(h, -K.one, 0, 0, K), 0, K)

        if g != [K.one] and g != f:
            factors = dup_gf_edf_zassenhaus(g, n, K) \
                + dup_gf_edf_zassenhaus(dmp_quo(f, g, 0, K), n, K)

    return _sort_factors(factors, multiple=False)


def dup_gf_ddf_shoup(f, K):
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

    >>> R, x = ring('x', FF(3))
    >>> f = R.to_dense(x**6 - x**5 + x**4 + x**3 - x)
    >>> dup_gf_ddf_shoup(f, R.domain)
    [([1 mod 3, 1 mod 3, 0 mod 3], 1), ([1 mod 3, 1 mod 3, 0 mod 3, 1 mod 3, 2 mod 3], 2)]

    References
    ==========

    * :cite:`Kaltofen1998subquadratic`
    * :cite:`Shoup1995factor`
    * :cite:`Gathen1992frobenious`

    """
    n = dmp_degree_in(f, 0, 0)
    k = math.ceil(math.sqrt(n//2))
    b = dup_gf_frobenius_monomial_base(f, K)
    h = dup_gf_frobenius_map([K.one, K.zero], f, b, K)
    # U[i] = x**(p**i)
    U = [[K.one, K.zero], h] + [K.zero]*(k - 1)

    for i in range(2, k + 1):
        U[i] = dup_gf_frobenius_map(U[i-1], f, b, K)

    h, U = U[k], U[:k]
    # V[i] = x**(p**(k*(i+1)))
    V = [h] + [K.zero]*(k - 1)

    for i in range(1, k):
        V[i] = dup_gf_compose_mod(V[i - 1], h, f, K)

    factors = []

    for i, v in enumerate(V):
        h, j = [K.one], k - 1

        for u in U:
            g = dmp_sub(v, u, 0, K)
            h = dmp_mul(h, g, 0, K)
            h = dmp_rem(h, f, 0, K)

        g = dmp_gcd(f, h, 0, K)
        f = dmp_quo(f, g, 0, K)

        for u in reversed(U):
            h = dmp_sub(v, u, 0, K)
            F = dmp_gcd(g, h, 0, K)

            if F != [K.one]:
                factors.append((F, k*(i + 1) - j))

            g, j = dmp_quo(g, F, 0, K), j - 1

    if f != [K.one]:
        factors.append((f, dmp_degree_in(f, 0, 0)))

    return factors


def dup_gf_edf_shoup(f, n, K):
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

    >>> R, x = ring('x', FF(2917))
    >>> f = R.to_dense(x**2 + 2837*x + 2277)
    >>> dup_gf_edf_shoup(f, 1, R.domain)
    [[1 mod 2917, 852 mod 2917], [1 mod 2917, 1985 mod 2917]]

    References
    ==========

    * :cite:`Shoup1991ffactor`
    * :cite:`Gathen1992frobenious`

    """
    p = K.mod
    N, q = dmp_degree_in(f, 0, 0), int(p)

    if not N:
        return []
    if N <= n:
        return [f]

    factors, x = [f], [K.one, K.zero]

    r = dup_gf_random(N - 1, K)

    if p == 2:
        h = dup_gf_pow_mod(x, q, f, K)
        H = dup_gf_trace_map(r, h, x, n - 1, f, K)[1]
        h1 = dmp_gcd(f, H, 0, K)
        h2 = dmp_quo(f, h1, 0, K)

        factors = dup_gf_edf_shoup(h1, n, K) + dup_gf_edf_shoup(h2, n, K)
    else:
        b = dup_gf_frobenius_monomial_base(f, K)
        H = _gf_trace_map(r, n, f, b, K)
        h = dup_gf_pow_mod(H, (q - 1)//2, f, K)

        h1 = dmp_gcd(f, h, 0, K)
        h2 = dmp_gcd(f, dmp_add_term(h, -K.one, 0, 0, K), 0, K)
        h3 = dmp_quo(f, dmp_mul(h1, h2, 0, K), 0, K)

        factors = dup_gf_edf_shoup(h1, n, K) \
            + dup_gf_edf_shoup(h2, n, K) \
            + dup_gf_edf_shoup(h3, n, K)

    return _sort_factors(factors, multiple=False)


def dup_gf_zassenhaus(f, K):
    """
    Factor a square-free ``f`` in ``GF(p)[x]`` for medium ``p``.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**2 + 4*x + 3)
    >>> dup_gf_zassenhaus(f, R.domain)
    [[1 mod 5, 1 mod 5], [1 mod 5, 3 mod 5]]

    """
    factors = []

    for factor, n in dup_gf_ddf_zassenhaus(f, K):
        factors += dup_gf_edf_zassenhaus(factor, n, K)

    return _sort_factors(factors, multiple=False)


def dup_gf_shoup(f, K):
    """
    Factor a square-free ``f`` in ``GF(p)[x]`` for large ``p``.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(x**2 + 4*x + 3)
    >>> dup_gf_shoup(f, R.domain)
    [[1 mod 5, 1 mod 5], [1 mod 5, 3 mod 5]]

    """
    factors = []

    for factor, n in dup_gf_ddf_shoup(f, K):
        factors += dup_gf_edf_shoup(factor, n, K)

    return _sort_factors(factors, multiple=False)


_factor_methods = {
    'berlekamp': dup_gf_berlekamp,  # ``p`` : small
    'zassenhaus': dup_gf_zassenhaus,  # ``p`` : medium
    'shoup': dup_gf_shoup,      # ``p`` : large
}


def dup_gf_factor_sqf(f, K):
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

    >>> R, x = ring('x', FF(5))
    >>> f = R.to_dense(3*x**2 + 2*x + 4)
    >>> dup_gf_factor_sqf(f, R.domain)
    (3 mod 5, [[1 mod 5, 1 mod 5], [1 mod 5, 3 mod 5]])

    References
    ==========

    * :cite:`Gathen1999modern`

    """
    lc = dmp_ground_LC(f, 0, K)
    f = dmp_ground_monic(f, 0, K)
    method = query('GF_FACTOR_METHOD')

    return lc, _factor_methods[method](f, K)
