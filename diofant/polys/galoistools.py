"""Dense univariate polynomials with coefficients in Galois fields."""

import math
import random

from ..ntheory import factorint
from .densebasic import dmp_degree_in, dmp_one_p, dmp_strip
from .polyconfig import query
from .polyutils import _sort_factors


def dup_gf_pow_mod(f, n, g, K):
    """
    Compute ``f**n`` in ``GF(q)[x]/(g)`` using repeated squaring.

    Given polynomials ``f`` and ``g`` in ``GF(q)[x]`` and a non-negative
    integer ``n``, efficiently computes ``f**n (mod g)`` i.e. the remainder
    of ``f**n`` from division by ``g``, using the repeated squaring algorithm.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (3*x**2 + 2*x + 4).to_dense()
    >>> g = (x + 1).to_dense()
    >>> dup_gf_pow_mod(f, 3, g, R.domain)
    []

    References
    ==========

    * :cite:`Gathen1999modern`, algorithm 4.8

    """
    if not n:
        return [K.one]
    elif n == 1:
        return dmp_rem(f, g, 0, K)
    elif n == 2:
        return dmp_rem(dmp_pow(f, 2, 0, K), g, 0, K)

    h = [K.one]

    while True:
        if n & 1:
            h = dmp_mul(h, f, 0, K)
            h = dmp_rem(h, g, 0, K)
            n -= 1

        n >>= 1

        if not n:
            break

        f = dmp_pow(f, 2, 0, K)
        f = dmp_rem(f, g, 0, K)

    return h


def dup_gf_compose_mod(g, h, f, K):
    """
    Compute polynomial composition ``g(h)`` in ``GF(q)[x]/(f)``.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> g = (3*x**2 + 2*x + 4).to_dense()
    >>> h = (2*x**2 + 2*x + 2).to_dense()
    >>> f = (4*x + 3).to_dense()
    >>> dup_gf_compose_mod(g, h, f, R.domain)
    [4 mod 5]

    """
    if not g:
        return []

    comp = [g[0]]

    for a in g[1:]:
        comp = dmp_mul(comp, h, 0, K)
        comp = dmp_add(comp, [a], 0, K)
        comp = dmp_rem(comp, f, 0, K)

    return comp


def dup_gf_trace_map(a, b, c, n, f, K):
    """
    Compute polynomial trace map in ``GF(q)[x]/(f)``.

    Given a polynomial ``f`` in ``GF(q)[x]``, polynomials ``a``, ``b``,
    ``c`` in the quotient ring ``GF(q)[x]/(f)`` such that ``b = c**t
    (mod f)`` for some positive power ``t`` of ``q``, and a positive
    integer ``n``, returns a mapping::

       a -> a**t**n, a + a**t + a**t**2 + ... + a**t**n (mod f)

    In factorization context, ``b = x**q mod f`` and ``c = x mod f``.
    This way we can efficiently compute trace polynomials in equal
    degree factorization routine, much faster than with other methods,
    like iterated Frobenius algorithm, for large degrees.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> a = (x + 2).to_dense()
    >>> b = (4*x + 4).to_dense()
    >>> c = (x + 1).to_dense()
    >>> f = (3*x**2 + 2*x + 4).to_dense()
    >>> dup_gf_trace_map(a, b, c, 4, f, R.domain)
    ([1 mod 5, 3 mod 5], [1 mod 5, 3 mod 5])

    References
    ==========

    * :cite:`Gathen1992ComputingFM`, algorithm 5.2

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


def dup_gf_random(n, K):
    """
    Generate a random polynomial in ``GF(q)[x]`` of degree ``n``.

    Examples
    ========

    >>> dup_gf_random(4, FF(5))  # doctest: +SKIP
    [1 mod 5, 4 mod 5, 4 mod 5, 2 mod 5, 1 mod 5]

    """
    return [K.one] + [K(random.randint(0, K.order - 1)) for i in range(n)]


def dup_gf_irreducible(n, K):
    """
    Generate random irreducible polynomial of degree ``n`` in ``GF(q)[x]``.

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
    >>> f = (x**10 + 4*x**9 + 2*x**8 + 2*x**7 + 3*x**6 +
    ...      2*x**5 + 4*x**4 + x**3 + 4*x**2 + 4).to_dense()
    >>> dup_gf_irred_p_ben_or(f, R.domain)
    True
    >>> f = (3*x**2 + 2*x + 4).to_dense()
    >>> dup_gf_irred_p_ben_or(f, R.domain)
    False

    References
    ==========

    * :cite:`Ben-Or1981ff`

    """
    n, q = dmp_degree_in(f, 0, 0), K.order

    if n <= 1:
        return True

    x = [K.one, K.zero]
    f = dmp_ground_monic(f, 0, K)

    H = h = dup_gf_pow_mod(x, q, f, K)

    for i in range(n//2):
        g = dmp_sub(h, x, 0, K)

        if dmp_one_p(dmp_gcd(f, g, 0, K), 0, K):
            h = dup_gf_compose_mod(h, H, f, K)
        else:
            return False

    return True


def dup_gf_irred_p_rabin(f, K):
    """
    Rabin's polynomial irreducibility test over finite fields.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (x**10 + 4*x**9 + 2*x**8 + 2*x**7 + 3*x**6 +
    ...      2*x**5 + 4*x**4 + x**3 + 4*x**2 + 4).to_dense()
    >>> dup_gf_irred_p_rabin(f, R.domain)
    True
    >>> f = (3*x**2 + 2*x + 4).to_dense()
    >>> dup_gf_irred_p_rabin(f, R.domain)
    False

    References
    ==========

    * :cite:`Gathen1999modern`, algorithm 14.36

    """
    n, q = dmp_degree_in(f, 0, 0), K.order

    if n <= 1:
        return True

    x = [K.one, K.zero]
    f = dmp_ground_monic(f, 0, K)

    indices = {n//d for d in factorint(n)}

    H = h = dup_gf_pow_mod(x, q, f, K)

    for i in range(1, n):
        if i in indices:
            g = dmp_sub(h, x, 0, K)

            if not dmp_one_p(dmp_gcd(f, g, 0, K), 0, K):
                return False

        h = dup_gf_compose_mod(h, H, f, K)

    return h == x


def dup_gf_irreducible_p(f, K):
    """Test irreducibility of a polynomial ``f`` in ``GF(q)[x]``."""
    ring = K.poly_ring('_0')
    f = ring.from_list(f)
    return f.is_irreducible


def dup_gf_primitive_p(f, K):
    """Test if ``f`` is a primitive polynomial over ``GF(p)``."""
    p = K.characteristic

    assert K.order == p

    if not dup_gf_irreducible_p(f, K):
        return False

    n = dmp_degree_in(f, 0, 0)
    t = [K.one] + [K.zero]*n

    for m in range(n, p**n - 1):
        r = dmp_rem(t, f, 0, K)
        if r == [K.one]:
            return False
        t = dmp_mul(r, [K.one, K.zero], 0, K)
    return True


def dup_gf_Qmatrix(f, K):
    """
    Calculate Berlekamp's ``Q`` matrix.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (3*x**2 + 2*x + 4).to_dense()
    >>> dup_gf_Qmatrix(f, R.domain)
    [[1 mod 5, 0 mod 5],
     [3 mod 5, 4 mod 5]]

    >>> f = (x**4 + 1).to_dense()
    >>> dup_gf_Qmatrix(f, R.domain)
    [[1 mod 5, 0 mod 5, 0 mod 5, 0 mod 5],
     [0 mod 5, 4 mod 5, 0 mod 5, 0 mod 5],
     [0 mod 5, 0 mod 5, 1 mod 5, 0 mod 5],
     [0 mod 5, 0 mod 5, 0 mod 5, 4 mod 5]]

    References
    ==========

    * :cite:`Geddes1992algorithms`, algorithm 8.5

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


def dup_gf_berlekamp(f, K):
    """
    Factor a square-free polynomial over finite fields of small order.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (x**4 + 1).to_dense()
    >>> dup_gf_berlekamp([1, 0, 0, 0, 1], R.domain)
    [[1 mod 5, 0 mod 5, 2 mod 5], [1 mod 5, 0 mod 5, 3 mod 5]]

    References
    ==========

    * :cite:`Geddes1992algorithms`, algorithm 8.4
    * :cite:`Knuth1985seminumerical`, section 4.6.2

    """
    from .solvers import RawMatrix
    Q = dup_gf_Qmatrix(f, K)
    Q = RawMatrix(Q) - RawMatrix.eye(len(Q))
    V = Q.T.nullspace()

    for i, v in enumerate(V):
        V[i] = dmp_strip(list(reversed(v)), 0)

    factors = [f]

    for v in V[1:]:
        for f in list(factors):
            for s in range(K.order):
                h = dmp_sub(v, [K(s)], 0, K)
                g = dmp_gcd(f, h, 0, K)

                if not dmp_one_p(g, 0, K) and g != f:
                    factors.remove(f)

                    f = dmp_quo(f, g, 0, K)
                    factors.extend([f, g])

                if len(factors) == len(V):
                    return _sort_factors(factors, multiple=False)

    return _sort_factors(factors, multiple=False)


def dup_gf_ddf_zassenhaus(f, K):
    """
    Cantor-Zassenhaus: Deterministic Distinct Degree Factorization.

    Given a monic square-free polynomial ``f`` in ``GF(q)[x]``, computes
    partial distinct degree factorization ``f_1 ... f_d`` of ``f`` where
    ``deg(f_i) != deg(f_j)`` for ``i != j``. The result is returned as a
    list of pairs ``(f_i, e_i)`` where ``deg(f_i) > 0`` and ``e_i > 0``
    is an argument to the equal degree factorization routine.

    Examples
    ========

    >>> R, x = ring('x', FF(11))
    >>> f = (x**15 - 1).to_dense()
    >>> dup_gf_ddf_zassenhaus(f, R.domain)
    [([1 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 10 mod 11], 1),
     ([1 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 0 mod 11, 1 mod 11, 0 mod 11,
       0 mod 11, 0 mod 11, 0 mod 11, 1 mod 11], 2)]

    To obtain factorization into irreducibles, use equal degree factorization
    procedure (EDF) with each of the factors.

    References
    ==========

    * :cite:`Gathen1999modern`, algorithm 14.3
    * :cite:`Geddes1992algorithms`, algorithm 8.8

    See Also
    ========

    dup_gf_edf_zassenhaus

    """
    factors, q = [], K.order
    g, x = [[K.one, K.zero]]*2

    for i in range(1, dmp_degree_in(f, 0, 0)//2 + 1):
        g = dup_gf_pow_mod(g, q, f, K)
        h = dmp_gcd(f, dmp_sub(g, x, 0, K), 0, K)

        if not dmp_one_p(h, 0, K):
            factors.append((h, i))

            f = dmp_quo(f, h, 0, K)
            g = dmp_rem(g, f, 0, K)

    if not dmp_one_p(f, 0, K):
        factors += [(f, dmp_degree_in(f, 0, 0))]

    return factors


def dup_gf_edf_zassenhaus(f, n, K):
    """
    Cantor-Zassenhaus: Probabilistic Equal Degree Factorization.

    Given a monic square-free polynomial ``f`` in ``GF(q)[x]`` and
    an integer ``n``, such that ``n`` divides ``deg(f)``, returns all
    irreducible factors ``f_1,...,f_d`` of ``f``, each of degree ``n``.
    EDF procedure gives complete factorization over Galois fields.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (x**3 + x**2 + x + 1).to_dense()
    >>> dup_gf_edf_zassenhaus(f, 1, R.domain)
    [[1 mod 5, 1 mod 5], [1 mod 5, 2 mod 5], [1 mod 5, 3 mod 5]]

    References
    ==========

    * :cite:`Geddes1992algorithms`, algorithm 8.9

    See Also
    ========

    dup_gf_ddf_zassenhaus

    """
    factors = [f]

    if dmp_degree_in(f, 0, 0) <= n:
        return factors

    p, q = K.characteristic, K.order
    N = dmp_degree_in(f, 0, 0) // n

    while len(factors) < N:
        r = dup_gf_random(2*n - 1, K)

        if p == 2:
            h = r

            for i in range(1, n):
                r = dup_gf_pow_mod(r, q, f, K)
                h = dmp_add(h, r, 0, K)
        else:
            h = dup_gf_pow_mod(r, (q**n - 1) // 2, f, K)
            h = dmp_sub(h, [K.one], 0, K)

        g = dmp_gcd(f, h, 0, K)

        if not dmp_one_p(g, 0, K) and g != f:
            factors = (dup_gf_edf_zassenhaus(g, n, K) +
                       dup_gf_edf_zassenhaus(dmp_quo(f, g, 0, K), n, K))

    return _sort_factors(factors, multiple=False)


def dup_gf_ddf_shoup(f, K):
    """
    Kaltofen-Shoup: Deterministic Distinct Degree Factorization.

    Given a monic square-free polynomial ``f`` in ``GF(q)[x]``, computes
    partial distinct degree factorization ``f_1 ... f_d`` of ``f`` where
    ``deg(f_i) != deg(f_j)`` for ``i != j``. The result is returned as a
    list of pairs ``(f_i, e_i)`` where ``deg(f_i) > 0`` and ``e_i > 0``
    is an argument to the equal degree factorization routine.

    Notes
    =====

    This algorithm is an improved version of Zassenhaus algorithm for
    large ``deg(f)`` and order ``q`` (especially for ``deg(f) ~ lg(q)``).

    Examples
    ========

    >>> R, x = ring('x', FF(3))
    >>> f = (x**6 - x**5 + x**4 + x**3 - x).to_dense()
    >>> dup_gf_ddf_shoup(f, R.domain)
    [([1 mod 3, 1 mod 3, 0 mod 3], 1), ([1 mod 3, 1 mod 3, 0 mod 3, 1 mod 3, 2 mod 3], 2)]

    References
    ==========

    * :cite:`Kaltofen1998subquadratic`, algorithm D
    * :cite:`Shoup1995factor`
    * :cite:`Gathen1992frobenious`

    See Also
    ========

    dup_gf_edf_shoup

    """
    n, q = dmp_degree_in(f, 0, 0), K.order
    k = math.ceil(math.sqrt(n//2))
    x = [K.one, K.zero]

    h = dup_gf_pow_mod(x, q, f, K)

    # U[i] = x**(q**i)
    U = [x, h] + [K.zero]*(k - 1)

    for i in range(2, k + 1):
        U[i] = dup_gf_compose_mod(U[i - 1], h, f, K)

    h, U = U[k], U[:k]
    # V[i] = x**(q**(k*(i+1)))
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

            if not dmp_one_p(F, 0, K):
                factors.append((F, k*(i + 1) - j))

            g, j = dmp_quo(g, F, 0, K), j - 1

    if not dmp_one_p(f, 0, K):
        factors.append((f, dmp_degree_in(f, 0, 0)))

    return factors


def dup_gf_edf_shoup(f, n, K):
    """
    Gathen-Shoup: Probabilistic Equal Degree Factorization.

    Given a monic square-free polynomial ``f`` in ``GF(q)[x]`` and
    an integer ``n``, such that ``n`` divides ``deg(f)``, returns all
    irreducible factors ``f_1,...,f_d`` of ``f``, each of degree ``n``.
    EDF procedure gives complete factorization over Galois fields.

    Notes
    =====

    This algorithm is an improved version of Zassenhaus algorithm for
    large ``deg(f)`` and order ``q`` (especially for ``deg(f) ~ lg(q)``).

    Examples
    ========

    >>> R, x = ring('x', FF(2917))
    >>> f = (x**2 + 2837*x + 2277).to_dense()
    >>> dup_gf_edf_shoup(f, 1, R.domain)
    [[1 mod 2917, 852 mod 2917], [1 mod 2917, 1985 mod 2917]]

    References
    ==========

    * :cite:`Shoup1991ffactor`
    * :cite:`Gathen1992ComputingFM`, algorithm 3.6

    See Also
    ========

    dup_gf_ddf_shoup

    """
    q, p = K.order, K.characteristic
    N = dmp_degree_in(f, 0, 0)

    if not N:
        return []
    if N <= n:
        return [f]

    factors, x = [f], [K.one, K.zero]

    r = dup_gf_random(N - 1, K)

    h = dup_gf_pow_mod(x, q, f, K)
    H = dup_gf_trace_map(r, h, x, n - 1, f, K)[1]

    if p == 2:
        h1 = dmp_gcd(f, H, 0, K)
        h2 = dmp_quo(f, h1, 0, K)

        factors = dup_gf_edf_shoup(h1, n, K) + dup_gf_edf_shoup(h2, n, K)
    else:
        h = dup_gf_pow_mod(H, (q - 1)//2, f, K)

        h1 = dmp_gcd(f, h, 0, K)
        h2 = dmp_gcd(f, dmp_sub(h, [K.one], 0, K), 0, K)
        h3 = dmp_quo(f, dmp_mul(h1, h2, 0, K), 0, K)

        factors = (dup_gf_edf_shoup(h1, n, K) + dup_gf_edf_shoup(h2, n, K) +
                   dup_gf_edf_shoup(h3, n, K))

    return _sort_factors(factors, multiple=False)


def dup_gf_zassenhaus(f, K):
    """
    Factor a square-free polynomial over finite fields of medium order.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (x**2 + 4*x + 3).to_dense()
    >>> dup_gf_zassenhaus(f, R.domain)
    [[1 mod 5, 1 mod 5], [1 mod 5, 3 mod 5]]

    """
    factors = []

    for factor, n in dup_gf_ddf_zassenhaus(f, K):
        factors += dup_gf_edf_zassenhaus(factor, n, K)

    return _sort_factors(factors, multiple=False)


def dup_gf_shoup(f, K):
    """
    Factor a square-free polynomial over finite fields of large order.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (x**2 + 4*x + 3).to_dense()
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
    Factor a square-free polynomial ``f`` in ``GF(q)[x]``.

    Returns its complete factorization into irreducibles::

                 f_1(x) f_2(x) ... f_d(x)

    where each ``f_i`` is a monic polynomial and ``gcd(f_i, f_j) == 1``,
    for ``i != j``.  The result is given as a list of factors of ``f``.

    Square-free factors of ``f`` can be factored into irreducibles over
    finite fields using three very different methods:

    Berlekamp
        efficient for very small values of order ``q`` (usually ``q < 25``)
    Cantor-Zassenhaus
        efficient on average input and with "typical" ``q``
    Shoup-Kaltofen-Gathen
        efficient with very large inputs and order

    If you want to use a specific factorization method - set
    ``GF_FACTOR_METHOD`` configuration option with one of ``"berlekamp"``,
    ``"zassenhaus"`` or ``"shoup"`` values.

    Examples
    ========

    >>> R, x = ring('x', FF(5))
    >>> f = (x**2 + 4*x + 3).to_dense()
    >>> dup_gf_factor_sqf(f, R.domain)
    [[1 mod 5, 1 mod 5], [1 mod 5, 3 mod 5]]

    References
    ==========

    * :cite:`Gathen1999modern`, chapter 14

    """
    method = query('GF_FACTOR_METHOD')

    return _factor_methods[method](f, K)


def dmp_ground_monic(f, u, K):
    """Divide all coefficients by ``LC(f)`` in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    f = f.monic()
    return f.to_dense()


def dmp_add(f, g, u, K):
    """Add dense polynomials in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return (f + g).to_dense()


def dmp_sub(f, g, u, K):
    """Subtract dense polynomials in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return (f - g).to_dense()


def dmp_mul(f, g, u, K):
    """Multiply dense polynomials in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return (f*g).to_dense()


def dmp_pow(f, n, u, K):
    """Raise ``f`` to the ``n``-th power in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    return (f**n).to_dense()


def dmp_div(f, g, u, K):
    """Polynomial division with remainder in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return tuple(map(lambda _: _.to_dense(), divmod(f, g)))


def dmp_rem(f, g, u, K):
    """Return polynomial remainder in ``K[X]``."""
    return dmp_div(f, g, u, K)[1]


def dmp_quo(f, g, u, K):
    """Return exact polynomial quotient in ``K[X]``."""
    return dmp_div(f, g, u, K)[0]


def dmp_gcd(f, g, u, K):
    """Computes polynomial GCD of `f` and `g` in `K[X]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return ring.gcd(f, g).to_dense()
