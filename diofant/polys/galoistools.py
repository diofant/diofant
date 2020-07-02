"""Dense univariate polynomials with coefficients in Galois fields."""

import random

from .polyconfig import query


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

    >>> f = dup_gf_irreducible(4, FF(5))
    >>> f  # doctest: +SKIP
    [1 mod 5, 2 mod 5, 4 mod 5, 4 mod 5, 3 mod 5]
    >>> dup_gf_irreducible_p(f, FF(5))
    True

    """
    while True:
        f = dup_gf_random(n, K)
        if dup_gf_irreducible_p(f, K):
            return f


def dup_gf_irreducible_p(f, K):
    """Test irreducibility of a polynomial ``f`` in ``GF(q)[x]``."""
    ring = K.poly_ring('_0')
    f = ring.from_list(f)
    return f.is_irreducible


def dup_gf_berlekamp(f, K):
    """Factor a square-free polynomial over finite fields of small order."""
    ring = K.poly_ring('_0')
    f = ring.from_list(f)
    r = ring._gf_berlekamp(f)
    return [_.to_dense() for _ in r]


def dup_gf_zassenhaus(f, K):
    """Factor a square-free polynomial over finite fields of medium order."""
    ring = K.poly_ring('_0')
    f = ring.from_list(f)
    r = ring._gf_zassenhaus(f)
    return [_.to_dense() for _ in r]


def dup_gf_shoup(f, K):
    """Factor a square-free polynomial over finite fields of large order."""
    ring = K.poly_ring('_0')
    f = ring.from_list(f)
    r = ring._gf_shoup(f)
    return [_.to_dense() for _ in r]


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
