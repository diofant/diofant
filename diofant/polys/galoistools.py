"""Dense univariate polynomials with coefficients in Galois fields."""

import random


def dup_gf_random(n, K, irreducible=False):
    """
    Generate a random polynomial in ``GF(q)[x]`` of degree ``n``.

    Examples
    ========

    >>> dup_gf_random(4, FF(5))  # doctest: +SKIP
    [1 mod 5, 4 mod 5, 4 mod 5, 2 mod 5, 1 mod 5]

    >>> f = dup_gf_random(4, FF(5), irreducible=True)
    >>> f  # doctest: +SKIP
    [1 mod 5, 2 mod 5, 4 mod 5, 4 mod 5, 3 mod 5]
    >>> dup_gf_irreducible_p(f, FF(5))
    True

    """
    while True:
        f = [K.one] + [K(random.randint(0, K.order - 1)) for i in range(n)]
        if not irreducible or dup_gf_irreducible_p(f, K):
            return f


def dup_gf_irreducible_p(f, K):
    """Test irreducibility of a polynomial ``f`` in ``GF(q)[x]``."""
    ring = K.poly_ring('_0')
    f = ring.from_list(f)
    return f.is_irreducible
