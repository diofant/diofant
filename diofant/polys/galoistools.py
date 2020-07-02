"""Dense univariate polynomials with coefficients in Galois fields."""

import random

from ..core import oo
from .polyconfig import query


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
    """Compute polynomial trace map in ``GF(q)[x]/(f)``."""
    ring = K.poly_ring('_0')
    a, b, c, f = map(ring.from_list, (a, b, c, f))
    return tuple(_.to_dense() for _ in ring._gf_trace_map(a, b, c, n, f))


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


def dup_gf_primitive_p(f, K):
    """Test if ``f`` is a primitive polynomial over ``GF(p)``."""
    p = K.characteristic

    assert K.order == p

    if not dup_gf_irreducible_p(f, K):
        return False

    n = dup_degree(f)
    t = [K.one] + [K.zero]*n

    for m in range(n, p**n - 1):
        r = dmp_rem(t, f, 0, K)
        if r == [K.one]:
            return False
        t = dmp_mul(r, [K.one, K.zero], 0, K)
    return True


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


def dmp_add(f, g, u, K):
    """Add dense polynomials in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return (f + g).to_dense()


def dmp_mul(f, g, u, K):
    """Multiply dense polynomials in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return (f*g).to_dense()


def dmp_div(f, g, u, K):
    """Polynomial division with remainder in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return tuple(map(lambda _: _.to_dense(), divmod(f, g)))


def dmp_rem(f, g, u, K):
    """Return polynomial remainder in ``K[X]``."""
    return dmp_div(f, g, u, K)[1]


def dup_degree(f):
    """Return the leading degree of ``f`` in ``K[x]``."""
    return -oo if not f else len(f) - 1
