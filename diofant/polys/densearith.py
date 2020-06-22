"""Arithmetics for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densebasic import dmp_degree_in, dmp_slice_in
from .polyconfig import query


def dmp_mul_ground(f, c, u, K):
    """Multiply ``f`` by a constant value in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    r = f.mul_ground(c)
    return r.to_dense()


def dup_lshift(f, n, K):
    """
    Efficiently multiply ``f`` by ``x**n`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dup_lshift(x**2 + 1, 2)
    x**4 + x**2

    """
    if not f:
        return f
    else:
        return f + [K.zero]*n


def dup_rshift(f, n, K):
    """
    Efficiently divide ``f`` by ``x**n`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dup_rshift(x**4 + x**2, 2)
    x**2 + 1
    >>> R.dup_rshift(x**4 + x**2 + 2, 2)
    x**2 + 1

    """
    return f[:-n]


def dmp_neg(f, u, K):
    """Negate a polynomial in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    return (-f).to_dense()


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


def dup_mul_karatsuba(f, g, K):
    """
    Multiply dense polynomials in ``K[x]`` using Karatsuba's algorithm.

    References
    ==========

    * :cite:`Hoeven02`

    """
    df = dmp_degree_in(f, 0, 0)
    dg = dmp_degree_in(g, 0, 0)

    n = max(df, dg) + 1

    n2 = n//2

    fl = dmp_slice_in(f, 0, n2, 0, 0, K)
    gl = dmp_slice_in(g, 0, n2, 0, 0, K)

    fh = dup_rshift(dmp_slice_in(f, n2, n, 0, 0, K), n2, K)
    gh = dup_rshift(dmp_slice_in(g, n2, n, 0, 0, K), n2, K)

    lo = dmp_mul(fl, gl, 0, K)
    hi = dmp_mul(fh, gh, 0, K)

    mid = dmp_mul(dmp_add(fl, fh, 0, K), dmp_add(gl, gh, 0, K), 0, K)
    mid = dmp_sub(mid, dmp_add(lo, hi, 0, K), 0, K)

    return dmp_add(dmp_add(lo, dup_lshift(mid, n2, K), 0, K),
                   dup_lshift(hi, 2*n2, K), 0, K)


def dmp_mul(f, g, u, K):
    """Multiply dense polynomials in ``K[X]``."""
    df = dmp_degree_in(f, 0, 0)
    dg = dmp_degree_in(g, 0, 0)

    if not u and max(df, dg) + 1 > query('KARATSUBA_CUTOFF'):
        return dup_mul_karatsuba(f, g, K)

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
