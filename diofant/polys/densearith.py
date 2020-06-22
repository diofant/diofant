"""Arithmetics for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densebasic import dmp_degree_in, dmp_ground, dmp_slice_in, dmp_strip
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

    lo = dup_mul(fl, gl, K)
    hi = dup_mul(fh, gh, K)

    mid = dup_mul(dmp_add(fl, fh, 0, K), dmp_add(gl, gh, 0, K), K)
    mid = dmp_sub(mid, dmp_add(lo, hi, 0, K), 0, K)

    return dmp_add(dmp_add(lo, dup_lshift(mid, n2, K), 0, K),
                   dup_lshift(hi, 2*n2, K), 0, K)


def dup_mul(f, g, K):
    """
    Multiply dense polynomials in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dmp_mul(x - 2, x + 2)
    x**2 - 4

    """
    if f == g:
        return dmp_pow(f, 2, 0, K)

    if not (f and g):
        return []

    df = dmp_degree_in(f, 0, 0)
    dg = dmp_degree_in(g, 0, 0)

    n = max(df, dg) + 1

    if n > query('KARATSUBA_CUTOFF'):
        return dup_mul_karatsuba(f, g, K)

    h = []

    for i in range(df + dg + 1):
        coeff = K.zero

        for j in range(max(0, i - dg), min(df, i) + 1):
            coeff += f[j]*g[i - j]

        h.append(coeff)

    return dmp_strip(h, 0)


def dmp_mul(f, g, u, K):
    """
    Multiply dense polynomials in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_mul(x*y + 1, x)
    x**2*y + x

    """
    if not u:
        return dup_mul(f, g, K)

    if f == g:
        return dmp_pow(f, 2, u, K)

    df = dmp_degree_in(f, 0, u)

    if df < 0:
        return f

    dg = dmp_degree_in(g, 0, u)

    if dg < 0:
        return g

    h, v = [], u - 1

    for i in range(df + dg + 1):
        coeff = dmp_ground(0, v)

        for j in range(max(0, i - dg), min(df, i) + 1):
            coeff = dmp_add(coeff, dmp_mul(f[j], g[i - j], v, K), v, K)

        h.append(coeff)

    return dmp_strip(h, u)


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
