"""Arithmetics for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densebasic import (dmp_degree_in, dmp_slice_in, dmp_strip, dmp_zero,
                         dmp_zero_p, dmp_zeros)
from .polyconfig import query


def dup_add_term(f, c, i, K):
    """
    Add ``c*x**i`` to ``f`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dmp_add_term(x**2 - 1, ZZ(2), 4)
    2*x**4 + x**2 - 1

    """
    if not c:
        return f

    n = len(f)
    m = n - i - 1

    if i == n - 1:
        return dmp_strip([f[0] + c] + f[1:], 0)
    else:
        if i >= n:
            return [c] + [K.zero]*(i - n) + f
        else:
            return f[:m] + [f[m] + c] + f[m + 1:]


def dmp_add_term(f, c, i, u, K):
    """
    Add ``c(x_2..x_u)*x_0**i`` to ``f`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_add_term(x*y + 1, R(2), 2)
    2*x**2 + x*y + 1

    """
    if not u:
        return dup_add_term(f, c, i, K)

    v = u - 1

    if dmp_zero_p(c, v):
        return f

    n = len(f)
    m = n - i - 1

    if i == n - 1:
        return dmp_strip([dmp_add(f[0], c, v, K)] + f[1:], u)
    else:
        if i >= n:
            return [c] + dmp_zeros(i - n, v, K) + f
        else:
            return f[:m] + [dmp_add(f[m], c, v, K)] + f[m + 1:]


def dup_mul_term(f, c, i, K):
    """
    Multiply ``f`` by ``c*x**i`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dmp_mul_term(x**2 - 1, ZZ(3), 2)
    3*x**4 - 3*x**2

    """
    if not c or not f:
        return []
    else:
        return [cf * c for cf in f] + [K.zero]*i


def dmp_mul_term(f, c, i, u, K):
    """
    Multiply ``f`` by ``c(x_2..x_u)*x_0**i`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_mul_term(x**2*y + x, 3*y, 2)
    3*x**4*y**2 + 3*x**3*y

    """
    if not u:
        return dup_mul_term(f, c, i, K)

    v = u - 1

    if dmp_zero_p(f, u):
        return f
    if dmp_zero_p(c, v):
        return dmp_zero(u)
    else:
        return [dmp_mul(cf, c, v, K) for cf in f] + dmp_zeros(i, v, K)


def dmp_mul_ground(f, c, u, K):
    """
    Multiply ``f`` by a constant value in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_mul_ground(2*x + 2*y, ZZ(3))
    6*x + 6*y

    """
    if not u:
        return dmp_strip([coeff * c for coeff in f], u)
    else:
        v = u - 1
        return [dmp_mul_ground(coeff, c, v, K) for coeff in f]


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
    """
    Negate a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_neg(x**2*y - x)
    -x**2*y + x

    """
    if not u:
        return [-coeff for coeff in f]
    else:
        v = u - 1
        return [dmp_neg(coeff, v, K) for coeff in f]


def dup_add(f, g, K):
    """
    Add dense polynomials in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dmp_add(x**2 - 1, x - 2)
    x**2 + x - 3

    """
    if not f:
        return g
    if not g:
        return f

    df = dmp_degree_in(f, 0, 0)
    dg = dmp_degree_in(g, 0, 0)

    if df == dg:
        return dmp_strip([a + b for a, b in zip(f, g)], 0)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [a + b for a, b in zip(f, g)]


def dmp_add(f, g, u, K):
    """
    Add dense polynomials in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_add(x**2 + y, x**2*y + x)
    x**2*y + x**2 + x + y

    """
    if not u:
        return dup_add(f, g, K)

    df = dmp_degree_in(f, 0, u)

    if df < 0:
        return g

    dg = dmp_degree_in(g, 0, u)

    if dg < 0:
        return f

    v = u - 1

    if df == dg:
        return dmp_strip([dmp_add(a, b, v, K) for a, b in zip(f, g)], u)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = g[:k], g[k:]

        return h + [dmp_add(a, b, v, K) for a, b in zip(f, g)]


def dup_sub(f, g, K):
    """
    Subtract dense polynomials in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dmp_sub(x**2 - 1, x - 2)
    x**2 - x + 1

    """
    if not f:
        return dmp_neg(g, 0, K)
    if not g:
        return f

    df = dmp_degree_in(f, 0, 0)
    dg = dmp_degree_in(g, 0, 0)

    if df == dg:
        return dmp_strip([a - b for a, b in zip(f, g)], 0)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = dmp_neg(g[:k], 0, K), g[k:]

        return h + [a - b for a, b in zip(f, g)]


def dmp_sub(f, g, u, K):
    """
    Subtract dense polynomials in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_sub(x**2 + y, x**2*y + x)
    -x**2*y + x**2 - x + y

    """
    if not u:
        return dup_sub(f, g, K)

    df = dmp_degree_in(f, 0, u)

    if df < 0:
        return dmp_neg(g, u, K)

    dg = dmp_degree_in(g, 0, u)

    if dg < 0:
        return f

    v = u - 1

    if df == dg:
        return dmp_strip([dmp_sub(a, b, v, K) for a, b in zip(f, g)], u)
    else:
        k = abs(df - dg)

        if df > dg:
            h, f = f[:k], f[k:]
        else:
            h, g = dmp_neg(g[:k], u, K), g[k:]

        return h + [dmp_sub(a, b, v, K) for a, b in zip(f, g)]


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

    mid = dup_mul(dup_add(fl, fh, K), dup_add(gl, gh, K), K)
    mid = dup_sub(mid, dup_add(lo, hi, K), K)

    return dup_add(dup_add(lo, dup_lshift(mid, n2, K), K),
                   dup_lshift(hi, 2*n2, K), K)


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
        coeff = dmp_zero(v)

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
