"""Arithmetics for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densebasic import (dmp_degree_in, dmp_one, dmp_one_p, dmp_slice_in,
                         dmp_strip, dmp_zero, dmp_zero_p, dmp_zeros)
from .polyconfig import query


def dup_add_term(f, c, i, K):
    """
    Add ``c*x**i`` to ``f`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

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

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_add_term(x*y + 1, 2, 2)
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

    >>> R, x = ring("x", ZZ)

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

    >>> R, x, y = ring("x y", ZZ)

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

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_mul_ground(2*x + 2*y, ZZ(3))
    6*x + 6*y

    """
    if not u:
        return dmp_strip([coeff * c for coeff in f], u)
    else:
        v = u - 1
        return [dmp_mul_ground(coeff, c, v, K) for coeff in f]


def dmp_quo_ground(f, c, u, K):
    """
    Quotient by a constant in ``K[X]``.

    Examples
    ========


    >>> R, x, y = ring("x y", ZZ)
    >>> R.dmp_quo_ground(2*x**2*y + 3*x, ZZ(2))
    x**2*y + x

    >>> R, x, y = ring("x y", QQ)
    >>> R.dmp_quo_ground(2*x**2*y + 3*x, QQ(2))
    x**2*y + 3/2*x

    """
    if not u:
        if not c:
            raise ZeroDivisionError('polynomial division')
        if not f:
            return f

        if K.is_Field:
            return [K.quo(coeff, c) for coeff in f]
        else:
            return [coeff // c for coeff in f]

    v = u - 1
    return [dmp_quo_ground(coeff, c, v, K) for coeff in f]


def dmp_exquo_ground(f, c, u, K):
    """
    Exact quotient by a constant in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", QQ)

    >>> R.dmp_exquo_ground(x**2*y + 2*x, QQ(2))
    1/2*x**2*y + x

    """
    if not u:
        if not c:
            raise ZeroDivisionError('polynomial division')
        if not f:
            return f

        return [K.exquo(coeff, c) for coeff in f]

    v = u - 1
    return [dmp_exquo_ground(coeff, c, v, K) for coeff in f]


def dup_lshift(f, n, K):
    """
    Efficiently multiply ``f`` by ``x**n`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

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

    >>> R, x = ring("x", ZZ)

    >>> R.dup_rshift(x**4 + x**2, 2)
    x**2 + 1
    >>> R.dup_rshift(x**4 + x**2 + 2, 2)
    x**2 + 1

    """
    return f[:-n]


def dmp_abs(f, u, K):
    """
    Make all coefficients positive in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_abs(x**2*y - x)
    x**2*y + x

    """
    if not u:
        return [abs(coeff) for coeff in f]
    else:
        v = u - 1
        return [dmp_abs(coeff, v, K) for coeff in f]


def dmp_neg(f, u, K):
    """
    Negate a polynomial in ``K[X]``.

    Examples
    =======

    >>> R, x, y = ring("x y", ZZ)

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

    >>> R, x = ring("x", ZZ)

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

    >>> R, x, y = ring("x y", ZZ)

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

    >>> R, x = ring("x", ZZ)

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

    >>> R, x, y = ring("x y", ZZ)

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


def dmp_add_mul(f, g, h, u, K):
    """
    Return ``f + g*h`` where ``f, g, h`` are in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_add_mul(x**2 + y, x, x + 2)
    2*x**2 + 2*x + y

    """
    return dmp_add(f, dmp_mul(g, h, u, K), u, K)


def dmp_sub_mul(f, g, h, u, K):
    """
    Return ``f - g*h`` where ``f, g, h`` are in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_sub_mul(x**2 + y, x, x + 2)
    -2*x + y

    """
    return dmp_sub(f, dmp_mul(g, h, u, K), u, K)


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

    >>> R, x = ring("x", ZZ)

    >>> R.dmp_mul(x - 2, x + 2)
    x**2 - 4

    """
    if f == g:
        return dup_sqr(f, K)

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

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_mul(x*y + 1, x)
    x**2*y + x

    """
    if not u:
        return dup_mul(f, g, K)

    if f == g:
        return dmp_sqr(f, u, K)

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


def dup_sqr(f, K):
    """
    Square dense polynomials in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dmp_sqr(x**2 + 1)
    x**4 + 2*x**2 + 1

    """
    df, h = len(f) - 1, []

    for i in range(2*df + 1):
        c = K.zero

        jmin = max(0, i - df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in range(jmin, jmax + 1):
            c += f[j]*f[i - j]

        c += c

        if n & 1:
            elem = f[jmax + 1]
            c += elem**2

        h.append(c)

    return dmp_strip(h, 0)


def dmp_sqr(f, u, K):
    """
    Square dense polynomials in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_sqr(x**2 + x*y + y**2)
    x**4 + 2*x**3*y + 3*x**2*y**2 + 2*x*y**3 + y**4

    """
    if not u:
        return dup_sqr(f, K)

    df = dmp_degree_in(f, 0, u)

    if df < 0:
        return f

    h, v = [], u - 1

    for i in range(2*df + 1):
        c = dmp_zero(v)

        jmin = max(0, i - df)
        jmax = min(i, df)

        n = jmax - jmin + 1

        jmax = jmin + n // 2 - 1

        for j in range(jmin, jmax + 1):
            c = dmp_add(c, dmp_mul(f[j], f[i - j], v, K), v, K)

        c = dmp_mul_ground(c, K(2), v, K)

        if n & 1:
            elem = dmp_sqr(f[jmax + 1], v, K)
            c = dmp_add(c, elem, v, K)

        h.append(c)

    return dmp_strip(h, u)


def dmp_pow(f, n, u, K):
    """
    Raise ``f`` to the ``n``-th power in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_pow(x*y + 1, 3)
    x**3*y**3 + 3*x**2*y**2 + 3*x*y + 1

    """
    if not n:
        return dmp_one(u, K)
    if n < 0:
        raise ValueError("can't raise polynomial to a negative power")
    if n == 1 or dmp_zero_p(f, u) or dmp_one_p(f, u, K):
        return f

    g = dmp_one(u, K)

    while True:
        n, m = n//2, n

        if m & 1:
            g = dmp_mul(g, f, u, K)

            if not n:
                break

        f = dmp_sqr(f, u, K)

    return g


def dmp_div(f, g, u, K):
    """
    Polynomial division with remainder in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)
    >>> R.dmp_div(x**2 + x*y, 2*x + 2)
    (0, x**2 + x*y)

    >>> R, x, y = ring("x y", QQ)
    >>> R.dmp_div(x**2 + x*y, 2*x + 2)
    (1/2*x + 1/2*y - 1/2, -y + 1)

    """
    ring = K.poly_ring(*["_%d" % i for i in range(u + 1)])
    f, g = map(ring.from_dense, (f, g))
    return tuple(map(ring.to_dense, divmod(f, g)))


def dmp_rem(f, g, u, K):
    """
    Return polynomial remainder in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)
    >>> R.dmp_rem(x**2 + x*y, 2*x + 2)
    x**2 + x*y

    >>> R, x, y = ring("x y", QQ)
    >>> R.dmp_rem(x**2 + x*y, 2*x + 2)
    -y + 1

    """
    return dmp_div(f, g, u, K)[1]


def dmp_quo(f, g, u, K):
    """
    Return exact polynomial quotient in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)
    >>> R.dmp_quo(x**2 + x*y, 2*x + 2)
    0

    >>> R, x, y = ring("x y", QQ)
    >>> R.dmp_quo(x**2 + x*y, 2*x + 2)
    1/2*x + 1/2*y - 1/2

    """
    return dmp_div(f, g, u, K)[0]


def dmp_max_norm(f, u, K):
    """
    Return maximum norm of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_max_norm(2*x*y - x - 3)
    3

    """
    if not u:
        return max(dmp_abs(f, 0, K), default=K.zero)

    v = u - 1
    return max(dmp_max_norm(c, v, K) for c in f)


def dmp_l1_norm(f, u, K):
    """
    Return l1 norm of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_l1_norm(2*x*y - x - 3)
    6

    """
    if not u:
        return sum(dmp_abs(f, u, K), K.zero)

    v = u - 1
    return sum(dmp_l1_norm(c, v, K) for c in f)


def dmp_expand(polys, u, K):
    """
    Multiply together several polynomials in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_expand([x**2 + y**2, x + 1])
    x**3 + x**2 + x*y**2 + y**2

    """
    if not polys:
        return dmp_one(u, K)

    f = polys[0]

    for g in polys[1:]:
        f = dmp_mul(f, g, u, K)

    return f
