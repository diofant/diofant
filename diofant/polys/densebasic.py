"""Basic tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

import random

from ..core import oo
from .monomials import Monomial


def dmp_LC(f, K):
    """
    Return leading coefficient of ``f``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dmp_LC(x**2 + 2*x + 3)
    1

    """
    if not f:
        return K.zero
    else:
        return f[0]


def dmp_TC(f, K):
    """
    Return trailing coefficient of ``f``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dmp_TC(x**2 + 2*x + 3)
    3

    """
    if not f:
        return K.zero
    else:
        return f[-1]


def dmp_ground_TC(f, u, K):
    """
    Return the ground trailing coefficient.

    Examples
    ========

    >>> R, x, y, z = ring('x y z', ZZ)

    >>> R.dmp_ground_TC(y + 2*z + 3)
    3

    """
    while u:
        f = dmp_TC(f, K)
        u -= 1

    return dmp_TC(f, K)


def dmp_degree_in(f, j, u):
    """
    Return the leading degree of ``f`` in ``x_j`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_degree_in(2*x + y**2 + 2*y + 3, 1)
    2

    """
    if not j:
        return -oo if dmp_zero_p(f, u) else len(f) - 1

    if j < 0 or j > u:
        raise IndexError(f'0 <= j <= {u} expected, got {j}')

    def degree_in(g, v, i, j):
        if i == j:
            return dmp_degree_in(g, 0, v)

        v, i = v - 1, i + 1

        return max(degree_in(c, v, i, j) for c in g)

    return degree_in(f, u, 0, j)


def dmp_strip(f, u):
    """
    Remove leading zeros from ``f`` in ``K[X]``.

    Examples
    ========

    >>> dmp_strip([[], [ZZ(0), ZZ(1), ZZ(2)], [ZZ(1)]], 1)
    [[0, 1, 2], [1]]

    """
    if not u:
        for i, c in enumerate(f):
            if c:
                return f[i:]
        return dmp_zero(u)

    v = u - 1

    for i, c in enumerate(f):
        if not dmp_zero_p(c, v):
            return f[i:]
    return dmp_zero(u)


def dup_reverse(f):
    """
    Compute ``x**n * f(1/x)``, i.e.: reverse ``f`` in ``K[x]``.

    Examples
    ========

    >>> dup_reverse([ZZ(1), ZZ(2), ZZ(3), ZZ(0)])
    [3, 2, 1]

    """
    return dmp_strip(list(reversed(f)), 0)


def dmp_to_tuple(f, u):
    """
    Convert ``f`` into a nested :class:`tuple`.

    This is needed for hashing.

    Examples
    ========

    >>> dmp_to_tuple([ZZ(1), ZZ(2), ZZ(3), ZZ(0)], 0)
    (1, 2, 3, 0)

    >>> dmp_to_tuple([[ZZ(1)], [ZZ(1), ZZ(2)]], 1)
    ((1,), (1, 2))

    """
    if not u:
        return tuple(f)

    v = u - 1
    return tuple(dmp_to_tuple(c, v) for c in f)


def dmp_convert(f, u, K0, K1):
    """
    Convert the ground domain of ``f`` from ``K0`` to ``K1``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> dmp_convert([[R(1)], [R(2)]], 1, R, ZZ)
    [[1], [2]]
    >>> dmp_convert([[ZZ(1)], [ZZ(2)]], 1, ZZ, R)
    [[1], [2]]

    """
    if K0 is not None and K0 == K1:
        return f

    if not u:
        r = [K1.convert(c, K0) for c in f]
    else:
        v = u - 1
        r = [dmp_convert(c, v, K0, K1) for c in f]

    return dmp_strip(r, u)


def dmp_zero_p(f, u):
    """
    Return ``True`` if ``f`` is zero in ``K[X]``.

    Examples
    ========

    >>> dmp_zero_p([[[[[]]]]], 4)
    True
    >>> dmp_zero_p([[[[[ZZ(1)]]]]], 4)
    False

    """
    while u:
        if len(f) != 1:
            return False

        f = f[0]
        u -= 1

    return not f


def dmp_zero(u):
    """
    Return a multivariate zero.

    Examples
    ========

    >>> dmp_zero(4)
    [[[[[]]]]]

    """
    r = []

    for i in range(u):
        r = [r]

    return r


def dmp_one_p(f, u, K):
    """
    Return ``True`` if ``f`` is one in ``K[X]``.

    Examples
    ========

    >>> dmp_one_p([[[ZZ(1)]]], 2, ZZ)
    True

    """
    while u:
        if len(f) != 1:
            return False
        f = f[0]
        u -= 1

    return f == [K.one]


def dmp_one(u, K):
    """
    Return a multivariate one over ``K``.

    Examples
    ========

    >>> dmp_one(2, ZZ)
    [[[1]]]

    """
    return dmp_ground(K.one, u)


def dmp_ground(c, u):
    """
    Return a multivariate constant.

    Examples
    ========

    >>> dmp_ground(ZZ(3), 5)
    [[[[[[3]]]]]]
    >>> dmp_ground(ZZ(1), -1)
    1

    """
    if not c:
        return dmp_zero(u)

    for i in range(u + 1):
        c = [c]

    return c


def dmp_from_dict(f, u, K):
    """
    Create a ``K[X]`` polynomial from a :class:`dict`.

    Examples
    ========

    >>> dmp_from_dict({(0,): ZZ(7), (2,): ZZ(5), (4,): ZZ(1)}, 0, ZZ)
    [1, 0, 5, 0, 7]

    >>> dmp_from_dict({(0, 0): ZZ(3), (0, 1): ZZ(2), (2, 1): ZZ(1)}, 1, ZZ)
    [[1, 0], [], [2, 3]]

    """
    if not f:
        return dmp_zero(u)
    elif not u:
        h = []
        n, = max(f)

        for k in range(n, -1, -1):
            h.append(f.get((k,), K.zero))

        return dmp_strip(h, 0)

    coeffs = {}

    for monom, coeff in f.items():
        head, tail = monom[0], monom[1:]

        if head in coeffs:
            coeffs[head][tail] = coeff
        else:
            coeffs[head] = {tail: coeff}

    n, v, h = max(coeffs), u - 1, []

    for k in range(n, -1, -1):
        coeff = coeffs.get(k)

        if coeff is not None:
            h.append(dmp_from_dict(coeff, v, K))
        else:
            h.append(dmp_zero(v))

    return dmp_strip(h, u)


def dmp_to_dict(f, u):
    """
    Convert a ``K[X]`` polynomial to a :class:`dict`.

    Examples
    ========

    >>> dmp_to_dict([[ZZ(1), ZZ(0)], [], [ZZ(2), ZZ(3)]], 1)
    {(0, 0): 3, (0, 1): 2, (2, 1): 1}

    """
    n, v, result = dmp_degree_in(f, 0, u), u - 1, {}

    if n == -oo:
        n = -1

    for k in range(n + 1):
        if f[n - k]:
            if u:
                h = dmp_to_dict(f[n - k], v)
                for exp, coeff in h.items():
                    result[Monomial((k,) + exp)] = coeff
            else:
                result[Monomial((k,))] = f[n - k]

    return result


def dmp_permute(f, P, u, K):
    """
    Return a polynomial in ``K[x_{P(1)},..,x_{P(n)}]``.

    Examples
    ========

    >>> dmp_permute([[[ZZ(2)], [ZZ(1), ZZ(0)]], []], [1, 0, 2], 2, ZZ)
    [[[2], []], [[1, 0], []]]

    """
    F, H = dmp_to_dict(f, u), {}

    for exp, coeff in F.items():
        new_exp = [0]*len(exp)

        for e, p in zip(exp, P):
            new_exp[p] = e

        H[tuple(new_exp)] = coeff

    return dmp_from_dict(H, u, K)


def dmp_terms_gcd(f, u, K):
    """Remove GCD of terms from ``f`` in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    G, f = f.terms_gcd()
    return G, f.to_dense()


def dmp_slice_in(f, m, n, j, u, K):
    """Take a continuous subsequence of terms of ``f`` in ``x_j`` in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    return f.slice(m, n, x=j).to_dense()


def dup_random(n, a, b, K, percent=None):
    """
    Return a polynomial of degree ``n`` with coefficients in ``[a, b]``.

    If ``percent`` is a natural number less than 100 then only approximately
    the given percentage of elements will be non-zero.

    Examples
    ========

    >>> dup_random(3, -10, 10, ZZ)  # doctest: +SKIP
    [-2, -8, 9, -4]

    """
    if percent is None:
        percent = 100//(b - a)
    percent = min(max(0, percent), 100)
    nz = ((n + 1)*percent)//100

    f = []
    while len(f) < n + 1:
        v = K.convert(random.randint(a, b))
        if v:
            f.append(v)

    if nz:
        f[-nz:] = [K.zero]*nz
        lt = f.pop(0)
        random.shuffle(f)
        f.insert(0, lt)

    return f
