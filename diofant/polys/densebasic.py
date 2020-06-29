"""Basic tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from ..core import oo
from .monomials import Monomial


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
        return dmp_ground(0, u)

    v = u - 1

    for i, c in enumerate(f):
        if not dmp_zero_p(c, v):
            return f[i:]
    return dmp_ground(0, u)


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
        r = []

        for i in range(u):
            r = [r]

        return r

    for i in range(u + 1):
        c = [c]

    return c


def dmp_from_dict(f, u, K):
    """Create a ``K[X]`` polynomial from a :class:`dict`."""
    if not f:
        return dmp_ground(0, u)
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
            h.append(dmp_ground(0, v))

    return dmp_strip(h, u)


def dmp_to_dict(f, u):
    """Convert a ``K[X]`` polynomial to a :class:`dict`."""
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
