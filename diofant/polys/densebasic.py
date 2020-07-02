"""Basic tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from ..core import oo
from .monomials import Monomial


def dmp_degree_in(f, j, u):
    """Return the leading degree of ``f`` in ``x_j`` in ``K[X]``."""
    assert not j
    return -oo if dmp_zero_p(f, u) else len(f) - 1


def dmp_strip(f, u):
    """Remove leading zeros from ``f`` in ``K[X]``."""
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


def dmp_zero_p(f, u):
    """Return ``True`` if ``f`` is zero in ``K[X]``."""
    while u:
        if len(f) != 1:
            return False

        f = f[0]
        u -= 1

    return not f


def dmp_ground(c, u):
    """Return a multivariate constant."""
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
