"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densearith import dmp_add, dmp_mul, dmp_mul_ground, dup_add, dup_mul
from .densebasic import dmp_zero_p


def dmp_diff_in(f, m, j, u, K):
    """``m``-th order derivative in ``x_j`` of a polynomial in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    return f.diff(x=j, m=m).to_dense()


def dmp_eval_in(f, a, j, u, K):
    """Evaluate a polynomial at ``x_j = a`` in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    r = f.eval(x=j, a=a)
    if ring.is_multivariate:
        r = r.to_dense()
    return r


def dmp_ground_monic(f, u, K):
    """Divide all coefficients by ``LC(f)`` in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    f = f.monic()
    return f.to_dense()


def dmp_ground_primitive(f, u, K):
    """Compute content and the primitive form of ``f`` in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    cont, p = f.primitive()
    return cont, p.to_dense()


def dup_transform(f, p, q, K):
    """
    Evaluate functional transformation ``q**n * f(p/q)`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> R.dup_transform(x**2 - 2*x + 1, x**2 + 1, x - 1)
    x**4 - 2*x**3 + 5*x**2 - 4*x + 4

    """
    if not f:
        return []

    n = len(f) - 1
    h, Q = [f[0]], [[K.one]]

    for i in range(n):
        Q.append(dup_mul(Q[-1], q, K))

    for c, q in zip(f[1:], Q[1:]):
        h = dup_mul(h, p, K)
        q = dmp_mul_ground(q, c, 0, K)
        h = dup_add(h, q, K)

    return h


def dmp_compose(f, g, u, K):
    """
    Evaluate functional composition ``f(g)`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dmp_compose(x*y + 2*x + y, y)
    y**2 + 3*y

    """
    if dmp_zero_p(f, u):
        return f

    h = [f[0]]

    for c in f[1:]:
        h = dmp_mul(h, g, u, K)
        h = dmp_add(h, [c], u, K)

    return h


def dmp_clear_denoms(f, u, K, convert=False):
    """Clear denominators."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    common, f = f.clear_denoms(convert=convert)
    return common, f.to_dense()
