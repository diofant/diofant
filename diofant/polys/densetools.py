"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densearith import (dmp_add, dmp_mul, dmp_mul_ground, dmp_neg, dmp_sub,
                         dup_add, dup_mul)
from .densebasic import dmp_ground, dmp_to_dict, dmp_zero, dmp_zero_p
from .polyerrors import DomainError


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


def dup_real_imag(f, K):
    """
    Return bivariate polynomials ``f1`` and ``f2``, such that ``f = f1 + f2*I``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dup_real_imag(x**3 + x**2 + x + 1)
    (x**3 + x**2 - 3*x*y**2 + x - y**2 + 1, 3*x**2*y + 2*x*y - y**3 + y)

    >>> R, x, y = ring('x y', QQ.algebraic_field(I))

    >>> R.dup_real_imag(x**2 + I*x - 1)
    (x**2 - y**2 - y - 1, 2*x*y + x)

    """
    if K.is_ComplexAlgebraicField:
        K0 = K.domain
        r1, i1 = dup_real_imag([_.real for _ in f], K0)
        r2, i2 = dup_real_imag([_.imag for _ in f], K0)
        return dmp_add(r1, dmp_neg(i2, 1, K0), 1, K0), dmp_add(r2, i1, 1, K0)
    elif not K.is_IntegerRing and not K.is_RationalField and not K.is_RealAlgebraicField:
        raise DomainError(f'computing real and imaginary parts is not supported over {K}')

    f1 = dmp_zero(1)
    f2 = dmp_zero(1)

    if not f:
        return f1, f2

    g = [[[K.one, K.zero]], [[K.one], []]]
    h = dmp_ground(f[0], 2)

    for c in f[1:]:
        h = dmp_mul(h, g, 2, K)
        h = dmp_add(h, [dmp_ground(c, 1)], 2, K)

    H = dmp_to_dict(h, 0)

    for (k,), h in H.items():
        m = k % 4

        if not m:
            f1 = dmp_add(f1, h, 1, K)
        elif m == 1:
            f2 = dmp_add(f2, h, 1, K)
        elif m == 2:
            f1 = dmp_sub(f1, h, 1, K)
        else:
            f2 = dmp_sub(f2, h, 1, K)

    return f1, f2


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
