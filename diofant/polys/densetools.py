"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""


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


def dmp_compose(f, g, u, K):
    """Evaluate functional composition ``f(g)`` in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    r = f.compose(0, g)
    return r.to_dense()


def dmp_clear_denoms(f, u, K, convert=False):
    """Clear denominators."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    common, f = f.clear_denoms(convert=convert)
    return common, f.to_dense()
