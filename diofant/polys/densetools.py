"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""


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
