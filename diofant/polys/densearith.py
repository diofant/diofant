"""Arithmetics for dense recursive polynomials in ``K[x]`` or ``K[X]``."""


def dmp_mul_ground(f, c, u, K):
    """Multiply ``f`` by a constant value in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    r = f.mul_ground(c)
    return r.to_dense()


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


def dmp_mul(f, g, u, K):
    """Multiply dense polynomials in ``K[X]``."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return (f*g).to_dense()


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
