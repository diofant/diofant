"""Efficient functions for generating orthogonal polynomials."""

from ..core import Dummy
from ..domains import QQ, ZZ
from .constructor import construct_domain
from .polytools import Poly, PurePoly


def _jacobi(n, a, b, K):
    """Low-level implementation of Jacobi polynomials."""
    ring = K.poly_ring('_0')
    x = ring._0

    j0 = ring.one
    if n < 1:
        return j0
    j1 = ((a + b + K(2))*x + (a - b))/K(2)

    for i in range(2, n + 1):
        den = K(i)*(a + b + i)*(a + b + K(2)*i - K(2))
        f0 = (a + b + K(2)*i - K.one) * (a*a - b*b) / (K(2)*den)
        f1 = (a + b + K(2)*i - K.one) * (a + b + K(2)*i - K(2)) * (a + b + K(2)*i) / (K(2)*den)
        f2 = (a + i - K.one)*(b + i - K.one)*(a + b + K(2)*i) / den
        j0, j1 = j1, j1*f0 + j1*x*f1 - j0*f2

    return j1


def jacobi_poly(n, a, b, x=None, **args):
    """Generates Jacobi polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError(f"can't generate Jacobi polynomial of degree {n}")

    K, v = construct_domain([a, b], field=True)
    poly = _jacobi(int(n), v[0], v[1], K)

    if x is not None:
        poly = Poly(poly, x, domain=K)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=K)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def _gegenbauer(n, a, K):
    """Low-level implementation of Gegenbauer polynomials."""
    ring = K.poly_ring('_0')
    x = ring._0

    g0 = ring.one
    if n < 1:
        return g0
    g1 = K(2)*a*x

    for i in range(2, n + 1):
        f1 = K(2) * (i + a - K.one) / i
        f2 = (i + K(2)*a - K(2)) / i
        g0, g1 = g1, g1*x*f1 - g0*f2

    return g1


def gegenbauer_poly(n, a, x=None, **args):
    """Generates Gegenbauer polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError(
            f"can't generate Gegenbauer polynomial of degree {n}")

    K, a = construct_domain(a, field=True)
    poly = _gegenbauer(int(n), a, K)

    if x is not None:
        poly = Poly(poly, x, domain=K)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=K)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def _chebyshevt(n, K):
    """Low-level implementation of Chebyshev polynomials of the 1st kind."""
    ring = K.poly_ring('_0')
    x = ring._0

    c0 = ring.one
    if n < 1:
        return c0
    c1 = x

    for i in range(2, n + 1):
        a = c1*x*K(2)
        c0, c1 = c1, a - c0

    return c1


def chebyshevt_poly(n, x=None, **args):
    """Generates Chebyshev polynomial of the first kind of degree `n` in `x`."""
    if n < 0:
        raise ValueError(
            f"can't generate 1st kind Chebyshev polynomial of degree {n}")

    poly = _chebyshevt(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, x, domain=ZZ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def _chebyshevu(n, K):
    """Low-level implementation of Chebyshev polynomials of the 2nd kind."""
    ring = K.poly_ring('_0')
    x = ring._0

    seq = [ring.one, K(2)*x]

    for i in range(2, n + 1):
        a = seq[-1]*x*K(2)
        seq.append(a - seq[-2])

    return seq[n]


def chebyshevu_poly(n, x=None, **args):
    """Generates Chebyshev polynomial of the second kind of degree `n` in `x`."""
    if n < 0:
        raise ValueError(
            f"can't generate 2nd kind Chebyshev polynomial of degree {n}")

    poly = _chebyshevu(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, x, domain=ZZ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def _hermite(n, K):
    """Low-level implementation of Hermite polynomials."""
    ring = K.poly_ring('_0')
    x = ring._0

    h0 = ring.one
    if n < 1:
        return h0
    h1 = K(2)*x

    for i in range(2, n + 1):
        a = h1*x
        b = h0*K(i - 1)

        h0, h1 = h1, (a - b)*K(2)

    return h1


def hermite_poly(n, x=None, **args):
    """Generates Hermite polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError(f"can't generate Hermite polynomial of degree {n}")

    poly = _hermite(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, x, domain=ZZ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def _legendre(n, K):
    """Low-level implementation of Legendre polynomials."""
    ring = K.poly_ring('_0')
    x = ring._0

    l0 = ring.one
    if n < 1:
        return l0
    l1 = x

    for i in range(2, n + 1):
        l0, l1 = l1, l1*x*K(2*i - 1, i) - l0*K(i - 1, i)

    return l1


def legendre_poly(n, x=None, **args):
    """Generates Legendre polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError(f"can't generate Legendre polynomial of degree {n}")

    poly = _legendre(int(n), QQ)

    if x is not None:
        poly = Poly(poly, x, domain=QQ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=QQ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def _laguerre(n, alpha, K):
    """Low-level implementation of Laguerre polynomials."""
    ring = K.poly_ring('_0')
    x = ring._0

    l0, l1 = ring.zero, ring.one

    for i in range(1, n + 1):
        l0, l1 = l1, l1*(-K.one/i*x + alpha/i + K(2*i - 1)/i) - l0*(alpha/i + K(i - 1)/i)

    return l1


def laguerre_poly(n, x=None, alpha=None, **args):
    """Generates Laguerre polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError(f"can't generate Laguerre polynomial of degree {n}")

    if alpha is not None:
        K, alpha = construct_domain(
            alpha, field=True)  # XXX: ground_field=True
    else:
        K, alpha = QQ, QQ(0)

    poly = _laguerre(int(n), alpha, K)

    if x is not None:
        poly = Poly(poly, x, domain=K)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=K)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def _spherical_bessel_fn(n, K):
    """Low-level implementation of fn(n, x)."""
    ring = K.poly_ring('_0')
    x = ring._0

    s0 = ring.one
    if n < 1:
        return s0*x
    s1 = x

    for i in range(2, n + 1):
        s0, s1 = s1, s1*x*K(2*i - 1) - s0

    return s1*x


def _spherical_bessel_fn_minus(n, K):
    """Low-level implementation of fn(-n, x)."""
    ring = K.poly_ring('_0')
    x = ring._0

    s0, s1 = x, ring.zero

    for i in range(2, n + 1):
        s0, s1 = s1, s1*x*K(3 - 2*i) - s0

    return s1


def spherical_bessel_fn(n, x=None, **args):
    """
    Coefficients for the spherical Bessel functions.

    Those are only needed in the jn() function.

    The coefficients are calculated from:

    fn(0, z) = 1/z
    fn(1, z) = 1/z**2
    fn(n-1, z) + fn(n+1, z) == (2*n+1)/z * fn(n, z)

    Examples
    ========

    >>> spherical_bessel_fn(1, z)
    z**(-2)
    >>> spherical_bessel_fn(2, z)
    -1/z + 3/z**3
    >>> spherical_bessel_fn(3, z)
    -6/z**2 + 15/z**4
    >>> spherical_bessel_fn(4, z)
    1/z - 45/z**3 + 105/z**5

    """
    if n < 0:
        poly = _spherical_bessel_fn_minus(-int(n), ZZ)
    else:
        poly = _spherical_bessel_fn(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, 1/x, domain=ZZ)
    else:
        poly = PurePoly(poly, 1/Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly
