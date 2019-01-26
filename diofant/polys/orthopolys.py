"""Efficient functions for generating orthogonal polynomials. """

from ..core import Dummy
from ..domains import QQ, ZZ
from .constructor import construct_domain
from .densearith import dmp_mul_ground, dup_add, dup_lshift, dup_mul, dup_sub
from .polytools import Poly, PurePoly


__all__ = ('jacobi_poly', 'chebyshevt_poly', 'chebyshevu_poly', 'hermite_poly',
           'legendre_poly', 'laguerre_poly', 'spherical_bessel_fn',
           'gegenbauer_poly')


def dup_jacobi(n, a, b, K):
    """Low-level implementation of Jacobi polynomials."""
    seq = [[K.one], [(a + b + K(2))/K(2), (a - b)/K(2)]]

    for i in range(2, n + 1):
        den = K(i)*(a + b + i)*(a + b + K(2)*i - K(2))
        f0 = (a + b + K(2)*i - K.one) * (a*a - b*b) / (K(2)*den)
        f1 = (a + b + K(2)*i - K.one) * (a + b + K(2)*i - K(2)) * (a + b + K(2)*i) / (K(2)*den)
        f2 = (a + i - K.one)*(b + i - K.one)*(a + b + K(2)*i) / den
        p0 = dmp_mul_ground(seq[-1], f0, 0, K)
        p1 = dmp_mul_ground(dup_lshift(seq[-1], 1, K), f1, 0, K)
        p2 = dmp_mul_ground(seq[-2], f2, 0, K)
        seq.append(dup_sub(dup_add(p0, p1, K), p2, K))

    return seq[n]


def jacobi_poly(n, a, b, x=None, **args):
    """Generates Jacobi polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError("can't generate Jacobi polynomial of degree %s" % n)

    K, v = construct_domain([a, b], field=True)
    poly = dup_jacobi(int(n), v[0], v[1], K)

    if x is not None:
        poly = Poly(poly, x, domain=K)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=K)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def dup_gegenbauer(n, a, K):
    """Low-level implementation of Gegenbauer polynomials."""
    seq = [[K.one], [K(2)*a, K.zero]]

    for i in range(2, n + 1):
        f1 = K(2) * (i + a - K.one) / i
        f2 = (i + K(2)*a - K(2)) / i
        p1 = dmp_mul_ground(dup_lshift(seq[-1], 1, K), f1, 0, K)
        p2 = dmp_mul_ground(seq[-2], f2, 0, K)
        seq.append(dup_sub(p1, p2, K))

    return seq[n]


def gegenbauer_poly(n, a, x=None, **args):
    """Generates Gegenbauer polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError(
            "can't generate Gegenbauer polynomial of degree %s" % n)

    K, a = construct_domain(a, field=True)
    poly = dup_gegenbauer(int(n), a, K)

    if x is not None:
        poly = Poly(poly, x, domain=K)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=K)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def dup_chebyshevt(n, K):
    """Low-level implementation of Chebyshev polynomials of the 1st kind."""
    seq = [[K.one], [K.one, K.zero]]

    for i in range(2, n + 1):
        a = dmp_mul_ground(dup_lshift(seq[-1], 1, K), K(2), 0, K)
        seq.append(dup_sub(a, seq[-2], K))

    return seq[n]


def chebyshevt_poly(n, x=None, **args):
    """Generates Chebyshev polynomial of the first kind of degree `n` in `x`."""
    if n < 0:
        raise ValueError(
            "can't generate 1st kind Chebyshev polynomial of degree %s" % n)

    poly = dup_chebyshevt(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, x, domain=ZZ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def dup_chebyshevu(n, K):
    """Low-level implementation of Chebyshev polynomials of the 2nd kind."""
    seq = [[K.one], [K(2), K.zero]]

    for i in range(2, n + 1):
        a = dmp_mul_ground(dup_lshift(seq[-1], 1, K), K(2), 0, K)
        seq.append(dup_sub(a, seq[-2], K))

    return seq[n]


def chebyshevu_poly(n, x=None, **args):
    """Generates Chebyshev polynomial of the second kind of degree `n` in `x`."""
    if n < 0:
        raise ValueError(
            "can't generate 2nd kind Chebyshev polynomial of degree %s" % n)

    poly = dup_chebyshevu(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, x, domain=ZZ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def dup_hermite(n, K):
    """Low-level implementation of Hermite polynomials."""
    seq = [[K.one], [K(2), K.zero]]

    for i in range(2, n + 1):
        a = dup_lshift(seq[-1], 1, K)
        b = dmp_mul_ground(seq[-2], K(i - 1), 0, K)

        c = dmp_mul_ground(dup_sub(a, b, K), K(2), 0, K)

        seq.append(c)

    return seq[n]


def hermite_poly(n, x=None, **args):
    """Generates Hermite polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError("can't generate Hermite polynomial of degree %s" % n)

    poly = dup_hermite(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, x, domain=ZZ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def dup_legendre(n, K):
    """Low-level implementation of Legendre polynomials."""
    seq = [[K.one], [K.one, K.zero]]

    for i in range(2, n + 1):
        a = dmp_mul_ground(dup_lshift(seq[-1], 1, K), K(2*i - 1, i), 0, K)
        b = dmp_mul_ground(seq[-2], K(i - 1, i), 0, K)

        seq.append(dup_sub(a, b, K))

    return seq[n]


def legendre_poly(n, x=None, **args):
    """Generates Legendre polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError("can't generate Legendre polynomial of degree %s" % n)

    poly = dup_legendre(int(n), QQ)

    if x is not None:
        poly = Poly(poly, x, domain=QQ)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=QQ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def dup_laguerre(n, alpha, K):
    """Low-level implementation of Laguerre polynomials."""
    seq = [[K.zero], [K.one]]

    for i in range(1, n + 1):
        a = dup_mul(seq[-1], [-K.one/i, alpha/i + K(2*i - 1)/i], K)
        b = dmp_mul_ground(seq[-2], alpha/i + K(i - 1)/i, 0, K)

        seq.append(dup_sub(a, b, K))

    return seq[-1]


def laguerre_poly(n, x=None, alpha=None, **args):
    """Generates Laguerre polynomial of degree `n` in `x`."""
    if n < 0:
        raise ValueError("can't generate Laguerre polynomial of degree %s" % n)

    if alpha is not None:
        K, alpha = construct_domain(
            alpha, field=True)  # XXX: ground_field=True
    else:
        K, alpha = QQ, QQ(0)

    poly = dup_laguerre(int(n), alpha, K)

    if x is not None:
        poly = Poly(poly, x, domain=K)
    else:
        poly = PurePoly(poly, Dummy('x'), domain=K)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def dup_spherical_bessel_fn(n, K):
    """Low-level implementation of fn(n, x)."""
    seq = [[K.one], [K.one, K.zero]]

    for i in range(2, n + 1):
        a = dmp_mul_ground(dup_lshift(seq[-1], 1, K), K(2*i - 1), 0, K)
        seq.append(dup_sub(a, seq[-2], K))

    return dup_lshift(seq[n], 1, K)


def dup_spherical_bessel_fn_minus(n, K):
    """Low-level implementation of fn(-n, x)."""
    seq = [[K.one, K.zero], [K.zero]]

    for i in range(2, n + 1):
        a = dmp_mul_ground(dup_lshift(seq[-1], 1, K), K(3 - 2*i), 0, K)
        seq.append(dup_sub(a, seq[-2], K))

    return seq[n]


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
        poly = dup_spherical_bessel_fn_minus(-int(n), ZZ)
    else:
        poly = dup_spherical_bessel_fn(int(n), ZZ)

    if x is not None:
        poly = Poly(poly, 1/x, domain=ZZ)
    else:
        poly = PurePoly(poly, 1/Dummy('x'), domain=ZZ)

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly
