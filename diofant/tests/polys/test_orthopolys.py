"""Tests for efficient functions for generating orthogonal polynomials."""

import pytest

from diofant import (ZZ, Rational, chebyshevt_poly, chebyshevu_poly,
                     hermite_poly, jacobi_poly, laguerre_poly, legendre_poly,
                     spherical_bessel_fn)
from diofant.abc import a, b, x
from diofant.polys.orthopolys import gegenbauer_poly


__all__ = ()


def test_jacobi_poly():
    pytest.raises(ValueError, lambda: jacobi_poly(-1, a, b, x))

    dom = ZZ.inject(a, b).field

    assert (jacobi_poly(1, a, b, x, polys=True) ==
            ((a/2 + b/2 + 1)*x + a/2 - b/2).as_poly(x, domain=dom))

    assert jacobi_poly(0, a, b, x) == 1
    assert jacobi_poly(1, a, b, x) == a/2 - b/2 + x*(a/2 + b/2 + 1)
    assert jacobi_poly(2, a, b, x) == (a**2/8 - a*b/4 - a/8 + b**2/8 - b/8 +
                                       x**2*(a**2/8 + a*b/4 + 7*a/8 +
                                             b**2/8 + 7*b/8 + Rational(3, 2)) +
                                       x*(a**2/4 + 3*a/4 - b**2/4 - 3*b/4) -
                                       Rational(1, 2))

    assert (jacobi_poly(1, a, b, polys=True) ==
            ((a/2 + b/2 + 1)*x + a/2 - b/2).as_poly(x, domain=dom))


def test_gegenbauer_poly():
    pytest.raises(ValueError, lambda: gegenbauer_poly(-1, a, x))

    dom = ZZ.inject(a).field

    assert gegenbauer_poly(
        1, a, x, polys=True) == (2*a*x).as_poly(x, domain=dom)

    assert gegenbauer_poly(0, a, x) == 1
    assert gegenbauer_poly(1, a, x) == 2*a*x
    assert gegenbauer_poly(2, a, x) == -a + x**2*(2*a**2 + 2*a)
    assert gegenbauer_poly(
        3, a, x) == x**3*(4*a**3/3 + 4*a**2 + 8*a/3) + x*(-2*a**2 - 2*a)

    assert gegenbauer_poly(1, Rational(1, 2), x) == x
    assert gegenbauer_poly(1, a, polys=True) == (2*a*x).as_poly(x, domain=dom)


def test_chebyshevt_poly():
    pytest.raises(ValueError, lambda: chebyshevt_poly(-1, x))

    assert chebyshevt_poly(1, x, polys=True) == x.as_poly()

    assert chebyshevt_poly(0, x) == 1
    assert chebyshevt_poly(1, x) == x
    assert chebyshevt_poly(2, x) == 2*x**2 - 1
    assert chebyshevt_poly(3, x) == 4*x**3 - 3*x
    assert chebyshevt_poly(4, x) == 8*x**4 - 8*x**2 + 1
    assert chebyshevt_poly(5, x) == 16*x**5 - 20*x**3 + 5*x
    assert chebyshevt_poly(6, x) == 32*x**6 - 48*x**4 + 18*x**2 - 1

    assert chebyshevt_poly(1, polys=True) == x.as_poly()


def test_chebyshevu_poly():
    pytest.raises(ValueError, lambda: chebyshevu_poly(-1, x))

    assert chebyshevu_poly(1, x, polys=True) == (2*x).as_poly()

    assert chebyshevu_poly(0, x) == 1
    assert chebyshevu_poly(1, x) == 2*x
    assert chebyshevu_poly(2, x) == 4*x**2 - 1
    assert chebyshevu_poly(3, x) == 8*x**3 - 4*x
    assert chebyshevu_poly(4, x) == 16*x**4 - 12*x**2 + 1
    assert chebyshevu_poly(5, x) == 32*x**5 - 32*x**3 + 6*x
    assert chebyshevu_poly(6, x) == 64*x**6 - 80*x**4 + 24*x**2 - 1

    assert chebyshevu_poly(1, polys=True) == (2*x).as_poly()


def test_hermite_poly():
    pytest.raises(ValueError, lambda: hermite_poly(-1, x))

    assert hermite_poly(1, x, polys=True) == (2*x).as_poly()

    assert hermite_poly(0, x) == 1
    assert hermite_poly(1, x) == 2*x
    assert hermite_poly(2, x) == 4*x**2 - 2
    assert hermite_poly(3, x) == 8*x**3 - 12*x
    assert hermite_poly(4, x) == 16*x**4 - 48*x**2 + 12
    assert hermite_poly(5, x) == 32*x**5 - 160*x**3 + 120*x
    assert hermite_poly(6, x) == 64*x**6 - 480*x**4 + 720*x**2 - 120

    assert hermite_poly(1, polys=True) == (2*x).as_poly()


def test_legendre_poly():
    pytest.raises(ValueError, lambda: legendre_poly(-1, x))

    assert legendre_poly(1, x, polys=True) == x.as_poly()

    assert legendre_poly(0, x) == 1
    assert legendre_poly(1, x) == x
    assert legendre_poly(2, x) == 3*x**2/2 - Rational(1, 2)
    assert legendre_poly(3, x) == 5*x**3/2 - 3*x/2
    assert legendre_poly(4, x) == 35*x**4/8 - 30*x**2/8 + Rational(3, 8)
    assert legendre_poly(5, x) == 63*x**5/8 - 70*x**3/8 + 15*x/8
    assert legendre_poly(6, x) == (231*x**6/16 - 315*x**4/16 +
                                   105*x**2/16 - Rational(5, 16))

    assert legendre_poly(1, polys=True) == x.as_poly()


def test_laguerre_poly():
    pytest.raises(ValueError, lambda: laguerre_poly(-1, x))

    assert laguerre_poly(1, x, polys=True) == (-x + 1).as_poly()

    assert laguerre_poly(0, x) == 1
    assert laguerre_poly(1, x) == -x + 1
    assert laguerre_poly(2, x) == x**2/2 - 2*x + 1
    assert laguerre_poly(3, x) == -x**3/6 + 3*x**2/2 - 3*x + 1
    assert laguerre_poly(4, x) == x**4/24 - 2*x**3/3 + 3*x**2 - 4*x + 1
    assert laguerre_poly(5, x) == (-x**5/120 + 5*x**4/24 - 5*x**3/3 +
                                   5*x**2 - 5*x + 1)
    assert laguerre_poly(6, x) == (x**6/720 - x**5/20 + 5*x**4/8 -
                                   10*x**3/3 + 15*x**2/2 - 6*x + 1)

    assert laguerre_poly(0, x, a) == 1
    assert laguerre_poly(1, x, a) == -x + a + 1
    assert laguerre_poly(2, x, a) == x**2/2 + (-a - 2)*x + a**2/2 + 3*a/2 + 1
    assert laguerre_poly(3, x, a) == (-x**3/6 + (a/2 + Rational(3, 2))*x**2 +
                                      (-a**2/2 - 5*a/2 - 3)*x + a**3/6 +
                                      a**2 + 11*a/6 + 1)

    assert laguerre_poly(1, x) == 1 - x
    assert laguerre_poly(1, polys=True) == (-x + 1).as_poly()


def test_spherical_bessel_fn():
    assert spherical_bessel_fn(1, x) == x**(-2)
    assert spherical_bessel_fn(1, polys=True)(1/x) == x**(-2)
