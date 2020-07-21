"""Tests for dense recursive polynomials' arithmetics."""

from diofant import ZZ, ring


__all__ = ()


def test_dmp_add_mul():
    R, x = ring('x', ZZ)

    assert (x**2 + 2*x + 3 +
            (3*x**2 + 2*x + 1)*(x + 2)) == 3*x**3 + 9*x**2 + 7*x + 5
    assert x**2 - 1 + (x - 2)*(x + 2) == 2*x**2 - 5

    R, x, y = ring('x y', ZZ)

    assert (x*y + 2*x + 3 +
            (3*x + 2*y + 1)*(x + 2)) == 3*x**2 + 3*x*y + 9*x + 4*y + 5
    assert x**2 + y + x*(x + 2) == 2*x**2 + y + 2*x


def test_dmp_sub_mul():
    R, x = ring('x', ZZ)

    assert (x**2 + 2*x + 3 - (3*x**2 + 2*x + 1)*(x + 2) ==
            -3*x**3 - 7*x**2 - 3*x + 1)
    assert x**2 - 1 - (x - 2)*(x + 2) == 3

    R, x, y = ring('x y', ZZ)

    assert (x*y + 2*x + 3 - (3*x + 2*y + 1)*(x + 2) ==
            -3*x**2 - x*y - 5*x - 4*y + 1)
    assert x**2 + y - x*(x + 2) == -2*x + y
