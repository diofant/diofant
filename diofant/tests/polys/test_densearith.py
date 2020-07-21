"""Tests for dense recursive polynomials' arithmetics."""

from diofant import QQ, ZZ, ring
from diofant.polys.specialpolys import f_polys


__all__ = ()


def test_dmp_add_term():
    R, x = ring('x', ZZ)

    f = R(0)

    assert f + 0 == 0
    assert f + 1 == 1
    assert f + x == x
    assert f + x**2 == x**2

    f = x**2 + x + 1

    assert f + 1 == x**2 + x + 2
    assert f + x == x**2 + 2*x + 1
    assert f + x**2 == 2*x**2 + x + 1

    assert f + x**3 == x**3 + x**2 + x + 1
    assert f + x**4 == x**4 + x**2 + x + 1
    assert f + x**5 == x**5 + x**2 + x + 1
    assert f + x**6 == x**6 + x**2 + x + 1

    assert f - x**2 == x + 1

    f = x**2 - 1

    assert f + 2*x**4 == 2*x**4 + x**2 - 1

    R, x, y, z = ring('x y z', ZZ)

    f = f_polys()[0]

    assert f + 0 == f

    f = x*y + 1

    assert f + 2*x**2 == 2*x**2 + x*y + 1

    R, x, y, z = ring('x y z', QQ)

    f = f.set_ring(R)/7

    assert f + 0 == f


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
