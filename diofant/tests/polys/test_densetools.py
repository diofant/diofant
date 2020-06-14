"""Tests for dense recursive polynomials' tools."""

from diofant import ZZ, ring


__all__ = ()


def test_dmp_compose():
    R, x = ring('x', ZZ)

    assert R.dmp_compose(R(0), R(0)) == 0
    assert R.dmp_compose(R(0), R(1)) == 0
    assert R.dmp_compose(R(0), x + 2) == 0

    assert R.dmp_compose(R(1), R(0)) == 1

    assert R.dmp_compose(x**2 + 2*x, R(0)) == 0
    assert R.dmp_compose(x**2 + 2*x + 1, R(0)) == 1

    assert R.dmp_compose(x**2 + 2*x + 1, R(1)) == 4
    assert R.dmp_compose(x**2 + 2*x + 1, R(7)) == 64

    assert R.dmp_compose(x**2 + 2*x + 1, x - 1) == x**2
    assert R.dmp_compose(x**2 + 2*x + 1, x + 1) == x**2 + 4*x + 4
    assert R.dmp_compose(x**2 + 2*x + 1, x**2 + 2*x + 1) == x**4 + 4*x**3 + 8*x**2 + 8*x + 4

    assert R.dmp_compose(x**2 + x, x - 1) == x**2 - x

    R, x, y, z = ring('x y z', ZZ)

    assert R.dmp_compose(R(0), R(0)) == 0
    assert R.dmp_compose(R(0), R(1)) == 0
    assert R.dmp_compose(R(1), R(0)) == 1
    assert R.dmp_compose(R(0), x + 2) == 0

    R, x, y = ring('x y', ZZ)

    assert R.dmp_compose(x**2 + 2*x, R(0)) == 0
    assert R.dmp_compose(x**2 + 2*x + 1, R(0)) == 1

    assert R.dmp_compose(x**2 + 2*x + 1, R(1)) == 4
    assert R.dmp_compose(x**2 + 2*x + 1, R(7)) == 64

    assert R.dmp_compose(x**2 + 2*x + 1, x - 1) == x**2
    assert R.dmp_compose(x**2 + 2*x + 1, x + 1) == x**2 + 4*x + 4

    assert R.dmp_compose(x**2 + 2*x + 1, x**2 + 2*x + 1) == x**4 + 4*x**3 + 8*x**2 + 8*x + 4
