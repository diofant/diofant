"""Tests for dense recursive polynomials' tools."""

from diofant import ZZ, ring


__all__ = ()


def test_dup_transform():
    R, x = ring('x', ZZ)

    assert R.dup_transform(R(0), R(0), x + 1) == 0
    assert R.dup_transform(R(0), R(1), x + 1) == 0
    assert R.dup_transform(R(0), x + 2, x + 1) == 0

    assert (R.dup_transform(6*x**4 - 5*x**3 + 4*x**2 - 3*x + 17,
                            x**2 - 3*x + 4, 2*x - 3) ==
            6*x**8 - 82*x**7 + 541*x**6 - 2205*x**5 + 6277*x**4 -
            12723*x**3 + 17191*x**2 - 13603*x + 4773)


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
