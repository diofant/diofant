"""Tests for dense recursive polynomials' tools."""

import pytest

from diofant import EX, QQ, ZZ, DomainError, I, ring, sqrt


__all__ = ()


def test_dup_real_imag():
    R, x, y = ring('x y', ZZ)

    assert R.dup_real_imag(R.zero) == (0, 0)
    assert R.dup_real_imag(R.one) == (1, 0)

    assert R.dup_real_imag(x + 1) == (x + 1, y)
    assert R.dup_real_imag(x + 2) == (x + 2, y)

    assert R.dup_real_imag(x**2 + 2*x + 3) == (x**2 - y**2 + 2*x + 3,
                                               2*x*y + 2*y)

    f = x**3 + x**2 + x + 1

    assert R.dup_real_imag(f) == (x**3 + x**2 - 3*x*y**2 + x - y**2 + 1,
                                  3*x**2*y + 2*x*y - y**3 + y)

    R, x, y = ring('x y', EX)
    pytest.raises(DomainError, lambda: R.dup_real_imag(x + 1))

    R = QQ.algebraic_field(I).inject('x', 'y')
    x, y = R.to_ground().gens

    f = R.x**4 + I*R.x**3 - R.x + 1
    r = x**4 - 6*x**2*y**2 - 3*x**2*y - x + y**4 + y**3 + 1
    i = 4*x**3*y + x**3 - 4*x*y**3 - 3*x*y**2 - y

    assert R.dup_real_imag(f) == (r, i)

    K = QQ.algebraic_field(sqrt(2))
    R = K.inject('x', 'y')
    x, y = R.gens

    f = R.x**2 + sqrt(2)*R.x - 1
    assert R.dup_real_imag(f) == (x**2 - y**2 + sqrt(2)*x - 1, 2*x*y + sqrt(2)*y)

    K = K.algebraic_field(I)
    R = K.inject('x', 'y')
    x, y = R.to_ground().gens

    f = R.x**2 + 2*sqrt(2)*I*R.x - 1 + I
    assert R.dup_real_imag(f) == (x**2 - y**2 - 2*sqrt(2)*y - 1,
                                  2*x*y + 2*sqrt(2)*x + 1)


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


def test_dup_decompose():
    R, x = ring('x', ZZ)

    assert R(1).decompose() == [1]

    assert x.decompose() == [x]
    assert (x**3).decompose() == [x**3]

    assert (x**4).decompose() == [x**2, x**2]
    assert (x**6).decompose() == [x**3, x**2]

    assert (7*x**4 + 1).decompose() == [7*x**2 + 1, x**2]
    assert (4*x**4 + 3*x**2 + 2).decompose() == [4*x**2 + 3*x + 2, x**2]

    f = x**12 + 20*x**10 + 150*x**8 + 500*x**6 + 625*x**4 - 2*x**3 - 10*x + 9

    assert f.decompose() == [x**4 - 2*x + 9, x**3 + 5*x]

    f = 2*x**12 + 40*x**10 + 300*x**8 + 1000*x**6 + 1250*x**4 - 4*x**3 - 20*x + 18

    assert f.decompose() == [2*x**4 - 4*x + 18, x**3 + 5*x]

    f = (x**12 + 20*x**10 - 8*x**9 + 150*x**8 - 120*x**7 + 524*x**6 -
         600*x**5 + 865*x**4 - 1034*x**3 + 600*x**2 - 170*x + 29)

    assert f.decompose() == [x**4 - 8*x**3 + 24*x**2 - 34*x + 29, x**3 + 5*x]

    Rt, t = ring('t', ZZ)
    R, x = ring('x', Rt)

    f = ((6*t**2 - 42)*x**4 + (48*t**2 + 96)*x**3 +
         (144*t**2 + 648*t + 288)*x**2 + (624*t**2 + 864*t + 384)*x +
         108*t**3 + 312*t**2 + 432*t + 192)

    assert f.decompose() == [f]
