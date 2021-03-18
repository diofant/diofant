"""Tests for real and complex root isolation and refinement algorithms."""

import math

import pytest

from diofant import (EX, QQ, ZZ, DomainError, I, RefinementFailed, ring, sqrt,
                     subsets)
from diofant.polys.rootisolation import RealInterval


__all__ = ()


def test__real_imag():
    R, x, y = ring('x y', ZZ)
    Rx = R.drop(y)
    _y = R.symbols[1]

    assert Rx._real_imag(Rx.zero, _y) == (0, 0)
    assert Rx._real_imag(Rx.one, _y) == (1, 0)

    assert Rx._real_imag(Rx.x + 1, _y) == (x + 1, y)
    assert Rx._real_imag(Rx.x + 2, _y) == (x + 2, y)

    assert Rx._real_imag(Rx.x**2 + 2*Rx.x + 3, _y) == (x**2 - y**2 + 2*x + 3,
                                                       2*x*y + 2*y)

    f = Rx.x**3 + Rx.x**2 + Rx.x + 1

    assert Rx._real_imag(f, _y) == (x**3 + x**2 - 3*x*y**2 + x - y**2 + 1,
                                    3*x**2*y + 2*x*y - y**3 + y)

    R, x, y = ring('x y', EX)
    Rx = R.drop(y)

    pytest.raises(DomainError, lambda: Rx._real_imag(Rx.x + 1))

    R = QQ.algebraic_field(I).inject('x', 'y')
    Rx = R.drop(1)
    _y = R.symbols[1]
    x, y = R.to_ground().gens

    f = Rx.x**2 + I*Rx.x - 1
    r = x**2 - y**2 - y - 1
    i = 2*x*y + x

    assert Rx._real_imag(f, _y) == (r, i)

    f = Rx.x**4 + I*Rx.x**3 - Rx.x + 1
    r = x**4 - 6*x**2*y**2 - 3*x**2*y - x + y**4 + y**3 + 1
    i = 4*x**3*y + x**3 - 4*x*y**3 - 3*x*y**2 - y

    assert Rx._real_imag(f, _y) == (r, i)

    K = QQ.algebraic_field(sqrt(2))
    R = K.inject('x', 'y')
    Rx = R.drop(1)
    _y = R.symbols[1]
    x, y = R.gens

    f = Rx.x**2 + sqrt(2)*Rx.x - 1
    assert Rx._real_imag(f, _y) == (x**2 - y**2 + sqrt(2)*x - 1, 2*x*y + sqrt(2)*y)

    K = K.algebraic_field(I)
    R = K.inject('x', 'y')
    Rx = R.drop(1)
    _y = R.symbols[1]
    x, y = R.to_ground().gens

    f = Rx.x**2 + 2*sqrt(2)*I*Rx.x - 1 + I
    assert Rx._real_imag(f, _y) == (x**2 - y**2 - 2*sqrt(2)*y - 1,
                                    2*x*y + 2*sqrt(2)*x + 1)


def test__transform():
    R, x = ring('x', ZZ)

    assert R._transform(R(0), R(0), x + 1) == 0
    assert R._transform(R(0), R(1), x + 1) == 0
    assert R._transform(R(0), x + 2, x + 1) == 0

    assert R._transform(x**2 - 2*x + 1, x**2 + 1,
                        x - 1) == x**4 - 2*x**3 + 5*x**2 - 4*x + 4

    assert (R._transform(6*x**4 - 5*x**3 + 4*x**2 - 3*x + 17,
                         x**2 - 3*x + 4, 2*x - 3) ==
            6*x**8 - 82*x**7 + 541*x**6 - 2205*x**5 + 6277*x**4 -
            12723*x**3 + 17191*x**2 - 13603*x + 4773)


def test__reverse():
    R, x = ring('x', ZZ)

    assert R._reverse(x**3 + 2*x**2 + 3) == 3*x**3 + 2*x + 1
    assert R._reverse(x**3 + 2*x**2 + 3*x) == 3*x**2 + 2*x + 1


def test_sturm():
    R, x = ring('x', QQ)

    assert R(5).sturm() == [1]
    assert x.sturm() == [x, 1]

    f = x**3 - 2*x**2 + 3*x - 5

    assert f.sturm() == [f, 3*x**2 - 4*x + 3,
                         -10*x/9 + QQ(13, 3), -QQ(3303, 100)]

    f = x**3 - 2*x**2 + x - 3

    assert f.sturm() == [f, 3*x**2 - 4*x + 1, 2*x/9 + QQ(25, 9), QQ(-2079, 4)]

    pytest.raises(DomainError, lambda: x.set_domain(ZZ).sturm())

    F = ZZ.frac_field('pi')
    pi = F.pi
    R, x = ring('x', F)

    f = (1024/(15625*pi**8)*x**5 - 4096/(625*pi**8)*x**4 + 32/(15625*pi**4)*x**3
         - 128/(625*pi**4)*x**2 + x/62500 - F((1, 625)))

    assert f.sturm() == [x**3 - 100*x**2 + pi**4/64*x - 25*pi**4/16,
                         3*x**2 - 200*x + pi**4/64,
                         (F((20000, 9)) - pi**4/96)*x + 25*pi**4/18,
                         (-3686400000000*pi**4 - 11520000*pi**8 -
                          9*pi**12)/(26214400000000 - 245760000*pi**4 +
                                     576*pi**8)]


def test__sign_variations():
    R, x = ring('x', ZZ)

    assert R._sign_variations(R(0)) == 0
    assert R._sign_variations(x) == 0
    assert R._sign_variations(x**2 + 2) == 0
    assert R._sign_variations(x*(x**2 + 3)) == 0
    assert R._sign_variations(x**4 + 4*x**2 + 5) == 0

    assert R._sign_variations(2 - x**2) == 1
    assert R._sign_variations(x*(3 - x**2)) == 1
    assert R._sign_variations((5 - x**2)*(x**2 + 1)) == 1

    assert R._sign_variations(-x**2 - 4*x - 5) == 0
    assert R._sign_variations((x - 5)*(x + 1)) == 1
    assert R._sign_variations((x - 1)*(x + 5)) == 1
    assert R._sign_variations(+x**2 - 4*x + 5) == 2
    assert R._sign_variations(-x**2 + 4*x - 5) == 2
    assert R._sign_variations((5 - x)*(x + 1)) == 1
    assert R._sign_variations((1 - x)*(x + 5)) == 1
    assert R._sign_variations(+x**2 + 4*x + 5) == 0

    assert R._sign_variations(-x**4 - 4*x**2 - 5) == 0
    assert R._sign_variations((x**2 - 5)*(x**2 + 1)) == 1
    assert R._sign_variations((x - 1)*(x + 1)*(x**2 + 5)) == 1
    assert R._sign_variations(+x**4 - 4*x**2 + 5) == 2
    assert R._sign_variations(-x**4 + 4*x**2 - 5) == 2
    assert R._sign_variations((5 - x**2)*(x**2 + 1)) == 1
    assert R._sign_variations((1 - x)*(x + 1)*(x**2 + 5)) == 1
    assert R._sign_variations(+x**4 + 4*x**2 + 5) == 0


def test__root_upper_bound():
    R, x = ring('x', ZZ)

    assert R._root_upper_bound(+x - 1) == 4
    assert R._root_upper_bound(-x - 1) is None

    R, x = ring('x', QQ)

    assert R._root_upper_bound(+x - 1) == 4
    assert R._root_upper_bound(-x - 1) is None

    assert R._root_upper_bound(+x/2 - 1) is None
    assert R._root_upper_bound(-x/2 - 1) is None


def test__step_refine_real_root():
    R, x = ring('x', ZZ)

    assert R._step_refine_real_root(x + 1,
                                    (-2, 0, 1, 1)) == (x + 2, (0, -2, 1, 2))


def test__inner_refine_real_root():
    R, x = ring('x', ZZ)

    f = 2 - x**2
    r = (1, QQ(3, 2))

    assert R._inner_refine_real_root(f, (1, 2, 1, 1), steps=1) == r


def test__refine_real_root():
    R, x = ring('x', ZZ)

    f = x**2 - 2

    assert R._refine_real_root(f, 1, 1, steps=1) == (1, 1)
    assert R._refine_real_root(f, 1, 1, steps=9) == (1, 1)

    pytest.raises(ValueError, lambda: R._refine_real_root(f, -2, 2))

    s, t = 1, 2

    assert R._refine_real_root(f, s, t, steps=0) == (1, 2)
    assert R._refine_real_root(f, s, t, steps=1) == (1, QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=2) == (QQ(4, 3), QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=3) == (QQ(7, 5), QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=4) == (QQ(7, 5), QQ(10, 7))
    assert R._refine_real_root(f, s, t, eps=QQ(1, 100)) == (QQ(24, 17), QQ(17, 12))

    s, t = 1, QQ(3, 2)

    assert R._refine_real_root(f, s, t, steps=0) == (1, QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=1) == (QQ(4, 3), QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=2) == (QQ(7, 5), QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=3) == (QQ(7, 5), QQ(10, 7))
    assert R._refine_real_root(f, s, t, steps=4) == (QQ(7, 5), QQ(17, 12))

    s, t = 1, QQ(5, 3)

    assert R._refine_real_root(f, s, t, steps=0) == (1, QQ(5, 3))
    assert R._refine_real_root(f, s, t, steps=1) == (1, QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=2) == (QQ(7, 5), QQ(3, 2))
    assert R._refine_real_root(f, s, t, steps=3) == (QQ(7, 5), QQ(13, 9))
    assert R._refine_real_root(f, s, t, steps=4) == (QQ(7, 5), QQ(27, 19))

    s, t = -1, -2

    assert R._refine_real_root(f, s, t, steps=0) == (-2, -1)
    assert R._refine_real_root(f, s, t, steps=1) == (-QQ(3, 2), -1)
    assert R._refine_real_root(f, s, t, steps=2) == (-QQ(3, 2), -QQ(4, 3))
    assert R._refine_real_root(f, s, t, steps=3) == (-QQ(3, 2), -QQ(7, 5))
    assert R._refine_real_root(f, s, t, steps=4) == (-QQ(10, 7), -QQ(7, 5))

    pytest.raises(RefinementFailed, lambda: R._refine_real_root(f, 0, 1))

    s, t, u, v, w = 1, 2, QQ(24, 17), QQ(17, 12), QQ(7, 5)

    assert R._refine_real_root(f, s, t, eps=QQ(1, 100)) == (u, v)
    assert R._refine_real_root(f, s, t, steps=6) == (u, v)

    assert R._refine_real_root(f, s, t, eps=QQ(1, 100), steps=5) == (w, v)
    assert R._refine_real_root(f, s, t, eps=QQ(1, 100), steps=6) == (u, v)
    assert R._refine_real_root(f, s, t, eps=QQ(1, 100), steps=7) == (u, v)

    s, t, u, v = -2, -1, QQ(-3, 2), QQ(-4, 3)

    assert R._refine_real_root(f, s, t, disjoint=-5) == (s, t)
    assert R._refine_real_root(f, s, t, disjoint=-v) == (s, t)
    assert R._refine_real_root(f, s, t, disjoint=v) == (u, v)

    s, t, u, v = 1, 2, QQ(4, 3), QQ(3, 2)

    assert R._refine_real_root(f, s, t, disjoint=5) == (s, t)
    assert R._refine_real_root(f, s, t, disjoint=-u) == (s, t)
    assert R._refine_real_root(f, s, t, disjoint=u) == (u, v)

    f = x**2 - 3

    assert R._refine_real_root(f, 1, 2,
                               eps=QQ(1, 100)) == (QQ(19, 11), QQ(26, 15))

    R, x = ring('x', QQ)

    f = (x - QQ(1, 2))*(x + QQ(1, 2))

    assert R._refine_real_root(f, 0, 1, steps=1) == (QQ(1, 2), QQ(1, 2))

    D, y = ring('y', ZZ)
    R, x = ring('x', D)

    f = x**2 + y*x - 1

    pytest.raises(DomainError, lambda: R._refine_real_root(f, 0, 1))


def test__isolate_real_roots_sqf():
    R, x = ring('x', ZZ)

    assert R._isolate_real_roots_sqf(R(0)) == []
    assert R._isolate_real_roots_sqf(R(5)) == []

    assert R._isolate_real_roots_sqf(x) == [(0, 0)]

    f = x*(x + 1)

    assert R._isolate_real_roots_sqf(f) == [(-1, -1), (0, 0)]
    assert R._isolate_real_roots_sqf(f, inf=+1) == []
    assert R._isolate_real_roots_sqf(f, sup=-1) == [(-1, -1)]
    assert R._isolate_real_roots_sqf(f, sup=-2) == []

    f = x*(x - 1)

    assert R._isolate_real_roots_sqf(f) == [(0, 0), (1, 1)]

    assert R._isolate_real_roots_sqf(x**4 + x + 1) == []

    i = [(-2, -1), (1, 2)]
    f = x**2 - 2

    assert R._isolate_real_roots_sqf(+f) == i
    assert R._isolate_real_roots_sqf(-f) == i

    for r in range(2, 7):
        for s in (1, 10, -1, -10):
            f = R(math.prod(x - s*_ for _ in range(1, r)))
            ans = sorted((s*_, s*_) for _ in range(1, r))
            assert R._isolate_real_roots_sqf(f) == ans

    assert R._isolate_real_roots_sqf(x**2 - 5) == [(-3, -2), (2, 3)]
    assert R._isolate_real_roots_sqf(x**3 - 5) == [(1, 2)]
    assert R._isolate_real_roots_sqf(x**4 - 5) == [(-2, -1), (1, 2)]
    assert R._isolate_real_roots_sqf(x**5 - 5) == [(1, 2)]
    assert R._isolate_real_roots_sqf(x**6 - 5) == [(-2, -1), (1, 2)]
    assert R._isolate_real_roots_sqf(x**7 - 5) == [(1, 2)]
    assert R._isolate_real_roots_sqf(x**8 - 5) == [(-2, -1), (1, 2)]
    assert R._isolate_real_roots_sqf(x**9 - 5) == [(1, 2)]

    for roots in subsets(range(1, 4)):
        f = R(math.prod(x - r for r in roots))
        ans = sorted((_, _) for _ in roots)
        assert R._isolate_real_roots_sqf(f) == ans

    assert R._isolate_real_roots_sqf((x - 3)*(x - 2)*(x - 1)*(x + 1)*(x + 2)*(x + 3)*(2*x + 1)) == \
        [(-3, -3), (-2, -2), (-1, -1), (-1, 0), (1, 1), (2, 2), (3, 3)]
    assert R._isolate_real_roots_sqf((x - 3)*(x - 2)*(x - 1)*(x + 1)*(x + 2)*(x + 3)*(2*x - 1)*(2*x + 1)) == \
        [(-3, -3), (-2, -2), (-1, -1), (-1, 0), (0, 1), (1, 1), (2, 2), (3, 3)]

    f = 9*x**2 - 2

    assert R._isolate_real_roots_sqf(f) == \
        [(-1, 0), (0, 1)]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 10)) == \
        [(QQ(-1, 2), QQ(-3, 7)), (QQ(3, 7), QQ(1, 2))]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 100)) == \
        [(QQ(-9, 19), QQ(-8, 17)), (QQ(8, 17), QQ(9, 19))]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 1000)) == \
        [(QQ(-33, 70), QQ(-8, 17)), (QQ(8, 17), QQ(33, 70))]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 10000)) == \
        [(QQ(-33, 70), QQ(-107, 227)), (QQ(107, 227), QQ(33, 70))]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 100000)) == \
        [(QQ(-305, 647), QQ(-272, 577)), (QQ(272, 577), QQ(305, 647))]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 1000000)) == \
        [(QQ(-1121, 2378), QQ(-272, 577)), (QQ(272, 577), QQ(1121, 2378))]

    f = (x - 2)*(x - 1)*(2*x - 1)*(10002*x - 1)*(10003*x - 1)

    assert R._isolate_real_roots_sqf(f) == \
        [(QQ(15, 150046), QQ(47, 470110)), (QQ(47, 470110), QQ(17, 170018)),
         (QQ(1, 2), QQ(1, 2)), (1, 1), (2, 2)]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 100000000000)) == \
        [(QQ(1, 10003), QQ(1, 10003)), (QQ(1, 10002), QQ(1, 10002)),
         (QQ(1, 2), QQ(1, 2)), (1, 1), (2, 2)]

    a, b, c, d = 10000090000001, 2000100003, 10000300007, 10000005000008
    f = (x - d)*(x + a)*(b*x + 1)*(c*x - 1)

    assert R._isolate_real_roots_sqf(f) == \
        [(-13194139533313, -8796093022209), (-1, 0), (0, 1),
         (8796093022209, 13194139533313)]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 100000000000)) == \
        [(-a, -a), (QQ(-7, 13958643719), QQ(-1, 2013265921)),
         (QQ(3, 30064771075), QQ(1, 9663676417)),
         (QQ(1328823874562668133568119, 132882321015),
          QQ(37336367728494399224248237, 3733634906029))]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 100000000000000000000000000000)) == \
        [(-a, -a), (-QQ(1, b), -QQ(1, b)), (QQ(1, c), QQ(1, c)), (d, d)]

    f = -2*(x - 2)*(x + 2)*(5*x**2 - 4*x - 20)

    assert R._isolate_real_roots_sqf(f) == \
        [(-2, -2), (-2, -1), (2, 2), (2, 3)]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 100)) == \
        [(-2, -2), (-QQ(23, 14), -QQ(18, 11)), (2, 2), (QQ(39, 16), QQ(22, 9))]

    f = x - 1

    assert R._isolate_real_roots_sqf(f, inf=2) == []
    assert R._isolate_real_roots_sqf(f, sup=0) == []
    assert R._isolate_real_roots_sqf(f) == [(1, 1)]
    assert R._isolate_real_roots_sqf(f, inf=1) == [(1, 1)]
    assert R._isolate_real_roots_sqf(f, sup=1) == [(1, 1)]
    assert R._isolate_real_roots_sqf(f, inf=1, sup=1) == [(1, 1)]

    f = x**2 - 2

    assert R._isolate_real_roots_sqf(f, inf=QQ(7, 4)) == []
    assert R._isolate_real_roots_sqf(f, inf=QQ(7, 5)) == [(QQ(7, 5), QQ(3, 2))]
    assert R._isolate_real_roots_sqf(f, sup=QQ(7, 5)) == [(-2, -1)]
    assert R._isolate_real_roots_sqf(f, sup=QQ(7, 4)) == [(-2, -1), (1, QQ(3, 2))]
    assert R._isolate_real_roots_sqf(f, sup=-QQ(7, 4)) == []
    assert R._isolate_real_roots_sqf(f, sup=-QQ(7, 5)) == [(-QQ(3, 2), -QQ(7, 5))]
    assert R._isolate_real_roots_sqf(f, inf=-QQ(7, 5)) == [(1, 2)]
    assert R._isolate_real_roots_sqf(f, inf=-QQ(7, 4)) == [(-QQ(3, 2), -1), (1, 2)]

    i = [(-2, -1), (1, 2)]

    assert R._isolate_real_roots_sqf(f, inf=-2) == i
    assert R._isolate_real_roots_sqf(f, sup=+2) == i
    assert R._isolate_real_roots_sqf(f, inf=-2, sup=2) == i
    assert R._isolate_real_roots_sqf(f, inf=+1) == [i[1]]
    assert R._isolate_real_roots_sqf(f, sup=-1) == [i[0]]

    f = (2*x - 3)*(x**2 - 3)*(x**2 - 2)

    assert R._isolate_real_roots_sqf(f) == \
        [(-2, -QQ(3, 2)), (-QQ(3, 2), -QQ(1, 1)), (1, QQ(3, 2)),
         (QQ(3, 2), QQ(3, 2)), (QQ(3, 2), 2)]

    f = 7*x**4 - 19*x**3 + 20*x**2 + 17*x + 20

    assert R._isolate_real_roots_sqf(f) == []

    R, x = ring('x', QQ)

    f = (6*x - 85)*(1028*x + 1)/3855

    assert R._isolate_real_roots_sqf(f) == [(-1, 0), (14, 15)]
    assert [_.as_tuple() for _ in R._isolate_real_roots_sqf(f, blackbox=True)] == [(-1, 0), (14, 15)]

    f = (2*x/5 - QQ(17, 3))*(4*x + QQ(1, 257))

    assert R._isolate_real_roots_sqf(f) == [(-1, 0), (14, 15)]

    R, x = ring('x', EX)

    pytest.raises(DomainError, lambda: R._isolate_real_roots_sqf(x + 3))

    R, x = ring('x', QQ.algebraic_field(I))

    f = (x - 1)*(x**3 + I*x - 2)

    assert R._isolate_real_roots_sqf(f) == [(1, 1)]
    assert R._isolate_real_roots_sqf(f, sup=0) == []

    f = (x**2 - 2)*(x**3 - x + I)

    assert R._isolate_real_roots_sqf(f) == [(QQ(-3, 2), QQ(-4, 3)), (QQ(4, 3), QQ(3, 2))]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 10), inf=0) == [(QQ(7, 5), QQ(10, 7))]

    assert R._isolate_real_roots_sqf(x) == [(0, 0)]
    assert R._isolate_real_roots_sqf(x - 1) == [(1, 1)]
    assert R._isolate_real_roots_sqf(x - I) == []

    f = (x + I)*(x - 1)

    assert [_.as_tuple() for _ in R._isolate_real_roots_sqf(f, blackbox=True)] == [(1, 1)]

    R, x = ring('x', QQ.algebraic_field(sqrt(2)))

    f = (-x**3 + sqrt(2)*x - 1)*(x**2 + 1)

    assert R._isolate_real_roots_sqf(f) == [(-2, -1)]
    assert R._isolate_real_roots_sqf(f, eps=QQ(1, 1000)) == [(QQ(-132, 91), QQ(-29, 20))]

    f = (x - sqrt(2))*(x + 2*sqrt(2))*(x - 7 + sqrt(2))*(x + 3*sqrt(2))*(x - 1)*(x + 1 - sqrt(2))

    assert R._isolate_real_roots_sqf(f) == [(-5, -4), (-3, -2), (0, 1),
                                            (1, 1), (1, 2), (5, 6)]

    R, x = ring('x', QQ.algebraic_field(sqrt(2), sqrt(3)))

    f = (x - sqrt(2))*(x - sqrt(3))*(x - 2*sqrt(6))*(x - sqrt(6))*(x**2 + 2)

    assert R._isolate_real_roots_sqf(f) == [(1, QQ(3, 2)), (QQ(3, 2), 2),
                                            (2, 3), (4, 5)]
    assert (R._isolate_real_roots_sqf(f, eps=QQ(1, 1000)) ==
            [(QQ(41, 29), QQ(58, 41)), (QQ(71, 41), QQ(97, 56)),
             (QQ(218, 89), QQ(49, 20)), (QQ(436, 89), QQ(485, 99))])


def test__isolate_real_roots():
    R, x = ring('x', ZZ)

    assert R._isolate_real_roots(R(0)) == []
    assert R._isolate_real_roots(R(1)) == []
    assert R._isolate_real_roots(R(3)) == []

    assert R._isolate_real_roots(x) == [((0, 0), 1)]
    assert R._isolate_real_roots(5*x) == [((0, 0), 1)]
    assert R._isolate_real_roots(7*x**4) == [((0, 0), 4)]
    assert R._isolate_real_roots(x**128) == [((0, 0), 128)]

    assert R._isolate_real_roots(x*(x + 1)) == [((-1, -1), 1), ((0, 0), 1)]
    assert R._isolate_real_roots(x*(x - 1)) == [((0, 0), 1), ((1, 1), 1)]

    assert R._isolate_real_roots(x**4 + x + 1) == []

    i = [((-2, -1), 1), ((1, 2), 1)]

    assert R._isolate_real_roots(+x**2 - 2) == i
    assert R._isolate_real_roots(-x**2 + 2) == i

    f = (2*x - 3)**4*(x**2 - 3)**2*(x**2 - 2)**3

    assert R._isolate_real_roots(f) == \
        [((-2, -QQ(3, 2)), 2), ((-QQ(3, 2), -QQ(1, 1)), 3), ((1, QQ(3, 2)), 3),
         ((QQ(3, 2), QQ(3, 2)), 4), ((QQ(5, 3), 2), 2)]

    f = (2*x - 3)*(x**2 - 3)*(x**2 - 2)

    assert R._isolate_real_roots(f) == \
        [((-2, -QQ(3, 2)), 1), ((-QQ(3, 2), -QQ(1, 1)), 1), ((1, QQ(3, 2)), 1),
         ((QQ(3, 2), QQ(3, 2)), 1), ((QQ(3, 2), 2), 1)]

    f = x - 1

    assert R._isolate_real_roots(f, inf=2) == []
    assert R._isolate_real_roots(f, sup=0) == []
    assert R._isolate_real_roots(f) == [((1, 1), 1)]
    assert R._isolate_real_roots(f, inf=1) == [((1, 1), 1)]
    assert R._isolate_real_roots(f, sup=1) == [((1, 1), 1)]
    assert R._isolate_real_roots(f, inf=1, sup=1) == [((1, 1), 1)]

    f = (x**2 - 2)**2

    assert R._isolate_real_roots(f, inf=QQ(7, 4)) == []
    assert R._isolate_real_roots(f, inf=QQ(7, 5)) == [((QQ(7, 5), QQ(3, 2)), 2)]
    assert R._isolate_real_roots(f, sup=QQ(7, 5)) == [((-2, -1), 2)]
    assert R._isolate_real_roots(f, sup=QQ(7, 4)) == [((-2, -1), 2), ((1, QQ(3, 2)), 2)]
    assert R._isolate_real_roots(f, sup=-QQ(7, 4)) == []
    assert R._isolate_real_roots(f, sup=-QQ(7, 5)) == [((-QQ(3, 2), -QQ(7, 5)), 2)]
    assert R._isolate_real_roots(f, inf=-QQ(7, 5)) == [((1, 2), 2)]
    assert R._isolate_real_roots(f, inf=-QQ(7, 4)) == [((-QQ(3, 2), -1), 2), ((1, 2), 2)]

    i = [((-2, -1), 2), ((1, 2), 2)]

    assert R._isolate_real_roots(f, inf=-2) == i
    assert R._isolate_real_roots(f, sup=+2) == i
    assert R._isolate_real_roots(f, inf=-2, sup=2) == i

    f = x**4*(x - 1)**3*(x**2 - 2)**2

    assert R._isolate_real_roots(f) == \
        [((-2, -1), 2), ((0, 0), 4), ((1, 1), 3), ((1, 2), 2)]

    f = x**45 - 45*x**44 + 990*x**43 - 1
    g = (x**46 - 15180*x**43 + 9366819*x**40 - 53524680*x**39 +
         260932815*x**38 - 1101716330*x**37 + 4076350421*x**36 -
         13340783196*x**35 + 38910617655*x**34 - 101766230790*x**33 +
         239877544005*x**32 - 511738760544*x**31 + 991493848554*x**30 -
         1749695026860*x**29 + 2818953098830*x**28 - 4154246671960*x**27 +
         5608233007146*x**26 - 6943526580276*x**25 + 7890371113950*x**24 -
         8233430727600*x**23 + 7890371113950*x**22 - 6943526580276*x**21 +
         5608233007146*x**20 - 4154246671960*x**19 + 2818953098830*x**18 -
         1749695026860*x**17 + 991493848554*x**16 - 511738760544*x**15 +
         239877544005*x**14 - 101766230790*x**13 + 38910617655*x**12 -
         13340783196*x**11 + 4076350421*x**10 - 1101716330*x**9 +
         260932815*x**8 - 53524680*x**7 + 9366819*x**6 - 1370754*x**5 +
         163185*x**4 - 15180*x**3 + 1035*x**2 - 47*x + 1)

    assert R._isolate_real_roots(f*g) == \
        [((0, QQ(1, 2)), 1), ((QQ(2, 3), QQ(3, 4)), 1), ((QQ(3, 4), 1), 1), ((6, 7), 1), ((24, 25), 1)]

    f = x**2 - 3

    assert R._isolate_real_roots(f) == [((-2, -1), 1), ((1, 2), 1)]
    assert R._isolate_real_roots(f, eps=QQ(1, 100)) == [((QQ(-26, 15), QQ(-19, 11)), 1), ((QQ(19, 11), QQ(26, 15)), 1)]

    f = x**4 - 4*x**2 + 4

    assert R._isolate_real_roots(f, inf=QQ(7, 4)) == []
    assert R._isolate_real_roots(f, inf=QQ(7, 5)) == [((QQ(7, 5), QQ(3, 2)), 2)]
    assert R._isolate_real_roots(f, sup=QQ(7, 4)) == [((-2, -1), 2), ((1, QQ(3, 2)), 2)]
    assert R._isolate_real_roots(f, sup=QQ(7, 5)) == [((-2, -1), 2)]

    f = (x**2 - 2)*(x**2 - 3)**7*(x + 1)*(7*x + 3)**3

    assert R._isolate_real_roots(f) == [((-2, -QQ(3, 2)), 7), ((-QQ(3, 2), -1), 1),
                                        ((-1, -1), 1), ((-1, 0), 3),
                                        ((1, QQ(3, 2)), 1), ((QQ(3, 2), 2), 7)]

    f = 7*x**4 - 19*x**3 + 20*x**2 + 17*x + 20

    assert R._isolate_real_roots(f) == []

    R, x = ring('x', QQ)

    f = (2*x/5 - QQ(17, 3))*(4*x + QQ(1, 257))

    assert R._isolate_real_roots(f) == [((-1, 0), 1), ((14, 15), 1)]

    assert R._isolate_real_roots(f, eps=QQ(1, 10)) == [((-QQ(1, 513), 0), 1), ((QQ(85, 6), QQ(85, 6)), 1)]
    assert R._isolate_real_roots(f, eps=QQ(1, 100)) == [((-QQ(1, 513), 0), 1), ((QQ(85, 6), QQ(85, 6)), 1)]
    assert R._isolate_real_roots(f, eps=QQ(1, 1000)) == [((-QQ(1, 1025), 0), 1), ((QQ(85, 6), QQ(85, 6)), 1)]
    assert R._isolate_real_roots(f, eps=QQ(1, 10000)) == [((-QQ(1, 1025), -QQ(65, 66881)), 1), ((QQ(85, 6), QQ(85, 6)), 1)]

    R, x = ring('x', EX)

    pytest.raises(DomainError, lambda: R._isolate_real_roots(x + 3))
    pytest.raises(DomainError, lambda: R._isolate_real_roots((x + 2)*(x + 3)**2))

    R, x = ring('x', QQ.algebraic_field(I))

    f = (x**2 - I)**2*(x - 2*I)**3

    assert R._isolate_real_roots(f) == []  # issue diofant/diofant#789
    assert R._isolate_real_roots(f*(x - 1)**3) == [((1, 1), 3)]

    f = x**4*(x - 1)**3*(x**2 - 2)**2

    assert R._isolate_real_roots(f) == \
        [((-2, -1), 2), ((0, 0), 4), ((1, 1), 3), ((QQ(4, 3), QQ(3, 2)), 2)]


def test__isolate_real_roots_pair():
    R, x = ring('x', ZZ)

    assert R._isolate_real_roots_pair(x*(x + 1), x) == \
        [((-1, -1), {0: 1}), ((0, 0), {0: 1, 1: 1})]
    assert R._isolate_real_roots_pair(x*(x - 1), x) == \
        [((0, 0), {0: 1, 1: 1}), ((1, 1), {0: 1})]

    f, g = (x**2 - 2)**2, x - 1

    assert R._isolate_real_roots_pair(f, g, inf=QQ(7, 4)) == []
    assert R._isolate_real_roots_pair(f, g, inf=QQ(7, 5)) == \
        [((QQ(7, 5), QQ(3, 2)), {0: 2})]
    assert R._isolate_real_roots_pair(f, g, sup=QQ(7, 5)) == \
        [((-2, -1), {0: 2}), ((1, 1), {1: 1})]
    assert R._isolate_real_roots_pair(f, g, sup=QQ(7, 4)) == \
        [((-2, -1), {0: 2}), ((1, 1), {1: 1}), ((1, QQ(3, 2)), {0: 2})]
    assert R._isolate_real_roots_pair(f, g, sup=-QQ(7, 4)) == []
    assert R._isolate_real_roots_pair(f, g, sup=-QQ(7, 5)) == \
        [((-QQ(3, 2), -QQ(7, 5)), {0: 2})]
    assert R._isolate_real_roots_pair(f, g, inf=-QQ(7, 5)) == \
        [((1, 1), {1: 1}), ((1, 2), {0: 2})]
    assert R._isolate_real_roots_pair(f, g, inf=-QQ(7, 4)) == \
        [((-QQ(3, 2), -1), {0: 2}), ((1, 1), {1: 1}), ((1, 2), {0: 2})]

    f, g = 2*x**2 - 1, x**2 - 2

    assert R._isolate_real_roots_pair(f, g) == \
        [((-2, -1), {1: 1}), ((-1, 0), {0: 1}),
         ((0, 1), {0: 1}), ((1, 2), {1: 1})]
    assert R._isolate_real_roots_pair(f, g, strict=True) == \
        [((-QQ(3, 2), -QQ(4, 3)), {1: 1}), ((-1, -QQ(2, 3)), {0: 1}),
         ((QQ(2, 3), 1), {0: 1}), ((QQ(4, 3), QQ(3, 2)), {1: 1})]

    f, g = x**2 - 2, (x - 1)*(x**2 - 2)

    assert R._isolate_real_roots_pair(f, g) == \
        [((-2, -1), {1: 1, 0: 1}), ((1, 1), {1: 1}), ((1, 2), {1: 1, 0: 1})]

    f, g = x*(x**2 - 2), x**2*(x - 1)*(x**2 - 2)

    assert R._isolate_real_roots_pair(f, g) == \
        [((-2, -1), {1: 1, 0: 1}), ((0, 0), {0: 1, 1: 2}),
         ((1, 1), {1: 1}), ((1, 2), {1: 1, 0: 1})]

    f, g = x**2*(x - 1)**3*(x**2 - 2)**2, x*(x - 1)**2*(x**2 + 2)
    _x = R.clone(domain=ZZ.field).x

    assert R._isolate_real_roots_pair(f, g) == \
        [((-2, -1), {0: 2}), ((0, 0), {0: 2, 1: 1}),
         ((1, 1), {0: 3, 1: 2}), ((1, 2), {0: 2})]
    assert R._isolate_real_roots_pair(f, g, basis=True) == \
        [((-2, -1), {0: 2}, _x**2 - 2), ((0, 0), {0: 2, 1: 1}, _x),
         ((1, 1), {0: 3, 1: 2}, _x - 1), ((1, 2), {0: 2}, _x**2 - 2)]

    f, g = x, R.zero

    assert R._isolate_real_roots_pair(f, g) == \
        R._isolate_real_roots_pair(g, f) == [((0, 0), {0: 1, 1: 1})]

    f *= x**2

    assert R._isolate_real_roots_pair(f, g) == \
        R._isolate_real_roots_pair(g, f) == [((0, 0), {0: 3, 1: 3})]

    R, x = ring('x', EX)

    pytest.raises(DomainError, lambda: R._isolate_real_roots_pair(x, x + 3))

    R, x = ring('x', ZZ)

    f, g = x**5 - 200, x**5 - 201

    assert R._isolate_real_roots_pair(f, g) == \
        [((QQ(75, 26), QQ(101, 35)), {0: 1}), ((QQ(309, 107), QQ(26, 9)), {1: 1})]

    R, x = ring('x', QQ)

    f, g = -x**5/200 + 1, -x**5/201 + 1

    assert R._isolate_real_roots_pair(f, g) == \
        [((QQ(75, 26), QQ(101, 35)), {0: 1}), ((QQ(309, 107), QQ(26, 9)), {1: 1})]


def test__count_real_roots():
    R, x = ring('x', ZZ)

    assert R._count_real_roots(R(0)) == 0
    assert R._count_real_roots(R(7)) == 0

    f = x - 1

    assert R._count_real_roots(f) == 1
    assert R._count_real_roots(f, inf=1) == 1
    assert R._count_real_roots(f, sup=0) == 0
    assert R._count_real_roots(f, sup=1) == 1
    assert R._count_real_roots(f, inf=0, sup=1) == 1
    assert R._count_real_roots(f, inf=0, sup=2) == 1
    assert R._count_real_roots(f, inf=1, sup=2) == 1

    f = x**2 - 2

    assert R._count_real_roots(f) == 2
    assert R._count_real_roots(f, sup=0) == 1
    assert R._count_real_roots(f, inf=-1, sup=1) == 0

    R, x = ring('x', QQ.algebraic_field(I))

    f = x**3 + I*x + 2

    assert R._count_real_roots(f) == 0

    f *= (x - 1)*(x + 1)

    assert R._count_real_roots(f) == 2


# parameters for test_dup_count_complex_roots_n(): n = 1..8
a, b = (-1, -1), (1, 1)
c, d = (+0, +0), (1, 1)


def test__count_complex_roots_1():
    R, x = ring('x', ZZ)

    f = x - 1

    assert R._count_complex_roots(f, a, b) == 1
    assert R._count_complex_roots(f, c, d) == 1

    f = -f

    assert R._count_complex_roots(f, a, b) == 1
    assert R._count_complex_roots(f, c, d) == 1

    f = x + 1

    assert R._count_complex_roots(f, a, b) == 1
    assert R._count_complex_roots(f, c, d) == 0

    R, x = ring('x', QQ)

    f = x - QQ(1, 2)

    assert R._count_complex_roots(f, c, d) == 1

    R, x = ring('x', EX)

    pytest.raises(DomainError, lambda: R._count_complex_roots(x))


def test__count_complex_roots_2():
    R, x = ring('x', ZZ)

    f = x*(x - 1)

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 2

    f = -f

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 2

    f = x*(x + 1)

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 1

    f = -f

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 1


def test__count_complex_roots_3():
    R, x = ring('x', ZZ)

    f = (x - 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 1

    f = x*(x - 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 2

    f = -f

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 2


def test__count_complex_roots_4():
    R, x = ring('x', ZZ)

    f = x**2 + 1

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 1

    f = x*(x**2 + 1)

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 2

    f = -f

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 2

    f = (x**2 + 1)*(x - 1)

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 2

    f *= x

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 3

    f = -f

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 3

    f = (x**2 + 1)*(x - 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 2

    f *= x

    assert R._count_complex_roots(f, a, b) == 5
    assert R._count_complex_roots(f, c, d) == 3

    f = -f

    assert R._count_complex_roots(f, a, b) == 5
    assert R._count_complex_roots(f, c, d) == 3


def test__count_complex_roots_5():
    R, x = ring('x', ZZ)

    f = (x + 1)**2 + 1

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 0

    f = ((x + 1)**2 + 1)*(x - 1)

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 1

    f *= x

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 2

    f = ((x + 1)**2 + 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 0

    f *= x

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 1

    f = ((x + 1)**2 + 1)*(x - 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 1

    f *= x

    assert R._count_complex_roots(f, a, b) == 5
    assert R._count_complex_roots(f, c, d) == 2


def test__count_complex_roots_6():
    R, x = ring('x', ZZ)

    f = (x - 1)**2 + 1

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 1

    f *= x - 1

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 2

    f *= x

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 3

    f = ((x - 1)**2 + 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 3
    assert R._count_complex_roots(f, c, d) == 1

    f *= x
    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 2

    f = ((x - 1)**2 + 1)*(x - 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 2

    f *= x

    assert R._count_complex_roots(f, a, b) == 5
    assert R._count_complex_roots(f, c, d) == 3


def test__count_complex_roots_7():
    R, x = ring('x', ZZ)

    f = ((x - 1)**2 + 1)*((x + 1)**2 + 1)

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 1

    f *= (x - 2)

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 1

    f = ((x - 1)**2 + 1)*((x + 1)**2 + 1)*(x**2 - 2)

    assert R._count_complex_roots(f, a, b) == 4
    assert R._count_complex_roots(f, c, d) == 1

    f = ((x - 1)**2 + 1)*((x + 1)**2 + 1)*(x - 1)

    assert R._count_complex_roots(f, a, b) == 5
    assert R._count_complex_roots(f, c, d) == 2

    f *= x

    assert R._count_complex_roots(f, a, b) == 6
    assert R._count_complex_roots(f, c, d) == 3

    f = ((x - 1)**2 + 1)*((x + 1)**2 + 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 5
    assert R._count_complex_roots(f, c, d) == 1

    f *= x

    assert R._count_complex_roots(f, a, b) == 6
    assert R._count_complex_roots(f, c, d) == 2

    f = ((x - 1)**2 + 1)*((x + 1)**2 + 1)*(x - 1)*(x + 1)

    assert R._count_complex_roots(f, a, b) == 6
    assert R._count_complex_roots(f, c, d) == 2

    f *= x

    assert R._count_complex_roots(f, a, b) == 7
    assert R._count_complex_roots(f, c, d) == 3

    f = ((x - 1)**2 + 1)*((x + 1)**2 + 1)*(x - 1)*(x + 1)*(x**2 + 1)

    assert R._count_complex_roots(f, a, b) == 8
    assert R._count_complex_roots(f, c, d) == 3


def test__count_complex_roots_8():
    R, x = ring('x', ZZ)

    f = ((x - 1)**2 + 1)*((x + 1)**2 + 1)*(x - 1)*(x + 1)*(x**2 + 1)*x

    assert R._count_complex_roots(f, a, b) == 9
    assert R._count_complex_roots(f, c, d) == 4

    f *= (x**2 - 2)

    assert R._count_complex_roots(f, a, b) == 9
    assert R._count_complex_roots(f, c, d) == 4


def test__count_complex_roots_9():
    R, x = ring('x', QQ.algebraic_field(sqrt(2)))

    f = -x**3 + sqrt(2)*x - 1

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 1

    R, x = ring('x', QQ.algebraic_field(sqrt(2)).algebraic_field(I))

    f = -x**3 + I*x**2 + sqrt(2)*x - 1

    assert R._count_complex_roots(f, a, b) == 2
    assert R._count_complex_roots(f, c, d) == 1


def test__count_complex_roots_implicit():
    R, x = ring('x', ZZ)

    f = (x**2 + 1)*(x - 1)*(x + 1)*x

    assert R._count_complex_roots(f) == 5

    assert R._count_complex_roots(f, sup=(0, 0)) == 3
    assert R._count_complex_roots(f, inf=(0, 0)) == 3

    assert R._count_complex_roots(f, inf=QQ(-2), sup=QQ(-1)) == 1


def test__count_complex_roots_exclude():
    R, x = ring('x', ZZ)

    f = (x**2 + 1)*(x - 1)*(x + 1)*x

    a, b = (-1, 0), (1, 1)

    assert R._count_complex_roots(f, a, b) == 4

    assert R._count_complex_roots(f, a, b, exclude=['S']) == 3
    assert R._count_complex_roots(f, a, b, exclude=['N']) == 3

    assert R._count_complex_roots(f, a, b, exclude=['S', 'N']) == 2

    assert R._count_complex_roots(f, a, b, exclude=['E']) == 4
    assert R._count_complex_roots(f, a, b, exclude=['W']) == 4

    assert R._count_complex_roots(f, a, b, exclude=['E', 'W']) == 4

    assert R._count_complex_roots(f, a, b, exclude=['N', 'S', 'E', 'W']) == 2

    assert R._count_complex_roots(f, a, b, exclude=['SW']) == 3
    assert R._count_complex_roots(f, a, b, exclude=['SE']) == 3

    assert R._count_complex_roots(f, a, b, exclude=['SW', 'SE']) == 2
    assert R._count_complex_roots(f, a, b, exclude=['SW', 'SE', 'S']) == 1
    assert R._count_complex_roots(f, a, b, exclude=['SW', 'SE', 'S', 'N']) == 0

    a, b = (0, 0), (1, 1)

    assert R._count_complex_roots(f, a, b, exclude=True) == 1

    R, x = ring('x', QQ.algebraic_field(I))

    f = x**4 + I*x**3 - x + 1

    assert R._count_complex_roots(f, inf=(0, 0), sup=(1, 1)) == 1

    r = R._isolate_complex_roots_sqf(f)

    assert r == [((QQ(-201, 100), QQ(-201, 100)), (0, 0)),
                 ((QQ(-201, 100), 0), (0, QQ(201, 100))),
                 ((0, QQ(-201, 100)), (QQ(201, 100), 0)),
                 ((0, 0), (QQ(201, 100), QQ(201, 100)))]
    assert all(R._count_complex_roots(f, inf=i, sup=s) == 1
               for i, s in r)


def test__isolate_complex_roots_sqf():
    R, x = ring('x', ZZ)

    f = x**2 - 2*x + 3

    assert R._isolate_complex_roots_sqf(f) == \
        [((0, -6), (6, 0)), ((0, 0), (6, 6))]
    assert [r.as_tuple() for r in R._isolate_complex_roots_sqf(f, blackbox=True)] == \
        [((0, -6), (6, 0)), ((0, 0), (6, 6))]

    assert R._isolate_complex_roots_sqf(f, inf=1, sup=3) == [((1, -3), (3, 0)), ((1, 0), (3, 3))]
    assert R._isolate_complex_roots_sqf(f, inf=(1, 0), sup=3) == [((1, 0), (3, 3))]
    assert R._isolate_complex_roots_sqf(f, inf=(1, QQ(-1, 2)), sup=3) == [((1, 0), (3, 3))]
    assert R._isolate_complex_roots_sqf(f, inf=(1, -3), sup=(3, -1)) == [((1, -3), (3, -1))]
    assert R._isolate_complex_roots_sqf(f, inf=0, sup=QQ(1, 6)) == []

    assert R._isolate_complex_roots_sqf(R.zero) == []

    pytest.raises(ValueError, lambda: R._isolate_complex_roots_sqf(f, inf=1, sup=1))

    assert R._isolate_complex_roots_sqf(f, eps=QQ(1, 10)) == \
        [((QQ(15, 16), -QQ(3, 2)), (QQ(33, 32), -QQ(45, 32))),
         ((QQ(15, 16), QQ(45, 32)), (QQ(33, 32), QQ(3, 2)))]
    assert R._isolate_complex_roots_sqf(f, eps=QQ(1, 100)) == \
        [((QQ(255, 256), -QQ(363, 256)), (QQ(513, 512), -QQ(723, 512))),
         ((QQ(255, 256), QQ(723, 512)), (QQ(513, 512), QQ(363, 256)))]

    f = 7*x**4 - 19*x**3 + 20*x**2 + 17*x + 20

    assert R._isolate_complex_roots_sqf(f) == \
        [((-QQ(40, 7), -QQ(40, 7)), (0, 0)), ((-QQ(40, 7), 0), (0, QQ(40, 7))),
         ((0, -QQ(40, 7)), (QQ(40, 7), 0)), ((0, 0), (QQ(40, 7), QQ(40, 7)))]
    assert R._isolate_complex_roots_sqf(f, eps=QQ(1, 10)) == \
        [((QQ(-25, 56), QQ(-5, 8)), (QQ(-5, 14), QQ(-15, 28))),
         ((QQ(-25, 56), QQ(15, 28)), (QQ(-5, 14), QQ(5, 8))),
         ((QQ(95, 56), QQ(-85, 56)), (QQ(25, 14), QQ(-10, 7))),
         ((QQ(95, 56), QQ(10, 7)), (QQ(25, 14), QQ(85, 56)))]

    R, x = ring('x', QQ)

    f = x**2/2 - 3*x/7 + 1

    assert R._isolate_complex_roots_sqf(f) == [((0, -4), (4, 0)), ((0, 0), (4, 4))]

    R, x = ring('x', EX)

    pytest.raises(DomainError,
                  lambda: R._isolate_complex_roots_sqf(x, inf=(-1, 0),
                                                       sup=(1, 1)))

    R, x = ring('x', QQ.algebraic_field(I))

    f = x**4 + I*x**3 - x + 1

    assert R._isolate_complex_roots_sqf(f, inf=(0, 0),
                                        sup=(1, 1)) == [((0, 0), (1, QQ(1, 2)))]
    assert R._isolate_complex_roots_sqf(f, inf=(0, 0), sup=(1, 1),
                                        eps=QQ(1, 100)) == [((QQ(79, 128), QQ(19, 64)),
                                                            (QQ(5, 8), QQ(39, 128)))]
    assert R._isolate_complex_roots_sqf(f, inf=(0, -1),
                                        sup=(1, 1)) == [((0, -1), (1, QQ(-1, 2))),
                                                        ((0, 0), (1, QQ(1, 2)))]
    assert R._isolate_complex_roots_sqf(f, inf=(0, -1), sup=(1, 1),
                                        eps=QQ(1, 100)) == [((QQ(79, 128), QQ(19, 64)),
                                                             (QQ(5, 8), QQ(39, 128))),
                                                            ((QQ(45, 64), QQ(-91, 128)),
                                                             (QQ(91, 128), QQ(-45, 64)))]

    f *= (x - 1)

    assert R._isolate_complex_roots_sqf(f) == [((QQ(-401, 100), QQ(-401, 100)), (0, 0)),
                                               ((QQ(-401, 100), 0), (0, QQ(401, 100))),
                                               ((0, QQ(-401, 100)), (QQ(401, 100), 0)),
                                               ((0, 0), (QQ(401, 100), QQ(401, 100)))]

    f = x**7 + I*x**4 - (2 + I)*x**3 - 3*x + 5

    assert R._isolate_complex_roots_sqf(f) == [((QQ(-1001, 100), 0), (0, QQ(1001, 100))),
                                               ((QQ(-1001, 400), QQ(-1001, 800)), (QQ(-1001, 800), 0)),
                                               ((QQ(-1001, 800), QQ(-1001, 800)), (0, 0)),
                                               ((0, QQ(-1001, 400)), (QQ(1001, 400), QQ(-1001, 800))),
                                               ((0, QQ(-1001, 800)), (QQ(1001, 400), 0)),
                                               ((0, 0), (QQ(1001, 400), QQ(1001, 800))),
                                               ((0, QQ(1001, 800)), (QQ(1001, 400), QQ(1001, 400)))]

    R, x = ring('x', QQ.algebraic_field(sqrt(2)))

    f = -x**3 + sqrt(2)*x - 1

    assert R._isolate_complex_roots_sqf(f) == [((0, QQ(-283, 100)), (QQ(283, 100), 0)),
                                               ((0, 0), (QQ(283, 100), QQ(283, 100)))]

    R, x = ring('x', QQ.algebraic_field(sqrt(2)).algebraic_field(I))

    f = -x**3 + I*x**2 + sqrt(2)*x - 1

    assert R._isolate_complex_roots_sqf(f) == [((QQ(-283, 100), 0), (0, QQ(283, 100))),
                                               ((0, QQ(-283, 100)), (QQ(283, 100), 0)),
                                               ((0, 0), (QQ(283, 100), QQ(283, 100)))]

    R, x = ring('x', EX)

    pytest.raises(DomainError, lambda: R._isolate_complex_roots_sqf(x))


@pytest.mark.timeout(300)
@pytest.mark.slow
@pytest.mark.skipif(isinstance(ZZ(42), int), reason='gmpy2 is not used')
def test__isolate_complex_roots_sqf_2():
    R, x = ring('x', ZZ)

    f = x**40 - 15*x**17 - 21*x**3 + 11

    res = R._isolate_complex_roots_sqf(f)
    ans = [((QQ(-21, 16), QQ(-21, 128)), (QQ(-63, 64), 0)),
           ((QQ(-21, 16), 0), (QQ(-63, 64), QQ(21, 128))),
           ((QQ(-21, 16), QQ(-21, 64)), (QQ(-63, 64), QQ(-21, 128))),
           ((QQ(-21, 16), QQ(21, 128)), (QQ(-63, 64), QQ(21, 64))),
           ((QQ(-21, 16), QQ(-21, 32)), (QQ(-63, 64), QQ(-21, 64))),
           ((QQ(-21, 16), QQ(21, 64)), (QQ(-63, 64), QQ(21, 32))),
           ((QQ(-63, 64), QQ(-21, 32)), (QQ(-21, 32), QQ(-21, 64))),
           ((QQ(-63, 64), QQ(21, 64)), (QQ(-21, 32), QQ(21, 32))),
           ((QQ(-63, 64), QQ(-105, 128)), (QQ(-21, 32), QQ(-21, 32))),
           ((QQ(-63, 64), QQ(21, 32)), (QQ(-21, 32), QQ(105, 128))),
           ((QQ(-63, 64), QQ(-63, 64)), (QQ(-21, 32), QQ(-105, 128))),
           ((QQ(-63, 64), QQ(105, 128)), (QQ(-21, 32), QQ(63, 64))),
           ((QQ(-21, 32), QQ(-105, 128)), (QQ(-21, 64), QQ(-21, 32))),
           ((QQ(-21, 32), QQ(21, 32)), (QQ(-21, 64), QQ(105, 128))),
           ((QQ(-21, 32), QQ(-63, 64)), (QQ(-21, 64), QQ(-105, 128))),
           ((QQ(-21, 32), QQ(105, 128)), (QQ(-21, 64), QQ(63, 64))),
           ((QQ(-21, 32), QQ(-21, 16)), (QQ(-21, 64), QQ(-63, 64))),
           ((QQ(-21, 32), QQ(63, 64)), (QQ(-21, 64), QQ(21, 16))),
           ((QQ(-21, 64), QQ(-21, 16)), (QQ(0, 1), QQ(-63, 64))),
           ((QQ(-21, 64), QQ(63, 64)), (QQ(0, 1), QQ(21, 16))),
           ((QQ(0, 1), QQ(-147, 128)), (QQ(21, 128), QQ(-63, 64))),
           ((QQ(0, 1), QQ(63, 64)), (QQ(21, 128), QQ(147, 128))),
           ((QQ(21, 128), QQ(-147, 128)), (QQ(21, 64), QQ(-63, 64))),
           ((QQ(21, 128), QQ(63, 64)), (QQ(21, 64), QQ(147, 128))),
           ((QQ(21, 64), QQ(-63, 64)), (QQ(63, 128), QQ(-105, 128))),
           ((QQ(21, 64), QQ(105, 128)), (QQ(63, 128), QQ(63, 64))),
           ((QQ(63, 128), QQ(-63, 64)), (QQ(21, 32), QQ(-105, 128))),
           ((QQ(63, 128), QQ(105, 128)), (QQ(21, 32), QQ(63, 64))),
           ((QQ(21, 32), QQ(-21, 64)), (QQ(63, 64), 0)),
           ((QQ(21, 32), 0), (QQ(63, 64), QQ(21, 64))),
           ((QQ(21, 32), QQ(-21, 32)), (QQ(21, 16), QQ(-21, 64))),
           ((QQ(21, 32), QQ(21, 64)), (QQ(21, 16), QQ(21, 32))),
           ((QQ(21, 32), QQ(-105, 128)), (QQ(63, 64), QQ(-21, 32))),
           ((QQ(21, 32), QQ(21, 32)), (QQ(63, 64), QQ(105, 128))),
           ((QQ(21, 32), QQ(-63, 64)), (QQ(63, 64), QQ(-105, 128))),
           ((QQ(21, 32), QQ(105, 128)), (QQ(63, 64), QQ(63, 64))),
           ((QQ(63, 64), QQ(-21, 64)), (QQ(21, 16), 0)),
           ((QQ(63, 64), 0), (QQ(21, 16), QQ(21, 64)))]

    assert res == ans


def test__isolate_all_roots_sqf():
    R, x = ring('x', ZZ)

    f = (4*x**3 - x**2 + 2*x + 5)*x

    assert R._isolate_all_roots_sqf(f) == \
        ([(-1, 0), (0, 0)],
         [((0, -QQ(5, 2)), (QQ(5, 2), 0)), ((0, 0), (QQ(5, 2), QQ(5, 2)))])

    assert R._isolate_all_roots_sqf(f, eps=QQ(1, 10)) == \
        ([(QQ(-7, 8), QQ(-6, 7)), (0, 0)],
         [((QQ(35, 64), -QQ(35, 32)), (QQ(5, 8), -QQ(65, 64))), ((QQ(35, 64), QQ(65, 64)), (QQ(5, 8), QQ(35, 32)))])

    R, x = ring('x', EX)

    pytest.raises(DomainError, lambda: R._isolate_all_roots_sqf(x, R))


def test__isolate_all_roots():
    R, x = ring('x', ZZ)

    f = (4*x**3 - x**2 + 2*x + 5)*x

    assert R._isolate_all_roots(f) == \
        ([((-1, 0), 1), ((0, 0), 1)],
         [(((0, -QQ(5, 2)), (QQ(5, 2), 0)), 1),
          (((0, 0), (QQ(5, 2), QQ(5, 2))), 1)])

    assert R._isolate_all_roots(f, eps=QQ(1, 10)) == \
        ([((QQ(-7, 8), QQ(-6, 7)), 1), ((0, 0), 1)],
         [(((QQ(35, 64), -QQ(35, 32)), (QQ(5, 8), -QQ(65, 64))), 1),
          (((QQ(35, 64), QQ(65, 64)), (QQ(5, 8), QQ(35, 32))), 1)])

    f = (x - 1)**2*(x + 1)**3

    pytest.raises(NotImplementedError, lambda: R._isolate_all_roots(f))

    D, y = ring('y', ZZ)
    R, x = ring('x', D)

    f = x**2 + y*x - 1

    pytest.raises(DomainError, lambda: R._isolate_all_roots(f))


def test_RealInterval():
    R, x = ring('x', ZZ)

    f = (x - 1)**2

    pytest.raises(ValueError, lambda: RealInterval((-2, 1), f))


def test_ComplexInterval():
    R, x = ring('x', QQ.algebraic_field(I))

    f = x**3 + x + I

    _, r1, r2 = R._isolate_complex_roots_sqf(f, blackbox=True)

    assert r1.is_disjoint(r2) is True
    assert r1.is_disjoint(r2, check_re_refinement=True) is False

    for i in range(4):
        r1, r2 = r1.refine(), r2.refine()

    assert r1.is_disjoint(r2, check_re_refinement=True) is True

    (u1, v1), (s1, t1) = r1.as_tuple()
    (u2, v2), (s2, t2) = r1.refine(vertical=True).as_tuple()

    assert v1 == v2 and t1 == t2
    assert u1 <= u2 < s2 < s1


def test_diofantissue_745():
    D, y = ring('y', ZZ)
    R, x = ring('x', D)

    pytest.raises(DomainError, lambda: R._count_real_roots(x**7 + y*x + 1))
