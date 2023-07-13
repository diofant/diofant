"""Tests for modular GCD algorithms."""

import pytest

from diofant import QQ, ZZ, ring, sqrt
from diofant.polys.modulargcd import (_chinese_remainder_reconstruction,
                                      _func_field_modgcd_m, _to_ANP_poly,
                                      _to_ZZ_poly)


__all__ = ()


def test_chinese_remainder():
    R, x, y = ring('x y', ZZ)
    p, q = 3, 5

    hp = x**3*y - x**2 - 1
    hq = -x**3*y - 2*x*y**2 + 2

    hpq = _chinese_remainder_reconstruction(hp, hq, p, q)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq

    _, z = ring('z', R)
    p, q = 3, 7

    hp = (x*y + 1)*z**2 + x
    hq = (x**2 - 3*y)*z + 2

    hpq = _chinese_remainder_reconstruction(hp, hq, p, q)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq


def test_to_ZZ_ANP_poly():
    A = QQ.algebraic_field(sqrt(2))
    R, x = ring('x', A)
    f = x*(sqrt(2) + 1)

    T, x_, z_ = ring('x_ z_', ZZ)
    f_ = x_*z_ + x_

    assert _to_ZZ_poly(f, T) == f_
    assert _to_ANP_poly(f_, R) == f

    R, x, t, s = ring('x t s', A)
    f = x*t**2 + x*s + sqrt(2)

    D, t_, s_ = ring('t_ s_', ZZ)
    T, x_, z_ = ring('x_ z_', D)
    f_ = (t_**2 + s_)*x_ + z_

    assert _to_ZZ_poly(f, T) == f_
    assert _to_ANP_poly(f_, R) == f


def test_modgcd_func_field():
    D, t = ring('t', ZZ)
    _, x, z = ring('x z', D)

    minpoly = (z**2*t**2 + z**2*t - 1).drop(0)
    f, g = x + 1, x - 1

    assert _func_field_modgcd_m(f, g, minpoly) == 1

    # First example from Monagan2004algebraic.
    m = z**2 - t
    f = 3*t*x**2 - (2*t**2 - 3*t)*x*z + 15*x + 15*z - 2*t**3
    g = 3*t*x**2*z + 15*x*z + (-2*t**3 + 3*t)*x - 2*t**2*z + 15

    assert _func_field_modgcd_m(f, g, m.drop(0)) == 3*t*x - 2*t**2*z + 15

    g = 3*t*x - 2*t**2*z + 15
    a = x + z
    b = x*z + 1

    assert _func_field_modgcd_m(a*g, b*g, m.drop(0)) == g % m
    assert _func_field_modgcd_m(a*g**2, b*g**2, m.drop(0)) == g**2 % m

    # issue diofant/diofant#850
    assert _func_field_modgcd_m(a*g**3, b*g**3, m.drop(0)) == g**3 % m


@pytest.mark.slow
def test_modgcd_func_field_slow():
    D, t = ring('t', ZZ)
    _, x, z = ring('x z', D)

    # Example from sec.4 of Monagan2004algebraic.
    m = (z**3 - (5 - t)*z**2 + (7 - t**2)*z - (9 - t**3))
    minpoly = m.drop(0)
    g = ((10 - 4*t)*x**2 - 5*x*z**2 + (4*t + 1)*x*z + (11 - 17*t**2 + 9*t)*x -
         19*z**2 + (-7*t + 6)*z + (-11*t**2 + 15*t + 3))
    a = ((18 + 10*t)*x**2 + 10*x*z**2 + (17*t + 2)*x*z +
         (2 + 17*t**2 + 8*t)*x + 6*z**2 + (17*t + 6)*z + (4*t**2 - 4*t + 2))
    b = ((-8 - 11*t)*x**2 - 14*x*z**2 + (8*t - 4)*x*z +
         (-17 - 5*t**2 + 19*t)*x - 11*z**2 + (17*t - 4)*z + (-14*t**2 - 19*t - 2))

    n = 4
    for k in range(n):
        f1 = g**k*a**(n - k)
        f2 = g**k*b**(n - k)
        assert _func_field_modgcd_m(f1, f2, minpoly) == (-g)**k % m
