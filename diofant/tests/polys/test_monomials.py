"""Tests for tools and arithmetics for monomials of distributed polynomials."""

import functools

import pytest

from diofant import Monomial, itermonomials
from diofant.abc import a, b, c, x, y, z


__all__ = ()


def test_monomials():
    assert set(itermonomials([], 0)) == {1}
    assert set(itermonomials([], 1)) == {1}
    assert set(itermonomials([], 2)) == {1}
    assert set(itermonomials([], 3)) == {1}

    assert set(itermonomials([x], 0)) == {1}
    assert set(itermonomials([x], 1)) == {1, x}
    assert set(itermonomials([x], 2)) == {1, x, x**2}
    assert set(itermonomials([x], 3)) == {1, x, x**2, x**3}

    assert set(itermonomials([x, y], 0)) == {1}
    assert set(itermonomials([x, y], 1)) == {1, x, y}
    assert set(itermonomials([x, y], 2)) == {1, x, y, x**2, y**2, x*y}
    assert set(itermonomials([x, y], 3)) == {1, x, y, x**2, x**3, y**2,
                                             y**3, x*y, x*y**2, y*x**2}


def test_monomial_ops():
    m1 = Monomial((3, 4, 1))
    assert m1*(1, 2, 0) == (4, 6, 1)
    assert m1/(1, 2, 0) == (2, 2, 1)
    assert m1**2 == (6, 8, 2)

    assert m1.gcd((1, 2, 0)) == (1, 2, 0)
    m2 = Monomial((1, 4, 1))
    assert m2.gcd((3, 2, 0)) == (1, 2, 0)
    m3 = Monomial((3, 4, 5))
    assert functools.reduce(Monomial.gcd, (m3, (0, 5, 1),
                                           (6, 3, 9))) == (0, 3, 1)

    assert m1.lcm((1, 2, 0)) == (3, 4, 1)
    m4 = Monomial((1, 4, 1))
    assert m4.lcm((3, 2, 0)) == (3, 4, 1)

    m5 = Monomial((1, 2, 3))
    assert m5.divides((4, 5, 6)) is True
    assert m5.divides((0, 5, 6)) is False
    m6 = Monomial((1, 2))
    assert m6.divides((3, 4)) is True
    assert m6.divides((0, 2)) is False

    pytest.raises(TypeError, lambda: m6.divides(2*x))


def test_Monomial():
    m = Monomial((3, 4, 1), (x, y, z))
    l = Monomial((3, 4, 1))
    n = Monomial((1, 2, 0), (x, y, z))

    assert m.as_expr() == x**3*y**4*z
    assert n.as_expr() == x**1*y**2

    assert m.as_expr(a, b, c) == a**3*b**4*c
    assert n.as_expr(a, b, c) == a**1*b**2

    pytest.raises(ValueError, lambda: l.as_expr())

    assert tuple(m) == (3, 4, 1)
    assert m.gens == (x, y, z)

    assert tuple(n) == (1, 2, 0)
    assert n.gens == (x, y, z)

    assert m == (3, 4, 1)
    assert n != (3, 4, 1)
    assert m != (1, 2, 0)
    assert n == (1, 2, 0)
    assert n != object()

    assert m != n
    assert hash(m) != hash(n)

    assert m[0] == m[-3] == 3
    assert m[1] == m[-2] == 4
    assert m[2] == m[-1] == 1

    assert n[0] == n[-3] == 1
    assert n[1] == n[-2] == 2
    assert n[2] == n[-1] == 0

    assert m[:2] == (3, 4)
    assert n[:2] == (1, 2)

    assert m*n == Monomial((4, 6, 1))
    assert m/n == Monomial((2, 2, 1))

    assert m*(1, 2, 0) == Monomial((4, 6, 1))
    assert m/(1, 2, 0) == Monomial((2, 2, 1))

    pytest.raises(TypeError, lambda: m*object())
    pytest.raises(TypeError, lambda: m/object())

    assert m.gcd(n) == Monomial((1, 2, 0))
    assert m.lcm(n) == Monomial((3, 4, 1))

    assert m.gcd((1, 2, 0)) == Monomial((1, 2, 0))
    assert m.lcm((1, 2, 0)) == Monomial((3, 4, 1))

    pytest.raises(TypeError, lambda: m.gcd(object()))
    pytest.raises(TypeError, lambda: m.lcm(object()))

    assert m**0 == Monomial((0, 0, 0))
    assert m**1 == m
    assert m**2 == Monomial((6, 8, 2))
    assert m**3 == Monomial((9, 12, 3))

    pytest.raises(ValueError, lambda: m**-3)

    assert m/Monomial((5, 2, 0)) == (-2, 2, 1)

    assert str(m) == 'x**3*y**4*z**1'
    assert str(l) == '(3, 4, 1)'
