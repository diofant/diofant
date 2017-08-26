"""Tests for tools and arithmetics for monomials of distributed polynomials. """

import pytest

from diofant.abc import a, b, c, x, y, z
from diofant.core import Integer
from diofant.domains import QQ, ZZ
from diofant.polys.monomials import (Monomial, itermonomials, monomial_count,
                                     monomial_div, monomial_divides,
                                     monomial_gcd, monomial_lcm, monomial_max,
                                     monomial_min, monomial_mul, monomial_pow,
                                     term_div)
from diofant.polys.polyerrors import ExactQuotientFailed


__all__ = ()


def test_monomials():
    assert itermonomials([], 0) == {Integer(1)}
    assert itermonomials([], 1) == {Integer(1)}
    assert itermonomials([], 2) == {Integer(1)}
    assert itermonomials([], 3) == {Integer(1)}

    assert itermonomials([x], 0) == {Integer(1)}
    assert itermonomials([x], 1) == {Integer(1), x}
    assert itermonomials([x], 2) == {Integer(1), x, x**2}
    assert itermonomials([x], 3) == {Integer(1), x, x**2, x**3}

    assert itermonomials([x, y], 0) == {Integer(1)}
    assert itermonomials([x, y], 1) == {Integer(1), x, y}
    assert itermonomials([x, y], 2) == {Integer(1), x, y, x**2, y**2, x*y}
    assert itermonomials([x, y], 3) == \
        {Integer(1), x, y, x**2, x**3, y**2, y**3, x*y, x*y**2, y*x**2}


def test_monomial_count():
    assert monomial_count(2, 2) == 6
    assert monomial_count(2, 3) == 10


def test_monomial_mul():
    assert monomial_mul((3, 4, 1), (1, 2, 0)) == (4, 6, 1)


def test_monomial_div():
    assert monomial_div((3, 4, 1), (1, 2, 0)) == (2, 2, 1)


def test_monomial_pow():
    assert monomial_pow((3, 4, 1), 2) == (6, 8, 2)


def test_monomial_gcd():
    assert monomial_gcd((3, 4, 1), (1, 2, 0)) == (1, 2, 0)


def test_monomial_lcm():
    assert monomial_lcm((3, 4, 1), (1, 2, 0)) == (3, 4, 1)


def test_monomial_max():
    assert monomial_max((3, 4, 5), (0, 5, 1), (6, 3, 9)) == (6, 5, 9)


def test_monomial_min():
    assert monomial_min((3, 4, 5), (0, 5, 1), (6, 3, 9)) == (0, 3, 1)


def test_monomial_divides():
    assert monomial_divides((1, 2, 3), (4, 5, 6)) is True
    assert monomial_divides((1, 2, 3), (0, 5, 6)) is False


def test_term_div():
    assert term_div(((3, 4, 1), 1), ((1, 2, 0), 1), QQ) == ((2, 2, 1), 1)
    assert term_div(((3, 4, 1), 1), ((1, 2, 0), 1), ZZ) == ((2, 2, 1), 1)
    assert term_div(((3, 4, 1), 1), ((1, 2, 2), 1), ZZ) is None
    assert term_div(((3, 4, 1), 1), ((1, 2, 2), 1), QQ) is None


def test_Monomial():
    m = Monomial((3, 4, 1), (x, y, z))
    m2 = Monomial((3, 4, 1))
    n = Monomial((1, 2, 0), (x, y, z))

    assert m.as_expr() == x**3*y**4*z
    assert n.as_expr() == x**1*y**2

    assert m.as_expr(a, b, c) == a**3*b**4*c
    assert n.as_expr(a, b, c) == a**1*b**2

    pytest.raises(ValueError, lambda: m2.as_expr())

    assert m.exponents == (3, 4, 1)
    assert m.gens == (x, y, z)

    assert n.exponents == (1, 2, 0)
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

    pytest.raises(ExactQuotientFailed, lambda: m/Monomial((5, 2, 0)))

    assert str(m) == "x**3*y**4*z**1"
    assert str(m2) == "Monomial((3, 4, 1))"
