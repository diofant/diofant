"""Symbolic optimization tests."""

import pytest

from diofant import (E, Eq, Rational, exp, maximize, minimize, nan, oo, sign,
                     sqrt)
from diofant.abc import t, w, x, y, z
from diofant.calculus.optimization import InfeasibleProblemError, simplex


__all__ = ()


def test_minimize():
    assert minimize(1) == (1, {})
    assert minimize((x - 2)**2) == (0, {x: 2})
    assert minimize((x - 2)**2, x) == (0, {x: 2})
    assert minimize(1/x, x) == (-oo, {x: 0})
    assert minimize(2*x**2 - 4*x + 5, x) == (3, {x: 1})
    assert minimize([2*x**2 - 4*x + 5, x > 0], x) == (3, {x: 1})
    assert minimize([2*x**2 - 4*x + 5, x > 1], x) is None
    assert minimize([2*x**2 - 4*x + 5, x >= 2], x) == (5, {x: 2})
    assert minimize([2*x**2 - 4*x + 5, x**2 <= 0], x) == (5, {x: 0})
    assert minimize([x**2 - 1,
                     (x - 1)*(x - 2)*(x - 3)*(4 - x) >= 0]) == (0, {x: 1})
    assert minimize([x**2, (x + 2)*(x - 1)*(1 - x) >= 0]) == (1, {x: 1})
    assert minimize(sign(x), x) == (-1, {x: 0})


def test_maximize():
    # issue sympy/sympy#4173
    assert maximize([x**(1/x), x > 0], x) == (exp(1/E), {x: E})

    # https://groups.google.com/forum/#!topic/sympy/tB2Sly4Gh_4
    assert maximize([12*x + 40*y, x + y <= 16, x + 3*y <= 36, x <= 10,
                     x >= 0, y >= 0], x, y) == (480, {x: 0, y: 12})


def test_minimize_linear():
    # see also sympy/sympy#20391
    assert minimize([-2*x - 3*y - 2*z, 2*x + y + z <= 4,
                     x + 2*y + z <= 7, z <= 5, x >= 0,
                     y >= 0, z >= 0], x, y, z) == (-11, {x: 0, y: 3, z: 1})
    assert minimize([-2*x - 3*y - 2*z, 2*x + y + z <= 4,
                     x + 2*y + z <= 7, z <= 5, x >= 0,
                     y >= 0, z >= 0], x, y, z) == (-11, {x: 0, y: 3, z: 1})
    assert minimize([-2*x - 3*y - 2*z, 2*x + y + z <= 4,
                     x + 2*y + z < 7, z <= 5, x >= 0,
                     y >= 0, z >= 0], x, y, z) is None
    assert minimize([-2*x - 3*y - 4*z, 3*x + 2*y + z <= 10,
                     2*x + 5*y + 3*z <= 15, x >= 0,
                     y >= 0, z >= 0], x, y, z) == (-20, {x: 0, y: 0, z: 5})
    assert maximize([12*x + 40*y, x + y <= 15, x + 3*y <= 36,
                     x <= 10, x >= 0, y >= 0], x, y) == (480, {x: 0, y: 12})
    assert minimize([-2*x - 3*y - 4*z, Eq(3*x + 2*y + z, 10),
                     Eq(2*x + 5*y + 3*z, 15), x >= 0, y >= 0, z >= 0],
                    x, y, z) == (Rational(-130, 7),
                                 {x: Rational(15, 7), y: 0,
                                  z: Rational(25, 7)})
    assert maximize([2*x + y, y - x <= 1, x - 2*y <= 2,
                     x >= 0, y >= 0], x, y) == (oo, {x: nan, y: nan})

    assert minimize([2*x + 3*y - z, 1 <= x + y + z,
                     x + y + z <= 2, 1 <= x - y + z,
                     x - y + z <= 2, Eq(x - y - z, 3)],
                    x, y, z) == (3, {x: 2, y: Rational(-1, 2),
                                     z: Rational(-1, 2)})

    assert minimize([-2*t + 4*w + 7*x + y + 5*z,
                     Eq(-t + w + 2*x + y + 2*z, 7),
                     Eq(-t + 2*w + 3*x + y + z, 6),
                     Eq(-t + w + x + 2*y + z, 4),
                     w >= 0, x >= 0, y >= 0, z >= 0],
                    t, w, x, y, z) == (19, {t: -1, w: 0,
                                            x: 1, y: 0, z: 2})
    assert minimize([-x - y, x + 2*y <= 8, 3*x + 2*y <= 12,
                     x + 3*y >= 6, x >= 0, y >= 0],
                    x, y) == (-5, {x: 2, y: 3})
    pytest.raises(InfeasibleProblemError,
                  lambda: minimize([-x - y, x + 2*y <= 8, 3*x + 2*y <= 12,
                                    x + 3*y >= 13, x >= 0, y >= 0],
                                   x, y))
    pytest.raises(InfeasibleProblemError,
                  lambda: minimize([-x - y, x + 2*y <= 8, 3*x + 2*y <= 12,
                                   x + 3*y >= 13, x >= 0, y >= 0], x, y))
    assert minimize([-x - y, 2*x + y >= 4,
                     Eq(x + 2*y, 6), x >= 0, y >= 0],
                    x, y) == (-6, {x: 6, y: 0})

    assert minimize([6*x + 3*y, x + y >= 1, 2*x - y >= 1,
                     3*y <= 2, x >= 0,
                     y >= 0], x, y) == (5, {x: Rational(2, 3),
                                            y: Rational(1, 3)})


@pytest.mark.xfail
def test_minimize_poly():
    assert minimize([x - 2*y, x**2 + y**2 <= 1], x, y)[0] == -sqrt(5)
    # at {x: 4/sqrt(5) - sqrt(5), y: 2/sqrt(5)}


@pytest.mark.xfail
def test_minimize_analytic():
    assert minimize(exp(x**2 + y**2), x, y) == (1, {x: 0, y: 0})


def test_simplex():
    pytest.raises(ValueError, lambda: simplex([1], [[1, 2], [2, 3]],
                                              [4, 7]))
    pytest.raises(ValueError, lambda: simplex([1, 3], [[1, 2],
                                                       [2, 3]], [4]))

    pytest.raises(InfeasibleProblemError,
                  lambda: simplex([1, 3], [[1, 2], [2, 3]], [-3, 1]))

    assert simplex([2, 3, 2], [[2, 1, 1], [1, 2, 1], [0, 0, 1]],
                   [4, 7, 5]) == (11, (0, 3, 1))
    assert simplex([2, 3, 4], [[3, 2, 1], [2, 5, 3]],
                   [10, 15]) == (20, (0, 0, 5))
    assert simplex([Rational(1, 2), 3, 1, 1],
                   [[1, 1, 1, 1],
                    [-2, -1, 1, 1],
                    [0, 1, 0, -1]],
                   [40, 10, 10]) == (90, (0, 25, 0, 15))

    assert simplex([2, 1], [[-1, 1], [1, -2]], [1, 2]) == (oo, (oo, oo))

    assert simplex([1, 2, 3, -2], [[3, -2, 1, 1],
                                   [-2, 1, 10, -1],
                                   [2, 0, 0, 1]],
                   [1, -2, 3]) == (-4, (0, 1, 0, 3))

    pytest.raises(InfeasibleProblemError,
                  lambda: simplex([1, 2, 3, -2],
                                  [[3, -2, 1, 1], [-2, 1, 10, -1],
                                   [2, 0, 0, 1], [0, -1, 2, 0]],
                                  [1, -2, 3, -3]))
