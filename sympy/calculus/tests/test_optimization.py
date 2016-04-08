from sympy.core import oo, E
from sympy.calculus.optimization import minimize, maximize, simplex
from sympy.functions import exp

from sympy.abc import x, y, z


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


def test_maximize():
    # issue sympy/sympy#4173
    assert maximize([x**(1/x), x > 0], x) == (exp(1/E), {x: E})


def test_simplex():
    assert simplex(-2*x - 3*y - 2*z, [2*x + y + z <= 4,
                                      x + 2*y + z <= 7,
                                      z <= 5]) == (11, [0, 3, 1])
    assert simplex(-2*x - 3*y - 4*z, [3*x + 2*y + z <= 10,
                                      2*x + 5*y + 3*z <= 15]) == (20, [0, 0, 5])
