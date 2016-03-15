from sympy.core import oo, E
from sympy.calculus import minimize, maximize
from sympy.functions import exp

from sympy.abc import x


def test_minimize():
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
