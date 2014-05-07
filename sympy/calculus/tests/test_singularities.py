from sympy import Symbol, exp, log, sign
from sympy.calculus.singularities import singularities

from sympy.abc import x


def test_singularities():
    assert singularities(x**2, x) == set()
    assert singularities(x/(x**2 + 3*x + 2), x) == {-2, -1}
    assert singularities(1/(1 + x), x) == {-1}
    assert singularities(sign(x), x) == {0}


def test_singularities_non_rational():
    assert singularities(exp(1/x), x) == {0}
    assert singularities(log((x - 2)**2), x) == {2}
    assert singularities(exp(1/x) + log(x + 1), x) == {-1, 0}
    assert singularities(exp(1/log(x + 1)), x) == {0}
