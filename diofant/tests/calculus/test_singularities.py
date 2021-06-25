from diofant import exp, log, sign
from diofant.abc import x
from diofant.calculus.singularities import singularities


__all__ = ()


def test_singularities():
    assert singularities(1, x) == set()
    assert singularities(x**2, x) == set()
    assert singularities(x/(x**2 + 3*x + 2), x) == {-2, -1}
    assert singularities(1/(1 + x), x) == {-1}
    assert singularities(1/(1 + x)**2, x) == {-1}
    assert singularities(sign(x), x) == {0}


def test_singularities_non_rational():
    assert singularities(exp(1/x), x) == {0}
    assert singularities(log((x - 2)**2), x) == {2}
    assert singularities(exp(1/x) + log(x + 1), x) == {-1, 0}
    assert singularities(exp(1/log(x + 1)), x) == {0}
