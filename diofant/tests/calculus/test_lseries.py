from diofant import (Derivative, E, I, Integer, O, cos, exp, gamma, log, pi,
                     sin, tanh)
from diofant.abc import x, y


__all__ = ()


def test_sin():
    e = sin(x).series(x, n=None)
    assert next(e) == x
    assert next(e) == -x**3/6
    assert next(e) == x**5/120


def test_cos():
    e = cos(x).series(x, n=None)
    assert next(e) == 1
    assert next(e) == -x**2/2
    assert next(e) == x**4/24


def test_exp():
    e = exp(x).series(x, n=None)
    assert next(e) == 1
    assert next(e) == x
    assert next(e) == x**2/2
    assert next(e) == x**3/6


def test_exp2():
    e = exp(cos(x)).series(x, n=None)
    assert next(e) == E
    assert next(e) == -E*x**2/2
    assert next(e) == E*x**4/6
    assert next(e) == -31*E*x**6/720


def test_simple():
    assert list(x.series(n=None)) == [x]
    assert list(Integer(1).series(x, n=None)) == [1]
    assert not next((x/(x + y)).series(y, n=None)).has(O)
    assert [t.doit() for t in Derivative(1 + x, x).series(x, n=None)] == [0, 1]

    # issue sympy/sympy#5183
    s = (x + 1/x).series(n=None)
    assert list(s) == [1/x, x]

    # issue sympy/sympy#5183
    s = (x + 1/x).series(n=None)
    assert next((x + x**2).series(n=None)) == x
    assert next(((1 + x)**7).series(x, n=None)) == 1
    assert next((sin(x + y)).series(x, n=3).series(y, n=None)) == x
    # it would be nice if all terms were grouped, but in the
    # following case that would mean that all the terms would have
    # to be known since, for example, every term has a constant in it.
    s = ((1 + x)**7).series(x, 1, n=None)
    assert [next(s) for i in range(2)] == [128, -448 + 448*x]


def test_tanh():
    # issue sympy/sympy#6999
    s = tanh(x).series(x, x0=1, n=None)
    assert next(s) == tanh(1)
    assert next(s) == x - (x - 1)*tanh(1)**2 - 1
    assert next(s) == -(x - 1)**2*tanh(1) + (x - 1)**2*tanh(1)**3
    assert next(s) == -(x - 1)**3*tanh(1)**4 - (x - 1)**3/3 + \
        4*(x - 1)**3*tanh(1)**2/3


def test_sympyissue_21859():
    e = -exp(-x*log(y)/2)*exp(-5*I*pi*x/6)*gamma(5*x/6 + 1)/(gamma(x/3 + 1)*gamma(x/2 + 1))
    s = e.series(x, n=None)

    assert next(s) == -1
    assert next(s) == x*log(y)/2 + 5*I*pi*x/6
