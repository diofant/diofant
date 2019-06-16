from diofant import (Function, I, Integer, Rational, cbrt, cos, exp, root,
                     sqrt, tan)
from diofant.abc import a, b, x, y


__all__ = ()


def test_add_eval():
    c = Integer(1)
    p = Integer(5)
    assert a*b + c + p == a*b + 6
    assert c + a + p == a + 6
    assert c + a - p == a + (-4)
    assert a + a == 2*a
    assert a + p + a == 2*a + 5
    assert c + p == 6
    assert b + a - b == a


def test_addmul_eval():
    c = Integer(1)
    p = Integer(5)
    assert c + a + b*c + a - p == 2*a + b + (-4)
    assert a*2 + p + a == a*2 + 5 + a
    assert a*2 + p + a == 3*a + 5
    assert a*2 + a == 3*a


def test_pow_eval():
    # XXX Pow does not fully support conversion of negative numbers
    #     to their complex equivalent

    assert sqrt(-1) == I

    assert sqrt(-4) == 2*I
    assert sqrt(+4) == 2
    assert cbrt(+8) == 2
    assert cbrt(-8) == 2*cbrt(-1)

    assert sqrt(-2) == I*sqrt(2)
    assert cbrt(-1) != I
    assert cbrt(-10) != I*cbrt(10)
    assert root(-2, 4) != root(2, 4)

    assert cbrt(64) == 4
    assert 64**Rational(2, 3) == 16
    assert 24/sqrt(64) == 3
    assert cbrt(-27) == 3*cbrt(-1)

    assert (cos(2) / tan(2))**2 == (cos(2) / tan(2))**2


def test_mulpow_eval():
    assert sqrt(50)/(sqrt(2)*x) == 5/x
    assert sqrt(27)/sqrt(3) == 3


def test_evalpow_bug():
    assert 1/(1/x) == x
    assert 1/(-1/x) == -x


def test_symbol_expand():
    f = x**4*y**4
    assert f == x**4*y**4
    assert f == f.expand()

    g = (x*y)**4
    assert g == f
    assert g.expand() == f
    assert g.expand() == g.expand().expand()


def test_function():
    f, l = map(Function, 'fl')
    assert exp(l(x))*l(x)/exp(l(x)) == l(x)
    assert exp(f(x))*f(x)/exp(f(x)) == f(x)
