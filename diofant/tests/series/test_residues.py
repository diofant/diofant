import pytest

from diofant import (Function, I, Mul, Rational, Symbol, cot, exp, factorial,
                     log, pi, residue, root, sin, sqrt, tan, tanh)
from diofant.abc import a, s, x, z


__all__ = ()


def test_basic1():
    assert residue(1/x, x, 0) == 1
    assert residue(-2/x, x, 0) == -2
    assert residue(81/x, x, 0) == 81
    assert residue(1/x**2, x, 0) == 0
    assert residue(0, x, 0) == 0
    assert residue(5, x, 0) == 0
    assert residue(x, x, 0) == 0
    assert residue(x**2, x, 0) == 0


def test_basic2():
    assert residue(1/x, x, 1) == 0
    assert residue(-2/x, x, 1) == 0
    assert residue(81/x, x, -1) == 0
    assert residue(1/x**2, x, 1) == 0
    assert residue(0, x, 1) == 0
    assert residue(5, x, 1) == 0
    assert residue(x, x, 1) == 0
    assert residue(x**2, x, 5) == 0


def test_f():
    f = Function('f')
    assert residue(f(x)/x**5, x, 0) == f(x).diff((x, 4)).subs({x: 0})/24


def test_functions():
    assert residue(1/sin(x), x, 0) == 1
    assert residue(2/sin(x), x, 0) == 2
    assert residue(1/sin(x)**2, x, 0) == 0
    assert residue(1/sin(x)**5, x, 0) == Rational(3, 8)


def test_expressions():
    assert residue(1/(x + 1), x, 0) == 0
    assert residue(1/(x + 1), x, -1) == 1
    assert residue(1/(x**2 + 1), x, -1) == 0
    assert residue(1/(x**2 + 1), x, I) == -I/2
    assert residue(1/(x**2 + 1), x, -I) == I/2
    assert residue(1/(x**4 + 1), x, 0) == 0
    assert residue(1/(x**4 + 1), x, exp(I*pi/4)) == -root(-1, 4)/4
    assert residue(1/(x**2 + a**2)**2, x, a*I) == -I/4/a**3


@pytest.mark.xfail
def test_expressions_failing():
    n = Symbol('n', integer=True, positive=True)
    assert residue(exp(z)/(z - pi*I/4*a)**n, z, I*pi*a) == \
        exp(I*pi*a/4)/factorial(n - 1)


def test_NotImplemented():
    pytest.raises(NotImplementedError, lambda: residue(exp(1/z), z, 0))


def test_bug():
    assert residue(2**(z)*(s + z)*(1 - s - z)/z**2, z, 0) == \
        1 + s*log(2) - s**2*log(2) - 2*s


def test_sympyissue_5654():
    assert residue(1/(x**2 + a**2)**2, x, a*I) == -I/(4*a**3)


def test_sympyissue_6499():
    assert residue(1/(exp(z) - 1), z, 0) == 1


def test_sympyissue_21177():
    e1 = cot(pi*x)/((x - 1)*(x - 2) + 1)
    e2 = cot(pi*x)/(x**2 - 3*x + 3)
    pt = Rational(3, 2) - sqrt(3)*I/2
    ans = -sqrt(3)*tanh(sqrt(3)*pi/2)/3

    assert residue(e1, x, pt) == ans
    assert residue(e2, x, pt) == ans


def test_sympyissue_21176():
    e = x**2*cot(pi*x)/(x**4 + 1)
    pt = -sqrt(2)/2 - sqrt(2)*I/2
    assert residue(e, x, pt) == sqrt(2)*I/Mul(2, -2 + 2*I, tan(sqrt(2)*pi/2 +
                                                               sqrt(2)*I*pi/2),
                                              evaluate=False)
