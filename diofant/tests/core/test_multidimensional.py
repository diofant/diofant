import pytest

from diofant import (Derivative, E, Function, Integer, Rational, diff, log, pi,
                     sin, vectorize, zeta, zoo)
from diofant.abc import x, y, z


__all__ = ()


def test_vectorize():
    vsin1 = vectorize(Integer(0))(sin)
    pytest.raises(TypeError, lambda: vectorize(Rational(1, 2))(sin))

    @vectorize()
    def vsin2(x):
        return sin(x)

    @vectorize(0, 'evaluate')
    def vsin3(x, evaluate=True):
        return sin(x, evaluate=evaluate)

    @vectorize(0, 'evaluate')
    def vzeta(z, a=None, evaluate=True):
        return zeta(z, a, evaluate=evaluate)

    @vectorize(0, 3, 'boo')
    def vlog(x, base=E, **kwargs):
        return log(x)/log(base)

    @vectorize(0, 1)
    def vdiff(f, y):
        return diff(f, y)

    f, g = map(Function, 'fg')

    assert vsin1([1, x, y]) == [sin(1), sin(x), sin(y)]
    assert vsin1((1, (x, y))) == [sin(1), [sin(x), sin(y)]]
    assert vsin1((1, x, y)) == [sin(1), sin(x), sin(y)]
    assert vsin2([1, x, y]) == [sin(1), sin(x), sin(y)]
    assert (vsin3([pi, 2*pi, pi/2], evaluate=[True, False]) ==
            [[0, sin(pi, evaluate=False)],
             [0, sin(2*pi, evaluate=False)],
             [1, sin(pi/2, evaluate=False)]])

    assert (vzeta([1, 2], evaluate=[True, False], a=2) ==
            [[zoo, zeta(1, 2, evaluate=False)],
             [-1 + pi**2/6, zeta(2, 2, evaluate=False)]])
    assert vlog([1, 2]) == [0, log(2)]

    assert (vdiff([f(x, y, z), g(x, y, z)], [x, y, z]) ==
            [[Derivative(f(x, y, z), x), Derivative(f(x, y, z), y),
              Derivative(f(x, y, z), z)],
             [Derivative(g(x, y, z), x), Derivative(g(x, y, z), y),
              Derivative(g(x, y, z), z)]])
