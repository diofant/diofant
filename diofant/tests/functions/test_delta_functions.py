import pytest

from diofant import (DiracDelta, Eq, Heaviside, I, Piecewise, Rational, Symbol,
                     adjoint, conjugate, nan, pi, sign, sqrt, symbols,
                     transpose)
from diofant.abc import x, y, z
from diofant.core.function import ArgumentIndexError


__all__ = ()


def test_DiracDelta():
    i = Symbol('i', nonzero=True)
    j = Symbol('j', positive=True)
    k = Symbol('k', negative=True)

    assert DiracDelta(1) == 0
    assert DiracDelta(5.1) == 0
    assert DiracDelta(-pi) == 0
    assert DiracDelta(5, 7) == 0
    assert DiracDelta(i) == 0
    assert DiracDelta(j) == 0
    assert DiracDelta(k) == 0
    assert DiracDelta(nan) == nan
    assert isinstance(DiracDelta(0), DiracDelta)
    assert isinstance(DiracDelta(x), DiracDelta)

    assert adjoint(DiracDelta(x)) == DiracDelta(x)
    assert adjoint(DiracDelta(x - y)) == DiracDelta(x - y)
    assert conjugate(DiracDelta(x)) == DiracDelta(x)
    assert conjugate(DiracDelta(x - y)) == DiracDelta(x - y)
    assert transpose(DiracDelta(x)) == DiracDelta(x)
    assert transpose(DiracDelta(x - y)) == DiracDelta(x - y)

    assert DiracDelta(x).diff(x) == DiracDelta(x, 1)
    assert DiracDelta(x, 1).diff(x) == DiracDelta(x, 2)

    assert DiracDelta(x).is_simple(x) is True
    assert DiracDelta(3*x).is_simple(x) is True
    assert DiracDelta(x**2).is_simple(x) is False
    assert DiracDelta(sqrt(x)).is_simple(x) is False
    assert DiracDelta(x).is_simple(y) is False

    assert DiracDelta(x*y).simplify(x) == DiracDelta(x)/abs(y)
    assert DiracDelta(x*y).simplify(y) == DiracDelta(y)/abs(x)
    assert DiracDelta(x**2*y).simplify(x) == DiracDelta(x**2*y)
    assert DiracDelta(y).simplify(x) == DiracDelta(y)
    assert DiracDelta((x - 1)*(x - 2)*(x - 3)).simplify(x) == \
        DiracDelta(x - 3)/2 + DiracDelta(x - 2) + DiracDelta(x - 1)/2

    pytest.raises(ArgumentIndexError, lambda: DiracDelta(x).fdiff(2))
    pytest.raises(ValueError, lambda: DiracDelta(x, -1))


def test_heaviside():
    x, y = symbols('x, y', extended_real=True)
    assert Heaviside(0) == 0.5
    assert Heaviside(-5) == 0
    assert Heaviside(1) == 1
    assert Heaviside(nan) == nan

    assert Heaviside(x).is_real

    assert adjoint(Heaviside(x)) == Heaviside(x)
    assert adjoint(Heaviside(x - y)) == Heaviside(x - y)
    assert conjugate(Heaviside(x)) == Heaviside(x)
    assert conjugate(Heaviside(x - y)) == Heaviside(x - y)
    assert transpose(Heaviside(x)) == Heaviside(x)
    assert transpose(Heaviside(x - y)) == Heaviside(x - y)

    assert Heaviside(x).diff(x) == DiracDelta(x)
    assert Heaviside(z + I).is_Function is True
    assert Heaviside(I*z).is_Function is True

    pytest.raises(ArgumentIndexError, lambda: Heaviside(x).fdiff(2))
    pytest.raises(ValueError, lambda: Heaviside(I))
    pytest.raises(ValueError, lambda: Heaviside(2 + 3*I))


def test_rewrite():
    x = Symbol('x', extended_real=True)
    assert Heaviside(x).rewrite(Piecewise) == \
        Piecewise((1, x > 0), (Rational(1, 2), Eq(x, 0)), (0, True))
    assert Heaviside(y).rewrite(Piecewise) == Heaviside(y)

    assert Heaviside(x).rewrite(sign) == (sign(x)+1)/2
    assert Heaviside(y).rewrite(sign) == Heaviside(y)
