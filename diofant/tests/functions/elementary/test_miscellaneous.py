import pytest

from diofant.abc import x, y, z
from diofant.core import (Eq, Function, I, Pow, Rational, Symbol, oo, symbols,
                          zoo)
from diofant.core.function import ArgumentIndexError
from diofant.functions import (Heaviside, Max, Min, Piecewise, cbrt, ceiling,
                               cos, floor, real_root, root, sin, sqrt)
from diofant.logic import true


__all__ = ()


def test_Min():
    n = Symbol('n', negative=True)
    n_ = Symbol('n_', negative=True)
    nn = Symbol('nn', nonnegative=True)
    nn_ = Symbol('nn_', nonnegative=True)
    p = Symbol('p', positive=True)
    p_ = Symbol('p_', positive=True)
    np = Symbol('np', nonpositive=True)
    np_ = Symbol('np_', nonpositive=True)

    assert Min(5, 4) == 4
    assert Min(-oo, -oo) == -oo
    assert Min(-oo, n) == -oo
    assert Min(n, -oo) == -oo
    assert Min(-oo, np) == -oo
    assert Min(np, -oo) == -oo
    assert Min(-oo, 0) == -oo
    assert Min(0, -oo) == -oo
    assert Min(-oo, nn) == -oo
    assert Min(nn, -oo) == -oo
    assert Min(-oo, p) == -oo
    assert Min(p, -oo) == -oo
    assert Min(-oo, oo) == -oo
    assert Min(oo, -oo) == -oo
    assert Min(n, n) == n
    assert Min(n, np) == Min(n, np)
    assert Min(np, n) == Min(np, n)
    assert Min(n, 0) == n
    assert Min(0, n) == n
    assert Min(n, nn) == n
    assert Min(nn, n) == n
    assert Min(n, p) == n
    assert Min(p, n) == n
    assert Min(n, oo) == n
    assert Min(oo, n) == n
    assert Min(np, np) == np
    assert Min(np, 0) == np
    assert Min(0, np) == np
    assert Min(np, nn) == np
    assert Min(nn, np) == np
    assert Min(np, p) == np
    assert Min(p, np) == np
    assert Min(np, oo) == np
    assert Min(oo, np) == np
    assert Min(0, 0) == 0
    assert Min(0, nn) == 0
    assert Min(nn, 0) == 0
    assert Min(0, p) == 0
    assert Min(p, 0) == 0
    assert Min(0, oo) == 0
    assert Min(oo, 0) == 0
    assert Min(nn, nn) == nn
    assert Min(nn, p) == Min(nn, p)
    assert Min(p, nn) == Min(p, nn)
    assert Min(nn, oo) == nn
    assert Min(oo, nn) == nn
    assert Min(p, p) == p
    assert Min(p, oo) == p
    assert Min(oo, p) == p
    assert Min(oo, oo) == oo

    assert isinstance(Min(n, n_), Min)
    assert isinstance(Min(nn, nn_), Min)
    assert isinstance(Min(np, np_), Min)
    assert isinstance(Min(p, p_), Min)

    # lists
    pytest.raises(ValueError, lambda: Min())
    assert Min(x, y) == Min(y, x)
    assert Min(x, y, z) == Min(z, y, x)
    assert Min(x, Min(y, z)) == Min(z, y, x)
    assert Min(x, Max(y, -oo)) == Min(x, y)
    assert Min(p, oo, n, p, p, p_) == n
    assert Min(p_, n_, p) == n_
    assert Min(n, oo, -7, p, p, 2) == Min(n, -7)
    assert Min(2, x, p, n, oo, n_, p, 2, -2, -2) == Min(-2, x, n, n_)
    assert Min(0, x, 1, y) == Min(0, x, y)
    assert Min(1000, 100, -100, x, p, n) == Min(n, x, -100)
    assert Min(cos(x), sin(x)) == Min(cos(x), sin(x))
    assert Min(cos(x), sin(x)).subs({x: 1}) == cos(1)
    assert Min(cos(x), sin(x)).subs({x: Rational(1, 2)}) == sin(Rational(1, 2))
    pytest.raises(ValueError, lambda: Min(cos(x), sin(x)).subs({x: I}))
    pytest.raises(ValueError, lambda: Min(I))
    pytest.raises(ValueError, lambda: Min(I, x))
    pytest.raises(ValueError, lambda: Min(zoo, x))

    assert Min(1, x).diff(x) == Heaviside(1 - x)
    assert Min(x, 1).diff(x) == Heaviside(1 - x)
    assert Min(0, -x, 1 - 2*x).diff(x) == -Heaviside(x + Min(0, -2*x + 1)) \
        - 2*Heaviside(2*x + Min(0, -x) - 1)

    pytest.raises(ArgumentIndexError, lambda: Min(1, x).fdiff(3))

    a, b = Symbol('a', extended_real=True), Symbol('b', extended_real=True)
    # a and b are both real, Min(a, b) should be real
    assert Min(a, b).is_extended_real

    # issue sympy/sympy#7619
    f = Function('f')
    assert Min(1, 2*Min(f(1), 2))  # doesn't fail

    # issue sympy/sympy#7233
    e = Min(0, x)
    assert e.evalf == e.n
    assert e.evalf(strict=False).args == (0, x)


def test_Max():
    n = Symbol('n', negative=True)
    n_ = Symbol('n_', negative=True)
    p = Symbol('p', positive=True)
    r = Symbol('r', extended_real=True)

    assert Max(5, 4) == 5

    # lists

    pytest.raises(ValueError, lambda: Max())
    assert Max(x, y) == Max(y, x)
    assert Max(x, y, z) == Max(z, y, x)
    assert Max(x, Max(y, z)) == Max(z, y, x)
    assert Max(x, Min(y, oo)) == Max(x, y)
    assert Max(n, -oo, n_, p, 2) == Max(p, 2)
    assert Max(n, -oo, n_, p) == p
    assert Max(2, x, p, n, -oo, -oo, n_, p, 2) == Max(2, x, p)
    assert Max(0, x, 1, y) == Max(1, x, y)
    assert Max(r, r + 1, r - 1) == 1 + r
    assert Max(1000, 100, -100, x, p, n) == Max(p, x, 1000)
    assert Max(cos(x), sin(x)) == Max(sin(x), cos(x))
    assert Max(cos(x), sin(x)).subs({x: 1}) == sin(1)
    assert Max(cos(x), sin(x)).subs({x: Rational(1, 2)}) == cos(Rational(1, 2))
    pytest.raises(ValueError, lambda: Max(cos(x), sin(x)).subs({x: I}))
    pytest.raises(ValueError, lambda: Max(I))
    pytest.raises(ValueError, lambda: Max(I, x))
    pytest.raises(ValueError, lambda: Max(zoo, 1))
    # interesting:
    # Max(n, -oo, n_,  p, 2) == Max(p, 2)
    # True
    # Max(n, -oo, n_,  p, 1000) == Max(p, 1000)
    # False

    assert Max(1, x).diff(x) == Heaviside(x - 1)
    assert Max(x, 1).diff(x) == Heaviside(x - 1)
    assert Max(x**2, 1 + x, 1).diff(x) == \
        2*x*Heaviside(x**2 - Max(1, x + 1)) \
        + Heaviside(x - Max(1, x**2) + 1)

    pytest.raises(ArgumentIndexError, lambda: Max(1, x).fdiff(3))

    a, b = Symbol('a', extended_real=True), Symbol('b', extended_real=True)
    # a and b are both real, Max(a, b) should be real
    assert Max(a, b).is_extended_real

    # issue sympy/sympy#7233
    e = Max(0, x)
    assert e.evalf == e.n
    assert e.evalf(strict=False).args == (0, x)


def test_sympyissue_8413():
    x = Symbol('x', extended_real=True)
    # we can't evaluate in general because non-reals are not
    # comparable: Min(floor(3.2 + I), 3.2 + I) -> ValueError
    assert Min(floor(x), x) == floor(x)
    assert Min(ceiling(x), x) == x
    assert Max(floor(x), x) == x
    assert Max(ceiling(x), x) == ceiling(x)


def test_root():
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True)

    assert root(2, 2) == sqrt(2)
    assert root(2, 1) == 2
    assert root(2, 3) == cbrt(2)
    assert root(2, 3) == cbrt(2)
    assert root(2, -5) == 2**Rational(4, 5)/2

    assert root(4, 2, evaluate=False) == Pow(4, Rational(1, 2), evaluate=False)
    assert root(4, 2, 1, evaluate=False) == -Pow(4, Rational(1, 2), evaluate=False)

    assert root(-2, 1) == -2

    assert root(-2, 2) == sqrt(2)*I
    assert root(-2, 1) == -2

    assert root(x, 2) == sqrt(x)
    assert root(x, 1) == x
    assert root(x, 3) == cbrt(x)
    assert root(x, 3) == cbrt(x)
    assert root(x, -5) == x**Rational(-1, 5)

    assert root(x, n) == x**(1/n)
    assert root(x, -n) == x**(-1/n)

    assert root(x, n, k) == x**(1/n)*(-1)**(2*k/n)


def test_real_root():
    assert real_root(-8, 3) == -2
    assert real_root(-16, 4) == root(-16, 4)
    r = root(-7, 4)
    assert real_root(r) == r
    r1 = root(-1, 3)
    r2 = r1**2
    r3 = root(-1, 4)
    assert real_root(r1 + r2 + r3) == -1 + r2 + r3
    assert real_root(root(-2, 3)) == -root(2, 3)
    assert real_root(-8., 3) == -2
    x = Symbol('x')
    n = Symbol('n')
    g = real_root(x, n)
    assert g.subs({x: -8, n: 3}) == -2
    assert g.subs({x: 8, n: 3}) == 2
    # give principle root if there is no real root -- if this is not desired
    # then maybe a Root class is needed to raise an error instead
    assert g.subs({x: I, n: 3}) == cbrt(I)
    assert g.subs({x: -8, n: 2}) == sqrt(-8)
    assert g.subs({x: I, n: 2}) == sqrt(I)


def test_rewrite_MaxMin_as_Heaviside():
    assert Max(0, x).rewrite(Heaviside) == x*Heaviside(x)
    assert Max(3, x).rewrite(Heaviside) == x*Heaviside(x - 3) + \
        3*Heaviside(-x + 3)
    assert Max(0, x+2, 2*x).rewrite(Heaviside) == \
        2*x*Heaviside(2*x)*Heaviside(x - 2) + \
        (x + 2)*Heaviside(-x + 2)*Heaviside(x + 2)

    assert Min(0, x).rewrite(Heaviside) == x*Heaviside(-x)
    assert Min(3, x).rewrite(Heaviside) == x*Heaviside(-x + 3) + \
        3*Heaviside(x - 3)
    assert Min(x, -x, -2).rewrite(Heaviside) == \
        x*Heaviside(-2*x)*Heaviside(-x - 2) - \
        x*Heaviside(2*x)*Heaviside(x - 2) \
        - 2*Heaviside(-x + 2)*Heaviside(x + 2)


def test_rewrite_as_Piecewise():
    x, y = symbols('x, y', real=True)
    assert (Max(x, y).rewrite(Piecewise) ==
            x*Piecewise((1, x - y > 0), (Rational(1, 2), Eq(x - y, 0)), (0, true)) +
            y*Piecewise((1, -x + y > 0), (Rational(1, 2), Eq(-x + y, 0)), (0, true)))


def test_sympyissue_11099():
    e = Min(x, y)
    s = {x: -2, y: 3}
    assert e.evalf(subs=s) == e.subs(s).evalf()
