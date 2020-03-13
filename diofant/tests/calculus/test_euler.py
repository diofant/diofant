"""Tests for Euler-Lagrange equations."""

import pytest

from diofant import Eq, Function, Symbol, cos, diff, sin
from diofant.calculus.euler import euler_equations as euler


__all__ = ()


def test_euler_interface():
    x = Function('x')
    y = Symbol('y')
    t = Symbol('t')
    pytest.raises(TypeError, lambda: euler())
    pytest.raises(TypeError, lambda: euler(diff(x(t), t)*y(t), [x(t), y]))
    pytest.raises(ValueError, lambda: euler(diff(x(t), t)*x(y), [x(t), x(y)]))
    pytest.raises(TypeError, lambda: euler(diff(x(t), t)**2, x(0)))
    pytest.raises(TypeError, lambda: euler(1, y))
    assert euler(diff(x(t), t)**2/2, {x(t)}) == [Eq(-diff(x(t), t, t), 0)]
    assert euler(diff(x(t), t)**2/2, x(t), {t}) == [Eq(-diff(x(t), t, t), 0)]


def test_euler_pendulum():
    x = Function('x')
    t = Symbol('t')
    L = diff(x(t), t)**2/2 + cos(x(t))
    assert euler(L, x(t), t) == [Eq(-sin(x(t)) - diff(x(t), t, t), 0)]


def test_euler_henonheiles():
    x = Function('x')
    y = Function('y')
    t = Symbol('t')
    L = sum(diff(z(t), t)**2/2 - z(t)**2/2 for z in [x, y])
    L += -x(t)**2*y(t) + y(t)**3/3
    assert euler(L, [x(t), y(t)], t) == [Eq(-2*x(t)*y(t) - x(t) -
                                            diff(x(t), t, t), 0),
                                         Eq(-x(t)**2 + y(t)**2 -
                                            y(t) - diff(y(t), t, t), 0)]


def test_euler_sineg():
    psi = Function('psi')
    t = Symbol('t')
    x = Symbol('x')
    L = diff(psi(t, x), t)**2/2 - diff(psi(t, x), x)**2/2 + cos(psi(t, x))
    assert euler(L, psi(t, x), [t, x]) == [Eq(-sin(psi(t, x)) -
                                              diff(psi(t, x), t, t) +
                                              diff(psi(t, x), x, x), 0)]


def test_euler_high_order():
    # an example from hep-th/0309038
    m = Symbol('m')
    k = Symbol('k')
    x = Function('x')
    y = Function('y')
    t = Symbol('t')
    L = (m*diff(x(t), t)**2/2 + m*diff(y(t), t)**2/2 -
         k*diff(x(t), t)*diff(y(t), t, t) + k*diff(y(t), t)*diff(x(t), t, t))
    assert euler(L, [x(t), y(t)]) == [Eq(2*k*diff(y(t), t, t, t) -
                                         m*diff(x(t), t, t), 0),
                                      Eq(-2*k*diff(x(t), t, t, t) -
                                         m*diff(y(t), t, t), 0)]

    w = Symbol('w')
    L = diff(x(t, w), t, w)**2/2
    assert euler(L) == [Eq(diff(x(t, w), t, t, w, w), 0)]
