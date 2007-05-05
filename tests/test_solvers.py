import sys
sys.path.append(".")

from sympy import Rational, Symbol, cos, solve, dsolve, Function, diff, \
        log, sin, exp
from sympy.core.functions import Derivative

import decimal

def test_linear():
    x = Symbol("x")
    assert solve(3*x-2, x) == Rational(2,3)

def test_quadratic():
    x = Symbol("x")
    assert solve(x**2-1, x) == [1, -1]
    assert solve(((x-1)*(x-2)).expand(), x) in [[1,2], [2,1]]
    assert solve(((x-1)*(x-1)).expand(), x) == [1]
    
def test_cubic():
    x = Symbol('x')
    f = x**3 - x**2 + x +1
  #  for root in solve(f, x):
  #      assert f.subs(x, root).evalf() < decimal.Decimal("1e-10")

def test_ODE_first_order():
    x = Symbol("x")
    f = Function(x)
    assert dsolve(3*diff(f, x) -1, f) == x/3+Symbol("C1")
    assert dsolve(x*diff(f, x) -1, f) == log(abs(x))+Symbol("C1")

def test_ODE_second_order():
    x = Symbol("x")
    f = Function(x)
    C1, C2 = Symbol("C1"), Symbol("C2")
    assert dsolve(Derivative(Derivative(f,x),x)+9*f, [f]) == \
        sin(3*x)*C1+cos(3*x)*C2

def test_ODE_1():
    class l(Function): pass

    r = Symbol("r")

    e = Derivative(l(r),r)/r+Derivative(Derivative(l(r),r),r)/2- \
        Derivative(l(r),r)**2/2

    sol = dsolve(e, [l(r)])
    assert (e.subs(l(r), sol).doit()).expand() == 0

    e = e*exp(-l(r))/exp(l(r))

    sol = dsolve(e, [l(r)])
    assert (e.subs(l(r), sol).doit()).expand() == 0
