import sys
sys.path.append(".")

from sympy import Rational, Symbol, cos, solve, dsolve, Function, Derivative, \
        log, sin

def test_linear():
    x = Symbol("x")
    assert solve(3*x-2, x) == Rational(2,3)

def test_quadratic():
    x = Symbol("x")
    assert solve(x**2-1, x) == [1, -1]
    assert solve(((x-1)*(x-2)).expand(), x) in [[1,2], [2,1]]
    assert solve(((x-1)*(x-1)).expand(), x) == [1]

def test_ODE_first_order():
    class f(Function):
        pass
    x = Symbol("x")
    assert dsolve(3*Derivative(f(x),x)-1, [f(x)]) == x/3+Symbol("C1")
    assert dsolve(x*Derivative(f(x),x)-1, [f(x)]) == log(abs(x))+Symbol("C1")

def test_ODE_second_order():
    class f(Function):
        pass
    x = Symbol("x")
    C1, C2 = Symbol("C1"), Symbol("C2")
    assert dsolve(Derivative(Derivative(f(x),x),x)+9*f(x), [f(x)]) == \
        sin(3*x)*C1+cos(3*x)*C2
