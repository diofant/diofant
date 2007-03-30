import sys
sys.path.append(".")

from sympy import Rational, Symbol, cos, solve, dsolve, Function, Derivative

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
