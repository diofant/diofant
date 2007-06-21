import sys
sys.path.append(".")

from sympy import Rational, Symbol, cos, solve, dsolve, Function, diff, \
        log, sin, exp, Matrix
from sympy.core.functions import Derivative
from sympy.modules.solvers import solve_linear_system, solve_linear_system_LU

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
  
def test_linear_system():
    x, y, z, u = Symbol('x'), Symbol('y'), Symbol('z'), Symbol('t')
    
    assert solve([x+5*y-2, -3*x+6*y-15], [x, y]) == {x: -3, y: 1} 
    n = Symbol("n")
    M = Matrix( [0,0,n*(n+1), (n+1)**2,0], [n+1,n+1,-2*n-1,-(n+1),0],
                [-1, 0, 1, 0, 0])
    assert solve_linear_system(M, [x,y,z,u]) \
           ==  {y: 0, z: (-u-u*n)/n, x: (-u-u*n)/n}

def test_linear_systemLU():
    x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
    n = Symbol('n')
    M = Matrix( [1,2,0,1],[1,3,2*n,1],[4,-1,n**2,1])
    assert solve_linear_system_LU(M, [x,y,z]) == {z: -3/(n**2+18*n),
                                                  x: 1-12*n/(n**2+18*n),
                                                  y: 6*n/(n**2+18*n)}

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
