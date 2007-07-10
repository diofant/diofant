import sys
sys.path.append(".")

from sympy import Rational, Symbol, symbols, sin, cos, solve, dsolve, Matrix
from sympy.core.functions import Function, Derivative, diff, exp, log
from sympy.modules.solvers import solve_linear_system, solve_linear_system_LU

def test_solve():
    x, y = symbols('x', 'y')

    assert solve(3*x-2, x) == [Rational(2,3)]
    assert solve(3*x == 2, x) == [Rational(2,3)]

    assert solve(x**2-1, x) == [-1, 1]
    assert solve(x**2 == 1, x) == [-1, 1]

def test_linear_system():
    x, y, z, t, n = symbols('x', 'y', 'z', 't', 'n')

    assert solve([x+5*y-2, -3*x+6*y-15], [x, y]) == {x: -3, y: 1}

    M = Matrix( [0,0,n*(n+1),(n+1)**2,0],
                [n+1,n+1,-2*n-1,-(n+1),0],
                [-1, 0, 1, 0, 0] )

    assert solve_linear_system(M, [x, y, z, t]) == \
           {y: 0, z: (-t-t*n)/n, x: (-t-t*n)/n}

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
