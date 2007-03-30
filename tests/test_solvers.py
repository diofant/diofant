import sys
sys.path.append(".")

from sympy import Rational, Symbol, cos, solve

def test_lienar():
    x = Symbol("x")
    assert solve(3*x-2, x) == Rational(2,3)
    assert solve(x**2-1, x) == [1, -1]
    assert solve(((x-1)*(x-2)).expand(), x) in [[1,2], [2,1]]
    assert solve(((x-1)*(x-1)).expand(), x) == [1]
