import sys
sys.path.append(".")

from sympy import Rational, Symbol, cos, solve

def test_lienar():
    x = Symbol("x")
    eq = 3*x-2
    solve(eq, x)
