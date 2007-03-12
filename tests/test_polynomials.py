
import sys
sys.path.append(".")

from sympy import Symbol
from sympy.modules.polynomials import *

def test_ispoly():
    x = Symbol("x")
    y = Symbol("y")
    assert not ispoly( x.sqrt(), x )
    assert ispoly( Rational(2), x)
    assert ispoly( x, x)
    assert ispoly( x**2, x)
    assert ispoly( x**2 + 3*x - 8, x)
    assert ispoly( x**2 + 3*x*y.sqrt() - 8, x)
    assert not ispoly( x**2 + 3*x*y.sqrt() - 8 , y)

    #assert Rational(1).ispoly(sin(x))
    #assert not exp(x).ispoly(sin(x))
