import sys
sys.path.append(".")

from sympy import Symbol
from sympy.modules.simplify import *

def test_ratsimp():
    x = Symbol("x")
    y = Symbol("y")
    e = 1/x+1/y
    assert e != (x+y)/(x*y)
    assert ratsimp(e) == (x+y)/(x*y)

    e = -x-y-(x+y)**(-1)*y**2+(x+y)**(-1)*x**2
    assert e != -2*y
    assert ratsimp(e) == -2*y

    e = x/(x+y)+y/(x+y)
    assert e != 1
    assert ratsimp(e) == 1