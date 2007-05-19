import sys
sys.path.append(".")

from sympy import Symbol, exp
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
    
    e = 1/(1+1/x)
    assert ratsimp(e) == x/(x+1)
    assert (x+1)*ratsimp(e)/x == 1
    assert ratsimp(exp(e)) == exp(x/(x+1))
    
def test_simplify():
    x = Symbol('x')
    y = Symbol('y')
    e = 1/x + 1/y
    assert e != (x+y)/(x*y)
    assert simplify(e) == (x+y)/(x*y)

    e = (4+4*x-2*(2+2*x))/(2+2*x)
    assert e != 0
    assert simplify(e) == 0
    
