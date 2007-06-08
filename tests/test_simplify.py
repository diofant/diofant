import sys
sys.path.append(".")

from sympy import Symbol, exp, cos, sin, tan, sec, csc, cot
from sympy.core.numbers import Rational
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

def test_trigsimp():
    x = Symbol('x')
    y = Symbol('y')
    assert trigsimp(5*cos(x)**2 + 5*sin(x)**2) == 5
    assert trigsimp(5*cos(x/2)**2 + 2*sin(x/2)**2) == 2 + 3*cos(x/2)**2
    assert trigsimp(1 + tan(x)**2) == sec(x)**2
    assert trigsimp(1 + cot(x)**2) == csc(x)**2

def test_simplify():
    x = Symbol('x')
    y = Symbol('y')
    e = 1/x + 1/y
    assert e != (x+y)/(x*y)
    assert simplify(e) == (x+y)/(x*y)

    e = (4+4*x-2*(2+2*x))/(2+2*x)
    assert e != 0
    assert simplify(e) == 0

def test_fraction():
    x, y = Symbol('x'), Symbol('y')

    assert fraction(Rational(1, 2)) == (1, 2)

    assert fraction(x) == (x, 1)
    assert fraction(1/x) == (1, x)
    assert fraction(x/y) == (x, y)

    assert fraction((x**2+1)/y) == (x**2+1, y)
    assert fraction(x*(y+1)/y**7) == (x*(y+1), y**7)

def test_together():
    x, y = Symbol('x'), Symbol('y')

    assert together(1/x) == 1/x

    assert together(1/x + 1) == (x+1)/x
    assert together(1/x + x) == (x**2+1)/x

    assert together(1/x + Rational(1, 2)) == (x+2)/(2*x)

    assert together(1/x + 2/y) == (2*x+y)/(y*x)
    assert together(1/(1 + 1/x)) == x/(1+x)
    assert together(x/(1 + 1/x)) == x**2/(1+x)
