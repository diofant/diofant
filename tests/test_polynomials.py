
import sys
sys.path.append(".")

import py

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

def test_coeff():
    x = Symbol("x")
    assert coeff(x**2, x, 1) == 0
    assert coeff(x**2, x, 2) == 1
    assert coeff(x**2, x, 2) != 0

    assert coeff(2*x+18*x**8, x, 1) == 2
    assert coeff(2*x+18*x**8, x, 4) == 0
    assert coeff(2*x+18*x**8, x, 8) == 18

def test_get_poly():
    x = Symbol("x")
    y = Symbol("y")
    assert get_poly(3*x**2,x) == [(3,2)]
    assert get_poly(2*x+3*x**2 - 5,x) == [(-5, 0), (2, 1), (3,2)]
    assert get_poly(2*x**100+3*x**2 - 5,x) == [(-5, 0), (3,2), (2, 100)]

    assert get_poly(y.sqrt()*x,x) == [(y.sqrt(),1)]
    assert get_poly(x**2 + 3*x*y.sqrt() - 8, x) == [(-8, 0), (3*y.sqrt(), 1), 
            (1, 2)]
    py.test.raises(PolynomialException, "get_poly(x.sqrt(),x)")

def test_poly():
    x = Symbol("x")
    y = Symbol("y")
    assert 3*x**2 == poly([(3,2)],x)
    assert 2*x+3*x**2 - 5 == poly([(-5, 0), (2, 1), (3,2)],x)
    assert 2*x**100+3*x**2 - 5 == poly([(-5, 0), (3,2), (2, 100)],x)
    assert 2*x**100+3*x**2 - 6 != poly([(-5, 0), (3,2), (2, 100)],x)

    assert y.sqrt()*x == poly([(y.sqrt(),1)],x)
    assert x**2 + 3*x*y.sqrt() - 8 == poly([(-8, 0), (3*y.sqrt(), 1), 
        (1, 2)],x)

def test_gcd():
    x = Symbol("x")
    assert gcd(x**2, x, x) == x
    assert gcd(3*x**2, x, x) == x
    assert gcd(3*x**2, 3*x, x) == 3*x
    assert gcd(x**2+2*x+1, x+1, x) == x+1
    assert gcd(x**2+2*x+2, x+1, x) == 1

def test_rep():
    assert rep(101,100) == (1,1)
    assert rep(300,100) == (0,3)
    assert rep(100,100) == (0,1)

    assert rep(100,10) == (0,0,1)

def test_sqf():
    x = Symbol("x")
    assert sqf(3*x**2, x) == 3*x**2
    assert sqf(x**2+2*x+1, x) == (x+1)**2
