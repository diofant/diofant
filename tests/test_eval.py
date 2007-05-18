import sys
sys.path.append(".")

import sympy as g
from sympy import *

def test_add_eval():
    a = Symbol("a")
    b = Symbol("b")
    
    e = a+b+a+b
    s1 = str(e)
    e.eval()
    s2 = str(e)
    assert s1 == s2

    c=g.Rational(1)
    p=g.Rational(5)
    
    e=a*b+c+p
    assert e == a*b+6
    
    e = c+a+p
    assert e == a+6
    
    e = c+a-p
    assert e == a+(-4)
    
    e = a+a
    assert e == 2*a
    
    e=a+p+a
    assert e == 2*a+5
    
    e = c+p
    assert e == g.Rational(6)
    
    e=b+a-b
    assert e == a
    
def test_addmul_eval():
    
    a = Symbol("a")
    b = Symbol("b")
    
    c=g.Rational(1)
    p=g.Rational(5)
    
    e = c+a+b*c+a-p
    assert e == 2*a+b+(-4)
    
    e=a*g.Rational(2)+p+a
    assert e == a*2+5+a
    
    e=a*g.Rational(2)+p+a
    assert e == 3*a+5
    
    e = a*g.Rational(2)+a
    assert e == 3*a
    
def test_pow_eval():
    
    assert (-1)**Rational(1,2) == I
    assert (-1)**Rational(1,3) == Rational(1,2)+Rational(1,2)*I*3**Rational(1,2)
    
    assert sqrt(-4) == 2*I
    assert sqrt( 4) == 2
    
    assert (8)**Rational(1,3) == 2
    
    assert sqrt(-2) == I*sqrt(2)
    assert (-1)**Rational(1,3) != I
    assert (-10)**Rational(1,3) != I*((10)**Rational(1,3))
    assert (-2)**Rational(1,4) != (2)**Rational(1,4)
    
def test_mulpow_eval():
    x = Symbol('x')
    
    assert sqrt(50)/(sqrt(2)*x) == 5/x
    assert sqrt(27)/sqrt(3) == 3
    
    
def test_evalpow_bug():
    x = Symbol("x")
    
    assert 1/(1/x) == x
    
    e=1/(-1/x)
    assert e==-x


def test_expbug():
    assert exp(-log(3))**(-1) == 3


def test_symbol_expand():
    x = Symbol('x')
    y = Symbol('y')

    f = x**4*y**4
    assert f == x**4*y**4
    assert f == f.expand()

    g = (x*y)**4
    assert f != g
    assert f == g.expand()
    assert g.expand() == g.expand().expand()
