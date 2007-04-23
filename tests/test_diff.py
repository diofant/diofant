import sys
sys.path.append(".")

import sympy as g
from sympy.core.functions import diff

def testdiff():
    a=g.Symbol("a")
    b=g.Symbol("b")
    c=g.Symbol("c")
    p=g.Rational(5)
    e=a*b+b**p
    assert e.diff(a) == b
    assert e.diff(b) == a+5*b**4
    assert e.diff(b).diff(a) == g.Rational(1)
    e=a*(b+c)
    assert e.diff(a) == b+c
    assert e.diff(b) == a
    assert e.diff(b).diff(a) == g.Rational(1)
    e=c**p
    assert diff(e, c,6) == g.Rational(0)
    assert diff(e, c,5) == g.Rational(120)
    e=c**g.Rational(2)
    assert e.diff(c) == 2*c
    e=(a*b*c)
    assert e.diff(c) == a*b

def testdiff2():
    n3=g.Rational(3)
    n2=g.Rational(2)
    n6=g.Rational(6)
    x=g.Symbol("x")
    c=g.Symbol("c")
    e=n3*(-n2+x**n2)*g.cos(x)+x*(-n6+x**n2)*g.sin(x)
    assert e == 3*((-2)+x**2)*g.cos(x)+x*((-6)+x**2)*g.sin(x)
    assert e.diff(x).expand() == x**3*g.cos(x)

    e=(x+1)**3
    assert e.diff(x) == 3*(x+1)**2
    e=x*(x+1)**3
    assert e.diff(x) == (x+1)**3+3*x*(x+1)**2
    e=(2*g.exp(x*x)*x)
    assert e.diff(x) == 2*g.exp(x*x)+4*x**2*g.exp(x*x)
