import sys
sys.path.append(".")

import sym as g

def testdiff():
    a=g.symbol("a")
    b=g.symbol("b")
    c=g.symbol("c")
    p=g.rational(5)
    e=a*b+b**p
    assert e.diff(a) == b
    assert e.diff(b) == a+5*b**4
    assert e.diff(b).diff(a) == g.rational(1)
    e=a*(b+c)
    assert e.diff(a) == b+c
    assert e.diff(b) == a
    assert e.diff(b).diff(a) == g.rational(1)
    e=c**p
    assert e.diffn(c,6) == g.rational(0)
    assert e.diffn(c,5) == g.rational(120)
    e=c**g.rational(2)
    assert e.diff(c) == 2*c
    e=(a*b*c)
    assert e.diff(c) == a*b

def testdiff2():
    n3=g.rational(3)
    n2=g.rational(2)
    n6=g.rational(6)
    x=g.symbol("x")
    c=g.symbol("c")
    e=n3*(-n2+x**n2)*g.cos(x)+x*(-n6+x**n2)*g.sin(x)
    assert e == 3*((-2)+x**2)*g.cos(x)+x*((-6)+x**2)*g.sin(x)
    assert e.diff(x).expandterms() == x**3*g.cos(x)

    e=(x+1)**3
    assert e.diff(x) == 3*(x+1)**2
    e=x*(x+1)**3
    assert e.diff(x) == (x+1)**3+3*x*(x+1)**2
    e=(2*g.exp(x*x)*x)
    assert e.diff(x) == 2*g.exp(x*x)+4*x**2*g.exp(x*x)
