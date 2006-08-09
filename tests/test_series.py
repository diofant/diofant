import sys
sys.path.append(".")

import py

import sym as g

def testseries():
    n3=g.rational(3)
    n2=g.rational(2)
    n6=g.rational(6)
    x=g.symbol("x")
    c=g.symbol("c")
    e=g.sin(x)
    assert str(e) == "sin(x)"
    assert str(e.series(x,0)) == "0"
    assert str(e.series(x,1)) == "x"
    assert str(e.series(x,2)) == "x"
    assert e.series(x,3) == x+(-g.rational(1)/6)*x**3
    assert e.series(x,4) == x+(-g.rational(1)/6)*x**3

    e=((g.exp(x)-1)/x)
    assert e.series(x,1) == g.rational(1)
    py.test.raises(g.pole_error, g.basic.series, e,x,0)

    #e=2*g.sin(x)*g.cos(x)
    #print
    #print e.series(x,5)
    #e=g.sin(2*x)
    #e=g.tan(2*x)
    #e=1/g.cos(x)
    #print e.series(x,8)

def xtestfind(self):
    a=g.symbol("a")
    b=g.symbol("b")
    c=g.symbol("c")
    p=g.rational(5)
    e=a*b+b**p
    assert e.find(b)
    assert not e.find(c)

def testsubs():
    n3=g.rational(3)
    n2=g.rational(2)
    n6=g.rational(6)
    x=g.symbol("x")
    c=g.symbol("c")
    e=x
    assert str(e) == "x"
    e=e.subs(x,n3)
    assert str(e) == "3"

    e=2*x
    assert e == 2*x
    e=e.subs(x,n3)
    assert str(e) == "6"

    e=(g.sin(x)**2).diff(x)
    assert e == 2*g.sin(x)*g.cos(x)
    e=e.subs(x,n3)
    assert e == 2*g.cos(n3)*g.sin(n3)

    e=(g.sin(x)**2).diff(x)
    assert e == 2*g.sin(x)*g.cos(x)
    e=e.subs(g.sin(x),g.cos(x))
    assert e == 2*g.cos(x)**2
