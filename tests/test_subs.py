import sys
sys.path.append(".")

import py

import sym as g

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

def test_lnexppow():
    x=g.symbol("x")
    w=g.symbol("dummy :)")
    e=(3**(1+x)+2**(1+x))/(3**x+2**x)
    e=e.eval()
    assert e.subs(2**x,w)!=e
    assert e.subs(g.exp(x*g.ln(g.rational(2))),w)!=e
