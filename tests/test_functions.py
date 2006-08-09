import sys
sys.path.append(".")

import sym as g

def testfunc():
    a=g.symbol("a")
    b=g.symbol("b")
    c=g.symbol("c")
    p=g.rational(5)
    e=a*b+g.sin(b**p)
    assert e == a*b+g.sin(b**5)
    assert e.diff(a) == b
    assert e.diff(b) == a+5*b**4*g.cos(b**5)
    e=g.tan(c)
    assert e == g.tan(c)
    assert e.diff(c) == g.cos(c)**(-2)
    e=c*g.ln(c)-c
    assert e == -c+c*g.ln(c)
    assert e.diff(c) == g.ln(c)
    e=g.ln(g.sin(c))
    assert e == g.ln(g.sin(c))
    assert e.diff(c) == g.sin(c)**(-1)*g.cos(c)
    assert e.diff(c) != g.cos(c)**(-1)*g.sin(c)
    assert e.diff(c) != g.sin(c)**(-2)*g.cos(c)
    assert e.diff(c) != g.sin(c)**(-3)*g.cos(c)
    t=g.rational(2)
    e=(t**a/g.ln(t))
    assert e == 2**a*g.ln(g.rational(2))**(-1)
    assert e.diff(a) == 2**a
