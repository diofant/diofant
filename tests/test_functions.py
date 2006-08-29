import sys
sys.path.append(".")

import sym as g
import sym as s

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

def testexpln():
    x=g.symbol("x")
    assert g.ln(g.exp(x))==x
    assert g.exp(g.ln(x))==x

def testlnexpansion():
    x=g.symbol("x")
    y=g.symbol("y")
    assert g.ln(x*y)==g.ln(x)+g.ln(y)
    assert g.ln(x**2)==2*g.ln(x)

def testlnhashingbug():
    x=s.symbol("y")
    assert x!=s.ln(s.ln(x))
    assert s.ln(x)!=s.ln(s.ln(s.ln(x)))

    e=1/s.ln(s.ln(x)+s.ln(s.ln(x)))
    e=e.eval()
    assert isinstance(e.a,s.ln)
    e=1/s.ln(s.ln(x)+s.ln(s.ln(s.ln(x))))
    e=e.eval()
    assert isinstance(e.a,s.ln)

def testexpbug():
    x=s.symbol("x")
    assert s.exp(1*s.ln(x))==x

def testexpexpand():
    x=s.symbol("x")
    e=s.exp(s.ln(s.rational(2))*(1+x)-s.ln(s.rational(2))*x)
    assert e.expand()==2
