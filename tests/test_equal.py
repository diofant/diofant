import sys
sys.path.append(".")

import sym as g

def testequal():
    b=g.symbol("b")
    a=g.symbol("a")
    c=g.symbol("c")
    e1=a+b
    e2=2*a*b
    e3=a**3*b**2
    e4=a*b+b*a
    assert not e1.isequal(e2)
    assert not e1==e2
    assert e1!=e2
    assert e2==e4
    assert e2!=e3
    assert not e2==e3

    x=g.symbol("x")
    e1=g.exp(x+1/x)
    y=g.symbol("x")
    e2=g.exp(y+1/y)
    assert e1==e2
    assert not e1!=e2
    y=g.symbol("y")
    e2=g.exp(y+1/y)
    assert not e1==e2
    assert e1!=e2

    e5=g.rational(3)+2*x-x-x
    assert e5==3
    assert 3==e5
    assert e5!=4
    assert 4!=e5
    assert e5!=3+x
    assert 3+x!=e5

def test_expevalbug():
    x=g.symbol("x")
    e1=g.exp(1*x)
    h1=e1.hash()
    e2=e1.eval()
    e3=g.exp(x)
    assert e1==e2
    assert e1==e3
    assert e2==e3
