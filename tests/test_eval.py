import sys
sys.path.append(".")

import sym as g

def testeval():
    a=g.symbol("a")
    b=g.symbol("b")
    c=g.rational(1)
    p=g.rational(5)
    e=(a*b+c+p).eval()
    assert e == a*b+6
    e=(c+a+p).eval()
    assert e == a+6
    e=(c+a-p).eval()
    assert e == a+(-4)
    e=(c+a+b*c+a-p).eval()
    assert e == 2*a+b+(-4)
    e=(a+a).eval()
    assert e == 2*a
    e=(a+p+a).eval()
    assert e == 2*a+5
    e=(a*g.rational(2)+p+a)
    assert e == a*2+5+a
    e=(a*g.rational(2)+p+a).eval()
    assert e == 3*a+5
    e=(c+p).eval()
    assert e == g.rational(6)
    e=(a*g.rational(2)+a).eval()
    assert e == 3*a
    e=b+a-b
    assert e == a

def testeval():
    a=g.symbol("a")
    b=g.symbol("b")
    e=a+b+a+b
    s1=str(e)
    e.eval()
    s2=str(e)
    assert s1 == s2
