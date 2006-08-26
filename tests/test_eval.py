import sys
sys.path.append(".")

import sym as g

def testeval():
    a=g.symbol("a")
    b=g.symbol("b")
    e=a+b+a+b
    s1=str(e)
    e.eval()
    s2=str(e)
    assert s1 == s2

    c=g.rational(1)
    p=g.rational(5)
    e=a*b+c+p
    assert e == a*b+6
    e=c+a+p
    assert e == a+6
    e=c+a-p
    assert e == a+(-4)
    e=c+a+b*c+a-p
    assert e == 2*a+b+(-4)
    e=a+a
    assert e == 2*a
    e=a+p+a
    assert e == 2*a+5
    e=a*g.rational(2)+p+a
    assert e == a*2+5+a
    e=a*g.rational(2)+p+a
    assert e == 3*a+5
    e=c+p
    assert e == g.rational(6)
    e=a*g.rational(2)+a
    assert e == 3*a
    e=b+a-b
    assert e == a

def test_evalpow_bug():
    x=g.symbol("x")
    e=1/(1/x)
    assert e==x
    e=1/(-1/x)
    print e.eval()
    assert e==-x
