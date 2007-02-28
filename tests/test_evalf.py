import sys
sys.path.append(".")

import sym as s

def eq(a,b):
    return abs(a-b)<0.0001

def testeval():
    e=s.log(3)/s.log(2)-1
    assert eq(e.evalf(),0.58496)

def test_bug1():
    x=s.Symbol('x')
    y=x*x
    assert eq(y.subs(x,s.Real(3.0)).evalf(),9)

def test_bug2():
    a = s.Real(4.)
    x = a + a
    x.eval()

