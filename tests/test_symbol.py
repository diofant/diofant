import sys
sys.path.append(".")

import sym as s

def testSymbol():
    a=s.Symbol("a")
    x1=s.Symbol("x")
    x2=s.Symbol("x")
    xdummy1=s.Symbol("x",dummy=True)
    xdummy2=s.Symbol("x",dummy=True)
    assert a!=x1
    assert a!=x2
    assert x1==x2
    assert x1!=xdummy1
    assert xdummy1!=xdummy2

    assert s.Symbol("x")==s.Symbol("x")
    assert s.Symbol("x",True)!=s.Symbol("x",True)
