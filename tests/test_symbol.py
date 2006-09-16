import sys
sys.path.append(".")

import sym as s

def testsymbol():
    a=s.symbol("a")
    x1=s.symbol("x")
    x2=s.symbol("x")
    xdummy1=s.symbol("x",dummy=True)
    xdummy2=s.symbol("x",dummy=True)
    assert a!=x1
    assert a!=x2
    assert x1==x2
    assert x1!=xdummy1
    assert xdummy1!=xdummy2

    assert s.symbol("x")==s.symbol("x")
    assert s.symbol("x",True)!=s.symbol("x",True)
