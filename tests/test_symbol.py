import sys
sys.path.append(".")

import sympy as s

def testSymbol():
    a=s.Symbol("a")
    x1=s.Symbol("x")
    x2=s.Symbol("x")
    xdummy1=s.Symbol("x", is_dummy=True)
    xdummy2=s.Symbol("x", is_dummy=True)
    assert a != x1
    assert a != x2
    assert x1 == x2
    assert x1 != xdummy1
    assert xdummy1 != xdummy2

    assert s.Symbol("x")==s.Symbol("x")
    assert s.Symbol("x", is_dummy=True)!=s.Symbol("x", is_dummy=True)
