import sys
sys.path.append(".")

from sympy import Rational, Symbol

def test_basics():
    x,y,a,b = [Symbol(Y) for Y in ["x","y","a","b"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = (x+y)**a
    assert e.match(p**q,[p,q]) == (x+y, a)
