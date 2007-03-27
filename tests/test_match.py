import sys
sys.path.append(".")

from sympy import Rational, Symbol

def test_symbol():
    x,y,a,b,c = [Symbol(Y) for Y in ["x","y","a","b","c"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = a
    assert e.match(b,[p]) == None
    assert e.match(a,[p]) == {}
    assert e.match(p,[p]) == {p: a}
    assert e.match(p,[p]) != {p: b}

def xtest_add():
    x,y,a,b,c = [Symbol(Y) for Y in ["x","y","a","b","c"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = (a+b)
    assert e.match(p+b,[p]) == (a,)

def test_basics():
    x,y,a,b,c = [Symbol(Y) for Y in ["x","y","a","b","c"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = (x+y)**a
    assert e.match(p**q,[p,q]) == (x+y, a)
    assert e.match(p**p,[p]) == None
    e = (x+y)**(x+y)
    assert e.match(p**p,[p]) == (x+y,)
    assert e.match(p**q,[p,q]) == (x+y, x+y)
    #e = (a+b)*(a+c)
    #assert e.match((p+b)*(p+c),[p]) == (a,)
