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

def test_add():
    x,y,a,b,c = [Symbol(Y) for Y in ["x","y","a","b","c"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = a+b
    assert e.match(p+b,[p]) == {p: a}
    assert e.match(p+a,[p]) == {p: b}
    e = a+b+c
    assert e.match(a+p+c,[p]) == {p: b}
    assert e.match(b+p+c,[p]) == {p: a}
    e = a+b+c+x
    assert e.match(a+p+x+c,[p]) == {p: b}
    assert e.match(b+p+c+x,[p]) == {p: a}
    #e = a+b+c+x
    #assert e.match(a+p+c,[p]) == {p: b+x}
    #assert e.match(b+p+c,[p]) == {p: a+x}

def test_basics():
    x,y,a,b,c = [Symbol(Y) for Y in ["x","y","a","b","c"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = (x+y)**a
    assert e.match(p**q,[p,q]) == (x+y, a)
    assert e.match(p**p,[p]) == None
    e = (x+y)**(x+y)
    assert e.match(p**p,[p]) == (x+y,)
    assert e.match(p**q,[p,q]) == (x+y, x+y)
    e = (a+b)*(a+c)
    assert e.match((p+b)*(p+c),[p]) == {p: a}
