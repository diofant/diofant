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
    
    e = x
    assert e.match(x,[b,c]) == {}

    e = Rational(5)
    assert e.match(c, [c]) == {c: 5}
    assert e.match(c, [b,c]) == {c: 5}
    assert e.match(e, [b]) == {}
    assert e.match(e, [b,c]) == {}
    assert e.match(e+1, [b,c]) == None

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
    assert e.match(b,[p]) == None
    assert e.match(b+p,[p]) == {p: a+c+x}
    assert e.match(a+p+c,[p]) == {p: b+x}
    assert e.match(b+p+c,[p]) == {p: a+x}
    e = 4*x+5
    assert e.match(4*x+c,[c]) == {c: 5}
    assert e.match(3*x+c,[c]) == None
    assert e.match(b*x+5,[b]) == {b: 4}
    assert e.match(b*x+c,[b,c]) == {b: 4, c: 5}
    e = 4*x+5*y+6
    assert e.match(a*x+b*y+c,[a,b,c]) == {a: 4, b: 5, c: 6}

def test_power():
    x,y,a,b,c = [Symbol(Y) for Y in ["x","y","a","b","c"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = (x+y)**a
    assert e.match(p**q,[p,q]) == {p: x+y, q: a}
    assert e.match(p**p,[p]) == None
    e = (x+y)**(x+y)
    assert e.match(p**p,[p]) == {p: x+y}
    assert e.match(p**q,[p,q]) == {p: x+y, q: x+y}

    #e = 3/(4*x+5)
    #assert e.match(a/(b*x+c),[a,b,c]) == {a: 3, b: 4, c: 5}

def test_mul():
    x,y,a,b,c = [Symbol(Y) for Y in ["x","y","a","b","c"]]
    p,q = [Symbol(Y) for Y in ["p","q"]]
    e = 4*x
    assert e.match(b*x,[b]) == {b: 4}
    assert e.match(b*x,[b,c]) == {b: 4}
    assert e.match(b*y,[b]) == None
    assert e.match(b*y,[a,b,c]) == None
    e = a*x*b*c
    assert e.match(p*x,[p]) == {p: a*b*c}
    assert e.match(c*p*x,[p]) == {p: a*b}
    e = (a+b)*(a+c)
    #assert e.match((p+b)*(p+c),[p]) == {p: a}
