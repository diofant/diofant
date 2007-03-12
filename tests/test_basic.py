import sys
sys.path.append(".")

import sympy as g
from sympy import Symbol, Rational, sin, exp

def dotest(s):
    x = g.Symbol("x")
    y = g.Symbol("y")
    l = [
    g.Rational(2),
    g.Real("1.3"), 
    x,
    y,
    pow(x,y)*y,
    5,
    5.5
    ]
    for x in l:
        for y in l:
            s(x,y)

def test_basic():
    def s(a,b):
        x = a
        x = +a
        x = -a
        x = a+b
        x = a-b
        x = a*b
        x = a/b
        x = a**b
    dotest(s)

def test_ibasic():
    def s(a,b):
        x = a
        x += b
        x = a
        x -= b
        x = a
        x *= b
        x = a
        x /= b
    dotest(s)

def test_ldegree():
    x=g.Symbol("x")
    assert (1/x**2+1+x+x**2).ldegree(x)==-2
    assert (1/x+1+x+x**2).ldegree(x)==-1
    assert (x**2+1/x).ldegree(x)==-1
    assert (1+x**2).ldegree(x)==0
    assert (x+1).ldegree(x)==0
    assert (x+x**2).ldegree(x)==1
    assert (x**2).ldegree(x)==2

def test_leadterm():
    x=g.Symbol("x")
    log=g.log
    assert (3+2*x**(log(3)/log(2)-1)).eval().leadterm(x)==(3,0)


def test_print_tree():
    x=g.Symbol("x") 
    y=g.Symbol("y") 

    e=(2*x-(7*x**2 - 2) + 3*y)
    e.print_tree()

def test_ispoly():
    x = Symbol("x")
    y = Symbol("y")
    assert not ( x.sqrt() ).ispoly(x)
    assert ( Rational(2) ).ispoly(x)
    assert ( x ).ispoly(x)
    assert ( x**2 ).ispoly(x)
    assert ( x**2 + 3*x - 8 ).ispoly(x)
    assert ( x**2 + 3*x*y.sqrt() - 8 ).ispoly(x)
    assert not ( x**2 + 3*x*y.sqrt() - 8 ).ispoly(y)

    #assert Rational(1).ispoly(sin(x))
    #assert not exp(x).ispoly(sin(x))
