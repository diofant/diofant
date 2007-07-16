import sys
sys.path.append(".")

import py

import sympy as g
from sympy import Basic, Symbol, Rational, Add, Mul, sin, exp

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
    assert (3+2*x**(log(3)/log(2)-1)).leadterm(x)==(3,0)

def test_print_tree():
    x=g.Symbol("x")
    y=g.Symbol("y")

    e=(2*x-(7*x**2 - 2) + 3*y)
    e.print_tree()

def test_atoms():
   x = Symbol('x')
   y = Symbol('y')
   assert (1+x).atoms() == [1,x]
   assert x.atoms() == [x]
   assert (1+2*g.cos(x)).atoms() == [1,2,x]
   assert (2*(x**(y**x))).atoms() == [2,x,y]
   assert g.Rational(1,2).atoms() == [g.Rational(1,2)]

   assert g.Rational(1,2).atoms(type=(g.core.numbers.Infinity)) == []

def test_has_class():
    expr = sin(4+exp(Symbol('x')))
    assert expr.has_class(sin)
    assert expr.has_class(Add)
    assert expr.has_class(Rational)
    assert expr.has_class(exp)
    assert expr.has_class(Symbol)
    assert not expr.has_class(Mul)

def test_assumptions():
    x = Basic(is_integer=True)

    assert x.is_integer == True
    assert x.is_real == True

    x = Basic(is_nonnegative=True)

    assert x.is_nonnegative == True
    assert x.is_negative == False
    assert x.is_positive == None

    x = Basic(is_nonpositive=True)

    assert x.is_nonpositive == True
    assert x.is_positive == False
    assert x.is_negative == None

    x = Basic(is_negative=True)

    assert x.is_negative == True
    assert x.is_positive == False
    assert x.is_nonnegative == False

    x = Basic(is_positive=True)

    assert x.is_positive == True
    assert x.is_negative == False
    assert x.is_nonpositive == False

    x = Basic(is_integer=True, is_nonnegative=True)

    assert x.is_nonnegative_integer == True
    assert x.is_negative_integer == False
    assert x.is_positive_integer == None

    x = Basic(is_integer=True, is_nonpositive=True)

    assert x.is_nonpositive_integer == True
    assert x.is_positive_integer == False
    assert x.is_negative_integer == None

    x = Basic(is_odd=True)

    assert x.is_odd == True
    assert x.is_even == False
    assert x.is_integer == True

    x = Basic(is_odd=False)

    assert x.is_odd == False
    assert x.is_even == None
    assert x.is_integer == None

    x = Basic(is_even=True)

    assert x.is_even == True
    assert x.is_odd == False
    assert x.is_integer == True

    x = Basic(is_even=False)

    assert x.is_even == False
    assert x.is_odd == None
    assert x.is_integer == None

    x = Basic(is_prime=True)

    assert x.is_prime == True
    assert x.is_integer == True
    assert x.is_positive == True
    assert x.is_negative == False
    assert x.is_nonpositive == False
    assert x.is_nonnegative == None

    x = Basic(is_prime=False)

    assert x.is_prime == False
    assert x.is_integer == None
    assert x.is_positive == None
    assert x.is_negative == None
    assert x.is_nonpositive == None
    assert x.is_nonnegative == None

    x = Basic(is_nni=True)

    assert x.is_integer == True
    assert x.is_nonnegative == True

    x = Basic(is_npi=True)

    assert x.is_integer == True
    assert x.is_nonpositive == True

    py.test.raises(AttributeError, "x.is_real = False")

