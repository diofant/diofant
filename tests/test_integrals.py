import sys
sys.path.append(".")

import py

from sym import Integral, Symbol, IntegralError, log

def test_basics():
    x=Symbol("x")
    t=Symbol("t")
    a=Symbol("a")
    e=(t+1)**2
    assert Integral(0,x,e,t).diff(x)==(1+x)**2
    assert Integral(0,x,e,t).diff(a)==0
    py.test.raises(IntegralError,"Integral(0,x,e,t).diff(t)")

    assert Integral(a,x,e,t).diff(x)==(1+x)**2
    assert Integral(a,x,e,t).diff(x)!=-(1+x)**2
    assert Integral(x,a,e,t).diff(x)==-(1+x)**2

    assert Integral(x,2*x,t**2,t).diff(x)==7*x**2

def test_integration():
    x=Symbol("x")
    t=Symbol("t")
    assert Integral(0,x, 0 ,t).doit()==0
    assert Integral(0,x, 3 ,t).doit()==3*x
    assert Integral(0,x, t ,t).doit()==x**2/2
    assert Integral(0,x, 3*t ,t).doit()==3*x**2/2
    assert Integral(0,x, 3*t**2 ,t).doit()==x**3
    assert Integral(1,x, -1/t**2 ,t).doit()==1/x-1
    assert Integral(1,x, 1/t ,t).doit()==log(x)

    assert Integral(0,x, t**2+5*t-8 ,t).doit()==x**3/3+5*x**2/2-8*x

    a=Symbol("a")
    b=Symbol("b")
    c=Symbol("c")
    assert Integral(0,x, a*t ,t).doit()==a*x**2/2
    assert Integral(0,x, a*t**4 ,t).doit()==a*x**5/5
    assert Integral(0,x, a*t**2 +b*t +c ,t).doit()==a*x**3/3+b*x**2/2+c*x
