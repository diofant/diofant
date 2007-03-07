import sys
sys.path.append(".")

import py

from sym import Integral, Symbol, IntegralError, log

def test_basics():
    x=Symbol("x")
    t=Symbol("t")
    a=Symbol("a")
    e=(t+1)**2
    assert Integral(e, (t,0,x)).diff(x)==(1+x)**2
    assert Integral(e, (t,0,x)).diff(a)==0
    py.test.raises(IntegralError,"Integral(e,(t,0,x)).diff(t)")

    assert Integral(e, (t,a,x)).diff(x)==(1+x)**2
    assert Integral(e, (t,a,x)).diff(x)!=-(1+x)**2
    assert Integral(e, (t,x,a)).diff(x)==-(1+x)**2

    assert Integral(t**2, (t,x,2*x)).diff(x)==7*x**2

def test_integration():
    x=Symbol("x")
    t=Symbol("t")
    assert Integral(0, (t,0,x)).doit()==0
    assert Integral(3, (t,0,x)).doit()==3*x
    assert Integral(t, (t,0,x)).doit()==x**2/2
    assert Integral(3*t, (t,0,x)).doit()==3*x**2/2
    assert Integral(3*t**2, (t,0,x)).doit()==x**3
    assert Integral(-1/t**2, (t,1,x)).doit()==1/x-1
    assert Integral(1/t, (t,1,x)).doit()==log(x)

    assert Integral(t**2+5*t-8, (t,0,x)).doit()==x**3/3+5*x**2/2-8*x

    a=Symbol("a")
    b=Symbol("b")
    c=Symbol("c")
    assert Integral(a*t, (t,0,x)).doit()==a*x**2/2
    assert Integral(a*t**4, (t,0,x)).doit()==a*x**5/5
    assert Integral(a*t**2+b*t+c, (t,0,x)).doit()==a*x**3/3+b*x**2/2+c*x
