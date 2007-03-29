import sys
sys.path.append(".")

import sympy as g
import sympy as s
from sympy import Symbol, log, arctan
from sympy.core.functions import Function
#from sympy.modules.derivatives import Derivative
from sympy.core.functions import Derivative

def test_func():
    a=g.Symbol("a")
    b=g.Symbol("b")
    c=g.Symbol("c")
    p=g.Rational(5)
    e=a*b+g.sin(b**p)
    assert e == a*b+g.sin(b**5)
    assert e.diff(a) == b
    assert e.diff(b) == a+5*b**4*g.cos(b**5)
    e=g.tan(c)
    assert e == g.tan(c)
    assert e.diff(c) in [g.cos(c)**(-2),1+g.sin(c)**2/g.cos(c)**2]
    e=c*g.log(c)-c
    assert e == -c+c*g.log(c)
    assert e.diff(c) == g.log(c)
    e=g.log(g.sin(c))
    assert e == g.log(g.sin(c))
    assert e.diff(c) == g.sin(c)**(-1)*g.cos(c)
    assert e.diff(c) != g.cos(c)**(-1)*g.sin(c)
    assert e.diff(c) != g.sin(c)**(-2)*g.cos(c)
    assert e.diff(c) != g.sin(c)**(-3)*g.cos(c)
    t=g.Rational(2)
    e=(t**a/g.log(t))
    assert e == 2**a*g.log(g.Rational(2))**(-1)
    assert e.diff(a) == 2**a

def test_log():
    assert g.log(2) > 0

def test_exp_log():
    x=g.Symbol("x")
    assert g.log(g.exp(x))==x
    assert g.exp(g.log(x))==x

def test_log_expansion():
    x=g.Symbol("x")
    y=g.Symbol("y")
    assert g.log(x*y)==g.log(x)+g.log(y)
    assert g.log(x**2)==2*g.log(x)

def test_log_hashing_bug():
    x=s.Symbol("y")
    assert x!=s.log(s.log(x))
    assert x.hash()!=s.log(s.log(x)).hash()
    assert s.log(x)!=s.log(s.log(s.log(x)))

    e=1/s.log(s.log(x)+s.log(s.log(x)))
    e=e.eval()
    assert isinstance(e.base,s.log)
    e=1/s.log(s.log(x)+s.log(s.log(s.log(x))))
    e=e.eval()
    assert isinstance(e.base,s.log)

    x=s.Symbol("x")
    e=s.log(s.log(x))
    assert isinstance(e,s.log)
    assert not isinstance(x,s.log)
    assert s.log(s.log(x)).hash() != x.hash()
    assert e!=x

def test_sign():
    assert s.sign(s.log(2)) == 1


def test_exp_bug():
    x=s.Symbol("x")
    assert s.exp(1*s.log(x))==x

def test_exp_expand():
    x=s.Symbol("x")
    e=s.exp(s.log(s.Rational(2))*(1+x)-s.log(s.Rational(2))*x)
    assert e.expand()==2

def test_pi():
    assert s.cos(s.pi)==-1
    assert s.cos(2*s.pi)==1
    assert s.sin(s.pi)==0
    assert s.sin(2*s.pi)==0

def test_bug1():
    x=Symbol("x")
    w=Symbol("w")
    e=(-log(w)).sqrt()
    assert e.subs(log(w),-x)!=-x.sqrt()
    assert e.subs(log(w),-x)==x.sqrt()

    e=(-5*log(w)).sqrt()
    assert e.subs(log(w),-x)==(5*x).sqrt()

def test_Derivative():
    x=Symbol("x")
    e=Derivative(log(x),x)
    assert e!=1/x
    assert e.doit()==1/x

def test_invtrig():
    x=Symbol("x")
    assert arctan(0) == 0
    assert arctan(x).diff(x) == 1/(1+x**2)

def test_general_function():
    class nu(Function):
        pass

    x=Symbol("x")
    y=Symbol("y")
    e=nu(x)
    edx=e.diff(x)
    edy=e.diff(y)
    edxdx=e.diff(x).diff(x)
    edxdy=e.diff(x).diff(y)
    assert e == nu(x)
    assert edx != nu(x)
    assert edx == Derivative(nu(x), x)
    assert edy == 0
    assert edxdx == Derivative(Derivative(nu(x), x), x)
    assert edxdy == 0

    #this works, but is semantically wrong, we need to settle on some interface
    #first
    assert nu(x**2).diff(x) == Derivative(nu(x**2), x**2) * 2*x
