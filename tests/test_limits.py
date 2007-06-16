import sys
sys.path.append(".")

from sympy.modules.limits import Limit
from sympy import log, limit, atan, Symbol, oo, pi, Rational
from sympy import sin, cos, exp

def testsets():
    from sympy.modules.limits import member,intersect,union
    x=[1,2,3]
    y=[2,4]
    assert member(2,y)
    assert not member(3,y)
    assert intersect(x,y)
    assert not intersect(x,[4,5])
    assert union(x,y)==[1,2,3,4]
    assert union(x,[1,3])==[1,2,3]
    assert union(x,[4,5])==[1,2,3,4,5]

def testcompare():
    from sympy.modules.limits import compare
    x=Symbol("y")
    assert compare(exp(x),x**5,x) == ">"
    assert compare(exp(x**2),exp(x)**2,x) == ">"
    assert compare(exp(x),exp(x+exp(-x)),x) == "="
    assert compare(exp(x+exp(-x)),exp(x),x) == "="
    assert compare(exp(x+exp(-x)),exp(-x),x) == "="
    assert compare(exp(exp(x)),exp(x+exp(-exp(x))),x) == ">"
    assert compare(exp(-x),x,x) ==  ">"
    assert compare(x,exp(-x),x) ==  "<"
    assert compare(exp(x+1/x),x,x) == ">"
    assert compare(exp(exp(x)),exp(x+exp(-exp(x))),x) == ">"
    assert compare(exp(-exp(x)),exp(x),x) == ">"
    assert compare(exp(exp(-exp(x))+x),exp(-exp(x)),x) == "<"

def eq(a,b):
    assert len(a)==len(b)
    for x,y in zip(a,b):
        assert x==y

def testmax():
    from sympy.modules.limits import max
    x=Symbol("y")
    eq(max([exp(x)],[x**5],x),  [exp(x)])
    eq(max([exp(-x)],[x],x),  [exp(-x)])


def testmrv():
    from sympy.modules.limits import mrv
    x=Symbol("y")
    eq(mrv(exp(x+1/x),x),[exp(x+1/x)])
    eq(mrv(-exp(1/x),x),[x])
    eq(mrv(x,x),[x])
    eq(mrv(exp(-x),x),[exp(-x)])
    eq(mrv(exp(x+exp(-x)),x),[exp(x+exp(-x)),exp(-x)])
    eq(mrv(exp(x+exp(-exp(x))),x),[exp(-exp(x))] )
    eq(mrv(exp(x+exp(-x**2)),x),[exp(-x**2)] )

def test_simple_limit_manual():
    "example 3.15"
    from sympy.modules.limits import mrv
    x = Symbol("y")
    f = (exp(1/x-exp(-x))-exp(1/x))/exp(-x)
    Omega=mrv(f,x)
    assert Omega==[exp(-x)]
    assert Omega!=[exp(x)]
    wexpr=Omega[0]
    w = Symbol("w")
    f2 = f.subs(wexpr,w)
    ser = f2.series(w,3)
    lterm=ser.leadterm(w)
    assert lterm[0]==-exp(1/x)
    assert lterm[1]==0
    Omega=mrv(lterm[0],x)
    assert Omega==[x]

def test_simple_limit_lessmanual():
    "example 3.15"
    x = Symbol("y")
    f = (exp(1/x-exp(-x))-exp(1/x))/exp(-x)
    from sympy.modules.limits import mrv_leadterm
    lterm = mrv_leadterm(f,x)
    assert lterm[0]==-exp(1/x)
    assert lterm[1]==0

def test_simple_limit_automatic():
    "example 3.15"
    x = Symbol("y")
    f = (exp(1/x-exp(-x))-exp(1/x))/exp(-x)
    assert limit(f,x,oo) == -1

def testlimitinf_lenmrveq1():
    x = Symbol("y")
    assert limit(x,x,oo) == oo
    assert limit(-x,x,oo) == oo
    assert limit(-x**2,x,oo) == oo
    assert limit(1/x,x,oo) == 0
    assert limit(1/x,x,oo) != 1
    assert limit(1+1/x,x,oo) == 1
    assert limit(1+1/x,x,oo) != 0
    assert limit(exp(x),x,oo) == oo
    assert limit(-exp(x),x,oo) == oo
    assert limit(-exp(1/x),x,oo) == -1
    assert limit(exp(x)/x,x,oo) == oo
    assert limit(exp(x)/x,x,oo) != 1
    assert limit(x+exp(-x),x,oo) == oo
    assert limit(x+exp(-x**2),x,oo) == oo
    assert limit(x+exp(-exp(x)),x,oo) == oo
    assert limit(1/x-exp(-x),x,oo) == 0
    assert limit(13+1/x-exp(-x),x,oo) == 13
    assert limit(x+1/x,x,oo) == oo

def testlimit():
    x = Symbol("y")
    e = exp(Rational(1))
    assert limit((exp(x)-1)/x,x,0) == 1
    assert limit(exp(x),x,0) == 1
    assert limit(exp(x),x,1) == e
    assert limit(exp(x),x,-1) == exp(Rational(-1))
    assert limit(log(x)*x,x,0) == 0

def testlimitinf_lenmrveq2():
    x = Symbol("y")
    assert limit(exp(x+exp(-x))-exp(x), x, oo) == 1
    assert limit(1/exp(-x+exp(-x))-exp(x), x, oo) == -1
    #example 8.19
    #e=(log(log(x)+log(log(x)))-log(log(x)))/log(log(x)+log(log(log(x))))*log(x)
    #assert limit(e,x,oo)==1

def xtestlonglimit1():
    "example 8.18"
    x = Symbol("y")
    h = exp(-x/(1+exp(-x)))
    e = (exp(h)*exp(-x/(1+h))*exp(exp(-x+h)))/h**2-exp(x)+x
    l = limit(e, x, oo)
    print "limit=",l
    assert l== 2
    assert l!= 1

def testlog():
    x=Symbol("x")
    assert limit(log(x),x,0)==-oo
    assert limit(log(x),x,oo)==oo

def testsubexp():
    x=Symbol("x")
    e=log(x)
    from sympy.modules.limits import subexp
    assert subexp(e,x)
    assert not subexp(e,exp(x))
    assert subexp(exp(exp(x)+log(x)),x)
    assert subexp(exp(exp(x)+log(x)),log(x))
    assert subexp(exp(exp(x)+log(x)),exp(x))
    assert not subexp(exp(exp(x)+log(x)),2*x)

def sqrt(x):
    return x**Rational(1,2)

def test_functions():
    x=Symbol("x")
    assert limit(sin(x)/x,x,0) == 1
    assert limit(cos(x)/sin(x),x,0) == oo
    assert limit(cos(1/x),x,oo) == 1
    assert limit(exp(x)*(sin(1/x+exp(-x))-sin(1/x)),x,oo) == 1

def test_sign():
    x=Symbol("x")
    from sympy.modules.limits import sign
    assert sign((1/(log(2)+log(x))).eval(),x)==1

def test_others():
    x=Symbol("x")
    a=Symbol("a")
    m=Symbol("m")
    n=Symbol("n")
    assert limit(sqrt(log(x+1))-sqrt(log(x)),x,oo)==0
    assert limit(((1+x)**a-1)/x,x,0)==a
    assert limit((x**(1/n)-1)/(x**(1/m)-1),x,1)==m/n
    #8.12 
    assert limit((3**x+5**x)**(1/x),x,oo)==5
    #this is a similar limit
    assert limit((5**x-3**x)**(1/x),x,oo)==5

def test_Limit():
    x=Symbol("x")
    e=Limit(x*log(x),x,0)
    assert e!=0
    assert e.doit()==0
    e2=Limit(x*log(x*x),x,0)
    assert e!=e2
    assert e.doit()==e2.doit()

def test_error():
    # make sure it exits greacefully if the limit can't be computed
    import py
    x = Symbol('x')
    py.test.raises( NotImplementedError, limit , x , x+1, 0)

def test_atan():
    x = Symbol("x")
    assert limit(atan(x), x, 0) == 0
    assert limit(atan(x), x, -oo) == -pi/2
    assert limit(atan(x), x, oo) == pi/2

def test_minusbug():
    x = Symbol("x")
    a=limit(log(1+exp(x))/x,x,-oo)
    b=limit(log(1+exp(-x))/(-x),x,oo)
    assert a==0
    assert a==b
