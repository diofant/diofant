import sys
sys.path.append(".")

import sympy as s

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
    x=s.Symbol("y")
    assert compare(s.exp(x),x**5,x) == ">"
    assert compare(s.exp(x**2),s.exp(x)**2,x) == ">"
    assert compare(s.exp(x),s.exp(x+s.exp(-x)),x) == "="
    assert compare(s.exp(x+s.exp(-x)),s.exp(x),x) == "="
    assert compare(s.exp(x+s.exp(-x)),s.exp(-x),x) == "="
    assert compare(s.exp(s.exp(x)),s.exp(x+s.exp(-s.exp(x))),x) == ">"
    assert compare(s.exp(-x),x,x) ==  ">"
    assert compare(x,s.exp(-x),x) ==  "<"
    assert compare(s.exp(x+1/x),x,x) == ">"
    assert compare(s.exp(s.exp(x)),s.exp(x+s.exp(-s.exp(x))),x) == ">"
    assert compare(s.exp(-s.exp(x)),s.exp(x),x) == ">"
    assert compare(s.exp(s.exp(-s.exp(x))+x),s.exp(-s.exp(x)),x) == "<"

def eq(a,b):
    assert len(a)==len(b)
    for x,y in zip(a,b):
        assert x==y

def testmax():
    from sympy.modules.limits import max
    x=s.Symbol("y")
    eq(max([s.exp(x)],[x**5],x),  [s.exp(x)])
    eq(max([s.exp(-x)],[x],x),  [s.exp(-x)])


def testmrv():
    from sympy.modules.limits import mrv
    x=s.Symbol("y")
    eq(mrv(s.exp(x+1/x),x),[s.exp(x+1/x)])
    eq(mrv(-s.exp(1/x),x),[x])
    eq(mrv(x,x),[x])
    eq(mrv(s.exp(-x),x),[s.exp(-x)])
    eq(mrv(s.exp(x+s.exp(-x)),x),[s.exp(x+s.exp(-x)),s.exp(-x)])
    eq(mrv(s.exp(x+s.exp(-s.exp(x))),x),[s.exp(-s.exp(x))] )
    eq(mrv(s.exp(x+s.exp(-x**2)),x),[s.exp(-x**2)] )

def test_simple_limit_manual():
    "example 3.15"
    from sympy.modules.limits import mrv
    x = s.Symbol("y")
    f = (s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    Omega=mrv(f,x)
    assert Omega==[s.exp(-x)]
    assert Omega!=[s.exp(x)]
    wexpr=Omega[0]
    w = s.Symbol("w")
    f2 = f.subs(wexpr,w)
    ser = f2.series(w,3)
    lterm=ser.leadterm(w)
    assert lterm[0]==-s.exp(1/x)
    assert lterm[1]==0
    Omega=mrv(lterm[0],x)
    assert Omega==[x]

def test_simple_limit_lessmanual():
    "example 3.15"
    x = s.Symbol("y")
    f = (s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    from sympy.modules.limits import mrv_leadterm
    lterm = mrv_leadterm(f,x)
    assert lterm[0]==-s.exp(1/x)
    assert lterm[1]==0

def test_simple_limit_automatic():
    "example 3.15"
    x = s.Symbol("y")
    f = (s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    assert s.limitinf(f,x) == -1

def testlimitinf_lenmrveq1():
    x = s.Symbol("y")
    assert s.limitinf(x,x) == s.infty
    assert s.limitinf(-x,x) == s.infty
    assert s.limitinf(-x**2,x) == s.infty
    assert s.limitinf(1/x,x) == 0
    assert s.limitinf(1/x,x) != 1
    assert s.limitinf(1+1/x,x) == 1
    assert s.limitinf(1+1/x,x) != 0
    assert s.limitinf(s.exp(x),x) == s.infty
    assert s.limitinf(-s.exp(x),x) == s.infty
    assert s.limitinf(-s.exp(1/x),x) == -1
    assert s.limitinf(s.exp(x)/x,x) == s.infty
    assert s.limitinf(s.exp(x)/x,x) != 1
    assert s.limitinf(x+s.exp(-x),x) == s.infty
    assert s.limitinf(x+s.exp(-x**2),x) == s.infty
    assert s.limitinf(x+s.exp(-s.exp(x)),x) == s.infty
    assert s.limitinf(1/x-s.exp(-x),x) == 0
    assert s.limitinf(13+1/x-s.exp(-x),x) == 13

    assert s.limitinf(x+1/x,x) == s.infty

def testlimit():
    x = s.Symbol("y")
    e = s.exp(s.Rational(1))
    assert s.limit((s.exp(x)-1)/x,x,0) == 1
    assert s.limit(s.exp(x),x,0) == 1
    assert s.limit(s.exp(x),x,1) == e
    assert s.limit(s.exp(x),x,-1) == s.exp(s.Rational(-1))
    assert s.limit(s.log(x)*x,x,0) == 0

def testlimitinf_lenmrveq2():
    x = s.Symbol("y")
    assert s.limitinf(s.exp(x+s.exp(-x))-s.exp(x),x) == 1
    assert s.limitinf(1/s.exp(-x+s.exp(-x))-s.exp(x),x) == -1
    #example 8.19
    e=(s.log(s.log(x)+s.log(s.log(x)))-s.log(s.log(x)))/s.log(s.log(x)+s.log(s.log(s.log(x)))) *s.log(x)
    assert s.limitinf(e,x)==1

def xtestlonglimit1():
    "example 8.18"
    x = s.Symbol("y")
    h = s.exp(-x/(1+s.exp(-x)))
    e = (s.exp(h)*s.exp(-x/(1+h))*s.exp(s.exp(-x+h)))/h**2-s.exp(x)+x
    l = s.limitinf(e,x)
    print "limit=",l
    assert l== 2
    assert l!= 1

def testlog():
    x=s.Symbol("x")
    e=s.log(x)
    assert s.limit(e,x,0)==s.infty

def testsubexp():
    x=s.Symbol("x")
    e=s.log(x)
    from sympy.modules.limits import subexp
    assert subexp(e,x)
    assert not subexp(e,s.exp(x))
    assert subexp(s.exp(s.exp(x)+s.log(x)),x)
    assert subexp(s.exp(s.exp(x)+s.log(x)),s.log(x))
    assert subexp(s.exp(s.exp(x)+s.log(x)),s.exp(x))
    assert not subexp(s.exp(s.exp(x)+s.log(x)),2*x)

def sqrt(x):
    return x**s.Rational(1,2)

def test_functions():
    from sympy import sin,cos,infty,Rational,exp,limit,limitinf
    x=s.Symbol("x")
    assert limit(sin(x)/x,x,0) == 1
    assert limit(cos(x)/sin(x),x,0) == infty
    assert limitinf(cos(1/x),x) == 1
    assert limitinf(exp(x)*(sin(1/x+exp(-x))-sin(1/x)),x) == 1

def test_sign():
    x=s.Symbol("x")
    from sympy.modules.limits import sign
    assert sign((1/(s.log(2)+s.log(x))).eval(),x)==1

def test_others():
    x=s.Symbol("x")
    a=s.Symbol("a")
    m=s.Symbol("m")
    n=s.Symbol("n")
    log=s.log
    assert s.limitinf(sqrt(log(x+1))-sqrt(log(x)),x)==0
    assert s.limit(((1+x)**a-1)/x,x,0)==a
    assert s.limit((x**(1/n)-1)/(x**(1/m)-1),x,1)==m/n
    #8.12 
    assert s.limitinf((3**x+5**x)**(1/x),x)==5
    #this is a similar limit
    assert s.limitinf((5**x-3**x)**(1/x),x)==5
