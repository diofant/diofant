import sys
sys.path.append(".")

import sym as s
import limits

def testsets():
    x=[1,2,3]
    y=[2,4]
    assert limits.member(2,y)
    assert not limits.member(3,y)
    assert limits.intersect(x,y)
    assert not limits.intersect(x,[4,5])
    assert limits.union(x,y)==[1,2,3,4]
    assert limits.union(x,[1,3])==[1,2,3]
    assert limits.union(x,[4,5])==[1,2,3,4,5]

def testcompare():
    x=s.symbol("y")
    assert limits.compare(s.exp(x),x**5,x) == ">"
    assert limits.compare(s.exp(x**2),s.exp(x)**2,x) == ">"
    assert limits.compare(s.exp(x),s.exp(x+s.exp(-x)),x) == "="
    assert limits.compare(s.exp(x+s.exp(-x)),s.exp(x),x) == "="
    assert limits.compare(s.exp(x+s.exp(-x)),s.exp(-x),x) == "="
    assert limits.compare(s.exp(s.exp(x)),s.exp(x+s.exp(-s.exp(x))),x) == ">"
    assert limits.compare(s.exp(-x),x,x) ==  ">"
    assert limits.compare(x,s.exp(-x),x) ==  "<"
    assert limits.compare(s.exp(x+1/x),x,x) == ">"
    assert limits.compare(s.exp(s.exp(x)),s.exp(x+s.exp(-s.exp(x))),x) == ">"
    assert limits.compare(s.exp(-s.exp(x)),s.exp(x),x) == ">"
    assert limits.compare(s.exp(s.exp(-s.exp(x))+x),s.exp(-s.exp(x)),x) == "<"

def eq(a,b):
    assert len(a)==len(b)
    for x,y in zip(a,b):
        assert x==y

def testmax():
    x=s.symbol("y")
    eq(limits.max([s.exp(x)],[x**5],x),  [s.exp(x)])
    eq(limits.max([s.exp(-x)],[x],x),  [s.exp(-x)])


def testmrv():
    x=s.symbol("y")
    eq(limits.mrv(s.exp(x+1/x),x),[s.exp(x+1/x)])
    eq(limits.mrv(-s.exp(1/x),x),[x])
    eq(limits.mrv(x,x),[x])
    eq(limits.mrv(s.exp(-x),x),[s.exp(-x)])
    eq(limits.mrv(s.exp(x+s.exp(-x)),x),[s.exp(x+s.exp(-x)),s.exp(-x)])
    eq(limits.mrv(s.exp(x+s.exp(-s.exp(x))),x),[s.exp(-s.exp(x))] )
    eq(limits.mrv(s.exp(x+s.exp(-x**2)),x),[s.exp(-x**2)] )

def test_simple_limit_manual():
    "example 3.15"
    x=s.symbol("y")
    f=(s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    Omega=limits.mrv(f,x)
    assert Omega==[s.exp(-x)]
    assert Omega!=[s.exp(x)]
    wexpr=Omega[0]
    w=s.symbol("w")
    f2=f.subs(wexpr,w)
    ser=f2.series(w,3)
    lterm=ser.leadterm(w)
    assert lterm[0]==-s.exp(1/x)
    assert lterm[1]==0
    Omega=limits.mrv(lterm[0],x)
    assert Omega==[x]

def test_simple_limit_lessmanual():
    "example 3.15"
    x=s.symbol("y")
    f=(s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    lterm=limits.mrvleadterm(f,x)
    assert lterm[0]==-s.exp(1/x)
    assert lterm[1]==0

def test_simple_limit_automatic():
    "example 3.15"
    x=s.symbol("y")
    f=(s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    assert limits.limitinf(f,x) == -1

def testlimitinf_lenmrveq1():
    x=s.symbol("y")
    assert limits.limitinf(x,x) == s.infty
    assert limits.limitinf(-x,x) == s.infty
    assert limits.limitinf(-x**2,x) == s.infty
    assert limits.limitinf(1/x,x) == 0
    assert limits.limitinf(1/x,x) != 1
    assert limits.limitinf(1+1/x,x) == 1
    assert limits.limitinf(1+1/x,x) != 0
    assert limits.limitinf(s.exp(x),x) == s.infty
    assert limits.limitinf(-s.exp(x),x) == s.infty
    assert limits.limitinf(-s.exp(1/x),x) == -1
    assert limits.limitinf(s.exp(x)/x,x) == s.infty
    assert limits.limitinf(s.exp(x)/x,x) != 1
    assert limits.limitinf(x+s.exp(-x),x) == s.infty
    assert limits.limitinf(x+s.exp(-x**2),x) == s.infty
    assert limits.limitinf(x+s.exp(-s.exp(x)),x) == s.infty
    assert limits.limitinf(1/x-s.exp(-x),x) == 0
    assert limits.limitinf(13+1/x-s.exp(-x),x) == 13

    assert limits.limitinf(x+1/x,x) == s.infty

def testlimit():
    x=s.symbol("y")
    e=s.exp(s.rational(1))
    assert limits.limit((s.exp(x)-1)/x,x,0) == 1
    assert limits.limit(s.exp(x),x,0) == 1
    assert limits.limit(s.exp(x),x,1) == e
    assert limits.limit(s.exp(x),x,-1) == s.exp(s.rational(-1))
    assert limits.limit(s.ln(x)*x,x,0) == 0

def testlimitinf_lenmrveq2():
    x=s.symbol("y")
    assert limits.limitinf(s.exp(x+s.exp(-x))-s.exp(x),x) == 1
    assert limits.limitinf(1/s.exp(-x+s.exp(-x))-s.exp(x),x) == -1
    #example 8.19
    e=(s.ln(s.ln(x)+s.ln(s.ln(x)))-s.ln(s.ln(x)))/s.ln(s.ln(x)+s.ln(s.ln(s.ln(x)))) *s.ln(x)
    assert limits.limitinf(e,x)==1

def xtestlonglimit1():
    "example 8.18"
    x=s.symbol("y")
    h=s.exp(-x/(1+s.exp(-x)))
    e=(s.exp(h)*s.exp(-x/(1+h))*s.exp(s.exp(-x+h)))/h**2-s.exp(x)+x
    l=limits.limitinf(e,x)
    print "limit=",l
    assert l== 2
    assert l!= 1

def testln():
    x=s.symbol("x")
    e=s.ln(x)
    assert limits.limit(e,x,0)==s.infty

def testsubexp():
    x=s.symbol("x")
    e=s.ln(x)
    assert limits.subexp(e,x)
    assert not limits.subexp(e,s.exp(x))
    assert limits.subexp(s.exp(s.exp(x)+s.ln(x)),x)
    assert limits.subexp(s.exp(s.exp(x)+s.ln(x)),s.ln(x))
    assert limits.subexp(s.exp(s.exp(x)+s.ln(x)),s.exp(x))
    assert not limits.subexp(s.exp(s.exp(x)+s.ln(x)),2*x)

def xtesttaylor():
    x=s.symbol("x")
    assert limits.taylor(s.exp,s.rational(0),x,3)==1+x+x**2/2+x**3/6
    assert limits.taylor(s.ln,s.rational(1),x+1,3)==x-x**2/2+x**3/3

def xtest_singularities():
    x=s.symbol("x")
    e=1/(1/x+1)
    assert essential_singularity(e.eval(),limits.series(e.eval(),x,5))==1+x

def xtestseries():
    x=s.symbol("x")
    e=1+x
    assert limits.series(e.eval(),x,5)==1+x
    e=1+x+x**2
    assert limits.series(e.eval(),x,5)==1+x+x**2
    e=s.ln(1+x)/x
    assert limits.series(e.eval(),x,3)==1-x/2+x**2/3
    e=1/(1/x+1)
    print limits.series(e.eval(),x,3)

def sqrt(x):
    return x**s.rational(1,2)

def test_functions():
    from sym import sin,cos,infty,rational,exp
    from limits import limit,limitinf
    x=s.symbol("x")
    assert limit(sin(x)/x,x,0) == 1
    assert limit(cos(x)/sin(x),x,0) == infty
    assert limitinf(cos(1/x),x) == 1
    assert limitinf(exp(x)*(sin(1/x+exp(-x))-sin(1/x)),x) == 1
