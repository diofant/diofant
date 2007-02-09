import sys
sys.path.append(".")

import sym as s

def testsets():
    from sym.modules.limits import member,intersect,union
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
    from sym.modules.limits import compare
    x=s.symbol("y")
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
    from sym.modules.limits import max
    x=s.symbol("y")
    eq(max([s.exp(x)],[x**5],x),  [s.exp(x)])
    eq(max([s.exp(-x)],[x],x),  [s.exp(-x)])


def testmrv():
    from sym.modules.limits import mrv
    x=s.symbol("y")
    eq(mrv(s.exp(x+1/x),x),[s.exp(x+1/x)])
    eq(mrv(-s.exp(1/x),x),[x])
    eq(mrv(x,x),[x])
    eq(mrv(s.exp(-x),x),[s.exp(-x)])
    eq(mrv(s.exp(x+s.exp(-x)),x),[s.exp(x+s.exp(-x)),s.exp(-x)])
    eq(mrv(s.exp(x+s.exp(-s.exp(x))),x),[s.exp(-s.exp(x))] )
    eq(mrv(s.exp(x+s.exp(-x**2)),x),[s.exp(-x**2)] )

def test_simple_limit_manual():
    "example 3.15"
    from sym.modules.limits import mrv
    x=s.symbol("y")
    f=(s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    Omega=mrv(f,x)
    assert Omega==[s.exp(-x)]
    assert Omega!=[s.exp(x)]
    wexpr=Omega[0]
    w=s.symbol("w")
    f2=f.subs(wexpr,w)
    ser=f2.series(w,3)
    lterm=ser.leadterm(w)
    assert lterm[0]==-s.exp(1/x)
    assert lterm[1]==0
    Omega=mrv(lterm[0],x)
    assert Omega==[x]

def test_simple_limit_lessmanual():
    "example 3.15"
    x=s.symbol("y")
    f=(s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    from sym.modules.limits import mrvleadterm
    lterm=mrvleadterm(f,x)
    assert lterm[0]==-s.exp(1/x)
    assert lterm[1]==0

def test_simple_limit_automatic():
    "example 3.15"
    x=s.symbol("y")
    f=(s.exp(1/x-s.exp(-x))-s.exp(1/x))/s.exp(-x)
    assert s.limitinf(f,x) == -1

def testlimitinf_lenmrveq1():
    x=s.symbol("y")
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
    x=s.symbol("y")
    e=s.exp(s.rational(1))
    assert s.limit((s.exp(x)-1)/x,x,0) == 1
    assert s.limit(s.exp(x),x,0) == 1
    assert s.limit(s.exp(x),x,1) == e
    assert s.limit(s.exp(x),x,-1) == s.exp(s.rational(-1))
    assert s.limit(s.ln(x)*x,x,0) == 0

def testlimitinf_lenmrveq2():
    x=s.symbol("y")
    assert s.limitinf(s.exp(x+s.exp(-x))-s.exp(x),x) == 1
    assert s.limitinf(1/s.exp(-x+s.exp(-x))-s.exp(x),x) == -1
    #example 8.19
    e=(s.ln(s.ln(x)+s.ln(s.ln(x)))-s.ln(s.ln(x)))/s.ln(s.ln(x)+s.ln(s.ln(s.ln(x)))) *s.ln(x)
    assert s.limitinf(e,x)==1

def xtestlonglimit1():
    "example 8.18"
    x=s.symbol("y")
    h=s.exp(-x/(1+s.exp(-x)))
    e=(s.exp(h)*s.exp(-x/(1+h))*s.exp(s.exp(-x+h)))/h**2-s.exp(x)+x
    l=s.limitinf(e,x)
    print "limit=",l
    assert l== 2
    assert l!= 1

def testln():
    x=s.symbol("x")
    e=s.ln(x)
    assert s.limit(e,x,0)==s.infty

def testsubexp():
    x=s.symbol("x")
    e=s.ln(x)
    from sym.modules.limits import subexp
    assert subexp(e,x)
    assert not subexp(e,s.exp(x))
    assert subexp(s.exp(s.exp(x)+s.ln(x)),x)
    assert subexp(s.exp(s.exp(x)+s.ln(x)),s.ln(x))
    assert subexp(s.exp(s.exp(x)+s.ln(x)),s.exp(x))
    assert not subexp(s.exp(s.exp(x)+s.ln(x)),2*x)

def sqrt(x):
    return x**s.rational(1,2)

def test_functions():
    from sym import sin,cos,infty,rational,exp,limit,limitinf
    x=s.symbol("x")
    assert limit(sin(x)/x,x,0) == 1
    assert limit(cos(x)/sin(x),x,0) == infty
    assert limitinf(cos(1/x),x) == 1
    assert limitinf(exp(x)*(sin(1/x+exp(-x))-sin(1/x)),x) == 1

def test_sign():
    x=s.symbol("x")
    from sym.modules.limits import sign
    assert sign((1/(s.ln(2)+s.ln(x))).eval(),x)==1

def test_others():
    x=s.symbol("x")
    a=s.symbol("a")
    m=s.symbol("m")
    n=s.symbol("n")
    ln=s.ln
    assert s.limitinf(sqrt(ln(x+1))-sqrt(ln(x)),x)==0
    assert s.limit(((1+x)**a-1)/x,x,0)==a
    assert s.limit((x**(1/n)-1)/(x**(1/m)-1),x,1)==m/n
    #8.12
    #assert s.limitinf((3**x-5**x)**(1/x),x)==5
