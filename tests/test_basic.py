import sys
sys.path.append(".")

import sym as g

def dotest(s):
    x=g.symbol("x")
    y=g.symbol("y")
    l=[
    g.rational(2),
    g.real("1.3"), 
    x,
    y,
    pow(x,y)*y,
    5,
    5.5
    ]
    for x in l:
        for y in l:
            s(x,y)

def testbasic():
    def s(a,b):
        x= a
        x= +a
        x= -a
        x= a+b
        x= a-b
        x= a*b
        x= a/b
        x= a**b
    dotest(s)

def testibasic():
    def s(a,b):
        x= a
        x+=b
        x= a
        x-=b
        x= a
        x*=b
        x= a
        x/=b
    dotest(s)

def test_ldegree():
    x=g.symbol("x")
    assert (1/x**2+1+x+x**2).ldegree(x)==-2
    assert (1/x+1+x+x**2).ldegree(x)==-1
    assert (x**2+1/x).ldegree(x)==-1
    assert (1+x**2).ldegree(x)==0
    assert (x+1).ldegree(x)==0
    assert (x+x**2).ldegree(x)==1
    assert (x**2).ldegree(x)==2

def test_leadterm():
    x=g.symbol("x")
    ln=g.ln
    assert (3+2*x**(ln(3)/ln(2)-1)).eval().leadterm(x)==(3,0)
