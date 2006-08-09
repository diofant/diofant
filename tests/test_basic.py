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
