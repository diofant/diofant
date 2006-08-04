#! /usr/bin/python

import sym as g

def doit():
    x=sym.symbol("x")
    s=sym.exp(x*x)
    e=s.series(x,14)
    #for a in e.args:
    #    print a,hash(a),a.hash()
    return e

def lim():
    x=g.symbol("x")
    o=g.symbol("o")
    #s=(g.ln(1/o+x)-x)/g.ln(1/o+g.ln(x))*(1/o)
    #s=(g.ln(1+x)-x-g.ln(o))/ (  o*(g.ln(1+g.ln(x))-g.ln(o)) )
    s=(g.exp(x)-1)/x
    s=s.eval()
    print s
    print s.series(x,3)

print "first:"
#lim()

a=g.symbol("a")
b=g.symbol("b")
#e=(a+b)**5
#print e
#print e.eval()

#e=1/g.cos(a)
#print e.series(a,8)
#e=b*a + -4 + b + a*b + 4 + (a+b)**2
#" becomes "2ab + b + (a+b)^2"
#print e
#print e.eval()
#print (b*a*b).eval()
print g.sin(a).series(a,10)
