#! /usr/bin/python

#This is a testing file, which I use for development. Just ignore it.

import sym as g

def doit():
    x=sym.Symbol("x")
    s=sym.exp(x*x)
    e=s.series(x,14)
    #for a in e.args:
    #    print a,hash(a),a.hash()
    return e

def lim():
    x=g.Symbol("x")
    o=g.Symbol("o")
    #s=(g.log(1/o+x)-x)/g.log(1/o+g.log(x))*(1/o)
    #s=(g.log(1+x)-x-g.log(o))/ (  o*(g.log(1+g.log(x))-g.log(o)) )
    s=(g.exp(x)-1)/x
    s=s.eval()
    print s
    print s.series(x,3)

print "first:"
#lim()

a=g.Symbol("a")
b=g.Symbol("b")
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
