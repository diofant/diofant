#a sandbox, for playing with new ideas and algorithms


import sys
sys.path.append("..")
sys.path.append(".")

from sym import exp,log,Symbol,infty,Rational,sin,cos,limit,I,pi,Mul
from sym import hashing, Integral, limitinf

x=Symbol("x") 
y=Symbol("y") 


def sqrt(x):
    return x.sqrt()

#print limit((sin(2*x)/x)**(1+x),x,0)
#e=(2*x-(7*x**2 - 2) + 3*y)
#print e.print_tree()
#print e.printtree()

#e= (1-(1-x**2))/(x**2*(1+(1-x)**Rational(1,2))) 
#print e
#print limit( e ,x,0)

#print e.expand().subs(x,0)
#print type(Rational(2)**Rational(2))
#e=(x+I*y)*(x-I*y)

#e=x+I*y
#print e.conjugate()
#print e.abs()

#p=1-x
#if (1-x).subs(x,1) == 0:
#    print "I found a root 0"
#if (1-x).subs(x,1) == 0.0:
#    print "I found a root 0.0"

e=(x)**2
w=Symbol("w")

import pdb
#pdb.run("limitinf(sqrt(log(x+1))-sqrt(log(x)),x)")
#print limitinf(sqrt(log(x+1))-sqrt(log(x)),x)

#e=sqrt(w-log(w))
#print e
#e2=e.subs(log(w),-x)
#print e2
#e3=e2.series(w,1)
#print e3
#e4=e3.subs(x,-log(w))
#print e4
#print e.series(w,1)

#x=Symbol("x") 
#w=Symbol("w")
#e=(-2)*x.sqrt()-Rational(1)/2*w/x.sqrt()
#print e
#print e.print_tree()

e=(-log(w)).sqrt()
print e
print e.subs(log(w),-x)
