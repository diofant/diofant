#a sandbox, for playing with new ideas and algorithms


import sys
sys.path.append("..")
sys.path.append(".")

from sym import exp,log,Symbol,infty,Rational,sin,cos,limit,I,pi,Mul
from sym import hashing, Integral, limitinf
from sym.modules import limits
limits.debug=True

x=Symbol("x") 
y=Symbol("y") 
w=Symbol("w") 


def sqrt(x):
    return x.sqrt()


#print limitinf((3**x-5**x)**(1/x),x),5
print limitinf((5**x-3**x)**(1/x),x),5
#e=(w**(-log(5)/log(3))-1/w)**(1/x)
#this later:
#print e
#print e.series(w,1)
#print e.series(w,1).subs(log(w),-log(3)*x)
#print e.series(w,1).subs(log(w),-log(3)*x).subs(w,0),"?=",5

#this first
#e=1/x*log(-w**(-log(5)/log(3))+1/w)
#e=log(w**(-log(5)/log(3))-1/w)
#print e
#print e.series(w,1)
