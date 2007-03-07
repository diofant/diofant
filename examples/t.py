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


def sqrt(x):
    return x.sqrt()


print limit((sin(2*x)/x)**(1+x),x,0),2

print limitinf((3**x-5**x)**(1/x),x),5
