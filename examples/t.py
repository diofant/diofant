#a sandbox, for playing with new ideas and algorithms


import sys
sys.path.append("..")
sys.path.append(".")

from sympy import exp,log,Symbol,infty,Rational,sin,cos,limit,I,pi,Mul
from sympy import hashing, Integral, limitinf
from sympy.modules import limits
limits.debug=True

x=Symbol("x") 
y=Symbol("y") 
w=Symbol("w") 


def sqrt(x):
    return x.sqrt()


print limitinf(sin(x)/x,x)
