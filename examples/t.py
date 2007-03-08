#a sandbox, for playing with new ideas and algorithms


import sys
sys.path.append("..")
sys.path.append(".")

from sympy import exp,log,Symbol,infty,Rational,sin,cos,limit,I,pi,Mul
from sympy import hashing, Integral, limitinf
from sympy.modules import limits
from sympy.core import basic
limits.debug=True
basic.outputType="pretty"

x=Symbol("x") 
y=Symbol("y") 
w=Symbol("w") 


def sqrt(x):
    return x.sqrt()

#import pdb
#pdb.run('print limitinf(sin(x)/x,x)')
print x**x
print x+y+x
print sin(x)**x
print sin(x)**cos(x)
print sin(x)/(cos(x)**2 * x**x +(2*y))

print sin(x**2+exp(x))
print exp(x).sqrt()
print exp(x).sqrt().sqrt()
