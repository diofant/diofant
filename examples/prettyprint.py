import sys
sys.path.append("..")

from sympy import exp,log,Symbol,infty,Rational,sin,cos,limit,I,pi,Mul
from sympy.core import basic

x=Symbol("x") 
y=Symbol("y") 

def p():
    print x**x
    print x+y+x
    print sin(x)**x
    print sin(x)**cos(x)
    print sin(x)/(cos(x)**2 * x**x +(2*y))

    print sin(x**2+exp(x))
    print exp(x).sqrt()
    print exp(x).sqrt().sqrt()

print "sympy print:"
p()
basic.outputType="pretty"
print "_"*70
print "pretty print:"
p()
