import sys
sys.path.append("..")

from sympy import exp,log,Symbol,Rational,sin,limit,limitinf

x=Symbol("x")
a=Symbol("a")
h=Symbol("h")

def sqrt(x):
    return x**Rational(1,2)

def sqrt3(x):
    return x**Rational(1,3)

def limitminf(f,x):
    return limitinf(f.subs(x,-x),x)

def show(computed, correct):
    print "computed:",computed,"correct:",correct

show( limitinf(sqrt(x**2-5*x+6)-x,x) , -Rational(5)/2 )
show( limitinf(x*(sqrt(x**2+1)-x),x) , Rational(1)/2 )
show( limitinf(x-sqrt3(x**3-1),x) , Rational(0) )
show( limitminf(log(1+exp(x))/x,x) , Rational(0) )
show( limitinf(log(1+exp(x))/x,x) , Rational(1) )
show( limit(sin(3*x)/x,x,0) , Rational(3) )
show( limit(sin(5*x)/sin(2*x),x,0) , Rational(5)/2 )
show( limitinf(((x-1)/(x+1))**x,x) , exp(-2))
