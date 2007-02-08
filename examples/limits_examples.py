import sys
sys.path.append("..")

from sym import exp,ln,symbol,rational,sin
from limits import limit,limitinf

x=symbol("x")
a=symbol("a")
h=symbol("h")

def sqrt(x):
    return x**rational(1,2)

def sqrt3(x):
    return x**rational(1,3)

def limitminf(f,x):
    return limitinf(f.subs(x,-x),x)

def show(computed, correct):
    print "computed:",computed,"correct:",correct.eval()

show( limitinf(sqrt(x**2-5*x+6)-x,x) , -rational(5)/2 )
show( limitinf(x*(sqrt(x**2+1)-x),x) , rational(1)/2 )
show( limitinf(x-sqrt3(x**3-1),x) , rational(0) )
show( limitminf(ln(1+exp(x))/x,x) , rational(0) )
show( limitinf(ln(1+exp(x))/x,x) , rational(1) )
show( limit(sin(3*x)/x,x,0) , rational(3) )
show( limit(sin(5*x)/sin(2*x),x,0) , rational(5)/2 )
show( limitinf(((x-1)/(x+1))**x,x) , exp(-2))
