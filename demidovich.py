from sym import exp,ln,symbol,infty,rational
from limits import limit,limitinf

x=symbol("x")
def sqrt(x):
    return x**rational(1,2)

def sqrt3(x):
    return x**rational(1,3)

def limitminf(f,x):
    return limitinf(f.subs(x,-x),x)

def testsimpleproblems():
    assert limitinf((x+1)*(x+2)*(x+3)/x**3,x)==1  #172
    #assert limitinf((2**(x+1)+3**(x+1))/(2**x+3**x),x)==3  #175
    assert limitinf(sqrt(x+1)-sqrt(x),x)==0  #179
    assert limitinf((2*x-3)*(3*x+5)*(4*x-6)/(3*x**3+x-1),x)==8  #Primjer 1
    assert limitinf(x/sqrt3(x**3+10),x)==1  #Primjer 2
    assert limitinf((x+1)**2/(x**2+1),x)==1  #181
    assert limitinf(1000*x/(x**2-1),x)==0  #182
    assert limitinf((x**2-5*x+1)/(3*x+7),x)==infty  #183
    assert limitinf((2*x**2-x+3)/(x**3-8*x+5),x)==0  #184
    assert limitminf(ln(1+exp(x))/x,x)==0  #267a
    assert limitinf(ln(1+exp(x))/x,x)==1  #267b
