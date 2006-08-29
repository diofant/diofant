from sym import exp,ln,symbol,infty,rational
from limits import limit,limitinf

x=symbol("x")
a=symbol("a")
h=symbol("h")

def sqrt(x):
    return x**rational(1,2)

def sqrt3(x):
    return x**rational(1,3)

def sqrt4(x):
    return x**rational(1,4)

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
    assert limitinf((2*x**2-3*x-4)/sqrt(x**4+1),x)==2  #186
    assert limitinf((2*x+3)/(x+sqrt3(x)),x)==2  #187
    assert limitinf(x**2/(10+x*sqrt(x)),x)==infty  #188
    assert limitinf(sqrt3(x**2+1)/(x+1),x)==0  #189
    assert limitinf(sqrt(x)/sqrt(x+sqrt(x+sqrt(x))),x)==1  #190
    #assert limit((x**2-(a+1)*x+a)/(x**3-a**3),x,a)==(a-1)/(3*a**2)  #196
    #assert limit(((x+h)**3-x**3)/h,h,0)==3*x**2  #197
    assert limit((1/(1-x)-3/(1-x**3)),x,1)==-1  #198
    assert limit((sqrt(1+x)-1)/(sqrt3(1+x)-1),x,0)==rational(3)/2  #Primer 4
    assert limit((sqrt(x)-1)/(x-1),x,1)==rational(1)/2  #199
    assert limit((sqrt(x)-8)/(sqrt3(x)-4),x,64)==3  #200
    assert limit((sqrt3(x)-1)/(sqrt4(x)-1),x,1)==rational(4)/3  #201
    assert limit((sqrt3(x**2)-2*sqrt3(x)+1)/(x-1)**2,x,1)==rational(1)/9  #202
    assert limit((sqrt(x)-sqrt(a))/(x-a),x,a)==1/(2*sqrt(a))  #Primer 5
    assert limit((sqrt(x)-1)/(sqrt3(x)-1),x,1)==rational(3)/2  #205
    assert limit((sqrt(1+x)-sqrt(1-x))/x,x,0)==1  #207
    assert limitinf(sqrt(x**2-5*x+6)-x,x)==-rational(5)/2  #213
    assert limitinf(x*(sqrt(x**2+1)-x),x)==rational(1)/2  #214
    assert limitinf(x-sqrt3(x**3-1),x)==0  #215
    assert limitminf(ln(1+exp(x))/x,x)==0  #267a
    assert limitinf(ln(1+exp(x))/x,x)==1  #267b

def xtestbug1():
    assert limitinf((2**(x+1)+3**(x+1))/(2**x+3**x),x)==3  #175
    #the rewritten expression is:
    #(w^(-1)+w^(-ln(3)*ln(2)^(-1)))^(-1)*(3*w^(-ln(3)*ln(2)^(-1))+2*w^(-1))
    #which seems to be correct, and the leading term is 3.
    #but the series() returns 2...
