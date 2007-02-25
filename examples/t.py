import sys
sys.path.append("..")

from sym import exp,ln,symbol,infty,rational,sin,cos,limit

x=symbol("x") 
y=symbol("y") 

#print limit((sin(2*x)/x)**(1+x),x,0)
#e=(2*x-(7*x**2 - 2) + 3*y)
#print e.printtree()
#print e.eval().printtree()

e= (1-(1-x**2))/(x**2*(1+(1-x)**rational(1,2))) 
print e
print limit( e ,x,0)

print e.expand().subs(x,0)
