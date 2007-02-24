import sys
sys.path.append("..")

from sym import exp,ln,symbol,infty,rational,sin,cos,limit

x=symbol("x") 
y=symbol("y") 

#print limit((sin(2*x)/x)**(1+x),x,0)
e=(2*x-(7*x**2 - 2) + 3*y)
print e.printtree()
print e.eval().printtree()
