#a sandbox, for playing with new ideas and algorithms


import sys
sys.path.append("..")

from sym import exp,log,Symbol,infty,Rational,sin,cos,limit,I,pi

x=Symbol("x") 
y=Symbol("y") 

#print limit((sin(2*x)/x)**(1+x),x,0)
e=(2*x-(7*x**2 - 2) + 3*y)
#print e.print_tree()
#print e.eval().printtree()

#e= (1-(1-x**2))/(x**2*(1+(1-x)**Rational(1,2))) 
#print e
#print limit( e ,x,0)

#print e.expand().subs(x,0)
#print type(Rational(2)**Rational(2))
e=(x+I*y)*(x-I*y)

e=x+I*y
print e.eval()
print e.eval().conjugate()
print e.abs()
