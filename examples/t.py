import sys
sys.path.append("..")

from sym import exp,log,Symbol,infty,Rational,sin,cos,limit

x=Symbol("x") 
y=Symbol("y") 

#print limit((sin(2*x)/x)**(1+x),x,0)
e=(2*x-(7*x**2 - 2) + 3*y)
print e.print_tree()
#print e.eval().printtree()

#e= (1-(1-x**2))/(x**2*(1+(1-x)**Rational(1,2))) 
#print e
#print limit( e ,x,0)

#print e.expand().subs(x,0)
