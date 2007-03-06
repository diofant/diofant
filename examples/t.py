#a sandbox, for playing with new ideas and algorithms


import sys
sys.path.append("..")

from sym import exp,log,Symbol,infty,Rational,sin,cos,limit,I,pi,Mul
from sym import hashing

x=Symbol("x") 
y=Symbol("y") 

#print limit((sin(2*x)/x)**(1+x),x,0)
#e=(2*x-(7*x**2 - 2) + 3*y)
#print e.print_tree()
#print e.eval().printtree()

#e= (1-(1-x**2))/(x**2*(1+(1-x)**Rational(1,2))) 
#print e
#print limit( e ,x,0)

#print e.expand().subs(x,0)
#print type(Rational(2)**Rational(2))
#e=(x+I*y)*(x-I*y)

#e=x+I*y
#print e.eval()
#print e.eval().conjugate()
#print e.abs()

#p=1-x
#if (1-x).subs(x,1) == 0:
#    print "I found a root 0"
#if (1-x).subs(x,1) == 0.0:
#    print "I found a root 0.0"

a1=cos(x)
a2=(1/sin(x)).eval()

b1=sin(x)
b2=(1/cos(x)).eval()

a= a1*a2
b= b1*b2
print a
print b
print a.eval()
print b.eval()
print a.hash(),b.hash(),a.eval().hash(),b.eval().hash()
print a==b
print a.eval().isequal(b.eval())
print a1.hash(),a2.hash()
print b1.hash(),b2.hash()

mhash=hashing.mhash()
mhash.addstr("<class 'sym.core.numbers.Rational'>")
mhash.addint(-1)
mhash.addint(1)
mhash=mhash.value

xhash=hashing.mhash()
xhash.addstr("<class 'sym.core.symbol.Symbol'>")
xhash.addstr("x")
xhash=xhash.value

a1hash=hashing.mhash()
a1hash.addstr("<class 'sym.modules.trigonometric.cos'>")
a1hash.addint(xhash)
a1hash=a1hash.value

b1hash=hashing.mhash()
b1hash.addstr("<class 'sym.modules.trigonometric.sin'>")
b1hash.addint(xhash)
b1hash=b1hash.value

b2hash=hashing.mhash()
b2hash.addstr("<class 'sym.core.power.Pow'>")
b2hash.add(a1hash)
b2hash.add(mhash)
b2hash=b2hash.value

a2hash=hashing.mhash()
a2hash.addstr("<class 'sym.core.power.Pow'>")
a2hash.add(b1hash)
a2hash.add(mhash)
a2hash=a2hash.value

m=hashing.mhash()
m.addstr("<class 'sym.core.addmul.Mul'>")
m.add(a1hash)
m.add(a2hash)

m2=hashing.mhash()
m2.addstr("<class 'sym.core.addmul.Mul'>")
m2.add(b1hash)
m2.add(b2hash)

print m.value
print m2.value

assert m.value==m2.value
