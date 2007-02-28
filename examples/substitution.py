import sys
sys.path.append("..")

import sym
x=sym.Symbol('x')
y=sym.Symbol('y')
e=1/sym.cos(x)
print e
print e.subs(sym.cos(x),y)
print e.subs(sym.cos(x),y).subs(y,x**2)
e=1/sym.log(x)
e=e.subs(x,sym.Real(2.71828))
print e
print e.eval()
print e.evalf()
