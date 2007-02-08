import sys
sys.path.append("..")

import sym
x=sym.symbol('x')
y=sym.symbol('y')
e=1/sym.cos(x)
print e
print e.subs(sym.cos(x),y)
print e.subs(sym.cos(x),y).subs(y,x**2)
e=1/sym.ln(x)
e=e.subs(x,sym.real(2.71828))
print e
print e.eval()
print e.evalf()
