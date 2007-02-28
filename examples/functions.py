import sys
sys.path.append("..")

import sym
a=sym.Symbol('a')
b=sym.Symbol('b')
e=sym.log((a+b)**5)
print e
print e.eval()
e=sym.exp(e)
print e
print e.eval()
e=sym.log(sym.exp((a+b)**5))
print e
print e.eval()
