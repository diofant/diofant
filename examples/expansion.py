import sys
sys.path.append("..")

import sym
a=sym.Symbol('a')
b=sym.Symbol('b')
e=(a+b)**5
print e
print e.eval()
print e.expand()
