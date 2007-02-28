import sys
sys.path.append("..")

import sym
x=sym.Symbol('x')
e=1/sym.cos(x)
print e.series(x,10)
e=1/sym.sin(x)
print e.series(x,4)
