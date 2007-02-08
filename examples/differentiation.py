import sys
sys.path.append("..")

import sym
a=sym.symbol('a')
b=sym.symbol('b')
e=(a+2*b)**5
print e
print e.eval()
print e.diff(a)
print e.diff(b)
print e.diff(b).diffn(a,3)
print e.expand().diff(b).diffn(a,3)
