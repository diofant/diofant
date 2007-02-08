import sys
sys.path.append("..")

import sym
a=sym.symbol('a')
b=sym.symbol('b')
e=(a+b)**5
print e
print e.eval()
print e.expand()
