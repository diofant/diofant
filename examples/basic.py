import sys
sys.path.append("..")

import sym
a=sym.symbol('a')
b=sym.symbol('b')
c=sym.symbol('c')
e=( a*b*b+2*b*a*b )**c
print e
print e.eval()
