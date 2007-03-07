import sys
sys.path.append("..")

import sym
a=sym.Symbol('a')
b=sym.Symbol('b')
c=sym.Symbol('c')
e=( a*b*b+2*b*a*b )**c
print e
