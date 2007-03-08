import sys
sys.path.append("..")

import sympy
x=sympy.Symbol('x')
e=1/sympy.cos(x)
print e.series(x,10)
e=1/sympy.sin(x)
print e.series(x,4)
