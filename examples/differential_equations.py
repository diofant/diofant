import sys
sys.path.append("..")

from sympy import Derivative, Symbol, Function

class l(Function): pass
class n(Function): pass

r = Symbol("r")

e = Derivative(l(r),r)/r-Derivative(Derivative(n(r),r),r)/2- \
    Derivative(n(r),r)**2/4+Derivative(n(r),r)*Derivative(l(r),r)/4

print e
print e.subs(n(r), -l(r))
