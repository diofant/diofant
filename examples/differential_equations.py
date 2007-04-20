import sys
sys.path.append("..")

from sympy import Derivative, Symbol, Function, exp, Rational, log, \
    dsolve

class l(Function): pass
class n(Function): pass

r = Symbol("r")

e = Derivative(l(r),r)/r-Derivative(Derivative(n(r),r),r)/2- \
    Derivative(n(r),r)**2/4+Derivative(n(r),r)*Derivative(l(r),r)/4

e = e.subs(n(r), -l(r))
sol = dsolve(e, [l(r)])
print e
print sol
print
print (e.subs(l(r), sol).doit()).expand()
