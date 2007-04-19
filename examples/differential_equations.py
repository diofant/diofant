import sys
sys.path.append("..")

from sympy import Derivative, Symbol, Function, exp, Rational, log

class l(Function): pass
class n(Function): pass

r = Symbol("r")

e = Derivative(l(r),r)/r-Derivative(Derivative(n(r),r),r)/2- \
    Derivative(n(r),r)**2/4+Derivative(n(r),r)*Derivative(l(r),r)/4

e = e.subs(n(r), -l(r))
t = r*exp(-l(r))

t2 = ( t.diffn(r,2)/t ).expand()

print e
a = Symbol("a", is_dummy = True)
tt = (a*t2).expand()
re = e.match(tt, [a])
#there is a bug in match(), it should actually return this:
re = {a: -Rational(1)/2}
assert ( t.diffn(r,2)*re[a]/t ).expand() == e

C1 = Symbol("C1")
C2 = Symbol("C2")
sol = -log(C1+C2/r)
print sol
print
print (e.subs(l(r), sol).doit()).expand()
