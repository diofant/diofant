import sys
sys.path.append("..")

from sympy import Derivative, Symbol, Function, exp, Rational, log, \
    dsolve

import relativity


def eq1():
    r = Symbol("r")
    e = relativity.Rmn.dd(0,0)
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    #this should be the same as 2
    #however, that's because we use the relation nu = -lam.
    #so this relation probably follows from eq.1 and 2. But maybe
    #the 3 and 4 needs to be implemented.
    print e
    #print dsolve(e, [relativity.lam(r)])

def eq2():
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

    #the same:
    e = relativity.Rmn.dd(1,1)
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    print e
    print dsolve(e, [relativity.lam(r)])

eq1()
