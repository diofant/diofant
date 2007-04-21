import sys
sys.path.append("..")

from sympy import Derivative, Symbol, Function, exp, Rational, log, \
    dsolve

import relativity


def eq1():
    r = Symbol("r")
    e = relativity.Rmn.dd(0,0)
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    print e
    print dsolve(e, [relativity.lam(r)])

def eq2():
    #class l(Function): pass
    #class n(Function): pass

    r = Symbol("r")

    #e = Derivative(l(r),r)/r-Derivative(Derivative(n(r),r),r)/2- \
    #    Derivative(n(r),r)**2/4+Derivative(n(r),r)*Derivative(l(r),r)/4

    #e = e.subs(n(r), -l(r))
    #sol = dsolve(e, [l(r)])
    #print e
    #print sol
    #print
    #print (e.subs(l(r), sol).doit()).expand()

    #the same:
    e = relativity.Rmn.dd(1,1)
    e = e.subs(relativity.nu(r), -relativity.lam(r))
    print e
    print dsolve(e, [relativity.lam(r)])

eq1()
#eq2()
