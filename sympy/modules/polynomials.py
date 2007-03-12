"""Module with some routines for polynomials"""

from sympy.core import Pow, Add, Mul, Rational, Number, Symbol

def ispoly(p,x):

    #if not isinstance(x, Symbol):
    #    d=Symbol("d",dummy=True)
    #    e=p.subs(x,d)
    #    if e.has(x): #this "x" is the variable in sin(x), for example
    #        return False
    #    print e,x
    #    return e.ispoly(d)
    if not p.has(x):
        return True
    if isinstance(p,Number):
        return True
    if p==x:
        return True
    if isinstance(p, Pow):
        if isinstance(p.exp, Rational) and p.exp.isinteger():
            if int(p.exp)>0 and ispoly( p.base, x):
                return True
    if isinstance(p,Add):
        a,b = p.getab()
        return ispoly(a, x) and ispoly(b, x)
    if isinstance(p,Mul):
        a,b = p.getab()
        return ispoly(a, x) and ispoly(b, x)
    return False
