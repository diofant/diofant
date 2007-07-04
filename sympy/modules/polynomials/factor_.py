"""Various algorithms for the factorization of polynomials."""

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import div_
from sympy.modules.polynomials import sqf_

def uv_int(f):
    """Find the factorization of an univariate integer polynomial.

    Using square-free factorization and Kronecker's algorithm."""
    a = sqf_.uv_int(f)
    result = []
    for i, p in enumerate(a):
        # Filter out constant Rational(1) factors
        if p.cl[0][1] != 0:
            # TODO: Check for rational roots first? (cheaper)
            for pp in kronecker(p):
                result += [pp]*(i+1)
    return result

def kronecker(f):
    """Recursive factorization of an univariate polynomial with integer
    coefficients using interpolation.
    """
    def lagrange_base(pos):
        """Compute the base polynomials used for interpolation."""
        l=[]
        for x in pos:
            l.append(Polynomial([[Rational(1),1],[-x,0]], f.var, f.order, f.coeff))
        b=[]
        for i, x in enumerate(pos):
            p = Polynomial([[Rational(1),0]], f.var, f.order, f.coeff)
            for ll in l[:i]+l[i+1:]:
                p *= ll
            c = Rational(1)
            for xx in pos[:i]+pos[i+1:]:
                c *= (x-xx)
            p.cl = map(lambda t:[t[0]/c]+t[1:], p.cl)
            b.append(p)
        return b

    def divisors(n):
        n = abs(n)
        r = []
        for i in range(1, n/2+1):
            if n % i == 0:
                r.append(i)
        r.append(n)
        return r

    def combine(divs):
        # Don't try negative divisors for first value.
        lst = map(lambda el: [el], divs[0])
        for choices in divs[1:]:
            lst2 = []
            for el in lst:
                for new in choices:
                    # Also try negative divisors:
                    lst2.append(el + [new])
                    lst2.append(el + [-new])
            lst = lst2
        return lst

    # Half the degree for a possible polynomial divisor g.
    deg = int(f.cl[0][1])/2
    # Points for interpolation
    pos = map(Rational, range(0-deg/2, deg+1-deg/2))
    # Reusable Lagrange base polynomials.
    base = lagrange_base(pos)
    # Evaluate at points.
    values = map(f, pos)
    # Look for roots, that is, zeros in values:
    # TODO: Move out of kronecker()!
    lfs = []
    for x, y in zip(pos, values):
        if y == 0:
            lfs.append(Polynomial(
                [[Rational(1),1], [-x, 0]], f.var, f.order, f.coeff))
    if len(lfs) > 0:
        ff = Polynomial([[Rational(1),0]], f.var, f.order, f.coeff)
        for lf in lfs:
            ff *= lf
        return lfs + kronecker(div_.mv_int(f, ff)[0][0])
    # All divisors of the values give possible values for g.
    divs = map(divisors, map(int, values))
    # Assemble all possible divisor combination
    combs = combine(divs)
    # Construct candidates for g.
    cands = []
    for comb in combs:
        cand = Polynomial(Rational(0), f.var, f.order, f.coeff) 
        for c,b in zip(comb, base):
            cand += c*b
        # Filter out constant and non-integer polynomials!
        if not (len(cand.cl) == 1 and cand.cl[0][1] == 0):
            if all(map(lambda t:t[0].is_integer, cand.cl)):
                cands.append(cand)

    # Get leading coefficient positive:
    for cand in cands:
        if cand.cl[0][0] < 0:
            cand.cl = map(lambda t:[t[0]*Rational(-1)] + t[1:], cand.cl)

    # Filter double entries:
    cands2 = []
    for cand in cands:
        if not cand in cands2:
            cands2.append(cand)
    cands = cands2

    # TODO: Too many candidates?
    # TODO: Use iterators instead of lists!

    for g in cands:
        q, r =  div_.mv_int(f,g)
        if r.cl[0][0] == 0:
            return kronecker(q[0]) + kronecker(g)
    else:
        # No divisor found, f irreducible.
        # TODO: Try again with smaller degree divisors?
        return [f]
