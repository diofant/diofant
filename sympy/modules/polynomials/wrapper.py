"""Module providing a user-friendly interface to the polynomial algorithms."""

from sympy.core.functions import diff
from sympy.modules.matrices import zero

from sympy.modules.polynomials.base import *
from sympy.modules.polynomials import div_
from sympy.modules.polynomials import gcd_
from sympy.modules.polynomials import groebner_
from sympy.modules.polynomials import lcm_
             
def div(f, g, var=None, order=None, coeff=None):
    """Polynomial division of f by g, returns quotients and remainder.

    Univariate and multivariate polynomials are possible. The argument
    g can also be a list of polynomials to be used as divisors. This
    algorithm doesn't stop when the leading terms don't divide, but
    instead tries to reduce smaller terms after that. When dealing
    with multiple variables, the monomial ordering in use has
    influence over the result. Optionally, a ring of coefficients can
    be indicated, to restrict computations to this domain.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> div(x**2+6*x+1, 3*x-1)
    (19/9+1/3*x, 28/9)
    >>> div(x**2+6*x+1, 3*x-1, coeff='int')
    (2, 3+x**2)
    >>> div(2*x**3*y**2 - x*y + y**3, [x-y, y**2], [x,y], 'lex')
    ([2*x**2*y**2+2*y**4-y+2*x*y**3, -1+y+2*y**3], 0)

    """
    f = Basic.sympify(f)
    if not isinstance(g, list):
        g = [g]
    g = map(lambda x: Basic.sympify(x), g)
    if not isinstance(var, list):
        var = [var]
    if len(var) > 0 and var[0] == None:
        var = merge_var(f.atoms(type=Symbol,
                        *[g_i.atoms(type=Symbol) for g_i in g]))
    f = Polynomial(f, var, order, coeff)
    g = map(lambda x: Polynomial(x, var, order, coeff), g)
    if coeff != None:
        if not coeff in coeff_rings:
            raise PolynomialException(
                "%s is not an implemented coefficient ring." % coeff)
        elif coeff == 'int':
            # TODO: Check if given polynomials have integer coeffs?
            q, r = div_.mv_int(f, g)
    else: # coeff == None
        q, r = div_.mv(f, g)
    if len(q) == 1:
        return q[0].basic, r.basic
    else:
        return map(lambda x: x.basic, q), r.basic

def gcd(f, g, var=None, order=None, coeff=None):
    """Greatest common divisor of two polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> gcd(4*x**2*y, 6*x*y**2)
    x*y
    >>> gcd(4*x**2*y, 6*x*y**2, coeff='int')
    2*x*y

    """
    f = Basic.sympify(f)
    g = Basic.sympify(g)    
    if isinstance(var, Symbol):
        var = [var]
    if var == None:
        var = merge_var(f.atoms(type=Symbol), g.atoms(type=Symbol))
    if len(var) == 0:
        if coeff == 'int':
            assert isinstance(f, Rational) and isinstance(g, Rational)
            assert f.is_integer and g.is_integer
            return abs(Rational(0).gcd(f.p, g.p))
        else:
            return Rational(1)
    elif len(var) == 1:
        if coeff == 'int':
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            return abs(Rational(0).gcd(cf.p, cg.p)) * \
                   gcd_.uv(f, g).basic
        else:
            r =  gcd_.uv(Polynomial(f, var, order, coeff),
                         Polynomial(g, var, order, coeff))
            return r.basic
    else:
        if coeff == 'int':
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            return abs(Rational(0).gcd(cf.p, cg.p)) * \
                   gcd_.mv(f, g).basic
        else:
            r =  gcd_.mv(Polynomial(f, var, order, coeff),
                         Polynomial(g, var, order, coeff))
            return r.basic

def groebner(f, var=None, order=None, coeff=None, reduced=True):
    """Computes the (reduced) Groebner base of given polynomials.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> groebner([x**2 + y**3, y**2-x], order='lex')
    [x-y**2, y**3+y**4]

    """
    if isinstance(f, Basic):
        f = [f]
    f = map(lambda p: Basic.sympify(p).expand(), f)
    # filter trivial or double entries
    ff = filter(lambda x: x!=0, f)
    if not ff: # Zero Ideal
        return [Rational(0)]
    f = []
    for p in ff:
        if not p in f:
            f.append(p)
    if isinstance(var, Symbol):
        var = [var]
    if var == None:
        var = merge_var(*map(lambda p: p.atoms(type=Symbol),f))
    f = map(lambda p: Polynomial(p, var, order, coeff), f)
    g = groebner_.groebner(f, reduced)
    return map(lambda p: p.basic, g)

def lcm(f, g, var=None, order=None, coeff=None):
    """Least common divisor of two polynomials.
    
    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> lcm(4*x**2*y, 6*x*y**2)
    x**2*y**2
    >>> lcm(4*x**2*y, 6*x*y**2, coeff='int')
    12*x**2*y**2
    """
    f = Basic.sympify(f)
    g = Basic.sympify(g)
    if isinstance(var, Symbol):
        var = [var]
    if var == None:
        var = merge_var(f.atoms(type=Symbol), g.atoms(type=Symbol))
    if len(var) == 0:
        if coeff == 'int':
            assert isinstance(f, Rational) and isinstance(g, Rational)
            assert f.is_integer and g.is_integer
            if f == Rational(0) or g == Rational(0):
                return Rational(0)
            return abs(f*g / Rational(0).gcd(f.p, g.p))
        else:
            return Rational(1)
    if len(var) == 1:
        if coeff == 'int':
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            if cf == Rational(0) or cg == Rational(0):
                return Rational(0)
            return abs(cf*cg / Rational(0).gcd(cf.p, cg.p)) * \
                   lcm_.uv(f, g).basic
        else:
            r = lcm_.uv(Polynomial(f, var, order, coeff),
                        Polynomial(g, var, order, coeff))
            return r.basic
    else:
        if coeff == 'int':
            f = Polynomial(f, var, order, coeff)
            cf = f.content()
            f.cl = map(lambda t:[t[0]/cf] + t[1:], f.cl)
            g = Polynomial(g, var, order, coeff)
            cg = g.content()
            g.cl = map(lambda t:[t[0]/cg] + t[1:], g.cl)
            if cf == Rational(0) or cg == Rational(0):
                return Rational(0)
            return abs(cf*cg / Rational(0).gcd(cf.p, cg.p)) * \
                   lcm_.mv(f, g).basic
        else:
            r = lcm_.mv(Polynomial(f, var, order, coeff),
                        Polynomial(g, var, order, coeff))
            return r.basic
            
class Ideal:
    """Describes a polynomial ideal over a field, with several variables.

    Try to avoid different variables and orders between the ideals, give
    them explicitly; the automatic handlers don't always act as expected.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> I = Ideal([x, y**2])
    >>> J = Ideal(x*y)
    >>> I == J
    False
    >>> I*J == I.intersect(J)
    False
    >>> I**2 == I*I
    True
    """
    def __init__(self, f=[Rational(0)], var=None, order='grevlex',
                 is_groebner=None):
        if not isinstance(f, list):
            f = [f]
        self.f = map(Basic.sympify, f)
        if var == None:
            var = []
            for p in self.f:
                var += filter(lambda x: not x in var, p.atoms(type=Symbol))
            var.sort(key=str)
        self.var = var
        self.order = order
        self.is_groebner = is_groebner
        for p in self.f:
            if not ispoly(p, self.var):
                raise PolynomialException('Non-polynomial as generator.')

    def __len__(self):
        return len(self.f)

    def __iter__(self):
        return self.f.__iter__()

    def __str__(self):
        return self.f.__str__()

    def __repr__(self):
        return 'Ideal(%s, %s, %s, %s)' % (repr(self.f), repr(self.var), 
               self.order, repr(self.is_groebner) )

    def __add__(self, other):
        """f is in I + J iff there are g in I, h in J such that f = g+h.
        """
        if not isinstance(other, Ideal):
            other = Ideal(other)
    
        var = self.var + filter(lambda x: not x in self.var, other.var)
        var.sort(key=str)

        if self.order == other.order:
            order = self.order
        else:
            order = None

        f = self.f + filter(lambda x: not x in self.f, other.f)
        
        return Ideal(f, var, order)

    def __radd__(self, other):
        return self.__add__(other)

    def __mul__(self, other):
        """f is in I*J iff f is a sum of products of elements of I and J.
        """
        if not isinstance(other, Ideal):
            other = Ideal(other)

        var = self.var + filter(lambda x: not x in self.var, other.var)
        var.sort(key=str)

        if self.order == other.order:
            order = self.order
        else:
            order = None

        f = []
        for p in self.f:
            for q in other.f:
                f.append(p*q)

        return Ideal(f, var, order)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        """Repeated product of self.
        """
        other = Basic.sympify(other)
        if not (isinstance(other, Number)
                and other.is_integer
                and other >= 0):
            raise PolynomialException('Illegal power of ideal.')

        if other == 0:
            return Ideal(Rational(1), self.var, self.order, True)
        elif other == 1:
            return self

        # avoid repeated elements
        f = self.f[:]
        I = Ideal(Rational(0), self.var, self.order)

        for p in self.f:
            I += p*Ideal(f,self.var, self.order)**(other-1)
            f.remove(p)

        return I

    def __mod__(self, other):
        if not isinstance(other, Ideal):
            other = Ideal(other)
        other.groebner()

        f = []
        for p in self.f:
            f.append(div(p, other.f, other.var, other.order)[-1])
        f = filter(lambda x: x!=0, f)
        return Ideal(f, other.var, other.order)

    def __contains__(self, other):
        other = Basic.sympify(other)
        if not ispoly(other, self.var):
            return False
        self.groebner()
        rem = div(other, self.f, self.var, self.order)[-1]
        if rem == 0:
            return True
        else:
            return False

    def __eq__(self, other):
        if not isinstance(other, Ideal):
            other = Ideal(other, self.var, self.order)
        # Try to save Groebner base computations.       
        if self.var != other.var or self.order != other.order:
            if self.is_groebner:
                s = self
                o = Ideal(other.f, self.var, self.order)
                o.groebner()
            elif other.is_groebner:
                o = other
                s = Ideal(self.f, other.var, other.order)
                s.groebner()
            else:
                s = Ideal(self.f)
                o = Ideal(other.f)
                s.groebner()
                o.groebner()
            return s.f == o.f
        else:
            self.groebner()
            other.groebner()
            return self.f == other.f

    def groebner(self, reduced=True):
        if not self.is_groebner:
            self.f = groebner(self.f, self.var, self.order, None, reduced)
            self.is_groebner = True

    def intersect(self, other):
        if not isinstance(other, Ideal):
            other = Ideal(other)
        t = Symbol('t', dummy=True)
        I = t*self + (1-t)*other
        I.order = '1-el'
        I.groebner()
        f = filter(lambda p: not t in p.atoms(type=Symbol), I.f)
        if self.var == other.var:
            var = self.var
        else:
            var = None
        if self.order == other.order:
            order = self.order
        else:
            order = None
        return Ideal(f, var, order)

def coeff(poly, x, n):
    """Returns the coefficient of x**n in the polynomial"""
    assert ispoly(poly,x)
    p = Polynomial(poly, x)
    t = filter(lambda t:t[1]==n,p.cl)
    if len(t) == 0:
        return Rational(0)
    else:
        return t[0][0]
#    return diff(poly, x,n).subs(x,0)/fact(n)

def sqf(p, x):
    """Calculates the square free decomposition of 'p'.
    """
    g = gcd(p, p.diff(x), x)
    if g == 1: return p
    a, b = div(p, g, x)
    assert b == 0
    return sqf(a, x) * g

def resultant(f, g, x, method='bezout'):
    """Computes resultant of two polynomials in one variable. This
       method can be used to verify if selected polynomials have
       common root withot factoring them or computing any GCD's.

       It can be also be used for variable elemination in case of
       bivariate polynomials. Just assume that one of the variables
       is a parameter and compute resultant with respect to the other
       one and you will get univariate polynomial in single variable.

       For now two algorithm have been implemented, based on
       Sylvester and Bezout matrices. The default is Bezout.

       TODO: Make Bezout approach run in O(s**2). Currently
             it is O(s**3), where s = max(deg(f), deg(g)).
    """

    fp = coeff_list(f, x)
    gp = coeff_list(g, x)

    m, n = int(fp[0][1]), int(gp[0][1])

    if method is 'sylvester':
        M = zero(m+n)

        for i in range(n):
            for coeff, j in fp:
                M[i, m-int(j)+i] = coeff

        for i in range(m):
            for coeff, j in gp:
                M[i+n, n-int(j)+i] = coeff

        return M.det()
    elif method is 'bezout':
        if n > m:
            s, fp, gp = n, gp, fp
        else:
            s = m

        p, q = [0]*(s+1), [0]*(s+1)

        for coeff, j in fp:
            p[int(j)] = coeff

        for coeff, j in gp:
            q[int(j)] = coeff

        M = zero(s)

        for i in range(s):
            for j in range(i, s):
                z = 1 + min(i, s-1-j)
                terms = [0] * z

                for k in range(z):
                    terms[k] = p[j+k+1]*q[i-k] - p[i-k]*q[j+k+1]

                M[i, j] = M[j, i] = Add(*terms)

        if (s*(s-1)/2) % 2 == 0:
            b = p[-1]**abs(n-m)
        else:
            b = -p[-1]**abs(n-m)

        return (1/b) * M.det()
    else:
        raise ValueError("Invalid method: '%s'" % method)

def collect(expr, syms):
    """Collect additive terms with respect to a list of variables in a linear
       multivariate polynomial. This function assumes the input expression is
       in an expanded form and will return None if this is not a linear
       polynomial or else a pair of the following form:

          ( { variable : coefficient }, free term )

       Example:
       >>> from sympy import *
       >>> from sympy.modules.polynomials import collect
       >>> x, y, z = Symbol('x'), Symbol('y'), Symbol('z')
       >>> collect(2*x + sin(z)*x + cos(z)*y + 1, [x, y])
       ({x: 2+sin(z), y: cos(z)}, 1)

    """
    if isinstance(expr, (Add, Mul)):
        content, tail = {}, 0

        if isinstance(expr, Mul):
            expr = [expr]

        for term in expr:
            coeff = 1

            while isinstance(term, Mul):
                a, term = term.getab()

                if isinstance(a, Symbol) and a in syms:
                    if (term.has_any(syms)):
                        return None
                    else:
                        coeff *= term

                        if a in content:
                            content[a] += coeff
                        else:
                            content[a] = coeff

                        break
                else:
                    if (a.has_any(syms)):
                        return None
                    else:
                        coeff *= a
            else:
                if isinstance(term, Symbol) and term in syms:
                    if term in content:
                        content[term] += coeff
                    else:
                        content[term] = coeff
                else:
                    tail += coeff * term

        return (content, tail)
    elif isinstance(expr, Symbol) and expr in syms:
        return ({expr : 1}, 0)
    elif isinstance(expr, Basic) and not expr.has_any(syms):
        return ({}, expr)
    else:
        return None

def roots(a, var=None):
    """Returns all rational roots of the equation a=0, where a is a
    polynomial."""

    def find_divisors(num):
        num = abs(int(num))
        r = []
        for x in range(1,num+1):
            if num % x == 0:
                r.append(x)
        return r

    def candidates(f, l):
        r = []
        for x in f:
            for y in l:
                r.append(Rational(y,x))
                r.append(-Rational(y,x))
        return r

    if var==None:
        var = a.atoms(type=Symbol)[0]

    r = []
    c = coeff_list(a)
    if len(c[0]) == 1:
        return []
    if c[-1][1] != 0:
        r.append(Rational(0))
    lastnum = c[-1][0]
    firstnum = c[0][0]
    if not lastnum.is_integer:
        m = lastnum.q
        lastnum *= m
        firstnum *= m
    if not firstnum.is_integer:
        m = firstnum.q
        lastnum *= m
        firstnum *= m
    last_d = find_divisors(lastnum)
    first_d = find_divisors(firstnum)
    for t in candidates(first_d, last_d):
        if a.subs(var, t) == 0:
            r.append(t)
    return r

def factor(a, var=None):
    """Factors the polynomial a.

    Example:
    >>> x = Symbol('x')
    >>> factor(x**6-1)
    (1+x)*(1+x**4+x**2)*(-1+x)

    Note: as you can see, it only factors out the rational roots, here the
    correct answer should be:
    (x - 1)*(x + 1)*(x**2 - x + 1)*(x**2 + x + 1)
    """
    if var==None:
        var = a.atoms(type=Symbol)[0]
    p = 1
    for r in roots(a, var): 
        p *= (var-r)
    if p == 1: 
        return a
    q,r = div(a, p, var)
    assert r == 0
    return factor(q, var)*p
