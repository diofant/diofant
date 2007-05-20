"""Module with some routines for polynomials"""

from sympy.core import Pow, Add, Mul, Rational, Number, Symbol, Basic
from sympy.core.functions import diff

class PolynomialException(Exception):
    pass

def ispoly(p, var=None):
    """
    Usage
    =====
      ispoly(p, var) -> Returns True if p is a polynomial in variable var. 
                        Returns False otherwise.
        
    Notes
    =====
        You can check wether it's a polynomial in several variables at
        once giving a tuple of symbols second argument
        (like ispoly(x**2 + y + 1, (x,y)) ).

    
    Examples
    ========
        >>> from sympy import *
        >>> from sympy.modules.polynomials import *
        >>> x = Symbol('x')
        >>> ispoly(x**2+x+1, x)
        True
        >>> y = Symbol('y')
        >>> ispoly(x**2 + y + 1, (x,y)) #polynomial in variables x and y
        True
        >>> ispoly(x**2 + exp(y) + 1, (x,y))
        False
        
    See also
    ========
       L{coeff_list}, L{coeff}
       
    """
    p = Basic.sympify(p)

    if var == None:
        var = p.atoms(type=Symbol)
    elif isinstance(var, Basic):
        # if the var argument is not a tuple or list
        var = [var] # so we can iterate over it

    if len(var) == 0: 
        return True # constant is polynomial.
    elif len(var) > 1:
        return ispoly(p, var[0]) and ispoly(p, var[1:])

    if not var[0] in p.atoms(type=Symbol): 
        return True # constant is polynomial.

    # Now we look for one variable, that is in the expression.
    if isinstance(p, Pow):
        if isinstance(p.exp, Number) \
           and p.exp.is_integer \
           and p.exp > 0:
            return ispoly(p.base, var[0])
        else:
            return False
    elif isinstance(p, (Add, Mul)):
        a, b = p.getab()
        return ispoly(a, var[0]) and ispoly(b, var[0])
    elif isinstance(p, Number):
        return True
    elif isinstance(p, Symbol):
        return True
    else:
        return False

def fact(n):
    """Returns n!"""
    if n == 0: return 1
    else: return fact(n-1)*n

def coeff(poly, x, n):
    """Returns the coefficient of x**n in the polynomial"""
    assert ispoly(poly,x)
    return diff(poly, x,n).subs(x,0)/fact(n)

def poly(p, var):
    """Returns a sympy polynomial from the representation "p" returned by
    coeff_list().
    """

    if isinstance(var, Symbol):
        var = [var]

    if len(p) == 0:
        return Rational(0)
    elif len(p[0]) != len(var) + 1:
        raise PolynomialException('Wrong number of variables given.')
    
    r = 0
    for item in p:
        c = item[0]
        for v in var:
            c *= v**item[var.index(v)+1]
        r += c
    return r

def rep(n, base):
    """Returns a representation of the integer 'n' in the base 'base'."""
    r = []
    while n!=0:
        r.append(n % base)
        n = (n - r[-1])/base
    return tuple(r)

def gcd(a, b, x):
    """Calculates a greatest common divisor of two polynomials.

    Currently using a heuristics algorithm.
    """

    def getcandidate(a, b, x, x0):
        n1 = a.subs(x, x0)
        n2 = b.subs(x, x0)
        n3 = n1.gcd(int(n1),int(n2))
        c = []
        for n, t in enumerate(rep(n3, x0)):
            if t != 0:
                c.append((t,n))
        return poly(c, x)

    #try some values of x0. If you find polynomials for which gcd doesn't
    #work, just find a number of x0, that works and add it to the end
    #of this list:
    for x0 in [100, 101]:
        c = getcandidate(a, b, x, x0)
        if div(a, c, x)[1] == 0 and div(b, c, x)[1] == 0: 
            return c

    raise PolynomialException("Can't calculate gcd for these polynomials")

def sqf(p, x):
    """Calculates the square free decomposition of 'p'.
    """
    g = gcd(p, p.diff(x), x)
    if g == 1: return p
    a, b = div(p, g, x)
    assert b == 0
    return sqf(a, x) * g

def div(f, g, x):
    """Expresses f = g*q + r, returns (q,r)

    All coefficients of 'f' and 'g' are assumed to be integers,
    and coefficients in 'q' and 'r' are then guaranteed to be integers.
    """
    fp = coeff_list(f, x)
    gp = coeff_list(g, x)
    q = 0
    while fp[0][1] >= gp[0][1] and fp[0][0]!=0:
        s1 = poly([fp[0]], x) / poly([gp[0]], x)
        if isinstance(s1, Mul):
            a,b = s1.getab()
            if isinstance(a, Number) and not a.is_integer:
                #the coefficient is rational but not integer, let's
                #put it in the remainder and we are done
                return q, f
        f = (f - g*s1).expand()
        fp = coeff_list(f, x)
        q+=s1
    return q, f
    
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

def coeff_list(p, var=None, order='lex'):
    """Return the list of coeffs and exponents.

    Currently, lexicographic ('lex') and graded lexicographic
    ('grlex') orders are implemented. The list of variables determines
    the order of the variables.
    The coefficients are assumed to be real numbers, that is, you can
    divide by them.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> from sympy.modules.trigonometric import sin
    >>> coeff_list(x+sin(y)*x**2, [x])
    [[sin(y), 2], [1, 1]]

    >>> coeff_list(6*y**3+7*y*x**2)
    [[7, 2, 1], [6, 0, 3]]

    >>> coeff_list(6*y**3+7*y*x**2, [y,x])
    [[6, 3, 0], [7, 1, 2]]
    
    """
    p = Basic.sympify(p)
    p = p.expand()

    if isinstance(var, Symbol):
        var = [var]
    if var == None:
        var = p.atoms(type=Symbol)
        var.sort()

    res = []
    
    if isinstance(p, Add):
        for a in p._args:
            res.append(*coeff_list(a, var, order))
    else:
        if not isinstance(p, Mul):
            p = Mul(Rational(1), p, evaluate=False) 
        item = [Rational(1)] + len(var)*[Rational(0)]
        for factor in p._args:
            # check if any of the var appear
            if filter(lambda x:x in var, factor.atoms(type=Symbol)):
                if isinstance(factor, Pow) \
                   and (factor.base in var) \
                   and isinstance(factor.exp, Number) \
                   and factor.exp.is_integer \
                   and factor.exp > 0:
                    item[var.index(factor.base)+1] += factor.exp
                elif isinstance(factor, Symbol):
                    item[var.index(factor)+1] += 1
                else:
                    raise PolynomialException('Not a polynomial!')
            else: # the factor is relativly constant
                item[0] *= factor
        res = [item]
    
    # sort list according to monomial ordering
    if order == 'lex':
        res.sort(key=lambda x: x[1:], reverse=True)
    elif order == 'grlex':
        res.sort(key=lambda x: [sum(x[1:])] + x[1:], reverse=True)
    else:
        raise PolynomialException(order + 'is not an implemented order.')

    # unify monomials
    result = []
    for item in res:
        filt = filter(lambda x: x[1:] == item[1:], result)
        if filt:
            result[ result.index(filt[0]) ][0] += item[0]
        else:
            result.append(item)
    
    return result

def div_mv(f, g_i, var=None, order='lex'):
    """Returns q_i and r such that f = g_1*q_1 +...+ g_n*q_n + r

    The g_i can be a single or a list of polynomials. The remainder r
    has a leading term that is not divisible by any of the leading
    terms of the g_i.
    Different monomial orderings are possible, see coeff_list() for
    details.

    Examples:
    >>> x = Symbol('x')
    >>> y = Symbol('y')
    >>> div_mv(x**3+2*x+5, x+1, [x])
    [3+x**2-x, 2]
    >>> div_mv(2*x**3*y**2 - x*y + y**3, x-y, [x,y])
    [2*x**2*y**2+2*y**4-y+2*x*y**3, 2*y**5+y**3-y**2]
    >>> div_mv(2*x**3*y**2 - x*y + y**3, x-y, [y,x])
    [-2*x**3*y-x**2+x-x*y-y**2-2*x**4, -x**2+x**3+2*x**5]
    >>> div_mv(2*x**3*y**2 - x*y + y**3, [x-y, y**2], [x,y])
    [2*x**2*y**2+2*y**4-y+2*x*y**3, -1+y+2*y**3, 0]
    """

    if not isinstance(g_i, list):
        g_i = [g_i]
    f = Basic.sympify(f)
    g_i = map(Basic.sympify, g_i)

    if var == None:
        var = f.atoms(type=Symbol)
        for g in g_i:
            var_g = g.atoms(type=Symbol)
            if var != var_g:
                raise PolynomialException('Bad variables.')

    f_cl = coeff_list(f, var, order)
    g_i_cl = map(lambda g: coeff_list(g, var, order), g_i)

    r = Rational(0)
    q_i = len(g_i)*[Rational(0)]

    while f != 0:
        for g, g_cl in zip(g_i, g_i_cl):
            # check if leading term of f is divisible by that of g
            if all([x>=y for x,y in zip(f_cl[0][1:],g_cl[0][1:])]):
                ff = poly([f_cl[0]], var) / poly([g_cl[0]], var)
                q_i[g_i.index(g)] += ff
                f = (f - ff * g).expand()
                f_cl = coeff_list(f, var, order)
                break
        else: # no division occured, add the leading term to remainder
            ff = poly([f_cl[0]], var)
            r += ff
            f -= ff
            #f_cl = coeff_list(f, var, order)
            f_cl = f_cl[1:]
            

    return q_i + [r]
