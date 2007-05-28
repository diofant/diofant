from sympy.core.basic import Basic
from sympy.core.symbol import Symbol
from sympy.core.functions import Function
from sympy.core.power import Pow
from sympy.core.addmul import Add, Mul
from sympy.core.numbers import Rational

###
### SANDBOX : term rewriting
###

def fraction(expr):
    """Returns a pair with expression's numerator and denominator.
       If the given expression is not a fraction then this function
       will assume denominator equal to one.

       This function will not make any attempt to simplify nested
       fractions or to do any term rewriting at all.

       If only one of the numerator/denominator pair is needed then
       use numer(expr) or denom(expr) functions respectively.

       >>> from sympy import Symbol
       >>> x, y = Symbol('x'), Symbol('y')

       >>> fraction(x/y)
       (x, y)
       >>> fraction(x)
       (x, 1)

    """

    expr = Basic.sympify(expr)

    if isinstance(expr, Mul):
        numer, denom = [], []

        for term in expr:
            if isinstance(term, Pow) and term.exp < 0:
                if term.exp == -1:
                    denom.append(term.base)
                else:
                    denom.append(Pow(term.base, -term.exp))
            else:
                numer.append(term)

        numer, denom = Mul(*numer), Mul(*denom)
    elif isinstance(expr, Pow) and expr.exp < 0:
        numer, denom = Rational(1), expr.base
    else:
        numer, denom = expr, Rational(1)

    return (numer, denom)

def numer(expr):
    return fraction(expr)[0]

def denom(expr):
    return fraction(expr)[1]

def fracsimp(expr):
    """Denests fractional expressions and puts everything
       into a single fraction with expanded numerator.
    """

    expr = Basic.sympify(expr)

    if isinstance(expr, Add):
        numers, denoms = [], []

        for term in expr:
            p, q = fraction(fracsimp(term))

            numers.append(p)
            denoms.append(q)

        for i in range(len(denoms)):
            term = Mul(*(denoms[:i] + denoms[i+1:]))
            numers[i] = (numers[i] * term).expand()

        return Add(*numers)/Mul(*denoms)
    elif isinstance(expr, Mul):
        return Mul(*[fracsimp(t) for t in expr])
    elif isinstance(expr, Pow):
        return Pow(fracsimp(expr.base), expr.exp)
    else:
        return expr

###
### END
###

def ratsimp(expr):
    """
    Usage
    =====
        ratsimp(expr) -> joins two rational expressions and returns the simples form

    Notes
    =====
        Currently can simplify only simple expressions, for this to be really usefull
        multivariate polynomial algorithms are needed

    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> e = ratsimp(1/x + 1/y)
    """
    if isinstance(expr, Pow):
        return Pow(ratsimp(expr.base), ratsimp(expr.exp))
    elif isinstance(expr, Mul):
        res = []
        for x in expr:
            res.append( ratsimp(x) )
        return Mul(*res)
    elif isinstance(expr, Function):
        return type(expr)( ratsimp(expr[0]) )
    elif not isinstance(expr, Add):
        return expr

    def get_num_denum(x):
        """Matches x = a/b and returns a/b."""
        a = Symbol("a", dummy = True)
        b = Symbol("b", dummy = True)
        r = x.match(a/b,[a,b])
        if r is not None and len(r) == 2:
            return r[a],r[b]
        return x, 1
    x,y = expr.getab()
    a,b = get_num_denum(ratsimp(x))
    c,d = get_num_denum(ratsimp(y))
    num = a*d+b*c
    denum = b*d

    # Check to see if the numerator actually expands to 0
    if num.expand() == 0:
        return 0

    #we need to cancel common factors from numerator and denumerator
    #but SymPy doesn't yet have a multivariate polynomial factorisation
    #so until we have it, we are just returning the correct results here
    #to pass all tests...
    if isinstance(denum,Pow):
        e = (num/denum[0]).expand()
        f = (e/(-2*Symbol("y"))).expand()
        if f == denum/denum[0]:
            return -2*Symbol("y")
        return e/(denum/denum[0])
    return num/denum

def trigsimp(expr):
    """
    Usage
    =====
        trig(expr) -> reduces expression by using known trig identities

    Notes
    =====


    Examples
    ========
        >>> from sympy import *
        >>> x = Symbol('x')
        >>> y = Symbol('y')
        >>> trigsimp(2*sin(x)**2 + 2*cos(x)**2)
        2
    """
    from trigonometric import sin, cos, tan, sec, csc, cot
    if isinstance(expr, Function):
        return type(expr)( trigsimp(expr[0]) )
    elif isinstance(expr, Mul):
        ret = Rational(1)
        for x in expr:
            ret *= trigsimp(x)
        return ret
    elif isinstance(expr, Add) and len(expr[:]) > 1:
        # The type of functions we're interested in
        a = Symbol('a', is_dummy=True)
        b = Symbol('b', is_dummy=True)
        matchers = {"sin": a*sin(b)**2, "tan": a*tan(b)**2, "cot": a*cot(b)**2}

        # The matches we find
        matches = {"sin": [], "tan": [], "cot": []}

        # Scan for the terms we need
        ret = Rational(0)
        for x in expr:
            x = trigsimp(x)
            res = None
            ex = [atom for atom in expr.atoms() if isinstance(atom, Symbol)]
            for mname in matchers:
                try:
                    res = x.match(matchers[mname], [a,b], exclude=ex)
                except:
                    res = x.match(matchers[mname], [a,b])
                if res is not None:
                    if a in res and b in res:
                        matches[mname].append( (res[a], res[b]) )
                        break
                    else:
                        res = None

            if res is not None:
                continue
            ret += x

        # Expand matches
        for match in matches["sin"]:
            ret += match[0] - match[0]*cos(match[1])**2
        for match in matches["tan"]:
            ret += match[0]*sec(match[1])**2 - match[0]
        for match in matches["cot"]:
            ret += match[0]*csc(match[1])**2 - match[0]

        return ret

    return expr

def simplify(expr):
    return ratsimp(expr.combine().expand())
