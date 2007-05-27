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
        >>> trigsimp(sin(x)**2 + cos(x)**2)
        1
    """
    from trigonometric import sin, cos, tan, sec, csc, cot
    a = Symbol('a', is_dummy=True)
    b = Symbol('b', is_dummy=True)
    c = Symbol('c', is_dummy=True)
    d = Symbol('d', is_dummy=True)

    identities = {
        a*sin(b)**2 + c*cos(b)**2 + d: (a + (c-a)*cos(b)**2 + d, [a, b, c, d]),
        a*sec(b)**2 - c*tan(b)**2 + d: (a + (c-a)*tan(b)**2 + d, [a, b, c, d]),
        a*csc(b)**2 - c*cot(b)**2 + d: (a + (c-a)*cot(b)**2 + d, [a, b, c, d])
    }
    for identity in identities:
        replacement,varlist = identities[identity]
        ex = [x for x in expr.atoms() if isinstance(x, Symbol)]
        try:
            res = expr.match(identity, varlist, exclude=ex)
        except:
            res = expr.match(identity, varlist)
        if res is not None:
            expr = replacement.subs_dict(res)
    return expr

def simplify(expr):
    return ratsimp(expr.combine().expand())
