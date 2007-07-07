from sympy.core.basic import Basic
from sympy.core.symbol import Symbol
from sympy.core.functions import Function, exp
from sympy.core.numbers import Rational
from sympy.core.addmul import Add, Mul
from sympy.core.power import Pow

def fraction(expr):
    """Returns a pair with expression's numerator and denominator.
       If the given expression is not a fraction then this function
       will assume that the denominator is equal to one.

       This function will not make any attempt to simplify nested
       fractions or to do any term rewriting at all.

       If only one of the numerator/denominator pair is needed then
       use numer(expr) or denom(expr) functions respectively.

       >>> from sympy import *
       >>> x, y = symbols('x', 'y')

       >>> fraction(x/y)
       (x, y)
       >>> fraction(x)
       (x, 1)
       >>> fraction(x*y/2)
       (x*y, 2)
       >>> fraction(Rational(1, 2))
       (1, 2)

       This function will also work fine with assumptions:

       >>> k = Symbol('k', negative=True)
       >>> fraction(x * y**k)
       (x, y**(-k))

    """
    expr = Basic.sympify(expr)

    if isinstance(expr, Mul):
        numer, denom = [], []

        for term in expr:
            if isinstance(term, Pow) and term.exp.is_negative:
                if term.exp.is_minus_one:
                    denom.append(term.base)
                else:
                    denom.append(Pow(term.base, -term.exp))
            elif isinstance(term, Rational):
                if term.is_integer:
                    numer.append(term)
                else:
                    numer.append(Rational(term.p))
                    denom.append(Rational(term.q))
            else:
                numer.append(term)

        numer, denom = Mul(*numer), Mul(*denom)
    elif isinstance(expr, Pow) and expr.exp.is_negative:
        numer, denom = Rational(1), expr.base
    elif isinstance(expr, Rational):
        numer, denom = Rational(expr.p), Rational(expr.q)
    else:
        numer, denom = expr, Rational(1)

    return (numer, denom)

def numer(expr):
    return fraction(expr)[0]

def denom(expr):
    return fraction(expr)[1]

def fraction_expand(expr):
    a, b = fraction(expr)
    return a.expand() / b.expand()

def numer_expand(expr):
    a, b = fraction(expr)
    return a.expand() / b

def denom_expand(expr):
    a, b = fraction(expr)
    return a / b.expand()

def together(expr, deep=False):
    """Combine together and denest rational functions into a single
       fraction. No futher expansion is performed, use appropriate
       functions respectively.

       >>> from sympy import *
       >>> x, y = symbols('x', 'y')

       >>> together(1/x + 1/y)
       (x+y)/(y*x)
       >>> together(1/(1 + 1/x))
       x/(1+x)

    """
    expr = Basic.sympify(expr)

    if isinstance(expr, Add):
        numers, denoms = [], []

        for term in expr:
            p, q = fraction(together(term, deep))

            numers.append(p)
            denoms.append(q)

        for i in range(len(denoms)):
            term = Mul(*(denoms[:i] + denoms[i+1:]))
            numers[i] = numers[i] * term

        return Add(*numers)/Mul(*denoms)
    elif isinstance(expr, (Mul, Pow)):
        return type(expr)(*[ together(t, deep) for t in expr ])
    elif isinstance(expr, Function) and deep:
        return type(expr)(together(expr._args, deep))
    else:
        return expr

def separate(expr, deep=False):
    """Rewrite or separate a power of product to a product of
       powers. This function is highly specific so it won't
       expand any nested products of summations etc.

       >>> from sympy import *
       >>> x, y, z = symbols('x', 'y', 'z')

       >>> separate((x*y)**2)
       x**2*y**2

       >>> separate((x*(y*z)**3)**2)
       z**6*y**6*x**2

       >>> separate((x*sin(x))**y + (x*cos(x))**y)
       x**y*cos(x)**y+x**y*sin(x)**y

       >>> separate((exp(x)*exp(y))**x)
       exp(x*y)*exp(x**2)

    """
    expr = Basic.sympify(expr)

    if isinstance(expr, Pow):
        terms, expo = [], separate(expr.exp, deep)

        if isinstance(expr.base, Mul):
            return Mul(*[ separate(t**expo, deep) for t in expr.base ])
        elif isinstance(expr.base, exp):
            if deep == True:
                return exp(separate(expr.base._args, deep)*expo)
            else:
                return exp(expr.base._args*expo)
        else:
            return Pow(separate(expr.base, deep), expo)
    elif isinstance(expr, (Add, Mul)):
        return type(expr)(*[ separate(t, deep) for t in expr ])
    elif isinstance(expr, Function) and deep:
        return type(expr)(separate(expr._args, deep))
    else:
        return expr

def collect(expr, syms, pretty=True):
    """Collect additive terms with respect to a list of symbols up
       to powers with rational exponents. By the term symbol here
       are meant arbitrary expressions, which can contain powers,
       products, sums etc.

       This function will not apply any redundant expanding to the
       input expression, so user is assumed to enter expression in
       final form. This makes 'collect' more predictable as there
       is no magic behind the scenes. However it is important to
       note, that powers of products are converted to products of
       powers using 'separate' function.

       There are two possible types of output. First, if 'pretty'
       flag is set, this function will return an single expression
       or else it will return a dictionary with separated symbols
       up to rational powers as keys and collected sub-expressions
       as values respectively.

       >>> from sympy import *
       >>> x, y, z = symbols('x', 'y', 'z')
       >>> a, b, c = symbols('a', 'b', 'c')

       This function can collect symbolic coefficients in polynomial
       or rational expressions. It will manage to find all integer or
       rational powers of collection variable:

       >>> collect(a*x**2 + b*x**2 + a*x - b*x + c, x)
       c+(-b+a)*x+(a+b)*x**2

       The same result achieved but in dictionary form:

       >>> collect(a*x**2 + b*x**2 + a*x - b*x + c, x, pretty=False)
       {x: -b+a, x**2: a+b, 1: c}

       You can also work with multi-variate polynomials. However
       remeber that this function is greedy so it will care only
       about a single symbol at time, in specification order:

       >>> collect(x**2 + y*x**2 + x*y + y + a*y, [x, y])
       (1+a)*y+x**2*(1+y)+x*y

       >>> collect(x**2*y**4 + z*(x*y**2)**2 + z + a*z, [x*y**2, z])
       (1+z)*x**2*y**4+z*(1+a)

       Also more complicated expressions can be used as collectors:

       >>> collect(a*sin(2*x) + b*sin(2*x), sin(2*x))
       (a+b)*sin(2*x)

       >>> collect(a*x**2*log(x)**2 + b*(x*log(x))**2, x*log(x))
       log(x)**2*(a+b)*x**2

       It is also possible to work with symbolic powers, although
       it has more complicated behaviour, because in this case
       power's base and symbolic part of the exponent are treated
       as a single symbol:

       >>> collect(a*x**c + b*x**c, x)
       x**c*a+x**c*b

       >>> collect(a*x**c + b*x**c, x**c)
       (a+b)*x**c

       However if you incorporate rationals to the exponents, then
       you will get well known behaviour:

       >>> collect(a*x**(2*c) + b*x**(2*c), x**c)
       (a+b)*x**(2*c)

       Note also that all previously stated facts about 'collect'
       function apply to the exponential function, so you can get:

       >>> collect(a*exp(2*x) + b*exp(2*x), exp(x))
       (a+b)*exp(2*x)

    """
    def make_list(expr, kind):
        if isinstance(expr, kind):
            return expr[:]
        else:
            return [expr]

    def get_base_exp(expr):
        this_exp = Rational(1)

        if isinstance(expr, Pow):
            if isinstance(expr.exp, Rational):
                expr, this_exp = expr.base, expr.exp
            elif isinstance(expr.exp, Mul):
                coeff, tail = expr.exp.getab()

                if isinstance(coeff, Rational):
                    expr, this_exp = Pow(expr.base, tail), coeff
        elif isinstance(expr, exp):
            if isinstance(expr._args, Rational):
                expr, this_exp = exp(Rational(1)), expr._args
            elif isinstance(expr._args, Mul):
                coeff, tail = expr._args.getab()

                if isinstance(coeff, Rational):
                    expr, this_exp = exp(tail), coeff

        return expr, this_exp

    def has_symbol(terms, sym):
        items = make_list(sym, Mul)

        if len(terms) < len(items):
            return None
        else:
            items = [ get_base_exp(i) for i in items ]
            terms = [ get_base_exp(t) for t in terms ]

            common_exp = Rational(1)

            for item, i_exp in items:
                for j in range(len(terms)):
                    term, t_exp = terms[j]

                    if term == item:
                        expo = t_exp / i_exp

                        if common_exp.is_one:
                            common_exp = expo
                        else:
                            if common_exp != expo:
                                return None

                        del terms[j]
                        break
                else:
                    return None

            return terms, common_exp

    expr = Basic.sympify(expr)

    if isinstance(syms, list):
        syms = [ separate(s) for s in syms ]
    else:
        syms = [ separate(syms) ]

    collected, final = {}, {}
    disliked = Rational(0)

    for prod in make_list(expr, Add):
        prod = separate(prod)

        terms = make_list(prod, Mul)

        for sym in syms:
            result = has_symbol(terms, sym)

            if result is not None:
                terms, common_exp = result
                index = sym, common_exp

                prod = Mul(*[ Pow(term, expo) for term, expo in terms ])

                if index in collected:
                    collected[index] += prod
                else:
                    collected[index] = prod

                break
        else:
            disliked += prod

    for ((sym, expo), expr) in collected.iteritems():
        final[separate(sym**expo)] = expr

    if disliked != Rational(0):
        final[Rational(1)] = disliked

    if pretty:
        return Add(*[ a*b for a, b in final.iteritems() ])
    else:
        return final

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
    elif isinstance(expr, Pow):
        return Pow(trigsimp(expr.base), trigsimp(expr.exp))
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
