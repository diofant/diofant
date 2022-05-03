import functools
import math
from collections import defaultdict

from ..core import (Add, Basic, Dummy, E, Expr, FunctionClass, I, Integer, Mul,
                    Pow, Rational, Wild, cacheit, count_ops, expand,
                    expand_mul, factor_terms, symbols)
from ..core.compatibility import iterable
from ..core.function import _mexpand
from ..core.strategies import greedy, identity
from ..core.sympify import sympify
from ..domains import ZZ
from ..functions import cos, cosh, cot, coth, exp, sin, sinh, tan, tanh
from ..functions.elementary.hyperbolic import HyperbolicFunction
from ..functions.elementary.trigonometric import TrigonometricFunction
from ..polys import Poly, cancel, factor, parallel_poly_from_expr
from ..polys.polyerrors import PolificationFailed
from ..polys.polytools import groebner
from .cse_main import cse


def trigsimp_groebner(expr, hints=[], quick=False, order='grlex',
                      polynomial=False):
    """
    Simplify trigonometric expressions using a groebner basis algorithm.

    This routine takes a fraction involving trigonometric or hyperbolic
    expressions, and tries to simplify it. The primary metric is the
    total degree. Some attempts are made to choose the simplest possible
    expression of the minimal degree, but this is non-rigorous, and also
    very slow (see the ``quick=True`` option).

    If ``polynomial`` is set to True, instead of simplifying numerator and
    denominator together, this function just brings numerator and denominator
    into a canonical form. This is much faster, but has potentially worse
    results. However, if the input is a polynomial, then the result is
    guaranteed to be an equivalent polynomial of minimal degree.

    The most important option is hints. Its entries can be any of the
    following:

    - a natural number
    - a function
    - an iterable of the form (func, var1, var2, ...)
    - anything else, interpreted as a generator

    A number is used to indicate that the search space should be increased.
    A function is used to indicate that said function is likely to occur in a
    simplified expression.
    An iterable is used indicate that func(var1 + var2 + ...) is likely to
    occur in a simplified .
    An additional generator also indicates that it is likely to occur.
    (See examples below).

    This routine carries out various computationally intensive algorithms.
    The option ``quick=True`` can be used to suppress one particularly slow
    step (at the expense of potentially more complicated results, but never at
    the expense of increased total degree).

    Examples
    ========

    Suppose you want to simplify ``sin(x)*cos(x)``. Naively, nothing happens:

    >>> ex = sin(x)*cos(x)
    >>> trigsimp_groebner(ex)
    sin(x)*cos(x)

    This is because ``trigsimp_groebner`` only looks for a simplification
    involving just ``sin(x)`` and ``cos(x)``. You can tell it to also try
    ``2*x`` by passing ``hints=[2]``:

    >>> trigsimp_groebner(ex, hints=[2])
    sin(2*x)/2
    >>> trigsimp_groebner(sin(x)**2 - cos(x)**2, hints=[2])
    -cos(2*x)

    Increasing the search space this way can quickly become expensive. A much
    faster way is to give a specific expression that is likely to occur:

    >>> trigsimp_groebner(ex, hints=[sin(2*x)])
    sin(2*x)/2

    Hyperbolic expressions are similarly supported:

    >>> trigsimp_groebner(sinh(2*x)/sinh(x))
    2*cosh(x)

    Note how no hints had to be passed, since the expression already involved
    ``2*x``.

    The tangent function is also supported. You can either pass ``tan`` in the
    hints, to indicate that than should be tried whenever cosine or sine are,
    or you can pass a specific generator:

    >>> trigsimp_groebner(sin(x)/cos(x), hints=[tan])
    tan(x)
    >>> trigsimp_groebner(sinh(x)/cosh(x), hints=[tanh(x)])
    tanh(x)

    Finally, you can use the iterable form to suggest that angle sum formulae
    should be tried:

    >>> ex = (tan(x) + tan(y))/(1 - tan(x)*tan(y))
    >>> trigsimp_groebner(ex, hints=[(tan, x, y)])
    tan(x + y)

    """
    # TODO
    #  - preprocess by replacing everything by funcs we can handle
    # - optionally use cot instead of tan
    # - more intelligent hinting.
    #     For example, if the ideal is small, and we have sin(x), sin(y),
    #     add sin(x + y) automatically... ?
    # - algebraic numbers ...
    # - expressions of lowest degree are not distinguished properly
    #   e.g. 1 - sin(x)**2
    # - we could try to order the generators intelligently, so as to influence
    #   which monomials appear in the quotient basis

    # THEORY
    # ------
    # Ratsimpmodprime above can be used to "simplify" a rational function
    # modulo a prime ideal. "Simplify" mainly means finding an equivalent
    # expression of lower total degree.
    #
    # We intend to use this to simplify trigonometric functions. To do that,
    # we need to decide (a) which ring to use, and (b) modulo which ideal to
    # simplify. In practice, (a) means settling on a list of "generators"
    # a, b, c, ..., such that the fraction we want to simplify is a rational
    # function in a, b, c, ..., with coefficients in ZZ (integers).
    # (2) means that we have to decide what relations to impose on the
    # generators. There are two practical problems:
    #   (1) The ideal has to be *prime* (a technical term).
    #   (2) The relations have to be polynomials in the generators.
    #
    # We typically have two kinds of generators:
    # - trigonometric expressions, like sin(x), cos(5*x), etc
    # - "everything else", like gamma(x), pi, etc.
    #
    # Since this function is trigsimp, we will concentrate on what to do with
    # trigonometric expressions. We can also simplify hyperbolic expressions,
    # but the extensions should be clear.
    #
    # One crucial point is that all *other* generators really should behave
    # like indeterminates. In particular if (say) "I" is one of them, then
    # in fact I**2 + 1 = 0 and we may and will compute non-sensical
    # expressions. However, we can work with a dummy and add the relation
    # I**2 + 1 = 0 to our ideal, then substitute back in the end.
    #
    # Now regarding trigonometric generators. We split them into groups,
    # according to the argument of the trigonometric functions. We want to
    # organise this in such a way that most trigonometric identities apply in
    # the same group. For example, given sin(x), cos(2*x) and cos(y), we would
    # group as [sin(x), cos(2*x)] and [cos(y)].
    #
    # Our prime ideal will be built in three steps:
    # (1) For each group, compute a "geometrically prime" ideal of relations.
    #     Geometrically prime means that it generates a prime ideal in
    #     CC[gens], not just ZZ[gens].
    # (2) Take the union of all the generators of the ideals for all groups.
    #     By the geometric primality condition, this is still prime.
    # (3) Add further inter-group relations which preserve primality.
    #
    # Step (1) works as follows. We will isolate common factors in the
    # argument, so that all our generators are of the form sin(n*x), cos(n*x)
    # or tan(n*x), with n an integer. Suppose first there are no tan terms.
    # The ideal [sin(x)**2 + cos(x)**2 - 1] is geometrically prime, since
    # X**2 + Y**2 - 1 is irreducible over CC.
    # Now, if we have a generator sin(n*x), than we can, using trig identities,
    # express sin(n*x) as a polynomial in sin(x) and cos(x). We can add this
    # relation to the ideal, preserving geometric primality, since the quotient
    # ring is unchanged.
    # Thus we have treated all sin and cos terms.
    # For tan(n*x), we add a relation tan(n*x)*cos(n*x) - sin(n*x) = 0.
    # (This requires of course that we already have relations for cos(n*x) and
    # sin(n*x).) It is not obvious, but it seems that this preserves geometric
    # primality.
    # XXX A real proof would be nice. HELP!
    #     Sketch that <S**2 + C**2 - 1, C*T - S> is a prime ideal of
    #     CC[S, C, T]:
    #     - it suffices to show that the projective closure in CP**3 is
    #       irreducible
    #     - using the half-angle substitutions, we can express sin(x), tan(x),
    #       cos(x) as rational functions in tan(x/2)
    #     - from this, we get a rational map from CP**1 to our curve
    #     - this is a morphism, hence the curve is prime
    #
    # Step (2) is trivial.
    #
    # Step (3) works by adding selected relations of the form
    # sin(x + y) - sin(x)*cos(y) - sin(y)*cos(x), etc. Geometric primality is
    # preserved by the same argument as before.

    def parse_hints(hints):
        """Split hints into (n, funcs, iterables, gens)."""
        n = 1
        funcs, iterables, gens = [], [], []
        for e in hints:
            if isinstance(e, (int, Integer)):
                n = e
            elif isinstance(e, FunctionClass):
                funcs.append(e)
            elif iterable(e):
                iterables.append((e[0], e[1:]))
                # XXX sin(x+2y)?
                # Note: we go through polys so e.g.
                # sin(-x) -> -sin(x) -> sin(x)
                gens.extend(parallel_poly_from_expr(
                    [e[0](x) for x in e[1:]] + [e[0](Add(*e[1:]))])[1].gens)
            else:
                gens.append(e)
        return n, funcs, iterables, gens

    def build_ideal(x, terms):
        """
        Build generators for our ideal. Terms is an iterable with elements of
        the form (fn, coeff), indicating that we have a generator fn(coeff*x).

        If any of the terms is trigonometric, sin(x) and cos(x) are guaranteed
        to appear in terms. Similarly for hyperbolic functions. For tan(n*x),
        sin(n*x) and cos(n*x) are guaranteed.

        """
        I = []
        y = Dummy('y')
        for fn, coeff in terms:
            for c, s, t, rel in (
                    [cos, sin, tan, cos(x)**2 + sin(x)**2 - 1],
                    [cosh, sinh, tanh, cosh(x)**2 - sinh(x)**2 - 1]):
                if coeff == 1 and fn in [c, s]:
                    I.append(rel)
                elif fn == t:
                    I.append(t(coeff*x)*c(coeff*x) - s(coeff*x))
                elif fn in [c, s]:
                    cn = fn(coeff*y).expand(trig=True).subs({y: x})
                    I.append(fn(coeff*x) - cn)
        return list(set(I))

    def analyse_gens(gens, hints):
        """
        Analyse the generators ``gens``, using the hints ``hints``.

        The meaning of ``hints`` is described in the main docstring.
        Return a new list of generators, and also the ideal we should
        work with.

        """
        # First parse the hints
        n, funcs, iterables, extragens = parse_hints(hints)

        # We just add the extragens to gens and analyse them as before
        gens = list(gens)
        gens.extend(extragens)

        # remove duplicates
        funcs = list(set(funcs))
        iterables = list(set(iterables))
        gens = list(set(gens))

        # all the functions we can do anything with
        allfuncs = {sin, cos, tan, sinh, cosh, tanh}
        # sin(3*x) -> ((3, x), sin)
        trigterms = [(g.args[0].as_coeff_mul(), g.func) for g in gens
                     if g.func in allfuncs]
        # Our list of new generators - start with anything that we cannot
        # work with (i.e. is not a trigonometric term)
        freegens = [g for g in gens if g.func not in allfuncs]
        newgens = []
        trigdict = defaultdict(list)
        for (coeff, var), fn in trigterms:
            trigdict[var].append((coeff, fn))
        res = []  # the ideal

        for key, val in trigdict.items():
            # We have now assembeled a dictionary. Its keys are common
            # arguments in trigonometric expressions, and values are lists of
            # pairs (fn, coeff). x0, (fn, coeff) in trigdict means that we
            # need to deal with fn(coeff*x0). We take the rational gcd of the
            # coeffs, call it ``gcd``. We then use x = x0/gcd as "base symbol",
            # all other arguments are integral multiples thereof.
            # We will build an ideal which works with sin(x), cos(x).
            # If hint tan is provided, also work with tan(x). Moreover, if
            # n > 1, also work with sin(k*x) for k <= n, and similarly for cos
            # (and tan if the hint is provided). Finally, any generators which
            # the ideal does not work with but we need to accomodate (either
            # because it was in expr or because it was provided as a hint)
            # we also build into the ideal.
            # This selection process is expressed in the list ``terms``.
            # build_ideal then generates the actual relations in our ideal,
            # from this list.
            fns = [x[1] for x in val]
            val = [x[0] for x in val]
            gcd = functools.reduce(math.gcd, val)
            terms = [(fn, v/gcd) for (fn, v) in zip(fns, val)]
            fs = set(funcs + fns)
            for c, s, t in ([cos, sin, tan], [cosh, sinh, tanh]):
                if any(x in fs for x in (c, s, t)):
                    fs.add(c)
                    fs.add(s)
            for fn in fs:
                for k in range(1, n + 1):
                    terms.append((fn, k))
            extra = []
            for fn, v in terms:
                if fn == tan:
                    extra.append((sin, v))
                    extra.append((cos, v))
                if fn in [sin, cos] and tan in fs:
                    extra.append((tan, v))
                if fn == tanh:
                    extra.append((sinh, v))
                    extra.append((cosh, v))
                if fn in [sinh, cosh] and tanh in fs:
                    extra.append((tanh, v))
            terms.extend(extra)
            x = gcd*Mul(*key)
            r = build_ideal(x, terms)
            res.extend(r)
            newgens.extend({fn(v*x) for fn, v in terms})

        # Add generators for compound expressions from iterables
        for fn, args in iterables:
            if fn == tan:
                # Tan expressions are recovered from sin and cos.
                iterables.extend([(sin, args), (cos, args)])
            elif fn == tanh:
                # Tanh expressions are recovered from sihn and cosh.
                iterables.extend([(sinh, args), (cosh, args)])
            else:
                dummys = symbols(f'd:{len(args):d}', cls=Dummy)
                expr = fn(Add(*dummys)).expand(trig=True).subs(list(zip(dummys, args)))
                res.append(fn(Add(*args)) - expr)

        if myI in gens:
            res.append(myI**2 + 1)
            freegens.remove(myI)
            newgens.append(myI)

        return res, freegens, newgens

    myI = Dummy('I')
    expr = expr.subs({I: myI})
    subs = [(myI, I)]

    num, denom = cancel(expr).as_numer_denom()
    try:
        (_, pdenom), opt = parallel_poly_from_expr([num, denom])
    except PolificationFailed:
        return expr
    ideal, freegens, gens = analyse_gens(opt.gens, hints)
    # NOTE we force the domain to be ZZ to stop polys from injecting generators
    #      (which is usually a sign of a bug in the way we build the ideal)
    if not gens:
        return expr
    G = groebner(ideal, order=order, gens=gens, domain=ZZ)

    # If our fraction is a polynomial in the free generators, simplify all
    # coefficients separately:

    from .ratsimp import ratsimpmodprime

    if freegens and pdenom.has_only_gens(*set(gens).intersection(pdenom.gens)):
        num = Poly(num, gens=gens+freegens).eject(*gens)
        res = []
        for monom, coeff in num.terms():
            ourgens = set(parallel_poly_from_expr([coeff, denom])[1].gens)
            # We compute the transitive closure of all generators that can
            # be reached from our generators through relations in the ideal.
            changed = True
            while changed:
                changed = False
                for p in ideal:
                    p = Poly(p)
                    if not ourgens.issuperset(p.gens) and \
                       not p.has_only_gens(*set(p.gens).difference(ourgens)):
                        changed = True
                        ourgens.update(p.exclude().gens)
            # NOTE preserve order!
            realgens = [x for x in gens if x in ourgens]
            # The generators of the ideal have now been (implicitely) split
            # into two groups: those involving ourgens and those that don't.
            # Since we took the transitive closure above, these two groups
            # live in subgrings generated by a *disjoint* set of variables.
            # Any sensible groebner basis algorithm will preserve this disjoint
            # structure (i.e. the elements of the groebner basis can be split
            # similarly), and and the two subsets of the groebner basis then
            # form groebner bases by themselves. (For the smaller generating
            # sets, of course.)
            ourG = [g.as_expr() for g in G.polys if
                    g.has_only_gens(*ourgens.intersection(g.gens))]
            res.append(Mul(*[a**b for a, b in zip(freegens, monom)]) *
                       ratsimpmodprime(coeff/denom, ourG, order=order,
                                       gens=realgens, quick=quick, domain=ZZ,
                                       polynomial=polynomial).subs(subs))
        return Add(*res)
        # NOTE The following is simpler and has less assumptions on the
        #      groebner basis algorithm. If the above turns out to be broken,
        #      use this.
        # return Add(*[Mul(*[a**b for a, b in zip(freegens, monom)]) *
        #              ratsimpmodprime(coeff/denom, list(G), order=order,
        #                              gens=gens, quick=quick, domain=ZZ)
        #              for monom, coeff in num.terms()])
    else:
        return ratsimpmodprime(
            expr, list(G), order=order, gens=freegens+gens,
            quick=quick, domain=ZZ, polynomial=polynomial).subs(subs)


_trigs = (TrigonometricFunction, HyperbolicFunction)


def trigsimp(expr, **opts):
    """
    Reduces expression by using known trig identities.

    Notes
    =====

    method:
    - Determine the method to use. Valid choices are 'matching' (default),
    'groebner', 'combined', and 'fu'. If 'matching', simplify the
    expression recursively by targeting common patterns. If 'groebner', apply
    an experimental groebner basis algorithm. In this case further options
    are forwarded to ``trigsimp_groebner``, please refer to its docstring.
    If 'combined', first run the groebner basis algorithm with small
    default parameters, then run the 'matching' algorithm. 'fu' runs the
    collection of trigonometric transformations described by Fu, et al.

    Examples
    ========

    >>> e = 2*sin(x)**2 + 2*cos(x)**2
    >>> trigsimp(e)
    2

    Simplification occurs wherever trigonometric functions are located.

    >>> trigsimp(log(e))
    log(2)

    Using `method="groebner"` (or `"combined"`) might lead to greater
    simplification.

    The old trigsimp routine can be accessed as with method 'old'.

    >>> t = 3*tanh(x)**7 - 2/coth(x)**7
    >>> trigsimp(t, method='old') == t
    True
    >>> trigsimp(t)
    tanh(x)**7

    See Also
    ========

    diofant.simplify.fu.fu

    """
    from ..matrices import MatrixExpr
    from .fu import fu

    expr = sympify(expr)

    try:
        return expr._eval_trigsimp(**opts)
    except AttributeError:
        pass

    old = opts.pop('old', False)
    if not old:
        opts.pop('deep', None)
        method = opts.pop('method', 'matching')
    else:
        method = 'old'

    def groebnersimp(ex, **opts):
        def traverse(e):
            if e.is_Atom:
                return e
            args = [traverse(x) for x in e.args]
            if e.is_Function or e.is_Pow:
                args = [trigsimp_groebner(x, **opts) for x in args]
            return e.func(*args)
        new = traverse(ex)
        if not isinstance(new, Expr) or isinstance(new, MatrixExpr):
            return new
        return trigsimp_groebner(new, **opts)

    trigsimpfunc = {
        'fu': (functools.partial(fu, **opts)),
        'matching': (futrig),
        'groebner': (functools.partial(groebnersimp, **opts)),
        'combined': (lambda x: futrig(groebnersimp(x,
                                                   polynomial=True, hints=[2, tan]))),
        'old': functools.partial(trigsimp_old, **opts)}[method]

    return trigsimpfunc(expr)


def exptrigsimp(expr, simplify=True):
    """
    Simplifies exponential / trigonometric / hyperbolic functions.
    When ``simplify`` is True (default) the expression obtained after the
    simplification step will be then be passed through simplify to
    precondition it so the final transformations will be applied.

    Examples
    ========

    >>> exptrigsimp(exp(z) + exp(-z))
    2*cosh(z)
    >>> exptrigsimp(cosh(z) - sinh(z))
    E**(-z)

    """
    from .fu import TR2i, hyper_as_trig
    from .simplify import bottom_up

    def exp_trig(e):
        # select the better of e, and e rewritten in terms of exp or trig
        # functions
        choices = [e]
        if e.has(*_trigs):
            choices.append(e.rewrite(exp))
        choices.append(e.rewrite(cos))
        return min(*choices, key=count_ops)
    newexpr = bottom_up(expr, exp_trig)

    if simplify:
        newexpr = newexpr.simplify()

    # conversion from exp to hyperbolic
    ex = {a for a in newexpr.atoms(Pow) if a.is_Exp} | newexpr.atoms(E)
    if ex:
        ex0 = {list(ex)[0]}
        ex = [ei for ei in ex if 1/ei not in ex]
        if not ex:
            ex = ex0

    # sinh and cosh
    for ei in ex:
        a = ei.exp if ei is not E else Integer(1)
        newexpr = newexpr.subs({ei + 1/ei: 2*cosh(a)})
        newexpr = newexpr.subs({ei - 1/ei: 2*sinh(a)})
        e2 = ei**-2
        if e2 in ex:
            a = e2.exp/2 if e2 is not E else Rational(1, 2)
            newexpr = newexpr.subs({(e2 + 1)*ei: 2*cosh(a)})
            newexpr = newexpr.subs({(e2 - 1)*ei: 2*sinh(a)})

    # exp ratios to tan and tanh
    for ei in ex:
        n, d = ei - 1, ei + 1
        et = n/d
        etinv = d/n  # not 1/et or else recursion errors arise
        a = ei.exp if ei.is_Exp else Integer(1)
        if a.is_Mul or a is I:
            c = a.as_coefficient(I)
            if c:
                t = I*tan(c/2)
                newexpr = newexpr.subs({etinv: 1/t})
                newexpr = newexpr.subs({et: t})
                continue
        t = tanh(a/2)
        newexpr = newexpr.subs({etinv: 1/t})
        newexpr = newexpr.subs({et: t})

    # sin/cos and sinh/cosh ratios to tan and tanh, respectively
    if newexpr.has(HyperbolicFunction):
        e, f = hyper_as_trig(newexpr)
        newexpr = f(TR2i(e))
    if newexpr.has(TrigonometricFunction):
        newexpr = TR2i(newexpr)

    # can we ever generate an I where there was none previously?
    if not (newexpr.has(I) and not expr.has(I)):
        expr = newexpr
    return expr


# -------------------- the old trigsimp routines ---------------------


def trigsimp_old(expr, **opts):
    """
    Reduces expression by using known trig identities.

    Notes
    =====

    deep:
    - Apply trigsimp inside all objects with arguments

    recursive:
    - Use common subexpression elimination (cse()) and apply
    trigsimp recursively (this is quite expensive if the
    expression is large)

    method:
    - Determine the method to use. Valid choices are 'matching' (default),
    'groebner', 'combined', 'fu' and 'futrig'. If 'matching', simplify the
    expression recursively by pattern matching. If 'groebner', apply an
    experimental groebner basis algorithm. In this case further options
    are forwarded to ``trigsimp_groebner``, please refer to its docstring.
    If 'combined', first run the groebner basis algorithm with small
    default parameters, then run the 'matching' algorithm. 'fu' runs the
    collection of trigonometric transformations described by Fu, et al.,
    while `futrig` runs a subset of Fu-transfoms that mimic the behavior
    of `trigsimp`.

    compare:
    - show input and output from `trigsimp` and `futrig` when different,
    but returns the `trigsimp` value.

    Examples
    ========

    >>> e = 2*sin(x)**2 + 2*cos(x)**2
    >>> trigsimp(e, old=True)
    2
    >>> trigsimp(log(e), old=True)
    log(2*sin(x)**2 + 2*cos(x)**2)
    >>> trigsimp(log(e), deep=True, old=True)
    log(2)

    Using `method="groebner"` (or `"combined"`) can sometimes lead to a lot
    more simplification:

    >>> e = (-sin(x) + 1)/cos(x) + cos(x)/(-sin(x) + 1)
    >>> trigsimp(e, old=True)
    (-sin(x) + 1)/cos(x) - cos(x)/(sin(x) - 1)
    >>> trigsimp(e, method='groebner', old=True)
    2/cos(x)

    >>> trigsimp(1/cot(x)**2, compare=True, old=True)
          futrig: tan(x)**2
    cot(x)**(-2)

    See Also
    ========

    diofant.simplify.fu.fu
    futrig

    """
    from .simplify import separatevars

    old = expr
    first = opts.pop('first', True)
    if first:
        if not expr.has(*_trigs):
            return expr

        trigsyms = set().union(*[t.free_symbols for t in expr.atoms(*_trigs)])
        if len(trigsyms) > 1:
            d = separatevars(expr)
            if d.is_Mul:
                d = separatevars(d, dict=True) or d
            if isinstance(d, dict):
                expr = 1
                for v in d.values():
                    # remove hollow factoring
                    was = v
                    v = expand_mul(v)
                    opts['first'] = False
                    vnew = trigsimp(v, **opts)
                    if vnew == v:
                        vnew = was
                    expr *= vnew
                old = expr
            else:
                if d.is_Add:
                    for s in trigsyms:
                        r, e = expr.as_independent(s)
                        if r:
                            opts['first'] = False
                            expr = r + trigsimp(e, **opts)
                            if not expr.is_Add:
                                break
                    old = expr

    recursive = opts.pop('recursive', False)
    deep = opts.pop('deep', False)
    method = opts.pop('method', 'matching')

    def groebnersimp(ex, deep, **opts):
        def traverse(e):
            if e.is_Atom:
                return e
            args = [traverse(x) for x in e.args]
            if e.is_Function or e.is_Pow:
                args = [trigsimp_groebner(x, **opts) for x in args]
            return e.func(*args)
        if deep:
            ex = traverse(ex)
        return trigsimp_groebner(ex, **opts)

    trigsimpfunc = {
        'matching': (_trigsimp),
        'groebner': (functools.partial(groebnersimp, **opts)),
        'combined': (lambda x, d: _trigsimp(groebnersimp(x,
                                                         d, polynomial=True, hints=[2, tan]),
                                            d))}[method]

    if recursive:
        w, g = cse(expr)
        g = trigsimpfunc(g[0], deep)

        for sub in reversed(w):
            g = g.subs({sub[0]: sub[1]})
            g = trigsimpfunc(g, deep)
        result = g
    else:
        result = trigsimpfunc(expr, deep)

    if opts.get('compare', False):
        f = futrig(old)
        if f != result:
            print('\tfutrig:', f)

    return result


def _dotrig(a, b):
    """Helper to tell whether ``a`` and ``b`` have the same sorts
    of symbols in them -- no need to test hyperbolic patterns against
    expressions that have no hyperbolics in them.

    """
    return a.func == b.func and (
        a.has(TrigonometricFunction) and b.has(TrigonometricFunction) or
        a.has(HyperbolicFunction) and b.has(HyperbolicFunction))


_trigpat = None


def _trigpats():
    global _trigpat
    a, b, c = symbols('a b c', cls=Wild)
    d = Wild('d', commutative=False)

    # for the simplifications like sinh/cosh -> tanh:
    # DO NOT REORDER THE FIRST 14 since these are assumed to be in this
    # order in _match_div_rewrite.
    matchers_division = (
        (a*sin(b)**c/cos(b)**c, a*tan(b)**c, sin(b), cos(b)),
        (a*tan(b)**c*cos(b)**c, a*sin(b)**c, sin(b), cos(b)),
        (a*cot(b)**c*sin(b)**c, a*cos(b)**c, sin(b), cos(b)),
        (a*tan(b)**c/sin(b)**c, a/cos(b)**c, sin(b), cos(b)),
        (a*cot(b)**c/cos(b)**c, a/sin(b)**c, sin(b), cos(b)),
        (a*cot(b)**c*tan(b)**c, a, sin(b), cos(b)),
        (a*(cos(b) + 1)**c*(cos(b) - 1)**c,
            a*(-sin(b)**2)**c, cos(b) + 1, cos(b) - 1),
        (a*(sin(b) + 1)**c*(sin(b) - 1)**c,
            a*(-cos(b)**2)**c, sin(b) + 1, sin(b) - 1),

        (a*sinh(b)**c/cosh(b)**c, a*tanh(b)**c, Integer(1), Integer(1)),
        (a*tanh(b)**c*cosh(b)**c, a*sinh(b)**c, Integer(1), Integer(1)),
        (a*coth(b)**c*sinh(b)**c, a*cosh(b)**c, Integer(1), Integer(1)),
        (a*tanh(b)**c/sinh(b)**c, a/cosh(b)**c, Integer(1), Integer(1)),
        (a*coth(b)**c/cosh(b)**c, a/sinh(b)**c, Integer(1), Integer(1)),
        (a*coth(b)**c*tanh(b)**c, a, Integer(1), Integer(1)),

        (c*(tanh(a) + tanh(b))/(1 + tanh(a)*tanh(b)),
            tanh(a + b)*c, Integer(1), Integer(1)),
    )

    matchers_add = (
        (c*sin(a)*cos(b) + c*cos(a)*sin(b) + d, sin(a + b)*c + d),
        (c*cos(a)*cos(b) - c*sin(a)*sin(b) + d, cos(a + b)*c + d),
        (c*sin(a)*cos(b) - c*cos(a)*sin(b) + d, sin(a - b)*c + d),
        (c*cos(a)*cos(b) + c*sin(a)*sin(b) + d, cos(a - b)*c + d),
        (c*sinh(a)*cosh(b) + c*sinh(b)*cosh(a) + d, sinh(a + b)*c + d),
        (c*cosh(a)*cosh(b) + c*sinh(a)*sinh(b) + d, cosh(a + b)*c + d),
    )

    # for cos(x)**2 + sin(x)**2 -> 1
    matchers_identity = (
        (a*sin(b)**2, a - a*cos(b)**2),
        (a*tan(b)**2, a*(1/cos(b))**2 - a),
        (a*cot(b)**2, a*(1/sin(b))**2 - a),
        (a*sin(b + c), a*(sin(b)*cos(c) + sin(c)*cos(b))),
        (a*cos(b + c), a*(cos(b)*cos(c) - sin(b)*sin(c))),
        (a*tan(b + c), a*((tan(b) + tan(c))/(1 - tan(b)*tan(c)))),

        (a*sinh(b)**2, a*cosh(b)**2 - a),
        (a*tanh(b)**2, a - a*(1/cosh(b))**2),
        (a*coth(b)**2, a + a*(1/sinh(b))**2),
        (a*sinh(b + c), a*(sinh(b)*cosh(c) + sinh(c)*cosh(b))),
        (a*cosh(b + c), a*(cosh(b)*cosh(c) + sinh(b)*sinh(c))),
        (a*tanh(b + c), a*((tanh(b) + tanh(c))/(1 + tanh(b)*tanh(c)))),

    )

    # Reduce any lingering artifacts, such as sin(x)**2 changing
    # to 1-cos(x)**2 when sin(x)**2 was "simpler"
    artifacts = (
        (a - a*cos(b)**2 + c, a*sin(b)**2 + c, cos),
        (a - a*(1/cos(b))**2 + c, -a*tan(b)**2 + c, cos),
        (a - a*(1/sin(b))**2 + c, -a*cot(b)**2 + c, sin),

        (a - a*cosh(b)**2 + c, -a*sinh(b)**2 + c, cosh),
        (a - a*(1/cosh(b))**2 + c, a*tanh(b)**2 + c, cosh),
        (a + a*(1/sinh(b))**2 + c, a*coth(b)**2 + c, sinh),

        # same as above but with noncommutative prefactor
        (a*d - a*d*cos(b)**2 + c, a*d*sin(b)**2 + c, cos),
        (a*d - a*d*(1/cos(b))**2 + c, -a*d*tan(b)**2 + c, cos),
        (a*d - a*d*(1/sin(b))**2 + c, -a*d*cot(b)**2 + c, sin),

        (a*d - a*d*cosh(b)**2 + c, -a*d*sinh(b)**2 + c, cosh),
        (a*d - a*d*(1/cosh(b))**2 + c, a*d*tanh(b)**2 + c, cosh),
        (a*d + a*d*(1/sinh(b))**2 + c, a*d*coth(b)**2 + c, sinh),
    )

    _trigpat = (a, b, c, d, matchers_division, matchers_add,
                matchers_identity, artifacts)
    return _trigpat


def _replace_mul_fpowxgpow(expr, f, g, rexp, h, rexph):
    """Helper for _match_div_rewrite.

    Replace f(b_)**c_*g(b_)**(rexp(c_)) with h(b)**rexph(c) if f(b_)
    and g(b_) are both positive or if c_ is an integer.

    """
    # assert expr.is_Mul and expr.is_commutative and f != g
    fargs = defaultdict(int)
    gargs = defaultdict(int)
    args = []
    for x in expr.args:
        if x.is_Pow or x.func in (f, g):
            b, e = x.as_base_exp()
            if b.is_positive or e.is_integer:
                if b.func == f:
                    fargs[b.args[0]] += e
                    continue
                if b.func == g:
                    gargs[b.args[0]] += e
                    continue
        args.append(x)
    common = set(fargs) & set(gargs)
    hit = False
    while common:
        key = common.pop()
        fe = fargs.pop(key)
        ge = gargs.pop(key)
        if fe == rexp(ge):
            args.append(h(key)**rexph(fe))
            hit = True
        else:
            fargs[key] = fe
            gargs[key] = ge
    if not hit:
        return expr
    while fargs:
        key, e = fargs.popitem()
        args.append(f(key)**e)
    while gargs:
        key, e = gargs.popitem()
        args.append(g(key)**e)
    return Mul(*args)


def _idn(x):
    return x


def _midn(x):
    return -x


def _one(x):
    return Integer(1)


def _match_div_rewrite(expr, i):
    """Helper for __trigsimp."""
    if i == 0:
        expr = _replace_mul_fpowxgpow(expr, sin, cos,
                                      _midn, tan, _idn)
    elif i == 1:
        expr = _replace_mul_fpowxgpow(expr, tan, cos,
                                      _idn, sin, _idn)
    elif i == 2:
        expr = _replace_mul_fpowxgpow(expr, cot, sin,
                                      _idn, cos, _idn)
    elif i == 3:
        expr = _replace_mul_fpowxgpow(expr, tan, sin,
                                      _midn, cos, _midn)
    elif i == 4:
        expr = _replace_mul_fpowxgpow(expr, cot, cos,
                                      _midn, sin, _midn)
    elif i == 5:
        expr = _replace_mul_fpowxgpow(expr, cot, tan,
                                      _idn, _one, _idn)
    # i in (6, 7) is skipped
    elif i == 8:
        expr = _replace_mul_fpowxgpow(expr, sinh, cosh,
                                      _midn, tanh, _idn)
    elif i == 9:
        expr = _replace_mul_fpowxgpow(expr, tanh, cosh,
                                      _idn, sinh, _idn)
    elif i == 10:
        expr = _replace_mul_fpowxgpow(expr, coth, sinh,
                                      _idn, cosh, _idn)
    elif i == 11:
        expr = _replace_mul_fpowxgpow(expr, tanh, sinh,
                                      _midn, cosh, _midn)
    elif i == 12:
        expr = _replace_mul_fpowxgpow(expr, coth, cosh,
                                      _midn, sinh, _midn)
    elif i == 13:
        expr = _replace_mul_fpowxgpow(expr, coth, tanh,
                                      _idn, _one, _idn)
    else:
        return
    return expr


def _trigsimp(expr, deep=False):
    # protect the cache from non-trig patterns; we only allow
    # trig patterns to enter the cache
    if expr.has(*_trigs):
        return __trigsimp(expr, deep)
    return expr


@cacheit
def __trigsimp(expr, deep=False):
    """Recursive helper for trigsimp."""
    from .fu import TR10i

    if _trigpat is None:
        _trigpats()
    a, b, c, d, matchers_division, matchers_add, \
        matchers_identity, artifacts = _trigpat

    if expr.is_Mul:
        # do some simplifications like sin/cos -> tan:
        if not expr.is_commutative:
            com, nc = expr.args_cnc()
            expr = _trigsimp(Mul._from_args(com), deep)*Mul._from_args(nc)
        else:
            for i, (pattern, simp, ok1, ok2) in enumerate(matchers_division):
                if not _dotrig(expr, pattern):
                    continue

                newexpr = _match_div_rewrite(expr, i)
                if newexpr is not None:
                    if newexpr != expr:
                        expr = newexpr
                        break
                    continue

                # use Diofant matching instead
                res = expr.match(pattern)
                if res and res.get(c, 0):
                    if not res[c].is_integer:
                        ok = ok1.subs(res)
                        if not ok.is_positive:
                            continue
                        ok = ok2.subs(res)
                        if not ok.is_positive:
                            continue
                    # if "a" contains any of trig or hyperbolic funcs with
                    # argument "b" then skip the simplification
                    if any(w.args[0] == res[b] for w in res[a].atoms(
                            TrigonometricFunction, HyperbolicFunction)):
                        continue
                    # simplify and finish:
                    expr = simp.subs(res)
                    break  # process below

    if expr.is_Add:
        args = []
        for term in expr.args:
            if not term.is_commutative:
                com, nc = term.args_cnc()
                nc = Mul._from_args(nc)
                term = Mul._from_args(com)
            else:
                nc = Integer(1)
            term = _trigsimp(term, deep)
            for pattern, result in matchers_identity:
                res = term.match(pattern)
                if res is not None:
                    term = result.subs(res)
                    break
            args.append(term*nc)
        if args != expr.args:
            expr = Add(*args)
            expr = min(expr, expand(expr), key=count_ops)
        if expr.is_Add:
            for pattern, result in matchers_add:
                if not _dotrig(expr, pattern):
                    continue
                expr = TR10i(expr)
                if expr.has(HyperbolicFunction):
                    res = expr.match(pattern)
                    # if "d" contains any trig or hyperbolic funcs with
                    # argument "a" or "b" then skip the simplification;
                    # this isn't perfect -- see tests
                    if res is None or not (a in res and b in res) or any(
                        w.args[0] in (res[a], res[b]) for w in res[d].atoms(
                            TrigonometricFunction, HyperbolicFunction)):
                        continue
                    expr = result.subs(res)
                    break

        # Reduce any lingering artifacts, such as sin(x)**2 changing
        # to 1 - cos(x)**2 when sin(x)**2 was "simpler"
        for pattern, result, ex in artifacts:
            if not _dotrig(expr, pattern):
                continue
            # Substitute a new wild that excludes some function(s)
            # to help influence a better match. This is because
            # sometimes, for example, 'a' would match sec(x)**2
            a_t = Wild('a', exclude=[ex])
            pattern = pattern.subs({a: a_t})
            result = result.subs({a: a_t})

            m = expr.match(pattern)
            was = None
            while m and was != expr:
                was = expr
                if m[a_t] == 0 or \
                        -m[a_t] in m[c].args or m[a_t] + m[c] == 0:
                    break
                if d in m and m[a_t]*m[d] + m[c] == 0:
                    break
                expr = result.subs(m)
                m = expr.match(pattern)
                m.setdefault(c, Integer(0))

    elif expr.is_Mul or expr.is_Pow or deep and expr.args:
        expr = expr.func(*[_trigsimp(a, deep) for a in expr.args])

    try:
        if not expr.has(*_trigs):
            raise TypeError
        e = {a for a in expr.atoms(Pow) if a.is_Exp}
        new = expr.rewrite(exp, deep=deep)
        if new == e:
            raise TypeError
        fnew = factor(new)
        if fnew != new:
            new = sorted([new, factor(new)], key=count_ops)[0]
        # if all exp that were introduced disappeared then accept it
        ne = {a for a in new.atoms(Pow) if a.is_Exp}
        if not ne - e:
            expr = new
    except TypeError:
        pass

    return expr
# ------------------- end of old trigsimp routines --------------------


def futrig(e, **kwargs):
    """Return simplified ``e`` using Fu-like transformations.
    This is not the "Fu" algorithm. This is called by default
    from ``trigsimp``. By default, hyperbolics subexpressions
    will be simplified, but this can be disabled by setting
    ``hyper=False``.

    Examples
    ========

    >>> trigsimp(1/tan(x)**2)
    tan(x)**(-2)

    >>> futrig(sinh(x)/tanh(x))
    cosh(x)

    """
    from .fu import hyper_as_trig
    from .simplify import bottom_up

    e = sympify(e)

    if not isinstance(e, Basic):
        return e

    if not e.args:
        return e

    old = e
    e = bottom_up(e, functools.partial(_futrig, **kwargs))

    if kwargs.pop('hyper', True) and e.has(HyperbolicFunction):
        e, f = hyper_as_trig(e)
        e = f(_futrig(e))

    if e != old and e.is_Mul and e.args[0].is_Rational:
        # redistribute leading coeff on 2-arg Add
        e = Mul(*e.as_coeff_Mul())
    return e


def _futrig(e, **kwargs):
    """Helper for futrig."""
    from ..utilities.iterables import _nodes
    from .fu import (TR1, TR2, TR3, TR5, TR6, TR8, TR10, TR11, TR12, TR14,
                     TR15, TR16, TR22, TR111, L, TR2i, TR10i, TRmorrie)

    if not e.has(TrigonometricFunction):
        return e

    if e.is_Mul:
        coeff, e = e.as_independent(TrigonometricFunction)
    else:
        coeff = Integer(1)

    def Lops(x):
        return L(x), x.count_ops(), _nodes(x), len(x.args), x.is_Add

    def trigs(x):
        return x.has(TrigonometricFunction)

    tree = [identity,
            (
                TR3,  # canonical angles
                TR1,  # sec-csc -> cos-sin
                TR12,  # expand tan of sum
                lambda x: _eapply(factor, x, trigs),
                TR2,  # tan-cot -> sin-cos
                [identity, lambda x: _eapply(_mexpand, x, trigs)],
                TR2i,  # sin-cos ratio -> tan
                lambda x: _eapply(lambda i: factor(i.normal()), x, trigs),
                TR14,  # factored identities
                TR5,  # sin-pow -> cos_pow
                TR10,  # sin-cos of sums -> sin-cos prod
                TR11, TR6,  # reduce double angles and rewrite cos pows
                lambda x: _eapply(factor, x, trigs),
                TR14,  # factored powers of identities
                [identity, lambda x: _eapply(_mexpand, x, trigs)],
                TRmorrie,
                TR10i,  # sin-cos products > sin-cos of sums
                [identity, TR8],  # sin-cos products -> sin-cos of sums
                [identity, lambda x: TR2i(TR2(x))],  # tan -> sin-cos -> tan
                [
                    lambda x: _eapply(expand_mul, TR5(x), trigs),
                    lambda x: _eapply(
                        expand_mul, TR15(x), trigs)],  # pos/neg powers of sin
                [
                    lambda x:  _eapply(expand_mul, TR6(x), trigs),
                    lambda x:  _eapply(
                        expand_mul, TR16(x), trigs)],  # pos/neg powers of cos
                TR111,  # tan, sin, cos to neg power -> cot, csc, sec
                [identity, TR2i],  # sin-cos ratio to tan
                [identity, lambda x: _eapply(
                    expand_mul, TR22(x), trigs)],  # tan-cot to sec-csc
                TR1, TR2, TR2i,
                [identity, lambda x: _eapply(
                    factor_terms, TR12(x), trigs)],  # expand tan of sum
            )]
    e = greedy(tree, objective=Lops)(e)

    return coeff*e


def _is_Expr(e):
    """_eapply helper to tell whether ``e`` and all its args
    are Exprs.

    """
    if not isinstance(e, Expr):
        return False
    return all(_is_Expr(i) for i in e.args)


def _eapply(func, e, cond=None):
    """Apply ``func`` to ``e`` if all args are Exprs else only
    apply it to those args that *are* Exprs.

    """
    if not isinstance(e, Expr):
        return e
    if _is_Expr(e) or not e.args:
        return func(e)
    return e.func(*[
        _eapply(func, ei) if (cond is None or cond(ei)) else ei
        for ei in e.args])
