r"""
The Gruntz Algorithm
====================

This section explains the basics of the algorithm [1]_ used for computing
limits.  Most of the time the :py:func:`~sympy.series.limits.limit` function
should just work.  However it is still useful to keep in mind how it is
implemented in case something does not work as expected.

First we define an ordering on functions of single variable `x` according
to how rapidly varying they at infinity.  Any two functions `f(x)` and
`g(x)` can be compared using the properties of:

    .. math::
        L = \lim\limits_{x\to\infty}\frac{\log|f(x)|}{\log|g(x)|}

We shall say that `f(x)` *dominates* `g(x)`, written `f(x) \succ
g(x)`, iff `L=\pm\infty`.  We also say that `f(x)` and `g(x)` are *of the
same comparability class* if neither `f(x) \succ g(x)` nor `g(x) \succ
f(x)` and shall denote it as `f(x) \asymp g(x)`.

It is easy to show the following examples:

* `e^{e^x} \succ e^{x^2} \succ e^x \succ x \succ 42`
* `2 \asymp 3 \asymp -5`
* `x \asymp x^2 \asymp x^3 \asymp -x`
* `e^x \asymp e^{-x} \asymp e^{2x} \asymp e^{x + e^{-x}}`
* `f(x) \asymp 1/f(x)`

Using these definitions yields the following strategy for
computing `\lim_{x \to \infty} f(x)`:

1. Given the function `f(x)`, we find the set of *most rapidly varying
   subexpressions* (MRV set) of it.  All items of this set belongs to the
   same comparability class.  Let's say it is `\{e^x, e^{2x}\}`.

2. Choose an expression `\omega` which is positive and tends to zero and
   which is in the same comparability class as any element of the MRV set.
   Such element always exists.  Then we rewrite the MRV set using `\omega`,
   in our case `\{\omega^{-1}, \omega^{-2}\}`, and substitute it into `f(x)`.

3. Let `f(\omega)` be the function which is obtained from `f(x)` after the
   rewrite step above.  Consider all expressions independent of `\omega` as
   constants and compute the leading term of the power series of `f(\omega)`
   around `\omega = 0^+`:

       .. math:: f(\omega) = c_0 \omega^{e_0} + c_1 \omega^{e_1} + \dots

   where `e_0 < e_1 < e_2 \dots`

4. If the leading exponent `e_0 > 0` then the limit is `0`.  If `e_0 < 0`,
   then the answer is `\pm\infty` (depends on sign of `c_0`).  Finally,
   if `e_0 = 0`, the limit is the limit of the leading coefficient `c_0`.

Notes
-----

This exposition glossed over several details.  For example, limits could be
computed recursively (steps 1 and 4).  Please address to the Gruntz thesis [1]_
for proof of the termination (pp. 52-60).

References
----------

.. [1] `Gruntz Thesis <http://www.cybertester.com/data/gruntz.pdf>`_
"""

from __future__ import print_function, division

from sympy.core import S, oo, Symbol, I, Dummy, Wild, Mul
from sympy.core.compatibility import default_sort_key
from sympy.functions import log, exp
from sympy.series.order import Order
from sympy.simplify import powsimp
from sympy import cacheit
from sympy.core.compatibility import reduce

from sympy.utilities.timeutils import timethis
timeit = timethis('gruntz')


def compare(a, b, x):
    """Returns "<" if a<b, "=" for a == b, ">" for a>b"""
    # log(exp(...)) must always be simplified here for termination
    la, lb = log(a), log(b)
    if a.func is exp:
        la = a.exp
    if b.func is exp:
        lb = b.exp

    c = limitinf(la/lb, x)
    if c.is_zero:
        return "<"
    elif c.is_infinite:
        return ">"
    else:
        return "="


class SubsSet(dict):
    """
    Stores (expr, dummy) pairs, and how to rewrite expr-s.

    The gruntz algorithm needs to rewrite certain expressions in term of a new
    variable w. We cannot use subs, because it is just too smart for us. For
    example::

        > Omega=[exp(exp(_p - exp(-_p))/(1 - 1/_p)), exp(exp(_p))]
        > O2=[exp(-exp(_p) + exp(-exp(-_p))*exp(_p)/(1 - 1/_p))/_w, 1/_w]
        > e = exp(exp(_p - exp(-_p))/(1 - 1/_p)) - exp(exp(_p))
        > e.subs(Omega[0],O2[0]).subs(Omega[1],O2[1])
        -1/w + exp(exp(p)*exp(-exp(-p))/(1 - 1/p))

    is really not what we want!

    So we do it the hard way and keep track of all the things we potentially
    want to substitute by dummy variables. Consider the expression::

        exp(x - exp(-x)) + exp(x) + x.

    The mrv set is {exp(x), exp(-x), exp(x - exp(-x))}.
    We introduce corresponding dummy variables d1, d2, d3 and rewrite::

        d3 + d1 + x.

    This class first of all keeps track of the mapping expr->variable, i.e.
    will at this stage be a dictionary::

        {exp(x): d1, exp(-x): d2, exp(x - exp(-x)): d3}.

    [It turns out to be more convenient this way round.]
    But sometimes expressions in the mrv set have other expressions from the
    mrv set as subexpressions, and we need to keep track of that as well. In
    this case, d3 is really exp(x - d2), so rewrites at this stage is::

        {d3: exp(x-d2)}.

    The function rewrite uses all this information to correctly rewrite our
    expression in terms of w. In this case w can be choosen to be exp(-x),
    i.e. d2. The correct rewriting then is::

        exp(-w)/w + 1/w + x.
    """
    def __init__(self):
        self.rewrites = {}

    def __repr__(self):
        return super(SubsSet, self).__repr__() + ', ' + self.rewrites.__repr__()

    def __getitem__(self, key):
        if not key in self:
            self[key] = Dummy()
        return dict.__getitem__(self, key)

    def do_subs(self, e):
        for expr, var in self.items():
            e = e.xreplace({var: expr})
        return e

    def meets(self, s2):
        """Tell whether or not self and s2 have non-empty intersection"""
        return set(self.keys()).intersection(list(s2.keys())) != set()

    def union(self, s2, exps=None):
        """Compute the union of self and s2, adjusting exps"""
        res = self.copy()
        tr = {}
        for expr, var in s2.items():
            if expr in self:
                if exps:
                    exps = exps.xreplace({var: res[expr]})
                tr[var] = res[expr]
            else:
                res[expr] = var
        for var, rewr in s2.rewrites.items():
            res.rewrites[var] = rewr.xreplace(tr)
        return res, exps

    def copy(self):
        r = SubsSet()
        r.rewrites = self.rewrites.copy()
        for expr, var in self.items():
            r[expr] = var
        return r


def mrv(e, x):
    """Returns a SubsSet of most rapidly varying (mrv) subexpressions of 'e',
       and e rewritten in terms of these"""
    e = powsimp(e, deep=True, combine='exp')
    if not e.has(x):
        return SubsSet(), e
    elif e == x:
        s = SubsSet()
        return s, s[x]
    elif e.is_Mul or e.is_Add:
        i, d = e.as_independent(x)  # throw away x-independent terms
        if d.func != e.func:
            s, expr = mrv(d, x)
            return s, e.func(i, expr)
        a, b = d.as_two_terms()
        s1, e1 = mrv(a, x)
        s2, e2 = mrv(b, x)
        return mrv_max1(s1, s2, e.func(i, e1, e2), x)
    elif e.is_Pow:
        b, e = e.as_base_exp()
        if e.has(x):
            return mrv(exp(e * log(b)), x)
        else:
            s, expr = mrv(b, x)
            return s, expr**e
    elif e.func is log:
        s, expr = mrv(e.args[0], x)
        return s, log(expr)
    elif e.func is exp:
        # We know from the theory of this algorithm that exp(log(...)) may always
        # be simplified here, and doing so is vital for termination.
        if e.args[0].func is log:
            return mrv(e.args[0].args[0], x)
        # if a product has an infinite factor the result will be
        # infinite if there is no zero, otherwise NaN; here, we
        # consider the result infinite if any factor is infinite
        li = limitinf(e.args[0], x)
        if any(_.is_infinite for _ in Mul.make_args(li)):
            s1 = SubsSet()
            e1 = s1[e]
            s2, e2 = mrv(e.args[0], x)
            su = s1.union(s2)[0]
            su.rewrites[e1] = exp(e2)
            return mrv_max3(s1, e1, s2, exp(e2), su, e1, x)
        else:
            s, expr = mrv(e.args[0], x)
            return s, exp(expr)
    elif e.is_Function:
        l = [mrv(a, x) for a in e.args]
        l2 = [s for (s, _) in l if s != SubsSet()]
        if len(l2) != 1:
            # e.g. something like BesselJ(x, x)
            raise NotImplementedError("MRV set computation for functions in"
                                      " several variables not implemented.")
        s, ss = l2[0], SubsSet()
        args = [ss.do_subs(x[1]) for x in l]
        return s, e.func(*args)
    elif e.is_Derivative:
        raise NotImplementedError("MRV set computation for derviatives"
                                  " not implemented yet.")
        return mrv(e.args[0], x)
    raise NotImplementedError(
        "Don't know how to calculate the mrv of '%s'" % e)


def mrv_max3(f, expsf, g, expsg, union, expsboth, x):
    """Computes the maximum of two sets of expressions f and g, which
    are in the same comparability class, i.e. max() compares (two elements of)
    f and g and returns either (f, expsf) [if f is larger], (g, expsg)
    [if g is larger] or (union, expsboth) [if f, g are of the same class].
    """
    if not isinstance(f, SubsSet):
        raise TypeError("f should be an instance of SubsSet")
    if not isinstance(g, SubsSet):
        raise TypeError("g should be an instance of SubsSet")
    if f == SubsSet():
        return g, expsg
    elif g == SubsSet():
        return f, expsf
    elif f.meets(g):
        return union, expsboth

    c = compare(list(f.keys())[0], list(g.keys())[0], x)
    if c == ">":
        return f, expsf
    elif c == "<":
        return g, expsg
    else:
        if c != "=":
            raise ValueError("c should be =")
        return union, expsboth


def mrv_max1(f, g, exps, x):
    """Computes the maximum of two sets of expressions f and g, which
    are in the same comparability class, i.e. mrv_max1() compares (two elements of)
    f and g and returns the set, which is in the higher comparability class
    of the union of both, if they have the same order of variation.
    Also returns exps, with the appropriate substitutions made.
    """
    u, b = f.union(g, exps)
    return mrv_max3(f, g.do_subs(exps), g, f.do_subs(exps),
                    u, b, x)


@cacheit
@timeit
def sign(e, x):
    """
    Returns a sign of an expression e(x) for x->oo.

    ::

        e >  0 for x sufficiently large ...  1
        e == 0 for x sufficiently large ...  0
        e <  0 for x sufficiently large ... -1

    The result of this function is currently undefined if e changes sign
    arbitarily often for arbitrarily large x (e.g. sin(x)).

    Note that this returns zero only if e is *constantly* zero
    for x sufficiently large. [If e is constant, of course, this is just
    the same thing as the sign of e.]
    """
    from sympy import sign as _sign

    if e.is_positive:
        return 1
    elif e.is_negative:
        return -1
    elif e.is_zero:
        return 0

    elif not e.has(x):
        return _sign(e)
    elif e == x:
        return 1
    elif e.is_Mul:
        a, b = e.as_two_terms()
        sa = sign(a, x)
        if not sa:
            return 0
        return sa * sign(b, x)
    elif e.func is exp:
        return 1
    elif e.is_Pow:
        s = sign(e.base, x)
        if s == 1:
            return 1
        if e.exp.is_Integer:
            return s**e.exp
    elif e.func is log:
        return sign(e.args[0] - 1, x)

    # if all else fails, do it the hard way
    c0, e0 = mrv_leadterm(e, x)
    return sign(c0, x)


@timeit
@cacheit
def limitinf(e, x):
    """Limit e(x) for x-> oo"""
    # rewrite e in terms of tractable functions only
    e = e.rewrite('tractable', deep=True)

    if not e.has(x):
        return e  # e is a constant
    if e.has(Order):
        e = e.expand().removeO()
    if not x.is_positive:
        # We make sure that x.is_positive is True so we
        # get all the correct mathematical behavior from the expression.
        # We need a fresh variable.
        p = Dummy('p', positive=True, finite=True)
        e = e.subs(x, p)
        x = p
    c0, e0 = mrv_leadterm(e, x)
    sig = sign(e0, x)
    if sig == 1:
        return S.Zero  # e0>0: lim f = 0
    elif sig == -1:  # e0<0: lim f = +-oo (the sign depends on the sign of c0)
        if c0.match(I*Wild("a", exclude=[I])):
            return c0*oo
        s = sign(c0, x)
        # the leading term shouldn't be 0:
        if s == 0:
            raise ValueError("Leading term should not be 0")
        return s*oo
    elif sig == 0:
        return limitinf(c0, x)  # e0=0: lim f = lim c0


def moveup2(s, x):
    r = SubsSet()
    for expr, var in s.items():
        r[expr.xreplace({x: exp(x)})] = var
    for var, expr in s.rewrites.items():
        r.rewrites[var] = s.rewrites[var].xreplace({x: exp(x)})
    return r


def moveup(l, x):
    return [e.xreplace({x: exp(x)}) for e in l]


@timeit
def calculate_series(e, x, logx=None):
    """ Calculates at least one term of the series of "e" in "x".

    This is a place that fails most often, so it is in its own function.
    """
    from sympy.polys import cancel

    for t in e.lseries(x, logx=logx):
        t = cancel(t)

        if t.simplify():
            break

    return t


@timeit
@cacheit
def mrv_leadterm(e, x):
    """Returns (c0, e0) for e."""
    Omega = SubsSet()
    if not e.has(x):
        return (e, S.Zero)
    if Omega == SubsSet():
        Omega, exps = mrv(e, x)
    if not Omega:
        # e really does not depend on x after simplification
        series = calculate_series(e, x)
        c0, e0 = series.leadterm(x)
        if e0 != 0:
            raise ValueError("e0 should be 0")
        return c0, e0
    if x in Omega:
        # move the whole omega up (exponentiate each term):
        Omega_up = moveup2(Omega, x)
        e_up = moveup([e], x)[0]
        exps_up = moveup([exps], x)[0]
        # NOTE: there is no need to move this down!
        e = e_up
        Omega = Omega_up
        exps = exps_up
    #
    # The positive dummy, w, is used here so log(w*2) etc. will expand;
    # a unique dummy is needed in this algorithm
    #
    # For limits of complex functions, the algorithm would have to be
    # improved, or just find limits of Re and Im components separately.
    #
    w = Dummy("w", extended_real=True, positive=True, finite=True)
    f, logw = rewrite(exps, Omega, x, w)
    series = calculate_series(f, w, logx=logw)
    return series.leadterm(w)


def build_expression_tree(Omega, rewrites):
    r""" Helper function for rewrite.

    We need to sort Omega (mrv set) so that we replace an expression before
    we replace any expression in terms of which it has to be rewritten::

        e1 ---> e2 ---> e3
                 \
                  -> e4

    Here we can do e1, e2, e3, e4 or e1, e2, e4, e3.
    To do this we assemble the nodes into a tree, and sort them by height.

    This function builds the tree, rewrites then sorts the nodes.
    """
    class Node:
        def ht(self):
            return reduce(lambda x, y: x + y,
                          [x.ht() for x in self.before], 1)
    nodes = {}
    for expr, v in Omega:
        n = Node()
        n.before = []
        n.var = v
        n.expr = expr
        nodes[v] = n
    for _, v in Omega:
        if v in rewrites:
            n = nodes[v]
            r = rewrites[v]
            for _, v2 in Omega:
                if r.has(v2):
                    n.before.append(nodes[v2])

    return nodes


@timeit
def rewrite(e, Omega, x, wsym):
    """e(x) ... the function
    Omega ... the mrv set
    wsym ... the symbol which is going to be used for w

    Returns the rewritten e in terms of w and log(w). See test_rewrite1()
    for examples and correct results.
    """
    from sympy import ilcm
    if not isinstance(Omega, SubsSet):
        raise TypeError("Omega should be an instance of SubsSet")
    if len(Omega) == 0:
        raise ValueError("Length can not be 0")
    # all items in Omega must be exponentials
    for t in Omega.keys():
        if t.func is not exp:
            raise ValueError("Value should be exp")
    rewrites = Omega.rewrites
    Omega = sorted(list(Omega.items()), key=default_sort_key)

    nodes = build_expression_tree(Omega, rewrites)
    Omega.sort(key=lambda x: nodes[x[1]].ht(), reverse=True)

    # make sure we know the sign of each exp() term; after the loop,
    # g is going to be the "w" - the simplest one in the mrv set
    for g, _ in Omega:
        sig = sign(g.args[0], x)
        if sig != 1 and sig != -1:
            raise NotImplementedError('Result depends on the sign of %s' % sig)
    if sig == 1:
        wsym = 1/wsym  # if g goes to oo, substitute 1/w
    # O2 is a list, which results by rewriting each item in Omega using "w"
    O2 = []
    denominators = []
    for f, var in Omega:
        c = limitinf(f.args[0]/g.args[0], x)
        if c.is_Rational:
            denominators.append(c.q)
        arg = f.args[0]
        if var in rewrites:
            if not rewrites[var].func is exp:
                raise ValueError("Value should be exp")
            arg = rewrites[var].args[0]
        O2.append((var, exp((arg - c*g.args[0]).expand())*wsym**c))

    # Remember that Omega contains subexpressions of "e". So now we find
    # them in "e" and substitute them for our rewriting, stored in O2

    # the following powsimp is necessary to automatically combine exponentials,
    # so that the .xreplace() below succeeds:
    # TODO this should not be necessary
    f = powsimp(e, deep=True, combine='exp')
    for a, b in O2:
        f = f.xreplace({a: b})

    for _, var in Omega:
        assert not f.has(var)

    # finally compute the logarithm of w (logw).
    logw = g.args[0]
    if sig == 1:
        logw = -logw  # log(w)->log(1/w)=-log(w)

    # Some parts of sympy have difficulty computing series expansions with
    # non-integral exponents. The following heuristic improves the situation:
    exponent = reduce(ilcm, denominators, 1)
    f = f.xreplace({wsym: wsym**exponent})
    logw /= exponent

    return f, logw


def gruntz(e, z, z0, dir="+"):
    """
    Compute the limit of e(z) at the point z0 using the Gruntz algorithm.

    z0 can be any expression, including oo and -oo.

    For dir="+" (default) it calculates the limit from the right
    (z->z0+) and for dir="-" the limit from the left (z->z0-). For infinite z0
    (oo or -oo), the dir argument doesn't matter.

    This algorithm is fully described in the module docstring in the gruntz.py
    file. It relies heavily on the series expansion. Most frequently, gruntz()
    is only used if the faster limit() function (which uses heuristics) fails.
    """
    if not isinstance(z, Symbol):
        raise NotImplementedError("Second argument must be a Symbol")

    # convert all limits to the limit z->oo; sign of z is handled in limitinf
    r = None
    if z0 == oo:
        r = limitinf(e, z)
    elif z0 == -oo:
        r = limitinf(e.subs(z, -z), z)
    else:
        if str(dir) == "-":
            e0 = e.subs(z, z0 - 1/z)
        elif str(dir) == "+":
            e0 = e.subs(z, z0 + 1/z)
        else:
            raise NotImplementedError("dir must be '+' or '-'")
        r = limitinf(e0, z)

    # This is a bit of a heuristic for nice results... we always rewrite
    # tractable functions in terms of familiar intractable ones.
    # It might be nicer to rewrite the exactly to what they were initially,
    # but that would take some work to implement.
    return r.rewrite('intractable', deep=True)
