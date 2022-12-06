"""
Implementation of the trigsimp algorithm by Fu et al.

The idea behind the ``fu`` algorithm is to use a sequence of rules, applied
in what is heuristically known to be a smart order, to select a simpler
expression that is equivalent to the input.

There are transform rules in which a single rule is applied to the
expression tree. The following are just mnemonic in nature; see the
docstrings for examples.

    TR0 - simplify expression
    TR1 - sec-csc to cos-sin
    TR2 - tan-cot to sin-cos ratio
    TR2i - sin-cos ratio to tan
    TR3 - angle canonicalization
    TR4 - functions at special angles
    TR5 - powers of sin to powers of cos
    TR6 - powers of cos to powers of sin
    TR7 - reduce cos power (increase angle)
    TR8 - expand products of sin-cos to sums
    TR9 - contract sums of sin-cos to products
    TR10 - separate sin-cos arguments
    TR10i - collect sin-cos arguments
    TR11 - reduce double angles
    TR12 - separate tan arguments
    TR12i - collect tan arguments
    TR13 - expand product of tan-cot
    TRmorrie - prod(cos(x*2**i), (i, 0, k - 1)) -> sin(2**k*x)/(2**k*sin(x))
    TR14 - factored powers of sin or cos to cos or sin power
    TR15 - negative powers of sin to cot power
    TR16 - negative powers of cos to tan power
    TR22 - tan-cot powers to negative powers of sec-csc functions
    TR111 - negative sin-cos-tan powers to csc-sec-cot

There are 4 combination transforms (CTR1 - CTR4) in which a sequence of
transformations are applied and the simplest expression is selected from
a few options.

Finally, there are the 2 rule lists (RL1 and RL2), which apply a
sequence of transformations and combined transformations, and the ``fu``
algorithm itself, which applies rules and rule lists and selects the
best expressions. There is also a function ``L`` which counts the number
of trigonometric functions that appear in the expression.

Other than TR0, re-writing of expressions is not done by the transformations.
e.g. TR10i finds pairs of terms in a sum that are in the form like
``cos(x)*cos(y) + sin(x)*sin(y)``. Such expression are targeted in a bottom-up
traversal of the expression, but no manipulation to make them appear is
attempted. For example,

    Set-up for examples below:

    >>> from time import time

>>> eq = cos(x + y)/cos(x)
>>> TR10i(eq.expand(trig=True))
-sin(x)*sin(y)/cos(x) + cos(y)

If the expression is put in "normal" form (with a common denominator) then
the transformation is successful:

>>> TR10i(_.normal())
cos(x + y)/cos(x)

TR11's behavior is similar. It rewrites double angles as smaller angles but
doesn't do any simplification of the result.

>>> TR11(sin(2)**a*cos(1)**(-a), 1)
(2*sin(1)*cos(1))**a*cos(1)**(-a)
>>> powsimp(_)
(2*sin(1))**a

The temptation is to try make these TR rules "smarter" but that should really
be done at a higher level; the TR rules should try maintain the "do one thing
well" principle.  There is one exception, however. In TR10i and TR9 terms are
recognized even when they are each multiplied by a common factor:

>>> fu(a*cos(x)*cos(y) + a*sin(x)*sin(y))
a*cos(x - y)

Factoring with ``factor_terms`` is used but it it "JIT"-like, being delayed
until it is deemed necessary. Furthermore, if the factoring does not
help with the simplification, it is not retained, so
``a*cos(x)*cos(y) + a*sin(x)*sin(z)`` does not become the factored
(but unsimplified in the trigonometric sense) expression:

>>> fu(a*cos(x)*cos(y) + a*sin(x)*sin(z))
a*sin(x)*sin(z) + a*cos(x)*cos(y)

In some cases factoring might be a good idea, but the user is left
to make that decision. For example:

>>> expr = ((15*sin(2*x) + 19*sin(x + y) + 17*sin(x + z) + 19*cos(x - z) +
...          25)*(20*sin(2*x) + 15*sin(x + y) + sin(y + z) + 14*cos(x - z) +
...               14*cos(y - z))*(9*sin(2*y) + 12*sin(y + z) +
...                               10*cos(x - y) + 2*cos(y - z) +
...                               18)).expand(trig=True).expand()

In the expanded state, there are nearly 1000 trig functions:

>>> L(expr)
932

If the expression where factored first, this would take time but the
resulting expression would be transformed very quickly:

>>> def clock(f, n=2):
...     t = time()
...     f()
...     return round(time() - t, n)
...
>>> clock(lambda: factor(expr))  # doctest: +SKIP
0.86
>>> clock(lambda: TR10i(expr), 3)  # doctest: +SKIP
0.016

If the unexpanded expression is used, the transformation takes longer but
not as long as it took to factor it and then transform it:

>>> clock(lambda: TR10i(expr), 2)  # doctest: +SKIP
0.28

So neither expansion nor factoring is used in ``TR10i``: if the
expression is already factored (or partially factored) then expansion
with ``trig=True`` would destroy what is already known and take
longer; if the expression is expanded, factoring may take longer than
simply applying the transformation itself.

Although the algorithms should be canonical, always giving the same
result, they may not yield the best result. This, in general, is
the nature of simplification where searching all possible transformation
paths is very expensive. Here is a simple example. There are 6 terms
in the following sum:

>>> expr = (sin(x)**2*cos(y)*cos(z) + sin(x)*sin(y)*cos(x)*cos(z) +
...         sin(x)*sin(z)*cos(x)*cos(y) + sin(y)*sin(z)*cos(x)**2 +
...         sin(y)*sin(z) + cos(y)*cos(z))
>>> args = expr.args

Serendipitously, fu gives the best result:

>>> fu(expr)
3*cos(y - z)/2 - cos(2*x + y + z)/2

But if different terms were combined, a less-optimal result might be
obtained, requiring some additional work to get better simplification,
but still less than optimal. The following shows an alternative form
of ``expr`` that resists optimal simplification once a given step
is taken since it leads to a dead end:

>>> TR9(-cos(x)**2*cos(y + z) + 3*cos(y - z)/2 +
...     cos(y + z)/2 + cos(-2*x + y + z)/4 - cos(2*x + y + z)/4)
sin(2*x)*sin(y + z)/2 - cos(x)**2*cos(y + z) + 3*cos(y - z)/2 + cos(y + z)/2

Here is a smaller expression that exhibits the same behavior:

>>> a = sin(x)*sin(z)*cos(x)*cos(y) + sin(x)*sin(y)*cos(x)*cos(z)
>>> TR10i(a)
sin(x)*sin(y + z)*cos(x)
>>> newa = _
>>> TR10i(expr - a)  # this combines two more of the remaining terms
sin(x)**2*cos(y)*cos(z) + sin(y)*sin(z)*cos(x)**2 + cos(y - z)
>>> TR10i(_ + newa) == _ + newa  # but now there is no more simplification
True

Without getting lucky or trying all possible pairings of arguments, the
final result may be less than optimal and impossible to find without
better heuristics or brute force trial of all possibilities.

Notes
=====

This work was started by Dimitar Vlahovski at the Technological School
"Electronic systems" (30.11.2011).

References
==========

Fu, Hongguang, Xiuqin Zhong, and Zhenbing Zeng. "Automated and readable
simplification of trigonometric expressions." Mathematical and computer
modelling 44.11 (2006): 1169-1177.
http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/DESTIME2006/DES_contribs/Fu/simplification.pdf

http://www.sosmath.com/trig/Trig5/trig5/pdf/pdf.html gives a formula sheet.

"""

from collections import defaultdict

from ..core import (Add, Dummy, Expr, I, Integer, Mul, Pow, Rational,
                    expand_mul, factor_terms, gcd_terms, pi)
from ..core.exprtools import Factors
from ..core.strategies import greedy, identity
from ..core.sympify import sympify
from ..functions import (binomial, cos, cosh, cot, coth, csc, sec, sin, sinh,
                         sqrt, tan, tanh)
from ..functions.elementary.hyperbolic import HyperbolicFunction
from ..functions.elementary.trigonometric import TrigonometricFunction
from ..ntheory import perfect_power
from ..polys.polytools import factor
from ..utilities import ordered
from .simplify import bottom_up


# ================== Fu-like tools ===========================


def TR0(rv):
    """Simplification of rational polynomials, trying to simplify
    the expression, e.g. combine things like 3*x + 2*x, etc....

    """
    # although it would be nice to use cancel, it doesn't work
    # with noncommutatives
    return rv.normal().factor().expand()


def TR1(rv):
    """Replace sec, csc with 1/cos, 1/sin

    Examples
    ========

    >>> TR1(2*csc(x) + sec(x))
    1/cos(x) + 2/sin(x)

    """

    def f(rv):
        if isinstance(rv, sec):
            a = rv.args[0]
            return 1/cos(a)
        elif isinstance(rv, csc):
            a = rv.args[0]
            return 1/sin(a)
        return rv

    return bottom_up(rv, f)


def TR2(rv):
    """Replace tan and cot with sin/cos and cos/sin

    Examples
    ========

    >>> TR2(tan(x))
    sin(x)/cos(x)
    >>> TR2(cot(x))
    cos(x)/sin(x)
    >>> TR2(tan(tan(x) - sin(x)/cos(x)))
    0

    """

    def f(rv):
        if isinstance(rv, tan):
            a = rv.args[0]
            return sin(a)/cos(a)
        elif isinstance(rv, cot):
            a = rv.args[0]
            return cos(a)/sin(a)
        return rv

    return bottom_up(rv, f)


def TR2i(rv, half=False):
    """Converts ratios involving sin and cos as follows::
        sin(x)/cos(x) -> tan(x)
        sin(x)/(cos(x) + 1) -> tan(x/2) if half=True

    Examples
    ========

    >>> TR2i(sin(x)/cos(x))
    tan(x)

    Powers of the numerator and denominator are also recognized

    >>> TR2i(sin(x)**2/(cos(x) + 1)**2, half=True)
    tan(x/2)**2

    The transformation does not take place unless assumptions allow
    (i.e. the base must be positive or the exponent must be an integer
    for both numerator and denominator)

    >>> TR2i(sin(x)**a/(cos(x) + 1)**a)
    (cos(x) + 1)**(-a)*sin(x)**a

    """

    def f(rv):
        if not rv.is_Mul:
            return rv

        n, d = rv.as_numer_denom()
        if n.is_Atom or d.is_Atom:
            return rv

        def ok(k, e):
            # initial filtering of factors
            return (
                (e.is_integer or k.is_positive) and (
                    k.func in (sin, cos) or (half and
                                             k.is_Add and
                                             len(k.args) >= 2 and
                                             any(any(isinstance(ai, cos) or ai.is_Pow and ai.base is cos
                                                     for ai in Mul.make_args(a)) for a in k.args))))

        n = n.as_powers_dict()
        ndone = [(k, n.pop(k)) for k in list(n) if not ok(k, n[k])]
        if not n:
            return rv

        d = d.as_powers_dict()
        ddone = [(k, d.pop(k)) for k in list(d) if not ok(k, d[k])]
        if not d:
            return rv

        # factoring if necessary

        def factorize(d, ddone):
            newk = []
            for k in d:
                if k.is_Add and len(k.args) > 1:
                    knew = factor(k) if half else factor_terms(k)
                    if knew != k:
                        newk.append((k, knew))
            if newk:
                for i, (k, knew) in enumerate(newk):
                    del d[k]
                    newk[i] = knew
                newk = Mul(*newk).as_powers_dict()
                for k in newk:
                    v = d[k] + newk[k]
                    if ok(k, v):
                        d[k] = v
                    else:
                        ddone.append((k, v))
        factorize(n, ndone)
        factorize(d, ddone)

        # joining
        t = []
        for k in n:
            if isinstance(k, sin):
                a = cos(k.args[0], evaluate=False)
                if a in d and d[a] == n[k]:
                    t.append(tan(k.args[0])**n[k])
                    n[k] = d[a] = None
                elif half:
                    a1 = 1 + a
                    if a1 in d and d[a1] == n[k]:
                        t.append((tan(k.args[0]/2))**n[k])
                        n[k] = d[a1] = None
            elif isinstance(k, cos):
                a = sin(k.args[0], evaluate=False)
                if a in d and d[a] == n[k]:
                    t.append(tan(k.args[0])**-n[k])
                    n[k] = d[a] = None
            elif (half and k.is_Add and k.args[0] == 1 and
                  isinstance(k.args[1], cos)):
                a = sin(k.args[1].args[0], evaluate=False)
                if a in d and d[a] == n[k] and (d[a].is_integer or
                                                a.is_positive):
                    t.append(tan(a.args[0]/2)**-n[k])
                    n[k] = d[a] = None

        if t:
            rv = Mul(*(t + [b**e for b, e in n.items() if e])) /\
                Mul(*[b**e for b, e in d.items() if e])
            rv *= Mul(*[b**e for b, e in ndone])/Mul(*[b**e for b, e in ddone])

        return rv

    return bottom_up(rv, f)


def TR3(rv):
    """Induced formula: example sin(-a) = -sin(a)

    Examples
    ========

    >>> TR3(cos(y - x*(y - x)))
    cos(x*(x - y) + y)
    >>> cos(pi/2 + x)
    -sin(x)
    >>> cos(30*pi/2 + x)
    -cos(x)

    """
    from .simplify import signsimp

    # Negative argument (already automatic for funcs like sin(-x) -> -sin(x)
    # but more complicated expressions can use it, too). Also, trig angles
    # between pi/4 and pi/2 are not reduced to an angle between 0 and pi/4.
    # The following are automatically handled:
    #   Argument of type: pi/2 +/- angle
    #   Argument of type: pi +/- angle
    #   Argument of type : 2k*pi +/- angle

    def f(rv):
        if not isinstance(rv, TrigonometricFunction):
            return rv
        rv = rv.func(signsimp(rv.args[0]))
        if (rv.args[0] - pi/4).is_positive is (pi/2 - rv.args[0]).is_positive is True:
            fmap = {cos: sin, sin: cos, tan: cot, cot: tan, sec: csc, csc: sec}
            rv = fmap[rv.func](pi/2 - rv.args[0])
        return rv

    return bottom_up(rv, f)


def TR4(rv):
    """Identify values of special angles.

        a=  0   pi/6        pi/4        pi/3        pi/2
    ----------------------------------------------------
    cos(a)  0   1/2         sqrt(2)/2   sqrt(3)/2   1
    sin(a)  1   sqrt(3)/2   sqrt(2)/2   1/2         0
    tan(a)  0   sqt(3)/3    1           sqrt(3)     --

    Examples
    ========

    >>> for s in (0, pi/6, pi/4, pi/3, pi/2):
    ...     print(f'{cos(s)} {sin(s)} {tan(s)} {cot(s)}')
    ...
    1 0 0 zoo
    sqrt(3)/2 1/2 sqrt(3)/3 sqrt(3)
    sqrt(2)/2 sqrt(2)/2 1 1
    1/2 sqrt(3)/2 sqrt(3) sqrt(3)/3
    0 1 zoo 0

    """
    # special values at 0, pi/6, pi/4, pi/3, pi/2 already handled
    return rv


def _TR56(rv, f, g, h, max, pow):
    """Helper for TR5 and TR6 to replace f**2 with h(g**2)

    Options
    =======

    max :   controls size of exponent that can appear on f
            e.g. if max=4 then f**4 will be changed to h(g**2)**2.
    pow :   controls whether the exponent must be a perfect power of 2
            e.g. if pow=True (and max >= 6) then f**6 will not be changed
            but f**8 will be changed to h(g**2)**4

    >>> def h(x):
    ...     return 1 - x
    >>> _TR56(sin(x)**3, sin, cos, h, 4, False)
    sin(x)**3
    >>> _TR56(sin(x)**6, sin, cos, h, 6, False)
    (-cos(x)**2 + 1)**3
    >>> _TR56(sin(x)**6, sin, cos, h, 6, True)
    sin(x)**6
    >>> _TR56(sin(x)**8, sin, cos, h, 10, True)
    (-cos(x)**2 + 1)**4

    """

    def _f(rv):
        # I'm not sure if this transformation should target all even powers
        # or only those expressible as powers of 2. Also, should it only
        # make the changes in powers that appear in sums -- making an isolated
        # change is not going to allow a simplification as far as I can tell.
        if not (rv.is_Pow and rv.base.func == f):
            return rv

        if rv.exp.is_negative:
            return rv
        if (rv.exp - max).is_positive:
            return rv
        if rv.exp == 2:
            return h(g(rv.base.args[0])**2)
        else:
            if rv.exp == 4:
                e = 2
            elif not pow:
                if rv.exp % 2:
                    return rv
                e = rv.exp//2
            else:
                p = perfect_power(rv.exp)
                if not p:
                    return rv
                e = rv.exp//2
            return h(g(rv.base.args[0])**2)**e

    return bottom_up(rv, _f)


def TR5(rv, max=4, pow=False):
    """Replacement of sin**2 with 1 - cos(x)**2.

    See _TR56 docstring for advanced use of ``max`` and ``pow``.

    Examples
    ========

    >>> TR5(sin(x)**2)
    -cos(x)**2 + 1
    >>> TR5(sin(x)**-2)  # unchanged
    sin(x)**(-2)
    >>> TR5(sin(x)**4)
    (-cos(x)**2 + 1)**2

    """
    return _TR56(rv, sin, cos, lambda x: 1 - x, max=max, pow=pow)


def TR6(rv, max=4, pow=False):
    """Replacement of cos**2 with 1 - sin(x)**2.

    See _TR56 docstring for advanced use of ``max`` and ``pow``.

    Examples
    ========

    >>> TR6(cos(x)**2)
    -sin(x)**2 + 1
    >>> TR6(cos(x)**-2)  # unchanged
    cos(x)**(-2)
    >>> TR6(cos(x)**4)
    (-sin(x)**2 + 1)**2

    """
    return _TR56(rv, cos, sin, lambda x: 1 - x, max=max, pow=pow)


def TR7(rv):
    """Lowering the degree of cos(x)**2

    Examples
    ========

    >>> TR7(cos(x)**2)
    cos(2*x)/2 + 1/2
    >>> TR7(cos(x)**2 + 1)
    cos(2*x)/2 + 3/2

    """

    def f(rv):
        if not (rv.is_Pow and rv.base.func == cos and rv.exp == 2):
            return rv
        return (1 + cos(2*rv.base.args[0]))/2

    return bottom_up(rv, f)


def TR8(rv, first=True):
    """Converting products of ``cos`` and/or ``sin`` to a sum or
    difference of ``cos`` and or ``sin`` terms.

    Examples
    ========

    >>> TR8(cos(2)*cos(3))
    cos(5)/2 + cos(1)/2
    >>> TR8(cos(2)*sin(3))
    sin(5)/2 + sin(1)/2
    >>> TR8(sin(2)*sin(3))
    -cos(5)/2 + cos(1)/2

    """

    def f(rv):
        if not (rv.is_Mul or rv.is_Pow and
                rv.base.func in (cos, sin) and
                (rv.exp.is_integer or rv.base.is_positive)):
            return rv

        if first:
            n, d = map(expand_mul, rv.as_numer_denom())
            newn = TR8(n, first=False)
            newd = TR8(d, first=False)
            if newn != n or newd != d:
                rv = gcd_terms(newn/newd)
                if rv.is_Mul and rv.args[0].is_Rational and \
                        len(rv.args) == 2 and rv.args[1].is_Add:
                    rv = Mul(*rv.as_coeff_Mul())
            return rv

        args = {cos: [], sin: [], None: []}
        for a in ordered(Mul.make_args(rv)):
            if a.func in (cos, sin):
                args[a.func].append(a.args[0])
            elif (a.is_Pow and a.exp.is_Integer and a.exp > 0 and
                    a.base.func in (cos, sin)):
                # XXX this is ok but pathological expression could be handled
                # more efficiently as in TRmorrie
                args[a.base.func].extend([a.base.args[0]]*a.exp)
            else:
                args[None].append(a)
        c = args[cos]
        s = args[sin]
        if not (c and s or len(c) > 1 or len(s) > 1):
            return rv

        args = args[None]
        n = min(len(c), len(s))
        for _ in range(n):
            a1 = s.pop()
            a2 = c.pop()
            args.append((sin(a1 + a2) + sin(a1 - a2))/2)
        while len(c) > 1:
            a1 = c.pop()
            a2 = c.pop()
            args.append((cos(a1 + a2) + cos(a1 - a2))/2)
        if c:
            args.append(cos(c.pop()))
        while len(s) > 1:
            a1 = s.pop()
            a2 = s.pop()
            args.append((-cos(a1 + a2) + cos(a1 - a2))/2)
        if s:
            args.append(sin(s.pop()))
        return TR8(expand_mul(Mul(*args)))

    return bottom_up(rv, f)


def TR9(rv):
    """Sum of ``cos`` or ``sin`` terms as a product of ``cos`` or ``sin``.

    Examples
    ========

    >>> TR9(cos(1) + cos(2))
    2*cos(1/2)*cos(3/2)
    >>> TR9(cos(1) + 2*sin(1) + 2*sin(2))
    cos(1) + 4*sin(3/2)*cos(1/2)

    If no change is made by TR9, no re-arrangement of the
    expression will be made. For example, though factoring
    of common term is attempted, if the factored expression
    wasn't changed, the original expression will be returned:

    >>> TR9(cos(3) + cos(3)*cos(2))
    cos(3) + cos(2)*cos(3)

    """

    def f(rv):
        if not rv.is_Add:
            return rv

        def do(rv, first=True):
            # cos(a)+/-cos(b) can be combined into a product of cosines and
            # sin(a)+/-sin(b) can be combined into a product of cosine and
            # sine.
            #
            # If there are more than two args, the pairs which "work" will
            # have a gcd extractable and the remaining two terms will have
            # the above structure -- all pairs must be checked to find the
            # ones that work. args that don't have a common set of symbols
            # are skipped since this doesn't lead to a simpler formula and
            # also has the arbitrariness of combining, for example, the x
            # and y term instead of the y and z term in something like
            # cos(x) + cos(y) + cos(z).

            if not rv.is_Add:
                return rv

            args = list(ordered(rv.args))
            if len(args) != 2:
                hit = False
                for i, ai in enumerate(args):
                    if ai is None:
                        continue
                    for j in range(i + 1, len(args)):
                        aj = args[j]
                        if aj is None:
                            continue
                        was = ai + aj
                        new = do(was)
                        if new != was:
                            args[i] = new  # update in place
                            args[j] = None
                            hit = True
                            break  # go to next i
                if hit:
                    rv = Add(*[_f for _f in args if _f])
                    if rv.is_Add:
                        rv = do(rv)

                return rv

            # two-arg Add
            split = trig_split(*args)
            if not split:
                return rv
            gcd, n1, n2, a, b, iscos = split

            # application of rule if possible
            if iscos:
                if n1 == n2:
                    return gcd*n1*2*cos((a + b)/2)*cos((a - b)/2)
                if n1 < 0:
                    a, b = b, a
                return -2*gcd*sin((a + b)/2)*sin((a - b)/2)
            else:
                if n1 == n2:
                    return gcd*n1*2*sin((a + b)/2)*cos((a - b)/2)
                if n1 < 0:
                    a, b = b, a
                return 2*gcd*cos((a + b)/2)*sin((a - b)/2)

        return process_common_addends(rv, do)  # DON'T sift by free symbols

    return bottom_up(rv, f)


def TR10(rv, first=True):
    """Separate sums in ``cos`` and ``sin``.

    Examples
    ========

    >>> TR10(cos(a + b))
    -sin(a)*sin(b) + cos(a)*cos(b)
    >>> TR10(sin(a + b))
    sin(a)*cos(b) + sin(b)*cos(a)
    >>> TR10(sin(a + b + c))
    (-sin(a)*sin(b) + cos(a)*cos(b))*sin(c) +
    (sin(a)*cos(b) + sin(b)*cos(a))*cos(c)

    """

    def f(rv):
        if rv.func not in (cos, sin):
            return rv

        f = rv.func
        arg = rv.args[0]
        if arg.is_Add:
            if first:
                args = list(ordered(arg.args))
            else:
                args = list(arg.args)
            a = args.pop()
            b = Add._from_args(args)
            if b.is_Add:
                if f == sin:
                    return sin(a)*TR10(cos(b), first=False) + \
                        cos(a)*TR10(sin(b), first=False)
                else:
                    return cos(a)*TR10(cos(b), first=False) - \
                        sin(a)*TR10(sin(b), first=False)
            else:
                if f == sin:
                    return sin(a)*cos(b) + cos(a)*sin(b)
                else:
                    return cos(a)*cos(b) - sin(a)*sin(b)
        return rv

    return bottom_up(rv, f)


def TR10i(rv):
    """Sum of products to function of sum.

    Examples
    ========

    >>> TR10i(cos(1)*cos(3) + sin(1)*sin(3))
    cos(2)
    >>> TR10i(cos(1)*sin(3) + sin(1)*cos(3) + cos(3))
    cos(3) + sin(4)
    >>> TR10i(sqrt(2)*cos(x)*x + sqrt(6)*sin(x)*x)
    2*sqrt(2)*x*sin(x + pi/6)

    """
    global _ROOT2, _ROOT3, _invROOT3
    if _ROOT2 is None:
        _roots()

    def f(rv):
        if not rv.is_Add:
            return rv

        def do(rv, first=True):
            # args which can be expressed as A*(cos(a)*cos(b)+/-sin(a)*sin(b))
            # or B*(cos(a)*sin(b)+/-cos(b)*sin(a)) can be combined into
            # A*f(a+/-b) where f is either sin or cos.
            #
            # If there are more than two args, the pairs which "work" will have
            # a gcd extractable and the remaining two terms will have the above
            # structure -- all pairs must be checked to find the ones that
            # work.

            if not rv.is_Add:
                return rv

            args = list(ordered(rv.args))
            if len(args) != 2:
                hit = False
                for i, ai in enumerate(args):
                    if ai is None:
                        continue
                    for j in range(i + 1, len(args)):
                        aj = args[j]
                        if aj is None:
                            continue
                        was = ai + aj
                        new = do(was)
                        if new != was:
                            args[i] = new  # update in place
                            args[j] = None
                            hit = True
                            break  # go to next i
                if hit:
                    rv = Add(*[_f for _f in args if _f])
                    if rv.is_Add:
                        rv = do(rv)

                return rv

            # two-arg Add
            split = trig_split(*args, two=True)
            if not split:
                return rv
            gcd, n1, n2, a, b, same = split

            # identify and get c1 to be cos then apply rule if possible
            if same:  # coscos, sinsin
                gcd = n1*gcd
                if n1 == n2:
                    return gcd*cos(a - b)
                return gcd*cos(a + b)
            else:  # cossin, cossin
                gcd = n1*gcd
                if n1 == n2:
                    return gcd*sin(a + b)
                return gcd*sin(b - a)

        rv = process_common_addends(
            rv, do, lambda x: tuple(ordered(x.free_symbols)))

        # need to check for inducible pairs in ratio of sqrt(3):1 that
        # appeared in different lists when sorting by coefficient
        while rv.is_Add:
            byrad = defaultdict(list)
            for a in rv.args:
                hit = 0
                if a.is_Mul:
                    for ai in a.args:
                        if ai.is_Pow and ai.exp == Rational(1, 2) and \
                                ai.base.is_Integer:
                            byrad[ai].append(a)
                            hit = 1
                            break
                if not hit:
                    byrad[Integer(1)].append(a)

            # no need to check all pairs -- just check for the onees
            # that have the right ratio
            args = []
            for a in byrad:
                for b in [_ROOT3*a, _invROOT3]:
                    if b in byrad:
                        for i in range(len(byrad[a])):
                            if byrad[a][i] is None:
                                continue
                            for j in range(len(byrad[b])):
                                if byrad[b][j] is None:
                                    continue
                                was = Add(byrad[a][i] + byrad[b][j])
                                new = do(was)
                                if new != was:
                                    args.append(new)
                                    byrad[a][i] = None
                                    byrad[b][j] = None
                                    break
            if args:
                rv = Add(*(args + [Add(*[_f for _f in v if _f])
                                   for v in byrad.values()]))
            else:
                rv = do(rv)  # final pass to resolve any new inducible pairs
                break

        return rv

    return bottom_up(rv, f)


def TR11(rv, base=None):
    """Function of double angle to product. The ``base`` argument can be used
    to indicate what is the un-doubled argument, e.g. if 3*pi/7 is the base
    then cosine and sine functions with argument 6*pi/7 will be replaced.

    Examples
    ========

    >>> TR11(sin(2*x))
    2*sin(x)*cos(x)
    >>> TR11(cos(2*x))
    -sin(x)**2 + cos(x)**2
    >>> TR11(sin(4*x))
    4*(-sin(x)**2 + cos(x)**2)*sin(x)*cos(x)
    >>> TR11(sin(4*x/3))
    4*(-sin(x/3)**2 + cos(x/3)**2)*sin(x/3)*cos(x/3)

    If the arguments are simply integers, no change is made
    unless a base is provided:

    >>> TR11(cos(2))
    cos(2)
    >>> TR11(cos(4), 2)
    -sin(2)**2 + cos(2)**2

    There is a subtle issue here in that autosimplification will convert
    some higher angles to lower angles

    >>> cos(6*pi/7) + cos(3*pi/7)
    -cos(pi/7) + cos(3*pi/7)

    The 6*pi/7 angle is now pi/7 but can be targeted with TR11 by supplying
    the 3*pi/7 base:

    >>> TR11(_, 3*pi/7)
    -sin(3*pi/7)**2 + cos(3*pi/7)**2 + cos(3*pi/7)

    """

    def f(rv):
        if rv.func not in (cos, sin):
            return rv

        if base:
            f = rv.func
            t = f(base*2)
            co = Integer(1)
            if t.is_Mul:
                co, t = t.as_coeff_Mul()
            if t.func not in (cos, sin):
                return rv
            if rv.args[0] == t.args[0]:
                c = cos(base)
                s = sin(base)
                if f is cos:
                    return (c**2 - s**2)/co
                else:
                    return 2*c*s/co
            return rv

        elif not rv.args[0].is_Number:
            # make a change if the leading coefficient's numerator is
            # divisible by 2
            c, m = rv.args[0].as_coeff_Mul(rational=True)
            if c.numerator % 2 == 0:
                arg = c.numerator//2*m/c.denominator
                c = TR11(cos(arg))
                s = TR11(sin(arg))
                if rv.func == sin:
                    rv = 2*s*c
                else:
                    rv = c**2 - s**2
        return rv

    return bottom_up(rv, f)


def TR12(rv, first=True):
    """Separate sums in ``tan``.

    Examples
    ========

    >>> TR12(tan(x + y))
    (tan(x) + tan(y))/(-tan(x)*tan(y) + 1)

    """

    def f(rv):
        if not rv.func == tan:
            return rv

        arg = rv.args[0]
        if arg.is_Add:
            if first:
                args = list(ordered(arg.args))
            else:
                args = list(arg.args)
            a = args.pop()
            b = Add._from_args(args)
            if b.is_Add:
                tb = TR12(tan(b), first=False)
            else:
                tb = tan(b)
            return (tan(a) + tb)/(1 - tan(a)*tb)
        return rv

    return bottom_up(rv, f)


def TR12i(rv):
    """Combine tan arguments as
    (tan(y) + tan(x))/(tan(x)*tan(y) - 1) -> -tan(x + y)

    Examples
    ========

    >>> ta, tb, tc = [tan(i) for i in (a, b, c)]
    >>> TR12i((ta + tb)/(-ta*tb + 1))
    tan(a + b)
    >>> TR12i((ta + tb)/(ta*tb - 1))
    -tan(a + b)
    >>> TR12i((-ta - tb)/(ta*tb - 1))
    tan(a + b)
    >>> eq = (ta + tb)/(-ta*tb + 1)**2*(-3*ta - 3*tc)/(2*(ta*tc - 1))
    >>> TR12i(eq.expand())
    -3*tan(a + b)*tan(a + c)/(2*(tan(a) + tan(b) - 1))

    """
    from ..polys import factor

    def f(rv):
        if not (rv.is_Add or rv.is_Mul or rv.is_Pow):
            return rv

        n, d = rv.as_numer_denom()
        if not d.args or not n.args:
            return rv

        dok = {}

        def ok(di):
            m = as_f_sign_1(di)
            if m:
                g, f, s = m
                if s == -1 and f.is_Mul and len(f.args) == 2 and \
                        all(isinstance(fi, tan) for fi in f.args):
                    return g, f

        d_args = list(Mul.make_args(d))
        for i, di in enumerate(d_args):
            m = ok(di)
            if m:
                g, t = m
                s = Add(*[_.args[0] for _ in t.args])
                dok[s] = Integer(1)
                d_args[i] = g
                continue
            if di.is_Add:
                di = factor(di)
                if di.is_Mul:
                    d_args.extend(di.args)
                    d_args[i] = Integer(1)
            elif di.is_Pow and (di.exp.is_integer or di.base.is_positive):
                m = ok(di.base)
                if m:
                    g, t = m
                    s = Add(*[_.args[0] for _ in t.args])
                    dok[s] = di.exp
                    d_args[i] = g**di.exp
                else:
                    di = factor(di)
                    if di.is_Mul:
                        d_args.extend(di.args)
                        d_args[i] = Integer(1)
        if not dok:
            return rv

        def ok2(ni):
            if ni.is_Add and len(ni.args) == 2:
                a, b = ni.args
                if isinstance(a, tan) and isinstance(b, tan):
                    return a, b
        n_args = list(Mul.make_args(factor_terms(n)))
        hit = False
        for i, ni in enumerate(n_args):
            m = ok2(ni)
            if not m:
                m = ok2(-ni)
                if m:
                    n_args[i] = Integer(-1)
                else:
                    if ni.is_Add:
                        ni = factor(ni)
                        if ni.is_Mul:
                            n_args.extend(ni.args)
                            n_args[i] = Integer(1)
                        continue
                    if ni.is_Pow and (
                            ni.exp.is_integer or ni.base.is_positive):
                        m = ok2(ni.base)
                        if m:
                            n_args[i] = Integer(1)
                        else:
                            ni = factor(ni)
                            if ni.is_Mul:
                                n_args.extend(ni.args)
                                n_args[i] = Integer(1)
                            continue
                    else:
                        continue
            else:
                n_args[i] = Integer(1)
            hit = True
            s = Add(*[_.args[0] for _ in m])
            ed = dok[s]
            newed = ed.extract_additively(1)
            if newed is not None:
                if newed:
                    dok[s] = newed
                else:
                    dok.pop(s)
            n_args[i] *= -tan(s)

        if hit:
            rv = Mul(*n_args)/Mul(*d_args)/Mul(*[(Add(*[
                tan(a) for a in i.args]) - 1)**e for i, e in dok.items()])

        return rv

    return bottom_up(rv, f)


def TR13(rv):
    """Change products of ``tan`` or ``cot``.

    Examples
    ========

    >>> TR13(tan(3)*tan(2))
    -tan(2)/tan(5) - tan(3)/tan(5) + 1
    >>> TR13(cot(3)*cot(2))
    cot(2)*cot(5) + 1 + cot(3)*cot(5)

    """

    def f(rv):
        if not rv.is_Mul:
            return rv

        # XXX handle products of powers? or let power-reducing handle it?
        args = {tan: [], cot: [], None: []}
        for a in ordered(Mul.make_args(rv)):
            if a.func in (tan, cot):
                args[a.func].append(a.args[0])
            else:
                args[None].append(a)
        t = args[tan]
        c = args[cot]
        if len(t) < 2 and len(c) < 2:
            return rv
        args = args[None]
        while len(t) > 1:
            t1 = t.pop()
            t2 = t.pop()
            args.append(1 - (tan(t1)/tan(t1 + t2) + tan(t2)/tan(t1 + t2)))
        if t:
            args.append(tan(t.pop()))
        while len(c) > 1:
            t1 = c.pop()
            t2 = c.pop()
            args.append(1 + cot(t1)*cot(t1 + t2) + cot(t2)*cot(t1 + t2))
        if c:
            args.append(cot(c.pop()))
        return Mul(*args)

    return bottom_up(rv, f)


def TRmorrie(rv):
    """Returns cos(x)*cos(2*x)*...*cos(2**(k-1)*x) -> sin(2**k*x)/(2**k*sin(x))

    Examples
    ========

    >>> TRmorrie(cos(x)*cos(2*x))
    sin(4*x)/(4*sin(x))
    >>> TRmorrie(7*Mul(*[cos(x) for x in range(10)]))
    7*sin(12)*sin(16)*cos(5)*cos(7)*cos(9)/(64*sin(1)*sin(3))

    Sometimes autosimplification will cause a power to be
    not recognized. e.g. in the following, cos(4*pi/7) automatically
    simplifies to -cos(3*pi/7) so only 2 of the 3 terms are
    recognized:

    >>> TRmorrie(cos(pi/7)*cos(2*pi/7)*cos(4*pi/7))
    -sin(3*pi/7)*cos(3*pi/7)/(4*sin(pi/7))

    A touch by TR8 resolves the expression to a Rational

    >>> TR8(_)
    -1/8

    In this case, if eq is unsimplified, the answer is obtained
    directly:

    >>> eq = cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9)
    >>> TRmorrie(eq)
    1/16

    But if angles are made canonical with TR3 then the answer
    is not simplified without further work:

    >>> TR3(eq)
    sin(pi/18)*cos(pi/9)*cos(2*pi/9)/2
    >>> TRmorrie(_)
    sin(pi/18)*sin(4*pi/9)/(8*sin(pi/9))
    >>> TR8(_)
    cos(7*pi/18)/(16*sin(pi/9))
    >>> TR3(_)
    1/16

    The original expression would have resolve to 1/16 directly with TR8,
    however:

    >>> TR8(eq)
    1/16

    References
    ==========

    https://en.wikipedia.org/wiki/Morrie%27s_law

    """

    def f(rv):
        if not rv.is_Mul:
            return rv

        args = defaultdict(list)
        coss = {}
        other = []
        for c in rv.args:
            b, e = c.as_base_exp()
            if e.is_Integer and isinstance(b, cos):
                co, a = b.args[0].as_coeff_Mul()
                args[a].append(co)
                coss[b] = e
            else:
                other.append(c)

        new = []
        for a in args:
            c = args[a]
            c.sort()
            no = []
            while c:
                k = 0
                cc = ci = c[0]
                while cc in c:
                    k += 1
                    cc *= 2
                if k > 1:
                    newarg = sin(2**k*ci*a)/2**k/sin(ci*a)
                    # see how many times this can be taken
                    take = None
                    ccs = []
                    for _ in range(k):
                        cc /= 2
                        key = cos(a*cc, evaluate=False)
                        ccs.append(cc)
                        take = min(coss[key], take or coss[key])
                    # update exponent counts
                    for _ in range(k):
                        cc = ccs.pop()
                        key = cos(a*cc, evaluate=False)
                        coss[key] -= take
                        if not coss[key]:
                            c.remove(cc)
                    new.append(newarg**take)
                else:
                    no.append(c.pop(0))
            c[:] = no

        if new:
            rv = Mul(*(new + other + [
                cos(k*a, evaluate=False) for a in args for k in args[a]]))

        return rv

    return bottom_up(rv, f)


def TR14(rv, first=True):
    """Convert factored powers of sin and cos identities into simpler
    expressions.

    Examples
    ========

    >>> TR14((cos(x) - 1)*(cos(x) + 1))
    -sin(x)**2
    >>> TR14((sin(x) - 1)*(sin(x) + 1))
    -cos(x)**2
    >>> p1 = (cos(x) + 1)*(cos(x) - 1)
    >>> p2 = (cos(y) - 1)*2*(cos(y) + 1)
    >>> p3 = (3*(cos(y) - 1))*(3*(cos(y) + 1))
    >>> TR14(p1*p2*p3*(x - 1))
    -18*(x - 1)*sin(x)**2*sin(y)**4

    """

    def f(rv):
        if not rv.is_Mul:
            return rv

        if first:
            # sort them by location in numerator and denominator
            # so the code below can just deal with positive exponents
            n, d = rv.as_numer_denom()
            if d != 1:
                newn = TR14(n, first=False)
                newd = TR14(d, first=False)
                if newn != n or newd != d:
                    rv = newn/newd
                return rv

        other = []
        process = []
        for a in rv.args:
            if a.is_Pow:
                b, e = a.as_base_exp()
                if not (e.is_integer or b.is_positive):
                    other.append(a)
                    continue
                a = b
            else:
                e = Integer(1)
            m = as_f_sign_1(a)
            if not m or m[1].func not in (cos, sin):
                if e == 1:
                    other.append(a)
                else:
                    other.append(a**e)
                continue
            g, f, si = m
            process.append((g, e.is_Number, e, f, si, a))

        # sort them to get like terms next to each other
        process = list(ordered(process))

        # keep track of whether there was any change
        nother = len(other)

        # access keys
        keys = (g, t, e, f, si, a) = list(range(6))

        while process:
            A = process.pop(0)
            if process:
                B = process[0]

                if A[e].is_Number and B[e].is_Number:
                    # both exponents are numbers
                    if A[f] == B[f]:
                        if A[si] != B[si]:
                            B = process.pop(0)
                            take = min(A[e], B[e])

                            # reinsert any remainder
                            # the B will likely sort after A so check it first
                            if B[e] != take:
                                rem = [B[i] for i in keys]
                                rem[e] -= take
                                process.insert(0, rem)
                            elif A[e] != take:
                                rem = [A[i] for i in keys]
                                rem[e] -= take
                                process.insert(0, rem)

                            if isinstance(A[f], cos):
                                t = sin
                            else:
                                t = cos
                            other.append((-A[g]*B[g]*t(A[f].args[0])**2)**take)
                            continue

                elif A[e] == B[e]:
                    # both exponents are equal symbols
                    if A[f] == B[f]:
                        if A[si] != B[si]:
                            B = process.pop(0)
                            take = A[e]
                            if isinstance(A[f], cos):
                                t = sin
                            else:
                                t = cos
                            other.append((-A[g]*B[g]*t(A[f].args[0])**2)**take)
                            continue

            # either we are done or neither condition above applied
            other.append(A[a]**A[e])

        if len(other) != nother:
            rv = Mul(*other)

        return rv

    return bottom_up(rv, f)


def TR15(rv, max=4, pow=False):
    """Convert sin(x)*-2 to 1 + cot(x)**2.

    See _TR56 docstring for advanced use of ``max`` and ``pow``.

    Examples
    ========

    >>> TR15(1 - 1/sin(x)**2)
    -cot(x)**2

    """

    def f(rv):
        if not (isinstance(rv, Pow) and isinstance(rv.base, sin)):
            return rv

        ia = 1/rv
        a = _TR56(ia, sin, cot, lambda x: 1 + x, max=max, pow=pow)
        if a != ia:
            rv = a
        return rv

    return bottom_up(rv, f)


def TR16(rv, max=4, pow=False):
    """Convert cos(x)*-2 to 1 + tan(x)**2.

    See _TR56 docstring for advanced use of ``max`` and ``pow``.

    Examples
    ========

    >>> TR16(1 - 1/cos(x)**2)
    -tan(x)**2

    """

    def f(rv):
        if not (isinstance(rv, Pow) and isinstance(rv.base, cos)):
            return rv

        ia = 1/rv
        a = _TR56(ia, cos, tan, lambda x: 1 + x, max=max, pow=pow)
        if a != ia:
            rv = a
        return rv

    return bottom_up(rv, f)


def TR111(rv):
    """Convert f(x)**-i to g(x)**i where either ``i`` is an integer
    or the base is positive and f, g are: tan, cot; sin, csc; or cos, sec.

    Examples
    ========

    >>> TR111(1 - 1/tan(x)**2)
    -cot(x)**2 + 1

    """

    def f(rv):
        if not (isinstance(rv, Pow) and
                (rv.base.is_positive or rv.exp.is_integer and rv.exp.is_negative)):
            return rv

        if isinstance(rv.base, tan):
            return cot(rv.base.args[0])**-rv.exp
        elif isinstance(rv.base, sin):
            return csc(rv.base.args[0])**-rv.exp
        elif isinstance(rv.base, cos):
            return sec(rv.base.args[0])**-rv.exp
        return rv

    return bottom_up(rv, f)


def TR22(rv, max=4, pow=False):
    """Convert tan(x)**2 to sec(x)**2 - 1 and cot(x)**2 to csc(x)**2 - 1.

    See _TR56 docstring for advanced use of ``max`` and ``pow``.

    Examples
    ========

    >>> TR22(1 + tan(x)**2)
    sec(x)**2
    >>> TR22(1 + cot(x)**2)
    csc(x)**2

    """

    def f(rv):
        if not (isinstance(rv, Pow) and rv.base.func in (cot, tan)):
            return rv

        rv = _TR56(rv, tan, sec, lambda x: x - 1, max=max, pow=pow)
        rv = _TR56(rv, cot, csc, lambda x: x - 1, max=max, pow=pow)
        return rv

    return bottom_up(rv, f)


def TRpower(rv):
    """Convert sin(x)**n and cos(x)**n with positive n to sums.

    Examples
    ========

    >>> TRpower(sin(x)**6)
    -15*cos(2*x)/32 + 3*cos(4*x)/16 - cos(6*x)/32 + 5/16
    >>> TRpower(sin(x)**3*cos(2*x)**4)
    (3*sin(x)/4 - sin(3*x)/4)*(cos(4*x)/2 + cos(8*x)/8 + 3/8)

    References
    ==========

    * https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Power-reduction_formulae

    """

    def f(rv):
        if not (isinstance(rv, Pow) and isinstance(rv.base, (sin, cos))):
            return rv
        b, n = rv.as_base_exp()
        x = b.args[0]
        if n.is_Integer and n.is_positive:
            if n.is_odd and isinstance(b, cos):
                rv = 2**(1-n)*Add(*[binomial(n, k)*cos((n - 2*k)*x)
                                    for k in range((n + 1)/2)])
            elif n.is_odd and isinstance(b, sin):
                rv = (2**(1-n)*(-1)**((n-1)/2) *
                      Add(*[binomial(n, k)*(-1)**k*sin((n - 2*k)*x)
                            for k in range((n + 1)/2)]))
            elif n.is_even and isinstance(b, cos):
                rv = 2**(1-n)*Add(*[binomial(n, k)*cos((n - 2*k)*x)
                                    for k in range(n/2)])
            elif n.is_even and isinstance(b, sin):
                rv = (2**(1-n)*(-1)**(n/2) *
                      Add(*[binomial(n, k)*(-1)**k*cos((n - 2*k)*x)
                            for k in range(n/2)]))
            if n.is_even:
                rv += 2**(-n)*binomial(n, n/2)
        return rv

    return bottom_up(rv, f)


def L(rv):
    """Return count of trigonometric functions in expression.

    Examples
    ========

    >>> L(cos(x)+sin(x))
    2

    """
    return Integer(rv.count(TrigonometricFunction))


# ============== end of basic Fu-like tools =====================

# tuples are chains  --  (f, g) -> lambda x: g(f(x))
# lists are choices  --  [f, g] -> lambda x: min(f(x), g(x), key=objective)

CTR1 = [(TR5, TR0), (TR6, TR0), identity]

CTR2 = (TR11, [(TR5, TR0), (TR6, TR0), TR0])

CTR3 = [(TRmorrie, TR8, TR0), (TRmorrie, TR8, TR10i, TR0), identity]

CTR4 = [(TR4, TR10i), identity]

RL1 = (TR4, TR3, TR4, TR12, TR4, TR13, TR4, TR0)


# XXX it's a little unclear how this one is to be implemented
# see Fu paper of reference, page 7. What is the Union symbol refering to?
# The diagram shows all these as one chain of transformations, but the
# text refers to them being applied independently. Also, a break
# if L starts to increase has not been implemented.
RL2 = [
    (TR4, TR3, TR10, TR4, TR3, TR11),
    (TR5, TR7, TR11, TR4),
    (CTR3, CTR1, TR9, CTR2, TR4, TR9, TR9, CTR4),
    identity,
]


def fu(rv, measure=lambda x: (L(x), x.count_ops())):
    """Attempt to simplify expression by using transformation rules given
    in the algorithm by Fu et al.

    :func:`fu` will try to minimize the objective function ``measure``.
    By default this first minimizes the number of trig terms and then minimizes
    the number of total operations.

    Examples
    ========

    >>> fu(sin(50)**2 + cos(50)**2 + sin(pi/6))
    3/2
    >>> fu(sqrt(6)*cos(x) + sqrt(2)*sin(x))
    2*sqrt(2)*sin(x + pi/3)

    CTR1 example

    >>> eq = sin(x)**4 - cos(y)**2 + sin(y)**2 + 2*cos(x)**2
    >>> fu(eq)
    cos(x)**4 - 2*cos(y)**2 + 2

    CTR2 example

    >>> fu(Rational(1, 2) - cos(2*x)/2)
    sin(x)**2

    CTR3 example

    >>> fu(sin(a)*(cos(b) - sin(b)) + cos(a)*(sin(b) + cos(b)))
    sqrt(2)*sin(a + b + pi/4)

    CTR4 example

    >>> fu(sqrt(3)*cos(x)/2 + sin(x)/2)
    sin(x + pi/3)

    Example 1

    >>> fu(1-sin(2*x)**2/4-sin(y)**2-cos(x)**4)
    -cos(x)**2 + cos(y)**2

    Example 2

    >>> fu(cos(4*pi/9))
    sin(pi/18)
    >>> fu(cos(pi/9)*cos(2*pi/9)*cos(3*pi/9)*cos(4*pi/9))
    1/16

    Example 3

    >>> fu(tan(7*pi/18)+tan(5*pi/18)-sqrt(3)*tan(5*pi/18)*tan(7*pi/18))
    -sqrt(3)

    Objective function example

    >>> fu(sin(x)/cos(x))  # default objective function
    tan(x)
    >>> fu(sin(x)/cos(x), measure=lambda x: -x.count_ops())  # maximize op count
    sin(x)/cos(x)

    References
    ==========
    http://rfdz.ph-noe.ac.at/fileadmin/Mathematik_Uploads/ACDCA/
    DESTIME2006/DES_contribs/Fu/simplification.pdf

    """
    fRL1 = greedy(RL1, measure)
    fRL2 = greedy(RL2, measure)

    was = rv
    rv = sympify(rv)
    if not isinstance(rv, Expr):
        return rv.func(*[fu(a, measure=measure) for a in rv.args])
    rv = TR1(rv)
    if rv.has(tan, cot):
        rv1 = fRL1(rv)
        if measure(rv1) < measure(rv):
            rv = rv1
        if rv.has(tan, cot):
            rv = TR2(rv)
    if rv.has(sin, cos):
        rv1 = fRL2(rv)
        rv2 = TR8(TRmorrie(rv1))
        rv = min([was, rv, rv1, rv2], key=measure)
    return min(TR2i(rv), rv, key=measure)


def process_common_addends(rv, do, key2=None, key1=True):
    """Apply ``do`` to addends of ``rv`` that (if key1=True) share at least
    a common absolute value of their coefficient and the value of ``key2`` when
    applied to the argument. If ``key1`` is False ``key2`` must be supplied and
    will be the only key applied.

    """
    # collect by absolute value of coefficient and key2
    absc = defaultdict(list)
    if key1:
        for a in rv.args:
            c, a = a.as_coeff_Mul()
            if c < 0:
                c = -c
                a = -a  # put the sign on `a`
            absc[(c, key2(a) if key2 else 1)].append(a)
    elif key2:
        for a in rv.args:
            absc[(Integer(1), key2(a))].append(a)
    else:
        raise ValueError('must have at least one key')

    args = []
    hit = False
    for k in absc:
        v = absc[k]
        c, _ = k
        if len(v) > 1:
            e = Add(*v, evaluate=False)
            new = do(e)
            if new != e:
                e = new
                hit = True
            args.append(c*e)
        else:
            args.append(c*v[0])
    if hit:
        rv = Add(*args)

    return rv


fufuncs = """
    TR0 TR1 TR2 TR3 TR4 TR5 TR6 TR7 TR8 TR9 TR10 TR10i TR11
    TR12 TR13 L TR2i TRmorrie TR12i
    TR14 TR15 TR16 TR111 TR22""".split()
FU = dict(zip(fufuncs, list(map(locals().get, fufuncs))))


def _roots():
    global _ROOT2, _ROOT3, _invROOT3
    _ROOT2, _ROOT3 = sqrt(2), sqrt(3)
    _invROOT3 = 1/_ROOT3


_ROOT2 = None


def trig_split(a, b, two=False):
    """Return the gcd, s1, s2, a1, a2, bool where

    If two is False (default) then::
        a + b = gcd*(s1*f(a1) + s2*f(a2)) where f = cos if bool else sin
    else:
        if bool, a + b was +/- cos(a1)*cos(a2) +/- sin(a1)*sin(a2) and equals
            n1*gcd*cos(a - b) if n1 == n2 else
            n1*gcd*cos(a + b)
        else a + b was +/- cos(a1)*sin(a2) +/- sin(a1)*cos(a2) and equals
            n1*gcd*sin(a + b) if n1 = n2 else
            n1*gcd*sin(b - a)

    Examples
    ========

    >>> trig_split(cos(x), cos(y))
    (1, 1, 1, x, y, True)
    >>> trig_split(2*cos(x), -2*cos(y))
    (2, 1, -1, x, y, True)
    >>> trig_split(cos(x)*sin(y), cos(y)*sin(y))
    (sin(y), 1, 1, x, y, True)

    >>> trig_split(cos(x), -sqrt(3)*sin(x), two=True)
    (2, 1, -1, x, pi/6, False)
    >>> trig_split(cos(x), sin(x), two=True)
    (sqrt(2), 1, 1, x, pi/4, False)
    >>> trig_split(cos(x), -sin(x), two=True)
    (sqrt(2), 1, -1, x, pi/4, False)
    >>> trig_split(sqrt(2)*cos(x), -sqrt(6)*sin(x), two=True)
    (2*sqrt(2), 1, -1, x, pi/6, False)
    >>> trig_split(-sqrt(6)*cos(x), -sqrt(2)*sin(x), two=True)
    (-2*sqrt(2), 1, 1, x, pi/3, False)
    >>> trig_split(cos(x)/sqrt(6), sin(x)/sqrt(2), two=True)
    (sqrt(6)/3, 1, 1, x, pi/6, False)
    >>> trig_split(-sqrt(6)*cos(x)*sin(y), -sqrt(2)*sin(x)*sin(y), two=True)
    (-2*sqrt(2)*sin(y), 1, 1, x, pi/3, False)

    >>> trig_split(cos(x), sin(x))
    >>> trig_split(cos(x), sin(z))
    >>> trig_split(2*cos(x), -sin(x))
    >>> trig_split(cos(x), -sqrt(3)*sin(x))
    >>> trig_split(cos(x)*cos(y), sin(x)*sin(z))
    >>> trig_split(cos(x)*cos(y), sin(x)*sin(y))
    >>> trig_split(-sqrt(6)*cos(x), sqrt(2)*sin(x)*sin(y), two=True)

    """
    global _ROOT2, _ROOT3, _invROOT3
    if _ROOT2 is None:
        _roots()

    a, b = map(Factors, (a, b))
    ua, ub = a.normal(b)
    gcd = a.gcd(b).as_expr()
    n1 = n2 = 1
    if -1 in ua.factors:
        ua = ua.quo(Integer(-1))
        n1 = -n1
    elif -1 in ub.factors:
        ub = ub.quo(Integer(-1))
        n2 = -n2
    a, b = (i.as_expr() for i in (ua, ub))

    def pow_cos_sin(a, two):
        """Return ``a`` as a tuple (r, c, s) such that
        ``a = (r or 1)*(c or 1)*(s or 1)``.

        Three arguments are returned (radical, c-factor, s-factor) as
        long as the conditions set by ``two`` are met; otherwise None is
        returned. If ``two`` is True there will be one or two non-None
        values in the tuple: c and s or c and r or s and r or s or c with c
        being a cosine function (if possible) else a sine, and s being a sine
        function (if possible) else oosine. If ``two`` is False then there
        will only be a c or s term in the tuple.

        ``two`` also require that either two cos and/or sin be present (with
        the condition that if the functions are the same the arguments are
        different or vice versa) or that a single cosine or a single sine
        be present with an optional radical.

        If the above conditions dictated by ``two`` are not met then None
        is returned.

        """
        c = s = None
        co = Integer(1)
        if a.is_Mul:
            co, a = a.as_coeff_Mul()
            if len(a.args) > 2 or not two:
                return
            if a.is_Mul:
                args = list(a.args)
            else:
                args = [a]
            a = args.pop(0)
            if isinstance(a, cos):
                c = a
            elif isinstance(a, sin):
                s = a
            elif a.is_Pow and a.exp == Rational(1, 2):  # autoeval doesn't allow -1/2
                co *= a
            else:
                return
            if args:
                b = args[0]
                if isinstance(b, cos):
                    if c:
                        s = b
                    else:
                        c = b
                elif isinstance(b, sin):
                    if s:
                        c = b
                    else:
                        s = b
                elif b.is_Pow and b.exp == Rational(1, 2):
                    co *= b
                else:
                    return
            return co if co != 1 else None, c, s
        elif isinstance(a, cos):
            c = a
        elif isinstance(a, sin):
            s = a
        if c is None and s is None:
            return
        co = co if co != 1 else None
        return co, c, s

    # get the parts
    m = pow_cos_sin(a, two)
    if m is None:
        return
    coa, ca, sa = m
    m = pow_cos_sin(b, two)
    if m is None:
        return
    cob, cb, sb = m

    # check them
    if (not ca) and cb or ca and isinstance(ca, sin):
        coa, ca, sa, cob, cb, sb = cob, cb, sb, coa, ca, sa
        n1, n2 = n2, n1
    if not two:  # need cos(x) and cos(y) or sin(x) and sin(y)
        c = ca or sa
        s = cb or sb
        if not isinstance(c, s.func):
            return
        return gcd, n1, n2, c.args[0], s.args[0], isinstance(c, cos)
    else:
        if not coa and not cob:
            if (ca and cb and sa and sb):
                if not isinstance(ca, sa.func) is isinstance(cb, sb.func):
                    return
                args = {j.args for j in (ca, sa)}
                if not all(i.args in args for i in (cb, sb)):
                    return
                return gcd, n1, n2, ca.args[0], sa.args[0], isinstance(ca, sa.func)
        if ca and sa or cb and sb or \
           two and (ca is None and sa is None or cb is None and sb is None):
            return
        c = ca or sa
        s = cb or sb
        if c.args != s.args:
            return
        if not coa:
            coa = Integer(1)
        if not cob:
            cob = Integer(1)
        if coa is cob:
            gcd *= _ROOT2
            return gcd, n1, n2, c.args[0], pi/4, False
        elif coa/cob == _ROOT3:
            gcd *= 2*cob
            return gcd, n1, n2, c.args[0], pi/3, False
        elif coa/cob == _invROOT3:
            gcd *= 2*coa
            return gcd, n1, n2, c.args[0], pi/6, False


def as_f_sign_1(e):
    """If ``e`` is a sum that can be written as ``g*(a + s)`` where
    ``s`` is ``+/-1``, return ``g``, ``a``, and ``s`` where ``a`` does
    not have a leading negative coefficient.

    Examples
    ========

    >>> as_f_sign_1(x + 1)
    (1, x, 1)
    >>> as_f_sign_1(x - 1)
    (1, x, -1)
    >>> as_f_sign_1(-x + 1)
    (-1, x, -1)
    >>> as_f_sign_1(-x - 1)
    (-1, x, 1)
    >>> as_f_sign_1(2*x + 2)
    (2, x, 1)

    """
    if not e.is_Add or len(e.args) != 2:
        return
    # exact match
    a, b = e.args
    if a in (-1, 1):
        g = Integer(1)
        if b.is_Mul and b.args[0].is_Number and b.args[0] < 0:
            a, b = -a, -b
            g = -g
        return g, b, a
    # gcd match
    a, b = map(Factors, e.args)
    ua, ub = a.normal(b)
    gcd = a.gcd(b).as_expr()
    if -1 in ua.factors:
        ua = ua.quo(Integer(-1))
        n1 = -1
        n2 = 1
    elif -1 in ub.factors:
        ub = ub.quo(Integer(-1))
        n1 = 1
        n2 = -1
    else:
        n1 = n2 = 1
    a, b = (i.as_expr() for i in (ua, ub))
    if a == 1:
        a, b = b, a
        n1, n2 = n2, n1
    if n1 == -1:
        gcd = -gcd
        n2 = -n2

    if b == 1:
        return gcd, a, n2


def _osborne(e, d):
    """Replace all hyperbolic functions with trig functions using
    the Osborne rule.

    Notes
    =====

    ``d`` is a dummy variable to prevent automatic evaluation
    of trigonometric/hyperbolic functions.


    References
    ==========

    https://en.wikipedia.org/wiki/Hyperbolic_function

    """

    def f(rv):
        if not isinstance(rv, HyperbolicFunction):
            return rv
        a = rv.args[0]
        a = a*d if not a.is_Add else Add._from_args([i*d for i in a.args])
        if isinstance(rv, sinh):
            return I*sin(a)
        elif isinstance(rv, cosh):
            return cos(a)
        elif isinstance(rv, tanh):
            return I*tan(a)
        elif isinstance(rv, coth):
            return cot(a)/I
        else:
            raise NotImplementedError(f'unhandled {rv.func}')

    return bottom_up(e, f)


def _osbornei(e, d):
    """Replace all trig functions with hyperbolic functions using
    the Osborne rule.

    Notes
    =====

    ``d`` is a dummy variable to prevent automatic evaluation
    of trigonometric/hyperbolic functions.

    References
    ==========

    https://en.wikipedia.org/wiki/Hyperbolic_function

    """

    def f(rv):
        if not isinstance(rv, TrigonometricFunction):
            return rv
        a = rv.args[0].xreplace({d: Integer(1)})
        if isinstance(rv, sin):
            return sinh(a)/I
        elif isinstance(rv, cos):
            return cosh(a)
        elif isinstance(rv, tan):
            return tanh(a)/I
        elif isinstance(rv, cot):
            return coth(a)*I
        elif isinstance(rv, sec):
            return 1/cosh(a)
        elif isinstance(rv, csc):
            return I/sinh(a)
        else:
            raise NotImplementedError(f'unhandled {rv.func}')

    return bottom_up(e, f)


def hyper_as_trig(rv):
    """Return an expression containing hyperbolic functions in terms
    of trigonometric functions. Any trigonometric functions initially
    present are replaced with Dummy symbols and the function to undo
    the masking and the conversion back to hyperbolics is also returned. It
    should always be true that::

        t, f = hyper_as_trig(expr)
        expr == f(t)

    Examples
    ========

    >>> eq = sinh(x)**2 + cosh(x)**2
    >>> t, f = hyper_as_trig(eq)
    >>> f(fu(t))
    cosh(2*x)

    References
    ==========

    https://en.wikipedia.org/wiki/Hyperbolic_function

    """
    from .radsimp import collect
    from .simplify import signsimp

    # mask off trig functions
    trigs = rv.atoms(TrigonometricFunction)
    reps = [(t, Dummy()) for t in trigs]
    masked = rv.xreplace(dict(reps))

    # get inversion substitutions in place
    reps = [(v, k) for k, v in reps]

    d = Dummy()

    return _osborne(masked, d), lambda x: collect(signsimp(
        _osbornei(x, d).xreplace(dict(reps))), I)


def sincos_to_sum(expr):
    """Convert products and powers of sin and cos to sums.

    Applied power reduction TRpower first, then expands products, and
    converts products to sums with TR8.

    Examples
    ========

    >>> sincos_to_sum(16*sin(x)**3*cos(2*x)**2)
    7*sin(x) - 5*sin(3*x) + 3*sin(5*x) - sin(7*x)

    """
    return TR8(expand_mul(TRpower(expr)))
