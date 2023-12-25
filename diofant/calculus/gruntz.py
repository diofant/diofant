r"""
The Gruntz Algorithm
====================

This section explains the basics of the algorithm :cite:`Gruntz1996limits` used for computing
limits.  Most of the time the :py:func:`~diofant.calculus.limits.limit` function
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

       .. math:: f(\omega) = c_0 \omega^{e_0} + c_1 \omega^{e_1} + \ldots

   where `e_0 < e_1 < e_2 \ldots`

4. If the leading exponent `e_0 > 0` then the limit is `0`.  If `e_0 < 0`,
   then the answer is `\pm\infty` (depends on sign of `c_0`).  Finally,
   if `e_0 = 0`, the limit is the limit of the leading coefficient `c_0`.

Notes
-----
This exposition glossed over several details.  For example, limits could be
computed recursively (steps 1 and 4).  Please address to the Gruntz thesis :cite:`Gruntz1996limits`
for proof of the termination (pp. 52-60).

"""

import functools

from ..core import Add, Dummy, Integer, Mul, cacheit, evaluate, oo
from ..core.function import UndefinedFunction
from ..utilities import ordered


def compare(a, b, x):
    r"""
    Determine order relation between two functons.

    Returns
    =======

    {1, 0, -1}
        Respectively, if `a(x) \succ b(x)`, `a(x) \asymp b(x)`
        or `b(x) \succ a(x)`.

    Examples
    ========

    >>> compare(exp(x), x**5, x)
    1

    """
    from ..functions import log

    # The log(exp(...)) must always be simplified here for termination.
    la = a.exp if a.is_Exp else log(a)
    lb = b.exp if b.is_Exp else log(b)

    c = limitinf(la/lb, x)
    if c.is_zero:
        return -1
    if c.is_infinite:
        return 1
    return 0


def mrv(e, x):
    """
    Calculate the MRV set of the expression.

    Examples
    ========

    >>> mrv(log(x - log(x))/log(x), x)
    {x}

    """
    from ..functions import log

    if not e.has(x):
        return set()
    if e == x:
        return {x}
    if e.is_Mul or e.is_Add:
        a, b = e.as_two_terms()
        return _mrv_max(mrv(a, x), mrv(b, x), x)
    if e.is_Exp:
        if e.exp == x:
            return {e}
        if any(a.is_infinite for a in Mul.make_args(limitinf(e.exp, x))):
            return _mrv_max({e}, mrv(e.exp, x), x)
        return mrv(e.exp, x)
    if e.is_Pow:
        return mrv(e.base, x)
    if isinstance(e, log):
        return mrv(e.args[0], x)
    if e.is_Function and not isinstance(e.func, UndefinedFunction):
        return functools.reduce(lambda a, b: _mrv_max(a, b, x),
                                [mrv(a, x) for a in e.args])
    raise NotImplementedError(f"Can't calculate the MRV of {e}.")


def _mrv_max(f, g, x):
    """Compute the maximum of two MRV sets."""
    for a, b in zip(f, g):
        if (c := compare(a, b, x)) in (1, -1):
            return f if c > 0 else g
        break
    return f | g


@cacheit
def signinf(e, x):
    r"""
    Determine sign of the expression at the infinity.

    Returns
    =======

    {1, 0, -1}
        One or minus one, if `e > 0` or `e < 0` for `x` sufficiently
        large and zero if `e` is *constantly* zero for `x\to\infty`.

    """
    from ..functions import sign

    if not e.has(x):
        return sign(e).simplify()
    if e == x:
        return Integer(1)
    if e.is_Mul:
        a, b = e.as_two_terms()
        return signinf(a, x)*signinf(b, x)
    if e.is_Pow and signinf(e.base, x) == 1:
        return Integer(1)

    c0, _ = leadterm(e, x)
    return signinf(c0, x)


@cacheit
def limitinf(e, x):
    """
    Compute the limit of the expression at the infinity.

    Examples
    ========

    >>> limitinf(exp(x)*(exp(1/x - exp(-x)) - exp(1/x)), x)
    -1

    """
    # Rewrite e in terms of tractable functions only:
    e = e.rewrite('tractable', wrt=x)

    if not e.has(x):
        # This is a bit of a heuristic for nice results.  We always rewrite
        # tractable functions in terms of familiar intractable ones.
        # TODO: It might be nicer to rewrite the exactly to what they were
        # initially, but that would take some work to implement.
        return e.rewrite('intractable')

    c0, e0 = leadterm(e, x)
    sig = signinf(e0, x)
    if sig == 1:
        return Integer(0)
    if sig == -1:
        return signinf(c0, x)*oo
    if sig == 0:
        return limitinf(c0, x)
    raise NotImplementedError(f'Result depends on the sign of {sig}.')


@cacheit
def leadterm(e, x):
    """
    Compute the leading term of the series.

    Returns
    =======

    tuple
        The leading term `c_0 w^{e_0}` of the series of `e` in terms
        of the most rapidly varying subexpression `w` in form of
        the pair ``(c0, e0)`` of Expr.

    Examples
    ========

    >>> leadterm(1/exp(-x + exp(-x)) - exp(x), x)
    (-1, 0)

    """
    from ..functions import exp, log

    if not e.has(x):
        return e, Integer(0)

    e = e.replace(lambda f: f.is_Pow and f.exp.has(x),
                  lambda f: exp(log(f.base)*f.exp))
    e = e.replace(lambda f: f.is_Mul and sum(a.is_Pow for a in f.args) > 1,
                  lambda f: Mul(exp(Add(*[a.exp for a in f.args if a.is_Exp])),
                                *[a for a in f.args if not a.is_Exp]))

    # The positive dummy, w, is used here so log(w*2) etc. will expand.
    # TODO: For limits of complex functions, the algorithm would have to
    # be improved, or just find limits of Re and Im components separately.
    w = Dummy('w', real=True, positive=True)
    e, logw = rewrite(e, x, w)

    lt = e.compute_leading_term(w, logx=logw)
    c0, e0 = lt.as_coeff_exponent(w)
    if c0.has(w):
        raise NotImplementedError(f'Cannot compute leadterm({e}, {x}). '
                                  'The coefficient should have been free of '
                                  f'{w}, but got {c0}.')
    return c0, e0


def rewrite(e, x, w):
    r"""
    Rewrites the expression in terms of the MRV subexpression.

    Parameters
    ==========

    e : Expr
        an expression
    x : Symbol
        variable of the `e`
    w : Symbol
        The symbol which is going to be used for substitution in place
        of the MRV in `x` subexpression.

    Returns
    =======

    tuple
        A pair: rewritten (in `w`) expression and `\log(w)`.

    Examples
    ========

    >>> rewrite(exp(x)*log(x), x, y)
    (log(x)/y, -x)

    """
    from ..functions import exp

    Omega = mrv(e, x)
    if not Omega:
        return e, None  # e really does not depend on x

    if x in Omega:
        # Moving up in the asymptotical scale:
        with evaluate(False):
            e = e.xreplace({x: exp(x)})
            Omega = {s.xreplace({x: exp(x)}) for s in Omega}

    Omega = list(ordered(Omega, keys=lambda a: -len(mrv(a, x))))

    for g in Omega:
        sig = signinf(g.exp, x)
        if sig not in (1, -1):
            raise NotImplementedError(f'Result depends on the sign of {sig}.')

    if sig == 1:
        w = 1/w  # if g goes to oo, substitute 1/w

    # Rewrite and substitute subexpressions in the Omega.
    for a in Omega:
        c = limitinf(a.exp/g.exp, x)
        b = exp(a.exp - c*g.exp)*w**c  # exponential must never be expanded here
        with evaluate(False):
            e = e.xreplace({a: b})

    return e, -sig*g.exp
