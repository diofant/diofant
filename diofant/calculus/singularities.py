from diofant.core import sympify, Add, Mul, Pow, PoleError
from diofant.series.limits import Limit
from diofant.functions import log, sign
from diofant.solvers import solve


def singularities(f, x):
    """Find singularities of real-valued function `f` with respect to `x`.

    Examples
    ========

    >>> from diofant import Symbol, exp, log
    >>> from diofant.abc import x

    >>> singularities(1/(1 + x), x) == {-1}
    True

    >>> singularities(exp(1/x) + log(x + 1), x) == {-1, 0}
    True

    >>> singularities(exp(1/log(x + 1)), x) == {0}
    True

    Notes
    =====

    Removable singularities are not supported now.

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Mathematical_singularity
    """
    f, x = sympify(f), sympify(x)
    guess, res = set(), set()

    assert x.is_Symbol

    if f.is_number:
        return set()
    elif f.is_polynomial(x):
        return set()
    elif f.func in (Add, Mul):
        guess = guess.union(*[singularities(a, x) for a in f.args])
    elif f.func is Pow:
        if f.exp.is_number and f.exp.is_negative:
            guess = {v for v in solve(f.base, x) if v.is_real}
        else:
            guess |= singularities(log(f.base)*f.exp, x)
    elif f.func in (log, sign) and len(f.args) == 1:
        guess |= singularities(f.args[0], x)
        guess |= {v for v in solve(f.args[0], x) if v.is_real}
    else:  # pragma: no cover
        raise NotImplementedError

    for s in guess:
        l = Limit(f, x, s, dir="real")
        try:
            r = l.doit()
            if r == l or f.subs(x, s) != r:  # pragma: no cover
                raise NotImplementedError
        except PoleError:
            res.add(s)

    return res
