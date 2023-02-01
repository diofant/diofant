from ..core import Add, Mul, PoleError, Pow, sympify
from ..sets import Reals


def singularities(f, x):
    """Find singularities of real-valued function `f` with respect to `x`.

    Examples
    ========

    >>> singularities(1/(1 + x), x)
    {-1}

    >>> singularities(exp(1/x) + log(x + 1), x)
    {-1, 0}

    >>> singularities(exp(1/log(x + 1)), x)
    {0}

    Notes
    =====

    Removable singularities are not supported now.

    References
    ==========

    * https://en.wikipedia.org/wiki/Mathematical_singularity

    """
    from ..functions import log, sign
    from ..solvers import solve
    from .limits import Limit

    f, x = sympify(f), sympify(x)
    guess, res = set(), set()

    assert x.is_Symbol

    if f.is_number:
        return set()
    if f.is_polynomial(x):
        return set()
    if f.func in (Add, Mul):
        guess = guess.union(*[singularities(a, x) for a in f.args])
    elif isinstance(f, Pow):
        if f.exp.is_number and f.exp.is_negative:
            guess = {s[x] for s in solve(f.base, x) if s[x].is_real}
        else:
            guess |= singularities(log(f.base)*f.exp, x)
    elif f.func in (log, sign) and len(f.args) == 1:
        guess |= singularities(f.args[0], x)
        guess |= {s[x] for s in solve(f.args[0], x) if s[x].is_real}
    else:
        raise NotImplementedError

    for s in guess:
        l = Limit(f, x, s, dir=Reals)
        try:
            r = l.doit()
            if r.is_infinite:
                raise PoleError
            raise NotImplementedError
        except PoleError:
            res.add(s)

    return res
