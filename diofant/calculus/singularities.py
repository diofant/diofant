from ..core import Add, Mul, PoleError, Pow, log, sympify
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
    from ..functions import Abs, cos, sign, sin
    from ..solvers import solve
    from .limits import limit

    f, x = sympify(f), sympify(x)

    assert x.is_Symbol

    if f.is_number:
        return set()
    if f.is_polynomial(x):
        return set()
    if f.func in (Add, Mul):
        res = set().union(*[singularities(a, x) for a in f.args])
    elif isinstance(f, Pow):
        if f.exp.is_number and f.exp.is_negative:
            res = {s[x] for s in solve(f.base, x) if s[x].is_real}
        else:
            res = singularities(log(f.base)*f.exp, x)
    elif f.func in (log, sign) and len(f.args) == 1:
        res = singularities(f.args[0], x)
        res |= {s[x] for s in solve(f.args[0], x) if s[x].is_real}
    elif f.func in (Abs, sin, cos):
        res = singularities(f.args[0], x)
    else:
        raise NotImplementedError

    for s in res.copy():
        try:
            l = limit(f, x, s, dir=Reals)
            if l.is_finite:
                res.remove(s)
        except PoleError:
            pass

    return res
