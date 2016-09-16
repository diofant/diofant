from diofant.core import oo, nan, diff, sympify
from diofant.sets import Interval
from diofant.core.compatibility import is_sequence
from diofant.series import limit
from diofant.functions import Min
from diofant.solvers import solve
from diofant.calculus import singularities


__all__ = ('minimize', 'maximize')


def minimize(f, *v):
    """Minimizes `f` with respect to given variables `v`.

    Examples
    ========

    >>> from diofant.calculus import minimize
    >>> from diofant.abc import x
    >>> minimize(x**2, x)
    (0, {x: 0})

    >>> minimize([x**2, x >= 1], x)
    (1, {x: 1})
    >>> minimize([-x**2, x >= -2, x <= 1], x)
    (-4, {x: -2})

    See Also
    ========

    maximize
    """
    f = set(map(sympify, f if is_sequence(f) else [f]))

    constr = {c for c in f if c.is_Relational}

    assert len(f - constr) == 1

    f = (f - constr).pop()

    if not v:
        v = f.free_symbols
        if not v:
            return f, dict()
        v = tuple(v)

    assert all(x.is_Symbol for x in v)

    if constr:
        dom = solve(constr, *v).as_set()
    else:
        dom = Interval(-oo, oo, True, True)**len(v)

    if len(v) == 1:
        return minimize_univariate(f, v[0], dom)
    else:  # pragma: no cover
        return NotImplementedError


def maximize(f, *v):
    """
    Maximizes `f` with respect to given variables `v`.

    See Also
    ========

    minimize
    """
    f = set(map(sympify, f if is_sequence(f) else [f]))

    fv, d = minimize([e if e.is_Relational else -e for e in f], *v)
    return -fv, d


def minimize_univariate(f, x, dom):
    extr = {}

    if dom.is_Union:
        for d in dom.args:
            fp, r = minimize_univariate(f, x, d)
            extr[r[x]] = fp
    elif dom.is_Interval:
        if not dom.left_open:
            extr[dom.start] = limit(f, x, dom.start)
        if not dom.right_open:
            extr[dom.end] = limit(f, x, dom.end, dir="-")
        for s in singularities(f, x):
            if s in dom:
                m = Min(limit(f, x, s), limit(f, x, s, dir="-"))
                if m is -oo:
                    return -oo, dict({x: s})
                else:
                    extr[s] = m

        for p in solve(diff(f, x), x):
            if p in dom:
                extr[p] = f.subs(x, p)
    elif dom.is_FiniteSet:
        for p in dom.args:
            extr[p] = f.subs(x, p)
    else:  # pragma: no cover
        raise NotImplementedError

    if extr:
        min, point = oo, nan
        for p, fp in sorted(extr.items()):
            if fp < min:
                point, min = p, fp
        return min, dict({x: point})
