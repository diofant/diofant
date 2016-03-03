from sympy.core import oo, diff, Dummy, sympify
from sympy.sets import Interval
from sympy.core.compatibility import is_sequence
from sympy.series import limit
from sympy.functions import Min
from sympy.solvers import solve
from sympy.calculus import singularities


__all__ = ['minimize', 'maximize']


def minimize(f, *v):
    """Minimizes `f` with respect to given variables `v`.

    Examples
    ========

    >>> from sympy.calculus import minimize
    >>> from sympy.abc import x
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

    constr = set([c for c in f if c.is_Relational])

    assert len(f - constr) == 1

    f = (f - constr).pop()

    if not v:
        v = f.free_symbols
        if not v:
            raise ValueError
        v = tuple(v)

    assert all(x.is_Symbol for x in v)

    if constr:
        dom = solve(constr, *v).as_set()
    else:
        dom = Interval(-oo, oo)**len(v)

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
    return (-fv, d)


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
                    return (-oo, dict({x: s}))
                else:
                    extr[s] = m

        for p in solve(diff(f, x), x):
            if p in dom:
                extr[p] = f.subs(x, p)
    elif dom.is_FiniteSet:
        for p in dom.args:
            extr[p] = f.subs(x, p)

    if extr:
        m = Min(*extr.values())
        for p, fp in extr.items():
            if fp == m:
                return (m, dict({x: p}))


def simplex(f, cs):
    """
    Simplex algorithm for linear programming.

    Examples
    ========

    >>> from sympy.calculus.optimization import simplex
    >>> from sympy.abc import x, y, z

    >>> simplex(-2*x - 3*y - 2*z, [2*x + y + z <= 4,
    ...                            x + 2*y + z <= 7,
    ...                            z <= 5])
    (11, [0, 3, 1])
    >>> simplex(-2*x - 3*y - 4*z, [3*x + 2*y + z <= 10,
    ...                            2*x + 5*y + 3*z <= 15])
    (20, [0, 0, 5])

    References
    ==========

    .. [1] https://en.wikipedia.org/wiki/Simplex_method
    """
    from sympy import Matrix, default_sort_key, Le

    syms = sorted(f.free_symbols, key=default_sort_key)
    n = len(syms)
    obj = [1] + [f.coeff(s) for s in syms]
    rows, cons = [], []
    for c in cs:
        assert c.func is Le and not c.rhs.has(*syms)
        rows.append([0] + [c.lhs.coeff(s) for s in syms])
        cons.append(c.rhs)

    def pivot_column(obj):
        low, idx = 0, 0
        for i in range(1, len(obj) - 1):
            if obj[i] < low:
                low, idx = obj[i], i
        return -1 if idx == 0 else idx

    def pivot_row(rows, col):
        rhs = [rows[i][-1] for i in range(len(rows))]
        lhs = [rows[i][col] for i in range(len(rows))]
        ratio, idx = oo, 0
        for i in range(len(rhs)):
            if lhs[i] != 0:
                r = rhs[i]/lhs[i]
                if r < ratio:
                    ratio, idx = r, i
        return idx

    # build full tableau
    for i in range(len(rows)):
        obj += [0]
        ident = [0 if r != i else 1 for r in range(len(rows))]
        rows[i] += ident + [cons[i]]
        rows[i] = Matrix(rows[i])
    obj = Matrix(obj + [0])

    # solve
    while min(obj[1:-1]) < 0:
        col = pivot_column(obj)
        row = pivot_row(rows, col)

        e = rows[row][col]
        rows[row] /= e
        for r in range(len(rows)):
            if r == row:
                continue
            rows[r] = rows[r] - rows[r][col]*rows[row]
        obj = obj - obj[col]*rows[row]

    ans = list(range(n))
    for i in range(1, n + 1):
        ans[i - 1] = 0
        for j in range(len(rows)):
            if rows[j][i] != 0:
                if ans[i - 1] != 0:
                    ans[i - 1] = 0
                    break
                ans[i - 1] = rows[j][-1]
    return (obj[-1], ans)
