from ..calculus import singularities
from ..core import Eq, Ge, Gt, Lt, Ne, S, diff, nan, oo, sympify
from ..core.compatibility import is_sequence, ordered
from ..functions import Min
from ..matrices import Matrix, eye, zeros
from ..series import limit
from ..sets import Interval
from ..solvers import reduce_inequalities, solve


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

    constraints = {c for c in f if c.is_Relational}

    assert len(f - constraints) == 1

    obj = (f - constraints).pop()

    if not v:
        v = obj.free_symbols
        if not v:
            return obj, {}
    v = list(ordered(v))
    dim = len(v)

    assert all(x.is_Symbol for x in v)

    # Canonicalize constraints, Ne -> pair Lt
    constraints |= {Lt(*c.args) for c in constraints if isinstance(c, Ne)}
    constraints |= {Lt(c.lts, c.gts) for c in constraints if isinstance(c, Ne)}
    constraints -= {c for c in constraints if isinstance(c, Ne)}

    # Gt/Ge -> Lt, Le
    constraints = {c.reversed if c.func in (Gt, Ge) else c
                   for c in constraints}

    # Now we have only Lt/Le/Eq
    constraints = list(ordered(c.func(c.lhs - c.rhs, 0)
                               for c in constraints))

    polys = [obj.as_poly(*v)] + [c.lhs.as_poly(*v) for c in constraints]
    is_polynomial = all(p is not None for p in polys)
    is_linear = is_polynomial and all(p.is_linear for p in polys)

    # Eliminate equalities, in the linear case for now
    elims = solve([c for c in constraints if isinstance(c, Eq)], *v)
    if elims and is_linear:
        elims = elims[0]
        res, sol = minimize([obj.subs(elims)] +
                            [c.subs(elims)
                             for c in constraints if not isinstance(c, Eq)],
                            *(set(v) - set(elims)))
        return res, {x: x.subs(elims).subs(sol) for x in v}

    if dim == 1:
        if constraints:
            dom = reduce_inequalities(constraints, *v).as_set()
        else:
            dom = Interval(-oo, oo, True, True)**len(v)
        return minimize_univariate(obj, v[0], dom)

    if is_linear:
        # Quick exit for strict forms
        if any(isinstance(c, Lt) for c in constraints):
            return

        # Transform to the standard form: maximize cᵀx with m⋅x≤b, x≥0.
        # We replace original vector of unrestricted variables v with
        # x of doubled size, so e.g. for the first component of v we
        # will have v₁ = x₁⁺ - x₁⁻, where x₁⁺≥0 and x₁⁻≥0.
        c = [-polys[0].coeff_monomial(x) for x in v]
        c.extend([-_ for _ in c])
        m = []
        for p in polys[1:]:
            r = [p.coeff_monomial(x) for x in v]
            m.extend(r)
            m.extend(-_ for _ in r)
        b = [-p.coeff_monomial(1) for p in polys[1:]]
        m = Matrix(m).reshape(len(b), len(c))

        res, sol = simplex(c, m, b)
        res -= polys[0].coeff_monomial(1)
        sol = map(lambda x, y: x - y, sol[:dim], sol[dim:])

        return -res, dict(zip(v, sol))

    raise NotImplementedError  # pragma: no cover


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
            p = p[x]
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


def simplex(c, m, b):
    """
    Simplex algorithm for linear programming.

    Find a vector x with nonnegative elements, that maximizes
    quantity `c^T x`, subject to the constraints `m x <= b`.

    Examples
    ========

    >>> simplex([2, 3, 4], [[3, 2, 1],
    ...                     [2, 5, 3]], [10, 15])
    (20, (0, 0, 5))

    References
    ==========

    .. [1] http://mathworld.wolfram.com/SimplexMethod.html
    """

    m = Matrix(m)

    if len(c) != m.cols or len(b) != m.rows:
        raise ValueError("The dimensions doesn't match")

    # build full tableau
    tableau = zeros(m.rows + 1, m.cols + m.rows + 2)
    tableau[-1, :-1] = Matrix([[1] + [-_ for _ in c] + [0]*m.rows])
    tableau[:-1, 1:m.cols + 1] = m
    tableau[:-1, m.cols + 1:-1] = eye(m.rows)
    tableau[:, -1] = Matrix(b + [0])

    if any(_.is_negative for _ in tableau[:-1, -1]):
        raise NotImplementedError("Phase I for simplex isn't implemented.")

    # Pivoting strategy use Bland's rule

    def pivot_col(obj):
        low, idx = 0, 0
        for i in range(1, len(obj) - 1):
            if obj[i] < low:
                low, idx = obj[i], i
        return -1 if idx == 0 else idx

    def pivot_row(lhs, rhs):
        ratio, idx = oo, 0
        for i in range(len(rhs)):
            if lhs[i] > 0:
                r = rhs[i]/lhs[i]
                if r < ratio:
                    ratio, idx = r, i
        return idx

    # Now solve

    while min(tableau[-1, 1:-1]) < 0:
        col = pivot_col(tableau[-1, :])
        row = pivot_row(tableau[0:-1, col], tableau[0:-1, -1])
        if tableau[row, col] <= 0:
            return oo, (oo,)*m.cols

        tableau[row, :] /= tableau[row, col]
        for r in range(tableau.rows - 1):
            if r == row:
                continue
            tableau[r, :] -= tableau[r, col]*tableau[row, :]
        tableau[-1, :] -= tableau[-1, col]*tableau[row, :]

    ans = [S.Zero]*m.cols
    for i in range(1, m.cols + 1):
        if tableau[-1, i] == 0:
            for j in range(tableau.rows - 1):
                if tableau[j, i] == 1:
                    ans[i - 1] = tableau[j, -1]
                    break

    return tableau[-1, -1], tuple(ans)
