from ..core import Integer, Lt, diff, nan, oo, sympify
from ..core.compatibility import is_sequence
from ..functions import Min
from ..matrices import eye, zeros
from ..series import limit
from ..sets import Interval
from ..solvers import reduce_inequalities, solve
from ..solvers.inequalities import canonicalize_inequalities
from ..utilities import ordered
from .singularities import singularities


__all__ = 'minimize', 'maximize'


def minimize(f, *v):
    """Minimizes `f` with respect to given variables `v`.

    Examples
    ========

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

    constraints = canonicalize_inequalities(constraints)

    if dim == 1:
        x = v[0]
        if constraints:
            constraints.extend([x - oo < 0, -oo - x < 0])
            dom = reduce_inequalities(constraints, x).as_set()
        else:
            dom = Interval(-oo, oo, True, True)**len(v)
        return minimize_univariate(obj, x, dom)

    polys = [obj.as_poly(*v)] + [c.lhs.as_poly(*v) for c in constraints]
    is_polynomial = all(p is not None for p in polys)
    is_linear = is_polynomial and all(p.is_linear for p in polys)

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
        m = [([+p.coeff_monomial(x) for x in v] +
              [-p.coeff_monomial(x) for x in v])
             for p in polys[1:]]
        b = [-p.coeff_monomial(1) for p in polys[1:]]

        res, sol = simplex(c, m, b)
        res -= polys[0].coeff_monomial(1)
        sol = map(lambda x, y: x - y, sol[:dim], sol[dim:])

        return -res, dict(zip(v, sol))

    raise NotImplementedError


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
            extr[dom.end] = limit(f, x, dom.end, dir='-')
        for s in singularities(f, x):
            if s in dom:
                m = Min(limit(f, x, s), limit(f, x, s, dir='-'))
                if m == -oo:
                    return -oo, dict({x: s})
                else:
                    extr[s] = m

        for p in solve(diff(f, x), x):
            p = p[x]
            if p in dom:
                extr[p] = f.subs({x: p})
    elif dom.is_FiniteSet:
        for p in dom.args:
            extr[p] = f.subs({x: p})
    else:
        raise NotImplementedError

    if extr:
        min, point = oo, nan
        for p, fp in sorted(extr.items()):
            if fp < min:
                point, min = p, fp
        return min, dict({x: point})


class InfeasibleProblem(Exception):
    pass


def simplex(c, m, b):
    """
    Simplex algorithm for linear programming.

    Find a vector x with nonnegative elements, that maximizes
    quantity `c^T x`, subject to the constraints `m x <= b`.

    Examples
    ========

    >>> simplex([2, 3, 4], [[3, 2, 1], [2, 5, 3]], [10, 15])
    (20, (0, 0, 5))

    References
    ==========

    * Paul R. Thie, Gerard E. Keough, An Introduction to Linear
      Programming and Game Theory, Third edition, 2008, Ch. 3.

    """
    rows, cols = len(b), len(c)

    if len(m) != rows or any(len(_) != cols for _ in m):
        raise ValueError("The dimensions doesn't match")

    m = sorted(m, key=lambda v: b[m.index(v)])
    b = sorted(b)

    # build full tableau
    tableau = zeros(rows + 1, cols + rows + 1)
    tableau[-1, :-1] = [[-_ for _ in c] + [0]*rows]
    tableau[:-1, :cols] = m
    tableau[:-1, cols:-1] = eye(rows)
    tableau[:, -1] = b + [0]

    def pivot_col(obj):
        # use Bland's rule
        for i in range(len(obj) - 1):  # pragma: no branch
            if obj[i] < 0:
                return i

    def pivot_row(lhs, rhs):
        ratio, idx = oo, 0
        for i in range(len(lhs)):
            if lhs[i] > 0:
                r = rhs[i]/lhs[i]
                if r < ratio:
                    ratio, idx = r, i
        return idx

    def solve_simplex(tableau, basis, phase1=False):
        while min(tableau[-1, :-1]) < 0:
            col = pivot_col(tableau[-1, :])
            row = pivot_row(tableau[:-1 - phase1, col], tableau[:, -1])

            if tableau[row, col] <= 0:
                return 1
            else:
                basis[row] = col

            tableau[row, :] /= tableau[row, col]
            for r in range(tableau.rows):
                if r != row:
                    tableau[r, :] -= tableau[r, col]*tableau[row, :]
        return 0

    # Now solve

    neg_idx = [b.index(_) for _ in b if _ < 0]
    nneg = len(neg_idx)
    basis = list(range(cols + nneg - 1, cols + nneg + rows - 1))

    if neg_idx:
        tableau = tableau.col_insert(-1, zeros(tableau.rows, nneg))
        tableau = tableau.row_insert(tableau.cols, zeros(1, tableau.cols))
        j = tableau.cols - nneg - 1
        for i in neg_idx:
            tableau[i, :] *= -1
            tableau[i, j] = 1
            tableau[-1, :-1 - nneg] -= tableau[i, :-1 - nneg]
            tableau[-1, -1] -= tableau[i, -1]
            j += 1

        status = solve_simplex(tableau, basis, phase1=True)
        assert status == 0

        if tableau[-1, -1].is_nonzero:
            raise InfeasibleProblem

        del tableau[-1, :]
        for i in range(nneg):
            del tableau[:, -2]

        for row in [_ for _ in range(rows) if basis[_] > cols + rows - 1]:
            for col in range(tableau.cols - 1):  # pragma: no branch
                if tableau[row, col] != 0:
                    break
            basis[row] = col
            tableau[row, :] /= tableau[row, col]
            for r in range(tableau.rows):
                if r != row:
                    tableau[r, :] -= tableau[r, col]*tableau[row, :]

    status = solve_simplex(tableau, basis)
    if status == 1:
        return oo, (oo,)*cols

    ans = [Integer(0)]*cols
    for c, b in enumerate(basis):
        if b < cols:
            ans[b] = tableau[:-1, -1][c]
    return tableau[-1, -1], tuple(ans)
