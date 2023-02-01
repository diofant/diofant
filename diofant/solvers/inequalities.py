"""Tools for solving inequalities and systems of inequalities."""

import collections
import itertools

from ..core import Dummy, Eq, Ge, Gt, Integer, Le, Lt, Ne, S, Symbol, oo
from ..core.compatibility import iterable
from ..core.relational import Relational
from ..functions import Abs, Max, Min, Piecewise, sign
from ..logic import And, Or, false, true
from ..matrices import Matrix, diag
from ..polys import PolificationFailedError, Poly, parallel_poly_from_expr
from ..polys.polyutils import _nsort
from ..sets import FiniteSet, Interval, Union
from ..utilities import filldedent, ordered


__all__ = 'reduce_inequalities',


def canonicalize_inequalities(eqs):
    """Canonicalize system of inequalities to have only Lt/Le."""
    eqs = set(eqs)

    # Canonicalize constraints, Ne -> pair Lt, Eq -> pair Le
    eqs |= {Lt(*e.args) for e in eqs if isinstance(e, Ne)}
    eqs |= {Lt(e.rhs, e.lhs) for e in eqs if isinstance(e, Ne)}
    eqs |= {Le(*e.args) for e in eqs if isinstance(e, Eq)}
    eqs |= {Le(e.rhs, e.lhs) for e in eqs if isinstance(e, Eq)}
    eqs -= {e for e in eqs if isinstance(e, (Ne, Eq))}

    # Gt/Ge -> Lt, Le
    eqs = {e.reversed if e.func in (Gt, Ge) else e for e in eqs}

    # Now we have only Lt/Le
    return list(ordered(e.func(e.lhs - e.rhs, 0) for e in eqs))


def fourier_motzkin(A, b, c, j):
    """
    Fourier-Motzkin elimination for `j`-th variable.

    Parameters
    ==========

    A : Matrix
        The coefficients of the system.
    b : Matrix
        The constant terms in the right hand side of relations.
    c : Matrix
        The vector of boolean elements, which determine the
        type of relation (1 for Le and 0 - for Lt).
    j : int
        The variable index.

    Example
    =======

    >>> A = Matrix([[-1, 0], [2, 4], [1, -2]])
    >>> b = Matrix([-1, 14, -1])
    >>> c = Matrix([1, 1, 1])
    >>> fourier_motzkin(A, b, c, 0)
    (Matrix([
    [0,  4],
    [0, -2]]), Matrix([
    [12],
    [-2]]), Matrix([
    [1],
    [1]]))

    References
    ==========

    * :cite:`Schrijver1998theory`, pp. 155–156.

    """
    m = A.rows
    rows = [[], [], []]
    D, d, k = [Matrix()]*3

    assert m == b.rows == c.rows
    assert all(_.is_comparable for _ in A)

    for i, a in enumerate(A[:, j]):
        rows[int(sign(a) + 1)].append(i)

    for p in itertools.chain(rows[1], itertools.product(*rows[::2])):
        if p in rows[1]:
            D = D.col_join(A[p, :])
            d = d.col_join(Matrix([b[p]]))
            k = k.col_join(Matrix([c[p]]))
        else:
            s, t = p
            D = D.col_join(A[t, j]*A[s, :] - A[s, j]*A[t, :])
            d = d.col_join(Matrix([A[t, j]*b[s] - A[s, j]*b[t]]))
            k = k.col_join(Matrix([c[s] and c[t]]))

    return D, d, k


def solve_linear_inequalities(eqs, *gens, **args):
    """
    Solve system of linear inequalities.

    Examples
    ========

    >>> solve_linear_inequalities([x >= 0, 2*x + 4*y <= 14, x - 2*y <= 1])
    (x >= 0) & (x <= 4) & (y >= x/2 - 1/2) & (y <= -x/2 + 7/2)

    """
    assert all(e.is_Relational for e in eqs)

    eqs = canonicalize_inequalities(eqs)

    polys, opt = parallel_poly_from_expr([e.lhs for e in eqs], *gens, **args)

    if not all(p.is_linear for p in polys):
        raise ValueError(f'Got non-linear inequality in {eqs}')

    gens = Matrix(opt.gens)
    A = Matrix([[p.coeff_monomial(x) for x in gens] for p in polys])
    b = Matrix([-p.coeff_monomial(1) for p in polys])
    c = Matrix([e.func is Le for e in eqs])
    res = []
    failed = []
    op_map = {(1, 1): Le, (1, 0): Lt, (0, 1): Ge, (0, 0): Gt}

    for i, g in reversed(list(enumerate(gens))):
        D, d, e = fourier_motzkin(A, b, c, i)

        if not D:
            failed.append(i)
            continue

        gens_g = gens.copy()
        gens_g[i] = 0

        for j, (r, x) in enumerate(zip(b - A*gens_g, c)):
            gc = A[j, i]
            if not gc:
                continue

            res.append(op_map[(int(gc > 0), int(x))](g, r/gc))

        A, b, c = D, d, e

    if not A.is_zero:
        i = failed.pop(0)
        g = gens[i]
        gens_g = gens.copy()
        gens_g[i] = 0
        strict = []
        non_strict = []

        for r, x in zip(diag(*A[:, i])**-1*(b - A*gens_g), c):
            if x:
                non_strict.append(r)
            else:
                strict.append(r)

        pos = int(A[0, i] > 0)
        other_op = Min if pos else Max

        if strict and non_strict:
            a, b = other_op(*non_strict), other_op(*strict)
            opn, ops = op_map[(pos, 1)], op_map[(pos, 0)]
            res.append(Or(And(opn(g, a), ops(a, b)), And(ops(g, b), opn(b, a))))
        else:
            both = non_strict + strict
            res.append(op_map[(pos, int(not strict))](g, other_op(*(both))))
    elif any(_ < 0 for _ in b):
        return false

    return And(*res)


def solve_poly_inequality(poly, rel):
    """
    Solve a polynomial inequality with rational coefficients.

    Examples
    ========

    >>> solve_poly_inequality(x.as_poly(), '==')
    [{0}]
    >>> solve_poly_inequality((x**2 - 1).as_poly(), '!=')
    [[-oo, -1), (-1, 1), (1, oo]]
    >>> solve_poly_inequality((x**2 - 1).as_poly(), '==')
    [{-1}, {1}]

    See Also
    ========

    solve_poly_inequalities

    """
    if not isinstance(poly, Poly):
        raise ValueError('`poly` should be a Poly instance')
    if rel not in {'>', '<', '>=', '<=', '==', '!='}:
        raise ValueError(f'Invalid relational operator symbol: {rel!r}')
    if poly.is_number:
        t = Relational(poly.as_expr(), 0, rel)
        if t == true:
            return [S.ExtendedReals]
        if t == false:
            return [S.EmptySet]
        raise NotImplementedError(f"Couldn't determine truth value of {t}")

    reals, intervals = poly.real_roots(multiple=False), []

    if rel == '==':
        for root, _ in reals:
            interval = Interval(root, root)
            intervals.append(interval)
    elif rel == '!=':
        left = -oo

        for right, _ in reals + [(oo, 1)]:
            interval = Interval(left, right, left.is_finite, right.is_finite)
            intervals.append(interval)
            left = right
    else:
        sign = +1 if poly.LC() > 0 else -1
        eq_sign, equal = None, False

        if rel == '>':
            eq_sign = +1
        elif rel == '<':
            eq_sign = -1
        elif rel == '>=':
            eq_sign, equal = +1, True
        else:
            eq_sign, equal = -1, True

        right, right_open = oo, False

        for left, multiplicity in reversed(reals):
            if multiplicity % 2:
                if sign == eq_sign:
                    intervals.insert(0, Interval(left, right, not equal and left.is_finite, right_open and right.is_finite))

                sign, right, right_open = -sign, left, not equal
            else:
                if sign == eq_sign and not equal:
                    intervals.insert(0, Interval(left, right, left.is_finite, right_open and right.is_finite))
                    right, right_open = left, True
                elif sign != eq_sign and equal:
                    intervals.insert(0, Interval(left, left))

        if sign == eq_sign:
            intervals.insert(0, Interval(-oo, right, False, right_open and right.is_finite))

    return intervals


def solve_poly_inequalities(polys):
    """
    Solve polynomial inequalities with rational coefficients.

    Examples
    ========

    >>> solve_poly_inequalities((((+x**2 - 3).as_poly(), '>'),
    ...                          ((-x**2 + 1).as_poly(), '>')))
    [-oo, -sqrt(3)) U (-1, 1) U (sqrt(3), oo]

    """
    return Union(*[solve_poly_inequality(*p) for p in polys])


def solve_rational_inequalities(eqs):
    """
    Solve a system of rational inequalities with rational coefficients.

    Examples
    ========

    >>> solve_rational_inequalities([[(((-x + 1).as_poly(),
    ...                                 Integer(1).as_poly(x)), '>='),
    ...                               (((-x + 1).as_poly(),
    ...                                 Integer(1).as_poly(x)), '<=')]])
    {1}

    >>> solve_rational_inequalities([[((x.as_poly(), Integer(1).as_poly(x)), '!='),
    ...                               (((-x + 1).as_poly(),
    ...                                 Integer(1).as_poly(x)), '>=')]])
    [-oo, 0) U (0, 1]

    See Also
    ========

    solve_poly_inequality

    """
    result = S.EmptySet

    for eq in eqs:
        global_intervals = [S.ExtendedReals]

        for (numer, denom), rel in eq:
            intervals = []

            for numer_interval in solve_poly_inequality(numer*denom, rel):
                for global_interval in global_intervals:
                    interval = numer_interval & global_interval

                    if interval is not S.EmptySet:
                        intervals.append(interval)

            global_intervals = intervals

            intervals = []

            for global_interval in global_intervals:
                for denom_interval in solve_poly_inequality(denom, '=='):
                    global_interval -= denom_interval

                if global_interval is not S.EmptySet:
                    intervals.append(global_interval)

            global_intervals = intervals

            if not global_intervals:
                break

            intervals = []

            for interval in global_intervals:
                intervals.append(interval)

            global_intervals = intervals

        for interval in global_intervals:
            result |= interval

    return result


def reduce_rational_inequalities(exprs, gen, relational=True):
    """
    Reduce a system of rational inequalities with rational coefficients.

    Examples
    ========

    >>> reduce_rational_inequalities([[x**2 <= 0]], x)
    Eq(x, 0)
    >>> reduce_rational_inequalities([[x + 2 > 0]], x)
    -2 < x
    >>> reduce_rational_inequalities([[(x + 2, '>')]], x)
    -2 < x
    >>> reduce_rational_inequalities([[x + 2]], x)
    Eq(x, -2)

    """
    exact = True
    eqs = []
    solution = S.ExtendedReals if exprs else S.EmptySet
    for _exprs in exprs:
        _eqs = []

        for expr in _exprs:
            if isinstance(expr, tuple):
                expr, rel = expr
            else:
                if expr.is_Relational:
                    expr, rel = expr.lhs - expr.rhs, expr.rel_op
                else:
                    rel = '=='

            if expr == true:
                numer, denom, rel = Integer(0), Integer(1), '=='
            elif expr == false:
                numer, denom, rel = Integer(1), Integer(1), '=='
            else:
                numer, denom = expr.together().as_numer_denom()

            (numer, denom), opt = parallel_poly_from_expr((numer, denom), gen)

            if not opt.domain.is_Exact:
                numer, denom, exact = numer.to_exact(), denom.to_exact(), False

            domain = opt.domain.get_exact()

            if not (domain.is_IntegerRing or domain.is_RationalField):
                expr = numer/denom
                expr = Relational(expr, 0, rel)
                solution &= solve_univariate_inequality(expr, gen, relational=False)
            else:
                _eqs.append(((numer, denom), rel))

        if _eqs:
            eqs.append(_eqs)

    if eqs:
        solution &= solve_rational_inequalities(eqs)

    if not exact:
        solution = solution.evalf()

    if relational:
        solution = solution.as_relational(gen)

    return solution


def reduce_piecewise_inequality(expr, rel, gen):
    """
    Reduce an inequality with nested piecewise functions.

    Examples
    ========

    >>> reduce_piecewise_inequality(abs(x - 5) - 3, '<', x)
    (2 < x) & (x < 8)
    >>> reduce_piecewise_inequality(abs(x + 2)*3 - 13, '<', x)
    (-19/3 < x) & (x < 7/3)

    >>> reduce_piecewise_inequality(Piecewise((1, x < 1),
    ...                                       (3, True)) - 1, '>', x)
    1 <= x

    See Also
    ========

    reduce_piecewise_inequalities

    """
    if gen.is_extended_real is False:
        raise TypeError(filldedent("""
            can't solve inequalities with piecewise
            functions containing non-real variables"""))

    def _bottom_up_scan(expr):
        exprs = []

        if expr.is_Add or expr.is_Mul:
            op = expr.func

            for arg in expr.args:
                _exprs = _bottom_up_scan(arg)

                if not exprs:
                    exprs = _exprs
                else:
                    args = []

                    for expr, conds in exprs:
                        for _expr, _conds in _exprs:
                            args.append((op(expr, _expr), conds + _conds))

                    exprs = args
        elif expr.is_Pow:
            n = expr.exp

            if not n.is_Integer:
                raise NotImplementedError('only integer powers are supported')

            _exprs = _bottom_up_scan(expr.base)

            for expr, conds in _exprs:
                exprs.append((expr**n, conds))
        elif isinstance(expr, Abs):
            _exprs = _bottom_up_scan(expr.args[0])

            for expr, conds in _exprs:
                exprs.append((+expr, conds + [Ge(expr, 0)]))
                exprs.append((-expr, conds + [Lt(expr, 0)]))
        elif isinstance(expr, Piecewise):
            for a in expr.args:
                _exprs = _bottom_up_scan(a.expr)

                for ex, conds in _exprs:
                    if a.cond != true:
                        exprs.append((ex, conds + [a.cond]))
                    else:
                        oconds = [c[1] for c in expr.args if c[1] != true]
                        exprs.append((ex, conds + [And(*[~c for c in oconds])]))
        else:
            exprs = [(expr, [])]

        return exprs

    exprs = _bottom_up_scan(expr)

    mapping = {'<': '>', '<=': '>='}
    inequalities = []

    for expr, conds in exprs:
        if rel not in mapping:
            expr = Relational(+expr, 0, rel)
        else:
            expr = Relational(-expr, 0, mapping[rel])

        inequalities.append([expr] + conds)

    return reduce_rational_inequalities(inequalities, gen)


def reduce_piecewise_inequalities(exprs, gen):
    """
    Reduce a system of inequalities with nested piecewise functions.

    Examples
    ========

    >>> reduce_piecewise_inequalities([(abs(3*x - 5) - 7, '<'),
    ...                                (abs(x + 25) - 13, '>')], x)
    (-2/3 < x) & (x < 4) & ((-12 < x) | (x < -38))
    >>> reduce_piecewise_inequalities([(abs(x - 4) + abs(3*x - 5) - 7, '<')], x)
    (1/2 < x) & (x < 4)

    See Also
    ========

    reduce_piecewise_inequality

    """
    return And(*[reduce_piecewise_inequality(expr, rel, gen)
                 for expr, rel in exprs])


def solve_univariate_inequality(expr, gen, relational=True):
    """
    Solves a real univariate inequality.

    Examples
    ========

    >>> solve_univariate_inequality(x**2 >= 4, x)
    (2 <= x) | (x <= -2)
    >>> solve_univariate_inequality(x**2 >= 4, x, relational=False)
    [-oo, -2] U [2, oo]

    """
    from ..simplify import simplify
    from .solvers import denoms, solve

    e = expr.lhs - expr.rhs
    solns = solve(e, gen)
    singularities = []
    for d in denoms(e):
        singularities.extend(solve(d, gen))
    solns = [s[gen] for s in solns]
    singularities = [s[gen] for s in singularities]

    include_x = expr.func(0, 0)

    def valid(x):
        v = e.subs({gen: x})
        try:
            r = expr.func(v, 0)
        except TypeError:
            r = false
        r = simplify(r)
        if r in (true, false):
            return r
        if v.is_comparable is False:
            return False
        raise NotImplementedError

    start = -oo
    sol_sets = [S.EmptySet]
    reals = _nsort(set(solns + singularities), separated=True)[0]
    for x in reals:
        end = x

        if end in [-oo, oo]:
            if valid(Integer(0)):
                sol_sets.append(Interval(start, oo, start in reals, end == oo))
                break

        if valid((start + end)/2 if start != -oo else end - 1):
            sol_sets.append(Interval(start, end, start.is_finite is not False, end.is_finite is not False))

        if x in singularities:
            singularities.remove(x)
        elif include_x:
            sol_sets.append(FiniteSet(x))

        start = end

    end = oo

    if valid(start + 1):
        sol_sets.append(Interval(start, end, True, end in reals))

    rv = Union(*sol_sets)

    rel_map = {Lt: Le, Gt: Ge}

    for t in [oo, -oo]:
        try:
            rel = rel_map.get(expr.func, expr.func)
            if rv.contains(t) is true and rel(e.limit(gen, t)) is false:
                rv -= FiniteSet(t)
        except TypeError:
            pass

    return rv if not relational else rv.as_relational(gen)


def _reduce_inequalities(inequalities, symbols):
    if len(symbols) > 1:
        try:
            return solve_linear_inequalities(inequalities, *symbols)
        except (PolificationFailedError, ValueError):
            pass

    rat_part = collections.defaultdict(list)
    pw_part = rat_part.copy()
    other = []

    for inequality in inequalities:
        if inequality == true:
            continue
        if inequality == false:
            return false

        expr, rel = inequality.lhs, inequality.rel_op  # rhs is 0

        # check for gens using atoms which is more strict than free_symbols to
        # guard against EX domain which won't be handled by
        # reduce_rational_inequalities
        gens = expr.atoms(Dummy, Symbol)

        if len(gens) == 1:
            gen = gens.pop()
        else:
            common = expr.free_symbols & set(symbols)
            if len(common) == 1:
                gen = common.pop()
                other.append(solve_univariate_inequality(Relational(expr, 0, rel), gen))
                continue
            raise NotImplementedError('Solving multivariate inequalities is '
                                      'implemented only for linear case yet.')

        if expr.is_rational_function(gen):
            rat_part[gen].append((expr, rel))
        else:
            components = set(expr.find(lambda u: u.has(gen) and
                                       (u.is_Function or u.is_Pow and
                                        not u.exp.is_Integer)))
            if components and all(isinstance(i, (Abs, Piecewise)) for i in components):
                pw_part[gen].append((expr, rel))
            else:
                other.append(solve_univariate_inequality(Relational(expr, 0, rel), gen))

    rat_reduced = []
    pw_reduced = []

    for gen, exprs in rat_part.items():
        rat_reduced.append(reduce_rational_inequalities([exprs], gen))

    for gen, exprs in pw_part.items():
        pw_reduced.append(reduce_piecewise_inequalities(exprs, gen))

    return And(*(rat_reduced + pw_reduced + other))


def reduce_inequalities(inequalities, symbols=[]):
    """
    Reduces a system of inequalities or equations.

    Examples
    ========

    >>> reduce_inequalities(0 <= x + 3, [])
    -3 <= x
    >>> reduce_inequalities(0 <= x + y*2 - 1, [x])
    -2*y + 1 <= x

    See Also
    ========

    diofant.solvers.solvers.solve : solve algebraic equations

    """
    if not iterable(inequalities):
        inequalities = [inequalities]

    # prefilter
    keep = []
    for i in inequalities:
        if isinstance(i, Relational):
            i = i.func(i.lhs.as_expr() - i.rhs.as_expr(), 0)
        elif i not in (True, False):
            i = Eq(i, 0)
        if i == true:
            continue
        if i == false:
            return false
        keep.append(i)
    inequalities = keep
    del keep

    gens = set().union(*[i.free_symbols for i in inequalities])

    if not iterable(symbols):
        symbols = [symbols]
    symbols = ordered(set(symbols) or gens)

    # make vanilla symbol real
    recast = {i: Dummy(i.name, extended_real=True)
              for i in gens if i.is_extended_real is None}
    inequalities = [i.xreplace(recast) for i in inequalities]
    symbols = ordered(i.xreplace(recast) for i in symbols)

    # solve system
    rv = _reduce_inequalities(inequalities, list(symbols))

    # restore original symbols and return
    return rv.xreplace({v: k for k, v in recast.items()})
