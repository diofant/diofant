"""Tools for solving inequalities and systems of inequalities. """

from collections import defaultdict

from ..core import Dummy, Eq, Ge, Integer, Lt, S, Symbol, oo
from ..core.compatibility import iterable
from ..core.relational import Relational
from ..functions import Abs, Piecewise
from ..logic import And, false, true
from ..polys import Poly, parallel_poly_from_expr
from ..polys.polyutils import _nsort
from ..sets import FiniteSet, Interval, Reals, Union
from ..utilities import filldedent


__all__ = 'reduce_inequalities',


def solve_poly_inequality(poly, rel):
    """
    Solve a polynomial inequality with rational coefficients.

    Examples
    ========

    >>> solve_poly_inequality(Poly(x), '==')
    [{0}]
    >>> solve_poly_inequality(Poly(x**2 - 1), '!=')
    [(-oo, -1), (-1, 1), (1, oo)]
    >>> solve_poly_inequality(Poly(x**2 - 1), '==')
    [{-1}, {1}]

    See Also
    ========

    solve_poly_inequalities

    """
    if not isinstance(poly, Poly):
        raise ValueError('`poly` should be a Poly instance')
    if rel not in {'>', '<', '>=', '<=', '==', '!='}:
        raise ValueError("Invalid relational operator symbol: %r" % rel)
    if poly.is_number:
        t = Relational(poly.as_expr(), 0, rel)
        if t == true:
            return [Reals]
        elif t == false:
            return [S.EmptySet]
        else:
            raise NotImplementedError("Couldn't determine truth value of %s" % t)

    reals, intervals = poly.real_roots(multiple=False), []

    if rel == '==':
        for root, _ in reals:
            interval = Interval(root, root)
            intervals.append(interval)
    elif rel == '!=':
        left = -oo

        for right, _ in reals + [(oo, 1)]:
            interval = Interval(left, right, True, True)
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

        right, right_open = oo, True

        for left, multiplicity in reversed(reals):
            if multiplicity % 2:
                if sign == eq_sign:
                    intervals.insert(0, Interval(left, right, not equal, right_open))

                sign, right, right_open = -sign, left, not equal
            else:
                if sign == eq_sign and not equal:
                    intervals.insert(0, Interval(left, right, True, right_open))
                    right, right_open = left, True
                elif sign != eq_sign and equal:
                    intervals.insert(0, Interval(left, left))

        if sign == eq_sign:
            intervals.insert(0, Interval(-oo, right, True, right_open))

    return intervals


def solve_poly_inequalities(polys):
    """
    Solve polynomial inequalities with rational coefficients.

    Examples
    ========

    >>> solve_poly_inequalities(((Poly(x**2 - 3), ">"),
    ...                          (Poly(-x**2 + 1), ">")))
    (-oo, -sqrt(3)) U (-1, 1) U (sqrt(3), oo)

    """
    return Union(*[solve_poly_inequality(*p) for p in polys])


def solve_rational_inequalities(eqs):
    """
    Solve a system of rational inequalities with rational coefficients.

    Examples
    ========

    >>> solve_rational_inequalities([[((Poly(-x + 1), Poly(1, x)), '>='),
    ...                               ((Poly(-x + 1), Poly(1, x)), '<=')]])
    {1}

    >>> solve_rational_inequalities([[((Poly(x), Poly(1, x)), '!='),
    ...                               ((Poly(-x + 1), Poly(1, x)), '>=')]])
    (-oo, 0) U (0, 1]

    See Also
    ========

    solve_poly_inequality

    """
    result = S.EmptySet

    for eq in eqs:
        global_intervals = [Reals]

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

        for interval in global_intervals:
            result |= interval

    return result


def reduce_rational_inequalities(exprs, gen, relational=True):
    """
    Reduce a system of rational inequalities with rational coefficients.

    Examples
    ========

    >>> x = Symbol('x', real=True)

    >>> reduce_rational_inequalities([[x**2 <= 0]], x)
    Eq(x, 0)
    >>> reduce_rational_inequalities([[x + 2 > 0]], x)
    -2 < x
    >>> reduce_rational_inequalities([[(x + 2, ">")]], x)
    -2 < x
    >>> reduce_rational_inequalities([[x + 2]], x)
    Eq(x, -2)

    """
    exact = True
    eqs = []
    solution = Reals if exprs else S.EmptySet
    for _exprs in exprs:
        _eqs = []

        for expr in _exprs:
            if isinstance(expr, tuple):
                expr, rel = expr
            else:
                if expr.is_Relational:
                    expr, rel = expr.lhs - expr.rhs, expr.rel_op
                else:
                    expr, rel = expr, '=='

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

    >>> x = Symbol('x', real=True)

    >>> reduce_piecewise_inequality(Abs(x - 5) - 3, '<', x)
    (2 < x) & (x < 8)
    >>> reduce_piecewise_inequality(Abs(x + 2)*3 - 13, '<', x)
    (-19/3 < x) & (x < 7/3)

    >>> reduce_piecewise_inequality(Piecewise((1, x < 1),
    ...                                       (3, True)) - 1, '>', x)
    1 <= x

    See Also
    ========

    reduce_piecewise_inequalities

    """
    if gen.is_extended_real is False:
        raise TypeError(filldedent('''
            can't solve inequalities with piecewise
            functions containing non-real variables'''))

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

            if not n.is_Integer:  # pragma: no cover
                raise NotImplementedError("only integer powers are supported")

            _exprs = _bottom_up_scan(expr.base)

            for expr, conds in _exprs:
                exprs.append((expr**n, conds))
        elif isinstance(expr, Abs):
            _exprs = _bottom_up_scan(expr.args[0])

            for expr, conds in _exprs:
                exprs.append(( expr, conds + [Ge(expr, 0)]))
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
            expr = Relational( expr, 0, rel)
        else:
            expr = Relational(-expr, 0, mapping[rel])

        inequalities.append([expr] + conds)

    return reduce_rational_inequalities(inequalities, gen)


def reduce_piecewise_inequalities(exprs, gen):
    """
    Reduce a system of inequalities with nested piecewise functions.

    Examples
    ========

    >>> x = Symbol('x', real=True)

    >>> reduce_piecewise_inequalities([(Abs(3*x - 5) - 7, '<'),
    ...                                (Abs(x + 25) - 13, '>')], x)
    (-2/3 < x) & (x < 4) & ((-12 < x) | (x < -38))
    >>> reduce_piecewise_inequalities([(Abs(x - 4) + Abs(3*x - 5) - 7, '<')], x)
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

    >>> x = Symbol('x', real=True)

    >>> solve_univariate_inequality(x**2 >= 4, x)
    (2 <= x) | (x <= -2)
    >>> solve_univariate_inequality(x**2 >= 4, x, relational=False)
    (-oo, -2] U [2, oo)

    """
    from ..simplify import simplify
    from .solvers import solve, denoms

    e = expr.lhs - expr.rhs
    parts = n, d = e.as_numer_denom()
    if all(i.is_polynomial(gen) for i in parts):
        solns = solve(n, gen, check=False)
        singularities = solve(d, gen, check=False)
    else:
        solns = solve(e, gen, check=False)
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
        elif v.is_comparable is False:
            return False
        else:
            raise NotImplementedError

    start = -oo
    sol_sets = [S.EmptySet]
    reals = _nsort(set(solns + singularities), separated=True)[0]
    for x in reals:
        end = x

        if end in [-oo, oo]:
            if valid(Integer(0)):
                sol_sets.append(Interval(start, oo, True, True))
                break

        if valid((start + end)/2 if start != -oo else end - 1):
            sol_sets.append(Interval(start, end, True, True))

        if x in singularities:
            singularities.remove(x)
        elif include_x:
            sol_sets.append(FiniteSet(x))

        start = end

    end = oo

    if valid(start + 1):
        sol_sets.append(Interval(start, end, True, True))

    rv = Union(*sol_sets)
    return rv if not relational else rv.as_relational(gen)


def _reduce_inequalities(inequalities, symbols):
    # helper for reduce_inequalities

    poly_part, pw_part = defaultdict(list), defaultdict(list)
    other = []

    for inequality in inequalities:
        if inequality == true:
            continue
        elif inequality == false:
            return false

        expr, rel = inequality.lhs, inequality.rel_op  # rhs is 0

        # check for gens using atoms which is more strict than free_symbols to
        # guard against EX domain which won't be handled by
        # reduce_rational_inequalities
        gens = expr.atoms(Dummy, Symbol)

        if len(gens) == 1:
            gen = gens.pop()
        else:
            common = expr.free_symbols & symbols
            if len(common) == 1:
                gen = common.pop()
                other.append(solve_univariate_inequality(Relational(expr, 0, rel), gen))
                continue
            else:
                raise NotImplementedError(filldedent('''
                    inequality has more than one
                    symbol of interest'''))

        if expr.is_polynomial(gen):
            poly_part[gen].append((expr, rel))
        else:
            components = set(expr.find(lambda u: u.has(gen) and
                                       (u.is_Function or u.is_Pow and
                                        not u.exp.is_Integer)))
            if components and all(isinstance(i, Abs) or isinstance(i, Piecewise) for i in components):
                pw_part[gen].append((expr, rel))
            else:
                other.append(solve_univariate_inequality(Relational(expr, 0, rel), gen))

    poly_reduced = []
    pw_reduced = []

    for gen, exprs in poly_part.items():
        poly_reduced.append(reduce_rational_inequalities([exprs], gen))

    for gen, exprs in pw_part.items():
        pw_reduced.append(reduce_piecewise_inequalities(exprs, gen))

    return And(*(poly_reduced + pw_reduced + other))


def reduce_inequalities(inequalities, symbols=[]):
    """
    Reduces a system of inequalities or equations.

    Examples
    ========

    >>> x = Symbol('x', real=True)
    >>> y = Symbol('y', real=True)

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
        elif i == false:
            return false
        keep.append(i)
    inequalities = keep
    del keep

    gens = set().union(*[i.free_symbols for i in inequalities])

    if not iterable(symbols):
        symbols = [symbols]
    symbols = set(symbols) or gens

    # make vanilla symbol real
    recast = {i: Dummy(i.name, extended_real=True)
              for i in gens if i.is_extended_real is None}
    inequalities = [i.xreplace(recast) for i in inequalities]
    symbols = {i.xreplace(recast) for i in symbols}

    # solve system
    rv = _reduce_inequalities(inequalities, symbols)

    # restore original symbols and return
    return rv.xreplace({v: k for k, v in recast.items()})
