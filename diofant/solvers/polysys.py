"""Solvers of systems of polynomial equations. """

import collections

from ..core import S
from ..matrices import Matrix
from ..polys import Poly, groebner, sring
from ..polys.polytools import parallel_poly_from_expr
from ..polys.polyerrors import (ComputationFailed, PolificationFailed,
                                CoercionFailed)
from ..polys.solvers import solve_lin_sys
from ..simplify import rcollect, simplify
from ..utilities import default_sort_key


__all__ = ('solve_linear_system', 'solve_poly_system')


class SolveFailed(Exception):
    """Raised when solver's conditions weren't met. """


def roots(p):
    r = collections.defaultdict(int)
    for v in p.all_roots():
        r[v] += 1
    return r


def solve_linear_system(system, *symbols, **flags):
    r"""Solve system of linear equations.

    Both under- and overdetermined systems are supported. The possible
    number of solutions is zero, one or infinite.

    Parameters
    ==========

    system : Matrix
        Nx(M+1) matrix, which means it has to be in augmented
        form.  This matrix will not be modified.
    \*symbols : list
        List of M Symbol's

    Returns
    =======

    solution: dict or None
        Respectively, this procedure will return None or
        a dictionary with solutions.  In the case of underdetermined
        systems, all arbitrary parameters are skipped.  This may
        cause a situation in which an empty dictionary is returned.
        In that case, all symbols can be assigned arbitrary values.

    Examples
    ========

    >>> from diofant.abc import x, y

    Solve the following system::

           x + 4 y ==  2
        -2 x +   y == 14

    >>> system = Matrix(( (1, 4, 2), (-2, 1, 14)))
    >>> solve_linear_system(system, x, y)
    {x: -6, y: 2}

    A degenerate system returns an empty dictionary.

    >>> system = Matrix(( (0,0,0), (0,0,0) ))
    >>> solve_linear_system(system, x, y)
    {}

    See Also
    ========

    diofant.matrices.matrices.MatrixBase.rref
    """

    eqs = system*Matrix(symbols + (-1,))
    domain, eqs = sring(eqs.transpose().tolist()[0], *symbols, field=True)

    res = solve_lin_sys(eqs, domain)
    if res is None:
        return

    for k in list(res.keys()):
        s = domain.symbols[domain.index(k)]
        res[s] = res[k].as_expr()
        del res[k]
        if flags.get('simplify', True):
            res[s] = simplify(res[s])

    return res


def solve_poly_system(seq, *gens, **args):
    """
    Solve a system of polynomial equations.

    Examples
    ========

    >>> from diofant.abc import x, y

    >>> solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y)
    [{x: 0, y: 0}, {x: 2, y: -sqrt(2)}, {x: 2, y: sqrt(2)}]
    """
    try:
        polys, opt = parallel_poly_from_expr(seq, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('solve_poly_system', len(seq), exc)

    if len(polys) == len(opt.gens) == 2:
        f, g = polys

        a, b = f.degree_list()
        c, d = g.degree_list()

        if a <= 2 and b <= 2 and c <= 2 and d <= 2:
            try:
                return solve_biquadratic(f, g, opt)
            except SolveFailed:
                pass

    return solve_generic(polys, opt)


def solve_biquadratic(f, g, opt):
    """
    Solve a system of two bivariate quadratic polynomial equations.

    Examples
    ========

    >>> from diofant.polys import Options, Poly
    >>> from diofant.abc import x, y

    >>> NewOption = Options((x, y), {'domain': 'ZZ'})

    >>> a = Poly(y**2 - 4 + x, y, x, domain='ZZ')
    >>> b = Poly(y*2 + 3*x - 7, y, x, domain='ZZ')
    >>> solve_biquadratic(a, b, NewOption)
    [{x: 1/3, y: 3}, {x: 41/27, y: 11/9}]

    >>> a = Poly(y + x**2 - 3, y, x, domain='ZZ')
    >>> b = Poly(-y + x - 4, y, x, domain='ZZ')
    >>> solve_biquadratic(a, b, NewOption)
    [{x: -sqrt(29)/2 + 7/2, y: -sqrt(29)/2 - 1/2},
     {x: sqrt(29)/2 + 7/2, y: -1/2 + sqrt(29)/2}]
    """
    G = groebner([f, g])

    if len(G) == 1 and G[0].is_ground:
        return []

    if len(G) != 2:
        raise SolveFailed

    p, q = G
    x, y = opt.gens

    p = Poly(p, x, expand=False)
    q = q.ltrim(-1)

    p_roots = [rcollect(expr, y) for expr in roots(p).keys()]
    q_roots = list(roots(q).keys())

    solutions = []

    for q_root in q_roots:
        for p_root in p_roots:
            solution = {x: p_root.subs(y, q_root), y: q_root}
            solutions.append(solution)

    return sorted(solutions, key=default_sort_key)


def solve_generic(polys, opt):
    """
    Solve a generic system of polynomial equations.

    Returns all possible solutions over C[x_1, x_2, ..., x_m] of a
    set F = { f_1, f_2, ..., f_n } of polynomial equations,  using
    Groebner basis approach. For now only zero-dimensional systems
    are supported, which means F can have at most a finite number
    of solutions.

    The algorithm works by the fact that, supposing G is the basis
    of F with respect to an elimination order  (here lexicographic
    order is used), G and F generate the same ideal, they have the
    same set of solutions. By the elimination property,  if G is a
    reduced, zero-dimensional Groebner basis, then there exists an
    univariate polynomial in G (in its last variable). This can be
    solved by computing its roots. Substituting all computed roots
    for the last (eliminated) variable in other elements of G, new
    polynomial system is generated. Applying the above procedure
    recursively, a finite number of solutions can be found.

    References
    ==========

    .. [Buchberger01] B. Buchberger, Groebner Bases: A Short
    Introduction for Systems Theorists, In: R. Moreno-Diaz,
    B. Buchberger, J.L. Freire, Proceedings of EUROCAST'01,
    February, 2001

    .. [Cox97] D. Cox, J. Little, D. O'Shea, Ideals, Varieties
    and Algorithms, Springer, Second Edition, 1997, pp. 112

    Examples
    ========

    >>> from diofant.polys import Poly, Options
    >>> from diofant.abc import x, y
    >>> NewOption = Options((x, y), {'domain': 'ZZ'})

    >>> a = Poly(x - y + 5, x, y, domain='ZZ')
    >>> b = Poly(x + y - 3, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [{x: -1, y: 4}]

    >>> a = Poly(x - 2*y + 5, x, y, domain='ZZ')
    >>> b = Poly(2*x - y - 3, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [{x: 11/3, y: 13/3}]

    >>> a = Poly(x**2 + y, x, y, domain='ZZ')
    >>> b = Poly(x + y*4, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [{x: 0, y: 0}, {x: 1/4, y: -1/16}]
    """

    def _solve_reduced_system(system, gens):
        """Recursively solves reduced polynomial systems. """

        basis = groebner(system, gens, polys=True)

        if len(basis) == 1 and basis[0].is_ground:
            return []

        if not basis.is_zero_dimensional:
            raise NotImplementedError("only zero-dimensional systems "
                                      "supported (finite number of solutions)")

        f = basis[-1]
        gens = f.gens
        gen = gens[-1]

        zeros = [k.doit() for k in roots(f.ltrim(gen)).keys()]

        if len(basis) == 1:
            return [{gen: zero} for zero in zeros]

        solutions = []

        for zero in zeros:
            new_system = []
            new_gens = gens[:-1]

            for b in basis[:-1]:
                eq = b.eval(gen, zero)

                if not eq.is_zero:
                    new_system.append(eq)

            for solution in _solve_reduced_system(new_system, new_gens):
                solution[gen] = zero
                solutions.append(solution)

        return solutions

    try:
        result = _solve_reduced_system(polys, opt.gens)
    except CoercionFailed:  # pragma: no cover
        raise NotImplementedError

    return sorted(result, key=default_sort_key)
