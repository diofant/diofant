"""Solvers of systems of polynomial equations. """

from diofant.core import S
from diofant.polys import Poly, groebner, roots
from diofant.polys.polytools import parallel_poly_from_expr
from diofant.polys.polyerrors import (ComputationFailed,
                                      PolificationFailed, CoercionFailed)
from diofant.simplify import rcollect
from diofant.utilities import default_sort_key, postfixes


class SolveFailed(Exception):
    """Raised when solver's conditions weren't met. """


def solve_poly_system(seq, *gens, **args):
    """
    Solve a system of polynomial equations.

    Examples
    ========

    >>> from diofant import solve_poly_system
    >>> from diofant.abc import x, y

    >>> solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y)
    [(0, 0), (2, -sqrt(2)), (2, sqrt(2))]

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
        return

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

    The ability of finding all solutions by this procedure depends
    on the root finding algorithms. If no solutions were found, it
    means only that roots() failed, but the system is solvable. To
    overcome this difficulty use numerical algorithms instead.

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
    >>> from diofant.solvers.polysys import solve_generic
    >>> from diofant.abc import x, y
    >>> NewOption = Options((x, y), {'domain': 'ZZ'})

    >>> a = Poly(x - y + 5, x, y, domain='ZZ')
    >>> b = Poly(x + y - 3, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [(-1, 4)]

    >>> a = Poly(x - 2*y + 5, x, y, domain='ZZ')
    >>> b = Poly(2*x - y - 3, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [(11/3, 13/3)]

    >>> a = Poly(x**2 + y, x, y, domain='ZZ')
    >>> b = Poly(x + y*4, x, y, domain='ZZ')
    >>> solve_generic([a, b], NewOption)
    [(0, 0), (1/4, -1/16)]
    """
    def _is_univariate(f):
        """Returns True if 'f' is univariate in its last variable. """
        for monom in f.monoms():
            if any(m > 0 for m in monom[:-1]):
                return False

        return True

    def _subs_root(f, gen, zero):
        """Replace generator with a root so that the result is nice. """
        p = f.as_expr({gen: zero})

        if f.degree(gen) >= 2:
            p = p.expand(deep=False)

        return p

    def _solve_reduced_system(system, gens):
        """Recursively solves reduced polynomial systems. """

        basis = groebner(system, gens, polys=True)

        if len(basis) == 1 and basis[0].is_ground:
            return

        univariate = list(filter(_is_univariate, basis))

        if len(univariate) == 1:
            f = univariate.pop()
        else:
            raise NotImplementedError("only zero-dimensional systems "
                                      "supported (finite number of solutions)")

        gens = f.gens
        gen = gens[-1]

        zeros = list(roots(f.ltrim(gen)).keys())

        if len(basis) == 1:
            return [ (zero,) for zero in zeros ]

        solutions = []

        for zero in zeros:
            new_system = []
            new_gens = gens[:-1]

            for b in basis[:-1]:
                eq = _subs_root(b, gen, zero)

                if eq is not S.Zero:
                    new_system.append(eq)

            for solution in _solve_reduced_system(new_system, new_gens):
                solutions.append(solution + (zero,))

        return solutions

    try:
        result = _solve_reduced_system(polys, opt.gens)
    except CoercionFailed:  # pragma: no cover
        raise NotImplementedError

    if result is not None:
        return sorted(result, key=default_sort_key)
    else:
        return
