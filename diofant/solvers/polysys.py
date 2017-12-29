"""Solvers of systems of polynomial equations. """

from ..matrices import Matrix
from ..polys import groebner, sring
from ..polys.polyerrors import (CoercionFailed, ComputationFailed,
                                PolificationFailed)
from ..polys.polytools import parallel_poly_from_expr
from ..polys.solvers import solve_lin_sys
from ..simplify import simplify
from ..utilities import default_sort_key, subsets


__all__ = ('solve_linear_system', 'solve_poly_system')


class SolveFailed(Exception):
    """Raised when solver's conditions weren't met. """


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

    >>> system = Matrix(((1, 4, 2), (-2, 1, 14)))
    >>> solve_linear_system(system, x, y)
    {x: -6, y: 2}

    A degenerate system returns an empty dictionary.

    >>> system = Matrix(((0, 0, 0), (0, 0, 0)))
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

    for k in list(res):
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

    return solve_generic(polys, opt)


def solve_generic(polys, opt):
    """
    Solve a generic system of polynomial equations.

    Returns all possible solutions over C[x_1, x_2, ..., x_m] of a
    set F = { f_1, f_2, ..., f_n } of polynomial equations,  using
    Gröbner basis approach. For now only zero-dimensional systems
    are supported, which means F can have at most a finite number
    of solutions.

    The algorithm works by the fact that, supposing G is the basis
    of F with respect to an elimination order  (here lexicographic
    order is used), G and F generate the same ideal, they have the
    same set of solutions. By the elimination property,  if G is a
    reduced, zero-dimensional Gröbner basis, then there exists an
    univariate polynomial in G (in its last variable). This can be
    solved by computing its roots. Substituting all computed roots
    for the last (eliminated) variable in other elements of G, new
    polynomial system is generated. Applying the above procedure
    recursively, a finite number of solutions can be found.

    References
    ==========

    .. [Buchberger01] B. Buchberger, Gröbner Bases: A Short
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
            solutions = []
            solved_syms = set()
            for syms in subsets(gens, min(len(system), len(basis))):
                res = _solve_reduced_system(system, syms)
                for r in res:
                    if not any(solved_syms & v.free_symbols
                               for v in r.values()):
                        solved_syms.update(syms)
                        solutions.append(r)
            return solutions

        f = basis[-1]
        gens = f.gens
        gen = gens[-1]

        zeros = {k.doit() for k in f.ltrim(gen).all_roots()}

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
