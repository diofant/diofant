"""Solvers of systems of polynomial equations. """

from ..domains import EX
from ..matrices import Matrix
from ..polys import groebner, poly, sring
from ..polys.polyerrors import ComputationFailed, PolificationFailed
from ..polys.polytools import parallel_poly_from_expr
from ..polys.solvers import solve_lin_sys
from ..simplify import simplify
from ..utilities import default_sort_key


__all__ = 'solve_linear_system', 'solve_poly_system'


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


def solve_poly_system(eqs, *gens, **args):
    """
    Solve a system of polynomial equations.

    Polynomial system may have finite number of solutions or
    infinitely many (positive-dimensional systems).

    References
    ==========

    * :cite:`Cox2015ideals`, p. 98

    Examples
    ========

    >>> solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y)
    [{x: 0, y: 0}, {x: 2, y: -sqrt(2)}, {x: 2, y: sqrt(2)}]

    >>> solve_poly_system([x*y], x, y)
    [{x: 0}, {y: 0}]

    """
    try:
        args['extension'] = False
        polys, opt = parallel_poly_from_expr(eqs, *gens, **args)
        polys = [p.to_exact() for p in polys]
    except PolificationFailed as exc:
        raise ComputationFailed('solve_poly_system', len(eqs), exc)

    def _solve_reduced_system(system, gens):
        """Recursively solves reduced polynomial systems."""

        basis = groebner(system, gens, polys=True, extension=False)
        dim = basis.dimension
        solutions = []

        if dim is None:
            return []

        elif dim > 0:
            max_iset = max(basis.independent_sets, key=len)
            new_gens = [g for g in gens if g not in max_iset]

            # After removing variables from the maximal set of independent
            # variables for the given ideal - the new ideal is of dimension
            # zero with the independent variables as parameters in the
            # coefficient domain.
            solutions.extend(_solve_reduced_system(system, new_gens))

            # Now we should examine cases when leading coefficient of
            # some polynomial in the system is zero.
            for p in basis.polys:
                lc = poly(p, *new_gens).LC(order=basis.order)
                for special in _solve_reduced_system(system + [lc], gens):
                    # This heuristics wipe out some redundant special
                    # solutions, which already there in solutions after
                    # solving the system with new set of generators.
                    if all(any((_.subs(s) - _).subs(special).simplify()
                               for _ in gens) for s in solutions):
                        solutions.insert(0, special)

        else:
            # By the elimination property, the last polynomial should
            # be univariate in the last variable.
            f = basis[-1]
            gen = gens[-1]

            zeros = {k.doit() for k in f.ltrim(gen).all_roots()}

            if len(basis) == 1:
                return [{gen: zero} for zero in zeros]

            new_basis = [b.set_domain(EX) for b in basis[:-1]]

            # Now substitute zeros for the last variable and
            # solve recursively new obtained zero-dimensional systems.
            for zero in zeros:
                new_system = []
                new_gens = gens[:-1]

                for b in new_basis:
                    eq = b.eval(gen, zero)

                    if not eq.is_zero:
                        new_system.append(eq)

                for solution in _solve_reduced_system(new_system, new_gens):
                    solution[gen] = zero
                    solutions.append(solution)

        return solutions

    result = _solve_reduced_system(polys, opt.gens)

    if not opt.domain.is_Exact:
        result = [{k: r[k].evalf(opt.domain.dps) for k in r} for r in result]

    return sorted(result, key=default_sort_key)
