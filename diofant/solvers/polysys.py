"""Solvers of systems of polynomial equations."""

import collections

from ..core import Integer, diff, expand_mul
from ..domains import EX
from ..matrices import Matrix
from ..polys import (LC, LT, ComputationFailedError, PolificationFailedError,
                     degree, groebner, parallel_poly_from_expr, real_roots,
                     subresultants)
from ..polys.solvers import solve_lin_sys
from ..simplify.simplify import simplify
from ..utilities import default_sort_key, numbered_symbols
from .utils import checksol


__all__ = ('solve_linear_system', 'solve_poly_system',
           'solve_surd_system')


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
    polys, _ = parallel_poly_from_expr([expand_mul(e) for e in eqs],
                                       *symbols, field=True)
    domain = polys[0].rep.ring
    polys = [_.rep for _ in polys]

    res = solve_lin_sys(polys, domain)
    if res is None:
        return

    for k in list(res):
        s = domain.symbols[domain.index(k)]
        res[s] = domain.to_expr(res[k])
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
    except PolificationFailedError as exc:
        raise ComputationFailedError('solve_poly_system', len(eqs), exc) from exc

    def _solve_reduced_system(system, gens):
        """Recursively solves reduced polynomial systems."""
        basis = groebner(system, *gens, polys=True, extension=False)
        dim = basis.dimension
        solutions = []

        if dim is None:
            return []

        if dim > 0:
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
                lc = p.as_poly(*new_gens).LC(order=basis.order)
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

            zeros = {k.doit() for k in f.exclude().all_roots()}

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


def solve_surd_system(eqs, *gens, **args):
    """
    Solve a system of algebraic equations.

    Examples
    ========

    >>> solve_surd_system([x + sqrt(x + 1) - 2])
    [{x: -sqrt(13)/2 + 5/2}]

    """
    eqs = list(eqs)

    if not gens:
        gens = set().union(*[_.free_symbols for _ in eqs])
        gens = sorted(gens, key=default_sort_key)
    else:
        gens = list(gens)

    aux = numbered_symbols('a')
    neqs = len(eqs)
    orig_eqs = eqs[:]
    ngens = len(gens)
    bases = collections.defaultdict(dict)

    def q_surd(e):
        return e.is_Pow and e.exp.is_Rational and not e.exp.is_Integer

    def tr_surd(e):
        n, d = e.exp.as_numer_denom()
        for v2, d2 in sorted(bases.get(e.base, {}).items(),
                             key=lambda _: -_[1]):
            if not d2 % d:
                return v2**(d2 // d)
        v = next(aux)
        bases[e.base][v] = d
        gens.append(v)
        eqs.append(v**d - e.base)
        return v**n

    for i in range(neqs):
        eqs[i] = eqs[i].replace(q_surd, tr_surd)

    denoms = []
    for i, e in enumerate(eqs):
        eqs[i], d = e.as_numer_denom()
        if not d.is_constant(*gens):
            denoms.insert(0, d)

    weaksols = solve_poly_system(eqs, *gens, **args)

    for i in range(len(weaksols) - 1, -1, -1):
        if any(checksol(_, weaksols[i], warn=True) for _ in denoms):
            del weaksols[i]
        elif any(checksol(_, weaksols[i], warn=True) is False for _ in orig_eqs):
            del weaksols[i]
        else:
            for g in gens[ngens:]:
                del weaksols[i][g]

    return weaksols


def solve_poly_system_cad(seq, gens, return_one_sample=True):
    """
    Solves a system of polynomial inequalities/equalities via
    cylindrical algebraic decomposition. Returns a sample point
    in one (if return_one_sample=True) or all cells over which
    the system holds. If the system is invalid over all cells, then
    we return None.

    Parameters
    ==========

    seq: a list/tuple/set
        Listing all the (in)equalities that are needed to be solved
    gens: generators
        Generators of the (in)equalities in seq for which we want the
        solutions
    return_one_sample: bool
        If True, returns a single satisfying point. If False, returns
        a sample point from each CAD cell over which the system holds.

    Returns
    =======

    List[Dict]
        a list of dicts with the returned sample points. Each dict
        is a point, with the keys being the variables. If the system
        is unsatisfiable, then an empty list is returned.

    Examples
    ========

    >>> solve_poly_system_cad([-x**2 - 1 > 0], [x])
    []

    >>> solve_poly_system_cad([x**2 - 1 > 0], [x])
    [{x: -2}]

    >>> solve_poly_system_cad([x**2 - 1 > 0], [x], False)
    [{x: -2}, {x: 2}]

    >>> solve_poly_system_cad([y*x**2 > 0, x + y < 1], [x, y], False)
    [{x: -1, y: 1/2}, {x: 1/4, y: 1/2}, {x: -1, y: 1}, {x: -2, y: 2}]

    """
    # prepare the atoms
    atoms = [(p.lhs - p.rhs).as_poly(*gens) for p in seq]

    sample_points = cylindrical_algebraic_decomposition(atoms, gens)
    valid_samples = []

    for sample in sample_points:
        if all(expr.subs(sample) for expr in seq):
            valid_samples.append(sample)
            if return_one_sample:
                break

    return valid_samples


# HONG PROJECTOR OPERATOR AND OPERATIONS USED FOR IT


def red(f, mvar):
    """
    The reductum of a function f, as treated as a univariate function
    of the main variable (mvar), is red(f) = f - lt(f), where lt is
    the leading term.

    Parameters
    ==========

    f: Expr or Poly
        A polynomial
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    Poly
        The reductum.

    Examples
    ========

    >>> red(x**3 + x**2 + 3*x, x)
    x**2 + 3*x

    """
    return f - LT(f, mvar)


def red_set(f, mvar):
    """
    The set of reducta of a function f is defined recursively.

    The ith level reducta of f, red^i(f), is defined recursively.
    red^0(f) = f
    red^i(f) = red(red^{i-1}(f))

    The reducta set RED(f) is defined as: {red^i(f) | 0 <= i <= deg(f)}.

    This function returns RED(f). Note the ith level reductum, if
    needed, can be accessed by indexing from the reducta set.

    Parameters
    ==========

    f: Expr or Poly
        A polynomial
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    Poly
        The reducta set.

    Examples
    ========

    >>> red_set(x**3 + x**2 + 3*x, x)
    [x**3 + x**2 + 3*x, x**2 + 3*x, 3*x, 0]

    """
    reds = []

    # handle the constant case here
    # because otherwise diofant says deg= -infty
    try:
        if f.is_number:
            return []
    except AttributeError:  # if its not a diofant Basic object, then also return []
        return []

    for i in range(degree(f, mvar) + 1):
        if i == 0:
            reds.append(f)
        else:
            reds.append(red(reds[i-1], mvar))

    return reds


def subresultant_polynomials(f, g, mvar):
    """
    Computes the subresultant polynomials themselves. It uses the
    subresultant PRS which is already built into SymPy.

    Assume without loss of generality that the degree of g is
    no more than the degree of f (in this function, we gracefully
    handle the opposite case by swapping them). Then, we can compute
    the subresultant polynomials from the subresultant PRS.

    The remainder r_i is the deg(r_{i-1})-th subresultant polynomial.
    Additionally, if deg(r_i) < deg(r_{i-1}) - 1, then the deg(r_i)
    subresultant polynomial is r_i * LC(r_i) ^ c_i, where
    c_i = deg(r_{i-1})-deg(r_i)-1).

    Parameters
    ==========

    f: Expr or Poly
        A polynomial
    g: Expr or Poly
        A polynomial
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    List
        The subresultants of f and g, in increasing degree.
        The 0th element is the resultant itself.

    Examples
    ========

    >>> subresultant_polynomials(x**2, x, x)
    [0, x]

    """
    # ensure deg(f) >= deg(g)
    if degree(f, mvar) < degree(g, mvar):
        f, g = g, f

    prs = subresultants(f, g, mvar)
    if len(prs) <= 1:
        return []

    subres_polys = [Integer(0)] * (degree(g, mvar) + 1)

    for i in reversed(range(2, len(prs))):
        subres_polys[degree(prs[i-1], mvar) - 1] = prs[i]

        if degree(prs[i], mvar) < degree(prs[i-1], mvar) - 1:
            degree_jump = degree(prs[i-1], mvar) - degree(prs[i], mvar) - 1
            subres_polys[degree(prs[i], mvar)] = prs[i] * LC(prs[i], mvar)**degree_jump

    # get last one
    subres_polys[-1] = prs[1] * LC(g, mvar) ** (degree(f, mvar) - degree(g, mvar) - 1)

    # try to expand to simplify
    for i, sp in enumerate(subres_polys):
        subres_polys[i] = sp.expand()

    return subres_polys


def subresultant_coefficients(f, g, mvar):
    """
    Computes the principal subresultant coefficients (PSC). Given the
    subresultant polynomials, in increasing degree, the ith PSC is
    the coefficient of mvar^i in the ith subresultant.

    Parameters
    ==========

    f: Expr or Poly
        A polynomial
    g: Expr or Poly
        A polynomial
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    List
        The principal subresultant coefficients (PSC) of f and g.

    Examples
    ========

    >>> subresultant_coefficients(x**2, x, x)
    [0, 1]

    """
    subres_polys = subresultant_polynomials(f, g, mvar)

    subres_coeffs = []

    for i, sp in enumerate(subres_polys):
        curr_coeff = sp.as_poly(mvar).coeff_monomial((i,))
        subres_coeffs.append(curr_coeff)

    return subres_coeffs


# HONG PROJECTION OPERATOR (1990)
# input: set F of k-variate polynomials
# output: set F' of (k-1)-variate polynomials such that a CAD of R^{k-1} can be lifted to R^k

def projone(F, mvar):
    r"""
    Computes the PROJ1 operator as defined in Hong 1990.

    Let F be a set of polynomials with a given mvar. Then,

    PROJ1 = \cup_{f \in F, g \in RED(f)} (ldcf(g) \cup PSC(g, D(g)))

    where RED is the reducta set, ldcf is the leading coefficient,
    PSC is the principal subresultant coefficient set, and D is the
    derivative operator.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    Set
        A set of projection factors.

    Examples
    ========

    >>> projone([y*x**2, x + 1], x)
    {0, 1, y, 2*y}

    """
    proj_set = set()
    for f in F:
        for g in red_set(f, mvar):
            proj_set.add(LC(g, mvar))
            proj_set.update(subresultant_coefficients(g, diff(g, mvar), mvar))

    return proj_set


def projtwo(F, mvar):
    r"""
    Computes the PROJ2* operator as defined in Hong (1990).

    This is an updated version of the PROJ2 operator from Collins.
    We will just call it PROJ2 here.

    Let F be a set of polynomials with a given mvar. Then,

    PROJ2 = \cup_{f,g \in F, f < g} \cup_{f' \in RED(f)} PSC(f', g)

    where RED is the reducta set, < indicates an arbitray "linear
    ordering" to not loop over redundant pairs, and PSC is the
    principal subresultant coefficient set.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    Set
        A set of projection factors.

    Examples
    ========

    >>> projtwo([y*x**2, x + 1], x)
    {1, y}

    """
    proj_set = set()
    for i, f in enumerate(F):
        # impose "linear ordering"
        for j in range(i+1, len(F)):
            g = F[j]
            for f_ in red_set(f, mvar):
                proj_set.update(subresultant_coefficients(f_, g, mvar))

    return proj_set


def hongproj(F, mvar):
    r"""
    The Hong projection operator, as defined in Hong(1990).
    PROJH takes a set of k-variate polynomials F, with an mvar, and
    returns a set of (k-1)-variate polynomials F, with the mvar
    eliminated. These projection factors satisfy the property that a
    CAD of R^{k-1} can be lifted to a CAD of R^k.

    The Hong projector, PROJH, is defined as:

    PROJH(F) = PROJ1(F) \cup PROJ2(F)

    PROJ1 = \cup_{f \in F, g \in RED(f)} (ldcf(g) \cup PSC(g, D(g)))

    PROJ2 = \cup_{f,g \in F, f < g} \cup_{f' \in RED(f)} PSC(f', g)

    where RED is the reducta set, ldcf is the leading coefficient,
    PSC is the principal subresultant coefficient set, < indicates
    an arbitray "linear ordering" to not loop over redundant pairs,
    and D is the derivative operator.

    We remove constants, and keep polynomials that are unique up to
    sign.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    mvar: a generator
        The "main variable".
        Polynomials are treated as univariate in the mvar.

    Returns
    =======

    Set
        The set of projection factors.

    Examples
    ========

    >>> hongproj([y*x**2, x + 1], x)
    {y, 2*y}

    """
    proj_factors = projone(F, mvar).union(projtwo(F, mvar))
    proj_factors_clean = set()

    for p in proj_factors:
        if p.is_number:
            continue
        if -p not in proj_factors_clean:
            proj_factors_clean.add(p)

    return proj_factors_clean


def cylindrical_algebraic_decomposition(F, gens):
    """
    Calculates a cylindrical algebraic decomposition adapted to F.
    Uses the Hong projection operator. Returns sample points which
    represent cells over which each f in F is sign-invariant. It
    projects iteratively down to lower-dimension spaces according to
    the list of generators given in gens, in their order.

    Parameters
    ==========

    F: a list/tuple/set
        A list of polyomials
    gens: a list of generators

    Returns
    =======

    List[Dict]
        a list of dicts with the returned sample points. Each dict
        is a point, with the keys being the variables. A sample point
        is returned from every cell made by the CAD algorithm.

    Examples
    ========

    >>> cylindrical_algebraic_decomposition([x**2 - 1], [x])
    [{x: -2}, {x: -1}, {x: 0}, {x: 1}, {x: 2}]

    """
    # Compute the projection sets
    projs_set = [F]
    for i in range(len(gens) - 1):
        projs_set.append(hongproj(projs_set[-1], gens[i]))

    # Lifting
    sample_points = [{}]

    for i in reversed(range(len(gens))):
        projs = projs_set[i]
        gen = gens[i]

        new_sample_points = []

        for i, point in enumerate(sample_points):
            roots = set()
            for proj in projs:
                subbed = proj.subs(point)
                subbed = subbed.expand()
                if not subbed.is_number:
                    roots.update(real_roots(subbed))
            # have to sort them overall now
            roots = sorted(roots)

            # Calculate sample points
            if not roots:
                samples = [0]
            elif len(roots) == 1:
                samples = [roots[0] - 1, roots[0], roots[0] + 1]
            else:
                samples = [roots[0] - 1]  # Point below the smallest root
                for r1, r2 in zip(roots, roots[1:]):
                    samples.extend([r1, (r1 + r2) / 2])
                samples.extend([roots[-1], roots[-1] + 1])  # Last root and point above it

            for value in samples:
                new_point = point.copy()
                new_point[gen] = value
                new_sample_points.append(new_point)

        sample_points = new_sample_points

    return sample_points
