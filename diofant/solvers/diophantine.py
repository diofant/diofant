import math

from ..core import (Add, Eq, Integer, Rational, Symbol, factor_terms,
                    integer_nthroot, oo, symbols, sympify)
from ..core.assumptions import check_assumptions
from ..core.compatibility import as_int, is_sequence
from ..core.function import _mexpand
from ..core.numbers import igcdex
from ..functions import floor, sign, sqrt
from ..matrices import Matrix
from ..ntheory import (divisors, factorint, is_square, isprime, multiplicity,
                       nextprime, perfect_power, sqrt_mod, square_factor)
from ..polys import GeneratorsNeeded, factor_list
from ..simplify import signsimp
from ..utilities import default_sort_key, filldedent, numbered_symbols
from .solvers import solve


__all__ = 'diophantine', 'classify_diop'

# these types are known (but not necessarily handled)
diop_known = {
    'binary_quadratic',
    'cubic_thue',
    'general_pythagorean',
    'general_sum_of_even_powers',
    'general_sum_of_squares',
    'homogeneous_general_quadratic',
    'homogeneous_ternary_quadratic',
    'homogeneous_ternary_quadratic_normal',
    'inhomogeneous_general_quadratic',
    'inhomogeneous_ternary_quadratic',
    'linear',
    'univariate'}


def _is_int(i):
    try:
        as_int(i)
        return True
    except ValueError:
        pass


def _sorted_tuple(*i):
    return tuple(sorted(i))


def _remove_gcd(*x):
    try:
        g = math.gcd(*x)
        return tuple(i//g for i in x)
    except (TypeError, ValueError):
        return x


def _rational_pq(a, b):
    # return `(numer, denom)` for a/b; sign in numer and gcd removed
    return _remove_gcd(sign(b)*a, abs(b))


def _nint_or_floor(p, q):
    # return nearest int to p/q; in case of tie return floor(p/q)
    w, r = divmod(p, q)
    if abs(r) <= abs(q)//2:
        return w
    return w + 1


def _odd(i):
    return i % 2 != 0


def _even(i):
    return i % 2 == 0


def diophantine(eq, param=symbols('t', integer=True), syms=None):
    """
    Simplify the solution procedure of diophantine equation ``eq`` by
    converting it into a product of terms which should equal zero.

    `(x + y)(x - y) = 0` and `x + y = 0` and `x - y = 0` are solved
    independently and combined. Each term is solved by calling
    ``diop_solve()``.

    Output of ``diophantine()`` is a set of tuples. The elements of the
    tuple are the solutions for each variable in the the equation and
    are arranged according to the alphabetic ordering of the variables.
    e.g. For an equation with two variables, `a` and `b`, the first
    element of the tuple is the solution for `a` and the second for `b`.

    Parameters
    ==========

    eq : Relational or Expr
        an equation (to be solved)
    t : Symbol, optional
        the parameter to be used in the solution.
    syms : list of Symbol's, optional
        which determines the order of the elements in the returned tuple.

    Examples
    ========

    >>> diophantine(x**2 - y**2)
    {(t_0, -t_0), (t_0, t_0)}

    >>> diophantine(x*(2*x + 3*y - z))
    {(0, n1, n2), (t_0, t_1, 2*t_0 + 3*t_1)}
    >>> diophantine(x**2 + 3*x*y + 4*x)
    {(0, n1), (3*t_0 - 4, -t_0)}

    See Also
    ========

    diofant.solvers.diophantine.diop_solve

    """
    if isinstance(eq, Eq):
        eq = eq.lhs - eq.rhs

    if eq == 0:
        return {(param,)}

    try:
        var = list(eq.expand(force=True).free_symbols)
        var.sort(key=default_sort_key)
        if syms:
            if not is_sequence(syms):
                raise TypeError('syms should be given as a sequence, e.g. a list')
            syms = [i for i in syms if i in var]
            if syms != var:
                map = dict(zip(syms, range(len(syms))))
                return {tuple(t[map[i]] for i in var)
                        for t in diophantine(eq, param)}
        n, d = eq.as_numer_denom()
        if not n.free_symbols:
            return set()
        if d.free_symbols:
            dsol = diophantine(d)
            good = diophantine(n) - dsol
            return {s for s in good if _mexpand(d.subs(zip(var, s)))}
        else:
            eq = n
        eq = factor_terms(eq)
        assert not eq.is_number
        eq = eq.as_independent(*var, as_Add=False)[1]
        p = eq.as_poly()
        assert not any(g.is_number for g in p.gens)
        eq = p.as_expr()
        assert eq.is_polynomial()
    except (GeneratorsNeeded, AssertionError, AttributeError) as exc:
        raise TypeError('Equation should be a polynomial with '
                        'Rational coefficients.') from exc

    try:
        # if we know that factoring should not be attempted, skip
        # the factoring step
        *_, t = classify_diop(eq)
        if t == 'general_sum_of_squares':
            # trying to factor such expressions will sometimes hang
            terms = [(eq, 1)]
        else:
            raise TypeError
    except (TypeError, NotImplementedError):
        terms = factor_list(eq)[1]

    sols = set()

    for term in terms:

        base, _ = term
        var_t, _, eq_type = classify_diop(base, _dict=False)
        _, base = signsimp(base, evaluate=False).as_coeff_Mul()
        solution = diop_solve(base, param)

        if eq_type in ['linear', 'homogeneous_ternary_quadratic',
                       'homogeneous_ternary_quadratic_normal',
                       'general_pythagorean']:
            sols.add(merge_solution(var, var_t, solution))

        elif eq_type in ['binary_quadratic', 'general_sum_of_squares',
                         'general_sum_of_even_powers', 'univariate']:
            for sol in solution:
                sols.add(merge_solution(var, var_t, sol))

        else:
            raise NotImplementedError(f'unhandled type: {eq_type}')

    # remove null merge results
    if () in sols:
        sols.remove(())

    null = tuple([0]*len(var))
    # if there is no solution, return trivial solution
    if not sols and eq.subs(zip(var, null)) == 0:
        sols.add(null)

    return {sympify(i) for i in sols}


def merge_solution(var, var_t, solution):
    """
    This is used to construct the full solution from the solutions of sub
    equations.

    For example when solving the equation `(x - y)(x^2 + y^2 - z^2) = 0`,
    solutions for each of the equations `x - y = 0` and `x^2 + y^2 - z^2` are
    found independently. Solutions for `x - y = 0` are `(x, y) = (t, t)`. But
    we should introduce a value for z when we output the solution for the
    original equation. This function converts `(t, t)` into `(t, t, n_{1})`
    where `n_{1}` is an integer parameter.

    """
    sol = []

    if None in solution:
        return ()

    solution = iter(solution)
    params = numbered_symbols('n', integer=True, start=1)
    for v in var:
        if v in var_t:
            sol.append(next(solution))
        else:
            sol.append(next(params))

    for val, symb in zip(sol, var):
        if check_assumptions(val, **symb._assumptions) is False:
            return ()

    return tuple(sol)


def diop_solve(eq, param=symbols('t', integer=True)):
    """
    Solves the diophantine equation ``eq``.

    Unlike ``diophantine()``, factoring of ``eq`` is not attempted. Uses
    ``classify_diop()`` to determine the type of the equation and calls
    the appropriate solver function.

    Parameters
    ==========

    eq : Expr
       an expression, which is assumed to be zero.
    t : Symbol, optional
       a parameter, to be used in the solution.

    Examples
    ========

    >>> from diofant.abc import w
    >>> diop_solve(2*x + 3*y - 5)
    (3*t_0 - 5, -2*t_0 + 5)
    >>> diop_solve(4*x + 3*y - 4*z + 5)
    (t_0, 8*t_0 + 4*t_1 + 5, 7*t_0 + 3*t_1 + 5)
    >>> diop_solve(x + 3*y - 4*z + w - 6)
    (t_0, t_0 + t_1, 6*t_0 + 5*t_1 + 4*t_2 - 6, 5*t_0 + 4*t_1 + 3*t_2 - 6)
    >>> diop_solve(x**2 + y**2 - 5)
    {(-2, -1), (-2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (2, -1), (2, 1)}

    See Also
    ========

    diofant.solvers.diophantine.diophantine

    """
    var, coeff, eq_type = classify_diop(eq, _dict=False)

    if eq_type == 'linear':
        return _diop_linear(var, coeff, param)

    elif eq_type == 'binary_quadratic':
        return _diop_quadratic(var, coeff, param)

    elif eq_type == 'homogeneous_ternary_quadratic':
        x_0, y_0, z_0 = _diop_ternary_quadratic(var, coeff)
        return _parametrize_ternary_quadratic((x_0, y_0, z_0), var, coeff)

    elif eq_type == 'homogeneous_ternary_quadratic_normal':
        x_0, y_0, z_0 = _diop_ternary_quadratic_normal(var, coeff)
        return _parametrize_ternary_quadratic((x_0, y_0, z_0), var, coeff)

    elif eq_type == 'general_pythagorean':
        return _diop_general_pythagorean(var, coeff, param)

    elif eq_type == 'univariate':
        l = solve(eq)
        s = set()

        for soln in l:
            soln = list(soln.values())[0]
            if isinstance(soln, Integer):
                s.add((soln,))
        return s

    elif eq_type == 'general_sum_of_squares':
        return _diop_general_sum_of_squares(var, -int(coeff[1]), limit=oo)

    elif eq_type == 'general_sum_of_even_powers':
        for k in coeff:
            if k.is_Pow and coeff[k]:
                p = k.exp
        return _diop_general_sum_of_even_powers(var, p, -int(coeff[1]), limit=oo)

    raise NotImplementedError(f'No solver has been written for {eq_type}.')


def classify_diop(eq, _dict=True):
    try:
        var = list(eq.free_symbols)
        assert var
    except (AttributeError, AssertionError) as exc:
        raise ValueError('equation should have 1 or more free symbols') from exc
    var.sort(key=default_sort_key)
    eq = eq.expand(force=True)
    coeff = eq.as_coefficients_dict()
    if not all(_is_int(c) for c in coeff.values()):
        raise TypeError('Coefficients should be Integers')

    diop_type = None
    total_degree = eq.as_poly().total_degree()
    homogeneous = 1 not in coeff
    if total_degree == 1:
        diop_type = 'linear'

    elif len(var) == 1:
        diop_type = 'univariate'

    elif total_degree == 2 and len(var) == 2:
        diop_type = 'binary_quadratic'

    elif total_degree == 2 and len(var) == 3 and homogeneous:
        if set(coeff) & set(var):
            diop_type = 'inhomogeneous_ternary_quadratic'
        else:
            nonzero = [k for k in coeff if coeff[k]]
            if len(nonzero) == 3 and all(i**2 in nonzero for i in var):
                diop_type = 'homogeneous_ternary_quadratic_normal'
            else:
                diop_type = 'homogeneous_ternary_quadratic'

    elif total_degree == 2 and len(var) >= 3:
        if set(coeff) & set(var):
            diop_type = 'inhomogeneous_general_quadratic'
        else:
            # there may be Pow keys like x**2 or Mul keys like x*y
            if any(k.is_Mul for k in coeff):  # cross terms
                if not homogeneous:
                    diop_type = 'inhomogeneous_general_quadratic'
                else:
                    diop_type = 'homogeneous_general_quadratic'
            else:  # all squares: x**2 + y**2 + ... + constant
                if all(coeff[k] == 1 for k in coeff if k != 1):
                    diop_type = 'general_sum_of_squares'
                elif all(is_square(abs(coeff[k])) for k in coeff):
                    if abs(sum(sign(coeff[k]) for k in coeff)) == len(var) - 2:
                        # all but one has the same sign
                        # e.g. 4*x**2 + y**2 - 4*z**2
                        diop_type = 'general_pythagorean'

    elif total_degree == 3 and len(var) == 2:
        diop_type = 'cubic_thue'

    elif (total_degree > 3 and total_degree % 2 == 0 and
          all(k.is_Pow and k.exp == total_degree for k in coeff if k != 1)):
        if all(coeff[k] == 1 for k in coeff if k != 1):
            diop_type = 'general_sum_of_even_powers'

    if diop_type is not None:
        return var, dict(coeff) if _dict else coeff, diop_type

    # new diop type instructions
    # --------------------------
    # if this error raises and the equation *can* be classified,
    #  * it should be identified in the if-block above
    #  * the type should be added to the diop_known
    # if a solver can be written for it,
    #  * a dedicated handler should be written (e.g. diop_linear)
    #  * it should be passed to that handler in diop_solve
    raise NotImplementedError(filldedent("""
        This equation is not yet recognized or else has not been
        simplified sufficiently to put it in a form recognized by
        diop_classify()."""))


classify_diop.__doc__ = """
    Helper routine used by diop_solve() to find the type of the ``eq`` etc.

    Parameters
    ==========

    eq : Expr
        an expression, which is assumed to be zero.

    Examples
    ========

    >>> classify_diop(4*x + 6*y - 4)
    ([x, y], {1: -4, x: 4, y: 6}, 'linear')
    >>> classify_diop(x + 3*y - 4*z + 5)
    ([x, y, z], {1: 5, x: 1, y: 3, z: -4}, 'linear')
    >>> classify_diop(x**2 + y**2 - x*y + x + 5)
    ([x, y], {1: 5, x: 1, x**2: 1, y**2: 1, x*y: -1}, 'binary_quadratic')

    Returns
    =======

    Returns a tuple containing the type of the diophantine equation along with
    the variables(free symbols) and their coefficients. Variables are returned
    as a list and coefficients are returned as a dict with the key being the
    respective term and the constant term is keyed to Integer(1).  The type
    is one of the following:

        * %s
""" % ('\n        * '.join(sorted(diop_known)))  # noqa: SFS101


def diop_linear(eq, param=symbols('t', integer=True)):
    """
    Solves linear diophantine equations.

    A linear diophantine equation is an equation of the form `a_{1}x_{1} +
    a_{2}x_{2} + .. + a_{n}x_{n} = 0` where `a_{1}, a_{2}, ..a_{n}` are
    integer constants and `x_{1}, x_{2}, ..x_{n}` are integer variables.

    Parameters
    ==========

    eq : Expr
       is a linear diophantine equation which is assumed to be zero.
    param : Symbol, optional
       is the parameter to be used in the solution.

    Examples
    ========

    >>> diop_linear(2*x - 3*y - 5)
    (3*t_0 - 5, 2*t_0 - 5)

    Here x = -3*t_0 - 5 and y = -2*t_0 - 5

    >>> diop_linear(2*x - 3*y - 4*z - 3)
    (t_0, 2*t_0 + 4*t_1 + 3, -t_0 - 3*t_1 - 3)

    See Also
    ========

    diofant.solvers.diophantine.diop_quadratic
    diofant.solvers.diophantine.diop_ternary_quadratic
    diofant.solvers.diophantine.diop_general_pythagorean
    diofant.solvers.diophantine.diop_general_sum_of_squares

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type == 'linear':
        return _diop_linear(var, coeff, param)


def _diop_linear(var, coeff, param):
    """
    Solves diophantine equations of the form:

    a_0*x_0 + a_1*x_1 + ... + a_n*x_n == c

    Note that no solution exists if gcd(a_0, ..., a_n) doesn't divide c.

    """
    if 1 in coeff:
        # negate coeff[] because input is of the form: ax + by + c ==  0
        #                              but is used as: ax + by     == -c
        c = -coeff[1]
    else:
        c = Integer(0)

    # Some solutions will have multiple free variables in their solutions.
    if param is None:
        params = [symbols('t')]*len(var)
    else:
        temp = str(param) + '_%i'
        params = [symbols(temp % i, integer=True) for i in range(len(var))]

    if len(var) == 1:
        q, r = divmod(c, coeff[var[0]])
        if not r:
            return q,
        else:
            return None,

    # base_solution_linear() can solve diophantine equations of the form:
    #
    # a*x + b*y == c
    #
    # We break down multivariate linear diophantine equations into a
    # series of bivariate linear diophantine equations which can then
    # be solved individually by base_solution_linear().
    #
    # Consider the following:
    #
    # a_0*x_0 + a_1*x_1 + a_2*x_2 == c
    #
    # which can be re-written as:
    #
    # a_0*x_0 + g_0*y_0 == c
    #
    # where
    #
    # g_0 == gcd(a_1, a_2)
    #
    # and
    #
    # y == (a_1*x_1)/g_0 + (a_2*x_2)/g_0
    #
    # This leaves us with two binary linear diophantine equations.
    # For the first equation:
    #
    # a == a_0
    # b == g_0
    # c == c
    #
    # For the second:
    #
    # a == a_1/g_0
    # b == a_2/g_0
    # c == the solution we find for y_0 in the first equation.
    #
    # The arrays A and B are the arrays of integers used for
    # 'a' and 'b' in each of the n-1 bivariate equations we solve.

    A = [coeff[v] for v in var]
    B = []
    if len(var) > 2:
        B.append(math.gcd(A[-2], A[-1]))
        A[-2] = A[-2] // B[0]
        A[-1] = A[-1] // B[0]
        for i in range(len(A) - 3, 0, -1):
            gcd = math.gcd(B[0], A[i])
            B[0] = B[0] // gcd
            A[i] = A[i] // gcd
            B.insert(0, gcd)
    B.append(A[-1])

    # Consider the trivariate linear equation:
    #
    # 4*x_0 + 6*x_1 + 3*x_2 == 2
    #
    # This can be re-written as:
    #
    # 4*x_0 + 3*y_0 == 2
    #
    # where
    #
    # y_0 == 2*x_1 + x_2
    # (Note that gcd(3, 6) == 3)
    #
    # The complete integral solution to this equation is:
    #
    # x_0 ==  2 + 3*t_0
    # y_0 == -2 - 4*t_0
    #
    # where 't_0' is any integer.
    #
    # Now that we have a solution for 'x_0', find 'x_1' and 'x_2':
    #
    # 2*x_1 + x_2 == -2 - 4*t_0
    #
    # We can then solve for '-2' and '-4' independently,
    # and combine the results:
    #
    # 2*x_1a + x_2a == -2
    # x_1a == 0 + t_0
    # x_2a == -2 - 2*t_0
    #
    # 2*x_1b + x_2b == -4*t_0
    # x_1b == 0*t_0 + t_1
    # x_2b == -4*t_0 - 2*t_1
    #
    # ==>
    #
    # x_1 == t_0 + t_1
    # x_2 == -2 - 6*t_0 - 2*t_1
    #
    # where 't_0' and 't_1' are any integers.
    #
    # Note that:
    #
    # 4*(2 + 3*t_0) + 6*(t_0 + t_1) + 3*(-2 - 6*t_0 - 2*t_1) == 2
    #
    # for any integral values of 't_0', 't_1'; as required.
    #
    # This method is generalized for many variables, below.

    solutions = []
    for i, Bi in enumerate(B):
        tot_x, tot_y = [], []

        for arg in Add.make_args(c):
            if arg.is_Integer:
                # example: 5 -> k = 5
                k, p = arg, Integer(1)
                pnew = params[0]
            else:  # arg is a Mul or Symbol
                # example: 3*t_1 -> k = 3
                # example: t_0 -> k = 1
                k, p = arg.as_coeff_Mul()
                pnew = params[params.index(p) + 1]

            sol = sol_x, sol_y = base_solution_linear(k, A[i], Bi, pnew)

            if p == 1:
                if None in sol:
                    return tuple([None]*len(var))
            else:
                # convert a + b*pnew -> a*p + b*pnew
                if isinstance(sol_x, Add):
                    sol_x = sol_x.args[0]*p + sol_x.args[1]
                if isinstance(sol_y, Add):
                    sol_y = sol_y.args[0]*p + sol_y.args[1]

            tot_x.append(sol_x)
            tot_y.append(sol_y)

        solutions.append(Add(*tot_x))
        c = Add(*tot_y)

    solutions.append(c)
    if param is None:
        # just keep the additive constant (i.e. replace t with 0)
        solutions = [i.as_coeff_Add()[0] for i in solutions]
    return tuple(solutions)


def base_solution_linear(c, a, b, t=None):
    """
    Return the base solution for the linear equation, `ax + by = c`.

    Used by ``diop_linear()`` to find the base solution of a linear
    Diophantine equation. If ``t`` is given then the parametrized solution is
    returned.

    ``base_solution_linear(c, a, b, t)``: ``a``, ``b``, ``c`` are coefficients
    in `ax + by = c` and ``t`` is the parameter to be used in the solution.

    Examples
    ========

    >>> base_solution_linear(5, 2, 3)  # equation 2*x + 3*y = 5
    (-5, 5)
    >>> base_solution_linear(0, 5, 7)  # equation 5*x + 7*y = 0
    (0, 0)
    >>> base_solution_linear(5, 2, 3, t)  # equation 2*x + 3*y = 5
    (3*t - 5, -2*t + 5)
    >>> base_solution_linear(0, 5, 7, t)  # equation 5*x + 7*y = 0
    (7*t, -5*t)

    """
    a, b, c = _remove_gcd(a, b, c)

    if c == 0:
        if t is not None:
            if b < 0:
                t = -t
            return b*t, -a*t
        else:
            return 0, 0
    else:
        x0, y0, d = igcdex(abs(a), abs(b))

        x0 *= sign(a)
        y0 *= sign(b)

        if divisible(c, d):
            if t is not None:
                if b < 0:
                    t = -t
                return c*x0 + b*t, c*y0 - a*t
            else:
                return c*x0, c*y0
        else:
            return None, None


def divisible(a, b):
    """Returns `True` if ``a`` is divisible by ``b`` and `False` otherwise."""
    return not a % b


def diop_quadratic(eq, param=symbols('t', integer=True)):
    """
    Solves quadratic diophantine equations.

    i.e. equations of the form `Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0`. Returns a
    set containing the tuples `(x, y)` which contains the solutions. If there
    are no solutions then `(None, None)` is returned.

    Parameters
    ==========

    eq : Expr
        should be a quadratic bivariate expression which is assumed
        to be zero.
    param : Symbol, optional
        is a parameter to be used in the solution.

    Examples
    ========

    >>> diop_quadratic(x**2 + y**2 + 2*x + 2*y + 2, t)
    {(-1, -1)}

    References
    ==========

    * Methods to solve Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0, [online],
      Available: https://web.archive.org/web/20181231080858/https://www.alpertron.com.ar/METHODS.HTM
    * Solving the equation ax^2+ bxy + cy^2 + dx + ey + f= 0, [online],
      Available: https://web.archive.org/web/20180831180321/http://www.jpr2718.org/ax2p.pdf

    See Also
    ========

    diofant.solvers.diophantine.diop_linear
    diofant.solvers.diophantine.diop_ternary_quadratic
    diofant.solvers.diophantine.diop_general_sum_of_squares
    diofant.solvers.diophantine.diop_general_pythagorean

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type == 'binary_quadratic':
        return _diop_quadratic(var, coeff, param)


def _diop_quadratic(var, coeff, t):

    x, y = var

    A = coeff[x**2]
    B = coeff[x*y]
    C = coeff[y**2]
    D = coeff[x]
    E = coeff[y]
    F = coeff[1]

    A, B, C, D, E, F = map(as_int, _remove_gcd(A, B, C, D, E, F))

    # (1) Simple-Hyperbolic case: A = C = 0, B != 0
    # In this case equation can be converted to (Bx + E)(By + D) = DE - BF
    # We consider two cases; DE - BF = 0 and DE - BF != 0
    # More details, https://web.archive.org/web/20181231080858/https://www.alpertron.com.ar/METHODS.HTM#SHyperb

    sol = set()
    discr = B**2 - 4*A*C
    if A == 0 and C == 0 and B != 0:

        if D*E - B*F == 0:
            q, r = divmod(E, B)
            if not r:
                sol.add((-q, t))
            q, r = divmod(D, B)
            if not r:
                sol.add((t, -q))
        else:
            div = divisors(D*E - B*F)
            div = div + [-term for term in div]
            for d in div:
                x0, r = divmod(d - E, B)
                if not r:
                    q, r = divmod(D*E - B*F, d)
                    assert not r
                    y0, r = divmod(q - D, B)
                    if not r:
                        sol.add((x0, y0))

    # (2) Parabolic case: B**2 - 4*A*C = 0
    # There are two subcases to be considered in this case.
    # sqrt(c)D - sqrt(a)E = 0 and sqrt(c)D - sqrt(a)E != 0
    # More Details, https://web.archive.org/web/20181231080858/https://www.alpertron.com.ar/METHODS.HTM#Parabol

    elif discr == 0:

        if A == 0:
            s = _diop_quadratic([y, x], coeff, t)
            for soln in s:
                sol.add((soln[1], soln[0]))

        else:
            g = sign(A)*math.gcd(A, C)
            a = A // g
            c = C // g
            e = sign(B/A)
            sqa = math.isqrt(a)
            sqc = math.isqrt(c)
            _c = e*sqc*D - sqa*E
            if not _c:
                z = symbols('z', extended_real=True)
                roots = solve(sqa*g*z**2 + D*z + sqa*F, z)
                for root in roots:
                    root = root[z]
                    if isinstance(root, Integer):
                        ans = diop_solve(sqa*x + e*sqc*y - root)
                        sol.add((ans[0], ans[1]))

            else:
                def solve_x(u):
                    return -e*sqc*g*_c*t**2 - (E + 2*e*sqc*g*u)*t - (e*sqc*g*u**2 + E*u + e*sqc*F) // _c

                def solve_y(u):
                    return sqa*g*_c*t**2 + (D + 2*sqa*g*u)*t + (sqa*g*u**2 + D*u + sqa*F) // _c

                for z0 in range(0, abs(_c)):
                    if (divisible(sqa*g*z0**2 + D*z0 + sqa*F, _c) and
                            divisible(e*sqc**g*z0**2 + E*z0 + e*sqc*F, _c)):
                        sol.add((solve_x(z0), solve_y(z0)))

    # (3) Method used when B**2 - 4*A*C is a square, is described in p. 6 of the below paper
    # by John P. Robertson, https://web.archive.org/web/20180831180321/http://www.jpr2718.org/ax2p.pdf

    elif is_square(discr):
        if A != 0:
            r = sqrt(discr)
            u, v = symbols('u, v', integer=True)
            eq = _mexpand(4*A*r*u*v + 4*A*D*(B*v + r*u + r*v - B*u) +
                          2*A*4*A*E*(u - v) + 4*A*r*4*A*F)

            solution = diop_solve(eq, t)

            for s0, t0 in solution:

                num = B*t0 + r*s0 + r*t0 - B*s0
                x_0 = num/(4*A*r)
                y_0 = (s0 - t0)/(2*r)
                if ((isinstance(s0, Symbol) or isinstance(t0, Symbol)) and
                        check_param(x_0, y_0, 4*A*r, t) != (None, None)):
                    ans = check_param(x_0, y_0, 4*A*r, t)
                    sol.add((ans[0], ans[1]))
                elif x_0.is_Integer and y_0.is_Integer and is_solution_quad(var, coeff, x_0, y_0):
                    sol.add((x_0, y_0))

        else:
            s = _diop_quadratic(var[::-1], coeff, t)  # Interchange x and y
            while s:
                sol.add(s.pop()[::-1])  # and solution <--

    # (4) B**2 - 4*A*C > 0 and B**2 - 4*A*C not a square or B**2 - 4*A*C < 0

    else:
        P, Q = _transformation_to_DN(var, coeff)
        D, N = _find_DN(var, coeff)
        solns_pell = diop_DN(D, N)

        if D < 0:
            for x0, y0 in solns_pell:
                for x in [-x0, x0]:
                    for y in [-y0, y0]:
                        s = P*Matrix([x, y]) + Q
                        try:
                            sol.add(tuple(as_int(_) for _ in s))
                        except ValueError:
                            pass
        else:
            # In this case equation can be transformed into a Pell equation

            solns_pell = set(solns_pell)
            for X, Y in list(solns_pell):
                solns_pell.add((-X, -Y))

            a = diop_DN(D, 1)
            T = a[0][0]
            U = a[0][1]

            if all(_is_int(_) for _ in P[:4] + Q[:2]):
                for r, s in solns_pell:
                    _a = (r + s*sqrt(D))*(T + U*sqrt(D))**t
                    _b = (r - s*sqrt(D))*(T - U*sqrt(D))**t
                    x_n = _mexpand((_a + _b)/2)
                    y_n = _mexpand((_a - _b)/(2*sqrt(D)))
                    s = P*Matrix([x_n, y_n]) + Q
                    sol.add(tuple(s))

            else:
                L = math.lcm(*[_.denominator for _ in P[:4] + Q[:2]])

                k = 1

                T_k = T
                U_k = U

                while (T_k - 1) % L != 0 or U_k % L != 0:
                    T_k, U_k = T_k*T + D*U_k*U, T_k*U + U_k*T
                    k += 1

                for X, Y in solns_pell:

                    for _ in range(k):
                        if all(_is_int(_) for _ in P*Matrix([X, Y]) + Q):
                            _a = (X + sqrt(D)*Y)*(T_k + sqrt(D)*U_k)**t
                            _b = (X - sqrt(D)*Y)*(T_k - sqrt(D)*U_k)**t
                            Xt = (_a + _b)/2
                            Yt = (_a - _b)/(2*sqrt(D))
                            s = P*Matrix([Xt, Yt]) + Q
                            sol.add(tuple(s))

                        X, Y = X*T + D*U*Y, X*U + Y*T

    return sol


def is_solution_quad(var, coeff, u, v):
    """
    Check whether `(u, v)` is solution to the quadratic binary diophantine
    equation with the variable list ``var`` and coefficient dictionary
    ``coeff``.

    Not intended for use by normal users.

    """
    reps = dict(zip(var, (u, v)))
    eq = Add(*[j*i.xreplace(reps) for i, j in coeff.items()])
    return _mexpand(eq) == 0


def diop_DN(D, N, t=symbols('t', integer=True)):
    """
    Solves the equation `x^2 - Dy^2 = N`.

    Mainly concerned with the case `D > 0, D` is not a perfect square,
    which is the same as the generalized Pell equation. The LMM
    algorithm is used to solve this equation.

    Returns
    =======

    A tuple of pairs, (`x, y)`, for each class of the solutions.
    Other solutions of the class can be constructed according to the
    values of ``D`` and ``N``.

    Parameters
    ==========

    D, N : Integer
        correspond to D and N in the equation.
    t : Symbol, optional
        is the parameter to be used in the solutions.

    Examples
    ========

    >>> diop_DN(13, -4)  # Solves equation x**2 - 13*y**2 = -4
    [(3, 1), (393, 109), (36, 10)]

    The output can be interpreted as follows: There are three fundamental
    solutions to the equation `x^2 - 13y^2 = -4` given by (3, 1), (393, 109)
    and (36, 10). Each tuple is in the form (x, y), i.e. solution (3, 1) means
    that `x = 3` and `y = 1`.

    >>> diop_DN(986, 1)  # Solves equation x**2 - 986*y**2 = 1
    [(49299, 1570)]

    See Also
    ========

    diofant.solvers.diophantine.find_DN
    diofant.solvers.diophantine.diop_bf_DN

    References
    ==========

    * Solving the generalized Pell equation x**2 - D*y**2 = N, John P.
      Robertson, July 31, 2004, Pages 16 - 17. [online], Available:
      https://web.archive.org/web/20180831180333/http://www.jpr2718.org/pell.pdf

    """
    if D < 0:
        if N == 0:
            return [(0, 0)]
        elif N < 0:
            return []
        else:  # N > 0
            sol = []
            for d in divisors(square_factor(N)):
                sols = cornacchia(1, -D, N // d**2)
                if sols:
                    for x, y in sols:
                        sol.append((d*x, d*y))
                        if D == -1:
                            sol.append((d*y, d*x))
            return sol

    elif D == 0:
        if N < 0:
            return []
        if N == 0:
            return [(0, t)]
        sN, _exact = integer_nthroot(N, 2)
        if _exact:
            return [(sN, t)]
        else:
            return []

    else:  # D > 0
        sD, _exact = integer_nthroot(D, 2)
        if _exact:
            if N == 0:
                return [(sD*t, t)]
            else:
                sol = []

                for y in range(floor(sign(N)*(N - 1)/(2*sD)) + 1):
                    try:
                        sq, _exact = integer_nthroot(D*y**2 + N, 2)
                    except ValueError:
                        _exact = False
                    if _exact:
                        sol.append((sq, y))

                return sol
        elif 1 < N**2 < D:
            return _special_diop_DN(D, N)
        else:
            if N == 0:
                return [(0, 0)]

            elif abs(N) == 1:

                pqa = PQa(0, 1, D)
                j = 0
                G = []
                B = []

                for i in pqa:  # pragma: no branch

                    a = i[2]
                    G.append(i[5])
                    B.append(i[4])

                    if j != 0 and a == 2*sD:
                        break
                    j += 1

                if _odd(j):

                    if N == -1:
                        x = G[j - 1]
                        y = B[j - 1]
                    else:
                        count = j
                        while count < 2*j - 1:
                            i = next(pqa)
                            G.append(i[5])
                            B.append(i[4])
                            count += 1

                        x = G[count]
                        y = B[count]
                else:
                    if N == 1:
                        x = G[j - 1]
                        y = B[j - 1]
                    else:
                        return []

                return [(x, y)]

            else:

                fs = []
                sol = []
                div = divisors(N)

                for d in div:
                    if divisible(N, d**2):
                        fs.append(d)

                for f in fs:
                    m = N // f**2

                    zs = sqrt_mod(D, abs(m), all_roots=True)
                    zs = [i for i in zs if i <= abs(m) // 2]

                    if abs(m) != 2:
                        zs = zs + [-i for i in zs if i]  # Omit duplicate zero

                    for z in zs:

                        pqa = PQa(z, abs(m), D)
                        j = 0
                        G = []
                        B = []

                        for i in pqa:  # pragma: no branch

                            G.append(i[5])
                            B.append(i[4])

                            if j != 0 and abs(i[1]) == 1:
                                r = G[j - 1]
                                s = B[j - 1]

                                if r**2 - D*s**2 == m:
                                    sol.append((f*r, f*s))

                                elif diop_DN(D, -1):
                                    a = diop_DN(D, -1)
                                    sol.append((f*(r*a[0][0] + a[0][1]*s*D), f*(r*a[0][1] + s*a[0][0])))

                                break

                            j += 1
                            if j == length(z, abs(m), D):
                                break

                return sol


def _special_diop_DN(D, N):
    """
    Solves the equation `x^2 - Dy^2 = N` for the special case where
    `1 < N**2 < D` and `D` is not a perfect square.

    References
    ==========

    * :cite:`Andreescu15`, Section 4.4.4.

    """
    assert 1 < N**2 < D
    assert not integer_nthroot(D, 2)[1]

    sqrt_D = sqrt(D)
    F = [(N, 1)]
    f = 2
    while True:
        f2 = f**2
        if f2 > abs(N):
            break
        n, r = divmod(N, f2)
        if r == 0:
            F.append((n, f))
        f += 1

    P = 0
    Q = 1
    G0, G1 = 0, 1
    B0, B1 = 1, 0

    solutions = []

    i = 0
    while True:
        a = floor((P + sqrt_D) / Q)
        P = a*Q - P
        Q = (D - P**2) // Q
        G2 = a*G1 + G0
        B2 = a*B1 + B0

        for n, f in F:
            if G2**2 - D*B2**2 == n:
                solutions.append((f*G2, f*B2))

        i += 1
        if Q == 1 and i % 2 == 0:
            break

        G0, G1 = G1, G2
        B0, B1 = B1, B2

    return solutions


def cornacchia(a, b, m):
    r"""
    Solves `ax^2 + by^2 = m` where `\gcd(a, b) = 1 = gcd(a, m)` and `a, b > 0`.

    Uses the algorithm due to Cornacchia. The method only finds primitive
    solutions, i.e. ones with `\gcd(x, y) = 1`. So this method can't be used to
    find the solutions of `x^2 + y^2 = 20` since the only solution to former is
    `(x, y) = (4, 2)` and it is not primitive. When `a = b`, only the
    solutions with `x \leq y` are found. For more details, see the References.

    Examples
    ========

    >>> cornacchia(2, 3, 35)
    {(2, 3), (4, 1)}
    >>> cornacchia(1, 1, 25)
    {(4, 3)}

    References
    ===========

    * A. Nitaj, "L'algorithme de Cornacchia"
    * Solving the diophantine equation ax**2 + by**2 = m by Cornacchia's
      method, [online], Available:
      http://www.numbertheory.org/php/cornacchia.html

    See Also
    ========

    diofant.utilities.iterables.signed_permutations

    """
    sols = set()

    a1 = igcdex(a, m)[0]
    v = sqrt_mod(-b*a1, m, all_roots=True)
    if not v:
        return

    for t in v:
        if t < m // 2:
            continue

        u, r = t, m

        while True:
            u, r = r, u % r
            if a*r**2 < m:
                break

        m1 = m - a*r**2

        if m1 % b == 0:
            m1 = m1 // b
            s, _exact = integer_nthroot(m1, 2)
            if _exact:
                if a == b and r < s:
                    r, s = s, r
                sols.add((int(r), int(s)))

    return sols


def PQa(P_0, Q_0, D):
    r"""
    Returns useful information needed to solve the Pell equation.

    There are six sequences of integers defined related to the continued
    fraction representation of `\frac{P + \sqrt{D}}{Q}`, namely {`P_{i}`},
    {`Q_{i}`}, {`a_{i}`},{`A_{i}`}, {`B_{i}`}, {`G_{i}`}. ``PQa()`` Returns
    these values as a 6-tuple in the same order as mentioned above.

    Parameters
    ==========

    P_0, Q_0, D : Integer
        integers corresponding to `P_{0}`, `Q_{0}` and `D` in the continued
        fraction `\frac{P_{0} + \sqrt{D}}{Q_{0}}`.  Also it's assumed
        that `P_{0}^2 == D mod(|Q_{0}|)` and `D` is square free.

    Examples
    ========

    >>> pqa = PQa(13, 4, 5)  # (13 + sqrt(5))/4
    >>> next(pqa)  # (P_0, Q_0, a_0, A_0, B_0, G_0)
    (13, 4, 3, 3, 1, -1)
    >>> next(pqa)  # (P_1, Q_1, a_1, A_1, B_1, G_1)
    (-1, 1, 1, 4, 1, 3)

    References
    ==========

    * Solving the generalized Pell equation x^2 - Dy^2 = N, John P.
      Robertson, July 31, 2004, Pages 4 - 8. https://web.archive.org/web/20180831180333/http://www.jpr2718.org/pell.pdf

    """
    A_i_2 = B_i_1 = 0
    A_i_1 = B_i_2 = 1

    G_i_2 = -P_0
    G_i_1 = Q_0

    P_i = P_0
    Q_i = Q_0

    while 1:
        a_i = floor((P_i + sqrt(D))/Q_i)
        A_i = a_i*A_i_1 + A_i_2
        B_i = a_i*B_i_1 + B_i_2
        G_i = a_i*G_i_1 + G_i_2

        yield P_i, Q_i, a_i, A_i, B_i, G_i

        A_i_1, A_i_2 = A_i, A_i_1
        B_i_1, B_i_2 = B_i, B_i_1
        G_i_1, G_i_2 = G_i, G_i_1

        P_i = a_i*Q_i - P_i
        Q_i = (D - P_i**2)/Q_i


def diop_bf_DN(D, N, t=symbols('t', integer=True)):
    r"""
    Uses brute force to solve the equation, `x^2 - Dy^2 = N`.

    Mainly concerned with the generalized Pell equation which is the case when
    `D > 0, D` is not a perfect square.
    Let `(t, u)` be the minimal positive solution of the equation
    `x^2 - Dy^2 = 1`. Then this method requires
    `\sqrt{\frac{\mid N \mid (t \pm 1)}{2D}}` to be small.

    Parameters
    ==========

    D, N : Integer
        correspond to D and N in the equation.
    t : Symbol, optional
        is the parameter to be used in the solutions.

    Examples
    ========

    >>> diop_bf_DN(13, -4)
    [(3, 1), (-3, 1), (36, 10)]
    >>> diop_bf_DN(986, 1)
    [(49299, 1570)]

    See Also
    ========

    diofant.solvers.diophantine.diop_DN

    References
    ==========

    * Solving the generalized Pell equation x**2 - D*y**2 = N, John P.
      Robertson, July 31, 2004, Page 15. https://web.archive.org/web/20180831180333/http://www.jpr2718.org/pell.pdf

    """
    D = as_int(D)
    N = as_int(N)

    sol = []
    a = diop_DN(D, 1)
    u = a[0][0]

    if abs(N) == 1:
        return diop_DN(D, N)

    elif N > 1:
        L1 = 0
        L2 = integer_nthroot(int(N*(u - 1)/(2*D)), 2)[0] + 1

    elif N < -1:
        L1, _exact = integer_nthroot(-int(N/D), 2)
        if not _exact:
            L1 += 1
        L2 = integer_nthroot(-int(N*(u + 1)/(2*D)), 2)[0] + 1

    else:  # N = 0
        if D < 0:
            return [(0, 0)]
        elif D == 0:
            return [(0, t)]
        else:
            sD, _exact = integer_nthroot(D, 2)
            if _exact:
                return [(sD*t, t), (-sD*t, t)]
            else:
                return [(0, 0)]

    for y in range(L1, L2):
        try:
            x, _exact = integer_nthroot(N + D*y**2, 2)
        except ValueError:
            _exact = False
        if _exact:
            sol.append((x, y))
            if not equivalent(x, y, -x, y, D, N):
                sol.append((-x, y))

    return sol


def equivalent(u, v, r, s, D, N):
    """
    Returns True if two solutions `(u, v)` and `(r, s)` of `x^2 - Dy^2 = N`
    belongs to the same equivalence class and False otherwise.

    Two solutions `(u, v)` and `(r, s)` to the above equation fall to the same
    equivalence class iff both `(ur - Dvs)` and `(us - vr)` are divisible by
    `N`. No test is performed to test whether `(u, v)` and
    `(r, s)` are actually solutions to the equation. User should take care of
    this.

    Parameters
    ==========

    u, v, r, s, D, N : Integer

    Examples
    ========

    >>> equivalent(18, 5, -18, -5, 13, -1)
    True
    >>> equivalent(3, 1, -18, 393, 109, -4)
    False

    References
    ==========

    * Solving the generalized Pell equation x**2 - D*y**2 = N, John P.
      Robertson, July 31, 2004, Page 12. https://web.archive.org/web/20180831180333/http://www.jpr2718.org/pell.pdf

    """
    return divisible(u*r - D*v*s, N) and divisible(u*s - v*r, N)


def length(P, Q, D):
    r"""
    Returns the (length of aperiodic part + length of periodic part) of
    continued fraction representation of `\frac{P + \sqrt{D}}{Q}`.

    It is important to remember that this does NOT return the length of the
    periodic part but the sum of the lengths of the two parts as mentioned
    above.

    Parameters
    ==========

    P, D, Q : Integer

    Examples
    ========

    >>> length(-2, 4, 5)  # (-2 + sqrt(5))/4
    3
    >>> length(-5, 4, 17)  # (-5 + sqrt(17))/4
    5

    See Also
    ========

    diofant.ntheory.continued_fraction.continued_fraction_periodic

    """
    from ..ntheory import continued_fraction_periodic
    v = continued_fraction_periodic(P, Q, D)
    if type(v[-1]) is list:
        rpt = len(v[-1])
        nonrpt = len(v) - 1
    else:
        rpt = 0
        nonrpt = len(v)
    return rpt + nonrpt


def transformation_to_DN(eq):
    """
    This function transforms general quadratic,
    `ax^2 + bxy + cy^2 + dx + ey + f = 0`
    to more easy to deal with `X^2 - DY^2 = N` form.

    This is used to solve the general quadratic equation by transforming it to
    the latter form.  This function returns a tuple (A, B) where A is a 2 X 2
    matrix and B is a 2 X 1 matrix such that,

    Transpose([x y]) =  A * Transpose([X Y]) + B

    Parameters
    ==========

    eq : Expr
        the quadratic expression to be transformed.

    Examples
    ========

    >>> A, B = transformation_to_DN(x**2 - 3*x*y - y**2 - 2*y + 1)
    >>> A
    Matrix([
    [1/26, 3/26],
    [   0, 1/13]])
    >>> B
    Matrix([
    [-6/13],
    [-4/13]])

    A, B  returned are such that Transpose((x y)) =  A * Transpose((X Y)) + B.
    Substituting these values for `x` and `y` and a bit of simplifying work
    will give an equation of the form `x^2 - Dy^2 = N`.

    >>> from diofant.abc import X, Y
    >>> u = (A*Matrix([X, Y]) + B)[0]  # Transformation for x
    >>> u
    X/26 + 3*Y/26 - 6/13
    >>> v = (A*Matrix([X, Y]) + B)[1]  # Transformation for y
    >>> v
    Y/13 - 4/13

    Next we will substitute these formulas for `x` and `y` and do
    ``simplify()``.

    >>> eq = simplify((x**2 - 3*x*y - y**2 - 2*y + 1).subs({x: u, y: v}))
    >>> eq
    X**2/676 - Y**2/52 + 17/13

    By multiplying the denominator appropriately, we can get a Pell equation
    in the standard form.

    >>> eq * 676
    X**2 - 13*Y**2 + 884

    If only the final equation is needed, ``find_DN()`` can be used.

    See Also
    ========

    diofant.solvers.diophantine.find_DN

    References
    ==========

    * Solving the equation ax^2 + bxy + cy^2 + dx + ey + f = 0,
      John P.Robertson, May 8, 2003, Page 7 - 11.
      https://web.archive.org/web/20180831180321/http://www.jpr2718.org/ax2p.pdf

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)
    if diop_type == 'binary_quadratic':
        return _transformation_to_DN(var, coeff)


def _transformation_to_DN(var, coeff):

    x, y = var

    a = coeff[x**2]
    b = coeff[x*y]
    c = coeff[y**2]
    d = coeff[x]
    e = coeff[y]
    f = coeff[1]

    a, b, c, d, e, f = map(as_int, _remove_gcd(a, b, c, d, e, f))

    X, Y = symbols('X, Y', integer=True)

    if b:
        B, C = _rational_pq(2*a, b)
        A, T = _rational_pq(a, B**2)

        # eq_1 = A*B*X**2 + B*(c*T - A*C**2)*Y**2 + d*T*X + (B*e*T - d*T*C)*Y + f*T*B
        coeff = {X**2: A*B, X*Y: 0, Y**2: B*(c*T - A*C**2), X: d*T, Y: B*e*T - d*T*C, 1: f*T*B}
        A_0, B_0 = _transformation_to_DN([X, Y], coeff)
        return Matrix(2, 2, [Integer(1)/B, -C/B, 0, 1])*A_0, Matrix(2, 2, [Integer(1)/B, -C/B, 0, 1])*B_0

    else:
        if d:
            B, C = _rational_pq(2*a, d)
            A, T = _rational_pq(a, B**2)

            # eq_2 = A*X**2 + c*T*Y**2 + e*T*Y + f*T - A*C**2
            coeff = {X**2: A, X*Y: 0, Y**2: c*T, X: 0, Y: e*T, 1: f*T - A*C**2}
            A_0, B_0 = _transformation_to_DN([X, Y], coeff)
            return Matrix(2, 2, [Integer(1)/B, 0, 0, 1])*A_0, Matrix(2, 2, [Integer(1)/B, 0, 0, 1])*B_0 + Matrix([-C/B, 0])

        else:
            if e:
                B, C = _rational_pq(2*c, e)
                A, T = _rational_pq(c, B**2)

                # eq_3 = a*T*X**2 + A*Y**2 + f*T - A*C**2
                coeff = {X**2: a*T, X*Y: 0, Y**2: A, X: 0, Y: 0, 1: f*T - A*C**2}
                A_0, B_0 = _transformation_to_DN([X, Y], coeff)
                return Matrix(2, 2, [1, 0, 0, Integer(1)/B])*A_0, Matrix(2, 2, [1, 0, 0, Integer(1)/B])*B_0 + Matrix([0, -C/B])

            else:
                # TODO: pre-simplification: Not necessary but may simplify
                # the equation.

                return Matrix(2, 2, [Integer(1)/a, 0, 0, 1]), Matrix([0, 0])


def find_DN(eq):
    """
    This function returns a tuple, `(D, N)` of the simplified form,
    `x^2 - Dy^2 = N`, corresponding to the general quadratic,
    `ax^2 + bxy + cy^2 + dx + ey + f = 0`.

    Solving the general quadratic is then equivalent to solving the equation
    `X^2 - DY^2 = N` and transforming the solutions by using the transformation
    matrices returned by ``transformation_to_DN()``.

    Parameters
    ==========

    eq : Expr
        is the quadratic expression to be transformed.

    Examples
    ========

    >>> find_DN(x**2 - 3*x*y - y**2 - 2*y + 1)
    (13, -884)

    Interpretation of the output is that we get `X^2 -13Y^2 = -884` after
    transforming `x^2 - 3xy - y^2 - 2y + 1` using the transformation returned
    by ``transformation_to_DN()``.

    See Also
    ========

    diofant.solvers.diophantine.transformation_to_DN

    References
    ==========

    * Solving the equation ax^2 + bxy + cy^2 + dx + ey + f = 0,
      John P.Robertson, May 8, 2003, Page 7 - 11.
      https://web.archive.org/web/20180831180321/http://www.jpr2718.org/ax2p.pdf

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)
    if diop_type == 'binary_quadratic':
        return _find_DN(var, coeff)


def _find_DN(var, coeff):

    x, y = var
    X, Y = symbols('X, Y', integer=True)
    A, B = _transformation_to_DN(var, coeff)

    u = (A*Matrix([X, Y]) + B)[0]
    v = (A*Matrix([X, Y]) + B)[1]
    eq = x**2*coeff[x**2] + x*y*coeff[x*y] + y**2*coeff[y**2] + x*coeff[x] + y*coeff[y] + coeff[1]
    simplified = _mexpand(eq.subs({x: u, y: v}))

    coeff = simplified.as_coefficients_dict()

    return -coeff[Y**2]/coeff[X**2], -coeff[1]/coeff[X**2]


def check_param(x, y, a, t):
    """
    If there is a number modulo ``a`` such that ``x`` and ``y`` are both
    integers, then return a parametric representation for ``x`` and ``y``
    else return (None, None).

    Here ``x`` and ``y`` are functions of ``t``.

    """
    from ..simplify.simplify import clear_coefficients

    if x.is_number and not x.is_Integer:
        return None, None

    if y.is_number and not y.is_Integer:
        return None, None

    m, n = symbols('m, n', integer=True)
    c, _ = (m*x + n*y).as_content_primitive()
    if a % c.denominator:
        return None, None

    # clear_coefficients(mx + b, R)[1] -> (R - b)/m
    eq = clear_coefficients(x, m)[1] - clear_coefficients(y, n)[1]
    _, eq = eq.as_content_primitive()

    return diop_solve(eq, t)


def diop_ternary_quadratic(eq):
    """
    Solves the general quadratic ternary form,
    `ax^2 + by^2 + cz^2 + fxy + gyz + hxz = 0`.

    Returns
    =======

    tuple
       which is a base solution for the above equation. If there are no
       solutions, ``(None, None, None)`` is returned.

    Parameters
    ===========

    eq : Expr
        should be an homogeneous expression of degree two in three variables
        and it is assumed to be zero.

    Examples
    ========

    >>> diop_ternary_quadratic(x**2 + 3*y**2 - z**2)
    (1, 0, 1)
    >>> diop_ternary_quadratic(4*x**2 + 5*y**2 - z**2)
    (1, 0, 2)
    >>> diop_ternary_quadratic(45*x**2 - 7*y**2 - 8*x*y - z**2)
    (28, 45, 105)
    >>> diop_ternary_quadratic(x**2 - 49*y**2 - z**2 + 13*z*y - 8*x*y)
    (9, 1, 5)

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type in ('homogeneous_ternary_quadratic',
                     'homogeneous_ternary_quadratic_normal'):
        return _diop_ternary_quadratic(var, coeff)


def _diop_ternary_quadratic(_var, coeff):

    x, y, z = _var
    var = [x, y, z]

    # Equations of the form B*x*y + C*z*x + E*y*z = 0 and At least two of the
    # coefficients A, B, C are non-zero.
    # There are infinitely many solutions for the equation.
    # Ex: (0, 0, t), (0, t, 0), (t, 0, 0)
    # Equation can be re-written as y*(B*x + E*z) = -C*x*z and we can find rather
    # unobviuos solutions. Set y = -C and B*x + E*z = x*z. The latter can be solved by
    # using methods for binary quadratic diophantine equations. Let's select the
    # solution which minimizes |x| + |z|

    if not any(coeff[i**2] for i in var):
        if coeff[x*z]:
            sols = diophantine(coeff[x*y]*x + coeff[y*z]*z - x*z)
            s = sols.pop()
            min_sum = abs(s[0]) + abs(s[1])

            for r in sols:
                if abs(r[0]) + abs(r[1]) < min_sum:
                    s = r
                    min_sum = abs(s[0]) + abs(s[1])

                x_0, y_0, z_0 = s[0], -coeff[x*z], s[1]

        else:
            var[0], var[1] = _var[1], _var[0]
            y_0, x_0, z_0 = _diop_ternary_quadratic(var, coeff)

        return _remove_gcd(x_0, y_0, z_0)

    if coeff[x**2] == 0:
        # If the coefficient of x is zero change the variables
        if coeff[y**2] == 0:
            var[0], var[2] = _var[2], _var[0]
            z_0, y_0, x_0 = _diop_ternary_quadratic(var, coeff)

        else:
            var[0], var[1] = _var[1], _var[0]
            y_0, x_0, z_0 = _diop_ternary_quadratic(var, coeff)

    else:
        if coeff[x*y] or coeff[x*z]:
            # Apply the transformation x --> X - (B*y + C*z)/(2*A)
            A = coeff[x**2]
            B = coeff[x*y]
            C = coeff[x*z]
            D = coeff[y**2]
            E = coeff[y*z]
            F = coeff[z**2]

            _coeff = {}

            _coeff[x**2] = 4*A**2
            _coeff[y**2] = 4*A*D - B**2
            _coeff[z**2] = 4*A*F - C**2
            _coeff[y*z] = 4*A*E - 2*B*C
            _coeff[x*y] = 0
            _coeff[x*z] = 0

            x_0, y_0, z_0 = _diop_ternary_quadratic(var, _coeff)

            if x_0 is None:
                return None, None, None

            p, q = _rational_pq(B*y_0 + C*z_0, 2*A)
            x_0, y_0, z_0 = x_0*q - p, y_0*q, z_0*q

        elif coeff[z*y] != 0:
            if coeff[y**2] == 0:
                if coeff[z**2] == 0:
                    # Equations of the form A*x**2 + E*yz = 0.
                    A = coeff[x**2]
                    E = coeff[y*z]

                    b, a = _rational_pq(-E, A)

                    x_0, y_0, z_0 = b, a, b

                else:
                    # Ax**2 + E*y*z + F*z**2  = 0
                    var[0], var[2] = _var[2], _var[0]
                    z_0, y_0, x_0 = _diop_ternary_quadratic(var, coeff)

            else:
                # A*x**2 + D*y**2 + E*y*z + F*z**2 = 0, C may be zero
                var[0], var[1] = _var[1], _var[0]
                y_0, x_0, z_0 = _diop_ternary_quadratic(var, coeff)

        else:
            # Ax**2 + D*y**2 + F*z**2 = 0, C may be zero
            x_0, y_0, z_0 = _diop_ternary_quadratic_normal(var, coeff)

    return _remove_gcd(x_0, y_0, z_0)


def transformation_to_normal(eq):
    """
    Returns the transformation Matrix that converts a general ternary
    quadratic equation `eq` (`ax^2 + by^2 + cz^2 + dxy + eyz + fxz`)
    to a form without cross terms: `ax^2 + by^2 + cz^2 = 0`. This is
    not used in solving ternary quadratics; it is only implemented for
    the sake of completeness.

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type in ('homogeneous_ternary_quadratic',
                     'homogeneous_ternary_quadratic_normal'):
        return _transformation_to_normal(var, coeff)


def _transformation_to_normal(var, coeff):

    _var = list(var)  # copy
    x, y, z = var

    if not any(coeff[i**2] for i in var):
        # https://math.stackexchange.com/questions/448051/transform-quadratic-ternary-form-to-normal-form/448065#448065
        a = coeff[x*y]
        b = coeff[y*z]
        c = coeff[x*z]
        swap = False
        if not a:  # b can't be 0 or else there aren't 3 vars
            swap = True
            a, b = b, a
        T = Matrix(((1, 1, -b/a), (1, -1, -c/a), (0, 0, 1)))
        if swap:
            T.row_swap(0, 1)
            T.col_swap(0, 1)
        return T

    if coeff[x**2] == 0:
        # If the coefficient of x is zero change the variables
        if coeff[y**2] == 0:
            _var[0], _var[2] = var[2], var[0]
            T = _transformation_to_normal(_var, coeff)
            T.row_swap(0, 2)
            T.col_swap(0, 2)
            return T

        else:
            _var[0], _var[1] = var[1], var[0]
            T = _transformation_to_normal(_var, coeff)
            T.row_swap(0, 1)
            T.col_swap(0, 1)
            return T

    # Apply the transformation x --> X - (B*Y + C*Z)/(2*A)
    if coeff[x*y] != 0 or coeff[x*z] != 0:
        A = coeff[x**2]
        B = coeff[x*y]
        C = coeff[x*z]
        D = coeff[y**2]
        E = coeff[y*z]
        F = coeff[z**2]

        _coeff = {x**2: 4*A**2, y**2: 4*A*D - B**2, z**2: 4*A*F - C**2,
                  y*z: 4*A*E - 2*B*C, x*y: 0, x*z: 0}

        T_0 = _transformation_to_normal(_var, _coeff)
        return Matrix(3, 3, [1, Integer(-B)/(2*A), Integer(-C)/(2*A), 0, 1, 0, 0, 0, 1])*T_0

    elif coeff[y*z] != 0:
        if coeff[y**2] == 0:
            if coeff[z**2] == 0:
                # Equations of the form A*x**2 + E*yz = 0.
                # Apply transformation y -> Y + Z ans z -> Y - Z
                return Matrix(3, 3, [1, 0, 0, 0, 1, 1, 0, 1, -1])

            else:
                # Ax**2 + E*y*z + F*z**2  = 0
                _var[0], _var[2] = var[2], var[0]
                T = _transformation_to_normal(_var, coeff)
                T.row_swap(0, 2)
                T.col_swap(0, 2)
                return T

        else:
            # A*x**2 + D*y**2 + E*y*z + F*z**2 = 0, F may be zero
            _var[0], _var[1] = var[1], var[0]
            T = _transformation_to_normal(_var, coeff)
            T.row_swap(0, 1)
            T.col_swap(0, 1)
            return T

    else:
        return Matrix.eye(3)


def parametrize_ternary_quadratic(eq):
    """
    Returns the parametrized general solution for the ternary quadratic
    equation ``eq`` which has the form
    `ax^2 + by^2 + cz^2 + fxy + gyz + hxz = 0`.

    Examples
    ========

    >>> parametrize_ternary_quadratic(x**2 + y**2 - z**2)
    (2*p*q, p**2 - q**2, p**2 + q**2)

    Here `p` and `q` are two co-prime integers.

    >>> parametrize_ternary_quadratic(3*x**2 + 2*y**2 - z**2 - 2*x*y + 5*y*z - 7*y*z)
    (2*p**2 - 2*p*q - q**2, 2*p**2 + 2*p*q - q**2, 2*p**2 - 2*p*q + 3*q**2)
    >>> parametrize_ternary_quadratic(124*x**2 - 30*y**2 - 7729*z**2)
    (-1410*p**2 - 363263*q**2, 2700*p**2 + 30916*p*q - 695610*q**2, -60*p**2 + 5400*p*q + 15458*q**2)

    References
    ==========

    * The algorithmic resolution of Diophantine equations, Nigel P. Smart,
      London Mathematical Society Student Texts 41, Cambridge University
      Press, Cambridge, 1998.

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type in ('homogeneous_ternary_quadratic',
                     'homogeneous_ternary_quadratic_normal'):
        x_0, y_0, z_0 = _diop_ternary_quadratic(var, coeff)
        return _parametrize_ternary_quadratic((x_0, y_0, z_0), var, coeff)


def _parametrize_ternary_quadratic(solution, _var, coeff):
    # called for a*x**2 + b*y**2 + c*z**2 + d*x*y + e*y*z + f*x*z = 0
    assert 1 not in coeff

    x_0, y_0, z_0 = solution

    v = list(_var)  # copy

    if x_0 is None:
        return None, None, None

    if solution.count(0) >= 2:
        # if there are 2 zeros the the equation reduces
        # to k*X**2 == 0 where X is x, y, or z so X must
        # be zero, too. So there is only the trivial solution.
        return None, None, None

    if x_0 == 0:
        v[0], v[1] = v[1], v[0]
        y_p, x_p, z_p = _parametrize_ternary_quadratic((y_0, x_0, z_0), v, coeff)
        return x_p, y_p, z_p

    x, y, z = v
    r, p, q = symbols('r, p, q', integer=True)

    eq = sum(k*v for k, v in coeff.items())
    eq_1 = _mexpand(eq.subs({x: r*x_0, y: r*y_0 + p, z: r*z_0 + q}))
    A, B = eq_1.as_independent(r, as_Add=True)

    x = A*x_0
    y = (A*y_0 - _mexpand(B/r*p))
    z = (A*z_0 - _mexpand(B/r*q))

    return x, y, z


def diop_ternary_quadratic_normal(eq):
    """
    Solves the quadratic ternary diophantine equation,
    `ax^2 + by^2 + cz^2 = 0`.

    Here the coefficients `a`, `b`, and `c` should be non zero. Otherwise the
    equation will be a quadratic binary or univariate equation. If solvable,
    returns a tuple `(x, y, z)` that satisifes the given equation. If the
    equation does not have integer solutions, `(None, None, None)` is returned.

    Examples
    ========

    >>> diop_ternary_quadratic_normal(x**2 + 3*y**2 - z**2)
    (1, 0, 1)
    >>> diop_ternary_quadratic_normal(4*x**2 + 5*y**2 - z**2)
    (1, 0, 2)
    >>> diop_ternary_quadratic_normal(34*x**2 - 3*y**2 - 301*z**2)
    (4, 9, 1)

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type == 'homogeneous_ternary_quadratic_normal':
        return _diop_ternary_quadratic_normal(var, coeff)


def _diop_ternary_quadratic_normal(var, coeff):

    x, y, z = var

    a = coeff[x**2]
    b = coeff[y**2]
    c = coeff[z**2]
    try:
        assert len([k for k in coeff if coeff[k]]) == 3
        assert all(coeff[i**2] for i in var)
    except AssertionError as exc:
        raise ValueError(filldedent("""
    coeff dict is not consistent with assumption of this routine:
    coefficients should be those of an expression in the form
    a*x**2 + b*y**2 + c*z**2 where a*b*c != 0.""")) from exc

    (sqf_of_a, sqf_of_b, sqf_of_c), (a_1, b_1, c_1), (a_2, b_2, c_2) = sqf_normal(a, b, c, steps=True)

    A = -a_2*c_2
    B = -b_2*c_2

    # If following two conditions are satisified then there are no solutions
    if A < 0 and B < 0:
        return None, None, None

    if any(_ is None for _ in [sqrt_mod(-b_2*c_2, a_2),
                               sqrt_mod(-c_2*a_2, b_2),
                               sqrt_mod(-a_2*b_2, c_2)]):
        return None, None, None

    z_0, x_0, y_0 = descent(A, B)

    z_0, q = _rational_pq(z_0, abs(c_2))
    x_0 *= q
    y_0 *= q

    x_0, y_0, z_0 = _remove_gcd(x_0, y_0, z_0)

    # Holzer reduction
    if sign(a) == sign(b):
        x_0, y_0, z_0 = holzer(x_0, y_0, z_0, abs(a_2), abs(b_2), abs(c_2))
    elif sign(a) == sign(c):
        x_0, z_0, y_0 = holzer(x_0, z_0, y_0, abs(a_2), abs(c_2), abs(b_2))
    else:
        y_0, z_0, x_0 = holzer(y_0, z_0, x_0, abs(b_2), abs(c_2), abs(a_2))

    x_0 = reconstruct(b_1, c_1, x_0)
    y_0 = reconstruct(a_1, c_1, y_0)
    z_0 = reconstruct(a_1, b_1, z_0)

    sq_lcm = math.lcm(sqf_of_a, sqf_of_b, sqf_of_c)

    x_0 = abs(x_0*sq_lcm//sqf_of_a)
    y_0 = abs(y_0*sq_lcm//sqf_of_b)
    z_0 = abs(z_0*sq_lcm//sqf_of_c)

    return _remove_gcd(x_0, y_0, z_0)


def sqf_normal(a, b, c, steps=False):
    """
    Return `a', b', c'`, the coefficients of the square-free normal
    form of `ax^2 + by^2 + cz^2 = 0`, where `a', b', c'` are pairwise
    prime.  If `steps` is True then also return three tuples:
    `sq`, `sqf`, and `(a', b', c')` where `sq` contains the square
    factors of `a`, `b` and `c` after removing the `gcd(a, b, c)`;
    `sqf` contains the values of `a`, `b` and `c` after removing
    both the `gcd(a, b, c)` and the square factors.

    The solutions for `ax^2 + by^2 + cz^2 = 0` can be
    recovered from the solutions of `a'x^2 + b'y^2 + c'z^2 = 0`.

    Examples
    ========

    >>> sqf_normal(2 * 3**2 * 5, 2 * 5 * 11, 2 * 7**2 * 11)
    (11, 1, 5)
    >>> sqf_normal(2 * 3**2 * 5, 2 * 5 * 11, 2 * 7**2 * 11, True)
    ((3, 1, 7), (5, 55, 11), (11, 1, 5))

    References
    ==========

    * Legendre's Theorem, Legrange's Descent,
      https://public.csusm.edu/aitken_html/notes/legendre.pdf

    See Also
    ========

    diofant.solvers.diophantine.reconstruct

    """
    ABC = _remove_gcd(a, b, c)
    sq = tuple(square_factor(i) for i in ABC)
    sqf = A, B, C = tuple(i//j**2 for i, j in zip(ABC, sq))
    pc = math.gcd(A, B)
    A //= pc
    B //= pc
    pa = math.gcd(B, C)
    B //= pa
    C //= pa
    pb = math.gcd(A, C)
    A //= pb
    B //= pb

    A *= pa
    B *= pb
    C *= pc

    if steps:
        return sq, sqf, (A, B, C)
    else:
        return A, B, C


def reconstruct(A, B, z):
    """
    Reconstruct the `z` value of an equivalent solution of `ax^2 + by^2 + cz^2`
    from the `z` value of a solution of the square-free normal form of the
    equation, `a'*x^2 + b'*y^2 + c'*z^2`, where `a'`, `b'` and `c'` are square
    free and `gcd(a', b', c') == 1`.

    """
    f = factorint(math.gcd(A, B))
    for p, e in f.items():
        if e != 1:
            raise ValueError('a and b should be square-free')
        z *= p
    return z


def ldescent(A, B):
    r"""
    Return a non-trivial solution to `w^2 = Ax^2 + By^2` using
    Lagrange's method; return None if there is no such solution.

    Here, `A \neq 0` and `B \neq 0` and `A` and `B` are square free. Output a
    tuple `(w_0, x_0, y_0)` which is a solution to the above equation.

    Examples
    ========

    >>> ldescent(1, 1)  # w^2 = x^2 + y^2
    (1, 1, 0)
    >>> ldescent(4, -7)  # w^2 = 4x^2 - 7y^2
    (2, -1, 0)

    This means that `x = -1, y = 0` and `w = 2` is a solution to the equation
    `w^2 = 4x^2 - 7y^2`

    >>> ldescent(5, -1)  # w^2 = 5x^2 - y^2
    (2, 1, -1)

    References
    ==========

    * The algorithmic resolution of Diophantine equations, Nigel P. Smart,
      London Mathematical Society Student Texts 41, Cambridge University
      Press, Cambridge, 1998.
    * Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
      Mathematics of Computation, Volume 00, Number 0,
      http://eprints.nottingham.ac.uk/60/1/kvxefz87.pdf

    """
    if abs(A) > abs(B):
        w, y, x = ldescent(B, A)
        return w, x, y

    if A == 1:
        return 1, 1, 0

    if B == 1:
        return 1, 0, 1

    if B == -1:  # and A == -1
        return

    r = sqrt_mod(A, B)

    Q = (r**2 - A) // B

    if Q == 0:
        B_0 = 1
        d = 0
    else:
        div = divisors(Q)
        B_0 = None

        for i in div:  # pragma: no cover
            sQ, _exact = integer_nthroot(abs(Q) // i, 2)
            if _exact:
                B_0, d = sign(Q)*i, sQ
                break

    assert B_0 is not None

    W, X, Y = ldescent(A, B_0)
    return _remove_gcd((-A*X + r*W), (r*X - W), Y*(B_0*d))


def descent(A, B):
    """
    Returns a non-trivial solution, (x, y, z), to `x^2 = Ay^2 + Bz^2`
    using Lagrange's descent method with lattice-reduction. `A` and `B`
    are assumed to be valid for such a solution to exist.

    This is faster than the normal Lagrange's descent algorithm because
    the Gaussian reduction is used.

    Examples
    ========

    >>> descent(3, 1)  # x**2 = 3*y**2 + z**2
    (1, 0, 1)

    `(x, y, z) = (1, 0, 1)` is a solution to the above equation.

    >>> descent(41, -113)
    (-16, -3, 1)

    References
    ==========

    * Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
      Mathematics of Computation, Volume 00, Number 0.

    """
    if abs(A) > abs(B):
        x, y, z = descent(B, A)
        return x, z, y

    if B == 1:
        return 1, 0, 1
    if A == 1:
        return 1, 1, 0
    if B == -A:
        return 0, 1, 1
    if B == A:
        x, z, y = descent(-1, A)
        return A*y, z, x

    w = sqrt_mod(A, B)
    x_0, z_0 = gaussian_reduce(w, A, B)

    t = (x_0**2 - A*z_0**2) // B
    t_2 = square_factor(t)
    t_1 = t // t_2**2

    x_1, z_1, y_1 = descent(A, t_1)

    return _remove_gcd(x_0*x_1 + A*z_0*z_1, z_0*x_1 + x_0*z_1, t_1*t_2*y_1)


def gaussian_reduce(w, a, b):
    r"""
    Returns a reduced solution `(x, z)` to the congruence
    `X^2 - aZ^2 \equiv 0 \ (mod \ b)` so that `x^2 + |a|z^2` is minimal.

    Here ``w`` is a solution of the congruence `x^2 \equiv a \ (mod \ b)`

    References
    ==========

    * Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
      Mathematics of Computation, Volume 00, Number 0.

    """
    u = (0, 1)
    v = (1, 0)

    if dot(u, v, w, a, b) < 0:
        v = (-v[0], -v[1])

    if norm(u, w, a, b) < norm(v, w, a, b):
        u, v = v, u

    while norm(u, w, a, b) > norm(v, w, a, b):
        k = dot(u, v, w, a, b) // dot(v, v, w, a, b)
        u, v = v, (u[0] - k*v[0], u[1] - k*v[1])

    u, v = v, u

    if dot(u, v, w, a, b) < dot(v, v, w, a, b)/2 or norm((u[0]-v[0], u[1]-v[1]), w, a, b) > norm(v, w, a, b):
        c = v
    else:
        c = (u[0] - v[0], u[1] - v[1])

    return c[0]*w + b*c[1], c[0]


def dot(u, v, w, a, b):
    r"""
    Returns a special dot product of the vectors `u = (u_{1}, u_{2})` and
    `v = (v_{1}, v_{2})` which is defined in order to reduce solution of
    the congruence equation `X^2 - aZ^2 \equiv 0 \ (mod \ b)`.

    """
    u_1, u_2 = u
    v_1, v_2 = v
    return (w*u_1 + b*u_2)*(w*v_1 + b*v_2) + abs(a)*u_1*v_1


def norm(u, w, a, b):
    r"""
    Returns the norm of the vector `u = (u_{1}, u_{2})` under the dot product
    defined by `u \cdot v = (wu_{1} + bu_{2})(w*v_{1} + bv_{2}) + |a|*u_{1}*v_{1}`
    where `u = (u_{1}, u_{2})` and `v = (v_{1}, v_{2})`.

    """
    u_1, u_2 = u
    return sqrt(dot((u_1, u_2), (u_1, u_2), w, a, b))


def holzer(x, y, z, a, b, c):
    r"""
    Simplify the solution `(x, y, z)` of the equation
    `ax^2 + by^2 = cz^2` with `a, b, c > 0` and `z^2 \geq \mid ab \mid` to
    a new reduced solution `(x', y', z')` such that `z'^2 \leq \mid ab \mid`.

    The algorithm is an interpretation of Mordell's reduction as described
    on page 8 of Cremona and Rusin's paper and the work of Mordell.

    References
    ==========

    * Efficient Solution of Rational Conices, J. E. Cremona and D. Rusin,
      Mathematics of Computation, Volume 00, Number 0.
    * Diophantine Equations, L. J. Mordell, page 48.

    """
    if _odd(c):
        k = 2*c
    else:
        k = c//2

    small = a*b*c
    step = 0
    while True:
        t1, t2, t3 = a*x**2, b*y**2, c*z**2
        # check that it's a solution
        if t1 + t2 != t3:
            if step == 0:
                raise ValueError('bad starting solution')
            break
        x_0, y_0, z_0 = x, y, z
        if max(t1, t2, t3) <= small:
            # Holzer condition
            break

        uv = u, v = base_solution_linear(k, y_0, -x_0)
        if None in uv:
            break

        p, q = -(a*u*x_0 + b*v*y_0), c*z_0
        r = Rational(p, q)
        if _even(c):
            w = _nint_or_floor(p, q)
            assert abs(w - r) <= Rational(1, 2)
        else:
            w = p//q  # floor
            if _odd(a*u + b*v + c*w):
                w += 1
            assert abs(w - r) <= 1

        A = a*u**2 + b*v**2 + c*w**2
        B = a*u*x_0 + b*v*y_0 + c*w*z_0
        x = Rational(x_0*A - 2*u*B, k)
        y = Rational(y_0*A - 2*v*B, k)
        z = Rational(z_0*A - 2*w*B, k)
        assert all(i.is_Integer for i in (x, y, z))
        step += 1

    return tuple(int(i) for i in (x_0, y_0, z_0))


def diop_general_pythagorean(eq, param=symbols('m', integer=True)):
    """
    Solves the general pythagorean equation,
    `a_{1}^2x_{1}^2 + a_{2}^2x_{2}^2 + . . . + a_{n}^2x_{n}^2 - a_{n + 1}^2x_{n + 1}^2 = 0`.

    Returns a tuple which contains a parametrized solution to the equation,
    sorted in the same order as the input variables.

    Parameters
    ==========

    eq : Expr
        is a general pythagorean equation which is assumed to be zero
    param : Symbol, optional
       is the base parameter used to construct other parameters by subscripting.

    Examples
    ========

    >>> from diofant.abc import e
    >>> diop_general_pythagorean(a**2 + b**2 + c**2 - d**2)
    (m1**2 + m2**2 - m3**2, 2*m1*m3, 2*m2*m3, m1**2 + m2**2 + m3**2)
    >>> diop_general_pythagorean(9*a**2 - 4*b**2 + 16*c**2 + 25*d**2 + e**2)
    (10*m1**2  + 10*m2**2  + 10*m3**2 - 10*m4**2, 15*m1**2  + 15*m2**2  + 15*m3**2  + 15*m4**2, 15*m1*m4, 12*m2*m4, 60*m3*m4)

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type == 'general_pythagorean':
        return _diop_general_pythagorean(var, coeff, param)


def _diop_general_pythagorean(var, coeff, t):

    if sign(coeff[var[0]**2]) + sign(coeff[var[1]**2]) + sign(coeff[var[2]**2]) < 0:
        for key in list(coeff):
            coeff[key] = -coeff[key]

    n = len(var)
    index = 0

    for i, v in enumerate(var):
        if sign(coeff[v**2]) == -1:
            index = i

    m = symbols(f'{t}1:{n:d}', integer=True)
    ith = sum(m_i**2 for m_i in m)
    L = [ith - 2*m[n - 2]**2]
    L.extend([2*m[i]*m[n-2] for i in range(n - 2)])
    sol = L[:index] + [ith] + L[index:]

    lcm = 1
    for i, v in enumerate(var):
        if i == index or (index > 0 and i == 0) or (index == 0 and i == 1):
            lcm = math.lcm(lcm, sqrt(abs(coeff[v**2])))
        else:
            s = sqrt(coeff[v**2])
            lcm = math.lcm(lcm, s if _odd(s) else s//2)

    for i, v in enumerate(var):
        sol[i] = (lcm*sol[i]) / sqrt(abs(coeff[v**2]))

    return tuple(sol)


def diop_general_sum_of_squares(eq, limit=1):
    r"""
    Solves the equation `x_{1}^2 + x_{2}^2 + . . . + x_{n}^2 - k = 0`.

    Returns at most ``limit`` number of solutions.

    Parameters
    ==========

    eq : Expr
        is an expression which is assumed to be zero. Also, ``eq`` should be in the form,
        `x_{1}^2 + x_{2}^2 + . . . + x_{n}^2 - k = 0`.
    limit : int, optional
        upper limit (the default is 1) for number of solutions returned.

    Notes
    =====

    When `n = 3` if `k = 4^a(8m + 7)` for some `a, m \in Z` then there will be
    no solutions.

    Examples
    ========

    >>> from diofant.abc import e
    >>> diop_general_sum_of_squares(a**2 + b**2 + c**2 + d**2 + e**2 - 2345)
    {(15, 22, 22, 24, 24)}

    References
    ==========

    * Representing an integer as a sum of three squares, [online],
      Available:
      https://proofwiki.org/wiki/Integer_as_Sum_of_Three_Squares

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type == 'general_sum_of_squares':
        return _diop_general_sum_of_squares(var, -int(coeff[1]), limit)


def _diop_general_sum_of_squares(var, k, limit=1):
    # solves Eq(sum(i**2 for i in var), k)
    n = len(var)
    if n < 3:
        raise ValueError('n must be greater than 2')

    s = set()

    if k < 0 or limit < 1:
        return s

    sign = [-1 if x.is_nonpositive else 1 for x in var]
    negs = sign.count(-1) != 0

    took = 0
    for t in sum_of_squares(k, n, zeros=True):
        if negs:
            s.add(tuple(sign[i]*j for i, j in enumerate(t)))
        else:
            s.add(t)
        took += 1
        if took == limit:
            break
    return s


def diop_general_sum_of_even_powers(eq, limit=1):
    """
    Solves the equation `x_{1}^e + x_{2}^e + . . . + x_{n}^e - k = 0`
    where `e` is an even, integer power.

    Returns at most ``limit`` number of solutions.

    Parameters
    ==========

    eq : Expr
        An expression which is assumed to be zero.  Also, ``eq`` should
        be in the form, `x_{1}^e + x_{2}^e + . . . + x_{n}^e - k = 0`.
    limit : Expr, optional
        Limit number of returned solutions.  Default is 1.

    Examples
    ========

    >>> diop_general_sum_of_even_powers(a**4 + b**4 - (2**4 + 3**4))
    {(2, 3)}

    See Also
    ========

    diofant.solvers.diophantine.power_representation

    """
    var, coeff, diop_type = classify_diop(eq, _dict=False)

    if diop_type == 'general_sum_of_even_powers':
        for k in coeff:
            if k.is_Pow and coeff[k]:
                p = k.exp
        return _diop_general_sum_of_even_powers(var, p, -coeff[1], limit)


def _diop_general_sum_of_even_powers(var, p, n, limit=1):
    # solves Eq(sum(i**2 for i in var), n)
    k = len(var)

    s = set()

    if n < 0 or limit < 1:
        return s

    sign = [-1 if x.is_nonpositive else 1 for x in var]
    negs = sign.count(-1) != 0

    took = 0
    for t in power_representation(n, p, k):
        if negs:
            s.add(tuple(sign[i]*j for i, j in enumerate(t)))
        else:
            s.add(t)
        took += 1
        if took == limit:
            break
    return s


# Functions below this comment can be more suitably grouped under
# an Additive number theory module rather than the Diophantine
# equation module.


def partition(n, k=None, zeros=False):
    """
    Returns a generator that can be used to generate partitions of an integer
    `n`.

    A partition of `n` is a set of positive integers which add up to `n`. For
    example, partitions of 3 are 3, 1 + 2, 1 + 1 + 1. A partition is returned
    as a tuple. If ``k`` equals None, then all possible partitions are returned
    irrespective of their size, otherwise only the partitions of size ``k`` are
    returned. If the ``zero`` parameter is set to True then a suitable
    number of zeros are added at the end of every partition of size less than
    ``k``.

    Parameters
    ==========

    n : int
        is a positive integer
    k : int, optional
        is the size of the partition which is also positive
        integer.  The default is ``None``.
    zeros : boolean, optional
        parameter is considered only if ``k`` is not ``None``. When the
        partitions are over, the last `next()` call throws the ``StopIteration``
        exception, so this function should always be used inside a
        try - except block.

    Examples
    ========

    >>> f = partition(5)
    >>> next(f)
    (1, 1, 1, 1, 1)
    >>> next(f)
    (1, 1, 1, 2)
    >>> g = partition(5, 3)
    >>> next(g)
    (1, 1, 3)
    >>> next(g)
    (1, 2, 2)
    >>> g = partition(5, 3, zeros=True)
    >>> next(g)
    (0, 0, 5)

    """
    from diofant.utilities.iterables import ordered_partitions
    if not zeros or k is None:
        for i in ordered_partitions(n, k):
            yield tuple(i)
    else:
        for m in range(1, k + 1):
            for i in ordered_partitions(n, m):
                i = tuple(i)
                yield (0,)*(k - len(i)) + i


def prime_as_sum_of_two_squares(p):
    """
    Represent a prime `p` as a unique sum of two squares; this can
    only be done if the prime is congruent to 1 mod 4.

    Examples
    ========

    >>> prime_as_sum_of_two_squares(7)  # can't be done
    >>> prime_as_sum_of_two_squares(5)
    (1, 2)

    References
    ==========

    * Representing a number as a sum of four squares, [online],
      Available: https://schorn.ch/lagrange.html

    See Also
    ========

    diofant.solvers.diophantine.sum_of_squares

    """
    if not p % 4 == 1:
        return

    if p % 8 == 5:
        b = 2
    else:
        b = 3

        while pow(b, (p - 1) // 2, p) == 1:
            b = nextprime(b)

    b = pow(b, (p - 1) // 4, p)
    a = p

    while b**2 > p:
        a, b = b, a % b

    return a % b, b


def sum_of_three_squares(n):
    r"""
    Returns a 3-tuple `(a, b, c)` such that `a^2 + b^2 + c^2 = n` and
    `a, b, c \geq 0`.

    Returns None if `n = 4^a(8m + 7)` for some `a, m \in Z`.

    Parameters
    ==========

    n : int
        a non-negative integer.

    Examples
    ========

    >>> sum_of_three_squares(44542)
    (18, 37, 207)

    References
    ==========

    * Representing a number as a sum of three squares, [online],
      Available: https://schorn.ch/lagrange.html

    See Also
    ========

    diofant.solvers.diophantine.sum_of_squares

    """
    special = {1: (1, 0, 0), 2: (1, 1, 0), 3: (1, 1, 1), 10: (1, 3, 0), 34: (3, 3, 4), 58: (3, 7, 0),
               85: (6, 7, 0), 130: (3, 11, 0), 214: (3, 6, 13), 226: (8, 9, 9), 370: (8, 9, 15),
               526: (6, 7, 21), 706: (15, 15, 16), 730: (1, 27, 0), 1414: (6, 17, 33), 1906: (13, 21, 36),
               2986: (21, 32, 39), 9634: (56, 57, 57)}

    v = 0

    if n == 0:
        return 0, 0, 0

    v = multiplicity(4, n)
    n //= 4**v

    if n % 8 == 7:
        return

    if n in special:
        x, y, z = special[n]
        return _sorted_tuple(2**v*x, 2**v*y, 2**v*z)

    s, _exact = integer_nthroot(n, 2)

    if _exact:
        return 2**v*s, 0, 0

    x = None

    if n % 8 == 3:
        s = s if _odd(s) else s - 1

        for x in range(s, -1, -2):  # pragma: no branch
            N = (n - x**2) // 2
            if isprime(N):
                y, z = prime_as_sum_of_two_squares(N)
                return _sorted_tuple(2**v*x, 2**v*(y + z), 2**v*abs(y - z))

    if n % 8 == 2 or n % 8 == 6:
        s = s if _odd(s) else s - 1
    else:
        s = s - 1 if _odd(s) else s

    for x in range(s, -1, -2):  # pragma: no branch
        N = n - x**2
        if isprime(N):
            y, z = prime_as_sum_of_two_squares(N)
            return _sorted_tuple(2**v*x, 2**v*y, 2**v*z)


def sum_of_four_squares(n):
    r"""
    Returns a 4-tuple `(a, b, c, d)` such that `a^2 + b^2 + c^2 + d^2 = n`.

    Here `a, b, c, d \geq 0`.

    Parameters
    ==========

    n : int
        is a non-negative integer.

    Examples
    ========

    >>> sum_of_four_squares(3456)
    (8, 8, 32, 48)
    >>> sum_of_four_squares(1294585930293)
    (0, 1234, 2161, 1137796)

    References
    ==========

    * Representing a number as a sum of four squares, [online],
      Available: https://schorn.ch/lagrange.html

    See Also
    ========

    diofant.solvers.diophantine.sum_of_squares

    """
    if n == 0:
        return 0, 0, 0, 0

    v = multiplicity(4, n)
    n //= 4**v

    if n % 8 == 7:
        d = 2
        n = n - 4
    elif n % 8 == 6 or n % 8 == 2:
        d = 1
        n = n - 1
    else:
        d = 0

    x, y, z = sum_of_three_squares(n)

    return _sorted_tuple(2**v*d, 2**v*x, 2**v*y, 2**v*z)


def power_representation(n, p, k, zeros=False):
    """
    Returns a generator for finding k-tuples of integers,
    `(n_{1}, n_{2}, ... n_{k})`, such that
    `n = n_{1}^p + n_{2}^p + ... n_{k}^p`.

    Parameters
    ==========

    n : int
        a non-negative integer
    k, p : int
        parameters to control representation ``n`` as a sum
        of ``k``, ``p``-th powers.
    zeros : boolean, optional
        if ``True`` (the default is ``False``), then the solutions
        will contain zeros.

    Examples
    ========

    Represent 1729 as a sum of two cubes:

    >>> f = power_representation(1729, 3, 2)
    >>> next(f)
    (9, 10)
    >>> next(f)
    (1, 12)

    If the flag `zeros` is True, the solution may contain tuples with
    zeros; any such solutions will be generated after the solutions
    without zeros:

    >>> list(power_representation(125, 2, 3, zeros=True))
    [(5, 6, 8), (3, 4, 10), (0, 5, 10), (0, 2, 11)]

    For even `p` the `permute_sign` function can be used to get all
    signed values:

    >>> from diofant.utilities.iterables import permute_signs
    >>> list(permute_signs((1, 12)))
    [(1, 12), (-1, 12), (1, -12), (-1, -12)]

    All possible signed permutations can also be obtained:

    >>> from diofant.utilities.iterables import signed_permutations
    >>> list(signed_permutations((1, 12)))
    [(1, 12), (-1, 12), (1, -12), (-1, -12), (12, 1), (-12, 1),
     (12, -1), (-12, -1)]

    """
    n, p, k = map(as_int, (n, p, k))

    if n < 0:
        if p % 2:
            for t in power_representation(-n, p, k, zeros):
                yield tuple(-i for i in t)
        return

    if p < 1 or k < 1:
        raise ValueError(filldedent(f"""
    Expecting positive integers for `(p, k)`, but got `({p}, {k})`"""))

    if n == 0:
        if zeros:
            yield (0,)*k
        return

    if k == 1:
        if p == 1:
            yield n,
        else:
            be = perfect_power(n)
            if be:
                b, e = be
                d, r = divmod(e, p)
                if not r:
                    yield b**d,
        return

    if p == 1:
        for t in partition(n, k, zeros=zeros):
            yield t
        return

    if p == 2:
        feasible = _can_do_sum_of_squares(n, k)
        if not feasible:
            return
        if not zeros and n > 33 and n >= k >= 5 and n - k in (13, 10, 7,
                                                              5, 4, 2, 1):
            # Todd G. Will, "When Is n^2 a Sum of k Squares?", [online].
            # Available: https://www.maa.org/sites/default/files/Will-MMz-201037918.pdf
            return
        if feasible == 2:  # it's prime and k == 2
            yield prime_as_sum_of_two_squares(n)
            return

    if k == 2 and p > 2:
        be = perfect_power(n)
        if be and be[1] % p == 0:
            return  # Fermat: a**n + b**n = c**n has no solution for n > 2

    if n >= k:
        a = integer_nthroot(n - (k - 1), p)[0]
        for t in pow_rep_recursive(a, k, n, [], p):
            yield tuple(reversed(t))

    if zeros:
        a = integer_nthroot(n, p)[0]
        for i in range(1, k):
            for t in pow_rep_recursive(a, i, n, [], p):
                yield tuple(reversed(t + (0,)*(k - i)))


sum_of_powers = power_representation


def pow_rep_recursive(n_i, k, n_remaining, terms, p):

    if k == 0 and n_remaining == 0:
        yield tuple(terms)
    else:
        if n_i >= 1 and k > 0:
            yield from pow_rep_recursive(n_i - 1, k, n_remaining, terms, p)
            residual = n_remaining - pow(n_i, p)
            if residual >= 0:
                yield from pow_rep_recursive(n_i, k - 1, residual, terms + [n_i], p)


def sum_of_squares(n, k, zeros=False):
    """Return a generator that yields the k-tuples of nonnegative
    values, the squares of which sum to n. If zeros is False (default)
    then the solution will not contain zeros. The nonnegative
    elements of a tuple are sorted.

    * If k == 1 and n is square, (n,) is returned.

    * If k == 2 then n can only be written as a sum of squares if
      every prime in the factorization of n that has the form
      4*k + 3 has an even multiplicity. If n is prime then
      it can only be written as a sum of two squares if it is
      in the form 4*k + 1.

    * if k == 3 then n can be written as a sum of squares if it does
      not have the form 4**m*(8*k + 7).

    * all integers can be written as the sum of 4 squares.

    * if k > 4 then n can be partitioned and each partition can
      be written as a sum of 4 squares; if n is not evenly divisible
      by 4 then n can be written as a sum of squares only if the
      an additional partition can be written as as sum of squares.
      For example, if k = 6 then n is partitioned into two parts,
      the first being written as a sum of 4 squares and the second
      being written as a sum of 2 squares -- which can only be
      done if the contition above for k = 2 can be met, so this will
      automatically reject certain partitions of n.

    Examples
    ========

    >>> list(sum_of_squares(25, 2))
    [(3, 4)]
    >>> list(sum_of_squares(25, 2, True))
    [(3, 4), (0, 5)]
    >>> list(sum_of_squares(25, 4))
    [(1, 2, 2, 4)]

    See Also
    ========

    diofant.utilities.iterables.signed_permutations

    """
    yield from power_representation(n, 2, k, zeros)


def _can_do_sum_of_squares(n, k):
    """Return True if n can be written as the sum of k squares,
    False if it can't, or 1 if k == 2 and n is prime (in which
    case it *can* be written as a sum of two squares). A False
    is returned only if it can't be written as k-squares, even
    if 0s are allowed.

    """
    if k < 1:
        return False
    if n < 0:
        return False
    if n == 0:
        return True
    if k == 1:
        return is_square(n)
    if k == 2:
        if n in (1, 2):
            return True
        if isprime(n):
            if n % 4 == 1:
                return 2  # signal that it was prime
            return False
        else:
            f = factorint(n)
            for p, m in f.items():
                # we can proceed iff no prime factor in the form 4*k + 3
                # has an odd multiplicity
                if (p % 4 == 3) and m % 2:
                    return False
            return True
    if k == 3:
        if (n//4**multiplicity(4, n)) % 8 == 7:
            return False
    # every number can be written as a sum of 4 squares; for k > 4 partitions
    # can be 0
    return True
