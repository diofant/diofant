"""
This module contain solvers for all kinds of equations,
algebraic or transcendental.
"""

import warnings
from collections import defaultdict
from types import GeneratorType

from ..core import (Add, Derivative, Dummy, Equality, Expr, Float, Function,
                    Ge, Integer, Lambda, Mul, Pow, S, Symbol, expand_log,
                    expand_mul, expand_multinomial, expand_power_exp, nfloat,
                    preorder_traversal, sympify)
from ..core.assumptions import check_assumptions
from ..core.compatibility import (default_sort_key, is_sequence, iterable,
                                  ordered)
from ..core.function import AppliedUndef
from ..core.relational import Relational
from ..functions import (Abs, Max, Min, Piecewise, acos, arg, asin, atan,
                         atan2, cos, exp, im, log, piecewise_fold, re, sin,
                         sqrt, tan)
from ..functions.elementary.trigonometric import (HyperbolicFunction,
                                                  TrigonometricFunction)
from ..integrals import Integral
from ..matrices import Matrix, zeros
from ..polys import Poly, RootOf, cancel, factor, roots, together
from ..polys.polyerrors import GeneratorsNeeded, PolynomialError
from ..simplify import (collect, denom, logcombine, nsimplify, posify,
                        powdenest, powsimp, simplify)
from ..simplify.fu import TR1
from ..simplify.sqrtdenest import unrad
from ..utilities import filldedent, subsets
from ..utilities.iterables import uniq
from .polysys import solve_linear_system, solve_poly_system


def denoms(eq, symbols=None):
    """Return (recursively) set of all denominators that appear in eq
    that contain any symbol in iterable ``symbols``; if ``symbols`` is
    None (default) then all denominators will be returned.

    Examples
    ========

    >>> from diofant.abc import x, y, z

    >>> denoms(x/y)
    {y}
    >>> denoms(x/(y*z))
    {y, z}
    >>> denoms(3/x + y/z)
    {x, z}
    >>> denoms(x/2 + y/z)
    {2, z}
    """

    pot = preorder_traversal(eq)
    dens = set()
    for p in pot:
        den = denom(p)
        if den is S.One:
            continue
        for d in Mul.make_args(den):
            dens.add(d)
    if not symbols:
        return dens
    rv = []
    for d in dens:
        free = d.free_symbols
        if any(s in free for s in symbols):
            rv.append(d)
    return set(rv)


def checksol(f, symbol, sol=None, **flags):
    """Checks whether sol is a solution of equation f == 0.

    Input can be either a single symbol and corresponding value
    or a dictionary of symbols and values. When given as a dictionary
    and flag ``simplify=True``, the values in the dictionary will be
    simplified. ``f`` can be a single equation or an iterable of equations.
    A solution must satisfy all equations in ``f`` to be considered valid;
    if a solution does not satisfy any equation, False is returned; if one or
    more checks are inconclusive (and none are False) then None
    is returned.

    Examples
    ========

    >>> from diofant import symbols
    >>> from diofant.solvers import checksol
    >>> x, y = symbols('x,y')
    >>> checksol(x**4 - 1, x, 1)
    True
    >>> checksol(x**4 - 1, x, 0)
    False
    >>> checksol(x**2 + y**2 - 5**2, {x: 3, y: 4})
    True

    To check if an expression is zero using checksol, pass it
    as ``f`` and send an empty dictionary for ``symbol``:

    >>> checksol(x**2 + x - x*(x + 1), {})
    True

    None is returned if checksol() could not conclude.

    flags:
        'numerical=True (default)'
           do a fast numerical check if ``f`` has only one symbol.
        'minimal=True (default is False)'
           a very fast, minimal testing.
        'warn=True (default is False)'
           show a warning if checksol() could not conclude.
        'simplify=True (default)'
           simplify solution before substituting into function and
           simplify the function before trying specific simplifications
        'force=True (default is False)'
           make positive all symbols without assumptions regarding sign.
    """
    minimal = flags.get('minimal', False)

    if sol is not None:
        sol = {symbol: sol}
    elif isinstance(symbol, dict):
        sol = symbol
    else:
        raise ValueError("Expecting (sym, val) or ({sym: val}, "
                         "None) but got (%s, %s)" % (symbol, sol))

    if sol and not f.has(*list(sol.keys())):
        # if f(y) == 0, x=3 does not set f(y) to zero...nor does it not
        if f.is_Number:
            return f.is_zero
        else:
            return

    illegal = {S.NaN,
               S.ComplexInfinity,
               S.Infinity,
               S.NegativeInfinity}
    if any(sympify(v).atoms() & illegal for k, v in sol.items()):
        return False

    was = f
    attempt = -1
    numerical = flags.get('numerical', True)
    while 1:
        attempt += 1
        if attempt == 0:
            val = f.subs(sol)
            if val.atoms() & illegal:
                return False
        elif attempt == 1:
            if val.free_symbols:
                if not val.is_constant(*list(sol.keys()), simplify=not minimal):
                    return False
                # there are free symbols -- simple expansion might work
                _, val = val.as_content_primitive()
                val = expand_mul(expand_multinomial(val))
        elif attempt == 2:
            if minimal:
                return
            if flags.get('simplify', True):
                for k in sol:
                    sol[k] = simplify(sol[k])
            # start over without the failed expanded form, possibly
            # with a simplified solution
            val = f.subs(sol)
            if flags.get('force', True):
                val, reps = posify(val)
                # expansion may work now, so try again and check
                exval = expand_mul(expand_multinomial(val))
                if exval.is_number or not exval.free_symbols:
                    # we can decide now
                    val = exval
        elif attempt == 3:
            val = powsimp(val)
        elif attempt == 4:
            val = cancel(val)
        elif attempt == 5:
            val = val.expand()
        elif attempt == 6:
            val = together(val)
        elif attempt == 7:
            val = powsimp(val)
        else:
            # if there are no radicals and no functions then this can't be
            # zero anymore -- can it?
            pot = preorder_traversal(expand_mul(val))
            seen = set()
            saw_pow_func = False
            for p in pot:
                if p in seen:
                    continue
                seen.add(p)
                if p.is_Pow and not p.exp.is_Integer:
                    saw_pow_func = True
                elif p.is_Function:
                    saw_pow_func = True
                if saw_pow_func:
                    break
            if saw_pow_func is False:
                return False
            if flags.get('force', True):
                # don't do a zero check with the positive assumptions in place
                val = val.subs(reps)
            break

        if val == was:
            continue
        elif val.is_Rational:
            return val == 0
        elif val.is_nonzero:
            return False
        if numerical and not val.free_symbols:
            return bool(abs(val.n(18).n(12, chop=True)) < 1e-9)
        was = val

    if flags.get('warn', False):
        warnings.warn("\n\tWarning: could not verify solution %s." % sol)


def solve(f, *symbols, **flags):
    r"""Algebraically solves equation or system of equations.

    Parameters
    ==========

    f : Expr, Equality or iterable of above
        All expressions are assumed to be equal to 0.

    \*symbols : tuple
        If none symbols given (empty tuple), free symbols
        of expressions will be used.

    \*\*flags : dict
        A dictionary of following parameters:

        exclude : iterable, optional
            Don't try to solve for any of the symbols in
            exclude.  Default is [].
        check : bool, optional
            If False, don't do any testing of solutions.  Default is
            True, i.e. the solutions are checked and those that doesn't
            satisfy given assumptions on symbols solved for or make any
            denominator zero - are automatically excluded.
        numerical : bool, optional
            If enabled (default), do a fast numerical check
            if ``f`` has only one symbol.
        minimal : bool, optional
            A very fast, minimal testing.  Default is False.
        warn : bool, optional
            Show a warning if :func:`~diofant.solvers.solvers.checksol`
            could not conclude.  Default is False.
        simplify : bool, optional
            Enable simplification (default) for all but polynomials of
            order 3 or greater before returning them and (if check is
            not False) use the general simplify function on the solutions
            and the expression obtained when they are substituted into the
            function which should be zero.
        force : bool, optional
            Make positive all symbols without assumptions regarding
            sign.  Default is False.
        rational : bool or None, optional
            If True (default), recast Floats as Rational.  If None,
            Floats will be recast as rationals but the answer will be
            recast as Floats.  If the flag is False then nothing
            will be done to the Floats.
        cubics, quartics, quintics : bool, optional
            Return explicit solutions (with radicals, which can be quite
            long) when, respectively, cubic, quartic or quintic expressions
            are encountered.  Default is True.  If False,
            :class:`~diofant.polys.rootoftools.RootOf` instances will
            be returned instead.

    Examples
    ========

    The output varies according to the input and can be seen by example::

        >>> from diofant import solve, Poly, Eq, Function, exp, I, sqrt
        >>> from diofant.abc import x, y, z, a, b
        >>> f = Function('f')

    * single expression and single symbol that is in the expression

        >>> solve(x - y, x)
        [{x: y}]
        >>> solve(x - 3, x)
        [{x: 3}]
        >>> solve(Eq(x, 3), x)
        [{x: 3}]
        >>> solve(Poly(x - 3), x)
        [{x: 3}]

    * single expression with no symbol that is in the expression

        >>> solve(3, x)
        []
        >>> solve(x - 3, y)
        []

    * single expression with no symbol given

          In this case, all free symbols will be selected as potential
          symbols to solve for. If the equation is univariate then a list
          of solutions is returned; otherwise -- as is the case when symbols are
          given as an iterable of length > 1 -- a list of mappings will be returned.

            >>> solve(x - 3)
            [{x: 3}]
            >>> solve(x**2 - y**2)
            [{x: -y}, {x: y}]
            >>> solve(z**2*x**2 - z**2*y**2)
            [{x: -y}, {x: y}, {z: 0}]
            >>> solve(z**2*x - z**2*y**2)
            [{x: y**2}, {z: 0}]

    * when an object other than a Symbol is given as a symbol, it is
      isolated algebraically and an implicit solution may be obtained.
      This is mostly provided as a convenience to save one from replacing
      the object with a Symbol and solving for that Symbol. It will only
      work if the specified object can be replaced with a Symbol using the
      subs method.

          >>> solve(f(x) - x, f(x))
          [{f(x): x}]
          >>> solve(f(x).diff(x) - f(x) - x, f(x).diff(x))
          [{Derivative(f(x), x): x + f(x)}]
          >>> solve(f(x).diff(x) - f(x) - x, f(x))
          [{f(x): -x + Derivative(f(x), x)}]

          >>> from diofant import Indexed, IndexedBase, Tuple, sqrt
          >>> A = IndexedBase('A')
          >>> eqs = Tuple(A[1] + A[2] - 3, A[1] - A[2] + 1)
          >>> solve(eqs, eqs.atoms(Indexed))
          [{A[1]: 1, A[2]: 2}]

        * It is possible to solve for anything that can be targeted with
          subs:

            >>> solve(x + 2 + sqrt(3), x + 2)
            [{x + 2: -sqrt(3)}]
            >>> solve((x + 2 + sqrt(3), x + 4 + y), y, x + 2)
            [{y: -2 + sqrt(3), x + 2: -sqrt(3)}]

        * Nothing heroic is done in this implicit solving so you may end up
          with a symbol still in the solution:

            >>> eqs = (x*y + 3*y + sqrt(3), x + 4 + y)
            >>> solve(eqs, y, x + 2)
            [{y: -sqrt(3)/(x + 3), x + 2: (-2*x - 6 + sqrt(3))/(x + 3)}]
            >>> solve(eqs, y*x, x)
            [{x: -y - 4, x*y: -3*y - sqrt(3)}]

        * if you attempt to solve for a number remember that the number
          you have obtained does not necessarily mean that the value is
          equivalent to the expression obtained:

            >>> solve(sqrt(2) - 1, 1, check=False)
            [{1: sqrt(2)}]
            >>> solve(x - y + 1, 1)  # /!\ -1 is targeted, too
            [{1: x/(y - 1)}]
            >>> [_[1].subs(z, -1) for _ in solve((x - y + 1).subs(-1, z), 1)]
            [-x + y]

        * To solve for a function within a derivative, use dsolve.

    * single expression and more than 1 symbol

        * when there is a linear solution

            >>> solve(x - y**2, x, y)
            [{x: y**2}]
            >>> solve(x**2 - y, x, y)
            [{y: x**2}]

        * if there is no linear solution then the first successful
          attempt for a nonlinear solution will be returned

            >>> solve(x**2 - y**2, x, y)
            [{x: -y}, {x: y}]
            >>> solve(x**2 - y**2/exp(x), x, y)
            [{x: 2*LambertW(y/2)}]
            >>> solve(x**2 - y**2/exp(x), y, x)
            [{y: -x*sqrt(E**x)}, {y: x*sqrt(E**x)}]

    * iterable of one or more of the above

        * when the system is linear

            * with a solution

                >>> solve([x - 3], x)
                [{x: 3}]
                >>> solve((x + 5*y - 2, -3*x + 6*y - 15), x, y)
                [{x: -3, y: 1}]
                >>> solve((x + 5*y - 2, -3*x + 6*y - 15), x, y, z)
                [{x: -3, y: 1}]
                >>> solve((x + 5*y - 2, -3*x + 6*y - z), z, x, y)
                [{x: -5*y + 2, z: 21*y - 6}]

            * without a solution

                >>> solve([x + 3, x - 3])
                []

        * if no symbols are given, all free symbols will be selected and a list
          of mappings returned

            >>> solve([x - 2, x**2 + y])
            [{x: 2, y: -4}]
            >>> solve([x - 2, x**2 + f(x)], {f(x), x})
            [{x: 2, f(x): -4}]

    Notes
    =====

    solve() with check=True (default) will run through the symbol tags to
    eliminate unwanted solutions.  If no assumptions are included all possible
    solutions will be returned.

        >>> from diofant import Symbol, solve
        >>> x = Symbol("x")
        >>> solve(x**2 - 1)
        [{x: -1}, {x: 1}]

    By using the positive tag only one solution will be returned:

        >>> pos = Symbol("pos", positive=True)
        >>> solve(pos**2 - 1)
        [{pos: 1}]

    When the solutions are checked, those that make any denominator zero
    are automatically excluded. If you do not want to exclude such solutions
    then use the check=False option:

        >>> from diofant import sin, limit
        >>> solve(sin(x)/x)  # 0 is excluded
        [{x: pi}]

    If check=False then a solution to the numerator being zero is found: x = 0.
    In this case, this is a spurious solution since sin(x)/x has the well known
    limit (without discontinuity) of 1 at x = 0:

        >>> solve(sin(x)/x, check=False)
        [{x: 0}, {x: pi}]

    In the following case, however, the limit exists and is equal to the the
    value of x = 0 that is excluded when check=True:

        >>> eq = x**2*(1/x - z**2/x)
        >>> solve(eq, x)
        []
        >>> solve(eq, x, check=False)
        [{x: 0}]
        >>> limit(eq, x, 0, '-')
        0
        >>> limit(eq, x, 0, '+')
        0

    See Also
    ========

    diofant.solvers.recurr.rsolve : solving recurrence equations
    diofant.solvers.ode.dsolve : solving differential equations
    diofant.solvers.inequalities.reduce_inequalities : solving inequalities
    """
    # keeping track of how f was passed since if it is a list
    # a dictionary of results will be returned.
    ###########################################################################

    def _sympified_list(w):
        return list(map(sympify, w if iterable(w) else [w]))
    bare_f = not iterable(f)
    ordered_symbols = (symbols and symbols[0] and
                       (isinstance(symbols[0], (Dummy, Symbol)) or
                        is_sequence(symbols[0], include=GeneratorType)))
    f, symbols = (_sympified_list(w) for w in [f, symbols])

    # preprocess equation(s)
    ###########################################################################
    for i, fi in enumerate(f):
        if isinstance(fi, Equality):
            if 'ImmutableMatrix' in [type(a).__name__ for a in fi.args]:
                f[i] = fi.lhs - fi.rhs
            else:
                f[i] = Add(fi.lhs, -fi.rhs, evaluate=False)
        elif isinstance(fi, Relational):
            raise ValueError("Only expressions or equalities "
                             "supported, got %s" % fi)
        elif isinstance(fi, Poly):
            f[i] = fi.as_expr()

        # rewrite hyperbolics in terms of exp
        f[i] = f[i].replace(lambda w: isinstance(w, HyperbolicFunction),
                            lambda w: w.rewrite(exp))

        # replace min/max:
        f[i] = f[i].replace(lambda w: isinstance(w, (Min, Max)),
                            lambda w: w.rewrite(Piecewise))

        # if we have a Matrix, we need to iterate over its elements again
        if f[i].is_Matrix:
            bare_f = False
            f.extend(list(f[i]))
            f[i] = S.Zero

        # if we can split it into real and imaginary parts then do so
        freei = f[i].free_symbols
        if freei and all(s.is_extended_real or s.is_imaginary for s in freei):
            fr, fi = f[i].as_real_imag()
            # accept as long as new re, im, arg or atan2 are not introduced
            had = f[i].atoms(re, im, arg, atan2)
            if fr and fi and fr != fi and not any(
                    i.atoms(re, im, arg, atan2) - had for i in (fr, fi)):
                if bare_f:
                    bare_f = False
                f[i: i + 1] = [fr, fi]

    # preprocess symbol(s)
    ###########################################################################
    if not symbols:
        # get symbols from equations
        symbols = set().union(*[fi.free_symbols for fi in f])
        if len(symbols) < len(f):
            for fi in f:
                pot = preorder_traversal(fi)
                for p in pot:
                    if not (p.is_number or p.is_Add or p.is_Mul) or \
                            isinstance(p, AppliedUndef):
                        symbols.add(p)
                        pot.skip()  # don't go any deeper
        symbols = list(symbols)
        # supply dummy symbols so solve(3) behaves like solve(3, x)
        for i in range(len(f) - len(symbols)):
            symbols.append(Dummy())

        ordered_symbols = False
    elif len(symbols) == 1 and iterable(symbols[0]):
        symbols = symbols[0]

    # remove symbols the user is not interested in
    exclude = flags.pop('exclude', set())
    if exclude:
        if isinstance(exclude, Expr):
            exclude = [exclude]
        exclude = set().union(*[e.free_symbols for e in sympify(exclude)])
    symbols = [s for s in symbols if s not in exclude]

    # real/imag handling -----------------------------
    w = Dummy('w')
    piece = Lambda(w, Piecewise((w, Ge(w, 0)), (-w, True)))
    for i, fi in enumerate(f):
        # Abs
        reps = []
        for a in fi.atoms(Abs):
            if not a.has(*symbols):
                continue
            if a.args[0].is_extended_real is None and a.args[0].is_imaginary is not True:
                raise NotImplementedError('solving %s when the argument '
                                          'is not real or imaginary.' % a)
            reps.append((a, piece(a.args[0]) if a.args[0].is_extended_real else
                         piece(a.args[0]*S.ImaginaryUnit)))
        fi = fi.subs(reps)

        # arg
        _arg = [a for a in fi.atoms(arg) if a.has(*symbols)]
        fi = fi.xreplace({a: atan(im(a.args[0])/re(a.args[0])) for a in _arg})

        # save changes
        f[i] = fi

    # see if re(s) or im(s) appear
    irf = []
    for s in symbols:
        if s.is_extended_real or s.is_imaginary:
            continue  # neither re(x) nor im(x) will appear
        # if re(s) or im(s) appear, the auxiliary equation must be present
        if any(fi.has(re(s), im(s)) for fi in f):
            irf.append((s, re(s) + S.ImaginaryUnit*im(s)))
    if irf:
        for s, rhs in irf:
            for i, fi in enumerate(f):
                f[i] = fi.xreplace({s: rhs})
            f.append(s - rhs)
            symbols.extend([re(s), im(s)])
        if bare_f:
            bare_f = False
    # end of real/imag handling  -----------------------------

    symbols = list(uniq(symbols))
    if not ordered_symbols:
        # we do this to make the results returned canonical in case f
        # contains a system of nonlinear equations; all other cases should
        # be unambiguous
        symbols = sorted(symbols, key=default_sort_key)

    # we can solve for non-symbol entities by replacing them with Dummy symbols
    symbols_new = []
    symbol_swapped = False
    for i, s in enumerate(symbols):
        if s.is_Symbol:
            s_new = s
        else:
            symbol_swapped = True
            s_new = Dummy('X%d' % i)
        symbols_new.append(s_new)

    if symbol_swapped:
        swap_sym = list(zip(symbols, symbols_new))
        f = [fi.subs(swap_sym) for fi in f]
        symbols = symbols_new
        swap_sym = {v: k for k, v in swap_sym}
    else:
        swap_sym = {}

    # this is needed in the next two events
    symset = set(symbols)

    # mask off any Object that we aren't going to invert: Derivative,
    # Integral, etc... so that solving for anything that they contain will
    # give an implicit solution
    seen = set()
    non_inverts = set()
    for fi in f:
        pot = preorder_traversal(fi)
        for p in pot:
            if not isinstance(p, Expr) or isinstance(p, Piecewise):
                pass
            elif (isinstance(p, bool) or not p.args or p in symset or
                  p.is_Add or p.is_Mul or p.is_Pow or p.is_Function or
                  isinstance(p, RootOf)) and p.func not in (re, im):
                continue
            elif p not in seen:
                seen.add(p)
                if p.free_symbols & symset:
                    non_inverts.add(p)
                else:
                    continue
            pot.skip()
    del seen
    non_inverts = {d: Dummy() for d in non_inverts}
    f = [fi.subs(non_inverts) for fi in f]

    non_inverts = [(v, k.subs(swap_sym)) for k, v in non_inverts.items()]

    # rationalize Floats
    floats = False
    if flags.get('rational', True) is not False:
        for i, fi in enumerate(f):
            if fi.has(Float):
                floats = True
                f[i] = nsimplify(fi, rational=True)

    # Any embedded piecewise functions need to be brought out to the
    # top level so that the appropriate strategy gets selected.
    # However, this is necessary only if one of the piecewise
    # functions depends on one of the symbols we are solving for.
    for i, fi in enumerate(f):
        if any(e.has(*symbols) for e in fi.atoms(Piecewise)):
            f[i] = piecewise_fold(fi)

    #
    # try to get a solution
    ###########################################################################
    if bare_f:
        solution = _solve(f[0], *symbols, **flags)
    else:
        solution = _solve_system(f, symbols, **flags)

    #
    # postprocessing
    ###########################################################################
    # Restore masked-off objects
    if non_inverts:

        def _do_dict(solution):
            return {k: v.subs(non_inverts) for k, v in solution.items()}
        for i in range(1):
            if type(solution) is dict:
                solution = _do_dict(solution)
                break
            elif solution and type(solution) is list:
                if type(solution[0]) is dict:
                    solution = [_do_dict(s) for s in solution]
                    break
                else:
                    solution = [v.subs(non_inverts) for v in solution]
                    break
            elif not solution:
                break
        else:  # pragma: no cover
            raise NotImplementedError(filldedent('''
                            no handling of %s was implemented''' % solution))

    # Restore original "symbols" if a dictionary is returned.
    # This is not necessary for
    #   - the single univariate equation case
    #     since the symbol will have been removed from the solution;
    #   - the nonlinear poly_system since that only supports zero-dimensional
    #     systems and those results come back as a list
    #
    # ** unless there were Derivatives with the symbols, but those were handled
    #    above.
    if symbol_swapped:
        symbols = [swap_sym[k] for k in symbols]
        if type(solution) is dict:
            solution = {swap_sym[k]: v.subs(swap_sym)
                        for k, v in solution.items()}
        elif solution and type(solution) is list and type(solution[0]) is dict:
            for i, sol in enumerate(solution):
                solution[i] = {swap_sym[k]: v.subs(swap_sym)
                               for k, v in sol.items()}

    # Get assumptions about symbols, to filter solutions.
    # Note that if assumptions about a solution can't be verified, it is still
    # returned.
    check = flags.get('check', True)

    # restore floats
    if floats and solution and flags.get('rational', None) is None:
        solution = nfloat(solution, exponent=False)

    if check and solution:  # assumption checking

        warn = flags.get('warn', False)
        got_None = []  # solutions for which one or more symbols gave None
        no_False = []  # solutions for which no symbols gave False
        if type(solution) is list:
            if type(solution[0]) is dict:
                for sol in solution:
                    a_None = False
                    for symb, val in sol.items():
                        test = check_assumptions(val, **symb._assumptions)
                        if test:
                            continue
                        if test is False:
                            break
                        a_None = True
                    else:
                        no_False.append(sol)
                        if a_None:
                            got_None.append(sol)
            else:  # list of expressions
                for sol in solution:
                    test = check_assumptions(sol, **symbols[0]._assumptions)
                    if test is False:
                        continue
                    no_False.append(sol)
                    if test is None:
                        got_None.append(sol)

        elif type(solution) is dict:
            a_None = False
            for symb, val in solution.items():
                test = check_assumptions(val, **symb._assumptions)
                if test:
                    continue
                if test is False:
                    no_False = None
                    break
                a_None = True
            else:
                no_False = solution
                if a_None:
                    got_None.append(solution)
        else:
            raise TypeError('Unrecognized solution')  # improve the checker

        solution = no_False
        if warn and got_None:
            warnings.warn(filldedent("""
                \tWarning: assumptions concerning following solution(s)
                can't be checked:""" + '\n\t' +
                                     ', '.join(str(s) for s in got_None)))

    #
    # done
    ###########################################################################

    if isinstance(solution, list):
        # Make sure that a list of solutions is ordered in a canonical way.
        solution.sort(key=default_sort_key)

    # return a list of mappings or []
    if not solution:
        solution = []
    else:
        if isinstance(solution, dict):
            solution = [solution]
        elif isinstance(solution[0], dict):
            pass
        else:
            if len(symbols) != 1:
                raise ValueError("Length should be 1")
            solution = [{symbols[0]: s} for s in solution]

    return solution


def _solve(f, *symbols, **flags):
    """Return a checked solution for f in terms of one or more of the
    symbols. A list should be returned except for the case when a linear
    undetermined-coefficients equation is encountered (in which case
    a dictionary is returned).

    If no method is implemented to solve the equation, a NotImplementedError
    will be raised. In the case that conversion of an expression to a Poly
    gives None a ValueError will be raised.
    """

    not_impl_msg = "No algorithms are implemented to solve equation %s"

    if len(symbols) != 1:
        soln = None
        free = f.free_symbols
        ex = free - set(symbols)
        if len(ex) != 1:
            ind, dep = f.as_independent(*symbols)
            ex = ind.free_symbols & dep.free_symbols
        # find first successful solution
        failed = []
        got_s = set()
        result = []
        for s in symbols:
            n, d = solve_linear(f, symbols=[s])
            if n.is_Symbol:
                # no need to check but we should simplify if desired
                if flags.get('simplify', True):
                    d = simplify(d)
                if got_s and any(ss in d.free_symbols for ss in got_s):
                    # sol depends on previously solved symbols: discard it
                    continue
                got_s.add(n)
                result.append({n: d})
            elif n and d:  # otherwise there was no solution for s
                failed.append(s)
        if not failed:
            return result
        for s in failed:
            try:
                soln = _solve(f, s, **flags)
                for sol in soln:
                    if got_s and any(ss in sol.free_symbols for ss in got_s):
                        # sol depends on previously solved symbols: discard it
                        continue
                    got_s.add(s)
                    result.append({s: sol})
            except NotImplementedError:
                continue
        if got_s:
            return result
        else:
            raise NotImplementedError(not_impl_msg % f)
    symbol = symbols[0]

    # /!\ capture this flag then set it to False so that no checking in
    # recursive calls will be done; only the final answer is checked
    checkdens = check = flags.pop('check', True)
    flags['check'] = False

    # build up solutions if f is a Mul
    if f.is_Mul:
        result = set()
        for m in f.args:
            soln = _solve(m, symbol, **flags)
            result.update(set(soln))
        result = list(result)
        if check:
            # all solutions have been checked but now we must
            # check that the solutions do not set denominators
            # in any factor to zero
            dens = denoms(f, symbols)
            result = [s for s in result if
                      all(not checksol(den, {symbol: s}, **flags) for den in
                          dens)]
        # set flags for quick exit at end
        check = False
        flags['simplify'] = False

    elif f.is_Piecewise:
        result = set()
        for n, (expr, cond) in enumerate(f.args):
            candidates = _solve(piecewise_fold(expr), *symbols, **flags)
            for candidate in candidates:
                if candidate in result:
                    continue
                try:
                    v = (cond == S.true) or cond.subs(symbol, candidate)
                except:
                    v = False
                if v != S.false:
                    # Only include solutions that do not match the condition
                    # of any previous pieces.
                    matches_other_piece = False
                    for other_n, (other_expr, other_cond) in enumerate(f.args):
                        if other_n == n:
                            break
                        if other_cond == S.false:
                            continue
                        try:
                            if other_cond.subs(symbol, candidate) == S.true:
                                matches_other_piece = True
                                break
                        except:
                            pass
                    if not matches_other_piece:
                        v = v == S.true or v.doit()
                        if isinstance(v, Relational):
                            v = v.canonical
                        result.add(Piecewise(
                            (candidate, v),
                            (S.NaN, True)
                        ))
        check = False
    else:
        # first see if it really depends on symbol and whether there
        # is a linear solution
        f_num, sol = solve_linear(f, symbols=symbols)
        if symbol not in f_num.free_symbols:
            return []
        elif f_num.is_Symbol:
            # no need to check but simplify if desired
            if flags.get('simplify', True):
                sol = simplify(sol)
            return [sol]

        result = False  # no solution was obtained
        msg = ''  # there is no failure message

        # Poly is generally robust enough to convert anything to
        # a polynomial and tell us the different generators that it
        # contains, so we will inspect the generators identified by
        # polys to figure out what to do.

        # try to identify a single generator that will allow us to solve this
        # as a polynomial, followed (perhaps) by a change of variables if the
        # generator is not a symbol

        try:
            poly = Poly(f_num)
            if poly is None:
                raise ValueError('could not convert %s to Poly' % f_num)
        except GeneratorsNeeded:
            simplified_f = simplify(f_num)
            if simplified_f != f_num:
                return _solve(simplified_f, symbol, **flags)
            raise ValueError('expression appears to be a constant')

        gens = [g for g in poly.gens if g.has(symbol)]

        def _as_base_q(x):
            """Return (b**e, q) for x = b**(p*e/q) where p/q is the leading
            Rational of the exponent of x, e.g. exp(-2*x/3) -> (exp(x), 3)
            """
            b, e = x.as_base_exp()
            if e.is_Rational:
                return b, e.q
            if not e.is_Mul:
                return x, 1
            c, ee = e.as_coeff_Mul()
            if c.is_Rational and c is not S.One:  # c could be a Float
                return b**ee, c.q
            return x, 1

        if len(gens) > 1:
            # If there is more than one generator, it could be that the
            # generators have the same base but different powers, e.g.
            #   >>> Poly(exp(x) + 1/exp(x))
            #   Poly(exp(-x) + exp(x), exp(-x), exp(x), domain='ZZ')
            #
            # If unrad was not disabled then there should be no rational
            # exponents appearing as in
            #   >>> Poly(sqrt(x) + sqrt(sqrt(x)))
            #   Poly(sqrt(x) + x**(1/4), sqrt(x), x**(1/4), domain='ZZ')

            bases, qs = list(zip(*[_as_base_q(g) for g in gens]))
            bases = set(bases)

            if len(bases) > 1 or not all(q == 1 for q in qs):
                funcs = {b for b in bases if b.is_Function}

                trig = {_ for _ in funcs if
                        isinstance(_, TrigonometricFunction)}
                other = funcs - trig
                if not other and len(funcs.intersection(trig)) > 1:
                    newf = TR1(f_num).rewrite(tan)
                    if newf != f_num:
                        result = _solve(newf, symbol, **flags)

                # just a simple case - see if replacement of single function
                # clears all symbol-dependent functions, e.g.
                # log(x) - log(log(x) - 1) - 3 can be solved even though it has
                # two generators.

                if result is False and funcs:
                    funcs = list(ordered(funcs))  # put shallowest function first
                    f1 = funcs[0]
                    t = Dummy('t')
                    # perform the substitution
                    ftry = f_num.subs(f1, t)

                    # if no Functions left, we can proceed with usual solve
                    if not ftry.has(symbol):
                        cv_sols = _solve(ftry, t, **flags)
                        cv_inv = _solve(t - f1, symbol, **flags)[0]
                        sols = []
                        for sol in cv_sols:
                            sols.append(cv_inv.subs(t, sol))
                        result = list(ordered(sols))

                if result is False:
                    msg = 'multiple generators %s' % gens

            else:
                # e.g. case where gens are exp(x), exp(-x)
                u = bases.pop()
                t = Dummy('t')
                inv = _solve(u - t, symbol, **flags)
                if isinstance(u, Pow):
                    # this will be resolved by factor in _tsolve but we might
                    # as well try a simple expansion here to get things in
                    # order so something like the following will work now without
                    # having to factor:
                    #
                    # >>> eq = (exp(I*(-x-2))+exp(I*(x+2)))
                    # >>> eq.subs(exp(x),y)  # fails
                    # exp(I*(-x - 2)) + exp(I*(x + 2))
                    # >>> eq.expand().subs(exp(x),y)  # works
                    # y**I*exp(2*I) + y**(-I)*exp(-2*I)
                    def _expand(p):
                        b, e = p.as_base_exp()
                        e = expand_mul(e)
                        return expand_power_exp(b**e)
                    ftry = f_num.replace(
                        lambda w: w.is_Pow,
                        _expand).subs(u, t)
                    if not ftry.has(symbol):
                        soln = _solve(ftry, t, **flags)
                        sols = []
                        for sol in soln:
                            for i in inv:
                                sols.append(i.subs(t, sol))
                        result = list(ordered(sols))

        elif len(gens) == 1:

            # There is only one generator that we are interested in, but
            # there may have been more than one generator identified by
            # polys (e.g. for symbols other than the one we are interested
            # in) so recast the poly in terms of our generator of interest.
            # Also use composite=True with f_num since Poly won't update
            # poly as documented in issue sympy/sympy#8810.

            poly = Poly(f_num, gens[0], composite=True)

            # if we aren't on the tsolve-pass, use roots
            if not flags.pop('tsolve', False):
                soln = None
                deg = poly.degree()
                flags['tsolve'] = True
                solvers = {k: flags.get(k, True) for k in
                           ('cubics', 'quartics', 'quintics')}
                soln = roots(poly, **solvers)
                if sum(soln.values()) < deg:
                    # e.g. roots(32*x**5 + 400*x**4 + 2032*x**3 +
                    #            5000*x**2 + 6250*x + 3189) -> {}
                    # so all_roots is used and RootOf instances are
                    # returned *unless* the system is multivariate
                    # or high-order EX domain.
                    try:
                        soln = poly.all_roots()
                    except NotImplementedError:  # pragma: no cover
                        if not flags.get('incomplete', True):
                                raise NotImplementedError(
                                    filldedent('''
    Neither high-order multivariate polynomials
    nor sorting of EX-domain polynomials is supported.
    If you want to see any results, pass keyword incomplete=True to
    solve; to see numerical values of roots
    for univariate expressions, use nroots.
    '''))
                        else:
                            pass
                else:
                    soln = list(soln.keys())

                if soln is not None:
                    u = poly.gen
                    if u != symbol:
                        try:
                            t = Dummy('t')
                            iv = _solve(u - t, symbol, **flags)
                            soln = list(ordered({i.subs(t, s) for i in iv for s in soln}))
                        except NotImplementedError:
                            # perhaps _tsolve can handle f_num
                            soln = None
                    else:
                        check = False  # only dens need to be checked
                    if soln is not None:
                        if len(soln) > 2:
                            # if the flag wasn't set then unset it since high-order
                            # results are quite long. Perhaps one could base this
                            # decision on a certain critical length of the
                            # roots. In addition, wester test M2 has an expression
                            # whose roots can be shown to be real with the
                            # unsimplified form of the solution whereas only one of
                            # the simplified forms appears to be real.
                            flags['simplify'] = flags.get('simplify', False)
                        result = soln

    # fallback if above fails
    # -----------------------
    if result is False:
        # try unrad
        if flags.pop('_unrad', True):
            u = unrad(f_num, symbol)
            if u:
                eq, cov = u
                if cov:
                    isym, ieq = cov
                    inv = _solve(ieq, symbol, **flags)[0]
                    rv = {inv.subs(isym, xi) for xi in _solve(eq, isym, **flags)}
                else:
                    rv = set(_solve(eq, symbol, **flags))
                if rv is not None:
                    result = list(ordered(rv))
                    # if the flag wasn't set then unset it since unrad results
                    # can be quite long or of very high order
                    flags['simplify'] = flags.get('simplify', False)
            else:
                pass  # for coverage

    # try _tsolve
    if result is False:
        flags.pop('tsolve', None)  # allow tsolve to be used on next pass
        try:
            soln = _tsolve(f_num, symbol, **flags)
            if soln is not None:
                result = soln
        except PolynomialError:
            pass
    # ----------- end of fallback ----------------------------

    if result is False:
        raise NotImplementedError('\n'.join([msg, not_impl_msg % f]))

    if flags.get('simplify', True):
        result = list(map(simplify, result))
        # we just simplified the solution so we now set the flag to
        # False so the simplification doesn't happen again in checksol()
        flags['simplify'] = False

    if checkdens:
        # reject any result that makes any denom. affirmatively 0;
        # if in doubt, keep it
        dens = denoms(f, symbols)
        result = [s for s in result if
                  all(not checksol(d, {symbol: s}, **flags)
                      for d in dens)]
    if check:
        # keep only results if the check is not False
        result = [r for r in result if
                  checksol(f_num, {symbol: r}, **flags) is not False]
    return result


def _solve_system(exprs, symbols, **flags):
    polys = []
    dens = set()
    failed = []
    result = False
    linear = False
    checkdens = check = flags.get('check', True)

    for j, g in enumerate(exprs):
        dens.update(denoms(g, symbols))
        i, d = _invert(g, *symbols)
        g = d - i
        g = g.as_numer_denom()[0]

        poly = g.as_poly(*symbols, extension=True)

        if poly is not None:
            polys.append(poly)
        else:
            failed.append(g)

    if not polys:
        solved_syms = []
    else:
        if all(p.is_linear for p in polys):
            n, m = len(polys), len(symbols)
            matrix = zeros(n, m + 1)

            for i, poly in enumerate(polys):
                for monom, coeff in poly.terms():
                    try:
                        j = monom.index(1)
                        matrix[i, j] = coeff
                    except ValueError:
                        matrix[i, m] = -coeff

            # returns a dictionary {symbols: values} or None
            result = solve_linear_system(matrix, *symbols, **flags)
            if failed:
                if result:
                    solved_syms = list(result.keys())
                else:
                    solved_syms = []
            else:
                linear = True

        else:
            result = solve_poly_system(polys, *symbols)
            solved_syms = list(set().union(*[{k for k in r.keys()}
                                             for r in result]))

    if result:
        if type(result) is dict:
            result = [result]
    else:
        result = [{}]

    if failed:
        # For each failed equation, see if we can solve for one of the
        # remaining symbols from that equation. If so, we update the
        # solution set and continue with the next failed equation,
        # repeating until we are done or we get an equation that can't
        # be solved.
        def _ok_syms(e, sort=False):
            rv = (e.free_symbols - solved_syms) & legal
            if sort:
                rv = list(rv)
                rv.sort(key=default_sort_key)
            return rv

        solved_syms = set(solved_syms)  # set of symbols we have solved for
        legal = set(symbols)  # what we are interested in

        # sort so equation with the fewest potential symbols is first
        for eq in ordered(failed, lambda _: len(_ok_syms(_))):
            u = Dummy()  # used in solution checking
            newresult = []
            bad_results = []
            got_s = set()
            hit = False
            for r in result:
                # update eq with everything that is known so far
                eq2 = eq.subs(r)
                # if check is True then we see if it satisfies this
                # equation, otherwise we just accept it
                if check and r:
                    b = checksol(u, u, eq2, minimal=True)
                    if b is not None:
                        # this solution is sufficient to know whether
                        # it is valid or not so we either accept or
                        # reject it, then continue
                        if b:
                            newresult.append(r)
                        else:
                            bad_results.append(r)
                        continue
                # search for a symbol amongst those available that
                # can be solved for
                ok_syms = _ok_syms(eq2, sort=True)
                if not ok_syms:
                    if r:
                        newresult.append(r)
                    break  # skip as it's independent of desired symbols
                for s in ok_syms:
                    soln = _solve(eq2, s, **flags)
                    # put each solution in r and append the now-expanded
                    # result in the new result list; use copy since the
                    # solution for s in being added in-place
                    for sol in soln:
                        if got_s and any(ss in sol.free_symbols for ss in got_s):
                            # sol depends on previously solved symbols: discard it
                            continue
                        rnew = r.copy()
                        for k, v in r.items():
                            rnew[k] = v.subs(s, sol)
                        # and add this new solution
                        rnew[s] = sol
                        newresult.append(rnew)
                    hit = True
                    got_s.add(s)
                if not hit:  # pragma: no cover
                    raise NotImplementedError('could not solve %s' % eq2)
            else:
                result = newresult
                for b in bad_results:
                    if b in result:
                        result.remove(b)

    default_simplify = bool(failed)  # rely on system-solvers to simplify
    if flags.get('simplify', default_simplify):
        for r in result:
            for k in r:
                r[k] = simplify(r[k])
        flags['simplify'] = False  # don't need to do so in checksol now

    if checkdens:
        result = [r for r in result
                  if not any(checksol(d, r, **flags) for d in dens)]

    if check and not linear:
        result = [r for r in result
                  if not any(checksol(e, r, **flags) is False for e in exprs)]

    result = [r for r in result if r]
    if linear and result:
        result = result[0]
    return result


def solve_linear(lhs, rhs=0, symbols=[], exclude=[]):
    r""" Return a tuple derived from f = lhs - rhs that is either:

        (numerator, denominator) of ``f``
            If this comes back as (0, 1) it means
            that ``f`` is independent of the symbols in ``symbols``, e.g::

                y*cos(x)**2 + y*sin(x)**2 - y = y*(0) = 0
                cos(x)**2 + sin(x)**2 = 1

            If it comes back as (0, 0) there is no solution to the equation
            amongst the symbols given.

            If the numerator is not zero then the function is guaranteed
            to be dependent on a symbol in ``symbols``.

        or

        (symbol, solution) where symbol appears linearly in the numerator of
        ``f``, is in ``symbols`` (if given) and is not in ``exclude`` (if given).

        No simplification is done to ``f`` other than and mul=True expansion,
        so the solution will correspond strictly to a unique solution.

    Examples
    ========

    >>> from diofant.abc import x, y, z

    These are linear in x and 1/x:

    >>> solve_linear(x + y**2)
    (x, -y**2)
    >>> solve_linear(1/x - y**2)
    (x, y**(-2))

    When not linear in x or y then the numerator and denominator are returned.

    >>> solve_linear(x**2/y**2 - 3)
    (x**2 - 3*y**2, y**2)

    If the numerator is a symbol then (0, 0) is returned if the solution for
    that symbol would have set any denominator to 0:

    >>> solve_linear(1/(1/x - 2))
    (0, 0)
    >>> 1/(1/x)  # to Diofant, this looks like x ...
    x
    >>> solve_linear(1/(1/x))  # so a solution is given
    (x, 0)

    If x is allowed to cancel, then this appears linear, but this sort of
    cancellation is not done so the solution will always satisfy the original
    expression without causing a division by zero error.

    >>> solve_linear(x**2*(1/x - z**2/x))
    (x**2*(-z**2 + 1), x)

    You can give a list of what you prefer for x candidates:

    >>> solve_linear(x + y + z, symbols=[y])
    (y, -x - z)

    You can also indicate what variables you don't want to consider:

    >>> solve_linear(x + y + z, exclude=[x, z])
    (y, -x - z)

    If only x was excluded then a solution for y or z might be obtained.
    """
    if isinstance(lhs, Equality):
        if rhs:
            raise ValueError(filldedent('''
            If lhs is an Equality, rhs must be 0 but was %s''' % rhs))
        rhs = lhs.rhs
        lhs = lhs.lhs
    dens = None
    eq = lhs - rhs
    n, d = eq.as_numer_denom()
    if not n:
        return S.Zero, S.One

    free = n.free_symbols
    if not symbols:
        symbols = free
    else:
        bad = [s for s in symbols if not s.is_Symbol]
        if bad:
            if len(bad) == 1:
                bad = bad[0]
            if len(symbols) == 1:
                eg = 'solve(%s, %s)' % (eq, symbols[0])
            else:
                eg = 'solve(%s, *%s)' % (eq, list(symbols))
            raise ValueError(filldedent('''
                solve_linear only handles symbols, not %s. To isolate
                non-symbols use solve, e.g. >>> %s <<<.
                             ''' % (bad, eg)))
        symbols = free.intersection(symbols)
    symbols = symbols.difference(exclude)

    # derivatives are easy to do but tricky to analyze to see if they are going
    # to disallow a linear solution, so for simplicity we just evaluate the
    # ones that have the symbols of interest
    derivs = defaultdict(list)
    for der in n.atoms(Derivative):
        csym = der.free_symbols & symbols
        for c in csym:
            derivs[c].append(der)

    if symbols:
        all_zero = True
        for xi in symbols:
            # if there are derivatives in this var, calculate them now
            if type(derivs[xi]) is list:
                derivs[xi] = {der: der.doit() for der in derivs[xi]}
            nn = n.subs(derivs[xi])
            dn = nn.diff(xi)
            if dn:
                all_zero = False
                if dn is S.NaN:
                    break
                if xi not in dn.free_symbols:
                    vi = -(nn.subs(xi, 0))/dn
                    if dens is None:
                        dens = denoms(eq, symbols)
                    if not any(checksol(di, {xi: vi}, minimal=True) is True
                               for di in dens):
                        # simplify any trivial integral
                        irep = [(i, i.doit()) for i in vi.atoms(Integral) if
                                i.function.is_number]
                        # do a slight bit of simplification
                        vi = expand_mul(vi.subs(irep))
                        if not d.has(xi) or not (d/xi).has(xi):
                            return xi, vi

        if all_zero:
            return S.Zero, S.One
    if n.is_Symbol:  # there was no valid solution
        n = d = S.Zero
    return n, d  # should we cancel now?


def minsolve_linear_system(system, *symbols, **flags):
    r"""
    Find a particular solution to a linear system.

    In particular, try to find a solution with the minimal possible number
    of non-zero variables. This is a very computationally hard problem.
    If ``quick=True``, a heuristic is used. Otherwise a naive algorithm with
    exponential complexity is used.
    """
    quick = flags.get('quick', False)
    # Check if there are any non-zero solutions at all
    s0 = solve_linear_system(system, *symbols, **flags)
    if not s0 or all(v == 0 for v in s0.values()):
        return s0
    if quick:
        # We just solve the system and try to heuristically find a nice
        # solution.
        s = solve_linear_system(system, *symbols)

        def update(determined, solution):
            delete = []
            for k, v in solution.items():
                solution[k] = v.subs(determined)
                if not solution[k].free_symbols:
                    delete.append(k)
                    determined[k] = solution[k]
            for k in delete:
                del solution[k]

        determined = {}
        update(determined, s)
        while s:
            # NOTE sort by default_sort_key to get deterministic result
            k = max((k for k in s.values()),
                    key=lambda x: (len(x.free_symbols), default_sort_key(x)))
            x = max(k.free_symbols, key=default_sort_key)
            if len(k.free_symbols) != 1:
                determined[x] = Integer(0)
            else:
                val = solve(k)[0][x]
                if val == 0 and all(v.subs(x, val) == 0 for v in s.values()):
                    determined[x] = Integer(1)
                else:
                    determined[x] = val
            update(determined, s)
        return determined
    else:
        # We try to select n variables which we want to be non-zero.
        # All others will be assumed zero. We try to solve the modified system.
        # If there is a non-trivial solution, just set the free variables to
        # one. If we do this for increasing n, trying all combinations of
        # variables, we will find an optimal solution.
        # We speed up slightly by starting at one less than the number of
        # variables the quick method manages.
        from itertools import combinations
        from ..utilities.misc import debug
        N = len(symbols)
        bestsol = minsolve_linear_system(system, *symbols, quick=True)
        n0 = len([x for x in bestsol.values() if x != 0])
        for n in range(n0 - 1, 1, -1):
            debug('minsolve: %s' % n)
            thissol = None
            for nonzeros in combinations(list(range(N)), n):
                subm = Matrix([system.col(i).T for i in nonzeros] + [system.col(-1).T]).T
                s = solve_linear_system(subm, *[symbols[i] for i in nonzeros])
                if s and not all(v == 0 for v in s.values()):
                    subs = [(symbols[v], Integer(1)) for v in nonzeros]
                    for k, v in s.items():
                        s[k] = v.subs(subs)
                    for sym in symbols:
                        if sym not in s:
                            if symbols.index(sym) in nonzeros:
                                s[sym] = Integer(1)
                            else:
                                s[sym] = Integer(0)
                    thissol = s
                    break
            if thissol is None:
                break
            bestsol = thissol
        return bestsol


def solve_undetermined_coeffs(equ, coeffs, sym, **flags):
    """Solve equation of a type p(x; a_1, ..., a_k) == q(x) where both
    p, q are univariate polynomials and f depends on k parameters.

    The result of this functions is a dictionary with symbolic
    values of those parameters with respect to coefficients in q.

    This functions accepts both Equations class instances and ordinary
    Diofant expressions. Specification of parameters and variable is
    obligatory for efficiency and simplicity reason.

    Examples
    ========

    >>> from diofant import Eq, Rational
    >>> from diofant.abc import a, b, c, x

    >>> solve_undetermined_coeffs(Eq(2*a*x + a+b, x), [a, b], x)
    {a: 1/2, b: -1/2}

    >>> solve_undetermined_coeffs(Eq(a*c*x + a+b, x), [a, b], x)
    {a: 1/c, b: -1/c}
    """
    if isinstance(equ, Equality):
        # got equation, so move all the
        # terms to the left hand side
        equ = equ.lhs - equ.rhs

    equ = cancel(equ).as_numer_denom()[0]

    system = list(collect(equ.expand(), sym, evaluate=False).values())

    if not any(equ.has(sym) for equ in system):
        # consecutive powers in the input expressions have
        # been successfully collected, so solve remaining
        # system using Gaussian elimination algorithm
        sol = solve(system, *coeffs, **flags)
        if sol:
            sol = sol[0]
        return sol
    else:
        return  # no solutions


# these are functions that have multiple inverse values per period
multi_inverses = {
    sin: lambda x: (asin(x), S.Pi - asin(x)),
    cos: lambda x: (acos(x), 2*S.Pi - acos(x)),
}


def _tsolve(eq, sym, **flags):
    """
    Helper for _solve that solves a transcendental equation with respect
    to the given symbol. Various equations containing powers and logarithms,
    can be solved.

    There is currently no guarantee that all solutions will be returned or
    that a real solution will be favored over a complex one.

    Either a list of potential solutions will be returned or None will be
    returned (in the case that no method was known to get a solution
    for the equation). All other errors (like the inability to cast an
    expression as a Poly) are unhandled.

    Examples
    ========

    >>> from diofant import log
    >>> from diofant.solvers.solvers import _tsolve as tsolve
    >>> from diofant.abc import x

    >>> tsolve(3**(2*x + 5) - 4, x)
    [-5/2 + log(2)/log(3), (-5*log(3)/2 + log(2) + I*pi)/log(3)]

    >>> tsolve(log(x) + 2*x, x)
    [LambertW(2)/2]
    """
    from .bivariate import bivariate_type, _solve_lambert, _filtered_gens

    if 'tsolve_saw' not in flags:
        flags['tsolve_saw'] = []
    if eq in flags['tsolve_saw']:
        return
    else:
        flags['tsolve_saw'].append(eq)

    rhs, lhs = _invert(eq, sym)

    if lhs == sym:
        return [rhs]
    try:
        if lhs.is_Add:
            # it's time to try factoring; powdenest is used
            # to try get powers in standard form for better factoring
            f = factor(powdenest(lhs - rhs))
            if f.is_Mul:
                return _solve(f, sym, **flags)
            if rhs:
                f = logcombine(lhs, force=flags.get('force', True))
                if f.count(log) != lhs.count(log):
                    if f.func is log:
                        return _solve(f.args[0] - exp(rhs), sym, **flags)
                    return _tsolve(f - rhs, sym)

        elif lhs.is_Pow:
            if lhs.exp.is_Integer:
                if lhs - rhs != eq:
                    return _solve(lhs - rhs, sym, **flags)
            elif sym not in lhs.exp.free_symbols:
                return _solve(lhs.base - rhs**(1/lhs.exp), sym, **flags)
            elif not rhs and sym in lhs.exp.free_symbols:
                # f(x)**g(x) only has solutions where f(x) == 0 and g(x) != 0 at
                # the same place
                sol_base = _solve(lhs.base, sym, **flags)
                return list(ordered(set(sol_base) -
                                    set(_solve(lhs.exp, sym, **flags))))
            elif (rhs is not S.Zero and
                  lhs.base.is_positive and
                  lhs.exp.is_extended_real):
                return _solve(lhs.exp*log(lhs.base) - log(rhs), sym, **flags)
            elif lhs.base == 0 and rhs == 1:
                return _solve(lhs.exp, sym, **flags)

        elif lhs.is_Mul and rhs.is_positive:
            llhs = expand_log(log(lhs))
            if llhs.is_Add:
                return _solve(llhs - log(rhs), sym, **flags)

        elif lhs.is_Function and len(lhs.args) == 1 and lhs.func in multi_inverses:
            # sin(x) = 1/3 -> x - asin(1/3) & x - (pi - asin(1/3))
            soln = []
            for i in multi_inverses[lhs.func](rhs):
                soln.extend(_solve(lhs.args[0] - i, sym, **flags))
            return list(ordered(soln))

        rewrite = lhs.rewrite(exp)
        if rewrite != lhs:
            return _solve(rewrite - rhs, sym, **flags)
    except NotImplementedError:
        pass

    # maybe it is a lambert pattern
    if flags.pop('bivariate', True):
        # lambert forms may need some help being recognized, e.g. changing
        # 2**(3*x) + x**3*log(2)**3 + 3*x**2*log(2)**2 + 3*x*log(2) + 1
        # to 2**(3*x) + (x*log(2) + 1)**3
        g = _filtered_gens(eq.as_poly(), sym)
        up_or_log = set()
        for gi in g:
            if gi.is_Pow and gi.base is S.Exp1 or gi.func is log:
                up_or_log.add(gi)
            elif gi.is_Pow:
                gisimp = powdenest(expand_power_exp(gi))
                if gisimp.is_Pow and sym in gisimp.exp.free_symbols:
                    up_or_log.add(gi)
        eq_down = expand_log(expand_power_exp(eq)).subs(
            dict(zip(up_or_log, [0]*len(up_or_log))))
        eq = expand_power_exp(factor(eq_down, deep=True) + (eq - eq_down))
        rhs, lhs = _invert(eq, sym)
        if lhs.has(sym):
            try:
                poly = lhs.as_poly()
                g = _filtered_gens(poly, sym)
                return _solve_lambert(lhs - rhs, sym, g)
            except NotImplementedError:
                # maybe it's a convoluted function
                if len(g) == 2:
                    try:
                        gpu = bivariate_type(lhs - rhs, *g)
                        if gpu is None:
                            raise NotImplementedError
                        g, p, u = gpu
                        flags['bivariate'] = False
                        inversion = _tsolve(g - u, sym, **flags)
                        if inversion:
                            sol = _solve(p, u, **flags)
                            return list(ordered({i.subs(u, s)
                                                 for i in inversion for s in sol}))
                    except NotImplementedError:
                        pass
                else:
                    pass

    if flags.pop('force', True):
        flags['force'] = False
        pos, reps = posify(lhs - rhs)
        for u, s in reps.items():
            if s == sym:
                break
        else:
            u = sym
        if pos.has(u):
            try:
                soln = _solve(pos, u, **flags)
                return list(ordered([s.subs(reps) for s in soln]))
            except NotImplementedError:
                pass


def _invert(eq, *symbols, **kwargs):
    """Return tuple (i, d) where ``i`` is independent of ``symbols`` and ``d``
    contains symbols. ``i`` and ``d`` are obtained after recursively using
    algebraic inversion until an uninvertible ``d`` remains. If there are no
    free symbols then ``d`` will be zero. Some (but not necessarily all)
    solutions to the expression ``i - d`` will be related to the solutions of
    the original expression.

    Examples
    ========

    >>> from diofant.solvers.solvers import _invert as invert
    >>> from diofant import sqrt, cos
    >>> from diofant.abc import x, y
    >>> invert(x - 3)
    (3, x)
    >>> invert(3)
    (3, 0)
    >>> invert(2*cos(x) - 1)
    (1/2, cos(x))
    >>> invert(sqrt(x) - 3)
    (3, sqrt(x))
    >>> invert(sqrt(x) + y, x)
    (-y, sqrt(x))
    >>> invert(sqrt(x) + y, y)
    (-sqrt(x), y)
    >>> invert(sqrt(x) + y, x, y)
    (0, sqrt(x) + y)

    If there is more than one symbol in a power's base and the exponent
    is not an Integer, then the principal root will be used for the
    inversion:

    >>> invert(sqrt(x + y) - 2)
    (4, x + y)
    >>> invert(sqrt(x + y) - 2)
    (4, x + y)

    If the exponent is an integer, setting ``integer_power`` to True
    will force the principal root to be selected:

    >>> invert(x**2 - 4, integer_power=True)
    (2, x)

    """
    eq = sympify(eq)
    free = eq.free_symbols
    if not symbols:
        symbols = free
    if not free & set(symbols):
        return eq, S.Zero

    dointpow = bool(kwargs.get('integer_power', False))

    lhs = eq
    rhs = S.Zero
    while True:
        was = lhs
        while True:
            indep, dep = lhs.as_independent(*symbols)

            # dep + indep == rhs
            if lhs.is_Add:
                # this indicates we have done it all
                if indep is S.Zero:
                    break

                lhs = dep
                rhs -= indep

            # dep * indep == rhs
            else:
                # this indicates we have done it all
                if indep is S.One:
                    break

                lhs = dep
                rhs /= indep

        # collect like-terms in symbols
        if lhs.is_Add:
            terms = {}
            for a in lhs.args:
                i, d = a.as_independent(*symbols)
                terms.setdefault(d, []).append(i)
            if any(len(v) > 1 for v in terms.values()):
                args = []
                for d, i in terms.items():
                    if len(i) > 1:
                        args.append(Add(*i)*d)
                    else:
                        args.append(i[0]*d)
                lhs = Add(*args)

        # if it's a two-term Add with rhs = 0 and two powers we can get the
        # dependent terms together, e.g. 3*f(x) + 2*g(x) -> f(x)/g(x) = -2/3
        if lhs.is_Add and not rhs and len(lhs.args) == 2 and \
                not lhs.is_polynomial(*symbols):
            a, b = ordered(lhs.args)
            ai, ad = a.as_independent(*symbols)
            bi, bd = b.as_independent(*symbols)
            if any(i.is_Pow for i in (ad, bd)):
                a_base, a_exp = ad.as_base_exp()
                b_base, b_exp = bd.as_base_exp()
                if a_base == b_base:
                    # a = -b
                    lhs = powsimp(powdenest(ad/bd))
                    rhs = -bi/ai
                else:
                    rat = ad/bd
                    _lhs = powsimp(ad/bd)
                    if _lhs != rat:
                        lhs = _lhs
                        rhs = -bi/ai
            if ai*bi is S.NegativeOne:
                if all(
                        isinstance(i, Function) for i in (ad, bd)) and \
                        ad.func == bd.func and len(ad.args) == len(bd.args):
                    if len(ad.args) == 1:
                        lhs = ad.args[0] - bd.args[0]
                    else:  # pragma: no cover
                        # should be able to solve
                        # f(x, y) == f(2, 3) -> x == 2
                        # f(x, x + y) == f(2, 3) -> x == 2 or x == 3 - y
                        raise NotImplementedError('equal function with more than 1 argument')

        elif lhs.is_Mul and any(a.is_Pow for a in lhs.args):
            lhs = powsimp(powdenest(lhs))

        if lhs.is_Function:
            if hasattr(lhs, 'inverse') and len(lhs.args) == 1:
                #                    -1
                # f(x) = g  ->  x = f  (g)
                #
                # /!\ inverse should not be defined if there are multiple values
                # for the function -- these are handled in _tsolve
                #
                rhs = lhs.inverse()(rhs)
                lhs = lhs.args[0]
            elif lhs.func is atan2:
                y, x = lhs.args
                lhs = 2*atan(y/(sqrt(x**2 + y**2) + x))

        if lhs.is_Pow and lhs.base is S.Exp1:
            rhs = log(rhs)
            lhs = lhs.exp

        if rhs and lhs.is_Pow and lhs.exp.is_Integer and lhs.exp < 0:
            lhs = 1/lhs
            rhs = 1/rhs

        # base**a = b -> base = b**(1/a) if
        #    a is an Integer and dointpow=True (this gives real branch of root)
        #    a is not an Integer and the equation is multivariate and the
        #      base has more than 1 symbol in it
        # The rationale for this is that right now the multi-system solvers
        # doesn't try to resolve generators to see, for example, if the whole
        # system is written in terms of sqrt(x + y) so it will just fail, so we
        # do that step here.
        if lhs.is_Pow and (
            lhs.exp.is_Integer and dointpow or not lhs.exp.is_Integer and
                len(symbols) > 1 and len(lhs.base.free_symbols & set(symbols)) > 1):
            rhs = rhs**(1/lhs.exp)
            lhs = lhs.base

        if lhs == was:
            break
    return rhs, lhs
