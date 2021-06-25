"""
This module contain solvers for all kinds of equations,
algebraic or transcendental.
"""

import warnings
from collections import defaultdict
from types import GeneratorType

from ..core import (Add, Dummy, E, Equality, Expr, Float, Function, Ge, I,
                    Integer, Lambda, Mul, Symbol, expand_log, expand_mul,
                    expand_power_exp, nan, nfloat, pi, preorder_traversal)
from ..core.assumptions import check_assumptions
from ..core.compatibility import is_sequence, iterable
from ..core.function import AppliedUndef
from ..core.logic import fuzzy_and
from ..core.relational import Relational
from ..core.sympify import sympify
from ..functions import (Abs, Max, Min, Piecewise, acos, arg, asin, atan,
                         atan2, cos, exp, im, log, piecewise_fold, re, sin,
                         sqrt, tan)
from ..functions.elementary.hyperbolic import HyperbolicFunction
from ..functions.elementary.trigonometric import TrigonometricFunction
from ..logic import false, true
from ..matrices import Matrix, zeros
from ..polys import Poly, RootOf, factor, roots
from ..polys.polyerrors import PolynomialError
from ..simplify.fu import TR1
from ..simplify.powsimp import powdenest, powsimp
from ..simplify.radsimp import denom
from ..simplify.simplify import logcombine, nsimplify, posify, simplify
from ..simplify.sqrtdenest import unrad
from ..utilities import default_sort_key, filldedent, ordered
from ..utilities.iterables import uniq
from .polysys import solve_linear_system, solve_poly_system, solve_surd_system
from .utils import checksol


__all__ = 'solve', 'solve_linear', 'minsolve_linear_system'


def denoms(eq, symbols=None):
    """Return (recursively) set of all denominators that appear in eq
    that contain any symbol in iterable ``symbols``; if ``symbols`` is
    None (default) then all denominators will be returned.

    Examples
    ========

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
        if den == 1:
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

        check : bool, optional
            If False, don't do any testing of solutions.  Default is
            True, i.e. the solutions are checked and those that doesn't
            satisfy given assumptions on symbols solved for or make any
            denominator zero - are automatically excluded.
        warn : bool, optional
            Show a warning if :func:`~diofant.solvers.utils.checksol`
            could not conclude.  Default is False.
        simplify : bool, optional
            Enable simplification (default) for all but polynomials of
            order 3 or greater before returning them and (if check is
            not False) use the general simplify function on the solutions
            and the expression obtained when they are substituted into the
            function which should be zero.
        rational : bool or None, optional
            If True, recast Floats as Rational.  If None (default),
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

    Single equation:

    >>> solve(x**2 - y**2)
    [{x: -y}, {x: y}]
    >>> solve(x**2 - 1)
    [{x: -1}, {x: 1}]

    We could restrict solutions by using assumptions:

    >>> p = Symbol('p', positive=True)
    >>> solve(p**2 - 1)
    [{p: 1}]

    Several equations:

    >>> solve((x + 5*y - 2, -3*x + 6*y - 15))
    [{x: -3, y: 1}]
    >>> solve((x + 5*y - 2, -3*x + 6*y - z))
    [{x: -5*z/21 + 4/7, y: z/21 + 2/7}]

    No solution:

    >>> solve([x + 3, x - 3])
    []

    Notes
    =====

    When an object other than a Symbol is given as a symbol, it is
    isolated algebraically and an implicit solution may be obtained.
    This is mostly provided as a convenience to save one from replacing
    the object with a Symbol and solving for that Symbol. It will only
    work if the specified object can be replaced with a Symbol using the
    subs method.

    >>> solve(f(x) - x, f(x))
    [{f(x): x}]
    >>> solve(f(x).diff(x) - f(x) - x, f(x).diff(x))
    [{Derivative(f(x), x): x + f(x)}]

    See Also
    ========

    diofant.solvers.recurr.rsolve : solving recurrence equations
    diofant.solvers.ode.dsolve : solving differential equations
    diofant.solvers.inequalities.reduce_inequalities : solving inequalities

    """
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
            if 'ImmutableMatrix' in (type(a).__name__ for a in fi.args):
                f[i] = fi.lhs - fi.rhs
            else:
                f[i] = Add(fi.lhs, -fi.rhs, evaluate=False)
        elif isinstance(fi, Relational):
            raise ValueError(f'Only expressions or equalities supported, got {fi}')
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
            f[i] = Integer(0)

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
                raise NotImplementedError(f'solving {a} when the argument '
                                          'is not real or imaginary.')
            reps.append((a, piece(a.args[0]) if a.args[0].is_extended_real else
                         piece(a.args[0]*I)))
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
            irf.append((s, re(s) + I*im(s)))
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
            s_new = Dummy(f'X{i:d}')
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
                pot.skip()
            elif (isinstance(p, bool) or not p.args or p in symset or
                  p.is_Add or p.is_Mul or p.is_Pow or p.is_Function or
                  isinstance(p, RootOf)) and p.func not in (re, im):
                pass
            elif p not in seen:
                seen.add(p)
                if p.free_symbols & symset:
                    non_inverts.add(p)
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

    # piecewise_fold might cancel denominators, so be sure to check them.
    piecewise_dens = set()

    # Any embedded piecewise functions need to be brought out to the
    # top level so that the appropriate strategy gets selected.
    # However, this is necessary only if one of the piecewise
    # functions depends on one of the symbols we are solving for.
    for i, fi in enumerate(f):
        if any(e.has(*symbols) for e in fi.atoms(Piecewise)):
            piecewise_dens |= denoms(fi, symbols)
            f[i] = piecewise_fold(fi)

    if all(_ == 0 for _ in f):
        return [{}]

    #
    # try to get a solution
    ###########################################################################
    if bare_f and len(symbols) == 1:
        solution = [{symbols[0]: s} for s in _solve(f[0], symbols[0], **flags)]
    else:
        solution = _solve_system(f, symbols, **flags)

    #
    # postprocessing
    ###########################################################################
    # Restore masked-off objects
    if non_inverts:
        solution = [{k: v.subs(non_inverts) for k, v in s.items()}
                    for s in solution]

    # Restore original "symbols" if a dictionary is returned.
    # This is not necessary for
    #   - the single univariate equation case
    #     since the symbol will have been removed from the solution;
    #
    # ** unless there were Derivatives with the symbols, but those were handled
    #    above.
    if symbol_swapped:
        symbols = [swap_sym[k] for k in symbols]
        if solution:
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
        def test_assumptions(sol):
            return fuzzy_and([check_assumptions(sol[sym], **sym._assumptions)
                              for sym in sol])

        solution = [s for s in solution if test_assumptions(s) is not False]

        warn = flags.get('warn', False)
        got_None = [s for s in solution if not test_assumptions(s)]
        if warn and got_None:
            warnings.warn(filldedent("""
                \tWarning: assumptions concerning following solution(s)
                can't be checked:""" + '\n\t' +
                                     ', '.join(str(s) for s in got_None)))

        solution = [s for s in solution if
                    all(not checksol(den, s, **flags) for den in piecewise_dens)]

    #
    # done
    ###########################################################################

    # Make sure that a list of solutions is ordered in a canonical way.
    solution.sort(key=default_sort_key)

    return solution


def _solve(f, symbol, **flags):
    """Return a checked solution for f in terms of one or more of the
    symbols. A list (possibly empty) should be returned.

    If no method is implemented to solve the equation, a NotImplementedError
    will be raised. In the case that conversion of an expression to a Poly
    gives None a ValueError will be raised.

    """
    not_impl_msg = 'No algorithms are implemented to solve equation %s'

    # /!\ capture this flag then set it to False so that no checking in
    # recursive calls will be done; only the final answer is checked
    flags['check'] = checkdens = check = flags.pop('check', True)

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
            dens = denoms(f, [symbol])
            result = [s for s in result if
                      all(not checksol(den, {symbol: s}, **flags) for den in
                          dens)]
        # set flags for quick exit at end; solutions for each
        # factor were already checked and simplified
        check = False
        flags['simplify'] = False

    elif f.is_Piecewise:
        result = set()
        for n, (expr, cond) in enumerate(f.args):
            candidates = _solve(piecewise_fold(expr), symbol, **flags)
            for candidate in candidates:
                if candidate in result:
                    continue
                try:
                    v = (cond == true) or cond.subs({symbol: candidate})
                except TypeError:
                    v = False
                if v != false:
                    # Only include solutions that do not match the condition
                    # of any previous pieces.
                    matches_other_piece = False
                    for other_n, (other_expr, other_cond) in enumerate(f.args):  # pragma: no branch
                        if other_n == n:
                            break
                        try:
                            if other_cond.subs({symbol: candidate}) == true:
                                matches_other_piece = True
                                break
                        except TypeError:
                            pass
                    if not matches_other_piece:
                        v = v == true or v.doit()
                        if isinstance(v, Relational):
                            v = v.canonical
                        result.add(Piecewise(
                            (candidate, v),
                            (nan, True)
                        ))
        check = False
        flags['simplify'] = False
    else:
        # first see if it really depends on symbol and whether there
        # is a linear solution
        f_num, sol = solve_linear(f, symbol)
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

        poly = Poly(f_num)
        gens = [g for g in poly.gens if g.has(symbol)]

        def _as_base_q(x):
            """Return (b**e, q) for x = b**(p*e/q) where p/q is the leading
            Rational of the exponent of x, e.g. exp(-2*x/3) -> (exp(x), 3)

            """
            b, e = x.as_base_exp()
            if e.is_Rational:
                return b, e.denominator
            if not e.is_Mul:
                return x, 1
            c, ee = e.as_coeff_Mul()
            if c.is_Rational and c != 1:  # c could be a Float
                return b**ee, c.denominator
            return x, 1

        if len(gens) > 1:
            # If there is more than one generator, it could be that the
            # generators have the same base but different powers, e.g.
            #   >>> (exp(x) + 1/exp(x)).as_poly()
            #   Poly(exp(-x) + exp(x), exp(-x), exp(x), domain='ZZ')
            #
            # If unrad was not disabled then there should be no rational
            # exponents appearing as in
            #   >>> (sqrt(x) + sqrt(sqrt(x))).as_poly()
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
                        # don't check the rewritten form --check
                        # solutions in the un-rewritten form below
                        flags['check'] = False
                        result = _solve(newf, symbol, **flags)
                        flags['check'] = check

                # just a simple case - see if replacement of single function
                # clears all symbol-dependent functions, e.g.
                # log(x) - log(log(x) - 1) - 3 can be solved even though it has
                # two generators.

                if result is False and funcs:
                    funcs = list(ordered(funcs))  # put shallowest function first
                    f1 = funcs[0]
                    t = Dummy('t')
                    # perform the substitution
                    ftry = f_num.subs({f1: t})

                    # if no Functions left, we can proceed with usual solve
                    if not ftry.has(symbol):
                        cv_sols = _solve(ftry, t, **flags)
                        cv_inv = _solve(t - f1, symbol, **flags)[0]
                        sols = []
                        for sol in cv_sols:
                            sols.append(cv_inv.subs({t: sol}))
                        result = list(ordered(sols))

                if result is False:
                    msg = f'multiple generators {gens}'

            else:
                # e.g. case where gens are exp(x), exp(-x)
                u = bases.pop()
                t = Dummy('t')
                inv = _solve(u - t, symbol, **flags)
                # this will be resolved by factor in _tsolve but we might
                # as well try a simple expansion here to get things in
                # order so something like the following will work now without
                # having to factor:
                #
                # >>> eq = (exp(I*(-x-2))+exp(I*(x+2)))
                # >>> eq.subs({exp(x): y})  # fails
                # exp(I*(-x - 2)) + exp(I*(x + 2))
                # >>> eq.expand().subs({exp(x): y})  # works
                # y**I*exp(2*I) + y**(-I)*exp(-2*I)

                def _expand(p):
                    b, e = p.as_base_exp()
                    e = expand_mul(e)
                    return expand_power_exp(b**e)
                ftry = f_num.replace(lambda w: w.is_Pow, _expand).subs({u: t})
                assert not ftry.has(symbol)
                soln = _solve(ftry, t, **flags)
                sols = []
                for sol in soln:
                    for i in inv:
                        sols.append(i.subs({t: sol}))
                result = list(ordered(sols))

        else:
            # There is only one generator that we are interested in, but
            # there may have been more than one generator identified by
            # polys (e.g. for symbols other than the one we are interested
            # in) so recast the poly in terms of our generator of interest.

            poly = Poly(f_num, gens[0], extension=False)

            # if we aren't on the tsolve-pass, use roots
            if not flags.pop('tsolve', False):
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
                    soln = poly.all_roots()
                else:
                    soln = list(soln)

                u = poly.gen
                if u != symbol:
                    try:
                        t = Dummy('t')
                        iv = _solve(u - t, symbol, **flags)
                        soln = list(ordered({i.subs({t: s}) for i in iv for s in soln}))
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
        u = unrad(f_num, symbol)
        if u:
            eq, cov = u
            if cov:
                isym, ieq = cov
                inv = _solve(ieq, symbol, **flags)[0]
                rv = {inv.subs({isym: xi}) for xi in _solve(eq, isym, **flags)}
            else:
                rv = set(_solve(eq, symbol, **flags))
            result = list(ordered(rv))
            # if the flag wasn't set then unset it since unrad results
            # can be quite long or of very high order
            flags['simplify'] = flags.get('simplify', False)

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
        dens = denoms(f, [symbol])
        result = [s for s in result if
                  all(not checksol(d, {symbol: s}, **flags)
                      for d in dens)]
    if check:
        # keep only results if the check is not False
        result = [r for r in result if
                  checksol(f_num, {symbol: r}, **flags) is not False]
    return result


def _solve_system(exprs, symbols, **flags):
    """Return a checked solution for list of exprs in terms of one or more
    of the symbols. A list of dict's (possibly empty) should be returned.

    """
    if len(symbols) != 1 and len(exprs) == 1:
        f = exprs[0]
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

    polys = []
    surds = []
    dens = set()
    failed = []
    result = [{}]
    solved_syms = []
    algebraic = False
    inversions = False
    checkdens = check = flags.get('check', True)

    for j, g in enumerate(exprs):
        dens.update(denoms(g, symbols))
        i, d = _invert(g, *symbols)
        g = d - i
        if exprs[j] not in (+g, -g):
            inversions = True
        g = g.as_numer_denom()[0]

        poly = g.as_poly(*symbols)

        if poly is not None:
            polys.append(poly)
        elif g.is_algebraic_expr(*symbols):
            surds.append(g)
        else:
            failed.append(g)

    if surds:
        result = solve_surd_system([_.as_expr() for _ in polys] +
                                   surds, *symbols)
        solved_syms = list(set().union(*[set(r) for r in result]))
    elif polys and all(p.is_linear for p in polys):
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
        solved_syms = list(result) if result else []
        result = [result] if result else [{}]
    elif polys:
        result = solve_poly_system(polys, *symbols)
        solved_syms = list(set().union(*[set(r) for r in result]))

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
                    b = checksol(u, {u: eq2}, minimal=True)
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
                    newresult.append(r)
                    break  # skip as it's independent of desired symbols
                for s in ok_syms:
                    soln = _solve(eq2, s, **flags)
                    # put each solution in r and append the now-expanded
                    # result in the new result list; use copy since the
                    # solution for s in being added in-place
                    for sol in soln:
                        if got_s and got_s & sol.free_symbols:
                            # sol depends on previously solved symbols: discard it
                            continue
                        rnew = r.copy()
                        for k, v in r.items():
                            rnew[k] = v.subs({s: sol})
                        # and add this new solution
                        rnew[s] = sol
                        newresult.append(rnew)
                    hit = True
                    got_s.add(s)
                if not hit:
                    raise NotImplementedError(f'could not solve {eq2}')
            else:
                result = newresult
                assert not any(b in bad_results for b in result)
    else:
        algebraic = True

    default_simplify = bool(failed)  # rely on system-solvers to simplify
    if flags.get('simplify', default_simplify):
        for r in result:
            for k in r:
                r[k] = simplify(r[k])
        flags['simplify'] = False  # don't need to do so in checksol now

    if checkdens:
        result = [r for r in result
                  if not any(checksol(d, r, **flags) for d in dens)]

    if check and (inversions or not algebraic):
        result = [r for r in result
                  if not any(checksol(e, r, **flags) is False for e in exprs)]

    return [r for r in result if r]


def solve_linear(f, x):
    r"""
    Solve equation ``f`` wrt variable ``x``.

    Returns
    =======

    tuple
        ``(x, solution)``, if there is a linear solution, ``(0, 1)`` if
        ``f`` is independent of the symbol ``x``, ``(0, 0)`` if solution set
        any denominator of ``f`` to zero or ``(numerator, denominator)``
        of ``f``, if it's a nonlinear expression wrt ``x``.

    Examples
    ========

    >>> solve_linear(1/x - y**2, x)
    (x, y**(-2))
    >>> solve_linear(x**2/y**2 - 3, x)
    (x**2 - 3*y**2, y**2)
    >>> solve_linear(y, x)
    (0, 1)
    >>> solve_linear(1/(1/x - 2), x)
    (0, 0)

    """
    if not x.is_Symbol:
        raise ValueError(f'{x} is not a Symbol')
    f = f.replace(lambda e: e.is_Derivative, lambda e: e.doit())
    n, d = res = f.as_numer_denom()
    poly = n.as_poly(x, extension=False)
    if poly is not None and poly.is_linear:
        a, b = n.expand().coeff(x, 1), n.expand().coeff(x, 0)
        if a != 0 and d.subs({x: -b/a}) != 0:
            res = (x, -b/a)
    if not n.simplify().has(x):
        res = Integer(0), Integer(1)
    if x == res[0] and any(checksol(_, {x: res[1]}) for _ in denoms(f, [x])):
        res = Integer(0), Integer(0)
    return res


def minsolve_linear_system(system, *symbols, **flags):
    r"""Find a particular solution to a linear system.

    In particular, try to find a solution with the minimal possible number
    of non-zero variables. This is a very computationally hard problem.

    Parameters
    ==========

    system : Matrix
        Nx(M+1) matrix, which means it has to be in augmented form.
    \*symbols : list
        List of M Symbol’s.
    \*\*flags : dict
        A dictionary of following parameters:

        quick : boolean, optional
            If True, a heuristic is used.  Otherwise (default) a naive
            algorithm with exponential complexity is used.

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
                if val == 0 and all(v.subs({x: val}) == 0 for v in s.values()):
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

        N = len(symbols)
        bestsol = minsolve_linear_system(system, *symbols, quick=True)
        n0 = len([x for x in bestsol.values() if x != 0])
        for n in range(n0 - 1, 1, -1):
            thissol = None
            for nonzeros in combinations(list(range(N)), n):
                subm = Matrix([system[:, i].T for i in nonzeros] + [system[:, -1].T]).T
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


# these are functions that have multiple inverse values per period
multi_inverses = {
    sin: lambda x: (asin(x), pi - asin(x)),
    cos: lambda x: (acos(x), 2*pi - acos(x)),
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

    >>> _tsolve(3**(2*x + 5) - 4, x)
    [-5/2 + log(2)/log(3), (-5*log(3)/2 + log(2) + I*pi)/log(3)]

    >>> _tsolve(log(x) + 2*x, x)
    [LambertW(2)/2]

    """
    from .bivariate import _filtered_gens, _solve_lambert, bivariate_type

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
                    if isinstance(f, log):
                        return _solve(f.args[0] - exp(rhs), sym, **flags)
                    else:
                        raise NotImplementedError

        elif lhs.is_Pow:
            if lhs.exp.is_Integer and lhs - rhs != eq:
                return _solve(lhs - rhs, sym, **flags)
            elif sym not in lhs.exp.free_symbols:
                return _solve(lhs.base - rhs**(1/lhs.exp), sym, **flags)
            elif not rhs and sym in lhs.exp.free_symbols:
                # f(x)**g(x) only has solutions where f(x) == 0 and g(x) != 0 at
                # the same place
                sol_base = _solve(lhs.base, sym, **flags)
                return list(ordered(set(sol_base) -
                                    set(_solve(lhs.exp, sym, **flags))))
            elif (rhs != 0 and
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
            if gi.is_Pow and gi.base is E or isinstance(gi, log):
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
                            return list(ordered({i.subs({u: s})
                                                 for i in inversion for s in sol}))
                        else:
                            raise NotImplementedError
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

    >>> _invert(x - 3)
    (3, x)
    >>> _invert(3)
    (3, 0)
    >>> _invert(2*cos(x) - 1)
    (1/2, cos(x))
    >>> _invert(sqrt(x) - 3)
    (3, sqrt(x))
    >>> _invert(sqrt(x) + y, x)
    (-y, sqrt(x))
    >>> _invert(sqrt(x) + y, y)
    (-sqrt(x), y)
    >>> _invert(sqrt(x) + y, x, y)
    (0, sqrt(x) + y)

    If there is more than one symbol in a power's base and the exponent
    is not an Integer, then the principal root will be used for the
    inversion:

    >>> _invert(sqrt(x + y) - 2)
    (4, x + y)
    >>> _invert(sqrt(x + y) - 2)
    (4, x + y)

    If the exponent is an integer, setting ``integer_power`` to True
    will force the principal root to be selected:

    >>> _invert(x**2 - 4, integer_power=True)
    (2, x)

    """
    eq = sympify(eq)
    free = eq.free_symbols
    if not symbols:
        symbols = free
    if not free & set(symbols):
        return eq, Integer(0)

    dointpow = bool(kwargs.get('integer_power', False))

    lhs = eq
    rhs = Integer(0)
    while True:
        was = lhs
        while True:
            indep, dep = lhs.as_independent(*symbols)

            # dep + indep == rhs
            if lhs.is_Add:
                # this indicates we have done it all
                if indep == 0:
                    break

                lhs = dep
                rhs -= indep

            # dep * indep == rhs
            else:
                # this indicates we have done it all
                if indep == 1:
                    break

                lhs = dep
                rhs /= indep

        # collect like-terms in symbols
        if lhs.is_Add:
            terms = defaultdict(list)
            for a in lhs.args:
                i, d = a.as_independent(*symbols)
                terms[d].append(i)
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
            if ai*bi == -1:
                if all(
                        isinstance(i, Function) for i in (ad, bd)) and \
                        ad.func == bd.func and len(ad.args) == len(bd.args):
                    if len(ad.args) == 1:
                        lhs = ad.args[0] - bd.args[0]
                    else:
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
            elif isinstance(lhs, atan2):
                y, x = lhs.args
                lhs = 2*atan(y/(sqrt(x**2 + y**2) + x))

        if lhs.is_Pow and lhs.base is E:
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
