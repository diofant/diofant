from ...core import (Basic, Dummy, Equality, Expr, Function, Integer, Tuple,
                     diff, nan, oo)
from ...core.relational import Relational
from ...logic import And, Not, Or, false, to_cnf, true
from ...logic.boolalg import Boolean
from ...sets import ExtendedReals
from ...utilities import default_sort_key
from .miscellaneous import Max, Min


class ExprCondPair(Tuple):
    """Represents an expression, condition pair."""

    def __new__(cls, expr, cond):
        return super().__new__(cls, expr, cond)

    @property
    def expr(self):
        """Returns the expression of this pair."""
        return self.args[0]

    @property
    def cond(self):
        """Returns the condition of this pair."""
        return self.args[1]

    @property
    def is_commutative(self):
        return self.expr.is_commutative


class Piecewise(Function):
    """
    Represents a piecewise function.

    Usage:

      Piecewise( (expr,cond), (expr,cond), ... )
        - Each argument is a 2-tuple defining an expression and condition
        - The conds are evaluated in turn returning the first that is True.
          If any of the evaluated conds are not determined explicitly False,
          e.g. x < 1, the function is returned in symbolic form.
        - Pairs where the cond is explicitly False, will be removed.
        - The last condition is explicitely True.  It's value is a default
          value for the function (nan if not specified othewise).

    Examples
    ========

      >>> f = x**2
      >>> g = log(x)
      >>> p = Piecewise((0, x < -1), (f, x <= 1), (g, True))
      >>> p.subs({x: 1})
      1
      >>> p.subs({x: 5})
      log(5)

    See Also
    ========

    diofant.functions.elementary.piecewise.piecewise_fold

    """

    nargs = None
    is_Piecewise = True

    def __new__(cls, *args, **options):
        newargs = []
        for ec in args:
            # ec could be a ExprCondPair or a tuple
            pair = ExprCondPair(*getattr(ec, 'args', ec))
            cond = pair.cond
            if cond == false:
                continue
            if not isinstance(cond, (bool, Relational, Boolean)):
                raise TypeError(
                    f'Cond {cond} is of type {type(cond)}, but must be a Relational,'
                    ' Boolean, or a built-in bool.')
            newargs.append(pair)
            if cond == true:
                break

        if cond != true:
            newargs.append(ExprCondPair(nan, true))

        if options.pop('evaluate', True):
            r = cls.eval(*newargs)
        else:
            r = None

        if r is None:
            return Expr.__new__(cls, *newargs, **options)
        return r

    @classmethod
    def eval(cls, *args):
        # Check for situations where we can evaluate the Piecewise object.
        # 1) Hit an unevaluable cond (e.g. x<1) -> keep object
        # 2) Hit a true condition -> return that expr
        # 3) Remove false conditions, if no conditions left -> raise ValueError
        all_conds_evaled = True    # Do all conds eval to a bool?
        piecewise_again = False    # Should we pass args to Piecewise again?
        non_false_ecpairs = []
        or1 = Or(*[cond for (_, cond) in args if cond != true])
        for expr, cond in args:
            # Check here if expr is a Piecewise and collapse if one of
            # the conds in expr matches cond. This allows the collapsing
            # of Piecewise((Piecewise(x,x<0),x<0)) to Piecewise((x,x<0)).
            # This is important when using piecewise_fold to simplify
            # multiple Piecewise instances having the same conds.
            # Eventually, this code should be able to collapse Piecewise's
            # having different intervals, but this will probably require
            # using the new assumptions.
            if isinstance(expr, Piecewise):
                or2 = Or(*[c for (_, c) in expr.args if c != true])
                for e, c in expr.args:
                    # Don't collapse if cond is "True" as this leads to
                    # incorrect simplifications with nested Piecewises.
                    if c == cond and (or1 == or2 or cond != true):
                        expr = e
                        piecewise_again = True
            cond_eval = cls.__eval_cond(cond)
            if cond_eval is None:
                all_conds_evaled = False
            elif cond_eval:
                if all_conds_evaled:
                    return expr
            if len(non_false_ecpairs) != 0:
                if non_false_ecpairs[-1].cond == cond:
                    continue
                if non_false_ecpairs[-1].expr == expr:
                    newcond = Or(cond, non_false_ecpairs[-1].cond)
                    if isinstance(newcond, (And, Or)):
                        newcond = to_cnf(newcond)
                    non_false_ecpairs[-1] = ExprCondPair(expr, newcond)
                    continue
            non_false_ecpairs.append(ExprCondPair(expr, cond))
        if len(non_false_ecpairs) != len(args) or piecewise_again:
            return cls(*non_false_ecpairs)

    def doit(self, **hints):
        """Evaluate this piecewise function."""
        newargs = []
        for e, c in self.args:
            if hints.get('deep', True):
                e = e.doit(**hints)
                c = c.doit(**hints)
            newargs.append((e, c))
        return self.func(*newargs)

    def _eval_as_leading_term(self, x):
        for e, c in self.args:  # pragma: no branch
            if c == true or c.subs({x: 0}) == true:
                return e.as_leading_term(x)

    def _eval_adjoint(self):
        return self.func(*[(e.adjoint(), c) for e, c in self.args])

    def _eval_conjugate(self):
        return self.func(*[(e.conjugate(), c) for e, c in self.args])

    def _eval_derivative(self, s):
        return self.func(*[(diff(e, s), c) for e, c in self.args])

    def _eval_evalf(self, prec):
        return self.func(*[(e.evalf(prec), c) for e, c in self.args])

    def _eval_integral(self, x):
        from ...integrals import integrate
        return self.func(*[(integrate(e, x), c) for e, c in self.args])

    def _eval_interval(self, x, a, b):
        """Evaluates the function along the sym in a given interval ab."""
        # FIXME: Currently complex intervals are not supported.  A possible
        # replacement algorithm, discussed in issue sympy/sympy#5227, can be found in the
        # following papers;
        #     http://portal.acm.org/citation.cfm?id=281649
        #     http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.70.4127&rep=rep1&type=pdf
        from .complexes import Abs

        if x.is_real is None:
            d = Dummy('d', real=True)
            return self.subs({x: d})._eval_interval(d, a, b)

        if self.has(Abs):
            return piecewise_fold(self.rewrite(Abs, Piecewise))._eval_interval(x, a, b)

        if a is None or b is None:
            # In this case, it is just simple substitution
            return piecewise_fold(super()._eval_interval(x, a, b))

        mul = 1
        if a == b:
            return Integer(0)
        if (a > b) == true:
            a, b, mul = b, a, -1
        elif (a <= b) != true:
            newargs = []
            for e, c in self.args:
                intervals = self._sort_expr_cond(
                    x, -oo, oo, c)
                values = []
                for lower, upper, expr in intervals:
                    if (a < lower) == true:
                        mid = lower
                        rep = b
                        val = e._eval_interval(x, mid, b)
                        val += self._eval_interval(x, a, mid)
                    elif (a > upper) == true:
                        mid = upper
                        rep = b
                        val = e._eval_interval(x, mid, b)
                        val += self._eval_interval(x, a, mid)
                    elif (a >= lower) == true and (a <= upper) == true:
                        rep = b
                        val = e._eval_interval(x, a, b)
                    elif (b < lower) == true:
                        mid = lower
                        rep = a
                        val = e._eval_interval(x, a, mid)
                        val += self._eval_interval(x, mid, b)
                    elif (b > upper) == true:
                        mid = upper
                        rep = a
                        val = e._eval_interval(x, a, mid)
                        val += self._eval_interval(x, mid, b)
                    elif ((b >= lower) == true) and ((b <= upper) == true):
                        rep = a
                        val = e._eval_interval(x, a, b)
                    else:
                        raise NotImplementedError(
                            """The evaluation of a Piecewise interval when both the lower
                            and the upper limit are symbolic is not yet implemented.""")
                    values.append(val)
                if len(set(values)) == 1:
                    c = c.subs({x: rep})
                    e = values[0]
                    newargs.append((e, c))
                else:
                    for i, vi in enumerate(values):
                        newargs.append((vi, (c == true and i == len(values) - 1) or
                                        And(rep >= intervals[i][0], rep <= intervals[i][1])))
            return self.func(*newargs)

        # Determine what intervals the expr,cond pairs affect.
        int_expr = self._sort_expr_cond(x, a, b)

        # Finally run through the intervals and sum the evaluation.
        ret_fun = 0
        for int_a, int_b, expr in int_expr:
            if isinstance(expr, Piecewise):
                # If we still have a Piecewise by now, _sort_expr_cond would
                # already have determined that its conditions are independent
                # of the integration variable, thus we just use substitution.
                ret_fun += piecewise_fold(
                    super(Piecewise, expr)._eval_interval(x, Max(a, int_a), Min(b, int_b)))
            else:
                ret_fun += expr._eval_interval(x, Max(a, int_a), Min(b, int_b))
        return mul * ret_fun

    def _sort_expr_cond(self, sym, a, b, targetcond=None):
        """Determine what intervals the expr, cond pairs affect.

        1) If cond is True, then log it as default
        1.1) Currently if cond can't be evaluated, throw NotImplementedError.
        2) For each inequality, if previous cond defines part of the interval
           update the new conds interval.
           -  eg x < 1, x < 3 -> [oo,1],[1,3] instead of [oo,1],[oo,3]
        3) Sort the intervals to make it easier to find correct exprs

        Under normal use, we return the expr,cond pairs in increasing order
        along the real axis corresponding to the symbol sym.  If targetcond
        is given, we return a list of (lowerbound, upperbound) pairs for
        this condition.

        """
        from ...solvers.inequalities import solve_univariate_inequality
        default = None
        int_expr = []
        expr_cond = []
        or_cond = False
        or_intervals = []
        independent_expr_cond = []
        if isinstance(targetcond, Relational) and targetcond.has(sym):
            targetcond = solve_univariate_inequality(targetcond, sym)
        for expr, cond in self.args:  # pragma: no branch
            if isinstance(cond, Relational) and cond.has(sym):
                cond = solve_univariate_inequality(cond, sym)
            if isinstance(cond, Or):
                for cond2 in sorted(cond.args, key=default_sort_key):
                    expr_cond.append((expr, cond2))
            else:
                expr_cond.append((expr, cond))
            if cond == true:
                break
        for expr, cond in expr_cond:  # pragma: no branch
            if cond == true:
                independent_expr_cond.append((expr, cond))
                default = self.func(*independent_expr_cond)
                break
            orig_cond = cond
            if sym not in cond.free_symbols:
                independent_expr_cond.append((expr, cond))
                continue
            if isinstance(cond, Equality):
                continue
            if isinstance(cond, And):
                lower = -oo
                upper = oo
                for cond2 in cond.args:
                    if sym not in [cond2.lts, cond2.gts]:
                        cond2 = solve_univariate_inequality(cond2, sym)
                    if cond2.lts == sym:
                        upper = Min(cond2.gts, upper)
                    elif cond2.gts == sym:
                        lower = Max(cond2.lts, lower)
                    else:
                        raise NotImplementedError(
                            'Unable to handle interval evaluation of expression.')
            else:
                lower, upper = cond.lts, cond.gts  # part 1: initialize with givens
                if cond.lts == sym:                # part 1a: expand the side ...
                    lower = -oo   # e.g. x <= 0 ---> -oo <= 0
                elif cond.gts == sym:            # part 1a: ... that can be expanded
                    upper = oo           # e.g. x >= 0 --->  oo >= 0
                else:
                    raise NotImplementedError(
                        'Unable to handle interval evaluation of expression.')

            # part 1b: Reduce (-)infinity to what was passed in.
            lower, upper = Max(a, lower), Min(b, upper)

            for n, ni in enumerate(int_expr):
                # Part 2: remove any interval overlap.  For any conflicts, the
                # iterval already there wins, and the incoming interval updates
                # its bounds accordingly.
                if self.__eval_cond(lower < ni[1]) and \
                        self.__eval_cond(lower >= ni[0]):
                    lower = ni[1]
                elif len(ni[1].free_symbols) and \
                        self.__eval_cond(lower >= ni[0]):
                    if self.__eval_cond(lower == ni[0]):
                        lower = ni[1]
                    else:
                        int_expr[n][1] = Min(lower, ni[1])
                elif len(ni[0].free_symbols) and \
                        self.__eval_cond(upper == ni[1]):
                    upper = Min(upper, ni[0])
                elif len(ni[1].free_symbols) and \
                        (lower >= ni[0]) != true and \
                        (ni[1] == Min(lower, upper)) is not True:
                    upper = Min(upper, ni[0])
                elif self.__eval_cond(upper > ni[0]) and \
                        self.__eval_cond(upper <= ni[1]):
                    upper = ni[0]
                elif len(ni[0].free_symbols) and \
                        self.__eval_cond(upper < ni[1]):
                    int_expr[n][0] = Max(upper, ni[0])

            if self.__eval_cond(lower >= upper) is not True:  # Is it still an interval?
                int_expr.append([lower, upper, expr])
            if orig_cond == targetcond:
                return [(lower, upper, None)]
            if isinstance(targetcond, Or) and cond in targetcond.args:
                or_cond = Or(or_cond, cond)
                or_intervals.append((lower, upper, None))
                if or_cond == targetcond:
                    or_intervals.sort(key=lambda x: x[0])
                    return or_intervals

        int_expr.sort(key=lambda x: x[1].sort_key(
        ) if x[1].is_number else (-oo).sort_key())
        int_expr.sort(key=lambda x: x[0].sort_key(
        ) if x[0].is_number else oo.sort_key())

        for n, ni in enumerate(int_expr):
            if len(ni[0].free_symbols) or len(ni[1].free_symbols):
                if isinstance(ni[1], Min) or ni[1] == b:
                    newval = Min(*ni[:-1])
                    if n > 0 and ni[0] == int_expr[n - 1][1]:
                        int_expr[n - 1][1] = newval
                    int_expr[n][0] = newval
                else:
                    newval = Max(*ni[:-1])
                    if n < len(int_expr) - 1 and ni[1] == int_expr[n + 1][0]:
                        int_expr[n + 1][0] = newval
                    int_expr[n][1] = newval

        # Add holes to list of intervals if there is a default value,
        # otherwise raise a ValueError.
        holes = []
        curr_low = a
        for int_a, int_b, expr in int_expr:
            if (curr_low < int_a) == true:
                holes.append([curr_low, Min(b, int_a), default])
            elif (curr_low >= int_a) != true:
                holes.append([curr_low, Min(b, int_a), default])
            curr_low = Min(b, int_b)
        if (curr_low < b) == true:
            holes.append([Min(b, curr_low), b, default])
        elif (curr_low >= b) != true:
            holes.append([Min(b, curr_low), b, default])

        if holes and default is not None:
            int_expr.extend(holes)
            if targetcond == true:
                return [(h[0], h[1], None) for h in holes]

        return int_expr

    def _eval_nseries(self, x, n, logx):
        args = [(ec.expr._eval_nseries(x, n, logx),
                 ec.cond.limit(x, 0)) for ec in self.args]
        return self.func(*args)

    def _eval_power(self, other):
        return self.func(*[(e**other, c) for e, c in self.args])

    def _eval_subs(self, old, new):
        """Piecewise conditions may contain bool which are not of Basic type."""
        args = list(self.args)
        for i, (e, c) in enumerate(args):  # pragma: no branch
            c = c._subs(old, new)
            if c != false:
                e = e._subs(old, new)
            args[i] = e, c
            if c == true:
                return self.func(*args)

    def _eval_transpose(self):
        return self.func(*[(e.transpose(), c) for e, c in self.args])

    def _eval_template_is_attr(self, is_attr, when_multiple=None):
        b = None
        for expr, _ in self.args:
            a = getattr(expr, is_attr)
            if a is None:
                return
            if b is None:
                b = a
            elif b is not a:
                return when_multiple
        return b

    def _eval_is_finite(self):
        return self._eval_template_is_attr('is_finite')

    def _eval_is_complex(self):
        return self._eval_template_is_attr('is_complex')

    def _eval_is_imaginary(self):
        return self._eval_template_is_attr('is_imaginary')

    def _eval_is_integer(self):
        return self._eval_template_is_attr('is_integer')

    def _eval_is_irrational(self):
        return self._eval_template_is_attr('is_irrational')

    def _eval_is_negative(self):
        return self._eval_template_is_attr('is_negative')

    def _eval_is_nonnegative(self):
        return self._eval_template_is_attr('is_nonnegative')

    def _eval_is_nonpositive(self):
        return self._eval_template_is_attr('is_nonpositive')

    def _eval_is_nonzero(self):
        return self._eval_template_is_attr('is_nonzero')

    def _eval_is_odd(self):
        return self._eval_template_is_attr('is_odd')

    def _eval_is_polar(self):
        return self._eval_template_is_attr('is_polar')

    def _eval_is_positive(self):
        return self._eval_template_is_attr('is_positive')

    def _eval_is_extended_real(self):
        return self._eval_template_is_attr('is_extended_real')

    def _eval_is_zero(self):
        return self._eval_template_is_attr('is_zero')

    @classmethod
    def __eval_cond(cls, cond):
        """Return the truth value of the condition."""
        if cond == true:
            return True
        if isinstance(cond, Equality):
            diff = cond.lhs - cond.rhs
            if diff.is_commutative:
                return diff.is_zero

    def as_expr_set_pairs(self):
        exp_sets = []
        U = ExtendedReals
        for expr, cond in self.args:
            cond_int = U.intersection(cond.as_set())
            U = U - cond_int
            exp_sets.append((expr, cond_int))
        return exp_sets

    def _eval_rewrite_as_tractable(self, *args, wrt=None, **kwargs):
        if wrt is not None:
            for a, c in args:  # pragma: no branch
                if c.free_symbols - {wrt}:
                    return
                if c.limit(wrt, oo):
                    return a


def piecewise_fold(expr):
    """
    Takes an expression containing a piecewise function and returns the
    expression in piecewise form.

    Examples
    ========

    >>> p = Piecewise((x, x < 1), (1, True))
    >>> piecewise_fold(x*p)
    Piecewise((x**2, x < 1), (x, true))

    See Also
    ========

    diofant.functions.elementary.piecewise.Piecewise

    """
    if not isinstance(expr, Basic) or not expr.has(Piecewise):
        return expr
    new_args = list(map(piecewise_fold, expr.args))
    if isinstance(expr, ExprCondPair):
        return ExprCondPair(*new_args)
    piecewise_args = []
    for n, arg in enumerate(new_args):
        if isinstance(arg, Piecewise):
            piecewise_args.append(n)
    if len(piecewise_args) > 0:
        n = piecewise_args[0]
        new_args = [(expr.func(*(new_args[:n] + [e] + new_args[n + 1:])), c)
                    for e, c in new_args[n].args]
        if isinstance(expr, Boolean):
            # If expr is Boolean, we must return some kind of PiecewiseBoolean.
            # This is constructed by means of Or, And and Not.
            # piecewise_fold(0 < Piecewise( (sin(x), x<0), (cos(x), True)))
            # can't return Piecewise((0 < sin(x), x < 0), (0 < cos(x), True))
            # but instead Or(And(x < 0, 0 < sin(x)), And(0 < cos(x), Not(x<0)))
            other = True
            rtn = False
            for e, c in new_args:
                rtn = Or(rtn, And(other, c, e))
                other = And(other, Not(c))
            if len(piecewise_args) > 1:
                return piecewise_fold(rtn)
            return rtn
        if len(piecewise_args) > 1:
            return piecewise_fold(Piecewise(*new_args))
        return Piecewise(*new_args)
    return expr.func(*new_args)
