from ..core import (Add, Dummy, Expr, Integer, Mul, Symbol, Tuple, cacheit,
                    expand_log, expand_power_base, nan, oo, sympify)
from ..core.compatibility import is_sequence
from ..utilities import default_sort_key
from ..utilities.iterables import uniq


class Order(Expr):
    r"""Represents the limiting behavior of function.

    The formal definition for order symbol `O(f(x))` (Big O) is
    that `g(x) \in O(f(x))` as `x\to a` iff

    .. math:: \lim\limits_{x \rightarrow a} \sup
              \left|\frac{g(x)}{f(x)}\right| < \infty

    Parameters
    ==========

    expr : Expr
        an expression
    args : sequence of Symbol's or pairs (Symbol, Expr), optional
        If only symbols are provided, i.e. no limit point are
        passed, then the limit point is assumed to be zero.  If no
        symbols are passed then all symbols in the expression are used.

    Examples
    ========

    The order of a function can be intuitively thought of representing all
    terms of powers greater than the one specified.  For example, `O(x^3)`
    corresponds to any terms proportional to `x^3, x^4,\ldots` and any
    higher power.  For a polynomial, this leaves terms proportional
    to `x^2`, `x` and constants.

    >>> 1 + x + x**2 + x**3 + x**4 + O(x**3)
    1 + x + x**2 + O(x**3)

    ``O(f(x))`` is automatically transformed to ``O(f(x).as_leading_term(x))``:

    >>> O(x + x**2)
    O(x)
    >>> O(cos(x))
    O(1)

    Some arithmetic operations:

    >>> O(x)*x
    O(x**2)
    >>> O(x) - O(x)
    O(x)

    The Big O symbol is a set, so we support membership test:

    >>> x in O(x)
    True
    >>> O(1) in O(1, x)
    True
    >>> O(1, x) in O(1)
    False
    >>> O(x) in O(1, x)
    True
    >>> O(x**2) in O(x)
    True

    Limit points other then zero and multivariate Big O are also supported:

    >>> O(x) == O(x, (x, 0))
    True
    >>> O(x + x**2, (x, oo))
    O(x**2, (x, oo))
    >>> O(cos(x), (x, pi/2))
    O(x - pi/2, (x, pi/2))

    >>> O(1 + x*y)
    O(1, x, y)
    >>> O(1 + x*y, (x, 0), (y, 0))
    O(1, x, y)
    >>> O(1 + x*y, (x, oo), (y, oo))
    O(x*y, (x, oo), (y, oo))

    References
    ==========

    * https://en.wikipedia.org/wiki/Big_O_notation

    """

    is_Order = True

    @cacheit
    def __new__(cls, expr, *args, **kwargs):
        expr = sympify(expr)

        if not args:
            if expr.is_Order:
                variables = expr.variables
                point = expr.point
            else:
                variables = list(expr.free_symbols)
                point = [Integer(0)]*len(variables)
        else:
            args = list(args if is_sequence(args) else [args])
            variables, point = [], []
            if is_sequence(args[0]):
                for a in args:
                    v, p = list(map(sympify, a))
                    variables.append(v)
                    point.append(p)
            else:
                variables = list(map(sympify, args))
                point = [Integer(0)]*len(variables)

        if not all(isinstance(v, (Dummy, Symbol)) for v in variables):
            raise TypeError(f'Variables are not symbols, got {variables}')

        if len(list(uniq(variables))) != len(variables):
            raise ValueError(f'Variables are supposed to be unique symbols, got {variables}')

        if expr.is_Order:
            expr_vp = dict(expr.args[1:])
            new_vp = dict(expr_vp)
            vp = dict(zip(variables, point))
            for v, p in vp.items():
                if v in new_vp:
                    if p != new_vp[v]:
                        raise NotImplementedError(
                            'Mixing Order at different points is not supported.')
                else:
                    new_vp[v] = p
            if set(expr_vp) == set(new_vp):
                return expr
            variables = list(new_vp)
            point = [new_vp[v] for v in variables]

        if expr is nan:
            return nan

        if any(x in p.free_symbols for x in variables for p in point):
            raise ValueError(f'Got {point} as a point.')

        if variables:
            if any(p != point[0] for p in point):
                raise NotImplementedError
            if point[0].as_coefficient(oo):
                s = {k: 1/Dummy() for k in variables}
                rs = {1/v: 1/k for k, v in s.items()}
            elif point[0] != 0:
                s = {k: Dummy() + point[0] for k in variables}
                rs = {v - point[0]: k - point[0] for k, v in s.items()}
            else:
                s = ()
                rs = ()

            expr = expr.subs(s)

            if expr.is_Add:
                from ..core import expand_multinomial
                expr = expand_multinomial(expr)

            if s:
                args = tuple(r[0] for r in rs.items())
            else:
                args = tuple(variables)

            if len(variables) > 1:
                # XXX: better way?  We need this expand() to
                # workaround e.g: expr = x*(x + y).
                # (x*(x + y)).as_leading_term(x, y) currently returns
                # x*y (wrong order term!).  That's why we want to deal with
                # expand()'ed expr (handled in "if expr.is_Add" branch below).
                expr = expr.expand()

            if expr.is_Add:
                lst = expr.extract_leading_order(args)
                expr = Add(*[f.expr for (e, f) in lst])

            elif expr:
                expr = expr.as_leading_term(*args)
                expr = expr.as_independent(*args, as_Add=False)[1]

                expr = expand_power_base(expr)
                expr = expand_log(expr)

                if len(args) == 1:
                    # The definition of O(f(x)) symbol explicitly stated that
                    # the argument of f(x) is irrelevant.  That's why we can
                    # combine some power exponents (only "on top" of the
                    # expression tree for f(x)), e.g.:
                    # x**p * (-x)**q -> x**(p+q) for real p, q.
                    x = args[0]
                    margs = list(Mul.make_args(
                        expr.as_independent(x, as_Add=False)[1]))

                    for i, t in enumerate(margs):
                        if t.is_Pow:
                            b, q = t.base, t.exp
                            if b in (x, -x) and q.is_extended_real and not q.has(x):
                                margs[i] = x**q
                            elif b.is_Pow and not b.exp.has(x):
                                b, r = b.base, b.exp
                                if b in (x, -x) and r.is_extended_real:
                                    margs[i] = x**(r*q)
                            elif b.is_Mul and b.args[0] == -1:
                                b = -b
                                if b.is_Pow and not b.exp.has(x):
                                    b, r = b.base, b.exp
                                    if b in (x, -x) and r.is_extended_real:
                                        margs[i] = x**(r*q)

                    expr = Mul(*margs)

            expr = expr.subs(rs)

        if expr == 0:
            return expr

        if expr.is_Order:
            expr = expr.expr

        if not expr.has(*variables):
            expr = Integer(1)

        # create Order instance:
        vp = dict(zip(variables, point))
        variables.sort(key=default_sort_key)
        point = [vp[v] for v in variables]
        args = (expr,) + Tuple(*zip(variables, point))
        obj = Expr.__new__(cls, *args)
        return obj

    def _eval_nseries(self, x, n, logx):
        return self

    @property
    def expr(self):
        return self.args[0]

    @property
    def variables(self):
        if self.args[1:]:
            return tuple(x[0] for x in self.args[1:])
        return ()

    @property
    def point(self):
        if self.args[1:]:
            return tuple(x[1] for x in self.args[1:])
        return ()

    @property
    def free_symbols(self):
        return self.expr.free_symbols | set(self.variables)

    def _eval_power(self, other):
        if other.is_Number and other.is_nonnegative:
            return self.func(self.expr**other, *self.args[1:])
        if other == O(1):
            return self

    def as_expr_variables(self, order_symbols):
        if order_symbols is None:
            order_symbols = self.args[1:]
        else:
            if (not all(o[1] == order_symbols[0][1] for o in order_symbols) and
                    not all(p == self.point[0] for p in self.point)):
                raise NotImplementedError('Order at points other than 0 '
                                          f'or oo not supported, got {self.point} as a point.')
            if order_symbols and order_symbols[0][1] != self.point[0]:
                raise NotImplementedError(
                    'Multiplying Order at different points is not supported.')
            order_symbols = dict(order_symbols)
            for s, p in dict(self.args[1:]).items():
                if s not in order_symbols:
                    order_symbols[s] = p
            order_symbols = sorted(order_symbols.items(), key=lambda x: default_sort_key(x[0]))
        return self.expr, tuple(order_symbols)

    def removeO(self):
        return Integer(0)

    def getO(self):
        return self

    @cacheit
    def contains(self, expr):
        """Membership test.

        Returns
        =======

        Boolean or None
            Return True if ``expr`` belongs to ``self``.  Return False if
            ``self`` belongs to ``expr``.  Return None if the inclusion
            relation cannot be determined.

        """
        from ..simplify import powsimp
        from .limits import Limit
        if expr == 0:
            return True
        if expr is nan:
            return False
        if expr.is_Order:
            if (not all(p == expr.point[0] for p in expr.point) and
                    not all(p == self.point[0] for p in self.point)):
                raise NotImplementedError('Order at points other than 0 '
                                          f'or oo not supported, got {self.point} as a point.')
            # self and/or expr is O(1):
            if any(not p for p in [expr.point, self.point]):
                point = self.point + expr.point
                if point:
                    point = point[0]
                else:
                    point = Integer(0)
            else:
                point = self.point[0]

            if expr.expr == self.expr:
                # O(1) + O(1), O(1) + O(1, x), etc.
                return all(x in self.args[1:] for x in expr.args[1:])
            if expr.expr.is_Add:
                return all(self.contains(x) for x in expr.expr.args)
            if self.expr.is_Add and point == 0:
                return any(self.func(x, *self.args[1:]).contains(expr)
                           for x in self.expr.args)
            if self.variables and expr.variables:
                common_symbols = tuple(s for s in self.variables if s in expr.variables)
            elif self.variables:
                common_symbols = self.variables
            else:
                common_symbols = expr.variables
            if not common_symbols:
                return
            r = None
            ratio = self.expr/expr.expr
            ratio = powsimp(ratio, deep=True, combine='exp')
            for s in common_symbols:
                l = Limit(ratio, s, point).doit(heuristics=False)
                if not isinstance(l, Limit):
                    l = l != 0
                else:
                    l = None
                if r is None:
                    r = l
                else:
                    if r != l:
                        return
            return r
        obj = self.func(expr, *self.args[1:])
        return self.contains(obj)

    def __contains__(self, other):
        result = self.contains(other)
        if result is None:
            raise TypeError('contains did not evaluate to a bool')
        return result

    def _eval_subs(self, old, new):
        if old in self.variables:
            newexpr = self.expr.subs({old: new})
            i = self.variables.index(old)
            newvars = list(self.variables)
            newpt = list(self.point)
            if new.is_Symbol:
                newvars[i] = new
            else:
                syms = new.free_symbols
                if len(syms) == 1 or old in syms:
                    if old in syms:
                        var = self.variables[i]
                    else:
                        var = syms.pop()
                    # First, try to substitute self.point in the "new"
                    # expr to see if this is a fixed point.
                    # E.g.  O(y).subs({y: sin(x)})
                    point = new.subs({var: self.point[i]})
                    if point != self.point[i]:
                        from ..solvers import solve
                        d = Dummy()
                        res = solve(old - new.subs({var: d}), d)
                        point = d.subs(res[0]).limit(old, self.point[i])
                    newvars[i] = var
                    newpt[i] = point
                else:
                    del newvars[i], newpt[i]
                    if not syms and new == self.point[i]:
                        newvars.extend(syms)
                        newpt.extend([Integer(0)]*len(syms))
            return Order(newexpr, *zip(newvars, newpt))

    def _eval_conjugate(self):
        expr = self.expr._eval_conjugate()
        if expr is not None:
            return self.func(expr, *self.args[1:])

    def _eval_transpose(self):
        expr = self.expr._eval_transpose()
        if expr is not None:
            return self.func(expr, *self.args[1:])

    def _eval_is_commutative(self):
        return self.expr.is_commutative


O = Order
