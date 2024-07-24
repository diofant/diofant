from ..core import (Add, Dummy, Expr, Integer, Mul, Symbol, cacheit,
                    expand_log, expand_power_base, nan, oo, sympify)


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
    var : Symbol, optional
        a variable of the expr.  If not privided, the expression is assumed
        to be univariate and it's variable is used.
    point : Expr, optional
        a limit point, default is zero.

    Examples
    ========

    The order of a function can be intuitively thought of representing all
    terms of powers greater than the one specified.  For example, `O(x^3)`
    corresponds to any terms proportional to `x^3, x^4,\ldots` and any
    higher power:

    >>> 1 + x + x**2 + x**3 + x**4 + O(x**3)
    1 + x + x**2 + O(x**3)

    ``O(f(x))`` is automatically transformed to ``O(f(x).as_leading_term(x))``:

    >>> O(x + x**2)
    O(x)
    >>> O(cos(x))
    O(1, x)

    Some arithmetic operations:

    >>> O(x)*x
    O(x**2)
    >>> O(x) - O(x)
    O(x)

    The Big O symbol is a set, so we support membership test:

    >>> x in O(x)
    True
    >>> O(x) in O(1, x)
    True
    >>> O(x**2) in O(x)
    True

    Limit points other then zero are also supported:

    >>> O(x) == O(x, x, 0)
    True
    >>> O(x + x**2, x, oo)
    O(x**2, x, oo)
    >>> O(cos(x), x, pi/2)
    O(x - pi/2, x, pi/2)

    References
    ==========

    * https://en.wikipedia.org/wiki/Big_O_notation

    """

    is_Order = True

    @cacheit
    def __new__(cls, expr, var=None, point=0, **kwargs):
        expr = sympify(expr)
        point = sympify(point)

        if not var:
            if expr.is_Order:
                var = expr.var
            else:
                free_symbols = expr.free_symbols
                if len(free_symbols) != 1:
                    raise ValueError
                var = free_symbols.pop()

        if not isinstance(var, (Dummy, Symbol)):
            raise TypeError(f'Variables are not symbols, got {var}')

        if expr is nan:
            return nan

        if expr.is_Order:
            if point != expr.point:
                raise NotImplementedError('Mixing Order at different points is not supported.')
            expr = expr.expr

        if var in point.free_symbols:
            raise ValueError(f'Got {point} as a point.')

        if point.as_coefficient(oo):
            s = {var: 1/Dummy()}
            rs = {1/v: 1/k for k, v in s.items()}
        elif point != 0:
            s = {var: Dummy() + point}
            rs = {v - point: k - point for k, v in s.items()}
        else:
            s = ()
            rs = ()

        expr = expr.subs(s)

        if s:
            x = tuple(rs)[0]
        else:
            x = var

        if expr.is_Add:
            lst = expr._extract_leading_order([x])
            expr = Add(*[f.expr for (e, f) in lst])

        elif expr:
            expr = expr.as_leading_term(x)
            _, expr = expr.as_independent(x, as_Add=False)

            expr = expand_power_base(expr)
            expr = expand_log(expr)

            # The definition of O(f(x)) symbol explicitly stated that
            # the argument of f(x) is irrelevant.  That's why we can
            # combine some power exponents (only "on top" of the
            # expression tree for f(x)), e.g.:
            # x**p * (-x)**q -> x**(p+q) for real p, q.
            margs = list(Mul.make_args(expr.as_independent(x, as_Add=False)[1]))

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

        if not expr.has(var):
            expr = Integer(1)

        args = (expr, var, point)
        obj = Expr.__new__(cls, *args)
        return obj

    def _eval_nseries(self, x, n, logx):
        return self

    @property
    def expr(self):
        return self.args[0]

    @property
    def var(self):
        return self.args[1]

    @property
    def point(self):
        return self.args[2]

    @property
    def free_symbols(self):
        return self.expr.free_symbols | {self.var}

    def _eval_power(self, other):
        if other.is_Number and other.is_nonnegative:
            return self.func(self.expr**other, self.var, self.point)

    def as_expr_variables(self, order_symbols):
        if order_symbols is None:
            order_symbols = self.args[1:]
        else:
            if order_symbols[1] != self.point:
                raise NotImplementedError('Multiplying Order at different points is not supported.')
            if order_symbols[0] != self.var:
                raise NotImplementedError
        return self.expr, order_symbols

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
            if self.point != expr.point:
                return
            if expr.expr == self.expr:
                return all(x in self.args[1:] for x in expr.args[1:])
            if self.var != expr.var:
                return
            ratio = self.expr/expr.expr
            ratio = powsimp(ratio, deep=True, combine='exp')
            l = Limit(ratio, self.var, self.point).doit(heuristics=False)
            if not isinstance(l, Limit):
                l = l != 0
            else:
                l = None
            return l
        obj = self.func(expr, *self.args[1:])
        return self.contains(obj)

    def __contains__(self, other):
        result = self.contains(other)
        if result is None:
            raise TypeError('contains did not evaluate to a bool')
        return result

    def _eval_subs(self, old, new):
        if old == self.var:
            newexpr = self.expr.subs({old: new})
            newvar = self.var
            newpt = self.point
            if new.is_Symbol:
                newvar = new
                newpt = self.point
            else:
                syms = new.free_symbols
                if len(syms) == 1 or old in syms:
                    if old in syms:
                        var = self.var
                    else:
                        var = syms.pop()
                    # First, try to substitute self.point in the "new"
                    # expr to see if this is a fixed point.
                    # E.g.  O(y).subs({y: sin(x)})
                    point = new.subs({var: self.point})
                    if point != self.point:
                        from ..solvers import solve
                        d = Dummy()
                        res = solve(old - new.subs({var: d}), d)
                        point = d.subs(res[0]).limit(old, self.point)
                    newvar = var
                    newpt = point
                if not syms and new == self.point:
                    newvar = old
                    newpt = new
            return Order(newexpr, newvar, newpt)

    def _eval_conjugate(self):
        expr = self.expr._eval_conjugate()
        if expr is not None:
            return self.func(expr, self.var, self.point)

    def _eval_transpose(self):
        expr = self.expr._eval_transpose()
        if expr is not None:
            return self.func(expr, self.var, self.point)

    def _eval_is_commutative(self):
        return self.expr.is_commutative


O = Order
