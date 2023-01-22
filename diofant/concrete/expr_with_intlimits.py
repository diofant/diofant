from ..core import oo
from .expr_with_limits import ExprWithLimits


class ReorderError(NotImplementedError):
    """Exception raised when trying to reorder dependent limits."""

    def __init__(self, expr, msg):
        """Initialize self."""
        super().__init__(f'{expr} could not be reordered: {msg}.')


class ExprWithIntLimits(ExprWithLimits):
    """Represents an expression with integer limits."""

    def __init__(self, function, *symbols, **assumptions):
        """Initialize self."""
        if not all(all(abs(_) == oo or (_.is_integer is not False)
                       for _ in l[1:]) for l in self.limits):
            raise ValueError('Limits must be integers or ±oo.')

    def change_index(self, var, trafo, newvar=None):
        r"""
        Change index of a Sum or Product.

        Perform a linear transformation `x \mapsto a x + b` on the index variable
        `x`. For `a` the only values allowed are `\pm 1`. A new variable to be used
        after the change of index can also be specified.

        Parameters
        ==========

        var : Symbol
            specifies the index variable `x` to transform.
        trafo : Expr
            The linear transformation in terms of ``var``.
        newvar : Symbol, optional
            Replacement symbol to be used instead of
            ``var`` in the final expression.

        Examples
        ========

        >>> from diofant.abc import i, j, l, u, v

        >>> s = Sum(x, (x, a, b))
        >>> s.doit()
        -a**2/2 + a/2 + b**2/2 + b/2

        >>> sn = s.change_index(x, x + 1, y)
        >>> sn
        Sum(y - 1, (y, a + 1, b + 1))
        >>> sn.doit()
        -a**2/2 + a/2 + b**2/2 + b/2

        >>> sn = s.change_index(x, -x, y)
        >>> sn
        Sum(-y, (y, -b, -a))
        >>> sn.doit()
        -a**2/2 + a/2 + b**2/2 + b/2

        >>> sn = s.change_index(x, x+u)
        >>> sn
        Sum(-u + x, (x, a + u, b + u))
        >>> sn.doit()
        -a**2/2 - a*u + a/2 + b**2/2 + b*u + b/2 - u*(-a + b + 1) + u
        >>> simplify(sn.doit())
        -a**2/2 + a/2 + b**2/2 + b/2

        >>> sn = s.change_index(x, -x - u, y)
        >>> sn
        Sum(-u - y, (y, -b - u, -a - u))
        >>> sn.doit()
        -a**2/2 - a*u + a/2 + b**2/2 + b*u + b/2 - u*(-a + b + 1) + u
        >>> simplify(sn.doit())
        -a**2/2 + a/2 + b**2/2 + b/2

        >>> p = Product(i*j**2, (i, a, b), (j, c, d))
        >>> p
        Product(i*j**2, (i, a, b), (j, c, d))
        >>> p2 = p.change_index(i, i+3, k)
        >>> p2
        Product(j**2*(k - 3), (k, a + 3, b + 3), (j, c, d))
        >>> p3 = p2.change_index(j, -j, l)
        >>> p3
        Product(l**2*(k - 3), (k, a + 3, b + 3), (l, -d, -c))

        When dealing with symbols only, we can make a
        general linear transformation:

        >>> sn = s.change_index(x, u*x+v, y)
        >>> sn
        Sum((-v + y)/u, (y, b*u + v, a*u + v))
        >>> sn.doit()
        -v*(a*u - b*u + 1)/u + (a**2*u**2/2 + a*u*v + a*u/2 - b**2*u**2/2 - b*u*v + b*u/2 + v)/u
        >>> simplify(sn.doit())
        a**2*u/2 + a/2 - b**2*u/2 + b/2

        However, the last result can be inconsistent with usual
        summation where the index increment is always 1. This is
        obvious as we get back the original value only for ``u``
        equal +1 or -1.

        See Also
        ========

        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.index
        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder_limit
        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder
        diofant.concrete.summations.Sum.reverse_order
        diofant.concrete.products.Product.reverse_order

        """
        if newvar is None:
            newvar = var

        limits = []
        for limit in self.limits:
            if limit[0] == var:
                p = trafo.as_poly(var)
                if p.degree() != 1:
                    raise ValueError('Index transformation is not linear')
                alpha = p.coeff_monomial(var)
                beta = p.coeff_monomial(1)
                if alpha.is_number:
                    if alpha == 1:
                        limits.append((newvar, alpha*limit[1] + beta, alpha*limit[2] + beta))
                    elif alpha == -1:
                        limits.append((newvar, alpha*limit[2] + beta, alpha*limit[1] + beta))
                    else:
                        raise ValueError('Linear transformation results in non-linear summation stepsize')
                else:
                    # Note that the case of alpha being symbolic can give issues if alpha < 0.
                    limits.append((newvar, alpha*limit[2] + beta, alpha*limit[1] + beta))
            else:
                limits.append(limit)

        function = self.function.subs({var: (var - beta)/alpha})
        function = function.subs({var: newvar})

        return self.func(function, *limits)

    def index(self, x):
        """
        Return the index of a dummy variable in the list of limits.

        Note that we start counting with 0 at the inner-most
        limits tuple.

        Parameters
        ==========

        x : Symbol
            a dummy variable

        Examples
        ========

        >>> Sum(x*y, (x, a, b), (y, c, d)).index(x)
        0
        >>> Sum(x*y, (x, a, b), (y, c, d)).index(y)
        1
        >>> Product(x*y, (x, a, b), (y, c, d)).index(x)
        0
        >>> Product(x*y, (x, a, b), (y, c, d)).index(y)
        1

        See Also
        ========

        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder_limit
        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder
        diofant.concrete.summations.Sum.reverse_order
        diofant.concrete.products.Product.reverse_order

        """
        variables = [limit[0] for limit in self.limits]

        if variables.count(x) != 1:
            raise ValueError(self, 'Number of instances of variable not equal to one')
        return variables.index(x)

    def reorder(self, *arg):
        r"""
        Reorder limits in a expression containing a Sum or a Product.

        Parameters
        ==========

        \*arg : list of tuples
            These tuples can contain numerical indices or index
            variable names or involve both.

        Examples
        ========

        >>> from diofant.abc import e, f

        >>> Sum(x*y, (x, a, b), (y, c, d)).reorder((x, y))
        Sum(x*y, (y, c, d), (x, a, b))

        >>> Sum(x*y*z, (x, a, b), (y, c, d), (z, e, f)).reorder((x, y), (x, z), (y, z))
        Sum(x*y*z, (z, e, f), (y, c, d), (x, a, b))

        >>> P = Product(x*y*z, (x, a, b), (y, c, d), (z, e, f))
        >>> P.reorder((x, y), (x, z), (y, z))
        Product(x*y*z, (z, e, f), (y, c, d), (x, a, b))

        We can also select the index variables by counting them, starting
        with the inner-most one:

        >>> Sum(x**2, (x, a, b), (x, c, d)).reorder((0, 1))
        Sum(x**2, (x, c, d), (x, a, b))

        And of course we can mix both schemes:

        >>> Sum(x*y, (x, a, b), (y, c, d)).reorder((y, x))
        Sum(x*y, (y, c, d), (x, a, b))
        >>> Sum(x*y, (x, a, b), (y, c, d)).reorder((y, 0))
        Sum(x*y, (y, c, d), (x, a, b))

        See Also
        ========

        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.index
        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder_limit
        diofant.concrete.summations.Sum.reverse_order
        diofant.concrete.products.Product.reverse_order

        """
        new_expr = self

        for r in arg:
            if len(r) != 2:
                raise ValueError(r, 'Invalid number of arguments')

            index1 = r[0]
            index2 = r[1]

            if not isinstance(r[0], int):
                index1 = self.index(r[0])
            if not isinstance(r[1], int):
                index2 = self.index(r[1])

            new_expr = new_expr.reorder_limit(index1, index2)

        return new_expr

    def reorder_limit(self, x, y):
        """
        Interchange two limit tuples of a Sum or Product expression.

        Parameters
        ==========

        x, y: int
            are integers corresponding to the index
            variables of the two limits which are to be interchanged.

        Examples
        ========

        >>> from diofant.abc import e, f

        >>> Sum(x*y*z, (x, a, b), (y, c, d), (z, e, f)).reorder_limit(0, 2)
        Sum(x*y*z, (z, e, f), (y, c, d), (x, a, b))
        >>> Sum(x**2, (x, a, b), (x, c, d)).reorder_limit(1, 0)
        Sum(x**2, (x, c, d), (x, a, b))

        >>> Product(x*y*z, (x, a, b), (y, c, d), (z, e, f)).reorder_limit(0, 2)
        Product(x*y*z, (z, e, f), (y, c, d), (x, a, b))

        See Also
        ========

        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.index
        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder
        diofant.concrete.summations.Sum.reverse_order
        diofant.concrete.products.Product.reverse_order

        """
        var = {limit[0] for limit in self.limits}
        limit_x = self.limits[x]
        limit_y = self.limits[y]

        if (len(set(limit_x[1].free_symbols).intersection(var)) == 0 and
                len(set(limit_x[2].free_symbols).intersection(var)) == 0 and
                len(set(limit_y[1].free_symbols).intersection(var)) == 0 and
                len(set(limit_y[2].free_symbols).intersection(var)) == 0):

            limits = []
            for i, limit in enumerate(self.limits):
                if i == x:
                    limits.append(limit_y)
                elif i == y:
                    limits.append(limit_x)
                else:
                    limits.append(limit)

            return type(self)(self.function, *limits)
        raise ReorderError(self, 'could not interchange the two limits specified')
