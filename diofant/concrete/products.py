from ..core import Integer, Mul, nan
from ..functions import exp, log
from ..polys import quo, roots
from ..simplify.powsimp import powsimp
from .expr_with_intlimits import ExprWithIntLimits


class Product(ExprWithIntLimits):
    r"""
    Represents an unevaluated products.

    ``Product`` represents a finite or infinite product, with the first
    argument being the general form of terms in the series (which
    usually depend on the bound variable ``symbol``), and the second
    argument being ``(symbol, start, end)``, with ``symbol`` taking
    all integer values from ``start`` through ``end`` (inclusive).

    Notes
    =====

    We follow the the analogue of the summation convention described by
    Karr :cite:`Karr1981summation`, adopted by the
    :class:`~diofant.concrete.summations.Sum`:

    .. math::

        \prod\limits_{i=m}^n f_i = \frac{1}{\prod\limits_{i=n+1}^{m-1}f_i}

    Examples
    ========

    >>> Product(k**2, (k, 1, m)).doit()
    factorial(m)**2

    Products with the lower limit being larger than the upper one:

    >>> Product(1/k, (k, 6, 1)).doit()
    120
    >>> Product(k, (k, 2, 5)).doit()
    120

    The empty product:

    >>> Product(k, (k, n, n-1)).doit()
    1

    See Also
    ========

    diofant.concrete.summations.Sum
    diofant.concrete.summations.summation
    product

    References
    ==========

    * https://en.wikipedia.org/wiki/Multiplication#Capital_pi_notation
    * https://en.wikipedia.org/wiki/Empty_product

    """

    def _eval_rewrite_as_Sum(self, *args):
        from .summations import Sum
        return exp(Sum(log(self.function), *self.limits))

    @property
    def term(self):
        return self.args[0]
    function = term

    def _eval_is_zero(self):
        if self.term.is_zero:
            return True

    def doit(self, **hints):
        f = self.function
        for index, limit in enumerate(self.limits):
            i, a, b = limit
            dif = b - a
            if dif.is_Integer and dif < 0:
                a, b = b + 1, a - 1
                f = 1 / f

            g = self._eval_product(f, (i, a, b))
            if g in (None, nan):
                return self.func(powsimp(f), *self.limits[index:])
            f = g

        if hints.get('deep', True):
            return f.doit(**hints)
        return powsimp(f)

    def _eval_adjoint(self):
        if self.is_commutative:
            return self.func(self.function.adjoint(), *self.limits)

    def _eval_conjugate(self):
        return self.func(self.function.conjugate(), *self.limits)

    def _eval_product(self, term, limits):
        from ..functions import KroneckerDelta, RisingFactorial
        from .delta import _has_simple_delta, deltaproduct
        from .summations import summation

        (k, a, n) = limits

        if k not in term.free_symbols:
            if (term - 1).is_zero:
                return Integer(1)
            return term**(n - a + 1)

        if a == n:
            return term.subs({k: a})

        if term.has(KroneckerDelta) and _has_simple_delta(term, limits[0]):
            return deltaproduct(term, limits)

        dif = n - a
        if dif.is_Integer:
            return Mul(*[term.subs({k: a + i}) for i in range(dif + 1)])

        if term.is_polynomial(k):
            poly = term.as_poly(k)

            A = B = Q = Integer(1)

            all_roots = roots(poly)

            M = 0
            for r, m in all_roots.items():
                M += m
                A *= RisingFactorial(a - r, n - a + 1)**m
                Q *= (k - r)**m

            if 0 < M < poly.degree():
                arg = quo(poly, Q.as_poly(k)).as_expr(k)
                B = self.func(arg, (k, a, n)).doit()

            if M > 0:
                return poly.LC()**(n - a + 1) * A * B

        elif term.is_Add:
            p, q = term.as_numer_denom()
            q = self._eval_product(q, (k, a, n))
            if q.is_Number:
                # There is expression, which couldn't change by
                # as_numer_denom(). E.g. n**(2/3) + 1 --> (n**(2/3) + 1, 1).
                p = sum(self._eval_product(i, (k, a, n)) for i in p.as_coeff_Add())
            else:
                p = self._eval_product(p, (k, a, n))
            return p / q

        elif term.is_Mul:
            exclude, include = [], []

            for t in term.args:
                p = self._eval_product(t, (k, a, n))

                if p is not None:
                    exclude.append(p)
                else:
                    include.append(t)

            if exclude:
                arg = term._new_rawargs(*include)
                A = Mul(*exclude)
                B = self.func(arg, (k, a, n)).doit()
                return A * B

        elif term.is_Pow:
            if not term.base.has(k):
                s = summation(term.exp, (k, a, n))

                return term.base**s
            if not term.exp.has(k):
                p = self._eval_product(term.base, (k, a, n))

                if p is not None:
                    return p**term.exp

    def _eval_simplify(self, ratio, measure):
        from ..simplify.simplify import product_simplify
        return product_simplify(self)

    def _eval_transpose(self):
        if self.is_commutative:
            return self.func(self.function.transpose(), *self.limits)

    def reverse_order(self, *indices):
        r"""
        Reverse the order of a limit in a Product.

        Parameters
        ==========

        \*indices : list
            The selectors in the argument ``indices`` specify some indices whose
            limits get reversed.  These selectors are either variable names or
            numerical indices counted starting from the inner-most limit tuple.

        Examples
        ========

        >>> P = Product(x, (x, a, b))
        >>> Pr = P.reverse_order(x)
        >>> Pr
        Product(1/x, (x, b + 1, a - 1))
        >>> Pr = Pr.doit()
        >>> Pr
        1/RisingFactorial(b + 1, a - b - 1)
        >>> simplify(Pr)
        gamma(b + 1)/gamma(a)
        >>> P = P.doit()
        >>> P
        RisingFactorial(a, -a + b + 1)
        >>> simplify(P)
        gamma(b + 1)/gamma(a)

        While one should prefer variable names when specifying which limits
        to reverse, the index counting notation comes in handy in case there
        are several symbols with the same name.

        >>> s = Sum(x*y, (x, a, b), (y, c, d))
        >>> s
        Sum(x*y, (x, a, b), (y, c, d))
        >>> s0 = s.reverse_order(0)
        >>> s0
        Sum(-x*y, (x, b + 1, a - 1), (y, c, d))
        >>> s1 = s0.reverse_order(1)
        >>> s1
        Sum(x*y, (x, b + 1, a - 1), (y, d + 1, c - 1))

        Of course we can mix both notations:

        >>> Sum(x*y, (x, a, b), (y, 2, 5)).reverse_order(x, 1)
        Sum(x*y, (x, b + 1, a - 1), (y, 6, 1))
        >>> Sum(x*y, (x, a, b), (y, 2, 5)).reverse_order(y, x)
        Sum(x*y, (x, b + 1, a - 1), (y, 6, 1))

        See Also
        ========

        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.index,
        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder_limit,
        diofant.concrete.expr_with_intlimits.ExprWithIntLimits.reorder

        References
        ==========

        * :cite:`Karr1981summation`

        """
        l_indices = list(indices)

        for i, indx in enumerate(l_indices):
            if not isinstance(indx, int):
                l_indices[i] = self.index(indx)

        e = 1
        limits = []
        for i, limit in enumerate(self.limits):
            l = limit
            if i in l_indices:
                e = -e
                l = (limit[0], limit[2] + 1, limit[1] - 1)
            limits.append(l)

        return Product(self.function ** e, *limits)


def product(*args, **kwargs):
    r"""
    Compute the product.

    The notation for symbols is similar to the notation used in Sum or
    Integral. product(f, (i, a, b)) computes the product of f with
    respect to i from a to b, i.e.,

    ::

                                     b
                                   _____
        product(f(n), (i, a, b)) = |   | f(n)
                                   |   |
                                   i = a

    If it cannot compute the product, it returns an unevaluated Product object.
    Repeated products can be computed by introducing additional symbols tuples::

    >>> i = symbols('i', integer=True)

    >>> product(i, (i, 1, k))
    factorial(k)
    >>> product(m, (i, 1, k))
    m**k
    >>> product(i, (i, 1, k), (k, 1, n))
    Product(factorial(k), (k, 1, n))

    """
    prod = Product(*args, **kwargs)

    if isinstance(prod, Product):
        return prod.doit(deep=False)
    return prod
