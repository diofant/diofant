from ..core import Derivative, Dummy, Eq, Function, Integer, Wild, nan, oo
from ..functions import Piecewise
from ..logic import false
from ..polys import PolynomialError, apart
from ..solvers import solve
from .expr_with_intlimits import ExprWithIntLimits
from .expr_with_limits import AddWithLimits
from .gosper import gosper_sum


class Sum(AddWithLimits, ExprWithIntLimits):
    r"""Represents unevaluated summation.

    ``Sum`` represents a finite or infinite series, with the first argument
    being the general form of terms in the series, and the second argument
    being ``(dummy_variable, start, end)``, with ``dummy_variable`` taking
    all integer values from ``start`` through ``end``. In accordance with
    long-standing mathematical convention, the end term is included in the
    summation.

    For finite sums (and sums with symbolic limits assumed to be finite) we
    follow the summation convention described by Karr :cite:`Karr1981summation`, especially
    definition 3 of section 1.4. The sum:

    .. math::

        \sum_{m \leq i < n} f(i)

    has *the obvious meaning* for `m < n`, namely:

    .. math::

        \sum_{m \leq i < n} f(i) = f(m) + f(m+1) + \ldots + f(n-2) + f(n-1)

    with the upper limit value `f(n)` excluded. The sum over an empty set is
    zero if and only if `m = n`:

    .. math::

        \sum_{m \leq i < n} f(i) = 0  \quad \mathrm{for} \quad  m = n

    Finally, for all other sums over empty sets we assume the following
    definition:

    .. math::

        \sum_{m \leq i < n} f(i) = - \sum_{n \leq i < m} f(i)  \quad \mathrm{for} \quad  m > n

    It is important to note that Karr defines all sums with the upper
    limit being exclusive. This is in contrast to the usual mathematical notation,
    but does not affect the summation convention. Indeed we have:

    .. math::

        \sum_{m \leq i < n} f(i) = \sum_{i = m}^{n - 1} f(i)

    where the difference in notation is intentional to emphasize the meaning,
    with limits typeset on the top being inclusive.

    Examples
    ========

    >>> from diofant.abc import i

    >>> Sum(k, (k, 1, m))
    Sum(k, (k, 1, m))
    >>> Sum(k, (k, 1, m)).doit()
    m**2/2 + m/2
    >>> Sum(k**2, (k, 1, m))
    Sum(k**2, (k, 1, m))
    >>> Sum(k**2, (k, 1, m)).doit()
    m**3/3 + m**2/2 + m/6
    >>> Sum(x**k, (k, 0, oo))
    Sum(x**k, (k, 0, oo))
    >>> Sum(x**k, (k, 0, oo)).doit()
    Piecewise((1/(-x + 1), Abs(x) < 1), (Sum(x**k, (k, 0, oo)), true))
    >>> Sum(x**k/factorial(k), (k, 0, oo)).doit()
    E**x

    Here are examples to do summation with symbolic indices.  You
    can use either Function of IndexedBase classes:

    >>> f = Function('f')

    >>> Sum(f(n), (n, 0, 3)).doit()
    f(0) + f(1) + f(2) + f(3)
    >>> Sum(f(n), (n, 0, oo)).doit()
    Sum(f(n), (n, 0, oo))
    >>> f = IndexedBase('f')
    >>> Sum(f[n]**2, (n, 0, 3)).doit()
    f[0]**2 + f[1]**2 + f[2]**2 + f[3]**2

    An example showing that the symbolic result of a summation is still
    valid for seemingly nonsensical values of the limits. Then the Karr
    convention allows us to give a perfectly valid interpretation to
    those sums by interchanging the limits according to the above rules:

    >>> s = Sum(i, (i, 1, n)).doit()
    >>> s
    n**2/2 + n/2
    >>> s.subs({n: -4})
    6
    >>> Sum(i, (i, 1, -4)).doit()
    6
    >>> Sum(-i, (i, -3, 0)).doit()
    6

    An explicit example of the Karr summation convention:

    >>> s1 = Sum(i**2, (i, m, m+n-1)).doit()
    >>> s1
    m**2*n + m*n**2 - m*n + n**3/3 - n**2/2 + n/6
    >>> s2 = Sum(i**2, (i, m+n, m-1)).doit()
    >>> s2
    -m**2*n - m*n**2 + m*n - n**3/3 + n**2/2 - n/6
    >>> s1 + s2
    0
    >>> s3 = Sum(i, (i, m, m-1)).doit()
    >>> s3
    0

    See Also
    ========

    diofant.concrete.summations.summation
    diofant.concrete.products.Product
    diofant.concrete.products.product

    References
    ==========

    * https://en.wikipedia.org/wiki/Summation#Capital-sigma_notation
    * https://en.wikipedia.org/wiki/Empty_sum

    """

    def __new__(cls, function, *symbols, **assumptions):
        obj = AddWithLimits.__new__(cls, function, *symbols, **assumptions)
        if not hasattr(obj, 'limits'):
            return obj
        if any(len(l) != 3 or None in l for l in obj.limits):
            raise ValueError('Sum requires values for lower and upper bounds.')

        return obj

    def _eval_is_zero(self):
        if self.function.is_zero:
            return True

    def doit(self, **hints):
        if hints.get('deep', True):
            f = self.function.doit(**hints)
        else:
            f = self.function

        for n, limit in enumerate(self.limits):
            i, a, b = limit
            dif = b - a
            if dif.is_integer and dif.is_negative:
                a, b = b + 1, a - 1
                f = -f

            newf = eval_sum(f, (i, a, b))
            if newf is None:
                if f == self.function:
                    return self
                else:
                    return self.func(f, *self.limits[n:])
            f = newf

        if hints.get('deep', True):
            # eval_sum could return partially unevaluated
            # result with Piecewise.  In this case we won't
            # doit() recursively.
            if not isinstance(f, Piecewise):
                return f.doit(**hints)

        return f

    def _eval_derivative(self, x):
        """
        Differentiate wrt x as long as x is not in the free symbols of any of
        the upper or lower limits.

        Sum(a*b*x, (x, 1, a)) can be differentiated wrt x or b but not `a`
        since the value of the sum is discontinuous in `a`. In a case
        involving a limit variable, the unevaluated derivative is returned.

        """
        # get limits and the function
        f, limits = self.function, list(self.limits)

        limit = limits.pop(-1)

        if limits:  # f is the argument to a Sum
            f = self.func(f, *limits)

        if len(limit) == 3:
            _, a, b = limit
            if x in a.free_symbols or x in b.free_symbols:
                return
            df = Derivative(f, x, evaluate=True)
            rv = self.func(df, limit)
            if limit[0] not in df.free_symbols:
                rv = rv.doit()
            return rv
        else:
            raise NotImplementedError('Lower and upper bound expected.')

    def _eval_simplify(self, ratio, measure):
        from ..simplify.simplify import sum_simplify
        return sum_simplify(self)

    def euler_maclaurin(self, m=0, n=0, eps=0, eval_integral=True):
        """
        Return an Euler-Maclaurin approximation of self, where m is the
        number of leading terms to sum directly and n is the number of
        terms in the tail.

        With m = n = 0, this is simply the corresponding integral
        plus a first-order endpoint correction.

        Returns (s, e) where s is the Euler-Maclaurin approximation
        and e is the estimated error (taken to be the magnitude of
        the first omitted term in the tail):

            >>> Sum(1/k, (k, 2, 5)).doit().evalf()
            1.28333333333333
            >>> s, e = Sum(1/k, (k, 2, 5)).euler_maclaurin()
            >>> s
            -log(2) + 7/20 + log(5)
            >>> print(sstr((s.evalf(), e.evalf()), full_prec=True))
            (1.26629073187415, 0.0175000000000000)

        The endpoints may be symbolic:

            >>> s, e = Sum(1/k, (k, a, b)).euler_maclaurin()
            >>> s
            -log(a) + log(b) + 1/(2*b) + 1/(2*a)
            >>> e
            Abs(1/(12*b**2) - 1/(12*a**2))

        If the function is a polynomial of degree at most 2n+1, the
        Euler-Maclaurin formula becomes exact (and e = 0 is returned):

            >>> Sum(k, (k, 2, b)).euler_maclaurin()
            (b**2/2 + b/2 - 1, 0)
            >>> Sum(k, (k, 2, b)).doit()
            b**2/2 + b/2 - 1

        With a nonzero `eps` specified, the summation is ended
        as soon as the remainder term is less than the epsilon.

        """
        from ..functions import bernoulli, factorial
        from ..integrals import Integral

        m = int(m)
        n = int(n)
        f = self.function
        if len(self.limits) != 1:
            raise ValueError('More than 1 limit')
        i, a, b = self.limits[0]
        if (a - b).is_positive:
            if a - b == 1:
                return Integer(0), Integer(0)
            a, b = b + 1, a - 1
            f = -f
        s = Integer(0)
        if m:
            if b.is_Integer and a.is_Integer:
                m = min(m, b - a + 1)
            if not eps or f.is_polynomial(i):
                for k in range(m):
                    s += f.subs({i: a + k})
            else:
                term = f.subs({i: a})
                if term:
                    test = abs(term.evalf(3)) < eps
                    if not (test == false):
                        # a symbolic Relational class, can't go further
                        return term, Integer(0)
                s += term
                for k in range(1, m):
                    term = f.subs({i: a + k})
                    if abs(term.evalf(3)) < eps and term != 0:
                        return s, abs(term)
                    s += term
            if b - a + 1 == m:
                return s, Integer(0)
            a += m
        x = Dummy('x')
        I = Integral(f.subs({i: x}), (x, a, b))
        if eval_integral:
            I = I.doit()
        s += I

        def fpoint(expr):
            if b is oo:
                return expr.subs({i: a}), 0
            return expr.subs({i: a}), expr.subs({i: b})
        fa, fb = fpoint(f)
        iterm = (fa + fb)/2
        g = f.diff(i)
        for k in range(1, n + 2):
            ga, gb = fpoint(g)
            term = bernoulli(2*k)/factorial(2*k)*(gb - ga)
            if eps and term and abs(term.evalf(3)) < eps:
                break
            if k <= n:
                s += term
                g = g.diff((i, 2), simplify=False)
        return s + iterm, abs(term)

    def reverse_order(self, *indices):
        r"""Reverse the order of a limit in a Sum.

        Parameters
        ==========

        \*indices : list
            The selectors in the argument ``indices`` specify some indices whose
            limits get reversed.  These selectors are either variable names or
            numerical indices counted starting from the inner-most limit tuple.

        Examples
        ========

        >>> Sum(x, (x, 0, 3)).reverse_order(x)
        Sum(-x, (x, 4, -1))
        >>> Sum(x*y, (x, 1, 5), (y, 0, 6)).reverse_order(x, y)
        Sum(x*y, (x, 6, 0), (y, 7, -1))
        >>> Sum(x, (x, a, b)).reverse_order(x)
        Sum(-x, (x, b + 1, a - 1))
        >>> Sum(x, (x, a, b)).reverse_order(0)
        Sum(-x, (x, b + 1, a - 1))

        While one should prefer variable names when specifying which limits
        to reverse, the index counting notation comes in handy in case there
        are several symbols with the same name.

        >>> s = Sum(x**2, (x, a, b), (x, c, d))
        >>> s
        Sum(x**2, (x, a, b), (x, c, d))
        >>> s0 = s.reverse_order(0)
        >>> s0
        Sum(-x**2, (x, b + 1, a - 1), (x, c, d))
        >>> s1 = s0.reverse_order(1)
        >>> s1
        Sum(x**2, (x, b + 1, a - 1), (x, d + 1, c - 1))

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

        return Sum(e * self.function, *limits)

    def findrecur(self, F=Function('F'), n=None):
        r"""Find a recurrence formula for the summand of the sum.

        Given a sum `f(n) = \sum_k F(n, k)`, where `F(n, k)` is
        doubly hypergeometric (that's, both `F(n + 1, k)/F(n, k)`
        and `F(n, k + 1)/F(n, k)` are rational functions of `n` and `k`),
        we find a recurrence for the summand `F(n, k)` of the form

            .. math:: \sum_{i=0}^I\sum_{j=0}^J a_{i,j}F(n - j, k - i) = 0

        Examples
        ========

        >>> s = Sum(factorial(n)/(factorial(k)*factorial(n - k)), (k, 0, oo))
        >>> s.findrecur()
        -F(n, k) + F(n - 1, k) + F(n - 1, k - 1)

        Notes
        =====

        We use Sister Celine's algorithm, see :cite:`Petkovsek1997AeqB`, Ch. 4.

        """
        from ..core import Mul, expand_func
        from ..functions import gamma
        from ..polys import factor, together
        from ..simplify import collect

        if len(self.variables) > 1:
            raise ValueError
        else:
            if self.limits[0][1:] != (0, oo):
                raise ValueError
            k = self.variables[0]

        if not n:
            try:
                n = (self.function.free_symbols - {k}).pop()
            except KeyError:
                raise ValueError

        a = Function('a')

        def f(i, j):
            return self.function.subs({n: i, k: j})

        I, J, step = 0, 1, 1
        y, x, sols = Integer(0), [], {}

        while not any(v for a, v in sols.items()):
            if step % 2 != 0:
                dy = sum(a(I, j)*f(n - j, k - I)/f(n, k) for j in range(J))
                dx = [a(I, j) for j in range(J)]
                I += 1
            else:
                dy = sum(a(i, J)*f(n - J, k - i)/f(n, k) for i in range(I))
                dx = [a(i, J) for i in range(I)]
                J += 1
            step += 1
            y += expand_func(dy.rewrite(gamma))
            x += dx

            t = together(y)
            numer = t.as_numer_denom()[0]
            numer = Mul(*[t for t in factor(numer).as_coeff_mul()[1] if t.has(a)])

            if not numer.is_rational_function(n, k):
                raise ValueError

            z = collect(numer, k)
            eq = z.as_poly(k).all_coeffs()
            sols = solve(eq, *x)[0]

        y = sum(a(i, j)*F(n - j, k - i) for i in range(I) for j in range(J))
        y = y.subs(sols).subs({_: 1 for _ in x})

        return y if y else None


def summation(f, *symbols, **kwargs):
    r"""Compute the summation of f with respect to symbols.

    The notation for symbols is similar to the notation used in Integral.
    summation(f, (i, a, b)) computes the sum of f with respect to i from a to b,
    i.e.,

    ::

                                    b
                                  ____
                                  \   `
        summation(f, (i, a, b)) =  )    f
                                  /___,
                                  i = a

    If it cannot compute the sum, it returns an unevaluated Sum object.
    Repeated sums can be computed by introducing additional symbols tuples::

    >>> i = symbols('i', integer=True)

    >>> summation(2*i - 1, (i, 1, n))
    n**2
    >>> summation(1/2**i, (i, 0, oo))
    2
    >>> summation(1/log(n)**n, (n, 2, oo))
    Sum(log(n)**(-n), (n, 2, oo))
    >>> summation(i, (i, 0, n), (n, 0, m))
    m**3/6 + m**2/2 + m/3

    >>> summation(x**n/factorial(n), (n, 0, oo))
    E**x

    See Also
    ========

    diofant.concrete.summations.Sum
    diofant.concrete.products.Product
    diofant.concrete.products.product

    """
    return Sum(f, *symbols, **kwargs).doit(deep=False)


def telescopic_direct(L, R, n, limits):
    """Returns the direct summation of the terms of a telescopic sum

    L is the term with lower index
    R is the term with higher index
    n difference between the indexes of L and R

    For example:

    >>> telescopic_direct(1/k, -1/(k+2), 2, (k, a, b))
    -1/(b + 2) - 1/(b + 1) + 1/(a + 1) + 1/a

    """
    (i, a, b) = limits
    s = 0
    for m in range(n):
        s += L.subs({i: a + m}) + R.subs({i: b - m})
    return s


def telescopic(L, R, limits):
    """Tries to perform the summation using the telescopic property

    return None if not possible

    """
    (i, a, b) = limits
    if L.is_Add or R.is_Add:
        return

    k = Wild('k')
    sol = (-R).match(L.subs({i: i + k}))
    if sol:
        s = sol[k]
    else:
        return

    if s.is_negative:
        return telescopic_direct(R, L, abs(s), (i, a, b))
    elif s.is_positive:
        return telescopic_direct(L, R, s, (i, a, b))


def eval_sum(f, limits):
    from ..functions import KroneckerDelta
    from .delta import _has_simple_delta, deltasummation

    (i, a, b) = limits
    if f == 0:
        return Integer(0)
    if i not in f.free_symbols:
        return f*(b - a + 1)
    if a == b:
        return f.subs({i: a})

    if f.has(KroneckerDelta) and _has_simple_delta(f, limits[0]):
        return deltasummation(f, limits)

    dif = b - a
    definite = dif.is_Integer
    # Doing it directly may be faster if there are very few terms.
    if definite and (dif < 100):
        return eval_sum_direct(f, (i, a, b))
    if isinstance(f, Piecewise):
        return
    # Try to do it symbolically. Even when the number of terms is known,
    # this can save time when b-a is big.
    # We should try to transform to partial fractions
    value = eval_sum_symbolic(f.expand(), (i, a, b))
    if value is not None:
        return value
    # Do it directly
    if definite:
        return eval_sum_direct(f, (i, a, b))


def eval_sum_direct(expr, limits):
    from ..core import Add
    (i, a, b) = limits

    dif = b - a
    return Add(*[expr.subs({i: a + j}) for j in range(dif + 1)])


def eval_sum_symbolic(f, limits):
    from ..functions import bernoulli, harmonic

    f_orig = f
    (i, a, b) = limits
    if not f.has(i):
        return f*(b - a + 1)

    # Linearity
    if f.is_Mul:
        L, R = f.as_two_terms()

        if not L.has(i):
            sR = eval_sum_symbolic(R, (i, a, b))
            if sR:
                return L*sR

        if not R.has(i):
            sL = eval_sum_symbolic(L, (i, a, b))
            if sL:
                return R*sL

        try:
            f = apart(f, i)  # see if it becomes an Add
        except PolynomialError:
            pass

    if f.is_Add:
        L, R = f.as_two_terms()
        lrsum = telescopic(L, R, (i, a, b))

        if lrsum:
            return lrsum

        lsum = eval_sum_symbolic(L, (i, a, b))
        rsum = eval_sum_symbolic(R, (i, a, b))

        if None not in (lsum, rsum):
            r = lsum + rsum
            if r is not nan:
                return r

    # Polynomial terms with Faulhaber's formula
    n = Wild('n')
    result = f.match(i**n)

    if result is not None:
        n = result[n]

        if n.is_Integer:
            if n >= 0:
                if (b is oo and a != -oo) or \
                   (a == -oo and b is not oo):
                    return oo
                return ((bernoulli(n + 1, b + 1) - bernoulli(n + 1, a))/(n + 1)).expand()
            elif a.is_Integer and a >= 1:
                if n == -1:
                    return harmonic(b) - harmonic(a - 1)
                else:
                    return harmonic(b, abs(n)) - harmonic(a - 1, abs(n))

    if not (a.has(oo, -oo) or
            b.has(oo, -oo)):
        # Geometric terms
        c1 = Wild('c1', exclude=[i])
        c2 = Wild('c2', exclude=[i])
        c3 = Wild('c3', exclude=[i])

        e = f.match(c1**(c2*i + c3))

        if e is not None:
            p = (c1**c3).subs(e)
            q = (c1**c2).subs(e)

            r = p*(q**a - q**(b + 1))/(1 - q)
            l = p*(b - a + 1)

            return Piecewise((l, Eq(q, 1)), (r, True))

    r = gosper_sum(f, (i, a, b))
    if r is not None and r.is_finite:
        return r

    return eval_sum_hyper(f_orig, (i, a, b))


def _eval_sum_hyper(f, i, a):
    """Returns (res, cond). Sums from a to oo."""
    from ..functions import hyper
    from ..polys import factor
    from ..simplify import fraction, hyperexpand, hypersimp, simplify

    if a != 0:
        return _eval_sum_hyper(f.subs({i: i + a}), i, 0)

    if f.subs({i: 0}) == 0:
        if simplify(f.subs({i: Dummy('i', integer=True, positive=True)})) == 0:
            return Integer(0), True
        return _eval_sum_hyper(f.subs({i: i + 1}), i, 0)

    hs = hypersimp(f, i)
    if hs is None:
        return

    numer, denom = fraction(factor(hs))
    top, topl = numer.as_coeff_mul(i)
    bot, botl = denom.as_coeff_mul(i)
    ab = [top, bot]
    factors = [topl, botl]
    params = [[], []]
    for k in range(2):
        for fac in factors[k]:
            mul = 1
            if fac.is_Pow:
                mul = fac.exp
                fac = fac.base
            p = fac.as_poly(i)
            if p.degree() != 1:
                return
            n, m = p.all_coeffs()
            ab[k] *= m**mul
            params[k] += [n/m]*mul

    # Add "1" to numerator parameters, to account for implicit n! in
    # hypergeometric series.
    ap = params[0] + [1]
    bq = params[1]
    x = ab[0]/ab[1]
    h = hyper(ap, bq, x)

    e = h
    try:
        e = hyperexpand(h)
    except PolynomialError:
        pass

    return f.limit(i, 0)*e, h.convergence_statement


def eval_sum_hyper(f, i_a_b):
    from ..logic import And

    i, a, b = i_a_b

    if (b - a).is_Integer:
        # We are never going to do better than doing the sum in the obvious way
        return

    old_sum = Sum(f, (i, a, b))

    if b != oo:
        if a == -oo:
            res = _eval_sum_hyper(f.subs({i: -i}), i, -b)
            if res is not None:
                return Piecewise(res, (old_sum, True))
        else:
            res1 = _eval_sum_hyper(f, i, a)
            res2 = _eval_sum_hyper(f, i, b + 1)
            if res1 is None or res2 is None:
                return
            (res1, cond1), (res2, cond2) = res1, res2
            cond = And(cond1, cond2)
            if cond == false:
                return
            return Piecewise((res1 - res2, cond), (old_sum, True))

    if a != -oo:
        # Now b == oo, a != -oo
        res = _eval_sum_hyper(f, i, a)
        if res is not None:
            r, c = res
            if c == false:
                f = f.subs({i: Dummy('i', integer=True, positive=True) + a})
                if f.is_nonnegative:
                    return oo
                else:
                    return
            return Piecewise(res, (old_sum, True))
