"""Hypergeometric and Meijer G-functions."""

import functools
import math

import mpmath

from ...core import (Derivative, Dummy, Expr, Function, I, Integer, Mod, Mul,
                     Ne, Rational, Tuple, oo, pi, zoo)
from ...core.function import ArgumentIndexError
from .. import (acosh, acoth, asin, asinh, atan, atanh, cos, cosh, exp, log,
                sin, sinh, sqrt)


class TupleArg(Tuple):
    """Arguments of the hyper/meijerg functions."""

    def limit(self, x, xlim, dir='+'):
        """Compute limit x->xlim."""
        return self.func(*[_.limit(x, xlim, dir) for _ in self.args])


# TODO should __new__ accept **options?
# TODO should constructors should check if parameters are sensible?


def _prep_tuple(v):
    """
    Turn an iterable argument V into a Tuple and unpolarify, since both
    hypergeometric and meijer g-functions are unbranched in their parameters.

    Examples
    ========

    >>> _prep_tuple([1, 2, 3])
    (1, 2, 3)
    >>> _prep_tuple((4, 5))
    (4, 5)
    >>> _prep_tuple((7, 8, 9))
    (7, 8, 9)

    """
    from .. import unpolarify
    return TupleArg(*[unpolarify(x) for x in v])


class TupleParametersBase(Function):
    """Base class that takes care of differentiation, when some of
    the arguments are actually tuples.

    """

    # This is not deduced automatically since there are Tuples as arguments.
    is_commutative = True

    def _eval_derivative(self, s):
        try:
            res = 0
            if self.args[0].has(s) or self.args[1].has(s):
                for i, p in enumerate(self._diffargs):
                    m = self._diffargs[i].diff(s)
                    if m != 0:
                        res += self.fdiff((1, i))*m
            return res + self.fdiff(3)*self.args[2].diff(s)
        except (ArgumentIndexError, NotImplementedError):
            return Derivative(self, s)

    @property
    def is_number(self):
        """Returns True if 'self' has no free symbols."""
        return not self.free_symbols


class hyper(TupleParametersBase):
    r"""
    The (generalized) hypergeometric function is defined by a series where
    the ratios of successive terms are a rational function of the summation
    index. When convergent, it is continued analytically to the largest
    possible domain.

    The hypergeometric function depends on two vectors of parameters, called
    the numerator parameters `a_p`, and the denominator parameters
    `b_q`. It also has an argument `z`. The series definition is

    .. math ::
        {}_pF_q\left(\begin{matrix} a_1, \ldots, a_p \\ b_1, \ldots, b_q \end{matrix}
                     \middle| z \right)
        = \sum_{n=0}^\infty \frac{(a_1)_n \ldots (a_p)_n}{(b_1)_n \ldots (b_q)_n}
                            \frac{z^n}{n!},

    where `(a)_n = (a)(a+1)\ldots(a+n-1)` denotes the rising factorial.

    If one of the `b_q` is a non-positive integer then the series is
    undefined unless one of the `a_p` is a larger (i.e. smaller in
    magnitude) non-positive integer. If none of the `b_q` is a
    non-positive integer and one of the `a_p` is a non-positive
    integer, then the series reduces to a polynomial. To simplify the
    following discussion, we assume that none of the `a_p` or
    `b_q` is a non-positive integer. For more details, see the
    references.

    The series converges for all `z` if `p \le q`, and thus
    defines an entire single-valued function in this case. If `p =
    q+1` the series converges for `|z| < 1`, and can be continued
    analytically into a half-plane. If `p > q+1` the series is
    divergent for all `z`.

    Note: The hypergeometric function constructor currently does *not* check
    if the parameters actually yield a well-defined function.

    Examples
    ========

    The parameters `a_p` and `b_q` can be passed as arbitrary
    iterables, for example:

    >>> hyper((1, 2, 3), [3, 4], x)
    hyper((1, 2, 3), (3, 4), x)

    There is also pretty printing (it looks better using unicode):

    >>> pprint(hyper((1, 2, 3), [3, 4], x), use_unicode=False)
      _
     |_  /1, 2, 3 |  \
     |   |        | x|
    3  2 \  3, 4  |  /

    The parameters must always be iterables, even if they are vectors of
    length one or zero:

    >>> hyper([1], [], x)
    hyper((1,), (), x)

    But of course they may be variables (but if they depend on x then you
    should not expect much implemented functionality):

    >>> hyper([n, a], [n**2], x)
    hyper((n, a), (n**2,), x)

    The hypergeometric function generalizes many named special functions.
    The function hyperexpand() tries to express a hypergeometric function
    using named special functions.
    For example:

    >>> hyperexpand(hyper([], [], x))
    E**x

    You can also use expand_func:

    >>> expand_func(x*hyper([1, 1], [2], -x))
    log(x + 1)

    More examples:

    >>> hyperexpand(hyper([], [Rational(1, 2)], -x**2/4))
    cos(x)
    >>> hyperexpand(x*hyper([Rational(1, 2), Rational(1, 2)], [Rational(3, 2)], x**2))
    asin(x)

    We can also sometimes hyperexpand parametric functions:

    >>> hyperexpand(hyper([-a], [], x))
    (-x + 1)**a

    See Also
    ========

    diofant.simplify.hyperexpand
    diofant.functions.special.gamma_functions.gamma
    diofant.functions.special.hyper.meijerg

    References
    ==========

    * Luke, Y. L. (1969), The Special Functions and Their Approximations,
      Volume 1
    * https://en.wikipedia.org/wiki/Generalized_hypergeometric_function

    """

    def __new__(cls, ap, bq, z):
        # TODO should we check convergence conditions?
        return Function.__new__(cls, _prep_tuple(ap), _prep_tuple(bq), z)

    @classmethod
    def eval(cls, ap, bq, z):
        from .. import unpolarify
        if len(ap) <= len(bq):
            nz = unpolarify(z)
            if z != nz:
                return hyper(ap, bq, nz)

    def fdiff(self, argindex=3):
        if argindex != 3:
            raise ArgumentIndexError(self, argindex)
        nap = Tuple(*[a + 1 for a in self.ap])
        nbq = Tuple(*[b + 1 for b in self.bq])
        fac = Mul(*self.ap)/Mul(*self.bq)
        return fac*hyper(nap, nbq, self.argument)

    def _eval_expand_func(self, **hints):
        from ...simplify import hyperexpand
        from .gamma_functions import gamma
        if len(self.ap) == 2 and len(self.bq) == 1 and self.argument == 1:
            a, b = self.ap
            c = self.bq[0]
            return gamma(c)*gamma(c - a - b)/gamma(c - a)/gamma(c - b)
        return hyperexpand(self)

    def _eval_rewrite_as_Sum(self, ap, bq, z):
        from ...concrete import Sum
        from .. import Piecewise, RisingFactorial, factorial
        n = Dummy('n', integer=True)
        rfap = Tuple(*[RisingFactorial(a, n) for a in ap])
        rfbq = Tuple(*[RisingFactorial(b, n) for b in bq])
        coeff = Mul(*rfap) / Mul(*rfbq)
        return Piecewise((Sum(coeff * z**n / factorial(n), (n, 0, oo)),
                          self.convergence_statement), (self, True))

    @property
    def argument(self):
        """Argument of the hypergeometric function."""
        return self.args[2]

    @property
    def ap(self):
        """Numerator parameters of the hypergeometric function."""
        return Tuple(*self.args[0])

    @property
    def bq(self):
        """Denominator parameters of the hypergeometric function."""
        return Tuple(*self.args[1])

    @property
    def _diffargs(self):
        return self.ap + self.bq

    @property
    def eta(self):
        """A quantity related to the convergence of the series."""
        return sum(self.ap) - sum(self.bq)

    @property
    def radius_of_convergence(self):
        """
        Compute the radius of convergence of the defining series.

        Note that even if this is not oo, the function may still be evaluated
        outside of the radius of convergence by analytic continuation. But if
        this is zero, then the function is not actually defined anywhere else.

        >>> hyper((1, 2), [3], z).radius_of_convergence
        1
        >>> hyper((1, 2, 3), [4], z).radius_of_convergence
        0
        >>> hyper((1, 2), (3, 4), z).radius_of_convergence
        oo

        """
        if any(a.is_integer and a.is_nonpositive for a in self.ap + self.bq):
            aints = [a for a in self.ap if a.is_Integer and a.is_nonpositive]
            bints = [a for a in self.bq if a.is_Integer and a.is_nonpositive]
            if len(aints) < len(bints):
                return Integer(0)
            popped = False
            for b in bints:
                cancelled = False
                while aints:
                    a = aints.pop()
                    if a >= b:
                        cancelled = True
                        break
                    popped = True
                if not cancelled:
                    return Integer(0)
            if aints or popped:
                # There are still non-positive numerator parameters.
                # This is a polynomial.
                return oo
        if len(self.ap) == len(self.bq) + 1:
            return Integer(1)
        elif len(self.ap) <= len(self.bq):
            return oo
        else:
            return Integer(0)

    @property
    def convergence_statement(self):
        """Return a condition on z under which the series converges."""
        from ...logic import And, Or
        from .. import re
        R = self.radius_of_convergence
        if R == 0:
            return False
        if R == oo:
            return True
        # The special functions and their approximations, page 44
        e = self.eta
        z = self.argument
        c1 = And(re(e) < 0, abs(z) <= 1)
        c2 = And(0 <= re(e), re(e) < 1, abs(z) <= 1, Ne(z, 1))
        c3 = And(re(e) >= 1, abs(z) < 1)
        return Or(c1, c2, c3)

    def _eval_simplify(self, ratio, measure):
        from ...simplify import hyperexpand
        return hyperexpand(self)

    def _eval_evalf(self, prec):
        z = self.argument._to_mpmath(prec)
        ap = [a._to_mpmath(prec) for a in self.ap]
        bp = [b._to_mpmath(prec) for b in self.bq]
        with mpmath.workprec(prec):
            res = mpmath.hyper(ap, bp, z, eliminate=False)
        return Expr._from_mpmath(res, prec)


class meijerg(TupleParametersBase):
    r"""
    The Meijer G-function is defined by a Mellin-Barnes type integral that
    resembles an inverse Mellin transform. It generalizes the hypergeometric
    functions.

    The Meijer G-function depends on four sets of parameters. There are
    "*numerator parameters*"
    `a_1, \ldots, a_n` and `a_{n+1}, \ldots, a_p`, and there are
    "*denominator parameters*"
    `b_1, \ldots, b_m` and `b_{m+1}, \ldots, b_q`.
    Confusingly, it is traditionally denoted as follows (note the position
    of `m`, `n`, `p`, `q`, and how they relate to the lengths of the four
    parameter vectors):

    .. math ::
        G_{p,q}^{m,n} \left(\begin{matrix}a_1, \ldots, a_n & a_{n+1}, \ldots, a_p \\
                                        b_1, \ldots, b_m & b_{m+1}, \ldots, b_q
                          \end{matrix} \middle| z \right).

    However, in diofant the four parameter vectors are always available
    separately (see examples), so that there is no need to keep track of the
    decorating sub- and super-scripts on the G symbol.

    The G function is defined as the following integral:

    .. math ::
         \frac{1}{2 \pi i} \int_L \frac{\prod_{j=1}^m \Gamma(b_j - s)
         \prod_{j=1}^n \Gamma(1 - a_j + s)}{\prod_{j=m+1}^q \Gamma(1- b_j +s)
         \prod_{j=n+1}^p \Gamma(a_j - s)} z^s \mathrm{d}s,

    where `\Gamma(z)` is the gamma function. There are three possible
    contours which we will not describe in detail here (see the references).
    If the integral converges along more than one of them the definitions
    agree. The contours all separate the poles of `\Gamma(1-a_j+s)`
    from the poles of `\Gamma(b_k-s)`, so in particular the G function
    is undefined if `a_j - b_k \in \mathbb{Z}_{>0}` for some
    `j \le n` and `k \le m`.

    The conditions under which one of the contours yields a convergent integral
    are complicated and we do not state them here, see the references.

    Note: Currently the Meijer G-function constructor does *not* check any
    convergence conditions.

    Examples
    ========

    You can pass the parameters either as four separate vectors:

    >>> pprint(meijerg([1, 2], [a, 4], [5], [], x), use_unicode=False)
     __1, 2 /1, 2  a, 4 |  \
    /__     |           | x|
    \_|4, 1 \ 5         |  /

    or as two nested vectors:

    >>> pprint(meijerg(([1, 2], [3, 4]), ([5], []), x), use_unicode=False)
     __1, 2 /1, 2  3, 4 |  \
    /__     |           | x|
    \_|4, 1 \ 5         |  /

    As with the hypergeometric function, the parameters may be passed as
    arbitrary iterables. Vectors of length zero and one also have to be
    passed as iterables. The parameters need not be constants, but if they
    depend on the argument then not much implemented functionality should be
    expected.

    All the subvectors of parameters are available:

    >>> g = meijerg([1], [2], [3], [4], x)
    >>> pprint(g, use_unicode=False)
     __1, 1 /1  2 |  \
    /__     |     | x|
    \_|2, 2 \3  4 |  /
    >>> g.an
    (1,)
    >>> g.ap
    (1, 2)
    >>> g.aother
    (2,)
    >>> g.bm
    (3,)
    >>> g.bq
    (3, 4)
    >>> g.bother
    (4,)

    The Meijer G-function generalizes the hypergeometric functions.
    In some cases it can be expressed in terms of hypergeometric functions,
    using Slater's theorem. For example:

    >>> hyperexpand(meijerg([a], [], [c], [b], x), allow_hyper=True)
    x**c*gamma(-a + c + 1)*hyper((-a + c + 1,),
                                 (-b + c + 1,), -x)/gamma(-b + c + 1)

    Thus the Meijer G-function also subsumes many named functions as special
    cases. You can use expand_func or hyperexpand to (try to) rewrite a
    Meijer G-function in terms of named special functions. For example:

    >>> expand_func(meijerg([[], []], [[0], []], -x))
    E**x
    >>> hyperexpand(meijerg([[], []], [[Rational(1, 2)], [0]], (x/2)**2))
    sin(x)/sqrt(pi)

    See Also
    ========

    diofant.functions.special.hyper.hyper
    diofant.simplify.hyperexpand

    References
    ==========

    * Luke, Y. L. (1969), The Special Functions and Their Approximations,
      Volume 1
    * https://en.wikipedia.org/wiki/Meijer_G-function

    """

    def __new__(cls, *args):
        if len(args) == 5:
            args = [(args[0], args[1]), (args[2], args[3]), args[4]]
        if len(args) != 3:
            raise TypeError("args must be either as, as', bs, bs', z or "
                            'as, bs, z')

        def tr(p):
            if len(p) != 2:
                raise TypeError('wrong argument')
            return TupleArg(_prep_tuple(p[0]), _prep_tuple(p[1]))

        arg0, arg1 = tr(args[0]), tr(args[1])
        if Tuple(arg0, arg1).has(oo, zoo,
                                 -oo):
            raise ValueError('G-function parameters must be finite')

        if any((a - b).is_integer and (a - b).is_positive
               for a in arg0[0] for b in arg1[0]):
            raise ValueError('no parameter a1, ..., an may differ from '
                             'any b1, ..., bm by a positive integer')

        # TODO should we check convergence conditions?
        return Function.__new__(cls, arg0, arg1, args[2])

    def fdiff(self, argindex=3):
        if argindex != 3:
            return self._diff_wrt_parameter(argindex[1])
        if len(self.an) >= 1:
            a = list(self.an)
            a[0] -= 1
            G = meijerg(a, self.aother, self.bm, self.bother, self.argument)
            return 1/self.argument * ((self.an[0] - 1)*self + G)
        elif len(self.bm) >= 1:
            b = list(self.bm)
            b[0] += 1
            G = meijerg(self.an, self.aother, b, self.bother, self.argument)
            return 1/self.argument * (self.bm[0]*self - G)
        else:
            return Integer(0)

    def _diff_wrt_parameter(self, idx):
        # Differentiation wrt a parameter can only be done in very special
        # cases. In particular, if we want to differentiate with respect to
        # `a`, all other gamma factors have to reduce to rational functions.
        #
        # Let MT denote mellin transform. Suppose T(-s) is the gamma factor
        # appearing in the definition of G. Then
        #
        #   MT(log(z)G(z)) = d/ds T(s) = d/da T(s) + ...
        #
        # Thus d/da G(z) = log(z)G(z) - ...
        # The ... can be evaluated as a G function under the above conditions,
        # the formula being most easily derived by using
        #
        # d  Gamma(s + n)    Gamma(s + n) / 1    1                1     \
        # -- ------------ =  ------------ | - + ----  + ... + --------- |
        # ds Gamma(s)        Gamma(s)     \ s   s + 1         s + n - 1 /
        #
        # which follows from the difference equation of the digamma function.
        # (There is a similar equation for -n instead of +n).

        # We first figure out how to pair the parameters.
        an = list(self.an)
        ap = list(self.aother)
        bm = list(self.bm)
        bq = list(self.bother)
        if idx < len(an):
            an.pop(idx)
        else:
            idx -= len(an)
            if idx < len(ap):
                ap.pop(idx)
            else:
                idx -= len(ap)
                if idx < len(bm):
                    bm.pop(idx)
                else:
                    bq.pop(idx - len(bm))
        pairs1 = []
        pairs2 = []
        for l1, l2, pairs in [(an, bq, pairs1), (ap, bm, pairs2)]:
            while l1:
                x = l1.pop()
                found = None
                for i, y in enumerate(l2):
                    if not Mod((x - y).simplify(), 1):
                        found = i
                        break
                if found is None:
                    raise NotImplementedError('Derivative not expressible '
                                              'as G-function?')
                y = l2[i]
                l2.pop(i)
                pairs.append((x, y))

        # Now build the result.
        res = log(self.argument)*self

        for a, b in pairs1:
            sign = 1
            n = a - b
            base = b
            if n < 0:
                sign = -1
                n = b - a
                base = a
            for k in range(n):
                res -= sign*meijerg(self.an + (base + k + 1,), self.aother,
                                    self.bm, self.bother + (base + k + 0,),
                                    self.argument)

        for a, b in pairs2:
            sign = 1
            n = b - a
            base = a
            if n < 0:
                sign = -1
                n = a - b
                base = b
            for k in range(n):
                res -= sign*meijerg(self.an, self.aother + (base + k + 1,),
                                    self.bm + (base + k + 0,), self.bother,
                                    self.argument)

        return res

    def get_period(self):
        """
        Return a number P such that G(x*exp(I*P)) == G(x).

        >>> meijerg([1], [], [], [], z).get_period()
        2*pi
        >>> meijerg([pi], [], [], [], z).get_period()
        oo
        >>> meijerg([1, 2], [], [], [], z).get_period()
        oo
        >>> meijerg([1, 1], [2], [1, Rational(1, 2), Rational(1, 3)], [1], z).get_period()
        12*pi

        """
        # This follows from slater's theorem.
        def compute(l):
            # first check that no two differ by an integer
            for i, b in enumerate(l):
                if not b.is_Rational:
                    return oo
                for j in range(i + 1, len(l)):
                    if not Mod((b - l[j]).simplify(), 1):
                        return oo
            return functools.reduce(math.lcm, (x.denominator for x in l), 1)
        beta = compute(self.bm)
        alpha = compute(self.an)
        p, q = len(self.ap), len(self.bq)
        if p == q:
            if beta == oo or alpha == oo:
                return oo
            return 2*pi*math.lcm(alpha, beta)
        elif p < q:
            return 2*pi*beta
        else:
            return 2*pi*alpha

    def _eval_expand_func(self, **hints):
        from ...simplify import hyperexpand
        return hyperexpand(self)

    def _eval_evalf(self, prec):
        # The default code is insufficient for polar arguments.
        # mpmath provides an optional argument "r", which evaluates
        # G(z**(1/r)). I am not sure what its intended use is, but we hijack it
        # here in the following way: to evaluate at a number z of |argument|
        # less than (say) n*pi, we put r=1/n, compute z' = root(z, n)
        # (carefully so as not to loose the branch information), and evaluate
        # G(z'**(1/r)) = G(z'**n) = G(z).
        from .. import ceiling, exp_polar
        z = self.argument
        znum = self.argument.evalf(prec, strict=False)
        if znum.has(exp_polar):
            znum, branch = znum.as_coeff_mul(exp_polar)
            if len(branch) != 1:
                return
            branch = branch[0].args[0]/I
        else:
            branch = Integer(0)
        n = ceiling(abs(branch/pi)) + 1
        znum = znum**(Integer(1)/n)*exp(I*branch / n)

        # Convert all args to mpf or mpc
        [z, r, ap, bq] = [arg._to_mpmath(prec)
                          for arg in [znum, 1/n, self.args[0], self.args[1]]]

        with mpmath.workprec(prec):
            v = mpmath.meijerg(ap, bq, z, r)

        return Expr._from_mpmath(v, prec)

    def integrand(self, s):
        """Get the defining integrand D(s)."""
        from .gamma_functions import gamma
        return self.argument**s \
            * Mul(*(gamma(b - s) for b in self.bm)) \
            * Mul(*(gamma(1 - a + s) for a in self.an)) \
            / Mul(*(gamma(1 - b + s) for b in self.bother)) \
            / Mul(*(gamma(a - s) for a in self.aother))

    @property
    def argument(self):
        """Argument of the Meijer G-function."""
        return self.args[2]

    @property
    def an(self):
        """First set of numerator parameters."""
        return Tuple(*self.args[0][0])

    @property
    def ap(self):
        """Combined numerator parameters."""
        return Tuple(*(self.args[0][0] + self.args[0][1]))

    @property
    def aother(self):
        """Second set of numerator parameters."""
        return Tuple(*self.args[0][1])

    @property
    def bm(self):
        """First set of denominator parameters."""
        return Tuple(*self.args[1][0])

    @property
    def bq(self):
        """Combined denominator parameters."""
        return Tuple(*(self.args[1][0] + self.args[1][1]))

    @property
    def bother(self):
        """Second set of denominator parameters."""
        return Tuple(*self.args[1][1])

    @property
    def _diffargs(self):
        return self.ap + self.bq

    @property
    def nu(self):
        """A quantity related to the convergence region of the integral,
        c.f. references.

        """
        return sum(self.bq) - sum(self.ap)

    @property
    def delta(self):
        """A quantity related to the convergence region of the integral,
        c.f. references.

        """
        return len(self.bm) + len(self.an) - Integer(len(self.ap) + len(self.bq))/2


class HyperRep(Function):
    """
    A base class for "hyper representation functions".

    This is used exclusively in hyperexpand(), but fits more logically here.

    pFq is branched at 1 if p == q+1. For use with slater-expansion, we want
    define an "analytic continuation" to all polar numbers, which is
    continuous on circles and on the ray t*exp_polar(I*pi). Moreover, we want
    a "nice" expression for the various cases.

    This base class contains the core logic, concrete derived classes only
    supply the actual functions.

    """

    @classmethod
    def eval(cls, *args):
        from .. import unpolarify
        newargs = tuple(map(unpolarify, args[:-1])) + args[-1:]
        if args != newargs:
            return cls(*newargs)

    @classmethod
    def _expr_small(cls, x):
        """An expression for F(x) which holds for |x| < 1."""
        raise NotImplementedError

    @classmethod
    def _expr_small_minus(cls, x):
        """An expression for F(-x) which holds for |x| < 1."""
        raise NotImplementedError

    @classmethod
    def _expr_big(cls, x, n):
        """An expression for F(exp_polar(2*I*pi*n)*x), |x| > 1."""
        raise NotImplementedError

    @classmethod
    def _expr_big_minus(cls, x, n):
        """An expression for F(exp_polar(2*I*pi*n + pi*I)*x), |x| > 1."""
        raise NotImplementedError

    def _eval_rewrite_as_nonrep(self, *args):
        from .. import Piecewise
        x, n = self.args[-1].extract_branch_factor(allow_half=True)
        minus = False
        newargs = self.args[:-1] + (x,)
        if not n.is_Integer:
            minus = True
            n -= Rational(1, 2)
        newerargs = newargs + (n,)
        if minus:
            small = self._expr_small_minus(*newargs)
            big = self._expr_big_minus(*newerargs)
        else:
            small = self._expr_small(*newargs)
            big = self._expr_big(*newerargs)

        if big == small:
            return small
        return Piecewise((big, abs(x) > 1), (small, True))

    def _eval_rewrite_as_nonrepsmall(self, *args):
        x, n = self.args[-1].extract_branch_factor(allow_half=True)
        args = self.args[:-1] + (x,)
        if not n.is_Integer:
            return self._expr_small_minus(*args)
        return self._expr_small(*args)


class HyperRep_power1(HyperRep):
    """Return a representative for hyper([-a], [], z) == (1 - z)**a."""

    @classmethod
    def _expr_small(cls, a, x):
        return (1 - x)**a

    @classmethod
    def _expr_small_minus(cls, a, x):
        return (1 + x)**a

    @classmethod
    def _expr_big(cls, a, x, n):
        if a.is_integer:
            return cls._expr_small(a, x)
        return (x - 1)**a*exp((2*n - 1)*pi*I*a)

    @classmethod
    def _expr_big_minus(cls, a, x, n):
        if a.is_integer:
            return cls._expr_small_minus(a, x)
        return (1 + x)**a*exp(2*n*pi*I*a)


class HyperRep_power2(HyperRep):
    """Return a representative for hyper([a, a - 1/2], [2*a], z)."""

    @classmethod
    def _expr_small(cls, a, x):
        return 2**(2*a - 1)*(1 + sqrt(1 - x))**(1 - 2*a)

    @classmethod
    def _expr_small_minus(cls, a, x):
        return 2**(2*a - 1)*(1 + sqrt(1 + x))**(1 - 2*a)

    @classmethod
    def _expr_big(cls, a, x, n):
        sgn = -1
        if n.is_odd:
            sgn = 1
            n -= 1
        return 2**(2*a - 1)*(1 + sgn*I*sqrt(x - 1))**(1 - 2*a) \
            * exp(-2*n*pi*I*a)

    @classmethod
    def _expr_big_minus(cls, a, x, n):
        sgn = 1
        if n.is_odd:
            sgn = -1
        return sgn*2**(2*a - 1)*(sqrt(1 + x) + sgn)**(1 - 2*a)*exp(-2*pi*I*a*n)


class HyperRep_log1(HyperRep):
    """Represent -z*hyper([1, 1], [2], z) == log(1 - z)."""

    @classmethod
    def _expr_small(cls, x):
        return log(1 - x)

    @classmethod
    def _expr_small_minus(cls, x):
        return log(1 + x)

    @classmethod
    def _expr_big(cls, x, n):
        return log(x - 1) + (2*n - 1)*pi*I

    @classmethod
    def _expr_big_minus(cls, x, n):
        return log(1 + x) + 2*n*pi*I


class HyperRep_atanh(HyperRep):
    """Represent hyper([1/2, 1], [3/2], z) == atanh(sqrt(z))/sqrt(z)."""

    @classmethod
    def _expr_small(cls, x):
        return atanh(sqrt(x))/sqrt(x)

    def _expr_small_minus(self, x):
        return atan(sqrt(x))/sqrt(x)

    def _expr_big(self, x, n):
        if n.is_even:
            return (acoth(sqrt(x)) + I*pi/2)/sqrt(x)
        else:
            return (acoth(sqrt(x)) - I*pi/2)/sqrt(x)

    def _expr_big_minus(self, x, n):
        if n.is_even:
            return atan(sqrt(x))/sqrt(x)
        else:
            return (atan(sqrt(x)) - pi)/sqrt(x)


class HyperRep_asin1(HyperRep):
    """Represent hyper([1/2, 1/2], [3/2], z) == asin(sqrt(z))/sqrt(z)."""

    @classmethod
    def _expr_small(cls, z):
        return asin(sqrt(z))/sqrt(z)

    @classmethod
    def _expr_small_minus(cls, z):
        return asinh(sqrt(z))/sqrt(z)

    @classmethod
    def _expr_big(cls, z, n):
        return Integer(-1)**n*((Rational(1, 2) - n)*pi/sqrt(z) + I*acosh(sqrt(z))/sqrt(z))

    @classmethod
    def _expr_big_minus(cls, z, n):
        return Integer(-1)**n*(asinh(sqrt(z))/sqrt(z) + n*pi*I/sqrt(z))


class HyperRep_asin2(HyperRep):
    """Represent hyper([1, 1], [3/2], z) == asin(sqrt(z))/sqrt(z)/sqrt(1-z)."""

    # TODO this can be nicer
    @classmethod
    def _expr_small(cls, z):
        return HyperRep_asin1._expr_small(z) \
            / HyperRep_power1._expr_small(Rational(1, 2), z)

    @classmethod
    def _expr_small_minus(cls, z):
        return HyperRep_asin1._expr_small_minus(z) \
            / HyperRep_power1._expr_small_minus(Rational(1, 2), z)

    @classmethod
    def _expr_big(cls, z, n):
        return HyperRep_asin1._expr_big(z, n) \
            / HyperRep_power1._expr_big(Rational(1, 2), z, n)

    @classmethod
    def _expr_big_minus(cls, z, n):
        return HyperRep_asin1._expr_big_minus(z, n) \
            / HyperRep_power1._expr_big_minus(Rational(1, 2), z, n)


class HyperRep_sqrts1(HyperRep):
    """Return a representative for hyper([-a, 1/2 - a], [1/2], z)."""

    @classmethod
    def _expr_small(cls, a, z):
        return ((1 - sqrt(z))**(2*a) + (1 + sqrt(z))**(2*a))/2

    @classmethod
    def _expr_small_minus(cls, a, z):
        return (1 + z)**a*cos(2*a*atan(sqrt(z)))

    @classmethod
    def _expr_big(cls, a, z, n):
        if n.is_even:
            return ((sqrt(z) + 1)**(2*a)*exp(2*pi*I*n*a) +
                    (sqrt(z) - 1)**(2*a)*exp(2*pi*I*(n - 1)*a))/2
        else:
            n -= 1
            return ((sqrt(z) - 1)**(2*a)*exp(2*pi*I*a*(n + 1)) +
                    (sqrt(z) + 1)**(2*a)*exp(2*pi*I*a*n))/2

    @classmethod
    def _expr_big_minus(cls, a, z, n):
        if n.is_even:
            return (1 + z)**a*exp(2*pi*I*n*a)*cos(2*a*atan(sqrt(z)))
        else:
            return (1 + z)**a*exp(2*pi*I*n*a)*cos(2*a*atan(sqrt(z)) - 2*pi*a)


class HyperRep_sqrts2(HyperRep):
    """Return a representative for
    sqrt(z)/2*[(1-sqrt(z))**2a - (1 + sqrt(z))**2a]
    == -2*z/(2*a+1) d/dz hyper([-a - 1/2, -a], [1/2], z)

    """

    @classmethod
    def _expr_small(cls, a, z):
        return sqrt(z)*((1 - sqrt(z))**(2*a) - (1 + sqrt(z))**(2*a))/2

    @classmethod
    def _expr_small_minus(cls, a, z):
        return sqrt(z)*(1 + z)**a*sin(2*a*atan(sqrt(z)))

    @classmethod
    def _expr_big(cls, a, z, n):
        if n.is_even:
            return sqrt(z)/2*((sqrt(z) - 1)**(2*a)*exp(2*pi*I*a*(n - 1)) -
                              (sqrt(z) + 1)**(2*a)*exp(2*pi*I*a*n))
        else:
            n -= 1
            return sqrt(z)/2*((sqrt(z) - 1)**(2*a)*exp(2*pi*I*a*(n + 1)) -
                              (sqrt(z) + 1)**(2*a)*exp(2*pi*I*a*n))

    def _expr_big_minus(self, a, z, n):
        if n.is_even:
            return (1 + z)**a*exp(2*pi*I*n*a)*sqrt(z)*sin(2*a*atan(sqrt(z)))
        else:
            return (1 + z)**a*exp(2*pi*I*n*a)*sqrt(z) \
                * sin(2*a*atan(sqrt(z)) - 2*pi*a)


class HyperRep_log2(HyperRep):
    """Represent log(1/2 + sqrt(1 - z)/2) == -z/4*hyper([3/2, 1, 1], [2, 2], z)."""

    @classmethod
    def _expr_small(cls, z):
        return log(Rational(1, 2) + sqrt(1 - z)/2)

    @classmethod
    def _expr_small_minus(cls, z):
        return log(Rational(1, 2) + sqrt(1 + z)/2)

    @classmethod
    def _expr_big(cls, z, n):
        if n.is_even:
            return (n - Rational(1, 2))*pi*I + log(sqrt(z)/2) + I*asin(1/sqrt(z))
        else:
            return (n - Rational(1, 2))*pi*I + log(sqrt(z)/2) - I*asin(1/sqrt(z))

    def _expr_big_minus(self, z, n):
        if n.is_even:
            return pi*I*n + log(Rational(1, 2) + sqrt(1 + z)/2)
        else:
            return pi*I*n + log(sqrt(1 + z)/2 - Rational(1, 2))


class HyperRep_cosasin(HyperRep):
    """Represent hyper([a, -a], [1/2], z) == cos(2*a*asin(sqrt(z)))."""

    # Note there are many alternative expressions, e.g. as powers of a sum of
    # square roots.

    @classmethod
    def _expr_small(cls, a, z):
        return cos(2*a*asin(sqrt(z)))

    @classmethod
    def _expr_small_minus(cls, a, z):
        return cosh(2*a*asinh(sqrt(z)))

    @classmethod
    def _expr_big(cls, a, z, n):
        return cosh(2*a*acosh(sqrt(z)) + a*pi*I*(2*n - 1))

    @classmethod
    def _expr_big_minus(cls, a, z, n):
        return cosh(2*a*asinh(sqrt(z)) + 2*a*pi*I*n)


class HyperRep_sinasin(HyperRep):
    """Represent 2*a*z*hyper([1 - a, 1 + a], [3/2], z)
    == sqrt(z)/sqrt(1-z)*sin(2*a*asin(sqrt(z)))

    """

    @classmethod
    def _expr_small(cls, a, z):
        return sqrt(z)/sqrt(1 - z)*sin(2*a*asin(sqrt(z)))

    @classmethod
    def _expr_small_minus(cls, a, z):
        return -sqrt(z)/sqrt(1 + z)*sinh(2*a*asinh(sqrt(z)))

    @classmethod
    def _expr_big(cls, a, z, n):
        return -1/sqrt(1 - 1/z)*sinh(2*a*acosh(sqrt(z)) + a*pi*I*(2*n - 1))

    @classmethod
    def _expr_big_minus(cls, a, z, n):
        return -1/sqrt(1 + 1/z)*sinh(2*a*asinh(sqrt(z)) + 2*a*pi*I*n)
