"""Riemann zeta and related function."""

from ...core import (Add, Dummy, Function, I, Integer, Rational, expand_mul,
                     oo, pi, zoo)
from ...core.function import ArgumentIndexError
from ..combinatorial.numbers import bernoulli, factorial, harmonic
from ..elementary.exponential import exp, exp_polar, log


###############################################################################
# #################### LERCH TRANSCENDENT ################################### #
###############################################################################


class lerchphi(Function):
    r"""
    Lerch transcendent (Lerch phi function).

    For `\operatorname{Re}(a) > 0`, `|z| < 1` and `s \in \mathbb{C}`, the
    Lerch transcendent is defined as

    .. math :: \Phi(z, s, a) = \sum_{n=0}^\infty \frac{z^n}{(n + a)^s},

    where the standard branch of the argument is used for `n + a`,
    and by analytic continuation for other values of the parameters.

    A commonly used related function is the Lerch zeta function, defined by

    .. math:: L(q, s, a) = \Phi(e^{2\pi i q}, s, a).

    **Analytic Continuation and Branching Behavior**

    It can be shown that

    .. math:: \Phi(z, s, a) = z\Phi(z, s, a+1) + a^{-s}.

    This provides the analytic continuation to `\operatorname{Re}(a) \le 0`.

    Assume now `\operatorname{Re}(a) > 0`. The integral representation

    .. math:: \Phi_0(z, s, a) = \int_0^\infty \frac{t^{s-1} e^{-at}}{1 - ze^{-t}}
                                \frac{\mathrm{d}t}{\Gamma(s)}

    provides an analytic continuation to `\mathbb{C} - [1, \infty)`.
    Finally, for `x \in (1, \infty)` we find

    .. math:: \lim_{\epsilon \to 0^+} \Phi_0(x + i\epsilon, s, a)
             -\lim_{\epsilon \to 0^+} \Phi_0(x - i\epsilon, s, a)
             = \frac{2\pi i \log^{s-1}{x}}{x^a \Gamma(s)},

    using the standard branch for both `\log{x}` and
    `\log{\log{x}}` (a branch of `\log{\log{x}}` is needed to
    evaluate `\log{x}^{s-1}`).
    This concludes the analytic continuation. The Lerch transcendent is thus
    branched at `z \in \{0, 1, \infty\}` and
    `a \in \mathbb{Z}_{\le 0}`. For fixed `z, a` outside these
    branch points, it is an entire function of `s`.

    See Also
    ========

    polylog, zeta

    References
    ==========

    * Bateman, H.; Erdélyi, A. (1953), Higher Transcendental Functions,
      Vol. I, New York: McGraw-Hill. Section 1.11.
    * https://dlmf.nist.gov/25.14
    * https://en.wikipedia.org/wiki/Lerch_transcendent

    Examples
    ========

    The Lerch transcendent is a fairly general function, for this reason it does
    not automatically evaluate to simpler functions. Use expand_func() to
    achieve this.

    If `z=1`, the Lerch transcendent reduces to the Hurwitz zeta function:

    >>> from diofant.abc import s
    >>> expand_func(lerchphi(1, s, a))
    zeta(s, a)

    More generally, if `z` is a root of unity, the Lerch transcendent
    reduces to a sum of Hurwitz zeta functions:

    >>> expand_func(lerchphi(-1, s, a))
    2**(-s)*zeta(s, a/2) - 2**(-s)*zeta(s, a/2 + 1/2)

    If `a=1`, the Lerch transcendent reduces to the polylogarithm:

    >>> expand_func(lerchphi(z, s, 1))
    polylog(s, z)/z

    More generally, if `a` is rational, the Lerch transcendent reduces
    to a sum of polylogarithms:

    >>> expand_func(lerchphi(z, s, Rational(1, 2)))
    2**(s - 1)*(polylog(s, sqrt(z))/sqrt(z) -
                polylog(s, sqrt(z)*exp_polar(I*pi))/sqrt(z))
    >>> expand_func(lerchphi(z, s, Rational(3, 2)))
    -2**s/z + 2**(s - 1)*(polylog(s, sqrt(z))/sqrt(z) -
                          polylog(s, sqrt(z)*exp_polar(I*pi))/sqrt(z))/z

    The derivatives with respect to `z` and `a` can be computed in
    closed form:

    >>> lerchphi(z, s, a).diff(z)
    (-a*lerchphi(z, s, a) + lerchphi(z, s - 1, a))/z
    >>> lerchphi(z, s, a).diff(a)
    -s*lerchphi(z, s + 1, a)

    """

    def _eval_expand_func(self, **hints):
        from .. import floor, unpolarify
        z, s, a = self.args
        if z == 1:
            return zeta(s, a)
        if s.is_Integer and s <= 0:
            t = Dummy('t')
            p = ((t + a)**(-s)).as_poly(t)
            start = 1/(1 - t)
            res = Integer(0)
            for c in p.all_coeffs():
                res += c*start
                start = t*start.diff(t)
            return res.subs({t: z})

        if a.is_Rational:
            # See section 18 of
            #   Kelly B. Roach.  Hypergeometric Function Representations.
            #   In: Proceedings of the 1997 International Symposium on Symbolic and
            #   Algebraic Computation, pages 205-211, New York, 1997. ACM.
            # TODO should something be polarified here?
            add = Integer(0)
            mul = Integer(1)
            # First reduce a to the interaval (0, 1]
            if a > 1:
                n = floor(a)
                if n == a:
                    n -= 1
                a -= n
                mul = z**(-n)
                add = Add(*[-z**(k - n)/(a + k)**s for k in range(n)])
            elif a <= 0:
                n = floor(-a) + 1
                a += n
                mul = z**n
                add = Add(*[z**(n - 1 - k)/(a - k - 1)**s for k in range(n)])

            m, n = Integer(a.numerator), Integer(a.denominator)
            zet = exp_polar(2*pi*I/n)
            root = z**(1/n)
            return add + mul*n**(s - 1)*Add(
                *[polylog(s, zet**k*root)._eval_expand_func(**hints)
                  / (unpolarify(zet)**k*root)**m for k in range(n)])

        # TODO use minpoly instead of ad-hoc methods when issue sympy/sympy#5888 is fixed
        if z.is_Exp and (z.exp/(pi*I)).is_Rational or z in [-1, I, -I]:
            # TODO reference?
            if z == -1:
                p, q = Integer(1), Integer(2)
            elif z == I:
                p, q = Integer(1), Integer(4)
            elif z == -I:
                p, q = Integer(-1), Integer(4)
            else:
                arg = z.exp/(2*pi*I)
                p, q = Integer(arg.numerator), Integer(arg.denominator)
            return Add(*[exp(2*pi*I*k*p/q)/q**s*zeta(s, (k + a)/q)
                         for k in range(q)])

        return lerchphi(z, s, a)

    def fdiff(self, argindex=1):
        z, s, a = self.args
        if argindex == 3:
            return -s*lerchphi(z, s + 1, a)
        if argindex == 1:
            return (lerchphi(z, s - 1, a) - a*lerchphi(z, s, a))/z
        raise ArgumentIndexError

    def _eval_rewrite_helper(self, z, s, a, target):
        res = self._eval_expand_func()
        if res.has(target):
            return res

    def _eval_rewrite_as_zeta(self, z, s, a):
        return self._eval_rewrite_helper(z, s, a, zeta)

    def _eval_rewrite_as_polylog(self, z, s, a):
        return self._eval_rewrite_helper(z, s, a, polylog)

###############################################################################
# #################### POLYLOGARITHM ######################################## #
###############################################################################


class polylog(Function):
    r"""
    Polylogarithm function.

    For `|z| < 1` and `s \in \mathbb{C}`, the polylogarithm is
    defined by

    .. math:: \operatorname{Li}_s(z) = \sum_{n=1}^\infty \frac{z^n}{n^s},

    where the standard branch of the argument is used for `n`. It admits
    an analytic continuation which is branched at `z=1` (notably not on the
    sheet of initial definition), `z=0` and `z=\infty`.

    The name polylogarithm comes from the fact that for `s=1`, the
    polylogarithm is related to the ordinary logarithm (see examples), and that

    .. math:: \operatorname{Li}_{s+1}(z) =
                    \int_0^z \frac{\operatorname{Li}_s(t)}{t} \mathrm{d}t.

    The polylogarithm is a special case of the Lerch transcendent:

    .. math:: \operatorname{Li}_{s}(z) = z \Phi(z, s, 1)

    See Also
    ========

    zeta, lerchphi

    Examples
    ========

    For `z \in \{0, 1, -1\}`, the polylogarithm is automatically expressed
    using other functions:

    >>> from diofant.abc import s
    >>> polylog(s, 0)
    0
    >>> polylog(s, 1)
    zeta(s)
    >>> polylog(s, -1)
    -dirichlet_eta(s)

    If `s` is a negative integer, `0` or `1`, the
    polylogarithm can be expressed using elementary functions. This can be
    done using expand_func():

    >>> expand_func(polylog(1, z))
    -log(-z + 1)
    >>> expand_func(polylog(0, z))
    z/(-z + 1)

    The derivative with respect to `z` can be computed in closed form:

    >>> polylog(s, z).diff(z)
    polylog(s - 1, z)/z

    The polylogarithm can be expressed in terms of the lerch transcendent:

    >>> polylog(s, z).rewrite(lerchphi)
    z*lerchphi(z, s, 1)

    References
    ==========

    * https://en.wikipedia.org/wiki/Polylogarithm
    * https://mathworld.wolfram.com/Polylogarithm.html

    """

    @classmethod
    def eval(cls, s, z):
        from .. import unpolarify
        if z == 1:
            return zeta(s)
        if z == -1:
            return -dirichlet_eta(s)
        if z == 0:
            return Integer(0)

        # branch handling
        if (1 - abs(z)).is_nonnegative:
            newz = unpolarify(z)
            if newz != z:
                return cls(s, newz)

    def fdiff(self, argindex=1):
        s, z = self.args
        if argindex == 2:
            return polylog(s - 1, z)/z
        raise ArgumentIndexError

    def _eval_rewrite_as_lerchphi(self, s, z):
        return z*lerchphi(z, s, 1)

    def _eval_expand_func(self, **hints):
        s, z = self.args
        if s == 1:
            return -log(1 - z)
        if s.is_Integer and s <= 0:
            u = Dummy('u')
            start = u/(1 - u)
            for _ in range(-s):
                start = u*start.diff(u)
            return expand_mul(start).subs({u: z})
        return polylog(s, z)

###############################################################################
# #################### HURWITZ GENERALIZED ZETA FUNCTION #################### #
###############################################################################


class zeta(Function):
    r"""
    Hurwitz zeta function (or Riemann zeta function).

    For `\operatorname{Re}(a) > 0` and `\operatorname{Re}(s) > 1`, this function is defined as

    .. math:: \zeta(s, a) = \sum_{n=0}^\infty \frac{1}{(n + a)^s},

    where the standard choice of argument for `n + a` is used. For fixed
    `a` with `\operatorname{Re}(a) > 0` the Hurwitz zeta function admits a
    meromorphic continuation to all of `\mathbb{C}`, it is an unbranched
    function with a simple pole at `s = 1`.

    Analytic continuation to other `a` is possible under some circumstances,
    but this is not typically done.

    The Hurwitz zeta function is a special case of the Lerch transcendent:

    .. math:: \zeta(s, a) = \Phi(1, s, a).

    This formula defines an analytic continuation for all possible values of
    `s` and `a` (also `\operatorname{Re}(a) < 0`), see the documentation of
    :class:`lerchphi` for a description of the branching behavior.

    If no value is passed for `a`, by this function assumes a default value
    of `a = 1`, yielding the Riemann zeta function.

    See Also
    ========

    dirichlet_eta, lerchphi, polylog

    References
    ==========

    * https://dlmf.nist.gov/25.11
    * https://en.wikipedia.org/wiki/Hurwitz_zeta_function

    Examples
    ========

    For `a = 1` the Hurwitz zeta function reduces to the famous Riemann
    zeta function:

    .. math:: \zeta(s, 1) = \zeta(s) = \sum_{n=1}^\infty \frac{1}{n^s}.

    >>> from diofant.abc import s
    >>> zeta(s, 1)
    zeta(s)
    >>> zeta(s)
    zeta(s)

    The Riemann zeta function can also be expressed using the Dirichlet eta
    function:

    >>> zeta(s).rewrite(dirichlet_eta)
    dirichlet_eta(s)/(-2**(-s + 1) + 1)

    The Riemann zeta function at positive even integer and negative odd integer
    values is related to the Bernoulli numbers:

    >>> zeta(2)
    pi**2/6
    >>> zeta(4)
    pi**4/90
    >>> zeta(-1)
    -1/12

    The specific formulae are:

    .. math:: \zeta(2n) = (-1)^{n+1} \frac{B_{2n} (2\pi)^{2n}}{2(2n)!}
    .. math:: \zeta(-n) = -\frac{B_{n+1}}{n+1}

    At negative even integers the Riemann zeta function is zero:

    >>> zeta(-4)
    0

    No closed-form expressions are known at positive odd integers, but
    numerical evaluation is possible:

    >>> zeta(3).evalf()
    1.20205690315959

    The derivative of `\zeta(s, a)` with respect to `a` is easily
    computed:

    >>> zeta(s, a).diff(a)
    -s*zeta(s + 1, a)

    However the derivative with respect to `s` has no useful closed form
    expression:

    >>> zeta(s, a).diff(s)
    Derivative(zeta(s, a), s)

    The Hurwitz zeta function can be expressed in terms of the Lerch transcendent,
    :class:`diofant.functions.special.zeta_functions.lerchphi`:

    >>> zeta(s, a).rewrite(lerchphi)
    lerchphi(1, s, a)

    """

    @classmethod
    def eval(cls, z, a=None):
        if a == 1:
            return cls(z)
        if a is None:
            a = Integer(1)
        if z.is_Number:
            if z is oo:
                return Integer(1)
            if z == 0:
                return Rational(1, 2) - a
            if z == 1:
                return zoo
            if z.is_Integer:
                if a.is_Integer:
                    if z.is_negative:
                        zeta = (-1)**z * bernoulli(-z + 1)/(-z + 1)
                    elif z.is_even:
                        B, F = bernoulli(z), factorial(z)
                        zeta = 2**(z - 1) * abs(B) * pi**z / F
                    else:
                        return

                    if a.is_negative:
                        return zeta + harmonic(abs(a), z)
                    return zeta - harmonic(a - 1, z)

    def _eval_rewrite_as_dirichlet_eta(self, s, a=1):
        if a != 1:
            return self
        s = self.args[0]
        return dirichlet_eta(s)/(1 - 2**(1 - s))

    def _eval_rewrite_as_lerchphi(self, s, a=1):
        return lerchphi(1, s, a)

    def fdiff(self, argindex=1):
        if len(self.args) == 2:
            s, a = self.args
        else:
            s, a = self.args + (1,)
        if argindex == 2:
            return -s*zeta(s + 1, a)
        raise ArgumentIndexError

    def _eval_rewrite_as_tractable(self, s, a=1, **kwargs):
        if len(self.args) == 1:
            return _zetas(exp(s))


class _zetas(Function):
    def _eval_rewrite_as_intractable(self, s, **kwargs):
        return zeta(log(s))

    def _eval_aseries(self, n, args0, x, logx):
        from ...calculus import Order
        point = args0[0]

        # Expansion at oo
        if point is oo:
            if n < 1:
                return Order(1, x)
            z = self.args[0]
            l = [(1/z)**log(p) for p in range(1, n)]
            o = Order(1/z**log(n), x)
            return (Add(*l))._eval_nseries(x, n, logx) + o

        # All other points are not handled
        return super()._eval_aseries(n, args0, x, logx)


class dirichlet_eta(Function):
    r"""
    Dirichlet eta function.

    For `\operatorname{Re}(s) > 0`, this function is defined as

    .. math:: \eta(s) = \sum_{n=1}^\infty \frac{(-1)^{n-1}}{n^s}.

    It admits a unique analytic continuation to all of `\mathbb{C}`.
    It is an entire, unbranched function.

    See Also
    ========

    zeta

    References
    ==========

    * https://en.wikipedia.org/wiki/Dirichlet_eta_function
    * https://mathworld.wolfram.com/DirichletEtaFunction.html

    Examples
    ========

    The Dirichlet eta function is closely related to the Riemann zeta function:

    >>> from diofant.abc import s
    >>> dirichlet_eta(s).rewrite(zeta)
    (-2**(-s + 1) + 1)*zeta(s)

    """

    @classmethod
    def eval(cls, s):
        if s == 1:
            return log(2)
        z = zeta(s)
        if not z.has(zeta):
            return (1 - 2**(1 - s))*z

    def _eval_rewrite_as_zeta(self, s):
        return (1 - 2**(1 - s)) * zeta(s)
