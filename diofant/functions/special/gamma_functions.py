from mpmath import mp, workprec

from ...core import (Add, Dummy, EulerGamma, Expr, Function, I, Integer, Pow,
                     Rational, oo, pi, zoo)
from ...core.function import ArgumentIndexError
from ..combinatorial.factorials import RisingFactorial, factorial, rf
from ..combinatorial.numbers import bernoulli, harmonic
from ..elementary.exponential import exp, log
from ..elementary.integers import ceiling, floor
from ..elementary.miscellaneous import sqrt
from .error_functions import erf
from .zeta_functions import zeta


###############################################################################
# ########################## COMPLETE GAMMA FUNCTION ######################## #
###############################################################################

class gamma(Function):
    r"""
    The gamma function

    .. math::
        \Gamma(x) := \int^{\infty}_{0} t^{x-1} e^{t} \mathrm{d}t.

    The ``gamma`` function implements the function which passes through the
    values of the factorial function, i.e. `\Gamma(n) = (n - 1)!` when n is
    an integer. More general, `\Gamma(z)` is defined in the whole complex
    plane except at the negative integers where there are simple poles.

    Examples
    ========

    Several special values are known:

    >>> gamma(1)
    1
    >>> gamma(4)
    6
    >>> gamma(Rational(3, 2))
    sqrt(pi)/2

    The Gamma function obeys the mirror symmetry:

    >>> conjugate(gamma(x))
    gamma(conjugate(x))

    Differentiation with respect to x is supported:

    >>> diff(gamma(x), x)
    gamma(x)*polygamma(0, x)

    Series expansion is also supported:

    >>> gamma(x).series(x, 0, 3)
    1/x - EulerGamma + x*(EulerGamma**2/2 + pi**2/12) + x**2*(-EulerGamma*pi**2/12 + polygamma(2, 1)/6 - EulerGamma**3/6) + O(x**3)

    We can numerically evaluate the gamma function to arbitrary precision
    on the whole complex plane:

    >>> gamma(pi).evalf(40)
    2.288037795340032417959588909060233922890
    >>> gamma(1+I).evalf(20)
    0.49801566811835604271 - 0.15494982830181068512*I

    See Also
    ========

    lowergamma: Lower incomplete gamma function.
    uppergamma: Upper incomplete gamma function.
    polygamma: Polygamma function.
    loggamma: Log Gamma function.
    digamma: Digamma function.
    trigamma: Trigamma function.
    diofant.functions.special.beta_functions.beta: Euler Beta function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Gamma_function
    * https://dlmf.nist.gov/5
    * https://mathworld.wolfram.com/GammaFunction.html
    * http://functions.wolfram.com/GammaBetaErf/Gamma/

    """

    unbranched = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return self.func(self.args[0])*polygamma(0, self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is oo:
                return oo
            elif arg.is_Integer:
                if arg.is_positive:
                    return factorial(arg - 1)
                else:
                    return zoo
            elif arg.is_Rational:
                if arg.denominator == 2:
                    n = abs(arg.numerator) // arg.denominator

                    if arg.is_positive:
                        k, coeff = n, Integer(1)
                    else:
                        n = k = n + 1

                        if n & 1 == 0:
                            coeff = Integer(1)
                        else:
                            coeff = Integer(-1)

                    for i in range(3, 2*k, 2):
                        coeff *= i

                    if arg.is_positive:
                        return coeff*sqrt(pi) / 2**n
                    else:
                        return 2**n*sqrt(pi) / coeff

        if arg.is_integer and arg.is_nonpositive:
            return zoo

    def _eval_expand_func(self, **hints):
        arg = self.args[0]
        if arg.is_Rational:
            if abs(arg.numerator) > arg.denominator:
                x = Dummy('x')
                n = arg.numerator // arg.denominator
                p = arg.numerator - n*arg.denominator
                return self.func(x + n)._eval_expand_func().subs({x: Rational(p, arg.denominator)})

        if arg.is_Add:
            coeff, tail = arg.as_coeff_add()
            if coeff and coeff.denominator != 1:
                intpart = floor(coeff)
                tail = (coeff - intpart,) + tail
                coeff = intpart
            tail = arg._new_rawargs(*tail, reeval=False)
            return self.func(tail)*RisingFactorial(tail, coeff)

        return self.func(*self.args)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_is_extended_real(self):
        x = self.args[0]
        if x.is_positive or x.is_noninteger:
            return True

    def _eval_is_positive(self):
        x = self.args[0]
        if x.is_positive:
            return True
        elif x.is_noninteger:
            return floor(x).is_even

    def _eval_rewrite_as_tractable(self, z, **kwargs):
        return exp(loggamma(z))

    def _eval_rewrite_as_factorial(self, z):
        return factorial(z - 1)

    def _eval_nseries(self, x, n, logx):
        x0 = self.args[0].limit(x, 0)
        if not (x0.is_Integer and x0 <= 0):
            return super()._eval_nseries(x, n, logx)
        t = self.args[0] - x0
        return (self.func(t + 1)/rf(self.args[0], -x0 + 1))._eval_nseries(x, n, logx)

    def _latex(self, printer, exp=None):
        aa = printer._print(self.args[0])
        if exp:
            return f'\\Gamma^{{{printer._print(exp)}}}{{\\left({aa} \\right)}}'
        else:
            return f'\\Gamma{{\\left({aa} \\right)}}'

    @staticmethod
    def _latex_no_arg(printer):
        return r'\Gamma'


###############################################################################
# ################ LOWER and UPPER INCOMPLETE GAMMA FUNCTIONS ############### #
###############################################################################

class lowergamma(Function):
    r"""
    The lower incomplete gamma function.

    It can be defined as the meromorphic continuation of

    .. math::
        \gamma(s, x) := \int_0^x t^{s-1} e^{-t} \mathrm{d}t = \Gamma(s) - \Gamma(s, x).

    This can be shown to be the same as

    .. math::
        \gamma(s, x) = \frac{x^s}{s} {}_1F_1\left({s \atop s+1} \middle| -x\right),

    where `{}_1F_1` is the (confluent) hypergeometric function.

    Examples
    ========

    >>> from diofant.abc import s
    >>> lowergamma(s, x)
    lowergamma(s, x)
    >>> lowergamma(3, x)
    2 - E**(-x)*x**2 - 2*E**(-x)*x - 2*E**(-x)
    >>> lowergamma(-Rational(1, 2), x)
    -2*sqrt(pi)*erf(sqrt(x)) - 2*E**(-x)/sqrt(x)

    See Also
    ========

    gamma: Gamma function.
    uppergamma: Upper incomplete gamma function.
    polygamma: Polygamma function.
    loggamma: Log Gamma function.
    digamma: Digamma function.
    trigamma: Trigamma function.
    diofant.functions.special.beta_functions.beta: Euler Beta function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Incomplete_gamma_function#Lower_incomplete_Gamma_function
    * Abramowitz, Milton; Stegun, Irene A., eds. (1965), Chapter 6, Section 5,
      Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables
    * https://dlmf.nist.gov/8
    * http://functions.wolfram.com/GammaBetaErf/Gamma2/
    * http://functions.wolfram.com/GammaBetaErf/Gamma3/

    """

    def fdiff(self, argindex=2):
        from .. import unpolarify
        from .hyper import meijerg
        if argindex == 2:
            a, z = self.args
            return exp(-unpolarify(z))*z**(a - 1)
        elif argindex == 1:
            a, z = self.args
            return gamma(a)*digamma(a) - log(z)*uppergamma(a, z) \
                - meijerg([], [1, 1], [0, 0, a], [], z)

        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, a, x):
        # For lack of a better place, we use this one to extract branching
        # information. The following can be
        # found in the literature (c/f references given above), albeit scattered:
        # 1) For fixed x != 0, lowergamma(s, x) is an entire function of s
        # 2) For fixed positive integers s, lowergamma(s, x) is an entire
        #    function of x.
        # 3) For fixed non-positive integers s,
        #    lowergamma(s, exp(I*2*pi*n)*x) =
        #              2*pi*I*n*(-1)**(-s)/factorial(-s) + lowergamma(s, x)
        #    (this follows from lowergamma(s, x).diff(x) = x**(s-1)*exp(-x)).
        # 4) For fixed non-integral s,
        #    lowergamma(s, x) = x**s*gamma(s)*lowergamma_unbranched(s, x),
        #    where lowergamma_unbranched(s, x) is an entire function (in fact
        #    of both s and x), i.e.
        #    lowergamma(s, exp(2*I*pi*n)*x) = exp(2*pi*I*n*a)*lowergamma(a, x)
        from .. import unpolarify
        nx, n = x.extract_branch_factor()
        if a.is_integer and a.is_positive:
            nx = unpolarify(x)
            if nx != x:
                return lowergamma(a, nx)
        elif a.is_integer and a.is_nonpositive:
            if n != 0:
                return 2*pi*I*n*(-1)**(-a)/factorial(-a) + lowergamma(a, nx)
        elif n != 0:
            return exp(2*pi*I*n*a)*lowergamma(a, nx)

        # Special values.
        if a.is_Number:
            # TODO this should be non-recursive
            if a == 1:
                return 1 - exp(-x)
            elif a == Rational(1, 2):
                return sqrt(pi)*erf(sqrt(x))
            elif a.is_Integer or (2*a).is_Integer:
                b = a - 1
                if b.is_positive:
                    return b*cls(b, x) - x**b * exp(-x)
                elif a == 0:
                    return zoo

                if not a.is_Integer:
                    return (cls(a + 1, x) + x**a * exp(-x))/a

    def _eval_evalf(self, prec):
        a = self.args[0]._to_mpmath(prec)
        z = self.args[1]._to_mpmath(prec)
        with workprec(prec):
            res = mp.gammainc(a, 0, z)
        return Expr._from_mpmath(res, prec)

    def _eval_conjugate(self):
        z = self.args[1]
        if z not in (0, -oo):
            return self.func(self.args[0].conjugate(), z.conjugate())

    def _eval_rewrite_as_uppergamma(self, s, x):
        return gamma(s) - uppergamma(s, x)

    def _eval_rewrite_as_tractable(self, s, x, **kwargs):
        return self.rewrite(uppergamma)

    def _eval_rewrite_as_expint(self, s, x):
        from .error_functions import expint
        if s.is_integer and s.is_nonpositive:
            return self
        return self.rewrite(uppergamma).rewrite(expint)

    @staticmethod
    def _latex_no_arg(printer):
        return r'\gamma'


class uppergamma(Function):
    r"""
    The upper incomplete gamma function.

    It can be defined as the meromorphic continuation of

    .. math::
        \Gamma(s, x) := \int_x^\infty t^{s-1} e^{-t} \mathrm{d}t = \Gamma(s) - \gamma(s, x).

    where `\gamma(s, x)` is the lower incomplete gamma function,
    :class:`lowergamma`. This can be shown to be the same as

    .. math::
        \Gamma(s, x) = \Gamma(s) - \frac{x^s}{s} {}_1F_1\left({s \atop s+1} \middle| -x\right),

    where `{}_1F_1` is the (confluent) hypergeometric function.

    The upper incomplete gamma function is also essentially equivalent to the
    generalized exponential integral:

    .. math::
        \operatorname{E}_{n}(x) = \int_{1}^{\infty}{\frac{e^{-xt}}{t^n} \, dt} = x^{n-1}\Gamma(1-n,x).

    Examples
    ========

    >>> from diofant.abc import s
    >>> uppergamma(s, x)
    uppergamma(s, x)
    >>> uppergamma(3, x)
    E**(-x)*x**2 + 2*E**(-x)*x + 2*E**(-x)
    >>> uppergamma(-Rational(1, 2), x)
    -2*sqrt(pi)*(-erf(sqrt(x)) + 1) + 2*E**(-x)/sqrt(x)
    >>> uppergamma(-2, x)
    expint(3, x)/x**2

    See Also
    ========

    gamma: Gamma function.
    lowergamma: Lower incomplete gamma function.
    polygamma: Polygamma function.
    loggamma: Log Gamma function.
    digamma: Digamma function.
    trigamma: Trigamma function.
    diofant.functions.special.beta_functions.beta: Euler Beta function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Incomplete_gamma_function#Upper_incomplete_Gamma_function
    * Abramowitz, Milton; Stegun, Irene A., eds. (1965), Chapter 6, Section 5,
      Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables
    * https://dlmf.nist.gov/8
    * http://functions.wolfram.com/GammaBetaErf/Gamma2/
    * http://functions.wolfram.com/GammaBetaErf/Gamma3/
    * https://en.wikipedia.org/wiki/Exponential_integral#Relation_with_other_functions

    """

    def fdiff(self, argindex=2):
        from .. import unpolarify
        from .hyper import meijerg
        if argindex == 2:
            a, z = self.args
            return -exp(-unpolarify(z))*z**(a - 1)
        elif argindex == 1:
            a, z = self.args
            return uppergamma(a, z)*log(z) + meijerg([], [1, 1], [0, 0, a], [], z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_evalf(self, prec):
        a = self.args[0]._to_mpmath(prec)
        z = self.args[1]._to_mpmath(prec)
        with workprec(prec):
            res = mp.gammainc(a, z, mp.inf)
        return Expr._from_mpmath(res, prec)

    @classmethod
    def eval(cls, a, z):
        from .. import unpolarify
        from .error_functions import expint
        if z.is_Number:
            if z is oo:
                return Integer(0)
            elif z == 0:
                # TODO: Holds only for Re(a) > 0:
                return gamma(a)

        # We extract branching information here. C/f lowergamma.
        nx, n = z.extract_branch_factor()
        if a.is_integer and a.is_positive:
            nx = unpolarify(z)
            if z != nx:
                return uppergamma(a, nx)
        elif a.is_integer and a.is_nonpositive:
            if n != 0:
                return -2*pi*I*n*(-1)**(-a)/factorial(-a) + uppergamma(a, nx)
        elif n != 0:
            return gamma(a)*(1 - exp(2*pi*I*n*a)) + exp(2*pi*I*n*a)*uppergamma(a, nx)

        # Special values.
        if a.is_Number:
            # TODO this should be non-recursive
            if a == 1:
                return exp(-z)
            elif a == Rational(1, 2):
                return sqrt(pi)*(1 - erf(sqrt(z)))  # TODO could use erfc...
            elif a.is_Integer or (2*a).is_Integer:
                b = a - 1
                if b.is_positive:
                    return b*cls(b, z) + z**b * exp(-z)
                elif b.is_Integer:
                    return expint(-b, z)*unpolarify(z)**(b + 1)
                else:
                    return (cls(a + 1, z) - z**a * exp(-z))/a

    def _eval_conjugate(self):
        z = self.args[1]
        if z not in (0, -oo):
            return self.func(self.args[0].conjugate(), z.conjugate())

    def _eval_rewrite_as_lowergamma(self, s, x):
        return gamma(s) - lowergamma(s, x)

    def _eval_rewrite_as_expint(self, s, x):
        from .error_functions import expint
        return expint(1 - s, x)*x**s


###############################################################################
# #################### POLYGAMMA and LOGGAMMA FUNCTIONS ##################### #
###############################################################################

class polygamma(Function):
    r"""
    The function ``polygamma(n, z)`` returns ``log(gamma(z)).diff(n + 1)``.

    It is a meromorphic function on `\mathbb{C}` and defined as the (n+1)-th
    derivative of the logarithm of the gamma function:

    .. math::
        \psi^{(n)} (z) := \frac{\mathrm{d}^{n+1}}{\mathrm{d} z^{n+1}} \log\Gamma(z).

    Examples
    ========

    Several special values are known:

    >>> polygamma(0, 1)
    -EulerGamma
    >>> polygamma(0, Rational(1, 2))
    -2*log(2) - EulerGamma
    >>> polygamma(0, Rational(1, 3))
    -3*log(3)/2 - sqrt(3)*pi/6 - EulerGamma
    >>> polygamma(0, Rational(1, 4))
    -3*log(2) - pi/2 - EulerGamma
    >>> polygamma(0, 2)
    -EulerGamma + 1
    >>> polygamma(0, 23)
    -EulerGamma + 19093197/5173168

    >>> polygamma(0, oo)
    oo
    >>> polygamma(0, -oo)
    oo
    >>> polygamma(0, I*oo)
    oo
    >>> polygamma(0, -I*oo)
    oo

    Differentiation with respect to x is supported:

    >>> diff(polygamma(0, x), x)
    polygamma(1, x)
    >>> diff(polygamma(0, x), (x, 2))
    polygamma(2, x)
    >>> diff(polygamma(0, x), (x, 3))
    polygamma(3, x)
    >>> diff(polygamma(1, x), x)
    polygamma(2, x)
    >>> diff(polygamma(1, x), (x, 2))
    polygamma(3, x)
    >>> diff(polygamma(2, x), x)
    polygamma(3, x)
    >>> diff(polygamma(2, x), (x, 2))
    polygamma(4, x)

    >>> diff(polygamma(n, x), x)
    polygamma(n + 1, x)
    >>> diff(polygamma(n, x), (x, 2))
    polygamma(n + 2, x)

    We can rewrite polygamma functions in terms of harmonic numbers:

    >>> polygamma(0, x).rewrite(harmonic)
    harmonic(x - 1) - EulerGamma
    >>> polygamma(2, x).rewrite(harmonic)
    2*harmonic(x - 1, 3) - 2*zeta(3)
    >>> ni = Symbol('n', integer=True)
    >>> polygamma(ni, x).rewrite(harmonic)
    (-1)**(n + 1)*(-harmonic(x - 1, n + 1) + zeta(n + 1))*factorial(n)

    See Also
    ========

    gamma: Gamma function.
    lowergamma: Lower incomplete gamma function.
    uppergamma: Upper incomplete gamma function.
    loggamma: Log Gamma function.
    digamma: Digamma function.
    trigamma: Trigamma function.
    diofant.functions.special.beta_functions.beta: Euler Beta function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Polygamma_function
    * https://mathworld.wolfram.com/PolygammaFunction.html
    * http://functions.wolfram.com/GammaBetaErf/PolyGamma/
    * http://functions.wolfram.com/GammaBetaErf/PolyGamma2/

    """

    def fdiff(self, argindex=2):
        if argindex == 2:
            n, z = self.args[:2]
            return polygamma(n + 1, z)
        else:
            raise ArgumentIndexError(self, argindex)

    def _eval_is_positive(self):
        if self.args[0].is_positive and self.args[1].is_positive:
            if self.args[0].is_odd:
                return True

    def _eval_aseries(self, n, args0, x, logx):
        from ...series import Order
        if args0[1] != oo or not \
                (self.args[0].is_Integer and self.args[0].is_nonnegative):
            return super()._eval_aseries(n, args0, x, logx)
        z = self.args[1]
        N = self.args[0]

        if N == 0:
            # digamma function series
            # Abramowitz & Stegun, p. 259, 6.3.18
            r = log(z) - 1/(2*z)
            o = None
            if n < 2:
                o = Order(1/z, x)
            else:
                m = ceiling((n + 1)//2)
                l = [bernoulli(2*k) / (2*k*z**(2*k)) for k in range(1, m)]
                r -= Add(*l)
                o = Order(1/z**(2*m), x)
            return r._eval_nseries(x, n, logx) + o
        else:
            # proper polygamma function
            # Abramowitz & Stegun, p. 260, 6.4.10
            # We return terms to order higher than O(x**n) on purpose
            # -- otherwise we would not be able to return any terms for
            #    quite a long time!
            fac = gamma(N)
            e0 = fac + N*fac/(2*z)
            m = ceiling((n + 1)//2)
            for k in range(1, m):
                fac = fac*(2*k + N - 1)*(2*k + N - 2) / ((2*k)*(2*k - 1))
                e0 += bernoulli(2*k)*fac/z**(2*k)
            o = Order(1/z**(2*m), x)
            if n == 0:
                o = Order(1/z, x)
            elif n == 1:
                o = Order(1/z**2, x)
            r = e0 + o
            return (-1 * (-1/z)**N * r)._eval_nseries(x, n, logx)

    @classmethod
    def eval(cls, n, z):
        from .. import unpolarify

        if n.is_integer:
            if n.is_nonnegative:
                nz = unpolarify(z)
                if z != nz:
                    return polygamma(n, nz)

            if n == -1:
                return loggamma(z)
            else:
                if z.is_Number:
                    if z is oo:
                        if n.is_Number:
                            if n == 0:
                                return oo
                            else:
                                return Integer(0)
                    elif z.is_Integer:
                        if z.is_nonpositive:
                            return zoo
                        else:
                            if n == 0:
                                return -EulerGamma + harmonic(z - 1, 1)
                            elif n.is_odd:
                                return (-1)**(n + 1)*factorial(n)*zeta(n + 1, z)

        if n == 0:
            if z.is_Rational:
                # TODO actually *any* n/m can be done, but that is messy
                lookup = {Rational(1, 2): -2*log(2) - EulerGamma,
                          Rational(1, 3): -pi/2/sqrt(3) - 3*log(3)/2 - EulerGamma,
                          Rational(1, 4): -pi/2 - 3*log(2) - EulerGamma,
                          Rational(3, 4): -3*log(2) - EulerGamma + pi/2,
                          Rational(2, 3): -3*log(3)/2 + pi/2/sqrt(3) - EulerGamma}
                if z > 0:
                    n = floor(z)
                    z0 = z - n
                    if z0 in lookup:
                        return lookup[z0] + Add(*[1/(z0 + k) for k in range(n)])
                else:  # z < 0
                    n = floor(1 - z)
                    z0 = z + n
                    if z0 in lookup:
                        return lookup[z0] - Add(*[1/(z0 - 1 - k) for k in range(n)])
            elif z in (oo, -oo):
                return oo
            else:
                t = z.as_coefficient(I)
                if t in (oo, -oo):
                    return oo

        # TODO n == 1 also can do some rational z

    def _eval_expand_func(self, **hints):
        n, z = self.args

        if n.is_Integer and n.is_nonnegative:
            if z.is_Add:
                coeff = z.args[0]
                if coeff.is_Integer:
                    e = -(n + 1)
                    if coeff > 0:
                        tail = Add(*[Pow(z - i, e)
                                     for i in range(1, int(coeff) + 1)])
                    else:
                        tail = -Add(*[Pow(z + i, e)
                                      for i in range(int(-coeff))])
                    return polygamma(n, z - coeff) + (-1)**n*factorial(n)*tail

            elif z.is_Mul:
                coeff, z = z.as_two_terms()
                if coeff.is_Integer and coeff.is_positive:
                    tail = [polygamma(n, z + Rational(i, coeff))
                            for i in range(int(coeff))]
                    if n == 0:
                        return Add(*tail)/coeff + log(coeff)
                    else:
                        return Add(*tail)/coeff**(n + 1)
                z *= coeff

        return polygamma(n, z)

    def _eval_rewrite_as_zeta(self, n, z):
        if (n - 1).is_nonnegative:
            return (-1)**(n + 1)*factorial(n)*zeta(n + 1, z)
        else:
            return self

    def _eval_rewrite_as_harmonic(self, n, z):
        if n.is_integer:
            if n == 0:
                return harmonic(z - 1) - EulerGamma
            else:
                return (-1)**(n+1) * factorial(n) * (zeta(n+1) - harmonic(z-1, n+1))

    def _eval_as_leading_term(self, x):
        from ...series import Order
        n, z = (a.as_leading_term(x) for a in self.args)
        o = Order(z, x)
        if n == 0 and o.contains(1/x):
            return o.getn() * log(x)
        else:
            return self.func(n, z)


class loggamma(Function):
    r"""
    The ``loggamma`` function implements the logarithm of the
    gamma function i.e, `\log\Gamma(x)`.

    Examples
    ========

    Several special values are known. For numerical integral
    arguments we have:

    >>> loggamma(-2)
    oo
    >>> loggamma(0)
    oo
    >>> loggamma(1)
    0
    >>> loggamma(2)
    0
    >>> loggamma(3)
    log(2)

    and for symbolic values:

    >>> n = Symbol('n', integer=True, positive=True)
    >>> loggamma(n)
    log(gamma(n))
    >>> loggamma(-n)
    oo

    for half-integral values:

    >>> loggamma(Rational(5, 2))
    log(3*sqrt(pi)/4)
    >>> loggamma(n/2)
    log(2**(-n + 1)*sqrt(pi)*gamma(n)/gamma(n/2 + 1/2))

    and general rational arguments:

    >>> L = loggamma(Rational(16, 3))
    >>> expand_func(L).doit()
    -5*log(3) + loggamma(1/3) + log(4) + log(7) + log(10) + log(13)
    >>> L = loggamma(Rational(19, 4))
    >>> expand_func(L).doit()
    -4*log(4) + loggamma(3/4) + log(3) + log(7) + log(11) + log(15)
    >>> L = loggamma(Rational(23, 7))
    >>> expand_func(L).doit()
    -3*log(7) + log(2) + loggamma(2/7) + log(9) + log(16)

    The loggamma function has the following limits towards infinity:

    >>> loggamma(oo)
    oo
    >>> loggamma(-oo)
    zoo

    The loggamma function obeys the mirror symmetry
    if `x \in \mathbb{C} \setminus \{-\infty, 0\}`:

    >>> c = Symbol('c', complex=True, real=False)
    >>> conjugate(loggamma(c))
    loggamma(conjugate(c))

    Differentiation with respect to x is supported:

    >>> diff(loggamma(x), x)
    polygamma(0, x)

    Series expansion is also supported:

    >>> loggamma(x).series(x, 0, 4)
    -log(x) - EulerGamma*x + pi**2*x**2/12 + x**3*polygamma(2, 1)/6 + O(x**4)

    We can numerically evaluate the gamma function to arbitrary precision
    on the whole complex plane:

    >>> loggamma(5).evalf(30)
    3.17805383034794561964694160130
    >>> loggamma(I).evalf(20)
    -0.65092319930185633889 - 1.8724366472624298171*I

    See Also
    ========

    gamma: Gamma function.
    lowergamma: Lower incomplete gamma function.
    uppergamma: Upper incomplete gamma function.
    polygamma: Polygamma function.
    digamma: Digamma function.
    trigamma: Trigamma function.
    diofant.functions.special.beta_functions.beta: Euler Beta function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Gamma_function
    * https://dlmf.nist.gov/5
    * https://mathworld.wolfram.com/LogGammaFunction.html
    * http://functions.wolfram.com/GammaBetaErf/LogGamma/

    """

    @classmethod
    def eval(cls, z):
        if z.is_integer:
            if z.is_nonpositive:
                return oo
            elif z.is_positive:
                return log(gamma(z))
        elif z.is_rational:
            p, q = z.as_numer_denom()
            # Half-integral values:
            if p.is_positive and q == 2:
                return log(sqrt(pi) * 2**(1 - p) * gamma(p) / gamma((p + 1) / 2))

        if z is oo:
            return oo
        elif abs(z) is oo:
            return zoo

    def _eval_expand_func(self, **hints):
        from ...concrete import Sum
        z = self.args[0]

        if z.is_Rational:
            p, q = z.as_numer_denom()
            # General rational arguments (u + p/q)
            # Split z as n + p/q with p < q
            n = p // q
            p = p - n*q
            assert p.is_positive
            assert q.is_positive
            assert p < q
            k = Dummy('k')
            if n.is_positive:
                return loggamma(p / q) - n*log(q) + Sum(log((k - 1)*q + p), (k, 1, n))
            elif n.is_negative:
                return loggamma(p / q) - n*log(q) + pi*I*n - Sum(log(k*q - p), (k, 1, -n))

        return self

    def _eval_nseries(self, x, n, logx=None):
        x0 = self.args[0].limit(x, 0)
        if x0 == 0:
            f = self._eval_rewrite_as_intractable(*self.args)
            return f._eval_nseries(x, n, logx)
        return super()._eval_nseries(x, n, logx)

    def _eval_aseries(self, n, args0, x, logx):
        from ...series import Order
        if args0[0] != oo:
            return super()._eval_aseries(n, args0, x, logx)
        z = self.args[0]
        m = min(n, ceiling(Rational(n + 1, 2)))
        r = log(z)*(z - Rational(1, 2)) - z + log(2*pi)/2
        l = [bernoulli(2*k) / (2*k*(2*k - 1)*z**(2*k - 1)) for k in range(1, m)]
        o = None
        if m == 0:
            o = Order(1, x)
        else:
            o = Order(1/z**(2*m - 1), x)
        # It is very inefficient to first add the order and then do the nseries
        return (r + Add(*l))._eval_nseries(x, n, logx) + o

    def _eval_rewrite_as_intractable(self, z, **kwargs):
        return log(gamma(z))

    def _eval_is_extended_real(self):
        if self.args[0].is_nonnegative:
            return True

    def _eval_conjugate(self):
        z = self.args[0]
        if (z.is_extended_real and z.is_nonpositive) is False:
            return self.func(z.conjugate())

    def fdiff(self, argindex=1):
        if argindex == 1:
            return polygamma(0, self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)


def digamma(x):
    r"""
    The digamma function is the first derivative of the loggamma function i.e,

    .. math::
        \psi(x) := \frac{\mathrm{d}}{\mathrm{d} z} \log\Gamma(z)
                = \frac{\Gamma'(z)}{\Gamma(z) }

    In this case, ``digamma(z) = polygamma(0, z)``.

    See Also
    ========

    gamma: Gamma function.
    lowergamma: Lower incomplete gamma function.
    uppergamma: Upper incomplete gamma function.
    polygamma: Polygamma function.
    loggamma: Log Gamma function.
    trigamma: Trigamma function.
    diofant.functions.special.beta_functions.beta: Euler Beta function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Digamma_function
    * https://mathworld.wolfram.com/DigammaFunction.html
    * http://functions.wolfram.com/GammaBetaErf/PolyGamma2/

    """
    return polygamma(0, x)


def trigamma(x):
    r"""
    The trigamma function is the second derivative of the loggamma function i.e,

    .. math::
        \psi^{(1)}(z) := \frac{\mathrm{d}^{2}}{\mathrm{d} z^{2}} \log\Gamma(z).

    In this case, ``trigamma(z) = polygamma(1, z)``.

    See Also
    ========

    gamma: Gamma function.
    lowergamma: Lower incomplete gamma function.
    uppergamma: Upper incomplete gamma function.
    polygamma: Polygamma function.
    loggamma: Log Gamma function.
    digamma: Digamma function.
    diofant.functions.special.beta_functions.beta: Euler Beta function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigamma_function
    * https://mathworld.wolfram.com/TrigammaFunction.html
    * http://functions.wolfram.com/GammaBetaErf/PolyGamma2/

    """
    return polygamma(1, x)
