from ...core import (Add, E, Function, I, Integer, Mul, Pow, expand_log, nan,
                     oo, pi, zoo)
from ...core.function import ArgumentIndexError, _coeff_isneg
from ...core.sympify import sympify
from ...ntheory import multiplicity, perfect_power
from .miscellaneous import sqrt


class exp_polar(Function):
    r"""
    Represent a 'polar number' (see g-function Sphinx documentation).

    ``exp_polar`` represents the function
    `Exp: \mathbb{C} \rightarrow \mathcal{S}`, sending the complex number
    `z = a + bi` to the polar number `r = exp(a), \theta = b`. It is one of
    the main functions to construct polar numbers.

    The main difference is that polar numbers don't "wrap around" at `2 \pi`:

    >>> exp(2*pi*I)
    1
    >>> exp_polar(2*pi*I)
    exp_polar(2*I*pi)

    apart from that they behave mostly like classical complex numbers:

    >>> exp_polar(2)*exp_polar(3)
    exp_polar(5)

    See also
    ========

    diofant.simplify.powsimp.powsimp
    diofant.functions.elementary.complexes.polar_lift
    diofant.functions.elementary.complexes.periodic_argument
    diofant.functions.elementary.complexes.principal_branch

    """

    is_polar = True
    is_comparable = False  # cannot be evalf'd

    unbranched = True

    def _eval_as_numer_denom(self):
        """
        Returns this with a positive exponent as a 2-tuple (a fraction).

        Examples
        ========

        >>> exp(-x).as_numer_denom()
        (1, E**x)
        >>> exp(x).as_numer_denom()
        (E**x, 1)

        """
        # this should be the same as Pow.as_numer_denom wrt
        # exponent handling
        exp = self.exp
        neg_exp = exp.is_negative
        if not neg_exp and not (-exp).is_negative:
            neg_exp = _coeff_isneg(exp)
        if neg_exp:
            return Integer(1), self.func(-exp)
        return self, Integer(1)

    @property
    def exp(self):
        """Returns the exponent of the function."""
        return self.args[0]

    def _eval_conjugate(self):
        return self.func(self.exp.conjugate())

    def _eval_is_finite(self):
        arg = self.exp
        if arg.is_infinite:
            if arg.is_positive:
                return False
            elif arg.is_negative:
                return True
        if arg.is_finite:
            return True

    def _eval_is_rational(self):
        if self.exp == 0:
            return True
        elif self.exp.is_rational and self.exp.is_nonzero:
            return False

    def _eval_is_zero(self):
        if self.exp.is_infinite and self.exp.is_negative:
            return True

    def _eval_expand_power_exp(self, **hints):
        arg = self.exp
        if arg.is_Add and arg.is_commutative:
            expr = 1
            for x in arg.args:
                expr *= self.func(x)
            return expr
        return self.func(arg)

    def _eval_Abs(self):
        from ...core import expand_mul
        return sqrt(expand_mul(self * self.conjugate()))

    def _eval_evalf(self, prec):
        """Careful! any evalf of polar numbers is flaky."""
        from .complexes import im, re
        i = im(self.exp)
        try:
            bad = (i <= -pi or i > pi)
        except TypeError:
            bad = True
        if bad:
            return self  # cannot evalf for this argument
        res = exp(self.exp).evalf(prec, strict=False)
        if i > 0 and im(res) < 0:
            # i ~ pi, but exp(I*i) evaluated to argument slightly bigger than pi
            return re(res)
        return res

    def _eval_power(self, other):
        return self.func(self.exp*other)

    def _eval_is_extended_real(self):
        if self.exp.is_extended_real:
            return True

    def as_base_exp(self):
        if self.exp == 0:
            return self, Integer(1)
        return self.func(1), Mul(*self.args)


def exp(arg, **kwargs):
    """
    The exponential function, `e^x`.

    See Also
    ========

    diofant.functions.elementary.exponential.log

    """
    return Pow(E, arg, **kwargs)


class log(Function):
    r"""
    The natural logarithm function `\ln(x)` or `\log(x)`.
    Logarithms are taken with the natural base, `e`. To get
    a logarithm of a different base ``b``, use ``log(x, b)``,
    which is essentially short-hand for ``log(x)/log(b)``.

    See Also
    ========

    diofant.functions.elementary.exponential.exp

    """

    def fdiff(self, argindex=1):
        """Returns the first derivative of the function."""
        if argindex == 1:
            return 1/self.args[0]
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        r"""Returns `e^x`, the inverse function of `\log(x)`."""
        return exp

    @classmethod
    def eval(cls, arg, base=None):
        from .complexes import unpolarify
        arg = sympify(arg)

        if base is not None:
            base = sympify(base)
            if base == 1:
                if arg == 1:
                    return nan
                else:
                    return zoo
            try:
                # handle extraction of powers of the base now
                # or else expand_log in Mul would have to handle this
                n = multiplicity(base, arg)
                if n:
                    den = base**n
                    if den.is_Integer:
                        return n + log(arg // den) / log(base)
                    else:
                        return n + log(arg / den) / log(base)
                else:
                    return log(arg)/log(base)
            except ValueError:
                pass
            if base is not E:
                return cls(arg)/cls(base)
            else:
                return cls(arg)

        if arg.is_Number:
            if arg == 0:
                return zoo
            elif arg == 1:
                return Integer(0)
            elif arg in (oo, -oo):
                return oo
            elif arg.is_Rational:
                if arg.denominator != 1:
                    return cls(arg.numerator) - cls(arg.denominator)

        if arg.is_Pow and arg.base is E and arg.exp.is_extended_real:
            return arg.exp
        elif isinstance(arg, exp_polar):
            return unpolarify(arg.exp)

        if arg.is_number:
            if arg.is_negative:
                return pi * I + cls(-arg)
            elif arg is zoo:
                return zoo
            elif arg is E:
                return Integer(1)

        # don't autoexpand Pow or Mul (see the issue sympy/sympy#3351):
        if not arg.is_Add:
            coeff = arg.as_coefficient(I)

            if coeff is not None:
                if coeff in (oo, -oo):
                    return oo
                elif coeff.is_Rational:
                    if coeff.is_nonnegative:
                        return +pi*I/2 + cls(+coeff)
                    else:
                        return -pi*I/2 + cls(-coeff)

    def as_base_exp(self):
        """Returns this function in the form (base, exponent)."""
        return self, Integer(1)

    def _eval_expand_log(self, deep=True, **hints):
        from ...concrete import Product, Sum
        from .complexes import unpolarify
        force = hints.get('force', False)
        if (len(self.args) == 2):
            return expand_log(self.func(*self.args), deep=deep, force=force)
        arg = self.args[0]
        if arg.is_Integer:
            # remove perfect powers
            p = perfect_power(int(arg))
            if p is not False:
                return p[1]*self.func(p[0])
        elif arg.is_Mul:
            expr = []
            nonpos = []
            for x in arg.args:
                if force or x.is_positive or x.is_polar:
                    a = self.func(x)
                    if isinstance(a, log):
                        expr.append(self.func(x)._eval_expand_log(**hints))
                    else:
                        expr.append(a)
                elif x.is_negative:
                    a = self.func(-x)
                    expr.append(a)
                    nonpos.append(Integer(-1))
                else:
                    nonpos.append(x)
            return Add(*expr) + log(Mul(*nonpos))
        elif arg.is_Pow:
            if force or (arg.exp.is_extended_real and arg.base.is_positive) or \
                    arg.base.is_polar:
                b = arg.base
                e = arg.exp
                a = self.func(b)
                if isinstance(a, log):
                    return unpolarify(e) * a._eval_expand_log(**hints)
                else:
                    return unpolarify(e) * a
        elif isinstance(arg, Product):
            if arg.function.is_positive:
                return Sum(log(arg.function), *arg.limits)

        return self.func(arg)

    def _eval_simplify(self, ratio, measure):
        from ...simplify import simplify
        if (len(self.args) == 2):
            return simplify(self.func(*self.args), ratio=ratio, measure=measure)
        expr = self.func(simplify(self.args[0], ratio=ratio, measure=measure))
        expr = expand_log(expr, deep=True)
        return min([expr, self], key=measure)

    def as_real_imag(self, deep=True, **hints):
        """
        Returns this function as a complex coordinate.

        Examples
        ========

        >>> log(x).as_real_imag()
        (log(Abs(x)), arg(x))
        >>> log(I).as_real_imag()
        (0, pi/2)
        >>> log(1 + I).as_real_imag()
        (log(sqrt(2)), pi/4)
        >>> log(I*x).as_real_imag()
        (log(Abs(x)), arg(I*x))

        """
        from .complexes import Abs, arg
        if deep:
            abs = Abs(self.args[0].expand(deep, **hints))
            arg = arg(self.args[0].expand(deep, **hints))
        else:
            abs = Abs(self.args[0])
            arg = arg(self.args[0])
        if hints.get('log', False):  # Expand the log
            hints['complex'] = False
            return log(abs).expand(deep, **hints), arg
        else:
            return log(abs), arg

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational and (self.args[0] - 1).is_nonzero:
                return False
        else:
            return s.is_rational

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if self.args[0].is_algebraic and (self.args[0] - 1).is_nonzero:
                return False
        else:
            return s.is_algebraic

    def _eval_is_extended_real(self):
        return self.args[0].is_positive

    def _eval_is_finite(self):
        arg = self.args[0]
        if arg.is_zero:
            return False
        elif arg.is_nonzero:
            return arg.is_finite

    def _eval_is_complex(self):
        arg = self.args[0]
        if arg.is_nonzero:
            return arg.is_complex

    def _eval_is_positive(self):
        return (self.args[0] - 1).is_positive

    def _eval_is_zero(self):
        return (self.args[0] - 1).is_zero

    def _eval_nseries(self, x, n, logx):
        from ...series import Order
        from .complexes import arg
        from .integers import floor
        if not logx:
            logx = log(x)
        arg_series = self.args[0].nseries(x, n=n, logx=logx)
        while arg_series.is_Order:
            n += 1
            arg_series = self.args[0].nseries(x, n=n, logx=logx)
        arg0 = arg_series.as_leading_term(x)
        c, e = arg0.as_coeff_exponent(x)
        t = (arg_series/arg0 - 1).cancel().nseries(x, n=n, logx=logx)
        # series of log(1 + t) in t
        log_series = term = t
        for i in range(1, n):
            term *= -i*t/(i + 1)
            term = term.nseries(x, n=n, logx=logx)
            log_series += term
        if t != 0:
            log_series += Order(t**n, x)
            # branch handling
            if c.is_negative:
                l = floor(arg(t.removeO()*c)/(2*pi)).limit(x, 0)
                if l.is_finite:
                    log_series += 2*I*pi*l
                else:
                    raise NotImplementedError
        return log_series + log(c) + e*logx

    def _eval_as_leading_term(self, x):
        arg = self.args[0].as_leading_term(x)
        if arg == 1:
            return (self.args[0] - 1).as_leading_term(x)
        return self.func(arg)


class LambertW(Function):
    r"""
    The Lambert W function `W(z)` is defined as the inverse
    function of `w \exp(w)`.

    In other words, the value of `W(z)` is such that `z = W(z) \exp(W(z))`
    for any complex number `z`.  The Lambert W function is a multivalued
    function with infinitely many branches `W_k(z)`, indexed by
    `k \in \mathbb{Z}`.  Each branch gives a different solution `w`
    of the equation `z = w \exp(w)`.

    The Lambert W function has two partially real branches: the
    principal branch (`k = 0`) is real for real `z > -1/e`, and the
    `k = -1` branch is real for `-1/e < z < 0`. All branches except
    `k = 0` have a logarithmic singularity at `z = 0`.

    Examples
    ========

    >>> LambertW(1.2)
    0.635564016364870
    >>> LambertW(1.2, -1).evalf()
    -1.34747534407696 - 4.41624341514535*I
    >>> LambertW(-1).is_real
    False

    References
    ==========

    * https://en.wikipedia.org/wiki/Lambert_W_function

    """

    @classmethod
    def eval(cls, x, k=None):
        if k == 0:
            return cls(x)
        elif k is None:
            k = Integer(0)

        if k == 0:
            if x == 0:
                return Integer(0)
            if x is E:
                return Integer(1)
            if x == -1/E:
                return Integer(-1)
            if x == -log(2)/2:
                return -log(2)
            if x is oo:
                return oo

        if k.is_nonzero:
            if x == 0:
                return -oo
        if k == -1:
            if x == -pi/2:
                return -I*pi/2
            elif x == -1/E:
                return Integer(-1)
            elif x == -2*exp(-2):
                return -Integer(2)

    def fdiff(self, argindex=1):
        """Return the first derivative of this function."""
        x = self.args[0]

        if len(self.args) == 1:
            if argindex == 1:
                return LambertW(x)/(x*(1 + LambertW(x)))
            else:
                raise ArgumentIndexError(self, argindex)
        else:
            k = self.args[1]
            if argindex == 1:
                return LambertW(x, k)/(x*(1 + LambertW(x, k)))
            else:
                raise ArgumentIndexError(self, argindex)

    def _eval_is_extended_real(self):
        x = self.args[0]
        if len(self.args) == 1:
            k = Integer(0)
        else:
            k = self.args[1]
        if k.is_zero:
            if (x + 1/E).is_positive:
                return True
            elif (x + 1/E).is_nonpositive:
                return False
        elif (k + 1).is_zero:
            if x.is_negative and (x + 1/E).is_nonnegative:
                return True
            elif x.is_nonpositive or (x + 1/E).is_positive:
                return False

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if self.args[0].is_nonzero and self.args[0].is_algebraic:
                return False
        else:
            return s.is_algebraic

    def _eval_nseries(self, x, n, logx):
        if len(self.args) == 1:
            from ...series import Order
            from .. import factorial
            x = self.args[0]
            o = Order(x**n, x)
            l = Integer(0)
            if n > 0:
                l += Add(*[Integer(-k)**(k - 1)*x**k/factorial(k)
                           for k in range(1, n)])
            return l + o
        return super()._eval_nseries(x, n=n, logx=logx)
