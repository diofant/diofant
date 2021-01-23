import typing

from ...core import Function, I, Integer, Rational, cacheit, nan, oo, pi, zoo
from ...core.function import ArgumentIndexError, _coeff_isneg
from ...core.sympify import sympify
from ..combinatorial.factorials import RisingFactorial, factorial
from .exponential import exp, log
from .miscellaneous import sqrt


def _rewrite_hyperbolics_as_exp(expr):
    expr = sympify(expr)
    return expr.xreplace({h: h.rewrite(exp)
                          for h in expr.atoms(HyperbolicFunction)})


###############################################################################
# ######################### HYPERBOLIC FUNCTIONS ############################ #
###############################################################################


class HyperbolicFunction(Function):
    """
    Base class for hyperbolic functions.

    See Also
    ========

    diofant.functions.elementary.hyperbolic.sinh
    diofant.functions.elementary.hyperbolic.cosh
    diofant.functions.elementary.hyperbolic.tanh
    diofant.functions.elementary.hyperbolic.coth

    """

    unbranched = True


class sinh(HyperbolicFunction):
    r"""
    The hyperbolic sine function, `\frac{e^x - e^{-x}}{2}`.

    * sinh(x) -> Returns the hyperbolic sine of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.cosh
    diofant.functions.elementary.hyperbolic.tanh
    diofant.functions.elementary.hyperbolic.asinh

    """

    def fdiff(self, argindex=1):
        """Returns the first derivative of this function."""
        if argindex == 1:
            return cosh(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return asinh

    @classmethod
    def eval(cls, arg):
        from .trigonometric import sin

        arg = sympify(arg)

        if arg.is_Number:
            if arg in (oo, -oo, 0):
                return arg
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is zoo:
                return nan

            i_coeff = arg.as_coefficient(I)

            if i_coeff is not None:
                return I * sin(i_coeff)
            else:
                if _coeff_isneg(arg):
                    return -cls(-arg)

            if arg.func == asinh:
                return arg.args[0]

            if arg.func == acosh:
                x = arg.args[0]
                return sqrt(x - 1) * sqrt(x + 1)

            if arg.func == atanh:
                x = arg.args[0]
                return x/sqrt(1 - x**2)

            if arg.func == acoth:
                x = arg.args[0]
                return 1/(sqrt(x - 1) * sqrt(x + 1))

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        """Returns the next term in the Taylor series expansion."""
        if n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            if len(previous_terms) >= 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n - 1))
            else:
                return x**n / factorial(n)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        """Returns this function as a complex coordinate."""
        from .trigonometric import cos, sin
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            else:
                return self, Integer(0)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return sinh(re)*cos(im), cosh(re)*sin(im)

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=deep, **hints)
        return re_part + im_part*I

    def _eval_expand_trig(self, **hints):
        arg = self.args[0]
        x = None
        if arg.is_Add:  # TODO, implement more if deep stuff here
            x, y = arg.as_two_terms()
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff != 1 and coeff.is_Integer and terms != 1:
                x = terms
                y = (coeff - 1)*x
        if x is not None:
            return (sinh(x)*cosh(y) + sinh(y)*cosh(x)).expand(trig=True)
        return sinh(arg)

    def _eval_rewrite_as_tractable(self, arg):
        return (exp(arg) - exp(-arg)) / 2

    def _eval_rewrite_as_exp(self, arg):
        return (exp(arg) - exp(-arg)) / 2

    def _eval_rewrite_as_cosh(self, arg):
        return -I*cosh(arg + pi*I/2)

    def _eval_rewrite_as_tanh(self, arg):
        tanh_half = tanh(arg/2)
        return 2*tanh_half/(1 - tanh_half**2)

    def _eval_rewrite_as_coth(self, arg):
        coth_half = coth(arg/2)
        return 2*coth_half/(coth_half**2 - 1)

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def _eval_is_finite(self):
        if self.args[0].is_imaginary:
            return True


class cosh(HyperbolicFunction):
    r"""
    The hyperbolic cosine function, `\frac{e^x + e^{-x}}{2}`.

    * cosh(x) -> Returns the hyperbolic cosine of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.sinh
    diofant.functions.elementary.hyperbolic.tanh
    diofant.functions.elementary.hyperbolic.acosh

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return sinh(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from .trigonometric import cos
        arg = sympify(arg)

        if arg.is_Number:
            if arg in (oo, -oo):
                return oo
            elif arg == 0:
                return Integer(1)
            elif arg.is_negative:
                return cls(-arg)
        else:
            if arg is zoo:
                return nan

            i_coeff = arg.as_coefficient(I)

            if i_coeff is not None:
                return cos(i_coeff)
            else:
                if _coeff_isneg(arg):
                    return cls(-arg)

            if arg.func == asinh:
                return sqrt(1 + arg.args[0]**2)

            if arg.func == acosh:
                return arg.args[0]

            if arg.func == atanh:
                return 1/sqrt(1 - arg.args[0]**2)

            if arg.func == acoth:
                x = arg.args[0]
                return x/(sqrt(x - 1) * sqrt(x + 1))

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 1:
            return Integer(0)
        else:
            x = sympify(x)

            if len(previous_terms) >= 2:
                p = previous_terms[-2]
                return p * x**2 / (n*(n - 1))
            else:
                return x**n/factorial(n)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        from .trigonometric import cos, sin
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            else:
                return self, Integer(0)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()

        return cosh(re)*cos(im), sinh(re)*sin(im)

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=deep, **hints)
        return re_part + im_part*I

    def _eval_expand_trig(self, deep=True, **hints):
        arg = self.args[0]
        x = None
        if arg.is_Add:  # TODO, implement more if deep stuff here
            x, y = arg.as_two_terms()
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff != 1 and coeff.is_Integer and terms != 1:
                x = terms
                y = (coeff - 1)*x
        if x is not None:
            return (cosh(x)*cosh(y) + sinh(x)*sinh(y)).expand(trig=True)
        return cosh(arg)

    def _eval_rewrite_as_tractable(self, arg):
        return (exp(arg) + exp(-arg)) / 2

    def _eval_rewrite_as_exp(self, arg):
        return (exp(arg) + exp(-arg)) / 2

    def _eval_rewrite_as_sinh(self, arg):
        return -I*sinh(arg + pi*I/2)

    def _eval_rewrite_as_tanh(self, arg):
        tanh_half = tanh(arg/2)**2
        return (1 + tanh_half)/(1 - tanh_half)

    def _eval_rewrite_as_coth(self, arg):
        coth_half = coth(arg/2)**2
        return (coth_half + 1)/(coth_half - 1)

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return Integer(1)
        else:
            return self.func(arg)

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def _eval_is_finite(self):
        if self.args[0].is_imaginary:
            return True


class tanh(HyperbolicFunction):
    r"""
    The hyperbolic tangent function, `\frac{\sinh(x)}{\cosh(x)}`.

    * tanh(x) -> Returns the hyperbolic tangent of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.sinh
    diofant.functions.elementary.hyperbolic.cosh
    diofant.functions.elementary.hyperbolic.atanh

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1 - tanh(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return atanh

    @classmethod
    def eval(cls, arg):
        from .trigonometric import tan
        arg = sympify(arg)

        if arg.is_Number:
            if arg is oo:
                return Integer(1)
            elif arg == -oo:
                return Integer(-1)
            elif arg == 0:
                return Integer(0)
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is zoo:
                return nan

            i_coeff = arg.as_coefficient(I)

            if i_coeff is not None:
                if _coeff_isneg(i_coeff):
                    return -I * tan(-i_coeff)
                return I * tan(i_coeff)
            else:
                if _coeff_isneg(arg):
                    return -cls(-arg)

            if arg.func == asinh:
                x = arg.args[0]
                return x/sqrt(1 + x**2)

            if arg.func == acosh:
                x = arg.args[0]
                return sqrt(x - 1) * sqrt(x + 1) / x

            if arg.func == atanh:
                return arg.args[0]

            if arg.func == acoth:
                return 1/arg.args[0]

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from .. import bernoulli
        if n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)

            a = 2**(n + 1)

            B = bernoulli(n + 1)
            F = factorial(n + 1)

            return a*(a - 1) * B/F * x**n

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        from .trigonometric import cos, sin
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            else:
                return self, Integer(0)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = sinh(re)**2 + cos(im)**2
        return sinh(re)*cosh(re)/denom, sin(im)*cos(im)/denom

    def _eval_rewrite_as_tractable(self, arg):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp - neg_exp)/(pos_exp + neg_exp)

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp - neg_exp)/(pos_exp + neg_exp)

    def _eval_rewrite_as_sinh(self, arg):
        return I*sinh(arg)/sinh(pi*I/2 - arg)

    def _eval_rewrite_as_cosh(self, arg):
        return I*cosh(pi*I/2 - arg)/cosh(arg)

    def _eval_rewrite_as_coth(self, arg):
        return 1/coth(arg)

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def _eval_is_finite(self):
        if self.args[0].is_extended_real:
            return True


class coth(HyperbolicFunction):
    r"""
    The hyperbolic cotangent function, `\frac{\cosh(x)}{\sinh(x)}`.

    * coth(x) -> Returns the hyperbolic cotangent of x

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1/sinh(self.args[0])**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return acoth

    @classmethod
    def eval(cls, arg):
        from .trigonometric import cot
        arg = sympify(arg)

        if arg.is_Number:
            if arg is oo:
                return Integer(1)
            elif arg == -oo:
                return Integer(-1)
            elif arg == 0:
                return zoo
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is zoo:
                return nan

            i_coeff = arg.as_coefficient(I)

            if i_coeff is not None:
                if _coeff_isneg(i_coeff):
                    return I * cot(-i_coeff)
                return -I * cot(i_coeff)
            else:
                if _coeff_isneg(arg):
                    return -cls(-arg)

            if arg.func == asinh:
                x = arg.args[0]
                return sqrt(1 + x**2)/x

            if arg.func == acosh:
                x = arg.args[0]
                return x/(sqrt(x - 1) * sqrt(x + 1))

            if arg.func == atanh:
                return 1/arg.args[0]

            if arg.func == acoth:
                return arg.args[0]

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from .. import bernoulli
        if n == 0:
            return 1 / sympify(x)
        elif n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)

            B = bernoulli(n + 1)
            F = factorial(n + 1)

            return 2**(n + 1) * B/F * x**n

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        from .trigonometric import cos, sin
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            else:
                return self, Integer(0)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        denom = sinh(re)**2 + sin(im)**2
        return sinh(re)*cosh(re)/denom, -sin(im)*cos(im)/denom

    def _eval_rewrite_as_tractable(self, arg):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp + neg_exp)/(pos_exp - neg_exp)

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = exp(-arg), exp(arg)
        return (pos_exp + neg_exp)/(pos_exp - neg_exp)

    def _eval_rewrite_as_sinh(self, arg):
        return -I*sinh(pi*I/2 - arg)/sinh(arg)

    def _eval_rewrite_as_cosh(self, arg):
        return -I*cosh(arg)/cosh(pi*I/2 - arg)

    def _eval_rewrite_as_tanh(self, arg):
        return 1/tanh(arg)

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return 1/arg
        else:
            return self.func(arg)


class ReciprocalHyperbolicFunction(HyperbolicFunction):
    """Base class for reciprocal functions of hyperbolic functions."""

    # To be defined in class
    _reciprocal_of = None
    _is_even: typing.Optional[bool] = None
    _is_odd: typing.Optional[bool] = None

    @classmethod
    def eval(cls, arg):
        if arg.could_extract_minus_sign():
            if cls._is_even:
                return cls(-arg)
            elif cls._is_odd:
                return -cls(-arg)

        t = cls._reciprocal_of.eval(arg)
        return 1/t if t is not None else t

    def _call_reciprocal(self, method_name, *args, **kwargs):
        # Calls method_name on _reciprocal_of
        o = self._reciprocal_of(self.args[0])
        return getattr(o, method_name)(*args, **kwargs)

    def _rewrite_reciprocal(self, method_name, arg):
        # Special handling for rewrite functions. If reciprocal rewrite returns
        # unmodified expression, then return None
        t = self._call_reciprocal(method_name, arg)
        assert t is not None and t != self._reciprocal_of(arg)
        return 1/t

    def _eval_rewrite_as_exp(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_exp', arg)

    def _eval_rewrite_as_tractable(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_tractable', arg)

    def _eval_rewrite_as_tanh(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_tanh', arg)

    def _eval_rewrite_as_coth(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_coth', arg)

    def as_real_imag(self, deep=True, **hints):
        return (1 / self._reciprocal_of(self.args[0])).as_real_imag(deep, **hints)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=True, **hints)
        return re_part + I*im_part

    def _eval_as_leading_term(self, x):
        return (1/self._reciprocal_of(self.args[0]))._eval_as_leading_term(x)

    def _eval_is_extended_real(self):
        return self._reciprocal_of(self.args[0]).is_extended_real

    def _eval_is_finite(self):
        return (1/self._reciprocal_of(self.args[0])).is_finite


class csch(ReciprocalHyperbolicFunction):
    r"""
    The hyperbolic cosecant function, `\frac{2}{e^x - e^{-x}}`

    * csch(x) -> Returns the hyperbolic cosecant of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.sinh
    diofant.functions.elementary.hyperbolic.cosh
    diofant.functions.elementary.hyperbolic.tanh
    diofant.functions.elementary.hyperbolic.sech
    diofant.functions.elementary.hyperbolic.asinh
    diofant.functions.elementary.hyperbolic.acosh

    """

    _reciprocal_of = sinh
    _is_odd = True

    def fdiff(self, argindex=1):
        """Returns the first derivative of this function."""
        if argindex == 1:
            return -coth(self.args[0]) * csch(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        """Returns the next term in the Taylor series expansion."""
        from .. import bernoulli
        if n == 0:
            return 1/sympify(x)
        elif n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)

            B = bernoulli(n + 1)
            F = factorial(n + 1)

            return 2 * (1 - 2**n) * B/F * x**n

    def _eval_rewrite_as_cosh(self, arg):
        return I / cosh(arg + I * pi / 2)


class sech(ReciprocalHyperbolicFunction):
    r"""
    The hyperbolic secant function, `\frac{2}{e^x + e^{-x}}`

    * sech(x) -> Returns the hyperbolic secant of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.sinh
    diofant.functions.elementary.hyperbolic.cosh
    diofant.functions.elementary.hyperbolic.tanh
    diofant.functions.elementary.hyperbolic.coth
    diofant.functions.elementary.hyperbolic.csch
    diofant.functions.elementary.hyperbolic.asinh
    diofant.functions.elementary.hyperbolic.acosh

    """

    _reciprocal_of = cosh
    _is_even = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return - tanh(self.args[0])*sech(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from ..combinatorial.numbers import euler
        if n < 0 or n % 2 == 1:
            return Integer(0)
        else:
            x = sympify(x)
            return euler(n) / factorial(n) * x**n

    def _eval_rewrite_as_sinh(self, arg):
        return I / sinh(arg + I * pi / 2)


###############################################################################
# ########################### HYPERBOLIC INVERSES ########################### #
###############################################################################

class asinh(Function):
    """
    The inverse hyperbolic sine function.

    * asinh(x) -> Returns the inverse hyperbolic sine of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.cosh
    diofant.functions.elementary.hyperbolic.tanh
    diofant.functions.elementary.hyperbolic.sinh

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sqrt(self.args[0]**2 + 1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from .trigonometric import asin
        arg = sympify(arg)

        if arg.is_Number:
            if arg in (oo, -oo, 0):
                return arg
            elif arg == 1:
                return log(sqrt(2) + 1)
            elif arg == -1:
                return log(sqrt(2) - 1)
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is zoo:
                return zoo

            i_coeff = arg.as_coefficient(I)

            if i_coeff is not None:
                return I * asin(i_coeff)
            else:
                if _coeff_isneg(arg):
                    return -cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            if len(previous_terms) >= 2 and n > 2:
                p = previous_terms[-2]
                return -p * (n - 2)**2/(n*(n - 1)) * x**2
            else:
                k = (n - 1) // 2
                R = RisingFactorial(Rational(1, 2), k)
                F = factorial(k)
                return (-1)**k * R / F * x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_rewrite_as_log(self, x):
        return log(x + sqrt(x**2 + 1))

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return sinh


class acosh(Function):
    """
    The inverse hyperbolic cosine function.

    * acosh(x) -> Returns the inverse hyperbolic cosine of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.asinh
    diofant.functions.elementary.hyperbolic.atanh
    diofant.functions.elementary.hyperbolic.cosh

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sqrt(self.args[0]**2 - 1)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        arg = sympify(arg)

        if arg.is_Number:
            if arg in (oo, -oo):
                return oo
            elif arg == 0:
                return pi*I / 2
            elif arg == 1:
                return Integer(0)
            elif arg == -1:
                return pi*I

        if arg.is_number:
            cst_table = {
                I: log(I*(1 + sqrt(2))),
                -I: log(-I*(1 + sqrt(2))),
                Rational(+1, 2): pi/3,
                Rational(-1, 2): 2*pi/3,
                sqrt(2)/2: pi/4,
                -sqrt(2)/2: 3*pi/4,
                1/sqrt(2): pi/4,
                -1/sqrt(2): 3*pi/4,
                sqrt(3)/2: pi/6,
                -sqrt(3)/2: 5*pi/6,
                (sqrt(3) - 1)/sqrt(2**3): 5*pi/12,
                -(sqrt(3) - 1)/sqrt(2**3): 7*pi/12,
                sqrt(2 + sqrt(2))/2: pi/8,
                -sqrt(2 + sqrt(2))/2: 7*pi/8,
                sqrt(2 - sqrt(2))/2: 3*pi/8,
                -sqrt(2 - sqrt(2))/2: 5*pi/8,
                (1 + sqrt(3))/(2*sqrt(2)): pi/12,
                -(1 + sqrt(3))/(2*sqrt(2)): 11*pi/12,
                (sqrt(5) + 1)/4: pi/5,
                -(sqrt(5) + 1)/4: 4*pi/5
            }

            if arg in cst_table:
                if arg.is_extended_real:
                    return cst_table[arg]*I
                return cst_table[arg]

        if arg.is_infinite:
            return oo

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return pi*I / 2
        elif n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            if len(previous_terms) >= 2 and n > 2:
                p = previous_terms[-2]
                return p * (n - 2)**2/(n*(n - 1)) * x**2
            else:
                k = (n - 1) // 2
                R = RisingFactorial(Rational(1, 2), k)
                F = factorial(k)
                return -R / F * I * x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return I*pi/2
        else:
            return self.func(arg)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return cosh

    def _eval_rewrite_as_log(self, x):
        return log(x + sqrt(x - 1)*sqrt(x + 1))

    def _eval_nseries(self, x, n, logx):
        x0 = self.args[0].limit(x, 0)
        if x0 == 1:
            return self._eval_rewrite_as_log(self.args[0])._eval_nseries(x, n, logx)
        else:
            return super()._eval_nseries(x, n, logx)


class atanh(Function):
    """
    The inverse hyperbolic tangent function.

    * atanh(x) -> Returns the inverse hyperbolic tangent of x

    See Also
    ========

    diofant.functions.elementary.hyperbolic.asinh
    diofant.functions.elementary.hyperbolic.acosh
    diofant.functions.elementary.hyperbolic.tanh

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from .trigonometric import atan
        arg = sympify(arg)

        if arg.is_Number:
            if arg == 0:
                return Integer(0)
            elif arg == 1:
                return oo
            elif arg == -1:
                return -oo
            elif arg is oo:
                return -I * atan(arg)
            elif arg == -oo:
                return I * atan(-arg)
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is zoo:
                return nan

            i_coeff = arg.as_coefficient(I)

            if i_coeff is not None:
                return I * atan(i_coeff)
            else:
                if _coeff_isneg(arg):
                    return -cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            return x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return tanh


class acoth(Function):
    """
    The inverse hyperbolic cotangent function.

    * acoth(x) -> Returns the inverse hyperbolic cotangent of x

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from .trigonometric import acot
        arg = sympify(arg)

        if arg.is_Number:
            if arg in (oo, -oo):
                return Integer(0)
            elif arg == 0:
                return pi*I / 2
            elif arg == 1:
                return oo
            elif arg == -1:
                return -oo
            elif arg.is_negative:
                return -cls(-arg)
        else:
            if arg is zoo:
                return 0

            i_coeff = arg.as_coefficient(I)

            if i_coeff is not None:
                return -I * acot(i_coeff)
            else:
                if _coeff_isneg(arg):
                    return -cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return pi*I / 2
        elif n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            return x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return I*pi/2
        else:
            return self.func(arg)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return coth
