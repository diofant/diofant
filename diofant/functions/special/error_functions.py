"""This module contains various functions that are special cases
of incomplete gamma functions. It should probably be renamed.
"""

from ...core import (Add, EulerGamma, Function, I, Integer, Pow, Rational,
                     cacheit, expand_mul, oo, pi, zoo)
from ...core.function import ArgumentIndexError
from ...core.sympify import sympify
from ..combinatorial.factorials import factorial
from ..elementary.complexes import polar_lift
from ..elementary.exponential import exp, log
from ..elementary.hyperbolic import cosh, sinh
from ..elementary.integers import floor
from ..elementary.miscellaneous import root, sqrt
from ..elementary.trigonometric import cos, sin
from .hyper import hyper, meijerg


# TODO series expansions
# TODO see the "Note:" in Ei

###############################################################################
# ############################## ERROR FUNCTION ############################# #
###############################################################################


class erf(Function):
    r"""
    The Gauss error function. This function is defined as:

    .. math ::
        \mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \mathrm{d}t.

    Examples
    ========

    Several special values are known:

    >>> erf(0)
    0
    >>> erf(oo)
    1
    >>> erf(-oo)
    -1
    >>> erf(I*oo)
    oo*I
    >>> erf(-I*oo)
    -oo*I

    In general one can pull out factors of -1 and I from the argument:

    >>> erf(-z)
    -erf(z)

    The error function obeys the mirror symmetry:

    >>> conjugate(erf(z))
    erf(conjugate(z))

    Differentiation with respect to z is supported:

    >>> diff(erf(z), z)
    2*E**(-z**2)/sqrt(pi)

    We can numerically evaluate the error function to arbitrary precision
    on the whole complex plane:

    >>> erf(4).evalf(30)
    0.999999984582742099719981147840

    >>> erf(-4*I).evalf(30)
    -1296959.73071763923152794095062*I

    See Also
    ========

    erfc: Complementary error function.
    erfi: Imaginary error function.
    erf2: Two-argument error function.
    erfinv: Inverse error function.
    erfcinv: Inverse Complementary error function.
    erf2inv: Inverse two-argument error function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Error_function
    * https://dlmf.nist.gov/7
    * https://mathworld.wolfram.com/Erf.html
    * http://functions.wolfram.com/GammaBetaErf/Erf

    """

    unbranched = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 2*exp(-self.args[0]**2)/sqrt(pi)
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return erfinv

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is oo:
                return Integer(1)
            if arg == -oo:
                return Integer(-1)
            if arg == 0:
                return Integer(0)

        if isinstance(arg, erfinv):
            return arg.args[0]

        if isinstance(arg, erfcinv):
            return 1 - arg.args[0]

        # Try to pull out factors of I
        t = arg.as_coefficient(I)
        if t in (oo, -oo):
            return arg

        # Try to pull out factors of -1
        if arg.could_extract_minus_sign():
            return -cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return Integer(0)
        x = sympify(x)
        k = floor(Rational(n - 1, 2))
        if len(previous_terms) >= 2:
            return -previous_terms[-2] * x**2 * (n - 2)/(n*k)
        return 2*(-1)**k * x**n/(n*factorial(k)*sqrt(pi))

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_is_extended_real(self):
        arg = self.args[0]
        if arg.is_extended_real:
            return True
        if arg.is_imaginary and arg.is_nonzero:
            return False

    def _eval_rewrite_as_uppergamma(self, z):
        from .gamma_functions import uppergamma
        return sqrt(z**2)/z*(1 - uppergamma(Rational(1, 2), z**2)/sqrt(pi))

    def _eval_rewrite_as_fresnels(self, z):
        arg = (1 - I)*z/sqrt(pi)
        return (1 + I)*(fresnelc(arg) - I*fresnels(arg))

    def _eval_rewrite_as_fresnelc(self, z):
        arg = (1 - I)*z/sqrt(pi)
        return (1 + I)*(fresnelc(arg) - I*fresnels(arg))

    def _eval_rewrite_as_meijerg(self, z):
        return z/sqrt(pi)*meijerg([Rational(1, 2)], [], [0], [-Rational(1, 2)], z**2)

    def _eval_rewrite_as_hyper(self, z):
        return 2*z/sqrt(pi)*hyper([Rational(1, 2)], [Rational(3, 2)], -z**2)

    def _eval_rewrite_as_expint(self, z):
        return sqrt(z**2)/z - z*expint(Rational(1, 2), z**2)/sqrt(pi)

    def _eval_rewrite_as_tractable(self, z, wrt=None, **kwargs):
        if wrt is not None and z.limit(wrt, oo) == -oo:
            return -1 + _erfs(-z)*exp(-z**2)
        return 1 - _erfs(z)*exp(-z**2)

    def _eval_rewrite_as_erfc(self, z):
        return 1 - erfc(z)

    def _eval_rewrite_as_erfi(self, z):
        return -I*erfi(I*z)

    def _eval_as_leading_term(self, x):
        from ...calculus import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return 2*x/sqrt(pi)
        return self.func(arg)

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            return self, Integer(0)
        if deep:
            x, y = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            x, y = self.args[0].as_real_imag()

        if x.is_zero:
            re = Integer(0)
            im = erfi(y)
        else:
            sq = -y**2/x**2
            re = (self.func(x + x*sqrt(sq)) + self.func(x - x*sqrt(sq)))/2
            im = x/(2*y)*sqrt(sq)*(self.func(x - x*sqrt(sq)) - self.func(x + x*sqrt(sq)))
        return re, im


class erfc(Function):
    r"""
    Complementary Error Function. The function is defined as:

    .. math ::
        \mathrm{erfc}(x) = \frac{2}{\sqrt{\pi}} \int_x^\infty e^{-t^2} \mathrm{d}t

    Examples
    ========

    Several special values are known:

    >>> erfc(0)
    1
    >>> erfc(oo)
    0
    >>> erfc(-oo)
    2
    >>> erfc(I*oo)
    -oo*I
    >>> erfc(-I*oo)
    oo*I

    The error function obeys the mirror symmetry:

    >>> conjugate(erfc(z))
    erfc(conjugate(z))

    Differentiation with respect to z is supported:

    >>> diff(erfc(z), z)
    -2*E**(-z**2)/sqrt(pi)

    It also follows

    >>> erfc(-z)
    -erfc(z) + 2

    We can numerically evaluate the complementary error function to arbitrary precision
    on the whole complex plane:

    >>> erfc(4).evalf(30)
    0.0000000154172579002800188521596734869

    >>> erfc(4*I).evalf(30)
    1.0 - 1296959.73071763923152794095062*I

    See Also
    ========

    erf: Gaussian error function.
    erfi: Imaginary error function.
    erf2: Two-argument error function.
    erfinv: Inverse error function.
    erfcinv: Inverse Complementary error function.
    erf2inv: Inverse two-argument error function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Error_function
    * https://dlmf.nist.gov/7
    * https://mathworld.wolfram.com/Erfc.html
    * http://functions.wolfram.com/GammaBetaErf/Erfc

    """

    unbranched = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -2*exp(-self.args[0]**2)/sqrt(pi)
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return erfcinv

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is oo:
                return Integer(0)
            if arg == 0:
                return Integer(1)

        if isinstance(arg, erfinv):
            return 1 - arg.args[0]

        if isinstance(arg, erfcinv):
            return arg.args[0]

        # Try to pull out factors of I
        t = arg.as_coefficient(I)
        if t in (oo, -oo):
            return -arg

        # Try to pull out factors of -1
        if arg.could_extract_minus_sign():
            return Integer(2) - cls(-arg)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return Integer(1)
        if n < 0 or n % 2 == 0:
            return Integer(0)
        x = sympify(x)
        k = floor(Rational(n - 1, 2))
        if len(previous_terms) >= 2:
            return -previous_terms[-2] * x**2 * (n - 2)/(n*k)
        return -2*(-1)**k * x**n/(n*factorial(k)*sqrt(pi))

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_is_extended_real(self):
        arg = self.args[0]
        if arg.is_extended_real:
            return True
        if arg.is_imaginary and arg.is_nonzero:
            return False

    def _eval_rewrite_as_tractable(self, z, **kwargs):
        return self.rewrite(erf).rewrite('tractable', **kwargs)

    def _eval_rewrite_as_erf(self, z):
        return 1 - erf(z)

    def _eval_rewrite_as_erfi(self, z):
        return 1 + I*erfi(I*z)

    def _eval_rewrite_as_fresnels(self, z):
        arg = (1 - I)*z/sqrt(pi)
        return 1 - (1 + I)*(fresnelc(arg) - I*fresnels(arg))

    def _eval_rewrite_as_fresnelc(self, z):
        arg = (1 - I)*z/sqrt(pi)
        return 1 - (1 + I)*(fresnelc(arg) - I*fresnels(arg))

    def _eval_rewrite_as_meijerg(self, z):
        return 1 - z/sqrt(pi)*meijerg([Rational(1, 2)], [], [0], [-Rational(1, 2)], z**2)

    def _eval_rewrite_as_hyper(self, z):
        return 1 - 2*z/sqrt(pi)*hyper([Rational(1, 2)], [Rational(3, 2)], -z**2)

    def _eval_rewrite_as_uppergamma(self, z):
        from .gamma_functions import uppergamma
        return 1 - sqrt(z**2)/z*(1 - uppergamma(Rational(1, 2), z**2)/sqrt(pi))

    def _eval_rewrite_as_expint(self, z):
        return 1 - sqrt(z**2)/z + z*expint(Rational(1, 2), z**2)/sqrt(pi)

    def _eval_as_leading_term(self, x):
        from ...calculus import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return Integer(1)
        return self.func(arg)

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            return self, Integer(0)
        if deep:
            x, y = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            x, y = self.args[0].as_real_imag()

        if x.is_zero:
            re = Integer(1)
            im = -erfi(y)
        else:
            sq = -y**2/x**2
            re = (self.func(x + x*sqrt(sq)) + self.func(x - x*sqrt(sq)))/2
            im = x/(2*y)*sqrt(sq)*(self.func(x - x*sqrt(sq)) - self.func(x + x*sqrt(sq)))
        return re, im


class erfi(Function):
    r"""
    Imaginary error function. The function erfi is defined as:

    .. math ::
        \mathrm{erfi}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{t^2} \mathrm{d}t

    Examples
    ========

    Several special values are known:

    >>> erfi(0)
    0
    >>> erfi(oo)
    oo
    >>> erfi(-oo)
    -oo
    >>> erfi(I*oo)
    I
    >>> erfi(-I*oo)
    -I

    In general one can pull out factors of -1 and I from the argument:

    >>> erfi(-z)
    -erfi(z)

    >>> conjugate(erfi(z))
    erfi(conjugate(z))

    Differentiation with respect to z is supported:

    >>> diff(erfi(z), z)
    2*E**(z**2)/sqrt(pi)

    We can numerically evaluate the imaginary error function to arbitrary precision
    on the whole complex plane:

    >>> erfi(2).evalf(30)
    18.5648024145755525987042919132

    >>> erfi(-2*I).evalf(30)
    -0.995322265018952734162069256367*I

    See Also
    ========

    erf: Gaussian error function.
    erfc: Complementary error function.
    erf2: Two-argument error function.
    erfinv: Inverse error function.
    erfcinv: Inverse Complementary error function.
    erf2inv: Inverse two-argument error function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Error_function
    * https://mathworld.wolfram.com/Erfi.html
    * http://functions.wolfram.com/GammaBetaErf/Erfi

    """

    unbranched = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 2*exp(self.args[0]**2)/sqrt(pi)
        raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, z):
        if z.is_Number:
            if z == 0:
                return Integer(0)
            if z is oo:
                return oo

        # Try to pull out factors of -1
        if z.could_extract_minus_sign():
            return -cls(-z)

        # Try to pull out factors of I
        nz = z.as_coefficient(I)
        if nz is not None:
            if nz is oo:
                return I
            if isinstance(nz, erfinv):
                return I*nz.args[0]
            if isinstance(nz, erfcinv):
                return I*(1 - nz.args[0])

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return Integer(0)
        x = sympify(x)
        k = floor(Rational(n - 1, 2))
        if len(previous_terms) >= 2:
            return previous_terms[-2] * x**2 * (n - 2)/(n*k)
        return 2 * x**n/(n*factorial(k)*sqrt(pi))

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def _eval_is_extended_real(self):
        arg = self.args[0]
        if arg.is_extended_real:
            return True
        if arg.is_imaginary and arg.is_nonzero:
            return False

    def _eval_rewrite_as_tractable(self, z, **kwargs):
        return self.rewrite(erf).rewrite('tractable', **kwargs)

    def _eval_rewrite_as_erf(self, z):
        return -I*erf(I*z)

    def _eval_rewrite_as_erfc(self, z):
        return I*erfc(I*z) - I

    def _eval_rewrite_as_fresnels(self, z):
        arg = (1 + I)*z/sqrt(pi)
        return (1 - I)*(fresnelc(arg) - I*fresnels(arg))

    def _eval_rewrite_as_fresnelc(self, z):
        arg = (1 + I)*z/sqrt(pi)
        return (1 - I)*(fresnelc(arg) - I*fresnels(arg))

    def _eval_rewrite_as_meijerg(self, z):
        return z/sqrt(pi)*meijerg([Rational(1, 2)], [], [0], [-Rational(1, 2)], -z**2)

    def _eval_rewrite_as_hyper(self, z):
        return 2*z/sqrt(pi)*hyper([Rational(1, 2)], [Rational(3, 2)], z**2)

    def _eval_rewrite_as_uppergamma(self, z):
        from .gamma_functions import uppergamma
        return sqrt(-z**2)/z*(uppergamma(Rational(1, 2), -z**2)/sqrt(pi) - 1)

    def _eval_rewrite_as_expint(self, z):
        return sqrt(-z**2)/z - z*expint(Rational(1, 2), -z**2)/sqrt(pi)

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            return self, Integer(0)
        if deep:
            x, y = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            x, y = self.args[0].as_real_imag()

        if x.is_zero:
            re = Integer(0)
            im = erf(y)
        else:
            sq = -y**2/x**2
            re = (self.func(x + x*sqrt(sq)) + self.func(x - x*sqrt(sq)))/2
            im = x/(2*y)*sqrt(sq)*(self.func(x - x*sqrt(sq)) -
                                   self.func(x + x*sqrt(sq)))
        return re, im


class erf2(Function):
    r"""
    Two-argument error function. This function is defined as:

    .. math ::
        \mathrm{erf2}(x, y) = \frac{2}{\sqrt{\pi}} \int_x^y e^{-t^2} \mathrm{d}t

    Examples
    ========

    Several special values are known:

    >>> erf2(0, 0)
    0
    >>> erf2(x, x)
    0
    >>> erf2(x, oo)
    -erf(x) + 1
    >>> erf2(x, -oo)
    -erf(x) - 1
    >>> erf2(oo, y)
    erf(y) - 1
    >>> erf2(-oo, y)
    erf(y) + 1

    In general one can pull out factors of -1:

    >>> erf2(-x, -y)
    -erf2(x, y)

    The error function obeys the mirror symmetry:

    >>> conjugate(erf2(x, y))
    erf2(conjugate(x), conjugate(y))

    Differentiation with respect to x, y is supported:

    >>> diff(erf2(x, y), x)
    -2*E**(-x**2)/sqrt(pi)
    >>> diff(erf2(x, y), y)
    2*E**(-y**2)/sqrt(pi)

    See Also
    ========

    erf: Gaussian error function.
    erfc: Complementary error function.
    erfi: Imaginary error function.
    erfinv: Inverse error function.
    erfcinv: Inverse Complementary error function.
    erf2inv: Inverse two-argument error function.

    References
    ==========

    * http://functions.wolfram.com/GammaBetaErf/Erf2/

    """

    def fdiff(self, argindex=1):
        x, y = self.args
        if argindex == 1:
            return -2*exp(-x**2)/sqrt(pi)
        if argindex == 2:
            return 2*exp(-y**2)/sqrt(pi)
        raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, x, y):
        I = oo
        N = -oo
        O = Integer(0)
        if x == y:
            return Integer(0)
        if (x is I or x is N or x is O) or (y is I or y is N or y is O):
            return erf(y) - erf(x)

        if isinstance(y, erf2inv) and y.args[0] == x:
            return y.args[1]

        # Try to pull out -1 factor
        sign_x = x.could_extract_minus_sign()
        sign_y = y.could_extract_minus_sign()
        if sign_x and sign_y:
            return -cls(-x, -y)
        if sign_x or sign_y:
            return erf(y)-erf(x)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate(), self.args[1].conjugate())

    def _eval_is_extended_real(self):
        x, y = self.args
        if y.is_extended_real:
            if x.is_extended_real:
                return True
            if x.is_imaginary and x.is_nonzero:
                return False

    def _eval_rewrite_as_erf(self, x, y):
        return erf(y) - erf(x)

    def _eval_rewrite_as_erfc(self, x, y):
        return erfc(x) - erfc(y)

    def _eval_rewrite_as_erfi(self, x, y):
        return I*(erfi(I*x)-erfi(I*y))

    def _eval_rewrite_as_fresnels(self, x, y):
        return erf(y).rewrite(fresnels) - erf(x).rewrite(fresnels)

    def _eval_rewrite_as_fresnelc(self, x, y):
        return erf(y).rewrite(fresnelc) - erf(x).rewrite(fresnelc)

    def _eval_rewrite_as_meijerg(self, x, y):
        return erf(y).rewrite(meijerg) - erf(x).rewrite(meijerg)

    def _eval_rewrite_as_hyper(self, x, y):
        return erf(y).rewrite(hyper) - erf(x).rewrite(hyper)

    def _eval_rewrite_as_uppergamma(self, x, y):
        from .gamma_functions import uppergamma
        return (sqrt(y**2)/y*(1 - uppergamma(Rational(1, 2), y**2)/sqrt(pi)) -
                sqrt(x**2)/x*(1 - uppergamma(Rational(1, 2), x**2)/sqrt(pi)))

    def _eval_rewrite_as_expint(self, x, y):
        return erf(y).rewrite(expint) - erf(x).rewrite(expint)


class erfinv(Function):
    r"""
    Inverse Error Function. The erfinv function is defined as:

    .. math ::
        \mathrm{erf}(x) = y \quad \Rightarrow \quad \mathrm{erfinv}(y) = x

    Examples
    ========

    Several special values are known:

    >>> erfinv(0)
    0
    >>> erfinv(1)
    oo

    Differentiation with respect to x is supported:

    >>> diff(erfinv(x), x)
    E**(erfinv(x)**2)*sqrt(pi)/2

    We can numerically evaluate the inverse error function to arbitrary precision
    on [-1, 1]:

    >>> erfinv(0.2)
    0.179143454621292

    See Also
    ========

    erf: Gaussian error function.
    erfc: Complementary error function.
    erfi: Imaginary error function.
    erf2: Two-argument error function.
    erfcinv: Inverse Complementary error function.
    erf2inv: Inverse two-argument error function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Error_function#Inverse_functions
    * http://functions.wolfram.com/GammaBetaErf/InverseErf/

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return sqrt(pi)*exp(self.func(self.args[0])**2)/2
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return erf

    @classmethod
    def eval(cls, z):
        if z == -1:
            return -oo
        if z == 0:
            return Integer(0)
        if z == 1:
            return oo

        if isinstance(z, erf) and z.args[0].is_extended_real:
            return z.args[0]

        # Try to pull out factors of -1
        nz = z.as_coefficient(-1)
        if isinstance(nz, erf) and nz.args[0].is_extended_real:
            return -nz.args[0]

    def _eval_rewrite_as_erfcinv(self, z):
        return erfcinv(1-z)


class erfcinv(Function):
    r"""
    Inverse Complementary Error Function. The erfcinv function is defined as:

    .. math ::
        \mathrm{erfc}(x) = y \quad \Rightarrow \quad \mathrm{erfcinv}(y) = x

    Examples
    ========

    Several special values are known:

    >>> erfcinv(1)
    0
    >>> erfcinv(0)
    oo

    Differentiation with respect to x is supported:

    >>> diff(erfcinv(x), x)
    -E**(erfcinv(x)**2)*sqrt(pi)/2

    See Also
    ========

    erf: Gaussian error function.
    erfc: Complementary error function.
    erfi: Imaginary error function.
    erf2: Two-argument error function.
    erfinv: Inverse error function.
    erf2inv: Inverse two-argument error function.

    References
    ==========

    * https://en.wikipedia.org/wiki/Error_function#Inverse_functions
    * http://functions.wolfram.com/GammaBetaErf/InverseErfc/

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -sqrt(pi)*exp(self.func(self.args[0])**2)/2
        raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return erfc

    @classmethod
    def eval(cls, z):
        if z == 0:
            return oo
        if z == 1:
            return Integer(0)
        if z == 2:
            return -oo

    def _eval_rewrite_as_erfinv(self, z):
        return erfinv(1-z)


class erf2inv(Function):
    r"""
    Two-argument Inverse error function. The erf2inv function is defined as:

    .. math ::
        \mathrm{erf2}(x, w) = y \quad \Rightarrow \quad \mathrm{erf2inv}(x, y) = w

    Examples
    ========

    Several special values are known:

    >>> erf2inv(0, 0)
    0
    >>> erf2inv(1, 0)
    1
    >>> erf2inv(0, 1)
    oo
    >>> erf2inv(0, y)
    erfinv(y)
    >>> erf2inv(oo, y)
    erfcinv(-y)

    Differentiation with respect to x and y is supported:

    >>> diff(erf2inv(x, y), x)
    E**(-x**2 + erf2inv(x, y)**2)
    >>> diff(erf2inv(x, y), y)
    E**(erf2inv(x, y)**2)*sqrt(pi)/2

    See Also
    ========

    erf: Gaussian error function.
    erfc: Complementary error function.
    erfi: Imaginary error function.
    erf2: Two-argument error function.
    erfinv: Inverse error function.
    erfcinv: Inverse complementary error function.

    References
    ==========

    * http://functions.wolfram.com/GammaBetaErf/InverseErf2/

    """

    def fdiff(self, argindex=1):
        x, y = self.args
        if argindex == 1:
            return exp(self.func(x, y)**2-x**2)
        if argindex == 2:
            return sqrt(pi)*exp(self.func(x, y)**2)/2
        raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, x, y):
        if x == 0 and y == 0:
            return Integer(0)
        if x == 0 and y == 1:
            return oo
        if x == 1 and y == 0:
            return Integer(1)
        if x == 0:
            return erfinv(y)
        if x is oo:
            return erfcinv(-y)
        if y == 0:
            return x
        if y is oo:
            return erfinv(x)


###############################################################################
# ################## EXPONENTIAL INTEGRALS ################################## #
###############################################################################

class Ei(Function):
    r"""
    The classical exponential integral.

    For use in Diofant, this function is defined as

    .. math:: \operatorname{Ei}(x) = \sum_{n=1}^\infty \frac{x^n}{n\, n!}
                                     + \log(x) + \gamma,

    where `\gamma` is the Euler-Mascheroni constant.

    If `x` is a polar number, this defines an analytic function on the
    Riemann surface of the logarithm. Otherwise this defines an analytic
    function in the cut plane `\mathbb{C} \setminus (-\infty, 0]`.

    **Background**

    The name *exponential integral* comes from the following statement:

    .. math:: \operatorname{Ei}(x) = \int_{-\infty}^x \frac{e^t}{t} \mathrm{d}t

    If the integral is interpreted as a Cauchy principal value, this statement
    holds for `x > 0` and `\operatorname{Ei}(x)` as defined above.

    Note that we carefully avoided defining `\operatorname{Ei}(x)` for
    negative real `x`. This is because above integral formula does not hold for
    any polar lift of such `x`, indeed all branches of
    `\operatorname{Ei}(x)` above the negative reals are imaginary.

    However, the following statement holds for all `x \in \mathbb{R}^*`:

    .. math:: \int_{-\infty}^x \frac{e^t}{t} \mathrm{d}t =
              \frac{\operatorname{Ei}\left(|x|e^{i \arg(x)}\right) +
                    \operatorname{Ei}\left(|x|e^{- i \arg(x)}\right)}{2},

    where the integral is again understood to be a principal value if
    `x > 0`, and `|x|e^{i \arg(x)}`,
    `|x|e^{- i \arg(x)}` denote two conjugate polar lifts of `x`.

    Examples
    ========

    The exponential integral in Diofant is strictly undefined for negative values
    of the argument. For convenience, exponential integrals with negative
    arguments are immediately converted into an expression that agrees with
    the classical integral definition:

    >>> Ei(-1)
    -I*pi + Ei(exp_polar(I*pi))

    This yields a real value:

    >>> Ei(-1).evalf(chop=True)
    -0.219383934395520

    On the other hand the analytic continuation is not real:

    >>> Ei(polar_lift(-1)).evalf(chop=True)
    -0.21938393439552 + 3.14159265358979*I

    The exponential integral has a logarithmic branch point at the origin:

    >>> Ei(x*exp_polar(2*I*pi))
    Ei(x) + 2*I*pi

    Differentiation is supported:

    >>> Ei(x).diff(x)
    E**x/x

    The exponential integral is related to many other special functions.
    For example:

    >>> Ei(x).rewrite(expint)
    -expint(1, x*exp_polar(I*pi)) - I*pi
    >>> Ei(x).rewrite(Shi)
    Chi(x) + Shi(x)

    See Also
    ========

    expint: Generalized exponential integral.
    E1: Special case of the generalized exponential integral.
    li: Logarithmic integral.
    Li: Offset logarithmic integral.
    Si: Sine integral.
    Ci: Cosine integral.
    Shi: Hyperbolic sine integral.
    Chi: Hyperbolic cosine integral.
    diofant.functions.special.gamma_functions.uppergamma: Upper incomplete gamma function.

    References
    ==========

    * https://dlmf.nist.gov/6.6
    * https://en.wikipedia.org/wiki/Exponential_integral
    * Abramowitz & Stegun, section 5: http://people.math.sfu.ca/~cbm/aands/page_228.htm

    """

    @classmethod
    def eval(cls, z):
        if z == 0:
            return -oo
        if z is oo:
            return oo
        if z == -oo:
            return Integer(0)

        if not z.is_polar and z.is_negative:
            # Note: is this a good idea?
            return Ei(polar_lift(z)) - pi*I
        nz, n = z.extract_branch_factor()
        if n:
            return Ei(nz) + 2*I*pi*n

    def fdiff(self, argindex=1):
        from .. import unpolarify
        arg = unpolarify(self.args[0])
        if argindex == 1:
            return exp(arg)/arg
        raise ArgumentIndexError(self, argindex)

    def _eval_evalf(self, prec):
        if (self.args[0]/polar_lift(-1)).is_positive:
            return Function._eval_evalf(self, prec) + (I*pi)._eval_evalf(prec)
        return Function._eval_evalf(self, prec)

    def _eval_rewrite_as_uppergamma(self, z):
        from .gamma_functions import uppergamma

        # XXX this does not currently work usefully because uppergamma
        #     immediately turns into expint
        return -uppergamma(0, polar_lift(-1)*z) - I*pi

    def _eval_rewrite_as_expint(self, z, **kwargs):
        return -expint(1, polar_lift(-1)*z) - I*pi

    def _eval_rewrite_as_li(self, z):
        if isinstance(z, log):
            return li(z.args[0])
        # TODO:
        # Actually it only holds that:
        #  Ei(z) = li(exp(z))
        # for -pi < imag(z) <= pi
        return li(exp(z))

    def _eval_rewrite_as_Si(self, z):
        return Shi(z) + Chi(z)
    _eval_rewrite_as_Ci = _eval_rewrite_as_Si
    _eval_rewrite_as_Chi = _eval_rewrite_as_Si
    _eval_rewrite_as_Shi = _eval_rewrite_as_Si

    def _eval_rewrite_as_tractable(self, z, **kwargs):
        return exp(z) * _eis(z)

    def _eval_nseries(self, x, n, logx):
        x0 = self.args[0].limit(x, 0)
        if x0 == 0:
            f = self._eval_rewrite_as_Si(*self.args)
            return f._eval_nseries(x, n, logx)
        return super()._eval_nseries(x, n, logx)


class expint(Function):
    r"""
    Generalized exponential integral.

    This function is defined as

    .. math:: \operatorname{E}_\nu(z) = z^{\nu - 1} \Gamma(1 - \nu, z),

    where `\Gamma(1 - \nu, z)` is the upper incomplete gamma function
    (``uppergamma``).

    Hence for `z` with positive real part we have

    .. math:: \operatorname{E}_\nu(z)
              =   \int_1^\infty \frac{e^{-zt}}{t^\nu} \mathrm{d}t,

    which explains the name.

    The representation as an incomplete gamma function provides an analytic
    continuation for `\operatorname{E}_\nu(z)`. If `\nu` is a
    non-positive integer the exponential integral is thus an unbranched
    function of `z`, otherwise there is a branch point at the origin.
    Refer to the incomplete gamma function documentation for details of the
    branching behavior.

    Examples
    ========

    >>> from diofant.abc import nu

    Differentiation is supported. Differentiation with respect to z explains
    further the name: for integral orders, the exponential integral is an
    iterated integral of the exponential function.

    >>> expint(nu, z).diff(z)
    -expint(nu - 1, z)

    Differentiation with respect to nu has no classical expression:

    >>> expint(nu, z).diff(nu)
    -z**(nu - 1)*meijerg(((), (1, 1)), ((0, 0, -nu + 1), ()), z)

    At non-postive integer orders, the exponential integral reduces to the
    exponential function:

    >>> expint(0, z)
    E**(-z)/z
    >>> expint(-1, z)
    E**(-z)/z + E**(-z)/z**2

    At half-integers it reduces to error functions:

    >>> expint(Rational(1, 2), z)
    -sqrt(pi)*erf(sqrt(z))/sqrt(z) + sqrt(pi)/sqrt(z)

    At positive integer orders it can be rewritten in terms of exponentials
    and expint(1, z). Use expand_func() to do this:

    >>> expand_func(expint(5, z))
    z**4*expint(1, z)/24 + E**(-z)*(-z**3 + z**2 - 2*z + 6)/24

    The generalized exponential integral is essentially equivalent to the
    incomplete gamma function:

    >>> expint(nu, z).rewrite(uppergamma)
    z**(nu - 1)*uppergamma(-nu + 1, z)

    As such it is branched at the origin:

    >>> expint(4, z*exp_polar(2*pi*I))
    I*pi*z**3/3 + expint(4, z)
    >>> expint(nu, z*exp_polar(2*pi*I))
    z**(nu - 1)*(E**(2*I*pi*nu) - 1)*gamma(-nu + 1) + expint(nu, z)

    See Also
    ========

    Ei: Another related function called exponential integral.
    E1: The classical case, returns expint(1, z).
    li: Logarithmic integral.
    Li: Offset logarithmic integral.
    Si: Sine integral.
    Ci: Cosine integral.
    Shi: Hyperbolic sine integral.
    Chi: Hyperbolic cosine integral.
    diofant.functions.special.gamma_functions.uppergamma

    References
    ==========

    * https://dlmf.nist.gov/8.19
    * http://functions.wolfram.com/GammaBetaErf/ExpIntegralE/
    * https://en.wikipedia.org/wiki/Exponential_integral

    """

    @classmethod
    def eval(cls, nu, z):
        from .. import exp, factorial, gamma, unpolarify, uppergamma
        nu2 = unpolarify(nu)
        if nu != nu2:
            return expint(nu2, z)
        if nu.is_Integer and nu <= 0 or (not nu.is_Integer and (2*nu).is_Integer):
            return unpolarify(expand_mul(z**(nu - 1)*uppergamma(1 - nu, z)))

        # Extract branching information. This can be deduced from what is
        # explained in lowergamma.eval().
        z, n = z.extract_branch_factor()
        if n == 0:
            return
        if nu.is_integer:
            if nu.is_positive:
                return expint(nu, z) \
                    - 2*pi*I*n*(-1)**(nu - 1)/factorial(nu - 1)*unpolarify(z)**(nu - 1)
        else:
            return (exp(2*I*pi*nu*n) - 1)*z**(nu - 1)*gamma(1 - nu) + expint(nu, z)

    def fdiff(self, argindex=1):
        from .hyper import meijerg
        nu, z = self.args
        if argindex == 1:
            return -z**(nu - 1)*meijerg([], [1, 1], [0, 0, 1 - nu], [], z)
        if argindex == 2:
            return -expint(nu - 1, z)
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_uppergamma(self, nu, z):
        from .gamma_functions import uppergamma
        return z**(nu - 1)*uppergamma(1 - nu, z)

    def _eval_rewrite_as_Ei(self, nu, z):
        from .. import exp, exp_polar, factorial, unpolarify
        if nu == 1:
            return -Ei(z*exp_polar(-I*pi)) - I*pi
        if nu.is_Integer and nu > 1:
            # DLMF, 8.19.7
            x = -unpolarify(z)
            return x**(nu - 1)/factorial(nu - 1)*E1(z).rewrite(Ei) + \
                exp(x)/factorial(nu - 1) * \
                Add(*[factorial(nu - k - 2)*x**k for k in range(nu - 1)])
        return self

    def _eval_expand_func(self, **hints):
        return self.rewrite(Ei).rewrite(expint, **hints)

    def _eval_rewrite_as_Si(self, nu, z):
        if nu != 1:
            return self
        return Shi(z) - Chi(z)
    _eval_rewrite_as_Ci = _eval_rewrite_as_Si
    _eval_rewrite_as_Chi = _eval_rewrite_as_Si
    _eval_rewrite_as_Shi = _eval_rewrite_as_Si

    def _eval_nseries(self, x, n, logx):
        nu = self.args[0]
        if not nu.has(x) and nu.is_Integer and nu.is_positive:
            f = self._eval_rewrite_as_Ei(*self.args)
            return f._eval_nseries(x, n, logx)
        return super()._eval_nseries(x, n, logx)


def E1(z):
    """
    Classical case of the generalized exponential integral.

    This is equivalent to ``expint(1, z)``.

    See Also
    ========

    Ei: Exponential integral.
    expint: Generalized exponential integral.
    li: Logarithmic integral.
    Li: Offset logarithmic integral.
    Si: Sine integral.
    Ci: Cosine integral.
    Shi: Hyperbolic sine integral.
    Chi: Hyperbolic cosine integral.

    """
    return expint(1, z)


class li(Function):
    r"""
    The classical logarithmic integral.

    For the use in Diofant, this function is defined as

    .. math:: \operatorname{li}(x) = \int_0^x \frac{1}{\log(t)} \mathrm{d}t \,.

    Examples
    ========

    Several special values are known:

    >>> li(0)
    0
    >>> li(1)
    -oo
    >>> li(oo)
    oo

    Differentiation with respect to z is supported:

    >>> diff(li(z), z)
    1/log(z)

    Defining the `li` function via an integral:


    The logarithmic integral can also be defined in terms of Ei:

    >>> li(z).rewrite(Ei)
    Ei(log(z))
    >>> diff(li(z).rewrite(Ei), z)
    1/log(z)

    We can numerically evaluate the logarithmic integral to arbitrary precision
    on the whole complex plane (except the singular points):

    >>> li(2).evalf(30)
    1.04516378011749278484458888919

    >>> li(2*I).evalf(30)
    1.0652795784357498247001125598 + 3.08346052231061726610939702133*I

    We can even compute Soldner's constant by the help of mpmath:

    >>> from mpmath import findroot
    >>> print(findroot(li, 2))
    1.45136923488338

    Further transformations include rewriting `li` in terms of
    the trigonometric integrals `Si`, `Ci`, `Shi` and `Chi`:

    >>> li(z).rewrite(Si)
    -log(I*log(z)) - log(1/log(z))/2 + log(log(z))/2 + Ci(I*log(z)) + Shi(log(z))
    >>> li(z).rewrite(Ci)
    -log(I*log(z)) - log(1/log(z))/2 + log(log(z))/2 + Ci(I*log(z)) + Shi(log(z))
    >>> li(z).rewrite(Shi)
    -log(1/log(z))/2 + log(log(z))/2 + Chi(log(z)) - Shi(log(z))
    >>> li(z).rewrite(Chi)
    -log(1/log(z))/2 + log(log(z))/2 + Chi(log(z)) - Shi(log(z))

    See Also
    ========

    Li: Offset logarithmic integral.
    Ei: Exponential integral.
    expint: Generalized exponential integral.
    E1: Special case of the generalized exponential integral.
    Si: Sine integral.
    Ci: Cosine integral.
    Shi: Hyperbolic sine integral.
    Chi: Hyperbolic cosine integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Logarithmic_integral
    * https://mathworld.wolfram.com/LogarithmicIntegral.html
    * https://dlmf.nist.gov/6
    * https://mathworld.wolfram.com/SoldnersConstant.html

    """

    @classmethod
    def eval(cls, z):
        if z == 0:
            return Integer(0)
        if z == 1:
            return -oo
        if z is oo:
            return oo

    def fdiff(self, argindex=1):
        arg = self.args[0]
        if argindex == 1:
            return 1/log(arg)
        raise ArgumentIndexError(self, argindex)

    def _eval_conjugate(self):
        z = self.args[0]
        # Exclude values on the branch cut (-oo, 0)
        if not (z.is_extended_real and z.is_negative):
            return self.func(z.conjugate())

    def _eval_rewrite_as_Li(self, z):
        return Li(z) + li(2)

    def _eval_rewrite_as_Ei(self, z):
        return Ei(log(z))

    def _eval_rewrite_as_uppergamma(self, z):
        from .gamma_functions import uppergamma
        return (-uppergamma(0, -log(z)) +
                (log(log(z)) - log(1/log(z)))/2 - log(-log(z)))

    def _eval_rewrite_as_Si(self, z):
        return (Ci(I*log(z)) - I*Si(I*log(z)) -
                (log(1/log(z)) - log(log(z)))/2 - log(I*log(z)))

    _eval_rewrite_as_Ci = _eval_rewrite_as_Si

    def _eval_rewrite_as_Shi(self, z):
        return (Chi(log(z)) - Shi(log(z)) - (log(1/log(z)) - log(log(z)))/2)

    _eval_rewrite_as_Chi = _eval_rewrite_as_Shi

    def _eval_rewrite_as_hyper(self, z):
        return (log(z)*hyper((1, 1), (2, 2), log(z)) +
                (log(log(z)) - log(1/log(z)))/2 + EulerGamma)

    def _eval_rewrite_as_meijerg(self, z):
        return (-log(-log(z)) - (log(1/log(z)) - log(log(z)))/2
                - meijerg(((), (1,)), ((0, 0), ()), -log(z)))

    def _eval_rewrite_as_tractable(self, z, **kwargs):
        return z * _eis(log(z))


class Li(Function):
    r"""
    The offset logarithmic integral.

    For the use in Diofant, this function is defined as

    .. math:: \operatorname{Li}(x) = \operatorname{li}(x) - \operatorname{li}(2)

    Examples
    ========

    The following special value is known:

    >>> Li(2)
    0

    Differentiation with respect to z is supported:

    >>> diff(Li(z), z)
    1/log(z)

    The shifted logarithmic integral can be written in terms of `li(z)`:

    >>> Li(z).rewrite(li)
    li(z) - li(2)

    We can numerically evaluate the logarithmic integral to arbitrary precision
    on the whole complex plane (except the singular points):

    >>> Li(2).evalf(30)
    0

    >>> Li(4).evalf(30)
    1.92242131492155809316615998938

    See Also
    ========

    li: Logarithmic integral.
    Ei: Exponential integral.
    expint: Generalized exponential integral.
    E1: Special case of the generalized exponential integral.
    Si: Sine integral.
    Ci: Cosine integral.
    Shi: Hyperbolic sine integral.
    Chi: Hyperbolic cosine integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Logarithmic_integral
    * https://mathworld.wolfram.com/LogarithmicIntegral.html
    * https://dlmf.nist.gov/6

    """

    @classmethod
    def eval(cls, z):
        if z is oo:
            return oo
        if z == 2:
            return Integer(0)

    def fdiff(self, argindex=1):
        arg = self.args[0]
        if argindex == 1:
            return 1/log(arg)
        raise ArgumentIndexError(self, argindex)

    def _eval_evalf(self, prec):
        return self.rewrite(li).evalf(prec)

    def _eval_rewrite_as_li(self, z):
        return li(z) - li(2)

    def _eval_rewrite_as_tractable(self, z, **kwargs):
        return self.rewrite(li).rewrite('tractable', **kwargs)

###############################################################################
# ################## TRIGONOMETRIC INTEGRALS ################################ #
###############################################################################


class TrigonometricIntegral(Function):
    """Base class for trigonometric integrals."""

    @classmethod
    def eval(cls, z):
        if z == 0:
            return cls._atzero
        if z is oo:
            return cls._atinf()
        if z == -oo:
            return cls._atneginf()

        nz = z.as_coefficient(polar_lift(I))
        if nz is None and cls._trigfunc(0) == 0:
            nz = z.as_coefficient(I)
        if nz is not None:
            return cls._Ifactor(nz, 1)
        nz = z.as_coefficient(polar_lift(-I))
        if nz is not None:
            return cls._Ifactor(nz, -1)

        nz = z.as_coefficient(polar_lift(-1))
        if nz is None and cls._trigfunc(0) == 0:
            nz = z.as_coefficient(-1)
        if nz is not None:
            return cls._minusfactor(nz)

        nz, n = z.extract_branch_factor()
        if n != 0 or nz != z:
            return 2*pi*I*n*cls._trigfunc(0) + cls(nz)

    def fdiff(self, argindex=1):
        from .. import unpolarify
        arg = unpolarify(self.args[0])
        if argindex == 1:
            return self._trigfunc(arg)/arg
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_Ei(self, z):
        return self._eval_rewrite_as_expint(z).rewrite(Ei)

    def _eval_nseries(self, x, n, logx):
        # NOTE this is fairly inefficient
        n += 1
        if self.args[0].subs({x: 0}) != 0:
            return super()._eval_nseries(x, n, logx)
        baseseries = self._trigfunc(x)._eval_nseries(x, n, logx)
        if self._trigfunc(0) != 0:
            baseseries -= 1
        baseseries = baseseries.replace(Pow, lambda t, n: t**n/n)
        if self._trigfunc(0) != 0:
            baseseries += EulerGamma + log(x)
        return baseseries.subs({x: self.args[0]})._eval_nseries(x, n, logx)


class Si(TrigonometricIntegral):
    r"""
    Sine integral.

    This function is defined by

    .. math:: \operatorname{Si}(z) = \int_0^z \frac{\sin{t}}{t} \mathrm{d}t.

    It is an entire function.

    Examples
    ========

    The sine integral is an antiderivative of sin(z)/z:

    >>> Si(z).diff(z)
    sin(z)/z

    It is unbranched:

    >>> Si(z*exp_polar(2*I*pi))
    Si(z)

    Sine integral behaves much like ordinary sine under multiplication by ``I``:

    >>> Si(I*z)
    I*Shi(z)
    >>> Si(-z)
    -Si(z)

    It can also be expressed in terms of exponential integrals, but beware
    that the latter is branched:

    >>> Si(z).rewrite(expint)
    -I*(-expint(1, z*exp_polar(-I*pi/2))/2 +
         expint(1, z*exp_polar(I*pi/2))/2) + pi/2

    See Also
    ========

    Ci: Cosine integral.
    Shi: Hyperbolic sine integral.
    Chi: Hyperbolic cosine integral.
    Ei: Exponential integral.
    expint: Generalized exponential integral.
    E1: Special case of the generalized exponential integral.
    li: Logarithmic integral.
    Li: Offset logarithmic integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_integral

    """

    _trigfunc = sin
    _atzero = Integer(0)

    @classmethod
    def _atinf(cls):
        return pi/2

    @classmethod
    def _atneginf(cls):
        return -pi/2

    @classmethod
    def _minusfactor(cls, z):
        return -Si(z)

    @classmethod
    def _Ifactor(cls, z, sign):
        return I*Shi(z)*sign

    def _eval_rewrite_as_expint(self, z):
        # XXX should we polarify z?
        return pi/2 + (E1(polar_lift(I)*z) - E1(polar_lift(-I)*z))/2/I


class Ci(TrigonometricIntegral):
    r"""
    Cosine integral.

    This function is defined for positive `x` by

    .. math:: \operatorname{Ci}(x) = \gamma + \log{x}
                         + \int_0^x \frac{\cos{t} - 1}{t} \mathrm{d}t
           = -\int_x^\infty \frac{\cos{t}}{t} \mathrm{d}t,

    where `\gamma` is the Euler-Mascheroni constant.

    We have

    .. math:: \operatorname{Ci}(z) =
        -\frac{\operatorname{E}_1\left(e^{i\pi/2} z\right)
               + \operatorname{E}_1\left(e^{-i \pi/2} z\right)}{2}

    which holds for all polar `z` and thus provides an analytic
    continuation to the Riemann surface of the logarithm.

    The formula also holds as stated
    for `z \in \mathbb{C}` with `\Re(z) > 0`.
    By lifting to the principal branch we obtain an analytic function on the
    cut complex plane.

    Examples
    ========

    The cosine integral is a primitive of `\cos(z)/z`:

    >>> Ci(z).diff(z)
    cos(z)/z

    It has a logarithmic branch point at the origin:

    >>> Ci(z*exp_polar(2*I*pi))
    Ci(z) + 2*I*pi

    The cosine integral behaves somewhat like ordinary `\cos` under multiplication by `i`:

    >>> Ci(polar_lift(I)*z)
    Chi(z) + I*pi/2
    >>> Ci(polar_lift(-1)*z)
    Ci(z) + I*pi

    It can also be expressed in terms of exponential integrals:

    >>> Ci(z).rewrite(expint)
    -expint(1, z*exp_polar(-I*pi/2))/2 - expint(1, z*exp_polar(I*pi/2))/2

    See Also
    ========

    Si: Sine integral.
    Shi: Hyperbolic sine integral.
    Chi: Hyperbolic cosine integral.
    Ei: Exponential integral.
    expint: Generalized exponential integral.
    E1: Special case of the generalized exponential integral.
    li: Logarithmic integral.
    Li: Offset logarithmic integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_integral

    """

    _trigfunc = cos
    _atzero = zoo

    @classmethod
    def _atinf(cls):
        return Integer(0)

    @classmethod
    def _atneginf(cls):
        return I*pi

    @classmethod
    def _minusfactor(cls, z):
        return Ci(z) + I*pi

    @classmethod
    def _Ifactor(cls, z, sign):
        return Chi(z) + I*pi/2*sign

    def _eval_rewrite_as_expint(self, z):
        return -(E1(polar_lift(I)*z) + E1(polar_lift(-I)*z))/2


class Shi(TrigonometricIntegral):
    r"""
    Sinh integral.

    This function is defined by

    .. math:: \operatorname{Shi}(z) = \int_0^z \frac{\sinh{t}}{t} \mathrm{d}t.

    It is an entire function.

    Examples
    ========

    The Sinh integral is a primitive of `\sinh(z)/z`:

    >>> Shi(z).diff(z)
    sinh(z)/z

    It is unbranched:

    >>> Shi(z*exp_polar(2*I*pi))
    Shi(z)

    The `\sinh` integral behaves much like ordinary `\sinh` under multiplication by `i`:

    >>> Shi(I*z)
    I*Si(z)
    >>> Shi(-z)
    -Shi(z)

    It can also be expressed in terms of exponential integrals, but beware
    that the latter is branched:

    >>> Shi(z).rewrite(expint)
    expint(1, z)/2 - expint(1, z*exp_polar(I*pi))/2 - I*pi/2

    See Also
    ========

    Si: Sine integral.
    Ci: Cosine integral.
    Chi: Hyperbolic cosine integral.
    Ei: Exponential integral.
    expint: Generalized exponential integral.
    E1: Special case of the generalized exponential integral.
    li: Logarithmic integral.
    Li: Offset logarithmic integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_integral

    """

    _trigfunc = sinh
    _atzero = Integer(0)

    @classmethod
    def _atinf(cls):
        return oo

    @classmethod
    def _atneginf(cls):
        return -oo

    @classmethod
    def _minusfactor(cls, z):
        return -Shi(z)

    @classmethod
    def _Ifactor(cls, z, sign):
        return I*Si(z)*sign

    def _eval_rewrite_as_expint(self, z):
        from .. import exp_polar

        # XXX should we polarify z?
        return (E1(z) - E1(exp_polar(I*pi)*z))/2 - I*pi/2


class Chi(TrigonometricIntegral):
    r"""
    Cosh integral.

    This function is defined for positive `x` by

    .. math:: \operatorname{Chi}(x) = \gamma + \log{x}
                         + \int_0^x \frac{\cosh{t} - 1}{t} \mathrm{d}t,

    where `\gamma` is the Euler-Mascheroni constant.

    We have

    .. math:: \operatorname{Chi}(z) = \operatorname{Ci}\left(e^{i \pi/2}z\right)
                         - i\frac{\pi}{2},

    which holds for all polar `z` and thus provides an analytic
    continuation to the Riemann surface of the logarithm.
    By lifting to the principal branch we obtain an analytic function on the
    cut complex plane.

    Examples
    ========

    The `\cosh` integral is a primitive of `\cosh(z)/z`:

    >>> Chi(z).diff(z)
    cosh(z)/z

    It has a logarithmic branch point at the origin:

    >>> Chi(z*exp_polar(2*I*pi))
    Chi(z) + 2*I*pi

    The `\cosh` integral behaves somewhat like ordinary `\cosh` under multiplication by `i`:

    >>> Chi(polar_lift(I)*z)
    Ci(z) + I*pi/2
    >>> Chi(polar_lift(-1)*z)
    Chi(z) + I*pi

    It can also be expressed in terms of exponential integrals:

    >>> Chi(z).rewrite(expint)
    -expint(1, z)/2 - expint(1, z*exp_polar(I*pi))/2 - I*pi/2

    See Also
    ========

    Si: Sine integral.
    Ci: Cosine integral.
    Shi: Hyperbolic sine integral.
    Ei: Exponential integral.
    expint: Generalized exponential integral.
    E1: Special case of the generalized exponential integral.
    li: Logarithmic integral.
    Li: Offset logarithmic integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_integral

    """

    _trigfunc = cosh
    _atzero = zoo

    @classmethod
    def _atinf(cls):
        return oo

    @classmethod
    def _atneginf(cls):
        return oo

    @classmethod
    def _minusfactor(cls, z):
        return Chi(z) + I*pi

    @classmethod
    def _Ifactor(cls, z, sign):
        return Ci(z) + I*pi/2*sign

    def _eval_rewrite_as_expint(self, z):
        from .. import exp_polar
        return -I*pi/2 - (E1(z) + E1(exp_polar(I*pi)*z))/2

    def _latex(self, printer, exp=None):
        if exp:
            return (r'\operatorname{Chi}^{%s}{\left (%s \right )}'  # noqa: SFS101
                    % (printer._print(exp), printer._print(self.args[0])))
        return (r'\operatorname{Chi}{\left (%s \right )}'  # noqa: SFS101
                % printer._print(self.args[0]))

    @staticmethod
    def _latex_no_arg(printer):
        return r'\operatorname{Chi}'


###############################################################################
# ################## FRESNEL INTEGRALS ###################################### #
###############################################################################

class FresnelIntegral(Function):
    """Base class for the Fresnel integrals."""

    unbranched = True

    @classmethod
    def eval(cls, z):
        # Value at zero
        if z == 0:
            return Integer(0)

        # Try to pull out factors of -1 and I
        prefact = Integer(1)
        newarg = z
        changed = False

        nz = newarg.as_coefficient(-1)
        if nz is not None:
            prefact = -prefact
            newarg = nz
            changed = True

        nz = newarg.as_coefficient(I)
        if nz is not None:
            prefact = cls._sign*I*prefact
            newarg = nz
            changed = True

        if changed:
            return prefact*cls(newarg)

        # Values at positive infinities signs
        # if any were extracted automatically
        if z is oo:
            return Rational(1, 2)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return self._trigfunc(pi*self.args[0]**2/2)
        raise ArgumentIndexError(self, argindex)

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.expand(deep, **hints), Integer(0)
            return self, Integer(0)
        if deep:
            x, y = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            x, y = self.args[0].as_real_imag()

        # Fresnel S
        # http://functions.wolfram.com/06.32.19.0003.01
        # http://functions.wolfram.com/06.32.19.0006.01
        # Fresnel C
        # http://functions.wolfram.com/06.33.19.0003.01
        # http://functions.wolfram.com/06.33.19.0006.01
        if x.is_zero:
            re, im = self.func(I*y).rewrite(erf).as_real_imag()
        else:
            sq = -y**2/x**2
            re = (self.func(x + x*sqrt(sq)) + self.func(x - x*sqrt(sq)))/2
            im = x/(2*y)*sqrt(sq)*(self.func(x - x*sqrt(sq)) -
                                   self.func(x + x*sqrt(sq)))
        return re, im


class fresnels(FresnelIntegral):
    r"""
    Fresnel integral S.

    This function is defined by

    .. math:: \operatorname{S}(z) = \int_0^z \sin{\frac{\pi}{2} t^2} \mathrm{d}t.

    It is an entire function.

    Examples
    ========

    Several special values are known:

    >>> fresnels(0)
    0
    >>> fresnels(oo)
    1/2
    >>> fresnels(-oo)
    -1/2
    >>> fresnels(I*oo)
    -I/2
    >>> fresnels(-I*oo)
    I/2

    In general one can pull out factors of -1 and `i` from the argument:

    >>> fresnels(-z)
    -fresnels(z)
    >>> fresnels(I*z)
    -I*fresnels(z)

    The Fresnel S integral obeys the mirror symmetry
    `\overline{S(z)} = S(\bar{z})`:

    >>> conjugate(fresnels(z))
    fresnels(conjugate(z))

    Differentiation with respect to `z` is supported:

    >>> diff(fresnels(z), z)
    sin(pi*z**2/2)

    Defining the Fresnel functions via an integral

    >>> integrate(sin(pi*z**2/2), z)
    3*fresnels(z)*gamma(3/4)/(4*gamma(7/4))
    >>> expand_func(integrate(sin(pi*z**2/2), z))
    fresnels(z)

    We can numerically evaluate the Fresnel integral to arbitrary precision
    on the whole complex plane:

    >>> fresnels(2).evalf(30)
    0.343415678363698242195300815958

    >>> fresnels(-2*I).evalf(30)
    0.343415678363698242195300815958*I

    See Also
    ========

    fresnelc: Fresnel cosine integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Fresnel_integral
    * https://dlmf.nist.gov/7
    * https://mathworld.wolfram.com/FresnelIntegrals.html
    * http://functions.wolfram.com/GammaBetaErf/FresnelS
    * The converging factors for the fresnel integrals
      by John W. Wrench Jr. and Vicki Alley

    """

    _trigfunc = sin
    _sign = -Integer(1)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0:
            return Integer(0)
        x = sympify(x)
        if len(previous_terms) >= 1:
            p = previous_terms[-1]
            return (-pi**2*x**4*(4*n - 1)/(8*n*(2*n + 1)*(4*n + 3))) * p
        return x**3 * (-x**4)**n * (Integer(2)**(-2*n - 1)*pi**(2*n + 1)) / ((4*n + 3)*factorial(2*n + 1))

    def _eval_rewrite_as_erf(self, z):
        return (1 + I)/4 * (erf((1 + I)/2*sqrt(pi)*z) - I*erf((1 - I)/2*sqrt(pi)*z))

    def _eval_rewrite_as_hyper(self, z):
        return pi*z**3/6 * hyper([Rational(3, 4)], [Rational(3, 2), Rational(7, 4)], -pi**2*z**4/16)

    def _eval_rewrite_as_meijerg(self, z):
        return (pi*z**Rational(9, 4) / (sqrt(2)*(z**2)**Rational(3, 4)*(-z)**Rational(3, 4))
                * meijerg([], [1], [Rational(3, 4)], [Rational(1, 4), 0], -pi**2*z**4/16))

    def _eval_aseries(self, n, args0, x, logx):
        from ...calculus import Order
        point = args0[0]

        # Expansion at oo
        if point is oo:
            z = self.args[0]

            # expansion of S(x) = S1(x*sqrt(pi/2)), see reference[5] page 1-8
            p = [(-1)**k * factorial(4*k + 1) /
                 (2**(2*k + 2) * z**(4*k + 3) * 2**(2*k)*factorial(2*k))
                 for k in range(n)]
            q = [1/(2*z)] + [(-1)**k * factorial(4*k - 1) /
                             (2**(2*k + 1) * z**(4*k + 1) * 2**(2*k - 1)*factorial(2*k - 1))
                             for k in range(1, n)]

            p = [-sqrt(2/pi)*t for t in p] + [Order(1/z**n, x)]
            q = [-sqrt(2/pi)*t for t in q] + [Order(1/z**n, x)]

            return Rational(1, 2) + (sin(z**2)*Add(*p) + cos(z**2)*Add(*q)).subs({x: sqrt(2/pi)*x})

        # All other points are not handled
        return super()._eval_aseries(n, args0, x, logx)


class fresnelc(FresnelIntegral):
    r"""
    Fresnel integral C.

    This function is defined by

    .. math:: \operatorname{C}(z) = \int_0^z \cos{\frac{\pi}{2} t^2} \mathrm{d}t.

    It is an entire function.

    Examples
    ========

    Several special values are known:

    >>> fresnelc(0)
    0
    >>> fresnelc(oo)
    1/2
    >>> fresnelc(-oo)
    -1/2
    >>> fresnelc(I*oo)
    I/2
    >>> fresnelc(-I*oo)
    -I/2

    In general one can pull out factors of -1 and `i` from the argument:

    >>> fresnelc(-z)
    -fresnelc(z)
    >>> fresnelc(I*z)
    I*fresnelc(z)

    The Fresnel C integral obeys the mirror symmetry
    `\overline{C(z)} = C(\bar{z})`:

    >>> conjugate(fresnelc(z))
    fresnelc(conjugate(z))

    Differentiation with respect to `z` is supported:

    >>> diff(fresnelc(z), z)
    cos(pi*z**2/2)

    Defining the Fresnel functions via an integral

    >>> integrate(cos(pi*z**2/2), z)
    fresnelc(z)*gamma(1/4)/(4*gamma(5/4))
    >>> expand_func(integrate(cos(pi*z**2/2), z))
    fresnelc(z)

    We can numerically evaluate the Fresnel integral to arbitrary precision
    on the whole complex plane:

    >>> fresnelc(2).evalf(30)
    0.488253406075340754500223503357

    >>> fresnelc(-2*I).evalf(30)
    -0.488253406075340754500223503357*I

    See Also
    ========

    fresnels: Fresnel sine integral.

    References
    ==========

    * https://en.wikipedia.org/wiki/Fresnel_integral
    * https://dlmf.nist.gov/7
    * https://mathworld.wolfram.com/FresnelIntegrals.html
    * http://functions.wolfram.com/GammaBetaErf/FresnelC
    * The converging factors for the fresnel integrals
      by John W. Wrench Jr. and Vicki Alley

    """

    _trigfunc = cos
    _sign = Integer(1)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0:
            return Integer(0)
        x = sympify(x)
        if len(previous_terms) >= 1:
            p = previous_terms[-1]
            return (-pi**2*x**4*(4*n - 3)/(8*n*(2*n - 1)*(4*n + 1))) * p
        return x * (-x**4)**n * (Integer(2)**(-2*n)*pi**(2*n)) / ((4*n + 1)*factorial(2*n))

    def _eval_rewrite_as_erf(self, z):
        return (1 - I)/4 * (erf((1 + I)/2*sqrt(pi)*z) + I*erf((1 - I)/2*sqrt(pi)*z))

    def _eval_rewrite_as_hyper(self, z):
        return z * hyper([Rational(1, 4)], [Rational(1, 2), Rational(5, 4)], -pi**2*z**4/16)

    def _eval_rewrite_as_meijerg(self, z):
        return (pi*z**Rational(3, 4) / (sqrt(2)*root(z**2, 4)*root(-z, 4))
                * meijerg([], [1], [Rational(1, 4)], [Rational(3, 4), 0], -pi**2*z**4/16))

    def _eval_aseries(self, n, args0, x, logx):
        from ...calculus import Order
        point = args0[0]

        # Expansion at oo
        if point is oo:
            z = self.args[0]

            # expansion of C(x) = C1(x*sqrt(pi/2)), see reference[5] page 1-8
            p = [(-1)**k * factorial(4*k + 1) /
                 (2**(2*k + 2) * z**(4*k + 3) * 2**(2*k)*factorial(2*k))
                 for k in range(n)]
            q = [1/(2*z)] + [(-1)**k * factorial(4*k - 1) /
                             (2**(2*k + 1) * z**(4*k + 1) * 2**(2*k - 1)*factorial(2*k - 1))
                             for k in range(1, n)]

            p = [-sqrt(2/pi)*t for t in p] + [Order(1/z**n, x)]
            q = [+sqrt(2/pi)*t for t in q] + [Order(1/z**n, x)]

            return Rational(1, 2) + (cos(z**2)*Add(*p) + sin(z**2)*Add(*q)).subs({x: sqrt(2/pi)*x})

        # All other points are not handled
        return super()._eval_aseries(n, args0, x, logx)


###############################################################################
# ################## HELPER FUNCTIONS ####################################### #
###############################################################################


class _erfs(Function):
    r"""
    Helper function to make the `\mathrm{erf}(z)` function
    tractable for the Gruntz algorithm.

    """

    @classmethod
    def eval(cls, z):
        if z.is_zero:
            return Integer(1)

    def _eval_aseries(self, n, args0, x, logx):
        from ...calculus import Order
        point = args0[0]

        # Expansion at oo
        if point is oo:
            z = self.args[0]
            l = [1/sqrt(pi)*factorial(2*k)*(-Integer(4))**(-k) /
                 factorial(k)*(1/z)**(2*k + 1) for k in range(n)]
            o = Order(1/z**(2*n + 1), x)
            # It is very inefficient to first add the order and then do the nseries
            return (Add(*l))._eval_nseries(x, n, logx) + o

        # All other points are not handled
        return super()._eval_aseries(n, args0, x, logx)

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = self.args[0]
            return -2/sqrt(pi) + 2*z*_erfs(z)
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_intractable(self, z):
        return (1 - erf(z))*exp(z**2)

    def _eval_evalf(self, prec):
        return self.rewrite('intractable').evalf(prec)


class _eis(Function):
    r"""
    Helper function to make the `\mathrm{Ei}(z)` and `\mathrm{li}(z)` functions
    tractable for the Gruntz algorithm.

    """

    def _eval_aseries(self, n, args0, x, logx):
        from ...calculus import Order
        if args0[0] != oo:
            return super()._eval_aseries(n, args0, x, logx)

        z = self.args[0]
        l = [factorial(k) * (1/z)**(k + 1) for k in range(n)]
        o = Order(1/z**(n + 1), x)
        # It is very inefficient to first add the order and then do the nseries
        return (Add(*l))._eval_nseries(x, n, logx) + o

    def fdiff(self, argindex=1):
        if argindex == 1:
            z = self.args[0]
            return 1/z - _eis(z)
        raise ArgumentIndexError(self, argindex)

    def _eval_rewrite_as_intractable(self, z):
        return exp(-z)*Ei(z)

    def _eval_nseries(self, x, n, logx):
        x0 = self.args[0].limit(x, 0)
        if x0 == 0:
            f = self._eval_rewrite_as_intractable(*self.args)
            return f._eval_nseries(x, n, logx)
        return super()._eval_nseries(x, n, logx)

    def _eval_evalf(self, prec):
        return self.rewrite('intractable').evalf(prec)
