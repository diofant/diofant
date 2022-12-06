from ...core import (Add, Derivative, Dummy, Eq, Expr, Function, I, Integer,
                     Mul, Rational, Symbol, Tuple, factor_terms, nan, oo, pi,
                     zoo)
from ...core.function import AppliedUndef, ArgumentIndexError
from ...core.sympify import sympify
from ...logic.boolalg import BooleanAtom
from .exponential import exp, exp_polar, log
from .miscellaneous import sqrt
from .piecewise import ExprCondPair, Piecewise
from .trigonometric import atan2


###############################################################################
# ####################### REAL and IMAGINARY PARTS ########################## #
###############################################################################


class re(Function):
    """Returns real part of expression.

    This function performs only elementary analysis and so it will fail to
    decompose properly more complicated expressions. If completely simplified
    result is needed then use Basic.as_real_imag() or perform complex
    expansion on instance of this function.

    Examples
    ========

    >>> re(2*E)
    2*E
    >>> re(2*I + 17)
    17
    >>> re(2*I)
    0
    >>> re(im(x) + x*I + 2)
    2

    See Also
    ========

    diofant.functions.elementary.complexes.im

    """

    is_extended_real = True
    unbranched = True  # implicitely works on the projection to C

    @classmethod
    def eval(cls, arg):
        if arg is zoo:
            return nan
        elif arg.is_extended_real:
            return arg
        elif arg.is_imaginary or (I*arg).is_extended_real:
            return Integer(0)
        elif arg.is_Function and isinstance(arg, conjugate):
            return re(arg.args[0])
        else:

            included, reverted, excluded = [], [], []
            args = Add.make_args(arg)
            for term in args:
                coeff = term.as_coefficient(I)

                if coeff is not None:
                    if not coeff.is_extended_real:
                        reverted.append(coeff)
                elif not term.has(I) and term.is_extended_real:
                    excluded.append(term)
                else:
                    # Try to do some advanced expansion.  If
                    # impossible, don't try to do re(arg) again
                    # (because this is what we are trying to do now).
                    real_imag = term.as_real_imag(ignore=arg)
                    if real_imag:
                        excluded.append(real_imag[0])
                    else:
                        included.append(term)

            if len(args) != len(included):
                a, b, c = (Add(*xs) for xs in [included, reverted, excluded])

                return cls(a) - im(b) + c

    def as_real_imag(self, deep=True, **hints):
        """Returns the real number with a zero imaginary part."""
        return self, Integer(0)

    def _eval_derivative(self, s):
        if s.is_extended_real or self.args[0].is_extended_real:
            return re(Derivative(self.args[0], s, evaluate=True))
        elif s.is_imaginary or self.args[0].is_imaginary:
            return -I*im(Derivative(self.args[0], s, evaluate=True))

    def _eval_rewrite_as_im(self, arg):
        return self.args[0] - I*im(self.args[0])

    def _eval_is_algebraic(self):
        return self.args[0].is_algebraic

    def _eval_is_real(self):
        if self.args[0].is_complex:
            return True


class im(Function):
    """Returns imaginary part of expression.

    This function performs only elementary analysis and so it will fail to
    decompose properly more complicated expressions. If completely simplified
    result is needed then use Basic.as_real_imag() or perform complex expansion
    on instance of this function.

    Examples
    ========

    >>> im(2*E)
    0
    >>> re(2*I + 17)
    17
    >>> im(x*I)
    re(x)
    >>> im(re(x) + y)
    im(y)

    See Also
    ========

    diofant.functions.elementary.complexes.re

    """

    is_extended_real = True
    unbranched = True  # implicitely works on the projection to C

    @classmethod
    def eval(cls, arg):
        if arg is zoo:
            return nan
        elif arg.is_extended_real:
            return Integer(0)
        elif arg.is_imaginary or (I*arg).is_extended_real:
            return -I * arg
        elif arg.is_Function and isinstance(arg, conjugate):
            return -im(arg.args[0])
        else:
            included, reverted, excluded = [], [], []
            args = Add.make_args(arg)
            for term in args:
                coeff = term.as_coefficient(I)

                if coeff is not None:
                    if not coeff.is_extended_real:
                        reverted.append(coeff)
                    else:
                        excluded.append(coeff)
                elif term.has(I) or not term.is_extended_real:
                    # Try to do some advanced expansion.  If
                    # impossible, don't try to do im(arg) again
                    # (because this is what we are trying to do now).
                    real_imag = term.as_real_imag(ignore=arg)
                    if real_imag:
                        excluded.append(real_imag[1])
                    else:
                        included.append(term)

            if len(args) != len(included):
                a, b, c = (Add(*xs) for xs in [included, reverted, excluded])

                return cls(a) + re(b) + c

    def as_real_imag(self, deep=True, **hints):
        """
        Return the imaginary part with a zero real part.

        Examples
        ========

        >>> im(2 + 3*I).as_real_imag()
        (3, 0)

        """
        return self, Integer(0)

    def _eval_derivative(self, s):
        if s.is_extended_real or self.args[0].is_extended_real:
            return im(Derivative(self.args[0], s, evaluate=True))
        elif s.is_imaginary or self.args[0].is_imaginary:
            return -I*re(Derivative(self.args[0], s, evaluate=True))

    def _eval_rewrite_as_re(self, arg):
        return -I*(self.args[0] - re(self.args[0]))

    def _eval_is_algebraic(self):
        return self.args[0].is_algebraic

    def _eval_is_real(self):
        if self.args[0].is_complex:
            return True


###############################################################################
# ############# SIGN, ABSOLUTE VALUE, ARGUMENT and CONJUGATION ############## #
###############################################################################

class sign(Function):
    """
    Returns the complex sign of an expression.

    For nonzero complex number z is an equivalent of z/abs(z).
    Else returns zero.

    Examples
    ========

    >>> sign(-1)
    -1
    >>> sign(0)
    0
    >>> sign(-3*I)
    -I
    >>> sign(1 + I)
    sign(1 + I)
    >>> _.evalf()
    0.707106781186548 + 0.707106781186548*I

    See Also
    ========

    Abs
    conjugate

    """

    is_complex = True

    def doit(self, **hints):
        if self.args[0].is_nonzero:
            return self.args[0] / Abs(self.args[0])
        return self

    @classmethod
    def eval(cls, arg):
        # handle what we can
        if arg.is_Mul:
            c, args = arg.as_coeff_mul()
            unk = []
            s = sign(c)
            for a in args:
                if a.is_negative:
                    s = -s
                elif a.is_positive:
                    pass
                elif a.is_imaginary and im(a).is_comparable:
                    s *= sign(a)
                else:
                    unk.append(a)
            if c == 1 and len(unk) == len(args):
                return
            return s * cls(arg._new_rawargs(*unk))
        if arg.is_zero:  # it may be an Expr that is zero
            return Integer(0)
        if arg.is_positive:
            return Integer(1)
        if arg.is_negative:
            return Integer(-1)
        if arg.is_Function:
            if isinstance(arg, sign):
                return arg
        if arg.is_imaginary:
            if arg.is_Pow and arg.exp == Rational(1, 2):
                # we catch this because non-trivial sqrt args are not expanded
                # e.g. sqrt(1-sqrt(2)) --x-->  to I*sqrt(sqrt(2) - 1)
                return I
            arg2 = -I * arg
            if arg2.is_positive:
                return I

    def _eval_Abs(self):
        if self.args[0].is_nonzero:
            return Integer(1)

    def _eval_conjugate(self):
        return sign(conjugate(self.args[0]))

    def _eval_derivative(self, s):
        if self.args[0].is_extended_real:
            from ..special.delta_functions import DiracDelta
            return 2 * Derivative(self.args[0], s, evaluate=True) \
                * DiracDelta(self.args[0])
        elif self.args[0].is_imaginary:
            from ..special.delta_functions import DiracDelta
            return 2 * Derivative(self.args[0], s, evaluate=True) \
                * DiracDelta(-I * self.args[0])

    def _eval_is_nonnegative(self):
        if self.args[0].is_nonnegative:
            return True

    def _eval_is_nonpositive(self):
        if self.args[0].is_nonpositive:
            return True

    def _eval_is_imaginary(self):
        return self.args[0].is_imaginary

    def _eval_is_integer(self):
        return self.args[0].is_extended_real

    def _eval_is_zero(self):
        return self.args[0].is_zero

    def _eval_power(self, other):
        if (
            self.args[0].is_extended_real and
            self.args[0].is_nonzero and
            other.is_integer and
            other.is_even
        ):
            return Integer(1)

    def _eval_rewrite_as_Piecewise(self, arg):
        if arg.is_extended_real:
            return Piecewise((1, arg > 0), (-1, arg < 0), (0, True))

    def _eval_rewrite_as_Heaviside(self, arg):
        from .. import Heaviside
        if arg.is_extended_real:
            return Heaviside(arg)*2-1

    def _eval_rewrite_as_exp(self, x):
        return exp(I*arg(x))

    def _eval_simplify(self, ratio, measure):
        return self.func(self.args[0].factor())

    def _eval_nseries(self, x, n, logx):
        direction = self.args[0].as_leading_term(x).as_coeff_exponent(x)[0]
        if direction.is_extended_real:
            return self.func(direction)
        else:
            return super()._eval_nseries(x, n, logx)


class Abs(Function):
    """Return the absolute value of the argument.

    This is an extension of the built-in function abs() to accept symbolic
    values.  If you pass a Diofant expression to the built-in abs(), it will
    pass it automatically to Abs().

    Examples
    ========

    >>> Abs(-1)
    1
    >>> x = Symbol('x', real=True)
    >>> abs(-x)  # The Python built-in
    Abs(x)
    >>> abs(x**2)
    x**2

    Note that the Python built-in will return either an Expr or int depending on
    the argument:

    >>> type(abs(-1))
    <... 'int'>
    >>> type(abs(Integer(-1)))
    <class 'diofant.core.numbers.One'>

    Abs will always return a diofant object.

    See Also
    ========

    diofant.functions.elementary.complexes.sign
    diofant.functions.elementary.complexes.conjugate

    """

    is_extended_real = True
    is_negative = False
    unbranched = True

    def fdiff(self, argindex=1):
        """
        Get the first derivative of the argument to Abs().

        Examples
        ========

        >>> abs(-x).fdiff()
        sign(x)

        """
        if argindex == 1:
            return sign(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from ...simplify import signsimp
        if hasattr(arg, '_eval_Abs'):
            obj = arg._eval_Abs()
            if obj is not None:
                return obj
        if not isinstance(arg, Expr):
            raise TypeError(f'Bad argument type for Abs(): {type(arg)}')
        # handle what we can
        arg = signsimp(arg, evaluate=False)
        if arg.is_Mul:
            known = []
            unk = []
            for t in arg.args:
                tnew = cls(t)
                if isinstance(tnew, cls):
                    unk.append(tnew.args[0])
                else:
                    known.append(tnew)
            known = Mul(*known)
            unk = cls(Mul(*unk), evaluate=False) if unk else Integer(1)
            return known*unk
        if arg.is_Pow:
            base, exponent = arg.as_base_exp()
            if base.is_extended_real:
                if exponent.is_integer:
                    if exponent.is_even:
                        return arg
                    if base == -1:
                        return Integer(1)
                    if isinstance(base, cls) and exponent == -1:
                        return arg
                    return Abs(base)**exponent
                if base.is_nonnegative:
                    return base**re(exponent)
                if base.is_negative:
                    return (-base)**re(exponent)*exp(-pi*im(exponent))
        if arg.is_zero:  # it may be an Expr that is zero
            return Integer(0)
        if arg.is_nonnegative:
            return arg
        if arg.is_nonpositive:
            return -arg
        if arg.is_imaginary:
            arg2 = -I * arg
            if arg2.is_nonnegative:
                return arg2
        if arg.is_Add:
            if any(a.is_infinite for a in arg.as_real_imag()):
                return oo
            if arg.is_extended_real is not True and arg.is_imaginary is None:
                if all(a.is_extended_real or a.is_imaginary or (I*a).is_extended_real for a in arg.args):
                    from ...core import expand_mul
                    return sqrt(expand_mul(arg*arg.conjugate()))
        if arg is zoo:
            return oo
        if arg.is_extended_real is not True and arg.is_imaginary is False:
            from ...core import expand_mul
            return sqrt(expand_mul(arg*arg.conjugate()))

    def _eval_is_integer(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_integer

    def _eval_is_nonzero(self):
        return self.args[0].is_nonzero

    def _eval_is_finite(self):
        if self.args[0].is_complex:
            return True

    def _eval_is_positive(self):
        return self.is_nonzero

    def _eval_is_rational(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_rational

    def _eval_is_even(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_even

    def _eval_is_odd(self):
        if self.args[0].is_extended_real:
            return self.args[0].is_odd

    def _eval_is_algebraic(self):
        return self.args[0].is_algebraic

    def _eval_power(self, other):
        if self.args[0].is_extended_real and other.is_integer:
            if other.is_even:
                return self.args[0]**other
            elif other != -1 and other.is_Integer:
                return self.args[0]**(other - 1)*self

    def _eval_nseries(self, x, n, logx):
        direction = self.args[0].as_leading_term(x).as_coeff_exponent(x)[0]
        s = self.args[0]._eval_nseries(x, n=n, logx=logx)
        when, lim = Eq(direction, 0), direction.limit(x, 0)
        if lim.equals(0) is False:
            return s/sign(lim)
        else:
            return Piecewise((lim, when), (s/sign(direction), True))

    def _eval_derivative(self, s):
        if self.args[0].is_extended_real or self.args[0].is_imaginary:
            return Derivative(self.args[0], s, evaluate=True) \
                * sign(conjugate(self.args[0]))
        return (re(self.args[0]) * Derivative(re(self.args[0]), s,
                                              evaluate=True) + im(self.args[0]) * Derivative(im(self.args[0]),
                                                                                             s, evaluate=True)) / Abs(self.args[0])

    def _eval_rewrite_as_Heaviside(self, arg):
        # Note this only holds for real arg (since Heaviside is not defined
        # for complex arguments).
        from .. import Heaviside
        if arg.is_extended_real:
            return arg*(Heaviside(arg) - Heaviside(-arg))

    def _eval_rewrite_as_Piecewise(self, arg):
        if arg.is_extended_real:
            return Piecewise((arg, arg >= 0), (-arg, True))

    def _eval_rewrite_as_sign(self, arg):
        return arg/sign(arg)

    def _eval_rewrite_as_tractable(self, arg, wrt=None, **kwargs):
        if wrt is not None and (s := sign(arg.limit(wrt, oo))) in (1, -1):
            return s*arg


class arg(Function):
    """Returns the argument (in radians) of a complex number.

    For a real number, the argument is always 0.

    Examples
    ========

    >>> arg(2.0)
    0
    >>> arg(I)
    pi/2
    >>> arg(sqrt(2) + I*sqrt(2))
    pi/4

    """

    is_real = True

    @classmethod
    def eval(cls, arg):
        if not arg.is_Atom:
            c, arg_ = factor_terms(arg).as_coeff_Mul()
            if arg_.is_Mul:
                arg_ = Mul(*[a if (sign(a) not in (-1, 1)) else
                             sign(a) for a in arg_.args])
            arg_ = sign(c)*arg_
        else:
            arg_ = arg
        if arg_.is_zero:
            return Integer(0)
        x, y = re(arg_), im(arg_)
        rv = atan2(y, x)
        if rv.is_number and not rv.atoms(AppliedUndef):
            return rv
        if arg_ != arg:
            return cls(arg_, evaluate=False)

    def _eval_derivative(self, s):
        x, y = re(self.args[0]), im(self.args[0])
        return (x * Derivative(y, s, evaluate=True) - y *
                Derivative(x, s, evaluate=True)) / (x**2 + y**2)

    def _eval_rewrite_as_atan2(self, arg):
        x, y = re(self.args[0]), im(self.args[0])
        return atan2(y, x)


class conjugate(Function):
    """Returns the complex conjugate of an argument.

    In mathematics, the complex conjugate of a complex number
    is given by changing the sign of the imaginary part.

    Thus, the conjugate of the complex number
    `a + i b` (where a and b are real numbers) is `a - i b`

    Examples
    ========

    >>> conjugate(2)
    2
    >>> conjugate(I)
    -I

    See Also
    ========

    diofant.functions.elementary.complexes.sign
    diofant.functions.elementary.complexes.Abs

    References
    ==========

    * https://en.wikipedia.org/wiki/Complex_conjugation

    """

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_conjugate()
        if obj is not None:
            return obj

    def _eval_Abs(self):
        return Abs(self.args[0], evaluate=True)

    def _eval_adjoint(self):
        return transpose(self.args[0])

    def _eval_conjugate(self):
        return self.args[0]

    def _eval_derivative(self, s):
        if s.is_extended_real:
            return conjugate(Derivative(self.args[0], s, evaluate=True))
        elif s.is_imaginary:
            return -conjugate(Derivative(self.args[0], s, evaluate=True))

    def _eval_transpose(self):
        return adjoint(self.args[0])

    def _eval_is_algebraic(self):
        return self.args[0].is_algebraic


class transpose(Function):
    """Linear map transposition."""

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_transpose()
        if obj is not None:
            return obj

    def _eval_adjoint(self):
        return conjugate(self.args[0])

    def _eval_conjugate(self):
        return adjoint(self.args[0])

    def _eval_transpose(self):
        return self.args[0]


class adjoint(Function):
    """Conjugate transpose or Hermite conjugation."""

    @classmethod
    def eval(cls, arg):
        obj = arg._eval_adjoint()
        if obj is not None:
            return obj
        obj = arg._eval_transpose()
        if obj is not None:
            return conjugate(obj)

    def _eval_adjoint(self):
        return self.args[0]

    def _eval_conjugate(self):
        return transpose(self.args[0])

    def _eval_transpose(self):
        return conjugate(self.args[0])

###############################################################################
# ############# HANDLING OF POLAR NUMBERS ################################### #
###############################################################################


class polar_lift(Function):
    """
    Lift argument to the Riemann surface of the logarithm, using the
    standard branch.

    >>> p = Symbol('p', polar=True)
    >>> polar_lift(4)
    4*exp_polar(0)
    >>> polar_lift(-4)
    4*exp_polar(I*pi)
    >>> polar_lift(-I)
    exp_polar(-I*pi/2)
    >>> polar_lift(I + 2)
    polar_lift(2 + I)

    >>> polar_lift(4*x)
    4*polar_lift(x)
    >>> polar_lift(4*p)
    4*p

    See Also
    ========

    diofant.functions.elementary.exponential.exp_polar
    diofant.functions.elementary.complexes.periodic_argument

    """

    is_polar = True
    is_comparable = False  # Cannot be evalf'd.

    @classmethod
    def eval(cls, arg):
        from .. import arg as argument
        from .exponential import exp_polar
        if arg.is_number and (arg.is_finite or arg.is_extended_real):
            ar = argument(arg)
            # In general we want to affirm that something is known,
            # e.g. `not ar.has(argument) and not ar.has(atan)`
            # but for now we will just be more restrictive and
            # see that it has evaluated to one of the known values.
            if ar in (0, pi/2, -pi/2, pi):
                return exp_polar(I*ar)*abs(arg)

        if arg.is_Mul:
            args = arg.args
        else:
            args = [arg]
        included = []
        excluded = []
        positive = []
        for arg in args:
            if arg.is_polar:
                included += [arg]
            elif arg.is_positive:
                positive += [arg]
            else:
                excluded += [arg]
        if len(excluded) < len(args):
            if excluded:
                return Mul(*(included + positive))*polar_lift(Mul(*excluded))
            elif included:
                return Mul(*(included + positive))
            else:
                return Mul(*positive)*exp_polar(0)

    def _eval_evalf(self, prec):
        """Careful! any evalf of polar numbers is flaky."""
        return self.args[0]._eval_evalf(prec)

    def _eval_Abs(self):
        return Abs(self.args[0], evaluate=True)


class periodic_argument(Function):
    """
    Represent the argument on a quotient of the Riemann surface of the
    logarithm. That is, given a period P, always return a value in
    (-P/2, P/2], by using exp(P*I) == 1.

    >>> unbranched_argument(exp(5*I*pi))
    pi
    >>> unbranched_argument(exp_polar(5*I*pi))
    5*pi
    >>> periodic_argument(exp_polar(5*I*pi), 2*pi)
    pi
    >>> periodic_argument(exp_polar(5*I*pi), 3*pi)
    -pi
    >>> periodic_argument(exp_polar(5*I*pi), pi)
    0

    See Also
    ========

    diofant.functions.elementary.exponential.exp_polar
    diofant.functions.elementary.complexes.polar_lift : Lift argument to the Riemann surface of the logarithm
    diofant.functions.elementary.complexes.principal_branch

    """

    @classmethod
    def _getunbranched(cls, ar):
        if ar.is_Mul:
            args = ar.args
        else:
            args = [ar]
        unbranched = 0
        for a in args:
            if not a.is_polar:
                unbranched += arg(a)
            elif isinstance(a, exp_polar):
                unbranched += a.exp.as_real_imag()[1]
            elif a.is_Pow:
                re, im = a.exp.as_real_imag()
                unbranched += re*unbranched_argument(
                    a.base) + im*log(abs(a.base))
            elif isinstance(a, polar_lift):
                unbranched += arg(a.args[0])
            else:
                return
        return unbranched

    @classmethod
    def eval(cls, ar, period):
        # Our strategy is to evaluate the argument on the Riemann surface of the
        # logarithm, and then reduce.
        # NOTE evidently this means it is a rather bad idea to use this with
        # period != 2*pi and non-polar numbers.
        from .integers import ceiling
        from .trigonometric import atan, atan2
        if not period.is_positive:
            return
        if period == oo and isinstance(ar, principal_branch):
            return periodic_argument(*ar.args)
        if isinstance(ar, polar_lift) and period >= 2*pi:
            return periodic_argument(ar.args[0], period)
        if ar.is_Mul:
            newargs = [x for x in ar.args if not x.is_positive]
            if len(newargs) != len(ar.args):
                return periodic_argument(Mul(*newargs), period)
        unbranched = cls._getunbranched(ar)
        if unbranched is None:
            return
        if unbranched.has(periodic_argument, atan2, arg, atan):
            return
        if period == oo:
            return unbranched
        else:
            n = ceiling(unbranched/period - Rational(1, 2))*period
            assert not n.has(ceiling)
            return unbranched - n

    def _eval_is_real(self):
        if self.args[1].is_real and self.args[1].is_positive:
            return True


def unbranched_argument(arg):
    return periodic_argument(arg, oo)


class principal_branch(Function):
    """
    Represent a polar number reduced to its principal branch on a quotient
    of the Riemann surface of the logarithm.

    This is a function of two arguments. The first argument is a polar
    number `z`, and the second one a positive real number of infinity, `p`.
    The result is "z mod exp_polar(I*p)".

    >>> principal_branch(z, oo)
    z
    >>> principal_branch(exp_polar(2*pi*I)*3, 2*pi)
    3*exp_polar(0)
    >>> principal_branch(exp_polar(2*pi*I)*3*z, 2*pi)
    3*principal_branch(z, 2*pi)

    See Also
    ========

    diofant.functions.elementary.exponential.exp_polar
    diofant.functions.elementary.complexes.polar_lift : Lift argument to the Riemann surface of the logarithm
    diofant.functions.elementary.complexes.periodic_argument

    """

    is_polar = True
    is_comparable = False  # cannot always be evalf'd

    @classmethod
    def eval(cls, x, period):
        if isinstance(x, polar_lift):
            return principal_branch(x.args[0], period)
        if period == oo:
            return x
        ub = periodic_argument(x, oo)
        barg = periodic_argument(x, period)
        if ub != barg and not ub.has(periodic_argument) \
                and not barg.has(periodic_argument):
            pl = polar_lift(x)

            def mr(expr):
                if not isinstance(expr, Symbol):
                    return polar_lift(expr)
                return expr
            pl = pl.replace(polar_lift, mr)
            if not pl.has(polar_lift):
                res = exp_polar(I*(barg - ub))*pl
                if not res.is_polar and not res.has(exp_polar):
                    res *= exp_polar(0)
                return res

        if not x.free_symbols:
            c, m = x, ()
        else:
            c, m = x.as_coeff_mul(*x.free_symbols)
        others = []
        for y in m:
            if y.is_positive:
                c *= y
            else:
                others += [y]
        m = tuple(others)
        arg = periodic_argument(c, period)
        if arg.has(periodic_argument):
            return
        if arg.is_number and (unbranched_argument(c) != arg or
                              (arg == 0 and m and c != 1)):
            if arg == 0:
                return abs(c)*principal_branch(Mul(*m), period)
            return principal_branch(exp_polar(I*arg)*Mul(*m), period)*abs(c)
        if arg.is_number and ((abs(arg) - period/2).is_negative or arg == period/2) \
                and not m:
            return exp_polar(arg*I)*abs(c)

    def _eval_evalf(self, prec):
        from .exponential import exp
        z, period = self.args
        p = periodic_argument(z, period)._eval_evalf(prec)
        if p is None or abs(p) > pi or p == -pi:
            return self  # Cannot evalf for this argument.
        return (abs(z)*exp(I*p))._eval_evalf(prec)


def _polarify(eq, lift, pause=False):
    from ...integrals import Integral

    if isinstance(eq, Tuple):
        return eq.func(*[_polarify(arg, lift, pause=False) for arg in eq.args])

    if eq.is_polar:
        return eq
    if isinstance(eq, BooleanAtom):
        return eq
    if eq.is_number and not pause:
        return polar_lift(eq)
    if isinstance(eq, (Dummy, Symbol)) and not pause and lift:
        return polar_lift(eq)
    elif eq.is_Atom:
        return eq
    elif eq.is_Add:
        r = eq.func(*[_polarify(arg, lift, pause=True) for arg in eq.args])
        if lift:
            return polar_lift(r)
        return r
    elif eq.is_Function:
        return eq.func(*[_polarify(arg, lift, pause=False) for arg in eq.args])
    elif eq.is_Exp:
        return eq.func(eq.base, _polarify(eq.exp, lift, pause=False))
    elif isinstance(eq, Integral):
        # Don't lift the integration variable
        func = _polarify(eq.function, lift, pause=pause)
        limits = []
        for limit in eq.args[1:]:
            var = _polarify(limit[0], lift=False, pause=pause)
            rest = tuple(_polarify(x, lift=lift, pause=pause) for x in limit[1:])
            limits.append((var,) + rest)
        return Integral(*((func,) + tuple(limits)))
    else:
        return eq.func(*[_polarify(arg, lift, pause=pause)
                         if isinstance(arg, Expr) else arg for arg in eq.args])


def polarify(eq, subs=True, lift=False):
    """
    Turn all numbers in eq into their polar equivalents (under the standard
    choice of argument).

    Note that no attempt is made to guess a formal convention of adding
    polar numbers, expressions like 1 + x will generally not be altered.

    Note also that this function does not promote exp(x) to exp_polar(x).

    If ``subs`` is True, all symbols which are not already polar will be
    substituted for polar dummies; in this case the function behaves much
    like posify.

    If ``lift`` is True, both addition statements and non-polar symbols are
    changed to their polar_lift()ed versions.
    Note that lift=True implies subs=False.

    >>> expr = (-x)**y
    >>> expr.expand()
    (-x)**y
    >>> polarify(expr)[0]
    (_x*exp_polar(I*pi))**_y
    >>> sorted(polarify(expr)[1].items(), key=default_sort_key)
    [(_x, x), (_y, y)]
    >>> polarify(expr)[0].expand()
    _x**_y*exp_polar(I*pi*_y)
    >>> polarify(x, lift=True)
    polar_lift(x)
    >>> polarify(x*(1+y), lift=True)
    polar_lift(x)*polar_lift(y + 1)

    Adds are treated carefully:

    >>> polarify(1 + sin((1 + I)*x))
    (sin(_x*polar_lift(1 + I)) + 1, {_x: x})

    """
    if lift:
        subs = False
    eq = _polarify(sympify(eq), lift)
    if not subs:
        return eq
    reps = {s: Dummy(s.name, polar=True) for s in eq.free_symbols}
    eq = eq.subs(reps)
    return eq, {r: s for s, r in reps.items()}


def _unpolarify(eq, exponents_only, pause=False):
    if isinstance(eq, bool) or eq.is_Atom:
        return eq

    if not pause:
        if isinstance(eq, exp_polar):
            return exp(_unpolarify(eq.exp, exponents_only))
        if isinstance(eq, principal_branch) and eq.args[1] == 2*pi:
            return _unpolarify(eq.args[0], exponents_only)
        if (
            eq.is_Add or eq.is_Mul or eq.is_Boolean or
            eq.is_Relational and (
                eq.rel_op in ('==', '!=') and 0 in eq.args or
                eq.rel_op not in ('==', '!='))
        ):
            return eq.func(*[_unpolarify(x, exponents_only) for x in eq.args])
        if isinstance(eq, polar_lift):
            return _unpolarify(eq.args[0], exponents_only)

    if eq.is_Pow and not eq.is_Exp:
        expo = _unpolarify(eq.exp, exponents_only)
        base = _unpolarify(eq.base, exponents_only,
                           not (expo.is_integer and not pause))
        return base**expo
    elif eq.is_Exp:
        return exp(_unpolarify(eq.exp, exponents_only, exponents_only))
    elif isinstance(eq, ExprCondPair):
        return eq.func(_unpolarify(eq.expr, exponents_only, exponents_only),
                       _unpolarify(eq.cond, exponents_only, exponents_only))

    if eq.is_Function and getattr(eq.func, 'unbranched', False):
        return eq.func(*[_unpolarify(x, exponents_only, exponents_only)
                         for x in eq.args])

    return eq.func(*[_unpolarify(x, exponents_only, True) for x in eq.args])


def unpolarify(eq, subs={}, exponents_only=False):
    """
    If p denotes the projection from the Riemann surface of the logarithm to
    the complex line, return a simplified version eq' of `eq` such that
    p(eq') == p(eq).
    Also apply the substitution subs in the end. (This is a convenience, since
    ``unpolarify``, in a certain sense, undoes polarify.)

    >>> unpolarify(polar_lift(I + 2))
    2 + I
    >>> unpolarify(sin(polar_lift(I + 7)))
    sin(7 + I)

    """
    if isinstance(eq, bool):
        return eq

    eq = sympify(eq)
    if subs != {}:
        return unpolarify(eq.subs(subs))
    changed = True
    pause = False
    if exponents_only:
        pause = True
    while changed:
        changed = False
        res = _unpolarify(eq, exponents_only, pause)
        if res != eq:
            changed = True
            eq = res
    # Finally, replacing Exp(0) by 1 is always correct.
    # So is polar_lift(0) -> 0.
    return res.subs({exp_polar(0): 1, polar_lift(0): 0})
