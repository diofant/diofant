from mpmath.libmp.libmpf import prec_to_dps

from ...core import (Add, Function, Ge, Gt, I, Integer, Le, Lt,
                     PrecisionExhausted)
from ...logic import false, true


###############################################################################
# ####################### FLOOR and CEILING FUNCTIONS ####################### #
###############################################################################


class RoundFunction(Function):
    """The base class for rounding functions."""

    @classmethod
    def eval(cls, arg):
        from .complexes import im
        if arg.is_integer:
            return arg
        if isinstance(arg, cls):
            return arg
        if arg.is_imaginary or (I*arg).is_extended_real:
            i = im(arg)
            if not i.has(I):
                return cls(i)*I
            return cls(arg, evaluate=False)

        v = cls._eval_number(arg)
        if v is not None:
            return v

        # Integral, numerical, symbolic part
        ipart = npart = spart = Integer(0)

        # Extract integral (or complex integral) terms
        terms = Add.make_args(arg)

        for t in terms:
            if t.is_integer or (t.is_imaginary and im(t).is_integer):
                ipart += t
            elif t.free_symbols:
                spart += t
            else:
                npart += t

        if not (npart or spart):
            return ipart

        # Evaluate npart numerically if independent of spart
        if npart and (not spart or npart.is_extended_real and
                      (spart.is_imaginary or (I*spart).is_extended_real) or
                      npart.is_imaginary and spart.is_extended_real):
            npart_int = None
            try:
                from ...core.evalf import DEFAULT_MAXPREC as TARGET
                prec = 10
                r, i = Integer(0), Integer(0)
                npart_re, npart_im = npart.as_real_imag()
                while prec < TARGET:
                    dps = prec_to_dps(prec)
                    r, i = npart_re.evalf(dps), npart_im.evalf(dps)
                    if ((not r or int(2**prec*abs(r)) > 2**prec*abs(int(r))) and
                            (not i or int(2**prec*abs(i)) > 2**prec*abs(int(i)))):
                        npart_int = cls(r) + cls(i)*I
                        break
                    prec += 10
                else:
                    raise PrecisionExhausted
            except PrecisionExhausted:
                npart_int = cls(r) + cls(i)*I
                if not npart.equals(npart_int):
                    npart_int = None

            if npart_int is not None:
                ipart += npart_int
                npart = Integer(0)

        spart += npart
        if not spart:
            return ipart
        if spart.is_imaginary or (I*spart).is_extended_real:
            return ipart + cls(im(spart), evaluate=False)*I
        return ipart + cls(spart, evaluate=False)

    def _eval_is_finite(self):
        return self.args[0].is_finite

    def _eval_is_extended_real(self):
        if self.args[0].is_extended_real:
            return True

    def _eval_is_integer(self):
        if self.args[0].is_real:
            return True


class floor(RoundFunction):
    """
    Floor is a univariate function which returns the largest integer
    value not greater than its argument. However this implementation
    generalizes floor to complex numbers.

    Examples
    ========

    >>> floor(17)
    17
    >>> floor(Rational(23, 10))
    2
    >>> floor(2*E)
    5
    >>> floor(-Float(0.567))
    -1
    >>> floor(-I/2)
    -I

    See Also
    ========

    diofant.functions.elementary.integers.ceiling

    References
    ==========

    * "Concrete mathematics" by Graham, pp. 87
    * https://mathworld.wolfram.com/FloorFunction.html

    """

    _dir = -1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            if arg.is_Rational:
                return Integer(arg.numerator // arg.denominator)
            if arg.is_Float:
                return Integer(int(arg.floor()))
            return arg
        if isinstance(arg, (floor, ceiling)):
            return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(Integer)[0]

    def _eval_nseries(self, x, n, logx):
        r = self.subs({x: 0})
        args = self.args[0]
        args0 = args.subs({x: 0})
        if args0 == r:
            direction = (args - args0).as_leading_term(x).as_coeff_exponent(x)[0]
            if direction.is_positive:
                return r
            return r - 1
        return r

    def __le__(self, other):
        if self.args[0] == other and other.is_extended_real:
            return true
        return Le(self, other, evaluate=False)

    def __gt__(self, other):
        if self.args[0] == other and other.is_extended_real:
            return false
        return Gt(self, other, evaluate=False)

    def _eval_as_leading_term(self, x):
        return self


class ceiling(RoundFunction):
    """
    Ceiling is a univariate function which returns the smallest integer
    value not less than its argument. Ceiling function is generalized
    in this implementation to complex numbers.

    Examples
    ========

    >>> ceiling(17)
    17
    >>> ceiling(Rational(23, 10))
    3
    >>> ceiling(2*E)
    6
    >>> ceiling(-Float(0.567))
    0
    >>> ceiling(I/2)
    I

    See Also
    ========

    diofant.functions.elementary.integers.floor

    References
    ==========

    * "Concrete mathematics" by Graham, pp. 87
    * https://mathworld.wolfram.com/CeilingFunction.html

    """

    _dir = 1

    @classmethod
    def _eval_number(cls, arg):
        if arg.is_Number:
            if arg.is_Rational:
                return -Integer(-arg.numerator // arg.denominator)
            if arg.is_Float:
                return Integer(int(arg.ceiling()))
            return arg
        if isinstance(arg, (ceiling, floor)):
            return arg
        if arg.is_NumberSymbol:
            return arg.approximation_interval(Integer)[1]

    def _eval_nseries(self, x, n, logx):
        r = self.subs({x: 0})
        args = self.args[0]
        args0 = args.subs({x: 0})
        if args0 == r:
            direction = (args - args0).as_leading_term(x).as_coeff_exponent(x)[0]
            if direction.is_positive:
                return r + 1
            return r
        return r

    def __lt__(self, other):
        if self.args[0] == other and other.is_extended_real:
            return false
        return Lt(self, other, evaluate=False)

    def __ge__(self, other):
        if self.args[0] == other and other.is_extended_real:
            return true
        return Ge(self, other, evaluate=False)
