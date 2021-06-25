import decimal
import fractions
import math
import numbers

import mpmath
import mpmath.libmp as mlib

from .cache import cacheit
from .compatibility import GROUND_TYPES, HAS_GMPY, as_int, gmpy
from .containers import Tuple
from .decorators import _sympifyit
from .expr import AtomicExpr, Expr
from .logic import fuzzy_not
from .singleton import S, SingletonWithManagedProperties
from .sympify import SympifyError, converter, sympify


rnd = mlib.round_nearest


def comp(z1, z2, tol=None):
    """Return a bool indicating whether the error between z1 and z2 is <= tol.

    If ``tol`` is None then True will be returned if there is a significant
    difference between the numbers: ``abs(z1 - z2)*10**p <= 1/2`` where ``p``
    is the lower of the precisions of the values. A comparison of strings will
    be made if ``z1`` is a Number and a) ``z2`` is a string or b) ``tol`` is ''
    and ``z2`` is a Number.

    When ``tol`` is a nonzero value, if z2 is non-zero and ``|z1| > 1``
    the error is normalized by ``|z1|``, so if you want to see if the
    absolute error between ``z1`` and ``z2`` is <= ``tol`` then call this
    as ``comp(z1 - z2, 0, tol)``.

    """
    if type(z2) is str:
        if not isinstance(z1, Number):
            raise ValueError('when z2 is a str z1 must be a Number')
        return str(z1) == z2
    if not z1:
        z1, z2 = z2, z1
    if not z1:
        return True
    if not tol:
        if tol is None:
            a, b = Float(z1), Float(z2)
            return int(abs(a - b)*10**mlib.prec_to_dps(
                min(a._prec, b._prec)))*2 <= 1
        elif all(getattr(i, 'is_Number', False) for i in (z1, z2)):
            return z1._prec == z2._prec and str(z1) == str(z2)
        raise ValueError('exact comparison requires two Numbers')
    diff = abs(z1 - z2)
    az1 = abs(z1)
    if z2 and az1 > 1:
        return diff/az1 <= tol
    else:
        return diff <= tol


def mpf_norm(mpf, prec):
    """Return the mpf tuple normalized appropriately for the indicated
    precision after doing a check to see if zero should be returned or
    not when the mantissa is 0. ``mlib.normalize`` always assumes that this
    is zero, but it may not be since the mantissa for mpf's values "+inf",
    "-inf" and "nan" have a mantissa of zero, too.

    Note: this is not intended to validate a given mpf tuple, so sending
    mpf tuples that were not created by mpmath may produce bad results. This
    is only a wrapper to ``mpmath.libmp.normalize`` which provides the check
    for non-zero mpfs that have a 0 for the mantissa.

    """
    sign, man, expt, bc = mpf
    if not man:
        # hack for normalize which does not do this;
        # it assumes that if man is zero the result is 0
        # (see issue sympy/sympy#6639)
        if not bc:
            return mlib.fzero
        else:
            # don't change anything; this should already
            # be a well formed mpf tuple
            return mpf
    rv = mlib.normalize(sign, man, expt, bc, prec, rnd)
    return rv


def _str_to_Decimal_dps(s):
    """Convert a string to pair of a Decimal instance and its precision."""
    try:
        num = decimal.Decimal(s)
    except decimal.InvalidOperation:
        raise ValueError(f'string-float not recognized: {s}')
    else:
        return num, len(num.as_tuple().digits)


def igcdex(a, b):
    """Returns x, y, g such that g = x*a + y*b = gcd(a, b).

    >>> igcdex(2, 3)
    (-1, 1, 1)
    >>> igcdex(10, 12)
    (-1, 1, 2)

    >>> igcdex(100, 2004)
    (-20, 1, 4)
    >>> x, y = _[:-1]
    >>> x*100 + y*2004
    4

    """
    if (not a) and (not b):
        return 0, 1, 0

    if not a:
        return 0, b//abs(b), abs(b)
    if not b:
        return a//abs(a), 0, abs(a)

    if a < 0:
        a, x_sign = -a, -1
    else:
        x_sign = 1

    if b < 0:
        b, y_sign = -b, -1
    else:
        y_sign = 1

    x, y, r, s = 1, 0, 0, 1

    while b:
        c, q = a % b, a // b
        a, b, r, s, x, y = b, c, x - q*r, y - q*s, r, s

    return x*x_sign, y*y_sign, a


def mod_inverse(a, m):
    """
    Return the number c such that, ( a * c ) % m == 1 where
    c has the same sign as a. If no such value exists, a
    ValueError is raised.

    Examples
    ========

    Suppose we wish to find multiplicative inverse x of
    3 modulo 11. This is the same as finding x such
    that 3 * x = 1 (mod 11). One value of x that satisfies
    this congruence is 4. Because 3 * 4 = 12 and 12 = 1 mod(11).
    This is the value return by mod_inverse:

    >>> mod_inverse(3, 11)
    4
    >>> mod_inverse(-3, 11)
    -4

    When there is a common factor between the numerators of
    ``a`` and ``m`` the inverse does not exist:

    >>> mod_inverse(2, 4)
    Traceback (most recent call last):
    ...
    ValueError: inverse of 2 mod 4 does not exist

    >>> mod_inverse(Integer(2)/7, Integer(5)/2)
    7/2

    References
    ==========

    * https://en.wikipedia.org/wiki/Modular_multiplicative_inverse
    * https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm

    """
    from ..utilities import filldedent
    c = None
    try:
        a, m = as_int(a), as_int(m)
        if m > 1:
            x, y, g = igcdex(a, m)
            if g == 1:
                c = x % m
            if a < 0:
                c -= m
    except ValueError:
        a, m = sympify(a), sympify(m)
        if not (a.is_number and m.is_number):
            raise TypeError(filldedent("""
                Expected numbers for arguments; symbolic `mod_inverse`
                is not implemented
                but symbolic expressions can be handled with the
                similar function,
                sympy.polys.polytools.invert"""))
        big = (m > 1)
        if not (big is S.true or big is S.false):
            raise ValueError(f'm > 1 did not evaluate; try to simplify {m}')
        elif big:
            c = 1/a
    if c is None:
        raise ValueError(f'inverse of {a} (mod {m}) does not exist')
    return c


def integer_digits(n, b):
    """Gives a list of the base `b` digits in the integer `n`."""
    n, b = map(as_int, (n, b))
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(n % b)
        n //= b
    return digits[::-1]


class Number(AtomicExpr):
    """
    Represents any kind of number in diofant.

    Floating point numbers are represented by the Float class.
    Integer numbers (of any size), together with rational numbers (again,
    there is no limit on their size) are represented by the Rational class.

    """

    is_commutative = True
    is_number = True
    is_Number = True

    # Used to make max(x._prec, y._prec) return x._prec when only x is a float
    _prec = -1

    def __new__(cls, *obj):
        if len(obj) == 1:
            obj = obj[0]

        if isinstance(obj, Number):
            return obj
        if isinstance(obj, numbers.Integral):
            return Integer(obj)
        if isinstance(obj, tuple) and len(obj) == 2:
            return Rational(*obj)
        if isinstance(obj, (float, mpmath.mpf, decimal.Decimal)):
            return Float(obj)
        if isinstance(obj, str):
            val = sympify(obj)
            if isinstance(val, Number):
                return val
            else:
                raise ValueError(f'String "{obj}" does not denote a Number')
        msg = 'expected str|int|float|Decimal|Number object but got %r'
        raise TypeError(msg % type(obj).__name__)

    def invert(self, other, *gens, **args):
        from ..polys.polytools import invert
        if getattr(other, 'is_number', True):
            return mod_inverse(self, other)
        return invert(self, other, *gens, **args)

    @_sympifyit('other', NotImplemented)
    def __divmod__(self, other):
        if not other:
            raise ZeroDivisionError('modulo by zero')
        elif self.is_Integer and other.is_Integer:
            return Tuple(*divmod(self.numerator, other.numerator))
        else:
            rat = self/other
            w = math.floor(rat)
            r = self - other*w
            return Tuple(w, r)

    @_sympifyit('other', NotImplemented)
    def __rdivmod__(self, other):
        return divmod(other, self)

    def __round__(self, *args):
        return round(float(self), *args)

    def _as_mpf_val(self, prec):  # pragma: no cover
        """Evaluation of mpf tuple accurate to at least prec bits."""
        raise NotImplementedError(f'{self.__class__.__name__} needs ._as_mpf_val() method')

    def _eval_evalf(self, prec):
        return Float._new(self._as_mpf_val(prec), prec)

    def _as_mpf_op(self, prec):
        prec = max(prec, self._prec)
        return self._as_mpf_val(prec), prec

    def __float__(self):
        return mlib.to_float(self._as_mpf_val(53))

    def _eval_conjugate(self):
        return self

    def _eval_subs(self, old, new):
        if old == -self:
            return -new
        return self  # there is no other possibility

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 1, 0, 'Number'

    @cacheit
    def sort_key(self, order=None):
        """Return a sort key."""
        return self.class_key(), (0, ()), (), self

    def __hash__(self):
        return super().__hash__()

    def is_constant(self, *wrt, **flags):
        """Return True if self is constant.

        See Also
        ========

        diofant.core.expr.Expr.is_constant

        """
        return True

    def as_coeff_mul(self, *deps, **kwargs):
        """Return the tuple (c, args) where self is written as a Mul.

        See Also
        ========

        diofant.core.expr.Expr.as_coeff_mul

        """
        # a -> c*t
        if self.is_Rational or not kwargs.pop('rational', True):
            return self, ()
        elif self.is_negative:
            return S.NegativeOne, (-self,)
        return S.One, (self,)

    def as_coeff_add(self, *deps):
        """Return the tuple (c, args) where self is written as an Add.

        See Also
        ========

        diofant.core.expr.Expr.as_coeff_add

        """
        # a -> c + t
        if self.is_Rational:
            return self, ()
        return S.Zero, (self,)

    def as_coeff_Mul(self, rational=False):
        """Efficiently extract the coefficient of a product."""
        if rational and not self.is_Rational:
            return S.One, self
        return (self, S.One) if self else (S.One, self)

    def as_coeff_Add(self, rational=False):
        """Efficiently extract the coefficient of a summation."""
        if not rational:
            return self, S.Zero
        return S.Zero, self

    def gcd(self, other):
        """Compute GCD of `self` and `other`."""
        from ..polys import gcd
        return gcd(self, other)

    def lcm(self, other):
        """Compute LCM of `self` and `other`."""
        from ..polys import lcm
        return lcm(self, other)

    def cofactors(self, other):
        """Compute GCD and cofactors of `self` and `other`."""
        from ..polys import cofactors
        return cofactors(self, other)


class Float(Number):
    """Represent a floating-point number of arbitrary precision.

    Examples
    ========

    >>> Float(3.5)
    3.50000000000000
    >>> Float(3)
    3.

    Creating Floats from strings (and Python ``int`` type) will
    give a minimum precision of 15 digits, but the precision
    will automatically increase to capture all digits entered.

    >>> Float(1)
    1.
    >>> Float(10**20)
    100000000000000000000.
    >>> Float('1e20')
    1.e+20

    However, *floating-point* numbers (Python ``float`` types) retain
    only 15 digits of precision:

    >>> Float(1e20)
    1.00000000000000e+20
    >>> Float(1.23456789123456789)
    1.23456789123457

    It may be preferable to enter high-precision decimal numbers
    as strings:

    >>> Float('1.23456789123456789')
    1.23456789123456789

    The desired number of digits can also be specified:

    >>> Float('1e-3', 3)
    0.00100
    >>> Float(100, 4)
    100.0

    Float can automatically count significant figures if decimal precision
    argument is omitted. (Auto-counting is only allowed for strings and ints).

    >>> Float('12e-3')
    0.012
    >>> Float(3)
    3.
    >>> Float('60.e2')  # 2 digits significant
    6.0e+3
    >>> Float('6000.')  # 4 digits significant
    6000.
    >>> Float('600e-2')  # 3 digits significant
    6.00

    Notes
    =====

    Floats are inexact by their nature unless their value is a binary-exact
    value.

    >>> approx, exact = Float(.1, 1), Float(.125, 1)

    For calculation purposes, you can change the precision of Float,
    but this will not increase the accuracy of the inexact value. The
    following is the most accurate 5-digit approximation of a value of 0.1
    that had only 1 digit of precision:

    >>> Float(approx, 5)
    0.099609

    Please note that you can't increase precision with evalf:

    >>> approx.evalf(5)
    Traceback (most recent call last):
    ...
    PrecisionExhausted: ...

    By contrast, 0.125 is exact in binary (as it is in base 10) and so it
    can be passed to Float constructor to obtain an arbitrary precision with
    matching accuracy:

    >>> Float(exact, 5)
    0.12500
    >>> Float(exact, 20)
    0.12500000000000000000

    Trying to make a high-precision Float from a float is not disallowed,
    but one must keep in mind that the *underlying float* (not the apparent
    decimal value) is being obtained with high precision. For example, 0.3
    does not have a finite binary representation. The closest rational is
    the fraction 5404319552844595/2**54. So if you try to obtain a Float of
    0.3 to 20 digits of precision you will not see the same thing as 0.3
    followed by 19 zeros:

    >>> Float(0.3, 20)
    0.29999999999999998890

    If you want a 20-digit value of the decimal 0.3 (not the floating point
    approximation of 0.3) you should send the 0.3 as a string. The underlying
    representation is still binary but a higher precision than Python's float
    is used:

    >>> Float('0.3', 20)
    0.30000000000000000000

    Although you can increase the precision of an existing Float using Float
    it will not increase the accuracy -- the underlying value is not changed:

    >>> def show(f):  # binary rep of Float
    ...     from diofant import Mul, Pow
    ...     s, m, e, b = f._mpf_
    ...     v = Mul(int(m), Pow(2, int(e), evaluate=False), evaluate=False)
    ...     print(f'{v} at prec={f._prec}')
    ...
    >>> t = Float('0.3', 3)
    >>> show(t)
    4915/2**14 at prec=13
    >>> show(Float(t, 20))  # higher prec, not higher accuracy
    4915/2**14 at prec=70
    >>> show(Float(t, 2))  # lower prec
    307/2**10 at prec=10

    """

    is_number = True

    is_extended_real = True

    is_Float = True

    def __new__(cls, num, dps=None):
        if dps is None:
            if isinstance(num, Float):
                return num
            elif isinstance(num, (str, numbers.Integral)):
                num, dps = _str_to_Decimal_dps(str(num))
            else:
                dps = 15

        prec = mlib.libmpf.dps_to_prec(dps)

        if isinstance(num, decimal.Decimal):
            _mpf_ = mlib.from_Decimal(num, prec, rnd)
        elif isinstance(num, Number):
            _mpf_ = num._as_mpf_val(prec)
        else:
            _mpf_ = mpmath.mpf(num, prec=prec, rounding=rnd)._mpf_

        # special cases
        if _mpf_ == mlib.fzero:
            pass  # we want a Float
        elif _mpf_ == mlib.fnan:
            return nan

        obj = Expr.__new__(cls)
        obj._mpf_ = _mpf_
        obj._prec = prec
        return obj

    @classmethod
    def _new(cls, _mpf_, _prec):
        # special cases
        if _mpf_ == mlib.fzero:
            return S.Zero  # XXX this is different from Float which gives 0.0
        elif _mpf_ == mlib.fnan:
            return nan

        obj = Expr.__new__(cls)
        obj._mpf_ = mpf_norm(_mpf_, _prec)
        obj._prec = _prec
        return obj

    def __getnewargs__(self):
        return self._mpf_, mlib.libmpf.prec_to_dps(self._prec)

    def _hashable_content(self):
        return self._mpf_, self._prec

    def floor(self):
        """Compute floor of self."""
        return Integer(mlib.to_int(mlib.mpf_floor(self._mpf_, self._prec)))

    def ceiling(self):
        """Compute ceiling of self."""
        return Integer(mlib.to_int(mlib.mpf_ceil(self._mpf_, self._prec)))

    @property
    def num(self):
        """Return mpmath representation."""
        return mpmath.mpf(self._mpf_)

    def _as_mpf_val(self, prec):
        return mpf_norm(self._mpf_, prec)

    def _as_mpf_op(self, prec):
        return self._mpf_, max(prec, self._prec)

    def _eval_is_finite(self):
        return self._mpf_ not in (mlib.finf, mlib.fninf)

    def _eval_is_integer(self):
        return self._mpf_ == mlib.fzero

    def _eval_is_positive(self):
        return mlib.mpf_gt(self._mpf_, mlib.fzero)

    def _eval_is_zero(self):
        return self._mpf_ == mlib.fzero

    def __bool__(self):
        return self.is_nonzero

    def __neg__(self):
        return Float._new(mlib.mpf_neg(self._mpf_), self._prec)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_add(self._mpf_, rhs, prec, rnd), prec)
        return Number.__add__(self, other)
    __radd__ = __add__

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_sub(self._mpf_, rhs, prec, rnd), prec)
        return Number.__sub__(self, other)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_mul(self._mpf_, rhs, prec, rnd), prec)
        return Number.__mul__(self, other)
    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        if isinstance(other, Number) and other != 0:
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_div(self._mpf_, rhs, prec, rnd), prec)
        return Number.__truediv__(self, other)

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if isinstance(other, Rational) and other.denominator != 1:
            # calculate mod with Rationals, *then* round the result
            return Float(Rational.__mod__(Rational(self), other),
                         mlib.prec_to_dps(self._prec))
        if isinstance(other, Float):
            r = self/other
            if r == int(r):
                prec = max(mlib.prec_to_dps(i) for i in (self._prec, other._prec))
                return Float(0, prec)
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_mod(self._mpf_, rhs, prec, rnd), prec)
        return Number.__mod__(self, other)

    @_sympifyit('other', NotImplemented)
    def __rmod__(self, other):
        if isinstance(other, Float):
            return other.__mod__(self)
        elif isinstance(other, Rational):
            # calculate mod with Rationals, *then* round the answer
            return Float(other.__mod__(Rational(self)),
                         mlib.prec_to_dps(self._prec))
        else:
            return NotImplemented

    def _eval_power(self, expt):
        """
        Expt is symbolic object but not equal to 0, 1.

        (-p)**r -> exp(r*log(-p)) -> exp(r*(log(p) + I*Pi)) ->
                  -> p**r*(sin(Pi*r) + cos(Pi*r)*I)

        """
        from .power import Pow
        if self == 0:
            if expt.is_positive:
                return S.Zero
            if expt.is_negative:
                return Float('inf')
        if isinstance(expt, Number):
            if isinstance(expt, Integer):
                prec = self._prec
                return Float._new(
                    mlib.mpf_pow_int(self._mpf_, expt.numerator, prec, rnd), prec)
            elif isinstance(expt, Rational) and \
                    expt.numerator == 1 and expt.denominator % 2 and self.is_negative:
                return Pow(S.NegativeOne, expt, evaluate=False)*(
                    -self)._eval_power(expt)
            expt, prec = expt._as_mpf_op(self._prec)
            mpfself = self._mpf_
            try:
                y = mlib.mpf_pow(mpfself, expt, prec, rnd)
                return Float._new(y, prec)
            except mlib.ComplexResult:
                re, im = mlib.mpc_pow(
                    (mpfself, mlib.fzero), (expt, mlib.fzero), prec, rnd)
                return Float._new(re, prec) + \
                    Float._new(im, prec)*I

    def __abs__(self):
        return Float._new(mlib.mpf_abs(self._mpf_), self._prec)

    def __int__(self):
        return int(mlib.to_int(self._mpf_))  # uses round_fast = round_down

    def __eq__(self, other):
        if isinstance(other, float):
            # coerce to Float at same precision
            o = Float(other)
            ompf = o._as_mpf_val(self._prec)
            return bool(mlib.mpf_eq(self._mpf_, ompf))
        try:
            other = sympify(other, strict=True)
        except SympifyError:
            return False    # diofant != other  -->  not ==
        if isinstance(other, NumberSymbol):
            if other.is_irrational:
                return False
            return other.__eq__(self)
        if isinstance(other, Float):
            return bool(mlib.mpf_eq(self._mpf_, other._mpf_))
        if isinstance(other, Number):
            # numbers should compare at the same precision;
            # all _as_mpf_val routines should be sure to abide
            # by the request to change the prec if necessary; if
            # they don't, the equality test will fail since it compares
            # the mpf tuples
            ompf = other._as_mpf_val(self._prec)
            return bool(mlib.mpf_eq(self._mpf_, ompf))
        return False    # Float != non-Number

    @_sympifyit('other', NotImplemented)
    def __gt__(self, other):
        if isinstance(other, NumberSymbol):
            return other.__lt__(self)
        if other.is_comparable:
            other = other.evalf(strict=False)
        if isinstance(other, Number) and other is not nan:
            return sympify(bool(mlib.mpf_gt(self._mpf_,
                                            other._as_mpf_val(self._prec))),
                           strict=True)
        return Expr.__gt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __ge__(self, other):
        if isinstance(other, NumberSymbol):
            return other.__le__(self)
        if other.is_comparable:
            other = other.evalf(strict=False)
        if isinstance(other, Number) and other is not nan:
            return sympify(bool(mlib.mpf_ge(self._mpf_,
                                            other._as_mpf_val(self._prec))),
                           strict=True)
        return Expr.__ge__(self, other)

    @_sympifyit('other', NotImplemented)
    def __lt__(self, other):
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        if other.is_extended_real and other.is_number:
            other = other.evalf(strict=False)
        if isinstance(other, Number) and other is not nan:
            return sympify(bool(mlib.mpf_lt(self._mpf_,
                                            other._as_mpf_val(self._prec))),
                           strict=True)
        return Expr.__lt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __le__(self, other):
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        if other.is_extended_real and other.is_number:
            other = other.evalf(strict=False)
        if isinstance(other, Number) and other is not nan:
            return sympify(bool(mlib.mpf_le(self._mpf_,
                                            other._as_mpf_val(self._prec))),
                           strict=True)
        return Expr.__le__(self, other)

    def __hash__(self):
        return super().__hash__()

    def epsilon_eq(self, other, epsilon='1e-15'):
        """Test approximate equality."""
        return abs(self - other) < Float(epsilon)

    def __format__(self, format_spec):
        return format(decimal.Decimal(str(self)), format_spec)


# Ground type for components of Rational
_int_dtype = gmpy.mpz if GROUND_TYPES == 'gmpy' else int


class Rational(Number):
    """Represents integers and rational numbers (p/q) of any size.

    Examples
    ========

    >>> Rational(3)
    3
    >>> Rational(1, 2)
    1/2

    Rational is unprejudiced in accepting input. If a float is passed, the
    underlying value of the binary representation will be returned:

    >>> Rational(.5)
    1/2
    >>> Rational(.2)
    3602879701896397/18014398509481984

    If the simpler representation of the float is desired then consider
    limiting the denominator to the desired value or convert the float to
    a string (which is roughly equivalent to limiting the denominator to
    10**12):

    >>> Rational(str(.2))
    1/5
    >>> Rational(.2).limit_denominator(10**12)
    1/5

    An arbitrarily precise Rational is obtained when a string literal is
    passed:

    >>> Rational('1.23')
    123/100
    >>> Rational('1e-2')
    1/100
    >>> Rational('.1')
    1/10

    The conversion of floats to expressions or simple fractions can
    be handled with nsimplify:

    >>> nsimplify(.3)  # numbers that have a simple form
    3/10

    But if the input does not reduce to a literal Rational, an error will
    be raised:

    >>> Rational(pi)
    Traceback (most recent call last):
    ...
    TypeError: invalid input: pi


    Low-level access numerator and denominator:

    >>> r = Rational(3, 4)
    >>> r
    3/4
    >>> r.numerator
    3
    >>> r.denominator
    4

    Note that these properties return integers (not Diofant Integers) so some care
    is needed when using them in expressions:

    >>> r.numerator/r.denominator
    0.75

    See Also
    ========

    diofant.core.sympify.sympify
    diofant.simplify.simplify.nsimplify

    """

    is_real = True
    is_integer = False
    is_rational = True
    is_number = True

    is_Rational = True

    @cacheit
    def __new__(cls, p, q=1):
        if q == 1:
            if isinstance(p, Rational):
                return p
            elif isinstance(p, Float):
                with mpmath.workprec(p._prec):
                    p, q = mlib.to_rational(p._mpf_)

        try:
            f = fractions.Fraction(p)/fractions.Fraction(q)
            p, q = f.numerator, f.denominator
        except ValueError:
            raise TypeError(f'invalid input: {p}, {q}')
        except ZeroDivisionError:
            if p == 0:
                return nan
            else:
                return zoo

        if q == 1:
            return Integer(p)
        if p == 1 and q == 2:
            return S.Half

        obj = Expr.__new__(cls)

        obj._numerator = _int_dtype(p)
        obj._denominator = _int_dtype(q)

        return obj

    def limit_denominator(self, max_denominator=1000000):
        """Closest Rational to self with denominator at most max_denominator.

        >>> Rational('3.141592653589793').limit_denominator(10)
        22/7
        >>> Rational('3.141592653589793').limit_denominator(100)
        311/99

        """
        f = fractions.Fraction(self.numerator, self.denominator)
        return Rational(f.limit_denominator(max_denominator))

    def __getnewargs__(self):
        return self.numerator, self.denominator

    def _hashable_content(self):
        return self.numerator, self.denominator

    def _eval_is_positive(self):
        return self.numerator > 0

    def _eval_is_zero(self):
        return self.numerator == 0

    def __bool__(self):
        return self.is_nonzero

    def __neg__(self):
        return Rational(-self.numerator, self.denominator)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Rational):
            n, d = other.numerator, other.denominator
            return Rational(self.numerator*d + self.denominator*n, self.denominator*d)
        else:
            return other.__radd__(self)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Rational):
            n, d = other.numerator, other.denominator
            return Rational(self.numerator*d - self.denominator*n, self.denominator*d)
        else:
            return other.__rsub__(self)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Rational):
            n, d = other.numerator, other.denominator
            return Rational(self.numerator*n, self.denominator*d)
        else:
            return other.__rmul__(self)

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        if isinstance(other, Rational):
            n, d = other.numerator, other.denominator
            return Rational(self.numerator*d, self.denominator*n)
        else:
            return other.__rtruediv__(self)

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if isinstance(other, Rational):
            n = (self.numerator*other.denominator) // (other.numerator*self.denominator)
            return Rational(self.numerator*other.denominator - n*other.numerator*self.denominator, self.denominator*other.denominator)
        else:
            return other.__rmod__(self)

    @_sympifyit('other', NotImplemented)
    def __rmod__(self, other):
        if isinstance(other, Rational):
            return Rational.__mod__(other, self)
        return Number.__rmod__(self, other)

    def _eval_power(self, expt):
        if isinstance(expt, Number):
            if expt.is_negative:
                # (3/4)**-2 -> (4/3)**2
                ne = -expt
                if (ne is S.One):
                    return Rational(self.denominator, self.numerator)
                if self.is_negative:
                    return -((S.NegativeOne)**((expt.numerator % expt.denominator) /
                                               Integer(expt.denominator)) *
                             Rational(self.denominator, -self.numerator)**ne)
                else:
                    return Rational(self.denominator, self.numerator)**ne
            if expt is oo:  # -oo already caught by test for negative
                if self.numerator > self.denominator:
                    # (3/2)**oo -> oo
                    return oo
                if self.numerator < -self.denominator:
                    # (-3/2)**oo -> oo + I*oo
                    return oo + oo*I
                return S.Zero
            if isinstance(expt, Float):
                return self._eval_evalf(expt._prec)**expt
            elif isinstance(expt, Integer):
                # (4/3)**2 -> 4**2 / 3**2
                return Rational(self.numerator**expt.numerator, self.denominator**expt.numerator)
            else:  # Rational
                if self.numerator != 1:
                    # (4/3)**(5/6) -> 4**(5/6)*3**(-5/6)
                    return Integer(self.numerator)**expt*Integer(self.denominator)**(-expt)
                # as the above caught negative self.numerator, now self is positive
                return Integer(self.denominator)**Rational(
                    expt.numerator*(expt.denominator - 1), expt.denominator) / \
                    Integer(self.denominator)**Integer(expt.numerator)

    def _as_mpf_val(self, prec):
        return mlib.from_rational(self.numerator, self.denominator, prec, rnd)

    def _mpmath_(self, prec, rnd):
        return mpmath.make_mpf(mlib.from_rational(self.numerator, self.denominator, prec, rnd))

    def __abs__(self):
        return Rational(abs(self.numerator), self.denominator)

    def __int__(self):
        p, q = self.numerator, self.denominator
        if p < 0:
            return -int(-p//q)
        return int(p//q)

    @_sympifyit('other', NotImplemented)
    def __eq__(self, other):
        if isinstance(other, NumberSymbol):
            if other.is_irrational:
                return False
            return other.__eq__(self)
        if isinstance(other, Number):
            if isinstance(other, Rational):
                # a Rational is always in reduced form so will never be 2/4
                # so we can just check equivalence of args
                return self.numerator == other.numerator and self.denominator == other.denominator
            if isinstance(other, Float):
                return mlib.mpf_eq(self._as_mpf_val(other._prec), other._mpf_)
        return False

    @_sympifyit('other', NotImplemented)
    def __gt__(self, other):
        if isinstance(other, NumberSymbol):
            return other.__lt__(self)
        expr = self
        if isinstance(other, Number):
            if isinstance(other, Rational):
                return sympify(bool(self.numerator*other.denominator > self.denominator*other.numerator),
                               strict=True)
            if isinstance(other, Float):
                return sympify(bool(mlib.mpf_gt(self._as_mpf_val(other._prec),
                                                other._mpf_)),
                               strict=True)
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.numerator), self.denominator*other
        return Expr.__gt__(expr, other)

    @_sympifyit('other', NotImplemented)
    def __ge__(self, other):
        if isinstance(other, NumberSymbol):
            return other.__le__(self)
        expr = self
        if isinstance(other, Number):
            if isinstance(other, Rational):
                return sympify(bool(self.numerator*other.denominator >= self.denominator*other.numerator),
                               strict=True)
            if isinstance(other, Float):
                return sympify(bool(mlib.mpf_ge(self._as_mpf_val(other._prec),
                                                other._mpf_)),
                               strict=True)
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.numerator), self.denominator*other
        return Expr.__ge__(expr, other)

    @_sympifyit('other', NotImplemented)
    def __lt__(self, other):
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        expr = self
        if isinstance(other, Number):
            if isinstance(other, Rational):
                return sympify(bool(self.numerator*other.denominator < self.denominator*other.numerator),
                               strict=True)
            if isinstance(other, Float):
                return sympify(bool(mlib.mpf_lt(self._as_mpf_val(other._prec),
                                                other._mpf_)),
                               strict=True)
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.numerator), self.denominator*other
        return Expr.__lt__(expr, other)

    @_sympifyit('other', NotImplemented)
    def __le__(self, other):
        expr = self
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        elif isinstance(other, Number):
            if isinstance(other, Rational):
                return sympify(bool(self.numerator*other.denominator <= self.denominator*other.numerator),
                               strict=True)
            if isinstance(other, Float):
                return sympify(bool(mlib.mpf_le(self._as_mpf_val(other._prec),
                                                other._mpf_)),
                               strict=True)
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.numerator), self.denominator*other
        return Expr.__le__(expr, other)

    def __hash__(self):
        return super().__hash__()

    def factors(self, limit=None, use_trial=True, use_rho=False,
                use_pm1=False, verbose=False, visual=False):
        """A wrapper to factorint which return factors of self that are
        smaller than limit (or cheap to compute). Special methods of
        factoring are disabled by default so that only trial division is used.

        """
        from ..ntheory import factorrat

        return factorrat(self, limit=limit, use_trial=use_trial,
                         use_rho=use_rho, use_pm1=use_pm1,
                         verbose=verbose, visual=visual).copy()

    @_sympifyit('other', NotImplemented)
    def gcd(self, other):
        """Compute GCD of `self` and `other`."""
        if isinstance(other, Rational):
            return Rational(
                Integer(math.gcd(self.numerator, other.numerator)),
                Integer(math.lcm(self.denominator, other.denominator)))
        return Number.gcd(self, other)

    @_sympifyit('other', NotImplemented)
    def lcm(self, other):
        """Compute LCM of `self` and `other`."""
        if isinstance(other, Rational):
            return Rational(
                self.numerator*other.numerator//math.gcd(self.numerator, other.numerator),
                math.gcd(self.denominator, other.denominator))
        return Number.lcm(self, other)

    def _eval_as_numer_denom(self):
        """Expression -> a/b -> a, b.

        See Also
        ========

        diofant.core.expr.Expr.as_numer_denom

        """
        return Integer(self.numerator), Integer(self.denominator)

    def as_content_primitive(self, radical=False):
        """Return the tuple (R, self/R) where R is the positive Rational
        extracted from self.

        Examples
        ========

        >>> Rational(-3, 2).as_content_primitive()
        (3/2, -1)

        See Also
        ========

        diofant.core.expr.Expr.as_content_primitive

        """
        if self:
            if self.is_positive:
                return self, S.One
            return -self, S.NegativeOne
        return S.One, self

    @property
    def numerator(self):
        return self._numerator

    @property
    def denominator(self):
        return self._denominator


numbers.Rational.register(Rational)


class Integer(Rational):
    """Represents integer numbers."""

    is_integer = True
    is_number = True

    is_Integer = True

    _denominator = _int_dtype(1)

    def _as_mpf_val(self, prec):
        return mlib.from_int(self.numerator, prec)

    def _mpmath_(self, prec, rnd):
        return mpmath.make_mpf(self._as_mpf_val(prec))

    @cacheit
    def __new__(cls, i):
        try:
            i = _int_dtype(i)
        except TypeError:
            raise TypeError('Integer can only work with integer expressions.')

        if i == 0:
            return S.Zero
        elif i == 1:
            return S.One
        elif i == -1:
            return S.NegativeOne
        else:
            obj = Expr.__new__(cls)
            obj._numerator = i
            return obj

    def __getnewargs__(self):
        return self.numerator,

    def __hash__(self):
        return hash(self.numerator)

    def __index__(self):
        return int(self.numerator)

    def __format__(self, format_spec):
        return int(self.numerator).__format__(format_spec)

    @_sympifyit('other', NotImplemented)
    def __eq__(self, other):
        if isinstance(other, Integer):
            return self.numerator == other.numerator
        return Rational.__eq__(self, other)

    ########################################

    def _eval_is_odd(self):
        return bool(self.numerator % 2)

    def _eval_power(self, expt):
        """
        Tries to do some simplifications on self**expt

        Returns None if no further simplifications can be done

        When exponent is a fraction (so we have for example a square root),
        we try to find a simpler representation by factoring the argument
        up to factors of 2**15, e.g.

          - sqrt(4) becomes 2
          - sqrt(-4) becomes 2*I
          - root(2**(3+7)*3**(6+7), 7) becomes 6*18**(3/7)

        Further simplification would require a special call to factorint on
        the argument which is not done here for sake of speed.

        """
        from ..ntheory import perfect_power
        from .power import Pow, integer_nthroot

        if expt is oo:
            if self.numerator > S.One:
                return oo
            # cases -1, 0, 1 are done in their respective classes
            return oo + I*oo
        if expt == -oo:
            return Rational(1, self)**oo
        if isinstance(expt, Float):
            # Rational knows how to exponentiate by a Float
            return super()._eval_power(expt)
        if not isinstance(expt, Rational):
            return
        if expt is S.Half and self.is_negative:
            # we extract I for this special case since everyone is doing so
            return I*Pow(-self, expt)
        if expt.is_negative:
            # invert base and change sign on exponent
            ne = -expt
            if self.is_negative:
                return -((S.NegativeOne)**((expt.numerator % expt.denominator) /
                                           Integer(expt.denominator))*Rational(1, -self)**ne)
            else:
                return Rational(1, self.numerator)**ne
        # see if base is a perfect root, sqrt(4) --> 2
        x, xexact = integer_nthroot(abs(self.numerator), expt.denominator)
        if xexact:
            # if it's a perfect root we've finished
            result = Integer(x**abs(expt.numerator))
            if self.is_negative:
                result *= S.NegativeOne**expt
            return result

        # The following is an algorithm where we collect perfect roots
        # from the factors of base.

        # if it's not an nth root, it still might be a perfect power
        b_pos = int(abs(self.numerator))
        p = perfect_power(b_pos)
        if p is not False:
            dict = {p[0]: p[1]}
        else:
            dict = Integer(self).factors(limit=2**15)

        # now process the dict of factors
        if self.is_negative:
            dict[-1] = 1
        out_int = 1  # integer part
        out_rad = 1  # extracted radicals
        sqr_int = 1
        sqr_gcd = 0
        sqr_dict = {}
        for prime, exponent in dict.items():
            exponent *= expt.numerator
            # remove multiples of expt.denominator: (2**12)**(1/10) -> 2*(2**2)**(1/10)
            div_e, div_m = divmod(exponent, expt.denominator)
            if div_e > 0:
                out_int *= prime**div_e
            if div_m > 0:
                # see if the reduced exponent shares a gcd with e.denominator
                # (2**2)**(1/10) -> 2**(1/5)
                g = math.gcd(div_m, expt.denominator)
                if g != 1:
                    out_rad *= Pow(prime, Rational(div_m//g, expt.denominator//g))
                else:
                    sqr_dict[prime] = div_m
        # identify gcd of remaining powers
        for p, ex in sqr_dict.items():
            if sqr_gcd == 0:
                sqr_gcd = ex
            else:
                sqr_gcd = math.gcd(sqr_gcd, ex)
                if sqr_gcd == 1:
                    break
        for k, v in sqr_dict.items():
            sqr_int *= k**(v//sqr_gcd)
        if sqr_int == self and out_int == 1 and out_rad == 1:
            result = None
        else:
            result = out_int*out_rad*Pow(sqr_int, Rational(sqr_gcd, expt.denominator))
        return result

    def _eval_is_prime(self):
        from ..ntheory import isprime

        return isprime(self)

    def _eval_is_composite(self):
        if self > 1:
            return fuzzy_not(self.is_prime)
        else:
            return False

    def __floordiv__(self, other):
        return Integer(self.numerator // Integer(other).numerator)

    def __rfloordiv__(self, other):
        return Integer(Integer(other).numerator // self.numerator)


numbers.Integral.register(Integer)


class RationalConstant(Rational):
    """
    Abstract base class for rationals with specific behaviors

    Derived classes must define class attributes p and q and should probably all
    be singletons.

    """

    def __new__(cls):
        return AtomicExpr.__new__(cls)


class IntegerConstant(Integer):
    """Abstract class for integer constants."""

    def __new__(cls):
        return AtomicExpr.__new__(cls)


class Zero(IntegerConstant, metaclass=SingletonWithManagedProperties):
    """The number zero.

    Zero is a singleton, and can be accessed by ``S.Zero``

    Examples
    ========

    >>> Integer(0) is S.Zero
    True
    >>> 1/S.Zero
    zoo

    References
    ==========

    * https://en.wikipedia.org/wiki/Zero

    """

    _numerator = _int_dtype(0)
    _denominator = _int_dtype(1)

    is_positive = False
    is_negative = False
    is_zero = True
    is_number = True
    is_imaginary = True

    def _eval_power(self, expt):
        if expt.is_positive:
            return self
        if expt.is_negative:
            return zoo
        if expt.is_extended_real is False:
            return nan
        # infinities are already handled with pos and neg
        # tests above; now throw away leading numbers on Mul
        # exponent
        coeff, terms = expt.as_coeff_Mul()
        if coeff.is_negative:
            return zoo**terms
        if coeff is not S.One:  # there is a Number to discard
            return self**terms


class One(IntegerConstant, metaclass=SingletonWithManagedProperties):
    """The number one.

    One is a singleton, and can be accessed by ``S.One``.

    Examples
    ========

    >>> Integer(1) is S.One
    True

    References
    ==========

    * https://en.wikipedia.org/wiki/1_%28number%29

    """

    is_number = True

    _numerator = _int_dtype(1)
    _denominator = _int_dtype(1)


class NegativeOne(IntegerConstant, metaclass=SingletonWithManagedProperties):
    """The number negative one.

    NegativeOne is a singleton, and can be accessed by ``S.NegativeOne``.

    Examples
    ========

    >>> Integer(-1) is S.NegativeOne
    True

    See Also
    ========

    One

    References
    ==========

    * https://en.wikipedia.org/wiki/%E2%88%921_%28number%29

    """

    is_number = True

    _numerator = _int_dtype(-1)
    _denominator = _int_dtype(1)

    def _eval_power(self, expt):
        from .add import Add
        if isinstance(expt, Number):
            if isinstance(expt, Float):
                return Float(-1.0)**expt
            elif expt in (oo, -oo):
                return nan
            elif expt is S.Half:
                return I
            else:
                assert isinstance(expt, Rational)
                if expt.denominator == 2:
                    return I**Integer(expt.numerator)
                i, r = divmod(expt.numerator, expt.denominator)
                if i:
                    return self**i*self**Rational(r, expt.denominator)
        if isinstance(expt, Add):
            # Handle (-1)**((-1)**n/2 + m/2)
            e2 = 2*expt
            if e2.is_even and e2.could_extract_minus_sign():
                e2 *= self
            assert e2.is_Add
            i, p = e2.as_two_terms()
            if p.is_Pow and p.base is S.NegativeOne and p.exp.is_integer:
                i = (i + 1)/2
                if i.is_even:
                    return self**p.exp


class Half(RationalConstant, metaclass=SingletonWithManagedProperties):
    """The rational number 1/2.

    Half is a singleton, and can be accessed by ``S.Half``.

    Examples
    ========

    >>> Rational(1, 2) is S.Half
    True

    References
    ==========

    * https://en.wikipedia.org/wiki/One_half

    """

    is_number = True

    _numerator = _int_dtype(1)
    _denominator = _int_dtype(2)


class Infinity(Number, metaclass=SingletonWithManagedProperties):
    r"""Positive infinite quantity.

    In real analysis the symbol `\infty` denotes an unbounded
    limit: `x\to\infty` means that `x` grows without bound.

    Infinity is often used not only to define a limit but as a value
    in the affinely extended real number system.  Points labeled `+\infty`
    and `-\infty` can be added to the topological space of the real numbers,
    producing the two-point compactification of the real numbers.  Adding
    algebraic properties to this gives us the extended real numbers.

    Infinity is a singleton, and can be accessed by ``oo``,
    or can be imported as ``oo``.

    Examples
    ========

    >>> 1 + oo
    oo
    >>> 42/oo
    0
    >>> limit(exp(x), x, oo)
    oo

    See Also
    ========

    NegativeInfinity, NaN

    References
    ==========

    * https://en.wikipedia.org/wiki/Infinity

    """

    is_commutative = True
    is_positive = True
    is_infinite = True
    is_number = True
    is_prime = False

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Number):
            if other in (-oo, nan, Float('-inf')):
                return nan
            elif other.is_Float:
                return Float('inf')
            else:
                return oo
        return NotImplemented
    __radd__ = __add__

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Number):
            if other in (oo, nan, Float('inf')):
                return nan
            elif other.is_Float:
                return Float('inf')
            else:
                return oo
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        if isinstance(other, Number):
            return (-self) + other
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Number):
            if other is S.Zero or other is nan:
                return nan
            elif other.is_Float:
                if other == 0:
                    return nan
                if other > 0:
                    return Float('inf')
                else:
                    return Float('-inf')
            else:
                if other > 0:
                    return oo
                else:
                    return -oo
        return NotImplemented
    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        if isinstance(other, Number):
            if other in (oo, -oo, nan, Float('inf'), Float('-inf')):
                return nan
            elif other.is_Float:
                if other.is_nonnegative:
                    return Float('inf')
                else:
                    return Float('-inf')
            else:
                if other > 0:
                    return oo
                elif other < 0:
                    return -oo
                else:
                    return zoo
        return NotImplemented

    def __neg__(self):
        return S.NegativeInfinity

    def _eval_power(self, expt):
        """
        ``expt`` is symbolic object but not equal to 0 or 1.

        ================ ======= ==============================
        Expression       Result  Notes
        ================ ======= ==============================
        ``oo ** nan``    ``nan``
        ``oo ** -p``     ``0``   ``p`` is number, ``oo``
        ================ ======= ==============================

        See Also
        ========
        Pow
        NaN
        NegativeInfinity

        """
        from ..functions import re

        if expt.is_positive:
            return oo
        if expt.is_negative:
            return S.Zero
        if expt.is_real is False and expt.is_number:
            expt_real = re(expt)
            if expt_real.is_positive:
                return zoo
            elif expt_real.is_negative:
                return S.Zero
            elif expt_real.is_zero:
                return nan

    def _as_mpf_val(self, prec):
        return mlib.finf

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        return other is oo

    @_sympifyit('other', NotImplemented)
    def __lt__(self, other):
        if other.is_extended_real:
            return S.false
        return Expr.__lt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __le__(self, other):
        if other.is_extended_real:
            if other.is_finite or other == -oo:
                return S.false
            elif other.is_nonpositive:
                return S.false
            elif other is oo:
                return S.true
        return Expr.__le__(self, other)

    @_sympifyit('other', NotImplemented)
    def __gt__(self, other):
        if other.is_extended_real:
            if other.is_finite or other == -oo:
                return S.true
            elif other.is_nonpositive:
                return S.true
            elif other == oo:
                return S.false
        return Expr.__gt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __ge__(self, other):
        if other.is_extended_real:
            return S.true
        return Expr.__ge__(self, other)

    def __mod__(self, other):
        return nan

    __rmod__ = __mod__


oo: Infinity = S.Infinity


class NegativeInfinity(Number, metaclass=SingletonWithManagedProperties):
    """Negative infinite quantity.

    NegativeInfinity is a singleton, and can be accessed by ``-oo``.

    See Also
    ========

    Infinity

    """

    is_commutative = True
    is_negative = True
    is_infinite = True
    is_number = True

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Number):
            if other in (oo, nan, Float('inf')):
                return nan
            elif other.is_Float:
                return Float('-inf')
            else:
                return -oo
        return NotImplemented
    __radd__ = __add__

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Number):
            if other in (-oo, nan, Float('-inf')):
                return nan
            elif other.is_Float:
                return Float('-inf')
            else:
                return -oo
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __rsub__(self, other):
        if isinstance(other, Number):
            return (-self) + other
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Number):
            if other is S.Zero or other is nan:
                return nan
            elif other.is_Float:
                if other.is_zero:
                    return nan
                elif other.is_positive:
                    return Float('-inf')
                else:
                    return Float('inf')
            else:
                if other.is_positive:
                    return -oo
                else:
                    return oo
        return NotImplemented
    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        if isinstance(other, Number):
            if other in (oo, -oo, nan, Float('-inf'), Float('inf')):
                return nan
            elif other.is_Float:
                if other.is_nonnegative:
                    return Float('-inf')
                else:
                    return Float('inf')
            else:
                if other > 0:
                    return -oo
                elif other < 0:
                    return oo
                else:
                    return zoo
        return NotImplemented

    def __neg__(self):
        return oo

    def _eval_power(self, expt):
        """
        ``expt`` is symbolic object but not equal to 0 or 1.

        ================ ======= ==============================
        Expression       Result  Notes
        ================ ======= ==============================
        ``(-oo) ** nan`` ``nan``
        ``(-oo) ** oo``  ``nan``
        ``(-oo) ** -oo`` ``nan``
        ``(-oo) ** e``   ``oo``  ``e`` is positive even integer
        ``(-oo) ** o``   ``-oo`` ``o`` is positive odd integer
        ================ ======= ==============================

        See Also
        ========

        Infinity
        Pow
        NaN

        """
        if expt.is_number:
            if expt in (oo, -oo, nan):
                return nan

            return S.NegativeOne**expt*oo**expt

    def _as_mpf_val(self, prec):
        return mlib.fninf

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        return other is -oo

    @_sympifyit('other', NotImplemented)
    def __lt__(self, other):
        if other.is_extended_real:
            if other.is_finite or other is oo:
                return S.true
            elif other.is_nonnegative:
                return S.true
            elif other == -oo:
                return S.false
        return Expr.__lt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __le__(self, other):
        if other.is_extended_real:
            return S.true
        return Expr.__le__(self, other)

    @_sympifyit('other', NotImplemented)
    def __gt__(self, other):
        if other.is_extended_real:
            return S.false
        return Expr.__gt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __ge__(self, other):
        if other.is_extended_real:
            if other.is_finite or other is oo:
                return S.false
            elif other.is_nonnegative:
                return S.false
            elif other == -oo:
                return S.true
        return Expr.__ge__(self, other)

    def __mod__(self, other):
        return nan

    __rmod__ = __mod__


class NaN(Number, metaclass=SingletonWithManagedProperties):
    """
    Not a Number.

    This serves as a place holder for numeric values that are indeterminate.
    Most operations on NaN, produce another NaN.  Most indeterminate forms,
    such as ``0/0`` or ``oo - oo` produce NaN.  Two exceptions are ``0**0``
    and ``oo**0``, which all produce ``1`` (this is consistent with Python's
    float).

    NaN is loosely related to floating point nan, which is defined in the
    IEEE 754 floating point standard, and corresponds to the Python
    ``float('nan')``.  Differences are noted below.

    NaN is mathematically not equal to anything else, even NaN itself.  This
    explains the initially counter-intuitive results with ``Eq`` and ``==`` in
    the examples below.

    NaN is not comparable so inequalities raise a TypeError.  This is in
    constrast with floating point nan where all inequalities are false.

    NaN is a singleton, and can be accessed by ``nan``.

    Examples
    ========

    >>> nan is nan
    True
    >>> oo - oo
    nan
    >>> nan + 1
    nan
    >>> Eq(nan, nan)   # mathematical equality
    false
    >>> nan == nan     # structural equality
    True

    References
    ==========

    * https://en.wikipedia.org/wiki/NaN

    """

    is_commutative = True
    is_comparable = False
    is_finite = False
    is_number = True

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        return self

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        return self

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        return self

    @_sympifyit('other', NotImplemented)
    def __truediv__(self, other):
        return self

    def _as_mpf_val(self, prec):
        return mlib.fnan

    def __hash__(self):
        return super().__hash__()

    def __eq__(self, other):
        # NaN is structurally equal to another NaN
        return other is nan

    def _eval_Eq(self, other):
        # NaN is not mathematically equal to anything, even NaN
        return S.false

    # Expr will _sympify and raise TypeError
    __gt__ = Expr.__gt__
    __ge__ = Expr.__ge__
    __lt__ = Expr.__lt__
    __le__ = Expr.__le__


nan = S.NaN


class ComplexInfinity(AtomicExpr, metaclass=SingletonWithManagedProperties):
    r"""Complex infinity.

    In complex analysis the symbol `\tilde\infty`, called "complex
    infinity", represents a quantity with infinite magnitude, but
    undetermined complex phase.

    ComplexInfinity is a singleton, and can be accessed by as ``zoo``.

    Examples
    ========

    >>> zoo + 42
    zoo
    >>> 42/zoo
    0
    >>> zoo + zoo
    nan
    >>> zoo*zoo
    zoo

    See Also
    ========

    Infinity

    """

    is_commutative = True
    is_infinite = True
    is_number = True
    is_prime = False
    is_extended_real = False

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    def __neg__(self):
        return self

    def _eval_power(self, expt):
        if expt.is_positive:
            return zoo
        elif expt.is_negative:
            return S.Zero


zoo: ComplexInfinity = S.ComplexInfinity


class NumberSymbol(AtomicExpr):
    """Base class for symbolic numbers."""

    is_commutative = True
    is_finite = True
    is_number = True

    is_NumberSymbol = True

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    def approximation_interval(self, number_cls):
        """Return an interval with number_cls endpoints that contains the
        value of NumberSymbol.  If not implemented, then return None.

        """
        return  # pragma: no cover

    def _eval_evalf(self, prec):
        return Float._new(self._as_mpf_val(prec), prec)

    @_sympifyit('other', NotImplemented)
    def __eq__(self, other):
        if self is other:
            return True
        if isinstance(other, Number) and self.is_irrational:
            return False

        return False    # NumberSymbol != non-(Number|self)

    @_sympifyit('other', NotImplemented)
    def __lt__(self, other):
        if self is other:
            return S.false
        return Expr.__lt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __le__(self, other):
        if self is other:
            return S.true
        return Expr.__le__(self, other)

    @_sympifyit('other', NotImplemented)
    def __gt__(self, other):
        r = sympify((-self) < (-other), strict=True)
        if r in (S.true, S.false):
            return r
        else:
            return Expr.__gt__(self, other)

    @_sympifyit('other', NotImplemented)
    def __ge__(self, other):
        r = sympify((-self) <= (-other), strict=True)
        if r in (S.true, S.false):
            return r
        else:
            return Expr.__ge__(self, other)

    def __int__(self):
        raise NotImplementedError

    def __hash__(self):
        return super().__hash__()


class Exp1(NumberSymbol, metaclass=SingletonWithManagedProperties):
    r"""The `e` constant.

    The transcendental number `e = 2.718281828\ldots` is the base of the
    natural logarithm and of the exponential function, `e = \exp(1)`.
    Sometimes called Euler's number or Napier's constant.

    Exp1 is a singleton, and can be imported as ``E``.

    Examples
    ========

    >>> E is exp(1)
    True
    >>> log(E)
    1

    References
    ==========

    * https://en.wikipedia.org/wiki/E_%28mathematical_constant%29

    """

    is_real = True
    is_positive = True
    is_number = True
    is_transcendental = True

    def __abs__(self):
        return self

    def __int__(self):
        return 2

    def _as_mpf_val(self, prec):
        return mpmath.e(prec)._mpf_

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return Integer(2), Integer(3)

    def _eval_power(self, arg):
        from ..functions.elementary.exponential import log
        from . import Add, Mul, Pow
        if arg.is_Number:
            if arg is oo:
                return oo
            elif arg == -oo:
                return S.Zero
        elif isinstance(arg, log):
            return arg.args[0]
        elif arg.is_Mul:
            Ioo = I*oo
            if arg in [Ioo, -Ioo]:
                return nan

            coeff = arg.coeff(pi*I)
            if coeff:
                if (2*coeff).is_integer:
                    if coeff.is_even:
                        return S.One
                    elif coeff.is_odd:
                        return S.NegativeOne
                    elif (coeff + S.Half).is_even:
                        return -I
                    elif (coeff + S.Half).is_odd:
                        return I

            # Warning: code in risch.py will be very sensitive to changes
            # in this (see DifferentialExtension).

            # look for a single log factor

            coeff, terms = arg.as_coeff_Mul()

            # but it can't be multiplied by oo
            if coeff in (oo, -oo):
                return

            coeffs, log_term = [coeff], None
            for term in Mul.make_args(terms):
                if isinstance(term, log):
                    if log_term is None:
                        log_term = term.args[0]
                    else:
                        return
                elif term.is_comparable:
                    coeffs.append(term)
                else:
                    return

            return log_term**Mul(*coeffs) if log_term else None
        elif arg.is_Add:
            out = []
            add = []
            for a in arg.args:
                if a is S.One:
                    add.append(a)
                    continue
                newa = self**a
                if newa.is_Pow and newa.base is self:
                    add.append(a)
                else:
                    out.append(newa)
            if out:
                return Mul(*out)*Pow(self, Add(*add), evaluate=False)
        elif arg.is_Matrix:
            return arg.exp()

    def _eval_rewrite_as_sin(self):
        from ..functions import sin
        return sin(I + pi/2) - I*sin(I)

    def _eval_rewrite_as_cos(self):
        from ..functions import cos
        return cos(I) + I*cos(I + pi/2)


E = S.Exp1


class Pi(NumberSymbol, metaclass=SingletonWithManagedProperties):
    r"""The `\pi` constant.

    The transcendental number `\pi = 3.141592654\ldots` represents the ratio
    of a circle's circumference to its diameter, the area of the unit circle,
    the half-period of trigonometric functions, and many other things
    in mathematics.

    Pi is a singleton, and can be imported as ``pi``.

    Examples
    ========

    >>> pi > 3
    true
    >>> pi.is_irrational
    True
    >>> sin(x + 2*pi)
    sin(x)
    >>> integrate(exp(-x**2), (x, -oo, oo))
    sqrt(pi)

    References
    ==========

    * https://en.wikipedia.org/wiki/Pi

    """

    is_real = True
    is_positive = True
    is_number = True
    is_transcendental = True

    def __abs__(self):
        return self

    def __int__(self):
        return 3

    def _as_mpf_val(self, prec):
        return mpmath.pi(prec)._mpf_

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return Integer(3), Integer(4)
        elif issubclass(number_cls, Rational):
            return Rational(223, 71), Rational(22, 7)


pi: Pi = S.Pi


class GoldenRatio(NumberSymbol, metaclass=SingletonWithManagedProperties):
    r"""The golden ratio, `\phi`.

    `\phi = \frac{1 + \sqrt{5}}{2}` is algebraic number.  Two quantities
    are in the golden ratio if their ratio is the same as the ratio of
    their sum to the larger of the two quantities, i.e. their maximum.

    Examples
    ========

    >>> GoldenRatio > 1
    true
    >>> GoldenRatio.expand(func=True)
    1/2 + sqrt(5)/2
    >>> GoldenRatio.is_irrational
    True

    References
    ==========

    * https://en.wikipedia.org/wiki/Golden_ratio

    """

    is_real = True
    is_positive = True
    is_irrational = True
    is_number = True
    is_algebraic = True

    def __int__(self):
        return 1

    def _as_mpf_val(self, prec):
        return mpmath.phi(prec)._mpf_

    def _eval_expand_func(self, **hints):
        from ..functions import sqrt
        return S.Half + S.Half*sqrt(5)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return S.One, Integer(2)


class EulerGamma(NumberSymbol, metaclass=SingletonWithManagedProperties):
    r"""The Euler-Mascheroni constant.

    `\gamma = 0.5772157\ldots` (also called Euler's constant) is a mathematical
    constant recurring in analysis and number theory.  It is defined as the
    limiting difference between the harmonic series and the
    natural logarithm:

    .. math:: \gamma = \lim\limits_{n\to\infty}
              \left(\sum\limits_{k=1}^n\frac{1}{k} - \ln n\right)

    Examples
    ========

    >>> EulerGamma.is_irrational
    >>> EulerGamma > 0
    true
    >>> EulerGamma > 1
    false

    References
    ==========

    * https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant

    """

    is_real = True
    is_positive = True
    is_number = True

    def __int__(self):
        return 0

    def _as_mpf_val(self, prec):
        return mpmath.euler(prec)._mpf_

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return S.Zero, S.One
        elif issubclass(number_cls, Rational):
            return S.Half, Rational(3, 5)


class Catalan(NumberSymbol, metaclass=SingletonWithManagedProperties):
    r"""Catalan's constant.

    `K = 0.91596559\ldots` is given by the infinite series

    .. math:: K = \sum_{k=0}^{\infty} \frac{(-1)^k}{(2k+1)^2}

    Examples
    ========

    >>> Catalan.is_irrational
    >>> Catalan > 0
    true
    >>> Catalan > 1
    false

    References
    ==========

    * https://en.wikipedia.org/wiki/Catalan%27s_constant

    """

    is_real = True
    is_positive = True
    is_number = True

    def __int__(self):
        return 0

    def _as_mpf_val(self, prec):
        return mpmath.catalan(prec)._mpf_

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return S.Zero, S.One
        elif issubclass(number_cls, Rational):
            return Rational(9, 10), S.One


class ImaginaryUnit(AtomicExpr, metaclass=SingletonWithManagedProperties):
    r"""The imaginary unit, `i = \sqrt{-1}`.

    I is a singleton, and can be imported as ``I``.

    Examples
    ========

    >>> sqrt(-1)
    I
    >>> I*I
    -1
    >>> 1/I
    -I

    References
    ==========

    * https://en.wikipedia.org/wiki/Imaginary_unit

    """

    is_commutative = True
    is_imaginary = True
    is_extended_real = False
    is_finite = True
    is_number = True
    is_algebraic = True
    is_transcendental = False
    is_real = False

    def __abs__(self):
        return S.One

    def _eval_evalf(self, prec):
        return self

    def _eval_conjugate(self):
        return -I

    def _eval_power(self, expt):
        """Helper for Pow constructor.

        b is I = sqrt(-1)
        e is symbolic object but not equal to 0, 1

        I**r -> (-1)**(r/2) -> exp(r/2*Pi*I) -> sin(Pi*r/2) + cos(Pi*r/2)*I, r is decimal
        I**0 mod 4 -> 1
        I**1 mod 4 -> I
        I**2 mod 4 -> -1
        I**3 mod 4 -> -I

        """
        if isinstance(expt, Number):
            if isinstance(expt, Integer):
                expt = expt.numerator % 4
                if expt == 0:
                    return S.One
                if expt == 1:
                    return I
                if expt == 2:
                    return -S.One
                return -I
            return S.NegativeOne**(expt*S.Half)

    def as_base_exp(self):
        return S.NegativeOne, S.Half


I = S.ImaginaryUnit


# Add sympify converters

converter[float] = converter[decimal.Decimal] = Float

try:
    import numpy as np

    def _sympify_numpy(x):
        return Float(mpmath.mpmathify(x))

    converter[np.floating] = _sympify_numpy
except ImportError:
    pass

converter[int] = Integer
converter[fractions.Fraction] = Rational


if HAS_GMPY:
    converter[gmpy.mpz] = Integer
    converter[gmpy.mpq] = Rational


def _sympify_mpmath(x):
    return Expr._from_mpmath(x, x.context.prec)


converter[mpmath.mpf] = _sympify_mpmath
converter[mpmath.mpc] = _sympify_mpmath


def _sympify_complex(a):
    return sympify(a.real) + I*sympify(a.imag)


converter[complex] = _sympify_complex
