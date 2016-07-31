import decimal
import fractions
import math
import re as regex
from collections import defaultdict

import mpmath
import mpmath.libmp as mlib
from mpmath.libmp import mpf_pow, mpf_pi, mpf_e, phi_fixed
from mpmath.ctx_mp import mpnumeric
from mpmath.libmp.libmpf import (finf as _mpf_inf, fninf as _mpf_ninf,
                                 fnan as _mpf_nan, fzero as _mpf_zero,
                                 _normalize as mpf_normalize, prec_to_dps)

from .containers import Tuple
from .sympify import converter, sympify, _sympify, SympifyError
from .singleton import S, Singleton
from .expr import Expr, AtomicExpr
from .decorators import _sympifyit
from .cache import cacheit, clear_cache
from .logic import fuzzy_not
from .compatibility import as_int, HAS_GMPY, DIOFANT_INTS

rnd = mlib.round_nearest

_LOG2 = math.log(2)


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
            if type(z2) is str and getattr(z1, 'is_Number', False):
                return str(z1) == z2
            a, b = Float(z1), Float(z2)
            return int(abs(a - b)*10**prec_to_dps(
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
    not when the mantissa is 0. ``mpf_normlize`` always assumes that this
    is zero, but it may not be since the mantissa for mpf's values "+inf",
    "-inf" and "nan" have a mantissa of zero, too.

    Note: this is not intended to validate a given mpf tuple, so sending
    mpf tuples that were not created by mpmath may produce bad results. This
    is only a wrapper to ``mpf_normalize`` which provides the check for non-
    zero mpfs that have a 0 for the mantissa.
    """
    sign, man, expt, bc = mpf
    if not man:
        # hack for mpf_normalize which does not do this;
        # it assumes that if man is zero the result is 0
        # (see issue 6639)
        if not bc:
            return _mpf_zero
        else:
            # don't change anything; this should already
            # be a well formed mpf tuple
            return mpf
    rv = mpf_normalize(sign, man, expt, bc, prec, rnd)
    return rv

# TODO: we should use the warnings module
_errdict = {"divide": False}


def seterr(divide=False):
    """
    Should diofant raise an exception on 0/0 or return a nan?

    divide == True .... raise an exception
    divide == False ... return nan
    """
    if _errdict["divide"] != divide:
        clear_cache()
        _errdict["divide"] = divide


def _decimal_to_Rational_prec(dec):
    """Convert an ordinary decimal instance to a Rational."""
    if not dec.is_finite():
        raise TypeError("dec must be finite, got %s." % dec)
    s, d, e = dec.as_tuple()
    prec = len(d)
    if e >= 0:  # it's an integer
        rv = Integer(int(dec))
    else:
        s = (-1)**s
        d = sum(di*10**i for i, di in enumerate(reversed(d)))
        rv = Rational(s*d, 10**-e)
    return rv, prec


def _literal_float(f):
    """Return True if n can be interpreted as a floating point number."""
    pat = r"[-+]?((\d*\.\d+)|(\d+\.?))(eE[-+]?\d+)?"
    return bool(regex.match(pat, f))


@cacheit
def igcd(*args):
    """Computes positive integer greatest common divisor.

    Examples
    ========

    >>> from diofant.core.numbers import igcd
    >>> igcd(2, 4)
    2
    >>> igcd(5, 10, 15)
    5
    """
    a = args[0]
    for b in args[1:]:
        a, b = as_int(a), abs(as_int(b))

        a = fractions.gcd(a, b)
        if a == 1:
            break
    return a


def ilcm(*args):
    """Computes integer least common multiple.

    Examples
    ========

    >>> from diofant.core.numbers import ilcm
    >>> ilcm(5, 10)
    10
    >>> ilcm(7, 3)
    21
    >>> ilcm(5, 10, 15)
    30

    """
    if 0 in args:
        return 0
    a = args[0]
    for b in args[1:]:
        a = a*b // igcd(a, b)
    return a


def igcdex(a, b):
    """Returns x, y, g such that g = x*a + y*b = gcd(a, b).

       >>> from diofant.core.numbers import igcdex
       >>> igcdex(2, 3)
       (-1, 1, 1)
       >>> igcdex(10, 12)
       (-1, 1, 2)

       >>> x, y, g = igcdex(100, 2004)
       >>> x, y, g
       (-20, 1, 4)
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
        (c, q) = (a % b, a // b)
        (a, b, r, s, x, y) = (b, c, x - q*r, y - q*s, r, s)

    return x*x_sign, y*y_sign, a


class Number(AtomicExpr):
    """
    Represents any kind of number in diofant.

    Floating point numbers are represented by the Float class.
    Integer numbers (of any size), together with rational numbers (again,
    there is no limit on their size) are represented by the Rational class.

    If you want to represent, for example, ``1+sqrt(2)``, then you need to do::

      Rational(1) + sqrt(Rational(2))
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
        if isinstance(obj, DIOFANT_INTS):
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
                raise ValueError('String "%s" does not denote a Number' % obj)
        msg = "expected str|int|float|Decimal|Number object but got %r"
        raise TypeError(msg % type(obj).__name__)

    def __divmod__(self, other):
        from .containers import Tuple
        from diofant.functions.elementary.complexes import sign

        try:
            other = Number(other)
        except TypeError:
            msg = "unsupported operand type(s) for divmod(): '%s' and '%s'"
            raise TypeError(msg % (type(self).__name__, type(other).__name__))
        if not other:
            raise ZeroDivisionError('modulo by zero')
        if self.is_Integer and other.is_Integer:
            return Tuple(*divmod(self.p, other.p))
        else:
            rat = self/other
        w = sign(rat)*int(abs(rat))  # = rat.floor()
        r = self - other*w
        return Tuple(w, r)

    def __rdivmod__(self, other):
        try:
            other = Number(other)
        except TypeError:
            msg = "unsupported operand type(s) for divmod(): '%s' and '%s'"
            raise TypeError(msg % (type(other).__name__, type(self).__name__))
        return divmod(other, self)

    def __round__(self, *args):
        return round(float(self), *args)

    def _as_mpf_val(self, prec):
        """Evaluation of mpf tuple accurate to at least prec bits."""
        raise NotImplementedError('%s needs ._as_mpf_val() method' %
            (self.__class__.__name__))

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

    def _eval_is_finite(self):
        return True

    @classmethod
    def class_key(cls):
        return 1, 0, 'Number'

    @cacheit
    def sort_key(self, order=None):
        return self.class_key(), (0, ()), (), self

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if self is S.Zero:
            return other
        if isinstance(other, Number):
            if other is S.NaN:
                return S.NaN
            elif other is S.Infinity:
                return S.Infinity
            elif other is S.NegativeInfinity:
                return S.NegativeInfinity
        return AtomicExpr.__add__(self, other)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Number):
            if other is S.NaN:
                return S.NaN
            elif other is S.Infinity:
                return S.NegativeInfinity
            elif other is S.NegativeInfinity:
                return S.Infinity
        return AtomicExpr.__sub__(self, other)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if self is S.One:
            return other
        if isinstance(other, Number):
            if other is S.NaN:
                return S.NaN
            elif other is S.Infinity:
                if self.is_zero:
                    return S.NaN
                elif self.is_positive:
                    return S.Infinity
                else:
                    return S.NegativeInfinity
            elif other is S.NegativeInfinity:
                if self.is_zero:
                    return S.NaN
                elif self.is_positive:
                    return S.NegativeInfinity
                else:
                    return S.Infinity
        elif isinstance(other, Tuple):
            return NotImplemented
        return AtomicExpr.__mul__(self, other)

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        if isinstance(other, Number):
            if other is S.NaN:
                return S.NaN
            elif other is S.Infinity or other is S.NegativeInfinity:
                return S.Zero
        return AtomicExpr.__div__(self, other)

    __truediv__ = __div__

    def __eq__(self, other):
        raise NotImplementedError('%s needs .__eq__() method' %
            (self.__class__.__name__))

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s < %s" % (self, other))
        raise NotImplementedError('%s needs .__lt__() method' %
            (self.__class__.__name__))

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s <= %s" % (self, other))
        raise NotImplementedError('%s needs .__le__() method' %
            (self.__class__.__name__))

    def __gt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s > %s" % (self, other))
        return _sympify(other).__lt__(self)

    def __ge__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s >= %s" % (self, other))
        return _sympify(other).__le__(self)

    def __hash__(self):
        return super(Number, self).__hash__()

    def is_constant(self, *wrt, **flags):
        return True

    def as_coeff_mul(self, *deps, **kwargs):
        # a -> c*t
        if self.is_Rational or not kwargs.pop('rational', True):
            return self, tuple()
        elif self.is_negative:
            return S.NegativeOne, (-self,)
        return S.One, (self,)

    def as_coeff_add(self, *deps):
        # a -> c + t
        if self.is_Rational:
            return self, tuple()
        return S.Zero, (self,)

    def as_coeff_Mul(self, rational=False):
        """Efficiently extract the coefficient of a product. """
        if rational and not self.is_Rational:
            return S.One, self
        return self, S.One

    def as_coeff_Add(self):
        """Efficiently extract the coefficient of a summation. """
        return self, S.Zero

    def gcd(self, other):
        """Compute GCD of `self` and `other`. """
        from diofant.polys import gcd
        return gcd(self, other)

    def lcm(self, other):
        """Compute LCM of `self` and `other`. """
        from diofant.polys import lcm
        return lcm(self, other)

    def cofactors(self, other):
        """Compute GCD and cofactors of `self` and `other`. """
        from diofant.polys import cofactors
        return cofactors(self, other)


class Float(Number):
    """Represent a floating-point number of arbitrary precision.

    Examples
    ========

    >>> from diofant import Float
    >>> Float(3.5)
    3.50000000000000
    >>> Float(3)
    3.00000000000000

    Creating Floats from strings (and Python ``int`` type) will
    give a minimum precision of 15 digits, but the precision
    will automatically increase to capture all digits entered.

    >>> Float(1)
    1.00000000000000
    >>> Float(10**20)
    100000000000000000000.
    >>> Float('1e20')
    100000000000000000000.

    However, *floating-point* numbers (Python ``float`` types) retain
    only 15 digits of precision:

    >>> Float(1e20)
    1.00000000000000e+20
    >>> Float(1.23456789123456789)
    1.23456789123457

    It may be preferable to enter high-precision decimal numbers
    as strings:

    Float('1.23456789123456789')
    1.23456789123456789

    The desired number of digits can also be specified:

    >>> Float('1e-3', 3)
    0.00100
    >>> Float(100, 4)
    100.0

    Float can automatically count significant figures if a null string
    is sent for the precision; space are also allowed in the string. (Auto-
    counting is only allowed for strings and ints).

    >>> Float('123 456 789 . 123 456', '')
    123456789.123456
    >>> Float('12e-3', '')
    0.012
    >>> Float(3, '')
    3.

    If a number is written in scientific notation, only the digits before the
    exponent are considered significant if a decimal appears, otherwise the
    "e" signifies only how to move the decimal:

    >>> Float('60.e2', '')  # 2 digits significant
    6.0e+3
    >>> Float('60e2', '')  # 4 digits significant
    6000.
    >>> Float('600e-2', '')  # 3 digits significant
    6.00

    Notes
    =====

    Floats are inexact by their nature unless their value is a binary-exact
    value.

    >>> approx, exact = Float(.1, 1), Float(.125, 1)

    For calculation purposes, evalf needs to be able to change the precision
    but this will not increase the accuracy of the inexact value. The
    following is the most accurate 5-digit approximation of a value of 0.1
    that had only 1 digit of precision:

    >>> approx.evalf(5)
    0.099609

    By contrast, 0.125 is exact in binary (as it is in base 10) and so it
    can be passed to Float or evalf to obtain an arbitrary precision with
    matching accuracy:

    >>> Float(exact, 5)
    0.12500
    >>> exact.evalf(20)
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

    >>> def show(f): # binary rep of Float
    ...     from diofant import Mul, Pow
    ...     s, m, e, b = f._mpf_
    ...     v = Mul(int(m), Pow(2, int(e), evaluate=False), evaluate=False)
    ...     print('%s at prec=%s' % (v, f._prec))
    ...
    >>> t = Float('0.3', 3)
    >>> show(t)
    4915/2**14 at prec=13
    >>> show(Float(t, 20)) # higher prec, not higher accuracy
    4915/2**14 at prec=70
    >>> show(Float(t, 2)) # lower prec
    307/2**10 at prec=10

    The same thing happens when evalf is used on a Float:

    >>> show(t.evalf(20))
    4915/2**14 at prec=70
    >>> show(t.evalf(2))
    307/2**10 at prec=10

    Finally, Floats can be instantiated with an mpf tuple (n, c, p) to
    produce the number (-1)**n*c*2**p:

    >>> n, c, p = 1, 5, 0
    >>> (-1)**n*c*2**p
    -5
    >>> Float((1, 5, 0))
    -5.00000000000000

    An actual mpf tuple also contains the number of bits in c as the last
    element of the tuple:

    >>> _._mpf_
    (1, 5, 0, 3)

    This is not needed for instantiation and is not the same thing as the
    precision. The mpf tuple and the precision are two separate quantities
    that Float tracks.

    """

    # A Float represents many real numbers,
    # both rational and irrational.
    is_rational = None
    is_irrational = None
    is_number = True

    is_extended_real = True

    is_Float = True

    def __new__(cls, num, prec=None):
        if isinstance(num, str):
            num = num.replace(' ', '')
            if num.startswith('.') and len(num) > 1:
                num = '0' + num
            elif num.startswith('-.') and len(num) > 2:
                num = '-0.' + num[2:]
        elif isinstance(num, float) and num == 0:
            num = '0'
        elif isinstance(num, (DIOFANT_INTS, Integer)):
            num = str(num)  # faster than mlib.from_int
        elif isinstance(num, mpmath.mpf):
            num = num._mpf_

        if prec is None:
            dps = 15
            if isinstance(num, Float):
                return num
            if isinstance(num, str) and _literal_float(num):
                try:
                    Num = decimal.Decimal(num)
                except decimal.InvalidOperation:
                    pass
                else:
                    isint = '.' not in num
                    num, dps = _decimal_to_Rational_prec(Num)
                    if num.is_Integer and isint:
                        dps = max(dps, len(str(num).lstrip('-')))
                    dps = max(15, dps)
        elif prec == '':
            if not isinstance(num, str):
                raise ValueError('The null string can only be used when '
                'the number to Float is passed as a string or an integer.')
            ok = None
            if _literal_float(num):
                try:
                    Num = decimal.Decimal(num)
                except decimal.InvalidOperation:
                    pass
                else:
                    isint = '.' not in num
                    num, dps = _decimal_to_Rational_prec(Num)
                    if num.is_Integer and isint:
                        dps = max(dps, len(str(num).lstrip('-')))
                    ok = True
            if ok is None:
                raise ValueError('string-float not recognized: %s' % num)
        else:
            dps = prec

        prec = mlib.libmpf.dps_to_prec(dps)
        if isinstance(num, float):
            _mpf_ = mlib.from_float(num, prec, rnd)
        elif isinstance(num, str):
            _mpf_ = mlib.from_str(num, prec, rnd)
        elif isinstance(num, decimal.Decimal):
            if num.is_finite():
                _mpf_ = mlib.from_str(str(num), prec, rnd)
            elif num.is_nan():
                _mpf_ = _mpf_nan
            elif num.is_infinite():
                if num > 0:
                    _mpf_ = _mpf_inf
                else:
                    _mpf_ = _mpf_ninf
            else:
                raise ValueError("unexpected decimal value %s" % str(num))
        elif isinstance(num, Rational):
            _mpf_ = mlib.from_rational(num.p, num.q, prec, rnd)
        elif isinstance(num, tuple) and len(num) in (3, 4):
            if type(num[1]) is str:
                # it's a hexadecimal (coming from a pickled object)
                # assume that it is in standard form
                num = list(num)
                num[1] = int(num[1], 16)
                _mpf_ = tuple(num)
            else:
                if not num[1] and len(num) == 4:
                    # handle normalization hack
                    return Float._new(num, prec)
                else:
                    _mpf_ = mpmath.mpf(
                        S.NegativeOne**num[0]*num[1]*2**num[2])._mpf_
        elif isinstance(num, Float):
            _mpf_ = num._mpf_
            if prec < num._prec:
                _mpf_ = mpf_norm(_mpf_, prec)
        else:
            _mpf_ = mpmath.mpf(num)._mpf_

        # special cases
        if _mpf_ == _mpf_zero:
            pass  # we want a Float
        elif _mpf_ == _mpf_nan:
            return S.NaN

        obj = Expr.__new__(cls)
        obj._mpf_ = _mpf_
        obj._prec = prec
        return obj

    @classmethod
    def _new(cls, _mpf_, _prec):
        # special cases
        if _mpf_ == _mpf_zero:
            return S.Zero  # XXX this is different from Float which gives 0.0
        elif _mpf_ == _mpf_nan:
            return S.NaN

        obj = Expr.__new__(cls)
        obj._mpf_ = mpf_norm(_mpf_, _prec)
        obj._prec = _prec
        return obj

    # mpz can't be pickled
    def __getnewargs__(self):
        return mlib.to_pickable(self._mpf_),

    def __getstate__(self):
        return {'_prec': self._prec}

    def _hashable_content(self):
        return self._mpf_, self._prec

    def floor(self):
        return Integer(int(mlib.to_int(
            mlib.mpf_floor(self._mpf_, self._prec))))

    def ceiling(self):
        return Integer(int(mlib.to_int(
            mlib.mpf_ceil(self._mpf_, self._prec))))

    @property
    def num(self):
        return mpmath.mpf(self._mpf_)

    def _as_mpf_val(self, prec):
        from diofant.utilities.misc import debug
        rv = mpf_norm(self._mpf_, prec)
        if rv != self._mpf_ and self._prec == prec:
            debug(self._mpf_, rv)
        return rv

    def _as_mpf_op(self, prec):
        return self._mpf_, max(prec, self._prec)

    def _eval_is_finite(self):
        if self._mpf_ in (_mpf_inf, _mpf_ninf):
            return False
        return True

    def _eval_is_infinite(self):
        if self._mpf_ in (_mpf_inf, _mpf_ninf):
            return True
        return False

    def _eval_is_integer(self):
        return self._mpf_ == _mpf_zero

    def _eval_is_negative(self):
        if self._mpf_ == _mpf_ninf:
            return True
        if self._mpf_ == _mpf_inf:
            return False
        return self.num < 0

    def _eval_is_positive(self):
        if self._mpf_ == _mpf_inf:
            return True
        if self._mpf_ == _mpf_ninf:
            return False
        return self.num > 0

    def _eval_is_zero(self):
        return self._mpf_ == _mpf_zero

    def __bool__(self):
        return self._mpf_ != _mpf_zero

    def __neg__(self):
        return Float._new(mlib.mpf_neg(self._mpf_), self._prec)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_add(self._mpf_, rhs, prec, rnd), prec)
        return Number.__add__(self, other)

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

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        if isinstance(other, Number) and other != 0:
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_div(self._mpf_, rhs, prec, rnd), prec)
        return Number.__div__(self, other)

    __truediv__ = __div__

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if isinstance(other, Rational) and other.q != 1:
            # calculate mod with Rationals, *then* round the result
            return Float(Rational.__mod__(Rational(self), other),
                prec_to_dps(self._prec))
        if isinstance(other, Float):
            r = self/other
            if r == int(r):
                prec = max([prec_to_dps(i)
                    for i in (self._prec, other._prec)])
                return Float(0, prec)
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_mod(self._mpf_, rhs, prec, rnd), prec)
        return Number.__mod__(self, other)

    @_sympifyit('other', NotImplemented)
    def __rmod__(self, other):
        if isinstance(other, Float):
            return other.__mod__(self)
        if isinstance(other, Number):
            rhs, prec = other._as_mpf_op(self._prec)
            return Float._new(mlib.mpf_mod(rhs, self._mpf_, prec, rnd), prec)
        return Number.__rmod__(self, other)

    def _eval_power(self, expt):
        """
        expt is symbolic object but not equal to 0, 1

        (-p)**r -> exp(r*log(-p)) -> exp(r*(log(p) + I*Pi)) ->
                  -> p**r*(sin(Pi*r) + cos(Pi*r)*I)
        """
        if self == 0:
            if expt.is_positive:
                return S.Zero
            if expt.is_negative:
                return Float('inf')
        if isinstance(expt, Number):
            if isinstance(expt, Integer):
                prec = self._prec
                return Float._new(
                    mlib.mpf_pow_int(self._mpf_, expt.p, prec, rnd), prec)
            elif isinstance(expt, Rational) and \
                    expt.p == 1 and expt.q % 2 and self.is_negative:
                return Pow(S.NegativeOne, expt, evaluate=False)*(
                    -self)._eval_power(expt)
            expt, prec = expt._as_mpf_op(self._prec)
            mpfself = self._mpf_
            try:
                y = mpf_pow(mpfself, expt, prec, rnd)
                return Float._new(y, prec)
            except mlib.ComplexResult:
                re, im = mlib.mpc_pow(
                    (mpfself, _mpf_zero), (expt, _mpf_zero), prec, rnd)
                return Float._new(re, prec) + \
                    Float._new(im, prec)*S.ImaginaryUnit

    def __abs__(self):
        return Float._new(mlib.mpf_abs(self._mpf_), self._prec)

    def __int__(self):
        if self._mpf_ == _mpf_zero:
            return 0
        return int(mlib.to_int(self._mpf_))  # uses round_fast = round_down

    def __eq__(self, other):
        if isinstance(other, float):
            # coerce to Float at same precision
            o = Float(other)
            try:
                ompf = o._as_mpf_val(self._prec)
            except ValueError:
                return False
            return bool(mlib.mpf_eq(self._mpf_, ompf))
        try:
            other = _sympify(other)
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

    def __gt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s > %s" % (self, other))
        if isinstance(other, NumberSymbol):
            return other.__le__(self)
        if other.is_comparable:
            other = other.evalf()
        if isinstance(other, Number) and other is not S.NaN:
            return _sympify(bool(
                mlib.mpf_gt(self._mpf_, other._as_mpf_val(self._prec))))
        return Expr.__gt__(self, other)

    def __ge__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s >= %s" % (self, other))
        if isinstance(other, NumberSymbol):
            return other.__lt__(self)
        if other.is_comparable:
            other = other.evalf()
        if isinstance(other, Number) and other is not S.NaN:
            return _sympify(bool(
                mlib.mpf_ge(self._mpf_, other._as_mpf_val(self._prec))))
        return Expr.__ge__(self, other)

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s < %s" % (self, other))
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        if other.is_extended_real and other.is_number:
            other = other.evalf()
        if isinstance(other, Number) and other is not S.NaN:
            return _sympify(bool(
                mlib.mpf_lt(self._mpf_, other._as_mpf_val(self._prec))))
        return Expr.__lt__(self, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s <= %s" % (self, other))
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        if other.is_extended_real and other.is_number:
            other = other.evalf()
        if isinstance(other, Number) and other is not S.NaN:
            return _sympify(bool(
                mlib.mpf_le(self._mpf_, other._as_mpf_val(self._prec))))
        return Expr.__le__(self, other)

    def __hash__(self):
        return super(Float, self).__hash__()

    def epsilon_eq(self, other, epsilon="1e-15"):
        return abs(self - other) < Float(epsilon)

    def __format__(self, format_spec):
        return format(decimal.Decimal(str(self)), format_spec)


# Add sympify converters
converter[float] = converter[decimal.Decimal] = Float

# this is here to work nicely in Sage
RealNumber = Float


class Rational(Number):
    """Represents integers and rational numbers (p/q) of any size.

    Examples
    ========

    >>> from diofant import Rational, nsimplify, sympify, pi
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

    >>> Rational("1.23")
    123/100
    >>> Rational('1e-2')
    1/100
    >>> Rational(".1")
    1/10

    The conversion of other types of strings can be handled by
    the sympify() function, and conversion of floats to expressions
    or simple fractions can be handled with nsimplify:

    >>> sympify('3**2/10')  # general expressions
    9/10
    >>> nsimplify(.3)  # numbers that have a simple form
    3/10

    But if the input does not reduce to a literal Rational, an error will
    be raised:

    >>> Rational(pi)
    Traceback (most recent call last):
    ...
    TypeError: invalid input: pi


    Low-level access numerator and denominator as .p and .q:

    >>> r = Rational(3, 4)
    >>> r
    3/4
    >>> r.p
    3
    >>> r.q
    4

    Note that p and q return integers (not Diofant Integers) so some care
    is needed when using them in expressions:

    >>> r.p/r.q
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
                p, q = float(p).as_integer_ratio()

        if isinstance(p, Rational):
            p = fractions.Fraction(p.p, p.q)
        if isinstance(q, Rational):
            q = fractions.Fraction(q.p, q.q)

        try:
            f = fractions.Fraction(p)/fractions.Fraction(q)
            p, q = f.numerator, f.denominator
        except ValueError:
            raise TypeError('invalid input: %s, %s' % (p, q))
        except ZeroDivisionError:
            if p == 0:
                if _errdict["divide"]:
                    raise ValueError("Indeterminate 0/0")
                else:
                    return S.NaN
            else:
                return S.ComplexInfinity

        if q == 1:
            return Integer(p)
        if p == 1 and q == 2:
            return S.Half

        obj = Expr.__new__(cls)
        obj.p = p
        obj.q = q

        return obj

    def limit_denominator(self, max_denominator=1000000):
        """Closest Rational to self with denominator at most max_denominator.

        >>> from diofant import Rational
        >>> Rational('3.141592653589793').limit_denominator(10)
        22/7
        >>> Rational('3.141592653589793').limit_denominator(100)
        311/99

        """
        f = fractions.Fraction(self.p, self.q)
        return Rational(f.limit_denominator(fractions.Fraction(int(max_denominator))))

    def __getnewargs__(self):
        return self.p, self.q

    def _hashable_content(self):
        return self.p, self.q

    def _eval_is_positive(self):
        return self.p > 0

    def _eval_is_zero(self):
        return self.p == 0

    def __neg__(self):
        return Rational(-self.p, self.q)

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Rational):
            return Rational(self.p*other.q + self.q*other.p, self.q*other.q)
        elif isinstance(other, Float):
            return other + self
        else:
            return Number.__add__(self, other)

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Rational):
            return Rational(self.p*other.q - self.q*other.p, self.q*other.q)
        elif isinstance(other, Float):
            return -other + self
        else:
            return Number.__sub__(self, other)

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Rational):
            return Rational(self.p*other.p, self.q*other.q)
        elif isinstance(other, Float):
            return other*self
        else:
            return Number.__mul__(self, other)

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        if isinstance(other, Rational):
            if self.p and other.p == S.Zero:
                return S.ComplexInfinity
            else:
                return Rational(self.p*other.q, self.q*other.p)
        elif isinstance(other, Float):
            return self*(1/other)
        else:
            return Number.__div__(self, other)

    __truediv__ = __div__

    @_sympifyit('other', NotImplemented)
    def __mod__(self, other):
        if isinstance(other, Rational):
            n = (self.p*other.q) // (other.p*self.q)
            return Rational(self.p*other.q - n*other.p*self.q, self.q*other.q)
        if isinstance(other, Float):
            # calculate mod with Rationals, *then* round the answer
            return Float(self.__mod__(Rational(other)),
                prec_to_dps(other._prec))
        return Number.__mod__(self, other)

    @_sympifyit('other', NotImplemented)
    def __rmod__(self, other):
        if isinstance(other, Rational):
            return Rational.__mod__(other, self)
        return Number.__rmod__(self, other)

    def _eval_power(self, expt):
        if isinstance(expt, Number):
            if isinstance(expt, Float):
                return self._eval_evalf(expt._prec)**expt
            if expt.is_negative:
                # (3/4)**-2 -> (4/3)**2
                ne = -expt
                if (ne is S.One):
                    return Rational(self.q, self.p)
                if self.is_negative:
                    if expt.q != 1:
                        return -(S.NegativeOne)**((expt.p % expt.q) /
                               Integer(expt.q))*Rational(self.q, -self.p)**ne
                    else:
                        return S.NegativeOne**ne*Rational(self.q, -self.p)**ne
                else:
                    return Rational(self.q, self.p)**ne
            if expt is S.Infinity:  # -oo already caught by test for negative
                if self.p > self.q:
                    # (3/2)**oo -> oo
                    return S.Infinity
                if self.p < -self.q:
                    # (-3/2)**oo -> oo + I*oo
                    return S.Infinity + S.Infinity*S.ImaginaryUnit
                return S.Zero
            if isinstance(expt, Integer):
                # (4/3)**2 -> 4**2 / 3**2
                return Rational(self.p**expt.p, self.q**expt.p)
            if isinstance(expt, Rational):
                if self.p != 1:
                    # (4/3)**(5/6) -> 4**(5/6)*3**(-5/6)
                    return Integer(self.p)**expt*Integer(self.q)**(-expt)
                # as the above caught negative self.p, now self is positive
                return Integer(self.q)**Rational(
                expt.p*(expt.q - 1), expt.q) / \
                    Integer(self.q)**Integer(expt.p)

        if self.is_negative and expt.is_even:
            return (-self)**expt

        return

    def _as_mpf_val(self, prec):
        return mlib.from_rational(self.p, self.q, prec, rnd)

    def _mpmath_(self, prec, rnd):
        return mpmath.make_mpf(mlib.from_rational(self.p, self.q, prec, rnd))

    def __abs__(self):
        return Rational(abs(self.p), self.q)

    def __int__(self):
        p, q = self.p, self.q
        if p < 0:
            return -(-p//q)
        return p//q

    def __eq__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # diofant != other  -->  not ==
        if isinstance(other, NumberSymbol):
            if other.is_irrational:
                return False
            return other.__eq__(self)
        if isinstance(other, Number):
            if isinstance(other, Rational):
                # a Rational is always in reduced form so will never be 2/4
                # so we can just check equivalence of args
                return self.p == other.p and self.q == other.q
            if isinstance(other, Float):
                return mlib.mpf_eq(self._as_mpf_val(other._prec), other._mpf_)
        return False

    def __gt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s > %s" % (self, other))
        if isinstance(other, NumberSymbol):
            return other.__le__(self)
        expr = self
        if isinstance(other, Number):
            if isinstance(other, Rational):
                return _sympify(bool(self.p*other.q > self.q*other.p))
            if isinstance(other, Float):
                return _sympify(bool(mlib.mpf_gt(
                    self._as_mpf_val(other._prec), other._mpf_)))
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.p), self.q*other
        return Expr.__gt__(expr, other)

    def __ge__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s >= %s" % (self, other))
        if isinstance(other, NumberSymbol):
            return other.__lt__(self)
        expr = self
        if isinstance(other, Number):
            if isinstance(other, Rational):
                return _sympify(bool(self.p*other.q >= self.q*other.p))
            if isinstance(other, Float):
                return _sympify(bool(mlib.mpf_ge(
                    self._as_mpf_val(other._prec), other._mpf_)))
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.p), self.q*other
        return Expr.__ge__(expr, other)

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s < %s" % (self, other))
        if isinstance(other, NumberSymbol):
            return other.__ge__(self)
        expr = self
        if isinstance(other, Number):
            if isinstance(other, Rational):
                return _sympify(bool(self.p*other.q < self.q*other.p))
            if isinstance(other, Float):
                return _sympify(bool(mlib.mpf_lt(
                    self._as_mpf_val(other._prec), other._mpf_)))
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.p), self.q*other
        return Expr.__lt__(expr, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s <= %s" % (self, other))
        expr = self
        if isinstance(other, NumberSymbol):
            return other.__gt__(self)
        elif isinstance(other, Number):
            if isinstance(other, Rational):
                return _sympify(bool(self.p*other.q <= self.q*other.p))
            if isinstance(other, Float):
                return _sympify(bool(mlib.mpf_le(
                    self._as_mpf_val(other._prec), other._mpf_)))
        elif other.is_number and other.is_extended_real:
            expr, other = Integer(self.p), self.q*other
        return Expr.__le__(expr, other)

    def __hash__(self):
        return super(Rational, self).__hash__()

    def factors(self, limit=None, use_trial=True, use_rho=False,
                use_pm1=False, verbose=False, visual=False):
        """A wrapper to factorint which return factors of self that are
        smaller than limit (or cheap to compute). Special methods of
        factoring are disabled by default so that only trial division is used.
        """
        from diofant.ntheory import factorint

        f = factorint(self.p, limit=limit, use_trial=use_trial,
                      use_rho=use_rho, use_pm1=use_pm1,
                      verbose=verbose).copy()
        f = defaultdict(int, f)
        for p, e in factorint(self.q, limit=limit,
                              use_trial=use_trial,
                              use_rho=use_rho,
                              use_pm1=use_pm1,
                              verbose=verbose).items():
            f[p] += -e

        if len(f) > 1 and 1 in f:
            del f[1]
        if not f:
            f = {1: 1}
        if not visual:
            return dict(f)
        else:
            if -1 in f:
                f.pop(-1)
                args = [S.NegativeOne]
            else:
                args = []
            args.extend([Pow(*i, evaluate=False)
                         for i in sorted(f.items())])
            return Mul(*args, evaluate=False)

    @_sympifyit('other', NotImplemented)
    def gcd(self, other):
        if isinstance(other, Rational):
            if other is S.Zero:
                return other
            return Rational(
                Integer(igcd(self.p, other.p)),
                Integer(ilcm(self.q, other.q)))
        return Number.gcd(self, other)

    @_sympifyit('other', NotImplemented)
    def lcm(self, other):
        if isinstance(other, Rational):
            return Rational(
                self.p*other.p//igcd(self.p, other.p),
                igcd(self.q, other.q))
        return Number.lcm(self, other)

    def as_numer_denom(self):
        return Integer(self.p), Integer(self.q)

    def as_content_primitive(self, radical=False):
        """Return the tuple (R, self/R) where R is the positive Rational
        extracted from self.

        Examples
        ========

        >>> from diofant import Rational
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


class Integer(Rational):

    q = 1
    is_integer = True
    is_number = True

    is_Integer = True

    def _as_mpf_val(self, prec):
        return mlib.from_int(self.p, prec)

    def _mpmath_(self, prec, rnd):
        return mpmath.make_mpf(self._as_mpf_val(prec))

    @cacheit
    def __new__(cls, i):
        if isinstance(i, str):
            i = i.replace(' ', '')
        # whereas we cannot, in general, make a Rational from an
        # arbitrary expression, we can make an Integer unambiguously
        # (except when a non-integer expression happens to round to
        # an integer). So we proceed by taking int() of the input and
        # let the int routines determine whether the expression can
        # be made into an int or whether an error should be raised.
        try:
            ival = int(i)
        except TypeError:
            raise TypeError(
                'Integer can only work with integer expressions.')

        # We only work with well-behaved integer types. This converts, for
        # example, numpy.int32 instances.
        if ival == 0:
            obj = S.Zero
        elif ival == 1:
            obj = S.One
        elif ival == -1:
            obj = S.NegativeOne
        else:
            obj = Expr.__new__(cls)
            obj.p = ival

        return obj

    def __getnewargs__(self):
        return self.p,

    # Arithmetic operations are here for efficiency
    def __int__(self):
        return self.p

    def __neg__(self):
        return Integer(-self.p)

    def __abs__(self):
        if self.p >= 0:
            return self
        else:
            return Integer(-self.p)

    def __divmod__(self, other):
        from .containers import Tuple
        if isinstance(other, Integer):
            return Tuple(*(divmod(self.p, other.p)))
        else:
            return Number.__divmod__(self, other)

    def __rdivmod__(self, other):
        from .containers import Tuple
        if isinstance(other, int):
            return Tuple(*(divmod(other, self.p)))
        else:
            try:
                other = Number(other)
            except TypeError:
                msg = "unsupported operand type(s) for divmod(): '%s' and '%s'"
                oname = type(other).__name__
                sname = type(self).__name__
                raise TypeError(msg % (oname, sname))
            return Number.__divmod__(other, self)

    # TODO make it decorator + bytecodehacks?
    def __add__(self, other):
        if isinstance(other, int):
            return Integer(self.p + other)
        elif isinstance(other, Integer):
            return Integer(self.p + other.p)
        return Rational.__add__(self, other)

    def __radd__(self, other):
        if isinstance(other, int):
            return Integer(other + self.p)
        return Rational.__add__(self, other)

    def __sub__(self, other):
        if isinstance(other, int):
            return Integer(self.p - other)
        elif isinstance(other, Integer):
            return Integer(self.p - other.p)
        return Rational.__sub__(self, other)

    def __rsub__(self, other):
        if isinstance(other, int):
            return Integer(other - self.p)
        return Rational.__rsub__(self, other)

    def __mul__(self, other):
        if isinstance(other, int):
            return Integer(self.p*other)
        elif isinstance(other, Integer):
            return Integer(self.p*other.p)
        return Rational.__mul__(self, other)

    def __rmul__(self, other):
        if isinstance(other, int):
            return Integer(other*self.p)
        return Rational.__mul__(self, other)

    def __mod__(self, other):
        if isinstance(other, int):
            return Integer(self.p % other)
        elif isinstance(other, Integer):
            return Integer(self.p % other.p)
        return Rational.__mod__(self, other)

    def __rmod__(self, other):
        if isinstance(other, int):
            return Integer(other % self.p)
        elif isinstance(other, Integer):
            return Integer(other.p % self.p)
        return Rational.__rmod__(self, other)

    def __eq__(self, other):
        if isinstance(other, int):
            return (self.p == other)
        elif isinstance(other, Integer):
            return (self.p == other.p)
        return Rational.__eq__(self, other)

    def __gt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s > %s" % (self, other))
        if isinstance(other, Integer):
            return _sympify(self.p > other.p)
        return Rational.__gt__(self, other)

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s < %s" % (self, other))
        if isinstance(other, Integer):
            return _sympify(self.p < other.p)
        return Rational.__lt__(self, other)

    def __ge__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s >= %s" % (self, other))
        if isinstance(other, Integer):
            return _sympify(self.p >= other.p)
        return Rational.__ge__(self, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s <= %s" % (self, other))
        if isinstance(other, Integer):
            return _sympify(self.p <= other.p)
        return Rational.__le__(self, other)

    def __hash__(self):
        return hash(self.p)

    def __index__(self):
        return self.p

    ########################################

    def _eval_is_odd(self):
        return bool(self.p % 2)

    def _eval_power(self, expt):
        """
        Tries to do some simplifications on self**expt

        Returns None if no further simplifications can be done

        When exponent is a fraction (so we have for example a square root),
        we try to find a simpler representation by factoring the argument
        up to factors of 2**15, e.g.

          - sqrt(4) becomes 2
          - sqrt(-4) becomes 2*I
          - (2**(3+7)*3**(6+7))**Rational(1,7) becomes 6*18**(3/7)

        Further simplification would require a special call to factorint on
        the argument which is not done here for sake of speed.

        """
        from diofant import perfect_power

        if expt is S.Infinity:
            if self.p > S.One:
                return S.Infinity
            # cases -1, 0, 1 are done in their respective classes
            return S.Infinity + S.ImaginaryUnit*S.Infinity
        if expt is S.NegativeInfinity:
            return Rational(1, self)**S.Infinity
        if not isinstance(expt, Number):
            # simplify when expt is even
            # (-2)**k --> 2**k
            if self.is_negative and expt.is_even:
                return (-self)**expt
        if isinstance(expt, Float):
            # Rational knows how to exponentiate by a Float
            return super(Integer, self)._eval_power(expt)
        if not isinstance(expt, Rational):
            return
        if expt is S.Half and self.is_negative:
            # we extract I for this special case since everyone is doing so
            return S.ImaginaryUnit*Pow(-self, expt)
        if expt.is_negative:
            # invert base and change sign on exponent
            ne = -expt
            if self.is_negative:
                if expt.q != 1:
                    return -(S.NegativeOne)**((expt.p % expt.q) /
                            Integer(expt.q))*Rational(1, -self)**ne
                else:
                    return (S.NegativeOne)**ne*Rational(1, -self)**ne
            else:
                return Rational(1, self.p)**ne
        # see if base is a perfect root, sqrt(4) --> 2
        x, xexact = integer_nthroot(abs(self.p), expt.q)
        if xexact:
            # if it's a perfect root we've finished
            result = Integer(x**abs(expt.p))
            if self.is_negative:
                result *= S.NegativeOne**expt
            return result

        # The following is an algorithm where we collect perfect roots
        # from the factors of base.

        # if it's not an nth root, it still might be a perfect power
        b_pos = int(abs(self.p))
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
            exponent *= expt.p
            # remove multiples of expt.q: (2**12)**(1/10) -> 2*(2**2)**(1/10)
            div_e, div_m = divmod(exponent, expt.q)
            if div_e > 0:
                out_int *= prime**div_e
            if div_m > 0:
                # see if the reduced exponent shares a gcd with e.q
                # (2**2)**(1/10) -> 2**(1/5)
                g = igcd(div_m, expt.q)
                if g != 1:
                    out_rad *= Pow(prime, Rational(div_m//g, expt.q//g))
                else:
                    sqr_dict[prime] = div_m
        # identify gcd of remaining powers
        for p, ex in sqr_dict.items():
            if sqr_gcd == 0:
                sqr_gcd = ex
            else:
                sqr_gcd = igcd(sqr_gcd, ex)
                if sqr_gcd == 1:
                    break
        for k, v in sqr_dict.items():
            sqr_int *= k**(v//sqr_gcd)
        if sqr_int == self and out_int == 1 and out_rad == 1:
            result = None
        else:
            result = out_int*out_rad*Pow(sqr_int, Rational(sqr_gcd, expt.q))
        return result

    def _eval_is_prime(self):
        from diofant.ntheory import isprime

        return isprime(self)

    def _eval_is_composite(self):
        if self > 1:
            return fuzzy_not(self.is_prime)
        else:
            return False

    def as_numer_denom(self):
        return self, S.One

    def __floordiv__(self, other):
        return Integer(self.p // Integer(other).p)

    def __rfloordiv__(self, other):
        return Integer(Integer(other).p // self.p)

# Add sympify converters
converter[int] = Integer


class AlgebraicNumber(Expr):
    """Class for representing algebraic numbers in Diofant. """

    is_AlgebraicNumber = True
    is_algebraic = True
    is_number = True

    def __new__(cls, expr, coeffs=Tuple(), alias=None, **args):
        """Construct a new algebraic number. """
        from diofant import Poly
        from diofant.polys.polyclasses import ANP, DMP
        from diofant.polys.numberfields import minimal_polynomial
        from diofant.core.symbol import Symbol

        expr = sympify(expr)

        if isinstance(expr, (tuple, Tuple)):
            minpoly, root = expr

            if not minpoly.is_Poly:
                minpoly = Poly(minpoly)
        elif expr.is_AlgebraicNumber:
            minpoly, root = expr.minpoly, expr.root
        else:
            minpoly, root = minimal_polynomial(
                expr, args.get('gen'), polys=True), expr

        dom = minpoly.get_domain()

        if coeffs != Tuple():
            if not isinstance(coeffs, ANP):
                rep = DMP.from_diofant_list(sympify(coeffs), 0, dom)
                scoeffs = Tuple(*coeffs)
            else:
                rep = DMP.from_list(coeffs.to_list(), 0, dom)
                scoeffs = Tuple(*coeffs.to_list())

            if rep.degree() >= minpoly.degree():
                rep = rep.rem(minpoly.rep)

            sargs = (root, scoeffs)

        else:
            rep = DMP.from_list([1, 0], 0, dom)

            if root.is_negative:
                rep = -rep

            sargs = (root, coeffs)

        if alias is not None:
            if not isinstance(alias, Symbol):
                alias = Symbol(alias)
            sargs = sargs + (alias,)

        obj = Expr.__new__(cls, *sargs)

        obj.rep = rep
        obj.root = root
        obj.alias = alias
        obj.minpoly = minpoly

        return obj

    def __hash__(self):
        return super(AlgebraicNumber, self).__hash__()

    def _eval_evalf(self, prec):
        return self.as_expr()._evalf(prec)

    @property
    def is_aliased(self):
        """Returns ``True`` if ``alias`` was set. """
        return self.alias is not None

    def as_poly(self, x=None):
        """Create a Poly instance from ``self``. """
        from diofant import Dummy, Poly, PurePoly
        if x is not None:
            return Poly.new(self.rep, x)
        else:
            if self.alias is not None:
                return Poly.new(self.rep, self.alias)
            else:
                return PurePoly.new(self.rep, Dummy('x'))

    def as_expr(self, x=None):
        """Create a Basic expression from ``self``. """
        return self.as_poly(x or self.root).as_expr().expand()

    def coeffs(self):
        """Returns all Diofant coefficients of an algebraic number. """
        return [ self.rep.domain.to_diofant(c) for c in self.rep.all_coeffs() ]

    def native_coeffs(self):
        """Returns all native coefficients of an algebraic number. """
        return self.rep.all_coeffs()

    def to_algebraic_integer(self):
        """Convert ``self`` to an algebraic integer. """
        from diofant import Poly
        f = self.minpoly

        if f.LC() == 1:
            return self

        coeff = f.LC()**(f.degree() - 1)
        poly = f.compose(Poly(f.gen/f.LC()))

        minpoly = poly*coeff
        root = f.LC()*self.root

        return AlgebraicNumber((minpoly, root), self.coeffs())

    def _eval_simplify(self, ratio, measure):
        from diofant.polys import RootOf, minpoly

        for r in [r for r in self.minpoly.all_roots() if r.func != RootOf]:
            if minpoly(self.root - r).is_Symbol:
                # use the matching root if it's simpler
                if measure(r) < ratio*measure(self.root):
                    return AlgebraicNumber(r)
        return self


class RationalConstant(Rational):
    """
    Abstract base class for rationals with specific behaviors

    Derived classes must define class attributes p and q and should probably all
    be singletons.
    """

    def __new__(cls):
        return AtomicExpr.__new__(cls)


class IntegerConstant(Integer):

    def __new__(cls):
        return AtomicExpr.__new__(cls)


class Zero(IntegerConstant, metaclass=Singleton):
    """The number zero.

    Zero is a singleton, and can be accessed by ``S.Zero``

    Examples
    ========

    >>> from diofant import S, Integer, zoo
    >>> Integer(0) is S.Zero
    True
    >>> 1/S.Zero
    zoo

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Zero
    """

    p = 0
    q = 1
    is_positive = False
    is_negative = False
    is_zero = True
    is_number = True
    is_imaginary = True

    @staticmethod
    def __abs__():
        return S.Zero

    @staticmethod
    def __neg__():
        return S.Zero

    def _eval_power(self, expt):
        if expt.is_positive:
            return self
        if expt.is_negative:
            return S.ComplexInfinity
        if expt.is_extended_real is False:
            return S.NaN
        # infinities are already handled with pos and neg
        # tests above; now throw away leading numbers on Mul
        # exponent
        coeff, terms = expt.as_coeff_Mul()
        if coeff.is_negative:
            return S.ComplexInfinity**terms
        if coeff is not S.One:  # there is a Number to discard
            return self**terms

    def __bool__(self):
        return False


class One(IntegerConstant, metaclass=Singleton):
    """The number one.

    One is a singleton, and can be accessed by ``S.One``.

    Examples
    ========

    >>> from diofant import S, Integer
    >>> Integer(1) is S.One
    True

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/1_%28number%29
    """
    is_number = True

    p = 1
    q = 1

    @staticmethod
    def __abs__():
        return S.One

    @staticmethod
    def __neg__():
        return S.NegativeOne

    def _eval_power(self, expt):
        return self

    @staticmethod
    def factors(limit=None, use_trial=True, use_rho=False, use_pm1=False,
                verbose=False, visual=False):
        if visual:
            return S.One
        return {1: 1}


class NegativeOne(IntegerConstant, metaclass=Singleton):
    """The number negative one.

    NegativeOne is a singleton, and can be accessed by ``S.NegativeOne``.

    Examples
    ========

    >>> from diofant import S, Integer
    >>> Integer(-1) is S.NegativeOne
    True

    See Also
    ========

    One

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/%E2%88%921_%28number%29

    """
    is_number = True

    p = -1
    q = 1

    @staticmethod
    def __abs__():
        return S.One

    @staticmethod
    def __neg__():
        return S.One

    def _eval_power(self, expt):
        if expt.is_odd:
            return S.NegativeOne
        if expt.is_even:
            return S.One
        if isinstance(expt, Number):
            if isinstance(expt, Float):
                return Float(-1.0)**expt
            if expt is S.NaN:
                return S.NaN
            if expt is S.Infinity or expt is S.NegativeInfinity:
                return S.NaN
            if expt is S.Half:
                return S.ImaginaryUnit
            if isinstance(expt, Rational):
                if expt.q == 2:
                    return S.ImaginaryUnit**Integer(expt.p)
                i, r = divmod(expt.p, expt.q)
                if i:
                    return self**i*self**Rational(r, expt.q)
        if isinstance(expt, Add):
            # Handle (-1)**((-1)**n/2 + m/2)
            e2 = 2*expt
            if e2.is_even:
                if e2.could_extract_minus_sign():
                    e2 *= self
            if e2.is_Add:
                i, p = e2.as_two_terms()
                if p.is_Pow and p.base is S.NegativeOne:
                    if p.exp.is_integer:
                        i = (i + 1)/2
                        if i.is_even:
                            return self**p.exp
                        elif i.is_odd:
                            return self**(p.exp + 1)
                        else:
                            return self**(p.exp + i)


class Half(RationalConstant, metaclass=Singleton):
    """The rational number 1/2.

    Half is a singleton, and can be accessed by ``S.Half``.

    Examples
    ========

    >>> from diofant import S, Rational
    >>> Rational(1, 2) is S.Half
    True

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/One_half
    """
    is_number = True

    p = 1
    q = 2

    @staticmethod
    def __abs__():
        return S.Half


class Infinity(Number, metaclass=Singleton):
    r"""Positive infinite quantity.

    In real analysis the symbol `\infty` denotes an unbounded
    limit: `x\to\infty` means that `x` grows without bound.

    Infinity is often used not only to define a limit but as a value
    in the affinely extended real number system.  Points labeled `+\infty`
    and `-\infty` can be added to the topological space of the real numbers,
    producing the two-point compactification of the real numbers.  Adding
    algebraic properties to this gives us the extended real numbers.

    Infinity is a singleton, and can be accessed by ``S.Infinity``,
    or can be imported as ``oo``.

    Examples
    ========

    >>> from diofant import oo, exp, limit, Symbol
    >>> 1 + oo
    oo
    >>> 42/oo
    0
    >>> x = Symbol('x')
    >>> limit(exp(x), x, oo)
    oo

    See Also
    ========

    NegativeInfinity, NaN

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Infinity
    """

    is_commutative = True
    is_positive = True
    is_infinite = True
    is_number = True
    is_prime = False

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    def _latex(self, printer):
        return r"\infty"

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Number):
            if other is S.NegativeInfinity or other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other == Float('-inf'):
                    return S.NaN
                else:
                    return Float('inf')
            else:
                return S.Infinity
        return NotImplemented
    __radd__ = __add__

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Number):
            if other is S.Infinity or other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other == Float('inf'):
                    return S.NaN
                else:
                    return Float('inf')
            else:
                return S.Infinity
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Number):
            if other is S.Zero or other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other == 0:
                    return S.NaN
                if other > 0:
                    return Float('inf')
                else:
                    return Float('-inf')
            else:
                if other > 0:
                    return S.Infinity
                else:
                    return S.NegativeInfinity
        return NotImplemented
    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        if isinstance(other, Number):
            if other is S.Infinity or \
                other is S.NegativeInfinity or \
                    other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other == Float('-inf') or \
                        other == Float('inf'):
                    return S.NaN
                elif other.is_nonnegative:
                    return Float('inf')
                else:
                    return Float('-inf')
            else:
                if other >= 0:
                    return S.Infinity
                else:
                    return S.NegativeInfinity
        return NotImplemented

    __truediv__ = __div__

    def __abs__(self):
        return S.Infinity

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
        from diofant.functions import re

        if expt.is_positive:
            return S.Infinity
        if expt.is_negative:
            return S.Zero
        if expt is S.NaN:
            return S.NaN
        elif expt is S.ComplexInfinity:
            return S.NaN
        if expt.is_real is False and expt.is_number:
            expt_real = re(expt)
            if expt_real.is_positive:
                return S.ComplexInfinity
            if expt_real.is_negative:
                return S.Zero
            if expt_real.is_zero:
                return S.NaN

            return self**expt.evalf()

    def _as_mpf_val(self, prec):
        return mlib.finf

    def __hash__(self):
        return super(Infinity, self).__hash__()

    def __eq__(self, other):
        return other is S.Infinity

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s < %s" % (self, other))
        if other.is_extended_real:
            return S.false
        return Expr.__lt__(self, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s <= %s" % (self, other))
        if other.is_extended_real:
            if other.is_finite or other is S.NegativeInfinity:
                return S.false
            elif other.is_nonpositive:
                return S.false
            elif other is S.Infinity:
                return S.true
        return Expr.__le__(self, other)

    def __gt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s > %s" % (self, other))
        if other.is_extended_real:
            if other.is_finite or other is S.NegativeInfinity:
                return S.true
            elif other.is_nonpositive:
                return S.true
            elif other is S.Infinity:
                return S.false
        return Expr.__gt__(self, other)

    def __ge__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s >= %s" % (self, other))
        if other.is_extended_real:
            return S.true
        return Expr.__ge__(self, other)

    def __mod__(self, other):
        return S.NaN

    __rmod__ = __mod__

oo = S.Infinity


class NegativeInfinity(Number, metaclass=Singleton):
    """Negative infinite quantity.

    NegativeInfinity is a singleton, and can be accessed
    by ``S.NegativeInfinity``.

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

    def _latex(self, printer):
        return r"-\infty"

    @_sympifyit('other', NotImplemented)
    def __add__(self, other):
        if isinstance(other, Number):
            if other is S.Infinity or other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other == Float('inf'):
                    return Float('nan')
                else:
                    return Float('-inf')
            else:
                return S.NegativeInfinity
        return NotImplemented
    __radd__ = __add__

    @_sympifyit('other', NotImplemented)
    def __sub__(self, other):
        if isinstance(other, Number):
            if other is S.NegativeInfinity or other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other == Float('-inf'):
                    return Float('nan')
                else:
                    return Float('-inf')
            else:
                return S.NegativeInfinity
        return NotImplemented

    @_sympifyit('other', NotImplemented)
    def __mul__(self, other):
        if isinstance(other, Number):
            if other is S.Zero or other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other is S.NaN or other.is_zero:
                    return S.NaN
                elif other.is_positive:
                    return Float('-inf')
                else:
                    return Float('inf')
            else:
                if other.is_positive:
                    return S.NegativeInfinity
                else:
                    return S.Infinity
        return NotImplemented
    __rmul__ = __mul__

    @_sympifyit('other', NotImplemented)
    def __div__(self, other):
        if isinstance(other, Number):
            if other is S.Infinity or \
                other is S.NegativeInfinity or \
                    other is S.NaN:
                return S.NaN
            elif other.is_Float:
                if other == Float('-inf') or \
                    other == Float('inf') or \
                        other is S.NaN:
                    return S.NaN
                elif other.is_nonnegative:
                    return Float('-inf')
                else:
                    return Float('inf')
            else:
                if other >= 0:
                    return S.NegativeInfinity
                else:
                    return S.Infinity
        return NotImplemented

    __truediv__ = __div__

    def __abs__(self):
        return S.Infinity

    def __neg__(self):
        return S.Infinity

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
            if expt is S.NaN or \
                expt is S.Infinity or \
                    expt is S.NegativeInfinity:
                return S.NaN

            if isinstance(expt, Integer) and expt.is_positive:
                if expt.is_odd:
                    return S.NegativeInfinity
                else:
                    return S.Infinity

            return S.NegativeOne**expt*S.Infinity**expt

    def _as_mpf_val(self, prec):
        return mlib.fninf

    def __hash__(self):
        return super(NegativeInfinity, self).__hash__()

    def __eq__(self, other):
        return other is S.NegativeInfinity

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s < %s" % (self, other))
        if other.is_extended_real:
            if other.is_finite or other is S.Infinity:
                return S.true
            elif other.is_nonnegative:
                return S.true
            elif other is S.NegativeInfinity:
                return S.false
        return Expr.__lt__(self, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s <= %s" % (self, other))
        if other.is_extended_real:
            return S.true
        return Expr.__le__(self, other)

    def __gt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s > %s" % (self, other))
        if other.is_extended_real:
            return S.false
        return Expr.__gt__(self, other)

    def __ge__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s >= %s" % (self, other))
        if other.is_extended_real:
            if other.is_finite or other is S.Infinity:
                return S.false
            elif other.is_nonnegative:
                return S.false
            elif other is S.NegativeInfinity:
                return S.true
        return Expr.__ge__(self, other)

    def __mod__(self, other):
        return S.NaN

    __rmod__ = __mod__


class NaN(Number, metaclass=Singleton):
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

    NaN is a singleton, and can be accessed by ``S.NaN``, or can be imported
    as ``nan``.

    Examples
    ========

    >>> from diofant import nan, S, oo, Eq
    >>> nan is S.NaN
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

    .. [1] http://en.wikipedia.org/wiki/NaN

    """
    is_commutative = True
    is_real = None
    is_rational = None
    is_algebraic = None
    is_transcendental = None
    is_integer = None
    is_comparable = False
    is_finite = None
    is_zero = None
    is_prime = None
    is_positive = None
    is_negative = None
    is_number = True

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    def _latex(self, printer):
        return r"\mathrm{NaN}"

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
    def __div__(self, other):
        return self

    __truediv__ = __div__

    def _as_mpf_val(self, prec):
        return _mpf_nan

    def __hash__(self):
        return super(NaN, self).__hash__()

    def __eq__(self, other):
        # NaN is structurally equal to another NaN
        return other is S.NaN

    def _eval_Eq(self, other):
        # NaN is not mathematically equal to anything, even NaN
        return S.false

    # Expr will _sympify and raise TypeError
    __gt__ = Expr.__gt__
    __ge__ = Expr.__ge__
    __lt__ = Expr.__lt__
    __le__ = Expr.__le__

nan = S.NaN


class ComplexInfinity(AtomicExpr, metaclass=Singleton):
    r"""Complex infinity.

    In complex analysis the symbol `\tilde\infty`, called "complex
    infinity", represents a quantity with infinite magnitude, but
    undetermined complex phase.

    ComplexInfinity is a singleton, and can be accessed by
    ``S.ComplexInfinity``, or can be imported as ``zoo``.

    Examples
    ========

    >>> from diofant import zoo, oo
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

    def _latex(self, printer):
        return r"\tilde{\infty}"

    @staticmethod
    def __abs__():
        return S.Infinity

    @staticmethod
    def __neg__():
        return S.ComplexInfinity

    def _eval_power(self, expt):
        if expt is S.ComplexInfinity:
            return S.NaN

        if isinstance(expt, Number):
            if expt is S.Zero:
                return S.NaN
            else:
                if expt.is_positive:
                    return S.ComplexInfinity
                else:
                    return S.Zero


zoo = S.ComplexInfinity


class NumberSymbol(AtomicExpr):

    is_commutative = True
    is_finite = True
    is_number = True

    is_NumberSymbol = True

    def __new__(cls):
        return AtomicExpr.__new__(cls)

    def approximation_interval(self, number_cls):
        """ Return an interval with number_cls endpoints that contains the
        value of NumberSymbol.  If not implemented, then return None.
        """
        return  # pragma: no cover

    def _eval_evalf(self, prec):
        return Float._new(self._as_mpf_val(prec), prec)

    def __eq__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            return False    # diofant != other  -->  not ==
        if self is other:
            return True
        if isinstance(other, Number) and self.is_irrational:
            return False

        return False    # NumberSymbol != non-(Number|self)

    def __lt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s < %s" % (self, other))
        if self is other:
            return S.false
        if isinstance(other, Number):
            approx = self.approximation_interval(other.__class__)
            if approx is not None:
                l, u = approx
                if other < l:
                    return S.false
                if other > u:
                    return S.true
            return _sympify(self.evalf() < other)
        if other.is_extended_real and other.is_number:
            other = other.evalf()
            return _sympify(self.evalf() < other)
        return Expr.__lt__(self, other)

    def __le__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s <= %s" % (self, other))
        if self is other:
            return S.true
        if other.is_extended_real and other.is_number:
            other = other.evalf()
        if isinstance(other, Number):
            return _sympify(self.evalf() <= other)
        return Expr.__le__(self, other)

    def __gt__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s > %s" % (self, other))
        r = _sympify((-self) < (-other))
        if r in (S.true, S.false):
            return r
        else:
            return Expr.__gt__(self, other)

    def __ge__(self, other):
        try:
            other = _sympify(other)
        except SympifyError:
            raise TypeError("Invalid comparison %s >= %s" % (self, other))
        r = _sympify((-self) <= (-other))
        if r in (S.true, S.false):
            return r
        else:
            return Expr.__ge__(self, other)

    def __int__(self):
        raise NotImplementedError  # pragma: no cover

    def __hash__(self):
        return super(NumberSymbol, self).__hash__()


class Exp1(NumberSymbol, metaclass=Singleton):
    r"""The `e` constant.

    The transcendental number `e = 2.718281828\dots` is the base of the
    natural logarithm and of the exponential function, `e = \exp(1)`.
    Sometimes called Euler's number or Napier's constant.

    Exp1 is a singleton, and can be accessed by ``S.Exp1``,
    or can be imported as ``E``.

    Examples
    ========

    >>> from diofant import exp, log, E
    >>> E is exp(1)
    True
    >>> log(E)
    1

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/E_%28mathematical_constant%29
    """

    is_real = True
    is_positive = True
    is_negative = False  # XXX Forces is_negative/is_nonnegative
    is_irrational = True
    is_number = True
    is_algebraic = False
    is_transcendental = True

    def _latex(self, printer):
        return r"e"

    @staticmethod
    def __abs__():
        return S.Exp1

    def __int__(self):
        return 2

    def _as_mpf_val(self, prec):
        return mpf_e(prec)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return Integer(2), Integer(3)

    def _eval_power(self, arg):
        from diofant.functions.elementary.exponential import log
        from diofant import Add, Mul, Pow
        if arg.is_Number:
            if arg is S.Infinity:
                return S.Infinity
            elif arg is S.NegativeInfinity:
                return S.Zero
        elif arg is S.ComplexInfinity:
            return S.NaN
        elif arg.func is log:
            return arg.args[0]
        elif arg.is_Mul:
            Ioo = S.ImaginaryUnit*S.Infinity
            if arg in [Ioo, -Ioo]:
                return S.NaN

            coeff = arg.coeff(S.Pi*S.ImaginaryUnit)
            if coeff:
                if (2*coeff).is_integer:
                    if coeff.is_even:
                        return S.One
                    elif coeff.is_odd:
                        return S.NegativeOne
                    elif (coeff + S.Half).is_even:
                        return -S.ImaginaryUnit
                    elif (coeff + S.Half).is_odd:
                        return S.ImaginaryUnit

            # Warning: code in risch.py will be very sensitive to changes
            # in this (see DifferentialExtension).

            # look for a single log factor

            coeff, terms = arg.as_coeff_Mul()

            # but it can't be multiplied by oo
            if coeff in [S.NegativeInfinity, S.Infinity]:
                return None

            coeffs, log_term = [coeff], None
            for term in Mul.make_args(terms):
                if term.func is log:
                    if log_term is None:
                        log_term = term.args[0]
                    else:
                        return None
                elif term.is_comparable:
                    coeffs.append(term)
                else:
                    return None

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
        from diofant import sin
        I = S.ImaginaryUnit
        return sin(I + S.Pi/2) - I*sin(I)

    def _eval_rewrite_as_cos(self):
        from diofant import cos
        I = S.ImaginaryUnit
        return cos(I) + I*cos(I + S.Pi/2)

E = S.Exp1


class Pi(NumberSymbol, metaclass=Singleton):
    r"""The `\pi` constant.

    The transcendental number `\pi = 3.141592654\dots` represents the ratio
    of a circle's circumference to its diameter, the area of the unit circle,
    the half-period of trigonometric functions, and many other things
    in mathematics.

    Pi is a singleton, and can be accessed by ``S.Pi``, or can
    be imported as ``pi``.

    Examples
    ========

    >>> from diofant import S, pi, oo, sin, exp, integrate, Symbol
    >>> S.Pi
    pi
    >>> pi > 3
    true
    >>> pi.is_irrational
    True
    >>> x = Symbol('x')
    >>> sin(x + 2*pi)
    sin(x)
    >>> integrate(exp(-x**2), (x, -oo, oo))
    sqrt(pi)

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Pi
    """

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = True
    is_number = True
    is_algebraic = False
    is_transcendental = True

    def _latex(self, printer):
        return r"\pi"

    @staticmethod
    def __abs__():
        return S.Pi

    def __int__(self):
        return 3

    def _as_mpf_val(self, prec):
        return mpf_pi(prec)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return Integer(3), Integer(4)
        elif issubclass(number_cls, Rational):
            return Rational(223, 71), Rational(22, 7)

pi = S.Pi


class GoldenRatio(NumberSymbol, metaclass=Singleton):
    r"""The golden ratio, `\phi`.

    `\phi = \frac{1 + \sqrt{5}}{2}` is algebraic number.  Two quantities
    are in the golden ratio if their ratio is the same as the ratio of
    their sum to the larger of the two quantities, i.e. their maximum.

    GoldenRatio is a singleton, and can be accessed by ``S.GoldenRatio``.

    Examples
    ========

    >>> from diofant import S
    >>> S.GoldenRatio > 1
    true
    >>> S.GoldenRatio.expand(func=True)
    1/2 + sqrt(5)/2
    >>> S.GoldenRatio.is_irrational
    True

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Golden_ratio
    """

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = True
    is_number = True
    is_algebraic = True
    is_transcendental = False

    def _latex(self, printer):
        return r"\phi"

    def __int__(self):
        return 1

    def _as_mpf_val(self, prec):
        # XXX track down why this has to be increased
        rv = mlib.from_man_exp(phi_fixed(prec + 10), -prec - 10)
        return mpf_norm(rv, prec)

    def _eval_expand_func(self, **hints):
        from diofant import sqrt
        return S.Half + S.Half*sqrt(5)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return S.One, Rational(2)


class EulerGamma(NumberSymbol, metaclass=Singleton):
    r"""The Euler-Mascheroni constant.

    `\gamma = 0.5772157\dots` (also called Euler's constant) is a mathematical
    constant recurring in analysis and number theory.  It is defined as the
    limiting difference between the harmonic series and the
    natural logarithm:

    .. math:: \gamma = \lim\limits_{n\to\infty}
              \left(\sum\limits_{k=1}^n\frac{1}{k} - \ln n\right)

    EulerGamma is a singleton, and can be accessed by ``S.EulerGamma``.

    Examples
    ========

    >>> from diofant import S
    >>> S.EulerGamma.is_irrational
    >>> S.EulerGamma > 0
    true
    >>> S.EulerGamma > 1
    false

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant
    """

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = None
    is_number = True

    def _latex(self, printer):
        return r"\gamma"

    def __int__(self):
        return 0

    def _as_mpf_val(self, prec):
        # XXX track down why this has to be increased
        v = mlib.libhyper.euler_fixed(prec + 10)
        rv = mlib.from_man_exp(v, -prec - 10)
        return mpf_norm(rv, prec)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return S.Zero, S.One
        elif issubclass(number_cls, Rational):
            return S.Half, Rational(3, 5)


class Catalan(NumberSymbol, metaclass=Singleton):
    r"""Catalan's constant.

    `K = 0.91596559\dots` is given by the infinite series

    .. math:: K = \sum_{k=0}^{\infty} \frac{(-1)^k}{(2k+1)^2}

    Catalan is a singleton, and can be accessed by ``S.Catalan``.

    Examples
    ========

    >>> from diofant import S
    >>> S.Catalan.is_irrational
    >>> S.Catalan > 0
    true
    >>> S.Catalan > 1
    false

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Catalan%27s_constant
    """

    is_real = True
    is_positive = True
    is_negative = False
    is_irrational = None
    is_number = True

    def __int__(self):
        return 0

    def _as_mpf_val(self, prec):
        # XXX track down why this has to be increased
        v = mlib.catalan_fixed(prec + 10)
        rv = mlib.from_man_exp(v, -prec - 10)
        return mpf_norm(rv, prec)

    def approximation_interval(self, number_cls):
        if issubclass(number_cls, Integer):
            return S.Zero, S.One
        elif issubclass(number_cls, Rational):
            return Rational(9, 10), S.One


class ImaginaryUnit(AtomicExpr, metaclass=Singleton):
    r"""The imaginary unit, `i = \sqrt{-1}`.

    I is a singleton, and can be accessed by ``S.I``, or can be
    imported as ``I``.

    Examples
    ========

    >>> from diofant import I, sqrt
    >>> sqrt(-1)
    I
    >>> I*I
    -1
    >>> 1/I
    -I

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Imaginary_unit
    """

    is_commutative = True
    is_imaginary = True
    is_extended_real = False
    is_finite = True
    is_number = True
    is_algebraic = True
    is_transcendental = False
    is_real = False

    def _latex(self, printer):
        return r"i"

    @staticmethod
    def __abs__():
        return S.One

    def _eval_evalf(self, prec):
        return self

    def _eval_conjugate(self):
        return -S.ImaginaryUnit

    def _eval_power(self, expt):
        """
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
                expt = expt.p % 4
                if expt == 0:
                    return S.One
                if expt == 1:
                    return S.ImaginaryUnit
                if expt == 2:
                    return -S.One
                return -S.ImaginaryUnit
            return (S.NegativeOne)**(expt*S.Half)
        return

    def as_base_exp(self):
        return S.NegativeOne, S.Half

I = S.ImaginaryUnit


def sympify_fractions(f):
    return Rational(f.numerator, f.denominator)

converter[fractions.Fraction] = sympify_fractions


try:
    if HAS_GMPY:
        import gmpy2 as gmpy
    else:
        raise ImportError

    def sympify_mpz(x):
        return Integer(int(x))

    def sympify_mpq(x):
        return Rational(int(x.numerator), int(x.denominator))

    converter[type(gmpy.mpz(1))] = sympify_mpz
    converter[type(gmpy.mpq(1, 2))] = sympify_mpq
except ImportError:
    pass


def sympify_mpmath(x):
    return Expr._from_mpmath(x, x.context.prec)

converter[mpnumeric] = sympify_mpmath


def sympify_complex(a):
    real, imag = list(map(sympify, (a.real, a.imag)))
    return real + S.ImaginaryUnit*imag

converter[complex] = sympify_complex

from .power import Pow, integer_nthroot
from .mul import Mul
Mul.identity = One()
from .add import Add
Add.identity = Zero()
