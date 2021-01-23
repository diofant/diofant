import typing

from ...core import (Add, Function, Integer, Rational, Symbol, cacheit,
                     expand_mul)
from ...core.function import ArgumentIndexError
from ...core.logic import fuzzy_and, fuzzy_not
from ...core.numbers import I, igcdex, nan, oo, pi, zoo
from ...core.sympify import sympify
from ...utilities import numbered_symbols
from ..combinatorial.factorials import RisingFactorial, factorial
from .exponential import exp, log
from .hyperbolic import acoth, asinh, atanh, cosh, coth, csch, sech, sinh, tanh
from .miscellaneous import sqrt


###############################################################################
# ######################## TRIGONOMETRIC FUNCTIONS ########################## #
###############################################################################


class TrigonometricFunction(Function):
    """Base class for trigonometric functions."""

    unbranched = True

    def _eval_expand_complex(self, deep=True, **hints):
        re_part, im_part = self.as_real_imag(deep=deep, **hints)
        return re_part + im_part*I

    def _as_real_imag(self, deep=True, **hints):
        if self.args[0].is_extended_real:
            if deep:
                hints['complex'] = False
                return self.args[0].expand(deep, **hints), Integer(0)
            else:
                return self.args[0], Integer(0)
        if deep:
            re, im = self.args[0].expand(deep, **hints).as_real_imag()
        else:
            re, im = self.args[0].as_real_imag()
        return re, im


def _peeloff_pi(arg):
    """
    Split ARG into two parts, a "rest" and a multiple of pi/2.
    This assumes ARG to be an Add.
    The multiple of pi returned in the second position is always a Rational.

    Examples
    ========

    >>> _peeloff_pi(x + pi/2)
    (x, pi/2)
    >>> _peeloff_pi(x + 2*pi/3 + pi*y)
    (x + pi*y + pi/6, pi/2)

    """
    for a in Add.make_args(arg):
        if a is pi:
            K = Integer(1)
            break
        elif a.is_Mul:
            K, p = a.as_two_terms()
            if p is pi and K.is_Rational:
                break
    else:
        return arg, Integer(0)

    m1 = (K % Rational(1, 2)) * pi
    m2 = K*pi - m1
    return arg - m2, m2


def _pi_coeff(arg, cycles=1):
    """
    When arg is a Number times pi (e.g. 3*pi/2) then return the Number
    normalized to be in the range [0, 2], else None.

    When an even multiple of pi is encountered, if it is multiplying
    something with known parity then the multiple is returned as 0 otherwise
    as 2.

    Examples
    ========

    >>> _pi_coeff(3*x*pi)
    3*x
    >>> _pi_coeff(11*pi/7)
    11/7
    >>> _pi_coeff(-11*pi/7)
    3/7
    >>> _pi_coeff(4*pi)
    0
    >>> _pi_coeff(5*pi)
    1
    >>> _pi_coeff(5.0*pi)
    1
    >>> _pi_coeff(5.5*pi)
    3/2
    >>> _pi_coeff(2 + pi)

    >>> _pi_coeff(2*Dummy(integer=True)*pi)
    2
    >>> _pi_coeff(2*Dummy(even=True)*pi)
    0

    """
    arg = sympify(arg)
    if arg is pi:
        return Integer(1)
    elif not arg:
        return Integer(0)
    elif arg.is_Mul:
        cx = arg.coeff(pi)
        if cx:
            c, x = cx.as_coeff_Mul()  # pi is not included as coeff
            if c.is_Float:
                # recast exact binary fractions to Rationals
                f = abs(c) % 1
                if f != 0:
                    p = -round(log(f, 2).evalf(strict=False))
                    m = 2**p
                    cm = c*m
                    i = int(cm)
                    if i == cm:
                        c = Rational(i, m)
                        cx = c*x
                else:
                    c = Rational(int(c))
                    cx = c*x
            if x.is_integer:
                c2 = c % 2
                if c2 == 1:
                    return x
                elif not c2:
                    if x.is_even is not None:  # known parity
                        return Integer(0)
                    return Integer(2)
                else:
                    return c2*x
            return cx


class sin(TrigonometricFunction):
    """
    The sine function.

    Returns the sine of x (measured in radians).

    Notes
    =====

    This function will evaluate automatically in the
    case x/pi is some rational number.  For example,
    if x is a multiple of pi, pi/2, pi/3, pi/4 and pi/6.

    Examples
    ========

    >>> sin(x**2).diff(x)
    2*x*cos(x**2)
    >>> sin(pi)
    0
    >>> sin(pi/2)
    1
    >>> sin(pi/6)
    1/2
    >>> sin(pi/12)
    -sqrt(2)/4 + sqrt(6)/4


    See Also
    ========

    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_functions
    * https://dlmf.nist.gov/4.14
    * http://functions.wolfram.com/ElementaryFunctions/Sin
    * https://mathworld.wolfram.com/TrigonometryAngles.html

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return cos(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg == 0:
                return Integer(0)
            elif arg in (oo, -oo):
                return

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        i_coeff = arg.as_coefficient(I)
        if i_coeff is not None:
            return I * sinh(i_coeff)

        pi_coeff = _pi_coeff(arg)
        if pi_coeff is not None:
            if pi_coeff.is_integer:
                return Integer(0)

            if (2*pi_coeff).is_odd:
                return (-1)**(pi_coeff - Rational(1, 2))

            if not pi_coeff.is_Rational:
                narg = pi_coeff*pi
                if narg != arg:
                    return cls(narg)
                return

            else:
                # https://github.com/sympy/sympy/issues/6048
                # transform a sine to a cosine, to avoid redundant code

                x = pi_coeff % 2
                if x > 1:
                    return -cls((x % 1)*pi)
                if 2*x > 1:
                    return cls((1 - x)*pi)
                narg = ((pi_coeff + Rational(3, 2)) % 2)*pi
                result = cos(narg)
                if not isinstance(result, cos):
                    return result
                if pi_coeff*pi != arg:
                    return cls(pi_coeff*pi)
                return

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                return sin(m)*cos(x) + cos(m)*sin(x)

        if isinstance(arg, asin):
            return arg.args[0]

        if isinstance(arg, atan):
            x = arg.args[0]
            return x / sqrt(1 + x**2)

        if isinstance(arg, atan2):
            y, x = arg.args
            return y / sqrt(x**2 + y**2)

        if isinstance(arg, acos):
            x = arg.args[0]
            return sqrt(1 - x**2)

        if isinstance(arg, acot):
            x = arg.args[0]
            return 1 / (sqrt(1 + 1 / x**2) * x)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)

            if len(previous_terms) >= 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n - 1))
            else:
                return (-1)**(n//2) * x**n/factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        return (exp(arg*I) - exp(-arg*I)) / (2*I)

    def _eval_rewrite_as_Pow(self, arg):
        if isinstance(arg, log):
            x = arg.args[0]
            return I*x**-I / 2 - I*x**I / 2

    def _eval_rewrite_as_cos(self, arg):
        return -cos(arg + pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = tan(arg/2)
        return 2*tan_half/(1 + tan_half**2)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)*cos(arg)/cos(arg)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = cot(arg/2)
        return 2*cot_half/(1 + cot_half**2)

    def _eval_rewrite_as_sqrt(self, arg):
        return self.rewrite(cos).rewrite(sqrt)

    def _eval_rewrite_as_csc(self, arg):
        return 1/csc(arg)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        re, im = self._as_real_imag(deep=deep, **hints)
        return sin(re)*cosh(im), cos(re)*sinh(im)

    def _eval_expand_trig(self, **hints):
        from .. import chebyshevt, chebyshevu
        arg = self.args[0]
        x = None
        if arg.is_Add:  # TODO, implement more if deep stuff here
            # TODO: Do this more efficiently for more than two terms
            x, y = arg.as_two_terms()
            sx = sin(x, evaluate=False)._eval_expand_trig()
            sy = sin(y, evaluate=False)._eval_expand_trig()
            cx = cos(x, evaluate=False)._eval_expand_trig()
            cy = cos(y, evaluate=False)._eval_expand_trig()
            return sx*cy + sy*cx
        else:
            n, x = arg.as_coeff_Mul(rational=True)
            if n.is_Integer:  # n will be positive because of .eval
                # canonicalization

                # See https://mathworld.wolfram.com/Multiple-AngleFormulas.html
                if n.is_odd:
                    return (-1)**((n - 1)/2)*chebyshevt(n, sin(x))
                else:
                    return expand_mul((-1)**(n/2 - 1)*cos(x)*chebyshevu(n -
                                                                        1, sin(x)), deep=False)
            pi_coeff = _pi_coeff(arg)
            if pi_coeff is not None:
                if pi_coeff.is_Rational:
                    return self.rewrite(sqrt)
        return sin(arg)

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_complex(self):
        if self.args[0].is_complex:
            return True

    def _eval_is_real(self):
        if self.args[0].is_real:
            return True

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_zero:
                return True
            elif s.args[0].is_rational and s.args[0].is_nonzero:
                return False
        else:
            return s.is_rational

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if self.args[0].is_algebraic and self.args[0].is_nonzero:
                return False
            pi_coeff = _pi_coeff(self.args[0])
            if pi_coeff is not None and pi_coeff.is_rational:
                return True
        else:
            return s.is_algebraic


class cos(TrigonometricFunction):
    """
    The cosine function.

    Returns the cosine of x (measured in radians).

    Notes
    =====

    See :func:`~diofant.functions.elementary.trigonometric.sin`
    for notes about automatic evaluation.

    Examples
    ========

    >>> cos(x**2).diff(x)
    -2*x*sin(x**2)
    >>> cos(pi)
    -1
    >>> cos(pi/2)
    0
    >>> cos(2*pi/3)
    -1/2
    >>> cos(pi/12)
    sqrt(2)/4 + sqrt(6)/4

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_functions
    * https://dlmf.nist.gov/4.14
    * http://functions.wolfram.com/ElementaryFunctions/Cos

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -sin(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        from .. import chebyshevt
        if arg.is_Number:
            if arg == 0:
                return Integer(1)
            elif arg in (oo, -oo):
                # In this cases, it is unclear if we should
                # return nan or leave un-evaluated.  One
                # useful test case is how "limit(sin(x)/x,x,oo)"
                # is handled.
                # See test_sin_cos_with_infinity() an
                # Test for issue sympy/sympy#3308
                # https://github.com/sympy/sympy/issues/5196
                # For now, we return un-evaluated.
                return

        if arg.could_extract_minus_sign():
            return cls(-arg)

        i_coeff = arg.as_coefficient(I)
        if i_coeff is not None:
            return cosh(i_coeff)

        pi_coeff = _pi_coeff(arg)
        if pi_coeff is not None:
            if pi_coeff.is_integer:
                return (-1)**pi_coeff

            if (2*pi_coeff).is_odd:
                return Integer(0)

            if not pi_coeff.is_Rational:
                narg = pi_coeff*pi
                if narg != arg:
                    return cls(narg)
                return

            else:

                # cosine formula #####################
                # https://github.com/sympy/sympy/issues/6048
                # explicit calculations are preformed for
                # cos(k pi/n) for n = 8,10,12,15,20,24,30,40,60,120
                # Some other exact values like cos(k pi/240) can be
                # calculated using a partial-fraction decomposition
                # by calling cos( X ).rewrite(sqrt)
                cst_table_some = {3: Rational(1, 2),
                                  5: (sqrt(5) + 1)/4}

                q = pi_coeff.denominator
                p = pi_coeff.numerator % (2*q)
                if p > q:
                    narg = (pi_coeff - 1)*pi
                    return -cls(narg)
                if 2*p > q:
                    narg = (1 - pi_coeff)*pi
                    return -cls(narg)

                # If nested sqrt's are worse than un-evaluation
                # you can require q to be in (1, 2, 3, 4, 6, 12)
                # q <= 12, q=15, q=20, q=24, q=30, q=40, q=60, q=120 return
                # expressions with 2 or fewer sqrt nestings.
                table2 = {
                    12: (3, 4),
                    20: (4, 5),
                    30: (5, 6),
                    15: (6, 10),
                    24: (6, 8),
                    40: (8, 10),
                    60: (20, 30),
                    120: (40, 60)
                }
                if q in table2:
                    a, b = p*pi/table2[q][0], p*pi/table2[q][1]
                    nvala, nvalb = cls(a), cls(b)
                    assert None not in (nvala, nvalb)
                    return nvala*nvalb + cls(pi/2 - a)*cls(pi/2 - b)

                if q > 12:
                    return

                if q in cst_table_some:
                    cts = cst_table_some[pi_coeff.denominator]
                    return chebyshevt(pi_coeff.numerator, cts).expand()

                if 0 == q % 2:
                    narg = (pi_coeff*2)*pi
                    nval = cls(narg)
                    assert nval is not None
                    x = (2*pi_coeff + 1)/2
                    sign_cos = (-1)**((-1 if x < 0 else 1)*int(abs(x)))
                    return sign_cos*sqrt( (1 + nval)/2 )
            return

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                return cos(m)*cos(x) - sin(m)*sin(x)

        if isinstance(arg, acos):
            return arg.args[0]

        if isinstance(arg, atan):
            x = arg.args[0]
            return 1 / sqrt(1 + x**2)

        if isinstance(arg, atan2):
            y, x = arg.args
            return x / sqrt(x**2 + y**2)

        if isinstance(arg, asin):
            x = arg.args[0]
            return sqrt(1 - x ** 2)

        if isinstance(arg, acot):
            x = arg.args[0]
            return 1 / sqrt(1 + 1 / x**2)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 1:
            return Integer(0)
        else:
            x = sympify(x)

            if len(previous_terms) >= 2:
                p = previous_terms[-2]
                return -p * x**2 / (n*(n - 1))
            else:
                return (-1)**(n//2)*x**n/factorial(n)

    def _eval_rewrite_as_exp(self, arg):
        return (exp(arg*I) + exp(-arg*I)) / 2

    def _eval_rewrite_as_Pow(self, arg):
        if isinstance(arg, log):
            x = arg.args[0]
            return x**I/2 + x**-I/2

    def _eval_rewrite_as_sin(self, arg):
        return sin(arg + pi/2)

    def _eval_rewrite_as_tan(self, arg):
        tan_half = tan(arg/2)**2
        return (1 - tan_half)/(1 + tan_half)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)*cos(arg)/sin(arg)

    def _eval_rewrite_as_cot(self, arg):
        cot_half = cot(arg/2)**2
        return (cot_half - 1)/(cot_half + 1)

    def _eval_rewrite_as_sqrt(self, arg):
        from .. import chebyshevt

        def migcdex(x):
            # recursive calcuation of gcd and linear combination
            # for a sequence of integers.
            # Given  (x1, x2, x3)
            # Returns (y1, y1, y3, g)
            # such that g is the gcd and x1*y1+x2*y2+x3*y3 - g = 0
            # Note, that this is only one such linear combination.
            assert len(x) != 1
            if len(x) == 2:
                return igcdex(x[0], x[-1])
            g = migcdex(x[1:])
            u, v, h = igcdex(x[0], g[-1])
            return tuple([u] + [v*i for i in g[0:-1] ] + [h])

        def ipartfrac(r, factors=None):
            from ...ntheory import factorint
            assert isinstance(r, Rational)
            n = r.denominator
            assert 2 < r.denominator*r.denominator
            if factors is None:
                a = [n//x**y for x, y in factorint(r.denominator).items()]
            else:
                a = [n//x for x in factors]
            if len(a) == 1:
                return [ r ]
            h = migcdex(a)
            ans = [r.numerator*Rational(i*j, r.denominator) for i, j in zip(h[:-1], a)]
            assert r == sum(ans)
            return ans
        pi_coeff = _pi_coeff(arg)
        if pi_coeff is None:
            return

        assert not pi_coeff.is_integer, 'should have been simplified already'

        if not pi_coeff.is_Rational:
            return

        cst_table_some = {
            3: Rational(1, 2),
            5: (sqrt(5) + 1)/4,
            17: sqrt((15 + sqrt(17))/32 + sqrt(2)*(sqrt(17 - sqrt(17)) +
                                                   sqrt(sqrt(2)*(-8*sqrt(17 + sqrt(17)) - (1 - sqrt(17))
                                                                 * sqrt(17 - sqrt(17))) + 6*sqrt(17) + 34))/32)
            # 65537 and 257 are the only other known Fermat primes
            # Please add if you would like them
        }

        def fermatCoords(n):
            assert n > 1 and n % 2
            primes = {p: 0 for p in cst_table_some}
            assert 1 not in primes
            for p_i in primes:
                while 0 == n % p_i:
                    n = n/p_i
                    primes[p_i] += 1
            if 1 != n:
                return False
            if max(primes.values()) > 1:
                return False
            return tuple(p for p in primes if primes[p] == 1)

        if pi_coeff.denominator in cst_table_some:
            return chebyshevt(pi_coeff.numerator, cst_table_some[pi_coeff.denominator]).expand()

        if 0 == pi_coeff.denominator % 2:  # recursively remove powers of 2
            narg = (pi_coeff*2)*pi
            nval = cos(narg)
            assert nval is not None
            nval = nval.rewrite(sqrt)
            x = (2*pi_coeff + 1)/2
            sign_cos = (-1)**((-1 if x < 0 else 1)*int(abs(x)))
            return sign_cos*sqrt( (1 + nval)/2 )

        FC = fermatCoords(pi_coeff.denominator)
        if FC:
            decomp = ipartfrac(pi_coeff, FC)
            X = [(x[1], x[0]*pi) for x in zip(decomp, numbered_symbols('z'))]
            pcls = cos(sum(x[0] for x in X))._eval_expand_trig().subs(X)
            return pcls.rewrite(sqrt)
        else:
            decomp = ipartfrac(pi_coeff)
            X = [(x[1], x[0]*pi) for x in zip(decomp, numbered_symbols('z'))]
            pcls = cos(sum(x[0] for x in X))._eval_expand_trig().subs(X)
            return pcls

    def _eval_rewrite_as_sec(self, arg):
        return 1/sec(arg)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        re, im = self._as_real_imag(deep=deep, **hints)
        return cos(re)*cosh(im), -sin(re)*sinh(im)

    def _eval_expand_trig(self, **hints):
        from .. import chebyshevt
        arg = self.args[0]
        x = None
        if arg.is_Add:  # TODO: Do this more efficiently for more than two terms
            x, y = arg.as_two_terms()
            sx = sin(x, evaluate=False)._eval_expand_trig()
            sy = sin(y, evaluate=False)._eval_expand_trig()
            cx = cos(x, evaluate=False)._eval_expand_trig()
            cy = cos(y, evaluate=False)._eval_expand_trig()
            return cx*cy - sx*sy
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff.is_Integer:
                return chebyshevt(coeff, cos(terms))
            pi_coeff = _pi_coeff(arg)
            if pi_coeff is not None:
                if pi_coeff.is_Rational:
                    return self.rewrite(sqrt)
        return cos(arg)

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return Integer(1)
        else:
            return self.func(arg)

    def _eval_is_real(self):
        if self.args[0].is_real:
            return True

    def _eval_is_complex(self):
        if self.args[0].is_complex:
            return True

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_zero:
                return True
            if s.args[0].is_rational and s.args[0].is_nonzero:
                return False
        else:
            return s.is_rational

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if self.args[0].is_algebraic and self.args[0].is_nonzero:
                return False
            pi_coeff = _pi_coeff(self.args[0])
            if pi_coeff is not None and pi_coeff.is_rational:
                return True
        else:
            return s.is_algebraic


class tan(TrigonometricFunction):
    """
    The tangent function.

    Returns the tangent of x (measured in radians).

    Notes
    =====

    See :func:`~diofant.functions.elementary.trigonometric.sin`
    for notes about automatic evaluation.

    Examples
    ========

    >>> tan(x**2).diff(x)
    2*x*(tan(x**2)**2 + 1)
    >>> tan(pi/8).expand()
    -1 + sqrt(2)

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_functions
    * https://dlmf.nist.gov/4.14
    * http://functions.wolfram.com/ElementaryFunctions/Tan

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1 + self**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return atan

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg == 0:
                return Integer(0)

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        i_coeff = arg.as_coefficient(I)
        if i_coeff is not None:
            return I * tanh(i_coeff)

        pi_coeff = _pi_coeff(arg, 2)
        if pi_coeff is not None:
            if pi_coeff.is_integer:
                return Integer(0)

            if not pi_coeff.is_Rational:
                narg = pi_coeff*pi
                if narg != arg:
                    return cls(narg)
                return

            else:
                if not pi_coeff.denominator % 2:
                    narg = pi_coeff*pi*2
                    cresult, sresult = cos(narg), cos(narg - pi/2)
                    if not isinstance(cresult, cos) \
                            and not isinstance(sresult, cos):
                        if sresult == 0:
                            return zoo
                        return (1 - cresult)/sresult
                table2 = {
                    12: (3, 4),
                    20: (4, 5),
                    30: (5, 6),
                    15: (6, 10),
                    24: (6, 8),
                    40: (8, 10),
                    60: (20, 30),
                    120: (40, 60)
                }
                q = pi_coeff.denominator
                p = pi_coeff.numerator % q
                if q in table2:
                    nvala, nvalb = cls(p*pi/table2[q][0]), cls(p*pi/table2[q][1])
                    assert None not in (nvala, nvalb)
                    return (nvala - nvalb)/(1 + nvala*nvalb)
                narg = ((pi_coeff + Rational(1, 2)) % 1 - Rational(1, 2))*pi
                # see cos() to specify which expressions should  be
                # expanded automatically in terms of radicals
                cresult, sresult = cos(narg), cos(narg - pi/2)
                if not isinstance(cresult, cos) and not isinstance(sresult, cos):
                    return sresult/cresult
                if narg != arg:
                    return cls(narg)

        if arg.is_Add:
            x, m = _peeloff_pi(arg)
            if m:
                tanm = tan(m)
                tanx = tan(x)
                if tanm is zoo:
                    return -cot(x)
                return (tanm + tanx)/(1 - tanm*tanx)

        if isinstance(arg, atan):
            return arg.args[0]

        if isinstance(arg, atan2):
            y, x = arg.args
            return y/x

        if isinstance(arg, asin):
            x = arg.args[0]
            return x / sqrt(1 - x**2)

        if isinstance(arg, acos):
            x = arg.args[0]
            return sqrt(1 - x**2) / x

        if isinstance(arg, acot):
            x = arg.args[0]
            return 1 / x

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from .. import bernoulli
        if n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)

            a, b = ((n - 1)//2), 2**(n + 1)

            B = bernoulli(n + 1)
            F = factorial(n + 1)

            return (-1)**a * b*(b - 1) * B/F * x**n

    def _eval_nseries(self, x, n, logx):
        i = self.args[0].limit(x, 0)*2/pi
        if i and i.is_Integer:
            return self.rewrite(cos)._eval_nseries(x, n=n, logx=logx)
        return Function._eval_nseries(self, x, n=n, logx=logx)

    def _eval_rewrite_as_Pow(self, arg):
        if isinstance(arg, log):
            x = arg.args[0]
            return I*(x**-I - x**I)/(x**-I + x**I)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        re, im = self._as_real_imag(deep=deep, **hints)
        if im:
            denom = cos(2*re) + cosh(2*im)
            return sin(2*re)/denom, sinh(2*im)/denom
        else:
            return self.func(re), Integer(0)

    def _eval_expand_trig(self, **hints):
        from .complexes import im, re
        arg = self.args[0]
        if arg.is_Add:
            from ...polys import symmetric_poly
            n = len(arg.args)
            TX = []
            for x in arg.args:
                tx = tan(x, evaluate=False)._eval_expand_trig()
                TX.append(tx)

            Yg = numbered_symbols('Y')
            Y = [ next(Yg) for i in range(n) ]

            p = [0, 0]
            for i in range(n + 1):
                p[1 - i % 2] += symmetric_poly(i, *Y)*(-1)**((i % 4)//2)
            return (p[0]/p[1]).subs(list(zip(Y, TX)))

        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff.is_Integer and coeff > 1:
                z = Symbol('dummy', extended_real=True)
                P = ((1 + I*z)**coeff).expand()
                return (im(P)/re(P)).subs({z: tan(terms)})
        return tan(arg)

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(neg_exp - pos_exp)/(neg_exp + pos_exp)

    def _eval_rewrite_as_sin(self, x):
        return 2*sin(x)**2/sin(2*x)

    def _eval_rewrite_as_cos(self, x):
        return -cos(x + pi/2)/cos(x)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)/cos(arg)

    def _eval_rewrite_as_cot(self, arg):
        return 1/cot(arg)

    def _eval_rewrite_as_sqrt(self, arg):
        y = self.rewrite(cos).rewrite(sqrt)
        if not y.has(cos):
            return y

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_real(self):
        if (2*self.args[0]/pi).is_noninteger:
            return True

    def _eval_is_finite(self):
        if self.args[0].is_imaginary:
            return True

    def _eval_is_nonzero(self):
        if (2*self.args[0]/pi).is_noninteger:
            return True
        elif self.args[0].is_imaginary and self.args[0].is_nonzero:
            return True

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_zero:
                return True
            elif s.args[0].is_rational and s.args[0].is_nonzero:
                return False
        else:
            return s.is_rational

    def _eval_is_algebraic(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_algebraic and s.args[0].is_nonzero:
                return False
            pi_coeff = _pi_coeff(s.args[0])
            if pi_coeff is not None and pi_coeff.is_rational:
                if (2*pi_coeff).is_noninteger:
                    return True
        else:
            return s.is_algebraic


class ReciprocalTrigonometricFunction(TrigonometricFunction):
    """Base class for reciprocal functions of trigonometric functions."""

    _reciprocal_of = None       # mandatory, to be defined in subclass

    # _is_even and _is_odd are used for correct evaluation of csc(-x), sec(-x)
    # TODO refactor into TrigonometricFunction common parts of
    # trigonometric functions eval() like even/odd, func(x+2*k*pi), etc.
    _is_even: typing.Optional[bool] = None  # optional, to be defined in subclass
    _is_odd: typing.Optional[bool] = None  # optional, to be defined in subclass

    @classmethod
    def eval(cls, arg):
        if arg.could_extract_minus_sign():
            if cls._is_even:
                return cls(-arg)
            elif cls._is_odd:  # pragma: no branch XXX: make _is_odd mandatory?
                return -cls(-arg)

        i_coeff = arg.as_coefficient(I)
        if i_coeff is not None:
            c, n = Integer(1), sech
            if cls._reciprocal_of == sin:
                c, n = -I, csch
            elif cls._reciprocal_of == tan:
                c, n = -I, coth
            return c*n(i_coeff)

        pi_coeff = _pi_coeff(arg)
        if pi_coeff is not None:
            narg = pi_coeff*pi
            if narg != arg:
                return cls(narg)
            elif pi_coeff.is_Number and pi_coeff > 1:
                s = Integer(1) if cls is cot else Integer(-1)
                return s*cls((pi_coeff - 1)*pi)

        t = cls._reciprocal_of.eval(arg)
        if hasattr(arg, 'inverse') and arg.inverse() == cls:
            return arg.args[0]

        if t is not None:
            t = 1/t
            if cls is cot:
                if cls in ((1/t).func, (-1/t).func):
                    t = t.rewrite(cls._reciprocal_of)
                elif cls._reciprocal_of in ((1/t).func, (-1/t).func):
                    t = t.rewrite(cls)
            return t

    def _call_reciprocal(self, method_name, *args, **kwargs):
        # Calls method_name on _reciprocal_of
        o = self._reciprocal_of(self.args[0])
        return getattr(o, method_name)(*args, **kwargs)

    def _calculate_reciprocal(self, method_name, *args, **kwargs):
        # If calling method_name on _reciprocal_of returns a value != None
        # then return the reciprocal of that value
        t = self._call_reciprocal(method_name, *args, **kwargs)
        return 1/t if t is not None else t

    def _rewrite_reciprocal(self, method_name, arg):
        # Special handling for rewrite functions. If reciprocal rewrite returns
        # unmodified expression, then return None
        t = self._call_reciprocal(method_name, arg)
        if t is not None and t != self._reciprocal_of(arg):
            return 1/t

    def _eval_rewrite_as_exp(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_exp', arg)

    def _eval_rewrite_as_Pow(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_Pow', arg)

    def _eval_rewrite_as_sin(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_sin', arg)

    def _eval_rewrite_as_cos(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_cos', arg)

    def _eval_rewrite_as_tan(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_tan', arg)

    def _eval_rewrite_as_sqrt(self, arg):
        return self._rewrite_reciprocal('_eval_rewrite_as_sqrt', arg)

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate())

    def as_real_imag(self, deep=True, **hints):
        return (1/self._reciprocal_of(self.args[0])).as_real_imag(deep,
                                                                  **hints)

    def _eval_expand_trig(self, **hints):
        return self._calculate_reciprocal('_eval_expand_trig', **hints)

    def _eval_is_extended_real(self):
        return (1/self._reciprocal_of(self.args[0])).is_extended_real

    def _eval_is_finite(self):
        return (1/self._reciprocal_of(self.args[0])).is_finite

    def _eval_is_rational(self):
        return (1/self._reciprocal_of(self.args[0])).is_rational

    def _eval_is_algebraic(self):
        return (1/self._reciprocal_of(self.args[0])).is_algebraic

    def _eval_nseries(self, x, n, logx):
        return (1/self._reciprocal_of(self.args[0]))._eval_nseries(x, n, logx)


class sec(ReciprocalTrigonometricFunction):
    """
    The secant function.

    Returns the secant of x (measured in radians).

    Notes
    =====

    See :func:`~diofant.functions.elementary.trigonometric.sin`
    for notes about automatic evaluation.

    Examples
    ========

    >>> sec(x**2).diff(x)
    2*x*tan(x**2)*sec(x**2)

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_functions
    * https://dlmf.nist.gov/4.14
    * http://functions.wolfram.com/ElementaryFunctions/Sec

    """

    _reciprocal_of = cos
    _is_even = True

    def _eval_rewrite_as_cot(self, arg):
        cot_half_sq = cot(arg/2)**2
        return (cot_half_sq + 1)/(cot_half_sq - 1)

    def _eval_rewrite_as_cos(self, arg):
        return 1/cos(arg)

    def _eval_rewrite_as_sincos(self, arg):
        return sin(arg)/(cos(arg)*sin(arg))

    def fdiff(self, argindex=1):
        if argindex == 1:
            return tan(self.args[0])*sec(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        # Reference Formula:
        # http://functions.wolfram.com/ElementaryFunctions/Sec/06/01/02/01/
        from ..combinatorial.numbers import euler
        if n < 0 or n % 2 == 1:
            return Integer(0)
        else:
            x = sympify(x)
            k = n//2
            return (-1)**k*euler(2*k)/factorial(2*k)*x**(2*k)


class csc(ReciprocalTrigonometricFunction):
    """
    The cosecant function.

    Returns the cosecant of x (measured in radians).

    Notes
    =====

    See :func:`~diofant.functions.elementary.trigonometric.sin`
    for notes about automatic evaluation.

    Examples
    ========

    >>> csc(x**2).diff(x)
    -2*x*cot(x**2)*csc(x**2)

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_functions
    * https://dlmf.nist.gov/4.14
    * http://functions.wolfram.com/ElementaryFunctions/Csc

    """

    _reciprocal_of = sin
    _is_odd = True

    def _eval_rewrite_as_sin(self, arg):
        return 1/sin(arg)

    def _eval_rewrite_as_sincos(self, arg):
        return cos(arg)/(sin(arg)*cos(arg))

    def _eval_rewrite_as_cot(self, arg):
        cot_half = cot(arg/2)
        return (1 + cot_half**2)/(2*cot_half)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -cot(self.args[0])*csc(self.args[0])
        else:
            raise ArgumentIndexError(self, argindex)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        from .. import bernoulli
        if n == 0:
            return 1/sympify(x)
        elif n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            k = n//2 + 1
            return ((-1)**(k - 1)*2*(2**(2*k - 1) - 1) *
                    bernoulli(2*k)*x**(2*k - 1)/factorial(2*k))


class cot(ReciprocalTrigonometricFunction):
    """
    The cotangent function.

    Returns the cotangent of x (measured in radians).

    Notes
    =====

    See :func:`~diofant.functions.elementary.trigonometric.sin`
    for notes about automatic evaluation.

    Examples
    ========

    >>> cot(x**2).diff(x)
    2*x*(-cot(x**2)**2 - 1)

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Trigonometric_functions
    * https://dlmf.nist.gov/4.14
    * http://functions.wolfram.com/ElementaryFunctions/Cot

    """

    _reciprocal_of = tan
    _is_odd = True

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1 - self**2
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return acot

    def _eval_rewrite_as_tan(self, arg):
        return 1/tan(arg)

    def _eval_rewrite_as_sincos(self, arg):
        return cos(arg)/sin(arg)

    def _eval_rewrite_as_exp(self, arg):
        neg_exp, pos_exp = exp(-arg*I), exp(arg*I)
        return I*(pos_exp + neg_exp)/(pos_exp - neg_exp)

    def _eval_rewrite_as_sin(self, x):
        return 2*sin(2*x)/sin(x)**2

    def as_real_imag(self, deep=True, **hints):
        re, im = self._as_real_imag(deep=deep, **hints)
        if im:
            denom = cos(2*re) - cosh(2*im)
            return -sin(2*re)/denom, -sinh(2*im)/denom
        else:
            return self.func(re), Integer(0)

    def _eval_expand_trig(self, **hints):
        from .complexes import im, re
        arg = self.args[0]
        if arg.is_Add:
            from ...polys import symmetric_poly
            n = len(arg.args)
            CX = []
            for x in arg.args:
                cx = cot(x, evaluate=False)._eval_expand_trig()
                CX.append(cx)

            Yg = numbered_symbols('Y')
            Y = [ next(Yg) for i in range(n) ]

            p = [0, 0]
            for i in range(n, -1, -1):
                p[(n - i) % 2] += symmetric_poly(i, *Y)*(-1)**(((n - i) % 4)//2)
            return (p[0]/p[1]).subs(list(zip(Y, CX)))
        else:
            coeff, terms = arg.as_coeff_Mul(rational=True)
            if coeff.is_Integer and coeff > 1:
                z = Symbol('dummy', real=True)
                P = ((z + I)**coeff).expand()
                return (re(P)/im(P)).subs({z: cot(terms)})
        return cot(arg)

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return 1/arg
        else:
            return self.func(arg)

    def _eval_subs(self, old, new):
        arg = self.args[0]
        argnew = arg.subs({old: new})
        if arg != argnew and (argnew/pi).is_integer:
            return zoo
        return cot(argnew)


###############################################################################
# ######################### TRIGONOMETRIC INVERSES ########################## #
###############################################################################


class InverseTrigonometricFunction(Function):
    """Base class for inverse trigonometric functions."""


class asin(InverseTrigonometricFunction):
    """
    The inverse sine function.

    Returns the arcsine of x in radians.

    Notes
    =====

    asin(x) will evaluate automatically in the cases oo, -oo, 0, 1,
    -1 and for some instances when the result is a rational multiple
    of pi (see the eval class method).

    Examples
    ========

    >>> asin(1)
    pi/2
    >>> asin(-1)
    -pi/2

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Inverse_trigonometric_functions
    * https://dlmf.nist.gov/4.23
    * http://functions.wolfram.com/ElementaryFunctions/ArcSin

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg in (oo, -oo):
                return -arg * I
            elif arg == 0:
                return Integer(0)
            elif arg == 1:
                return pi / 2
            elif arg == -1:
                return -pi / 2

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        if arg.is_number:
            cst_table = {
                sqrt(3)/2: 3,
                -sqrt(3)/2: -3,
                sqrt(2)/2: 4,
                -sqrt(2)/2: -4,
                1/sqrt(2): 4,
                -1/sqrt(2): -4,
                sqrt((5 - sqrt(5))/8): 5,
                -sqrt((5 - sqrt(5))/8): -5,
                Rational(+1, 2): 6,
                Rational(-1, 2): -6,
                sqrt(2 - sqrt(2))/2: 8,
                -sqrt(2 - sqrt(2))/2: -8,
                (sqrt(5) - 1)/4: 10,
                (1 - sqrt(5))/4: -10,
                (sqrt(3) - 1)/sqrt(2**3): 12,
                (1 - sqrt(3))/sqrt(2**3): -12,
                (sqrt(5) + 1)/4: Rational(10, 3),
                -(sqrt(5) + 1)/4: -Rational(10, 3)
            }

            if arg in cst_table:
                return pi / cst_table[arg]

        i_coeff = arg.as_coefficient(I)
        if i_coeff is not None:
            return I * asinh(i_coeff)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
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
                return R / F * x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_rewrite_as_acos(self, x):
        return pi/2 - acos(x)

    def _eval_rewrite_as_atan(self, x):
        return 2*atan(x/(1 + sqrt(1 - x**2)))

    def _eval_rewrite_as_log(self, x):
        return -I*log(I*x + sqrt(1 - x**2))

    def _eval_rewrite_as_acot(self, arg):
        return 2*acot((1 + sqrt(1 - arg**2))/arg)

    def _eval_rewrite_as_asec(self, arg):
        return pi/2 - asec(1/arg)

    def _eval_rewrite_as_acsc(self, arg):
        return acsc(1/arg)

    def _eval_is_complex(self):
        if self.args[0].is_complex:
            return True

    def _eval_is_extended_real(self):
        x = self.args[0]
        return fuzzy_and([x.is_real, (1 - abs(x)).is_nonnegative])

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_zero:
                return True
            elif s.args[0].is_rational and s.args[0].is_nonzero:
                return False
        else:
            return s.is_rational

    def _eval_is_positive(self):
        if self.args[0].is_positive:
            return (self.args[0] - 1).is_negative
        if self.args[0].is_negative:
            return fuzzy_not((self.args[0] + 1).is_positive)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return sin


class acos(InverseTrigonometricFunction):
    """
    The inverse cosine function.

    Returns the arc cosine of x (measured in radians).

    Notes
    =====

    ``acos(x)`` will evaluate automatically in the cases
    ``oo``, ``-oo``, ``0``, ``1``, ``-1``.

    ``acos(zoo)`` evaluates to ``zoo``
    (see note in :py:class`diofant.functions.elementary.trigonometric.asec`)

    Examples
    ========

    >>> acos(1)
    0
    >>> acos(0)
    pi/2
    >>> acos(oo)
    oo*I

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Inverse_trigonometric_functions
    * https://dlmf.nist.gov/4.23
    * http://functions.wolfram.com/ElementaryFunctions/ArcCos

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1/sqrt(1 - self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg in (oo, -oo):
                return arg * I
            elif arg == 0:
                return pi / 2
            elif arg == 1:
                return Integer(0)
            elif arg == -1:
                return pi

        if arg is zoo:
            return zoo

        if arg.is_number:
            cst_table = {
                Rational(+1, 2): pi/3,
                Rational(-1, 2): 2*pi/3,
                sqrt(2)/2: pi/4,
                -sqrt(2)/2: 3*pi/4,
                1/sqrt(2): pi/4,
                -1/sqrt(2): 3*pi/4,
                sqrt(3)/2: pi/6,
                -sqrt(3)/2: 5*pi/6,
            }

            if arg in cst_table:
                return cst_table[arg]

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return pi / 2
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
                return -R / F * x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if (s.args[0] - 1).is_zero:
                return True
            elif s.args[0].is_rational and (s.args[0] - 1).is_nonzero:
                return False
        else:
            return s.is_rational

    def _eval_is_positive(self):
        x = self.args[0]
        return fuzzy_and([x.is_real, (1 - abs(x)).is_positive])

    def _eval_is_extended_real(self):
        x = self.args[0]
        return fuzzy_and([x.is_real, (1 - abs(x)).is_nonnegative])

    def _eval_nseries(self, x, n, logx):
        return self._eval_rewrite_as_log(self.args[0])._eval_nseries(x, n, logx)

    def _eval_rewrite_as_log(self, x):
        return pi/2 + I * \
            log(I * x + sqrt(1 - x**2))

    def _eval_rewrite_as_asin(self, x):
        return pi/2 - asin(x)

    def _eval_rewrite_as_atan(self, x):
        return atan(sqrt(1 - x**2)/x) + (pi/2)*(1 - x*sqrt(1/x**2))

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return cos

    def _eval_rewrite_as_acot(self, arg):
        return pi/2 - 2*acot((1 + sqrt(1 - arg**2))/arg)

    def _eval_rewrite_as_asec(self, arg):
        return asec(1/arg)

    def _eval_rewrite_as_acsc(self, arg):
        return pi/2 - acsc(1/arg)

    def _eval_conjugate(self):
        z = self.args[0]
        r = self.func(self.args[0].conjugate())
        if z.is_extended_real is False:
            return r
        elif z.is_extended_real and (z + 1).is_nonnegative and (z - 1).is_nonpositive:
            return r


class atan(InverseTrigonometricFunction):
    """
    The inverse tangent function.

    Returns the arc tangent of x (measured in radians).

    Notes
    =====

    atan(x) will evaluate automatically in the cases
    oo, -oo, 0, 1, -1.

    Examples
    ========

    >>> atan(0)
    0
    >>> atan(1)
    pi/4
    >>> atan(oo)
    pi/2

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Inverse_trigonometric_functions
    * https://dlmf.nist.gov/4.23
    * http://functions.wolfram.com/ElementaryFunctions/ArcTan

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(1 + self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg is oo:
                return pi / 2
            elif arg == -oo:
                return -pi / 2
            elif arg == 0:
                return Integer(0)
            elif arg == 1:
                return pi / 4
            elif arg == -1:
                return -pi / 4
        if arg.could_extract_minus_sign():
            return -cls(-arg)

        if arg.is_number:
            cst_table = {
                sqrt(3)/3: 6,
                -sqrt(3)/3: -6,
                1/sqrt(3): 6,
                -1/sqrt(3): -6,
                sqrt(3): 3,
                -sqrt(3): -3,
                (1 + sqrt(2)): Rational(8, 3),
                -(1 + sqrt(2)): Rational(8, 3),
                (sqrt(2) - 1): 8,
                (1 - sqrt(2)): -8,
                sqrt((5 + 2*sqrt(5))): Rational(5, 2),
                -sqrt((5 + 2*sqrt(5))): -Rational(5, 2),
                (2 - sqrt(3)): 12,
                -(2 - sqrt(3)): -12
            }

            if arg in cst_table:
                return pi / cst_table[arg]

        i_coeff = arg.as_coefficient(I)
        if i_coeff is not None:
            return I * atanh(i_coeff)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            return (-1)**((n - 1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_zero:
                return True
            elif s.args[0].is_rational and s.args[0].is_nonzero:
                return False
        else:
            return s.is_rational

    def _eval_is_positive(self):
        return self.args[0].is_positive

    def _eval_is_extended_real(self):
        return self.args[0].is_extended_real

    def _eval_rewrite_as_log(self, x):
        return I/2 * (log(
            (Integer(1) - I * x)/(Integer(1) + I * x)))

    def _eval_aseries(self, n, args0, x, logx):
        if args0[0] == oo:
            return (pi/2 - atan(1/self.args[0]))._eval_nseries(x, n, logx)
        else:
            return super()._eval_aseries(n, args0, x, logx)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return tan

    def _eval_rewrite_as_asin(self, arg):
        return sqrt(arg**2)/arg*(pi/2 - asin(1/sqrt(1 + arg**2)))

    def _eval_rewrite_as_acos(self, arg):
        return sqrt(arg**2)/arg*acos(1/sqrt(1 + arg**2))

    def _eval_rewrite_as_acot(self, arg):
        return acot(1/arg)

    def _eval_rewrite_as_asec(self, arg):
        return sqrt(arg**2)/arg*asec(sqrt(1 + arg**2))

    def _eval_rewrite_as_acsc(self, arg):
        return sqrt(arg**2)/arg*(pi/2 - acsc(sqrt(1 + arg**2)))


class acot(InverseTrigonometricFunction):
    """
    The inverse cotangent function.

    Returns the arc cotangent of x (measured in radians).  This
    function has a branch cut discontinuity in the complex plane
    running from `-i` to `i`.

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://dlmf.nist.gov/4.23
    * http://functions.wolfram.com/ElementaryFunctions/ArcCot
    * https://en.wikipedia.org/wiki/Inverse_trigonometric_functions

    """

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1 / (1 + self.args[0]**2)
        else:
            raise ArgumentIndexError(self, argindex)

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg in (oo, -oo):
                return Integer(0)
            elif arg == 0:
                return pi / 2
            elif arg == 1:
                return pi / 4
            elif arg == -1:
                return -pi / 4

        if arg.could_extract_minus_sign():
            return -cls(-arg)

        if arg.is_number:
            cst_table = {
                sqrt(3)/3: 3,
                -sqrt(3)/3: -3,
                1/sqrt(3): 3,
                -1/sqrt(3): -3,
                sqrt(3): 6,
                -sqrt(3): -6,
                (1 + sqrt(2)): 8,
                -(1 + sqrt(2)): -8,
                (1 - sqrt(2)): -Rational(8, 3),
                (sqrt(2) - 1): Rational(8, 3),
                sqrt(5 + 2*sqrt(5)): 10,
                -sqrt(5 + 2*sqrt(5)): -10,
                (2 + sqrt(3)): 12,
                -(2 + sqrt(3)): -12,
                (2 - sqrt(3)): Rational(12, 5),
                -(2 - sqrt(3)): -Rational(12, 5),
            }

            if arg in cst_table:
                return pi / cst_table[arg]

        i_coeff = arg.as_coefficient(I)
        if i_coeff is not None:
            return -I * acoth(i_coeff)

    @staticmethod
    @cacheit
    def taylor_term(n, x, *previous_terms):
        if n == 0:
            return pi / 2  # FIX THIS
        elif n < 0 or n % 2 == 0:
            return Integer(0)
        else:
            x = sympify(x)
            return (-1)**((n + 1)//2) * x**n / n

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)

        if x in arg.free_symbols and Order(1, x).contains(arg):
            return arg
        else:
            return self.func(arg)

    def _eval_is_rational(self):
        s = self.func(*self.args)
        if s.func == self.func:
            if s.args[0].is_rational:
                return False
        else:
            return s.is_rational

    def _eval_is_positive(self):
        return self.args[0].is_nonnegative

    def _eval_is_extended_real(self):
        return self.args[0].is_extended_real

    def _eval_aseries(self, n, args0, x, logx):
        if args0[0] == oo:
            return (pi/2 - acot(1/self.args[0]))._eval_nseries(x, n, logx)
        else:
            return super()._eval_aseries(n, args0, x, logx)

    def _eval_rewrite_as_log(self, x):
        return I/2 * \
            (log((x - I)/(x + I)))

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return cot

    def _eval_rewrite_as_asin(self, arg):
        return (arg*sqrt(1/arg**2) *
                (pi/2 - asin(sqrt(-arg**2)/sqrt(-arg**2 - 1))))

    def _eval_rewrite_as_acos(self, arg):
        return arg*sqrt(1/arg**2)*acos(sqrt(-arg**2)/sqrt(-arg**2 - 1))

    def _eval_rewrite_as_atan(self, arg):
        return atan(1/arg)

    def _eval_rewrite_as_asec(self, arg):
        return arg*sqrt(1/arg**2)*asec(sqrt((1 + arg**2)/arg**2))

    def _eval_rewrite_as_acsc(self, arg):
        return arg*sqrt(1/arg**2)*(pi/2 - acsc(sqrt((1 + arg**2)/arg**2)))


class asec(InverseTrigonometricFunction):
    r"""
    The inverse secant function.

    Returns the arc secant of x (measured in radians).

    Notes
    =====

    ``asec(x)`` will evaluate automatically in the cases
    ``oo``, ``-oo``, ``0``, ``1``, ``-1``.

    ``asec(x)`` has branch cut in the interval [-1, 1]. For complex arguments,
    it can be defined as

    .. math::
        sec^{-1}(z) = -i*(log(\sqrt{1 - z^2} + 1) / z)

    At ``x = 0``, for positive branch cut, the limit evaluates to ``zoo``. For
    negative branch cut, the limit

    .. math::
        \lim_{z \to 0}-i*(log(-\sqrt{1 - z^2} + 1) / z)

    simplifies to `-i*log(z/2 + O(z^3))` which ultimately evaluates to
    ``zoo``.

    As ``asex(x)`` = ``asec(1/x)``, a similar argument can be given for
    ``acos(x)``.

    Examples
    ========

    >>> asec(1)
    0
    >>> asec(-1)
    pi

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Inverse_trigonometric_functions
    * https://dlmf.nist.gov/4.23
    * http://functions.wolfram.com/ElementaryFunctions/ArcSec
    * https://reference.wolfram.com/language/ref/ArcSec.html

    """

    @classmethod
    def eval(cls, arg):
        if arg.is_zero:
            return zoo
        if arg.is_Number:
            if arg == 1:
                return Integer(0)
            elif arg == -1:
                return pi
        if arg in [oo, -oo, zoo]:
            return pi/2

    def fdiff(self, argindex=1):
        if argindex == 1:
            return 1/(self.args[0]**2*sqrt(1 - 1/self.args[0]**2))
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return sec

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)
        if Order(1, x).contains(arg):
            return log(arg)
        else:
            return self.func(arg)

    def _eval_rewrite_as_log(self, arg):
        return pi/2 + I*log(I/arg + sqrt(1 - 1/arg**2))

    def _eval_rewrite_as_asin(self, arg):
        return pi/2 - asin(1/arg)

    def _eval_rewrite_as_acos(self, arg):
        return acos(1/arg)

    def _eval_rewrite_as_atan(self, arg):
        return sqrt(arg**2)/arg*(-pi/2 + 2*atan(arg + sqrt(arg**2 - 1)))

    def _eval_rewrite_as_acot(self, arg):
        return sqrt(arg**2)/arg*(-pi/2 + 2*acot(arg - sqrt(arg**2 - 1)))

    def _eval_rewrite_as_acsc(self, arg):
        return pi/2 - acsc(arg)


class acsc(InverseTrigonometricFunction):
    """
    The inverse cosecant function.

    Returns the arc cosecant of x (measured in radians).

    Notes
    =====

    acsc(x) will evaluate automatically in the cases
    oo, -oo, 0, 1, -1.

    Examples
    ========

    >>> acsc(1)
    pi/2
    >>> acsc(-1)
    -pi/2

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot
    diofant.functions.elementary.trigonometric.atan2

    References
    ==========

    * https://en.wikipedia.org/wiki/Inverse_trigonometric_functions
    * https://dlmf.nist.gov/4.23
    * http://functions.wolfram.com/ElementaryFunctions/ArcCsc

    """

    @classmethod
    def eval(cls, arg):
        if arg.is_Number:
            if arg == 1:
                return pi/2
            elif arg == -1:
                return -pi/2
        if arg in [oo, -oo, zoo]:
            return Integer(0)

    def fdiff(self, argindex=1):
        if argindex == 1:
            return -1/(self.args[0]**2*sqrt(1 - 1/self.args[0]**2))
        else:
            raise ArgumentIndexError(self, argindex)

    def inverse(self, argindex=1):
        """Returns the inverse of this function."""
        return csc

    def _eval_as_leading_term(self, x):
        from ...series import Order
        arg = self.args[0].as_leading_term(x)
        if Order(1, x).contains(arg):
            return log(arg)
        else:
            return self.func(arg)

    def _eval_rewrite_as_log(self, arg):
        return -I*log(I/arg + sqrt(1 - 1/arg**2))

    def _eval_rewrite_as_asin(self, arg):
        return asin(1/arg)

    def _eval_rewrite_as_acos(self, arg):
        return pi/2 - acos(1/arg)

    def _eval_rewrite_as_atan(self, arg):
        return sqrt(arg**2)/arg*(pi/2 - atan(sqrt(arg**2 - 1)))

    def _eval_rewrite_as_acot(self, arg):
        return sqrt(arg**2)/arg*(pi/2 - acot(1/sqrt(arg**2 - 1)))

    def _eval_rewrite_as_asec(self, arg):
        return pi/2 - asec(arg)


class atan2(InverseTrigonometricFunction):
    r"""
    The function ``atan2(y, x)`` computes `\operatorname{atan}(y/x)` taking
    two arguments `y` and `x`.  Signs of both `y` and `x` are considered to
    determine the appropriate quadrant of `\operatorname{atan}(y/x)`.
    The range is `(-\pi, \pi]`. The complete definition reads as follows:

    .. math::

        \operatorname{atan2}(y, x) =
        \begin{cases}
          \arctan\left(\frac y x\right) & \qquad x > 0 \\
          \arctan\left(\frac y x\right) + \pi& \qquad y \ge 0 , x < 0 \\
          \arctan\left(\frac y x\right) - \pi& \qquad y < 0 , x < 0 \\
          +\frac{\pi}{2} & \qquad y > 0 , x = 0 \\
          -\frac{\pi}{2} & \qquad y < 0 , x = 0 \\
          \text{undefined} & \qquad y = 0, x = 0
        \end{cases}

    Attention: Note the role reversal of both arguments. The `y`-coordinate
    is the first argument and the `x`-coordinate the second.

    Examples
    ========

    Going counter-clock wise around the origin we find the
    following angles:

    >>> atan2(0, 1)
    0
    >>> atan2(1, 1)
    pi/4
    >>> atan2(1, 0)
    pi/2
    >>> atan2(1, -1)
    3*pi/4
    >>> atan2(0, -1)
    pi
    >>> atan2(-1, -1)
    -3*pi/4
    >>> atan2(-1, 0)
    -pi/2
    >>> atan2(-1, 1)
    -pi/4

    which are all correct. Compare this to the results of the ordinary
    `\operatorname{atan}` function for the point `(x, y) = (-1, 1)`

    >>> atan(Integer(1) / -1)
    -pi/4
    >>> atan2(1, -1)
    3*pi/4

    where only the `\operatorname{atan2}` function returns what we expect.
    We can differentiate the function with respect to both arguments:

    >>> diff(atan2(y, x), x)
    -y/(x**2 + y**2)

    >>> diff(atan2(y, x), y)
    x/(x**2 + y**2)

    We can express the `\operatorname{atan2}` function in terms of
    complex logarithms:

    >>> atan2(y, x).rewrite(log)
    -I*log((x + I*y)/sqrt(x**2 + y**2))

    and in terms of `\operatorname(atan)`:

    >>> atan2(y, x).rewrite(atan)
    2*atan(y/(x + sqrt(x**2 + y**2)))

    but note that this form is undefined on the negative real axis.

    See Also
    ========

    diofant.functions.elementary.trigonometric.sin
    diofant.functions.elementary.trigonometric.csc
    diofant.functions.elementary.trigonometric.cos
    diofant.functions.elementary.trigonometric.sec
    diofant.functions.elementary.trigonometric.tan
    diofant.functions.elementary.trigonometric.cot
    diofant.functions.elementary.trigonometric.asin
    diofant.functions.elementary.trigonometric.acsc
    diofant.functions.elementary.trigonometric.acos
    diofant.functions.elementary.trigonometric.asec
    diofant.functions.elementary.trigonometric.atan
    diofant.functions.elementary.trigonometric.acot

    References
    ==========

    * https://en.wikipedia.org/wiki/Inverse_trigonometric_functions
    * https://en.wikipedia.org/wiki/Atan2
    * http://functions.wolfram.com/ElementaryFunctions/ArcTan2

    """

    @classmethod
    def eval(cls, y, x):
        from .. import Heaviside
        from .complexes import im, re
        if x == -oo:
            if y.is_zero:
                # Special case y = 0 because we define Heaviside(0) = 1/2
                return pi
            return 2*pi*(Heaviside(re(y))) - pi
        elif x is oo:
            return Integer(0)
        elif x.is_imaginary and y.is_imaginary and x.is_number and y.is_number:
            x = im(x)
            y = im(y)

        if x.is_extended_real and y.is_extended_real:
            if x.is_positive:
                return atan(y / x)
            elif x.is_negative:
                if y.is_negative:
                    return atan(y / x) - pi
                elif y.is_nonnegative:
                    return atan(y / x) + pi
            elif x.is_zero:
                if y.is_positive:
                    return pi/2
                elif y.is_negative:
                    return -pi/2
                elif y.is_zero:
                    return nan
        if y.is_zero and x.is_extended_real and x.is_nonzero:
            return pi*(1 - Heaviside(x))
        if x.is_number and y.is_number:
            return -I*log((x + I*y)/sqrt(x**2 + y**2))

    def _eval_rewrite_as_log(self, y, x):
        return -I*log((x + I*y) / sqrt(x**2 + y**2))

    def _eval_rewrite_as_atan(self, y, x):
        return 2*atan(y / (sqrt(x**2 + y**2) + x))

    def _eval_rewrite_as_arg(self, y, x):
        from .complexes import arg
        if x.is_extended_real and y.is_extended_real:
            return arg(x + y*I)
        n = x + I*y
        d = x**2 + y**2
        return arg(n/sqrt(d)) - I*log(abs(n)/sqrt(abs(d)))

    def _eval_is_real(self):
        y, x = self.args
        if x.is_real and y.is_real:
            if x.is_nonzero and y.is_nonzero:
                return True

    def _eval_conjugate(self):
        return self.func(self.args[0].conjugate(), self.args[1].conjugate())

    def fdiff(self, argindex):
        y, x = self.args
        if argindex == 1:
            # Diff wrt y
            return x/(x**2 + y**2)
        elif argindex == 2:
            # Diff wrt x
            return -y/(x**2 + y**2)
        else:
            raise ArgumentIndexError(self, argindex)
