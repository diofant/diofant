import math

from mpmath.libmp import sqrtrem as mpmath_sqrtrem

from ..logic import true
from ..utilities import sift
from .add import Add
from .cache import cacheit
from .compatibility import as_int
from .evalf import PrecisionExhausted
from .evaluate import global_evaluate
from .expr import Expr
from .function import (_coeff_isneg, expand_complex, expand_mul,
                       expand_multinomial)
from .logic import fuzzy_or
from .mul import Mul, _keep_coeff
from .numbers import E, I, Integer, Rational, nan, oo, pi, zoo
from .symbol import Dummy, symbols
from .sympify import sympify


def integer_nthroot(y, n):
    """
    Return a tuple containing x = floor(y**(1/n))
    and a boolean indicating whether the result is exact (that is,
    whether x**n == y).

    >>> integer_nthroot(16, 2)
    (4, True)
    >>> integer_nthroot(26, 2)
    (5, False)

    """
    y, n = int(y), int(n)
    if y < 0:
        raise ValueError('y must be nonnegative')
    if n < 1:
        raise ValueError('n must be positive')
    if y in (0, 1):
        return y, True
    if n == 1:
        return y, True
    if n == 2:
        x, rem = mpmath_sqrtrem(y)
        return int(x), not rem
    if n > y:
        return 1, False
    # Get initial estimate for Newton's method. Care must be taken to
    # avoid overflow
    try:
        guess = int(y**(1./n) + 0.5)
    except OverflowError:
        exp = math.log(y, 2)/n
        if exp > 53:
            shift = int(exp - 53)
            guess = int(2.0**(exp - shift) + 1) << shift
        else:
            guess = int(2.0**exp)
    if guess > 2**50:
        # Newton iteration
        xprev, x = -1, guess
        while 1:
            t = x**(n - 1)
            xprev, x = x, ((n - 1)*x + y//t)//n
            if abs(x - xprev) < 2:
                break
    else:
        x = guess
    # Compensate
    t = x**n
    while t < y:
        x += 1
        t = x**n
    while t > y:
        x -= 1
        t = x**n
    return x, t == y


class Pow(Expr):
    """
    Defines the expression x**y as "x raised to a power y".

    For complex numbers `x` and `y`, ``Pow`` gives the principal
    value of `exp(y*log(x))`.

    Singleton definitions involving (0, 1, -1, oo, -oo, I, -I):

    +--------------+---------+-----------------------------------------------+
    | expr         | value   | reason                                        |
    +==============+=========+===============================================+
    | z**0         | 1       | Although arguments over 0**0 exist, see [2].  |
    +--------------+---------+-----------------------------------------------+
    | z**1         | z       |                                               |
    +--------------+---------+-----------------------------------------------+
    | (-oo)**(-1)  | 0       |                                               |
    +--------------+---------+-----------------------------------------------+
    | (-1)**-1     | -1      |                                               |
    +--------------+---------+-----------------------------------------------+
    | 0**-1        | zoo     | This is not strictly true, as 0**-1 may be    |
    |              |         | undefined, but is convenient in some contexts |
    |              |         | where the base is assumed to be positive.     |
    +--------------+---------+-----------------------------------------------+
    | 1**-1        | 1       |                                               |
    +--------------+---------+-----------------------------------------------+
    | oo**-1       | 0       |                                               |
    +--------------+---------+-----------------------------------------------+
    | 0**oo        | 0       | Because for all complex numbers z near        |
    |              |         | 0, z**oo -> 0.                                |
    +--------------+---------+-----------------------------------------------+
    | 0**-oo       | zoo     | This is not strictly true, as 0**oo may be    |
    |              |         | oscillating between positive and negative     |
    |              |         | values or rotating in the complex plane.      |
    |              |         | It is convenient, however, when the base      |
    |              |         | is positive.                                  |
    +--------------+---------+-----------------------------------------------+
    | 1**oo        | nan     | Because there are various cases where         |
    | 1**-oo       |         | lim(x(t),t)=1, lim(y(t),t)=oo (or -oo),       |
    |              |         | but lim(x(t)**y(t), t) != 1.  See [3].        |
    +--------------+---------+-----------------------------------------------+
    | z**zoo       | nan     | No limit for z**t for t -> zoo.               |
    +--------------+---------+-----------------------------------------------+
    | (-1)**oo     | nan     | Because of oscillations in the limit.         |
    | (-1)**(-oo)  |         |                                               |
    +--------------+---------+-----------------------------------------------+
    | oo**oo       | oo      |                                               |
    +--------------+---------+-----------------------------------------------+
    | oo**-oo      | 0       |                                               |
    +--------------+---------+-----------------------------------------------+
    | (-oo)**oo    | nan     |                                               |
    | (-oo)**-oo   |         |                                               |
    +--------------+---------+-----------------------------------------------+
    | oo**I        | nan     | oo**e could probably be best thought of as    |
    | (-oo)**I     |         | the limit of x**e for real x as x tends to    |
    |              |         | oo. If e is I, then the limit does not exist  |
    |              |         | and nan is used to indicate that.             |
    +--------------+---------+-----------------------------------------------+
    | oo**(1+I)    | zoo     | If the real part of e is positive, then the   |
    | (-oo)**(1+I) |         | limit of abs(x**e) is oo. So the limit value  |
    |              |         | is zoo.                                       |
    +--------------+---------+-----------------------------------------------+
    | oo**(-1+I)   | 0       | If the real part of e is negative, then the   |
    | -oo**(-1+I)  |         | limit is 0.                                   |
    +--------------+---------+-----------------------------------------------+

    Because symbolic computations are more flexible that floating point
    calculations and we prefer to never return an incorrect answer,
    we choose not to conform to all IEEE 754 conventions.  This helps
    us avoid extra test-case code in the calculation of limits.

    See Also
    ========

    diofant.core.numbers.Infinity
    diofant.core.numbers.NaN

    References
    ==========

    * https://en.wikipedia.org/wiki/Exponentiation
    * https://en.wikipedia.org/wiki/Zero_to_the_power_of_zero
    * https://en.wikipedia.org/wiki/Indeterminate_forms

    """

    is_Pow = True

    @cacheit
    def __new__(cls, b, e, evaluate=None):
        if evaluate is None:
            evaluate = global_evaluate[0]
        from ..functions.elementary.exponential import exp_polar

        b = sympify(b, strict=True)
        e = sympify(e, strict=True)
        if evaluate:
            if nan in (b, e):
                return nan
            if e is Integer(0):
                return Integer(1)
            if e is Integer(1):
                return b
            if e is zoo:
                return nan
            if e.is_integer and _coeff_isneg(b):
                if e.is_even:
                    b = -b
                elif e.is_odd:
                    return -Pow(-b, e)
            if b is Integer(1):
                if abs(e).is_infinite:
                    return nan
                return Integer(1)
            # recognize base as E
            if not e.is_Atom and b is not E and not isinstance(b, exp_polar):
                from ..functions import im, log, sign
                from ..simplify import denom, numer
                from .exprtools import factor_terms
                c, ex = factor_terms(e, sign=False).as_coeff_Mul()
                den = denom(ex)
                if isinstance(den, log) and den.args[0] == b:
                    return E**(c*numer(ex))
                if den.is_Add:
                    s = sign(im(b))
                    if s.is_Number and s and den == \
                            log(-factor_terms(b, sign=False)) + s*I*pi:
                        return E**(c*numer(ex))

            obj = b._eval_power(e)
            if obj is not None:
                return obj
        obj = Expr.__new__(cls, b, e)
        if b is E:
            obj.is_Exp = True
        return obj

    def _eval_is_commutative(self):
        return self.base.is_commutative and self.exp.is_commutative

    @property
    def base(self):
        """Returns base of the power expression."""
        return self.args[0]

    @property
    def exp(self):
        """Returns exponent of the power expression."""
        return self.args[1]

    @classmethod
    def class_key(cls):
        """Nice order of classes."""
        return 4, 2, cls.__name__

    def _eval_power(self, other):
        from ..functions import Abs, arg, exp, floor, im, log, re, sign
        b, e = self.as_base_exp()
        if b is nan:
            return (b**e)**other  # let __new__ handle it

        s = None
        if other.is_integer:
            s = 1
        elif b.is_polar:  # e.g. exp_polar, besselj, var('p', polar=True)...
            s = 1
        elif e.is_extended_real is not None:
            # helper functions ===========================
            def _half(e):
                """Return True if the exponent has a literal 2 as the
                denominator, else None.

                """
                if getattr(e, 'denominator', None) == 2:
                    return True
                n, d = e.as_numer_denom()
                if n.is_integer and d == 2:
                    return True

            def _n2(e):
                """Return ``e`` evaluated to a Number with 2 significant
                digits, else None.

                """
                try:
                    rv = e.evalf(2)
                    if rv.is_Number:
                        return rv
                except PrecisionExhausted:  # pragma: no cover
                    pass

            # ===================================================
            if e.is_extended_real:
                # we need _half(other) with constant floor or
                # floor(Rational(1, 2) - e*arg(b)/2/pi) == 0

                # handle -1 as special case
                if e == -1:
                    # floor arg. is 1/2 + arg(b)/2/pi
                    if _half(other):
                        if b.is_negative is True:
                            return (-1)**other*Pow(-b, e*other)
                        if b.is_extended_real is False:
                            return Pow(b.conjugate()/Abs(b)**2, other)
                elif e.is_even:
                    if b.is_extended_real:
                        b = abs(b)
                    if b.is_imaginary:
                        b = abs(im(b))*I

                if (abs(e) < 1) == true or (e == 1):
                    s = 1  # floor = 0
                elif b.is_nonnegative:
                    s = 1  # floor = 0
                elif re(b).is_nonnegative and (abs(e) < 2) == true:
                    s = 1  # floor = 0
                elif im(b).is_nonzero and (abs(e) == 2):
                    s = 1  # floor = 0
                elif b.is_imaginary and (abs(e) == 2):
                    s = 1  # floor = 0
                elif _half(other):
                    s = exp(2*pi*I*other*floor(
                        Rational(1, 2) - e*arg(b)/(2*pi)))
                    if s.is_extended_real and _n2(sign(s) - s) == 0:
                        s = sign(s)
                    else:
                        s = None
            else:
                # e.is_extended_real is False requires:
                #     _half(other) with constant floor or
                #     floor(Rational(1, 2) - im(e*log(b))/2/pi) == 0
                s = exp(2*I*pi*other*floor(Rational(1, 2) - im(e*log(b))/2/pi))
                # be careful to test that s is -1 or 1 b/c sign(I) == I:
                # so check that s is real
                if s.is_extended_real and _n2(sign(s) - s) == 0:
                    s = sign(s)
                else:
                    s = None

        if s is not None:
            return s*Pow(b, e*other)

    def _eval_is_positive(self):
        b, e = self.base, self.exp

        if b.is_nonnegative and b == e:
            return True
        if b.is_positive and (e.is_real or e.is_positive):
            return True
        if b.is_negative and e.is_integer and (b.is_finite or e.is_nonnegative):
            return e.is_even
        if b.is_nonpositive and e.is_odd and (b.is_finite or e.is_nonnegative):
            return False
        if b in {I, -I} and e.is_imaginary:
            return True

    def _eval_is_nonnegative(self):
        b, e = self.base, self.exp

        if b.is_imaginary and e.is_nonnegative:
            m = e % 4
            if m.is_integer:
                return m.is_zero

    def _eval_is_negative(self):
        b, e = self.base, self.exp

        if b.is_negative:
            if e.is_odd and (b.is_finite or e.is_positive):
                return True
            if e.is_even:
                return False
        elif b.is_positive:
            if e.is_extended_real:
                return False
        elif b.is_nonnegative:
            if e.is_nonnegative:
                return False
        elif b.is_nonpositive:
            if e.is_even:
                return False
        elif b.is_extended_real:
            if e.is_even:
                return False

    def _eval_is_zero(self):
        b, e = self.base, self.exp

        if b.is_zero:
            if e.is_positive:
                return True
            if e.is_nonpositive:
                return False
        elif b.is_nonzero:
            if e.is_finite:
                return False
            if e.is_infinite:
                if (1 - abs(b)).is_positive:
                    return e.is_positive
                if (1 - abs(b)).is_negative:
                    return e.is_negative

    def _eval_is_integer(self):
        b, e = self.base, self.exp

        if b.is_rational:
            if b.is_integer is False and e.is_positive:
                return False  # rat**nonneg
        if b.is_integer and e.is_integer:
            if b is Integer(-1):
                return True
            if e.is_nonnegative or e.is_positive:
                return True
        if b.is_integer and e.is_negative and (e.is_finite or e.is_integer):
            if (b - 1).is_nonzero and (b + 1).is_nonzero:
                return False
        if b.is_Number and e.is_Number:
            check = self.func(*self.args)
            if check.is_Integer:
                return True

    def _eval_is_extended_real(self):
        from ..functions import arg, log
        from .mul import Mul

        b, e = self.base, self.exp

        if b is E:
            if e.is_extended_real:
                return True
            if e.is_imaginary:
                return (2*I*e/pi).is_even

        if b.is_extended_real is None:
            if b.func == self.func and b.is_Exp and b.exp.is_imaginary:
                return e.is_imaginary
        if e.is_extended_real is None:
            return

        if b.is_extended_real and e.is_extended_real:
            if b.is_positive:
                return True
            if b.is_nonnegative:
                if e.is_nonnegative:
                    return True
            else:
                if e.is_integer:
                    if b.is_nonzero or e.is_nonnegative:
                        return True
                elif b.is_negative:
                    if e.is_rational and e.is_noninteger:
                        return False

        if b.is_nonzero and e.is_negative:
            return (b**-e).is_extended_real

        if b.is_imaginary:
            if e.is_integer:
                if e.is_even:
                    if b.is_nonzero or e.is_nonnegative:
                        return True
                elif e.is_odd:
                    return False
            elif e.is_imaginary and log(b).is_imaginary:
                return True
            elif e.is_Add:
                c, a = e.as_coeff_Add()
                if c and c.is_Integer:
                    return Mul(b**c, b**a, evaluate=False).is_extended_real
            elif b in (-I, I) and (e/2).is_noninteger:
                return False
            return

        if b.is_extended_real and e.is_imaginary:
            if b is Integer(-1):
                return True
            c = e.coeff(I)
            if c in (1, -1):
                if b == 2:
                    return False

        if b.is_extended_real is False:  # we already know it's not imag
            i = arg(b)*e/pi
            return i.is_integer

    def _eval_is_complex(self):
        from ..functions import log

        b, e = self.base, self.exp

        if b.is_complex:
            exp = log(b)*e
            return fuzzy_or([exp.is_complex, exp.is_negative])

    def _eval_is_imaginary(self):
        from ..functions import arg, log

        b, e = self.base, self.exp

        if b.is_imaginary:
            if e.is_integer:
                return e.is_odd

        if e.is_imaginary and e.is_nonzero:
            if log(b).is_imaginary:
                return False

        if b.is_real and e.is_real:
            if b.is_positive:
                return False
            if e.is_integer:
                return False
            if (2*e).is_integer:
                return b.is_negative

        if b.is_real is False:  # we already know it's not imag
            return (2*arg(b)*e/pi).is_odd

    def _eval_is_odd(self):
        b, e = self.base, self.exp

        if e.is_integer:
            if e.is_positive:
                return b.is_odd
            if e.is_nonnegative and b.is_odd:
                return True
            if b is Integer(-1):
                return True

    def _eval_is_finite(self):
        b, e = self.base, self.exp

        if e.is_negative:
            if b.is_zero:
                return False
        if b.is_finite and e.is_finite:
            if e.is_nonnegative or b.is_nonzero:
                return True

    def _eval_is_polar(self):
        b, e = self.base, self.exp

        if b.is_polar and e.is_commutative:
            return True

    def _eval_subs(self, old, new):
        from ..functions import log
        from .symbol import Symbol

        def _check(ct1, ct2, old):
            """Return bool, pow where, if bool is True, then the exponent of
            Pow `old` will combine with `pow` so the substitution is valid,
            otherwise bool will be False,

            cti are the coefficient and terms of an exponent of self or old
            In this _eval_subs routine a change like (b**(2*x)).subs({b**x: y})
            will give y**2 since (b**x)**2 == b**(2*x); if that equality does
            not hold then the substitution should not occur so `bool` will be
            False.

            """
            coeff1, terms1 = ct1
            coeff2, terms2 = ct2
            if terms1 == terms2:
                pow = coeff1/coeff2
                try:
                    pow = as_int(pow)
                    combines = True
                except ValueError:
                    combines = Pow._eval_power(
                        Pow(*old.as_base_exp(), evaluate=False),
                        pow)
                    if isinstance(combines, Pow):
                        combines = combines.base is old.base
                    else:
                        combines = False
                return combines, pow
            return False, None

        if old == self.base:
            return new**self.exp._subs(old, new)

        if old.func is self.func and self.exp == old.exp:
            l = log(self.base, old.base)
            if l.is_Number:
                return Pow(new, l)

        if isinstance(old, self.func) and self.base == old.base:
            if self.exp.is_Add is False:
                ct2 = old.exp.as_independent(Symbol, as_Add=False)
                ct1 = (self.exp/ct2[1], ct2[1])
                ok, pow = _check(ct1, ct2, old)
                if ok:
                    # issue sympy/sympy#5180: (x**(6*y)).subs({x**(3*y):z})->z**2
                    return self.func(new, pow)
            else:  # b**(6*x+a).subs({b**(3*x): y}) -> y**2 * b**a
                # exp(exp(x) + exp(x**2)).subs({exp(exp(x)): w}) -> w * exp(exp(x**2))
                oarg = old.exp
                new_l = []
                o_al = []
                ct2 = oarg.as_coeff_mul()
                for a in self.exp.args:
                    newa = a._subs(old, new)
                    ct1 = newa.as_coeff_mul()
                    ok, pow = _check(ct1, ct2, old)
                    if ok:
                        new_l.append(new**pow)
                        continue
                    o_al.append(newa)
                if new_l:
                    new_l.append(Pow(self.base, Add(*o_al), evaluate=False))
                    return Mul(*new_l)

        if old.is_Exp and self.exp.is_extended_real and self.base.is_positive:
            ct1 = old.exp.as_independent(Symbol, as_Add=False)
            ct2 = (self.exp*log(self.base)).as_independent(
                Symbol, as_Add=False)
            ok, pow = _check(ct1, ct2, old)
            if ok:
                return self.func(new, pow)  # (2**x).subs({exp(x*log(2)): z}) -> z

    def as_base_exp(self):
        """Return base and exp of self.

        If base is 1/Integer, then return Integer, -exp. If this extra
        processing is not needed, the base and exp properties will
        give the raw arguments

        Examples
        ========

        >>> p = Pow(Rational(1, 2), 2, evaluate=False)
        >>> p.as_base_exp()
        (2, -2)
        >>> p.args
        (1/2, 2)

        """
        b, e = self.base, self.exp
        if b.is_Rational and b.numerator == 1 and b.denominator != 1:
            return Integer(b.denominator), -e
        return b, e

    def _eval_adjoint(self):
        from ..functions.elementary.complexes import adjoint
        i, p = self.exp.is_integer, self.base.is_positive
        if i:
            return adjoint(self.base)**self.exp
        if p:
            return self.base**adjoint(self.exp)

    def _eval_conjugate(self):
        if self.is_extended_real:
            return self
        from ..functions.elementary.complexes import conjugate as c
        i, p = self.exp.is_integer, self.base.is_positive
        if i:
            return c(self.base)**self.exp
        if p:
            return self.base**c(self.exp)
        if i is False and p is False:
            expanded = expand_complex(self)
            assert expanded != self
            return c(expanded)

    def _eval_transpose(self):
        from ..functions.elementary.complexes import transpose
        i, p = self.exp.is_integer, self.base.is_complex
        if p:
            return self.base**self.exp
        if i:
            return transpose(self.base)**self.exp

    def _eval_expand_power_exp(self, **hints):
        """a**(n+m) -> a**n*a**m."""
        b = self.base
        e = self.exp
        if e.is_Add and e.is_commutative:
            expr = []
            for x in e.args:
                expr.append(self.func(self.base, x))
            return Mul(*expr)
        return self.func(b, e)

    def _eval_expand_power_base(self, **hints):
        """(a*b)**n -> a**n * b**n."""
        force = hints.get('force', False)

        b = self.base
        e = self.exp
        if not b.is_Mul:
            return self

        cargs, nc = b.args_cnc(split_1=False)

        # expand each term - this is top-level-only
        # expansion but we have to watch out for things
        # that don't have an _eval_expand method
        if nc:
            nc = [i._eval_expand_power_base(**hints)
                  if hasattr(i, '_eval_expand_power_base') else i
                  for i in nc]

            if e.is_Integer:
                if e.is_positive:
                    rv = Mul(*nc*e)
                else:
                    rv = 1/Mul(*nc*-e)
                if cargs:
                    rv *= Mul(*cargs)**e
                return rv

            if not cargs:
                return self.func(Mul(*nc), e, evaluate=False)

            nc = [Mul(*nc)]

        # sift the commutative bases
        def pred(x):
            if x is I:
                return I
            polar = x.is_polar
            if polar:
                return True
            if polar is None:
                return fuzzy_or([x.is_nonnegative, (1/x).is_nonnegative])
        sifted = sift(cargs, pred)
        nonneg = sifted[True]
        other = sifted[None]
        neg = sifted[False]
        imag = sifted[I]
        if imag:
            i = len(imag) % 4
            if i == 0:
                pass
            elif i == 1:
                other.append(I)
            elif i == 2:
                if neg:
                    nonn = -neg.pop()
                    if nonn is not Integer(1):
                        nonneg.append(nonn)
                else:
                    neg.append(Integer(-1))
            else:
                if neg:
                    nonn = -neg.pop()
                    if nonn is not Integer(1):
                        nonneg.append(nonn)
                else:
                    neg.append(Integer(-1))
                other.append(I)
            del imag

        # bring out the bases that can be separated from the base

        if force or e.is_integer:
            # treat all commutatives the same and put nc in other
            cargs = nonneg + neg + other
            other = nc
        else:
            # this is just like what is happening automatically, except
            # that now we are doing it for an arbitrary exponent for which
            # no automatic expansion is done

            assert not e.is_Integer

            # handle negatives by making them all positive and putting
            # the residual -1 in other
            if len(neg) > 1:
                o = Integer(1)
                if not other and neg[0].is_Number:
                    o *= neg.pop(0)
                if len(neg) % 2:
                    o = -o
                for n in neg:
                    nonneg.append(-n)
                if o is not Integer(1):
                    other.append(o)
            elif neg and other:
                if neg[0].is_Number and neg[0] is not Integer(-1):
                    other.append(Integer(-1))
                    nonneg.append(-neg[0])
                else:
                    other.extend(neg)
            else:
                other.extend(neg)
            del neg

            cargs = nonneg
            other += nc

        rv = Integer(1)
        if cargs:
            rv *= Mul(*[self.func(b, e, evaluate=False) for b in cargs])
        if other:
            rv *= self.func(Mul(*other), e, evaluate=False)
        return rv

    def _eval_expand_multinomial(self, **hints):
        """(a+b+..) ** n -> a**n + n*a**(n-1)*b + .., n is nonzero integer."""
        base, exp = self.base, self.exp
        result = self

        if exp.is_Rational and exp.numerator > 0 and base.is_Add:
            if not exp.is_Integer:
                n = Integer(exp.numerator // exp.denominator)

                if not n:
                    return result
                radical, result = self.func(base, exp - n), []

                expanded_base_n = self.func(base, n)
                if expanded_base_n.is_Pow:
                    expanded_base_n = \
                        expanded_base_n._eval_expand_multinomial()
                for term in Add.make_args(expanded_base_n):
                    result.append(term*radical)

                return Add(*result)

            n = int(exp)

            if base.is_commutative:
                order_terms, other_terms = [], []

                for b in base.args:
                    if b.is_Order:
                        order_terms.append(b)
                    else:
                        other_terms.append(b)

                if order_terms:
                    # (f(x) + O(x^n))^m -> f(x)^m + m*f(x)^{m-1} *O(x^n)
                    f = Add(*other_terms)
                    o = Add(*order_terms)

                    if n == 2:
                        return expand_multinomial(f**n, deep=False) + n*f*o
                    g = expand_multinomial(f**(n - 1), deep=False)
                    return expand_mul(f*g, deep=False) + n*g*o

                if base.is_number:
                    # Efficiently expand expressions of the form (a + b*I)**n
                    # where 'a' and 'b' are real numbers and 'n' is integer.
                    a, b = base.as_real_imag()

                    if a.is_Rational and b.is_Rational:
                        if not a.is_Integer:
                            if not b.is_Integer:
                                k = self.func(a.denominator * b.denominator, n)
                                a, b = a.numerator*b.denominator, a.denominator*b.numerator
                            else:
                                k = self.func(a.denominator, n)
                                a, b = a.numerator, a.denominator*b
                        elif not b.is_Integer:
                            k = self.func(b.denominator, n)
                            a, b = a*b.denominator, b.numerator
                        else:
                            k = 1

                        a, b, c, d = int(a), int(b), 1, 0

                        while n:
                            if n & 1:
                                c, d = a*c - b*d, b*c + a*d
                                n -= 1
                            a, b = a*a - b*b, 2*a*b
                            n //= 2

                        if k == 1:
                            return c + I*d
                        return Integer(c)/k + I*d/k

                p = other_terms
                # (x+y)**3 -> x**3 + 3*x**2*y + 3*x*y**2 + y**3
                # in this particular example:
                # p = [x,y]; n = 3
                # so now it's easy to get the correct result -- we get the
                # coefficients first:
                from ..ntheory import multinomial_coefficients
                from ..polys import Poly
                expansion_dict = multinomial_coefficients(len(p), n)
                # in our example: {(3, 0): 1, (1, 2): 3, (0, 3): 1, (2, 1): 3}
                # and now construct the expression.
                return Poly(expansion_dict, *p).as_expr()
            if n == 2:
                return Add(*[f*g for f in base.args for g in base.args])
            multi = (base**(n - 1))._eval_expand_multinomial()
            assert multi.is_Add
            return Add(*[f*g for f in base.args for g in multi.args])
        if (exp.is_Rational and exp.numerator < 0 and base.is_Add and
                abs(exp.numerator) > exp.denominator):
            return 1 / self.func(base, -exp)._eval_expand_multinomial()
        if exp.is_Add and base.is_Number:
            #  a + b      a  b
            # n      --> n  n  , where n, a, b are Numbers

            coeff, tail = Integer(1), Integer(0)
            for term in exp.args:
                if term.is_Number:
                    coeff *= self.func(base, term)
                else:
                    tail += term

            return coeff * self.func(base, tail)
        return result

    def as_real_imag(self, deep=True, **hints):
        """Returns real and imaginary parts of self

        See Also
        ========

        diofant.core.expr.Expr.as_real_imag

        """
        from ..functions import arg, cos, sin

        if self.exp.is_Integer:
            exp = self.exp
            re, im = self.base.as_real_imag(deep=deep)
            if not im:
                return self, Integer(0)
            a, b = symbols('a b', cls=Dummy)
            if exp >= 0:
                if re.is_Number and im.is_Number:
                    # We can be more efficient in this case
                    expr = expand_multinomial(self.base**exp)
                    return expr.as_real_imag()

                expr = ((a + b)**exp).as_poly()  # a = re, b = im; expr = (a + b*I)**exp
            else:
                mag = re**2 + im**2
                re, im = re/mag, -im/mag
                if re.is_Number and im.is_Number:
                    # We can be more efficient in this case
                    expr = expand_multinomial((re + im*I)**-exp)
                    return expr.as_real_imag()

                expr = ((a + b)**-exp).as_poly()

            # Terms with even b powers will be real
            r = [i for i in expr.terms() if not i[0][1] % 2]
            re_part = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])
            # Terms with odd b powers will be imaginary
            r = [i for i in expr.terms() if i[0][1] % 4 == 1]
            im_part1 = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])
            r = [i for i in expr.terms() if i[0][1] % 4 == 3]
            im_part3 = Add(*[cc*a**aa*b**bb for (aa, bb), cc in r])

            return (re_part.subs({a: re, b: I*im}),
                    im_part1.subs({a: re, b: im}) + im_part3.subs({a: re, b: -im}))

        if self.exp.is_Rational:
            re, im = self.base.as_real_imag(deep=deep)

            if im.is_zero and self.exp is Rational(1, 2):
                if re.is_nonnegative:
                    return self, Integer(0)
                if re.is_nonpositive:
                    return Integer(0), (-self.base)**self.exp

            # XXX: This is not totally correct since for x**(p/q) with
            #      x being imaginary there are actually q roots, but
            #      only a single one is returned from here.
            r = self.func(self.func(re, 2) + self.func(im, 2), Rational(1, 2))
            t = arg(re + I*im)

            rp, tp = self.func(r, self.exp), t*self.exp

            return rp*cos(tp), rp*sin(tp)
        if self.is_Exp:
            from ..functions import exp
            re, im = self.exp.as_real_imag()
            if deep:
                re = re.expand(deep, **hints)
                im = im.expand(deep, **hints)
            c, s = cos(im), sin(im)
            return exp(re)*c, exp(re)*s
        from ..functions import im, re
        if deep:
            hints['complex'] = False

            expanded = self.expand(deep, **hints)
            if hints.get('ignore') != expanded:
                return re(expanded), im(expanded)
        else:
            return re(self), im(self)

    def _eval_derivative(self, s):
        from ..functions import log
        dbase = self.base.diff(s)
        dexp = self.exp.diff(s)
        return self * (dexp * log(self.base) + dbase * self.exp/self.base)

    def _eval_evalf(self, prec):
        base, exp = self.as_base_exp()
        base = base._evalf(prec)
        if not exp.is_Integer:
            exp = exp._evalf(prec)
        return self.func(base, exp)

    def _eval_is_polynomial(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return bool(self.base._eval_is_polynomial(syms) and
                        self.exp.is_Integer and (self.exp >= 0))
        return True

    def _eval_is_rational(self):
        p = self.func(*self.as_base_exp())  # in case it's unevaluated
        if not p.is_Pow:
            return p.is_rational
        b, e = p.base, p.exp

        if e.is_Rational and b.is_Rational:
            # we didn't check that e is not an Integer
            # because Rational**Integer autosimplifies
            return False
        if e.is_integer:
            if b.is_rational:
                if b.is_nonzero or e.is_nonnegative:
                    return True
                if b == e:  # always rational, even for 0**0
                    return True
            elif b.is_irrational:
                if e.is_zero:
                    return True
        if b is E:
            if e.is_rational and e.is_nonzero:
                return False

    def _eval_is_algebraic(self):
        b, e = self.base, self.exp

        if b.is_zero or (b - 1).is_zero:
            return True
        if b is E:
            s = self.doit()
            if s.func == self.func:
                if e.is_nonzero:
                    if e.is_algebraic:
                        return False
                    if (e/pi).is_rational:
                        return False
                    if (e/(I*pi)).is_rational:
                        return True
            else:
                return s.is_algebraic
        elif e.is_rational and e.is_nonzero:
            if b.is_nonzero or e.is_nonnegative:
                return b.is_algebraic
        elif b.is_algebraic and e.is_algebraic:
            if (b.is_nonzero and (b - 1).is_nonzero) or b.is_irrational:
                return e.is_rational

    def _eval_is_rational_function(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return self.base._eval_is_rational_function(syms) and \
                self.exp.is_Integer
        return True

    def _eval_is_algebraic_expr(self, syms):
        if self.exp.has(*syms):
            return False

        if self.base.has(*syms):
            return self.base._eval_is_algebraic_expr(syms) and \
                self.exp.is_Rational
        return True

    def _eval_as_numer_denom(self):
        """Expression -> a/b -> a, b.

        See Also
        ========

        diofant.core.expr.Expr.as_numer_denom

        """
        if not self.is_commutative:
            return self, Integer(1)
        base, exp = self.as_base_exp()
        if base is Integer(1):
            return self, Integer(1)
        n, d = base.as_numer_denom()
        neg_exp = exp.is_negative
        if not neg_exp and not (-exp).is_negative:
            neg_exp = _coeff_isneg(exp)
        int_exp = exp.is_integer
        # the denominator cannot be separated from the numerator if
        # its sign is unknown unless the exponent is an integer, e.g.
        # sqrt(a/b) != sqrt(a)/sqrt(b) when a=1 and b=-1. But if the
        # denominator is negative the numerator and denominator can
        # be negated and the denominator (now positive) separated.
        if not (d.is_extended_real or int_exp):
            n = base
            d = Integer(1)
        dnonpos = d.is_nonpositive
        if dnonpos:
            n, d = -n, -d
        elif dnonpos is None and not int_exp:
            n = base
            d = Integer(1)
        if neg_exp:
            n, d = d, n
            exp = -exp
        if d is Integer(1):
            return self.func(n, exp), Integer(1)
        if n is Integer(1):
            return Integer(1), self.func(d, exp)
        return self.func(n, exp), self.func(d, exp)

    def _matches(self, expr, repl_dict={}):
        """Helper method for match().

        See Also
        ========

        diofant.core.basic.Basic.matches

        """
        expr = sympify(expr, strict=True)

        # special case, pattern = 1 and expr.exp can match to 0
        if expr is Integer(1):
            d = repl_dict.copy()
            d = self.exp._matches(Integer(0), d)
            if d is not None:
                return d

        # make sure the expression to be matched is an Expr
        if not isinstance(expr, Expr):
            return

        b, e = expr.as_base_exp()

        # special case number
        sb, se = self.as_base_exp()
        if sb.is_Symbol and se.is_Integer and expr:
            if e.is_rational:
                return sb._matches(b**(e/se), repl_dict)
            return sb._matches(expr**(1/se), repl_dict)

        d = repl_dict.copy()
        d = self.base._matches(b, d)
        if d is None:
            return

        d = self.exp.xreplace(d)._matches(e, d)
        if d is None:
            return Expr._matches(self, expr, repl_dict)
        return d

    def _eval_nseries(self, x, n, logx):
        from ..calculus import Order, limit
        from ..functions import arg, exp, floor, log
        from ..simplify import powsimp
        if self.is_Exp:
            e_series = self.exp.nseries(x, n=n, logx=logx)
            if e_series.is_Order:
                return 1 + e_series
            e0 = limit(e_series.removeO(), x, 0)
            if e0 in (-oo, oo):
                return self
            t = e_series - e0
            exp_series = term = exp(e0)
            # series of exp(e0 + t) in t
            for i in range(1, n):
                term *= t/i
                term = term.nseries(x, n=n, logx=logx)
                exp_series += term
            exp_series += Order(t**n, x)
            return powsimp(exp_series, deep=True, combine='exp')
        if self.exp.has(x):
            return exp(self.exp*log(self.base)).nseries(x, n=n, logx=logx)

        b_series = self.base.nseries(x, n=n, logx=logx)
        while b_series.is_Order:
            n += 1
            b_series = self.base.nseries(x, n=n, logx=logx)
        b0 = b_series.as_leading_term(x)
        t = expand_mul(expand_multinomial(b_series/b0 - 1).cancel())
        if t.is_Add:
            t = t.func(*[i for i in t.args if i.limit(x, 0).is_finite])
        c, e = b0.as_coeff_exponent(x)
        if self.exp is oo:
            if e != 0:
                sig = -e
            else:
                sig = abs(c) - 1 if c != 1 else t.removeO()
            if sig.is_positive:
                return oo
            if sig.is_negative:
                return Integer(0)
            raise NotImplementedError
        pow_series = term = Integer(1)
        # series of (1 + t)**e in t
        for i in range(1, n):
            term *= (self.exp - i + 1)*t/i
            term = term.nseries(x, n=n, logx=logx)
            pow_series += term
        factor = b0**self.exp
        if t != 0 and not (self.exp.is_Integer and self.exp >= 0 and n > self.exp):
            pow_series += Order(t**n, x)
            # branch handling
            if c.is_negative:
                l = floor(arg(t.removeO()*c)/(2*pi)).limit(x, 0)
                assert l.is_finite
                factor *= exp(2*pi*I*self.exp*l)
        pow_series = expand_mul(factor*pow_series)
        return powsimp(pow_series, deep=True, combine='exp')

    def _eval_as_leading_term(self, x):
        from ..calculus import Order
        from ..functions import exp, log
        if not self.exp.has(x):
            return self.func(self.base.as_leading_term(x), self.exp)
        if self.is_Exp:
            if self.exp.is_Mul:
                k, arg = self.exp.as_independent(x)
            else:
                k, arg = Integer(1), self.exp
            if arg.is_Add:
                return Mul(*[exp(k*f).as_leading_term(x) for f in arg.args])
            arg = self.exp.as_leading_term(x)
            if Order(1, x).contains(arg):
                return Integer(1)
            return exp(arg)
        return exp(self.exp*log(self.base)).as_leading_term(x)

    def _eval_rewrite_as_sin(self, base, exp):
        from ..functions import sin
        if self.is_Exp:
            return sin(I*self.exp + pi/2) - I*sin(I*self.exp)

    def _eval_rewrite_as_cos(self, base, exp):
        from ..functions import cos
        if self.is_Exp:
            return cos(I*self.exp) + I*cos(I*self.exp + pi/2)

    def _eval_rewrite_as_tanh(self, base, exp):
        from ..functions import tanh
        if self.is_Exp:
            return (1 + tanh(self.exp/2))/(1 - tanh(self.exp/2))

    def as_content_primitive(self, radical=False):
        """Return the tuple (R, self/R) where R is the positive Rational
        extracted from self.

        Examples
        ========

        >>> sqrt(4 + 4*sqrt(2)).as_content_primitive()
        (2, sqrt(1 + sqrt(2)))
        >>> sqrt(3 + 3*sqrt(2)).as_content_primitive()
        (1, sqrt(3)*sqrt(1 + sqrt(2)))

        >>> ((2*x + 2)**2).as_content_primitive()
        (4, (x + 1)**2)
        >>> (4**((1 + y)/2)).as_content_primitive()
        (2, 4**(y/2))
        >>> (3**((1 + y)/2)).as_content_primitive()
        (1, 3**((y + 1)/2))
        >>> (3**((5 + y)/2)).as_content_primitive()
        (9, 3**((y + 1)/2))
        >>> eq = 3**(2 + 2*x)
        >>> powsimp(eq) == eq
        True
        >>> eq.as_content_primitive()
        (9, 3**(2*x))
        >>> powsimp(Mul(*_))
        3**(2*x + 2)

        >>> eq = (2 + 2*x)**y
        >>> s = expand_power_base(eq)
        >>> s.is_Mul, s
        (False, (2*x + 2)**y)
        >>> eq.as_content_primitive()
        (1, (2*(x + 1))**y)
        >>> s = expand_power_base(_[1])
        >>> s.is_Mul, s
        (True, 2**y*(x + 1)**y)

        See Also
        ========

        diofant.core.expr.Expr.as_content_primitive

        """
        b, e = self.as_base_exp()
        b = _keep_coeff(*b.as_content_primitive(radical=radical))
        ce, pe = e.as_content_primitive(radical=radical)
        if b.is_Rational:
            # e
            # = ce*pe
            # = ce*(h + t)
            # = ce*h + ce*t
            # => self
            # = b**(ce*h)*b**(ce*t)
            # = b**(cehp/cehq)*b**(ce*t)
            # = b**(iceh+r/cehq)*b**(ce*t)
            # = b**iceh*b**(r/cehq)*b**(ce*t)
            # = b**iceh*b**(ce*t + r/cehq)
            h, t = pe.as_coeff_Add()
            if h.is_Rational:
                ceh = ce*h
                c = self.func(b, ceh)
                r = Integer(0)
                if not c.is_Rational:
                    iceh, r = divmod(ceh.numerator, ceh.denominator)
                    c = self.func(b, iceh)
                return c, self.func(b, _keep_coeff(ce, t + r/ce/ceh.denominator))
        e = _keep_coeff(ce, pe)
        # b**e = (h*t)**e = h**e*t**e = c*m*t**e
        if e.is_Rational and b.is_Mul:
            h, t = b.as_content_primitive(radical=radical)  # h is positive
            c, m = self.func(h, e).as_coeff_Mul()  # so c is positive
            m, me = m.as_base_exp()
            if m is Integer(1) or me == e:  # probably always true
                # return the following, not return c, m*Pow(t, e)
                # which would change Pow into Mul; we let diofant
                # decide what to do by using the unevaluated Mul, e.g
                # should it stay as sqrt(2 + 2*sqrt(5)) or become
                # sqrt(2)*sqrt(1 + sqrt(5))
                return c, self.func(_keep_coeff(m, t), e)
        return Integer(1), self.func(b, e)
