import decimal
import fractions

import mpmath
import pytest
from mpmath.libmp.libmpf import finf, fninf

from diofant import (Catalan, E, EulerGamma, Float, Ge, GoldenRatio, Gt, I,
                     Integer, Le, Lt, Mul, Number, Pow, Rational, Symbol, cbrt,
                     comp, cos, exp, factorial, false, integer_digits,
                     integer_nthroot, latex, log, mod_inverse, nan, nextprime,
                     oo, pi, root, sin, sqrt, true, zoo)
from diofant.abc import x
from diofant.core.cache import clear_cache
from diofant.core.numbers import igcdex, mpf_norm


__all__ = ()


def same_and_same_prec(a, b):
    # stricter matching for Floats
    return a == b and a._prec == b._prec


def test_mod():
    x = Rational(1, 2)
    y = Rational(3, 4)
    z = Rational(5, 18043)

    assert x % x == 0
    assert x % y == Rational(1, 2)
    assert x % z == Rational(3, 36086)
    assert y % x == Rational(1, 4)
    assert y % y == 0
    assert y % z == Rational(9, 72172)
    assert z % x == Rational(5, 18043)
    assert z % y == Rational(5, 18043)
    assert z % z == 0

    a = Float(2.6)

    assert (a % .2) == 0
    assert (a % 2).round(15) == 0.6
    assert (a % 0.5).round(15) == 0.1

    a = Rational(7, 2)
    assert (a % pi) == a - pi

    p = Symbol('p', infinite=True)

    assert zoo % 0 == nan
    assert oo % oo == nan
    assert zoo % oo == nan
    assert 5 % oo == nan
    assert p % 5 == nan

    # In these two tests, if the precision of m does
    # not match the precision of the ans, then it is
    # likely that the change made now gives an answer
    # with degraded accuracy.
    r = Rational(500, 41)
    f = Float('.36', 3)
    m = r % f
    ans = Float(r % Rational(f), 3)
    assert m == ans
    assert m._prec == ans._prec
    f = Float('8.36', 3)
    m = f % r
    ans = Float(Rational(f) % r, 3)
    assert m == ans
    assert m._prec == ans._prec

    assert Integer(0) % float(1) == 0

    # No rounding required since these numbers can be represented
    # exactly.
    assert Rational(3, 4) % Float(1.1) == 0.75
    assert Float(1.5) % Rational(5, 4) == 0.25
    assert Float('2.75') % Float('1.5') == Float('1.25')
    assert 2.75 % Float('1.5') == Float('1.25')

    a = Integer(7)
    b = Integer(4)

    assert type(a % b) == Integer
    assert a % b == 3
    assert Integer(1) % Rational(2, 3) == Rational(1, 3)
    assert Rational(7, 5) % Integer(1) == Rational(2, 5)
    assert Integer(2) % 1.5 == 0.5

    assert Integer(10) % Integer(3) == Integer(1)
    assert Integer(10) % 4 == 2
    assert 15 % Integer(4) == 3

    assert (Float(1) % zoo) is nan


def test_divmod():
    assert divmod(Integer(12), Integer(8)) == (1, 4)
    assert divmod(-Integer(12), Integer(8)) == (-2, 4)
    assert divmod(Integer(0), Integer(1)) == (0, 0)
    pytest.raises(ZeroDivisionError, lambda: divmod(Integer(0), Integer(0)))
    pytest.raises(ZeroDivisionError, lambda: divmod(Integer(1), Integer(0)))
    assert divmod(Integer(12), 8) == (1, 4)
    assert divmod(12, Integer(8)) == (1, 4)

    assert divmod(Integer(2), Rational(3, 2)) == (1, Rational(1, 2))
    assert divmod(Rational(3, 2), Integer(2)) == (0, Rational(3, 2))
    assert divmod(Integer(2), Float(3.5)) == (0, 2)
    assert divmod(Float(3.5), Integer(2)) == (1, Float(1.5))
    assert divmod(Integer(2), Rational(1, 3)) == (6, 0)
    assert divmod(Rational(1, 3), Integer(2)) == (0, Rational(1, 3))
    assert divmod(Integer(2), Float(0.1)) == (20, 0)
    assert divmod(Float(0.1), Integer(2)) == (0, Float(0.1))
    assert divmod(Integer(2), 2) == (1, 0)
    assert divmod(2, Integer(2)) == (1, 0)
    assert divmod(Integer(2), 1.5) == (1, Float(0.5))
    assert divmod(1.5, Integer(2)) == (0, Float(1.5))
    assert divmod(0.3, Integer(2)) == (0, Float(0.3))
    assert divmod(Rational(3, 2), Float(3.5)) == (0, Rational(3, 2))
    assert divmod(Float(3.5), Rational(3, 2)) == (2, Float(0.5))
    assert divmod(Rational(3, 2), Rational(1, 3)) == (4, Float(0.16666666666666666))
    assert divmod(Rational(1, 3), Rational(3, 2)) == (0, Rational(1, 3))
    assert divmod(Rational(3, 2), Float(0.1)) == (15, 0)
    assert divmod(Float(0.1), Rational(3, 2)) == (0, Float(0.1))
    assert divmod(Rational(3, 2), 2) == (0, Rational(3, 2))
    assert divmod(2, Rational(3, 2)) == (1, Float(0.5))
    assert divmod(Rational(3, 2), 1.5) == (1, 0)
    assert divmod(1.5, Rational(3, 2)) == (1, 0)
    assert divmod(Rational(3, 2), 0.3) == (5, 0)
    assert divmod(0.3, Rational(3, 2)) == (0, Float(0.3))
    assert divmod(Rational(1, 3), Float(3.5)) == (0, Rational(1, 3))
    assert divmod(Float(3.5), Float(0.1)) == (35, 0)
    assert divmod(Float(0.1), Float(3.5)) == (0, Float(0.1))
    assert divmod(Float(3.5), 2) == (1, Float(1.5))
    assert divmod(2, Float(3.5)) == (0, 2)
    assert divmod(Float(3.5), 1.5) == (2, Float(0.5))
    assert divmod(1.5, Float(3.5)) == (0, Float(1.5))
    assert divmod(0.3, Float(3.5)) == (0, Float(0.3))
    assert divmod(Float(0.1), Rational(1, 3)) == (0, Float(0.1))
    assert divmod(Rational(1, 3), 2) == (0, Rational(1, 3))
    assert divmod(2, Rational(1, 3)) == (6, 0)
    assert divmod(Rational(1, 3), 1.5) == (0, Rational(1, 3))
    assert divmod(0.3, Rational(1, 3)) == (0, Float(0.3))
    assert divmod(Float(0.1), 2) == (0, Float(0.1))
    assert divmod(2, Float(0.1)) == (20, 0)
    assert divmod(Float(0.1), 1.5) == (0, Float(0.1))
    assert divmod(1.5, Float(0.1)) == (15, 0)
    assert divmod(Float(0.1), 0.3) == (0, Float(0.1))

    assert str(divmod(Integer(2), 0.3)) == '(6, 0.2)'
    assert str(divmod(Float(3.5), Rational(1, 3))) == '(10, 0.166666666666667)'
    assert str(divmod(Float(3.5), 0.3)) == '(11, 0.2)'
    assert str(divmod(Rational(1, 3), Float(0.1))) == '(3, 0.0333333333333333)'
    assert str(divmod(1.5, Rational(1, 3))) == '(4, 0.166666666666667)'
    assert str(divmod(Rational(1, 3), 0.3)) == '(1, 0.0333333333333333)'
    assert str(divmod(0.3, Float(0.1))) == '(2, 0.1)'

    assert divmod(-3, Integer(2)) == (-2, 1)
    assert divmod(Integer(-3), Integer(2)) == (-2, 1)
    assert divmod(Integer(-3), 2) == (-2, 1)

    pytest.raises(ZeroDivisionError, lambda: divmod(oo, 0))

    # issue sympy/sympy#15561
    assert divmod(Integer(4), Float(-2.1)) == divmod(4, -2.1)
    f = fractions.Fraction(-31, 10)
    assert divmod(Integer(4), Rational(f)) == divmod(4, f)


def test_igcdex():
    assert igcdex(0, 0) == (0, 1, 0)
    assert igcdex(-2, 0) == (-1, 0, 2)
    assert igcdex(0, -2) == (0, -1, 2)
    assert igcdex(2, 3) == (-1, 1, 1)
    assert igcdex(10, 12) == (-1, 1, 2)
    assert igcdex(100, 2004) == (-20, 1, 4)
    assert igcdex(100, -2004) == (-20, -1, 4)


def _strictly_equal(a, b):
    return (a.numerator, a.denominator, type(a.numerator), type(a.denominator)) == \
           (b.numerator, b.denominator, type(b.numerator), type(b.denominator))


def _test_rational_new(cls):
    """Tests that are common between Integer and Rational."""
    for a, b in ((0, Integer(0)), (1, Integer(1)), (-1, Integer(-1)),
                 # These look odd, but are similar to int():
                 ('1', Integer(1)), ('-1', Integer(-1))):
        clear_cache()
        assert cls(a) is b

    i = Integer(10)
    assert _strictly_equal(i, cls('10'))
    assert _strictly_equal(i, cls('10'))
    assert _strictly_equal(i, cls(int(10)))
    assert _strictly_equal(i, cls(i))

    pytest.raises(TypeError, lambda: cls(Symbol('x')))


def test_Integer_new():
    _test_rational_new(Integer)

    assert _strictly_equal(Integer(0.9), Integer(0))
    assert _strictly_equal(Integer(10.5), Integer(10))
    pytest.raises(ValueError, lambda: Integer('10.5'))
    assert Integer(Rational('1.' + '9'*20)) == 1


def test_Rational_new():
    _test_rational_new(Rational)

    n1 = Rational(1, 2)
    assert n1 == Rational(Integer(1), 2)
    assert n1 == Rational(Integer(1), Integer(2))
    assert n1 == Rational(1, Integer(2))
    assert n1 == Rational(Rational(1, 2))
    assert 1 == Rational(n1, n1)
    assert Rational(3, 2) == Rational(Rational(1, 2), Rational(1, 3))
    assert Rational(3, 1) == Rational(1, Rational(1, 3))
    n3_4 = Rational(3, 4)
    assert Rational('3/4') == n3_4
    assert -Rational('-3/4') == n3_4
    assert Rational('.76').limit_denominator(4) == n3_4
    assert Rational(19, 25).limit_denominator(4) == n3_4
    assert Rational('19/25').limit_denominator(4) == n3_4
    assert Rational(1.0, 3) == Rational(1, 3)
    assert Rational(1, 3.0) == Rational(1, 3)
    assert Rational(Float(0.5)) == Rational(1, 2)
    assert Rational(-1, 0) == zoo
    assert Rational(+1, 0) == zoo
    pytest.raises(TypeError, lambda: Rational('3**3'))
    pytest.raises(TypeError, lambda: Rational('1/2 + 2/3'))

    assert Rational(2, 4).numerator == 1
    assert Rational(2, 4).denominator == 2

    # handle fractions.Fraction instances
    assert Rational(fractions.Fraction(1, 2)) == Rational(1, 2)

    assert Rational(pi.evalf(100)).evalf(100) == pi.evalf(100)


def test_Number_new():
    # Expected behavior on numbers and strings
    assert Number(1) is Integer(1)
    assert Number(2).__class__ is Integer
    assert Number(-622).__class__ is Integer
    assert Number(5, 3).__class__ is Rational
    assert Number(5.3).__class__ is Float
    assert Number('1') is Integer(1)
    assert Number('2').__class__ is Integer
    assert Number('-622').__class__ is Integer
    assert Number('5/3').__class__ is Rational
    assert Number('5.3').__class__ is Float
    pytest.raises(ValueError, lambda: Number('cos'))
    pytest.raises(TypeError, lambda: Number(cos))
    a = Rational(3, 5)
    assert Number(a) is a  # Check idempotence on Numbers


def test_Rational_cmp():
    n1 = Rational(1, 4)
    n2 = Rational(1, 3)
    n3 = Rational(2, 4)
    n4 = Rational(2, -4)
    n5 = Integer(0)
    n6 = Integer(1)
    n7 = Integer(3)
    n8 = Integer(-3)

    assert n8 < n5
    assert n5 < n6
    assert n6 < n7
    assert n8 < n7
    assert n7 > n8
    assert (n1 + 1)**n2 < 2
    assert ((n1 + n6)/n7) < 1

    assert n4 < n3
    assert n2 < n3
    assert n1 < n2
    assert n3 > n1
    assert n3 >= n1
    assert Integer(-1) <= 0
    assert Integer(-1) < 0

    pytest.raises(TypeError, lambda: n1 < nan)
    pytest.raises(TypeError, lambda: n1 <= nan)
    pytest.raises(TypeError, lambda: n1 > nan)
    pytest.raises(TypeError, lambda: n1 >= nan)


def test_Float():
    def eq(a, b):
        t = Float('1.0E-15')
        return -t < a - b < t

    a = Float(2) ** Float(3)
    assert eq(a, Float(8))
    assert eq((pi ** -1).evalf(), Float('0.31830988618379067'))
    a = Float(2) ** Float(4)
    assert eq(a, Float(16))
    assert (Float(.3) == Float(.5)) is False
    x_dec = Float((0, 5404319552844595, -52, 53))
    x2_dec = Float(x_dec)
    assert x_dec == x2_dec == Float(1.2)

    assert Float((0, int(0), -123, -1)) == Float('nan')
    assert Float((0, int(0), -456, -2)) == Float('inf') == Float('+inf')
    assert Float((1, int(0), -789, -3)) == Float('-inf')

    assert Float((0, int(7), 1, 3)) == Float('14.0', dps=15)

    assert Float('+inf').is_finite is False
    assert Float('+inf').is_negative is False
    assert Float('+inf').is_positive is True
    assert Float('+inf').is_infinite is True
    assert Float('+inf').is_nonzero is True

    assert Float('-inf').is_finite is False
    assert Float('-inf').is_negative is True
    assert Float('-inf').is_positive is False
    assert Float('-inf').is_infinite is True
    assert Float('-inf').is_nonzero is True

    assert Float('0.0').is_finite is True
    assert Float('0.0').is_negative is False
    assert Float('0.0').is_positive is False
    assert Float('0.0').is_infinite is False
    assert Float('0.0').is_zero is True

    # rationality properties
    assert Float(1).is_rational is None
    assert Float(1).is_irrational is None
    assert sqrt(2).evalf(15).is_rational is None
    assert sqrt(2).evalf(15).is_irrational is None

    # do not automatically evalf
    def teq(a):
        assert (a.evalf(strict=False) == a) is False
        assert (a.evalf(strict=False) != a) is True
        assert (a == a.evalf(strict=False)) is False
        assert (a != a.evalf(strict=False)) is True

    teq(pi)
    teq(2*pi)
    teq(cos(0.1, evaluate=False))

    # integer
    i = 12345678901234567890
    assert same_and_same_prec(Float(12), Float('12'))
    assert same_and_same_prec(Float(Integer(i)), Float(i))
    assert same_and_same_prec(Float(i), Float(str(i), 20))
    assert same_and_same_prec(Float(str(i)), Float(i))
    assert same_and_same_prec(Float(i), Float(i))

    # inexact floats (repeating binary = denom not multiple of 2)
    # cannot have precision greater than 15
    assert Float(.125, 22) == .125
    assert Float(2.0, 22) == 2
    assert float(Float('.12500000000000001')) == .125
    assert Float(.12500000000000001) == Float('0.125', dps=15)

    # allow auto precision detection
    assert Float('.1') == Float(.1, 1)
    assert Float('.125') == Float(.125, 3)
    assert Float('.100') == Float(.1, 3)
    assert Float('2.0') == Float('2', 2)

    pytest.raises(ValueError, lambda: Float('12.3d-4'))
    pytest.raises(ValueError, lambda: Float('.'))
    pytest.raises(ValueError, lambda: Float('-.'))

    assert Float(12.3) == Float('12.300000000000001', dps=15)

    zero = Float('0.0')
    assert Float('-0') == zero
    assert Float('.0') == zero
    assert Float('-.0') == zero
    assert Float('-0.0') == zero
    assert Float(0.0) == zero
    assert Float(0) == zero
    assert Float(0) == Float('0')
    assert Float(1) == Float(1.0)
    assert Float(Integer(0)) == zero
    assert Float(Integer(1)) == Float(1.0)

    assert Float(decimal.Decimal('0.1'), 3) == Float('.1', 3)
    assert Float(decimal.Decimal('nan')) == nan
    assert Float(decimal.Decimal('+Infinity')) == +oo
    assert Float(decimal.Decimal('-Infinity')) == -oo

    assert f'{Float(4.236622):.3f}' == '4.237'
    assert f'{Float(pi.evalf(40), 40):.35f}' == '3.14159265358979323846264338327950288'

    assert Float(+oo) == Float('+inf')
    assert Float(-oo) == Float('-inf')

    assert Float(0)**2 is Integer(0)

    t = Symbol('t', real=False)
    assert Float(0)**t == Pow(Float(0), t, evaluate=False)


def test_Float_eval():
    a = Float(3.2)
    assert (a**2).is_Float


def test_sympyissue_5206():
    a = Float(0.1, 10)
    b = Float('0.1', 10)

    assert a - a == 0
    assert a + (-a) == 0
    assert 0 + a - a == 0
    assert 0 + a + (-a) == 0

    assert b - b == 0
    assert b + (-b) == 0
    assert 0 + b - b == 0
    assert 0 + b + (-b) == 0


def test_Infinity():
    assert oo != 1
    assert 1*oo == oo
    assert 1 != oo
    assert oo != -oo
    assert oo != Symbol('x')**3
    assert oo + 1 == oo
    assert 2 + oo == oo
    assert 3*oo + 2 == oo
    assert Rational(1, 2)**oo == 0
    assert Rational(1, 2)**(-oo) == oo
    assert oo**zoo == nan
    assert -oo*3 == -oo
    assert oo + oo == oo
    assert -oo + oo*(-5) == -oo
    assert 1/oo == 0
    assert 1/(-oo) == 0
    assert 8/oo == 0
    assert Integer(8)/oo == 0
    assert Integer(8)/(-oo) == 0
    assert oo % 2 == nan
    assert 2 % oo == nan
    assert oo/oo == nan
    assert oo/-oo == nan
    assert -oo/oo == nan
    assert -oo/-oo == nan
    assert oo - oo == nan
    assert oo - -oo == oo
    assert -oo - oo == -oo
    assert -oo - -oo == nan
    assert oo + -oo == nan
    assert -oo + oo == nan
    assert oo + oo == oo
    assert -oo + oo == nan
    assert oo + -oo == nan
    assert -oo + -oo == -oo
    assert oo*oo == oo
    assert -oo*oo == -oo
    assert oo*-oo == -oo
    assert -oo*-oo == oo
    assert oo/0 == zoo
    assert -oo/0 == zoo
    assert 0/oo == 0
    assert 0/-oo == 0
    assert oo*0 == nan
    assert -oo*0 == nan
    assert 0*oo == nan
    assert 0*-oo == nan
    assert oo + 0 == oo
    assert -oo + 0 == -oo
    assert 0 + oo == oo
    assert 0 + -oo == -oo
    assert oo - 0 == oo
    assert -oo - 0 == -oo
    assert 0 - oo == -oo
    assert 0 - -oo == oo
    assert oo/2 == oo
    assert -oo/2 == -oo
    assert oo/-2 == -oo
    assert -oo/-2 == oo
    assert oo*2 == oo
    assert -oo*2 == -oo
    assert -oo*Float(0.0) == nan
    assert oo*-2 == -oo
    assert 2/oo == 0
    assert 2/-oo == 0
    assert -2/oo == 0
    assert -2/-oo == 0
    assert 2*oo == oo
    assert 2*-oo == -oo
    assert -2*oo == -oo
    assert -2*-oo == oo
    assert 2 + oo == oo
    assert 2 - oo == -oo
    assert -2 + oo == oo
    assert -2 - oo == -oo
    assert 2 + -oo == -oo
    assert 2 - -oo == oo
    assert -2 + -oo == -oo
    assert -2 - -oo == oo
    assert 2 + oo == oo
    assert 2 - oo == -oo
    assert oo/I == -oo*I
    assert -oo/I == oo*I
    assert oo*float(1) == Float('inf')
    assert (oo*float(1)).is_Float
    assert -oo*float(1) == Float('-inf')
    assert (-oo*float(1)).is_Float
    assert oo/float(1) == Float('inf')
    assert (oo/float(1)).is_Float
    assert -oo/float(1) == Float('-inf')
    assert (-oo/float(1)).is_Float
    assert oo*float(-1) == Float('-inf')
    assert (oo*float(-1)).is_Float
    assert -oo*float(-1) == Float('inf')
    assert (-oo*float(-1)).is_Float
    assert oo/float(-1) == Float('-inf')
    assert (oo/float(-1)).is_Float
    assert -oo/float(-1) == Float('inf')
    assert (-oo/float(-1)).is_Float
    assert oo + float(1) == Float('inf')
    assert (oo + float(1)).is_Float
    assert -oo + float(1) == Float('-inf')
    assert (-oo + float(1)).is_Float
    assert oo - float(1) == Float('inf')
    assert (oo - float(1)).is_Float
    assert -oo - float(1) == Float('-inf')
    assert (-oo - float(1)).is_Float
    assert float(1)*oo == Float('inf')
    assert (float(1)*oo).is_Float
    assert float(1)*-oo == Float('-inf')
    assert (float(1)*-oo).is_Float
    assert float(1)/oo == 0
    assert float(1)/-oo == 0
    assert float(-1)*oo == Float('-inf')
    assert (float(-1)*oo).is_Float
    assert float(-1)*-oo == Float('inf')
    assert (float(-1)*-oo).is_Float
    assert float(-1)/oo == 0
    assert float(-1)/-oo == 0
    assert float(1) + oo == Float('inf')
    assert float(1) + -oo == Float('-inf')
    assert float(1) - oo == Float('-inf')
    assert float(1) - -oo == Float('inf')
    assert oo*float(0) == nan

    assert Float('nan') == nan
    assert nan*1.0 == nan
    assert -1.0*nan == nan
    assert nan*oo == nan
    assert nan*-oo == nan
    assert nan/oo == nan
    assert nan/-oo == nan
    assert nan + oo == nan
    assert nan + -oo == nan
    assert nan - oo == nan
    assert nan - -oo == nan
    assert -oo * 0 == nan

    assert oo*nan == nan
    assert -oo*nan == nan
    assert oo/nan == nan
    assert -oo/nan == nan
    assert oo + nan == nan
    assert -oo + nan == nan
    assert oo - nan == nan
    assert -oo - nan == nan
    assert 0 * oo == nan
    assert oo.is_Rational is False
    assert isinstance(oo, Rational) is False

    assert 1/oo == 0
    assert -1/oo == 0
    assert 1/-oo == 0
    assert -1/-oo == 0
    assert 1*oo == oo
    assert -1*oo == -oo
    assert 1*-oo == -oo
    assert -1*-oo == oo
    assert 1/nan == nan
    assert 1 - -oo == oo
    assert 1 + nan == nan
    assert 1 - nan == nan
    assert nan - 1 == nan
    assert nan/1 == nan
    assert -oo - 1 == -oo

    e = I + cos(1)**2 + sin(1)**2 - 1
    assert oo**e == Pow(oo, e, evaluate=False)

    # issue sympy/sympy#10020
    assert oo**I is nan
    assert oo**(1 + I) is zoo
    assert oo**(-1 + I) is Integer(0)
    assert (-oo)**I is nan
    assert (-oo)**(-1 + I) is Integer(0)
    t = Symbol('t', real=False)
    assert oo**t == Pow(oo, t, evaluate=False)
    assert (-oo)**t == Pow(-oo, t, evaluate=False)

    # issue sympy/sympy#7742
    assert -oo % 1 == nan


def test_Infinity_2():
    assert oo*x != oo
    assert oo*(pi - 1) == oo
    assert oo*(1 - pi) == -oo

    assert (-oo)*x != -oo
    assert (-oo)*(pi - 1) == -oo
    assert (-oo)*(1 - pi) == oo

    assert (-1)**nan is nan
    assert oo - Float('inf') is nan
    assert oo + Float('-inf') is nan
    assert oo*0 is nan
    assert oo/Float('inf') is nan
    assert oo/Float('-inf') is nan
    assert oo**nan is nan
    assert -oo + Float('inf') is nan
    assert -oo - Float('-inf') is nan
    assert -oo*nan is nan
    assert -oo*0 is nan
    assert -oo/Float('inf') is nan
    assert -oo/Float('-inf') is nan
    assert -oo/nan is nan
    assert abs(-oo) == oo
    assert all((-oo)**i is nan for i in (oo, -oo, nan))
    assert (-oo)**3 == -oo
    assert (-oo)**2 == oo
    assert abs(zoo) == oo


def test_Mul_Infinity_Zero():
    assert 0*Float('inf') == nan
    assert 0*Float('-inf') == nan
    assert 0*Float('inf') == nan
    assert 0*Float('-inf') == nan
    assert Float('inf')*0 == nan
    assert Float('-inf')*0 == nan
    assert Float('inf')*0 == nan
    assert Float('-inf')*0 == nan
    assert Float(0)*Float('inf') == nan
    assert Float(0)*Float('-inf') == nan
    assert Float(0)*Float('inf') == nan
    assert Float(0)*Float('-inf') == nan
    assert Float('inf')*Float(0) == nan
    assert Float('-inf')*Float(0) == nan
    assert Float('inf')*Float(0) == nan
    assert Float('-inf')*Float(0) == nan


def test_Div_By_Zero():
    assert 1/Integer(0) == zoo
    assert 1/Float(0) == Float('inf')
    assert 0/Integer(0) == nan
    assert 0/Float(0) == nan
    assert Rational(0, 0) == nan
    assert Float(0)/0 == nan
    assert Float(0)/Float(0) == nan
    assert -1/Integer(0) == zoo
    assert -1/Float(0) == Float('-inf')


def test_Infinity_inequations():
    assert oo > pi
    assert not oo < pi
    assert exp(-3) < oo

    assert Float('+inf') > pi
    assert not Float('+inf') < pi
    assert exp(-3) < Float('+inf')

    pytest.raises(TypeError, lambda: oo < I)
    pytest.raises(TypeError, lambda: oo <= I)
    pytest.raises(TypeError, lambda: oo > I)
    pytest.raises(TypeError, lambda: oo >= I)
    pytest.raises(TypeError, lambda: -oo < I)
    pytest.raises(TypeError, lambda: -oo <= I)
    pytest.raises(TypeError, lambda: -oo > I)
    pytest.raises(TypeError, lambda: -oo >= I)

    pytest.raises(TypeError, lambda: I < oo)
    pytest.raises(TypeError, lambda: I <= oo)
    pytest.raises(TypeError, lambda: I > oo)
    pytest.raises(TypeError, lambda: I >= oo)
    pytest.raises(TypeError, lambda: I < -oo)
    pytest.raises(TypeError, lambda: I <= -oo)
    pytest.raises(TypeError, lambda: I > -oo)
    pytest.raises(TypeError, lambda: I >= -oo)

    assert oo > -oo
    assert oo >= -oo
    assert (oo < -oo) is false
    assert (oo <= -oo) is false
    assert -oo < oo
    assert -oo <= oo
    assert (-oo > oo) is false
    assert (-oo >= oo) is false

    # pylint: disable=comparison-with-itself

    assert (oo < oo) is false  # issue sympy/sympy#7775
    assert (oo > oo) is false
    assert (-oo > -oo) is false
    assert (-oo < -oo) is false
    assert oo >= oo
    assert oo <= oo
    assert -oo >= -oo
    assert -oo <= -oo

    b = Symbol('b', real=True)
    assert (x < oo) == Lt(x, oo)  # issue sympy/sympy#7775
    assert -oo < b < oo
    assert -oo <= b <= oo
    assert oo > b
    assert oo >= b
    assert (oo < b) is false
    assert (oo <= b) is false
    assert (-oo > b) is false
    assert (-oo >= b) is false
    assert -oo < b
    assert -oo <= b
    assert (oo < x) == Lt(oo, x)
    assert (oo > x) == Gt(oo, x)
    assert (oo <= x) == Le(oo, x)
    assert (oo >= x) == Ge(oo, x)
    assert (-oo < x) == Lt(-oo, x)
    assert (-oo > x) == Gt(-oo, x)
    assert (-oo <= x) == Le(-oo, x)
    assert (-oo >= x) == Ge(-oo, x)


def test_NaN():
    assert nan == nan  # pylint: disable=comparison-with-itself
    assert nan != 1
    assert 1*nan == nan
    assert 1 != nan
    assert nan == -nan
    assert oo != Symbol('x')**3
    assert nan + 1 == nan
    assert 2 + nan == nan
    assert Integer(2) + nan == nan
    assert 3*nan + 2 == nan
    assert -nan*3 == nan
    assert nan + nan == nan
    assert -nan + nan*(-5) == nan
    assert 1/nan == nan
    assert 1/(-nan) == nan
    assert 8/nan == nan
    assert Integer(8)/nan == nan
    pytest.raises(TypeError, lambda: nan > 0)
    pytest.raises(TypeError, lambda: nan < 0)
    pytest.raises(TypeError, lambda: nan >= 0)
    pytest.raises(TypeError, lambda: nan <= 0)
    pytest.raises(TypeError, lambda: 0 < nan)
    pytest.raises(TypeError, lambda: 0 > nan)
    pytest.raises(TypeError, lambda: 0 <= nan)
    pytest.raises(TypeError, lambda: 0 >= nan)
    assert 1 + nan == nan
    assert 1 - nan == nan
    assert 1*nan == nan
    assert 1/nan == nan
    assert nan - 1 == nan
    assert nan*1 == nan
    assert nan + 1 == nan
    assert nan/1 == nan
    assert nan**0 == nan
    assert 1**nan == nan
    # test Pow._eval_power's handling of NaN
    assert Pow(nan, 0, evaluate=False)**2 == nan
    assert Integer(0)/Integer(0) == nan


def test_special_numbers():
    assert isinstance(nan, Number) is True
    assert isinstance(+oo, Number) is True
    assert isinstance(-oo, Number) is True

    assert nan.is_number is True
    assert (+oo).is_number is True
    assert (-oo).is_number is True
    assert zoo.is_number is True

    assert isinstance(nan, Rational) is False
    assert isinstance(+oo, Rational) is False
    assert isinstance(-oo, Rational) is False

    assert nan.is_rational is not True
    assert (+oo).is_rational is not True
    assert (-oo).is_rational is not True


def test_powers():
    assert integer_nthroot(1, 2) == (1, True)
    assert integer_nthroot(1, 5) == (1, True)
    assert integer_nthroot(2, 1) == (2, True)
    assert integer_nthroot(2, 2) == (1, False)
    assert integer_nthroot(2, 5) == (1, False)
    assert integer_nthroot(4, 2) == (2, True)
    assert integer_nthroot(123**25, 25) == (123, True)
    assert integer_nthroot(123**25 + 1, 25) == (123, False)
    assert integer_nthroot(123**25 - 1, 25) == (122, False)
    assert integer_nthroot(1, 1) == (1, True)
    assert integer_nthroot(0, 1) == (0, True)
    assert integer_nthroot(0, 3) == (0, True)
    assert integer_nthroot(10000, 1) == (10000, True)
    assert integer_nthroot(4, 2) == (2, True)
    assert integer_nthroot(16, 2) == (4, True)
    assert integer_nthroot(26, 2) == (5, False)
    assert integer_nthroot(1234567**7, 7) == (1234567, True)
    assert integer_nthroot(1234567**7 + 1, 7) == (1234567, False)
    assert integer_nthroot(1234567**7 - 1, 7) == (1234566, False)
    b = 25**1000
    assert integer_nthroot(b, 1000) == (25, True)
    assert integer_nthroot(b + 1, 1000) == (25, False)
    assert integer_nthroot(b - 1, 1000) == (24, False)
    c = 10**400
    c2 = c**2
    assert integer_nthroot(c2, 2) == (c, True)
    assert integer_nthroot(c2 + 1, 2) == (c, False)
    assert integer_nthroot(c2 - 1, 2) == (c - 1, False)
    assert integer_nthroot(2, 10**10) == (1, False)

    p, r = integer_nthroot(int(factorial(10000)), 100)
    assert p % (10**10) == 5322420655
    assert not r

    # Test that this is fast
    assert integer_nthroot(2, 10**10) == (1, False)

    pytest.raises(ValueError, lambda: integer_nthroot(-1, 2))
    pytest.raises(ValueError, lambda: integer_nthroot(2, 0))

    # output should be int if possible
    assert type(integer_nthroot(2**61, 2)[0]) is int


def test_integer_nthroot_overflow():
    assert integer_nthroot(10**(50*50), 50) == (10**50, True)
    assert integer_nthroot(10**100000, 10000) == (10**10, True)


def test_powers_Integer():
    # check infinity
    assert (+1) ** oo == nan
    assert (-1) ** oo == nan
    assert (+2) ** oo == oo
    assert (-2) ** oo == oo + oo * I
    assert 0 ** oo == 0

    # check Nan
    assert (+1) ** nan == nan
    assert (-1) ** nan == nan

    # check for exact roots
    assert (-1) ** Rational(6, 5) == - root(-1, 5)
    assert sqrt(4) == 2
    assert sqrt(-4) == I * 2
    assert root(16, 4) == 2
    assert root(-16, 4) == 2 * root(-1, 4)
    assert 9 ** Rational(3, 2) == 27
    assert (-9) ** Rational(3, 2) == -27*I
    assert 27 ** Rational(2, 3) == 9
    assert (-27) ** Rational(2, 3) == 9 * ((-1) ** Rational(2, 3))
    assert (-2) ** Rational(-2, 1) == Rational(1, 4)

    # not exact roots
    assert sqrt(-3) == I*sqrt(3)
    assert 3 ** Rational(3, 2) == 3 * sqrt(3)
    assert (-3) ** Rational(3, 2) == - 3 * sqrt(-3)
    assert (-3) ** Rational(5, 2) == 9 * I * sqrt(3)
    assert (-3) ** Rational(7, 2) == - I * 27 * sqrt(3)
    assert 2 ** Rational(3, 2) == 2 * sqrt(2)
    assert 2 ** Rational(-3, 2) == sqrt(2) / 4
    assert 81 ** Rational(2, 3) == 9 * 3 ** Rational(2, 3)
    assert (-81) ** Rational(2, 3) == 9 * (-3) ** Rational(2, 3)
    assert (-3) ** Rational(-7, 3) == -(-1)**Rational(2, 3)*3**Rational(2, 3)/27
    assert (-3) ** Rational(-2, 3) == -cbrt(-1)*cbrt(3)/3

    # join roots
    assert sqrt(6) + sqrt(24) == 3*sqrt(6)
    assert sqrt(2) * sqrt(3) == sqrt(6)

    # separate symbols & constansts
    assert sqrt(49 * x) == 7 * sqrt(x)
    assert sqrt((3 - sqrt(pi)) ** 2) == 3 - sqrt(pi)

    # check that it is fast for big numbers
    assert (2**64 + 1) ** Rational(4, 3)
    assert (2**64 + 1) ** Rational(17, 25)

    # negative rational power and negative base
    assert (-3) ** Rational(-7, 3) == -(-1)**Rational(2, 3)*3**Rational(2, 3)/27
    assert (-3) ** Rational(-2, 3) == -cbrt(-1)*cbrt(3)/3

    assert Integer(1).factors(visual=True) == 1
    assert Integer(1234).factors() == {617: 1, 2: 1}
    assert Rational(2*3, 3*5*7).factors() == {2: 1, 5: -1, 7: -1}

    # test that eval_power factors numbers bigger than
    # the current limit in factor_trial_division (2**15)
    n = nextprime(2**15)
    assert sqrt(n**2) == n
    assert sqrt(n**3) == n*sqrt(n)
    assert sqrt(4*n) == 2*sqrt(n)

    # check that factors of base with powers sharing gcd with power are removed
    assert root(2**4*3, 6) == 2**Rational(2, 3)*root(3, 6)
    assert (2**4*3)**Rational(5, 6) == 8*cbrt(2)*3**Rational(5, 6)

    # check that bases sharing a gcd are exptracted
    assert cbrt(2)*root(3, 4)*root(6, 5) == \
        2**Rational(8, 15)*3**Rational(9, 20)
    assert sqrt(8)*cbrt(24)*root(6, 5) == \
        4*2**Rational(7, 10)*3**Rational(8, 15)
    assert sqrt(8)*cbrt(-24)*root(-6, 5) == \
        4*(-3)**Rational(8, 15)*2**Rational(7, 10)
    assert cbrt(2)*2**Rational(8, 9) == 2*2**Rational(2, 9)
    assert 2**Rational(2, 3)*cbrt(6) == 2*cbrt(3)
    assert 2**Rational(2, 3)*6**Rational(8, 9) == \
        2*2**Rational(5, 9)*3**Rational(8, 9)
    assert (-2)**Rational(2, 3)*cbrt(-4) == -2*cbrt(2)
    assert 3*Pow(3, 2, evaluate=False) == 3**3
    assert 3*Pow(3, Rational(-1, 3), evaluate=False) == 3**Rational(2, 3)
    assert (-2)**Rational(1, 3)*(-3)**Rational(1, 4)*(-5)**Rational(5, 6) == \
        -(-1)**Rational(5, 12)*cbrt(2)*root(3, 4) * \
        5**Rational(5, 6)

    assert (-2)**Symbol('', even=True) == 2**Symbol('', even=True)
    assert (-1)**Float(.5) == 1.0*I

    n = Symbol('n', integer=True)
    e = (-1)**n/2 + Rational(5, 2)
    assert (-1)**e == Pow(-1, e, evaluate=False)


def test_powers_Rational():
    # check infinity
    assert Rational(1, 2) ** oo == 0
    assert Rational(3, 2) ** oo == oo
    assert Rational(-1, 2) ** oo == 0
    assert Rational(-3, 2) ** oo == oo + oo*I

    # check Nan
    assert Rational(3, 4) ** nan == nan
    assert Rational(-2, 3) ** nan == nan

    # exact roots on numerator
    assert sqrt(Rational(4, 3)) == 2 * sqrt(3) / 3
    assert Rational(4, 3) ** Rational(3, 2) == 8 * sqrt(3) / 9
    assert sqrt(Rational(-4, 3)) == I * 2 * sqrt(3) / 3
    assert Rational(-4, 3) ** Rational(3, 2) == - I * 8 * sqrt(3) / 9
    assert cbrt(Rational(27, 2)) == 3 * (2 ** Rational(2, 3)) / 2
    assert Rational(5**3, 8**3) ** Rational(4, 3) == Rational(5**4, 8**4)

    # exact root on denominator
    assert sqrt(Rational(1, 4)) == Rational(1, 2)
    assert sqrt(Rational(1, -4)) == I/2
    assert sqrt(Rational(3, 4)) == sqrt(3)/2
    assert sqrt(Rational(3, -4)) == I*sqrt(3)/2
    assert cbrt(Rational(5, 27)) == cbrt(5)/3

    # not exact roots
    assert sqrt(Rational(1, 2)) == sqrt(2)/2
    assert sqrt(Rational(-4, 7)) == I * sqrt(Rational(4, 7))
    assert Rational(-3, 2)**Rational(-7, 3) == \
        -4*(-1)**Rational(2, 3)*cbrt(2)*3**Rational(2, 3)/27
    assert Rational(-3, 2)**Rational(-2, 3) == \
        -cbrt(-1)*2**Rational(2, 3)*cbrt(3)/3

    # negative integer power and negative rational base
    assert Rational(-2, 3) ** Rational(-2, 1) == Rational(9, 4)

    a = Rational(1, 10)
    assert a**Float(a, 2) == Float(a, 2)**Float(a, 2)
    assert Rational(-2, 3)**Symbol('', even=True) == \
        Rational(2, 3)**Symbol('', even=True)


def test_powers_Float():
    assert str((Rational(-1, 10)**Rational(3, 10)).evalf()) == str(Float(-.1)**(.3))


def test_abs1():
    assert Rational(1, 6) != Rational(-1, 6)
    assert abs(Rational(1, 6)) == abs(Rational(-1, 6))


def test_accept_int():
    assert Float(4) == 4


def test_dont_accept_str():
    assert Float('0.2') != '0.2'
    assert not Float('0.2') == '0.2'  # pylint: disable=unneeded-not


def test_int():
    a = Integer(5)
    assert int(a) == 5
    a = Rational(9, 10)
    assert int(a) == int(-a) == 0
    assert 1/(-1)**Rational(2, 3) == -cbrt(-1)
    assert int(pi) == 3
    assert int(E) == 2
    assert int(GoldenRatio) == 1

    a = Integer(2**100)
    assert int(a) == a


def test_real_bug():
    assert str(2.0*x*x) in ['(2.0*x)*x', '2.0*x**2', '2.00000000000000*x**2']
    assert str(2.1*x*x) != '(2.0*x)*x'


def test_bug_sqrt():
    assert ((sqrt(2) + 1)*(sqrt(2) - 1)).expand() == 1


def test_pi_Pi():
    # Test that pi (instance) is imported, but Pi (class) is not.
    with pytest.raises(ImportError):
        # pylint: disable=unused-import,no-name-in-module
        from diofant import Pi  # noqa: F401


def test_no_len():
    # there should be no len for numbers
    pytest.raises(TypeError, lambda: len(Rational(2, 3)))
    pytest.raises(TypeError, lambda: len(Integer(2)))


def test_sympyissue_3321():
    assert sqrt(Rational(1, 5)) == sqrt(Rational(1, 5))
    assert 5 * sqrt(Rational(1, 5)) == sqrt(5)


def test_sympyissue_3692():
    assert root(-1, 6).expand(complex=True) == I/2 + sqrt(3)/2
    assert root(-5, 6).expand(complex=True) == root(5, 6)*I/2 + root(5, 6)*sqrt(3)/2
    assert root(-64, 6).expand(complex=True) == I + sqrt(3)


def test_sympyissue_3423():
    assert sqrt(x - 1).as_base_exp() == (x - 1, Rational(1, 2))
    assert sqrt(x - 1) != I*sqrt(1 - x)


def test_Integer_factors():
    def F(i):
        return Integer(i).factors()

    assert F(1) == {}
    assert F(2) == {2: 1}
    assert F(3) == {3: 1}
    assert F(4) == {2: 2}
    assert F(5) == {5: 1}
    assert F(6) == {2: 1, 3: 1}
    assert F(7) == {7: 1}
    assert F(8) == {2: 3}
    assert F(9) == {3: 2}
    assert F(10) == {2: 1, 5: 1}
    assert F(11) == {11: 1}
    assert F(12) == {2: 2, 3: 1}
    assert F(13) == {13: 1}
    assert F(14) == {2: 1, 7: 1}
    assert F(15) == {3: 1, 5: 1}
    assert F(16) == {2: 4}
    assert F(17) == {17: 1}
    assert F(18) == {2: 1, 3: 2}
    assert F(19) == {19: 1}
    assert F(20) == {2: 2, 5: 1}
    assert F(21) == {3: 1, 7: 1}
    assert F(22) == {2: 1, 11: 1}
    assert F(23) == {23: 1}
    assert F(24) == {2: 3, 3: 1}
    assert F(25) == {5: 2}
    assert F(26) == {2: 1, 13: 1}
    assert F(27) == {3: 3}
    assert F(28) == {2: 2, 7: 1}
    assert F(29) == {29: 1}
    assert F(30) == {2: 1, 3: 1, 5: 1}
    assert F(31) == {31: 1}
    assert F(32) == {2: 5}
    assert F(33) == {3: 1, 11: 1}
    assert F(34) == {2: 1, 17: 1}
    assert F(35) == {5: 1, 7: 1}
    assert F(36) == {2: 2, 3: 2}
    assert F(37) == {37: 1}
    assert F(38) == {2: 1, 19: 1}
    assert F(39) == {3: 1, 13: 1}
    assert F(40) == {2: 3, 5: 1}
    assert F(41) == {41: 1}
    assert F(42) == {2: 1, 3: 1, 7: 1}
    assert F(43) == {43: 1}
    assert F(44) == {2: 2, 11: 1}
    assert F(45) == {3: 2, 5: 1}
    assert F(46) == {2: 1, 23: 1}
    assert F(47) == {47: 1}
    assert F(48) == {2: 4, 3: 1}
    assert F(49) == {7: 2}
    assert F(50) == {2: 1, 5: 2}
    assert F(51) == {3: 1, 17: 1}


def test_Rational_factors():
    def F(p, q, visual=None):
        return Rational(p, q).factors(visual=visual)

    assert F(2, 3) == {2: 1, 3: -1}
    assert F(2, 9) == {2: 1, 3: -2}
    assert F(2, 15) == {2: 1, 3: -1, 5: -1}
    assert F(6, 10) == {3: 1, 5: -1}


def test_sympyissue_4107():
    assert pi*(E + 10) + pi*(-E - 10) != 0
    assert pi*(E + 10**10) + pi*(-E - 10**10) != 0
    assert pi*(E + 10**20) + pi*(-E - 10**20) != 0
    assert pi*(E + 10**80) + pi*(-E - 10**80) != 0

    assert (pi*(E + 10) + pi*(-E - 10)).expand() == 0
    assert (pi*(E + 10**10) + pi*(-E - 10**10)).expand() == 0
    assert (pi*(E + 10**20) + pi*(-E - 10**20)).expand() == 0
    assert (pi*(E + 10**80) + pi*(-E - 10**80)).expand() == 0


def test_IntegerInteger():
    a = Integer(4)
    b = Integer(a)

    assert a == b


def test_Rational_gcd_lcm_cofactors():
    assert Integer(4).gcd(2) == 2
    assert Integer(4).lcm(2) == 4
    assert Integer(4).gcd(Integer(2)) == 2
    assert Integer(4).lcm(Integer(2)) == 4

    assert Integer(4).gcd(3) == 1
    assert Integer(4).lcm(3) == 12
    assert Integer(4).gcd(Integer(3)) == 1
    assert Integer(4).lcm(Integer(3)) == 12

    assert Rational(4, 3).gcd(2) == Rational(2, 3)
    assert Rational(4, 3).lcm(2) == 4
    assert Rational(4, 3).gcd(Integer(2)) == Rational(2, 3)
    assert Rational(4, 3).lcm(Integer(2)) == 4

    assert Integer(4).gcd(Rational(2, 9)) == Rational(2, 9)
    assert Integer(4).lcm(Rational(2, 9)) == 4

    assert Rational(4, 3).gcd(Rational(2, 9)) == Rational(2, 9)
    assert Rational(4, 3).lcm(Rational(2, 9)) == Rational(4, 3)
    assert Rational(4, 5).gcd(Rational(2, 9)) == Rational(2, 45)
    assert Rational(4, 5).lcm(Rational(2, 9)) == 4

    assert Integer(4).cofactors(2) == (2, 2, 1)
    assert Integer(4).cofactors(Integer(2)) == (2, 2, 1)

    assert Integer(4).gcd(Float(2.0)) == 1
    assert Integer(4).lcm(Float(2.0)) == Float(8.0)
    assert Integer(4).cofactors(Float(2.0)) == (1, 4, Float(2.0))

    assert Rational(1, 2).gcd(Float(2.0)) == 1
    assert Rational(1, 2).lcm(Float(2.0)) == Float(1.0)
    assert Rational(1, 2).cofactors(Float(2.0)) == (1, Rational(1, 2), Float(2.0))


def test_Float_gcd_lcm_cofactors():
    assert Float(2.0).gcd(Integer(4)) == 1
    assert Float(2.0).lcm(Integer(4)) == Float(8.0)
    assert Float(2.0).cofactors(Integer(4)) == (1, Float(2.0), 4)

    assert Float(2.0).gcd(Rational(1, 2)) == 1
    assert Float(2.0).lcm(Rational(1, 2)) == Float(1.0)
    assert Float(2.0).cofactors(Rational(1, 2)) == \
        (1, Float(2.0), Rational(1, 2))


def test_sympyissue_4611():
    assert abs(pi._evalf(50) - 3.14159265358979) < 1e-10
    assert abs(E._evalf(50) - 2.71828182845905) < 1e-10
    assert abs(Catalan._evalf(50) - 0.915965594177219) < 1e-10
    assert abs(EulerGamma._evalf(50) - 0.577215664901533) < 1e-10
    assert abs(GoldenRatio._evalf(50) - 1.61803398874989) < 1e-10
    assert (pi + x).evalf(strict=False) == pi.evalf() + x
    assert (E + x).evalf(strict=False) == E.evalf() + x
    assert (Catalan + x).evalf(strict=False) == Catalan.evalf() + x
    assert (EulerGamma + x).evalf(strict=False) == EulerGamma.evalf() + x
    assert (GoldenRatio + x).evalf(strict=False) == GoldenRatio.evalf() + x


def test_conversion_to_mpmath():
    assert mpmath.mpmathify(Integer(1)) == mpmath.mpf(1)
    assert mpmath.mpmathify(Rational(1, 2)) == mpmath.mpf(0.5)
    assert mpmath.mpmathify(Float('1.23', 15)) == mpmath.mpf('1.23')

    assert mpmath.mpf(Rational(1, 3)) == mpmath.mpf('0.33333333333333331')


def test_relational():
    # real
    x = Float(.1)
    assert (x != cos) is True
    assert (x == cos) is False

    # rational
    x = Rational(1, 3)
    assert (x != cos) is True
    assert (x == cos) is False

    # integer defers to rational so these tests are omitted

    # number symbol
    x = pi
    assert (x != cos) is True
    assert (x == cos) is False

    r = Symbol('r', extended_real=True)
    assert (oo > r) == Gt(oo, r)
    assert (oo <= r) == Le(oo, r)
    assert (oo >= r) is true

    assert (Float(3.0) >= pi) is false
    assert (Float(3.0) <= pi) is true


def test_Integer_as_index():
    assert 'hello'[Integer(2):] == 'llo'


def test_Rational_int():
    assert int(+Rational(7, 5)) == 1
    assert int(+Rational(1, 2)) == 0
    assert int(-Rational(1, 2)) == 0
    assert int(-Rational(7, 5)) == -1


def test_zoo():
    b = Symbol('b', finite=True)
    nz = Symbol('nz', nonzero=True)
    p = Symbol('p', positive=True)
    n = Symbol('n', negative=True)
    im = Symbol('i', imaginary=True)
    c = Symbol('c', complex=True)
    pb = Symbol('pb', positive=True, finite=True)
    nb = Symbol('nb', negative=True, finite=True)
    imb = Symbol('ib', imaginary=True, finite=True)
    for i in [I, oo, -oo, Integer(0), Integer(1), pi, Rational(1, 2), Integer(3), log(3),
              b, nz, p, n, im, pb, nb, imb, c]:
        if i.is_finite and (i.is_extended_real or i.is_imaginary):
            assert i + zoo is zoo
            assert i - zoo is zoo
            assert zoo + i is zoo
            assert zoo - i is zoo
        elif i.is_finite is not False:
            assert (i + zoo).is_Add
            assert (i - zoo).is_Add
            assert (zoo + i).is_Add
            assert (zoo - i).is_Add
        else:
            assert (i + zoo) is nan
            assert (i - zoo) is nan
            assert (zoo + i) is nan
            assert (zoo - i) is nan

        if i.is_nonzero and (i.is_extended_real or i.is_imaginary):
            assert i*zoo is zoo
            assert zoo*i is zoo
        elif i.is_zero:
            assert i*zoo is nan
            assert zoo*i is nan
        else:
            assert (i*zoo).is_Mul
            assert (zoo*i).is_Mul

        if (1/i).is_nonzero and (i.is_extended_real or i.is_imaginary):
            assert zoo/i is zoo
        elif (1/i).is_zero:
            assert zoo/i is nan
        elif i.is_zero:
            assert zoo/i is zoo
        else:
            assert (zoo/i).is_Mul

    assert (I*oo).is_Mul  # allow directed infinity
    assert zoo + zoo is nan
    assert zoo * zoo is zoo
    assert zoo - zoo is nan
    assert zoo/zoo is nan
    assert zoo**zoo is nan
    assert zoo**0 is Integer(1)
    assert zoo**2 is zoo
    assert 1/zoo is Integer(0)

    assert Mul.flatten([Integer(-1), oo, Integer(0)]) == ([nan], [], None)


def test_sympyissue_4122():
    x = Symbol('x', nonpositive=True)
    assert (oo + x).is_Add
    x = Symbol('x', finite=True)
    assert (oo + x).is_Add  # x could be imaginary
    x = Symbol('x', nonnegative=True)
    assert oo + x == oo
    x = Symbol('x', real=True)
    assert oo + x == oo

    # similarily for negative infinity
    x = Symbol('x', nonnegative=True)
    assert (-oo + x).is_Add
    x = Symbol('x', finite=True)
    assert (-oo + x).is_Add
    x = Symbol('x', nonpositive=True)
    assert -oo + x == -oo
    x = Symbol('x', real=True)
    assert -oo + x == -oo


def test_GoldenRatio_expand():
    assert GoldenRatio.expand(func=True) == Rational(1, 2) + sqrt(5)/2


def test_as_content_primitive():
    assert Integer(0).as_content_primitive() == (1, 0)
    assert Rational(1, 2).as_content_primitive() == (Rational(1, 2), 1)
    assert Rational(-1, 2).as_content_primitive() == (Rational(1, 2), -1)
    assert Integer(3).as_content_primitive() == (3, 1)
    assert Float(3.1).as_content_primitive() == (1, 3.1)


def test_hashing_diofant_integers():
    # Test for issue sympy/sympy#5072
    assert {Integer(3)} == {int(3)}
    assert hash(Integer(4)) == hash(int(4))


def test_sympyissue_4172():
    assert int((E**100).round()) == 26881171418161354484126255515800135873611119
    assert int((pi**100).round()) == 51878483143196131920862615246303013562686760680406
    assert int((1/EulerGamma**100).round()) == 734833795660954410469466


def test_Catalan_EulerGamma_prec():
    n = GoldenRatio
    f = Float(n.evalf(), 5)
    assert f._mpf_ == (0, int(212079), -17, 18)
    assert f._prec == 20
    assert n._as_mpf_val(20) == f._mpf_

    n = EulerGamma
    f = Float(n.evalf(), 5)
    assert f._mpf_ == (0, int(302627), -19, 19)
    assert f._prec == 20
    assert n._as_mpf_val(20) == f._mpf_


def test_Float_eq():
    assert Float(.12, 3) != Float(.12, 4)
    assert Float(.12, 3) == .12
    assert 0.12 == Float(.12, 3)
    assert Float('.12', 22) != .12


def test_int_NumberSymbols():
    assert [int(i) for i in [pi, EulerGamma, E,
                             GoldenRatio, Catalan]] == [3, 0, 2, 1, 0]


def test_approximation_interval():
    assert E.approximation_interval(Integer) == (2, 3)
    assert E.approximation_interval(Float) is None

    assert GoldenRatio.approximation_interval(Integer) == (1, 2)
    assert GoldenRatio.approximation_interval(Float) is None

    assert EulerGamma.approximation_interval(Integer) == (0, 1)
    assert EulerGamma.approximation_interval(Rational) == (Rational(1, 2),
                                                           Rational(3, 5))
    assert EulerGamma.approximation_interval(Float) is None

    assert Catalan.approximation_interval(Integer) == (0, 1)
    assert Catalan.approximation_interval(Rational) == (Rational(9, 10), 1)
    assert Catalan.approximation_interval(Float) is None

    assert pi.approximation_interval(Integer) == (3, 4)
    assert pi.approximation_interval(Rational) == (Rational(223, 71),
                                                   Rational(22, 7))
    assert pi.approximation_interval(Float) is None


def test_sympyissue_6640():
    # fnan is not included because Float no longer returns fnan,
    # but otherwise, the same sort of test could apply
    assert Float(finf).is_nonzero is True
    assert Float(fninf).is_nonzero is True
    assert bool(Float(0)) is False


def test_sympyissue_6349():
    assert Float('23.e3')._prec == 10
    assert Float('23e3')._prec == 10
    assert Float('23000')._prec == 20
    assert Float('-23000')._prec == 20


def test_mpf_norm():
    assert mpf_norm((1, 0, 1, 0), 10) == mpmath.mpf('0')._mpf_
    assert Float._new((1, 0, 1, 0), 10)._mpf_ == mpmath.mpf('0')._mpf_


def test_latex():
    assert latex(pi) == r'\pi'
    assert latex(E) == r'e'
    assert latex(GoldenRatio) == r'\phi'
    assert latex(EulerGamma) == r'\gamma'
    assert latex(oo) == r'\infty'
    assert latex(-oo) == r'-\infty'
    assert latex(zoo) == r'\tilde{\infty}'
    assert latex(nan) == r'\mathrm{NaN}'
    assert latex(I) == r'i'


def test_Float_idempotence():
    x = Float('1.23')
    y = Float(x)
    z = Float(x, 15)
    assert same_and_same_prec(y, x)
    assert not same_and_same_prec(z, x)
    x = Float(10**20)
    y = Float(x)
    z = Float(x, 15)
    assert same_and_same_prec(y, x)
    assert not same_and_same_prec(z, x)


def test_comp():
    # sqrt(2) = 1.414213 5623730950...
    a = sqrt(2).evalf(7)
    assert comp(a, 1.41421346) is False
    assert comp(a, 1.41421347)
    assert comp(a, 1.41421366)
    assert comp(a, 1.41421367) is False
    assert comp(sqrt(2).evalf(2), '1.4')
    assert comp(sqrt(2).evalf(2), Float(1.4, 2), '')
    pytest.raises(ValueError, lambda: comp(sqrt(2).evalf(2), 1.4, ''))
    assert comp(sqrt(2).evalf(2), Float(1.4, 3)) is False
    pytest.raises(ValueError, lambda: comp('123', '123'))


def test_sympyissue_10063():
    assert 2**Float(3) == Float(8)


def test_invert_numbers():
    assert Integer(2).invert(5) == 3
    assert Integer(2).invert(Rational(5, 2)) == Rational(1, 2)
    assert Integer(2).invert(5.) == 3
    assert Integer(2).invert(Integer(5)) == 3
    assert Integer(2.).invert(5) == 3
    assert sqrt(2).invert(5) == 1/sqrt(2)
    assert sqrt(2).invert(sqrt(3)) == 1/sqrt(2)


def test_mod_inverse():
    assert mod_inverse(3, 11) == 4
    assert mod_inverse(5, 11) == 9
    assert mod_inverse(21124921, 521512) == 7713
    assert mod_inverse(124215421, 5125) == 2981
    assert mod_inverse(214, 12515) == 1579
    assert mod_inverse(5823991, 3299) == 1442
    assert mod_inverse(123, 44) == 39
    assert mod_inverse(2, 5) == 3
    assert mod_inverse(-2, 5) == -3
    assert Integer(2).invert(x) == Rational(1, 2)
    pytest.raises(TypeError, lambda: mod_inverse(2, x))
    pytest.raises(ValueError, lambda: mod_inverse(2, Rational(1, 2)))
    pytest.raises(ValueError, lambda: mod_inverse(2, cos(1)**2 + sin(1)**2))
    pytest.raises(ValueError, lambda: mod_inverse(2, 1))
    pytest.raises(ValueError, lambda: mod_inverse(3, 6))


def test_integer_digits():
    assert integer_digits(0, 2) == integer_digits(0, 3) == [0]
    assert integer_digits(10, 3) == [1, 0, 1]
    assert integer_digits(874881, 712) == [1, 516, 545]


def test_sympyissue_13081():
    r = Rational(905502432259640373, 288230376151711744)
    assert (pi < r) is true
    assert (r > pi) is true
    r2 = Rational(472202503979844695356573871761845338575143343779448489867569357017941709222155070092152068445390137810467671349,
                  150306725297525326584926758194517569752043683130132471725266622178061377607334940381676735896625196994043838464)
    assert (r2 < pi) is true
    assert (r2 > pi) is false


def test_comparisons_with_unknown_type():
    class Foo:
        # Class that is unaware of Basic, and relies on both classes returning
        # the NotImplemented singleton for equivalence to evaluate to False, and
        # the other comparisons to raise a TypeError.
        pass

    ni, nf, nr = Integer(3), Float(1.0), Rational(1, 3)
    foo = Foo()

    for n in ni, nf, nr, oo:
        assert n != foo
        assert foo != n
        assert not n == foo
        assert not foo == n
        pytest.raises(TypeError, lambda: n < foo)
        pytest.raises(TypeError, lambda: foo > n)
        pytest.raises(TypeError, lambda: n > foo)
        pytest.raises(TypeError, lambda: foo < n)
        pytest.raises(TypeError, lambda: n <= foo)
        pytest.raises(TypeError, lambda: foo >= n)
        pytest.raises(TypeError, lambda: n >= foo)
        pytest.raises(TypeError, lambda: foo <= n)

    class Bar:
        # Class that considers itself greater than any instance of Number except
        # Infinity, and relies on the NotImplemented singleton for symmetric
        # relations.

        def __eq__(self, other):
            if isinstance(other, Number):
                return False
            return NotImplemented

        def __lt__(self, other):
            if other is oo:
                return True
            if isinstance(other, Number):
                return False
            return NotImplemented

        def __le__(self, other):
            return self < other or self == other

        def __gt__(self, other):
            return not self <= other

        def __ge__(self, other):
            return not self < other

    bar = Bar()

    for n in ni, nf, nr:
        assert n != bar
        assert bar != n
        assert not n == bar
        assert not bar == n
        assert n < bar
        assert bar > n
        assert not n > bar
        assert not bar < n
        assert n <= bar
        assert bar >= n
        assert not n >= bar
        assert not bar <= n

    assert oo != bar
    assert bar != oo
    assert not oo == bar
    assert not bar == oo
    assert not oo < bar
    assert not bar > oo
    assert oo > bar
    assert bar < oo
    assert not oo <= bar
    assert not bar >= oo
    assert oo >= bar
    assert bar <= oo


def test_sympyissue_24543():
    assert Rational('0.5', '100') == Rational(1, 200)
