import pytest

from diofant import (E, I, O, Rational, Symbol, acosh, acoth, asinh, atanh,
                     cos, cosh, cot, coth, csch, exp, log, nan, oo, pi, sec,
                     sech, sin, sinh, sqrt, symbols, tan, tanh, zoo)
from diofant.abc import x, y
from diofant.core.function import ArgumentIndexError
from diofant.functions.elementary.hyperbolic import \
    ReciprocalHyperbolicFunction


__all__ = ()


def test_sinh():
    k = Symbol('k', integer=True)

    assert sinh(nan) == nan
    assert sinh(zoo) == nan

    assert sinh(oo) == oo
    assert sinh(-oo) == -oo

    assert sinh(0) == 0

    assert sinh(1) == sinh(1)
    assert sinh(-1) == -sinh(1)

    assert sinh(x) == sinh(x)
    assert sinh(-x) == -sinh(x)

    assert sinh(pi) == sinh(pi)
    assert sinh(-pi) == -sinh(pi)

    assert sinh(2**1024 * E) == sinh(2**1024 * E)
    assert sinh(-2**1024 * E) == -sinh(2**1024 * E)

    assert sinh(pi*I) == 0
    assert sinh(-pi*I) == 0
    assert sinh(2*pi*I) == 0
    assert sinh(-2*pi*I) == 0
    assert sinh(-3*10**73*pi*I) == 0
    assert sinh(7*10**103*pi*I) == 0

    assert sinh(pi*I/2) == I
    assert sinh(-pi*I/2) == -I
    assert sinh(5*pi*I/2) == I
    assert sinh(7*pi*I/2) == -I

    assert sinh(pi*I/3) == sqrt(3)*I/2
    assert sinh(-2*pi*I/3) == -sqrt(3)*I/2

    assert sinh(pi*I/4) == sqrt(2)*I/2
    assert sinh(-pi*I/4) == -sqrt(2)*I/2
    assert sinh(17*pi*I/4) == sqrt(2)*I/2
    assert sinh(-3*pi*I/4) == -sqrt(2)*I/2

    assert sinh(pi*I/6) == I/2
    assert sinh(-pi*I/6) == -I/2
    assert sinh(7*pi*I/6) == -I/2
    assert sinh(-5*pi*I/6) == -I/2

    assert sinh(pi*I/105) == sin(pi/105)*I
    assert sinh(-pi*I/105) == -sin(pi/105)*I

    assert sinh(2 + 3*I) == sinh(2 + 3*I)

    assert sinh(x*I) == sin(x)*I

    assert sinh(k*pi*I) == 0
    assert sinh(17*k*pi*I) == 0

    assert sinh(k*pi*I/2) == sin(k*pi/2)*I

    r = Symbol('r', extended_real=True)
    assert sinh(r).is_extended_real
    assert sinh(x).is_extended_real is None
    i = Symbol('i', imaginary=True)
    assert sinh(i).is_finite
    assert sinh(x).is_finite is None

    pytest.raises(ArgumentIndexError, lambda: sinh(x).fdiff(2))

    a, b = symbols('a b', extended_real=True)
    z = a + b*I
    for deep in [True, False]:
        assert sinh(z).as_real_imag(deep=deep) == (sinh(a)*cos(b),
                                                   cosh(a)*sin(b))
        assert sinh(a).as_real_imag(deep=deep) == (sinh(a), 0)


def test_sinh_series():
    assert sinh(x).series(x, 0, 10) == \
        x + x**3/6 + x**5/120 + x**7/5040 + x**9/362880 + O(x**10)

    assert sinh(x).taylor_term(1, x) == x
    assert sinh(x).taylor_term(3, x) == x**3/6
    assert sinh(x).taylor_term(3, x, *(x, 0)) == sinh(x).taylor_term(3, x)


def test_cosh():
    k = Symbol('k', integer=True)

    assert cosh(nan) == nan
    assert cosh(zoo) == nan

    assert cosh(oo) == oo
    assert cosh(-oo) == oo

    assert cosh(0) == 1

    assert cosh(1) == cosh(1)
    assert cosh(-1) == cosh(1)

    assert cosh(x) == cosh(x)
    assert cosh(-x) == cosh(x)

    assert cosh(pi*I) == cos(pi)
    assert cosh(-pi*I) == cos(pi)

    assert cosh(2**1024 * E) == cosh(2**1024 * E)
    assert cosh(-2**1024 * E) == cosh(2**1024 * E)

    assert cosh(pi*I/2) == 0
    assert cosh(-pi*I/2) == 0
    assert cosh((-3*10**73 + 1)*pi*I/2) == 0
    assert cosh((7*10**103 + 1)*pi*I/2) == 0

    assert cosh(pi*I) == -1
    assert cosh(-pi*I) == -1
    assert cosh(5*pi*I) == -1
    assert cosh(8*pi*I) == 1

    assert cosh(pi*I/3) == Rational(1, 2)
    assert cosh(-2*pi*I/3) == Rational(-1, 2)

    assert cosh(pi*I/4) == sqrt(2)/2
    assert cosh(-pi*I/4) == sqrt(2)/2
    assert cosh(11*pi*I/4) == -sqrt(2)/2
    assert cosh(-3*pi*I/4) == -sqrt(2)/2

    assert cosh(pi*I/6) == sqrt(3)/2
    assert cosh(-pi*I/6) == sqrt(3)/2
    assert cosh(7*pi*I/6) == -sqrt(3)/2
    assert cosh(-5*pi*I/6) == -sqrt(3)/2

    assert cosh(pi*I/105) == cos(pi/105)
    assert cosh(-pi*I/105) == cos(pi/105)

    assert cosh(2 + 3*I) == cosh(2 + 3*I)

    assert cosh(x*I) == cos(x)

    assert cosh(k*pi*I) == cos(k*pi)
    assert cosh(17*k*pi*I) == cos(17*k*pi)

    assert cosh(k*pi) == cosh(k*pi)

    r = Symbol('r', extended_real=True)
    assert cosh(r).is_extended_real
    assert cosh(x).is_extended_real is None
    i = Symbol('i', imaginary=True)
    assert cosh(i).is_finite
    assert cosh(x).is_finite is None

    pytest.raises(ArgumentIndexError, lambda: cosh(x).fdiff(2))

    a, b = symbols('a b', extended_real=True)
    z = a + b*I
    for deep in [True, False]:
        assert cosh(z).as_real_imag(deep=deep) == (cosh(a)*cos(b),
                                                   sinh(a)*sin(b))
        assert cosh(a).as_real_imag(deep=deep) == (cosh(a), 0)


def test_cosh_series():
    assert cosh(x).series(x, 0, 10) == \
        1 + x**2/2 + x**4/24 + x**6/720 + x**8/40320 + O(x**10)

    assert cosh(x).taylor_term(0, x) == 1
    assert cosh(x).taylor_term(2, x) == x**2/2
    assert cosh(x).taylor_term(2, x, *(1, 0)) == cosh(x).taylor_term(2, x)


def test_tanh():
    k = Symbol('k', integer=True)

    assert tanh(nan) == nan
    assert tanh(zoo) == nan

    assert tanh(oo) == 1
    assert tanh(-oo) == -1

    assert tanh(0) == 0

    assert tanh(1) == tanh(1)
    assert tanh(-1) == -tanh(1)

    assert tanh(x) == tanh(x)
    assert tanh(-x) == -tanh(x)

    assert tanh(pi) == tanh(pi)
    assert tanh(-pi) == -tanh(pi)

    assert tanh(2**1024 * E) == tanh(2**1024 * E)
    assert tanh(-2**1024 * E) == -tanh(2**1024 * E)

    assert tanh(pi*I) == 0
    assert tanh(-pi*I) == 0
    assert tanh(2*pi*I) == 0
    assert tanh(-2*pi*I) == 0
    assert tanh(-3*10**73*pi*I) == 0
    assert tanh(7*10**103*pi*I) == 0

    assert tanh(pi*I/2) == tanh(pi*I/2)
    assert tanh(-pi*I/2) == -tanh(pi*I/2)
    assert tanh(5*pi*I/2) == tanh(5*pi*I/2)
    assert tanh(7*pi*I/2) == tanh(7*pi*I/2)

    assert tanh(pi*I/3) == sqrt(3)*I
    assert tanh(-2*pi*I/3) == sqrt(3)*I

    assert tanh(pi*I/4) == I
    assert tanh(-pi*I/4) == -I
    assert tanh(17*pi*I/4) == I
    assert tanh(-3*pi*I/4) == I

    assert tanh(pi*I/6) == I/sqrt(3)
    assert tanh(-pi*I/6) == -I/sqrt(3)
    assert tanh(7*pi*I/6) == I/sqrt(3)
    assert tanh(-5*pi*I/6) == I/sqrt(3)

    assert tanh(pi*I/105) == tan(pi/105)*I
    assert tanh(-pi*I/105) == -tan(pi/105)*I

    assert tanh(2 + 3*I) == tanh(2 + 3*I)

    assert tanh(x*I) == tan(x)*I

    assert tanh(k*pi*I) == 0
    assert tanh(17*k*pi*I) == 0

    assert tanh(k*pi*I/2) == tan(k*pi/2)*I

    r = Symbol('r', extended_real=True)
    assert tanh(r).is_extended_real
    assert tanh(x).is_extended_real is None
    assert tanh(r).is_finite
    assert tanh(x).is_finite is None

    pytest.raises(ArgumentIndexError, lambda: tanh(x).fdiff(2))

    a, b = symbols('a b', extended_real=True)
    z = a + b*I
    for deep in [True, False]:
        d = sinh(a)**2 + cos(b)**2
        assert tanh(z).as_real_imag(deep=deep) == (sinh(a)*cosh(a)/d,
                                                   sin(b)*cos(b)/d)
        assert tanh(a).as_real_imag(deep=deep) == (tanh(a), 0)


def test_tanh_series():
    assert tanh(x).series(x, 0, 10) == \
        x - x**3/3 + 2*x**5/15 - 17*x**7/315 + 62*x**9/2835 + O(x**10)


def test_coth():
    k = Symbol('k', integer=True)

    assert coth(nan) == nan
    assert coth(zoo) == nan

    assert coth(oo) == 1
    assert coth(-oo) == -1

    assert coth(0) == coth(0)
    assert coth(0) == zoo
    assert coth(1) == coth(1)
    assert coth(-1) == -coth(1)

    assert coth(x) == coth(x)
    assert coth(-x) == -coth(x)

    assert coth(pi*I) == -I*cot(pi)
    assert coth(-pi*I) == cot(pi)*I

    assert coth(2**1024 * E) == coth(2**1024 * E)
    assert coth(-2**1024 * E) == -coth(2**1024 * E)

    assert coth(pi*I) == -I*cot(pi)
    assert coth(-pi*I) == I*cot(pi)
    assert coth(2*pi*I) == -I*cot(2*pi)
    assert coth(-2*pi*I) == I*cot(2*pi)
    assert coth(-3*10**73*pi*I) == I*cot(3*10**73*pi)
    assert coth(7*10**103*pi*I) == -I*cot(7*10**103*pi)

    assert coth(pi*I/2) == 0
    assert coth(-pi*I/2) == 0
    assert coth(5*pi*I/2) == 0
    assert coth(7*pi*I/2) == 0

    assert coth(pi*I/3) == -I/sqrt(3)
    assert coth(-2*pi*I/3) == -I/sqrt(3)

    assert coth(pi*I/4) == -I
    assert coth(-pi*I/4) == I
    assert coth(17*pi*I/4) == -I
    assert coth(-3*pi*I/4) == -I

    assert coth(pi*I/6) == -sqrt(3)*I
    assert coth(-pi*I/6) == sqrt(3)*I
    assert coth(7*pi*I/6) == -sqrt(3)*I
    assert coth(-5*pi*I/6) == -sqrt(3)*I

    assert coth(pi*I/105) == -cot(pi/105)*I
    assert coth(-pi*I/105) == cot(pi/105)*I

    assert coth(2 + 3*I) == coth(2 + 3*I)

    assert coth(x*I) == -cot(x)*I

    assert coth(k*pi*I) == -cot(k*pi)*I
    assert coth(17*k*pi*I) == -cot(17*k*pi)*I

    assert coth(k*pi*I) == -cot(k*pi)*I

    pytest.raises(ArgumentIndexError, lambda: coth(x).fdiff(2))

    a, b = symbols('a b', extended_real=True)
    z = a + b*I
    for deep in [True, False]:
        d = sinh(a)**2 + sin(b)**2
        assert coth(z).as_real_imag(deep=deep) == (sinh(a)*cosh(a)/d,
                                                   -sin(b)*cos(b)/d)
        assert coth(a).as_real_imag(deep=deep) == (coth(a), 0)


def test_coth_series():
    assert coth(x).series(x, 0, 8) == \
        1/x + x/3 - x**3/45 + 2*x**5/945 - x**7/4725 + O(x**8)


def test_csch():
    k = Symbol('k', integer=True)
    n = Symbol('n', positive=True)

    assert csch(nan) == nan
    assert csch(zoo) == nan

    assert csch(oo) == 0
    assert csch(-oo) == 0

    assert csch(0) == zoo

    assert csch(-1) == -csch(1)

    assert csch(-x) == -csch(x)
    assert csch(-pi) == -csch(pi)
    assert csch(-2**1024 * E) == -csch(2**1024 * E)

    assert csch(pi*I) == zoo
    assert csch(-pi*I) == zoo
    assert csch(2*pi*I) == zoo
    assert csch(-2*pi*I) == zoo
    assert csch(-3*10**73*pi*I) == zoo
    assert csch(7*10**103*pi*I) == zoo

    assert csch(pi*I/2) == -I
    assert csch(-pi*I/2) == I
    assert csch(5*pi*I/2) == -I
    assert csch(7*pi*I/2) == I

    assert csch(pi*I/3) == -2/sqrt(3)*I
    assert csch(-2*pi*I/3) == 2/sqrt(3)*I

    assert csch(pi*I/4) == -sqrt(2)*I
    assert csch(-pi*I/4) == sqrt(2)*I
    assert csch(7*pi*I/4) == sqrt(2)*I
    assert csch(-3*pi*I/4) == sqrt(2)*I

    assert csch(pi*I/6) == -2*I
    assert csch(-pi*I/6) == 2*I
    assert csch(7*pi*I/6) == 2*I
    assert csch(-7*pi*I/6) == -2*I
    assert csch(-5*pi*I/6) == 2*I

    assert csch(pi*I/105) == -1/sin(pi/105)*I
    assert csch(-pi*I/105) == 1/sin(pi/105)*I

    assert csch(x*I) == -1/sin(x)*I

    assert csch(k*pi*I) == zoo
    assert csch(17*k*pi*I) == zoo

    assert csch(k*pi*I/2) == -1/sin(k*pi/2)*I

    assert csch(n).is_extended_real is True
    assert csch(n).is_finite is None

    pytest.raises(ArgumentIndexError, lambda: csch(x).fdiff(2))


def test_csch_series():
    assert (csch(x).series(x, 0, 10) ==
            1/x - x/6 + 7*x**3/360 - 31*x**5/15120 + 127*x**7/604800 -
            73*x**9/3421440 + O(x**10))


def test_reciprocal():
    class FakeSech(ReciprocalHyperbolicFunction):
        _reciprocal_of = cosh

    assert FakeSech(-x) == 1/cosh(x)


def test_sech():
    k = Symbol('k', integer=True)
    n = Symbol('n', positive=True)

    assert sech(nan) == nan
    assert sech(zoo) == nan

    assert sech(oo) == 0
    assert sech(-oo) == 0

    assert sech(0) == 1

    assert sech(-1) == sech(1)
    assert sech(-x) == sech(x)

    assert sech(pi*I) == sec(pi)

    assert sech(-pi*I) == sec(pi)
    assert sech(-2**1024 * E) == sech(2**1024 * E)

    assert sech(pi*I/2) == zoo
    assert sech(-pi*I/2) == zoo
    assert sech((-3*10**73 + 1)*pi*I/2) == zoo
    assert sech((7*10**103 + 1)*pi*I/2) == zoo

    assert sech(pi*I) == -1
    assert sech(-pi*I) == -1
    assert sech(5*pi*I) == -1
    assert sech(8*pi*I) == 1

    assert sech(pi*I/3) == 2
    assert sech(-2*pi*I/3) == -2

    assert sech(pi*I/4) == sqrt(2)
    assert sech(-pi*I/4) == sqrt(2)
    assert sech(5*pi*I/4) == -sqrt(2)
    assert sech(-5*pi*I/4) == -sqrt(2)

    assert sech(pi*I/6) == 2/sqrt(3)
    assert sech(-pi*I/6) == 2/sqrt(3)
    assert sech(7*pi*I/6) == -2/sqrt(3)
    assert sech(-5*pi*I/6) == -2/sqrt(3)

    assert sech(pi*I/105) == 1/cos(pi/105)
    assert sech(-pi*I/105) == 1/cos(pi/105)

    assert sech(x*I) == 1/cos(x)

    assert sech(k*pi*I) == 1/cos(k*pi)
    assert sech(17*k*pi*I) == 1/cos(17*k*pi)

    assert sech(n).is_extended_real is True
    assert csch(n).is_finite is None

    pytest.raises(ArgumentIndexError, lambda: sech(x).fdiff(2))


def test_sech_series():
    assert sech(x).series(x, 0, 10) == \
        1 - x**2/2 + 5*x**4/24 - 61*x**6/720 + 277*x**8/8064 + O(x**10)


def test_asinh():
    assert asinh(x) == asinh(x)
    assert asinh(-x) == -asinh(x)
    assert asinh(-2) == -asinh(2)
    assert asinh(nan) == nan
    assert asinh( 0) == 0
    assert asinh(+1) == log(sqrt(2) + 1)

    assert asinh(-1) == log(sqrt(2) - 1)
    assert asinh(I) == pi*I/2
    assert asinh(-I) == -pi*I/2
    assert asinh(I/2) == pi*I/6
    assert asinh(-I/2) == -pi*I/6

    assert asinh(oo) == oo
    assert asinh(-oo) == -oo

    assert asinh(I*oo) == oo
    assert asinh(-I*oo) == -oo

    assert asinh(zoo) == zoo

    assert asinh(I*(sqrt(3) - 1)/2**Rational(3, 2)) == pi*I/12
    assert asinh(-I*(sqrt(3) - 1)/2**Rational(3, 2)) == -pi*I/12

    assert asinh(I*(sqrt(5) - 1)/4) == pi*I/10
    assert asinh(-I*(sqrt(5) - 1)/4) == -pi*I/10

    assert asinh(I*(sqrt(5) + 1)/4) == 3*pi*I/10
    assert asinh(-I*(sqrt(5) + 1)/4) == -3*pi*I/10

    pytest.raises(ArgumentIndexError, lambda: asinh(x).fdiff(2))


def test_asinh_rewrite():
    assert asinh(x).rewrite(log) == log(x + sqrt(x**2 + 1))


def test_asinh_series():
    assert asinh(x).series(x, 0, 8) == \
        x - x**3/6 + 3*x**5/40 - 5*x**7/112 + O(x**8)
    t5 = asinh(x).taylor_term(5, x)
    assert t5 == 3*x**5/40
    assert asinh(x).taylor_term(7, x, t5, 0) == -5*x**7/112


def test_acosh():
    # TODO please write more tests  -- see issue sympy/sympy#3751
    # From http://functions.wolfram.com/ElementaryFunctions/ArcCosh/03/01/
    # at specific points

    assert acosh(-x) == acosh(-x)

    assert acosh(1) == 0
    assert acosh(-1) == pi*I
    assert acosh(0) == I*pi/2
    assert acosh(Rational(1, 2)) == I*pi/3
    assert acosh(Rational(-1, 2)) == 2*pi*I/3

    assert acosh(+oo) == oo
    assert acosh(-oo) == oo
    assert acosh(+I*oo) == oo
    assert acosh(-I*oo) == oo

    assert acosh(zoo) == oo

    assert acosh(I) == log(I*(1 + sqrt(2)))
    assert acosh(-I) == log(-I*(1 + sqrt(2)))
    assert acosh((sqrt(3) - 1)/(2*sqrt(2))) == 5*pi*I/12
    assert acosh(-(sqrt(3) - 1)/(2*sqrt(2))) == 7*pi*I/12
    assert acosh(sqrt(2)/2) == I*pi/4
    assert acosh(-sqrt(2)/2) == 3*I*pi/4
    assert acosh(sqrt(3)/2) == I*pi/6
    assert acosh(-sqrt(3)/2) == 5*I*pi/6
    assert acosh(sqrt(2 + sqrt(2))/2) == I*pi/8
    assert acosh(-sqrt(2 + sqrt(2))/2) == 7*I*pi/8
    assert acosh(sqrt(2 - sqrt(2))/2) == 3*I*pi/8
    assert acosh(-sqrt(2 - sqrt(2))/2) == 5*I*pi/8
    assert acosh((1 + sqrt(3))/(2*sqrt(2))) == I*pi/12
    assert acosh(-(1 + sqrt(3))/(2*sqrt(2))) == 11*I*pi/12
    assert acosh((sqrt(5) + 1)/4) == I*pi/5
    assert acosh(-(sqrt(5) + 1)/4) == 4*I*pi/5

    assert str(acosh(5*I).evalf(6)) == '2.31244 + 1.5708*I'
    assert str(acosh(-5*I).evalf(6)) == '2.31244 - 1.5708*I'

    pytest.raises(ArgumentIndexError, lambda: acosh(x).fdiff(2))


def test_acosh_series():
    assert acosh(x).series(x, 0, 8) == \
        -I*x + pi*I/2 - I*x**3/6 - 3*I*x**5/40 - 5*I*x**7/112 + O(x**8)
    t5 = acosh(x).taylor_term(5, x)
    assert t5 == - 3*I*x**5/40
    assert acosh(x).taylor_term(7, x, t5, 0) == - 5*I*x**7/112


# TODO please write more tests -- see issue sympy/sympy#3751
def test_atanh():
    # TODO please write more tests  -- see issue sympy/sympy#3751
    # From http://functions.wolfram.com/ElementaryFunctions/ArcTanh/03/01/
    # at specific points

    # at specific points
    assert atanh(0) == 0
    assert atanh(I) == I*pi/4
    assert atanh(-I) == -I*pi/4
    assert atanh(1) == oo
    assert atanh(-1) == -oo
    assert atanh(-2) == -atanh(2)

    # at infinites
    assert atanh(I*oo) == I*pi/2
    assert atanh(-I*oo) == -I*pi/2

    assert atanh(zoo) == nan

    # properties
    assert atanh(-x) == -atanh(x)

    assert atanh(I/sqrt(3)) == I*pi/6
    assert atanh(-I/sqrt(3)) == -I*pi/6
    assert atanh(I*sqrt(3)) == I*pi/3
    assert atanh(-I*sqrt(3)) == -I*pi/3
    assert atanh(I*(1 + sqrt(2))) == 3*pi*I/8
    assert atanh(I*(sqrt(2) - 1)) == pi*I/8
    assert atanh(I*(1 - sqrt(2))) == -pi*I/8
    assert atanh(-I*(1 + sqrt(2))) == -3*pi*I/8
    assert atanh(I*sqrt(5 + 2*sqrt(5))) == 2*I*pi/5
    assert atanh(-I*sqrt(5 + 2*sqrt(5))) == -2*I*pi/5
    assert atanh(I*(2 - sqrt(3))) == pi*I/12
    assert atanh(I*(sqrt(3) - 2)) == -pi*I/12
    assert atanh(oo) == -I*pi/2
    assert atanh(-oo) == I*pi/2

    pytest.raises(ArgumentIndexError, lambda: atanh(x).fdiff(2))


def test_atanh_series():
    assert atanh(x).series(x, 0, 10) == \
        x + x**3/3 + x**5/5 + x**7/7 + x**9/9 + O(x**10)


def test_acoth():
    # TODO please write more tests  -- see issue sympy/sympy#3751
    # From http://functions.wolfram.com/ElementaryFunctions/ArcCoth/03/01/
    # at specific points

    # at specific points
    assert acoth(0) == I*pi/2
    assert acoth(I) == -I*pi/4
    assert acoth(-I) == I*pi/4
    assert acoth(1) == oo
    assert acoth(-1) == -oo
    assert acoth(-2) == -acoth(2)

    # at infinites
    assert acoth(oo) == 0
    assert acoth(-oo) == 0
    assert acoth(I*oo) == 0
    assert acoth(-I*oo) == 0
    assert acoth(zoo) == 0

    # properties
    assert acoth(-x) == -acoth(x)

    assert acoth(I/sqrt(3)) == -I*pi/3
    assert acoth(-I/sqrt(3)) == I*pi/3
    assert acoth(I*sqrt(3)) == -I*pi/6
    assert acoth(-I*sqrt(3)) == I*pi/6
    assert acoth(I*(1 + sqrt(2))) == -pi*I/8
    assert acoth(-I*(sqrt(2) + 1)) == pi*I/8
    assert acoth(I*(1 - sqrt(2))) == 3*pi*I/8
    assert acoth(I*(sqrt(2) - 1)) == -3*pi*I/8
    assert acoth(I*sqrt(5 + 2*sqrt(5))) == -I*pi/10
    assert acoth(-I*sqrt(5 + 2*sqrt(5))) == I*pi/10
    assert acoth(I*(2 + sqrt(3))) == -pi*I/12
    assert acoth(-I*(2 + sqrt(3))) == pi*I/12
    assert acoth(I*(2 - sqrt(3))) == -5*pi*I/12
    assert acoth(I*(sqrt(3) - 2)) == 5*pi*I/12

    pytest.raises(ArgumentIndexError, lambda: acoth(x).fdiff(2))


def test_acoth_series():
    assert acoth(x).series(x, 0, 10) == \
        I*pi/2 + x + x**3/3 + x**5/5 + x**7/7 + x**9/9 + O(x**10)


def test_inverses():
    assert sinh(x).inverse() == asinh
    pytest.raises(AttributeError, lambda: cosh(x).inverse())
    assert tanh(x).inverse() == atanh
    assert coth(x).inverse() == acoth
    assert asinh(x).inverse() == sinh
    assert acosh(x).inverse() == cosh
    assert atanh(x).inverse() == tanh
    assert acoth(x).inverse() == coth


def test_leading_term():
    assert cosh(x).as_leading_term(x) == 1
    assert coth(x).as_leading_term(x) == 1/x
    assert acosh(x).as_leading_term(x) == I*pi/2
    assert acoth(x).as_leading_term(x) == I*pi/2
    for func in [sinh, tanh, asinh, atanh]:
        assert func(x).as_leading_term(x) == x
    for func in [sinh, cosh, tanh, coth, asinh, acosh, atanh, acoth]:
        for arg in (1/x, Rational(1, 2)):
            eq = func(arg)
            assert eq.as_leading_term(x) == eq
    for func in [csch, sech]:
        eq = func(Rational(1, 2))
        assert eq.as_leading_term(x) == eq
    assert csch(x).as_leading_term(x) == 1/x


def test_complex():
    a, b = symbols('a,b', extended_real=True)
    z = a + b*I
    for func in [sinh, cosh, tanh, coth, sech, csch]:
        assert func(z).conjugate() == func(a - b*I)
    for deep in [True, False]:
        assert sinh(z).expand(
            complex=True, deep=deep) == sinh(a)*cos(b) + I*cosh(a)*sin(b)
        assert cosh(z).expand(
            complex=True, deep=deep) == cosh(a)*cos(b) + I*sinh(a)*sin(b)
        assert tanh(z).expand(complex=True, deep=deep) == sinh(a)*cosh(
            a)/(cos(b)**2 + sinh(a)**2) + I*sin(b)*cos(b)/(cos(b)**2 + sinh(a)**2)
        assert coth(z).expand(complex=True, deep=deep) == sinh(a)*cosh(
            a)/(sin(b)**2 + sinh(a)**2) - I*sin(b)*cos(b)/(sin(b)**2 + sinh(a)**2)
        assert csch(z).expand(complex=True, deep=deep) == cos(b) * sinh(a) / (sin(b)**2
                                                                              * cosh(a)**2 + cos(b)**2 * sinh(a)**2) - I*sin(b) * cosh(a) / (sin(b)**2
                                                                                                                                             * cosh(a)**2 + cos(b)**2 * sinh(a)**2)
        assert sech(z).expand(complex=True, deep=deep) == cos(b) * cosh(a) / (sin(b)**2
                                                                              * sinh(a)**2 + cos(b)**2 * cosh(a)**2) - I*sin(b) * sinh(a) / (sin(b)**2
                                                                                                                                             * sinh(a)**2 + cos(b)**2 * cosh(a)**2)


def test_complex_2899():
    a, b = symbols('a,b', extended_real=True)
    for deep in [True, False]:
        for func in [sinh, cosh, tanh, coth]:
            assert func(a).expand(complex=True, deep=deep) == func(a)


def test_simplifications():
    assert sinh(asinh(x)) == x
    assert sinh(acosh(x)) == sqrt(x - 1) * sqrt(x + 1)
    assert sinh(atanh(x)) == x/sqrt(1 - x**2)
    assert sinh(acoth(x)) == 1/(sqrt(x - 1) * sqrt(x + 1))

    assert cosh(asinh(x)) == sqrt(1 + x**2)
    assert cosh(acosh(x)) == x
    assert cosh(atanh(x)) == 1/sqrt(1 - x**2)
    assert cosh(acoth(x)) == x/(sqrt(x - 1) * sqrt(x + 1))

    assert tanh(asinh(x)) == x/sqrt(1 + x**2)
    assert tanh(acosh(x)) == sqrt(x - 1) * sqrt(x + 1) / x
    assert tanh(atanh(x)) == x
    assert tanh(acoth(x)) == 1/x

    assert coth(asinh(x)) == sqrt(1 + x**2)/x
    assert coth(acosh(x)) == x/(sqrt(x - 1) * sqrt(x + 1))
    assert coth(atanh(x)) == 1/x
    assert coth(acoth(x)) == x

    assert csch(asinh(x)) == 1/x
    assert csch(acosh(x)) == 1/(sqrt(x - 1) * sqrt(x + 1))
    assert csch(atanh(x)) == sqrt(1 - x**2)/x
    assert csch(acoth(x)) == sqrt(x - 1) * sqrt(x + 1)

    assert sech(asinh(x)) == 1/sqrt(1 + x**2)
    assert sech(acosh(x)) == 1/x
    assert sech(atanh(x)) == sqrt(1 - x**2)
    assert sech(acoth(x)) == sqrt(x - 1) * sqrt(x + 1)/x


def test_sympyissue_4136():
    assert cosh(asinh(Rational(3, 2))) == sqrt(Rational(13, 4))


def test_sinh_rewrite():
    assert sinh(x).rewrite(exp) == (exp(x) - exp(-x))/2 \
        == sinh(x).rewrite('tractable')
    assert sinh(x).rewrite(cosh) == -I*cosh(x + I*pi/2)
    assert sinh(x).rewrite(tanh) == 2*tanh(x/2)/(1 - tanh(x/2)**2)
    assert sinh(x).rewrite(coth) == 2*coth(x/2)/(coth(x/2)**2 - 1)


def test_cosh_rewrite():
    assert cosh(x).rewrite(exp) == (exp(x) + exp(-x))/2 \
        == cosh(x).rewrite('tractable')
    assert cosh(x).rewrite(sinh) == -I*sinh(x + I*pi/2)
    assert cosh(x).rewrite(tanh) == (1 + tanh(x/2)**2)/(1 - tanh(x/2)**2)
    assert cosh(x).rewrite(coth) == (coth(x/2)**2 + 1)/(coth(x/2)**2 - 1)


def test_tanh_rewrite():
    assert tanh(x).rewrite(exp) == (exp(x) - exp(-x))/(exp(x) + exp(-x)) \
        == tanh(x).rewrite('tractable')
    assert tanh(x).rewrite(sinh) == I*sinh(x)/sinh(I*pi/2 - x)
    assert tanh(x).rewrite(cosh) == I*cosh(I*pi/2 - x)/cosh(x)
    assert tanh(x).rewrite(coth) == 1/coth(x)


def test_coth_rewrite():
    assert coth(x).rewrite(exp) == (exp(x) + exp(-x))/(exp(x) - exp(-x)) \
        == coth(x).rewrite('tractable')
    assert coth(x).rewrite(sinh) == -I*sinh(I*pi/2 - x)/sinh(x)
    assert coth(x).rewrite(cosh) == -I*cosh(x)/cosh(I*pi/2 - x)
    assert coth(x).rewrite(tanh) == 1/tanh(x)


def test_csch_rewrite():
    assert csch(x).rewrite(exp) == 1 / (exp(x)/2 - exp(-x)/2) \
        == csch(x).rewrite('tractable')
    assert csch(x).rewrite(cosh) == I/cosh(x + I*pi/2)
    assert csch(x).rewrite(tanh) == (1 - tanh(x/2)**2)/(2*tanh(x/2))
    assert csch(x).rewrite(coth) == (coth(x/2)**2 - 1)/(2*coth(x/2))


def test_sech_rewrite():
    assert sech(x).rewrite(exp) == 1 / (exp(x)/2 + exp(-x)/2) \
        == sech(x).rewrite('tractable')
    assert sech(x).rewrite(sinh) == I/sinh(x + I*pi/2)
    assert sech(x).rewrite(tanh) == (1 - tanh(x/2)**2)/(1 + tanh(x/2)**2)
    assert sech(x).rewrite(coth) == (coth(x/2)**2 - 1)/(coth(x/2)**2 + 1)


def test_derivs():
    assert coth(x).diff(x) == -sinh(x)**(-2)
    assert sinh(x).diff(x) == cosh(x)
    assert cosh(x).diff(x) == sinh(x)
    assert tanh(x).diff(x) == -tanh(x)**2 + 1
    assert csch(x).diff(x) == -coth(x)*csch(x)
    assert sech(x).diff(x) == -tanh(x)*sech(x)
    assert acoth(x).diff(x) == 1/(-x**2 + 1)
    assert asinh(x).diff(x) == 1/sqrt(x**2 + 1)
    assert acosh(x).diff(x) == 1/sqrt(x**2 - 1)
    assert atanh(x).diff(x) == 1/(-x**2 + 1)


def test_sinh_expansion():
    assert sinh(x+y).expand(trig=True) == sinh(x)*cosh(y) + cosh(x)*sinh(y)
    assert sinh(2*x).expand(trig=True) == 2*sinh(x)*cosh(x)
    assert sinh(3*x).expand(trig=True).expand() == \
        sinh(x)**3 + 3*sinh(x)*cosh(x)**2


def test_cosh_expansion():
    assert cosh(x+y).expand(trig=True) == cosh(x)*cosh(y) + sinh(x)*sinh(y)
    assert cosh(2*x).expand(trig=True) == cosh(x)**2 + sinh(x)**2
    assert cosh(3*x).expand(trig=True).expand() == \
        3*sinh(x)**2*cosh(x) + cosh(x)**3
