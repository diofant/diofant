import pytest

from diofant import (Basic, Float, I, Integer, Number, O, Pow, Rational,
                     Symbol, cbrt, cos, exp, log, nan, oo, pi, root, sin, sqrt,
                     symbols, true, zoo)
from diofant.abc import a, b, c, x, y
from diofant.tests.core.test_evalf import NS


__all__ = ()


def test_rational():
    a = Rational(1, 5)

    r = sqrt(5)/5
    assert sqrt(a) == r
    assert 2*sqrt(a) == 2*r

    r = a*sqrt(a)
    assert a**Rational(3, 2) == r
    assert 2*a**Rational(3, 2) == 2*r

    r = a**5*a**Rational(2, 3)
    assert a**Rational(17, 3) == r
    assert 2 * a**Rational(17, 3) == 2*r


def test_large_rational():
    e = cbrt(Rational(123712**12 - 1, 7) + Rational(1, 7))
    assert e == 234232585392159195136*cbrt(Rational(1, 7))


def test_negative_real():
    def feq(a, b):
        return abs(a - b) < 1E-10

    assert feq(1/Float(-0.5), -Integer(2))


def test_expand():
    assert (2**(-1 - x)).expand() == Rational(1, 2)*2**(-x)


def test_sqrt():
    # issue sympy/sympy#7638
    e = 1 + I/5
    assert sqrt(e**5) == e**Rational(5, 2)
    assert sqrt(e**6) == e**3
    r = symbols('r', extended_real=True)
    assert sqrt((1 + I*r)**6) != (1 + I*r)**3


def test_sympyissue_3449():
    # test if powers are simplified correctly, see also sympy/sympy#3995
    assert cbrt(x)**2 == x**Rational(2, 3)
    assert (x**3)**Rational(2, 5) == Pow(x**3, Rational(2, 5), evaluate=False)

    a = Symbol('a', extended_real=True)
    b = Symbol('b', extended_real=True)
    assert (a**2)**b == (abs(a)**b)**2
    assert sqrt(1/a) != 1/sqrt(a)  # e.g. for a = -1
    assert cbrt(a**3) != a
    assert (x**a)**b != x**(a*b)  # e.g. x = -1, a=2, b=1/2
    assert (x**.5)**b == x**(.5*b)
    assert (x**.5)**.5 == x**.25
    assert (x**2.5)**.5 != x**1.25  # e.g. for x = 5*I

    k = Symbol('k', integer=True)
    m = Symbol('m', integer=True)
    assert (x**k)**m == x**(k*m)
    assert Number(5)**Rational(2, 3) == cbrt(Number(25))

    assert (x**.5)**2 == x**1.0
    assert (x**2)**k == (x**k)**2 == x**(2*k)

    a = Symbol('a', positive=True)
    assert (a**3)**Rational(2, 5) == a**Rational(6, 5)
    assert (a**2)**b == (a**b)**2
    assert (a**Rational(2, 3))**x == (a**(2*x/3)) != (a**x)**Rational(2, 3)

    assert sqrt(x - 1).subs({x: 5}) == 2


def test_sympyissue_3866():
    # pylint: disable=nonexistent-operator
    assert -(-sqrt(sqrt(5) - 1)) == sqrt(sqrt(5) - 1)


def test_negative_one():
    x = Symbol('x', complex=True)
    y = Symbol('y', complex=True)
    assert 1/x**y == x**(-y)


def test_sympyissue_4362():
    neg = Symbol('neg', negative=True)
    nonneg = Symbol('nonneg', nonnegative=True)
    any = Symbol('any')
    num, den = sqrt(1/neg).as_numer_denom()
    assert num == sqrt(-1)
    assert den == sqrt(-neg)
    num, den = sqrt(1/nonneg).as_numer_denom()
    assert num == 1
    assert den == sqrt(nonneg)
    num, den = sqrt(1/any).as_numer_denom()
    assert num == sqrt(1/any)
    assert den == 1

    def eqn(num, den, pow):
        return (num/den)**pow
    npos = 1
    nneg = -1
    dpos = 2 - sqrt(3)
    dneg = 1 - sqrt(3)
    assert dpos > 0 > dneg
    assert npos > 0 > nneg
    # pos or neg integer
    eq = eqn(npos, dpos, 2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (1, dpos**2)
    eq = eqn(npos, dneg, 2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (1, dneg**2)
    eq = eqn(nneg, dpos, 2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (1, dpos**2)
    eq = eqn(nneg, dneg, 2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (1, dneg**2)
    eq = eqn(npos, dpos, -2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (dpos**2, 1)
    eq = eqn(npos, dneg, -2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (dneg**2, 1)
    eq = eqn(nneg, dpos, -2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (dpos**2, 1)
    eq = eqn(nneg, dneg, -2)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (dneg**2, 1)
    # pos or neg rational
    pow = Rational(1, 2)
    eq = eqn(npos, dpos, pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (npos**pow, dpos**pow)
    eq = eqn(npos, dneg, pow)
    assert eq.is_Pow is False
    assert eq.as_numer_denom() == ((-npos)**pow, (-dneg)**pow)
    eq = eqn(nneg, dpos, pow)
    assert not eq.is_Pow or eq.as_numer_denom() == (nneg**pow, dpos**pow)
    eq = eqn(nneg, dneg, pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == ((-nneg)**pow, (-dneg)**pow)
    eq = eqn(npos, dpos, -pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (dpos**pow, npos**pow)
    eq = eqn(npos, dneg, -pow)
    assert eq.is_Pow is False
    assert eq.as_numer_denom() == (-(-npos)**pow*(-dneg)**pow, npos)
    eq = eqn(nneg, dpos, -pow)
    assert not eq.is_Pow or eq.as_numer_denom() == (dpos**pow, nneg**pow)
    eq = eqn(nneg, dneg, -pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == ((-dneg)**pow, (-nneg)**pow)
    # unknown exponent
    pow = 2*any
    eq = eqn(npos, dpos, pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (npos**pow, dpos**pow)
    eq = eqn(npos, dneg, pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == ((-npos)**pow, (-dneg)**pow)
    eq = eqn(nneg, dpos, pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (nneg**pow, dpos**pow)
    eq = eqn(nneg, dneg, pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == ((-nneg)**pow, (-dneg)**pow)
    eq = eqn(npos, dpos, -pow)
    assert eq.as_numer_denom() == (dpos**pow, npos**pow)
    eq = eqn(npos, dneg, -pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == ((-dneg)**pow, (-npos)**pow)
    eq = eqn(nneg, dpos, -pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == (dpos**pow, nneg**pow)
    eq = eqn(nneg, dneg, -pow)
    assert eq.is_Pow
    assert eq.as_numer_denom() == ((-dneg)**pow, (-nneg)**pow)

    assert ((1/(1 + x/3))**-1).as_numer_denom() == (3 + x, 3)
    notp = Symbol('notp', positive=False)  # not positive does not imply real
    b = (1 + x/notp)**-2
    assert (b**(-y)).as_numer_denom() == (1, b**y)
    assert (b**-1).as_numer_denom() == ((notp + x)**2, notp**2)
    nonp = Symbol('nonp', nonpositive=True)
    assert (((1 + x/nonp)**-2)**-1).as_numer_denom() == ((-nonp - x)**2,
                                                         nonp**2)

    n = Symbol('n', negative=True)
    assert (x**n).as_numer_denom() == (1, x**-n)
    assert sqrt(1/n).as_numer_denom() == (I, sqrt(-n))
    n = Symbol('0 or neg', nonpositive=True)
    # if x and n are split up without negating each term and n is negative
    # then the answer might be wrong; if n is 0 it won't matter since
    # 1/oo and 1/zoo are both zero as is sqrt(0)/sqrt(-x) unless x is also
    # zero (in which case the negative sign doesn't matter):
    # 1/sqrt(1/-1) = -I but sqrt(-1)/sqrt(1) = I
    assert (1/sqrt(x/n)).as_numer_denom() == (sqrt(-n), sqrt(-x))
    c = Symbol('c', complex=True)
    e = sqrt(1/c)
    assert e.as_numer_denom() == (e, 1)
    i = Symbol('i', integer=True)
    assert ((1 + x/y)**i).as_numer_denom() == ((x + y)**i, y**i)


def test_Pow_signs():
    # Cf. issues sympy/sympy#4595 and sympy/sympy#5250
    n = Symbol('n', even=True)
    assert (3 - y)**2 != (y - 3)**2
    assert (3 - y)**n != (y - 3)**n
    assert (-3 + y - x)**2 != (3 - y + x)**2
    assert (y - 3)**3 != -(3 - y)**3


def test_power_with_noncommutative_mul_as_base():
    x = Symbol('x', commutative=False)
    y = Symbol('y', commutative=False)
    assert (x*y)**3 != x**3*y**3
    assert (2*x*y)**3 == 8*(x*y)**3


def test_zero():
    assert 0**x != 0
    assert 0**(2*x) == 0**x
    assert 0**(1.0*x) == 0**x
    assert 0**(2.0*x) == 0**x
    assert (0**(2 - x)).as_base_exp() == (0, 2 - x)
    assert 0**(x - 2) != oo**(2 - x)
    assert 0**(2*x*y) == 0**(x*y)
    assert 0**(-2*x*y) == zoo**(x*y)
    assert 0**I == nan
    i = Symbol('i', imaginary=True, nonzero=True)
    assert 0**i == nan


def test_pow_as_base_exp():
    assert (oo**(2 - x)).as_base_exp() == (oo, 2 - x)
    assert (oo**(x - 2)).as_base_exp() == (oo, x - 2)
    p = Rational(1, 2)**x
    assert p.base, p.exp == p.as_base_exp() == (2, -x)
    # issue sympy/sympy#8344:
    assert Pow(1, 2, evaluate=False).as_base_exp() == (1, 2)
    # issue sympy/sympy#21396
    assert I.as_base_exp() == (-1, Rational(1, 2))
    assert sqrt(I).as_base_exp() == (-1, Rational(1, 4))


def test_sympyissue_6100():
    assert x**1.0 != x
    assert x != x**1.0
    assert true != x**1.0
    assert x**1.0 is not True
    assert x is not True
    assert x*y != (x*y)**1.0
    assert (x**1.0)**1.0 != x
    assert (x**1.0)**2.0 == x**2
    b = Basic()
    assert Pow(b, 1.0, evaluate=False) != b
    # if the following gets distributed as a Mul (x**1.0*y**1.0 then
    # __eq__ methods could be added to Symbol and Pow to detect the
    # power-of-1.0 case.
    assert isinstance((x*y)**1.0, Pow)


def test_sympyissue_6208():
    assert sqrt(33**(9*I/10)) == -33**(9*I/20)
    assert root((6*I)**(2*I), 3).as_base_exp()[1] == Rational(1, 3)  # != 2*I/3
    assert root((6*I)**(I/3), 3).as_base_exp()[1] == I/9
    assert sqrt(exp(3*I)) == exp(3*I/2)
    assert sqrt(-sqrt(3)*(1 + 2*I)) == sqrt(sqrt(3))*sqrt(-1 - 2*I)
    assert sqrt(exp(5*I)) == -exp(5*I/2)
    assert root(exp(5*I), 3).exp == Rational(1, 3)


def test_sympyissue_6990():
    assert (sqrt(a + b*x + x**2)).series(x, 0, 3).removeO() == \
        b*x/(2*sqrt(a)) + x**2*(1/(2*sqrt(a)) -
                                b**2/(8*a**Rational(3, 2))) + sqrt(a)


def test_sympyissue_6068():
    assert sqrt(sin(x)).series(x, 0, 8) == \
        sqrt(x) - x**Rational(5, 2)/12 + x**Rational(9, 2)/1440 - \
        x**Rational(13, 2)/24192 + O(x**8)
    assert sqrt(sin(x)).series(x, 0, 10) == \
        sqrt(x) - x**Rational(5, 2)/12 + x**Rational(9, 2)/1440 - \
        x**Rational(13, 2)/24192 - 67*x**Rational(17, 2)/29030400 + O(x**10)
    assert sqrt(sin(x**3)).series(x, 0, 19) == \
        x**Rational(3, 2) - x**Rational(15, 2)/12 + x**Rational(27, 2)/1440 + O(x**19)
    assert sqrt(sin(x**3)).series(x, 0, 20) == \
        x**Rational(3, 2) - x**Rational(15, 2)/12 + x**Rational(27, 2)/1440 - \
        x**Rational(39, 2)/24192 + O(x**20)


def test_sympyissue_6782():
    assert sqrt(sin(x**3)).series(x, 0, 7) == x**Rational(3, 2) + O(x**7)
    assert sqrt(sin(x**4)).series(x, 0, 3) == x**2 + O(x**3)


def test_sympyissue_6653():
    assert (1 / sqrt(1 + sin(x**2))).series(x, 0, 3) == 1 - x**2/2 + O(x**3)


def test_sympyissue_6429():
    f = (c**2 + x)**(0.5)
    assert f.series(x, x0=0, n=1) == (c**2)**0.5 + O(x)
    assert f.taylor_term(0, x) == (c**2)**0.5
    assert f.taylor_term(1, x) == 0.5*x*(c**2)**(-0.5)
    assert f.taylor_term(2, x) == -0.125*x**2*(c**2)**(-1.5)


def test_sympyissue_7638():
    f = pi/log(sqrt(2))
    assert ((1 + I)**(I*f/2))**0.3 == (1 + I)**(0.15*I*f)
    # if 1/3 -> 1.0/3 this should fail since it cannot be shown that the
    # sign will be +/-1; for the previous "small arg" case, it didn't matter
    # that this could not be proved
    assert (1 + I)**(4*I*f) == cbrt((1 + I)**(12*I*f))

    assert cbrt((1 + I)**(I*(1 + 7*f))).exp == Rational(1, 3)
    r = symbols('r', extended_real=True)
    assert sqrt(r**2) == abs(r)
    assert cbrt(r**3) != r
    assert sqrt(Pow(2*I, Rational(5, 2))) != (2*I)**Rational(5, 4)
    p = symbols('p', positive=True)
    assert cbrt(p**2) == p**Rational(2, 3)
    assert NS(((0.2 + 0.7*I)**(0.7 + 1.0*I))**(0.5 - 0.1*I), 1) == '0.4 + 0.2*I'
    assert sqrt(1/(1 + I)) == sqrt((1 - I)/2)  # or 1/sqrt(1 + I)
    e = 1/(1 - sqrt(2))
    assert sqrt(e) == I/sqrt(-1 + sqrt(2))
    assert e**Rational(-1, 2) == -I*sqrt(-1 + sqrt(2))
    assert sqrt((cos(1)**2 + sin(1)**2 - 1)**(3 + I)).exp == Rational(1, 2)
    assert sqrt(r**Rational(4, 3)) != r**Rational(2, 3)
    assert sqrt((p + I)**Rational(4, 3)) == (p + I)**Rational(2, 3)
    assert sqrt((p - p**2*I)**2) == p - p**2*I
    assert sqrt((p + r*I)**2) != p + r*I


def test_sympyissue_8582():
    assert 1**oo is nan
    assert 1**(-oo) is nan
    assert 1**zoo is nan
    assert 1**(oo + I) is nan
    assert 1**(1 + I*oo) is nan
    assert 1**(oo + I*oo) is nan


def test_sympyissue_8650():
    n = Symbol('n', integer=True, nonnegative=True)
    assert (n**n).is_positive is True
    x = 5*n+5
    assert (x**(5*(n+1))).is_positive is True


@pytest.mark.slow
def test_sympyissue_12578():
    s = root(1 - ((x - 1/x)/2)**(-4), 8)
    assert s.series(x, n=17) == (1 - 2*x**4 - 8*x**6 - 34*x**8 -
                                 152*x**10 - 714*x**12 - 3472*x**14 -
                                 17318*x**16 + O(x**17))
    d10 = s.diff((x, 10))
    assert d10.limit(x, 0) == -551577600


def test_sympyissue_13914():
    assert x**zoo is nan


def test_sympyissue_18470():
    assert nan**0 is nan


def test_sympyissue_18499():
    assert (1/oo)**(-oo) is zoo
