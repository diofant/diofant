import pytest

from diofant import (EulerGamma, GoldenRatio, I, Integer, Product, Rational,
                     Sum, Symbol, bell, bernoulli, binomial, catalan, cos, cot,
                     diff, digamma, euler, expand_func, factorial, fibonacci,
                     gamma, genocchi, harmonic, hyper, im, limit, log, lucas,
                     nan, oo, pi, polygamma, re, simplify, sin, sqrt, sstr,
                     symbols, trigamma, zeta, zoo)
from diofant.abc import k, m, n, x
from diofant.functions.combinatorial.numbers import stirling


__all__ = ()


def test_bernoulli():
    assert bernoulli(0) == 1
    assert bernoulli(1) == Rational(-1, 2)
    assert bernoulli(2) == Rational(1, 6)
    assert bernoulli(3) == 0
    assert bernoulli(4) == Rational(-1, 30)
    assert bernoulli(5) == 0
    assert bernoulli(6) == Rational(1, 42)
    assert bernoulli(7) == 0
    assert bernoulli(8) == Rational(-1, 30)
    assert bernoulli(10) == Rational(5, 66)
    assert bernoulli(1000001) == 0

    assert bernoulli(0, x) == 1
    assert bernoulli(1, x) == x - Rational(1, 2)
    assert bernoulli(2, x) == x**2 - x + Rational(1, 6)
    assert bernoulli(3, x) == x**3 - (3*x**2)/2 + x/2

    # Should be fast; computed with mpmath
    b = bernoulli(1000)
    assert b.numerator % 10**10 == 7950421099
    assert b.denominator == 342999030

    b = bernoulli(10**6, evaluate=False).evalf()
    assert str(b) == '-2.23799235765713e+4767529'

    # Issue sympy/sympy#8527
    l = Symbol('l', integer=True)
    m = Symbol('m', integer=True, nonnegative=True)
    n = Symbol('n', integer=True, positive=True)
    assert isinstance(bernoulli(2 * l + 1), bernoulli)
    assert isinstance(bernoulli(2 * m + 1), bernoulli)
    assert bernoulli(2 * n + 1) == 0

    pytest.raises(ValueError, lambda: bernoulli(-20))
    assert isinstance(bernoulli(l, x), bernoulli)


def test_fibonacci():
    assert [fibonacci(n) for n in range(-3, 5)] == [2, -1, 1, 0, 1, 1, 2, 3]
    assert fibonacci(100) == 354224848179261915075
    assert [lucas(n) for n in range(-3, 5)] == [-4, 3, -1, 2, 1, 3, 4, 7]
    assert lucas(100) == 792070839848372253127
    assert lucas(x) == lucas(x, evaluate=False)

    assert fibonacci(1, x) == 1
    assert fibonacci(2, x) == x
    assert fibonacci(3, x) == x**2 + 1
    assert fibonacci(4, x) == x**3 + 2*x

    assert fibonacci(x).rewrite(sqrt) == (GoldenRatio**x - cos(pi*x)/GoldenRatio**x)/sqrt(5)
    assert fibonacci(x).rewrite('tractable') == fibonacci(x).rewrite(sqrt)

    pytest.raises(ValueError, lambda: fibonacci(-2, x))

    n = Symbol('n', integer=True)
    assert isinstance(fibonacci(n, x).rewrite(sqrt), fibonacci)


def test_bell():
    assert [bell(n) for n in range(8)] == [1, 1, 2, 5, 15, 52, 203, 877]

    assert bell(0, x) == 1
    assert bell(1, x) == x
    assert bell(2, x) == x**2 + x
    assert bell(5, x) == x**5 + 10*x**4 + 25*x**3 + 15*x**2 + x

    X = symbols('x:6')
    # X = (x0, x1, .. x5)
    # at the same time: X[1] = x1, X[2] = x2 for standard readablity.
    # but we must supply zero-based indexed object X[1:] = (x1, .. x5)

    assert bell(6, 2, X[1:]) == 6*X[5]*X[1] + 15*X[4]*X[2] + 10*X[3]**2
    assert bell(
        6, 3, X[1:]) == 15*X[4]*X[1]**2 + 60*X[3]*X[2]*X[1] + 15*X[2]**3

    X = (1, 10, 100, 1000, 10000)
    assert bell(6, 2, X) == (6 + 15 + 10)*10000

    X = (1, 2, 3, 3, 5)
    assert bell(6, 2, X) == 6*5 + 15*3*2 + 10*3**2

    X = (1, 2, 3, 5)
    assert bell(6, 3, X) == 15*5 + 60*3*2 + 15*2**3

    # Dobinski's formula
    n = Symbol('n', integer=True, nonnegative=True)
    # For large numbers, this is too slow
    # For nonintegers, there are significant precision errors
    for i in [0, 2, 3, 7, 13, 42, 55]:
        assert bell(i).evalf() == bell(n).rewrite(Sum).evalf(subs={n: i}, strict=False)

    # For negative numbers, the formula does not hold
    m = Symbol('m', integer=True)
    assert bell(-1).evalf() == bell(m).rewrite(Sum).evalf(subs={m: -1})

    assert isinstance(bell(m, x).rewrite(Sum), bell)


def test_harmonic():
    assert harmonic(n, 0) == n
    assert harmonic(n).evalf() == harmonic(n)
    assert harmonic(n, 1) == harmonic(n)
    assert harmonic(1, n).evalf() == harmonic(1, n)

    assert harmonic(0, 1) == 0
    assert harmonic(1, 1) == 1
    assert harmonic(2, 1) == Rational(3, 2)
    assert harmonic(3, 1) == Rational(11, 6)
    assert harmonic(4, 1) == Rational(25, 12)
    assert harmonic(0, 2) == 0
    assert harmonic(1, 2) == 1
    assert harmonic(2, 2) == Rational(5, 4)
    assert harmonic(3, 2) == Rational(49, 36)
    assert harmonic(4, 2) == Rational(205, 144)
    assert harmonic(0, 3) == 0
    assert harmonic(1, 3) == 1
    assert harmonic(2, 3) == Rational(9, 8)
    assert harmonic(3, 3) == Rational(251, 216)
    assert harmonic(4, 3) == Rational(2035, 1728)

    assert harmonic(oo, -1) == nan
    assert harmonic(oo, 0) == oo
    assert harmonic(oo, Rational(1, 2)) == oo
    assert harmonic(oo, 1) == oo
    assert harmonic(oo, 2) == (pi**2)/6
    assert harmonic(oo, 3) == zeta(3)
    assert harmonic(oo, x) == harmonic(oo, x, evaluate=False)


def test_harmonic_rational():
    ne = Integer(6)
    no = Integer(5)
    pe = Integer(8)
    po = Integer(9)
    qe = Integer(10)
    qo = Integer(13)

    Heee = harmonic(ne + pe/qe)
    Aeee = (-log(10) + 2*(Rational(-1, 4) + sqrt(5)/4)*log(sqrt(-sqrt(5)/8 + Rational(5, 8)))
            + 2*(-sqrt(5)/4 - Rational(1, 4))*log(sqrt(sqrt(5)/8 + Rational(5, 8)))
            + pi*(Rational(1, 4) + sqrt(5)/4)/(2*sqrt(-sqrt(5)/8 + Rational(5, 8)))
            + Rational(13944145, 4720968))

    Heeo = harmonic(ne + pe/qo)
    Aeeo = (-log(26) + 2*log(sin(3*pi/13))*cos(4*pi/13) + 2*log(sin(2*pi/13))*cos(32*pi/13)
            + 2*log(sin(5*pi/13))*cos(80*pi/13) - 2*log(sin(6*pi/13))*cos(5*pi/13)
            - 2*log(sin(4*pi/13))*cos(pi/13) + pi*cot(5*pi/13)/2 - 2*log(sin(pi/13))*cos(3*pi/13)
            + Rational(2422020029, 702257080))

    Heoe = harmonic(ne + po/qe)
    Aeoe = (-log(20) + 2*(Rational(1, 4) + sqrt(5)/4)*log(Rational(-1, 4) + sqrt(5)/4)
            + 2*(Rational(-1, 4) + sqrt(5)/4)*log(sqrt(-sqrt(5)/8 + Rational(5, 8)))
            + 2*(-sqrt(5)/4 - Rational(1, 4))*log(sqrt(sqrt(5)/8 + Rational(5, 8)))
            + 2*(-sqrt(5)/4 + Rational(1, 4))*log(Rational(1, 4) + sqrt(5)/4)
            + Rational(11818877030, 4286604231) + pi*(sqrt(5)/8 + Rational(5, 8))/sqrt(-sqrt(5)/8 + Rational(5, 8)))

    Heoo = harmonic(ne + po/qo)
    Aeoo = (-log(26) + 2*log(sin(3*pi/13))*cos(54*pi/13) + 2*log(sin(4*pi/13))*cos(6*pi/13)
            + 2*log(sin(6*pi/13))*cos(108*pi/13) - 2*log(sin(5*pi/13))*cos(pi/13)
            - 2*log(sin(pi/13))*cos(5*pi/13) + pi*cot(4*pi/13)/2
            - 2*log(sin(2*pi/13))*cos(3*pi/13) + Rational(11669332571, 3628714320))

    Hoee = harmonic(no + pe/qe)
    Aoee = (-log(10) + 2*(Rational(-1, 4) + sqrt(5)/4)*log(sqrt(-sqrt(5)/8 + Rational(5, 8)))
            + 2*(-sqrt(5)/4 - Rational(1, 4))*log(sqrt(sqrt(5)/8 + Rational(5, 8)))
            + pi*(Rational(1, 4) + sqrt(5)/4)/(2*sqrt(-sqrt(5)/8 + Rational(5, 8)))
            + Rational(779405, 277704))

    Hoeo = harmonic(no + pe/qo)
    Aoeo = (-log(26) + 2*log(sin(3*pi/13))*cos(4*pi/13) + 2*log(sin(2*pi/13))*cos(32*pi/13)
            + 2*log(sin(5*pi/13))*cos(80*pi/13) - 2*log(sin(6*pi/13))*cos(5*pi/13)
            - 2*log(sin(4*pi/13))*cos(pi/13) + pi*cot(5*pi/13)/2
            - 2*log(sin(pi/13))*cos(3*pi/13) + Rational(53857323, 16331560))

    Hooe = harmonic(no + po/qe)
    Aooe = (-log(20) + 2*(Rational(1, 4) + sqrt(5)/4)*log(Rational(-1, 4) + sqrt(5)/4)
            + 2*(Rational(-1, 4) + sqrt(5)/4)*log(sqrt(-sqrt(5)/8 + Rational(5, 8)))
            + 2*(-sqrt(5)/4 - Rational(1, 4))*log(sqrt(sqrt(5)/8 + Rational(5, 8)))
            + 2*(-sqrt(5)/4 + Rational(1, 4))*log(Rational(1, 4) + sqrt(5)/4)
            + Rational(486853480, 186374097) + pi*(sqrt(5)/8 + Rational(5, 8))/sqrt(-sqrt(5)/8 + Rational(5, 8)))

    Hooo = harmonic(no + po/qo)
    Aooo = (-log(26) + 2*log(sin(3*pi/13))*cos(54*pi/13) + 2*log(sin(4*pi/13))*cos(6*pi/13)
            + 2*log(sin(6*pi/13))*cos(108*pi/13) - 2*log(sin(5*pi/13))*cos(pi/13)
            - 2*log(sin(pi/13))*cos(5*pi/13) + pi*cot(4*pi/13)/2
            - 2*log(sin(2*pi/13))*cos(3*pi/13) + Rational(383693479, 125128080))

    H = [Heee, Heeo, Heoe, Heoo, Hoee, Hoeo, Hooe, Hooo]
    A = [Aeee, Aeeo, Aeoe, Aeoo, Aoee, Aoeo, Aooe, Aooo]

    for h, a in zip(H, A):
        e = expand_func(h).doit()
        assert simplify(e/a) == 1
        assert h.evalf() == a.evalf()


def test_harmonic_evalf():
    assert str(harmonic(1.5).evalf(10)) == '1.280372306'
    assert str(harmonic(1.5, 2).evalf(10)) == '1.154576311'  # issue sympy/sympy#7443


def test_harmonic_rewrite_polygamma():
    assert harmonic(n).rewrite(digamma) == polygamma(0, n + 1) + EulerGamma
    assert harmonic(n).rewrite(trigamma) == polygamma(0, n + 1) + EulerGamma
    assert harmonic(n).rewrite(polygamma) == polygamma(0, n + 1) + EulerGamma

    assert harmonic(n, 3).rewrite(polygamma) == polygamma(2, n + 1)/2 - polygamma(2, 1)/2
    assert harmonic(n, m).rewrite(polygamma) == (-1)**m*(polygamma(m - 1, 1) - polygamma(m - 1, n + 1))/factorial(m - 1)

    assert expand_func(harmonic(n+4)) == harmonic(n) + 1/(n + 4) + 1/(n + 3) + 1/(n + 2) + 1/(n + 1)
    assert expand_func(harmonic(n-4)) == harmonic(n) - 1/(n - 1) - 1/(n - 2) - 1/(n - 3) - 1/n

    assert harmonic(n, m).rewrite('tractable') == harmonic(n, m).rewrite(polygamma).rewrite(gamma).rewrite('tractable')

    assert isinstance(expand_func(harmonic(n, 2)), harmonic)

    assert expand_func(harmonic(n + Rational(1, 2))) == expand_func(harmonic(n + Rational(1, 2)))
    assert expand_func(harmonic(Rational(-1, 2))) == harmonic(Rational(-1, 2))
    assert expand_func(harmonic(x)) == harmonic(x)


def test_harmonic_limit():
    m = Symbol('m', positive=True)
    assert limit(harmonic(n, m + 1), n, oo) == zeta(m + 1)

    assert limit(harmonic(n, 2), n, oo) == pi**2/6
    assert limit(harmonic(n, 3), n, oo) == -polygamma(2, 1)/2


def test_harmonic_rewrite_sum():
    assert harmonic(n).rewrite(Sum) == Sum(1/k, (k, 1, n))
    assert harmonic(n, m).rewrite(Sum) == Sum(k**(-m), (k, 1, n))


def test_euler():
    assert euler(0) == 1
    assert euler(1) == 0
    assert euler(2) == -1
    assert euler(3) == 0
    assert euler(4) == 5
    assert euler(6) == -61
    assert euler(8) == 1385

    assert euler(20, evaluate=False) != 370371188237525

    n = Symbol('n', integer=True)
    assert euler(n) != -1
    assert euler(n).subs({n: 2}) == -1

    assert euler(20).evalf() == 370371188237525.0
    assert euler(20, evaluate=False).evalf() == 370371188237525.0
    assert euler(Rational(1, 2)).evalf() == euler(Rational(1, 2))

    assert euler(n).rewrite(Sum) == euler(n)
    # XXX: Not sure what the guy who wrote this test was trying to do with the _j and _k stuff
    assert euler(2*n + 1).rewrite(Sum) == 0


def test_catalan():
    n = Symbol('n', integer=True)
    m = Symbol('n', integer=True, positive=True)

    catalans = [1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796, 58786]
    for i, c in enumerate(catalans):
        assert catalan(i) == c
        assert catalan(n).rewrite(factorial).subs({n: i}) == c
        assert catalan(n).rewrite(Product).subs({n: i}).doit() == c

    assert catalan(x) == catalan(x)
    assert catalan(2*x).rewrite(binomial) == binomial(4*x, 2*x)/(2*x + 1)
    assert catalan(Rational(1, 2)).rewrite(gamma) == 8/(3*pi)
    assert catalan(Rational(1, 2)).rewrite(factorial).rewrite(gamma) ==\
        8 / (3 * pi)
    assert catalan(3*x).rewrite(gamma) == 4**(
        3*x)*gamma(3*x + Rational(1, 2))/(sqrt(pi)*gamma(3*x + 2))
    assert catalan(x).rewrite(hyper) == hyper((-x + 1, -x), (2,), 1)

    assert catalan(n).rewrite(factorial) == factorial(2*n) / (factorial(n + 1)
                                                              * factorial(n))
    assert isinstance(catalan(n).rewrite(Product), catalan)
    assert isinstance(catalan(m).rewrite(Product), Product)

    assert diff(catalan(x), x) == (polygamma(
        0, x + Rational(1, 2)) - polygamma(0, x + 2) + log(4))*catalan(x)

    assert catalan(x).evalf() == catalan(x)
    c = catalan(Rational(1, 2)).evalf()
    assert str(c) == '0.848826363156775'
    c = catalan(I).evalf(3)
    assert sstr((re(c), im(c))) == '(0.398, -0.0209)'

    # issue sympy/sympy#8601
    n = Symbol('n', integer=True, negative=True)

    assert catalan(n - 1) == 0
    assert catalan(Rational(-1, 2)) == zoo
    assert catalan(-1) == Rational(-1, 2)
    c1 = catalan(-5.6).evalf(strict=False)
    assert str(c1) == '6.93334070531408e-5'
    c2 = catalan(-35.4).evalf(strict=False)
    assert str(c2) == '-4.14189164517449e-24'


def test_genocchi():
    genocchis = [1, -1, 0, 1, 0, -3, 0, 17]
    for n, g in enumerate(genocchis):
        assert genocchi(n + 1) == g

    assert genocchi(Symbol('z', zero=True) + 1) == 1
    pytest.raises(ValueError, lambda: genocchi(Rational(1, 2)))

    m = Symbol('m', integer=True)
    n = Symbol('n', integer=True, positive=True)

    assert genocchi(2*n + 1) == 0
    assert genocchi(m) == genocchi(m)
    assert genocchi(n).rewrite(bernoulli) == 2 * (1 - 2 ** n) * bernoulli(n)
    assert genocchi(x).rewrite(bernoulli) == genocchi(x)

    assert genocchi(2*n).is_odd
    assert genocchi(2*n + 1, evaluate=False).is_odd is None
    assert genocchi(4*n).is_positive
    # This should work for 4 * n - 2, but fails due to some variation of issue
    # 8632 ((4*n-2).is_positive returns None)
    assert genocchi(4*n + 2).is_negative
    assert genocchi(2*n + 1, evaluate=False).is_negative is None

    assert genocchi(m).is_odd is None
    assert genocchi(m).is_negative is None


def test_stirling():
    pytest.raises(ValueError, lambda: stirling(-1, -1))
    pytest.raises(ValueError, lambda: stirling(10, 1, kind=3))
    for kind in [1, 2]:
        assert stirling(0, 0, kind=kind) == 1
    assert [stirling(6, i, kind=1) for i in range(7)] == [0, 120, 274, 225,
                                                          85, 15, 1]
    assert [stirling(4, i, kind=1) for i in range(5)] == [0, 6, 11, 6, 1]
    assert [stirling(4, i, signed=True) for i in range(5)] == [0, -6, 11, -6, 1]
    assert [stirling(10, i) for i in range(12)] == [0, 1, 511, 9330, 34105,
                                                    42525, 22827, 5880, 750,
                                                    45, 1, 0]
    assert stirling(5, 3, 2) == 7


def test_sympyissue_8496():
    pytest.raises(TypeError, lambda: catalan(n, k))
    pytest.raises(TypeError, lambda: euler(n, k))
