"""
Tests from Michael Wester's 1999 paper.

References
==========

* http://www.math.unm.edu/~wester/cas/book/Wester.pdf
* http://math.unm.edu/~wester/cas_review.html
"""
import itertools
import math

import pytest

from diofant import (ZZ, And, BlockMatrix, Complement, Derivative, DiracDelta,
                     E, EulerGamma, FiniteSet, Float, Function, GoldenRatio,
                     GramSchmidt, Heaviside, I, Integral, Lambda, LambertW,
                     LaplaceTransform, Le, Lt, Matrix, MatrixSymbol, Max, Mul,
                     N, O, Or, Piecewise, Product, Rational, Subs, Sum, Symbol,
                     ZeroMatrix, acos, acot, apart, asin, asinh,
                     assoc_legendre, atan, atanh, bernoulli, besselj, binomial,
                     block_collapse, cancel, cbrt, ceiling, chebyshevt,
                     combsimp, conjugate)
from diofant import continued_fraction_convergents as cf_c
from diofant import continued_fraction_iterator as cf_i
from diofant import continued_fraction_periodic as cf_p
from diofant import continued_fraction_reduce as cf_r
from diofant import (cos, cosh, cot, diff, dsolve, elliptic_e, elliptic_f, erf,
                     exp, expand, expand_func, eye, factor, factorial,
                     factorial2, factorint, fibonacci, floor,
                     fourier_transform, gamma, gcd, groebner, hessian, hyper,
                     hyperexpand, im, integrate, inverse_laplace_transform,
                     laplace_transform, legendre_poly, limit, log, logcombine,
                     maximize, mellin_transform, minimize, nan, npartitions,
                     oo, pi, polygamma, polylog, powdenest, powsimp,
                     primerange, primitive, primitive_root, product, radsimp,
                     re, reduce_inequalities, residue, resultant, rf, root,
                     rsolve, sec, sign, simplify, sin, sinh, solve, sqrt,
                     sqrtdenest, summation, symbols, tan, tanh, totient,
                     trigsimp, trunc, wronskian, zeta, zoo)
from diofant.abc import a, b, c, d, s, t, v, w, x, y, z
from diofant.functions.combinatorial.numbers import stirling
from diofant.utilities.iterables import partitions


__all__ = ()

i, j, k, l, m, n = symbols('i j k l m n', integer=True)
f, g, h = symbols('f g h', cls=Function)


# B. Set Theory

def test_B1():
    assert (FiniteSet(i, j, j, k, k, k) | FiniteSet(l, k, j) |
            FiniteSet(j, m, j)) == FiniteSet(i, j, k, l, m)


def test_B2():
    a, b, c = FiniteSet(j), FiniteSet(m), FiniteSet(j, k)
    d, e = FiniteSet(i), FiniteSet(j, k, l)

    assert (FiniteSet(i, j, j, k, k, k) & FiniteSet(l, k, j) &
            FiniteSet(j, m, j)) == a | (b & (c | (d & e)))


def test_B3():
    assert (FiniteSet(i, j, k, l, m) - FiniteSet(j) ==
            Complement(FiniteSet(i, k, l, m), FiniteSet(j)))


def test_B4():
    assert (FiniteSet(*(FiniteSet(i, j)*FiniteSet(k, l))) ==
            FiniteSet((i, k), (i, l), (j, k), (j, l)))


# C. Numbers

def test_C1():
    assert (factorial(50) ==
            30414093201713378043612608166064768844377641568960512000000000000)


def test_C2():
    assert (factorint(factorial(50)) == {2: 47, 3: 22, 5: 12, 7: 8,
                                         11: 4, 13: 3, 17: 2, 19: 2, 23: 2, 29: 1, 31: 1, 37: 1,
                                         41: 1, 43: 1, 47: 1})


def test_C3():
    assert (factorial2(10), factorial2(9)) == (3840, 945)


def test_C7():
    assert log(32768, 8) == 5


def test_C8():
    # Modular multiplicative inverse. Would be nice if divmod could do this.
    assert ZZ.invert(5, 7) == 3
    assert ZZ.invert(5, 6) == 5


def test_C9():
    assert math.gcd(1776, 1554, 5698) == 74


def test_C10():
    x = 0
    for n in range(2, 11):
        x += Rational(1, n)
    assert x == Rational(4861, 2520)


def test_C12():
    assert Rational(7, 11) * Rational(22, 7) == 2


def test_C13():
    test = Rational(10, 7) * root(1 + Rational(29, 1000), 3)
    good = root(3, 3)
    assert test == good


def test_C14():
    assert sqrtdenest(sqrt(2*sqrt(3) + 4)) == 1 + sqrt(3)


def test_C15():
    test = sqrtdenest(sqrt(14 + 3*sqrt(3 + 2*sqrt(5 - 12*sqrt(3 - 2*sqrt(2))))))
    good = sqrt(2) + 3
    assert test == good


def test_C16():
    test = sqrtdenest(sqrt(10 + 2*sqrt(6) + 2*sqrt(10) + 2*sqrt(15)))
    good = sqrt(2) + sqrt(3) + sqrt(5)
    assert test == good


def test_C17():
    test = radsimp((sqrt(3) + sqrt(2)) / (sqrt(3) - sqrt(2)))
    good = 5 + 2*sqrt(6)
    assert test == good


def test_C18():
    assert simplify((sqrt(-2 + sqrt(-5)) * sqrt(-2 - sqrt(-5))).expand(complex=True)) == 3


@pytest.mark.xfail
def test_C19():
    assert simplify(root(90 + 34*sqrt(7), 3)) == 3 + sqrt(7)


@pytest.mark.xfail
def test_C20():
    inside = 135 + 78*sqrt(3)
    test = (root(inside, 3)**2 + 3)*sqrt(3)/root(inside, 3)
    assert simplify(test) == 12


@pytest.mark.xfail
def test_C21():
    assert simplify(root(41 + 29*sqrt(2), 5)) == 1 + sqrt(2)


@pytest.mark.xfail
def test_C22():
    test = ((6 - 4*sqrt(2))*log(3 - 2*sqrt(2)) +
            (3 - 2*sqrt(2))*log(17 - 12*sqrt(2)) +
            32 - 24*sqrt(2))/(48*sqrt(2) - 72)
    assert simplify(test - sqrt(2)/3 + log(sqrt(2) - 1)/3) == 0


def test_C23():
    assert 2 * oo - 3 == oo


# D. Numerical Analysis

def test_D1():
    assert 0.0 / sqrt(2) == 0.0


def test_D2():
    assert str(exp(-1000000).evalf()) == '3.29683147808856e-434295'


def test_D3():
    assert exp(pi*sqrt(163)).evalf(50).num.ae(262537412640768744)


def test_D4():
    assert floor(Rational(-5, 3)) == -2
    assert ceiling(Rational(-5, 3)) == -1


# F. Combinatorial Theory.

def test_F1():
    assert rf(x, 3) == x*(1 + x)*(2 + x)


def test_F2():
    assert expand_func(binomial(n, 3)) == n*(n - 1)*(n - 2)/6


@pytest.mark.xfail
def test_F3():
    assert combsimp(2**n*factorial(n)*factorial2(2*n - 1)) == factorial(2*n)


@pytest.mark.xfail
def test_F4():
    assert combsimp(2**n*factorial(n)*product(2*k - 1, (k, 1, n))) == factorial(2*n)


def test_F5():
    assert combsimp(gamma((2*n + 1)/2)/sqrt(pi)/factorial(n) -
                    factorial(2*n)/2**(2*n)/factorial(n)**2) == 0


def test_F6():
    partTest = [p.copy() for p in partitions(4)]
    partDesired = [{4: 1}, {1: 1, 3: 1}, {2: 2}, {1: 2, 2: 1}, {1: 4}]
    assert partTest == partDesired


def test_F7():
    assert npartitions(4) == 5


def test_F8():
    assert stirling(5, 2, signed=True) == -50  # if signed, then kind=1


def test_F9():
    assert totient(1776) == 576

# G. Number Theory


def test_G1():
    assert list(primerange(999983, 1000004)) == [999983, 1000003]


def test_G2():
    assert primitive_root(191) == 19


# ... G14 Modular equations are not implemented.


def test_G15():
    assert Rational(sqrt(3).evalf()).limit_denominator(15) == Rational(26, 15)
    assert list(itertools.takewhile(lambda x: x.denominator <= 15,
                                    cf_c(cf_i(sqrt(3)))))[-1] == Rational(26, 15)


def test_G16():
    assert list(itertools.islice(cf_i(pi), 10)) == [3, 7, 15, 1, 292,
                                                    1, 1, 1, 2, 1]


def test_G17():
    assert cf_p(0, 1, 23) == [4, [1, 3, 1, 8]]


def test_G18():
    assert cf_p(1, 2, 5) == [[1]]
    assert cf_r([[1]]) == Rational(1, 2) + sqrt(5)/2


def test_G20():
    s = symbols('s', integer=True, positive=True)
    # Wester erroneously has this as -s + sqrt(s**2 + 1)
    assert cf_r([[2*s]]) == s + sqrt(s**2 + 1)


# H. Algebra

def test_H1():
    assert simplify(2*2**n) == simplify(2**(n + 1))
    assert powdenest(2*2**n) == simplify(2**(n + 1))


def test_H2():
    assert powsimp(4 * 2**n) == 2**(n + 2)


def test_H3():
    assert (-1)**(n*(n + 1)) == 1


def test_H4():
    assert primitive(6*x - 10) == (2, 3*x - 5)
    assert factor(6*x - 10) == Mul(2, 3*x - 5, evaluate=False)


p1 = 64*x**34 - 21*x**47 - 126*x**8 - 46*x**5 - 16*x**60 - 81
p2 = 72*x**60 - 25*x**25 - 19*x**23 - 22*x**39 - 83*x**52 + 54*x**10 + 81
q = 34*x**19 - 25*x**16 + 70*x**7 + 20*x**3 - 91*x - 86


def test_H5():
    assert gcd(p1, p2, x) == 1


def test_H6():
    assert gcd(expand(p1 * q), expand(p2 * q)) == q


def test_H7():
    p1 = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    p2 = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    assert gcd(p1, p2, x, y, z) == 1


def test_H8():
    p1 = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    p2 = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    q = 11*x**12*y**7*z**13 - 23*x**2*y**8*z**10 + 47*x**17*y**5*z**8
    assert gcd(p1 * q, p2 * q, x, y, z) == q


def test_H9():
    p1 = 2*x**(n + 4) - x**(n + 2)
    p2 = 4*x**(n + 1) + 3*x**n
    assert gcd(p1, p2) == x**n


def test_H10():
    p1 = 3*x**4 + 3*x**3 + x**2 - x - 2
    p2 = x**3 - 3*x**2 + x + 5
    assert resultant(p1, p2, x) == 0
    assert factor(p1) == (x + 1)*(3*x**3 + x - 2)
    assert factor(p2) == (x + 1)*(x**2 - 4*x + 5)


def test_H11():
    assert resultant(p1 * q, p2 * q, x) == 0


def test_H12():
    num = x**2 - 4
    den = x**2 + 4*x + 4
    assert cancel(num/den) == (x - 2)/(x + 2)
    assert simplify(num/den) == (x - 2)/(x + 2)


@pytest.mark.xfail
def test_H13():
    assert simplify((exp(x) - 1) / (exp(x/2) + 1)) == exp(x/2) - 1


def test_H14():
    p = (x + 1) ** 20
    ep = expand(p)
    assert ep == (1 + 20*x + 190*x**2 + 1140*x**3 + 4845*x**4 + 15504*x**5
                  + 38760*x**6 + 77520*x**7 + 125970*x**8 + 167960*x**9 + 184756*x**10
                  + 167960*x**11 + 125970*x**12 + 77520*x**13 + 38760*x**14 + 15504*x**15
                  + 4845*x**16 + 1140*x**17 + 190*x**18 + 20*x**19 + x**20)
    dep = diff(ep, x)
    assert dep == (20 + 380*x + 3420*x**2 + 19380*x**3 + 77520*x**4
                   + 232560*x**5 + 542640*x**6 + 1007760*x**7 + 1511640*x**8 + 1847560*x**9
                   + 1847560*x**10 + 1511640*x**11 + 1007760*x**12 + 542640*x**13
                   + 232560*x**14 + 77520*x**15 + 19380*x**16 + 3420*x**17 + 380*x**18
                   + 20*x**19)
    assert factor(dep) == 20*(1 + x)**19


def test_H15():
    assert simplify((Mul(*[x - r[x] for r in solve(x**3 + x**2 - 7)]))) == x**3 + x**2 - 7


def test_H16():
    assert (factor(x**15 - 1) ==
            (x - 1)*(x**2 + x + 1)*(x**4 + x**3 + x**2 + x + 1) *
            (x**8 - x**7 + x**5 - x**4 + x**3 - x + 1))
    assert (factor(x**100 - 1) ==
            (x - 1)*(x + 1)*(x**2 + 1)*(x**4 - x**3 + x**2 - x + 1) *
            (x**4 + x**3 + x**2 + x + 1)*(x**8 - x**6 + x**4 - x**2 + 1) *
            (x**20 - x**15 + x**10 - x**5 + 1) *
            (x**20 + x**15 + x**10 + x**5 + 1) *
            (x**40 - x**30 + x**20 - x**10 + 1))


def test_H17():
    assert simplify(factor(expand(p1 * p2)) - p1*p2) == 0


def test_H18():
    e = 4*x**4 + 8*x**3 + 77*x**2 + 18*x + 153
    assert (factor(e, extension=I) ==
            Mul(4, x - 3*I/2, x + 3*I/2, x + 1 - 4*I, x + 1 + 4*I,
                evaluate=False))


def test_H19():
    # The idea is to let a**2 == 2, then solve 1/(a-1). Answer is a+1")
    assert (a - 1).as_poly().invert((a**2 - 2).as_poly()) == a + 1


def test_H20():
    f = x**3 + (sqrt(2) - 2)*x**2 - (2*sqrt(2) + 3)*x - 3*sqrt(2)
    g = x**2 - 2
    r = (x**2 - 2*x - 3)/(x - sqrt(2))
    assert cancel(f/g, extension=True) == cancel(f/g, extension=sqrt(2)) == r


def test_H22():
    f = x**4 - 3*x**2 + 1
    assert factor(f) == (x**2 - x - 1)*(x**2 + x - 1)
    assert factor(f, modulus=5) == (x + 3)**2 * (x + 2)**2


def test_H23():
    f = x**11 + x + 1
    s = (x**2 + x + 1)*(x**9 + 65536*x**8 + x**6 + 65536*x**5 +
                        x**3 + 65536*x**2 + 1)
    assert factor(f, modulus=65537) == s
    assert trunc(expand(s), 65537) == f


def test_H24():
    phi = GoldenRatio.expand(func=True)
    assert (factor(x**4 - 3*x**2 + 1, extension=phi) ==
            (x - phi)*(x + 1 - phi)*(x - 1 + phi)*(x + phi))


def test_H25():
    e = (x - 2*y**2 + 3*z**3) ** 20
    assert factor(expand(e)) == e


def test_H26():
    g = expand((sin(x) - 2*cos(y)**2 + 3*tan(z)**3)**20)
    assert factor(g, expand=False) == (-sin(x) + 2*cos(y)**2 - 3*tan(z)**3)**20


def test_H27():
    f = 24*x*y**19*z**8 - 47*x**17*y**5*z**8 + 6*x**15*y**9*z**2 - 3*x**22 + 5
    g = 34*x**5*y**8*z**13 + 20*x**7*y**7*z**7 + 12*x**9*y**16*z**4 + 80*y**14*z
    h = -2*z*y**7 \
        * (6*x**9*y**9*z**3 + 10*x**7*z**6 + 17*y*x**5*z**12 + 40*y**7) \
        * (3*x**22 + 47*x**17*y**5*z**8 - 6*x**15*y**9*z**2 - 24*x*y**19*z**8 - 5)
    assert factor(expand(f*g)) == h


def test_H28():
    f = expand((1 - c**2)**5 * (1 - s**2)**5 * (c**2 + s**2)**10)
    g = groebner([f, c**2 + s**2 - 1])
    assert g[-1] == c**20 - 5*c**18 + 10*c**16 - 10*c**14 + 5*c**12 - c**10
    assert factor(g[-1]) == c**10 * (c - 1)**5 * (c + 1)**5
    assert solve([f, c**2 + s**2 - 1]) == [{c: -1, s: 0}, {c: 0, s: -1},
                                           {c: 0, s: 1}, {c: 1, s: 0}]


@pytest.mark.xfail
def test_H29():
    assert factor(4*x**2 - 21*x*y + 20*y**2, modulus=3) == (x + y)*(x - y)


def test_H30():
    test = factor(x**3 + y**3, extension=sqrt(-3))
    answer = (x + y)*(x + y*(-Rational(1, 2) - sqrt(3)/2*I))*(x + y*(-Rational(1, 2) + sqrt(3)/2*I))
    assert answer == test


def test_H31():
    f = (x**2 + 2*x + 3)/(x**3 + 4*x**2 + 5*x + 2)
    g = 2/(x + 1)**2 - 2/(x + 1) + 3/(x + 2)
    assert apart(f) == g
    assert cancel(g) == f


# I. Trigonometry

@pytest.mark.xfail
def test_I1():
    assert simplify(tan(7*pi/10) + sqrt(1 + 2/sqrt(5))) == 0


@pytest.mark.xfail
def test_I2():
    assert simplify(sqrt((1 + cos(6))/2) + cos(3)) == 0


def test_I3():
    assert cos(n*pi) + sin((4*n - 1)*pi/2) == (-1)**n - 1


def test_I4():
    assert cos(pi*cos(n*pi)) + sin(pi/2*cos(n*pi)) == (-1)**n - 1


@pytest.mark.xfail
def test_I5():
    assert simplify(sin((n**5/5 + n**4/2 + n**3/3 - n/30)*pi)) == 0


def test_I7():
    assert simplify(cos(3*x)/cos(x) - cos(x)**2 + 3*sin(x)**2) == 0


def test_I8():
    assert simplify(cos(3*x)/cos(x) - 2*cos(2*x) + 1) == 0


def test_I9():
    assert simplify(cos(3*x)/cos(x) - cos(x)**2 + 3*sin(x)**2) == 0


def test_I10():
    assert trigsimp((tan(x)**2 + 1 - cos(x)**-2) / (sin(x)**2 + cos(x)**2 - 1)) == nan


# J. Special functions.

def test_J1():
    assert bernoulli(16) == Rational(-3617, 510)


def test_J2():
    assert diff(elliptic_e(x, y**2), y) == (elliptic_e(x, y**2) - elliptic_f(x, y**2))/y


def test_J4():
    assert gamma(Rational(-1, 2)) == -2*sqrt(pi)


def test_J5():
    assert polygamma(0, Rational(1, 3)) == -EulerGamma - pi/2*sqrt(Rational(1, 3)) - Rational(3, 2)*log(3)


def test_J6():
    assert besselj(2, 1 + I).evalf() == Float('0.04157988694396212', 15) + Float('0.24739764151330632', 15)*I


def test_J7():
    assert simplify(besselj(Rational(-5, 2), pi/2)) == 12/(pi**2)


def test_J8():
    p = besselj(Rational(3, 2), z)
    q = (sin(z)/z - cos(z))/sqrt(pi*z/2)
    assert simplify(expand_func(p) - q) == 0


def test_J9():
    assert besselj(0, z).diff(z) == - besselj(1, z)


def test_J10():
    mu, nu = symbols('mu, nu', integer=True)
    assert assoc_legendre(nu, mu, 0) == 2**mu*sqrt(pi)/gamma((nu - mu)/2 + 1)/gamma((-nu - mu + 1)/2)


def test_J11():
    assert simplify(assoc_legendre(3, 1, x)) == simplify(-Rational(3, 2)*sqrt(1 - x**2)*(5*x**2 - 1))


@pytest.mark.slow
def test_J12():
    assert simplify(chebyshevt(1008, x) - 2*x*chebyshevt(1007, x) + chebyshevt(1006, x)) == 0


def test_J13():
    a = symbols('a', integer=True, negative=False)
    assert chebyshevt(a, -1) == (-1)**a


def test_J14():
    p = hyper([Rational(1, 2), Rational(1, 2)], [Rational(3, 2)], z**2)
    assert hyperexpand(p) == asin(z)/z


@pytest.mark.xfail
def test_J17():
    assert integrate(f((x + 2)/5)*DiracDelta((x - 2)/3) -
                     g(x)*diff(DiracDelta(x - 1), x),
                     (x, 0, 3)) == 3*f(Rational(4, 5)) + g(x).diff(x).subs(x, 1)


# K. The Complex Domain

def test_K1():
    z1, z2 = symbols('z1, z2', complex=True)
    assert re(z1 + I*z2) == -im(z2) + re(z1)
    assert im(z1 + I*z2) == im(z1) + re(z2)


def test_K2():
    assert abs(3 - sqrt(7) + I*sqrt(6*sqrt(7) - 15)) == 1


@pytest.mark.xfail
def test_K3():
    a, b = symbols('a, b', real=True)
    assert simplify(abs(1/(a + I/a + I*b)) - 1/sqrt(a**2 + (1/a + b)**2)) == 0


def test_K4():
    assert log(3 + 4*I).expand(complex=True) == log(5) + I*atan(Rational(4, 3))


def test_K5():
    x, y = symbols('x, y', real=True)
    assert tan(x + I*y).expand(complex=True) == (sin(2*x)/(cos(2*x) +
                                                           cosh(2*y)) + I*sinh(2*y)/(cos(2*x) + cosh(2*y)))


def test_K6():
    assert sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z)) == sqrt(x*y)/sqrt(x)
    assert sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z)) != sqrt(y)


def test_K7():
    y = symbols('y', real=True, negative=False)
    expr = sqrt(x*y*abs(z)**2)/(sqrt(x)*abs(z))
    sexpr = simplify(expr)
    assert sexpr == sqrt(y)


def test_K8():
    assert simplify(sqrt(1/z) - 1/sqrt(z)) != 0


@pytest.mark.xfail
def test_K8_fail():
    z = symbols('z', negative=False)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) == 0


def test_K9():
    z = symbols('z', real=True, positive=True)
    assert simplify(sqrt(1/z) - 1/sqrt(z)) == 0


def test_K10():
    z = symbols('z', real=True, negative=True)
    assert simplify(sqrt(1/z) + 1/sqrt(z)) == 0


# L. Determining Zero Equivalence

def test_L1():
    assert sqrt(997) - root(997**3, 6) == 0


def test_L2():
    assert sqrt(999983) - root(999983**3, 6) == 0


def test_L3():
    assert simplify((root(2, 3) + root(4, 3))**3 - 6*(root(2, 3) + root(4, 3)) - 6) == 0


def test_L4():
    assert trigsimp(cos(x)**3 + cos(x)*sin(x)**2 - cos(x)) == 0


@pytest.mark.xfail
def test_L5():
    assert simplify(log(tan(x/2 + pi/4)) - asinh(tan(x))) == 0


def test_L6():
    assert (log(tan(x/2 + pi/4)) - asinh(tan(x))).diff(x).subs({x: 0}) == 0


@pytest.mark.xfail
def test_L7():
    assert simplify(log((2*sqrt(x) + 1)/(sqrt(4*x + 4*sqrt(x) + 1)))) == 0


@pytest.mark.xfail
def test_L8():
    assert simplify((4*x + 4*sqrt(x) + 1)**(sqrt(x)/(2*sqrt(x) + 1))
                    * (2*sqrt(x) + 1)**(1/(2*sqrt(x) + 1)) - 2*sqrt(x) - 1) == 0


@pytest.mark.xfail
def test_L9():
    assert simplify(2**(1 - z)*gamma(z)*zeta(z)*cos(z*pi/2) - pi**z*zeta(1 - z)) == 0


# M. Equations

def test_M2():
    # The roots of this equation should all be real. Note that this
    # doesn't test that they are correct.
    sol = solve(3*x**3 - 18*x**2 + 33*x - 19, x)
    assert all(s[x].expand(complex=True).is_extended_real for s in sol)


def test_M5():
    assert (solve(x**6 - 9*x**4 - 4*x**3 + 27*x**2 - 36*x - 23, x) ==
            [{x: r} for r in (x**6 - 9*x**4 - 4*x**3 + 27*x**2 -
                              36*x - 23).as_poly().all_roots()])


def test_M6():
    assert ({s[x] for s in solve(x**7 - 1, x)} ==
            {cos(n*2*pi/7) + I*sin(n*2*pi/7) for n in range(7)})
    # The paper asks for exp terms, but sin's and cos's may be acceptable;
    # if the results are simplified, exp terms appear for all but
    # -sin(pi/14) - I*cos(pi/14) and -sin(pi/14) + I*cos(pi/14) which
    # will simplify if you apply the transformation foo.rewrite(exp).expand()


def test_M7():
    sol = solve(x**8 - 8*x**7 + 34*x**6 - 92*x**5 + 175*x**4 - 236*x**3 +
                226*x**2 - 140*x + 46, x)
    assert [s[x].simplify() for s in sol] == [
        1 - sqrt(-6 - 2*I*sqrt(3 + 4*sqrt(3)))/2,
        1 + sqrt(-6 - 2*I*sqrt(3 + 4*sqrt(3)))/2,
        1 - sqrt(-6 + 2*I*sqrt(3 + 4*sqrt(3)))/2,
        1 + sqrt(-6 + 2*I*sqrt(3 + 4*sqrt(3)))/2,
        1 - sqrt(-6 + 2*sqrt(-3 + 4*sqrt(3)))/2,
        1 + sqrt(-6 + 2*sqrt(-3 + 4*sqrt(3)))/2,
        1 - sqrt(-6 - 2*sqrt(-3 + 4*sqrt(3)))/2,
        1 + sqrt(-6 - 2*sqrt(-3 + 4*sqrt(3)))/2]


def test_M8():
    # There should be an infinite number of solutions.
    assert solve(exp(2*x) + 2*exp(x) + 1 - z, x) == [{x: log(-sqrt(z) - 1)},
                                                     {x: log(+sqrt(z) - 1)}]


def test_M10():
    assert solve(exp(x) - x, x) == [{x: -LambertW(-1)}]


@pytest.mark.xfail
def test_M11():
    assert solve(x**x - x, x) != [{x: 1}]


def test_M12():
    # TODO: x = [-1, 2*(+/-asinh(1)*I + n*pi}, 3*(pi/6 + n*pi/3)]
    assert solve((x + 1)*(sin(x)**2 + 1)**2*cos(3*x)**3, x) == [
        {x: -1}, {x: pi/6}, {x: pi/2},
        {x: -I*log(1 + sqrt(2))}, {x: I*log(1 + sqrt(2))},
        {x: pi - I*log(1 + sqrt(2))}, {x: pi + I*log(1 + sqrt(2))},
    ]


def test_M13():
    assert solve(sin(x) - cos(x), x) == [{x: -3*pi/4}, {x: pi/4}]


def test_M14():
    assert solve(tan(x) - 1, x) == [{x: pi/4}]


def test_M15():
    assert solve(sin(x) - Rational(1, 2)) == [{x: pi/6}, {x: 5*pi/6}]


def test_M16():
    assert solve(sin(x) - tan(x), x) == [{x: 0}, {x: -pi},
                                         {x: pi}, {x: 2*pi}]


@pytest.mark.xfail
def test_M17():
    assert solve(asin(x) - atan(x), x) == [{x: 0}]


@pytest.mark.xfail
def test_M18():
    assert solve(acos(x) - atan(x), x) != []  # == [{x: sqrt((sqrt(5) - 1)/2)}]


def test_M19():
    assert solve((x - 2)/root(x, 3), x) == [{x: 2}]


def test_M20():
    assert solve(sqrt(x**2 + 1) - x + 2, x) == []


def test_M21():
    assert solve(x + sqrt(x) - 2) == [{x: 1}]


def test_M22():
    assert solve(2*sqrt(x) + 3*root(x, 4) - 2) == [{x: Rational(1, 16)}]


def test_M23():
    assert solve(x - 1/sqrt(1 + x**2)) == [{x: -I*sqrt(Rational(1, 2) + sqrt(5)/2)},
                                           {x: sqrt(Rational(-1, 2) + sqrt(5)/2)}]


def test_M24():
    solution = solve(1 - binomial(m, 2)*2**k, k)
    answer = log(2/(m*(m - 1)), 2)
    assert solution[0][k].expand() == answer.expand()


def test_M25():
    a, b, c, d = symbols(':d', positive=True)
    assert solve(a*b**x - c*d**x, x)[0][x].expand() == (log(c/a)/log(b/d)).expand()


def test_M26():
    assert solve(sqrt(log(x)) - log(sqrt(x))) == [{x: 1}, {x: exp(4)}]


@pytest.mark.xfail
def test_M28():
    assert solve(5*x + exp((x - 5)/2) - 8*x**3, x) != []


def test_M29():
    x = symbols('x', real=True)
    assert solve(abs(x - 1) - 2) == [{x: -1}, {x: 3}]


def test_M30():
    x = symbols('x', real=True)
    assert solve(abs(2*x + 5) - abs(x - 2), x) == [{x: -7}, {x: -1}]


def test_M31():
    x = symbols('x', real=True)
    assert solve(1 - abs(x) - Max(-x - 2, x - 2), x) == [{x: -Rational(3, 2)},
                                                         {x: Rational(3, 2)}]


@pytest.mark.xfail
def test_M32():
    assert solve(Max(2 - x**2, x) - Max(-x, (x**3)/9), x) == [{x: -1}, {x: 3}]


@pytest.mark.xfail
def test_M34():
    assert solve((1 + I)*z + (2 - I)*conjugate(z) + 3*I, z) == [{z: 2 + 3*I}]


def test_M35():
    x, y = symbols('x y', real=True)
    assert solve((3*x - 2*y - I*y + 3*I).as_real_imag()) == [{y: 3, x: 2}]


@pytest.mark.xfail
def test_M36():
    assert solve(f(x)**2 + f(x) - 2, x) != []


def test_M37():
    assert solve([x + y + z - 6, 2*x + y + 2*z - 10, x + 3*y + z - 10]) == [{x: -z + 4, y: 2}]


@pytest.mark.slow
def test_M38():
    [k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15,
     k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29,
     k30, k31, k32, k33, k34, k35, k36, k37, k38, k39, k40, k41, k42, k43,
     k44, k45, k46, k47, k48, k49] = variables = symbols('k1:50')
    del k43
    system = [
        -b*k8/a + c*k8/a, -b*k11/a + c*k11/a, -b*k10/a + c*k10/a + k2, -k3 - b*k9/a + c*k9/a,
        -b*k14/a + c*k14/a, -b*k15/a + c*k15/a, -b*k18/a + c*k18/a - k2, -b*k17/a + c*k17/a,
        -b*k16/a + c*k16/a + k4, -b*k13/a + c*k13/a - b*k21/a + c*k21/a + b*k5/a - c*k5/a,
        b*k44/a - c*k44/a, -b*k45/a + c*k45/a, -b*k20/a + c*k20/a, -b*k44/a + c*k44/a,
        b*k46/a - c*k46/a, b**2*k47/a**2 - 2*b*c*k47/a**2 + c**2*k47/a**2, k3, -k4,
        -b*k12/a + c*k12/a - a*k6/b + c*k6/b, -b*k19/a + c*k19/a + a*k7/c - b*k7/c,
        b*k45/a - c*k45/a, -b*k46/a + c*k46/a, -k48 + c*k48/a + c*k48/b - c**2*k48/(a*b),
        -k49 + b*k49/a + b*k49/c - b**2*k49/(a*c), a*k1/b - c*k1/b, a*k4/b - c*k4/b,
        a*k3/b - c*k3/b + k9, -k10 + a*k2/b - c*k2/b, a*k7/b - c*k7/b, -k9, k11,
        b*k12/a - c*k12/a + a*k6/b - c*k6/b, a*k15/b - c*k15/b, k10 + a*k18/b - c*k18/b,
        -k11 + a*k17/b - c*k17/b, a*k16/b - c*k16/b, -a*k13/b + c*k13/b + a*k21/b - c*k21/b + a*k5/b - c*k5/b,
        -a*k44/b + c*k44/b, a*k45/b - c*k45/b, a*k14/c - b*k14/c + a*k20/b - c*k20/b,
        a*k44/b - c*k44/b, -a*k46/b + c*k46/b, -k47 + c*k47/a + c*k47/b - c**2*k47/(a*b),
        a*k19/b - c*k19/b, -a*k45/b + c*k45/b, a*k46/b - c*k46/b, a**2*k48/b**2 - 2*a*c*k48/b**2 + c**2*k48/b**2,
        -k49 + a*k49/b + a*k49/c - a**2*k49/(b*c), k16, -k17, -a*k1/c + b*k1/c,
        -k16 - a*k4/c + b*k4/c, -a*k3/c + b*k3/c, k18 - a*k2/c + b*k2/c, b*k19/a - c*k19/a - a*k7/c + b*k7/c,
        -a*k6/c + b*k6/c, -a*k8/c + b*k8/c, -a*k11/c + b*k11/c + k17, -a*k10/c + b*k10/c - k18,
        -a*k9/c + b*k9/c, -a*k14/c + b*k14/c - a*k20/b + c*k20/b, -a*k13/c + b*k13/c + a*k21/c - b*k21/c - a*k5/c + b*k5/c,
        a*k44/c - b*k44/c, -a*k45/c + b*k45/c, -a*k44/c + b*k44/c, a*k46/c - b*k46/c,
        -k47 + b*k47/a + b*k47/c - b**2*k47/(a*c), -a*k12/c + b*k12/c, a*k45/c - b*k45/c,
        -a*k46/c + b*k46/c, -k48 + a*k48/b + a*k48/c - a**2*k48/(b*c),
        a**2*k49/c**2 - 2*a*b*k49/c**2 + b**2*k49/c**2, k8, k11, -k15, k10 - k18,
        -k17, k9, -k16, -k29, k14 - k32, -k21 + k23 - k31, -k24 - k30, -k35, k44,
        -k45, k36, k13 - k23 + k39, -k20 + k38, k25 + k37, b*k26/a - c*k26/a - k34 + k42,
        -2*k44, k45, k46, b*k47/a - c*k47/a, k41, k44, -k46, -b*k47/a + c*k47/a,
        k12 + k24, -k19 - k25, -a*k27/b + c*k27/b - k33, k45, -k46, -a*k48/b + c*k48/b,
        a*k28/c - b*k28/c + k40, -k45, k46, a*k48/b - c*k48/b, a*k49/c - b*k49/c,
        -a*k49/c + b*k49/c, -k1, -k4, -k3, k15, k18 - k2, k17, k16, k22, k25 - k7,
        k24 + k30, k21 + k23 - k31, k28, -k44, k45, -k30 - k6, k20 + k32, k27 + b*k33/a - c*k33/a,
        k44, -k46, -b*k47/a + c*k47/a, -k36, k31 - k39 - k5, -k32 - k38, k19 - k37,
        k26 - a*k34/b + c*k34/b - k42, k44, -2*k45, k46, a*k48/b - c*k48/b,
        a*k35/c - b*k35/c - k41, -k44, k46, b*k47/a - c*k47/a, -a*k49/c + b*k49/c,
        -k40, k45, -k46, -a*k48/b + c*k48/b, a*k49/c - b*k49/c, k1, k4, k3, -k8,
        -k11, -k10 + k2, -k9, k37 + k7, -k14 - k38, -k22, -k25 - k37, -k24 + k6,
        -k13 - k23 + k39, -k28 + b*k40/a - c*k40/a, k44, -k45, -k27, -k44, k46,
        b*k47/a - c*k47/a, k29, k32 + k38, k31 - k39 + k5, -k12 + k30, k35 - a*k41/b + c*k41/b,
        -k44, k45, -k26 + k34 + a*k42/c - b*k42/c, k44, k45, -2*k46, -b*k47/a + c*k47/a,
        -a*k48/b + c*k48/b, a*k49/c - b*k49/c, k33, -k45, k46, a*k48/b - c*k48/b,
        -a*k49/c + b*k49/c
    ]
    solution = [{
        k49: 0, k48: 0, k47: 0, k46: 0, k45: 0, k44: 0, k41: 0, k40: 0,
        k38: 0, k37: 0, k36: 0, k35: 0, k33: 0, k32: 0, k30: 0, k29: 0,
        k28: 0, k27: 0, k25: 0, k24: 0, k22: 0, k21: 0, k20: 0, k19: 0,
        k18: 0, k17: 0, k16: 0, k15: 0, k14: 0, k13: 0, k12: 0, k11: 0,
        k10: 0, k9:  0, k8:  0, k7:  0, k6:  0, k5:  0, k4:  0, k3:  0,
        k2:  0, k1:  0,
        k34: b/c*k42, k31: k39, k26: a/c*k42, k23: k39
    }]
    assert solve(system, *variables) == solution


@pytest.mark.slow
def test_M39():
    r0, r1, r2, r3, r4 = (6*z**5 - 6*z**4 - 9*z**3 -
                          7*z**2 - 3*z - 1).as_poly().all_roots()
    sol = [{x: -1, y: 1, z: 1}, {x: 1, y: 1, z: 1}]
    for r in [r1, r2, r3, r4, r0]:
        sol.extend([{x: -sqrt(3)*sqrt(r*(-12*r**3 + 12*r**2 + 30*r - 7))/3,
                     y: -9*r**4 + 12*r**3 + 21*r**2/2 + 6*r + 3/2, z: r},
                    {x: sqrt(3)*sqrt(r*(-12*r**3 + 12*r**2 + 30*r - 7))/3,
                     y: -9*r**4 + 12*r**3 + 21*r**2/2 + 6*r + 3/2, z: r}])
    sol.extend([{x: -sqrt(-1 - sqrt(2)*I), y: sqrt(2)*I,
                 z: Rational(1, 3) - sqrt(2)*I/3},
                {x: sqrt(-1 - sqrt(2)*I), y: sqrt(2)*I,
                 z: Rational(1, 3) - sqrt(2)*I/3},
                {x: -sqrt(-1 + sqrt(2)*I), y: -sqrt(2)*I,
                 z: Rational(1, 3) + sqrt(2)*I/3},
                {x: sqrt(-1 + sqrt(2)*I), y: -sqrt(2)*I,
                 z: Rational(1, 3) + sqrt(2)*I/3}])

    assert solve((x**2*y + 3*y*z - 4, -3*x**2*z + 2*y**2 + 1,
                  2*y*z**2 - z**2 - 1), check=False) == sol


# N. Inequalities

def test_N1():
    assert E**pi > pi**E


def test_N2():
    x = symbols('x', real=True)
    assert reduce_inequalities(x**4 - x + 1 > 0)
    assert reduce_inequalities(x**4 - x + 1 > 1) == Or(Lt(1, x), x < 0)


def test_N9():
    x = Symbol('x', real=True)
    assert reduce_inequalities(abs(x - 1) > 2) == Or(And(Lt(-oo, x), Lt(x, -1)),
                                                     And(Lt(3, x), Lt(x, oo)))


def test_N10():
    x = Symbol('x', real=True)
    p = (x - 1)*(x - 2)*(x - 3)*(x - 4)*(x - 5)
    assert (reduce_inequalities(expand(p) < 0) ==
            Or(And(Lt(-oo, x), Lt(x, 1)), And(Lt(2, x), Lt(x, 3)),
               And(Lt(4, x), Lt(x, 5))))


def test_N11():
    x = Symbol('x', real=True)
    assert reduce_inequalities(6/(x - 3) <= 3) == Or(And(Le(5, x), Lt(x, oo)),
                                                     And(Lt(-oo, x), Lt(x, 3)))


def test_N12():
    x = Symbol('x', real=True)
    assert reduce_inequalities(sqrt(x) < 2) == And(-oo < x, x < 4)


@pytest.mark.xfail
def test_N13():
    x = Symbol('x', real=True)
    reduce_inequalities(sin(x) < 2)


def test_N14():
    x = Symbol('x', real=True)
    assert reduce_inequalities(sin(x) < 1) == Or(pi/2 < x, x < pi/2)


@pytest.mark.xfail
def test_N15():
    r, t = symbols('r t', real=True)
    reduce_inequalities(abs(2*r*(cos(t) - 1) + 1) <= 1, r)


@pytest.mark.xfail
def test_N16():
    r, t = symbols('r t', real=True)
    reduce_inequalities((r**2)*((cos(t) - 4)**2)*sin(t)**2 < 9, r)


def test_N17():
    assert reduce_inequalities([x + y > 0,
                                x - y < 0]) == ((y > 0) & (x < y) & (x > -y))


# O. Vector Analysis

def test_O1():
    M = Matrix((1 + I, -2, 3*I))
    assert sqrt(expand(M.dot(M.H))) == sqrt(15)


def test_O2():
    assert Matrix((2, 2, -3)).cross(Matrix((1, 3, 1))) == Matrix([[11],
                                                                  [-5],
                                                                  [4]])


def test_O10():
    L = [Matrix([2, 3, 5]), Matrix([3, 6, 2]), Matrix([8, 3, 6])]
    assert GramSchmidt(L) == [Matrix([
                              [2],
                              [3],
                              [5]]),
                              Matrix([
                                  [Rational(23, 19)],
                                  [Rational(63, 19)],
                                  [Rational(-47, 19)]]),
                              Matrix([
                                  [Rational(1692, 353)],
                                  [Rational(-1551, 706)],
                                  [Rational(-423, 706)]])]


# P. Matrix Theory

def test_P2():
    M = Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    del M[1, :]
    del M[:, 2]
    assert M == Matrix([[1, 2],
                        [7, 8]])


def test_P3():
    A = Matrix([[11, 12, 13, 14],
                [21, 22, 23, 24],
                [31, 32, 33, 34],
                [41, 42, 43, 44]])
    A11 = A[0:3, 1:4]
    A12 = A[(0, 1, 3), (2, 0, 3)]
    A21 = A
    A221 = A[0:2, 2:4]
    A222 = A[(3, 0), (2, 1)]
    A22 = A221.col_join(A222)
    B = A11.row_join(A12)
    B = B.col_join(A21.row_join(A22))
    assert B == Matrix([[12, 13, 14, 13, 11, 14],
                        [22, 23, 24, 23, 21, 24],
                        [32, 33, 34, 43, 41, 44],
                        [11, 12, 13, 14, 13, 14],
                        [21, 22, 23, 24, 23, 24],
                        [31, 32, 33, 34, 43, 42],
                        [41, 42, 43, 44, 13, 12]])


@pytest.mark.xfail
def test_P5():
    M = Matrix([[7, 11],
                [3, 8]])
    assert M % 2 == Matrix([[1, 1],
                            [1, 0]])


def test_P5_workaround():
    M = Matrix([[7, 11],
                [3, 8]])
    assert M.applyfunc(lambda i: i % 2) == Matrix([[1, 1],
                                                   [1, 0]])


def test_P6():
    M = Matrix([[cos(x), sin(x)],
                [-sin(x), cos(x)]])
    assert M.diff((x, 2)) == Matrix([[-cos(x), -sin(x)],
                                     [sin(x), -cos(x)]])


def test_P7():
    M = Matrix([[x, y]])*(
        z*Matrix([[1, 3, 5],
                  [2, 4, 6]]) + Matrix([[7, -9, 11],
                                        [-8, 10, -12]]))
    assert M == Matrix([[x*(z + 7) + y*(2*z - 8), x*(3*z - 9) + y*(4*z + 10),
                         x*(5*z + 11) + y*(6*z - 12)]])


@pytest.mark.xfail
def test_P8():
    M = Matrix([[1, -2*I],
                [-3*I, 4]])
    assert M.norm(ord=oo) == 7


def test_P9():
    a, b, c = symbols('a b c', real=True)
    M = Matrix([[a/(b*c), 1/c, 1/b],
                [1/c, b/(a*c), 1/a],
                [1/b, 1/a, c/(a*b)]])
    assert factor(M.norm('fro')) == (a**2 + b**2 + c**2)/(abs(a)*abs(b)*abs(c))


def test_P10():
    M = Matrix([[1, 2 + 3*I],
                [f(4 - 5*I), 6]])
    assert M.H == Matrix([[1, f(4 - 5*I).conjugate()],
                          [2 - 3*I, 6]])


def test_P11():
    assert (simplify(Matrix([[x, y], [1, x*y]]).inv()) ==
            Matrix([[x, -1], [-1/y, x/y]])/(x**2 - 1))


def test_P12():
    A11 = MatrixSymbol('A11', n, n)
    A12 = MatrixSymbol('A12', n, n)
    A22 = MatrixSymbol('A22', n, n)
    B = BlockMatrix([[A11, A12],
                     [ZeroMatrix(n, n), A22]])
    assert block_collapse(B.inverse()) == BlockMatrix([[A11.inverse(), (-1)*A11.inverse()*A12*A22.inverse()],
                                                       [ZeroMatrix(n, n), A22.inverse()]])


def test_P13():
    M = Matrix([[1,     x - 2,                         x - 3],
                [x - 1, x**2 - 3*x + 6,       x**2 - 3*x - 2],
                [x - 2, x**2 - 8,       2*(x**2) - 12*x + 14]])
    L, U, _ = M.LUdecomposition()
    assert simplify(L) == Matrix([[1,     0,     0],
                                  [x - 1, 1,     0],
                                  [x - 2, x - 3, 1]])
    assert simplify(U) == Matrix([[1, x - 2, x - 3],
                                  [0,     4, x - 5],
                                  [0,     0, x - 7]])


def test_P14():
    M = Matrix([[1, 2, 3, 1, 3],
                [3, 2, 1, 1, 7],
                [0, 2, 4, 1, 1],
                [1, 1, 1, 1, 4]])
    R, _ = M.rref()
    assert R == Matrix([[1, 0, -1, 0,  2],
                        [0, 1,  2, 0, -1],
                        [0, 0,  0, 1,  3],
                        [0, 0,  0, 0,  0]])


def test_P15():
    M = Matrix([[-1, 3,  7, -5],
                [4, -2,  1,  3],
                [2,  4, 15, -7]])
    assert M.rank() == 2


def test_P16():
    M = Matrix([[2*sqrt(2), 8],
                [6*sqrt(6), 24*sqrt(3)]])
    assert M.rank() == 1


def test_P17():
    t = symbols('t', real=True)
    M = Matrix([[sin(2*t), cos(2*t)],
                [2*(1 - (cos(t)**2))*cos(t), (1 - 2*(sin(t)**2))*sin(t)]])
    assert M.rank() == 1


def test_P18():
    M = Matrix([[1,  0, -2, 0],
                [-2, 1,  0, 3],
                [-1, 2, -6, 6]])
    assert M.nullspace() == [Matrix([[2],
                                     [4],
                                     [1],
                                     [0]]),
                             Matrix([[0],
                                     [-3],
                                     [0],
                                     [1]])]


def test_P19():
    M = Matrix([[1,    1,    1,    1],
                [w,    x,    y,    z],
                [w**2, x**2, y**2, z**2],
                [w**3, x**3, y**3, z**3]])
    assert M.det() == (w**3*x**2*y - w**3*x**2*z - w**3*x*y**2 + w**3*x*z**2
                       + w**3*y**2*z - w**3*y*z**2 - w**2*x**3*y + w**2*x**3*z
                       + w**2*x*y**3 - w**2*x*z**3 - w**2*y**3*z + w**2*y*z**3
                       + w*x**3*y**2 - w*x**3*z**2 - w*x**2*y**3 + w*x**2*z**3
                       + w*y**3*z**2 - w*y**2*z**3 - x**3*y**2*z + x**3*y*z**2
                       + x**2*y**3*z - x**2*y*z**3 - x*y**3*z**2 + x*y**2*z**3
                       )


def test_P21():
    M = Matrix([[5, -3, -7],
                [-2, 1,  2],
                [2, -3, -4]])
    assert M.charpoly(x).as_expr() == x**3 - 2*x**2 - 5*x + 6


@pytest.mark.slow
def test_P22():
    # Wester test calculates eigenvalues for a diagonal matrix of dimension 100
    # This currently takes forever with diofant:
    # M = (2 - x)*eye(100);
    # assert M.eigenvals() == {-x + 2: 100}
    # So we will speed-up the test checking only for dimension=12
    d = 12
    M = (2 - x)*eye(d)
    assert M.eigenvals() == {-x + 2: d}


def test_P23():
    M = Matrix([
        [2, 1, 0, 0, 0],
        [1, 2, 1, 0, 0],
        [0, 1, 2, 1, 0],
        [0, 0, 1, 2, 1],
        [0, 0, 0, 1, 2]])
    assert M.eigenvals() == {1: 1, 2: 1, 3: 1, sqrt(3) + 2: 1, -sqrt(3) + 2: 1}


def test_P24():
    M = Matrix([[611,  196,  -192,  407,   -8,  -52,  -49,   29],
                [196,  899,   113, -192,  -71,  -43,   -8,  -44],
                [-192,  113,  899,  196,   61,   49,    8,   52],
                [ 407, -192,  196,  611,    8,   44,   59,  -23],
                [-8,   -71,    61,    8,  411, -599,  208,  208],
                [-52,   -43,   49,   44, -599,  411,  208,  208],
                [-49,    -8,    8,   59,  208,  208,   99, -911],
                [29,   -44,    52,  -23,  208,  208, -911,   99]])
    assert M.eigenvals() == {0: 1, 10*sqrt(10405): 1, 100*sqrt(26) + 510: 1,
                             1000: 2, -100*sqrt(26) + 510: 1,
                             -10*sqrt(10405): 1, 1020: 1}


def test_P25():
    MF = N(Matrix([[ 611,  196, -192,  407,   -8,  -52,  -49,   29],
                   [ 196,  899,  113, -192,  -71,  -43,   -8,  -44],
                   [-192,  113,  899,  196,   61,   49,    8,   52],
                   [ 407, -192,  196,  611,    8,   44,   59,  -23],
                   [-8,   -71,    61,    8,  411, -599,  208,  208],
                   [-52,   -43,   49,   44, -599,  411,  208,  208],
                   [-49,    -8,    8,   59,  208,  208,   99, -911],
                   [ 29,   -44,   52,  -23,  208,  208, -911,   99]]))
    assert (Matrix(sorted(MF.eigenvals())) - Matrix(
            [-1020.0490184299969, 0.0, 0.09804864072151699, 1000.0,
             1019.9019513592784, 1020.0, 1020.0490184299969])).norm() < 1e-13


def test_P26():
    a0, a1, a2, a3, a4 = symbols('a0 a1 a2 a3 a4')
    M = Matrix([[-a4, -a3, -a2, -a1, -a0,  0,  0,  0,  0],
                [  1,   0,   0,   0,   0,  0,  0,  0,  0],
                [  0,   1,   0,   0,   0,  0,  0,  0,  0],
                [  0,   0,   1,   0,   0,  0,  0,  0,  0],
                [  0,   0,   0,   1,   0,  0,  0,  0,  0],
                [  0,   0,   0,   0,   0, -1, -1,  0,  0],
                [  0,   0,   0,   0,   0,  1,  0,  0,  0],
                [  0,   0,   0,   0,   0,  0,  1, -1, -1],
                [  0,   0,   0,   0,   0,  0,  0,  1,  0]])
    assert M.eigenvals() == {
        Rational(-1, 2) - sqrt(3)*I/2: 2,
        Rational(-1, 2) + sqrt(3)*I/2: 2}


def test_P27():
    M = Matrix([[a,  0, 0, 0, 0],
                [0,  0, 0, 0, 1],
                [0,  0, a, 0, 0],
                [0,  0, 0, a, 0],
                [0, -2, 0, 0, 2]])
    assert M.eigenvects() == [(a, 3, [Matrix([[1],
                                              [0],
                                              [0],
                                              [0],
                                              [0]]),
                                      Matrix([[0],
                                              [0],
                                              [1],
                                              [0],
                                              [0]]),
                                      Matrix([[0],
                                              [0],
                                              [0],
                                              [1],
                                              [0]])]),
                              (1 - I, 1, [Matrix([[          0],
                                                  [-1/(-1 + I)],
                                                  [          0],
                                                  [          0],
                                                  [          1]])]),
                              (1 + I, 1, [Matrix([[          0],
                                                  [-1/(-1 - I)],
                                                  [          0],
                                                  [          0],
                                                  [          1]])])]


def test_P30():
    M = Matrix([[1,  0,  0,  1, -1],
                [0,  1, -2,  3, -3],
                [0,  0, -1,  2, -2],
                [1, -1,  1,  0,  1],
                [1, -1,  1, -1,  2]])
    J, _ = M.jordan_form()
    assert J == Matrix([[-1, 0, 0, 0, 0],
                        [0,  1, 1, 0, 0],
                        [0,  0, 1, 0, 0],
                        [0,  0, 0, 1, 1],
                        [0,  0, 0, 0, 1]])


def test_P32():
    M = Matrix([[1, -2],
                [2, 1]]).as_immutable()
    assert exp(M).rewrite(cos).simplify() == Matrix([[E*cos(2), -E*sin(2)],
                                                     [E*sin(2),  E*cos(2)]])


def test_P33():
    M = Matrix([[0,    1,      0,   0],
                [0,    0,      0, 2*w],
                [0,    0,      0,   1],
                [0, -2*w, 3*w**2,   0]]).as_immutable()
    assert exp(M*t).rewrite(cos).expand() == Matrix([
        [1, -3*t + 4*sin(t*w)/w,  6*t*w - 6*sin(t*w), -2*cos(t*w)/w + 2/w],
        [0,      4*cos(t*w) - 3, -6*w*cos(t*w) + 6*w,          2*sin(t*w)],
        [0,  2*cos(t*w)/w - 2/w,     -3*cos(t*w) + 4,          sin(t*w)/w],
        [0,         -2*sin(t*w),        3*w*sin(t*w),            cos(t*w)]])


@pytest.mark.xfail
def test_P34():
    a, b, c = symbols('a b c', real=True)
    M = Matrix([[a, 1, 0, 0, 0, 0],
                [0, a, 0, 0, 0, 0],
                [0, 0, b, 0, 0, 0],
                [0, 0, 0, c, 1, 0],
                [0, 0, 0, 0, c, 1],
                [0, 0, 0, 0, 0, c]])
    assert sin(M) == Matrix([[sin(a), cos(a), 0, 0, 0, 0],
                             [0, sin(a), 0, 0, 0, 0],
                             [0, 0, sin(b), 0, 0, 0],
                             [0, 0, 0, sin(c), cos(c), -sin(c)/2],
                             [0, 0, 0, 0, sin(c), cos(c)],
                             [0, 0, 0, 0, 0, sin(c)]])


@pytest.mark.xfail
def test_P35():
    M = pi/2*Matrix([[2, 1, 1],
                     [2, 3, 2],
                     [1, 1, 2]])
    assert sin(M) == eye(3)


@pytest.mark.xfail
def test_P36():
    M = Matrix([[10, 7],
                [7, 17]])
    assert sqrt(M) == Matrix([[3, 1],
                              [1, 4]])


def test_P37():
    M = Matrix([[1, 1, 0],
                [0, 1, 0],
                [0, 0, 1]])
    assert M**Rational(1, 2) == Matrix([[1, 1/2, 0],
                                        [0, 1, 0],
                                        [0, 0, 1]])


def test_P38():
    M = Matrix([[0, 1, 0],
                [0, 0, 0],
                [0, 0, 0]])
    assert all(e in (nan, zoo) for e in M**Rational(1, 2))


def test_P40():
    r, t = symbols('r t', real=True)
    M = Matrix([r*cos(t), r*sin(t)])
    assert M.jacobian(Matrix([r, t])) == Matrix([[cos(t), -r*sin(t)],
                                                 [sin(t),  r*cos(t)]])


def test_P41():
    r, t = symbols('r t', real=True)
    assert hessian(r**2*sin(t), (r, t)) == Matrix([[  2*sin(t),   2*r*cos(t)],
                                                   [2*r*cos(t), -r**2*sin(t)]])


def test_P42():
    assert wronskian([cos(x), sin(x)], x).simplify() == 1


def test_P43():
    def __my_jacobian(M, Y):
        return Matrix([M.diff(v).T for v in Y]).T
    r, t = symbols('r t', real=True)
    M = Matrix([r*cos(t), r*sin(t)])
    assert __my_jacobian(M, [r, t]) == Matrix([[cos(t), -r*sin(t)],
                                               [sin(t),  r*cos(t)]])


def test_P44():
    def __my_hessian(f, Y):
        V = Matrix([diff(f, v) for v in Y])
        return Matrix([V.T.diff(v) for v in Y])
    r, t = symbols('r t', real=True)
    assert __my_hessian(r**2*sin(t), (r, t)) == Matrix([
        [  2*sin(t),   2*r*cos(t)],
        [2*r*cos(t), -r**2*sin(t)]])


def test_P45():
    def __my_wronskian(Y, v):
        M = Matrix([Matrix(Y).T.diff((x, n)) for n in range(len(Y))])
        return M.det()
    assert __my_wronskian([cos(x), sin(x)], x).simplify() == 1


# R. Sums

@pytest.mark.xfail
def test_R1():
    i, j, n = symbols('i j n', integer=True, positive=True)
    xn = MatrixSymbol('xn', n, 1)
    summation((xn[i, 0] - Sum(xn[j, 0], (j, 0, n - 1))/n)**2, (i, 0, n - 1))


@pytest.mark.xfail
def test_R2():
    m, b = symbols('m b', real=True)
    i, n = symbols('i n', integer=True, positive=True)
    xn = MatrixSymbol('xn', n, 1)
    yn = MatrixSymbol('yn', n, 1)
    f = Sum((yn[i, 0] - m*xn[i, 0] - b)**2, (i, 0, n - 1))
    f1 = diff(f, m)
    f2 = diff(f, b)
    assert solve((f1, f2), m, b) != []


def test_R3():
    n, k = symbols('n k', integer=True, positive=True)
    sk = ((-1)**k) * (binomial(2*n, k))**2
    T = summation(sk, (k, 0, oo))
    assert simplify(T - (-1)**n*binomial(2*n, n)) == 0


def test_R6():
    n, k = symbols('n k', integer=True, positive=True)
    gn = MatrixSymbol('gn', n + 2, 1)
    T = summation(gn[k, 0] - gn[k - 1, 0], (k, 1, n + 1))
    assert T == -gn[0, 0] + gn[n + 1, 0]


def test_R7():
    n, k = symbols('n k', integer=True, positive=True)
    T = Sum(k**3, (k, 1, n)).doit()
    assert T.factor() == n**2*(n + 1)**2/4


@pytest.mark.xfail
def test_R8():
    n, k = symbols('n k', integer=True, positive=True)
    T = summation(k**2*binomial(n, k), (k, 1, n))
    assert simplify(T - n*(n + 1)*2**(n - 2)) == 0


def test_R9():
    n, k = symbols('n k', integer=True, positive=True)
    Sm = summation(binomial(n, k - 1)/k, (k, 1, n + 1))
    assert Sm.simplify() == (2**(n + 1) - 1)/(n + 1)


def test_R10():
    n, m, r, k = symbols('n m r k', integer=True, positive=True)
    T = summation(binomial(n, k)*binomial(m, r - k), (k, 0, r))
    assert (T - binomial(m + n, r)).rewrite(factorial) == 0


@pytest.mark.xfail
def test_R11():
    n, k = symbols('n k', integer=True, positive=True)
    sk = binomial(n, k)*fibonacci(k)
    T = summation(sk, (k, 0, n))
    assert simplify(T - fibonacci(2*n)) == 0


@pytest.mark.xfail
def test_R12():
    n, k = symbols('n k', integer=True, positive=True)
    T = summation(fibonacci(k)**2, (k, 0, n))
    assert simplify(T - fibonacci(n)*fibonacci(n + 1)) == 0


@pytest.mark.xfail
def test_R13():
    n, k = symbols('n k', integer=True, positive=True)
    T = summation(sin(k*x), (k, 1, n))
    assert simplify(T - (cot(x/2)/2 - cos(x*(2*n + 1)/2)/(2*sin(x/2)))) == 0


@pytest.mark.xfail
def test_R14():
    n, k = symbols('n k', integer=True, positive=True)
    T = summation(sin((2*k - 1)*x), (k, 1, n))
    assert simplify(T - sin(n*x)**2/sin(x)) == 0


def test_R16():
    k = symbols('k', integer=True, positive=True)
    Sm = Sum(1/k**2 + 1/k**3, (k, 1, oo))
    assert Sm.doit() == zeta(3) + pi**2/6


def test_R17():
    k = symbols('k', integer=True, positive=True)
    assert abs(float(Sum(1/k**2 + 1/k**3, (k, 1, oo)))
               - 2.8469909700078206) < 1e-15


def test_R18():
    k = symbols('k', integer=True, positive=True)
    assert summation(1/(2**k*k**2), (k, 1, oo)) == polylog(2, Rational(1, 2))


@pytest.mark.slow
@pytest.mark.xfail
def test_R19():
    k = symbols('k', integer=True, positive=True)
    s = summation(1/((3*k + 1)*(3*k + 2)*(3*k + 3)), (k, 0, oo))
    assert simplify(s + log(3)/4 - sqrt(3)*pi/12) == 0


@pytest.mark.xfail
def test_R20():
    n, k = symbols('n k', integer=True, positive=True)
    T = summation(binomial(n, 4*k), (k, 0, oo))
    assert simplify(T - 2**(n/2)*cos(pi*n/4)/2 - 2**(n - 1)/2) == 0


@pytest.mark.xfail
def test_R21():
    k = symbols('k', integer=True, positive=True)
    T = summation(1/(sqrt(k*(k + 1)) * (sqrt(k) + sqrt(k + 1))), (k, 1, oo))
    assert simplify(T) == 1


def test_R24():
    m, k = symbols('m k', integer=True, positive=True)
    Sm = Sum(Product(k/(2*k - 1), (k, 1, m)), (m, 2, oo))
    assert Sm.doit() == pi/2


# S. Products

def test_S1():
    k = symbols('k', integer=True, positive=True)
    Pr = product(gamma(k/3), (k, 1, 8))
    assert Pr.simplify() == 640*sqrt(3)*pi**3/6561


def test_S2():
    n, k = symbols('n k', integer=True, positive=True)
    assert Product(k, (k, 1, n)).doit() == factorial(n)


def test_S3():
    n, k = symbols('n k', integer=True, positive=True)
    assert product(x**k, (k, 1, n)).simplify() == x**(n*(n + 1)/2)


def test_S4():
    n, k = symbols('n k', integer=True, positive=True)
    assert product(1 + 1/k, (k, 1, n - 1)).simplify() == n


def test_S5():
    n, k = symbols('n k', integer=True, positive=True)
    assert (Product((2*k - 1)/(2*k), (k, 1, n)).doit().combsimp() ==
            factorial(n - Rational(1, 2))/(sqrt(pi)*factorial(n)))


@pytest.mark.xfail(reason='https://github.com/sympy/sympy/issues/7133')
def test_S6():
    n, k = symbols('n k', integer=True, positive=True)
    assert simplify(product(x**2 - 2*x*cos(k*pi/n) + 1, (k, 1, n - 1))
                    - (x**(2*n) - 1)/(x**2 - 1)) == 0


@pytest.mark.xfail
def test_S7():
    k = symbols('k', integer=True, positive=True)
    T = product((k**3 - 1)/(k**3 + 1), (k, 2, oo))
    assert simplify(T) == Rational(2, 3)


@pytest.mark.xfail
def test_S8():
    k = symbols('k', integer=True, positive=True)
    T = product(1 - 1/(2*k)**2, (k, 1, oo))
    assert simplify(T) == 2/pi


@pytest.mark.xfail(reason='https://github.com/sympy/sympy/issues/7133')
def test_S9():
    k = symbols('k', integer=True, positive=True)
    T = product(1 + (-1)**(k + 1)/(2*k - 1), (k, 1, oo))
    assert simplify(T) == sqrt(2)


@pytest.mark.xfail(reason='https://github.com/sympy/sympy/issues/7137')
def test_S10():
    k = symbols('k', integer=True, positive=True)
    T = product((k*(k + 1) + 1 + I)/(k*(k + 1) + 1 - I), (k, 0, oo))
    assert simplify(T) == -1


# T. Limits

def test_T1():
    assert limit((1 + 1/n)**n, n, oo) == E
    assert limit((1 - cos(x))/x**2, x, 0) == Rational(1, 2)


def test_T2():
    assert limit((3**x + 5**x)**(1/x), x, oo) == 5


@pytest.mark.xfail
def test_T3():
    assert limit(log(x)/(log(x) + sin(x)), x, oo) == 1


def test_T4():
    assert limit((exp(x*exp(-x)/(exp(-x) + exp(-2*x**2/(x + 1))))
                  - exp(x))/x, x, oo) == -exp(2)


def test_T5():
    assert limit(x*log(x)*log(x*exp(x) - x**2)**2/log(log(x**2
                                                          + 2*exp(exp(3*x**3*log(x))))), x, oo) == Rational(1, 3)


def test_T6():
    assert limit(1/n * factorial(n)**(1/n), n, oo) == exp(-1)


def test_T7():
    limit(1/n * gamma(n + 1)**(1/n), n, oo)


def test_T8():
    a, z = symbols('a z', real=True, positive=True)
    assert limit(gamma(z + a)/gamma(z)*exp(-a*log(z)), z, oo) == 1


@pytest.mark.xfail
def test_T9():
    z, k = symbols('z k', real=True, positive=True)
    assert limit(hyper((1, k), (1,), z/k), k, oo) == exp(z)


@pytest.mark.xfail
def test_T10():
    assert limit(zeta(x) - 1/(x - 1), x, 1) == EulerGamma


@pytest.mark.xfail
def test_T11():
    n, k = symbols('n k', integer=True, positive=True)
    assert limit(n**x/(x*product((1 + x/k), (k, 1, n))), n, oo) == gamma(x)


def test_T12():
    x, t = symbols('x t', real=True)
    assert limit(x*integrate(exp(-t**2), (t, 0, x))/(1 - exp(-x**2)),
                 x, 0) == 1


def test_T13():
    assert [limit(x/abs(x), x, 0, dir=1),
            limit(x/abs(x), x, 0)] == [-1, 1]


def test_T14():
    assert limit(atan(-log(x)), x, 0) == pi/2


# U. Calculus

def test_U1():
    x = symbols('x', real=True)
    assert diff(abs(x), x) == sign(x)


def test_U2():
    f = Lambda(x, Piecewise((-x, x < 0), (x, True)))
    assert diff(f(x), x) == Piecewise((-1, x < 0), (1, True))


def test_U3():
    f = Lambda(x, Piecewise((x**2 - 1, x == 1), (x**3, x != 1)))
    f1 = Lambda(x, diff(f(x), x))
    assert f1(x) == 3*x**2
    assert f1(1) == 3


def test_U6():
    h = Function('h')
    T = integrate(f(y), y, h(x), g(x))
    T.diff(x)


@pytest.mark.xfail
def test_U7():
    p, t = symbols('p t', real=True)
    diff(f(p, t))


def test_U8():
    x, y = symbols('x y', real=True)
    eq = cos(x*y) + x
    eq = eq.subs({y: f(x)})
    #  If Diofant had implicit_diff() function this hack could be avoided
    assert (solve((f(x) - eq).diff(x),
                  f(x).diff(x))[0][f(x).diff(x)].subs({f(x): y}) ==
            (-y*sin(x*y) + 1)/(x*sin(x*y) + 1))


def test_U9():
    su = diff(f(x, y), x) + diff(f(x, y), y)
    s2 = su.subs({f(x, y): g(x**2 + y**2)})
    s3 = s2.doit()
    assert s3 == (2*x*Subs(Derivative(g(z), z), (z, x**2 + y**2)) +
                  2*y*Subs(Derivative(g(z), z), (z, x**2 + y**2)))


def test_U10():
    assert residue((z**3 + 5)/((z**4 - 1)*(z + 1)), z, -1) == Rational(-9, 4)


def test_U13():
    assert minimize(x**4 - x + 1, x)[0] == -3*cbrt(2)/8 + 1


@pytest.mark.xfail
def test_U14():
    f = 1/(x**2 + y**2 + 1)
    assert [minimize(f, x, y)[0], maximize(f, x, y)[0]] == [0, 1]


def test_U17():
    x1, x2, x3, x4 = symbols('x1:5')
    assert minimize([4*x1 - x2 + 2*x3 - 2*x4,
                     2*x1 + x2 + x3 + x4 <= 10,
                     x1 - 2*x2 - x3 + x4 >= 4,
                     x1 + x2 + 3*x3 - x4 >= 4,
                     x1 >= 0, x2 >= 0, x3 >= 0,
                     x4 >= 0], x1, x2,
                    x3, x4) == (4, {x1: 2, x2: 0, x3: 2, x4: 4})


# V. Indefinite Integrals

@pytest.mark.xfail
def test_V1():
    x = symbols('x', real=True)
    assert integrate(abs(x), x) == x*abs(x)/2


def test_V2():
    assert (integrate(Piecewise((-x, x < 0), (x, x >= 0)), x) ==
            Piecewise((-x**2/2, x < 0), (x**2/2, x >= 0)))


def test_V3():
    assert integrate(1/(x**3 + 2), x).diff().simplify() == 1/(x**3 + 2)


@pytest.mark.xfail
def test_V4():
    assert simplify(integrate(2**x/sqrt(1 + 4**x), x) - asinh(2**x)/log(2)) == 0


@pytest.mark.slow
def test_V5():
    assert simplify(integrate((3*x - 5)**2/sqrt(2*x - 1)**7, x) -
                    (-41 + 80*x - 45*x**2)/(5*sqrt(2*x - 1)**5)) == 0


@pytest.mark.xfail
def test_V6():
    assert simplify(integrate(1/(2*exp(m*x) - 5*exp(-m*x)), x) +
                    atanh(exp(m*x)*sqrt(Rational(2, 5)))/(m*sqrt(10))) == 0


def test_V7():
    r1 = integrate(sinh(x)**4/cosh(x)**2)
    assert r1.simplify() == -3*x/2 + sinh(x)**3/(2*cosh(x)) + 3*tanh(x)/2


def test_V10():
    assert integrate(1/(3 + 3*cos(x) + 4*sin(x)), x) == log(tan(x/2) + Rational(3, 4))/4


def test_V11():
    r1 = integrate(1/(4 + 3*cos(x) + 4*sin(x)), x)
    r2 = factor(r1)
    assert logcombine(r2, force=True) == log(cbrt((tan(x/2) + 1)/(tan(x/2) + 7)))


def test_V12():
    r1 = integrate(1/(5 + 3*cos(x) + 4*sin(x)), x)
    assert r1 == -1/(tan(x/2) + 2)


@pytest.mark.slow
@pytest.mark.xfail
def test_V13():
    r1 = integrate(1/(6 + 3*cos(x) + 4*sin(x)), x)
    assert simplify(r1 - 2*sqrt(11)*atan(sqrt(11)*(3*tan(x/2) + 4)/11)/11) == 0


def test_V15():
    r1 = integrate(x*acot(x/y), x)
    assert simplify(r1 - (x*y + (x**2 + y**2)*acot(x/y))/2) == 0


@pytest.mark.slow
@pytest.mark.xfail
def test_V17():
    r1 = integrate((diff(f(x), x)*g(x)
                    - f(x)*diff(g(x), x))/(f(x)**2 - g(x)**2), x)
    assert simplify(r1 - log((f(x) - g(x))/(f(x) + g(x)))/2) == 0


# W. Definite Integrals

@pytest.mark.xfail
@pytest.mark.slow
def test_W3():
    assert simplify(integrate(sqrt(x + 1/x - 2), (x, 0, 1)) - Rational(4, 3)) == 0


@pytest.mark.xfail
def test_W4():
    assert integrate(sqrt(x + 1/x - 2), (x, 1, 2)) == -2*sqrt(2)/3 + Rational(4, 3)


@pytest.mark.xfail
def test_W5():
    assert integrate(sqrt(x + 1/x - 2), (x, 0, 2)) == -2*sqrt(2)/3 + Rational(8, 3)


def test_W7():
    a = symbols('a', real=True, positive=True)
    r1 = integrate(cos(x)/(x**2 + a**2), (x, -oo, oo))
    assert r1.simplify() == pi*exp(-a)/a


@pytest.mark.xfail
@pytest.mark.slow
def test_W10():
    r1 = integrate(x/(1 + x + x**2 + x**4), (x, -oo, oo))
    assert not isinstance(r1, Integral)


@pytest.mark.xfail
@pytest.mark.slow
def test_W11():
    assert simplify(integrate(sqrt(1 - x**2)/(1 + x**2), (x, -1, 1)) -
                    pi*(-1 + sqrt(2))) == 0


def test_W12():
    p = symbols('p', real=True, positive=True)
    q = symbols('q', real=True)
    r1 = integrate(x*exp(-p*x**2 + 2*q*x), (x, -oo, oo))
    assert simplify(r1) == sqrt(pi)*q*exp(q**2/p)/sqrt(p)**3


def test_W14():
    assert integrate(sin(x)/x*exp(2*I*x), (x, -oo, oo)) == 0


@pytest.mark.xfail
@pytest.mark.slow
def test_W15():
    assert simplify(integrate(log(gamma(x))*cos(6*pi*x), (x, 0, 1)) -
                    Rational(1, 12)) == 0


def test_W16():
    assert integrate((1 + x)**3*legendre_poly(1, x)*legendre_poly(2, x),
                     (x, -1, 1)) == Rational(36, 35)


def test_W17():
    a, b = symbols('a b', real=True, positive=True)
    assert integrate(exp(-a*x)*besselj(0, b*x),
                     (x, 0, oo)) == 1/(b*sqrt(a**2/b**2 + 1))


def test_W18():
    assert integrate((besselj(1, x)/x)**2, (x, 0, oo)) == 4/(3*pi)


@pytest.mark.xfail
def test_W20():
    assert not integrate(x**2*polylog(3, 1/(x + 1)), (x, 0, 1)).has(Integral)


def test_W21():
    assert abs(N(integrate(x**2*polylog(3, 1/(x + 1)), (x, 0, 1)))
               - 0.210882859565594) < 1e-15


def test_W22():
    t, u = symbols('t u', real=True)
    s = Lambda(x, Piecewise((1, And(x >= 1, x <= 2)), (0, True)))
    assert (integrate(s(t)*cos(t), (t, 0, u)) ==
            Piecewise((sin(u) - sin(1), And(u <= 2, u >= 1)),
                      (0, And(u <= 1, u >= -oo)),
                      (-sin(1) + sin(2), True)))


def test_W23():
    a, b = symbols('a b', real=True, positive=True)
    r1 = integrate(integrate(x/(x**2 + y**2), (x, a, b)), (y, -oo, oo))
    assert r1.simplify() == pi*(-a + b)


@pytest.mark.xfail
@pytest.mark.slow
def test_W23b():
    a, b = symbols('a b', real=True, positive=True)
    r2 = integrate(integrate(x/(x**2 + y**2), (y, -oo, oo)), (x, a, b))
    assert simplify(r2 - pi*(-a + b)) == 0


@pytest.mark.xfail
@pytest.mark.slow
def test_W24():
    x, y = symbols('x y', real=True)
    r1 = integrate(integrate(sqrt(x**2 + y**2), (x, 0, 1)), (y, 0, 1))
    assert simplify(r1 - (sqrt(2) + asinh(1))/3) == 0


@pytest.mark.slow
def test_W26():
    x, y = symbols('x y', real=True)
    assert integrate(integrate(abs(y - x**2), (y, 0, 2)),
                     (x, -1, 1)) == Rational(46, 15)


def test_W27():
    assert integrate(integrate(integrate(1, (z, 0, c*(1 - x/a - y/b))),
                               (y, 0, b*(1 - x/a))),
                     (x, 0, a)) == a*b*c/6


# X. Series

def test_X1():
    assert ((1/sqrt(1 - (v/c)**2)).series(v, n=8) ==
            5*v**6/(16*c**6) + 3*v**4/(8*c**4) + v**2/(2*c**2) + 1 + O(v**8))


def test_X2():
    s1 = (1/sqrt(1 - (v/c)**2)).series(v, n=8)
    assert (1/s1**2).series(v, x0=0, n=8) == -v**2/c**2 + 1 + O(v**8)


def test_X3():
    s1 = (sin(x).series()/cos(x).series()).series()
    s2 = tan(x).series()
    assert s1 == s2 == x + x**3/3 + 2*x**5/15 + O(x**6)


def test_X4():
    s1 = log(sin(x)/x).series()
    s2 = log((sin(x)/x).series()).series()
    assert s1 == s2 == -x**2/6 - x**4/180 + O(x**6)


def test_X5():
    assert ((diff(f(a*x), x) + g(b*x) +
             integrate(h(c*y), (y, 0, x))).series(x, x0=d, n=2) ==
            (Integral(h(c*y), (y, 0, x)) + g(b*d) +
             b*(-d + x)*Subs(Derivative(g(y), y), (y, b*d)) +
             a*Subs(Derivative(f(y), y), (y, a*x)) + O((-d + x)**2, (x, d))))


def test_X6():
    # Taylor series of nonscalar objects (noncommutative multiplication)
    # expected result => (B A - A B) t^2/2 + O(t^3)   [Stanly Steinberg]
    a, b = symbols('a b', commutative=False)
    assert ((exp((a + b)*x) - exp(a*x)*exp(b*x)).series(x, n=3) ==
            x**2*(b*a/2 - a*b/2) + O(x**3))


@pytest.mark.slow
def test_X7():
    # => sum( Bernoulli[k]/k! x^(k - 2), k = 1..infinity )
    #    = 1/x^2 - 1/(2 x) + 1/12 - x^2/720 + x^4/30240 + O(x^6)
    #    [Levinson and Redheffer, p. 173]
    assert ((1/(x*(exp(x) - 1))).series(x, n=7) == x**(-2) - 1/(2*x) +
            Rational(1, 12) - x**2/720 + x**4/30240 - x**6/1209600 + O(x**7))


def test_X8():
    # Puiseux series (terms with fractional degree):
    # => 1/sqrt(x - 3/2 pi) + (x - 3/2 pi)^(3/2) / 12 + O([x - 3/2 pi]^(7/2))

    # see issue sympy/sympy#7167
    assert (sqrt(sec(x)).series(x, x0=pi*3/2, n=4) ==
            1/sqrt(x - 3*pi/2) + sqrt(x - 3*pi/2)**3/12 +
            sqrt(x - 3*pi/2)**7/160 + O((x - 3*pi/2)**4, (x, 3*pi/2)))


def test_X9():
    assert ((x**x).series(x, n=4) == 1 + x*log(x) + x**2*log(x)**2/2 +
            x**3*log(x)**3/6 + x**4*log(x)**4) + O(x**4)


def test_X10():
    assert ((log(sinh(z)) + log(cosh(z + w))).series(z, n=2) ==
            log(cosh(w)) + log(z) + z*sinh(w)/cosh(w) + O(z**2))


def test_X11():
    assert (log(sinh(z)*cosh(z + w)).series(z, n=2) ==
            log(cosh(w)) + log(z) + z*sinh(w)/cosh(w) + O(z**2))


def test_X13():
    assert sqrt(2*x**2 + 1).series(x, x0=oo, n=1) == sqrt(2)*x + O(1/x, (x, oo))


@pytest.mark.xfail
def test_X14():
    (1/2**(2*n)*binomial(2*n, n)).series(n, x0=oo, n=1)


@pytest.mark.xfail(reason='https://github.com/sympy/sympy/issues/7164')
def test_X15():
    x, t = symbols('x t', real=True)
    e1 = integrate(exp(-t)/t, (t, x, oo))
    e1.series(x, x0=oo, n=5)


@pytest.mark.xfail(reason='https://github.com/diofant/diofant/pull/158')
def test_X16():
    cos(x + y).series(x + y, x0=0, n=4)


# Y. Transforms

def test_Y1():
    t = symbols('t', real=True, positive=True)
    w = symbols('w', real=True)
    F, _, _ = laplace_transform(cos((w - 1)*t), t, s)
    assert F == s/(s**2 + (w - 1)**2)


def test_Y2():
    t = symbols('t', real=True, positive=True)
    w = symbols('w', real=True)
    f = inverse_laplace_transform(s/(s**2 + (w - 1)**2), s, t)
    assert f == cos(t*abs(w - 1))


@pytest.mark.slow
def test_Y3():
    t = symbols('t', real=True, positive=True)
    w = symbols('w', real=True)
    F, _, _ = laplace_transform(sinh(w*t)*cosh(w*t), t, s)
    assert F == w/(s**2 - 4*w**2)


def test_Y4():
    t = symbols('t', real=True, positive=True)
    F, _, _ = laplace_transform(erf(3/sqrt(t)), t, s)
    assert F == (1 - exp(-6*sqrt(s)))/s


@pytest.mark.xfail
def test_Y5():
    t = symbols('t', real=True, positive=True)
    y = Function('y')
    F, _, _ = laplace_transform(diff(y(t), t, 2) + y(t) -
                                4*(Heaviside(t - 1) - Heaviside(t - 2)), t, s)
    assert (F == s**2*LaplaceTransform(y(t), t, s) - s
            + LaplaceTransform(y(t), t, s) - 4*exp(-s)/s + 4*exp(-2*s)/s)


@pytest.mark.xfail
def test_Y8():
    assert fourier_transform(1, x, z) != 0


def test_Y9():
    assert (fourier_transform(exp(-9*x**2), x, z) ==
            sqrt(pi)*exp(-pi**2*z**2/9)/3)


def test_Y10():
    assert (fourier_transform(abs(x)*exp(-3*abs(x)), x, z) ==
            2*((-4*pi**2*z**2 + 9)/(16*pi**4*z**4 + 72*pi**2*z**2 + 81)))


def test_Y12():
    F, _, _ = mellin_transform(besselj(3, x)/x**3, x, s)
    assert simplify(F - (2**(s - 4)*gamma(s/2))/gamma(4 - s/2)) == 0


# Z. Ordinary Difference and Differential Equations

def test_Z1():
    assert (rsolve(f(n + 2) - 2*f(n + 1) + f(n) - 2,
                   init={f(0): 1, f(1): m}) == [{f: Lambda(n, n**2 +
                                                           n*(m - 2) + 1)}])


def test_Z2():
    assert rsolve(f(n) - (5*f(n - 1) - 6*f(n - 2)),
                  init={f(0): 0, f(1): 1}) == [{f: Lambda(n, -2**n + 3**n)}]


@pytest.mark.slow
def test_Z3():
    # => f(n) = Fibonacci[n + 1]   [Cohen, p. 83]
    # recurrence solution is correct, Wester expects it to be simplified to
    # fibonacci(n+1), but that is quite hard
    assert (rsolve(f(n) - (f(n - 1) + f(n - 2)),
                   init={f(1): 1, f(2): 2})
            == [{f: Lambda(n, 2**(-n)*((1 + sqrt(5))**n*(sqrt(5) + 5) +
                           (-sqrt(5) + 1)**n*(-sqrt(5) + 5))/10)}])


@pytest.mark.xfail
def test_Z4():
    rsolve(f(n) - ((1 + c - c**(n-1) - c**(n+1))/(1 - c**n)*f(n - 1)
                   - c*(1 - c**(n-2))/(1 - c**(n-1))*f(n - 2) + 1),
           init={f(1): 1, f(2): (2 + 2*c + c**2)/(1 + c)})


def test_Z5():
    eq = Derivative(f(x), (x, 2)) + 4*f(x) - sin(2*x)
    sol = dsolve(eq, f(x), init={f(0): 0, f(x).diff(x).subs({x: 0}): 0})
    assert solve(sol, f(x))[0][f(x)] == -x*cos(2*x)/4 + sin(2*x)/8


@pytest.mark.xfail
def test_Z6():
    t = symbols('t', real=True, positive=True)
    eq = Derivative(f(t), t, 2) + 4*f(t) - sin(2*t)
    F, _, _ = laplace_transform(eq, t, s)
    assert (F == s**2*LaplaceTransform(f(t), t, s) +
            4*LaplaceTransform(f(t), t, s) - 2/(s**2 + 4))
