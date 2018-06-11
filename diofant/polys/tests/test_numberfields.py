"""Tests for computational algebraic number field theory. """

import pytest

from diofant import (Add, GoldenRatio, I, Integer, Mul, Poly, PurePoly,
                     Rational, Tuple, cbrt, cos, exp, exp_polar, expand,
                     expand_multinomial, nsimplify, oo, pi, root, sin, solve,
                     sqrt)
from diofant.abc import x, y, z
from diofant.domains import QQ
from diofant.polys.numberfields import (AlgebraicNumber, field_isomorphism,
                                        is_isomorphism_possible,
                                        minimal_polynomial, primitive_element)
from diofant.polys.polyclasses import DMP
from diofant.polys.polyerrors import CoercionFailed, NotAlgebraic
from diofant.polys.polytools import degree
from diofant.polys.rootoftools import RootOf


__all__ = ()


def test_minimal_polynomial():
    assert minimal_polynomial(-7)(x) == x + 7
    assert minimal_polynomial(-1)(x) == x + 1
    assert minimal_polynomial( 0)(x) == x
    assert minimal_polynomial( 1)(x) == x - 1
    assert minimal_polynomial( 7)(x) == x - 7

    assert minimal_polynomial(Rational(1, 3), method='groebner')(x) == 3*x - 1

    pytest.raises(NotAlgebraic,
                  lambda: minimal_polynomial(pi, method='groebner'))
    pytest.raises(NotAlgebraic,
                  lambda: minimal_polynomial(sin(sqrt(2)), method='groebner'))
    pytest.raises(NotAlgebraic,
                  lambda: minimal_polynomial(2**pi, method='groebner'))

    pytest.raises(ValueError, lambda: minimal_polynomial(1, method='spam'))

    assert minimal_polynomial(sqrt(2))(x) == x**2 - 2
    assert minimal_polynomial(sqrt(5))(x) == x**2 - 5
    assert minimal_polynomial(sqrt(6))(x) == x**2 - 6

    assert minimal_polynomial(2*sqrt(2))(x) == x**2 - 8
    assert minimal_polynomial(3*sqrt(5))(x) == x**2 - 45
    assert minimal_polynomial(4*sqrt(6))(x) == x**2 - 96

    assert minimal_polynomial(2*sqrt(2) + 3)(x) == x**2 - 6*x + 1
    assert minimal_polynomial(3*sqrt(5) + 6)(x) == x**2 - 12*x - 9
    assert minimal_polynomial(4*sqrt(6) + 7)(x) == x**2 - 14*x - 47

    assert minimal_polynomial(2*sqrt(2) - 3)(x) == x**2 + 6*x + 1
    assert minimal_polynomial(3*sqrt(5) - 6)(x) == x**2 + 12*x - 9
    assert minimal_polynomial(4*sqrt(6) - 7)(x) == x**2 + 14*x - 47

    assert minimal_polynomial(sqrt(1 + sqrt(6)))(x) == x**4 - 2*x**2 - 5
    assert minimal_polynomial(sqrt(I + sqrt(6)))(x) == x**8 - 10*x**4 + 49

    assert minimal_polynomial(2*I + sqrt(2 + I))(x) == x**4 + 4*x**2 + 8*x + 37

    assert minimal_polynomial(sqrt(2) + sqrt(3))(x) == x**4 - 10*x**2 + 1
    assert minimal_polynomial(sqrt(2) + sqrt(3) + sqrt(6))(x) == x**4 - 22*x**2 - 48*x - 23

    e = 1/sqrt(sqrt(1 + sqrt(3)) - 4)
    assert minimal_polynomial(e) == minimal_polynomial(e, method='groebner')
    assert minimal_polynomial(e)(x) == (222*x**8 + 240*x**6 +
                                        94*x**4 + 16*x**2 + 1)

    a = 1 - 9*sqrt(2) + 7*sqrt(3)

    assert minimal_polynomial(1/a)(x) == 392*x**4 - 1232*x**3 + 612*x**2 + 4*x - 1
    assert minimal_polynomial(1/sqrt(a))(x) == 392*x**8 - 1232*x**6 + 612*x**4 + 4*x**2 - 1

    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(oo))
    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(2**y))
    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(sin(1)))

    assert minimal_polynomial(sqrt(2))(x) == x**2 - 2

    assert minimal_polynomial(sqrt(2)) == PurePoly(x**2 - 2)
    assert minimal_polynomial(sqrt(2), method='groebner') == PurePoly(x**2 - 2)

    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(3))

    assert minimal_polynomial(a)(x) == x**2 - 2
    assert minimal_polynomial(b)(x) == x**2 - 3

    assert minimal_polynomial(a) == PurePoly(x**2 - 2)
    assert minimal_polynomial(b) == PurePoly(x**2 - 3)

    assert minimal_polynomial(sqrt(a)) == PurePoly(x**4 - 2)
    assert minimal_polynomial(a + 1) == PurePoly(x**2 - 2*x - 1)
    assert minimal_polynomial(sqrt(a/2 + 17))(x) == 2*x**4 - 68*x**2 + 577
    assert minimal_polynomial(sqrt(b/2 + 17))(x) == 4*x**4 - 136*x**2 + 1153

    # issue diofant/diofant#431
    theta = AlgebraicNumber(sqrt(2), (Rational(1, 2), 17))
    assert minimal_polynomial(theta)(x) == 2*x**2 - 68*x + 577

    theta = AlgebraicNumber(RootOf(x**7 + x - 1, 3), (1, 2, 0, 0, 1))
    ans = minimal_polynomial(theta)(x)
    assert ans == (x**7 - 7*x**6 + 19*x**5 - 27*x**4 + 63*x**3 -
                   115*x**2 + 82*x - 147)
    assert minimal_polynomial(theta.as_expr(), method='groebner')(x) == ans
    theta = AlgebraicNumber(RootOf(x**5 + 5*x - 1, 2), (1, -1, 1))
    ans = (x**30 - 15*x**28 - 10*x**27 + 135*x**26 + 330*x**25 - 705*x**24 -
           150*x**23 + 3165*x**22 - 6850*x**21 + 7182*x**20 + 3900*x**19 +
           4435*x**18 + 11970*x**17 - 259725*x**16 - 18002*x**15 +
           808215*x**14 - 200310*x**13 - 647115*x**12 + 299280*x**11 -
           1999332*x**10 + 910120*x**9 + 2273040*x**8 - 5560320*x**7 +
           5302000*x**6 - 2405376*x**5 + 1016640*x**4 - 804480*x**3 +
           257280*x**2 - 53760*x + 1280)
    assert minimal_polynomial(sqrt(theta) + root(theta, 3))(x) == ans
    theta = sqrt(1 + 1/(AlgebraicNumber(RootOf(x**3 + 4*x - 15, 1),
                                        (1, 0, 1)) +
                        1/(sqrt(3) +
                           AlgebraicNumber(RootOf(x**3 - x + 1, 0),
                                           (1, 2, -1)))))
    ans = (2262264837876687263*x**36 - 38939909597855051866*x**34 +
           315720420314462950715*x**32 - 1601958657418182606114*x**30 +
           5699493671077371036494*x**28 - 15096777696140985506150*x**26 +
           30847690820556462893974*x**24 - 49706549068200640994022*x**22 +
           64013601241426223813103*x**20 - 66358713088213594372990*x**18 +
           55482571280904904971976*x**16 - 37309340229165533529076*x**14 +
           20016999328983554519040*x**12 - 8446273798231518826782*x**10 +
           2738866994867366499481*x**8 - 657825125060873756424*x**6 +
           110036313740049140508*x**4 - 11416087328869938298*x**2 +
           551322649782053543)
    assert minimal_polynomial(theta)(x) == ans

    a, b = sqrt(2)/3 + 7, AlgebraicNumber(sqrt(2)/3 + 7)

    f = 81*x**8 - 2268*x**6 - 4536*x**5 + 22644*x**4 + 63216*x**3 - \
        31608*x**2 - 189648*x + 141358

    assert minimal_polynomial(sqrt(a) + sqrt(sqrt(a)))(x) == f
    assert minimal_polynomial(sqrt(b) + sqrt(sqrt(b)))(x) == f

    assert minimal_polynomial(a**Rational(3, 2))(x) == 729*x**4 - 506898*x**2 + 84604519

    a = AlgebraicNumber(RootOf(x**3 + x - 1, 0))
    assert minimal_polynomial(1/a**2)(x) == x**3 - x**2 - 2*x - 1

    # issue sympy/sympy#5994
    eq = (-1/(800*sqrt(Rational(-1, 240) + 1/(18000*cbrt(Rational(-1, 17280000) +
                                                         sqrt(15)*I/28800000)) + 2*cbrt(Rational(-1, 17280000) +
                                                                                        sqrt(15)*I/28800000))))
    assert minimal_polynomial(eq)(x) == 8000*x**2 - 1

    ex = 1 + sqrt(2) + sqrt(3)
    mp = minimal_polynomial(ex)(x)
    assert mp == x**4 - 4*x**3 - 4*x**2 + 16*x - 8

    ex = 1/(1 + sqrt(2) + sqrt(3))
    mp = minimal_polynomial(ex)(x)
    assert mp == 8*x**4 - 16*x**3 + 4*x**2 + 4*x - 1

    p = cbrt(expand((1 + sqrt(2) - 2*sqrt(3) + sqrt(7))**3))
    mp = minimal_polynomial(p)(x)
    assert mp == x**8 - 8*x**7 - 56*x**6 + 448*x**5 + 480*x**4 - 5056*x**3 + 1984*x**2 + 7424*x - 3008
    p = expand((1 + sqrt(2) - 2*sqrt(3) + sqrt(7))**3)
    mp = minimal_polynomial(p)(x)
    assert mp == x**8 - 512*x**7 - 118208*x**6 + 31131136*x**5 + 647362560*x**4 - 56026611712*x**3 + 116994310144*x**2 + 404854931456*x - 27216576512

    assert minimal_polynomial(-sqrt(5)/2 - Rational(1, 2) + (-sqrt(5)/2 - Rational(1, 2))**2)(x) == x - 1
    a = 1 + sqrt(2)
    assert minimal_polynomial((a*sqrt(2) + a)**3)(x) == x**2 - 198*x + 1

    p = 1/(1 + sqrt(2) + sqrt(3))
    assert minimal_polynomial(p, method='groebner')(x) == 8*x**4 - 16*x**3 + 4*x**2 + 4*x - 1

    p = 2/(1 + sqrt(2) + sqrt(3))
    assert minimal_polynomial(p, method='groebner')(x) == x**4 - 4*x**3 + 2*x**2 + 4*x - 2

    assert minimal_polynomial(1 + sqrt(2)*I, method='groebner')(x) == x**2 - 2*x + 3
    assert minimal_polynomial(1/(1 + sqrt(2)) + 1, method='groebner')(x) == x**2 - 2
    assert minimal_polynomial(sqrt(2)*I + I*(1 + sqrt(2)),
                              method='groebner')(x) == x**4 + 18*x**2 + 49

    assert minimal_polynomial(exp_polar(0))(x) == x - 1


def test_minimal_polynomial_hi_prec():
    p = 1/sqrt(1 - 9*sqrt(2) + 7*sqrt(3) + Rational(1, 10)**30)
    mp = minimal_polynomial(p)(x)
    # checked with Wolfram Alpha
    assert mp.coeff(x**6) == -1232000000000000000000000000001223999999999999999999999999999987999999999999999999999999999996000000000000000000000000000000


def test_minimal_polynomial_sq():
    p = expand_multinomial((1 + 5*sqrt(2) + 2*sqrt(3))**3)
    mp = minimal_polynomial(cbrt(p))(x)
    assert mp == x**4 - 4*x**3 - 118*x**2 + 244*x + 1321
    p = expand_multinomial((1 + sqrt(2) - 2*sqrt(3) + sqrt(7))**3)
    mp = minimal_polynomial(cbrt(p))(x)
    assert mp == x**8 - 8*x**7 - 56*x**6 + 448*x**5 + 480*x**4 - 5056*x**3 + 1984*x**2 + 7424*x - 3008
    p = Add(*[sqrt(i) for i in range(1, 12)])
    mp = minimal_polynomial(p)(x)
    assert mp.subs({x: 0}) == -71965773323122507776


def test_minpoly_compose():
    # issue sympy/sympy#6868
    eq = (-1/(800*sqrt(Rational(-1, 240) + 1/(18000*cbrt(Rational(-1, 17280000) +
                                                         sqrt(15)*I/28800000)) + 2*cbrt(Rational(-1, 17280000) +
                                                                                        sqrt(15)*I/28800000))))
    mp = minimal_polynomial(eq + 3)(x)
    assert mp == 8000*x**2 - 48000*x + 71999

    # issue sympy/sympy#5888
    assert minimal_polynomial(exp(I*pi/8))(x) == x**8 + 1

    mp = minimal_polynomial(sin(pi/7) + sqrt(2))(x)
    assert mp == 4096*x**12 - 63488*x**10 + 351488*x**8 - 826496*x**6 + \
        770912*x**4 - 268432*x**2 + 28561
    mp = minimal_polynomial(cos(pi/7) + sqrt(2))(x)
    assert mp == 64*x**6 - 64*x**5 - 432*x**4 + 304*x**3 + 712*x**2 - \
        232*x - 239
    mp = minimal_polynomial(exp(I*pi/7) + sqrt(2))(x)
    assert mp == x**12 - 2*x**11 - 9*x**10 + 16*x**9 + 43*x**8 - 70*x**7 - 97*x**6 + 126*x**5 + 211*x**4 - 212*x**3 - 37*x**2 + 142*x + 127

    mp = minimal_polynomial(cos(pi/9))(x)
    assert mp == 8*x**3 - 6*x - 1

    mp = minimal_polynomial(sin(pi/7) + sqrt(2))(x)
    assert mp == 4096*x**12 - 63488*x**10 + 351488*x**8 - 826496*x**6 + \
        770912*x**4 - 268432*x**2 + 28561
    mp = minimal_polynomial(cos(pi/7) + sqrt(2))(x)
    assert mp == 64*x**6 - 64*x**5 - 432*x**4 + 304*x**3 + 712*x**2 - \
        232*x - 239
    mp = minimal_polynomial(exp(I*pi/7) + sqrt(2))(x)
    assert mp == x**12 - 2*x**11 - 9*x**10 + 16*x**9 + 43*x**8 - 70*x**7 - 97*x**6 + 126*x**5 + 211*x**4 - 212*x**3 - 37*x**2 + 142*x + 127

    mp = minimal_polynomial(exp(2*I*pi/7))(x)
    assert mp == x**6 + x**5 + x**4 + x**3 + x**2 + x + 1
    mp = minimal_polynomial(exp(2*I*pi/15))(x)
    assert mp == x**8 - x**7 + x**5 - x**4 + x**3 - x + 1
    mp = minimal_polynomial(cos(2*pi/7))(x)
    assert mp == 8*x**3 + 4*x**2 - 4*x - 1
    mp = minimal_polynomial(sin(2*pi/7))(x)
    ex = (5*cos(2*pi/7) - 7)/(9*cos(pi/7) - 5*cos(3*pi/7))
    mp = minimal_polynomial(ex)(x)
    assert mp == x**3 + 2*x**2 - x - 1
    assert minimal_polynomial(-1/(2*cos(pi/7)))(x) == x**3 + 2*x**2 - x - 1
    assert minimal_polynomial(sin(2*pi/15))(x) == \
        256*x**8 - 448*x**6 + 224*x**4 - 32*x**2 + 1
    assert minimal_polynomial(sin(5*pi/14))(x) == 8*x**3 - 4*x**2 - 4*x + 1
    assert minimal_polynomial(cos(pi/15))(x) == 16*x**4 + 8*x**3 - 16*x**2 - 8*x + 1
    assert minimal_polynomial(cos(pi/17))(x) == (256*x**8 - 128*x**7 -
                                                 448*x**6 + 192*x**5 +
                                                 240*x**4 - 80*x**3 -
                                                 40*x**2 + 8*x + 1)
    assert minimal_polynomial(cos(2*pi/21))(x) == (64*x**6 - 32*x**5 - 96*x**4 +
                                                   48*x**3 + 32*x**2 - 16*x + 1)

    ex = RootOf(x**3 + x*4 + 1, 0)
    mp = minimal_polynomial(ex)(x)
    assert mp == x**3 + 4*x + 1
    mp = minimal_polynomial(ex + 1)(x)
    assert mp == x**3 - 3*x**2 + 7*x - 4
    assert minimal_polynomial(exp(I*pi/3))(x) == x**2 - x + 1
    assert minimal_polynomial(exp(I*pi/4))(x) == x**4 + 1
    assert minimal_polynomial(exp(I*pi/6))(x) == x**4 - x**2 + 1
    assert minimal_polynomial(exp(I*pi/9))(x) == x**6 - x**3 + 1
    assert minimal_polynomial(exp(I*pi/10))(x) == x**8 - x**6 + x**4 - x**2 + 1
    assert minimal_polynomial(exp(I*pi/18))(x) == x**12 - x**6 + 1
    assert minimal_polynomial(sin(pi/9))(x) == 64*x**6 - 96*x**4 + 36*x**2 - 3
    assert minimal_polynomial(sin(pi/11))(x) == 1024*x**10 - 2816*x**8 + \
        2816*x**6 - 1232*x**4 + 220*x**2 - 11

    ex = cbrt(2)*exp(Rational(2, 3)*I*pi)
    assert minimal_polynomial(ex)(x) == x**3 - 2

    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(cos(pi*sqrt(2))))
    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(sin(pi*sqrt(2))))
    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(exp(I*pi*sqrt(2))))

    # issue sympy/sympy#5934
    ex = 1/(-36000 - 7200*sqrt(5) + (12*sqrt(10)*sqrt(sqrt(5) + 5) +
                                     24*sqrt(10)*sqrt(-sqrt(5) + 5))**2) + 1
    pytest.raises(ZeroDivisionError, lambda: minimal_polynomial(ex))

    ex = sqrt(1 + cbrt(2)) + sqrt(1 + root(2, 4)) + sqrt(2)
    mp = minimal_polynomial(ex)(x)
    assert degree(mp) == 48 and mp.subs({x: 0}) == -16630256576

    mp = minimal_polynomial(sin(pi/27))(x)
    assert mp == (262144*x**18 - 1179648*x**16 + 2211840*x**14 -
                  2236416*x**12 + 1317888*x**10 - 456192*x**8 +
                  88704*x**6 - 8640*x**4 + 324*x**2 - 3)


def test_minpoly_sympyissue_7113():
    # see discussion in https://github.com/sympy/sympy/pull/2234
    r = nsimplify(pi, tolerance=0.000000001)
    mp = minimal_polynomial(r)(x)
    assert mp == 1768292677839237920489538677417507171630859375*x**109 - \
        2734577732179183863586489182929671773182898498218854181690460140337930774573792597743853652058046464


def test_minpoly_sympyissue_7574():
    ex = -cbrt(-1) + (-1)**Rational(2, 3)
    assert minimal_polynomial(ex)(x) == x + 1


def test_primitive_element():
    assert primitive_element([sqrt(2)]) == (PurePoly(x**2 - 2), [1], [[1, 0]])
    assert (primitive_element([sqrt(2), sqrt(3)]) ==
            (PurePoly(x**4 - 10*x**2 + 1), [1, 1], [[QQ(+1, 2), 0, -QQ(9, 2), 0],
                                                    [QQ(-1, 2), 0, QQ(11, 2), 0]]))

    pytest.raises(ValueError, lambda: primitive_element([]))

    # issue sympy/sympy#13849
    assert (primitive_element([sqrt(2), sqrt(2) + sqrt(5)]) ==
            (PurePoly(x**4 - 76*x**2 + 4), [1, 2], [[QQ(1, 12), 0, QQ(-37, 6), 0],
                                                    [QQ(-1, 24), 0, QQ(43, 12), 0]]))

    # issue sympy/sympy#14117
    assert (primitive_element([I*sqrt(2*sqrt(2) + 3), I*sqrt(-2*sqrt(2) + 3), I]) ==
            (PurePoly(x**4 + 54*x**2 + 81), [1, 2, 4], [[QQ(1, 3), 0], [QQ(1, 27), 0, 2, 0],
                                                        [QQ(-1, 54), 0, QQ(-5, 6), 0]]))


def test_field_isomorphism():
    a = QQ.algebraic_field(sqrt(2))
    b = QQ.algebraic_field(sqrt(3))
    c = QQ.algebraic_field(sqrt(7))
    d = QQ.algebraic_field(sqrt(2) + sqrt(3))
    e = QQ.algebraic_field(sqrt(2) + sqrt(3) + sqrt(7))

    assert field_isomorphism(a, a) == [1, 0]
    assert field_isomorphism(a, b) is None
    assert field_isomorphism(a, c) is None
    assert field_isomorphism(a, d) == [QQ(1, 2), 0, -QQ(9, 2), 0]
    assert field_isomorphism(a, e) == [QQ(1, 80), 0, -QQ(1, 2), 0, QQ(59, 20), 0]

    assert field_isomorphism(b, a) is None
    assert field_isomorphism(b, b) == [1, 0]
    assert field_isomorphism(b, c) is None
    assert field_isomorphism(b, d) == [-QQ(1, 2), 0, QQ(11, 2), 0]
    assert field_isomorphism(b, e) == [-QQ(3, 640), 0, QQ(67, 320), 0, -QQ(297, 160), 0, QQ(313, 80), 0]

    assert field_isomorphism(c, a) is None
    assert field_isomorphism(c, b) is None
    assert field_isomorphism(c, c) == [1, 0]
    assert field_isomorphism(c, d) is None
    assert field_isomorphism(c, e) == [QQ(3, 640), 0, -QQ(71, 320), 0, QQ(377, 160), 0, -QQ(469, 80), 0]

    assert field_isomorphism(d, a) is None
    assert field_isomorphism(d, b) is None
    assert field_isomorphism(d, c) is None
    assert field_isomorphism(d, d) == [1, 0]
    assert field_isomorphism(d, e) == [-QQ(3, 640), 0, QQ(71, 320), 0, -QQ(377, 160), 0, QQ(549, 80), 0]

    assert field_isomorphism(e, a) is None
    assert field_isomorphism(e, b) is None
    assert field_isomorphism(e, c) is None
    assert field_isomorphism(e, d) is None
    assert field_isomorphism(e, e) == [1, 0]

    f = QQ.algebraic_field(3*sqrt(2) + 8*sqrt(7) - 5)

    assert field_isomorphism(f, e) == [QQ(3, 80), 0, -QQ(139, 80), 0, QQ(347, 20), 0, -QQ(761, 20), -5]

    assert field_isomorphism(QQ.algebraic_field(3), QQ.algebraic_field(sqrt(2))) == [1, 0]

    assert field_isomorphism(QQ.algebraic_field(+I*sqrt(3)), QQ.algebraic_field(I*sqrt(3)/2)) == [+1, 0]
    assert field_isomorphism(QQ.algebraic_field(-I*sqrt(3)), QQ.algebraic_field(I*sqrt(3)/2)) == [-1, 0]

    assert field_isomorphism(QQ.algebraic_field(+I*sqrt(3)), QQ.algebraic_field(-I*sqrt(3)/2)) == [-1, 0]
    assert field_isomorphism(QQ.algebraic_field(-I*sqrt(3)), QQ.algebraic_field(-I*sqrt(3)/2)) == [+1, 0]

    assert field_isomorphism(QQ.algebraic_field(+2*I*sqrt(3)/7), QQ.algebraic_field(5*I*sqrt(3)/3)) == [QQ(+2, 5), 0]
    assert field_isomorphism(QQ.algebraic_field(-2*I*sqrt(3)/7), QQ.algebraic_field(5*I*sqrt(3)/3)) == [QQ(-2, 5), 0]

    assert field_isomorphism(QQ.algebraic_field(+2*I*sqrt(3)/7), QQ.algebraic_field(-5*I*sqrt(3)/3)) == [QQ(-2, 5), 0]
    assert field_isomorphism(QQ.algebraic_field(-2*I*sqrt(3)/7), QQ.algebraic_field(-5*I*sqrt(3)/3)) == [QQ(+2, 5), 0]

    assert field_isomorphism(QQ.algebraic_field(+2*I*sqrt(3)/7 + 27), QQ.algebraic_field(5*I*sqrt(3)/3)) == [QQ(+2, 5), 189]
    assert field_isomorphism(QQ.algebraic_field(-2*I*sqrt(3)/7 + 27), QQ.algebraic_field(5*I*sqrt(3)/3)) == [QQ(-2, 5), 189]

    assert field_isomorphism(QQ.algebraic_field(+2*I*sqrt(3)/7 + 27), QQ.algebraic_field(-5*I*sqrt(3)/3)) == [QQ(-2, 5), 189]
    assert field_isomorphism(QQ.algebraic_field(-2*I*sqrt(3)/7 + 27), QQ.algebraic_field(-5*I*sqrt(3)/3)) == [QQ(+2, 5), 189]

    p = QQ.algebraic_field(+sqrt(2) + sqrt(3))
    q = QQ.algebraic_field(-sqrt(2) + sqrt(3))
    r = QQ.algebraic_field(+sqrt(2) - sqrt(3))
    s = QQ.algebraic_field(-sqrt(2) - sqrt(3))
    c = QQ.algebraic_field(cbrt(2))

    pos_coeffs = [+QQ(1, 2), 0, -QQ(9, 2), 0]
    neg_coeffs = [-QQ(1, 2), 0, +QQ(9, 2), 0]

    a = QQ.algebraic_field(sqrt(2))

    assert is_isomorphism_possible(a, c) is False
    assert field_isomorphism(a, c) is None

    assert is_isomorphism_possible(a, p) is True
    assert is_isomorphism_possible(a, q) is True
    assert is_isomorphism_possible(a, r) is True
    assert is_isomorphism_possible(a, s) is True

    assert field_isomorphism(a, p, fast=True) == pos_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_coeffs
    assert field_isomorphism(a, s, fast=True) == neg_coeffs

    assert field_isomorphism(a, p, fast=False) == pos_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_coeffs
    assert field_isomorphism(a, s, fast=False) == neg_coeffs

    a = QQ.algebraic_field(-sqrt(2))

    assert is_isomorphism_possible(a, p) is True
    assert is_isomorphism_possible(a, q) is True
    assert is_isomorphism_possible(a, r) is True
    assert is_isomorphism_possible(a, s) is True

    assert field_isomorphism(a, p, fast=True) == neg_coeffs
    assert field_isomorphism(a, q, fast=True) == pos_coeffs
    assert field_isomorphism(a, r, fast=True) == neg_coeffs
    assert field_isomorphism(a, s, fast=True) == pos_coeffs

    assert field_isomorphism(a, p, fast=False) == neg_coeffs
    assert field_isomorphism(a, q, fast=False) == pos_coeffs
    assert field_isomorphism(a, r, fast=False) == neg_coeffs
    assert field_isomorphism(a, s, fast=False) == pos_coeffs

    pos_coeffs = [+QQ(1, 2), 0, -QQ(11, 2), 0]
    neg_coeffs = [-QQ(1, 2), 0, +QQ(11, 2), 0]

    a = QQ.algebraic_field(sqrt(3))

    assert is_isomorphism_possible(a, p) is True
    assert is_isomorphism_possible(a, q) is True
    assert is_isomorphism_possible(a, r) is True
    assert is_isomorphism_possible(a, s) is True

    assert field_isomorphism(a, p, fast=True) == neg_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_coeffs
    assert field_isomorphism(a, s, fast=True) == pos_coeffs

    assert field_isomorphism(a, p, fast=False) == neg_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_coeffs
    assert field_isomorphism(a, s, fast=False) == pos_coeffs

    a = QQ.algebraic_field(-sqrt(3))

    assert is_isomorphism_possible(a, p) is True
    assert is_isomorphism_possible(a, q) is True
    assert is_isomorphism_possible(a, r) is True
    assert is_isomorphism_possible(a, s) is True

    assert field_isomorphism(a, p, fast=True) == pos_coeffs
    assert field_isomorphism(a, q, fast=True) == pos_coeffs
    assert field_isomorphism(a, r, fast=True) == neg_coeffs
    assert field_isomorphism(a, s, fast=True) == neg_coeffs

    assert field_isomorphism(a, p, fast=False) == pos_coeffs
    assert field_isomorphism(a, q, fast=False) == pos_coeffs
    assert field_isomorphism(a, r, fast=False) == neg_coeffs
    assert field_isomorphism(a, s, fast=False) == neg_coeffs

    pos_coeffs = [+QQ(3, 2), 0, -QQ(33, 2), -8]
    neg_coeffs = [-QQ(3, 2), 0, +QQ(33, 2), -8]

    a = QQ.algebraic_field(3*sqrt(3) - 8)

    assert is_isomorphism_possible(a, p) is True
    assert is_isomorphism_possible(a, q) is True
    assert is_isomorphism_possible(a, r) is True
    assert is_isomorphism_possible(a, s) is True

    assert field_isomorphism(a, p, fast=True) == neg_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_coeffs
    assert field_isomorphism(a, s, fast=True) == pos_coeffs

    assert field_isomorphism(a, p, fast=False) == neg_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_coeffs
    assert field_isomorphism(a, s, fast=False) == pos_coeffs

    a = QQ.algebraic_field(3*sqrt(2) + 2*sqrt(3) + 1)

    pos_1_coeffs = [+QQ(1, 2), 0, -QQ(5, 2), 1]
    neg_1_coeffs = [-QQ(1, 2), 0, +QQ(5, 2), 1]
    pos_5_coeffs = [+QQ(5, 2), 0, -QQ(49, 2), 1]
    neg_5_coeffs = [-QQ(5, 2), 0, +QQ(49, 2), 1]

    assert is_isomorphism_possible(a, p) is True
    assert is_isomorphism_possible(a, q) is True
    assert is_isomorphism_possible(a, r) is True
    assert is_isomorphism_possible(a, s) is True

    assert field_isomorphism(a, p, fast=True) == pos_1_coeffs
    assert field_isomorphism(a, q, fast=True) == neg_5_coeffs
    assert field_isomorphism(a, r, fast=True) == pos_5_coeffs
    assert field_isomorphism(a, s, fast=True) == neg_1_coeffs

    assert field_isomorphism(a, p, fast=False) == pos_1_coeffs
    assert field_isomorphism(a, q, fast=False) == neg_5_coeffs
    assert field_isomorphism(a, r, fast=False) == pos_5_coeffs
    assert field_isomorphism(a, s, fast=False) == neg_1_coeffs

    a = QQ.algebraic_field(sqrt(2))
    b = QQ.algebraic_field(sqrt(3))
    c = QQ.algebraic_field(sqrt(7))

    assert is_isomorphism_possible(a, b) is True
    assert is_isomorphism_possible(b, a) is True

    assert is_isomorphism_possible(c, p) is False

    assert field_isomorphism(a, b, fast=True) is None
    assert field_isomorphism(b, a, fast=True) is None

    assert field_isomorphism(a, b, fast=False) is None
    assert field_isomorphism(b, a, fast=False) is None

    pytest.raises(ValueError, lambda: field_isomorphism(1, 2))


def test_to_number_field():
    A = QQ.algebraic_field(sqrt(2))
    assert A.convert(sqrt(2)) == A([1, 0])
    B = A.algebraic_field(sqrt(3))
    assert B.convert(sqrt(2) + sqrt(3)) == B([1, 0])

    a = AlgebraicNumber(sqrt(2) + sqrt(3), [Rational(1, 2), Integer(0), -Rational(9, 2), Integer(0)])

    assert B.convert(sqrt(2)) == B(a.coeffs())

    pytest.raises(CoercionFailed, lambda: QQ.algebraic_field(sqrt(3)).convert(sqrt(2)))

    # issue sympy/sympy#5649
    assert AlgebraicNumber(1).rep.rep == QQ.algebraic_field(1).convert(1).rep
    assert AlgebraicNumber(sqrt(2)).rep.rep == A.convert(sqrt(2)).rep

    p = x**6 - 6*x**4 - 6*x**3 + 12*x**2 - 36*x + 1
    r0, r1 = p.as_poly(x).all_roots()[:2]
    a = AlgebraicNumber(r0, [Rational(-96, 755), Rational(-54, 755),
                             Rational(128, 151), Rational(936, 755),
                             Rational(-1003, 755), Rational(2184, 755)])
    A = QQ.algebraic_field(r0)
    assert A.convert(r1) == A(a.coeffs())


def test_AlgebraicNumber():
    minpoly, root = x**2 - 2, sqrt(2)

    a = AlgebraicNumber(root, gen=x)

    assert a.rep == DMP([QQ(1), QQ(0)], QQ)
    assert a.root == root
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.coeffs() == [Integer(1), Integer(0)]
    assert a.native_coeffs() == [QQ(1), QQ(0)]

    a = AlgebraicNumber(root, gen=x)

    assert a.rep == DMP([QQ(1), QQ(0)], QQ)
    assert a.root == root
    assert a.minpoly == minpoly
    assert a.is_number

    a = AlgebraicNumber(root, gen=x)

    assert a.rep == DMP([QQ(1), QQ(0)], QQ)
    assert a.root == root
    assert a.minpoly == minpoly
    assert a.is_number

    assert AlgebraicNumber(sqrt(2), []).rep == DMP([], QQ)

    assert AlgebraicNumber(sqrt(2), [8]).rep == DMP([QQ(8)], QQ)
    assert AlgebraicNumber(sqrt(2), [Rational(8, 3)]).rep == DMP([QQ(8, 3)], QQ)

    assert AlgebraicNumber(sqrt(2), [7, 3]).rep == DMP([QQ(7), QQ(3)], QQ)
    assert AlgebraicNumber(
        sqrt(2), [Rational(7, 9), Rational(3, 2)]).rep == DMP([QQ(7, 9), QQ(3, 2)], QQ)

    assert AlgebraicNumber(sqrt(2), [1, 2, 3]).rep == DMP([QQ(2), QQ(5)], QQ)

    a = AlgebraicNumber(AlgebraicNumber(root, gen=x), [1, 2])

    assert a.rep == DMP([QQ(1), QQ(2)], QQ)
    assert a.root == root
    assert a.minpoly == minpoly
    assert a.is_number

    assert a.coeffs() == [Integer(1), Integer(2)]
    assert a.native_coeffs() == [QQ(1), QQ(2)]

    a = AlgebraicNumber((minpoly, root), [1, 2])

    assert a.rep == DMP([QQ(1), QQ(2)], QQ)
    assert a.root == root
    assert a.minpoly == minpoly
    assert a.is_number

    a = AlgebraicNumber((Poly(minpoly), root), [1, 2])

    assert a.rep == DMP([QQ(1), QQ(2)], QQ)
    assert a.root == root
    assert a.minpoly == minpoly
    assert a.is_number

    assert AlgebraicNumber( sqrt(3)).rep == DMP([ QQ(1), QQ(0)], QQ)
    assert AlgebraicNumber(-sqrt(3)).rep == DMP([ QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(2))

    assert a == b

    c = AlgebraicNumber(sqrt(2), gen=x)

    assert a == b
    assert a == c

    a = AlgebraicNumber(sqrt(2), [1, 2])
    b = AlgebraicNumber(sqrt(2), [1, 3])

    assert a != b and a != sqrt(2) + 3

    assert (a == x) is False and (a != x) is True

    a = AlgebraicNumber(sqrt(2), [1, 0])

    assert a.as_poly(x) == Poly(x)

    assert a.as_expr() == sqrt(2)
    assert a.as_expr(x) == x

    a = AlgebraicNumber(sqrt(2), [2, 3])

    p = a.as_poly()

    assert p == Poly(2*p.gen + 3)

    assert a.as_poly(x) == Poly(2*x + 3)

    assert a.as_expr() == 2*sqrt(2) + 3
    assert a.as_expr(x) == 2*x + 3

    A = QQ.algebraic_field(AlgebraicNumber(sqrt(2)))
    a = A([1, 0])
    b = A.convert(sqrt(2))
    assert a == b

    a = AlgebraicNumber(sqrt(2), [1, 2, 3])
    assert a.args == (sqrt(2), Tuple(2, 5))

    pytest.raises(ValueError, lambda: AlgebraicNumber(RootOf(x**3 + y*x + 1,
                                                             x, 0)))

    a = AlgebraicNumber(RootOf(x**3 + 2*x - 1, 1))
    assert a.free_symbols == set()

    # integer powers:
    assert a**0 == 1
    assert a**2 == AlgebraicNumber(a, (1, 0, 0))
    assert a**5 == AlgebraicNumber(a, (1, 0, 0, 0, 0, 0))
    assert a**110 == AlgebraicNumber(a, ([1] + [0]*110))
    assert (a**pi).is_Pow

    b = AlgebraicNumber(sqrt(3))
    assert b + 1 == AlgebraicNumber(sqrt(3), (1, 1))
    assert (b + 1) + b == AlgebraicNumber(sqrt(3), (2, 1))

    assert (2*b + 1)**3 == 30*b + 37
    assert 1/b == b/3

    b = AlgebraicNumber(RootOf(x**7 - x + 1, 1), (1, 2, -1))
    t = AlgebraicNumber(b.root)
    assert b**7 == 490*t**6 - 119*t**5 - 196*t**4 - 203*t**3 - 265*t**2 + 637*t - 198

    b = AlgebraicNumber(sqrt(2), (1, 0))
    c = b + 1
    assert c**2 == 2*b + 3
    assert c**5 == 29*b + 41
    assert c**-2 == 3 - 2*b
    assert c**-11 == 5741*b - 8119

    # arithmetics
    assert a**3 == -2*a + 1 == a*(-2) + 1 == 1 + (-2)*a == 1 - 2*a
    assert a**5 == a**2 + 4*a - 2
    assert a**4 == -2*a**2 + a == a - 2*a**2
    assert a**110 == (-2489094528619081871*a**2 + 3737645722703173544*a -
                      1182958048412500088)

    assert a + a == 2*a
    assert 2*a - a == a
    assert Integer(1) - a == (-a) + 1

    assert (a + pi).is_Add
    assert (pi + a).is_Add
    assert (a - pi).is_Add
    assert (pi - a).is_Add
    assert (a*pi).is_Mul
    assert (pi*a).is_Mul

    a = AlgebraicNumber(sqrt(2), (5, 7))
    b = AlgebraicNumber(sqrt(2), (-4, -6))
    assert a*b == AlgebraicNumber(sqrt(2), (-58, -82))

    # different extensions
    a = AlgebraicNumber(sqrt(2))
    b = AlgebraicNumber(sqrt(3))
    assert a + b == Add(a, +b, evaluate=False)
    assert a * b == Mul(a, +b, evaluate=False)
    assert a - b == Add(a, -b, evaluate=False)


def test_to_algebraic_integer():
    a = AlgebraicNumber(sqrt(3), gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 3
    assert a.root == sqrt(3)
    assert a.rep == DMP([QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(2*sqrt(3), gen=x).to_algebraic_integer()
    assert a.minpoly == x**2 - 12
    assert a.root == 2*sqrt(3)
    assert a.rep == DMP([QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(3)/2, gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 12
    assert a.root == 2*sqrt(3)
    assert a.rep == DMP([QQ(1), QQ(0)], QQ)

    a = AlgebraicNumber(sqrt(3)/2, [Rational(7, 19), 3], gen=x).to_algebraic_integer()

    assert a.minpoly == x**2 - 12
    assert a.root == 2*sqrt(3)
    assert a.rep == DMP([QQ(7, 19), QQ(3)], QQ)


def test_minpoly_fraction_field():
    assert minimal_polynomial(1/x)(y) == x*y - 1
    assert minimal_polynomial(1/(x + 1))(y) == x*y + y - 1

    assert minimal_polynomial(sqrt(x))(y) == y**2 - x
    assert minimal_polynomial(sqrt(x), method='groebner')(y) == y**2 - x

    assert minimal_polynomial(sqrt(x + 1))(y) == y**2 - x - 1
    assert minimal_polynomial(sqrt(x)/x)(y) == x*y**2 - 1
    assert minimal_polynomial(sqrt(2)*sqrt(x))(y) == y**2 - 2 * x

    assert minimal_polynomial(sqrt(2) + sqrt(x))(y) == \
        y**4 - 2*x*y**2 - 4*y**2 + x**2 - 4*x + 4
    assert minimal_polynomial(sqrt(2) + sqrt(x), method='groebner')(y) == \
        y**4 - 2*x*y**2 - 4*y**2 + x**2 - 4*x + 4

    assert minimal_polynomial(cbrt(x))(y) == y**3 - x
    assert minimal_polynomial(cbrt(x) + sqrt(x))(y) == \
        y**6 - 3*x*y**4 - 2*x*y**3 + 3*x**2*y**2 - 6*x**2*y - x**3 + x**2

    assert minimal_polynomial(sqrt(x)/z)(y) == z**2*y**2 - x
    assert minimal_polynomial(sqrt(x)/(z + 1))(y) == z**2*y**2 + 2*z*y**2 + y**2 - x

    assert minimal_polynomial(1/x) == PurePoly(x*y - 1, y)
    assert minimal_polynomial(1/(x + 1)) == PurePoly((x + 1)*y - 1, y)
    assert minimal_polynomial(sqrt(x)) == PurePoly(y**2 - x, y)
    assert minimal_polynomial(sqrt(x) / z) == PurePoly(z**2*y**2 - x, y)

    # this is (sqrt(1 + x**3)/x).integrate(x).diff(x) - sqrt(1 + x**3)/x
    a = sqrt(x)/sqrt(1 + x**(-3)) - sqrt(x**3 + 1)/x + 1/(x**Rational(5, 2) *
                                                          (1 + x**(-3))**Rational(3, 2)) + 1/(x**Rational(11, 2)*(1 + x**(-3))**Rational(3, 2))

    assert minimal_polynomial(a)(y) == y

    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(exp(x)))


@pytest.mark.slow
def test_minpoly_fraction_field_slow():
    assert minimal_polynomial(minimal_polynomial(sqrt(root(x, 5) - 1))(y).subs(y, sqrt(root(x, 5) - 1)))(z) == z


def test_minpoly_domain():
    F = QQ.algebraic_field(sqrt(2))

    assert minimal_polynomial(sqrt(2), domain=F)(x) == x - sqrt(2)
    assert minimal_polynomial(sqrt(8), domain=F)(x) == x - 2*sqrt(2)
    assert minimal_polynomial(sqrt(Rational(3, 2)), domain=F)(x) == 2*x**2 - 3

    pytest.raises(NotAlgebraic, lambda: minimal_polynomial(y, domain=QQ))

    # issue sympy/sympy#14494

    F = QQ.algebraic_field(I)
    assert minimal_polynomial(I, domain=F)(x) == x - I

    F = QQ.algebraic_field(sqrt(3)*I)
    assert minimal_polynomial(exp(I*pi/3), domain=F)(x) == 2*x - sqrt(3)*I - 1


def test_sympyissue_11553():
    eqs = (x + y + 1, x + GoldenRatio)
    assert solve(eqs, x, y) == [{x: -GoldenRatio, y: -1 + GoldenRatio}]
