"""Tests for algorithms for computing symbolic roots of polynomials."""

import itertools

import mpmath
import pytest

from diofant import (EX, QQ, ZZ, I, Integer, Interval, Piecewise,
                     PolynomialError, Rational, RootOf, Symbol, Wild, acos,
                     cbrt, cos, cyclotomic_poly, exp, im, legendre_poly,
                     nroots, pi, powsimp, re, root, roots, sin, sqrt, symbols)
from diofant.abc import a, b, c, d, e, q, x, y, z
from diofant.polys.polyroots import (preprocess_roots, root_factors,
                                     roots_binomial, roots_cubic,
                                     roots_cyclotomic, roots_linear,
                                     roots_quadratic, roots_quartic,
                                     roots_quintic)
from diofant.polys.polyutils import _nsort
from diofant.utilities.randtest import verify_numerically


__all__ = ()


def test_roots_linear():
    assert roots_linear((2*x + 1).as_poly()) == [-Rational(1, 2)]


def test_roots_quadratic():
    assert roots_quadratic((2*x**2).as_poly()) == [0, 0]
    assert roots_quadratic((2*x**2 + 3*x).as_poly()) == [-Rational(3, 2), 0]
    assert roots_quadratic((2*x**2 + 3).as_poly()) == [-I*sqrt(6)/2, I*sqrt(6)/2]
    assert roots_quadratic((2*x**2 + 4*x + 3).as_poly()) == [-1 - I*sqrt(2)/2, -1 + I*sqrt(2)/2]

    f = x**2 + (2*a*e + 2*c*e)/(a - c)*x + (d - b + a*e**2 - c*e**2)/(a - c)
    assert (roots_quadratic(f.as_poly(x)) ==
            [-e*(a + c)/(a - c) - sqrt((a*b + 4*a*c*e**2 -
                                        a*d - b*c + c*d)/(a - c)**2),
             -e*(a + c)/(a - c) + sqrt((a*b + 4*a*c*e**2 -
                                        a*d - b*c + c*d)/(a - c)**2)])

    # check for simplification
    f = (y*x**2 - 2*x - 2*y).as_poly(x)
    assert roots_quadratic(f) == [-sqrt((2*y**2 + 1)/y**2) + 1/y,
                                  sqrt((2*y**2 + 1)/y**2) + 1/y]
    f = (x**2 + (-y**2 - 2)*x + y**2 + 1).as_poly(x)
    assert roots_quadratic(f) == [y**2/2 - sqrt(y**4)/2 + 1,
                                  y**2/2 + sqrt(y**4)/2 + 1]

    f = (sqrt(2)*x**2 - 1).as_poly(x)
    r = roots_quadratic(f)
    assert r == _nsort(r)

    # issue sympy/sympy#8255
    f = (-24*x**2 - 180*x + 264).as_poly()
    assert [w.evalf(2) for w in f.all_roots(radicals=True)] == \
           [w.evalf(2) for w in f.all_roots(radicals=False)]
    for _a, _b, _c in itertools.product((-2, 2), (-2, 2), (0, -1)):
        f = (_a*x**2 + _b*x + _c).as_poly()
        roots = roots_quadratic(f)
        assert roots == _nsort(roots)


def test_sympyissue_8285():
    roots = ((4*x**8 - 1).as_poly()*(x**2 + 1).as_poly()).all_roots()
    assert roots == _nsort(roots)
    f = (x**4 + 5*x**2 + 6).as_poly()
    ro = [RootOf(f, i) for i in range(4)]
    roots = (x**4 + 5*x**2 + 6).as_poly().all_roots()
    assert roots == ro
    assert roots == _nsort(roots)
    # more than 2 complex roots from which to identify the
    # imaginary ones
    roots = (2*x**8 - 1).as_poly().all_roots()
    assert roots == _nsort(roots)
    assert len((2*x**10 - 1).as_poly().all_roots()) == 10  # doesn't fail


def test_sympyissue_8289():
    roots = ((x**2 + 2).as_poly()*(x**4 + 2).as_poly()).all_roots()
    assert roots == _nsort(roots)
    roots = (x**6 + 3*x**3 + 2).as_poly().all_roots()
    assert roots == _nsort(roots)
    roots = (x**6 - x + 1).as_poly().all_roots()
    assert roots == _nsort(roots)
    # all imaginary roots
    roots = (x**4 + 4*x**2 + 4).as_poly().all_roots()
    assert roots == _nsort(roots)


def test_sympyissue_14293():
    roots = (x**8 + 2*x**6 + 37*x**4 - 36*x**2 + 324).as_poly().all_roots()
    assert roots == _nsort(roots)


def test_roots_cubic():
    assert roots_cubic((2*x**3).as_poly()) == [0, 0, 0]
    assert roots_cubic((x**3 - 3*x**2 + 3*x - 1).as_poly()) == [1, 1, 1]

    assert roots_cubic((x**3 + 1).as_poly()) == \
        [-1, Rational(1, 2) - I*sqrt(3)/2, Rational(1, 2) + I*sqrt(3)/2]
    assert roots_cubic((2*x**3 - 3*x**2 - 3*x - 1).as_poly())[0] == \
        Rational(1, 2) + cbrt(3)/2 + 3**Rational(2, 3)/2
    eq = -x**3 + 2*x**2 + 3*x - 2
    assert roots(eq, trig=True, multiple=True) == \
        roots_cubic(eq.as_poly(), trig=True) == [
        Rational(2, 3) + 2*sqrt(13)*cos(acos(8*sqrt(13)/169)/3)/3,
        -2*sqrt(13)*sin(-acos(8*sqrt(13)/169)/3 + pi/6)/3 + Rational(2, 3),
        -2*sqrt(13)*cos(-acos(8*sqrt(13)/169)/3 + pi/3)/3 + Rational(2, 3),
    ]
    res = roots_cubic((x**3 + 2*a/27).as_poly(x))
    assert res == [-root(2, 3)*root(a, 3)/3,
                   -root(2, 3)*root(a, 3)*(-Rational(1, 2) + sqrt(3)*I/2)/3,
                   -root(2, 3)*root(a, 3)*(-Rational(1, 2) - sqrt(3)*I/2)/3]
    res = roots_cubic((x**3 - 2*a/27).as_poly(x))
    assert res == [root(2, 3)*root(a, 3)/3,
                   root(2, 3)*root(a, 3)*(-Rational(1, 2) + sqrt(3)*I/2)/3,
                   root(2, 3)*root(a, 3)*(-Rational(1, 2) - sqrt(3)*I/2)/3]

    # issue sympy/sympy#8438
    p = -3*x**3 - 2*x**2 + x*y + 1
    croots = roots_cubic(p.as_poly(x), x)
    z = -Rational(3, 2) - 7*I/2  # this will fail in code given in commit msg
    post = [r.subs({y: z}) for r in croots]
    assert set(post) == set(roots_cubic(p.subs({y: z}).as_poly(x)))
    # /!\ if p is not made an expression, this is *very* slow
    assert all(p.subs({y: z, x: i}).evalf(2, chop=True) == 0 for i in post)


def test_roots_quartic():
    assert roots_quartic((x**4).as_poly()) == [0, 0, 0, 0]
    assert roots_quartic((x**4 + x**3).as_poly()) in [
        [-1, 0, 0, 0],
        [0, -1, 0, 0],
        [0, 0, -1, 0],
        [0, 0, 0, -1]
    ]
    assert roots_quartic((x**4 - x**3).as_poly()) in [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ]

    lhs = roots_quartic((x**4 + x).as_poly())
    rhs = [Rational(1, 2) + I*sqrt(3)/2, Rational(1, 2) - I*sqrt(3)/2, 0, -1]

    assert sorted(lhs, key=hash) == sorted(rhs, key=hash)

    # test of all branches of roots quartic
    for i, (a, b, c, d) in enumerate([(1, 2, 3, 0),
                                      (3, -7, -9, 9),
                                      (1, 2, 3, 4),
                                      (1, 2, 3, 4),
                                      (-7, -3, 3, -6),
                                      (-3, 5, -6, -4),
                                      (6, -5, -10, -3)]):
        if i == 2:
            c = -a*(a**2/Integer(8) - b/Integer(2))
        elif i == 3:
            d = a*(a*(3*a**2/Integer(256) - b/Integer(16)) + c/Integer(4))
        eq = x**4 + a*x**3 + b*x**2 + c*x + d
        ans = roots_quartic(eq.as_poly())
        assert all(eq.subs({x: ai}).evalf(chop=True) == 0 for ai in ans)

    # not all symbolic quartics are unresolvable
    eq = (q*x + q/4 + x**4 + x**3 + 2*x**2 - Rational(1, 3)).as_poly(x)
    sol = roots_quartic(eq)
    assert all(verify_numerically(eq.subs({x: i}), 0) for i in sol)
    z = symbols('z', negative=True)
    eq = x**4 + 2*x**3 + 3*x**2 + x*(z + 11) + 5
    zans = roots_quartic(eq.as_poly(x))
    assert all(verify_numerically(eq.subs({x: i, z: -1}), 0) for i in zans)
    # but some are (see also issue sympy/sympy#4989)
    # it's ok if the solution is not Piecewise, but the tests below should pass
    eq = (y*x**4 + x**3 - x + z).as_poly(x)
    ans = roots_quartic(eq)
    assert all(type(i) == Piecewise for i in ans)
    reps = ({y: -Rational(1, 3), z: -Rational(1, 4)},  # 4 real
            {y: -Rational(1, 3), z: -Rational(1, 2)},  # 2 real
            {y: -Rational(1, 3), z: -2})  # 0 real
    for rep in reps:
        sol = roots_quartic(eq.subs(rep).as_poly(x))
        assert all(verify_numerically(w.subs(rep) - s, 0) for w, s in zip(ans, sol))


def test_roots_quintic():
    assert roots_quintic((x**5 + x**4 + 1).as_poly()) == []
    assert roots_quintic((x**5 - 6*x + 3).as_poly()) == []
    assert roots_quintic((6*x**5 + 9*x**3 - 10*x**2 - 9*x).as_poly()) == []


def test_roots_cyclotomic():
    assert roots_cyclotomic(cyclotomic_poly(1, x, polys=True)) == [1]
    assert roots_cyclotomic(cyclotomic_poly(2, x, polys=True)) == [-1]
    assert roots_cyclotomic(cyclotomic_poly(
        3, x, polys=True)) == [-Rational(1, 2) - I*sqrt(3)/2, -Rational(1, 2) + I*sqrt(3)/2]
    assert roots_cyclotomic(cyclotomic_poly(4, x, polys=True)) == [-I, I]
    assert roots_cyclotomic(cyclotomic_poly(
        6, x, polys=True)) == [Rational(1, 2) - I*sqrt(3)/2, Rational(1, 2) + I*sqrt(3)/2]

    assert roots_cyclotomic(cyclotomic_poly(7, x, polys=True)) == [
        -cos(pi/7) - I*sin(pi/7),
        -cos(pi/7) + I*sin(pi/7),
        -cos(3*pi/7) - I*sin(3*pi/7),
        -cos(3*pi/7) + I*sin(3*pi/7),
        cos(2*pi/7) - I*sin(2*pi/7),
        cos(2*pi/7) + I*sin(2*pi/7),
    ]

    assert roots_cyclotomic(cyclotomic_poly(8, x, polys=True)) == [
        -sqrt(2)/2 - I*sqrt(2)/2,
        -sqrt(2)/2 + I*sqrt(2)/2,
        sqrt(2)/2 - I*sqrt(2)/2,
        sqrt(2)/2 + I*sqrt(2)/2,
    ]

    assert roots_cyclotomic(cyclotomic_poly(12, x, polys=True)) == [
        -sqrt(3)/2 - I/2,
        -sqrt(3)/2 + I/2,
        sqrt(3)/2 - I/2,
        sqrt(3)/2 + I/2,
    ]

    assert roots_cyclotomic(
        cyclotomic_poly(1, x, polys=True), factor=True) == [1]
    assert roots_cyclotomic(
        cyclotomic_poly(2, x, polys=True), factor=True) == [-1]

    assert roots_cyclotomic(cyclotomic_poly(3, x, polys=True), factor=True) == \
        [-root(-1, 3), -1 + root(-1, 3)]
    assert roots_cyclotomic(cyclotomic_poly(4, x, polys=True), factor=True) == \
        [-I, I]
    assert roots_cyclotomic(cyclotomic_poly(5, x, polys=True), factor=True) == \
        [-root(-1, 5), -root(-1, 5)**3, root(-1, 5)**2, -1 - root(-1, 5)**2 + root(-1, 5) + root(-1, 5)**3]

    assert roots_cyclotomic(cyclotomic_poly(6, x, polys=True), factor=True) == \
        [1 - root(-1, 3), root(-1, 3)]


def test_roots_binomial():
    assert roots_binomial((5*x).as_poly()) == [0]
    assert roots_binomial((5*x**4).as_poly()) == [0, 0, 0, 0]
    assert roots_binomial((5*x + 2).as_poly()) == [-Rational(2, 5)]

    A = 10**Rational(3, 4)/10

    assert roots_binomial((5*x**4 + 2).as_poly()) == \
        [-A - A*I, -A + A*I, A - A*I, A + A*I]

    a1 = Symbol('a1', nonnegative=True)
    b1 = Symbol('b1', nonnegative=True)

    r0 = roots_quadratic((a1*x**2 + b1).as_poly(x))
    r1 = roots_binomial((a1*x**2 + b1).as_poly(x))

    assert powsimp(r0[0]) == powsimp(r1[0])
    assert powsimp(r0[1]) == powsimp(r1[1])
    for a, b, s, n in itertools.product((1, 2), (1, 2), (-1, 1), (2, 3, 4, 5)):
        if a == b and a != 1:  # a == b == 1 is sufficient
            continue
        p = (a*x**n + s*b).as_poly()
        ans = roots_binomial(p)
        assert ans == _nsort(ans)

    # issue sympy/sympy#8813
    assert roots((2*x**3 - 16*y**3).as_poly(x)) == {
        2*y*(-Rational(1, 2) - sqrt(3)*I/2): 1,
        2*y: 1,
        2*y*(-Rational(1, 2) + sqrt(3)*I/2): 1}

    p = (exp(I*x/3)**4 + exp(I*x/3)).as_poly(exp(I*x/3))
    assert roots(p) == roots(x**4 + x)


def test_roots_preprocessing():
    f = a*y*x**2 + y - b

    coeff, poly = preprocess_roots(f.as_poly(x))

    assert coeff == 1
    assert poly == (a*y*x**2 + y - b).as_poly(x)

    f = c**3*x**3 + c**2*x**2 + c*x + a

    coeff, poly = preprocess_roots(f.as_poly(x))

    assert coeff == 1/c
    assert poly == (x**3 + x**2 + x + a).as_poly(x)

    f = c**3*x**3 + c**2*x**2 + a

    coeff, poly = preprocess_roots(f.as_poly(x))

    assert coeff == 1/c
    assert poly == (x**3 + x**2 + a).as_poly(x)

    f = c**3*x**3 + c*x + a

    coeff, poly = preprocess_roots(f.as_poly(x))

    assert coeff == 1/c
    assert poly == (x**3 + x + a).as_poly(x)

    f = c**3*x**3 + a

    coeff, poly = preprocess_roots(f.as_poly(x))

    assert coeff == 1/c
    assert poly == (x**3 + a).as_poly(x)

    E, F, J, L = symbols('E,F,J,L')

    f = -21601054687500000000*E**8*J**8/L**16 + \
        508232812500000000*F*x*E**7*J**7/L**14 - \
        4269543750000000*E**6*F**2*J**6*x**2/L**12 + \
        16194716250000*E**5*F**3*J**5*x**3/L**10 - \
        27633173750*E**4*F**4*J**4*x**4/L**8 + \
        14840215*E**3*F**5*J**3*x**5/L**6 + \
        54794*E**2*F**6*J**2*x**6/(5*L**4) - \
        1153*E*J*F**7*x**7/(80*L**2) + \
        633*F**8*x**8/160000

    coeff, poly = preprocess_roots(f.as_poly(x))

    assert coeff == 20*E*J/(F*L**2)
    assert poly == 633*x**8 - 115300*x**7 + 4383520*x**6 + 296804300*x**5 - 27633173750*x**4 + \
        809735812500*x**3 - 10673859375000*x**2 + 63529101562500*x - 135006591796875

    f = J**8 + 7*E*x**2*L**16 + 5*F*x*E**5*J**7*L**2
    coeff, poly = preprocess_roots(f.as_poly(x))
    assert coeff == 1
    assert poly == f.as_poly(x)

    f = (-y**2 + x**2*exp(x)).as_poly(y, domain=ZZ.inject(x, exp(x)))
    g = (y**2 - exp(x)).as_poly(y, domain=ZZ.inject(exp(x)))

    assert preprocess_roots(f) == (x, g)


def test_roots0():
    assert roots(1, x) == {}
    assert roots(x, x) == {0: 1}
    assert roots(x**9, x) == {0: 9}
    assert roots(((x - 2)*(x + 3)*(x - 4)).expand(), x) == {-3: 1, 2: 1, 4: 1}

    assert roots(x**2 - 2*x + 1, x, auto=False) == {1: 2}

    assert roots(2*x + 1, x) == {Rational(-1, 2): 1}
    assert roots((2*x + 1)**2, x) == {Rational(-1, 2): 2}
    assert roots((2*x + 1)**5, x) == {Rational(-1, 2): 5}
    assert roots((2*x + 1)**10, x) == {Rational(-1, 2): 10}

    assert roots(x**4 - 1, x) == {I: 1, 1: 1, -1: 1, -I: 1}
    assert roots((x**4 - 1)**2, x) == {I: 2, 1: 2, -1: 2, -I: 2}

    assert roots(((2*x - 3)**2).expand(), x) == {+Rational(3, 2): 2}
    assert roots(((2*x + 3)**2).expand(), x) == {-Rational(3, 2): 2}

    assert roots(((2*x - 3)**3).expand(), x) == {+Rational(3, 2): 3}
    assert roots(((2*x + 3)**3).expand(), x) == {-Rational(3, 2): 3}

    assert roots(((2*x - 3)**5).expand(), x) == {+Rational(3, 2): 5}
    assert roots(((2*x + 3)**5).expand(), x) == {-Rational(3, 2): 5}

    assert roots(((a*x - b)**5).expand(), x) == {+b/a: 5}
    assert roots(((a*x + b)**5).expand(), x) == {-b/a: 5}

    assert roots(x**2 + (-a - 1)*x + a, x) == {a: 1, 1: 1}

    assert roots(x**4 - 2*x**2 + 1, x) == {1: 2, -1: 2}

    assert roots(x**6 - 4*x**4 + 4*x**3 - x**2, x) == \
        {1: 2, -1 - sqrt(2): 1, 0: 2, -1 + sqrt(2): 1}

    assert roots(x**8 - 1, x) == {
        sqrt(2)/2 + I*sqrt(2)/2: 1,
        sqrt(2)/2 - I*sqrt(2)/2: 1,
        -sqrt(2)/2 + I*sqrt(2)/2: 1,
        -sqrt(2)/2 - I*sqrt(2)/2: 1,
        1: 1, -1: 1, I: 1, -I: 1
    }

    f = -2016*x**2 - 5616*x**3 - 2056*x**4 + 3324*x**5 + 2176*x**6 - \
        224*x**7 - 384*x**8 - 64*x**9

    assert roots(f) == {0: 2, -2: 2, 2: 1, -Rational(7, 2): 1, -Rational(3, 2): 1, -Rational(1, 2): 1, Rational(3, 2): 1}

    assert roots((a + b + c)*x - (a + b + c + d), x) == {(a + b + c + d)/(a + b + c): 1}

    assert roots(x**3 + x**2 - x + 1, x, cubics=False) == {}
    assert roots(((x - 2)*(
        x + 3)*(x - 4)).expand(), x, cubics=False) == {-3: 1, 2: 1, 4: 1}
    assert roots(((x - 2)*(x + 3)*(x - 4)*(x - 5)).expand(), x, cubics=False) == \
        {-3: 1, 2: 1, 4: 1, 5: 1}
    assert roots(x**3 + 2*x**2 + 4*x + 8, x) == {-2: 1, -2*I: 1, 2*I: 1}
    assert roots(x**3 + 2*x**2 + 4*x + 8, x, cubics=True) == \
        {-2*I: 1, 2*I: 1, -2: 1}
    assert roots((x**2 - x)*(x**3 + 2*x**2 + 4*x + 8), x) == \
        {1: 1, 0: 1, -2: 1, -2*I: 1, 2*I: 1}

    r1_2, r1_3 = Rational(1, 2), Rational(1, 3)

    x0 = (3*sqrt(33) + 19)**r1_3
    x1 = 4/x0/3
    x2 = x0/3
    x3 = sqrt(3)*I/2
    x4 = x3 - r1_2
    x5 = -x3 - r1_2
    assert roots(x**3 + x**2 - x + 1, x, cubics=True) == {
        -x1 - x2 - r1_3: 1,
        -x1/x4 - x2*x4 - r1_3: 1,
        -x1/x5 - x2*x5 - r1_3: 1,
    }

    f = (x**2 + 2*x + 3).subs({x: 2*x**2 + 3*x}).subs({x: 5*x - 4})

    r13_20, r1_20 = [Rational(*r) for r in ((13, 20), (1, 20))]

    s2 = sqrt(2)
    assert roots(f, x) == {
        r13_20 + r1_20*sqrt(1 - 8*I*s2): 1,
        r13_20 - r1_20*sqrt(1 - 8*I*s2): 1,
        r13_20 + r1_20*sqrt(1 + 8*I*s2): 1,
        r13_20 - r1_20*sqrt(1 + 8*I*s2): 1,
    }

    f = x**4 + x**3 + x**2 + x + 1

    r1_4, r1_8, r5_8 = [Rational(*r) for r in ((1, 4), (1, 8), (5, 8))]

    assert roots(f, x) == {
        -r1_4 + r1_4*5**r1_2 + I*(r5_8 + r1_8*5**r1_2)**r1_2: 1,
        -r1_4 + r1_4*5**r1_2 - I*(r5_8 + r1_8*5**r1_2)**r1_2: 1,
        -r1_4 - r1_4*5**r1_2 + I*(r5_8 - r1_8*5**r1_2)**r1_2: 1,
        -r1_4 - r1_4*5**r1_2 - I*(r5_8 - r1_8*5**r1_2)**r1_2: 1,
    }

    f = z**3 + (-2 - y)*z**2 + (1 + 2*y - 2*x**2)*z - y + 2*x**2

    assert roots(f, z) == {
        1: 1,
        Rational(1, 2) + y/2 + sqrt(1 - 2*y + y**2 + 8*x**2)/2: 1,
        Rational(1, 2) + y/2 - sqrt(1 - 2*y + y**2 + 8*x**2)/2: 1,
    }

    assert roots(a*b*c*x**3 + 2*x**2 + 4*x + 8, x, cubics=False) == {}
    assert roots(a*b*c*x**3 + 2*x**2 + 4*x + 8, x, cubics=True) != {}

    assert roots(x**4 - 1, x, filter='Z') == {1: 1, -1: 1}
    assert roots(x**4 - 1, x, filter='R') == {1: 1, -1: 1}
    assert roots(x**4 - 1, x, filter='I') == {I: 1, -I: 1}

    pytest.raises(ValueError, lambda: roots(x**4 - 1, x, filter='spam'))

    assert roots((x - 1)*(x + 1), x) == {1: 1, -1: 1}
    assert roots(
        (x - 1)*(x + 1), x, predicate=lambda r: r.is_positive) == {1: 1}

    assert roots(x**4 - 1, x, filter='Z', multiple=True) == [-1, 1]
    assert roots(x**4 - 1, x, filter='I', multiple=True) == [I, -I]

    assert roots(x**3, x, multiple=True) == [0, 0, 0]
    assert roots(1234, x, multiple=True) == []

    f = x**6 - x**5 + x**4 - x**3 + x**2 - x + 1

    assert roots(f) == {
        -I*sin(pi/7) + cos(pi/7): 1,
        -I*sin(2*pi/7) - cos(2*pi/7): 1,
        -I*sin(3*pi/7) + cos(3*pi/7): 1,
        I*sin(pi/7) + cos(pi/7): 1,
        I*sin(2*pi/7) - cos(2*pi/7): 1,
        I*sin(3*pi/7) + cos(3*pi/7): 1,
    }

    g = ((x**2 + 1)*f**2).expand()

    assert roots(g) == {
        -I*sin(pi/7) + cos(pi/7): 2,
        -I*sin(2*pi/7) - cos(2*pi/7): 2,
        -I*sin(3*pi/7) + cos(3*pi/7): 2,
        I*sin(pi/7) + cos(pi/7): 2,
        I*sin(2*pi/7) - cos(2*pi/7): 2,
        I*sin(3*pi/7) + cos(3*pi/7): 2,
        -I: 1, I: 1,
    }

    r = roots(x**3 + 40*x + 64)
    real_root = [rx for rx in r if rx.is_extended_real][0]
    cr = 108 + 6*sqrt(1074)
    assert real_root == -2*root(cr, 3)/3 + 20/root(cr, 3)

    eq = ((7 + 5*sqrt(2))*x**3 + (-6 - 4*sqrt(2))*x**2 + (-sqrt(2) - 1)*x + 2).as_poly(x, domain=EX)
    assert roots(eq) == {-1 + sqrt(2): 1, -2 + 2*sqrt(2): 1, -sqrt(2) + 1: 1}

    eq = (41*x**5 + 29*sqrt(2)*x**5 - 153*x**4 - 108*sqrt(2)*x**4 +
          175*x**3 + 125*sqrt(2)*x**3 - 45*x**2 - 30*sqrt(2)*x**2 - 26*sqrt(2)*x -
          26*x + 24).as_poly(x, domain=EX)
    assert roots(eq) == {-sqrt(2) + 1: 1, -2 + 2*sqrt(2): 1, -1 + sqrt(2): 1,
                         -4 + 4*sqrt(2): 1, -3 + 3*sqrt(2): 1}

    eq = (x**3 - 2*x**2 + 6*sqrt(2)*x**2 - 8*sqrt(2)*x + 23*x - 14 +
          14*sqrt(2)).as_poly(x, domain=EX)
    assert roots(eq) == {-2*sqrt(2) + 2: 1, -2*sqrt(2) + 1: 1, -2*sqrt(2) - 1: 1}

    assert roots(((x + sqrt(2))**3 - 7).as_poly(x, domain=EX)) == \
        {-sqrt(2) - root(7, 3)/2 - sqrt(3)*root(7, 3)*I/2: 1,
         -sqrt(2) - root(7, 3)/2 + sqrt(3)*root(7, 3)*I/2: 1,
         -sqrt(2) + root(7, 3): 1}

    pytest.raises(PolynomialError, lambda: roots(x*y, x, y))


def test_roots1():
    assert roots(1) == {}
    assert roots(1, multiple=True) == []
    q = Symbol('q', real=True)
    assert roots(x**3 - q, x) == {cbrt(q): 1,
                                  -cbrt(q)/2 - sqrt(3)*I*cbrt(q)/2: 1,
                                  -cbrt(q)/2 + sqrt(3)*I*cbrt(q)/2: 1}
    assert roots_cubic((x**3 - 1).as_poly()) == [1, Rational(-1, 2) + sqrt(3)*I/2,
                                                 Rational(-1, 2) - sqrt(3)*I/2]

    assert roots([1, x, y]) == {-x/2 - sqrt(x**2 - 4*y)/2: 1,
                                -x/2 + sqrt(x**2 - 4*y)/2: 1}
    pytest.raises(ValueError, lambda: roots([1, x, y], z))


def test_roots_slow():
    """Just test that calculating these roots does not hang."""
    a, b, c, d, x = symbols('a,b,c,d,x')

    f1 = x**2*c + (a/b) + x*c*d - a
    f2 = x**2*(a + b*(c - d)*a) + x*a*b*c/(b*d - d) + (a*d - c/d)

    assert list(roots(f1, x).values()) == [1, 1]
    assert list(roots(f2, x).values()) == [1, 1]

    zz, yy, xx, zy, zx, yx, k = symbols('zz,yy,xx,zy,zx,yx,k')

    e1 = (zz - k)*(yy - k)*(xx - k) + zy*yx*zx + zx - zy - yx
    e2 = (zz - k)*yx*yx + zx*(yy - k)*zx + zy*zy*(xx - k)

    assert list(roots(e1 - e2, k).values()) == [1, 1, 1]

    f = x**3 + 2*x**2 + 8
    R = list(roots(f))

    assert not any(i for i in [f.subs({x: ri}).evalf(chop=True) for ri in R])


def test_roots_inexact():
    R1 = roots(x**2 + x + 1, x, multiple=True)
    R2 = roots(x**2 + x + 1.0, x, multiple=True)

    for r1, r2 in zip(R1, R2):
        assert abs(r1 - r2) < 1e-12

    f = x**4 + 3.0*sqrt(2.0)*x**3 - (78.0 + 24.0*sqrt(3.0))*x**2 \
        + 144.0*(2*sqrt(3.0) + 9.0)

    R1 = roots(f, multiple=True)
    R2 = (-12.7530479110482, -3.85012393732929,
          4.89897948556636, 7.46155167569183)

    for r1, r2 in zip(R1, R2):
        assert abs(r1 - r2) < 1e-10


def test_roots_preprocessed():
    E, F, J, L = symbols('E,F,J,L')

    f = -21601054687500000000*E**8*J**8/L**16 + \
        508232812500000000*F*x*E**7*J**7/L**14 - \
        4269543750000000*E**6*F**2*J**6*x**2/L**12 + \
        16194716250000*E**5*F**3*J**5*x**3/L**10 - \
        27633173750*E**4*F**4*J**4*x**4/L**8 + \
        14840215*E**3*F**5*J**3*x**5/L**6 + \
        54794*E**2*F**6*J**2*x**6/(5*L**4) - \
        1153*E*J*F**7*x**7/(80*L**2) + \
        633*F**8*x**8/160000

    assert roots(f, x) == {}

    R1 = roots(f.evalf(strict=False), x, multiple=True)
    R2 = [-1304.88375606366, 97.1168816800648, 186.946430171876, 245.526792947065,
          503.441004174773, 791.549343830097, 1273.16678129348, 1850.10650616851]

    w = Wild('w')
    p = w*E*J/(F*L**2)

    assert len(R1) == len(R2)

    for r1, r2 in zip(R1, R2):
        match = r1.match(p)
        assert match is not None and abs(match[w] - r2) < 1e-10


def test_roots_mixed():
    f = -1936 - 5056*x - 7592*x**2 + 2704*x**3 - 49*x**4

    _re, _im = [], []
    p = f.as_poly()
    for r in p.all_roots():
        c, (r,) = r.as_coeff_mul()
        if r.is_real:
            r = r.interval
            _re.append((c*QQ.to_expr(r.a), c*QQ.to_expr(r.b)))
        else:
            r = r.interval
            _im.append((c*QQ.to_expr(r.ax) + c*I*QQ.to_expr(r.ay),
                        c*QQ.to_expr(r.bx) + c*I*QQ.to_expr(r.by)))

    _nroots = nroots(f)
    _sroots = roots(f, multiple=True)

    _re = [Interval(a, b) for (a, b) in _re]
    _im = [Interval(re(a), re(b))*Interval(im(a), im(b)) for (a, b) in _im]

    _intervals = _re + _im
    _sroots = [r.evalf() for r in _sroots]

    _nroots = sorted(_nroots, key=lambda x: x.sort_key())
    _sroots = sorted(_sroots, key=lambda x: x.sort_key())

    for _roots in (_nroots, _sroots):
        for i, r in zip(_intervals, _roots):
            if r.is_extended_real:
                assert r in i
            else:
                assert (re(r), im(r)) in i


def test_root_factors():
    assert root_factors(Integer(1).as_poly(x)) == [Integer(1).as_poly(x)]
    assert root_factors(x.as_poly()) == [x.as_poly()]

    assert root_factors(x**2 - 1, x) == [x + 1, x - 1]
    assert root_factors(x**2 - y, x) == [x - sqrt(y), x + sqrt(y)]

    assert root_factors((x**4 - 1)**2) == \
        [x + 1, x + 1, x - 1, x - 1, x - I, x - I, x + I, x + I]

    assert root_factors((x**4 - 1).as_poly(), filter='Z') == \
        [(x + 1).as_poly(), (x - 1).as_poly(), (x**2 + 1).as_poly()]
    assert root_factors(8*x**2 + 12*x**4 + 6*x**6 + x**8, x, filter='Q') == \
        [x, x, x**6 + 6*x**4 + 12*x**2 + 8]

    pytest.raises(ValueError, lambda: root_factors((x*y).as_poly()))


def test_nroots1():
    n = 64
    p = legendre_poly(n, x, polys=True)

    pytest.raises(mpmath.mp.NoConvergence, lambda: p.nroots(n=3, maxsteps=5))

    roots = p.nroots(n=3)
    # The order of roots matters. They are ordered from smallest to the
    # largest.
    assert [str(r) for r in roots] == \
        ['-0.999', '-0.996', '-0.991', '-0.983', '-0.973', '-0.961',
         '-0.946', '-0.930', '-0.911', '-0.889', '-0.866', '-0.841',
         '-0.813', '-0.784', '-0.753', '-0.720', '-0.685', '-0.649',
         '-0.611', '-0.572', '-0.531', '-0.489', '-0.446', '-0.402',
         '-0.357', '-0.311', '-0.265', '-0.217', '-0.170', '-0.121',
         '-0.0730', '-0.0243', '0.0243', '0.0730', '0.121', '0.170',
         '0.217', '0.265', '0.311', '0.357', '0.402', '0.446', '0.489',
         '0.531', '0.572', '0.611', '0.649', '0.685', '0.720', '0.753',
         '0.784', '0.813', '0.841', '0.866', '0.889', '0.911', '0.930',
         '0.946', '0.961', '0.973', '0.983', '0.991', '0.996', '0.999']


def test_nroots2():
    p = (x**5 + 3*x + 1).as_poly()

    roots = p.nroots(n=3)
    # The order of roots matters. The roots are ordered by their real
    # components (if they agree, then by their imaginary components),
    # with real roots appearing first.
    assert [str(r) for r in roots] == \
        ['-0.332', '-0.839 - 0.944*I', '-0.839 + 0.944*I',
         '1.01 - 0.937*I', '1.01 + 0.937*I']

    roots = p.nroots(n=5)
    assert [str(r) for r in roots] == \
        ['-0.33199', '-0.83907 - 0.94385*I', '-0.83907 + 0.94385*I',
         '1.0051 - 0.93726*I', '1.0051 + 0.93726*I']


def test_roots_composite():
    assert len(roots((y**3 + y**2*sqrt(x) + y + x).as_poly(y, composite=True))) == 3


def test_sympyissue_7724():
    e = x**4*I + x**2 + I
    r1, r2 = roots(e, x), e.as_poly(x).all_roots()
    assert len(r1) == 4
    assert {_.evalf() for _ in r1} == {_.evalf() for _ in r2}


def test_sympyissue_14291():
    p = (((x - 1)**2 + 1)*((x - 1)**2 + 2)*(x - 1)).as_poly()
    assert set(p.all_roots()) == {1, 1 - I, 1 + I, 1 - I*sqrt(2), 1 + sqrt(2)*I}


def test_sympyissue_16589():
    e = x**4 - 8*sqrt(2)*x**3 + 4*x**3 - 64*sqrt(2)*x**2 + 1024*x
    rs = roots(e, x)
    assert 0 in rs
    assert all(not e.evalf(chop=True, subs={x: r}) for r in rs)


def test_sympyissue_21263():
    e = x**3 + 3*x**2 + 3*x + y + 1
    r = roots(e, x)
    assert r == {-root(y, 3) - 1: 1,
                 -root(y, 3)*(-Rational(1, 2) - sqrt(3)*I/2) - 1: 1,
                 -root(y, 3)*(-Rational(1, 2) + sqrt(3)*I/2) - 1: 1}
