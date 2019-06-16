import pytest

from diofant import (Dummy, Eq, I, Integer, Integral, Rational, RootOf, cbrt,
                     cos, cosh, posify, root, solve, sqrt, sqrtdenest, sstr,
                     symbols, sympify)
from diofant.abc import F, a, b, c, d, t, x, y
from diofant.simplify.sqrtdenest import _subsets as subsets
from diofant.simplify.sqrtdenest import unrad


__all__ = ()

r2, r3, r5, r6, r7, r10, r15, r29 = [sqrt(x) for x in [2, 3, 5, 6, 7, 10,
                                                       15, 29]]


def NS(e, n=15, **options):
    return sstr(sympify(e).evalf(n, **options), full_prec=True)


def test_sqrtdenest():
    d = {sqrt(5 + 2 * r6): r2 + r3,
         sqrt(5. + 2 * r6): sqrt(5. + 2 * r6),
         sqrt(5. + 4*sqrt(5 + 2 * r6)): sqrt(5.0 + 4*r2 + 4*r3),
         sqrt(r2): sqrt(r2),
         sqrt(5 + r7): sqrt(5 + r7),
         sqrt(3 + sqrt(5 + 2*r7)):
         3*r2*root(5 + 2*r7, 4)/(2*sqrt(6 + 3*r7)) +
         r2*sqrt(6 + 3*r7)/(2*root(5 + 2*r7, 4)),
         sqrt(3 + 2*r3): 3**Rational(3, 4)*(r6/2 + 3*r2/2)/3}
    for i in d:
        assert sqrtdenest(i) == d[i]


def test_sqrtdenest2():
    assert sqrtdenest(sqrt(16 - 2*r29 + 2*sqrt(55 - 10*r29))) == \
        r5 + sqrt(11 - 2*r29)
    e = sqrt(-r5 + sqrt(-2*r29 + 2*sqrt(-10*r29 + 55) + 16))
    assert sqrtdenest(e) == root(-2*r29 + 11, 4)
    r = sqrt(1 + r7)
    assert sqrtdenest(sqrt(1 + r)) == sqrt(1 + r)
    e = sqrt(((1 + sqrt(1 + 2*sqrt(3 + r2 + r5)))**2).expand())
    assert sqrtdenest(e) == 1 + sqrt(1 + 2*sqrt(r2 + r5 + 3))

    assert sqrtdenest(sqrt(5*r3 + 6*r2)) == sqrt(2)*root(3, 4) + root(27, 4)

    assert sqrtdenest(sqrt(((1 + r5 + sqrt(1 + r3))**2).expand())) == \
        1 + r5 + sqrt(1 + r3)

    assert sqrtdenest(sqrt(((1 + r5 + r7 + sqrt(1 + r3))**2).expand())) == \
        1 + sqrt(1 + r3) + r5 + r7

    e = sqrt(((1 + cos(2) + cos(3) + sqrt(1 + r3))**2).expand())
    assert sqrtdenest(e) == cos(3) + cos(2) + 1 + sqrt(1 + r3)

    e = sqrt(-2*r10 + 2*r2*sqrt(-2*r10 + 11) + 14)
    assert sqrtdenest(e) == sqrt(-2*r10 - 2*r2 + 4*r5 + 14)

    # check that the result is not more complicated than the input
    z = sqrt(-2*r29 + cos(2) + 2*sqrt(-10*r29 + 55) + 16)
    assert sqrtdenest(z) == z

    assert sqrtdenest(sqrt(r6 + sqrt(15))) == sqrt(r6 + sqrt(15))

    z = sqrt(15 - 2*sqrt(31) + 2*sqrt(55 - 10*r29))
    assert sqrtdenest(z) == z


def test_sqrtdenest_rec():
    assert sqrtdenest(sqrt(-4*sqrt(14) - 2*r6 + 4*sqrt(21) + 33)) == \
        -r2 + r3 + 2*r7
    assert sqrtdenest(sqrt(-28*r7 - 14*r5 + 4*sqrt(35) + 82)) == \
        -7 + r5 + 2*r7
    assert sqrtdenest(sqrt(6*r2/11 + 2*sqrt(22)/11 + 6*sqrt(11)/11 + 2)) == \
        sqrt(11)*(r2 + 3 + sqrt(11))/11
    assert sqrtdenest(sqrt(468*r3 + 3024*r2 + 2912*r6 + 19735)) == \
        9*r3 + 26 + 56*r6
    z = sqrt(-490*r3 - 98*sqrt(115) - 98*sqrt(345) - 2107)
    assert sqrtdenest(z) == sqrt(-1)*(7*r5 + 7*r15 + 7*sqrt(23))
    z = sqrt(-4*sqrt(14) - 2*r6 + 4*sqrt(21) + 34)
    assert sqrtdenest(z) == z
    assert sqrtdenest(sqrt(-8*r2 - 2*r5 + 18)) == -r10 + 1 + r2 + r5
    assert sqrtdenest(sqrt(8*r2 + 2*r5 - 18)) == \
        sqrt(-1)*(-r10 + 1 + r2 + r5)
    assert sqrtdenest(sqrt(8*r2/3 + 14*r5/3 + Rational(154, 9))) == \
        -r10/3 + r2 + r5 + 3
    assert sqrtdenest(sqrt(sqrt(2*r6 + 5) + sqrt(2*r7 + 8))) == \
        sqrt(1 + r2 + r3 + r7)
    assert sqrtdenest(sqrt(4*r15 + 8*r5 + 12*r3 + 24)) == 1 + r3 + r5 + r15

    w = 1 + r2 + r3 + r5 + r7
    assert sqrtdenest(sqrt((w**2).expand())) == w
    z = sqrt((w**2).expand() + 1)
    assert sqrtdenest(z) == z

    z = sqrt(2*r10 + 6*r2 + 4*r5 + 12 + 10*r15 + 30*r3)
    assert sqrtdenest(z) == z


def test_sympyissue_6241():
    z = sqrt( -320 + 32*sqrt(5) + 64*r15)
    assert sqrtdenest(z) == z


def test_sqrtdenest3():
    z = sqrt(13 - 2*r10 + 2*r2*sqrt(-2*r10 + 11))
    assert sqrtdenest(z) == -1 + r2 + r10
    assert sqrtdenest(z, max_iter=1) == -1 + sqrt(2) + sqrt(10)
    n = sqrt(2*r6/7 + 2*r7/7 + 2*sqrt(42)/7 + 2)
    d = sqrt(16 - 2*r29 + 2*sqrt(55 - 10*r29))
    assert sqrtdenest(n/d).equals(
        r7*(1 + r6 + r7)/(7*(sqrt(-2*r29 + 11) + r5)))
    z = sqrt(sqrt(r2 + 2) + 2)
    assert sqrtdenest(z) == z
    assert sqrtdenest(sqrt(-2*r10 + 4*r2*sqrt(-2*r10 + 11) + 20)) == \
        sqrt(-2*r10 - 4*r2 + 8*r5 + 20)
    assert sqrtdenest(sqrt((112 + 70*r2) + (46 + 34*r2)*r5)) == \
        r10 + 5 + 4*r2 + 3*r5
    z = sqrt(5 + sqrt(2*r6 + 5)*sqrt(-2*r29 + 2*sqrt(-10*r29 + 55) + 16))
    r = sqrt(-2*r29 + 11)
    assert sqrtdenest(z) == sqrt(r2*r + r3*r + r10 + r15 + 5)


def test_sqrtdenest4():
    # see Denest_en.pdf in https://github.com/sympy/sympy/issues/3192
    z = sqrt(8 - r2*sqrt(5 - r5) - sqrt(3)*(1 + r5))
    z1 = sqrtdenest(z)
    c = sqrt(-r5 + 5)
    z1 = ((-r15*c - r3*c + c + r5*c - r6 - r2 + r10 + sqrt(30))/4).expand()
    assert sqrtdenest(z) == z1

    z = sqrt(2*r2*sqrt(r2 + 2) + 5*r2 + 4*sqrt(r2 + 2) + 8)
    assert sqrtdenest(z) == r2 + sqrt(r2 + 2) + 2

    w = 2 + r2 + r3 + (1 + r3)*sqrt(2 + r2 + 5*r3)
    z = sqrt((w**2).expand())
    assert sqrtdenest(z) == w.expand()


def test_sqrt_symbolic_denest():
    z = sqrt(((1 + sqrt(sqrt(2 + x) + 3))**2).expand())
    assert sqrtdenest(z) == sqrt((1 + sqrt(sqrt(2 + x) + 3))**2)
    z = sqrt(((1 + sqrt(sqrt(2 + cos(1)) + 3))**2).expand())
    assert sqrtdenest(z) == 1 + sqrt(sqrt(2 + cos(1)) + 3)
    z = ((1 + cos(2))**4 + 1).expand()
    assert sqrtdenest(z) == z
    z = sqrt(((1 + sqrt(sqrt(2 + cos(3*x)) + 3))**2 + 1).expand())
    assert sqrtdenest(z) == z
    c = cos(3)
    c2 = c**2
    assert sqrtdenest(sqrt(2*sqrt(1 + r3)*c + c2 + 1 + r3*c2)) == \
        -1 - sqrt(1 + r3)*c
    ra = sqrt(1 + r3)
    z = sqrt(20*ra*sqrt(3 + 3*r3) + 12*r3*ra*sqrt(3 + 3*r3) + 64*r3 + 112)
    assert sqrtdenest(z) == z


def test_sympyissue_5857():
    z = sqrt(1/(4*r3 + 7) + 1)
    ans = (r2 + r6)/(r3 + 2)
    assert sqrtdenest(z) == ans
    assert sqrtdenest(1 + z) == 1 + ans
    assert sqrtdenest(Integral(z + 1, (x, 1, 2))) == \
        Integral(1 + ans, (x, 1, 2))
    assert sqrtdenest(x + sqrt(y)) == x + sqrt(y)
    ans = (r2 + r6)/(r3 + 2)
    assert sqrtdenest(z) == ans
    assert sqrtdenest(1 + z) == 1 + ans
    assert sqrtdenest(Integral(z + 1, (x, 1, 2))) == \
        Integral(1 + ans, (x, 1, 2))
    assert sqrtdenest(x + sqrt(y)) == x + sqrt(y)


def test_subsets():
    assert subsets(1) == [[1]]
    assert subsets(4) == [
        [1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 0], [1, 0, 1, 0],
        [0, 1, 1, 0], [1, 1, 1, 0], [0, 0, 0, 1], [1, 0, 0, 1], [0, 1, 0, 1],
        [1, 1, 0, 1], [0, 0, 1, 1], [1, 0, 1, 1], [0, 1, 1, 1], [1, 1, 1, 1]]


def test_sympyissue_5653():
    assert sqrtdenest(
        sqrt(2 + sqrt(2 + sqrt(2)))) == sqrt(2 + sqrt(2 + sqrt(2)))


def test_unrad0():
    s = symbols('s', cls=Dummy)

    # checkers to deal with possibility of answer coming
    # back with a sign change (cf issue sympy/sympy#5203)
    def check(rv, ans):
        assert bool(rv[1]) == bool(ans[1])
        if ans[1]:
            return s_check(rv, ans)
        e = rv[0].expand()
        a = ans[0].expand()
        return e in [a, -a] and rv[1] == ans[1]

    def s_check(rv, ans):
        # get the dummy
        rv = list(rv)
        d = rv[0].atoms(Dummy)
        reps = list(zip(d, [s]*len(d)))
        # replace s with this dummy
        rv = (rv[0].subs(reps).expand(), [rv[1][0].subs(reps), rv[1][1].subs(reps)])
        ans = (ans[0].subs(reps).expand(), [ans[1][0].subs(reps), ans[1][1].subs(reps)])
        return str(rv[0]) in [str(ans[0]), str(-ans[0])] and \
            str(rv[1]) == str(ans[1])

    assert check(unrad(sqrt(x) - root(x + 1, 3)*sqrt(x + 2) + 2),
                 (s**10 + 8*s**8 + 24*s**6 - 12*s**5 - 22*s**4 - 160*s**3 - 212*s**2 -
                  192*s - 56, [s, s**2 - x]))


@pytest.mark.slow
def test_unrad1():
    pytest.raises(NotImplementedError, lambda:
                  unrad(sqrt(x) + sqrt(x + 1) + sqrt(1 - sqrt(x)) + 3))
    pytest.raises(NotImplementedError, lambda:
                  unrad(sqrt(x) + cbrt(x + 1) + 2*sqrt(y)))

    s = symbols('s', cls=Dummy)

    # checkers to deal with possibility of answer coming
    # back with a sign change (cf issue sympy/sympy#5203)
    def check(rv, ans):
        assert bool(rv[1]) == bool(ans[1])
        if ans[1]:
            return s_check(rv, ans)
        e = rv[0].expand()
        a = ans[0].expand()
        return e in [a, -a] and rv[1] == ans[1]

    def s_check(rv, ans):
        # get the dummy
        rv = list(rv)
        d = rv[0].atoms(Dummy)
        reps = list(zip(d, [s]*len(d)))
        # replace s with this dummy
        rv = (rv[0].subs(reps).expand(), [rv[1][0].subs(reps), rv[1][1].subs(reps)])
        ans = (ans[0].subs(reps).expand(), [ans[1][0].subs(reps), ans[1][1].subs(reps)])
        return str(rv[0]) in [str(ans[0]), str(-ans[0])] and \
            str(rv[1]) == str(ans[1])

    assert check(unrad(sqrt(x)),
                 (x, []))
    assert check(unrad(sqrt(x) + 1),
                 (x - 1, []))
    assert check(unrad(sqrt(x) + root(x, 3) + 2),
                 (s**3 + s**2 + 2, [s, s**6 - x]))
    assert check(unrad(sqrt(x)*root(x, 3) + 2),
                 (x**5 - 64, []))
    assert check(unrad(sqrt(x) + cbrt(x + 1)),
                 (x**3 - (x + 1)**2, []))
    assert check(unrad(sqrt(x) + sqrt(x + 1) + sqrt(2*x)),
                 (-2*sqrt(2)*x - 2*x + 1, []))
    assert check(unrad(sqrt(x) + sqrt(x + 1) + 2),
                 (16*x - 9, []))
    assert check(unrad(sqrt(x) + sqrt(x + 1) + sqrt(1 - x)),
                 (5*x**2 - 4*x, []))
    assert check(unrad(a*sqrt(x) + b*sqrt(x) + c*sqrt(y) + d*sqrt(y)),
                 ((a*sqrt(x) + b*sqrt(x))**2 - (c*sqrt(y) + d*sqrt(y))**2, []))
    assert check(unrad(sqrt(x) + sqrt(1 - x)),
                 (2*x - 1, []))
    assert check(unrad(sqrt(x) + sqrt(1 - x) - 3),
                 (x**2 - x + 16, []))
    assert check(unrad(sqrt(x) + sqrt(1 - x) + sqrt(2 + x)),
                 (5*x**2 - 2*x + 1, []))
    assert unrad(sqrt(x) + sqrt(1 - x) + sqrt(2 + x) - 3) in [
        (25*x**4 + 376*x**3 + 1256*x**2 - 2272*x + 784, []),
        (25*x**8 - 476*x**6 + 2534*x**4 - 1468*x**2 + 169, [])]
    assert unrad(sqrt(x) + sqrt(1 - x) + sqrt(2 + x) - sqrt(1 - 2*x)) == \
        (41*x**4 + 40*x**3 + 232*x**2 - 160*x + 16, [])  # orig root at 0.487
    assert check(unrad(sqrt(x) + sqrt(x + 1)), (Integer(1), []))

    eq = sqrt(x) + sqrt(x + 1) + sqrt(1 - sqrt(x))
    assert check(unrad(eq),
                 (16*x**2 - 9*x, []))
    assert {s[x] for s in solve(eq, check=False)} == {0, Rational(9, 16)}
    assert solve(eq) == []
    # but this one really does have those solutions
    assert ({s[x] for s in solve(sqrt(x) - sqrt(x + 1) +
                                 sqrt(1 - sqrt(x)))} ==
            {0, Rational(9, 16)})

    assert check(unrad(sqrt(x) + root(x + 1, 3) + 2*sqrt(y), y),
                 (2*sqrt(x)*cbrt(x + 1) + x - 4*y +
                     (x + 1)**Rational(2, 3), []))
    assert check(unrad(sqrt(x/(1 - x)) + cbrt(x + 1)),
                 (x**5 - x**4 - x**3 + 2*x**2 + x - 1, []))
    assert check(unrad(sqrt(x/(1 - x)) + 2*sqrt(y), y),
                 (4*x*y + x - 4*y, []))
    assert check(unrad(sqrt(x)*sqrt(1 - x) + 2, x),
                 (x**2 - x + 4, []))

    # http://tutorial.math.lamar.edu/
    #        Classes/Alg/SolveRadicalEqns.aspx#Solve_Rad_Ex2_a
    assert solve(Eq(x, sqrt(x + 6))) == [{x: 3}]
    assert solve(Eq(x + sqrt(x - 4), 4)) == [{x: 4}]
    assert solve(Eq(1, x + sqrt(2*x - 3))) == []
    assert {s[x] for s in solve(Eq(sqrt(5*x + 6) - 2, x))} == {-1, 2}
    assert {s[x] for s in solve(Eq(sqrt(2*x - 1) - sqrt(x - 4), 2))} == {5, 13}
    assert solve(Eq(sqrt(x + 7) + 2, sqrt(3 - x))) == [{x: -6}]
    # http://www.purplemath.com/modules/solverad.htm
    assert solve(cbrt(2*x - 5) - 3) == [{x: 16}]
    assert {s[x] for s in solve(x + 1 - root(x**4 + 4*x**3 - x, 4))} == {-Rational(1, 2),
                                                                         -Rational(1, 3)}
    assert {s[x] for s in solve(sqrt(2*x**2 - 7) - (3 - x))} == {-8, 2}
    assert solve(sqrt(2*x + 9) - sqrt(x + 1) - sqrt(x + 4)) == [{x: 0}]
    assert solve(sqrt(x + 4) + sqrt(2*x - 1) - 3*sqrt(x - 1)) == [{x: 5}]
    assert solve(sqrt(x)*sqrt(x - 7) - 12) == [{x: 16}]
    assert solve(sqrt(x - 3) + sqrt(x) - 3) == [{x: 4}]
    assert solve(sqrt(9*x**2 + 4) - (3*x + 2)) == [{x: 0}]
    assert solve(sqrt(x) - 2 - 5) == [{x: 49}]
    assert solve(sqrt(x - 3) - sqrt(x) - 3) == []
    assert solve(sqrt(x - 1) - x + 7) == [{x: 10}]
    assert solve(sqrt(x - 2) - 5) == [{x: 27}]
    assert solve(sqrt(17*x - sqrt(x**2 - 5)) - 7) == [{x: 3}]
    assert solve(sqrt(x) - sqrt(x - 1) + sqrt(sqrt(x))) == []

    # don't posify the expression in unrad and do use _mexpand
    z = sqrt(2*x + 1)/sqrt(x) - sqrt(2 + 1/x)
    p = posify(z)[0]
    assert solve(p) == []
    assert solve(z) == []
    assert solve(z + 6*I) == [{x: -Rational(1, 11)}]
    assert solve(p + 6*I) == []
    # issue sympy/sympy#8622
    assert unrad((root(x + 1, 5) - root(x, 3))) == (
        x**5 - x**3 - 3*x**2 - 3*x - 1, [])
    # issue sympy/sympy#8679
    assert check(unrad(x + root(x, 3) + root(x, 3)**2 + sqrt(y), x),
                 (s**3 + s**2 + s + sqrt(y), [s, s**3 - x]))

    # for coverage
    assert check(unrad(sqrt(x) + root(x, 3) + y),
                 (s**3 + s**2 + y, [s, s**6 - x]))
    assert solve(sqrt(x) + root(x, 3) - 2) == [{x: 1}]
    pytest.raises(NotImplementedError, lambda:
                  solve(sqrt(x) + root(x, 3) + root(x + 1, 5) - 2))
    # fails through a different code path
    pytest.raises(NotImplementedError, lambda: solve(-sqrt(2) + cosh(x)/x))
    # unrad some
    e = root(x + 1, 3) + root(x, 3)
    assert unrad(e) == (2*x + 1, [])
    eq = (sqrt(x) + sqrt(x + 1) + sqrt(1 - x) - 6*sqrt(5)/5)
    assert check(unrad(eq),
                 (15625*x**4 + 173000*x**3 + 355600*x**2 - 817920*x + 331776, []))
    assert check(unrad(root(x, 4) + root(x, 4)**3 - 1),
                 (s**3 + s - 1, [s, s**4 - x]))
    assert check(unrad(root(x, 2) + root(x, 2)**3 - 1),
                 (x**3 + 2*x**2 + x - 1, []))
    assert unrad(x**0.5) is None
    assert check(unrad(t + root(x + y, 5) + root(x + y, 5)**3),
                 (s**3 + s + t, [s, s**5 - x - y]))
    assert check(unrad(x + root(x + y, 5) + root(x + y, 5)**3, y),
                 (s**3 + s + x, [s, s**5 - x - y]))
    assert check(unrad(x + root(x + y, 5) + root(x + y, 5)**3, x),
                 (s**5 + s**3 + s - y, [s, s**5 - x - y]))
    assert check(unrad(root(x - 1, 3) + root(x + 1, 5) + root(2, 5)),
                 (s**5 + 5*root(2, 5)*s**4 + s**3 + 10*2**Rational(2, 5)*s**3 +
                  10*2**Rational(3, 5)*s**2 + 5*2**Rational(4, 5)*s + 4, [s, s**3 - x + 1]))
    pytest.raises(NotImplementedError,
                  lambda: unrad((root(x, 2) + root(x, 3) +
                                 root(x, 4)).subs({x: x**5 - x + 1})))

    # the simplify flag should be reset to False for unrad results;
    # if it's not then this next test will take a long time
    assert solve(root(x, 3) + root(x, 5) - 2) == [{x: 1}]
    eq = (sqrt(x) + sqrt(x + 1) + sqrt(1 - x) - 6*sqrt(5)/5)
    assert check(unrad(eq),
                 ((5*x - 4)*(3125*x**3 + 37100*x**2 + 100800*x - 82944), []))
    ans = [{x: Rational(4, 5)},
           {x: Rational(-1484, 375) + 172564/(140625*cbrt(114*sqrt(12657)/78125 +
                                                          Rational(12459439, 52734375))) +
               4*cbrt(114*sqrt(12657)/78125 +
                      Rational(12459439, 52734375))}]
    assert solve(eq) == ans
    # duplicate radical handling
    assert check(unrad(sqrt(x + root(x + 1, 3)) - root(x + 1, 3) - 2),
                 (s**3 - s**2 - 3*s - 5, [s, s**3 - x - 1]))
    # cov post-processing
    e = root(x**2 + 1, 3) - root(x**2 - 1, 5) - 2
    assert check(unrad(e),
                 (s**5 - 10*s**4 + 39*s**3 - 80*s**2 + 80*s - 30,
                  [s, s**3 - x**2 - 1]))

    e = sqrt(x + root(x + 1, 2)) - root(x + 1, 3) - 2
    assert check(unrad(e),
                 (s**6 - 2*s**5 - 7*s**4 - 3*s**3 + 26*s**2 + 40*s + 25,
                  [s, s**3 - x - 1]))
    assert check(unrad(e, _reverse=True),
                 (s**6 - 14*s**5 + 73*s**4 - 187*s**3 + 276*s**2 - 228*s + 89,
                  [s, s**2 - x - sqrt(x + 1)]))
    # this one needs r0, r1 reversal to work
    assert check(unrad(sqrt(x + sqrt(root(x, 3) - 1)) - root(x, 6) - 2),
                 (s**12 - 2*s**8 - 8*s**7 - 8*s**6 + s**4 + 8*s**3 + 23*s**2 +
                  32*s + 17, [s, s**6 - x]))

    # is this needed?
    # assert unrad(root(cosh(x), 3)/x*root(x + 1, 5) - 1) == (
    #    x**15 - x**3*cosh(x)**5 - 3*x**2*cosh(x)**5 - 3*x*cosh(x)**5 - cosh(x)**5, [])
    pytest.raises(NotImplementedError, lambda:
                  unrad(sqrt(cosh(x)/x) + root(x + 1, 3)*sqrt(x) - 1))
    assert unrad((x+y)**(2*y/3) + cbrt(x+y) + 1) is None
    assert check(unrad((x+y)**(2*y/3) + cbrt(x+y) + 1, x),
                 (s**(2*y) + s + 1, [s, s**3 - x - y]))

    # This tests two things: that if full unrad is attempted and fails
    # the solution should still be found; also it tests that the use of
    # composite
    assert len(solve(sqrt(y)*x + x**3 - 1, x)) == 3
    assert len(solve(-512*y**3 + 1344*cbrt(x + 2)*y**2 -
                     1176*(x + 2)**Rational(2, 3)*y -
                     169*x + 686, y, _unrad=False)) == 3

    # watch out for when the cov doesn't involve the symbol of interest
    eq = -x + (7*y/8 - cbrt(27*x/2 + 27*sqrt(x**2)/2)/3)**3 - 1
    assert solve(eq, y) == [
        {y: RootOf(-768*x + 343*y**3 - 588*cbrt(4)*y**2*cbrt(x + sqrt(x**2)) +
                   672*cbrt(2)*y*cbrt(x + sqrt(x**2))**2 - 256*sqrt(x**2) - 512, y, 0)},
        {y: RootOf(-768*x + 343*y**3 - 588*cbrt(4)*y**2*cbrt(x + sqrt(x**2)) +
                   672*cbrt(2)*y*cbrt(x + sqrt(x**2))**2 - 256*sqrt(x**2) - 512, y, 1)},
        {y: RootOf(-768*x + 343*y**3 - 588*cbrt(4)*y**2*cbrt(x + sqrt(x**2)) +
                   672*cbrt(2)*y*cbrt(x + sqrt(x**2))**2 - 256*sqrt(x**2) - 512, y, 2)}]

    eq = root(x + 1, 3) - (root(x, 3) + root(x, 5))
    assert check(unrad(eq),
                 (3*s**13 + 3*s**11 + s**9 - 1, [s, s**15 - x]))
    assert check(unrad(eq - 2),
                 (3*s**13 + 3*s**11 + 6*s**10 + s**9 + 12*s**8 + 6*s**6 + 12*s**5 +
                  12*s**3 + 7, [s, s**15 - x]))
    assert check(unrad(root(x, 3) - root(x + 1, 4)/2 + root(x + 2, 3)),
                 (4096*s**13 + 960*s**12 + 48*s**11 - s**10 - 1728*s**4,
                  [s, s**4 - x - 1]))  # orig expr has two real roots: -1, -.389
    assert check(unrad(root(x, 3) + root(x + 1, 4) - root(x + 2, 3)/2),
                 (343*s**13 + 2904*s**12 + 1344*s**11 + 512*s**10 - 1323*s**9 -
                  3024*s**8 - 1728*s**7 + 1701*s**5 + 216*s**4 - 729*s, [s, s**4 - x -
                                                                         1]))  # orig expr has one real root: -0.048
    assert check(unrad(root(x, 3)/2 - root(x + 1, 4) + root(x + 2, 3)),
                 (729*s**13 - 216*s**12 + 1728*s**11 - 512*s**10 + 1701*s**9 -
                  3024*s**8 + 1344*s**7 + 1323*s**5 - 2904*s**4 + 343*s, [s, s**4 - x -
                                                                          1]))  # orig expr has 2 real roots: -0.91, -0.15

    # orig expr has 1 real root: 19.53
    assert check(unrad(root(x, 3)/2 - root(x + 1, 4) + root(x + 2, 3) - 2),
                 (729*s**13 + 1242*s**12 + 18496*s**10 + 129701*s**9 + 388602*s**8 +
                  453312*s**7 - 612864*s**6 - 3337173*s**5 - 6332418*s**4 - 7134912*s**3
                  - 5064768*s**2 - 2111913*s - 398034, [s, s**4 - x - 1]))

    ans = solve(sqrt(x) + sqrt(x + 1) -
                sqrt(1 - x) - sqrt(2 + x))
    assert len(ans) == 1 and NS(ans[0][x])[:4] == '0.73'
    # the fence optimization problem
    # https://github.com/sympy/sympy/issues/4793#issuecomment-36994519
    eq = F - (2*x + 2*y + sqrt(x**2 + y**2))
    ans = 2*F/7 - sqrt(2)*F/14
    X = solve(eq, x, check=False)
    for xi in reversed(X):  # reverse since currently, ans is the 2nd one
        Y = solve((x*y).subs(xi).diff(y), y,
                  simplify=False, check=False)
        if any((a[y] - ans).expand().is_zero for a in Y):
            break
    else:
        assert None  # no answer was found
    assert (solve(sqrt(x + 1) + root(x, 3) - 2) ==
            [{x: (-11/(9*cbrt(Rational(47, 54) + sqrt(93)/6)) +
                  Rational(1, 3) + cbrt(Rational(47, 54) +
                                        sqrt(93)/6))**3}])
    assert (solve(sqrt(sqrt(x + 1)) + cbrt(x) - 2) ==
            [{x: (-sqrt(-2*cbrt(Rational(-1, 16) + sqrt(6913)/16) +
                        6/cbrt(Rational(-1, 16) + sqrt(6913)/16) +
                        Rational(17, 2) +
                        121/(4*sqrt(-6/cbrt(Rational(-1, 16) +
                                            sqrt(6913)/16) +
                                    2*cbrt(Rational(-1, 16) +
                                           sqrt(6913)/16) +
                                    Rational(17, 4))))/2 +
                  sqrt(-6/cbrt(Rational(-1, 16) + sqrt(6913)/16) +
                       2*cbrt(Rational(-1, 16) + sqrt(6913)/16) +
                       Rational(17, 4))/2 + Rational(9, 4))**3}])
    assert (solve(sqrt(x) + root(sqrt(x) + 1, 3) - 2) ==
            [{x: (-cbrt(Rational(81, 2) + 3*sqrt(741)/2)/3 +
                  (Rational(81, 2) + 3*sqrt(741)/2)**Rational(-1, 3) + 2)**2}])
    eq = (-x + (Rational(1, 2) - sqrt(3)*I/2)*cbrt(3*x**3/2 - x*(3*x**2 - 34)/2 +
                                                   sqrt((-3*x**3 + x*(3*x**2 - 34) + 90)**2/4 - Rational(39304, 27)) -
                                                   45) + 34/(3*(Rational(1, 2) - sqrt(3)*I/2)*cbrt(3*x**3/2 -
                                                                                                   x*(3*x**2 - 34)/2 + sqrt((-3*x**3 + x*(3*x**2 - 34) + 90)**2/4 -
                                                                                                                            Rational(39304, 27)) - 45)))
    assert check(unrad(eq),
                 (s**6 - sqrt(3)*s**6*I + 102*cbrt(12)*s**4 +
                  102*2**Rational(2, 3)*3**Rational(5, 6)*s**4*I + 1620*s**3 - 1620*sqrt(3)*s**3*I -
                  13872*cbrt(18)*s**2 + 471648 - 471648*sqrt(3)*I, [s, s**3 - 306*x
                                                                    - sqrt(3)*sqrt(31212*x**2 - 165240*x + 61484) + 810]))
    assert solve(eq, x, check=False) != []  # not other code errors


@pytest.mark.xfail
def test_unrad1_fail():
    assert solve(sqrt(x + root(x, 3)) + root(x - y, 5), y) != []


@pytest.mark.slow
def test_unrad_slow():
    # this has roots with multiplicity > 1; there should be no
    # repeats in roots obtained, however
    eq = (sqrt(1 + sqrt(1 - 4*x**2)) - x*((1 + sqrt(1 + 2*sqrt(1 - 4*x**2)))))
    assert solve(eq) == [{x: Rational(1, 2)}]


def test_unrad2():
    assert solve(root(x**3 - 3*x**2, 3) + 1 - x) == []
    assert solve(root(x + 1, 3) + root(x**2 - 2, 5) + 1) == []
