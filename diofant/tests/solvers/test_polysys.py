"""Tests for solvers of systems of polynomial equations."""

import pytest

from diofant import (ComputationFailedError, I, Integer, Matrix, Mul,
                     PolynomialError, Rational, RootOf, flatten, ordered, root,
                     sqrt, symbols)
from diofant.abc import c, n, s, t, x, y, z
from diofant.solvers.polysys import (cylindrical_algebraic_decomposition,
                                     get_nice_roots, hongproj, projone,
                                     projtwo, red, red_set,
                                     solve_linear_system, solve_poly_system,
                                     solve_poly_system_cad, solve_surd_system,
                                     subresultant_coefficients,
                                     subresultant_polynomials)


__all__ = ()


def test_solve_linear_system():
    assert solve_linear_system(Matrix([[0, 0]]), x) == {}

    M = Matrix([[0, 0, n*(n + 1), (n + 1)**2, 0],
                [n + 1, n + 1, -2*n - 1, -(n + 1), 0],
                [-1, 0, 1, 0, 0]])

    assert solve_linear_system(M, x, y, z, t) == {x: -t - t/n,
                                                  z: -t - t/n, y: 0}


def test_solve_poly_system():
    assert solve_poly_system([x - 1], x) == [{x: 1}]

    pytest.raises(ComputationFailedError, lambda: solve_poly_system([0, 1]))

    assert solve_poly_system([y - x, y - x - 1], x, y) == []

    assert solve_poly_system([x - y + 5, x + y - 3], x, y) == [{x: -1, y: 4}]
    assert solve_poly_system([x - 2*y + 5, 2*x - y - 3], x, y) == [{x: Rational(11, 3), y: Rational(13, 3)}]
    assert solve_poly_system([x**2 + y, x + y*4], x, y) == [{x: 0, y: 0}, {x: Rational(1, 4), y: Rational(-1, 16)}]

    assert solve_poly_system([y - x**2, y + x**2], x, y) == [{x: 0, y: 0}]

    assert (solve_poly_system([2*x - 3, 3*y/2 - 2*x, z - 5*y], x, y, z) ==
            [{x: Rational(3, 2), y: 2, z: 10}])

    assert (solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y) ==
            [{x: 0, y: 0}, {x: 2, y: -sqrt(2)}, {x: 2, y: sqrt(2)}])

    assert (solve_poly_system([x*y - 2*y, 2*y**2 - x**3], x, y) ==
            [{x: 0, y: 0}, {x: 2, y: -2}, {x: 2, y: 2}])

    assert (solve_poly_system([y - x**2, y + x**2 + 1], x, y) ==
            [{x: -I/sqrt(2), y: Rational(-1, 2)}, {x: I/sqrt(2), y: Rational(-1, 2)}])

    f_1 = x**2 + y + z - 1
    f_2 = x + y**2 + z - 1
    f_3 = x + y + z**2 - 1

    a, b = sqrt(2) - 1, -sqrt(2) - 1

    assert (solve_poly_system([f_1, f_2, f_3], x, y, z) ==
            [{x: 0, y: 0, z: 1}, {x: 0, y: 1, z: 0}, {x: 1, y: 0, z: 0},
             {x: a, y: a, z: a}, {x: b, y: b, z: b}])

    solution = [{x: 1, y: -1}, {x: 1, y: 1}]

    assert solve_poly_system([(x**2 - y**2).as_poly(), (x - 1).as_poly()]) == solution
    assert solve_poly_system([x**2 - y**2, x - 1], x, y) == solution
    assert solve_poly_system([x**2 - y**2, x - 1]) == solution

    assert (solve_poly_system([x + x*y - 3, y + x*y - 4], x, y) ==
            [{x: -3, y: -2}, {x: 1, y: 2}])

    assert (solve_poly_system([x**3 - y**3], x, y) ==
            [{x: y}, {x: y*(-1/2 - sqrt(3)*I/2)}, {x: y*(-1/2 + sqrt(3)*I/2)}])

    pytest.raises(PolynomialError, lambda: solve_poly_system([1/x], x))

    assert (solve_poly_system([x**6 + x - 1], x) ==
            [{x: r} for r in (x**6 + x - 1).as_poly().all_roots()])

    # Arnold's problem on two walking old women
    eqs = (4*n + 9*t - y, n*(12 - x) - 9*t, -4*n + t*(12 - x))
    res = solve_poly_system(eqs, n, t, x, y)
    assert res == [{n: 0, t: 0, y: 0}, {n: -y/2, t: y/3, x: 18}, {n: y/10, t: y/15, x: 6}]
    assert [_ for _ in res if 12 > _.get(x, 0) > 0] == [{n: y/10, t: y/15, x: 6}]
    # Now add redundant equation
    eqs = (n*(12 - x) + t*(12 - x) - y,
           4*n + 9*t - y, n*(12 - x) - 9*t, -4*n + t*(12 - x))
    res = solve_poly_system(eqs, n, x, y, t)
    assert res == [{n: 0, t: 0, y: 0}, {n: -3*t/2, x: 18, y: 3*t},
                   {n: 3*t/2, x: 6, y: 15*t}]

    assert solve_poly_system(eqs[1:], n, t, y, x) != [{n: 0, t: 0, y: 0}]


@pytest.mark.slow
def test_solve_poly_system2():
    assert solve_poly_system((x, y)) == [{x: 0, y: 0}]
    assert (solve_poly_system((x**3 + y**2,)) ==
            [{x: r} for r in (x**3 + y**2).as_poly(x).all_roots()])
    assert solve_poly_system((x, y, z)) == [{x: 0, y: 0, z: 0}]
    assert solve_poly_system((x, y, z), x, y, z, t) == [{x: 0, y: 0, z: 0}]
    assert solve_poly_system((x*y - z, y*z - x,
                              x*y - y)) == [{x: 0, y: 0, z: 0},
                                            {x: 1, y: -1, z: -1},
                                            {x: 1, y: 1, z: 1}]

    assert solve_poly_system((x + y, x - y)) == [{x: 0, y: 0}]
    assert solve_poly_system((x + y, 2*x + 2*y)) == [{x: -y}]
    assert solve_poly_system((x**2 + y**2,)) == [{x: -I*y}, {x: I*y}]
    assert (solve_poly_system((x**3*y**2 - 1,)) ==
            [{x: r} for r in (x**3*y**2 - 1).as_poly(x).all_roots()])
    assert (solve_poly_system((x**3 - y**3,)) ==
            [{x: y}, {x: y*(Rational(-1, 2) - sqrt(3)*I/2)},
             {x: y*(Rational(-1, 2) + sqrt(3)*I/2)}])
    assert solve_poly_system((y - x, y - x - 1)) == []
    assert (solve_poly_system((x*y - z**2 - z, x**2 + x - y*z,
                               x*z - y**2 - y)) ==
            [{x: -z/2 - sqrt(-3*z**2 - 2*z + 1)/2 - Rational(1, 2),
              y: -z/2 + sqrt(Mul(-1, z + 1, 3*z - 1,
                                 evaluate=False))/2 - Rational(1, 2)},
             {x: -z/2 + sqrt(-3*z**2 - 2*z + 1)/2 - Rational(1, 2),
              y: -z/2 - sqrt(Mul(-1, z + 1, 3*z - 1,
                                 evaluate=False))/2 - Rational(1, 2)},
             {x: 0, y: 0, z: 0}])

    assert solve_poly_system((x*y*z,)) == [{x: 0}, {y: 0}, {z: 0}]
    assert solve_poly_system((x**2 - 1, (x - 1)*y, (x + 1)*z)) == [{x: -1, y: 0},
                                                                   {x: 1, z: 0}]
    assert (solve_poly_system((x**2 + y**2 + z**2,
                               x + y - z, y + z**2)) ==
            [{x: -1, y: Rational(1, 2) - sqrt(3)*I/2,
              z: Rational(-1, 2) - sqrt(3)*I/2},
             {x: -1, y: Rational(1, 2) + sqrt(3)*I/2,
              z: Rational(-1, 2) + sqrt(3)*I/2}, {x: 0, y: 0, z: 0}])

    assert (solve_poly_system((x*z - 2*y + 1, y*z - 1 + z,
                               y*z + x*y*z + z)) ==
            [{x: Rational(-3, 2) - sqrt(7)*I/2,
              y: Rational(1, 4) - sqrt(7)*I/4,
              z: Rational(5, 8) + sqrt(7)*I/8},
             {x: Rational(-3, 2) + sqrt(7)*I/2,
              y: Rational(1, 4) + sqrt(7)*I/4, z: Rational(5, 8) - sqrt(7)*I/8}])

    assert solve_poly_system((x**3*y*z - x*z**2, x*y**2*z - x*y*z,
                              x**2*y**2 - z)) == [{x: 0, z: 0},
                                                  {x: -sqrt(z), y: 1},
                                                  {x: sqrt(z), y: 1},
                                                  {y: 0, z: 0}]
    assert (solve_poly_system((x*y**2 - z - z**2, x**2*y - y, y**2 - z**2)) ==
            [{y: 0, z: 0}, {x: -1, z: Rational(-1, 2), y: Rational(-1, 2)},
             {x: -1, z: Rational(-1, 2), y: Rational(1, 2)}])

    assert solve_poly_system((z*x - y - x + x*y, y*z - z + x**2 + y*x**2, x - x**2 + y, z)) == [{x: 0, y: 0, z: 0}]

    assert solve_poly_system((x*y - x*z + y**2,
                              y*z - x**2 + x**2*y,
                              x - x*y + y)) == [{x: 0, y: 0},
                                                {x: Rational(1, 2) - sqrt(3)*I/2,
                                                 y: Rational(1, 2) + sqrt(3)*I/2,
                                                 z: Rational(-1, 2) + sqrt(3)*I/2},
                                                {x: Rational(1, 2) + sqrt(3)*I/2,
                                                 y: Rational(1, 2) - sqrt(3)*I/2,
                                                 z: Rational(-1, 2) - sqrt(3)*I/2}]

    assert (solve_poly_system((y*z + x**2 + z, x*y*z + x*z - y**3,
                               x*z + y**2)) ==
            [{x: 0, y: 0, z: 0}, {x: Rational(1, 2),
                                  y: Rational(-1, 2), z: Rational(-1, 2)},
             {x: Rational(-1, 4) - sqrt(3)*I/4, y: Rational(-1, 2),
              z: Rational(1, 4) - sqrt(3)*I/4},
             {x: Rational(-1, 4) + sqrt(3)*I/4, y: Rational(-1, 2),
              z: Rational(1, 4) + sqrt(3)*I/4}])
    assert solve_poly_system((x**2 + z**2*y + y*z, y**2 - z*x + x,
                              x*y + z**2 - 1)) == [{x: 0, y: 0, z: -1},
                                                   {x: 0, y: 0, z: 1}]
    assert (solve_poly_system((x + y**2*z - 2*y**2 + 4*y - 2*z - 1,
                               -x + y**2*z - 1)) ==
            [{x: (z**3 - 2*z**2*sqrt(z**2/(z**2 - 2*z + 1)) - z**2 +
                  2*z*sqrt(z**2/(z**2 - 2*z + 1)) + 3*z - 1)/(z**2 - 2*z + 1),
              y: sqrt(z**2/(z - 1)**2) - 1/(z - 1)},
             {x: (z**3 + 2*z**2*sqrt(z**2/(z**2 - 2*z + 1)) - z**2 -
                  2*z*sqrt(z**2/(z**2 - 2*z + 1)) + 3*z - 1)/(z**2 - 2*z + 1),
              y: -sqrt(z**2/(z - 1)**2) - 1/(z - 1)}, {x: 0, y: 1, z: 1}])
    assert solve_poly_system((x, y - 1, z)) == [{x: 0, y: 1, z: 0}]

    V = A31, A32, A21, B1, B2, B3, C3, C2 = symbols('A31 A32 A21 B1 B2 B3 C3 C2')
    S = (C2 - A21, C3 - A31 - A32, B1 + B2 + B3 - 1,
         B2*C2 + B3*C3 - Rational(1, 2), B2*C2**2 + B3*C3**2 - Rational(1, 3),
         B3*A32*C2 - Rational(1, 6))
    assert (solve_poly_system(S, *V) ==
            [{A21: C2, A31: C3*(3*C2**2 - 3*C2 + C3)/(C2*(3*C2 - 2)),
              A32: C3*(C2 - C3)/(C2*(3*C2 - 2)),
              B1: (6*C2*C3 - 3*C2 - 3*C3 + 2)/(6*C2*C3),
              B2: Mul(-1, 3*C3 - 2, evaluate=False)/(6*C2*(C2 - C3)),
              B3: (3*C2 - 2)/(6*C3*(C2 - C3))},
             {A21: Rational(2, 3), A31: -1/(4*B3), A32: 1/(4*B3),
              B1: -B3 + Rational(1, 4), B2: Rational(3, 4),
              C2: Rational(2, 3), C3: 0},
             {A21: Rational(2, 3),
              A31: (8*B3 - 3)/(12*B3), A32: 1/(4*B3), B1: Rational(1, 4),
              B2: -B3 + Rational(3, 4), C2: Rational(2, 3), C3: Rational(2, 3)}])

    V = ax, bx, cx, gx, jx, lx, mx, nx, q = symbols('ax bx cx gx jx lx mx nx q')
    S = (ax*q - lx*q - mx, ax - gx*q - lx, bx*q**2 + cx*q - jx*q - nx,
         q*(-ax*q + lx*q + mx), q*(-ax + gx*q + lx))
    assert solve_poly_system(S, *V) == [{ax: (lx*q + mx)/q,
                                         bx: Mul(-1, cx*q - jx*q - nx,
                                                 evaluate=False)/q**2,
                                         gx: mx/q**2}, {ax: lx, mx: 0,
                                                        nx: 0, q: 0}]


def test_solve_biquadratic():
    x0, y0, x1, y1, r = symbols('x0 y0 x1 y1 r')

    f_1 = (x - 1)**2 + (y - 1)**2 - r**2
    f_2 = (x - 2)**2 + (y - 2)**2 - r**2

    assert (solve_poly_system([f_1, f_2], x, y) ==
            [{x: Rational(3, 2) - sqrt(-1 + 2*r**2)/2,
              y: Rational(3, 2) + sqrt(-1 + 2*r**2)/2},
             {x: Rational(3, 2) + sqrt(-1 + 2*r**2)/2,
              y: Rational(3, 2) - sqrt(-1 + 2*r**2)/2}])

    f_1 = (x - 1)**2 + (y - 2)**2 - r**2
    f_2 = (x - 1)**2 + (y - 1)**2 - r**2

    assert (solve_poly_system([f_1, f_2], x, y) ==
            [{x: 1 - sqrt((2*r - 1)*(2*r + 1))/2, y: Rational(3, 2)},
             {x: 1 + sqrt((2*r - 1)*(2*r + 1))/2, y: Rational(3, 2)}])

    def query(expr):
        return expr.is_Pow and expr.exp is Rational(1, 2)

    f_1 = (x - 1)**2 + (y - 2)**2 - r**2
    f_2 = (x - x1)**2 + (y - 1)**2 - r**2

    result = [tuple(r.values()) for r in solve_poly_system([f_1, f_2], x, y)]

    assert len(result) == 2
    assert all(len(r) == 2 for r in result)
    assert all(r.count(query) == 1 for r in flatten(result))

    f_1 = (x - x0)**2 + (y - y0)**2 - r**2
    f_2 = (x - x1)**2 + (y - y1)**2 - r**2

    result = [tuple(r.values()) for r in solve_poly_system([f_1, f_2], x, y)]

    assert len(result) == 2
    assert all(len(r) == 2 for r in result)
    assert all(len(r.find(query)) == 1 for r in flatten(result))

    eqs = [y**2 - 4 + x, y*2 + 3*x - 7]
    assert solve_poly_system(eqs, x, y) == [{x: Rational(11, 9),
                                             y: Rational(5, 3)},
                                            {x: 3, y: -1}]
    eqs = [y + x**2 - 3, -y + x - 4]
    assert solve_poly_system(eqs, x, y) == [{x: Rational(-1, 2) + sqrt(29)/2,
                                             y: Rational(-9, 2) + sqrt(29)/2},
                                            {x: -sqrt(29)/2 - Rational(1, 2),
                                             y: Rational(-9, 2) - sqrt(29)/2}]

    # issue sympy/sympy#6785
    roots = solve_poly_system([((x - 5)**2/250000 +
                                (y - Rational(5, 10))**2/250000) - 1, x],
                              x, y)
    assert roots == [{x: 0, y: Rational(1, 2) + 15*sqrt(1111)},
                     {x: 0, y: Rational(1, 2) - 15*sqrt(1111)}]

    roots = solve_poly_system([((x - 5)**2/250000 +
                                (y - 5.0/10)**2/250000) - 1, x], x, y)
    assert len(roots) == 2
    assert roots[0][x] == roots[1][x] == 0
    assert roots[0][y].epsilon_eq(-499.474999374969, 1e12)
    assert roots[1][y].epsilon_eq(+500.474999374969, 1e12)


def test_sympyissue_12345():
    eqs = (x**2 - y - sqrt(2), x**2 + x*y - y**2)
    r0, r1, r2, r3 = (y**4 - 3*y**3 + y**2*(-3*sqrt(2) + 1) +
                      2*sqrt(2)*y + 2).as_poly(y).all_roots()
    sol = [{x: sqrt(2)*r0**3/2 - 3*sqrt(2)*r0**2/2 - 2*r0 + sqrt(2)*r0/2 + 1,
            y: r0},
           {x: sqrt(2)*r1**3/2 - 3*sqrt(2)*r1**2/2 - 2*r1 + sqrt(2)*r1/2 + 1,
            y: r1},
           {x: sqrt(2)*r2**3/2 - 3*sqrt(2)*r2**2/2 - 2*r2 + sqrt(2)*r2/2 + 1,
            y: r2},
           {x: sqrt(2)*r3**3/2 - 3*sqrt(2)*r3**2/2 - 2*r3 + sqrt(2)*r3/2 + 1,
            y: r3}]
    assert solve_poly_system(eqs, x, y) == sol


@pytest.mark.slow
def test_sympyissue_16038():
    sys1 = [(2*x - 8*y)**2 + (5*x - 10*y)**2 + (10*x - 7*y)**2 - 437,
            (7*y - 10*z)**2 + (8*y - z)**2 + (10*y - 9*z)**2 - 474,
            (-10*x + 10*z)**2 + (-5*x + 9*z)**2 + (-2*x + z)**2 - 885]

    sys2 = [(2.0*x - 8*y)**2 + (5*x - 10*y)**2 + (10*x - 7*y)**2 - 437,
            (7*y - 10*z)**2 + (8*y - z)**2 + (10*y - 9*z)**2 - 474,
            (-10*x + 10*z)**2 + (-5*x + 9*z)**2 + (-2*x + z)**2 - 885]

    sols1 = solve_poly_system(sys1, x, y, z)
    sols2 = solve_poly_system(sys2, x, y, z)

    assert len(sols1) == len(sols2) == 8
    assert {x: -1, y: -2, z: -3} in sols1
    assert {x: +1, y: +2, z: +3} in sols2
    assert (list(ordered([{k: v.evalf(7) for k, v in _.items()} for _ in sols1])) ==
            list(ordered([{k: v.evalf(7) for k, v in _.items()} for _ in sols2])))


def test_solve_surd_system():
    eqs = [x + sqrt(x + 1) - 2]
    res = [{x: -sqrt(13)/2 + Rational(5, 2)}]
    assert solve_surd_system(eqs) == solve_surd_system(eqs, x) == res

    eqs = [root(x, 3) + root(x, 2) - 2]
    res = [{x: 1}]
    assert solve_surd_system(eqs) == res

    eqs = [x - (-x + 1)/(x - 1)]
    res = [{x: -1}]
    assert solve_surd_system(eqs) == res

    eqs = [root(x, 4) + root(x, 3) + sqrt(x)]
    res = [{x: 0}]
    assert solve_surd_system(eqs) == res

    a0, a1 = symbols('a:2')

    eqs = [x + sqrt(x + sqrt(x + 1)) - 2]
    res = [{x: -RootOf(a1**3 + a1**2 - 4*a1 + 1, 2) + 2}]
    assert solve_surd_system(eqs) == res

    eqs = [sqrt(x) + y + 2, y*x - 1]
    _, r1, r2 = (a0**3 + 2*a0**2 + 1).as_poly().all_roots()
    res = [{x: r1**2, y: -2 - r1}, {x: r2**2, y: -2 - r2}]
    assert solve_surd_system(eqs) == res

    eqs = [sqrt(17*x - sqrt(x**2 - 5)) - y]
    res = [{x: 17*y**2/288 - sqrt(y**4 - 1440)/288},
           {x: 17*y**2/288 + sqrt(y**4 - 1440)/288}]
    assert solve_surd_system(eqs, x) == res


def test_sympyissue_21999():
    assert solve_poly_system([x - 1], x, y) == [{x: 1}]
    assert solve_poly_system([y - 1], x, y) == [{y: 1}]


def test_sympyissue_26682():
    eqns = [c**2 + s**2 - 1,
            -1.98079646822393*c - 0.887785747630113*s - 0.15634896910398]
    assert solve_poly_system(eqns, s, c) == [{c: -0.47366192584758887,
                                              s: 0.88070675028771817},
                                             {c: 0.34220436670067722,
                                              s: -0.93962554850907942}]


# NEW TESTS FOR CAD AND RELATED FUNCTIONS

def test_red():
    # simple univar degree-one example
    assert red(x, x) == 0
    # univar degree-two
    assert red(x**2 + x + 1, x) == x+1
    # bivar
    assert red(x*y + x**2 * y**2, x) == x*y


def test_red_set():
    assert not red_set(1, x)
    assert red_set(x, x) == [x, 0]

    assert red_set(x**3 + x**2 + x + 1, x) == [x**3 + x**2 + x + 1,
                                               x**2 + x + 1,
                                               x + 1,
                                               1]

    assert red_set(y*x + y, x) == [y*x + y, y]


def test_subresultant_polynomials():
    # edge cases
    assert subresultant_polynomials(x, 0, x) == []
    assert subresultant_polynomials(x, 1, x) == [1]

    # simple monic univariate examples
    assert subresultant_polynomials(x**2+1, x**2-1, x) == [4, -2, -1+x**2]
    assert subresultant_polynomials(x**3+1, x**2-1, x) == [0, 1+x, -1+x**2]

    # should be order-invariant
    assert subresultant_polynomials(x**3+1, x**2-1, x) == subresultant_polynomials(x**3+1, x**2-1, x)

    # battery of univariate examples with diff degrees and coefficients
    # all checked in Mathematica
    fs = [
        2*x**5 - 3*x**4 + x**3 - 7*x + 5,
        -x**5 + 2*x**4 - 5*x**2 + x - 4,
        4*x**4 - x**3 + 2*x**2 - x + 3,
        x**3 - x**2 + x - 1,
        5*x**2 + 3*x - 2
    ]
    gs = [
        x**5 + 4*x**4 - x**3 + 2*x**2 - 3*x + 6,
        3*x**3 - x + 2,
        -2*x**3 + 3*x - 5,
        2*x**2 + x - 3,
        -x + 1
    ]
    answers = [
        [
            45695124,
            692022 - 809988*x,
            1349 - 743*x - 901*x**2,
            397 - 487*x + 43*x**2 - 24*x**3,
            7 + x + 4*x**2 - 3*x**3 + 11*x**4,
            6 - 3*x + 2*x**2 - x**3 + 4*x**4 + x**5
        ],
        [
            31514,
            862 - 1469*x,
            102 + 12*x + 99*x**2,
            6 - 3*x + 9*x**3
        ],
        [
            21650,
            -730 - 130*x,
            22 - 50*x + 32*x**2,
            -5 + 3*x - 2*x**3
        ],
        [
            0,
            -13 + 13*x,
            -3 + x + 2*x**2
        ],
        [
            6,
            1 - x
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_polynomials(fs[0], gs[0], x) == answers[0]
    assert subresultant_polynomials(fs[1], gs[1], x) == answers[1]
    assert subresultant_polynomials(fs[2], gs[2], x) == answers[2]
    assert subresultant_polynomials(fs[3], gs[3], x) == answers[3]
    assert subresultant_polynomials(fs[4], gs[4], x) == answers[4]

    # a battery of bivariate tests with varying deg and coeffs
    # checked with mathematica
    fs = [
        2*x**3 + y**2 + 3*x*y - 4,
        x**3 + 4*x**2*y + y - 1,
        3*x**2 + 2*x*y**2 - y + 6,
        x**3 - y**3 + x*y,
        y*x**2+1
    ]
    gs = [
        x**2 + 2*y**3 + 5*x*y,
        x**2 + y**4 - x,
        x + y**3 - 2*x*y + 1,
        x**2 - 2*y**2 + 3*x*y**2,
        y**2*x**2 - 1
    ]
    answers = [
        [
            16 + 52*y**2 + 1000*y**3 - 254*y**4 - 232*y**5 + 360*y**6 - 48*y**7 + 32*y**9,
            -4 + 3*x*y + y**2 + 50*x*y**2 - 4*x*y**3 + 20*y**4,
            x**2 + 5*x*y + 2*y**3
        ],
        [
            -5*y + 5*y**2 + 3*y**4 + 5*y**5 - 8*y**6 + 4*y**9 + 16*y**10 + y**12,
            -1 + x + y + 4*x*y - y**4 - x*y**4 - 4*y**5,
            -x + x**2 + y**4
        ],
        [
            9 - 25*y + 26*y**2 + 6*y**3 - 2*y**5 + 7*y**6,
            1 + x - 2*x*y + y**3
        ],
        [
            -2*y**4 - 8*y**5 - 4*y**6 + 27*y**9,
            x*y + 2*x*y**2 - y**3 - 6*y**4 + 9*x*y**4,
            x**2 - 2*y**2 + 3*x*y**2
        ],
        [
            y**2 + 2*y**3 + y**4,
            -y - y**2,
            x**2 - 1/y**2
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_polynomials(fs[0], gs[0], x) == answers[0]
    assert subresultant_polynomials(fs[1], gs[1], x) == answers[1]
    assert subresultant_polynomials(fs[2], gs[2], x) == answers[2]
    assert subresultant_polynomials(fs[3], gs[3], x) == answers[3]
    assert subresultant_polynomials(fs[4], gs[4], x) == answers[4]


def test_subresultant_coefficients():
    # edge cases
    assert not subresultant_coefficients(x, 0, x)
    assert subresultant_coefficients(x, 1, x) == [1]

    # simple monic univariate examples
    assert subresultant_coefficients(x**2+1, x**2-1, x) == [4, 0, 1]
    assert subresultant_coefficients(x**3+1, x**2-1, x) == [0, 1, 1]

    # should be order-invariant
    assert subresultant_coefficients(x**3+1, x**2-1, x) == subresultant_coefficients(x**3+1, x**2-1, x)

    # battery of univariate examples with diff degrees and coefficients
    # all checked in Mathematica
    fs = [
        2*x**5 - 3*x**4 + x**3 - 7*x + 5,
        -x**5 + 2*x**4 - 5*x**2 + x - 4,
        4*x**4 - x**3 + 2*x**2 - x + 3,
        x**3 - x**2 + x - 1,
        5*x**2 + 3*x - 2
    ]
    gs = [
        x**5 + 4*x**4 - x**3 + 2*x**2 - 3*x + 6,
        3*x**3 - x + 2,
        -2*x**3 + 3*x - 5,
        2*x**2 + x - 3,
        -x + 1
    ]
    answers = [
        [
            45695124,
            -809988,
            -901,
            -24,
            11,
            1
        ],
        [
            31514,
            -1469,
            99,
            9
        ],
        [
            21650,
            -130,
            32,
            -2
        ],
        [
            0,
            13,
            2
        ],
        [
            6,
            -1
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_coefficients(fs[0], gs[0], x) == answers[0]
    assert subresultant_coefficients(fs[1], gs[1], x) == answers[1]
    assert subresultant_coefficients(fs[2], gs[2], x) == answers[2]
    assert subresultant_coefficients(fs[3], gs[3], x) == answers[3]
    assert subresultant_coefficients(fs[4], gs[4], x) == answers[4]

    # a battery of bivariate tests with varying deg and coeffs
    # checked with mathematica
    fs = [
        2*x**3 + y**2 + 3*x*y - 4,
        x**3 + 4*x**2*y + y - 1,
        3*x**2 + 2*x*y**2 - y + 6,
        x**3 - y**3 + x*y,
        y*x**2+1
    ]
    gs = [
        x**2 + 2*y**3 + 5*x*y,
        x**2 + y**4 - x,
        x + y**3 - 2*x*y + 1,
        x**2 - 2*y**2 + 3*x*y**2,
        y**2*x**2 - 1
    ]
    answers = [
        [
            16 + 52*y**2 + 1000*y**3 - 254*y**4 - 232*y**5 + 360*y**6 - 48*y**7 + 32*y**9,
            3*y + 50*y**2 - 4*y**3,
            1
        ],
        [
            -5*y + 5*y**2 + 3*y**4 + 5*y**5 - 8*y**6 + 4*y**9 + 16*y**10 + y**12,
            1 + 4*y - y**4,
            1
        ],
        [
            9 - 25*y + 26*y**2 + 6*y**3 - 2*y**5 + 7*y**6,
            1 - 2*y
        ],
        [
            -2*y**4 - 8*y**5 - 4*y**6 + 27*y**9,
            y + 2*y**2 + 9*y**4,
            1
        ],
        [
            y**2 + 2*y**3 + y**4,
            0,
            1
        ]
    ]
    # unrolled it so we can see exactly which example fails
    assert subresultant_coefficients(fs[0], gs[0], x) == answers[0]
    assert subresultant_coefficients(fs[1], gs[1], x) == answers[1]
    assert subresultant_coefficients(fs[2], gs[2], x) == answers[2]
    assert subresultant_coefficients(fs[3], gs[3], x) == answers[3]
    assert subresultant_coefficients(fs[4], gs[4], x) == answers[4]


def test_get_nice_roots():
    # constants have no roots
    assert get_nice_roots(3) == []
    assert get_nice_roots(Integer(3).as_poly(x)) == []

    # if not implemented, just solve numerically
    # eg if coefficient is algebraic
    # the answer here can be solved with basic algebra
    assert get_nice_roots(sqrt(2) * x**2 - 1)[1].evalf() == sqrt(1 / sqrt(2)).evalf()

    # if roots are RootOf, then they should be numeric
    assert get_nice_roots(x**5 + x**2 - 1)[0] == RootOf(x**5 + x**2 - 1, 0).evalf()

    # the algebraic roots should stay algebraic
    # bc of the multiplication, we get the roots from x^2 - 1 of +- sqrt(2)
    assert get_nice_roots((x**2 - 2) * (x**5 - x**2 - 1)) == \
        [-sqrt(2), RootOf(x**5 - x**2 - 1, 0).evalf(), sqrt(2)]


def test_projone():
    # simple example: work it out manually by looping through
    # for x^2
    #   red_set = [x^2,0,0]
    #   LC = 1
    #   subres coeffs of x^2 and 2x (deriv) = [0,2]
    #   so add [0,1,2] to the proj factors
    # for x
    #   red_set = [x,0]
    #   LC = 1
    #   subres coeffs of x and 1 (deriv) = [1]
    #   so add [1] to the proj factors
    # hence return {0,1,2}
    assert projone([x**2, x], x) == {0, 1, 2}

    # slighly harder bivar example: loop through
    # for y*x**2
    #   red_set = [y*x^2]
    #   LC = y
    #   subres coeffs of y*x^2 and 2yx (deriv) = [0,2y]
    #   so add [0,y,2y] to the proj factors
    # for y**2*x
    #   red_set = [y**2*x]
    #   LC = y**2
    #   subres coeffs of y**2*x and y**2 (deriv) = [y**2]
    #   so add [y**2, y**2] to the proj factors
    # hence return {0,y,2y,y**2}
    assert projone([y*x**2, y**2*x], x) == {0, y**2, y, 2*y}


def test_projtwo():
    # simple example: loop through (pairs)
    # for the pair (f,g) = (x^2, x)
    #   red_set(f) = [x^2]
    #   for f_ = x^2
    #       subres coeffs of x^2 and x = [0,1]
    #   so add [0,1] to the proj factors
    # hence return {0,1}
    assert projtwo([x**2, x], x) == {0, 1}

    # slighly harder bivar example: loop through (pairs)
    # for the pair (f,g) = (y*x^2, y**2*x)
    #   red_set(f) = [y*x^2]
    #   for f_ = y*x^2
    #       subres coeffs of y*x^2 and y**2*x = [0,y**2]
    #   so add [0,y**2] to the proj factors
    # hence return {0,y**2}
    assert projtwo([y*x**2, y**2*x], x) == {0, y**2}


def test_hongproj():
    # the same two examples as in test_projone and test_projtwo
    # the tests here are really just testing the 'cleaning up'

    # this should be empty as they are all constants
    assert hongproj([x**2, x], x) == set()

    assert hongproj([y*x**2, y**2*x], x) == {y**2, y, 2*y}


def test_cylindrical_algebraic_decomposition():
    # simple univar example
    # x^2-1 has roots +-1, which implies these cells
    assert cylindrical_algebraic_decomposition([x**2-1], [x]) ==\
        [{x: val} for val in [-2, -1, 0, 1, 2]]

    # harder univar example
    # the collection of roots are (-1, 0, 1, 5**(1/3))
    # have to do the sympify thing to do algebraic comparisons
    # is there an easier way??
    assert cylindrical_algebraic_decomposition([x**2-1,
                                                x,
                                                x**3-5], [x]) ==\
        [{x: val} for val in
         [-2, -1, -Integer(1)/2, 0, Integer(1)/2, 1,
          (1 + 5**(Integer(1)/3))/2,
          5**(Integer(1)/3),
          5**(Integer(1)/3) + 1]]

    # multivar example
    # projecting on x gets [-y^2-1, y], only root is 0
    # then lift on y=-1, y=0, y=1 -- easy but tedious
    assert cylindrical_algebraic_decomposition([x + y, x*y - 1], [x, y]) == [
        {y: -1, x: -2},
        {y: -1, x: -1},
        {y: -1, x: 0},
        {y: -1, x: 1},
        {y: -1, x: 2},
        {y: 0, x: -1},
        {y: 0, x: 0},
        {y: 0, x: 1},
        {y: 1, x: -2},
        {y: 1, x: -1},
        {y: 1, x: 0},
        {y: 1, x: 1},
        {y: 1, x: 2}]


def test_solve_poly_system_cad():
    # no solution as this quadratic lies below x axis
    assert not solve_poly_system_cad([-x**2 - 1 >= 0], [x])

    # testing utility
    # solves a system and then subs the sample points back in
    def solve_and_sub(ineqs, vars, return_one_sample=True):
        soln = solve_poly_system_cad(ineqs, vars, return_one_sample)

        for sample in soln:
            if not all(ineq.subs(sample) for ineq in ineqs):
                return False
        return True

    assert solve_and_sub([x**2 - 1 >= 3], [x]) is True
    assert solve_and_sub([x**2 - 1 >= 3], [x], False) is True

    # harder example
    assert solve_and_sub([x**2 * y**2 - 1 > 0, x <= 0.2,
                          x + y >= 1], [x, y]) is True
