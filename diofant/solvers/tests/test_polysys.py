"""Tests for solvers of systems of polynomial equations. """

import pytest

from diofant import I, Matrix, Poly, Rational, flatten, sqrt, symbols
from diofant.abc import n, t, x, y, z
from diofant.polys import ComputationFailed, PolynomialError, RootOf
from diofant.solvers.polysys import solve_linear_system, solve_poly_system


__all__ = ()


def test_solve_linear_system():
    M = Matrix([[0, 0, n*(n + 1), (n + 1)**2, 0],
                [n + 1, n + 1, -2*n - 1, -(n + 1), 0],
                [-1, 0, 1, 0, 0]])

    assert solve_linear_system(M, x, y, z, t) == {x: -t - t/n,
                                                  z: -t - t/n, y: 0}


def test_solve_poly_system():
    assert solve_poly_system([x - 1], x) == [{x: 1}]

    pytest.raises(ComputationFailed, lambda: solve_poly_system([0, 1]))

    assert solve_poly_system([y - x, y - x - 1], x, y) == []

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

    assert solve_poly_system([Poly(x**2 - y**2), Poly(x - 1)]) == solution
    assert solve_poly_system([x**2 - y**2, x - 1], x, y) == solution
    assert solve_poly_system([x**2 - y**2, x - 1]) == solution

    assert (solve_poly_system([x + x*y - 3, y + x*y - 4], x, y) ==
            [{x: -3, y: -2}, {x: 1, y: 2}])

    assert (solve_poly_system([x**3 - y**3], x, y) ==
            [{x: y}, {x: y*(-1/2 - sqrt(3)*I/2)}, {x: y*(-1/2 + sqrt(3)*I/2)}])

    pytest.raises(PolynomialError, lambda: solve_poly_system([1/x], x))

    assert (solve_poly_system([x**6 + x - 1], x) ==
            [{x: RootOf(x**6 + x - 1, 0)}, {x: RootOf(x**6 + x - 1, 1)},
             {x: RootOf(x**6 + x - 1, 2)}, {x: RootOf(x**6 + x - 1, 3)},
             {x: RootOf(x**6 + x - 1, 4)}, {x: RootOf(x**6 + x - 1, 5)}])

    # Arnold's problem on two walking old women
    eqs = (4*n + 9*t - y, n*(12 - x) - 9*t, -4*n + t*(12 - x))
    res = solve_poly_system(eqs, n, t, x, y)
    assert res == [{n: 0, t: 0, y: 0}, {n: -y/2, t: y/3, x: 18}, {n: y/10, t: y/15, x: 6}]
    assert [_ for _ in res if 12 > _.get(x, 0) > 0] == [{n: y/10, t: y/15, x: 6}]


@pytest.mark.xfail
@pytest.mark.slow
def test_arnold_xfail():
    eqs = (n*(12 - x) + t*(12 - x) - y,
           4*n + 9*t - y, n*(12 - x) - 9*t, -4*n + t*(12 - x))
    res = solve_poly_system(eqs, n, x, y, t)
    assert res == [{n: 0, t: 0, y: 0}, {n: -3*t/2, x: 18, y: 3*t},
                   {n: 3*t/2, x: 6, y: 15*t}]

    assert solve_poly_system(eqs[1:], n, t, y, x) != [{n: 0, t: 0, y: 0}]


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
            [{x: 1 - sqrt(((2*r - 1)*(2*r + 1)))/2, y: Rational(3, 2)},
             {x: 1 + sqrt(((2*r - 1)*(2*r + 1)))/2, y: Rational(3, 2)}])

    def query(expr):
        return expr.is_Pow and expr.exp is Rational(1, 2)

    f_1 = (x - 1 )**2 + (y - 2)**2 - r**2
    f_2 = (x - x1)**2 + (y - 1)**2 - r**2

    result = [tuple(r.values()) for r in solve_poly_system([f_1, f_2], x, y)]

    assert len(result) == 2 and all(len(r) == 2 for r in result)
    assert all(r.count(query) == 1 for r in flatten(result))

    f_1 = (x - x0)**2 + (y - y0)**2 - r**2
    f_2 = (x - x1)**2 + (y - y1)**2 - r**2

    result = [tuple(r.values()) for r in solve_poly_system([f_1, f_2], x, y)]

    assert len(result) == 2 and all(len(r) == 2 for r in result)
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


def test_solve_sympyissue_6785():
    roots = solve_poly_system([((x - 5)**2/250000 +
                                (y - Rational(5, 10))**2/250000) - 1, x],
                              x, y)
    assert roots == [{x: 0, y: Rational(1, 2) + 15*sqrt(1111)},
                     {x: 0, y: Rational(1, 2) - 15*sqrt(1111)}]


@pytest.mark.xfail
def test_solve_sympyissue_6785_RR():
    roots = solve_poly_system([((x - 5)**2/250000 +
                                (y - 5.0/10)**2/250000) - 1, x], x, y)
    # TODO: does this really have to be so complicated?!
    assert len(roots) == 2
    assert roots[0][x] == 0
    assert roots[0][y].epsilon_eq(-499.474999374969, 1e12)
    assert roots[1][x] == 0
    assert roots[1][y].epsilon_eq(500.474999374969, 1e12)


def test_sympyissue_12345():
    eqs = (x**2 - y - sqrt(2), x**2 + x*y - y**2)
    r0, r1, r2, r3 = Poly(y**4 - 3*y**3 + y**2*(-3*sqrt(2) + 1) +
                          2*sqrt(2)*y + 2, y).all_roots()
    sol = [{x: sqrt(2)*r0**3/2 - 3*sqrt(2)*r0**2/2 - 2*r0 + sqrt(2)*r0/2 + 1,
            y: r0},
           {x: sqrt(2)*r1**3/2 - 3*sqrt(2)*r1**2/2 - 2*r1 + sqrt(2)*r1/2 + 1,
            y: r1},
           {x: sqrt(2)*r2**3/2 - 3*sqrt(2)*r2**2/2 - 2*r2 + sqrt(2)*r2/2 + 1,
            y: r2},
           {x: sqrt(2)*r3**3/2 - 3*sqrt(2)*r3**2/2 - 2*r3 + sqrt(2)*r3/2 + 1,
            y: r3}]
    assert solve_poly_system(eqs, x, y) == sol
