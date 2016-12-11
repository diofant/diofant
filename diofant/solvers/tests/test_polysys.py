"""Tests for solvers of systems of polynomial equations. """

import pytest

from diofant import flatten, I, Integer, Poly, QQ, Rational, S, sqrt, symbols
from diofant.polys import PolynomialError, ComputationFailed
from diofant.solvers.polysys import solve_poly_system

from diofant.abc import x, y, z


def test_solve_poly_system():
    assert solve_poly_system([x - 1], x) == [{x: 1}]

    pytest.raises(ComputationFailed, lambda: solve_poly_system([0, 1]))

    assert solve_poly_system([y - x, y - x - 1], x, y) is None

    assert solve_poly_system([y - x**2, y + x**2], x, y) == [{x: 0, y: 0}]

    assert (solve_poly_system([2*x - 3, 3*y/2 - 2*x, z - 5*y], x, y, z) ==
            [{x: Rational(3, 2), y: 2, z: 10}])

    assert (solve_poly_system([x*y - 2*y, 2*y**2 - x**2], x, y) ==
            [{x: 0, y: 0}, {x: 2, y: -sqrt(2)}, {x: 2, y: sqrt(2)}])

    assert (solve_poly_system([x*y - 2*y, 2*y**2 - x**3], x, y) ==
            [{x: 0, y: 0}, {x: 2, y: -2}, {x: 2, y: 2}])

    assert (solve_poly_system([y - x**2, y + x**2 + 1], x, y) ==
            [{x: -I/sqrt(2), y: -S.Half}, {x: I/sqrt(2), y: -S.Half}])

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

    pytest.raises(NotImplementedError,
                  lambda: solve_poly_system([x**3 - y**3], x, y))
    pytest.raises(PolynomialError, lambda: solve_poly_system([1/x], x))


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
        return expr.is_Pow and expr.exp is S.Half

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


def test_solve_issue_3686():
    roots = solve_poly_system([((x - 5)**2/250000 +
                                (y - Rational(5, 10))**2/250000) - 1, x], x, y)
    assert roots == [{x: 0, y: Rational(1, 2) + 15*sqrt(1111)},
                     {x: 0, y: Rational(1, 2) - 15*sqrt(1111)}]


@pytest.mark.xfail
def test_solve_issue_3686_1():
    roots = solve_poly_system([((x - 5)**2/250000 + (y - 5.0/10)**2/250000) - 1, x], x, y)
    # TODO: does this really have to be so complicated?!
    assert len(roots) == 2
    assert roots[0][x] == 0
    assert roots[0][y].epsilon_eq(-499.474999374969, 1e12)
    assert roots[1][x] == 0
    assert roots[1][y].epsilon_eq(500.474999374969, 1e12)
