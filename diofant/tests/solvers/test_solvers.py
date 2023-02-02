import pytest

from diofant import (Derivative, E, Eq, Float, Function, I, Indexed,
                     IndexedBase, Integer, Integral, LambertW, Matrix, Max,
                     Mul, Piecewise, Pow, Rational, RootOf, Symbol, Tuple,
                     Wild, acos, arg, asin, atan, atan2, cbrt, cos, cosh, diff,
                     erf, erfc, erfcinv, erfinv, exp, expand_log, im, log, nan,
                     nfloat, ordered, pi, posify, re, real_root, root, sec,
                     sech, simplify, sin, sinh, solve, sqrt, sstr, symbols,
                     tan, tanh)
from diofant.abc import (F, a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q,
                         r, t, x, y, z)
from diofant.solvers.bivariate import _filtered_gens, _lambert, _solve_lambert
from diofant.solvers.solvers import (_invert, minsolve_linear_system,
                                     solve_linear)
from diofant.solvers.utils import checksol
from diofant.tests.core.test_evalf import NS
from diofant.utilities.randtest import verify_numerically as tn


__all__ = ()


def test_swap_back():
    f, g = map(Function, 'fg')
    fx, gx = f(x), g(x)
    assert (solve([fx + y - 2, fx - gx - 5], fx, y, gx) ==
            [{fx: gx + 5, y: -gx - 3}])
    assert solve(fx + gx**2*x - y, [fx, gx]) == [{fx: y - gx**2*x}]
    assert solve([f(1) - 2, x + 2]) == [{x: -2, f(1): 2}]


def test_solve_args():
    # equation container, issue sympy/sympy#5113
    ans = [{x: -3, y: 1}]
    eqs = (x + 5*y - 2, -3*x + 6*y - 15)
    assert all(solve(container(eqs), x, y) == ans
               for container in (tuple, list, set, frozenset))
    assert solve(Tuple(*eqs), x, y) == ans
    # implicit symbol to solve for
    assert solve(x**2 - 4) == [{x: -2}, {x: 2}]
    assert solve([x + y - 3, x - y - 5]) == [{x: 4, y: -1}]
    # no symbol to solve for
    assert solve(42) == []
    assert solve([1, 2]) == []
    # duplicate symbols removed
    assert solve((x - 3, y + 2), x, y, x) == [{x: 3, y: -2}]
    # unordered symbols
    # only 1
    assert solve(y - 3, {y}) == [{y: 3}]
    # more than 1
    assert solve(y - 3, {x, y}) == [{y: 3}]
    # multiple symbols: take the first linear solution
    assert solve(x + y - 3, [x, y]) == [{x: 3 - y}]
    # failing undetermined system
    assert (solve(a*x + b**2/(x + 4) - 3*x - 4/x, a, b) ==
            [{a: (-b**2*x + 3*x**3 + 12*x**2 + 4*x + 16)/(x**2*(x + 4))}])
    # failed single equation
    assert solve(1/(1/x - y + exp(y))) == []
    pytest.raises(NotImplementedError,
                  lambda: solve(exp(x) + sin(x) + exp(y) + sin(y)))
    # failed system
    # --  when no symbols given, 1 fails
    assert solve([y, exp(x) + x]) == [{x: -LambertW(1), y: 0}]
    #     both fail
    assert solve((exp(x) - x, exp(y) - y)) == [{x: -LambertW(-1),
                                                y: -LambertW(-1)}]
    # --  when symbols given
    assert solve([y, exp(x) + x], x, y) == [{x: -LambertW(1), y: 0}]
    # symbol is a number
    assert solve(x**2 - pi, pi) == [{pi: x**2}]
    # no equations
    assert solve([], [x]) == [{}]
    # overdetermined system
    # - nonlinear
    assert solve([(x + y)**2 - 4, x + y - 2]) == [{x: -y + 2}]
    # - linear
    assert solve((x + y - 2, 2*x + 2*y - 4)) == [{x: -y + 2}]
    # iterable with one equation
    assert solve([x - 3], x) == [{x: 3}]

    pytest.raises(ValueError, lambda: solve([x**2 * y**2 <= x**2 * y,
                                             x**2 * y**2 > x**2 * y]))

    assert solve(Eq(x, 3), x) == [{x: 3}]
    assert solve((x - 3).as_poly(), x) == [{x: 3}]

    assert solve([x - 2, x**2 + y]) == [{x: 2, y: -4}]


def test_solve_max():
    x = symbols('x', real=True)
    assert solve(1 - abs(x) - Max(-x - 2, x - 2), x) == [{x: -Rational(3, 2)},
                                                         {x: Rational(3, 2)}]


def test_solve_polynomial1():
    assert solve(x - y, x) == [{x: +y}]
    assert solve(x + y, x) == [{x: -y}]
    assert solve(3*x - 2, x) == [{x: Rational(2, 3)}]
    assert solve(Eq(3*x, 2), x) == [{x: Rational(2, 3)}]

    assert solve(x**2 - 1, x) == [{x: -1}, {x: 1}]
    assert solve(Eq(x**2, 1), x) == [{x: -1}, {x: 1}]

    assert solve(x - y**3, x) == [{x: y**3}]
    rx = root(x, 3)
    assert (solve(x - y**3, y) ==
            [{y: rx}, {y: -rx/2 - sqrt(3)*I*rx/2}, {y: -rx/2 + sqrt(3)*I*rx/2}])
    a11, a12, a21, a22, b1, b2 = symbols('a11,a12,a21,a22,b1,b2')

    assert (solve([a11*x + a12*y - b1, a21*x + a22*y - b2], x, y) ==
            [{x: (a22*b1 - a12*b2)/(a11*a22 - a12*a21),
              y: (a11*b2 - a21*b1)/(a11*a22 - a12*a21)}])

    solution = [{y: 0, x: 0}]

    assert solve((x - y, x + y), x, y) == solution
    assert solve((x - y, x + y), (x, y)) == solution
    assert solve((x - y, x + y), [x, y]) == solution

    assert solve(x**3 - 15*x - 4) == [{x: 4}, {x: -2 - sqrt(3)},
                                      {x: -2 + sqrt(3)}]

    assert solve((x**2 - 1)**2 - a, x) == [{x: -sqrt(-sqrt(a) + 1)},
                                           {x: +sqrt(-sqrt(a) + 1)},
                                           {x: -sqrt(sqrt(a) + 1)},
                                           {x: +sqrt(sqrt(a) + 1)}]

    assert (solve(x**3 - x + a, x, cubics=False) ==
            [{x: r} for r in (x**3 - x + a).as_poly(x).all_roots()])
    assert (solve(x**3 + 3*x**2 + x - 1, cubics=False) ==
            [{x: -1}, {x: -1 + sqrt(2)}, {x: -sqrt(2) - 1}])

    assert solve(x - y**2, x, y) == [{x: y**2}]
    assert solve(x**2 - y, y, x) == [{y: x**2}]

    assert solve(x*y + y, x) == [{x: -1}]
    assert solve(x*exp(y) + y, x) == [{x: -y/exp(y)}]

    assert solve(x**3 + 2*x + 3, x) == [{x: -1},
                                        {x: Rational(1, 2) - sqrt(11)*I/2},
                                        {x: Rational(1, 2) + sqrt(11)*I/2}]
    assert solve([y - x, x - 5]) == [{x: 5, y: 5}]

    assert solve(x**2 + p*x + q, x) == [{x: -p/2 - sqrt(p**2 - 4*q)/2},
                                        {x: -p/2 + sqrt(p**2 - 4*q)/2}]

    assert (solve([y**2 - x**3 + 1, x*y]) ==
            [{x: 0, y: -I}, {x: 0, y: I}, {x: 1, y: 0},
             {x: -Rational(1, 2) - sqrt(3)*I/2, y: 0},
             {x: -Rational(1, 2) + sqrt(3)*I/2, y: 0}])

    assert solve(x**5 + 5*x**4 + 10*x**3 + 10*x**2 + 5*x + 1, x) == [{x: -1}]


def test_solve_polynomial2():
    assert solve(4, x) == []
    assert solve(x, x) == [{x: 0}]
    assert solve(2*x, x) == [{x: 0}]
    assert solve(x - 3, y) == []
    assert solve(x - 3, x) == [{x: 3}]
    assert solve(x - 3) == [{x: 3}]
    assert solve(x + sqrt(2), x) == [{x: -sqrt(2)}]
    assert solve(x + root(2, 4), x) == [{x: -root(2, 4)}]
    assert solve(x**2 + 1, x) == [{x: -I}, {x: I}]
    assert solve(x**2 - 1, x) == [{x: -1}, {x: 1}]
    assert solve([x**2 - 3, y - 1]) == [{x: -sqrt(3), y: 1},
                                        {x: sqrt(3), y: 1}]
    assert solve(x**4 - 1, x) == [{x: -1}, {x: 1}, {x: -I}, {x: I}]
    assert (solve([x**2 + y - 2, y**2 - 4], x, y) ==
            [{x: -2, y: -2}, {x: 0, y: 2}, {x: 2, y: -2}])

    assert solve(z**2*x**2 - z**2*y**2) == [{x: -y}, {x: y}, {z: 0}]
    assert solve(z**2*x - z**2*y**2) == [{x: y**2}, {z: 0}]
    assert solve(z**2*x - z**2*y**2, simplify=False) == [{x: y**2}, {z: 0}]


def test_solve_polynomial3():
    # issue sympy/sympy#12114
    terms = (1 + a*b + d*e, 1 + a*c + d*f, 1 + b*c + e*f,
             g - a**2 - d**2, g - b**2 - e**2, g - c**2 - f**2)
    s = solve(terms, [a, b, c, d, e, f, g])
    assert s == [{a: -sqrt(-f**2 - 1), b: -sqrt(-f**2 - 1),
                  c: -sqrt(-f**2 - 1), d: f, e: f, g: -1},
                 {a: sqrt(-f**2 - 1), b: sqrt(-f**2 - 1),
                  c: sqrt(-f**2 - 1), d: f, e: f, g: -1},
                 {a: -sqrt(3)*f/2 - sqrt(-f**2 + 2)/2,
                  b: sqrt(3)*f/2 - sqrt(-f**2 + 2)/2, c: sqrt(-f**2 + 2),
                  d: -f/2 + sqrt(-3*f**2 + 6)/2,
                  e: -f/2 - sqrt(3)*sqrt(-f**2 + 2)/2, g: 2},
                 {a: -sqrt(3)*f/2 + sqrt(-f**2 + 2)/2,
                  b: sqrt(3)*f/2 + sqrt(-f**2 + 2)/2, c: -sqrt(-f**2 + 2),
                  d: -f/2 - sqrt(-3*f**2 + 6)/2,
                  e: -f/2 + sqrt(3)*sqrt(-f**2 + 2)/2, g: 2},
                 {a: sqrt(3)*f/2 - sqrt(-f**2 + 2)/2,
                  b: -sqrt(3)*f/2 - sqrt(-f**2 + 2)/2, c: sqrt(-f**2 + 2),
                  d: -f/2 - sqrt(-3*f**2 + 6)/2,
                  e: -f/2 + sqrt(3)*sqrt(-f**2 + 2)/2, g: 2},
                 {a: sqrt(3)*f/2 + sqrt(-f**2 + 2)/2,
                  b: -sqrt(3)*f/2 + sqrt(-f**2 + 2)/2, c: -sqrt(-f**2 + 2),
                  d: -f/2 + sqrt(-3*f**2 + 6)/2,
                  e: -f/2 - sqrt(3)*sqrt(-f**2 + 2)/2, g: 2}]


def test_solve_polynomial_cv_1a():
    # Test for solving on equations that can be converted to a polynomial
    # equation using the change of variable y -> x**Rational(p, q).
    assert solve(sqrt(x) - 1, x) == [{x: 1}]
    assert solve(sqrt(x) - 2, x) == [{x: 4}]
    assert solve(root(x, 4) - 2, x) == [{x: 16}]
    assert solve(cbrt(x) - 3, x) == [{x: 27}]
    assert solve(sqrt(x) + cbrt(x) + root(x, 4), x) == [{x: 0}]


def test_solve_polynomial_cv_1b():
    assert solve(4*x*(1 - a*sqrt(x)), x) == [{x: 0}, {x: 1/a**2}]
    assert solve(x*(root(x, 3) - 3)) == [{x: 0}, {x: 27}]


def test_solve_polynomial_cv_2():
    # Test for solving on equations that can be converted to a polynomial
    # equation multiplying both sides of the equation by x**m.
    assert (solve(x + 1/x - 1, x) in
            [[{x: Rational(1, 2) + I*sqrt(3)/2},
              {x: Rational(1, 2) - I*sqrt(3)/2}],
             [{x: Rational(1, 2) - I*sqrt(3)/2},
              {x: Rational(1, 2) + I*sqrt(3)/2}]])


def test_solve_qubics():
    assert (solve(x**3 - x + 1) ==
            [{x: -1/((-Rational(1, 2) -
                      sqrt(3)*I/2)*cbrt(3*sqrt(69)/2 +
                                        Rational(27, 2))) -
                 (-Rational(1, 2) -
                  sqrt(3)*I/2)*cbrt(3*sqrt(69)/2 +
                                    Rational(27, 2))/3},
             {x: Mul(-1, -Rational(1, 2) + sqrt(3)*I/2,
                     evaluate=False)*cbrt(3*sqrt(69)/2 +
                                          Rational(27, 2))/3 -
                 1/((-Rational(1, 2) +
                     sqrt(3)*I/2)*cbrt(3*sqrt(69)/2 +
                                       Rational(27, 2)))},
             {x: -cbrt(3*sqrt(69)/2 + Rational(27, 2))/3 -
                 1/cbrt(3*sqrt(69)/2 + Rational(27, 2))}])
    assert (solve(x**3 - x + 1, cubics=False) ==
            [{x: r} for r in (x**3 - x + 1).as_poly().all_roots()])


def test_quintics_1():
    f = x**5 - 110*x**3 - 55*x**2 + 2310*x + 979
    assert all(tn(f.subs(root).evalf(strict=False), 0)
               for root in solve(f, check=False))

    f = x**5 - 15*x**3 - 5*x**2 + 10*x + 20
    assert all(root[x].func == RootOf for root in solve(f))

    # if one uses solve to get the roots of a polynomial that has a RootOf
    # solution, make sure that the use of nfloat during the solve process
    # doesn't fail. Note: if you want numerical solutions to a polynomial
    # it is *much* faster to use nroots to get them than to solve the
    # equation only to get RootOf solutions which are then numerically
    # evaluated. So for eq = x**5 + 3*x + 7 do eq.as_poly().nroots() rather
    # than [i[x].evalf() for i in solve(eq)] to get the numerical roots of eq.
    assert (nfloat(solve(x**5 + 3*x**3 + 7)[0][x], exponent=False) ==
            RootOf(x**5 + 3*x**3 + 7, 0).evalf())


def test_highorder_poly():
    # just testing that the uniq generator is unpacked
    sol = solve(x**6 - 2*x + 2)
    assert all(isinstance(i[x], RootOf) for i in sol)
    assert len(sol) == 6


def test_quintics_2():
    f = x**5 + 15*x + 12
    assert all(tn(f.subs(root).evalf(strict=False), 0)
               for root in solve(f, check=False))

    f = x**5 - 15*x**3 - 5*x**2 + 10*x + 20
    assert all(root[x].func == RootOf for root in solve(f))


def test_solve_rational():
    assert solve((x - y**3)/((y**2)*sqrt(1 - y**2)), x) == [{x: y**3}]

    eq = x**2*(1/x - z**2/x)
    assert solve(eq, x) == []
    assert solve(eq, x, check=False) == [{x: 0}]

    assert solve(x + 1/x + y, x) == [{x: -y/2 - sqrt(y**2 - 4)/2},
                                     {x: -y/2 + sqrt(y**2 - 4)/2}]

    assert solve((x + 1)/(x**2 + 2)) == [{x: -1}]


def test_solve_nonlinear():
    assert solve(x**2 - y**2, x, y) == [{x: -y}, {x: y}]
    assert solve(x**2 - y**2) == [{x: -y}, {x: y}]
    assert solve(x**2 - y**2/exp(x), x, y) == [{x: 2*LambertW(y/2)}]
    assert solve(x**2 - y**2/exp(x), y, x) == [{y: -sqrt(E**x*x**2)},
                                               {y: sqrt(E**x*x**2)}]

    # issue sympy/sympy#5132
    assert (solve([exp(x) - sin(y), 1/y - 3], [x, y]) ==
            [{x: log(sin(Rational(1, 3))), y: Rational(1, 3)}])
    assert (solve((sqrt(x**2 + y**2) - sqrt(10), x + y - 4)) ==
            [{x: 1, y: 3}, {x: 3, y: 1}])
    eqs = [exp(x)**2 - sin(y) + z, 1/exp(y) - 3]
    assert solve(eqs) == [{x: log(-sqrt(-z - sin(log(3)))),
                           y: -log(3)},
                          {x: log(sqrt(-z - sin(log(3)))),
                           y: -log(3)}]

    # issue sympy/sympy#6056
    assert {simplify(w[x]) for w in solve(tanh(x - 1)*tanh(x + 1) + 1)} == {
        -log(2)/2 + log(1 - I),
        -log(2)/2 + log(-1 - I),
        -log(2)/2 + log(1 + I),
        -log(2)/2 + log(-1 + I), }

    # issue sympy/sympy#6248
    f = Function('f')
    assert solve(f(x + 1) - f(2*x - 1)) == [{x: 2}]
    assert solve(log(x + 1) - log(2*x - 1)) == [{x: 2}]


def test_solve_matrix():
    # issue sympy/sympy#3870
    A = Matrix(2, 2, [a, b, c, d])
    B = Matrix(2, 2, [0, 2, -3, 0])
    C = Matrix(2, 2, [1, 2, 3, 4])

    assert solve(A*B - C, [a, b, c, d]) == [{a: 1, b: -Rational(1, 3), c: 2, d: -1}]
    assert solve([A*B - C], [a, b, c, d]) == [{a: 1, b: -Rational(1, 3), c: 2, d: -1}]
    assert solve(Eq(A*B, C), [a, b, c, d]) == [{a: 1, b: -Rational(1, 3), c: 2, d: -1}]


def test_sympyissue_8666():
    assert solve(Eq(x**2 - 1/(x**2 - 4), 4 - 1/(x**2 - 4)), x) == []
    assert solve(Eq(x + 1/x, 1/x), x) == []


def test_sympyissue_7228():
    assert solve(4**(2*(x**2) + 2*x) - 8, x) == [{x: -Rational(3, 2)},
                                                 {x: Rational(1, 2)}]


def test_sympyissue_7190():
    assert solve(log(x-3) + log(x+3), x) == [{x: sqrt(10)}]


def test_linear_system():
    assert solve([x - 1, x - y, x - 2*y, y - 1], [x, y]) == []

    assert solve([x - 1, x - y, x - 2*y, x - 1], [x, y]) == []
    assert solve([x - 1, x - 1, x - y, x - 2*y], [x, y]) == []

    assert solve([x + 5*y - 2, -3*x + 6*y - 15], x, y) == [{x: -3, y: 1}]
    assert solve((x + 5*y - 2, -3*x + 6*y - 15), x, y, z) == [{x: -3, y: 1}]

    assert solve([x + y + z + t, -z - t], x, y, z, t) == [{x: -y, z: -t}]
    assert solve((x + 5*y - 2, -3*x + 6*y - z), z, x, y) == [{x: -5*y + 2,
                                                              z: 21*y - 6}]

    assert solve([x + 3, x - 3]) == []

    # issues sympy/sympy/5114 and sympy/sympy#6611
    # See that it doesn't hang; this solves in about 2 seconds.
    # Also check that the solution is relatively small.
    # Note: the system in issue sympy/sympy#6611 solves in about 5 seconds and has
    # an op-count of 138336 (with simplify=False).
    eqs = Matrix([
        [b - c/d + r/d], [c*(1/g + 1/e + 1/d) - f/g - r/d],
        [-c/g + f*(1/j + 1/i + 1/g) - h/i], [-f/i + h*(1/m + 1/l + 1/i) - k/m],
        [-h/m + k*(1/p + 1/o + 1/m) - n/p], [-k/p + n*(1/q + 1/p)]])
    v = Matrix([f, h, k, n, b, c])
    ans = solve(list(eqs), list(v), simplify=False)[0]
    # If time is taken to simplify then then 3270 below becomes
    # 3093 and the time is about 50 seconds instead of 2.
    assert sum(s.count_ops() for s in ans.values()) <= 3270


def test_linear_system_function():
    a = Function('a')
    assert solve([a(0, 0) + a(0, 1) + a(1, 0) + a(1, 1), -a(1, 0) - a(1, 1)],
                 a(0, 0), a(0, 1), a(1, 0), a(1, 1)) == [{a(1, 0): -a(1, 1), a(0, 0): -a(0, 1)}]


def test_solve_radicals():
    eq = root(x**3 - 3*x**2, 3) + 1 - x
    assert solve(eq) == []
    assert solve(eq, check=False) == [{x: Rational(1, 3)}]

    eq = root(x, 3) - root(x, 5) + Rational(1, 7)
    assert solve(eq) == [{x: RootOf(7*x**5 - 7*x**3 + 1, 1)**15},
                         {x: RootOf(7*x**5 - 7*x**3 + 1, 2)**15}]

    # XXX is this correct?
    sol = solve(eq, check=False)
    assert abs(real_root(eq.subs(sol[0])).evalf(2, strict=False)).epsilon_eq(0)

    assert solve(sqrt(x) + 1, x) == []
    assert solve(cbrt(x) + sqrt(x) + 1, x) == []
    assert solve(4*x*(1 - sqrt(x))) == [{x: 0}, {x: 1}]

    # issue diofant/diofant#925
    assert solve((sqrt(x) + 1)/(cbrt(x) + sqrt(x) + 1)) == []


# Note: multiple solutions exist for some of these equations, so the tests
# should be expected to break if the implementation of the solver changes
# in such a way that a different branch is chosen


def test_solve_transcendental():
    a, b = symbols('a, b')
    assert solve(sin(x)/x) == [{x: pi}]  # 0 is excluded
    assert solve(sin(x)/x, check=False) == [{x: 0}, {x: pi}]

    assert solve(exp(x) - 3, x) == [{x: log(3)}]
    assert solve((a*x + b)*(exp(x) - 3), x) == [{x: -b/a}, {x: log(3)}]
    assert solve(cos(x) - y, x) == [{x: -acos(y) + 2*pi}, {x: acos(y)}]
    assert solve(2*cos(x) - y, x) == [{x: -acos(y/2) + 2*pi}, {x: acos(y/2)}]
    assert solve(Eq(cos(x), sin(x)), x) == [{x: -3*pi/4}, {x: pi/4}]

    assert (solve(exp(x) + exp(-x) - y, x) in
            [[{x: log(y/2 - sqrt(y**2 - 4)/2)},
              {x: log(y/2 + sqrt(y**2 - 4)/2)}],
             [{x: log(y - sqrt(y**2 - 4)) - log(2)},
              {x: log(y + sqrt(y**2 - 4)) - log(2)}],
             [{x: log(y/2 - sqrt((y - 2)*(y + 2))/2)},
              {x: log(y/2 + sqrt((y - 2)*(y + 2))/2)}]])
    assert solve(exp(x) - 3, x) == [{x: log(3)}]
    assert solve(Eq(exp(x), 3), x) == [{x: log(3)}]
    assert solve(log(x) - 3, x) == [{x: exp(3)}]
    assert solve(sqrt(3*x) - 4, x) == [{x: Rational(16, 3)}]
    assert solve(3**(x + 2), x) == []
    assert solve(3**(2 - x), x) == []
    assert solve(x + 2**x, x) == [{x: -LambertW(log(2))/log(2)}]
    ans = solve(3*x + 5 + 2**(-5*x + 3), x)
    assert len(ans) == 1
    assert ans[0][x].expand() == -Rational(5, 3) + LambertW(-10240*root(2, 3)*log(2)/3)/(5*log(2))
    assert (solve(5*x - 1 + 3*exp(2 - 7*x), x) ==
            [{x: Rational(1, 5) + LambertW(-21*exp(Rational(3, 5))/5)/7}])
    assert (solve(2*x + 5 + log(3*x - 2), x) ==
            [{x: Rational(2, 3) + LambertW(2*exp(-Rational(19, 3))/3)/2}])
    assert solve(3*x + log(4*x), x) == [{x: LambertW(Rational(3, 4))/3}]
    assert solve((2*x + 8)*(8 + exp(x))) == [{x: -4}, {x: log(8) + pi*I}]
    eq = 2*exp(3*x + 4) - 3
    ans = solve(eq, x)  # this generated a failure in flatten
    assert len(ans) == 3
    assert all(eq.subs(a).evalf(chop=True) == 0 for a in ans)
    assert solve(2*log(3*x + 4) - 3, x) == [{x: (exp(Rational(3, 2)) - 4)/3}]
    assert solve(exp(x) + 1, x) == [{x: pi*I}]

    eq = 2*(3*x + 4)**5 - 6*7**(3*x + 9)
    result = solve(eq, x)
    ans = [{x: (log(2401) + 5*LambertW(-log(7**(7*root(3, 5)/5))))/(3*log(7))/-1}]
    assert result == ans
    # it works if expanded, too
    assert solve(eq.expand(), x) == result

    assert solve(z*cos(x) - y, x) == [{x: -acos(y/z) + 2*pi},
                                      {x: acos(y/z)}]
    assert solve(z*cos(2*x) - y, x) == [{x: -acos(y/z)/2 + pi},
                                        {x: acos(y/z)/2}]
    assert (solve(z*cos(sin(x)) - y, x) ==
            [{x: asin(acos(y/z) - 2*pi) + pi}, {x: -asin(acos(y/z)) + pi},
             {x: -asin(acos(y/z) - 2*pi)}, {x: asin(acos(y/z))}])

    assert solve(z*cos(x), x) == [{x: pi/2}, {x: 3*pi/2}]

    # issue sympy/sympy#4508
    assert solve(y - b*x/(a + x), x) in [[{x: -a*y/(y - b)}],
                                         [{x: a*y/(b - y)}]]
    assert solve(y - b*exp(a/x), x) == [{x: a/log(y/b)}]
    # issue sympy/sympy#4507
    assert solve(y - b/(1 + a*x), x) in [[{x: (b - y)/(a*y)}],
                                         [{x: -((y - b)/(a*y))}]]
    # issue sympy/sympy#4506
    assert solve(y - a*x**b, x) == [{x: (y/a)**(1/b)}]
    # issue sympy/sympy#4505
    assert solve(z**x - y, x) == [{x: log(y)/log(z)}]
    # issue sympy/sympy#4504
    assert solve(2**x - 10, x) == [{x: log(10)/log(2)}]
    # issue sympy/sympy#6744
    assert solve(x*y) == [{x: 0}, {y: 0}]
    assert solve([x*y]) == [{x: 0}, {y: 0}]
    assert solve(x**y - 1) == [{x: 1}, {y: 0}]
    assert solve([x**y - 1]) == [{x: 1}, {y: 0}]
    assert solve(x*y*(x**2 - y**2)) == [{x: 0}, {x: -y}, {x: y}, {y: 0}]
    assert solve([x*y*(x**2 - y**2)], check=False) == [{x: 0}, {x: -y}, {x: y}, {y: 0}]
    # issue sympy/sympy#4739
    assert solve(exp(log(5)*x) - 2**x, x) == [{x: 0}]

    # misc
    # make sure that the right variables is picked up in tsolve
    pytest.raises(NotImplementedError, lambda: solve((exp(x) + 1)**x - 2))

    # shouldn't generate a GeneratorsNeededError error in _tsolve when the NaN is generated
    # for eq_down. Actual answers, as determined numerically are approx. +/- 0.83
    pytest.raises(NotImplementedError, lambda:
                  solve(sinh(x)*sinh(sinh(x)) + cosh(x)*cosh(sinh(x)) - 3))

    # watch out for recursive loop in tsolve
    pytest.raises(NotImplementedError, lambda: solve((x + 2)**y*x - 3, x))

    # issue sympy/sympy#7245
    assert solve(sin(sqrt(x))) == [{x: 0}, {x: pi**2}]

    # issue sympy/sympy#7602
    a, b = symbols('a, b', extended_real=True, negative=False)
    assert sstr(solve(Eq(a, 0.5 - cos(pi*b)/2), b)) == \
        '[{b: -0.318309886183791*acos(-2.0*a + 1.0) + 2.0}, {b: 0.318309886183791*acos(-2.0*a + 1.0)}]'

    expr = root(x, 3) - root(x, 5)
    expr1 = root(x, 3, 1) - root(x, 5, 1)
    v = expr1.subs({x: -3})
    eq = Eq(expr, v)
    eq1 = Eq(expr1, v)
    assert solve(eq, check=False) == [{x: _**15}
                                      for _ in (-x**5 + x**3 +
                                                v).as_poly(x, extension=False).all_roots()]
    for s, v in zip((x.subs(_) for _ in solve(eq1, check=False)),
                    (_**15 for _ in ((-1)**Rational(2, 3)*x**5 -
                                     (-1)**Rational(2, 5)*x**3 -
                                     v).as_poly(x, extension=False).all_roots())):
        assert simplify(s - v) == 0

    assert solve(+3**x - 10) == [{x: log(10)/log(3)}]
    assert solve(-3**x + 10) == [{x: log(10)/log(3)}]

    assert solve(a*x**b - y, x) == [{x: (y/a)**(1/b)}]


def test_solve_for_exprs():
    f = Function('f')
    df = f(x).diff(x)

    assert solve(f(x) - x, f(x)) == [{f(x): x}]
    assert (solve(f(x).diff(x) - f(x) - x, f(x).diff(x)) ==
            [{f(x).diff(x): x + f(x)}])
    assert solve(f(x).diff(x) - f(x) - x, f(x)) == [{f(x): -x + df}]

    A = IndexedBase('A')
    eqs = Tuple(A[1] + A[2] - 3, A[1] - A[2] + 1)
    assert solve(eqs, eqs.atoms(Indexed)) == [{A[1]: 1, A[2]: 2}]

    assert solve(x + 2 + sqrt(3), x + 2) == [{x + 2: -sqrt(3)}]
    assert (solve((x + 2 + sqrt(3), x + 4 + y), y, x + 2) ==
            [{y: -2 + sqrt(3), x + 2: -sqrt(3)}])

    eqs = (x*y + 3*y + sqrt(3), x + 4 + y)
    assert solve(eqs, y, x + 2) == [{y: -sqrt(3)/(x + 3),
                                     x + 2: (-2*x - 6 +
                                             sqrt(3))/(x + 3)}]
    assert solve(eqs, y*x, x) == [{x: -y - 4, x*y: -3*y - sqrt(3)}]

    assert solve(sqrt(2) - 1, 1, check=False) == [{1: sqrt(2)}]
    assert solve(x - y + 1, 1) == [{1: x/(y - 1)}]  # /!\ -1 is targeted, too
    assert [_[1].subs({z: -1})
            for _ in solve((x - y + 1).subs({-1: z}), 1)] == [-x + y]

    assert solve([x - 2, x**2 + f(x)], {f(x), x}) == [{x: 2, f(x): -4}]

    eq = f(x).diff(x) + f(x).diff(x, x) - f(x).diff(x)/f(x)
    assert solve(eq, f(x).diff(x, x)) == [{f(x).diff(x, x):
                                           (-f(x) + 1)*f(x).diff(x)/f(x)}]


def test_solve_for_functions_derivatives():
    t = Symbol('t')
    x = Function('x')(t)
    y = Function('y')(t)
    a11, a12, a21, a22, b1, b2 = symbols('a11,a12,a21,a22,b1,b2')

    soln = solve([a11*x + a12*y - b1, a21*x + a22*y - b2], x, y)
    assert soln == [{x: (a22*b1 - a12*b2)/(a11*a22 - a12*a21),
                     y: (a11*b2 - a21*b1)/(a11*a22 - a12*a21)}]

    assert solve(x - 1, x) == [{x: 1}]
    assert solve(3*x - 2, x) == [{x: Rational(2, 3)}]

    soln = solve([a11*x.diff(t) + a12*y.diff(t) - b1,
                  a21*x.diff(t) + a22*y.diff(t) - b2],
                 x.diff(t), y.diff(t))
    assert soln == [{y.diff(t): (a11*b2 - a21*b1)/(a11*a22 - a12*a21),
                     x.diff(t): (a22*b1 - a12*b2)/(a11*a22 - a12*a21)}]

    assert solve(x.diff(t) - 1, x.diff(t)) == [{x.diff(t): 1}]
    assert solve(3*x.diff(t) - 2, x.diff(t)) == [{x.diff(t): Rational(2, 3)}]

    eqns = {3*x - 1, 2*y - 4}
    assert solve(eqns, {x, y}) == [{x: Rational(1, 3), y: 2}]
    x = Symbol('x')
    f = Function('f')
    F = x**2 + f(x)**2 - 4*x - 1
    assert solve(F.diff(x), diff(f(x), x)) == [{diff(f(x), x): (-x + 2)/f(x)}]

    # Mixed cased with a Symbol and a Function
    x = Symbol('x')
    y = Function('y')(t)

    soln = solve([a11*x + a12*y.diff(t) - b1, a21*x +
                  a22*y.diff(t) - b2], x, y.diff(t))
    assert soln == [{y.diff(t): (a11*b2 - a21*b1)/(a11*a22 - a12*a21),
                     x: (a22*b1 - a12*b2)/(a11*a22 - a12*a21)}]


def test_sympyissue_3725():
    f = Function('f')
    F = x**2 + f(x)**2 - 4*x - 1
    e = F.diff(x)
    assert solve(e, f(x).diff(x)) in [[{f(x).diff(x): (2 - x)/f(x)}],
                                      [{f(x).diff(x): -((x - 2)/f(x))}]]


def test_sympyissue_3870():
    A = Matrix(2, 2, [a, b, c, d])
    B = Matrix(2, 2, [0, 2, -3, 0])
    C = Matrix(2, 2, [1, 2, 3, 4])

    assert solve([A*B - B*A], [a, b, c, d]) == [{a: d, b: -Rational(2, 3)*c}]
    assert solve([A*C - C*A], [a, b, c, d]) == [{a: d - c, b: Rational(2, 3)*c}]
    assert solve([A*B - B*A, A*C - C*A], [a, b, c, d]) == [{a: d, b: 0, c: 0}]

    assert solve([Eq(A*B, B*A)], [a, b, c, d]) == [{a: d, b: -Rational(2, 3)*c}]
    assert solve([Eq(A*C, C*A)], [a, b, c, d]) == [{a: d - c, b: Rational(2, 3)*c}]
    assert solve([Eq(A*B, B*A), Eq(A*C, C*A)], [a, b, c, d]) == [{a: d, b: 0, c: 0}]


def test_solve_linear():
    w = Wild('w')
    assert solve_linear(Integer(0), x) == (0, 1)
    assert solve_linear(3*x - y, x) == (x, y/3)
    assert solve_linear(3*x - y, y) == (y, 3*x)
    assert solve_linear(x**2/y - 1, y) == (y, x**2)
    assert solve_linear(w - x, x) == (x, w)
    assert solve_linear(cos(x)**2 + sin(x)**2 + 2 + y,
                        y) == (y, -2 - cos(x)**2 - sin(x)**2)
    assert solve_linear(cos(x)**2 + sin(x)**2 + 2 + y, x) == (0, 1)
    assert solve_linear(x - 3, x) == (x, 3)
    assert solve_linear(1/(1/x - 2), x) == (0, 0)
    assert solve_linear((x + 1)*exp(-x), x) == (x, -1)
    assert solve_linear((x + 1)*exp(x), x) == ((x + 1)*exp(x), 1)
    assert solve_linear(x*exp(-x**2), x) == (x, 0)
    assert solve_linear(0**x - 1, x) == (0**x - 1, 1)
    pytest.raises(ValueError, lambda: solve_linear(x - 3, Integer(3)))
    assert solve_linear(x + y**2, x) == (x, -y**2)
    assert solve_linear(1/x - y**2, x) == (x, 1/y**2)
    assert solve_linear(x**2/y**2 - 3, x) == (x**2 - 3*y**2, y**2)
    assert solve_linear(1/(1/x), x) == (x, 0)
    assert solve_linear(x**2*(1/x - z**2/x), x) == (x**2*(-z**2 + 1), x)
    assert solve_linear(x + y + z, y) == (y, -x - z)


def test_sympyissue_4793():
    assert solve(1/x) == []
    assert solve(x*(1 - 5/x)) == [{x: 5}]
    assert solve(x + sqrt(x) - 2) == [{x: 1}]
    assert solve(-(1 + x)/(2 + x)**2 + 1/(2 + x)) == []
    assert solve(-x**2 - 2*x + (x + 1)**2 - 1) == []
    assert solve((x/(x + 1) + 3)**(-2)) == []
    assert solve(x/sqrt(x**2 + 1), x) == [{x: 0}]
    assert solve(exp(x) - y, x) == [{x: log(y)}]
    assert solve(exp(x)) == []
    assert solve(x**2 + x + sin(y)**2 + cos(y)**2 - 1, x) in [[{x: 0}, {x: -1}], [{x: -1}, {x: 0}]]
    eq = 4*3**(5*x + 2) - 7
    ans = solve(eq, x)
    assert len(ans) == 5
    assert all(eq.subs(a).evalf(chop=True) == 0 for a in ans)
    assert solve(log(x**2) - y**2/exp(x),
                 x, y) == [{y: -sqrt(exp(x)*log(x**2))},
                           {y: sqrt(exp(x)*log(x**2))}]
    assert solve(x**2*z**2 - z**2*y**2) == [{x: -y}, {x: y}, {z: 0}]
    assert solve((x - 1)/(1 + 1/(x - 1))) == []
    assert solve(x**(y*z) - x, x) == [{x: 1}]
    pytest.raises(NotImplementedError, lambda: solve(log(x) - exp(x), x))
    pytest.raises(NotImplementedError, lambda: solve(2**x - exp(x) - 3))


def test_sympypull_1964():
    # issue sympy/sympy#5171
    assert solve(sqrt(x)) == solve(sqrt(x**3)) == [{x: 0}]
    assert solve(sqrt(x - 1)) == [{x: 1}]
    # issue sympy/sympy#4462
    a = Symbol('a')
    assert solve(-3*a/sqrt(x), x) == []
    # issue sympy/sympy#4486
    assert solve(2*x/(x + 2) - 1, x) == [{x: 2}]
    # issue sympy/sympy#4496
    assert solve((x**2/(7 - x)).diff(x)) == [{x: 0}, {x: 14}]
    # issue sympy/sympy#4695
    f = Function('f')
    assert solve((3 - 5*x/f(x))*f(x), f(x)) == [{f(x): 5*x/3}]
    # issue sympy/sympy#4497
    assert solve(1/root(5 + x, 5) - 9, x) == [{x: Rational(-295244, 59049)}]

    assert solve(sqrt(x) + sqrt(sqrt(x)) - 4) == [{x: (Rational(-1, 2) + sqrt(17)/2)**4}]
    assert (solve((sqrt(exp(x)) + sqrt(exp(-x)) - 4).as_poly()) in
            [[{x: log((-sqrt(3) + 2)**2)}, {x: log((sqrt(3) + 2)**2)}],
             [{x: 2*log(-sqrt(3) + 2)}, {x: 2*log(sqrt(3) + 2)}],
             [{x: log(-4*sqrt(3) + 7)}, {x: log(4*sqrt(3) + 7)}]])
    assert solve((exp(x) + exp(-x) - 4).as_poly()) == [{x: log(-sqrt(3) + 2)},
                                                       {x: log(+sqrt(3) + 2)}]
    assert (solve(x**y + x**(2*y) - 1, x) ==
            [{x: (Rational(-1, 2) + sqrt(5)/2)**(1/y)},
             {x: (Rational(-1, 2) - sqrt(5)/2)**(1/y)}])

    assert solve(exp(x/y)*exp(-z/y) - 2, y) == [{y: (x - z)/log(2)}]
    assert solve(x**z*y**z - 2, z) in [[{z: log(2)/(log(x) + log(y))}],
                                       [{z: log(2)/(log(x*y))}]]
    # if you do inversion too soon then multiple roots (as for the following)
    # will be missed, e.g. if exp(3*x) = exp(3) -> 3*x = 3
    assert (solve(exp(3*x) - exp(3), x) in
            [[{x: 1}, {x: log(E*(Rational(-1, 2) - sqrt(3)*I/2))},
              {x: log(E*(Rational(-1, 2) + sqrt(3)*I/2))}],
             [{x: 1}, {x: log(-E/2 - sqrt(3)*E*I/2)},
              {x: log(-E/2 + sqrt(3)*E*I/2)}]])

    # coverage test
    p = Symbol('p', positive=True)
    assert solve((1/p + 1)**(p + 1)) == []


def test_sympyissue_5197():
    x = Symbol('x', extended_real=True)
    assert solve(x**2 + 1, x) == []
    n = Symbol('n', integer=True, positive=True)
    assert solve((n - 1)*(n + 2)*(2*n - 1), n) == [{n: 1}]
    x = Symbol('x', positive=True)
    y = Symbol('y')
    assert solve([x + 5*y - 2, -3*x + 6*y - 15], x, y) == []
    y = Symbol('y', positive=True)
    # The solution following should not contain {y: -x*exp(x/2)}
    assert solve(x**2 - y**2/exp(x), y, x) == [{y: x*exp(x/2)}]
    assert solve(x**2 - y**2/exp(x), x, y) == [{x: 2*LambertW(y/2)}]
    x, y, z = symbols('x y z', positive=True)
    assert solve(z**2*x**2 - z**2*y**2/exp(x), y, x, z) == [{y: x*exp(x/2)}]


def test_checking():
    assert solve(x*(x - y/x), x, check=False) == [{x: 0}, {x: -sqrt(y)},
                                                  {x: sqrt(y)}]
    assert solve(x*(x - y/x), x) == [{x: -sqrt(y)}, {x: sqrt(y)}]
    # {x: 0, y: 4} sets denominator to 0 in the following so system should return None
    assert solve((1/(1/x + 2), 1/(y - 3) - 1)) == []
    # 0 sets denominator of 1/x to zero so None is returned
    assert solve(1/(1/x + 2)) == []


def test_sympyissue_4671():
    assert (solve((2**exp(y**2/x) + 2)/(x**2 + 15), y) ==
            [{y: -sqrt(x)*sqrt(-log(log(2)) + log(log(2) + I*pi))},
             {y: sqrt(x)*sqrt(-log(log(2)) + log(log(2) + I*pi))}])
    assert (solve(sqrt(x**2 - 1) - 2) in
            ([{x: sqrt(5)}, {x: -sqrt(5)}], [{x: -sqrt(5)}, {x: sqrt(5)}]))


def test_sympyissue_4467():
    assert solve((a**2 + 1)*(sin(a*x) + cos(a*x)), x) == [{x: -pi/(4*a)},
                                                          {x: 3*pi/(4*a)}]


def test_sympyissue_4463():
    C1, C2 = symbols('C1 C2')
    f = Function('f')

    assert solve(C1 + C2/x**2 - exp(-f(x)),
                 f(x)) == [{f(x): log(x**2/(C1*x**2 + C2))}]
    assert (solve(1 - log(a + 4*x**2), x) in
            ([{x: -sqrt(-a + E)/2}, {x: sqrt(-a + E)/2}],
             [{x: sqrt(-a + E)/2}, {x: -sqrt(-a + E)/2}]))
    assert (solve(log(a**(-3) - x**2)/a, x) in
            ([{x: -sqrt(-1 + a**(-3))}, {x: sqrt(-1 + a**(-3))}],
             [{x: sqrt(-1 + a**(-3))}, {x: -sqrt(-1 + a**(-3))}]))
    assert (solve(1 - log(a + 4*x**2), x) in
            ([{x: -sqrt(-a + E)/2}, {x: sqrt(-a + E)/2}],
             [{x: sqrt(-a + E)/2}, {x: -sqrt(-a + E)/2}]))
    assert solve(3 - (sinh(a*x) + cosh(a*x)), x) == [{x: log(3)/a}]
    assert (solve(3 - (sinh(a*x) + cosh(a*x)**2), x) ==
            [{x: (log(-1 + sqrt(2)) + I*pi)/a}, {x: (log(2 + sqrt(5)) + I*pi)/a},
             {x: log(-2 + sqrt(5))/a}, {x: log(1 + sqrt(2))/a}])
    assert solve(atan(x) - 1) == [{x: tan(1)}]
    assert solve(-a*x + 2*x*log(x), x) == [{x: exp(a/2)}]
    assert solve(a/x + exp(x/2), x) == [{x: 2*LambertW(-a/2)}]
    assert solve(x**x) == []
    assert solve(x**x - 2) == [{x: exp(LambertW(log(2)))}]
    assert solve(((x - 3)*(x - 2))**((x - 3)*(x - 4))) == [{x: 2}]
    assert solve((a/x + exp(x/2)).diff(x), x) == [{x: 4*LambertW(sqrt(2)*sqrt(a)/4)}]


def test_sympyissue_5132():
    assert (solve([r - x**2 - y**2, tan(t) - y/x], [x, y]) ==
            [{x: -sqrt(r*sin(t)**2)/tan(t), y: -sqrt(r*sin(t)**2)},
             {x: sqrt(r*sin(t)**2)/tan(t), y: sqrt(r*sin(t)**2)}])
    assert (solve([exp(x) - sin(y), 1/exp(y) - 3], [x, y]) ==
            [{x: log(-sin(log(3))), y: -log(3)}])
    assert (solve([exp(x) - sin(y), y**2 - 4]) ==
            [{x: log(-sin(2)), y: -2}, {x: log(sin(2)), y: 2}])
    eqs = [exp(x)**2 - sin(y) + z**2, 1/exp(y) - 3]
    assert solve(eqs, [x, y]) == [{x: log(-sqrt(-z**2 - sin(log(3)))),
                                   y: -log(3)},
                                  {x: log(sqrt(-z**2 - sin(log(3)))),
                                   y: -log(3)}]
    assert solve(eqs, x, z) == [{x: log(-sqrt(-z**2 + sin(y)))},
                                {x: log(sqrt(-z**2 + sin(y)))}]
    assert (solve(eqs, x, y) ==
            [{x: log(-sqrt(-z**2 - sin(log(3)))), y: -log(3)},
             {x: log(sqrt(-z**2 - sin(log(3)))), y: -log(3)}])
    assert (solve(eqs, y, z) ==
            [{y: -log(3), z: -sqrt(-exp(2*x) - sin(log(3)))},
             {y: -log(3), z: +sqrt(-exp(2*x) - sin(log(3)))}])
    eqs = [exp(x)**2 - sin(y) + z, 1/exp(y) - 3]
    assert solve(eqs, x, z) == [{x: log(-sqrt(-z + sin(y)))},
                                {x: log(sqrt(-z + sin(y)))}]
    assert (solve(eqs, x, y) ==
            [{x: log(-sqrt(-z - sin(log(3)))), y: -log(3)},
             {x: log(+sqrt(-z - sin(log(3)))), y: -log(3)}])
    assert solve(eqs, z, y) == [{z: -exp(2*x) - sin(log(3)), y: -log(3)}]
    assert (solve((sqrt(x**2 + y**2) - sqrt(10), x + y - 4), x, y) ==
            [{x: 1, y: 3}, {x: 3, y: 1}])


@pytest.mark.slow
def test_sympyissue_5335():
    lam, a0, conc = symbols('lam a0 conc')
    eqs = [lam + 2*y - a0*(1 - x/2)*x - 0.005*x/2*x,
           a0*(1 - x/2)*x - 1*y - 0.743436700916726*y,
           x + y - conc]
    sym = [x, y, a0]
    # there are 4 solutions but only two are valid
    assert len(solve(eqs, sym, simplify=False, check=False)) == 2


def test_sympyissue_5767():
    assert solve([x**2 + y + 4], [x]) == [{x: -sqrt(-y - 4)},
                                          {x: +sqrt(-y - 4)}]


def test_polysys():
    assert (solve([x**2 + 2/y - 2, x + y - 3], [x, y]) ==
            [{x: 1, y: 2}, {x: 1 + sqrt(5), y: 2 - sqrt(5)},
             {x: 1 - sqrt(5), y: 2 + sqrt(5)}])
    assert solve([x**2 + y - 2, x**2 + y]) == []
    assert (solve([x**2 + y - 3, x - y - 4], (x, y)) ==
            [{x: -Rational(1, 2) + sqrt(29)/2, y: -Rational(9, 2) + sqrt(29)/2},
             {x: -sqrt(29)/2 - Rational(1, 2), y: -Rational(9, 2) - sqrt(29)/2}])


def test__invert():
    assert _invert(x - 2) == (2, x)
    assert _invert(2) == (2, 0)
    assert _invert(exp(1/x) - 3, x) == (1/log(3), x)
    assert _invert(exp(1/x + a/x) - 3, x) == ((a + 1)/log(3), x)
    assert _invert(a, x) == (a, 0)


@pytest.mark.slow
def test_sympyissue_5114():
    # there is no 'a' in the equation set but this is how the
    # problem was originally posed
    syms = a, b, c, f, h, k, n
    eqs = [b + r/d - c/d, c*(1/d + 1/e + 1/g) - f/g - r/d,
           f*(1/g + 1/i + 1/j) - c/g - h/i,
           h*(1/i + 1/l + 1/m) - f/i - k/m,
           k*(1/m + 1/o + 1/p) - h/m - n/p,
           n*(1/p + 1/q) - k/p]
    assert len(solve(eqs, syms, check=False, simplify=False)) == 1


def test_sympyissue_5849():
    I1, I2, I3, I4, I5, I6 = symbols('I1:7')
    dI1, dI4, dQ2, dQ4, Q2, Q4 = symbols('dI1,dI4,dQ2,dQ4,Q2,Q4')

    e = (I1 - I2 - I3,
         I3 - I4 - I5,
         I4 + I5 - I6,
         -I1 + I2 + I6,
         -2*I1 - 2*I3 - 2*I5 - 3*I6 - dI1/2 + 12,
         -I4 + dQ4,
         -I2 + dQ2,
         2*I3 + 2*I5 + 3*I6 - Q2,
         I4 - 2*I5 + 2*Q4 + dI4)
    e = tuple(_.subs({I3: I6}) for _ in e)

    ans = [{dQ4: I3 - I5,
            dI1: -4*I2 - 8*I3 - 4*I5 - 6*I6 + 24,
            I4: I3 - I5,
            dQ2: I2,
            Q2: 2*I3 + 2*I5 + 3*I6,
            I1: I2 + I3,
            Q4: -I3/2 + 3*I5/2 - dI4/2}]
    ans = [{k: v.subs({I3: I6}) for k, v in ans[0].items()}]
    syms = I1, I4, Q2, Q4, dI1, dI4, dQ2, dQ4
    assert solve(e, *syms) == ans
    assert [_.subs(ans[0]) for _ in e] == [0]*9


def test_sympyissue_5901():
    f, g, h = map(Function, 'fgh')
    a = Symbol('a')
    D = Derivative(f(x), x)
    G = Derivative(g(a), a)
    assert solve(f(x) + f(x).diff(x), f(x)) == [{f(x): -D}]
    assert solve(f(x) - 3, f(x)) == [{f(x): 3}]
    assert solve(f(x) - 3*f(x).diff(x), f(x)) == [{f(x): 3*D}]
    assert solve([f(x) - 3*f(x).diff(x)], f(x)) == [{f(x): 3*D}]
    assert (solve([f(x) - 3*f(x).diff(x), f(x)**2 - y + 4], f(x), y) ==
            [{f(x): 3*D, y: 9*D**2 + 4}])
    assert solve(-f(a)**2*g(a)**2 + f(a)**2*h(a)**2 + g(a).diff(a),
                 h(a), g(a)) == [{g(a): -sqrt(h(a)**2 + G/f(a)**2)},
                                 {g(a): sqrt(h(a)**2 + G/f(a)**2)}]
    eqs = [f(x)**2 + g(x) - 2*f(x).diff(x), g(x)**2 - 4]
    assert solve(eqs, f(x), g(x)) == [{f(x): -sqrt(2)*sqrt(D - 1), g(x): 2},
                                      {f(x): sqrt(2)*sqrt(D - 1), g(x): 2},
                                      {f(x): -sqrt(2)*sqrt(D + 1), g(x): -2},
                                      {f(x): sqrt(2)*sqrt(D + 1), g(x): -2}]

    # the underlying problem was in solve_linear that was not masking off
    # anything but a Mul or Add; it now raises an error if it gets anything
    # but a symbol and solve handles the substitutions necessary so solve_linear
    # won't make this error
    pytest.raises(ValueError,
                  lambda: solve_linear(f(x) + f(x).diff(x), f(x)))
    assert solve_linear(f(x) + f(x).diff(x), x) == \
        (f(x) + Derivative(f(x), x), 1)
    assert solve_linear(f(x) + Integral(x, (x, y)), x) == \
        (f(x) + Integral(x, (x, y)), 1)
    assert solve_linear(f(x) + Integral(x, (x, y)) + x, x) == \
        (x + f(x) + Integral(x, (x, y)), 1)
    assert solve_linear(f(y) + Integral(x, (x, y)) + x, x) == \
        (x, -f(y) - Integral(x, (x, y)))
    assert solve_linear(x - f(x)/a + (f(x) - 1)/a, x) == \
        (x, 1/a)
    assert solve_linear(x + Derivative(2*x, x), x) == \
        (x, -2)
    assert solve_linear(x + Integral(x, y), x) == \
        (x + Integral(x, y), 1)
    assert solve_linear(x + Integral(x, y) - 2, x) == \
        (x + Integral(x, y) - 2, 1)

    assert solve(x + exp(x)**2, exp(x)) == [{exp(x): -sqrt(-x)},
                                            {exp(x): +sqrt(-x)}]


def test_sympyissue_5912():
    assert (solve(x**2 - x - 0.1, rational=True) ==
            [{x: Rational(1, 2) + sqrt(35)/10},
             {x: Rational(1, 2) - sqrt(35)/10}])
    ans = solve(x**2 - x - 0.1, rational=False)
    assert len(ans) == 2
    assert all(a[x].is_Number for a in ans)
    ans = solve(x**2 - x - 0.1)
    assert len(ans) == 2
    assert all(a[x].is_Number for a in ans)


def test_float_handling():
    def test(e1, e2):
        return len(e1.atoms(Float)) == len(e2.atoms(Float))
    assert solve(x - 0.5, rational=True)[0][x].is_Rational
    assert solve(x - 0.5, rational=False)[0][x].is_Float
    assert solve(x - Rational(1, 2), rational=False)[0][x].is_Rational
    assert solve(x - 0.5, rational=None)[0][x].is_Float
    assert solve(x - Rational(1, 2), rational=None)[0][x].is_Rational
    assert test(nfloat(1 + 2*x), 1.0 + 2.0*x)
    for contain in [list, tuple, set]:
        ans = nfloat(contain([1 + 2*x]))
        assert type(ans) is contain
        assert test(list(ans)[0], 1.0 + 2.0*x)
    k, v = list(nfloat({2*x: [1 + 2*x]}).items())[0]
    assert test(k, 2*x)
    assert test(v[0], 1.0 + 2.0*x)
    assert test(nfloat(cos(2*x)), cos(2.0*x))
    assert test(nfloat(3*x**2), 3.0*x**2)
    assert test(nfloat(3*x**2, exponent=True), 3.0*x**2.0)
    assert test(nfloat(exp(2*x)), exp(2.0*x))
    assert test(nfloat(x/3), x/3.0)
    assert test(nfloat(x**4 + 2*x + cos(Rational(1, 3)) + 1),
                x**4 + 2.0*x + 1.94495694631474)
    # don't call nfloat if there is no solution
    tot = 100 + c + z + t
    assert solve(((.7 + c)/tot - .6, (.2 + z)/tot - .3, t/tot - .1)) == []


def test_check_assumptions():
    x = symbols('x', positive=True)
    assert solve(x**2 - 1) == [{x: 1}]

    with pytest.warns(UserWarning) as warn:
        assert solve(x**2 - y, x, warn=True) == [{x: -sqrt(y)}, {x: sqrt(y)}]
    assert len(warn) == 1
    assert warn[0].message.args[0][:112] == """
                        Warning: assumptions concerning following
solution(s)                 can't be checked:"""


def test_sympyissue_6056():
    assert solve(tanh(x + 3)*tanh(x - 3) - 1) == []
    assert {simplify(w[x]) for w in solve((tanh(x + 3)*tanh(x - 3) + 1)**2)} == {
        -log(2)/2 + log(1 - I),
        -log(2)/2 + log(-1 - I),
        -log(2)/2 + log(1 + I),
        -log(2)/2 + log(-1 + I), }


def test_sympyissue_6060():
    absxm3 = Piecewise(
        (x - 3, 0 <= x - 3),
        (3 - x, 0 > x - 3)
    )
    y = Symbol('y', positive=True)
    assert solve(absxm3 - y, x) == [{x: -y + 3}, {x: y + 3}]


def test_checksol():
    pytest.raises(ValueError, lambda: checksol(x**4 - 1, 1))
    assert checksol(x*(x - y/x), {x: 1}, force=False) is False

    sol = {y: sqrt(x)}
    with pytest.warns(UserWarning) as warn:
        assert checksol(sqrt(y**2), sol, warn=True, force=False) is None
    assert len(warn) == 1
    assert warn[0].message.args[0] == f"""
\tWarning: could not verify solution {sol}."""

    eq = r - x**2 - y**2
    dict_var_soln = {y: - sqrt(r) / sqrt(tan(t)**2 + 1),
                     x: -sqrt(r)*tan(t)/sqrt(tan(t)**2 + 1)}
    assert checksol(eq, dict_var_soln) is True


def test_exclude():
    R, C, Ri, Vout, V1, Vminus, Vplus, s = \
        symbols('R, C, Ri, Vout, V1, Vminus, Vplus, s')
    Rf = symbols('Rf', positive=True)  # to eliminate Rf = 0 soln
    eqs = [C*V1*s + Vplus*(-2*C*s - 1/R),
           Vminus*(-1/Ri - 1/Rf) + Vout/Rf,
           C*Vplus*s + V1*(-C*s - 1/R) + Vout/R,
           -Vminus + Vplus]
    assert solve(eqs, Rf, Ri, V1, Vminus, Vout, Vplus) == [
        {
            Rf: Ri*(C*R*s + 1)**2/(C*R*s),
            Vminus: Vplus,
            V1: 2*Vplus + Vplus/(C*R*s),
            Vout: C*R*Vplus*s + 3*Vplus + Vplus/(C*R*s)},
        {
            Vplus: 0,
            Vminus: 0,
            V1: 0,
            Vout: 0},
    ]

    # TODO: Investingate why currently solution [0] is preferred over [1].
    assert solve(eqs, R, Rf, Ri, V1, Vminus, Vout) in [[{
        Vminus: Vplus,
        V1: Vout/2 + Vplus/2 + sqrt((Vout - 5*Vplus)*(Vout - Vplus))/2,
        R: (Vout - 3*Vplus - sqrt(Vout**2 - 6*Vout*Vplus + 5*Vplus**2))/(2*C*Vplus*s),
        Rf: Ri*(Vout - Vplus)/Vplus,
    }, {
        Vminus: Vplus,
        V1: Vout/2 + Vplus/2 - sqrt((Vout - 5*Vplus)*(Vout - Vplus))/2,
        R: (Vout - 3*Vplus + sqrt(Vout**2 - 6*Vout*Vplus + 5*Vplus**2))/(2*C*Vplus*s),
        Rf: Ri*(Vout - Vplus)/Vplus,
    }], [{
        Vminus: Vplus,
        Vout: (V1**2 - V1*Vplus - Vplus**2)/(V1 - 2*Vplus),
        Rf: Ri*(V1 - Vplus)**2/(Vplus*(V1 - 2*Vplus)),
        R: Vplus/(C*s*(V1 - 2*Vplus)),
    }]]


def test_high_order_roots():
    s = x**5 + 4*x**3 + 3*x**2 + Rational(7, 4)
    assert solve(s) == [{x: RootOf(s*4, i)} for i in range(5)]


def test_minsolve_linear_system():
    def count(dic):
        return len([x for x in dic.values() if x == 0])

    m = Matrix([[0, 0, 1, 1, 1, 0], [1, 1, 0, 1, 1, 0]])
    v = (a, t, x, y, z)
    assert count(minsolve_linear_system(m, *v, quick=True)) == 3
    assert count(minsolve_linear_system(m, *v)) == 3

    m = Matrix([[0, 1, 1, 1, 0], [1, 0, 1, 1, 0]])
    v = (a, x, y, z)
    assert count(minsolve_linear_system(m, *v, quick=True)) == 1
    assert count(minsolve_linear_system(m, *v)) == 2


def test_real_roots():
    # cf. issue sympy/sympy#6650
    x = Symbol('x', extended_real=True)
    assert len(solve(x**5 + x**3 + 1)) == 1


def test_sympyissue_6528():
    eqs = [
        327600995*x**2 - 37869137*x + 1809975124*y**2 - 9998905626,
        895613949*x**2 - 273830224*x*y + 530506983*y**2 - 10000000000]
    # two expressions encountered are > 1400 ops long so if this hangs
    # it is likely because simplification is being done
    assert len(solve(eqs, y, x, check=False)) == 4


def test_overdetermined():
    x = symbols('x', extended_real=True)
    eqs = [abs(4*x - 7) - 5, abs(3 - 8*x) - 1]
    assert solve(eqs, x) == [{x: Rational(1, 2)}]
    assert solve(eqs, x) == [{x: Rational(1, 2)}]
    assert solve(eqs, x, check=False) == [{x: Rational(1, 2)}, {x: 3}]


def test_sympyissue_6605():
    x = symbols('x')
    assert solve(4**(x/2) - 2**(x/3)) == [{x: 0}, {x: 3*I*pi/log(2)}]
    # while the first one passed, this one failed
    assert solve(5**(x/2) - 2**(x/3)) == [{x: 0}]
    x = symbols('x', real=True, nonzero=True)
    b = sqrt(6)*sqrt(log(2))/sqrt(log(5))
    assert [expand_log(s[x]) for s in solve(5**(x/2) - 2**(3/x))] == [-b, b]


def test_sympyissue_6644():
    eq = -sqrt((m - q)**2 + (-m/(2*q) + Rational(1, 2))**2) + sqrt((-m**2/2 - sqrt(
        4*m**4 - 4*m**2 + 8*m + 1)/4 - Rational(1, 4))**2 + (m**2/2 - m - sqrt(
            4*m**4 - 4*m**2 + 8*m + 1)/4 - Rational(1, 4))**2)
    sol = solve(eq, q, simplify=False, check=False)
    assert len(sol) == 5


def test_sympyissue_6752():
    assert solve([a**2 + a, a - b], [a, b]) == [{a: -1, b: -1},
                                                {a: 0, b: 0}]
    assert solve([a**2 + a*c, a - b], [a, b]) == [{a: 0, b: 0},
                                                  {a: -c, b: -c}]


def test_sympyissue_6792():
    assert (solve(x*(x - 1)**2*(x + 1)*(x**6 - x + 1)) ==
            [{x: -1}, {x: 0}, {x: 1}] +
            [{x: r} for r in (x**6 - x + 1).as_poly().all_roots()])


def test_sympyissues_6819_6820_6821():
    # issue sympy/sympy#6821
    x, y = symbols('x y', extended_real=True)
    assert solve(abs(x + 3) - 2*abs(x - 3)) == [{x: 1}, {x: 9}]
    assert solve(abs(x - 7) - 8) == [{x: -1}, {x: 15}]

    # issue sympy/sympy#7145
    assert solve(2*abs(x) - abs(x - 1)) == [{x: -1},
                                            {x: Rational(1, 3)}]

    i = symbols('i', imaginary=True)
    x = symbols('x')
    assert solve(abs(i) - 3) == [{i: -3*I}, {i: 3*I}]
    pytest.raises(NotImplementedError, lambda: solve(abs(x) - 3))

    w = symbols('w', integer=True)
    assert solve(2*x**w - 4*y**w, w) == solve((x/y)**w - 2, w)

    x, y = symbols('x y', extended_real=True)
    # issue sympy/sympy#2642
    assert solve(x*(1 + I)) == [{x: 0}]

    x, y = symbols('x y', imaginary=True)
    assert solve(x + y*I + 3 + 2*I) == [{x: -2*I, y: 3*I}]

    x = symbols('x', extended_real=True)
    assert solve(x + y + 3 + 2*I) == [{x: -3, y: -2*I}]


def test_sympyissue_6248():
    assert solve(2**x + 4**x) == [{x: I*pi/log(2)}]


def test_lambert_multivariate():
    for _ in range(7):
        assert _filtered_gens((x + 1/x + exp(x) + y).as_poly(), x) == {x, exp(x)}
        assert _filtered_gens((x + 1/x + exp(x)).as_poly(), x) == {exp(x), x}
        assert _filtered_gens((x + log(x) + 1/x + exp(x)).as_poly(),
                              x) == {exp(x), log(x), x}
        assert _filtered_gens((exp(I*x) - 1/x + log(x)/exp(I*x) + 2*x).as_poly(),
                              x) == {exp(I*x), x, log(x)}
    assert not _lambert(x, x)
    assert solve((x**2 - 2*x + 1).subs({x: log(x) + 3*x})) == [{x: LambertW(3*E)/3}]
    assert (solve((x**2 - 2*x + 1).subs({x: (log(x) + 3*x)**2 - 1})) ==
            [{x: LambertW(3*exp(-sqrt(2)))/3}, {x: LambertW(3*exp(sqrt(2)))/3}])
    assert (solve((x**2 - 2*x - 2).subs({x: log(x) + 3*x})) ==
            [{x: LambertW(3*exp(1 + sqrt(3)))/3},
             {x: LambertW(3*exp(-sqrt(3) + 1))/3}])
    assert solve(x*log(x) + 3*x + 1, x) == [{x: exp(-3 + LambertW(-exp(3)))}]
    eq = (x*exp(x) - 3).subs({x: x*exp(x)})
    assert solve(eq) == [{x: LambertW(3*exp(-LambertW(3)))}]
    # coverage test
    pytest.raises(NotImplementedError, lambda: solve(x - sin(x)*log(y - x), x))

    _13 = Rational(1, 3)
    _56 = Rational(5, 6)
    _53 = Rational(5, 3)
    assert (solve(3*log(a**(3*x + 5)) + a**(3*x + 5), x) ==
            [{x: log(exp(-LambertW(_13)/3)/a**_53)/log(a)},
             {x: log(exp(-LambertW(_13)/3)*(-1 - sqrt(3)*I)/(2*a**_53))/log(a)},
             {x: log(exp(-LambertW(_13)/3)*(-1 + sqrt(3)*I)/(2*a**_53))/log(a)}])
    p = symbols('p', positive=True)
    assert (solve(3*log(p**(3*x + 5)) + p**(3*x + 5), x) ==
            [{x: (-5*log(p) + log(LambertW(_13)) + log(3))/(3*log(p))},
             {x: log((-3**_13 - 3**_56*I)*LambertW(_13)**_13/(2*p**_53))/log(p)},
             {x: log((-3**_13 + 3**_56*I)*LambertW(_13)**_13/(2*p**_53))/log(p)}])

    # check collection
    assert (solve(3*log(a**(3*x + 5)) + b*log(a**(3*x + 5)) + a**(3*x + 5), x) ==
            [{x: log(exp(-LambertW(1/(b + 3))/3)/a**_53)/log(a)},
             {x: log(exp(-LambertW(1/(b + 3))/3)*(-1 - sqrt(3)*I)/(2*a**_53))/log(a)},
             {x: log(exp(-LambertW(1/(b + 3))/3)*(-1 + sqrt(3)*I)/(2*a**_53))/log(a)}])

    eq = 4*2**(2*p + 3) - 2*p - 3
    assert _solve_lambert(eq, p, _filtered_gens(eq.as_poly(), p)) == [
        -Rational(3, 2) - LambertW(-4*log(2))/(2*log(2))]

    # issue sympy/sympy#4271
    assert (solve((a/x + exp(x/2)).diff((x, 2)), x) ==
            [{x: 6*LambertW(root(-1, 3)*root(a, 3)/3)}])

    assert (solve((log(x) + x).subs({x: x**2 + 1})) ==
            [{x: -I*sqrt(-LambertW(1) + 1)}, {x: sqrt(-1 + LambertW(1))}])

    # these only give one of the solutions (see XFAIL below)
    assert solve(x**3 - 3**x, x) == [{x: -3/log(3)*LambertW(-log(3)/3)},
                                     {x: -3*LambertW(-log(3)/3, -1)/log(3)}]
    #     replacing 3 with 2 in the above solution gives 2
    assert solve(x**2 - 2**x, x) == [{x: 2}, {x: -2*LambertW(-log(2)/2, -1)/log(2)}]
    assert solve(-x**2 + 2**x, x) == [{x: 2}, {x: -2*LambertW(-log(2)/2, -1)/log(2)}]
    assert (solve(3**cos(x) - cos(x)**3) ==
            [{x: acos(-3*LambertW(-log(3)/3)/log(3))},
             {x: acos(-3*LambertW(-log(3)/3, -1)/log(3))}])

    # issue sympy/sympy#5673
    eq = -x + exp(exp(LambertW(log(x)))*LambertW(log(x)))
    assert checksol(eq, {x: 2}) is True


@pytest.mark.xfail
@pytest.mark.slow
def test_other_lambert():
    solve(3*sin(x) - x*sin(3), x)  # == [{x: 3}]
    # assert {s[x] for s in solve(3*log(x) - x*log(3))} == {3, -3*LambertW(-log(3)/3)/log(3)}
    # a = Rational(6, 5)
    # assert {s[x] for s in solve(x**a - a**x)} == {a, -a*LambertW(-log(a)/a)/log(a)}
    # assert ({s[x] for s in solve(3**cos(x) - cos(x)**3)} ==
    #         {acos(3), acos(-3*LambertW(-log(3)/3)/log(3))})
    # assert {s[x] for s in solve(x**2 - 2**x)} == {2, -2/log(2)*LambertW(log(2)/2)}


def test_rewrite_trig():
    assert solve(sin(x) + tan(x)) == [{x: 0}, {x: -pi}, {x: pi}, {x: 2*pi}]
    assert (solve(sin(x) + sec(x)) ==
            [{x: -2*atan(Rational(-1, 2) + sqrt(4 + (1 - sqrt(3)*I)**2)/2 + sqrt(3)*I/2)},
             {x: 2*atan(Rational(1, 2) - sqrt(3)*I/2 + sqrt(4 + (1 - sqrt(3)*I)**2)/2)},
             {x: 2*atan(Rational(1, 2) - sqrt(4 + (1 + sqrt(3)*I)**2)/2 + sqrt(3)*I/2)},
             {x: 2*atan(Rational(1, 2) + sqrt(4 + (1 + sqrt(3)*I)**2)/2 + sqrt(3)*I/2)}])
    assert solve(sinh(x) + tanh(x)) == [{x: 0}, {x: I*pi}]

    # issue sympy/sympy#6157
    assert solve(2*sin(x) - cos(x), x) == [{x: -2*atan(2 + sqrt(5))},
                                           {x: -2*atan(-sqrt(5) + 2)}]


def test_rewrite_trigh():
    assert solve(sinh(x) + sech(x)) == [{x: log(sqrt(-2 + sqrt(5))) + I*pi},
                                        {x: log(-I*sqrt(2 + sqrt(5)))},
                                        {x: log(I*sqrt(2 + sqrt(5)))},
                                        {x: log(sqrt(-2 + sqrt(5)))}]


def test_uselogcombine():
    eq = z - log(x) + log(y/(x*(-1 + y**2/x**2)))
    assert solve(eq, x) == [{x: sqrt(y*(-exp(z) + y))}, {x: -sqrt(-y*(exp(z) - y))}]
    assert (solve(log(x + 3) + log(1 + 3/x) - 3) in
            [[{x: -3 + sqrt(-12 + exp(3))*exp(Rational(3, 2))/2 + exp(3)/2},
              {x: -sqrt(-12 + exp(3))*exp(Rational(3, 2))/2 - 3 + exp(3)/2}],
             [{x: -3 + sqrt(-36 + (-exp(3) + 6)**2)/2 + exp(3)/2},
              {x: -3 - sqrt(-36 + (-exp(3) + 6)**2)/2 + exp(3)/2}]])
    assert solve(log(exp(2*x) + 1) + log(-tanh(x) + 1) - log(2)) == []


def test_atan2():
    assert solve(atan2(x, 2) - pi/3, x) == [{x: 2*sqrt(3)}]


def test_errorinverses():
    assert solve(erf(x) - y, x) == [{x: erfinv(y)}]
    assert solve(erfinv(x) - y, x) == [{x: erf(y)}]
    assert solve(erfc(x) - y, x) == [{x: erfcinv(y)}]
    assert solve(erfcinv(x) - y, x) == [{x: erfc(y)}]


def test_sympyissue_2725():
    R = Symbol('R')
    eq = sqrt(2)*R*sqrt(1/(R + 1)) + (R + 1)*(sqrt(2)*sqrt(1/(R + 1)) - 1)
    sol = solve(eq, R)
    assert sol == [{R: Rational(5, 3) + (-Rational(1, 2) -
                                         sqrt(3)*I/2)*cbrt(Rational(251, 27) +
                                                           sqrt(111)*I/9) +
                    40/(9*((-Rational(1, 2) -
                            sqrt(3)*I/2)*cbrt(Rational(251, 27) +
                                              sqrt(111)*I/9)))},
                   {R: Rational(5, 3) + 40/(9*cbrt(Rational(251, 27) +
                                                   sqrt(111)*I/9)) +
                    cbrt(Rational(251, 27) + sqrt(111)*I/9)}]


def test_piecewise():
    x = symbols('x')
    # if no symbol is given the piecewise detection must still work
    assert solve(Piecewise((x - 2, x > 2), (2 - x, True)) - 3) == [{x: -1}, {x: 5}]

    assert solve(abs(y)*x - 1, x) == [{x: 1/abs(y)}]

    # issue sympy/sympy#6989
    f = Function('f')
    assert (solve(Eq(-f(x), Piecewise((1, x > 0), (0, True))), f(x)) ==
            [{f(x): Piecewise((-1, x > 0), (0, True))}])

    # issue sympy/sympy#6060
    absxm3 = Piecewise(
        (x - 3, 0 <= x - 3),
        (3 - x, 0 > x - 3)
    )
    assert (solve(absxm3 - y, x) ==
            [{x: Piecewise((-y + 3, y > 0), (nan, True))},
             {x: Piecewise((y + 3, 0 <= y), (nan, True))}])

    # issue sympy/sympy#6821
    x = symbols('x', extended_real=True)
    assert solve([abs(x) - 2, arg(x) - pi], x) == [{x: -2}, {x: 2}]

    # issue sympy/sympy#8692
    assert (solve(Eq(abs(x + 1) + abs(x**2 - 7), 9), x) ==
            [{x: -Rational(1, 2) + sqrt(61)/2},
             {x: -sqrt(69)/2 + Rational(1, 2)}])


def test_real_imag_splitting():
    a, b = symbols('a b', extended_real=True)
    assert solve(sqrt(a**2 + b**2) - 3, a) == [{a: -sqrt(-b**2 + 9)},
                                               {a: +sqrt(-b**2 + 9)}]
    a, b = symbols('a b', imaginary=True)
    x = symbols('x')
    assert solve(sqrt(a**2 + b**2) - 3, a) == []

    assert solve([re(x) - 1, im(x) - 2], x) == [{re(x): 1,
                                                 x: 1 + 2*I,
                                                 im(x): 2}]
    # check for 'dict' handling of solution
    eq = sqrt(re(x)**2 + im(x)**2) - 3
    assert solve(eq) == solve(eq, x)

    x, y = symbols('x y', extended_real=True)
    assert solve(x + y*I + 3) == [{y: 0, x: -3}]
    assert solve([x + y*I + 3, y]) == [{x: -3, y: 0}]


def test_sympyissue_7110():
    y = -2*x**3 + 4*x**2 - 2*x + 5
    assert any(i[x].is_real for i in solve(y))


def test_sympyissue_7895():
    r = symbols('r', extended_real=True)
    assert solve(sqrt(r) - 2) == [{r: 4}]


def test_sympyissue_2777():
    # the equations represent two circles
    x, y = symbols('x y', extended_real=True)
    e1, e2 = sqrt(x**2 + y**2) - 10, sqrt(y**2 + (-x + 10)**2) - 3
    a, b = Rational(191, 20), 3*sqrt(391)/20
    ans = [{x: a, y: -b}, {x: a, y: b}]
    assert solve((e1, e2), (x, y)) == ans
    assert solve((e1, e2/(x - a)), (x, y)) == []
    # make the 2nd circle's radius be -3
    e2 += 6
    assert solve((e1, e2), (x, y)) == []
    assert solve((e1, e2), (x, y), check=False) == ans


def test_sympyissue_7322():
    number = 5.62527e-35
    assert solve(x - number, x)[0][x] == number


def test_sympyissue_8587():
    f = Piecewise((2*x**2, (0 < x) & (x < 1)), (2, True))
    assert solve(f - 1) == [{x: 1/sqrt(2)}]


def test_high_order_multivariate():
    assert len(solve(a*x**3 - x + 1, x)) == 3
    assert len(solve(a*x**4 - x + 1, x)) == 4
    assert (solve(a*x**5 - x + 1, x) ==
            [{x: r} for r in (a*x**5 - x + 1).as_poly(x).all_roots()])

    # result checking must always consider the denominator and RootOf
    # must be checked, too
    d = x**5 - x + 1
    assert solve(d*(1 + 1/d)) == [{x: r} for r in (d + 1).as_poly().all_roots()]
    d = x - 1
    assert solve(d*(2 + 1/d)) == [{x: Rational(1, 2)}]


def test_base_0_exp_0():
    assert solve(0**x - 1) == [{x: 0}]
    assert solve(0**(x - 2) - 1) == [{x: 2}]
    e = x*(Pow(Pow(x, 0, evaluate=False), -1, evaluate=False) - x)
    assert solve(e) == [{x: 0}, {x: 1}]


def test_sympyissue_8755():
    # This tests two things: that if full unrad is attempted and fails
    # the solution should still be found; also it tests the use of
    # keyword `composite`.
    assert len(solve(sqrt(y)*x + x**3 - 1, x)) == 3
    assert len(solve(-512*y**3 + 1344*cbrt(x + 2)*y**2 -
                     1176*(x + 2)**Rational(2, 3)*y - 169*x + 686, y, _unrad=False)) == 3


@pytest.mark.slow
def test_sympyissue_8828():
    x1 = 0
    y1 = -620
    r1 = 920
    x2 = 126
    y2 = 276
    x3 = 51
    y3 = 205
    r3 = 104
    v = x, y, z

    f1 = (x - x1)**2 + (y - y1)**2 - (r1 - z)**2
    f2 = (x2 - x)**2 + (y2 - y)**2 - z**2
    f3 = (x - x3)**2 + (y - y3)**2 - (r3 - z)**2
    F = f1, f2, f3

    g1 = sqrt((x - x1)**2 + (y - y1)**2) + z - r1
    g2 = f2
    g3 = sqrt((x - x3)**2 + (y - y3)**2) + z - r3
    G = g1, g2, g3

    A = solve(F, v)
    B = solve(G, v)

    p, q = ({tuple(i.evalf(2) for i in ordered(j)) for j in R} for R in [A, B])
    assert p == q


def test_sympyissue_10391():
    assert solve((2*x + 8)*exp(-6*x), x) == [{x: -4}]


def test_sympyissue_11538():
    eqs = (x - y**3 + 4, x + y + 4 + 4*E)
    assert solve(eqs, x, y) == [{x: -RootOf(y**3 + y + 4*E, i) - 4*E - 4,
                                 y: +RootOf(y**3 + y + 4*E, i)}
                                for i in range(3)]


def test_sympyissue_12180():
    e1, e2 = x - y*b, x*a - y
    assert (solve(e1, [x, y]) ==
            solve(e1, x, y) == [{x: y*b}])
    assert (solve(e2, [x, y]) ==
            solve(e2, x, y) == [{x: y/a}])


def test_issue_427():
    assert solve([1 + y, x - y], x) == []
    assert solve([x - y, y - 3], x) == []


def test_sympyissue_14645():
    eq = x*y - x - y
    ans = [{x: y/(y - 1)}]
    assert solve([eq], [x, y]) == ans
    assert solve([eq, eq], [x, y]) == ans


def test_sympyissue_14721():
    assert solve([h, -1 + (-k + 1)**2/b**2 + (-h - 1)**2/a**2,
                  -1 + (-k + 1)**2/b**2 + (-h + 1)**2/a**2, k + 2],
                 h, k, a, b) == [{a: -sqrt(b**2/(b**2 - 9)), h: 0, k: -2},
                                 {a: sqrt(b**2/(b**2 - 9)), h: 0, k: -2}]


def test_sympyissue_14791():
    assert solve(exp(log(5)*x) - exp(log(2)*x), x) == [{x: 0}]


def test_sympyissue_14779():
    x = symbols('x', real=True)
    assert solve(sqrt(x**4 - 130*x**2 + 1089) +
                 sqrt(x**4 - 130*x**2 + 3969) - 96*abs(x)/x) == [{x: sqrt(130)}]


@pytest.mark.slow
def test_issue_709():
    eqs = [sqrt(x) + y + 2, root(y, 3)*x - 1]
    a1 = symbols('a1')
    r5, r6 = (a1**7 + 4*a1**4 + 4*a1 - 1).as_poly().all_roots()[-2:]
    res = [{x: 4 + 4*r5**3 + r5**6, y: r5**3},
           {x: 4 + r6**6 + 4*r6**3, y: r6**3}]
    assert solve(eqs) == res


@pytest.mark.slow
def test_unrad1():
    eq = sqrt(x) + sqrt(x + 1) + sqrt(1 - sqrt(x))

    assert {s[x] for s in solve(eq, check=False)} == {0, Rational(9, 16)}
    assert solve(eq) == []
    # but this one really does have those solutions
    assert ({s[x] for s in solve(sqrt(x) - sqrt(x + 1) +
                                 sqrt(1 - sqrt(x)))} ==
            {0, Rational(9, 16)})

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

    # for coverage
    assert solve(sqrt(x) + root(x, 3) - 2) == [{x: 1}]
    pytest.raises(NotImplementedError, lambda:
                  solve(sqrt(x) + root(x, 3) + root(x + 1, 5) - 2))
    # fails through a different code path
    pytest.raises(NotImplementedError, lambda: solve(-sqrt(2) + cosh(x)/x))

    # the simplify flag should be reset to False for unrad results;
    # if it's not then this next test will take a long time
    assert solve(root(x, 3) + root(x, 5) - 2) == [{x: 1}]
    eq = sqrt(x) + sqrt(x + 1) + sqrt(1 - x) - 6*sqrt(5)/5
    ans = [{x: Rational(4, 5)},
           {x: Rational(-1484, 375) + 172564/(140625*cbrt(114*sqrt(12657)/78125 +
                                                          Rational(12459439, 52734375))) +
               4*cbrt(114*sqrt(12657)/78125 +
                      Rational(12459439, 52734375))}]
    assert solve(eq) == ans

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

    ans = solve(sqrt(x) + sqrt(x + 1) -
                sqrt(1 - x) - sqrt(2 + x))
    assert len(ans) == 1
    assert NS(ans[0][x])[:4] == '0.73'
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
    assert solve(eq, x, check=False) != []  # not other code errors


@pytest.mark.xfail
def test_unrad1_fail():
    assert solve(sqrt(x + root(x, 3)) + root(x - y, 5), y) != []


@pytest.mark.slow
def test_unrad_slow():
    # this has roots with multiplicity > 1; there should be no
    # repeats in roots obtained, however
    eq = (sqrt(1 + sqrt(1 - 4*x**2)) - x*(1 + sqrt(1 + 2*sqrt(1 - 4*x**2))))
    assert solve(eq) == [{x: Rational(1, 2)}]


def test_unrad2():
    assert solve(root(x**3 - 3*x**2, 3) + 1 - x) == []
    assert solve(root(x + 1, 3) + root(x**2 - 2, 5) + 1) == []


def test_sympyissue_20610():
    assert solve([x + y, sqrt(2)], [x, y]) == []


def test_sympyissue_21167():
    assert solve(cbrt(x - 1) + cbrt(x) + cbrt(x + 1)) == []


def test_sympyissue_21766():
    assert solve([z + y/x, - (z + y/x)]) == [{x: -y/z}, {y: 0, z: 0}]


def test_sympyissue_21852():
    assert solve(2*x + sqrt(2*x**2) - 21) == [{x: 21 - 21*sqrt(2)/2}]


def test_sympyissue_21882():
    eqs = [-a*k + 4*a/3 + b + 2*c/9 + 5*d/6 + 5*f/6, d/2 - f*k + 4*f/3,
           -d*k + d + f/6, 13*a/18 + 13*b/18 + 13*c/18, a + b/2 - c*k + 20*c/9,
           a/6 - b*k + b + c/18, a + 5*b/3 + c/3, 4*a/3 + 2*b/3 + 2*c, -g]
    sols = [{a: 0, b: 0, c: 0, d: 0, f: 0, g: 0},
            {a: 0, b: 0, c: 0, d: -f, g: 0, k: Rational(5, 6)},
            {a: -2*c, b: c, d: 0, f: 0, g: 0, k: Rational(13, 18)}]
    assert solve(eqs) == sols


def test_sympyissue_21890():
    assert solve([4*x**3*y**4 - 2*y,
                  4*x**4*y**3 - 2*x]) == [{x: root(4, 3)/(2*y)},
                                          {x: (-root(4, 3)/4 -
                                               root(4, 3)*sqrt(3)*I/4)/y},
                                          {x: (-root(4, 3)/4 +
                                               root(4, 3)*sqrt(3)*I/4)/y},
                                          {x: 0, y: 0}]


def test_sympyissue_21905():
    f = 0.07*x*y + 10.0/y + 30.0/x
    eqs = [f.diff(x), f.diff(y)]
    assert (solve(eqs) ==
            [{x: Float('10.873803730028921', dps=15),
              y: Float('3.624601243342974', dps=15)},
             {x: Float('-5.4369018650144607', dps=15) - Float('9.4169902659710321', dps=15)*I,
              y: Float('-1.812300621671487', dps=15) - Float('3.1389967553236771', dps=15)*I},
             {x: Float('-5.4369018650144607', dps=15) + Float('9.4169902659710321', dps=15)*I,
              y: Float('-1.812300621671487', dps=15) + Float('3.1389967553236771', dps=15)*I}])


def test_sympyissue_21984():
    ka = 10**5
    C0 = 10**-10
    kw = 10**-14  # that is constant
    H, OH, HA, A = symbols('H OH HA A')

    eqs = [(H*A/HA) - ka, H*OH - kw, A + OH - H, HA + A - C0]

    res = [{H: Float('-100000.0000000001', dps=15),
            OH: Float('-9.9999999999999901e-20', dps=15),
            HA: Float('100000.0000000002', dps=15),
            A: Float('-100000.0000000001', dps=15)},
           {H: Float('1.0005409447144861e-7', dps=15),
            OH: Float('9.9954094471448725e-8', dps=15),
            HA: Float('1.0005409447134852e-22', dps=15),
            A: Float('9.9999999999899941e-11', dps=15)},
           {H: Float('-9.7518697933531939e-8', dps=15),
            OH: Float('-9.7618697933532039e-8', dps=15),
            HA: Float('-9.7518697933631848e-23', dps=15),
            A: Float('1.0000000000009752e-10', dps=15)}]
    assert solve(eqs, [H, OH, HA, A]) == res


def test_sympyissue_22051():
    lamda, mu = symbols('lamda mu')
    eqs = [(x*y*(lamda + 30000*mu) + 6*log(2*y))/x,
           lamda*x + 30000*mu*x + 6*log(x)/y, x*y - 10, 30000*x*y - 1550]
    assert solve(eqs) == []
    eqs2 = [y*(exp(x) + 1), x + exp(x)]
    assert solve(eqs2) == [{x: -LambertW(1), y: 0}]


def test_sympyissue_22058():
    assert solve(-sqrt(t)*x**2 + 2*x + sqrt(t),
                 x) == [{x: -sqrt(1 + 1/t) + 1/sqrt(t)},
                        {x: sqrt(1 + 1/t) + 1/sqrt(t)}]


def test_sympyissue_22248():
    g = 9.81
    y0 = 100
    v0 = 55
    m = 80
    c = 15
    y = y0 + (m/c)*(v0 + ((m*g)/c))*(1 - exp(- (c/m)*x)) - ((m*g)/c)*x

    assert solve(y) == [{x: Float('-1.4164130909148258', dps=15)},
                        {x: Float('11.61083847106101', dps=15)}]


def test_sympyissue_22837():
    assert solve(Eq(0, (4 - 4*x + x**2)/(4*a**2)), x) == [{x: 2}]


def test_issue_1195():
    assert solve(0.9**x - 0.35) == [{x: 9.964094404284511}]


@pytest.mark.slow
@pytest.mark.timeout(2100)
def test_sympyissue_23637():
    p, q, s, t = 4*y + 3*x + 1, -2*y - x + 1, -y - x + 1, 2*y + 2*x + 6
    rea1 = 1.0*10**(-4)*(x + 2)*p**3 - 1.3*q*s*t**2
    rea2 = 1.0*10**(-4)*(y + 1)*p**4 - 2.99*q**2*s*t**2
    sol = [{x: -11.0, y: 8.0}, {x: -8.474434848793706, y: 9.4724262173802511},
           {x: -4.0297080219649217, y: 0.99590437023147171},
           {x: -3.814341662315357, y: 0.84581641948719222},
           {x: -3.0, y: 2.0}, {x: -2.0091163509386085, y: -0.98835947935546686},
           {x: -2.0, y: -1.0}, {x: 1.3235849477227735, y: -0.32419504840343982},
           {x: 3.0040159362898193, y: -1.0015924793400075}, {x: 5.0, y: -4.0}]
    assert solve([rea1, rea2], x, y) == sol


def test_sympyissue_23855():
    assert solve([x - 1], x, x) == [{x: 1}]
