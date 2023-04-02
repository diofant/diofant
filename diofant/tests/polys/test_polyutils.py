"""Tests for useful utilities for higher level polynomial classes."""

import pytest

from diofant import (ZZ, GeneratorsNeededError, I, Integer, Integral, Mul,
                     PolynomialError, Rational, cos, erf, exp, factor,
                     integrate, pi, sin, sqrt, symbols)
from diofant.abc import p, q, t, x, y, z
from diofant.polys.polyoptions import build_options
from diofant.polys.polyutils import (_nsort, _parallel_dict_from_expr,
                                     _sort_factors, _sort_gens, _unify_gens)


__all__ = ()

A, B = symbols('A,B', commutative=False)


def test__nsort():
    # issue sympy/sympy#6137
    r = ([Rational(3, 2) + sqrt(Rational(-14, 3) - 2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3) -
                                4/sqrt(Rational(-7, 3) + 61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) +
                                       2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) -
                                61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)))/2 -
          sqrt(-7/3 + 61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) +
               2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3))/2,
          Rational(3, 2) - sqrt(Rational(-7, 3) + 61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) +
                                2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3))/2 -
          sqrt(Rational(-14, 3) - 2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3) -
               4/sqrt(Rational(-7, 3) + 61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) +
                      2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) -
               61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)))/2, Rational(3, 2) +
          sqrt(Rational(-14, 3) - 2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3) +
               4/sqrt(Rational(-7, 3) + 61/(18*(Rational(-415, 216) +
                                                13*I/12)**Rational(1, 3)) + 2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) -
               61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)))/2 +
          sqrt(Rational(-7, 3) + 61/(18*(Rational(-415, 216)
                                         + 13*I/12)**Rational(1, 3)) + 2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3))/2,
          Rational(3, 2) + sqrt(Rational(-7, 3) + 61/(18*(Rational(-415, 216) +
                                                          13*I/12)**Rational(1, 3)) + 2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3))/2 -
          sqrt(Rational(-14, 3) - 2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3) +
               4/sqrt(Rational(-7, 3) + 61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) +
                      2*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)) -
               61/(18*(Rational(-415, 216) + 13*I/12)**Rational(1, 3)))/2])
    ans = [r[1], r[0], r[-1], r[-2]]
    assert _nsort(r) == ans
    assert len(_nsort(r, separated=True)[0]) == 0
    b, c, a = exp(-1000), exp(-999), exp(-1001)
    assert _nsort((b, c, a)) == [a, b, c]
    d = symbols('d', extended_real=True)
    assert _nsort((d,)) == [d]
    assert _nsort((d,), separated=True) == [[d], []]
    c = symbols('c', complex=True, real=False)
    assert _nsort((c,)) == [c]
    assert _nsort((c,), separated=True) == [[], [c]]
    assert _nsort((I, Rational(1)), separated=True) == ([Rational(1)], [I])


def test__sort_gens():
    assert not _sort_gens([])

    assert _sort_gens([x]) == (x,)
    assert _sort_gens([p]) == (p,)
    assert _sort_gens([q]) == (q,)

    assert _sort_gens([x, p]) == (x, p)
    assert _sort_gens([p, x]) == (x, p)
    assert _sort_gens([q, p]) == (p, q)

    assert _sort_gens([q, p, x]) == (x, p, q)

    assert _sort_gens([x, p, q], wrt=x) == (x, p, q)
    assert _sort_gens([x, p, q], wrt=p) == (p, x, q)
    assert _sort_gens([x, p, q], wrt=q) == (q, x, p)

    assert _sort_gens([x, p, q], wrt='x') == (x, p, q)
    assert _sort_gens([x, p, q], wrt='p') == (p, x, q)
    assert _sort_gens([x, p, q], wrt='q') == (q, x, p)

    assert _sort_gens([x, p, q], wrt='x,q') == (x, q, p)
    assert _sort_gens([x, p, q], wrt='q,x') == (q, x, p)
    assert _sort_gens([x, p, q], wrt='p,q') == (p, q, x)
    assert _sort_gens([x, p, q], wrt='q,p') == (q, p, x)

    assert _sort_gens([x, p, q], wrt='x, q') == (x, q, p)
    assert _sort_gens([x, p, q], wrt='q, x') == (q, x, p)
    assert _sort_gens([x, p, q], wrt='p, q') == (p, q, x)
    assert _sort_gens([x, p, q], wrt='q, p') == (q, p, x)

    assert _sort_gens([x, p, q], wrt=[x, 'q']) == (x, q, p)
    assert _sort_gens([x, p, q], wrt=[q, 'x']) == (q, x, p)
    assert _sort_gens([x, p, q], wrt=[p, 'q']) == (p, q, x)
    assert _sort_gens([x, p, q], wrt=[q, 'p']) == (q, p, x)

    assert _sort_gens([x, p, q], wrt=['x', 'q']) == (x, q, p)
    assert _sort_gens([x, p, q], wrt=['q', 'x']) == (q, x, p)
    assert _sort_gens([x, p, q], wrt=['p', 'q']) == (p, q, x)
    assert _sort_gens([x, p, q], wrt=['q', 'p']) == (q, p, x)

    assert _sort_gens([x, p, q], sort='x > p > q') == (x, p, q)
    assert _sort_gens([x, p, q], sort='p > x > q') == (p, x, q)
    assert _sort_gens([x, p, q], sort='p > q > x') == (p, q, x)

    assert _sort_gens([x, p, q], wrt='x', sort='q > p') == (x, q, p)
    assert _sort_gens([x, p, q], wrt='p', sort='q > x') == (p, q, x)
    assert _sort_gens([x, p, q], wrt='q', sort='p > x') == (q, p, x)

    X = symbols('x0:3 x10:13 x20:23')

    assert _sort_gens(X) == X


def test__unify_gens():
    assert not _unify_gens([], [])

    assert _unify_gens([x], [x]) == (x,)
    assert _unify_gens([y], [y]) == (y,)

    assert _unify_gens([x, y], [x]) == (x, y)
    assert _unify_gens([x], [x, y]) == (x, y)

    assert _unify_gens([x, y], [x, y]) == (x, y)
    assert _unify_gens([y, x], [y, x]) == (y, x)

    assert _unify_gens([x], [y]) == (x, y)
    assert _unify_gens([y], [x]) == (y, x)

    assert _unify_gens([x], [y, x]) == (y, x)
    assert _unify_gens([y, x], [x]) == (y, x)

    assert _unify_gens([x, y, z], [x, y, z]) == (x, y, z)
    assert _unify_gens([z, y, x], [x, y, z]) == (z, y, x)
    assert _unify_gens([x, y, z], [z, y, x]) == (x, y, z)
    assert _unify_gens([z, y, x], [z, y, x]) == (z, y, x)

    assert _unify_gens([x, y, z], [t, x, p, q, z]) == (t, x, y, p, q, z)


def test__sort_factors():
    assert _sort_factors([], multiple=True) == []
    assert _sort_factors([], multiple=False) == []

    F = [[1, 2, 3], [1, 2], [1]]
    G = [[1], [1, 2], [1, 2, 3]]

    assert _sort_factors(F, multiple=False) == G

    F = [[1, 2], [1, 2, 3], [1, 2], [1]]
    G = [[1], [1, 2], [1, 2], [1, 2, 3]]

    assert _sort_factors(F, multiple=False) == G

    F = [[2, 2], [1, 2, 3], [1, 2], [1]]
    G = [[1], [1, 2], [2, 2], [1, 2, 3]]

    assert _sort_factors(F, multiple=False) == G

    F = [([1, 2, 3], 1), ([1, 2], 1), ([1], 1)]
    G = [([1], 1), ([1, 2], 1), ([1, 2, 3], 1)]

    assert _sort_factors(F, multiple=True) == G

    F = [([1, 2], 1), ([1, 2, 3], 1), ([1, 2], 1), ([1], 1)]
    G = [([1], 1), ([1, 2], 1), ([1, 2], 1), ([1, 2, 3], 1)]

    assert _sort_factors(F, multiple=True) == G

    F = [([2, 2], 1), ([1, 2, 3], 1), ([1, 2], 1), ([1], 1)]
    G = [([1], 1), ([1, 2], 1), ([2, 2], 1), ([1, 2, 3], 1)]

    assert _sort_factors(F, multiple=True) == G

    F = [([2, 2], 1), ([1, 2, 3], 1), ([1, 2], 2), ([1], 1)]
    G = [([1], 1), ([2, 2], 1), ([1, 2], 2), ([1, 2, 3], 1)]

    assert _sort_factors(F, multiple=True) == G


def test__dict_from_expr_if_gens():
    opt = build_options([x], {})
    assert _parallel_dict_from_expr([Integer(17)],
                                    opt) == ([{(0,): 17}], opt)
    assert _parallel_dict_from_expr([Integer(-17)],
                                    opt) == ([{(0,): -17}], opt)
    assert _parallel_dict_from_expr([17*x], opt) == ([{(1,): 17}], opt)
    assert _parallel_dict_from_expr([17*x**7], opt) == ([{(7,): 17}], opt)
    assert _parallel_dict_from_expr([x + 2*y + 3*z],
                                    opt) == ([{(1,): 1,
                                               (0,): 2*y + 3*z}], opt)
    assert _parallel_dict_from_expr([x*y + 2*x*z + 3*y*z],
                                    opt) == ([{(1,): y + 2*z,
                                               (0,): 3*y*z}], opt)
    assert _parallel_dict_from_expr([2**y*x],
                                    opt) == ([{(1,): 2**y}], opt)

    opt = build_options([x, y], {})
    assert _parallel_dict_from_expr([Integer(17)],
                                    opt) == ([{(0, 0): 17}], opt)
    assert _parallel_dict_from_expr([Integer(-17)],
                                    opt) == ([{(0, 0): -17}], opt)
    assert _parallel_dict_from_expr([17*x],
                                    opt) == ([{(1, 0): 17}], opt)
    assert _parallel_dict_from_expr([17*x**7*y],
                                    opt) == ([{(7, 1): 17}], opt)
    assert _parallel_dict_from_expr([x + 2*y + 3*z],
                                    opt) == ([{(1, 0): 1, (0, 1): 2,
                                               (0, 0): 3*z}], opt)
    assert _parallel_dict_from_expr([x*y + 2*x*z + 3*y*z],
                                    opt) == ([{(1, 1): 1, (1, 0): 2*z,
                                               (0, 1): 3*z}], opt)
    pytest.raises(PolynomialError,
                  lambda: _parallel_dict_from_expr([2**y*x], opt))

    opt = build_options([x, y, z], {})
    assert _parallel_dict_from_expr([Integer(17)],
                                    opt) == ([{(0, 0, 0): 17}], opt)
    assert _parallel_dict_from_expr([Integer(-17)],
                                    opt) == ([{(0, 0, 0): -17}], opt)
    assert _parallel_dict_from_expr([17*x],
                                    opt) == ([{(1, 0, 0): 17}], opt)
    assert _parallel_dict_from_expr([17*x**7*y*z**12],
                                    opt) == ([{(7, 1, 12): 17}], opt)
    assert _parallel_dict_from_expr([x + 2*y + 3*z],
                                    opt) == ([{(1, 0, 0): 1,
                                               (0, 1, 0): 2,
                                               (0, 0, 1): 3}], opt)
    assert _parallel_dict_from_expr([x*y + 2*x*z + 3*y*z],
                                    opt) == ([{(1, 1, 0): 1,
                                               (1, 0, 1): 2,
                                               (0, 1, 1): 3}], opt)

    opt = build_options([], {})
    assert _parallel_dict_from_expr([Integral(x, (x, 1, 2)) +
                                     x], opt) == ([{(0, 1): 1, (1, 0): 1}],
                                                  build_options([x, Integral(x, (x, 1, 2))], {}))


def test__dict_from_expr_no_gens():
    opt = build_options([], {})
    pytest.raises(GeneratorsNeededError,
                  lambda: _parallel_dict_from_expr([Integer(17)], opt))

    assert _parallel_dict_from_expr([x], opt) == ([{(1,): 1}], build_options([x], {}))
    assert _parallel_dict_from_expr([y], opt) == ([{(1,): 1}], build_options([y], {}))

    assert _parallel_dict_from_expr([x*y], opt) == ([{(1, 1): 1}], build_options([x, y], {}))
    assert _parallel_dict_from_expr([x + y], opt) == ([{(1, 0): 1, (0, 1): 1}],
                                                      build_options([x, y], {}))

    assert _parallel_dict_from_expr([sqrt(2)], opt) == ([{(1,): 1}], build_options([sqrt(2)], {}))

    f = cos(x)*sin(x) + cos(x)*sin(y) + cos(y)*sin(x) + cos(y)*sin(y)

    assert _parallel_dict_from_expr([f], opt) == ([{(0, 1, 0, 1): 1, (0, 1, 1, 0): 1,
                                                    (1, 0, 0, 1): 1,
                                                    (1, 0, 1, 0): 1}],
                                                  build_options([cos(x), cos(y), sin(x), sin(y)], {}))

    opt = build_options([], {'greedy': False})
    pytest.raises(GeneratorsNeededError,
                  lambda: _parallel_dict_from_expr([sqrt(2)], opt))

    opt = build_options([], {'domain': ZZ.inject(x)})
    assert _parallel_dict_from_expr([x*y], opt) == ([{(1,): x}], build_options([y], {'domain': ZZ.inject(x)}))
    opt = build_options([], {'domain': ZZ.inject(y)})
    assert _parallel_dict_from_expr([x*y], opt) == ([{(1,): y}], build_options([x], {'domain': ZZ.inject(y)}))

    opt = build_options([], {'extension': None})
    assert _parallel_dict_from_expr([3*sqrt(2)*pi*x*y],
                                    opt) == ([{(1, 1, 1, 1): 3}],
                                             build_options([x, y, pi, sqrt(2)], {'extension': None}))
    opt = build_options([], {'extension': True})
    assert _parallel_dict_from_expr([3*sqrt(2)*pi*x*y],
                                    opt) == ([{(1, 1, 1): 3*sqrt(2)}],
                                             build_options([x, y, pi], {'extension': True}))


def test__parallel_dict_from_expr_if_gens():
    opt = build_options([x], {})
    assert _parallel_dict_from_expr([x + 2*y + 3*z, Integer(7)], opt) == \
        ([{(1,): 1, (0,): 2*y + 3*z}, {(0,): 7}], opt)
    assert _parallel_dict_from_expr([Mul(x, x**2, evaluate=False)], opt) == \
        ([{(3,): 1}], opt)

    pytest.raises(PolynomialError, lambda: _parallel_dict_from_expr([A*x], opt))


def test__parallel_dict_from_expr_no_gens():
    opt = build_options([], {})
    assert _parallel_dict_from_expr([Mul(x, x**2, evaluate=False)], opt) == \
        ([{(3,): 1}], build_options([x], {}))
    assert _parallel_dict_from_expr([x*y, Integer(3)], opt) == \
        ([{(1, 1): 1}, {(0, 0): 3}], build_options([x, y], {}))
    assert _parallel_dict_from_expr([x*y, 2*z, Integer(3)], opt) == \
        ([{(1, 1, 0): 1}, {(0, 0, 1): 2}, {(0, 0, 0): 3}], build_options([x, y, z], {}))


def test__parallel_dict_from_expr():
    opt = build_options([], {})
    assert _parallel_dict_from_expr([x - 1, x**2 - 2], opt) == ([{(0,): -1, (1,): 1},
                                                                {(0,): -2,
                                                                 (2,): 1}], build_options([x], {}))
    assert _parallel_dict_from_expr([x - 1], opt) == ([{(0,): -1, (1,): 1}], build_options([x], {}))

    pytest.raises(PolynomialError, lambda: _parallel_dict_from_expr([A*B - B*A], opt))


def test_sympyissue_7383():
    x, z, R, a = symbols('x z R a')
    r = sqrt(x**2 + z**2)
    u = erf(a*r/sqrt(2))/r
    Ec = u.diff(z, z).subs({x: sqrt(R*R - z*z)})
    assert integrate(Ec, (z, -R, R)).simplify() == \
        -2*sqrt(2)*R*a**3*exp(-R**2*a**2/2)/(3*sqrt(pi))


def test_sympyissue_10161():
    x = symbols('x', real=True)
    h = (2*x*(-2*x + abs(x))*(x**2 - 1)/abs(x**2 - 1)
         + (x/abs(x) - 2)*abs(x**2 - 1))
    assert (h - factor(h)).simplify() == 0
