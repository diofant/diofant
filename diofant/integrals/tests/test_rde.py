"""Most of these tests come from the examples in Bronstein's book."""

import pytest

from diofant import I, Poly, Rational, oo, symbols
from diofant.abc import n, t, x, z
from diofant.integrals.rde import (bound_degree, cancel_exp, cancel_primitive,
                                   no_cancel_equal, normal_denom, order_at,
                                   order_at_oo, rischDE, solve_poly_rde, spde,
                                   special_denom, weak_normalizer)
from diofant.integrals.risch import (DifferentialExtension,
                                     NonElementaryIntegralException)


__all__ = ()

t0, t1, t2, k = symbols('t:3 k')


def test_order_at():
    a = Poly(t**4, t)
    b = Poly((t**2 + 1)**3*t, t)
    c = Poly((t**2 + 1)**6*t, t)
    d = Poly((t**2 + 1)**10*t**10, t)
    e = Poly((t**2 + 1)**100*t**37, t)
    p1 = Poly(t, t)
    p2 = Poly(1 + t**2, t)
    assert order_at(a, p1, t) == 4
    assert order_at(b, p1, t) == 1
    assert order_at(c, p1, t) == 1
    assert order_at(d, p1, t) == 10
    assert order_at(e, p1, t) == 37
    assert order_at(a, p2, t) == 0
    assert order_at(b, p2, t) == 3
    assert order_at(c, p2, t) == 6
    assert order_at(d, p1, t) == 10
    assert order_at(e, p2, t) == 100
    assert order_at(Poly(0, t), Poly(t, t), t) == oo
    assert order_at_oo(Poly(t**2 - 1, t), Poly(t + 1), t) == \
        order_at_oo(Poly(t - 1, t), Poly(1, t), t) == -1
    assert order_at_oo(Poly(0, t), Poly(1, t), t) == oo


def test_weak_normalizer():
    a = Poly((1 + x)*t**5 + 4*t**4 + (-1 - 3*x)*t**3 - 4*t**2 + (-2 + 2*x)*t, t)
    d = Poly(t**4 - 3*t**2 + 2, t)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    r = weak_normalizer(a, d, DE, z)
    assert r == (Poly(t**5 - t**4 - 4*t**3 + 4*t**2 + 4*t - 4, t),
                 (Poly((1 + x)*t**2 + x*t, t), Poly(t + 1, t)))
    assert weak_normalizer(r[1][0], r[1][1], DE) == (Poly(1, t), r[1])
    r = weak_normalizer(Poly(1 + t**2), Poly(t**2 - 1, t), DE, z)
    assert r == (Poly(t**4 - 2*t**2 + 1, t), (Poly(-3*t**2 + 1, t), Poly(t**2 - 1, t)))
    assert weak_normalizer(r[1][0], r[1][1], DE, z) == (Poly(1, t), r[1])
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2)]})
    r = weak_normalizer(Poly(1 + t**2), Poly(t, t), DE, z)
    assert r == (Poly(t, t), (Poly(0, t), Poly(1, t)))
    assert weak_normalizer(r[1][0], r[1][1], DE, z) == (Poly(1, t), r[1])


def test_normal_denom():
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    pytest.raises(NonElementaryIntegralException, lambda: normal_denom(Poly(1, x), Poly(1, x),
                                                                       Poly(1, x), Poly(x, x), DE))
    fa, fd = Poly(t**2 + 1, t), Poly(1, t)
    ga, gd = Poly(1, t), Poly(t**2, t)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert normal_denom(fa, fd, ga, gd, DE) == \
        (Poly(t, t), (Poly(t**3 - t**2 + t - 1, t), Poly(1, t)), (Poly(1, t),
                                                                  Poly(1, t)), Poly(t, t))


def test_special_denom():
    # TODO: add more tests here
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert special_denom(Poly(1, t), Poly(t**2, t), Poly(1, t), Poly(t**2 - 1, t),
                         Poly(t, t), DE) == \
        (Poly(1, t), Poly(t**2 - 1, t), Poly(t**2 - 1, t), Poly(t, t))
#    assert special_denom(Poly(1, t), Poly(2*x, t), Poly((1 + 2*x)*t, t), DE) == 1

    # issue sympy/sympy#3940
    # Note, this isn't a very good test, because the denominator is just 1,
    # but at least it tests the exp cancellation case
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(-2*x*t0, t0),
                                                Poly(I*k*t1, t1)]})
    DE.decrement_level()
    assert special_denom(Poly(1, t0), Poly(I*k, t0), Poly(1, t0), Poly(t0, t0),
                         Poly(1, t0), DE) == \
        (Poly(1, t0), Poly(I*k, t0), Poly(t0, t0), Poly(1, t0))


@pytest.mark.xfail
def test_bound_degree_fail():
    # Primitive
    DE = DifferentialExtension(extension={'D': [Poly(1, x),
                                                Poly(t0/x**2, t0), Poly(1/x, t)]})
    assert bound_degree(Poly(t**2, t), Poly(-(1/x**2*t**2 + 1/x), t),
                        Poly((2*x - 1)*t**4 + (t0 + x)/x*t**3 - (t0 + 4*x**2)/2*x*t**2 + x*t,
                             t), DE) == 3


def test_bound_degree():
    # Base
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert bound_degree(Poly(1, x), Poly(-2*x, x), Poly(1, x), DE) == 0

    # Primitive (see above test_bound_degree_fail)
    # TODO: Add test for when the degree bound becomes larger after limited_integrate
    # TODO: Add test for db == da - 1 case

    # Exp
    # TODO: Add tests
    # TODO: Add test for when the degree becomes larger after parametric_log_deriv()

    # Nonlinear
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert bound_degree(Poly(t, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t), DE) == 0


def test_spde():
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    pytest.raises(NonElementaryIntegralException, lambda: spde(Poly(t, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t), 0, DE))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert spde(Poly(t**2 + x*t*2 + x**2, t), Poly(t**2/x**2 + (2/x - 1)*t, t),
                Poly(t**2/x**2 + (2/x - 1)*t, t), 0, DE) == \
        (Poly(0, t), Poly(0, t), 0, Poly(0, t), Poly(1, t))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t0/x**2, t0), Poly(1/x, t)]})
    assert spde(Poly(t**2, t), Poly(-t**2/x**2 - 1/x, t),
                Poly((2*x - 1)*t**4 + (t0 + x)/x*t**3 - (t0 + 4*x**2)/(2*x)*t**2 + x*t, t), 3, DE) == \
        (Poly(0, t), Poly(0, t), 0, Poly(0, t),
         Poly(t0*t**2/2 + x**2*t**2 - x**2*t, t))
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert spde(Poly(x**2 + x + 1, x), Poly(-2*x - 1, x), Poly(x**5/2 +
                                                               3*x**4/4 + x**3 - x**2 + 1, x), 4, DE) == \
        (Poly(0, x), Poly(x/2 - Rational(1, 4), x), 2, Poly(x**2 + x + 1, x), Poly(5*x/4, x))
    assert spde(Poly(x**2 + x + 1, x), Poly(-2*x - 1, x), Poly(x**5/2 +
                                                               3*x**4/4 + x**3 - x**2 + 1, x), n, DE) == \
        (Poly(0, x), Poly(x/2 - Rational(1, 4), x), -2 + n, Poly(x**2 + x + 1, x), Poly(5*x/4, x))
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1, t)]})
    pytest.raises(NonElementaryIntegralException, lambda: spde(Poly((t - 1)*(t**2 + 1)**2, t), Poly((t - 1)*(t**2 + 1), t), Poly(1, t), 0, DE))
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert spde(Poly(x**2 - x, x), Poly(1, x), Poly(9*x**4 - 10*x**3 + 2*x**2, x), 4, DE) == (Poly(0, x), Poly(0, x), 0, Poly(0, x), Poly(3*x**3 - 2*x**2, x))
    assert spde(Poly(x**2 - x, x), Poly(x**2 - 5*x + 3, x), Poly(x**7 - x**6 - 2*x**4 + 3*x**3 - x**2, x), 5, DE) == \
        (Poly(1, x), Poly(x + 1, x), 1, Poly(x**4 - x**3, x), Poly(x**3 - x**2, x))


def test_solve_poly_rde_no_cancel():
    # deg(b) large
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1 + t**2, t)]})
    assert solve_poly_rde(Poly(t**2 + 1, t), Poly(t**3 + (x + 1)*t**2 + t + x + 2, t),
                          oo, DE) == Poly(t + x, t)
    # deg(b) small
    DE = DifferentialExtension(extension={'D': [Poly(1, x)]})
    assert solve_poly_rde(Poly(0, x), Poly(x/2 - Rational(1, 4), x), oo, DE) == \
        Poly(x**2/4 - x/4, x)
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert solve_poly_rde(Poly(2, t), Poly(t**2 + 2*t + 3, t), 1, DE) == \
        Poly(t + 1, t, x)
    # deg(b) == deg(D) - 1
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t**2 + 1, t)]})
    assert no_cancel_equal(Poly(1 - t, t),
                           Poly(t**3 + t**2 - 2*x*t - 2*x, t), oo, DE) == \
        (Poly(t**2, t), 1, Poly((-2 - 2*x)*t - 2*x, t))


def test_solve_poly_rde_cancel():
    # exp
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    assert cancel_exp(Poly(2*x, t), Poly(2*x, t), 0, DE) == \
        Poly(1, t)
    assert cancel_exp(Poly(2*x, t), Poly((1 + 2*x)*t, t), 1, DE) == \
        Poly(t, t)
    # TODO: Add more exp tests, including tests that require is_deriv_in_field()

    # primitive
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(1/x, t)]})

    # If the DecrementLevel context manager is working correctly, this shouldn't
    # cause any problems with the further tests.
    pytest.raises(NonElementaryIntegralException, lambda: cancel_primitive(Poly(1, t), Poly(t, t), oo, DE))

    assert cancel_primitive(Poly(1, t), Poly(t + 1/x, t), 2, DE) == \
        Poly(t, t)
    assert cancel_primitive(Poly(4*x, t), Poly(4*x*t**2 + 2*t/x, t), 3, DE) == \
        Poly(t**2, t)

    # TODO: Add more primitive tests, including tests that require is_deriv_in_field()


def test_rischDE():
    # TODO: Add more tests for rischDE, including ones from the text
    DE = DifferentialExtension(extension={'D': [Poly(1, x), Poly(t, t)]})
    DE.decrement_level()
    assert rischDE(Poly(-2*x, x), Poly(1, x), Poly(1 - 2*x - 2*x**2, x),
                   Poly(1, x), DE) == \
        (Poly(x + 1, x), Poly(1, x))
