import pytest

from diofant import (Derivative, Dummy, Float, I, O, Piecewise, Rational,
                     RisingFactorial, Sum, Tuple, besseli, cos, exp, exp_polar,
                     expand_func, factorial, false, gamma, hyper, limit, log,
                     meijerg, oo, pi, polar_lift, sqrt, symbols)
from diofant.abc import a, b, c, d, k, l, s, x, z
from diofant.functions.special.hyper import (HyperRep, HyperRep_asin1,
                                             HyperRep_asin2, HyperRep_atanh,
                                             HyperRep_cosasin, HyperRep_log1,
                                             HyperRep_log2, HyperRep_power1,
                                             HyperRep_power2, HyperRep_sinasin,
                                             HyperRep_sqrts1, HyperRep_sqrts2)
from diofant.utilities.randtest import random_complex_number as randcplx
from diofant.utilities.randtest import verify_derivative_numerically as td
from diofant.utilities.randtest import verify_numerically as tn


__all__ = ()


def test_TupleParametersBase():
    # test that our implementation of the chain rule works
    p = hyper((), (), z**2)
    assert p.diff(z) == p*2*z


def test_hyper():
    pytest.raises(TypeError, lambda: hyper(1, 2, z))

    assert hyper((1, 2), (1,), z) == hyper(Tuple(1, 2), Tuple(1), z)

    h = hyper((1, 2), (3, 4, 5), z)
    assert h.ap == Tuple(1, 2)
    assert h.bq == Tuple(3, 4, 5)
    assert h.argument == z
    assert h.is_commutative is True

    # just a few checks to make sure that all arguments go where they should
    assert tn(hyper(Tuple(), Tuple(), z), exp(z), z)
    assert tn(z*hyper((1, 1), Tuple(2), -z), log(1 + z), z)

    # differentiation
    h = hyper(
        (randcplx(), randcplx(), randcplx()), (randcplx(), randcplx()), z)
    assert td(h, z)

    a1, a2, b1, b2, b3 = symbols('a1:3, b1:4')
    assert hyper((a1, a2), (b1, b2, b3), z).diff(z) == \
        a1*a2/(b1*b2*b3) * hyper((a1 + 1, a2 + 1), (b1 + 1, b2 + 1, b3 + 1), z)

    # differentiation wrt parameters is not supported
    assert hyper([z], [], z).diff(z) == Derivative(hyper([z], [], z), z)

    # hyper is unbranched wrt parameters
    assert hyper([polar_lift(z)], [polar_lift(k)], polar_lift(x)) == \
        hyper([z], [k], polar_lift(x))

    assert hyper((1, 2, 3), [3, 4], 1).is_number
    assert not hyper((1, 2, 3), [3, x], 1).is_number


def test_expand_func():
    # evaluation at 1 of Gauss' hypergeometric function:
    a1, b1, c1 = randcplx(), randcplx(), randcplx() + 5
    assert expand_func(hyper([a, b], [c], 1)) == \
        gamma(c)*gamma(-a - b + c)/(gamma(-a + c)*gamma(-b + c))
    assert abs(expand_func(hyper([a1, b1], [c1], 1))
               - hyper([a1, b1], [c1], 1)).evalf(strict=False) < 1e-10

    # hyperexpand wrapper for hyper:
    assert expand_func(hyper([], [], z)) == exp(z)
    assert expand_func(hyper([1, 2, 3], [], z)) == hyper([1, 2, 3], [], z)
    assert expand_func(meijerg([[1, 1], []], [[1], [0]], z)) == log(z + 1)
    assert expand_func(meijerg([[1, 1], []], [[], []], z)) == \
        meijerg([[1, 1], []], [[], []], z)


def replace_dummy(expr, sym):
    dum = expr.atoms(Dummy)
    if not dum:
        return expr
    assert len(dum) == 1
    return expr.xreplace({dum.pop(): sym})


def test_hyper_rewrite_sum():
    _k = Dummy('k')
    assert replace_dummy(hyper((1, 2), (1, 3), x).rewrite(Sum), _k) == \
        Sum(x**_k / factorial(_k) * RisingFactorial(2, _k) /
            RisingFactorial(3, _k), (_k, 0, oo))

    assert hyper((1, 2, 3), (-1, 3), z).rewrite(Sum) == \
        hyper((1, 2, 3), (-1, 3), z)


def test_radius_of_convergence():
    assert hyper((1, 2), [3], z).radius_of_convergence == 1
    assert hyper((1, 2), [3, 4], z).radius_of_convergence == oo
    assert hyper((1, 2, 3), [4], z).radius_of_convergence == 0
    assert hyper((0, 1, 2), [4], z).radius_of_convergence == oo
    assert hyper((-1, 1, 2), [-4], z).radius_of_convergence == 0
    assert hyper((-1, -2, 2), [-1], z).radius_of_convergence == oo
    assert hyper((-1, 2), [-1, -2], z).radius_of_convergence == 0
    assert hyper([-1, 1, 3], [-2, 2], z).radius_of_convergence == 1
    assert hyper([-1, 1], [-2, 2], z).radius_of_convergence == oo
    assert hyper([-1, 1, 3], [-2], z).radius_of_convergence == 0
    assert hyper((-1, 2, 3, 4), [], z).radius_of_convergence == oo
    assert hyper((-4, 2), [-3], z).radius_of_convergence == 0

    assert hyper([1, 1], [3], 1).convergence_statement
    assert hyper([1, 1], [2], 1).convergence_statement == false
    assert hyper([1, 1], [2], -1).convergence_statement
    assert hyper([1, 1], [1], -1).convergence_statement == false


def test_meijer():
    pytest.raises(TypeError, lambda: meijerg(1, z))
    pytest.raises(TypeError, lambda: meijerg(((1,), (2,)), (3,), (4,), z))
    pytest.raises(TypeError, lambda: meijerg((1, 2, 3), (4, 5), z))

    assert meijerg(((1, 2), (3,)), ((4,), (5,)), z) == \
        meijerg(Tuple(1, 2), Tuple(3), Tuple(4), Tuple(5), z)

    g = meijerg((1, 2), (3, 4, 5), (6, 7, 8, 9), (10, 11, 12, 13, 14), z)
    assert g.an == Tuple(1, 2)
    assert g.ap == Tuple(1, 2, 3, 4, 5)
    assert g.aother == Tuple(3, 4, 5)
    assert g.bm == Tuple(6, 7, 8, 9)
    assert g.bq == Tuple(6, 7, 8, 9, 10, 11, 12, 13, 14)
    assert g.bother == Tuple(10, 11, 12, 13, 14)
    assert g.argument == z
    assert g.nu == 75
    assert g.delta == -1
    assert g.is_commutative is True

    assert meijerg([1, 2], [3], [4], [5], z).delta == Rational(1, 2)

    # just a few checks to make sure that all arguments go where they should
    assert tn(meijerg(Tuple(), Tuple(), Tuple(0), Tuple(), -z), exp(z), z)
    assert tn(sqrt(pi)*meijerg(Tuple(), Tuple(),
                               Tuple(0), Tuple(Rational(1, 2)), z**2/4), cos(z), z)
    assert tn(meijerg(Tuple(1, 1), Tuple(), Tuple(1), Tuple(0), z),
              log(1 + z), z)

    # test exceptions
    pytest.raises(ValueError, lambda: meijerg(((3, 1), (2,)),
                                              ((oo,), (2, 0)), x))
    pytest.raises(ValueError, lambda: meijerg(((3, 1), (2,)),
                                              ((1,), (2, 0)), x))

    # differentiation
    g = meijerg((randcplx(),), (randcplx() + 2*I,), Tuple(),
                (randcplx(), randcplx()), z)
    assert td(g, z)

    g = meijerg(Tuple(), (randcplx(),), Tuple(),
                (randcplx(), randcplx()), z)
    assert td(g, z)

    g = meijerg(Tuple(), Tuple(), Tuple(randcplx()),
                Tuple(randcplx(), randcplx()), z)
    assert td(g, z)

    a1, a2, b1, b2, c1, c2, d1, d2 = symbols('a1:3, b1:3, c1:3, d1:3')
    assert meijerg((a1, a2), (b1, b2), (c1, c2), (d1, d2), z).diff(z) == \
        (meijerg((a1 - 1, a2), (b1, b2), (c1, c2), (d1, d2), z)
         + (a1 - 1)*meijerg((a1, a2), (b1, b2), (c1, c2), (d1, d2), z))/z

    assert meijerg([z, z], [], [], [], z).diff(z) == \
        Derivative(meijerg([z, z], [], [], [], z), z)

    # meijerg is unbranched wrt parameters
    assert meijerg([polar_lift(a1)], [polar_lift(a2)], [polar_lift(b1)],
                   [polar_lift(b2)], polar_lift(z)) == meijerg([a1], [a2],
                                                               [b1], [b2],
                                                               polar_lift(z))

    # integrand
    assert meijerg([a], [b], [c], [d], z).integrand(s) == \
        z**s*gamma(c - s)*gamma(-a + s + 1)/(gamma(b - s)*gamma(-d + s + 1))

    assert meijerg([[], []], [[Rational(1, 2)], [0]], 1).is_number
    assert not meijerg([[], []], [[x], [0]], 1).is_number


def test_meijerg_derivative():
    assert meijerg([], [1, 1], [0, 0, x], [], z).diff(x) == \
        log(z)*meijerg([], [1, 1], [0, 0, x], [], z) \
        + 2*meijerg([], [1, 1, 1], [0, 0, x, 0], [], z)

    y = randcplx()
    a = 5  # mpmath chokes with non-real numbers, and Mod1 with floats
    assert td(meijerg([x], [], [], [], y), x)
    assert td(meijerg([x**2], [], [], [], y), x)
    assert td(meijerg([], [x], [], [], y), x)
    assert td(meijerg([], [], [x], [], y), x)
    assert td(meijerg([], [], [], [x], y), x)
    assert td(meijerg([x], [a], [a + 1], [], y), x)
    assert td(meijerg([x], [a + 1], [a], [], y), x)
    assert td(meijerg([x, a], [], [], [a + 1], y), x)
    assert td(meijerg([x, a + 1], [], [], [a], y), x)

    b = Rational(3, 2)
    assert td(meijerg([a + 2], [b], [b - 3, x], [a], y), x)
    assert td(meijerg([x], [2, b], [1, b + 1], [], y), x)


def test_meijerg_period():
    assert meijerg([], [1], [0], [], x).get_period() == 2*pi
    assert meijerg([1], [], [], [0], x).get_period() == 2*pi
    assert meijerg([], [], [0], [], x).get_period() == 2*pi  # exp(x)
    assert meijerg(
        [], [], [0], [Rational(1, 2)], x).get_period() == 2*pi  # cos(sqrt(x))
    assert meijerg(
        [], [], [Rational(1, 2)], [0], x).get_period() == 4*pi  # sin(sqrt(x))
    assert meijerg([1, 1], [], [1], [0], x).get_period() == oo  # log(1 + x)


def test_hyper_unpolarify():
    a = exp_polar(2*pi*I)*x
    b = x
    assert hyper([], [], a).argument == b
    assert hyper([0], [], a).argument == a
    assert hyper([0], [0], a).argument == b
    assert hyper([0, 1], [0], a).argument == a


def test_hyperrep():
    # First test the base class works.
    a, b, c, d, z = symbols('a b c d z')

    class MyRep(HyperRep):
        @classmethod
        def _expr_small(cls, x):
            return a

        @classmethod
        def _expr_small_minus(cls, x):
            return b

        @classmethod
        def _expr_big(cls, x, n):
            return c*n

        @classmethod
        def _expr_big_minus(cls, x, n):
            return d*n
    assert MyRep(z).rewrite('nonrep') == Piecewise((0, abs(z) > 1), (a, True))
    assert MyRep(exp_polar(I*pi)*z).rewrite('nonrep') == \
        Piecewise((0, abs(z) > 1), (b, True))
    assert MyRep(exp_polar(2*I*pi)*z).rewrite('nonrep') == \
        Piecewise((c, abs(z) > 1), (a, True))
    assert MyRep(exp_polar(3*I*pi)*z).rewrite('nonrep') == \
        Piecewise((d, abs(z) > 1), (b, True))
    assert MyRep(exp_polar(4*I*pi)*z).rewrite('nonrep') == \
        Piecewise((2*c, abs(z) > 1), (a, True))
    assert MyRep(exp_polar(5*I*pi)*z).rewrite('nonrep') == \
        Piecewise((2*d, abs(z) > 1), (b, True))
    assert MyRep(z).rewrite('nonrepsmall') == a
    assert MyRep(exp_polar(I*pi)*z).rewrite('nonrepsmall') == b

    def t(func, hyp, z):
        """Test that func is a valid representation of hyp."""
        # First test that func agrees with hyp for small z
        if not tn(func.rewrite('nonrepsmall'), hyp, z,
                  a=Rational(-1, 2), b=Rational(-1, 2), c=Rational(1, 2), d=Rational(1, 2)):
            return False
        # Next check that the two small representations agree.
        if not tn(
            func.rewrite('nonrepsmall').subs({z: exp_polar(I*pi)*z}).replace(exp_polar, exp),
            func.subs({z: exp_polar(I*pi)*z}).rewrite('nonrepsmall'),
                z, a=Rational(-1, 2), b=Rational(-1, 2), c=Rational(1, 2), d=Rational(1, 2)):
            return False
        # Next check continuity along exp_polar(I*pi)*t
        expr = func.subs({z: exp_polar(I*pi)*z}).rewrite('nonrep')
        if abs(expr.subs({z: 1 + 1e-15}) - expr.subs({z: 1 - 1e-15})).evalf(strict=False) > 1e-10:
            return False
        # Finally check continuity of the big reps.

        def dosubs(func, a, b):
            rv = func.subs({z: exp_polar(a)*z}).rewrite('nonrep')
            return rv.subs({z: exp_polar(b)*z}).replace(exp_polar, exp)
        for n in [0, 1, 2, 3, 4, -1, -2, -3, -4]:
            expr1 = dosubs(func, 2*I*pi*n, I*pi/2)
            expr2 = dosubs(func, 2*I*pi*n + I*pi, -I*pi/2)
            if not tn(expr1, expr2, z):
                return False
            expr1 = dosubs(func, 2*I*pi*(n + 1), -I*pi/2)
            expr2 = dosubs(func, 2*I*pi*n + I*pi, I*pi/2)
            if not tn(expr1, expr2, z):
                return False
        return True

    # Now test the various representatives.
    a = Rational(1, 3)
    assert t(HyperRep_atanh(z), hyper([Rational(1, 2), 1], [Rational(3, 2)], z), z)
    assert t(HyperRep_power1(a, z), hyper([-a], [], z), z)
    assert t(HyperRep_power2(a, z), hyper([a, a - Rational(1, 2)], [2*a], z), z)
    assert t(HyperRep_log1(z), -z*hyper([1, 1], [2], z), z)
    assert t(HyperRep_asin1(z), hyper([Rational(1, 2), Rational(1, 2)], [Rational(3, 2)], z), z)
    assert t(HyperRep_asin2(z), hyper([1, 1], [Rational(3, 2)], z), z)
    assert t(HyperRep_sqrts1(a, z), hyper([-a, Rational(1, 2) - a], [Rational(1, 2)], z), z)
    assert t(HyperRep_sqrts2(a, z),
             -2*z/(2*a + 1)*hyper([-a - Rational(1, 2), -a], [Rational(1, 2)], z).diff(z), z)
    assert t(HyperRep_log2(z), -z/4*hyper([Rational(3, 2), 1, 1], [2, 2], z), z)
    assert t(HyperRep_cosasin(a, z), hyper([-a, a], [Rational(1, 2)], z), z)
    assert t(HyperRep_sinasin(a, z), 2*a*z*hyper([1 - a, 1 + a], [Rational(3, 2)], z), z)


def test_meijerg_eval():
    a = randcplx()
    arg = x*exp_polar(k*pi*I)
    expr1 = pi*meijerg([[], [(a + 1)/2]], [[a/2], [-a/2, (a + 1)/2]], arg**2/4)
    expr2 = besseli(a, arg)

    # Test that the two expressions agree for all arguments.
    for x_ in [0.5, 1.5]:
        for k_ in [0.0, 0.1, 0.3, 0.5, 0.8, 1, 5.751, 15.3]:
            assert abs((expr1 - expr2).evalf(subs={x: x_, k: k_}, strict=False)) < 1e-10
            assert abs((expr1 - expr2).evalf(subs={x: x_, k: -k_}, strict=False)) < 1e-10

    # Test continuity independently
    eps = 1e-13
    expr2 = expr1.subs({k: l})
    for x_ in [0.5, 1.5]:
        for k_ in [0.5, Rational(1, 3), 0.25, 0.75, Rational(2, 3), 1.0, 1.5]:
            assert abs((expr1 - expr2).evalf(
                       subs={x: x_, k: k_ + eps, l: k_ - eps})) < 1e-10
            assert abs((expr1 - expr2).evalf(
                       subs={x: x_, k: -k_ + eps, l: -k_ - eps})) < 1e-10

    expr = (meijerg(((0.5,), ()), ((0.5, 0, 0.5), ()), exp_polar(-I*pi)/4)
            + meijerg(((0.5,), ()), ((0.5, 0, 0.5), ()), exp_polar(I*pi)/4)) \
        / (2*sqrt(pi))
    assert (expr - pi/exp(1)).evalf(chop=True) == 0


def test_limits():
    k, x = symbols('k, x')
    assert hyper((1,), (Rational(4, 3), Rational(5, 3)), k**2).series(k) == \
        hyper((1,), (Rational(4, 3), Rational(5, 3)), 0) + \
        9*k**2*hyper((2,), (Rational(7, 3), Rational(8, 3)), 0)/20 + \
        81*k**4*hyper((3,), (Rational(10, 3), Rational(11, 3)), 0)/1120 + \
        O(k**6)  # issue sympy/sympy#6350
    assert limit(meijerg((), (), (1,), (0,), -x), x, 0) == \
        meijerg(((), ()), ((1,), (0,)), 0)  # issue sympy/sympy#6052


def test_evalf():
    assert hyper((-1, 1), (-1,), 1).evalf() == Float('2.0')

    e = meijerg([], [], [], [], (exp_polar(-I*pi))*cos(exp_polar(-I*pi)))
    assert e.evalf() == e
