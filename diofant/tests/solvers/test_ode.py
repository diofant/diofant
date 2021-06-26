import pytest

from diofant import (Derivative, Dummy, E, Ei, Eq, Function, I, Integer,
                     Integral, LambertW, Matrix, Mul, O, Piecewise, Rational,
                     RootOf, Subs, Symbol, acos, acosh, asin, asinh, atan,
                     cbrt, cos, diff, dsolve, erf, erfi, exp, log, pi, root,
                     simplify, sin, sinh, sqrt, sstr, symbols, tan, variations)
from diofant.abc import A
from diofant.solvers.deutils import ode_order
from diofant.solvers.ode import (_lie_group_remove, _linear_coeff_match,
                                 _undetermined_coefficients_match, checkinfsol,
                                 checkodesol, checksysodesol, classify_ode,
                                 classify_sysode, constant_renumber,
                                 constantsimp, get_numbered_constants,
                                 homogeneous_order, infinitesimals,
                                 ode_nth_linear_euler_eq_homogeneous,
                                 ode_sol_simplicity, odesimp, solve_init)


__all__ = ()

C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10 = symbols('C0:11')
x, y, z, t = symbols('x y z t', real=True)
f, g, h = symbols('f g h', cls=Function)

# Note: the tests below may fail (but still be correct) if ODE solver,
# the integral engine, solve(), or even simplify() changes. Also, in
# differently formatted solutions, the arbitrary constants might not be
# equal.  Using specific hints in tests can help to avoid this.

# Tests of order higher than 1 should run the solutions through
# constant_renumber because it will normalize it (constant_renumber causes
# dsolve() to return different results on different machines)


def test_linear_2eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    k, l, m, n = symbols('k, l, m, n', integer=True)
    t, a = symbols('t, a')
    x0, y0 = symbols('x0, y0', cls=Function)
    eq1 = [Eq(diff(x(t), t), 9*y(t)), Eq(diff(y(t), t), 12*x(t))]
    sol1 = [Eq(x(t), sqrt(3)*exp(6*sqrt(3)*t)*C2/2 - sqrt(3)*exp(-6*sqrt(3)*t)*C1/2),
            Eq(y(t), exp(6*sqrt(3)*t)*C2 + exp(-6*sqrt(3)*t)*C1)]
    assert dsolve(eq1) == sol1

    eq2 = [Eq(diff(x(t), t), 2*x(t) + 4*y(t)),
           Eq(diff(y(t), t), 12*x(t) + 41*y(t))]
    sol2 = [Eq(x(t), -4*exp(t*(-sqrt(1713) + 43)/2)*C1/(Rational(-39, 2) + sqrt(1713)/2) - 4*exp(t*(sqrt(1713) + 43)/2)*C2/(-sqrt(1713)/2 - Rational(39, 2))),
            Eq(y(t), exp(t*(-sqrt(1713) + 43)/2)*C1 + exp(t*(sqrt(1713) + 43)/2)*C2)]
    assert dsolve(eq2) == sol2

    eq3 = [Eq(diff(x(t), t), x(t) + y(t)), Eq(diff(y(t), t), -2*x(t) + 2*y(t))]
    sol3 = [Eq(x(t), -exp(t*(3 - sqrt(7)*I)/2)*C1/(Rational(-1, 2) + sqrt(7)*I/2) - exp(t*(3 + sqrt(7)*I)/2)*C2/(Rational(-1, 2) - sqrt(7)*I/2)),
            Eq(y(t), exp(t*(3 - sqrt(7)*I)/2)*C1 + exp(t*(3 + sqrt(7)*I)/2)*C2)]
    assert dsolve(eq3) == sol3

    eq4 = [Eq(diff(x(t), t), x(t) + y(t) + 9),
           Eq(diff(y(t), t), 2*x(t) + 5*y(t) + 23)]
    sol4 = [Eq(x(t), -exp(t*(-sqrt(6) + 3))*C1/(-2 + sqrt(6)) - exp(t*(sqrt(6) + 3))*C2/(-sqrt(6) - 2) + (-69 + 32*sqrt(6))/((-18 + 6*sqrt(6))*(-2 + sqrt(6))) + (-32*sqrt(6) - 69)/((-18 - 6*sqrt(6))*(-sqrt(6) - 2))),
            Eq(y(t), exp(t*(-sqrt(6) + 3))*C1 + exp(t*(sqrt(6) + 3))*C2 - (-32*sqrt(6) - 69)/(-18 - 6*sqrt(6)) - (-69 + 32*sqrt(6))/(-18 + 6*sqrt(6)))]
    assert dsolve(eq4) == sol4

    eq4_1 = [Eq(diff(x(t), t), x(t) + y(t) + 9),
             Eq(42*diff(y(t), t), 2*x(t) + 5*y(t) + 23)]
    sol4_1 = [Eq(x(t), -exp(t*(-sqrt(1705) + 47)/84)*C1/(Rational(37, 84) + sqrt(1705)/84) - exp(t*(sqrt(1705) + 47)/84)*C2/(-sqrt(1705)/84 + Rational(37, 84)) + (-3294060 - 55524*sqrt(1705))/((-6731340 - 143220*sqrt(1705))*(-sqrt(1705)/84 + Rational(37, 84))) + (-3294060 + 55524*sqrt(1705))/((-6731340 + 143220*sqrt(1705))*(Rational(37, 84) + sqrt(1705)/84))),
              Eq(y(t), exp(t*(-sqrt(1705) + 47)/84)*C1 + exp(t*(sqrt(1705) + 47)/84)*C2 - (-3294060 + 55524*sqrt(1705))/(-6731340 + 143220*sqrt(1705)) - (-3294060 - 55524*sqrt(1705))/(-6731340 - 143220*sqrt(1705)))]
    assert dsolve(eq4_1) == sol4_1

    eq5 = [Eq(diff(x(t), t), x(t) + y(t) + 81),
           Eq(diff(y(t), t), -2*x(t) + y(t) + 23)]
    sol5 = [Eq(x(t), sqrt(2)*exp(t*(1 - sqrt(2)*I))*I*C1/2 - sqrt(2)*exp(t*(1 + sqrt(2)*I))*I*C2/2 - sqrt(2)*I*(185 - 58*sqrt(2)*I)/12 + sqrt(2)*I*(185 + 58*sqrt(2)*I)/12),
            Eq(y(t), exp(t*(1 - sqrt(2)*I))*C1 + exp(t*(1 + sqrt(2)*I))*C2 - Rational(185, 3))]
    assert dsolve(eq5) == sol5

    eq6 = (Eq(diff(x(t), t), 5*t*x(t) + 2*y(t)), Eq(diff(y(t), t), 2*x(t) + 5*t*y(t)))
    sol6 = [Eq(x(t), (C1*exp(Integral(2, t)) + C2*exp(-Integral(2, t)))*exp(Integral(5*t, t))),
            Eq(y(t), (C1*exp(Integral(2, t)) - C2*exp(-Integral(2, t)))*exp(Integral(5*t, t)))]
    assert dsolve(eq6) == sol6

    eq7 = (Eq(diff(x(t), t), 5*t*x(t) + t**2*y(t)), Eq(diff(y(t), t), -t**2*x(t) + 5*t*y(t)))
    sol7 = [Eq(x(t), (C1*cos(Integral(t**2, t)) + C2*sin(Integral(t**2, t)))*exp(Integral(5*t, t))),
            Eq(y(t), (-C1*sin(Integral(t**2, t)) + C2*cos(Integral(t**2, t)))*exp(Integral(5*t, t)))]
    assert dsolve(eq7) == sol7

    eq8 = [Eq(diff(x(t), t), 5*t*x(t) + t**2*y(t)),
           Eq(diff(y(t), t), -t**2*x(t) + (5*t + 9*t**2)*y(t))]
    sol8 = [Eq(x(t), E**Integral(5*t, t)*(-exp((-sqrt(77) + 9)*Integral(t**2, t)/2)*C1/(Rational(-9, 2) + sqrt(77)/2) - exp((sqrt(77) + 9)*Integral(t**2, t)/2)*C2/(Rational(-9, 2) - sqrt(77)/2))),
            Eq(y(t), E**Integral(5*t, t)*(exp((-sqrt(77) + 9)*Integral(t**2, t)/2)*C1 + exp((sqrt(77) + 9)*Integral(t**2, t)/2)*C2))]
    assert dsolve(eq8) == sol8

    eq10 = (Eq(diff(x(t), t), 5*t*x(t) + t**2*y(t)), Eq(diff(y(t), t), (1-t**2)*x(t) + (5*t+9*t**2)*y(t)))
    sol10 = [Eq(x(t), C1*x0(t) + C2*x0(t)*Integral(t**2*exp(Integral(5*t, t))*exp(Integral(9*t**2 + 5*t, t))/x0(t)**2, t)),
             Eq(y(t), C1*y0(t) + C2*(y0(t)*Integral(t**2*exp(Integral(5*t, t))*exp(Integral(9*t**2 + 5*t, t))/x0(t)**2, t) +
                                     exp(Integral(5*t, t))*exp(Integral(9*t**2 + 5*t, t))/x0(t)))]
    assert dsolve(eq10) == sol10

    eq11 = (Eq(diff(x(t), t), f(t)*x(t) + g(t)*y(t)), Eq(diff(y(t), t), a*(f(t) + a*h(t))*x(t) + a*(g(t) - h(t))*y(t)))
    sol11 = [Eq(x(t), exp(-Integral(-a*g(t) - f(t), t))*(C1 - Integral(-exp(-a*Integral(h(t), t))*exp(Integral(-a*g(t) - f(t), t))*C2*g(t), t))),
             Eq(y(t), exp(-Integral(-a*g(t) - f(t), t))*a*(C1 - Integral(-exp(-a*Integral(h(t), t))*exp(Integral(-a*g(t) - f(t), t))*C2*g(t), t)) + exp(-a*Integral(g(t) - (a*g(t) - a*h(t))/a, t))*C1)]
    assert dsolve(eq11) == sol11


def test_linear_2eq_order1_nonhomog_linear():
    e = [Eq(diff(f(x), x), f(x) + g(x) + 5*x),
         Eq(diff(g(x), x), f(x) - g(x))]
    s = [Eq(f(x), -exp(sqrt(2)*x)*C2/(-sqrt(2) + 1) + (10*x - 5*sqrt(2))/Mul(8, 1 + sqrt(2), evaluate=False) + (10*x + 5*sqrt(2))/Mul(8, -sqrt(2) + 1, evaluate=False) - exp(-sqrt(2)*x)*C1/(1 + sqrt(2))), Eq(g(x), exp(sqrt(2)*x)*C2 - 5*x/2 + exp(-sqrt(2)*x)*C1)]
    assert dsolve(e) == s


def test_linear_2eq_order1_nonhomog():
    # Note: once implemented, add some tests esp. with resonance
    e = [Eq(diff(f(x), x), f(x) + exp(x)),
         Eq(diff(g(x), x), f(x) + g(x) + x*exp(x))]
    s = [Eq(f(x), exp(x)*C2 + exp(x)*x), Eq(g(x), exp(x)*C1 + exp(x)*C2*x + exp(x)*x**2)]
    assert dsolve(e) == s


def test_linear_2eq_order1_type2_degen():
    e = [Eq(diff(f(x), x), f(x) + 5),
         Eq(diff(g(x), x), f(x) + 7)]
    s1 = [Eq(f(x), exp(x)*C2 - 5), Eq(g(x), exp(x)*C2 + C1 + 2*x - 5)]
    s = dsolve(e)
    assert checksysodesol(e, s) == (True, [0, 0])
    assert s == s1


def test_dsolve_linear_2eq_order1_diag_triangular():
    e = [Eq(diff(f(x), x), f(x)),
         Eq(diff(g(x), x), g(x))]
    s1 = [Eq(f(x), C1*exp(x)), Eq(g(x), C2*exp(x))]
    s = dsolve(e)
    assert checksysodesol(e, s) == (True, [0, 0])
    assert s == s1

    e = [Eq(diff(f(x), x), 2*f(x)),
         Eq(diff(g(x), x), 3*f(x) + 7*g(x))]
    s1 = [Eq(f(x), -5*exp(2*x)*C1/3), Eq(g(x), exp(7*x)*C2 + exp(2*x)*C1)]
    s = dsolve(e)
    assert checksysodesol(e, s) == (True, [0, 0])
    assert s == s1


def test_sysode_linear_2eq_order1_type1_D_lt_0():
    e = [Eq(diff(f(x), x), -9*I*f(x) - 4*g(x)),
         Eq(diff(g(x), x), -4*I*g(x))]
    s = dsolve(e)
    assert checksysodesol(e, s) == (True, [0, 0])


def test_sysode_linear_2eq_order1_type1_D_lt_0_b_eq_0():
    e = [Eq(diff(f(x), x), -9*I*f(x)),
         Eq(diff(g(x), x), -4*I*g(x))]
    s1 = [Eq(f(x), C1*exp(-9*I*x)), Eq(g(x), C2*exp(-4*I*x))]
    s = dsolve(e)
    assert checksysodesol(e, s) == (True, [0, 0])
    assert s == s1


def test_linear_2eq_order2():
    x, y, z = symbols('x, y, z', cls=Function)
    k, l, m, n = symbols('k, l, m, n', integer=True)
    t, l = symbols('t, l')
    x0, y0 = symbols('x0, y0', cls=Function)

    eq1 = (Eq(diff(x(t), t, t), 5*x(t) + 43*y(t)), Eq(diff(y(t), t, t), x(t) + 9*y(t)))
    sol1 = [Eq(x(t), 43*C1*exp(t*RootOf(l**4 - 14*l**2 + 2, 0)) + 43*C2*exp(t*RootOf(l**4 - 14*l**2 + 2, 1)) +
               43*C3*exp(t*RootOf(l**4 - 14*l**2 + 2, 2)) + 43*C4*exp(t*RootOf(l**4 - 14*l**2 + 2, 3))),
            Eq(y(t), C1*(RootOf(l**4 - 14*l**2 + 2, 0)**2 - 5)*exp(t*RootOf(l**4 - 14*l**2 + 2, 0)) +
               C2*(RootOf(l**4 - 14*l**2 + 2, 1)**2 - 5)*exp(t*RootOf(l**4 - 14*l**2 + 2, 1)) +
               C3*(RootOf(l**4 - 14*l**2 + 2, 2)**2 - 5)*exp(t*RootOf(l**4 - 14*l**2 + 2, 2)) +
               C4*(RootOf(l**4 - 14*l**2 + 2, 3)**2 - 5)*exp(t*RootOf(l**4 - 14*l**2 + 2, 3)))]
    assert dsolve(eq1) == sol1

    eq1_1 = (Eq(diff(x(t), t, t), 5*x(t) - 4*y(t)), Eq(diff(y(t), t, t), x(t) + 9*y(t)))
    sol1_1 = [Eq(x(t), 2*exp(sqrt(7)*t)*C1*(-4*t + 4*sqrt(7)) - 8*exp(sqrt(7)*t)*C3*t + 2*exp(-sqrt(7)*t)*C2*(-4*t - 4*sqrt(7)) - 8*exp(-sqrt(7)*t)*C4*t), Eq(y(t), 4*exp(sqrt(7)*t)*C1*t + exp(sqrt(7)*t)*C3*(4*t + 4*sqrt(7)) + 4*exp(-sqrt(7)*t)*C2*t + exp(-sqrt(7)*t)*C4*(4*t - 4*sqrt(7)))]
    assert dsolve(eq1_1) == sol1_1

    eq1_2 = (Eq(diff(x(t), t, t), 5*x(t)), Eq(diff(y(t), t, t), x(t) + 5*y(t)))
    sol1_2 = [Eq(x(t), 2*sqrt(5)*exp(sqrt(5)*t)*C1 + 2*sqrt(5)*exp(-sqrt(5)*t)*C2), Eq(y(t), exp(sqrt(5)*t)*C1*t + exp(sqrt(5)*t)*C3 - exp(-sqrt(5)*t)*C2*t + exp(-sqrt(5)*t)*C4)]
    assert dsolve(eq1_2) == sol1_2

    eq1_3 = (Eq(diff(x(t), t, t), 5*x(t) - 4*y(t)), Eq(diff(y(t), t, t), 5*y(t)))
    sol1_3 = [Eq(x(t), -4*exp(sqrt(5)*t)*C1*t + exp(sqrt(5)*t)*C3 + 4*exp(-sqrt(5)*t)*C2*t + exp(-sqrt(5)*t)*C4), Eq(y(t), 2*sqrt(5)*exp(sqrt(5)*t)*C1 + 2*sqrt(5)*exp(-sqrt(5)*t)*C2)]
    assert dsolve(eq1_3) == sol1_3

    eq2 = (Eq(diff(x(t), t, t), 8*x(t)+3*y(t)+31), Eq(diff(y(t), t, t), 9*x(t)+7*y(t)+12))
    sol2 = [Eq(x(t), 3*C1*exp(t*RootOf(l**4 - 15*l**2 + 29, 0)) + 3*C2*exp(t*RootOf(l**4 - 15*l**2 + 29, 1)) +
               3*C3*exp(t*RootOf(l**4 - 15*l**2 + 29, 2)) + 3*C4*exp(t*RootOf(l**4 - 15*l**2 + 29, 3)) - 181/29),
            Eq(y(t), C1*(RootOf(l**4 - 15*l**2 + 29, 0)**2 - 8)*exp(t*RootOf(l**4 - 15*l**2 + 29, 0)) +
               C2*(RootOf(l**4 - 15*l**2 + 29, 1)**2 - 8)*exp(t*RootOf(l**4 - 15*l**2 + 29, 1)) +
               C3*(RootOf(l**4 - 15*l**2 + 29, 2)**2 - 8)*exp(t*RootOf(l**4 - 15*l**2 + 29, 2)) +
               C4*(RootOf(l**4 - 15*l**2 + 29, 3)**2 - 8)*exp(t*RootOf(l**4 - 15*l**2 + 29, 3)) + 183/29)]
    assert dsolve(eq2) == sol2

    eq3 = (Eq(diff(x(t), t, t) - 9*diff(y(t), t) + 7*x(t), 0), Eq(diff(y(t), t, t) + 9*diff(x(t), t) + 7*y(t), 0))
    sol3 = [Eq(x(t), C1*cos(t*(9/2 + sqrt(109)/2)) + C2*sin(t*(9/2 + sqrt(109)/2)) + C3*cos(t*(-sqrt(109)/2 + 9/2)) +
               C4*sin(t*(-sqrt(109)/2 + 9/2))), Eq(y(t), -C1*sin(t*(9/2 + sqrt(109)/2)) + C2*cos(t*(9/2 + sqrt(109)/2)) -
                                                   C3*sin(t*(-sqrt(109)/2 + 9/2)) + C4*cos(t*(-sqrt(109)/2 + 9/2)))]
    assert dsolve(eq3) == sol3

    eq4 = (Eq(diff(x(t), t, t), 9*t*diff(y(t), t)-9*y(t)), Eq(diff(y(t), t, t), 7*t*diff(x(t), t)-7*x(t)))
    sol4 = [Eq(x(t), C3*t + t*Integral((9*C1*exp(3*sqrt(7)*t**2/2) + 9*C2*exp(-3*sqrt(7)*t**2/2))/t**2, t)),
            Eq(y(t), C4*t + t*Integral((3*sqrt(7)*C1*exp(3*sqrt(7)*t**2/2) - 3*sqrt(7)*C2*exp(-3*sqrt(7)*t**2/2))/t**2, t))]
    assert dsolve(eq4) == sol4

    eq4_1 = (Eq(diff(x(t), t, t), -9*t*diff(y(t), t) + 9*y(t)), Eq(diff(y(t), t, t), 7*t*diff(x(t), t) - 7*x(t)))
    sol4_1 = [Eq(x(t), C3*t + t*Integral((-9*C1*cos(3*sqrt(7)*t**2/2) - 9*C2*sin(3*sqrt(7)*t**2/2))/t**2, t)), Eq(y(t), C4*t + t*Integral((-3*sqrt(7)*C1*sin(3*sqrt(7)*t**2/2) + 3*sqrt(7)*C2*cos(3*sqrt(7)*t**2/2))/t**2, t))]
    assert dsolve(eq4_1) == sol4_1

    eq5 = (Eq(diff(x(t), t, t), (log(t)+t**2)*diff(x(t), t)+(log(t)+t**2)*3*diff(y(t), t)), Eq(diff(y(t), t, t),
                                                                                               (log(t)+t**2)*2*diff(x(t), t)+(log(t)+t**2)*9*diff(y(t), t)))
    sol5 = [Eq(x(t), -sqrt(22)*(C1*Integral(exp((-sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) + C2 -
                                C3*Integral(exp((sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) - C4 -
                                (sqrt(22) + 5)*(C1*Integral(exp((-sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) + C2) +
                                (-sqrt(22) + 5)*(C3*Integral(exp((sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) + C4))/88),
            Eq(y(t), -sqrt(22)*(C1*Integral(exp((-sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) +
                                C2 - C3*Integral(exp((sqrt(22) + 5)*Integral(t**2 + log(t), t)), t) - C4)/44)]
    assert dsolve(eq5) == sol5

    eq6 = (Eq(diff(x(t), t, t), log(t)*t*diff(y(t), t) - log(t)*y(t)), Eq(diff(y(t), t, t), log(t)*t*diff(x(t), t) - log(t)*x(t)))
    sol6 = [Eq(x(t), C3*t + t*Integral((C1*exp(Integral(t*log(t), t)) +
                                        C2*exp(-Integral(t*log(t), t)))/t**2, t)), Eq(y(t), C4*t + t*Integral((C1*exp(Integral(t*log(t), t)) -
                                                                                                               C2*exp(-Integral(t*log(t), t)))/t**2, t))]
    assert dsolve(eq6) == sol6

    eq6_1 = (Eq(diff(x(t), t, t), -log(t)*t*diff(y(t), t) + log(t)*y(t)), Eq(diff(y(t), t, t), log(t)*t*diff(x(t), t) - log(t)*x(t)))
    sol6_1 = [Eq(x(t), C3*t + t*Integral((-C1*cos(Integral(t*log(t), t)) - C2*sin(Integral(t*log(t), t)))/t**2, t)), Eq(y(t), C4*t + t*Integral((-C1*sin(Integral(t*log(t), t)) + C2*cos(Integral(t*log(t), t)))/t**2, t))]
    assert dsolve(eq6_1) == sol6_1

    eq7 = (Eq(diff(x(t), t, t), log(t)*(t*diff(x(t), t) - x(t)) + exp(t)*(t*diff(y(t), t) - y(t))),
           Eq(diff(y(t), t, t), (t**2)*(t*diff(x(t), t) - x(t)) + (t)*(t*diff(y(t), t) - y(t))))
    sol7 = [Eq(x(t), C3*t + t*Integral((C1*x0(t) + C2*x0(t)*Integral(t*exp(t)*exp(Integral(t**2, t)) *
                                                                     exp(Integral(t*log(t), t))/x0(t)**2, t))/t**2, t)), Eq(y(t), C4*t + t*Integral((C1*y0(t) +
                                                                                                                                                     C2*(y0(t)*Integral(t*exp(t)*exp(Integral(t**2, t))*exp(Integral(t*log(t), t))/x0(t)**2, t) +
                                                                                                                                                         exp(Integral(t**2, t))*exp(Integral(t*log(t), t))/x0(t)))/t**2, t))]
    assert dsolve(eq7) == sol7

    eq8 = (Eq(diff(x(t), t, t), t*(4*x(t) + 9*y(t))), Eq(diff(y(t), t, t), t*(12*x(t) - 6*y(t))))
    sol8 = ('[Eq(x(t), -sqrt(133)*((-sqrt(133) - 1)*(C2*(133*t**8/24 - t**3/6 + sqrt(133)*t**3/2 + 1) + '
            'C1*t*(sqrt(133)*t**4/6 - t**3/12 + 1) + O(t**6)) - (-1 + sqrt(133))*(C2*(-sqrt(133)*t**3/6 - t**3/6 + 1) + '
            'C1*t*(-sqrt(133)*t**3/12 - t**3/12 + 1) + O(t**6)) - 4*C2*(133*t**8/24 - t**3/6 + sqrt(133)*t**3/2 + 1) + '
            '4*C2*(-sqrt(133)*t**3/6 - t**3/6 + 1) - 4*C1*t*(sqrt(133)*t**4/6 - t**3/12 + 1) + '
            '4*C1*t*(-sqrt(133)*t**3/12 - t**3/12 + 1) + O(t**6))/3192), Eq(y(t), -sqrt(133)*(-C2*(133*t**8/24 - t**3/6 + '
            'sqrt(133)*t**3/2 + 1) + C2*(-sqrt(133)*t**3/6 - t**3/6 + 1) - C1*t*(sqrt(133)*t**4/6 - t**3/12 + 1) + '
            'C1*t*(-sqrt(133)*t**3/12 - t**3/12 + 1) + O(t**6))/266)]')
    assert sstr(dsolve(eq8)) == sol8

    eq9 = (Eq(diff(x(t), t, t), t*(4*diff(x(t), t) + 9*diff(y(t), t))), Eq(diff(y(t), t, t), t*(12*diff(x(t), t) - 6*diff(y(t), t))))
    sol9 = [Eq(x(t), -sqrt(133)*(4*C1*Integral(exp((-sqrt(133) - 1)*Integral(t, t)), t) + 4*C2 -
                                 4*C3*Integral(exp((-1 + sqrt(133))*Integral(t, t)), t) - 4*C4 - (-1 + sqrt(133))*(C1*Integral(exp((-sqrt(133) -
                                                                                                                                    1)*Integral(t, t)), t) + C2) + (-sqrt(133) - 1)*(C3*Integral(exp((-1 + sqrt(133))*Integral(t, t)), t) +
                                                                                                                                                                                     C4))/3192), Eq(y(t), -sqrt(133)*(C1*Integral(exp((-sqrt(133) - 1)*Integral(t, t)), t) + C2 -
                                                                                                                                                                                                                      C3*Integral(exp((-1 + sqrt(133))*Integral(t, t)), t) - C4)/266)]
    assert dsolve(eq9) == sol9

    eq10 = (t**2*diff(x(t), t, t) + 3*t*diff(x(t), t) + 4*t*diff(y(t), t) + 12*x(t) + 9*y(t),
            t**2*diff(y(t), t, t) + 2*t*diff(x(t), t) - 5*t*diff(y(t), t) + 15*x(t) + 8*y(t))
    sol10 = [Eq(x(t), -C1*(-2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                               5*sqrt(70771857)/36)**Rational(1, 3)) + 13 + 2*sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                             4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 +
                                                                                                                                                                   346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))))*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                  4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                   5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) -
                C2*(-2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                    13 - 2*sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                      5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                        5*sqrt(70771857)/36)**Rational(1, 3))))*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                           2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                               4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) - C3*t**(1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              5*sqrt(70771857)/36)**Rational(1, 3)))/2)*(2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            5*sqrt(70771857)/36)**Rational(1, 3)) + 13 + 2*sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        5*sqrt(70771857)/36)**Rational(1, 3)))) - C4*t**(-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       5*sqrt(70771857)/36)**Rational(1, 3)))/2 + 1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       5*sqrt(70771857)/36)**Rational(1, 3))/2)*(-2*sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 5*sqrt(70771857)/36)**Rational(1, 3))) + 2*sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             5*sqrt(70771857)/36)**Rational(1, 3)) + 13)), Eq(y(t), C1*(-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 14 + (-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          5*sqrt(70771857)/36)**Rational(1, 3)))/2)**2 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 5*sqrt(70771857)/36)**Rational(1, 3))))*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 + sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) + C2*(-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 14 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    5*sqrt(70771857)/36)**Rational(1, 3))) + (-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      5*sqrt(70771857)/36)**Rational(1, 3)))/2)**2)*exp((-sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           5*sqrt(70771857)/36)**Rational(1, 3))/2 + 1 - sqrt(-284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) - 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 5*sqrt(70771857)/36)**Rational(1, 3)))/2)*log(t)) + C3*t**(1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2 + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   5*sqrt(70771857)/36)**Rational(1, 3)))/2)*(sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               5*sqrt(70771857)/36)**Rational(1, 3)) + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    5*sqrt(70771857)/36)**Rational(1, 3))) + 14 + (1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        5*sqrt(70771857)/36)**Rational(1, 3))/2 + sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 + 346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               5*sqrt(70771857)/36)**Rational(1, 3)))/2)**2) + C4*t**(-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)))/2 + 1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2)*(-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   8 + 346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))) + (-sqrt(-2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3) + 8 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 284/sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)))/2 + 1 + sqrt(-346/(3*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3))/2)**2 + sqrt(-346/(3*(Rational(4333, 4) +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  5*sqrt(70771857)/36)**Rational(1, 3)) + 4 + 2*(Rational(4333, 4) + 5*sqrt(70771857)/36)**Rational(1, 3)) + 14))]
    assert dsolve(eq10) == sol10

    eq11 = [x(t).diff(t, t) - 2*x(t) + 3*y(t) + 1,
            y(t).diff(t, t) + 2*x(t) - 3*y(t) - 2]
    sol11 = [Eq(x(t), (E**(sqrt(5)*t)*C1 - 3*C3*t - 3*C4 +
                       3*t**2/10 + Rational(8, 25) + E**(-sqrt(5)*t)*C2)),
             Eq(y(t), (-E**(sqrt(5)*t)*C1 - 2*C3*t - 2*C4 + t**2/5 -
                       Rational(8, 25) - E**(-sqrt(5)*t)*C2))]
    assert dsolve(eq11) == sol11
    eq12 = [x(t).diff(t, t) - 2*x(t) + 2*y(t) + 2,
            y(t).diff(t, t) - 2*x(t) + 2*y(t) + 1]
    sol12 = [Eq(x(t), -2*C1*t**3 - 2*C2*t**2 + C3*t + C4 - t**4/12 - t**2),
             Eq(y(t), (-2*C1*t**3 + 6*C1*t - 2*C2*t**2 + 2*C2 +
                       C3*t + C4 - t**4/12 - t**2/2))]
    assert dsolve(eq12) == sol12


def test_linear_3eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')
    eq1 = [Eq(diff(x(t), t), 21*x(t)), Eq(diff(y(t), t), 17*x(t) + 3*y(t)),
           Eq(diff(z(t), t), 5*x(t)+7*y(t)+9*z(t))]
    sol1 = [Eq(x(t), 216*exp(21*t)*C3/209),
            Eq(y(t), 204*exp(21*t)*C3/209 - 6*exp(3*t)*C1/7),
            Eq(z(t), exp(21*t)*C3 + exp(9*t)*C2 + exp(3*t)*C1)]
    assert dsolve(eq1) == sol1

    eq2 = [Eq(diff(x(t), t), 3*y(t) - 11*z(t)),
           Eq(diff(y(t), t), 7*z(t) - 3*x(t)),
           Eq(diff(z(t), t), 11*x(t)-7*y(t))]
    sol2 = [Eq(x(t), exp(sqrt(179)*I*t)*C3*(Rational(-21, 170) + 11*sqrt(179)*I/170) + 7*C1/3 + exp(-sqrt(179)*I*t)*C2*(Rational(-21, 170) - 11*sqrt(179)*I/170)),
            Eq(y(t), -sqrt(179)*exp(sqrt(179)*I*t)*I*C3*(7 - 33*sqrt(179)*I/179)/170 + 11*C1/3 + sqrt(179)*exp(-sqrt(179)*I*t)*I*C2*(7 + 33*sqrt(179)*I/179)/170), Eq(z(t), exp(sqrt(179)*I*t)*C3 + C1 + exp(-sqrt(179)*I*t)*C2)]
    assert dsolve(eq2) == sol2

    eq3 = [Eq(3*diff(x(t), t), 4*5*(y(t) - z(t))),
           Eq(4*diff(y(t), t), 3*5*(z(t) - x(t))),
           Eq(5*diff(z(t), t), 3*4*(x(t) - y(t)))]
    sol3 = [Eq(x(t), exp(5*sqrt(2)*I*t)*C3*(-1 + 4*sqrt(2)*I/3) + C1 + exp(-5*sqrt(2)*I*t)*C2*(-1 - 4*sqrt(2)*I/3)),
            Eq(y(t), -sqrt(2)*exp(5*sqrt(2)*I*t)*I*C3*(Rational(15, 4) - 5*sqrt(2)*I/2)/5 + C1 + sqrt(2)*exp(-5*sqrt(2)*I*t)*I*C2*(Rational(15, 4) + 5*sqrt(2)*I/2)/5),
            Eq(z(t), exp(5*sqrt(2)*I*t)*C3 + C1 + exp(-5*sqrt(2)*I*t)*C2)]
    assert dsolve(eq3) == sol3

    f = t**3 + log(t)
    g = t**2 + sin(t)
    eq4 = [Eq(diff(x(t), t), (4*f + g)*x(t) - f*y(t) - 2*f*z(t)),
           Eq(diff(y(t), t), 2*f*x(t) + (f + g)*y(t) - 2*f*z(t)),
           Eq(diff(z(t), t), 5*f*x(t) + f*y(t) + (-3*f + g)*z(t))]
    sol4 = [Eq(x(t), E**Integral(t**2 + sin(t), t)*(exp(sqrt(3)*I*Integral(-t**3 - log(t), t))*C3*(-2/(-4 - sqrt(3)*I) + (2 + 4/(-4 - sqrt(3)*I))/((-4 - sqrt(3)*I)*(Rational(-27, 19) - 17*sqrt(3)*I/19))) + exp(-2*Integral(-t**3 - log(t), t))*C1 + exp(-sqrt(3)*I*Integral(-t**3 - log(t), t))*C2*((2 + 4/(-4 + sqrt(3)*I))/((-4 + sqrt(3)*I)*(Rational(-27, 19) + 17*sqrt(3)*I/19)) - 2/(-4 + sqrt(3)*I)))),
            Eq(y(t), E**Integral(t**2 + sin(t), t)*(-exp(sqrt(3)*I*Integral(-t**3 - log(t), t))*C3*(2 + 4/(-4 - sqrt(3)*I))/(Rational(-27, 19) - 17*sqrt(3)*I/19) - exp(-sqrt(3)*I*Integral(-t**3 - log(t), t))*C2*(2 + 4/(-4 + sqrt(3)*I))/(Rational(-27, 19) + 17*sqrt(3)*I/19))),
            Eq(z(t), E**Integral(t**2 + sin(t), t)*(exp(sqrt(3)*I*Integral(-t**3 - log(t), t))*C3 + exp(-2*Integral(-t**3 - log(t), t))*C1 + exp(-sqrt(3)*I*Integral(-t**3 - log(t), t))*C2))]
    assert dsolve(eq4) == sol4

    eq5 = [Eq(diff(x(t), t), 4*x(t) - z(t)),
           Eq(diff(y(t), t), 2*x(t) + 2*y(t) - z(t)),
           Eq(diff(z(t), t), 3*x(t) + y(t))]
    sol5 = [Eq(x(t), -2*exp(2*t)*C1 + C2*(-2*exp(2*t)*t - 2*exp(2*t)) + C3*(-exp(2*t)*t**2 - 2*exp(2*t)*t - exp(2*t))),
            Eq(y(t), -2*exp(2*t)*C1 + C2*(-2*exp(2*t)*t - 2*exp(2*t)) + C3*(-exp(2*t)*t**2 - 2*exp(2*t)*t + exp(2*t))),
            Eq(z(t), -4*exp(2*t)*C1 + C2*(-4*exp(2*t)*t - 2*exp(2*t)) + C3*(-2*exp(2*t)*t**2 - 2*exp(2*t)*t))]
    assert dsolve(eq5) == sol5

    eq6 = [Eq(diff(x(t), t), 4*x(t) - y(t) - 2*z(t)),
           Eq(diff(y(t), t), 2*x(t) + y(t) - 2*z(t)),
           Eq(diff(z(t), t), 5*x(t) - 3*z(t))]
    sol6 = [Eq(x(t), exp(2*t)*C1 + exp(I*t)*C3*(Mul(-1, (-2 + 4/(4 - I))/((Rational(25, 17) - 15*I/17)*(4 - I)), evaluate=False) + 2/(4 - I)) + exp(-I*t)*C2*(2/(4 + I) - (-2 + 4/(4 + I))/((Rational(25, 17) + 15*I/17)*(4 + I)))),
            Eq(y(t), -exp(I*t)*C3*(-2 + 4/(4 - I))/(Rational(25, 17) - 15*I/17) - exp(-I*t)*C2*(-2 + 4/(4 + I))/(Rational(25, 17) + 15*I/17)),
            Eq(z(t), exp(2*t)*C1 + exp(I*t)*C3 + exp(-I*t)*C2)]
    assert dsolve(eq6) == sol6


def test_linear_3eq_order1_nonhomog():
    e = [Eq(diff(f(x), x), -9*f(x) - 4*g(x)),
         Eq(diff(g(x), x), -4*g(x)),
         Eq(diff(h(x), x), h(x) + exp(x))]
    s = [Eq(f(x), -4*exp(-4*x)*C2/5 + exp(-9*x)*C1),
         Eq(g(x), exp(-4*x)*C2), Eq(h(x), E**x*C3 + E**x*x)]
    assert dsolve(e) == s


def test_linear_3eq_order1_diagonal():
    # older code made assumptions about coefficients being nonzero
    e = [Eq(diff(f(x), x), f(x)),
         Eq(diff(g(x), x), g(x)),
         Eq(diff(h(x), x), h(x))]
    s1 = [Eq(f(x), C1*exp(x)), Eq(g(x), C2*exp(x)), Eq(h(x), C3*exp(x))]
    s = dsolve(e)
    assert s == s1


def test_nonlinear_2eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')
    eq1 = (Eq(diff(x(t), t), x(t)*y(t)**3), Eq(diff(y(t), t), y(t)**5))
    sol1 = [
        Eq(x(t), C1*exp((-1/(4*C2 + 4*t))**(-Rational(1, 4)))),
        Eq(y(t), -root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), C1*exp(-1/root(-1/(4*C2 + 4*t), 4))),
        Eq(y(t), root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), C1*exp(-I/root(-1/(4*C2 + 4*t), 4))),
        Eq(y(t), -I*root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), C1*exp(I/root(-1/(4*C2 + 4*t), 4))),
        Eq(y(t), I*root(-1/(4*C2 + 4*t), 4))]
    assert dsolve(eq1) == sol1

    eq2 = (Eq(diff(x(t), t), exp(3*x(t))*y(t)**3), Eq(diff(y(t), t), y(t)**5))
    sol2 = [
        Eq(x(t), -log(C1 - 3/root(-1/(4*C2 + 4*t), 4))/3),
        Eq(y(t), -root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), -log(C1 + 3/root(-1/(4*C2 + 4*t), 4))/3),
        Eq(y(t), root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), -log(C1 + 3*I/root(-1/(4*C2 + 4*t), 4))/3),
        Eq(y(t), -I*root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), -log(C1 - 3*I/root(-1/(4*C2 + 4*t), 4))/3),
        Eq(y(t), I*root(-1/(4*C2 + 4*t), 4))]
    assert dsolve(eq2) == sol2

    eq2_1 = (Eq(diff(x(t), t), exp(0*x(t))*x(t)*y(t)**2), Eq(diff(y(t), t), x(t)*y(t)))
    sol2_1 = [Eq(x(t), -exp(2*C1*(C2 + t))*C1/(exp(2*C1*(C2 + t)) - 1) + C1), Eq(y(t), -sqrt(2)*sqrt(-exp(2*C1*(C2 + t))*C1/(exp(2*C1*(C2 + t)) - 1))), Eq(x(t), -exp(2*C1*(C2 + t))*C1/(exp(2*C1*(C2 + t)) - 1) + C1), Eq(y(t), sqrt(2)*sqrt(-exp(2*C1*(C2 + t))*C1/(exp(2*C1*(C2 + t)) - 1)))]
    assert dsolve(eq2_1) == sol2_1

    eq3 = (Eq(diff(x(t), t), y(t)*x(t)), Eq(diff(y(t), t), x(t)**3))
    tt = Rational(2, 3)
    sol3 = [
        Eq(x(t), 6**tt/(6*(-sinh(sqrt(C1)*(C2 + t)/2)/sqrt(C1))**tt)),
        Eq(y(t), sqrt(C1 + C1/sinh(sqrt(C1)*(C2 + t)/2)**2)/3)]
    assert dsolve(eq3) == sol3

    eq4 = (Eq(diff(x(t), t), x(t)*y(t)*sin(t)**2), Eq(diff(y(t), t), y(t)**2*sin(t)**2))
    sol4 = {Eq(x(t), -2*exp(C1)/(C2*exp(C1) + t - sin(2*t)/2)), Eq(y(t), -2/(C1 + t - sin(2*t)/2))}
    assert dsolve(eq4) == sol4

    eq5 = (Eq(x(t), t*diff(x(t), t)+diff(x(t), t)*diff(y(t), t)), Eq(y(t), t*diff(y(t), t)+diff(y(t), t)**2))
    sol5 = {Eq(x(t), C1*C2 + C1*t), Eq(y(t), C2**2 + C2*t)}
    assert dsolve(eq5, [x(t), y(t)]) == sol5

    eq6 = (Eq(diff(x(t), t), x(t)**2*y(t)**3), Eq(diff(y(t), t), y(t)**5))
    sol6 = [
        Eq(x(t), 1/(C1 - 1/root(-1/(4*C2 + 4*t), 4))),
        Eq(y(t), -root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), 1/(C1 + (-1/(4*C2 + 4*t))**(-Rational(1, 4)))),
        Eq(y(t), root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), 1/(C1 + I/root(-1/(4*C2 + 4*t), 4))),
        Eq(y(t), -I*root(-1/(4*C2 + 4*t), 4)),
        Eq(x(t), 1/(C1 - I/root(-1/(4*C2 + 4*t), 4))),
        Eq(y(t), I*root(-1/(4*C2 + 4*t), 4))]
    assert dsolve(eq6) == sol6


def test_checksysodesol():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')
    eq = (Eq(diff(x(t), t), 9*y(t)), Eq(diff(y(t), t), 12*x(t)))
    sol = [Eq(x(t), 9*C1*exp(-6*sqrt(3)*t) + 9*C2*exp(6*sqrt(3)*t)),
           Eq(y(t), -6*sqrt(3)*C1*exp(-6*sqrt(3)*t) + 6*sqrt(3)*C2*exp(6*sqrt(3)*t))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (diff(x(t), t) - 9*y(t), diff(y(t), t) - 12*x(t))
    assert checksysodesol(eq, sol) == (True, [0, 0])
    assert checksysodesol(eq, sol, func=(x(t), y(t))) == (True, [0, 0])

    pytest.raises(ValueError, lambda: checksysodesol(eq, sol, func=(x(t), y(t, C0))))
    pytest.raises(ValueError, lambda: checksysodesol(eq, [Eq(x(t), y(t)), Eq(y(t), 0)]))
    pytest.raises(ValueError, lambda: checksysodesol(eq, (sol[0],)))

    assert checksysodesol(eq, [sol[0], Eq(sol[1].lhs, 0)])[0] is False

    # with change lhs to rhs:
    assert checksysodesol(eq, [_.reversed for _ in sol]) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), 2*x(t) + 4*y(t)), Eq(diff(y(t), t), 12*x(t) + 41*y(t)))
    sol = [Eq(x(t), 4*C1*exp(t*(-sqrt(1713)/2 + Rational(43, 2))) + 4*C2*exp(t*(sqrt(1713)/2 +
                                                                                Rational(43, 2)))), Eq(y(t), C1*(-sqrt(1713)/2 + Rational(39, 2))*exp(t*(-sqrt(1713)/2 +
                                                                                                                                                         Rational(43, 2))) + C2*(Rational(39, 2) + sqrt(1713)/2)*exp(t*(sqrt(1713)/2 + Rational(43, 2))))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), x(t) + y(t)), Eq(diff(y(t), t), -2*x(t) + 2*y(t)))
    sol = [Eq(x(t), (C1*sin(sqrt(7)*t/2) + C2*cos(sqrt(7)*t/2))*exp(3*t/2)),
           Eq(y(t), ((C1/2 - sqrt(7)*C2/2)*sin(sqrt(7)*t/2) + (sqrt(7)*C1/2 +
                                                               C2/2)*cos(sqrt(7)*t/2))*exp(3*t/2))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), x(t) + y(t) + 9), Eq(diff(y(t), t), 2*x(t) + 5*y(t) + 23))
    sol = [Eq(x(t), C1*exp(t*(-sqrt(6) + 3)) + C2*exp(t*(sqrt(6) + 3)) -
              Rational(22, 3)), Eq(y(t), C1*(-sqrt(6) + 2)*exp(t*(-sqrt(6) + 3)) + C2*(2 +
                                                                                       sqrt(6))*exp(t*(sqrt(6) + 3)) - Rational(5, 3))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), x(t) + y(t) + 81), Eq(diff(y(t), t), -2*x(t) + y(t) + 23))
    sol = [Eq(x(t), (C1*sin(sqrt(2)*t) + C2*cos(sqrt(2)*t))*exp(t) - Rational(58, 3)),
           Eq(y(t), (sqrt(2)*C1*cos(sqrt(2)*t) - sqrt(2)*C2*sin(sqrt(2)*t))*exp(t) - Rational(185, 3))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), 5*t*x(t) + 2*y(t)), Eq(diff(y(t), t), 2*x(t) + 5*t*y(t)))
    sol = [Eq(x(t), (C1*exp((Integral(2, t).doit())) + C2*exp(-(Integral(2, t)).doit())) *
              exp((Integral(5*t, t)).doit())), Eq(y(t), (C1*exp((Integral(2, t)).doit()) -
                                                         C2*exp(-(Integral(2, t)).doit()))*exp((Integral(5*t, t)).doit()))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), 5*t*x(t) + t**2*y(t)), Eq(diff(y(t), t), -t**2*x(t) + 5*t*y(t)))
    sol = [Eq(x(t), (C1*cos((Integral(t**2, t)).doit()) + C2*sin((Integral(t**2, t)).doit())) *
              exp((Integral(5*t, t)).doit())), Eq(y(t), (-C1*sin((Integral(t**2, t)).doit()) +
                                                         C2*cos((Integral(t**2, t)).doit()))*exp((Integral(5*t, t)).doit()))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), 5*t*x(t) + t**2*y(t)), Eq(diff(y(t), t), -t**2*x(t) + (5*t+9*t**2)*y(t)))
    sol = [Eq(x(t), (C1*exp((-sqrt(77)/2 + Rational(9, 2))*(Integral(t**2, t)).doit()) +
                     C2*exp((sqrt(77)/2 + Rational(9, 2))*(Integral(t**2, t)).doit()))*exp((Integral(5*t, t)).doit())),
           Eq(y(t), (C1*(-sqrt(77)/2 + Rational(9, 2))*exp((-sqrt(77)/2 + Rational(9, 2))*(Integral(t**2, t)).doit()) +
                     C2*(sqrt(77)/2 + Rational(9, 2))*exp((sqrt(77)/2 + Rational(9, 2))*(Integral(t**2, t)).doit()))*exp((Integral(5*t, t)).doit()))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t, t), 5*x(t) + 43*y(t)), Eq(diff(y(t), t, t), x(t) + 9*y(t)))
    root0 = -sqrt(-sqrt(47) + 7)
    root1 = sqrt(-sqrt(47) + 7)
    root2 = -sqrt(sqrt(47) + 7)
    root3 = sqrt(sqrt(47) + 7)
    sol = [Eq(x(t), 43*C1*exp(t*root0) + 43*C2*exp(t*root1) + 43*C3*exp(t*root2) + 43*C4*exp(t*root3)),
           Eq(y(t), C1*(root0**2 - 5)*exp(t*root0) + C2*(root1**2 - 5)*exp(t*root1) +
              C3*(root2**2 - 5)*exp(t*root2) + C4*(root3**2 - 5)*exp(t*root3))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t, t), 8*x(t)+3*y(t)+31), Eq(diff(y(t), t, t), 9*x(t)+7*y(t)+12))
    root0 = -sqrt(-sqrt(109)/2 + Rational(15, 2))
    root1 = sqrt(-sqrt(109)/2 + Rational(15, 2))
    root2 = -sqrt(sqrt(109)/2 + Rational(15, 2))
    root3 = sqrt(sqrt(109)/2 + Rational(15, 2))
    sol = [Eq(x(t), 3*C1*exp(t*root0) + 3*C2*exp(t*root1) + 3*C3*exp(t*root2) + 3*C4*exp(t*root3) - Rational(181, 29)),
           Eq(y(t), C1*(root0**2 - 8)*exp(t*root0) + C2*(root1**2 - 8)*exp(t*root1) +
              C3*(root2**2 - 8)*exp(t*root2) + C4*(root3**2 - 8)*exp(t*root3) + Rational(183, 29))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t, t) - 9*diff(y(t), t) + 7*x(t), 0), Eq(diff(y(t), t, t) + 9*diff(x(t), t) + 7*y(t), 0))
    sol = [Eq(x(t), C1*cos(t*(Rational(9, 2) + sqrt(109)/2)) + C2*sin(t*(Rational(9, 2) + sqrt(109)/2)) +
              C3*cos(t*(-sqrt(109)/2 + Rational(9, 2))) + C4*sin(t*(-sqrt(109)/2 + Rational(9, 2)))), Eq(y(t), -C1*sin(t*(Rational(9, 2) + sqrt(109)/2))
                                                                                                         + C2*cos(t*(Rational(9, 2) + sqrt(109)/2)) - C3*sin(t*(-sqrt(109)/2 + Rational(9, 2))) + C4*cos(t*(-sqrt(109)/2 + Rational(9, 2))))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t, t), 9*t*diff(y(t), t)-9*y(t)), Eq(diff(y(t), t, t), 7*t*diff(x(t), t)-7*x(t)))
    I1 = sqrt(6)*root(7, 4)*sqrt(pi)*erfi(sqrt(6)*root(7, 4)*t/2)/2 - exp(3*sqrt(7)*t**2/2)/t
    I2 = -sqrt(6)*root(7, 4)*sqrt(pi)*erf(sqrt(6)*root(7, 4)*t/2)/2 - exp(-3*sqrt(7)*t**2/2)/t
    sol = [Eq(x(t), C3*t + t*(9*C1*I1 + 9*C2*I2)), Eq(y(t), C4*t + t*(3*sqrt(7)*C1*I1 - 3*sqrt(7)*C2*I2))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), 21*x(t)), Eq(diff(y(t), t), 17*x(t)+3*y(t)), Eq(diff(z(t), t), 5*x(t)+7*y(t)+9*z(t)))
    sol = [Eq(x(t), C1*exp(21*t)), Eq(y(t), 17*C1*exp(21*t)/18 + C2*exp(3*t)),
           Eq(z(t), 209*C1*exp(21*t)/216 - 7*C2*exp(3*t)/6 + C3*exp(9*t))]
    assert checksysodesol(eq, sol) == (True, [0, 0, 0])

    eq = (Eq(diff(x(t), t), 3*y(t)-11*z(t)), Eq(diff(y(t), t), 7*z(t)-3*x(t)), Eq(diff(z(t), t), 11*x(t)-7*y(t)))
    sol = [Eq(x(t), 7*C0 + sqrt(179)*C1*cos(sqrt(179)*t) + (77*C1/3 + 130*C2/3)*sin(sqrt(179)*t)),
           Eq(y(t), 11*C0 + sqrt(179)*C2*cos(sqrt(179)*t) + (-58*C1/3 - 77*C2/3)*sin(sqrt(179)*t)),
           Eq(z(t), 3*C0 + sqrt(179)*(-7*C1/3 - 11*C2/3)*cos(sqrt(179)*t) + (11*C1 - 7*C2)*sin(sqrt(179)*t))]
    assert checksysodesol(eq, sol) == (True, [0, 0, 0])

    eq = (Eq(3*diff(x(t), t), 4*5*(y(t)-z(t))), Eq(4*diff(y(t), t), 3*5*(z(t)-x(t))), Eq(5*diff(z(t), t), 3*4*(x(t)-y(t))))
    sol = [Eq(x(t), C0 + 5*sqrt(2)*C1*cos(5*sqrt(2)*t) + (12*C1/5 + 164*C2/15)*sin(5*sqrt(2)*t)),
           Eq(y(t), C0 + 5*sqrt(2)*C2*cos(5*sqrt(2)*t) + (-51*C1/10 - 12*C2/5)*sin(5*sqrt(2)*t)),
           Eq(z(t), C0 + 5*sqrt(2)*(-9*C1/25 - 16*C2/25)*cos(5*sqrt(2)*t) + (12*C1/5 - 12*C2/5)*sin(5*sqrt(2)*t))]
    assert checksysodesol(eq, sol) == (True, [0, 0, 0])

    eq = (Eq(diff(x(t), t), 4*x(t) - z(t)), Eq(diff(y(t), t), 2*x(t)+2*y(t)-z(t)), Eq(diff(z(t), t), 3*x(t)+y(t)))
    sol = [Eq(x(t), C1*exp(2*t) + C2*t*exp(2*t) + C2*exp(2*t) + C3*t**2*exp(2*t)/2 + C3*t*exp(2*t) + C3*exp(2*t)),
           Eq(y(t), C1*exp(2*t) + C2*t*exp(2*t) + C2*exp(2*t) + C3*t**2*exp(2*t)/2 + C3*t*exp(2*t)),
           Eq(z(t), 2*C1*exp(2*t) + 2*C2*t*exp(2*t) + C2*exp(2*t) + C3*t**2*exp(2*t) + C3*t*exp(2*t) + C3*exp(2*t))]
    assert checksysodesol(eq, sol) == (True, [0, 0, 0])

    eq = (Eq(diff(x(t), t), 4*x(t) - y(t) - 2*z(t)), Eq(diff(y(t), t), 2*x(t) + y(t) - 2*z(t)), Eq(diff(z(t), t), 5*x(t)-3*z(t)))
    sol = [Eq(x(t), C1*exp(2*t) + C2*(-sin(t) + 3*cos(t)) + C3*(3*sin(t) + cos(t))),
           Eq(y(t), C2*(-sin(t) + 3*cos(t)) + C3*(3*sin(t) + cos(t))), Eq(z(t), C1*exp(2*t) + 5*C2*cos(t) + 5*C3*sin(t))]
    assert checksysodesol(eq, sol) == (True, [0, 0, 0])

    eq = (Eq(diff(x(t), t), x(t)*y(t)**3), Eq(diff(y(t), t), y(t)**5))
    sol = [Eq(x(t), C1*exp((-1/(4*C2 + 4*t))**(-Rational(1, 4)))), Eq(y(t), -root(-1/(4*C2 + 4*t), 4)),
           Eq(x(t), C1*exp(-1/root(-1/(4*C2 + 4*t), 4))), Eq(y(t), root(-1/(4*C2 + 4*t), 4)),
           Eq(x(t), C1*exp(-I/root(-1/(4*C2 + 4*t), 4))), Eq(y(t), -I*root(-1/(4*C2 + 4*t), 4)),
           Eq(x(t), C1*exp(I/root(-1/(4*C2 + 4*t), 4))), Eq(y(t), I*root(-1/(4*C2 + 4*t), 4))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(diff(x(t), t), exp(3*x(t))*y(t)**3), Eq(diff(y(t), t), y(t)**5))
    sol = [Eq(x(t), -log(C1 - 3/root(-1/(4*C2 + 4*t), 4))/3), Eq(y(t), -root(-1/(4*C2 + 4*t), 4)),
           Eq(x(t), -log(C1 + 3/root(-1/(4*C2 + 4*t), 4))/3), Eq(y(t), root(-1/(4*C2 + 4*t), 4)),
           Eq(x(t), -log(C1 + 3*I/root(-1/(4*C2 + 4*t), 4))/3), Eq(y(t), -I*root(-1/(4*C2 + 4*t), 4)),
           Eq(x(t), -log(C1 - 3*I/root(-1/(4*C2 + 4*t), 4))/3), Eq(y(t), I*root(-1/(4*C2 + 4*t), 4))]
    assert checksysodesol(eq, sol) == (True, [0, 0])

    eq = (Eq(x(t), t*diff(x(t), t)+diff(x(t), t)*diff(y(t), t)), Eq(y(t), t*diff(y(t), t)+diff(y(t), t)**2))
    sol = {Eq(x(t), C1*C2 + C1*t), Eq(y(t), C2**2 + C2*t)}
    assert checksysodesol(eq, sol) == (True, [0, 0])


def test_nonlinear_3eq_order1():
    x, y, z = symbols('x, y, z', cls=Function)
    t = Symbol('t')
    eq1 = (4*diff(x(t), t) + 2*y(t)*z(t), 3*diff(y(t), t) - z(t)*x(t), 5*diff(z(t), t) - x(t)*y(t))
    sol1 = '[Eq(x(t), Eq(Integral(4/(sqrt(-3*C1 + C2 - 4*_y**2)*sqrt(5*C1 - C2 - 4*_y**2)), (_y, x(t))), C3 + Integral(-sqrt(15)/15, t))), Eq(y(t), Eq(Integral(3/(sqrt(-C1 + 5*C2 - 6*_y**2)*sqrt(C1 - 4*C2 + 3*_y**2)), (_y, y(t))), C3 + Integral(sqrt(5)/10, t))), Eq(z(t), Eq(Integral(5/(sqrt(-3*C1 + C2 - 10*_y**2)*sqrt(4*C1 - C2 + 5*_y**2)), (_y, z(t))), C3 + Integral(sqrt(3)/6, t)))]'
    assert sstr(dsolve(eq1)) == sol1

    eq2 = (4*diff(x(t), t) + 2*y(t)*z(t)*sin(t), 3*diff(y(t), t) - z(t)*x(t)*sin(t), 5*diff(z(t), t) - x(t)*y(t)*sin(t))
    sol2 = '[Eq(x(t), Eq(Integral(3/(sqrt(-C1 + 5*C2 - 6*_y**2)*sqrt(C1 - 4*C2 + 3*_y**2)), (_y, x(t))), C3 + Integral(-sqrt(5)*sin(t)/10, t))), Eq(y(t), Eq(Integral(4/(sqrt(-3*C1 + C2 - 4*_y**2)*sqrt(5*C1 - C2 - 4*_y**2)), (_y, y(t))), C3 + Integral(sqrt(15)*sin(t)/15, t))), Eq(z(t), Eq(Integral(5/(sqrt(-3*C1 + C2 - 10*_y**2)*sqrt(4*C1 - C2 + 5*_y**2)), (_y, z(t))), C3 + Integral(-sqrt(3)*sin(t)/6, t)))]'
    assert sstr(dsolve(eq2)) == sol2


def test_checkodesol():
    # For the most part, checkodesol is well tested in the tests below.
    # These tests only handle cases not checked below.
    pytest.raises(ValueError, lambda: checkodesol(f(x, y).diff(x), Eq(f(x, y), x)))
    pytest.raises(ValueError, lambda: checkodesol(f(x).diff(x), Eq(f(x, y),
                                                                   x), f(x, y)))
    assert checkodesol(f(x).diff(x), Eq(f(x, y), x)) == \
        (False, -f(x).diff(x) + f(x, y).diff(x) - 1)
    assert checkodesol(f(x).diff(x), Eq(f(x), x)) is not True
    assert checkodesol(f(x).diff(x), Eq(f(x), x)) == (False, 1)
    sol1 = Eq(f(x)**5 + 11*f(x) - 2*f(x) + x, 0)
    assert checkodesol(diff(sol1.lhs, x), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, x)*exp(f(x)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, (x, 2)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, (x, 2))*exp(f(x)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, (x, 3)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, (x, 3))*exp(f(x)), sol1) == (True, 0)
    assert checkodesol(diff(sol1.lhs, (x, 3)), Eq(f(x), x*log(x))) == \
        (False, 60*x**4*((log(x) + 1)**2 + log(x))*(
            log(x) + 1)*log(x)**2 - 5*x**4*log(x)**4 - 9)
    assert checkodesol(diff(exp(f(x)) + x, x)*x, Eq(exp(f(x)) + x, 0)) == \
        (True, 0)
    assert checkodesol(diff(exp(f(x)) + x, x)*x, Eq(exp(f(x)) + x, 0),
                       solve_for_func=False) == (True, 0)
    assert checkodesol(f(x).diff((x, 2)), [Eq(f(x), C1 + C2*x),
                                           Eq(f(x), C2 + C1*x), Eq(f(x), C1*x + C2*x**2)]) == \
        [(True, 0), (True, 0), (False, 2*C2)]
    assert checkodesol(f(x).diff((x, 2)), {Eq(f(x), C1 + C2*x),
                                           Eq(f(x), C2 + C1*x), Eq(f(x), C1*x + C2*x**2)}) == \
        {(True, 0), (True, 0), (False, 2*C2)}
    assert checkodesol(f(x).diff(x) - 1/f(x)/2, Eq(f(x)**2, x)) == \
        [(True, 0), (True, 0)]
    assert checkodesol(f(x).diff(x) - f(x), Eq(C1*exp(x), f(x))) == (True, 0)
    # Based on test_1st_homogeneous_coeff_ode2_eq3sol.  Make sure that
    # checkodesol tries back substituting f(x) when it can.
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    assert not checkodesol(eq3, sol3)[1].has(f(x))


def test_dsolve_options():
    eq = x*f(x).diff(x) + f(x)
    a = dsolve(eq, hint='all')
    b = dsolve(eq, hint='all', simplify=False)
    c = dsolve(eq, hint='all_Integral')
    keys = ['1st_exact', '1st_exact_Integral', '1st_homogeneous_coeff_best',
            '1st_homogeneous_coeff_subs_dep_div_indep',
            '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
            '1st_homogeneous_coeff_subs_indep_div_dep',
            '1st_homogeneous_coeff_subs_indep_div_dep_Integral', '1st_linear',
            '1st_linear_Integral', 'almost_linear', 'almost_linear_Integral',
            'best', 'best_hint', 'default', 'lie_group',
            'nth_linear_euler_eq_homogeneous', 'order',
            'separable', 'separable_Integral']
    Integral_keys = ['1st_exact_Integral',
                     '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
                     '1st_homogeneous_coeff_subs_indep_div_dep_Integral', '1st_linear_Integral',
                     'almost_linear_Integral', 'best', 'best_hint', 'default',
                     'nth_linear_euler_eq_homogeneous',
                     'order', 'separable_Integral']
    assert sorted(a) == keys
    assert a['order'] == ode_order(eq, f(x))
    assert a['best'] == Eq(f(x), C1/x)
    assert dsolve(eq, hint='best') == Eq(f(x), C1/x)
    assert a['default'] == 'separable'
    assert a['best_hint'] == 'separable'
    assert not a['1st_exact'].has(Integral)
    assert not a['separable'].has(Integral)
    assert not a['1st_homogeneous_coeff_best'].has(Integral)
    assert not a['1st_homogeneous_coeff_subs_dep_div_indep'].has(Integral)
    assert not a['1st_homogeneous_coeff_subs_indep_div_dep'].has(Integral)
    assert not a['1st_linear'].has(Integral)
    assert a['1st_linear_Integral'].has(Integral)
    assert a['1st_exact_Integral'].has(Integral)
    assert a['1st_homogeneous_coeff_subs_dep_div_indep_Integral'].has(Integral)
    assert a['1st_homogeneous_coeff_subs_indep_div_dep_Integral'].has(Integral)
    assert a['separable_Integral'].has(Integral)
    assert sorted(b) == keys
    assert b['order'] == ode_order(eq, f(x))
    assert b['best'] == Eq(f(x), C1/x)
    assert dsolve(eq, hint='best', simplify=False) == Eq(f(x), C1/x)
    assert b['default'] == 'separable'
    assert b['best_hint'] == '1st_linear'
    assert a['separable'] != b['separable']
    assert a['1st_homogeneous_coeff_subs_dep_div_indep'] != \
        b['1st_homogeneous_coeff_subs_dep_div_indep']
    assert a['1st_homogeneous_coeff_subs_indep_div_dep'] != \
        b['1st_homogeneous_coeff_subs_indep_div_dep']
    assert not b['1st_exact'].has(Integral)
    assert not b['separable'].has(Integral)
    assert not b['1st_homogeneous_coeff_best'].has(Integral)
    assert not b['1st_homogeneous_coeff_subs_dep_div_indep'].has(Integral)
    assert not b['1st_homogeneous_coeff_subs_indep_div_dep'].has(Integral)
    assert not b['1st_linear'].has(Integral)
    assert b['1st_linear_Integral'].has(Integral)
    assert b['1st_exact_Integral'].has(Integral)
    assert b['1st_homogeneous_coeff_subs_dep_div_indep_Integral'].has(Integral)
    assert b['1st_homogeneous_coeff_subs_indep_div_dep_Integral'].has(Integral)
    assert b['separable_Integral'].has(Integral)
    assert sorted(c) == Integral_keys
    pytest.raises(ValueError, lambda: dsolve(eq, hint='notarealhint'))
    pytest.raises(ValueError, lambda: dsolve(eq, hint='Liouville'))
    assert dsolve(f(x).diff(x) - 1/f(x)**2, hint='all')['best'] == \
        dsolve(f(x).diff(x) - 1/f(x)**2, hint='best')
    assert dsolve(f(x) + f(x).diff(x) + sin(x).diff(x) + 1, f(x),
                  hint='1st_linear_Integral') == \
        Eq(f(x), (C1 + Integral((-sin(x).diff(x) - 1) *
                                exp(Integral(1, x)), x))*exp(-Integral(1, x)))


def test_classify_ode():
    assert classify_ode(f(x).diff((x, 2)), f(x)) == \
        ('nth_linear_constant_coeff_homogeneous', 'Liouville',
         '2nd_power_series_ordinary', 'Liouville_Integral')
    assert classify_ode(f(x), f(x)) == ()
    assert classify_ode(Eq(f(x).diff(x), 0), f(x)) == ('separable',
                                                       '1st_linear', '1st_homogeneous_coeff_best',
                                                       '1st_homogeneous_coeff_subs_indep_div_dep',
                                                       '1st_homogeneous_coeff_subs_dep_div_indep',
                                                       '1st_power_series', 'lie_group',
                                                       'nth_linear_constant_coeff_homogeneous',
                                                       'separable_Integral',
                                                       '1st_linear_Integral',
                                                       '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
                                                       '1st_homogeneous_coeff_subs_dep_div_indep_Integral')
    assert classify_ode(f(x).diff(x)**2, f(x)) == ('lie_group',)
    # issue sympy/sympy#4749: f(x) should be cleared from highest derivative before classifying
    a = classify_ode(Eq(f(x).diff(x) + f(x), x), f(x))
    b = classify_ode(f(x).diff(x)*f(x) + f(x)*f(x) - x*f(x), f(x))
    c = classify_ode(f(x).diff(x)/f(x) + f(x)/f(x) - x/f(x), f(x))
    assert a == ('1st_linear',
                 'Bernoulli',
                 'almost_linear',
                 '1st_power_series', 'lie_group',
                 'nth_linear_constant_coeff_undetermined_coefficients',
                 'nth_linear_constant_coeff_variation_of_parameters',
                 '1st_linear_Integral',
                 'Bernoulli_Integral',
                 'almost_linear_Integral',
                 'nth_linear_constant_coeff_variation_of_parameters_Integral')
    assert b == c != ()
    assert classify_ode(
        2*x*f(x)*f(x).diff(x) + (1 + x)*f(x)**2 - exp(x), f(x)
    ) == ('Bernoulli', 'almost_linear', 'lie_group',
          'Bernoulli_Integral', 'almost_linear_Integral')
    assert 'Riccati_special_minus2' in \
        classify_ode(2*f(x).diff(x) + f(x)**2 - f(x)/x + 3*x**(-2), f(x))
    pytest.raises(ValueError,
                  lambda: classify_ode(x + f(x, y).diff(x).diff(y), f(x, y)))
    # issue sympy/sympy#5176
    k = Symbol('k')
    assert classify_ode(f(x).diff(x)/(k*f(x) + k*x*f(x)) + 2*f(x)/(k*f(x) +
                                                                   k*x*f(x)) + x*f(x).diff(x)/(k*f(x) + k*x*f(x)) + z, f(x)) == \
        ('separable', '1st_exact', '1st_power_series', 'lie_group',
         'separable_Integral', '1st_exact_Integral')
    # preprocessing
    ans = ('separable', '1st_exact', '1st_linear', 'Bernoulli',
           '1st_homogeneous_coeff_best',
           '1st_homogeneous_coeff_subs_indep_div_dep',
           '1st_homogeneous_coeff_subs_dep_div_indep',
           '1st_power_series', 'lie_group',
           'nth_linear_constant_coeff_undetermined_coefficients',
           'nth_linear_constant_coeff_variation_of_parameters',
           'separable_Integral', '1st_exact_Integral',
           '1st_linear_Integral',
           'Bernoulli_Integral',
           '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
           '1st_homogeneous_coeff_subs_dep_div_indep_Integral',
           'nth_linear_constant_coeff_variation_of_parameters_Integral')
    #     w/o f(x) given
    assert classify_ode(diff(f(x) + x, x) + diff(f(x), x)) == ans
    #     w/ f(x) and prep=True
    assert classify_ode(diff(f(x) + x, x) + diff(f(x), x), f(x),
                        prep=True) == ans

    assert classify_ode(Eq(2*x**3*f(x).diff(x), 0), f(x)) == \
        ('separable', '1st_linear', '1st_power_series', 'lie_group',
         'separable_Integral', '1st_linear_Integral')
    assert classify_ode(Eq(2*f(x)**3*f(x).diff(x), 0), f(x)) == \
        ('separable', '1st_power_series', 'lie_group', 'separable_Integral')

    pytest.raises(ValueError, lambda: classify_ode(Derivative(f(x), x) +
                                                   Derivative(g(x), x)))

    # issue sympy/sympy4825
    pytest.raises(ValueError, lambda: dsolve(f(x, y).diff(x) - y*f(x, y), f(x)))
    assert classify_ode(f(x, y).diff(x) - y*f(x, y), f(x), dict=True) == \
        {'default': None, 'order': 0}
    # See also issue sympy/sympy#3793, test Z13.
    pytest.raises(ValueError, lambda: dsolve(f(x).diff(x), f(y)))
    assert classify_ode(f(x).diff(x), f(y), dict=True) == \
        {'default': None, 'order': 0}


def test_classify_ode_init():
    # Dummy
    eq = f(x).diff(x, x) - f(x)

    # Not f(0) or f'(0)
    init = {x: 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    ############################
    # f(0) type (AppliedUndef) #
    ############################

    # Wrong function
    init = {g(0): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Contains x
    init = {f(x): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Too many args
    init = {f(0, 0): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # point contains f
    # XXX: Should be NotImplementedError
    init = {f(0): f(1)}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Does not raise
    init = {f(0): 1}
    classify_ode(eq, f(x), init=init)

    #####################
    # f'(0) type (Subs) #
    #####################

    # Wrong function
    init = {g(x).diff(x).subs({x: 0}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Contains x
    init = {f(y).diff(y).subs({y: x}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Wrong variable
    init = {f(y).diff(y).subs({y: 0}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Too many args
    init = {f(x, y).diff(x).subs({x: 0}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Derivative wrt wrong vars
    init = {Derivative(f(x), x, y).subs({x: 0}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # point contains f
    # XXX: Should be NotImplementedError
    init = {f(x).diff(x).subs({x: 0}): f(0)}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Does not raise
    init = {f(x).diff(x).subs({x: 0}): 1}
    classify_ode(eq, f(x), init=init)

    ###########################
    # f'(y) type (Derivative) #
    ###########################

    # Wrong function
    init = {g(x).diff(x).subs({x: y}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Contains x
    init = {f(y).diff(y).subs({y: x}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Too many args
    init = {f(x, y).diff(x).subs({x: y}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Derivative wrt wrong vars
    init = {Derivative(f(x), x, z).subs({x: y}): 1}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # point contains f
    # XXX: Should be NotImplementedError
    init = {f(x).diff(x).subs({x: y}): f(0)}
    pytest.raises(ValueError, lambda: classify_ode(eq, f(x), init=init))

    # Does not raise
    init = {f(x).diff(x).subs({x: y}): 1}
    classify_ode(eq, f(x), init=init)


def test_classify_sysode():
    # Here x is assumed to be x(t) and y as y(t) for simplicity.
    # Similarly diff(x,t) and diff(y,y) is assumed to be x1 and y1 respectively.
    k, l, m, n = symbols('k, l, m, n', integer=True)
    k1, k2, k3, l1, l2, l3, m1, m2, m3 = symbols('k1, k2, k3, l1, l2, l3, m1, m2, m3', integer=True)
    P, Q, R, p, q, r = symbols('P, Q, R, p, q, r', cls=Function)
    P1, P2, P3, Q1, Q2, R1, R2 = symbols('P1, P2, P3, Q1, Q2, R1, R2', cls=Function)
    x, y, z = symbols('x, y, z', cls=Function)
    t = symbols('t')
    x1 = diff(x(t), t)
    y1 = diff(y(t), t)
    z1 = diff(z(t), t)
    x2 = diff(x(t), t, t)
    y2 = diff(y(t), t, t)

    eq1 = (Eq(diff(x(t), t), 5*t*x(t) + 2*y(t)), Eq(diff(y(t), t), 2*x(t) + 5*t*y(t)))
    sol1 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): -5*t, (1, x(t), 1): 0, (0, x(t), 1): 1,
                                                (1, y(t), 0): -5*t, (1, x(t), 0): -2, (0, y(t), 1): 0, (0, y(t), 0): -2, (1, y(t), 1): 1},
            'type_of_equation': 'type3', 'func': [x(t), y(t)], 'is_linear': True, 'eq': [-5*t*x(t) - 2*y(t) +
                                                                                         Derivative(x(t), t), -5*t*y(t) - 2*x(t) + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq1) == sol1

    eq2 = (Eq(x2, k*x(t) - l*y1), Eq(y2, l*x1 + k*y(t)))
    sol2 = {'order': {y(t): 2, x(t): 2}, 'type_of_equation': 'type3', 'is_linear': True, 'eq':
            [-k*x(t) + l*Derivative(y(t), t) + Derivative(x(t), t, t), -k*y(t) - l*Derivative(x(t), t) +
             Derivative(y(t), t, t)], 'no_of_equation': 2, 'func_coeff': {(0, y(t), 0): 0, (0, x(t), 2): 1,
                                                                          (1, y(t), 1): 0, (1, y(t), 2): 1, (1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): -k, (1, x(t), 1):
                                                                          -l, (0, x(t), 1): 0, (0, y(t), 1): l, (1, x(t), 0): 0, (1, y(t), 0): -k}, 'func': [x(t), y(t)]}
    assert classify_sysode(eq2) == sol2

    eq3 = (Eq(x2+4*x1+3*y1+9*x(t)+7*y(t), 11*exp(I*t)), Eq(y2+5*x1+8*y1+3*x(t)+12*y(t), 2*exp(I*t)))
    sol3 = {'no_of_equation': 2, 'func_coeff': {(1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): 9,
                                                (1, x(t), 1): 5, (0, x(t), 1): 4, (0, y(t), 1): 3, (1, x(t), 0): 3, (1, y(t), 0): 12, (0, y(t), 0): 7,
                                                (0, x(t), 2): 1, (1, y(t), 2): 1, (1, y(t), 1): 8}, 'type_of_equation': 'type4', 'func': [x(t), y(t)],
            'is_linear': True, 'eq': [9*x(t) + 7*y(t) - 11*exp(I*t) + 4*Derivative(x(t), t) + 3*Derivative(y(t), t) +
                                      Derivative(x(t), t, t), 3*x(t) + 12*y(t) - 2*exp(I*t) + 5*Derivative(x(t), t) + 8*Derivative(y(t), t) +
                                      Derivative(y(t), t, t)], 'order': {y(t): 2, x(t): 2}}
    assert classify_sysode(eq3) == sol3

    eq4 = (Eq((4*t**2 + 7*t + 1)**2*x2, 5*x(t) + 35*y(t)), Eq((4*t**2 + 7*t + 1)**2*y2, x(t) + 9*y(t)))
    sol4 = {'no_of_equation': 2, 'func_coeff': {(1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): -5,
                                                (1, x(t), 1): 0, (0, x(t), 1): 0, (0, y(t), 1): 0, (1, x(t), 0): -1, (1, y(t), 0): -9, (0, y(t), 0): -35,
                                                (0, x(t), 2): 16*t**4 + 56*t**3 + 57*t**2 + 14*t + 1, (1, y(t), 2): 16*t**4 + 56*t**3 + 57*t**2 + 14*t + 1,
                                                (1, y(t), 1): 0}, 'type_of_equation': 'type10', 'func': [x(t), y(t)], 'is_linear': True,
            'eq': [(4*t**2 + 7*t + 1)**2*Derivative(x(t), t, t) - 5*x(t) - 35*y(t), (4*t**2 + 7*t + 1)**2*Derivative(y(t), t, t)
                   - x(t) - 9*y(t)], 'order': {y(t): 2, x(t): 2}}
    assert classify_sysode(eq4) == sol4

    eq5 = (Eq(diff(x(t), t), x(t) + y(t) + 9), Eq(diff(y(t), t), 2*x(t) + 5*y(t) + 23))
    sol5 = {'eq': [-x(t) - y(t) + Derivative(x(t), t) - 9, -2*x(t) - 5*y(t) + Derivative(y(t), t) - 23], 'func': [x(t), y(t)], 'func_coeff': {(0, x(t), 0): -1, (0, x(t), 1): 1, (0, y(t), 0): -1, (0, y(t), 1): 0, (1, x(t), 0): -2, (1, x(t), 1): 0, (1, y(t), 0): -5, (1, y(t), 1): 1}, 'is_linear': True, 'no_of_equation': 2, 'order': {x(t): 1, y(t): 1}, 'type_of_equation': 'type1'}
    assert classify_sysode(eq5) == sol5

    eq6 = (Eq(x1, exp(k*x(t))*P(x(t), y(t))), Eq(y1, r(y(t))*P(x(t), y(t))))
    sol6 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): 0, (1, x(t), 1): 0, (0, x(t), 1): 1, (1, y(t), 0): 0,
                                                (1, x(t), 0): 0, (0, y(t), 1): 0, (0, y(t), 0): 0, (1, y(t), 1): 1}, 'type_of_equation': 'type2', 'func':
            [x(t), y(t)], 'is_linear': False, 'eq': [-P(x(t), y(t))*exp(k*x(t)) + Derivative(x(t), t), -P(x(t),
                                                                                                          y(t))*r(y(t)) + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq6) == sol6

    eq7 = (Eq(x1, x(t)**2+y(t)/x(t)), Eq(y1, x(t)/y(t)))
    sol7 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): 0, (1, x(t), 1): 0, (0, x(t), 1): 1, (1, y(t), 0): 0,
                                                (1, x(t), 0): -1/y(t), (0, y(t), 1): 0, (0, y(t), 0): -1/x(t), (1, y(t), 1): 1}, 'type_of_equation': 'type3',
            'func': [x(t), y(t)], 'is_linear': False, 'eq': [-x(t)**2 + Derivative(x(t), t) - y(t)/x(t), -x(t)/y(t) +
                                                             Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq7) == sol7

    eq8 = (Eq(x1, P1(x(t))*Q1(y(t))*R(x(t), y(t), t)), Eq(y1, P1(x(t))*Q1(y(t))*R(x(t), y(t), t)))
    sol8 = {'func': [x(t), y(t)], 'is_linear': False, 'type_of_equation': 'type4', 'eq':
            [-P1(x(t))*Q1(y(t))*R(x(t), y(t), t) + Derivative(x(t), t), -P1(x(t))*Q1(y(t))*R(x(t), y(t), t) +
             Derivative(y(t), t)], 'func_coeff': {(0, y(t), 1): 0, (1, y(t), 1): 1, (1, x(t), 1): 0, (0, y(t), 0): 0,
                                                  (1, x(t), 0): 0, (0, x(t), 0): 0, (1, y(t), 0): 0, (0, x(t), 1): 1}, 'order': {y(t): 1, x(t): 1}, 'no_of_equation': 2}
    assert classify_sysode(eq8) == sol8

    eq9 = (Eq(x1, 3*y(t)-11*z(t)), Eq(y1, 7*z(t)-3*x(t)), Eq(z1, 11*x(t)-7*y(t)))
    sol9 = {'no_of_equation': 3, 'func_coeff': {(1, y(t), 0): 0, (2, y(t), 1): 0, (2, z(t), 1): 1,
                                                (0, x(t), 0): 0, (2, x(t), 1): 0, (1, x(t), 1): 0, (2, y(t), 0): 7, (0, x(t), 1): 1, (1, z(t), 1): 0,
                                                (0, y(t), 1): 0, (1, x(t), 0): 3, (0, z(t), 0): 11, (0, y(t), 0): -3, (1, z(t), 0): -7, (0, z(t), 1): 0,
                                                (2, x(t), 0): -11, (2, z(t), 0): 0, (1, y(t), 1): 1}, 'type_of_equation': 'type1', 'func': [x(t), y(t), z(t)],
            'is_linear': True, 'eq': [-3*y(t) + 11*z(t) + Derivative(x(t), t), 3*x(t) - 7*z(t) + Derivative(y(t), t),
                                      -11*x(t) + 7*y(t) + Derivative(z(t), t)], 'order': {z(t): 1, y(t): 1, x(t): 1}}
    assert classify_sysode(eq9) == sol9

    eq10 = (x2 + log(t)*(t*x1 - x(t)) + exp(t)*(t*y1 - y(t)), y2 + (t**2)*(t*x1 - x(t)) + (t)*(t*y1 - y(t)))
    sol10 = {'no_of_equation': 2, 'func_coeff': {(1, x(t), 2): 0, (0, y(t), 2): 0, (0, x(t), 0): -log(t),
                                                 (1, x(t), 1): t**3, (0, x(t), 1): t*log(t), (0, y(t), 1): t*exp(t), (1, x(t), 0): -t**2, (1, y(t), 0): -t,
                                                 (0, y(t), 0): -exp(t), (0, x(t), 2): 1, (1, y(t), 2): 1, (1, y(t), 1): t**2}, 'type_of_equation': 'type11',
             'func': [x(t), y(t)], 'is_linear': True, 'eq': [(t*Derivative(x(t), t) - x(t))*log(t) + (t*Derivative(y(t), t) -
                                                                                                      y(t))*exp(t) + Derivative(x(t), t, t), t**2*(t*Derivative(x(t), t) - x(t)) + t*(t*Derivative(y(t), t) - y(t))
                                                             + Derivative(y(t), t, t)], 'order': {y(t): 2, x(t): 2}}
    assert classify_sysode(eq10) == sol10

    eq11 = (Eq(x1, x(t)*y(t)**3), Eq(y1, y(t)**5))
    sol11 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): -y(t)**3, (1, x(t), 1): 0, (0, x(t), 1): 1,
                                                 (1, y(t), 0): 0, (1, x(t), 0): 0, (0, y(t), 1): 0, (0, y(t), 0): 0, (1, y(t), 1): 1}, 'type_of_equation':
             'type1', 'func': [x(t), y(t)], 'is_linear': False, 'eq': [-x(t)*y(t)**3 + Derivative(x(t), t),
                                                                       -y(t)**5 + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq11) == sol11

    eq12 = (Eq(x1, y(t)), Eq(y1, x(t)))
    sol12 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): 0, (1, x(t), 1): 0, (0, x(t), 1): 1, (1, y(t), 0): 0,
                                                 (1, x(t), 0): -1, (0, y(t), 1): 0, (0, y(t), 0): -1, (1, y(t), 1): 1}, 'type_of_equation': 'type1', 'func':
             [x(t), y(t)], 'is_linear': True, 'eq': [-y(t) + Derivative(x(t), t), -x(t) + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq12) == sol12

    eq13 = (Eq(x1, x(t)*y(t)*sin(t)**2), Eq(y1, y(t)**2*sin(t)**2))
    sol13 = {'no_of_equation': 2, 'func_coeff': {(0, x(t), 0): -y(t)*sin(t)**2, (1, x(t), 1): 0, (0, x(t), 1): 1,
                                                 (1, y(t), 0): 0, (1, x(t), 0): 0, (0, y(t), 1): 0, (0, y(t), 0): -x(t)*sin(t)**2, (1, y(t), 1): 1},
             'type_of_equation': 'type4', 'func': [x(t), y(t)], 'is_linear': False, 'eq': [-x(t)*y(t)*sin(t)**2 +
                                                                                           Derivative(x(t), t), -y(t)**2*sin(t)**2 + Derivative(y(t), t)], 'order': {y(t): 1, x(t): 1}}
    assert classify_sysode(eq13) == sol13

    eq14 = (Eq(x1, 21*x(t)), Eq(y1, 17*x(t)+3*y(t)), Eq(z1, 5*x(t)+7*y(t)+9*z(t)))
    sol14 = {'no_of_equation': 3, 'func_coeff': {(1, y(t), 0): -3, (2, y(t), 1): 0, (2, z(t), 1): 1,
                                                 (0, x(t), 0): -21, (2, x(t), 1): 0, (1, x(t), 1): 0, (2, y(t), 0): -7, (0, x(t), 1): 1, (1, z(t), 1): 0,
                                                 (0, y(t), 1): 0, (1, x(t), 0): -17, (0, z(t), 0): 0, (0, y(t), 0): 0, (1, z(t), 0): 0, (0, z(t), 1): 0,
                                                 (2, x(t), 0): -5, (2, z(t), 0): -9, (1, y(t), 1): 1}, 'type_of_equation': 'type1', 'func': [x(t), y(t), z(t)],
             'is_linear': True, 'eq': [-21*x(t) + Derivative(x(t), t), -17*x(t) - 3*y(t) + Derivative(y(t), t), -5*x(t) -
                                       7*y(t) - 9*z(t) + Derivative(z(t), t)], 'order': {z(t): 1, y(t): 1, x(t): 1}}
    assert classify_sysode(eq14) == sol14

    eq15 = (Eq(x1, 4*x(t)+5*y(t)+2*z(t)), Eq(y1, x(t)+13*y(t)+9*z(t)), Eq(z1, 32*x(t)+41*y(t)+11*z(t)))
    sol15 = {'no_of_equation': 3, 'func_coeff': {(1, y(t), 0): -13, (2, y(t), 1): 0, (2, z(t), 1): 1,
                                                 (0, x(t), 0): -4, (2, x(t), 1): 0, (1, x(t), 1): 0, (2, y(t), 0): -41, (0, x(t), 1): 1, (1, z(t), 1): 0,
                                                 (0, y(t), 1): 0, (1, x(t), 0): -1, (0, z(t), 0): -2, (0, y(t), 0): -5, (1, z(t), 0): -9, (0, z(t), 1): 0,
                                                 (2, x(t), 0): -32, (2, z(t), 0): -11, (1, y(t), 1): 1}, 'type_of_equation': 'type1', 'func':
             [x(t), y(t), z(t)], 'is_linear': True, 'eq': [-4*x(t) - 5*y(t) - 2*z(t) + Derivative(x(t), t), -x(t) - 13*y(t) -
                                                           9*z(t) + Derivative(y(t), t), -32*x(t) - 41*y(t) - 11*z(t) + Derivative(z(t), t)], 'order': {z(t): 1, y(t): 1, x(t): 1}}
    assert classify_sysode(eq15) == sol15

    eq16 = (Eq(3*x1, 4*5*(y(t)-z(t))), Eq(4*y1, 3*5*(z(t)-x(t))), Eq(5*z1, 3*4*(x(t)-y(t))))
    sol16 = {'no_of_equation': 3, 'func_coeff': {(1, y(t), 0): 0, (2, y(t), 1): 0, (2, z(t), 1): 5,
                                                 (0, x(t), 0): 0, (2, x(t), 1): 0, (1, x(t), 1): 0, (2, y(t), 0): 12, (0, x(t), 1): 3, (1, z(t), 1): 0,
                                                 (0, y(t), 1): 0, (1, x(t), 0): 15, (0, z(t), 0): 20, (0, y(t), 0): -20, (1, z(t), 0): -15, (0, z(t), 1): 0,
                                                 (2, x(t), 0): -12, (2, z(t), 0): 0, (1, y(t), 1): 4}, 'type_of_equation': 'type1', 'func': [x(t), y(t), z(t)],
             'is_linear': True, 'eq': [-20*y(t) + 20*z(t) + 3*Derivative(x(t), t), 15*x(t) - 15*z(t) + 4*Derivative(y(t), t),
                                       -12*x(t) + 12*y(t) + 5*Derivative(z(t), t)], 'order': {z(t): 1, y(t): 1, x(t): 1}}
    assert classify_sysode(eq16) == sol16

    eq17 = (4*diff(x(t), t) + 2*y(t)*z(t), diff(y(t), t) - z(t)*x(t), 5*diff(z(t), t) - x(t)*y(t))
    assert classify_sysode(eq17)['type_of_equation'] is None

    eq18 = (x(t).diff(t) - t*x(t) + 1, y(t).diff(t) + y(t))
    assert classify_sysode(eq18)['type_of_equation'] is None

    eq18_1 = (x(t).diff(t) - t*x(t) + t, y(t).diff(t) + y(t))
    assert classify_sysode(eq18_1)['type_of_equation'] is None


def test_solve_init():
    # Basic tests that things work from dsolve.
    assert dsolve(f(x).diff(x) - f(x), f(x), init={f(0): 1}) == Eq(f(x), exp(x))
    assert dsolve(f(x).diff(x) - f(x), f(x), init={f(x).diff(x).subs({x: 0}): 1}) == Eq(f(x), exp(x))
    assert dsolve(f(x).diff(x, x) + f(x), f(x),
                  init={f(0): 1,
                        f(x).diff(x).subs({x: 0}): 1}) == Eq(f(x), sin(x) + cos(x))
    assert (dsolve([f(x).diff(x) - f(x) + g(x), g(x).diff(x) - g(x) - f(x)],
                   [f(x), g(x)], init={f(0): 1, g(0): 0}) ==
            [Eq(f(x), E**(x*(1 - I))/2 + E**(x*(1 + I))/2),
             Eq(g(x), E**(x*(1 - I))*I/2 - E**(x*(1 + I))*I/2)])

    assert (dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x), f(x),
                   init={f(0): 1}, hint='1st_exact', simplify=False) ==
            Eq(x*cos(f(x)) + f(x)**3/3, Rational(1, 3)))
    assert (dsolve(cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x), f(x),
                   init={f(0): 1}, hint='1st_exact', simplify=True) ==
            Eq(x*cos(f(x)) + f(x)**3/3, Rational(1, 3)))

    assert (dsolve(x + f(x)*Derivative(f(x), x), init={f(1): 0}) ==
            [Eq(f(x), -sqrt(-x**2 + 1)), Eq(f(x), sqrt(-x**2 + 1))])
    assert (dsolve(x + f(x)*Derivative(f(x), x), init={Subs(f(x), (x, 1)): 0}) ==
            [Eq(f(x), -sqrt(-x**2 + 1)), Eq(f(x), sqrt(-x**2 + 1))])

    assert solve_init([Eq(f(x), C1*exp(x))], [f(x)], [C1], {f(0): 1}) == {C1: 1}
    assert solve_init([Eq(f(x), C1*sin(x) + C2*cos(x))], [f(x)], [C1, C2],
                      {f(0): 1, f(pi/2): 1}) == {C1: 1, C2: 1}

    assert solve_init([Eq(f(x), C1*sin(x) + C2*cos(x))], [f(x)], [C1, C2],
                      {f(0): 1, f(x).diff(x).subs({x: 0}): 1}) == {C1: 1, C2: 1}

    # XXX: Ought to be ValueError
    pytest.raises(NotImplementedError, lambda: solve_init([Eq(f(x), C1*sin(x) + C2*cos(x))], [f(x)], [C1, C2], {f(0): 1, f(pi): 1}))

    EI, q, L = symbols('EI q L')

    # eq = Eq(EI*diff(f(x), x, 4), q)
    sols = [Eq(f(x), C1 + C2*x + C3*x**2 + C4*x**3 + q*x**4/(24*EI))]
    funcs = [f(x)]
    constants = [C1, C2, C3, C4]
    # Test both cases, Derivative (the default from f(x).diff(x).subs({x: L})),
    # and Subs
    init1 = {f(0): 0,
             f(x).diff(x).subs({x: 0}): 0,
             f(L).diff((L, 2)): 0,
             f(L).diff((L, 3)): 0}
    init2 = {f(0): 0,
             f(x).diff(x).subs({x: 0}): 0,
             Subs(f(x).diff((x, 2)), (x, L)): 0,
             Subs(f(x).diff((x, 3)), (x, L)): 0}

    solved_constants1 = solve_init(sols, funcs, constants, init1)
    solved_constants2 = solve_init(sols, funcs, constants, init2)
    assert solved_constants1 == solved_constants2 == {
        C1: 0,
        C2: 0,
        C3: L**2*q/(4*EI),
        C4: -L*q/(6*EI)}

    # Under-specified case:
    assert solve_init([Eq(f(x), E**(2*x)*C2 + E**(-2*x)*C1)],
                      [f(x)], {C1, C2}, {f(0): 1}) == {C1: 1 - C2}
    assert solve_init([Eq(f(x), C1*sin(x) + C2*cos(x))],
                      [f(x)], [C1, C2], {f(0): 1}) == {C2: 1}


def test_ode_order():
    f = Function('f')
    g = Function('g')
    x = Symbol('x')
    assert ode_order(3*x*exp(f(x)), f(x)) == 0
    assert ode_order(x*diff(f(x), x) + 3*x*f(x) - sin(x)/x, f(x)) == 1
    assert ode_order(x**2*f(x).diff(x, x) + x*diff(f(x), x) - f(x), f(x)) == 2
    assert ode_order(diff(x*exp(f(x)), x, x), f(x)) == 2
    assert ode_order(diff(x*diff(x*exp(f(x)), x, x), x), f(x)) == 3
    assert ode_order(diff(f(x), x, x), g(x)) == 0
    assert ode_order(diff(f(x), x, x)*diff(g(x), x), f(x)) == 2
    assert ode_order(diff(f(x), x, x)*diff(g(x), x), g(x)) == 1
    assert ode_order(diff(x*diff(x*exp(f(x)), x, x), x), g(x)) == 0
    # issue sympy/sympy#5835: ode_order has to also work for unevaluated derivatives
    # (ie, without using doit()).
    assert ode_order(Derivative(x*f(x), x), f(x)) == 1
    assert ode_order(x*sin(Derivative(x*f(x)**2, x, x)), f(x)) == 2
    assert ode_order(Derivative(x*Derivative(x*exp(f(x)), x, x), x), g(x)) == 0
    assert ode_order(Derivative(f(x), x, x), g(x)) == 0
    assert ode_order(Derivative(x*exp(f(x)), x, x), f(x)) == 2
    assert ode_order(Derivative(f(x), x, x)*Derivative(g(x), x), g(x)) == 1
    assert ode_order(Derivative(x*Derivative(f(x), x, x), x), f(x)) == 3
    assert ode_order(
        x*sin(Derivative(x*Derivative(f(x), x)**2, x, x)), f(x)) == 3


# In all tests below, checkodesol has the order option set to prevent
# superfluous calls to ode_order(), and the solve_for_func flag set to False
# because dsolve() already tries to solve for the function, unless the
# simplify=False option is set.
def test_old_ode_tests():
    # These are simple tests from the old ode module
    eq1 = Eq(f(x).diff(x), 0)
    eq2 = Eq(3*f(x).diff(x) - 5, 0)
    eq3 = Eq(3*f(x).diff(x), 5)
    eq4 = Eq(9*f(x).diff(x, x) + f(x), 0)
    eq5 = Eq(9*f(x).diff(x, x), f(x))
    # Type: a(x)f'(x)+b(x)*f(x)+c(x)=0
    eq6 = Eq(x**2*f(x).diff(x) + 3*x*f(x) - sin(x)/x, 0)
    eq7 = Eq(f(x).diff(x, x) - 3*diff(f(x), x) + 2*f(x), 0)
    # Type: 2nd order, constant coefficients (two real different roots)
    eq8 = Eq(f(x).diff(x, x) - 4*diff(f(x), x) + 4*f(x), 0)
    # Type: 2nd order, constant coefficients (two real equal roots)
    eq9 = Eq(f(x).diff(x, x) + 2*diff(f(x), x) + 3*f(x), 0)
    # Type: 2nd order, constant coefficients (two complex roots)
    eq10 = Eq(3*f(x).diff(x) - 1, 0)
    eq11 = Eq(x*f(x).diff(x) - 1, 0)
    sol1 = Eq(f(x), C1)
    sol2 = Eq(f(x), C1 + 5*x/3)
    sol3 = Eq(f(x), C1 + 5*x/3)
    sol4 = Eq(f(x), C1*sin(x/3) + C2*cos(x/3))
    sol5 = Eq(f(x), C1*exp(-x/3) + C2*exp(x/3))
    sol6 = Eq(f(x), (C1 - cos(x))/x**3)
    sol7 = Eq(f(x), (C1 + C2*exp(x))*exp(x))
    sol8 = Eq(f(x), (C1 + C2*x)*exp(2*x))
    sol9 = Eq(f(x), (C1*sin(x*sqrt(2)) + C2*cos(x*sqrt(2)))*exp(-x))
    sol10 = Eq(f(x), C1 + x/3)
    sol11 = Eq(f(x), C1 + log(x))
    assert dsolve(eq1) == sol1
    assert dsolve(eq1.lhs) == sol1
    assert dsolve(eq2) == sol2
    assert dsolve(eq3) == sol3
    assert dsolve(eq4) == sol4
    assert dsolve(eq5) == sol5
    assert dsolve(eq6) == sol6
    assert dsolve(eq7) == sol7
    assert dsolve(eq8) == sol8
    assert dsolve(eq9) == sol9
    assert dsolve(eq10) == sol10
    assert dsolve(eq11) == sol11
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=2, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=1, solve_for_func=False)[0]
    assert checkodesol(eq11, sol11, order=1, solve_for_func=False)[0]


@pytest.mark.slow
def test_1st_linear():
    # Type: first order linear form f'(x)+p(x)f(x)=q(x)
    eq = Eq(f(x).diff(x) + x*f(x), x**2)
    sol = Eq(f(x), (C1 + x*exp(x**2/2)
                    - sqrt(2)*sqrt(pi)*erfi(sqrt(2)*x/2)/2)*exp(-x**2/2))
    assert dsolve(eq, hint='1st_linear') == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_Bernoulli():
    # Type: Bernoulli, f'(x) + p(x)*f(x) == q(x)*f(x)**n
    eq = Eq(x*f(x).diff(x) + f(x) - f(x)**2, 0)
    sol = dsolve(eq, f(x), hint='Bernoulli')
    assert sol == Eq(f(x), 1/(x*(C1 + 1/x)))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_Riccati_special_minus2():
    # Type: Riccati special alpha = -2, a*dy/dx + b*y**2 + c*y/x +d/x**2
    eq = 2*f(x).diff(x) + f(x)**2 - f(x)/x + 3*x**(-2)
    sol = dsolve(eq, f(x), hint='Riccati_special_minus2')
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_1st_exact1():
    # Type: Exact differential equation, p(x,f) + q(x,f)*f' == 0,
    # where dp/df == dq/dx
    eq1 = sin(x)*cos(f(x)) + cos(x)*sin(f(x))*f(x).diff(x)
    eq2 = (2*x*f(x) + 1)/f(x) + (f(x) - x)/f(x)**2*f(x).diff(x)
    eq3 = 2*x + f(x)*cos(x) + (2*f(x) + sin(x) - sin(f(x)))*f(x).diff(x)
    eq4 = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
    eq5 = 2*x*f(x) + (x**2 + f(x)**2)*f(x).diff(x)
    sol1 = [Eq(f(x), -acos(C1/cos(x)) + 2*pi), Eq(f(x), acos(C1/cos(x)))]
    sol2 = Eq(f(x), exp(C1 - x**2 + LambertW(-x*exp(-C1 + x**2))))
    sol2b = Eq(log(f(x)) + x/f(x) + x**2, C1)
    sol3 = Eq(f(x)*sin(x) + cos(f(x)) + x**2 + f(x)**2, C1)
    sol4 = Eq(x*cos(f(x)) + f(x)**3/3, C1)
    sol5 = Eq(x**2*f(x) + f(x)**3/3, C1)
    assert dsolve(eq1, f(x), hint='1st_exact') == sol1
    assert dsolve(eq2, f(x), hint='1st_exact') == sol2
    assert dsolve(eq3, f(x), hint='1st_exact') == sol3
    assert dsolve(eq4, hint='1st_exact') == sol4
    assert dsolve(eq5, hint='1st_exact', simplify=False) == sol5
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    # issue sympy/sympy#5080 blocks the testing of this solution
    # assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2b, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]


def test_separable1():
    # test_separable1-5 are from Ordinary Differential Equations, Tenenbaum and
    # Pollard, pg. 55
    eq1 = f(x).diff(x) - f(x)
    eq2 = x*f(x).diff(x) - f(x)
    eq3 = f(x).diff(x) + sin(x)
    eq4 = f(x)**2 + 1 - (x**2 + 1)*f(x).diff(x)
    eq5 = f(x).diff(x)/tan(x) - f(x) - 2
    sol1 = Eq(f(x), C1*exp(x))
    sol2 = Eq(f(x), C1*x)
    sol3 = Eq(f(x), C1 + cos(x))
    sol4 = Eq(atan(f(x)), C1 + atan(x))
    sol5 = Eq(f(x), -2 + C1*sqrt(1 + tan(x)**2))
    assert dsolve(eq1, hint='separable') == sol1
    assert dsolve(eq2, hint='separable') == sol2
    assert dsolve(eq3, hint='separable') == sol3
    assert dsolve(eq4, hint='separable', simplify=False) == sol4
    assert dsolve(eq5, hint='separable') == simplify(sol5).expand()
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]


def test_separable2():
    a = Symbol('a')
    eq6 = f(x)*x**2*f(x).diff(x) - f(x)**3 - 2*x**2*f(x).diff(x)
    eq7 = f(x)**2 - 1 - (2*f(x) + x*f(x))*f(x).diff(x)
    eq8 = x*log(x)*f(x).diff(x) + sqrt(1 + f(x)**2)
    eq9 = exp(x + 1)*tan(f(x)) + cos(f(x))*f(x).diff(x)
    eq10 = (x*cos(f(x)) + x**2*sin(f(x))*f(x).diff(x) -
            a**2*sin(f(x))*f(x).diff(x))
    # solve() messes this one up a little bit, so lets test _Integral here
    # We have to test strings with _Integral because y is a dummy variable.
    sol6str = ('Eq(Integral((_y - 2)/_y**3, (_y, f(x))), '
               'C1 + Integral(x**(-2), x))')
    sol7 = Eq(-log(-1 + f(x)**2)/2, C1 - log(2 + x))
    sol8 = Eq(asinh(f(x)), C1 - log(log(x)))
    # integrate cannot handle the integral on the lhs (cos/tan)
    sol9str = ('Eq(Integral(cos(_y)/tan(_y), (_y, f(x))), '
               'C1 + Integral(-E*E**x, x))')
    sol10 = Eq(-log(-1 + sin(f(x))**2)/2, C1 - log(x**2 - a**2)/2)
    assert str(dsolve(eq6, hint='separable_Integral')) == sol6str
    assert dsolve(eq7, hint='separable', simplify=False) == sol7
    assert dsolve(eq8, hint='separable', simplify=False) == sol8
    assert str(dsolve(eq9, hint='separable_Integral')) == sol9str
    assert dsolve(eq10, hint='separable', simplify=False) == sol10
    assert checkodesol(eq7, sol7, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=1, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=1, solve_for_func=False)[0]


def test_separable3():
    eq11 = f(x).diff(x) - f(x)*tan(x)
    eq12 = (x - 1)*cos(f(x))*f(x).diff(x) - 2*x*sin(f(x))
    eq13 = f(x).diff(x) - f(x)*log(f(x))/tan(x)
    sol11 = Eq(f(x), C1*sqrt(1 + tan(x)**2))
    sol12 = Eq(log(-1 + cos(f(x))**2)/2, C1 + 2*x + 2*log(x - 1))
    sol13 = Eq(log(log(f(x))), C1 + log(cos(x)**2 - 1)/2)
    assert dsolve(eq11, hint='separable') == simplify(sol11)
    assert dsolve(eq12, hint='separable', simplify=False) == sol12
    assert dsolve(eq13, hint='separable', simplify=False) == sol13
    assert checkodesol(eq11, sol11, order=1, solve_for_func=False)[0]
    assert checkodesol(eq13, sol13, order=1, solve_for_func=False)[0]


def test_separable4():
    # This has a slow integral (1/((1 + y**2)*atan(y))), so we isolate it.
    eq14 = x*f(x).diff(x) + (1 + f(x)**2)*atan(f(x))
    sol14 = Eq(log(atan(f(x))), C1 - log(x))
    assert dsolve(eq14, hint='separable', simplify=False) == sol14
    assert checkodesol(eq14, sol14, order=1, solve_for_func=False)[0]


def test_separable5():
    eq15 = f(x).diff(x) + x*(f(x) + 1)
    eq16 = exp(f(x)**2)*(x**2 + 2*x + 1) + (x*f(x) + f(x))*f(x).diff(x)
    eq17 = f(x).diff(x) + f(x)
    eq18 = sin(x)*cos(2*f(x)) + cos(x)*sin(2*f(x))*f(x).diff(x)
    eq19 = (1 - x)*f(x).diff(x) - x*(f(x) + 1)
    eq20 = f(x)*diff(f(x), x) + x - 3*x*f(x)**2
    eq21 = f(x).diff(x) - exp(x + f(x))
    sol15 = Eq(f(x), -1 + C1*exp(-x**2/2))
    sol16 = Eq(-exp(-f(x)**2)/2, C1 - x - x**2/2)
    sol17 = Eq(f(x), C1*exp(-x))
    sol18 = Eq(-log(-1 + sin(2*f(x))**2)/4, C1 + log(-1 + sin(x)**2)/2)
    sol19 = Eq(f(x), (C1*exp(-x) - x + 1)/(x - 1))
    sol20 = Eq(log(-1 + 3*f(x)**2)/6, C1 + x**2/2)
    sol21 = Eq(-exp(-f(x)), C1 + exp(x))
    assert dsolve(eq15, hint='separable') == sol15
    assert dsolve(eq16, hint='separable', simplify=False) == sol16
    assert dsolve(eq17, hint='separable') == sol17
    assert dsolve(eq18, hint='separable', simplify=False) == sol18
    assert dsolve(eq19, hint='separable') == sol19
    assert dsolve(eq20, hint='separable', simplify=False) == sol20
    assert dsolve(eq21, hint='separable', simplify=False) == sol21
    assert checkodesol(eq15, sol15, order=1, solve_for_func=False)[0]
    assert checkodesol(eq16, sol16, order=1, solve_for_func=False)[0]
    assert checkodesol(eq17, sol17, order=1, solve_for_func=False)[0]
    assert checkodesol(eq18, sol18, order=1, solve_for_func=False)[0]
    assert checkodesol(eq19, sol19, order=1, solve_for_func=False)[0]
    assert checkodesol(eq20, sol20, order=1, solve_for_func=False)[0]
    assert checkodesol(eq21, sol21, order=1, solve_for_func=False)[0]


def test_separable_1_5_checkodesol():
    eq12 = (x - 1)*cos(f(x))*f(x).diff(x) - 2*x*sin(f(x))
    sol12 = Eq(-log(1 - cos(f(x))**2)/2, C1 - 2*x - 2*log(1 - x))
    assert checkodesol(eq12, sol12, order=1, solve_for_func=False)[0]


def test_separable6():
    eq1 = f(x).diff(x)*(1 - sin(f(x)))
    sol1 = Eq(f(x) + cos(f(x)), C1)
    assert dsolve(eq1) == sol1


def test_homogeneous_order():
    assert homogeneous_order(exp(y/x) + tan(y/x), x, y) == 0
    assert homogeneous_order(x**2 + sin(x)*cos(y), x, y) is None
    assert homogeneous_order(x - y - x*sin(y/x), x, y) == 1
    assert homogeneous_order((x*y + sqrt(x**4 + y**4) + x**2*(log(x) - log(y))) /
                             (pi*x**Rational(2, 3)*sqrt(y)**3), x, y) == Rational(-1, 6)
    assert homogeneous_order(y/x*cos(y/x) - x/y*sin(y/x) + cos(y/x), x, y) == 0
    assert homogeneous_order(f(x), x, f(x)) == 1
    assert homogeneous_order(f(x)**2, x, f(x)) == 2
    assert homogeneous_order(x*y*z, x, y) == 2
    assert homogeneous_order(x*y*z, x, y, z) == 3
    assert homogeneous_order(x**2*f(x)/sqrt(x**2 + f(x)**2), f(x)) is None
    assert homogeneous_order(f(x, y)**2, x, f(x, y), y) == 2
    assert homogeneous_order(f(x, y)**2, x, f(x), y) is None
    assert homogeneous_order(f(x, y)**2, x, f(x, y)) is None
    assert homogeneous_order(f(y, x)**2, x, y, f(x, y)) is None
    assert homogeneous_order(f(y), f(x), x) is None
    assert homogeneous_order(-f(x)/x + 1/sin(f(x) / x), f(x), x) == 0
    assert homogeneous_order(log(1/y) + log(x**2), x, y) is None
    assert homogeneous_order(log(1/y) + log(x), x, y) == 0
    assert homogeneous_order(log(x/y), x, y) == 0
    assert homogeneous_order(2*log(1/y) + 2*log(x), x, y) == 0
    a = Symbol('a')
    assert homogeneous_order(a*log(1/y) + a*log(x), x, y) == 0
    assert homogeneous_order(f(x).diff(x), x, y) is None
    assert homogeneous_order(-f(x).diff(x) + x, x, y) is None
    assert homogeneous_order(O(x), x, y) is None
    assert homogeneous_order(x + O(x**2), x, y) is None
    assert homogeneous_order(x**pi, x) == pi
    assert homogeneous_order(x**x, x) is None
    pytest.raises(ValueError, lambda: homogeneous_order(x*y))


@pytest.mark.slow
def test_1st_homogeneous_coeff_ode():
    # Type: First order homogeneous, y'=f(y/x)
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    eq8 = x + f(x) - (x - f(x))*f(x).diff(x)
    sol1 = Eq(log(x), C1 - log(f(x)*sin(f(x)/x)/x))
    sol2 = Eq(log(x), log(C1) + log(cos(f(x)/x) - 1)/2 - log(cos(f(x)/x) + 1)/2)
    sol3 = Eq(f(x), -exp(C1)*LambertW(-x*exp(-C1 + 1)))
    sol4 = Eq(log(f(x)), C1 - 2*exp(x/f(x)))
    sol5 = Eq(f(x), exp(2*C1 + LambertW(-2*x**4*exp(-4*C1))/2)/x)
    sol6 = Eq(log(x),
              C1 + exp(-f(x)/x)*sin(f(x)/x)/2 + exp(-f(x)/x)*cos(f(x)/x)/2)
    sol7 = Eq(log(f(x)), C1 - 2*sqrt(-x/f(x) + 1))
    sol8 = Eq(log(x), C1 - log(sqrt(1 + f(x)**2/x**2)) + atan(f(x)/x))
    assert dsolve(eq1, hint='1st_homogeneous_coeff_subs_dep_div_indep') == \
        sol1
    # indep_div_dep actually has a simpler solution for eq2,
    # but it runs too slow
    assert dsolve(eq2, hint='1st_homogeneous_coeff_subs_dep_div_indep',
                  simplify=False) == sol2
    assert dsolve(eq3, hint='1st_homogeneous_coeff_best') == sol3
    assert dsolve(eq4, hint='1st_homogeneous_coeff_best') == sol4
    assert dsolve(eq5, hint='1st_homogeneous_coeff_best') == sol5
    assert dsolve(eq6, hint='1st_homogeneous_coeff_subs_dep_div_indep') == \
        sol6
    assert dsolve(eq7, hint='1st_homogeneous_coeff_best') == sol7
    assert dsolve(eq8, hint='1st_homogeneous_coeff_best') == sol8
    # checks are below


@pytest.mark.slow
def test_1st_homogeneous_coeff_ode_check134568():
    # These are the checkodesols from test_homogeneous_coeff_ode().
    eq1 = f(x)/x*cos(f(x)/x) - (x/f(x)*sin(f(x)/x) + cos(f(x)/x))*f(x).diff(x)
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    eq4 = 2*f(x)*exp(x/f(x)) + f(x)*f(x).diff(x) - 2*x*exp(x/f(x))*f(x).diff(x)
    eq5 = 2*x**2*f(x) + f(x)**3 + (x*f(x)**2 - 2*x**3)*f(x).diff(x)
    eq6 = x*exp(f(x)/x) - f(x)*sin(f(x)/x) + x*sin(f(x)/x)*f(x).diff(x)
    eq8 = x + f(x) - (x - f(x))*f(x).diff(x)
    sol1 = Eq(f(x)*sin(f(x)/x), C1)
    sol4 = Eq(log(C1*f(x)) + 2*exp(x/f(x)), 0)
    sol3 = Eq(-f(x)/(1 + log(x/f(x))), C1)
    sol5 = Eq(log(C1*x*sqrt(1/x)*sqrt(f(x))) + x**2/(2*f(x)**2), 0)
    sol6 = Eq(-exp(-f(x)/x)*sin(f(x)/x)/2 + log(C1*x) -
              cos(f(x)/x)*exp(-f(x)/x)/2, 0)
    sol8 = Eq(-atan(f(x)/x) + log(C1*x*sqrt(1 + f(x)**2/x**2)), 0)
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=1, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=1, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=1, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check2():
    eq2 = x*f(x).diff(x) - f(x) - x*sin(f(x)/x)
    sol2 = Eq(x/tan(f(x)/(2*x)), C1)
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check3():
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    sol3a = Eq(f(x), x*exp(1 - LambertW(C1*x)))
    assert checkodesol(eq3, sol3a, solve_for_func=True)[0]


@pytest.mark.xfail
def test_1st_homogeneous_coeff_ode_check3_fail():
    eq3 = f(x) + (x*log(f(x)/x) - 2*x)*diff(f(x), x)
    sol3b = Eq(f(x), C1*LambertW(C2*x))
    assert checkodesol(eq3, sol3b, solve_for_func=True)[0]


def test_1st_homogeneous_coeff_ode_check7():
    eq7 = (x + sqrt(f(x)**2 - x*f(x)))*f(x).diff(x) - f(x)
    sol7 = Eq(log(C1*f(x)) + 2*sqrt(1 - x/f(x)), 0)
    assert checkodesol(eq7, sol7, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode2():
    eq1 = f(x).diff(x) - f(x)/x + 1/sin(f(x)/x)
    eq2 = x**2 + f(x)**2 - 2*x*f(x)*f(x).diff(x)
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol1 = [Eq(f(x), x*(-acos(C1 + log(x)) + 2*pi)), Eq(f(x), x*acos(C1 + log(x)))]
    sol2 = Eq(log(f(x)), log(C1) + log(x/f(x)) - log(x**2/f(x)**2 - 1))
    sol3 = Eq(f(x), log((1/(C1 - log(x)))**x))
    # specific hints are applied for speed reasons
    assert dsolve(eq1, hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol1
    assert dsolve(eq2, hint='1st_homogeneous_coeff_best', simplify=False) == sol2
    assert dsolve(eq3, hint='1st_homogeneous_coeff_subs_dep_div_indep') == sol3
    assert checkodesol(eq1, sol1, order=1, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=1, solve_for_func=False)[0]
    # test for eq3 is in test_1st_homogeneous_coeff_ode2_check3 below


def test_1st_homogeneous_coeff_ode2_check3():
    eq3 = x*exp(f(x)/x) + f(x) - x*f(x).diff(x)
    sol3 = Eq(f(x), log(log(C1/x)**(-x)))
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode_check9():
    _u2 = Dummy('u2')
    __a = Dummy('a')
    eq9 = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol9 = Eq(-Integral(-1/(-(1 - sqrt(1 - _u2**2))*_u2 + _u2), (_u2, __a,
                                                                 x/f(x))) + log(C1*f(x)), 0)
    assert checkodesol(eq9, sol9, order=1, solve_for_func=False)[0]


def test_1st_homogeneous_coeff_ode3():
    # The standard integration engine cannot handle one of the integrals
    # involved (see issue sympy/sympy#4551).  meijerg code comes up with an answer, but in
    # unconventional form.
    # checkodesol fails for this equation, so its test is in
    # test_1st_homogeneous_coeff_ode_check9 above. It has to compare string
    # expressions because u2 is a dummy variable.
    eq = f(x)**2 + (x*sqrt(f(x)**2 - x**2) - x*f(x))*f(x).diff(x)
    sol = Eq(log(f(x)), C1 - Piecewise(
            (-acosh(f(x)/x), abs(f(x)**2)/x**2 > 1),
            (I*asin(f(x)/x), True)))
    assert dsolve(eq, hint='1st_homogeneous_coeff_subs_indep_div_dep') == sol


def test_1st_homogeneous_coeff_corner_case():
    eq1 = f(x).diff(x) - f(x)/x
    c1 = classify_ode(eq1, f(x))
    eq2 = x*f(x).diff(x) - f(x)
    c2 = classify_ode(eq2, f(x))
    sdi = '1st_homogeneous_coeff_subs_dep_div_indep'
    sid = '1st_homogeneous_coeff_subs_indep_div_dep'
    assert sid not in c1 and sdi not in c1
    assert sid not in c2 and sdi not in c2


@pytest.mark.slow
def test_nth_linear_constant_coeff_homogeneous():
    # From Exercise 20, in Ordinary Differential Equations,
    #                      Tenenbaum and Pollard, pg. 220
    a = Symbol('a', positive=True)
    k = Symbol('k', extended_real=True)
    eq1 = f(x).diff((x, 2)) + 2*f(x).diff(x)
    eq2 = f(x).diff((x, 2)) - 3*f(x).diff(x) + 2*f(x)
    eq3 = f(x).diff((x, 2)) - f(x)
    eq4 = f(x).diff((x, 3)) + f(x).diff((x, 2)) - 6*f(x).diff(x)
    eq5 = 6*f(x).diff((x, 2)) - 11*f(x).diff(x) + 4*f(x)
    eq6 = Eq(f(x).diff((x, 2)) + 2*f(x).diff(x) - f(x), 0)
    eq7 = diff(f(x), (x, 3)) + diff(f(x), (x, 2)) - 10*diff(f(x), x) - 6*f(x)
    eq8 = f(x).diff((x, 4)) - f(x).diff((x, 3)) - 4*f(x).diff((x, 2)) + \
        4*f(x).diff(x)
    eq9 = f(x).diff((x, 4)) + 4*f(x).diff((x, 3)) + f(x).diff((x, 2)) - \
        4*f(x).diff(x) - 2*f(x)
    eq10 = f(x).diff((x, 4)) - a**2*f(x)
    eq11 = f(x).diff((x, 2)) - 2*k*f(x).diff(x) - 2*f(x)
    eq12 = f(x).diff((x, 2)) + 4*k*f(x).diff(x) - 12*k**2*f(x)
    eq13 = f(x).diff((x, 4))
    eq14 = f(x).diff((x, 2)) + 4*f(x).diff(x) + 4*f(x)
    eq15 = 3*f(x).diff((x, 3)) + 5*f(x).diff((x, 2)) + f(x).diff(x) - f(x)
    eq16 = f(x).diff((x, 3)) - 6*f(x).diff((x, 2)) + 12*f(x).diff(x) - 8*f(x)
    eq17 = f(x).diff((x, 2)) - 2*a*f(x).diff(x) + a**2*f(x)
    eq18 = f(x).diff((x, 4)) + 3*f(x).diff((x, 3))
    eq19 = f(x).diff((x, 4)) - 2*f(x).diff((x, 2))
    eq20 = f(x).diff((x, 4)) + 2*f(x).diff((x, 3)) - 11*f(x).diff((x, 2)) - \
        12*f(x).diff(x) + 36*f(x)
    eq21 = 36*f(x).diff((x, 4)) - 37*f(x).diff((x, 2)) + 4*f(x).diff(x) + 5*f(x)
    eq22 = f(x).diff((x, 4)) - 8*f(x).diff((x, 2)) + 16*f(x)
    eq23 = f(x).diff((x, 2)) - 2*f(x).diff(x) + 5*f(x)
    eq24 = f(x).diff((x, 2)) - f(x).diff(x) + f(x)
    eq25 = f(x).diff((x, 4)) + 5*f(x).diff((x, 2)) + 6*f(x)
    eq26 = f(x).diff((x, 2)) - 4*f(x).diff(x) + 20*f(x)
    eq27 = f(x).diff((x, 4)) + 4*f(x).diff((x, 2)) + 4*f(x)
    eq28 = f(x).diff((x, 3)) + 8*f(x)
    eq29 = f(x).diff((x, 4)) + 4*f(x).diff((x, 2))
    eq30 = f(x).diff((x, 5)) + 2*f(x).diff((x, 3)) + f(x).diff(x)
    sol1 = Eq(f(x), C1 + C2*exp(-2*x))
    sol2 = Eq(f(x), (C1 + C2*exp(x))*exp(x))
    sol3 = Eq(f(x), C1*exp(x) + C2*exp(-x))
    sol4 = Eq(f(x), C1 + C2*exp(-3*x) + C3*exp(2*x))
    sol5 = Eq(f(x), C1*exp(x/2) + C2*exp(4*x/3))
    sol6 = Eq(f(x), C1*exp(x*(-1 + sqrt(2))) + C2*exp(x*(-sqrt(2) - 1)))
    sol7 = Eq(f(x),
              C1*exp(3*x) + C2*exp(x*(-2 - sqrt(2))) + C3*exp(x*(-2 + sqrt(2))))
    sol8 = Eq(f(x), C1 + C2*exp(x) + C3*exp(-2*x) + C4*exp(2*x))
    sol9 = Eq(f(x),
              C1*exp(x) + C2*exp(-x) + C3*exp(x*(-2 + sqrt(2))) +
              C4*exp(x*(-2 - sqrt(2))))
    sol10 = Eq(f(x),
               C1*sin(x*sqrt(a)) + C2*cos(x*sqrt(a)) + C3*exp(x*sqrt(a)) +
               C4*exp(-x*sqrt(a)))
    sol11 = Eq(f(x),
               C1*exp(x*(k - sqrt(k**2 + 2))) + C2*exp(x*(k + sqrt(k**2 + 2))))
    sol12 = Eq(f(x), E**(-2*x*(k + 2*abs(k)))*C1 + E**(-2*x*(k - 2*abs(k)))*C2)
    sol13 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*x**3)
    sol14 = Eq(f(x), (C1 + C2*x)*exp(-2*x))
    sol15 = Eq(f(x), (C1 + C2*x)*exp(-x) + C3*exp(x/3))
    sol16 = Eq(f(x), (C1 + C2*x + C3*x**2)*exp(2*x))
    sol17 = Eq(f(x), (C1 + C2*x)*exp(a*x))
    sol18 = Eq(f(x), C1 + C2*x + C3*x**2 + C4*exp(-3*x))
    sol19 = Eq(f(x), C1 + C2*x + C3*exp(x*sqrt(2)) + C4*exp(-x*sqrt(2)))
    sol20 = Eq(f(x), (C1 + C2*x)*exp(-3*x) + (C3 + C4*x)*exp(2*x))
    sol21 = Eq(f(x), C1*exp(x/2) + C2*exp(-x) + C3*exp(-x/3) + C4*exp(5*x/6))
    sol22 = Eq(f(x), (C1 + C2*x)*exp(-2*x) + (C3 + C4*x)*exp(2*x))
    sol23 = Eq(f(x), (C1*sin(2*x) + C2*cos(2*x))*exp(x))
    sol24 = Eq(f(x), (C1*sin(x*sqrt(3)/2) + C2*cos(x*sqrt(3)/2))*exp(x/2))
    sol25 = Eq(f(x),
               C1*cos(x*sqrt(3)) + C2*sin(x*sqrt(3)) + C3*sin(x*sqrt(2)) +
               C4*cos(x*sqrt(2)))
    sol26 = Eq(f(x), (C1*sin(4*x) + C2*cos(4*x))*exp(2*x))
    sol27 = Eq(f(x), (C1 + C2*x)*sin(x*sqrt(2)) + (C3 + C4*x)*cos(x*sqrt(2)))
    sol28 = Eq(f(x),
               (C1*sin(x*sqrt(3)) + C2*cos(x*sqrt(3)))*exp(x) + C3*exp(-2*x))
    sol29 = Eq(f(x), C1 + C2*sin(2*x) + C3*cos(2*x) + C4*x)
    sol30 = Eq(f(x), C1 + (C2 + C3*x)*sin(x) + (C4 + C5*x)*cos(x))
    sol1s = constant_renumber(sol1, 'C', 1, 2)
    sol2s = constant_renumber(sol2, 'C', 1, 2)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 3)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    sol6s = constant_renumber(sol6, 'C', 1, 2)
    sol7s = constant_renumber(sol7, 'C', 1, 3)
    sol8s = constant_renumber(sol8, 'C', 1, 4)
    sol9s = constant_renumber(sol9, 'C', 1, 4)
    sol10s = constant_renumber(sol10, 'C', 1, 4)
    sol11s = constant_renumber(sol11, 'C', 1, 2)
    sol12s = constant_renumber(sol12, 'C', 1, 2)
    sol13s = constant_renumber(sol13, 'C', 1, 4)
    sol14s = constant_renumber(sol14, 'C', 1, 2)
    sol15s = constant_renumber(sol15, 'C', 1, 3)
    sol16s = constant_renumber(sol16, 'C', 1, 3)
    sol17s = constant_renumber(sol17, 'C', 1, 2)
    sol18s = constant_renumber(sol18, 'C', 1, 4)
    sol19s = constant_renumber(sol19, 'C', 1, 4)
    sol20s = constant_renumber(sol20, 'C', 1, 4)
    sol21s = constant_renumber(sol21, 'C', 1, 4)
    sol22s = constant_renumber(sol22, 'C', 1, 4)
    sol23s = constant_renumber(sol23, 'C', 1, 2)
    sol24s = constant_renumber(sol24, 'C', 1, 2)
    sol25s = constant_renumber(sol25, 'C', 1, 4)
    sol26s = constant_renumber(sol26, 'C', 1, 2)
    sol27s = constant_renumber(sol27, 'C', 1, 4)
    sol28s = constant_renumber(sol28, 'C', 1, 3)
    sol29s = constant_renumber(sol29, 'C', 1, 4)
    sol30s = constant_renumber(sol30, 'C', 1, 5)
    assert dsolve(eq1) in (sol1, sol1s)
    assert dsolve(eq2) in (sol2, sol2s)
    assert dsolve(eq3) in (sol3, sol3s)
    assert dsolve(eq4) in (sol4, sol4s)
    assert dsolve(eq5) in (sol5, sol5s)
    assert dsolve(eq6) in (sol6, sol6s)
    assert dsolve(eq7) in (sol7, sol7s)
    assert dsolve(eq8) in (sol8, sol8s)
    assert dsolve(eq9) in (sol9, sol9s)
    assert dsolve(eq10) in (sol10, sol10s)
    assert dsolve(eq11) in (sol11, sol11s)
    assert dsolve(eq12) in (sol12, sol12s)
    assert dsolve(eq13) in (sol13, sol13s)
    assert dsolve(eq14) in (sol14, sol14s)
    assert dsolve(eq15) in (sol15, sol15s)
    assert dsolve(eq16) in (sol16, sol16s)
    assert dsolve(eq17) in (sol17, sol17s)
    assert dsolve(eq18) in (sol18, sol18s)
    assert dsolve(eq19) in (sol19, sol19s)
    assert dsolve(eq20) in (sol20, sol20s)
    assert dsolve(eq21) in (sol21, sol21s)
    assert dsolve(eq22) in (sol22, sol22s)
    assert dsolve(eq23) in (sol23, sol23s)
    assert dsolve(eq24) in (sol24, sol24s)
    assert dsolve(eq25) in (sol25, sol25s)
    assert dsolve(eq26) in (sol26, sol26s)
    assert dsolve(eq27) in (sol27, sol27s)
    assert dsolve(eq28) in (sol28, sol28s)
    assert dsolve(eq29) in (sol29, sol29s)
    assert dsolve(eq30) in (sol30, sol30s)
    assert checkodesol(eq1, sol1, order=2, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=2, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=2, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=3, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=3, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=4, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=4, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=4, solve_for_func=False)[0]
    assert checkodesol(eq11, sol11, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, sol12, order=2, solve_for_func=False)[0]
    assert checkodesol(eq13, sol13, order=4, solve_for_func=False)[0]
    assert checkodesol(eq14, sol14, order=2, solve_for_func=False)[0]
    assert checkodesol(eq15, sol15, order=3, solve_for_func=False)[0]
    assert checkodesol(eq16, sol16, order=3, solve_for_func=False)[0]
    assert checkodesol(eq17, sol17, order=2, solve_for_func=False)[0]
    assert checkodesol(eq18, sol18, order=4, solve_for_func=False)[0]
    assert checkodesol(eq19, sol19, order=4, solve_for_func=False)[0]
    assert checkodesol(eq20, sol20, order=4, solve_for_func=False)[0]
    assert checkodesol(eq21, sol21, order=4, solve_for_func=False)[0]
    assert checkodesol(eq22, sol22, order=4, solve_for_func=False)[0]
    assert checkodesol(eq23, sol23, order=2, solve_for_func=False)[0]
    assert checkodesol(eq24, sol24, order=2, solve_for_func=False)[0]
    assert checkodesol(eq25, sol25, order=4, solve_for_func=False)[0]
    assert checkodesol(eq26, sol26, order=2, solve_for_func=False)[0]
    assert checkodesol(eq27, sol27, order=4, solve_for_func=False)[0]
    assert checkodesol(eq28, sol28, order=3, solve_for_func=False)[0]
    assert checkodesol(eq29, sol29, order=4, solve_for_func=False)[0]
    assert checkodesol(eq30, sol30, order=5, solve_for_func=False)[0]


def test_nth_linear_constant_coeff_homogeneous_RootOf():
    c = [C1, C2, C3, C4, C5]
    eq = f(x).diff((x, 5)) + 11*f(x).diff(x) - 2*f(x)
    sol = Eq(f(x), sum(exp(x*RootOf(x**5 + 11*x - 2, i))*c[i]
                       for i in range(5)))
    assert dsolve(eq) == sol
    assert checkodesol(eq, sol, order=5, solve_for_func=False)[0]


@pytest.mark.slow
def test_sympyissue_15520():
    c = [C1, C2, C3, C4, C5]
    eq = f(x).diff((x, 5)) + sqrt(3)*f(x).diff(x) - 2*f(x)
    sol = Eq(f(x), sum(exp(x*RootOf(x**5 + sqrt(3)*x - 2, i))*c[i]
                       for i in range(5)))
    assert dsolve(eq) == sol


def test_undetermined_coefficients_match():
    assert _undetermined_coefficients_match(g(x), x) == {'test': False}
    assert _undetermined_coefficients_match(sin(2*x + sqrt(5)), x) == \
        {'test': True, 'trialset':
            {cos(2*x + sqrt(5)), sin(2*x + sqrt(5))}}
    assert _undetermined_coefficients_match(sin(x)*cos(x), x) == \
        {'test': False}
    s = {cos(x), x*cos(x), x**2*cos(x), x**2*sin(x), x*sin(x), sin(x)}
    assert _undetermined_coefficients_match(sin(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': s}
    assert _undetermined_coefficients_match(
        sin(x)*x**2 + sin(x)*x + sin(x), x) == {'test': True, 'trialset': s}
    assert _undetermined_coefficients_match(
        exp(2*x)*sin(x)*(x**2 + x + 1), x
    ) == {
        'test': True, 'trialset': {exp(2*x)*sin(x), x**2*exp(2*x)*sin(x),
                                   cos(x)*exp(2*x), x**2*cos(x)*exp(2*x), x*cos(x)*exp(2*x),
                                   x*exp(2*x)*sin(x)}}
    assert _undetermined_coefficients_match(1/sin(x), x) == {'test': False}
    assert _undetermined_coefficients_match(log(x), x) == {'test': False}
    assert _undetermined_coefficients_match(2**(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': {2**x, x*2**x, x**2*2**x}}
    assert _undetermined_coefficients_match(x**y, x) == {'test': False}
    assert _undetermined_coefficients_match(exp(x)*exp(2*x + 1), x) == \
        {'test': True, 'trialset': {exp(1 + 3*x)}}
    assert _undetermined_coefficients_match(sin(x)*(x**2 + x + 1), x) == \
        {'test': True, 'trialset': {x*cos(x), x*sin(x), x**2*cos(x),
                                    x**2*sin(x), cos(x), sin(x)}}
    assert _undetermined_coefficients_match(sin(x)*(x + sin(x)), x) == \
        {'test': False}
    assert _undetermined_coefficients_match(sin(x)*(x + sin(2*x)), x) == \
        {'test': False}
    assert _undetermined_coefficients_match(sin(x)*tan(x), x) == \
        {'test': False}
    assert _undetermined_coefficients_match(
        x**2*sin(x)*exp(x) + x*sin(x) + x, x
    ) == {
        'test': True, 'trialset': {x**2*cos(x)*exp(x), x, cos(x), 1,
                                   exp(x)*sin(x), sin(x), x*exp(x)*sin(x), x*cos(x), x*cos(x)*exp(x),
                                   x*sin(x), cos(x)*exp(x), x**2*exp(x)*sin(x)}}
    assert _undetermined_coefficients_match(4*x*sin(x - 2), x) == {
        'trialset': {x*cos(x - 2), x*sin(x - 2), cos(x - 2), sin(x - 2)},
        'test': True,
    }
    assert _undetermined_coefficients_match(2**x*x, x) == \
        {'test': True, 'trialset': {2**x, x*2**x}}
    assert _undetermined_coefficients_match(2**x*exp(2*x), x) == \
        {'test': True, 'trialset': {2**x*exp(2*x)}}
    assert _undetermined_coefficients_match(exp(-x)/x, x) == \
        {'test': False}
    # Below are from Ordinary Differential Equations,
    #                Tenenbaum and Pollard, pg. 231
    assert _undetermined_coefficients_match(Integer(4), x) == \
        {'test': True, 'trialset': {1}}
    assert _undetermined_coefficients_match(12*exp(x), x) == \
        {'test': True, 'trialset': {exp(x)}}
    assert _undetermined_coefficients_match(exp(I*x), x) == \
        {'test': True, 'trialset': {exp(I*x)}}
    assert _undetermined_coefficients_match(sin(x), x) == \
        {'test': True, 'trialset': {cos(x), sin(x)}}
    assert _undetermined_coefficients_match(cos(x), x) == \
        {'test': True, 'trialset': {cos(x), sin(x)}}
    assert _undetermined_coefficients_match(8 + 6*exp(x) + 2*sin(x), x) == \
        {'test': True, 'trialset': {1, cos(x), sin(x), exp(x)}}
    assert _undetermined_coefficients_match(x**2, x) == \
        {'test': True, 'trialset': {1, x, x**2}}
    assert _undetermined_coefficients_match(9*x*exp(x) + exp(-x), x) == \
        {'test': True, 'trialset': {x*exp(x), exp(x), exp(-x)}}
    assert _undetermined_coefficients_match(2*exp(2*x)*sin(x), x) == \
        {'test': True, 'trialset': {exp(2*x)*sin(x), cos(x)*exp(2*x)}}
    assert _undetermined_coefficients_match(x - sin(x), x) == \
        {'test': True, 'trialset': {1, x, cos(x), sin(x)}}
    assert _undetermined_coefficients_match(x**2 + 2*x, x) == \
        {'test': True, 'trialset': {1, x, x**2}}
    assert _undetermined_coefficients_match(4*x*sin(x), x) == \
        {'test': True, 'trialset': {x*cos(x), x*sin(x), cos(x), sin(x)}}
    assert _undetermined_coefficients_match(x*sin(2*x), x) == \
        {'test': True, 'trialset':
            {x*cos(2*x), x*sin(2*x), cos(2*x), sin(2*x)}}
    assert _undetermined_coefficients_match(x**2*exp(-x), x) == \
        {'test': True, 'trialset': {x*exp(-x), x**2*exp(-x), exp(-x)}}
    assert _undetermined_coefficients_match(2*exp(-x) - x**2*exp(-x), x) == \
        {'test': True, 'trialset': {x*exp(-x), x**2*exp(-x), exp(-x)}}
    assert _undetermined_coefficients_match(exp(-2*x) + x**2, x) == \
        {'test': True, 'trialset': {1, x, x**2, exp(-2*x)}}
    assert _undetermined_coefficients_match(x*exp(-x), x) == \
        {'test': True, 'trialset': {x*exp(-x), exp(-x)}}
    assert _undetermined_coefficients_match(x + exp(2*x), x) == \
        {'test': True, 'trialset': {1, x, exp(2*x)}}
    assert _undetermined_coefficients_match(sin(x) + exp(-x), x) == \
        {'test': True, 'trialset': {cos(x), sin(x), exp(-x)}}
    assert _undetermined_coefficients_match(exp(x), x) == \
        {'test': True, 'trialset': {exp(x)}}
    # converted from sin(x)**2
    assert _undetermined_coefficients_match(Rational(1, 2) - cos(2*x)/2, x) == \
        {'test': True, 'trialset': {1, cos(2*x), sin(2*x)}}
    # converted from exp(2*x)*sin(x)**2
    assert _undetermined_coefficients_match(
        exp(2*x)*(Rational(1, 2) + cos(2*x)/2), x
    ) == {
        'test': True, 'trialset': {exp(2*x)*sin(2*x), cos(2*x)*exp(2*x),
                                   exp(2*x)}}
    assert _undetermined_coefficients_match(2*x + sin(x) + cos(x), x) == \
        {'test': True, 'trialset': {1, x, cos(x), sin(x)}}
    # converted from sin(2*x)*sin(x)
    assert _undetermined_coefficients_match(cos(x)/2 - cos(3*x)/2, x) == \
        {'test': True, 'trialset': {cos(x), cos(3*x), sin(x), sin(3*x)}}
    assert _undetermined_coefficients_match(cos(x**2), x) == {'test': False}
    assert _undetermined_coefficients_match(2**(x**2), x) == {'test': False}


@pytest.mark.slow
def test_nth_linear_constant_coeff_undetermined_coefficients():
    hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    g = exp(-x)
    f2 = f(x).diff((x, 2))
    c = 3*f(x).diff((x, 3)) + 5*f2 + f(x).diff(x) - f(x) - x
    eq1 = c - x*g
    eq2 = c - g
    # 3-27 below are from Ordinary Differential Equations,
    #                     Tenenbaum and Pollard, pg. 231
    eq3 = f2 + 3*f(x).diff(x) + 2*f(x) - 4
    eq4 = f2 + 3*f(x).diff(x) + 2*f(x) - 12*exp(x)
    eq5 = f2 + 3*f(x).diff(x) + 2*f(x) - exp(I*x)
    eq6 = f2 + 3*f(x).diff(x) + 2*f(x) - sin(x)
    eq7 = f2 + 3*f(x).diff(x) + 2*f(x) - cos(x)
    eq8 = f2 + 3*f(x).diff(x) + 2*f(x) - (8 + 6*exp(x) + 2*sin(x))
    eq9 = f2 + f(x).diff(x) + f(x) - x**2
    eq10 = f2 - 2*f(x).diff(x) - 8*f(x) - 9*x*exp(x) - 10*exp(-x)
    eq11 = f2 - 3*f(x).diff(x) - 2*exp(2*x)*sin(x)
    eq12 = f(x).diff((x, 4)) - 2*f2 + f(x) - x + sin(x)
    eq13 = f2 + f(x).diff(x) - x**2 - 2*x
    eq14 = f2 + f(x).diff(x) - x - sin(2*x)
    eq15 = f2 + f(x) - 4*x*sin(x)
    eq16 = f2 + 4*f(x) - x*sin(2*x)
    eq17 = f2 + 2*f(x).diff(x) + f(x) - x**2*exp(-x)
    eq18 = f(x).diff((x, 3)) + 3*f2 + 3*f(x).diff(x) + f(x) - 2*exp(-x) + \
        x**2*exp(-x)
    eq19 = f2 + 3*f(x).diff(x) + 2*f(x) - exp(-2*x) - x**2
    eq20 = f2 - 3*f(x).diff(x) + 2*f(x) - x*exp(-x)
    eq21 = f2 + f(x).diff(x) - 6*f(x) - x - exp(2*x)
    eq22 = f2 + f(x) - sin(x) - exp(-x)
    eq23 = f(x).diff((x, 3)) - 3*f2 + 3*f(x).diff(x) - f(x) - exp(x)
    # sin(x)**2
    eq24 = f2 + f(x) - Rational(1, 2) - cos(2*x)/2
    # exp(2*x)*sin(x)**2
    eq25 = f(x).diff((x, 3)) - f(x).diff(x) - exp(2*x)*(Rational(1, 2) - cos(2*x)/2)
    eq26 = (f(x).diff((x, 5)) + 2*f(x).diff((x, 3)) + f(x).diff(x) - 2*x -
            sin(x) - cos(x))
    # sin(2*x)*sin(x), skip 3127 for now, match bug
    eq27 = f2 + f(x) - cos(x)/2 + cos(3*x)/2
    eq28 = f(x).diff(x) - 1
    sol1 = Eq(f(x),
              -1 - x + (C1 + C2*x - 3*x**2/32 - x**3/24)*exp(-x) + C3*exp(x/3))
    sol2 = Eq(f(x), -1 - x + (C1 + C2*x - x**2/8)*exp(-x) + C3*exp(x/3))
    sol3 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol4 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), C1*exp(-2*x) + C2*exp(-x) + exp(I*x)/10 - 3*I*exp(I*x)/10)
    sol6 = Eq(f(x), -3*cos(x)/10 + sin(x)/10 + C1*exp(-x) + C2*exp(-2*x))
    sol7 = Eq(f(x), cos(x)/10 + 3*sin(x)/10 + C1*exp(-x) + C2*exp(-2*x))
    sol8 = Eq(f(x),
              4 - 3*cos(x)/5 + sin(x)/5 + exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol9 = Eq(f(x),
              -2*x + x**2 + (C1*sin(x*sqrt(3)/2) + C2*cos(x*sqrt(3)/2))*exp(-x/2))
    sol10 = Eq(f(x), -x*exp(x) - 2*exp(-x) + C1*exp(-2*x) + C2*exp(4*x))
    sol11 = Eq(f(x), C1 + C2*exp(3*x) - (3*sin(x) + cos(x))*exp(2*x)/5)
    sol12 = Eq(f(x), x - sin(x)/4 + (C1 + C2*x)*exp(-x) + (C3 + C4*x)*exp(x))
    sol13 = Eq(f(x), C1 + x**3/3 + C2*exp(-x))
    sol14 = Eq(f(x), C1 - x - sin(2*x)/5 - cos(2*x)/10 + x**2/2 + C2*exp(-x))
    sol15 = Eq(f(x), (C1 + x)*sin(x) + (C2 - x**2)*cos(x))
    sol16 = Eq(f(x), (C1 + x/16)*sin(2*x) + (C2 - x**2/8)*cos(2*x))
    sol17 = Eq(f(x), (C1 + C2*x + x**4/12)*exp(-x))
    sol18 = Eq(f(x), (C1 + C2*x + C3*x**2 - x**5/60 + x**3/3)*exp(-x))
    sol19 = Eq(f(x), Rational(7, 4) - 3*x/2 + x**2/2 + C1*exp(-x) + (C2 - x)*exp(-2*x))
    sol20 = Eq(f(x), C1*exp(x) + C2*exp(2*x) + (6*x + 5)*exp(-x)/36)
    sol21 = Eq(f(x), -Rational(1, 36) - x/6 + C1*exp(-3*x) + (C2 + x/5)*exp(2*x))
    sol22 = Eq(f(x), C1*sin(x) + (C2 - x/2)*cos(x) + exp(-x)/2)
    sol23 = Eq(f(x), (C1 + C2*x + C3*x**2 + x**3/6)*exp(x))
    sol24 = Eq(f(x), Rational(1, 2) - cos(2*x)/6 + C1*sin(x) + C2*cos(x))
    sol25 = Eq(f(x), C1 + C2*exp(-x) + C3*exp(x) +
               (-21*sin(2*x) + 27*cos(2*x) + 130)*exp(2*x)/1560)
    sol26 = Eq(f(x),
               C1 + (C2 + C3*x - x**2/8)*sin(x) + (C4 + C5*x + x**2/8)*cos(x) + x**2)
    sol27 = Eq(f(x), cos(3*x)/16 + C1*cos(x) + (C2 + x/4)*sin(x))
    sol28 = Eq(f(x), C1 + x)
    sol1s = constant_renumber(sol1, 'C', 1, 3)
    sol2s = constant_renumber(sol2, 'C', 1, 3)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 2)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    sol6s = constant_renumber(sol6, 'C', 1, 2)
    sol7s = constant_renumber(sol7, 'C', 1, 2)
    sol8s = constant_renumber(sol8, 'C', 1, 2)
    sol9s = constant_renumber(sol9, 'C', 1, 2)
    sol10s = constant_renumber(sol10, 'C', 1, 2)
    sol11s = constant_renumber(sol11, 'C', 1, 2)
    sol12s = constant_renumber(sol12, 'C', 1, 2)
    sol13s = constant_renumber(sol13, 'C', 1, 4)
    sol14s = constant_renumber(sol14, 'C', 1, 2)
    sol15s = constant_renumber(sol15, 'C', 1, 2)
    sol16s = constant_renumber(sol16, 'C', 1, 2)
    sol17s = constant_renumber(sol17, 'C', 1, 2)
    sol18s = constant_renumber(sol18, 'C', 1, 3)
    sol19s = constant_renumber(sol19, 'C', 1, 2)
    sol20s = constant_renumber(sol20, 'C', 1, 2)
    sol21s = constant_renumber(sol21, 'C', 1, 2)
    sol22s = constant_renumber(sol22, 'C', 1, 2)
    sol23s = constant_renumber(sol23, 'C', 1, 3)
    sol24s = constant_renumber(sol24, 'C', 1, 2)
    sol25s = constant_renumber(sol25, 'C', 1, 3)
    sol26s = constant_renumber(sol26, 'C', 1, 5)
    sol27s = constant_renumber(sol27, 'C', 1, 2)
    assert dsolve(eq1, hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, hint=hint) in (sol2, sol2s)
    assert dsolve(eq3, hint=hint) in (sol3, sol3s)
    assert dsolve(eq4, hint=hint) in (sol4, sol4s)
    assert dsolve(eq5, hint=hint) in (sol5, sol5s)
    assert dsolve(eq6, hint=hint) in (sol6, sol6s)
    assert dsolve(eq7, hint=hint) in (sol7, sol7s)
    assert dsolve(eq8, hint=hint) in (sol8, sol8s)
    assert dsolve(eq9, hint=hint) in (sol9, sol9s)
    assert dsolve(eq10, hint=hint) in (sol10, sol10s)
    assert dsolve(eq11, hint=hint) in (sol11, sol11s)
    assert dsolve(eq12, hint=hint) in (sol12, sol12s)
    assert dsolve(eq13, hint=hint) in (sol13, sol13s)
    assert dsolve(eq14, hint=hint) in (sol14, sol14s)
    assert dsolve(eq15, hint=hint) in (sol15, sol15s)
    assert dsolve(eq16, hint=hint) in (sol16, sol16s)
    assert dsolve(eq17, hint=hint) in (sol17, sol17s)
    assert dsolve(eq18, hint=hint) in (sol18, sol18s)
    assert dsolve(eq19, hint=hint) in (sol19, sol19s)
    assert dsolve(eq20, hint=hint) in (sol20, sol20s)
    assert dsolve(eq21, hint=hint) in (sol21, sol21s)
    assert dsolve(eq22, hint=hint) in (sol22, sol22s)
    assert dsolve(eq23, hint=hint) in (sol23, sol23s)
    assert dsolve(eq24, hint=hint) in (sol24, sol24s)
    assert dsolve(eq25, hint=hint) in (sol25, sol25s)
    assert dsolve(eq26, hint=hint) in (sol26, sol26s)
    assert dsolve(eq27, hint=hint) in (sol27, sol27s)
    assert dsolve(eq28, hint=hint) == sol28
    assert checkodesol(eq1, sol1, order=3, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=3, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=2, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=2, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=2, solve_for_func=False)[0]
    assert checkodesol(eq11, sol11, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, sol12, order=4, solve_for_func=False)[0]
    assert checkodesol(eq13, sol13, order=2, solve_for_func=False)[0]
    assert checkodesol(eq14, sol14, order=2, solve_for_func=False)[0]
    assert checkodesol(eq15, sol15, order=2, solve_for_func=False)[0]
    assert checkodesol(eq16, sol16, order=2, solve_for_func=False)[0]
    assert checkodesol(eq17, sol17, order=2, solve_for_func=False)[0]
    assert checkodesol(eq18, sol18, order=3, solve_for_func=False)[0]
    assert checkodesol(eq19, sol19, order=2, solve_for_func=False)[0]
    assert checkodesol(eq20, sol20, order=2, solve_for_func=False)[0]
    assert checkodesol(eq21, sol21, order=2, solve_for_func=False)[0]
    assert checkodesol(eq22, sol22, order=2, solve_for_func=False)[0]
    assert checkodesol(eq23, sol23, order=3, solve_for_func=False)[0]
    assert checkodesol(eq24, sol24, order=2, solve_for_func=False)[0]
    assert checkodesol(eq25, sol25, order=3, solve_for_func=False)[0]
    assert checkodesol(eq26, sol26, order=5, solve_for_func=False)[0]
    assert checkodesol(eq27, sol27, order=2, solve_for_func=False)[0]
    assert checkodesol(eq28, sol28, order=1, solve_for_func=False)[0]


def test_sympyissue_5787():
    # This test case is to show the classification of imaginary constants under
    # nth_linear_constant_coeff_undetermined_coefficients
    eq = Eq(diff(f(x), x), I*f(x) + Rational(1, 2) - I)
    out_hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    assert out_hint in classify_ode(eq)


@pytest.mark.xfail
def test_nth_linear_constant_coeff_undetermined_coefficients_imaginary_exp():
    # Equivalent to eq26 in
    # test_nth_linear_constant_coeff_undetermined_coefficients above.
    # This fails because the algorithm for undetermined coefficients
    # doesn't know to multiply exp(I*x) by sufficient x because it is linearly
    # dependent on sin(x) and cos(x).
    hint = 'nth_linear_constant_coeff_undetermined_coefficients'
    eq26a = f(x).diff(x, 5) + 2*f(x).diff(x, 3) + f(x).diff(x) - 2*x - exp(I*x)
    # sol26 = Eq(f(x),
    #            C1 + (C2 + C3*x - x**2/8)*sin(x) + (C4 + C5*x + x**2/8)*cos(x) + x**2)
    dsolve(eq26a, hint=hint)  # == sol26
    # assert checkodesol(eq26a, sol26, order=5, solve_for_func=False)[0]


@pytest.mark.slow
def test_nth_linear_constant_coeff_variation_of_parameters():
    hint = 'nth_linear_constant_coeff_variation_of_parameters'
    g = exp(-x)
    f2 = f(x).diff((x, 2))
    c = 3*f(x).diff((x, 3)) + 5*f2 + f(x).diff(x) - f(x) - x
    eq1 = c - x*g
    eq2 = c - g
    eq3 = f(x).diff(x) - 1
    eq4 = f2 + 3*f(x).diff(x) + 2*f(x) - 4
    eq5 = f2 + 3*f(x).diff(x) + 2*f(x) - 12*exp(x)
    eq6 = f2 - 2*f(x).diff(x) - 8*f(x) - 9*x*exp(x) - 10*exp(-x)
    eq7 = f2 + 2*f(x).diff(x) + f(x) - x**2*exp(-x)
    eq8 = f2 - 3*f(x).diff(x) + 2*f(x) - x*exp(-x)
    eq9 = f(x).diff((x, 3)) - 3*f2 + 3*f(x).diff(x) - f(x) - exp(x)
    eq10 = f2 + 2*f(x).diff(x) + f(x) - exp(-x)/x
    eq11 = f2 + f(x) - 1/sin(x)*1/cos(x)
    eq12 = f(x).diff((x, 4)) - 1/x
    sol1 = Eq(f(x),
              -1 - x + (C1 + C2*x - 3*x**2/32 - x**3/24)*exp(-x) + C3*exp(x/3))
    sol2 = Eq(f(x), -1 - x + (C1 + C2*x - x**2/8)*exp(-x) + C3*exp(x/3))
    sol3 = Eq(f(x), C1 + x)
    sol4 = Eq(f(x), 2 + C1*exp(-x) + C2*exp(-2*x))
    sol5 = Eq(f(x), 2*exp(x) + C1*exp(-x) + C2*exp(-2*x))
    sol6 = Eq(f(x), -x*exp(x) - 2*exp(-x) + C1*exp(-2*x) + C2*exp(4*x))
    sol7 = Eq(f(x), (C1 + C2*x + x**4/12)*exp(-x))
    sol8 = Eq(f(x), C1*exp(x) + C2*exp(2*x) + (6*x + 5)*exp(-x)/36)
    sol9 = Eq(f(x), (C1 + C2*x + C3*x**2 + x**3/6)*exp(x))
    sol10 = Eq(f(x), (C1 + x*(C2 + log(x)))*exp(-x))
    sol11 = Eq(f(x), cos(x)*(C2 - Integral(1/cos(x), x)) + sin(x)*(C1 +
                                                                   Integral(1/sin(x), x)))
    sol12 = Eq(f(x), C1 + C2*x + x**3*(C3 + log(x)/6) + C4*x**2)
    sol1s = constant_renumber(sol1, 'C', 1, 3)
    sol2s = constant_renumber(sol2, 'C', 1, 3)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 2)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    sol6s = constant_renumber(sol6, 'C', 1, 2)
    sol7s = constant_renumber(sol7, 'C', 1, 2)
    sol8s = constant_renumber(sol8, 'C', 1, 2)
    sol9s = constant_renumber(sol9, 'C', 1, 3)
    sol10s = constant_renumber(sol10, 'C', 1, 2)
    sol11s = constant_renumber(sol11, 'C', 1, 2)
    sol12s = constant_renumber(sol12, 'C', 1, 4)
    assert dsolve(eq1, hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, hint=hint) in (sol2, sol2s)
    assert dsolve(eq3, hint=hint) in (sol3, sol3s)
    assert dsolve(eq4, hint=hint) in (sol4, sol4s)
    assert dsolve(eq5, hint=hint) in (sol5, sol5s)
    assert dsolve(eq6, hint=hint) in (sol6, sol6s)
    assert dsolve(eq7, hint=hint) in (sol7, sol7s)
    assert dsolve(eq8, hint=hint) in (sol8, sol8s)
    assert dsolve(eq9, hint=hint) in (sol9, sol9s)
    assert dsolve(eq10, hint=hint) in (sol10, sol10s)
    assert dsolve(eq11, hint=hint + '_Integral') in (sol11, sol11s)
    assert dsolve(eq12, hint=hint) in (sol12, sol12s)
    assert checkodesol(eq1, sol1, order=3, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=3, solve_for_func=False)[0]
    assert checkodesol(eq3, sol3, order=1, solve_for_func=False)[0]
    assert checkodesol(eq4, sol4, order=2, solve_for_func=False)[0]
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    assert checkodesol(eq6, sol6, order=2, solve_for_func=False)[0]
    assert checkodesol(eq7, sol7, order=2, solve_for_func=False)[0]
    assert checkodesol(eq8, sol8, order=2, solve_for_func=False)[0]
    assert checkodesol(eq9, sol9, order=3, solve_for_func=False)[0]
    assert checkodesol(eq10, sol10, order=2, solve_for_func=False)[0]
    assert checkodesol(eq12, sol12, order=4, solve_for_func=False)[0]


def test_nth_linear_constant_coeff_variation_of_parameters_coverage():
    hint = 'nth_linear_constant_coeff_variation_of_parameters'
    eq = (f(x).diff((x, 3)) - 3*f(x).diff((x, 2)) + 3*f(x).diff(x) -
          f(x) - exp(x)*log(x))
    sol = Eq(f(x), exp(x)*(C1 + C2*x + C3*x**2 + x**3*(6*log(x) - 11)/36))
    assert dsolve(eq, f(x), hint=hint) == sol


def test_nth_linear_constant_coeff_variation_of_parameters_simplify_False():
    # solve_variation_of_parameters shouldn't attempt to simplify the
    # Wronskian if simplify=False.  If wronskian() ever gets good enough
    # to simplify the result itself, this test might fail.
    hint = 'nth_linear_constant_coeff_variation_of_parameters'
    assert dsolve(f(x).diff((x, 5)) + 2*f(x).diff((x, 3)) + f(x).diff(x) -
                  2*x - exp(I*x), f(x), hint + '_Integral', simplify=False) != \
        dsolve(f(x).diff((x, 5)) + 2*f(x).diff((x, 3)) + f(x).diff(x) -
               2*x - exp(I*x), f(x), hint + '_Integral', simplify=True)


def test_Liouville_ODE():
    hint = 'Liouville'
    # The first part here used to be test_ODE_1() from test_solvers.py
    eq1 = diff(f(x), x)/x + diff(f(x), x, x)/2 - diff(f(x), x)**2/2
    eq1a = diff(x*exp(-f(x)), x, x)
    # compare to test_unexpanded_Liouville_ODE() below
    eq2 = (eq1*exp(-f(x))/exp(f(x))).expand()
    eq3 = diff(f(x), x, x) + 1/f(x)*(diff(f(x), x))**2 + 1/x*diff(f(x), x)
    eq4 = x*diff(f(x), x, x) + x/f(x)*diff(f(x), x)**2 + x*diff(f(x), x)
    eq5 = Eq((x*exp(f(x))).diff(x, x), 0)
    sol1 = Eq(f(x), log(x/(C1 + C2*x)))
    sol1a = Eq(C1 + C2/x - exp(-f(x)), 0)
    sol2 = sol1
    sol3 = {Eq(f(x), -sqrt(C1 + C2*log(x))),
            Eq(f(x), sqrt(C1 + C2*log(x)))}
    sol4 = {Eq(f(x), sqrt(C1 + C2*exp(x))*exp(-x/2)),
            Eq(f(x), -sqrt(C1 + C2*exp(x))*exp(-x/2))}
    sol5 = Eq(f(x), log(C1 + C2/x))
    sol1s = constant_renumber(sol1, 'C', 1, 2)
    sol2s = constant_renumber(sol2, 'C', 1, 2)
    sol3s = constant_renumber(sol3, 'C', 1, 2)
    sol4s = constant_renumber(sol4, 'C', 1, 2)
    sol5s = constant_renumber(sol5, 'C', 1, 2)
    assert dsolve(eq1, hint=hint) in (sol1, sol1s)
    assert dsolve(eq1a, hint=hint) in (sol1, sol1s)
    assert dsolve(eq2, hint=hint) in (sol2, sol2s)
    assert set(dsolve(eq3, hint=hint)) in (sol3, sol3s)
    assert set(dsolve(eq4, hint=hint)) in (sol4, sol4s)
    assert dsolve(eq5, hint=hint) in (sol5, sol5s)
    assert checkodesol(eq1, sol1, order=2, solve_for_func=False)[0]
    assert checkodesol(eq1a, sol1a, order=2, solve_for_func=False)[0]
    assert checkodesol(eq2, sol2, order=2, solve_for_func=False)[0]
    assert all(i[0] for i in checkodesol(eq3, sol3, order=2,
                                         solve_for_func=False))
    assert all(i[0] for i in checkodesol(eq4, sol4, order=2,
                                         solve_for_func=False))
    assert checkodesol(eq5, sol5, order=2, solve_for_func=False)[0]
    not_Liouville1 = classify_ode(diff(f(x), x)/x + f(x)*diff(f(x), x, x)/2 -
                                  diff(f(x), x)**2/2, f(x))
    not_Liouville2 = classify_ode(diff(f(x), x)/x + diff(f(x), x, x)/2 -
                                  x*diff(f(x), x)**2/2, f(x))
    assert hint not in not_Liouville1
    assert hint not in not_Liouville2
    assert hint + '_Integral' not in not_Liouville1
    assert hint + '_Integral' not in not_Liouville2


def test_unexpanded_Liouville_ODE():
    # This is the same as eq1 from test_Liouville_ODE() above.
    eq1 = diff(f(x), x)/x + diff(f(x), x, x)/2 - diff(f(x), x)**2/2
    eq2 = eq1*exp(-f(x))/exp(f(x))
    sol2 = Eq(f(x), log(x/(C1 + C2*x)))
    sol2s = constant_renumber(sol2, 'C', 1, 2)
    assert dsolve(eq2) in (sol2, sol2s)
    assert checkodesol(eq2, sol2, order=2, solve_for_func=False)[0]


def test_sympyissue_4785():
    eq = x + A*(x + diff(f(x), x) + f(x)) + diff(f(x), x) + f(x) + 2
    assert classify_ode(eq, f(x)) == ('1st_linear', 'almost_linear',
                                      '1st_power_series', 'lie_group',
                                      'nth_linear_constant_coeff_undetermined_coefficients',
                                      'nth_linear_constant_coeff_variation_of_parameters',
                                      '1st_linear_Integral', 'almost_linear_Integral',
                                      'nth_linear_constant_coeff_variation_of_parameters_Integral')
    # issue sympy/sympy#4864
    eq = (x**2 + f(x)**2)*f(x).diff(x) - 2*x*f(x)
    assert classify_ode(eq, f(x)) == ('1st_exact',
                                      '1st_homogeneous_coeff_best',
                                      '1st_homogeneous_coeff_subs_indep_div_dep',
                                      '1st_homogeneous_coeff_subs_dep_div_indep',
                                      '1st_power_series',
                                      'lie_group', '1st_exact_Integral',
                                      '1st_homogeneous_coeff_subs_indep_div_dep_Integral',
                                      '1st_homogeneous_coeff_subs_dep_div_indep_Integral')


def test_constant_renumber_order_sympyissue_5308():
    assert constant_renumber(C1*x + C2*y, 'C', 1, 2) == \
        constant_renumber(C1*y + C2*x, 'C', 1, 2) == \
        C1*x + C2*y
    e = C1*(C2 + x)*(C3 + y)
    for a, b, c in variations([C1, C2, C3], 3):
        assert constant_renumber(a*(b + x)*(c + y), 'C', 1, 3) == e


def test_sympyissue_5112_5430():
    assert homogeneous_order(-log(x) + acosh(x), x) is None
    assert homogeneous_order(y - log(x), x, y) is None


def test_nth_order_linear_euler_eq_homogeneous():
    x, t, a, b, c = symbols('x t a b c')
    y = Function('y')
    our_hint = 'nth_linear_euler_eq_homogeneous'

    eq = diff(f(t), (t, 4))*t**4 - 13*diff(f(t), (t, 2))*t**2 + 36*f(t)
    assert our_hint in classify_ode(eq)

    eq = a*y(t) + b*t*diff(y(t), t) + c*t**2*diff(y(t), (t, 2))
    assert our_hint in classify_ode(eq)

    eq = Eq(-3*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0)
    sol = C1 + C2*x**Rational(5, 2)
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(3*f(x) - 5*diff(f(x), x)*x + 2*x**2*diff(f(x), x, x), 0)
    sol = C1*sqrt(x) + C2*x**3
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(4*f(x) + 5*diff(f(x), x)*x + x**2*diff(f(x), x, x), 0)
    sol = (C1 + C2*log(x))/x**2
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(6*f(x) - 6*diff(f(x), x)*x + 1*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0)
    sol = dsolve(eq, f(x), hint=our_hint)
    sol = C1/x**2 + C2*x + C3*x**3
    sols = constant_renumber(sol, 'C', 1, 4)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(-125*f(x) + 61*diff(f(x), x)*x - 12*x**2*diff(f(x), x, x) + x**3*diff(f(x), x, x, x), 0)
    sol = x**5*(C1 + C2*log(x) + C3*log(x)**2)
    sols = [sol, constant_renumber(sol, 'C', 1, 4)]
    sols += [sols[-1].expand()]
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in sols
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = t**2*diff(y(t), (t, 2)) + t*diff(y(t), t) - 9*y(t)
    sol = C1*t**3 + C2*t**-3
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, y(t), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = x**3*f(x).diff((x, 3)) + f(x)
    sol1 = dsolve(eq, f(x), hint=our_hint)
    _x = sol1.atoms(RootOf).pop().poly.gen
    r0, r1, r2 = (_x**3 - 3*_x**2 + 2*_x + 1).as_poly().all_roots()
    sol0 = Eq(f(x), C1*x**r0 + C2*x**r1 + C3*x**r2)
    assert sol0 == sol1

    eq = x**2*f(x).diff((x, 2)) + f(x)
    sol = Eq(f(x), sqrt(x)*(C1*sin(sqrt(3)*log(x)/2) + C2*cos(sqrt(3)*log(x)/2)))
    assert dsolve(eq, f(x), hint=our_hint) == sol


def test_nth_order_linear_euler_eq_nonhomogeneous_undetermined_coefficients():
    x, t = symbols('x t')
    a, b, c, d = symbols('a b c d', integer=True)
    our_hint = 'nth_linear_euler_eq_nonhomogeneous_undetermined_coefficients'

    eq = x**4*diff(f(x), (x, 4)) - 13*x**2*diff(f(x), (x, 2)) + 36*f(x) + x
    assert our_hint in classify_ode(eq, f(x))

    eq = a*x**2*diff(f(x), (x, 2)) + b*x*diff(f(x), x) + c*f(x) + d*log(x)
    assert our_hint in classify_ode(eq, f(x))

    eq = Eq(x**2*diff(f(x), x, x) + x*diff(f(x), x), 1)
    sol = C1 + C2*log(x) + log(x)**2/2
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq, f(x))
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(x**2*diff(f(x), x, x) - 2*x*diff(f(x), x) + 2*f(x), x**3)
    sol = x*(C1 + C2*x + Rational(1, 2)*x**2)
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq, f(x))
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(x**2*diff(f(x), x, x) - x*diff(f(x), x) - 3*f(x), log(x)/x)
    sol = C1/x + C2*x**3 - Rational(1, 16)*log(x)/x - Rational(1, 8)*log(x)**2/x
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq, f(x))
    assert dsolve(eq, f(x), hint=our_hint).rhs.expand() in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(x**2*diff(f(x), x, x) + 3*x*diff(f(x), x) - 8*f(x), log(x)**3 - log(x))
    sol = C1/x**4 + C2*x**2 - Rational(1, 8)*log(x)**3 - Rational(3, 32)*log(x)**2 - Rational(1, 64)*log(x) - Rational(7, 256)
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs.expand() in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(x**3*diff(f(x), x, x, x) - 3*x**2*diff(f(x), x, x) + 6*x*diff(f(x), x) - 6*f(x), log(x))
    sol = C1*x + C2*x**2 + C3*x**3 - Rational(1, 6)*log(x) - Rational(11, 36)
    sols = constant_renumber(sol, 'C', 1, 3)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs.expand() in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]


def test_nth_order_linear_euler_eq_nonhomogeneous_variation_of_parameters():
    x, t = symbols('x, t')
    a, b, c, d = symbols('a, b, c, d', integer=True)
    our_hint = 'nth_linear_euler_eq_nonhomogeneous_variation_of_parameters'

    eq = Eq(x**2*diff(f(x), (x, 2)) - 8*x*diff(f(x), x) + 12*f(x), x**2)
    assert our_hint in classify_ode(eq, f(x))

    eq = Eq(a*x**3*diff(f(x), (x, 3)) + b*x**2*diff(f(x), (x, 2)) + c*x*diff(f(x), x) + d*f(x), x*log(x))
    assert our_hint in classify_ode(eq, f(x))

    eq = Eq(x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x), x**4)
    sol = C1*x + C2*x**2 + x**4/6
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs.expand() in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(3*x**2*diff(f(x), x, x) + 6*x*diff(f(x), x) - 6*f(x), x**3*exp(x))
    sol = C1/x**2 + C2*x + x*exp(x)/3 - 4*exp(x)/3 + 8*exp(x)/(3*x) - 8*exp(x)/(3*x**2)
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs.expand() in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = Eq(x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x), x**4*exp(x))
    sol = C1*x + C2*x**2 + x**2*exp(x) - 2*x*exp(x)
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs.expand() in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = x**2*Derivative(f(x), x, x) - 2*x*Derivative(f(x), x) + 2*f(x) - log(x)
    sol = C1*x + C2*x**2 + log(x)/2 + 3/4
    sols = constant_renumber(sol, 'C', 1, 2)
    assert our_hint in classify_ode(eq)
    assert dsolve(eq, f(x), hint=our_hint).rhs in (sol, sols)
    assert checkodesol(eq, sol, order=2, solve_for_func=False)[0]

    eq = x**2*f(x).diff((x, 2)) + f(x) - 1
    sol = Eq(f(x), sqrt(x)*(C1*sin(sqrt(3)*log(x)/2) + C2*cos(sqrt(3)*log(x)/2) + sin(sqrt(3)*log(x)/2)*Integral(2*sqrt(3)*cos(sqrt(3)*log(x)/2)/(3*x**Rational(3, 2)), x) - cos(sqrt(3)*log(x)/2)*Integral(2*sqrt(3)*sin(sqrt(3)*log(x)/2)/(3*x**Rational(3, 2)), x)))
    assert dsolve(eq, f(x), hint=our_hint + '_Integral') == sol


def test_sympyissue_5095():
    f = Function('f')
    pytest.raises(ValueError, lambda: dsolve(f(x).diff(x)**2, f(x), 'separable'))
    pytest.raises(ValueError, lambda: dsolve(f(x).diff(x)**2, f(x), 'fdsjf'))


def test_almost_linear():
    A = Symbol('A', positive=True)
    our_hint = 'almost_linear'
    f = Function('f')
    d = f(x).diff(x)
    eq = x**2*f(x)**2*d + f(x)**3 + 1
    sol = dsolve(eq, f(x), hint='almost_linear')
    assert sol[0].rhs == cbrt(C1*exp(3/x) - 1)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x*f(x)*d + 2*x*f(x)**2 + 1
    sol = dsolve(eq, f(x), hint='almost_linear')
    assert sol[0].rhs == -sqrt(C1 - 2*Ei(4*x))*exp(-2*x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x*d + x*f(x) + 1
    sol = dsolve(eq, f(x), hint='almost_linear')
    assert sol.rhs == (C1 - Ei(x))*exp(-x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]
    assert our_hint in classify_ode(eq, f(x))

    eq = x*exp(f(x))*d + exp(f(x)) + 3*x
    sol = dsolve(eq, f(x), hint='almost_linear')
    assert sol.rhs == log(C1/x - 3*x/2)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x + A*(x + diff(f(x), x) + f(x)) + diff(f(x), x) + f(x) + 2
    sol = dsolve(eq, f(x), hint='almost_linear')
    assert sol.rhs == (C1 + Piecewise(
        (x, Eq(A + 1, 0)), ((-A*x + A - x - 1)*exp(x)/(A + 1), True)))*exp(-x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_exact_enhancement():
    f = Function('f')(x)
    df = Derivative(f, x)
    eq = f/x**2 + ((f*x - 1)/x)*df
    sol = dsolve(eq, f)
    assert sol == [Eq(f, -sqrt(C1*x**2 + 1)/abs(x) + 1/x),
                   Eq(f, sqrt(C1*x**2 + 1)/abs(x) + 1/x)]

    eq = (x*f - 1) + df*(x**2 - x*f)
    rhs = [sol.rhs for sol in dsolve(eq, f)]
    assert rhs[0] == x - sqrt(C1 + x**2 - 2*log(x))
    assert rhs[1] == x + sqrt(C1 + x**2 - 2*log(x))

    eq = (x + 2)*sin(f) + df*x*cos(f)
    rhs = [sol.rhs for sol in dsolve(eq, f)]
    assert rhs == [
        -acos(-sqrt(C1*exp(-2*x)/x**4 + 1)) + 2*pi,
        -acos(sqrt(C1*exp(-2*x)/x**4 + 1)) + 2*pi,
        acos(-sqrt(C1*exp(-2*x)/x**4 + 1)),
        acos(sqrt(C1*exp(-2*x)/x**4 + 1))]


def test_separable_reduced():
    f = Function('f')
    x = Symbol('x')  # BUG: if x is real, a more complex solution is returned!
    df = f(x).diff(x)
    eq = (x / f(x))*df + tan(x**2*f(x) / (x**2*f(x) - 1))
    assert classify_ode(eq) == ('separable_reduced', 'lie_group',
                                'separable_reduced_Integral')

    eq = x * df + f(x) * (1 / (x**2*f(x) - 1))
    assert classify_ode(eq) == ('separable_reduced', 'lie_group',
                                'separable_reduced_Integral')
    sol = dsolve(eq, hint='separable_reduced', simplify=False)
    assert sol.lhs == log(x**2*f(x))/3 + log(x**2*f(x) - Rational(3, 2))/6
    assert sol.rhs == C1 + log(x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    # this is the equation that does not like x to be real
    eq = f(x).diff(x) + (f(x) / (x**4*f(x) - x))
    assert classify_ode(eq) == ('separable_reduced', 'lie_group',
                                'separable_reduced_Integral')
    # generates PolynomialError in solve attempt
    sol = dsolve(eq, hint='separable_reduced')
    assert sol.lhs - sol.rhs == \
        log(x**3*f(x))/4 + log(x**3*f(x) - Rational(4, 3))/12 - C1 - log(x)
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]

    eq = x*df + f(x)*(x**2*f(x))
    sol = dsolve(eq, hint='separable_reduced', simplify=False)
    assert sol == Eq(log(x**2*f(x))/2 - log(x**2*f(x) - 2)/2, C1 + log(x))
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_homogeneous_function():
    f = Function('f')
    eq1 = tan(x + f(x))
    eq2 = sin((3*x)/(4*f(x)))
    eq3 = cos(3*x/4*f(x))
    eq4 = log((3*x + 4*f(x))/(5*f(x) + 7*x))
    eq5 = exp((2*x**2)/(3*f(x)**2))
    eq6 = log((3*x + 4*f(x))/(5*f(x) + 7*x) + exp((2*x**2)/(3*f(x)**2)))
    eq7 = sin((3*x)/(5*f(x) + x**2))
    assert homogeneous_order(eq1, x, f(x)) is None
    assert homogeneous_order(eq2, x, f(x)) == 0
    assert homogeneous_order(eq3, x, f(x)) is None
    assert homogeneous_order(eq4, x, f(x)) == 0
    assert homogeneous_order(eq5, x, f(x)) == 0
    assert homogeneous_order(eq6, x, f(x)) == 0
    assert homogeneous_order(eq7, x, f(x)) is None


def test_linear_coeff_match():
    n, d = z*(2*x + 3*f(x) + 5), z*(7*x + 9*f(x) + 11)
    rat = n/d
    eq1 = sin(rat) + cos(rat.expand())
    eq2 = rat
    eq3 = log(sin(rat))
    ans = (4, -Rational(13, 3))
    assert _linear_coeff_match(eq1, f(x)) == ans
    assert _linear_coeff_match(eq2, f(x)) == ans
    assert _linear_coeff_match(eq3, f(x)) == ans

    # no c
    eq4 = (3*x)/f(x)
    # not x and f(x)
    eq5 = (3*x + 2)/x
    # denom will be zero
    eq6 = (3*x + 2*f(x) + 1)/(3*x + 2*f(x) + 5)
    # not rational coefficient
    eq7 = (3*x + 2*f(x) + sqrt(2))/(3*x + 2*f(x) + 5)
    assert _linear_coeff_match(eq4, f(x)) is None
    assert _linear_coeff_match(eq5, f(x)) is None
    assert _linear_coeff_match(eq6, f(x)) is None
    assert _linear_coeff_match(eq7, f(x)) is None


def test_linear_coefficients():
    f = Function('f')
    sol = Eq(f(x), (C1 - 3*x**2 - 18*x)/Mul(2, x**2 + 6*x + 9, evaluate=False))
    eq = f(x).diff(x) + (3 + 2*f(x))/(x + 3)
    assert dsolve(eq, hint='linear_coefficients') == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_constantsimp_take_problem():
    c = exp(C1) + 2
    assert len((constantsimp(exp(C1) + c + c*x, [C1])).as_poly().gens) == 2


def test_sympyissue_6879():
    f = Function('f')
    eq = Eq(Derivative(f(x), (x, 2)) - 2*Derivative(f(x), x) + f(x), sin(x))
    sol = (C1 + C2*x)*exp(x) + cos(x)/2
    assert dsolve(eq).rhs == sol
    assert checkodesol(eq, sol, order=1, solve_for_func=False)[0]


def test_sympyissue_6989():
    f = Function('f')
    k = Symbol('k')
    assert dsolve(f(x).diff(x) - x*exp(-k*x), f(x)) == \
        Eq(f(x), C1 + Piecewise(
            (x**2/2, Eq(k**3, 0)),
            ((-k**2*x - k)*exp(-k*x)/k**3, True)
        ))
    eq = -f(x).diff(x) + x*exp(-k*x)
    sol = dsolve(eq, f(x))
    actual_sol = Eq(f(x), Piecewise((C1 + x**2/2, Eq(k**3, 0)),
                                    (C1 - x*exp(-k*x)/k - exp(-k*x)/k**2, True)
                                    ))
    errstr = str(eq) + ' : ' + str(sol) + ' == ' + str(actual_sol)
    assert sol == actual_sol, errstr


def test_heuristic1():
    y, a, b, c, a4, a3, a2, a1, a0 = symbols('y a b c a4 a3 a2 a1 a0')
    f = Function('f')
    xi = Function('xi')
    eta = Function('eta')
    df = f(x).diff(x)
    eq = Eq(df, x**2*f(x))
    eq1 = f(x).diff(x) + a*f(x) - c*exp(b*x)
    eq2 = f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)
    eq3 = (1 + 2*x)*df + 2 - 4*exp(-f(x))
    eq4 = f(x).diff(x) - (a4*x**4 + a3*x**3 + a2*x**2 + a1*x + a0)**Rational(-1, 2)
    eq5 = x**2*df - f(x) + x**2*exp(x - (1/x))
    eqlist = [eq, eq1, eq2, eq3, eq4, eq5]

    i = infinitesimals(eq, hint='abaco1_simple')
    assert i == [{eta(x, f(x)): exp(x**3/3), xi(x, f(x)): 0},
                 {eta(x, f(x)): f(x), xi(x, f(x)): 0},
                 {eta(x, f(x)): 0, xi(x, f(x)): x**(-2)}]
    i1 = infinitesimals(eq1, hint='abaco1_simple')
    assert i1 == [{eta(x, f(x)): exp(-a*x), xi(x, f(x)): 0}]
    i2 = infinitesimals(eq2, hint='abaco1_simple')
    assert i2 == [{eta(x, f(x)): exp(-x**2), xi(x, f(x)): 0}]
    i3 = infinitesimals(eq3, hint='abaco1_simple')
    assert i3 == [{eta(x, f(x)): 0, xi(x, f(x)): 2*x + 1},
                  {eta(x, f(x)): 0, xi(x, f(x)): 1/(exp(f(x)) - 2)}]
    i4 = infinitesimals(eq4, hint='abaco1_simple')
    assert i4 == [{eta(x, f(x)): 1, xi(x, f(x)): 0},
                  {eta(x, f(x)): 0,
                   xi(x, f(x)): sqrt(a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)}]
    i5 = infinitesimals(eq5, hint='abaco1_simple')
    assert i5 == [{xi(x, f(x)): 0, eta(x, f(x)): exp(-1/x)}]

    ilist = [i, i1, i2, i3, i4, i5]
    for eq, i in (zip(eqlist, ilist)):
        check = checkinfsol(eq, i)
        assert check[0]

    eq = f(x).diff(x) - x**2*f(x)
    assert infinitesimals(eq) == [{eta(x, f(x)): exp(x**3/3), xi(x, f(x)): 0}]


def test_sympyissue_6247():
    eq = x**2*f(x)**2 + x*Derivative(f(x), x)
    sol = dsolve(eq, hint='separable_reduced')
    assert checkodesol(eq, sol, order=1)[0]
    eq = f(x).diff(x, x) + 4*f(x)
    sol = dsolve(eq, f(x), simplify=False)
    assert sol == Eq(f(x), C1*sin(2*x) + C2*cos(2*x))


def test_heuristic2():
    xi = Function('xi')
    eta = Function('eta')
    df = f(x).diff(x)

    # This ODE can be solved by the Lie Group method, when there are
    # better assumptions
    eq = df - (f(x)/x)*(x*log(x**2/f(x)) + 2)
    i = infinitesimals(eq, hint='abaco1_product')
    assert i == [{eta(x, f(x)): f(x)*exp(-x), xi(x, f(x)): 0}]
    assert checkinfsol(eq, i)[0]


def test_heuristic3():
    xi = Function('xi')
    eta = Function('eta')
    a, b = symbols('a b')
    df = f(x).diff(x)

    eq = x**2*df + x*f(x) + f(x)**2 + x**2
    i = infinitesimals(eq, hint='bivariate')
    assert i == [{eta(x, f(x)): f(x), xi(x, f(x)): x}]
    assert checkinfsol(eq, i)[0]

    eq = x**2*(-f(x)**2 + df) - a*x**2*f(x) + 2 - a*x
    i = infinitesimals(eq, hint='bivariate')
    assert checkinfsol(eq, i)[0]


def test_heuristic_4():
    y, a = symbols('y a')

    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    i = infinitesimals(eq, hint='chi')
    assert checkinfsol(eq, i)[0]


def test_heuristic_function_sum():
    xi = Function('xi')
    eta = Function('eta')
    eq = f(x).diff(x) - (3*(1 + x**2/f(x)**2)*atan(f(x)/x) + (1 - 2*f(x))/x +
                         (1 - 3*f(x))*(x/f(x)**2))
    i = infinitesimals(eq, hint='function_sum')
    assert i == [{eta(x, f(x)): f(x)**(-2) + x**(-2), xi(x, f(x)): 0}]
    assert checkinfsol(eq, i)[0]


def test_heuristic_abaco2_similar():
    xi = Function('xi')
    eta = Function('eta')
    F = Function('F')
    a, b = symbols('a b')
    eq = f(x).diff(x) - F(a*x + b*f(x))
    i = infinitesimals(eq, hint='abaco2_similar')
    assert i == [{eta(x, f(x)): -a/b, xi(x, f(x)): 1}]
    assert checkinfsol(eq, i)[0]

    eq = f(x).diff(x) - (f(x)**2 / (sin(f(x) - x) - x**2 + 2*x*f(x)))
    i = infinitesimals(eq, hint='abaco2_similar')
    assert i == [{eta(x, f(x)): f(x)**2, xi(x, f(x)): f(x)**2}]
    assert checkinfsol(eq, i)[0]


def test_heuristic_abaco2_unique_unknown():
    xi = Function('xi')
    eta = Function('eta')
    F = Function('F')
    a, b = symbols('a b')
    x = Symbol('x', positive=True)

    eq = f(x).diff(x) - x**(a - 1)*(f(x)**(1 - b))*F(x**a/a + f(x)**b/b)
    i = infinitesimals(eq, hint='abaco2_unique_unknown')
    assert i == [{eta(x, f(x)): -f(x)*f(x)**(-b), xi(x, f(x)): x*x**(-a)}]
    assert checkinfsol(eq, i)[0]

    eq = f(x).diff(x) + tan(F(x**2 + f(x)**2) + atan(x/f(x)))
    i = infinitesimals(eq, hint='abaco2_unique_unknown')
    assert i == [{eta(x, f(x)): x, xi(x, f(x)): -f(x)}]
    assert checkinfsol(eq, i)[0]

    eq = (x*f(x).diff(x) + f(x) + 2*x)**2 - 4*x*f(x) - 4*x**2 - 4*a
    i = infinitesimals(eq, hint='abaco2_unique_unknown')
    assert checkinfsol(eq, i)[0]


def test_heuristic_linear():
    a, b, m, n = symbols('a b m n')

    eq = x**(n*(m + 1) - m)*(f(x).diff(x)) - a*f(x)**n - b*x**(n*(m + 1))
    i = infinitesimals(eq, hint='linear')
    assert checkinfsol(eq, i)[0]


@pytest.mark.xfail
def test_kamke():
    a, b, alpha, c = symbols('a b alpha c')
    eq = x**2*(a*f(x)**2+(f(x).diff(x))) + b*x**alpha + c
    i = infinitesimals(eq, hint='sum_function')
    assert checkinfsol(eq, i)[0]


def test_series():
    C1 = Symbol('C1')
    eq = f(x).diff(x) - f(x)
    assert dsolve(eq, hint='1st_power_series') == Eq(f(x),
                                                     C1 + C1*x + C1*x**2/2 + C1*x**3/6 + C1*x**4/24 +
                                                     C1*x**5/120 + O(x**6))
    eq = f(x).diff(x) - x*f(x)
    assert dsolve(eq, hint='1st_power_series') == Eq(f(x),
                                                     C1*x**4/8 + C1*x**2/2 + C1 + O(x**6))
    eq = f(x).diff(x) - sin(x*f(x))
    sol = Eq(f(x), (x - 2)**2*(1 + sin(4))*cos(4) + (x - 2)*sin(4) + 2 + O(x**3))
    assert dsolve(eq, hint='1st_power_series', init={f(2): 2}, n=3) == sol


def test_lie_group():
    C1 = Symbol('C1')
    x = Symbol('x')  # assuming x is real generates an error!
    a, b, c = symbols('a b c')
    eq = f(x).diff(x)**2
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol)[0]

    eq = Eq(f(x).diff(x), x**2*f(x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), C1*exp(x**3/3))
    assert checkodesol(eq, sol)[0]

    eq = f(x).diff(x) + a*f(x) - c*exp(b*x)
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol)[0]

    eq = f(x).diff(x) + 2*x*f(x) - x*exp(-x**2)
    sol = dsolve(eq, f(x), hint='lie_group')
    actual_sol = Eq(f(x), (C1 + x**2/2)*exp(-x**2))
    errstr = str(eq)+' : '+str(sol)+' == '+str(actual_sol)
    assert sol == actual_sol, errstr
    assert checkodesol(eq, sol)[0]

    eq = (1 + 2*x)*(f(x).diff(x)) + 2 - 4*exp(-f(x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), log(C1/(2*x + 1) + 2))
    assert checkodesol(eq, sol)[0]

    eq = x**2*(f(x).diff(x)) - f(x) + x**2*exp(x - (1/x))
    sol = dsolve(eq, f(x), hint='lie_group')
    assert checkodesol(eq, sol)[0]

    eq = x**2*f(x)**2 + x*Derivative(f(x), x)
    sol = dsolve(eq, f(x), hint='lie_group')
    assert sol == Eq(f(x), 2/(C1 + x**2))
    assert checkodesol(eq, sol)[0]

    # issue diofant/diofant#309
    assert dsolve(f(x).diff(x)**2 - 1, f(x)) == [Eq(f(x), C1 - x),
                                                 Eq(f(x), C1 + x)]


def test_user_infinitesimals():
    x = Symbol('x')  # assuming x is real generates an error
    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    sol = dsolve(eq, hint='lie_group', xi=sqrt(f(x) - 1)/sqrt(f(x) + 1),
                 eta=0)
    actual_sol = Eq(f(x), (C1 + x**2)/(C1 - x**2))
    errstr = str(eq)+' : '+str(sol)+' == '+str(actual_sol)
    assert sol == actual_sol, errstr
    pytest.raises(ValueError, lambda: dsolve(eq, hint='lie_group', xi=0, eta=f(x)))


def test_sympyissue_7081():
    eq = x*(f(x).diff(x)) + 1 - f(x)**2
    assert dsolve(eq) == Eq(f(x), (C1 + x**2)/(C1 - x**2))


def test_2nd_power_series_ordinary():
    C1, C2 = symbols('C1 C2')
    eq = f(x).diff((x, 2)) - x*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x),
                            C2*(x**3/6 + 1) + C1*x*(x**3/12 + 1) + O(x**6))
    assert dsolve(eq, x0=-2) == Eq(f(x),
                                   C2*((x + 2)**4/6 + (x + 2)**3/6 - (x + 2)**2 + 1)
                                   + C1*(x + (x + 2)**4/12 - (x + 2)**3/3 + 2)
                                   + O(x**6))
    assert dsolve(eq, n=2) == Eq(f(x), C2*x + C1 + O(x**2))

    eq = (1 + x**2)*(f(x).diff((x, 2))) + 2*x*(f(x).diff(x)) - 2*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x), C2*(-x**4/3 + x**2 + 1) + C1*x
                            + O(x**6))

    eq = f(x).diff((x, 2)) + x*(f(x).diff(x)) + f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x), C2*(
        x**4/8 - x**2/2 + 1) + C1*x*(-x**2/3 + 1) + O(x**6))
    # for coverage
    assert classify_ode(f(x).diff((x, 2)) +
                        sin(x)*(f(x).diff(x)) + f(x)) == ()

    eq = f(x).diff((x, 2)) + f(x).diff(x) - x*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq) == Eq(f(x), C2*(
        -x**4/24 + x**3/6 + 1) + C1*x*(x**3/24 + x**2/6 - x/2
                                       + 1) + O(x**6))

    eq = f(x).diff((x, 2)) + x*f(x)
    assert classify_ode(eq) == ('2nd_power_series_ordinary',)
    assert dsolve(eq, n=7) == Eq(f(x), C2*(
        x**6/180 - x**3/6 + 1) + C1*x*(-x**3/12 + 1) + O(x**7))


def test_2nd_power_series_regular():
    C1, C2 = symbols('C1 C2')
    eq = x**2*(f(x).diff((x, 2))) - 3*x*(f(x).diff(x)) + (4*x + 4)*f(x)
    assert dsolve(eq) == Eq(f(x), C1*x**2*(-16*x**3/9 +
                                           4*x**2 - 4*x + 1) + O(x**6))

    eq = 4*x**2*(f(x).diff((x, 2))) - 8*x**2*(f(x).diff(x)) + (4*x**2 +
                                                               1)*f(x)
    assert dsolve(eq) == Eq(f(x), C1*sqrt(x)*(
        x**4/24 + x**3/6 + x**2/2 + x + 1) + O(x**6))

    eq = x**2*(f(x).diff((x, 2))) - x**2*(f(x).diff(x)) + (
        x**2 - 2)*f(x)
    assert dsolve(eq) == Eq(f(x), C1*(-x**6/720 - 3*x**5/80 - x**4/8 +
                                      x**2/2 + x/2 + 1)/x + C2*x**2*(-x**3/60 + x**2/20 + x/2 + 1)
                            + O(x**6))

    eq = x**2*(f(x).diff((x, 2))) + x*(f(x).diff(x)) + (x**2 - Rational(1, 4))*f(x)
    assert dsolve(eq) == Eq(f(x), C1*(x**4/24 - x**2/2 + 1)/sqrt(x) +
                            C2*sqrt(x)*(x**4/120 - x**2/6 + 1) + O(x**6))

    eq = x*(f(x).diff((x, 2))) - f(x).diff(x) + 4*x**3*f(x)
    assert dsolve(eq) == Eq(f(x), C2*(-x**4/2 + 1) + C1*x**2 + O(x**6))

    eq = x**3*(f(x).diff((x, 2))) - 3*x*(f(x).diff(x)) + (4*x + 4)*f(x)
    assert '2nd_power_series_regular' not in classify_ode(eq, init={f(0): 1})
    eq = x**3*(f(x).diff((x, 2))) - 3*x*2*(f(x).diff(x)) + 4*f(x)
    assert '2nd_power_series_regular' not in classify_ode(eq, init={f(0): 1})


def test_sympyissue_7093():
    sol = [Eq(f(x), C1 - 2*x*sqrt(x**3)/5), Eq(f(x), C1 + 2*x*sqrt(x**3)/5)]
    eq = Derivative(f(x), x)**2 - x**3
    assert dsolve(eq) == sol
    assert checkodesol(eq, sol) == [(True, 0), (True, 0)]


def test_dsolve_linsystem_symbol():
    eps = Symbol('epsilon', positive=True)
    eq1 = (Eq(diff(f(x), x), -eps*g(x)), Eq(diff(g(x), x), eps*f(x)))
    sol1 = [Eq(f(x), exp(I*eps*x)*I*C2 - exp(-I*eps*x)*I*C1),
            Eq(g(x), exp(I*eps*x)*C2 + exp(-I*eps*x)*C1)]
    s = dsolve(eq1)
    assert checksysodesol(eq1, s) == (True, [0, 0])
    assert s == sol1


def test_C1_function_9239():
    c1, c2 = symbols('C1, C2', cls=Function)
    t = Symbol('t')
    eq = (Eq(diff(c1(t), t), 9*c2(t)), Eq(diff(c2(t), t), 12*c1(t)))
    sol = [Eq(c1(t), sqrt(3)*exp(6*sqrt(3)*t)*C4/2 - sqrt(3)*exp(-6*sqrt(3)*t)*C3/2),
           Eq(c2(t), exp(6*sqrt(3)*t)*C4 + exp(-6*sqrt(3)*t)*C3)]
    assert dsolve(eq) == sol


def test_sympyissue_11290():
    eq = cos(f(x)) - (x*sin(f(x)) - f(x)**2)*f(x).diff(x)
    s0 = dsolve(eq, f(x), simplify=False, hint='1st_exact')
    s1 = dsolve(eq, f(x), simplify=False, hint='1st_exact_Integral')
    assert (str(s1) ==
            'Eq(Subs(Integral(-x*sin(_y) + _y**2 '
            '- Integral(-sin(_y), x), _y) + '
            'Integral(cos(_y), x), (_y, f(x))), C1)')
    assert s1.doit() == s0


def test_sympyissue_7138():
    eqs = [Eq(f(x).diff(x), f(x) - 1), Eq(g(x).diff(x), f(x) + 2*g(x) - 3)]
    assert dsolve(eqs) == [Eq(f(x), -exp(x)*C1 + 1),
                           Eq(g(x), exp(2*x)*C2 + exp(x)*C1 + 1)]


def test_sympyissue_10379():
    t, y = symbols('t,y')
    sol = dsolve(f(t).diff(t) - (1 - 51.05*y*f(t)))
    ans = Eq(f(t), (0.019588638589618 + 0.019588638589618*E**(-1.0*y*(-1.0*C1 + 51.05*t)))/y)
    assert str(sol) == str(ans)


def test_sympyissue_10867():
    v = Eq(g(x).diff(x).diff(x), (x-2)**2 + (x-3)**3)
    ans = Eq(g(x), C1 + C2*x + x**5/20 - 2*x**4/3 + 23*x**3/6 - 23*x**2/2)
    assert dsolve(v, g(x)) == ans


def test_sympyissue_15407():
    eqs = [Eq(Derivative(f(t), t), -(x + y)*f(t)),
           Eq(Derivative(g(t), t), x*f(t)), Eq(Derivative(h(t), t), y*f(t))]
    ans = [Eq(f(t), -exp(-t*(x + y))*C3*(x + y)/y),
           Eq(g(t), C1 + exp(-t*(x + y))*C3*x/y),
           Eq(h(t), C2 + exp(-t*(x + y))*C3)]
    assert dsolve(eqs) == ans


def test_sympyissue_15311():
    eqn = sqrt(2) * f(x).diff((x, 3)) + f(x).diff(x)
    assert dsolve(eqn) == Eq(f(x), C1 + C2*sin(root(8, 4)*x/2) + C3*cos(root(8, 4)*x/2))


def test_sympyissue_15474():
    a, b = symbols('a b')
    eqs = [Eq(f(t).diff(t), a*f(t)), Eq(g(t).diff(t), b*g(t))]
    ans = [Eq(f(t), C1*exp(a*t)), Eq(g(t), C2*exp(b*t))]
    assert dsolve(eqs) == ans
    eqs = [Eq(f(t).diff(t), -a*g(t)), Eq(g(t).diff(t), a*f(t))]
    ans = [Eq(f(t), exp(I*a*t)*I*C2 - exp(-I*a*t)*I*C1),
           Eq(g(t), exp(I*a*t)*C2 + exp(-I*a*t)*C1)]
    assert dsolve(eqs) == ans


def test_sympyissue_15574():
    f1, f2, f3, f4 = symbols('f1 f2 f3 f4', cls=Function)
    eqs = [Eq(f(x).diff(x), f(x)) for f in (f1, f2, f3, f4)]
    ans = [Eq(f(x), c*exp(x)) for f, c in zip((f1, f2, f3, f4), (C1, C2, C3, C4))]
    assert dsolve(eqs[:2]) == ans[:2]
    assert dsolve(eqs[:3]) == ans[:3]
    assert dsolve(eqs[:4]) == ans[:4]


def test__lie_group_remove():
    eq = x**2*y
    assert _lie_group_remove(eq) == x**2*y
    eq = f(x**2*y)
    assert _lie_group_remove(eq) == x**2*y
    eq = y**2*x + f(x**3)
    assert _lie_group_remove(eq) == x*y**2
    eq = (f(x**3) + y)*x**4
    assert _lie_group_remove(eq) == x**4*y


def test_ode_sol_simplicity():
    eq1 = Eq(f(x)/tan(f(x)/(2*x)), C1)
    eq2 = Eq(f(x)/tan(f(x)/(2*x) + f(x)), C2)
    assert [ode_sol_simplicity(eq, f(x)) for eq in [eq1, eq2]] == [28, 35]


def test_get_numbered_constants():
    pytest.raises(ValueError, lambda: get_numbered_constants('spam'))


def test_dsolve_interface():
    pytest.raises(ValueError, lambda: dsolve((f(t).diff(t) - g(t), g(t).diff((t, 2)) + t)))
    pytest.raises(ValueError, lambda: dsolve((f(t).diff(t) - g(t),)))
    pytest.raises(ValueError, lambda: classify_sysode((f(t).diff(t) - g(t),)))
    pytest.raises(ValueError, lambda: classify_sysode((f(t).diff(t) - g(t), g(t).diff(t) + 1), (f(t),)))
    pytest.raises(ValueError, lambda: classify_sysode((f(t).diff(t) - g(t), g(t).diff(t) + 1), (f(t), g(t, x))))
    pytest.raises(ValueError, lambda: ode_nth_linear_euler_eq_homogeneous(2*x**2*f(x).diff((x, 2)) - 3*x*f(x).diff(x),
                                                                          f(x), 2,
                                                                          match={-1: Integer(0), 0: Integer(0),
                                                                                 1: -3*x, 2: 2*x**2},
                                                                          returns='spam'))

    # issue diofant/diofant#293
    eqs = (-f(x) + Derivative(f(x), x) + Derivative(g(x), x),
           g(x) + Derivative(f(x), x) - Derivative(g(x), x))
    sol = [Eq(f(x), -exp(x*(1 - I)/2)*I*C1 + exp(x*(1 + I)/2)*I*C2),
           Eq(g(x), +exp(x*(1 - I)/2)*C1 + exp(x*(1 + I)/2)*C2)]
    assert dsolve(eqs) == sol


def test_odesimp():
    sol = dsolve(x*f(x).diff(x) - f(x) - x*sin(f(x)/x), f(x),
                 hint='1st_homogeneous_coeff_subs_indep_div_dep_Integral',
                 simplify=False)
    sol1 = odesimp(sol, f(x), 1, {C1}, hint='1st_homogeneous_coeff_subs_indep_div_dep')
    ssol = Eq(f(x), 2*x*atan(C1*x))
    assert sol1 == ssol
    sol2 = odesimp(Eq(sol.rhs, sol.lhs), f(x), 1, {C1}, hint='1st_homogeneous_coeff_subs_indep_div_dep')
    assert sol2 == ssol
    pytest.raises(TypeError, lambda: odesimp(sol.lhs - sol.rhs,
                                             f(x), 1, {C1},
                                             hint='1st_homogeneous_coeff_subs_indep_div_dep'))


def test_sympyissue_16635():
    eqs = [Eq(f(t).diff(t), f(t) - g(t) + 15*t - 10),
           Eq(g(t).diff(t), f(t) - g(t) - 15*t - 5)]
    sol = [Eq(f(t), -2*C1 + C2*(-2*t - 1) - 10*t**3 + 5*t**2/2 - 15*t/2 - (-2*t - 1)*(15*t**2/2 - 5*t/2)),
           Eq(g(t), -2*C1 + C2*(-2*t + 1) - 10*t**3 + 5*t**2/2 - 15*t/2 - (-2*t + 1)*(15*t**2/2 - 5*t/2))]
    assert dsolve(eqs) == sol


def test_sympyissue_14312():
    a, b = symbols('a b')
    eqs1 = (Eq(f(t).diff(t), b*g(t)), Eq(g(t).diff(t), -(b + a)*g(t)),
            Eq(h(t).diff(t), a*g(t)))
    sol1 = dsolve(eqs1)
    assert checksysodesol(eqs1, sol1) == (True, [0, 0, 0])

    eqs2 = (Eq(h(t).diff(t), a*g(t)), Eq(f(t).diff(t), b*g(t)),
            Eq(g(t).diff(t), -(b + a)*g(t)))
    sol2 = dsolve(eqs2)
    assert checksysodesol(eqs2, sol2) == (True, [0, 0, 0])


def test_sympyissue_8859():
    eqs = [Eq(f(t).diff(t), f(t) + 3*t), Eq(g(t).diff(t), g(t))]
    sol = [Eq(f(t), -3 - 3*t + C1*exp(t)), Eq(g(t), C2*exp(t))]
    assert dsolve(eqs) == sol


def test_sympyissue_9204():
    m, q = symbols('m q', real=True)
    El = Matrix(symbols('e1:4', real=True))
    B = Matrix(symbols('b1:4', real=True))

    V = Matrix([f(t), g(t), h(t)])
    F = m*V.diff(t) - q*El - q*V.cross(B)

    assert all(_.has(*El) for _ in dsolve(F))

    El[0] = El[1] = 0
    B[0] = B[1] = 0
    e3 = El[2]
    b3 = B[2]
    F = m*V.diff(t) - q*El - q*V.cross(B)

    sol2 = [Eq(f(t), -exp(I*b3*q*t/m)*I*C3 + exp(-I*b3*q*t/m)*I*C2),
            Eq(g(t), exp(I*b3*q*t/m)*C3 + exp(-I*b3*q*t/m)*C2),
            Eq(h(t), C1 + e3*q*t/m)]

    assert dsolve(F) == sol2
