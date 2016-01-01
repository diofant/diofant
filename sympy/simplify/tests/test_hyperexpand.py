from random import randrange

import pytest

from sympy.simplify.hyperexpand import (ShiftA, ShiftB, UnShiftA, UnShiftB,
                       MeijerShiftA, MeijerShiftB, MeijerShiftC, MeijerShiftD,
                       MeijerUnShiftA, MeijerUnShiftB, MeijerUnShiftC,
                       MeijerUnShiftD,
                       ReduceOrder, reduce_order, apply_operators,
                       devise_plan, make_derivative_operator, Formula,
                       hyperexpand, Hyper_Function, G_Function,
                       reduce_order_meijer,
                       build_hypergeometric_formula)
from sympy import hyper, I, S, meijerg, Piecewise
from sympy.abc import z, a, b, c
from sympy.utilities.randtest import verify_numerically as tn
from sympy import (cos, sin, log, exp, asin, lowergamma, atanh, besseli,
                   gamma, sqrt, pi, erf, exp_polar)


def test_branch_bug():
    assert hyperexpand(hyper((-S(1)/3, S(1)/2), (S(2)/3, S(3)/2), -z)) == \
        -z**S('1/3')*lowergamma(exp_polar(I*pi)/3, z)/5 \
        + sqrt(pi)*erf(sqrt(z))/(5*sqrt(z))
    assert hyperexpand(meijerg([S(7)/6, 1], [], [S(2)/3], [S(1)/6, 0], z)) == \
        2*z**S('2/3')*(2*sqrt(pi)*erf(sqrt(z))/sqrt(z) - 2*lowergamma(
                       S(2)/3, z)/z**S('2/3'))*gamma(S(2)/3)/gamma(S(5)/3)


def test_hyperexpand():
    # Luke, Y. L. (1969), The Special Functions and Their Approximations,
    # Volume 1, section 6.2

    assert hyperexpand(hyper([], [], z)) == exp(z)
    assert hyperexpand(hyper([1, 1], [2], -z)*z) == log(1 + z)
    assert hyperexpand(hyper([], [S.Half], -z**2/4)) == cos(z)
    assert hyperexpand(z*hyper([], [S('3/2')], -z**2/4)) == sin(z)
    assert hyperexpand(hyper([S('1/2'), S('1/2')], [S('3/2')], z**2)*z) \
        == asin(z)


def can_do(ap, bq, numerical=True, div=1, lowerplane=False):
    from sympy import exp_polar, exp
    r = hyperexpand(hyper(ap, bq, z))
    if r.has(hyper):
        return False
    if not numerical:
        return True
    repl = {}
    for n, a in enumerate(r.free_symbols - {z}):
        repl[a] = randcplx(n)/div
    [a, b, c, d] = [2, -1, 3, 1]
    if lowerplane:
        [a, b, c, d] = [2, -2, 3, -1]
    return tn(
        hyper(ap, bq, z).subs(repl),
        r.replace(exp_polar, exp).subs(repl),
        z, a=a, b=b, c=c, d=d)


def test_roach():
    # Kelly B. Roach.  Meijer G Function Representations.
    # Section "Gallery"
    assert can_do([S(1)/2], [S(9)/2])
    assert can_do([], [1, S(5)/2, 4])
    assert can_do([-S.Half, 1, 2], [3, 4])
    assert can_do([S(1)/3], [-S(2)/3, -S(1)/2, S(1)/2, 1])
    assert can_do([-S(3)/2, -S(1)/2], [-S(5)/2, 1])
    assert can_do([-S(3)/2, ], [-S(1)/2, S(1)/2])  # shine-integral
    assert can_do([-S(3)/2, -S(1)/2], [2])  # elliptic integrals


@pytest.mark.xfail
def test_roach_fail():
    assert can_do([-S(1)/2, 1], [S(1)/4, S(1)/2, S(3)/4])  # PFDD
    assert can_do([S(3)/2], [S(5)/2, 5])  # struve function
    assert can_do([-S(1)/2, S(1)/2, 1], [S(3)/2, S(5)/2])  # polylog, pfdd
    assert can_do([1, 2, 3], [S(1)/2, 4])  # XXX ?
    assert can_do([S(1)/2], [-S(1)/3, -S(1)/2, -S(2)/3])  # PFDD ?

# For the long table tests, see end of file


def test_polynomial():
    from sympy import oo
    assert hyperexpand(hyper([], [-1], z)) == oo
    assert hyperexpand(hyper([-2], [-1], z)) == oo
    assert hyperexpand(hyper([0, 0], [-1], z)) == 1
    assert can_do([-5, -2, randcplx(), randcplx()], [-10, randcplx()])


def test_hyperexpand_bases():
    assert hyperexpand(hyper([2], [a], z)) == \
        a + z**(-a + 1)*(-a**2 + 3*a + z*(a - 1) - 2)*exp(z)* \
        lowergamma(a - 1, z) - 1
    # TODO [a+1, a-S.Half], [2*a]
    assert hyperexpand(hyper([1, 2], [3], z)) == -2/z - 2*log(-z + 1)/z**2
    assert hyperexpand(hyper([S.Half, 2], [S(3)/2], z)) == \
        -1/(2*z - 2) + atanh(sqrt(z))/sqrt(z)/2
    assert hyperexpand(hyper([S(1)/2, S(1)/2], [S(5)/2], z)) == \
        (-3*z + 3)/4/(z*sqrt(-z + 1)) \
        + (6*z - 3)*asin(sqrt(z))/(4*z**(S(3)/2))
    assert hyperexpand(hyper([1, 2], [S(3)/2], z)) == -1/(2*z - 2) \
        - asin(sqrt(z))/(sqrt(z)*(2*z - 2)*sqrt(-z + 1))
    assert hyperexpand(hyper([-S.Half - 1, 1, 2], [S.Half, 3], z)) == \
        sqrt(z)*(6*z/7 - S(6)/5)*atanh(sqrt(z)) \
        + (-30*z**2 + 32*z - 6)/35/z - 6*log(-z + 1)/(35*z**2)
    assert hyperexpand(hyper([1 + S.Half, 1, 1], [2, 2], z)) == \
        -4*log(sqrt(-z + 1)/2 + S(1)/2)/z
    # TODO hyperexpand(hyper([a], [2*a + 1], z))
    # TODO [S.Half, a], [S(3)/2, a+1]
    assert hyperexpand(hyper([2], [b, 1], z)) == \
        z**(-b/2 + S(1)/2)*besseli(b - 1, 2*sqrt(z))*gamma(b) \
        + z**(-b/2 + 1)*besseli(b, 2*sqrt(z))*gamma(b)
    # TODO [a], [a - S.Half, 2*a]


def test_hyperexpand_parametric():
    assert hyperexpand(hyper([a, S(1)/2 + a], [S(1)/2], z)) \
        == (1 + sqrt(z))**(-2*a)/2 + (1 - sqrt(z))**(-2*a)/2
    assert hyperexpand(hyper([a, -S(1)/2 + a], [2*a], z)) \
        == 2**(2*a - 1)*((-z + 1)**(S(1)/2) + 1)**(-2*a + 1)


def test_shifted_sum():
    from sympy import simplify
    assert simplify(hyperexpand(z**4*hyper([2], [3, S('3/2')], -z**2))) \
        == z*sin(2*z) + (-z**2 + S.Half)*cos(2*z) - S.Half


def _randrat():
    """ Steer clear of integers. """
    return S(randrange(25) + 10)/50


def randcplx(offset=-1):
    """ Polys is not good with real coefficients. """
    return _randrat() + I*_randrat() + I*(1 + offset)


@pytest.mark.slow
def test_formulae():
    from sympy.simplify.hyperexpand import FormulaCollection
    formulae = FormulaCollection().formulae
    for formula in formulae:
        h = formula.func(formula.z)
        rep = {}
        for n, sym in enumerate(formula.symbols):
            rep[sym] = randcplx(n)

        # NOTE hyperexpand returns truly branched functions. We know we are
        #      on the main sheet, but numerical evaluation can still go wrong
        #      (e.g. if exp_polar cannot be evalf'd).
        #      Just replace all exp_polar by exp, this usually works.

        # first test if the closed-form is actually correct
        h = h.subs(rep)
        closed_form = formula.closed_form.subs(rep).rewrite('nonrepsmall')
        z = formula.z
        assert tn(h, closed_form.replace(exp_polar, exp), z)

        # now test the computed matrix
        cl = (formula.C * formula.B)[0].subs(rep).rewrite('nonrepsmall')
        assert tn(closed_form.replace(
            exp_polar, exp), cl.replace(exp_polar, exp), z)
        deriv1 = z*formula.B.applyfunc(lambda t: t.rewrite(
            'nonrepsmall')).diff(z)
        deriv2 = formula.M * formula.B
        for d1, d2 in zip(deriv1, deriv2):
            assert tn(d1.subs(rep).replace(exp_polar, exp),
                      d2.subs(rep).rewrite('nonrepsmall').replace(exp_polar, exp), z)


def test_meijerg_formulae():
    from sympy.simplify.hyperexpand import MeijerFormulaCollection
    formulae = MeijerFormulaCollection().formulae
    for sig in formulae:
        for formula in formulae[sig]:
            g = meijerg(formula.func.an, formula.func.ap,
                        formula.func.bm, formula.func.bq,
                        formula.z)
            rep = {}
            for sym in formula.symbols:
                rep[sym] = randcplx()

            # first test if the closed-form is actually correct
            g = g.subs(rep)
            closed_form = formula.closed_form.subs(rep)
            z = formula.z
            assert tn(g, closed_form, z)

            # now test the computed matrix
            cl = (formula.C * formula.B)[0].subs(rep)
            assert tn(closed_form, cl, z)
            deriv1 = z*formula.B.diff(z)
            deriv2 = formula.M * formula.B
            for d1, d2 in zip(deriv1, deriv2):
                assert tn(d1.subs(rep), d2.subs(rep), z)


def op(f):
    return z*f.diff(z)


def test_plan():
    assert devise_plan(Hyper_Function([0], ()),
            Hyper_Function([0], ()), z) == []
    with pytest.raises(ValueError):
        devise_plan(Hyper_Function([1], ()), Hyper_Function((), ()), z)
    with pytest.raises(ValueError):
        devise_plan(Hyper_Function([2], [1]), Hyper_Function([2], [2]), z)
    with pytest.raises(ValueError):
        devise_plan(Hyper_Function([2], []), Hyper_Function([S("1/2")], []), z)

    # We cannot use pi/(10000 + n) because polys is insanely slow.
    a1, a2, b1 = (randcplx(n) for n in range(3))
    b1 += 2*I
    h = hyper([a1, a2], [b1], z)

    h2 = hyper((a1 + 1, a2), [b1], z)
    assert tn(apply_operators(h,
        devise_plan(Hyper_Function((a1 + 1, a2), [b1]),
            Hyper_Function((a1, a2), [b1]), z), op),
        h2, z)

    h2 = hyper((a1 + 1, a2 - 1), [b1], z)
    assert tn(apply_operators(h,
        devise_plan(Hyper_Function((a1 + 1, a2 - 1), [b1]),
            Hyper_Function((a1, a2), [b1]), z), op),
        h2, z)


def test_plan_derivatives():
    a1, a2, a3 = 1, 2, S('1/2')
    b1, b2 = 3, S('5/2')
    h = Hyper_Function((a1, a2, a3), (b1, b2))
    h2 = Hyper_Function((a1 + 1, a2 + 1, a3 + 2), (b1 + 1, b2 + 1))
    ops = devise_plan(h2, h, z)
    f = Formula(h, z, h(z), [])
    deriv = make_derivative_operator(f.M, z)
    assert tn((apply_operators(f.C, ops, deriv)*f.B)[0], h2(z), z)

    h2 = Hyper_Function((a1, a2 - 1, a3 - 2), (b1 - 1, b2 - 1))
    ops = devise_plan(h2, h, z)
    assert tn((apply_operators(f.C, ops, deriv)*f.B)[0], h2(z), z)


def test_reduction_operators():
    a1, a2, b1 = (randcplx(n) for n in range(3))
    h = hyper([a1], [b1], z)

    assert ReduceOrder(2, 0) is None
    assert ReduceOrder(2, -1) is None
    assert ReduceOrder(1, S('1/2')) is None

    h2 = hyper((a1, a2), (b1, a2), z)
    assert tn(ReduceOrder(a2, a2).apply(h, op), h2, z)

    h2 = hyper((a1, a2 + 1), (b1, a2), z)
    assert tn(ReduceOrder(a2 + 1, a2).apply(h, op), h2, z)

    h2 = hyper((a2 + 4, a1), (b1, a2), z)
    assert tn(ReduceOrder(a2 + 4, a2).apply(h, op), h2, z)

    # test several step order reduction
    ap = (a2 + 4, a1, b1 + 1)
    bq = (a2, b1, b1)
    func, ops = reduce_order(Hyper_Function(ap, bq))
    assert func.ap == (a1,)
    assert func.bq == (b1,)
    assert tn(apply_operators(h, ops, op), hyper(ap, bq, z), z)


def test_shift_operators():
    a1, a2, b1, b2, b3 = (randcplx(n) for n in range(5))
    h = hyper((a1, a2), (b1, b2, b3), z)

    pytest.raises(ValueError, lambda: ShiftA(0))
    pytest.raises(ValueError, lambda: ShiftB(1))

    assert tn(ShiftA(a1).apply(h, op), hyper((a1 + 1, a2), (b1, b2, b3), z), z)
    assert tn(ShiftA(a2).apply(h, op), hyper((a1, a2 + 1), (b1, b2, b3), z), z)
    assert tn(ShiftB(b1).apply(h, op), hyper((a1, a2), (b1 - 1, b2, b3), z), z)
    assert tn(ShiftB(b2).apply(h, op), hyper((a1, a2), (b1, b2 - 1, b3), z), z)
    assert tn(ShiftB(b3).apply(h, op), hyper((a1, a2), (b1, b2, b3 - 1), z), z)


def test_ushift_operators():
    a1, a2, b1, b2, b3 = (randcplx(n) for n in range(5))
    h = hyper((a1, a2), (b1, b2, b3), z)

    pytest.raises(ValueError, lambda: UnShiftA((1,), (), 0, z))
    pytest.raises(ValueError, lambda: UnShiftB((), (-1,), 0, z))
    pytest.raises(ValueError, lambda: UnShiftA((1,), (0, -1, 1), 0, z))
    pytest.raises(ValueError, lambda: UnShiftB((0, 1), (1,), 0, z))

    s = UnShiftA((a1, a2), (b1, b2, b3), 0, z)
    assert tn(s.apply(h, op), hyper((a1 - 1, a2), (b1, b2, b3), z), z)
    s = UnShiftA((a1, a2), (b1, b2, b3), 1, z)
    assert tn(s.apply(h, op), hyper((a1, a2 - 1), (b1, b2, b3), z), z)

    s = UnShiftB((a1, a2), (b1, b2, b3), 0, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1 + 1, b2, b3), z), z)
    s = UnShiftB((a1, a2), (b1, b2, b3), 1, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1, b2 + 1, b3), z), z)
    s = UnShiftB((a1, a2), (b1, b2, b3), 2, z)
    assert tn(s.apply(h, op), hyper((a1, a2), (b1, b2, b3 + 1), z), z)


def can_do_meijer(a1, a2, b1, b2, numeric=True):
    """
    This helper function tries to hyperexpand() the meijer g-function
    corresponding to the parameters a1, a2, b1, b2.
    It returns False if this expansion still contains g-functions.
    If numeric is True, it also tests the so-obtained formula numerically
    (at random values) and returns False if the test fails.
    Else it returns True.
    """
    from sympy import unpolarify, expand
    r = hyperexpand(meijerg(a1, a2, b1, b2, z))
    if r.has(meijerg):
        return False
    # NOTE hyperexpand() returns a truly branched function, whereas numerical
    #      evaluation only works on the main branch. Since we are evaluating on
    #      the main branch, this should not be a problem, but expressions like
    #      exp_polar(I*pi/2*x)**a are evaluated incorrectly. We thus have to get
    #      rid of them. The expand heuristically does this...
    r = unpolarify(expand(r, force=True, power_base=True, power_exp=False,
                          mul=False, log=False, multinomial=False, basic=False))

    if not numeric:
        return True

    repl = {}
    for n, a in enumerate(meijerg(a1, a2, b1, b2, z).free_symbols - {z}):
        repl[a] = randcplx(n)
    return tn(meijerg(a1, a2, b1, b2, z).subs(repl), r.subs(repl), z)


@pytest.mark.slow
def test_meijerg_expand():
    from sympy import combsimp, simplify
    # from mpmath docs
    assert hyperexpand(meijerg([[], []], [[0], []], -z)) == exp(z)

    assert hyperexpand(meijerg([[1, 1], []], [[1], [0]], z)) == \
        log(z + 1)
    assert hyperexpand(meijerg([[1, 1], []], [[1], [1]], z)) == \
        z/(z + 1)
    assert hyperexpand(meijerg([[], []], [[S(1)/2], [0]], (z/2)**2)) \
        == sin(z)/sqrt(pi)
    assert hyperexpand(meijerg([[], []], [[0], [S(1)/2]], (z/2)**2)) \
        == cos(z)/sqrt(pi)
    assert can_do_meijer([], [a], [a - 1, a - S.Half], [])
    assert can_do_meijer([], [], [a/2], [-a/2], False)  # branches...
    assert can_do_meijer([a], [b], [a], [b, a - 1])

    # wikipedia
    assert hyperexpand(meijerg([1], [], [], [0], z)) == \
        Piecewise((0, abs(z) < 1), (1, abs(1/z) < 1),
                 (meijerg([1], [], [], [0], z), True))
    assert hyperexpand(meijerg([], [1], [0], [], z)) == \
        Piecewise((1, abs(z) < 1), (0, abs(1/z) < 1),
                 (meijerg([], [1], [0], [], z), True))

    # The Special Functions and their Approximations
    assert can_do_meijer([], [], [a + b/2], [a, a - b/2, a + S.Half])
    assert can_do_meijer(
        [], [], [a], [b], False)  # branches only agree for small z
    assert can_do_meijer([], [S.Half], [a], [-a])
    assert can_do_meijer([], [], [a, b], [])
    assert can_do_meijer([], [], [a, b], [])
    assert can_do_meijer([], [], [a, a + S.Half], [b, b + S.Half])
    assert can_do_meijer([], [], [a, -a], [0, S.Half], False)  # dito
    assert can_do_meijer([], [], [a, a + S.Half, b, b + S.Half], [])
    assert can_do_meijer([S.Half], [], [0], [a, -a])
    assert can_do_meijer([S.Half], [], [a], [0, -a], False)  # dito
    assert can_do_meijer([], [a - S.Half], [a, b], [a - S.Half], False)
    assert can_do_meijer([], [a + S.Half], [a + b, a - b, a], [], False)
    assert can_do_meijer([a + S.Half], [], [b, 2*a - b, a], [], False)

    # This for example is actually zero.
    assert can_do_meijer([], [], [], [a, b])

    # Testing a bug:
    assert hyperexpand(meijerg([0, 2], [], [], [-1, 1], z)) == \
        Piecewise((0, abs(z) < 1),
                  (z/2 - 1/(2*z), abs(1/z) < 1),
                  (meijerg([0, 2], [], [], [-1, 1], z), True))

    # Test that the simplest possible answer is returned:
    assert combsimp(simplify(hyperexpand(
        meijerg([1], [1 - a], [-a/2, -a/2 + S(1)/2], [], 1/z)))) == \
        -2*sqrt(pi)*(sqrt(z + 1) + 1)**a/a

    # Test that hyper is returned
    assert hyperexpand(meijerg([1], [], [a], [0, 0], z)) == hyper(
        (a,), (a + 1, a + 1), z*exp_polar(I*pi))*z**a*gamma(a)/gamma(a + 1)**2


def test_meijerg_lookup():
    from sympy import uppergamma, Si, Ci
    assert hyperexpand(meijerg([a], [], [b, a], [], z)) == \
        z**b*exp(z)*gamma(-a + b + 1)*uppergamma(a - b, z)
    assert hyperexpand(meijerg([0], [], [0, 0], [], z)) == \
        exp(z)*uppergamma(0, z)
    assert can_do_meijer([a], [], [b, a + 1], [])
    assert can_do_meijer([a], [], [b + 2, a], [])
    assert can_do_meijer([a], [], [b - 2, a], [])

    assert hyperexpand(meijerg([a], [], [a, a, a - S(1)/2], [], z)) == \
        -sqrt(pi)*z**(a - S(1)/2)*(2*cos(2*sqrt(z))*(Si(2*sqrt(z)) - pi/2)
                                   - 2*sin(2*sqrt(z))*Ci(2*sqrt(z))) == \
        hyperexpand(meijerg([a], [], [a, a - S(1)/2, a], [], z)) == \
        hyperexpand(meijerg([a], [], [a - S(1)/2, a, a], [], z))
    assert can_do_meijer([a - 1], [], [a + 2, a - S(3)/2, a + 1], [])


@pytest.mark.xfail
def test_meijerg_expand_fail():
    # These basically test hyper([], [1/2 - a, 1/2 + 1, 1/2], z),
    # which is *very* messy. But since the meijer g actually yields a
    # sum of bessel functions, things can sometimes be simplified a lot and
    # are then put into tables...
    assert can_do_meijer([], [], [a + S.Half], [a, a - b/2, a + b/2])
    assert can_do_meijer([], [], [0, S.Half], [a, -a])
    assert can_do_meijer([], [], [3*a - S.Half, a, -a - S.Half], [a - S.Half])
    assert can_do_meijer([], [], [0, a - S.Half, -a - S.Half], [S.Half])
    assert can_do_meijer([], [], [a, b + S(1)/2, b], [2*b - a])
    assert can_do_meijer([], [], [a, b + S(1)/2, b, 2*b - a])
    assert can_do_meijer([S.Half], [], [-a, a], [0])


@pytest.mark.slow
def test_meijerg():
    # carefully set up the parameters.
    # NOTE: this used to fail sometimes. I believe it is fixed, but if you
    #       hit an inexplicable test failure here, please let me know the seed.
    a1, a2 = (randcplx(n) - 5*I - n*I for n in range(2))
    b1, b2 = (randcplx(n) + 5*I + n*I for n in range(2))
    b3, b4, b5, a3, a4, a5 = (randcplx() for n in range(6))
    g = meijerg([a1], [a3, a4], [b1], [b3, b4], z)

    assert ReduceOrder.meijer_minus(3, 4) is None
    assert ReduceOrder.meijer_plus(4, 3) is None

    g2 = meijerg([a1, a2], [a3, a4], [b1], [b3, b4, a2], z)
    assert tn(ReduceOrder.meijer_plus(a2, a2).apply(g, op), g2, z)

    g2 = meijerg([a1, a2], [a3, a4], [b1], [b3, b4, a2 + 1], z)
    assert tn(ReduceOrder.meijer_plus(a2, a2 + 1).apply(g, op), g2, z)

    g2 = meijerg([a1, a2 - 1], [a3, a4], [b1], [b3, b4, a2 + 2], z)
    assert tn(ReduceOrder.meijer_plus(a2 - 1, a2 + 2).apply(g, op), g2, z)

    g2 = meijerg([a1], [a3, a4, b2 - 1], [b1, b2 + 2], [b3, b4], z)
    assert tn(ReduceOrder.meijer_minus(
        b2 + 2, b2 - 1).apply(g, op), g2, z, tol=1e-6)

    # test several-step reduction
    an = [a1, a2]
    bq = [b3, b4, a2 + 1]
    ap = [a3, a4, b2 - 1]
    bm = [b1, b2 + 1]
    niq, ops = reduce_order_meijer(G_Function(an, ap, bm, bq))
    assert niq.an == (a1,)
    assert set(niq.ap) == {a3, a4}
    assert niq.bm == (b1,)
    assert set(niq.bq) == {b3, b4}
    assert tn(apply_operators(g, ops, op), meijerg(an, ap, bm, bq, z), z)


def test_meijerg_shift_operators():
    # carefully set up the parameters. XXX this still fails sometimes
    a1, a2, a3, a4, a5, b1, b2, b3, b4, b5 = (randcplx(n) for n in range(10))
    g = meijerg([a1], [a3, a4], [b1], [b3, b4], z)

    assert tn(MeijerShiftA(b1).apply(g, op),
              meijerg([a1], [a3, a4], [b1 + 1], [b3, b4], z), z)
    assert tn(MeijerShiftB(a1).apply(g, op),
              meijerg([a1 - 1], [a3, a4], [b1], [b3, b4], z), z)
    assert tn(MeijerShiftC(b3).apply(g, op),
              meijerg([a1], [a3, a4], [b1], [b3 + 1, b4], z), z)
    assert tn(MeijerShiftD(a3).apply(g, op),
              meijerg([a1], [a3 - 1, a4], [b1], [b3, b4], z), z)

    s = MeijerUnShiftA([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(
        s.apply(g, op), meijerg([a1], [a3, a4], [b1 - 1], [b3, b4], z), z)

    s = MeijerUnShiftC([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(
        s.apply(g, op), meijerg([a1], [a3, a4], [b1], [b3 - 1, b4], z), z)

    s = MeijerUnShiftB([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(
        s.apply(g, op), meijerg([a1 + 1], [a3, a4], [b1], [b3, b4], z), z)

    s = MeijerUnShiftD([a1], [a3, a4], [b1], [b3, b4], 0, z)
    assert tn(
        s.apply(g, op), meijerg([a1], [a3 + 1, a4], [b1], [b3, b4], z), z)


@pytest.mark.slow
def test_meijerg_confluence():
    def t(m, a, b):
        from sympy import sympify, Piecewise
        a, b = sympify([a, b])
        m_ = m
        m = hyperexpand(m)
        if not m == Piecewise((a, abs(z) < 1), (b, abs(1/z) < 1), (m_, True)):
            return False
        if not (m.args[0].args[0] == a and m.args[1].args[0] == b):
            return False
        z0 = randcplx()/10
        if abs(m.subs(z, z0).n() - a.subs(z, z0).n()).n() > 1e-10:
            return False
        if abs(m.subs(z, 1/z0).n() - b.subs(z, 1/z0).n()).n() > 1e-10:
            return False
        return True

    assert t(meijerg([], [1, 1], [0, 0], [], z), -log(z), 0)
    assert t(meijerg(
        [], [3, 1], [0, 0], [], z), -z**2/4 + z - log(z)/2 - S(3)/4, 0)
    assert t(meijerg([], [3, 1], [-1, 0], [], z),
             z**2/12 - z/2 + log(z)/2 + S(1)/4 + 1/(6*z), 0)
    assert t(meijerg([], [1, 1, 1, 1], [0, 0, 0, 0], [], z), -log(z)**3/6, 0)
    assert t(meijerg([1, 1], [], [], [0, 0], z), 0, -log(1/z))
    assert t(meijerg([1, 1], [2, 2], [1, 1], [0, 0], z),
             -z*log(z) + 2*z, -log(1/z) + 2)
    assert t(meijerg([S(1)/2], [1, 1], [0, 0], [S(3)/2], z), log(z)/2 - 1, 0)

    def u(an, ap, bm, bq):
        m = meijerg(an, ap, bm, bq, z)
        m2 = hyperexpand(m, allow_hyper=True)
        if m2.has(meijerg) and not (m2.is_Piecewise and len(m2.args) == 3):
            return False
        return tn(m, m2, z)
    assert u([], [1], [0, 0], [])
    assert u([1, 1], [], [], [0])
    assert u([1, 1], [2, 2, 5], [1, 1, 6], [0, 0])
    assert u([1, 1], [2, 2, 5], [1, 1, 6], [0])


def test_lerchphi():
    from sympy import combsimp, exp_polar, polylog, log, lerchphi
    assert hyperexpand(hyper([1, a], [a + 1], z)/a) == lerchphi(z, 1, a)
    assert hyperexpand(
        hyper([1, a, a], [a + 1, a + 1], z)/a**2) == lerchphi(z, 2, a)
    assert hyperexpand(hyper([1, a, a, a], [a + 1, a + 1, a + 1], z)/a**3) == \
        lerchphi(z, 3, a)
    assert hyperexpand(hyper([1] + [a]*10, [a + 1]*10, z)/a**10) == \
        lerchphi(z, 10, a)
    assert combsimp(hyperexpand(meijerg([0, 1 - a], [], [0],
        [-a], exp_polar(-I*pi)*z))) == lerchphi(z, 1, a)
    assert combsimp(hyperexpand(meijerg([0, 1 - a, 1 - a], [], [0],
        [-a, -a], exp_polar(-I*pi)*z))) == lerchphi(z, 2, a)
    assert combsimp(hyperexpand(meijerg([0, 1 - a, 1 - a, 1 - a], [], [0],
        [-a, -a, -a], exp_polar(-I*pi)*z))) == lerchphi(z, 3, a)

    assert hyperexpand(z*hyper([1, 1], [2], z)) == -log(1 + -z)
    assert hyperexpand(z*hyper([1, 1, 1], [2, 2], z)) == polylog(2, z)
    assert hyperexpand(z*hyper([1, 1, 1, 1], [2, 2, 2], z)) == polylog(3, z)

    assert hyperexpand(hyper([1, a, 1 + S(1)/2], [a + 1, S(1)/2], z)) == \
        -2*a/(z - 1) + (-2*a**2 + a)*lerchphi(z, 1, a)

    # Now numerical tests. These make sure reductions etc are carried out
    # correctly

    # a rational function (polylog at negative integer order)
    assert can_do([2, 2, 2], [1, 1])

    # NOTE these contain log(1-x) etc ... better make sure we have |z| < 1
    # reduction of order for polylog
    assert can_do([1, 1, 1, b + 5], [2, 2, b], div=10)

    # reduction of order for lerchphi
    # XXX lerchphi in mpmath is flaky
    assert can_do(
        [1, a, a, a, b + 5], [a + 1, a + 1, a + 1, b], numerical=False)

    # test a bug
    from sympy import Abs
    assert hyperexpand(hyper([S(1)/2, S(1)/2, S(1)/2, 1],
                             [S(3)/2, S(3)/2, S(3)/2], S(1)/4)) == \
        Abs(-polylog(3, exp_polar(I*pi)/2) + polylog(3, S(1)/2))


def test_partial_simp():
    # First test that hypergeometric function formulae work.
    a, b, c, d, e = (randcplx() for _ in range(5))
    for func in [Hyper_Function([a, b, c], [d, e]),
            Hyper_Function([], [a, b, c, d, e])]:
        f = build_hypergeometric_formula(func)
        z = f.z
        assert f.closed_form == func(z)
        deriv1 = f.B.diff(z)*z
        deriv2 = f.M*f.B
        for func1, func2 in zip(deriv1, deriv2):
            assert tn(func1, func2, z)

    # Now test that formulae are partially simplified.
    from sympy.abc import a, b, z
    assert hyperexpand(hyper([3, a], [1, b], z)) == \
        (-a*b/2 + a*z/2 + 2*a)*hyper([a + 1], [b], z) \
        + (a*b/2 - 2*a + 1)*hyper([a], [b], z)
    assert tn(
        hyperexpand(hyper([3, d], [1, e], z)), hyper([3, d], [1, e], z), z)
    assert hyperexpand(hyper([3], [1, a, b], z)) == \
        hyper((), (a, b), z) \
        + z*hyper((), (a + 1, b), z)/(2*a) \
        - z*(b - 4)*hyper((), (a + 1, b + 1), z)/(2*a*b)
    assert tn(
        hyperexpand(hyper([3], [1, d, e], z)), hyper([3], [1, d, e], z), z)


def test_hyperexpand_special():
    assert hyperexpand(hyper([a, b], [c], 1)) == \
        gamma(c)*gamma(c - a - b)/gamma(c - a)/gamma(c - b)
    assert hyperexpand(hyper([a, b], [1 + a - b], -1)) == \
        gamma(1 + a/2)*gamma(1 + a - b)/gamma(1 + a)/gamma(1 + a/2 - b)
    assert hyperexpand(hyper([a, b], [1 + b - a], -1)) == \
        gamma(1 + b/2)*gamma(1 + b - a)/gamma(1 + b)/gamma(1 + b/2 - a)
    assert hyperexpand(meijerg([1 - z - a/2], [1 - z + a/2], [b/2], [-b/2], 1)) == \
        gamma(1 - 2*z)*gamma(z + a/2 + b/2)/gamma(1 - z + a/2 - b/2) \
        /gamma(1 - z - a/2 + b/2)/gamma(1 - z + a/2 + b/2)
    assert hyperexpand(hyper([a], [b], 0)) == 1
    assert hyper([a], [b], 0) != 0


def test_Mod1_behavior():
    from sympy import Symbol, simplify, lowergamma
    n = Symbol('n', integer=True)
    # Note: this should not hang.
    assert simplify(hyperexpand(meijerg([1], [], [n + 1], [0], z))) == \
        lowergamma(n + 1, z)


@pytest.mark.slow
def test_prudnikov_misc():
    assert can_do([1, (3 + I)/2, (3 - I)/2], [S(3)/2, 2])
    assert can_do([S.Half, a - 1], [S(3)/2, a + 1], lowerplane=True)
    assert can_do([], [b + 1])
    assert can_do([a], [a - 1, b + 1])

    assert can_do([a], [a - S.Half, 2*a])
    assert can_do([a], [a - S.Half, 2*a + 1])
    assert can_do([a], [a - S.Half, 2*a - 1])
    assert can_do([a], [a + S.Half, 2*a])
    assert can_do([a], [a + S.Half, 2*a + 1])
    assert can_do([a], [a + S.Half, 2*a - 1])
    assert can_do([S.Half], [b, 2 - b])
    assert can_do([S.Half], [b, 3 - b])
    assert can_do([1], [2, b])

    assert can_do([a, a + S.Half], [2*a, b, 2*a - b + 1])
    assert can_do([a, a + S.Half], [S.Half, 2*a, 2*a + S.Half])
    assert can_do([a], [a + 1], lowerplane=True)  # lowergamma


@pytest.mark.slow
def test_prudnikov_1():
    # A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
    # Integrals and Series: More Special Functions, Vol. 3,.
    # Gordon and Breach Science Publisher

    # 7.3.1
    assert can_do([a, -a], [S.Half])
    assert can_do([a, 1 - a], [S.Half])
    assert can_do([a, 1 - a], [S(3)/2])
    assert can_do([a, 2 - a], [S.Half])
    assert can_do([a, 2 - a], [S(3)/2])
    assert can_do([a, 2 - a], [S(3)/2])
    assert can_do([a, a + S(1)/2], [2*a - 1])
    assert can_do([a, a + S(1)/2], [2*a])
    assert can_do([a, a + S(1)/2], [2*a + 1])
    assert can_do([a, a + S(1)/2], [S(1)/2])
    assert can_do([a, a + S(1)/2], [S(3)/2])
    assert can_do([a, a/2 + 1], [a/2])
    assert can_do([1, b], [2])
    assert can_do([1, b], [b + 1], numerical=False)  # Lerch Phi
             # NOTE: branches are complicated for |z| > 1

    assert can_do([a], [2*a])
    assert can_do([a], [2*a + 1])
    assert can_do([a], [2*a - 1])


@pytest.mark.slow
def test_prudnikov_2():
    h = S.Half
    assert can_do([-h, -h], [h])
    assert can_do([-h, h], [3*h])
    assert can_do([-h, h], [5*h])
    assert can_do([-h, h], [7*h])
    assert can_do([-h, 1], [h])

    for p in [-h, h]:
        for n in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]:
            for m in [-h, h, 3*h, 5*h, 7*h]:
                assert can_do([p, n], [m])
        for n in [1, 2, 3, 4]:
            for m in [1, 2, 3, 4]:
                assert can_do([p, n], [m])


@pytest.mark.slow
def test_prudnikov_3():
    h = S.Half
    assert can_do([S(1)/4, S(3)/4], [h])
    assert can_do([S(1)/4, S(3)/4], [3*h])
    assert can_do([S(1)/3, S(2)/3], [3*h])
    assert can_do([S(3)/4, S(5)/4], [h])
    assert can_do([S(3)/4, S(5)/4], [3*h])

    for p in [1, 2, 3, 4]:
        for n in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4, 9*h]:
            for m in [1, 3*h, 2, 5*h, 3, 7*h, 4]:
                assert can_do([p, m], [n])


@pytest.mark.slow
def test_prudnikov_4():
    h = S.Half
    for p in [3*h, 5*h, 7*h]:
        for n in [-h, h, 3*h, 5*h, 7*h]:
            for m in [3*h, 2, 5*h, 3, 7*h, 4]:
                assert can_do([p, m], [n])
        for n in [1, 2, 3, 4]:
            for m in [2, 3, 4]:
                assert can_do([p, m], [n])


@pytest.mark.slow
def test_prudnikov_5():
    h = S.Half

    for p in [1, 2, 3]:
        for q in range(p, 4):
            for r in [1, 2, 3]:
                for s in range(r, 4):
                    assert can_do([-h, p, q], [r, s])

    for p in [h, 1, 3*h, 2, 5*h, 3]:
        for q in [h, 3*h, 5*h]:
            for r in [h, 3*h, 5*h]:
                for s in [h, 3*h, 5*h]:
                    if s <= q and s <= r:
                        assert can_do([-h, p, q], [r, s])

    for p in [h, 1, 3*h, 2, 5*h, 3]:
        for q in [1, 2, 3]:
            for r in [h, 3*h, 5*h]:
                for s in [1, 2, 3]:
                    assert can_do([-h, p, q], [r, s])


@pytest.mark.slow
def test_prudnikov_6():
    h = S.Half

    for m in [3*h, 5*h]:
        for n in [1, 2, 3]:
            for q in [h, 1, 2]:
                for p in [1, 2, 3]:
                    assert can_do([h, q, p], [m, n])
            for q in [1, 2, 3]:
                for p in [3*h, 5*h]:
                    assert can_do([h, q, p], [m, n])

    for q in [1, 2]:
        for p in [1, 2, 3]:
            for m in [1, 2, 3]:
                for n in [1, 2, 3]:
                    assert can_do([h, q, p], [m, n])

    assert can_do([h, h, 5*h], [3*h, 3*h])
    assert can_do([h, 1, 5*h], [3*h, 3*h])
    assert can_do([h, 2, 2], [1, 3])

    # pages 435 to 457 contain more PFDD and stuff like this


@pytest.mark.slow
def test_prudnikov_7():
    assert can_do([3], [6])

    h = S.Half
    for n in [h, 3*h, 5*h, 7*h]:
        assert can_do([-h], [n])
    for m in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]:  # HERE
        for n in [-h, h, 3*h, 5*h, 7*h, 1, 2, 3, 4]:
            assert can_do([m], [n])


@pytest.mark.slow
def test_prudnikov_8():
    h = S.Half

    # 7.12.2
    for a in [1, 2, 3]:
        for b in [1, 2, 3]:
            for c in range(1, a + 1):
                for d in [h, 1, 3*h, 2, 5*h, 3]:
                    assert can_do([a, b], [c, d])
        for b in [3*h, 5*h]:
            for c in [h, 1, 3*h, 2, 5*h, 3]:
                for d in [1, 2, 3]:
                    assert can_do([a, b], [c, d])

    for a in [-h, h, 3*h, 5*h]:
        for b in [1, 2, 3]:
            for c in [h, 1, 3*h, 2, 5*h, 3]:
                for d in [1, 2, 3]:
                    assert can_do([a, b], [c, d])
        for b in [h, 3*h, 5*h]:
            for c in [h, 3*h, 5*h, 3]:
                for d in [h, 1, 3*h, 2, 5*h, 3]:
                    if c <= b:
                        assert can_do([a, b], [c, d])


@pytest.mark.slow
def test_prudnikov_9():
    # 7.13.1 [we have a general formula ... so this is a bit pointless]
    for i in range(9):
        assert can_do([], [(S(i) + 1)/2])
    for i in range(5):
        assert can_do([], [-(2*S(i) + 1)/2])


@pytest.mark.slow
def test_prudnikov_10():
    # 7.14.2
    h = S.Half
    for p in [-h, h, 1, 3*h, 2, 5*h, 3, 7*h, 4]:
        for m in [1, 2, 3, 4]:
            for n in range(m, 5):
                assert can_do([p], [m, n])

    for p in [1, 2, 3, 4]:
        for n in [h, 3*h, 5*h, 7*h]:
            for m in [1, 2, 3, 4]:
                assert can_do([p], [n, m])

    for p in [3*h, 5*h, 7*h]:
        for m in [h, 1, 2, 5*h, 3, 7*h, 4]:
            assert can_do([p], [h, m])
            assert can_do([p], [3*h, m])

    for m in [h, 1, 2, 5*h, 3, 7*h, 4]:
        assert can_do([7*h], [5*h, m])

    assert can_do([-S(1)/2], [S(1)/2, S(1)/2])  # shine-integral shi


@pytest.mark.slow
def test_prudnikov_11():
    # 7.15
    assert can_do([a, a + S.Half], [2*a, b, 2*a - b])
    assert can_do([a, a + S.Half], [S(3)/2, 2*a, 2*a - S(1)/2])

    assert can_do([S(1)/4, S(3)/4], [S(1)/2, S(1)/2, 1])
    assert can_do([S(5)/4, S(3)/4], [S(3)/2, S(1)/2, 2])
    assert can_do([S(5)/4, S(3)/4], [S(3)/2, S(3)/2, 1])
    assert can_do([S(5)/4, S(7)/4], [S(3)/2, S(5)/2, 2])

    assert can_do([1, 1], [S(3)/2, 2, 2])  # cosh-integral chi


@pytest.mark.slow
def test_prudnikov_12():
    # 7.16
    assert can_do(
        [], [a, a + S.Half, 2*a], False)  # branches only agree for some z!
    assert can_do([], [a, a + S.Half, 2*a + 1], False)  # dito
    assert can_do([], [S.Half, a, a + S.Half])
    assert can_do([], [S(3)/2, a, a + S.Half])

    assert can_do([], [S(1)/4, S(1)/2, S(3)/4])
    assert can_do([], [S(1)/2, S(1)/2, 1])
    assert can_do([], [S(1)/2, S(3)/2, 1])
    assert can_do([], [S(3)/4, S(3)/2, S(5)/4])
    assert can_do([], [1, 1, S(3)/2])
    assert can_do([], [1, 2, S(3)/2])
    assert can_do([], [1, S(3)/2, S(3)/2])
    assert can_do([], [S(5)/4, S(3)/2, S(7)/4])
    assert can_do([], [2, S(3)/2, S(3)/2])


@pytest.mark.slow
def test_prudnikov_2F1():
    h = S.Half
    # Elliptic integrals
    for p in [-h, h]:
        for m in [h, 3*h, 5*h, 7*h]:
            for n in [1, 2, 3, 4]:
                assert can_do([p, m], [n])


@pytest.mark.xfail
def test_prudnikov_fail_2F1():
    assert can_do([a, b], [b + 1])  # incomplete beta function
    assert can_do([-1, b], [c])    # Poly. also -2, -3 etc

    # TODO polys

    # Legendre functions:
    assert can_do([a, b], [a + b + S.Half])
    assert can_do([a, b], [a + b - S.Half])
    assert can_do([a, b], [a + b + S(3)/2])
    assert can_do([a, b], [(a + b + 1)/2])
    assert can_do([a, b], [(a + b)/2 + 1])
    assert can_do([a, b], [a - b + 1])
    assert can_do([a, b], [a - b + 2])
    assert can_do([a, b], [2*b])
    assert can_do([a, b], [S.Half])
    assert can_do([a, b], [S(3)/2])
    assert can_do([a, 1 - a], [c])
    assert can_do([a, 2 - a], [c])
    assert can_do([a, 3 - a], [c])
    assert can_do([a, a + S(1)/2], [c])
    assert can_do([1, b], [c])
    assert can_do([1, b], [S(3)/2])

    assert can_do([S(1)/4, S(3)/4], [1])

    # PFDD
    o = S(1)
    assert can_do([o/8, 1], [o/8*9])
    assert can_do([o/6, 1], [o/6*7])
    assert can_do([o/6, 1], [o/6*13])
    assert can_do([o/5, 1], [o/5*6])
    assert can_do([o/5, 1], [o/5*11])
    assert can_do([o/4, 1], [o/4*5])
    assert can_do([o/4, 1], [o/4*9])
    assert can_do([o/3, 1], [o/3*4])
    assert can_do([o/3, 1], [o/3*7])
    assert can_do([o/8*3, 1], [o/8*11])
    assert can_do([o/5*2, 1], [o/5*7])
    assert can_do([o/5*2, 1], [o/5*12])
    assert can_do([o/5*3, 1], [o/5*8])
    assert can_do([o/5*3, 1], [o/5*13])
    assert can_do([o/8*5, 1], [o/8*13])
    assert can_do([o/4*3, 1], [o/4*7])
    assert can_do([o/4*3, 1], [o/4*11])
    assert can_do([o/3*2, 1], [o/3*5])
    assert can_do([o/3*2, 1], [o/3*8])
    assert can_do([o/5*4, 1], [o/5*9])
    assert can_do([o/5*4, 1], [o/5*14])
    assert can_do([o/6*5, 1], [o/6*11])
    assert can_do([o/6*5, 1], [o/6*17])
    assert can_do([o/8*7, 1], [o/8*15])


@pytest.mark.xfail
def test_prudnikov_fail_3F2():
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(1)/3, S(2)/3])
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(2)/3, S(4)/3])
    assert can_do([a, a + S(1)/3, a + S(2)/3], [S(4)/3, S(5)/3])

    # page 421
    assert can_do([a, a + S(1)/3, a + S(2)/3], [3*a/2, (3*a + 1)/2])

    # pages 422 ...
    assert can_do([-S.Half, S.Half, S.Half], [1, 1])  # elliptic integrals
    assert can_do([-S.Half, S.Half, 1], [S(3)/2, S(3)/2])
    # TODO LOTS more

    # PFDD
    assert can_do([S(1)/8, S(3)/8, 1], [S(9)/8, S(11)/8])
    assert can_do([S(1)/8, S(5)/8, 1], [S(9)/8, S(13)/8])
    assert can_do([S(1)/8, S(7)/8, 1], [S(9)/8, S(15)/8])
    assert can_do([S(1)/6, S(1)/3, 1], [S(7)/6, S(4)/3])
    assert can_do([S(1)/6, S(2)/3, 1], [S(7)/6, S(5)/3])
    assert can_do([S(1)/6, S(2)/3, 1], [S(5)/3, S(13)/6])
    assert can_do([S.Half, 1, 1], [S(1)/4, S(3)/4])
    # LOTS more


@pytest.mark.xfail
def test_prudnikov_fail_other():
    # 7.11.2

    # 7.12.1
    assert can_do([1, a], [b, 1 - 2*a + b])  # ???

    # 7.14.2
    assert can_do([-S(1)/2], [S(1)/2, 1])  # struve
    assert can_do([1], [S(1)/2, S(1)/2])  # struve
    assert can_do([S(1)/4], [S(1)/2, S(5)/4])  # PFDD
    assert can_do([S(3)/4], [S(3)/2, S(7)/4])  # PFDD
    assert can_do([1], [S(1)/4, S(3)/4])  # PFDD
    assert can_do([1], [S(3)/4, S(5)/4])  # PFDD
    assert can_do([1], [S(5)/4, S(7)/4])  # PFDD
    # TODO LOTS more

    # 7.15.2
    assert can_do([S(1)/2, 1], [S(3)/4, S(5)/4, S(3)/2])  # PFDD
    assert can_do([S(1)/2, 1], [S(7)/4, S(5)/4, S(3)/2])  # PFDD

    # 7.16.1
    assert can_do([], [S(1)/3, S(2/3)])  # PFDD
    assert can_do([], [S(2)/3, S(4/3)])  # PFDD
    assert can_do([], [S(5)/3, S(4/3)])  # PFDD

    # XXX this does not *evaluate* right??
    assert can_do([], [a, a + S.Half, 2*a - 1])


def test_bug():
    h = hyper([-1, 1], [z], -1)
    assert hyperexpand(h) == (z + 1)/z
