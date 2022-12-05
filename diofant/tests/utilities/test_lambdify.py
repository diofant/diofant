import itertools
import math

import mpmath
import pytest

import diofant
from diofant import (ITE, And, Float, Function, I, Integral, Lambda, Matrix,
                     Max, Min, Not, Or, Piecewise, Rational, cos, exp, false,
                     lambdify, oo, pi, sin, sqrt, symbols, true)
from diofant.abc import t, w, x, y, z
from diofant.external import import_module
from diofant.printing.lambdarepr import LambdaPrinter
from diofant.utilities.decorator import conserve_mpmath_dps
from diofant.utilities.lambdify import (MATH_TRANSLATIONS, MPMATH_TRANSLATIONS,
                                        NUMPY_TRANSLATIONS, _get_namespace,
                                        implemented_function, lambdastr)


__all__ = ()

MutableDenseMatrix = Matrix

numpy = import_module('numpy')
with_numpy = pytest.mark.skipif(numpy is None,
                                reason="Couldn't import numpy.")

# ================= Test different arguments =======================


def test_no_args():
    f = lambdify([], 1)
    pytest.raises(TypeError, lambda: f(-1))
    assert f() == 1


def test_single_arg():
    f = lambdify(x, 2*x)
    assert f(1) == 2


def test_list_args():
    f = lambdify([x, y], x + y)
    assert f(1, 2) == 3


def test_nested_args():
    # issue sympy/sympy#2790
    assert lambdify((x, (y, z)), x + y)(1, (2, 4)) == 3
    assert lambdify((x, (y, (w, z))), w + x + y + z)(1, (2, (3, 4))) == 10
    assert lambdify(x, x + 1, dummify=False)(1) == 2


def test_str_args():
    f = lambdify('x,y,z', 'z,y,x')
    assert f(3, 2, 1) == (1, 2, 3)
    assert f(1.0, 2.0, 3.0) == (3.0, 2.0, 1.0)
    # make sure correct number of args required
    pytest.raises(TypeError, lambda: f(0))


def test_own_namespace():
    def myfunc(x):
        return 1
    f = lambdify(x, sin(x), {'sin': myfunc})
    assert f(0.1) == 1
    assert f(100) == 1


def test_own_module():
    f = lambdify(x, sin(x), math)
    assert f(0) == 0.0


def test_bad_args():
    # no vargs given
    pytest.raises(TypeError, lambda: lambdify(1))
    # same with vector exprs
    pytest.raises(TypeError, lambda: lambdify([1, 2]))
    # reserved name
    pytest.raises(ValueError, lambda: lambdify((('__flatten_args__',),), 1))

    pytest.raises(NameError, lambda: lambdify(x, 1, 'spam'))


def test__get_namespace():
    pytest.raises(TypeError, lambda: _get_namespace(1))


def test_lambdastr():
    assert lambdastr(x, x**2) == 'lambda x: (x**2)'
    assert lambdastr(x, None, dummify=True).find('None') > 0


def test_atoms():
    # Non-Symbol atoms should not be pulled out from the expression namespace
    f = lambdify(x, pi + x, {'pi': 3.14})
    assert f(0) == 3.14
    f = lambdify(x, I + x, {'I': 1j})
    assert f(1) == 1 + 1j

# ================= Test different modules =========================

# high precision output of sin(0.2*pi) is used to detect if precision is lost unwanted


@conserve_mpmath_dps
def test_diofant_lambda():
    mpmath.mp.dps = 50
    sin02 = mpmath.mpf('0.19866933079506121545941262711838975037020672954020')
    f = lambdify(x, sin(x), 'diofant')
    assert f(x) == sin(x)
    prec = 1e-15
    assert -prec < f(Rational(1, 5)).evalf() - Float(str(sin02)) < prec


@conserve_mpmath_dps
def test_math_lambda():
    mpmath.mp.dps = 50
    sin02 = mpmath.mpf('0.19866933079506121545941262711838975037020672954020')
    f = lambdify(x, sin(x), 'math')
    prec = 1e-15
    assert -prec < f(0.2) - sin02 < prec

    # if this succeeds, it can't be a python math function
    pytest.raises(TypeError, lambda: f(x))


@conserve_mpmath_dps
def test_mpmath_lambda():
    mpmath.mp.dps = 50
    sin02 = mpmath.mpf('0.19866933079506121545941262711838975037020672954020')
    f = lambdify(x, sin(x), 'mpmath')
    prec = 1e-49  # mpmath precision is around 50 decimal places
    assert -prec < f(mpmath.mpf('0.2')) - sin02 < prec

    # if this succeeds, it can't be a mpmath function
    pytest.raises(TypeError, lambda: f(x))


@conserve_mpmath_dps
def test_number_precision():
    mpmath.mp.dps = 50
    sin02 = mpmath.mpf('0.19866933079506121545941262711838975037020672954020')
    f = lambdify(x, sin02, 'mpmath')
    prec = 1e-49  # mpmath precision is around 50 decimal places
    assert -prec < f(0) - sin02 < prec


@conserve_mpmath_dps
def test_mpmath_precision():
    mpmath.mp.dps = 100
    assert str(lambdify((), pi.evalf(100), 'mpmath')()) == str(pi.evalf(100))


# ================= Test Translations ==============================
# We can only check if all translated functions are valid. It has to be checked
# by hand if they are complete.


def test_math_transl():
    for sym, mat in MATH_TRANSLATIONS.items():
        assert sym in diofant.__dict__
        assert mat in math.__dict__


def test_mpmath_transl():
    for sym, mat in MPMATH_TRANSLATIONS.items():
        assert sym in diofant.__dict__ or sym == 'Matrix'
        assert mat in mpmath.__dict__


@with_numpy
def test_numpy_transl():
    for sym, nump in NUMPY_TRANSLATIONS.items():
        assert sym in diofant.__dict__
        assert nump in numpy.__dict__


@with_numpy
def test_numpy_translation_abs():
    f = lambdify(x, abs(x), 'numpy')
    assert f(-1) == 1
    assert f(1) == 1


# ================= Test some functions ============================


def test_exponentiation():
    f = lambdify(x, x**2)
    assert f(-1) == 1
    assert f(0) == 0
    assert f(1) == 1
    assert f(-2) == 4
    assert f(2) == 4
    assert f(2.5) == 6.25


def test_sqrt():
    f = lambdify(x, sqrt(x))
    assert f(0) == 0.0
    assert f(1) == 1.0
    assert f(4) == 2.0
    assert abs(f(2) - 1.414) < 0.001
    assert f(6.25) == 2.5


def test_trig():
    f = lambdify([x], [cos(x), sin(x)], 'math')
    d = f(pi)
    prec = 1e-11
    assert -prec < d[0] + 1 < prec
    assert -prec < d[1] < prec
    d = f(3.14159)
    prec = 1e-5
    assert -prec < d[0] + 1 < prec
    assert -prec < d[1] < prec

# ================= Test vectors ===================================


def test_vector_simple():
    f = lambdify((x, y, z), (z, y, x))
    assert f(3, 2, 1) == (1, 2, 3)
    assert f(1.0, 2.0, 3.0) == (3.0, 2.0, 1.0)
    # make sure correct number of args required
    pytest.raises(TypeError, lambda: f(0))


def test_vector_discontinuous():
    f = lambdify(x, (-1/x, 1/x))
    pytest.raises(ZeroDivisionError, lambda: f(0))
    assert f(1) == (-1.0, 1.0)
    assert f(2) == (-0.5, 0.5)
    assert f(-2) == (0.5, -0.5)


def test_trig_symbolic():
    f = lambdify([x], [cos(x), sin(x)], 'math')
    d = f(pi)
    assert abs(d[0] + 1) < 0.0001
    assert abs(d[1] - 0) < 0.0001


def test_trig_float():
    f = lambdify([x], [cos(x), sin(x)])
    d = f(3.14159)
    assert abs(d[0] + 1) < 0.0001
    assert abs(d[1] - 0) < 0.0001


def test_docs():
    f = lambdify(x, x**2)
    assert f(2) == 4
    f = lambdify([x, y, z], [z, y, x])
    assert f(1, 2, 3) == [3, 2, 1]
    f = lambdify(x, sqrt(x))
    assert f(4) == 2.0
    f = lambdify((x, y), sin(x*y)**2)
    assert f(0, 5) == 0


def test_math():
    f = lambdify((x, y), sin(x), modules='math')
    assert f(0, 5) == 0


def test_sin():
    f = lambdify(x, sin(x)**2)
    assert isinstance(f(2), float)
    f = lambdify(x, sin(x)**2, modules='math')
    assert isinstance(f(2), float)


def test_matrix():
    A = Matrix([[x, x*y], [sin(z) + 4, x**z]])
    sol = Matrix([[1, 2], [sin(3) + 4, 1]])
    f = lambdify((x, y, z), A, modules='diofant')
    assert f(1, 2, 3) == sol
    f = lambdify((x, y, z), (A, [A]), modules='diofant')
    assert f(1, 2, 3) == (sol, [sol])
    J = Matrix((x, x + y)).jacobian((x, y))
    v = Matrix((x, y))
    sol = Matrix([[1, 0], [1, 1]])
    assert lambdify(v, J, modules='diofant')(1, 2) == sol
    assert lambdify(v.T, J, modules='diofant')(1, 2) == sol


@with_numpy
def test_numpy_matrix():
    A = Matrix([[x, x*y], [sin(z) + 4, x**z]])
    sol_arr = numpy.array([[1, 2], [numpy.sin(3) + 4, 1]])
    # Lambdify array first, to ensure return to array as default
    f = lambdify((x, y, z), A, ['numpy'])
    numpy.testing.assert_allclose(f(1, 2, 3), sol_arr)
    # Check that the types are arrays and matrices
    assert isinstance(f(1, 2, 3), numpy.ndarray)


@with_numpy
def test_numpy_transpose():
    A = Matrix([[1, x], [0, 1]])
    f = lambdify(x, A.T, modules='numpy')
    numpy.testing.assert_array_equal(f(2), numpy.array([[1, 0], [2, 1]]))


@with_numpy
def test_numpy_inverse():
    A = Matrix([[1, x], [0, 1]])
    f = lambdify(x, A**-1, modules='numpy')
    numpy.testing.assert_array_equal(f(2), numpy.array([[1, -2], [0,  1]]))


@with_numpy
def test_numpy_old_matrix():
    A = Matrix([[x, x*y], [sin(z) + 4, x**z]])
    sol_arr = numpy.array([[1, 2], [numpy.sin(3) + 4, 1]])
    f = lambdify((x, y, z), A, [{'ImmutableMatrix': numpy.array}, 'numpy'])
    numpy.testing.assert_allclose(f(1, 2, 3), sol_arr)
    assert isinstance(f(1, 2, 3), numpy.ndarray)


@with_numpy
@pytest.mark.filterwarnings('ignore::RuntimeWarning')
def test_sympyissue_11306():
    p = Piecewise((1/x, y < -1), (x, y <= 1), (1/x, True))
    lambdify([x, y], p, modules='numpy')(0, 1)


@with_numpy
def test_numpy_piecewise():
    pieces = Piecewise((x, x < 3), (x**2, x > 5), (0, True))
    f = lambdify(x, pieces, modules='numpy')
    numpy.testing.assert_array_equal(f(numpy.arange(10)),
                                     numpy.array([0, 1, 2, 0, 0, 0, 36, 49, 64, 81]))
    # If we evaluate somewhere all conditions are False, we should get back NaN
    nodef_func = lambdify(x, Piecewise((x, x > 0), (-x, x < 0)))
    numpy.testing.assert_array_equal(nodef_func(numpy.array([-1, 0, 1])),
                                     numpy.array([1, numpy.nan, 1]))


@with_numpy
def test_numpy_logical_ops():
    and_func = lambdify((x, y), And(x, y), modules='numpy')
    or_func = lambdify((x, y), Or(x, y), modules='numpy')
    not_func = lambdify(x, Not(x), modules='numpy')
    arr1 = numpy.array([True, True])
    arr2 = numpy.array([False, True])
    numpy.testing.assert_array_equal(and_func(arr1, arr2), numpy.array([False, True]))
    numpy.testing.assert_array_equal(or_func(arr1, arr2), numpy.array([True, True]))
    numpy.testing.assert_array_equal(not_func(arr2), numpy.array([True, False]))


@with_numpy
def test_numpy_matmul():
    xmat = Matrix([[x, y], [z, 1+z]])
    ymat = Matrix([[x**2], [abs(x)]])
    mat_func = lambdify((x, y, z), xmat*ymat, modules='numpy')
    numpy.testing.assert_array_equal(mat_func(0.5, 3, 4), numpy.array([[1.625], [3.5]]))
    numpy.testing.assert_array_equal(mat_func(-0.5, 3, 4), numpy.array([[1.375], [3.5]]))
    # Multiple matrices chained together in multiplication
    f = lambdify((x, y, z), xmat*xmat*xmat, modules='numpy')
    numpy.testing.assert_array_equal(f(0.5, 3, 4), numpy.array([[72.125, 119.25],
                                                                [159, 251]]))


def test_integral():
    f = Lambda(x, exp(-x**2))
    l = lambdify(x, Integral(f(x), (x, -oo, oo)), modules='diofant')
    assert l(x) == Integral(exp(-x**2), (x, -oo, oo))

# ================= Test symbolic ==================================


def test_sym_single_arg():
    f = lambdify(x, x * y)
    assert f(z) == z * y


def test_sym_list_args():
    f = lambdify([x, y], x + y + z)
    assert f(1, 2) == 3 + z


def test_sym_integral():
    f = Lambda(x, exp(-x**2))
    l = lambdify(x, Integral(f(x), (x, -oo, oo)), modules='diofant')
    assert l(y).doit() == sqrt(pi)


def test_namespace_order():
    # lambdify had a bug, such that module dictionaries or cached module
    # dictionaries would pull earlier namespaces into themselves.
    # Because the module dictionaries form the namespace of the
    # generated lambda, this meant that the behavior of a previously
    # generated lambda function could change as a result of later calls
    # to lambdify.
    n1 = {'f': lambda x: 'first f'}
    n2 = {'f': lambda x: 'second f',
          'g': lambda x: 'function g'}
    f = diofant.Function('f')
    g = diofant.Function('g')
    if1 = lambdify(x, f(x), modules=(n1, 'diofant'))
    assert if1(1) == 'first f'
    if2 = lambdify(x, g(x), modules=(n2, 'diofant'))
    assert if2(1) == 'function g'
    # previously gave 'second f'
    assert if1(1) == 'first f'


def test_imps():
    # Here we check if the default returned functions are anonymous - in
    # the sense that we can have more than one function with the same name
    f = implemented_function('f', lambda x: 2*x)
    g = implemented_function('f', math.sqrt)
    l1 = lambdify(x, f(x))
    l2 = lambdify(x, g(x))
    assert str(f(x)) == str(g(x))
    assert l1(3) == 6
    assert l2(3) == math.sqrt(3)
    # check that we can pass in a Function as input
    func = diofant.Function('myfunc')
    assert not hasattr(func, '_imp_')
    my_f = implemented_function(func, lambda x: 2*x)
    assert hasattr(func, '_imp_')
    assert hasattr(my_f, '_imp_')
    # Error for functions with same name and different implementation
    f2 = implemented_function('f', lambda x: x + 101)
    pytest.raises(ValueError, lambda: lambdify(x, f(f2(x))))


def test_imps_errors():
    # Test errors that implemented functions can return, and still be
    # able to form expressions.  See issue sympy/sympy#10810.
    for val, error_class in itertools.product((0, 0., 2, 2.0),
                                              (AttributeError, TypeError,
                                               ValueError)):

        def myfunc(a):
            if a == 0:
                raise error_class
            return 1

        f = implemented_function('f', myfunc)
        expr = f(val)
        assert expr == f(val)


def test_imps_wrong_args():
    pytest.raises(ValueError, lambda: implemented_function(sin, lambda x: x))


def test_lambdify_imps():
    # Test lambdify with implemented functions
    # first test basic (diofant) lambdify
    f = diofant.cos
    assert lambdify(x, f(x))(0) == 1
    assert lambdify(x, 1 + f(x))(0) == 2
    assert lambdify((x, y), y + f(x))(0, 1) == 2
    # make an implemented function and test
    f = implemented_function('f', lambda x: x + 100)
    assert lambdify(x, f(x))(0) == 100
    assert lambdify(x, 1 + f(x))(0) == 101
    assert lambdify((x, y), y + f(x))(0, 1) == 101
    # Can also handle tuples, lists, dicts as expressions
    lam = lambdify(x, (f(x), x))
    assert lam(3) == (103, 3)
    lam = lambdify(x, [f(x), x])
    assert lam(3) == [103, 3]
    lam = lambdify(x, [f(x), (f(x), x)])
    assert lam(3) == [103, (103, 3)]
    lam = lambdify(x, {f(x): x})
    assert lam(3) == {103: 3}
    lam = lambdify(x, {f(x): x})
    assert lam(3) == {103: 3}
    lam = lambdify(x, {x: f(x)})
    assert lam(3) == {3: 103}
    # Check that imp preferred to other namespaces by default
    d = {'f': lambda x: x + 99}
    lam = lambdify(x, f(x), d)
    assert lam(3) == 103
    # Unless flag passed
    lam = lambdify(x, f(x), d, use_imps=False)
    assert lam(3) == 102


def test_dummification():
    F = Function('F')
    G = Function('G')
    # "\alpha" is not a valid python variable name
    # lambdify should sub in a dummy for it, and return
    # without a syntax error
    alpha = symbols(r'\alpha')
    some_expr = 2 * F(t)**2 / G(t)
    lam = lambdify((F(t), G(t)), some_expr)
    assert lam(3, 9) == 2
    lam = lambdify(sin(t), 2 * sin(t)**2)
    assert lam(F(t)) == 2 * F(t)**2
    # Test that \alpha was properly dummified
    lam = lambdify((alpha, t), 2*alpha + t)
    assert lam(2, 1) == 5
    pytest.raises(SyntaxError, lambda: lambdify(F(t) * G(t), F(t) * G(t) + 5))
    pytest.raises(SyntaxError, lambda: lambdify(2 * F(t), 2 * F(t) + 5))
    pytest.raises(SyntaxError, lambda: lambdify(2 * F(t), 4 * F(t) + 5))


def test_python_keywords():
    # Test for issue sympy/sympy#7452. The automatic dummification should ensure use of
    # Python reserved keywords as symbol names will create valid lambda
    # functions. This is an additional regression test.
    python_if = symbols('if')
    expr = python_if / 2
    f = lambdify(python_if, expr)
    assert f(4.0) == 2.0


def test_lambdify_docstring():
    func = lambdify((w, x, y, z), w + x + y + z)
    assert func.__doc__ == (
        'Created with lambdify. Signature:\n\n'
        'func(w, x, y, z)\n\n'
        'Expression:\n\n'
        'w + x + y + z')
    syms = symbols('a1:26')
    func = lambdify(syms, sum(syms))
    assert func.__doc__ == (
        'Created with lambdify. Signature:\n\n'
        'func(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15,\n'
        '        a16, a17, a18, a19, a20, a21, a22, a23, a24, a25)\n\n'
        'Expression:\n\n'
        'a1 + a10 + a11 + a12 + a13 + a14 + a15 + a16 + a17 + a18 + a19 + a2 + a20 +...')


# ================= Test special printers ==========================


def test_special_printers():
    class IntervalPrinter(LambdaPrinter):
        """Use ``lambda`` printer but print numbers as ``mpi`` intervals."""

        def _print_Integer(self, expr):
            return f"mpi('{super()._print_Integer(expr)}')"

        def _print_Rational(self, expr):
            return f"mpi('{super()._print_Rational(expr)}')"

    def intervalrepr(expr):
        return IntervalPrinter().doprint(expr)

    expr = diofant.sqrt(diofant.sqrt(2) + diofant.sqrt(3)) + diofant.Rational(1, 2)

    func0 = lambdify((), expr, modules='mpmath', printer=intervalrepr)
    func1 = lambdify((), expr, modules='mpmath', printer=IntervalPrinter)
    func2 = lambdify((), expr, modules='mpmath', printer=IntervalPrinter())

    mpi = type(mpmath.mpi(1, 2))

    assert isinstance(func0(), mpi)
    assert isinstance(func1(), mpi)
    assert isinstance(func2(), mpi)


def test_true_false():
    # We want exact is comparison here, not just ==
    assert lambdify([], true)() is True
    assert lambdify([], false)() is False


def test_ITE():
    assert lambdify((x, y, z), ITE(x, y, z))(True, 5, 3) == 5
    assert lambdify((x, y, z), ITE(x, y, z))(False, 5, 3) == 3


def test_Min_Max():
    # see sympy/sympy#10375
    assert lambdify((x, y, z), Min(x, y, z))(1, 2, 3) == 1
    assert lambdify((x, y, z), Max(x, y, z))(1, 2, 3) == 3


def test_sympyissue_12092():
    f = implemented_function('f', lambda x: x**2)
    assert f(f(2)).evalf() == Float(16)


def test_sympyissue_23224():
    f = lambdify([], (1,))
    assert f() == (1,)
