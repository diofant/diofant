# This testfile tests Diofant <-> NumPy compatibility

# Don't test any Diofant features here. Just pure interaction with NumPy.
# Always write regular Diofant tests for anything, that can be tested in pure
# Python (without numpy). Here we test everything, that a user may need when
# using Diofant with NumPy

import mpmath
import pytest

import diofant
from diofant import (Float, Integer, Matrix, MatrixSymbol, Rational, Symbol,
                     lambdify, list2numpy, matrix2numpy, sin, symarray,
                     symbols, sympify)
from diofant.abc import x, y, z
from diofant.matrices.expressions.matexpr import MatrixElement
from diofant.utilities.decorator import conserve_mpmath_dps
from diofant.utilities.lambdify import NUMPY_TRANSLATIONS


__all__ = ()

numpy = pytest.importorskip('numpy')


# first, systematically check, that all operations are implemented and don't
# raise an exception


def test_systematic_basic():
    def s(diofant_object, numpy_array):
        # pylint: disable=pointless-statement
        diofant_object + numpy_array
        numpy_array + diofant_object
        diofant_object - numpy_array
        numpy_array - diofant_object
        diofant_object * numpy_array
        numpy_array * diofant_object
        diofant_object / numpy_array
        numpy_array / diofant_object
        diofant_object ** numpy_array
        numpy_array ** diofant_object
    x = Symbol('x')
    y = Symbol('y')
    diofant_objs = [
        Rational(2, 3),
        Float('1.3'),
        x,
        y,
        pow(x, y)*y,
        Integer(5),
        Float(5.5),
    ]
    numpy_objs = [
        numpy.array([1]),
        numpy.array([3, 8, -1]),
        numpy.array([x, x**2, Integer(5)]),
        numpy.array([x/y*sin(y), 5, Integer(5)]),
    ]
    for x in diofant_objs:
        for y in numpy_objs:
            s(x, y)


# now some random tests, that test particular problems and that also
# check that the results of the operations are correct

def test_basics():
    one = Integer(1)
    zero = Integer(0)
    assert numpy.array(1) == numpy.array(one)
    assert numpy.array([one]) == numpy.array([one])
    assert numpy.array([x]) == numpy.array([x])
    assert numpy.array(x) == numpy.array(Symbol('x'))
    assert numpy.array(one + x) == numpy.array(1 + x)

    X = numpy.array([one, zero, zero])
    assert (X == numpy.array([one, zero, zero])).all()
    assert (X == numpy.array([one, 0, 0])).all()


def test_arrays():
    X = numpy.array([Symbol('a') + Rational(1, 2)])
    Y = X + X
    assert Y == numpy.array([1 + 2*Symbol('a')])
    Y = Y + 1
    assert Y == numpy.array([2 + 2*Symbol('a')])
    Y = X - X
    assert Y == numpy.array([0])


def test_conversion1():
    a = list2numpy([x**2, x])
    # looks like an array?
    assert isinstance(a, numpy.ndarray)
    assert a[0] == x**2
    assert a[1] == x
    assert len(a) == 2
    # yes, it's the array


def test_conversion2():
    a = 2*list2numpy([x**2, x])
    b = list2numpy([2*x**2, 2*x])
    assert (a == b).all()

    X = list2numpy([Symbol('a') + Rational(1, 2)])
    Y = X + X
    assert Y == numpy.array([1 + 2*Symbol('a')])
    Y = Y + 1
    assert Y == numpy.array([2 + 2*Symbol('a')])
    Y = X - X
    assert Y == numpy.array([0])


def test_list2numpy():
    assert (numpy.array([x**2, x]) == list2numpy([x**2, x])).all()


def test_Matrix1():
    m = Matrix([[x, x**2], [5, 2/x]])
    assert (numpy.array(m.subs({x: 2})) == numpy.array([[2, 4], [5, 1]])).all()
    m = Matrix([[sin(x), x**2], [5, 2/x]])
    assert (numpy.array(m.subs({x: 2})) == numpy.array([[sin(2), 4], [5, 1]])).all()


def test_Matrix2():
    a = numpy.array([[2, 4], [5, 1]])
    assert Matrix(a) == Matrix([[2, 4], [5, 1]])
    assert Matrix(a) != Matrix([[2, 4], [5, 2]])
    a = numpy.array([[sin(2), 4], [5, 1]])
    assert Matrix(a) == Matrix([[sin(2), 4], [5, 1]])
    assert Matrix(a) != Matrix([[sin(0), 4], [5, 1]])


def test_Matrix_sum():
    M = Matrix([[1, 2, 3], [x, y, x], [2*y, -50, z*x]])
    m = numpy.array([[2, 3, 4], [x, 5, 6], [x, y, z**2]])
    assert M + m == Matrix([[3, 5, 7], [2*x, y + 5, x + 6], [2*y + x, y - 50, z*x + z**2]])
    assert m + M == Matrix([[3, 5, 7], [2*x, y + 5, x + 6], [2*y + x, y - 50, z*x + z**2]])
    assert M + m == M.add(m)


def test_Matrix_mul():
    M = Matrix([[1, 2, 3], [x, y, x]])
    m = numpy.array([[2, 4], [x, 6], [x, z**2]])
    assert M*m == Matrix([
        [         2 + 5*x,        16 + 3*z**2],
        [2*x + x*y + x**2, 4*x + 6*y + x*z**2],
    ])

    assert m*M == Matrix([
        [   2 + 4*x,      4 + 4*y,      6 + 4*x],
        [       7*x,    2*x + 6*y,          9*x],
        [x + x*z**2, 2*x + y*z**2, 3*x + x*z**2],
    ])
    a = numpy.array([2])
    assert a[0] * M == 2 * M
    assert M * a[0] == 2 * M


def test_Matrix_numpy_array():
    class MatArray:
        def __array__(self):
            return numpy.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    matarr = MatArray()
    assert Matrix(matarr) == Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])


def test_matrix2numpy():
    a = matrix2numpy(Matrix([[1, x**2], [3*sin(x), 0]]))
    assert isinstance(a, numpy.ndarray)
    assert a.shape == (2, 2)
    assert a[0, 0] == 1
    assert a[0, 1] == x**2
    assert a[1, 0] == 3*sin(x)
    assert a[1, 1] == 0


def test_matrix2numpy_conversion():
    a = Matrix([[1, 2, sin(x)], [x**2, x, Rational(1, 2)]])
    b = numpy.array([[1, 2, sin(x)], [x**2, x, Rational(1, 2)]])
    assert (matrix2numpy(a) == b).all()
    assert matrix2numpy(a).dtype == numpy.dtype('object')

    c = matrix2numpy(Matrix([[1, 2], [10, 20]]), dtype='int8')
    d = matrix2numpy(Matrix([[1, 2], [10, 20]]), dtype='float64')
    assert c.dtype == numpy.dtype('int8')
    assert d.dtype == numpy.dtype('float64')


def test_sympyissue_3728():
    assert (Rational(1, 2)*numpy.array([2*x, 0]) == numpy.array([x, 0])).all()
    assert (Rational(1, 2) + numpy.array(
        [2*x, 0]) == numpy.array([2*x + Rational(1, 2), Rational(1, 2)])).all()
    assert (Float('0.5')*numpy.array([2*x, 0]) == numpy.array([Float('1.0')*x, 0])).all()
    assert (Float('0.5') + numpy.array(
        [2*x, 0]) == numpy.array([2*x + Float('0.5'), Float('0.5')])).all()


@conserve_mpmath_dps
def test_lambdify():
    mpmath.mp.dps = 16
    sin02 = mpmath.mpf('0.198669330795061215459412627')
    f = lambdify(x, sin(x), 'numpy')
    prec = 1e-15
    assert -prec < f(0.2) - sin02 < prec
    with pytest.raises(TypeError):
        f(x)


def test_lambdify_numpy_matrix():
    f = lambdify(x, Matrix([[x, 2*x], [1, 2]]), [{'ImmutableMatrix': numpy.array}, 'numpy'])
    assert (f(1) == numpy.array([[1, 2], [1, 2]])).all()


def test_lambdify_matrix_multi_input():
    M = Matrix([[x**2, x*y, x*z],
                [y*x, y**2, y*z],
                [z*x, z*y, z**2]])
    f = lambdify((x, y, z), M, [{'ImmutableMatrix': numpy.array}, 'numpy'])

    xh, yh, zh = 1.0, 2.0, 3.0
    expected = numpy.array([[xh**2, xh*yh, xh*zh],
                            [yh*xh, yh**2, yh*zh],
                            [zh*xh, zh*yh, zh**2]])
    actual = f(xh, yh, zh)
    assert numpy.allclose(actual, expected)


def test_lambdify_transl():
    for sym, mat in NUMPY_TRANSLATIONS.items():
        assert sym in dir(diofant)
        assert mat in dir(numpy)


def test_symarray():
    """Test creation of numpy arrays of diofant symbols."""
    syms = symbols('_0,_1,_2')
    s1 = symarray('', 3)
    s2 = symarray('', 3)
    numpy.testing.assert_array_equal(s1, numpy.array(syms, dtype=object))
    assert s1[0] == s2[0]

    a = symarray('a', 3)
    b = symarray('b', 3)
    assert not a[0] == b[0]

    asyms = symbols('a_0,a_1,a_2')
    numpy.testing.assert_array_equal(a, numpy.array(asyms, dtype=object))

    # Multidimensional checks
    a2d = symarray('a', (2, 3))
    assert a2d.shape == (2, 3)
    a00, a12 = symbols('a_0_0,a_1_2')
    assert a2d[0, 0] == a00
    assert a2d[1, 2] == a12

    a3d = symarray('a', (2, 3, 2))
    assert a3d.shape == (2, 3, 2)
    a000, a120, a121 = symbols('a_0_0_0,a_1_2_0,a_1_2_1')
    assert a3d[0, 0, 0] == a000
    assert a3d[1, 2, 0] == a120
    assert a3d[1, 2, 1] == a121


def test_vectorize():
    assert (numpy.vectorize(
        sin)([1, 2, 3]) == numpy.array([sin(1), sin(2), sin(3)])).all()


def test_array_coeersion():
    A = MatrixSymbol('A', 2, 2)
    assert numpy.array(A)[1][0] == MatrixElement(A, 1, 0)


def test_from_ndarray():
    # See sympy/sympy#7465
    assert Matrix(numpy.array([1, 2, 3])) == Matrix([1, 2, 3])
    assert Matrix(numpy.array([[1, 2, 3]])) == Matrix([[1, 2, 3]])
    assert Matrix(numpy.array([[1, 2, 3], [4, 5, 6]])) == \
        Matrix([[1, 2, 3], [4, 5, 6]])
    assert Matrix(numpy.array([x, y, z])) == Matrix([x, y, z])
    pytest.raises(NotImplementedError, lambda: Matrix(numpy.array([[
        [1, 2], [3, 4]], [[5, 6], [7, 8]]])))


def test_sympify():
    assert sympify(numpy.float128(1.1)) == Float(1.1)
