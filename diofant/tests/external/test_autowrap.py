import os
import tempfile

import pytest

import diofant
from diofant import Eq, Idx, IndexedBase, sin, symbols
from diofant.abc import a, b, c
from diofant.external import import_module
from diofant.utilities.autowrap import (CodeWrapError, F2PyCodeWrapper,
                                        autowrap, ufuncify)


__all__ = ()

numpy = import_module('numpy', min_module_version='1.6.1')
with_numpy = pytest.mark.skipif(numpy is None,
                                reason="Couldn't import numpy.")

Cython = import_module('Cython', min_module_version='0.15.1')
with_cython = pytest.mark.skipif(Cython is None,
                                 reason="Couldn't import Cython.")

f2py = import_module('numpy.f2py', import__kwargs={'fromlist': ['f2py']})
with_f2py = pytest.mark.skipif(f2py is None, reason="Couldn't run f2py.")

f2pyworks = False
if f2py:
    try:
        autowrap(symbols('x'), 'f95', 'f2py')
    except (CodeWrapError, ImportError, OSError):
        f2pyworks = False
    else:
        f2pyworks = True

n, m, d = symbols('n m d', integer=True)
A, B, C = symbols('A B C', cls=IndexedBase)
i = Idx('i', m)
j = Idx('j', n)
k = Idx('k', d)


#
# test runners used by several language-backend combinations
#

def runtest_autowrap_twice(language, backend):
    f = autowrap((((a + b)/c)**5).expand(), language, backend)
    g = autowrap((((a + b)/c)**4).expand(), language, backend)

    # check that autowrap updates the module name.  Else, g gives the same as f
    assert f(1, -2, 1) == -1.0
    assert g(1, -2, 1) == 1.0


@with_numpy
def runtest_autowrap_trace(language, backend):
    trace = autowrap(A[i, i], language, backend)
    assert trace(numpy.eye(100)) == 100


@with_numpy
def runtest_autowrap_matrix_vector(language, backend):
    x, y = symbols('x y', cls=IndexedBase)
    expr = Eq(y[i], A[i, j]*x[j])
    mv = autowrap(expr, language, backend)

    # compare with numpy's dot product
    M = numpy.random.rand(10, 20)
    x = numpy.random.rand(20)
    y = numpy.dot(M, x)
    assert numpy.sum(numpy.abs(y - mv(M, x))) < 1e-13


@with_numpy
def runtest_autowrap_matrix_matrix(language, backend):
    expr = Eq(C[i, j], A[i, k]*B[k, j])
    matmat = autowrap(expr, language, backend)

    # compare with numpy's dot product
    M1 = numpy.random.rand(10, 20)
    M2 = numpy.random.rand(20, 15)
    M3 = numpy.dot(M1, M2)
    assert numpy.sum(numpy.abs(M3 - matmat(M1, M2))) < 1e-13


@with_numpy
def runtest_ufuncify(language, backend):
    a, b, c = symbols('a b c')

    if backend == 'numpy':
        pytest.raises(ValueError, lambda: ufuncify((a, b), Eq(b, a + b), backend=backend))
        pytest.raises(ValueError, lambda: ufuncify((a, b, c), [Eq(b, a), Eq(c, a**2)], backend=backend))

    helpers = [('helper', a*b, (a, b)), ('spam', sin(a), (a,))]

    fabc = ufuncify((a, b, c), a*b + c, backend=backend)
    fabc_2 = ufuncify((a, b, c), a*b + c, backend=backend, helpers=helpers)
    facb = ufuncify((a, c, b), a*b + c, backend=backend)
    grid = numpy.linspace(-2, 2, 50)
    b = numpy.linspace(-5, 4, 50)
    c = numpy.linspace(-1, 1, 50)
    expected = grid*b + c
    numpy.testing.assert_allclose(fabc(grid, b, c), expected)
    numpy.testing.assert_allclose(fabc_2(grid, b, c), expected)
    numpy.testing.assert_allclose(facb(grid, c, b), expected)


def runtest_autowrap_helpers(language, backend):
    # issue sympy/sympy#10274
    expr = (a - b + c)**13
    tmp = tempfile.mkdtemp()
    f = autowrap(expr, language, backend, tempdir=tmp, helpers=[('helper', a - b + c, (a, b, c))])
    assert f(1, 1, 1) == 1

    for file in os.listdir(tmp):
        if file.startswith('wrapped_code_') and file.endswith('.c'):
            with open(tmp + '/' + file, encoding='utf-8') as f:
                assert f.read() == ('/******************************************************************************\n'
                                    ' *' + ('Code generated with diofant ' + diofant.__version__).center(76) + '*\n'
                                    ' *                                                                            *\n'
                                    ' *         See https://diofant.readthedocs.io/ for more information.          *\n'
                                    ' *                                                                            *\n'
                                    " *                      This file is part of 'autowrap'                       *\n"
                                    ' ******************************************************************************/\n'
                                    '#include ' + '"' + file[:-1] + 'h"' + '\n'
                                    '#include <math.h>\n'
                                    '\n'
                                    'double helper(double a, double b, double c) {\n'
                                    '\n'
                                    '   double helper_result;\n'
                                    '   helper_result = a - b + c;\n'
                                    '   return helper_result;\n'
                                    '\n'
                                    '}\n'
                                    '\n'
                                    'double autofunc(double a, double b, double c) {\n'
                                    '\n'
                                    '   double autofunc_result;\n'
                                    '   autofunc_result = pow(helper(a, b, c), 13);\n'
                                    '   return autofunc_result;\n'
                                    '\n'
                                    '}\n')

#
# tests of language-backend combinations
#

# f2py


@with_f2py
def test_wrap_twice_f95_f2py():
    runtest_autowrap_twice('f95', 'f2py')

    f = ufuncify(a, 2*a, backend='f2py')
    f2 = ufuncify((a,), 2*a, backend='f2py')
    f3 = ufuncify((a,), 2*a, backend='f2py', language='f95')

    grid = numpy.linspace(-2, 2, 50)
    numpy.testing.assert_allclose(f(grid), f2(grid))
    numpy.testing.assert_allclose(f(grid), f3(grid))


@with_f2py
def test_autowrap_trace_f95_f2py():
    runtest_autowrap_trace('f95', 'f2py')


@with_f2py
def test_autowrap_matrix_vector_f95_f2py():
    runtest_autowrap_matrix_vector('f95', 'f2py')


@with_f2py
def test_autowrap_matrix_matrix_f95_f2py():
    runtest_autowrap_matrix_matrix('f95', 'f2py')


@with_f2py
def test_ufuncify_f95_f2py():
    runtest_ufuncify('f95', 'f2py')


@with_numpy
def test_ufuncify_f95_numpy():
    runtest_ufuncify('f95', 'numpy')


@with_f2py
def test_autowrap_verbose(capsys):
    f = autowrap((((a + b)/c)**5).expand(), backend='f2py', verbose=True)

    assert capsys.readouterr().out.find('running build') == 2
    assert f(1, -2, 1) == -1.0


@with_f2py
def test_autowrap_command_err(monkeypatch):
    monkeypatch.setattr(F2PyCodeWrapper, 'command', ['/bin/false'])
    pytest.raises(CodeWrapError, lambda: autowrap(1, backend='f2py'))


# Cython

@with_cython
def test_wrap_twice_c_cython():
    runtest_autowrap_twice('C', 'cython')


@with_cython
@with_numpy
def test_autowrap_trace_C_Cython():
    runtest_autowrap_trace('C', 'cython')


@with_cython
@with_numpy
def test_autowrap_matrix_vector_C_cython():
    runtest_autowrap_matrix_vector('C', 'cython')


@with_cython
@with_numpy
def test_autowrap_matrix_matrix_C_cython():
    runtest_autowrap_matrix_matrix('C', 'cython')


@with_cython
@with_numpy
def test_ufuncify_C_Cython():
    runtest_ufuncify('C', 'cython')


@with_cython
def test_autowrap_helpers_C_cython():
    runtest_autowrap_helpers('C', 'cython')


# Numpy

@with_cython
@with_numpy
def test_ufuncify_numpy():
    # This test doesn't use Cython, but if Cython works, then there is a valid
    # C compiler, which is needed.
    runtest_ufuncify('C', 'numpy')
