# Tests that require installed backends go into
# diofant/test_external/test_autowrap

import io
import os
import shutil
import tempfile

import pytest

from diofant import Eq
from diofant.abc import x, y, z
from diofant.utilities.autowrap import (CodeWrapper, CythonCodeWrapper,
                                        UfuncifyCodeWrapper, autowrap,
                                        binary_function)
from diofant.utilities.codegen import (CCodeGen, CodeGenArgumentListError,
                                       make_routine)


__all__ = ()


def get_string(dump_fn, routines, prefix='file'):
    # Wrapper for dump_fn. dump_fn writes its results to a stream object and
    # this wrapper returns the contents of that stream as a string. This
    # auxiliary function is used by many tests below.
    #
    # The header and the empty lines are not generator to facilitate the
    # testing of the output.
    output = io.StringIO()
    dump_fn(routines, output, prefix)
    source = output.getvalue()
    output.close()
    return source


def test_cython_wrapper_scalar_function():
    expr = (x + y)*z
    routine = make_routine('test', expr)
    code_gen = CythonCodeWrapper(CCodeGen())
    source = get_string(code_gen.dump_pyx, [routine])
    expected = (
        "cdef extern from 'file.h':\n"
        '    double test(double x, double y, double z)\n'
        '\n'
        'def test_c(double x, double y, double z):\n'
        '\n'
        '    return test(x, y, z)')
    assert source == expected


def test_cython_wrapper_outarg():
    code_gen = CythonCodeWrapper(CCodeGen())

    routine = make_routine('test', Eq(z, x + y))
    source = get_string(code_gen.dump_pyx, [routine])
    expected = (
        "cdef extern from 'file.h':\n"
        '    void test(double x, double y, double *z)\n'
        '\n'
        'def test_c(double x, double y):\n'
        '\n'
        '    cdef double z = 0\n'
        '    test(x, y, &z)\n'
        '    return z')
    assert source == expected


def test_cython_wrapper_inoutarg():
    code_gen = CythonCodeWrapper(CCodeGen())
    routine = make_routine('test', Eq(z, x + y + z))
    source = get_string(code_gen.dump_pyx, [routine])
    expected = (
        "cdef extern from 'file.h':\n"
        '    void test(double x, double y, double *z)\n'
        '\n'
        'def test_c(double x, double y, double z):\n'
        '\n'
        '    test(x, y, &z)\n'
        '    return z')
    assert source == expected


def test_autowrap_dummy():
    # Uses DummyWrapper to test that codegen works as expected

    tempdir = tempfile.mkstemp()[-1]
    os.unlink(tempdir)
    f = autowrap(x + y, backend='dummy', tempdir=tempdir)
    assert f() == str(x + y)
    assert f.args == 'x, y'
    assert f.returns == 'nameless'
    f = autowrap(Eq(z, x + y), backend='dummy')
    assert f() == str(x + y)
    assert f.args == 'x, y'
    assert f.returns == 'z'
    f = autowrap(Eq(z, x + y + z), backend='dummy')
    assert f() == str(x + y + z)
    assert f.args == 'x, y, z'
    assert f.returns == 'z'

    e = x + y

    pytest.raises(ValueError, lambda: autowrap(e, backend='spam'))
    pytest.raises(ValueError, lambda: autowrap(e, backend='spam', language='C'))
    pytest.raises(ValueError, lambda: autowrap(e, language='spam'))


def test_autowrap_args():
    pytest.raises(CodeGenArgumentListError,
                  lambda: autowrap(Eq(z, x + y), backend='dummy', args=(x,)))
    f = autowrap(Eq(z, x + y), backend='dummy', args=(y, x))
    assert f() == str(x + y)
    assert f.args == 'y, x'
    assert f.returns == 'z'

    pytest.raises(CodeGenArgumentListError,
                  lambda: autowrap(Eq(z, x + y + z), backend='dummy', args=(x, y)))
    f = autowrap(Eq(z, x + y + z), backend='dummy', args=(y, x, z))
    assert f() == str(x + y + z)
    assert f.args == 'y, x, z'
    assert f.returns == 'z'

    f = autowrap(Eq(z, x + y + z), backend='dummy', args=(y, x, z))
    assert f() == str(x + y + z)
    assert f.args == 'y, x, z'
    assert f.returns == 'z'


def test_autowrap_store_files():
    tmp = tempfile.mkdtemp()
    try:
        f = autowrap(x + y, backend='dummy', tempdir=tmp)
        assert f() == str(x + y)
        assert os.access(tmp, os.F_OK)
    finally:
        shutil.rmtree(tmp)


def test_binary_function():
    f = binary_function('f', x + y, backend='dummy')
    assert f._imp_() == str(x + y)


def test_ufuncify_source():
    code_wrapper = UfuncifyCodeWrapper(CCodeGen('ufuncify'))
    routine = make_routine('test', x + y + z)
    source = get_string(code_wrapper.dump_c, [routine])
    expected = """\
#include "Python.h"
#include "math.h"
#include "numpy/ndarraytypes.h"
#include "numpy/ufuncobject.h"
#include "numpy/halffloat.h"
#include "file.h"

static PyMethodDef wrapper_module_{num}Methods[] = {{
        {{NULL, NULL, 0, NULL}}
}};

static void test_ufunc(char **args, npy_intp *dimensions, npy_intp* steps, void* data)
{{
    npy_intp i;
    npy_intp n = dimensions[0];
    char *in0 = args[0];
    char *in1 = args[1];
    char *in2 = args[2];
    char *out1 = args[3];
    npy_intp in0_step = steps[0];
    npy_intp in1_step = steps[1];
    npy_intp in2_step = steps[2];
    npy_intp out1_step = steps[3];
    for (i = 0; i < n; i++) {{
        *((double *)out1) = test(*(double *)in0, *(double *)in1, *(double *)in2);
        in0 += in0_step;
        in1 += in1_step;
        in2 += in2_step;
        out1 += out1_step;
    }}
}}
PyUFuncGenericFunction test_funcs[1] = {{&test_ufunc}};
static char test_types[4] = {{NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE}};
static void *test_data[1] = {{NULL}};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {{
    PyModuleDef_HEAD_INIT,
    "wrapper_module_{num}",
    NULL,
    -1,
    wrapper_module_{num}Methods,
    NULL,
    NULL,
    NULL,
    NULL
}};

PyMODINIT_FUNC PyInit_wrapper_module_{num}(void)
{{
    PyObject *m, *d;
    PyObject *ufunc0;
    m = PyModule_Create(&moduledef);
    if (!m) {{
        return NULL;
    }}
    import_array();
    import_umath();
    d = PyModule_GetDict(m);
    ufunc0 = PyUFunc_FromFuncAndData(test_funcs, test_data, test_types, 1, 3, 1,
            PyUFunc_None, "wrapper_module_{num}", "Created in Diofant with Ufuncify", 0);
    PyDict_SetItemString(d, "test", ufunc0);
    Py_DECREF(ufunc0);
    return m;
}}
#else
PyMODINIT_FUNC initwrapper_module_{num}(void)
{{
    PyObject *m, *d;
    PyObject *ufunc0;
    m = Py_InitModule("wrapper_module_{num}", wrapper_module_{num}Methods);
    if (m == NULL) {{
        return;
    }}
    import_array();
    import_umath();
    d = PyModule_GetDict(m);
    ufunc0 = PyUFunc_FromFuncAndData(test_funcs, test_data, test_types, 1, 3, 1,
            PyUFunc_None, "wrapper_module_{num}", "Created in Diofant with Ufuncify", 0);
    PyDict_SetItemString(d, "test", ufunc0);
    Py_DECREF(ufunc0);
}}
#endif""".format(num=CodeWrapper._module_counter)  # noqa: SFS201
    assert source == expected
