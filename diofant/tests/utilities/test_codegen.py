import io

import pytest

from diofant import (Catalan, Dummy, Eq, Equality, Function, Idx, IndexedBase,
                     Integral, Lambda, Matrix, MatrixSymbol, acos, asin, atan,
                     atan2, besseli, ceiling, cos, cosh, erf, floor, ln, log,
                     pi, sin, sinh, sqrt, symbols, tan, tanh)
from diofant.abc import B, C, X, a, b, c, d, t, x, y, z
from diofant.utilities.codegen import (CCodeGen, CodeGenArgumentListError,
                                       CodeGenError, FCodeGen, InOutArgument,
                                       InputArgument, OutputArgument, Result,
                                       Routine, codegen, default_datatypes,
                                       make_routine)
from diofant.utilities.lambdify import implemented_function


__all__ = ()


def get_string(dump_fn, routines, prefix='file', header=False, empty=False):
    # Wrapper for dump_fn. dump_fn writes its results to a stream object and
    # this wrapper returns the contents of that stream as a string. This
    # auxiliary function is used by many tests below.
    #
    # The header and the empty lines are not generated to facilitate the
    # testing of the output.
    output = io.StringIO()
    dump_fn(routines, output, prefix, header, empty)
    source = output.getvalue()
    output.close()
    return source


def test_for_bad_arguments():
    pytest.raises(ValueError, lambda: make_routine('test', x, language='foo'))

    cg = CCodeGen()
    pytest.raises(CodeGenError, lambda: cg.write([x], prefix='test'))

    pytest.raises(CodeGenError, lambda: cg.routine('test', Eq(sin(y), x),
                                                   argument_sequence=None))

    pytest.raises(TypeError, lambda: Result([x]))
    pytest.raises(TypeError, lambda: InOutArgument(sin(y), y, y*sin(x)))
    pytest.raises(TypeError, lambda: InOutArgument(y, y, y*sin(x), 'spam'))
    pytest.raises(TypeError, lambda: InOutArgument(y, y, y*sin(x), dimensions='spam'))

    r = cg.routine('test', x, argument_sequence=None)
    pytest.raises(ValueError, lambda: Routine(r.name, r.arguments + ['spam'],
                                              r.results, r.local_vars,
                                              r.global_vars))
    pytest.raises(ValueError, lambda: Routine(r.name, r.arguments,
                                              r.results + ['spam'],
                                              r.local_vars, r.global_vars))
    pytest.raises(ValueError, lambda: Routine(r.name, r.arguments, [Result(y)],
                                              r.local_vars, r.global_vars))


def test_low_level():
    a = InOutArgument(y, y, y*sin(x), default_datatypes['float'])
    assert a.name == y
    assert a.expr == y*sin(x)
    assert a.get_datatype('C') == 'double'
    pytest.raises(CodeGenError, lambda: a.get_datatype('spam'))


def test_Routine_argument_order():
    expr = (x + y)*z
    pytest.raises(CodeGenArgumentListError, lambda: make_routine('test', expr,
                                                                 argument_sequence=[z, x]))
    pytest.raises(CodeGenArgumentListError, lambda: make_routine('test', Eq(a,
                                                                            expr), argument_sequence=[z, x, y]))
    r = make_routine('test', Eq(a, expr), argument_sequence=[z, x, a, y])
    assert [arg.name for arg in r.arguments] == [z, x, a, y]
    assert [type(arg) for arg in r.arguments] == [
        InputArgument, InputArgument, OutputArgument, InputArgument]
    r = make_routine('test', Eq(z, expr), argument_sequence=[z, x, y])
    assert [type(arg) for arg in r.arguments] == [
        InOutArgument, InputArgument, InputArgument]

    A, B = map(IndexedBase, ['A', 'B'])
    m = symbols('m', integer=True)
    i = Idx('i', m)
    r = make_routine('test', Eq(A[i], B[i]), argument_sequence=[B, A, m])
    assert [arg.name for arg in r.arguments] == [B.label, A.label, m]

    expr = Integral(x*y*z, (x, 1, 2), (y, 1, 3))
    r = make_routine('test', Eq(a, expr), argument_sequence=[z, x, a, y])
    assert [arg.name for arg in r.arguments] == [z, x, a, y]


def test_empty_c_code():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [])
    assert source == '#include "file.h"\n#include <math.h>\n'
    code_gen = CCodeGen(preprocessor_statements='')
    source = get_string(code_gen.dump_c, [])
    assert source == '#include "file.h"\n'


def test_empty_c_code_with_comment():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [], header=True)
    assert source[:82] == (
        '/******************************************************************************\n *'
    )
    #   "                    Code generated with diofant x.y.z                    "
    assert source[158:] == (                                                              '*\n'
                                                                                          ' *                                                                            *\n'
                                                                                          ' *         See https://diofant.readthedocs.io/ for more information.          *\n'
                                                                                          ' *                                                                            *\n'
                                                                                          " *                       This file is part of 'project'                       *\n"
                                                                                          ' ******************************************************************************/\n'
                                                                                          '#include "file.h"\n'
                                                                                          '#include <math.h>\n')


def test_empty_c_header():
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_h, [])
    assert source == '#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n#endif\n'


def test_simple_c_code():
    expr = (x + y)*z
    routine = make_routine('test', expr)
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'double test(double x, double y, double z) {\n'
        '   double test_result;\n'
        '   test_result = z*(x + y);\n'
        '   return test_result;\n'
        '}\n'
    )
    assert source == expected


def test_c_code_reserved_words():
    if_, typedef_, while_ = symbols('if, typedef, while')
    expr = (if_ + typedef_) * while_
    routine = make_routine('test', expr)
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'double test(double if_, double typedef_, double while_) {\n'
        '   double test_result;\n'
        '   test_result = while_*(if_ + typedef_);\n'
        '   return test_result;\n'
        '}\n'
    )
    assert source == expected


def test_numbersymbol_c_code():
    routine = make_routine('test', pi**Catalan)
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'double test() {\n'
        '   double test_result;\n'
        '   double const Catalan = 0.915965594177219;\n'
        '   test_result = pow(M_PI, Catalan);\n'
        '   return test_result;\n'
        '}\n'
    )
    assert source == expected


def test_c_code_argument_order():
    expr = x + y
    routine = make_routine('test', expr, argument_sequence=[z, x, y])
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'double test(double z, double x, double y) {\n'
        '   double test_result;\n'
        '   test_result = x + y;\n'
        '   return test_result;\n'
        '}\n'
    )
    assert source == expected

    p = MatrixSymbol('p', 3, 1)
    routine = make_routine('test', expr, argument_sequence=[p, x, y])
    source = get_string(code_gen.dump_c, [routine])
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'double test(double *p, double x, double y) {\n'
        '   double test_result;\n'
        '   test_result = x + y;\n'
        '   return test_result;\n'
        '}\n'
    )
    assert source == expected


def test_simple_c_header():
    expr = (x + y)*z
    routine = make_routine('test', expr)
    code_gen = CCodeGen()
    source = get_string(code_gen.dump_h, [routine])
    expected = (
        '#ifndef PROJECT__FILE__H\n'
        '#define PROJECT__FILE__H\n'
        'double test(double x, double y, double z);\n'
        '#endif\n'
    )
    assert source == expected


def test_simple_c_codegen():
    expr = (x + y)*z
    result = codegen(('test', expr), 'C', 'file', header=False, empty=False)
    expected = [
        ('file.c',
         '#include "file.h"\n'
         '#include <math.h>\n'
         'double test(double x, double y, double z) {\n'
         '   double test_result;\n'
         '   test_result = z*(x + y);\n'
         '   return test_result;\n'
         '}\n'),
        ('file.h',
         '#ifndef PROJECT__FILE__H\n'
         '#define PROJECT__FILE__H\n'
         'double test(double x, double y, double z);\n'
         '#endif\n')
    ]
    assert result == expected


def test_multiple_results_c():
    expr1 = (x + y)*z
    expr2 = (x - y)*z
    routine = make_routine(
        'test',
        [expr1, expr2]
    )
    code_gen = CCodeGen()
    pytest.raises(CodeGenError, lambda: get_string(code_gen.dump_h, [routine]))


def test_no_results_c():
    pytest.raises(ValueError, lambda: make_routine('test', []))


def test_ansi_math1_codegen():
    # not included: log10
    name_expr = [
        ('test_fabs', abs(x)),
        ('test_acos', acos(x)),
        ('test_asin', asin(x)),
        ('test_atan', atan(x)),
        ('test_ceil', ceiling(x)),
        ('test_cos', cos(x)),
        ('test_cosh', cosh(x)),
        ('test_floor', floor(x)),
        ('test_log', log(x)),
        ('test_ln', ln(x)),
        ('test_sin', sin(x)),
        ('test_sinh', sinh(x)),
        ('test_sqrt', sqrt(x)),
        ('test_tan', tan(x)),
        ('test_tanh', tanh(x)),
    ]
    result = codegen(name_expr, 'C', 'file', header=False, empty=False)
    assert result[0][0] == 'file.c'
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test_fabs(double x) {\n   double test_fabs_result;\n   test_fabs_result = fabs(x);\n   return test_fabs_result;\n}\n'
        'double test_acos(double x) {\n   double test_acos_result;\n   test_acos_result = acos(x);\n   return test_acos_result;\n}\n'
        'double test_asin(double x) {\n   double test_asin_result;\n   test_asin_result = asin(x);\n   return test_asin_result;\n}\n'
        'double test_atan(double x) {\n   double test_atan_result;\n   test_atan_result = atan(x);\n   return test_atan_result;\n}\n'
        'double test_ceil(double x) {\n   double test_ceil_result;\n   test_ceil_result = ceil(x);\n   return test_ceil_result;\n}\n'
        'double test_cos(double x) {\n   double test_cos_result;\n   test_cos_result = cos(x);\n   return test_cos_result;\n}\n'
        'double test_cosh(double x) {\n   double test_cosh_result;\n   test_cosh_result = cosh(x);\n   return test_cosh_result;\n}\n'
        'double test_floor(double x) {\n   double test_floor_result;\n   test_floor_result = floor(x);\n   return test_floor_result;\n}\n'
        'double test_log(double x) {\n   double test_log_result;\n   test_log_result = log(x);\n   return test_log_result;\n}\n'
        'double test_ln(double x) {\n   double test_ln_result;\n   test_ln_result = log(x);\n   return test_ln_result;\n}\n'
        'double test_sin(double x) {\n   double test_sin_result;\n   test_sin_result = sin(x);\n   return test_sin_result;\n}\n'
        'double test_sinh(double x) {\n   double test_sinh_result;\n   test_sinh_result = sinh(x);\n   return test_sinh_result;\n}\n'
        'double test_sqrt(double x) {\n   double test_sqrt_result;\n   test_sqrt_result = sqrt(x);\n   return test_sqrt_result;\n}\n'
        'double test_tan(double x) {\n   double test_tan_result;\n   test_tan_result = tan(x);\n   return test_tan_result;\n}\n'
        'double test_tanh(double x) {\n   double test_tanh_result;\n   test_tanh_result = tanh(x);\n   return test_tanh_result;\n}\n'
    )
    assert result[1][0] == 'file.h'
    assert result[1][1] == (
        '#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n'
        'double test_fabs(double x);\ndouble test_acos(double x);\n'
        'double test_asin(double x);\ndouble test_atan(double x);\n'
        'double test_ceil(double x);\ndouble test_cos(double x);\n'
        'double test_cosh(double x);\ndouble test_floor(double x);\n'
        'double test_log(double x);\ndouble test_ln(double x);\n'
        'double test_sin(double x);\ndouble test_sinh(double x);\n'
        'double test_sqrt(double x);\ndouble test_tan(double x);\n'
        'double test_tanh(double x);\n#endif\n'
    )


def test_ansi_math2_codegen():
    # not included: frexp, ldexp, modf, fmod
    name_expr = [
        ('test_atan2', atan2(x, y)),
        ('test_pow', x**y),
    ]
    result = codegen(name_expr, 'C', 'file', header=False, empty=False)
    assert result[0][0] == 'file.c'
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test_atan2(double x, double y) {\n   double test_atan2_result;\n   test_atan2_result = atan2(x, y);\n   return test_atan2_result;\n}\n'
        'double test_pow(double x, double y) {\n   double test_pow_result;\n   test_pow_result = pow(x, y);\n   return test_pow_result;\n}\n'
    )
    assert result[1][0] == 'file.h'
    assert result[1][1] == (
        '#ifndef PROJECT__FILE__H\n#define PROJECT__FILE__H\n'
        'double test_atan2(double x, double y);\n'
        'double test_pow(double x, double y);\n'
        '#endif\n'
    )


def test_complicated_codegen():
    name_expr = [
        ('test1', ((sin(x) + cos(y) + tan(z))**7).expand()),
        ('test2', cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))),
    ]
    result = codegen(name_expr, 'C', 'file', header=False, empty=False)
    assert result[0][0] == 'file.c'
    assert result[0][1] == (
        '#include "file.h"\n#include <math.h>\n'
        'double test1(double x, double y, double z) {\n'
        '   double test1_result;\n'
        '   test1_result = '
        'pow(sin(x), 7) + '
        '7*pow(sin(x), 6)*cos(y) + '
        '7*pow(sin(x), 6)*tan(z) + '
        '21*pow(sin(x), 5)*pow(cos(y), 2) + '
        '42*pow(sin(x), 5)*cos(y)*tan(z) + '
        '21*pow(sin(x), 5)*pow(tan(z), 2) + '
        '35*pow(sin(x), 4)*pow(cos(y), 3) + '
        '105*pow(sin(x), 4)*pow(cos(y), 2)*tan(z) + '
        '105*pow(sin(x), 4)*cos(y)*pow(tan(z), 2) + '
        '35*pow(sin(x), 4)*pow(tan(z), 3) + '
        '35*pow(sin(x), 3)*pow(cos(y), 4) + '
        '140*pow(sin(x), 3)*pow(cos(y), 3)*tan(z) + '
        '210*pow(sin(x), 3)*pow(cos(y), 2)*pow(tan(z), 2) + '
        '140*pow(sin(x), 3)*cos(y)*pow(tan(z), 3) + '
        '35*pow(sin(x), 3)*pow(tan(z), 4) + '
        '21*pow(sin(x), 2)*pow(cos(y), 5) + '
        '105*pow(sin(x), 2)*pow(cos(y), 4)*tan(z) + '
        '210*pow(sin(x), 2)*pow(cos(y), 3)*pow(tan(z), 2) + '
        '210*pow(sin(x), 2)*pow(cos(y), 2)*pow(tan(z), 3) + '
        '105*pow(sin(x), 2)*cos(y)*pow(tan(z), 4) + '
        '21*pow(sin(x), 2)*pow(tan(z), 5) + '
        '7*sin(x)*pow(cos(y), 6) + '
        '42*sin(x)*pow(cos(y), 5)*tan(z) + '
        '105*sin(x)*pow(cos(y), 4)*pow(tan(z), 2) + '
        '140*sin(x)*pow(cos(y), 3)*pow(tan(z), 3) + '
        '105*sin(x)*pow(cos(y), 2)*pow(tan(z), 4) + '
        '42*sin(x)*cos(y)*pow(tan(z), 5) + '
        '7*sin(x)*pow(tan(z), 6) + '
        'pow(cos(y), 7) + '
        '7*pow(cos(y), 6)*tan(z) + '
        '21*pow(cos(y), 5)*pow(tan(z), 2) + '
        '35*pow(cos(y), 4)*pow(tan(z), 3) + '
        '35*pow(cos(y), 3)*pow(tan(z), 4) + '
        '21*pow(cos(y), 2)*pow(tan(z), 5) + '
        '7*cos(y)*pow(tan(z), 6) + '
        'pow(tan(z), 7);\n'
        '   return test1_result;\n'
        '}\n'
        'double test2(double x, double y, double z) {\n'
        '   double test2_result;\n'
        '   test2_result = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))));\n'
        '   return test2_result;\n'
        '}\n'
    )
    assert result[1][0] == 'file.h'
    assert result[1][1] == (
        '#ifndef PROJECT__FILE__H\n'
        '#define PROJECT__FILE__H\n'
        'double test1(double x, double y, double z);\n'
        'double test2(double x, double y, double z);\n'
        '#endif\n'
    )


def test_loops_c():
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y[i], A[i, j]*x[j])), 'C', 'file', header=False, empty=False)

    assert f1 == 'file.c'
    expected = (
        '#include "file.h"\n'
        '#include <math.h>\n'
        'void matrix_vector(double *A, int m, int n, double *x, double *y) {\n'
        '   for (int i=0; i<m; i++){\n'
        '      y[i] = 0;\n'
        '   }\n'
        '   for (int i=0; i<m; i++){\n'
        '      for (int j=0; j<n; j++){\n'
        '         y[i] = %(rhs)s + y[i];\n'
        '      }\n'
        '   }\n'
        '}\n'
    )

    assert code in (expected % {'rhs': f'A[{i * n + j}]*x[j]'},
                    expected % {'rhs': f'A[{j + i * n}]*x[j]'},
                    expected % {'rhs': f'x[j]*A[{i * n + j}]'},
                    expected % {'rhs': f'x[j]*A[{j + i * n}]'})
    assert f2 == 'file.h'
    assert interface == (
        '#ifndef PROJECT__FILE__H\n'
        '#define PROJECT__FILE__H\n'
        'void matrix_vector(double *A, int m, int n, double *x, double *y);\n'
        '#endif\n'
    )


def test_dummy_loops_c():
    i, m = symbols('i m', integer=True, cls=Dummy)
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx(i, m)
    ino = i.label.dummy_index
    mno = m.dummy_index
    expected = f"""#include "file.h"
#include <math.h>
void test_dummies(int m_{mno}, double *x, double *y) {{
   for (int i_{ino}=0; i_{ino}<m_{mno}; i_{ino}++){{
      y[i_{ino}] = x[i_{ino}];
   }}
}}
"""
    r = make_routine('test_dummies', Eq(y[i], x[i]))
    c = CCodeGen()
    code = get_string(c.dump_c, [r])
    assert code == expected


def test_partial_loops_c():
    # check that loop boundaries are determined by Idx, and array strides
    # determined by shape of IndexedBase object.
    n, m, o, p = symbols('n m o p', integer=True)
    A = IndexedBase('A', shape=(m, p))
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', (o, m - 5))  # Note: bounds are inclusive
    j = Idx('j', n)          # dimension n corresponds to bounds (0, n - 1)

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y[i], A[i, j]*x[j])), 'C', 'file', header=False, empty=False)

    assert f1 == 'file.c'

    upperi = m - 4
    rhs = f'x[j]*A[{j + i * p}]'
    expected = f"""#include "file.h"
#include <math.h>
void matrix_vector(double *A, int m, int n, int o, int p, double *x, double *y) {{
   for (int i=o; i<{upperi}; i++){{
      y[i] = 0;
   }}
   for (int i=o; i<{upperi}; i++){{
      for (int j=0; j<n; j++){{
         y[i] = {rhs} + y[i];
      }}
   }}
}}
"""

    assert code == expected
    assert f2 == 'file.h'
    assert interface == """#ifndef PROJECT__FILE__H
#define PROJECT__FILE__H
void matrix_vector(double *A, int m, int n, int o, int p, double *x, double *y);
#endif
"""


def test_output_arg_c():
    r = make_routine('foo', [Equality(y, sin(x)), cos(x)])
    c = CCodeGen()
    result = c.write([r], 'test', header=False, empty=False)
    assert result[0][0] == 'test.c'
    expected = (
        '#include "test.h"\n'
        '#include <math.h>\n'
        'double foo(double x, double *y) {\n'
        '   (*y) = sin(x);\n'
        '   double foo_result;\n'
        '   foo_result = cos(x);\n'
        '   return foo_result;\n'
        '}\n'
    )
    assert result[0][1] == expected


def test_output_arg_c_reserved_words():
    if_, while_ = symbols('if, while')
    r = make_routine('foo', [Equality(while_, sin(if_)), cos(if_)])
    c = CCodeGen()
    result = c.write([r], 'test', header=False, empty=False)
    assert result[0][0] == 'test.c'
    expected = (
        '#include "test.h"\n'
        '#include <math.h>\n'
        'double foo(double if_, double *while_) {\n'
        '   (*while_) = sin(if_);\n'
        '   double foo_result;\n'
        '   foo_result = cos(if_);\n'
        '   return foo_result;\n'
        '}\n'
    )
    assert result[0][1] == expected


def test_ccode_results_named_ordered():
    A = MatrixSymbol('A', 1, 3)
    expr1 = Equality(A, Matrix([[1, 2, x]]))
    expr2 = Equality(C, (x + y)*z)
    expr3 = Equality(B, 2*x)
    name_expr = ('test', [expr1, expr2, expr3])
    result = codegen(name_expr, 'c', 'test', header=False, empty=False,
                     argument_sequence=(x, C, z, y, A, B))
    source = result[0][1]
    expected = (
        '#include "test.h"\n'
        '#include <math.h>\n'
        'void test(double x, double *C, double z, double y, double *A, double *B) {\n'
        '   (*C) = z*(x + y);\n'
        '   A[0] = 1;\n'
        '   A[1] = 2;\n'
        '   A[2] = x;\n'
        '   (*B) = 2*x;\n'
        '}\n'
    )
    assert source == expected


def test_ccode_matrixsymbol_slice():
    A = MatrixSymbol('A', 5, 3)
    B = MatrixSymbol('B', 1, 3)
    C = MatrixSymbol('C', 1, 3)
    D = MatrixSymbol('D', 5, 1)
    name_expr = ('test', [Equality(B, A[0, :]),
                          Equality(C, A[1, :]),
                          Equality(D, A[:, 2])])
    result = codegen(name_expr, 'c', 'test', header=False, empty=False)
    source = result[0][1]
    expected = (
        '#include "test.h"\n'
        '#include <math.h>\n'
        'void test(double *A, double *B, double *C, double *D) {\n'
        '   B[0] = A[0];\n'
        '   B[1] = A[1];\n'
        '   B[2] = A[2];\n'
        '   C[0] = A[3];\n'
        '   C[1] = A[4];\n'
        '   C[2] = A[5];\n'
        '   D[0] = A[2];\n'
        '   D[1] = A[5];\n'
        '   D[2] = A[8];\n'
        '   D[3] = A[11];\n'
        '   D[4] = A[14];\n'
        '}\n'
    )
    assert source == expected


def test_empty_f_code():
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [])
    assert source == ''


def test_empty_f_code_with_header():
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [], header=True)
    assert source[:82] == (
        '!******************************************************************************\n!*'
    )
    #   "                    Code generated with diofant x.y.z                    "
    assert source[158:] == (                                                              '*\n'
                                                                                          '!*                                                                            *\n'
                                                                                          '!*         See https://diofant.readthedocs.io/ for more information.          *\n'
                                                                                          '!*                                                                            *\n'
                                                                                          "!*                       This file is part of 'project'                       *\n"
                                                                                          '!******************************************************************************\n')


def test_empty_f_header():
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_h, [])
    assert source == ''


def test_simple_f_code():
    expr = (x + y)*z
    routine = make_routine('test', expr)
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        'REAL*8 function test(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'test = z*(x + y)\n'
        'end function\n'
    )
    assert source == expected


def test_numbersymbol_f_code():
    routine = make_routine('test', pi**Catalan)
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        'REAL*8 function test()\n'
        'implicit none\n'
        'REAL*8, parameter :: Catalan = 0.915965594177219d0\n'
        'REAL*8, parameter :: pi = 3.14159265358979d0\n'
        'test = pi**Catalan\n'
        'end function\n'
    )
    assert source == expected


def test_erf_f_code():
    routine = make_routine('test', erf(x) - erf(-2 * x))
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        'REAL*8 function test(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test = erf(x) + erf(2.0d0*x)\n'
        'end function\n'
    )
    assert source == expected, source


def test_not_fortran_f_code():
    routine = make_routine('test', besseli(1, x))
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        'REAL*8 function test(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8 :: besseli\n'
        'test = besseli(1.0, x)\n'
        'end function\n'
    )
    assert source == expected, source

    f = Function('f')
    routine = make_routine('test', f(x).diff(x))
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        'REAL*8 function test(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8 :: Derivative(f(x), x)\n'
        'test = Derivative(f(x), x)\n'
        'end function\n'
    )
    assert source == expected, source


def test_f_code_argument_order():
    expr = x + y
    routine = make_routine('test', expr, argument_sequence=[z, x, y])
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = (
        'REAL*8 function test(z, x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: z\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'test = x + y\n'
        'end function\n'
    )
    assert source == expected


def test_simple_f_header():
    expr = (x + y)*z
    routine = make_routine('test', expr)
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_h, [routine])
    expected = (
        'interface\n'
        'REAL*8 function test(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'end function\n'
        'end interface\n'
    )
    assert source == expected


def test_simple_f_codegen():
    expr = (x + y)*z
    result = codegen(
        ('test', expr), 'F95', 'file', header=False, empty=False)
    expected = [
        ('file.f90',
         'REAL*8 function test(x, y, z)\n'
         'implicit none\n'
         'REAL*8, intent(in) :: x\n'
         'REAL*8, intent(in) :: y\n'
         'REAL*8, intent(in) :: z\n'
         'test = z*(x + y)\n'
         'end function\n'),
        ('file.h',
         'interface\n'
         'REAL*8 function test(x, y, z)\n'
         'implicit none\n'
         'REAL*8, intent(in) :: x\n'
         'REAL*8, intent(in) :: y\n'
         'REAL*8, intent(in) :: z\n'
         'end function\n'
         'end interface\n')
    ]
    assert result == expected


def test_multiple_results_f():
    expr1 = (x + y)*z
    expr2 = (x - y)*z
    routine = make_routine(
        'test',
        [expr1, expr2]
    )
    code_gen = FCodeGen()
    pytest.raises(CodeGenError, lambda: get_string(code_gen.dump_h, [routine]))


def test_no_results_f():
    pytest.raises(ValueError, lambda: make_routine('test', []))


def test_intrinsic_math_codegen():
    # not included: log10
    name_expr = [
        ('test_abs', abs(x)),
        ('test_acos', acos(x)),
        ('test_asin', asin(x)),
        ('test_atan', atan(x)),
        ('test_cos', cos(x)),
        ('test_cosh', cosh(x)),
        ('test_log', log(x)),
        ('test_ln', ln(x)),
        ('test_sin', sin(x)),
        ('test_sinh', sinh(x)),
        ('test_sqrt', sqrt(x)),
        ('test_tan', tan(x)),
        ('test_tanh', tanh(x)),
    ]
    result = codegen(name_expr, 'F95', 'file', header=False, empty=False)
    assert result[0][0] == 'file.f90'
    expected = (
        'REAL*8 function test_abs(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_abs = Abs(x)\n'
        'end function\n'
        'REAL*8 function test_acos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_acos = acos(x)\n'
        'end function\n'
        'REAL*8 function test_asin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_asin = asin(x)\n'
        'end function\n'
        'REAL*8 function test_atan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_atan = atan(x)\n'
        'end function\n'
        'REAL*8 function test_cos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_cos = cos(x)\n'
        'end function\n'
        'REAL*8 function test_cosh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_cosh = cosh(x)\n'
        'end function\n'
        'REAL*8 function test_log(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_log = log(x)\n'
        'end function\n'
        'REAL*8 function test_ln(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_ln = log(x)\n'
        'end function\n'
        'REAL*8 function test_sin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_sin = sin(x)\n'
        'end function\n'
        'REAL*8 function test_sinh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_sinh = sinh(x)\n'
        'end function\n'
        'REAL*8 function test_sqrt(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_sqrt = sqrt(x)\n'
        'end function\n'
        'REAL*8 function test_tan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_tan = tan(x)\n'
        'end function\n'
        'REAL*8 function test_tanh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'test_tanh = tanh(x)\n'
        'end function\n'
    )
    assert result[0][1] == expected

    assert result[1][0] == 'file.h'
    expected = (
        'interface\n'
        'REAL*8 function test_abs(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_acos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_asin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_atan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_cos(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_cosh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_log(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_ln(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_sin(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_sinh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_sqrt(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_tan(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_tanh(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'end function\n'
        'end interface\n'
    )
    assert result[1][1] == expected


def test_intrinsic_math2_codegen():
    # not included: frexp, ldexp, modf, fmod
    name_expr = [
        ('test_atan2', atan2(x, y)),
        ('test_pow', x**y),
    ]
    result = codegen(name_expr, 'F95', 'file', header=False, empty=False)
    assert result[0][0] == 'file.f90'
    expected = (
        'REAL*8 function test_atan2(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'test_atan2 = atan2(x, y)\n'
        'end function\n'
        'REAL*8 function test_pow(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'test_pow = x**y\n'
        'end function\n'
    )
    assert result[0][1] == expected

    assert result[1][0] == 'file.h'
    expected = (
        'interface\n'
        'REAL*8 function test_atan2(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test_pow(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'end function\n'
        'end interface\n'
    )
    assert result[1][1] == expected


def test_complicated_codegen_f95():
    name_expr = [
        ('test1', ((sin(x) + cos(y) + tan(z))**7).expand()),
        ('test2', cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))),
    ]
    result = codegen(name_expr, 'F95', 'file', header=False, empty=False)
    assert result[0][0] == 'file.f90'
    expected = (
        'REAL*8 function test1(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'test1 = sin(x)**7 + 7*sin(x)**6*cos(y) + 7*sin(x)**6*tan(z) + 21*sin(x) &\n'
        '      **5*cos(y)**2 + 42*sin(x)**5*cos(y)*tan(z) + 21*sin(x)**5*tan(z) &\n'
        '      **2 + 35*sin(x)**4*cos(y)**3 + 105*sin(x)**4*cos(y)**2*tan(z) + &\n'
        '      105*sin(x)**4*cos(y)*tan(z)**2 + 35*sin(x)**4*tan(z)**3 + 35*sin( &\n'
        '      x)**3*cos(y)**4 + 140*sin(x)**3*cos(y)**3*tan(z) + 210*sin(x)**3* &\n'
        '      cos(y)**2*tan(z)**2 + 140*sin(x)**3*cos(y)*tan(z)**3 + 35*sin(x) &\n'
        '      **3*tan(z)**4 + 21*sin(x)**2*cos(y)**5 + 105*sin(x)**2*cos(y)**4* &\n'
        '      tan(z) + 210*sin(x)**2*cos(y)**3*tan(z)**2 + 210*sin(x)**2*cos(y) &\n'
        '      **2*tan(z)**3 + 105*sin(x)**2*cos(y)*tan(z)**4 + 21*sin(x)**2*tan &\n'
        '      (z)**5 + 7*sin(x)*cos(y)**6 + 42*sin(x)*cos(y)**5*tan(z) + 105* &\n'
        '      sin(x)*cos(y)**4*tan(z)**2 + 140*sin(x)*cos(y)**3*tan(z)**3 + 105 &\n'
        '      *sin(x)*cos(y)**2*tan(z)**4 + 42*sin(x)*cos(y)*tan(z)**5 + 7*sin( &\n'
        '      x)*tan(z)**6 + cos(y)**7 + 7*cos(y)**6*tan(z) + 21*cos(y)**5*tan( &\n'
        '      z)**2 + 35*cos(y)**4*tan(z)**3 + 35*cos(y)**3*tan(z)**4 + 21*cos( &\n'
        '      y)**2*tan(z)**5 + 7*cos(y)*tan(z)**6 + tan(z)**7\n'
        'end function\n'
        'REAL*8 function test2(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'test2 = cos(cos(cos(cos(cos(cos(cos(cos(x + y + z))))))))\n'
        'end function\n'
    )
    assert result[0][1] == expected
    assert result[1][0] == 'file.h'
    expected = (
        'interface\n'
        'REAL*8 function test1(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'end function\n'
        'end interface\n'
        'interface\n'
        'REAL*8 function test2(x, y, z)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'end function\n'
        'end interface\n'
    )
    assert result[1][1] == expected


def test_loops():
    n, m = symbols('n,m', integer=True)
    A, x, y = map(IndexedBase, 'Axy')
    i = Idx('i', m)
    j = Idx('j', n)

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y[i], A[i, j]*x[j])), 'F95', 'file', header=False, empty=False)

    assert f1 == 'file.f90'
    expected = (
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(out), dimension(1:m) :: y\n'
        'INTEGER*4 :: i\n'
        'INTEGER*4 :: j\n'
        'do i = 1, m\n'
        '   y(i) = 0\n'
        'end do\n'
        'do i = 1, m\n'
        '   do j = 1, n\n'
        '      y(i) = %(rhs)s + y(i)\n'
        '   end do\n'
        'end do\n'
        'end subroutine\n'
    )

    assert code in (expected % {'rhs': 'A(i, j)*x(j)'},
                    expected % {'rhs': 'x(j)*A(i, j)'})
    assert f2 == 'file.h'
    assert interface == (
        'interface\n'
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(out), dimension(1:m) :: y\n'
        'end subroutine\n'
        'end interface\n'
    )


def test_dummy_loops_f95():
    i, m = symbols('i m', integer=True, cls=Dummy)
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx(i, m)
    icount = i.label.dummy_index
    mcount = m.dummy_index
    expected = f"""subroutine test_dummies(m_{mcount}, x, y)
implicit none
INTEGER*4, intent(in) :: m_{mcount}
REAL*8, intent(in), dimension(1:m_{mcount}) :: x
REAL*8, intent(out), dimension(1:m_{mcount}) :: y
INTEGER*4 :: i_{icount}
do i_{icount} = 1, m_{mcount}
   y(i_{icount}) = x(i_{icount})
end do
end subroutine
"""
    r = make_routine('test_dummies', Eq(y[i], x[i]))
    c = FCodeGen()
    code = get_string(c.dump_f95, [r])
    assert code == expected


def test_loops_InOut():
    i, j, n, m = symbols('i,j,n,m', integer=True)
    A = IndexedBase('A')[Idx(i, m), Idx(j, n)]
    x = IndexedBase('x')[Idx(j, n)]
    y = IndexedBase('y')[Idx(i, m)]

    (f1, code), (f2, interface) = codegen(
        ('matrix_vector', Eq(y, y + A*x)), 'F95', 'file', header=False, empty=False)

    assert f1 == 'file.f90'
    expected = (
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(inout), dimension(1:m) :: y\n'
        'INTEGER*4 :: i\n'
        'INTEGER*4 :: j\n'
        'do i = 1, m\n'
        '   do j = 1, n\n'
        '      y(i) = %(rhs)s + y(i)\n'
        '   end do\n'
        'end do\n'
        'end subroutine\n'
    )

    assert code in (expected % {'rhs': 'A(i, j)*x(j)'},
                    expected % {'rhs': 'x(j)*A(i, j)'})
    assert f2 == 'file.h'
    assert interface == (
        'interface\n'
        'subroutine matrix_vector(A, m, n, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'INTEGER*4, intent(in) :: n\n'
        'REAL*8, intent(in), dimension(1:m, 1:n) :: A\n'
        'REAL*8, intent(in), dimension(1:n) :: x\n'
        'REAL*8, intent(inout), dimension(1:m) :: y\n'
        'end subroutine\n'
        'end interface\n'
    )


def test_partial_loops_f():
    # check that loop boundaries are determined by Idx, and array strides
    # determined by shape of IndexedBase object.
    n, m, o, p = symbols('n m o p', integer=True)
    A = IndexedBase('A', shape=(m, p))
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', (o, m - 5))  # Note: bounds are inclusive
    j = Idx('j', n)          # dimension n corresponds to bounds (0, n - 1)

    (_, code), _ = codegen(('matrix_vector', Eq(y[i], A[i, j]*x[j])),
                           'F95', 'file', header=False, empty=False)

    rhs = 'x(j)*A(i, j)'
    iup = str(m - 4)
    ilow = str(1 + o)
    iup_ilow = str(m - 4 - o)
    expected = f"""subroutine matrix_vector(A, m, n, o, p, x, y)
implicit none
INTEGER*4, intent(in) :: m
INTEGER*4, intent(in) :: n
INTEGER*4, intent(in) :: o
INTEGER*4, intent(in) :: p
REAL*8, intent(in), dimension(1:m, 1:p) :: A
REAL*8, intent(in), dimension(1:n) :: x
REAL*8, intent(out), dimension(1:{iup_ilow}) :: y
INTEGER*4 :: i
INTEGER*4 :: j
do i = {ilow}, {iup}
   y(i) = 0
end do
do i = {ilow}, {iup}
   do j = 1, n
      y(i) = {rhs} + y(i)
   end do
end do
end subroutine
"""

    assert code == expected


def test_output_arg_f():
    r = make_routine('foo', [Equality(y, sin(x)), cos(x)])
    c = FCodeGen()
    result = c.write([r], 'test', header=False, empty=False)
    assert result[0][0] == 'test.f90'
    assert result[0][1] == (
        'REAL*8 function foo(x, y)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(out) :: y\n'
        'y = sin(x)\n'
        'foo = cos(x)\n'
        'end function\n'
    )


def test_inline_function():
    n, m = symbols('n m', integer=True)
    x, y = map(IndexedBase, 'xy')
    i = Idx('i', m)
    p = FCodeGen()
    func = implemented_function('func', Lambda(n, n*(n + 1)))
    routine = make_routine('test_inline', Eq(y[i], func(x[i])))
    code = get_string(p.dump_f95, [routine])
    expected = (
        'subroutine test_inline(m, x, y)\n'
        'implicit none\n'
        'INTEGER*4, intent(in) :: m\n'
        'REAL*8, intent(in), dimension(1:m) :: x\n'
        'REAL*8, intent(out), dimension(1:m) :: y\n'
        'INTEGER*4 :: i\n'
        'do i = 1, m\n'
        '   y(i) = %s*%s\n'
        'end do\n'
        'end subroutine\n'
    )
    args = ('x(i)', '(x(i) + 1)')
    assert code in (expected % args, expected % args[::-1])


def test_f_code_call_signature_wrap():
    # Issue sympy/sympy#7934
    x = symbols('x:20')
    expr = 0
    for sym in x:
        expr += sym
    routine = make_routine('test', expr)
    code_gen = FCodeGen()
    source = get_string(code_gen.dump_f95, [routine])
    expected = """\
REAL*8 function test(x0, x1, x10, x11, x12, x13, x14, x15, x16, x17, x18, &
      x19, x2, x3, x4, x5, x6, x7, x8, x9)
implicit none
REAL*8, intent(in) :: x0
REAL*8, intent(in) :: x1
REAL*8, intent(in) :: x10
REAL*8, intent(in) :: x11
REAL*8, intent(in) :: x12
REAL*8, intent(in) :: x13
REAL*8, intent(in) :: x14
REAL*8, intent(in) :: x15
REAL*8, intent(in) :: x16
REAL*8, intent(in) :: x17
REAL*8, intent(in) :: x18
REAL*8, intent(in) :: x19
REAL*8, intent(in) :: x2
REAL*8, intent(in) :: x3
REAL*8, intent(in) :: x4
REAL*8, intent(in) :: x5
REAL*8, intent(in) :: x6
REAL*8, intent(in) :: x7
REAL*8, intent(in) :: x8
REAL*8, intent(in) :: x9
test = x0 + x1 + x10 + x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 + &
      x19 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
end function
"""
    assert source == expected


def test_check_case():
    pytest.raises(CodeGenError, lambda: codegen(('test', x*X), 'f95', 'prefix'))


def test_check_case_false_positive():
    # The upper case/lower case exception should not be triggered by Diofant
    # objects that differ only because of assumptions.  (It may be useful to
    # have a check for that as well, but here we only want to test against
    # false positives with respect to case checking.)
    x2 = symbols('x', my_assumption=True)
    try:
        codegen(('test', x*x2), 'f95', 'prefix')
    except CodeGenError as exc:
        if exc.args[0].startswith('Fortran ignores case.'):
            raise AssertionError('This exception should not '
                                 'be raised!') from exc


def test_c_fortran_omit_routine_name():
    name_expr = [('foo', 2*x)]
    result = codegen(name_expr, 'F95', header=False, empty=False)
    expresult = codegen(name_expr, 'F95', 'foo', header=False, empty=False)
    assert result[0][1] == expresult[0][1]

    name_expr = ('foo', x*y)
    result = codegen(name_expr, 'F95', header=False, empty=False)
    expresult = codegen(name_expr, 'F95', 'foo', header=False, empty=False)
    assert result[0][1] == expresult[0][1]

    name_expr = ('foo', Matrix([[x, y], [x+y, x-y]]))
    result = codegen(name_expr, 'C', header=False, empty=False)
    expresult = codegen(name_expr, 'C', 'foo', header=False, empty=False)
    assert result[0][1] == expresult[0][1]


def test_fcode_matrix_output():
    e1 = x + y
    e2 = Matrix([[x, y], [z, 16]])
    name_expr = ('test', (e1, e2))
    result = codegen(name_expr, 'f95', 'test', header=False, empty=False)
    source = result[0][1]
    expected = (
        'REAL*8 function test(x, y, z, out_%(hash)s)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(in) :: z\n'
        'REAL*8, intent(out), dimension(1:2, 1:2) :: out_%(hash)s\n'
        'out_%(hash)s(1, 1) = x\n'
        'out_%(hash)s(2, 1) = z\n'
        'out_%(hash)s(1, 2) = y\n'
        'out_%(hash)s(2, 2) = 16\n'
        'test = x + y\n'
        'end function\n'
    )
    # look for the magic number
    a = source.splitlines()[5]
    b = a.split('_')
    out = b[1]
    expected = expected % {'hash': out}
    assert source == expected


def test_fcode_results_named_ordered():
    A = MatrixSymbol('A', 1, 3)
    expr1 = Equality(A, Matrix([[1, 2, x]]))
    expr2 = Equality(C, (x + y)*z)
    expr3 = Equality(B, 2*x)
    name_expr = ('test', [expr1, expr2, expr3])
    result = codegen(name_expr, 'f95', 'test', header=False, empty=False,
                     argument_sequence=(x, z, y, C, A, B))
    source = result[0][1]
    expected = (
        'subroutine test(x, z, y, C, A, B)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'REAL*8, intent(in) :: z\n'
        'REAL*8, intent(in) :: y\n'
        'REAL*8, intent(out) :: C\n'
        'REAL*8, intent(out) :: B\n'
        'REAL*8, intent(out), dimension(1:1, 1:3) :: A\n'
        'C = z*(x + y)\n'
        'A(1, 1) = 1\n'
        'A(1, 2) = 2\n'
        'A(1, 3) = x\n'
        'B = 2*x\n'
        'end subroutine\n'
    )
    assert source == expected


def test_fcode_matrixsymbol_slice():
    A = MatrixSymbol('A', 2, 3)
    B = MatrixSymbol('B', 1, 3)
    C = MatrixSymbol('C', 1, 3)
    D = MatrixSymbol('D', 2, 1)
    name_expr = ('test', [Equality(B, A[0, :]),
                          Equality(C, A[1, :]),
                          Equality(D, A[:, 2])])
    result = codegen(name_expr, 'f95', 'test', header=False, empty=False)
    source = result[0][1]
    expected = (
        'subroutine test(A, B, C, D)\n'
        'implicit none\n'
        'REAL*8, intent(in), dimension(1:2, 1:3) :: A\n'
        'REAL*8, intent(out), dimension(1:1, 1:3) :: B\n'
        'REAL*8, intent(out), dimension(1:1, 1:3) :: C\n'
        'REAL*8, intent(out), dimension(1:2, 1:1) :: D\n'
        'B(1, 1) = A(1, 1)\n'
        'B(1, 2) = A(1, 2)\n'
        'B(1, 3) = A(1, 3)\n'
        'C(1, 1) = A(2, 1)\n'
        'C(1, 2) = A(2, 2)\n'
        'C(1, 3) = A(2, 3)\n'
        'D(1, 1) = A(1, 3)\n'
        'D(2, 1) = A(2, 3)\n'
        'end subroutine\n'
    )
    assert source == expected


def test_fcode_matrixsymbol_slice_autoname():
    # see issue sympy/sympy#8093
    A = MatrixSymbol('A', 2, 3)
    name_expr = ('test', A[:, 1])
    result = codegen(name_expr, 'f95', 'test', header=False, empty=False)
    source = result[0][1]
    expected = (
        'subroutine test(A, out_%(hash)s)\n'
        'implicit none\n'
        'REAL*8, intent(in), dimension(1:2, 1:3) :: A\n'
        'REAL*8, intent(out), dimension(1:2, 1:1) :: out_%(hash)s\n'
        'out_%(hash)s(1, 1) = A(1, 2)\n'
        'out_%(hash)s(2, 1) = A(2, 2)\n'
        'end subroutine\n'
    )
    # look for the magic number
    a = source.splitlines()[3]
    b = a.split('_')
    out = b[1]
    expected = expected % {'hash': out}
    assert source == expected


def test_global_vars():
    result = codegen(('f', x*y), 'F95', header=False, empty=False,
                     global_vars=(y,))
    source = result[0][1]
    expected = (
        'REAL*8 function f(x)\n'
        'implicit none\n'
        'REAL*8, intent(in) :: x\n'
        'f = x*y\n'
        'end function\n'
    )
    assert source == expected

    result = codegen(('f', x*y+z), 'C', header=False, empty=False,
                     global_vars=(z, t))
    source = result[0][1]
    expected = (
        '#include "f.h"\n'
        '#include <math.h>\n'
        'double f(double x, double y) {\n'
        '   double f_result;\n'
        '   f_result = x*y + z;\n'
        '   return f_result;\n'
        '}\n'
    )
    assert source == expected


def test_ccode_cse():
    cg = CCodeGen(cse=True)
    e = MatrixSymbol('e', 3, 1)

    pytest.raises(ValueError, lambda: cg.routine('test', [], None))
    pytest.raises(CodeGenError, lambda: cg.routine('test', [e], None))

    routines = [cg.routine('test', [Equality(e, Matrix([[a*b], [a*b + c*d], [a*b*c*d]]))], None)]
    result = cg.write(routines, prefix='test', to_files=False, header=False, empty=False)
    source = result[0][1]
    expected = (
        '#include "test.h"\n'
        '#include <math.h>\n'
        'void test(double a, double b, double c, double d, double *e) {\n'
        '   const double x0 = a*b;\n'
        '   const double x1 = c*d;\n'
        '   e[0] = x0;\n'
        '   e[1] = x0 + x1;\n'
        '   e[2] = x0*x1;\n'
        '}\n'
    )
    assert source == expected

    routines = [cg.routine('test', Equality(e, Matrix([[a*b], [a*b + c*d], [a*b*c*d]])), None)]
    result = cg.write(routines, prefix='test', to_files=False, header=False, empty=False)
    source = result[0][1]
    assert source == expected

    routines = [cg.routine('test', Matrix([[a*b], [a*b + c*d], [a*b*c*d]]), None)]
    result = cg.write(routines, prefix='test', to_files=False, header=False, empty=False)
    source = result[0][1]
    expected = (
        '#include "test.h"\n'
        '#include <math.h>\n'
        'void test(double a, double b, double c, double d, double *out_%(hash)s) {\n'
        '   const double x0 = a*b;\n'
        '   const double x1 = c*d;\n'
        '   out_%(hash)s[0] = x0;\n'
        '   out_%(hash)s[1] = x0 + x1;\n'
        '   out_%(hash)s[2] = x0*x1;\n'
        '}\n'
    )
    # look for the magic number
    out = source.splitlines()[5].split('_')[1].split('[')[0]
    expected = expected % {'hash': out}
    assert source == expected
