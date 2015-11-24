from sympy.core import pi, oo, symbols, Rational, Integer, GoldenRatio, EulerGamma, Catalan, Lambda, Dummy
from sympy.functions import Piecewise, sin, cos, Abs, exp, ceiling, sqrt
from sympy.utilities.pytest import raises
from sympy.printing.jscode import JavascriptCodePrinter
from sympy.utilities.lambdify import implemented_function
from sympy.tensor import IndexedBase, Idx
from sympy.matrices import Matrix, MatrixSymbol

from sympy import jscode

x, y, z = symbols('x,y,z')


def test_printmethod():
    assert jscode(Abs(x)) == "Math.abs(x)"


def test_jscode_sqrt():
    assert jscode(sqrt(x)) == "Math.sqrt(x)"
    assert jscode(x**0.5) == "Math.sqrt(x)"
    assert jscode(sqrt(x)) == "Math.sqrt(x)"


def test_jscode_Pow():
    g = implemented_function('g', Lambda(x, 2*x))
    assert jscode(x**3) == "Math.pow(x, 3)"
    assert jscode(x**(y**3)) == "Math.pow(x, Math.pow(y, 3))"
    assert jscode(1/(g(x)*3.5)**(x - y**x)/(x**2 + y)) == \
        "Math.pow(3.5*2*x, -x + Math.pow(y, x))/(Math.pow(x, 2) + y)"
    assert jscode(x**-1.0) == '1/x'


def test_jscode_constants_mathh():
    assert jscode(exp(1)) == "Math.E"
    assert jscode(pi) == "Math.PI"
    assert jscode(oo) == "Number.POSITIVE_INFINITY"
    assert jscode(-oo) == "Number.NEGATIVE_INFINITY"


def test_jscode_constants_other():
    assert jscode(
        2*GoldenRatio) == "var GoldenRatio = 1.61803398874989;\n2*GoldenRatio"
    assert jscode(2*Catalan) == "var Catalan = 0.915965594177219;\n2*Catalan"
    assert jscode(
        2*EulerGamma) == "var EulerGamma = 0.577215664901533;\n2*EulerGamma"


def test_jscode_Rational():
    assert jscode(Rational(3, 7)) == "3/7"
    assert jscode(Rational(18, 9)) == "2"
    assert jscode(Rational(3, -7)) == "-3/7"
    assert jscode(Rational(-3, -7)) == "3/7"


def test_jscode_Integer():
    assert jscode(Integer(67)) == "67"
    assert jscode(Integer(-1)) == "-1"


def test_jscode_functions():
    assert jscode(sin(x) ** cos(x)) == "Math.pow(Math.sin(x), Math.cos(x))"


def test_jscode_inline_function():
    x = symbols('x')
    g = implemented_function('g', Lambda(x, 2*x))
    assert jscode(g(x)) == "2*x"
    g = implemented_function('g', Lambda(x, 2*x/Catalan))
    assert jscode(g(x)) == "var Catalan = %s;\n2*x/Catalan" % Catalan.n()
    A = IndexedBase('A')
    i = Idx('i', symbols('n', integer=True))
    g = implemented_function('g', Lambda(x, x*(1 + x)*(2 + x)))
    assert jscode(g(A[i]), assign_to=A[i]) == (
        "for (var i=0; i<n; i++){\n"
        "   A[i] = (A[i] + 1)*(A[i] + 2)*A[i];\n"
        "}"
    )


def test_jscode_exceptions():
    assert jscode(ceiling(x)) == "Math.ceil(x)"
    assert jscode(Abs(x)) == "Math.abs(x)"


def test_jscode_boolean():
    assert jscode(x & y) == "x && y"
    assert jscode(x | y) == "x || y"
    assert jscode(~x) == "!x"
    assert jscode(x & y & z) == "x && y && z"
    assert jscode(x | y | z) == "x || y || z"
    assert jscode((x & y) | z) == "z || x && y"
    assert jscode((x | y) & z) == "z && (x || y)"


def test_jscode_Piecewise():
    expr = Piecewise((x, x < 1), (x**2, True))
    p = jscode(expr)
    s = \
"""\
((x < 1) ? (
   x
)
: (
   Math.pow(x, 2)
))\
"""
    assert p == s
    assert jscode(expr, assign_to="c") == (
    "if (x < 1) {\n"
    "   c = x;\n"
    "}\n"
    "else {\n"
    "   c = Math.pow(x, 2);\n"
    "}")
    # Check that Piecewise without a True (default) condition error
    expr = Piecewise((x, x < 1), (x**2, x > 1), (sin(x), x > 0))
    raises(ValueError, lambda: jscode(expr))


def test_jscode_Piecewise_deep():
    p = jscode(2*Piecewise((x, x < 1), (x**2, True)))
    s = \
"""\
2*((x < 1) ? (
   x
)
: (
   Math.pow(x, 2)
))\
"""
    assert p == s


def test_jscode_settings():
    raises(TypeError, lambda: jscode(sin(x), method="garbage"))


def test_jscode_Indexed():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o = symbols('n m o', integer=True)
    i, j, k = Idx('i', n), Idx('j', m), Idx('k', o)
    p = JavascriptCodePrinter()
    p._not_c = set()

    x = IndexedBase('x')[j]
    assert p._print_Indexed(x) == 'x[j]'
    A = IndexedBase('A')[i, j]
    assert p._print_Indexed(A) == 'A[%s]' % (m*i+j)
    B = IndexedBase('B')[i, j, k]
    assert p._print_Indexed(B) == 'B[%s]' % (i*o*m+j*o+k)

    assert p._not_c == set()


def test_jscode_loops_matrix_vector():
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)

    s = (
        'for (var i=0; i<m; i++){\n'
        '   y[i] = 0;\n'
        '}\n'
        'for (var i=0; i<m; i++){\n'
        '   for (var j=0; j<n; j++){\n'
        '      y[i] = x[j]*A[n*i + j] + y[i];\n'
        '   }\n'
        '}'
    )
    c = jscode(A[i, j]*x[j], assign_to=y[i])
    assert c == s


def test_dummy_loops():
    i, m = symbols('i m', integer=True, cls=Dummy)
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx(i, m)

    expected = (
        'for (var i_%(icount)i=0; i_%(icount)i<m_%(mcount)i; i_%(icount)i++){\n'
        '   y[i_%(icount)i] = x[i_%(icount)i];\n'
        '}'
    ) % {'icount': i.label.dummy_index, 'mcount': m.dummy_index}
    code = jscode(x[i], assign_to=y[i])
    assert code == expected


def test_jscode_loops_add():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m = symbols('n m', integer=True)
    A = IndexedBase('A')
    x = IndexedBase('x')
    y = IndexedBase('y')
    z = IndexedBase('z')
    i = Idx('i', m)
    j = Idx('j', n)

    s = (
        'for (var i=0; i<m; i++){\n'
        '   y[i] = x[i] + z[i];\n'
        '}\n'
        'for (var i=0; i<m; i++){\n'
        '   for (var j=0; j<n; j++){\n'
        '      y[i] = x[j]*A[n*i + j] + y[i];\n'
        '   }\n'
        '}'
    )
    c = jscode(A[i, j]*x[j] + x[i] + z[i], assign_to=y[i])
    assert c == s


def test_jscode_loops_multiple_contractions():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    a = IndexedBase('a')
    b = IndexedBase('b')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)
    k = Idx('k', o)
    l = Idx('l', p)

    s = (
        'for (var i=0; i<m; i++){\n'
        '   y[i] = 0;\n'
        '}\n'
        'for (var i=0; i<m; i++){\n'
        '   for (var j=0; j<n; j++){\n'
        '      for (var k=0; k<o; k++){\n'
        '         for (var l=0; l<p; l++){\n'
        '            y[i] = y[i] + b[%s]*a[%s];\n' % (j*o*p + k*p + l, i*n*o*p + j*o*p + k*p + l) +
        '         }\n'
        '      }\n'
        '   }\n'
        '}'
    )
    c = jscode(b[j, k, l]*a[i, j, k, l], assign_to=y[i])
    assert c == s


def test_jscode_loops_addfactor():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    a = IndexedBase('a')
    b = IndexedBase('b')
    c = IndexedBase('c')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)
    k = Idx('k', o)
    l = Idx('l', p)

    s = (
        'for (var i=0; i<m; i++){\n'
        '   y[i] = 0;\n'
        '}\n'
        'for (var i=0; i<m; i++){\n'
        '   for (var j=0; j<n; j++){\n'
        '      for (var k=0; k<o; k++){\n'
        '         for (var l=0; l<p; l++){\n'
        '            y[i] = (a[%s] + b[%s])*c[%s] + y[i];\n' % (i*n*o*p + j*o*p + k*p + l, i*n*o*p + j*o*p + k*p + l, j*o*p + k*p + l) +
        '         }\n'
        '      }\n'
        '   }\n'
        '}'
    )
    c = jscode((a[i, j, k, l] + b[i, j, k, l])*c[j, k, l], assign_to=y[i])
    assert c == s


def test_jscode_loops_multiple_terms():
    from sympy.tensor import IndexedBase, Idx
    from sympy import symbols
    n, m, o, p = symbols('n m o p', integer=True)
    a = IndexedBase('a')
    b = IndexedBase('b')
    c = IndexedBase('c')
    y = IndexedBase('y')
    i = Idx('i', m)
    j = Idx('j', n)
    k = Idx('k', o)

    s0 = (
        'for (var i=0; i<m; i++){\n'
        '   y[i] = 0;\n'
        '}\n'
    )
    s1 = (
        'for (var i=0; i<m; i++){\n'
        '   for (var j=0; j<n; j++){\n'
        '      for (var k=0; k<o; k++){\n'
        '         y[i] = b[j]*b[k]*c[%s] + y[i];\n' % (i*n*o + j*o + k) +
        '      }\n'
        '   }\n'
        '}\n'
    )
    s2 = (
        'for (var i=0; i<m; i++){\n'
        '   for (var k=0; k<o; k++){\n'
        '      y[i] = b[k]*a[%s] + y[i];\n' % (i*o + k) +
        '   }\n'
        '}\n'
    )
    s3 = (
        'for (var i=0; i<m; i++){\n'
        '   for (var j=0; j<n; j++){\n'
        '      y[i] = b[j]*a[%s] + y[i];\n' % (i*n + j) +
        '   }\n'
        '}\n'
    )
    c = jscode(
        b[j]*a[i, j] + b[k]*a[i, k] + b[j]*b[k]*c[i, j, k], assign_to=y[i])
    assert (c == s0 + s1 + s2 + s3[:-1] or
            c == s0 + s1 + s3 + s2[:-1] or
            c == s0 + s2 + s1 + s3[:-1] or
            c == s0 + s2 + s3 + s1[:-1] or
            c == s0 + s3 + s1 + s2[:-1] or
            c == s0 + s3 + s2 + s1[:-1])


def test_Matrix_printing():
    # Test returning a Matrix
    mat = Matrix([x*y, Piecewise((2 + x, y>0), (y, True)), sin(z)])
    A = MatrixSymbol('A', 3, 1)
    assert jscode(mat, A) == (
        "A[0] = x*y;\n"
        "if (y > 0) {\n"
        "   A[1] = x + 2;\n"
        "}\n"
        "else {\n"
        "   A[1] = y;\n"
        "}\n"
        "A[2] = Math.sin(z);")
    # Test using MatrixElements in expressions
    expr = Piecewise((2*A[2, 0], x > 0), (A[2, 0], True)) + sin(A[1, 0]) + A[0, 0]
    assert jscode(expr) == (
        "((x > 0) ? (\n"
        "   2*A[2]\n"
        ")\n"
        ": (\n"
        "   A[2]\n"
        ")) + Math.sin(A[1]) + A[0]")
    # Test using MatrixElements in a Matrix
    q = MatrixSymbol('q', 5, 1)
    M = MatrixSymbol('M', 3, 3)
    m = Matrix([[sin(q[1,0]), 0, cos(q[2,0])],
        [q[1,0] + q[2,0], q[3, 0], 5],
        [2*q[4, 0]/q[1,0], sqrt(q[0,0]) + 4, 0]])
    assert jscode(m, M) == (
        "M[0] = Math.sin(q[1]);\n"
        "M[1] = 0;\n"
        "M[2] = Math.cos(q[2]);\n"
        "M[3] = q[1] + q[2];\n"
        "M[4] = q[3];\n"
        "M[5] = 5;\n"
        "M[6] = 2*q[4]*1/q[1];\n"
        "M[7] = 4 + Math.sqrt(q[0]);\n"
        "M[8] = 0;")
