import fractions
import re

import mpmath
import pytest

from diofant import (Add, Float, Function, I, Integer, Lambda, Matrix, Mul, Or,
                     Poly, Pow, Range, Rational, Symbol, SympifyError, Tuple,
                     Xor, evaluate, exp, false, pi, sin, sqrt, sympify, true)
from diofant.abc import _clash, _clash1, _clash2, x, y
from diofant.core.compatibility import HAS_GMPY
from diofant.core.decorators import _sympifyit
from diofant.geometry import Line, Point
from diofant.utilities.decorator import conserve_mpmath_dps


__all__ = ()


def test_sympyissue_3538():
    v = sympify('exp(x)')
    assert v == exp(x)
    assert isinstance(v, Pow)
    assert str(type(v)) == str(type(exp(x)))


def test_sympify1():
    assert sympify(None) is None
    with pytest.raises(SympifyError) as ex:
        sympify(None, strict=True)
    assert str(ex).find('SympifyError') >= 0
    assert sympify('x') == Symbol('x')
    assert sympify('   x') == Symbol('x')
    assert sympify('   x   ') == Symbol('x')
    # issue sympy/sympy#4877
    n1 = Rational(1, 2)
    assert sympify('--.5') == n1
    assert sympify('-1/2') == -n1
    assert sympify('-+--.5') == -n1
    # options to make reals into rationals
    assert sympify('2/2.6', rational=True) == Rational(10, 13)
    assert sympify('2.6/2', rational=True) == Rational(13, 10)
    assert sympify('2.6e2/17', rational=True) == Rational(260, 17)
    assert sympify('2.6e+2/17', rational=True) == Rational(260, 17)
    assert sympify('2.6e-2/17', rational=True) == Rational(26, 17000)
    assert sympify('2.1+3/4', rational=True) == \
        Rational(21, 10) + Rational(3, 4)
    assert sympify('2.234456', rational=True) == Rational(279307, 125000)
    assert sympify('2.234456e23', rational=True) == 223445600000000000000000
    assert sympify('2.234456e-23', rational=True) == \
        Rational(279307, 12500000000000000000000000000)
    assert sympify('-2.234456e-23', rational=True) == \
        Rational(-279307, 12500000000000000000000000000)
    assert sympify('12345678901/17', rational=True) == \
        Rational(12345678901, 17)
    assert sympify('1/.3 + x', rational=True) == Rational(10, 3) + x
    # make sure longs in fractions work
    assert sympify('222222222222/11111111111') == \
        Rational(222222222222, 11111111111)
    # ... or from high precision reals
    assert sympify('.1234567890123456', rational=True) == \
        Rational(19290123283179, 156250000000000)


def test_sympify_Fraction():
    value = sympify(fractions.Fraction(101, 127))
    assert value == Rational(101, 127) and type(value) is Rational


def test_sympify_gmpy():
    if HAS_GMPY:
        import gmpy2 as gmpy

        value = sympify(gmpy.mpz(1000001))
        assert value == Integer(1000001) and type(value) is Integer

        value = sympify(gmpy.mpq(101, 127))
        assert value == Rational(101, 127) and type(value) is Rational


@conserve_mpmath_dps
def test_sympify_mpmath():
    value = sympify(mpmath.mpf(1.0))
    assert value == Float(1.0) and type(value) is Float

    mpmath.mp.dps = 12
    assert sympify(
        mpmath.pi).epsilon_eq(Float('3.14159265359'), Float('1e-12')) is true
    assert sympify(
        mpmath.pi).epsilon_eq(Float('3.14159265359'), Float('1e-13')) is false

    mpmath.mp.dps = 6
    assert sympify(
        mpmath.pi).epsilon_eq(Float('3.14159'), Float('1e-5')) is true
    assert sympify(
        mpmath.pi).epsilon_eq(Float('3.14159'), Float('1e-6')) is false

    assert sympify(mpmath.mpc(1.0 + 2.0j)) == Float(1.0) + Float(2.0)*I


def test_sympify2():
    class A:
        def _diofant_(self):
            return Symbol('x')**3

    a = A()

    assert sympify(a, strict=True) == x**3
    assert sympify(a) == x**3
    assert a == x**3


def test_sympify3():
    assert sympify('x**3') == x**3
    assert sympify('x^3') == x**3
    assert sympify('x^3', convert_xor=False) == Xor(x, 3)
    assert sympify('1/2') == Rational(1, 2)

    pytest.raises(SympifyError, lambda: sympify('x**3', strict=True))
    pytest.raises(SympifyError, lambda: sympify('1/2', strict=True))


def test_sympify_keywords():
    pytest.raises(SympifyError, lambda: sympify('if'))
    pytest.raises(SympifyError, lambda: sympify('for'))
    pytest.raises(SympifyError, lambda: sympify('while'))
    pytest.raises(SympifyError, lambda: sympify('lambda'))


def test_sympify_float():
    assert sympify('1e-64') != 0
    assert sympify('1e-20000') != 0


def test_sympify_bool():
    assert sympify(True) is true
    assert sympify(False) is false


def test_sympyify_iterables():
    ans = [Rational(3, 10), Rational(1, 5)]
    assert sympify(['.3', '.2'], rational=True) == ans
    assert sympify({'.3', '.2'}, rational=True) == set(ans)
    assert sympify(('.3', '.2'), rational=True) == Tuple(*ans)
    assert sympify({x: 0, y: 1}) == {x: 0, y: 1}
    assert sympify(['1', '2', ['3', '4']]) == [Integer(1), Integer(2), [Integer(3), Integer(4)]]


def test_sympify4():
    class A:
        def _diofant_(self):
            return Symbol('x')

    a = A()

    assert sympify(a, strict=True)**3 == x**3
    assert sympify(a)**3 == x**3
    assert a == x


def test_sympify5():
    class A:
        def __str__(self):
            raise TypeError

    with pytest.raises(SympifyError) as err:
        sympify(A())
    assert re.match(r"^Sympify of expression '<diofant\.tests\.core\.test_sympify"
                    r"\.test_sympify5\.<locals>\.A object at 0x[0-9a-f]+>' failed,"
                    ' because of exception being raised:\nTypeError: $', str(err.value))


def test_sympify_text():
    assert sympify('some') == Symbol('some')
    assert sympify('core') == Symbol('core')

    assert sympify('True') is True
    assert sympify('False') is False

    assert sympify('Poly') == Poly
    assert sympify('sin') == sin


def test_sympify_function():
    assert sympify('factor(x**2-1, x)') == -(1 - x)*(x + 1)
    assert sympify('sin(pi/2)*cos(pi)') == -Integer(1)


def test_sympify_poly():
    p = (x**2 + x + 1).as_poly()

    assert sympify(p, strict=True) is p
    assert sympify(p) is p


def test_sympyissue_3595():
    assert sympify('a_') == Symbol('a_')
    assert sympify('_a') == Symbol('_a')


def test_lambda():
    x = Symbol('x')
    assert sympify('lambda: 1') == Lambda((), 1)
    assert sympify('lambda x: x') == Lambda(x, x)
    assert sympify('lambda x: 2*x') == Lambda(x, 2*x)
    assert sympify('lambda x, y: 2*x+y') == Lambda([x, y], 2*x + y)


def test_lambda_raises():
    pytest.raises(NotImplementedError, lambda: sympify('lambda *args: args'))  # args argument error
    pytest.raises(NotImplementedError, lambda: sympify('lambda **kwargs: kwargs'))  # kwargs argument error
    pytest.raises(SympifyError, lambda: sympify('lambda x = 1: x'))    # Keyword argument error
    with pytest.raises(SympifyError):
        sympify('lambda: 1', strict=True)


def test_sympify_raises():
    pytest.raises(SympifyError, lambda: sympify('fx)'))


def test_sympify_strict():
    x = Symbol('x')
    f = Function('f')

    # positive sympify
    assert sympify(x, strict=True) is x
    assert sympify(f, strict=True) is f
    assert sympify(1, strict=True) == Integer(1)
    assert sympify(0.5, strict=True) == Float('0.5')
    assert sympify(1 + 1j, strict=True) == 1.0 + I*1.0

    class A:
        def _diofant_(self):
            return Integer(5)

    a = A()
    assert sympify(a, strict=True) == Integer(5)

    # negative sympify
    pytest.raises(SympifyError, lambda: sympify('1', strict=True))
    pytest.raises(SympifyError, lambda: sympify([1, 2, 3], strict=True))


def test_sympifyit():
    x = Symbol('x')
    y = Symbol('y')

    @_sympifyit('b', NotImplemented)
    def add(a, b):
        return a + b

    assert add(x, 1) == x + 1
    assert add(x, 0.5) == x + Float('0.5')
    assert add(x, y) == x + y

    assert add(x, '1') == NotImplemented

    @_sympifyit('b')
    def add_raises(a, b):
        return a + b

    assert add_raises(x, 1) == x + 1
    assert add_raises(x, 0.5) == x + Float('0.5')
    assert add_raises(x, y) == x + y

    pytest.raises(SympifyError, lambda: add_raises(x, '1'))

    with pytest.raises(LookupError):
        @_sympifyit('x', NotImplemented)
        def spam():
            return


def test_int_float():
    class F1dot1:
        def __float__(self):
            return 1.1

    class F1dot1b:
        """
        This class is still a float, even though it also implements __int__().
        """

        def __float__(self):
            return 1.1

        def __int__(self):
            return 1

    class F1dot1c:
        """
        This class is still a float, because it implements _diofant_()
        """

        def __float__(self):
            return 1.1

        def __int__(self):
            return 1

        def _diofant_(self):
            return Float(1.1)

    class I5:
        def __int__(self):
            return 5

    class I5b:
        """
        This class implements both __int__() and __float__(), so it will be
        treated as Float in Diofant. One could change this behavior, by using
        float(a) == int(a), but deciding that integer-valued floats represent
        exact numbers is arbitrary and often not correct, so we do not do it.
        If, in the future, we decide to do it anyway, the tests for I5b need to
        be changed.
        """

        def __float__(self):
            return 5.0

        def __int__(self):
            return 5

    class I5c:
        """
        This class implements both __int__() and __float__(), but also
        a _diofant_() method, so it will be Integer.
        """

        def __float__(self):
            return 5.0

        def __int__(self):
            return 5

        def _diofant_(self):
            return Integer(5)

    i5 = I5()
    i5b = I5b()
    i5c = I5c()
    f1_1 = F1dot1()
    f1_1b = F1dot1b()
    f1_1c = F1dot1c()
    assert sympify(i5) == 5
    assert isinstance(sympify(i5), Integer)
    assert sympify(i5b) == 5
    assert isinstance(sympify(i5b), Float)
    assert sympify(i5c) == 5
    assert isinstance(sympify(i5c), Integer)
    assert abs(sympify(f1_1) - 1.1) < 1e-5
    assert abs(sympify(f1_1b) - 1.1) < 1e-5
    assert abs(sympify(f1_1c) - 1.1) < 1e-5

    assert sympify(i5, strict=True) == 5
    assert isinstance(sympify(i5, strict=True), Integer)
    assert sympify(i5b, strict=True) == 5
    assert isinstance(sympify(i5b, strict=True), Float)
    assert sympify(i5c, strict=True) == 5
    assert isinstance(sympify(i5c, strict=True), Integer)
    assert abs(sympify(f1_1, strict=True) - 1.1) < 1e-5
    assert abs(sympify(f1_1b, strict=True) - 1.1) < 1e-5
    assert abs(sympify(f1_1c, strict=True) - 1.1) < 1e-5


def test_evaluate_false():
    cases = {
        '2 + 3': Add(2, 3, evaluate=False),
        '2**2 / 3': Mul(Pow(2, 2, evaluate=False), Pow(3, -1, evaluate=False), evaluate=False),
        '2 + 3 * 5': Add(2, Mul(3, 5, evaluate=False), evaluate=False),
        '2 - 3 * 5': Add(2, -Mul(3, 5, evaluate=False), evaluate=False),
        '1 / 3': Mul(1, Pow(3, -1, evaluate=False), evaluate=False),
        'True | False': Or(True, False, evaluate=False),
        '1 + 2 + 3 + 5*3 + integrate(x)': Add(1, 2, 3, Mul(5, 3, evaluate=False), x**2/2, evaluate=False),
        '2 * 4 * 6 + 8': Add(Mul(2, 4, 6, evaluate=False), 8, evaluate=False),
    }
    for case, result in cases.items():
        assert sympify(case, evaluate=False) == result


def test_sympyissue_4133():
    a = sympify('Integer(4)')

    assert a == Integer(4)
    assert a.is_Integer


def test_sympyissue_3982():
    a = [3, 2.0]
    assert sympify(a) == [Integer(3), Float(2.0)]
    assert sympify(tuple(a)) == Tuple(Integer(3), Float(2.0))
    assert sympify(set(a)) == {Integer(3), Float(2.0)}


def test_S_sympify():
    assert Rational(1, 2) == sympify(1)/2
    assert sqrt(-2) == sqrt(2)*I


def test_sympyissue_4788():
    assert repr(sympify(1.0 + 0J)) == repr(Float(1.0)) == repr(Float(1.0))


def test_sympyissue_4798_None():
    assert sympify(None) is None


def test_sympyissue_3218():
    assert sympify('x+\ny') == x + y


def test_sympyissue_4988_builtins():
    C = Symbol('C')
    vars = {}
    vars['C'] = C
    exp1 = sympify('C')
    assert exp1 == C  # Make sure it did not get mixed up with diofant.C

    exp2 = sympify('C', vars)
    assert exp2 == C  # Make sure it did not get mixed up with diofant.C


def test_geometry():
    p = sympify(Point(0, 1))
    assert p == Point(0, 1) and isinstance(p, Point)
    L = sympify(Line(p, (1, 0)))
    assert L == Line((0, 1), (1, 0)) and isinstance(L, Line)


def test_sympyissue_6540_6552():
    assert sympify('[[1/3,2], (2/5,)]') == [[Rational(1, 3), 2], (Rational(2, 5),)]
    assert sympify('[[2/6,2], (2/4,)]') == [[Rational(1, 3), 2], (Rational(1, 2),)]
    assert sympify('[[[2*(1)]]]') == [[[2]]]
    assert sympify('Matrix([2*(1)])') == Matrix([2])


def test_sympyissue_6046():
    assert str(sympify('Q & C', locals=_clash1)) == 'C & Q'
    assert str(sympify('pi(x)', locals=_clash2)) == 'pi(x)'
    assert str(sympify('pi(C, Q)', locals=_clash)) == 'pi(C, Q)'
    locals = {}
    exec('from diofant.abc import S, O', locals)
    assert str(sympify('O&S', locals)) == 'O & S'


def test_sympyissue_8821_highprec_from_str():
    s = str(pi.evalf(128))
    p = sympify(s)
    assert abs(sin(p)) < 1e-127


def test_Range():
    assert sympify(range(10)) == Range(10)
    assert sympify(range(10), strict=True) == Range(10)


def test_sympyissue_10773():
    with evaluate(False):
        ans = Mul(Integer(-10), Pow(Integer(5), Integer(-1)))
    assert sympify('-10/5', evaluate=False) == ans

    with evaluate(False):
        ans = Mul(Integer(-10), Pow(Integer(-5), Integer(-1)))
    assert sympify('-10/-5', evaluate=False) == ans
