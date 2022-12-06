import pytest

from diofant import (Add, Basic, Derivative, DiracDelta, Dummy, E, Float,
                     Function, Ge, Gt, Heaviside, I, Integer, Integral, Le, Lt,
                     Max, Mul, Number, NumberSymbol, O, Piecewise, Pow,
                     Rational, Si, Subs, Sum, Symbol, Tuple, Wild,
                     WildFunction, apart, cancel, cbrt, collect, combsimp, cos,
                     default_sort_key, diff, exp, exp_polar, expand,
                     expand_multinomial, factor, factorial, false, gamma, log,
                     lucas, nan, nsimplify, oo, pi, posify, powsimp, radsimp,
                     ratsimp, root, simplify, sin, sqrt, symbols, sympify, tan,
                     tanh, together, trigsimp, true, zoo)
from diofant.abc import a, b, c, n, r, t, u, x, y, z
from diofant.core.function import AppliedUndef
from diofant.solvers.utils import checksol


__all__ = ()


class DummyNumber:
    """
    Minimal implementation of a number that works with Diofant.

    If one has a Number class (e.g. Sage Integer, or some other custom class)
    that one wants to work well with Diofant, one has to implement at least the
    methods of this class DummyNumber, resp. its subclasses I5 and F1dot1.

    Basically, one just needs to implement either __int__() or __float__() and
    then one needs to make sure that the class works with Python integers and
    with itself.
    """

    def __radd__(self, other):
        if isinstance(other, (int, float)):
            return other + self.number
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, (int, float, DummyNumber)):
            return self.number + other
        return NotImplemented

    def __rsub__(self, other):
        if isinstance(other, (int, float)):
            return other - self.number
        return NotImplemented

    def __sub__(self, other):
        if isinstance(other, (int, float, DummyNumber)):
            return self.number - other
        return NotImplemented

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return other * self.number
        return NotImplemented

    def __mul__(self, other):
        if isinstance(other, (int, float, DummyNumber)):
            return self.number * other
        return NotImplemented

    def __rtruediv__(self, other):
        if isinstance(other, (int, float)):
            return other / self.number
        return NotImplemented

    def __truediv__(self, other):
        if isinstance(other, (int, float, DummyNumber)):
            return self.number / other
        return NotImplemented

    def __rpow__(self, other):
        if isinstance(other, (int, float)):
            return other ** self.number
        return NotImplemented

    def __pow__(self, other):
        if isinstance(other, (int, float, DummyNumber)):
            return self.number ** other
        return NotImplemented

    def __pos__(self):
        return self.number

    def __neg__(self):
        return - self.number


class I5(DummyNumber):
    number = 5

    def __int__(self):
        return self.number


class F1dot1(DummyNumber):
    number = 1.1

    def __float__(self):
        return self.number


i5 = I5()
f1_1 = F1dot1()

# basic diofant objects
basic_objs = [
    Integer(2),
    Float('1.3'),
    x,
    y,
    pow(x, y)*y,
]

# all supported objects
all_objs = basic_objs + [
    5,
    5.5,
    i5,
    f1_1
]


def dotest(s):
    for x in all_objs:
        for y in all_objs:
            s(x, y)
    return True


def test_basic():
    def j(a, b):
        # pylint: disable=pointless-statement
        a
        +a
        -a
        a + b
        a - b
        a*b
        a/b
        a**b
    assert dotest(j)


def test_ibasic():
    def s(a, b):
        # pylint: disable=pointless-statement
        x = a
        x += b
        x = a
        x -= b
        x = a
        x *= b
        x = a
        x /= b
    assert dotest(s)


def test_relational():
    assert (pi < 3) is false
    assert (pi <= 3) is false
    assert (pi > 3) is true
    assert (pi >= 3) is true
    assert (-pi < 3) is true
    assert (-pi <= 3) is true
    assert (-pi > 3) is false
    assert (-pi >= 3) is false
    r = Symbol('r', extended_real=True)
    assert (r - 2 < r - 3) is false
    assert Lt(x + I, x + I + 2).func == Lt  # issue sympy/sympy#8288


def test_relational_assumptions():
    m1 = Symbol('m1', nonnegative=False)
    m2 = Symbol('m2', positive=False)
    m3 = Symbol('m3', nonpositive=False)
    m4 = Symbol('m4', negative=False)
    assert (m1 < 0) == Lt(m1, 0)
    assert (m2 <= 0) == Le(m2, 0)
    assert (m3 > 0) == Gt(m3, 0)
    assert (m4 >= 0) == Ge(m4, 0)
    m1 = Symbol('m1', nonnegative=False, extended_real=True)
    m2 = Symbol('m2', positive=False, extended_real=True)
    m3 = Symbol('m3', nonpositive=False, extended_real=True)
    m4 = Symbol('m4', negative=False, extended_real=True)
    assert (m1 < 0) is true
    assert (m2 <= 0) is true
    assert (m3 > 0) is true
    assert (m4 >= 0) is true
    m1 = Symbol('m1', negative=True)
    m2 = Symbol('m2', nonpositive=True)
    m3 = Symbol('m3', positive=True)
    m4 = Symbol('m4', nonnegative=True)
    assert (m1 < 0) is true
    assert (m2 <= 0) is true
    assert (m3 > 0) is true
    assert (m4 >= 0) is true
    m1 = Symbol('m1', negative=False, extended_real=True)
    m2 = Symbol('m2', nonpositive=False, extended_real=True)
    m3 = Symbol('m3', positive=False, extended_real=True)
    m4 = Symbol('m4', nonnegative=False, extended_real=True)
    assert (m1 < 0) is false
    assert (m2 <= 0) is false
    assert (m3 > 0) is false
    assert (m4 >= 0) is false


def test_relational_noncommutative():
    A, B = symbols('A,B', commutative=False)
    assert (A < B) == Lt(A, B)
    assert (A <= B) == Le(A, B)
    assert (A > B) == Gt(A, B)
    assert (A >= B) == Ge(A, B)


def test_basic_nostr():
    for obj in basic_objs:
        pytest.raises(TypeError, lambda: obj + '1')
        pytest.raises(TypeError, lambda: obj - '1')
        if obj == 2:
            assert obj * '1' == '11'
        else:
            pytest.raises(TypeError, lambda: obj * '1')
        pytest.raises(TypeError, lambda: obj / '1')
        pytest.raises(TypeError, lambda: obj ** '1')


def test_series0():
    # issue sympy/sympy#7231
    f = Function('f')
    ans1 = f(x).series(x, a)
    _xi_1 = ans1.atoms(Dummy).pop()
    res = (f(a) + (-a + x)*Subs(Derivative(f(_xi_1), _xi_1), (_xi_1, a)) +
           (-a + x)**2*Subs(Derivative(f(_xi_1), _xi_1, _xi_1), (_xi_1, a))/2 +
           (-a + x)**3*Subs(Derivative(f(_xi_1), _xi_1, _xi_1, _xi_1), (_xi_1, a))/6 +
           (-a + x)**4*Subs(Derivative(f(_xi_1), _xi_1, _xi_1, _xi_1, _xi_1), (_xi_1, a))/24 +
           (-a + x)**5*Subs(Derivative(f(_xi_1), _xi_1, _xi_1, _xi_1, _xi_1, _xi_1),
                            (_xi_1, a))/120 + O((-a + x)**6, (x, a)))
    assert res == ans1
    ans2 = f(x).series(x, a)
    assert res == ans2


def test_series_expansion_for_uniform_order():
    assert (1/x + y + x).series(x, 0, 0) == 1/x + O(1, x)
    assert (1/x + y + x).series(x, 0, 1) == 1/x + y + O(x)
    assert (1/x + 1 + x).series(x, 0, 0) == 1/x + O(1, x)
    assert (1/x + 1 + x).series(x, 0, 1) == 1/x + 1 + O(x)
    assert (1/x + x).series(x, 0, 0) == 1/x + O(1, x)
    assert (1/x + y + y*x + x).series(x, 0, 0) == 1/x + O(1, x)
    assert (1/x + y + y*x + x).series(x, 0, 1) == 1/x + y + O(x)


def test_as_leading_term():
    pytest.raises(ValueError, lambda: (1/x + x).as_leading_term(1))
    assert (3 + 2*x**(log(3)/log(2) - 1)).as_leading_term(x) == 3
    assert (1/x**2 + 1 + x + x**2).as_leading_term(x) == 1/x**2
    assert (1/x + 1 + x + x**2).as_leading_term(x) == 1/x
    assert (x**2 + 1/x).as_leading_term(x) == 1/x
    assert (1 + x**2).as_leading_term(x) == 1
    assert (x + 1).as_leading_term(x) == 1
    assert (x + x**2).as_leading_term(x) == x
    assert (x**2).as_leading_term(x) == x**2
    assert (x + oo).as_leading_term(x) == oo

    assert (x*cos(1)*cos(1 + sin(1)) + sin(1 + sin(1))).as_leading_term(x) == \
        sin(1 + sin(1))

    assert (2 + pi + x).as_leading_term(x) == 2 + pi
    assert (2*x + pi*x + x**2).as_leading_term(x) == (2 + pi)*x

    # see issue sympy/sympy#6843
    n = Symbol('n', integer=True, positive=True)
    r = -n**3/(2*n**2 + 4*n + 2) - n**2/(n**2 + 2*n + 1) + \
        n**2/(n + 1) - n/(2*n**2 + 4*n + 2) + n/(n*x + x) + 2*n/(n + 1) - \
        1 + 1/(n*x + x) + 1/(n + 1) - 1/x
    assert r.as_leading_term(x).cancel() == n/2

    # see issue sympy/sympy#9075
    assert (6**(x + 1)).as_leading_term(x) == 6
    assert (6**(x + n)).as_leading_term(x) == 6**n

    # issue sympy/sympy#17847
    assert (1 - cos(x)).as_leading_term(x) == x**2/2
    assert (1 - cos(x) + x**6).as_leading_term(x) == x**2/2
    assert (1 + cos(x) + x**6).as_leading_term(x) == 2
    assert (sin(x) - x + x**7).as_leading_term(x) == -x**3/6


def test_as_leading_term_stub():
    class Foo(Function):
        pass
    assert Foo(1/x).as_leading_term(x) == Foo(1/x)
    assert Foo(1).as_leading_term(x) == Foo(1)
    pytest.raises(NotImplementedError, lambda: Foo(x).as_leading_term(x))


def test_atoms():
    assert x.atoms() == {x}
    assert (1 + x).atoms() == {x, 1}

    assert (1 + 2*cos(x)).atoms(Symbol) == {x}
    assert (1 + 2*cos(x)).atoms(Symbol, Number) == {1, 2, x}

    assert (2*(x**(y**x))).atoms() == {2, x, y}

    assert Rational(1, 2).atoms() == {Rational(1, 2)}
    assert Rational(1, 2).atoms(Symbol) == set()

    assert sin(oo).atoms(oo) == {oo}

    assert Integer(0).as_poly(x).atoms() == {0, x}
    assert Integer(1).as_poly(x).atoms() == {1, x}

    assert x.as_poly().atoms() == {x}
    assert x.as_poly(x, y).atoms() == {x, y}
    assert (x + y).as_poly().atoms() == {x, y}
    assert (x + y).as_poly(x, y, z).atoms() == {x, y, z}
    assert (x + y*t).as_poly(x, y, z).atoms() == {t, x, y, z}

    assert (I*pi).atoms(NumberSymbol) == {pi}
    assert (I*pi).atoms(NumberSymbol, I) == {pi, I}
    assert (I*pi).atoms(I, NumberSymbol) == {pi, I}

    assert exp(exp(x)).atoms(Pow) == {exp(exp(x)), exp(x)}
    assert (1 + x*(2 + y) + exp(3 + z)).atoms(Add) == {1 + x*(2 + y) + exp(3 + z),
                                                       2 + y, 3 + z}

    # issue sympy/sympy#6132
    f = Function('f')
    e = (f(x) + sin(x) + 2)
    assert e.atoms(AppliedUndef) == {f(x)}
    assert e.atoms(AppliedUndef, Function) == {f(x), sin(x)}
    assert e.atoms(Function) == {f(x), sin(x)}
    assert e.atoms(AppliedUndef, Number) == {f(x), 2}
    assert e.atoms(Function, Number) == {2, sin(x), f(x)}


def test_is_polynomial():
    k = Symbol('k', nonnegative=True, integer=True)

    assert Integer(2).is_polynomial(x, y, z) is True
    assert pi.is_polynomial(x, y, z) is True

    assert x.is_polynomial(x) is True
    assert x.is_polynomial(y) is True

    assert (x**2).is_polynomial(x) is True
    assert (x**2).is_polynomial(y) is True

    assert (x**(-2)).is_polynomial(x) is False
    assert (x**(-2)).is_polynomial(y) is True

    assert (2**x).is_polynomial(x) is False
    assert (2**x).is_polynomial(y) is True

    assert (x**k).is_polynomial(x) is False
    assert (x**k).is_polynomial(k) is False
    assert (x**x).is_polynomial(x) is False
    assert (k**k).is_polynomial(k) is False
    assert (k**x).is_polynomial(k) is False

    assert (x**(-k)).is_polynomial(x) is False
    assert ((2*x)**k).is_polynomial(x) is False

    assert (x**2 + 3*x - 8).is_polynomial(x) is True
    assert (x**2 + 3*x - 8).is_polynomial(y) is True

    assert (x**2 + 3*x - 8).is_polynomial() is True

    assert sqrt(x).is_polynomial(x) is False
    assert (sqrt(x)**3).is_polynomial(x) is False

    assert (x**2 + 3*x*sqrt(y) - 8).is_polynomial(x) is True
    assert (x**2 + 3*x*sqrt(y) - 8).is_polynomial(y) is False

    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(2)).is_polynomial() is True
    assert ((x**2)*(y**2) + x*(y**2) + y*x + exp(x)).is_polynomial() is False

    assert (
        (x**2)*(y**2) + x*(y**2) + y*x + exp(2)).is_polynomial(x, y) is True
    assert (
        (x**2)*(y**2) + x*(y**2) + y*x + exp(x)).is_polynomial(x, y) is False


def test_is_rational_function():
    assert Integer(1).is_rational_function() is True
    assert Integer(1).is_rational_function(x) is True

    assert Rational(17, 54).is_rational_function() is True
    assert Rational(17, 54).is_rational_function(x) is True

    assert (12/x).is_rational_function() is True
    assert (12/x).is_rational_function(x) is True

    assert (x/y).is_rational_function() is True
    assert (x/y).is_rational_function(x) is True
    assert (x/y).is_rational_function(x, y) is True

    assert (x**2 + 1/x/y).is_rational_function() is True
    assert (x**2 + 1/x/y).is_rational_function(x) is True
    assert (x**2 + 1/x/y).is_rational_function(x, y) is True

    assert (sin(y)/x).is_rational_function() is False
    assert (sin(y)/x).is_rational_function(y) is False
    assert (sin(y)/x).is_rational_function(x) is True
    assert (sin(y)/x).is_rational_function(x, y) is False

    assert nan.is_rational_function() is False
    assert (+oo).is_rational_function() is False
    assert (-oo).is_rational_function() is False
    assert zoo.is_rational_function() is False


def test_is_algebraic_expr():
    assert sqrt(3).is_algebraic_expr(x) is True
    assert sqrt(3).is_algebraic_expr() is True

    eq = cbrt((1 + x**2)/(1 - y**2))
    assert eq.is_algebraic_expr(x) is True
    assert eq.is_algebraic_expr(y) is True

    assert (sqrt(x) + y**Rational(2, 3)).is_algebraic_expr(x) is True
    assert (sqrt(x) + y**Rational(2, 3)).is_algebraic_expr(y) is True
    assert (sqrt(x) + y**Rational(2, 3)).is_algebraic_expr() is True

    assert (cos(y)/sqrt(x)).is_algebraic_expr() is False
    assert (cos(y)/sqrt(x)).is_algebraic_expr(x) is True
    assert (cos(y)/sqrt(x)).is_algebraic_expr(y) is False
    assert (cos(y)/sqrt(x)).is_algebraic_expr(x, y) is False

    a = sqrt(exp(x)**2 + 2*exp(x) + 1)/(exp(x) + 1)
    assert a.is_algebraic_expr(x) is False


def test_SAGE1():
    # see https://github.com/sympy/sympy/issues/3346
    class MyInt:
        def _diofant_(self):
            return Integer(5)
    m = MyInt()
    e = Integer(2)*m
    assert e == 10

    pytest.raises(TypeError, lambda: Integer(2)*MyInt)


def test_SAGE2():
    class MyInt:
        def __int__(self):
            return 5
    assert sympify(MyInt()) == 5
    e = Integer(2)*MyInt()
    assert e == 10

    pytest.raises(TypeError, lambda: Integer(2)*MyInt)


def test_SAGE3():
    class MySymbol:
        def __rmul__(self, other):
            return 'mys', other, self

    o = MySymbol()
    e = x*o

    assert e == ('mys', x, o)


def test_len():
    e = x*y
    assert len(e.args) == 2
    e = x + y + z
    assert len(e.args) == 3


def test_doit():
    a = Integral(x**2, x)

    assert isinstance(a.doit(), Integral) is False

    assert isinstance(a.doit(integrals=True), Integral) is False
    assert isinstance(a.doit(integrals=False), Integral) is True

    assert (2*Integral(x, x)).doit() == x**2


def test_args():
    assert (x*y).args in ((x, y), (y, x))
    assert (x + y).args in ((x, y), (y, x))
    assert (x*y + 1).args in ((x*y, 1), (1, x*y))
    assert sin(x*y).args == (x*y,)
    assert sin(x*y).args[0] == x*y
    assert (x**y).args == (x, y)
    assert (x**y).args[0] == x
    assert (x**y).args[1] == y


def test_sympyissue_3757():
    A, B, C = symbols('A,B,C', commutative=False)
    assert A*B - B*A != 0
    assert (A*(A + B)*B).expand() == A**2*B + A*B**2
    assert (A*(A + B + C)*B).expand() == A**2*B + A*B**2 + A*C*B


def test_as_numer_denom():
    assert nan.as_numer_denom() == (nan, 1)
    assert oo.as_numer_denom() == (oo, 1)
    assert (-oo).as_numer_denom() == (-oo, 1)
    assert zoo.as_numer_denom() == (zoo, 1)
    assert (-zoo).as_numer_denom() == (zoo, 1)

    assert x.as_numer_denom() == (x, 1)
    assert (1/x).as_numer_denom() == (1, x)
    assert (x/y).as_numer_denom() == (x, y)
    assert (x/2).as_numer_denom() == (x, 2)
    assert (x*y/z).as_numer_denom() == (x*y, z)
    assert (x/(y*z)).as_numer_denom() == (x, y*z)
    assert Rational(1, 2).as_numer_denom() == (1, 2)
    assert (1/y**2).as_numer_denom() == (1, y**2)
    assert (x/y**2).as_numer_denom() == (x, y**2)
    assert ((x**2 + 1)/y).as_numer_denom() == (x**2 + 1, y)
    assert (x*(y + 1)/y**7).as_numer_denom() == (x*(y + 1), y**7)
    assert (x**-2).as_numer_denom() == (1, x**2)
    assert (a/x + b/2/x + c/3/x).as_numer_denom() == \
        (6*a + 3*b + 2*c, 6*x)
    assert (a/x + b/2/x + c/3/y).as_numer_denom() == \
        (2*c*x + y*(6*a + 3*b), 6*x*y)
    assert (a/x + b/2/x + c/.5/x).as_numer_denom() == \
        (2*a + b + 4.0*c, 2*x)
    # this should take no more than a few seconds
    assert int(log(Add(*[Dummy()/i/x for i in range(1, 705)]
                       ).as_numer_denom()[1]/x).evalf(4)) == 705
    for i in (oo, -oo, zoo):
        assert (i + x/3).as_numer_denom() == \
            (x + i, 3)
    assert (oo + x/3 + y/4).as_numer_denom() == \
        (4*x + 3*y + oo, 12)
    assert (oo*x + zoo*y).as_numer_denom() == \
        (zoo*y + oo*x, 1)

    A, B, C = symbols('A,B,C', commutative=False)

    assert (A*B*C**-1).as_numer_denom() == (A*B*C**-1, 1)
    assert (A*B*C**-1/x).as_numer_denom() == (A*B*C**-1, x)
    assert (C**-1*A*B).as_numer_denom() == (C**-1*A*B, 1)
    assert (C**-1*A*B/x).as_numer_denom() == (C**-1*A*B, x)
    assert ((A*B*C)**-1).as_numer_denom() == ((A*B*C)**-1, 1)
    assert ((A*B*C)**-1/x).as_numer_denom() == ((A*B*C)**-1, x)


def test_as_independent():
    assert (2*x*sin(x) + y + x).as_independent(x) == (y, x + 2*x*sin(x))
    assert (2*x*sin(x) + y + x).as_independent(y) == (x + 2*x*sin(x), y)

    assert (2*x*sin(x) + y + x).as_independent(x, y) == (0, y + x + 2*x*sin(x))

    assert (x*sin(x)*cos(y)).as_independent(x) == (cos(y), x*sin(x))
    assert (x*sin(x)*cos(y)).as_independent(y) == (x*sin(x), cos(y))

    assert (x*sin(x)*cos(y)).as_independent(x, y) == (1, x*sin(x)*cos(y))

    assert (sin(x)).as_independent(x) == (1, sin(x))
    assert (sin(x)).as_independent(y) == (sin(x), 1)

    assert (2*sin(x)).as_independent(x) == (2, sin(x))
    assert (2*sin(x)).as_independent(y) == (2*sin(x), 1)

    # issue sympy/sympy#4903 = 4865b
    n1, n2, n3 = symbols('n1 n2 n3', commutative=False)
    assert (n1 + n1*n2).as_independent(n2) == (n1, n1*n2)
    assert (n2*n1 + n1*n2).as_independent(n2) == (0, n1*n2 + n2*n1)
    assert (n1*n2*n1).as_independent(n2) == (n1, n2*n1)
    assert (n1*n2*n1).as_independent(n1) == (1, n1*n2*n1)

    assert (3*x).as_independent(x, as_Add=True) == (0, 3*x)
    assert (3*x).as_independent(x, as_Add=False) == (3, x)
    assert (3 + x).as_independent(x, as_Add=True) == (3, x)
    assert (3 + x).as_independent(x, as_Add=False) == (1, 3 + x)

    # issue sympy/sympy#5479
    assert (3*x).as_independent(Symbol) == (3, x)

    # issue sympy/sympy#5648
    assert (n1*x*y).as_independent(x) == (n1*y, x)
    assert ((x + n1)*(x - y)).as_independent(x) == (1, (x + n1)*(x - y))
    assert ((x + n1)*(x - y)).as_independent(y) == (x + n1, x - y)
    assert (DiracDelta(x - z)*DiracDelta(x - y)).as_independent(x) \
        == (1, DiracDelta(x - z)*DiracDelta(x - y))
    assert (x*y*n1*n2*n3).as_independent(n2) == (x*y*n1, n2*n3)
    assert (x*y*n1*n2*n3).as_independent(n1) == (x*y, n1*n2*n3)
    assert (x*y*n1*n2*n3).as_independent(n3) == (x*y*n1*n2, n3)
    assert (DiracDelta(x - z)*DiracDelta(y - z)*DiracDelta(x - t)).as_independent(y) == \
           (DiracDelta(x - z)*DiracDelta(x - t), DiracDelta(y - z))

    # issue sympy/sympy#5784
    assert (x + Integral(x, (x, 1, 2))).as_independent(x, strict=True) == \
           (Integral(x, (x, 1, 2)), x)

    assert (x*y).as_independent(z, as_Add=True) == (x*y, 0)


@pytest.mark.xfail
def test_call_2():
    f = Function('f')
    assert (2*f)(x) == 2*f(x)  # UndefinedFunction does not subclass Expr yet


def test_replace():
    f = log(sin(x)) + tan(sin(x**2))

    assert f.replace(sin, cos) == log(cos(x)) + tan(cos(x**2))
    pytest.raises(TypeError, lambda: f.replace(sin, x))
    assert f.replace(
        sin, lambda a: sin(2*a)) == log(sin(2*x)) + tan(sin(2*x**2))

    a = Wild('a')
    b = Wild('b')

    assert f.replace(sin(a), cos(a)) == log(cos(x)) + tan(cos(x**2))
    assert f.replace(
        sin(a), lambda a: sin(2*a)) == log(sin(2*x)) + tan(sin(2*x**2))
    pytest.raises(TypeError, lambda: f.replace(sin(a), None))
    # test exact
    assert (2*x).replace(a*x + b, b - a, exact=True) == 2*x
    assert (2*x).replace(a*x + b, b - a) == 2/x
    assert (2*x).replace(a*x + b, lambda a, b: b - a, exact=True) == 2*x
    assert (2*x).replace(a*x + b, lambda a, b: b - a) == 2/x
    assert (2*x - 1).replace(a*x + b, lambda a, b: b - a, exact=True) == -3

    g = 2*sin(x**3)

    assert g.replace(
        lambda expr: expr.is_Number, lambda expr: expr**2) == 4*sin(x**9)

    assert sin(x).replace(cos, sin) == sin(x)

    def cond(x):
        return x.is_Mul

    assert (x**2 + O(x**3)).replace(Pow, lambda b, e: b**e/e) == x**2/2 + O(x**3)
    assert (x*(x*y + 3)).replace(lambda x: x.is_Mul, lambda x: 2 + x) == \
        x*(x*y + 5) + 2
    e = (x*y + 1)*(2*x*y + 1) + 1
    assert x.replace(x, y) == y
    assert (x + 1).replace(1, 2) == x + 2
    pytest.raises(TypeError, lambda: e.replace(cond, x))

    # https://groups.google.com/forum/#!topic/sympy/8wCgeC95tz0
    n1, n2, n3 = symbols('n1:4', commutative=False)
    f = Function('f')
    assert (n1*f(n2)).replace(f, lambda x: x) == n1*n2
    assert (n3*f(n2)).replace(f, lambda x: x) == n3*n2

    # for test coverage
    pytest.raises(TypeError, lambda: x.replace(None, y))

    assert (2*x + y).replace(a*x + b, b - a, exact=True) == y - 2


def test_find():
    expr = (x + y + 2 + sin(3*x))

    assert expr.find(lambda u: u.is_Integer) == {2: 1, 3: 1}
    assert expr.find(lambda u: u.is_Symbol) == {x: 2, y: 1}

    assert expr.find(Integer) == {2: 1, 3: 1}
    assert expr.find(Symbol) == {x: 2, y: 1}

    a = Wild('a')

    expr = sin(sin(x)) + sin(x) + cos(x) + x

    assert expr.find(lambda u: type(u) is sin) == {sin(x): 2, sin(sin(x)): 1}
    assert expr.find(sin(a)) == {sin(x): 2, sin(sin(x)): 1}
    assert expr.find(sin) == {sin(x): 2, sin(sin(x)): 1}


def test_count():
    expr = (x + y + 2 + sin(3*x))

    assert expr.count(lambda u: u.is_Integer) == 2
    assert expr.count(lambda u: u.is_Symbol) == 3

    assert expr.count(Integer) == 2
    assert expr.count(Symbol) == 3
    assert expr.count(2) == 1

    a = Wild('a')

    assert expr.count(sin) == 1
    assert expr.count(sin(a)) == 1
    assert expr.count(lambda u: type(u) is sin) == 1


def test_has_basics():
    f = Function('f')
    g = Function('g')
    p = Wild('p')

    assert sin(x).has(x)
    assert sin(x).has(sin)
    assert not sin(x).has(y)
    assert not sin(x).has(cos)
    assert f(x).has(x)
    assert f(x).has(f)
    assert not f(x).has(y)
    assert not f(x).has(g)

    assert f(x).diff(x).has(x)
    assert f(x).diff(x).has(f)
    assert f(x).diff(x).has(Derivative)
    assert not f(x).diff(x).has(y)
    assert not f(x).diff(x).has(g)
    assert not f(x).diff(x).has(sin)

    assert (x**2).has(Symbol)
    assert not (x**2).has(Wild)
    assert (2*p).has(Wild)

    assert not x.has()


def test_has_multiple():
    f = x**2*y + sin(2**t + log(z))

    assert f.has(x)
    assert f.has(y)
    assert f.has(z)
    assert f.has(t)

    assert not f.has(u)

    assert f.has(x, y, z, t)
    assert f.has(x, y, z, t, u)

    i = Integer(4400)

    assert not i.has(x)

    assert (i*x**i).has(x)
    assert not (i*y**i).has(x)
    assert (i*y**i).has(x, y)
    assert not (i*y**i).has(x, z)


def test_has_piecewise():
    f = (x*y + 3/y)**(3 + 2)
    g = Function('g')
    h = Function('h')
    p = Piecewise((g(x), x < -1), (1, x <= 1), (f, True))

    assert p.has(x)
    assert p.has(y)
    assert not p.has(z)
    assert p.has(1)
    assert p.has(3)
    assert not p.has(4)
    assert p.has(f)
    assert p.has(g)
    assert not p.has(h)


def test_has_iterative():
    A, B, C = symbols('A,B,C', commutative=False)
    f = x*gamma(x)*sin(x)*exp(x*y)*A*B*C*cos(x*A*B)

    assert f.has(x)
    assert f.has(x*y)
    assert f.has(x*sin(x))
    assert not f.has(x*sin(y))
    assert f.has(x*A)
    assert f.has(x*A*B)
    assert not f.has(x*A*C)
    assert f.has(x*A*B*C)
    assert not f.has(x*A*C*B)
    assert f.has(x*sin(x)*A*B*C)
    assert not f.has(x*sin(x)*A*C*B)
    assert not f.has(x*sin(y)*A*B*C)
    assert f.has(x*gamma(x))
    assert not f.has(x + sin(x))

    assert (x & y & z).has(x & z)


def test_has_integrals():
    f = Integral(x**2 + sin(x*y*z), (x, 0, x + y + z))

    assert f.has(x + y)
    assert f.has(x + z)
    assert f.has(y + z)

    assert f.has(x*y)
    assert f.has(x*z)
    assert f.has(y*z)

    assert not f.has(2*x + y)
    assert not f.has(2*x*y)


def test_has_tuple():
    f = Function('f')
    g = Function('g')
    h = Function('h')

    assert Tuple(x, y).has(x)
    assert not Tuple(x, y).has(z)
    assert Tuple(f(x), g(x)).has(x)
    assert not Tuple(f(x), g(x)).has(y)
    assert Tuple(f(x), g(x)).has(f)
    assert Tuple(f(x), g(x)).has(f(x))
    assert not Tuple(f, g).has(x)
    assert Tuple(f, g).has(f)
    assert not Tuple(f, g).has(h)
    assert Tuple(True).has(True) is True  # .has(1) will also be True


def test_has_polys():
    poly = (x**2 + x*y*sin(z)).as_poly(x, y, t)

    assert poly.has(x)
    assert poly.has(x, y, z)
    assert poly.has(x, y, z, t)

    # the following morphs from Add to Mul during processing
    assert Add(0, (x + y)/z/-2,
               evaluate=False).as_numer_denom() == (-x - y, 2*z)


def test_as_poly_as_expr():
    f = x**2 + 2*x*y

    assert f.as_poly().as_expr() == f
    assert f.as_poly(x, y).as_expr() == f

    assert (f + sin(x)).as_poly(x, y) is None

    p = f.as_poly(x, y)

    assert p.as_poly() == p


def test_nonzero():
    assert bool(Integer(0)) is False
    assert bool(Integer(1)) is True
    assert bool(x) is True
    assert bool(x + y) is True
    assert bool(x - x) is False
    assert bool(x*y) is True
    assert bool(x*1) is True
    assert bool(x*0) is False


def test_is_number():
    assert Float(3.14).is_number is True
    assert Integer(737).is_number is True
    assert Rational(3, 2).is_number is True
    assert Integer(8).is_number is True
    assert x.is_number is False
    assert (2*x).is_number is False
    assert (x + y).is_number is False
    assert log(2).is_number is True
    assert log(x).is_number is False
    assert (2 + log(2)).is_number is True
    assert (8 + log(2)).is_number is True
    assert (2 + log(x)).is_number is False
    assert (8 + log(2) + x).is_number is False
    assert (1 + x**2/x - x).is_number is True
    assert Tuple(1).is_number is False
    assert Add(2, x).is_number is False
    assert Mul(3, 4).is_number is True
    assert Pow(log(2), 2).is_number is True
    assert oo.is_number is True
    g = WildFunction('g')
    assert g.is_number is False
    assert (2*g).is_number is False
    assert (x**2).subs({x: 3}).is_number is True

    # test extensibility of .is_number
    # on subinstances of Basic
    class A(Basic):
        pass
    a = A()
    assert a.is_number is False


def test_as_coeff_add():
    assert Integer(2).as_coeff_add() == (2, ())
    assert Float(3.0).as_coeff_add() == (0, (Float(3.0),))
    assert Float(-3.0).as_coeff_add() == (0, (Float(-3.0),))
    assert x.as_coeff_add() == (0, (x,))
    assert (x - 1).as_coeff_add() == (-1, (x,))
    assert (x + 1).as_coeff_add() == (1, (x,))
    assert (x + 2).as_coeff_add() == (2, (x,))
    assert (x + y).as_coeff_add(y) == (x, (y,))
    assert (3*x).as_coeff_add(y) == (3*x, ())
    # don't do expansion
    e = (x + y)**2
    assert e.as_coeff_add(y) == (0, (e,))


def test_as_coeff_mul():
    assert Integer(2).as_coeff_mul() == (2, ())
    assert Float(3.0).as_coeff_mul() == (1, (Float(3.0),))
    assert Float(-3.0).as_coeff_mul() == (-1, (Float(3.0),))
    assert Float(-3.0).as_coeff_mul(rational=False) == (-Float(3.0), ())
    assert x.as_coeff_mul() == (1, (x,))
    assert (-x).as_coeff_mul() == (-1, (x,))
    assert (2*x).as_coeff_mul() == (2, (x,))
    assert (x*y).as_coeff_mul(y) == (x, (y,))
    assert (3 + x).as_coeff_mul() == (1, (3 + x,))
    assert (3 + x).as_coeff_mul(y) == (3 + x, ())
    # don't do expansion
    e = exp(x + y)
    assert e.as_coeff_mul(y) == (1, (e,))
    e = 2**(x + y)
    assert e.as_coeff_mul(y) == (1, (e,))
    assert (1.1*x).as_coeff_mul(rational=False) == (1.1, (x,))
    assert (1.1*x).as_coeff_mul() == (1, (1.1, x))
    assert (-oo*x).as_coeff_mul(rational=True) == (-1, (oo, x))


def test_as_coeff_exponent():
    assert (3*x**4).as_coeff_exponent(x) == (3, 4)
    assert (2*x**3).as_coeff_exponent(x) == (2, 3)
    assert (4*x**2).as_coeff_exponent(x) == (4, 2)
    assert (6*x**1).as_coeff_exponent(x) == (6, 1)
    assert (3*x**0).as_coeff_exponent(x) == (3, 0)
    assert (2*x**0).as_coeff_exponent(x) == (2, 0)
    assert (1*x**0).as_coeff_exponent(x) == (1, 0)
    assert (0*x**0).as_coeff_exponent(x) == (0, 0)
    assert (-1*x**0).as_coeff_exponent(x) == (-1, 0)
    assert (-2*x**0).as_coeff_exponent(x) == (-2, 0)
    assert (2*x**3 + pi*x**3).as_coeff_exponent(x) == (2 + pi, 3)
    assert (x*log(2)/(2*x + pi*x)).as_coeff_exponent(x) == (log(2)/(2 + pi), 0)
    assert ((-x)**-pi).as_coeff_exponent(x) == ((-1)**-pi, -pi)
    # issue sympy/sympy#4784
    D = Derivative
    f = Function('f')
    fx = D(f(x), x)
    assert fx.as_coeff_exponent(f(x)) == (fx, 0)


def test_extractions():
    assert ((x*y)**3).extract_multiplicatively(1) == (x*y)**3
    assert ((x*y)**3).extract_multiplicatively(x**2 * y) == x*y**2
    assert ((x*y)**3).extract_multiplicatively(x**4 * y) is None
    assert (2*x).extract_multiplicatively(2) == x
    assert (2*x).extract_multiplicatively(3) is None
    assert (2*x).extract_multiplicatively(-1) is None
    assert (Rational(1, 2)*x).extract_multiplicatively(3) == x/6
    assert (sqrt(x)).extract_multiplicatively(x) is None
    assert (sqrt(x)).extract_multiplicatively(1/x) is None
    assert x.extract_multiplicatively(-x) is None
    assert oo.extract_multiplicatively(2) is oo
    assert oo.extract_multiplicatively(-1) is None
    assert (-oo).extract_multiplicatively(2) == -oo
    assert (-oo).extract_multiplicatively(-2) is oo
    assert (-oo).extract_multiplicatively(0) is None

    assert ((x*y)**3).extract_additively(1) is None
    assert (x + 1).extract_additively(x) == 1
    assert (x + 1).extract_additively(2*x) is None
    assert (x + 1).extract_additively(-x) is None
    assert (-x + 1).extract_additively(2*x) is None
    assert (2*x + 3).extract_additively(x) == x + 3
    assert (2*x + 3).extract_additively(2) == 2*x + 1
    assert (2*x + 3).extract_additively(3) == 2*x
    assert (2*x + 3).extract_additively(-2) is None
    assert (2*x + 3).extract_additively(3*x) is None
    assert (2*x + 3).extract_additively(2*x) == 3
    assert x.extract_additively(0) == x
    assert Integer(2).extract_additively(x) is None
    assert Float(2.).extract_additively(2) == 0
    assert (2*x + 3).extract_additively(x + 1) == x + 2
    assert (2*x + 3).extract_additively(y + 1) is None
    assert (2*x - 3).extract_additively(x + 1) is None
    assert (2*x - 3).extract_additively(y + z) is None
    assert ((a + 1)*x*4 + y).extract_additively(x).expand() == \
        4*a*x + 3*x + y
    assert ((a + 1)*x*4 + 3*y).extract_additively(x + 2*y).expand() == \
        4*a*x + 3*x + y
    assert (y*(x + 1)).extract_additively(x + 1) is None
    assert ((y + 1)*(x + 1) + 3).extract_additively(x + 1) == \
        y*(x + 1) + 3
    assert ((x + y)*(x + 1) + x + y + 3).extract_additively(x + y) == \
        x*(x + y) + 3
    assert (x + y + 2*((x + y)*(x + 1)) + 3).extract_additively((x + y)*(x + 1)) == \
        x + y + (x + 1)*(x + y) + 3
    assert ((y + 1)*(x + 2*y + 1) + 3).extract_additively(y + 1) == \
        (x + 2*y)*(y + 1) + 3

    n = Symbol('n', integer=True)
    assert (Integer(-3)).could_extract_minus_sign() is True
    assert (-n*x + x).could_extract_minus_sign() != \
        (n*x - x).could_extract_minus_sign()
    assert (x - y).could_extract_minus_sign() != \
        (-x + y).could_extract_minus_sign()
    assert (1 - x - y).could_extract_minus_sign() is True
    assert (1 - x + y).could_extract_minus_sign() is False
    assert ((-x - x*y)/y).could_extract_minus_sign() is True
    assert (-(x + x*y)/y).could_extract_minus_sign() is True
    assert ((x + x*y)/(-y)).could_extract_minus_sign() is True
    assert ((x + x*y)/y).could_extract_minus_sign() is False
    assert (x*(-x - x**3)).could_extract_minus_sign() is True
    assert ((-x - y)/(x + y)).could_extract_minus_sign() is True
    # The results of each of these will vary on different machines, e.g.
    # the first one might be False and the other (then) is true or vice versa,
    # so both are included.
    assert ((-x - y)/(x - y)).could_extract_minus_sign() is False or \
           ((-x - y)/(y - x)).could_extract_minus_sign() is False
    assert (x - y).could_extract_minus_sign() is False
    assert (-x + y).could_extract_minus_sign() is True

    # issue sympy/sympy#5843
    e = 1 + x
    assert (2*e).extract_multiplicatively(e) == 2
    assert (4*e).extract_multiplicatively(2*e) == 2
    assert ((3*e)*(2*e)).extract_multiplicatively(e) == 6*e


def test_nan_extractions():
    for r in (1, 0, I, nan):
        assert nan.extract_additively(r) is None
        assert nan.extract_multiplicatively(r) is None


def test_coeff():
    assert (x + 1).coeff(x + 1) == 1
    assert (3*x).coeff(0) == 0
    assert (z*(1 + x)*x**2).coeff(1 + x) == z*x**2
    assert (1 + 2*x*x**(1 + x)).coeff(x*x**(1 + x)) == 2
    assert (1 + 2*x**(y + z)).coeff(x**(y + z)) == 2
    assert (3 + 2*x + 4*x**2).coeff(1) == 0
    assert (3 + 2*x + 4*x**2).coeff(-1) == 0
    assert (3 + 2*x + 4*x**2).coeff(x) == 2
    assert (3 + 2*x + 4*x**2).coeff(x**2) == 4
    assert (3 + 2*x + 4*x**2).coeff(x**3) == 0

    assert (-x/8 + x*y).coeff(x) == -Rational(1, 8) + y
    assert (-x/8 + x*y).coeff(-x) == Rational(1, 8)
    assert (4*x).coeff(2*x) == 0
    assert (2*x).coeff(2*x) == 1
    assert (-oo*x).coeff(x*oo) == -1
    assert (10*x).coeff(x, 0) == 0
    assert (10*x).coeff(10*x, 0) == 0

    n1, n2 = symbols('n1 n2', commutative=False)
    assert (n1*n2).coeff(n1) == 1
    assert (n1*n2).coeff(n2) == n1
    assert (n1*n2 + x*n1).coeff(n1) == 1  # 1*n1*(n2+x)
    assert (n2*n1 + x*n1).coeff(n1) == n2 + x
    assert (n2*n1 + x*n1**2).coeff(n1) == n2
    assert (n1**x).coeff(n1) == 0
    assert (n1*n2 + n2*n1).coeff(n1) == 0
    assert (2*(n1 + n2)*n2).coeff(n1 + n2, right=1) == n2
    assert (2*(n1 + n2)*n2).coeff(n1 + n2, right=0) == 2
    assert (3*n1*n2 + 3*x*n1).coeff(n1) == 3

    f = Function('f')
    assert (2*f(x) + 3*f(x).diff(x)).coeff(f(x)) == 2

    expr = z*(x + y)**2
    expr2 = z*(x + y)**2 + z*(2*x + 2*y)**2
    assert expr.coeff(z) == (x + y)**2
    assert expr.coeff(x + y) == 0
    assert expr2.coeff(z) == (x + y)**2 + (2*x + 2*y)**2

    assert (x + y + 3*z).coeff(1) == x + y
    assert (-x + 2*y).coeff(-1) == x
    assert (x - 2*y).coeff(-1) == 2*y
    assert (3 + 2*x + 4*x**2).coeff(1) == 0
    assert (-x - 2*y).coeff(2) == -y
    assert (x + sqrt(2)*x).coeff(sqrt(2)) == x
    assert (3 + 2*x + 4*x**2).coeff(x) == 2
    assert (3 + 2*x + 4*x**2).coeff(x**2) == 4
    assert (3 + 2*x + 4*x**2).coeff(x**3) == 0
    assert (z*(x + y)**2).coeff((x + y)**2) == z
    assert (z*(x + y)**2).coeff(x + y) == 0
    assert (2 + 2*x + (x + 1)*y).coeff(x + 1) == y

    assert (x + 2*y + 3).coeff(1) == x
    assert (x + 2*y + 3).coeff(x, 0) == 2*y + 3
    assert (x**2 + 2*y + 3*x).coeff(x**2, 0) == 2*y + 3*x
    assert x.coeff(0, 0) == 0
    assert x.coeff(x, 0) == 0

    n, m, o = symbols('n m o', commutative=False)
    assert n.coeff(n) == 1
    assert y.coeff(n) == 0
    assert (3*n).coeff(n) == 3
    assert (2 + n).coeff(x*m) == 0
    assert (2*x*n*m).coeff(x) == 2*n*m
    assert (2 + n).coeff(x*m*n + y) == 0
    assert (2*x*n*m).coeff(3*n) == 0
    assert (n*m + m*n*m).coeff(n) == 1 + m
    assert (n*m + m*n*m).coeff(n, right=True) == m  # = (1 + m)*n*m
    assert (n*m + m*n).coeff(n) == 0
    assert (n*m + o*m*n).coeff(m*n) == o
    assert (n*m + o*m*n).coeff(m*n, right=1) == 1
    assert (n*m + n*m*n).coeff(n*m, right=1) == 1 + n  # = n*m*(n + 1)

    assert (x*y).coeff(z, 0) == x*y


def test_coeff2():
    psi = Function('psi')
    g = 1/r**2 * (2*r*psi(r).diff((r, 1)) + r**2 * psi(r).diff((r, 2)))
    g = g.expand()
    assert g.coeff(psi(r).diff(r)) == 2/r


def test_coeff2_0():
    psi = Function('psi')
    g = 1/r**2 * (2*r*psi(r).diff((r, 1)) + r**2 * psi(r).diff((r, 2)))
    g = g.expand()

    assert g.coeff(psi(r).diff((r, 2))) == 1


def test_coeff_expand():
    expr = z*(x + y)**2
    expr2 = z*(x + y)**2 + z*(2*x + 2*y)**2
    assert expr.coeff(z) == (x + y)**2
    assert expr2.coeff(z) == (x + y)**2 + (2*x + 2*y)**2


def test_integrate():
    assert x.integrate(x) == x**2/2
    assert x.integrate((x, 0, 1)) == Rational(1, 2)


def test_as_base_exp():
    assert x.as_base_exp() == (x, 1)
    assert (x*y*z).as_base_exp() == (x*y*z, 1)
    assert (x + y + z).as_base_exp() == (x + y + z, 1)
    assert ((x + y)**z).as_base_exp() == (x + y, z)


def test_sympyissue_4963():
    assert hasattr(Mul(x, y), 'is_commutative')
    assert hasattr(Mul(x, y, evaluate=False), 'is_commutative')
    assert hasattr(Pow(x, y), 'is_commutative')
    assert hasattr(Pow(x, y, evaluate=False), 'is_commutative')
    expr = Mul(Pow(2, 2, evaluate=False), 3, evaluate=False) + 1
    assert hasattr(expr, 'is_commutative')


def test_action_verbs():
    assert nsimplify(1/(exp(3*pi*x/5) + 1)) == \
        (1/(exp(3*pi*x/5) + 1)).nsimplify()
    assert ratsimp(1/x + 1/y) == (1/x + 1/y).ratsimp()
    assert trigsimp(log(x), deep=True) == (log(x)).trigsimp(deep=True)
    assert radsimp(1/(2 + sqrt(2))) == (1/(2 + sqrt(2))).radsimp()
    assert radsimp(1/(a + b*sqrt(c)), symbolic=False) == \
        (1/(a + b*sqrt(c))).radsimp(symbolic=False)
    assert powsimp(x**y*x**z*y**z, combine='all') == \
        (x**y*x**z*y**z).powsimp(combine='all')
    assert (x**t*y**t).powsimp(force=True) == (x*y)**t
    assert simplify(x**y*x**z*y**z) == (x**y*x**z*y**z).simplify()
    assert together(1/x + 1/y) == (1/x + 1/y).together()
    assert collect(a*x**2 + b*x**2 + a*x - b*x + c, x) == \
        (a*x**2 + b*x**2 + a*x - b*x + c).collect(x)
    assert apart(y/(y + 2)/(y + 1), y) == (y/(y + 2)/(y + 1)).apart(y)
    assert combsimp(y/(x + 2)/(x + 1)) == (y/(x + 2)/(x + 1)).combsimp()
    assert factor(x**2 + 5*x + 6) == (x**2 + 5*x + 6).factor()
    assert cancel((x**2 + 5*x + 6)/(x + 2)) == ((x**2 + 5*x + 6)/(x + 2)).cancel()


def test_as_powers_dict():
    assert x.as_powers_dict() == {x: 1}
    assert (x**y*z).as_powers_dict() == {x: y, z: 1}
    assert Mul(2, 2, evaluate=False).as_powers_dict() == {2: 2}
    assert (x*y).as_powers_dict()[z] == 0
    assert (x + y).as_powers_dict()[z] == 0


def test_as_coefficients_dict():
    check = [Integer(1), x, y, x*y, 1]
    assert [Add(3*x, 2*x, y, 3).as_coefficients_dict()[i] for i in check] == \
        [3, 5, 1, 0, 3]
    assert [Add(3*x, 2*x, y, 3, evaluate=False).as_coefficients_dict()[i]
            for i in check] == [3, 5, 1, 0, 3]
    assert [(3*x*y).as_coefficients_dict()[i] for i in check] == \
        [0, 0, 0, 3, 0]
    assert (3.0*x*y).as_coefficients_dict()[3.0*x*y] == 1


def test_args_cnc():
    A = symbols('A', commutative=False)
    assert (x + A).args_cnc() == \
        [[], [x + A]]
    assert (x + a).args_cnc() == \
        [[a + x], []]
    assert (x*a).args_cnc() == \
        [[a, x], []]
    assert (x*y*A*(A + 1)).args_cnc(cset=True) == \
        [{x, y}, [A, 1 + A]]
    assert Mul(x, x, evaluate=False).args_cnc(cset=True, warn=False) == \
        [{x}, []]
    assert Mul(x, x**2, evaluate=False).args_cnc(cset=True, warn=False) == \
        [{x, x**2}, []]
    pytest.raises(ValueError, lambda: Mul(x, x, evaluate=False).args_cnc(cset=True))
    assert Mul(x, y, x, evaluate=False).args_cnc() == \
        [[x, y, x], []]
    # always split -1 from leading number
    assert (-1.*x).args_cnc() == [[-1, 1.0, x], []]


def test_new_rawargs():
    n = Symbol('n', commutative=False)
    a = x + n
    assert a.is_commutative is False
    assert a._new_rawargs(x).is_commutative
    assert a._new_rawargs(x, y).is_commutative
    assert a._new_rawargs(x, n).is_commutative is False
    assert a._new_rawargs(x, y, n).is_commutative is False
    m = x*n
    assert m.is_commutative is False
    assert m._new_rawargs(x).is_commutative
    assert m._new_rawargs(n).is_commutative is False
    assert m._new_rawargs(x, y).is_commutative
    assert m._new_rawargs(x, n).is_commutative is False
    assert m._new_rawargs(x, y, n).is_commutative is False

    assert m._new_rawargs(x, n, reeval=False).is_commutative is False
    assert m._new_rawargs(Integer(1)) is Integer(1)


def test_sympyissue_5226():
    assert Add(evaluate=False) == 0
    assert Mul(evaluate=False) == 1
    assert Mul(x + y, evaluate=False).is_Add


def test_free_symbols():
    # free_symbols should return the free symbols of an object
    assert Integer(1).free_symbols == set()
    assert x.free_symbols == {x}
    assert Integral(x, (x, 1, y)).free_symbols == {y}
    assert (-Integral(x, (x, 1, y))).free_symbols == {y}


def test_sympyissue_5300():
    x = Symbol('x', commutative=False)
    assert x*sqrt(2)/sqrt(6) == x*sqrt(3)/3


def test_as_coeff_Mul():
    assert Integer(0).as_coeff_Mul() == (1, 0)
    assert Integer(3).as_coeff_Mul() == (3, 1)
    assert Rational(3, 4).as_coeff_Mul() == (Rational(3, 4), 1)
    assert Float(5.0).as_coeff_Mul() == (Float(5.0), 1)

    assert (Integer(3)*x).as_coeff_Mul() == (3, x)
    assert (Rational(3, 4)*x).as_coeff_Mul() == (Rational(3, 4), x)
    assert (Float(5.0)*x).as_coeff_Mul() == (Float(5.0), x)

    assert (Integer(3)*x*y).as_coeff_Mul() == (3, x*y)
    assert (Rational(3, 4)*x*y).as_coeff_Mul() == (Rational(3, 4), x*y)
    assert (Float(5.0)*x*y).as_coeff_Mul() == (Float(5.0), x*y)

    assert x.as_coeff_Mul() == (1, x)
    assert (x*y).as_coeff_Mul() == (1, x*y)
    assert (-oo*x).as_coeff_Mul(rational=True) == (-1, oo*x)


def test_as_coeff_Add():
    assert Integer(3).as_coeff_Add() == (3, 0)
    assert Rational(3, 4).as_coeff_Add() == (Rational(3, 4), 0)
    assert Float(5.0).as_coeff_Add() == (Float(5.0), 0)

    assert (Integer(3) + x).as_coeff_Add() == (3, x)
    assert (Rational(3, 4) + x).as_coeff_Add() == (Rational(3, 4), x)
    assert (Float(5.0) + x).as_coeff_Add() == (Float(5.0), x)

    assert (Integer(3) + x + y).as_coeff_Add() == (3, x + y)
    assert (Rational(3, 4) + x + y).as_coeff_Add() == (Rational(3, 4), x + y)
    assert (Float(5.0) + x + y).as_coeff_Add() == (Float(5.0), x + y)

    assert x.as_coeff_Add() == (0, x)
    assert (x*y).as_coeff_Add() == (0, x*y)


def test_expr_sorting():
    f, g = symbols('f,g', cls=Function)

    exprs = [1/x**2, 1/x, sqrt(sqrt(x)), sqrt(x), x, sqrt(x)**3, x**2]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [x, 2*x, 2*x**2, 2*x**3, x**n, 2*x**n, sin(x), sin(x)**n,
             sin(x**2), cos(x), cos(x**2), tan(x)]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [x + 1, x**2 + x + 1, x**3 + x**2 + x + 1]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [Integer(4), x - 3*I/2, x + 3*I/2, x - 4*I + 1, x + 4*I + 1]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [f(1), f(2), f(3), f(1, 2, 3), g(1), g(2), g(3), g(1, 2, 3)]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [exp(x), f(x), g(x), sin(x), cos(x), factorial(x)]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [Tuple(x, y), Tuple(x, z), Tuple(x, y, z)]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [[3], [1, 2]]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [[1, 2], [2, 3]]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [[1, 2], [1, 2, 3]]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [{x: -y}, {x: y}]
    assert sorted(exprs, key=default_sort_key) == exprs

    exprs = [{1}, {1, 2}]
    assert sorted(exprs, key=default_sort_key) == exprs

    a, b = exprs = [Dummy('x'), Dummy('x')]
    assert sorted([b, a], key=default_sort_key) == exprs


def test_as_ordered_factors():
    f, g = symbols('f,g', cls=Function)

    assert x.as_ordered_factors() == [x]
    assert (2*x*x**n*sin(x)*cos(x)).as_ordered_factors() \
        == [2, x, x**n, sin(x), cos(x)]

    args = [f(1), f(2), f(3), f(1, 2, 3), g(1), g(2), g(3), g(1, 2, 3)]
    expr = Mul(*args)

    assert expr.as_ordered_factors() == args

    A, B = symbols('A,B', commutative=False)

    assert (A*B).as_ordered_factors() == [A, B]
    assert (B*A).as_ordered_factors() == [B, A]


def test_as_ordered_terms():
    f, g = symbols('f,g', cls=Function)

    assert x.as_ordered_terms() == [x]
    assert (sin(x)**2*cos(x) + sin(x)*cos(x)**2 + 1).as_ordered_terms() \
        == [sin(x)**2*cos(x), sin(x)*cos(x)**2, 1]

    args = [f(1), f(2), f(3), f(1, 2, 3), g(1), g(2), g(3), g(1, 2, 3)]
    expr = Add(*args)

    assert expr.as_ordered_terms() == args

    assert (1 + 4*sqrt(3)*pi*x).as_ordered_terms() == [4*pi*x*sqrt(3), 1]

    assert (+2 + 3*I).as_ordered_terms() == [2, 3*I]
    assert (-2 + 3*I).as_ordered_terms() == [-2, 3*I]
    assert (+2 - 3*I).as_ordered_terms() == [2, -3*I]
    assert (-2 - 3*I).as_ordered_terms() == [-2, -3*I]

    assert (+4 + 3*I).as_ordered_terms() == [4, 3*I]
    assert (-4 + 3*I).as_ordered_terms() == [-4, 3*I]
    assert (+4 - 3*I).as_ordered_terms() == [4, -3*I]
    assert (-4 - 3*I).as_ordered_terms() == [-4, -3*I]

    assert (x*y).as_ordered_terms(data=True) == ([(x*y, ((1.0, 0.0),
                                                         (1, 1), ()))], [x, y])

    f = x**2*y**2 + x*y**4 + y + 2

    assert f.as_ordered_terms(order='lex') == [x**2*y**2, x*y**4, y, 2]
    assert f.as_ordered_terms(order='grlex') == [x*y**4, x**2*y**2, y, 2]
    assert f.as_ordered_terms(order='rev-lex') == [2, y, x*y**4, x**2*y**2]
    assert f.as_ordered_terms(order='rev-grlex') == [2, y, x**2*y**2, x*y**4]


def test__eval_interval():
    # issue sympy/sympy#4199
    # first subs and limit gives nan
    a = x/y
    assert a._eval_interval(x, 0, oo)._eval_interval(y, oo, 0) is nan
    # second subs and limit gives nan
    assert a._eval_interval(x, 0, oo)._eval_interval(y, 0, oo) is nan
    # difference gives nan
    a = x - y
    assert a._eval_interval(x, 1, oo)._eval_interval(y, oo, 1) is nan
    pytest.raises(ValueError, lambda: x._eval_interval(x, None, None))
    a = -y*Heaviside(x - y)
    assert a._eval_interval(x, -oo, oo) == -y
    assert a._eval_interval(x, oo, -oo) == y

    # Test that limit is used when zoo is returned
    assert Si(1/x)._eval_interval(x, 0, 1) == -pi/2 + Si(1)


def test_primitive():
    assert (3*(x + 1)**2).primitive() == (3, (x + 1)**2)
    assert (6*x + 2).primitive() == (2, 3*x + 1)
    assert (x/2 + 3).primitive() == (Rational(1, 2), x + 6)
    eq = (6*x + 2)*(x/2 + 3)
    assert eq.primitive()[0] == 1
    eq = (2 + 2*x)**2
    assert eq.primitive()[0] == 1
    assert (4.0*x).primitive() == (1, 4.0*x)
    assert (4.0*x + y/2).primitive() == (Rational(1, 2), 8.0*x + y)
    assert (-2*x).primitive() == (2, -x)
    assert Add(5*z/7, 0.5*x, 3*y/2, evaluate=False).primitive() == \
        (Rational(1, 14), 7.0*x + 21*y + 10*z)
    for i in (oo, -oo, zoo):
        assert (i + x/3).primitive() == \
            (Rational(1, 3), i + x)
    assert (oo + 2*x/3 + 4*y/7).primitive() == \
        (Rational(1, 21), 14*x + 12*y + oo)
    assert Integer(0).primitive() == (1, 0)


def test_is_constant():
    assert Sum(x, (x, 1, 10)).is_constant() is True
    assert Sum(x, (x, 1, n)).is_constant() is False
    assert Sum(x, (x, 1, n)).is_constant(y) is True
    assert Sum(x, (x, 1, n)).is_constant(n) is False
    assert Sum(x, (x, 1, n)).is_constant(x) is True
    eq = a*cos(x)**2 + a*sin(x)**2 - a
    assert eq.is_constant() is True
    assert eq.subs({x: pi, a: 2}) == eq.subs({x: pi, a: 3}) == 0
    assert x.is_constant() is False
    assert x.is_constant(y) is True

    assert checksol(x, {x: Sum(x, (x, 1, n))}) is False
    assert checksol(x, {x: Sum(x, (x, 1, n))}) is False
    f = Function('f')
    assert checksol(x, {x: f(x)}) is False

    p = symbols('p', positive=True)
    assert Pow(x, 0, evaluate=False).is_constant() is True  # == 1
    assert Pow(0, x, evaluate=False).is_constant() is False  # == 0 or 1
    assert Pow(0, p, evaluate=False).is_constant() is True  # == 1
    assert (2**x).is_constant() is False
    assert Pow(2, 3, evaluate=False).is_constant() is True

    z1, z2 = symbols('z1 z2', zero=True)
    assert (z1 + 2*z2).is_constant() is True

    e = factorial(x) % x
    assert e.subs({x: x - 1}).is_constant() is False

    e = a*sin(x)**2 + a*cos(x)**2
    assert e.is_constant(x) is True
    assert e.is_constant(a) is False

    assert Integer(2).is_constant() is True

    for _ in range(5):
        assert Piecewise((x, (x < 0)), (0, True)).is_constant() is not True
        assert Piecewise((1, (x < 0)), (0, True)).is_constant() is not True


def test_equals():
    assert (-3 - sqrt(5) + (-sqrt(10)/2 - sqrt(2)/2)**2).equals(0)
    assert (x**2 - 1).equals((x + 1)*(x - 1))
    assert (cos(x)**2 + sin(x)**2).equals(1)
    assert (a*cos(x)**2 + a*sin(x)**2).equals(a)
    r = sqrt(2)
    assert (-1/(r + r*x) + 1/r/(1 + x)).equals(0)
    assert factorial(x + 1).equals((x + 1)*factorial(x))
    assert sqrt(3).equals(2*sqrt(3)) is False
    assert (sqrt(5)*sqrt(3)).equals(sqrt(3)) is False
    assert (sqrt(5) + sqrt(3)).equals(0) is False
    assert (sqrt(5) + pi).equals(0) is False
    eq = -(-1)**Rational(3, 4)*root(6, 4) + root(-6, 4)*I
    if eq != 0:  # if canonicalization makes this zero, skip the test
        assert eq.equals(0)
    assert sqrt(x).equals(0) is False

    # from integrate(x*sqrt(1 + 2*x), x);
    # diff is zero only when assumptions allow
    i = 2*sqrt(2)*x**Rational(5, 2)*(1 + 1/(2*x))**Rational(5, 2)/5 + \
        2*sqrt(2)*x**Rational(3, 2)*(1 + 1/(2*x))**Rational(5, 2)/(-6 - 3/x)
    ans = sqrt(2*x + 1)*(6*x**2 + x - 1)/15
    diff = i - ans
    assert diff.equals(0) is False
    assert diff.subs({x: Rational(-1, 4)}) == 7*sqrt(2)/120
    # there are regions for x for which the expression is True, for
    # example, when x < -1/2 or x > 0 the expression is zero
    p = Symbol('p', positive=True)
    assert diff.subs({x: p}).equals(0) is True
    assert diff.subs({x: -1}).equals(0) is True

    # prove via minimal_polynomial or self-consistency
    eq = sqrt(1 + sqrt(3)) + sqrt(3 + 3*sqrt(3)) - sqrt(10 + 6*sqrt(3))
    assert eq.equals(0)
    q = cbrt(3) + 3
    p = cbrt(expand(q**3))
    assert (p - q).equals(0)

    # issue sympy/sympy#6829
    # eq = q*x + q/4 + x**4 + x**3 + 2*x**2 - Rational(1, 3)
    # z = eq.subs(solve(eq, x)[0])
    q, x0, x1, x2 = symbols('q, x:3')
    z = q*x2 + q/4 + x2**4 + x2**3 + 2*x2**2 - Rational(1, 3)
    z = z.subs(((x2, (-x1/2 - sqrt(x0 - Rational(13, 6) +
                                   (2*q - Rational(7, 4))/x1)/2 -
                      Rational(1, 4))),
                (x1, sqrt(-x0 - Rational(13, 12))),
                (x0, 2*cbrt(-(q - Rational(7, 8))**2/8 -
                            Rational(2197, 13824)))))
    z = expand_multinomial(z)
    assert z.equals(0)


def test_random():
    assert posify(x)[0]._random() is not None
    assert lucas(n)._random(2, -2, 0, -1, 1) is None

    # issue sympy/sympy#8662
    assert Piecewise((Max(x, y), z))._random() is None

    assert sqrt(2)._random(2) == Float('1.4141', dps=2)


def test_round():
    assert Float('0.1249999').round(2) == 0.12
    d20 = 12345678901234567890
    ans = Integer(d20).round(2)
    assert ans.is_Float
    assert ans == d20
    ans = Integer(d20).round(-2)
    assert ans.is_Float
    assert ans == 12345678901234567900
    assert Rational(1, 7).round(4) == 0.1429
    assert Float(.1349).round(2) == 0.13
    n = Integer(12345)
    ans = n.round()
    assert ans.is_Float
    assert ans == n
    ans = n.round(1)
    assert ans.is_Float
    assert ans == n
    ans = n.round(4)
    assert ans.is_Float
    assert ans == n
    assert n.round(-1) == 12350

    r = n.round(-4)
    assert r == 10000
    # in fact, it should equal many values since __eq__
    # compares at equal precision
    assert all(r == i for i in range(9984, 10049))

    assert n.round(-5) == 0

    assert (pi + sqrt(2)).round(2) == 4.56
    assert (10*(pi + sqrt(2))).round(-1) == 50
    pytest.raises(TypeError, lambda: round(x + 2, 2))
    assert Float(2.3).round(1) == 2.3
    e = Float(12.345).round(2)
    assert e == round(12.345, 2)
    assert type(e) is Float

    assert (Float(.3, 3) + 2*pi).round() == 7
    assert (Float(.3, 3) + 2*pi*100).round() == 629
    assert (Float(.03, 3) + 2*pi/100).round(5) == 0.09283
    assert (Float(.03, 3) + 2*pi/100).round(4) == 0.0928
    assert (pi + 2*E*I).round() == 3 + 5*I

    assert Integer(0).round() == 0

    a = (Add(1, Float('1.' + '9'*27), evaluate=0))
    assert a.round(10) == Float('3.0000000000')
    assert a.round(25) == Float('3.0000000000000000000000000')
    assert a.round(26) == Float('3.00000000000000000000000000')
    assert a.round(27) == Float('2.999999999999999999999999999')
    assert a.round(30) == Float('2.999999999999999999999999999')

    pytest.raises(TypeError, x.round)

    # exact magnitude of 10
    assert str(Integer(1).round()) == '1.'
    assert str(Integer(100).round()) == '100.'

    # applied to real and imaginary portions
    assert (2*pi + E*I).round() == 6 + 3*I
    assert (2*pi + I/10).round() == 6
    assert (pi/10 + 2*I).round() == 2*I
    # the lhs re and im parts are Float with dps of 2
    # and those on the right have dps of 15 so they won't compare
    # equal unless we use string or compare components (which will
    # then coerce the floats to the same precision) or re-create
    # the floats
    assert str((pi/10 + E*I).round(2)) == '0.31 + 2.72*I'
    assert (pi/10 + E*I).round(2).as_real_imag() == (0.31, 2.72)
    assert (pi/10 + E*I).round(2) == Float(0.31, 2) + I*Float(2.72, 3)

    # issue sympy/sympy#6914
    assert (I**(I + 3)).round(3) == Float('-0.208')*I

    # issue sympy/sympy#8720
    assert Float(-123.6).round() == -124.
    assert Float(-1.5).round() == -2.
    assert Float(-100.5).round() == -101.
    assert (Float(-1.5) - Float(10.5)*I).round() == -2.0 - 11.0*I

    # issue sympy/sympy#7961
    assert str(Float(0.006).round(2)) == '0.01'
    assert str(Float(0.00106).round(4)) == '0.0011'

    # issue sympy/sympy#8147
    assert nan.round() == nan
    assert (+oo).round() == +oo
    assert (-oo).round() == -oo
    assert zoo.round() == zoo


def test_round_exception_nostr():
    # Don't use the string form of the expression in the round exception, as
    # it's too slow
    s = Symbol('bad')
    try:
        s.round()
    except TypeError as e:
        assert 'bad' not in str(e)
    else:
        # Did not raise
        raise AssertionError('Did not raise')


def test_extract_branch_factor():
    assert exp_polar(2.0*I*pi).extract_branch_factor() == (1, 1)
    assert exp_polar(3*pi*I + x).extract_branch_factor() == (exp_polar(x + I*pi), 1)
    e = y*exp_polar(-5*pi*I)*exp_polar(3*pi*I + 2*pi*x)
    assert e.extract_branch_factor() == (y*exp_polar(2*pi*x), -1)


def test_identity_removal():
    assert Add.make_args(x + 0) == (x,)
    assert Mul.make_args(x*1) == (x,)


def test_float_0():
    assert Float(0.0) + 1 == Float(1.0)


@pytest.mark.xfail
def test_float_0_fail():
    assert Float(0.0)*x == Float(0.0)
    assert (x + Float(0.0)).is_Add


def test_invert():
    assert (x**2 - 1).invert(2*x - 1) == Rational(-4, 3)


def test_sympyissue_6325():
    ans = (b**2 + z**2 - (b*(a + b*t) + z*(c + t*z))**2/(
        (a + b*t)**2 + (c + t*z)**2))/sqrt((a + b*t)**2 + (c + t*z)**2)
    e = sqrt((a + b*t)**2 + (c + z*t)**2)
    assert diff(e, (t, 2)) == ans
    assert e.diff((t, 2)) == ans
    assert diff(e, (t, 2), simplify=False) != ans


def test_sympyissue_7426():
    f1 = a % c
    f2 = x % z
    assert f1.equals(f2) is False


def test_sympyissue_11122():
    p = Symbol('p', positive=False)
    assert (p > 0) is false


def test_pow_rewrite():
    assert (2**x).rewrite(sin) == 2**x
    assert (2**x).rewrite(tanh) == 2**x


@pytest.mark.xfail
@pytest.mark.timeout(20)
def test_sympyissue_13645():
    # NB: this does work if r and M are real
    r, r_m = symbols('r r_m', positive=True)
    th = symbols('th', extended_real=True)
    a, M = symbols('a M', extended_real=True)
    kappa, gamma = symbols('kappa gamma', extended_real=True)

    l = sqrt(M/r_m**3)*(r_m**4 + r_m**2*a**2 - 2*M*r_m*a**2 -
                        a*sqrt(M*r_m)*(r_m**2-a**2))
    l /= r_m**2 - 3*M*r_m + 2*a*sqrt(M*r_m)

    Delta = r**2 - 2*M*r + a**2
    Sigma = r**2 + a**2*cos(th)**2
    A = (r**2 + a**2)**2 - Delta*a**2*sin(th)**2

    tmp3 = sqrt(1 + l**2*Sigma**2*Delta/(sin(th))**2)
    ln_h = log((1 + tmp3) / (Sigma*Delta/A)) - tmp3
    hm1 = exp(ln_h) - 1

    # not hangs
    assert (hm1*(gamma-1)/(kappa*gamma))**(1/(gamma - 1)) is not None


def test_sympyissue_21334():
    e = exp(-x**2/(x + 1) + x) - exp(x/(x + 1)) + O(y)

    assert e.as_leading_term(y) == 0


def test_sympyissue_22583():
    f = Function('f')
    assert (1/f(x) + 1).is_polynomial(f(x)) is False
