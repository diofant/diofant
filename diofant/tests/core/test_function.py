import inspect
from random import random

import pytest

from diofant import (Derivative, Dummy, E, Eq, Expr, FiniteSet, Float,
                     Function, I, Integer, Lambda, O, PoleError, Rational,
                     RootOf, S, Subs, Sum, Symbol, Tuple, acos, cbrt, ceiling,
                     cos, diff, exp, expand, expint, floor, im, log, loggamma,
                     nan, nfloat, oo, pi, polygamma, re, sin, solve, sqrt,
                     symbols, zoo)
from diofant.abc import a, b, t, w, x, y, z
from diofant.core.basic import _aresame
from diofant.core.cache import clear_cache
from diofant.core.function import (ArgumentIndexError, UndefinedFunction,
                                   _mexpand)


__all__ = ()

f, g, h = symbols('f g h', cls=Function)


def test_f_expand_complex():
    x = Symbol('x', extended_real=True)

    assert f(x).expand(complex=True) == I*im(f(x)) + re(f(x))
    assert exp(x).expand(complex=True) == exp(x)
    assert exp(I*x).expand(complex=True) == cos(x) + I*sin(x)
    assert exp(z).expand(complex=True) == cos(im(z))*exp(re(z)) + \
        I*sin(im(z))*exp(re(z))


def test_bug1():
    e = sqrt(-log(w))
    assert e.subs({log(w): -x}) == sqrt(x)

    e = sqrt(-5*log(w))
    assert e.subs({log(w): -x}) == sqrt(5*x)


def test_general_function():
    nu = Function('nu')

    e = nu()
    assert e.diff(x) == 0
    assert e.diff(y) == 0

    e = nu(x)
    edx = e.diff(x)
    edy = e.diff(y)
    edxdx = e.diff(x).diff(x)
    edxdy = e.diff(x).diff(y)
    assert e == nu(x)
    assert edx != nu(x)
    assert edx == diff(nu(x), x)
    assert edy == 0
    assert edxdx == diff(diff(nu(x), x), x)
    assert edxdy == 0


def test_derivative_subs_bug():
    e = diff(g(x), x)
    assert e.subs({g(x): +f(x)}) != e
    assert e.subs({g(x): +f(x)}) == Derivative(f(x), x)
    assert e.subs({g(x): -f(x)}) == Derivative(-f(x), x)

    assert e.subs({x: y}) == Derivative(g(y), y)


def test_derivative_subs_self_bug():
    d = diff(f(x), x)

    assert d.subs({d: y}) == y


def test_derivative_linearity():
    assert diff(-f(x), x) == -diff(f(x), x)
    assert diff(8*f(x), x) == 8*diff(f(x), x)
    assert diff(8*f(x), x) != 7*diff(f(x), x)
    assert diff(8*f(x)*x, x) == 8*f(x) + 8*x*diff(f(x), x)
    assert diff(8*f(x)*y*x, x) == 8*y*f(x) + 8*y*x*diff(f(x), x)


def test_derivative_evaluate():
    assert Derivative(sin(x), x) != diff(sin(x), x)
    assert Derivative(sin(x), x).doit() == diff(sin(x), x)

    assert Derivative(Derivative(f(x), x), x) == diff(f(x), x, x)
    assert Derivative(sin(x), (x, 0)) == sin(x)


def test_diff_symbols():
    assert diff(f(x, y, z), x, y, z) == Derivative(f(x, y, z), x, y, z)
    assert diff(f(x, y, z), x, x, x) == Derivative(f(x, y, z), x, x, x)
    assert diff(f(x, y, z), (x, 3)) == Derivative(f(x, y, z), (x, 3))

    # issue sympy/sympy#5028
    assert [diff(-z + x/y, sym) for sym in (z, x, y)] == [-1, 1/y, -x/y**2]
    assert diff(f(x, y, z), x, y, (z, 2)) == Derivative(f(x, y, z), x, y, z, z)
    assert diff(f(x, y, z), x, y, (z, 2), evaluate=False) == \
        Derivative(f(x, y, z), x, y, z, z)
    assert Derivative(f(x, y, z), x, y, z)._eval_derivative(z) == \
        Derivative(f(x, y, z), x, y, z, z)
    assert Derivative(Derivative(f(x, y, z), x), y)._eval_derivative(z) == \
        Derivative(f(x, y, z), x, y, z)


def test_Function():
    class MyFunc(Function):
        @classmethod
        def eval(cls, x):  # one arg
            return

    assert MyFunc.nargs == FiniteSet(1)
    assert MyFunc(x).nargs == FiniteSet(1)
    pytest.raises(TypeError, lambda: MyFunc(x, y).nargs)

    class MyFunc2(Function):
        @classmethod
        def eval(cls, *x):  # star args
            return

    assert MyFunc2.nargs == S.Naturals0
    assert MyFunc2(x).nargs == S.Naturals0


def test_nargs():
    f = Function('f')
    assert f.nargs == S.Naturals0
    assert f(1).nargs == S.Naturals0
    assert Function('f', nargs=2)(1, 2).nargs == FiniteSet(2)
    assert sin.nargs == FiniteSet(1)
    assert sin(2).nargs == FiniteSet(1)
    assert log.nargs == FiniteSet(1, 2)
    assert log(2).nargs == FiniteSet(1, 2)
    assert Function('f', nargs=2).nargs == FiniteSet(2)
    assert Function('f', nargs=0).nargs == FiniteSet(0)
    assert Function('f', nargs=(0, 1)).nargs == FiniteSet(0, 1)
    assert Function('f', nargs=None).nargs == S.Naturals0
    pytest.raises(ValueError, lambda: Function('f', nargs=()))

    p = Function('g', nargs=(1, 2))(1)
    assert p.args == (1,)
    assert p.nargs == FiniteSet(1, 2)

    pytest.raises(ValueError, lambda: sin(x, no_such_option=1))


def test_Lambda():
    e = Lambda(x, x**2)
    assert e(4) == 16
    assert e(x) == x**2
    assert e(y) == y**2

    assert Lambda(x, x**2) == Lambda(x, x**2)
    assert Lambda(x, x**2) == Lambda(y, y**2)
    assert Lambda(x, x**2) != Lambda(y, y**2 + 1)
    assert Lambda((x, y), x**y) == Lambda((y, x), y**x)
    assert Lambda((x, y), x**y) != Lambda((x, y), y**x)

    assert Lambda((x, y), x**y)(x, y) == x**y
    assert Lambda((x, y), x**y)(3, 3) == 3**3
    assert Lambda((x, y), x**y)(x, 3) == x**3
    assert Lambda((x, y), x**y)(3, y) == 3**y
    assert Lambda(x, f(x))(x) == f(x)
    assert Lambda(x, x**2)(e(x)) == x**4
    assert e(e(x)) == x**4

    assert Lambda((x, y), x + y).nargs == FiniteSet(2)

    p = x, y, z, t
    assert Lambda(p, t*(x + y + z))(*p) == t * (x + y + z)

    assert Lambda(x, 2*x) + Lambda(y, 2*y) == 2*Lambda(x, 2*x)
    assert Lambda(x, 2*x) not in [Lambda(x, x)]
    pytest.raises(TypeError, lambda: Lambda(1, x))
    assert Lambda(x, 1)(1) is Integer(1)

    assert (2*x).canonical_variables == {}
    assert Lambda(x, 2*x).canonical_variables == {x: Symbol('0_')}
    x_ = Symbol('x_')
    assert Lambda(x_, 2*x_).canonical_variables == {x_: Symbol('0__')}


def test_IdentityFunction():
    assert Lambda(x, x) is Lambda(y, y) is S.IdentityFunction
    assert Lambda(x, 2*x) is not S.IdentityFunction
    assert Lambda((x, y), x) is not S.IdentityFunction


def test_Lambda_symbols():
    assert Lambda(x, 2*x).free_symbols == set()
    assert Lambda(x, x*y).free_symbols == {y}


def test_Lambda_arguments():
    pytest.raises(TypeError, lambda: Lambda(x, 2*x)(x, y))
    pytest.raises(TypeError, lambda: Lambda((x, y), x + y)(x))


def test_Lambda_equality():
    assert Lambda(x, 2*x) == Lambda(y, 2*y)
    # although variables are casts as Dummies, the expressions
    # should still compare equal
    assert Lambda((x, y), 2*x) == Lambda((x, y), 2*x)
    assert Lambda(x, 2*x) != Lambda((x, y), 2*x)
    assert Lambda(x, 2*x) != 2*x


def test_Subs():
    assert Subs(x, (x, 0)) == Subs(y, (y, 0))
    assert Subs(x, (x, 0)).subs({x: 1}) == Subs(x, (x, 0))
    assert Subs(y, (x, 0)).subs({y: 1}) == Subs(1, (x, 0))
    assert Subs(f(x), (x, 0)).doit() == f(0)
    assert Subs(f(x**2), (x**2, 0)).doit() == f(0)
    assert Subs(f(x, y, z), (x, 0), (y, 1), (z, 1)) != \
        Subs(f(x, y, z), (x, 0), (y, 0), (z, 1))
    assert Subs(f(x, y), (x, 0), (y, 1), (z, 1)) == \
        Subs(f(x, y), (x, 0), (y, 1), (z, 2))
    assert Subs(f(x, y), (x, 0), (y, 1), (z, 1)) != \
        Subs(f(x, y) + z, (x, 0), (y, 1), (z, 0))
    assert Subs(f(x, y), (x, 0), (y, 1)).doit() == f(0, 1)
    assert Subs(Subs(f(x, y), (x, 0)), (y, 1)).doit() == f(0, 1)
    pytest.raises(ValueError, lambda: Subs(f(x, y), x))
    pytest.raises(ValueError, lambda: Subs(f(x, y), (x, 0), (x, 0), (y, 1)))

    assert len(Subs(f(x, y), (x, 0), (y, 1)).variables) == 2
    assert Subs(f(x, y), (x, 0), (y, 1)).point == Tuple(0, 1)

    assert Subs(f(x), (x, 0)) == Subs(f(y), (y, 0))
    assert Subs(f(x, y), (x, 0), (y, 1)) == Subs(f(x, y), (y, 1), (x, 0))
    assert Subs(f(x)*y, (x, 0), (y, 1)) == Subs(f(y)*x, (y, 0), (x, 1))
    assert Subs(f(x)*y, (x, 1), (y, 1)) == Subs(f(y)*x, (x, 1), (y, 1))

    assert Subs(f(x), (x, 0)).subs({x: 1}).doit() == f(0)
    assert Subs(f(x), (x, y)).subs({y: 0}) == Subs(f(x), (x, 0))
    assert Subs(y*f(x), (x, y)).subs({y: 2}) == Subs(2*f(x), (x, 2))
    assert (2 * Subs(f(x), (x, 0))).subs({Subs(f(x), (x, 0)): y}) == 2*y

    assert Subs(f(x), (x, 0)).free_symbols == set()
    assert Subs(f(x, y), (x, z)).free_symbols == {y, z}

    assert Subs(f(x).diff(x), (x, 0)).doit(), Subs(f(x).diff(x), (x, 0))
    assert Subs(1 + f(x).diff(x), (x, 0)).doit(), 1 + Subs(f(x).diff(x), (x, 0))
    assert Subs(y*f(x, y).diff(x), (x, 0), (y, 2)).doit() == \
        2*Subs(Derivative(f(x, 2), x), (x, 0))
    assert Subs(y**2*f(x), (x, 0)).diff(y) == 2*y*f(0)

    e = Subs(y**2*f(x), (x, y))
    assert e.diff(y) == e.doit().diff(y) == y**2*Derivative(f(y), y) + 2*y*f(y)

    assert Subs(f(x), (x, 0)) + Subs(f(x), (x, 0)) == 2*Subs(f(x), (x, 0))
    e1 = Subs(z*f(x), (x, 1))
    e2 = Subs(z*f(y), (y, 1))
    assert e1 + e2 == 2*e1
    assert hash(e1) == hash(e2)
    assert Subs(z*f(x + 1), (x, 1)) not in (e1, e2)
    assert Derivative(f(x), x).subs({x: g(x)}) == Derivative(f(g(x)), g(x))
    assert Derivative(f(x), x).subs({x: x + y}) == Subs(Derivative(f(x), x),
                                                        (x, x + y))
    assert Subs(f(x)*cos(y) + z, (x, 0), (y, pi/3)).evalf(2, strict=False) == \
        Subs(f(x)*cos(y) + z, (x, 0), (y, pi/3)).evalf(2, strict=False) == \
        z + Rational('1/2').evalf(2)*f(0)

    assert f(x).diff(x).subs({x: 0}).subs({x: y}) == f(x).diff(x).subs({x: 0})
    assert (x*f(x).diff(x).subs({x: 0})).subs({x: y}) == y*f(x).diff(x).subs({x: 0})


def test_expand_function():
    assert expand(x + y) == x + y
    assert expand(x + y, complex=True) == I*im(x) + I*im(y) + re(x) + re(y)
    assert expand((x + y)**11, modulus=11) == x**11 + y**11


def test_function_comparable():
    assert sin(x).is_comparable is False
    assert cos(x).is_comparable is False

    assert sin(Float('0.1')).is_comparable is True
    assert cos(Float('0.1')).is_comparable is True

    assert sin(E).is_comparable is True
    assert cos(E).is_comparable is True

    assert sin(Rational(1, 3)).is_comparable is True
    assert cos(Rational(1, 3)).is_comparable is True


def test_function_comparable_infinities():
    assert sin(+oo).is_comparable is False
    assert sin(-oo).is_comparable is False
    assert sin(zoo).is_comparable is False
    assert sin(nan).is_comparable is False


def test_deriv1():
    # These all requre derivatives evaluated at a point (issue sympy/sympy#4719) to work.
    # See issue sympy/sympy#4624
    assert f(2*x).diff(x) == 2*Subs(Derivative(f(x), x), (x, 2*x))
    assert (f(x)**3).diff(x) == 3*f(x)**2*f(x).diff(x)
    assert (f(2*x)**3).diff(x) == 6*f(2*x)**2*Subs(Derivative(f(x), x), (x, 2*x))

    assert f(2 + x).diff(x) == Subs(Derivative(f(x), x), (x, x + 2))
    assert f(2 + 3*x).diff(x) == 3*Subs(Derivative(f(x), x), (x, 3*x + 2))
    assert f(3*sin(x)).diff(x) == 3*cos(x)*Subs(Derivative(f(x), x), (x, 3*sin(x)))

    # See issue sympy/sympy#8510
    assert f(x, x + z).diff(x) == Subs(Derivative(f(y, x + z), y), (y, x)) \
        + Subs(Derivative(f(x, y), y), (y, x + z))
    assert f(x, x**2).diff(x) == Subs(Derivative(f(y, x**2), y), (y, x)) \
        + 2*x*Subs(Derivative(f(x, y), y), (y, x**2))


def test_deriv2():
    assert (x**3).diff(x) == 3*x**2
    assert (x**3).diff(x, evaluate=False) != 3*x**2
    assert (x**3).diff(x, evaluate=False) == Derivative(x**3, x)

    assert diff(x**3, x) == 3*x**2
    assert diff(x**3, x, evaluate=False) != 3*x**2
    assert diff(x**3, x, evaluate=False) == Derivative(x**3, x)


def test_func_deriv():
    assert f(x).diff(x) == Derivative(f(x), x)
    # issue sympy/sympy#4534
    assert f(x, y).diff(x, y) - f(x, y).diff(y, x) == 0
    assert Derivative(f(x, y), x, y).args[1:] == (x, y)
    assert Derivative(f(x, y), y, x).args[1:] == (y, x)
    assert (Derivative(f(x, y), x, y) - Derivative(f(x, y), y, x)).doit() == 0


def test_suppressed_evaluation():
    a = sin(0, evaluate=False)
    assert a != 0
    assert isinstance(a, sin)
    assert a.args == (0,)


def test_function_evalf():
    def eq(a, b, eps):
        return abs(a - b) < eps
    assert eq(sin(1).evalf(15), Float('0.841470984807897'), 1e-13)
    assert eq(
        sin(2).evalf(25), Float('0.9092974268256816953960199', 25), 1e-23)
    assert eq(sin(1 + I).evalf(
        15), Float('1.29845758141598') + Float('0.634963914784736')*I, 1e-13)
    assert eq(exp(1 + I).evalf(15), Float(
        '1.46869393991588') + Float('2.28735528717884239')*I, 1e-13)
    assert eq(exp(-0.5 + 1.5*I).evalf(15, strict=False), Float(
        '0.0429042815937374') + Float('0.605011292285002')*I, 1e-13)
    assert eq(log(pi + sqrt(2)*I).evalf(
        15), Float('1.23699044022052') + Float('0.422985442737893')*I, 1e-13)
    assert eq(cos(100).evalf(15), Float('0.86231887228768'), 1e-13)


def test_extensibility_eval():
    class MyFunc(Function):
        @classmethod
        def eval(cls, *args):
            return 0, 0, 0
    assert MyFunc(0) == (0, 0, 0)


def test_function_signature():
    assert str(inspect.signature(sin)) == '(arg)'
    assert str(inspect.signature(log)) == '(arg, base=None)'


def test_function_non_commutative():
    x = Symbol('x', commutative=False)
    assert f(x).is_commutative is False
    assert sin(x).is_commutative is False
    assert exp(x).is_commutative is False
    assert log(x).is_commutative is False
    assert f(x).is_complex is False
    assert sin(x).is_complex is False
    assert exp(x).is_complex is False
    assert log(x).is_complex is False


def test_function_complex():
    x = Symbol('x', complex=True)
    assert f(x).is_commutative is True
    assert sin(x).is_commutative is True
    assert exp(x).is_commutative is True
    assert log(x).is_commutative is True
    assert sin(x).is_complex is True
    assert exp(x).is_complex is True
    assert log(x).is_complex is None  # could be zero
    n = Symbol('n', complex=True, nonzero=True)
    z = Symbol('z', zero=True)
    assert log(n).is_complex is True
    assert log(z).is_complex is False


def test_function__eval_nseries():
    n = Symbol('n')

    assert sin(x)._eval_nseries(x, 2, None) == x + O(x**3)
    assert sin(x + 1)._eval_nseries(x, 2, None) == x*cos(1) + sin(1) + O(x**2)
    assert sin(pi*(1 - x))._eval_nseries(x, 2, None) == pi*x + O(x**3)
    assert acos(1 - x**2)._eval_nseries(x, 2, None) == sqrt(2)*x + O(x**2)
    assert polygamma(n, x + 1)._eval_nseries(x, 2, None) == \
        polygamma(n, 1) + polygamma(n + 1, 1)*x + O(x**2)
    pytest.raises(PoleError, lambda: sin(1/x)._eval_nseries(x, 2, None))
    assert acos(1 - x)._eval_nseries(x, 4, None) == sqrt(2)*sqrt(x) + \
        sqrt(2)*x**Rational(3, 2)/12 + O(x**2)
    assert acos(1 + x)._eval_nseries(x, 4, None) == sqrt(2)*I*sqrt(x) - \
        sqrt(2)*I*x**(3/2)/12 + O(x**2)  # XXX: wrong, branch cuts
    assert loggamma(1/x)._eval_nseries(x, 0, None) == \
        log(x)/2 - log(x)/x - 1/x + O(1, x)
    assert loggamma(log(1/x)).series(x, n=1, logx=y) == loggamma(-y)

    # issue sympy/sympy#6725:
    assert expint(Rational(3, 2), -x)._eval_nseries(x, 8, None) == \
        2 - 2*I*sqrt(pi)*sqrt(x) - 2*x - x**2/3 - x**3/15 - x**4/84 + O(x**5)
    assert sin(sqrt(x))._eval_nseries(x, 6, None) == \
        sqrt(x) - x**Rational(3, 2)/6 + x**Rational(5, 2)/120 + O(x**Rational(7, 2))


def test_doit():
    n = Symbol('n', integer=True)
    f = Sum(2 * n * x, (n, 1, 3))
    d = Derivative(f, x)
    assert d.doit() == 12
    assert d.doit(deep=False) == Sum(2*n, (n, 1, 3))


def test_evalf_default():
    assert type(sin(4.0)) == Float
    assert type(re(sin(I + 1.0))) == Float
    assert type(im(sin(I + 1.0))) == Float
    assert type(sin(4)) == sin
    assert type(polygamma(2.0, 4.0)) == Float
    assert type(sin(Rational(1, 4))) == sin


def test_diff_args():
    # issue sympy/sympy#5399
    pytest.raises(ValueError, lambda: x.diff(Integer(4)))
    pytest.raises(ValueError, lambda: x.diff(x, Integer(4)))
    pytest.raises(ValueError, lambda: x.diff(Integer(4), x))
    pytest.raises(ValueError, (x*y).diff)


def test_derivative_numerically():
    z0 = random() + I*random()
    assert abs(Derivative(sin(x), x).doit_numerically(z0) - cos(z0)) < 1e-15


def test_fdiff_argument_index_error():
    class MyFunc(Function):
        nargs = 1  # define since there is no eval routine

        def fdiff(self, argindex=1):
            raise ArgumentIndexError
    mf = MyFunc(x)
    assert mf.diff(x) == Derivative(mf, x)
    pytest.raises(TypeError, lambda: MyFunc(x, x))

    pytest.raises(ArgumentIndexError, lambda: f(x).fdiff(2))

    with pytest.raises(ArgumentIndexError) as err:
        sin(x).fdiff(2)
    assert str(err.value) == ('Invalid operation with argument number 2 '
                              'for Function sin(x)')


def test_deriv_wrt_function():
    x = f(t)
    xd = diff(x, t)
    xdd = diff(xd, t)
    y = g(t)
    yd = diff(y, t)

    assert diff(x, t) == xd
    assert diff(2 * x + 4, t) == 2 * xd
    assert diff(2 * x + 4 + y, t) == 2 * xd + yd
    assert diff(2 * x + 4 + y * x, t) == 2 * xd + x * yd + xd * y
    assert diff(2 * x + 4 + y * x, x) == 2 + y
    assert (diff(4 * x**2 + 3 * x + x * y, t) == 3 * xd + x * yd + xd * y +
            8 * x * xd)
    assert (diff(4 * x**2 + 3 * xd + x * y, t) == 3 * xdd + x * yd + xd * y +
            8 * x * xd)
    assert diff(4 * x**2 + 3 * xd + x * y, xd) == 3
    assert diff(4 * x**2 + 3 * xd + x * y, xdd) == 0
    assert diff(sin(x), t) == xd * cos(x)
    assert diff(exp(x), t) == xd * exp(x)
    assert diff(sqrt(x), t) == xd / (2 * sqrt(x))


def test_diff_wrt_value():
    assert Expr()._diff_wrt is False
    assert x._diff_wrt is True
    assert f(x)._diff_wrt is True
    assert Derivative(f(x), x)._diff_wrt is True
    assert Derivative(x**2, x)._diff_wrt is False


def test_diff_wrt():
    fx = f(x)
    dfx = diff(f(x), x)
    ddfx = diff(f(x), x, x)

    assert diff(sin(fx) + fx**2, fx) == cos(fx) + 2*fx
    assert diff(sin(dfx) + dfx**2, dfx) == cos(dfx) + 2*dfx
    assert diff(sin(ddfx) + ddfx**2, ddfx) == cos(ddfx) + 2*ddfx
    assert diff(fx**2, dfx) == 0
    assert diff(fx**2, ddfx) == 0
    assert diff(dfx**2, fx) == 0
    assert diff(dfx**2, ddfx) == 0
    assert diff(ddfx**2, dfx) == 0

    assert diff(fx*dfx*ddfx, fx) == dfx*ddfx
    assert diff(fx*dfx*ddfx, dfx) == fx*ddfx
    assert diff(fx*dfx*ddfx, ddfx) == fx*dfx

    assert diff(f(x), x).diff(f(x)) == 0
    assert (sin(f(x)) - cos(diff(f(x), x))).diff(f(x)) == cos(f(x))

    assert diff(sin(fx), fx, x) == diff(sin(fx), x, fx)

    # Chain rule cases
    assert f(g(x)).diff(x) == \
        Subs(Derivative(f(x), x), (x, g(x)))*Derivative(g(x), x)
    assert diff(f(g(x), h(x)), x) == \
        Subs(Derivative(f(y, h(x)), y), (y, g(x)))*Derivative(g(x), x) + \
        Subs(Derivative(f(g(x), y), y), (y, h(x)))*Derivative(h(x), x)
    assert f(sin(x)).diff(x) == Subs(Derivative(f(x), x), (x, sin(x)))*cos(x)

    assert diff(f(g(x)), g(x)) == Subs(Derivative(f(x), x), (x, g(x)))


def test_diff_wrt_func_subs():
    assert f(g(x)).diff(x).subs({g: Lambda(x, 2*x)}).doit() == f(2*x).diff(x)


def test_diff_wrt_not_allowed():
    pytest.raises(ValueError, lambda: diff(sin(x**2), x**2))
    pytest.raises(ValueError, lambda: diff(exp(x*y), x*y))
    pytest.raises(ValueError, lambda: diff(1 + x, 1 + x))


def test_klein_gordon_lagrangian():
    m = Symbol('m')
    phi = f(x, t)

    L = -(diff(phi, t)**2 - diff(phi, x)**2 - m**2*phi**2)/2
    eqna = Eq(
        diff(L, phi) - diff(L, diff(phi, x), x) - diff(L, diff(phi, t), t), 0)
    eqnb = Eq(diff(phi, t, t) - diff(phi, x, x) + m**2*phi, 0)
    assert eqna == eqnb


def test_sho_lagrangian():
    m = Symbol('m')
    k = Symbol('k')
    x = f(t)

    L = m*diff(x, t)**2/2 - k*x**2/2
    eqna = Eq(diff(L, x), diff(L, diff(x, t), t))
    eqnb = Eq(-k*x, m*diff(x, t, t))
    assert eqna == eqnb

    assert diff(L, x, t) == diff(L, t, x)
    assert diff(L, diff(x, t), t) == m*diff(x, (t, 2))
    assert diff(L, t, diff(x, t)) == -k*x + m*diff(x, (t, 2))


def test_straight_line():
    F = f(x)
    Fd = F.diff(x)
    L = sqrt(1 + Fd**2)
    assert diff(L, F) == 0
    assert diff(L, Fd) == Fd/sqrt(1 + Fd**2)


def test_sort_variable():
    vsort = Derivative._sort_variables

    assert vsort((x, y, z)) == [x, y, z]
    assert vsort((h(x), g(x), f(x))) == [f(x), g(x), h(x)]
    assert vsort((z, y, x, h(x), g(x), f(x))) == [x, y, z, f(x), g(x), h(x)]
    assert vsort((x, f(x), y, f(y))) == [x, f(x), y, f(y)]
    assert vsort((y, x, g(x), f(x), z, h(x), y, x)) == \
        [x, y, f(x), g(x), z, h(x), x, y]
    assert vsort((z, y, f(x), x, f(x), g(x))) == [y, z, f(x), x, f(x), g(x)]
    assert vsort((z, y, f(x), x, f(x), g(x), z, z, y, x)) == \
        [y, z, f(x), x, f(x), g(x), x, y, z, z]


def test_unhandled():
    class MyExpr(Expr):
        def _eval_derivative(self, s):
            if not s.name.startswith('xi'):
                return self

    expr = MyExpr(x, y, z)
    assert diff(expr, x, y, f(x), z) == Derivative(expr, f(x), z)
    assert diff(expr, f(x), x) == Derivative(expr, f(x), x)


def test_nfloat():
    x = Symbol('x')
    eq = x**Rational(4, 3) + 4*cbrt(x)/3
    assert _aresame(nfloat(eq), x**Rational(4, 3) + (4.0/3)*cbrt(x))
    assert _aresame(nfloat(eq, exponent=True), x**(4.0/3) + (4.0/3)*x**(1.0/3))
    eq = x**Rational(4, 3) + 4*x**(x/3)/3
    assert _aresame(nfloat(eq), x**Rational(4, 3) + (4.0/3)*x**(x/3))
    big = 12345678901234567890
    # specify precision to match value used in nfloat
    Float_big = Float(big, 15)
    assert _aresame(nfloat(big), Float_big)
    assert _aresame(nfloat(big*x), Float_big*x)
    assert _aresame(nfloat(x**big, exponent=True), x**Float_big)
    assert nfloat({x: sqrt(2)}) == {x: nfloat(sqrt(2))}
    assert nfloat({sqrt(2): x}) == {sqrt(2): x}
    assert nfloat(cos(x + sqrt(2))) == cos(x + nfloat(sqrt(2)))

    # issue sympy/sympy#6342
    lamda = Symbol('lamda')
    f = x*lamda + lamda**3*(x/2 + Rational(1, 2)) + lamda**2 + Rational(1, 4)
    assert not any(a[lamda].free_symbols
                   for a in solve(f.subs({x: -0.139})))

    # issue sympy/sympy#6632
    assert nfloat(-100000*sqrt(2500000001) + 5000000001) == \
        9.99999999800000e-11

    # issue sympy/sympy#7122
    eq = cos(3*x**4 + y)*RootOf(x**5 + 3*x**3 + 1, 0)
    assert str(nfloat(eq, exponent=False, n=1)) == '-0.7*cos(3.0*x**4 + y)'


def test_sympyissue_7068():
    f = Function('f')
    y1 = Dummy('y')
    y2 = Dummy('y')
    func1 = f(a + y1 * b)
    func2 = f(a + y2 * b)
    func1_y = func1.diff(y1)
    func2_y = func2.diff(y2)
    assert func1_y != func2_y
    z1 = Subs(f(a), (a, y1))
    z2 = Subs(f(a), (a, y2))
    assert z1 != z2


def test_sympyissue_7687():
    f = Function('f')(x)
    ff = Function('f')(x)
    match_with_cache = f.match(ff)
    assert isinstance(f, type(ff))
    clear_cache()
    ff = Function('f')(x)
    assert isinstance(f, type(ff))
    assert match_with_cache == f.match(ff)


def test_sympyissue_7688():
    f = Function('f')  # actually an UndefinedFunction
    clear_cache()

    class A(UndefinedFunction):
        pass

    a = A('f')
    assert isinstance(a, type(f))


def test_mexpand():
    assert _mexpand(None) is None
    assert _mexpand(1) is Integer(1)
    assert _mexpand(x*(x + 1)**2) == (x*(x + 1)**2).expand()


def test_diff_series():
    # issue sympy/sympy#11313
    # test Derivative series & as_leading_term
    assert Derivative(x**3 + x**4, x).as_leading_term(x).doit() == 3*x**2
    s = Derivative(sin(x), x).series(x, n=3)
    assert s == Derivative(-x**3/6, x) + Derivative(x, x) + O(x**3)
    assert s.doit() == 1 - x**2/2 + O(x**3)


def test_sympyissue_12005():
    e1 = Subs(Derivative(f(x), x), (x, x))
    assert e1.diff(x) == Derivative(f(x), x, x)
    e2 = Subs(Derivative(f(x), x), (x, x**2 + 1))
    assert e2.diff(x) == 2*x*Subs(Derivative(f(x), x, x), (x, x**2 + 1))
    e3 = Subs(Derivative(f(x) + y**2 - y, y), (y, y**2))
    assert e3.diff(y) == 4*y
    e4 = Subs(Derivative(f(x + y), y), (y, x**2))
    assert e4.diff(y) == 0
    e5 = Subs(Derivative(f(x), x), (y, y), (z, z))
    assert e5.diff(x) == Derivative(f(x), x, x)
    assert f(g(x)).diff(g(x), g(x)) == Subs(Derivative(f(y), y, y), (y, g(x)))


def test_sympyissue_13098():
    assert floor(log(Float('9.9999999000000006'), 10)) == 0
    assert floor(log(Float('9.9999999899999992'), 10)) == 0
    assert floor(log(Float(('15.' + '9'*54 + '001'), dps=56), 2)) == 3
    assert floor(log(Float('16.0'), 2)) == 4
    assert floor(log(Float('99.' + '9'*23 + '007', dps=25), 10)) == 1
    assert floor(log(Float('999.99999000000003'), 10)) == 2
    assert floor(log(Float('999.999999'), 10)) == 2
    x = Float('9.'+'9'*20)
    assert floor(log(x, 10)) == 0
    assert ceiling(log(x, 10)) == 1
    assert floor(log(20 - x, 10)) == 1
    assert ceiling(log(20 - x, 10)) == 2


def test_sympyissue_6938():
    f = Function('ceil')
    assert isinstance(f(0.3), Function)
    f = Function('sin')
    assert isinstance(f(0.3), Function)
    assert isinstance(f(pi).evalf(), Function)
    assert isinstance(f(x).evalf(subs={x: 1.2}), Function)


def test_sympyissue_21773():
    assert Subs(0, (y, 1)) * Subs(0, (z, 1)) == Subs(0, (x, 1))**2
