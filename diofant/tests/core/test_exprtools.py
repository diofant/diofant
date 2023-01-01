"""Tests for tools for manipulating of large commutative expressions."""

import pytest

from diofant import (Add, Basic, Dict, Expr, Float, I, Integer, Integral,
                     Interval, Mul, O, Rational, Sum, Symbol, Tuple, cbrt,
                     collect, cos, exp, oo, root, simplify, sin, sqrt, symbols)
from diofant.abc import a, b, t, x, y, z
from diofant.core.coreerrors import NonCommutativeExpressionError
from diofant.core.exprtools import (Factors, Term, _gcd_terms, decompose_power,
                                    factor_nc, factor_terms, gcd_terms)
from diofant.core.function import _mexpand
from diofant.core.mul import _keep_coeff
from diofant.simplify.cse_opts import sub_pre


__all__ = ()


def test_decompose_power():
    assert decompose_power(x) == (x, 1)
    assert decompose_power(x**2) == (x, 2)
    assert decompose_power(x**(2*y)) == (x**y, 2)
    assert decompose_power(x**(2*y/3)) == (x**(y/3), 2)


def test_Factors():
    assert Factors() == Factors({}) == Factors(Integer(1))
    assert Factors(Integer(1)) == Factors(Factors(Integer(1)))
    assert Factors().as_expr() == 1
    assert Factors({x: 2, y: 3, sin(x): 4}).as_expr() == x**2*y**3*sin(x)**4
    assert Factors(+oo) == Factors({oo: 1})
    assert Factors(-oo) == Factors({oo: 1, -1: 1})

    f1 = Factors({oo: 1})
    f2 = Factors({oo: 1})
    assert hash(f1) == hash(f2)

    a = Factors({x: 5, y: 3, z: 7})
    b = Factors({      y: 4, z: 3, t: 10})

    assert a.mul(b) == a*b == Factors({x: 5, y: 7, z: 10, t: 10})

    assert a.div(b) == divmod(a, b) == \
        (Factors({x: 5, z: 4}), Factors({y: 1, t: 10}))
    assert a.quo(b) == a/b == Factors({x: 5, z: 4})
    assert a.rem(b) == a % b == Factors({y: 1, t: 10})

    assert a.pow(3) == a**3 == Factors({x: 15, y: 9, z: 21})
    assert b.pow(3) == b**3 == Factors({y: 12, z: 9, t: 30})

    pytest.raises(ValueError, lambda: a.pow(3.1))
    pytest.raises(ValueError, lambda: a.pow(Factors(3.1)))

    assert a.pow(0) == Factors()

    assert a.gcd(b) == Factors({y: 3, z: 3})
    assert a.lcm(b) == a.lcm(b.as_expr()) == Factors({x: 5, y: 4, z: 7, t: 10})

    a = Factors({x: 4, y: 7, t: 7})
    b = Factors({z: 1, t: 3})

    assert a.normal(b) == (Factors({x: 4, y: 7, t: 4}), Factors({z: 1}))

    assert Factors(sqrt(2)*x).as_expr() == sqrt(2)*x

    assert Factors(-I)*I == Factors()
    assert Factors({Integer(-1): Integer(3)})*Factors({Integer(-1): Integer(1), I: Integer(5)}) == \
        Factors(I)

    assert Factors(Integer(2)**x).div(Integer(3)**x) == \
        (Factors({Integer(2): x}), Factors({Integer(3): x}))
    assert Factors(2**(2*x + 2)).div(Integer(8)) == \
        (Factors({Integer(2): 2*x + 2}), Factors({Integer(8): Integer(1)}))

    # coverage
    # /!\ things break if this is not True
    assert Factors({Integer(-1): Rational(3, 2)}) == Factors({I: 1, -1: 1})
    assert Factors({I: Integer(1), Integer(-1): Rational(1, 3)}).as_expr() == I*cbrt(-1)

    assert Factors(-1.) == Factors({Integer(-1): Integer(1), Float(1.): 1})
    assert Factors(-2.) == Factors({Integer(-1): Integer(1), Float(2.): 1})
    assert Factors((-2.)**x) == Factors({Float(-2.): x})
    assert Factors(Integer(-2)) == Factors({Integer(-1): Integer(1), Integer(2): 1})
    assert Factors(Rational(1, 2)) == Factors({Integer(2): -1})
    assert Factors(Rational(3, 2)) == Factors({Integer(3): 1, Integer(2): Integer(-1)})
    assert Factors({I: Integer(1)}) == Factors(I)
    assert Factors({-1.0: 2, I: 1}) == Factors({Float(1.0): 1, I: 1})
    assert Factors({-1: -Rational(3, 2)}).as_expr() == I
    A = symbols('A', commutative=False)
    assert Factors(2*A**2) == Factors({Integer(2): 1, A**2: 1})
    assert Factors(I) == Factors({I: 1})
    assert Factors(x).normal(Integer(2)) == (Factors(x), Factors(Integer(2)))
    assert Factors(x).normal(Integer(0)) == (Factors(), Factors(Integer(0)))
    pytest.raises(ZeroDivisionError, lambda: Factors(x).div(Integer(0)))
    assert Factors(x).mul(Integer(2)) == Factors(2*x)
    assert Factors(x).mul(Integer(0)).is_zero
    assert Factors(x).mul(1/x).is_one
    assert Factors(x**sqrt(8)).as_expr() == x**(2*sqrt(2))
    assert Factors(x)**Factors(Integer(2)) == Factors(x**2)
    assert Factors(x).gcd(Integer(0)) == Factors(x)
    assert Factors(x).lcm(Integer(0)).is_zero
    assert Factors(Integer(0)).div(x) == (Factors(Integer(0)), Factors())
    assert Factors(x).div(x) == (Factors(), Factors())
    assert Factors({x: .2})/Factors({x: .2}) == Factors()
    assert Factors(x) != Factors()
    assert Factors(x) == x
    assert Factors(Integer(0)).normal(x) == (Factors(Integer(0)), Factors())
    n, d = x**(2 + y), x**2
    f = Factors(n)
    assert f.div(d) == f.normal(d) == (Factors(x**y), Factors())
    assert f.gcd(d) == Factors()
    d = x**y
    assert f.div(d) == f.normal(d) == (Factors(x**2), Factors())
    assert f.gcd(d) == Factors(d)
    n = d = 2**x
    f = Factors(n)
    assert f.div(d) == f.normal(d) == (Factors(), Factors())
    assert f.gcd(d) == Factors(d)
    n, d = 2**x, 2**y
    f = Factors(n)
    assert f.div(d) == f.normal(d) == (Factors({Integer(2): x}), Factors({Integer(2): y}))
    assert f.gcd(d) == Factors()

    assert f.div(f) == (Factors(), Factors())

    # extraction of constant only
    n = x**(x + 3)
    assert Factors(n).normal(x**-3) == (Factors({x: x + 6}), Factors({}))
    assert Factors(n).normal(x**3) == (Factors({x: x}), Factors({}))
    assert Factors(n).normal(x**4) == (Factors({x: x}), Factors({x: 1}))
    assert Factors(n).normal(x**(y - 3)) == \
        (Factors({x: x + 6}), Factors({x: y}))
    assert Factors(n).normal(x**(y + 3)) == (Factors({x: x}), Factors({x: y}))
    assert Factors(n).normal(x**(y + 4)) == \
        (Factors({x: x}), Factors({x: y + 1}))

    assert Factors(n).div(x**-3) == (Factors({x: x + 6}), Factors({}))
    assert Factors(n).div(x**3) == (Factors({x: x}), Factors({}))
    assert Factors(n).div(x**4) == (Factors({x: x}), Factors({x: 1}))
    assert Factors(n).div(x**(y - 3)) == \
        (Factors({x: x + 6}), Factors({x: y}))
    assert Factors(n).div(x**(y + 3)) == (Factors({x: x}), Factors({x: y}))
    assert Factors(n).div(x**(y + 4)) == \
        (Factors({x: x}), Factors({x: y + 1}))

    assert Factors({I: I}).as_expr() == (-1)**(I/2)
    assert Factors({-1: Rational(4, 3)}).as_expr() == -cbrt(-1)


def test_Term():
    a = Term(4*x*y**2/z/t**3)
    b = Term(2*x**3*y**5/t**3)

    assert a == Term(4, Factors({x: 1, y: 2}), Factors({z: 1, t: 3}))
    assert b == Term(2, Factors({x: 3, y: 5}), Factors({t: 3}))

    assert a.as_expr() == 4*x*y**2/z/t**3
    assert b.as_expr() == 2*x**3*y**5/t**3

    assert a.inv() == \
        Term(Rational(1, 4), Factors({z: 1, t: 3}), Factors({x: 1, y: 2}))
    assert b.inv() == Term(Rational(1, 2), Factors({t: 3}), Factors({x: 3, y: 5}))

    assert a.mul(b) == a*b == \
        Term(8, Factors({x: 4, y: 7}), Factors({z: 1, t: 6}))
    assert a.quo(b) == a/b == Term(2, Factors({}), Factors({x: 2, y: 3, z: 1}))

    assert a.pow(3) == a**3 == \
        Term(64, Factors({x: 3, y: 6}), Factors({z: 3, t: 9}))
    assert b.pow(3) == b**3 == Term(8, Factors({x: 9, y: 15}), Factors({t: 9}))

    assert a.pow(-3) == a**(-3) == \
        Term(Rational(1, 64), Factors({z: 3, t: 9}), Factors({x: 3, y: 6}))
    assert b.pow(-3) == b**(-3) == \
        Term(Rational(1, 8), Factors({t: 9}), Factors({x: 9, y: 15}))

    assert a.gcd(b) == Term(2, Factors({x: 1, y: 2}), Factors({t: 3}))
    assert a.lcm(b) == Term(4, Factors({x: 3, y: 5}), Factors({z: 1, t: 3}))

    a = Term(4*x*y**2/z/t**3)
    b = Term(2*x**3*y**5*t**7)

    assert a.mul(b) == Term(8, Factors({x: 4, y: 7, t: 4}), Factors({z: 1}))

    assert Term((2*x + 2)**3) == Term(8, Factors({x + 1: 3}), Factors({}))
    assert Term((2*x + 2)*(3*x + 6)**2) == \
        Term(18, Factors({x + 1: 1, x + 2: 2}), Factors({}))

    A = Symbol('A', commutative=False)
    pytest.raises(NonCommutativeExpressionError, lambda: Term(A))

    f1, f2 = Factors({x: 2}), Factors()
    assert Term(2, numer=f1) == Term(2, f1, f2)
    assert Term(2, denom=f1) == Term(2, f2, f1)

    pytest.raises(TypeError, lambda: a*2)
    pytest.raises(TypeError, lambda: a/3)
    pytest.raises(TypeError, lambda: a**3.1)


def test_gcd_terms():
    f = 2*(x + 1)*(x + 4)/(5*x**2 + 5) + (2*x + 2)*(x + 5)/(x**2 + 1)/5 + \
        (2*x + 2)*(x + 6)/(5*x**2 + 5)

    assert _gcd_terms(f) == (Rational(6, 5)*((1 + x)/(1 + x**2)), 5 + x, 1)
    assert _gcd_terms(Add.make_args(f)) == \
        (Rational(6, 5)*((1 + x)/(1 + x**2)), 5 + x, 1)

    newf = Rational(6, 5)*((1 + x)*(5 + x)/(1 + x**2))
    assert gcd_terms(f) == newf
    args = Add.make_args(f)
    # non-Basic sequences of terms treated as terms of Add
    assert gcd_terms(list(args)) == newf
    assert gcd_terms(tuple(args)) == newf
    assert gcd_terms(set(args)) == newf
    # but a Basic sequence is treated as a container
    assert gcd_terms(Tuple(*args)) != newf
    assert gcd_terms(Basic(Tuple(1, 3*y + 3*x*y), Tuple(1, 3))) == \
        Basic((1, 3*y*(x + 1)), (1, 3))
    # but we shouldn't change keys of a dictionary or some may be lost
    assert gcd_terms(Dict((x*(1 + y), 2), (x + x*y, y + x*y))) == \
        Dict({x*(y + 1): 2, x + x*y: y*(1 + x)})

    assert gcd_terms((2*x + 2)**3 + (2*x + 2)**2) == 4*(x + 1)**2*(2*x + 3)

    assert gcd_terms(0) == 0
    assert gcd_terms(1) == 1
    assert gcd_terms(x) == x
    assert gcd_terms(2 + 2*x) == Mul(2, 1 + x, evaluate=False)
    arg = x*(2*x + 4*y)
    garg = 2*x*(x + 2*y)
    assert gcd_terms(arg) == garg
    assert gcd_terms(sin(arg)) == sin(garg)

    # issue sympy/sympy#6139-like
    alpha, alpha1, alpha2, alpha3 = symbols('alpha:4')
    a = alpha**2 - alpha*x**2 + alpha + x**3 - x*(alpha + 1)
    rep = {alpha: (1 + sqrt(5))/2 + alpha1*x + alpha2*x**2 + alpha3*x**3}
    s = (a/(x - alpha)).subs(rep).series(x, 0, 1)
    assert simplify(collect(s, x)) == -sqrt(5)/2 - Rational(3, 2) + O(x)

    # issue sympy/sympy#5917
    assert _gcd_terms([Integer(0), Integer(0)]) == (0, 0, 1)
    assert _gcd_terms([2*x + 4]) == (2, x + 2, 1)

    eq = x/(x + 1/x)
    assert gcd_terms(eq, fraction=False) == eq


def test_factor_terms():
    A = Symbol('A', commutative=False)
    assert factor_terms(9*(x + x*y + 1) + (3*x + 3)**(2 + 2*x)) == \
        9*x*y + 9*x + _keep_coeff(Integer(3), x + 1)**_keep_coeff(Integer(2), x + 1) + 9
    assert factor_terms(9*(x + x*y + 1) + 3**(2 + 2*x)) == \
        _keep_coeff(Integer(9), 3**(2*x) + x*y + x + 1)
    assert factor_terms(3**(2 + 2*x) + a*3**(2 + 2*x)) == \
        9*3**(2*x)*(a + 1)
    assert factor_terms(x + x*A) == \
        x*(1 + A)
    assert factor_terms(sin(x + x*A)) == \
        sin(x*(1 + A))
    assert factor_terms((3*x + 3)**((2 + 2*x)/3)) == \
        _keep_coeff(Integer(3), x + 1)**_keep_coeff(Rational(2, 3), x + 1)
    assert factor_terms(x + (x*y + x)**(3*x + 3)) == \
        x + (x*(y + 1))**_keep_coeff(Integer(3), x + 1)
    assert factor_terms(a*(x + x*y) + b*(x*2 + y*x*2)) == \
        x*(a + 2*b)*(y + 1)
    i = Integral(x, (x, 0, oo))
    assert factor_terms(i) == i

    # check radical extraction
    eq = sqrt(2) + sqrt(10)
    assert factor_terms(eq) == eq
    assert factor_terms(eq, radical=True) == sqrt(2)*(1 + sqrt(5))
    eq = root(-6, 3) + root(6, 3)
    assert factor_terms(eq, radical=True) == cbrt(6)*(1 + cbrt(-1))

    eq = [x + x*y]
    ans = [x*(y + 1)]
    for c in [list, tuple, set]:
        assert factor_terms(c(eq)) == c(ans)
    assert factor_terms(Tuple(x + x*y)) == Tuple(x*(y + 1))
    assert factor_terms(Interval(0, 1)) == Interval(0, 1)
    e = 1/sqrt(a/2 + 1)
    assert factor_terms(e, clear=False) == 1/sqrt(a/2 + 1)
    assert factor_terms(e, clear=True) == sqrt(2)/sqrt(a + 2)

    eq = x/(x + 1/x) + 1/(x**2 + 1)
    assert factor_terms(eq, fraction=False) == eq
    assert factor_terms(eq, fraction=True) == 1

    assert factor_terms((1/(x**3 + x**2) + 2/x**2)*y) == \
        y*(2 + 1/(x + 1))/x**2

    # if not True, then processesing for this in factor_terms is not necessary
    assert gcd_terms(-x - y) == -x - y
    assert factor_terms(-x - y) == Mul(-1, x + y, evaluate=False)

    # if not True, then "special" processesing in factor_terms is not necessary
    assert gcd_terms(exp(Mul(-1, x + 1))) == exp(-x - 1)
    e = exp(-x - 2) + x
    assert factor_terms(e) == exp(Mul(-1, x + 2, evaluate=False)) + x
    assert factor_terms(e, sign=False) == e
    assert factor_terms(exp(-4*x - 2) - x) == -x + exp(Mul(-2, 2*x + 1, evaluate=False))


def test_xreplace():
    e = Mul(2, 1 + x, evaluate=False)
    assert e.xreplace({}) == e
    assert e.xreplace({y: x}) == e


def test_factor_nc():
    x, y = symbols('x,y')
    k = symbols('k', integer=True)
    n, m, o = symbols('n,m,o', commutative=False)

    # mul and multinomial expansion is needed
    e = x*(1 + y)**2
    assert _mexpand(e) == x + x*2*y + x*y**2

    def factor_nc_test(e):
        ex = _mexpand(e)
        assert ex.is_Add
        f = factor_nc(ex)
        assert not f.is_Add
        assert _mexpand(f) == ex

    factor_nc_test(x*(1 + y))
    factor_nc_test(n*(x + 1))
    factor_nc_test(n*(x + m))
    factor_nc_test((x + m)*n)
    factor_nc_test(n*m*(x*o + n*o*m)*n)
    s = Sum(x, (x, 1, 2))
    factor_nc_test(x*(1 + s))
    factor_nc_test(x*(1 + s)*s)
    factor_nc_test(x*(1 + sin(s)))
    factor_nc_test((1 + n)**2)

    factor_nc_test((x + n)*(x + m)*(x + y))
    factor_nc_test(x*(n*m + 1))
    factor_nc_test(x*(n*m + x))
    factor_nc_test(x*(x*n*m + 1))
    factor_nc_test(x*n*(x*m + 1))
    factor_nc_test(x*(m*n + x*n*m))
    factor_nc_test(n*(1 - m)*n**2)

    factor_nc_test((n + m)**2)
    factor_nc_test((n - m)*(n + m)**2)
    factor_nc_test((n + m)**2*(n - m))
    factor_nc_test((m - n)*(n + m)**2*(n - m))

    assert factor_nc(n*(n + n*m)) == n**2*(1 + m)
    assert factor_nc(m*(m*n + n*m*n**2)) == m*(m + n*m*n)*n
    eq = m*sin(n) - sin(n)*m
    assert factor_nc(eq) == eq

    eq = (sin(n) + x)*(cos(n) + x)
    assert factor_nc(eq.expand()) == eq

    # issue sympy/sympy#6534
    assert (2*n + 2*m).factor() == 2*(n + m)

    # issue sympy/sympy#6701
    assert factor_nc(n**k + n**(k + 1)) == n**k*(1 + n)
    assert factor_nc((m*n)**k + (m*n)**(k + 1)) == (1 + m*n)*(m*n)**k

    # issue sympy/sympy#6918
    assert factor_nc(-n*(2*x**2 + 2*x)) == -2*n*x*(x + 1)

    assert factor_nc(1 + Mul(Expr(), Expr(), evaluate=False)) == 1 + Expr()**2


def test_sympyissue_6360():
    a, b = symbols('a b')
    apb = a + b
    eq = apb + apb**2*(-2*a - 2*b)
    assert factor_terms(sub_pre(eq)) == a + b - 2*(a + b)**3


def test_sympyissue_7903():
    a = symbols(r'a', extended_real=True)
    t = exp(I*cos(a)) + exp(-I*sin(a))
    assert t.simplify()
