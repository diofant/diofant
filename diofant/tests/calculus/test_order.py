import pytest

from diofant import (Add, Derivative, Function, I, Integer, Integral, O,
                     Rational, Symbol, conjugate, cos, digamma, exp, expand,
                     factorial, log, nan, oo, pi, sin, sqrt, transpose)
from diofant.abc import w, x, y, z


__all__ = ()


def test_caching_bug():
    # needs to be a first test, so that all caches are clean
    # cache it
    O(w)
    # and test that this won't raise an exception
    O(w**(-1/x/log(3)*log(5)), w)


def test_free_symbols():
    assert O(x).free_symbols == {x}
    assert O(1, x).free_symbols == {x}


def test_simple_1():
    o = Integer(0)
    assert O(2*x) == O(x)
    assert O(x)*3 == O(x)
    assert -28*O(x) == O(x)
    assert O(O(x)) == O(x)
    assert O(exp(x)) == O(1, x)
    assert O(exp(1/x)).expr == exp(1/x)
    assert O(x*exp(1/x)).expr == x*exp(1/x)
    assert O(x**(o/3), x).expr == x**(o/3)
    assert O(x**(5*o/3), x).expr == x**(5*o/3)
    assert O(x**2 + x + y, x) == O(1, x)
    assert O(x**2 + x + y, y) == O(1, y)
    pytest.raises(TypeError, lambda: O(x, 2 - x))
    pytest.raises(ValueError, lambda: O(x, x, x**2))
    pytest.raises(ValueError, lambda: O(x*y))

    assert O(x**2).is_commutative


def test_simple_2():
    assert O(2*x)*x == O(x**2)
    assert O(2*x)/x == O(1, x)
    assert O(2*x)*x*exp(1/x) == O(x**2*exp(1/x))
    assert (O(2*x)*x*exp(1/x)/log(x)**3).expr == x**2*exp(1/x)*log(x)**-3


def test_simple_3():
    assert O(x) + x == O(x)
    assert O(x) + 2 == 2 + O(x)
    assert O(x) + x**2 == O(x)
    assert O(x) + 1/x == 1/x + O(x)
    assert O(1/x) + 1/x**2 == 1/x**2 + O(1/x)
    assert O(x) + exp(1/x) == O(x) + exp(1/x)


def test_simple_4():
    assert O(x)**2 == O(x**2)


def test_simple_5():
    assert O(x) + O(x**2) == O(x)
    assert O(x) + O(x**-2) == O(x**-2)
    assert O(x) + O(1/x) == O(1/x)


def test_simple_6():
    assert O(x) - O(x) == O(x)
    assert O(x) + O(x**2) == O(x)
    assert O(x) + O(exp(1/x)) == O(exp(1/x))
    assert O(x**3) + O(exp(2/x)) == O(exp(2/x))
    assert O(x**-3) + O(exp(2/x)) == O(exp(2/x))


def test_simple_8():
    assert O(sqrt(-x)) == O(sqrt(x))
    assert O(x**2*sqrt(x)) == O(x**Rational(5, 2))
    assert O(x**3*sqrt(-(-x)**3)) == O(x**Rational(9, 2))
    assert O(x**Rational(3, 2)*sqrt((-x)**3)) == O(x**3)
    assert O(x*(-2*x)**(I/2)) == O(x*(-x)**(I/2))
    assert O(sqrt((-x)**I)) == O(sqrt((-x)**I), evaluate=False)
    assert O(sqrt(-x**I)) == O(sqrt(-x**I), evaluate=False)


def test_as_expr_variables():
    assert O(x).as_expr_variables(None) == (x, (x, 0))
    assert O(x).as_expr_variables((x, 0)) == (x, (x, 0))


def test_contains():
    assert O(1, x).contains(O(1, x))
    assert O(x).contains(O(x))
    assert O(x).contains(O(x**2))
    assert not O(x**2).contains(O(x))
    assert not O(x).contains(O(1/x))
    assert not O(1/x).contains(O(exp(1/x)))
    assert not O(x).contains(O(exp(1/x)))
    assert O(1/x).contains(O(x))
    assert O(exp(1/x)).contains(O(x))
    assert O(exp(1/x)).contains(O(1/x))
    assert O(exp(1/x)).contains(O(exp(1/x)))
    assert O(exp(2/x)).contains(O(exp(1/x)))
    assert not O(exp(1/x)).contains(O(exp(2/x)))

    assert O(x).contains(O(y)) is None

    assert O(sin(1/x**2)).contains(O(cos(1/x**2))) is None
    assert O(cos(1/x**2)).contains(O(sin(1/x**2))) is None

    q = Symbol('q', positive=True)
    assert O(x**8).contains(x**(q + 7)) is None
    assert O(x**8).contains(x**(q + 8))

    pytest.raises(TypeError, lambda: O(x**y, x) in O(x**z, x))
    assert O(1/x) not in O(x)
    assert O(x) in O(1/x)


def test_add_1():
    assert O(x + x) == O(x)
    assert O(3*x - 2*x**2) == O(x)
    assert O(1 + x) == O(1, x)
    assert O(1 + 1/x) == O(1/x)
    assert O(log(x) + 1/log(x)) == O(log(x))
    assert O(exp(1/x) + x) == O(exp(1/x))
    assert O(exp(1/x) + 1/x**20) == O(exp(1/x))


def test_log_args():
    assert O(log(x)) + O(log(2*x)) == O(log(x))
    assert O(log(x)) + O(log(x**3)) == O(log(x))


def test_sympyissue_3468():
    y = Symbol('y', negative=True)
    z = Symbol('z', complex=True)

    # check that Order does not modify assumptions about symbols
    O(x)
    O(y)
    O(z)

    assert x.is_positive is None
    assert y.is_positive is False
    assert z.is_positive is None


def test_leading_order():
    assert (x + 1 + 1/x**5)._extract_leading_order(x) == ((1/x**5, O(1/x**5)),)
    assert (1 + 1/x)._extract_leading_order(x) == ((1/x, O(1/x)),)
    assert (1 + x)._extract_leading_order(x) == ((1, O(1, x)),)
    assert (1 + x**2)._extract_leading_order(x) == ((1, O(1, x)),)
    assert (2 + x**2)._extract_leading_order(x) == ((2, O(1, x)),)
    assert (x + x**2)._extract_leading_order(x) == ((x, O(x)),)


def test_leading_order2():
    assert set((2 + pi + x**2)._extract_leading_order(x)) == {(pi, O(1, x)),
                                                              (2, O(1, x))}
    assert set((2*x + pi*x + x**2)._extract_leading_order(x)) == {(2*x, O(x)),
                                                                  (x*pi, O(x))}


def test_order_leadterm():
    assert O(x**2).as_leading_term(x) == O(x**2)


def test_order_symbols():
    e = x*y*sin(x)*Integral(x, (x, 1, 2))
    assert O(e, x) == O(x**2)


def test_nan():
    assert O(nan, x) == nan
    assert not O(x).contains(nan)


def test_O1():
    assert O(1, x) * x == O(x)
    assert O(1, y) * x == O(1, y)


def test_getn():
    # other lines are tested incidentally by the suite
    assert O(x).getn() == 1
    assert O(x/log(x)).getn() == 1
    assert O(x**2/log(x)**2).getn() == 2
    assert O(x*log(x)).getn() == 1
    pytest.raises(NotImplementedError, (O(x) + O(y)).getn)
    pytest.raises(NotImplementedError, O(x**y*log(x)**z, x, 0).getn)
    pytest.raises(NotImplementedError, O(x**pi*log(x), x, 0).getn)

    f = Function('f')
    pytest.raises(NotImplementedError, O(f(x)).getn)


def test_diff():
    assert O(1, x).diff(x) == Derivative(O(1, x), x)
    assert O(x**2).diff(x) == Derivative(O(x**2), x)


def test_getO():
    assert x.getO() is None
    assert x.removeO() == x
    assert O(x).getO() == O(x)
    assert O(x).removeO() == 0
    assert (z + O(x) + O(y)).getO() == O(x) + O(y)
    assert (z + O(x) + O(y)).removeO() == z
    pytest.raises(NotImplementedError, (O(x) + O(y)).getn)


def test_leading_term():
    assert O(1/digamma(1/x)) == O(1/log(x))


def test_eval():
    assert O(x).subs({O(x): 1}) == 1
    assert O(x).subs({x: y}) == O(y)
    assert O(x).subs({y: x}) == O(x)
    assert O(x).subs({x: x + y}) == O(x + y, x, -y)


def test_pow():
    assert (1/O(x)).is_Pow


def test_sympyissue_4855():
    assert 1/O(x) != O(1/x)
    assert 1/O(x, x, oo) != O(1/x, x, oo)

    f = Function('f')
    assert 1/O(f(x)) != O(1/x)


def test_order_conjugate_transpose():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', imaginary=True)
    assert conjugate(O(x)) == O(conjugate(x))
    assert conjugate(O(y)) == O(conjugate(y))
    assert conjugate(O(x**2)) == O(conjugate(x)**2)
    assert conjugate(O(y**2)) == O(conjugate(y)**2)
    assert conjugate(O(z)) == conjugate(O(z), evaluate=False)
    assert transpose(O(x)) == O(transpose(x))
    assert transpose(O(y)) == O(transpose(y))
    assert transpose(O(x**2)) == O(transpose(x)**2)
    assert transpose(O(y**2)) == O(transpose(y)**2)
    assert transpose(O(z)) == transpose(O(z), evaluate=False)


def test_order_noncommutative():
    A = Symbol('A', commutative=False)
    assert O(A + A*x, x) == O(1, x)
    assert (A + A*x)*O(x) == O(x)
    assert (A*x)*O(x) == O(x**2, x)
    assert expand((1 + O(x))*A*A*x) == A*A*x + O(x**2, x)
    assert expand((A*A + O(x))*x) == A*A*x + O(x**2, x)
    assert expand((A + O(x))*A*x) == A*A*x + O(x**2, x)


def test_sympyissue_6753():
    assert (1 + x**2)**10000*O(x) == O(x)


def test_sympyissue_7872():
    assert O(x**3).subs({x: exp(-x**2)}) == O(exp(-3*x**2), x, -oo)


def test_order_at_infinity():
    assert O(1 + x, x, oo) == O(x, x, oo)
    assert O(3*x, x, oo) == O(x, x, oo)
    assert O(x, x, oo)*3 == O(x, x, oo)
    assert -28*O(x, x, oo) == O(x, x, oo)
    assert O(3, x, oo) == O(1, x, oo)
    assert O(x**2 + x + y, x, oo) == O(x**2, x, oo)
    assert O(x**2 + x + y, y, oo) == O(y, y, oo)

    assert O(2*x, x, oo)*x == O(x**2, x, oo)
    assert O(2*x, x, oo)/x == O(1, x, oo)
    assert O(2*x, x, oo)*x*exp(1/x) == O(x**2*exp(1/x), x, oo)
    assert O(2*x, x, oo)*x*exp(1/x)/log(x)**3 == O(x**2*exp(1/x)*log(x)**-3, x, oo)

    assert O(x, x, oo) + 1/x == 1/x + O(x, x, oo) == O(x, x, oo)
    assert O(x, x, oo) + 1 == 1 + O(x, x, oo) == O(x, x, oo)
    assert O(x, x, oo) + x == x + O(x, x, oo) == O(x, x, oo)
    assert O(x, x, oo) + x**2 == x**2 + O(x, x, oo)
    assert O(1/x, x, oo) + 1/x**2 == 1/x**2 + O(1/x, x, oo) == O(1/x, x, oo)
    assert O(x, x, oo) + exp(1/x) == exp(1/x) + O(x, x, oo)

    assert O(x, x, oo)**2 == O(x**2, x, oo)

    assert O(x, x, oo) + O(x**2, x, oo) == O(x**2, x, oo)
    assert O(x, x, oo) + O(x**-2, x, oo) == O(x, x, oo)
    assert O(x, x, oo) + O(1/x, x, oo) == O(x, x, oo)

    assert O(x, x, oo) - O(x, x, oo) == O(x, x, oo)
    assert O(x, x, oo) + O(1, x, oo) == O(x, x, oo)
    assert O(x, x, oo) + O(x**2, x, oo) == O(x**2, x, oo)
    assert O(1/x, x, oo) + O(1, x, oo) == O(1, x, oo)
    assert O(x, x, oo) + O(exp(1/x), x, oo) == O(x, x, oo)
    assert O(x**3, x, oo) + O(exp(2/x), x, oo) == O(x**3, x, oo)
    assert O(x**-3, x, oo) + O(exp(2/x), x, oo) == O(exp(2/x), x, oo)

    # issue sympy/sympy#7207
    assert O(exp(x), x, oo).expr == O(2*exp(x), x, oo).expr == exp(x)
    assert O(y**x, x, oo).expr == O(2*y**x, x, oo).expr == y**x

    # issue sympy/sympy#9917
    assert O(x*sin(x) + 1, x, oo) != O(x*sin(x), x, oo)

    # issue sympy/sympy#15539
    assert O(x**-6, x, -oo) == O(x**(-6), x, -oo, evaluate=False)


def test_mixing_order_at_zero_and_infinity():
    assert (O(x, x, 0) + O(x, x, oo)).is_Add
    assert O(x, x, 0) + O(x, x, oo) == O(x, x, oo) + O(x, x, 0)
    assert O(O(x, x, oo), point=oo) == O(x, x, oo)

    # not supported (yet)
    pytest.raises(NotImplementedError, lambda: O(x, x, 0)*O(x, x, oo))
    pytest.raises(NotImplementedError, lambda: O(x, x, oo)*O(x, x, 0))
    pytest.raises(NotImplementedError, lambda: O(O(x), x, oo))


def test_order_at_some_point():
    assert O(x, x, 1) == O(1, x, 1)
    assert O(2*x - 2, x, 1) == O(x - 1, x, 1)
    assert O(-x + 1, x, 1) == O(x - 1, x, 1)
    assert O(x - 1, x, 1)**2 == O((x - 1)**2, x, 1)
    assert O(x - 2, x, 2) - O(x - 2, x, 2) == O(x - 2, x, 2)


def test_order_subs_limits():
    # issue sympy/sympy#3333
    assert (1 + O(x)).subs({x: 1/x}) == 1 + O(1/x, x, oo)
    assert (1 + O(x)).limit(x, 0) == 1
    # issue sympy/sympy#5769
    assert ((x + O(x**2))/x).limit(x, 0) == 1

    assert O(x**2).subs({x: y - 1}) == O((y - 1)**2, y, 1)
    assert O(10*x**2, x, 2).subs({x: y - 1}) == O(1, y, 3)

    assert O(1/x, x, oo).subs({x: +I*x}) == O(1/x, x, -I*oo)
    assert O(1/x, x, oo).subs({x: -I*x}) == O(1/x, x, +I*oo)


def test_sympyissue_9351():
    assert exp(x).series(x, 10, 1) == exp(10) + O(x - 10, x, 10)


def test_sympyissue_7599():
    n = Symbol('n', integer=True)
    assert O(x**n, x) + O(x**2) == Add(O(x**2), O(x**n, x), evaluate=False)


def test_sympyissue_22836():
    assert O(2**x + factorial(x), x, oo) == O(factorial(x), x, oo)
    assert O(2**x + factorial(x) + x**x, x, oo) == O((1/x)**(-x), x, oo)
    assert O(x + factorial(x), x, oo) == O(factorial(x), x, oo)
