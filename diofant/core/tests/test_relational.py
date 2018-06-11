import random
from operator import ge, gt, le, lt

import pytest

from diofant import (And, Float, Function, I, Implies, Integer, Not, Or,
                     Rational, Symbol, Wild, Xor, ceiling, false, floor, nan,
                     oo, pi, simplify, sqrt, true, zoo)
from diofant.abc import t, w, x, y, z
from diofant.core.relational import (Eq, Equality, Ge, GreaterThan, Gt, Le,
                                     LessThan, Lt, Ne, Rel, Relational,
                                     StrictGreaterThan, StrictLessThan,
                                     Unequality)
from diofant.core.relational import _Inequality as Inequality
from diofant.sets.sets import FiniteSet, Interval


__all__ = ()


def test_rel_ne():
    assert Relational(x, y, '!=') == Ne(x, y)


def test_rel_subs():
    e = Relational(x, y, '==')
    e = e.subs(x, z)

    assert isinstance(e, Equality)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '>=')
    e = e.subs(x, z)

    assert isinstance(e, GreaterThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '<=')
    e = e.subs(x, z)

    assert isinstance(e, LessThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '>')
    e = e.subs(x, z)

    assert isinstance(e, StrictGreaterThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Relational(x, y, '<')
    e = e.subs(x, z)

    assert isinstance(e, StrictLessThan)
    assert e.lhs == z
    assert e.rhs == y

    e = Eq(x, 0)
    assert e.subs(x, 0) is true
    assert e.subs(x, 1) is false


def test_wrappers():
    e = x + x**2

    res = Relational(y, e, '==')
    assert Rel(y, x + x**2, '==') == res
    assert Eq(y, x + x**2) == res

    res = Relational(y, e, '<')
    assert Lt(y, x + x**2) == res

    res = Relational(y, e, '<=')
    assert Le(y, x + x**2) == res

    res = Relational(y, e, '>')
    assert Gt(y, x + x**2) == res

    res = Relational(y, e, '>=')
    assert Ge(y, x + x**2) == res

    res = Relational(y, e, '!=')
    assert Ne(y, x + x**2) == res


def test_Eq():
    assert Eq(x**2) == Eq(x**2, 0)
    assert Eq(x**2) != Eq(x**2, 1)

    assert Eq(x, x)  # issue sympy/sympy#5719

    # issue sympy/sympy#6116
    p = Symbol('p', positive=True)
    assert Eq(p, 0) is false


def test_rel_Infinity():
    # NOTE: All of these are actually handled by diofant.core.Number, and do
    # not create Relational objects.
    assert (oo > oo) is false
    assert (oo > -oo) is true
    assert (oo > 1) is true
    assert (oo < oo) is false
    assert (oo < -oo) is false
    assert (oo < 1) is false
    assert (oo >= oo) is true
    assert (oo >= -oo) is true
    assert (oo >= 1) is true
    assert (oo <= oo) is true
    assert (oo <= -oo) is false
    assert (oo <= 1) is false
    assert (-oo > oo) is false
    assert (-oo > -oo) is false
    assert (-oo > 1) is false
    assert (-oo < oo) is true
    assert (-oo < -oo) is false
    assert (-oo < 1) is true
    assert (-oo >= oo) is false
    assert (-oo >= -oo) is true
    assert (-oo >= 1) is false
    assert (-oo <= oo) is true
    assert (-oo <= -oo) is true
    assert (-oo <= 1) is true


def test_bool():
    assert Eq(0, 0) is true
    assert Eq(1, 0) is false
    assert Ne(0, 0) is false
    assert Ne(1, 0) is true
    assert Lt(0, 1) is true
    assert Lt(1, 0) is false
    assert Le(0, 1) is true
    assert Le(1, 0) is false
    assert Le(0, 0) is true
    assert Gt(1, 0) is true
    assert Gt(0, 1) is false
    assert Ge(1, 0) is true
    assert Ge(0, 1) is false
    assert Ge(1, 1) is true
    assert Eq(I, 2) is false
    assert Ne(I, 2) is true
    pytest.raises(TypeError, lambda: Gt(I, 2))
    pytest.raises(TypeError, lambda: Ge(I, 2))
    pytest.raises(TypeError, lambda: Lt(I, 2))
    pytest.raises(TypeError, lambda: Le(I, 2))
    a = Float('.000000000000000000001', '')
    b = Float('.0000000000000000000001', '')
    assert Eq(pi + a, pi + b) is false


def test_rich_cmp():
    assert (x < y) == Lt(x, y)
    assert (x <= y) == Le(x, y)
    assert (x > y) == Gt(x, y)
    assert (x >= y) == Ge(x, y)


def test_doit():
    p = Symbol('p', positive=True)
    n = Symbol('n', negative=True)
    np = Symbol('np', nonpositive=True)
    nn = Symbol('nn', nonnegative=True)

    assert Gt(p, 0).doit() is true
    assert Gt(p, 1).doit() == Gt(p, 1)
    assert Ge(p, 0).doit() is true
    assert Le(p, 0).doit() is false
    assert Lt(n, 0).doit() is true
    assert Le(np, 0).doit() is true
    assert Gt(nn, 0).doit() == Gt(nn, 0)
    assert Lt(nn, 0).doit() is false

    assert Eq(x, 0).doit() == Eq(x, 0)


def test_new_relational():
    assert Eq(x) == Relational(x, 0)       # None ==> Equality
    assert Eq(x) == Relational(x, 0, '==')
    assert Eq(x) == Relational(x, 0, 'eq')
    assert Eq(x) == Equality(x, 0)
    assert Eq(x, -1) == Relational(x, -1)       # None ==> Equality
    assert Eq(x, -1) == Relational(x, -1, '==')
    assert Eq(x, -1) == Relational(x, -1, 'eq')
    assert Eq(x, -1) == Equality(x, -1)
    assert Eq(x) != Relational(x, 1)       # None ==> Equality
    assert Eq(x) != Relational(x, 1, '==')
    assert Eq(x) != Relational(x, 1, 'eq')
    assert Eq(x) != Equality(x, 1)
    assert Eq(x, -1) != Relational(x, 1)       # None ==> Equality
    assert Eq(x, -1) != Relational(x, 1, '==')
    assert Eq(x, -1) != Relational(x, 1, 'eq')
    assert Eq(x, -1) != Equality(x, 1)

    assert Ne(x, 0) == Relational(x, 0, '!=')
    assert Ne(x, 0) == Relational(x, 0, '<>')
    assert Ne(x, 0) == Relational(x, 0, 'ne')
    assert Ne(x, 0) == Unequality(x, 0)
    assert Ne(x, 0) != Relational(x, 1, '!=')
    assert Ne(x, 0) != Relational(x, 1, '<>')
    assert Ne(x, 0) != Relational(x, 1, 'ne')
    assert Ne(x, 0) != Unequality(x, 1)

    assert Ge(x, 0) == Relational(x, 0, '>=')
    assert Ge(x, 0) == Relational(x, 0, 'ge')
    assert Ge(x, 0) == GreaterThan(x, 0)
    assert Ge(x, 1) != Relational(x, 0, '>=')
    assert Ge(x, 1) != Relational(x, 0, 'ge')
    assert Ge(x, 1) != GreaterThan(x, 0)
    assert (x >= 1) == Relational(x, 1, '>=')
    assert (x >= 1) == Relational(x, 1, 'ge')
    assert (x >= 1) == GreaterThan(x, 1)
    assert (x >= 0) != Relational(x, 1, '>=')
    assert (x >= 0) != Relational(x, 1, 'ge')
    assert (x >= 0) != GreaterThan(x, 1)

    assert Le(x, 0) == Relational(x, 0, '<=')
    assert Le(x, 0) == Relational(x, 0, 'le')
    assert Le(x, 0) == LessThan(x, 0)
    assert Le(x, 1) != Relational(x, 0, '<=')
    assert Le(x, 1) != Relational(x, 0, 'le')
    assert Le(x, 1) != LessThan(x, 0)
    assert (x <= 1) == Relational(x, 1, '<=')
    assert (x <= 1) == Relational(x, 1, 'le')
    assert (x <= 1) == LessThan(x, 1)
    assert (x <= 0) != Relational(x, 1, '<=')
    assert (x <= 0) != Relational(x, 1, 'le')
    assert (x <= 0) != LessThan(x, 1)

    assert Gt(x, 0) == Relational(x, 0, '>')
    assert Gt(x, 0) == Relational(x, 0, 'gt')
    assert Gt(x, 0) == StrictGreaterThan(x, 0)
    assert Gt(x, 1) != Relational(x, 0, '>')
    assert Gt(x, 1) != Relational(x, 0, 'gt')
    assert Gt(x, 1) != StrictGreaterThan(x, 0)
    assert (x > 1) == Relational(x, 1, '>')
    assert (x > 1) == Relational(x, 1, 'gt')
    assert (x > 1) == StrictGreaterThan(x, 1)
    assert (x > 0) != Relational(x, 1, '>')
    assert (x > 0) != Relational(x, 1, 'gt')
    assert (x > 0) != StrictGreaterThan(x, 1)

    assert Lt(x, 0) == Relational(x, 0, '<')
    assert Lt(x, 0) == Relational(x, 0, 'lt')
    assert Lt(x, 0) == StrictLessThan(x, 0)
    assert Lt(x, 1) != Relational(x, 0, '<')
    assert Lt(x, 1) != Relational(x, 0, 'lt')
    assert Lt(x, 1) != StrictLessThan(x, 0)
    assert (x < 1) == Relational(x, 1, '<')
    assert (x < 1) == Relational(x, 1, 'lt')
    assert (x < 1) == StrictLessThan(x, 1)
    assert (x < 0) != Relational(x, 1, '<')
    assert (x < 0) != Relational(x, 1, 'lt')
    assert (x < 0) != StrictLessThan(x, 1)

    # finally, some fuzz testing
    for i in range(100):
        while 1:
            strtype, length = (chr, 65535) if random.randint(0, 1) else (chr, 255)
            relation_type = strtype(random.randint(0, length))
            if random.randint(0, 1):
                relation_type += strtype(random.randint(0, length))
            if relation_type not in ('==', 'eq', '!=', '<>', 'ne', '>=', 'ge',
                                     '<=', 'le', '>', 'gt', '<', 'lt'):
                break

        pytest.raises(ValueError, lambda: Relational(x, 1, relation_type))

    assert all(Relational(x, 0, op).rel_op == '==' for op in ('eq', '=='))
    assert all(Relational(x, 0, op).rel_op == '!=' for op in ('ne', '<>', '!='))
    assert all(Relational(x, 0, op).rel_op == '>' for op in ('gt', '>'))
    assert all(Relational(x, 0, op).rel_op == '<' for op in ('lt', '<'))
    assert all(Relational(x, 0, op).rel_op == '>=' for op in ('ge', '>='))
    assert all(Relational(x, 0, op).rel_op == '<=' for op in ('le', '<='))


def test_relational_bool_output():
    # https://github.com/sympy/sympy/issues/5931
    pytest.raises(TypeError, lambda: bool(x > 3))
    pytest.raises(TypeError, lambda: bool(x >= 3))
    pytest.raises(TypeError, lambda: bool(x < 3))
    pytest.raises(TypeError, lambda: bool(x <= 3))
    pytest.raises(TypeError, lambda: bool(Eq(x, 3)))
    pytest.raises(TypeError, lambda: bool(Ne(x, 3)))


def test_relational_logic_symbols():
    # See issue sympy/sympy#6204
    assert (x < y) & (z < t) == And(x < y, z < t)
    assert (x < y) | (z < t) == Or(x < y, z < t)
    assert ~(x < y) == Not(x < y)
    assert (x < y) >> (z < t) == Implies(x < y, z < t)
    assert (x < y) << (z < t) == Implies(z < t, x < y)
    assert (x < y) ^ (z < t) == Xor(x < y, z < t)

    assert isinstance((x < y) & (z < t), And)
    assert isinstance((x < y) | (z < t), Or)
    assert isinstance(~(x < y), GreaterThan)
    assert isinstance((x < y) >> (z < t), Implies)
    assert isinstance((x < y) << (z < t), Implies)
    assert isinstance((x < y) ^ (z < t), (Or, Xor))


def test_univariate_relational_as_set():
    assert (x > 0).as_set() == Interval(0, oo, True, True)
    assert (x >= 0).as_set() == Interval(0, oo, False, True)
    assert (x < 0).as_set() == Interval(-oo, 0, True, True)
    assert (x <= 0).as_set() == Interval(-oo, 0, True)
    assert Eq(x, 0).as_set() == FiniteSet(0)
    assert Ne(x, 0).as_set() == Interval(-oo, 0, True, True) + \
        Interval(0, oo, True, True)

    assert (x**2 >= 4).as_set() == (Interval(-oo, -2, True) +
                                    Interval(2, oo, False, True))


@pytest.mark.xfail
def test_multivariate_relational_as_set():
    assert (x*y >= 0).as_set() == Interval(0, oo)*Interval(0, oo) + \
        Interval(-oo, 0)*Interval(-oo, 0)


def test_Not():
    assert Not(Equality(x, y)) == Unequality(x, y)
    assert Not(Unequality(x, y)) == Equality(x, y)
    assert Not(StrictGreaterThan(x, y)) == LessThan(x, y)
    assert Not(StrictLessThan(x, y)) == GreaterThan(x, y)
    assert Not(GreaterThan(x, y)) == StrictLessThan(x, y)
    assert Not(LessThan(x, y)) == StrictGreaterThan(x, y)


def test_evaluate():
    assert str(Eq(x, x, evaluate=False)) == 'Eq(x, x)'
    assert Eq(x, x, evaluate=False).doit() == true
    assert str(Ne(x, x, evaluate=False)) == 'Ne(x, x)'
    assert Ne(x, x, evaluate=False).doit() == false

    assert str(Ge(x, x, evaluate=False)) == 'x >= x'
    assert str(Le(x, x, evaluate=False)) == 'x <= x'
    assert str(Gt(x, x, evaluate=False)) == 'x > x'
    assert str(Lt(x, x, evaluate=False)) == 'x < x'


def assert_all_ineq_raise_TypeError(a, b):
    pytest.raises(TypeError, lambda: a > b)
    pytest.raises(TypeError, lambda: a >= b)
    pytest.raises(TypeError, lambda: a < b)
    pytest.raises(TypeError, lambda: a <= b)
    pytest.raises(TypeError, lambda: b > a)
    pytest.raises(TypeError, lambda: b >= a)
    pytest.raises(TypeError, lambda: b < a)
    pytest.raises(TypeError, lambda: b <= a)


def assert_all_ineq_give_class_Inequality(a, b):
    """All inequality operations on `a` and `b` result in class Inequality."""
    assert isinstance(a > b,  Inequality)
    assert isinstance(a >= b, Inequality)
    assert isinstance(a < b,  Inequality)
    assert isinstance(a <= b, Inequality)
    assert isinstance(b > a,  Inequality)
    assert isinstance(b >= a, Inequality)
    assert isinstance(b < a,  Inequality)
    assert isinstance(b <= a, Inequality)


def test_imaginary_compare_raises_TypeError():
    # See issue sympy/sympy#5724
    assert_all_ineq_raise_TypeError(I, x)


def test_complex_compare_not_real():
    # two cases which are not real
    y = Symbol('y', imaginary=True, nonzero=True)
    z = Symbol('z', complex=True, extended_real=False)
    for a in (y, z):
        assert_all_ineq_raise_TypeError(2, a)
    # some cases which should remain un-evaluated
    x = Symbol('x', extended_real=True)
    z = Symbol('z', complex=True)
    for a in (x, z, t):
        assert_all_ineq_give_class_Inequality(2, a)


def test_imaginary_and_inf_compare_raises_TypeError():
    # See pull request sympy/sympy#7835
    y = Symbol('y', imaginary=True, nonzero=True)
    assert_all_ineq_raise_TypeError(oo, y)
    assert_all_ineq_raise_TypeError(-oo, y)


def test_complex_pure_imag_not_ordered():
    pytest.raises(TypeError, lambda: 2*I < 3*I)

    # more generally
    x = Symbol('x', extended_real=True, nonzero=True)
    y = Symbol('y', imaginary=True)
    z = Symbol('z', complex=True)
    assert_all_ineq_raise_TypeError(I, y)

    t = I*x   # an imaginary number, should raise errors
    assert_all_ineq_raise_TypeError(2, t)

    t = -I*y   # a real number, so no errors
    assert_all_ineq_give_class_Inequality(2, t)

    t = I*z   # unknown, should be unevaluated
    assert_all_ineq_give_class_Inequality(2, t)


def test_x_minus_y_not_same_as_x_lt_y():
    """
    A consequence of pull request sympy/sympy#7792 is that `x - y < 0` and `x < y`
    are not synonymous.
    """
    x = I + 2
    y = I + 3
    pytest.raises(TypeError, lambda: x < y)
    assert x - y < 0

    ineq = Lt(x, y, evaluate=False)
    pytest.raises(TypeError, lambda: ineq.doit())
    assert ineq.lhs - ineq.rhs < 0

    t = Symbol('t', imaginary=True, nonzero=True)
    x = 2 + t
    y = 3 + t
    ineq = Lt(x, y, evaluate=False)
    pytest.raises(TypeError, lambda: ineq.doit())
    assert ineq.lhs - ineq.rhs < 0

    # this one should give error either way
    x = I + 2
    y = 2*I + 3
    pytest.raises(TypeError, lambda: x < y)
    pytest.raises(TypeError, lambda: x - y < 0)


def test_nan_equality_exceptions():
    # See issue sympy/sympy#7774
    assert Equality(nan, nan) is false
    assert Unequality(nan, nan) is true

    # See issue sympy/sympy#7773
    A = (x, Integer(0), Rational(1, 3), pi, oo, -oo)
    assert Equality(nan, random.choice(A)) is false
    assert Equality(random.choice(A), nan) is false
    assert Unequality(nan, random.choice(A)) is true
    assert Unequality(random.choice(A), nan) is true


def test_nan_inequality_raise_errors():
    # See discussion in pull request sympy/sympy#7776.  We test inequalities with
    # a set including examples of various classes.
    for q in (x, Integer(0), Integer(10), Rational(1, 3), pi, Float(1.3), oo, -oo, nan):
        assert_all_ineq_raise_TypeError(q, nan)


def test_nan_complex_inequalities():
    # Comparisons of NaN with non-real raise errors, we're not too
    # fussy whether its the NaN error or complex error.
    for r in (I, zoo, Symbol('z', imaginary=True)):
        assert_all_ineq_raise_TypeError(r, nan)


def test_complex_infinity_inequalities():
    pytest.raises(TypeError, lambda: zoo > 0)
    pytest.raises(TypeError, lambda: zoo >= 0)
    pytest.raises(TypeError, lambda: zoo < 0)
    pytest.raises(TypeError, lambda: zoo <= 0)


def test_inequalities_symbol_name_same():
    """Using the operator and functional forms should give same results."""
    # We test all combinations from a set
    # FIXME: could replace with random selection after test passes
    A = (x, y, Integer(0), Rational(1, 3), pi, oo, -oo)
    for a in A:
        for b in A:
            assert Gt(a, b) == (a > b)
            assert Lt(a, b) == (a < b)
            assert Ge(a, b) == (a >= b)
            assert Le(a, b) == (a <= b)

    for b in (y, Integer(0), Rational(1, 3), pi, oo, -oo):
        assert Gt(x, b, evaluate=False) == (x > b)
        assert Lt(x, b, evaluate=False) == (x < b)
        assert Ge(x, b, evaluate=False) == (x >= b)
        assert Le(x, b, evaluate=False) == (x <= b)

    for b in (y, Integer(0), Rational(1, 3), pi, oo, -oo):
        assert Gt(b, x, evaluate=False) == (b > x)
        assert Lt(b, x, evaluate=False) == (b < x)
        assert Ge(b, x, evaluate=False) == (b >= x)
        assert Le(b, x, evaluate=False) == (b <= x)


def test_inequalities_symbol_name_same_complex():
    """Using the operator and functional forms should give same results.
    With complex non-real numbers, both should raise errors.
    """
    # FIXME: could replace with random selection after test passes
    for a in (x, Integer(0), Rational(1, 3), pi, oo):
        pytest.raises(TypeError, lambda: Gt(a, I))
        pytest.raises(TypeError, lambda: a > I)
        pytest.raises(TypeError, lambda: Lt(a, I))
        pytest.raises(TypeError, lambda: a < I)
        pytest.raises(TypeError, lambda: Ge(a, I))
        pytest.raises(TypeError, lambda: a >= I)
        pytest.raises(TypeError, lambda: Le(a, I))
        pytest.raises(TypeError, lambda: a <= I)


def test_inequalities_cant_sympify_other():
    # see issue sympy/sympy#7833

    bar = "foo"

    for a in (x, Integer(0), Rational(1, 3), pi, I, zoo, oo, -oo, nan):
        for op in (lt, gt, le, ge):
            pytest.raises(TypeError, lambda: op(a, bar))


def test_ineq_avoid_wild_symbol_flip():
    p = Wild('p')
    assert Gt(x, p) == Gt(x, p, evaluate=False)
    assert (x < p) == Lt(x, p, evaluate=False)  # issue sympy/sympy#7951
    # Previously failed as 'p > x':
    e = Lt(x, y).subs({y: p})
    assert e == Lt(x, p, evaluate=False)
    # Previously failed as 'p <= x':
    e = Ge(x, p).doit()
    assert e == Ge(x, p, evaluate=False)


def test_sympyissue_8245():
    a = Rational(6506833320952669167898688709329, 5070602400912917605986812821504)
    q = a.n(10)
    assert (a == q) is True
    assert (a != q) is False
    assert (a > q) is false
    assert (a < q) is false
    assert (a >= q) is true
    assert (a <= q) is true

    a = sqrt(2)
    r = Rational(str(a.n(30)))
    assert (r == a) is False
    assert (r != a) is True
    assert (r > a) is true
    assert (r < a) is false
    assert (r >= a) is true
    assert (r <= a) is false
    a = sqrt(2)
    r = Rational(str(a.n(29)))
    assert (r == a) is False
    assert (r != a) is True
    assert (r > a) is false
    assert (r < a) is true
    assert (r >= a) is false
    assert (r <= a) is true


def test_sympyissue_8449():
    p = Symbol('p', nonnegative=True)
    assert Lt(-oo, p)
    assert Ge(-oo, p) is false
    assert Gt(oo, -p)
    assert Le(oo, -p) is false


def test_simplify():
    assert simplify(x*(y + 1) - x*y - x + 1 < x) == (x > 1)
    assert simplify(Integer(1) < -x) == (x < -1)

    # issue sympy/sympy#10304
    d = -(3*2**pi)**(1/pi) + 2*3**(1/pi)
    assert d.is_real
    assert simplify(Eq(1 + I*d, 0)) is False
    assert simplify(Ne(1 + I*d, 0)) is True


def test_equals():
    f = Function('f')
    assert Eq(x, 1).equals(1) is not True
    assert Eq(x, 1).equals(Eq(x*(y + 1) - x*y - x + 1, x))
    assert Eq(x, y).equals(x < y, True) is False
    assert Eq(x, f(1)).equals(Eq(x, f(2)), True) == f(1) - f(2)
    assert Eq(f(1), y).equals(Eq(f(2), y), True) == f(1) - f(2)
    assert Eq(x, f(1)).equals(Eq(f(2), x), True) == f(1) - f(2)
    assert Eq(f(1), x).equals(Eq(x, f(2)), True) == f(1) - f(2)
    assert Eq(w, x).equals(Eq(y, z), True) is False
    assert Eq(f(1), f(2)).equals(Eq(f(3), f(4)), True) == f(1) - f(3) + f(4) - f(2)
    assert Eq(f(1), f(2)).equals(Eq(f(3), f(4))) is None
    assert (x < y).equals(y > x, True) is True
    assert (x < y).equals(y >= x, True) is False
    assert (x < y).equals(z < y, True) is False
    assert (x < y).equals(x < z, True) is False
    assert (x < f(1)).equals(x < f(2), True) == f(1) - f(2)
    assert (f(1) < x).equals(f(2) < x, True) == f(1) - f(2)


def test_reversed():
    assert (x < y).reversed == (y > x)
    assert (x <= y).reversed == (y >= x)
    assert Eq(x, y, evaluate=False).reversed == Eq(y, x, evaluate=False)
    assert Ne(x, y, evaluate=False).reversed == Ne(y, x, evaluate=False)
    assert (x >= y).reversed == (y <= x)
    assert (x > y).reversed == (y < x)


def test_canonical():
    one = Integer(1)

    def unchanged(v):
        c = v.canonical
        return v.is_Relational and c.is_Relational and v == c

    def isreversed(v):
        return v.canonical == v.reversed

    assert unchanged(x < one)
    assert unchanged(x <= one)
    assert isreversed(Eq(one, x, evaluate=False))
    assert unchanged(Eq(x, one, evaluate=False))
    assert isreversed(Ne(one, x, evaluate=False))
    assert unchanged(Ne(x, one, evaluate=False))
    assert unchanged(x >= one)
    assert unchanged(x > one)

    assert unchanged(x < y)
    assert unchanged(x <= y)
    assert isreversed(Eq(y, x, evaluate=False))
    assert unchanged(Eq(x, y, evaluate=False))
    assert isreversed(Ne(y, x, evaluate=False))
    assert unchanged(Ne(x, y, evaluate=False))
    assert isreversed(x >= y)
    assert isreversed(x > y)
    assert (-x < 1).canonical == (x > -1)
    assert isreversed(-x > y)


@pytest.mark.xfail
def test_sympyissue_8444():
    x = Symbol('x', extended_real=True)
    assert (x <= oo) == (x >= -oo) == true

    x = Symbol('x', real=True)
    assert x >= floor(x)
    assert (x < floor(x)) is false
    assert Gt(x, floor(x)) == Gt(x, floor(x), evaluate=False)
    assert Ge(x, floor(x)) == Ge(x, floor(x), evaluate=False)
    assert x <= ceiling(x)
    assert (x > ceiling(x)) is false
    assert Lt(x, ceiling(x)) == Lt(x, ceiling(x), evaluate=False)
    assert Le(x, ceiling(x)) == Le(x, ceiling(x), evaluate=False)
    i = Symbol('i', integer=True)
    assert (i > floor(i)) is false
    assert (i < ceiling(i)) is false


def test_sympyissue_10633():
    assert Eq(True, False) is false
    assert Eq(False, True) is false
    assert Eq(True, True) is true
    assert Eq(False, False) is true
