import pytest

from sympy import I, sqrt, log, exp, sin, asin
from sympy.core import Symbol, S, Rational, Integer, Dummy, Wild, Pow
from sympy.core.facts import InconsistentAssumptions
from sympy import simplify


def test_symbol_unset():
    x = Symbol('x', extended_real=True, integer=True)
    assert x.is_extended_real is True
    assert x.is_integer is True
    assert x.is_imaginary is None
    assert x.is_noninteger is False
    assert x.is_number is False


def test_zero():
    z = Integer(0)
    assert z.is_commutative is True
    assert z.is_integer is True
    assert z.is_rational is True
    assert z.is_algebraic is True
    assert z.is_transcendental is False
    assert z.is_extended_real is True
    assert z.is_complex is True
    assert z.is_noninteger is False
    assert z.is_irrational is False
    assert z.is_imaginary is True
    assert z.is_positive is False
    assert z.is_negative is False
    assert z.is_nonpositive is True
    assert z.is_nonnegative is True
    assert z.is_even is True
    assert z.is_odd is False
    assert z.is_finite is True
    assert z.is_infinite is False
    assert z.is_comparable is True
    assert z.is_prime is False
    assert z.is_composite is False
    assert z.is_number is True


def test_one():
    z = Integer(1)
    assert z.is_commutative is True
    assert z.is_integer is True
    assert z.is_rational is True
    assert z.is_algebraic is True
    assert z.is_transcendental is False
    assert z.is_extended_real is True
    assert z.is_complex is True
    assert z.is_noninteger is False
    assert z.is_irrational is False
    assert z.is_imaginary is False
    assert z.is_positive is True
    assert z.is_negative is False
    assert z.is_nonpositive is False
    assert z.is_nonnegative is True
    assert z.is_even is False
    assert z.is_odd is True
    assert z.is_finite is True
    assert z.is_infinite is False
    assert z.is_comparable is True
    assert z.is_prime is False
    assert z.is_number is True
    assert z.is_composite is False  # issue 8807


def test_negativeone():
    z = Integer(-1)
    assert z.is_commutative is True
    assert z.is_integer is True
    assert z.is_rational is True
    assert z.is_algebraic is True
    assert z.is_transcendental is False
    assert z.is_extended_real is True
    assert z.is_complex is True
    assert z.is_noninteger is False
    assert z.is_irrational is False
    assert z.is_imaginary is False
    assert z.is_positive is False
    assert z.is_negative is True
    assert z.is_nonpositive is True
    assert z.is_nonnegative is False
    assert z.is_even is False
    assert z.is_odd is True
    assert z.is_finite is True
    assert z.is_infinite is False
    assert z.is_comparable is True
    assert z.is_prime is False
    assert z.is_composite is False
    assert z.is_number is True


def test_infinity():
    oo = S.Infinity

    assert oo.is_commutative is True
    assert oo.is_integer is False
    assert oo.is_rational is False
    assert oo.is_algebraic is False
    assert oo.is_transcendental is False
    assert oo.is_extended_real is True
    assert oo.is_complex is False
    assert oo.is_noninteger is False
    assert oo.is_irrational is False
    assert oo.is_imaginary is False
    assert oo.is_positive is True
    assert oo.is_negative is False
    assert oo.is_nonpositive is False
    assert oo.is_nonnegative is True
    assert oo.is_even is False
    assert oo.is_odd is False
    assert oo.is_finite is False
    assert oo.is_infinite is True
    assert oo.is_comparable is True
    assert oo.is_prime is False
    assert oo.is_composite is False
    assert oo.is_number is True


def test_neg_infinity():
    mm = S.NegativeInfinity

    assert mm.is_commutative is True
    assert mm.is_integer is False
    assert mm.is_rational is False
    assert mm.is_algebraic is False
    assert mm.is_transcendental is False
    assert mm.is_extended_real is True
    assert mm.is_complex is False
    assert mm.is_noninteger is False
    assert mm.is_irrational is False
    assert mm.is_imaginary is False
    assert mm.is_positive is False
    assert mm.is_negative is True
    assert mm.is_nonpositive is True
    assert mm.is_nonnegative is False
    assert mm.is_even is False
    assert mm.is_odd is False
    assert mm.is_finite is False
    assert mm.is_infinite is True
    assert mm.is_comparable is True
    assert mm.is_prime is False
    assert mm.is_composite is False
    assert mm.is_number is True


def test_zoo():
    zoo = S.ComplexInfinity
    assert zoo.is_complex is False
    assert zoo.is_real is False
    assert zoo.is_prime is False
    assert zoo.is_infinite


def test_nan():
    nan = S.NaN

    assert nan.is_commutative is True
    assert nan.is_integer is None
    assert nan.is_rational is None
    assert nan.is_algebraic is None
    assert nan.is_transcendental is None
    assert nan.is_extended_real is None
    assert nan.is_complex is None
    assert nan.is_noninteger is None
    assert nan.is_irrational is None
    assert nan.is_imaginary is None
    assert nan.is_positive is None
    assert nan.is_negative is None
    assert nan.is_nonpositive is None
    assert nan.is_nonnegative is None
    assert nan.is_even is None
    assert nan.is_odd is None
    assert nan.is_finite is None
    assert nan.is_infinite is None
    assert nan.is_comparable is False
    assert nan.is_prime is None
    assert nan.is_composite is None
    assert nan.is_number is True


def test_pos_rational():
    r = Rational(3, 4)
    assert r.is_commutative is True
    assert r.is_integer is False
    assert r.is_rational is True
    assert r.is_algebraic is True
    assert r.is_transcendental is False
    assert r.is_extended_real is True
    assert r.is_complex is True
    assert r.is_noninteger is True
    assert r.is_irrational is False
    assert r.is_imaginary is False
    assert r.is_positive is True
    assert r.is_negative is False
    assert r.is_nonpositive is False
    assert r.is_nonnegative is True
    assert r.is_even is False
    assert r.is_odd is False
    assert r.is_finite is True
    assert r.is_infinite is False
    assert r.is_comparable is True
    assert r.is_prime is False
    assert r.is_composite is False

    r = Rational(1, 4)
    assert r.is_nonpositive is False
    assert r.is_positive is True
    assert r.is_negative is False
    assert r.is_nonnegative is True
    r = Rational(5, 4)
    assert r.is_negative is False
    assert r.is_positive is True
    assert r.is_nonpositive is False
    assert r.is_nonnegative is True
    r = Rational(5, 3)
    assert r.is_nonnegative is True
    assert r.is_positive is True
    assert r.is_negative is False
    assert r.is_nonpositive is False


def test_neg_rational():
    r = Rational(-3, 4)
    assert r.is_positive is False
    assert r.is_nonpositive is True
    assert r.is_negative is True
    assert r.is_nonnegative is False
    r = Rational(-1, 4)
    assert r.is_nonpositive is True
    assert r.is_positive is False
    assert r.is_negative is True
    assert r.is_nonnegative is False
    r = Rational(-5, 4)
    assert r.is_negative is True
    assert r.is_positive is False
    assert r.is_nonpositive is True
    assert r.is_nonnegative is False
    r = Rational(-5, 3)
    assert r.is_nonnegative is False
    assert r.is_positive is False
    assert r.is_negative is True
    assert r.is_nonpositive is True


def test_pi():
    z = S.Pi
    assert z.is_commutative is True
    assert z.is_integer is False
    assert z.is_rational is False
    assert z.is_algebraic is False
    assert z.is_transcendental is True
    assert z.is_extended_real is True
    assert z.is_complex is True
    assert z.is_noninteger is True
    assert z.is_irrational is True
    assert z.is_imaginary is False
    assert z.is_positive is True
    assert z.is_negative is False
    assert z.is_nonpositive is False
    assert z.is_nonnegative is True
    assert z.is_even is False
    assert z.is_odd is False
    assert z.is_finite is True
    assert z.is_infinite is False
    assert z.is_comparable is True
    assert z.is_prime is False
    assert z.is_composite is False


def test_E():
    z = S.Exp1
    assert z.is_commutative is True
    assert z.is_integer is False
    assert z.is_rational is False
    assert z.is_algebraic is False
    assert z.is_transcendental is True
    assert z.is_extended_real is True
    assert z.is_complex is True
    assert z.is_noninteger is True
    assert z.is_irrational is True
    assert z.is_imaginary is False
    assert z.is_positive is True
    assert z.is_negative is False
    assert z.is_nonpositive is False
    assert z.is_nonnegative is True
    assert z.is_even is False
    assert z.is_odd is False
    assert z.is_finite is True
    assert z.is_infinite is False
    assert z.is_comparable is True
    assert z.is_prime is False
    assert z.is_composite is False


def test_I():
    z = S.ImaginaryUnit
    assert z.is_commutative is True
    assert z.is_integer is False
    assert z.is_rational is False
    assert z.is_algebraic is True
    assert z.is_transcendental is False
    assert z.is_extended_real is False
    assert z.is_complex is True
    assert z.is_noninteger is False
    assert z.is_irrational is False
    assert z.is_imaginary is True
    assert z.is_positive is False
    assert z.is_negative is False
    assert z.is_nonpositive is False
    assert z.is_nonnegative is False
    assert z.is_even is False
    assert z.is_odd is False
    assert z.is_finite is True
    assert z.is_infinite is False
    assert z.is_comparable is False
    assert z.is_prime is False
    assert z.is_composite is False


def test_symbol_real():
    # issue 3848
    a = Symbol('a', extended_real=False)

    assert a.is_extended_real is False
    assert a.is_integer is False
    assert a.is_negative is False
    assert a.is_positive is False
    assert a.is_nonnegative is False
    assert a.is_nonpositive is False
    assert a.is_zero is False


def test_symbol_zero():
    x = Symbol('x', zero=True)
    assert x.is_positive is False
    assert x.is_nonpositive
    assert x.is_negative is False
    assert x.is_nonnegative
    assert x.is_zero is True
    assert x.is_nonzero is False
    assert x.is_finite is True


def test_symbol_positive():
    x = Symbol('x', positive=True)
    assert x.is_positive is True
    assert x.is_nonpositive is False
    assert x.is_negative is False
    assert x.is_nonnegative is True
    assert x.is_zero is False
    assert x.is_nonzero is True


def test_neg_symbol_positive():
    x = -Symbol('x', positive=True)
    assert x.is_positive is False
    assert x.is_nonpositive is True
    assert x.is_negative is True
    assert x.is_nonnegative is False
    assert x.is_zero is False
    assert x.is_nonzero is True


def test_symbol_nonpositive():
    x = Symbol('x', nonpositive=True)
    assert x.is_positive is False
    assert x.is_nonpositive is True
    assert x.is_negative is None
    assert x.is_nonnegative is None
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_neg_symbol_nonpositive():
    x = -Symbol('x', nonpositive=True)
    assert x.is_positive is None
    assert x.is_nonpositive is None
    assert x.is_negative is False
    assert x.is_nonnegative is True
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_symbol_falsepositive():
    x = Symbol('x', positive=False)
    assert x.is_positive is False
    assert x.is_nonpositive is None
    assert x.is_negative is None
    assert x.is_nonnegative is None
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_symbol_falsepositive_mul():
    # Explicit handling of arg.is_positive=False was added to Mul._eval_is_positive
    x = 2*Symbol('x', positive=False)
    assert x.is_positive is False  # This was None before
    assert x.is_nonpositive is None
    assert x.is_negative is None
    assert x.is_nonnegative is None
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_neg_symbol_falsepositive():
    x = -Symbol('x', positive=False)
    assert x.is_positive is None
    assert x.is_nonpositive is None
    assert x.is_negative is False
    assert x.is_nonnegative is None
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_neg_symbol_falsenegative():
    # Explicit handling of arg.is_negative=False was added to Mul._eval_is_positive
    x = -Symbol('x', negative=False)
    assert x.is_positive is False  # This was None before
    assert x.is_nonpositive is None
    assert x.is_negative is None
    assert x.is_nonnegative is None
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_symbol_falsepositive_real():
    x = Symbol('x', positive=False, extended_real=True)
    assert x.is_positive is False
    assert x.is_nonpositive is True
    assert x.is_negative is None
    assert x.is_nonnegative is None
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_neg_symbol_falsepositive_real():
    x = -Symbol('x', positive=False, extended_real=True)
    assert x.is_positive is None
    assert x.is_nonpositive is None
    assert x.is_negative is False
    assert x.is_nonnegative is True
    assert x.is_zero is None
    assert x.is_nonzero is None


def test_symbol_falsenonnegative():
    x = Symbol('x', nonnegative=False)
    assert x.is_positive is False
    assert x.is_nonpositive is None
    assert x.is_negative is None
    assert x.is_nonnegative is False
    assert x.is_zero is False
    assert x.is_nonzero is True


@pytest.mark.xfail
def test_neg_symbol_falsenonnegative():
    x = -Symbol('x', nonnegative=False)
    assert x.is_positive is None
    assert x.is_nonpositive is False  # this currently returns None
    assert x.is_negative is False  # this currently returns None
    assert x.is_nonnegative is None
    assert x.is_zero is False  # this currently returns None
    assert x.is_nonzero is True  # this currently returns None


def test_symbol_falsenonnegative_real():
    x = Symbol('x', nonnegative=False, extended_real=True)
    assert x.is_positive is False
    assert x.is_nonpositive is True
    assert x.is_negative is True
    assert x.is_nonnegative is False
    assert x.is_zero is False
    assert x.is_nonzero is True


def test_neg_symbol_falsenonnegative_real():
    x = -Symbol('x', nonnegative=False, extended_real=True)
    assert x.is_positive is True
    assert x.is_nonpositive is False
    assert x.is_negative is False
    assert x.is_nonnegative is True
    assert x.is_zero is False
    assert x.is_nonzero is True


def test_prime():
    assert S(-1).is_prime is False
    assert S(-2).is_prime is False
    assert S(-4).is_prime is False
    assert S(0).is_prime is False
    assert S(1).is_prime is False
    assert S(2).is_prime is True
    assert S(17).is_prime is True
    assert S(4).is_prime is False


def test_composite():
    assert S(-1).is_composite is False
    assert S(-2).is_composite is False
    assert S(-4).is_composite is False
    assert S(0).is_composite is False
    assert S(2).is_composite is False
    assert S(17).is_composite is False
    assert S(4).is_composite is True


def test_prime_symbol():
    x = Symbol('x', prime=True)
    assert x.is_prime is True
    assert x.is_integer is True
    assert x.is_positive is True
    assert x.is_negative is False
    assert x.is_nonpositive is False
    assert x.is_nonnegative is True

    x = Symbol('x', prime=False)
    assert x.is_prime is False
    assert x.is_integer is None
    assert x.is_positive is None
    assert x.is_negative is None
    assert x.is_nonpositive is None
    assert x.is_nonnegative is None


def test_symbol_noncommutative():
    x = Symbol('x', commutative=True)
    assert x.is_complex is None

    x = Symbol('x', commutative=False)
    assert x.is_integer is False
    assert x.is_rational is False
    assert x.is_algebraic is False
    assert x.is_irrational is False
    assert x.is_extended_real is False
    assert x.is_complex is False


def test_other_symbol():
    x = Symbol('x', integer=True)
    assert x.is_integer is True
    assert x.is_extended_real is True

    x = Symbol('x', integer=True, nonnegative=True)
    assert x.is_integer is True
    assert x.is_nonnegative is True
    assert x.is_negative is False
    assert x.is_positive is None

    x = Symbol('x', integer=True, nonpositive=True)
    assert x.is_integer is True
    assert x.is_nonpositive is True
    assert x.is_positive is False
    assert x.is_negative is None

    x = Symbol('x', odd=True)
    assert x.is_odd is True
    assert x.is_even is False
    assert x.is_integer is True

    x = Symbol('x', odd=False)
    assert x.is_odd is False
    assert x.is_even is None
    assert x.is_integer is None

    x = Symbol('x', even=True)
    assert x.is_even is True
    assert x.is_odd is False
    assert x.is_integer is True

    x = Symbol('x', even=False)
    assert x.is_even is False
    assert x.is_odd is None
    assert x.is_integer is None

    x = Symbol('x', integer=True, nonnegative=True)
    assert x.is_integer is True
    assert x.is_nonnegative is True

    x = Symbol('x', integer=True, nonpositive=True)
    assert x.is_integer is True
    assert x.is_nonpositive is True

    with pytest.raises(AttributeError):
        x.is_extended_real = False

    x = Symbol('x', algebraic=True)
    assert x.is_transcendental is False
    x = Symbol('x', transcendental=True)
    assert x.is_algebraic is False
    assert x.is_rational is False
    assert x.is_integer is False


def test_issue_3825():
    """catch: hash instability"""
    x = Symbol("x")
    y = Symbol("y")
    a1 = x + y
    a2 = y + x
    a2.is_comparable

    h1 = hash(a1)
    h2 = hash(a2)
    assert h1 == h2


def test_issue_4822():
    z = (-1)**Rational(1, 3)*(1 - I*sqrt(3))
    assert z.is_extended_real in [True, None]


def test_hash_vs_typeinfo():
    """seemingly different typeinfo, but in fact equal"""

    # the following two are semantically equal
    x1 = Symbol('x', even=True)
    x2 = Symbol('x', integer=True, odd=False)

    assert hash(x1) == hash(x2)
    assert x1 == x2


def test_hash_vs_typeinfo_2():
    """different typeinfo should mean !eq"""
    # the following two are semantically different
    x = Symbol('x')
    x1 = Symbol('x', even=True)

    assert x != x1
    assert hash(x) != hash(x1)  # This might fail with very low probability


def test_hash_vs_eq():
    """catch: different hash for equal objects"""
    a = 1 + S.Pi    # important: do not fold it into a Number instance
    ha = hash(a)  # it should be Add/Mul/... to trigger the bug

    a.is_positive   # this uses .evalf() and deduces it is positive
    assert a.is_positive is True

    # be sure that hash stayed the same
    assert ha == hash(a)

    # now b should be the same expression
    b = a.expand(trig=True)
    hb = hash(b)

    assert a == b
    assert ha == hb


def test_Add_is_pos_neg():
    # these cover lines not covered by the rest of tests in core
    n = Symbol('n', negative=True, infinite=True)
    nn = Symbol('n', nonnegative=True, infinite=True)
    np = Symbol('n', nonpositive=True, infinite=True)
    p = Symbol('p', positive=True, infinite=True)
    r = Dummy(extended_real=True, finite=False)
    x = Symbol('x')
    xf = Symbol('xb', finite=True, real=True)
    assert (n + p).is_positive is None
    assert (n + x).is_positive is None
    assert (p + x).is_positive is None
    assert (n + p).is_negative is None
    assert (n + x).is_negative is None
    assert (p + x).is_negative is None

    assert (n + xf).is_positive is False
    assert (p + xf).is_positive is True
    assert (n + xf).is_negative is True
    assert (p + xf).is_negative is False

    assert (x - S.Infinity).is_negative is None  # issue 7798
    # issue 8046, 16.2
    assert (p + nn).is_positive
    assert (n + np).is_negative
    assert (p + r).is_positive is None


def test_Add_is_imaginary():
    nn = Dummy(nonnegative=True)
    assert (I*nn + I).is_imaginary  # issue 8046, 17


def test_Add_is_algebraic():
    a = Symbol('a', algebraic=True)
    b = Symbol('a', algebraic=True)
    na = Symbol('na', algebraic=False)
    nb = Symbol('nb', algebraic=False)
    x = Symbol('x')
    assert (a + b).is_algebraic
    assert (na + nb).is_algebraic is None
    assert (a + na).is_algebraic is False
    assert (a + x).is_algebraic is None
    assert (na + x).is_algebraic is None


def test_Mul_is_algebraic():
    a = Symbol('a', algebraic=True)
    b = Symbol('a', algebraic=True)
    na = Symbol('na', algebraic=False)
    an = Symbol('an', algebraic=True, nonzero=True)
    nb = Symbol('nb', algebraic=False)
    x = Symbol('x')
    assert (a*b).is_algebraic
    assert (na*nb).is_algebraic is None
    assert (a*na).is_algebraic is None
    assert (an*na).is_algebraic is False
    assert (a*x).is_algebraic is None
    assert (na*x).is_algebraic is None


def test_Pow_is_algebraic():
    e = Symbol('e', algebraic=True)

    assert Pow(1, e, evaluate=False).is_algebraic
    assert Pow(0, e, evaluate=False).is_algebraic

    a = Symbol('a', algebraic=True)
    na = Symbol('na', algebraic=False)
    ia = Symbol('ia', algebraic=True, irrational=True)
    ib = Symbol('ib', algebraic=True, irrational=True)
    r = Symbol('r', rational=True, nonzero=True)
    x = Symbol('x')
    assert (a**r).is_algebraic
    assert (a**x).is_algebraic is None
    assert (na**r).is_algebraic is False
    assert (ia**r).is_algebraic
    assert (ia**ib).is_algebraic is False

    assert (a**e).is_algebraic is None

    # Gelfond-Schneider constant:
    assert Pow(2, sqrt(2), evaluate=False).is_algebraic is False

    assert Pow(S.GoldenRatio, sqrt(3), evaluate=False).is_algebraic is False

    # sympy/sympy#8649
    t = Symbol('t', real=True, transcendental=True)
    n = Symbol('n', integer=True)
    assert (t**n).is_algebraic is None
    assert (t**n).is_integer is None


def test_Mul_is_infinite():
    x = Symbol('x')
    f = Symbol('f', finite=True)
    i = Symbol('i', infinite=True)
    z = Dummy(zero=True)
    nzf = Dummy(finite=True, zero=False)
    from sympy import Mul
    assert (x*f).is_finite is None
    assert (x*i).is_finite is None
    assert (f*i).is_finite is False
    assert (x*f*i).is_finite is None
    assert (z*i).is_finite is False
    assert (nzf*i).is_finite is False
    assert (z*f).is_finite is True
    assert Mul(0, f, evaluate=False).is_finite is True
    assert Mul(0, i, evaluate=False).is_finite is False

    assert (x*f).is_infinite is None
    assert (x*i).is_infinite is None
    assert (f*i).is_infinite is None
    assert (x*f*i).is_infinite is None
    assert (z*i).is_infinite is S.NaN.is_infinite
    assert (nzf*i).is_infinite is True
    assert (z*f).is_infinite is False
    assert Mul(0, f, evaluate=False).is_infinite is False
    assert Mul(0, i, evaluate=False).is_infinite is S.NaN.is_infinite


def test_special_is_rational():
    i = Symbol('i', integer=True)
    i2 = Symbol('i2', integer=True)
    ni = Symbol('ni', integer=True, nonzero=True)
    r = Symbol('r', rational=True)
    rn = Symbol('r', rational=True, nonzero=True)
    nr = Symbol('nr', irrational=True)
    x = Symbol('x')
    assert sqrt(3).is_rational is False
    assert (3 + sqrt(3)).is_rational is False
    assert (3*sqrt(3)).is_rational is False
    assert exp(3).is_rational is False
    assert exp(ni).is_rational is False
    assert exp(rn).is_rational is False
    assert exp(x).is_rational is None
    assert exp(log(3), evaluate=False).is_rational is True
    assert log(exp(3), evaluate=False).is_rational is True
    assert log(3).is_rational is False
    assert log(ni + 1).is_rational is False
    assert log(rn + 1).is_rational is False
    assert log(x).is_rational is None
    assert (sqrt(3) + sqrt(5)).is_rational is None
    assert (sqrt(3) + S.Pi).is_rational is False
    assert (x**i).is_rational is None
    assert (i**i).is_rational is True
    assert (i**i2).is_rational is None
    assert (r**i).is_rational is None
    assert (r**r).is_rational is None
    assert (r**x).is_rational is None
    assert (nr**i).is_rational is None  # issue 8598
    assert (nr**Symbol('z', zero=True)).is_rational
    assert sin(1).is_rational is False
    assert sin(ni).is_rational is False
    assert sin(rn).is_rational is False
    assert sin(x).is_rational is None
    assert asin(r).is_rational is False
    assert sin(asin(3), evaluate=False).is_rational is True


@pytest.mark.xfail
def test_issue_6275():
    x = Symbol('x')
    # both zero or both Muls...but neither "change would be very appreciated.
    # This is similar to x/x => 1 even though if x = 0, it is really nan.
    assert isinstance(x*0, type(0*S.Infinity))
    if 0*S.Infinity is S.NaN:
        b = Symbol('b', finite=None)
        assert (b*0).is_zero is None


def test_sanitize_assumptions():
    # issue 6666
    for cls in (Symbol, Dummy, Wild):
        x = cls('x', extended_real=1, positive=0)
        assert x.is_extended_real is True
        assert x.is_positive is False
        assert cls('', extended_real=True, positive=None).is_positive is None
        pytest.raises(ValueError, lambda: cls('', commutative=None))
    pytest.raises(ValueError, lambda: Symbol._sanitize(dict(commutative=None)))


def test_special_assumptions():
    e = -3 - sqrt(5) + (-sqrt(10)/2 - sqrt(2)/2)**2
    assert simplify(e < 0) is S.false
    assert simplify(e > 0) is S.false
    assert (e == 0) is False  # it's not a literal 0
    assert e.equals(0) is True


def test_inconsistent():
    # cf. issues 5795 and 5545
    pytest.raises(InconsistentAssumptions, lambda: Symbol('x', extended_real=True,
           commutative=False))


def test_issue_6631():
    assert ((-1)**(I)).is_extended_real is True
    assert ((-1)**(I*2)).is_extended_real is True
    assert ((-1)**(I/2)).is_extended_real is True
    assert ((-1)**(I*S.Pi)).is_extended_real is True
    assert (I**(I + 2)).is_extended_real is True


def test_issue_2730():
    assert (1/(1 + I)).is_extended_real is False


def test_issue_4149():
    assert (3 + I).is_complex
    assert (3 + I).is_imaginary is False
    assert (3*I + S.Pi*I).is_imaginary
    y = Symbol('y', extended_real=True)
    assert (3*I + S.Pi*I + y*I).is_imaginary is True
    p = Symbol('p', positive=True)
    assert (3*I + S.Pi*I + p*I).is_imaginary
    n = Symbol('n', negative=True)
    assert (-3*I - S.Pi*I + n*I).is_imaginary

    i = Symbol('i', imaginary=True)
    assert ([(i**a).is_imaginary for a in range(4)] ==
            [False, True, False, True])

    # tests from the PR #7887:
    e = S("-sqrt(3)*I/2 + 0.866025403784439*I")
    assert e.is_extended_real is False
    assert e.is_imaginary


def test_issue_2920():
    n = Symbol('n', negative=True)
    assert sqrt(n).is_imaginary


def test_issue_7899():
    x = Symbol('x', extended_real=True)
    assert (I*x).is_extended_real is None
    assert ((x - I)*(x - 1)).is_zero is None
    assert ((x - I)*(x - 1)).is_extended_real is None


@pytest.mark.xfail
def test_issue_7993():
    x = Dummy(integer=True)
    y = Dummy(noninteger=True)
    assert (x - y).is_zero is False


def test_issue_8075():
    pytest.raises(InconsistentAssumptions, lambda: Dummy(zero=True, finite=False))
    pytest.raises(InconsistentAssumptions, lambda: Dummy(zero=True, infinite=True))


def test_issue_8642():
    x = Symbol('x', extended_real=True, integer=False)
    assert (x*2).is_integer is None


def test_issue_9165():
    z = Symbol('z', zero=True)
    f = Symbol('f', finite=False)
    assert 0/z == S.NaN
    assert 0*(1/z) == S.NaN
    assert 0*f == S.NaN
