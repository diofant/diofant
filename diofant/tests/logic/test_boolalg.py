import itertools

import pytest

from diofant import (ITE, And, EmptySet, Equality, Equivalent, Implies,
                     Integer, Interval, Nand, Nor, Not, Or, S, Unequality,
                     Union, Xor, false, oo, simplify_logic, sqrt, to_cnf,
                     to_dnf, to_nnf, true)
from diofant.abc import a, b, c, d, w, x, y, z
from diofant.logic.boolalg import (BooleanAtom, BooleanFunction, _POSform,
                                   _SOPform, is_cnf, is_dnf, is_literal,
                                   is_nnf, to_int_repr)


__all__ = ()


def test_overloading():
    assert a & b == And(a, b)
    assert a | b == Or(a, b)
    assert a >> b == Implies(a, b)
    assert ~a == Not(a)
    assert a ^ b == Xor(a, b)


def test_And():
    assert And() is true
    assert And(a) == a
    assert And(True) is true
    assert And(False) is false
    assert And(True, True) is true
    assert And(True, False) is false
    assert And(True, False, evaluate=False) is not false
    assert And(False, False) is false
    assert And(True, a) == a
    assert And(False, a) is false
    assert And(True, True, True) is true
    assert And(True, True, a) == a
    assert And(True, False, a) is false
    assert And(2, a) == a
    assert And(2, 3) is true
    assert And(a < 1, a >= 1) is false

    e = a > 1
    assert And(e, e.canonical) == e.canonical

    g, l, ge, le = a > b, b < a, a >= b, b <= a
    assert And(g, l, ge, le) == And(l, le)


def test_Or():
    assert Or() is false
    assert Or(a) == a
    assert Or(True) is true
    assert Or(False) is false
    assert Or(True, True) is true
    assert Or(True, False) is true
    assert Or(False, False) is false
    assert Or(True, a) is true
    assert Or(False, a) == a
    assert Or(True, False, False) is true
    assert Or(True, False, a) is true
    assert Or(False, False, a) == a
    assert Or(2, a) is true
    assert Or(a < 1, a >= 1) is true

    e = a > 1
    assert Or(e, e.canonical) == e

    g, l, ge, le = a > b, b < a, a >= b, b <= a
    assert Or(g, l, ge, le) == Or(g, ge)


def test_Xor():
    assert Xor() is false
    assert Xor(a) == a
    assert Xor(a, a) is false
    assert Xor(True, a, a) is true
    assert Xor(a, a, a, a, a) == a
    assert Xor(True, False, False, a, b) == ~Xor(a, b)
    assert Xor(True) is true
    assert Xor(False) is false
    assert Xor(True, True) is false
    assert Xor(True, False) is true
    assert Xor(False, False) is false
    assert Xor(True, a) == ~a
    assert Xor(False, a) == a
    assert Xor(True, False, False) is true
    assert Xor(True, False, a) == ~a
    assert Xor(False, False, a) == a
    assert isinstance(Xor(a, b), Xor)
    assert Xor(a, b, Xor(c, d)) == Xor(a, b, c, d)
    assert Xor(a, b, Xor(b, c)) == Xor(a, c)
    assert Xor(a < 1, a >= 1, b) == Xor(0, 1, b) == Xor(1, 0, b)

    e = a > 1
    assert Xor(e, e.canonical) == Xor(0, 0) == Xor(1, 1)

    e = Integer(1) < a
    assert e != e.canonical
    assert Xor(e, e.canonical) is false

    assert Xor(a > 1, b > c) == Xor(a > 1, b > c, evaluate=False)


def test_Not():
    pytest.raises(TypeError, lambda: Not(True, False))
    assert Not(True) is false
    assert Not(False) is true
    assert Not(0) is true
    assert Not(1) is false
    assert Not(2) is false
    assert Not(Unequality(a, b)) == Equality(a, b)


def test_Nand():
    assert Nand() is false
    assert Nand(a) == ~a
    assert Nand(True) is false
    assert Nand(False) is true
    assert Nand(True, True) is false
    assert Nand(True, False) is true
    assert Nand(False, False) is true
    assert Nand(True, a) == ~a
    assert Nand(False, a) is true
    assert Nand(True, True, True) is false
    assert Nand(True, True, a) == ~a
    assert Nand(True, False, a) is true


def test_Nor():
    assert Nor() is true
    assert Nor(a) == ~a
    assert Nor(True) is false
    assert Nor(False) is true
    assert Nor(True, True) is false
    assert Nor(True, False) is false
    assert Nor(False, False) is true
    assert Nor(True, a) is false
    assert Nor(False, a) == ~a
    assert Nor(True, True, True) is false
    assert Nor(True, True, a) is false
    assert Nor(True, False, a) is false


def test_Implies():
    pytest.raises(ValueError, lambda: Implies(a, b, c))

    assert Implies(True, True) is true
    assert Implies(True, False) is false
    assert Implies(False, True) is true
    assert Implies(False, False) is true
    assert Implies(0, a) is true
    assert Implies(1, 1) is true
    assert Implies(1, 0) is false
    assert a >> b == b << a
    assert (a < 1) >> (a >= 1) == (a >= 1)
    assert (a < 1) >> (1 > a) is true
    assert (a < 1) >> (Integer(1) > a) is true
    assert a >> a is true
    assert ((a < 1) >> (b >= 1)) == Implies(a < 1, b >= 1, evaluate=False)


def test_Equivalent():
    assert Equivalent(a, b) == Equivalent(b, a) == Equivalent(a, b, a)
    assert Equivalent() is true
    assert Equivalent(a, a) == Equivalent(a) is true
    assert Equivalent(True, True) == Equivalent(False, False) is true
    assert Equivalent(True, False) == Equivalent(False, True) is false
    assert Equivalent(a, True) == a
    assert Equivalent(a, False) == ~a
    assert Equivalent(a, b, True) == a & b
    assert Equivalent(a, b, False) == ~a & ~b
    assert Equivalent(1, a) == a
    assert Equivalent(0, a) == Not(a)
    assert Equivalent(a, Equivalent(b, c)) != Equivalent(Equivalent(a, b), c)
    assert Equivalent(a < 1, a >= 1) is false
    assert Equivalent(a < 1, a >= 1, 0) is false
    assert Equivalent(a < 1, a >= 1, 1) is false
    assert Equivalent(a < 1, 1 > a) == Equivalent(1, 1) == Equivalent(0, 0)
    assert Equivalent(a < 1, Integer(1) > a) == Equivalent(1, 1) == Equivalent(0, 0)
    assert Equivalent(a < 1, b >= 1) == Equivalent(b >= 1, a < 1, evaluate=False)


def test_equals():
    assert (~(a | b)).equals(~a & ~b) is True
    assert Equivalent(a, b).equals((a >> b) & (b >> a)) is True
    assert ((a | ~b) & (~a | b)).equals((~a & ~b) | (a & b)) is True
    assert (a >> b).equals(~a >> ~b) is False
    assert (a >> (b >> a)).equals(a >> (c >> a)) is False
    pytest.raises(NotImplementedError, lambda: (a & (a < b)).equals(a & (b > a)))


def test_simplification():
    set1 = [[0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 0]]
    set2 = [[0, 0, 0], [0, 1, 0], [1, 0, 1], [1, 1, 1]]
    assert _SOPform([x, y, z], set1) == (~x & z) | (~z & x)
    assert ~_SOPform([x, y, z], set2) == ~((~x & ~z) | (x & z))
    assert _POSform([x, y, z], set1 + set2) is true
    assert _SOPform([x, y, z], set1 + set2) is true
    assert _SOPform([w, x, y, z], set1 + set2) is true

    minterms = [[0, 0, 0, 1], [0, 0, 1, 1], [0, 1, 1, 1], [1, 0, 1, 1],
                [1, 1, 1, 1]]
    dontcares = [[0, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 1]]
    assert _SOPform([w, x, y, z], minterms, dontcares) == (~w & z) | (y & z)
    assert _POSform([w, x, y, z], minterms, dontcares) == (~w | y) & z

    # test simplification
    assert (a & (b | c)).simplify() == a & (b | c)
    assert ((a & b) | (a & c)).simplify() == a & (b | c)
    assert (a >> b).simplify() == ~a | b
    assert to_dnf(Equivalent(a, b), simplify=True) == (a & b) | (~a & ~b)
    assert (Equality(a, 2) & c).simplify() == Equality(a, 2) & c
    assert (Equality(a, 2) & a).simplify() == Equality(a, 2) & a
    assert (Equality(a, b) & c).simplify() == Equality(a, b) & c
    assert ((Equality(a, 3) & b) | (Equality(a, 3) & c)).simplify() == Equality(a, 3) & (b | c)

    e = a & (x**2 - x)
    assert e.simplify() == a & x*(x - 1)
    assert simplify_logic(e, deep=False) == e

    pytest.raises(ValueError, lambda: simplify_logic(a & (b | c), form='spam'))

    e = x & y ^ z | (z ^ x)
    res = [(x & ~z) | (z & ~x) | (z & ~y), (x & ~y) | (x & ~z) | (z & ~x)]
    assert to_dnf(e, simplify=True) in res
    assert _SOPform([z, y, x], [[0, 0, 1], [0, 1, 1],
                                [1, 0, 0], [1, 0, 1], [1, 1, 0]]) == res[1]

    # check input
    ans = _SOPform([x, y], [[1, 0]])
    assert _SOPform([x, y], [[1, 0]]) == ans
    assert _POSform([x, y], [[1, 0]]) == ans

    pytest.raises(ValueError, lambda: _SOPform([x], [[1]], [[1]]))
    assert _SOPform([x], [[1]], [[0]]) is true
    assert _SOPform([x], [[0]], [[1]]) is true
    assert _SOPform([x], [], []) is false

    pytest.raises(ValueError, lambda: _POSform([x], [[1]], [[1]]))
    assert _POSform([x], [[1]], [[0]]) is true
    assert _POSform([x], [[0]], [[1]]) is true
    assert _POSform([x], [], []) is false

    # check working of simplify
    assert ((a & b) | (a & c)).simplify() == a & (b | c)
    assert (x & ~x).simplify() is false
    assert (x | ~x).simplify() is true


def test_bool_symbol():
    assert And(a, True) == a
    assert And(a, True, True) == a
    assert And(a, False) is false
    assert And(a, True, False) is false
    assert Or(a, True) is true
    assert Or(a, False) == a


def test_is_boolean():
    assert true.is_Boolean
    assert (a & b).is_Boolean
    assert (a | b).is_Boolean
    assert (~a).is_Boolean
    assert (a ^ b).is_Boolean


def test_subs():
    assert (a & b).subs({a: True}) == b
    assert (a & b).subs({a: False}) is false
    assert (a & b).subs({b: True}) == a
    assert (a & b).subs({b: False}) is false
    assert (a & b).subs({a: True, b: True}) is true
    assert (a | b).subs({a: True}) is true
    assert (a | b).subs({a: False}) == b
    assert (a | b).subs({b: True}) is true
    assert (a | b).subs({b: False}) == a
    assert (a | b).subs({a: True, b: True}) is true


def test_logic_commutatity():
    assert a & b == b & a
    assert a | b == b | a


def test_logic_associativity():
    assert (a & b) & c == a & (b & c)
    assert (a | b) | c == a | (b | c)


def test_double_negation():
    assert ~(~a) == a


def test_to_nnf():
    assert to_nnf(true) is true
    assert to_nnf(false) is false
    assert to_nnf(a) == a
    assert to_nnf(~a) == ~a

    class Foo(BooleanFunction):
        pass

    pytest.raises(ValueError, lambda: to_nnf(~Foo(a)))

    assert to_nnf(a | ~a | b) is true
    assert to_nnf(a & ~a & b) is false
    assert to_nnf(a >> b) == ~a | b
    assert to_nnf(Implies(a, b, evaluate=False)) == ~a | b
    assert to_nnf(a >> (c >> ~b)) == (~b | ~c) | ~a
    assert to_nnf(Equivalent(a, b)) == (a | ~b) & (b | ~a)
    assert to_nnf(Equivalent(a, b, c)) == (~a | b) & (~b | c) & (~c | a)
    assert to_nnf(Equivalent(a, b, c, d)) == (~a | b) & (~b | c) & (~c | d) & (~d | a)
    assert to_nnf(a ^ b ^ c) == (a | b | c) & (~a | ~b | c) & (a | ~b | ~c) & (~a | b | ~c)
    assert to_nnf(ITE(a, b, c)) == (~a | b) & (a | c)
    assert to_nnf(~(a | b | c)) == ~a & ~b & ~c
    assert to_nnf(~(a & b & c)) == ~a | ~b | ~c
    assert to_nnf(~(a >> b)) == a & ~b
    assert to_nnf(~(Equivalent(a, b, c))) == (a | b | c) & (~a | ~b | ~c)
    assert to_nnf(~(a ^ b ^ c)) == (~a | b | c) & (a | ~b | c) & (a | b | ~c) & (~a | ~b | ~c)
    assert to_nnf(~(ITE(a, b, c))) == (~a | ~b) & (a | ~c)
    assert to_nnf((a >> b) ^ (b >> a)) == (a & ~b) | (~a & b)
    assert to_nnf((a >> b) ^ (b >> a), False) == (~a | ~b | a | b) & ((a & ~b) | (~a & b))


def test_to_cnf():
    assert to_cnf(~(b | c)) == ~b & ~c
    assert to_cnf((a & b) | c) == (a | c) & (b | c)
    assert to_cnf(a >> b) == ~a | b
    assert to_cnf(a >> (b & c)) == (~a | b) & (~a | c)
    assert to_cnf(a & (b | c) | ~a & (b | c), True) == b | c
    assert to_cnf(a | (~b & ~c)) == (a | ~b) & (a | ~c)

    assert to_cnf(Equivalent(a, b)) == (a | ~b) & (b | ~a)
    assert to_cnf(Equivalent(a, b & c)) == (~a | b) & (~a | c) & (~b | ~c | a)
    assert to_cnf(Equivalent(a, b | c), True) == (~b | a) & (~c | a) & (b | c | ~a)
    assert to_cnf(~(a | b) | c) == (~a | c) & (~b | c)
    assert to_cnf((a & b) | c) == (a | c) & (b | c)


def test_to_dnf():
    assert to_dnf(true) == true
    assert to_dnf(~b & ~c) == ~b & ~c
    assert to_dnf(~(b | c)) == ~b & ~c
    assert to_dnf(a & (b | c)) == (a & b) | (a & c)
    assert to_dnf((~a | b) & c) == (b & c) | (c & ~a)
    assert to_dnf(a >> b) == ~a | b
    assert to_dnf(a >> (b & c)) == (~a) | (b & c)

    assert to_dnf(Equivalent(a, b), True) == (a & b) | (~a & ~b)
    assert to_dnf(Equivalent(a, b & c), True) == ((a & b & c) | (~a & ~b) | (~a & ~c))


def test_to_int_repr():
    assert to_int_repr([x | y, z | x], [x, y, z]) == [{1, 2}, {1, 3}]
    assert to_int_repr([x | y, z | ~x], [x, y, z]) == [{1, 2}, {-1, 3}]


def test_is_nnf():
    assert is_nnf(true) is True
    assert is_nnf(a) is True
    assert is_nnf(~a) is True
    assert is_nnf(a & b) is True
    assert is_nnf((a & b) | (~a & a) | (~b & b) | (~a & ~b), False) is True
    assert is_nnf((a | b) & (~a | ~b)) is True
    assert is_nnf(~(a | b)) is False
    assert is_nnf(a ^ b) is False
    assert is_nnf((a & b) | (~a & a) | (~b & b) | (~a & ~b), True) is False


def test_is_cnf():
    assert is_cnf(x) is True
    assert is_cnf(x | y | z) is True
    assert is_cnf(x & y & z) is True
    assert is_cnf((x | y) & z) is True
    assert is_cnf((x & y) | z) is False
    assert is_cnf(x & y & (z | ~(x ^ y))) is False
    assert is_cnf(x & y & (z | (x ^ y))) is False


def test_is_dnf():
    assert is_dnf(x) is True
    assert is_dnf(x | y | z) is True
    assert is_dnf(x & y & z) is True
    assert is_dnf((x & y) | z) is True
    assert is_dnf((x | y) & z) is False


def test_ITE():
    pytest.raises(ValueError, lambda: ITE(a, b))

    assert ITE(True, False, True) is false
    assert ITE(True, True, False) is true
    assert ITE(False, True, False) is false
    assert ITE(False, False, True) is true
    assert isinstance(ITE(a, b, c), ITE)

    assert ITE(True, b, c) == b
    assert ITE(False, b, c) == c

    assert ITE(a, b, b) == b

    assert ITE(c, False, True) == ~c
    assert ITE(c, True, False) == c

    assert ITE(a > 0, a**2, a).diff(a) == ITE(a > 0, 2*a, 1)


def test_is_literal():
    assert is_literal(True) is True
    assert is_literal(False) is True
    assert is_literal(a) is True
    assert is_literal(~a) is True
    assert is_literal(a | b) is False


def test_logic_operators():
    assert (True & a) == (a & True) == a
    assert (False & a) == (a & False) == false
    assert (True | a) == (a | True) == true
    assert (False | a) == (a | False) == a
    assert (True >> a) == (a << True) == a
    assert (False >> a) == (a << False) == true
    assert (a >> True) == (True << a) == true
    assert (a >> False) == (False << a) == ~a
    assert (True ^ a) == (a ^ True) == ~a
    assert (False ^ a) == (a ^ False) == a


def test_true_false():
    # pylint: disable=singleton-comparison,comparison-with-itself,unneeded-not
    assert true is true
    assert false is false
    assert true is not True
    assert false is not False
    assert true
    assert not false
    assert true == True  # noqa: E712
    assert false == False  # noqa: E712
    assert not true == False  # noqa: E712
    assert not false == True  # noqa: E712
    assert not true == false

    assert hash(true) == hash(True)
    assert hash(false) == hash(False)
    assert len({true, True}) == len({false, False}) == 1

    assert int(true) == 1
    assert int(false) == 0

    assert isinstance(true, BooleanAtom)
    assert isinstance(false, BooleanAtom)

    assert not isinstance(true, bool)
    assert not isinstance(false, bool)

    assert ~true is false
    assert Not(True) is false
    assert ~false is true
    assert Not(False) is true

    for T, F in itertools.product([True, true], [False, false]):
        assert And(T, F) is false
        assert And(F, T) is false
        assert And(F, F) is false
        assert And(T, T) is true
        assert And(T, x) == x
        assert And(F, x) is false
        if not (T is True and F is False):
            assert T & F is false
            assert F & T is false
        if F is not False:
            assert F & F is false
        if T is not True:
            assert T & T is true

        assert Or(T, F) is true
        assert Or(F, T) is true
        assert Or(F, F) is false
        assert Or(T, T) is true
        assert Or(T, x) is true
        assert Or(F, x) == x
        if not (T is True and F is False):
            assert T | F is true
            assert F | T is true
        if F is not False:
            assert F | F is false
        if T is not True:
            assert T | T is true

        assert Xor(T, F) is true
        assert Xor(F, T) is true
        assert Xor(F, F) is false
        assert Xor(T, T) is false
        assert Xor(T, x) == ~x
        assert Xor(F, x) == x
        if not (T is True and F is False):
            assert T ^ F is true
            assert F ^ T is true
        if F is not False:
            assert F ^ F is false
        if T is not True:
            assert T ^ T is false

        assert Nand(T, F) is true
        assert Nand(F, T) is true
        assert Nand(F, F) is true
        assert Nand(T, T) is false
        assert Nand(T, x) == ~x
        assert Nand(F, x) is true

        assert Nor(T, F) is false
        assert Nor(F, T) is false
        assert Nor(F, F) is true
        assert Nor(T, T) is false
        assert Nor(T, x) is false
        assert Nor(F, x) == ~x

        assert Implies(T, F) is false
        assert Implies(F, T) is true
        assert Implies(F, F) is true
        assert Implies(T, T) is true
        assert Implies(T, x) == x
        assert Implies(F, x) is true
        assert Implies(x, T) is true
        assert Implies(x, F) == ~x
        if not (T is True and F is False):
            assert T >> F is false
            assert F << T is false
            assert F >> T is true
            assert T << F is true
        if F is not False:
            assert F >> F is true
            assert F << F is true
        if T is not True:
            assert T >> T is true
            assert T << T is true

        assert Equivalent(T, F) is false
        assert Equivalent(F, T) is false
        assert Equivalent(F, F) is true
        assert Equivalent(T, T) is true
        assert Equivalent(T, x) == x
        assert Equivalent(F, x) == ~x
        assert Equivalent(x, T) == x
        assert Equivalent(x, F) == ~x

        assert ITE(T, T, T) is true
        assert ITE(T, T, F) is true
        assert ITE(T, F, T) is false
        assert ITE(T, F, F) is false
        assert ITE(F, T, T) is true
        assert ITE(F, T, F) is false
        assert ITE(F, F, T) is true
        assert ITE(F, F, F) is false


def test_bool_as_set():
    assert ((x <= 2) & (x >= -2)).as_set() == Interval(-2, 2)
    assert ((x >= 2) | (x <= -2)).as_set() == (Interval(-oo, -2) + Interval(2, oo, False))
    assert Not(x > 2, evaluate=False).as_set() == Interval(-oo, 2, True)
    assert true.as_set() == S.UniversalSet
    assert false.as_set() == EmptySet()


@pytest.mark.xfail
def test_multivariate_bool_as_set():
    ((x >= 0) & (y >= 0)).as_set()  # == Interval(0, oo)*Interval(0, oo)
    ((x >= 0) | (y >= 0)).as_set()  # == Reals*Reals - Interval(-oo, 0, True, True)*Interval(-oo, 0, True, True)


def test_all_or_nothing():
    args = x >= -oo, x <= oo
    v = And(*args)
    if isinstance(v, And):
        assert len(v.args) == len(args) - args.count(true)
    else:
        assert v
    v = Or(*args)
    if isinstance(v, Or):
        assert len(v.args) == 2
    else:
        assert v


def test_canonical_atoms():
    assert true.canonical == true
    assert false.canonical == false


def test_sympyissue_8777():
    assert ((x > 2) & (x < oo)).as_set() == Interval(2, oo, True, True)
    assert ((x >= 1) & (x < oo)).as_set() == Interval(1, oo, False, True)
    assert (x < oo).as_set() == Interval(-oo, oo, False, True)
    assert (x > -oo).as_set() == Interval(-oo, oo, True)


def test_sympyissue_8975():
    assert (((-oo < x) & (x <= -2)) | ((2 <= x) & (x < oo))).as_set() == Interval(-oo, -2, True) + Interval(2, oo, False, True)


def test_sympyissue_10240():
    assert (~((x > 2) & (x < 3))).as_set() == Union(Interval(-oo, 2, True), Interval(3, oo, False, True))


def test_sympyissue_10641():
    assert str(((x < sqrt(3)) | x).evalf(2)) == 'x | (x < 1.7)'
    assert str(((x < sqrt(3)) & x).evalf(2)) == 'x & (x < 1.7)'


def test_sympyissue_12522():
    assert Equality(1, 1).simplify() is true
    assert true.simplify() is true
    assert false.simplify() is false
