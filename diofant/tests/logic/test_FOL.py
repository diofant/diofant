"""Tests for the First Order Logic module."""

import operator

import pytest

from diofant import And, Or, Symbol, to_cnf, to_dnf, true
from diofant.abc import w, x, y, z
from diofant.logic.FOL import (FOL_KB, AppliedPredicate,
                               AppliedUndefinedBooleanFunction, Constant,
                               Exists, ForAll, Predicate, Quantifier,
                               UndefinedBooleanFunction, entails, fol_true,
                               mgu, resolve, standardize, to_pnf, to_snf)


def test_Predicate():
    a = Predicate('a')
    b = Predicate('b')
    p = a(x, y)
    q = b(z)
    assert isinstance(a, Predicate)
    assert isinstance(p, AppliedPredicate)
    assert p.name == 'a'
    assert p.func == a
    pytest.raises(ValueError, a)
    assert p != a(y, x)
    assert p != a(x, y, z)
    assert p & q == And(p, b(z), evaluate=False)
    assert Or(p) == p


def test_UndefinedBooleanFunction():
    f = UndefinedBooleanFunction('f')
    g = UndefinedBooleanFunction('g')
    fx = f(x)
    gyz = g(y, z)
    assert isinstance(f, UndefinedBooleanFunction)
    assert isinstance(fx, AppliedUndefinedBooleanFunction)
    assert fx.name == 'f'
    assert fx.func == f
    pytest.raises(ValueError, f)
    assert f(x, y) != f(y, x)
    assert f(x, y) != f(x, y, z)
    assert fx | gyz == Or(f(x), g(y, z), evaluate=False)
    assert And(fx) == f(x)


def test_Constant():
    a = Constant('a')
    b = Constant('b')
    assert a == Constant('a')
    assert a != b
    assert a.name == 'a'
    assert Constant(a) == a


def test_ForAll():
    a = Predicate('a')
    b = Predicate('b')
    assert ForAll(x, ForAll(y, a(x) & b(y))) == ForAll((x, y), a(x) & b(y))
    pytest.raises(ValueError, lambda: ForAll(x, a(x) | ForAll(x, b(x))))
    assert ForAll((), a(x)) == a(x)
    assert ForAll(y, a(x)) == a(x)
    assert ForAll((x, y), a(x)) == ForAll(x, a(x))
    assert ForAll((x, y, z), a(x) >> b(y, z)).vars == (x, y, z)
    assert ForAll((x, y, z), a(x) >> b(y, z)).expr == a(x) >> b(y, z)
    assert ForAll((x, y), ForAll(z, a(x, y) ^ b(y, z))).vars == (x, y, z)
    assert ForAll((x, y), ForAll(z, a(x, y) ^ b(y, z))).expr == a(x, y) ^ b(y, z)


def test_Exists():
    a = Predicate('a')
    b = Predicate('b')
    assert Exists(x, Exists(y, a(x) & b(y))) == Exists((x, y), a(x) & b(y))
    pytest.raises(ValueError, lambda: Exists(x, a(x) | Exists(x, b(x))))
    assert Exists((), a(x)) == a(x)
    assert Exists(y, a(x)) == a(x)
    assert Exists((x, y), a(x)) == Exists(x, a(x))
    assert Exists((x, y, z), a(x) >> b(y, z)).vars == (x, y, z)
    assert Exists((x, y, z), a(x) >> b(y, z)).expr == a(x) >> b(y, z)
    assert Exists((x, y), Exists(z, a(x, y) ^ b(y, z))).vars == (x, y, z)
    assert Exists((x, y), Exists(z, a(x, y) ^ b(y, z))).expr == a(x, y) ^ b(y, z)
    assert str(Exists((x, y), a(x) >> b(y, z))) == 'Exists((x, y), Implies(a(x), b(y, z)))'


def test_fol_true():
    p = Predicate('p')
    f = UndefinedBooleanFunction('f')
    _p = {0: False, 'default': True}
    _f = {0: 1, 'default': 0}
    _fn = {0: None, 'default': 0}
    assert fol_true(p(0)) is None
    assert fol_true(p(f(0)), {p: _p, f: _f}) is True
    assert fol_true(p(f(0)), {p: _p, f: _fn}) is True
    assert fol_true(p(f(x)), {x: 1, p: _p, f: _f}) is False
    assert fol_true(ForAll(x, p(x)), {p: _p}) is None
    assert fol_true(ForAll(x, Exists(y, p(x) & p(y))),
                    {x: [0, 1], p: _p}) is None

    pytest.raises(ValueError, lambda: fol_true(p(0), {object(): _f}))
    pytest.raises(ValueError, lambda: fol_true(Quantifier(x, x)))

    gt = Predicate('gt')
    lt = Predicate('lt')
    eq = Predicate('eq')
    add1 = UndefinedBooleanFunction('add1')

    def _add1(x):
        return (x + 1) % 6

    domain = range(6)

    assert fol_true(lt(add1(x), x), {x: 1, lt: operator.lt, add1: _add1}) is False
    assert fol_true(gt(x, y) & gt(y, z) >> gt(x, z),
                    {x: 3, y: 2, z: 1, gt: operator.gt}) is True
    assert fol_true(ForAll(x, gt(add1(x), x)),
                    {x: domain, gt: operator.gt, add1: _add1}) is False
    assert fol_true(Exists(x, eq(_add1(x), x),),
                    {x: domain, eq: operator.eq, add1: _add1}) is False
    assert fol_true(ForAll((x, y, z), (gt(x, y) & gt(y, z)) >> gt(x, z)),
                    {x: domain, y: domain, z: domain, gt: operator.gt}) is True


def test_standardize():
    x0 = Symbol('x0')
    x1 = Symbol('x1')
    y0 = Symbol('y0')
    p = Predicate('p')
    q = Predicate('q')
    r = Predicate('r')
    assert standardize(ForAll(x, p(x)) >> ForAll(x, q(x))) == \
        ForAll(x, p(x)) >> ForAll(x0, q(x0))
    assert standardize(ForAll(x, p(x)) & ForAll(x, q(x))) == \
        ForAll(x, p(x)) & ForAll(x, q(x))
    assert standardize(Exists(x, p(x)) | Exists(x, q(x))) == \
        Exists(x, p(x)) | Exists(x, q(x))
    assert standardize(ForAll(x, p(x) & q(x)) | ForAll(x, p(x) >> q(x))) == \
        ForAll(x, p(x) & q(x)) | ForAll(x0, p(x0) >> q(x0))
    assert standardize(Exists(x, p(x) & q(x)) & Exists(x, p(x) >> q(x))) == \
        Exists(x, p(x) & q(x)) & Exists(x0, p(x0) >> q(x0))
    assert standardize(Exists(x, p(x) >> ForAll(y, r(y))) |
                       Exists((x, y), p(x, y))) == \
        Exists(x, p(x) >> ForAll(y0, r(y0))) | Exists((x, y), p(x, y))
    assert standardize(ForAll(x, Exists(y, p(x) & q(y))) | Exists(x, p(x)) &
                       ForAll(x, q(x))) == ((Exists(x0, p(x0)) &
                                             ForAll(x, q(x))) |
                                            ForAll(x1, Exists(y, p(x1) & q(y))))


def test_to_pnf():
    p = Predicate('p')
    q = Predicate('q')
    assert to_pnf(True) is true
    assert to_pnf(p(x, y)) == p(x, y)
    assert to_pnf(~ForAll(x, p(x))) == Exists(x, ~p(x))
    assert to_pnf(~Exists(x, p(x))) == ForAll(x, ~p(x))
    assert to_pnf(Exists(x, ForAll(y, p(x, y)))) == \
        Exists(x, ForAll(y, p(x, y)))
    assert to_pnf(ForAll((x, y), ~(p(x) >> q(y)))) == \
        ForAll((x, y), p(x) & ~q(y))
    assert to_pnf(Exists(x, p(x)) & Exists(y, q(y))) == \
        Exists((x, y), p(x) & q(y))
    assert to_pnf(ForAll(x, p(x)) >> ForAll(y, q(y))) == \
        ForAll(y, Exists(x, ~p(x) | q(y)))
    assert to_pnf(ForAll(x, p(x)) & ForAll(x, q(x))) == \
        ForAll(x, p(x) & q(x))
    assert to_pnf(Exists(x, p(x)) | Exists(x, q(x))) == \
        Exists(x, p(x) | q(x))
    assert to_pnf(ForAll(x, p(x)) >> ForAll(x, q(x)), iter([y])) == \
        ForAll((x), Exists((y), ~p(y) | q(x)))

    pytest.raises(ValueError, lambda: to_pnf(~Quantifier(x, x)))


def test_to_snf():
    p = Predicate('p')
    q = Predicate('q')
    f0 = UndefinedBooleanFunction('f0')
    f1 = UndefinedBooleanFunction('f1')
    f = UndefinedBooleanFunction('f')
    c0 = Constant('c0')
    c = Constant('c')
    assert to_snf(ForAll((x, y), p(x) >> q(y))) == ~p(x) | q(y)
    assert to_snf(Exists(x, ForAll(y, p(x) | q(y)))) == p(c0) | q(y)
    assert to_snf(ForAll(x, Exists(y, p(x) & q(y)))) == p(x) & q(f0(x))
    assert to_snf(ForAll(w, Exists(x, ForAll(y, Exists(z, p(w, x) >> q(y, z)))))) == \
        ~p(w, f0(w)) | q(y, f1(w, y))
    assert to_snf(Exists(x, ForAll(y, Exists(z, p(x) | p(y) | p(z)))),
                  functions=iter([f]), constants=iter([c])) == p(c) | p(y) | p(f(y))

    pytest.raises(ValueError, lambda: to_snf(Quantifier(x, x)))


def test_to_cnf():
    p = Predicate('p')
    q = Predicate('q')
    assert to_cnf(p(x) | q(x) | p(y) | q(y)) == p(x) | q(x) | p(y) | q(y)
    assert to_cnf((p(x) & q(x)) | (p(y) & q(y))) == \
        (p(x) | p(y)) & (p(x) | q(y)) & (q(x) | p(y)) & (q(x) | q(y))


def test_to_dnf():
    p = Predicate('p')
    q = Predicate('q')
    assert to_dnf(p(x) & q(x) & p(y) & q(y)) == p(x) & q(x) & p(y) & q(y)
    assert to_dnf((p(x) | q(x)) & (p(y) | q(y))) == \
        (p(x) & p(y)) | (p(x) & q(y)) | (q(x) & p(y)) | (q(x) & q(y))


def test_mgu():
    p = Predicate('p')
    q = Predicate('q')
    f = UndefinedBooleanFunction('f')
    g = UndefinedBooleanFunction('g')
    a = Constant('a')
    b = Constant('b')
    assert mgu(p(x), q(x)) is False
    assert mgu(p(x), p(x, y)) is False
    assert mgu(p(x), p(x)) == {true: true}
    assert mgu(p(x, x), p(y, f(y))) is False
    assert mgu(p(x, y), p(f(x), z)) is False
    assert mgu(p(f(a)), p(f(b))) is False
    assert mgu(p(x, f(a)), p(y, x)) == {x: f(a), y: f(a)}
    assert mgu(p(a, x, f(g(z))), p(z, f(y), f(y))) == {x: f(g(a)), z: a, y: g(a)}


def test_resolution():
    p = Predicate('p')
    q = Predicate('q')
    a = Constant('a')
    assert resolve((p(x) | q(x)) & (~p(x) | q(x))) is True
    assert resolve(p(x) >> q(x) & p(a), ~q(a)) is False
    pytest.raises(ValueError, lambda: resolve(p(x) | a))


def test_entails():
    parent = Predicate('parent')
    male = Predicate('male')
    father = Predicate('father')
    sibling = Predicate('sibling')
    tom = Constant('tom')
    john = Constant('john')
    fred = Constant('fred')
    clauses = [
        parent(x, y) & male(x) >> father(x, y),
        father(x, y) & father(x, z) >> sibling(y, z),
        parent(tom, john),
        male(tom),
        parent(tom, fred)
    ]
    assert entails(sibling(john, fred), clauses) is True
    assert entails(~sibling(john, fred), clauses) is False

    dog = Predicate('dog')
    animal_lover = Predicate('animal_lover')
    kills = Predicate('kills')
    owns = Predicate('owns')
    cat = Predicate('cat')
    animal = Predicate('animal')
    d = Constant('d')
    jack = Constant('jack')
    tuna = Constant('tuna')
    curiosity = Constant('curiosity')
    clauses = [
        ForAll((x, y), (dog(y) & owns(x, y)) >> animal_lover(x)),
        ForAll((x, y), (animal_lover(x) & animal(y)) >> ~kills(x, y)),
        ForAll(x, cat(x) >> animal(x)),
        ForAll(x, cat(x) >> (kills(jack, x) | kills(curiosity, x))),
        dog(d),
        cat(tuna),
        owns(jack, d)
    ]
    assert entails(kills(jack, tuna), clauses) is False
    assert entails(kills(curiosity, tuna), clauses) is True


def test_FOL_KB():
    kb = FOL_KB()
    man = Predicate('man')
    mortal = Predicate('mortal')
    socrates = Constant('socrates')
    kb.tell(man(x) >> mortal(x))
    kb.tell(man(socrates))
    assert kb.ask(mortal(socrates))
    assert not kb.ask(~mortal(socrates))

    kb = FOL_KB()
    dangerous = Predicate('dangerous')
    predator = Predicate('predator')
    lion = Predicate('lion')
    human = Predicate('human')
    jack = Constant('jack')
    leo = Constant('leo')
    simba = Constant('simba')
    kb.tell((predator(x) & human(y)) >> dangerous(x, y))
    kb.tell(lion(x) >> predator(x))
    kb.tell(human(jack))
    kb.tell(lion(leo))
    kb.tell(lion(simba))
    assert kb.ask(predator(simba))
    assert kb.ask(dangerous(leo, jack))
    assert not kb.ask(dangerous(leo, simba))
    assert kb.ask(lion(x), all_answers=True) == [lion(leo), lion(simba)]

    kb = FOL_KB()
    knows = Predicate('knows')
    p = Constant('p')
    q = Constant('q')
    r = Constant('r')
    s = Constant('s')
    kb.tell((knows(x, y) & knows(y, z)) >> knows(x, z))
    kb.tell(knows(p, q))
    kb.tell(knows(q, r))
    assert kb.ask(knows(p, q))
    assert kb.ask(knows(q, r))
    assert kb.ask(knows(p, r))
    assert kb.ask(knows(x, y), all_answers=True) == \
        [knows(p, q), knows(p, r), knows(q, r)]
    kb.tell(knows(x, y) >> knows(y, x))
    assert kb.ask(knows(r, p))
    assert kb.ask(knows(x, y), all_answers=True) == \
        [knows(p, p), knows(p, q), knows(p, r), knows(q, p), knows(q, q),
         knows(q, r), knows(r, p), knows(r, q), knows(r, r)]
    kb.tell(knows(r, s))
    assert kb.ask(knows(p, s))
    assert kb.ask(knows(s, p))

    kb = FOL_KB()
    p = Predicate('p')
    q = Predicate('q')
    pytest.raises(ValueError, lambda: kb.tell(p(x) >> ~q(x)))
    pytest.raises(ValueError, lambda: kb.tell(~p(x) >> q(x)))
    pytest.raises(ValueError, lambda: kb.tell(~p(x)))

    kb = FOL_KB()
    r = Predicate('r')
    pytest.raises(ValueError, lambda: kb.tell((p(x) | q(x)) >> r(x)))
