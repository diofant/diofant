""" Test file for SymPy First Order Logic """

import operator

import pytest

from diofant.abc import W, X, Y, Z
from diofant.core import Symbol
from diofant.logic.boolalg import (And, Implies, Not, Or, Xor, to_cnf, to_dnf,
                                   true)
from diofant.logic.FOL import (FOL_KB, AppliedFunction, AppliedPredicate,
                               Constant, Exists, ForAll, Function, Predicate,
                               Quantifier, entails, fol_true, mgu, resolve,
                               standardize, to_pnf, to_snf)


def test_Predicate():
    A = Predicate('A')
    B = Predicate('B')
    P = A(X, Y)
    Q = B(Z)
    assert isinstance(A, Predicate)
    assert isinstance(P, AppliedPredicate)
    assert P.name == 'A'
    assert P.func == A
    pytest.raises(ValueError, lambda: A())
    assert not A(X, Y) == A(Y, X)
    assert not A(X, Y) == A(X, Y, Z)
    assert P & Q == And(A(X, Y), B(Z))
    assert Or(P) == A(X, Y)


def test_Function():
    f = Function('f')
    g = Function('g')
    fx = f(X)
    gyz = g(Y, Z)
    assert isinstance(f, Function)
    assert isinstance(fx, AppliedFunction)
    assert fx.name == 'f'
    assert fx.func == f
    pytest.raises(ValueError, lambda: f())
    assert not f(X, Y) == f(Y, X)
    assert not f(X, Y) == f(X, Y, Z)
    assert fx | gyz == Or(f(X), g(Y, Z))
    assert And(fx) == f(X)


def test_Constant():
    A = Constant('A')
    B = Constant('B')
    assert A == Constant('A')
    assert not A == B
    assert A.name == 'A'
    assert Constant(A) == A


def test_ForAll():
    A = Predicate('A')
    B = Predicate('B')
    assert ForAll(X, ForAll(Y, A(X) & B(Y))) == ForAll((X, Y), A(X) & B(Y))
    pytest.raises(ValueError, lambda: ForAll(X, A(X) | ForAll(X, B(X))))
    ForAll((), A(X)) == A(X)
    ForAll(Y, A(X)) == A(X)
    ForAll((X, Y), A(X)) == ForAll(X, A(X))
    ForAll((X, Y, Z), A(X) >> B(Y, Z)).vars == (X, Y, Z)
    ForAll((X, Y, Z), A(X) >> B(Y, Z)).expr == Implies(A(X), B(Y, Z))
    ForAll((X, Y), ForAll(Z, A(X, Y) ^ B(Y, Z))).vars == (X, Y, Z)
    ForAll((X, Y), ForAll(Z, A(X, Y) ^ B(Y, Z))).expr == Xor(A(X, Y), B(Y, Z))


def test_Exists():
    A = Predicate('A')
    B = Predicate('B')
    assert Exists(X, Exists(Y, A(X) & B(Y))) == Exists((X, Y), A(X) & B(Y))
    pytest.raises(ValueError, lambda: Exists(X, A(X) | Exists(X, B(X))))
    Exists((), A(X)) == A(X)
    Exists(Y, A(X)) == A(X)
    Exists((X, Y), A(X)) == Exists(X, A(X))
    Exists((X, Y, Z), A(X) >> B(Y, Z)).vars == (X, Y, Z)
    Exists((X, Y, Z), A(X) >> B(Y, Z)).expr == Implies(A(X), B(Y, Z))
    Exists((X, Y), Exists(Z, A(X, Y) ^ B(Y, Z))).vars == (X, Y, Z)
    Exists((X, Y), Exists(Z, A(X, Y) ^ B(Y, Z))).expr == Xor(A(X, Y), B(Y, Z))


def test_fol_true():
    P = Predicate('P')
    f = Function('f')
    _P = {0: False, 'default': True}
    _f = {0: 1, 'default': 0}
    _fn = {0: None, 'default': 0}
    assert fol_true(P(0)) is None
    assert fol_true(P(f(0)), {P: _P, f: _f}) is True
    assert fol_true(P(f(0)), {P: _P, f: _fn}) is True
    assert fol_true(P(f(X)), {X: 1, P: _P, f: _f}) is False
    assert fol_true(ForAll(X, P(X)), {P: _P}) is None
    assert fol_true(ForAll(X, Exists(Y, P(X) & P(Y))),
                    {X: [0, 1], P: _P}) is None

    pytest.raises(ValueError, lambda: fol_true(P(0), {object(): _f}))
    pytest.raises(ValueError, lambda: fol_true(Quantifier(X, X)))

    GT = Predicate('GT')
    LT = Predicate('LT')
    EQ = Predicate('EQ')
    add1 = Function('add1')

    def _add1(x):
        return (x + 1) % 6

    domain = range(6)
    assert fol_true(LT(add1(X), X), {X: 1, LT: operator.lt, add1: _add1}) is False
    assert fol_true(GT(X, Y) & GT(Y, Z) >> GT(X, Z),
                    {X: 3, Y: 2, Z: 1, GT: operator.gt}) is True
    assert fol_true(ForAll(X, GT(add1(X), X)),
                    {X: domain, GT: operator.gt, add1: _add1}) is False
    assert fol_true(Exists(X, EQ(_add1(X), X),),
                    {X: domain, EQ: operator.eq, add1: _add1}) is False
    assert fol_true(ForAll((X, Y, Z), (GT(X, Y) & GT(Y, Z)) >> GT(X, Z)),
                    {X: domain, Y: domain, Z: domain, GT: operator.gt}) is True


def test_standardize():
    X0 = Symbol('X0')
    X1 = Symbol('X1')
    Y0 = Symbol('Y0')
    P = Predicate('P')
    Q = Predicate('Q')
    R = Predicate('R')
    assert standardize(ForAll(X, P(X)) >> ForAll(X, Q(X))) == \
        ForAll(X, P(X)) >> ForAll(X0, Q(X0))
    assert standardize(ForAll(X, P(X)) & ForAll(X, Q(X))) == \
        ForAll(X, P(X)) & ForAll(X, Q(X))
    assert standardize(Exists(X, P(X)) | Exists(X, Q(X))) == \
        Exists(X, P(X)) | Exists(X, Q(X))
    assert standardize(ForAll(X, P(X) & Q(X)) | ForAll(X, P(X) >> Q(X))) == \
        ForAll(X, P(X) & Q(X)) | ForAll(X0, P(X0) >> Q(X0))
    assert standardize(Exists(X, P(X) & Q(X)) & Exists(X, P(X) >> Q(X))) == \
        Exists(X, P(X) & Q(X)) & Exists(X0, P(X0) >> Q(X0))
    assert standardize(Exists(X, P(X) >> ForAll(Y, R(Y))) |
                       Exists((X, Y), P(X, Y))) == \
        Exists(X, P(X) >> ForAll(Y0, R(Y0))) | Exists((X, Y), P(X, Y))
    assert standardize(ForAll(X, Exists(Y, P(X) & Q(Y))) | Exists(X, P(X)) &
                       ForAll(X, Q(X))) == Or(And(Exists((X0), P(X0)),
                                                  ForAll((X), Q(X))),
                                              ForAll((X1),
                                                     Exists((Y), And(P(X1),
                                                                     Q(Y)))))


def test_to_pnf():
    P = Predicate('P')
    Q = Predicate('Q')
    assert to_pnf(True) is true
    assert to_pnf(P(X, Y)) == P(X, Y)
    assert to_pnf(~ForAll(X, P(X))) == Exists(X, ~P(X))
    assert to_pnf(~Exists(X, P(X))) == ForAll(X, ~P(X))
    assert to_pnf(Exists(X, ForAll(Y, P(X, Y)))) == \
        Exists(X, ForAll(Y, P(X, Y)))
    assert to_pnf(ForAll((X, Y), ~(P(X) >> Q(Y)))) == \
        ForAll((X, Y), P(X) & ~Q(Y))
    assert to_pnf(Exists(X, P(X)) & Exists(Y, Q(Y))) == \
        Exists((X, Y), P(X) & Q(Y))
    assert to_pnf(ForAll(X, P(X)) >> ForAll(Y, Q(Y))) == \
        ForAll(Y, Exists(X, ~P(X) | Q(Y)))
    assert to_pnf(ForAll(X, P(X)) & ForAll(X, Q(X))) == \
        ForAll(X, P(X) & Q(X))
    assert to_pnf(Exists(X, P(X)) | Exists(X, Q(X))) == \
        Exists(X, P(X) | Q(X))
    assert to_pnf(ForAll(X, P(X)) >> ForAll(X, Q(X)), iter([Y])) == \
        ForAll((X), Exists((Y), ~P(Y) | Q(X)))

    pytest.raises(ValueError, lambda: to_pnf(Not(Quantifier(X, X))))


def test_to_snf():
    P = Predicate('P')
    Q = Predicate('Q')
    f0 = Function('f0')
    f1 = Function('f1')
    F = Function('F')
    c0 = Constant('c0')
    C = Constant('C')
    assert to_snf(ForAll((X, Y), P(X) >> Q(Y))) == ~P(X) | Q(Y)
    assert to_snf(Exists(X, ForAll(Y, P(X) | Q(Y)))) == P(c0) | Q(Y)
    assert to_snf(ForAll(X, Exists(Y, P(X) & Q(Y)))) == P(X) & Q(f0(X))
    assert to_snf(ForAll(W, Exists(X, ForAll(Y, Exists(Z, P(W, X) >> Q(Y, Z)))))) == \
        ~P(W, f0(W)) | Q(Y, f1(W, Y))
    assert to_snf(Exists(X, ForAll(Y, Exists(Z, P(X) | P(Y) | P(Z)))),
                  functions=iter([F]), constants=iter([C])) == P(C) | P(Y) | P(F(Y))

    pytest.raises(ValueError, lambda: to_snf(Quantifier(X, X)))


def test_to_cnf():
    P = Predicate('P')
    Q = Predicate('Q')
    assert to_cnf(P(X) | Q(X) | P(Y) | Q(Y)) == P(X) | Q(X) | P(Y) | Q(Y)
    assert to_cnf((P(X) & Q(X)) | (P(Y) & Q(Y))) == \
        (P(X) | P(Y)) & (P(X) | Q(Y)) & (Q(X) | P(Y)) & (Q(X) | Q(Y))


def test_to_dnf():
    P = Predicate('P')
    Q = Predicate('Q')
    assert to_dnf(P(X) & Q(X) & P(Y) & Q(Y)) == P(X) & Q(X) & P(Y) & Q(Y)
    assert to_dnf((P(X) | Q(X)) & (P(Y) | Q(Y))) == \
        (P(X) & P(Y)) | (P(X) & Q(Y)) | (Q(X) & P(Y)) | (Q(X) & Q(Y))


def test_mgu():
    P = Predicate('P')
    Q = Predicate('Q')
    f = Function('f')
    g = Function('g')
    a = Constant('a')
    b = Constant('b')
    assert mgu(P(X), Q(X)) is False
    assert mgu(P(X), P(X, Y)) is False
    assert mgu(P(X), P(X)) == {true: true}
    assert mgu(P(X, X), P(Y, f(Y))) is False
    assert mgu(P(X, Y), P(f(X), Z)) is False
    assert mgu(P(f(a)), P(f(b))) is False
    assert mgu(P(X, f(a)), P(Y, X)) == {X: f(a), Y: f(a)}
    assert mgu(P(a, X, f(g(Z))), P(Z, f(Y), f(Y))) == {X: f(g(a)), Z: a, Y: g(a)}


def test_resolution():
    P = Predicate('P')
    Q = Predicate('Q')
    a = Constant('a')
    assert resolve(And(Or(P(X), Q(X)), Or(~P(X), Q(X)))) is True
    assert resolve(And(P(X) >> Q(X), P(a), ~Q(a))) is False
    pytest.raises(ValueError, lambda: resolve(P(X) | a))


def test_entails():
    Parent = Predicate('Parent')
    Male = Predicate('Male')
    Father = Predicate('Father')
    Sibling = Predicate('Sibling')
    Tom = Constant('Tom')
    John = Constant('John')
    Fred = Constant('Fred')
    clauses = [
        Parent(X, Y) & Male(X) >> Father(X, Y),
        Father(X, Y) & Father(X, Z) >> Sibling(Y, Z),
        Parent(Tom, John),
        Male(Tom),
        Parent(Tom, Fred)
    ]
    assert entails(Sibling(John, Fred), clauses) is True
    assert entails(~Sibling(John, Fred), clauses) is False

    Dog = Predicate('Dog')
    AnimalLover = Predicate('AnimalLover')
    Kills = Predicate('Kills')
    Owns = Predicate('Owns')
    Cat = Predicate('Cat')
    Animal = Predicate('Animal')
    D = Constant('D')
    Jack = Constant('Jack')
    Tuna = Constant('Tuna')
    Curiosity = Constant('Curiosity')
    clauses = [
        ForAll((X, Y), (Dog(Y) & Owns(X, Y)) >> AnimalLover(X)),
        ForAll((X, Y), (AnimalLover(X) & Animal(Y)) >> ~Kills(X, Y)),
        ForAll(X, Cat(X) >> Animal(X)),
        ForAll(X, Cat(X) >> (Kills(Jack, X) | Kills(Curiosity, X))),
        Dog(D),
        Cat(Tuna),
        Owns(Jack, D)
    ]
    assert entails(Kills(Jack, Tuna), clauses) is False
    assert entails(Kills(Curiosity, Tuna), clauses) is True


def test_FOL_KB():
    KB = FOL_KB()
    Man = Predicate('Man')
    Mortal = Predicate('Mortal')
    Socrates = Constant('Socrates')
    KB.tell(Man(X) >> Mortal(X))
    KB.tell(Man(Socrates))
    assert KB.ask(Mortal(Socrates))
    assert not KB.ask(~Mortal(Socrates))

    KB = FOL_KB()
    Dangerous = Predicate('Dangerous')
    Predator = Predicate('Predator')
    Lion = Predicate('Lion')
    Human = Predicate('Human')
    Jack = Constant('Jack')
    Leo = Constant('Leo')
    Simba = Constant('Simba')
    KB.tell((Predator(X) & Human(Y)) >> Dangerous(X, Y))
    KB.tell(Lion(X) >> Predator(X))
    KB.tell(Human(Jack))
    KB.tell(Lion(Leo))
    KB.tell(Lion(Simba))
    assert KB.ask(Predator(Simba))
    assert KB.ask(Dangerous(Leo, Jack))
    assert not KB.ask(Dangerous(Leo, Simba))
    assert KB.ask(Lion(X), all_answers=True) == [Lion(Leo), Lion(Simba)]

    KB = FOL_KB()
    Knows = Predicate('Knows')
    P = Constant('P')
    Q = Constant('Q')
    R = Constant('R')
    S = Constant('S')
    KB.tell((Knows(X, Y) & Knows(Y, Z)) >> Knows(X, Z))
    KB.tell(Knows(P, Q))
    KB.tell(Knows(Q, R))
    assert KB.ask(Knows(P, Q))
    assert KB.ask(Knows(Q, R))
    assert KB.ask(Knows(P, R))
    assert KB.ask(Knows(X, Y), all_answers=True) == \
        [Knows(P, Q), Knows(P, R), Knows(Q, R)]
    KB.tell(Knows(X, Y) >> Knows(Y, X))
    assert KB.ask(Knows(R, P))
    assert KB.ask(Knows(X, Y), all_answers=True) == \
        [Knows(P, P), Knows(P, Q), Knows(P, R), Knows(Q, P), Knows(Q, Q),
         Knows(Q, R), Knows(R, P), Knows(R, Q), Knows(R, R)]
    KB.tell(Knows(R, S))
    assert KB.ask(Knows(P, S))
    assert KB.ask(Knows(S, P))

    KB = FOL_KB()
    P = Predicate('P')
    Q = Predicate('Q')
    pytest.raises(ValueError, lambda: KB.tell(P(X) >> ~Q(X)))
    pytest.raises(ValueError, lambda: KB.tell(~P(X) >> Q(X)))
    pytest.raises(ValueError, lambda: KB.tell(~P(X)))

    KB = FOL_KB()
    R = Predicate('R')
    pytest.raises(ValueError, lambda: KB.tell((P(X) | Q(X)) >> R(X)))
