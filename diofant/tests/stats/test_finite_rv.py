import pytest

from diofant import (And, Dict, Eq, FiniteSet, Integer, Matrix, Or, Rational,
                     Symbol, Tuple, binomial, cos, pi, simplify, sqrt, symbols,
                     sympify)
from diofant.abc import p, x
from diofant.stats import (Bernoulli, Binomial, Coin, Die, DiscreteUniform, E,
                           FiniteRV, Hypergeometric, P, Rademacher, cdf,
                           cmoment, correlation, covariance, density, moment,
                           pspace, sample, skewness, smoment, variance, where)
from diofant.stats.frv_types import DieDistribution


__all__ = ()


def BayesTest(A, B):
    assert P(A, B) == P(And(A, B)) / P(B)
    assert P(A, B) == P(B, A) * P(A) / P(B)


def test_discreteuniform():
    # Symbolic
    a, b, c = symbols('a b c')
    X = DiscreteUniform('X', [a, b, c])
    assert X.pspace.distribution.pdf(a) == Rational(1, 3)
    assert X.pspace.distribution.pdf(p) == 0

    assert E(X) == (a + b + c)/3
    assert simplify(variance(X)
                    - ((a**2 + b**2 + c**2)/3 - (a/3 + b/3 + c/3)**2)) == 0
    assert P(Eq(X, a)) == P(Eq(X, b)) == P(Eq(X, c)) == Rational(1, 3)

    Y = DiscreteUniform('Y', range(-5, 5))

    # Numeric
    assert E(Y) == -Rational(1, 2)
    assert variance(Y) == Rational(33, 4)

    for x in range(-5, 5):
        assert P(Eq(Y, x)) == Rational(1, 10)
        assert P(Y <= x) == Rational(x + 6, 10)
        assert P(Y >= x) == Rational(5 - x, 10)

    assert dict(density(Die('D', 6)).items()) == \
        dict(density(DiscreteUniform('U', range(1, 7))).items())


def test_dice():
    # TODO: Make iid method!
    X, Y, Z = Die('X', 6), Die('Y', 6), Die('Z', 6)
    a, b = symbols('a b')

    assert E(X) == 3 + Rational(1, 2)
    assert variance(X) == Rational(35, 12)
    assert E(X + Y) == 7
    assert E(X + X) == 7
    assert E(a*X + b) == a*E(X) + b
    assert variance(X + Y) == variance(X) + variance(Y) == cmoment(X + Y, 2)
    assert variance(X + X) == 4 * variance(X) == cmoment(X + X, 2)
    assert cmoment(X, 0) == 1
    assert cmoment(4*X, 3) == 64*cmoment(X, 3)
    assert covariance(X, Y) == 0
    assert covariance(X, X + Y) == variance(X)
    assert density(Eq(cos(X*pi), 1))[True] == Rational(1, 2)
    assert correlation(X, Y) == 0
    assert correlation(X, Y) == correlation(Y, X)
    assert smoment(X + Y, 3) == skewness(X + Y)
    assert smoment(X, 0) == 1
    assert P(X > 3) == Rational(1, 2)
    assert P(2*X > 6) == Rational(1, 2)
    assert P(X > Y) == Rational(5, 12)
    assert P(Eq(X, Y)) == P(Eq(X, 1))

    assert E(X, X > 3) == 5 == moment(X, 1, 0, X > 3)
    assert E(X, Y > 3) == E(X) == moment(X, 1, 0, Y > 3)
    assert E(X + Y, Eq(X, Y)) == E(2*X)
    assert moment(X, 0) == 1
    assert moment(5*X, 2) == 25*moment(X, 2)

    assert P(X > 3, X > 3) == 1
    assert P(X > Y, Eq(Y, 6)) == 0
    assert P(Eq(X + Y, 12)) == Rational(1, 36)
    assert P(Eq(X + Y, 12), Eq(X, 6)) == Rational(1, 6)

    assert density(X + Y) == density(Y + Z) != density(X + X)
    d = density(2*X + Y**Z)
    assert d[22] == Rational(1, 108) and d[4100] == Rational(1, 216) and 3130 not in d

    assert pspace(X).domain.as_boolean() == Or(
        *[Eq(X.symbol, i) for i in [1, 2, 3, 4, 5, 6]])

    assert where(X > 3).set == FiniteSet(4, 5, 6)

    X = Die('X', 2)
    x = X.symbol

    assert X.pspace.compute_cdf(X) == {1: Rational(1, 2), 2: 1}
    assert X.pspace.sorted_cdf(X) == [(1, Rational(1, 2)), (2, 1)]

    assert X.pspace.compute_density(X)(1) == Rational(1, 2)
    assert X.pspace.compute_density(X)(0) == 0
    assert X.pspace.compute_density(X)(8) == 0

    assert X.pspace.density == x


def test_given():
    X = Die('X', 6)
    assert density(X, X > 5) == {6: 1}
    assert where(X > 2, X > 5).as_boolean() == Eq(X.symbol, 6)
    assert sample(X, X > 5) == 6


def test_domains():
    X, Y = Die('x', 6), Die('y', 6)
    x, y = X.symbol, Y.symbol
    # Domains
    d = where(X > Y)
    assert d.condition == (x > y)
    d = where(And(X > Y, Y > 3))
    assert d.as_boolean() == Or(And(Eq(x, 5), Eq(y, 4)), And(Eq(x, 6),
                                                             Eq(y, 5)), And(Eq(x, 6), Eq(y, 4)))
    assert len(d.elements) == 3

    assert len(pspace(X + Y).domain.elements) == 36

    Z = Die('x', 4)

    pytest.raises(ValueError, lambda: P(X > Z))  # Two domains with same internal symbol

    assert pspace(X + Y).domain.set == FiniteSet(1, 2, 3, 4, 5, 6)**2

    assert where(X > 3).set == FiniteSet(4, 5, 6)
    assert X.pspace.domain.dict == FiniteSet(
        *[Dict({X.symbol: i}) for i in range(1, 7)])

    assert where(X > Y).dict == FiniteSet(*[Dict({X.symbol: i, Y.symbol: j})
                                            for i in range(1, 7) for j in range(1, 7) if i > j])

    X, Y = Die('x', 7), Die('y', 3)
    x, y = X.symbol, Y.symbol
    pytest.raises(ValueError, lambda: X.pspace.where(Y < 3))
    cset = X.pspace.where(X < 3)
    assert ((x, 1),) in cset
    assert ((x, 3),) not in cset

    assert X.pspace.where(True) == X.pspace.domain


def test_dice_bayes():
    X, Y, Z = Die('X', 6), Die('Y', 6), Die('Z', 6)

    BayesTest(X > 3, X + Y < 5)
    BayesTest(Eq(X - Y, Z), Z > Y)
    BayesTest(X > 3, X > 2)


def test_die_args():
    pytest.raises(ValueError, lambda: Die('X', -1))  # issue sympy/sympy#8105: negative sides.
    pytest.raises(ValueError, lambda: Die('X', 0))
    pytest.raises(ValueError, lambda: Die('X', 1.5))  # issue sympy/sympy#8103: non integer sides.

    k = Symbol('k')
    sym_die = Die('X', k)
    pytest.raises(ValueError, lambda: density(sym_die).dict)


def test_bernoulli():
    p, a, b = symbols('p a b')
    X = Bernoulli('B', p, a, b)

    assert E(X) == a*p + b*(-p + 1)
    assert density(X)[a] == p
    assert density(X)[b] == 1 - p

    X = Bernoulli('B', p, 1, 0)

    assert E(X) == p
    assert simplify(variance(X)) == p*(1 - p)
    assert E(a*X + b) == a*E(X) + b
    assert simplify(variance(a*X + b)) == simplify(a**2 * variance(X))


def test_cdf():
    D = Die('D', 6)
    o = Integer(1)

    assert cdf(
        D) == sympify({1: o/6, 2: o/3, 3: o/2, 4: 2*o/3, 5: 5*o/6, 6: o})


def test_coins():
    C, D = Coin('C'), Coin('D')
    H, T = symbols('H, T')
    assert P(Eq(C, D)) == Rational(1, 2)
    assert density(Tuple(C, D)) == {(H, H): Rational(1, 4), (H, T): Rational(1, 4),
                                    (T, H): Rational(1, 4), (T, T): Rational(1, 4)}
    assert dict(density(C).items()) == {H: Rational(1, 2), T: Rational(1, 2)}

    F = Coin('F', Rational(1, 10))
    assert P(Eq(F, H)) == Rational(1, 10)

    d = pspace(C).domain

    assert d.as_boolean() == Or(Eq(C.symbol, H), Eq(C.symbol, T))

    pytest.raises(ValueError, lambda: P(C > D))  # Can't intelligently compare H to T


def test_binomial_verify_parameters():
    pytest.raises(ValueError, lambda: Binomial('b', .2, .5))
    pytest.raises(ValueError, lambda: Binomial('b', 3, 1.5))


def test_binomial_numeric():
    nvals = range(5)
    pvals = [0, Rational(1, 4), Rational(1, 2), Rational(3, 4), 1]

    for n in nvals:
        for p in pvals:
            X = Binomial('X', n, p)
            assert E(X) == n*p
            assert variance(X) == n*p*(1 - p)
            if n > 0 and 0 < p < 1:
                assert skewness(X) == (1 - 2*p)/sqrt(n*p*(1 - p))
            for k in range(n + 1):
                assert P(Eq(X, k)) == binomial(n, k)*p**k*(1 - p)**(n - k)


@pytest.mark.slow
def test_binomial_symbolic():
    n = 10  # Because we're using for loops, can't do symbolic n
    p = symbols('p', positive=True)
    X = Binomial('X', n, p)
    assert simplify(E(X)) == n*p == simplify(moment(X, 1))
    assert simplify(variance(X)) == n*p*(1 - p) == simplify(cmoment(X, 2))
    assert simplify(skewness(X) - (1-2*p)/sqrt(n*p*(1-p))) == 0

    # Test ability to change success/failure winnings
    H, T = symbols('H T')
    Y = Binomial('Y', n, p, succ=H, fail=T)
    assert simplify(E(Y) - (n*(H*p + T*(1 - p)))) == 0


def test_hypergeometric_numeric():
    for N in range(1, 5):
        for m in range(N + 1):
            for n in range(1, N + 1):
                X = Hypergeometric('X', N, m, n)
                N, m, n = map(sympify, (N, m, n))
                assert sum(density(X).values()) == 1
                assert E(X) == n * m / N
                if N > 1:
                    assert variance(X) == n*(m/N)*(N - m)/N*(N - n)/(N - 1)
                # Only test for skewness when defined
                if N > 2 and 0 < m < N and n < N:
                    assert skewness(X) == simplify((N - 2*m)*sqrt(N - 1)*(N - 2*n)
                                                   / (sqrt(n*m*(N - m)*(N - n))*(N - 2)))


def test_rademacher():
    X = Rademacher('X')

    assert E(X) == 0
    assert variance(X) == 1
    assert density(X)[-1] == Rational(1, 2)
    assert density(X)[1] == Rational(1, 2)


def test_FiniteRV():
    F = FiniteRV('F', {1: Rational(1, 2), 2: Rational(1, 4), 3: Rational(1, 4)})

    assert dict(density(F).items()) == {1: Rational(1, 2), 2: Rational(1, 4), 3: Rational(1, 4)}
    assert P(F >= 2) == Rational(1, 2)

    assert pspace(F).domain.as_boolean() == Or(
        *[Eq(F.symbol, i) for i in [1, 2, 3]])


def test_density_call():
    x = Bernoulli('x', p)
    d = density(x)
    assert d(0) == 1 - p
    assert d(5) == 0

    assert 0 in d
    assert 5 not in d
    assert d(0) == d[0]


def test_DieDistribution():
    X = DieDistribution(6)
    assert X.pdf(Rational(1, 2)) == 0
    assert X.pdf(x).subs({x: 1}).doit() == Rational(1, 6)
    assert X.pdf(x).subs({x: 7}).doit() == 0
    assert X.pdf(x).subs({x: -1}).doit() == 0
    assert X.pdf(x).subs({x: Rational(1, 3)}).doit() == 0
    pytest.raises(TypeError, lambda: X.pdf(x).subs({x: Matrix([0, 0])}))
    pytest.raises(ValueError, lambda: X.pdf(x**2 - 1))
