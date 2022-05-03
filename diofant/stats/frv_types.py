"""
Finite Discrete Random Variables - Prebuilt variable types

Contains
========
FiniteRV
DiscreteUniform
Die
Bernoulli
Coin
Binomial
Hypergeometric
"""

from ..concrete import Sum
from ..core import Dict, Dummy, Expr, Integer, Rational, cacheit, sympify
from ..core.compatibility import as_int
from ..core.logic import fuzzy_and, fuzzy_not
from ..functions import KroneckerDelta, binomial
from .frv import SingleFiniteDistribution, SingleFinitePSpace


__all__ = ('FiniteRV', 'DiscreteUniform', 'Die', 'Bernoulli', 'Coin',
           'Binomial', 'Hypergeometric')


def rv(name, cls, *args):
    density = cls(*args)
    return SingleFinitePSpace(name, density).value


class FiniteDistributionHandmade(SingleFiniteDistribution):
    @property
    def dict(self):
        return self.args[0]

    def __new__(cls, density):
        density = Dict(density)
        return Expr.__new__(cls, density)


def FiniteRV(name, density):
    """
    Create a Finite Random Variable given a dict representing the density.

    Returns a RandomSymbol.

    >>> from diofant.stats import E, P

    >>> density = {0: .1, 1: .2, 2: .3, 3: .4}
    >>> X = FiniteRV('X', density)

    >>> E(X)
    2.00000000000000
    >>> P(X >= 2)
    0.700000000000000

    """
    return rv(name, FiniteDistributionHandmade, density)


class DiscreteUniformDistribution(SingleFiniteDistribution):
    @property
    def p(self):
        return Rational(1, len(self.args))

    @property  # type: ignore[misc]
    @cacheit
    def dict(self):
        return {k: self.p for k in self.set}

    @property
    def set(self):
        return self.args

    def pdf(self, x):  # pylint: disable=invalid-overridden-method
        if x in self.args:
            return self.p
        else:
            return Integer(0)


def DiscreteUniform(name, items):
    """
    Create a Finite Random Variable representing a uniform distribution over
    the input set.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> X = DiscreteUniform('X', (a, b, c))  # equally likely over a, b, c
    >>> density(X).dict
    {a: 1/3, b: 1/3, c: 1/3}

    >>> Y = DiscreteUniform('Y', list(range(5)))  # distribution over a range
    >>> density(Y).dict
    {0: 1/5, 1: 1/5, 2: 1/5, 3: 1/5, 4: 1/5}

    """
    return rv(name, DiscreteUniformDistribution, *items)


class DieDistribution(SingleFiniteDistribution):
    _argnames = 'sides',

    def __new__(cls, sides):
        sides_sym = sympify(sides)
        if fuzzy_not(fuzzy_and((sides_sym.is_integer, sides_sym.is_positive))):
            raise ValueError("'sides' must be a positive integer.")
        return super().__new__(cls, sides)

    @property  # type: ignore[misc]
    @cacheit
    def dict(self):
        as_int(self.sides)
        return super().dict

    @property
    def set(self):
        return list(map(Integer, range(1, self.sides + 1)))

    def pdf(self, x):  # pylint: disable=invalid-overridden-method
        x = sympify(x)
        if x.is_number:
            if x.is_Integer and 1 <= x <= self.sides:
                return Rational(1, self.sides)
            return Integer(0)
        if x.is_Symbol:
            i = Dummy('i', integer=True, positive=True)
            return Sum(KroneckerDelta(x, i)/self.sides, (i, 1, self.sides))
        raise ValueError("'x' expected as an argument of type 'number' or 'symbol', "
                         f'not {type(x)}')


def Die(name, sides=6):
    """
    Create a Finite Random Variable representing a fair die.

    Returns a RandomSymbol.

    >>> from diofant.stats import density

    >>> D6 = Die('D6', 6)  # Six sided Die
    >>> density(D6).dict
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}

    >>> D4 = Die('D4', 4)  # Four sided Die
    >>> density(D4).dict
    {1: 1/4, 2: 1/4, 3: 1/4, 4: 1/4}

    """
    return rv(name, DieDistribution, sides)


class BernoulliDistribution(SingleFiniteDistribution):
    _argnames = ('p', 'succ', 'fail')

    @property  # type: ignore[misc]
    @cacheit
    def dict(self):
        return {self.succ: self.p, self.fail: 1 - self.p}


def Bernoulli(name, p, succ=1, fail=0):
    """
    Create a Finite Random Variable representing a Bernoulli process.

    Returns a RandomSymbol

    >>> from diofant.stats import density

    >>> X = Bernoulli('X', Rational(3, 4))  # 1-0 Bernoulli variable, probability = 3/4
    >>> density(X).dict
    {0: 1/4, 1: 3/4}

    >>> X = Bernoulli('X', Rational(1, 2), 'Heads', 'Tails')  # A fair coin toss
    >>> density(X).dict
    {Heads: 1/2, Tails: 1/2}

    """
    return rv(name, BernoulliDistribution, p, succ, fail)


def Coin(name, p=Rational(1, 2)):
    """
    Create a Finite Random Variable representing a Coin toss.

    Probability p is the chance of getting "Heads." Half by default

    Returns a RandomSymbol.

    >>> from diofant.stats import density

    >>> H, T = Symbol('H'), Symbol('T')

    >>> C = Coin('C')  # A fair coin toss
    >>> density(C).dict
    {H: 1/2, T: 1/2}

    >>> C2 = Coin('C2', Rational(3, 5))  # An unfair coin
    >>> density(C2).dict
    {H: 3/5, T: 2/5}

    """
    return rv(name, BernoulliDistribution, p, 'H', 'T')


class BinomialDistribution(SingleFiniteDistribution):
    _argnames = ('n', 'p', 'succ', 'fail')

    def __new__(cls, *args):
        n = args[BinomialDistribution._argnames.index('n')]
        p = args[BinomialDistribution._argnames.index('p')]
        n_sym = sympify(n)
        p_sym = sympify(p)

        if fuzzy_not(fuzzy_and((n_sym.is_integer, n_sym.is_nonnegative))):
            raise ValueError(f"'n' must be positive integer. n = {n!s}.")
        if fuzzy_not(fuzzy_and((p_sym.is_nonnegative, (p_sym - 1).is_nonpositive))):
            raise ValueError(f"'p' must be: 0 <= p <= 1 . p = {p!s}")
        return super().__new__(cls, *args)

    @property
    def n(self):
        return self.args[self._argnames.index('n')]

    @property  # type: ignore[misc]
    @cacheit
    def dict(self):
        n, p, succ, fail = self.n, self.p, self.succ, self.fail
        n = as_int(n)
        return {k*succ + (n - k)*fail:
                binomial(n, k) * p**k * (1 - p)**(n - k) for k in range(n + 1)}


def Binomial(name, n, p, succ=1, fail=0):
    """
    Create a Finite Random Variable representing a binomial distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> X = Binomial('X', 4, Rational(1, 2))  # Four "coin flips"
    >>> density(X).dict
    {0: 1/16, 1: 1/4, 2: 3/8, 3: 1/4, 4: 1/16}

    """
    return rv(name, BinomialDistribution, n, p, succ, fail)


class HypergeometricDistribution(SingleFiniteDistribution):
    _argnames = ('N', 'm', 'n')

    @property  # type: ignore[misc]
    @cacheit
    def dict(self):
        N, m, n = self.N, self.m, self.n
        N, m, n = list(map(sympify, (N, m, n)))
        density = {sympify(k):
                   Rational(binomial(m, k) * binomial(N - m, n - k),
                            binomial(N, n))
                   for k in range(max(0, n + m - N), min(m, n) + 1)}
        return density

    @property
    def n(self):
        return self.args[self._argnames.index('n')]


def Hypergeometric(name, N, m, n):
    """
    Create a Finite Random Variable representing a hypergeometric distribution.

    Returns a RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> X = Hypergeometric('X', 10, 5, 3)  # 10 marbles, 5 white (success), 3 draws
    >>> density(X).dict
    {0: 1/12, 1: 5/12, 2: 5/12, 3: 1/12}

    """
    return rv(name, HypergeometricDistribution, N, m, n)


class RademacherDistribution(SingleFiniteDistribution):
    @property  # type: ignore[misc]
    @cacheit
    def dict(self):
        return {-1: Rational(1, 2), 1: Rational(1, 2)}


def Rademacher(name):
    """
    Create a Finite Random Variable representing a Rademacher distribution.

    Return a RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> X = Rademacher('X')
    >>> density(X).dict
    {-1: 1/2, 1: 1/2}

    See Also
    ========

    diofant.stats.Bernoulli

    References
    ==========

    * https://en.wikipedia.org/wiki/Rademacher_distribution

    """
    return rv(name, RademacherDistribution)
