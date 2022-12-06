"""
Finite Discrete Random Variables Module

See Also
========
diofant.stats.frv_types
diofant.stats.rv
diofant.stats.crv

"""

import random
from itertools import product

from ..core import Dict, Eq, Expr, Lambda, Mul, Symbol, Tuple, cacheit, sympify
from ..functions import Piecewise
from ..logic import And, Or
from ..sets import FiniteSet
from .rv import (ConditionalDomain, NamedArgsMixin, ProductDomain,
                 ProductPSpace, PSpace, RandomDomain, SinglePSpace,
                 random_symbols, rv_subs, sumsets)


class FiniteDensity(dict):
    """Finite probabibility density."""

    def __call__(self, item):
        item = sympify(item)
        if item in self:
            return self[item]
        else:
            return 0

    @property
    def dict(self):
        return dict(self)


class FiniteDomain(RandomDomain):
    """
    A domain with discrete finite support

    Represented using a FiniteSet.

    """

    is_Finite = True

    @property
    def dict(self):
        return FiniteSet(*[Dict(dict(el)) for el in self.elements])

    def as_boolean(self):
        return Or(*[And(*[Eq(sym, val) for sym, val in item]) for item in self])

    def __iter__(self):
        raise NotImplementedError


class SingleFiniteDomain(FiniteDomain):
    """
    A FiniteDomain over a single symbol/set

    Example: The possibilities of a *single* die roll.

    """

    def __new__(cls, symbol, set):
        if not isinstance(set, FiniteSet):
            set = FiniteSet(*set)
        return Expr.__new__(cls, symbol, set)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def symbols(self):
        return FiniteSet(self.symbol)

    @property
    def set(self):
        return self.args[1]

    @property
    def elements(self):
        return FiniteSet(*[frozenset(((self.symbol, elem), )) for elem in self.set])

    def __iter__(self):
        return (frozenset(((self.symbol, elem),)) for elem in self.set)

    def __contains__(self, other):
        sym, val = tuple(other)[0]
        return sym == self.symbol and val in self.set


class ProductFiniteDomain(ProductDomain, FiniteDomain):
    """
    A Finite domain consisting of several other FiniteDomains

    Example: The possibilities of the rolls of three independent dice

    """

    def __iter__(self):
        proditer = product(*self.domains)
        return (sumsets(items) for items in proditer)

    @property
    def elements(self):
        return FiniteSet(*self)


class ConditionalFiniteDomain(ConditionalDomain, ProductFiniteDomain):
    """
    A FiniteDomain that has been restricted by a condition

    Example: The possibilities of a die roll under the condition that the
    roll is even.

    """

    def __new__(cls, domain, condition):
        if condition is True:
            return domain
        cond = rv_subs(condition)
        # Check that we aren't passed a condition like die1 == z
        # where 'z' is a symbol that we don't know about
        # We will never be able to test this equality through iteration
        if not cond.free_symbols.issubset(domain.free_symbols):
            raise ValueError(f'Condition "{condition}" contains foreign symbols'
                             f' \n{tuple(cond.free_symbols - domain.free_symbols)}.\n'
                             'Will be unable to iterate using this condition')

        return Expr.__new__(cls, domain, cond)

    def _test(self, elem):
        val = self.condition.xreplace(dict(elem))
        if val in [True, False]:
            return val
        elif val.is_Equality:
            return val.lhs == val.rhs
        raise ValueError(f'Undeciable if {val!s}')

    def __contains__(self, other):
        return other in self.fulldomain and self._test(other)

    def __iter__(self):
        return (elem for elem in self.fulldomain if self._test(elem))

    @property
    def set(self):
        if self.fulldomain.__class__ is SingleFiniteDomain:
            return FiniteSet(*[elem for elem in self.fulldomain.set
                               if frozenset(((self.fulldomain.symbol, elem),)) in self])
        else:
            raise NotImplementedError(
                'Not implemented on multi-dimensional conditional domain')

    def as_boolean(self):
        return FiniteDomain.as_boolean(self)


class SingleFiniteDistribution(Expr, NamedArgsMixin):
    """Base class for finite distributions."""

    def __new__(cls, *args):
        args = list(map(sympify, args))
        return Expr.__new__(cls, *args)

    @property  # type: ignore[misc]
    @cacheit
    def dict(self):
        return {k: self.pdf(k) for k in self.set}

    @property
    def pdf(self):
        x = Symbol('x')
        return Lambda(x, Piecewise(*(
            [(v, Eq(k, x)) for k, v in self.dict.items()] + [(0, True)])))

    @property
    def set(self):
        return list(self.dict)

    values = property(lambda self: self.dict.values)
    items = property(lambda self: self.dict.items)
    __iter__ = property(lambda self: self.dict.__iter__)
    __getitem__ = property(lambda self: self.dict.__getitem__)

    __call__ = pdf

    def __contains__(self, other):
        return other in self.set


#####################
# Probability Space #
#####################


class FinitePSpace(PSpace):
    """
    A Finite Probability Space

    Represents the probabilities of a finite number of events.

    """

    is_Finite = True

    @property
    def domain(self):
        return self.args[0]

    @property
    def density(self):
        return self.args[0]

    def __new__(cls, domain, density):
        density = {sympify(key): sympify(val)
                   for key, val in density.items()}
        public_density = Dict(density)

        obj = PSpace.__new__(cls, domain, public_density)
        obj._density = density
        return obj

    def prob_of(self, elem):
        return self._density.get(elem, 0)

    def where(self, condition):
        return ConditionalFiniteDomain(self.domain, condition)

    def compute_density(self, expr):
        expr = expr.xreplace({rs: rs.symbol for rs in self.values})
        d = FiniteDensity()
        for elem in self.domain:
            val = expr.xreplace(dict(elem))
            prob = self.prob_of(elem)
            d[val] = d.get(val, 0) + prob
        return d

    @cacheit
    def compute_cdf(self, expr):
        d = self.compute_density(expr)
        cum_prob = 0
        cdf = []
        for key in sorted(d):
            prob = d[key]
            cum_prob += prob
            cdf.append((key, cum_prob))

        return dict(cdf)

    @cacheit
    def sorted_cdf(self, expr, python_float=False):
        cdf = self.compute_cdf(expr)
        items = list(cdf.items())
        sorted_items = sorted(items, key=lambda val_cumprob: val_cumprob[1])
        if python_float:
            sorted_items = [(v, float(cum_prob))
                            for v, cum_prob in sorted_items]
        return sorted_items

    def integrate(self, expr, rvs=None):
        rvs = rvs or self.values
        expr = expr.xreplace({rs: rs.symbol for rs in rvs})
        return sum(expr.xreplace(dict(elem)) * self.prob_of(elem)
                   for elem in self.domain)

    def probability(self, condition):
        cond_symbols = frozenset(rs.symbol for rs in random_symbols(condition))
        assert cond_symbols.issubset(self.symbols)
        return sum(self.prob_of(elem) for elem in self.where(condition))

    def conditional_space(self, condition):
        domain = self.where(condition)
        prob = self.probability(condition)
        density = {key: val / prob
                   for key, val in self._density.items() if key in domain}
        return FinitePSpace(domain, density)

    def sample(self):
        """
        Internal sample method

        Returns dictionary mapping RandomSymbol to realization value.

        """
        expr = Tuple(*self.values)
        cdf = self.sorted_cdf(expr, python_float=True)

        x = random.uniform(0, 1)
        # Find first occurrence with cumulative probability less than x
        # This should be replaced with binary search
        for value, cum_prob in cdf:  # pragma: no branch
            if x < cum_prob:
                # return dictionary mapping RandomSymbols to values
                return dict(zip(expr, value))


class SingleFinitePSpace(SinglePSpace, FinitePSpace):
    """
    A single finite probability space

    Represents the probabilities of a set of random events that can be
    attributed to a single variable/symbol.

    This class is implemented by many of the standard FiniteRV types such as
    Die, Bernoulli, Coin, etc....

    """

    @property
    def domain(self):
        return SingleFiniteDomain(self.symbol, self.distribution.set)

    @property  # type: ignore[misc]
    @cacheit
    def _density(self):
        return {frozenset(((self.symbol, val),)): prob
                for val, prob in self.distribution.dict.items()}


class ProductFinitePSpace(ProductPSpace, FinitePSpace):
    """A collection of several independent finite probability spaces."""

    @property
    def domain(self):
        return ProductFiniteDomain(*[space.domain for space in self.spaces])

    @property  # type: ignore[misc]
    @cacheit
    def _density(self):
        proditer = product(*[iter(space._density.items())
                             for space in self.spaces])
        d = {}
        for items in proditer:
            elems, probs = list(zip(*items))
            elem = sumsets(elems)
            prob = Mul(*probs)
            d[elem] = d.get(elem, 0) + prob
        return Dict(d)

    @property  # type: ignore[misc]
    @cacheit
    def density(self):
        return self._density
