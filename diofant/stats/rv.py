"""
Main Random Variables Module.

Defines abstract random variable type.
Contains interfaces for probability space object (PSpace) as well as standard
operators, P, E, sample, density, where

See Also
========
diofant.stats.crv
diofant.stats.frv
diofant.stats.rv_interface

"""

from __future__ import annotations

from ..abc import x
from ..core import (Add, Eq, Equality, Expr, Integer, Lambda, Symbol, Tuple,
                    oo, sympify)
from ..core.logic import fuzzy_or
from ..core.relational import Relational
from ..functions import DiracDelta
from ..logic.boolalg import Boolean, false, true
from ..sets import FiniteSet, ProductSet
from ..solvers import solve
from ..utilities import lambdify


class RandomDomain(Expr):
    """
    Represents a set of variables and the values which they can take

    See Also
    ========
    diofant.stats.crv.ContinuousDomain
    diofant.stats.frv.FiniteDomain

    """

    is_ProductDomain = False
    is_Finite = False
    is_Continuous = False

    def __new__(cls, symbols, *args):
        symbols = FiniteSet(*symbols)
        return Expr.__new__(cls, symbols, *args)

    @property
    def set(self):
        return self.args[1]

    def __contains__(self, other):
        raise NotImplementedError

    def integrate(self, expr):
        raise NotImplementedError


class SingleDomain(RandomDomain):
    """
    A single variable and its domain

    See Also
    ========

    diofant.stats.crv.SingleContinuousDomain
    diofant.stats.frv.SingleFiniteDomain

    """

    def __new__(cls, symbol, set):
        assert symbol.is_Symbol
        return Expr.__new__(cls, symbol, set)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def symbols(self):
        return FiniteSet(self.symbol)


class ConditionalDomain(RandomDomain):
    """
    A RandomDomain with an attached condition

    See Also
    ========

    diofant.stats.crv.ConditionalContinuousDomain
    diofant.stats.frv.ConditionalFiniteDomain

    """

    def __new__(cls, fulldomain, condition):
        condition = condition.xreplace({rs: rs.symbol
                                        for rs in random_symbols(condition)})
        return Expr.__new__(cls, fulldomain, condition)

    @property
    def symbols(self):
        return self.fulldomain.symbols

    @property
    def fulldomain(self):
        return self.args[0]

    @property
    def condition(self):
        return self.args[1]

    @property
    def set(self):
        raise NotImplementedError('Set of Conditional Domain not Implemented')


class PSpace(Expr):
    """
    A Probability Space

    Probability Spaces encode processes that equal different values
    probabilistically. These underly Random Symbols which occur in Diofant
    expressions and contain the mechanics to evaluate statistical statements.

    See Also
    ========
    diofant.stats.crv.ContinuousPSpace
    diofant.stats.frv.FinitePSpace

    """

    is_Finite: bool | None = None
    is_Continuous: bool | None = None

    @property
    def values(self):
        return frozenset(RandomSymbol(self, sym) for sym in self.domain.symbols)

    @property
    def symbols(self):
        return self.domain.symbols

    def where(self, condition):
        raise NotImplementedError

    def compute_density(self, expr):
        raise NotImplementedError

    def sample(self):
        raise NotImplementedError

    def probability(self, condition):
        raise NotImplementedError

    def integrate(self, expr):
        raise NotImplementedError


class SinglePSpace(PSpace):
    """
    Represents the probabilities of a set of random events that can be
    attributed to a single variable/symbol.

    """

    def __new__(cls, s, distribution):
        if isinstance(s, str):
            s = Symbol(s)
        if not isinstance(s, Symbol):
            raise TypeError('s should have been string or Symbol')
        return Expr.__new__(cls, s, distribution)

    @property
    def value(self):
        return RandomSymbol(self, self.symbol)

    @property
    def symbol(self):
        return self.args[0]

    @property
    def distribution(self):
        return self.args[1]


class RandomSymbol(Expr):
    """
    Random Symbols represent ProbabilitySpaces in Diofant Expressions
    In principle they can take on any value that their symbol can take on
    within the associated PSpace with probability determined by the PSpace
    Density.

    Random Symbols contain pspace and symbol properties.
    The pspace property points to the represented Probability Space
    The symbol is a standard Diofant Symbol that is used in that probability space
    for example in defining a density.

    You can form normal Diofant expressions using RandomSymbols and operate on
    those expressions with the Functions

    E - Expectation of a random expression
    P - Probability of a condition
    density - Probability Density of an expression
    given - A new random expression (with new random symbols) given a condition

    An object of the RandomSymbol type should almost never be created by the
    user. They tend to be created instead by the PSpace class's value method.
    Traditionally a user doesn't even do this but instead calls one of the
    convenience functions Normal, Exponential, Coin, Die, FiniteRV, etc....

    """

    def __new__(cls, pspace, symbol):
        if not isinstance(symbol, Symbol):
            raise TypeError('symbol should be of type Symbol')
        if not isinstance(pspace, PSpace):
            raise TypeError('pspace variable should be of type PSpace')
        return Expr.__new__(cls, pspace, symbol)

    is_finite = True
    is_Symbol = True
    is_Atom = True

    _diff_wrt = True

    pspace = property(lambda self: self.args[0])
    symbol = property(lambda self: self.args[1])
    name = property(lambda self: self.symbol.name)

    def _eval_is_positive(self):
        return self.symbol.is_positive

    def _eval_is_integer(self):
        return self.symbol.is_integer

    def _eval_is_extended_real(self):
        return fuzzy_or([self.symbol.is_extended_real,
                         self.pspace.is_extended_real])

    def _eval_is_commutative(self):
        return self.symbol.is_commutative

    def _hashable_content(self):
        return self.pspace, self.symbol

    @property
    def free_symbols(self):
        return {self}


class ProductPSpace(PSpace):
    """
    A probability space resulting from the merger of two independent probability
    spaces.

    Often created using the function, pspace

    """

    def __new__(cls, *spaces):
        rs_space_dict = {}
        for space in spaces:
            for value in space.values:
                rs_space_dict[value] = space

        symbols = FiniteSet(*[val.symbol for val in rs_space_dict])

        # Overlapping symbols
        if len(symbols) < sum(len(space.symbols) for space in spaces):
            raise ValueError('Overlapping Random Variables')

        new_cls = cls
        if all(space.is_Finite for space in spaces):
            from .frv import ProductFinitePSpace
            new_cls = ProductFinitePSpace
        if all(space.is_Continuous for space in spaces):
            from .crv import ProductContinuousPSpace
            new_cls = ProductContinuousPSpace

        obj = Expr.__new__(new_cls, *FiniteSet(*spaces))

        return obj

    @property
    def rs_space_dict(self):
        d = {}
        for space in self.spaces:
            for value in space.values:
                d[value] = space
        return d

    @property
    def symbols(self):
        return FiniteSet(*[val.symbol for val in self.rs_space_dict])

    @property
    def spaces(self):
        return FiniteSet(*self.args)

    @property
    def values(self):
        return sumsets(space.values for space in self.spaces)

    def integrate(self, expr, rvs=None, **kwargs):
        rvs = rvs or self.values
        rvs = frozenset(rvs)
        for space in self.spaces:
            expr = space.integrate(expr, rvs & space.values, **kwargs)
        return expr

    @property
    def domain(self):
        return ProductDomain(*[space.domain for space in self.spaces])

    @property
    def density(self):
        raise NotImplementedError('Density not available for ProductSpaces')

    def sample(self):
        return {k: v for space in self.spaces
                for k, v in space.sample().items()}


class ProductDomain(RandomDomain):
    """
    A domain resulting from the merger of two independent domains

    See Also
    ========

    diofant.stats.crv.ProductContinuousDomain
    diofant.stats.frv.ProductFiniteDomain

    """

    is_ProductDomain = True

    def __new__(cls, *domains):
        # Flatten any product of products
        domains2 = []
        for domain in domains:
            if not domain.is_ProductDomain:
                domains2.append(domain)
            else:
                domains2.extend(domain.domains)
        domains2 = FiniteSet(*domains2)

        new_cls = cls
        if all(domain.is_Finite for domain in domains2):
            from .frv import ProductFiniteDomain
            new_cls = ProductFiniteDomain
        if all(domain.is_Continuous for domain in domains2):
            from .crv import ProductContinuousDomain
            new_cls = ProductContinuousDomain

        return Expr.__new__(new_cls, *domains2)

    @property
    def symbols(self):
        return FiniteSet(*[sym for domain in self.domains
                           for sym in domain.symbols])

    @property
    def domains(self):
        return self.args

    @property
    def set(self):
        return ProductSet(domain.set for domain in self.domains)

    def __contains__(self, other):
        # Split event into each subdomain
        for domain in self.domains:
            # Collect the parts of this event which associate to this domain
            elem = frozenset(item for item in other
                             if domain.symbols.contains(item[0]) == true)
            # Test this sub-event
            if elem not in domain:
                return False
        # All subevents passed
        return True


def random_symbols(expr):
    """Returns all RandomSymbols within a Diofant Expression."""
    try:
        return list(expr.atoms(RandomSymbol))
    except AttributeError:
        return []


def pspace(expr):
    """
    Returns the underlying Probability Space of a random expression.

    For internal use.

    Examples
    ========

    >>> from diofant.stats import Normal
    >>> X = Normal('X', 0, 1)
    >>> pspace(2*X + 1) == X.pspace
    True

    """
    expr = sympify(expr)
    rvs = random_symbols(expr)
    if not rvs:
        raise ValueError(f'Expression containing Random Variable expected, not {expr}')
    # If only one space present
    if all(rv.pspace == rvs[0].pspace for rv in rvs):
        return rvs[0].pspace
    # Otherwise make a product space
    return ProductPSpace(*[rv.pspace for rv in rvs])


def sumsets(sets):
    """Union of sets."""
    return frozenset().union(*sets)


def rs_swap(a, b):
    """
    Build a dictionary to swap RandomSymbols based on their underlying symbol.

    i.e.
    if    ``X = ('x', pspace1)``
    and   ``Y = ('x', pspace2)``
    then ``X`` and ``Y`` match and the key, value pair
    ``{X:Y}`` will appear in the result

    Inputs: collections a and b of random variables which share common symbols
    Output: dict mapping RVs in a to RVs in b

    """
    d = {}
    for rsa in a:
        d[rsa] = [rsb for rsb in b if rsa.symbol == rsb.symbol][0]
    return d


def given(expr, condition=None, **kwargs):
    r"""Conditional Random Expression.

    From a random expression and a condition on that expression creates a new
    probability space from the condition and returns the same expression on that
    conditional probability space.

    Examples
    ========

    >>> from diofant.stats import Die
    >>> X = Die('X', 6)
    >>> Y = given(X, X > 3)
    >>> density(Y).dict
    {4: 1/3, 5: 1/3, 6: 1/3}

    Following convention, if the condition is a random symbol then that symbol
    is considered fixed.

    >>> from diofant.stats import Normal

    >>> X = Normal('X', 0, 1)
    >>> Y = Normal('Y', 0, 1)
    >>> pprint(density(X + Y, Y)(z), use_unicode=False)
                    2
           -(-Y + z)
           -----------
      ___       2
    \/ 2 *E
    ------------------
             ____
         2*\/ pi

    """
    if not random_symbols(condition) or pspace_independent(expr, condition):
        return expr

    if isinstance(condition, RandomSymbol):
        condition = Eq(condition, condition.symbol)

    condsymbols = random_symbols(condition)
    if (isinstance(condition, Equality) and len(condsymbols) == 1 and
            not isinstance(pspace(expr).domain, ConditionalDomain)):
        rv = tuple(condsymbols)[0]
        results = solve(condition, rv)
        return sum(expr.subs(res) for res in results)

    # Get full probability space of both the expression and the condition
    fullspace = pspace(Tuple(expr, condition))
    # Build new space given the condition
    space = fullspace.conditional_space(condition, **kwargs)
    # Dictionary to swap out RandomSymbols in expr with new RandomSymbols
    # That point to the new conditional space
    swapdict = rs_swap(fullspace.values, space.values)
    # Swap random variables in the expression
    expr = expr.xreplace(swapdict)
    return expr


def expectation(expr, condition=None, numsamples=None, evaluate=True, **kwargs):
    """
    Returns the expected value of a random expression

    Parameters
    ==========

    expr : Expr containing RandomSymbols
        The expression of which you want to compute the expectation value
    given : Expr containing RandomSymbols
        A conditional expression. E(X, X>0) is expectation of X given X > 0
    numsamples : int
        Enables sampling and approximates the expectation with this many samples
    evalf : Bool (defaults to True)
        If sampling return a number rather than a complex expression
    evaluate : Bool (defaults to True)
        In case of continuous systems return unevaluated integral

    Examples
    ========

    >>> from diofant.stats import Die, E
    >>> X = Die('X', 6)
    >>> E(X)
    7/2
    >>> E(2*X + 1)
    8

    >>> E(X, X > 3)  # Expectation of X given that it is above 3
    5

    """
    if not random_symbols(expr):  # expr isn't random?
        return expr
    if numsamples:  # Computing by monte carlo sampling?
        return sampling_E(expr, condition, numsamples=numsamples, **kwargs)

    # Create new expr and recompute E
    if condition is not None:  # If there is a condition
        return expectation(given(expr, condition), evaluate=evaluate)

    # A few known statements for efficiency

    if expr.is_Add:  # We know that E is Linear
        return Add(*[expectation(arg, evaluate=evaluate)
                     for arg in expr.args])

    # Otherwise case is simple, pass work off to the ProbabilitySpace
    result = pspace(expr).integrate(expr)
    if evaluate and hasattr(result, 'doit'):
        return result.doit(**kwargs)
    else:
        return result


def probability(condition, given_condition=None, numsamples=None,
                evaluate=True, **kwargs):
    """
    Probability that a condition is true, optionally given a second condition

    Parameters
    ==========

    condition : Combination of Relationals containing RandomSymbols
        The condition of which you want to compute the probability
    given_condition : Combination of Relationals containing RandomSymbols
        A conditional expression. P(X > 1, X > 0) is expectation of X > 1
        given X > 0
    numsamples : int
        Enables sampling and approximates the probability with this many samples
    evaluate : Bool (defaults to True)
        In case of continuous systems return unevaluated integral

    Examples
    ========

    >>> from diofant.stats import Die, P
    >>> X, Y = Die('X', 6), Die('Y', 6)
    >>> P(X > 3)
    1/2
    >>> P(Eq(X, 5), X > 2)  # Probability that X == 5 given that X > 2
    1/4
    >>> P(X > Y)
    5/12

    """
    condition = sympify(condition)
    given_condition = sympify(given_condition)

    if given_condition is not None and \
            not isinstance(given_condition, (Relational, Boolean)):
        raise ValueError(f'{given_condition} is not a relational or combination of relationals')
    if given_condition == false:
        return Integer(0)
    if not isinstance(condition, (Relational, Boolean)):
        raise ValueError(f'{condition} is not a relational or combination of relationals')
    if condition == true:
        return Integer(1)
    if condition == false:
        return Integer(0)

    if numsamples:
        return sampling_P(condition, given_condition, numsamples=numsamples,
                          **kwargs)
    if given_condition is not None:  # If there is a condition
        # Recompute on new conditional expr
        return probability(given(condition, given_condition, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    result = pspace(condition).probability(condition, **kwargs)
    if evaluate and hasattr(result, 'doit'):
        return result.doit()
    else:
        return result


class Density(Expr):
    """Probability density."""

    expr = property(lambda self: self.args[0])

    @property
    def condition(self):
        if len(self.args) > 1:
            return self.args[1]

    def doit(self, **kwargs):
        evaluate = kwargs.pop('evaluate', True)

        expr, condition = self.expr, self.condition
        if condition is not None:
            # Recompute on new conditional expr
            expr = given(expr, condition, **kwargs)
        if not random_symbols(expr):
            return Lambda(x, DiracDelta(x - expr))
        if (isinstance(expr, RandomSymbol) and
                hasattr(expr.pspace, 'distribution') and
                isinstance(pspace(expr), SinglePSpace)):
            return expr.pspace.distribution
        result = pspace(expr).compute_density(expr, **kwargs)

        if evaluate and hasattr(result, 'doit'):
            return result.doit()
        else:
            return result


def density(expr, condition=None, evaluate=True, numsamples=None, **kwargs):
    """
    Probability density of a random expression, optionally given a second
    condition.

    This density will take on different forms for different types of
    probability spaces. Discrete variables produce Dicts. Continuous
    variables produce Lambdas.

    Parameters
    ==========

    expr : Expr containing RandomSymbols
        The expression of which you want to compute the density value
    condition : Relational containing RandomSymbols
        A conditional expression. density(X > 1, X > 0) is density of X > 1
        given X > 0
    numsamples : int
        Enables sampling and approximates the density with this many samples

    Examples
    ========

    >>> from diofant.stats import Die, Normal

    >>> D = Die('D', 6)
    >>> X = Normal(x, 0, 1)

    >>> density(D).dict
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}
    >>> density(2*D).dict
    {2: 1/6, 4: 1/6, 6: 1/6, 8: 1/6, 10: 1/6, 12: 1/6}
    >>> density(X)(x)
    sqrt(2)*E**(-x**2/2)/(2*sqrt(pi))

    """
    if numsamples:
        return sampling_density(expr, condition, numsamples=numsamples,
                                **kwargs)

    kwargs['evaluate'] = evaluate

    return Density(expr, condition).doit(**kwargs)


def cdf(expr, condition=None, evaluate=True, **kwargs):
    """
    Cumulative Distribution Function of a random expression.

    optionally given a second condition

    This density will take on different forms for different types of
    probability spaces.
    Discrete variables produce Dicts.
    Continuous variables produce Lambdas.

    Examples
    ========

    >>> from diofant.stats import Die, Normal

    >>> D = Die('D', 6)
    >>> X = Normal('X', 0, 1)

    >>> density(D).dict
    {1: 1/6, 2: 1/6, 3: 1/6, 4: 1/6, 5: 1/6, 6: 1/6}
    >>> cdf(D)
    {1: 1/6, 2: 1/3, 3: 1/2, 4: 2/3, 5: 5/6, 6: 1}
    >>> cdf(3*D, D > 2)
    {9: 1/4, 12: 1/2, 15: 3/4, 18: 1}

    >>> cdf(X)
    Lambda(_z, erf(sqrt(2)*_z/2)/2 + 1/2)

    """
    if condition is not None:  # If there is a condition
        # Recompute on new conditional expr
        return cdf(given(expr, condition, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    result = pspace(expr).compute_cdf(expr, **kwargs)

    if evaluate and hasattr(result, 'doit'):
        return result.doit()
    else:
        return result


def where(condition, given_condition=None, **kwargs):
    """
    Returns the domain where a condition is True.

    Examples
    ========

    >>> from diofant.stats import Die, Normal

    >>> D1, D2 = Die('a', 6), Die('b', 6)
    >>> a, b = D1.symbol, D2.symbol
    >>> X = Normal('x', 0, 1)

    >>> where(X**2 < 1)
    Domain: (-1 < x) & (x < 1)

    >>> where(X**2 < 1).set
    (-1, 1)

    >>> where(And(D1 <= D2, D2 < 3))
    Domain: (Eq(a, 1) & Eq(b, 1)) | (Eq(a, 1) & Eq(b, 2)) | (Eq(a, 2) & Eq(b, 2))

    """
    if given_condition is not None:  # If there is a condition
        # Recompute on new conditional expr
        return where(given(condition, given_condition, **kwargs), **kwargs)

    # Otherwise pass work off to the ProbabilitySpace
    return pspace(condition).where(condition, **kwargs)


def sample(expr, condition=None, **kwargs):
    """
    A realization of the random expression

    Examples
    ========

    >>> from diofant.stats import Die
    >>> X, Y, Z = Die('X', 6), Die('Y', 6), Die('Z', 6)

    >>> die_roll = sample(X + Y + Z)  # A random realization of three dice

    """
    return next(sample_iter(expr, condition, numsamples=1))


def sample_iter(expr, condition=None, numsamples=oo, **kwargs):
    """
    Returns an iterator of realizations from the expression given a condition

    expr: Random expression to be realized
    condition: A conditional expression (optional)
    numsamples: Length of the iterator (defaults to infinity)

    Examples
    ========

    >>> from diofant.stats import Normal
    >>> X = Normal('X', 0, 1)
    >>> expr = X*X + 3
    >>> iterator = sample_iter(expr, numsamples=3)
    >>> list(iterator)  # doctest: +SKIP
    [12, 4, 7]

    See Also
    ========

    diofant.stats.sample
    diofant.stats.rv.sampling_P
    diofant.stats.rv.sampling_E

    """
    if condition is not None:
        ps = pspace(Tuple(expr, condition))
    else:
        ps = pspace(expr)

    rvs = list(ps.values)
    fn = lambdify(rvs, expr, **kwargs)
    if condition is not None:
        given_fn = lambdify(rvs, condition, **kwargs)

    # Check that lambdify can handle the expression
    # Some operations like Sum can prove difficult
    try:
        d = ps.sample()  # a dictionary that maps RVs to values
        args = [d[rv] for rv in rvs]
        fn(*args)
        if condition is not None:
            given_fn(*args)
    except (TypeError, ValueError) as exc:
        raise TypeError('Expr/condition too complex for lambdify') from exc

    def return_generator():
        count = 0
        while count < numsamples:
            d = ps.sample()  # a dictionary that maps RVs to values
            args = [d[rv] for rv in rvs]

            if condition is not None:  # Check that these values satisfy the condition
                gd = given_fn(*args)
                if gd not in (True, False):
                    raise ValueError(
                        'Conditions must not contain free symbols')
                if not gd:  # If the values don't satisfy then try again
                    continue

            yield fn(*args)
            count += 1
    return return_generator()


def sampling_P(condition, given_condition=None, numsamples=1,
               evalf=True, **kwargs):
    """
    Sampling version of P

    See Also
    ========
    diofant.stats.P
    diofant.stats.rv.sampling_E
    diofant.stats.rv.sampling_density

    """
    count_true = 0
    count_false = 0

    samples = sample_iter(condition, given_condition,
                          numsamples=numsamples, **kwargs)

    for x in samples:
        if x:
            count_true += 1
        else:
            count_false += 1

    result = Integer(count_true) / numsamples
    return result.evalf()


def sampling_E(expr, given_condition=None, numsamples=1,
               evalf=True, **kwargs):
    """
    Sampling version of E

    See Also
    ========
    diofant.stats.P
    diofant.stats.rv.sampling_P
    diofant.stats.rv.sampling_density

    """
    samples = sample_iter(expr, given_condition,
                          numsamples=numsamples, **kwargs)

    result = Add(*list(samples)) / numsamples
    return result.evalf(strict=False)


def sampling_density(expr, given_condition=None, numsamples=1, **kwargs):
    """
    Sampling version of density

    See Also
    ========
    diofant.stats.density
    diofant.stats.rv.sampling_P
    diofant.stats.rv.sampling_E

    """
    results = {}
    for result in sample_iter(expr, given_condition,
                              numsamples=numsamples, **kwargs):
        results[result] = results.get(result, 0) + 1
    return results


def dependent(a, b):
    """
    Dependence of two random expressions

    Two expressions are independent if knowledge of one does not change
    computations on the other.

    Examples
    ========

    >>> from diofant.stats import Normal

    >>> X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    >>> dependent(X, Y)
    False
    >>> dependent(2*X + Y, -Y)
    True
    >>> X, Y = given(Tuple(X, Y), Eq(X + Y, 3))
    >>> dependent(X, Y)
    True

    See Also
    ========
    diofant.stats.rv.independent

    """
    if pspace_independent(a, b):
        return False

    z = Symbol('z', extended_real=True)
    # Dependent if density is unchanged when one is given information about
    # the other
    return (density(a, Eq(b, z)) != density(a) or
            density(b, Eq(a, z)) != density(b))


def independent(a, b):
    """
    Independence of two random expressions

    Two expressions are independent if knowledge of one does not change
    computations on the other.

    Examples
    ========

    >>> from diofant.stats import Normal

    >>> X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    >>> independent(X, Y)
    True
    >>> independent(2*X + Y, -Y)
    False
    >>> X, Y = given(Tuple(X, Y), Eq(X + Y, 3))
    >>> independent(X, Y)
    False

    See Also
    ========
    diofant.stats.rv.dependent

    """
    return not dependent(a, b)


def pspace_independent(a, b):
    """
    Tests for independence between a and b by checking if their PSpaces have
    overlapping symbols. This is a sufficient but not necessary condition for
    independence and is intended to be used internally.

    Notes
    =====

    pspace_independent(a, b) implies independent(a, b)
    independent(a, b) does not imply pspace_independent(a, b)

    """
    a_symbols = set(pspace(b).symbols)
    b_symbols = set(pspace(a).symbols)

    if len(a_symbols.intersection(b_symbols)) == 0:
        return True


def rv_subs(expr):
    """Given a random expression replace all random variables with their symbols."""
    symbols = random_symbols(expr)
    if not symbols:
        return expr
    swapdict = {rv: rv.symbol for rv in symbols}
    return expr.xreplace(swapdict)


class NamedArgsMixin:
    """Helper class for named arguments."""

    _argnames: tuple[str, ...] = ()

    def __getattr__(self, attr):
        try:
            return self.args[self._argnames.index(attr)]
        except ValueError as exc:
            raise AttributeError(f"'{type(self).__name__}' object "
                                 f"has not attribute '{attr}'") from exc


def _value_check(condition, message):
    """
    Check a condition on input value.

    Raises ValueError with message if condition is not True

    """
    if condition != true:
        raise ValueError(message)
