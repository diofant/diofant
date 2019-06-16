"""
Diofant statistics module

Introduces a random variable type into the Diofant language.

Random variables may be declared using prebuilt functions such as
Normal, Exponential, Coin, Die, etc...  or built with functions like FiniteRV.

Queries on random expressions can be made using the functions

========================= =============================
    Expression                    Meaning
------------------------- -----------------------------
 ``P(condition)``          Probability
 ``E(expression)``         Expected value
 ``variance(expression)``  Variance
 ``density(expression)``   Probability Density Function
 ``sample(expression)``    Produce a realization
 ``where(condition)``      Where the condition is true
========================= =============================

Examples
========

>>> from diofant.stats import P, E, variance, Die, Normal
>>> X, Y = Die('X', 6), Die('Y', 6) # Define two six sided dice
>>> Z = Normal('Z', 0, 1) # Declare a Normal random variable with mean 0, std 1
>>> P(X>3) # Probability X is greater than 3
1/2
>>> E(X+Y) # Expectation of the sum of two dice
7
>>> variance(X+Y) # Variance of the sum of two dice
35/6
>>> simplify(P(Z>1)) # Probability of Z being greater than 1
-erf(sqrt(2)/2)/2 + 1/2

"""

from .crv_types import (Arcsin, Benini, Beta, BetaPrime, Cauchy, Chi,
                        ChiNoncentral, ChiSquared, ContinuousRV, Dagum, Erlang,
                        Exponential, FDistribution, FisherZ, Frechet, Gamma,
                        GammaInverse, Kumaraswamy, Laplace, Logistic,
                        LogNormal, Maxwell, Nakagami, Normal, Pareto,
                        QuadraticU, RaisedCosine, Rayleigh, StudentT,
                        Triangular, Uniform, UniformSum, VonMises, Weibull,
                        WignerSemicircle)
from .drv_types import Geometric, Poisson
from .frv_types import (Bernoulli, Binomial, Coin, Die, DiscreteUniform,
                        FiniteRV, Hypergeometric, Rademacher)
from .rv_interface import (E, P, cdf, cmoment, correlation, covariance,
                           density, dependent, given, independent, moment,
                           pspace, random_symbols, sample, sample_iter,
                           sampling_density, skewness, smoment, std, variance,
                           where)
