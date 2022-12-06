"""
Continuous Random Variables - Prebuilt variables

Contains
========
Arcsin
Benini
Beta
BetaPrime
Cauchy
Chi
ChiNoncentral
ChiSquared
Dagum
Erlang
Exponential
FDistribution
FisherZ
Frechet
Gamma
GammaInverse
Kumaraswamy
Laplace
Logistic
LogNormal
Maxwell
Nakagami
Normal
Pareto
QuadraticU
RaisedCosine
Rayleigh
StudentT
Triangular
Uniform
UniformSum
VonMises
Weibull
WignerSemicircle
"""

import random

from ..concrete import Sum
from ..core import Dummy, Eq, Expr, Lambda, Rational, oo, pi, sympify
from ..functions import Abs, Piecewise, besseli
from ..functions import beta as beta_fn
from ..functions import binomial, cos, exp, factorial, floor, gamma, log, sqrt
from ..logic import And
from ..sets import Interval
from .crv import (ContinuousDistributionHandmade, SingleContinuousDistribution,
                  SingleContinuousPSpace)
from .rv import _value_check


__all__ = ('ContinuousRV',
           'Arcsin',
           'Benini',
           'Beta',
           'BetaPrime',
           'Cauchy',
           'Chi',
           'ChiNoncentral',
           'ChiSquared',
           'Dagum',
           'Erlang',
           'Exponential',
           'FDistribution',
           'FisherZ',
           'Frechet',
           'Gamma',
           'GammaInverse',
           'Kumaraswamy',
           'Laplace',
           'Logistic',
           'LogNormal',
           'Maxwell',
           'Nakagami',
           'Normal',
           'Pareto',
           'QuadraticU',
           'RaisedCosine',
           'Rayleigh',
           'StudentT',
           'Triangular',
           'Uniform',
           'UniformSum',
           'VonMises',
           'Weibull',
           'WignerSemicircle')


def ContinuousRV(symbol, density, set=Interval(-oo, oo, True, True)):
    """
    Create a Continuous Random Variable given the following:

    -- a symbol
    -- a probability density function
    -- set on which the pdf is valid (defaults to entire real line)

    Returns a RandomSymbol.

    Many common continuous random variable types are already implemented.
    This function should be necessary only very rarely.

    Examples
    ========

    >>> from diofant.stats import E, P

    >>> pdf = sqrt(2)*exp(-x**2/2)/(2*sqrt(pi))  # Normal distribution
    >>> X = ContinuousRV(x, pdf)

    >>> E(X)
    0
    >>> P(X > 0)
    1/2

    """
    pdf = Lambda(symbol, density)
    dist = ContinuousDistributionHandmade(pdf, set)
    return SingleContinuousPSpace(symbol, dist).value


def rv(symbol, cls, args):
    args = list(map(sympify, args))
    dist = cls(*args)
    dist.check(*args)
    return SingleContinuousPSpace(symbol, dist).value

########################################
# Continuous Probability Distributions #
########################################

# ------------------------------------------------------------------------------
# Arcsin distribution ----------------------------------------------------------


class ArcsinDistribution(SingleContinuousDistribution):
    _argnames = ('a', 'b')

    def pdf(self, x):
        return 1/(pi*sqrt((x - self.a)*(self.b - x)))


def Arcsin(name, a=0, b=1):
    r"""
    Create a Continuous Random Variable with an arcsin distribution.

    The density of the arcsin distribution is given by

    .. math::
        f(x) := \frac{1}{\pi\sqrt{(x-a)(b-x)}}

    with `x \in [a,b]`. It must hold that `-\infty < a < b < \infty`.

    Parameters
    ==========

    a : Real number, the left interval boundary
    b : Real number, the right interval boundary

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> a, b = symbols('a b', real=True)

    >>> X = Arcsin('x', a, b)
    >>> density(X)(z)
    1/(pi*sqrt((-a + z)*(b - z)))

    References
    ==========

    * https://en.wikipedia.org/wiki/Arcsine_distribution

    """
    return rv(name, ArcsinDistribution, (a, b))

# ------------------------------------------------------------------------------
# Benini distribution ----------------------------------------------------------


class BeniniDistribution(SingleContinuousDistribution):
    _argnames = ('alpha', 'beta', 'sigma')

    @property
    def set(self):
        return Interval(self.sigma, oo, False, True)

    def pdf(self, x):
        alpha, beta, sigma = self.alpha, self.beta, self.sigma
        return (exp(-alpha*log(x/sigma) - beta*log(x/sigma)**2)
                * (alpha/x + 2*beta*log(x/sigma)/x))


def Benini(name, alpha, beta, sigma):
    r"""
    Create a Continuous Random Variable with a Benini distribution.

    The density of the Benini distribution is given by

    .. math::
        f(x) := e^{-\alpha\log{\frac{x}{\sigma}}
                -\beta\log^2\left[{\frac{x}{\sigma}}\right]}
                \left(\frac{\alpha}{x}+\frac{2\beta\log{\frac{x}{\sigma}}}{x}\right)

    This is a heavy-tailed distribution and is also known as the log-Rayleigh
    distribution.

    Parameters
    ==========

    alpha : Real number, `\alpha > 0`, a shape
    beta : Real number, `\beta > 0`, a shape
    sigma : Real number, `\sigma > 0`, a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> alpha = Symbol('alpha', positive=True)
    >>> beta = Symbol('beta', positive=True)
    >>> sigma = Symbol('sigma', positive=True)

    >>> X = Benini('x', alpha, beta, sigma)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
                /  z  \           2/  z  \ /                  /  z  \\
     - alpha*log|-----| - beta*log |-----| |        2*beta*log|-----||
                \sigma/            \sigma/ |alpha             \sigma/|
    E                                     *|----- + -----------------|
                                           \  z             z        /

    References
    ==========

    * https://en.wikipedia.org/wiki/Benini_distribution
    * https://reference.wolfram.com/legacy/v8/ref/BeniniDistribution.html

    """
    return rv(name, BeniniDistribution, (alpha, beta, sigma))

# ------------------------------------------------------------------------------
# Beta distribution ------------------------------------------------------------


class BetaDistribution(SingleContinuousDistribution):
    _argnames = ('alpha', 'beta')

    set = Interval(0, 1)

    @staticmethod
    def check(alpha, beta):
        _value_check(alpha > 0, 'Alpha must be positive')
        _value_check(beta > 0, 'Beta must be positive')

    def pdf(self, x):
        alpha, beta = self.alpha, self.beta
        return x**(alpha - 1) * (1 - x)**(beta - 1) / beta_fn(alpha, beta)

    def sample(self):
        return random.betavariate(self.alpha, self.beta)


def Beta(name, alpha, beta):
    r"""
    Create a Continuous Random Variable with a Beta distribution.

    The density of the Beta distribution is given by

    .. math::
        f(x) := \frac{x^{\alpha-1}(1-x)^{\beta-1}} {\mathrm{B}(\alpha,\beta)}

    with `x \in [0,1]`.

    Parameters
    ==========

    alpha : Real number, `\alpha > 0`, a shape
    beta : Real number, `\beta > 0`, a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, density

    >>> alpha = Symbol('alpha', positive=True)
    >>> beta = Symbol('beta', positive=True)

    >>> X = Beta('x', alpha, beta)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
     alpha - 1         beta - 1
    z         *(-z + 1)
    ---------------------------
         beta(alpha, beta)

    >>> expand_func(simplify(E(X, meijerg=True)))
    alpha/(alpha + beta)

    References
    ==========

    * https://en.wikipedia.org/wiki/Beta_distribution
    * https://mathworld.wolfram.com/BetaDistribution.html

    """
    return rv(name, BetaDistribution, (alpha, beta))

# ------------------------------------------------------------------------------
# Beta prime distribution ------------------------------------------------------


class BetaPrimeDistribution(SingleContinuousDistribution):
    _argnames = ('alpha', 'beta')

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        alpha, beta = self.alpha, self.beta
        return x**(alpha - 1)*(1 + x)**(-alpha - beta)/beta_fn(alpha, beta)


def BetaPrime(name, alpha, beta):
    r"""
    Create a continuous random variable with a Beta prime distribution.

    The density of the Beta prime distribution is given by

    .. math::
        f(x) := \frac{x^{\alpha-1} (1+x)^{-\alpha -\beta}}{B(\alpha,\beta)}

    with `x > 0`.

    Parameters
    ==========

    alpha : Real number, `\alpha > 0`, a shape
    beta : Real number, `\beta > 0`, a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> alpha = Symbol('alpha', positive=True)
    >>> beta = Symbol('beta', positive=True)

    >>> X = BetaPrime('x', alpha, beta)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
     alpha - 1        -alpha - beta
    z         *(z + 1)
    -------------------------------
           beta(alpha, beta)

    References
    ==========

    * https://en.wikipedia.org/wiki/Beta_prime_distribution
    * https://mathworld.wolfram.com/BetaPrimeDistribution.html

    """
    return rv(name, BetaPrimeDistribution, (alpha, beta))

# ------------------------------------------------------------------------------
# Cauchy distribution ----------------------------------------------------------


class CauchyDistribution(SingleContinuousDistribution):
    _argnames = ('x0', 'gamma')

    def pdf(self, x):
        return 1/(pi*self.gamma*(1 + ((x - self.x0)/self.gamma)**2))


def Cauchy(name, x0, gamma):
    r"""
    Create a continuous random variable with a Cauchy distribution.

    The density of the Cauchy distribution is given by

    .. math::
        f(x) := \frac{1}{\pi} \arctan\left(\frac{x-x_0}{\gamma}\right)
                +\frac{1}{2}

    Parameters
    ==========

    x0 : Real number, the location
    gamma : Real number, `\gamma > 0`, the scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> x0 = Symbol('x0')
    >>> gamma = Symbol('gamma', positive=True)

    >>> X = Cauchy('x', x0, gamma)

    >>> density(X)(z)
    1/(pi*gamma*(1 + (-x0 + z)**2/gamma**2))

    References
    ==========

    * https://en.wikipedia.org/wiki/Cauchy_distribution
    * https://mathworld.wolfram.com/CauchyDistribution.html

    """
    return rv(name, CauchyDistribution, (x0, gamma))

# ------------------------------------------------------------------------------
# Chi distribution -------------------------------------------------------------


class ChiDistribution(SingleContinuousDistribution):
    _argnames = 'k',

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        return 2**(1 - self.k/2)*x**(self.k - 1)*exp(-x**2/2)/gamma(self.k/2)


def Chi(name, k):
    r"""
    Create a continuous random variable with a Chi distribution.

    The density of the Chi distribution is given by

    .. math::
        f(x) := \frac{2^{1-k/2}x^{k-1}e^{-x^2/2}}{\Gamma(k/2)}

    with `x \geq 0`.

    Parameters
    ==========

    k : A positive Integer, `k > 0`, the number of degrees of freedom

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import Chi, density

    >>> X = Chi('x', k)

    >>> density(X)(z)
    2**(-k/2 + 1)*E**(-z**2/2)*z**(k - 1)/gamma(k/2)

    References
    ==========

    * https://en.wikipedia.org/wiki/Chi_distribution
    * https://mathworld.wolfram.com/ChiDistribution.html

    """
    return rv(name, ChiDistribution, (k,))

# ------------------------------------------------------------------------------
# Non-central Chi distribution -------------------------------------------------


class ChiNoncentralDistribution(SingleContinuousDistribution):
    _argnames = ('k', 'l')

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        k, l = self.k, self.l
        return exp(-(x**2+l**2)/2)*x**k*l / (l*x)**(k/2) * besseli(k/2-1, l*x)


def ChiNoncentral(name, k, l):
    r"""
    Create a continuous random variable with a non-central Chi distribution.

    The density of the non-central Chi distribution is given by

    .. math::
        f(x) := \frac{e^{-(x^2+\lambda^2)/2} x^k\lambda}
                {(\lambda x)^{k/2}} I_{k/2-1}(\lambda x)

    with `x \geq 0`. Here, `I_\nu (x)` is the
    :ref:`modified Bessel function of the first kind <besseli>`.

    Parameters
    ==========

    k : A positive Integer, `k > 0`, the number of degrees of freedom
    l : Shift parameter

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> l = Symbol('l')

    >>> X = ChiNoncentral('x', k, l)

    >>> density(X)(z)
    E**(-l**2/2 - z**2/2)*l*z**k*(l*z)**(-k/2)*besseli(k/2 - 1, l*z)

    References
    ==========

    * https://en.wikipedia.org/wiki/Noncentral_chi_distribution

    """
    return rv(name, ChiNoncentralDistribution, (k, l))

# ------------------------------------------------------------------------------
# Chi squared distribution -----------------------------------------------------


class ChiSquaredDistribution(SingleContinuousDistribution):
    _argnames = 'k',

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        k = self.k
        return 1/(2**(k/2)*gamma(k/2))*x**(k/2 - 1)*exp(-x/2)


def ChiSquared(name, k):
    r"""
    Create a continuous random variable with a Chi-squared distribution.

    The density of the Chi-squared distribution is given by

    .. math::
        f(x) := \frac{1}{2^{\frac{k}{2}}\Gamma\left(\frac{k}{2}\right)}
                x^{\frac{k}{2}-1} e^{-\frac{x}{2}}

    with `x \geq 0`.

    Parameters
    ==========

    k : A positive Integer, `k > 0`, the number of degrees of freedom

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, density, variance

    >>> k = Symbol('k', integer=True, positive=True)

    >>> X = ChiSquared('x', k)

    >>> density(X)(z)
    2**(-k/2)*E**(-z/2)*z**(k/2 - 1)/gamma(k/2)

    >>> combsimp(E(X))
    k

    >>> simplify(expand_func(variance(X)))
    2*k

    References
    ==========

    * https://en.wikipedia.org/wiki/Chi_squared_distribution
    * https://mathworld.wolfram.com/Chi-SquaredDistribution.html

    """
    return rv(name, ChiSquaredDistribution, (k, ))

# ------------------------------------------------------------------------------
# Dagum distribution -----------------------------------------------------------


class DagumDistribution(SingleContinuousDistribution):
    _argnames = ('p', 'a', 'b')

    def pdf(self, x):
        p, a, b = self.p, self.a, self.b
        return a*p/x*((x/b)**(a*p)/(((x/b)**a + 1)**(p + 1)))


def Dagum(name, p, a, b):
    r"""
    Create a continuous random variable with a Dagum distribution.

    The density of the Dagum distribution is given by

    .. math::
        f(x) := \frac{a p}{x} \left( \frac{\left(\tfrac{x}{b}\right)^{a p}}
                {\left(\left(\tfrac{x}{b}\right)^a + 1 \right)^{p+1}} \right)

    with `x > 0`.

    Parameters
    ==========

    p : Real number, `p > 0`, a shape
    a : Real number, `a > 0`, a shape
    b : Real number, `b > 0`, a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> p = Symbol('p', positive=True)
    >>> b = Symbol('b', positive=True)
    >>> a = Symbol('a', positive=True)

    >>> X = Dagum('x', p, a, b)

    >>> density(X)(z)
    a*p*(z/b)**(a*p)*((z/b)**a + 1)**(-p - 1)/z

    References
    ==========

    * https://en.wikipedia.org/wiki/Dagum_distribution

    """
    return rv(name, DagumDistribution, (p, a, b))

# ------------------------------------------------------------------------------
# Erlang distribution ----------------------------------------------------------


def Erlang(name, k, l):
    r"""
    Create a continuous random variable with an Erlang distribution.

    The density of the Erlang distribution is given by

    .. math::
        f(x) := \frac{\lambda^k x^{k-1} e^{-\lambda x}}{(k-1)!}

    with `x \in [0,\infty]`.

    Parameters
    ==========

    k : Integer
    l : Real number, `\lambda > 0`, the rate

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, cdf, density, variance

    >>> k = Symbol('k', integer=True, positive=True)
    >>> l = Symbol('l', positive=True)

    >>> X = Erlang('x', k, l)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
     -l*z  k  k - 1
    E    *l *z
    ---------------
        gamma(k)

    >>> C = cdf(X, meijerg=True)(z)
    >>> pprint(C, use_unicode=False)
    /  k*lowergamma(k, 0)   k*lowergamma(k, l*z)
    |- ------------------ + --------------------  for z >= 0
    <     gamma(k + 1)          gamma(k + 1)
    |
    \                     0                       otherwise

    >>> simplify(E(X))
    k/l

    >>> simplify(variance(X))
    k/l**2

    References
    ==========

    * https://en.wikipedia.org/wiki/Erlang_distribution
    * https://mathworld.wolfram.com/ErlangDistribution.html

    """
    return rv(name, GammaDistribution, (k, 1/l))

# ------------------------------------------------------------------------------
# Exponential distribution -----------------------------------------------------


class ExponentialDistribution(SingleContinuousDistribution):
    _argnames = 'rate',

    set = Interval(0, oo, False, True)

    @staticmethod
    def check(rate):
        _value_check(rate > 0, 'Rate must be positive.')

    def pdf(self, x):
        return self.rate * exp(-self.rate*x)

    def sample(self):
        return random.expovariate(self.rate)


def Exponential(name, rate):
    r"""
    Create a continuous random variable with an Exponential distribution.

    The density of the exponential distribution is given by

    .. math::
        f(x) := \lambda \exp(-\lambda x)

    with `x > 0`. Note that the expected value is `1/\lambda`.

    Parameters
    ==========

    rate : A positive Real number, `\lambda > 0`, the rate (or inverse scale/inverse mean)

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, cdf, density, skewness, std, variance

    >>> l = Symbol('lambda', positive=True)

    >>> X = Exponential('x', l)

    >>> density(X)(z)
    E**(-lambda*z)*lambda

    >>> cdf(X)(z)
    Piecewise((1 - E**(-lambda*z), z >= 0), (0, true))

    >>> E(X)
    1/lambda

    >>> variance(X)
    lambda**(-2)

    >>> skewness(X)
    2

    >>> X = Exponential('x', 10)

    >>> density(X)(z)
    10*E**(-10*z)

    >>> E(X)
    1/10

    >>> std(X)
    1/10

    References
    ==========

    * https://en.wikipedia.org/wiki/Exponential_distribution
    * https://mathworld.wolfram.com/ExponentialDistribution.html

    """
    return rv(name, ExponentialDistribution, (rate, ))

# ------------------------------------------------------------------------------
# F distribution ---------------------------------------------------------------


class FDistributionDistribution(SingleContinuousDistribution):
    _argnames = ('d1', 'd2')

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        d1, d2 = self.d1, self.d2
        return (sqrt((d1*x)**d1*d2**d2 / (d1*x+d2)**(d1+d2))
                / (x * beta_fn(d1/2, d2/2)))


def FDistribution(name, d1, d2):
    r"""
    Create a continuous random variable with a F distribution.

    The density of the F distribution is given by

    .. math::
        f(x) := \frac{\sqrt{\frac{(d_1 x)^{d_1} d_2^{d_2}}
                {(d_1 x + d_2)^{d_1 + d_2}}}}
                {x \mathrm{B} \left(\frac{d_1}{2}, \frac{d_2}{2}\right)}

    with `x > 0`.

    .. TODO - What do these parameters mean?

    Parameters
    ==========

    d1 : `d_1 > 0` a parameter
    d2 : `d_2 > 0` a parameter

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> d1 = Symbol('d1', positive=True)
    >>> d2 = Symbol('d2', positive=True)

    >>> X = FDistribution('x', d1, d2)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
      d2
      --    ______________________________
      2    /       d1            -d1 - d2
    d2  *\/  (d1*z)  *(d1*z + d2)
    --------------------------------------
                      /d1  d2\
                z*beta|--, --|
                      \2   2 /

    References
    ==========

    * https://en.wikipedia.org/wiki/F-distribution
    * https://mathworld.wolfram.com/F-Distribution.html

    """
    return rv(name, FDistributionDistribution, (d1, d2))

# ------------------------------------------------------------------------------
# Fisher Z distribution --------------------------------------------------------


class FisherZDistribution(SingleContinuousDistribution):
    _argnames = ('d1', 'd2')

    def pdf(self, x):
        d1, d2 = self.d1, self.d2
        return (2*d1**(d1/2)*d2**(d2/2) / beta_fn(d1/2, d2/2) *
                exp(d1*x) / (d1*exp(2*x)+d2)**((d1+d2)/2))


def FisherZ(name, d1, d2):
    r"""
    Create a Continuous Random Variable with an Fisher's Z distribution.

    The density of the Fisher's Z distribution is given by

    .. math::
        f(x) := \frac{2d_1^{d_1/2} d_2^{d_2/2}} {\mathrm{B}(d_1/2, d_2/2)}
                \frac{e^{d_1z}}{\left(d_1e^{2z}+d_2\right)^{\left(d_1+d_2\right)/2}}


    .. TODO - What is the difference between these degrees of freedom?

    Parameters
    ==========

    d1 : `d_1 > 0`, degree of freedom
    d2 : `d_2 > 0`, degree of freedom

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> d1 = Symbol('d1', positive=True)
    >>> d2 = Symbol('d2', positive=True)

    >>> X = FisherZ('x', d1, d2)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
                                      d1   d2
              d1   d2               - -- - --
              --   --                 2    2
       d1*z   2    2  / 2*z        \
    2*E    *d1  *d2  *\E   *d1 + d2/
    -----------------------------------------
                       /d1  d2\
                   beta|--, --|
                       \2   2 /

    References
    ==========

    * https://en.wikipedia.org/wiki/Fisher%27s_z-distribution
    * https://mathworld.wolfram.com/Fishersz-Distribution.html

    """
    return rv(name, FisherZDistribution, (d1, d2))

# ------------------------------------------------------------------------------
# Frechet distribution ---------------------------------------------------------


class FrechetDistribution(SingleContinuousDistribution):
    _argnames = ('a', 's', 'm')

    set = Interval(0, oo, False, True)

    def __new__(cls, a, s=1, m=0):
        a, s, m = list(map(sympify, (a, s, m)))
        return Expr.__new__(cls, a, s, m)

    def pdf(self, x):
        a, s, m = self.a, self.s, self.m
        return a/s * ((x-m)/s)**(-1-a) * exp(-((x-m)/s)**(-a))


def Frechet(name, a, s=1, m=0):
    r"""
    Create a continuous random variable with a Frechet distribution.

    The density of the Frechet distribution is given by

    .. math::
        f(x) := \frac{\alpha}{s} \left(\frac{x-m}{s}\right)^{-1-\alpha}
                 e^{-(\frac{x-m}{s})^{-\alpha}}

    with `x \geq m`.

    Parameters
    ==========

    a : Real number, `a \in \left(0, \infty\right)` the shape
    s : Real number, `s \in \left(0, \infty\right)` the scale
    m : Real number, `m \in \left(-\infty, \infty\right)` the minimum

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> a, s = symbols('a s', positive=True)
    >>> m = Symbol('m', real=True)

    >>> X = Frechet('x', a, s, m)

    >>> density(X)(z)
    E**(-((-m + z)/s)**(-a))*a*((-m + z)/s)**(-a - 1)/s

    References
    ==========

    * https://en.wikipedia.org/wiki/Fr%C3%A9chet_distribution

    """
    return rv(name, FrechetDistribution, (a, s, m))

# ------------------------------------------------------------------------------
# Gamma distribution -----------------------------------------------------------


class GammaDistribution(SingleContinuousDistribution):
    _argnames = ('k', 'theta')

    set = Interval(0, oo, False, True)

    @staticmethod
    def check(k, theta):
        _value_check(k > 0, 'k must be positive')
        _value_check(theta > 0, 'Theta must be positive')

    def pdf(self, x):
        k, theta = self.k, self.theta
        return x**(k - 1) * exp(-x/theta) / (gamma(k)*theta**k)

    def sample(self):
        return random.gammavariate(self.k, self.theta)


def Gamma(name, k, theta):
    r"""
    Create a continuous random variable with a Gamma distribution.

    The density of the Gamma distribution is given by

    .. math::
        f(x) := \frac{1}{\Gamma(k) \theta^k} x^{k - 1} e^{-\frac{x}{\theta}}

    with `x \in [0,1]`.

    Parameters
    ==========

    k : Real number, `k > 0`, a shape
    theta : Real number, `\theta > 0`, a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, cdf, density, variance

    >>> k = Symbol('k', positive=True)
    >>> theta = Symbol('theta', positive=True)

    >>> X = Gamma('x', k, theta)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
      -z
     -----
     theta      -k  k - 1
    E     *theta  *z
    ---------------------
           gamma(k)

    >>> C = cdf(X, meijerg=True)(z)
    >>> pprint(C, use_unicode=False)
    /                                   /     z  \
    |                       k*lowergamma|k, -----|
    |  k*lowergamma(k, 0)               \   theta/
    <- ------------------ + ----------------------  for z >= 0
    |     gamma(k + 1)           gamma(k + 1)
    |
    \                      0                        otherwise

    >>> E(X)
    theta*gamma(k + 1)/gamma(k)

    >>> V = simplify(variance(X))
    >>> pprint(V, use_unicode=False)
           2
    k*theta

    References
    ==========

    * https://en.wikipedia.org/wiki/Gamma_distribution
    * https://mathworld.wolfram.com/GammaDistribution.html

    """
    return rv(name, GammaDistribution, (k, theta))

# ------------------------------------------------------------------------------
# Inverse Gamma distribution ---------------------------------------------------


class GammaInverseDistribution(SingleContinuousDistribution):
    _argnames = ('a', 'b')

    set = Interval(0, oo, False, True)

    @staticmethod
    def check(a, b):
        _value_check(a > 0, 'alpha must be positive')
        _value_check(b > 0, 'beta must be positive')

    def pdf(self, x):
        a, b = self.a, self.b
        return b**a/gamma(a) * x**(-a-1) * exp(-b/x)


def GammaInverse(name, a, b):
    r"""
    Create a continuous random variable with an inverse Gamma distribution.

    The density of the inverse Gamma distribution is given by

    .. math::
        f(x) := \frac{\beta^\alpha}{\Gamma(\alpha)} x^{-\alpha - 1}
                \exp\left(\frac{-\beta}{x}\right)

    with `x > 0`.

    Parameters
    ==========

    a : Real number, `a > 0` a shape
    b : Real number, `b > 0` a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> a = Symbol('a', positive=True)
    >>> b = Symbol('b', positive=True)

    >>> X = GammaInverse('x', a, b)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
     -b
     ---
      z   a  -a - 1
    E   *b *z
    ---------------
        gamma(a)

    References
    ==========

    * https://en.wikipedia.org/wiki/Inverse-gamma_distribution

    """
    return rv(name, GammaInverseDistribution, (a, b))

# ------------------------------------------------------------------------------
# Kumaraswamy distribution -----------------------------------------------------


class KumaraswamyDistribution(SingleContinuousDistribution):
    _argnames = ('a', 'b')

    set = Interval(0, oo, False, True)

    @staticmethod
    def check(a, b):
        _value_check(a > 0, 'a must be positive')
        _value_check(b > 0, 'b must be positive')

    def pdf(self, x):
        a, b = self.a, self.b
        return a * b * x**(a-1) * (1-x**a)**(b-1)


def Kumaraswamy(name, a, b):
    r"""
    Create a Continuous Random Variable with a Kumaraswamy distribution.

    The density of the Kumaraswamy distribution is given by

    .. math::
        f(x) := a b x^{a-1} (1-x^a)^{b-1}

    with `x \in [0,1]`.

    Parameters
    ==========

    a : Real number, `a > 0` a shape
    b : Real number, `b > 0` a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> a = Symbol('a', positive=True)
    >>> b = Symbol('b', positive=True)

    >>> X = Kumaraswamy('x', a, b)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
                         b - 1
         a - 1 /   a    \
    a*b*z     *\- z  + 1/


    References
    ==========

    * https://en.wikipedia.org/wiki/Kumaraswamy_distribution

    """
    return rv(name, KumaraswamyDistribution, (a, b))

# ------------------------------------------------------------------------------
# Laplace distribution ---------------------------------------------------------


class LaplaceDistribution(SingleContinuousDistribution):
    _argnames = ('mu', 'b')

    def pdf(self, x):
        mu, b = self.mu, self.b
        return 1/(2*b)*exp(-Abs(x - mu)/b)


def Laplace(name, mu, b):
    r"""
    Create a continuous random variable with a Laplace distribution.

    The density of the Laplace distribution is given by

    .. math::
        f(x) := \frac{1}{2 b} \exp \left(-\frac{|x-\mu|}b \right)

    Parameters
    ==========

    mu : Real number, the location (mean)
    b : Real number, `b > 0`, a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> mu = Symbol('mu')
    >>> b = Symbol('b', positive=True)

    >>> X = Laplace('x', mu, b)

    >>> density(X)(z)
    E**(-Abs(mu - z)/b)/(2*b)

    References
    ==========

    * https://en.wikipedia.org/wiki/Laplace_distribution
    * https://mathworld.wolfram.com/LaplaceDistribution.html

    """
    return rv(name, LaplaceDistribution, (mu, b))

# ------------------------------------------------------------------------------
# Logistic distribution --------------------------------------------------------


class LogisticDistribution(SingleContinuousDistribution):
    _argnames = ('mu', 's')

    def pdf(self, x):
        mu, s = self.mu, self.s
        return exp(-(x - mu)/s)/(s*(1 + exp(-(x - mu)/s))**2)


def Logistic(name, mu, s):
    r"""
    Create a continuous random variable with a logistic distribution.

    The density of the logistic distribution is given by

    .. math::
        f(x) := \frac{e^{-(x-\mu)/s}} {s\left(1+e^{-(x-\mu)/s}\right)^2}

    Parameters
    ==========

    mu : Real number, the location (mean)
    s : Real number, `s > 0` a scale

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> mu = Symbol('mu', real=True)
    >>> s = Symbol('s', positive=True)

    >>> X = Logistic('x', mu, s)

    >>> density(X)(z)
    E**((mu - z)/s)/(s*(E**((mu - z)/s) + 1)**2)

    References
    ==========

    * https://en.wikipedia.org/wiki/Logistic_distribution
    * https://mathworld.wolfram.com/LogisticDistribution.html

    """
    return rv(name, LogisticDistribution, (mu, s))

# ------------------------------------------------------------------------------
# Log Normal distribution ------------------------------------------------------


class LogNormalDistribution(SingleContinuousDistribution):
    _argnames = ('mean', 'std')

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        mean, std = self.mean, self.std
        return exp(-(log(x) - mean)**2 / (2*std**2)) / (x*sqrt(2*pi)*std)

    def sample(self):
        return random.lognormvariate(self.mean, self.std)


def LogNormal(name, mean, std):
    r"""
    Create a continuous random variable with a log-normal distribution.

    The density of the log-normal distribution is given by

    .. math::
        f(x) := \frac{1}{x\sqrt{2\pi\sigma^2}}
                e^{-\frac{\left(\ln x-\mu\right)^2}{2\sigma^2}}

    with `x \geq 0`.

    Parameters
    ==========

    mu : Real number, the log-scale
    sigma : Real number, `\sigma^2 > 0` a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> mu = Symbol('mu', real=True)
    >>> sigma = Symbol('sigma', positive=True)

    >>> X = LogNormal('x', mu, sigma)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
                          2
           -(-mu + log(z))
           -----------------
                       2
      ___       2*sigma
    \/ 2 *E
    ------------------------
            ____
        2*\/ pi *sigma*z

    >>> X = LogNormal('x', 0, 1)  # Mean 0, standard deviation 1

    >>> density(X)(z)
    sqrt(2)*E**(-log(z)**2/2)/(2*sqrt(pi)*z)

    References
    ==========

    * https://en.wikipedia.org/wiki/Lognormal
    * https://mathworld.wolfram.com/LogNormalDistribution.html

    """
    return rv(name, LogNormalDistribution, (mean, std))

# ------------------------------------------------------------------------------
# Maxwell distribution ---------------------------------------------------------


class MaxwellDistribution(SingleContinuousDistribution):
    _argnames = 'a',

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        a = self.a
        return sqrt(2/pi)*x**2*exp(-x**2/(2*a**2))/a**3


def Maxwell(name, a):
    r"""
    Create a continuous random variable with a Maxwell distribution.

    The density of the Maxwell distribution is given by

    .. math::
        f(x) := \sqrt{\frac{2}{\pi}} \frac{x^2 e^{-x^2/(2a^2)}}{a^3}

    with `x \geq 0`.

    .. TODO - what does the parameter mean?

    Parameters
    ==========

    a : Real number, `a > 0`

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, density, variance

    >>> a = Symbol('a', positive=True)

    >>> X = Maxwell('x', a)

    >>> density(X)(z)
    sqrt(2)*E**(-z**2/(2*a**2))*z**2/(sqrt(pi)*a**3)

    >>> E(X)
    2*sqrt(2)*a/sqrt(pi)

    >>> simplify(variance(X))
    a**2*(-8 + 3*pi)/pi

    References
    ==========

    * https://en.wikipedia.org/wiki/Maxwell_distribution
    * https://mathworld.wolfram.com/MaxwellDistribution.html

    """
    return rv(name, MaxwellDistribution, (a, ))

# ------------------------------------------------------------------------------
# Nakagami distribution --------------------------------------------------------


class NakagamiDistribution(SingleContinuousDistribution):
    _argnames = ('mu', 'omega')

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        mu, omega = self.mu, self.omega
        return 2*mu**mu/(gamma(mu)*omega**mu)*x**(2*mu - 1)*exp(-mu/omega*x**2)


def Nakagami(name, mu, omega):
    r"""
    Create a continuous random variable with a Nakagami distribution.

    The density of the Nakagami distribution is given by

    .. math::
        f(x) := \frac{2\mu^\mu}{\Gamma(\mu)\omega^\mu} x^{2\mu-1}
                \exp\left(-\frac{\mu}{\omega}x^2 \right)

    with `x > 0`.

    Parameters
    ==========

    mu : Real number, `\mu \geq \frac{1}{2}` a shape
    omega : Real number, `\omega > 0`, the spread

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, density, variance

    >>> mu = Symbol('mu', positive=True)
    >>> omega = Symbol('omega', positive=True)

    >>> X = Nakagami('x', mu, omega)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
            2
       -mu*z
       -------
        omega    mu      -mu  2*mu - 1
    2*E       *mu  *omega   *z
    ----------------------------------
                gamma(mu)

    >>> simplify(E(X, meijerg=True))
    sqrt(mu)*sqrt(omega)*gamma(mu + 1/2)/gamma(mu + 1)

    >>> V = simplify(variance(X, meijerg=True))
    >>> pprint(V, use_unicode=False)
                        2
             omega*gamma (mu + 1/2)
    omega - -----------------------
            gamma(mu)*gamma(mu + 1)

    References
    ==========

    * https://en.wikipedia.org/wiki/Nakagami_distribution

    """
    return rv(name, NakagamiDistribution, (mu, omega))

# ------------------------------------------------------------------------------
# Normal distribution ----------------------------------------------------------


class NormalDistribution(SingleContinuousDistribution):
    _argnames = ('mean', 'std')

    @staticmethod
    def check(mean, std):
        _value_check(std > 0, 'Standard deviation must be positive')

    def pdf(self, x):
        return exp(-(x - self.mean)**2 / (2*self.std**2)) / (sqrt(2*pi)*self.std)

    def sample(self):
        return random.normalvariate(self.mean, self.std)


def Normal(name, mean, std):
    r"""
    Create a continuous random variable with a Normal distribution.

    The density of the Normal distribution is given by

    .. math::
        f(x) := \frac{1}{\sigma\sqrt{2\pi}} e^{ -\frac{(x-\mu)^2}{2\sigma^2} }

    Parameters
    ==========

    mu : Real number, the mean
    sigma : Real number, `\sigma^2 > 0` the variance

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, cdf, density, skewness, std

    >>> mu = Symbol('mu')
    >>> sigma = Symbol('sigma', positive=True)

    >>> X = Normal('x', mu, sigma)

    >>> density(X)(z)
    sqrt(2)*E**(-(-mu + z)**2/(2*sigma**2))/(2*sqrt(pi)*sigma)

    >>> C = simplify(expand(cdf(X)))(z)  # it needs a little more help...
    >>> pprint(C, use_unicode=False)
         /  ___         \
         |\/ 2 *(mu - z)|
      erf|--------------|
         \    2*sigma   /   1
    - ------------------- + -
               2            2

    >>> simplify(skewness(X))
    0

    >>> X = Normal('x', 0, 1)  # Mean 0, standard deviation 1
    >>> density(X)(z)
    sqrt(2)*E**(-z**2/2)/(2*sqrt(pi))

    >>> E(2*X + 1)
    1

    >>> simplify(std(2*X + 1))
    2

    References
    ==========

    * https://en.wikipedia.org/wiki/Normal_distribution
    * https://mathworld.wolfram.com/NormalDistributionFunction.html

    """
    return rv(name, NormalDistribution, (mean, std))

# ------------------------------------------------------------------------------
# Pareto distribution ----------------------------------------------------------


class ParetoDistribution(SingleContinuousDistribution):
    _argnames = ('xm', 'alpha')

    @property
    def set(self):
        return Interval(self.xm, oo, False, True)

    @staticmethod
    def check(xm, alpha):
        _value_check(xm > 0, 'Xm must be positive')
        _value_check(alpha > 0, 'Alpha must be positive')

    def pdf(self, x):
        xm, alpha = self.xm, self.alpha
        return alpha * xm**alpha / x**(alpha + 1)

    def sample(self):
        return random.paretovariate(self.alpha)


def Pareto(name, xm, alpha):
    r"""
    Create a continuous random variable with the Pareto distribution.

    The density of the Pareto distribution is given by

    .. math::
        f(x) := \frac{\alpha\,x_m^\alpha}{x^{\alpha+1}}

    with `x \in [x_m,\infty]`.

    Parameters
    ==========

    xm : Real number, `x_m > 0`, a scale
    alpha : Real number, `\alpha > 0`, a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> xm = Symbol('xm', positive=True)
    >>> beta = Symbol('beta', positive=True)

    >>> X = Pareto('x', xm, beta)

    >>> density(X)(z)
    beta*xm**beta*z**(-beta - 1)

    References
    ==========

    * https://en.wikipedia.org/wiki/Pareto_distribution
    * https://mathworld.wolfram.com/ParetoDistribution.html

    """
    return rv(name, ParetoDistribution, (xm, alpha))

# ------------------------------------------------------------------------------
# QuadraticU distribution ------------------------------------------------------


class QuadraticUDistribution(SingleContinuousDistribution):
    _argnames = ('a', 'b')

    @property
    def set(self):
        return Interval(self.a, self.b)

    def pdf(self, x):
        a, b = self.a, self.b
        alpha = 12/(b - a)**3
        beta = (a + b)/2
        return Piecewise((alpha*(x - beta)**2, And(a <= x, x <= b)),
                         (0, True))


def QuadraticU(name, a, b):
    r"""
    Create a Continuous Random Variable with a U-quadratic distribution.

    The density of the U-quadratic distribution is given by

    .. math::
        f(x) := \alpha (x-\beta)^2

    with `x \in [a,b]`.

    Parameters
    ==========

    a : Real number
    b : Real number, `a < b`

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> a, b = symbols('a b', real=True)

    >>> X = QuadraticU('x', a, b)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
    /                2
    |   /  a   b    \
    |12*|- - - - + z|
    |   \  2   2    /
    <-----------------  for And(a <= z, z <= b)
    |            3
    |    (-a + b)
    |
    \        0                 otherwise

    References
    ==========

    * https://en.wikipedia.org/wiki/U-quadratic_distribution

    """
    return rv(name, QuadraticUDistribution, (a, b))

# ------------------------------------------------------------------------------
# RaisedCosine distribution ----------------------------------------------------


class RaisedCosineDistribution(SingleContinuousDistribution):
    _argnames = ('mu', 's')

    @property
    def set(self):
        return Interval(self.mu - self.s, self.mu + self.s)

    @staticmethod
    def check(mu, s):
        _value_check(s > 0, 's must be positive')

    def pdf(self, x):
        mu, s = self.mu, self.s
        return Piecewise(
            ((1 + cos(pi*(x - mu)/s))/(2*s), And(mu - s <= x, x <= mu + s)),
            (0, True))


def RaisedCosine(name, mu, s):
    r"""
    Create a Continuous Random Variable with a raised cosine distribution.

    The density of the raised cosine distribution is given by

    .. math::
        f(x) := \frac{1}{2s}\left(1+\cos\left(\frac{x-\mu}{s}\pi\right)\right)

    with `x \in [\mu-s,\mu+s]`.

    Parameters
    ==========

    mu : Real number
    s : Real number, `s > 0`

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> mu = Symbol('mu', real=True)
    >>> s = Symbol('s', positive=True)

    >>> X = RaisedCosine('x', mu, s)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
    /   /pi*(-mu + z)\
    |cos|------------| + 1
    |   \     s      /
    <---------------------  for And(z <= mu + s, mu - s <= z)
    |         2*s
    |
    \          0                        otherwise

    References
    ==========

    * https://en.wikipedia.org/wiki/Raised_cosine_distribution

    """
    return rv(name, RaisedCosineDistribution, (mu, s))

# ------------------------------------------------------------------------------
# Rayleigh distribution --------------------------------------------------------


class RayleighDistribution(SingleContinuousDistribution):
    _argnames = 'sigma',

    set = Interval(0, oo, False, True)

    def pdf(self, x):
        sigma = self.sigma
        return x/sigma**2*exp(-x**2/(2*sigma**2))


def Rayleigh(name, sigma):
    r"""
    Create a continuous random variable with a Rayleigh distribution.

    The density of the Rayleigh distribution is given by

    .. math ::
        f(x) := \frac{x}{\sigma^2} e^{-x^2/2\sigma^2}

    with `x > 0`.

    Parameters
    ==========

    sigma : Real number, `\sigma > 0`

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, density, variance

    >>> sigma = Symbol('sigma', positive=True)

    >>> X = Rayleigh('x', sigma)

    >>> density(X)(z)
    E**(-z**2/(2*sigma**2))*z/sigma**2

    >>> E(X)
    sqrt(2)*sqrt(pi)*sigma/2

    >>> variance(X)
    -pi*sigma**2/2 + 2*sigma**2

    References
    ==========

    * https://en.wikipedia.org/wiki/Rayleigh_distribution
    * https://mathworld.wolfram.com/RayleighDistribution.html

    """
    return rv(name, RayleighDistribution, (sigma, ))

# ------------------------------------------------------------------------------
# StudentT distribution --------------------------------------------------------


class StudentTDistribution(SingleContinuousDistribution):
    _argnames = 'nu',

    def pdf(self, x):
        nu = self.nu
        return 1/(sqrt(nu)*beta_fn(Rational(1, 2), nu/2))*(1 + x**2/nu)**(-(nu + 1)/2)


def StudentT(name, nu):
    r"""
    Create a continuous random variable with a student's t distribution.

    The density of the student's t distribution is given by

    .. math::
        f(x) := \frac{\Gamma \left(\frac{\nu+1}{2} \right)}
                {\sqrt{\nu\pi}\Gamma \left(\frac{\nu}{2} \right)}
                \left(1+\frac{x^2}{\nu} \right)^{-\frac{\nu+1}{2}}

    Parameters
    ==========

    nu : Real number, `\nu > 0`, the degrees of freedom

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> nu = Symbol('nu', positive=True)

    >>> X = StudentT('x', nu)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
                nu   1
              - -- - -
                2    2
      /     2\
      |    z |
      |1 + --|
      \    nu/
    --------------------
      ____     /     nu\
    \/ nu *beta|1/2, --|
               \     2 /

    References
    ==========

    * https://en.wikipedia.org/wiki/Student_t-distribution
    * https://mathworld.wolfram.com/Studentst-Distribution.html

    """
    return rv(name, StudentTDistribution, (nu, ))

# ------------------------------------------------------------------------------
# Triangular distribution ------------------------------------------------------


class TriangularDistribution(SingleContinuousDistribution):
    _argnames = ('a', 'b', 'c')

    def pdf(self, x):
        a, b, c = self.a, self.b, self.c
        return Piecewise(
            (2*(x - a)/((b - a)*(c - a)), And(a <= x, x < c)),
            (2/(b - a), Eq(x, c)),
            (2*(b - x)/((b - a)*(b - c)), And(c < x, x <= b)),
            (0, True))


def Triangular(name, a, b, c):
    r"""
    Create a continuous random variable with a triangular distribution.

    The density of the triangular distribution is given by

    .. math::
        f(x) := \begin{cases}
                  0 & \mathrm{for\ } x < a, \\
                  \frac{2(x-a)}{(b-a)(c-a)} & \mathrm{for\ } a \le x < c, \\
                  \frac{2}{b-a} & \mathrm{for\ } x = c, \\
                  \frac{2(b-x)}{(b-a)(b-c)} & \mathrm{for\ } c < x \le b, \\
                  0 & \mathrm{for\ } b < x.
                \end{cases}

    Parameters
    ==========

    a : Real number, `a \in \left(-\infty, \infty\right)`
    b : Real number, `a < b`
    c : Real number, `a \leq c \leq b`

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> X = Triangular('x', a, b, c)

    >>> pprint(density(X)(z), use_unicode=False)
    /    -2*a + 2*z
    |-----------------  for And(a <= z, z < c)
    |(-a + b)*(-a + c)
    |
    |       2
    |     ------              for z = c
    <     -a + b
    |
    |   2*b - 2*z
    |----------------   for And(z <= b, c < z)
    |(-a + b)*(b - c)
    |
    \        0                otherwise

    References
    ==========

    * https://en.wikipedia.org/wiki/Triangular_distribution
    * https://mathworld.wolfram.com/TriangularDistribution.html

    """
    return rv(name, TriangularDistribution, (a, b, c))

# ------------------------------------------------------------------------------
# Uniform distribution ---------------------------------------------------------


class UniformDistribution(SingleContinuousDistribution):
    _argnames = ('left', 'right')

    def pdf(self, x):
        left, right = self.left, self.right
        return Piecewise(
            (1/(right - left), And(left <= x, x <= right)),
            (0, True))

    def compute_cdf(self, **kwargs):
        from ..functions import Min
        z = Dummy('z', real=True)
        result = SingleContinuousDistribution.compute_cdf(self, **kwargs)(z)
        reps = {
            Min(z, self.right): z,
            Min(z, self.left, self.right): self.left,
            Min(z, self.left): self.left}
        result = result.subs(reps)
        return Lambda(z, result)

    def expectation(self, expr, var, **kwargs):
        from ..functions import Max, Min
        kwargs['evaluate'] = True
        result = SingleContinuousDistribution.expectation(self, expr, var, **kwargs)
        result = result.subs({Max(self.left, self.right): self.right,
                              Min(self.left, self.right): self.left})
        return result

    def sample(self):
        return random.uniform(self.left, self.right)


def Uniform(name, left, right):
    r"""
    Create a continuous random variable with a uniform distribution.

    The density of the uniform distribution is given by

    .. math::
        f(x) := \begin{cases}
                  \frac{1}{b - a} & \text{for } x \in [a,b]  \\
                  0               & \text{otherwise}
                \end{cases}

    with `x \in [a,b]`.

    Parameters
    ==========

    a : Real number, `-\infty < a` the left boundary
    b : Real number, `a < b < \infty` the right boundary

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, cdf, density, variance

    >>> a = Symbol('a', negative=True)
    >>> b = Symbol('b', positive=True)

    >>> X = Uniform('x', a, b)

    >>> density(X)(z)
    Piecewise((1/(-a + b), (a <= z) & (z <= b)), (0, true))

    >>> cdf(X)(z)
    -a/(-a + b) + z/(-a + b)

    >>> simplify(E(X))
    a/2 + b/2

    >>> simplify(variance(X))
    a**2/12 - a*b/6 + b**2/12

    References
    ==========

    * https://en.wikipedia.org/wiki/Uniform_distribution_%28continuous%29
    * https://mathworld.wolfram.com/UniformDistribution.html

    """
    return rv(name, UniformDistribution, (left, right))

# ------------------------------------------------------------------------------
# UniformSum distribution ------------------------------------------------------


class UniformSumDistribution(SingleContinuousDistribution):
    _argnames = 'n',

    @property
    def n(self):
        return self.args[self._argnames.index('n')]

    @property
    def set(self):
        return Interval(0, self.n)

    def pdf(self, x):
        n = self.n
        k = Dummy('k')
        return 1/factorial(
            n - 1)*Sum((-1)**k*binomial(n, k)*(x - k)**(n - 1), (k, 0, floor(x)))


def UniformSum(name, n):
    r"""
    Create a continuous random variable with an Irwin-Hall distribution.

    The probability distribution function depends on a single parameter
    `n` which is an integer.

    The density of the Irwin-Hall distribution is given by

    .. math ::
        f(x) := \frac{1}{(n-1)!}\sum_{k=0}^{\lfloor x\rfloor}(-1)^k
                \binom{n}{k}(x-k)^{n-1}

    Parameters
    ==========

    n : A positive Integer, `n > 0`

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> X = UniformSum('x', n)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
    floor(z)
      ___
      \  `
       \         k        n - 1 /n\
        )    (-1) *(z - k)     *| |
       /                        \k/
      /__,
     k = 0
    -------------------------------
                (n - 1)!

    References
    ==========

    * https://en.wikipedia.org/wiki/Uniform_sum_distribution
    * https://mathworld.wolfram.com/UniformSumDistribution.html

    """
    return rv(name, UniformSumDistribution, (n, ))

# ------------------------------------------------------------------------------
# VonMises distribution --------------------------------------------------------


class VonMisesDistribution(SingleContinuousDistribution):
    _argnames = ('mu', 'k')

    set = Interval(0, 2*pi)

    @staticmethod
    def check(mu, k):
        _value_check(k > 0, 'k must be positive')

    def pdf(self, x):
        mu, k = self.mu, self.k
        return exp(k*cos(x-mu)) / (2*pi*besseli(0, k))


def VonMises(name, mu, k):
    r"""
    Create a Continuous Random Variable with a von Mises distribution.

    The density of the von Mises distribution is given by

    .. math::
        f(x) := \frac{e^{\kappa\cos(x-\mu)}}{2\pi I_0(\kappa)}

    with `x \in [0,2\pi]`.

    Parameters
    ==========

    mu : Real number, measure of location
    k : Real number, measure of concentration

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import density

    >>> mu = Symbol('mu')
    >>> k = Symbol('k', positive=True)

    >>> X = VonMises('x', mu, k)

    >>> D = density(X)(z)
    >>> pprint(D, use_unicode=False)
       k*cos(mu - z)
      E
    ------------------
    2*pi*besseli(0, k)

    References
    ==========

    * https://en.wikipedia.org/wiki/Von_Mises_distribution
    * https://mathworld.wolfram.com/vonMisesDistribution.html

    """
    return rv(name, VonMisesDistribution, (mu, k))

# ------------------------------------------------------------------------------
# Weibull distribution ---------------------------------------------------------


class WeibullDistribution(SingleContinuousDistribution):
    _argnames = ('alpha', 'beta')

    set = Interval(0, oo, False, True)

    @staticmethod
    def check(alpha, beta):
        _value_check(alpha > 0, 'Alpha must be positive')
        _value_check(beta > 0, 'Beta must be positive')

    def pdf(self, x):
        alpha, beta = self.alpha, self.beta
        return beta * (x/alpha)**(beta - 1) * exp(-(x/alpha)**beta) / alpha

    def sample(self):
        return random.weibullvariate(self.alpha, self.beta)


def Weibull(name, alpha, beta):
    r"""
    Create a continuous random variable with a Weibull distribution.

    The density of the Weibull distribution is given by

    .. math::
        f(x) := \begin{cases}
                  \frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}
                  e^{-(x/\lambda)^{k}} & x\geq0\\
                  0 & x<0
                \end{cases}

    Parameters
    ==========

    lambda : Real number, `\lambda > 0` a scale
    k : Real number, `k > 0` a shape

    Returns
    =======

    A RandomSymbol.

    Examples
    ========

    >>> from diofant.stats import E, density, variance

    >>> l = Symbol('lambda', positive=True)
    >>> k = Symbol('k', positive=True, real=True)

    >>> X = Weibull('x', l, k)

    >>> density(X)(z)
    E**(-(z/lambda)**k)*k*(z/lambda)**(k - 1)/lambda

    >>> simplify(E(X))
    lambda*gamma(1 + 1/k)

    >>> simplify(variance(X))
    lambda**2*(-gamma(1 + 1/k)**2 + gamma(1 + 2/k))

    References
    ==========

    * https://en.wikipedia.org/wiki/Weibull_distribution
    * https://mathworld.wolfram.com/WeibullDistribution.html

    """
    return rv(name, WeibullDistribution, (alpha, beta))

# ------------------------------------------------------------------------------
# Wigner semicircle distribution -----------------------------------------------


class WignerSemicircleDistribution(SingleContinuousDistribution):
    _argnames = 'R',

    @property
    def set(self):
        return Interval(-self.R, self.R)

    def pdf(self, x):
        R = self.R
        return 2/(pi*R**2)*sqrt(R**2 - x**2)


def WignerSemicircle(name, R):
    r"""
    Create a continuous random variable with a Wigner semicircle distribution.

    The density of the Wigner semicircle distribution is given by

    .. math::
        f(x) := \frac2{\pi R^2}\,\sqrt{R^2-x^2}

    with `x \in [-R,R]`.

    Parameters
    ==========

    R : Real number, `R > 0`, the radius

    Returns
    =======

    A `RandomSymbol`.

    Examples
    ========

    >>> from diofant.stats import E, density

    >>> R = Symbol('R', positive=True)

    >>> X = WignerSemicircle('x', R)

    >>> density(X)(z)
    2*sqrt(R**2 - z**2)/(pi*R**2)

    >>> E(X)
    0

    References
    ==========

    * https://en.wikipedia.org/wiki/Wigner_semicircle_distribution
    * https://mathworld.wolfram.com/WignersSemicircleLaw.html

    """
    return rv(name, WignerSemicircleDistribution, (R,))
