from diofant.stats import Poisson, Beta
from diofant.stats.rv import pspace, ProductPSpace, density
from diofant.stats.drv_types import PoissonDistribution
from diofant import Symbol, Eq


def test_density():
    x = Symbol('x')
    l = Symbol('l', positive=True)
    rate = Beta(l, 2, 3)
    X = Poisson(x, rate)
    assert isinstance(pspace(X), ProductPSpace)
    assert density(X, Eq(rate, rate.symbol)) == PoissonDistribution(l)
