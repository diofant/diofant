from diofant import Eq, Symbol
from diofant.stats import Beta, Poisson
from diofant.stats.drv_types import PoissonDistribution
from diofant.stats.rv import ProductPSpace, density, pspace


__all__ = ()


def test_density():
    x = Symbol('x')
    l = Symbol('l', positive=True)
    rate = Beta(l, 2, 3)
    X = Poisson(x, rate)
    assert isinstance(pspace(X), ProductPSpace)
    assert density(X, Eq(rate, rate.symbol)) == PoissonDistribution(l)
