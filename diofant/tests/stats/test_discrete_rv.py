from diofant import Rational, Sum, exp, factorial
from diofant.abc import x, z
from diofant.stats import E, density, variance
from diofant.stats.drv_types import (Geometric, GeometricDistribution, Poisson,
                                     PoissonDistribution)


__all__ = ()


def test_PoissonDistribution():
    l = 3
    p = PoissonDistribution(l)
    assert abs(p.cdf(10).evalf() - 1) < .001
    assert p.expectation(x, x) == l
    assert p.expectation(x**2, x) - p.expectation(x, x)**2 == l


def test_Poisson():
    l = 3
    x = Poisson('x', l)
    assert E(x) == l
    assert variance(x) == l
    assert density(x) == PoissonDistribution(l)
    assert isinstance(E(x, evaluate=False), Sum)
    assert isinstance(E(2*x, evaluate=False), Sum)
    assert density(x)(z) == 3**z/(exp(3)*factorial(z))


def test_GeometricDistribution():
    p = Rational(1, 5)
    d = GeometricDistribution(p)
    assert d.expectation(x, x) == 1/p
    assert d.expectation(x**2, x) - d.expectation(x, x)**2 == (1-p)/p**2
    assert abs(d.cdf(20000).evalf() - 1) < .001

    X = Geometric("x", p)
    assert density(X)(z) == Rational(4, 5)**(z - 1)/5
