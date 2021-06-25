import pytest

from diofant import (And, Basic, Dict, DiracDelta, Eq, Interval, Rational, Sum,
                     Symbol, Tuple, cos, integrate, oo, sin, symbols)
from diofant.abc import x, z
from diofant.stats import (Die, E, Exponential, Normal, P, density, dependent,
                           given, independent, pspace, random_symbols, sample,
                           variance, where)
from diofant.stats.rv import (Density, NamedArgsMixin, ProductDomain,
                              ProductPSpace, RandomSymbol, rs_swap,
                              sample_iter)


__all__ = ()


def test_where():
    X, Y = Die('X'), Die('Y')
    Z = Normal('Z', 0, 1)

    assert where(Z**2 <= 1).set == Interval(-1, 1)
    assert where(
        Z**2 <= 1).as_boolean() == Interval(-1, 1).as_relational(Z.symbol)
    assert where(And(X > Y, Y > 4)).as_boolean() == And(
        Eq(X.symbol, 6), Eq(Y.symbol, 5))

    assert len(where(X < 3).set) == 2
    assert 1 in where(X < 3).set

    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    assert where(And(X**2 <= 1, X >= 0)).set == Interval(0, 1)
    XX = given(X, And(X**2 <= 1, X >= 0))
    assert XX.pspace.domain.set == Interval(0, 1)
    assert XX.pspace.domain.as_boolean() == \
        And(0 <= X.symbol, X.symbol**2 <= 1, -oo < X.symbol, X.symbol < oo)

    with pytest.raises(TypeError):
        XX = given(X, X + 3)


def test_random_symbols():
    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)

    assert set(random_symbols(2*X + 1)) == {X}
    assert set(random_symbols(2*X + Y)) == {X, Y}
    assert set(random_symbols(2*X + Y.symbol)) == {X}
    assert set(random_symbols(2)) == set()

    assert X.is_commutative


def test_pspace():
    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)

    pytest.raises(ValueError, lambda: pspace(5 + 3))
    pytest.raises(ValueError, lambda: pspace(x < 1))
    assert pspace(X) == X.pspace
    assert pspace(2*X + 1) == X.pspace
    assert pspace(2*X + Y) == ProductPSpace(Y.pspace, X.pspace)


def test_rs_swap():
    X = Normal('x', 0, 1)
    Y = Exponential('y', 1)

    XX = Normal('x', 0, 2)
    YY = Normal('y', 0, 3)

    expr = 2*X + Y
    assert expr.subs(rs_swap((X, Y), (YY, XX))) == 2*XX + YY


def test_RandomSymbol():
    pytest.raises(TypeError, lambda: RandomSymbol(0, 1))
    pytest.raises(TypeError, lambda: RandomSymbol(0, Symbol('x')))

    X = Normal('x', 0, 1)
    Y = Normal('x', 0, 2)
    assert X.symbol == Y.symbol
    assert X != Y

    assert X.name == X.symbol.name

    X = Normal('lambda', 0, 1)  # make sure we can use protected terms
    X = Normal('Lambda', 0, 1)  # make sure we can use Diofant terms

    pytest.raises(TypeError, lambda: Normal(1, 0, 1))


def test_RandomSymbol_diff():
    X = Normal('x', 0, 1)
    assert (2*X).diff(X)


def test_overlap():
    X = Normal('x', 0, 1)
    Y = Normal('x', 0, 2)

    pytest.raises(ValueError, lambda: P(X > Y))


def test_ProductPSpace():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    px = X.pspace
    py = Y.pspace
    assert pspace(X + Y) == ProductPSpace(px, py)
    assert pspace(X + Y) == ProductPSpace(py, px)

    X = Die('X', 2)
    Y = Die('Y', 2)

    assert (pspace(X + Y).density ==
            Dict((frozenset({('X', 1), ('Y', 2)}), Rational(1, 4)),
                 (frozenset({('X', 1), ('Y', 1)}), Rational(1, 4)),
                 (frozenset({('X', 2), ('Y', 1)}), Rational(1, 4)),
                 (frozenset({('X', 2), ('Y', 2)}), Rational(1, 4))))
    d = pspace(X + Y).domain
    assert ((X.symbol, 1), (Y.symbol, 2)) in d
    assert ((X.symbol, 0), (Y.symbol, 2)) not in d

    Z = Die('Z', 2)
    d1 = pspace(X + Y).domain
    assert ProductDomain(d1, Z.pspace.domain) == pspace(X + Y + Z).domain


def test_E():
    assert E(5) == 5


def test_Sample():
    X = Die('X', 6)
    Y = Normal('Y', 0, 1)

    assert sample(X) in [1, 2, 3, 4, 5, 6]
    assert sample(X + Y).is_Float

    assert P(X + Y > 0, Y < 0, numsamples=10).is_number
    assert P(X > 10, numsamples=10).is_number
    assert E(X + Y, numsamples=10).is_number
    assert variance(X + Y, numsamples=10).is_number

    pytest.raises(TypeError, lambda: P(Y > z, numsamples=5))

    assert P(sin(Y) <= 1, numsamples=10, modules=['math']) == 1
    assert P(sin(Y) <= 1, cos(Y) < 1, numsamples=10, modules=['math']) == 1

    assert all(i in range(1, 7) for i in density(X, numsamples=10))
    assert all(i in range(4, 7) for i in density(X, X > 3, numsamples=10))

    # Make sure this doesn't raise an error
    Y = Normal('Y', 0, 1)
    E(Sum(1/z**Y, (z, 1, oo)), Y > 2, numsamples=3, modules='mpmath')


def test_given():
    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    A = given(X, True)
    B = given(X, Y > 2)

    assert X == A == B


def test_dependence():
    X, Y = Die('X'), Die('Y')
    assert independent(X, 2*Y)
    assert not dependent(X, 2*Y)

    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    assert independent(X, Y)
    assert dependent(X, 2*X)

    # Create a dependency
    XX, YY = given(Tuple(X, Y), Eq(X + Y, 3))
    assert dependent(XX, YY)


@pytest.mark.xfail
def test_dependent_finite():
    X, Y = Die('X'), Die('Y')
    # Dependence testing requires symbolic conditions which currently break
    # finite random variables
    assert dependent(X, Y + X)

    XX, YY = given(Tuple(X, Y), X + Y > 5)  # Create a dependency
    assert dependent(XX, YY)


def test_normality():
    X, Y = Normal('X', 0, 1), Normal('Y', 0, 1)
    x, z = symbols('x, z', real=True)
    dens = density(X - Y, Eq(X + Y, z))

    assert integrate(dens(x), (x, -oo, oo)) == 1


def test_Density():
    X = Die('X', 6)
    d = Density(X)
    assert d.doit() == density(X)
    r = Rational(1, 6)
    assert density(2*X).dict == {2: r, 4: r, 6: r, 8: r, 10: r, 12: r}

    X = Normal('X', 0, 1)
    Y = Normal('Y', 0, 1)
    assert str(density(X - Y, Y)(z)) == 'sqrt(2)*E**(-(Y + z)**2/2)/(2*sqrt(pi))'


def test_NamedArgsMixin():
    class Foo(Basic, NamedArgsMixin):
        _argnames = 'foo', 'bar'

    a = Foo(1, 2)

    assert a.foo == 1
    assert a.bar == 2

    pytest.raises(AttributeError, lambda: a.baz)

    class Bar(Basic, NamedArgsMixin):
        pass

    pytest.raises(AttributeError, lambda: Bar(1, 2).foo)


def test_density_constant():
    assert density(3)(2) == 0
    assert density(3)(3) == DiracDelta(0)


def test_real():
    x = Normal('x', 0, 1)
    assert x.is_extended_real


def test_exponential():
    # issue sympy/sympy#10052:
    X = Exponential('X', 3)
    assert P(X < oo) == 1
    assert P(X > oo) == 0
    assert P(X < 2, X > oo) == 0
    assert P(X < oo, X > oo) == 0
    assert P(X < oo, X > 2) == 1
    assert P(X < 3, X == 2) == 0
    pytest.raises(ValueError, lambda: P(1))
    pytest.raises(ValueError, lambda: P(X < 1, 2))


def test_sample_iter():
    X = Normal('X', 0, 1)
    expr = X*X + 3
    assert all(_ > 3 for _ in sample_iter(expr, numsamples=3))
    assert all(_ > 5 for _ in sample_iter(expr, numsamples=3, condition=X > 2))
    pytest.raises(ValueError, lambda: tuple(sample_iter(expr, numsamples=3,
                                                        condition=X)))


def test_sympyissue_8129():
    X = Exponential('X', 4)
    assert P(X >= X) == 1
    assert P(X > X) == 0
    assert P(X > X+1) == 0
