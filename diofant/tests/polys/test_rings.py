"""Test sparse polynomials."""

import functools
import math
import operator

import pytest

from diofant import (EX, FF, QQ, RR, ZZ, CoercionFailed, ExactQuotientFailed,
                     GeneratorsError, GeneratorsNeeded,
                     PolynomialDivisionFailed, PolynomialRing, Rational,
                     Symbol, field, grlex, lex, pi, ring, sin, sqrt, symbols)
from diofant.abc import t, x, y, z
from diofant.polys.rings import PolyElement
from diofant.polys.specialpolys import f_polys


__all__ = ()


def test_PolynomialRing___init__():
    assert len(PolynomialRing(ZZ, 'x,y,z').gens) == 3
    assert len(ZZ.inject(x).gens) == 1
    assert len(ZZ.inject('x', 'y', 'z').gens) == 3
    assert len(ZZ.inject(x, y, z).gens) == 3

    pytest.raises(GeneratorsNeeded, lambda: ZZ.inject())
    pytest.raises(GeneratorsError, lambda: ZZ.inject(0))

    assert ZZ.inject(t).poly_ring('x').domain == ZZ.inject(t)
    assert PolynomialRing('ZZ[t]', 'x').domain == ZZ.inject(t)

    pytest.raises(GeneratorsError, lambda: ZZ.inject('x').poly_ring('x'))

    _lex = Symbol('lex')

    assert PolynomialRing(ZZ, 'x').order == lex
    assert PolynomialRing(ZZ, 'x', _lex).order == lex
    assert PolynomialRing(ZZ, 'x', 'lex').order == lex

    R1 = ZZ.inject('x', 'y')
    R2 = ZZ.inject('x', 'y')
    R3 = ZZ.inject('x', 'y', 'z')

    assert R1.x == R1.gens[0]
    assert R1.y == R1.gens[1]
    assert R1.x == R2.x
    assert R1.y == R2.y
    assert R1.x != R3.x
    assert R1.y != R3.y

    R4 = ZZ.inject('gens')

    assert type(R4.gens) is tuple

    pytest.raises(GeneratorsError, lambda: PolynomialRing(ZZ, {1: 2}))
    pytest.raises(GeneratorsError, lambda: PolynomialRing(ZZ, ['x', ['y']]))

    R, _ = ring(x, ZZ)

    assert R.ngens == 1 and R.domain == ZZ


def test_PolynomialRing___hash__():
    R, x, y, z = ring('x y z', QQ)

    assert hash(R)


def test_PolynomialRing___eq__():
    assert ring('x y z', QQ)[0] == ring('x y z', QQ)[0]
    assert ring('x y z', QQ)[0] is ring('x y z', QQ)[0]

    assert ring('x y z', QQ)[0] != ring('x y z', ZZ)[0]
    assert ring('x y z', QQ)[0] is not ring('x y z', ZZ)[0]

    assert ring('x y z', ZZ)[0] != ring('x y z', QQ)[0]
    assert ring('x y z', ZZ)[0] is not ring('x y z', QQ)[0]

    assert ring('x y z', QQ)[0] != ring('x y', QQ)[0]
    assert ring('x y z', QQ)[0] is not ring('x y', QQ)[0]

    assert ring('x y', QQ)[0] != ring('x y z', QQ)[0]
    assert ring('x y', QQ)[0] is not ring('x y z', QQ)[0]


def test_PolynomialRing__call__():
    R, x = ring('x', ZZ)

    assert R({2: 1, 1: 0, 0: -1}) == x**2 - 1
    assert R([-1, 0, 1]) == x**2 - 1
    assert R([((2,), 1), ((0,), -1)]) == x**2 - 1

    D, t = ring('t', ZZ)
    R, x = ring('x', D)

    assert R(t) == t

    R, x, y, z = ring('x y z', QQ)

    assert R(7) == R({(0, 0, 0): 7}) == 7
    assert R(7*x*y*z) == 7*x*y*z

    f = x**2 + 2*x*y + 3*x + 4*z**2 + 5*z + 6

    assert R({(2, 0, 0): 1, (1, 1, 0): 2, (1, 0, 0): 3, (0, 0, 2): 4, (0, 0, 1): 5, (0, 0, 0): 6}) == f
    assert R([((2, 0, 0), 1), ((1, 1, 0), 2), ((1, 0, 0), 3), ((0, 0, 2), 4), ((0, 0, 1), 5), ((0, 0, 0), 6)]) == f


def test_PolynomialRing_domain_new():
    R1 = ZZ.frac_field(t).poly_ring(y)
    R2 = R1.poly_ring(x, z)

    assert R2.domain_new(R1.ring.domain.one) == 1


def test_PolynomialRing_drop():
    R, x, y, z = ring('x y z', ZZ)

    assert R.drop(x) == ZZ.inject('y', 'z')
    assert R.drop(y) == ZZ.inject('x', 'z')
    assert R.drop(z) == ZZ.inject('x', 'y')

    assert R.drop(0) == ZZ.inject('y', 'z')
    assert R.drop(0).drop(0) == ZZ.inject('z')
    assert R.drop(0).drop(0).drop(0) == ZZ

    assert R.drop(1) == ZZ.inject('x', 'z')

    assert R.drop(2) == ZZ.inject('x', 'y')
    assert R.drop(2).drop(1) == ZZ.inject('x')
    assert R.drop(2).drop(1).drop(0) == ZZ

    pytest.raises(ValueError, lambda: R.drop(3))
    pytest.raises(ValueError, lambda: R.drop(x).drop(y))


def test_PolynomialRing_index():
    R, x, y, z = ring('x y z', ZZ)

    assert R.index(0) == 0
    assert R.index(-1) == 2

    pytest.raises(ValueError, lambda: R.index(100))

    assert R.index(x) == 0
    assert R.index(y) == 1

    pytest.raises(ValueError, lambda: R.index(x + y))

    assert R.index('x') == 0
    assert R.index('z') == 2

    pytest.raises(ValueError, lambda: R.index('t'))

    assert R.index(Symbol('x')) == 0
    assert R.index(Symbol('z')) == 2

    pytest.raises(ValueError, lambda: R.index(Symbol('t')))


def test_PolynomialRing_is_():
    R = QQ.inject('x')

    assert R.is_univariate is True
    assert R.is_multivariate is False

    R = QQ.inject('x', 'y', 'z')

    assert R.is_univariate is False
    assert R.is_multivariate is True


def test_PolynomialRing_add():
    R, x = ring('x', ZZ)

    F = [x**2 + 2*i + 3 for i in range(4)]

    assert functools.reduce(operator.add, F) == 4*x**2 + 24


def test_PolynomialRing_mul():
    R, x = ring('x', ZZ)

    F = [x**2 + 2*i + 3 for i in range(4)]

    assert functools.reduce(operator.mul, F) == (x**8 + 24*x**6 +
                                                 206*x**4 + 744*x**2 + 945)


def test_PolynomialRing_to_ground():
    R, x = ring('x', ZZ)

    pytest.raises(ValueError, lambda: R.to_ground())

    R2, x, y = ring('x y', ZZ)

    assert R2.eject(x) == ZZ.inject('x').poly_ring('y')
    assert R2.eject(x, y) == R2


def test_PolyElement___hash__():
    R, x, y, z = ring('x y z', QQ)

    assert hash(x*y*z)


def test_PolyElement___eq__():
    R, x, y = ring('x y', ZZ)

    assert ((x*y + 5*x*y) == 6) is False
    assert ((x*y + 5*x*y) == 6*x*y) is True
    assert (6 == (x*y + 5*x*y)) is False
    assert (6*x*y == (x*y + 5*x*y)) is True

    assert ((x*y - x*y) == 0) is True
    assert (0 == (x*y - x*y)) is True

    assert ((x*y - x*y) == 1) is False
    assert (1 == (x*y - x*y)) is False

    assert ((x*y - x*y) == 1) is False
    assert (1 == (x*y - x*y)) is False

    assert ((x*y + 5*x*y) != 6) is True
    assert ((x*y + 5*x*y) != 6*x*y) is False
    assert (6 != (x*y + 5*x*y)) is True
    assert (6*x*y != (x*y + 5*x*y)) is False

    assert ((x*y - x*y) != 0) is False
    assert (0 != (x*y - x*y)) is False

    assert ((x*y - x*y) != 1) is True
    assert (1 != (x*y - x*y)) is True

    Rt, t = ring('t', ZZ)
    R, x, y = ring('x y', Rt)

    assert (t**3*x//x == t**3) is True
    assert (t**3*x//x == t**4) is False


def test_PolyElement_copy():
    R, x, y, z = ring('x y z', ZZ)

    f = x*y + 3*z
    g = f.copy()

    assert f == g

    g[(1, 1, 1)] = 7

    assert f != g


def test_PolyElement_inject():
    D, y, z = ring('y z', ZZ)
    R, x = ring('x', D)

    p = x*y*z + 1

    R2 = ZZ.inject('x', 'y', 'z')
    R3 = ZZ.inject('y', 'z', 'x')

    p2 = p.inject()

    assert p2.ring == R2
    assert p2 == R2.x*R2.y*R2.z + 1

    p3 = p.inject(front=True)

    assert p3.ring == R3
    assert p3 == R3.x*R3.y*R3.z + 1


def test_PolyElement_set_domain():
    R, x, y = ring('x y', ZZ)

    f = x + y

    assert f.set_domain(ZZ) is f

    g = f.set_domain(QQ)
    R2 = g.ring

    assert g is not f
    assert R2.to_expr(g) == R.to_expr(f)
    assert R2.domain is QQ


def test_PolyElement_items():
    R, x, y = ring('x y', ZZ)

    f = x*y**7 + 2*x**2*y**3

    assert set(f.keys()) == {(2, 3), (1, 7)}
    assert set(f.values()) == {1, 2}
    assert set(f.items()) == {((2, 3), 2), ((1, 7), 1)}

    R, x, y, z = ring('x y z', ZZ)

    f = x*y + 3*z

    assert list(f.items()) == [((1, 1, 0), 1), ((0, 0, 1), 3)]

    R, x, y, z = ring('x y z', QQ)

    f = x**2/3 + y**3/4 + z**4/5

    assert set(f.keys()) == {(2, 0, 0), (0, 3, 0), (0, 0, 4)}
    assert set(f.values()) == {QQ(1, 3), QQ(1, 4), QQ(1, 5)}
    assert set(f.items()) == {((2, 0, 0), QQ(1, 3)), ((0, 3, 0), QQ(1, 4)),
                              ((0, 0, 4), QQ(1, 5))}


def test_PolynomialRing_from_list():
    R, x = ring('x', ZZ)

    assert R.from_list([]) == R(0)

    f = [ZZ(8), ZZ(0), ZZ(0), ZZ(0), ZZ(0), ZZ(2), ZZ(0), ZZ(0), ZZ(3)]
    g = 3*x**8 + 2*x**5 + 8

    assert R.from_list(f) == g

    f = [ZZ(7), ZZ(0), ZZ(5), ZZ(0), ZZ(1)]
    g = x**4 + 5*x**2 + 7

    assert R.from_list(f) == g

    R, x, y = ring('x y', ZZ)
    R1, z = ring('z', R)

    f = [R(8), R(0), R(0), R(2), R(0), R(3)]
    g = 3*z**5 + 2*z**3 + 8

    assert R1.from_list(f) == g


def test_PolyElement_as_expr():
    R, x, y, z = ring('x y z', ZZ)

    f = 3*x**2*y - x*y*z + 7*z**3 + 1

    x, y, z = R.symbols
    g = 3*x**2*y - x*y*z + 7*z**3 + 1

    assert f != g
    assert R.to_expr(f) == g

    x, y, z = symbols('x y z')
    g = 3*x**2*y - x*y*z + 7*z**3 + 1

    assert f != g


def test_PolyElement_from_expr():
    x, y, z = symbols('x y z')
    R, X, Y, Z = ring((x, y, z), ZZ)

    f = R.convert(1)

    assert f == 1 and isinstance(f, R.dtype)

    f = R.convert(x)

    assert f == X and isinstance(f, R.dtype)

    f = R.convert(x*y*z)

    assert f == X*Y*Z and isinstance(f, R.dtype)

    f = R.convert(x*y*z + x*y + x)

    assert f == X*Y*Z + X*Y + X and isinstance(f, R.dtype)

    f = R.convert(x**3*y*z + x**2*y**7 + 1)

    assert f == X**3*Y*Z + X**2*Y**7 + 1 and isinstance(f, R.dtype)

    pytest.raises(CoercionFailed, lambda: R.convert(1/x))
    pytest.raises(CoercionFailed, lambda: R.convert(2**x))
    pytest.raises(CoercionFailed, lambda: R.convert(7*x + sqrt(2)))

    R, X, Y = ring((2**x, y), ZZ)

    f = R.convert(2**(2*x) + 1)

    assert f == X**2 + 1


def test_PolyElement_degree():
    R, x = ring('x', ZZ)

    assert R(0).degree() == -math.inf
    assert R(1).degree() == 0
    assert x.degree() == 1
    assert (x**4 + 1).degree() == 4
    assert (x**3 + 2*x**2 + 3).degree() == 3
    assert (x**3 + x**2 + 2*x).degree() == 3

    R, x, y = ring('x y', ZZ)

    assert R(0).degree() == -math.inf
    assert R(1).degree() == 0
    assert (2*x + 1).degree() == 1
    assert (2*x + y**2 + 2*y + 3).degree() == 1
    assert (2*x + y**2 + 2*y + 3).degree(y) == 2

    R, x, y, z = ring('x y z', ZZ)

    assert R(0).degree() == -math.inf
    assert R(1).degree() == 0
    assert (x + 1).degree() == 1
    assert (2*y**3 + z).degree() == 0
    assert (x*y**3 + z).degree() == 1
    assert (x**5*y**3 + z).degree() == 5

    assert (x**5*y**3 + z).degree(0) == 5
    assert (x**5*y**3 + z).degree(-3) == 5

    pytest.raises(ValueError, lambda: (x**5*y**3 + z).degree(100))

    assert R(0).degree(x) == -math.inf
    assert R(1).degree(x) == 0
    assert (x + 1).degree(x) == 1
    assert (2*y**3 + z).degree(x) == 0
    assert (x*y**3 + z).degree(x) == 1
    assert (7*x**5*y**3 + z).degree(x) == 5

    assert R(0).degree(y) == -math.inf
    assert R(1).degree(y) == 0
    assert (x + 1).degree(y) == 0
    assert (2*y**3 + z).degree(y) == 3
    assert (x*y**3 + z).degree(y) == 3
    assert (7*x**5*y**3 + z).degree(y) == 3

    assert R(0).degree(z) == -math.inf
    assert R(1).degree(z) == 0
    assert (x + 1).degree(z) == 0
    assert (2*y**3 + z).degree(z) == 1
    assert (x*y**3 + z).degree(z) == 1
    assert (7*x**5*y**3 + z).degree(z) == 1

    f = f_polys()[4]

    assert f.degree(x) == 9
    assert f.degree(y) == 12
    assert f.degree(z) == 8

    R, x, y, z, t = ring('x y z t', ZZ)

    f = f_polys()[6]

    assert f.degree(x) == 4
    assert f.degree(y) == 4
    assert f.degree(z) == 6
    assert f.degree(t) == 3


def test_PolyElement_tail_degree():
    R, x, y, z = ring('x y z', ZZ)

    assert R(0).tail_degree() == -math.inf
    assert R(1).tail_degree() == 0
    assert (x + 1).tail_degree() == 0
    assert (2*y**3 + x**3*z).tail_degree() == 0
    assert (x*y**3 + x**3*z).tail_degree() == 1
    assert (x**5*y**3 + x**3*z).tail_degree() == 3
    assert (x**2*y + x**3*z**2).tail_degree() == 2

    assert R(0).tail_degree(x) == -math.inf
    assert R(1).tail_degree(x) == 0
    assert (x + 1).tail_degree(x) == 0
    assert (2*y**3 + x**3*z).tail_degree(x) == 0
    assert (x*y**3 + x**3*z).tail_degree(x) == 1
    assert (7*x**5*y**3 + x**3*z).tail_degree(x) == 3
    assert (x**2*y + x**3*z**2).tail_degree(x) == 2

    assert R(0).tail_degree(y) == -math.inf
    assert R(1).tail_degree(y) == 0
    assert (x + 1).tail_degree(y) == 0
    assert (2*y**3 + x**3*z).tail_degree(y) == 0
    assert (x*y**3 + x**3*z).tail_degree(y) == 0
    assert (7*x**5*y**3 + x**3*z).tail_degree(y) == 0
    assert (x**2*y + x**3*z**2).tail_degree(y) == 0

    assert R(0).tail_degree(z) == -math.inf
    assert R(1).tail_degree(z) == 0
    assert (x + 1).tail_degree(z) == 0
    assert (2*y**3 + x**3*z).tail_degree(z) == 0
    assert (x*y**3 + x**3*z).tail_degree(z) == 0
    assert (7*x**5*y**3 + x**3*z).tail_degree(z) == 0
    assert (x**2*y + x**3*z**2).tail_degree(z) == 0


def test_PolyElement_degree_list():
    R, x, y = ring('x y', ZZ)

    assert [(x + y**2 + 2*y + 3).degree(_) for _ in (x, y)] == [1, 2]

    R, x, y, z = ring('x y z', ZZ)

    assert [R(0).degree(_) for _ in (x, y, z)] == [-math.inf]*3
    assert [R(1).degree(_) for _ in (x, y, z)] == [0, 0, 0]
    assert [(x**2*y + x**3*z**2).degree(_) for _ in (x, y, z)] == [3, 1, 2]

    assert [f_polys()[0].degree(_) for _ in (x, y, z)] == [2, 2, 2]
    assert [f_polys()[1].degree(_) for _ in (x, y, z)] == [3, 3, 3]
    assert [f_polys()[2].degree(_) for _ in (x, y, z)] == [5, 3, 3]
    assert [f_polys()[3].degree(_) for _ in (x, y, z)] == [5, 4, 7]
    assert [f_polys()[4].degree(_) for _ in (x, y, z)] == [9, 12, 8]
    assert [f_polys()[5].degree(_) for _ in (x, y, z)] == [3, 3, 3]

    R, x, y, z, t = ring('x y z t', ZZ)

    assert [R(0).degree(_) for _ in (x, y, z, t)] == [-math.inf]*4
    assert [R(1).degree(_) for _ in (x, y, z, t)] == [0, 0, 0, 0]

    assert [f_polys()[6].degree(_) for _ in (x, y, z, t)] == [4, 4, 6, 3]


def test_PolyEelemet_total_degree():
    R, x, y, z = ring('x y z', ZZ)

    assert R(0).total_degree() == -math.inf
    assert (x**2*y + x**3*z**2 + 1).total_degree() == 5
    assert (x**2 + z**3).total_degree() == 3
    assert (x*y*z + z**4).total_degree() == 4
    assert (x**3 + x + 1).total_degree() == 3


def test_PolyElement_coeff():
    R, x = ring('x', ZZ)

    assert R(0)[1] == 0
    assert R(1)[1] == 1
    assert (x + 2)[1] == 2
    assert (3*x**2 + 1)[1] == 1
    assert (2*x**3 + 3*x**2 + 4*x + 5)[1] == 5
    assert (x**2 + 2*x + 3)[1] == 3

    R, x, y = ring('x y', ZZ)

    f = R(0)

    assert f[1] == 0
    assert f.eject(y)[1] == 0

    f = 2*x*y**2 + 3*x*y + 4*x + 5

    assert f[1] == 5
    assert f.eject(y)[1] == 5

    R, x, y, z = ring('x y z', ZZ)

    f = R(0)

    assert f[1] == 0
    assert f.eject(y, z)[1] == 0

    f = 2*x*y + 3*x*z + 4*x + 5

    assert f[1] == 5
    assert f.eject(y, z)[1] == 5

    f = y + 2*z + 3

    assert f[1] == 3

    f = 3*x**2*y - x*y*z + 7*z**3 + 23

    assert f[1] == 23

    pytest.raises(ValueError, lambda: f[3])

    assert f[x] == 0
    assert f[y] == 0
    assert f[z] == 0

    assert f[x**2*y] == 3
    assert f[x*y*z] == -1
    assert f[(1, 1, 1)] == -1
    assert f[z**3] == 7
    assert f[(0, 0, 3)] == 7

    pytest.raises(ValueError, lambda: f[3*x**2*y])
    pytest.raises(ValueError, lambda: f[-x*y*z])
    pytest.raises(ValueError, lambda: f[7*z**3])
    pytest.raises(ValueError, lambda: f[x + y])

    f = 2*x + 3*x*y + 4*z + 5

    assert f[1] == R.domain(5)


def test_PolyElement_LC():
    R, x = ring('x', ZZ)

    assert R(0).LC == 0
    assert R(1).LC == 1
    assert (2*x**3 + 3*x**2 + 4*x + 5).LC == 2
    assert (3*x**2 + 1).LC == 3
    assert (x + 2).LC == 1
    assert (x**2 + 2*x + 3).LC == 1

    R, x, y = ring('x y', ZZ)

    assert R(0).LC == 0
    assert R(0).eject(-1).LC == 0
    assert (2*x*y**2 + 3*x*y + 4*x + 5).LC == 2

    R1 = R.eject(-1).domain

    f = x**2*y**2 + x**2*y - 1

    assert f.eject(-1).LC == R1.y**2 + R1.y

    f = 2*x*y**2 + 3*x*y + 4*x + 5

    assert f.eject(-1).LC == 2*R1.y**2 + 3*R1.y + 4

    R, x, y = ring('x y', QQ)

    assert R(0).LC == QQ(0)
    assert (x/2).LC == QQ(1, 2)
    assert (x*y/4 + x/2).LC == QQ(1, 4)

    R, x, y, z = ring('x y z', ZZ)

    R2 = R.eject(-1).domain

    assert R(0).LC == 0
    assert (y + 2*z + 3).LC == 1
    assert (2*x*y + 3*x*z + 4*x + 5).LC == 2

    f = x**2*y**2 + x**2*y - 1

    assert f.eject(-1).LC == 1

    f = x*y*z - y**2*z**2

    assert f.eject(-1).LC == R2.z

    R12 = R.drop(x)

    assert R(0).eject(y, z).LC == 0

    f = 2*x*y + 3*x*z + 4*x + 5

    assert f.eject(y, z).LC == 2*R12.y + 3*R12.z + 4


def test_PolyElement_LM():
    R, x, y = ring('x y', ZZ)

    f = x**2*y**2 + x**2*y - 1

    assert f.eject(-1).LM == (2,)

    R, x, y = ring('x y', QQ)

    assert R(0).LM == (0, 0)
    assert (x/2).LM == (1, 0)
    assert (x*y/4 + x/2).LM == (1, 1)
    assert (3*x*y + y**2).LM == (1, 1)

    R, x, y, z = ring('x y z', ZZ)

    f = x**2*y**2 + x**2*y - 1

    assert f.eject(-1).LM == (2, 2)

    f = x*y*z - y**2*z**2

    assert f.eject(-1).LM == (1, 1)


def test_PolyElement_LT():
    R, x, y = ring('x y', QQ)

    assert R(0).LT == ((0, 0), QQ(0))
    assert (x/2).LT == ((1, 0), QQ(1, 2))
    assert (x*y/4 + x/2).LT == ((1, 1), QQ(1, 4))


def test_PolyElement_leading_term():
    R, x, y = ring('x y', QQ)

    assert R(0).leading_term() == 0
    assert (x/2).leading_term() == x/2
    assert (x*y/4 + x/2).leading_term() == x*y/4


def test_PolyElement_all_coeffs():
    R, x = ring('x', ZZ)

    assert R.zero.all_coeffs() == [0]
    assert (3*x**2 + 2*x + 1).all_coeffs() == [1, 2, 3]
    assert (7*x**4 + 2*x + 1).all_coeffs() == [1, 2, 0, 0, 7]


def test_PolyElement__abs__():
    R, x = ring('x', ZZ)

    assert abs(R(0)) == 0
    assert abs(x**2 - 1) == x**2 + 1
    assert abs(R(1)) == 1
    assert abs(R(-7)) == 7
    assert abs(-x**2 + 2*x + 3) == x**2 + 2*x + 3
    assert abs(R(-1)) == 1

    R, x = ring('x', QQ)

    assert abs(R(0)) == 0
    assert abs(R(QQ(+1, 2))) == QQ(1, 2)
    assert abs(R(QQ(-7, 3))) == QQ(7, 3)
    assert abs(-x**2/7 + 2*x/7 + QQ(3, 7)) == x**2/7 + 2*x/7 + QQ(3, 7)
    assert abs(R(QQ(-1, 2))) == QQ(1, 2)

    R, x, y = ring('x y', ZZ)

    assert abs(x**2*y - x) == x**2*y + x

    R, x, y, z = ring('x y z', ZZ)

    assert abs(R(0)) == 0
    assert abs(R(1)) == 1
    assert abs(R(-7)) == 7

    R, x, y, z = ring('x y z', QQ)

    assert abs(R(0)) == 0
    assert abs(R(QQ(1, 2))) == QQ(1, 2)
    assert abs(R(QQ(-7, 9))) == QQ(7, 9)


def test_PolyElement__neg__():
    R, x = ring('x', ZZ)

    assert -R(0) == 0
    assert -(x**2 - 1) == 1 - x**2
    assert -R(1) == -1
    assert -R(-7) == 7
    assert -(-x**2 + 2*x + 3) == x**2 - 2*x - 3
    assert -R(-1) == 1

    R, x = ring('x', QQ)

    assert -R(0) == 0
    assert -R(QQ(1, 2)) == QQ(-1, 2)
    assert -R(QQ(-7, 9)) == QQ(7, 9)
    assert -(-x**2/7 + 2*x/7 + QQ(3, 7)) == x**2/7 - 2*x/7 - QQ(3, 7)
    assert -R(QQ(-1, 2)) == QQ(1, 2)

    R, x, y = ring('x y', ZZ)

    assert -(x**2*y - x) == -x**2*y + x

    R, x, y, z = ring('x y z', ZZ)

    assert -R(0) == 0
    assert -R(1) == -1
    assert -R(-7) == 7

    R, x, y, z = ring('x y z', QQ)

    assert -R(0) == 0
    assert -R(QQ(1, 9)) == QQ(-1, 9)
    assert -R(QQ(-7, 9)) == QQ(7, 9)


def test_PolyElement___add__():
    R, x = ring('x', ZZ)

    assert R(1) + R(0) == 1
    assert R(1) + R(1) == 2
    assert R(1) + R(2) == 3

    assert (x + 2) + R(1) == x + 3
    assert R(1) + (x + 2) == x + 3

    assert (x**2 + 2*x + 3) + (8*x**2 + 9*x + 10) == 9*x**2 + 11*x + 13

    assert (x**2 - 1) + (x - 2) == x**2 + x - 3

    f = R(0)

    assert f + 0 == 0
    assert f + 1 == 1
    assert f + x == x
    assert f + x**2 == x**2

    f = x**2 + x + 1

    assert f + 1 == x**2 + x + 2
    assert f + x == x**2 + 2*x + 1
    assert f + x**2 == 2*x**2 + x + 1

    assert f + x**3 == x**3 + x**2 + x + 1
    assert f + x**4 == x**4 + x**2 + x + 1
    assert f + x**5 == x**5 + x**2 + x + 1
    assert f + x**6 == x**6 + x**2 + x + 1

    assert f - x**2 == x + 1

    f = x**2 - 1

    assert f + 2*x**4 == 2*x**4 + x**2 - 1

    R, x = ring('x', QQ)

    assert R(0) + R(0) == 0
    assert R(QQ(1, 2)) + R(0) == QQ(1, 2)
    assert R(0) + R(QQ(1, 2)) == QQ(1, 2)
    assert R(QQ(1, 4)) + R(QQ(1, 4)) == QQ(1, 2)
    assert R(QQ(1, 4)) + R(QQ(1, 2)) == QQ(3, 4)

    assert (x/2 + QQ(2, 3)) + R(1) == x/2 + QQ(5, 3)
    assert R(1) + (x/2 + QQ(2, 3)) == x/2 + QQ(5, 3)

    assert ((x**2/7 + 2*x/7 + QQ(3, 7)) +
            (8*x**2/7 + 9*x/7 + QQ(10, 7))) == 9*x**2/7 + 11*x/7 + QQ(13, 7)

    R, x, y = ring('x y', ZZ)

    pytest.raises(CoercionFailed, lambda: R.convert(EX(pi)))

    p = x**4 + 2*y
    m = (1, 2)
    p1 = p._iadd_term((m, 5))

    assert p is p1 and p1 == x**4 + 5*x*y**2 + 2*y

    p2 = p._iadd_term(((0, 1), 2))

    assert p == p2 and p2 == x**4 + 5*x*y**2 + 4*y

    p3 = p._iadd_term(((0, 1), -4))

    assert p == p3 and p3 == x**4 + 5*x*y**2

    p = x
    p1 = p._iadd_term((m, 5))

    assert p is not p1 and p1 == 5*x*y**2 + x

    assert (x + y - 1) + 1 == x + y

    f, g = (x + y)**2, (x - y)**2

    assert f + g == 2*x**2 + 2*y**2

    f, g = x**2 + y, x**2*y + x

    assert f + g == x**2*y + x**2 + x + y

    R, x, y, z = ring('x y z', ZZ)

    assert R(0) + R(0) == 0
    assert R(1) + R(0) == 1
    assert R(0) + R(1) == 1
    assert R(2) + R(1) == 3
    assert R(1) + R(2) == 3

    p1 = x**4 + 2*y
    p2 = y + z
    m = (1, 2, 3)

    p = p1._iadd_poly_term(p2, (m, 3))

    assert p is p1 and p == x**4 + 3*x*y**3*z**3 + 3*x*y**2*z**4 + 2*y

    f = f_polys()[0]

    assert f + 0 == f

    f = x*y + 1

    assert f + 2*x**2 == 2*x**2 + x*y + 1

    R, x, y, z = ring('x y z', QQ)

    assert R(0) + R(0) == 0
    assert R(QQ(1, 2)) + R(0) == QQ(1, 2)
    assert R(0) + R(QQ(1, 2)) == QQ(1, 2)
    assert R(QQ(2, 7)) + R(QQ(1, 7)) == QQ(3, 7)
    assert R(QQ(1, 7)) + R(QQ(2, 7)) == QQ(3, 7)

    f = f.set_ring(R)/7

    assert f + 0 == f

    Rt, t = ring('t', ZZ)
    Ruv, u, v = ring('u v', ZZ)
    Rxyz, x, y, z = ring('x y z', Ruv)

    assert dict(+x) == dict(x)

    assert dict(x + 3*y) == {(1, 0, 0): 1, (0, 1, 0): 3}

    assert dict(u + x) == dict(x + u) == {(1, 0, 0): 1, (0, 0, 0): u}
    assert dict(u + x*y) == dict(x*y + u) == {(1, 1, 0): 1, (0, 0, 0): u}
    assert dict(u + x*y + z) == dict(x*y + z + u) == {(1, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): u}

    assert dict(u*x + x) == dict(x + u*x) == {(1, 0, 0): u + 1}
    assert dict(u*x + x*y) == dict(x*y + u*x) == {(1, 1, 0): 1, (1, 0, 0): u}
    assert dict(u*x + x*y + z) == dict(x*y + z + u*x) == {(1, 1, 0): 1, (0, 0, 1): 1, (1, 0, 0): u}

    pytest.raises(TypeError, lambda: t + x)
    pytest.raises(TypeError, lambda: x + t)
    pytest.raises(TypeError, lambda: t + u)
    pytest.raises(TypeError, lambda: u + t)

    Fuv, u, v = field('u v', ZZ)
    Rxyz, x, y, z = ring('x y z', Fuv)

    assert u + (x - u) == x
    assert dict(u + x) == dict(x + u) == {(1, 0, 0): 1, (0, 0, 0): u}

    Rxyz, x, y, z = ring('x y z', EX)

    assert dict(EX(pi) + x*y*z) == dict(x*y*z + EX(pi)) == {(1, 1, 1): EX(1), (0, 0, 0): EX(pi)}


def test_PolyElement___sub__():
    R, x = ring('x', ZZ)

    assert R(0) - R(0) == 0
    assert R(1) - R(0) == +1
    assert R(0) - R(1) == -1
    assert R(1) - R(1) == 0
    assert R(1) - R(2) == -1

    assert (x + 2) - R(1) == +x + 1
    assert R(1) - (x + 2) == -x - 1

    assert (3*x**2 + 2*x + 1) - (8*x**2 + 9*x + 10) == -5*x**2 - 7*x - 9

    assert (x**2 - 1) - (x - 2) == x**2 - x + 1

    R, x = ring('x', QQ)

    assert R(0) - R(0) == 0
    assert R(QQ(1, 2)) - R(0) == QQ(+1, 2)
    assert R(0) - R(QQ(1, 2)) == QQ(-1, 2)
    assert R(QQ(1, 3)) - R(QQ(1, 3)) == 0
    assert R(QQ(1, 3)) - R(QQ(2, 3)) == QQ(-1, 3)

    assert (x/7 + QQ(2, 7)) - R(1) == +x/7 - QQ(5, 7)
    assert R(1) - (x/7 + QQ(2, 7)) == -x/7 + QQ(5, 7)

    assert ((3*x**2/7 + 2*x/7 + QQ(1, 7)) -
            (8*x**2/7 + 9*x/7 + QQ(10, 7))) == -5*x**2/7 - x - QQ(9, 7)

    R, x, y = ring('x y', ZZ)

    assert (x + y + 1) - 1 == x + y

    f, g = x + y**2, x*y + y**2

    assert f - g == x - x*y

    f, g = 4, x + y

    assert f - g == -x - y + 4

    f, g = x**2 - 2, y**2

    assert f - g == x**2 - y**2 - 2
    assert g - f == 2 + y**2 - x**2

    f, g = x**2 + y, x**2*y + x

    assert f - g == -x**2*y + x**2 - x + y

    R, x, y, z = ring('x y z', ZZ)

    assert R(0) - R(0) == 0
    assert R(1) - R(0) == +1
    assert R(0) - R(1) == -1
    assert R(2) - R(1) == +1
    assert R(1) - R(2) == -1

    R, x, y, z = ring('x y z', QQ)

    assert R(0) - R(0) == 0
    assert R(QQ(1, 2)) - R(0) == QQ(+1, 2)
    assert R(0) - R(QQ(1, 2)) == QQ(-1, 2)
    assert R(QQ(2, 7)) - R(QQ(1, 7)) == QQ(+1, 7)
    assert R(QQ(1, 7)) - R(QQ(2, 7)) == QQ(-1, 7)

    Rt, t = ring('t', ZZ)
    Ruv, u, v = ring('u v', ZZ)
    Rxyz, x, y, z = ring('x y z', Ruv)

    assert u - x == -x + u
    assert (x + u) - 2*u == x - u

    assert dict(x - 3*y) == {(1, 0, 0): 1, (0, 1, 0): -3}

    assert dict(-u + x) == dict(x - u) == {(1, 0, 0): 1, (0, 0, 0): -u}
    assert dict(-u + x*y) == dict(x*y - u) == {(1, 1, 0): 1, (0, 0, 0): -u}
    assert dict(-u + x*y + z) == dict(x*y + z - u) == {(1, 1, 0): 1, (0, 0, 1): 1, (0, 0, 0): -u}

    assert dict(-u*x + x) == dict(x - u*x) == {(1, 0, 0): -u + 1}
    assert dict(-u*x + x*y) == dict(x*y - u*x) == {(1, 1, 0): 1, (1, 0, 0): -u}
    assert dict(-u*x + x*y + z) == dict(x*y + z - u*x) == {(1, 1, 0): 1, (0, 0, 1): 1, (1, 0, 0): -u}

    pytest.raises(TypeError, lambda: t - x)
    pytest.raises(TypeError, lambda: x - t)
    pytest.raises(TypeError, lambda: t - u)
    pytest.raises(TypeError, lambda: u - t)

    Fuv, u, v = field('u v', ZZ)
    Rxyz, x, y, z = ring('x y z', Fuv)

    assert dict(-u + x) == dict(x - u) == {(1, 0, 0): 1, (0, 0, 0): -u}

    Rxyz, x, y, z = ring('x y z', EX)

    assert dict(-EX(pi) + x*y*z) == dict(x*y*z - EX(pi)) == {(1, 1, 1): EX(1), (0, 0, 0): -EX(pi)}


def test_PolyElement___mul__():
    R, x = ring('x', FF(7))

    assert (2*x + 1)*(3*x + 4) == 6*x**2 + 4*x + 4

    R, x = ring('x', ZZ)

    f = R(0)

    assert f*R(0) == 0
    assert f*R(1) == 0
    assert f*ZZ(2) == 0

    f = R(1)

    assert f*R(0) == 0
    assert f*R(1) == 1

    assert R(5)*R(7) == 35

    f = x**2 + 2*x - 1

    assert f*ZZ(3) == 3*x**2 + 6*x - 3

    f = x**2 + 2*x + 3

    assert f*ZZ(0) == 0
    assert f*ZZ(2) == 2*x**2 + 4*x + 6

    assert (x - 2)*(x + 2) == x**2 - 4

    f = 3*x**5 + 6*x**2 + x + 2
    g = 4*x**3 + x
    h = 12*x**8 + 3*x**6 + 24*x**5 + 4*x**4 + 14*x**3 + x**2 + 2*x

    assert f*g == h
    assert g*f == h

    f = 2*x**4 + x + 7
    h = 4*x**8 + 4*x**5 + 28*x**4 + x**2 + 14*x + 49

    assert f*f == h

    p1 = 1 + 2*x - 3*x**2 + x**5
    p2 = 1 - x + 5*x**2 + x**3

    assert p1**25 * p2**25 == (p1*p2)**25
    assert p1**25 * p2**20 == p1**5 * (p1*p2)**20

    R, x = ring('x', QQ)

    assert R(0)*R(0) == 0
    assert R(0)*R(QQ(1, 2)) == 0
    assert R(QQ(1, 2))*R(0) == 0
    assert R(QQ(1, 2))*R(QQ(4, 7)) == QQ(2, 7)
    assert R(QQ(5, 7))*R(QQ(3, 7)) == QQ(15, 49)

    R, x, y = ring('x y', FF(5))

    assert (2*x + 1)*(3*x + 4) == x**2 + x + 4

    R, x, y = ring('x y', ZZ)

    assert (x*y + 1)*x == x**2*y + x

    assert (2*x + 2*y)*ZZ(3) == 6*x + 6*y

    R, x, y = ring('x y', QQ)

    f, g = x + y, x - y

    assert f*g == x**2 - y**2

    f, g = 4, x + y

    assert f*g == 4*x + 4*y

    R, x, y, z = ring('x y z', ZZ)

    assert R(0)*R(0) == 0
    assert R(1)*R(0) == 0
    assert R(0)*R(1) == 0
    assert R(2)*R(1) == 2
    assert R(1)*R(2) == 2

    f = f_polys()[0]

    assert (f*ZZ(2) ==
            2*x**2*y*z**2 + 4*x**2*y*z + 6*x**2*y + 4*x**2 + 6*x +
            8*y**2*z**2 + 10*y**2*z + 12*y**2 + 2*y*z**2 + 4*y*z + 2*y + 2)

    R, x, y, z = ring('x y z', QQ)

    assert R(0)*R(0) == 0
    assert R(QQ(1, 2))*R(0) == 0
    assert R(0)*R(QQ(1, 2)) == 0
    assert R(QQ(2, 7))*R(QQ(1, 3)) == QQ(2, 21)
    assert R(QQ(1, 7))*R(QQ(2, 3)) == QQ(2, 21)

    f = f.set_ring(R)/7

    assert (f*QQ(1, 2) ==
            x**2*y*z**2/14 + x**2*y*z/7 + 3*x**2*y/14 + x**2/7 + 3*x/14 +
            2*y**2*z**2/7 + 5*y**2*z/14 + 3*y**2/7 + y*z**2/14 + y*z/7 +
            y/14 + QQ(1, 14))

    Rt, t = ring('t', ZZ)
    Ruv, u, v = ring('u v', ZZ)
    Rxyz, x, y, z = ring('x y z', Ruv)

    assert dict(u*x) == dict(x*u) == {(1, 0, 0): u}

    assert dict(2*u*x + z) == dict(x*2*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*x + z) == dict(2*x*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*x + z) == dict(x*2*u + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*x*2 + z) == dict(x*u*2 + z) == {(1, 0, 0): 2*u, (0, 0, 1): 1}

    assert dict(2*u*x*y + z) == dict(x*y*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*x*y + z) == dict(2*x*y*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*x*y + z) == dict(x*y*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*x*y*2 + z) == dict(x*y*u*2 + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}

    assert dict(2*u*y*x + z) == dict(y*x*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*2*y*x + z) == dict(2*y*x*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(2*u*y*x + z) == dict(y*x*2*u + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}
    assert dict(u*y*x*2 + z) == dict(y*x*u*2 + z) == {(1, 1, 0): 2*u, (0, 0, 1): 1}

    assert dict(3*u*(x + y) + z) == dict((x + y)*3*u + z) == {(1, 0, 0): 3*u, (0, 1, 0): 3*u, (0, 0, 1): 1}

    pytest.raises(TypeError, lambda: t*x + z)
    pytest.raises(TypeError, lambda: x*t + z)
    pytest.raises(TypeError, lambda: t*u + z)
    pytest.raises(TypeError, lambda: u*t + z)

    Fuv, u, v = field('u v', ZZ)
    Rxyz, x, y, z = ring('x y z', Fuv)

    assert dict(u*x) == dict(x*u) == {(1, 0, 0): u}

    R, x, y, z = ring('x y z', EX)

    assert dict(EX(pi)*x*y*z) == dict(x*y*z*EX(pi)) == {(1, 1, 1): EX(pi)}
    assert (x + 2*y)*R(0) == 0


def test_PolyElement_misc_ariths():
    R, x = ring('x', ZZ)

    f, g, h = x**2 + 2*x + 3, 3*x**2 + 2*x + 1, x + 2

    assert f + g*h == 3*x**3 + 9*x**2 + 7*x + 5
    assert f - g*h == -3*x**3 - 7*x**2 - 3*x + 1

    f, g, h = x**2 - 1, x - 2, x + 2

    assert f + g*h == 2*x**2 - 5
    assert f - g*h == 3

    R, x, y = ring('x y', ZZ)

    f, g, h = x*y + 2*x + 3, 3*x + 2*y + 1, x + 2

    assert f + g*h == 3*x**2 + 3*x*y + 9*x + 4*y + 5
    assert f - g*h == -3*x**2 - x*y - 5*x - 4*y + 1

    f, g, h = x**2 + y, x, x + 2

    assert f + g*h == 2*x**2 + y + 2*x
    assert f - g*h == -2*x + y


def test_PolyElement___floordiv__truediv__():
    R, x = ring('x', ZZ)

    f, g = x**2 + 2*x + 3, R(0)

    pytest.raises(ZeroDivisionError, lambda: divmod(f, g))
    pytest.raises(ZeroDivisionError, lambda: f % g)
    pytest.raises(ZeroDivisionError, lambda: f // g)
    pytest.raises(ZeroDivisionError, lambda: f.exquo(g))

    f, g = x**2 + 1, 2*x - 4
    q, r = R(0), x**2 + 1

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    assert divmod(R.zero, f) == (0, 0)

    f, g = 3*x**3 + x**2 + x + 5, 5*x**2 - 3*x + 1
    q, r = R(0), f

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = 5*x**4 + 4*x**3 + 3*x**2 + 2*x + 1, x**2 + 2*x + 3
    q, r = 5*x**2 - 6*x, 20*x + 1

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = 5*x**5 + 4*x**4 + 3*x**3 + 2*x**2 + x, x**4 + 2*x**3 + 9
    q, r = 5*x - 6, 15*x**3 + 2*x**2 - 44*x + 54

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    R, x = ring('x', QQ)

    f, g = x**2 + 2*x + 3, R(0)

    pytest.raises(ZeroDivisionError, lambda: divmod(f, g))

    f, g = x**2 + 1, 2*x - 4
    q, r = x/2 + 1, R(5)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = 3*x**3 + x**2 + x + 5, 5*x**2 - 3*x + 1
    q, r = 3*x/5 + QQ(14, 25), 52*x/25 + QQ(111, 25)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    R, x = ring('x', RR)

    pytest.raises(PolynomialDivisionFailed,
                  lambda: divmod(R(2.0), R(-1.8438812457236466e-19)))

    R, x, y = ring('x y', ZZ)

    f, g = x*y + 2*x + 3, R(0)

    pytest.raises(ZeroDivisionError, lambda: divmod(f, g))
    pytest.raises(ZeroDivisionError, lambda: f % g)
    pytest.raises(ZeroDivisionError, lambda: f // g)
    pytest.raises(ZeroDivisionError, lambda: f.exquo(g))

    f, g = x**2 - y**2, x - y
    q, r = x + y, R(0)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    assert f.exquo(g) == q

    f, g = x**2 + y**2, x - y
    q, r = x + y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, -x + y
    q, r = -x - y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, 2*x - 2*y
    q, r = R(0), f

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + x*y, 2*x + 2

    assert divmod(f, g) == (0, f)

    R, x, y = ring('x y', QQ)

    f, g = x*y + 2*x + 3, R(0)

    pytest.raises(ZeroDivisionError, lambda: divmod(f, g))

    f, g = x**2 - y**2, x - y
    q, r = x + y, R(0)

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q
    assert f.exquo(g) == q

    f, g = x**2 + y**2, x - y
    q, r = x + y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, -x + y
    q, r = -x - y, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    f, g = x**2 + y**2, 2*x - 2*y
    q, r = x/2 + y/2, 2*y**2

    assert divmod(f, g) == (q, r)
    assert f % g == r
    assert f // g == q

    pytest.raises(ExactQuotientFailed, lambda: f.exquo(g))

    pytest.raises(ZeroDivisionError, lambda: f.quo_ground(0))
    pytest.raises(ZeroDivisionError, lambda: f.quo_term(((1, 1), 0)))
    pytest.raises(ZeroDivisionError, lambda: f.exquo_ground(0))

    assert R.zero.exquo_ground(2) == 0

    assert R.zero.quo_term(((1, 0), 1)) == 0
    assert g.quo_term((R.zero_monom, 2)) == x - y
    assert f.quo_term(((1, 0), 2)) == x/2

    R, x, y, z = ring('x y z', ZZ)

    f, g, q, r = R(1), 2*x + 1, 0, 1

    assert divmod(f, g) == (q, r)
    assert f // g == q
    assert f % g == r

    assert (2*x**2 - 4)/2 == x**2 - 2
    assert (2*x**2 - 3)/2 == x**2 - 2

    assert (x**2 - 1)//x == x
    assert (x**2 - x)//x == x - 1

    assert (x**2 - 1)//(2*x) == 0
    assert (x**2 - x)//(x - 1) == x

    assert len((x**2/3 + y**3/4 + z**4/5).items()) == 0

    R, x, y, z = ring('x y z', QQ)

    assert len((x**2/3 + y**3/4 + z**4/5).items()) == 3

    pytest.raises(ZeroDivisionError, lambda: x/0)

    Rt, t = ring('t', ZZ)
    Ruv, u, v = ring('u v', ZZ)
    Rxyz, x, y, z = ring('x y z', Ruv)

    assert dict((u**2*x + u)/u) == {(1, 0, 0): u, (0, 0, 0): 1}

    pytest.raises(TypeError, lambda: u//(u**2*x + u))
    pytest.raises(TypeError, lambda: t//x)
    pytest.raises(TypeError, lambda: x//t)
    pytest.raises(TypeError, lambda: t//u)
    pytest.raises(TypeError, lambda: u//t)

    assert divmod(x, u) == (0, x)
    assert x % u == x
    assert x // u == 0

    pytest.raises(TypeError, lambda: divmod(u, x))
    pytest.raises(TypeError, lambda: u % x)
    pytest.raises(TypeError, lambda: u // x)
    pytest.raises(TypeError, lambda: divmod(u, t))
    pytest.raises(TypeError, lambda: u % t)
    pytest.raises(TypeError, lambda: u // t)
    pytest.raises(TypeError, lambda: divmod(u, sqrt(2)))
    pytest.raises(TypeError, lambda: u % sqrt(2))
    pytest.raises(TypeError, lambda: u // sqrt(2))


def test_PolyElement_quo_term():
    R, x = ring('x', ZZ)

    assert R(0).quo_term(((3,), ZZ(1))) == 0
    assert (x**3).quo_term(((3,), ZZ(1))) == 1

    f = x**4 + 2*x**3 + 3*x**2 + 4*x + 5

    assert R(0).quo_term(((5,), ZZ(1))) == 0

    assert f.quo_term(((1,), ZZ(1))) == x**3 + 2*x**2 + 3*x + 4
    assert f.quo_term(((3,), ZZ(1))) == x + 2
    assert f.quo_term(((5,), ZZ(1))) == 0

    f = x**4 + x**2

    assert f.quo_term(((2,), ZZ(1))) == x**2 + 1

    f += 2

    assert f.quo_term(((2,), ZZ(1))) == x**2 + 1

    R, x, y, z = ring('x y z', ZZ)

    f = x**3*y**4*z

    assert f.quo_term(((1, 2, 0), 1)) == x**2*y**2*z
    assert f.quo_term(((1, 2, 2), 1)) == 0

    R, x, y, z = ring('x y z', QQ)

    f = x**3*y**4*z

    assert f.quo_term(((1, 2, 0), 1)) == x**2*y**2*z
    assert f.quo_term(((1, 2, 2), 1)) == 0


def test_PolyElement_mul_term():
    R, x = ring('x', ZZ)

    f = R(0)

    assert f*R([((3,), ZZ(2))]) == 0

    f = x + 1

    assert f*R([((3,), 0)]) == 0

    f = x**2 + 2*x + 3

    assert f*R([((0,), ZZ(2))]) == 2*x**2 + 4*x + 6
    assert f*R([((1,), ZZ(2))]) == 2*x**3 + 4*x**2 + 6*x
    assert f*R([((2,), ZZ(2))]) == 2*x**4 + 4*x**3 + 6*x**2
    assert f*R([((3,), ZZ(2))]) == 2*x**5 + 4*x**4 + 6*x**3

    assert (x**2 - 1)*R([((2,), ZZ(3))]) == 3*x**4 - 3*x**2

    R, x, y = ring('x y', ZZ)

    assert R(0)*R([((3, 0), ZZ(2))]) == 0
    assert R(1)*R([((3, 0), ZZ(0))]) == 0

    f = x*y + 2*x + 3

    assert f*R([((2, 0), ZZ(2))]) == 2*x**3*y + 4*x**3 + 6*x**2
    assert (x**2*y + x)*R([((2, 1), 3)]) == 3*x**4*y**2 + 3*x**3*y

    R, x, y = ring('x y', QQ)

    assert R(0)*R([((3, 0), QQ(2, 3))]) == 0
    assert R(QQ(1, 2))*R([((3, 0), QQ(0))]) == 0

    f = x*y/5 + 2*x/5 + QQ(3, 5)

    assert f*R([((2, 0), QQ(2, 3))]) == 2*x**3*y/15 + 4*x**3/15 + 2*x**2/5


def test_PolyElement_mul_monom():
    R, x = ring('x', ZZ)

    assert R(0)*R([((3,), ZZ(1))]) == 0
    assert R(1)*R([((3,), ZZ(1))]) == x**3

    f = x**4 + 2*x**3 + 3*x**2 + 4*x + 5

    assert f*R([((1,), ZZ(1))]) == x**5 + 2*x**4 + 3*x**3 + 4*x**2 + 5*x
    assert f*R([((2,), ZZ(1))]) == x**6 + 2*x**5 + 3*x**4 + 4*x**3 + 5*x**2

    f = x**2 + 1

    assert f*R([((2,), ZZ(1))]) == x**4 + x**2


def test_PolyElement___pow__():
    R, x = ring('x', FF(5))

    f, g = 3*x**2 + 2*x + 4, x + 1

    assert f**3 == 2*x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 2*x**2 + x + 4
    assert pow(f, 3, g) == 0

    R, x = ring('x', FF(7))

    assert (3*x + 4)**2 == 2*x**2 + 3*x + 2

    R, x = ring('x', FF(11))

    f = x + 1

    assert f**11 == x**11 + 1

    f = x**10 + 1

    assert f**11 == x**110 + 1

    f = x**4 + x + 8

    assert f**0 == 1
    assert f**1 == f
    assert f**2 == x**8 + 2*x**5 + 5*x**4 + x**2 + 5*x + 9
    assert f**5 == (x**20 + 5*x**17 + 7*x**16 + 10*x**14 + 6*x**13 +
                    2*x**12 + 10*x**11 + 9*x**10 + 6*x**9 + 10*x**8 +
                    6*x**7 + 6*x**6 + 5*x**4 + 2*x**3 + 5*x**2 + 9*x + 10)
    assert f**8 == (x**32 + 8*x**29 + 9*x**28 + 6*x**26 + 8*x**25 +
                    10*x**24 + x**23 + 2*x**22 + 5*x**21 + 10*x**20 +
                    7*x**19 + 7*x**18 + 9*x**17 + x**16 + 2*x**15 +
                    6*x**12 + 2*x**11 + 5*x**10 + 2*x**9 + 5*x**8 +
                    7*x**7 + 7*x**6 + 9*x**5 + 10*x**4 + 10*x**3 +
                    7*x**2 + 5*x + 5)
    assert f**45 == (x**180 + x**177 + 8*x**176 + 4*x**147 + 4*x**144 +
                     10*x**143 + 10*x**136 + 10*x**133 + 3*x**132 +
                     6*x**114 + 6*x**111 + 4*x**110 + 8*x**103 + 8*x**100 +
                     9*x**99 + 10*x**92 + 10*x**89 + 3*x**88 + 4*x**81 +
                     4*x**78 + 10*x**77 + 8*x**70 + 8*x**67 + 9*x**66 +
                     9*x**59 + 9*x**56 + 6*x**55 + 3*x**48 + 3*x**45 +
                     2*x**44 + 10*x**37 + 10*x**34 + 3*x**33 + 10*x**26 +
                     10*x**23 + 3*x**22 + 2*x**15 + 2*x**12 + 5*x**11 +
                     4*x**4 + 4*x + 10)

    g = 2*x**2 + 7

    assert pow(f, 0, g) == 1
    assert pow(f, 1, g) == x + 1
    assert pow(f, 2, g) == 2*x + 3
    assert pow(f, 5, g) == 7*x + 8
    assert pow(f, 8, g) == x + 5
    assert pow(f, 45, g) == 5*x + 4

    R, x = ring('x', FF(8))

    f = x + 1

    assert f**0 == 1
    assert f**1 == f
    assert f**2 == f*f == x**2 + 1
    assert f**3 == f*f*f == x**3 + x**2 + x + 1
    assert f**4 == f*f*f*f == x**4 + 1
    assert f**5 == f*f*f*f*f == x**5 + x**4 + x + 1
    assert f**8 == functools.reduce(operator.mul, [f]*8) == x**8 + 1

    F9 = FF(3, [2, 2, 1])
    R, x = ring('x', F9)

    f = x + F9(4)

    assert f**3 == f*f*f == x**3 + F9(8)

    R, x = ring('x', ZZ)

    assert R(0)**0 == 1
    assert R(0)**1 == 0
    assert R(0)**2 == 0
    assert R(0)**7 == 0

    assert R(2)**2 == 4
    assert (x + 2)**2 == x**2 + 4*x + 4

    assert (2*x**4 + x + 7)**2 == 4*x**8 + 4*x**5 + 28*x**4 + x**2 + 14*x + 49

    pytest.raises(ValueError, lambda: R(1)**object())
    pytest.raises(ValueError, lambda: R(1)**-1)

    assert R(1)**0 == 1
    assert R(1)**1 == 1
    assert R(1)**7 == 1

    assert R(3)**0 == 1
    assert R(3)**1 == 3
    assert R(3)**7 == 2187

    assert (x**2 + 1)**2 == x**4 + 2*x**2 + 1

    f = 2*x**4 + x + 7

    assert f**0 == 1
    assert f**1 == f
    assert f**2 == 4*x**8 + 4*x**5 + 28*x**4 + x**2 + 14*x + 49
    assert f**3 == (8*x**12 + 12*x**9 + 84*x**8 + 6*x**6 +
                    84*x**5 + 294*x**4 + x**3 + 21*x**2 + 147*x + 343)

    assert (x - 2)**3 == x**3 - 6*x**2 + 12*x - 8

    f = -11200*x**4 - 2604*x**2 + 49
    g = 15735193600000000*x**16 + 14633730048000000*x**14 + 4828147466240000*x**12 \
        + 598976863027200*x**10 + 3130812416256*x**8 - 2620523775744*x**6 \
        + 92413760096*x**4 - 1225431984*x**2 + 5764801

    assert f**4 == f._pow_generic(4) == f._pow_multinomial(4) == g

    R, x = ring('x', ZZ, grlex)

    assert R.zero**0 == 1

    f = 2*x + 3

    assert f**0 == 1
    assert f**1 == f

    assert f**2 == f._pow_generic(2) == f._pow_multinomial(2) == 4*x**2 + 12*x + 9
    assert f**3 == f._pow_generic(3) == f._pow_multinomial(3) == 8*x**3 + 36*x**2 + 54*x + 27
    assert f**4 == f._pow_generic(4) == f._pow_multinomial(4) == 16*x**4 + 96*x**3 + 216*x**2 + 216*x + 81
    assert f**5 == f._pow_generic(5) == f._pow_multinomial(5) == 32*x**5 + 240*x**4 + 720*x**3 + 1080*x**2 + 810*x + 243

    pytest.raises(ValueError, lambda: f**-2)

    f = x**2 - 2*x + x**3 + 1

    assert f**5 == (x**15 + 5*x**14 - 25*x**12 + 5*x**11 + 71*x**10 -
                    60*x**9 - 85*x**8 + 170*x**7 - 60*x**6 - 112*x**5 +
                    170*x**4 - 115*x**3 + 45*x**2 - 10*x + 1)

    R, x = ring('x', QQ)

    assert R(0)**0 == 1
    assert R(0)**2 == 0

    assert R(1)**0 == 1
    assert R(1)**1 == 1
    assert R(1)**7 == 1

    assert R(QQ(3, 7))**0 == 1
    assert R(QQ(3, 7))**1 == QQ(3, 7)
    assert R(QQ(2, 3))**2 == QQ(4, 9)
    assert R(QQ(3, 7))**7 == QQ(2187, 823543)

    assert (x/3 + QQ(2, 3))**2 == x**2/9 + 4*x/9 + QQ(4, 9)

    R, x, y = ring('x y', FF(7))

    assert (3*x + 4)**2 == 2*x**2 + 3*x + 2

    R, x, y = ring('x y', ZZ)

    assert R(0)**0 == 1
    assert R(0)**1 == 0
    assert R(0)**7 == 0

    assert R(1)**0 == 1
    assert R(1)**1 == 1
    assert R(1)**7 == 1

    f = x + y**2

    assert f**2 == x**2 + 2*x*y**2 + y**4
    assert f**3 == x**3 + 3*x**2*y**2 + 3*x*y**4 + y**6

    assert (x**2 + x*y +
            y**2)**2 == x**4 + 2*x**3*y + 3*x**2*y**2 + 2*x*y**3 + y**4

    pytest.raises(ValueError, lambda: R(1)**-1)

    R, x, y = ring('x y', QQ)

    assert R(0)**0 == 1

    assert R(QQ(3, 7))**0 == 1
    assert R(QQ(3, 7))**1 == QQ(3, 7)
    assert R(QQ(3, 7))**7 == QQ(2187, 823543)

    # issues sympy/sympy#20485 and sympy/sympy#20487
    pytest.raises(ValueError, lambda: x**QQ(3, 2))
    pytest.raises(ValueError, lambda: (x + y)**QQ(3, 2))

    R, x, y, z = ring('x y z', ZZ)

    assert R(0)**2 == 0
    assert R(2)**2 == 4

    R, x, y, z = ring('x y z', ZZ, grlex)

    f = x**3*y - 2*x*y**2 - 3*z + 1
    g = x**6*y**2 - 4*x**4*y**3 - 6*x**3*y*z + 2*x**3*y + 4*x**2*y**4 + 12*x*y**2*z - 4*x*y**2 + 9*z**2 - 6*z + 1

    assert f**2 == f._pow_generic(2) == f._pow_multinomial(2) == g

    f = x**3*y - 2*x*y**2 - 3*z + x + y + 1

    assert f**4 == f._pow_generic(4)

    R, x, y, z = ring('x y z', QQ)

    assert R(0)**2 == 0
    assert R(QQ(2, 3))**2 == QQ(4, 9)


def test_PolyElement_div():
    _, u = ring('u', ZZ)
    R, x = ring('x', ZZ)

    pytest.raises(ValueError, lambda: u.div([x]))
    pytest.raises(ZeroDivisionError, lambda: R.one.div([R.zero]))

    assert R.zero.div([x**2 + 1]) == ([0], 0)

    R, x = ring('x', ZZ, grlex)

    f = x**3 - 12*x**2 - 42
    g = x - 3

    q = x**2 - 9*x - 27
    r = -123

    assert f.div([g]) == ([q], r)

    f = x**2 + 2*x + 2

    assert f.div([R(1)]) == ([f], 0)

    R, x = ring('x', QQ, grlex)

    f = x**2 + 2*x + 2

    assert f.div([R(2)]) == ([x**2/2 + x + 1], 0)

    R, x = ring('x', RR)

    pytest.raises(PolynomialDivisionFailed,
                  lambda: divmod(R(2.0), R(-1.8438812457236466e-19)))

    R, x, y = ring('x y', ZZ)

    f = x**3
    qv = [x - y**2, x - y]

    assert f.div(qv) == ([x**2 + x*y**2 + y**4, 0], y**6)

    f = x**2*y
    qv = x**2 - y, x*y - 1

    assert f.div(qv) == ([y, 0], y**2)
    assert f.div(list(reversed(qv))) == ([x, 0], x)

    qv = x - y**2, y**3 - 1

    assert f.div(qv) == ([x*y + y**3, y**2], y**2)
    assert f.div(list(reversed(qv))) == ([x, x*y + 1], y**2)

    R, x, y = ring('x y', ZZ, grlex)

    f = 4*x**2*y - 2*x*y + 4*x - 2*y + 8

    assert f.div([R(2)]) == ([2*x**2*y - x*y + 2*x - y + 4], 0)
    assert f.div([2*y]) == ([2*x**2 - x - 1], 4*x + 8)

    f = x - 1
    g = y - 1

    assert f.div([g]) == ([0], f)

    f = x*y**2 + 1
    G = [x*y + 1, y + 1]

    Q = [y, -1]
    r = 2

    assert f.div(G) == (Q, r)

    f = x**2*y + x*y**2 + y**2
    G = [x*y - 1, y**2 - 1]

    Q = [x + y, 1]
    r = x + y + 1

    assert f.div(G) == (Q, r)

    G = [y**2 - 1, x*y - 1]

    Q = [x + 1, x]
    r = 2*x + 1

    assert f.div(G) == (Q, r)


def test_PolyElement_quo_ground():
    R, x = ring('x', ZZ)

    pytest.raises(ZeroDivisionError,
                  lambda: (x**2 + 2*x + 3).quo_ground(ZZ(0)))

    assert R(0).quo_ground(ZZ(3)) == 0

    f = 6*x**2 + 2*x + 8

    assert f.quo_ground(ZZ(1)) == f
    assert f.quo_ground(ZZ(2)) == 3*x**2 + x + 4
    assert f.quo_ground(ZZ(3)) == 2*x**2 + 2

    f = 3*x**2 + 2

    assert f.quo_ground(ZZ(2)) == x**2 + 1

    R, x = ring('x', QQ)

    f = 3*x**2 + 2

    assert f.quo_ground(QQ(2)) == 3*x**2/2 + 1

    f = 6*x**2 + 2*x + 8

    assert f.quo_ground(QQ(1)) == f
    assert f.quo_ground(QQ(2)) == 3*x**2 + x + 4
    assert f.quo_ground(QQ(7)) == 6*x**2/7 + 2*x/7 + QQ(8, 7)

    R, x, y = ring('x y', ZZ)

    f = 6*x**2 + 2*x + 8

    assert f.quo_ground(ZZ(1)) == f
    assert f.quo_ground(ZZ(2)) == 3*x**2 + x + 4
    assert f.quo_ground(ZZ(3)) == 2*x**2 + 2

    f = 2*x**2*y + 3*x

    assert f.quo_ground(ZZ(2)) == x**2*y + x

    R, x, y = ring('x y', QQ)

    f = 2*x**2*y + 3*x

    assert f.quo_ground(QQ(2)) == x**2*y + 3*x/2


def test_PolyElement_exquo_ground():
    R, x = ring('x', ZZ)

    pytest.raises(ZeroDivisionError,
                  lambda: (x**2 + 2*x + 3).exquo_ground(ZZ(0)))
    pytest.raises(ExactQuotientFailed,
                  lambda: (x**2 + 2*x + 3).exquo_ground(ZZ(3)))

    assert R(0).exquo_ground(ZZ(3)) == 0

    f = 6*x**2 + 2*x + 8

    assert f.exquo_ground(ZZ(1)) == f
    assert f.exquo_ground(ZZ(2)) == 3*x**2 + x + 4

    R, x = ring('x', QQ)

    f = x**2 + 2

    assert f.exquo_ground(QQ(2)) == x**2/2 + 1

    f = 6*x**2 + 2*x + 8

    assert f.exquo_ground(QQ(1)) == f
    assert f.exquo_ground(QQ(2)) == 3*x**2 + x + 4
    assert f.exquo_ground(QQ(7)) == 6*x**2/7 + 2*x/7 + QQ(8, 7)

    R, x, y = ring('x y', ZZ)

    f = 6*x**2 + 2*x + 8

    assert f.exquo_ground(ZZ(1)) == f
    assert f.exquo_ground(ZZ(2)) == 3*x**2 + x + 4

    R, x, y = ring('x y', QQ)

    assert (x**2*y + 2*x).exquo_ground(QQ(2)) == x**2*y/2 + x


def test_PolyElement_trunc_ground():
    R, x = ring('x', ZZ)

    assert (2*x**3 + 3*x**2 + 5*x + 7).trunc_ground(ZZ(3)) == -x**3 - x + 1
    assert (x**5 + 2*x**4 + 3*x**3 + 4*x**2 +
            5*x + 6).trunc_ground(ZZ(3)) == x**5 - x**4 + x**2 - x
    assert (6*x**5 + 5*x**4 + 4*x**3 + 3*x**2 +
            2*x + 1).trunc_ground(ZZ(3)) == -x**4 + x**3 - x + 1

    R, x = ring('x', QQ)

    assert (x**5 + 2*x**4 + 3*x**3 + 4*x**2 +
            5*x + 6).trunc_ground(ZZ(3)) == x**5 + 2*x**4 + x**2 + 2*x

    R, x, y = ring('x y', ZZ)

    f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3

    assert f.trunc_ground(ZZ(3)) == -x**2 - x*y - y

    R, x, y, z = ring('x y z', ZZ)

    f = f_polys()[0]

    assert f.trunc_ground(ZZ(3)) == (x**2*y*z**2 - x**2*y*z - x**2 +
                                     y**2*z**2 - y**2*z + y*z**2 - y*z + y + 1)


def test_PolyElement_monic():
    R, x = ring('x', ZZ)

    assert (2*x + 2).monic() == x + 1

    pytest.raises(ExactQuotientFailed, lambda: (2*x + 1).monic())

    assert (3*x**2 + 6*x + 9).monic() == x**2 + 2*x + 3

    pytest.raises(ExactQuotientFailed, lambda: (3*x**2 + 4*x + 5).monic())

    R, x = ring('x', QQ)

    assert R(0).monic() == 0
    assert R(1).monic() == 1
    assert (2*x + 1).monic() == x + QQ(1, 2)
    assert (7*x**2 + x + 21).monic() == x**2 + x/7 + 3
    assert (3*x**2 + 4*x + 2).monic() == x**2 + 4*x/3 + QQ(2, 3)
    assert (x - 1).monic() == x - 1
    assert (2*x - 1).monic() == x - QQ(1, 2)

    R, x, y = ring('x y', ZZ)

    assert (3*x**2 + 6*x + 9).monic() == x**2 + 2*x + 3

    pytest.raises(ExactQuotientFailed, lambda: (3*x**2 + 4*x + 5).monic())

    f = 3*x**2*y + 6*x**2 + 3*x*y + 9*y + 3

    assert f.monic() == x**2*y + 2*x**2 + x*y + 3*y + 1

    R, x, y = ring('x y', QQ)

    assert R(0).monic() == 0
    assert R(1).monic() == 1
    assert (7*x**2 + x + 21).monic() == x**2 + x/7 + 3

    f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3

    assert f.monic() == x**2*y + 8/3*x**2 + 5/3*x*y + 2*x + 2/3*y + 1


def test_PolyElement_content():
    R, x = ring('x', ZZ)

    assert R(0).content() == 0
    assert R(+1).content() == 1
    assert R(-1).content() == -1
    assert (x + 1).content() == 1
    assert (2*x + 2).content() == 2
    assert (x**2 + 2*x + 1).content() == 1
    assert (2*x**2 + 4*x + 2).content() == 2
    assert (6*x**2 + 8*x + 12).content() == 2

    R, x = ring('x', QQ)

    assert (6*x**2 + 8*x + 12).content() == 2

    assert (2*x/3 + QQ(4, 9)).content() == QQ(2, 9)
    assert (2*x/3 + QQ(4, 5)).content() == QQ(2, 15)

    R, x, y = ring('x y', ZZ)

    assert R(0).content() == 0
    assert R(+1).content() == 1
    assert R(-1).content() == -1
    assert (x + 1).content() == 1
    assert (2*x + 2).content() == 2
    assert (x**2 + 2*x + 1).content() == 1
    assert (2*x**2 + 4*x + 2).content() == 2

    R = R.eject(y)

    assert R(-2).content() == -2

    f, g, F = (3*y**2 + 2*y + 1).eject(y), R(1), R(0)

    for i in range(5):
        g *= f
        F += R.x**i*g

    assert F.content() == f.inject().drop(x)

    f = 2*x*y + 6*x + 4*y + 12

    assert f.content() == 2
    assert f.eject(y).content() == (2*y + 6).drop(x)

    R, x, y = ring('x y', QQ)

    assert R(0).content() == 0
    assert (2*x/3 + QQ(4, 9)).content() == QQ(2, 9)
    assert (2*x/3 + QQ(4, 5)).content() == QQ(2, 15)

    f = 2*x*y + 6*x + 4*y + 12

    assert f.content() == 2

    R, x, y, z = ring('x y z', ZZ)

    f = f_polys()[0]

    assert f.content() == 1
    assert (2*f).content() == 2

    f = f_polys()[1]

    assert f.content() == 1
    assert (3*f).content() == 3

    f = f_polys()[2]

    assert f.content() == 1
    assert (4*f).content() == 4

    f = f_polys()[3]

    assert f.content() == 1
    assert (5*f).content() == 5

    f = f_polys()[4]

    assert f.content() == -1
    assert (6*f).content() == -6

    assert f.eject(y, z).content() == -1

    f = f_polys()[5]

    assert f.content() == -1
    assert (7*f).content() == -7

    assert f.eject(y, z).content() == -1

    R, x, y, z, t = ring('x y z t', ZZ)

    f = f_polys()[6]

    assert f.content() == 1
    assert (8*f).content() == 8

    assert f.eject(y, z, t).content() == 1


def test_PolyElement_primitive():
    R, x = ring('x', ZZ)

    assert R(0).primitive() == (0, 0)  # issue sympy/sympy#21319
    assert R(1).primitive() == (1, 1)
    assert (x + 1).primitive() == (1, x + 1)
    assert (2*x + 2).primitive() == (2, x + 1)
    assert (x**2 + 2*x + 1).primitive() == (1, x**2 + 2*x + 1)
    assert (2*x**2 + 4*x + 2).primitive() == (2, x**2 + 2*x + 1)
    assert (6*x**2 + 8*x + 12).primitive() == (2, 3*x**2 + 4*x + 6)

    R, x = ring('x', QQ)

    assert R(0).primitive() == (0, 0)
    assert R(1).primitive() == (1, 1)
    assert (x + 1).primitive() == (1, x + 1)
    assert (2*x + 2).primitive() == (2, x + 1)
    assert (x**2 + 2*x + 1).primitive() == (1, x**2 + 2*x + 1)
    assert (2*x**2 + 4*x + 2).primitive() == (2, x**2 + 2*x + 1)
    assert (6*x**2 + 8*x + 12).primitive() == (2, 3*x**2 + 4*x + 6)

    assert (2*x/3 + QQ(4, 9)).primitive() == (QQ(2, 9), 3*x + 2)
    assert (2*x/3 + QQ(4, 5)).primitive() == (QQ(2, 15), 5*x + 6)

    R, x, y = ring('x y', ZZ)

    assert R(0).primitive() == (0, 0)
    assert R(2).primitive() == (2, 1)

    R = R.eject(y)

    assert R(0).primitive() == (0, 0)
    assert R(1).primitive() == (1, 1)

    f, g, F = (3*y**2 + 2*y + 1).eject(y), R(1), R(0)

    for i in range(5):
        g *= f
        F += R.x**i*g

    assert F.primitive() == (f.inject().drop(x), F // f)

    f = 2*x*y + 6*x + 4*y + 12

    assert f.primitive() == (2, x*y + 3*x + 2*y + 6)
    assert f.eject(y).primitive() == ((2*y + 6).drop(x), (x + 2).eject(y))

    R, x, y = ring('x y', QQ)

    assert R(0).primitive() == (0, 0)
    assert R(2).primitive() == (2, 1)

    assert (2*x/3 + QQ(4, 9)).primitive() == (QQ(2, 9), 3*x + 2)
    assert (2*x/3 + QQ(4, 5)).primitive() == (QQ(2, 15), 5*x + 6)
    assert (-3*x/4 + y + QQ(11, 8)).primitive() == (QQ(-1, 8), 6*x - 8*y - 11)

    f = 2*x*y + 6*x + 4*y + 12

    assert f.primitive() == (2, x*y + 3*x + 2*y + 6)

    R, x, y, z = ring('x y z', ZZ)

    f = f_polys()[0]

    assert f.primitive() == (1, f)
    assert (2*f).primitive() == (2, f)

    f = f_polys()[1]

    assert f.primitive() == (1, f)
    assert (3*f).primitive() == (3, f)

    f = f_polys()[2]

    assert f.primitive() == (1, f)
    assert (4*f).primitive() == (4, f)

    f = f_polys()[3]

    assert f.primitive() == (1, f)
    assert (5*f).primitive() == (5, f)

    f = f_polys()[4]

    assert f.primitive() == (-1, -f)
    assert (6*f).primitive() == (-6, -f)

    assert f.eject(y, z).primitive() == (-1, -f.eject(y, z))

    f = f_polys()[5]

    assert f.primitive() == (-1, -f)
    assert (7*f).primitive() == (-7, -f)

    assert f.eject(y, z).primitive() == (-1, -f.eject(y, z))

    R, x, y, z = ring('x y z', QQ)

    assert (3*x + 2).primitive() == (1, 3*x + 2)
    assert (4*x + 2).primitive() == (2, 2*x + 1)

    R, x, y, z, t = ring('x y z t', ZZ)

    f = f_polys()[6]

    assert f.primitive() == (1, f)
    assert (8*f).primitive() == (8, f)

    assert f.eject(y, z, t).primitive() == (1, f.eject(y, z, t))


def test_PolyElement_clear_denoms():
    R, x = ring('x', QQ)

    assert R(0).clear_denoms() == (1, 0)
    assert R(1).clear_denoms() == (1, 1)
    assert R(7).clear_denoms() == (1, 7)
    assert R(QQ(7, 3)).clear_denoms() == (3, 7)

    f = 3*x**2 + x

    assert f.clear_denoms() == (1, 3*x**2 + x)
    assert f.clear_denoms(convert=True) == (1, (3*x**2 + x).set_domain(ZZ))

    f = x**2 + x/2

    assert f.clear_denoms() == (2, 2*x**2 + x)
    assert f.clear_denoms(convert=True) == (2, (2*x**2 + x).set_domain(ZZ))

    f = x/2 + QQ(1, 3)

    assert f.clear_denoms() == (6, 3*x + 2)
    assert f.clear_denoms(convert=True) == (6, (3*x + 2).set_domain(ZZ))

    R, x = ring('x', EX)

    assert (3*x/2 + Rational(9, 4)).clear_denoms() == (4, 6*x + 9)

    assert R(7).clear_denoms() == (1, 7)

    a = EX.to_expr(EX('a'))

    assert (sin(a)/a*x).clear_denoms() == (a, x*sin(a))

    R, x, y = ring('x y', QQ)

    assert R(0).clear_denoms() == (1, 0)
    assert R(1).clear_denoms() == (1, 1)
    assert R(7).clear_denoms() == (1, 7)

    assert R(QQ(7, 3)).clear_denoms() == (3, 7)

    f = 3*x**2 + x

    assert f.clear_denoms() == (1, 3*x**2 + x)
    assert f.clear_denoms(convert=True) == (1, (3*x**2 + x).set_domain(ZZ))

    f = x**2 + x/2

    assert f.clear_denoms() == (2, 2*x**2 + x)
    assert f.clear_denoms(convert=True) == (2, (2*x**2 + x).set_domain(ZZ))

    f = x/2 + y/3 + 1

    assert f.clear_denoms() == (6, 3*x + 2*y + 6)
    assert f.clear_denoms(convert=True) == (6, (3*x + 2*y + 6).set_domain(ZZ))

    F = [x - 17824537287975195925064602467992950991718052713078834557692023531499318507213727406844943097*y**7/413954288007559433755329699713866804710749652268151059918115348815925474842910720000
           - 4882321164854282623427463828745855894130208215961904469205260756604820743234704900167747753*y**6/12936071500236232304854053116058337647210926633379720622441104650497671088840960000
           - 36398103304520066098365558157422127347455927422509913596393052633155821154626830576085097433*y**5/25872143000472464609708106232116675294421853266759441244882209300995342177681920000
           - 168108082231614049052707339295479262031324376786405372698857619250210703675982492356828810819*y**4/58212321751063045371843239022262519412449169850208742800984970927239519899784320000
           - 5694176899498574510667890423110567593477487855183144378347226247962949388653159751849449037*y**3/1617008937529529038106756639507292205901365829172465077805138081312208886105120000
           - 154482622347268833757819824809033388503591365487934245386958884099214649755244381307907779*y**2/60637835157357338929003373981523457721301218593967440417692678049207833228942000
           - 2452813096069528207645703151222478123259511586701148682951852876484544822947007791153163*y/2425513406294293557160134959260938308852048743758697616707707121968313329157680
           - QQ(34305265428126440542854669008203683099323146152358231964773310260498715579162112959703, 202126117191191129763344579938411525737670728646558134725642260164026110763140),
         y**8 + 693749860237914515552*y**7/67859264524169150569
              + 27761407182086143225024*y**6/610733380717522355121
              + 7785127652157884044288*y**5/67859264524169150569
              + 36567075214771261409792*y**4/203577793572507451707
              + 36336335165196147384320*y**3/203577793572507451707
              + 7452455676042754048000*y**2/67859264524169150569
              + 2593331082514399232000*y/67859264524169150569
              + QQ(390399197427343360000, 67859264524169150569)]

    G = [(3725588592068034903797967297424801242396746870413359539263038139343329273586196480000*x -
          160420835591776763325581422211936558925462474417709511019228211783493866564923546661604487873*y**7 -
          1406108495478033395547109582678806497509499966197028487131115097902188374051595011248311352864*y**6 -
          5241326875850889518164640374668786338033653548841427557880599579174438246266263602956254030352*y**5 -
          10758917262823299139373269714910672770004760114329943852726887632013485035262879510837043892416*y**4 -
          13119383576444715672578819534846747735372132018341964647712009275306635391456880068261130581248*y**3 -
          9491412317016197146080450036267011389660653495578680036574753839055748080962214787557853941760*y**2 -
          3767520915562795326943800040277726397326609797172964377014046018280260848046603967211258368000*y -
          632314652371226552085897259159210286886724229880266931574701654721512325555116066073245696000).set_domain(ZZ),
         (610733380717522355121*y**8 + 6243748742141230639968*y**7 +
          27761407182086143225024*y**6 + 70066148869420956398592*y**5 +
          109701225644313784229376*y**4 + 109009005495588442152960*y**3 +
          67072101084384786432000*y**2 + 23339979742629593088000*y +
          3513592776846090240000).set_domain(ZZ)]

    assert [f.clear_denoms()[1].set_domain(ZZ) for f in F] == G
    assert [f.clear_denoms(convert=True)[1] for f in F] == G

    R, x, y = ring('x y', EX)

    assert (3*x/2 + Rational(9, 4)).clear_denoms() == (4, 6*x + 9)
    assert R(7).clear_denoms() == (1, 7)

    a = EX.to_expr(EX('a'))

    assert (sin(a)/a*y).clear_denoms() == (a, y*sin(a))


def test_PolyElement_terms_gcd():
    R, x = ring('x', ZZ)

    assert R(0).terms_gcd() == ((0,), 0)
    assert (x**2 + 1).terms_gcd() == ((0,), x**2 + 1)
    assert (x**3 + x).terms_gcd() == ((1,), x**2 + 1)
    assert (x**5 + 3*x**4 + x**3 + 4*x**2 +
            2*x).terms_gcd() == ((1,), x**4 + 3*x**3 + x**2 + 4*x + 2)
    assert (x**4 + x**2).terms_gcd() == ((2,), x**2 + 1)

    R, x, y = ring('x y', ZZ)

    assert R(0).terms_gcd() == ((0, 0), 0)
    assert R(1).terms_gcd() == ((0, 0), 1)
    assert (x**3 + x).terms_gcd() == ((1, 0), x**2 + 1)
    assert (x**2*y + 1).terms_gcd() == ((0, 0), x**2*y + 1)
    assert (x**6*y**2 + x**3*y).terms_gcd() == ((3, 1), x**3*y + 1)
    assert (x**3*y + x**2*y**2).terms_gcd() == ((2, 1), x + y)
    assert (x**2 - x*y - 2*y**2).terms_gcd() == ((0, 0), x**2 - x*y - 2*y**2)


def test_PolyElement_max_norm():
    R, x = ring('x', ZZ)

    assert R(0).max_norm() == 0
    assert R(1).max_norm() == 1
    assert (-x**2 + 2*x - 3).max_norm() == 3

    assert (x**3 + 4*x**2 + 2*x + 3).max_norm() == 4

    R, x, y = ring('x y', ZZ)

    assert R(-1).max_norm() == 1
    assert R(+0).max_norm() == 0
    assert R(+1).max_norm() == 1

    assert (x**3 + 4*x**2 + 2*x + 3).max_norm() == 4
    assert (-x**2 + 2*x - 3).max_norm() == 3

    assert (2*x*y - x - 3).max_norm() == 3

    R, x, y, z = ring('x y z', ZZ)

    assert R(0).max_norm() == 0
    assert R(1).max_norm() == 1

    assert f_polys()[0].max_norm() == 6


def test_PolyElement_l1_norm():
    R, x = ring('x', ZZ)

    assert R(-1).l1_norm() == 1
    assert R(+0).l1_norm() == 0
    assert R(+1).l1_norm() == 1

    assert (2*x**3 - 3*x**2 + 1).l1_norm() == 6
    assert (x**3 + 4*x**2 + 2*x + 3).l1_norm() == 10

    R, x, y = ring('x y', ZZ)

    assert R(-1).l1_norm() == 1
    assert R(+0).l1_norm() == 0
    assert R(+1).l1_norm() == 1

    assert (x**3 + 4*x**2 + 2*x + 3).l1_norm() == 10
    assert (-x**2 + 2*x - 3).l1_norm() == 6
    assert (2*x*y - x - 3).l1_norm() == 6

    R, x, y, z = ring('x y z', ZZ)

    assert R(0).l1_norm() == 0
    assert R(1).l1_norm() == 1

    assert f_polys()[0].l1_norm() == 31


def test_PolyElement_diff():
    R, x = ring('x', FF(3))

    f = 2*x**10 + x**9 + 2*x**8 + 2*x**6 + x**5 + 1

    assert f.diff() == 2*x**9 + x**7 + 2*x**4
    assert f.diff(m=2) == x**6 + 2*x**3
    assert f.diff(m=3) == 0

    assert f.diff(m=0) == f
    assert f.diff(m=2) == f.diff().diff()
    assert f.diff(m=3) == f.diff().diff().diff()

    R, x = ring('x', FF(5))

    assert (3*x**2 + 2*x + 4).diff() == x + 2

    R, x = ring('x', FF(11))

    assert R.zero.diff() == 0
    assert R(7).diff() == 0
    assert (7*x + 3).diff() == R(7)
    assert (7*x**2 + 3*x + 1).diff() == 3*x + 3
    assert (x**11 + 1).diff() == 0

    R, x = ring('x', ZZ)

    assert R(0).diff() == 0
    assert R(7).diff() == 0
    assert (2*x + 7).diff() == 2
    assert (x**2 + 2*x + 1).diff() == 2*x + 2
    assert (x**3 + 2*x**2 + 3*x + 4).diff() == 3*x**2 + 4*x + 3
    assert (x**4 - x**3 + 2).diff() == 4*x**3 - 3*x**2
    assert (x**3 + 2*x**2 + 3*x + 4).diff(m=2) == 6*x + 4

    f = 17*x**10 + 34*x**9 + 56*x**8 - 345*x**7 + 23*x**6 + 76*x**5 + 12*x**2 + 3*x + 7

    assert f.diff(m=0) == f
    assert f.diff(m=2) == f.diff().diff()
    assert f.diff(m=3) == f.diff().diff().diff()

    R, x, y = ring('x y', ZZ)

    assert R(0).diff() == 0

    f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    assert f.diff() == y**2 + 2*y + 3
    assert f.diff(m=2) == 0
    assert f.diff(y) == 2*x*y + 2*x + 4*y + 3

    f = x + x**2*y**3

    assert f.diff(x, 0) == f
    assert f.diff(x) == 2*x*y**3 + 1
    assert f.diff(y) == 3*x**2*y**2
    assert f.diff(x, 2) == 2*y**3
    assert f.diff(x, 3) == 0

    R, x, y, z = ring('x y z', ZZ)

    assert R(0).diff() == 0
    assert (y + 2).diff() == 0
    assert x.diff() == 1
    assert (3*x**2 + x).diff() == 6*x + 1

    R, x, y, z, t = ring('x y z t', FF(23))

    f = f_polys()[6].set_ring(R)

    assert f.diff(m=0) == f
    assert f.diff(m=2) == f.diff().diff()
    assert f.diff(m=3) == f.diff().diff().diff()

    R, x, y, z, t = ring('x y z t', ZZ)

    f = f_polys()[6]

    assert f.diff(m=0) == f
    assert f.diff(m=2) == f.diff().diff()
    assert f.diff(m=3) == f.diff().diff().diff()

    assert (f.diff(y, m=2) ==
            -5076*x*y**2 - 282*x*y - 54*y*z**3*t**2 + 54*y*t**2 - 2*z**3*t**2 + 2*t**2)
    assert f.diff(y, m=3) == -10152*x*y - 282*x - 54*z**3*t**2 + 54*t**2
    assert (f.diff(z, m=2) ==
            270*x**3*z*t**2 + 846*x*y*z - 54*y**3*z*t**2 - 6*y**2*z*t**2 +
            90*z**4*t**2 + 24*z**2*t**3 - 18*z*t**2)
    assert (f.diff(z, m=3) ==
            270*x**3*t**2 + 846*x*y - 54*y**3*t**2 - 6*y**2*t**2 +
            360*z**3*t**2 + 48*z*t**3 - 18*t**2)

    R, *X = ring('x:11', QQ)

    f = 288*X[0]**8*X[1]**6*X[4]**3*X[10]**2/5 + 8*X[0]**2*X[2]**3*X[4]**3 + 2*X[0]**2 - 2*X[1]**2

    assert f.diff(X[0]) == 2304*X[0]**7*X[1]**6*X[4]**3*X[10]**2/5 + 16*X[0]*X[2]**3*X[4]**3 + 4*X[0]
    assert f.diff(X[4]) == 864*X[0]**8*X[1]**6*X[4]**2*X[10]**2/5 + 24*X[0]**2*X[2]**3*X[4]**2
    assert f.diff(X[10]) == 576*X[0]**8*X[1]**6*X[4]**3*X[10]/5

    pytest.raises(ValueError, lambda: f.diff(x='spam', m=2))


def test_PolyElement_integrate():
    R, x = ring('x', QQ)

    assert R(0).integrate() == 0
    assert R(0).integrate(m=2) == 0

    assert R(1).integrate() == x
    assert R(1).integrate(m=2) == x**2/2

    assert (x**2 + 2*x + 3).integrate(m=0) == x**2 + 2*x + 3
    assert (x**2 + 2*x + 3).integrate() == x**3/3 + x**2 + 3*x
    assert (x**2 + 2*x + 3).integrate(m=2) == x**4/12 + x**3/3 + 3*x**2/2
    assert (x**2 + 2*x + 3).integrate(m=3) == x**5/60 + x**4/12 + x**3/2

    assert (x**2 + 2*x).integrate() == x**3/3 + x**2
    assert (x**2 + 2*x).integrate(m=2) == x**4/12 + x**3/3

    assert (17*x**29).integrate(m=3) == 17*x**32/29760

    assert (17*x**29 + x**5/2).integrate(m=3) == 17*x**32/29760 + x**8/672

    f = x + 1

    assert f.integrate() == x**2/2 + x
    assert f.integrate(x=0) == x**2/2 + x
    assert f.integrate(x=x) == x**2/2 + x

    R, x, y = ring('x y', QQ)

    assert (x**2 + 2*x + 3).integrate(m=0) == x**2 + 2*x + 3
    assert (x**2 + 2*x + 3).integrate(m=1) == x**3/3 + x**2 + 3*x
    assert (x**2 + 2*x + 3).integrate(m=2) == x**4/12 + x**3/3 + 3*x**2/2
    assert (x**2 + 2*x + 3).integrate(m=3) == x**5/60 + x**4/12 + x**3/2

    assert (x + 2*y).integrate() == x**2/2 + 2*x*y
    assert (x + 2*y).integrate(x=y) == x*y + y**2
    assert (x + 2*y).integrate(m=2) == x**3/6 + x**2*y

    f = x*y + 1

    assert f.integrate(x=x) == x**2*y/2 + x
    assert f.integrate(x=y) == x*y**2/2 + y
    assert f.integrate(x=x, m=2) == x**3*y/6 + x**2/2
    assert f.integrate(x=y, m=2) == x*y**3/6 + y**2/2

    assert f.integrate(x=x).integrate(x=y) == x**2*y**2/4 + x*y
    assert f.integrate(x=y).integrate(x=x) == x**2*y**2/4 + x*y

    R, x, y, z = ring('x y z', QQ)

    assert R(0).integrate() == 0
    assert R(0).integrate(m=2) == 0

    assert R(1).integrate() == x
    assert R(1).integrate(m=2) == x**2/2

    R, x, y, z, t = ring('x y z t', QQ)

    f = f_polys()[6].set_ring(R)

    assert (f.integrate(x=y, m=2) ==
            705*x**4*y**3/2 + 45*x**3*y**2*z**3*t**2/2 - 45*x**3*y**2*t**2/2 -
            141*x*y**6/10 - 47*x*y**5/20 + 47*x*y**3*z**3/2 + 47*x*y**3*z*t/3 -
            9*y**5*z**3*t**2/20 + 9*y**5*t**2/20 - y**4*z**3*t**2/12 + y**4*t**2/12 +
            3*y**2*z**6*t**2/2 + y**2*z**4*t**3 - 3*y**2*z**3*t**2/2 - y**2*z*t**3)
    assert (f.integrate(x=y, m=3) ==
            705*x**4*y**4/8 + 15*x**3*y**3*z**3*t**2/2 - 15*x**3*y**3*t**2/2 -
            141*x*y**7/70 - 47*x*y**6/120 + 47*x*y**4*z**3/8 + 47*x*y**4*z*t/12 -
            3*y**6*z**3*t**2/40 + 3*y**6*t**2/40 - y**5*z**3*t**2/60 + y**5*t**2/60 +
            y**3*z**6*t**2/2 + y**3*z**4*t**3/3 - y**3*z**3*t**2/2 - y**3*z*t**3/3)
    assert (f.integrate(x=z, m=2) ==
            2115*x**4*y*z**2/2 + 9*x**3*z**5*t**2/4 - 45*x**3*z**2*t**2/2 -
            423*x*y**4*z**2/2 - 47*x*y**3*z**2/2 + 141*x*y*z**5/20 + 47*x*y*z**3*t/3 -
            9*y**3*z**5*t**2/20 + 9*y**3*z**2*t**2/2 - y**2*z**5*t**2/20 +
            y**2*z**2*t**2/2 + 3*z**8*t**2/56 + z**6*t**3/15 - 3*z**5*t**2/20 -
            z**3*t**3/3)
    assert (f.integrate(x=z, m=3) ==
            705*x**4*y*z**3/2 + 3*x**3*z**6*t**2/8 - 15*x**3*z**3*t**2/2 -
            141*x*y**4*z**3/2 - 47*x*y**3*z**3/6 + 47*x*y*z**6/40 + 47*x*y*z**4*t/12 -
            3*y**3*z**6*t**2/40 + 3*y**3*z**3*t**2/2 - y**2*z**6*t**2/120 +
            y**2*z**3*t**2/6 + z**9*t**2/168 + z**7*t**3/105 - z**6*t**2/40 -
            z**4*t**3/12)


def test_PolyElement___call__():
    R, x = ring('x', ZZ)

    f = 3*x + 1

    assert f(0) == 1
    assert f(1) == 4

    pytest.raises(ValueError, lambda: f())
    pytest.raises(ValueError, lambda: f(0, 1))

    pytest.raises(CoercionFailed, lambda: f(QQ(1, 7)))

    R, x, y = ring('x y', ZZ)

    f = 3*x + y**2 + 1

    assert f(0, 0) == 1
    assert f(1, 7) == 53

    assert f(0) == (y**2 + 1).drop(x)
    assert f(1) == (y**2 + 4).drop(x)

    pytest.raises(ValueError, lambda: f())
    pytest.raises(ValueError, lambda: f(0, 1, 2))

    pytest.raises(CoercionFailed, lambda: f(1, QQ(1, 7)))
    pytest.raises(CoercionFailed, lambda: f(QQ(1, 7), 1))
    pytest.raises(CoercionFailed, lambda: f(QQ(1, 7), QQ(1, 7)))


def test_PolyElement_eval():
    R, x = ring('x', FF(5))

    assert (3*x**2 + 2*x + 4).eval(x, 2) == 0

    R, x = ring('x', FF(11))

    assert R.zero.eval(x, 4) == 0
    assert R.zero.eval(x, 27) == 0
    assert R(7).eval(x, 4) == 7
    assert R(7).eval(x, 27) == 7

    f = x**8 + 3*x**6 + 2*x**5 + 4*x**4 + 3*x**3 + x**2 + 2*x

    assert f.eval(x, 0) == 0
    assert f.eval(x, 4) == 9
    assert f.eval(x, 27) == 5

    f = 4*x**8 + 4*x**5 + 6*x**4 + x**2 + 3*x + 5

    assert f.eval(x, 0) == 5
    assert f.eval(x, 4) == 3
    assert f.eval(x, 27) == 9

    R, x = ring('x', ZZ)

    assert R(0).eval(a=7) == 0
    assert R(0).eval(a=3) == 0
    assert (x + 2).eval() == 2

    assert (x**2 + x + 1).diff().eval(a=1) == 3

    f = x**2 + 2*x + 3

    assert f.eval(a=7) == 66
    assert f.eval(a=2) == 11

    f = x**3 + 4*x**2 + 2*x + 3
    r = f.eval(x, 0)

    assert r == 3 and not isinstance(r, PolyElement)

    pytest.raises(CoercionFailed, lambda: f.eval(x, QQ(1, 7)))

    R, x, y = ring('x y', ZZ)

    assert R(0).eval(a=3) == 0
    assert R(0).eval([(y, 1)]) == 0
    assert (y + 2).eval() == (y + 2).drop(x)

    assert (3*x*y + 2*x + y + 2).eval(a=3) == (10*y + 8).drop(x)

    f = 2*x*y + 3*x + y + 2

    assert f.eval(a=2) == (5*y + 8).drop(x)
    assert f.eval(x=1, a=2) == (7*x + 4).drop(y)
    assert f.eval([(x, 2), (y, 2)]) == 18

    f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    assert f.diff().eval(a=2) == (y**2 + 2*y + 3).drop(x)
    assert f.diff(x=1).eval(x=1, a=2) == (6*x + 11).drop(y)

    R, x, y, z = ring('x y z', ZZ)

    assert R(0).eval(a=3) == 0
    assert R(0).eval([(z, 1)]) == 0
    assert R(0).eval([(y, 1), (z, 2)]) == 0
    assert R(1).eval(a=3) == 1
    assert (z + 2).eval(a=3) == (z + 2).drop(x)
    assert (3*x*z + 2*x + z + 2).eval(a=3) == (10*z + 8).drop(x)

    f = 45*x**3 - 9*y**3 - y**2 + 3*z**3 + 10*z

    assert f.eval(x=z, a=-2) == (45*x**3 - 9*y**3 - y**2 - 44).drop(z)

    f = (x*y)**3 + 4*(x*y)**2 + 2*x*y + 3
    r = f.eval(x, 0)

    assert r == 3 and isinstance(r, R.drop(x).dtype)

    r = f.eval([(x, 0), (y, 0)])

    assert r == 3 and isinstance(r, R.drop(x, y).dtype)

    r = f.eval(y, 0)

    assert r == 3 and isinstance(r, R.drop(y).dtype)

    r = f.eval([(y, 0), (x, 0)])

    assert r == 3 and isinstance(r, R.drop(y, x).dtype)

    r = f.eval([(x, 0), (y, 0), (z, 0)])

    assert r == 3 and not isinstance(r, PolyElement)

    pytest.raises(CoercionFailed, lambda: f.eval([(x, 1), (y, QQ(1, 7))]))
    pytest.raises(CoercionFailed, lambda: f.eval([(x, QQ(1, 7)), (y, 1)]))
    pytest.raises(CoercionFailed, lambda: f.eval([(x, QQ(1, 7)), (y, QQ(1, 7))]))

    f = f_polys()[0]

    assert f.eval([(x, 1), (y, -17), (z, 8)]) == 84496
    assert f.eval([(y, -17), (z, 8)]) == (-1409*x**2 + 3*x + 85902).drop(-1).drop(-1)
    assert f.eval([(z, 8)]) == (83*x**2*y + 2*x**2 + 3*x + 302*y**2 + 81*y + 1).drop(z)

    f = f_polys()[1]

    assert (f.eval([(y, -17), (z, 8)]) ==
            (-136*x**3 + 15699*x**2 + 9166*x - 27144).drop(-1).drop(-1))

    f = f_polys()[2]

    assert (f.eval([(y, -12), (z, 3)]) ==
            (-1377*x**5 - 702*x**3 - 1224*x**2 - 624).drop(-1).drop(-1))

    f = f_polys()[3]

    assert (f.eval([(y, -12), (z, 3)]) ==
            (144*x**5 + 82*x**4 - 5181*x**3 - 28872*x**2 -
             14868*x - 540).drop(-1).drop(-1))

    f = f_polys()[4]

    assert (f.eval([(y, 25), (z, -1)]) ==
            (152587890625*x**9 + 9765625*x**8 - 59605407714843750*x**7 -
             3839159765625*x**6 - 1562475*x**5 + 9536712644531250*x**4 +
             610349546750*x**3 - 4*x**2 + 24414375000*x + 1562520).drop(-1).drop(-1))

    f = f_polys()[5]

    assert (f.eval([(y, 25), (z, -1)]) ==
            (-x**3 - 78*x**2 - 2028*x - 17576).drop(-1).drop(-1))

    R, x, y, z, t = ring('x y z t', ZZ)

    f = f_polys()[6]

    assert (f.eval(x=y, a=-2) ==
            (-4230*x**4 + 45*x**3*z**3*t**2 - 45*x**3*t**2 - 282*x*z**3 -
             188*x*z*t - 6392*x + 3*z**6*t**2 + 2*z**4*t**3 + 65*z**3*t**2 -
             2*z*t**3 - 68*t**2).drop(y))
    assert (f.eval(x=y, a=7) ==
            (14805*x**4 + 45*x**3*z**3*t**2 - 45*x**3*t**2 + 987*x*z**3 +
             658*x*z*t - 1031744*x + 3*z**6*t**2 + 2*z**4*t**3 -
             3139*z**3*t**2 - 2*z*t**3 + 3136*t**2).drop(y))
    assert (f.diff(x=y, m=2).eval(x=y, a=7) ==
            (-250698*x - 380*z**3*t**2 + 380*t**2).drop(y))

    assert (f.eval(x=z, a=-2) ==
            (2115*x**4*y - 405*x**3*t**2 - 423*x*y**4 - 47*x*y**3 - 188*x*y*t -
             1128*x*y + 81*y**3*t**2 + 9*y**2*t**2 + 36*t**3 + 216*t**2).drop(z))
    assert (f.eval(x=z, a=7) ==
            (2115*x**4*y + 15390*x**3*t**2 - 423*x*y**4 - 47*x*y**3 + 658*x*y*t +
             48363*x*y - 3078*y**3*t**2 - 342*y**2*t**2 + 4788*t**3 +
             351918*t**2).drop(z))
    assert (f.eval([(y, 0), (z, 2), (t, 4)]) ==
            (5040*x**3 + 4480).drop(-1).drop(-1).drop(-1))


def test_PolyElement_compose():
    R, x = ring('x', FF(5))

    assert (3*x**2 + 2*x +
            4).compose(x, 2*x**2 + 2*x + 2) == 2*x**4 + 4*x**3 + 3*x

    R, x = ring('x', FF(11))

    assert R.zero.compose(x, x) == 0
    assert R.one.compose(x, R.zero) == 1
    assert x.compose(x, R.zero) == 0
    assert x.compose(x, x) == x

    f = x**2 + x + 1
    g = x**3 + 2

    assert f.compose(x, g) == x**6 + 5*x**3 + 7

    R, x = ring('x', ZZ)

    assert R(0).compose(x, 0) == 0
    assert R(0).compose(x, 1) == 0
    assert R(0).compose(x, x + 2) == 0
    assert R(0).compose(x, -x) == 0
    assert R(0).compose(x, x + 1) == 0
    assert R(1).compose(x, 0) == 1
    assert R(1).compose(x, -x) == 1
    assert R(1).compose(x, x + 1) == 1

    assert (x**2 + 2*x).compose(x, 0) == 0

    f = x**2 + 2*x + 1

    assert f.compose(x, 0) == 1

    assert f.compose(x, 1) == 4
    assert f.compose(x, 7) == 64

    assert f.compose(x, x - 1) == x**2
    assert f.compose(x, x + 1) == x**2 + 4*x + 4
    assert f.compose(x, x**2 + 2*x + 1) == x**4 + 4*x**3 + 8*x**2 + 8*x + 4

    assert (x**2 + x).compose(x, x - 1) == x**2 - x

    f = x**2 - 2*x + 1

    assert f.compose(x, 2*x) == 4*x**2 - 4*x + 1
    assert f.compose(x, x + 2) == x**2 + 2*x + 1

    f = x**3 + 4*x**2 + 2*x + 3
    r = f.compose(x, 0)

    assert r == 3 and isinstance(r, R.dtype)

    assert f.compose(x, x) == f
    assert f.compose(x, x**2) == x**6 + 4*x**4 + 2*x**2 + 3

    pytest.raises(CoercionFailed, lambda: f.compose(x, QQ(1, 7)))

    f = x**3 + 2*x**2 - 4*x + 2

    assert f.compose(x, -x) == -x**3 + 2*x**2 + 4*x + 2

    f = x**4 + 2*x**3 + 3*x**2 + 4*x + 5

    assert f.compose(x, -x) == x**4 - 2*x**3 + 3*x**2 - 4*x + 5
    assert f.compose(x, -7*x) == 2401*x**4 - 686*x**3 + 147*x**2 - 28*x + 5
    assert f.compose(x, x + 1) == x**4 + 6*x**3 + 15*x**2 + 20*x + 15
    assert f.compose(x, x + 7) == x**4 + 30*x**3 + 339*x**2 + 1712*x + 3267

    f = x**5 + 2*x**4 + 3*x**3 + 4*x**2 + 5*x + 6

    assert f.compose(x, -x) == -x**5 + 2*x**4 - 3*x**3 + 4*x**2 - 5*x + 6

    R, x, y = ring('x y', ZZ)

    assert (x**2 + 2*x).compose(x, 0) == 0

    f = x**2 + 2*x + 1

    assert f.compose(x, 0) == 1

    assert f.compose(x, 1) == 4
    assert f.compose(x, 7) == 64

    assert f.compose(x, x - 1) == x**2
    assert f.compose(x, x + 1) == x**2 + 4*x + 4

    assert f.compose(x, x**2 + 2*x + 1) == x**4 + 4*x**3 + 8*x**2 + 8*x + 4

    assert (x*y + 2*x + y).compose(x, y) == y**2 + 3*y

    R, x, y, z = ring('x y z', ZZ)

    assert R(0).compose(x, 0) == 0
    assert R(0).compose(x, 1) == 0
    assert R(0).compose(x, x + 2) == 0
    assert R(1).compose(x, 0) == 1

    f = x**3 + 4*x**2 + 2*x + 3
    r = f.compose(x, 0)

    assert r == 3 and isinstance(r, R.dtype)

    r = f.compose([(x, 0), (y, 0)])

    assert r == 3 and isinstance(r, R.dtype)

    r = f.compose({x: 0, y: 0})

    assert r == 3 and isinstance(r, R.dtype)

    pytest.raises(ValueError, lambda: f.compose('spam'))

    r = (x**3 + 4*x**2 + 2*x*y*z + 3).compose(x, y*z**2 - 1)
    q = (y*z**2 - 1)**3 + 4*(y*z**2 - 1)**2 + 2*(y*z**2 - 1)*y*z + 3

    assert r == q and isinstance(r, R.dtype)


def test_PolyElement_is_():
    R, x = ring('x', ZZ)

    assert R(0).is_ground is True
    assert R(1).is_ground is True
    assert R(2).is_ground is True

    f = x**16 + x**14 - x**10 - x**8 - x**6 + x**2

    assert f.is_cyclotomic is False
    assert (f + 1).is_cyclotomic

    assert R.is_normal(f) is True

    R, x, y = ring('x y', ZZ)

    assert R(0).is_ground is True
    assert R(1).is_ground is True
    assert R(2).is_ground is True

    assert R.zero.is_homogeneous is True
    assert (x**2 + x*y).is_homogeneous is True
    assert (x**3 + x*y).is_homogeneous is False

    R, x, y, z = ring('x y z', QQ)

    assert R(0).is_ground is True
    assert R(1).is_ground is True
    assert R(2).is_ground is True
    assert (3*y).is_ground is False

    assert (x - x).is_generator is False
    assert (x - x).is_ground
    assert (x - x).is_monomial
    assert (x - x).is_term

    assert (x - x + 1).is_generator is False
    assert (x - x + 1).is_ground
    assert (x - x + 1).is_monomial
    assert (x - x + 1).is_term

    assert x.is_generator
    assert x.is_ground is False
    assert x.is_monomial
    assert x.is_term

    assert (x*y).is_generator is False
    assert (x*y).is_ground is False
    assert (x*y).is_monomial
    assert (x*y).is_term

    assert (3*x).is_generator is False
    assert (3*x).is_ground is False
    assert (3*x).is_monomial is False
    assert (3*x).is_term

    assert (3*x + 1).is_generator is False
    assert (3*x + 1).is_ground is False
    assert (3*x + 1).is_monomial is False
    assert (3*x + 1).is_term is False

    assert (x + y + z + 1).is_linear
    assert (x*y*z + 1).is_linear is False

    assert (x*y + z + 1).is_quadratic
    assert (x*y*z + 1).is_quadratic is False

    pytest.raises(AttributeError, lambda: x.is_cyclotomic)

    R, x, y, z, w, t = ring('x y z w t', ZZ)


def test_PolyElement_drop():
    R, x, y, z = ring('x y z', ZZ)

    assert R(1).drop(0).ring == ZZ.inject('y', 'z')
    assert R(1).drop(0).drop(0).ring == ZZ.inject('z')
    assert isinstance(R(1).drop(0).drop(0).drop(0), R.dtype) is False

    pytest.raises(ValueError, lambda: z.drop(0).drop(0).drop(0))
    pytest.raises(ValueError, lambda: x.drop(0))

    f = z**2*x + 2*z*y + x*z + 1
    R2 = R.eject(z)
    D = R2.domain

    assert f.eject(z) == D.z**2*R2.x + 2*D.z*R2.y + D.z*R2.x + 1

    R12 = R.eject(y, z)
    D = R12.domain

    assert f.eject(y, z) == R12.x*(D.z**2 + D.z) + 2*D.y*D.z + 1

    R02 = R.eject(x, z)
    D = R02.domain

    assert f.eject(x, z) == R02.y*2*D.z + D.x*D.z**2 + D.x*D.z + 1

    R3 = R.drop(y, z)

    assert R3 == ZZ.inject('x')

    pytest.raises(ValueError, lambda: R3.x.eject(R3.x))


def test_PolyElement_decompose():
    R, x = ring('x', ZZ)

    f = x**12 + 20*x**10 + 150*x**8 + 500*x**6 + 625*x**4 - 2*x**3 - 10*x + 9
    g = x**4 - 2*x + 9
    h = x**3 + 5*x

    assert g.compose(x, h) == f
    assert f.decompose() == [g, h]

    f = x**4 - 2*x**3 + x**2
    g = x**2
    h = x**2 - x

    assert g.compose(x, h) == f
    assert f.decompose() == [g, h]

    assert R(1).decompose() == [1]

    assert x.decompose() == [x]
    assert (x**3).decompose() == [x**3]
    assert (x**4).decompose() == [x**2, x**2]
    assert (x**6).decompose() == [x**3, x**2]

    assert (7*x**4 + 1).decompose() == [7*x**2 + 1, x**2]
    assert (4*x**4 + 3*x**2 + 2).decompose() == [4*x**2 + 3*x + 2, x**2]

    f = 2*x**12 + 40*x**10 + 300*x**8 + 1000*x**6 + 1250*x**4 - 4*x**3 - 20*x + 18

    assert f.decompose() == [2*x**4 - 4*x + 18, x**3 + 5*x]

    f = (x**12 + 20*x**10 - 8*x**9 + 150*x**8 - 120*x**7 + 524*x**6 -
         600*x**5 + 865*x**4 - 1034*x**3 + 600*x**2 - 170*x + 29)

    assert f.decompose() == [x**4 - 8*x**3 + 24*x**2 - 34*x + 29, x**3 + 5*x]

    Rt, t = ring('t', ZZ)
    R, x = ring('x', Rt)

    f = ((6*t**2 - 42)*x**4 + (48*t**2 + 96)*x**3 +
         (144*t**2 + 648*t + 288)*x**2 + (624*t**2 + 864*t + 384)*x +
         108*t**3 + 312*t**2 + 432*t + 192)

    assert f.decompose() == [f]


def test_PolyElement_shift():
    _, x = ring('x', ZZ)

    assert (x**2 - 2*x + 1).shift(2) == x**2 + 2*x + 1


def test_PolyElement_sturm():
    F, t = field('t', ZZ)
    _, x = ring('x', F)

    f = 1024/(15625*t**8)*x**5 - 4096/(625*t**8)*x**4 + 32/(15625*t**4)*x**3 - 128/(625*t**4)*x**2 + F(1)/62500*x - F(1)/625

    assert f.sturm() == [
        x**3 - 100*x**2 + t**4/64*x - 25*t**4/16,
        3*x**2 - 200*x + t**4/64,
        (-t**4/96 + F(20000)/9)*x + 25*t**4/18,
        (-9*t**12 - 11520000*t**8 - 3686400000000*t**4)/(576*t**8 - 245760000*t**4 + 26214400000000),
    ]


def test_PolyElement_slice():
    R, x = ring('x', ZZ)

    f = x**3 + 2*x**2 + 3*x + 4

    assert f.slice(0, 0) == f.slice(0, 0, x) == 0
    assert f.slice(0, 1) == f.slice(0, 1, x) == 4
    assert f.slice(0, 2) == f.slice(0, 2, x) == 3*x + 4
    assert f.slice(0, 3) == f.slice(0, 3, x) == 2*x**2 + 3*x + 4

    assert f.slice(0, 4) == f.slice(0, 4, x) == f
    assert f.slice(0, 9) == f

    assert f.slice(1, 0) == 0
    assert f.slice(1, 1) == 0
    assert f.slice(1, 2) == 3*x
    assert f.slice(1, 3) == 2*x**2 + 3*x
    assert f.slice(1, 4) == x**3 + 2*x**2 + 3*x

    assert (x + 2).slice(0, 3) == x + 2

    R, x, y = ring('x y', ZZ)

    f = x + 2*y**2 + 3*y + 4

    assert f.slice(1, 2) == f
    assert f.slice(2, 1) == 2*y**2 + 3*y + 5


def test_cache():
    R1 = QQ.inject(-sqrt(2))
    R2 = QQ.inject(-2*sqrt(2))

    assert R1 != R2


def test__gf_random():
    R, x = ring('x', FF(11))

    for n in range(8):
        f = R._gf_random(n, irreducible=True)

        assert f.monic() == f
        assert f.is_irreducible is True


def test_PolyElement__setitem__():
    R, x = ring('x', ZZ)

    p = x + 1
    p[(2,)] = 3

    assert p == 3*x**2 + x + 1

    p_set = {p}
    p_set.pop()

    with pytest.raises(RuntimeError):
        p[(3,)] = 1
