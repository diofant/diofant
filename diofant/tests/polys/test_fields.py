"""Test sparse rational functions."""

import pytest

from diofant import (CC, QQ, ZZ, CoercionFailed, I, Rational, field, ring,
                     sqrt, symbols)
from diofant.polys.fields import FracElement


__all__ = ()


def test_FractionField___init__():
    F1 = ZZ.inject('x', 'y').field
    F2 = ZZ.inject('x', 'y').field
    F3 = ZZ.inject('x', 'y', 'z').field

    assert F1.x == F1.gens[0]
    assert F1.y == F1.gens[1]
    assert F1.x == F2.x
    assert F1.y == F2.y
    assert F1.x != F3.x
    assert F1.y != F3.y

    F4 = ZZ.inject('gens').field
    assert type(F4.gens) is tuple


def test_FractionField___hash__():
    F, x, y, z = field('x y z', QQ)
    assert hash(F)


def test_FractionField___eq__():
    assert field('x y z', QQ)[0] == field('x y z', QQ)[0]
    assert field('x y z', QQ)[0] is field('x y z', QQ)[0]

    assert field('x y z', QQ)[0] != field('x y z', ZZ)[0]
    assert field('x y z', QQ)[0] is not field('x y z', ZZ)[0]

    assert field('x y z', ZZ)[0] != field('x y z', QQ)[0]
    assert field('x y z', ZZ)[0] is not field('x y z', QQ)[0]

    assert field('x y z', QQ)[0] != field('x y', QQ)[0]
    assert field('x y z', QQ)[0] is not field('x y', QQ)[0]

    assert field('x y', QQ)[0] != field('x y z', QQ)[0]
    assert field('x y', QQ)[0] is not field('x y z', QQ)[0]


def test_FractionField_methods():
    F = ZZ.inject('x').field

    assert F.domain_new(2) == ZZ(2)

    x = symbols('x')
    assert F(x**2 + x) == F.x**2 + F.x


def test_FracElement___hash__():
    F, x, y, z = field('x y z', QQ)
    assert hash(x*y/z)


def test_FracElement_copy():
    F, x, y, z = field('x y z', ZZ)

    f = x*y/3*z
    g = f.copy()

    assert f == g
    g.numerator[(1, 1, 1)] = 7
    assert f != g


def test_FracElement_as_expr():
    F, x, y, z = field('x y z', ZZ)
    f = (3*x**2*y - x*y*z)/(7*z**3 + 1)

    X, Y, Z = F.symbols
    g = (3*X**2*Y - X*Y*Z)/(7*Z**3 + 1)

    assert f != g
    assert F.to_expr(f) == g

    X, Y, Z = symbols('x y z')
    g = (3*X**2*Y - X*Y*Z)/(7*Z**3 + 1)

    assert f != g


def test_FracElement_from_expr():
    x, y, z = symbols('x y z')
    F, X, Y, Z = field((x, y, z), ZZ)

    f = F.convert(1)
    assert f == 1 and isinstance(f, F.dtype)

    f = F.convert(Rational(3, 7))
    assert f == F(3)/7 and isinstance(f, F.dtype)

    f = F.convert(x)
    assert f == X and isinstance(f, F.dtype)

    f = F.convert(Rational(3, 7)*x)
    assert f == 3*X/7 and isinstance(f, F.dtype)

    f = F.convert(1/x)
    assert f == 1/X and isinstance(f, F.dtype)

    f = F.convert(x*y*z)
    assert f == X*Y*Z and isinstance(f, F.dtype)

    f = F.convert(x*y/z)
    assert f == X*Y/Z and isinstance(f, F.dtype)

    f = F.convert(x*y*z + x*y + x)
    assert f == X*Y*Z + X*Y + X and isinstance(f, F.dtype)

    f = F.convert((x*y*z + x*y + x)/(x*y + 7))
    assert f == (X*Y*Z + X*Y + X)/(X*Y + 7) and isinstance(f, F.dtype)

    f = F.convert(x**3*y*z + x**2*y**7 + 1)
    assert f == X**3*Y*Z + X**2*Y**7 + 1 and isinstance(f, F.dtype)

    pytest.raises(CoercionFailed, lambda: F.convert(2**x))
    pytest.raises(CoercionFailed, lambda: F.convert(7*x + sqrt(2)))

    F,  X, Y = field((2**x, y), ZZ)
    f = F.convert(2**(2*x) + 1)
    assert f == X**2 + 1

    # issue sympy/sympy#20985
    F, X = field(x, CC)

    assert F.convert(I/x) == F.convert(CC(0, 1))/X


def test_FracElement_to_poly():
    F, x, y = field('x y', ZZ)
    pytest.raises(ValueError, lambda: (x/y).to_poly())


def test_FracElement__pos_neg__():
    F,  x, y = field('x y', QQ)

    f = (7*x - 9)/y
    g = (-7*x + 9)/y

    assert +f == f
    assert +g == g
    assert -f == g
    assert -g == f


def test_FracElement___add__():
    F,  x, y = field('x y', QQ)

    f, g = 1/x, 1/y
    assert f + g == g + f == (x + y)/(x*y)

    z = symbols('z')
    pytest.raises(TypeError, lambda: x + z)

    assert x + F.ring.gens[0] == F.ring.gens[0] + x == 2*x

    F,  x, y = field('x y', ZZ)
    assert x + 3 == 3 + x
    assert x + QQ(3, 7) == QQ(3, 7) + x == (7*x + 3)/7

    Fuv,  u, v = field('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Fuv)

    f = (u*v + x)/(y + u*v)
    assert dict(f.numerator) == {(1, 0, 0, 0): 1, (0, 0, 0, 0): u*v}
    assert dict(f.denominator) == {(0, 1, 0, 0): 1, (0, 0, 0, 0): u*v}

    Ruv,  u, v = ring('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Ruv)

    f = (u*v + x)/(y + u*v)
    assert dict(f.numerator) == {(1, 0, 0, 0): 1, (0, 0, 0, 0): u*v}
    assert dict(f.denominator) == {(0, 1, 0, 0): 1, (0, 0, 0, 0): u*v}


def test_FracElement___sub__():
    F,  x, y = field('x y', QQ)

    f, g = 1/x, 1/y
    assert f - g == (-x + y)/(x*y)

    assert x - F.ring.gens[0] == F.ring.gens[0] - x == 0

    F,  x, y = field('x y', ZZ)
    assert x - 3 == -(3 - x)
    assert x - QQ(3, 7) == -(QQ(3, 7) - x) == (7*x - 3)/7

    Fuv,  u, v = field('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Fuv)

    f = (u*v - x)/(y - u*v)
    assert dict(f.numerator) == {(1, 0, 0, 0): -1, (0, 0, 0, 0): u*v}
    assert dict(f.denominator) == {(0, 1, 0, 0): 1, (0, 0, 0, 0): -u*v}

    Ruv,  u, v = ring('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Ruv)

    f = (u*v - x)/(y - u*v)
    assert dict(f.numerator) == {(1, 0, 0, 0): -1, (0, 0, 0, 0): u*v}
    assert dict(f.denominator) == {(0, 1, 0, 0): 1, (0, 0, 0, 0): -u*v}

    Fuv,  u, v = field('u v', ZZ)
    Rxyz,  x, y, z = ring('x y z', Fuv)

    f = u - x
    assert dict(f) == {(0, 0, 0): u, (1, 0, 0): -Fuv.one}


def test_FracElement___mul__():
    F,  x, y = field('x y', QQ)

    f, g = 1/x, 1/y
    assert f*g == g*f == 1/(x*y)

    assert x*F.ring.gens[0] == F.ring.gens[0]*x == x**2

    F,  x, y = field('x y', ZZ)
    assert x*3 == 3*x
    assert x*QQ(3, 7) == QQ(3, 7)*x == 3*x/7

    Fuv,  u, v = field('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Fuv)

    f = ((u + 1)*x*y + 1)/((v - 1)*z - t*u*v - 1)
    assert dict(f.numerator) == {(1, 1, 0, 0): u + 1, (0, 0, 0, 0): 1}
    assert dict(f.denominator) == {(0, 0, 1, 0): v - 1, (0, 0, 0, 1): -u*v, (0, 0, 0, 0): -1}

    Ruv,  u, v = ring('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Ruv)

    f = ((u + 1)*x*y + 1)/((v - 1)*z - t*u*v - 1)
    assert dict(f.numerator) == {(1, 1, 0, 0): u + 1, (0, 0, 0, 0): 1}
    assert dict(f.denominator) == {(0, 0, 1, 0): v - 1, (0, 0, 0, 1): -u*v, (0, 0, 0, 0): -1}


def test_FracElement___truediv__():
    F,  x, y = field('x y', QQ)

    f, g = 1/x, 1/y
    assert f/g == y/x

    assert x/F.ring.gens[0] == F.ring.gens[0]/x == 1

    F,  x, y = field('x y', ZZ)
    assert x*3 == 3*x
    assert x/QQ(3, 7) == (QQ(3, 7)/x)**-1 == 7*x/3

    pytest.raises(ZeroDivisionError, lambda: x/0)
    pytest.raises(ZeroDivisionError, lambda: 1/(x - x))
    pytest.raises(ZeroDivisionError, lambda: x/(x - x))

    Fuv,  u, v = field('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Fuv)

    f = (u*v)/(x*y)
    assert dict(f.numerator) == {(0, 0, 0, 0): u*v}
    assert dict(f.denominator) == {(1, 1, 0, 0): 1}

    g = (x*y)/(u*v)
    assert dict(g.numerator) == {(1, 1, 0, 0): 1}
    assert dict(g.denominator) == {(0, 0, 0, 0): u*v}

    Ruv,  u, v = ring('u v', ZZ)
    Fxyzt,  x, y, z, t = field('x y z t', Ruv)

    f = (u*v)/(x*y)
    assert dict(f.numerator) == {(0, 0, 0, 0): u*v}
    assert dict(f.denominator) == {(1, 1, 0, 0): 1}

    g = (x*y)/(u*v)
    assert dict(g.numerator) == {(1, 1, 0, 0): 1}
    assert dict(g.denominator) == {(0, 0, 0, 0): u*v}

    Fuv,  u, v = field('u v', ZZ)
    Rxyz,  x, y, z = ring('x y z', Fuv)

    pytest.raises(TypeError, lambda: u/x)


def test_FracElement___pow__():
    F,  x, y = field('x y', QQ)

    f, g = 1/x, 1/y

    assert f**3 == 1/x**3
    assert g**3 == 1/y**3

    assert (f*g)**3 == 1/(x**3*y**3)
    assert (f*g)**-3 == (x*y)**3

    pytest.raises(ZeroDivisionError, lambda: (x - x)**-3)


def test_FracElement_diff():
    F,  x, y, z = field('x y z', ZZ)

    assert ((x**2 + y)/(z + 1)).diff(x) == 2*x/(z + 1)

    F,  x, y = field('x y', QQ.algebraic_field(I))

    assert ((x - y)/x).diff(x) == y/x**2


def test_FracElement___call__():
    F,  x, y, z = field('x y z', ZZ)
    f = (x**2 + 3*y)/z

    pytest.raises(ValueError, lambda: f(1, 1, 1, 1))

    r = f(1, 1, 1)
    assert r == 4 and not isinstance(r, FracElement)
    pytest.raises(ZeroDivisionError, lambda: f(1, 1, 0))

    Fz = ZZ.inject('z').field
    assert f(1, 1) == 4/Fz.z


def test_FracElement_eval():
    F,  x, y, z = field('x y z', ZZ)
    Fyz = field('y z', ZZ)[0]
    f = (x**2 + 3*y)/z

    assert f.eval(x, 0) == 3*Fyz.y/Fyz.z
    pytest.raises(ZeroDivisionError, lambda: f.eval(z, 0))


def test_FracElement_compose():
    F, x, y, z = field('x y z', QQ)

    f = x**3

    assert f.compose([(x, x/(y + z)),
                      (y, z/x)]) == x**3/(y**3 + 3*y**2*z + 3*y*z**2 + z**3)

    # issue sympy/sympy#20484
    assert f.compose(x, x/(y + z)) == x**3/(y**3 + 3*y**2*z + 3*y*z**2 + z**3)


def test_cache():
    F1 = QQ.frac_field(-sqrt(2))
    F2 = QQ.frac_field(-2*sqrt(2))

    assert F1 != F2
