import pytest

from diofant import ImmutableMatrix as Matrix
from diofant import cos, pi, simplify, sin, symbols
from diofant.vector import (BaseDyadic, CoordSysCartesian, Dyadic, DyadicAdd,
                            DyadicMul, DyadicZero, Vector, express)


__all__ = ()

A = CoordSysCartesian('A')


def test_dyadic():
    a, b = symbols('a, b')
    assert Dyadic.zero != 0
    assert isinstance(Dyadic.zero, DyadicZero)
    assert BaseDyadic(A.i, A.j) != BaseDyadic(A.j, A.i)
    assert (BaseDyadic(Vector.zero, A.i) ==
            BaseDyadic(A.i, Vector.zero) == Dyadic.zero)
    pytest.raises(TypeError, lambda: BaseDyadic(A.i, a))
    pytest.raises(TypeError, lambda: BaseDyadic(a, A.j))

    d1 = A.i | A.i
    d2 = A.j | A.j
    d3 = A.i | A.j

    assert isinstance(d1, BaseDyadic)
    d_mul = a*d1
    assert isinstance(d_mul, DyadicMul)
    assert d_mul.base_dyadic == d1
    assert d_mul.measure_number == a
    assert isinstance(a*d1 + b*d3, DyadicAdd)
    assert d1 == A.i.outer(A.i)
    assert d3 == A.i.outer(A.j)
    v1 = a*A.i - A.k
    v2 = A.i + b*A.j
    assert v1 | v2 == v1.outer(v2) == (a * (A.i | A.i) + (a*b) * (A.i | A.j)
                                       - (A.k | A.i) - b * (A.k | A.j))
    pytest.raises(TypeError, lambda: A.i.outer(A.x))

    assert d1 * 0 == Dyadic.zero
    assert d1 != Dyadic.zero
    assert d1 * 2 == 2 * (A.i | A.i)
    assert d1 / 2. == 0.5 * d1
    pytest.raises(TypeError, lambda: d1/d2)

    assert d1.dot(0 * d1) == Vector.zero
    assert d1 & d2 == Dyadic.zero
    assert d1.dot(A.i) == A.i == d1 & A.i
    pytest.raises(TypeError, lambda: d1.dot(a))

    assert d1.cross(Vector.zero) == Dyadic.zero
    assert d1.cross(A.i) == Dyadic.zero
    assert d1 ^ A.j == d1.cross(A.j)
    assert d1.cross(A.k) == - A.i | A.j
    assert d2.cross(A.i) == - A.j | A.k == d2 ^ A.i
    pytest.raises(TypeError, lambda: d1.cross(a))

    assert A.i ^ d1 == Dyadic.zero
    assert A.j.cross(d1) == - A.k | A.i == A.j ^ d1
    assert Vector.zero.cross(d1) == Dyadic.zero
    assert A.k ^ d1 == A.j | A.i
    assert A.i.dot(d1) == A.i & d1 == A.i
    assert A.j.dot(d1) == Vector.zero
    assert Vector.zero.dot(d1) == Vector.zero
    assert A.j & d2 == A.j

    assert d1.dot(d3) == d1 & d3 == A.i | A.j == d3
    assert d3 & d1 == Dyadic.zero

    q = symbols('q')
    B = A.orient_new_axis('B', q, A.k)
    assert express(d1, B) == express(d1, B, B)
    assert express(d1, B) == ((cos(q)**2) * (B.i | B.i) + (-sin(q) * cos(q)) *
                              (B.i | B.j) + (-sin(q) * cos(q)) * (B.j | B.i) + (sin(q)**2) *
                              (B.j | B.j))
    assert express(d1, B, A) == (cos(q)) * (B.i | A.i) + (-sin(q)) * (B.j | A.i)
    assert express(d1, A, B) == (cos(q)) * (A.i | B.i) + (-sin(q)) * (A.i | B.j)
    assert d1.to_matrix(A) == Matrix([[1, 0, 0], [0, 0, 0], [0, 0, 0]])
    assert d1.to_matrix(A, B) == Matrix([[cos(q), -sin(q), 0],
                                         [0, 0, 0],
                                         [0, 0, 0]])
    assert d3.to_matrix(A) == Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    a, b, c, d, e, f = symbols('a, b, c, d, e, f')
    v1 = a * A.i + b * A.j + c * A.k
    v2 = d * A.i + e * A.j + f * A.k
    d4 = v1.outer(v2)
    assert d4.to_matrix(A) == Matrix([[a * d, a * e, a * f],
                                      [b * d, b * e, b * f],
                                      [c * d, c * e, c * f]])
    d5 = v1.outer(v1)
    C = A.orient_new_axis('C', q, A.i)
    for expected, actual in zip(C.rotation_matrix(A) * d5.to_matrix(A) *
                                C.rotation_matrix(A).T, d5.to_matrix(C)):
        assert (expected - actual).simplify() == 0


def test_dyadic_simplify():
    x, y, k, m, s, A = symbols('x, y, k, m, s, A')
    N = CoordSysCartesian('N')

    dy = N.i | N.i
    test1 = (1 / x + 1 / y) * dy
    assert (N.i & test1 & N.i) != (x + y) / (x * y)
    test1 = test1.simplify()
    assert test1.simplify() == simplify(test1)
    assert (N.i & test1 & N.i) == (x + y) / (x * y)

    test2 = (A**2 * s**4 / (4 * pi * k * m**3)) * dy
    test2 = test2.simplify()
    assert (N.i & test2 & N.i) == (A**2 * s**4 / (4 * pi * k * m**3))

    test3 = ((4 + 4 * x - 2 * (2 + 2 * x)) / (2 + 2 * x)) * dy
    test3 = test3.simplify()
    assert (N.i & test3 & N.i) == 0

    test4 = ((-4 * x * y**2 - 2 * y**3 - 2 * x**2 * y) / (x + y)**2) * dy
    test4 = test4.simplify()
    assert (N.i & test4 & N.i) == -2 * y
