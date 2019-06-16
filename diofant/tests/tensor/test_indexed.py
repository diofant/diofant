import pytest

from diofant import Idx, Indexed, IndexedBase
from diofant.abc import a, b, c, x
from diofant.core import Basic, Dict, Symbol, Tuple, oo, symbols
from diofant.core.compatibility import iterable
from diofant.tensor.indexed import IndexException


__all__ = ()


def test_Idx_construction():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i) != Idx(i, 1)
    assert Idx(i, a) == Idx(i, (0, a - 1))
    assert Idx(i, oo) == Idx(i, (0, oo))

    pytest.raises(TypeError, lambda: Idx(x))
    pytest.raises(TypeError, lambda: Idx(0.5))
    pytest.raises(TypeError, lambda: Idx(i, x))
    pytest.raises(TypeError, lambda: Idx(i, 0.5))
    pytest.raises(TypeError, lambda: Idx(i, (x, 5)))
    pytest.raises(TypeError, lambda: Idx(i, (2, x)))
    pytest.raises(TypeError, lambda: Idx(i, (2, 3.5)))
    pytest.raises(ValueError, lambda: Idx(i, (1, 2, 3)))
    pytest.raises(TypeError, lambda: Idx(i, Basic()))


def test_Idx_properties():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i).is_integer


def test_Idx_bounds():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i).lower is None
    assert Idx(i).upper is None
    assert Idx(i, a).lower == 0
    assert Idx(i, a).upper == a - 1
    assert Idx(i, 5).lower == 0
    assert Idx(i, 5).upper == 4
    assert Idx(i, oo).lower == 0
    assert Idx(i, oo).upper == oo
    assert Idx(i, (a, b)).lower == a
    assert Idx(i, (a, b)).upper == b
    assert Idx(i, (1, 5)).lower == 1
    assert Idx(i, (1, 5)).upper == 5
    assert Idx(i, (-oo, oo)).lower == -oo
    assert Idx(i, (-oo, oo)).upper == oo


def test_Idx_fixed_bounds():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(2).lower is None
    assert Idx(2).upper is None
    assert Idx(2, a).lower == 0
    assert Idx(2, a).upper == a - 1
    assert Idx(2, 5).lower == 0
    assert Idx(2, 5).upper == 4
    assert Idx(2, oo).lower == 0
    assert Idx(2, oo).upper == oo
    assert Idx(2, (a, b)).lower == a
    assert Idx(2, (a, b)).upper == b
    assert Idx(2, (1, 5)).lower == 1
    assert Idx(2, (1, 5)).upper == 5
    assert Idx(2, (-oo, oo)).lower == -oo
    assert Idx(2, (-oo, oo)).upper == oo


def test_Idx_func_args():
    i, a, b = symbols('i a b', integer=True)
    ii = Idx(i)
    assert ii.func(*ii.args) == ii
    ii = Idx(i, a)
    assert ii.func(*ii.args) == ii
    ii = Idx(i, (a, b))
    assert ii.func(*ii.args) == ii


def test_Idx_subs():
    i, a, b = symbols('i a b', integer=True)
    assert Idx(i, a).subs({a: b}) == Idx(i, b)
    assert Idx(i, a).subs({i: b}) == Idx(b, a)

    assert Idx(i).subs({i: 2}) == Idx(2)
    assert Idx(i, a).subs({a: 2}) == Idx(i, 2)
    assert Idx(i, (a, b)).subs({i: 2}) == Idx(2, (a, b))


def test_IndexedBase_sugar():
    i, j = symbols('i j', integer=True)
    A1 = Indexed(a, i, j)
    A2 = IndexedBase(a)
    assert A1 == A2[i, j]
    assert A1 == A2[(i, j)]
    assert A1 == A2[[i, j]]
    assert A1 == A2[Tuple(i, j)]
    assert all(a.is_Integer for a in A2[1, 0].args[1:])


def test_IndexedBase_subs():
    i, j, k = symbols('i j k', integer=True)
    A = IndexedBase(a)
    B = IndexedBase(b)
    C = IndexedBase(c)
    assert A[i] == B[i].subs({b: a})
    assert isinstance(C[1].subs({C: Dict({1: 2})}), type(A[1]))


def test_IndexedBase_shape():
    i, j, m, n = symbols('i j m n', integer=True)
    a = IndexedBase('a', shape=(m, m))
    b = IndexedBase('a', shape=(m, n))
    assert b.shape == Tuple(m, n)
    assert a[i, j] != b[i, j]
    assert a[i, j] == b[i, j].subs({n: m})
    assert b.func(*b.args) == b
    assert b[i, j].func(*b[i, j].args) == b[i, j]
    pytest.raises(IndexException, lambda: b[i])
    pytest.raises(IndexException, lambda: b[i, i, j])

    F = IndexedBase("F", shape=m)
    assert F.shape == Tuple(m)
    assert F[i].subs({i: j}) == F[j]
    pytest.raises(IndexException, lambda: F[i, j])


def test_Indexed_constructor():
    i, j = symbols('i j', integer=True)
    A = Indexed('A', i, j)
    assert A == Indexed(Symbol('A'), i, j)
    assert A == Indexed(IndexedBase('A'), i, j)
    pytest.raises(TypeError, lambda: Indexed(A, i, j))
    pytest.raises(IndexException, lambda: Indexed("A"))


def test_Indexed_func_args():
    i, j = symbols('i j', integer=True)
    A = Indexed(a, i, j)
    assert A == A.func(*A.args)


def test_Indexed_subs():
    i, j, k = symbols('i j k', integer=True)
    A = IndexedBase(a)
    B = IndexedBase(b)
    assert A[i, j] == B[i, j].subs({b: a})
    assert A[i, j] == A[i, k].subs({k: j})


def test_Indexed_properties():
    i, j = symbols('i j', integer=True)
    A = Indexed('A', i, j)
    assert A.rank == 2
    assert A.indices == (i, j)
    assert A.base == IndexedBase('A')
    assert A.ranges == [None, None]
    pytest.raises(IndexException, lambda: A.shape)

    n, m = symbols('n m', integer=True)
    assert Indexed('A', Idx(
        i, m), Idx(j, n)).ranges == [Tuple(0, m - 1), Tuple(0, n - 1)]
    assert Indexed('A', Idx(i, m), Idx(j, n)).shape == Tuple(m, n)
    pytest.raises(IndexException, lambda: Indexed("A", Idx(i, m), Idx(j)).shape)


def test_Indexed_shape_precedence():
    i, j = symbols('i j', integer=True)
    o, p = symbols('o p', integer=True)
    n, m = symbols('n m', integer=True)
    a = IndexedBase('a', shape=(o, p))
    assert a.shape == Tuple(o, p)
    assert Indexed(
        a, Idx(i, m), Idx(j, n)).ranges == [Tuple(0, m - 1), Tuple(0, n - 1)]
    assert Indexed(a, Idx(i, m), Idx(j, n)).shape == Tuple(o, p)
    assert Indexed(
        a, Idx(i, m), Idx(j)).ranges == [Tuple(0, m - 1), Tuple(None, None)]
    assert Indexed(a, Idx(i, m), Idx(j)).shape == Tuple(o, p)


def test_complex_indices():
    i, j = symbols('i j', integer=True)
    A = Indexed('A', i, i + j)
    assert A.rank == 2
    assert A.indices == (i, i + j)


def test_not_interable():
    i, j = symbols('i j', integer=True)
    A = Indexed('A', i, i + j)
    assert not iterable(A)


def test_Indexed_coeff():
    N = Symbol('N', integer=True)
    len_y = N
    i = Idx('i', len_y-1)
    y = IndexedBase('y', shape=(len_y,))
    a = (1/y[i+1]*y[i]).coeff(y[i])
    b = (y[i]/y[i+1]).coeff(y[i])
    assert a == b


def test_indexed_is_constant():
    A = IndexedBase("A")
    i, j, k = symbols("i,j,k")
    assert not A[i].is_constant()
    assert A[i].is_constant(j)
    assert not A[1 + 2*i, k].is_constant()
    assert not A[1 + 2*i, k].is_constant(i)
    assert A[1 + 2*i, k].is_constant(j)
    assert not A[1 + 2*i, k].is_constant(k)
