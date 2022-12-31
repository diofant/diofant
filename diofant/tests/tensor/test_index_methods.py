import pytest

from diofant import (Function, Idx, IndexedBase, Pow, exp,
                     get_contraction_structure, get_indices)
from diofant.abc import a, x, y
from diofant.tensor.index_methods import IndexConformanceExceptionError


__all__ = ()


def test_trivial_indices():
    assert get_indices(x) == (set(), {})
    assert get_indices(x*y) == (set(), {})
    assert get_indices(x + y) == (set(), {})
    assert get_indices(x**y) == (set(), {})
    assert get_indices(None) == (set(), {})


def test_get_indices_Indexed():
    x = IndexedBase('x')
    i, j = Idx('i'), Idx('j')
    assert get_indices(x[i, j]) == ({i, j}, {})
    assert get_indices(x[j, i]) == ({j, i}, {})


def test_get_indices_Idx():
    f = Function('f')
    i, j = Idx('i'), Idx('j')
    assert get_indices(f(i)*j) == ({i, j}, {})
    assert get_indices(f(j, i)) == ({j, i}, {})
    assert get_indices(f(i)*i) == (set(), {})


def test_get_indices_mul():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    assert get_indices(x[j]*y[i]) == ({i, j}, {})
    assert get_indices(x[i]*y[j]) == ({i, j}, {})


def test_get_indices_exceptions():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    pytest.raises(IndexConformanceExceptionError, lambda: get_indices(x[i] + y[j]))


def test_scalar_broadcast():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i = Idx('i')
    assert get_indices(x[i] + y[i, i]) == ({i}, {})


def test_get_indices_add():
    x = IndexedBase('x')
    y = IndexedBase('y')
    A = IndexedBase('A')
    i, j, k = Idx('i'), Idx('j'), Idx('k')
    assert get_indices(x[i] + 2*y[i]) == ({i, }, {})
    assert get_indices(y[i] + 2*A[i, j]*x[j]) == ({i, }, {})
    assert get_indices(y[i] + 2*(x[i] + A[i, j]*x[j])) == ({i, }, {})
    assert get_indices(y[i] + x[i]*(A[j, j] + 1)) == ({i, }, {})
    assert get_indices(
        y[i] + x[i]*x[j]*(y[j] + A[j, k]*x[k])) == ({i, }, {})


def test_get_indices_Pow():
    x = IndexedBase('x')
    y = IndexedBase('y')
    A = IndexedBase('A')
    i, j, k = Idx('i'), Idx('j'), Idx('k')
    assert get_indices(Pow(x[i], y[j])) == ({i, j}, {})
    assert get_indices(Pow(x[i, k], y[j, k])) == ({i, j, k}, {})
    assert get_indices(Pow(A[i, k], y[k] + A[k, j]*x[j])) == ({i, k}, {})
    assert get_indices(Pow(2, x[i])) == get_indices(exp(x[i]))

    # test of a design decision, this may change:
    assert get_indices(Pow(x[i], 2)) == ({i, }, {})


def test_get_contraction_structure_basic():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    f = Function('f')
    assert get_contraction_structure(x[i]*y[j]) == {None: {x[i]*y[j]}}
    assert get_contraction_structure(x[i] + y[j]) == {None: {x[i], y[j]}}
    assert get_contraction_structure(x[i]*y[i]) == {(i,): {x[i]*y[i]}}
    assert get_contraction_structure(
        1 + x[i]*y[i]) == {None: {1}, (i,): {x[i]*y[i]}}
    assert get_contraction_structure(x[i]**y[i]) == {None: {x[i]**y[i]}}
    assert (get_contraction_structure(f(x[i, i])) ==
            {None: {f(x[i, i])}, f(x[i, i]): [{(i,): {x[i, i]}}]})


def test_get_contraction_structure_complex():
    x = IndexedBase('x')
    y = IndexedBase('y')
    A = IndexedBase('A')
    i, j, k = Idx('i'), Idx('j'), Idx('k')
    expr1 = y[i] + A[i, j]*x[j]
    d1 = {None: {y[i]}, (j,): {A[i, j]*x[j]}}
    assert get_contraction_structure(expr1) == d1
    expr2 = expr1*A[k, i] + x[k]
    d2 = {None: {x[k]}, (i,): {expr1*A[k, i]}, expr1*A[k, i]: [d1]}
    assert get_contraction_structure(expr2) == d2


def test_contraction_structure_simple_Pow():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    ii_jj = x[i, i]**y[j, j]
    assert get_contraction_structure(ii_jj) == {
        None: {ii_jj},
        ii_jj: [
            {(i,): {x[i, i]}},
            {(j,): {y[j, j]}}
        ]
    }


def test_contraction_structure_Mul_and_Pow():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j, k = Idx('i'), Idx('j'), Idx('k')

    i_ji = x[i]**(y[j]*x[i])
    assert get_contraction_structure(i_ji) == {None: {i_ji}}
    ij_i = (x[i]*y[j])**(y[i])
    assert get_contraction_structure(ij_i) == {None: {ij_i}}
    j_ij_i = x[j]*(x[i]*y[j])**(y[i])
    assert get_contraction_structure(j_ij_i) == {(j,): {j_ij_i}}
    j_i_ji = x[j]*x[i]**(y[j]*x[i])
    assert get_contraction_structure(j_i_ji) == {(j,): {j_i_ji}}
    ij_exp_kki = x[i]*y[j]*exp(y[i]*y[k, k])
    result = get_contraction_structure(ij_exp_kki)
    expected = {
        (i,): {ij_exp_kki},
        ij_exp_kki: [{
                     None: {exp(y[i]*y[k, k])},
                     exp(y[i]*y[k, k]): [{
                         None: {y[i]*y[k, k]},
                         y[i]*y[k, k]: [{(k,): {y[k, k]}}]
                     }]}
                     ]
    }
    assert result == expected


def test_contraction_structure_Add_in_Pow():
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')
    s_ii_jj_s = (1 + x[i, i])**(1 + y[j, j])
    expected = {
        None: {s_ii_jj_s},
        s_ii_jj_s: [
            {None: {1}, (i,): {x[i, i]}},
            {None: {1}, (j,): {y[j, j]}}
        ]
    }
    result = get_contraction_structure(s_ii_jj_s)
    assert result == expected


def test_contraction_structure_Pow_in_Pow():
    x = IndexedBase('x')
    y = IndexedBase('y')
    z = IndexedBase('z')
    i, j, k = Idx('i'), Idx('j'), Idx('k')
    ii_jj_kk = x[i, i]**y[j, j]**z[k, k]
    expected = {
        None: {ii_jj_kk},
        ii_jj_kk: [
            {(i,): {x[i, i]}},
            {
                None: {y[j, j]**z[k, k]},
                y[j, j]**z[k, k]: [
                    {(j,): {y[j, j]}},
                    {(k,): {z[k, k]}}
                ]
            }
        ]
    }
    assert get_contraction_structure(ii_jj_kk) == expected


def test_ufunc_support():
    f = Function('f')
    g = Function('g')
    x = IndexedBase('x')
    y = IndexedBase('y')
    i, j = Idx('i'), Idx('j')

    assert get_indices(f(x[i])) == ({i}, {})
    assert get_indices(f(x[i], y[j])) == ({i, j}, {})
    assert get_indices(f(y[i])*g(x[i])) == (set(), {})
    assert get_indices(f(a, x[i])) == ({i}, {})
    assert get_indices(f(a, y[i], x[j])*g(x[i])) == ({j}, {})
    assert get_indices(g(f(x[i]))) == ({i}, {})

    assert get_contraction_structure(f(x[i])) == {None: {f(x[i])}}
    assert get_contraction_structure(
        f(y[i])*g(x[i])) == {(i,): {f(y[i])*g(x[i])}}
    assert get_contraction_structure(
        f(y[i])*g(f(x[i]))) == {(i,): {f(y[i])*g(f(x[i]))}}
    assert get_contraction_structure(
        f(x[j], y[i])*g(x[i])) == {(i,): {f(x[j], y[i])*g(x[i])}}
