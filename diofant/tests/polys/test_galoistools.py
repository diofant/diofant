import pytest

from diofant import FF, nextprime, pi, ring
from diofant.polys.galoistools import (dup_gf_berlekamp, dup_gf_compose_mod,
                                       dup_gf_ddf_shoup, dup_gf_ddf_zassenhaus,
                                       dup_gf_edf_shoup, dup_gf_edf_zassenhaus,
                                       dup_gf_frobenius_map,
                                       dup_gf_frobenius_monomial_base,
                                       dup_gf_irreducible,
                                       dup_gf_irreducible_p, dup_gf_pow_mod,
                                       dup_gf_Qbasis, dup_gf_Qmatrix,
                                       dup_gf_trace_map)
from diofant.polys.polyconfig import using


__all__ = ()


def test_gf_powering():
    R, x = ring('x', FF(11))

    f = R.to_dense(x**4 + x + 8)
    g = R.to_dense(2*x**2 + 7)

    assert dup_gf_pow_mod(f, 0, g, R.domain) == [1]
    assert dup_gf_pow_mod(f, 1, g, R.domain) == [1, 1]
    assert dup_gf_pow_mod(f, 2, g, R.domain) == [2, 3]
    assert dup_gf_pow_mod(f, 5, g, R.domain) == [7, 8]
    assert dup_gf_pow_mod(f, 8, g, R.domain) == [1, 5]
    assert dup_gf_pow_mod(f, 45, g, R.domain) == [5, 4]


def test_dup_gf_compose_mod():
    R, x = ring('x', FF(11))

    g = []
    h = R.to_dense(x)
    f = h.copy()

    assert dup_gf_compose_mod(g, h, h, R.domain) == []

    f = R.to_dense(x**4 + x**3 + 4*x**2 + 9*x + 1)
    g = R.to_dense(x**2 + x + 1)
    h = R.to_dense(x**3 + 2)

    assert dup_gf_compose_mod(g, h, f, R.domain) == [3, 9, 6, 10]


def test_dup_gf_trace_map():
    R, x = ring('x', FF(11))

    f = R.to_dense(x**4 + x**3 + 4*x**2 + 9*x + 1)
    a = R.to_dense(x**2 + x + 1)
    c = R.to_dense(x)
    b = dup_gf_pow_mod(c, 11, f, R.domain)

    assert dup_gf_trace_map(a, b, c, 0, f, R.domain) == ([1, 1, 1], [1, 1, 1])
    assert dup_gf_trace_map(a, b, c, 1, f, R.domain) == ([5, 2, 10, 3], [5, 3, 0, 4])
    assert dup_gf_trace_map(a, b, c, 2, f, R.domain) == ([5, 9, 5, 3], [10, 1, 5, 7])
    assert dup_gf_trace_map(a, b, c, 3, f, R.domain) == ([1, 10, 6, 0], [7])
    assert dup_gf_trace_map(a, b, c, 4, f, R.domain) == ([1, 1, 1], [1, 1, 8])
    assert dup_gf_trace_map(a, b, c, 5, f, R.domain) == ([5, 2, 10, 3], [5, 3, 0, 0])
    assert dup_gf_trace_map(a, b, c, 11, f, R.domain) == ([1, 10, 6, 0], [10])


def test_dup_gf_irreducible():
    F11 = FF(11)

    for n in range(8):
        assert dup_gf_irreducible_p(dup_gf_irreducible(n, F11), F11) is True


def test_dup_gf_irreducible_p():
    R, x = ring('x', FF(11))

    f = R.to_dense(R(7))
    g = R.to_dense(7*x + 3)
    h = R.to_dense(7*x**2 + 3*x + 1)

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert dup_gf_irreducible_p(f, R.domain) is True
            assert dup_gf_irreducible_p(g, R.domain) is True
            assert dup_gf_irreducible_p(h, R.domain) is False

    with using(gf_irred_method='other'):
        pytest.raises(KeyError, lambda: dup_gf_irreducible_p(f, R.domain))

    R, x = ring('x', FF(13))

    f = R.to_dense(2*x**4 + 3*x**3 + 4*x**2 + 5*x + 6)
    g = R.to_dense(2*x**4 + 3*x**3 + 4*x**2 + 5*x + 8)

    with using(gf_irred_method='ben-or'):
        assert dup_gf_irreducible_p(f, R.domain) is False
        assert dup_gf_irreducible_p(g, R.domain) is True

    R, x = ring('x', FF(17))

    f = (x**10 + 9*x**9 + 9*x**8 + 13*x**7 + 16*x**6 + 15*x**5 +
         6*x**4 + 7*x**3 + 7*x**2 + 7*x + 10)
    g = (x**10 + 7*x**9 + 16*x**8 + 7*x**7 + 15*x**6 + 13*x**5 + 13*x**4 +
         11*x**3 + 16*x**2 + 10*x + 9)
    h = f*g
    f, g, h = map(R.to_dense, (f, g, h))

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert dup_gf_irreducible_p(f, R.domain) is True
            assert dup_gf_irreducible_p(g, R.domain) is True
            assert dup_gf_irreducible_p(h, R.domain) is False


def test_dup_gf_frobenius_map():
    R, x = ring('x', FF(3))

    f = R.to_dense(2*x**9 + x**7 + 2*x**5 + 2*x**4 + 2*x**2 + 2*x + 2)
    g = R.to_dense(x**9 + x**8 + 2*x**6 + x**4 + 2*x**2 + 1)
    b = dup_gf_frobenius_monomial_base(g, R.domain)
    h = dup_gf_frobenius_map(f, g, b, R.domain)
    h1 = dup_gf_pow_mod(f, R.domain.mod, g, R.domain)

    assert h == h1


def test_dup_gf_berlekamp():
    R, x = ring('x', FF(11))

    f = R.to_dense(x**6 + 8*x**5 + x**4 + 8*x**3 + 10*x**2 + 8*x + 1)

    Q = [[1, 0, 0, 0, 0, 0],
         [3, 5, 8, 8, 6, 5],
         [3, 6, 6, 1, 10, 0],
         [9, 4, 10, 3, 7, 9],
         [7, 8, 10, 0, 0, 8],
         [8, 10, 7, 8, 10, 8]]

    V = [[1, 0, 0, 0, 0, 0],
         [0, 1, 1, 1, 1, 0],
         [0, 0, 7, 9, 0, 1]]

    assert dup_gf_Qmatrix(f, R.domain) == Q
    assert dup_gf_Qbasis(Q, R.domain) == V

    assert dup_gf_berlekamp(f, R.domain) == [[1, 1], [1, 5, 3], [1, 2, 3, 4]]

    R, x = ring('x', FF(13))

    f = R.to_dense(x**8 + x**6 + 10*x**4 + 10*x**3 + 8*x**2 + 2*x + 8)
    Q = [[1, 0, 0, 0, 0, 0, 0, 0],
         [2, 1, 7, 11, 10, 12, 5, 11],
         [3, 6, 4, 3, 0, 4, 7, 2],
         [4, 3, 6, 5, 1, 6, 2, 3],
         [2, 11, 8, 8, 3, 1, 3, 11],
         [6, 11, 8, 6, 2, 7, 10, 9],
         [5, 11, 7, 10, 0, 11, 7, 12],
         [3, 3, 12, 5, 0, 11, 9, 12]]
    V = [[1, 0, 0, 0, 0, 0, 0, 0],
         [0, 5, 5, 0, 9, 5, 1, 0],
         [0, 9, 11, 9, 10, 12, 0, 1]]

    assert dup_gf_Qmatrix(f, R.domain) == Q
    assert dup_gf_Qbasis(Q, R.domain) == V

    assert dup_gf_berlekamp(f, R.domain) == [[1, 3], [1, 8, 4, 12], [1, 2, 3, 4, 6]]


def test_dup_gf_ddf():
    R, x = ring('x', FF(11))

    f = R.to_dense(x**15 - 1)
    g = [([1, 0, 0, 0, 0, 10], 1),
         ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], 2)]

    assert dup_gf_ddf_zassenhaus(f, R.domain) == g
    assert dup_gf_ddf_shoup(f, R.domain) == g

    R, x = ring('x', FF(2))

    f = R.to_dense(x**63 + 1)
    g = [([1, 1], 1),
         ([1, 1, 1], 2),
         ([1, 1, 1, 1, 1, 1, 1], 3),
         ([1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1], 6)]

    assert dup_gf_ddf_zassenhaus(f, R.domain) == g
    assert dup_gf_ddf_shoup(f, R.domain) == g

    R, x = ring('x', FF(3))

    f = R.to_dense(x**6 - x**5 + x**4 + x**3 - x)
    g = [([1, 1, 0], 1),
         ([1, 1, 0, 1, 2], 2)]

    assert dup_gf_ddf_zassenhaus(f, R.domain) == g
    assert dup_gf_ddf_shoup(f, R.domain) == g

    R, x = ring('x', FF(809))

    f = R.to_dense(x**10 + 2*x**9 + 5*x**8 + 26*x**7 + 677*x**6 + 436*x**5 +
                   791*x**4 + 325*x**3 + 456*x**2 + 24*x + 577)
    g = [([1, 701], 1),
         ([1, 110, 559, 532, 694, 151, 110, 70, 735, 122], 9)]

    assert dup_gf_ddf_zassenhaus(f, R.domain) == g
    assert dup_gf_ddf_shoup(f, R.domain) == g

    R, x = ring('x', FF(nextprime(int((2**15*pi)))))

    f = R.to_dense(x**15 + x + 1)
    g = [([1, 22730, 68144], 2),
         ([1, 64876, 83977, 10787, 12561, 68608, 52650, 88001, 84356], 4),
         ([1, 15347, 95022, 84569, 94508, 92335], 5)]

    assert dup_gf_ddf_zassenhaus(f, R.domain) == g
    assert dup_gf_ddf_shoup(f, R.domain) == g


def test_dup_gf_edf():
    R, x = ring('x', FF(3))

    f = R.to_dense(x**4 + x**3 + x + 2)
    g = [[1, 0, 1], [1, 1, 2]]

    assert dup_gf_edf_zassenhaus(f, 2, R.domain) == g
    assert dup_gf_edf_shoup(f, 2, R.domain) == g
