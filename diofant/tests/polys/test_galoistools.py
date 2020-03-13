import pytest

from diofant import FF, ring
from diofant.polys.galoistools import (dup_gf_compose_mod, dup_gf_irreducible,
                                       dup_gf_irreducible_p, dup_gf_pow_mod,
                                       dup_gf_primitive_p, dup_gf_trace_map)
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

    F9 = FF(3, [1, 2, 2])
    R, x = ring('x', F9)

    f = R.to_dense(x**3 + F9(8)*x**2 + F9(8)*x + F9(4))

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert dup_gf_irreducible_p(f, R.domain) is False

    F27 = FF(3, [1, 2, 0, 1])
    R, x = ring('x', F27)

    f = R.to_dense(x**3 + F27(8)*x**2 + F27(19)*x + F27(24))

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert dup_gf_irreducible_p(f, R.domain) is True


def test_dup_gf_primitive_p():
    R, x = ring('x', FF(11))

    f = R.to_dense(x**2 + 2*x)

    assert dup_gf_primitive_p(f, R.domain) is False

    f = R.to_dense(x**2 + 7*x + 5)

    assert dup_gf_primitive_p(f, R.domain) is False

    f = R.to_dense(x**2 + 2*x + 6)

    assert dup_gf_irreducible_p(f, R.domain) is True
    assert dup_gf_primitive_p(f, R.domain) is True
