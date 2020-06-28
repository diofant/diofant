from diofant import FF, ring
from diofant.polys.galoistools import (dup_gf_compose_mod, dup_gf_irreducible,
                                       dup_gf_primitive_p, dup_gf_trace_map)


__all__ = ()


def test_dup_gf_compose_mod():
    R, x = ring('x', FF(11))

    g = []
    h = x.to_dense()
    f = h.copy()

    assert dup_gf_compose_mod(g, h, h, R.domain) == []

    f = (x**4 + x**3 + 4*x**2 + 9*x + 1).to_dense()
    g = (x**2 + x + 1).to_dense()
    h = (x**3 + 2).to_dense()

    assert dup_gf_compose_mod(g, h, f, R.domain) == [3, 9, 6, 10]


def test_dup_gf_trace_map():
    R, x = ring('x', FF(11))

    f = x**4 + x**3 + 4*x**2 + 9*x + 1
    a = x**2 + x + 1
    c = x
    b = pow(c, 11, f)

    f, a, b, c = map(lambda _: _.to_dense(), (f, a, b, c))

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
        f = dup_gf_irreducible(n, F11)
        R, x = ring('x', F11)
        f = R.from_list(f)

        assert f.is_irreducible is True


def test_dup_gf_primitive_p():
    R, x = ring('x', FF(11))

    f = (x**2 + 2*x).to_dense()

    assert dup_gf_primitive_p(f, R.domain) is False

    f = (x**2 + 7*x + 5).to_dense()

    assert dup_gf_primitive_p(f, R.domain) is False

    f = x**2 + 2*x + 6

    assert f.is_irreducible is True

    f = f.to_dense()

    assert dup_gf_primitive_p(f, R.domain) is True
