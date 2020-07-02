from diofant import FF, ring
from diofant.polys.galoistools import dup_gf_irreducible, dup_gf_primitive_p


__all__ = ()


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
