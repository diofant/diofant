from diofant import FF, ring
from diofant.polys.galoistools import dup_gf_irreducible


__all__ = ()


def test_dup_gf_irreducible():
    F11 = FF(11)

    for n in range(8):
        f = dup_gf_irreducible(n, F11)
        R, x = ring('x', F11)
        f = R.from_list(f)

        assert f.is_irreducible is True
