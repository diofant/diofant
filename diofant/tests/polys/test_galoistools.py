from diofant import FF, ring
from diofant.polys.galoistools import dup_gf_random


__all__ = ()


def test_dup_gf_random():
    F11 = FF(11)

    for n in range(8):
        f = dup_gf_random(n, F11, irreducible=True)
        R, x = ring('x', F11)
        f = R.from_list(f)

        assert f.is_irreducible is True
