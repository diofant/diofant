from diofant import ZZ, ring
from diofant.polys.factorization_alg_field import (_distinct_prime_divisors,
                                                   _sqf_p)


def test__distinct_prime_divisors():
    s = [20, 6, 7]
    assert _distinct_prime_divisors(s, ZZ) == [5, 3, 7]

    s = [15, 18, 11]
    assert _distinct_prime_divisors(s, ZZ) == [5, 2, 11]


def test__sqf_p():
    *_, z = ring('x z', ZZ)

    assert _sqf_p(z**2, (z**2 - 2).drop(0), 2) is True
