import pytest

from diofant import nextprime, pi
from diofant.domains import ZZ
from diofant.polys import polyconfig as config
from diofant.polys.galoistools import (gf_add, gf_add_ground, gf_berlekamp,
                                       gf_compose_mod, gf_ddf_shoup,
                                       gf_ddf_zassenhaus, gf_div, gf_edf_shoup,
                                       gf_edf_zassenhaus, gf_factor_sqf,
                                       gf_frobenius_map,
                                       gf_frobenius_monomial_base,
                                       gf_from_dict, gf_gcd, gf_irred_p_ben_or,
                                       gf_irred_p_rabin, gf_irreducible,
                                       gf_irreducible_p, gf_monic, gf_mul,
                                       gf_mul_ground, gf_pow_mod, gf_Qbasis,
                                       gf_Qmatrix, gf_quo, gf_rem, gf_sqr,
                                       gf_sub, gf_sub_ground, gf_trace_map)


__all__ = ()


def test_gf_from_dict():
    f = {11: 12, 6: 2, 0: 25}
    g = [1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 3]

    assert gf_from_dict(f, 11, ZZ) == g

    f = {11: -5, 4: 0, 3: 1, 0: 12}
    g = [6, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1]

    assert gf_from_dict(f, 11, ZZ) == g


def test_gf_monic():
    assert gf_monic([], 11, ZZ) == (0, [])

    assert gf_monic([1], 11, ZZ) == (1, [1])
    assert gf_monic([2], 11, ZZ) == (2, [1])

    assert gf_monic([1, 2, 3, 4], 11, ZZ) == (1, [1, 2, 3, 4])
    assert gf_monic([2, 3, 4, 5], 11, ZZ) == (2, [1, 7, 2, 8])


def test_gf_arith():
    assert gf_add_ground([], 0, 11, ZZ) == []
    assert gf_sub_ground([], 0, 11, ZZ) == []

    assert gf_add_ground([], 3, 11, ZZ) == [3]
    assert gf_sub_ground([], 3, 11, ZZ) == [8]

    assert gf_add_ground([1], 3, 11, ZZ) == [4]
    assert gf_sub_ground([1], 3, 11, ZZ) == [9]

    assert gf_add_ground([8], 3, 11, ZZ) == []
    assert gf_sub_ground([3], 3, 11, ZZ) == []

    assert gf_add_ground([1, 2, 3], 3, 11, ZZ) == [1, 2, 6]
    assert gf_sub_ground([1, 2, 3], 3, 11, ZZ) == [1, 2, 0]

    assert gf_mul_ground([], 0, 11, ZZ) == []
    assert gf_mul_ground([], 1, 11, ZZ) == []

    assert gf_mul_ground([1], 0, 11, ZZ) == []
    assert gf_mul_ground([1], 1, 11, ZZ) == [1]

    assert gf_mul_ground([1, 2, 3], 0, 11, ZZ) == []
    assert gf_mul_ground([1, 2, 3], 1, 11, ZZ) == [1, 2, 3]
    assert gf_mul_ground([1, 2, 3], 7, 11, ZZ) == [7, 3, 10]

    assert gf_add([], [], 11, ZZ) == []
    assert gf_add([1], [], 11, ZZ) == [1]
    assert gf_add([], [1], 11, ZZ) == [1]
    assert gf_add([1], [1], 11, ZZ) == [2]
    assert gf_add([1], [2], 11, ZZ) == [3]

    assert gf_add([1, 2], [1], 11, ZZ) == [1, 3]
    assert gf_add([1], [1, 2], 11, ZZ) == [1, 3]

    assert gf_add([1, 2, 3], [8, 9, 10], 11, ZZ) == [9, 0, 2]

    assert gf_sub([], [], 11, ZZ) == []
    assert gf_sub([1], [], 11, ZZ) == [1]
    assert gf_sub([], [1], 11, ZZ) == [10]
    assert gf_sub([1], [1], 11, ZZ) == []
    assert gf_sub([1], [2], 11, ZZ) == [10]

    assert gf_sub([1, 2], [1], 11, ZZ) == [1, 1]
    assert gf_sub([1], [1, 2], 11, ZZ) == [10, 10]

    assert gf_sub([3, 2, 1], [8, 9, 10], 11, ZZ) == [6, 4, 2]

    assert gf_mul([], [], 11, ZZ) == []
    assert gf_mul([], [1], 11, ZZ) == []
    assert gf_mul([1], [], 11, ZZ) == []
    assert gf_mul([1], [1], 11, ZZ) == [1]
    assert gf_mul([5], [7], 11, ZZ) == [2]

    assert gf_mul([3, 0, 0, 6, 1, 2], [4, 0, 1, 0], 11, ZZ) == [1, 0,
                                                                3, 2, 4, 3, 1, 2, 0]
    assert gf_mul([4, 0, 1, 0], [3, 0, 0, 6, 1, 2], 11, ZZ) == [1, 0,
                                                                3, 2, 4, 3, 1, 2, 0]

    assert gf_mul([2, 0, 0, 1, 7], [2, 0, 0, 1, 7], 11, ZZ) == [4, 0,
                                                                0, 4, 6, 0, 1, 3, 5]

    assert gf_sqr([], 11, ZZ) == []
    assert gf_sqr([2], 11, ZZ) == [4]
    assert gf_sqr([1, 2], 11, ZZ) == [1, 4, 4]

    assert gf_sqr([2, 0, 0, 1, 7], 11, ZZ) == [4, 0, 0, 4, 6, 0, 1, 3, 5]


def test_gf_division():
    pytest.raises(ZeroDivisionError, lambda: gf_div([1, 2, 3], [], 11, ZZ))
    pytest.raises(ZeroDivisionError, lambda: gf_rem([1, 2, 3], [], 11, ZZ))
    pytest.raises(ZeroDivisionError, lambda: gf_quo([1, 2, 3], [], 11, ZZ))
    pytest.raises(ZeroDivisionError, lambda: gf_quo([1, 2, 3], [], 11, ZZ))

    assert gf_div([1], [1, 2, 3], 7, ZZ) == ([], [1])
    assert gf_rem([1], [1, 2, 3], 7, ZZ) == [1]
    assert gf_quo([1], [1, 2, 3], 7, ZZ) == []

    f = [5, 4, 3, 2, 1, 0]
    g = [1, 2, 3]
    q = [5, 1, 0, 6]
    r = [3, 3]

    assert gf_div(f, g, 7, ZZ) == (q, r)
    assert gf_rem(f, g, 7, ZZ) == r
    assert gf_quo(f, g, 7, ZZ) == q

    f = [5, 4, 3, 2, 1, 0]
    g = [1, 2, 3, 0]
    q = [5, 1, 0]
    r = [6, 1, 0]

    assert gf_div(f, g, 7, ZZ) == (q, r)
    assert gf_rem(f, g, 7, ZZ) == r
    assert gf_quo(f, g, 7, ZZ) == q

    assert gf_quo([1, 2, 1], [1, 1], 11, ZZ) == [1, 1]


def test_gf_powering():
    assert gf_pow_mod([1, 0, 0, 1, 8], 0, [2, 0, 7], 11, ZZ) == [1]
    assert gf_pow_mod([1, 0, 0, 1, 8], 1, [2, 0, 7], 11, ZZ) == [1, 1]
    assert gf_pow_mod([1, 0, 0, 1, 8], 2, [2, 0, 7], 11, ZZ) == [2, 3]
    assert gf_pow_mod([1, 0, 0, 1, 8], 5, [2, 0, 7], 11, ZZ) == [7, 8]
    assert gf_pow_mod([1, 0, 0, 1, 8], 8, [2, 0, 7], 11, ZZ) == [1, 5]
    assert gf_pow_mod([1, 0, 0, 1, 8], 45, [2, 0, 7], 11, ZZ) == [5, 4]


def test_gf_gcd():
    assert gf_gcd([], [], 11, ZZ) == []
    assert gf_gcd([2], [], 11, ZZ) == [1]
    assert gf_gcd([], [2], 11, ZZ) == [1]
    assert gf_gcd([2], [2], 11, ZZ) == [1]

    assert gf_gcd([], [1, 0], 11, ZZ) == [1, 0]
    assert gf_gcd([1, 0], [], 11, ZZ) == [1, 0]

    assert gf_gcd([3, 0], [3, 0], 11, ZZ) == [1, 0]
    assert gf_gcd([1, 8, 7], [1, 7, 1, 7], 11, ZZ) == [1, 7]


def test_gf_compose_mod():
    assert gf_compose_mod([], [1, 0], [1, 0], 11, ZZ) == []

    f = [1, 1, 4, 9, 1]
    g = [1, 1, 1]
    h = [1, 0, 0, 2]

    assert gf_compose_mod(g, h, f, 11, ZZ) == [3, 9, 6, 10]


def test_gf_trace_map():
    f = [1, 1, 4, 9, 1]
    a = [1, 1, 1]
    c = [1, 0]
    b = gf_pow_mod(c, 11, f, 11, ZZ)

    assert gf_trace_map(a, b, c, 0, f, 11, ZZ) == ([1, 1, 1], [1, 1, 1])
    assert gf_trace_map(a, b, c, 1, f, 11, ZZ) == ([5, 2, 10, 3], [5, 3, 0, 4])
    assert gf_trace_map(a, b, c, 2, f, 11, ZZ) == ([5, 9, 5, 3], [10, 1, 5, 7])
    assert gf_trace_map(a, b, c, 3, f, 11, ZZ) == ([1, 10, 6, 0], [7])
    assert gf_trace_map(a, b, c, 4, f, 11, ZZ) == ([1, 1, 1], [1, 1, 8])
    assert gf_trace_map(a, b, c, 5, f, 11, ZZ) == ([5, 2, 10, 3], [5, 3, 0, 0])
    assert gf_trace_map(a, b, c, 11, f, 11, ZZ) == ([1, 10, 6, 0], [10])


def test_gf_irreducible():
    assert gf_irreducible_p(gf_irreducible(1, 11, ZZ), 11, ZZ) is True
    assert gf_irreducible_p(gf_irreducible(2, 11, ZZ), 11, ZZ) is True
    assert gf_irreducible_p(gf_irreducible(3, 11, ZZ), 11, ZZ) is True
    assert gf_irreducible_p(gf_irreducible(4, 11, ZZ), 11, ZZ) is True
    assert gf_irreducible_p(gf_irreducible(5, 11, ZZ), 11, ZZ) is True
    assert gf_irreducible_p(gf_irreducible(6, 11, ZZ), 11, ZZ) is True
    assert gf_irreducible_p(gf_irreducible(7, 11, ZZ), 11, ZZ) is True


def test_gf_irreducible_p():
    assert gf_irred_p_ben_or([7], 11, ZZ) is True
    assert gf_irred_p_ben_or([7, 3], 11, ZZ) is True
    assert gf_irred_p_ben_or([7, 3, 1], 11, ZZ) is False

    assert gf_irred_p_rabin([7], 11, ZZ) is True
    assert gf_irred_p_rabin([7, 3], 11, ZZ) is True
    assert gf_irred_p_rabin([7, 3, 1], 11, ZZ) is False

    assert gf_irred_p_ben_or([2, 3, 4, 5, 6], 13, ZZ) is False
    assert gf_irred_p_ben_or([2, 3, 4, 5, 8], 13, ZZ) is True

    with config.using(gf_irred_method='ben-or'):
        assert gf_irreducible_p([7], 11, ZZ) is True
        assert gf_irreducible_p([7, 3], 11, ZZ) is True
        assert gf_irreducible_p([7, 3, 1], 11, ZZ) is False

    with config.using(gf_irred_method='rabin'):
        assert gf_irreducible_p([7], 11, ZZ) is True
        assert gf_irreducible_p([7, 3], 11, ZZ) is True
        assert gf_irreducible_p([7, 3, 1], 11, ZZ) is False

    with config.using(gf_irred_method='other'):
        pytest.raises(KeyError, lambda: gf_irreducible_p([7], 11, ZZ))

    f = [1, 9, 9, 13, 16, 15, 6, 7, 7, 7, 10]
    g = [1, 7, 16, 7, 15, 13, 13, 11, 16, 10, 9]

    h = gf_mul(f, g, 17, ZZ)

    assert gf_irred_p_ben_or(f, 17, ZZ) is True
    assert gf_irred_p_ben_or(g, 17, ZZ) is True

    assert gf_irred_p_ben_or(h, 17, ZZ) is False

    assert gf_irred_p_rabin(f, 17, ZZ) is True
    assert gf_irred_p_rabin(g, 17, ZZ) is True

    assert gf_irred_p_rabin(h, 17, ZZ) is False


def test_gf_frobenius_map():
    f = [2, 0, 1, 0, 2, 2, 0, 2, 2, 2]
    g = [1, 1, 0, 2, 0, 1, 0, 2, 0, 1]
    p = 3
    b = gf_frobenius_monomial_base(g, p, ZZ)
    h = gf_frobenius_map(f, g, b, p, ZZ)
    h1 = gf_pow_mod(f, p, g, p, ZZ)
    assert h == h1


def test_gf_berlekamp():
    f = [1, 8, 1, 8, 10, 8, 1]

    Q = [[1, 0, 0, 0, 0, 0],
         [3, 5, 8, 8, 6, 5],
         [3, 6, 6, 1, 10, 0],
         [9, 4, 10, 3, 7, 9],
         [7, 8, 10, 0, 0, 8],
         [8, 10, 7, 8, 10, 8]]

    V = [[1, 0, 0, 0, 0, 0],
         [0, 1, 1, 1, 1, 0],
         [0, 0, 7, 9, 0, 1]]

    assert gf_Qmatrix(f, 11, ZZ) == Q
    assert gf_Qbasis(Q, 11, ZZ) == V

    assert gf_berlekamp(f, 11, ZZ) == [[1, 1], [1, 5, 3], [1, 2, 3, 4]]

    f = [1, 0, 1, 0, 10, 10, 8, 2, 8]

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

    assert gf_Qmatrix(f, 13, ZZ) == Q
    assert gf_Qbasis(Q, 13, ZZ) == V

    assert gf_berlekamp(f, 13, ZZ) == [[1, 3], [1, 8, 4, 12], [1, 2, 3, 4, 6]]


def test_gf_ddf():
    f = gf_from_dict({15: ZZ(1), 0: ZZ(-1)}, 11, ZZ)
    g = [([1, 0, 0, 0, 0, 10], 1),
         ([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1], 2)]

    assert gf_ddf_zassenhaus(f, 11, ZZ) == g
    assert gf_ddf_shoup(f, 11, ZZ) == g

    f = gf_from_dict({63: ZZ(1), 0: ZZ(1)}, 2, ZZ)
    g = [([1, 1], 1),
         ([1, 1, 1], 2),
         ([1, 1, 1, 1, 1, 1, 1], 3),
         ([1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0,
           0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1], 6)]

    assert gf_ddf_zassenhaus(f, 2, ZZ) == g
    assert gf_ddf_shoup(f, 2, ZZ) == g

    f = gf_from_dict({6: ZZ(1), 5: ZZ(-1), 4: ZZ(1), 3: ZZ(1), 1: ZZ(-1)}, 3, ZZ)
    g = [([1, 1, 0], 1),
         ([1, 1, 0, 1, 2], 2)]

    assert gf_ddf_zassenhaus(f, 3, ZZ) == g
    assert gf_ddf_shoup(f, 3, ZZ) == g

    f = [1, 2, 5, 26, 677, 436, 791, 325, 456, 24, 577]
    g = [([1, 701], 1),
         ([1, 110, 559, 532, 694, 151, 110, 70, 735, 122], 9)]

    assert gf_ddf_zassenhaus(f, 809, ZZ) == g
    assert gf_ddf_shoup(f, 809, ZZ) == g

    p = ZZ(nextprime(int((2**15*pi))))
    f = gf_from_dict({15: 1, 1: 1, 0: 1}, p, ZZ)
    g = [([1, 22730, 68144], 2),
         ([1, 64876, 83977, 10787, 12561, 68608, 52650, 88001, 84356], 4),
         ([1, 15347, 95022, 84569, 94508, 92335], 5)]

    assert gf_ddf_zassenhaus(f, p, ZZ) == g
    assert gf_ddf_shoup(f, p, ZZ) == g


def test_gf_edf():
    f = [1, 1, 0, 1, 2]
    g = [[1, 0, 1], [1, 1, 2]]

    assert gf_edf_zassenhaus(f, 2, 3, ZZ) == g
    assert gf_edf_shoup(f, 2, 3, ZZ) == g


def test_gf_factor_sqf():
    assert gf_factor_sqf([], 11, ZZ) == (0, [])
    assert gf_factor_sqf([1], 11, ZZ) == (1, [])
    assert gf_factor_sqf([1, 1], 11, ZZ) == (1, [[1, 1]])
    assert gf_factor_sqf([2, 3], 11, ZZ) == (2, [[1, 7]])

    with config.using(gf_factor_method='berlekamp'):
        assert gf_factor_sqf([], 11, ZZ) == (0, [])
        assert gf_factor_sqf([1], 11, ZZ) == (1, [])
        assert gf_factor_sqf([1, 1], 11, ZZ) == (1, [[1, 1]])
        assert gf_factor_sqf([1, 0], 11, ZZ) == (1, [[1, 0]])

    with config.using(gf_factor_method='zassenhaus'):
        assert gf_factor_sqf([], 11, ZZ) == (0, [])
        assert gf_factor_sqf([1], 11, ZZ) == (1, [])
        assert gf_factor_sqf([1, 1], 11, ZZ) == (1, [[1, 1]])
        assert gf_factor_sqf([1, 0], 11, ZZ) == (1, [[1, 0]])

    with config.using(gf_factor_method='shoup'):
        assert gf_factor_sqf([], 11, ZZ) == (0, [])
        assert gf_factor_sqf([1], 11, ZZ) == (1, [])
        assert gf_factor_sqf([1, 1], 11, ZZ) == (1, [[1, 1]])
        assert gf_factor_sqf([1, 0], 11, ZZ) == (1, [[1, 0]])

    f, p = [1, 0, 0, 1, 0], 2
    g = (1, [[1, 0],
             [1, 1],
             [1, 1, 1]])

    with config.using(gf_factor_method='berlekamp'):
        assert gf_factor_sqf(f, p, ZZ) == g

    with config.using(gf_factor_method='zassenhaus'):
        assert gf_factor_sqf(f, p, ZZ) == g

    with config.using(gf_factor_method='shoup'):
        assert gf_factor_sqf(f, p, ZZ) == g

    # Gathen polynomials: x**n + x + 1 (mod p > 2**n * pi)

    p = ZZ(nextprime(int((2**15*pi))))
    f = gf_from_dict({15: 1, 1: 1, 0: 1}, p, ZZ)
    g = (1, [[1, 22730, 68144],
             [1, 81553, 77449, 86810, 4724],
             [1, 86276, 56779, 14859, 31575],
             [1, 15347, 95022, 84569, 94508, 92335]])

    with config.using(gf_factor_method='zassenhaus'):
        assert gf_factor_sqf(f, p, ZZ) == g

    with config.using(gf_factor_method='shoup'):
        assert gf_factor_sqf(f, p, ZZ) == g

    # Shoup polynomials: f = a_0 x**n + a_1 x**(n-1) + ... + a_n
    # (mod p > 2**(n-2) * pi), where a_n = a_{n-1}**2 + 1, a_0 = 1

    p = ZZ(nextprime(int((2**4*pi))))
    f = [1, 2, 5, 26, 41, 39, 38]

    g = (1, [[1, 44, 26],
             [1, 11, 25, 18, 30]])

    with config.using(gf_factor_method='zassenhaus'):
        assert gf_factor_sqf(f, p, ZZ) == g

    with config.using(gf_factor_method='shoup'):
        assert gf_factor_sqf(f, p, ZZ) == g
