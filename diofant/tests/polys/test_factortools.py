"""Tools for polynomial factorization routines in characteristic zero. """

import pytest

from diofant import I, nextprime, pi, sin, sqrt
from diofant.domains import EX, FF, QQ, RR, ZZ
from diofant.polys import polyconfig as config
from diofant.polys.factortools import dmp_zz_diophantine
from diofant.polys.polyerrors import DomainError
from diofant.polys.rings import ring
from diofant.polys.specialpolys import f_polys, w_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = f_polys()
w_1, w_2 = w_polys()


def test_dmp_trial_division():
    R, x = ring("x", ZZ)
    assert R.dmp_trial_division(x**5 + 8*x**4 + 25*x**3 + 38*x**2 + 28*x + 8, (x + 1, x + 2)) == [(x + 1, 2), (x + 2, 3)]

    R, x, y = ring("x,y", ZZ)
    assert R.dmp_trial_division(x**5 + 8*x**4 + 25*x**3 + 38*x**2 + 28*x + 8, (x + 1, x + 2)) == [(x + 1, 2), (x + 2, 3)]


def test_dmp_zz_mignotte_bound():
    R, x = ring("x", ZZ)
    assert R.dmp_zz_mignotte_bound(2*x**2 + 3*x + 4) == 32

    R, x, y = ring("x,y", ZZ)
    assert R.dmp_zz_mignotte_bound(2*x**2 + 3*x + 4) == 32


def test_dup_zz_hensel_step():
    R, x = ring("x", ZZ)

    f = x**4 - 1
    g = x**3 + 2*x**2 - x - 2
    h = x - 2
    s = -2
    t = 2*x**2 - 2*x - 1

    G, H, S, T = R.dup_zz_hensel_step(5, f, g, h, s, t)

    assert G == x**3 + 7*x**2 - x - 7
    assert H == x - 7
    assert S == 8
    assert T == -8*x**2 - 12*x - 1


def test_dup_zz_hensel_lift():
    R, x = ring("x", ZZ)

    f = x**4 - 1
    F = [x - 1, x - 2, x + 2, x + 1]

    assert R.dup_zz_hensel_lift(ZZ(5), f, F, 4) == \
        [x - 1, x - 182, x + 182, x + 1]


def test_dup_zz_irreducible_p():
    R, x = ring("x", ZZ)

    assert R.dup_zz_irreducible_p(x) is None

    assert R.dup_zz_irreducible_p(3*x**4 + 2*x**3 + 6*x**2 + 8*x + 7) is None
    assert R.dup_zz_irreducible_p(3*x**4 + 2*x**3 + 6*x**2 + 8*x + 4) is None

    assert R.dup_zz_irreducible_p(3*x**4 + 2*x**3 + 6*x**2 + 8*x + 10) is True
    assert R.dup_zz_irreducible_p(3*x**4 + 2*x**3 + 6*x**2 + 8*x + 14) is True


def test_dup_cyclotomic_p():
    R, x = ring("x", ZZ)

    assert R.dup_cyclotomic_p(x - 1) is True
    assert R.dup_cyclotomic_p(x + 1) is True
    assert R.dup_cyclotomic_p(x**2 + x + 1) is True
    assert R.dup_cyclotomic_p(x**2 + 1) is True
    assert R.dup_cyclotomic_p(x**2 + 1, irreducible=True) is True
    assert R.dup_cyclotomic_p(x**4 + x**3 + x**2 + x + 1) is True
    assert R.dup_cyclotomic_p(x**2 - x + 1) is True
    assert R.dup_cyclotomic_p(x**6 + x**5 + x**4 + x**3 + x**2 + x + 1) is True
    assert R.dup_cyclotomic_p(x**4 + 1) is True
    assert R.dup_cyclotomic_p(x**6 + x**3 + 1) is True

    assert R.dup_cyclotomic_p(0) is False
    assert R.dup_cyclotomic_p(1) is False
    assert R.dup_cyclotomic_p(x) is False
    assert R.dup_cyclotomic_p(x + 2) is False
    assert R.dup_cyclotomic_p(3*x + 1) is False
    assert R.dup_cyclotomic_p(x**2 - 1) is False

    f = x**16 + x**14 - x**10 + x**8 - x**6 + x**2 + 1
    assert R.dup_cyclotomic_p(f) is False

    g = x**16 + x**14 - x**10 - x**8 - x**6 + x**2 + 1
    assert R.dup_cyclotomic_p(g) is True

    R, x = ring("x", QQ)
    assert R.dup_cyclotomic_p(x**2 + x + 1) is True
    assert R.dup_cyclotomic_p(x**2/2 + x + 1) is False

    R, x = ring("x", ZZ.poly_ring("y"))
    assert R.dup_cyclotomic_p(x**2 + x + 1) is False


def test_dup_zz_cyclotomic_poly():
    R, x = ring("x", ZZ)

    assert R.dup_zz_cyclotomic_poly(1) == x - 1
    assert R.dup_zz_cyclotomic_poly(2) == x + 1
    assert R.dup_zz_cyclotomic_poly(3) == x**2 + x + 1
    assert R.dup_zz_cyclotomic_poly(4) == x**2 + 1
    assert R.dup_zz_cyclotomic_poly(5) == x**4 + x**3 + x**2 + x + 1
    assert R.dup_zz_cyclotomic_poly(6) == x**2 - x + 1
    assert R.dup_zz_cyclotomic_poly(7) == x**6 + x**5 + x**4 + x**3 + x**2 + x + 1
    assert R.dup_zz_cyclotomic_poly(8) == x**4 + 1
    assert R.dup_zz_cyclotomic_poly(9) == x**6 + x**3 + 1


def test_dup_zz_cyclotomic_factor():
    R, x = ring("x", ZZ)

    assert R.dup_zz_cyclotomic_factor(0) is None
    assert R.dup_zz_cyclotomic_factor(1) is None

    assert R.dup_zz_cyclotomic_factor(2*x**10 - 1) is None
    assert R.dup_zz_cyclotomic_factor(x**10 - 3) is None
    assert R.dup_zz_cyclotomic_factor(x**10 + x**5 - 1) is None

    assert R.dup_zz_cyclotomic_factor(x + 1) == [x + 1]
    assert R.dup_zz_cyclotomic_factor(x - 1) == [x - 1]

    assert R.dup_zz_cyclotomic_factor(x**2 + 1) == [x**2 + 1]
    assert R.dup_zz_cyclotomic_factor(x**2 - 1) == [x - 1, x + 1]

    assert R.dup_zz_cyclotomic_factor(x**27 + 1) == \
        [x + 1, x**2 - x + 1, x**6 - x**3 + 1, x**18 - x**9 + 1]
    assert R.dup_zz_cyclotomic_factor(x**27 - 1) == \
        [x - 1, x**2 + x + 1, x**6 + x**3 + 1, x**18 + x**9 + 1]


def test_dup_zz_factor():
    R, x = ring("x", ZZ)

    assert R.dup_zz_factor(0) == (0, [])
    assert R.dup_zz_factor(7) == (7, [])
    assert R.dup_zz_factor(-7) == (-7, [])

    assert R.dup_zz_factor_sqf(0) == (0, [])
    assert R.dup_zz_factor_sqf(7) == (7, [])
    assert R.dup_zz_factor_sqf(-7) == (-7, [])

    assert R.dup_zz_factor(2*x + 4) == (2, [(x + 2, 1)])
    assert R.dup_zz_factor_sqf(2*x + 4) == (2, [x + 2])

    f = x**4 + x + 1

    for i in range(20):
        assert R.dup_zz_factor(f) == (1, [(f, 1)])

    assert R.dup_zz_factor(x**2 + 2*x + 2) == (1, [(x**2 + 2*x + 2, 1)])

    with config.using(use_irreducible_in_factor=True):
        assert R.dup_zz_factor(x**2 + 2*x + 2) == (1, [(x**2 + 2*x + 2, 1)])

    assert R.dup_zz_factor(18*x**2 + 12*x + 2) == (2, [(3*x + 1, 2)])

    with config.using(use_irreducible_in_factor=True):
        assert R.dup_zz_factor(18*x**2 + 12*x + 2) == (2, [(3*x + 1, 2)])

    assert R.dup_zz_factor(-9*x**2 + 1) == \
        (-1, [(3*x - 1, 1),
              (3*x + 1, 1)])

    with config.using(use_irreducible_in_factor=True):
        assert R.dup_zz_factor_sqf(3*x**4 + 2*x**3 +
                                   6*x**2 + 8*x + 10) == (1, [3*x**4 + 2*x**3 +
                                                              6*x**2 + 8*x + 10])

    assert R.dup_zz_factor_sqf(-9*x**2 + 1) == (-1, [3*x - 1, 3*x + 1])

    with config.using(use_irreducible_in_factor=True):
        assert R.dup_zz_factor_sqf(-9*x**2 + 1) == (-1, [3*x - 1, 3*x + 1])

    with config.using(use_cyclotomic_factor=False):
        assert R.dup_zz_factor_sqf(-9*x**2 + 1) == (-1, [3*x - 1, 3*x + 1])

    assert R.dup_zz_factor(x**3 - 6*x**2 + 11*x - 6) == \
        (1, [(x - 3, 1),
             (x - 2, 1),
             (x - 1, 1)])

    assert R.dup_zz_factor_sqf(x**3 - 6*x**2 + 11*x - 6) == \
        (1, [x - 3,
             x - 2,
             x - 1])

    assert R.dup_zz_factor(3*x**3 + 10*x**2 + 13*x + 10) == \
        (1, [(x + 2, 1),
             (3*x**2 + 4*x + 5, 1)])

    assert R.dup_zz_factor_sqf(3*x**3 + 10*x**2 + 13*x + 10) == \
        (1, [x + 2,
             3*x**2 + 4*x + 5])

    assert R.dup_zz_factor(-x**6 + x**2) == \
        (-1, [(x - 1, 1),
              (x + 1, 1),
              (x, 2),
              (x**2 + 1, 1)])

    f = 1080*x**8 + 5184*x**7 + 2099*x**6 + 744*x**5 + 2736*x**4 - 648*x**3 + 129*x**2 - 324

    assert R.dup_zz_factor(f) == \
        (1, [(5*x**4 + 24*x**3 + 9*x**2 + 12, 1),
             (216*x**4 + 31*x**2 - 27, 1)])

    f = -29802322387695312500000000000000000000*x**25 \
        + 2980232238769531250000000000000000*x**20 \
        + 1743435859680175781250000000000*x**15 \
        + 114142894744873046875000000*x**10 \
        - 210106372833251953125*x**5 \
        + 95367431640625

    assert R.dup_zz_factor(f) == \
        (-95367431640625, [(5*x - 1, 1),
                           (100*x**2 + 10*x - 1, 2),
                           (625*x**4 + 125*x**3 + 25*x**2 + 5*x + 1, 1),
                           (10000*x**4 - 3000*x**3 + 400*x**2 - 20*x + 1, 2),
                           (10000*x**4 + 2000*x**3 + 400*x**2 + 30*x + 1, 2)])

    f = x**10 - 1

    config.setup('USE_CYCLOTOMIC_FACTOR', True)
    F_0 = R.dup_zz_factor(f)

    config.setup('USE_CYCLOTOMIC_FACTOR', False)
    F_1 = R.dup_zz_factor(f)

    assert F_0 == F_1 == \
        (1, [(x - 1, 1),
             (x + 1, 1),
             (x**4 - x**3 + x**2 - x + 1, 1),
             (x**4 + x**3 + x**2 + x + 1, 1)])

    config.setup('USE_CYCLOTOMIC_FACTOR')

    f = x**10 + 1

    config.setup('USE_CYCLOTOMIC_FACTOR', True)
    F_0 = R.dup_zz_factor(f)

    config.setup('USE_CYCLOTOMIC_FACTOR', False)
    F_1 = R.dup_zz_factor(f)

    assert F_0 == F_1 == \
        (1, [(x**2 + 1, 1),
             (x**8 - x**6 + x**4 - x**2 + 1, 1)])

    config.setup('USE_CYCLOTOMIC_FACTOR')


def test_dmp_zz_wang():
    R,  x, y, z = ring("x,y,z", ZZ)
    UV, _x = ring("x", ZZ)

    p = ZZ(nextprime(R.dmp_zz_mignotte_bound(w_1)))
    assert p == 6291469

    t_1, k_1, e_1 = y, 1, ZZ(-14)
    t_2, k_2, e_2 = z, 2, ZZ(3)
    t_3, k_3, e_3 = y + z, 2, ZZ(-11)
    t_4, k_4, e_4 = y - z, 1, ZZ(-17)

    T = [t_1, t_2, t_3, t_4]
    K = [k_1, k_2, k_3, k_4]
    E = [e_1, e_2, e_3, e_4]

    T = zip([t.drop(x) for t in T], K)

    A = [ZZ(-14), ZZ(3)]

    S = R.dmp_eval_tail(w_1, A)
    cs, s = UV.dmp_ground_primitive(S)

    assert cs == 1 and s == S == \
        1036728*_x**6 + 915552*_x**5 + 55748*_x**4 + 105621*_x**3 - 17304*_x**2 - 26841*_x - 644

    assert R.dmp_zz_wang_non_divisors(E, cs, ZZ(4)) == [7, 3, 11, 17]
    assert s.is_squarefree and UV.dmp_degree_in(s, 0) == R.dmp_degree_in(w_1, 0)

    _, H = UV.dup_zz_factor_sqf(s)

    h_1 = 44*_x**2 + 42*_x + 1
    h_2 = 126*_x**2 - 9*_x + 28
    h_3 = 187*_x**2 - 23

    assert H == [h_1, h_2, h_3]

    LC = [lc.drop(x) for lc in [-4*y - 4*z, -y*z**2, y**2 - z**2]]

    assert R.dmp_zz_wang_lead_coeffs(w_1, T, cs, E, H, A) == (w_1, H, LC)

    # H_1 = [44*x**2 + 42*x + 1, 126*x**2 - 9*x + 28, 187*x**2 - 23]
    # H_2 = [-4*x**2*y - 12*x**2 - 3*x*y + 1, -9*x**2*y - 9*x - 2*y, x**2*y**2 - 9*x**2 + y - 9]
    # H_3 = [-4*x**2*y - 12*x**2 - 3*x*y + 1, -9*x**2*y - 9*x - 2*y, x**2*y**2 - 9*x**2 + y - 9]

    # c_1 = -70686*x**5 - 5863*x**4 - 17826*x**3 + 2009*x**2 + 5031*x + 74
    # c_2 = 9*x**5*y**4 + 12*x**5*y**3 - 45*x**5*y**2 - 108*x**5*y - 324*x**5 + 18*x**4*y**3 - 216*x**4*y**2 - 810*x**4*y + 2*x**3*y**4 + 9*x**3*y**3 - 252*x**3*y**2 - 288*x**3*y - 945*x**3 - 30*x**2*y**2 - 414*x**2*y + 2*x*y**3 - 54*x*y**2 - 3*x*y + 81*x + 12*y
    # c_3 = -36*x**4*y**2 - 108*x**4*y - 27*x**3*y**2 - 36*x**3*y - 108*x**3 - 8*x**2*y**2 - 42*x**2*y - 6*x*y**2 + 9*x + 2*y

    # TODO
    # assert R.dmp_zz_diophantine(H_1, c_1, [], 5, p) == [-3*x, -2, 1]
    # assert R.dmp_zz_diophantine(H_2, c_2, [ZZ(-14)], 5, p) == [-x*y, -3*x, -6]
    # assert R.dmp_zz_diophantine(H_3, c_3, [ZZ(-14)], 5, p) == [0, 0, -1]

    factors = R.dmp_zz_wang_hensel_lifting(w_1, H, LC, A, p)
    assert R.dmp_expand(factors) == w_1


def test_dmp_zz_diophantine():
    F = [[[[47], []], [[9, 0, 0, -9]]],
         [[[45]], [[]], [[]], [[-9], [-1], [], [3, 0, -6, 0]]]]
    c = [[[-270, 0, 0, 270]], [[]], [[94, 0], []],
         [[54, 0, 0, -54], [6, 0, 0, -6], [], [-18, 0, 54, 18, 0, -54, 0]]]
    A = [-2, 0]
    d = 6
    p = 2345258188817
    u = 2
    K = ZZ

    r = dmp_zz_diophantine(F, c, A, d, p, u, K)
    assert r == [[[[-6, 0, 0, 6]]], [[[2, 0]]]]


def test_sympyissue_6355():
    # This tests a bug in the Wang algorithm that occured only with a very
    # specific set of random numbers.
    random_sequence = [-1, -1, 0, 0, 0, 0, -1, -1, 0, -1, 3, -1, 3, 3, 3, 3, -1, 3]

    R, x, y, z = ring("x,y,z", ZZ)
    f = 2*x**2 + y*z - y - z**2 + z

    assert R.dmp_zz_wang(f, seed=random_sequence) == [f]


def test_dmp_zz_factor():
    R, x = ring("x", ZZ)
    assert R.dmp_zz_factor(0) == (0, [])
    assert R.dmp_zz_factor(7) == (7, [])
    assert R.dmp_zz_factor(-7) == (-7, [])

    assert R.dmp_zz_factor(x**2 - 9) == (1, [(x - 3, 1), (x + 3, 1)])

    R, x, y = ring("x,y", ZZ)
    assert R.dmp_zz_factor(0) == (0, [])
    assert R.dmp_zz_factor(7) == (7, [])
    assert R.dmp_zz_factor(-7) == (-7, [])

    assert R.dmp_zz_factor(x) == (1, [(x, 1)])
    assert R.dmp_zz_factor(4*x) == (4, [(x, 1)])
    assert R.dmp_zz_factor(4*x + 2) == (2, [(2*x + 1, 1)])
    assert R.dmp_zz_factor(x*y + 1) == (1, [(x*y + 1, 1)])
    assert R.dmp_zz_factor(y**2 + 1) == (1, [(y**2 + 1, 1)])
    assert R.dmp_zz_factor(y**2 - 1) == (1, [(y - 1, 1), (y + 1, 1)])

    assert R.dmp_zz_factor(x**2*y**2 + 6*x**2*y + 9*x**2 - 1) == (1, [(x*y + 3*x - 1, 1), (x*y + 3*x + 1, 1)])
    assert R.dmp_zz_factor(x**2*y**2 - 9) == (1, [(x*y - 3, 1), (x*y + 3, 1)])

    R, x, y, z = ring("x,y,z", ZZ)
    assert R.dmp_zz_factor(x**2*y**2*z**2 - 9) == \
        (1, [(x*y*z - 3, 1),
             (x*y*z + 3, 1)])

    R, x, y, z, u = ring("x,y,z,u", ZZ)
    assert R.dmp_zz_factor(x**2*y**2*z**2*u**2 - 9) == \
        (1, [(x*y*z*u - 3, 1),
             (x*y*z*u + 3, 1)])

    R, x, y, z = ring("x,y,z", ZZ)
    assert R.dmp_zz_factor(f_1) == \
        (1, [(x + y*z + 20, 1),
             (x*z + y + 30, 1),
             (x*y + z + 10, 1)])

    assert R.dmp_zz_factor(f_2) == \
        (1, [(x**2*y**2 + x**2*z**2 + y + 90, 1),
             (x**3*y + x**3*z + z - 11, 1)])

    assert R.dmp_zz_factor(f_3) == \
        (1, [(x**2*y**2 + x*z**4 + x + z, 1),
             (x**3 + x*y*z + y**2 + y*z**3, 1)])

    assert R.dmp_zz_factor(f_4) == \
        (-1, [(x*y**3 + z**2, 1),
              (x**2*z + y**4*z**2 + 5, 1),
              (x**3*y - z**2 - 3, 1),
              (x**3*y**4 + z**2, 1)])

    assert R.dmp_zz_factor(f_5) == \
        (-1, [(x + y - z, 3)])

    R, x, y, z, t = ring("x,y,z,t", ZZ)
    assert R.dmp_zz_factor(f_6) == \
        (1, [(47*x*y + z**3*t**2 - t**2, 1),
             (45*x**3 - 9*y**3 - y**2 + 3*z**3 + 2*z*t, 1)])

    R, x, y, z = ring("x,y,z", ZZ)
    assert R.dmp_zz_factor(w_1) == (1, [(4*x**2*y + 4*x**2*z + x*y*z - 1, 1),
                                        (x**2*y*z**2 + 3*x*z + 2*y, 1),
                                        (x**2*y**2 - x**2*z**2 + y - z**2, 1)])

    R, x, y = ring("x,y", ZZ)
    f = -12*x**16*y + 240*x**12*y**3 - 768*x**10*y**4 + 1080*x**8*y**5 - 768*x**6*y**6 + 240*x**4*y**7 - 12*y**9

    assert R.dmp_zz_factor(f) == \
        (-12, [(y, 1),
               (x**2 - y, 6),
               (x**4 + 6*x**2*y + y**2, 1)])


def test_dmp_ext_factor():
    R, x = ring("x", QQ.algebraic_field(I))

    assert R(0).factor_list() == (0, [])

    f = x + 1

    assert f.factor_list() == (1, [(f, 1)])

    g = 2*x + 2

    assert g.factor_list() == (2, [(f, 1)])

    f = 7*x**4 + 1
    g = x**4 + QQ(1, 7)

    assert f.factor_list() == (7, [(g, 1)])

    f = x**4 + 1

    assert f.factor_list() == (1, [(x**2 - I, 1), (x**2 + I, 1)])

    f = 4*x**2 + 9

    assert f.factor_list() == (4, [(x - 3*I/2, 1), (x + 3*I/2, 1)])

    f = 4*x**4 + 8*x**3 + 77*x**2 + 18*x + 153

    assert f.factor_list() == (4, [(x - 3*I/2, 1), (x + 1 + 4*I, 1),
                                   (x + 1 - 4*I, 1), (x + 3*I/2, 1)])

    R, x = ring("x", QQ.algebraic_field(sqrt(2)))

    f = x**4 + 1

    assert f.factor_list() == (1, [(x**2 - sqrt(2)*x + 1, 1),
                                   (x**2 + sqrt(2)*x + 1, 1)])

    f = x**2 + 2*sqrt(2)*x + 2

    assert f.factor_list() == (1, [(x + sqrt(2), 2)])
    assert (f**3).factor_list() == (1, [(x + sqrt(2), 6)])

    f *= 2

    assert f.factor_list() == (2, [(x + sqrt(2), 2)])
    assert (f**3).factor_list() == (8, [(x + sqrt(2), 6)])

    R,  x, y = ring("x,y", QQ.algebraic_field(sqrt(2)))

    assert R(0).factor_list() == (0, [])

    f = x + 1

    assert f.factor_list() == (1, [(f, 1)])

    g = 2*x + 2

    assert g.factor_list() == (2, [(f, 1)])

    f = x**2 - 2*y**2

    assert f.factor_list() == (1, [(x - sqrt(2)*y, 1), (x + sqrt(2)*y, 1)])

    f = 2*x**2 - 4*y**2

    assert f.factor_list() == (2, [(x - sqrt(2)*y, 1), (x + sqrt(2)*y, 1)])

    R,  x = ring("x", QQ.algebraic_field(I))
    f = x**2 + 1
    assert f.factor_list() == (1, [(x - I, 1), (x + I, 1)])


def test_sympyissue_5786():
    R,  x, y, z, t = ring("x, y, z, t", QQ.algebraic_field(I))

    f = (z - I*t)*(x - I*y)
    assert f.factor_list() == (1, [(z - I*t, 1), (x - I*y, 1)])

    f = (z - I*t)**2*(x - I*y)
    assert f.factor_list() == (1, [(z - I*t, 2), (x - I*y, 1)])

    f = (z - I*t)*(x - I*y)**3
    assert f.factor_list() == (1, [(z - I*t, 1), (x - I*y, 3)])


def test_dmp_factor_list():
    R, x = ring("x", ZZ)
    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    R, x = ring("x", QQ)
    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    R, x = ring("x", ZZ.poly_ring('t'))
    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    R, x = ring("x", QQ.poly_ring('t'))
    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    R, x = ring("x", ZZ)

    assert (x**2 + 2*x + 1).factor_list() == (1, [(x + 1, 2)])
    # issue sympy/sympy#8037
    assert (6*x**2 - 5*x - 6).factor_list() == (1, [(2*x - 3, 1), (3*x + 2, 1)])

    R, x = ring("x", QQ)
    assert (x**2/2 + x + QQ(1, 2)).factor_list() == (QQ(1, 2), [(x + 1, 2)])

    R, x = ring("x", FF(2))
    assert (x**2 + 1).factor_list() == (1, [(x + 1, 2)])

    R, x = ring("x", RR)
    assert (1.0*x**2 + 2.0*x + 1.0).factor_list() == (1.0, [(1.0*x + 1.0, 2)])
    assert (2.0*x**2 + 4.0*x + 2.0).factor_list() == (2.0, [(1.0*x + 1.0, 2)])

    f = 6.7225336055071*x**2 - 10.6463972754741*x - 0.33469524022264
    coeff, factors = f.factor_list()
    assert coeff == 1.0 and len(factors) == 1 and factors[0][0].almosteq(f, 1e-10) and factors[0][1] == 1

    # issue diofant/diofant#238
    f = 0.1*x**2 + 1.1*x + 1.0
    assert f.factor_list() == (10.0, [(0.1*x + 0.1, 1), (0.1*x + 1.0, 1)])
    f = 0.25 + 1.0*x + 1.0*x**2
    assert f.factor_list() == (4.0, [(0.25 + 0.5*x, 2)])

    Rt, t = ring("t", ZZ)
    R, x = ring("x", Rt)

    f = 4*t*x**2 + 4*t**2*x

    assert f.factor_list() == (4*t, [(x, 1), (x + t, 1)])

    Rt, t = ring("t", QQ)
    R, x = ring("x", Rt)

    f = t*x**2/2 + t**2*x/2

    assert f.factor_list() == (t/2, [(x, 1), (x + t, 1)])

    R, x = ring("x", QQ.algebraic_field(I))

    f = x**4 + 2*x**2

    assert f.factor_list() == (1, [(x, 2), (x**2 + 2, 1)])

    R, x = ring("x", EX)
    pytest.raises(DomainError, lambda: R(EX(sin(1))).factor_list())

    R, x, y = ring("x,y", ZZ)
    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    R, x, y = ring("x,y", QQ)
    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    Rt, t = ring("t", ZZ)
    R, x, y = ring("x,y", Rt)
    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    Rt, t = ring("t", QQ)
    R, x, y = ring("x,y", Rt)
    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    R, *X = ring("x:200", ZZ)

    f, g = X[0]**2 + 2*X[0] + 1, X[0] + 1
    assert f.factor_list() == (1, [(g, 2)])

    f, g = X[-1]**2 + 2*X[-1] + 1, X[-1] + 1
    assert f.factor_list() == (1, [(g, 2)])

    R, x = ring("x", ZZ)
    assert (x**2 + 2*x + 1).factor_list() == (1, [(x + 1, 2)])
    R, x = ring("x", QQ)
    assert (x**2/2 + x + QQ(1, 2)).factor_list() == (QQ(1, 2), [(x + 1, 2)])

    R, x, y = ring("x,y", ZZ)
    assert (x**2 + 2*x + 1).factor_list() == (1, [(x + 1, 2)])
    R, x, y = ring("x,y", QQ)
    assert (x**2/2 + x + QQ(1, 2)).factor_list() == (QQ(1, 2), [(x + 1, 2)])

    R, x, y = ring("x,y", ZZ)
    f = 4*x**2*y + 4*x*y**2

    assert f.factor_list() == (4, [(y, 1), (x, 1), (x + y, 1)])

    R,  x, y = ring("x,y", QQ)
    f = x**2*y/2 + x*y**2/2

    assert f.factor_list() == (QQ(1, 2), [(y, 1), (x, 1), (x + y, 1)])

    R,  x, y = ring("x,y", RR)
    f = 2.0*x**2 - 8.0*y**2

    assert f.factor_list() == (2.0, [(1.0*x - 2.0*y, 1), (1.0*x + 2.0*y, 1)])

    f = 6.7225336055071*x**2*y**2 - 10.6463972754741*x*y - 0.33469524022264
    coeff, factors = f.factor_list()
    assert coeff == 1.0 and len(factors) == 1 and factors[0][0].almosteq(f, 1e-10) and factors[0][1] == 1

    # issue diofant/diofant#238
    R,  x, y, z = ring("x,y,z", RR)
    f = x*y + x*z + 0.1*y + 0.1*z
    assert f.factor_list() == (10.0, [(0.1*y + 0.1*z, 1), (x + 0.1, 1)])
    f = 0.25*x**2 + 1.0*x*y*z + 1.0*y**2*z**2
    assert f.factor_list() == (4.0, [(0.25*x + 0.5*y*z, 2)])

    Rt, t = ring("t", ZZ)
    R, x, y = ring("x,y", Rt)
    f = 4*t*x**2 + 4*t**2*x

    assert f.factor_list() == (4*t, [(x, 1), (x + t, 1)])

    Rt, t = ring("t", QQ)
    R, x, y = ring("x,y", Rt)
    f = t*x**2/2 + t**2*x/2

    assert f.factor_list() == (t/2, [(x, 1), (x + t, 1)])

    R, x, y = ring("x,y", FF(2))
    pytest.raises(NotImplementedError, lambda: (x**2 + y**2).factor_list())

    R, x, y = ring("x,y", EX)
    pytest.raises(DomainError, lambda: R(EX(sin(1))).factor_list())

    R, x, y = ring('x, y', QQ.algebraic_field(I))
    f, r = x**2 + y**2, (1, [(x - I*y, 1), (x + I*y, 1)])

    assert R.dmp_factor_list(f) == r

    with config.using(aa_factor_method='trager'):
        assert R.dmp_factor_list(f) == r


def test_gf_factor():
    R, x = ring('x', FF(11))

    assert R(0).factor_list() == (0, [])
    assert R(1).factor_list() == (1, [])
    assert (x + 1).factor_list() == (1, [(x + 1, 1)])

    assert (5*x**3 + 2*x**2 + 7*x +
            2).factor_list() == (5, [(x + 2, 1), (x + 8, 2)])

    f = x**6 + 8*x**5 + x**4 + 8*x**3 + 10*x**2 + 8*x + 1
    g = (1, [(x + 1, 1),
             (x**2 + 5*x + 3, 1),
             (x**3 + 2*x**2 + 3*x + 4, 1)])

    with config.using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    f = x**3 + 5*x**2 + 8*x + 4

    g = (1, [(x + 1, 1), (x + 2, 2)])

    with config.using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    f = x**9 + x**8 + 10*x**7 + x**6 + 10*x**4 + 10*x**3 + 10*x**2
    g = (1, [(x, 2), (x**2 + 9*x + 5, 1),
             (x**5 + 3*x**4 + 8*x**2 + 5*x + 2, 1)])

    with config.using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    f = x**32 + 1

    g = (1, [(x**16 + 3*x**8 + 10, 1),
             (x**16 + 8*x**8 + 10, 1)])

    with config.using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    f = 8*x**32 + 5
    g = (8, [(x + 3, 1),
             (x + 8, 1),
             (x**2 + 9, 1),
             (x**2 + 2*x + 2, 1),
             (x**2 + 9*x + 2, 1),
             (x**4 + 5*x**2 + 7, 1),
             (x**4 + 6*x**2 + 7, 1),
             (x**8 + x**4 + 6, 1),
             (x**8 + 10*x**4 + 6, 1)])

    with config.using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    f = 8*x**63 + 5
    g = (8, [(x + 7, 1),
             (x**2 + 4*x + 5, 1),
             (x**3 + 6*x**2 + 8*x + 2, 1),
             (x**3 + 9*x**2 + 9*x + 2, 1),
             (x**6 + 9*x**3 + 4, 1),
             (x**6 + 2*x**5 + 8*x**3 + 4*x**2 + 6*x + 4, 1),
             (x**6 + 2*x**5 + 3*x**4 + 8*x**3 + 6*x + 4, 1),
             (x**6 + 2*x**5 + 6*x**4 + 8*x**2 + 4*x + 4, 1),
             (x**6 + 3*x**5 + 3*x**4 + x**3 + 6*x**2 + 8*x + 4, 1),
             (x**6 + 5*x**5 + 6*x**4 + 8*x**2 + 6*x + 4, 1),
             (x**6 + 6*x**5 + 2*x**4 + 7*x**3 + 9*x**2 + 8*x + 4, 1),
             (x**6 + 10*x**5 + 4*x**4 + 7*x**3 + 10*x**2 + 7*x + 4, 1),
             (x**6 + 10*x**5 + 10*x**4 + x**3 + 4*x**2 + 9*x + 4, 1)])

    with config.using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='other'):
        pytest.raises(KeyError, lambda: (x + 1).factor_list())

    R, x = ring('x', FF(2))

    f = x**4 + x
    g = (1, [(x, 1),
             (x + 1, 1),
             (x**2 + x + 1, 1)])

    with config.using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    # Gathen polynomials: x**n + x + 1 (mod p > 2**n * pi)

    p = ZZ(nextprime(int((2**15*pi))))
    R, x = ring('x', FF(p))

    f = x**15 + x + 1
    g = (1, [(x**2 + 22730*x + 68144, 1),
             (x**4 + 81553*x**3 + 77449*x**2 + 86810*x + 4724, 1),
             (x**4 + 86276*x**3 + 56779*x**2 + 14859*x + 31575, 1),
             (x**5 + 15347*x**4 + 95022*x**3 + 84569*x**2 + 94508*x + 92335, 1)])

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g

    # Shoup polynomials: f = a_0 x**n + a_1 x**(n-1) + ... + a_n
    # (mod p > 2**(n-2) * pi), where a_n = a_{n-1}**2 + 1, a_0 = 1

    p = ZZ(nextprime(int((2**4*pi))))
    R, x = ring('x', FF(p))

    f = x**6 + 2*x**5 + 5*x**4 + 26*x**3 + 41*x**2 + 39*x + 38

    g = (1, [(x**2 + 44*x + 26, 1),
             (x**4 + 11*x**3 + 25*x**2 + 18*x + 30, 1)])

    with config.using(gf_factor_method='zassenhaus'):
        assert f.factor_list() == g

    with config.using(gf_factor_method='shoup'):
        assert f.factor_list() == g


def test_PolyElement_is_irreducible():
    R, x = ring("x", ZZ)

    assert (x**2 + x + 1).is_irreducible is True
    assert (x**2 + 2*x + 1).is_irreducible is False
    assert (x**2 - 1).is_irreducible is False

    R, x, y = ring("x,y", ZZ)

    assert R(2).is_irreducible is True
    assert (x**2 + x + 1).is_irreducible is True
    assert (x**2 + 2*x + 1).is_irreducible is False
    assert ((x - 2*y)*(x + y)).is_irreducible is False
    assert (x**2 + y**2).is_irreducible is True
