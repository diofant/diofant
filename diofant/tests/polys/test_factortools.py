"""Tools for polynomial factorization routines in characteristic zero."""

import functools
import operator

import pytest

from diofant import (EX, FF, QQ, RR, ZZ, DomainError, ExtraneousFactors, I,
                     nextprime, pi, ring, sin, sqrt)
from diofant.config import using
from diofant.polys.specialpolys import f_polys, w_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = f_polys()
w_1, w_2 = w_polys()


def test__trial_division():
    R, x = ring('x', ZZ)

    assert R._trial_division(x**5 + 8*x**4 + 25*x**3 + 38*x**2 + 28*x +
                             8, (x + 1, x + 2)) == [(x + 1, 2), (x + 2, 3)]

    R, x, y = ring('x y', ZZ)

    assert R._trial_division(x**5 + 8*x**4 + 25*x**3 + 38*x**2 + 28*x +
                             8, (x + 1, x + 2)) == [(x + 1, 2), (x + 2, 3)]


def test__zz_mignotte_bound():
    R, x = ring('x', ZZ)

    assert R._zz_mignotte_bound(2*x**2 + 3*x + 4) == 32

    R, x, y = ring('x y', ZZ)

    assert R._zz_mignotte_bound(2*x**2 + 3*x + 4) == 32


def test__zz_hensel_step():
    R, x = ring('x', ZZ)

    f = x**4 - 1
    g = x**3 + 2*x**2 - x - 2
    h = x - 2
    s = R(-2)
    t = 2*x**2 - 2*x - 1

    G, H, S, T = R._zz_hensel_step(5, f, g, h, s, t)

    assert G == x**3 + 7*x**2 - x - 7
    assert H == x - 7
    assert S == 8
    assert T == -8*x**2 - 12*x - 1


def test__zz_hensel_lift():
    R, x = ring('x', ZZ)

    f = x**4 - 1
    F = [x - 1, x - 2, x + 2, x + 1]

    assert R._zz_hensel_lift(ZZ(5), f, F, 4) == [x - 1, x - 182,
                                                 x + 182, x + 1]


def test__cyclotomic_p():
    R, x = ring('x', ZZ)

    assert (x - 1).is_cyclotomic is True
    assert (x + 1).is_cyclotomic is True
    assert (x**2 + x + 1).is_cyclotomic is True

    f = x**2 + 1

    assert f.is_cyclotomic is True
    assert R._cyclotomic_p(f, irreducible=True) is True

    assert (x**4 + x**3 + x**2 + x + 1).is_cyclotomic is True
    assert (x**2 - x + 1).is_cyclotomic is True
    assert (x**6 + x**5 + x**4 + x**3 + x**2 + x + 1).is_cyclotomic is True
    assert (x**4 + 1).is_cyclotomic is True
    assert (x**6 + x**3 + 1).is_cyclotomic is True

    assert R(0).is_cyclotomic is False
    assert R(1).is_cyclotomic is False
    assert x.is_cyclotomic is False
    assert (x + 2).is_cyclotomic is False
    assert (3*x + 1).is_cyclotomic is False
    assert (x**2 - 1).is_cyclotomic is False

    f = x**16 + x**14 - x**10 - x**6 + x**2 + 1

    assert (f + x**8).is_cyclotomic is False
    assert (f - x**8).is_cyclotomic is True

    R, x = ring('x', QQ)

    assert (x**2 + x + 1).is_cyclotomic is True
    assert (x**2/2 + x + 1).is_cyclotomic is False

    R, x = ring('x', ZZ.inject('y'))

    assert (x**2 + x + 1).is_cyclotomic is False


def test__zz_cyclotomic_poly():
    R, x = ring('x', ZZ)

    assert R._zz_cyclotomic_poly(1) == x - 1
    assert R._zz_cyclotomic_poly(2) == x + 1
    assert R._zz_cyclotomic_poly(3) == x**2 + x + 1
    assert R._zz_cyclotomic_poly(4) == x**2 + 1
    assert R._zz_cyclotomic_poly(5) == x**4 + x**3 + x**2 + x + 1
    assert R._zz_cyclotomic_poly(6) == x**2 - x + 1
    assert R._zz_cyclotomic_poly(7) == (x**6 + x**5 + x**4 + x**3 +
                                        x**2 + x + 1)
    assert R._zz_cyclotomic_poly(8) == x**4 + 1
    assert R._zz_cyclotomic_poly(9) == x**6 + x**3 + 1


def test__zz_cyclotomic_factor():
    R, x = ring('x', ZZ)

    assert R._zz_cyclotomic_factor(R(0)) is None
    assert R._zz_cyclotomic_factor(R(1)) is None

    assert R._zz_cyclotomic_factor(2*x**10 - 1) is None
    assert R._zz_cyclotomic_factor(x**10 - 3) is None
    assert R._zz_cyclotomic_factor(x**10 + x**5 - 1) is None

    assert R._zz_cyclotomic_factor(x + 1) == [x + 1]
    assert R._zz_cyclotomic_factor(x - 1) == [x - 1]

    assert R._zz_cyclotomic_factor(x**2 + 1) == [x**2 + 1]
    assert R._zz_cyclotomic_factor(x**2 - 1) == [x - 1, x + 1]

    assert R._zz_cyclotomic_factor(x**27 + 1) == [x + 1, x**2 - x + 1,
                                                  x**6 - x**3 + 1,
                                                  x**18 - x**9 + 1]
    assert R._zz_cyclotomic_factor(x**27 - 1) == [x - 1, x**2 + x + 1,
                                                  x**6 + x**3 + 1,
                                                  x**18 + x**9 + 1]


def test_dup_zz_factor():
    R, x = ring('x', ZZ)

    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])
    assert R(-7).factor_list() == (-7, [])

    assert R._zz_factor_sqf(R(0)) == (0, [])
    assert R._zz_factor_sqf(R(7)) == (7, [])
    assert R._zz_factor_sqf(R(-7)) == (-7, [])

    assert (2*x + 4).factor_list() == (2, [(x + 2, 1)])
    assert R._zz_factor_sqf(2*x + 4) == (2, [x + 2])

    f = x**4 + x + 1

    for i in range(20):
        assert f.factor_list() == (1, [(f, 1)])

    f = x**5 - x**3 - x**2 + 1

    assert f.factor_list() == (1, [(x + 1, 1), (x - 1, 2), (x**2 + x + 1, 1)])

    for test in (True, False):
        with using(use_irreducible_in_factor=test):
            assert (x**2 + 2*x + 2).factor_list() == (1, [(x**2 + 2*x + 2, 1)])

            assert (18*x**2 + 12*x + 2).factor_list() == (2, [(3*x + 1, 2)])

            f = -9*x**2 + 1

            assert R._zz_factor_sqf(f) == (-1, [3*x - 1, 3*x + 1])
            assert f.factor_list() == (-1, [(3*x - 1, 1), (3*x + 1, 1)])

            assert R._zz_factor_sqf(3*x**4 + 2*x**3 + 6*x**2 +
                                    8*x + 10) == (1, [3*x**4 + 2*x**3 +
                                                      6*x**2 + 8*x + 10])

    with using(use_cyclotomic_factor=False):
        assert R._zz_factor_sqf(-9*x**2 + 1) == (-1, [3*x - 1, 3*x + 1])

    assert (x**3 - 6*x**2 + 11*x - 6).factor_list() == (1, [(x - 3, 1),
                                                            (x - 2, 1),
                                                            (x - 1, 1)])

    assert R._zz_factor_sqf(x**3 - 6*x**2 + 11*x - 6) == (1, [x - 3, x - 2,
                                                              x - 1])

    assert (3*x**3 + 10*x**2 + 13*x +
            10).factor_list() == (1, [(x + 2, 1), (3*x**2 + 4*x + 5, 1)])

    assert R._zz_factor_sqf(3*x**3 + 10*x**2 +
                            13*x + 10) == (1, [x + 2, 3*x**2 + 4*x + 5])

    assert (-x**6 + x**2).factor_list() == (-1, [(x, 2), (x - 1, 1), (x + 1, 1),
                                                 (x**2 + 1, 1)])

    f = (1080*x**8 + 5184*x**7 + 2099*x**6 + 744*x**5 + 2736*x**4 -
         648*x**3 + 129*x**2 - 324)

    assert f.factor_list() == (1, [(216*x**4 + 31*x**2 - 27, 1),
                                   (5*x**4 + 24*x**3 + 9*x**2 + 12, 1)])

    f = (-29802322387695312500000000000000000000*x**25 +
         2980232238769531250000000000000000*x**20 +
         1743435859680175781250000000000*x**15 +
         114142894744873046875000000*x**10 - 210106372833251953125*x**5 +
         + 95367431640625)

    assert (f.factor_list() ==
            (-95367431640625,
             [(5*x - 1, 1), (100*x**2 + 10*x - 1, 2),
              (625*x**4 + 125*x**3 + 25*x**2 + 5*x + 1, 1),
              (10000*x**4 - 3000*x**3 + 400*x**2 - 20*x + 1, 2),
              (10000*x**4 + 2000*x**3 + 400*x**2 + 30*x + 1, 2)]))

    f = x**10 - 1

    for test in (True, False):
        with using(use_cyclotomic_factor=test):
            f = x**10 - 1

            assert f.factor_list() == (1, [(x - 1, 1), (x + 1, 1),
                                           (x**4 - x**3 + x**2 - x + 1, 1),
                                           (x**4 + x**3 + x**2 + x + 1, 1)])

            f = x**10 + 1

            assert f.factor_list() == (1, [(x**2 + 1, 1),
                                           (x**8 - x**6 + x**4 - x**2 + 1, 1)])


def test__zz_wang():
    R, x, y, z = ring('x y z', ZZ)
    UV, _x = ring('x', ZZ)

    p = ZZ(nextprime(R._zz_mignotte_bound(w_1)))

    assert p == 6291469

    t_1, k_1, e_1 = y, 1, ZZ(-14)
    t_2, k_2, e_2 = z, 2, ZZ(3)
    t_3, k_3, e_3 = y + z, 2, ZZ(-11)
    t_4, k_4, e_4 = y - z, 1, ZZ(-17)

    T = [t_1, t_2, t_3, t_4]
    K = [k_1, k_2, k_3, k_4]
    E = [e_1, e_2, e_3, e_4]

    T = list(zip([t.drop(x) for t in T], K))

    A = [ZZ(-14), ZZ(3)]

    S = w_1.eval([(y, A[0]), (z, A[1])])
    cs, s = S.primitive()

    assert cs == 1 and s == S == (1036728*_x**6 + 915552*_x**5 + 55748*_x**4 +
                                  105621*_x**3 - 17304*_x**2 - 26841*_x - 644)

    assert R._zz_wang_non_divisors(E, cs, ZZ(4)) == [7, 3, 11, 17]
    assert s.is_squarefree and s.degree() == w_1.degree()

    _, H = UV._zz_factor_sqf(s)

    h_1 = 187*_x**2 - 23
    h_2 = 44*_x**2 + 42*_x + 1
    h_3 = 126*_x**2 - 9*_x + 28

    LC = [lc.drop(x) for lc in [y**2 - z**2, -4*y - 4*z, -y*z**2]]
    factors = R._zz_wang_hensel_lifting(w_1, H, LC, A, p)

    assert H == [h_1, h_2, h_3]
    assert R._zz_wang_lead_coeffs(w_1, T, cs, E, H, A) == (w_1, H, LC)
    assert functools.reduce(operator.mul, factors) == w_1

    # coverage tests
    f = x**6 + 5*x**4*y - 5*x**2*y**2 - y**3

    assert R._zz_wang(f, mod=4, seed=1) == [x**2 - y, x**4 + 6*x**2*y + y**2]

    # This tests a bug in the Wang algorithm that occured only with a very
    # specific set of random numbers; issue sympy/sympy#6355.
    random_sequence = [-1, -1, 0, 0, 0, 0, -1, -1, 0, -1, 3, -1, 3, 3, 3,
                       3, -1, 3]

    R, x, y, z = ring('x y z', ZZ)

    f = 2*x**2 + y*z - y - z**2 + z

    assert R._zz_wang(f, seed=random_sequence) == [f]

    with using(eez_restart_if_needed=False):
        pytest.raises(ExtraneousFactors,
                      lambda: R._zz_wang(f, seed=random_sequence))


def test__zz_diophantine():
    R, x, y = ring('x y', ZZ)

    H_1 = [44*x**2 + 42*x + 1, 126*x**2 - 9*x + 28, 187*x**2 - 23]
    H_2 = [-4*x**2*y - 12*x**2 - 3*x*y + 1, -9*x**2*y - 9*x - 2*y,
           x**2*y**2 - 9*x**2 + y - 9]
    H_3 = [-4*x**2*y - 12*x**2 - 3*x*y + 1, -9*x**2*y - 9*x - 2*y,
           x**2*y**2 - 9*x**2 + y - 9]
    c_1 = -70686*x**5 - 5863*x**4 - 17826*x**3 + 2009*x**2 + 5031*x + 74
    c_2 = (9*x**5*y**4 + 12*x**5*y**3 - 45*x**5*y**2 - 108*x**5*y -
           324*x**5 + 18*x**4*y**3 - 216*x**4*y**2 - 810*x**4*y +
           2*x**3*y**4 + 9*x**3*y**3 - 252*x**3*y**2 - 288*x**3*y -
           945*x**3 - 30*x**2*y**2 - 414*x**2*y + 2*x*y**3 -
           54*x*y**2 - 3*x*y + 81*x + 12*y)
    c_3 = (-36*x**4*y**2 - 108*x**4*y - 27*x**3*y**2 - 36*x**3*y -
           108*x**3 - 8*x**2*y**2 - 42*x**2*y - 6*x*y**2 + 9*x + 2*y)
    p = 6291469

    assert R._zz_diophantine(H_1, c_1, [ZZ(0)], 5, p) == [-3*x, -2, 1]
    assert R._zz_diophantine(H_2, c_2, [ZZ(-14)], 5, p) == [-x*y, -3*x, -6]
    assert R._zz_diophantine(H_3, c_3, [ZZ(-14)], 5, p) == [0, 0, -1]

    R, x, y, z = ring('x y z', ZZ)

    F = [47*x*y + 9*z**3 - 9, 45*x**3 - 9*y**3 - y**2 + 3*z**3 - 6*z]
    c = (-270*x**3*z**3 + 270*x**3 + 94*x*y*z + 54*y**3*z**3 - 54*y**3 +
         6*y**2*z**3 - 6*y**2 - 18*z**6 + 54*z**4 + 18*z**3 - 54*z)
    p = 2345258188817

    assert R._zz_diophantine(F, c, [ZZ(-2), ZZ(0)], 6,
                             p) == [-6*z**3 + 6, 2*z]


def test_dmp_zz_factor():
    R, x = ring('x', ZZ)

    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])
    assert R(-7).factor_list() == (-7, [])

    assert (x**2 - 9).factor_list() == (1, [(x - 3, 1), (x + 3, 1)])

    R, x, y = ring('x y', ZZ)

    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])
    assert R(-7).factor_list() == (-7, [])

    assert x.factor_list() == (1, [(x, 1)])
    assert (4*x).factor_list() == (4, [(x, 1)])
    assert (4*x + 2).factor_list() == (2, [(2*x + 1, 1)])
    assert (x*y + 1).factor_list() == (1, [(x*y + 1, 1)])
    assert (y**2 + 1).factor_list() == (1, [(y**2 + 1, 1)])
    assert (y**2 - 1).factor_list() == (1, [(y - 1, 1), (y + 1, 1)])

    assert (x**2*y**2 + 6*x**2*y +
            9*x**2 - 1).factor_list() == (1, [(x*y + 3*x - 1, 1),
                                              (x*y + 3*x + 1, 1)])
    assert (x**2*y**2 - 9).factor_list() == (1, [(x*y - 3, 1), (x*y + 3, 1)])

    f = (-12*x**16*y + 240*x**12*y**3 - 768*x**10*y**4 + 1080*x**8*y**5 -
         768*x**6*y**6 + 240*x**4*y**7 - 12*y**9)

    assert f.factor_list() == (-12, [(y, 1), (x**2 - y, 6),
                                     (x**4 + 6*x**2*y + y**2, 1)])

    R, x, y, z = ring('x y z', ZZ)

    assert (x**2*y**2*z**2 - 9).factor_list() == (1, [(x*y*z - 3, 1),
                                                      (x*y*z + 3, 1)])

    assert f_1.factor_list() == (1, [(x*y + z + 10, 1), (x + y*z + 20, 1),
                                     (x*z + y + 30, 1)])

    assert f_2.factor_list() == (1, [(x**3*y + x**3*z + z - 11, 1),
                                     (x**2*y**2 + x**2*z**2 + y + 90, 1)])

    assert f_3.factor_list() == (1, [(x**2*y**2 + x*z**4 + x + z, 1),
                                     (x**3 + x*y*z + y**2 + y*z**3, 1)])

    assert f_4.factor_list() == (-1, [(x*y**3 + z**2, 1), (x**3*y**4 + z**2, 1),
                                      (x**3*y - z**2 - 3, 1),
                                      (x**2*z + y**4*z**2 + 5, 1)])

    assert f_5.factor_list() == (-1, [(x + y - z, 3)])

    assert w_1.factor_list() == (1, [(x**2*y*z**2 + 3*x*z + 2*y, 1),
                                     (4*x**2*y + 4*x**2*z + x*y*z - 1, 1),
                                     (x**2*y**2 - x**2*z**2 + y - z**2, 1)])

    R, x, y, z, t = ring('x y z t', ZZ)

    assert (x**2*y**2*z**2*t**2 - 9).factor_list() == (1, [(x*y*z*t - 3, 1),
                                                           (x*y*z*t + 3, 1)])

    assert f_6.factor_list() == (1, [(47*x*y + z**3*t**2 - t**2, 1),
                                     (45*x**3 - 9*y**3 - y**2 + 3*z**3 +
                                      2*z*t, 1)])


@pytest.mark.parametrize('method', ('modular', 'trager'))
def test_dmp_ext_factor(method):
    with using(aa_factor_method=method):
        R, x = ring('x', QQ.algebraic_field(I))

        assert R(0).factor_list() == (0, [])
        assert (x + 1).factor_list() == (1, [(x + 1, 1)])
        assert (2*x + 2).factor_list() == (2, [(x + 1, 1)])
        assert (7*x**4 + 1).factor_list() == (7, [(x**4 + QQ(1, 7), 1)])
        assert (x**4 + 1).factor_list() == (1, [(x**2 - I, 1), (x**2 + I, 1)])
        assert (4*x**2 + 9).factor_list() == (4, [(x - 3*I/2, 1), (x + 3*I/2, 1)])
        assert (4*x**4 + 8*x**3 + 77*x**2 + 18*x +
                153).factor_list() == (4, [(x - 3*I/2, 1), (x + 1 + 4*I, 1),
                                           (x + 1 - 4*I, 1), (x + 3*I/2, 1)])
        assert (x**2 + 1).factor_list() == (1, [(x - I, 1), (x + I, 1)])

        R, x = ring('x', QQ.algebraic_field(sqrt(2)))

        assert (x**4 + 1).factor_list() == (1, [(x**2 - sqrt(2)*x + 1, 1),
                                                (x**2 + sqrt(2)*x + 1, 1)])

        f = x**2 + 2*sqrt(2)*x + 2

        assert f.factor_list() == (1, [(x + sqrt(2), 2)])
        assert (f**3).factor_list() == (1, [(x + sqrt(2), 6)])

        f *= 2

        assert f.factor_list() == (2, [(x + sqrt(2), 2)])
        assert (f**3).factor_list() == (8, [(x + sqrt(2), 6)])

        R, x, y = ring('x y', QQ.algebraic_field(sqrt(2)))

        assert R(0).factor_list() == (0, [])
        assert (x + 1).factor_list() == (1, [(x + 1, 1)])
        assert (2*x + 2).factor_list() == (2, [(x + 1, 1)])
        assert (x**2 - 2*y**2).factor_list() == (1, [(x - sqrt(2)*y, 1),
                                                     (x + sqrt(2)*y, 1)])
        assert (2*x**2 - 4*y**2).factor_list() == (2, [(x - sqrt(2)*y, 1),
                                                       (x + sqrt(2)*y, 1)])


def test_sympyissue_5786():
    R, x, y, z, t = ring('x y z t', QQ.algebraic_field(I))

    f, g = z - I*t, x - I*y

    assert (f*g).factor_list() == (1, [(f, 1), (g, 1)])
    assert (f**2*g).factor_list() == (1, [(g, 1), (f, 2)])
    assert (f*g**3).factor_list() == (1, [(f, 1), (g, 3)])


def test_factor_list():
    R, x = ring('x', FF(2))

    assert (x**2 + 1).factor_list() == (1, [(x + 1, 2)])

    R, x = ring('x', ZZ)

    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    assert (x**2 + 2*x + 1).factor_list() == (1, [(x + 1, 2)])

    # issue sympy/sympy#8037
    assert (6*x**2 - 5*x - 6).factor_list() == (1, [(2*x - 3, 1), (3*x + 2, 1)])

    R, x = ring('x', QQ)

    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    assert (x**2/2 + x + QQ(1, 2)).factor_list() == (QQ(1, 2), [(x + 1, 2)])

    R, x = ring('x', QQ.algebraic_field(I))

    f = x**4 + 2*x**2

    assert f.factor_list() == (1, [(x, 2), (x**2 + 2, 1)])

    R, x = ring('x', RR)

    assert (1.0*x**2 + 2.0*x + 1.0).factor_list() == (1.0, [(1.0*x + 1.0, 2)])
    assert (2.0*x**2 + 4.0*x + 2.0).factor_list() == (2.0, [(1.0*x + 1.0, 2)])

    f = 6.7225336055071*x**2 - 10.6463972754741*x - 0.33469524022264

    assert f.factor_list() == (1.0, [(f, 1)])

    # issue diofant/diofant#238
    f = 0.1*x**2 + 1.1*x + 1.0

    assert f.factor_list() == (10.0, [(0.1*x + 0.1, 1), (0.1*x + 1.0, 1)])

    f = 0.25 + 1.0*x + 1.0*x**2

    assert f.factor_list() == (4.0, [(0.25 + 0.5*x, 2)])

    Rt, t = ring('t', ZZ)
    R, x = ring('x', Rt)

    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    assert (4*t*x**2 + 4*t**2*x).factor_list() == (4*t, [(x, 1), (x + t, 1)])

    Rt, t = ring('t', QQ)
    R, x = ring('x', Rt)

    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    assert (t*x**2/2 + t**2*x/2).factor_list() == (t/2, [(x, 1), (x + t, 1)])

    R, x = ring('x', EX)

    pytest.raises(DomainError, lambda: R(EX(sin(1))).factor_list())

    R, x, y = ring('x y', FF(2))

    pytest.raises(NotImplementedError, lambda: (x**2 + y**2).factor_list())

    R, x, y = ring('x y', ZZ)

    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    assert (x**2 + 2*x + 1).factor_list() == (1, [(x + 1, 2)])
    assert (4*x**2*y + 4*x*y**2).factor_list() == (4, [(y, 1), (x, 1),
                                                       (x + y, 1)])

    R, x, y = ring('x y', QQ)

    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    assert (x**2/2 + x + QQ(1, 2)).factor_list() == (QQ(1, 2), [(x + 1, 2)])
    assert (x**2*y/2 + x*y**2/2).factor_list() == (QQ(1, 2), [(y, 1), (x, 1),
                                                              (x + y, 1)])

    R, x, y = ring('x y', QQ.algebraic_field(I))

    f, r = x**2 + y**2, (1, [(x - I*y, 1), (x + I*y, 1)])

    for method in ('trager', 'modular'):
        with using(aa_factor_method=method):
            assert f.factor_list() == r

    R, x, y = ring('x y', RR)

    f = 2.0*x**2 - 8.0*y**2

    assert f.factor_list() == (2.0, [(1.0*x - 2.0*y, 1), (1.0*x + 2.0*y, 1)])

    f = 6.7225336055071*x**2*y**2 - 10.6463972754741*x*y - 0.33469524022264

    assert f.factor_list() == (1.0, [(f, 1)])

    Rt, t = ring('t', ZZ)
    R, x, y = ring('x y', Rt)

    assert R(0).factor_list() == (0, [])
    assert R(7).factor_list() == (7, [])

    assert (4*t*x**2 + 4*t**2*x).factor_list() == (4*t, [(x, 1), (x + t, 1)])

    Rt, t = ring('t', QQ)
    R, x, y = ring('x y', Rt)

    assert R(0).factor_list() == (0, [])
    assert R(QQ(1, 7)).factor_list() == (QQ(1, 7), [])

    assert (t*x**2/2 + t**2*x/2).factor_list() == (t/2, [(x, 1), (x + t, 1)])

    R, x, y = ring('x y', EX)

    pytest.raises(DomainError, lambda: R(EX(sin(1))).factor_list())

    # issue diofant/diofant#238
    R, x, y, z = ring('x y z', RR)

    f = x*y + x*z + 0.1*y + 0.1*z

    assert f.factor_list() == (10.0, [(x + 0.1, 1), (0.1*y + 0.1*z, 1)])

    f = 0.25*x**2 + 1.0*x*y*z + 1.0*y**2*z**2

    assert f.factor_list() == (4.0, [(0.25*x + 0.5*y*z, 2)])

    R, *X = ring('x:200', ZZ)

    f, g = X[0]**2 + 2*X[0] + 1, X[0] + 1

    assert f.factor_list() == (1, [(g, 2)])

    f, g = X[-1]**2 + 2*X[-1] + 1, X[-1] + 1

    assert f.factor_list() == (1, [(g, 2)])


def test_gf_factor():
    R, x = ring('x', FF(2))

    f = x**4 + x
    g = (1, [(x, 1),
             (x + 1, 1),
             (x**2 + x + 1, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**18 + x**17 + x**16 + x**14 + x**12 + x**11 + x**8 + x**5 + x**3 + 1
    g = (1, [(x + 1, 4), (x**4 + x**3 + 1, 1),
             (x**10 + x**8 + x**7 + x**5 + 1, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**63 + 1
    g = (1, [(x + 1, 1), (x**2 + x + 1, 1), (x**3 + x + 1, 1),
             (x**6 + x + 1, 1), (x**3 + x**2 + 1, 1),
             (x**6 + x**3 + 1, 1), (x**6 + x**5 + 1, 1),
             (x**6 + x**4 + x**2 + x + 1, 1), (x**6 + x**5 + x**2 + x + 1, 1),
             (x**6 + x**4 + x**3 + x + 1, 1), (x**6 + x**5 + x**4 + x + 1, 1),
             (x**6 + x**5 + x**3 + x**2 + 1, 1),
             (x**6 + x**5 + x**4 + x**2 + 1, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = (x**28 + x**27 + x**26 + x**25 + x**24 + x**20 + x**19 + x**17 +
         x**16 + x**15 + x**14 + x**13 + x**12 + x**11 + x**9 + x**8 +
         x**5 + x**4 + x**2 + x)
    g = (1, [(x, 1), (x + 1, 2), (x**5 + x**4 + x**3 + x + 1, 1),
             (x**10 + x**9 + x**8 + x**7 + 1, 1),
             (x**10 + x**9 + x**8 + x**5 + x**4 + x**2 + 1, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    R, x = ring('x', FF(3))

    f = x**6 - x**5 + x**4 + x**3 - x
    g = (1, [(x, 1), (x + 1, 1), (x**2 + 1, 1), (x**2 + x + 2, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**4 + x**3 + x + 2
    g = (1, [(x**2 + 1, 1), (x**2 + x + 2, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    R, x = ring('x', FF(11))

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert R(0).factor_list() == (0, [])
            assert R(1).factor_list() == (1, [])
            assert x.factor_list() == (1, [(x, 1)])
            assert (x + 1).factor_list() == (1, [(x + 1, 1)])
            assert (2*x + 3).factor_list() == (2, [(x + 7, 1)])

    assert (5*x**3 + 2*x**2 + 7*x +
            2).factor_list() == (5, [(x + 2, 1), (x + 8, 2)])

    f = x**6 + 8*x**5 + x**4 + 8*x**3 + 10*x**2 + 8*x + 1
    g = (1, [(x + 1, 1),
             (x**2 + 5*x + 3, 1),
             (x**3 + 2*x**2 + 3*x + 4, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**3 + 5*x**2 + 8*x + 4
    g = (1, [(x + 1, 1), (x + 2, 2)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**9 + x**8 + 10*x**7 + x**6 + 10*x**4 + 10*x**3 + 10*x**2
    g = (1, [(x, 2), (x**2 + 9*x + 5, 1),
             (x**5 + 3*x**4 + 8*x**2 + 5*x + 2, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**32 + 1
    g = (1, [(x**16 + 3*x**8 + 10, 1),
             (x**16 + 8*x**8 + 10, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = 8*x**32 + 5
    g = (8, [(x + 3, 1), (x + 8, 1), (x**2 + 9, 1), (x**2 + 2*x + 2, 1),
             (x**2 + 9*x + 2, 1), (x**8 + x**4 + 6, 1), (x**8 + 10*x**4 + 6, 1),
             (x**4 + 5*x**2 + 7, 1), (x**4 + 6*x**2 + 7, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = 8*x**63 + 5
    g = (8, [(x + 7, 1), (x**6 + 9*x**3 + 4, 1), (x**2 + 4*x + 5, 1),
             (x**3 + 6*x**2 + 8*x + 2, 1), (x**3 + 9*x**2 + 9*x + 2, 1),
             (x**6 + 2*x**5 + 6*x**4 + 8*x**2 + 4*x + 4, 1),
             (x**6 + 2*x**5 + 8*x**3 + 4*x**2 + 6*x + 4, 1),
             (x**6 + 5*x**5 + 6*x**4 + 8*x**2 + 6*x + 4, 1),
             (x**6 + 2*x**5 + 3*x**4 + 8*x**3 + 6*x + 4, 1),
             (x**6 + 10*x**5 + 4*x**4 + 7*x**3 + 10*x**2 + 7*x + 4, 1),
             (x**6 + 3*x**5 + 3*x**4 + x**3 + 6*x**2 + 8*x + 4, 1),
             (x**6 + 6*x**5 + 2*x**4 + 7*x**3 + 9*x**2 + 8*x + 4, 1),
             (x**6 + 10*x**5 + 10*x**4 + x**3 + 4*x**2 + 9*x + 4, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**15 - 1
    g = (1, [(x + 2, 1), (x + 6, 1), (x + 7, 1), (x + 8, 1), (x + 10, 1),
             (x**2 + x + 1, 1), (x**2 + 5*x + 3, 1), (x**2 + 9*x + 4, 1),
             (x**2 + 4*x + 5, 1), (x**2 + 3*x + 9, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    with using(gf_factor_method='other'):
        pytest.raises(KeyError, lambda: (x + 1).factor_list())

    R, x = ring('x', FF(13))

    f = x**8 + x**6 + 10*x**4 + 10*x**3 + 8*x**2 + 2*x + 8
    g = (1, [(x + 3, 1), (x**3 + 8*x**2 + 4*x + 12, 1),
             (x**4 + 2*x**3 + 3*x**2 + 4*x + 6, 1)])

    with using(gf_factor_method='berlekamp'):
        assert f.factor_list() == g

    R, x = ring('x', FF(809))

    f = (x**10 + 2*x**9 + 5*x**8 + 26*x**7 + 677*x**6 + 436*x**5 +
         791*x**4 + 325*x**3 + 456*x**2 + 24*x + 577)
    g = (1, [(x + 701, 1), (x**9 + 110*x**8 + 559*x**7 + 532*x**6 +
                            694*x**5 + 151*x**4 + 110*x**3 + 70*x**2 +
                            735*x + 122, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    # Gathen polynomials: x**n + x + 1 (mod p > 2**n * pi)

    R, x = ring('x', FF(nextprime(2**15*pi)))

    f = x**15 + x + 1
    g = (1, [(x**2 + 22730*x + 68144, 1),
             (x**4 + 81553*x**3 + 77449*x**2 + 86810*x + 4724, 1),
             (x**4 + 86276*x**3 + 56779*x**2 + 14859*x + 31575, 1),
             (x**5 + 15347*x**4 + 95022*x**3 + 84569*x**2 + 94508*x + 92335, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    # Shoup polynomials: f = a_0 x**n + a_1 x**(n-1) + ... + a_n
    # (mod p > 2**(n-2) * pi), where a_n = a_{n-1}**2 + 1, a_0 = 1

    R, x = ring('x', FF(nextprime(2**4*pi)))

    f = x**6 + 2*x**5 + 5*x**4 + 26*x**3 + 41*x**2 + 39*x + 38
    g = (1, [(x**2 + 44*x + 26, 1),
             (x**4 + 11*x**3 + 25*x**2 + 18*x + 30, 1)])

    for method in ('zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    F8 = FF(2, [1, 1, 0, 1])
    R, x = ring('x', F8)

    f = x**10 + x**9 + F8(2)*x**8 + F8(2)*x**7 + F8(5)*x**6 + F8(3)*x**5
    g = (F8(1), [(x, 5), (x + F8(3), 1), (x + F8(6), 1),
                 (x**3 + F8(4)*x**2 + x + F8(3), 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    F9 = FF(3, [2, 2, 1])
    R, x = ring('x', F9)

    f = x**5 + F9(2)*x**4 + F9(6)*x**3 + F9(8)*x**2 + F9(5)*x + F9(4)
    g = (1, [(x + F9(8), 1), (x**2 + 2*x + F9(4), 1),
             (x**2 + F9(4)*x + F9(4), 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g


def test_PolyElement_is_irreducible():
    R, x = ring('x', FF(5))

    f = (x**10 + 4*x**9 + 2*x**8 + 2*x**7 + 3*x**6 +
         2*x**5 + 4*x**4 + x**3 + 4*x**2 + 4)
    g = 3*x**2 + 2*x + 4

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert f.is_irreducible is True
            assert g.is_irreducible is False

    R, x = ring('x', FF(11))

    f = R(7)
    g = 7*x + 3
    h = 7*x**2 + 3*x + 1

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert f.is_irreducible is True
            assert g.is_irreducible is True
            assert h.is_irreducible is False

    with using(gf_irred_method='other'):
        pytest.raises(KeyError, lambda: f.is_irreducible)

    R, x = ring('x', FF(13))

    f = 2*x**4 + 3*x**3 + 4*x**2 + 5*x + 6
    g = 2*x**4 + 3*x**3 + 4*x**2 + 5*x + 8

    with using(gf_irred_method='ben-or'):
        assert f.is_irreducible is False
        assert g.is_irreducible is True

    R, x = ring('x', FF(17))

    f = (x**10 + 9*x**9 + 9*x**8 + 13*x**7 + 16*x**6 + 15*x**5 +
         6*x**4 + 7*x**3 + 7*x**2 + 7*x + 10)
    g = (x**10 + 7*x**9 + 16*x**8 + 7*x**7 + 15*x**6 + 13*x**5 + 13*x**4 +
         11*x**3 + 16*x**2 + 10*x + 9)
    h = f*g

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert f.is_irreducible is True
            assert g.is_irreducible is True
            assert h.is_irreducible is False

    F9 = FF(3, [2, 2, 1])
    R, x = ring('x', F9)

    f = x**3 + F9(8)*x**2 + F9(8)*x + F9(4)

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert f.is_irreducible is False

    F27 = FF(3, [1, 0, 2, 1])
    R, x = ring('x', F27)

    f = x**3 + F27(8)*x**2 + F27(19)*x + F27(24)

    for method in ('ben-or', 'rabin'):
        with using(gf_irred_method=method):
            assert f.is_irreducible is True

    R, x = ring('x', ZZ)

    assert x.is_irreducible is True
    assert (x**2 + x + 1).is_irreducible is True
    assert (x**2 + 2*x + 1).is_irreducible is False
    assert (x**2 - 1).is_irreducible is False

    f = 3*x**4 + 2*x**3 + 6*x**2 + 8*x

    assert (f + 7).is_irreducible is True
    assert (f + 4).is_irreducible is True
    assert (f + 10).is_irreducible is True
    assert (f + 14).is_irreducible is True

    R, x, y = ring('x y', ZZ)

    assert R(2).is_irreducible is True
    assert (x**2 + x + 1).is_irreducible is True
    assert (x**2 + 2*x + 1).is_irreducible is False
    assert ((x - 2*y)*(x + y)).is_irreducible is False
    assert (x**2 + y**2).is_irreducible is True

    R, x, y, z = ring('x y z', QQ)

    assert (x**2 + x + 1).is_irreducible
    assert (x**2 + 2*x + 1).is_irreducible is False


@pytest.mark.timeout(50)
def test_sympyissue_16620():
    R, x = ring('x', FF(2))

    f = x**17 + 1
    g = (1, [(x + 1, 1),
             (x**8 + x**5 + x**4 + x**3 + 1, 1),
             (x**8 + x**7 + x**6 + x**4 + x**2 + x + 1, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g

    f = x**31 + 1
    g = (1, [(x + 1, 1), (x**5 + x**2 + 1, 1), (x**5 + x**3 + 1, 1),
             (x**5 + x**3 + x**2 + x + 1, 1), (x**5 + x**4 + x**2 + x + 1, 1),
             (x**5 + x**4 + x**3 + x + 1, 1), (x**5 + x**4 + x**3 + x**2 + 1, 1)])

    for method in ('berlekamp', 'zassenhaus', 'shoup'):
        with using(gf_factor_method=method):
            assert f.factor_list() == g


def test__gf_trace_map():
    R, x = ring('x', FF(5))

    a = x + 2
    b = 4*x + 4
    c = x + 1
    f = 3*x**2 + 2*x + 4

    assert R._gf_trace_map(a, b, c, 4, f) == (x + 3, x + 3)

    R, x = ring('x', FF(11))

    f = x**4 + x**3 + 4*x**2 + 9*x + 1
    a = x**2 + x + 1
    c = x
    b = pow(c, 11, f)

    assert R._gf_trace_map(a, b, c, 0, f) == (x**2 + x + 1, x**2 + x + 1)
    assert R._gf_trace_map(a, b, c, 1, f) == (5*x**3 + 2*x**2 + 10*x + 3,
                                              5*x**3 + 3*x**2 + 4)
    assert R._gf_trace_map(a, b, c, 2, f) == (5*x**3 + 9*x**2 + 5*x + 3,
                                              10*x**3 + x**2 + 5*x + 7)
    assert R._gf_trace_map(a, b, c, 3, f) == (x**3 + 10*x**2 + 6*x, 7)
    assert R._gf_trace_map(a, b, c, 4, f) == (x**2 + x + 1, x**2 + x + 8)
    assert R._gf_trace_map(a, b, c, 5, f) == (5*x**3 + 2*x**2 + 10*x + 3,
                                              5*x**3 + 3*x**2)
    assert R._gf_trace_map(a, b, c, 11, f) == (x**3 + 10*x**2 + 6*x, 10)
