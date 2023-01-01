"""Tests for Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences."""

import pytest

from diofant import (CC, FF, QQ, RR, ZZ, I, NotInvertibleError, field, ring,
                     sqrt)
from diofant.config import using
from diofant.polys.specialpolys import f_polys


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = f_polys()


def test_gcdex():
    R, x = ring('x', FF(11))

    assert R.zero.gcdex(R(2)) == (0, 6, 1)
    assert R(2).gcdex(R(2)) == (0, 6, 1)

    assert R.zero.gcdex(3*x) == (0, 4, x)

    assert (3*x).gcdex(3*x) == (0, 4, x)

    assert (x**2 + 8*x + 7).gcdex(x**3 + 7*x**2 + x + 7) == (5*x + 6, 6, x + 7)

    R, x = ring('x', QQ)

    f = x**4 - 2*x**3 - 6*x**2 + 12*x + 15
    g = x**3 + x**2 - 4*x - 4

    s = -x/5 + QQ(3, 5)
    t = x**2/5 - 6*x/5 + 2
    h = x + 1

    assert f.half_gcdex(g) == (s, h)
    assert f.gcdex(g) == (s, t, h)

    f = x**4 + 4*x**3 - x + 1
    g = x**3 - x + 1

    s, t, h = f.gcdex(g)
    S, T, H = g.gcdex(f)

    assert s*f + t*g == h
    assert S*g + T*f == H

    f = 2*x
    g = x**2 - 16

    s = x/32
    t = -QQ(1, 16)
    h = 1

    assert f.half_gcdex(g) == (s, h)
    assert f.gcdex(g) == (s, t, h)


def test_dup_invert():
    R, x = ring('x', QQ)

    assert R.invert(2*x, x**2 - 16) == x/32
    assert R.invert(x**2 - 1, 2*x - 1) == QQ(-4, 3)

    pytest.raises(NotInvertibleError, lambda: R.invert(x**2 - 1, x - 1))


def test_dmp_prem():
    _, x = ring('x', FF(7))

    f = x**2 + x + 3
    g = 2*x + 2

    assert f.prem(g) == 5  # issue sympy/sympy#20397

    _, x = ring('x', ZZ)

    f = 3*x**3 + x**2 + x + 5
    g = 5*x**2 - 3*x + 1

    r = 52*x + 111

    assert f.prem(g) == r

    pytest.raises(ZeroDivisionError, lambda: f.prem(0))

    f = x**2 + 1
    g = 2*x - 4
    r = 20

    assert f.prem(g) == r

    _, x = ring('x', QQ)

    f = 3*x**3 + x**2 + x + 5
    g = 5*x**2 - 3*x + 1

    r = 52*x + 111

    assert g.prem(f) == g
    assert f.prem(g) == r

    _, x, y = ring('x y', ZZ)

    f = x**2 - y**2
    g = x - y

    assert f.prem(g) == 0

    f = x**2 + y**2
    g = x - y

    r = 2*y**2

    assert f.prem(g) == r

    pytest.raises(ZeroDivisionError, lambda: f.prem(0))

    g = 2*x - 2*y

    r = 8*y**2

    assert g.prem(f) == g
    assert f.prem(g) == r

    f = x**2 + x*y
    g = 2*x + 2

    r = -4*y + 4

    assert f.prem(g) == r


def test_PolyElement_subresultants():
    R, x = ring('x', ZZ)

    for check in (True, False):
        with using(use_collins_resultant=check):
            assert R(0).resultant(R(0)) == 0
            assert R(0).resultant(R(0), includePRS=True) == (0, [])
            assert R(1).resultant(R(0)) == 0
            assert R(1).subresultants(R(0)) == [1]
            assert R(0).resultant(R(1)) == 0
            assert R(0).resultant(R(1), includePRS=True) == (0, [1])

            f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
            g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

            a = 15*x**4 - 3*x**2 + 9
            b = 65*x**2 + 125*x - 245
            c = 9326*x - 12300
            d = R(260708)

            assert f.subresultants(g) == [f, g, a, b, c, d]
            assert f.resultant(g) == d.drop(x)

            f = x**2 - 2*x + 1
            g = x**2 - 1

            a = 2*x - 2

            assert f.subresultants(g) == [f, g, a]
            assert f.resultant(g) == 0

            f = x**2 + 1
            g = x**2 - 1

            a = -2

            assert f.subresultants(g) == [f, g, a]
            assert f.resultant(g) == 4
            assert f.resultant(g, includePRS=True) == (4, [x**2 + 1, x**2 - 1, -2])

            f = x**2 - 1
            g = x**3 - x**2 + 2

            assert f.resultant(g) == 0

            f = 3*x**3 - x
            g = 5*x**2 + 1

            assert f.resultant(g) == 64

            f = x**2 - 2*x + 7
            g = x**3 - x + 5

            assert f.resultant(g) == 265

            f = x**3 - 6*x**2 + 11*x - 6
            g = x**3 - 15*x**2 + 74*x - 120

            assert f.resultant(g) == -8640

            f = x**3 - 6*x**2 + 11*x - 6
            g = x**3 - 10*x**2 + 29*x - 20

            assert f.resultant(g) == 0

            f = x**3 - 1
            g = x**3 + 2*x**2 + 2*x - 1

            assert f.resultant(g) == 16

            f = x**8 - 2
            g = x - 1

            assert f.resultant(g) == -1

            # issue sympy/sympy#10666
            f = x**3 - 7*x + 7
            g = x

            assert f.resultant(g) == -g.resultant(f) == -7

    Rt, t = ring('t', ZZ)
    _, x = ring('x', Rt)

    f = x**6 - 5*x**4 + 5*x**2 + 4
    g = -6*t*x**5 + x**4 + 20*t*x**3 - 3*x**2 - 10*t*x + 6

    assert f.resultant(g) == 2930944*t**6 + 2198208*t**4 + 549552*t**2 + 45796
    assert (x - 1).resultant(x + 1, includePRS=True) == (2, [x - 1, x + 1, 2])

    R, x, y = ring('x y', ZZ)

    for check in (True, False):
        with using(use_collins_resultant=check):
            assert R(0).resultant(R(0)) == 0
            assert R(0).resultant(R(0), includePRS=True) == (0, [])

            assert R(0).resultant(R(1)) == 0
            assert R(1).resultant(R(0)) == 0
            assert R(1).subresultants(R(0)) == [1]
            assert R(0).resultant(R(1), includePRS=True) == (0, [1])

            f = x + y + 2
            g = 2*x*y + x + 3

            assert f.resultant(g) == (-2*y**2 - 5*y + 1).drop(x)

            f = 3*x**2*y - y**3 - 4
            g = x**2 + x*y**3 - 9

            a = 3*x*y**4 + y**3 - 27*y + 4
            b = (-3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 +
                 729*y**2 - 216*y + 16)

            r = b.drop(x)
            rr = (r, [3*x**2*y - y**3 - 4, x**2 + x*y**3 - 9, 3*x*y**4 +
                      y**3 - 27*y + 4, -3*y**10 - 12*y**7 + y**6 -
                      54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16])

            assert f.subresultants(g) == [f, g, a, b]

            assert f.resultant(g) == r
            assert f.resultant(g, includePRS=True) == rr

            f = -x**3 + 5
            g = 3*x**2*y + x**2

            a = 45*y**2 + 30*y + 5
            b = 675*y**3 + 675*y**2 + 225*y + 25

            r = b.drop(x)

            assert f.subresultants(g) == [f, g, a]
            assert f.resultant(g) == r
            assert f.resultant(g, includePRS=True)[0] == r

            f = x + y
            g = x**2 - x*y + 1

            assert f.resultant(g) == (1 + 2*y**2).drop(x)

            g += 1

            assert f.resultant(g) == (2 + 2*y**2).drop(x)

    R, x, y = ring('x y', QQ)

    for check in (True, False):
        with using(use_collins_resultant=check):
            assert R(0).resultant(R(0)) == 0
            assert R(0).resultant(R(1)) == 0

            f = x + y
            g = x**2 - x*y + 1

            assert f.resultant(g) == (1 + 2*y**2).drop(x)

            f = x/2 + y + QQ(2, 3)
            g = 2*x*y + x + 3

            assert f.resultant(g) == (-2*y**2 - 7*y/3 + QQ(5, 6)).drop(x)

            f = 3*x**2*y - y**3 - 4
            g = x**2 + x*y**3 - 9

            assert f.resultant(g) == (-3*y**10 - 12*y**7 + y**6 - 54*y**4 +
                                      8*y**3 + 729*y**2 - 216*y + 16).drop(x)

            f = -x**3 + 5
            g = 3*x**2*y + x**2

            assert f.resultant(g) == (675*y**3 + 675*y**2 + 225*y + 25).drop(x)

    R, x, y, z, u, v = ring('x y z u v', ZZ)

    for check in (True, False):
        with using(use_collins_resultant=check):
            f = 6*x**2 - 3*x*y - 2*x*z + y*z
            g = x**2 - x*u - x*v + u*v

            r = (y**2*z**2 - 3*y**2*z*u - 3*y**2*z*v + 9*y**2*u*v -
                 2*y*z**2*u - 2*y*z**2*v + 6*y*z*u**2 + 12*y*z*u*v +
                 6*y*z*v**2 - 18*y*u**2*v - 18*y*u*v**2 + 4*z**2*u*v -
                 12*z*u**2*v - 12*z*u*v**2 + 36*u**2*v**2)

            assert f.resultant(g) == r.drop(x)

    R, x, y, z, u, v = ring('x y z u v', QQ)

    for check in (True, False):
        with using(use_collins_resultant=check):
            f = x**2 - x*y/2 - x*z/3 + y*z/6
            g = x**2 - x*u - x*v + u*v

            r = (y**2*z**2/36 - y**2*z*u/12 - y**2*z*v/12 + y**2*u*v/4 -
                 y*z**2*u/18 - y*z**2*v/18 + y*z*u**2/6 + y*z*u*v/3 +
                 y*z*v**2/6 - y*u**2*v/2 - y*u*v**2/2 + z**2*u*v/9 -
                 z*u**2*v/3 - z*u*v**2/3 + u**2*v**2)

            assert f.resultant(g) == r.drop(x)


def test_PolyElement_discriminant():
    R, x = ring('x', ZZ)

    assert R(0).discriminant() == 0
    assert x.discriminant() == 1

    assert (x**3 + 3*x**2 + 9*x - 13).discriminant() == -11664
    assert (5*x**5 + x**3 + 2).discriminant() == 31252160
    assert (x**4 + 2*x**3 + 6*x**2 - 22*x + 13).discriminant() == 0
    assert (12*x**7 + 15*x**4 + 30*x**3 + x**2 + 1).discriminant() == -220289699947514112

    assert (x**2 + 2*x + 3).discriminant() == -8

    R, x, y = ring('x y', ZZ)

    assert R(0).discriminant() == 0
    assert y.discriminant() == 0

    assert (x**3 + 3*x**2 + 9*x - 13).discriminant() == -11664
    assert (5*x**5 + x**3 + 2).discriminant() == 31252160
    assert (x**4 + 2*x**3 + 6*x**2 - 22*x + 13).discriminant() == 0
    assert (12*x**7 + 15*x**4 + 30*x**3 + x**2 + 1).discriminant() == -220289699947514112

    assert (x**2*y + 2*y).discriminant() == (-8*y**2).drop(x)
    assert (x*y**2 + 2*x).discriminant() == 1

    R, x, y, z = ring('x y z', ZZ)

    assert (x*y + z).discriminant() == 1

    R, x, y, z, u = ring('x y z u', ZZ)

    assert (x**2*y + x*z + u).discriminant() == (-4*y*u + z**2).drop(x)

    R, x, y, z, u, v = ring('x y z u v', ZZ)

    assert (x**3*y + x**2*z + x*u + v).discriminant() == \
        (-27*y**2*v**2 + 18*y*z*u*v - 4*y*u**3 - 4*z**3*v + z**2*u**2).drop(x)

    F, a, b, c = ring('a b c', ZZ)
    _, x = ring('x', F)

    f, g = a*x**2 + b*x + c, b**2 - 4*a*c

    assert f.discriminant() == g


def test_dmp_gcd():
    R, x = ring('x', FF(5))

    f = 3*x**2 + 2*x + 4
    g = 2*x**2 + 2*x + 3

    assert f.cofactors(g) == (x + 3, 3*x + 3, 2*x + 1)

    R, x = ring('x', FF(11))

    assert R(0).cofactors(R(0)) == (0, 0, 0)
    assert R(2).cofactors(R(0)) == (1, 2, 0)
    assert R(0).cofactors(R(2)) == (1, 0, 2)
    assert R(2).cofactors(R(2)) == (1, 2, 2)

    assert R(0).cofactors(x) == (x, 0, 1)
    assert x.cofactors(R(0)) == (x, 1, 0)

    assert (3*x).cofactors(3*x) == (x, 3, 3)
    assert (x**2 + 8*x + 7).cofactors(x**3 + 7*x**2 + x + 7) == (x + 7, x + 1,
                                                                 x**2 + 1)

    R, x = ring('x', ZZ)

    for test in (True, False):
        for method in ('prs', 'modgcd'):
            with using(use_heu_gcd=test, fallback_gcd_zz_method=method):
                assert R(0).cofactors(R(0)) == (0, 0, 0)
                assert R(0).cofactors(x) == (x, 0, 1)
                assert x.cofactors(R(0)) == (x, 1, 0)
                assert R(0).cofactors(-x) == (x, 0, -1)
                assert (-x).cofactors(R(0)) == (x, -1, 0)
                assert (2*x).cofactors(R(2)) == (2, x, 1)
                assert R(2).cofactors(R(0)) == (2, 1, 0)
                assert R(-2).cofactors(R(0)) == (2, -1, 0)
                assert R(0).cofactors(R(-2)) == (2, 0, -1)
                assert R(0).cofactors(2*x + 4) == (2*x + 4, 0, 1)
                assert (2*x + 4).cofactors(R(0)) == (2*x + 4, 1, 0)
                assert R(2).cofactors(R(2)) == (2, 1, 1)
                assert R(-2).cofactors(R(2)) == (2, -1, 1)
                assert R(2).cofactors(R(-2)) == (2, 1, -1)
                assert R(-2).cofactors(R(-2)) == (2, -1, -1)
                assert (x**2 + 2*x + 1).cofactors(R(1)) == (1, x**2 + 2*x + 1, 1)
                assert (x**2 + 2*x + 1).cofactors(R(2)) == (1, x**2 + 2*x + 1, 2)
                assert (2*x**2 + 4*x + 2).cofactors(R(2)) == (2, x**2 + 2*x + 1, 1)
                assert R(2).cofactors(2*x**2 + 4*x + 2) == (2, 1, x**2 + 2*x + 1)
                assert (2*x**2 + 4*x + 2).cofactors(x + 1) == (x + 1, 2*x + 2, 1)
                assert (x + 1).cofactors(2*x**2 + 4*x + 2) == (x + 1, 1, 2*x + 2)
                assert (x - 31).cofactors(x) == (1, x - 31, x)

                f, g = 2*x + 2, 6*x**2 - 6

                assert f.cofactors(g) == (2*x + 2, 1, 3*x - 3)

                f, g = [1000000000000*x + 998549000000]*2

                assert f.cofactors(g) == (f, 1, 1)

                f, g = 999530000000*x + 1000000000000, 999530000000*x + 999999000000

                assert f.cofactors(g) == (1000000, 999530*x + 1000000, 999530*x + 999999)

                f = x**2 - 1
                g = x**2 - 3*x + 2

                assert f.cofactors(g) == (x - 1, x + 1, x - 2)

                f = x**4 + 8*x**3 + 21*x**2 + 22*x + 8
                g = x**3 + 6*x**2 + 11*x + 6

                assert f.cofactors(g) == (x**2 + 3*x + 2, x**2 + 5*x + 4, x + 3)

                f = x**4 - 4
                g = x**4 + 4*x**2 + 4

                assert f.cofactors(g) == (x**2 + 2, x**2 - 2, x**2 + 2)

                f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
                g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

                assert f.cofactors(g) == (1, f, g)

                f = (-352518131239247345597970242177235495263669787845475025293906825864749649589178600387510272*x**49 +
                     46818041807522713962450042363465092040687472354933295397472942006618953623327997952*x**42 +
                     378182690892293941192071663536490788434899030680411695933646320291525827756032*x**35 +
                     112806468807371824947796775491032386836656074179286744191026149539708928*x**28 -
                     12278371209708240950316872681744825481125965781519138077173235712*x**21 +
                     289127344604779611146960547954288113529690984687482920704*x**14 +
                     19007977035740498977629742919480623972236450681*x**7 +
                     311973482284542371301330321821976049)

                h = (365431878023781158602430064717380211405897160759702125019136*x**21 +
                     197599133478719444145775798221171663643171734081650688*x**14 -
                     9504116979659010018253915765478924103928886144*x**7 -
                     311973482284542371301330321821976049)
                cff = (-964661685087874498642420170752*x**28 + 649736296036977287118848*x**21 +
                       658473216967637120*x**14 - 30463679113*x**7 - 1)
                cfg = (-47268422569305850433478588366848*x**27 + 30940259392972115602096128*x**20 +
                       18261628279718027904*x**13 - 426497272383*x**6)

                assert f.cofactors(f.diff()) == (h, cff, cfg)

                f = 1317378933230047068160*x + 2945748836994210856960
                g = 120352542776360960*x + 269116466014453760

                H, cff, cfg = 120352542776360960*x + 269116466014453760, 10946, 1

                assert f.cofactors(g) == (H, cff, cfg)

                with using(heu_gcd_max=0):
                    assert f.cofactors(g) == (H, cff, cfg)

    R, x = ring('x', QQ)

    for test in (True, False):
        with using(use_heu_gcd=test, fallback_gcd_zz_method='prs'):
            f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
            g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

            assert f.cofactors(g) == (1, f, g)

            assert R(0).cofactors(R(0)) == (0, 0, 0)

            f, g = x**2/2 + x + QQ(1, 2), x/2 + QQ(1, 2)

            assert f.cofactors(g) == (x + 1, g, QQ(1, 2))

            f, g = x**2 - 1, x**2 - 3*x + 2

            assert f.cofactors(g) == (x - 1, x + 1, x - 2)

    R, x = ring('x', QQ.algebraic_field(sqrt(2)))

    for method in ('modgcd', 'prs'):
        with using(gcd_aa_method=method):
            f, g = 2*x, R(2)

            assert f.cofactors(g) == (1, f, g)

            f, g = 2*x, R(sqrt(2))

            assert f.cofactors(g) == (1, f, g)

            f, g = 2*x + 2, 6*x**2 - 6

            assert f.cofactors(g) == (x + 1, 2, 6*x - 6)

    R, x = ring('x', QQ.algebraic_field(sqrt(2)**(-1)*sqrt(3)))

    for method in ('modgcd', 'prs'):
        with using(gcd_aa_method=method):
            f, g = x + 1, x - 1

            assert f.cofactors(g) == (1, f, g)

    R, x = ring('x', CC)

    f, g = x**2 - 1, x**3 - 3*x + 2

    assert f.cofactors(g) == (x - 1, x + 1, x**2 + x - 2)

    R, x, y = ring('x y', ZZ)

    for test in (True, False):
        for method in ('prs', 'modgcd'):
            with using(use_heu_gcd=test, fallback_gcd_zz_method=method):
                assert R(0).cofactors(R(0)) == (0, 0, 0)
                assert R(2).cofactors(R(0)) == (2, 1, 0)
                assert R(-2).cofactors(R(0)) == (2, -1, 0)
                assert R(0).cofactors(R(-2)) == (2, 0, -1)
                assert R(0).cofactors(2*x + 4) == (2*x + 4, 0, 1)
                assert (2*x).cofactors(R(2)) == (2, x, 1)
                assert (2*x + 4).cofactors(R(0)) == (2*x + 4, 1, 0)
                assert R(2).cofactors(R(2)) == (2, 1, 1)
                assert R(-2).cofactors(R(2)) == (2, -1, 1)
                assert R(2).cofactors(R(-2)) == (2, 1, -1)
                assert R(-2).cofactors(R(-2)) == (2, -1, -1)
                assert (x**2 + 2*x + 1).cofactors(R(1)) == (1, x**2 + 2*x + 1, 1)
                assert (x**2 + 2*x + 1).cofactors(R(2)) == (1, x**2 + 2*x + 1, 2)
                assert (2*x**2 + 4*x + 2).cofactors(R(2)) == (2, x**2 + 2*x + 1, 1)
                assert R(2).cofactors(2*x**2 + 4*x + 2) == (2, 1, x**2 + 2*x + 1)

                f, g = 2*x**2 + 4*x + 2, x + 1

                assert f.cofactors(g) == (g, 2*x + 2, 1)
                assert g.cofactors(f) == (g, 1, 2*x + 2)

                with using(heu_gcd_max=0):
                    assert f.cofactors(g) == (g, 2*x + 2, 1)

                f = x**4 + 8*x**3 + 21*x**2 + 22*x + 8
                g = x**3 + 6*x**2 + 11*x + 6

                assert f.cofactors(g) == (x**2 + 3*x + 2, x**2 + 5*x + 4, x + 3)

                f, g = x + 2*y, x + y

                assert f.cofactors(g) == (1, f, g)

                f, g = x**2 + 2*x*y + y**2, x**2 + x*y

                assert f.cofactors(g) == (x + y, x + y, x)

                f, g = x**2 + 2*x*y + y**2, x**3 + y**3

                assert f.cofactors(g) == (x + y, x + y, x**2 - x*y + y**2)

                f, g = x*y**2 + 2*x*y + x, x*y**3 + x

                assert f.cofactors(g) == (x*y + x, y + 1, y**2 - y + 1)

                f, g = x**2*y**2 + x**2*y + 1, x*y**2 + x*y + 1

                assert f.cofactors(g) == (1, f, g)

                f = 2*x*y**2 + 4*x*y + 2*x + y**2 + 2*y + 1
                g = 2*x*y**3 + 2*x + y**3 + 1

                assert f.cofactors(g) == (2*x*y + 2*x + y + 1, y + 1, y**2 - y + 1)

                f = 2*x**2 + 4*x*y - 2*x - 4*y
                g = x**2 + x - 2

                assert f.cofactors(g) == (x - 1, 2*x + 4*y, x + 2)

                f = 2*x**2 + 2*x*y - 3*x - 3*y
                g = 4*x*y - 2*x + 4*y**2 - 2*y

                assert f.cofactors(g) == (x + y, 2*x - 3, 4*y - 2)

                f = (-17434367009167300000000000000000000000000000000000000000000000000000000*x**4*y -
                     250501827896299135568887342575961783764139560000000000000000000000000000000000000000000*x**3*y -
                     2440935909299672540738135183426056447877858000000000000000000000000000000*x**3 -
                     1349729941723537919695626818065131519270095220127010623905326719279566297660000000000000000000000000000*x**2*y -
                     26304033868956978374552886858060487282904504027042515077682955951658838800000000000000000*x**2 -
                     3232215785736369696036755035364398565076440134133908303058376297547504030528179314849416971379040931276000000000000000*x*y -
                     94485916261760032526508027937078714464844205539023800247528621905831259414691631156161537919255129011800*x -
                     2902585888465621357542575571971656665554321652262249362701116665830760628936600958940851960635161420991047110815678789984677193092993*y -
                     113133324167442997472440652189550843502029192913459268196939183295294085146407870078840385860571627108778756267503630290)
                g = (10000000000000000000000000000*x**2 + 71841388839807267676152024786000000000000000*x +
                     129029628760809605749020969023932901278290735413660734705971)

                h = (-1743436700916730000000000000000000000000000*x**2*y -
                     12525091394814956778444367128798089188206978000000000000000*x*y -
                     244093590929967254073813518342605644787785800*x -
                     22495499028725631994927113634418779135935898997901327211111875586270479483*y -
                     876801128965234839118530545935732755107147297241756982389990)

                assert f.cofactors(g) == (g, h, 1)

    R, x, y = ring('x y', QQ)

    for test in (True, False):
        with using(use_heu_gcd=test, fallback_gcd_zz_method='prs'):
            f, g = x**2/2 + x + QQ(1, 2), x/2 + QQ(1, 2)

            assert f.cofactors(g) == (x + 1, g, QQ(1, 2))
            assert g.cofactors(f) == (x + 1, QQ(1, 2), g)

            assert f.gcd(g) == x + 1
            with using(fallback_gcd_zz_method='modgcd'):
                assert f.gcd(g) == x + 1

            assert R(0).cofactors(R(0)) == (0, 0, 0)
            assert R(0).cofactors(g) == (x + 1, 0, QQ(1, 2))

            f, g = x**2/4 + x*y + y**2, x**2/2 + x*y

            assert f.cofactors(g) == (x + 2*y, x/4 + y/2, x/2)

            f, g = x**2/2 + x*y + y**2/2, x**2 + x*y

            assert f.cofactors(g) == (x + y, x/2 + y/2, x)

    R, x, y = ring('x y', QQ.algebraic_field(sqrt(2)))

    for method in ('modgcd', 'prs'):
        with using(gcd_aa_method=method):
            f, g = (x + sqrt(2)*y)**2, x + sqrt(2)*y

            assert f.cofactors(g) == (g, g, 1)

            f, g = x + sqrt(2)*y, x + y

            assert f.cofactors(g) == (1, f, g)

            f, g = x*y + sqrt(2)*y**2, sqrt(2)*y

            assert f.cofactors(g) == (y, x + sqrt(2)*y, sqrt(2))

            f, g = x**2 + 2*sqrt(2)*x*y + 2*y**2, x + sqrt(2)*y

            assert f.cofactors(g) == (g, g, 1)

    R, x, y = ring('x y', RR)

    for test in (True, False):
        with using(use_heu_gcd=test, fallback_gcd_zz_method='prs'):
            f, g = 2.1*x*y**2 - 2.2*x*y + 2.1*x, 1.0*x**3
            h = 1.0*x

            assert f.cofactors(g) == (h, 2.1*y**2 - 2.2*y + 2.1, 1.0*x**2)

            f, g = 2.1*x*y**2 - 2.1*x*y + 2.1*x, 2.1*x**3

            assert f.cofactors(g) == (h, f//h, g//h)
            assert g.cofactors(f) == (h, g//h, f//h)

    R, x, y = ring('x y', CC)

    f, g = x**2 - y, x**3 - y*x + 2

    assert f.cofactors(g) == (1, f, g)

    R, x, y, z = ring('x y z', ZZ)

    for test in (True, False):
        for method in ('prs', 'modgcd'):
            with using(use_heu_gcd=test, fallback_gcd_zz_method=method):
                f, g = x - y*z, x - y*z

                assert f.cofactors(g) == (x - y*z, 1, 1)

                f, g, h = R.fateman_poly_F_1()
                H, cff, cfg = f.cofactors(g)

                assert H == h
                assert H*cff == f
                assert H*cfg == g

                f, g, h = R.fateman_poly_F_2()
                H, cff, cfg = f.cofactors(g)

                assert H == h
                assert H*cff == f
                assert H*cfg == g

                f, g, h = R.fateman_poly_F_3()
                H, cff, cfg = f.cofactors(g)

                assert H == h
                assert H*cff == f
                assert H*cfg == g

    R, x, y, z = ring('x y z', QQ.algebraic_field(sqrt(2), sqrt(3)))

    with using(gcd_aa_method='modgcd'):
        h = x**2*y**7 + sqrt(6)/21*z
        f, g = h*(27*y**3 + 1), h*(y + x)

        assert f.cofactors(g) == (h, 27*y**3 + 1, x + y)

        h = x**13*y**3 + x**10/2 + 1/sqrt(2)
        f, g = h*(x + 1), h*sqrt(2)/sqrt(3)

        assert f.cofactors(g) == (h, x + 1, sqrt(2)/sqrt(3))

        h = x**4*y**9 + sqrt(6)/22*z
        f, g = h*(21*y**3 + 1), h*(y + x)

        assert f.cofactors(g) == (x**4*y**9 + sqrt(6)/22*z, 21*y**3 + 1, x + y)

        h = x**4*y**3 + sqrt(6)/22*z
        f, g = h*(11*y**3 + 1), h*(y + x)

        assert f.cofactors(g) == (x**4*y**3 + sqrt(6)/22*z, 11*y**3 + 1, x + y)

        h = x**2*y**3 + 1111*sqrt(6)/12*z
        a, b = 11*y**3 + 2, (y + x - 1)*h

        assert (h*a).cofactors(h*b) == (h, a, b)

        a, b = 12*y + 2*x - 1, (y + x - 1)*h

        assert (h*a).cofactors(h*b) == (h, a, b)

    R, x, y, z = ring('x y z', QQ.algebraic_field(I))

    for method in ('prs', 'modgcd'):
        with using(gcd_aa_method=method):
            f, g = R.one, I*z

            assert f.cofactors(g) == (1, f, g)
            assert g.cofactors(f) == (1, g, f)

    R, x, y, z, u = ring('x y z u', ZZ)

    for test in (True, False):
        for method in ('prs', 'modgcd'):
            with using(use_heu_gcd=test, fallback_gcd_zz_method=method):
                f, g = u**2 + 2*u + 1, 2*u + 2

                assert f.cofactors(g) == (u + 1, u + 1, 2)

                f, g = z**2*u**2 + 2*z**2*u + z**2 + z*u + z, u**2 + 2*u + 1
                h, cff, cfg = u + 1, z**2*u + z**2 + z, u + 1

                assert f.cofactors(g) == (h, cff, cfg)
                assert g.cofactors(f) == (h, cfg, cff)

                f, g = x + y + z, -x - y - z - u

                assert f.cofactors(g) == (1, f, g)

                f, g, h = R.fateman_poly_F_3()
                H, cff, cfg = f.cofactors(g)

                assert H == h
                assert H*cff == f
                assert H*cfg == g

                f, g, h = (1199999999999991*x**17 - y, 2*y - 19989798798 + x**211,
                           12*x*y**7 + x**4 - 1)

                for _ in range(10):
                    assert (f*h).cofactors(g*h) == (h, f, g)

    R, x, y, z, u, _ = ring('x y z u v', ZZ)

    for test in (True, False):
        with using(use_heu_gcd=test, fallback_gcd_zz_method='modgcd'):
            f, g, h = R.fateman_poly_F_1()
            H, cff, cfg = f.cofactors(g)

            assert H == h
            assert H*cff == f
            assert H*cfg == g

            f, g, h = R.fateman_poly_F_3()
            H, cff, cfg = f.cofactors(g)

            assert H == h
            assert H*cff == f
            assert H*cfg == g

    R, x, y, z, u, _, a, b = ring('x y z u v a b', ZZ)

    for test in (True, False):
        with using(use_heu_gcd=test, fallback_gcd_zz_method='modgcd'):
            f, g, h = R.fateman_poly_F_1()
            H, cff, cfg = f.cofactors(g)

            assert H == h
            assert H*cff == f
            assert H*cfg == g

    R, x, y, z, u, _, a, b, *_ = ring('x y z u v a b c d', ZZ)

    for test in (True, False):
        with using(use_heu_gcd=test, fallback_gcd_zz_method='modgcd'):
            f, g, h = R.fateman_poly_F_1()
            H, cff, cfg = f.cofactors(g)

            assert H == h
            assert H*cff == f
            assert H*cfg == g

    F, x = field('x', QQ)
    R, _ = ring('t', F)

    assert R(x).gcd(R(0)) == 1


def test_PolyElement_lcm():
    R, x = ring('x', FF(5))

    assert (3*x**2 + 2*x + 4).lcm(2*x**2 + 2*x + 3) == x**3 + 2*x**2 + 4

    R, x = ring('x', FF(11))

    assert R.zero.lcm(R(2)) == 0
    assert R(2).lcm(R(2)) == 1

    assert R.zero.lcm(x) == 0

    assert (3*x).lcm(3*x) == x
    assert (x**2 + 8*x + 7).lcm(x**3 + 7*x**2 + x + 7) == (x**4 + 8*x**3 +
                                                           8*x**2 + 8*x + 7)

    R, x = ring('x', ZZ)

    assert R(2).lcm(R(6)) == 6

    assert (2*x**3).lcm(6*x) == 6*x**3
    assert (2*x**3).lcm(3*x) == 6*x**3

    assert (x**2 + x).lcm(x) == x**2 + x
    assert (x**2 + x).lcm(2*x) == 2*x**2 + 2*x
    assert (x**2 + 2*x).lcm(x) == x**2 + 2*x
    assert (2*x**2 + x).lcm(x) == 2*x**2 + x
    assert (2*x**2 + x).lcm(2*x) == 4*x**2 + 2*x
    assert (x**2 - 1).lcm(x**2 - 3*x + 2) == x**3 - 2*x**2 - x + 2

    R, x = ring('x', QQ)

    f = (x**2 + 7*x/2 + 3)/2
    g = x**2/2 + x
    h = x**3 + 7/2*x**2 + 3*x

    assert f.lcm(g) == h

    R, x, y = ring('x y', ZZ)

    assert R(2).lcm(R(6)) == 6
    assert x.lcm(y) == x*y

    assert (2*x**3).lcm(6*x*y**2) == 6*x**3*y**2
    assert (2*x**3).lcm(3*x*y**2) == 6*x**3*y**2

    assert (x**2*y).lcm(x*y**2) == x**2*y**2

    f = 2*x*y**5 - 3*x*y**4 - 2*x*y**3 + 3*x*y**2
    g = y**5 - 2*y**3 + y
    h = 2*x*y**7 - 3*x*y**6 - 4*x*y**5 + 6*x*y**4 + 2*x*y**3 - 3*x*y**2

    assert f.lcm(g) == h

    f = x**3 - 3*x**2*y - 9*x*y**2 - 5*y**3
    g = x**4 + 6*x**3*y + 12*x**2*y**2 + 10*x*y**3 + 3*y**4
    h = x**5 + x**4*y - 18*x**3*y**2 - 50*x**2*y**3 - 47*x*y**4 - 15*y**5

    assert f.lcm(g) == h

    f = x**2 + 2*x*y + y**2
    g = x**2 + x*y
    h = x**3 + 2*x**2*y + x*y**2

    assert f.lcm(g) == h

    R, x, y = ring('x y', QQ)

    f = 2*x*y - x**2/2 + QQ(1, 3)
    g = 3*x**3 - x*y**2 - QQ(1, 2)
    h = (x**5 - 4*x**4*y - x**3*y**2/3 - 2*x**3/3 + 4*x**2*y**3/3 -
         x**2/6 + 2*x*y**2/9 + 2*x*y/3 + QQ(1, 9))

    assert f.lcm(g) == h

    f = x**2/4 + x*y + y**2
    g = x**2/2 + x*y
    h = x**3 + 4*x**2*y + 4*x*y**2

    assert f.lcm(g) == h


def test_PolyElement_cancel():
    R, x = ring('x', ZZ)

    f = 2*x**2 - 2
    g = x**2 - 2*x + 1

    p = 2*x + 2
    q = x - 1

    assert f.cancel(g) == (p, q)
    assert f.cancel(g, include=False) == (1, 1, p, q)

    f = -x - 2
    g = 3*x - 4

    F = x + 2
    G = -3*x + 4

    assert f.cancel(g) == (f, g)
    assert F.cancel(G) == (f, g)

    assert R(0).cancel(R(0)) == (0, 0)
    assert R(0).cancel(R(0), include=False) == (1, 1, 0, 0)

    assert x.cancel(R(0)) == (1, 0)
    assert x.cancel(R(0), include=False) == (1, 1, 1, 0)

    assert R(0).cancel(x) == (0, 1)
    assert R(0).cancel(x, include=False) == (1, 1, 0, 1)

    f = R(0)
    g = x
    one = 1

    assert f.cancel(g, include=True) == (f, one)

    R, x = ring('x', QQ)

    assert (x**2/4 - 1).cancel(x/2 - 1) == (x + 2, 2)

    Fx, x = field('x', ZZ)
    _, t = ring('t', Fx)

    f = (-x**2 - 4)/4*t
    g = t**2 + (x**2 + 2)/2

    assert f.cancel(g) == ((-x**2 - 4)*t, 4*t**2 + 2*x**2 + 4)

    R, x, y = ring('x y', ZZ)

    f = 2*x**2 - 2
    g = x**2 - 2*x + 1

    p = 2*x + 2
    q = x - 1

    assert f.cancel(g) == (p, q)
    assert f.cancel(g, include=False) == (1, 1, p, q)

    assert R(0).cancel(R(0)) == (0, 0)
    assert R(0).cancel(R(0), include=False) == (1, 1, 0, 0)

    assert y.cancel(R(0)) == (1, 0)
    assert y.cancel(R(0), include=False) == (1, 1, 1, 0)

    assert R(0).cancel(y) == (0, 1)
    assert R(0).cancel(y, include=False) == (1, 1, 0, 1)

    assert (y**2 - x**2).cancel(y - x) == (x + y, 1)

    f = 2*x**3 + 4*x**2 + 2*x
    g = 3*x**2 + 3*x
    F = 2*x + 2
    G = 3

    assert f.cancel(g) == (F, G)

    assert (-f).cancel(g) == (-F, G)
    assert f.cancel(-g) == (-F, G)

    R, x, y = ring('x y', QQ)

    f = x**3/2 + x**2 + x/2
    g = x**2/3 + x/3
    F = 3*x + 3
    G = 2

    assert f.cancel(g) == (F, G)

    assert (-f).cancel(g) == (-F, G)
    assert f.cancel(-g) == (-F, G)


def test_sympyissue_10996():
    _, x, y, z = ring('x y z', ZZ)

    f = 12*x**6*y**7*z**3 - 3*x**4*y**9*z**3 + 12*x**3*y**5*z**4
    g = (-48*x**7*y**8*z**3 + 12*x**5*y**10*z**3 - 48*x**5*y**7*z**2 +
         36*x**4*y**7*z - 48*x**4*y**6*z**4 + 12*x**3*y**9*z**2 -
         48*x**3*y**4 - 9*x**2*y**9*z - 48*x**2*y**5*z**3 + 12*x*y**6 +
         36*x*y**5*z**2 - 48*y**2*z)

    H, cff, cfg = f.cofactors(g)

    assert H == 12*x**3*y**4 - 3*x*y**6 + 12*y**2*z
    assert H*cff == f
    assert H*cfg == g


def test_sympyissue_21460():
    R = ZZ.inject('x')

    r = R.gcd(R(4), R(6))
    assert type(r) is R.dtype
    assert r == 2

    R = QQ.inject('x')

    r = R.gcd(R(4), R(6))
    assert type(r) is R.dtype
    assert r == 1


@pytest.mark.slow
def test_sympyissue_23479():
    xs = ' '.join(f'x{i}' for i in range(13))
    _, x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12 = ring(xs, ZZ)

    p1 = (917280*x0*x4*x12**933 - 5640390*x0*x4*x12**905 + 14354685*x0*x4*x12**877 -
          19360440*x0*x4*x12**849 + 14597910*x0*x4*x12**821 -
          5834970*x0*x4*x12**793 + 965925*x0*x4*x12**765 -
          970200*x1*x2*x5*x11**493 + 5904990*x1*x2*x5*x11**479 -
          14883885*x1*x2*x5*x11**465 + 19889640*x1*x2*x5*x11**451 -
          14862510*x1*x2*x5*x11**437 + 5887890*x1*x2*x5*x11**423 -
          965925*x1*x2*x5*x11**409 - 34300*x1*x2*x7*x11**507 +
          224420*x1*x2*x7*x11**493 - 607600*x1*x2*x7*x11**479 +
          872200*x1*x2*x7*x11**465 - 700700*x1*x2*x7*x11**451 +
          298900*x1*x2*x7*x11**437 - 52920*x1*x2*x7*x11**423 -
          12005*x1*x2*x10*x11**521 + 72030*x1*x2*x10*x11**507 -
          180075*x1*x2*x10*x11**493 + 240100*x1*x2*x10*x11**479 -
          180075*x1*x2*x10*x11**465 + 72030*x1*x2*x10*x11**451 -
          12005*x1*x2*x10*x11**437 - 929160*x1*x3*x6*x11**489 +
          5699790*x1*x3*x6*x11**475 - 14473485*x1*x3*x6*x11**461 +
          19479240*x1*x3*x6*x11**447 - 14657310*x1*x3*x6*x11**433 +
          5846850*x1*x3*x6*x11**419 - 965925*x1*x3*x6*x11**405 -
          2178*x1*x3*x8*x11**517 + 13068*x1*x3*x8*x11**503 -
          32670*x1*x3*x8*x11**489 + 43560*x1*x3*x8*x11**475 -
          32670*x1*x3*x8*x11**461 + 13068*x1*x3*x8*x11**447 -
          2178*x1*x3*x8*x11**433 - 4356*x1*x3*x9*x11**503 +
          33660*x1*x3*x9*x11**489 - 102960*x1*x3*x9*x11**475 +
          162360*x1*x3*x9*x11**461 - 140580*x1*x3*x9*x11**447 +
          63756*x1*x3*x9*x11**433 - 11880*x1*x3*x9*x11**419)

    p2 = (3889620*x0*x12**1129 - 23337720*x0*x12**1101 + 58344300*x0*x12**1073 -
          77792400*x0*x12**1045 + 58344300*x0*x12**1017 - 23337720*x0*x12**989 +
          3889620*x0*x12**961)

    assert p1.gcd(p2) == 1
