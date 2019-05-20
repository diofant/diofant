"""Tests for Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences. """

import pytest

from diofant.domains import CC, FF, QQ, RR, ZZ
from diofant.functions import sqrt
from diofant.polys.polyconfig import using
from diofant.polys.polyerrors import HeuristicGCDFailed, NotInvertible
from diofant.polys.rings import ring
from diofant.polys.specialpolys import (dmp_fateman_poly_F_1,
                                        dmp_fateman_poly_F_2,
                                        dmp_fateman_poly_F_3, f_polys)


__all__ = ()

f_0, f_1, f_2, f_3, f_4, f_5, f_6 = f_polys()


def test_dup_gcdex():
    R, x = ring("x", QQ)

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

    R, x = ring("x", FF(11))

    assert R.zero.gcdex(R(2)) == (0, 6, 1)
    assert R(2).gcdex(R(2)) == (0, 6, 1)

    assert R.zero.gcdex(3*x) == (0, 4, x)

    assert (3*x).gcdex(3*x) == (0, 4, x)

    assert (x**2 + 8*x + 7).gcdex(x**3 + 7*x**2 + x + 7) == (5*x + 6, 6, x + 7)


def test_dup_invert():
    R, x = ring("x", QQ)
    assert R.dup_invert(2*x, x**2 - 16) == x/32
    pytest.raises(NotInvertible, lambda: R.dup_invert(x**2 - 1, x - 1))


def test_dmp_prem():
    R, x = ring('x', ZZ)

    f = 3*x**3 + x**2 + x + 5
    g = 5*x**2 - 3*x + 1

    r = 52*x + 111

    assert f.prem(g) == r

    pytest.raises(ZeroDivisionError, lambda: f.prem(0))

    f = x**2 + 1
    g = 2*x - 4
    r = 20

    assert f.prem(g) == r

    R, x = ring('x', QQ)

    f = 3*x**3 + x**2 + x + 5
    g = 5*x**2 - 3*x + 1

    r = 52*x + 111

    assert g.prem(f) == g
    assert f.prem(g) == r

    R, x, y = ring('x y', ZZ)

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
    R, x = ring("x", ZZ)

    assert R(0).resultant(R(0)) == 0

    assert R(1).resultant(R(0)) == 0
    assert R(0).resultant(R(1)) == 0

    f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    a = 15*x**4 - 3*x**2 + 9
    b = 65*x**2 + 125*x - 245
    c = 9326*x - 12300
    d = 260708

    assert f.subresultants(g) == [f, g, a, b, c, d]
    assert f.resultant(g) == R.dmp_LC(d)

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

    assert R.dmp_inner_subresultants(0, 0) == ([], [])
    assert R.dmp_inner_subresultants(0, 1) == ([1], [1])

    R, x, y = ring("x,y", ZZ)

    assert R(0).resultant(R(0)) == 0
    assert R(0).resultant(R(0), includePRS=True)[0] == 0
    assert R.dmp_zz_collins_resultant(0, 0) == 0
    assert R.dmp_qq_collins_resultant(0, 0) == 0

    assert R(1).resultant(R(0)) == 0

    assert R(0).resultant(R(1)) == 0
    assert R(0).resultant(R(1), includePRS=True)[0] == 0
    assert R.dmp_zz_collins_resultant(0, 1) == 0
    assert R.dmp_qq_collins_resultant(0, 1) == 0

    assert R.dmp_inner_subresultants(0, 0) == ([], [])
    assert R.dmp_inner_subresultants(0, 1) == ([R.one], [[1]])

    f = 3*x**2*y - y**3 - 4
    g = x**2 + x*y**3 - 9

    a = 3*x*y**4 + y**3 - 27*y + 4
    b = -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16

    r = R.dmp_LC(b)
    rr = (r, [3*x**2*y - y**3 - 4, x**2 + x*y**3 - 9, 3*x*y**4 + y**3 - 27*y + 4,
              -3*y**10 - 12*y**7 + y**6 - 54*y**4 + 8*y**3 + 729*y**2 - 216*y + 16])

    assert f.subresultants(g) == [f, g, a, b]

    assert f.resultant(g) == r
    assert f.resultant(g, includePRS=True) == rr
    assert R.dmp_zz_collins_resultant(f, g) == r
    assert R.dmp_qq_collins_resultant(f, g) == r

    f = -x**3 + 5
    g = 3*x**2*y + x**2

    a = 45*y**2 + 30*y + 5
    b = 675*y**3 + 675*y**2 + 225*y + 25

    r = R.dmp_LC(b)

    assert f.subresultants(g) == [f, g, a]
    assert f.resultant(g) == r
    assert f.resultant(g, includePRS=True)[0] == r
    assert R.dmp_zz_collins_resultant(f, g) == r
    assert R.dmp_qq_collins_resultant(f, g) == r

    R, x, y, z, u, v = ring("x,y,z,u,v", ZZ)

    f = 6*x**2 - 3*x*y - 2*x*z + y*z
    g = x**2 - x*u - x*v + u*v

    r = y**2*z**2 - 3*y**2*z*u - 3*y**2*z*v + 9*y**2*u*v - 2*y*z**2*u \
        - 2*y*z**2*v + 6*y*z*u**2 + 12*y*z*u*v + 6*y*z*v**2 - 18*y*u**2*v \
        - 18*y*u*v**2 + 4*z**2*u*v - 12*z*u**2*v - 12*z*u*v**2 + 36*u**2*v**2

    assert R.dmp_zz_collins_resultant(f, g) == r.drop(x)

    R,  x, y, z, u, v = ring("x,y,z,u,v", QQ)

    f = x**2 - x*y/2 - x*z/3 + y*z/6
    g = x**2 - x*u - x*v + u*v

    r = y**2*z**2/36 - y**2*z*u/12 - y**2*z*v/12 + y**2*u*v/4 \
        - y*z**2*u/18 - y*z**2*v/18 + y*z*u**2/6 + y*z*u*v/3 \
        + y*z*v**2/6 - y*u**2*v/2 - y*u*v**2/2 + z**2*u*v/9 \
        - z*u**2*v/3 - z*u*v**2/3 + u**2*v**2

    assert R.dmp_qq_collins_resultant(f, g) == r.drop(x)

    Rt, t = ring("t", ZZ)
    Rx, x = ring("x", Rt)

    f = x**6 - 5*x**4 + 5*x**2 + 4
    g = -6*t*x**5 + x**4 + 20*t*x**3 - 3*x**2 - 10*t*x + 6

    assert f.resultant(g) == 2930944*t**6 + 2198208*t**4 + 549552*t**2 + 45796

    assert (x - 1).resultant(x + 1, includePRS=True) == (2, [x - 1, x + 1, 2])

    R, x, y = ring("x,y", ZZ)

    f = x + y
    g = x**2 - x*y + 1

    assert f.resultant(g) == (1 + 2*y**2).drop(x)

    g += 1
    with using(use_collins_resultant=True):
        assert f.resultant(g) == (2 + 2*y**2).drop(x)

    R, x, y = ring("x,y", QQ)

    f = x + y
    g = x**2 - x*y + 1

    with using(use_collins_resultant=True):
        assert f.resultant(g) == (1 + 2*y**2).drop(x)

    R, x, y = ring("x y", ZZ)

    f = x + y + 2
    g = 2*x*y + x + 3
    assert R.dmp_zz_collins_resultant(f, g) == (-2*y**2 - 5*y + 1).drop(x)


def test_PolyElement_discriminant():
    R, x = ring("x", ZZ)

    assert R(0).discriminant() == 0
    assert x.discriminant() == 1

    assert (x**3 + 3*x**2 + 9*x - 13).discriminant() == -11664
    assert (5*x**5 + x**3 + 2).discriminant() == 31252160
    assert (x**4 + 2*x**3 + 6*x**2 - 22*x + 13).discriminant() == 0
    assert (12*x**7 + 15*x**4 + 30*x**3 + x**2 + 1).discriminant() == -220289699947514112

    assert (x**2 + 2*x + 3).discriminant() == -8

    R, x, y = ring("x,y", ZZ)

    assert R(0).discriminant() == 0
    assert y.discriminant() == 0

    assert (x**3 + 3*x**2 + 9*x - 13).discriminant() == -11664
    assert (5*x**5 + x**3 + 2).discriminant() == 31252160
    assert (x**4 + 2*x**3 + 6*x**2 - 22*x + 13).discriminant() == 0
    assert (12*x**7 + 15*x**4 + 30*x**3 + x**2 + 1).discriminant() == -220289699947514112

    assert (x**2*y + 2*y).discriminant() == (-8*y**2).drop(x)
    assert (x*y**2 + 2*x).discriminant() == 1

    R, x, y, z = ring("x,y,z", ZZ)
    assert (x*y + z).discriminant() == 1

    R, x, y, z, u = ring("x,y,z,u", ZZ)
    assert (x**2*y + x*z + u).discriminant() == (-4*y*u + z**2).drop(x)

    R, x, y, z, u, v = ring("x,y,z,u,v", ZZ)
    assert (x**3*y + x**2*z + x*u + v).discriminant() == \
        (-27*y**2*v**2 + 18*y*z*u*v - 4*y*u**3 - 4*z**3*v + z**2*u**2).drop(x)


def test_dmp_gcd():
    R, x = ring("x", ZZ)

    f, g = 0, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (0, 0, 0)

    f, g = 2, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, 1, 0)

    f, g = -2, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, -1, 0)

    f, g = 0, -2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, 0, -1)

    f, g = 0, 2*x + 4
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2*x + 4, 0, 1)

    f, g = 2*x + 4, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2*x + 4, 1, 0)

    f, g = 2, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, 1, 1)

    f, g = -2, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, -1, 1)

    f, g = 2, -2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, 1, -1)

    f, g = -2, -2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, -1, -1)

    f, g = x**2 + 2*x + 1, 1
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (1, x**2 + 2*x + 1, 1)

    f, g = x**2 + 2*x + 1, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (1, x**2 + 2*x + 1, 2)

    f, g = 2*x**2 + 4*x + 2, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, x**2 + 2*x + 1, 1)

    f, g = 2, 2*x**2 + 4*x + 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (2, 1, x**2 + 2*x + 1)

    f, g = 2*x**2 + 4*x + 2, x + 1
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (x + 1, 2*x + 2, 1)

    f, g = x + 1, 2*x**2 + 4*x + 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (x + 1, 1, 2*x + 2)

    f, g = x - 31, x
    assert R.dmp_zz_heu_gcd(f, g) == R.dup_rr_prs_gcd(f, g) == (1, f, g)

    f = x**4 + 8*x**3 + 21*x**2 + 22*x + 8
    g = x**3 + 6*x**2 + 11*x + 6

    h = x**2 + 3*x + 2

    cff = x**2 + 5*x + 4
    cfg = x + 3

    assert R.dmp_zz_heu_gcd(f, g) == (h, cff, cfg)
    assert R.dup_rr_prs_gcd(f, g) == (h, cff, cfg)

    f = x**4 - 4
    g = x**4 + 4*x**2 + 4

    h = x**2 + 2

    cff = x**2 - 2
    cfg = x**2 + 2

    assert R.dmp_zz_heu_gcd(f, g) == (h, cff, cfg)
    assert R.dup_rr_prs_gcd(f, g) == (h, cff, cfg)

    f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    h = 1

    cff = f
    cfg = g

    assert R.dmp_zz_heu_gcd(f, g) == (h, cff, cfg)
    assert R.dup_rr_prs_gcd(f, g) == (h, cff, cfg)

    R, x = ring("x", QQ)

    f = x**8 + x**6 - 3*x**4 - 3*x**3 + 8*x**2 + 2*x - 5
    g = 3*x**6 + 5*x**4 - 4*x**2 - 9*x + 21

    h = 1

    cff = f
    cfg = g

    assert R.dmp_qq_heu_gcd(f, g) == (h, cff, cfg)
    assert R.dup_ff_prs_gcd(f, g) == (h, cff, cfg)

    assert R.dup_ff_prs_gcd(R.zero, R.zero) == ([], [], [])

    R, x = ring("x", ZZ)

    f = - 352518131239247345597970242177235495263669787845475025293906825864749649589178600387510272*x**49 \
        + 46818041807522713962450042363465092040687472354933295397472942006618953623327997952*x**42 \
        + 378182690892293941192071663536490788434899030680411695933646320291525827756032*x**35 \
        + 112806468807371824947796775491032386836656074179286744191026149539708928*x**28 \
        - 12278371209708240950316872681744825481125965781519138077173235712*x**21 \
        + 289127344604779611146960547954288113529690984687482920704*x**14 \
        + 19007977035740498977629742919480623972236450681*x**7 \
        + 311973482284542371301330321821976049

    g = 365431878023781158602430064717380211405897160759702125019136*x**21 \
        + 197599133478719444145775798221171663643171734081650688*x**14 \
        - 9504116979659010018253915765478924103928886144*x**7 \
        - 311973482284542371301330321821976049

    assert R.dmp_zz_heu_gcd(f, R.dmp_diff_in(f, 1, 0))[0] == g
    assert R.dup_rr_prs_gcd(f, R.dmp_diff_in(f, 1, 0))[0] == g

    R, x = ring("x", QQ)

    f = x**2/2 + x + QQ(1, 2)
    g = x/2 + QQ(1, 2)

    h = x + 1

    assert R.dmp_qq_heu_gcd(f, g) == (h, g, QQ(1, 2))
    assert R.dup_ff_prs_gcd(f, g) == (h, g, QQ(1, 2))

    R, x = ring("x", ZZ)

    f = 1317378933230047068160*x + 2945748836994210856960
    g = 120352542776360960*x + 269116466014453760

    h = 120352542776360960*x + 269116466014453760
    cff = 10946
    cfg = 1

    assert R.dmp_zz_heu_gcd(f, g) == (h, cff, cfg)

    with using(heu_gcd_max=0):
        pytest.raises(HeuristicGCDFailed, lambda: R.dmp_zz_heu_gcd(f, g))

    R, x = ring("x", CC)
    f, g = (x**2 - 1, x**3 - 3*x + 2)
    assert R.dmp_inner_gcd(f, g) == (1, f, g)

    R, x, y = ring("x,y", CC)
    f, g = (x**2 - y, x**3 - y*x + 2)
    assert R.dmp_inner_gcd(f, g) == (1, f, g)

    R,  x, y = ring("x,y", ZZ)

    f, g = 0, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (0, 0, 0)

    f, g = 2, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, 1, 0)

    f, g = -2, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, -1, 0)

    f, g = 0, -2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, 0, -1)

    f, g = 0, 2*x + 4
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2*x + 4, 0, 1)

    f, g = 2*x + 4, 0
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2*x + 4, 1, 0)

    f, g = 2, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, 1, 1)

    f, g = -2, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, -1, 1)

    f, g = 2, -2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, 1, -1)

    f, g = -2, -2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, -1, -1)

    f, g = x**2 + 2*x + 1, 1
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (1, x**2 + 2*x + 1, 1)

    f, g = x**2 + 2*x + 1, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (1, x**2 + 2*x + 1, 2)
    with using(use_simplify_gcd=0):
        assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (1, x**2 + 2*x + 1, 2)

    f, g = 2*x**2 + 4*x + 2, 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, x**2 + 2*x + 1, 1)

    f, g = 2, 2*x**2 + 4*x + 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (2, 1, x**2 + 2*x + 1)

    f, g = 2*x**2 + 4*x + 2, x + 1
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (x + 1, 2*x + 2, 1)

    f, g = x + 1, 2*x**2 + 4*x + 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (x + 1, 1, 2*x + 2)

    with using(heu_gcd_max=0):
        pytest.raises(HeuristicGCDFailed, lambda: R.dmp_zz_heu_gcd(f, g))

    f = x**2 + 2*x*y + y**2
    g = x**2 + x*y

    assert R.dmp_rr_prs_gcd(f, g) == (x + y, x + y, x)

    R, x, y, z, u = ring("x,y,z,u", ZZ)

    f, g = u**2 + 2*u + 1, 2*u + 2
    assert R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) == (u + 1, u + 1, 2)

    f, g = z**2*u**2 + 2*z**2*u + z**2 + z*u + z, u**2 + 2*u + 1
    h, cff, cfg = u + 1, z**2*u + z**2 + z, u + 1

    assert R.dmp_zz_heu_gcd(f, g) == (h, cff, cfg)
    assert R.dmp_rr_prs_gcd(f, g) == (h, cff, cfg)

    assert R.dmp_zz_heu_gcd(g, f) == (h, cfg, cff)
    assert R.dmp_rr_prs_gcd(g, f) == (h, cfg, cff)

    R, x, y, z = ring("x,y,z", ZZ)

    f, g, h = map(R.from_dense, dmp_fateman_poly_F_1(2, ZZ))
    H, cff, cfg = R.dmp_zz_heu_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    H, cff, cfg = R.dmp_rr_prs_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    R, x, y, z, u, v = ring("x,y,z,u,v", ZZ)

    f, g, h = map(R.from_dense, dmp_fateman_poly_F_1(4, ZZ))
    H, cff, cfg = R.dmp_zz_heu_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    R, x, y, z, u, v, a, b = ring("x,y,z,u,v,a,b", ZZ)

    f, g, h = map(R.from_dense, dmp_fateman_poly_F_1(6, ZZ))
    H, cff, cfg = R.dmp_zz_heu_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    R, x, y, z, u, v, a, b, c, d = ring("x,y,z,u,v,a,b,c,d", ZZ)

    f, g, h = map(R.from_dense, dmp_fateman_poly_F_1(8, ZZ))
    H, cff, cfg = R.dmp_zz_heu_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    R, x, y, z = ring("x,y,z", ZZ)

    f, g, h = map(R.from_dense, dmp_fateman_poly_F_2(2, ZZ))
    H, cff, cfg = R.dmp_zz_heu_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    H, cff, cfg = R.dmp_rr_prs_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    f, g, h = map(R.from_dense, dmp_fateman_poly_F_3(2, ZZ))
    H, cff, cfg = R.dmp_zz_heu_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    H, cff, cfg = R.dmp_rr_prs_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    R, x, y, z, u, v = ring("x,y,z,u,v", ZZ)

    f, g, h = map(R.from_dense, dmp_fateman_poly_F_3(4, ZZ))
    H, cff, cfg = R.dmp_inner_gcd(f, g)

    assert H == h and R.dmp_mul(H, cff) == f \
        and R.dmp_mul(H, cfg) == g

    R, x, y = ring("x,y", QQ)

    f = x**2/2 + x + QQ(1, 2)
    g = x/2 + QQ(1, 2)

    h = x + 1

    assert R.dmp_qq_heu_gcd(f, g) == (h, g, QQ(1, 2))
    assert R.dmp_ff_prs_gcd(f, g) == (h, g, QQ(1, 2))
    with using(use_simplify_gcd=0):
        assert R.dmp_qq_heu_gcd(f, g) == (h, g, QQ(1, 2))
        assert R.dmp_ff_prs_gcd(f, g) == (h, g, QQ(1, 2))

    assert R.dmp_ff_prs_gcd(R.zero, R.zero) == (0, 0, 0)
    assert R.dmp_qq_heu_gcd(R.zero, R.zero) == (0, 0, 0)
    assert R.dmp_ff_prs_gcd(R.zero, g) == (x + 1, R.zero, QQ(1, 2))
    assert R.dmp_qq_heu_gcd(R.zero, g) == (x + 1, R.zero, QQ(1, 2))

    R, x, y = ring("x,y", RR)

    f = 2.1*x*y**2 - 2.2*x*y + 2.1*x
    g = 1.0*x**3

    assert R.dmp_ff_prs_gcd(f, g) == \
        (1.0*x, 2.1*y**2 - 2.2*y + 2.1, 1.0*x**2)

    R, x, y = ring("x,y", ZZ)

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

    assert (R.dmp_zz_heu_gcd(f, g) == R.dmp_rr_prs_gcd(f, g) ==
            (g,
             -1743436700916730000000000000000000000000000*x**2*y -
             12525091394814956778444367128798089188206978000000000000000*x*y -
             244093590929967254073813518342605644787785800*x -
             22495499028725631994927113634418779135935898997901327211111875586270479483*y -
             876801128965234839118530545935732755107147297241756982389990, 1))

    R, x = ring("x", ZZ)

    f, g = x**2 - 1, x**2 - 3*x + 2
    assert R.dmp_gcd(f, g) == x - 1

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        R.dmp_gcd(f, g) == x - 1

    R, x = ring("x", QQ)

    f, g = x**2/2 + x + QQ(1, 2), x/2 + QQ(1, 2)

    assert R.dmp_gcd(f, g) == x + 1
    with using(use_heu_gcd=False):
        R.dmp_gcd(f, g) == x + 1

    R, x, y = ring("x,y", QQ.algebraic_field(sqrt(2)))

    f, g = (x + sqrt(2)*y)**2, x + sqrt(2)*y

    assert R.dmp_gcd(f, g) == g
    with using(gcd_aa_method='modgcd'):
        assert R.dmp_gcd(f, g) == g


def test_PolyElement_lcm():
    R, x = ring("x", ZZ)

    assert R(2).lcm(R(6)) == 6

    assert (2*x**3).lcm(6*x) == 6*x**3
    assert (2*x**3).lcm(3*x) == 6*x**3

    assert (x**2 + x).lcm(x) == x**2 + x
    assert (x**2 + x).lcm(2*x) == 2*x**2 + 2*x
    assert (x**2 + 2*x).lcm(x) == x**2 + 2*x
    assert (2*x**2 + x).lcm(x) == 2*x**2 + x
    assert (2*x**2 + x).lcm(2*x) == 4*x**2 + 2*x
    assert (x**2 - 1).lcm(x**2 - 3*x + 2) == x**3 - 2*x**2 - x + 2

    R,  x, y = ring("x,y", ZZ)

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

    R, x = ring('x', QQ)

    f = (x**2 + 7*x/2 + 3)/2
    g = x**2/2 + x
    h = x**3 + 7/2*x**2 + 3*x

    assert f.lcm(g) == h

    R,  x, y = ring("x,y", QQ)

    f = 2*x*y - x**2/2 + QQ(1, 3)
    g = 3*x**3 - x*y**2 - QQ(1, 2)
    h = (x**5 - 4*x**4*y - x**3*y**2/3 - 2*x**3/3 + 4*x**2*y**3/3 -
         x**2/6 + 2*x*y**2/9 + 2*x*y/3 + QQ(1, 9))

    assert f.lcm(g) == h

    f = x**2/4 + x*y + y**2
    g = x**2/2 + x*y
    h = x**3 + 4*x**2*y + 4*x*y**2

    assert f.lcm(g) == h

    R, x = ring("x", FF(11))

    assert R.zero.lcm(R(2)) == 0
    assert R(2).lcm(R(2)) == 1

    assert R.zero.lcm(x) == 0

    assert (3*x).lcm(3*x) == x
    assert (x**2 + 8*x + 7).lcm(x**3 + 7*x**2 + x + 7) == (x**4 + 8*x**3 +
                                                           8*x**2 + 8*x + 7)

    R, x = ring("x", FF(5))

    assert (3*x**2 + 2*x + 4).lcm(2*x**2 + 2*x + 3) == x**3 + 2*x**2 + 4


def test_dmp_content():
    R,  x, y = ring("x,y", ZZ)

    assert R.dmp_content(-2) == 2

    f, g, F = 3*y**2 + 2*y + 1, 1, 0

    for i in range(5):
        g *= f
        F += x**i*g

    assert R.dmp_content(F) == f.drop(x)

    R,  x, y, z = ring("x,y,z", ZZ)

    assert R.dmp_content(f_4) == 1
    assert R.dmp_content(f_5) == 1

    R,  x, y, z, t = ring("x,y,z,t", ZZ)
    assert R.dmp_content(f_6) == 1


def test_dmp_primitive():
    R,  x, y = ring("x,y", ZZ)

    assert R.dmp_primitive(0) == (0, 0)
    assert R.dmp_primitive(1) == (1, 1)

    f, g, F = 3*y**2 + 2*y + 1, 1, 0

    for i in range(5):
        g *= f
        F += x**i*g

    assert R.dmp_primitive(F) == (f.drop(x), F // f)

    R,  x, y, z = ring("x,y,z", ZZ)

    cont, f = R.dmp_primitive(f_4)
    assert cont == 1 and f == f_4
    cont, f = R.dmp_primitive(f_5)
    assert cont == 1 and f == f_5

    R,  x, y, z, t = ring("x,y,z,t", ZZ)

    cont, f = R.dmp_primitive(f_6)
    assert cont == 1 and f == f_6


def test_PolyElement_cancel():
    R, x = ring("x", ZZ)

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

    R, x, y = ring("x,y", ZZ)

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

    R, x = ring('x', QQ)

    assert (x**2/4 - 1).cancel(x/2 - 1) == (x + 2, 2)


def test_dmp_zz_modular_resultant():
    R, x, y = ring("x y", ZZ)
    R1 = R.drop(x)

    f = x + y + 2
    g = 2*x*y + x + 3

    assert R.dmp_zz_modular_resultant(f, g, 5) == -2*R1.y**2 + 1
