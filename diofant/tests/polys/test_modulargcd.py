from diofant import QQ, ZZ, ring, sqrt
from diofant.polys.modulargcd import (_chinese_remainder_reconstruction,
                                      _func_field_modgcd_m, _to_ANP_poly,
                                      _to_ZZ_poly)
from diofant.polys.polyconfig import using


__all__ = ()


def test_modgcd_univariate_integers():
    R, x = ring('x', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g = R.zero, R.zero

        assert f.cofactors(g) == (0, 0, 0)

        f, g = R.zero, x

        assert f.cofactors(g) == (x, 0, 1)
        assert g.cofactors(f) == (x, 1, 0)

        f, g = R.zero, -x

        assert f.cofactors(g) == (x, 0, -1)
        assert g.cofactors(f) == (x, -1, 0)

        f, g = 2*x, R(2)

        assert f.cofactors(g) == (2, x, 1)

        f, g = 2*x + 2, 6*x**2 - 6

        assert f.cofactors(g) == (2*x + 2, 1, 3*x - 3)

        f, g = [1000000000000*x + 998549000000]*2

        assert f.cofactors(g) == (f, 1, 1)

        f, g = 999530000000*x + 1000000000000, 999530000000*x + 999999000000

        assert f.cofactors(g) == (1000000, 999530*x + 1000000, 999530*x + 999999)

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
        g = (365431878023781158602430064717380211405897160759702125019136*x**21 +
             197599133478719444145775798221171663643171734081650688*x**14 -
             9504116979659010018253915765478924103928886144*x**7 -
             311973482284542371301330321821976049)

        assert f.gcd(f.diff(x)) == g

        f = 1317378933230047068160*x + 2945748836994210856960
        g = 120352542776360960*x + 269116466014453760

        assert f.cofactors(g) == (120352542776360960*x + 269116466014453760,
                                  10946, 1)


def test_modgcd_bivariate_integers():
    R, x, y = ring('x,y', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g = R.zero, R.zero

        assert f.cofactors(g) == (0, 0, 0)

        f, g = 2*x, R(2)

        assert f.cofactors(g) == (2, x, 1)

        f, g = x + 2*y, x + y

        assert f.cofactors(g) == (1, f, g)

        f, g = x**2 + 2*x*y + y**2, x**3 + y**3

        assert f.cofactors(g) == (x + y, x + y, x**2 - x*y + y**2)

        f, g = x*y**2 + 2*x*y + x, x*y**3 + x

        assert f.cofactors(g) == (x*y + x, y + 1, y**2 - y + 1)

        f, g = x**2*y**2 + x**2*y + 1, x*y**2 + x*y + 1

        assert f.cofactors(g) == (1, f, g)

        f = 2*x*y**2 + 4*x*y + 2*x + y**2 + 2*y + 1
        g = 2*x*y**3 + 2*x + y**3 + 1

        assert f.cofactors(g) == (2*x*y + 2*x + y + 1, y + 1, y**2 - y + 1)

        f, g = 2*x**2 + 4*x + 2, x + 1

        assert f.cofactors(g) == (x + 1, 2*x + 2, 1)

        f, g = x + 1, 2*x**2 + 4*x + 2

        assert f.cofactors(g) == (x + 1, 1, 2*x + 2)

        f = 2*x**2 + 4*x*y - 2*x - 4*y
        g = x**2 + x - 2

        assert f.cofactors(g) == (x - 1, 2*x + 4*y, x + 2)

        f = 2*x**2 + 2*x*y - 3*x - 3*y
        g = 4*x*y - 2*x + 4*y**2 - 2*y

        assert f.cofactors(g) == (x + y, 2*x - 3, 4*y - 2)


def test_chinese_remainder():
    R, x, y = ring('x, y', ZZ)
    p, q = 3, 5

    hp = x**3*y - x**2 - 1
    hq = -x**3*y - 2*x*y**2 + 2

    hpq = _chinese_remainder_reconstruction(hp, hq, p, q)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq

    T, z = ring('z', R)
    p, q = 3, 7

    hp = (x*y + 1)*z**2 + x
    hq = (x**2 - 3*y)*z + 2

    hpq = _chinese_remainder_reconstruction(hp, hq, p, q)

    assert hpq.trunc_ground(p) == hp
    assert hpq.trunc_ground(q) == hq


def test_modgcd_multivariate_integers():
    R, x, y = ring('x,y', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g = R.zero, R.zero

        assert f.cofactors(g) == (0, 0, 0)

        f, g = 2*x**2 + 4*x + 2, x + 1

        assert f.cofactors(g) == (x + 1, 2*x + 2, 1)

        f, g = x + 1, 2*x**2 + 4*x + 2

        assert f.cofactors(g) == (x + 1, 1, 2*x + 2)

        f = 2*x**2 + 2*x*y - 3*x - 3*y
        g = 4*x*y - 2*x + 4*y**2 - 2*y

        assert f.cofactors(g) == (x + y, 2*x - 3, 4*y - 2)

        f, g = x*y**2 + 2*x*y + x, x*y**3 + x

        assert f.cofactors(g) == (x*y + x, y + 1, y**2 - y + 1)

        f, g = x**2*y**2 + x**2*y + 1, x*y**2 + x*y + 1

        assert f.cofactors(g) == (1, f, g)

        f = x**4 + 8*x**3 + 21*x**2 + 22*x + 8
        g = x**3 + 6*x**2 + 11*x + 6

        assert f.cofactors(g) == (x**2 + 3*x + 2, x**2 + 5*x + 4, x + 3)

    R, x, y, z, u = ring('x,y,z,u', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g = x + y + z, -x - y - z - u

        assert f.cofactors(g) == (1, f, g)

        f, g = u**2 + 2*u + 1, 2*u + 2

        assert f.cofactors(g) == (u + 1, u + 1, 2)

        f, g = z**2*u**2 + 2*z**2*u + z**2 + z*u + z, u**2 + 2*u + 1
        h, cff, cfg = u + 1, z**2*u + z**2 + z, u + 1

        assert f.cofactors(g) == (h, cff, cfg)
        assert g.cofactors(f) == (h, cfg, cff)

        f, g, h = (1199999999999991*x**17 - y, 2*y - 19989798798 + x**211,
                   12*x*y**7 + x**4 - 1)

        for i in range(10):
            assert (f*h).cofactors(g*h) == (h, f, g)

    R, x, y, z = ring('x,y,z', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g = x - y*z, x - y*z

        assert f.cofactors(g) == (x - y*z, 1, 1)

        f, g, h = R.fateman_poly_F_1()
        H, cff, cfg = f.cofactors(g)

        assert H == h and H*cff == f and H*cfg == g

        f, g, h = R.fateman_poly_F_2()
        H, cff, cfg = f.cofactors(g)

        assert H == h and H*cff == f and H*cfg == g

        f, g, h = R.fateman_poly_F_3()
        H, cff, cfg = f.cofactors(g)

        assert H == h and H*cff == f and H*cfg == g

    R, x, y, z, u, v = ring('x,y,z,u,v', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g, h = R.fateman_poly_F_1()
        H, cff, cfg = f.cofactors(g)

        assert H == h and H*cff == f and H*cfg == g

    R, x, y, z, u, v, a, b = ring('x,y,z,u,v,a,b', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g, h = R.fateman_poly_F_1()
        H, cff, cfg = f.cofactors(g)

        assert H == h and H*cff == f and H*cfg == g

    R, x, y, z, u, v, a, b, c, d = ring('x,y,z,u,v,a,b,c,d', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g, h = R.fateman_poly_F_1()
        H, cff, cfg = f.cofactors(g)

        assert H == h and H*cff == f and H*cfg == g

    R, x, y, z, t = ring('x,y,z,t', ZZ)

    with using(use_heu_gcd=False, fallback_gcd_zz_method='modgcd'):
        f, g, h = R.fateman_poly_F_3()
        H, cff, cfg = f.cofactors(g)

        assert H == h and H*cff == f and H*cfg == g


def test_to_ZZ_ANP_poly():
    A = QQ.algebraic_field(sqrt(2))
    R, x = ring('x', A)
    f = x*(sqrt(2) + 1)

    T, x_, z_ = ring('x_, z_', ZZ)
    f_ = x_*z_ + x_

    assert _to_ZZ_poly(f, T) == f_
    assert _to_ANP_poly(f_, R) == f

    R, x, t, s = ring('x, t, s', A)
    f = x*t**2 + x*s + sqrt(2)

    D, t_, s_ = ring('t_, s_', ZZ)
    T, x_, z_ = ring('x_, z_', D)
    f_ = (t_**2 + s_)*x_ + z_

    assert _to_ZZ_poly(f, T) == f_
    assert _to_ANP_poly(f_, R) == f


def test_modgcd_algebraic_field():
    A = QQ.algebraic_field(sqrt(2))
    R, x = ring('x', A)

    with using(gcd_aa_method='modgcd'):
        f, g = 2*x, R(2)

        assert f.cofactors(g) == (1, f, g)

        f, g = 2*x, R(sqrt(2))

        assert f.cofactors(g) == (1, f, g)

        f, g = 2*x + 2, 6*x**2 - 6

        assert f.cofactors(g) == (x + 1, 2, 6*x - 6)

    R, x, y = ring('x, y', A)

    with using(gcd_aa_method='modgcd'):
        f, g = x + sqrt(2)*y, x + y

        assert f.cofactors(g) == (1, f, g)

        f, g = x*y + sqrt(2)*y**2, sqrt(2)*y

        assert f.cofactors(g) == (y, x + sqrt(2)*y, sqrt(2))

        f, g = x**2 + 2*sqrt(2)*x*y + 2*y**2, x + sqrt(2)*y

        assert f.cofactors(g) == (g, g, 1)

    A = QQ.algebraic_field(sqrt(2), sqrt(3))
    R, x, y, z = ring('x, y, z', A)

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

    A = QQ.algebraic_field(sqrt(2)**(-1)*sqrt(3))
    R, x = ring('x', A)

    with using(gcd_aa_method='modgcd'):
        f, g = x + 1, x - 1

        assert f.cofactors(g) == (1, f, g)


def test_modgcd_algebraic_field_random():
    A = QQ.algebraic_field(sqrt(2), sqrt(3))
    R, x, y, z = ring('x, y, z', A)

    with using(gcd_aa_method='modgcd'):
        h = x**2*y**3 + 1111*sqrt(6)/12*z
        a, b = 11*y**3 + 2, (y + x - 1)*h

        assert (h*a).cofactors(h*b) == (h, a, b)

        a, b = 12*y + 2*x - 1, (y + x - 1)*h

        assert (h*a).cofactors(h*b) == (h, a, b)


def test_modgcd_func_field():
    D, t = ring('t', ZZ)
    R, x, z = ring('x, z', D)

    minpoly = (z**2*t**2 + z**2*t - 1).drop(0)
    f, g = x + 1, x - 1

    assert _func_field_modgcd_m(f, g, minpoly) == R.one

    # First example from Monagan2004algebraic.
    m = z**2 - t
    f = 3*t*x**2 - (2*t**2 - 3*t)*x*z + 15*x + 15*z - 2*t**3
    g = 3*t*x**2*z + 15*x*z + (-2*t**3 + 3*t)*x - 2*t**2*z + 15

    assert _func_field_modgcd_m(f, g, m.drop(0)) == 3*t*x - 2*t**2*z + 15

    g = 3*t*x - 2*t**2*z + 15
    a = x + z
    b = x*z + 1

    assert _func_field_modgcd_m(a*g, b*g, m.drop(0)) == g % m
    assert _func_field_modgcd_m(a*g**2, b*g**2, m.drop(0)) == g**2 % m

    # issue diofant/diofant#850
    assert _func_field_modgcd_m(a*g**3, b*g**3, m.drop(0)) == g**3 % m
