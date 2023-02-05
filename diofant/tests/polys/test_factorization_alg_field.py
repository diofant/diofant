import random

import pytest

from diofant import QQ, ZZ, I, ring, root, sqrt
from diofant.polys.factorization_alg_field import (_distinct_prime_divisors,
                                                   _sqf_p, efactor)


def test__distinct_prime_divisors():
    s = [20, 6, 7]
    assert _distinct_prime_divisors(s, ZZ) == [5, 3, 7]

    s = [15, 18, 11]
    assert _distinct_prime_divisors(s, ZZ) == [5, 2, 11]


def test__sqf_p():
    *_, z = ring('x z', ZZ)

    assert _sqf_p(z**2, (z**2 - 2).drop(0), 2) is True


def test_efactor_1():
    R, x, y = ring('x y', QQ.algebraic_field(sqrt(2)))

    f = x**2 + sqrt(2)*y
    assert efactor(f) == (1, [(f, 1)])

    f1 = x + y
    f2 = x**2 + sqrt(2)*y
    f = f1 * f2

    assert efactor(f) == (1, [(f1, 1), (f2, 1)])
    assert efactor(f, save=False) == (1, [(f1, 1), (f2, 1)])

    f1 = x + 2**10*y
    f2 = x**2 + sqrt(2)*y
    f = f1 * f2

    assert efactor(f) == (1, [(f1, 1), (f2, 1)])

    R, x, y, z = ring('x y z', QQ.algebraic_field(sqrt(2)))

    f1 = x - sqrt(2)*z
    f = f1**2

    assert efactor(f) == (1, [(f1, 2)])

    A3 = QQ.algebraic_field(sqrt(3))
    R, x, y, z = ring('x y z', A3)

    f1 = z + 1
    f2 = 3*x*y**2/4 + sqrt(3)
    f3 = sqrt(3)*y*x**2 + 2*y + z
    f = f1 * f2**2 * f3

    lc2 = f2.LC
    lc3 = f3.LC

    assert efactor(f) == (lc2**2*lc3, [(f1, 1), (f2.quo_ground(lc2), 2),
                                       (f3.quo_ground(lc3), 1)])

    R, x, y = ring('x y', QQ.algebraic_field(I))

    f1 = x - I*y
    f2 = x + I*y
    f = f1*f2

    assert efactor(f) == (1, [(f1, 1), (f2, 1)])

    f1 = x*(y + 1) + 1
    f2 = x*(y + I) + 1
    f3 = x**2*(y - I) + 1
    f = f1*f2*f3

    assert efactor(f) == (1, [(f1, 1), (f2, 1), (f3, 1)])

    lc = R.domain(-2)
    f1 = x*(y - 3*I) + lc**(-1)
    f2 = x*(y + 2) + 1
    f3 = x*(y + I) + 1
    f = lc*f1*f2*f3

    assert efactor(f) == (lc, [(f1, 1), (f2, 1), (f3, 1)])

    a = sqrt(2)*(1 + I)/2
    A = QQ.algebraic_field(a)
    R, x, y = ring('x y', A)

    f1 = x - a**3*y
    f2 = x - a*y
    f3 = x + a**3*y
    f4 = x + a*y
    f = x**4 + y**4

    assert f1*f2*f3*f4 == f
    assert efactor(f) == (1, [(f1, 1), (f2, 1), (f3, 1), (f4, 1)])

    R, x, y, z = ring('x y z', QQ.algebraic_field(root(2, 5)))

    f1 = y
    f2 = x - y
    f3 = x**2 + root(2, 5)*y*z
    f = f1*f2*f3

    assert efactor(f) == (1, [(f1, 1), (f2, 1), (f3, 1)])

    R, x, y, z = ring('x y z', QQ.algebraic_field(sqrt(2), root(3, 3)))

    lc = R.domain.from_expr(root(3, 3))
    f1 = x - z*root(3, 3)**2/3
    f2 = x**2 + 2*y + sqrt(2)
    f = lc*f1*f2

    assert efactor(f) == (lc, [(f1, 1), (f2, 1)])

    a = (-1 + sqrt(5))/4 - I*sqrt((sqrt(5) + 5)/8)
    A = QQ.algebraic_field(a)
    a = A.unit
    R, x, y, z = ring('x y z', A)

    f1 = x**2 + 2*(a**3 + a**2 + a + 1)*x + a*z**2 + a**3*y + 12*a**2
    f2 = x**2 - 2*a*x - (a**3 + a**2 + a + 1)*z**2 + a**2*y + 12*a**3
    f = f1*f2

    assert efactor(f) == (1, [(f1, 1), (f2, 1)])

    R, x, y = ring('x y', QQ.algebraic_field(root(1, 5, 3)))
    A = R.domain

    a = A([QQ(-19125, 42722), QQ(23337, 21361), QQ(46350, 21361), QQ(17315, 21361)])
    b = A([QQ(-17355, 85444), QQ(-15120, 21361), QQ(-7870, 21361), QQ(45917, 85444)])
    c = A([QQ(5, 521), QQ(130, 521), QQ(650, 521), QQ(104, 521)])
    d = A([QQ(16196, 21361), QQ(-6645, 21361), QQ(-20200, 21361), QQ(-29909, 42722)])

    e = a*y + b*y**2 + x*c + x*y*d + x**2

    r = efactor(e)

    e1 = A([QQ(75, 82), QQ(-10, 41), QQ(-25, 41), QQ(-47, 41)])*y + x
    e2 = (x + A([QQ(-163, 1042), QQ(-35, 521), QQ(-175, 521),
                 QQ(465, 1042)])*y + A([QQ(5, 521), QQ(130, 521),
                                        QQ(650, 521), QQ(104, 521)]))

    assert r == (1, [(e1, 1), (e2, 1)])


def test_efactor_random():
    A3 = QQ.algebraic_field(sqrt(3))
    _, x, y, z, t = ring('x y z t', A3)

    f1 = x*y - sqrt(3)
    f2 = z*t + 1
    f3 = x**2 + 1
    f4 = x**2 + z*t
    f = f1*f2*f3*f4

    for seed in [0, 1, 2, 6]:
        random.seed(seed)
        assert efactor(f) == (1, [(f1, 1), (f2, 1), (f3, 1), (f4, 1)])


@pytest.mark.slow
def test_efactor_wang():
    a = (-1 + sqrt(5))/4 - I*sqrt((sqrt(5) + 5)/8)
    A = QQ.algebraic_field(a)
    a = A.unit
    _, x, y, z = ring('x y z', A)

    f1 = x**2 - 2*a**2*x + a**3*z**2 - (a**3 + a**2 + a + 1)*y + 12*a
    f2 = x**2 - 2*a**3*x + a**2*z**2 + a*y - 12*(a**3 + a**2 + a + 1)
    f3 = x**2 + 2*(a**3 + a**2 + a + 1)*x + a*z**2 + a**3*y + 12*a**2
    f4 = x**2 - 2*a*x - (a**3 + a**2 + a + 1)*z**2 + a**2*y + 12*a**3
    f = f1*f2*f3*f4

    assert efactor(f) == (1, [(f1, 1), (f2, 1), (f3, 1), (f4, 1)])


@pytest.mark.timeout(60)
def test_sympyissue_19196():
    A = QQ.algebraic_field(sqrt(2), root(3, 3))
    _, x, y, z = ring('x y z', A)

    f1, f2 = x - z/root(3, 3), x**2 + 2*y + sqrt(2)
    f = f1*f2
    assert f.factor_list() == (1, [(f1, 1), (f2, 1)])
