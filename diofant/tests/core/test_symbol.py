import pytest

from diofant import (Dummy, Float, GreaterThan, I, Integer, LessThan, Rational,
                     StrictGreaterThan, StrictLessThan, Symbol, Wild, pi, sstr,
                     symbols, sympify)


__all__ = ()


def test_Symbol():
    a = Symbol('a')
    x1 = Symbol('x')
    x2 = Symbol('x')
    xdummy1 = Dummy('x')
    xdummy2 = Dummy('x')

    assert a != x1
    assert a != x2
    assert x1 == x2
    assert x1 != xdummy1
    assert xdummy1 != xdummy2

    assert Symbol('x') == Symbol('x')
    assert Dummy('x') != Dummy('x')
    d = symbols('d', cls=Dummy)
    assert isinstance(d, Dummy)
    c, d = symbols('c,d', cls=Dummy)
    assert isinstance(c, Dummy)
    assert isinstance(d, Dummy)
    pytest.raises(TypeError, Symbol)
    pytest.raises(TypeError, lambda: Symbol(1))


def test_Dummy():
    assert Dummy() != Dummy()
    Dummy._count = 0
    d1 = Dummy()
    Dummy._count = 0
    assert d1 == Dummy()


def test_as_dummy():
    x = Symbol('x')
    x1 = x.as_dummy()
    assert x1 != x
    assert x1 != x.as_dummy()

    x = Symbol('x', commutative=False)
    x1 = x.as_dummy()
    assert x1 != x
    assert x1.is_commutative is False


def test_lt_gt():
    x, y = Symbol('x'), Symbol('y')

    assert (x >= y) == GreaterThan(x, y)
    assert (x >= 0) == GreaterThan(x, 0)
    assert (x <= y) == LessThan(x, y)
    assert (x <= 0) == LessThan(x, 0)

    assert (0 <= x) == GreaterThan(x, 0)
    assert (0 >= x) == LessThan(x, 0)
    assert (Integer(0) >= x) == GreaterThan(0, x)
    assert (Integer(0) <= x) == LessThan(0, x)

    assert (x > y) == StrictGreaterThan(x, y)
    assert (x > 0) == StrictGreaterThan(x, 0)
    assert (x < y) == StrictLessThan(x, y)
    assert (x < 0) == StrictLessThan(x, 0)

    assert (0 < x) == StrictGreaterThan(x, 0)
    assert (0 > x) == StrictLessThan(x, 0)
    assert (Integer(0) > x) == StrictGreaterThan(0, x)
    assert (Integer(0) < x) == StrictLessThan(0, x)

    e = x**2 + 4*x + 1
    assert (e >= 0) == GreaterThan(e, 0)
    assert (0 <= e) == GreaterThan(e, 0)
    assert (e > 0) == StrictGreaterThan(e, 0)
    assert (0 < e) == StrictGreaterThan(e, 0)

    assert (e <= 0) == LessThan(e, 0)
    assert (0 >= e) == LessThan(e, 0)
    assert (e < 0) == StrictLessThan(e, 0)
    assert (0 > e) == StrictLessThan(e, 0)

    assert (Integer(0) >= e) == GreaterThan(0, e)
    assert (Integer(0) <= e) == LessThan(0, e)
    assert (Integer(0) < e) == StrictLessThan(0, e)
    assert (Integer(0) > e) == StrictGreaterThan(0, e)


def test_no_len():
    # there should be no len for numbers
    x = Symbol('x')
    pytest.raises(TypeError, lambda: len(x))


def test_ineq_unequal():
    x, y, z = symbols('x,y,z')

    e = (
        Integer(-1) >= x, Integer(-1) >= y, Integer(-1) >= z,
        Integer(-1) > x, Integer(-1) > y, Integer(-1) > z,
        Integer(-1) <= x, Integer(-1) <= y, Integer(-1) <= z,
        Integer(-1) < x, Integer(-1) < y, Integer(-1) < z,
        Integer(0) >= x, Integer(0) >= y, Integer(0) >= z,
        Integer(0) > x, Integer(0) > y, Integer(0) > z,
        Integer(0) <= x, Integer(0) <= y, Integer(0) <= z,
        Integer(0) < x, Integer(0) < y, Integer(0) < z,
        Rational(3, 7) >= x, Rational(3, 7) >= y, Rational(3, 7) >= z,
        Rational(3, 7) > x, Rational(3, 7) > y, Rational(3, 7) > z,
        Rational(3, 7) <= x, Rational(3, 7) <= y, Rational(3, 7) <= z,
        Rational(3, 7) < x, Rational(3, 7) < y, Rational(3, 7) < z,
        Float(1.5) >= x, Float(1.5) >= y, Float(1.5) >= z,
        Float(1.5) > x, Float(1.5) > y, Float(1.5) > z,
        Float(1.5) <= x, Float(1.5) <= y, Float(1.5) <= z,
        Float(1.5) < x, Float(1.5) < y, Float(1.5) < z,
        Integer(2) >= x, Integer(2) >= y, Integer(2) >= z,
        Integer(2) > x, Integer(2) > y, Integer(2) > z,
        Integer(2) <= x, Integer(2) <= y, Integer(2) <= z,
        Integer(2) < x, Integer(2) < y, Integer(2) < z,
        x >= -1, y >= -1, z >= -1,
        x > -1, y > -1, z > -1,
        x <= -1, y <= -1, z <= -1,
        x < -1, y < -1, z < -1,
        x >= 0, y >= 0, z >= 0,
        x > 0, y > 0, z > 0,
        x <= 0, y <= 0, z <= 0,
        x < 0, y < 0, z < 0,
        x >= 1.5, y >= 1.5, z >= 1.5,
        x > 1.5, y > 1.5, z > 1.5,
        x <= 1.5, y <= 1.5, z <= 1.5,
        x < 1.5, y < 1.5, z < 1.5,
        x >= 2, y >= 2, z >= 2,
        x > 2, y > 2, z > 2,
        x <= 2, y <= 2, z <= 2,
        x < 2, y < 2, z < 2,

        x >= y, x >= z, y >= x, y >= z, z >= x, z >= y,
        x > y, x > z, y > x, y > z, z > x, z > y,
        x <= y, x <= z, y <= x, y <= z, z <= x, z <= y,
        x < y, x < z, y < x, y < z, z < x, z < y,

        x - pi >= y + z, y - pi >= x + z, z - pi >= x + y,
        x - pi > y + z, y - pi > x + z, z - pi > x + y,
        x - pi <= y + z, y - pi <= x + z, z - pi <= x + y,
        x - pi < y + z, y - pi < x + z, z - pi < x + y,
        True, False)

    left_e = e[:-1]
    for i, e1 in enumerate(left_e):
        for e2 in e[i + 1:]:
            assert e1 != e2


def test_Wild_properties():
    # these tests only include Atoms
    x = Symbol('x')
    y = Symbol('y')
    p = Symbol('p', positive=True)
    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True, positive=True)

    given_patterns = [x, y, p, k, -k, n, -n, sympify(-3), sympify(3),
                      pi, Rational(3, 2), I]

    def integerp(k):
        return k.is_integer

    def positivep(k):
        return k.is_positive

    def symbolp(k):
        return k.is_Symbol

    def realp(k):
        return k.is_extended_real

    S = Wild('S', properties=[symbolp])
    R = Wild('R', properties=[realp])
    Y = Wild('Y', exclude=[x, p, k, n])
    P = Wild('P', properties=[positivep])
    K = Wild('K', properties=[integerp])
    N = Wild('N', properties=[positivep, integerp])

    given_wildcards = [S, R, Y, P, K, N]

    goodmatch = {
        S: (x, y, p, k, n),
        R: (p, k, -k, n, -n, -3, 3, pi, Rational(3, 2)),
        Y: (y, -3, 3, pi, Rational(3, 2), I),
        P: (p, n, 3, pi, Rational(3, 2)),
        K: (k, -k, n, -n, -3, 3),
        N: (n, 3)}

    for A in given_wildcards:
        for pat in given_patterns:
            d = pat.match(A)
            if pat in goodmatch[A]:
                assert d[A] in goodmatch[A]
            else:
                assert d is None


def test_symbols():
    x = Symbol('x')
    y = Symbol('y')
    z = Symbol('z')

    assert symbols('x') == x
    assert symbols('x ') == x
    assert symbols(' x ') == x
    assert symbols('x,') == (x,)
    assert symbols('x, ') == (x,)
    assert symbols('x ,') == (x,)

    assert symbols('x , y') == (x, y)

    assert symbols('x,y,z') == (x, y, z)
    assert symbols('x y z') == (x, y, z)

    assert symbols('x,y,z,') == (x, y, z)
    assert symbols('x y z ') == (x, y, z)

    xyz = Symbol('xyz')
    abc = Symbol('abc')

    assert symbols('xyz') == xyz
    assert symbols('xyz,') == (xyz,)
    assert symbols('xyz,abc') == (xyz, abc)

    assert symbols(('xyz',)) == (xyz,)
    assert symbols(('xyz,',)) == ((xyz,),)
    assert symbols(('x,y,z,',)) == ((x, y, z),)
    assert symbols(('xyz', 'abc')) == (xyz, abc)
    assert symbols(('xyz,abc',)) == ((xyz, abc),)
    assert symbols(('xyz,abc', 'x,y,z')) == ((xyz, abc), (x, y, z))

    assert symbols(('x', 'y', 'z')) == (x, y, z)
    assert symbols(['x', 'y', 'z']) == [x, y, z]
    assert symbols({'x', 'y', 'z'}) == {x, y, z}

    pytest.raises(ValueError, lambda: symbols(''))
    pytest.raises(ValueError, lambda: symbols(','))
    pytest.raises(ValueError, lambda: symbols('x,,y,,z'))
    pytest.raises(ValueError, lambda: symbols(('x', '', 'y', '', 'z')))

    a, b = symbols('x,y', extended_real=True)
    assert a.is_extended_real
    assert b.is_extended_real

    x0 = Symbol('x0')
    x1 = Symbol('x1')
    x2 = Symbol('x2')

    y0 = Symbol('y0')
    y1 = Symbol('y1')

    assert symbols('x0:0') == ()
    assert symbols('x0:1') == (x0,)
    assert symbols('x0:2') == (x0, x1)
    assert symbols('x0:3') == (x0, x1, x2)

    assert symbols('x:0') == ()
    assert symbols('x:1') == (x0,)
    assert symbols('x:2') == (x0, x1)
    assert symbols('x:3') == (x0, x1, x2)

    assert symbols('x1:1') == ()
    assert symbols('x1:2') == (x1,)
    assert symbols('x1:3') == (x1, x2)

    assert symbols('x1:3,x,y,z') == (x1, x2, x, y, z)

    assert symbols('x:3,y:2') == (x0, x1, x2, y0, y1)
    assert symbols(('x:3', 'y:2')) == ((x0, x1, x2), (y0, y1))

    a = Symbol('a')
    b = Symbol('b')
    c = Symbol('c')
    d = Symbol('d')

    assert symbols('x:z') == (x, y, z)
    assert symbols('a:d,x:z') == (a, b, c, d, x, y, z)
    assert symbols(('a:d', 'x:z')) == ((a, b, c, d), (x, y, z))

    aa = Symbol('aa')
    ab = Symbol('ab')
    ac = Symbol('ac')
    ad = Symbol('ad')

    assert symbols('aa:d') == (aa, ab, ac, ad)
    assert symbols('aa:d,x:z') == (aa, ab, ac, ad, x, y, z)
    assert symbols(('aa:d', 'x:z')) == ((aa, ab, ac, ad), (x, y, z))

    # issue sympy/sympy#6675
    def sym(s):
        return sstr(symbols(s))
    assert sym('a0:4') == '(a0, a1, a2, a3)'
    assert sym('a2:4,b1:3') == '(a2, a3, b1, b2)'
    assert sym('a1(2:4)') == '(a12, a13)'
    assert sym('a0:2.0:2') == '(a0.0, a0.1, a1.0, a1.1)'
    assert sym('aa:cz') == '(aaz, abz, acz)'
    assert sym('aa:c0:2') == '(aa0, aa1, ab0, ab1, ac0, ac1)'
    assert sym('aa:ba:b') == '(aaa, aab, aba, abb)'
    assert sym('a:3b') == '(a0b, a1b, a2b)'
    assert sym('a-1:3b') == '(a-1b, a-2b)'

    c = chr(0)

    assert sym(r'a:2\,:2' + c) == f'(a0,0{c}, a0,1{c}, a1,0{c}, a1,1{c})'
    assert sym('x(:a:3)') == '(x(a0), x(a1), x(a2))'
    assert sym('x(:c):1') == '(xa0, xb0, xc0)'
    assert sym('x((:a)):3') == '(x(a)0, x(a)1, x(a)2)'
    assert sym('x(:a:3') == '(x(a0, x(a1, x(a2)'
    assert sym(':2') == '(0, 1)'
    assert sym(':b') == '(a, b)'
    assert sym(':b:2') == '(a0, a1, b0, b1)'
    assert sym(':2:2') == '(00, 01, 10, 11)'
    assert sym(':b:b') == '(aa, ab, ba, bb)'

    pytest.raises(ValueError, lambda: symbols(':'))
    pytest.raises(ValueError, lambda: symbols('a:'))
    pytest.raises(ValueError, lambda: symbols('::'))
    pytest.raises(ValueError, lambda: symbols('a::'))
    pytest.raises(ValueError, lambda: symbols(':a:'))
    pytest.raises(ValueError, lambda: symbols('::a'))


def test_sympyissue_9057():
    from diofant import beta

    beta(2, 3)  # not raises

    beta = Symbol('beta')
    pytest.raises(TypeError, lambda: beta(2))
    pytest.raises(TypeError, lambda: beta(2.5))
    pytest.raises(TypeError, lambda: beta(2, 3))
