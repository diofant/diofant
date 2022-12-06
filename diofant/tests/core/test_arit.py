import itertools

import pytest

from diofant import (Add, Dummy, E, Float, I, Integer, Mod, Mul, O, Pow,
                     Rational, Symbol, cbrt, comp, cos, exp, factorial, im,
                     log, nan, oo, pi, polar_lift, re, root, sign, sin, sqrt,
                     symbols, sympify, tan, zoo)
from diofant.abc import a, c, x, y, z
from diofant.utilities.randtest import verify_numerically


__all__ = ()


b = Symbol('b', positive=True)


def same_and_same_prec(a, b):
    # stricter matching for Floats
    return a == b and a._prec == b._prec


def test_bug1():
    assert re(x) != x
    x.series(x, 0, 1)
    assert re(x) != x


def test_Symbol():
    e = a*b
    assert e == a*b
    assert a*b*b == a*b**2
    assert a*b*b + c == c + a*b**2
    assert a*b*b - c == -c + a*b**2

    x = Symbol('x', complex=True, extended_real=False)
    assert x.is_imaginary is None  # could be I or 1 + I
    x = Symbol('x', complex=True, imaginary=False)
    assert x.is_extended_real is None  # could be 1 or 1 + I
    x = Symbol('x', real=True)
    assert x.is_complex
    x = Symbol('x', imaginary=True)
    assert x.is_complex
    x = Symbol('x', extended_real=False, imaginary=False)
    assert x.is_complex is None  # might be a non-number


def test_arit0():
    p = 5
    e = a*b
    assert e == a*b
    e = a*b + b*a
    assert e == 2*a*b
    e = a*b + b*a + a*b + p*b*a
    assert e == 8*a*b
    e = a*b + b*a + a*b + p*b*a + a
    assert e == a + 8*a*b
    e = a + a
    assert e == 2*a
    e = a + b + a
    assert e == b + 2*a
    e = a + b*b + a + b*b
    assert e == 2*a + 2*b**2
    e = a + 2 + b*b + a + b*b + p
    assert e == 7 + 2*a + 2*b**2
    e = (a + b*b + a + b*b)*p
    assert e == 5*(2*a + 2*b**2)
    e = (a*b*c + c*b*a + b*a*c)*p
    assert e == 15*a*b*c
    e = (a*b*c + c*b*a + b*a*c)*p - 15*a*b*c
    assert e == 0
    e = 50*(a - a)
    assert e == 0
    e = b*a - b - a*b + b
    assert e == 0
    e = a*b + c**p
    assert e == a*b + c**5
    e = a/b
    assert e == a*b**(-1)
    e = a*2*2
    assert e == 4*a
    e = 2 + a*2/2
    assert e == 2 + a
    e = 2 - a - 2
    assert e == -a
    e = 2*a*2
    assert e == 4*a
    e = 2/a/2
    assert e == a**(-1)
    e = 2**a**2
    assert e == 2**(a**2)
    e = -(1 + a)
    assert e == -1 - a
    e = Rational(1, 2)*(1 + a)
    assert e == Rational(1, 2) + a/2


def test_div():
    e = a/b
    assert e == a*b**(-1)
    e = a/b + c/2
    assert e == a*b**(-1) + c/2
    e = (1 - b)/(b - 1)
    assert e == (1 + -b)*((-1) + b)**(-1)


def test_pow():
    n1 = Integer(1)
    n2 = Integer(2)
    n5 = Integer(5)
    e = a*a
    assert e == a**2
    e = a*a*a
    assert e == a**3
    e = a*a*a*a**6
    assert e == a**9
    e = a*a*a*a**6 - a**9
    assert e == 0
    e = a**(b - b)
    assert e == 1
    e = (a - a)**b
    assert e == 0
    e = (a + 1 - a)**b
    assert e == 1

    e = (a + b + c)**n2
    assert e == (a + b + c)**2
    assert e.expand() == 2*b*c + 2*a*c + 2*a*b + a**2 + c**2 + b**2

    e = (a + b)**n2
    assert e == (a + b)**2
    assert e.expand() == 2*a*b + a**2 + b**2

    e = (a + b)**(n1/n2)
    assert e == sqrt(a + b)
    assert e.expand() == sqrt(a + b)

    n = n5**(n1/n2)
    assert n == sqrt(5)
    e = n*a*b - n*b*a
    assert e == 0
    e = n*a*b + n*b*a
    assert e == 2*a*b*sqrt(5)
    assert e.diff(a) == 2*b*sqrt(5)
    assert e.diff(a) == 2*b*sqrt(5)
    e = a/b**2
    assert e == a*b**(-2)

    assert sqrt(2*(1 + sqrt(2))) == (2*(1 + 2**Rational(1, 2)))**Rational(1, 2)

    assert ((x*y)**3).expand() == y**3 * x**3
    assert ((x*y)**-3).expand() == y**-3 * x**-3

    assert (x**5*(3*x)**3).expand() == 27 * x**8
    assert (x**5*(-3*x)**3).expand() == -27 * x**8
    assert (x**5*(3*x)**(-3)).expand() == Rational(1, 27) * x**2
    assert (x**5*(-3*x)**(-3)).expand() == -Rational(1, 27) * x**2

    # expand_power_exp
    assert (x**(y**(x + exp(x + y)) + z)).expand(deep=False) == \
        x**z*x**(y**(x + exp(x + y)))
    assert (x**(y**(x + exp(x + y)) + z)).expand() == \
        x**z*x**(y**x*y**(exp(x)*exp(y)))

    n = Symbol('n', even=False)
    k = Symbol('k', even=True)
    o = Symbol('o', odd=True)

    assert (-1)**x == (-1)**x
    assert (-1)**n == (-1)**n
    assert (-2)**k == 2**k
    assert (-2*x)**k == (2*x)**k  # we choose not to auto expand this
    assert (-2*x)**o == -(2*x)**o  # but we do handle coefficient sign
    assert (-1)**k == 1


def test_pow2():
    # x**(2*y) is always (x**y)**2 but is only (x**2)**y if
    #                                  x.is_positive or y.is_integer
    # let x = 1 to see why the following are not true.
    assert (-x)**Rational(2, 3) != x**Rational(2, 3)
    assert (-x)**Rational(5, 7) != -x**Rational(5, 7)
    assert cbrt((-x)**2) != (cbrt(-x))**2
    assert sqrt(x**2) != x


def test_pow3():
    assert sqrt(2)**3 == 2 * sqrt(2)
    assert sqrt(2)**3 == sqrt(8)


def test_pow_E():
    assert 2**(y/log(2)) == E**y
    assert 2**(y/log(2)/3) == E**(y/3)
    assert 3**(1/log(-3)) != E
    assert (3 + 2*I)**(1/(log(-3 - 2*I) + I*pi)) == E
    assert (4 + 2*I)**(1/(log(-4 - 2*I) + I*pi)) == E
    assert (3 + 2*I)**(1/(log(-3 - 2*I, 3)/2 + I*pi/log(3)/2)) == 9
    assert (3 + 2*I)**(1/(log(3 + 2*I, 3)/2)) == 9
    # every time tests are run they will affirm with a different random
    # value that this identity holds
    while 1:
        b = x._random()
        _, i = b.as_real_imag()
        if i:
            break
    assert verify_numerically(b**(1/(log(-b) + sign(i)*I*pi).evalf(strict=False)), E)


def test_sympyissue_3516():
    assert root(4, 4) == sqrt(2)


def test_pow_im():
    for m in (-2, -1, 2):
        for d in (3, 4, 5):
            b = m*I
            for i in range(1, 4*d + 1):
                e = Rational(i, d)
                assert (b**e - (b**e).evalf()).evalf(2, chop=True, strict=False) == 0

    e = Rational(7, 3)
    assert (2*x*I)**e == 4*cbrt(2)*(I*x)**e  # same as Wolfram Alpha
    im = symbols('im', imaginary=True)
    assert (2*im*I)**e == 4*cbrt(2)*(I*im)**e

    args = [I, I, I, I, 2]
    e = Rational(1, 3)
    ans = 2**e
    assert Mul(*args, evaluate=False)**e == ans
    assert Mul(*args)**e == ans
    args = [I, I, I, 2]
    e = Rational(1, 3)
    ans = 2**e*(-I)**e
    assert Mul(*args, evaluate=False)**e == ans
    assert Mul(*args)**e == ans
    args.append(-3)
    ans = (6*I)**e
    assert Mul(*args, evaluate=False)**e == ans
    assert Mul(*args)**e == ans
    args.append(-1)
    ans = (-6*I)**e
    assert Mul(*args, evaluate=False)**e == ans
    assert Mul(*args)**e == ans

    args = [I, I, 2]
    e = Rational(1, 3)
    ans = (-2)**e
    assert Mul(*args, evaluate=False)**e == ans
    assert Mul(*args)**e == ans
    args.append(-3)
    ans = 6**e
    assert Mul(*args, evaluate=False)**e == ans
    assert Mul(*args)**e == ans
    args.append(-1)
    ans = (-6)**e
    assert Mul(*args, evaluate=False)**e == ans
    assert Mul(*args)**e == ans
    assert Mul(Pow(-1, Rational(3, 2), evaluate=False), I, I) == I
    assert Mul(I*Pow(I, Rational(1, 2), evaluate=False)) == (-1)**Rational(3, 4)


def test_real_mul():
    assert Float(0) * pi * x == Float(0)
    assert set((Float(1) * pi * x).args) == {Float(1), pi, x}


def test_ncmul():
    A = Symbol('A', commutative=False)
    B = Symbol('B', commutative=False)
    C = Symbol('C', commutative=False)
    assert A*B != B*A
    assert A*B*C != C*B*A
    assert A*b*B*3*C == 3*b*A*B*C
    assert A*b*B*3*C != 3*b*B*A*C
    assert A*b*B*3*C == 3*A*B*C*b

    assert A + B == B + A
    assert (A + B)*C != C*(A + B)

    assert C*(A + B)*C != C*C*(A + B)

    assert A*A == A**2
    assert (A + B)*(A + B) == (A + B)**2

    assert A**-1 * A == 1
    assert A/A == 1
    assert A/(A**2) == 1/A

    assert A/(1 + A) == A/(1 + A)

    assert set((A + B + 2*(A + B)).args) == \
        {A, B, 2*(A + B)}


def test_ncpow():
    x = Symbol('x', commutative=False)
    y = Symbol('y', commutative=False)
    z = Symbol('z', commutative=False)
    b = Symbol('b')

    assert (x**2)*(y**2) != (y**2)*(x**2)
    assert (x**-2)*y != y*(x**2)
    assert 2**x*2**y != 2**(x + y)
    assert 2**x*2**y*2**z != 2**(x + y + z)
    assert 2**x*2**(2*x) == 2**(3*x)
    assert 2**x*2**(2*x)*2**x == 2**(4*x)
    assert exp(x)*exp(y) != exp(y)*exp(x)
    assert exp(x)*exp(y)*exp(z) != exp(y)*exp(x)*exp(z)
    assert exp(x)*exp(y)*exp(z) != exp(x + y + z)
    assert x**a*x**b != x**(a + b)
    assert x**a*x**b*x**c != x**(a + b + c)
    assert x**3*x**4 == x**7
    assert x**3*x**4*x**2 == x**9
    assert x**a*x**(4*a) == x**(5*a)
    assert x**a*x**(4*a)*x**a == x**(6*a)


def test_powerbug():
    assert x**1 != (-x)**1
    assert x**2 == (-x)**2
    assert x**3 != (-x)**3
    assert x**4 == (-x)**4
    assert x**5 != (-x)**5
    assert x**6 == (-x)**6

    assert x**128 == (-x)**128
    assert x**129 != (-x)**129

    assert (2*x)**2 == (-2*x)**2


def test_Mul_doesnt_expand_exp():
    assert exp(x)*exp(y) == exp(x)*exp(y)
    assert 2**x*2**y == 2**x*2**y
    assert x**2*x**3 == x**5
    assert 2**x*3**x == 6**x
    assert x**y*x**(2*y) == x**(3*y)
    assert sqrt(2)*sqrt(2) == 2
    assert 2**x*2**(2*x) == 2**(3*x)
    assert sqrt(2)*root(2, 4)*5**Rational(3, 4) == 10**Rational(3, 4)
    assert (x**(-log(5)/log(3))*x)/(x*x**(-log(5)/log(3))) == sympify(1)


def test_Add_Mul_is_integer():
    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True)

    assert (2*k).is_integer is True
    assert (-k).is_integer is True
    assert (k/3).is_integer is None
    assert (x*k*n).is_integer is None

    assert (k + n).is_integer is True
    assert (k + x).is_integer is None
    assert (k + n*x).is_integer is None
    assert (k + n/3).is_integer is None

    assert ((1 + sqrt(3))*(-sqrt(3) + 1)).is_integer is not False
    assert (1 + (1 + sqrt(3))*(-sqrt(3) + 1)).is_integer is not False


def test_Add_Mul_is_finite():
    x = Symbol('x', extended_real=True, finite=False)
    y = Symbol('y', real=True)
    z = Symbol('z', real=True)

    assert y.is_finite is True
    assert (x*y).is_finite is False
    assert (1024*y).is_finite is True
    assert (y*exp(x)).is_finite is not True
    assert (y*z).is_finite is True
    assert (x*y*exp(x)).is_finite is not True

    assert (y - 67).is_finite is True
    assert (y + exp(x)).is_finite is not True
    assert (1 + x).is_finite is False
    assert (1 + x**2 + (1 + x)*(1 - x)).is_finite is None
    assert (sqrt(2)*(1 + x)).is_finite is False
    assert (sqrt(2)*(1 + x)*(1 - x)).is_finite is False


def test_Mul_is_even_odd():
    x = Symbol('x', integer=True)
    y = Symbol('y', integer=True)

    k = Symbol('k', odd=True)
    n = Symbol('n', odd=True)
    m = Symbol('m', even=True)

    assert (2*x).is_even is True
    assert (2*x).is_odd is False

    assert (3*x).is_even is None
    assert (3*x).is_odd is None

    assert (k/3).is_integer is None
    assert (k/3).is_even is None
    assert (k/3).is_odd is None

    assert (2*n).is_even is True
    assert (2*n).is_odd is False

    assert (2*m).is_even is True
    assert (2*m).is_odd is False

    assert (-n).is_even is False
    assert (-n).is_odd is True

    assert (k*n).is_even is False
    assert (k*n).is_odd is True

    assert (k*m).is_even is True
    assert (k*m).is_odd is False

    assert (k*n*m).is_even is True
    assert (k*n*m).is_odd is False

    assert (k*m*x).is_even is True
    assert (k*m*x).is_odd is False

    # issue sympy/sympy#6791:
    assert (x/2).is_integer is None
    assert (k/2).is_integer is False
    assert (m/2).is_integer is True

    assert (x*y).is_even is None
    assert (x*x).is_even is None
    assert (x*(x + k)).is_even is True
    assert (x*(x + m)).is_even is None

    assert (x*y).is_odd is None
    assert (x*x).is_odd is None
    assert (x*(x + k)).is_odd is False
    assert (x*(x + m)).is_odd is None

    assert (x*y*(y + m)).is_even is None
    assert (x*y*(y + m)).is_odd is None


def test_even_odd_in_ternary_integer_product():
    # Tests that oddness inference is independent of term ordering.
    # We try to force a different order by modifying symbol names.

    # issues sympy/sympy#9127 and diofant/diofant#1003

    x = Symbol('x', integer=True)
    y = Symbol('y', integer=True)
    k = Symbol('k', odd=True)

    assert (x*y*(y + k)).is_even is True
    assert (y*x*(x + k)).is_even is True

    assert (x*y*(y + k)).is_odd is False
    assert (y*x*(x + k)).is_odd is False


def test_Mul_is_rational():
    n = Symbol('n', integer=True)
    m = Symbol('m', integer=True, nonzero=True)

    assert (n/m).is_rational is True
    assert (x/pi).is_rational is None
    assert (x/n).is_rational is None
    assert (m/pi).is_rational is False

    r = Symbol('r', rational=True)
    assert (pi*r).is_rational is None

    z = Symbol('z', zero=True)
    i = Symbol('i', imaginary=True)
    assert (z*i).is_rational


def test_Add_is_rational():
    n = Symbol('n', rational=True)
    m = Symbol('m', rational=True)

    assert (n + m).is_rational is True
    assert (x + pi).is_rational is None
    assert (x + n).is_rational is None
    assert (n + pi).is_rational is False


def test_Add_is_even_odd():
    x = Symbol('x', integer=True)

    k = Symbol('k', odd=True)
    n = Symbol('n', odd=True)
    m = Symbol('m', even=True)

    assert (k + 7).is_even is True
    assert (k + 7).is_odd is False

    assert (-k + 7).is_even is True
    assert (-k + 7).is_odd is False

    assert (k - 12).is_even is False
    assert (k - 12).is_odd is True

    assert (-k - 12).is_even is False
    assert (-k - 12).is_odd is True

    assert (k + n).is_even is True
    assert (k + n).is_odd is False

    assert (k + m).is_even is False
    assert (k + m).is_odd is True

    assert (k + n + m).is_even is True
    assert (k + n + m).is_odd is False

    assert (k + n + x + m).is_even is None
    assert (k + n + x + m).is_odd is None


def test_Mul_is_negative_positive():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=False, complex=True)
    z = Symbol('z', zero=True)

    e = 2*z
    assert e.is_Mul
    assert e.is_positive is False
    assert e.is_negative is False
    assert (x*y).is_positive is None

    neg = Symbol('neg', negative=True)
    pos = Symbol('pos', positive=True)
    nneg = Symbol('nneg', nonnegative=True)
    npos = Symbol('npos', nonpositive=True)

    assert neg.is_negative is True
    assert (-neg).is_negative is False
    assert (2*neg).is_negative is True

    assert (2*pos)._eval_is_negative() is False
    assert (2*pos).is_negative is False

    assert pos.is_negative is False
    assert (-pos).is_negative is True
    assert (2*pos).is_negative is False

    assert (pos*neg).is_negative is True
    assert (2*pos*neg).is_negative is True
    assert (-pos*neg).is_negative is False
    assert (pos*neg*y).is_negative is False     # y.is_extended_real=F;  !real -> !neg

    assert nneg.is_negative is False
    assert (-nneg).is_negative is None
    assert (2*nneg).is_negative is False

    assert npos.is_negative is None
    assert (-npos).is_negative is False
    assert (2*npos).is_negative is None

    assert (nneg*npos).is_negative is None

    assert (neg*nneg).is_negative is None
    assert (neg*npos).is_negative is False

    assert (pos*nneg).is_negative is False
    assert (pos*npos).is_negative is None

    assert (npos*neg*nneg).is_negative is False
    assert (npos*pos*nneg).is_negative is None

    assert (-npos*neg*nneg).is_negative is None
    assert (-npos*pos*nneg).is_negative is False

    assert (17*npos*neg*nneg).is_negative is False
    assert (17*npos*pos*nneg).is_negative is None

    assert (neg*npos*pos*nneg).is_negative is False

    assert (x*neg).is_negative is None
    assert (nneg*npos*pos*x*neg).is_negative is None

    assert neg.is_positive is False
    assert (-neg).is_positive is True
    assert (2*neg).is_positive is False

    assert pos.is_positive is True
    assert (-pos).is_positive is False
    assert (2*pos).is_positive is True

    assert (pos*neg).is_positive is False
    assert (2*pos*neg).is_positive is False
    assert (-pos*neg).is_positive is True
    assert (-pos*neg*y).is_positive is False    # y.is_extended_real=F;  !real -> !neg

    assert nneg.is_positive is None
    assert (-nneg).is_positive is False
    assert (2*nneg).is_positive is None

    assert npos.is_positive is False
    assert (-npos).is_positive is None
    assert (2*npos).is_positive is False

    assert (nneg*npos).is_positive is False

    assert (neg*nneg).is_positive is False
    assert (neg*npos).is_positive is None

    assert (pos*nneg).is_positive is None
    assert (pos*npos).is_positive is False

    assert (npos*neg*nneg).is_positive is None
    assert (npos*pos*nneg).is_positive is False

    assert (-npos*neg*nneg).is_positive is False
    assert (-npos*pos*nneg).is_positive is None

    assert (17*npos*neg*nneg).is_positive is None
    assert (17*npos*pos*nneg).is_positive is False

    assert (neg*npos*pos*nneg).is_positive is None

    assert (x*neg).is_positive is None
    assert (nneg*npos*pos*x*neg).is_positive is None


def test_Mul_is_negative_positive_2():
    a = Symbol('a', nonnegative=True)
    af = Symbol('af', nonnegative=True, finite=True)
    b = Symbol('b', nonnegative=True)
    bf = Symbol('bf', nonnegative=True, finite=True)
    c = Symbol('c', nonpositive=True)
    cf = Symbol('cf', nonpositive=True, finite=True)
    d = Symbol('d', nonpositive=True)
    df = Symbol('df', nonpositive=True, finite=True)

    assert (af*bf).is_nonnegative is True
    assert (a*b).is_negative is False
    assert (a*b).is_zero is None
    assert (a*b).is_positive is None

    assert (cf*df).is_nonnegative is True
    assert (c*d).is_negative is False
    assert (c*d).is_zero is None
    assert (c*d).is_positive is None

    assert (af*cf).is_nonpositive is True
    assert (a*c).is_positive is False
    assert (a*c).is_zero is None
    assert (a*c).is_negative is None


def test_Mul_is_nonpositive_nonnegative():
    x = Symbol('x', extended_real=True)

    k = Symbol('k', negative=True)
    kf = Symbol('kf', negative=True, finite=True)
    n = Symbol('n', positive=True)
    nf = Symbol('nf', positive=True, finite=True)
    u = Symbol('u', nonnegative=True)
    v = Symbol('v', nonpositive=True)
    uf = Symbol('uf', nonnegative=True, finite=True)
    vf = Symbol('vf', nonpositive=True, finite=True)

    assert k.is_nonpositive is True
    assert (-k).is_nonpositive is False
    assert (2*k).is_nonpositive is True

    assert n.is_nonpositive is False
    assert (-n).is_nonpositive is True
    assert (2*n).is_nonpositive is False

    assert (n*k).is_nonpositive is True
    assert (2*n*k).is_nonpositive is True
    assert (-n*k).is_nonpositive is False

    assert u.is_nonpositive is None
    assert (-u).is_nonpositive is True
    assert (2*u).is_nonpositive is None

    assert v.is_nonpositive is True
    assert (-v).is_nonpositive is None
    assert (2*v).is_nonpositive is True

    assert (uf*vf).is_nonpositive is True

    assert (k*u).is_nonpositive is True
    assert (k*v).is_nonpositive is None

    assert (n*u).is_nonpositive is None
    assert (n*v).is_nonpositive is True

    assert (v*k*u).is_nonpositive is None
    assert (vf*nf*uf).is_nonpositive is True

    assert (-vf*kf*uf).is_nonpositive is True
    assert (-v*n*u).is_nonpositive is None

    assert (17*v*k*u).is_nonpositive is None
    assert (17*vf*nf*uf).is_nonpositive is True

    assert (k*v*n*u).is_nonpositive is None

    assert (x*k).is_nonpositive is None
    assert (u*v*n*x*k).is_nonpositive is None

    assert k.is_nonnegative is False
    assert (-k).is_nonnegative is True
    assert (2*k).is_nonnegative is False

    assert n.is_nonnegative is True
    assert (-n).is_nonnegative is False
    assert (2*n).is_nonnegative is True

    assert (n*k).is_nonnegative is False
    assert (2*n*k).is_nonnegative is False
    assert (-n*k).is_nonnegative is True

    assert u.is_nonnegative is True
    assert (-u).is_nonnegative is None
    assert (2*u).is_nonnegative is True

    assert v.is_nonnegative is None
    assert (-v).is_nonnegative is True
    assert (2*v).is_nonnegative is None

    assert (u*v).is_nonnegative is None

    assert (k*u).is_nonnegative is None
    assert (k*v).is_nonnegative is True

    assert (n*u).is_nonnegative is True
    assert (n*v).is_nonnegative is None

    assert (vf*kf*uf).is_nonnegative is True
    assert (v*n*u).is_nonnegative is None

    assert (-v*k*u).is_nonnegative is None
    assert (-vf*nf*uf).is_nonnegative is True

    assert (17*vf*kf*uf).is_nonnegative is True
    assert (17*v*nf*u).is_nonnegative is None

    assert (kf*vf*nf*uf).is_nonnegative is True

    assert (x*k).is_nonnegative is None
    assert (u*v*n*x*k).is_nonnegative is None


def test_Add_is_negative_positive():
    x = Symbol('x', extended_real=True)

    k = Symbol('k', negative=True)
    n = Symbol('n', positive=True)
    u = Symbol('u', nonnegative=True)
    v = Symbol('v', nonpositive=True)

    assert (k - 2).is_negative is True
    assert (k + 17).is_negative is None
    assert (-k - 5).is_negative is None
    assert (-k + 123).is_negative is False

    assert (k - n).is_negative is True
    assert (k + n).is_negative is None
    assert (-k - n).is_negative is None
    assert (-k + n).is_negative is False

    assert (k - n - 2).is_negative is True
    assert (k + n + 17).is_negative is None
    assert (-k - n - 5).is_negative is None
    assert (-k + n + 123).is_negative is False

    assert (-2*k + 123*n + 17).is_negative is False

    assert (k + u).is_negative is None
    assert (k + v).is_negative is True
    assert (n + u).is_negative is False
    assert (n + v).is_negative is None

    assert (u - v).is_negative is False
    assert (u + v).is_negative is None
    assert (-u - v).is_negative is None
    assert (-u + v).is_negative is None

    assert (u - v + n + 2).is_negative is False
    assert (u + v + n + 2).is_negative is None
    assert (-u - v + n + 2).is_negative is None
    assert (-u + v + n + 2).is_negative is None

    assert (k + x).is_negative is None
    assert (k + x - n).is_negative is None

    assert (k - 2).is_positive is False
    assert (k + 17).is_positive is None
    assert (-k - 5).is_positive is None
    assert (-k + 123).is_positive is True

    assert (k - n).is_positive is False
    assert (k + n).is_positive is None
    assert (-k - n).is_positive is None
    assert (-k + n).is_positive is True

    assert (k - n - 2).is_positive is False
    assert (k + n + 17).is_positive is None
    assert (-k - n - 5).is_positive is None
    assert (-k + n + 123).is_positive is True

    assert (-2*k + 123*n + 17).is_positive is True

    assert (k + u).is_positive is None
    assert (k + v).is_positive is False
    assert (n + u).is_positive is True
    assert (n + v).is_positive is None

    assert (u - v).is_positive is None
    assert (u + v).is_positive is None
    assert (-u - v).is_positive is None
    assert (-u + v).is_positive is False

    assert (u - v - n - 2).is_positive is None
    assert (u + v - n - 2).is_positive is None
    assert (-u - v - n - 2).is_positive is None
    assert (-u + v - n - 2).is_positive is False

    assert (n + x).is_positive is None
    assert (n + x - k).is_positive is None

    z = (-3 - sqrt(5) + (-sqrt(10)/2 - sqrt(2)/2)**2)
    assert z.is_zero
    z = sqrt(1 + sqrt(3)) + sqrt(3 + 3*sqrt(3)) - sqrt(10 + 6*sqrt(3))
    assert z.is_zero


def test_Add_is_nonpositive_nonnegative():
    x = Symbol('x', extended_real=True)

    k = Symbol('k', negative=True)
    n = Symbol('n', positive=True)
    u = Symbol('u', nonnegative=True)
    v = Symbol('v', nonpositive=True)

    assert (u - 2).is_nonpositive is None
    assert (u + 17).is_nonpositive is False
    assert (-u - 5).is_nonpositive is True
    assert (-u + 123).is_nonpositive is None

    assert (u - v).is_nonpositive is None
    assert (u + v).is_nonpositive is None
    assert (-u - v).is_nonpositive is None
    assert (-u + v).is_nonpositive is True

    assert (u - v - 2).is_nonpositive is None
    assert (u + v + 17).is_nonpositive is None
    assert (-u - v - 5).is_nonpositive is None
    assert (-u + v - 123).is_nonpositive is True

    assert (-2*u + 123*v - 17).is_nonpositive is True

    assert (k + u).is_nonpositive is None
    assert (k + v).is_nonpositive is True
    assert (n + u).is_nonpositive is False
    assert (n + v).is_nonpositive is None

    assert (k - n).is_nonpositive is True
    assert (k + n).is_nonpositive is None
    assert (-k - n).is_nonpositive is None
    assert (-k + n).is_nonpositive is False

    assert (k - n + u + 2).is_nonpositive is None
    assert (k + n + u + 2).is_nonpositive is None
    assert (-k - n + u + 2).is_nonpositive is None
    assert (-k + n + u + 2).is_nonpositive is False

    assert (u + x).is_nonpositive is None
    assert (v - x - n).is_nonpositive is None

    assert (u - 2).is_nonnegative is None
    assert (u + 17).is_nonnegative is True
    assert (-u - 5).is_nonnegative is False
    assert (-u + 123).is_nonnegative is None

    assert (u - v).is_nonnegative is True
    assert (u + v).is_nonnegative is None
    assert (-u - v).is_nonnegative is None
    assert (-u + v).is_nonnegative is None

    assert (u - v + 2).is_nonnegative is True
    assert (u + v + 17).is_nonnegative is None
    assert (-u - v - 5).is_nonnegative is None
    assert (-u + v - 123).is_nonnegative is False

    assert (2*u - 123*v + 17).is_nonnegative is True

    assert (k + u).is_nonnegative is None
    assert (k + v).is_nonnegative is False
    assert (n + u).is_nonnegative is True
    assert (n + v).is_nonnegative is None

    assert (k - n).is_nonnegative is False
    assert (k + n).is_nonnegative is None
    assert (-k - n).is_nonnegative is None
    assert (-k + n).is_nonnegative is True

    assert (k - n - u - 2).is_nonnegative is False
    assert (k + n - u - 2).is_nonnegative is None
    assert (-k - n - u - 2).is_nonnegative is None
    assert (-k + n - u - 2).is_nonnegative is None

    assert (u - x).is_nonnegative is None
    assert (v + x + n).is_nonnegative is None


def test_Pow_is_integer():
    x = Symbol('x')

    k = Symbol('k', integer=True)
    n = Symbol('n', integer=True, nonnegative=True)
    m = Symbol('m', integer=True, positive=True)

    assert (k**2).is_integer is True
    assert (k**(-2)).is_integer is None
    assert ((m + 1)**(-2)).is_integer is False
    assert (m**(-1)).is_integer is None  # issue sympy/sympy#8580

    assert (2**k).is_integer is None
    assert (2**(-k)).is_integer is None

    assert (2**n).is_integer is True
    assert (2**(-n)).is_integer is None

    assert (2**m).is_integer is True
    assert (2**(-m)).is_integer is False

    assert (x**2).is_integer is None
    assert (2**x).is_integer is None

    assert (k**n).is_integer is True
    assert (k**(-n)).is_integer is None

    assert (k**x).is_integer is None
    assert (x**k).is_integer is None

    assert (k**(n*m)).is_integer is True
    assert (k**(-n*m)).is_integer is None

    assert sqrt(3).is_integer is False
    assert sqrt(.3).is_integer is False
    assert Pow(3, 2, evaluate=False).is_integer is True
    assert Pow(3, 0, evaluate=False).is_integer is True
    assert Pow(3, -2, evaluate=False).is_integer is False
    assert Pow(Rational(1, 2), 3, evaluate=False).is_integer is False
    # decided by re-evaluating
    assert Pow(3, Rational(1, 2), evaluate=False).is_integer is False
    assert Pow(3, Rational(1, 2), evaluate=False).is_integer is False
    assert Pow(4, Rational(1, 2), evaluate=False).is_integer is True
    assert Pow(Rational(1, 2), -2, evaluate=False).is_integer is True

    assert ((-1)**k).is_integer

    # issue sympy/sympy#23287
    assert (x**2/2).is_integer is None

    x = Symbol('x', extended_real=True, integer=False)
    assert (x**2).is_integer is None  # issue sympy/sympy#8641


def test_Pow_is_real():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True, positive=True)

    assert (x**2).is_extended_real is True
    assert (x**3).is_extended_real is True
    assert (x**x).is_extended_real is None
    assert (y**x).is_extended_real is True

    assert cbrt(x).is_extended_real is None
    assert cbrt(y).is_extended_real is True

    assert sqrt(-1 - sqrt(2)).is_extended_real is False

    i = Symbol('i', imaginary=True)
    ni = Symbol('ni', imaginary=True, nonzero=True)
    assert (i**i).is_extended_real is None
    assert (I**i).is_extended_real is True
    assert ((-I)**i).is_extended_real is True
    assert (2**i).is_extended_real is None  # (2**(pi/log(2) * I)) is real, 2**I is not
    assert (2**I).is_extended_real is False
    assert (2**-I).is_extended_real is False
    assert (i**2).is_extended_real is True
    assert (i**3).is_extended_real is False
    assert (i**x).is_extended_real is None  # could be (-I)**(2/3)
    e = Symbol('e', even=True)
    o = Symbol('o', odd=True)
    k = Symbol('k', integer=True)
    assert (i**(e**2)).is_extended_real is True
    assert (i**o).is_extended_real is False
    assert (i**k).is_extended_real is None
    assert (i**(4*k)).is_extended_real is None
    assert (ni**(4*k)).is_extended_real is True
    assert (x**i).is_extended_real is None
    assert (i**(Rational(1, 2) + x)).is_extended_real is None
    assert Pow(I, 2, evaluate=False).is_extended_real

    x = Symbol('x', nonnegative=True)
    y = Symbol('y', nonnegative=True)
    assert im(x**y).expand(complex=True) is Integer(0)
    assert (x**y).is_extended_real is True
    i = Symbol('i', imaginary=True)
    assert (exp(i)**I).is_extended_real is True
    assert log(exp(i)).is_imaginary is None  # i could be 2*pi*I
    c = Symbol('c', complex=True)
    assert log(c).is_extended_real is None  # c could be 0 or 2, too
    assert log(exp(c)).is_extended_real is None  # log(0), log(E), ...
    n = Symbol('n', negative=False)
    assert log(n).is_extended_real is None
    n = Symbol('n', nonnegative=True)
    assert log(n).is_extended_real is None

    assert sqrt(-I).is_extended_real is False  # issue sympy/sympy#7843

    # issue sympy/sympy#6631
    assert ((-1)**I).is_extended_real is True
    assert ((-1)**(I*2)).is_extended_real is True
    assert ((-1)**(I/2)).is_extended_real is True
    assert ((-1)**(I*pi)).is_extended_real is True
    assert (I**(I + 2)).is_extended_real is True


def test_real_Pow():
    k = Symbol('k', integer=True, nonzero=True)
    assert (k**(I*pi/log(k))).is_extended_real


def test_Pow_is_finite():
    x = Symbol('x', extended_real=True)
    p = Symbol('p', positive=True)
    n = Symbol('n', negative=True)
    y = Symbol('y', real=True)

    assert (x**2).is_finite is None  # x could be oo
    assert (x**x).is_finite is None  # ditto
    assert (p**x).is_finite is None  # ditto
    assert (n**x).is_finite is None  # ditto
    assert (1/pi).is_finite
    assert (y**2).is_finite is True
    assert (y**x).is_finite is None
    assert (y**exp(x)).is_finite is None
    assert (1/y).is_finite is None  # if zero, no, otherwise yes
    assert (1/exp(x)).is_finite is None  # x could be -oo


def test_Pow_is_even_odd():
    k = Symbol('k', even=True)
    n = Symbol('n', odd=True)
    m = Symbol('m', integer=True, nonnegative=True)
    p = Symbol('p', integer=True, positive=True)

    assert ((-1)**n).is_odd
    assert ((-1)**k).is_odd
    assert ((-1)**(m - p)).is_odd

    assert (k**2).is_even is True
    assert (n**2).is_even is False
    assert (2**k).is_even is None
    assert (x**2).is_even is None

    assert (k**m).is_even is None
    assert (n**m).is_even is False

    assert (k**p).is_even is True
    assert (n**p).is_even is False

    assert (m**k).is_even is None
    assert (p**k).is_even is None

    assert (m**n).is_even is None
    assert (p**n).is_even is None

    assert (k**x).is_even is None
    assert (n**x).is_even is None

    assert (k**2).is_odd is False
    assert (n**2).is_odd is True
    assert (3**k).is_odd is None

    assert (k**m).is_odd is None
    assert (n**m).is_odd is True

    assert (k**p).is_odd is False
    assert (n**p).is_odd is True

    assert (m**k).is_odd is None
    assert (p**k).is_odd is None

    assert (m**n).is_odd is None
    assert (p**n).is_odd is None

    assert (k**x).is_odd is None
    assert (n**x).is_odd is None


def test_Pow_is_negative_positive():
    er = Symbol('er', extended_real=True)
    r = Symbol('r', real=True)
    p = Symbol('p', positive=True)

    k = Symbol('k', integer=True, positive=True)
    n = Symbol('n', even=True)
    m = Symbol('m', odd=True)

    assert (2**p).is_positive is True
    assert (2**r).is_positive is True
    assert ((-2)**er).is_positive is None
    assert ((-2)**n).is_positive is True
    assert ((-2)**m).is_positive is False

    assert (k**2).is_positive is True
    assert (k**(-2)).is_positive is True

    assert (k**p).is_positive is True
    assert (k**r).is_positive is True
    assert ((-k)**er).is_positive is None
    assert ((-k)**n).is_positive is True
    assert ((-k)**m).is_positive is False

    assert (2**er).is_negative is False
    assert ((-2)**er).is_negative is None
    assert ((-2)**n).is_negative is False
    assert ((-2)**m).is_negative is True

    assert (k**2).is_negative is False
    assert (k**(-2)).is_negative is False

    assert (k**er).is_negative is False
    assert ((-k)**er).is_negative is None
    assert ((-k)**n).is_negative is False
    assert ((-k)**m).is_negative is True

    assert (2**x).is_positive is None
    assert (2**x).is_negative is None

    s = Symbol('s', nonpositive=True)
    assert (s**n).is_negative is False
    assert (s**m).is_positive is None
    n = Symbol('n', even=True, nonnegative=True)
    m = Symbol('m', odd=True, nonnegative=True)
    assert ((-p)**n).is_positive is True
    assert ((-p)**m).is_positive is False
    assert (s**m).is_positive is False
    assert ((-p)**(n + 1)).is_negative is True
    s = Symbol('s', nonpositive=True, finite=True)
    assert ((s - 1)**n).is_positive is True
    assert ((s - 1)**m).is_positive is False
    assert (s**m).is_positive is False
    assert ((s - 1)**m).is_negative is True

    i = Symbol('i', imaginary=True)
    assert (i**4).is_positive is None  # issue diofant/diofant#956


def test_Pow_is_zero():
    z = Symbol('z', zero=True)
    e = z**2
    assert e.is_zero
    assert e.is_positive is False
    assert e.is_negative is False

    assert Pow(0, 0, evaluate=False).is_nonzero
    assert Pow(0, 3, evaluate=False).is_zero
    assert Pow(0, oo, evaluate=False).is_zero
    assert Pow(0, -3, evaluate=False).is_nonzero
    assert Pow(0, -oo, evaluate=False).is_nonzero
    assert Pow(2, 2, evaluate=False).is_nonzero

    a = Symbol('a', zero=False)
    assert Pow(a, 3).is_nonzero is True  # issue sympy/sympy#7965

    assert Pow(2, oo, evaluate=False).is_nonzero
    assert Pow(2, -oo, evaluate=False).is_zero
    assert Pow(Rational(1, 2), oo, evaluate=False).is_zero
    assert Pow(Rational(1, 2), -oo, evaluate=False).is_nonzero

    n = Symbol('n', nonzero=True)
    assert Pow(n, oo).is_zero is None

    e = Symbol('e', nonpositive=True)
    assert Pow(0, e).is_zero is False


def test_Pow_is_nonpositive_nonnegative():
    x = Symbol('x', extended_real=True)
    r = Symbol('r', real=True)
    p = Symbol('p', positive=True)

    k = Symbol('k', integer=True, nonnegative=True)
    l = Symbol('l', integer=True, positive=True)
    n = Symbol('n', even=True)
    m = Symbol('m', odd=True)

    assert (x**(4*k)).is_nonnegative is True
    assert (2**x).is_nonnegative is True
    assert ((-2)**x).is_nonnegative is None
    assert ((-2)**n).is_nonnegative is True
    assert ((-2)**m).is_nonnegative is False

    assert (k**2).is_nonnegative is True
    assert (k**(-2)).is_nonnegative is None
    assert (k**k).is_nonnegative is True

    assert (k**x).is_nonnegative is None    # NOTE (0**x).is_extended_real = U
    assert (l**p).is_nonnegative is True
    assert (l**p).is_positive is True
    assert (l**r).is_nonnegative is True
    assert (l**r).is_positive is True
    assert ((-k)**x).is_nonnegative is None
    assert ((-k)**n).is_nonnegative is None
    assert ((-k)**m).is_nonnegative is None

    assert (2**p).is_nonpositive is False
    assert (2**r).is_nonpositive is False
    assert ((-2)**x).is_nonpositive is None
    assert ((-2)**n).is_nonpositive is False
    assert ((-2)**m).is_nonpositive is True

    assert (k**2).is_nonpositive is None
    assert (k**(-2)).is_nonpositive is None

    assert (k**x).is_nonpositive is None
    assert ((-k)**x).is_nonpositive is None
    assert ((-k)**n).is_nonpositive is None
    assert ((-k)**m).is_nonpositive is None

    assert (x**2).is_nonnegative is True
    i = symbols('i', imaginary=True)
    ni = symbols('ni', imaginary=True, nonzero=True)
    assert (ni**2).is_nonpositive is True
    assert (ni**4).is_nonpositive is False
    assert (i**3).is_nonpositive is False
    assert (I**i).is_nonnegative is True
    assert (exp(I)**i).is_nonnegative is True


def test_Mul_is_imaginary_real():
    r = Symbol('r', real=True)
    p = Symbol('p', positive=True, real=True)
    i = Symbol('i', imaginary=True)
    ii = Symbol('ii', imaginary=True)
    ni = Symbol('ni', imaginary=True, nonzero=True)
    nii = Symbol('nii', imaginary=True, nonzero=True)

    assert I.is_imaginary is True
    assert I.is_extended_real is False
    assert (-I).is_imaginary is True
    assert (-I).is_extended_real is False
    assert (3*I).is_imaginary is True
    assert (3*I).is_extended_real is False
    assert (I*I).is_imaginary is False
    assert (I*I).is_extended_real is True

    e = (p + p*I)
    j = Symbol('j', integer=True, zero=False)
    assert (e**j).is_extended_real is None
    assert (e**(2*j)).is_extended_real is None
    assert (e**j).is_imaginary is None
    assert (e**(2*j)).is_imaginary is None

    assert (e**-1).is_imaginary is False
    assert (e**2).is_imaginary
    assert (e**3).is_imaginary is False
    assert (e**4).is_imaginary is False
    assert (e**5).is_imaginary is False
    assert (e**-1).is_extended_real is False
    assert (e**2).is_extended_real is False
    assert (e**3).is_extended_real is False
    assert (e**4).is_extended_real
    assert (e**5).is_extended_real is False
    assert (e**3).is_complex

    assert (r*i).is_imaginary is True
    assert (r*i).is_extended_real is None

    assert (x*i).is_imaginary is None
    assert (x*i).is_extended_real is None

    assert (ni*nii).is_imaginary is False
    assert (i*ii).is_extended_real is True

    assert (r*i*ii).is_imaginary is None
    assert (r*i*ii).is_extended_real is True

    # Github's issue sympy/sympy#5874:
    nr = Symbol('nr', extended_real=False, complex=True)
    a = Symbol('a', extended_real=True, nonzero=True)
    b = Symbol('b', extended_real=True)
    assert (i*nr).is_extended_real is None
    assert (a*nr).is_extended_real is False
    assert (b*nr).is_extended_real is None


def test_Add_is_comparable():
    assert (x + y).is_comparable is False
    assert (x + 1).is_comparable is False
    assert (Rational(1, 3) - sqrt(8)).is_comparable is True


def test_Mul_is_comparable():
    assert (x*y).is_comparable is False
    assert (x*2).is_comparable is False
    assert (sqrt(2)*Rational(1, 3)).is_comparable is True


def test_Pow_is_comparable():
    assert (x**y).is_comparable is False
    assert (x**2).is_comparable is False
    assert (sqrt(Rational(1, 3))).is_comparable is True


def test_Add_is_positive_2():
    e = Rational(1, 3) - sqrt(8)
    assert e.is_positive is False
    assert e.is_negative is True

    e = pi - 1
    assert e.is_positive is True
    assert e.is_negative is False


def test_Add_is_irrational():
    i = Symbol('i', irrational=True)

    assert i.is_irrational is True
    assert i.is_rational is False

    assert (i + 1).is_irrational is True
    assert (i + 1).is_rational is False


@pytest.mark.xfail
def test_sympyissue_3531_fail():
    class MightyNumeric(tuple):
        def __rtruediv__(self, other):
            return 'something'

    assert sympify(1)/MightyNumeric((1, 2)) == 'something'


def test_sympyissue_3531():
    class Foo:
        def __init__(self):
            self.field = 1.0

        def __mul__(self, other):
            self.field = self.field * other

        def __rmul__(self, other):
            self.field = other * self.field
    f = Foo()
    assert f*x == x*f


def test_bug3():
    b = Symbol('b', positive=True)
    e = 2*a + b
    f = b + 2*a
    assert e == f


def test_suppressed_evaluation():
    a = Add(0, 3, 2, evaluate=False)
    b = Mul(1, 3, 2, evaluate=False)
    c = Pow(3, 2, evaluate=False)
    assert a != 6
    assert isinstance(a, Add)
    assert a.args == (0, 3, 2)
    assert b != 6
    assert isinstance(b, Mul)
    assert b.args == (1, 3, 2)
    assert c != 9
    assert isinstance(c, Pow)
    assert c.args == (3, 2)


def test_Add_as_coeff_mul():
    # issue sympy/sympy#5524.  These should all be (1, self)
    assert (x + 1).as_coeff_mul() == (1, (x + 1,))
    assert (x + 2).as_coeff_mul() == (1, (x + 2,))
    assert (x + 3).as_coeff_mul() == (1, (x + 3,))

    assert (x - 1).as_coeff_mul() == (1, (x - 1,))
    assert (x - 2).as_coeff_mul() == (1, (x - 2,))
    assert (x - 3).as_coeff_mul() == (1, (x - 3,))

    n = Symbol('n', integer=True)
    assert (n + 1).as_coeff_mul() == (1, (n + 1,))
    assert (n + 2).as_coeff_mul() == (1, (n + 2,))
    assert (n + 3).as_coeff_mul() == (1, (n + 3,))

    assert (n - 1).as_coeff_mul() == (1, (n - 1,))
    assert (n - 2).as_coeff_mul() == (1, (n - 2,))
    assert (n - 3).as_coeff_mul() == (1, (n - 3,))


def test_Pow_as_coeff_mul_doesnt_expand():
    assert exp(x + y).as_coeff_mul() == (1, (exp(x + y),))
    assert exp(x + exp(x + y)) != exp(x + exp(x)*exp(y))


def test_sympyissue_3514():
    assert sqrt(Rational(1, 2)) * sqrt(6) == 2 * sqrt(3)/2
    assert Rational(1, 2)*sqrt(6)*sqrt(2) == sqrt(3)
    assert sqrt(6)/2*sqrt(2) == sqrt(3)
    assert sqrt(6)*sqrt(2)/2 == sqrt(3)


def test_make_args():
    assert Add.make_args(x) == (x,)
    assert Mul.make_args(x) == (x,)

    assert Add.make_args(x*y*z) == (x*y*z,)
    assert Mul.make_args(x*y*z) == (x*y*z).args

    assert Add.make_args(x + y + z) == (x + y + z).args
    assert Mul.make_args(x + y + z) == (x + y + z,)

    assert Add.make_args((x + y)**z) == ((x + y)**z,)
    assert Mul.make_args((x + y)**z) == ((x + y)**z,)


def test_sympyissue_5126():
    assert (-2)**x*(-3)**x != 6**x
    i = Symbol('i', integer=1)
    assert (-2)**i*(-3)**i == 6**i


def test_Rational_as_content_primitive():
    c, p = Integer(1), Integer(0)
    assert (c*p).as_content_primitive() == (c, p)
    c, p = Rational(1, 2), Integer(1)
    assert (c*p).as_content_primitive() == (c, p)


def test_Add_as_content_primitive():
    assert (x + 2).as_content_primitive() == (1, x + 2)

    assert (3*x + 2).as_content_primitive() == (1, 3*x + 2)
    assert (3*x + 3).as_content_primitive() == (3, x + 1)
    assert (3*x + 6).as_content_primitive() == (3, x + 2)

    assert (3*x + 2*y).as_content_primitive() == (1, 3*x + 2*y)
    assert (3*x + 3*y).as_content_primitive() == (3, x + y)
    assert (3*x + 6*y).as_content_primitive() == (3, x + 2*y)

    assert (3/x + 2*x*y*z**2).as_content_primitive() == (1, 3/x + 2*x*y*z**2)
    assert (3/x + 3*x*y*z**2).as_content_primitive() == (3, 1/x + x*y*z**2)
    assert (3/x + 6*x*y*z**2).as_content_primitive() == (3, 1/x + 2*x*y*z**2)

    assert (2*x/3 + 4*y/9).as_content_primitive() == \
        (Rational(2, 9), 3*x + 2*y)
    assert (2*x/3 + 2.5*y).as_content_primitive() == \
        (Rational(1, 3), 2*x + 7.5*y)

    # the coefficient may sort to a position other than 0
    p = 3 + x + y
    assert (2*p).expand().as_content_primitive() == (2, p)
    assert (2.0*p).expand().as_content_primitive() == (1, 2.*p)
    p *= -1
    assert (2*p).expand().as_content_primitive() == (2, p)


def test_Mul_as_content_primitive():
    assert (2*x).as_content_primitive() == (2, x)
    assert (x*(2 + 2*x)).as_content_primitive() == (2, x*(1 + x))
    assert (x*(2 + 2*y)*(3*x + 3)**2).as_content_primitive() == \
        (18, x*(1 + y)*(x + 1)**2)
    assert ((2 + 2*x)**2*(3 + 6*x) + Rational(1, 2)).as_content_primitive() == \
        (Rational(1, 2), 24*(x + 1)**2*(2*x + 1) + 1)


def test_Pow_as_content_primitive():
    assert (x**y).as_content_primitive() == (1, x**y)
    assert ((2*x + 2)**y).as_content_primitive() == \
        (1, (Mul(2, (x + 1), evaluate=False))**y)
    assert ((2*x + 2)**3).as_content_primitive() == (8, (x + 1)**3)
    assert (2**(Float(0.1) + x)).as_content_primitive() == (1, 2**(Float(0.1) + x))


def test_Pow_as_numer_denom():
    # issue sympy/sympy#10095
    assert ((1/(2*E))**oo).as_numer_denom() == (1, (2*E)**oo)
    assert ((2*E)**oo).as_numer_denom() == ((2*E)**oo, 1)
    e = Pow(1, oo, evaluate=False)
    assert e.as_numer_denom() == (e, 1)


def test_sympyissue_5460():
    u = Mul(2, (1 + x), evaluate=False)
    assert (2 + u).args == (2, u)


def test_product_irrational():
    assert (I*pi).is_irrational is False
    # The following used to be deduced from the above bug:
    assert (I*pi).is_positive is False


def test_sympyissue_5919():
    assert (x/(y*(1 + y))).expand() == x/(y**2 + y)


def test_Mod():
    assert isinstance(Mod(x, 1), Mod)
    assert pi % pi == 0
    assert Mod(5, 3) == 2
    assert Mod(-5, 3) == 1
    assert Mod(5, -3) == -1
    assert Mod(-5, -3) == -2
    assert type(Mod(3.2, 2, evaluate=False)) == Mod
    assert 5 % x == Mod(5, x)
    assert x % 5 == Mod(x, 5)
    assert x % y == Mod(x, y)
    assert (x % y).subs({x: 5, y: 3}) == 2

    # Float handling
    point3 = Float(3.3) % 1
    assert (x - 3.3) % 1 == Mod(1.*x + 1 - point3, 1)
    assert Mod(-3.3, 1) == 1 - point3
    assert Mod(0.7, 1) == Float(0.7)
    e = Mod(1.3, 1)
    assert comp(e, .3)
    assert e.is_Float
    e = Mod(1.3, .7)
    assert comp(e, .6)
    assert e.is_Float
    e = Mod(1.3, Rational(7, 10))
    assert comp(e, .6)
    assert e.is_Float
    e = Mod(Rational(13, 10), 0.7)
    assert comp(e, .6)
    assert e.is_Float
    e = Mod(Rational(13, 10), Rational(7, 10))
    assert comp(e, .6)
    assert e.is_Rational

    # check that sign is right
    r2 = sqrt(2)
    r3 = sqrt(3)
    for i in [-r3, -r2, r2, r3]:
        for j in [-r3, -r2, r2, r3]:
            assert verify_numerically(i % j, i.evalf() % j.evalf())
    for _x in range(4):
        for _y in range(9):
            reps = [(x, _x), (y, _y)]
            assert Mod(3*x + y, 9).subs(reps) == (3*_x + _y) % 9

    # denesting
    #   easy case
    assert Mod(Mod(x, y), y) == Mod(x, y)
    #   in case someone attempts more denesting
    for i in [-3, -2, 2, 3]:
        for j in [-3, -2, 2, 3]:
            for k in range(3):
                assert Mod(Mod(k, i), j) == (k % i) % j

    # known difference
    assert Mod(5*sqrt(2), sqrt(5)) == 5*sqrt(2) - 3*sqrt(5)
    p = symbols('p', positive=True)
    assert Mod(p + 1, p + 3) == p + 1
    assert Mod(x + 1, x + 3) == Mod(x + 1, x + 3, evaluate=False)
    n = symbols('n', negative=True)
    assert Mod(n - 3, n - 1) == -2
    assert Mod(n - 2*p, n - p) == -p
    assert Mod(p - 2*n, p - n) == -n

    # handling sums
    assert (x + 3) % 1 == Mod(x, 1)
    assert (x + 3.0) % 1 == Mod(1.*x, 1)
    assert (x - Rational(33, 10)) % 1 == Mod(x + Rational(7, 10), 1)

    a = Mod(.6*x + y, .3*y)
    b = Mod(0.1*y + 0.6*x, 0.3*y)
    # Test that a, b are equal, with 1e-14 accuracy in coefficients
    eps = 1e-14
    assert abs((a.args[0] - b.args[0]).subs({x: 1, y: 1})) < eps
    assert abs((a.args[1] - b.args[1]).subs({x: 1, y: 1})) < eps

    assert (x + 1) % x == 1 % x
    assert (x + y) % x == y % x
    assert (x + y + 2) % x == (y + 2) % x
    assert (a + 3*x + 1) % (2*x) == Mod(a + x + 1, 2*x)
    assert (12*x + 18*y) % (3*x) == 3*Mod(6*y, x)

    # gcd extraction
    assert (-3*x) % (-2*y) == -Mod(3*x, 2*y)
    assert (.6*pi) % (.3*x*pi) == 0.3*pi*Mod(2, x)
    assert (.6*pi) % (.31*x*pi) == pi*Mod(0.6, 0.31*x)
    assert (6*pi) % (.3*x*pi) == 0.3*pi*Mod(20, x)
    assert (6*pi) % (.31*x*pi) == pi*Mod(6, 0.31*x)
    assert (6*pi) % (.42*x*pi) == pi*Mod(6, 0.42*x)
    assert (12*x) % (2*y) == 2*Mod(6*x, y)
    assert (12*x) % (3*5*y) == 3*Mod(4*x, 5*y)
    assert (12*x) % (15*x*y) == 3*x*Mod(4, 5*y)
    assert (-2*pi) % (3*pi) == pi
    assert (2*x + 2) % (x + 1) == 0
    assert (x*(x + 1)) % (x + 1) == (x + 1)*Mod(x, 1)
    assert Mod(5.0*x, 0.1*y) == 0.1*Mod(50*x, y)
    i = Symbol('i', integer=True)
    assert (3*i*x) % (2*i*y) == i*Mod(3*x, 2*y)
    assert Mod(4*i, 4) == 0

    # issue sympy/sympy#8677
    n = Symbol('n', integer=True, positive=True)
    assert (factorial(n) % n).equals(0) is not False

    # symbolic with known parity
    n = Symbol('n', even=True)
    assert Mod(n, 2) == 0
    n = Symbol('n', odd=True)
    assert Mod(n, 2) == 1

    # issue diofant/diofant#312
    assert Mod(-x, 2*x) == x

    # issue sympy/sympy#10963
    assert (x**6000 % 400).args[1] == 400


def test_Mod_is_integer():
    p = Symbol('p', integer=True)
    q1 = Symbol('q1', integer=True)
    q2 = Symbol('q2', integer=True, nonzero=True)
    assert Mod(x, y).is_integer is None
    assert Mod(p, q1).is_integer is None
    assert Mod(x, q2).is_integer is None
    assert Mod(p, q2).is_integer


def test_Mod_is_nonposneg():
    n = Symbol('n', integer=True)
    k = Symbol('k', integer=True, positive=True)
    assert (n % 3).is_nonnegative
    assert Mod(n, -3).is_nonpositive
    assert Mod(n, k).is_nonnegative
    assert Mod(n, -k).is_nonpositive
    assert Mod(k, n).is_nonnegative is None


def test_sympyissue_6001():
    A = Symbol('A', commutative=False)
    eq = A + A**2
    # it doesn't matter whether it's True or False; they should
    # just all be the same
    assert eq.is_commutative == (eq + 1).is_commutative

    B = Symbol('B', commutative=False)
    # Although commutative terms could cancel we return True
    # meaning there are non-commutative symbols; aftersubstitution
    # that definition can change, e.g. (A*B).subs({B: A**-1}) -> 1
    assert (sqrt(2)*A).is_commutative is False
    assert (sqrt(2)*A*B).is_commutative is False


def test_polar():
    p = Symbol('p', polar=True)
    assert p.is_polar
    assert x.is_polar is None
    assert Integer(1).is_polar is None
    assert (p**x).is_polar is True
    assert (x**p).is_polar is None
    assert ((2*p)**x).is_polar is True
    assert (2*p).is_polar is True
    assert (-2*p).is_polar is not True
    assert (polar_lift(-2)*p).is_polar is True

    q = Symbol('q', polar=True)
    assert (p*q)**2 == p**2 * q**2
    assert (2*q)**2 == 4 * q**2
    assert ((p*q)**x).expand() == p**x * q**x


def test_sympyissue_6040():
    # pylint: disable=unneeded-not
    a, b = Pow(1, 2, evaluate=False), 1
    assert a != b
    assert b != a
    assert not a == b
    assert not b == a


def test_sympyissue_6077():
    assert x**2.0/x == x**1.0
    assert x/x**2.0 == x**-1.0
    assert x*x**2.0 == x**3.0
    assert x**1.5*x**2.5 == x**4.0

    assert 2**(2.0*x)/2**x == 2**(1.0*x)
    assert 2**x/2**(2.0*x) == 2**(-1.0*x)
    assert 2**x*2**(2.0*x) == 2**(3.0*x)
    assert 2**(1.5*x)*2**(2.5*x) == 2**(4.0*x)


def test_mul_flatten_oo_zoo():
    p = symbols('p', positive=True)
    n, m = symbols('n,m', negative=True)
    x_im = symbols('x_im', imaginary=True)
    assert n*oo == -oo
    assert n*m*oo == oo
    assert p*oo == oo
    assert x_im*oo != I*oo  # i could be +/- 3*I -> +/-oo

    assert zoo*2*zoo is zoo

    # issue sympy/sympy#18507
    assert Mul(zoo, zoo, 0) is nan


def test_add_flatten():
    # see https://github.com/sympy/sympy/issues/2633#issuecomment-29545524
    a = oo + I*oo
    b = oo - I*oo
    assert a + b == nan
    assert a - b == nan
    assert (1/a).simplify() == (1/b).simplify() == 0

    a = Pow(2, 3, evaluate=False)
    assert a + a == 16

    assert zoo + 1 + zoo is nan


def test_issue_31():
    assert sin(x + O(x**2)) - sin(x + O(x**2)) == \
        Add(-sin(x + O(x**2)), sin(x + O(x**2)), evaluate=False)
    assert sin(O(x))/sin(O(x)) == Mul(1/sin(O(x)), sin(O(x)), evaluate=False)


def test_sympyissue_5160_6087_6089_6090():
    # issue sympy/sympy#6087
    assert ((-2*x*y**y)**3.2).evalf(2, strict=False) == (2**3.2*(-x*y**y)**3.2).evalf(2, strict=False)
    # issue sympy/sympy#6089
    A, B, C = symbols('A,B,C', commutative=False)
    assert (2.*B*C)**3 == 8.0*(B*C)**3
    assert (-2.*B*C)**3 == -8.0*(B*C)**3
    assert (-2*B*C)**2 == 4*(B*C)**2
    # issue sympy/sympy#5160
    assert sqrt(-1.0*x) == 1.0*sqrt(-x)
    assert sqrt(1.0*x) == 1.0*sqrt(x)
    # issue sympy/sympy#6090
    assert (-2*x*y*A*B)**2 == 4*x**2*y**2*(A*B)**2


def test_float_int():
    assert int(float(sqrt(10))) == int(sqrt(10))
    assert int(pi**1000) % 10 == 2
    assert int(Float('1.123456789012345678901234567890e20')) == \
        int(112345678901234567890)
    assert int(Float('1.123456789012345678901234567890e25')) == \
        int(11234567890123456789012345)
    # decimal forces float so it's not an exact integer ending in 000000
    assert int(Float('1.123456789012345678901234567890e35')) == \
        112345678901234567890123456789000192
    assert int(Float('123456789012345678901234567890e5')) == \
        12345678901234567890123456789000192
    assert Integer(Float('1.123456789012345678901234567890e20')) == \
        112345678901234567890
    assert Integer(Float('1.123456789012345678901234567890e25')) == \
        11234567890123456789012345
    # decimal forces float so it's not an exact integer ending in 000000
    assert Integer(Float('1.123456789012345678901234567890e35')) == \
        112345678901234567890123456789000192
    assert Integer(Float('123456789012345678901234567890e5')) == \
        12345678901234567890123456789000192
    assert same_and_same_prec(Float('123000e-2'), Float('1230.00'))
    assert same_and_same_prec(Float('123000e2'), Float('123000.e2'))

    assert int(1 + Rational('.9999999999999999999999999')) == 1
    assert int(pi/1e20) == 0
    assert int(1 + pi/1e20) == 1
    assert int(Add(1.2, -2, evaluate=False)) == int(1.2 - 2)
    assert int(Add(1.2, +2, evaluate=False)) == int(1.2 + 2)
    assert int(Add(1 + Float('.99999999999999999'), evaluate=False)) == 1
    pytest.raises(TypeError, lambda: float(x))
    pytest.raises(TypeError, lambda: float(sqrt(-1)))

    assert int(12345678901234567890 + cos(1)**2 + sin(1)**2) == \
        12345678901234567891


def test_sympyissue_6611a():
    assert Mul.flatten([cbrt(3),
                        Pow(-Rational(1, 9), Rational(2, 3), evaluate=False)]) == \
        ([Rational(1, 3), (-1)**Rational(2, 3)], [], None)


def test_denest_add_mul():
    # when working with evaluated expressions make sure they denest
    eq = x + 1
    eq = Add(eq, 2, evaluate=False)
    eq = Add(eq, 2, evaluate=False)
    assert Add(*eq.args) == x + 5
    eq = x*2
    eq = Mul(eq, 2, evaluate=False)
    eq = Mul(eq, 2, evaluate=False)
    assert Mul(*eq.args) == 8*x
    # but don't let them denest unecessarily
    eq = Mul(-2, x - 2, evaluate=False)
    assert 2*eq == Mul(-4, x - 2, evaluate=False)
    assert -eq == Mul(2, x - 2, evaluate=False)


def test_mul_coeff():
    # It is important that all Numbers be removed from the seq;
    # This can be tricky when powers combine to produce those numbers
    p = exp(I*pi/3)
    assert p**2*x*p*y*p*x*p**2 == x**2*y


def test_mul_zero_detection():
    nz = Dummy(real=True, nonzero=True)
    r = Dummy(extended_real=True)
    c = Dummy(real=False, complex=True)
    c2 = Dummy(real=False, complex=True)
    i = Dummy(imaginary=True)
    ni = Dummy(imaginary=True, nonzero=True)
    e = nz*r*c
    assert e.is_imaginary is None
    assert e.is_extended_real is None
    e = nz*c
    assert e.is_imaginary is None
    assert e.is_extended_real is False
    e = nz*ni*c
    assert e.is_imaginary is False
    assert e.is_extended_real is None
    # check for more than one complex; it is important to use
    # uniquely named Symbols to ensure that two factors appear
    # e.g. if the symbols have the same name they just become
    # a single factor, a power.
    e = nz*i*c*c2
    assert e.is_imaginary is None
    assert e.is_extended_real is None

    # _eval_is_extended_real and _eval_is_zero both employ trapping of the
    # zero value so args should be tested in both directions and
    # TO AVOID GETTING THE CACHED RESULT, Dummy MUST BE USED

    # real is unknonwn
    def test(z, b, e):
        if z.is_zero and b.is_finite:
            assert e.is_extended_real
            assert e.is_zero
        else:
            assert e.is_extended_real is None
            if b.is_finite:
                if z.is_zero:
                    assert e.is_zero
                else:
                    assert e.is_zero is None
            elif b.is_finite is False:
                if z.is_zero is None:
                    assert e.is_zero is None
                else:
                    assert e.is_nonzero

    for iz, ib in itertools.product(*[[True, False, None]]*2):
        z = Dummy('z', nonzero=iz)
        b = Dummy('f', finite=ib)
        e = Mul(z, b, evaluate=False)
        test(z, b, e)
        z = Dummy('nz', nonzero=iz)
        b = Dummy('f', finite=ib)
        e = Mul(b, z, evaluate=False)
        test(z, b, e)

    # real is True
    def test2(z, b, e):
        if z.is_zero and not b.is_finite:
            assert e.is_extended_real is None
        elif not z.is_finite:
            return e.is_extended_real is None
        else:
            assert e.is_extended_real

    for iz, ib in itertools.product(*[[True, False, None]]*2):
        z = Dummy('z', nonzero=iz, extended_real=True)
        b = Dummy('b', finite=ib, extended_real=True)
        e = Mul(z, b, evaluate=False)
        test2(z, b, e)
        z = Dummy('z', nonzero=iz, extended_real=True)
        b = Dummy('b', finite=ib, extended_real=True)
        e = Mul(b, z, evaluate=False)
        test2(z, b, e)


def test_sympyissue_8274():
    z = sqrt(1 + sqrt(3)) + sqrt(3 + 3*sqrt(3)) - sqrt(10 + 6*sqrt(3))
    assert z.is_positive is False  # it's 0
    z = (-cbrt(2)*(3*sqrt(93) + 29)**2 -
         4*(3*sqrt(93) + 29)**Rational(4, 3) +
         12*sqrt(93)*cbrt(3*sqrt(93) + 29) +
         116*cbrt(3*sqrt(93) + 29) +
         174*cbrt(2)*sqrt(93) + 1678*cbrt(2))
    assert z.is_positive is False  # it's 0
    z = 2*(-3*tan(19*pi/90) + sqrt(3))*cos(11*pi/90)*cos(19*pi/90) - \
        sqrt(3)*(-3 + 4*cos(19*pi/90)**2)
    assert z.is_positive is not True  # it's zero and it shouldn't hang


def test_sympyissue_8354():
    z = (9*(3*sqrt(93) + 29)**Rational(2, 3)*(cbrt(3*sqrt(93) +
                                                   29)*(-2**Rational(2, 3)*cbrt(3*sqrt(93) +
                                                                                29) - 2) - 2*cbrt(2))**3 +
         72*(3*sqrt(93) + 29)**Rational(2, 3)*(81*sqrt(93) + 783) +
         (162*sqrt(93) + 1566)*(cbrt(3*sqrt(93) + 29) *
                                (-2**Rational(2, 3)*cbrt(3*sqrt(93) + 29) - 2) -
                                2*cbrt(2))**2)
    assert z.is_positive is False  # it's 0 (and a single _mexpand isn't enough)


def test_sympyissue_9832():
    x = Symbol('x', extended_real=True)
    assert (x**2 - oo).is_negative is None


def test_sympyissue_10728():
    A, B = symbols('A B', commutative=False)
    assert (A + B).is_commutative is None
    assert (A + B).is_zero is None


def test_sympyissue_18509():
    e = 2**oo / pi**oo

    assert e != oo
    assert e == Mul(oo, pi**-oo, evaluate=False)

    e = 2**oo / (E + 1)**oo

    assert e != oo
    assert e == Mul(oo, (E + 1)**-oo, evaluate=False)


def test_sympyissue_16971():
    a = Symbol('a', extended_real=True)
    b = Symbol('b', extended_real=True)

    assert (a + b).is_extended_real is None
    assert (a - b).is_extended_real is None


def test_issue_849():
    a = Symbol('a', extended_real=True)
    b = Symbol('b', extended_real=True)

    # issue sympy/sympy#16971
    assert (a + b).is_extended_real is None
    assert (a - b).is_extended_real is None

    assert (a*b).is_extended_real is None


def test_issue_1004():
    assert Pow(Dummy(negative=True), -3,
               evaluate=False).is_negative is not True
    assert Pow(-oo, -3, evaluate=False).is_negative is not True
