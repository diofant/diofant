import pytest

from diofant import (Abs, Derivative, DiracDelta, E, Eq, Expr, Function,
                     Heaviside, I, Integer, Integral, Interval, Matrix, Ne, O,
                     Piecewise, Rational, Subs, Symbol, adjoint, arg, atan2,
                     cbrt, comp, conjugate, cos, erf, exp, exp_polar, expand,
                     gamma, im, log, nan, oo, periodic_argument, pi,
                     polar_lift, polarify, principal_branch, re, root, sign,
                     simplify, sin, sqrt, symbols, tanh, transpose,
                     unbranched_argument, unpolarify, uppergamma, zoo)
from diofant.abc import x, y, z
from diofant.core.function import ArgumentIndexError


__all__ = ()


def N_equals(a, b):
    """Check whether two complex numbers are numerically close."""
    return comp(a.evalf(), b.evalf(), 1.e-6)


def test_re():
    a, b = symbols('a,b', extended_real=True)

    r = Symbol('r', extended_real=True)
    i = Symbol('i', imaginary=True)

    assert re(nan) == nan

    assert re(oo) == oo
    assert re(-oo) == -oo

    assert re(0) == 0

    assert re(1) == 1
    assert re(-1) == -1

    assert re(E) == E
    assert re(-E) == -E

    x = Symbol('x')
    assert re(x) == re(x)
    assert re(x*I) == -im(x)
    assert re(r*I) == 0
    assert re(r) == r
    assert re(i*I) == I * i
    assert re(i) == 0

    y = Symbol('y')
    assert re(x + y) == re(x + y)
    assert re(x + r) == re(x) + r

    assert re(re(x)) == re(x)

    assert re(2 + I) == 2
    assert re(x + I) == re(x)

    assert re(x + y*I) == re(x) - im(y)
    assert re(x + r*I) == re(x)

    assert re(log(2*I)) == log(2)

    assert re((2 + I)**2).expand(complex=True) == 3

    assert re(conjugate(x)) == re(x)
    assert conjugate(re(x)) == re(x)

    assert re(x).as_real_imag() == (re(x), 0)

    assert re(i*r*x).diff(r) == re(i*x)
    assert re(i*r*x).diff(i) == I*r*im(x)

    assert re(sqrt(a + b*I)) == root(a**2 + b**2, 4)*cos(arg(a + I*b)/2)
    assert re(a * (2 + b*I)) == 2*a

    assert re((1 + sqrt(a + b*I))/2) == root(a**2 + b**2, 4)*cos(arg(a + I*b)/2)/2 + Rational(1, 2)

    assert re(x).rewrite(im) == x - I*im(x)  # issue sympy/sympy#10897
    assert (x + re(y)).rewrite(re, im) == x + y - I*im(y)

    a = Symbol('a', algebraic=True)
    t = Symbol('t', transcendental=True)
    assert re(a).is_algebraic
    assert re(x).is_algebraic is None
    assert re(t).is_algebraic is False

    assert re(zoo) == nan

    # issue sympy/sympy#4757
    x = Symbol('x', extended_real=True)
    y = Symbol('y', imaginary=True)
    f = Function('f')
    assert re(f(x)).diff(x) == re(f(x).diff(x))
    assert re(f(y)).diff(y) == -I*im(f(y).diff(y))
    assert re(f(z)).diff(z) == Derivative(re(f(z)), z)


def test_im():
    a, b = symbols('a,b', extended_real=True)

    r = Symbol('r', extended_real=True)
    i = Symbol('i', imaginary=True)

    assert im(nan) == nan

    assert im(oo*I) == oo
    assert im(-oo*I) == -oo

    assert im(0) == 0

    assert im(1) == 0
    assert im(-1) == 0

    assert im(E*I) == E
    assert im(-E*I) == -E

    x = Symbol('x')
    assert im(x) == im(x)
    assert im(x*I) == re(x)
    assert im(r*I) == r
    assert im(r) == 0
    assert im(i*I) == 0
    assert im(i) == -I * i

    y = Symbol('x')
    assert im(x + y) == im(x + y)
    assert im(x + r) == im(x)
    assert im(x + r*I) == im(x) + r

    assert im(im(x)*I) == im(x)

    assert im(2 + I) == 1
    assert im(x + I) == im(x) + 1

    assert im(x + y*I) == im(x) + re(y)
    assert im(x + r*I) == im(x) + r

    assert im(log(2*I)) == pi/2

    assert im((2 + I)**2).expand(complex=True) == 4

    assert im(conjugate(x)) == -im(x)
    assert conjugate(im(x)) == im(x)

    assert im(x).as_real_imag() == (im(x), 0)

    assert im(i*r*x).diff(r) == im(i*x)
    assert im(i*r*x).diff(i) == -I * re(r*x)

    assert im(sqrt(a + b*I)) == root(a**2 + b**2, 4)*sin(arg(a + I*b)/2)
    assert im(a * (2 + b*I)) == a*b

    assert im((1 + sqrt(a + b*I))/2) == root(a**2 + b**2, 4)*sin(arg(a + I*b)/2)/2

    assert im(x).rewrite(re) == -I*(x - re(x))  # sympy/sympy#10897
    assert (x + im(y)).rewrite(im, re) == x - I*(y - re(y))

    a = Symbol('a', algebraic=True)
    t = Symbol('t', transcendental=True)
    assert re(a).is_algebraic
    assert re(x).is_algebraic is None
    assert re(t).is_algebraic is False

    assert re(zoo) == nan

    # issue sympy/sympy#4757
    x = Symbol('x', extended_real=True)
    y = Symbol('y', imaginary=True)
    f = Function('f')
    assert im(f(x)).diff(x) == im(f(x).diff(x))
    assert im(f(y)).diff(y) == -I*re(f(y).diff(y))
    assert im(f(z)).diff(z) == Derivative(im(f(z)), z)


def test_sign():
    assert sign(1.2) == 1
    assert sign(-1.2) == -1
    assert sign(3*I) == I
    assert sign(-3*I) == -I
    assert sign(0) == 0
    assert sign(nan) == nan
    assert sign(2 + 2*I).doit() == sqrt(2)*(2 + 2*I)/4
    assert sign(2 + 3*I).simplify() == sign(2 + 3*I)
    assert sign(2 + 2*I).simplify() == sign(1 + I)
    assert sign(im(sqrt(1 - sqrt(3)))) == 1
    assert sign(sqrt(1 - sqrt(3))) == I

    x = Symbol('x')
    assert sign(x).is_finite is True
    assert sign(x).is_complex is True
    assert sign(x).is_imaginary is None
    assert sign(x).is_integer is None
    assert sign(x).is_extended_real is None
    assert sign(x).is_zero is None
    assert sign(x).doit() == sign(x)
    assert sign(1.2*x) == sign(x)
    assert sign(2*x) == sign(x)
    assert sign(I*x) == I*sign(x)
    assert sign(-2*I*x) == -I*sign(x)
    assert sign(conjugate(x)) == conjugate(sign(x))

    p = Symbol('p', positive=True)
    n = Symbol('n', negative=True)
    m = Symbol('m', negative=True)
    assert sign(2*p*x) == sign(x)
    assert sign(n*x) == -sign(x)
    assert sign(n*m*x) == sign(x)

    x = Symbol('x', imaginary=True)
    xn = Symbol('xn', imaginary=True, nonzero=True)
    assert sign(x).is_imaginary is True
    assert sign(x).is_integer is None
    assert sign(x).is_extended_real is None
    assert sign(x).is_zero is None
    assert sign(x).diff(x) == 2*DiracDelta(-I*x)
    assert sign(xn).doit() == xn/abs(xn)
    assert conjugate(sign(x)) == -sign(x)

    x = Symbol('x', extended_real=True)
    assert sign(x).is_imaginary is None
    assert sign(x).is_integer is True
    assert sign(x).is_extended_real is True
    assert sign(x).is_zero is None
    assert sign(x).diff(x) == 2*DiracDelta(x)
    assert sign(x).doit() == sign(x)
    assert conjugate(sign(x)) == sign(x)

    assert sign(sin(x)).series(x) == 1
    y = Symbol('y')
    assert sign(x*y).series(x).removeO() == sign(y)
    assert sign(I + x).series(x, n=2) == I + x*Subs(sign(x + I).diff(x),
                                                    (x, 0)) + O(x**2)

    x = Symbol('x', nonzero=True)
    assert sign(x).is_imaginary is None
    assert sign(x).is_integer is None
    assert sign(x).is_extended_real is None
    assert sign(x).is_nonzero is True
    assert sign(x).doit() == x/abs(x)
    assert sign(abs(x)) == 1
    assert abs(sign(x)) == 1

    x = Symbol('x', positive=True)
    assert sign(x).is_imaginary is False
    assert sign(x).is_integer is True
    assert sign(x).is_extended_real is True
    assert sign(x).is_nonzero is True
    assert sign(x).doit() == x/abs(x)
    assert sign(abs(x)) == 1
    assert abs(sign(x)) == 1

    x = 0
    assert sign(x).is_imaginary is True
    assert sign(x).is_integer is True
    assert sign(x).is_extended_real is True
    assert sign(x).is_zero is True
    assert sign(x).doit() == 0
    assert sign(abs(x)) == 0
    assert abs(sign(x)) == 0

    nz = Symbol('nz', nonzero=True, integer=True)
    assert sign(nz).is_imaginary is False
    assert sign(nz).is_integer is True
    assert sign(nz).is_extended_real is True
    assert sign(nz).is_nonzero is True
    assert sign(nz)**2 == 1
    assert (sign(nz)**3).args == (sign(nz), 3)

    assert sign(Symbol('x', nonnegative=True)).is_nonnegative
    assert sign(Symbol('x', nonnegative=True)).is_nonpositive is None
    assert sign(Symbol('x', nonpositive=True)).is_nonnegative is None
    assert sign(Symbol('x', nonpositive=True)).is_nonpositive
    assert sign(Symbol('x', extended_real=True)).is_nonnegative is None
    assert sign(Symbol('x', extended_real=True)).is_nonpositive is None
    assert sign(Symbol('x', extended_real=True, zero=False)).is_nonpositive is None

    x, y = Symbol('x', extended_real=True), Symbol('y')
    assert sign(x).rewrite(Piecewise) == \
        Piecewise((1, x > 0), (-1, x < 0), (0, True))
    assert sign(y).rewrite(Piecewise) == sign(y)
    assert sign(x).rewrite(Heaviside) == 2*Heaviside(x)-1
    assert sign(y).rewrite(Heaviside) == sign(y)
    assert sign(1 + I).rewrite(exp) == exp(I*pi/4)

    # evaluate what can be evaluated
    assert sign(exp_polar(I*pi)*pi) is Integer(-1)

    eq = -sqrt(10 + 6*sqrt(3)) + sqrt(1 + sqrt(3)) + sqrt(3 + 3*sqrt(3))
    # if there is a fast way to know when and when you cannot prove an
    # expression like this is zero then the equality to zero is ok
    assert sign(eq) == 0
    q = 1 + sqrt(2) - 2*sqrt(3) + 1331*sqrt(6)
    p = cbrt(expand(q**3))
    d = p - q
    assert sign(d) == 0

    assert abs(sign(z)) == Abs(sign(z), evaluate=False)


def test_as_real_imag():
    n = pi**1000
    # the special code for working out the real
    # and complex parts of a power with Integer exponent
    # should not run if there is no imaginary part, hence
    # this should not hang
    assert n.as_real_imag() == (n, 0)

    # issue sympy/sympy#6261
    assert sqrt(x).as_real_imag() == \
        (root(re(x)**2 + im(x)**2, 4)*cos(arg(re(x) + I*im(x))/2),
         root(re(x)**2 + im(x)**2, 4)*sin(arg(re(x) + I*im(x))/2))

    # issue sympy/sympy#3853
    a, b = symbols('a,b', extended_real=True)
    assert (((1 + sqrt(a + b*I))/2).as_real_imag() ==
            (root(a**2 + b**2, 4)*cos(arg(a + I*b)/2)/2 + Rational(1, 2),
             root(a**2 + b**2, 4)*sin(arg(a + I*b)/2)/2))

    assert sqrt(a**2).as_real_imag() == (sqrt(a**2), 0)
    i = symbols('i', imaginary=True)
    assert sqrt(i**2).as_real_imag() == (0, abs(i))


@pytest.mark.xfail
def test_sign_sympyissue_6167():
    n = pi**1000
    i = int(n)
    assert (n - i).round() == 1  # doesn't hang
    assert sign(n - i) == 1
    # perhaps it's not possible to get the sign right when
    # only 1 digit is being requested for this situtation;
    # 2 digits works
    assert (n - x).evalf(1, subs={x: i}, maxn=400) > 0
    assert (n - x).evalf(2, subs={x: i}, maxn=400) > 0


def test_Abs():
    pytest.raises(TypeError, lambda: Abs(Interval(2, 3)))  # issue sympy/sympy#8717

    x, y = symbols('x,y')
    assert sign(sign(x)) == sign(x)
    assert isinstance(sign(x*y), sign)
    assert Abs(0) == 0
    assert Abs(1) == 1
    assert Abs(-1) == 1
    assert abs(I) == 1
    assert abs(-I) == 1
    assert abs(nan) == nan
    assert abs(zoo) == oo
    assert abs(+I * pi) == pi
    assert abs(-I * pi) == pi
    assert abs(+I * x) == abs(x)
    assert abs(-I * x) == abs(x)
    assert abs(-2*x) == 2*abs(x)
    assert abs(-2.0*x) == 2.0*abs(x)
    assert abs(2*pi*x*y) == 2*pi*abs(x*y)
    assert abs(conjugate(x)) == abs(x)
    assert conjugate(abs(x)) == abs(x)

    a = cos(1)**2 + sin(1)**2 - 1
    assert abs(a*x).series(x).simplify() == 0

    a = Symbol('a', positive=True)
    assert abs(2*pi*x*a) == 2*pi*a*abs(x)
    assert abs(2*pi*I*x*a) == 2*pi*a*abs(x)

    x = Symbol('x', extended_real=True)
    n = Symbol('n', integer=True)
    assert abs((-1)**n) == 1
    assert x**(2*n) == abs(x)**(2*n)
    assert abs(x).diff(x) == sign(x)
    assert abs(-x).fdiff() == sign(x)
    assert abs(x) == Abs(x)  # Python built-in
    assert abs(x)**3 == x**2*abs(x)
    assert abs(x)**4 == x**4
    assert (abs(x)**(3*n)).args == (abs(x), 3*n)  # leave symbolic odd unchanged
    assert (1/abs(x)).args == (abs(x), -1)
    assert 1/abs(x)**3 == 1/(x**2*abs(x))
    assert abs(x)**-3 == abs(x)/x**4
    assert abs(x**3) == x**2*abs(x)
    assert abs(x**pi) == Abs(x**pi, evaluate=False)

    x = Symbol('x', imaginary=True)
    assert abs(x).diff(x) == -sign(x)

    pytest.raises(ArgumentIndexError, lambda: abs(z).fdiff(2))

    eq = -sqrt(10 + 6*sqrt(3)) + sqrt(1 + sqrt(3)) + sqrt(3 + 3*sqrt(3))
    # if there is a fast way to know when you can and when you cannot prove an
    # expression like this is zero then the equality to zero is ok
    assert abs(eq) == 0
    q = 1 + sqrt(2) - 2*sqrt(3) + 1331*sqrt(6)
    p = cbrt(expand(q**3))
    d = p - q
    assert abs(d) == 0

    assert abs(4*exp(pi*I/4)) == 4
    assert abs(3**(2 + I)) == 9
    assert abs((-3)**(1 - I)) == 3*exp(pi)

    assert abs(+oo) is oo
    assert abs(-oo) is oo
    assert abs(oo + I) is oo
    assert abs(oo + I*oo) is oo

    a = Symbol('a', algebraic=True)
    t = Symbol('t', transcendental=True)
    x = Symbol('x')
    assert re(a).is_algebraic
    assert re(x).is_algebraic is None
    assert re(t).is_algebraic is False

    assert abs(sign(z)) == Abs(sign(z), evaluate=False)

    # issue sympy/sympy#4757
    x = Symbol('x', extended_real=True)
    y = Symbol('y', imaginary=True)
    f = Function('f')
    assert abs(f(x)).diff(x).subs({f(x): 1 + I*x}).doit() == x/sqrt(1 + x**2)
    assert abs(f(y)).diff(y).subs({f(y): 1 + y}).doit() == -y/sqrt(1 - y**2)


def test_Abs_rewrite():
    x = Symbol('x', extended_real=True)
    a = abs(x).rewrite(Heaviside).expand()
    assert a == x*Heaviside(x) - x*Heaviside(-x)
    for i in [-2, -1, 0, 1, 2]:
        assert a.subs({x: i}) == abs(i)
    y = Symbol('y')
    assert abs(y).rewrite(Heaviside) == abs(y)

    x, y = Symbol('x', extended_real=True), Symbol('y')
    assert abs(x).rewrite(Piecewise) == Piecewise((x, x >= 0), (-x, True))
    assert abs(y).rewrite(Piecewise) == abs(y)
    assert abs(y).rewrite(sign) == y/sign(y)


def test_Abs_real():
    # test some properties of abs that only apply
    # to real numbers
    x = Symbol('x', complex=True)
    assert sqrt(x**2) != abs(x)
    assert abs(x**2) != x**2

    x = Symbol('x', extended_real=True)
    assert sqrt(x**2) == abs(x)
    assert abs(x**2) == x**2

    # if the symbol is zero, the following will still apply
    nn = Symbol('nn', nonnegative=True, extended_real=True)
    np = Symbol('np', nonpositive=True, extended_real=True)
    assert abs(nn) == nn
    assert abs(np) == -np


def test_Abs_properties():
    x, z = symbols('x, z')
    assert abs(x).is_extended_real is True
    assert abs(x).is_rational is None
    assert abs(x).is_positive is None
    assert abs(x).is_nonnegative is True
    assert abs(x).is_finite is None

    z = Symbol('z', complex=True, zero=False)
    assert abs(z).is_extended_real is True
    assert abs(z).is_rational is None
    assert abs(z).is_positive is True
    assert abs(z).is_nonzero is True
    assert abs(z).is_finite is True

    p = Symbol('p', positive=True)
    assert abs(p).is_extended_real is True
    assert abs(p).is_rational is None
    assert abs(p).is_positive is True
    assert abs(p).is_nonzero is True

    q = Symbol('q', rational=True)
    assert abs(q).is_rational is True
    assert abs(q).is_integer is None
    assert abs(q).is_positive is None
    assert abs(q).is_nonnegative is True
    assert abs(q).is_finite is True

    i = Symbol('i', integer=True)
    assert abs(i).is_integer is True
    assert abs(i).is_positive is None
    assert abs(i).is_nonnegative is True
    assert abs(i).is_finite is True

    e = Symbol('n', even=True)
    ne = Symbol('ne', extended_real=True, even=False)
    assert abs(e).is_even
    assert abs(ne).is_even is False
    assert abs(i).is_even is None

    o = Symbol('n', odd=True)
    no = Symbol('no', extended_real=True, odd=False)
    assert abs(o).is_odd
    assert abs(no).is_odd is False
    assert abs(i).is_odd is None


def test_abs():
    # this tests that abs calls Abs; don't rename to
    # test_Abs since that test is already above
    a = Symbol('a', positive=True)
    assert abs(I*(1 + a)**2) == (1 + a)**2


def test_arg():
    assert arg(0) == 0
    assert arg(1) == 0
    assert arg(-1) == pi
    assert arg(I) == pi/2
    assert arg(-I) == -pi/2
    assert arg(1 + I) == pi/4
    assert arg(-1 + I) == 3*pi/4
    assert arg(1 - I) == -pi/4
    f = Function('f')
    assert not arg(f(0) + I*f(1)).atoms(re)

    p = Symbol('p', positive=True)
    assert arg(p) == 0

    n = Symbol('n', negative=True)
    assert arg(n) == pi

    x = Symbol('x')
    assert conjugate(arg(x)) == arg(x)

    e = p + I*p**2
    assert arg(e) == arg(1 + p*I)
    # make sure sign doesn't swap
    e = -2*p + 4*I*p**2
    assert arg(e) == arg(-1 + 2*p*I)
    # make sure sign isn't lost
    x = symbols('x', extended_real=True)  # could be zero
    e = x + I*x
    assert arg(e) == arg(x*(1 + I))
    assert arg(e/p) == arg(x*(1 + I))
    e = p*cos(p) + I*log(p)*exp(p)
    assert arg(e).args[0] == e
    # keep it simple -- let the user do more advanced cancellation
    e = (p + 1) + I*(p**2 - 1)
    assert arg(e).args[0] == e


def test_arg_rewrite():
    assert arg(1 + I) == atan2(1, 1)

    x = Symbol('x', extended_real=True)
    y = Symbol('y', extended_real=True)
    assert arg(x + I*y).rewrite(atan2) == atan2(y, x)


def test_adjoint():
    x, y = symbols('x y')
    assert adjoint(adjoint(x)) == x
    assert adjoint(x + y) == adjoint(x) + adjoint(y)
    assert adjoint(x - y) == adjoint(x) - adjoint(y)
    assert adjoint(x * y) == adjoint(x) * adjoint(y)
    assert adjoint(x / y) == adjoint(x) / adjoint(y)
    assert adjoint(-x) == -adjoint(x)

    x, y = symbols('x y', commutative=False)
    assert adjoint(adjoint(x)) == x
    assert adjoint(x + y) == adjoint(x) + adjoint(y)
    assert adjoint(x - y) == adjoint(x) - adjoint(y)
    assert adjoint(x * y) == adjoint(y) * adjoint(x)
    assert adjoint(x / y) == 1 / adjoint(y) * adjoint(x)
    assert adjoint(-x) == -adjoint(x)


def test_conjugate():
    a = Symbol('a', extended_real=True)
    b = Symbol('b', imaginary=True)
    assert conjugate(a) == a
    assert conjugate(I*a) == -I*a
    assert conjugate(b) == -b
    assert conjugate(I*b) == I*b
    assert conjugate(a*b) == -a*b
    assert conjugate(I*a*b) == I*a*b

    x = Symbol('x')
    y = Symbol('y')
    assert conjugate(conjugate(x)) == x
    assert conjugate(x + y) == conjugate(x) + conjugate(y)
    assert conjugate(x - y) == conjugate(x) - conjugate(y)
    assert conjugate(x * y) == conjugate(x) * conjugate(y)
    assert conjugate(x / y) == conjugate(x) / conjugate(y)
    assert conjugate(-x) == -conjugate(x)

    a = Symbol('a', algebraic=True)
    t = Symbol('t', transcendental=True)
    assert re(a).is_algebraic
    assert re(x).is_algebraic is None
    assert re(t).is_algebraic is False

    assert conjugate(z).diff(z) == Derivative(conjugate(z), z)

    # issue sympy/sympy#4754
    x = Symbol('x', extended_real=True)
    y = Symbol('y', imaginary=True)
    f = Function('f')
    assert (f(x).conjugate()).diff(x) == (f(x).diff(x)).conjugate()
    assert (f(y).conjugate()).diff(y) == -(f(y).diff(y)).conjugate()


def test_conjugate_transpose():
    x = Symbol('x')
    assert conjugate(transpose(x)) == adjoint(x)
    assert transpose(conjugate(x)) == adjoint(x)
    assert adjoint(transpose(x)) == conjugate(x)
    assert transpose(adjoint(x)) == conjugate(x)
    assert adjoint(conjugate(x)) == transpose(x)
    assert conjugate(adjoint(x)) == transpose(x)

    class Symmetric(Expr):
        def _eval_adjoint(self):
            return

        def _eval_conjugate(self):
            return

        def _eval_transpose(self):
            return self
    x = Symmetric()
    assert conjugate(x) == adjoint(x)
    assert transpose(x) == x


def test_transpose():
    a = Symbol('a', complex=True)
    assert transpose(a) == a
    assert transpose(I*a) == I*a

    x, y = symbols('x y')
    assert transpose(transpose(x)) == x
    assert transpose(x + y) == transpose(x) + transpose(y)
    assert transpose(x - y) == transpose(x) - transpose(y)
    assert transpose(x * y) == transpose(x) * transpose(y)
    assert transpose(x / y) == transpose(x) / transpose(y)
    assert transpose(-x) == -transpose(x)

    x, y = symbols('x y', commutative=False)
    assert transpose(transpose(x)) == x
    assert transpose(x + y) == transpose(x) + transpose(y)
    assert transpose(x - y) == transpose(x) - transpose(y)
    assert transpose(x * y) == transpose(y) * transpose(x)
    assert transpose(x / y) == 1 / transpose(y) * transpose(x)
    assert transpose(-x) == -transpose(x)


def test_polarify():
    z = Symbol('z', polar=True)
    f = Function('f')
    ES = {}

    assert polarify(-1) == (polar_lift(-1), ES)
    assert polarify(1 + I) == (polar_lift(1 + I), ES)

    assert polarify(exp(x), subs=False) == exp(x)
    assert polarify(1 + x, subs=False) == 1 + x
    assert polarify(f(I) + x, subs=False) == f(polar_lift(I)) + x

    assert polarify(x, lift=True) == polar_lift(x)
    assert polarify(z, lift=True) == z
    assert polarify(f(x), lift=True) == f(polar_lift(x))
    assert polarify(1 + x, lift=True) == polar_lift(1 + x)
    assert polarify(1 + f(x), lift=True) == polar_lift(1 + f(polar_lift(x)))

    newex, subs = polarify(f(x) + z)
    assert newex.subs(subs) == f(x) + z

    mu = Symbol('mu')
    sigma = Symbol('sigma', positive=True)

    # Make sure polarify(lift=True) doesn't try to lift the integration
    # variable
    assert polarify(
        Integral(sqrt(2)*x*exp(-(-mu + x)**2/(2*sigma**2))/(2*sqrt(pi)*sigma),
                 (x, -oo, oo)), lift=True) == Integral(sqrt(2)*(sigma*exp_polar(0))**exp_polar(I*pi) *
                                                       exp((sigma*exp_polar(0))**(2*exp_polar(I*pi))*exp_polar(I*pi)*polar_lift(-mu + x) **
                                                           (2*exp_polar(0))/2)*exp_polar(0)*polar_lift(x)/(2*sqrt(pi)), (x, -oo, oo))


def test_unpolarify():
    p = exp_polar(7*I) + 1
    u = exp(7*I) + 1

    assert unpolarify(1) == 1
    assert unpolarify(p) == u
    assert unpolarify(p**2) == u**2
    assert unpolarify(p**x) == p**x
    assert unpolarify(p*x) == u*x
    assert unpolarify(p + x) == u + x
    assert unpolarify(sqrt(sin(p))) == sqrt(sin(u))

    # Test reduction to principal branch 2*pi.
    t = principal_branch(x, 2*pi)
    assert unpolarify(t) == x
    assert unpolarify(sqrt(t)) == sqrt(t)

    # Test exponents_only.
    assert unpolarify(p**p, exponents_only=True) == p**u
    assert unpolarify(uppergamma(x, p**p)) == uppergamma(x, p**u)

    # Test functions.
    assert unpolarify(sin(p)) == sin(u)
    assert unpolarify(tanh(p)) == tanh(u)
    assert unpolarify(gamma(p)) == gamma(u)
    assert unpolarify(erf(p)) == erf(u)
    assert unpolarify(uppergamma(x, p)) == uppergamma(x, p)

    assert unpolarify(uppergamma(sin(p), sin(p + exp_polar(0)))) == \
        uppergamma(sin(u), sin(u + 1))
    assert unpolarify(uppergamma(polar_lift(0), 2*exp_polar(0))) == \
        uppergamma(0, 2)

    assert unpolarify(Eq(p, 0)) == Eq(u, 0)
    assert unpolarify(Ne(p, 0)) == Ne(u, 0)
    assert unpolarify(polar_lift(x) > 0) == (x > 0)

    # Test bools
    assert unpolarify(True) is True


def test_sympyissue_4035():
    assert abs(x).expand(trig=True) == abs(x)
    assert sign(x).expand(trig=True) == sign(x)
    assert arg(x).expand(trig=True) == arg(x)


def test_sympyissue_3206():
    assert abs(abs(x)) == abs(x)


def test_sympyissue_4757():
    x = Symbol('x', extended_real=True)
    y = Symbol('y', imaginary=True)
    f = Function('f')
    assert arg(f(x)).diff(x).subs({f(x): 1 + I*x**2}).doit() == 2*x/(1 + x**4)
    assert arg(f(y)).diff(y).subs({f(y): I + y**2}).doit() == 2*y/(1 + y**4)


def test_periodic_argument():
    p = Symbol('p', positive=True)

    assert unbranched_argument(2 + I) == periodic_argument(2 + I, oo)
    assert unbranched_argument(1 + x) == periodic_argument(1 + x, oo)
    assert N_equals(unbranched_argument((1 + I)**2), pi/2)
    assert N_equals(unbranched_argument((1 - I)**2), -pi/2)
    assert N_equals(periodic_argument((1 + I)**2, 3*pi), pi/2)
    assert N_equals(periodic_argument((1 - I)**2, 3*pi), -pi/2)

    assert unbranched_argument(principal_branch(x, pi)) == \
        periodic_argument(x, pi)

    assert unbranched_argument(polar_lift(2 + I)) == unbranched_argument(2 + I)
    assert periodic_argument(polar_lift(2 + I), 2*pi) == \
        periodic_argument(2 + I, 2*pi)
    assert periodic_argument(polar_lift(2 + I), 3*pi) == \
        periodic_argument(2 + I, 3*pi)
    assert periodic_argument(polar_lift(2 + I), pi) == \
        periodic_argument(polar_lift(2 + I), pi)

    assert unbranched_argument(polar_lift(1 + I)) == pi/4
    assert periodic_argument(2*p, p) == periodic_argument(p, p)
    assert periodic_argument(pi*p, p) == periodic_argument(p, p)

    assert abs(polar_lift(1 + I)) == abs(1 + I)

    assert periodic_argument(x, pi).is_real is True
    assert periodic_argument(x, oo, evaluate=False).is_real is None

    a = Symbol('a', polar=True)
    assert periodic_argument(a, oo) == periodic_argument(a, oo, evaluate=False)


@pytest.mark.xfail
@pytest.mark.slow
def test_principal_branch_fail():
    # TODO XXX why does abs(x)._eval_evalf() not fall back to global evalf?
    assert N_equals(principal_branch((1 + I)**2, pi/2), 0)


def test_principal_branch():
    p = Symbol('p', positive=True)
    neg = Symbol('x', negative=True)

    assert principal_branch(polar_lift(x), p) == principal_branch(x, p)
    assert principal_branch(polar_lift(2 + I), p) == principal_branch(2 + I, p)
    assert principal_branch(2*x, p) == 2*principal_branch(x, p)
    assert principal_branch(1, pi) == exp_polar(0)
    assert principal_branch(-1, 2*pi) == exp_polar(I*pi)
    assert principal_branch(-1, pi) == exp_polar(0)
    assert principal_branch(exp_polar(3*pi*I)*x, 2*pi) == \
        principal_branch(exp_polar(I*pi)*x, 2*pi)
    assert principal_branch(neg*exp_polar(pi*I), 2*pi) == neg*exp_polar(-I*pi)

    assert N_equals(principal_branch((1 + I)**2, 2*pi), 2*I)
    assert N_equals(principal_branch((1 + I)**2, 3*pi), 2*I)
    assert N_equals(principal_branch((1 + I)**2, 1*pi), 2*I)

    # test argument sanitization
    assert isinstance(principal_branch(x, I), principal_branch)
    assert isinstance(principal_branch(x, -4), principal_branch)
    assert isinstance(principal_branch(x, -oo), principal_branch)
    assert isinstance(principal_branch(x, zoo), principal_branch)

    assert (principal_branch((4 + I)**2, 2*pi).evalf() ==
            principal_branch((4 + I)**2, 2*pi))

    assert principal_branch(exp_polar(-I*pi)*polar_lift(-1 - I),
                            2*pi) == sqrt(2)*exp_polar(I*pi/4)


@pytest.mark.xfail
def test_sympyissue_6167_6151():
    n = pi**1000
    i = int(n)
    assert sign(n - i) == 1
    assert abs(n - i) == n - i
    eps = pi**-1500
    big = pi**1000
    one = cos(x)**2 + sin(x)**2
    e = big*one - big + eps
    assert sign(simplify(e)) == 1
    for xi in (111, 11, 1, Rational(1, 10)):
        assert sign(e.subs({x: xi})) == 1


def test_sympyissue_11413():
    V = Matrix([[x], [y], [z]])
    U = V.normalized()
    r = sqrt(abs(x)**2 + abs(y)**2 + abs(z)**2)
    assert U == Matrix([[x/r], [y/r], [z/r]])
    assert U.norm().simplify() == 1
