import pytest

from diofant import (Add, Basic, E, Eq, Float, Function, GoldenRatio, I,
                     Integer, Integral, Lt, Matrix, MatrixSymbol, Mul, Number,
                     Piecewise, Rational, Sum, Symbol, acos, asin, atan,
                     besseli, besselj, besselsimp, binomial, cancel, cbrt,
                     combsimp, cos, cosh, cosine_transform, count_ops, diff,
                     erf, exp, exp_polar, expand, expand_multinomial,
                     expand_power_exp, factor, factorial, gamma, hyper,
                     hypersimp, integrate, ln, log, logcombine, nsimplify, oo,
                     pi, posify, root, separatevars, sign, signsimp, simplify,
                     sin, sinh, solve, sqrt, sqrtdenest, sstr, symbols, tan,
                     trigsimp, true, zoo)
from diofant.abc import (R, a, b, c, d, e, f, g, h, i, k, m, n, r, s, t, w, x,
                         y, z)
from diofant.core.mul import _keep_coeff
from diofant.simplify.simplify import clear_coefficients, nthroot


__all__ = ()


def test_factorial_simplify():
    # There are more tests in test_factorials.py. These are just to
    # ensure that simplify() calls factorial_simplify correctly
    assert simplify(factorial(x)/x) == factorial(x - 1)
    assert simplify(factorial(factorial(x))) == factorial(factorial(x))


def test_simplify_expr():
    A = Symbol('A')
    f = Function('f')

    assert all(simplify(tmp) == tmp for tmp in [I, E, oo, x, -x, -oo, -E, -I])

    e = 1/x + 1/y
    assert e != (x + y)/(x*y)
    assert simplify(e) == (x + y)/(x*y)

    e = A**2*s**4/(4*pi*k*m**3)
    assert simplify(e) == e

    e = (4 + 4*x - 2*(2 + 2*x))/(2 + 2*x)
    assert simplify(e) == 0

    e = (-4*x*y**2 - 2*y**3 - 2*x**2*y)/(x + y)**2
    assert simplify(e) == -2*y

    e = -x - y - (x + y)**(-1)*y**2 + (x + y)**(-1)*x**2
    assert simplify(e) == -2*y

    e = (x + x*y)/x
    assert simplify(e) == 1 + y

    e = (f(x) + y*f(x))/f(x)
    assert simplify(e) == 1 + y

    e = (2 * (1/n - cos(n * pi)/n))/pi
    assert simplify(e) == (-cos(pi*n) + 1)/(pi*n)*2

    e = integrate(1/(x**3 + 1), x).diff(x)
    assert simplify(e) == 1/(x**3 + 1)

    e = integrate(x/(x**2 + 3*x + 1), x).diff(x)
    assert simplify(e) == x/(x**2 + 3*x + 1)

    f = Symbol('f')
    A = Matrix([[2*k - m*w**2, -k], [-k, k - m*w**2]]).inv()
    assert simplify((A*Matrix([0, f]))[1]) == \
        -f*(2*k - m*w**2)/(k**2 - (k - m*w**2)*(2*k - m*w**2))

    f = -x + y/(z + t) + z*x/(z + t) + z*a/(z + t) + t*x/(z + t)
    assert simplify(f) == (y + a*z)/(z + t)

    A, B = symbols('A,B', commutative=False)

    assert simplify(A*B - B*A) == A*B - B*A
    assert simplify(A/(1 + y/x)) == x*A/(x + y)
    assert simplify(A*(1/x + 1/y)) == A/x + A/y  # (x + y)*A/(x*y)

    assert simplify(log(2) + log(3)) == log(6)
    assert simplify(log(2*x) - log(2)) == log(x)

    assert simplify(hyper([], [], x)) == exp(x)


def test_sympyissue_3557():
    f_1 = x*a + y*b + z*c - 1
    f_2 = x*d + y*e + z*f - 1
    f_3 = x*g + y*h + z*i - 1

    solutions = solve([f_1, f_2, f_3], x, y, z, simplify=False)

    assert simplify(solutions[0][y]) == \
        (a*i + c*d + f*g - a*f - c*g - d*i) / \
        (a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g)


def test_simplify_other():
    assert simplify(sin(x)**2 + cos(x)**2) == 1
    assert simplify(gamma(x + 1)/gamma(x)) == x
    assert simplify(sin(x)**2 + cos(x)**2 + factorial(x)/gamma(x)) == 1 + x
    assert simplify(
        Eq(sin(x)**2 + cos(x)**2, factorial(x)/gamma(x))) == Eq(x, 1)
    nc = symbols('nc', commutative=False)
    assert simplify(x + x*nc) == x*(1 + nc)
    # issue sympy/sympy#6123
    # f = exp(-I*(k*sqrt(t) + x/(2*sqrt(t)))**2)
    # ans = integrate(f, (k, -oo, oo), conds='none')
    ans = (I*(-pi*x*exp(-3*I*pi/4 + I*x**2/(4*t))*erf(x*exp(-3*I*pi/4) /
                                                      (2*sqrt(t)))/(2*sqrt(t)) + pi*x*exp(-3*I*pi/4 + I*x**2/(4*t)) /
              (2*sqrt(t)))*exp(-I*x**2/(4*t))/(sqrt(pi)*x) - I*sqrt(pi) *
           (-erf(x*exp(I*pi/4)/(2*sqrt(t))) + 1)*exp(I*pi/4)/(2*sqrt(t)))
    assert simplify(expand_power_exp(ans)) == -(-1)**Rational(3, 4)*sqrt(pi)/sqrt(t)
    # issue sympy/sympy#6370
    assert simplify(2**(2 + x)/4) == 2**x


def test_simplify_complex():
    cosAsExp = cos(x)._eval_rewrite_as_exp(x)
    tanAsExp = tan(x)._eval_rewrite_as_exp(x)
    assert simplify(cosAsExp*tanAsExp).expand() == (
        sin(x))._eval_rewrite_as_exp(x).expand()  # issue sympy/sympy#4341


def test_simplify_ratio():
    # roots of x**3-3*x+5
    roots = [(Rational(1, 2) - sqrt(3)*I/2) *
             cbrt(sqrt(21)/2 + Rational(5, 2)) +
             1/((Rational(1, 2) - sqrt(3)*I/2)*cbrt(sqrt(21)/2 +
                                                    Rational(5, 2))),
             1/((Rational(1, 2) + sqrt(3)*I/2) *
                cbrt(sqrt(21)/2 + Rational(5, 2))) +
             (Rational(1, 2) + sqrt(3)*I/2)*cbrt(sqrt(21)/2 +
                                                 Rational(5, 2)),
             -cbrt(sqrt(21)/2 + Rational(5, 2)) -
             1/cbrt(sqrt(21)/2 + Rational(5, 2))]

    for root in roots:
        assert count_ops(simplify(root, ratio=1)) <= count_ops(root)
        # If ratio=oo, simplify() is always applied:
        assert simplify(root, ratio=oo) is not root


def test_simplify_measure():
    def measure1(expr):
        return len(str(expr))

    def measure2(expr):
        return -count_ops(expr)  # Return the most complicated result

    expr = (x + 1)/(x + sin(x)**2 + cos(x)**2)
    assert measure1(simplify(expr, measure=measure1)) <= measure1(expr)
    assert measure2(simplify(expr, measure=measure2)) <= measure2(expr)

    expr2 = Eq(sin(x)**2 + cos(x)**2, 1)
    assert measure1(simplify(expr2, measure=measure1)) <= measure1(expr2)
    assert measure2(simplify(expr2, measure=measure2)) <= measure2(expr2)


def test_sympyissue_4407():
    assert simplify(exp(-Rational(1, 2)) + exp(-Rational(3, 2))) == \
        (1 + E)*exp(-Rational(3, 2))


def test_sympyissue_5652():
    assert simplify(E + exp(-E)) == exp(-E) + E
    n = symbols('n', commutative=False)
    assert simplify(n + n**(-n)) == n + n**(-n)


def test_simplify_fail1():
    e = (x + y)**2/(-4*x*y**2 - 2*y**3 - 2*x**2*y)
    assert simplify(e) == 1 / (-2*y)


def test_nthroot():
    assert nthroot(90 + 34*sqrt(7), 3) == sqrt(7) + 3
    q = 1 + sqrt(2) - 2*sqrt(3) + sqrt(6) + sqrt(7)
    assert nthroot(expand_multinomial(q**3), 3) == q
    assert nthroot(41 + 29*sqrt(2), 5) == 1 + sqrt(2)
    assert nthroot(-41 - 29*sqrt(2), 5) == -1 - sqrt(2)
    expr = 1320*sqrt(10) + 4216 + 2576*sqrt(6) + 1640*sqrt(15)
    assert nthroot(expr, 5) == 1 + sqrt(6) + sqrt(15)
    q = 1 + sqrt(2) + sqrt(3) + sqrt(5)
    assert expand_multinomial(nthroot(expand_multinomial(q**5), 5)) == q
    q = 1 + sqrt(2) + 7*sqrt(6) + 2*sqrt(10)
    assert nthroot(expand_multinomial(q**5), 5, 8) == q
    q = 1 + sqrt(2) - 2*sqrt(3) + 1171*sqrt(6)
    assert nthroot(expand_multinomial(q**3), 3) == q
    assert nthroot(expand_multinomial(q**6), 6) == q


@pytest.mark.slow
def test_nthroot1():
    q = 1 + sqrt(2) + sqrt(3) + Rational(1, 10)**20
    p = expand_multinomial(q**5)
    assert nthroot(p, 5) == q
    q = 1 + sqrt(2) + sqrt(3) + Rational(1, 10)**30
    p = expand_multinomial(q**5)
    assert nthroot(p, 5) == q


def test_separatevars():
    n = Symbol('n')
    assert separatevars(2*n*x*z + 2*x*y*z) == 2*x*z*(n + y)
    assert separatevars(x*z + x*y*z) == x*z*(1 + y)
    assert separatevars(pi*x*z + pi*x*y*z) == pi*x*z*(1 + y)
    assert separatevars(x*y**2*sin(x) + x*sin(x)*sin(y)) == \
        x*(sin(y) + y**2)*sin(x)
    assert separatevars(x*exp(x + y) + x*exp(x)) == x*(1 + exp(y))*exp(x)
    assert separatevars((x*(y + 1))**z).is_Pow  # != x**z*(1 + y)**z
    assert separatevars(1 + x + y + x*y) == (x + 1)*(y + 1)
    assert separatevars(y/pi*exp(-(z - x)/cos(n))) == \
        y*exp(x/cos(n))*exp(-z/cos(n))/pi
    assert separatevars((x + y)*(x - y) + y**2 + 2*x + 1) == (x + 1)**2
    # issue sympy/sympy#4858
    p = Symbol('p', positive=True)
    assert separatevars(sqrt(p**2 + x*p**2)) == p*sqrt(1 + x)
    assert separatevars(sqrt(y*(p**2 + x*p**2))) == p*sqrt(y*(1 + x))
    assert separatevars(sqrt(y*(p**2 + x*p**2)), force=True) == \
        p*sqrt(y)*sqrt(1 + x)
    # issue sympy/sympy#4865
    assert separatevars(sqrt(x*y)).is_Pow
    assert separatevars(sqrt(x*y), force=True) == sqrt(x)*sqrt(y)
    # issue sympy/sympy#4957
    # any type sequence for symbols is fine
    assert separatevars(((2*x + 2)*y), dict=True, symbols=()) == \
        {'coeff': 1, x: 2*x + 2, y: y}
    # separable
    assert separatevars(((2*x + 2)*y), dict=True, symbols=[x]) == \
        {'coeff': y, x: 2*x + 2}
    assert separatevars(((2*x + 2)*y), dict=True, symbols=[]) == \
        {'coeff': 1, x: 2*x + 2, y: y}
    assert separatevars(((2*x + 2)*y), dict=True) == \
        {'coeff': 1, x: 2*x + 2, y: y}
    assert separatevars(((2*x + 2)*y), dict=True, symbols=None) == \
        {'coeff': y*(2*x + 2)}
    # not separable
    assert separatevars(3, dict=True) is None
    assert separatevars(2*x + y, dict=True, symbols=()) is None
    assert separatevars(2*x + y, dict=True) is None
    assert separatevars(2*x + y, dict=True, symbols=None) == {'coeff': 2*x + y}
    # issue sympy/sympy#4808
    n, m = symbols('n,m', commutative=False)
    assert separatevars(m + n*m) == (1 + n)*m
    assert separatevars(x + x*n) == x*(1 + n)
    # issue sympy/sympy#4910
    f = Function('f')
    assert separatevars(f(x) + x*f(x)) == f(x) + x*f(x)
    # a noncommutable object present
    eq = x*(1 + hyper((), (), y*z))
    assert separatevars(eq) == eq
    pytest.raises(ValueError, lambda: separatevars(3, [x**y], dict=True))


def test_separatevars_advanced_factor():
    x, y = symbols('x,y')
    assert (separatevars(1 + log(x)*log(y) + log(x) + log(y)) ==
            (log(x) + 1)*(log(y) + 1))
    assert (separatevars(1 + x - log(z) - x*log(z) - exp(y)*log(z) -
                         x*exp(y)*log(z) + x*exp(y) +
                         exp(y)) == -((exp(y) + 1) *
                                      (x + 1)*(log(z) - 1)))
    x, y = symbols('x,y', positive=True)
    assert (separatevars(1 + log(x**log(y)) + log(x*y)) ==
            (log(x) + 1)*(log(y) + 1))


def test_hypersimp():
    n, k = symbols('n,k', integer=True)

    assert hypersimp(factorial(k), k) == k + 1
    assert hypersimp(factorial(k**2), k) is None

    assert hypersimp(1/factorial(k), k) == 1/(k + 1)

    assert hypersimp(2**k/factorial(k)**2, k) == 2/(k + 1)**2

    assert hypersimp(binomial(n, k), k) == (n - k)/(k + 1)
    assert hypersimp(binomial(n + 1, k), k) == (n - k + 1)/(k + 1)

    term = (4*k + 1)*factorial(k)/factorial(2*k + 1)
    assert hypersimp(term, k) == ((4*k + 5)/(3 + 14*k + 8*k**2))/2

    term = 1/((2*k - 1)*factorial(2*k + 1))
    assert hypersimp(term, k) == (k - Rational(1, 2))/((k + 1)*(2*k + 1)*(2*k + 3))

    term = binomial(n, k)*(-1)**k/factorial(k)
    assert hypersimp(term, k) == (k - n)/(k + 1)**2

    assert hypersimp(2**(I*k) * 2**k, k) == 2**(1 + I)


def test_nsimplify():
    assert nsimplify(0) == 0
    assert nsimplify(-1) == -1
    assert nsimplify(1) == 1
    assert nsimplify(1 + x) == 1 + x
    assert nsimplify(2.7) == Rational(27, 10)
    assert nsimplify(1 - GoldenRatio) == (1 - sqrt(5))/2
    assert nsimplify((1 + sqrt(5))/4, [GoldenRatio]) == GoldenRatio/2
    assert nsimplify(2/GoldenRatio, [GoldenRatio]) == 2*GoldenRatio - 2
    assert nsimplify(exp(5*pi*I/3, evaluate=False)) == Rational(1, 2) - sqrt(3)*I/2
    assert nsimplify(sin(3*pi/5, evaluate=False)) == sqrt(sqrt(5)/8 +
                                                          Rational(5, 8))
    assert nsimplify(sqrt(atan('1', evaluate=False))*(2 + I), [pi]) == \
        sqrt(pi) + sqrt(pi)/2*I
    assert nsimplify(2 + exp(2*atan('1/4')*I)) == Rational(49, 17) + 8*I/17
    assert nsimplify(pi, tolerance=0.01) == Rational(22, 7)
    assert nsimplify(pi, tolerance=0.001) == Rational(355, 113)
    assert nsimplify(0.33333, tolerance=1e-4) == Rational(1, 3)
    assert nsimplify(2.0**(1/3.), tolerance=0.001) == Rational(635, 504)
    assert nsimplify(2.0**(1/3.), tolerance=0.001, full=True) == cbrt(2)
    assert nsimplify(x + .5, rational=True) == Rational(1, 2) + x
    assert nsimplify(1/.3 + x, rational=True) == Rational(10, 3) + x
    assert nsimplify(log(3).evalf(), rational=True) == Rational(109861228866811,
                                                                100000000000000)
    assert nsimplify(Float(0.272198261287950), [pi, log(2)]) == pi*log(2)/8
    assert nsimplify(Float(0.272198261287950).evalf(3), [pi, log(2)]) == \
        -pi/4 - log(2) + Rational(7, 4)
    assert nsimplify(x/7.0) == x/7
    assert nsimplify(pi/1e2) == pi/100
    assert nsimplify(pi/1e2, rational=False) == pi/100.0
    assert nsimplify(pi/1e-7) == 10000000*pi
    assert not nsimplify(
        factor(-3.0*z**2*(z**2)**(-2.5) + 3*(z**2)**(-1.5))).atoms(Float)
    e = x**0.0
    assert e.is_Pow
    assert nsimplify(x**0.0) == 1
    assert nsimplify(3.333333, tolerance=0.1, rational=True) == Rational(10, 3)
    assert nsimplify(3.333333, tolerance=0.01, rational=True) == Rational(10, 3)
    assert nsimplify(3.666666, tolerance=0.1, rational=True) == Rational(11, 3)
    assert nsimplify(3.666666, tolerance=0.01, rational=True) == Rational(11, 3)
    assert nsimplify(33, tolerance=10, rational=True) == 33
    assert nsimplify(33.33, tolerance=10, rational=True) == 30
    assert nsimplify(37.76, tolerance=10, rational=True) == 40
    assert nsimplify(-203.1) == -Rational(2031, 10)
    assert nsimplify(+.2, tolerance=0) == Rational(+1, 5)
    assert nsimplify(-.2, tolerance=0) == Rational(-1, 5)
    assert nsimplify(.2222, tolerance=0) == Rational(1111, 5000)
    assert nsimplify(-.2222, tolerance=0) == -Rational(1111, 5000)
    # issue sympy/sympy#7211, PR sympy/sympy#4112
    assert nsimplify(Float(2e-8)) == Rational(1, 50000000)
    # issue sympy/sympy#7322 direct test
    assert nsimplify(1e-42, rational=True) != 0
    # issue sympy/sympy#10336
    inf = Float('inf')
    infs = (-oo, oo, inf, -inf)
    for i in infs:
        ans = sign(i)*oo
        assert nsimplify(i) == ans
        assert nsimplify(i + x) == x + ans

    assert nsimplify(Sum(1/n**2, (n, 1, oo)), [pi]) == pi**2/6


def test_sympyissue_9448():
    expr = (1/(1 - (-1)**Rational(2, 3) - cbrt(-1)) +
            1/(1 + (-1)**Rational(2, 3) + cbrt(-1)))
    assert nsimplify(expr) == Rational(1, 2)


def test_extract_minus_sign():
    assert simplify(-x/-y) == x/y
    assert simplify(-x/y) == -x/y
    assert simplify(x/y) == x/y
    assert simplify(x/-y) == -x/y
    assert simplify(-x/0) == zoo*x
    assert simplify(Rational(-5, 0)) == zoo
    assert simplify(-a*x/(-y - b)) == a*x/(b + y)


def test_diff():
    f = Function('f')
    g = Function('g')
    assert simplify(g(x).diff(x)*f(x).diff(x) - f(x).diff(x)*g(x).diff(x)) == 0
    assert simplify(2*f(x)*f(x).diff(x) - diff(f(x)**2, x)) == 0
    assert simplify(diff(1/f(x), x) + f(x).diff(x)/f(x)**2) == 0
    assert simplify(f(x).diff(x, y) - f(x).diff(y, x)) == 0


def test_logcombine_1():
    z, w = symbols('z,w', positive=True)
    b = Symbol('b', real=True)
    assert logcombine(log(x) + 2*log(y)) == log(x) + 2*log(y)
    assert logcombine(log(x) + 2*log(y), force=True) == log(x*y**2)
    assert logcombine(a*log(w) + log(z)) == a*log(w) + log(z)
    assert logcombine(b*log(z) + b*log(x)) == log(z**b) + b*log(x)
    assert logcombine(b*log(z) - log(w)) == log(z**b/w)
    assert logcombine(log(x)*log(z)) == log(x)*log(z)
    assert logcombine(log(w)*log(x)) == log(w)*log(x)
    assert logcombine(cos(-2*log(z) + b*log(w))) in [cos(log(w**b/z**2)),
                                                     cos(log(z**2/w**b))]
    assert logcombine(log(log(x) - log(y)) - log(z), force=True) == \
        log(log(x/y)/z)
    assert logcombine((2 + I)*log(x), force=True) == (2 + I)*log(x)
    assert logcombine((x**2 + log(x) - log(y))/(x*y), force=True) == \
        (x**2 + log(x/y))/(x*y)
    # the following could also give log(z*x**log(y**2)), what we
    # are testing is that a canonical result is obtained
    assert logcombine(log(x)*2*log(y) + log(z), force=True) == \
        log(z*y**log(x**2))
    assert logcombine((x*y + sqrt(x**4 + y**4) + log(x) -
                       log(y))/(pi*x**Rational(2, 3) *
                                sqrt(y)**3), force=True) == (
        x*y + sqrt(x**4 + y**4) + log(x/y))/(pi*x**Rational(2, 3) *
                                             y**Rational(3, 2))
    assert logcombine(gamma(-log(x/y))*acos(-log(x/y)), force=True) == \
        acos(-log(x/y))*gamma(-log(x/y))

    assert logcombine(2*log(z)*log(w)*log(x) + log(z) + log(w)) == \
        log(z**log(w**2))*log(x) + log(w*z)
    assert logcombine(3*log(w) + 3*log(z)) == log(w**3*z**3)
    assert logcombine(x*(y + 1) + log(2) + log(3)) == x*(y + 1) + log(6)
    assert logcombine((x + y)*log(w) + (-x - y)*log(3)) == (x + y)*log(w/3)


def test_logcombine_complex_coeff():
    i = Integral((sin(x**2) + cos(x**3))/x, x)
    assert logcombine(i, force=True) == i
    assert logcombine(i + 2*log(x), force=True) == i + log(x**2)


def test_posify():
    assert sstr(posify(
        x +
        Symbol('p', positive=True) +
        Symbol('n', negative=True))) == '(n + p + _x, {_x: x})'

    eq, rep = posify(1/x)
    assert log(eq).expand().subs(rep) == -log(x)
    assert sstr(posify([x, 1 + x])) == '([_x, _x + 1], {_x: x})'

    p = symbols('p', positive=True)
    n = symbols('n', negative=True)
    orig = [x, n, p]
    modified, reps = posify(orig)
    assert sstr(modified) == '[_x, n, p]'
    assert [w.subs(reps) for w in modified] == orig

    assert sstr(Integral(posify(1/x + y)[0], (y, 1, 3)).expand()) == \
        'Integral(1/_x, (y, 1, 3)) + Integral(_y, (y, 1, 3))'
    assert sstr(Sum(posify(1/x**n)[0], (n, 1, 3)).expand()) == \
        'Sum(_x**(-n), (n, 1, 3))'


def test_sympyissue_4194():
    # simplify should call cancel
    f = Function('f')
    assert simplify((4*x + 6*f(y))/(2*x + 3*f(y))) == 2


def test_as_content_primitive():
    # although the _as_content_primitive methods do not alter the underlying structure,
    # the as_content_primitive function will touch up the expression and join
    # bases that would otherwise have not been joined.
    assert (x*(2 + 2*x)*(3*x + 3)**2).as_content_primitive() == \
        (18, x*(x + 1)**3)
    assert (2 + 2*x + 2*y*(3 + 3*y)).as_content_primitive() == \
        (2, x + 3*y*(y + 1) + 1)
    assert ((2 + 6*x)**2).as_content_primitive() == (4, (3*x + 1)**2)
    assert ((2 + 6*x)**(2*y)).as_content_primitive() == \
        (1, (_keep_coeff(Integer(2), (3*x + 1)))**(2*y))
    assert (5 + 10*x + 2*y*(3 + 3*y)).as_content_primitive() == \
        (1, 10*x + 6*y*(y + 1) + 5)
    assert (5*(x*(1 + y)) + 2*x*(3 + 3*y)).as_content_primitive() == \
        (11, x*(y + 1))
    assert ((5*(x*(1 + y)) + 2*x*(3 + 3*y))**2).as_content_primitive() == \
        (121, x**2*(y + 1)**2)
    assert (y**2).as_content_primitive() == (1, y**2)
    assert oo.as_content_primitive() == (1, oo)
    eq = x**(2 + y)
    assert (eq).as_content_primitive() == (1, eq)
    assert (Rational(1, 2)**(2 + x)).as_content_primitive() == (Rational(1, 4), 2**-x)
    assert (Rational(-1, 2)**(2 + x)).as_content_primitive() == \
           (Rational(1, 4), Rational(-1, 2)**x)
    assert (Rational(-1, 2)**(2 + x)).as_content_primitive() == \
           (Rational(1, 4), Rational(-1, 2)**x)
    assert (4**((1 + y)/2)).as_content_primitive() == (2, 4**(y/2))
    assert (3**((1 + y)/2)).as_content_primitive() == \
           (1, 3**(Mul(Rational(1, 2), 1 + y, evaluate=False)))
    assert (5**Rational(3, 4)).as_content_primitive() == (1, 5**Rational(3, 4))
    assert (5**Rational(7, 4)).as_content_primitive() == (5, 5**Rational(3, 4))
    assert Add(5*z/7, 0.5*x, 3*y/2, evaluate=False).as_content_primitive() == \
              (Rational(1, 14), 7.0*x + 21*y + 10*z)
    assert (2**Rational(3, 4) + root(2, 4)*sqrt(3)).as_content_primitive(radical=True) == \
           (1, root(2, 4)*(sqrt(2) + sqrt(3)))


def test_signsimp():
    e = x*(-x + 1) + x*(x - 1)
    assert signsimp(Eq(e, 0)) is true
    assert abs(x - 1) == abs(1 - x)


def test_besselsimp():
    assert besselsimp(exp(-I*pi*y/2)*besseli(y, z*exp_polar(I*pi/2))) == \
        besselj(y, z)
    assert besselsimp(exp(-I*pi*a/2)*besseli(a, 2*sqrt(x)*exp_polar(I*pi/2))) == \
        besselj(a, 2*sqrt(x))
    assert besselsimp(sqrt(2)*sqrt(pi)*root(x, 4)*exp(I*pi/4)*exp(-I*pi*a/2) *
                      besseli(-Rational(1, 2), sqrt(x)*exp_polar(I*pi/2)) *
                      besseli(a, sqrt(x)*exp_polar(I*pi/2))/2) == \
        besselj(a, sqrt(x)) * cos(sqrt(x))
    assert besselsimp(besseli(Rational(-1, 2), z)) == \
        sqrt(2)*cosh(z)/(sqrt(pi)*sqrt(z))
    assert besselsimp(besseli(a, z*exp_polar(-I*pi/2))) == \
        exp(-I*pi*a/2)*besselj(a, z)
    assert cosine_transform(1/t*sin(a/t), t, y) == \
        sqrt(2)*sqrt(pi)*besselj(0, 2*sqrt(a)*sqrt(y))/2


def test_Piecewise():
    e1 = x*(x + y) - y*(x + y)
    e2 = sin(x)**2 + cos(x)**2
    e3 = expand((x + y)*y/x)
    s1 = simplify(e1)
    s2 = simplify(e2)
    s3 = simplify(e3)
    assert simplify(Piecewise((e1, x < e2), (e3, True))) == \
        Piecewise((s1, x < s2), (s3, True))


def test_polymorphism():
    class A(Basic):
        def _eval_simplify(self, **kwargs):
            return 1

    a = A(5, 2)
    assert simplify(a) == 1


def test_sympyissue_6811():
    eq = (x + 2*y)*(2*x + 2)
    assert simplify(eq) == (x + 1)*(x + 2*y)*2
    # reject the 2-arg Mul -- these are a headache for test writing
    assert simplify(eq.expand()) == 2*x**2 + 4*x*y + 2*x + 4*y


@pytest.mark.xfail
def test_sympyissue_6811_fail():
    eq = 4*(-19*sin(x)*y + 5*sin(3*x)*y + 15*cos(2*x)*z - 21*z)*t/(9*cos(x) - 5*cos(3*x))
    assert trigsimp(eq) == -2*(2*sin(x)*y + 3*z)*t/cos(x)


def test_sympyissue_6920():
    e = [cos(x) + I*sin(x), cos(x) - I*sin(x),
         cosh(x) - sinh(x), cosh(x) + sinh(x)]
    ok = [exp(I*x), exp(-I*x), exp(-x), exp(x)]
    # wrap in f to show that the change happens wherever ei occurs
    f = Function('f')
    assert [simplify(f(ei)).args[0] for ei in e] == ok


def test_sympyissue_7001():
    assert (simplify(-(r*Piecewise((4*pi/3, r <= R),
                                   (-8*pi*R**3/(3*r**3), True)) +
                       2*Piecewise((4*pi*r/3, r <= R),
                                   (4*pi*R**3/(3*r**2), True)))/(4*pi*r)) ==
            Piecewise((-1, r <= R), (0, True)))


def test_inequality_no_auto_simplify():
    # no simplify on creation but can be simplified
    lhs = cos(x)**2 + sin(x)**2
    rhs = 2
    e = Lt(lhs, rhs)
    assert e == Lt(lhs, rhs, evaluate=False)
    assert simplify(e)


def test_sympyissue_9398():
    assert cancel(1e-14) != 0
    assert cancel(1e-14*I) != 0

    assert simplify(1e-14) != 0
    assert simplify(1e-14*I) != 0

    assert (I*Number(1.)*Number(10)**Number(-14)).simplify() != 0

    assert cancel(1e-20) != 0
    assert cancel(1e-20*I) != 0

    assert simplify(1e-20) != 0
    assert simplify(1e-20*I) != 0

    assert cancel(1e-100) != 0
    assert cancel(1e-100*I) != 0

    assert simplify(1e-100) != 0
    assert simplify(1e-100*I) != 0

    f = Float('1e-1000', 15)
    assert cancel(f) != 0
    assert cancel(f*I) != 0

    assert simplify(f) != 0
    assert simplify(f*I) != 0


def test_sympyissue_6249():
    A = MatrixSymbol('A', 3, 3)
    B = MatrixSymbol('B', 3, 3)
    p = A*B - B*A
    assert cancel(p) == p
    assert combsimp(p) == p
    assert factor(p) == p
    assert separatevars(p) == p
    assert sqrtdenest(p) == p

    M = MatrixSymbol('M', 2, 1)
    assert simplify(M[0]/2) == M[0]/2


def test_clear_coefficients():
    assert clear_coefficients(4*y*(6*x + 3)) == (y*(2*x + 1), 0)
    assert clear_coefficients(4*y*(6*x + 3) - 2) == (y*(2*x + 1), Rational(1, 6))
    assert clear_coefficients(4*y*(6*x + 3) - 2, x) == (y*(2*x + 1), x/12 + Rational(1, 6))
    assert clear_coefficients(sqrt(2) - 2) == (sqrt(2), 2)
    assert clear_coefficients(4*sqrt(2) - 2) == (sqrt(2), Rational(1, 2))
    assert clear_coefficients(Integer(3), x) == (0, x - 3)
    assert clear_coefficients(oo, x) == (oo, x)
    assert clear_coefficients(-pi, x) == (pi, -x)
    assert clear_coefficients(2 - pi/3, x) == (pi, -3*x + 6)


def test_sympyissue_9296():
    q = symbols('q_1:5')
    dq = symbols('dq_1:5')
    a = (dq[0]*(0.1*(0.01*sin(q[2]) + 0.01*sin(q[1] + q[2]))*cos(q[1] + q[2]) +
                0.01*(0.1*sin(q[2]) + 0.1*sin(q[1] + q[2]))*cos(q[1] + q[2]) +
                0.05*(-0.05*cos(q[1]) - 0.025)*sin(q[1]) +
                0.01*(-0.1*cos(q[2]) - 0.1*cos(q[1] + q[2]))*sin(q[1] + q[2]) +
                0.1*(-0.01*cos(q[2]) - 0.01*cos(q[1] + q[2]))*sin(q[1] + q[2]) +
                0.0025*sin(q[1])*cos(q[1]) - 0.00125*sin(q[1])))
    b = dq[1]*(-0.0025*sin(q[1]) + 0.002*sin(q[2])*cos(q[1] + q[2]) -
               0.002*sin(q[1] + q[2])*cos(q[2]))
    r = simplify(a + b).replace(lambda x: x.is_Float and abs(x) < 1e-15,
                                lambda x: 0)
    assert r == (-Float('0.0045000000000000005', dps=15)*dq[0]*sin(q[1]) -
                 Float('0.0045000000000000005', dps=15)*dq[1]*sin(q[1]))


def test_sympyissue_9630():
    psi = (-0.999999972295856*sin(13.3579685223169*x/y) + 1.0*cos(13.3579685223169*x/y) +
           0.999999972295856*sinh(13.3579685223169*x/y) - 1.0*cosh(13.3579685223169*x/y))
    assert simplify(psi) == psi


def test_sympyissue_12792():
    expr = (0.25*y*sin(x)/(0.25*cos(x)**2 - 1.0*cos(x) + 1)
            + (z - 0.5)*(-0.25*y*sin(x)*cos(x)**2
                         / (0.0625*cos(x)**4 - 0.5*cos(x)**3 + 1.5*cos(x)**2 - 2.0*cos(x) + 1)
                         + 0.5*y*sin(x)*cos(x)/(0.0625*cos(x)**4 - 0.5*cos(x)**3 + 1.5*cos(x)**2
                                                - 2.0*cos(x) + 1) + 0.25*cos(x)**3
                         / (0.0625*cos(x)**4 - 0.5*cos(x)**3 + 1.5*cos(x)**2 - 2.0*cos(x) + 1)
                         - 1.0*cos(x)**2/(0.0625*cos(x)**4 - 0.5*cos(x)**3 + 1.5*cos(x)**2 - 2.0*cos(x) + 1)
                         + 1.0*cos(x)/(0.0625*cos(x)**4 - 0.5*cos(x)**3 + 1.5*cos(x)**2 - 2.0*cos(x) + 1)
                         - 1/(0.25*cos(x)**2 - 1.0*cos(x) + 1)) - 0.25*cos(x)
            / (0.25*cos(x)**2 - 1.0*cos(x) + 1) + 0.5/(0.25*cos(x)**2 - 1.0*cos(x) + 1))

    expr_simp = simplify(expr)
    assert expr_simp.equals(expr)


def test_sympyissue_12506():
    expr = 1.0 * cos(x) + 2.0 * cos(asin(0.5 * sin(x)))
    expr = expr.diff((x, 2))
    expr_simp = expr.simplify()
    assert expr_simp.equals(expr)


def test_sympyissue_13115():
    q_1, q_2, q_3 = symbols('q_1:4')
    Mq = Matrix([[(1.0*cos(q_2) + 0.5*cos(q_2 + q_3))**2*sin(q_1)**2 +
                  (1.0*cos(q_2) + 0.5*cos(q_2 + q_3))**2*cos(q_1)**2 +
                  0.25*sin(q_1)**2*cos(q_2)**2 + 0.25*cos(q_1)**2*cos(q_2)**2, 0, 0],
                 [0, (-1.0*sin(q_2) - 0.5*sin(q_2 + q_3))**2*sin(q_1)**2 +
                  (-1.0*sin(q_2) - 0.5*sin(q_2 + q_3))**2*cos(q_1)**2 +
                  (-1.0*cos(q_2) - 0.5*cos(q_2 + q_3))**2 +
                  0.25*sin(q_1)**2*sin(q_2)**2 + 0.25*sin(q_2)**2*cos(q_1)**2 +
                  0.25*cos(q_2)**2, -0.5*(-1.0*sin(q_2) - 0.5*sin(q_2 + q_3))*sin(q_1)**2*sin(q_2 + q_3) -
                  0.5*(-1.0*sin(q_2) - 0.5*sin(q_2 + q_3))*sin(q_2 + q_3)*cos(q_1)**2 -
                  0.5*(-1.0*cos(q_2) - 0.5*cos(q_2 + q_3))*cos(q_2 + q_3)],
                 [0, -0.5*(-1.0*sin(q_2) - 0.5*sin(q_2 + q_3))*sin(q_1)**2*sin(q_2 + q_3) -
                  0.5*(-1.0*sin(q_2) - 0.5*sin(q_2 + q_3))*sin(q_2 + q_3)*cos(q_1)**2 -
                  0.5*(-1.0*cos(q_2) - 0.5*cos(q_2 + q_3))*cos(q_2 + q_3),
                  0.25*sin(q_1)**2*sin(q_2 + q_3)**2 + 0.25*sin(q_2 + q_3)**2*cos(q_1)**2 +
                  0.25*cos(q_2 + q_3)**2]])

    Mqs = simplify(Mq)
    assert Mqs.subs({q_1: 0, q_2: 0, q_3: 0}) == Matrix([[2.5, 0, 0],
                                                         [0, 2.5, 0.75],
                                                         [0, 0.75, 0.25]])


@pytest.mark.xfail
def test_simplify_algebraic_numbers():
    e = (3 + 4*I)**Rational(3, 2)
    assert simplify(e) == 2 + 11*I  # issue sympy/sympy#4401


@pytest.mark.timeout(20)
def test_sympyissue_21641():
    assert (simplify(65712362363534280139543*ln(Rational(49, 50)) /
                     2441406250000000000) ==
            -65712362363534280139543*log(50)/2441406250000000000 +
            65712362363534280139543*log(49)/2441406250000000000)
