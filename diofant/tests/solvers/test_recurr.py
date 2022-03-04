import pytest

from diofant import (Eq, Function, I, Lambda, binomial, cos, factorial, gamma,
                     pi, rf, rsolve, sin, sqrt, symbols)
from diofant.abc import a, b
from diofant.solvers.recurr import rsolve_hyper, rsolve_poly, rsolve_ratio


__all__ = ()

f, g = symbols('f,g', cls=Function)
n, k = symbols('n,k', integer=True)
C0, C1, C2 = symbols('C:3')


def test_poly():
    assert rsolve_poly([-1, -1, 1], 0, n) == (0, [])
    assert rsolve_poly([-1, -1, 1], 1, n) == (-1, [])
    assert rsolve_poly([-n**2, n, -1, 1], 1, n) is None

    assert rsolve(-f(n) + (n + 1)*f(n + 1) -
                  n) == [{f: Lambda(n, C0/factorial(n) + 1)}]
    assert rsolve(-(n + 1)*f(n) + n*f(n + 1) -
                  1) == [{f: Lambda(n, C0*n - 1)}]
    assert rsolve(-(4*n + 2)*f(n) + f(n + 1) - 4*n -
                  1) == [{f: Lambda(n, 4**n*C0*gamma(n + 1/2)/sqrt(pi) - 1)}]

    assert (rsolve(-f(n) + f(n + 1) - n**5 - n**3) ==
            [{f: Lambda(n, C0 + n**2*(n**4 - 3*n**3 + 4*n**2 - 3*n + 1)/6)}])

    assert rsolve_poly([1, 1], sqrt(n), n) is None

    assert rsolve_poly([-2, -1, 1],
                       -2*n**4 - (n + 1)**4 + (n + 2)**4, n) == (n**4, [])
    assert rsolve(n*f(n) - f(n + 1) - n**3 +
                  (n + 1)**2) == [{f: Lambda(n, (C0*factorial(n) + n**3)/n)}]


def test_ratio():
    assert rsolve_ratio([-2*n**3 + n**2 + 2*n - 1, 2*n**3 + n**2 - 6*n,
                         -2*n**3 - 11*n**2 - 18*n - 9, 2*n**3 + 13*n**2 +
                         22*n + 8], 0, n) == (C2*(2*n - 3)/(n**2 - 1)/2, [C2])
    assert rsolve_ratio([1, 1], sqrt(n), n) is None
    assert rsolve_ratio([-n**3, n + 1], n, n) is None


def test_hyper():
    assert rsolve((n**2 - 2)*f(n) - (2*n + 1)*f(n + 1) +
                  f(n + 2)) == [{f: Lambda(n, C0*rf(-sqrt(2), n) +
                                           C1*rf(+sqrt(2), n))}]

    assert rsolve((n**2 - k)*f(n) - (2*n + 1)*f(n + 1) +
                  f(n + 2)) == [{f: Lambda(n, C1*rf(sqrt(k), n) +
                                           C0*rf(-sqrt(k), n))}]

    assert rsolve(2*n*(n + 1)*f(n) - (n**2 + 3*n - 2)*f(n + 1) +
                  (n - 1)*f(n + 2)) == [{f: Lambda(n, C1*factorial(n) +
                                                   C0*2**n)}]

    assert rsolve_hyper([n + 2, -(2*n + 3)*(17*n**2 + 51*n + 39), n + 1],
                        0, n) == (0, [])

    assert rsolve_hyper([-n - 1, -1, 1], 0, n) == (0, [])

    assert rsolve(-n - f(n) + f(n + 1)) == [{f: Lambda(n, C0 + n*(n - 1)/2)}]
    assert rsolve(-1 - n - f(n) +
                  f(n + 1)) == [{f: Lambda(n, C0 + n*(n + 1)/2)}]
    assert rsolve(-3*(n + n**2) - f(n) +
                  f(n + 1)) == [{f: Lambda(n, C0 + n**3 - n)}]

    assert rsolve(-n - factorial(n) - f(n) + f(n + 1)) is None

    assert rsolve(-a*f(n) + f(n + 1)) == [{f: Lambda(n, C0*a**n)}]
    assert rsolve(-a*f(n) +
                  f(n + 2)) == [{f: Lambda(n, a**(n/2)*((-1)**n*C1 + C0))}]

    assert (rsolve(f(n) + f(n + 1) + f(n + 2)) ==
            [{f: Lambda(n, 2**-n*(C0*(-1 - sqrt(3)*I)**n +
                                  C1*(-1 + sqrt(3)*I)**n))}])

    assert rsolve_hyper([1, -2*n/a - 2/a, 1], 0, n) == (0, [])

    assert rsolve_hyper([1, 1], sqrt(n), n) is None
    assert rsolve_hyper([1, 1], n + sqrt(n), n) is None


@pytest.mark.slow
def test_bulk():
    funcs = [n, n + 1, n**2, n**3, n**4, n + n**2,
             27*n + 52*n**2 - 3*n**3 + 12*n**4 - 52*n**5]
    coeffs = [[-2, 1], [-2, -1, 1], [-1, 1, 1, -1, 1],
              [-n, 1], [n**2 - n + 12, 1]]
    for p in funcs:
        for c in coeffs:
            q = sum(c[i]*p.subs({n: n + i}) for i in range(len(c)))
            if p.is_polynomial(n):
                assert rsolve_poly(c, q, n)[0] == p
            if p.is_hypergeometric(n) and len(c) <= 3:
                assert rsolve_hyper(c, q, n)[0].subs({C0: 0, C1: 0,
                                                      C2: 0}).expand() == p


def test_rsolve():
    eq = f(n + 2) - f(n + 1) - f(n)
    res = [{f: Lambda(n, 2**(-n)*(C0*(1 + sqrt(5))**n +
                                  C1*(-sqrt(5) + 1)**n))}]

    assert rsolve(eq) == res

    res = [{k: v.subs({C0: sqrt(5), C1: -sqrt(5)}).simplify()
            for k, v in r.items()} for r in res]

    assert rsolve(eq, init={f(0): 0, f(1): 5}) == res
    assert rsolve(f(n) - f(n - 1) - f(n - 2), init={f(0): 0, f(1): 5}) == res
    assert rsolve(Eq(f(n), f(n - 1) + f(n - 2)), init={f(0): 0, f(1): 5}) == res

    eq = (n - 1)*f(n + 2) - (n**2 + 3*n - 2)*f(n + 1) + 2*n*(n + 1)*f(n)
    res = [{f: Lambda(n, C1*factorial(n) + C0*2**n)}]

    assert rsolve(eq) == res

    res = [{f: Lambda(n, -3*factorial(n) + 3*2**n)}]

    assert rsolve(eq, init={f(0): 0, f(1): 3}) == res

    eq = f(n) - f(n - 1) - 2

    assert rsolve(eq, f(n)) == [{f: Lambda(n, C0 + 2*n)}]
    assert rsolve(eq) == [{f: Lambda(n, C0 + 2*n)}]
    assert rsolve(eq, init={f(0): 0}) == [{f: Lambda(n, 2*n)}]
    assert rsolve(eq, init={f(0): 1}) == [{f: Lambda(n, 2*n + 1)}]
    assert rsolve(eq, init={f(0): 0, f(1): 1}) is None

    eq = 3*f(n - 1) - f(n) - 1

    assert rsolve(eq) == [{f: Lambda(n, 3**n*C0 + 1/2)}]
    assert rsolve(eq, init={f(0): 0}) == [{f: Lambda(n, -3**n/2 + 1/2)}]
    assert rsolve(eq, init={f(0): 1}) == [{f: Lambda(n, 3**n/2 + 1/2)}]
    assert rsolve(eq, init={f(0): 2}) == [{f: Lambda(n, 3*3**n/2 + 1/2)}]

    assert rsolve(f(n) - 1/n*f(n - 1),
                  f(n)) == [{f: Lambda(n, C0/factorial(n))}]
    assert rsolve(f(n) - 1/n*f(n - 1) - 1, f(n)) is None

    eq = 2*f(n - 1) + (1 - n)*f(n)/n

    assert rsolve(eq) == [{f: Lambda(n, 2**n*C0*n)}]
    assert rsolve([eq]) == [{f: Lambda(n, 2**n*C0*n)}]
    assert rsolve(eq, init={f(1): 1}) == [{f: Lambda(n, 2**(n - 1)*n)}]
    assert rsolve(eq, init={f(1): 2},
                  simplify=False) == [{f: Lambda(n, 2**(n - 1)*n*2)}]
    assert rsolve(eq, init={f(1): 2}) == [{f: Lambda(n, 2**n*n)}]
    assert rsolve(eq, init={f(1): 3}) == [{f: Lambda(n, 3*2**n*n/2)}]

    eq = (n - 1)*(n - 2)*f(n + 2) - (n + 1)*(n + 2)*f(n)

    assert rsolve(eq) == [{f: Lambda(n,
                                     n*(n - 2)*(n - 1)*((-1)**n*C1 + C0))}]
    assert rsolve(eq, init={f(3): 6,
                            f(4): 24}) == [{f: Lambda((n), n*(n - 1)*(n - 2))}]
    assert (rsolve(eq, init={f(3): 6, f(4): -24}) ==
            [{f: Lambda(n, (-1)**(n + 1)*n*(n - 2)*(n - 1))}])

    assert rsolve(Eq(f(n + 1), a*f(n)),
                  init={f(1): a}) == [{f: Lambda(n, a**n)}]

    assert (rsolve(f(n) - a*f(n - 2),
                   init={f(1): sqrt(a)*(a + b), f(2): a*(a - b)}) ==
            [{f: Lambda(n, a**(n/2)*(-(-1)**n*b + a))}])

    eq = (-16*n**2 + 32*n - 12)*f(n - 1) + (4*n**2 - 12*n + 9)*f(n)

    assert (rsolve(eq, init={f(1): binomial(2*n + 1, 3)}) ==
            [{f: Lambda(n, 4**n*n*(2*n - 1)*gamma(n + 3/2)/(3*gamma(n - 1/2)))}])

    assert (rsolve(f(n) + a*(f(n + 1) + f(n - 1))/2) ==
            [{f: Lambda(n, a**-n*(C0*(-a*sqrt(-1 + a**-2) - 1)**n +
                                  C1*(a*sqrt(-1 + a**-2) - 1)**n))}])


def test_rsolve_raises():
    pytest.raises(ValueError, lambda: rsolve(f(n) - f(k + 1), f(n)))
    pytest.raises(ValueError, lambda: rsolve(f(n) - sqrt(n)*f(n + 1)))
    pytest.raises(ValueError,
                  lambda: rsolve(f(n) - f(n + 1), f(n), init={g(0): 0}))
    pytest.raises(NotImplementedError, lambda: rsolve(f(n) - sqrt(n)))


def test_sympyissue_6844():
    eq = f(n + 2) - f(n + 1) + f(n)/4

    assert rsolve(eq) == [{f: Lambda(n, 2**(-n)*(C0 + C1*n))}]
    assert rsolve(eq, init={f(0): 0,
                            f(1): 1}) == [{f: Lambda(n, 2**(1 - n)*n)}]


def test_issue_294():
    eq = f(n) - f(n - 1) - 2*f(n - 2) - 2*n

    assert rsolve(eq) == [{f: Lambda(n, (-1)**n*C0 + 2**n*C1 - n - 5/2)}]
    # issue sympy/sympy#11261
    assert rsolve(eq, init={f(0): -1,
                            f(1): 1}) == [{f: Lambda(n, -(-1)**n/2 +
                                                     2**(n + 1) - n - 5/2)}]
    # issue sympy/sympy#7055
    assert rsolve(-2*f(n) + f(n + 1) +
                  n - 1) == [{f: Lambda(n, 2**n*C0 + n)}]


def test_sympyissue_8697():
    assert rsolve(f(n + 3) - f(n + 2) - f(n + 1) +
                  f(n)) == [{f: Lambda(n, (-1)**n*C1 + C0 + C2*n)}]
    assert (rsolve(f(n + 3) + 3*f(n + 2) + 3*f(n + 1) + f(n)) ==
            [{f: Lambda(n, (-1)**n*(C0 + C1*n + C2*n**2))}])

    assert rsolve(f(n) - 2*f(n - 3) + 5*f(n - 2) - 4*f(n - 1),
                  init={f(0): 1, f(1): 3,
                        f(2): 8}) == [{f: Lambda(n, 3*2**n - n - 2)}]

    # From issue thread (but not related to the problem, fixed before):
    assert rsolve(f(n) - 2*f(n - 1) - n,
                  init={f(0): 1}) == [{f: Lambda(n, 3*2**n - n - 2)}]
    assert (rsolve(f(n + 2) - 5*f(n + 1) + 6*f(n) - n) ==
            [{f: Lambda(n, 2**n*C0 + 3**n*C1 + n/2 + 3/4)}])


def test_issue_451():
    assert rsolve(f(n) - 2*f(n - 1) - 3**n,
                  init={f(0): 1}) == [{f: Lambda(n, -2**(n + 1) +
                                                 3**(n + 1))}]


def test_issue_456():
    assert rsolve(f(n) - 2*f(n - 1) - 3**n*n,
                  init={f(0): 1}) == [{f: Lambda(n, 7*2**n +
                                                 3**(n + 1)*(n - 2))}]


def test_sympyissue_13629():
    assert rsolve(f(n + 1) - (f(n) + (n + 1)**2),
                  init={f(0): 0}) == [{f: Lambda(n, n*(2*n**2 + 3*n + 1)/6)}]


def test_sympyissue_15553():
    assert rsolve(Eq(f(n + 1), 2*f(n) +
                     n**2 + 1)) == [{f: Lambda(n, 2**n*C0 - n**2 - 2*n - 4)}]
    assert rsolve(Eq(f(n + 1), 2*f(n) + n**2 + 1),
                  init={f(1): 0}) == [{f: Lambda(n, 7*2**n/2 - n**2 -
                                                 2*n - 4)}]


def test_issue_922():
    assert rsolve(-2*n/3 + f(n) - f(n - 1) + 2*(n - 1)**3/3 + 2*(n - 1)**2/3,
                  init={f(0): 0}) == [{f: Lambda(n, n*(-3*n**3 + 2*n**2 +
                                                       9*n + 4)/18)}]


def test_issue_923():
    assert rsolve(4*f(n) + 4*f(n + 1) +
                  f(n + 2)) == [{f: Lambda(n, (-2)**n*(C0 + C1*n))}]


def test_sympyissue_17982():
    assert (rsolve(f(n + 3) + 10*f(n + 2) + 32*f(n + 1) + 32*f(n)) ==
            [{f: Lambda(n, (-2)**n*C0 + (-4)**n*C1 + (-4)**n*C2*n)}])


def test_sympyissue_18751():
    r = symbols('r', real=True, positive=True)
    theta = symbols('theta', real=True)

    eq = f(n) - 2*r*cos(theta)*f(n - 1) + r**2*f(n - 2)
    res = [{f: Lambda(n, r**n*(C0*(cos(theta) - I*abs(sin(theta)))**n +
                      C1*(cos(theta) + I*abs(sin(theta)))**n))}]

    assert rsolve(eq) == res


def test_sympyissue_19630():
    eq = f(n + 3) - 3*f(n + 1) + 2*f(n)
    res = [{f: Lambda(n, (-2)**n*C1 + C0 + C2*n)}]

    assert rsolve(eq) == res

    res0 = [{f: Lambda(n, (-2)**n + 2*n)}]

    assert rsolve(eq, init={f(1): 0, f(2): 8, f(3): -2}) == res0
