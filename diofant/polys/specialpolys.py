"""Functions for generating interesting polynomials, e.g. for benchmarking."""

from ..core import Add, Dummy, Integer, Mul, symbols
from ..core.sympify import sympify
from ..domains import ZZ
from ..ntheory import nextprime
from ..utilities import subsets
from . import polytools, rings


__all__ = ('swinnerton_dyer_poly', 'cyclotomic_poly', 'symmetric_poly',
           'random_poly', 'interpolating_poly')


def swinnerton_dyer_poly(n, x=None, **args):
    """Generates n-th Swinnerton-Dyer polynomial in `x`."""
    from ..functions import sqrt
    from .numberfields import minimal_polynomial
    if n <= 0:
        raise ValueError(
            f"can't generate Swinnerton-Dyer polynomial of order {n}")

    if x is not None:
        x = sympify(x)
    else:
        x = Dummy('x')

    if n == 1:
        ex = x**2 - 2
    elif n == 2:
        ex = x**4 - 10*x**2 + 1
    elif n == 3:
        ex = x**8 - 40*x**6 + 352*x**4 - 960*x**2 + 576
    else:
        p = 2
        a = [sqrt(2)]
        for i in range(2, n + 1):
            p = nextprime(p)
            a.append(sqrt(p))
        ex = minimal_polynomial(Add(*a))(x)

    if not args.get('polys', False):
        return ex
    else:
        return polytools.PurePoly(ex, x)


def cyclotomic_poly(n, x=None, **args):
    """Generates cyclotomic polynomial of order `n` in `x`."""
    if n <= 0:
        raise ValueError(
            f"can't generate cyclotomic polynomial of order {n}")

    R = ZZ.inject('_0')
    poly = R._zz_cyclotomic_poly(int(n)).all_coeffs()

    if x is not None:
        poly = polytools.Poly(poly, x)
    else:
        poly = polytools.PurePoly(poly, Dummy('x'))

    if not args.get('polys', False):
        return poly.as_expr()
    else:
        return poly


def symmetric_poly(n, *gens, **args):
    """Generates symmetric polynomial of order `n`."""
    if n < 0 or n > len(gens) or not gens:
        raise ValueError(f"can't generate symmetric polynomial of order {n} for {gens}")
    elif not n:
        poly = Integer(1)
    else:
        poly = Add(*[Mul(*s) for s in subsets(gens, int(n))])

    if not args.get('polys', False):
        return poly
    else:
        return poly.as_poly(*gens)


def random_poly(x, n, inf, sup, domain=ZZ, polys=False, percent=None):
    """Return a polynomial of degree ``n`` with coefficients in ``[inf, sup]``."""
    ring = domain.inject(x)
    poly = polytools.Poly(dict(ring._random(n, inf, sup, percent)), x, domain=domain)

    if not polys:
        return poly.as_expr()
    else:
        return poly


def interpolating_poly(n, x, X='x', Y='y'):
    """Construct Lagrange interpolating polynomial for ``n`` data points."""
    if isinstance(X, str):
        X = symbols(f'{X}:{n}')

    if isinstance(Y, str):
        Y = symbols(f'{Y}:{n}')

    coeffs = []

    for i in range(n):
        numer = []
        denom = []

        for j in range(n):
            if i == j:
                continue

            numer.append(x - X[j])
            denom.append(X[i] - X[j])

        numer = Mul(*numer)
        denom = Mul(*denom)

        coeffs.append(numer/denom)

    return Add(*[coeff*y for coeff, y in zip(coeffs, Y)])


class _test_polys:
    def fateman_poly_F_1(self):
        """Fateman's GCD benchmark: trivial GCD."""
        gens = self.gens
        x, y = gens[0:2]

        return ((1 + sum(gens))*(2 + sum(gens)),
                (1 + sum(_**2 for _ in gens))*(-3*y*x**2 + y**2 - 1),
                self.one)

    def fateman_poly_F_2(self):
        """Fateman's GCD benchmark: linearly dense quartic inputs."""
        gens = self.gens
        x = gens[0]

        D = (1 + sum(gens))**2

        return (D*(-2 + x - sum(gens[1:]))**2,
                D*(2 + sum(gens))**2, D)

    def fateman_poly_F_3(self):
        """Fateman's GCD benchmark: sparse inputs (deg f ~ vars f)."""
        gens = self.gens
        x = gens[0]
        n = self.ngens

        D = (1 + sum(_**n for _ in gens))**2

        return (D*(-2 + x**n - sum(_**n for _ in gens[1:]))**2,
                D*(+2 + sum(_**n for _ in gens))**2, D)


# A few useful polynomials from Wang's paper ('78).


def _f_0():
    R, x, y, z = rings.ring('x y z', ZZ)
    return x**2*y*z**2 + 2*x**2*y*z + 3*x**2*y + 2*x**2 + 3*x + 4*y**2*z**2 + 5*y**2*z + 6*y**2 + y*z**2 + 2*y*z + y + 1


def _f_1():
    R, x, y, z = rings.ring('x y z', ZZ)
    return x**3*y*z + x**2*y**2*z**2 + x**2*y**2 + 20*x**2*y*z + 30*x**2*y + x**2*z**2 + 10*x**2*z + x*y**3*z + 30*x*y**2*z + 20*x*y**2 + x*y*z**3 + 10*x*y*z**2 + x*y*z + 610*x*y + 20*x*z**2 + 230*x*z + 300*x + y**2*z**2 + 10*y**2*z + 30*y*z**2 + 320*y*z + 200*y + 600*z + 6000


def _f_2():
    R, x, y, z = rings.ring('x y z', ZZ)
    return x**5*y**3 + x**5*y**2*z + x**5*y*z**2 + x**5*z**3 + x**3*y**2 + x**3*y*z + 90*x**3*y + 90*x**3*z + x**2*y**2*z - 11*x**2*y**2 + x**2*z**3 - 11*x**2*z**2 + y*z - 11*y + 90*z - 990


def _f_3():
    R, x, y, z = rings.ring('x y z', ZZ)
    return x**5*y**2 + x**4*z**4 + x**4 + x**3*y**3*z + x**3*z + x**2*y**4 + x**2*y**3*z**3 + x**2*y*z**5 + x**2*y*z + x*y**2*z**4 + x*y**2 + x*y*z**7 + x*y*z**3 + x*y*z**2 + y**2*z + y*z**4


def _f_4():
    R, x, y, z = rings.ring('x y z', ZZ)
    return -x**9*y**8*z - x**8*y**5*z**3 - x**7*y**12*z**2 - 5*x**7*y**8 - x**6*y**9*z**4 + x**6*y**7*z**3 + 3*x**6*y**7*z - 5*x**6*y**5*z**2 - x**6*y**4*z**3 + x**5*y**4*z**5 + 3*x**5*y**4*z**3 - x**5*y*z**5 + x**4*y**11*z**4 + 3*x**4*y**11*z**2 - x**4*y**8*z**4 + 5*x**4*y**7*z**2 + 15*x**4*y**7 - 5*x**4*y**4*z**2 + x**3*y**8*z**6 + 3*x**3*y**8*z**4 - x**3*y**5*z**6 + 5*x**3*y**4*z**4 + 15*x**3*y**4*z**2 + x**3*y**3*z**5 + 3*x**3*y**3*z**3 - 5*x**3*y*z**4 + x**2*z**7 + 3*x**2*z**5 + x*y**7*z**6 + 3*x*y**7*z**4 + 5*x*y**3*z**4 + 15*x*y**3*z**2 + y**4*z**8 + 3*y**4*z**6 + 5*z**6 + 15*z**4


def _f_5():
    R, x, y, z = rings.ring('x y z', ZZ)
    return -x**3 - 3*x**2*y + 3*x**2*z - 3*x*y**2 + 6*x*y*z - 3*x*z**2 - y**3 + 3*y**2*z - 3*y*z**2 + z**3


def _f_6():
    R, x, y, z, t = rings.ring('x y z t', ZZ)
    return 2115*x**4*y + 45*x**3*z**3*t**2 - 45*x**3*t**2 - 423*x*y**4 - 47*x*y**3 + 141*x*y*z**3 + 94*x*y*z*t - 9*y**3*z**3*t**2 + 9*y**3*t**2 - y**2*z**3*t**2 + y**2*t**2 + 3*z**6*t**2 + 2*z**4*t**3 - 3*z**3*t**2 - 2*z*t**3


def _w_1():
    R, x, y, z = rings.ring('x y z', ZZ)
    return 4*x**6*y**4*z**2 + 4*x**6*y**3*z**3 - 4*x**6*y**2*z**4 - 4*x**6*y*z**5 + x**5*y**4*z**3 + 12*x**5*y**3*z - x**5*y**2*z**5 + 12*x**5*y**2*z**2 - 12*x**5*y*z**3 - 12*x**5*z**4 + 8*x**4*y**4 + 6*x**4*y**3*z**2 + 8*x**4*y**3*z - 4*x**4*y**2*z**4 + 4*x**4*y**2*z**3 - 8*x**4*y**2*z**2 - 4*x**4*y*z**5 - 2*x**4*y*z**4 - 8*x**4*y*z**3 + 2*x**3*y**4*z + x**3*y**3*z**3 - x**3*y**2*z**5 - 2*x**3*y**2*z**3 + 9*x**3*y**2*z - 12*x**3*y*z**3 + 12*x**3*y*z**2 - 12*x**3*z**4 + 3*x**3*z**3 + 6*x**2*y**3 - 6*x**2*y**2*z**2 + 8*x**2*y**2*z - 2*x**2*y*z**4 - 8*x**2*y*z**3 + 2*x**2*y*z**2 + 2*x*y**3*z - 2*x*y**2*z**3 - 3*x*y*z + 3*x*z**3 - 2*y**2 + 2*y*z**2


def _w_2():
    R, x, y = rings.ring('x y', ZZ)
    return 24*x**8*y**3 + 48*x**8*y**2 + 24*x**7*y**5 - 72*x**7*y**2 + 25*x**6*y**4 + 2*x**6*y**3 + 4*x**6*y + 8*x**6 + x**5*y**6 + x**5*y**3 - 12*x**5 + x**4*y**5 - x**4*y**4 - 2*x**4*y**3 + 292*x**4*y**2 - x**3*y**6 + 3*x**3*y**3 - x**2*y**5 + 12*x**2*y**3 + 48*x**2 - 12*y**3


def f_polys():
    return _f_0(), _f_1(), _f_2(), _f_3(), _f_4(), _f_5(), _f_6()


def w_polys():
    return _w_1(), _w_2()
