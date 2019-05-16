"""High-level polynomials manipulation functions. """

import itertools

from ..core import Add, Integer, Mul
from ..utilities import numbered_symbols
from .polyerrors import (ComputationFailed, MultivariatePolynomialError,
                         PolificationFailed)
from .polyoptions import allowed_flags
from .polytools import Poly, parallel_poly_from_expr, poly_from_expr
from .specialpolys import interpolating_poly, symmetric_poly


__all__ = 'symmetrize', 'horner', 'interpolate', 'viete'


def symmetrize(F, *gens, **args):
    """
    Rewrite a polynomial in terms of elementary symmetric polynomials.

    A symmetric polynomial is a multivariate polynomial that remains invariant
    under any variable permutation, i.e., if ``f = f(x_1, x_2, ..., x_n)``,
    then ``f = f(x_{i_1}, x_{i_2}, ..., x_{i_n})``, where
    ``(i_1, i_2, ..., i_n)`` is a permutation of ``(1, 2, ..., n)`` (an
    element of the group ``S_n``).

    Returns a tuple of symmetric polynomials ``(f1, f2, ..., fn)`` such that
    ``f = f1 + f2 + ... + fn``.

    Examples
    ========

    >>> symmetrize(x**2 + y**2)
    (-2*x*y + (x + y)**2, 0)

    >>> symmetrize(x**2 + y**2, formal=True)
    (s1**2 - 2*s2, 0, [(s1, x + y), (s2, x*y)])

    >>> symmetrize(x**2 - y**2)
    (-2*x*y + (x + y)**2, -2*y**2)

    >>> symmetrize(x**2 - y**2, formal=True)
    (s1**2 - 2*s2, -2*y**2, [(s1, x + y), (s2, x*y)])

    """
    allowed_flags(args, ['formal', 'symbols'])

    iterable = True

    if not hasattr(F, '__iter__'):
        iterable = False
        F = [F]

    try:
        F, opt = parallel_poly_from_expr(F, *gens, **args)
    except PolificationFailed as exc:
        result = []

        for expr in exc.exprs:
            assert expr.is_Number
            result.append((expr, Integer(0)))

        if not iterable:
            result, = result

        if not exc.opt.formal:
            return result
        else:
            if iterable:
                return result, []
            else:
                return result + ([],)

    polys, symbols = [], opt.symbols
    gens, dom = opt.gens, opt.domain

    for i in range(len(gens)):
        poly = symmetric_poly(i + 1, gens, polys=True)
        polys.append((next(symbols), poly.set_domain(dom)))

    indices = range(len(gens) - 1)
    weights = range(len(gens), 0, -1)

    result = []

    for f in F:
        symmetric = []

        if not f.is_homogeneous:
            symmetric.append(f.TC())
            f -= f.TC()

        while f:
            _height, _monom, _coeff = -1, None, None

            for i, (monom, coeff) in enumerate(f.terms()):
                if all(monom[i] >= monom[i + 1] for i in indices):
                    height = max(n*m for n, m in zip(weights, monom))

                    if height > _height:
                        _height, _monom, _coeff = height, monom, coeff

            if _height != -1:
                monom, coeff = _monom, _coeff
            else:
                break

            exponents = []

            for m1, m2 in zip(monom, monom[1:] + (0,)):
                exponents.append(m1 - m2)

            term = [s**n for (s, _), n in zip(polys, exponents)]
            poly = [p**n for (_, p), n in zip(polys, exponents)]

            symmetric.append(Mul(coeff, *term))
            product = poly[0]*coeff

            for p in poly[1:]:
                product *= p

            f -= product

        result.append((Add(*symmetric), f.as_expr()))

    polys = [(s, p.as_expr()) for s, p in polys]

    if not opt.formal:
        for i, (sym, non_sym) in enumerate(result):
            result[i] = (sym.subs(polys), non_sym)

    if not iterable:
        result, = result

    if not opt.formal:
        return result
    else:
        if iterable:
            return result, polys
        else:
            return result + (polys,)


def horner(f, *gens, **args):
    """
    Rewrite a polynomial in Horner form.

    Among other applications, evaluation of a polynomial at a point is optimal
    when it is applied using the Horner scheme.

    Examples
    ========

    >>> from diofant.abc import e

    >>> horner(9*x**4 + 8*x**3 + 7*x**2 + 6*x + 5)
    x*(x*(x*(9*x + 8) + 7) + 6) + 5

    >>> horner(a*x**4 + b*x**3 + c*x**2 + d*x + e)
    e + x*(d + x*(c + x*(a*x + b)))

    >>> f = 4*x**2*y**2 + 2*x**2*y + 2*x*y**2 + x*y

    >>> horner(f, wrt=x)
    x*(x*y*(4*y + 2) + y*(2*y + 1))

    >>> horner(f, wrt=y)
    y*(x*y*(4*x + 2) + x*(2*x + 1))

    References
    ==========

    * https://en.wikipedia.org/wiki/Horner_scheme

    """
    allowed_flags(args, [])

    try:
        F, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        return exc.expr

    form, gen = Integer(0), F.gen

    if F.is_univariate:
        for coeff in F.all_coeffs():
            form = form*gen + coeff
    else:
        F, gens = Poly(F, gen), gens[1:]

        for coeff in F.all_coeffs():
            form = form*gen + horner(coeff, *gens, **args)

    return form


def interpolate(data, x):
    """
    Construct an interpolating polynomial for the data points.

    Examples
    ========

    A list is interpreted as though it were paired with a range starting
    from 1:

    >>> interpolate([1, 4, 9, 16], x)
    x**2

    This can be made explicit by giving a list of coordinates:

    >>> interpolate([(1, 1), (2, 4), (3, 9)], x)
    x**2

    The (x, y) coordinates can also be given as keys and values of a
    dictionary (and the points need not be equispaced):

    >>> interpolate([(-1, 2), (1, 2), (2, 5)], x)
    x**2 + 1
    >>> interpolate({-1: 2, 1: 2, 2: 5}, x)
    x**2 + 1

    """
    n = len(data)

    if isinstance(data, dict):
        X, Y = list(zip(*data.items()))
    else:
        if isinstance(data[0], tuple):
            X, Y = list(zip(*data))
        else:
            X = list(range(1, n + 1))
            Y = list(data)

    poly = interpolating_poly(n, x, X, Y)

    return poly.expand()


def viete(f, roots=None, *gens, **args):
    """
    Generate Viete's formulas for ``f``.

    Examples
    ========

    >>> r1, r2 = symbols('r1:3')

    >>> viete(a*x**2 + b*x + c, [r1, r2], x)
    [(r1 + r2, -b/a), (r1*r2, c/a)]

    """
    allowed_flags(args, [])

    try:
        f, opt = poly_from_expr(f, *gens, **args)
    except PolificationFailed as exc:
        raise ComputationFailed('viete', 1, exc)

    if f.is_multivariate:
        raise MultivariatePolynomialError(
            "multivariate polynomials are not allowed")

    n = f.degree()

    if n < 1:
        raise ValueError(
            "can't derive Viete's formulas for a constant polynomial")

    if roots is None:
        roots = numbered_symbols('r', start=1)

    roots = list(itertools.islice(roots, n))

    if n != len(roots):
        raise ValueError("required %s roots, got %s" % (n, len(roots)))

    lc, coeffs = f.LC(), f.all_coeffs()
    result, sign = [], -1

    for i, coeff in enumerate(coeffs[1:]):
        poly = symmetric_poly(i + 1, roots)
        coeff = sign*(coeff/lc)
        result.append((poly, coeff))
        sign = -sign

    return result
