"""Computational algebraic field theory."""

import functools
import math

import mpmath

from ..config import query
from ..core import (Add, Dummy, E, GoldenRatio, I, Integer, Mul, Rational,
                    cacheit, pi)
from ..core.exprtools import Factors
from ..core.function import _mexpand, count_ops
from ..core.sympify import sympify
from ..domains import QQ, AlgebraicField
from ..functions import Abs, conjugate, cos, exp_polar, im, re, root, sin, sqrt
from ..ntheory import divisors, factorint
from ..simplify.radsimp import _split_gcd
from ..simplify.simplify import _is_sum_surds
from ..utilities import lambdify, numbered_symbols, sift
from ..utilities.iterables import uniq
from .orthopolys import chebyshevt_poly
from .polyerrors import NotAlgebraic
from .polytools import (Poly, PurePoly, degree, factor_list, groebner, lcm,
                        parallel_poly_from_expr, resultant)
from .rootoftools import RootOf
from .specialpolys import cyclotomic_poly


__all__ = 'minimal_polynomial', 'primitive_element', 'field_isomorphism'


def _choose_factor(factors, x, v, dom=QQ, prec=200, bound=5):
    """
    Return a factor having root ``v``
    It is assumed that one of the factors has root ``v``.

    """
    if isinstance(factors[0], tuple):
        factors = [f[0] for f in factors]
    if len(factors) == 1:
        return factors[0]

    points = {x: v}
    symbols = dom.symbols if hasattr(dom, 'symbols') else []
    t = QQ(1, 10)

    for n in range(bound**len(symbols)):
        prec1 = 10
        n_temp = n
        for s in symbols:
            points[s] = n_temp % bound
            n_temp = n_temp // bound

        while True:
            candidates = []
            eps = t**(prec1 // 2)
            for f in factors:
                if abs(f.as_expr().evalf(prec1, points, strict=False)) < eps:
                    candidates.append(f)
            if candidates:
                factors = candidates
            if len(factors) == 1:
                return factors[0]
            if prec1 > prec:
                break
            prec1 *= 2
    else:
        raise NotImplementedError(f'multiple candidates for the minimal polynomial of {v}')


def _separate_sq(p):
    """
    Helper function for ``_minimal_polynomial_sq``.

    It selects a rational ``g`` such that the polynomial ``p``
    consists of a sum of terms whose surds squared have gcd equal to ``g``
    and a sum of terms with surds squared prime with ``g``;
    then it takes the field norm to eliminate ``sqrt(g)``

    See simplify.simplify.split_surds and polytools.sqf_norm.

    Examples
    ========

    >>> p = -x + sqrt(2) + sqrt(3) + sqrt(7)
    >>> p = _separate_sq(p)
    >>> p
    -x**2 + 2*sqrt(3)*x + 2*sqrt(7)*x - 2*sqrt(21) - 8
    >>> p = _separate_sq(p)
    >>> p
    -x**4 + 4*sqrt(7)*x**3 - 32*x**2 + 8*sqrt(7)*x + 20
    >>> p = _separate_sq(p)
    >>> p
    -x**8 + 48*x**6 - 536*x**4 + 1728*x**2 - 400

    """
    def is_sqrt(expr):
        return expr.is_Pow and expr.exp == Rational(1, 2)

    p = p.doit()

    # p = c1*sqrt(q1) + ... + cn*sqrt(qn) -> a = [(c1, q1), .., (cn, qn)]
    a = []
    for y in p.args:
        if not y.is_Mul:
            if is_sqrt(y):
                a.append((Integer(1), y**2))
            elif y.is_Atom:
                a.append((y, Integer(1)))
            else:
                raise NotImplementedError
        else:
            sifted = sift(y.args, is_sqrt)
            a.append((Mul(*sifted[False]), Mul(*sifted[True])**2))
    a.sort(key=lambda z: z[1])
    if a[-1][1] == 1:
        # there are no surds
        return p
    surds = [z for y, z in a]
    for i in range(len(surds)):  # pragma: no branch
        if surds[i] != 1:
            break
    g, b1, b2 = _split_gcd(*surds[i:])
    a1 = []
    a2 = []
    for y, z in a:
        if z in b1:
            a1.append(y*sqrt(z))
        else:
            a2.append(y*sqrt(z))
    p1 = Add(*a1)
    p2 = Add(*a2)
    return _mexpand(p1**2) - _mexpand(p2**2)


def _minimal_polynomial_sq(p, n, x):
    """
    Returns the minimal polynomial for the ``nth-root`` of a sum of surds
    or ``None`` if it fails.

    Parameters
    ==========

    p : sum of surds
    n : positive integer
    x : variable of the returned polynomial

    Examples
    ========

    >>> q = 1 + sqrt(2) + sqrt(3)
    >>> _minimal_polynomial_sq(q, 3, x)
    x**12 - 4*x**9 - 4*x**6 + 16*x**3 - 8

    """
    p = sympify(p)
    n = sympify(n)
    assert n.is_Integer and n > 1 and _is_sum_surds(p)
    pn = root(p, n)
    # eliminate the square roots
    p -= x
    while 1:
        p1 = _separate_sq(p)
        if p1 is p:
            p = p1.subs({x: x**n})
            break
        else:
            p = p1

    # by construction `p` has root `pn`
    # the minimal polynomial is the factor vanishing in x = pn
    factors = factor_list(p)[1]

    return _choose_factor(factors, x, pn)


def _minpoly_op_algebraic_element(op, ex1, ex2, x, dom, mp1=None, mp2=None):
    """
    Return the minimal polynomial for ``op(ex1, ex2)``.

    Parameters
    ==========

    op : operation ``Add`` or ``Mul``
    ex1, ex2 : expressions for the algebraic elements
    x : indeterminate of the polynomials
    dom: ground domain
    mp1, mp2 : minimal polynomials for ``ex1`` and ``ex2`` or None

    Examples
    ========

    >>> p1 = sqrt(sqrt(2) + 1)
    >>> p2 = sqrt(sqrt(2) - 1)
    >>> _minpoly_op_algebraic_element(Mul, p1, p2, x, QQ)
    x - 1
    >>> q1 = sqrt(y)
    >>> q2 = 1 / y
    >>> _minpoly_op_algebraic_element(Add, q1, q2, x, QQ.inject(y).field)
    x**2*y**2 - 2*x*y - y**3 + 1

    References
    ==========

    * https://en.wikipedia.org/wiki/Resultant
    * I.M. Isaacs, Proc. Amer. Math. Soc. 25 (1970), 638
      "Degrees of sums in a separable field extension".

    """
    y = Dummy(str(x))
    if mp1 is None:
        mp1 = _minpoly_compose(ex1, x, dom)
    if mp2 is None:
        mp2 = _minpoly_compose(ex2, y, dom)
    else:
        mp2 = mp2.subs({x: y})

    if op is Add:
        # mp1a = mp1.subs({x: x - y})
        (p1, p2), _ = parallel_poly_from_expr((mp1, x - y), x, y)
        r = p1.compose(p2)
        mp1a = r.as_expr()

    elif op is Mul:
        mp1a = _muly(mp1, x, y)
    else:
        raise NotImplementedError('option not available')

    r = resultant(mp1a, mp2, gens=[y, x])

    deg1 = degree(mp1, x)
    deg2 = degree(mp2, y)
    if op is Mul and deg1 == 1 or deg2 == 1:
        # if deg1 = 1, then mp1 = x - a; mp1a = x - y - a;
        # r = mp2(x - a), so that `r` is irreducible
        return r

    r = r.as_poly(x, domain=dom)
    _, factors = r.factor_list()
    res = _choose_factor(factors, x, op(ex1, ex2), dom)
    return res.as_expr()


def _invertx(p, x):
    """Returns ``expand_mul(x**degree(p, x)*p.subs({x: 1/x}))``."""
    (p1,) = parallel_poly_from_expr((p,), x)[0]

    n = degree(p1)
    a = [c * x**(n - i) for (i,), c in p1.terms()]
    return Add(*a)


def _muly(p, x, y):
    """Returns ``_mexpand(y**deg*p.subs({x:x / y}))``."""
    (p1,) = parallel_poly_from_expr((p,), x)[0]

    n = degree(p1)
    a = [c * x**i * y**(n - i) for (i,), c in p1.terms()]
    return Add(*a)


def _minpoly_pow(ex, pw, x, dom):
    """
    Returns ``minimal_polynomial(ex**pw)``

    Parameters
    ==========

    ex : algebraic element
    pw : rational number
    x : indeterminate of the polynomial
    dom: ground domain

    Examples
    ========

    >>> p = sqrt(1 + sqrt(2))
    >>> _minpoly_pow(p, 2, x, QQ)
    x**2 - 2*x - 1
    >>> minimal_polynomial(p**2)(x)
    x**2 - 2*x - 1
    >>> _minpoly_pow(y, Rational(1, 3), x, QQ.inject(y).field)
    x**3 - y
    >>> minimal_polynomial(cbrt(y))(x)
    x**3 - y

    """
    pw = sympify(pw)
    mp = _minpoly_compose(ex, x, dom)
    if not pw.is_rational:
        raise NotAlgebraic(f"{ex} doesn't seem to be an algebraic element")
    if pw < 0:
        if mp == x:
            raise ZeroDivisionError(f'{ex} is zero')
        mp = _invertx(mp, x)
        if pw == -1:
            return mp
        pw = -pw
        ex = 1/ex

    y = Dummy(str(x))
    mp = mp.subs({x: y})
    n, d = pw.as_numer_denom()
    res = resultant(mp, x**d - y**n, gens=[y]).as_poly(x, domain=dom)
    _, factors = res.factor_list()
    res = _choose_factor(factors, x, ex**pw, dom)
    return res.as_expr()


def _minpoly_add(x, dom, *a):
    """Returns ``minimal_polynomial(Add(*a), dom)``."""
    mp = _minpoly_op_algebraic_element(Add, a[0], a[1], x, dom)
    p = a[0] + a[1]
    for px in a[2:]:
        mp = _minpoly_op_algebraic_element(Add, p, px, x, dom, mp1=mp)
        p = p + px
    return mp


def _minpoly_mul(x, dom, *a):
    """Returns ``minimal_polynomial(Mul(*a), dom)``."""
    mp = _minpoly_op_algebraic_element(Mul, a[0], a[1], x, dom)
    p = a[0] * a[1]
    for px in a[2:]:
        mp = _minpoly_op_algebraic_element(Mul, p, px, x, dom, mp1=mp)
        p = p * px
    return mp


def _minpoly_sin(ex, x):
    """
    Returns the minimal polynomial of ``sin(ex)``
    see https://mathworld.wolfram.com/TrigonometryAngles.html

    """
    c, a = ex.args[0].as_coeff_Mul()
    if a is pi:
        n = c.denominator
        q = sympify(n)
        if q.is_prime:
            # for a = pi*p/q with q odd prime, using chebyshevt
            # write sin(q*a) = mp(sin(a))*sin(a);
            # the roots of mp(x) are sin(pi*p/q) for p = 1,..., q - 1
            a = chebyshevt_poly(n, polys=True).all_coeffs()
            return Add(*[x**(n - i - 1)*a[n - i] for i in range(n)])
        if c.numerator == 1:
            if q == 9:
                return 64*x**6 - 96*x**4 + 36*x**2 - 3

        if n % 2 == 1:
            # for a = pi*p/q with q odd, use
            # sin(q*a) = 0 to see that the minimal polynomial must be
            # a factor of chebyshevt_poly(n)
            a = chebyshevt_poly(n, polys=True).all_coeffs()
            a = [x**(n - i)*a[n - i] for i in range(n + 1)]
            r = Add(*a)
            _, factors = factor_list(r)
            res = _choose_factor(factors, x, ex)
            return res

        expr = sqrt((1 - cos(2*c*pi))/2)
        return _minpoly_compose(expr, x, QQ)

    raise NotAlgebraic(f"{ex} doesn't seem to be an algebraic element")


def _minpoly_cos(ex, x):
    """
    Returns the minimal polynomial of ``cos(ex)``
    see https://mathworld.wolfram.com/TrigonometryAngles.html

    """
    c, a = ex.args[0].as_coeff_Mul()
    if a is pi:
        if c.numerator == 1:
            if c.denominator == 7:
                return 8*x**3 - 4*x**2 - 4*x + 1
            elif c.denominator == 9:
                return 8*x**3 - 6*x - 1
        elif c.numerator == 2:
            q = sympify(c.denominator)
            if q.is_prime:
                s = _minpoly_sin(ex, x)
                return _mexpand(s.subs({x: sqrt((1 - x)/2)}))

        # for a = pi*p/q, cos(q*a) =T_q(cos(a)) = (-1)**p
        n = int(c.denominator)
        a = chebyshevt_poly(n, polys=True).all_coeffs()
        a = [x**(n - i)*a[n - i] for i in range(n + 1)]
        r = Add(*a) - (-1)**c.numerator
        _, factors = factor_list(r)
        return _choose_factor(factors, x, ex)

    raise NotAlgebraic(f"{ex} doesn't seem to be an algebraic element")


def _minpoly_exp(ex, x):
    """Returns the minimal polynomial of ``exp(ex)``."""
    c, a = ex.exp.as_coeff_Mul()
    q = sympify(c.denominator)
    if a == I*pi:
        if c.numerator == 1 or c.numerator == -1:
            if q == 3:
                return x**2 - x + 1
            if q == 4:
                return x**4 + 1
            if q == 6:
                return x**4 - x**2 + 1
            if q == 8:
                return x**8 + 1
            if q == 9:
                return x**6 - x**3 + 1
            if q == 10:
                return x**8 - x**6 + x**4 - x**2 + 1
            if q.is_prime:
                s = 0
                for i in range(q):
                    s += (-x)**i
                return s

        # x**(2*q) = product(factors)
        factors = [cyclotomic_poly(i, x) for i in divisors(2*q)]
        return _choose_factor(factors, x, ex)
    raise NotAlgebraic(f"{ex} doesn't seem to be an algebraic element")


def _minpoly_rootof(ex, x):
    """Returns the minimal polynomial of a ``RootOf`` object."""
    domain = ex.poly.domain
    if domain.is_IntegerRing:
        return ex.poly(x)
    else:
        return ex.poly.sqf_norm()[-1](x)


def _minpoly_compose(ex, x, dom):
    """
    Computes the minimal polynomial of an algebraic element
    using operations on minimal polynomials

    Examples
    ========

    >>> minimal_polynomial(sqrt(2) + 3*Rational(1, 3), method='compose')(x)
    x**2 - 2*x - 1
    >>> minimal_polynomial(sqrt(y) + 1/y, method='compose')(x)
    x**2*y**2 - 2*x*y - y**3 + 1

    """
    if ex.is_Rational:
        return ex.denominator*x - ex.numerator
    if ex is I:
        return x**2 + 1
    if ex is GoldenRatio:
        return x**2 - x - 1
    if ex == exp_polar(0):
        return x - 1
    if hasattr(dom, 'symbols') and ex in dom.symbols:
        return x - ex

    if dom.is_RationalField and _is_sum_surds(ex):
        # eliminate the square roots
        ex -= x
        while 1:
            ex1 = _separate_sq(ex)
            if ex1 is ex:
                return ex
            else:
                ex = ex1

    if ex.is_Add:
        res = _minpoly_add(x, dom, *sorted(ex.args, key=count_ops, reverse=True))
    elif ex.is_Mul:
        f = Factors(ex).factors
        r = sift(f.items(), lambda itx: itx[0].is_Rational and itx[1].is_Rational)
        if r[True] and dom == QQ:
            ex1 = Mul(*[bx**ex for bx, ex in r[False] + r[None]])
            r1 = r[True]
            dens = [y.denominator for _, y in r1]
            lcmdens = functools.reduce(lcm, dens, 1)
            nums = [base**(y.numerator*lcmdens // y.denominator) for base, y in r1]
            ex2 = Mul(*nums)
            mp1 = minimal_polynomial(ex1)(x)
            # use the fact that in Diofant canonicalization products of integers
            # raised to rational powers are organized in relatively prime
            # bases, and that in ``base**(n/d)`` a perfect power is
            # simplified with the root
            mp2 = ex2.denominator*x**lcmdens - ex2.numerator
            ex2 = Mul(*[bx**ex for bx, ex in r1])
            res = _minpoly_op_algebraic_element(Mul, ex1, ex2, x, dom, mp1=mp1, mp2=mp2)
        else:
            res = _minpoly_mul(x, dom, *sorted(ex.args, key=count_ops, reverse=True))
    elif ex.is_Pow:
        if ex.base is E:
            res = _minpoly_exp(ex, x)
        else:
            res = _minpoly_pow(ex.base, ex.exp, x, dom)
    elif isinstance(ex, sin):
        res = _minpoly_sin(ex, x)
    elif isinstance(ex, cos):
        res = _minpoly_cos(ex, x)
    elif isinstance(ex, RootOf) and ex.poly.domain.is_Numerical:
        res = _minpoly_rootof(ex, x)
    elif isinstance(ex, conjugate):
        res = _minpoly_compose(ex.args[0], x, dom)
    elif isinstance(ex, Abs):
        res = _minpoly_compose(sqrt(ex.args[0]*ex.args[0].conjugate()), x, dom)
    elif isinstance(ex, re):
        res = _minpoly_compose((ex.args[0] + ex.args[0].conjugate())/2, x, dom)
    elif isinstance(ex, im):
        res = _minpoly_compose((ex.args[0] - ex.args[0].conjugate())/2/I, x, dom)
    else:
        raise NotAlgebraic(f"{ex} doesn't seem to be an algebraic element")
    return res


@cacheit
def minimal_polynomial(ex, method=None, **args):
    """
    Computes the minimal polynomial of an algebraic element.

    Parameters
    ==========

    ex : algebraic element expression
    method : str, optional
        If ``compose``, the minimal polynomial of the subexpressions
        of ``ex`` are computed, then the arithmetic operations on them are
        performed using the resultant and factorization.  If ``groebner``,
        a bottom-up algorithm, using Gröbner bases is used.
        Defaults are determined by :func:`~diofant.config.setup`.
    domain : Domain, optional
        If no ground domain is given, it will be generated automatically
        from the expression.

    Examples
    ========

    >>> minimal_polynomial(sqrt(2))(x)
    x**2 - 2
    >>> minimal_polynomial(sqrt(2), domain=QQ.algebraic_field(sqrt(2)))(x)
    x - sqrt(2)
    >>> minimal_polynomial(sqrt(2) + sqrt(3))(x)
    x**4 - 10*x**2 + 1
    >>> minimal_polynomial(solve(x**3 + x + 3)[0][x])(x)
    x**3 + x + 3
    >>> minimal_polynomial(sqrt(y))(x)
    x**2 - y

    """
    if method is None:
        method = query('minpoly_method')
    _minpoly_methods = {'compose': _minpoly_compose, 'groebner': minpoly_groebner}
    try:
        _minpoly = _minpoly_methods[method]
    except KeyError:
        raise ValueError(f"'{method}' is not a valid algorithm for computing minimal "
                         ' polynomial')

    ex = sympify(ex)
    if ex.is_number:
        # not sure if it's always needed but try it for numbers (issue sympy/sympy#8354)
        ex = _mexpand(ex, recursive=True)

    x = Dummy('x')
    domain = args.get('domain',
                      QQ.inject(*ex.free_symbols).field if ex.free_symbols else QQ)

    result = _minpoly(ex, x, domain)
    _, factors = factor_list(result, x, domain=domain)
    result = _choose_factor(factors, x, ex, dom=domain)
    result = result.primitive()[1]

    return PurePoly(result, x, domain=domain)


def minpoly_groebner(ex, x, domain):
    """
    Computes the minimal polynomial of an algebraic number
    using Gröbner bases

    Examples
    ========

    >>> minimal_polynomial(sqrt(2) + 1, method='groebner')(x)
    x**2 - 2*x - 1

    References
    ==========

    * :cite:`Adams1994intro`

    """
    generator = numbered_symbols('a', cls=Dummy)
    mapping, symbols = {}, {}

    def update_mapping(ex, exp, base=None):
        if ex in mapping:
            return symbols[ex]

        a = next(generator)
        symbols[ex] = a

        if base is not None:
            mapping[ex] = a**exp + base
        else:
            mapping[ex] = exp.as_expr(a)

        return a

    def bottom_up_scan(ex):
        if ex.is_Atom:
            if ex is I:
                return update_mapping(ex, 2, 1)
            elif ex is GoldenRatio:
                return bottom_up_scan(ex.expand(func=True))
            elif ex.is_Rational:
                return ex
            elif ex.is_Symbol:
                return ex
        elif ex.is_Add or ex.is_Mul:
            return ex.func(*[bottom_up_scan(g) for g in ex.args])
        elif ex.is_Pow:
            if ex.exp.is_Rational:
                base, exp = ex.base, ex.exp
                if exp.is_nonnegative:
                    if exp.is_noninteger:
                        base, exp = base**exp.numerator, Rational(1, exp.denominator)
                    base = bottom_up_scan(base)
                else:
                    bmp = PurePoly(minpoly_groebner(1/base, x, domain=domain), x)
                    base, exp = update_mapping(1/base, bmp), -exp
                return update_mapping(ex, exp.denominator, -base**exp.numerator)
        elif isinstance(ex, RootOf) and ex.poly.domain.is_Numerical:
            if ex.poly.domain.is_IntegerRing:
                return update_mapping(ex, ex.poly)
            else:
                return update_mapping(ex, ex.poly.sqf_norm()[-1])
        elif isinstance(ex, conjugate):
            return update_mapping(ex, minimal_polynomial(ex.args[0], domain=domain,
                                                         method='groebner'))
        elif isinstance(ex, Abs):
            return bottom_up_scan(sqrt(ex.args[0]*ex.args[0].conjugate()))
        elif isinstance(ex, re):
            return bottom_up_scan((ex.args[0] + ex.args[0].conjugate())/2)
        elif isinstance(ex, im):
            return bottom_up_scan((ex.args[0] - ex.args[0].conjugate())/2/I)

        raise NotAlgebraic(f"{ex} doesn't seem to be an algebraic number")

    if ex.is_Pow and ex.exp.is_negative:
        n, d = Integer(1), bottom_up_scan(1/ex)
    else:
        n, d = bottom_up_scan(ex), Integer(1)

    F = [d*x - n] + list(mapping.values())
    G = groebner(F, *(list(symbols.values()) + [x]), order='lex', domain=domain)

    return G[-1]  # by construction G[-1] has root `ex`


def primitive_element(extension, **args):
    """Construct a common number field for all extensions.

    References
    ==========

    * :cite:`Yokoyama1989primitive`
    * :cite:`Arno1996alg`

    """
    if not extension:
        raise ValueError("can't compute primitive element for empty extension")

    extension = list(uniq(extension))

    x = Dummy('x')
    domain = args.get('domain', QQ)
    F = [minimal_polynomial(e, domain=domain) for e in extension]
    Y = [p.gen for p in F]

    for u in range(1, (len(F) - 1)*math.prod(f.degree() for f in F) + 1):
        coeffs = [u**n for n in range(len(Y))]
        f = x - sum(c*y for c, y in zip(coeffs, Y))

        *H, g = groebner(F + [f], *(Y + [x]), domain=domain)

        for i, (h, y) in enumerate(zip(H, Y)):
            H[i] = (y - h).eject(*Y).retract(field=True)
            if not (H[i].domain.is_RationalField or H[i].domain.is_AlgebraicField):
                break  # G is not a triangular set
            else:
                H[i] = H[i].set_domain(domain)
        else:
            g = g.eject(*Y).set_domain(domain)
            break
    else:
        if len(F) == 1:
            g, coeffs, H = F[0].replace(x), [Integer(1)], [x.as_poly(domain=domain)]
        else:  # pragma: no cover
            raise RuntimeError('run out of coefficient configurations')

    _, factors = factor_list(g, domain=domain)
    t = sum(c*e for c, e in zip(coeffs, extension))
    g = _choose_factor(factors, x, t, dom=domain)

    H = [h.rem(g).rep.all_coeffs() for h in H]

    _, g = PurePoly(g).clear_denoms(convert=True)

    if g.LC() != 1:
        H = [[c/g.LC()**n for n, c in enumerate(h)] for h in H]
        coeffs = [c*g.LC() for c in coeffs]
        g = (g.compose((g.gen/g.LC()).as_poly())*g.LC()**g.degree()//g.LC()).retract()

    return g, list(coeffs), H


def field_isomorphism_pslq(a, b):
    """Construct field isomorphism using PSLQ algorithm."""
    if not all(_.domain.is_RationalField and _.ext.is_real for _ in (a, b)):
        raise NotImplementedError("PSLQ doesn't support complex coefficients")

    f = a.minpoly
    x = f.gen

    g = b.minpoly.replace(x)
    m = g.degree()

    a, b = a.ext, b.ext

    for n in mpmath.libmp.libintmath.giant_steps(32, 256):  # pragma: no branch
        with mpmath.workdps(n):
            A, B = lambdify((), [a, b], 'mpmath')()
            basis = [A] + [B**i for i in reversed(range(m))]
            coeffs = mpmath.pslq(basis, maxcoeff=10**10, maxsteps=10**3)

        if coeffs:
            assert coeffs[0]  # basis[1:] elements are linearly independent

            h = -Poly(coeffs[1:], x, field=True).quo_ground(coeffs[0])

            if f.compose(h).rem(g).is_zero:
                return h.rep.all_coeffs()
        else:
            break


def field_isomorphism_factor(a, b):
    """Construct field isomorphism via factorization."""
    p = a.minpoly.set_domain(b)
    _, factors = p.factor_list()

    for f, _ in factors:
        if f.degree() == 1:
            root = -f.rep[(0,)]/f.rep[(1,)]

            if (a.ext - b.to_expr(root)).evalf(chop=True) == 0:
                return root.rep.all_coeffs()


def field_isomorphism(a, b, **args):
    """Construct an isomorphism between two number fields."""
    if not all(isinstance(_, AlgebraicField) for _ in (a, b)):
        raise ValueError(f'Arguments should be algebraic fields, got {a} and {b}')

    if a == b:
        return a.unit.rep.all_coeffs()

    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if a.domain == b.domain:
        if m % n:
            return
        elif a.domain.is_RationalField:
            da = a.minpoly.discriminant()
            db = b.minpoly.discriminant()
            k = m // n

            for p, q in factorint(da).items():
                if q % 2 and db % (p**k):
                    return

    if args.get('fast', True):
        try:
            result = field_isomorphism_pslq(a, b)

            if result is not None:
                return result
        except NotImplementedError:
            pass

    return field_isomorphism_factor(a, b)
