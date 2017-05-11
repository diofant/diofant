"""Computational algebraic field theory. """

from functools import reduce

from mpmath import pslq, mp

from ..core import (S, Rational, AlgebraicNumber, Add, Mul, sympify,
                    Dummy, expand_mul, I, pi, GoldenRatio, Integer)
from ..functions import cos, sin, sqrt
from .polytools import (Poly, PurePoly, sqf_norm, invert, factor_list,
                        groebner, resultant, degree, poly_from_expr,
                        parallel_poly_from_expr, lcm)
from .polyerrors import (IsomorphismFailed, CoercionFailed, NotAlgebraic,
                         GeneratorsError)
from .rootoftools import RootOf
from .specialpolys import cyclotomic_poly
from .polyutils import dict_from_expr, expr_from_dict
from ..domains import ZZ, QQ
from .orthopolys import dup_chebyshevt
from .rings import ring
from .ring_series import rs_compose_add
from ..printing.lambdarepr import LambdaPrinter
from ..utilities import numbered_symbols, variations, lambdify, public, sift
from ..core.exprtools import Factors
from ..core.function import _mexpand
from ..simplify.radsimp import _split_gcd
from ..simplify.simplify import _is_sum_surds
from ..ntheory import sieve, divisors


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
                if abs(f.as_expr().evalf(prec1, points)) < eps:
                    candidates.append(f)
            if candidates:
                factors = candidates
            if len(factors) == 1:
                return factors[0]
            if prec1 > prec:
                break
            prec1 *= 2

    raise NotImplementedError("multiple candidates for the minimal polynomial of %s" % v)


def _separate_sq(p):
    """
    helper function for ``_minimal_polynomial_sq``

    It selects a rational ``g`` such that the polynomial ``p``
    consists of a sum of terms whose surds squared have gcd equal to ``g``
    and a sum of terms with surds squared prime with ``g``;
    then it takes the field norm to eliminate ``sqrt(g)``

    See simplify.simplify.split_surds and polytools.sqf_norm.

    Examples
    ========

    >>> from diofant import sqrt
    >>> from diofant.abc import x

    >>> p= -x + sqrt(2) + sqrt(3) + sqrt(7)
    >>> p = _separate_sq(p); p
    -x**2 + 2*sqrt(3)*x + 2*sqrt(7)*x - 2*sqrt(21) - 8
    >>> p = _separate_sq(p); p
    -x**4 + 4*sqrt(7)*x**3 - 32*x**2 + 8*sqrt(7)*x + 20
    >>> p = _separate_sq(p); p
    -x**8 + 48*x**6 - 536*x**4 + 1728*x**2 - 400

    """
    from ..utilities import sift

    def is_sqrt(expr):
        return expr.is_Pow and expr.exp is S.Half

    # p = c1*sqrt(q1) + ... + cn*sqrt(qn) -> a = [(c1, q1), .., (cn, qn)]
    a = []
    for y in p.args:
        if not y.is_Mul:
            if is_sqrt(y):
                a.append((S.One, y**2))
            elif y.is_Atom:
                a.append((y, S.One))
            elif y.is_Pow and y.exp.is_integer:
                a.append((y, S.One))
            else:  # pragma: no cover
                raise NotImplementedError
        else:
            sifted = sift(y.args, is_sqrt)
            a.append((Mul(*sifted[False]), Mul(*sifted[True])**2))
    a.sort(key=lambda z: z[1])
    if a[-1][1] is S.One:
        # there are no surds
        return p
    surds = [z for y, z in a]
    for i in range(len(surds)):
        if surds[i] != 1:
            break
    g, b1, b2 = _split_gcd(*surds[i:])
    a1 = []
    a2 = []
    for y, z in a:
        if z in b1:
            a1.append(y*z**S.Half)
        else:
            a2.append(y*z**S.Half)
    p1 = Add(*a1)
    p2 = Add(*a2)
    p = _mexpand(p1**2) - _mexpand(p2**2)
    return p


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

    >>> from diofant import sqrt
    >>> from diofant.abc import x

    >>> q = 1 + sqrt(2) + sqrt(3)
    >>> _minimal_polynomial_sq(q, 3, x)
    x**12 - 4*x**9 - 4*x**6 + 16*x**3 - 8
    """
    p = sympify(p)
    n = sympify(n)
    r = _is_sum_surds(p)
    if not n.is_Integer or not n > 0 or not _is_sum_surds(p):
        return
    pn = p**Rational(1, n)
    # eliminate the square roots
    p -= x
    while 1:
        p1 = _separate_sq(p)
        if p1 is p:
            p = p1.subs({x: x**n})
            break
        else:
            p = p1

    # _separate_sq eliminates field extensions in a minimal way, so that
    # if n = 1 then `p = constant*(minimal_polynomial(p))`
    # if n > 1 it contains the minimal polynomial as a factor.
    if n == 1:
        p1 = Poly(p)
        if p.coeff(x**p1.degree(x)) < 0:
            p = -p
        p = p.primitive()[1]
        return p
    # by construction `p` has root `pn`
    # the minimal polynomial is the factor vanishing in x = pn
    factors = factor_list(p)[1]

    result = _choose_factor(factors, x, pn)
    return result


def _minpoly_op_algebraic_element(op, ex1, ex2, x, dom, mp1=None, mp2=None):
    """
    return the minimal polynomial for ``op(ex1, ex2)``

    Parameters
    ==========

    op : operation ``Add`` or ``Mul``
    ex1, ex2 : expressions for the algebraic elements
    x : indeterminate of the polynomials
    dom: ground domain
    mp1, mp2 : minimal polynomials for ``ex1`` and ``ex2`` or None

    Examples
    ========

    >>> from diofant import sqrt, Add, Mul, QQ
    >>> from diofant.abc import x, y

    >>> p1 = sqrt(sqrt(2) + 1)
    >>> p2 = sqrt(sqrt(2) - 1)
    >>> _minpoly_op_algebraic_element(Mul, p1, p2, x, QQ)
    x - 1
    >>> q1 = sqrt(y)
    >>> q2 = 1 / y
    >>> _minpoly_op_algebraic_element(Add, q1, q2, x, QQ.frac_field(y))
    x**2*y**2 - 2*x*y - y**3 + 1

    References
    ==========

    .. [1] http://en.wikipedia.org/wiki/Resultant
    .. [2] I.M. Isaacs, Proc. Amer. Math. Soc. 25 (1970), 638
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
        if dom == QQ:
            R, X = ring('X', QQ)
            p1 = R(dict_from_expr(mp1)[0])
            p2 = R(dict_from_expr(mp2)[0])
        else:
            (p1, p2), _ = parallel_poly_from_expr((mp1, x - y), x, y)
            r = p1.compose(p2)
            mp1a = r.as_expr()

    elif op is Mul:
        mp1a = _muly(mp1, x, y)
    else:  # pragma: no cover
        raise NotImplementedError('option not available')

    if op is Mul or dom != QQ:
        r = resultant(mp1a, mp2, gens=[y, x])
    else:
        r = rs_compose_add(p1, p2)
        r = expr_from_dict(r.as_expr_dict(), x)

    deg1 = degree(mp1, x)
    deg2 = degree(mp2, y)
    if op is Mul and deg1 == 1 or deg2 == 1:
        # if deg1 = 1, then mp1 = x - a; mp1a = x - y - a;
        # r = mp2(x - a), so that `r` is irreducible
        return r

    r = Poly(r, x, domain=dom)
    _, factors = r.factor_list()
    res = _choose_factor(factors, x, op(ex1, ex2), dom)
    return res.as_expr()


def _invertx(p, x):
    """
    Returns ``expand_mul(x**degree(p, x)*p.subs(x, 1/x))``
    """
    p1 = poly_from_expr(p, x)[0]

    n = degree(p1)
    a = [c * x**(n - i) for (i,), c in p1.terms()]
    return Add(*a)


def _muly(p, x, y):
    """
    Returns ``_mexpand(y**deg*p.subs({x:x / y}))``
    """
    p1 = poly_from_expr(p, x)[0]

    n = degree(p1)
    a = [c * x**i * y**(n - i) for (i,), c in p1.terms()]
    return Add(*a)


def _minpoly_pow(ex, pw, x, dom, mp=None):
    """
    Returns ``minpoly(ex**pw, x)``

    Parameters
    ==========

    ex : algebraic element
    pw : rational number
    x : indeterminate of the polynomial
    dom: ground domain
    mp : minimal polynomial of ``p``

    Examples
    ========

    >>> from diofant import sqrt, QQ, Rational
    >>> from diofant.abc import x, y

    >>> p = sqrt(1 + sqrt(2))
    >>> _minpoly_pow(p, 2, x, QQ)
    x**2 - 2*x - 1
    >>> minpoly(p**2, x)
    x**2 - 2*x - 1
    >>> _minpoly_pow(y, Rational(1, 3), x, QQ.frac_field(y))
    x**3 - y
    >>> minpoly(y**Rational(1, 3), x)
    x**3 - y
    """
    pw = sympify(pw)
    if not mp:
        mp = _minpoly_compose(ex, x, dom)
    if not pw.is_rational:
        raise NotAlgebraic("%s doesn't seem to be an algebraic element" % ex)
    if pw < 0:
        if mp == x:
            raise ZeroDivisionError('%s is zero' % ex)
        mp = _invertx(mp, x)
        if pw == -1:
            return mp
        pw = -pw
        ex = 1/ex

    y = Dummy(str(x))
    mp = mp.subs({x: y})
    n, d = pw.as_numer_denom()
    res = Poly(resultant(mp, x**d - y**n, gens=[y]), x, domain=dom)
    _, factors = res.factor_list()
    res = _choose_factor(factors, x, ex**pw, dom)
    return res.as_expr()


def _minpoly_add(x, dom, *a):
    """
    returns ``minpoly(Add(*a), dom, x)``
    """
    mp = _minpoly_op_algebraic_element(Add, a[0], a[1], x, dom)
    p = a[0] + a[1]
    for px in a[2:]:
        mp = _minpoly_op_algebraic_element(Add, p, px, x, dom, mp1=mp)
        p = p + px
    return mp


def _minpoly_mul(x, dom, *a):
    """
    returns ``minpoly(Mul(*a), dom, x)``
    """
    mp = _minpoly_op_algebraic_element(Mul, a[0], a[1], x, dom)
    p = a[0] * a[1]
    for px in a[2:]:
        mp = _minpoly_op_algebraic_element(Mul, p, px, x, dom, mp1=mp)
        p = p * px
    return mp


def _minpoly_sin(ex, x):
    """
    Returns the minimal polynomial of ``sin(ex)``
    see http://mathworld.wolfram.com/TrigonometryAngles.html
    """
    c, a = ex.args[0].as_coeff_Mul()
    if a is pi:
        if c.is_rational:
            n = c.q
            q = sympify(n)
            if q.is_prime:
                # for a = pi*p/q with q odd prime, using chebyshevt
                # write sin(q*a) = mp(sin(a))*sin(a);
                # the roots of mp(x) are sin(pi*p/q) for p = 1,..., q - 1
                a = dup_chebyshevt(n, ZZ)
                return Add(*[x**(n - i - 1)*a[i] for i in range(n)])
            if c.p == 1:
                if q == 9:
                    return 64*x**6 - 96*x**4 + 36*x**2 - 3

            if n % 2 == 1:
                # for a = pi*p/q with q odd, use
                # sin(q*a) = 0 to see that the minimal polynomial must be
                # a factor of dup_chebyshevt(n, ZZ)
                a = dup_chebyshevt(n, ZZ)
                a = [x**(n - i)*a[i] for i in range(n + 1)]
                r = Add(*a)
                _, factors = factor_list(r)
                res = _choose_factor(factors, x, ex)
                return res

            expr = ((1 - cos(2*c*pi))/2)**S.Half
            res = _minpoly_compose(expr, x, QQ)
            return res

    raise NotAlgebraic("%s doesn't seem to be an algebraic element" % ex)


def _minpoly_cos(ex, x):
    """
    Returns the minimal polynomial of ``cos(ex)``
    see http://mathworld.wolfram.com/TrigonometryAngles.html
    """
    c, a = ex.args[0].as_coeff_Mul()
    if a is pi:
        if c.is_rational:
            if c.p == 1:
                if c.q == 7:
                    return 8*x**3 - 4*x**2 - 4*x + 1
                if c.q == 9:
                    return 8*x**3 - 6*x - 1
            elif c.p == 2:
                q = sympify(c.q)
                if q.is_prime:
                    s = _minpoly_sin(ex, x)
                    return _mexpand(s.subs({x: sqrt((1 - x)/2)}))

            # for a = pi*p/q, cos(q*a) =T_q(cos(a)) = (-1)**p
            n = int(c.q)
            a = dup_chebyshevt(n, ZZ)
            a = [x**(n - i)*a[i] for i in range(n + 1)]
            r = Add(*a) - (-1)**c.p
            _, factors = factor_list(r)
            res = _choose_factor(factors, x, ex)
            return res

    raise NotAlgebraic("%s doesn't seem to be an algebraic element" % ex)


def _minpoly_exp(ex, x):
    """
    Returns the minimal polynomial of ``exp(ex)``
    """
    c, a = ex.exp.as_coeff_Mul()
    p = sympify(c.p)
    q = sympify(c.q)
    if a == I*pi:
        if c.is_rational:
            if c.p == 1 or c.p == -1:
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
            mp = _choose_factor(factors, x, ex)
            return mp
        else:
            raise NotAlgebraic("%s doesn't seem to be an algebraic element" % ex)
    raise NotAlgebraic("%s doesn't seem to be an algebraic element" % ex)


def _minpoly_rootof(ex, x):
    """
    Returns the minimal polynomial of a ``RootOf`` object.
    """
    p = ex.expr
    p = p.subs({ex.poly.gens[0]: x})
    _, factors = factor_list(p, x)
    result = _choose_factor(factors, x, ex)
    return result


def _minpoly_compose(ex, x, dom):
    """
    Computes the minimal polynomial of an algebraic element
    using operations on minimal polynomials

    Examples
    ========

    >>> from diofant import minimal_polynomial, sqrt, Rational
    >>> from diofant.abc import x, y
    >>> minimal_polynomial(sqrt(2) + 3*Rational(1, 3), x, compose=True)
    x**2 - 2*x - 1
    >>> minimal_polynomial(sqrt(y) + 1/y, x, compose=True)
    x**2*y**2 - 2*x*y - y**3 + 1

    """
    if ex.is_Rational:
        return ex.q*x - ex.p
    if ex is I:
        return x**2 + 1
    if ex is GoldenRatio:
        return x**2 - x - 1
    if hasattr(dom, 'symbols') and ex in dom.symbols:
        return x - ex

    if dom.is_QQ and _is_sum_surds(ex):
        # eliminate the square roots
        ex -= x
        while 1:
            ex1 = _separate_sq(ex)
            if ex1 is ex:
                return ex
            else:
                ex = ex1

    if ex.is_Add:
        res = _minpoly_add(x, dom, *ex.args)
    elif ex.is_Mul:
        f = Factors(ex).factors
        r = sift(f.items(), lambda itx: itx[0].is_Rational and itx[1].is_Rational)
        if r[True] and dom == QQ:
            ex1 = Mul(*[bx**ex for bx, ex in r[False] + r[None]])
            r1 = r[True]
            dens = [y.q for _, y in r1]
            lcmdens = reduce(lcm, dens, 1)
            nums = [base**(y.p*lcmdens // y.q) for base, y in r1]
            ex2 = Mul(*nums)
            mp1 = minimal_polynomial(ex1, x)
            # use the fact that in Diofant canonicalization products of integers
            # raised to rational powers are organized in relatively prime
            # bases, and that in ``base**(n/d)`` a perfect power is
            # simplified with the root
            mp2 = ex2.q*x**lcmdens - ex2.p
            ex2 = ex2**Rational(1, lcmdens)
            res = _minpoly_op_algebraic_element(Mul, ex1, ex2, x, dom, mp1=mp1, mp2=mp2)
        else:
            res = _minpoly_mul(x, dom, *ex.args)
    elif ex.is_Pow:
        if ex.base is S.Exp1:
            res = _minpoly_exp(ex, x)
        else:
            res = _minpoly_pow(ex.base, ex.exp, x, dom)
    elif ex.__class__ is sin:
        res = _minpoly_sin(ex, x)
    elif ex.__class__ is cos:
        res = _minpoly_cos(ex, x)
    elif ex.__class__ is RootOf:
        res = _minpoly_rootof(ex, x)
    else:
        raise NotAlgebraic("%s doesn't seem to be an algebraic element" % ex)
    return res


@public
def minimal_polynomial(ex, x=None, **args):
    """
    Computes the minimal polynomial of an algebraic element.

    Parameters
    ==========

    ex : algebraic element expression
    x : independent variable of the minimal polynomial
    compose : boolean, optional
        If ``True`` (default), the minimal polynomial of the subexpressions
        of ``ex`` are computed, then the arithmetic operations on them are
        performed using the resultant and factorization.  Else a bottom-up
        algorithm is used with ``groebner``.  The default algorithm
        stalls less frequently.
    polys : boolean, optional
        if ``True`` returns a ``Poly`` object (the default is ``False``).
    domain : Domain, optional
        If no ground domain is given, it will be generated automatically
        from the expression.

    Examples
    ========

    >>> from diofant import minimal_polynomial, sqrt, solve, QQ
    >>> from diofant.abc import x, y

    >>> minimal_polynomial(sqrt(2), x)
    x**2 - 2
    >>> minimal_polynomial(sqrt(2), x, domain=QQ.algebraic_field(sqrt(2)))
    x - sqrt(2)
    >>> minimal_polynomial(sqrt(2) + sqrt(3), x)
    x**4 - 10*x**2 + 1
    >>> minimal_polynomial(solve(x**3 + x + 3)[0][x], x)
    x**3 + x + 3
    >>> minimal_polynomial(sqrt(y), x)
    x**2 - y

    """
    from .polytools import degree
    from ..domains import FractionField
    from ..core import preorder_traversal

    compose = args.get('compose', True)
    polys = args.get('polys', False)
    dom = args.get('domain', None)

    ex = sympify(ex)
    if ex.is_number:
        # not sure if it's always needed but try it for numbers (issue sympy/sympy#8354)
        ex = _mexpand(ex, recursive=True)
    if any(e.is_AlgebraicNumber for e in preorder_traversal(ex)):
        compose = False

    if x is not None:
        x, cls = sympify(x), Poly
    else:
        x, cls = Dummy('x'), PurePoly

    if not dom:
        dom = FractionField(QQ, list(ex.free_symbols)) if ex.free_symbols else QQ
    if hasattr(dom, 'symbols') and x in dom.symbols:
        raise GeneratorsError("the variable %s is an element of the ground domain %s" % (x, dom))

    if compose:
        result = _minpoly_compose(ex, x, dom)
        result = result.primitive()[1]
        c = result.coeff(x**degree(result, x))
        if c.is_negative:
            result = expand_mul(-result)
        return cls(result, x, field=True) if polys else result.collect(x)

    if not dom.is_QQ:
        raise NotImplementedError("groebner method only works for QQ")

    result = _minpoly_groebner(ex, x)
    return cls(result, x, field=True) if polys else result.collect(x)


def _minpoly_groebner(ex, x):
    """
    Computes the minimal polynomial of an algebraic number
    using Gröbner bases

    Examples
    ========

    >>> from diofant import minimal_polynomial, sqrt, Rational
    >>> from diofant.abc import x
    >>> minimal_polynomial(sqrt(2) + 3*Rational(1, 3), x, compose=False)
    x**2 - 2*x - 1

    References
    ==========

    .. [1] [Adams94]_
    """
    from .polytools import degree
    from ..core import expand_multinomial

    generator = numbered_symbols('a', cls=Dummy)
    mapping, symbols, replace = {}, {}, []

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
            if ex is S.ImaginaryUnit:
                return update_mapping(ex, 2, 1)
            elif ex.is_Rational:
                return ex
        elif ex.is_Add or ex.is_Mul:
            return ex.func(*[bottom_up_scan(g) for g in ex.args])
        elif ex.is_Pow:
            if ex.exp.is_Rational:
                base, exp = ex.base, ex.exp
                if exp.is_nonnegative:
                    if exp.is_noninteger:
                        base, exp = base**exp.p, Rational(1, exp.q)
                    base = bottom_up_scan(base)
                else:
                    bmp = PurePoly(_minpoly_groebner(1/base, x), x)
                    base, exp = update_mapping(1/base, bmp), -exp
                return update_mapping(ex, exp.q, -base**exp.p)
        elif ex.is_AlgebraicNumber:
            base = update_mapping(ex.root, ex.minpoly)
            res = Integer(0)
            for exp, coeff in ex.rep.terms():
                exp = Integer(exp[0])
                if exp:
                    res += coeff*update_mapping(base**exp, 1, -base**exp)
                else:
                    res += coeff
            return res
        elif isinstance(ex, RootOf) and ex.poly.domain.is_ZZ:
            return update_mapping(ex, ex.poly)

        raise NotAlgebraic("%s doesn't seem to be an algebraic number" % ex)

    if ex.is_Pow and ex.exp.is_negative:
        n, d = S.One, bottom_up_scan(1/ex)
    else:
        n, d = bottom_up_scan(ex), S.One

    F = [d*x - n] + list(mapping.values())
    G = groebner(F, list(symbols.values()) + [x], order='lex')

    _, factors = factor_list(G[-1])  # by construction G[-1] has root `ex`
    return _choose_factor(factors, x, ex)


minpoly = minimal_polynomial
__all__.append('minpoly')


def _coeffs_generator(n):
    """Generate coefficients for `primitive_element()`. """
    for coeffs in variations([1, -1], n, repetition=True):
        yield list(coeffs)


@public
def primitive_element(extension, x=None, **args):
    """Construct a common number field for all extensions. """
    if not extension:
        raise ValueError("can't compute primitive element for empty extension")

    if x is not None:
        x, cls = sympify(x), Poly
    else:
        x, cls = Dummy('x'), PurePoly
    if not args.get('ex', False):
        extension = [ AlgebraicNumber(ext, gen=x) for ext in extension ]

        g, coeffs = extension[0].minpoly.replace(x), [1]

        for ext in extension[1:]:
            s, _, g = sqf_norm(g, x, extension=ext)
            coeffs = [ s*c for c in coeffs ] + [1]

        if not args.get('polys', False):
            return g.as_expr(), coeffs
        else:
            return cls(g), coeffs

    generator = numbered_symbols('y', cls=Dummy)

    F, Y = [], []

    for ext in extension:
        y = next(generator)

        if ext.is_Poly:
            if ext.is_univariate:
                f = ext.as_expr(y)
            else:
                raise ValueError("expected minimal polynomial, got %s" % ext)
        else:
            f = minpoly(ext, y)

        F.append(f)
        Y.append(y)

    coeffs_generator = args.get('coeffs', _coeffs_generator)

    for coeffs in coeffs_generator(len(Y)):
        f = x - sum(c*y for c, y in zip(coeffs, Y))
        G = groebner(F + [f], Y + [x], order='lex', field=True)

        H, g = G[:-1], cls(G[-1], x, domain='QQ')

        for i, (h, y) in enumerate(zip(H, Y)):
            try:
                H[i] = Poly(y - h, x,
                            domain='QQ').all_coeffs()  # XXX: composite=False
            except CoercionFailed:  # pragma: no cover
                break  # G is not a triangular set
        else:
            break
    else:  # pragma: no cover
        raise RuntimeError("run out of coefficient configurations")

    _, g = g.clear_denoms()

    if not args.get('polys', False):
        return g.as_expr(), coeffs, H
    else:
        return g, coeffs, H


def is_isomorphism_possible(a, b):
    """Returns `True` if there is a chance for isomorphism. """
    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if m % n != 0:
        return False

    if n == m:
        return True

    da = a.minpoly.discriminant()
    db = b.minpoly.discriminant()

    i, k, half = 1, m//n, db//2

    while True:
        p = sieve[i]
        P = p**k

        if P > half:
            break

        if ((da % p) % 2) and not (db % P):
            return False

        i += 1

    return True


def field_isomorphism_pslq(a, b):
    """Construct field isomorphism using PSLQ algorithm. """
    if not a.root.is_extended_real or not b.root.is_extended_real:
        raise NotImplementedError("PSLQ doesn't support complex coefficients")

    f = a.minpoly
    g = b.minpoly.replace(f.gen)

    n, m, prev = 100, b.minpoly.degree(), None

    for i in range(1, 5):
        A = a.root.evalf(n)
        B = b.root.evalf(n)

        basis = [1, B] + [ B**i for i in range(2, m) ] + [A]

        dps, mp.dps = mp.dps, n
        coeffs = pslq(basis, maxcoeff=int(1e10), maxsteps=1000)
        mp.dps = dps

        if coeffs is None:
            break

        if coeffs != prev:
            prev = coeffs
        else:
            break

        coeffs = [sympify(c)/coeffs[-1] for c in coeffs[:-1]]

        while not coeffs[-1]:
            coeffs.pop()

        coeffs = list(reversed(coeffs))
        h = Poly(coeffs, f.gen, domain='QQ')

        if f.compose(h).rem(g).is_zero:
            d, approx = len(coeffs) - 1, 0

            for i, coeff in enumerate(coeffs):
                approx += coeff*B**(d - i)

            if A*approx < 0:
                return [ -c for c in coeffs ]
            else:
                return coeffs
        elif f.compose(-h).rem(g).is_zero:
            return [ -c for c in coeffs ]
        else:
            n *= 2

    return


def field_isomorphism_factor(a, b):
    """Construct field isomorphism via factorization. """
    _, factors = factor_list(a.minpoly, extension=b)

    for f, _ in factors:
        if f.degree() == 1:
            coeffs = f.rep.TC().to_diofant_list()
            d, terms = len(coeffs) - 1, []

            for i, coeff in enumerate(coeffs):
                terms.append(coeff*b.root**(d - i))

            root = Add(*terms)

            if (a.root - root).evalf(chop=True) == 0:
                return coeffs

            if (a.root + root).evalf(chop=True) == 0:
                return [ -c for c in coeffs ]
    else:
        return


@public
def field_isomorphism(a, b, **args):
    """Construct an isomorphism between two number fields. """
    a, b = sympify(a), sympify(b)

    if not a.is_AlgebraicNumber:
        a = AlgebraicNumber(a)

    if not b.is_AlgebraicNumber:
        b = AlgebraicNumber(b)

    if a == b:
        return a.coeffs()

    n = a.minpoly.degree()
    m = b.minpoly.degree()

    if n == 1:
        return [a.root]

    if m % n != 0:
        return

    if args.get('fast', True):
        try:
            result = field_isomorphism_pslq(a, b)

            if result is not None:
                return result
        except NotImplementedError:
            pass

    return field_isomorphism_factor(a, b)


@public
def to_number_field(extension, theta=None, **args):
    """Express `extension` in the field generated by `theta`. """
    gen = args.get('gen')

    if hasattr(extension, '__iter__'):
        extension = list(extension)
    else:
        extension = [extension]

    if len(extension) == 1 and type(extension[0]) is tuple:
        return AlgebraicNumber(extension[0])

    minpoly, coeffs = primitive_element(extension, gen, polys=True)
    root = sum(coeff*ext for coeff, ext in zip(coeffs, extension))

    if theta is None:
        return AlgebraicNumber((minpoly, root))
    else:
        theta = sympify(theta)

        if not theta.is_AlgebraicNumber:
            theta = AlgebraicNumber(theta, gen=gen)

        coeffs = field_isomorphism(root, theta)

        if coeffs is not None:
            return AlgebraicNumber(theta, coeffs)
        else:
            raise IsomorphismFailed(
                "%s is not in a subfield of %s" % (root, theta.root))


class IntervalPrinter(LambdaPrinter):
    """Use ``lambda`` printer but print numbers as ``mpi`` intervals. """

    def _print_Integer(self, expr):
        return "mpi('%s')" % super(IntervalPrinter, self)._print_Integer(expr)

    def _print_Rational(self, expr):
        return "mpi('%s')" % super(IntervalPrinter, self)._print_Rational(expr)

    def _print_Pow(self, expr):
        return super(IntervalPrinter, self)._print_Pow(expr, rational=True)


@public
def isolate(alg, eps=None, fast=False):
    """Give a rational isolating interval for an algebraic number. """
    alg = sympify(alg)

    if alg.is_Rational:
        return alg, alg
    elif not alg.is_extended_real:
        raise NotImplementedError(
            "complex algebraic numbers are not supported")

    func = lambdify((), alg, modules="mpmath", printer=IntervalPrinter())

    poly = minpoly(alg, polys=True)
    intervals = poly.intervals(sqf=True)

    dps, done = mp.dps, False

    try:
        while not done:
            alg = func()

            for a, b in intervals:
                if a <= alg.a and alg.b <= b:
                    done = True
                    break
            else:
                mp.dps *= 2
    finally:
        mp.dps = dps

    if eps is not None:
        a, b = poly.refine_root(a, b, eps=eps, fast=fast)

    return a, b
