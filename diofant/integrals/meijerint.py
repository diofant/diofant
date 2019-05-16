"""
Integrate functions by rewriting them as Meijer G-functions.

There are three user-visible functions that can be used by other parts of the
diofant library to solve various integration problems:

- meijerint_indefinite
- meijerint_definite
- meijerint_inversion

They can be used to compute, respectively, indefinite integrals, definite
integrals over intervals of the real line, and inverse laplace-type integrals
(from c-I*oo to c+I*oo). See the respective docstrings for details.

The main references for this are:

[L] Luke, Y. L. (1969), The Special Functions and Their Approximations,
    Volume 1

[R] Kelly B. Roach.  Meijer G Function Representations.
    In: Proceedings of the 1997 International Symposium on Symbolic and
    Algebraic Computation, pages 205-211, New York, 1997. ACM.

[P] A. P. Prudnikov, Yu. A. Brychkov and O. I. Marichev (1990).
    Integrals and Series: More Special Functions, Vol. 3,.
    Gordon and Breach Science Publisher
"""

from collections import defaultdict

from ..core import (Add, Dummy, E, Eq, Expr, Function, I, Integer, Mul, Ne,
                    Pow, Rational, Tuple, Wild, cacheit, expand, expand_mul,
                    expand_power_base, factor_terms, nan, oo, pi, symbols,
                    sympify, zoo)
from ..core.compatibility import default_sort_key, ordered
from ..functions import Heaviside, Piecewise, cos, meijerg, piecewise_fold, sin
from ..functions.elementary.hyperbolic import (HyperbolicFunction,
                                               _rewrite_hyperbolics_as_exp)
from ..logic import And, Not, Or, false
from ..logic.boolalg import BooleanAtom
from ..simplify import collect, hyperexpand, powdenest
from ..simplify.fu import sincos_to_sum
from ..utilities.iterables import multiset_partitions
from ..utilities.misc import debug as _debug


# keep this at top for easy reference
z = Dummy('z')


def _has(res, *f):
    # return True if res has f; in the case of Piecewise
    # only return True if *all* pieces have f
    res = piecewise_fold(res)
    if getattr(res, 'is_Piecewise', False):
        return all(_has(i, *f) for i in res.args)
    return res.has(*f)


def _create_lookup_table(table):
    """Add formulae for the function -> meijerg lookup table."""
    def wild(n):
        return Wild(n, exclude=[z])
    p, q, a, b, c = list(map(wild, 'pqabc'))
    n = Wild('n', properties=[lambda x: x.is_Integer and x > 0])
    t = p*z**q

    def add(formula, an, ap, bm, bq, arg=t, fac=Integer(1), cond=True, hint=True):
        table[_mytype(formula, z)].append((formula,
                                           [(fac, meijerg(an, ap, bm, bq, arg))], cond, hint))

    def addi(formula, inst, cond, hint=True):
        table[_mytype(formula, z)].append((formula, inst, cond, hint))

    def constant(a):
        return [(a, meijerg([1], [], [], [0], z)),
                (a, meijerg([], [1], [0], [], z))]
    table[()] = [(a, constant(a), True, True)]

    # [P], Section 8.

    from ..logic import Not
    from ..functions import unpolarify

    class IsNonPositiveInteger(Function):

        @classmethod
        def eval(cls, arg):
            arg = unpolarify(arg)
            if arg.is_Integer is True:
                return arg <= 0

    # Section 8.4.2
    from ..functions import (gamma, cos, exp, re, sin, sqrt, sinh, cosh,
                             factorial, log, erf, erfc, erfi, polar_lift)
    # TODO this needs more polar_lift (c/f entry for exp)
    add(Heaviside(t - b)*(t - b)**(a - 1), [a], [], [], [0], t/b,
        gamma(a)*b**(a - 1), And(b > 0))
    add(Heaviside(b - t)*(b - t)**(a - 1), [], [a], [0], [], t/b,
        gamma(a)*b**(a - 1), And(b > 0))
    add(Heaviside(z - (b/p)**(1/q))*(t - b)**(a - 1), [a], [], [], [0], t/b,
        gamma(a)*b**(a - 1), And(b > 0))
    add(Heaviside((b/p)**(1/q) - z)*(b - t)**(a - 1), [], [a], [0], [], t/b,
        gamma(a)*b**(a - 1), And(b > 0))
    add((b + t)**(-a), [1 - a], [], [0], [], t/b, b**(-a)/gamma(a),
        hint=Not(IsNonPositiveInteger(a)))
    add(abs(b - t)**(-a), [1 - a], [(1 - a)/2], [0], [(1 - a)/2], t/b,
        2*sin(pi*a/2)*gamma(1 - a)*abs(b)**(-a), re(a) < 1)
    add((t**a - b**a)/(t - b), [0, a], [], [0, a], [], t/b,
        b**(a - 1)*sin(a*pi)/pi)

    add((z**a - b**a)/(z - b), [0, a], [], [0, a], [], z/b,
        b**(a - 1)*sin(a*pi)/pi)

    # 12
    def A1(r, sign, nu):
        return pi**(-Rational(1, 2))*(-sign*nu/2)**(1 - 2*r)

    def tmpadd(r, sgn):
        # XXX the a**2 is bad for matching
        add((sqrt(a**2 + t) + sgn*a)**b/(a**2 + t)**r,
            [(1 + b)/2, 1 - 2*r + b/2], [],
            [(b - sgn*b)/2], [(b + sgn*b)/2], t/a**2,
            a**(b - 2*r)*A1(r, sgn, b))
    tmpadd(0, 1)
    tmpadd(0, -1)
    tmpadd(Rational(1, 2), 1)
    tmpadd(Rational(1, 2), -1)

    # 13
    def tmpadd(r, sgn):
        add((sqrt(a + p*z**q) + sgn*sqrt(p)*z**(q/2))**b/(a + p*z**q)**r,
            [1 - r + sgn*b/2], [1 - r - sgn*b/2], [0, Rational(1, 2)], [],
            p*z**q/a, a**(b/2 - r)*A1(r, sgn, b))
    tmpadd(0, 1)
    tmpadd(0, -1)
    tmpadd(Rational(1, 2), 1)
    tmpadd(Rational(1, 2), -1)
    # (those after look obscure)

    # Section 8.4.3
    add(exp(t), [], [], [0], [], t/polar_lift(-1))

    # TODO can do sin^n, sinh^n by expansion ... where?
    # 8.4.4 (hyperbolic functions)
    add(sinh(t), [], [1], [Rational(1, 2)], [1, 0], t**2/4, pi**Rational(3, 2))
    add(cosh(t), [], [Rational(1, 2)], [0], [Rational(1, 2), Rational(1, 2)], t**2/4, pi**Rational(3, 2))

    # Section 8.4.5
    # TODO can do t + a. but can also do by expansion... (XXX not really)
    add(sin(t), [], [], [Rational(1, 2)], [0], t**2/4, sqrt(pi))
    add(cos(t), [], [], [0], [Rational(1, 2)], t**2/4, sqrt(pi))

    # Section 8.5.5
    def make_log1(subs):
        N = subs[n]
        return [((-1)**N*factorial(N),
                 meijerg([], [1]*(N + 1), [0]*(N + 1), [], t))]

    def make_log2(subs):
        N = subs[n]
        return [(factorial(N),
                 meijerg([1]*(N + 1), [], [], [0]*(N + 1), t))]
    # TODO these only hold for positive p, and can be made more general
    #      but who uses log(x)*Heaviside(a-x) anyway ...
    # TODO also it would be nice to derive them recursively ...
    addi(log(t)**n*Heaviside(1 - t), make_log1, True)
    addi(log(t)**n*Heaviside(t - 1), make_log2, True)

    def make_log3(subs):
        return make_log1(subs) + make_log2(subs)
    addi(log(t)**n, make_log3, True)
    addi(log(t + a),
         constant(log(a)) + [(Integer(1), meijerg([1, 1], [], [1], [0], t/a))],
         True)
    addi(log(abs(t - a)), constant(log(abs(a))) +
         [(pi, meijerg([1, 1], [Rational(1, 2)], [1], [0, Rational(1, 2)], t/a))],
         True)
    # TODO log(x)/(x+a) and log(x)/(x-1) can also be done. should they
    #      be derivable?
    # TODO further formulae in this section seem obscure

    # Sections 8.4.9-10
    # TODO

    # Section 8.4.11
    from ..functions import (Ei, expint, Si, Ci, Shi, Chi,
                             fresnels, fresnelc)
    addi(Ei(t),
         constant(-I*pi) + [(Integer(-1), meijerg([], [1], [0, 0], [],
                                                  t*polar_lift(-1)))],
         True)

    # Section 8.4.12
    add(Si(t), [1], [], [Rational(1, 2)], [0, 0], t**2/4, sqrt(pi)/2)
    add(Ci(t), [], [1], [0, 0], [Rational(1, 2)], t**2/4, -sqrt(pi)/2)

    # Section 8.4.13
    add(Shi(t), [Rational(1, 2)], [], [0], [Rational(-1, 2), Rational(-1, 2)], polar_lift(-1)*t**2/4,
        t*sqrt(pi)/4)
    add(Chi(t), [], [Rational(1, 2), 1], [0, 0], [Rational(1, 2), Rational(1, 2)], t**2/4, -
        pi**Rational(3, 2)/2)

    # generalized exponential integral
    add(expint(a, t), [], [a], [a - 1, 0], [], t)

    # Section 8.4.14
    add(erf(t), [1], [], [Rational(1, 2)], [0], t**2, 1/sqrt(pi))
    # TODO exp(-x)*erf(I*x) does not work
    add(erfc(t), [], [1], [0, Rational(1, 2)], [], t**2, 1/sqrt(pi))
    # This formula for erfi(z) yields a wrong(?) minus sign
    # add(erfi(t), [1], [], [Rational(1, 2)], [0], -t**2, I/sqrt(pi))
    add(erfi(t), [Rational(1, 2)], [], [0], [-Rational(1, 2)], -t**2, t/sqrt(pi))

    # Fresnel Integrals
    add(fresnels(t), [1], [], [Rational(3, 4)], [0, Rational(1, 4)], pi**2*t**4/16, Rational(1, 2))
    add(fresnelc(t), [1], [], [Rational(1, 4)], [0, Rational(3, 4)], pi**2*t**4/16, Rational(1, 2))

    # ##### bessel-type functions #####
    from ..functions import besselj, bessely, besseli, besselk

    # Section 8.4.19
    add(besselj(a, t), [], [], [a/2], [-a/2], t**2/4)

    # all of the following are derivable
    # add(sin(t)*besselj(a, t), [Rational(1, 4), Rational(3, 4)], [], [(1+a)/2],
    #     [-a/2, a/2, (1-a)/2], t**2, 1/sqrt(2))
    # add(cos(t)*besselj(a, t), [Rational(1, 4), Rational(3, 4)], [], [a/2],
    #     [-a/2, (1+a)/2, (1-a)/2], t**2, 1/sqrt(2))
    # add(besselj(a, t)**2, [Rational(1, 2)], [], [a], [-a, 0], t**2, 1/sqrt(pi))
    # add(besselj(a, t)*besselj(b, t), [0, Rational(1, 2)], [], [(a + b)/2],
    #    [-(a+b)/2, (a - b)/2, (b - a)/2], t**2, 1/sqrt(pi))

    # Section 8.4.20
    add(bessely(a, t), [], [-(a + 1)/2], [a/2, -a/2], [-(a + 1)/2], t**2/4)

    # TODO all of the following should be derivable
    # add(sin(t)*bessely(a, t), [Rational(1, 4), Rational(3, 4)], [(1 - a - 1)/2],
    #     [(1 + a)/2, (1 - a)/2], [(1 - a - 1)/2, (1 - 1 - a)/2, (1 - 1 + a)/2],
    #     t**2, 1/sqrt(2))
    # add(cos(t)*bessely(a, t), [Rational(1, 4), Rational(3, 4)], [(0 - a - 1)/2],
    #     [(0 + a)/2, (0 - a)/2], [(0 - a - 1)/2, (1 - 0 - a)/2, (1 - 0 + a)/2],
    #     t**2, 1/sqrt(2))
    # add(besselj(a, t)*bessely(b, t), [0, Rational(1, 2)], [(a - b - 1)/2],
    #     [(a + b)/2, (a - b)/2], [(a - b - 1)/2, -(a + b)/2, (b - a)/2],
    #     t**2, 1/sqrt(pi))
    # addi(bessely(a, t)**2,
    #      [(2/sqrt(pi), meijerg([], [Rational(1, 2), Rational(1, 2) - a], [0, a, -a],
    #                            [Rational(1, 2) - a], t**2)),
    #       (1/sqrt(pi), meijerg([Rational(1, 2)], [], [a], [-a, 0], t**2))],
    #      True)
    # addi(bessely(a, t)*bessely(b, t),
    #      [(2/sqrt(pi), meijerg([], [0, Rational(1, 2), (1 - a - b)/2],
    #                            [(a + b)/2, (a - b)/2, (b - a)/2, -(a + b)/2],
    #                            [(1 - a - b)/2], t**2)),
    #       (1/sqrt(pi), meijerg([0, Rational(1, 2)], [], [(a + b)/2],
    #                            [-(a + b)/2, (a - b)/2, (b - a)/2], t**2))],
    #      True)

    # Section 8.4.21 ?
    # Section 8.4.22
    add(besseli(a, t), [], [(1 + a)/2], [a/2], [-a/2, (1 + a)/2], t**2/4, pi)
    # TODO many more formulas. should all be derivable

    # Section 8.4.23
    add(besselk(a, t), [], [], [a/2, -a/2], [], t**2/4, Rational(1, 2))
    # TODO many more formulas. should all be derivable

    # Complete elliptic integrals K(z) and E(z)
    from ..functions import elliptic_k, elliptic_e
    add(elliptic_k(t), [Rational(1, 2), Rational(1, 2)], [], [0], [0], -t, Rational(1, 2))
    add(elliptic_e(t), [Rational(1, 2), 3*Rational(1, 2)], [], [0], [0], -t, -Rational(1, 2)/2)


####################################################################
# First some helper functions.
####################################################################


def _mytype(f, x):
    """Create a hashable entity describing the type of f."""
    if x not in f.free_symbols:
        return ()
    elif f.is_Function:
        return type(f),
    else:
        types = [_mytype(a, x) for a in f.args]
        res = []
        for t in types:
            res += list(t)
        res.sort(key=default_sort_key)
        return tuple(res)


class _CoeffExpValueError(ValueError):
    """Exception raised by _get_coeff_exp, for internal use only."""

    pass


def _get_coeff_exp(expr, x):
    """
    When expr is known to be of the form c*x**b, with c and/or b possibly 1,
    return c, b.

    >>> _get_coeff_exp(a*x**b, x)
    (a, b)
    >>> _get_coeff_exp(x, x)
    (1, 1)
    >>> _get_coeff_exp(2*x, x)
    (2, 1)
    >>> _get_coeff_exp(x**3, x)
    (1, 3)

    """
    from ..simplify import powsimp
    (c, m) = expand_power_base(powsimp(expr)).as_coeff_mul(x)
    if not m:
        return c, Integer(0)
    [m] = m
    if m.is_Pow:
        if m.base != x:
            raise _CoeffExpValueError('expr not of form a*x**b')
        return c, m.exp
    elif m == x:
        return c, Integer(1)
    else:
        raise _CoeffExpValueError('expr not of form a*x**b: %s' % expr)


def _exponents(expr, x):
    """
    Find the exponents of ``x`` (not including zero) in ``expr``.

    >>> _exponents(x, x)
    {1}
    >>> _exponents(x**2, x)
    {2}
    >>> _exponents(x**2 + x, x)
    {1, 2}
    >>> _exponents(x**3*sin(x + x**y) + 1/x, x)
    {-1, 1, 3, y}

    """
    def _exponents_(expr, x, res):
        if expr == x:
            res.update([1])
            return
        if expr.is_Pow and expr.base == x:
            res.update([expr.exp])
            return
        for arg in expr.args:
            _exponents_(arg, x, res)
    res = set()
    _exponents_(expr, x, res)
    return res


def _functions(expr, x):
    """Find the types of functions in expr, to estimate the complexity."""
    return ({e.func for e in expr.atoms(Function) if x in e.free_symbols} |
            {e.func for e in expr.atoms(Pow) if e.base is E and x in e.free_symbols})


def _find_splitting_points(expr, x):
    """
    Find numbers a such that a linear substitution x -> x + a would
    (hopefully) simplify expr.

    >>> _find_splitting_points(x, x)
    {0}
    >>> _find_splitting_points((x-1)**3, x)
    {1}
    >>> _find_splitting_points(sin(x+3)*x, x)
    {-3, 0}

    """
    p, q = [Wild(n, exclude=[x]) for n in 'pq']

    def compute_innermost(expr, res):
        if not isinstance(expr, Expr):
            return
        m = expr.match(p*x + q)
        if m and m[p] != 0:
            res.add(-m[q]/m[p])
            return
        if expr.is_Atom:
            return
        for arg in expr.args:
            compute_innermost(arg, res)
    innermost = set()
    compute_innermost(expr, innermost)
    return innermost


def _split_mul(f, x):
    """
    Split expression ``f`` into fac, po, g, where fac is a constant factor,
    po = x**s for some s independent of s, and g is "the rest".

    >>> from diofant.abc import s
    >>> _split_mul((3*x)**s*sin(x**2)*x, x)
    (3**s, x*x**s, sin(x**2))

    """
    from ..functions import polarify, unpolarify
    fac = Integer(1)
    po = Integer(1)
    g = Integer(1)
    f = expand_power_base(f)

    args = Mul.make_args(f)
    for a in args:
        if a == x:
            po *= x
        elif x not in a.free_symbols:
            fac *= a
        else:
            if a.is_Pow and x not in a.exp.free_symbols:
                c, t = a.base.as_coeff_mul(x)
                if t != (x,):
                    c, t = expand_mul(a.base).as_coeff_mul(x)
                if t == (x,):
                    po *= x**a.exp
                    fac *= unpolarify(polarify(c**a.exp, subs=False))
                    continue
            g *= a

    return fac, po, g


def _mul_args(f):
    """
    Return a list ``L`` such that Mul(*L) == f.

    If f is not a Mul or Pow, L=[f].
    If f=g**n for an integer n, L=[g]*n.
    If f is a Mul, L comes from applying _mul_args to all factors of f.

    """
    args = Mul.make_args(f)
    gs = []
    for g in args:
        if g.is_Pow and g.exp.is_Integer:
            n = g.exp
            base = g.base
            if n < 0:
                n = -n
                base = 1/base
            gs += [base]*n
        else:
            gs.append(g)
    return gs


def _mul_as_two_parts(f):
    """
    Find all the ways to split f into a product of two terms.
    Return None on failure.

    Although the order is canonical from multiset_partitions, this is
    not necessarily the best order to process the terms. For example,
    if the case of len(gs) == 2 is removed and multiset is allowed to
    sort the terms, some tests fail.

    >>> list(ordered(_mul_as_two_parts(x*sin(x)*exp(x))))
    [(x, E**x*sin(x)), (E**x*x, sin(x)), (x*sin(x), E**x)]

    """

    gs = _mul_args(f)
    if len(gs) < 2:
        return
    if len(gs) == 2:
        if ((gs[0].is_Pow and gs[0].base is E) and
                (not gs[1].is_Pow or gs[1].base is not E)):
            gs = [gs[1], gs[0]]
        return [tuple(gs)]
    return [(Mul(*x), Mul(*y)) for (x, y) in multiset_partitions(gs, 2)]


def _inflate_g(g, n):
    """ Return C, h such that h is a G function of argument z**n and
    g = C*h.

    """
    # TODO should this be a method of meijerg?
    # See: [L, page 150, equation (5)]
    def inflate(params, n):
        """(a1, .., ak) -> (a1/n, (a1+1)/n, ..., (ak + n-1)/n)."""
        res = []
        for a in params:
            for i in range(n):
                res.append((a + i)/n)
        return res
    v = Integer(len(g.ap) - len(g.bq))
    C = n**(1 + g.nu + v/2)
    C /= (2*pi)**((n - 1)*g.delta)
    return C, meijerg(inflate(g.an, n), inflate(g.aother, n),
                      inflate(g.bm, n), inflate(g.bother, n),
                      g.argument**n * n**(n*v))


def _flip_g(g):
    """ Turn the G function into one of inverse argument
    (i.e. G(1/x) -> G'(x))

    """
    # See [L], section 5.2
    def tr(l):
        return [1 - a for a in l]
    return meijerg(tr(g.bm), tr(g.bother), tr(g.an), tr(g.aother), 1/g.argument)


def _inflate_fox_h(g, a):
    r"""
    Let d denote the integrand in the definition of the G function ``g``.
    Consider the function H which is defined in the same way, but with
    integrand d/Gamma(a*s) (contour conventions as usual).

    If a is rational, the function H can be written as C*G, for a constant C
    and a G-function G.

    This function returns C, G.

    """
    if a < 0:
        return _inflate_fox_h(_flip_g(g), -a)
    p = Integer(a.numerator)
    q = Integer(a.denominator)
    # We use the substitution s->qs, i.e. inflate g by q. We are left with an
    # extra factor of Gamma(p*s), for which we use Gauss' multiplication
    # theorem.
    D, g = _inflate_g(g, q)
    z = g.argument
    D /= (2*pi)**((1 - p)/2)*p**(-Rational(1, 2))
    z /= p**p
    bs = [(n + 1)/p for n in range(p)]
    return D, meijerg(g.an, g.aother, g.bm, list(g.bother) + bs, z)


_dummies = {}


def _dummy(name, token, expr, **kwargs):
    """
    Return a dummy. This will return the same dummy if the same token+name is
    requested more than once, and it is not already in expr.
    This is for being cache-friendly.

    """
    d = _dummy_(name, token, **kwargs)
    if d in expr.free_symbols:
        return Dummy(name, **kwargs)
    return d


def _dummy_(name, token, **kwargs):
    """
    Return a dummy associated to name and token. Same effect as declaring
    it globally.

    """
    global _dummies
    if not (name, token) in _dummies:
        _dummies[(name, token)] = Dummy(name, **kwargs)
    return _dummies[(name, token)]


def _is_analytic(f, x):
    """ Check if f(x), when expressed using G functions on the positive reals,
    will in fact agree with the G functions almost everywhere

    """
    from ..functions import Heaviside, Abs
    return not any(x in expr.free_symbols for expr in f.atoms(Heaviside, Abs))


def _condsimp(cond):
    """
    Do naive simplifications on ``cond``.

    Note that this routine is completely ad-hoc, simplification rules being
    added as need arises rather than following any logical pattern.

    >>> _condsimp(Or(x < y, z, Eq(x, y)))
    z | (x <= y)
    >>> _condsimp(Or(x <= y, And(x < y, z)))
    x <= y

    """
    from ..functions import (unbranched_argument, exp_polar,
                             periodic_argument, polar_lift)
    from ..logic.boolalg import BooleanFunction
    if not isinstance(cond, BooleanFunction):
        return cond
    cond = cond.func(*list(map(_condsimp, cond.args)))
    change = True
    p, q, r = symbols('p q r', cls=Wild)
    rules = [
        (Or(p < q, Eq(p, q)), p <= q),
        # The next two obviously are instances of a general pattern, but it is
        # easier to spell out the few cases we care about.
        (And(abs(unbranched_argument(p)) <= pi,
             abs(unbranched_argument(exp_polar(-2*pi*I)*p)) <= pi),
         Eq(unbranched_argument(exp_polar(-I*pi)*p), 0)),
        (And(abs(unbranched_argument(p)) <= pi/2,
             abs(unbranched_argument(exp_polar(-pi*I)*p)) <= pi/2),
         Eq(unbranched_argument(exp_polar(-I*pi/2)*p), 0)),
        (Or(p <= q, And(p < q, r)), p <= q)
    ]
    while change:
        change = False
        for fro, to in rules:
            if fro.func != cond.func:
                continue
            for n, arg in enumerate(cond.args):
                if r in fro.args[0].free_symbols:
                    m = arg.match(fro.args[1])
                    num = 1
                else:
                    num = 0
                    m = arg.match(fro.args[0])
                if not m:
                    continue
                otherargs = [x.subs(m) for x in fro.args[:num] + fro.args[num + 1:]]
                otherlist = [n]
                for arg2 in otherargs:
                    for k, arg3 in enumerate(cond.args):
                        if k in otherlist:
                            continue
                        if arg2 == arg3:
                            otherlist += [k]
                            break
                        if isinstance(arg3, And) and arg2.args[1] == r and \
                                isinstance(arg2, And) and arg2.args[0] in arg3.args:
                            otherlist += [k]
                            break
                        if isinstance(arg3, And) and arg2.args[0] == r and \
                                isinstance(arg2, And) and arg2.args[1] in arg3.args:
                            otherlist += [k]
                            break
                if len(otherlist) != len(otherargs) + 1:
                    continue
                newargs = [arg for (k, arg) in enumerate(cond.args)
                           if k not in otherlist] + [to.subs(m)]
                cond = cond.func(*newargs)
                change = True
                break
    # final tweak

    def repl_eq(orig):
        if orig.lhs == 0:
            expr = orig.rhs
        elif orig.rhs == 0:
            expr = orig.lhs
        else:
            return orig
        m = expr.match(unbranched_argument(polar_lift(p)**q))
        if not m:
            if isinstance(expr, periodic_argument) and not expr.args[0].is_polar \
                    and expr.args[1] == oo:
                return expr.args[0] > 0
            return orig
        return m[p] > 0
    return cond.replace(
        lambda expr: expr.is_Relational and expr.rel_op == '==',
        repl_eq)


def _eval_cond(cond):
    """Re-evaluate the conditions."""
    if isinstance(cond, bool):
        return cond
    return _condsimp(cond.doit())

####################################################################
# Now the "backbone" functions to do actual integration.
####################################################################


def _my_principal_branch(expr, period, full_pb=False):
    """ Bring expr nearer to its principal branch by removing
    superfluous factors.

    This function does *not* guarantee to yield the principal branch,
    to avoid introducing opaque principal_branch() objects,
    unless full_pb=True.

    """
    from ..functions import principal_branch
    res = principal_branch(expr, period)
    if not full_pb:
        res = res.replace(principal_branch, lambda x, y: x)
    return res


def _rewrite_saxena_1(fac, po, g, x):
    """
    Rewrite the integral fac*po*g dx, from zero to infinity, as
    integral fac*G, where G has argument a*x. Note po=x**s.
    Return fac, G.

    """
    _, s = _get_coeff_exp(po, x)
    a, b = _get_coeff_exp(g.argument, x)
    period = g.get_period()
    a = _my_principal_branch(a, period)

    # We substitute t = x**b.
    C = fac/(abs(b)*a**((s + 1)/b - 1))
    # Absorb a factor of (at)**((1 + s)/b - 1).

    def tr(l):
        return [a + (1 + s)/b - 1 for a in l]
    return C, meijerg(tr(g.an), tr(g.aother), tr(g.bm), tr(g.bother),
                      a*x)


def _check_antecedents_1(g, x, helper=False):
    r"""
    Return a condition under which the mellin transform of g exists.
    Any power of x has already been absorbed into the G function,
    so this is just int_0^\infty g dx.

    See [L, section 5.6.1]. (Note that s=1.)

    If ``helper`` is True, only check if the MT exists at infinity, i.e. if
    int_1^\infty g dx exists.

    """
    # NOTE if you update these conditions, please update the documentation as well
    from ..functions import ceiling, re, unbranched_argument as arg
    delta = g.delta
    eta, _ = _get_coeff_exp(g.argument, x)
    m, n, p, q = sympify([len(g.bm), len(g.an), len(g.ap), len(g.bq)])

    if p > q:
        def tr(l):
            return [1 - x for x in l]
        return _check_antecedents_1(meijerg(tr(g.bm), tr(g.bother),
                                            tr(g.an), tr(g.aother), x/eta),
                                    x)

    tmp = []
    for b in g.bm:
        tmp += [-re(b) < 1]
    for a in g.an:
        tmp += [1 < 1 - re(a)]
    cond_3 = And(*tmp)

    for b in g.bother:
        tmp += [-re(b) < 1]
    for a in g.aother:
        tmp += [1 < 1 - re(a)]
    cond_3_star = And(*tmp)

    cond_4 = (-re(g.nu) + (q + 1 - p)/2 > q - p)

    def debug(*msg):
        _debug(*msg)

    debug('Checking antecedents for 1 function:')
    debug('  delta=%s, eta=%s, m=%s, n=%s, p=%s, q=%s'
          % (delta, eta, m, n, p, q))
    debug('  ap = %s, %s' % (list(g.an), list(g.aother)))
    debug('  bq = %s, %s' % (list(g.bm), list(g.bother)))
    debug('  cond_3=%s, cond_3*=%s, cond_4=%s' % (cond_3, cond_3_star, cond_4))

    conds = []

    # case 1
    case1 = []
    tmp1 = [1 <= n, p < q, 1 <= m]
    tmp2 = [1 <= p, 1 <= m, Eq(q, p + 1), Not(And(Eq(n, 0), Eq(m, p + 1)))]
    tmp3 = [1 <= p, Eq(q, p)]
    for k in range(ceiling(delta/2) + 1):
        tmp3 += [Ne(abs(arg(eta)), (delta - 2*k)*pi)]
    tmp = [delta > 0, abs(arg(eta)) < delta*pi]
    extra = [Ne(eta, 0), cond_3]
    if helper:
        extra = []
    for t in [tmp1, tmp2, tmp3]:
        case1 += [And(*(t + tmp + extra))]
    conds += case1
    debug('  case 1:', case1)

    # case 2
    extra = [cond_3]
    if helper:
        extra = []
    case2 = [And(Eq(n, 0), p + 1 <= m, m <= q,
                 abs(arg(eta)) < delta*pi, *extra)]
    conds += case2
    debug('  case 2:', case2)

    # case 3
    extra = [cond_3, cond_4]
    if helper:
        extra = []
    case3 = [And(p < q, 1 <= m, delta > 0, Eq(abs(arg(eta)), delta*pi),
                 *extra)]
    case3 += [And(p <= q - 2, Eq(delta, 0), Eq(abs(arg(eta)), 0), *extra)]
    conds += case3
    debug('  case 3:', case3)

    # TODO altered cases 4-7

    # extra case from wofram functions site:
    # (reproduced verbatim from Prudnikov, section 2.24.2)
    # http://functions.wolfram.com/HypergeometricFunctions/MeijerG/21/02/01/
    case_extra = []
    case_extra += [Eq(p, q), Eq(delta, 0), Eq(arg(eta), 0), Ne(eta, 0)]
    if not helper:
        case_extra += [cond_3]
    s = []
    for a, b in zip(g.ap, g.bq):
        s += [b - a]
    case_extra += [re(Add(*s)) < 0]
    case_extra = And(*case_extra)
    conds += [case_extra]
    debug('  extra case:', [case_extra])

    case_extra_2 = [And(delta > 0, abs(arg(eta)) < delta*pi)]
    if not helper:
        case_extra_2 += [cond_3]
    case_extra_2 = And(*case_extra_2)
    conds += [case_extra_2]
    debug('  second extra case:', [case_extra_2])

    # TODO This leaves only one case from the three listed by Prudnikov.
    #      Investigate if these indeed cover everything; if so, remove the rest.

    return Or(*conds)


def _int0oo_1(g, x):
    r"""
    Evaluate int_0^\infty g dx using G functions,
    assuming the necessary conditions are fulfilled.

    >>> _int0oo_1(meijerg([a], [b], [c], [d], x*y), x)
    gamma(-a)*gamma(c + 1)/(y*gamma(-d)*gamma(b + 1))

    """
    # See [L, section 5.6.1]. Note that s=1.
    from ..functions import gamma, unpolarify
    from ..simplify import combsimp
    eta, _ = _get_coeff_exp(g.argument, x)
    res = 1/eta
    # XXX TODO we should reduce order first
    for b in g.bm:
        res *= gamma(b + 1)
    for a in g.an:
        res *= gamma(1 - a - 1)
    for b in g.bother:
        res /= gamma(1 - b - 1)
    for a in g.aother:
        res /= gamma(a + 1)
    return combsimp(unpolarify(res))


def _rewrite_saxena(fac, po, g1, g2, x, full_pb=False):
    """
    Rewrite the integral fac*po*g1*g2 from 0 to oo in terms of G functions
    with argument c*x.
    Return C, f1, f2 such that integral C f1 f2 from 0 to infinity equals
    integral fac po g1 g2 from 0 to infinity.

    >>> from diofant.abc import s
    >>> g1 = meijerg([], [], [0], [], s*t)
    >>> g2 = meijerg([], [], [m/2], [-m/2], t**2/4)
    >>> r = _rewrite_saxena(1, t**0, g1, g2, t)
    >>> r[0]
    s/(4*sqrt(pi))
    >>> r[1]
    meijerg(((), ()), ((-1/2, 0), ()), s**2*t/4)
    >>> r[2]
    meijerg(((), ()), ((m/2,), (-m/2,)), t/4)

    """
    from ..core.numbers import ilcm

    def pb(g):
        a, b = _get_coeff_exp(g.argument, x)
        per = g.get_period()
        return meijerg(g.an, g.aother, g.bm, g.bother,
                       _my_principal_branch(a, per, full_pb)*x**b)

    _, s = _get_coeff_exp(po, x)
    _, b1 = _get_coeff_exp(g1.argument, x)
    _, b2 = _get_coeff_exp(g2.argument, x)
    if b1.is_negative:
        b1 = -b1
        g1 = _flip_g(g1)
    if b2.is_negative:
        b2 = -b2
        g2 = _flip_g(g2)
    if not b1.is_Rational or not b2.is_Rational:
        return
    m1, n1 = b1.numerator, b1.denominator
    m2, n2 = b2.numerator, b2.denominator
    tau = ilcm(m1*n2, m2*n1)
    r1 = tau//(m1*n2)
    r2 = tau//(m2*n1)

    C1, g1 = _inflate_g(g1, r1)
    C2, g2 = _inflate_g(g2, r2)
    g1 = pb(g1)
    g2 = pb(g2)

    fac *= C1*C2
    a1, b = _get_coeff_exp(g1.argument, x)
    a2, _ = _get_coeff_exp(g2.argument, x)

    # arbitrarily tack on the x**s part to g1
    # TODO should we try both?
    exp = (s + 1)/b - 1
    fac = fac/(abs(b) * a1**exp)

    def tr(l):
        return [a + exp for a in l]
    g1 = meijerg(tr(g1.an), tr(g1.aother), tr(g1.bm), tr(g1.bother), a1*x)
    g2 = meijerg(g2.an, g2.aother, g2.bm, g2.bother, a2*x)

    return powdenest(fac, polar=True), g1, g2


def _check_antecedents(g1, g2, x):
    """Return a condition under which the integral theorem applies."""
    from ..functions import (re, cos, exp, sin, sign, unpolarify,
                             arg as arg_, unbranched_argument as arg)
    #  Yes, this is madness.
    # XXX TODO this is a testing *nightmare*
    # NOTE if you update these conditions, please update the documentation as well

    # The following conditions are found in
    # [P], Section 2.24.1
    #
    # They are also reproduced (verbatim!) at
    # http://functions.wolfram.com/HypergeometricFunctions/MeijerG/21/02/03/
    #
    # Note: k=l=r=alpha=1
    sigma, _ = _get_coeff_exp(g1.argument, x)
    omega, _ = _get_coeff_exp(g2.argument, x)
    s, t, u, v = sympify([len(g1.bm), len(g1.an), len(g1.ap), len(g1.bq)])
    m, n, p, q = sympify([len(g2.bm), len(g2.an), len(g2.ap), len(g2.bq)])
    bstar = s + t - (u + v)/2
    cstar = m + n - (p + q)/2
    rho = g1.nu + (u - v)/2 + 1
    mu = g2.nu + (p - q)/2 + 1
    phi = q - p - (v - u)
    eta = 1 - (v - u) - mu - rho
    psi = (pi*(q - m - n) + abs(arg(omega)))/(q - p)
    theta = (pi*(v - s - t) + abs(arg(sigma)))/(v - u)

    _debug('Checking antecedents:')
    _debug('  sigma=%s, s=%s, t=%s, u=%s, v=%s, b*=%s, rho=%s'
           % (sigma, s, t, u, v, bstar, rho))
    _debug('  omega=%s, m=%s, n=%s, p=%s, q=%s, c*=%s, mu=%s,'
           % (omega, m, n, p, q, cstar, mu))
    _debug('  phi=%s, eta=%s, psi=%s, theta=%s' % (phi, eta, psi, theta))

    def _c1():
        for g in [g1, g2]:
            for i in g.an:
                for j in g.bm:
                    diff = i - j
                    if diff.is_integer and diff.is_positive:
                        return False
        return True
    c1 = _c1()
    c2 = And(*[re(1 + i + j) > 0 for i in g1.bm for j in g2.bm])
    c3 = And(*[re(1 + i + j) < 1 + 1 for i in g1.an for j in g2.an])
    c4 = And(*[(p - q)*re(1 + i - 1) - re(mu) > -Rational(3, 2) for i in g1.an])
    c5 = And(*[(p - q)*re(1 + i) - re(mu) > -Rational(3, 2) for i in g1.bm])
    c6 = And(*[(u - v)*re(1 + i - 1) - re(rho) > -Rational(3, 2) for i in g2.an])
    c7 = And(*[(u - v)*re(1 + i) - re(rho) > -Rational(3, 2) for i in g2.bm])
    c8 = (abs(phi) + 2*re((rho - 1)*(q - p) + (v - u)*(q - p) + (mu -
                                                                 1)*(v - u)) > 0)
    c9 = (abs(phi) - 2*re((rho - 1)*(q - p) + (v - u)*(q - p) + (mu -
                                                                 1)*(v - u)) > 0)
    c10 = (abs(arg(sigma)) < bstar*pi)
    c11 = Eq(abs(arg(sigma)), bstar*pi)
    c12 = (abs(arg(omega)) < cstar*pi)
    c13 = Eq(abs(arg(omega)), cstar*pi)

    # The following condition is *not* implemented as stated on the wolfram
    # function site. In the book of Prudnikov there is an additional part
    # (the And involving re()). However, I only have this book in russian, and
    # I don't read any russian. The following condition is what other people
    # have told me it means.
    # Worryingly, it is different from the condition implemented in REDUCE.
    # The REDUCE implementation:
    #   https://reduce-algebra.svn.sourceforge.net/svnroot/reduce-algebra/trunk/packages/defint/definta.red
    #   (search for tst14)
    # The Wolfram alpha version:
    #   http://functions.wolfram.com/HypergeometricFunctions/MeijerG/21/02/03/03/0014/
    z0 = exp(-(bstar + cstar)*pi*I)
    zos = unpolarify(z0*omega/sigma)
    zso = unpolarify(z0*sigma/omega)
    if zos == 1/zso:
        c14 = And(Eq(phi, 0), bstar + cstar <= 1,
                  Or(Ne(zos, 1), re(mu + rho + v - u) < 1,
                     re(mu + rho + q - p) < 1))
    else:
        c14 = And(Eq(phi, 0), bstar - 1 + cstar <= 0,
                  Or(And(Ne(zos, 1), abs(arg_(1 - zos)) < pi),
                     And(re(mu + rho + v - u) < 1, Eq(zos, 1))))

        def _cond():
            """
            Note: if `zso` is 1 then tmp will be NaN.  This raises a
            TypeError on `NaN < pi`.  Previously this gave `False` so
            this behavior has been hardcoded here but someone should
            check if this NaN is more serious! This NaN is triggered by
            test_meijerint() in test_meijerint.py:
            `meijerint_definite(exp(x), x, 0, I)`

            """
            tmp = abs(arg_(1 - zso))
            return False if tmp is nan else tmp < pi
        c14_alt = And(Eq(phi, 0), cstar - 1 + bstar <= 0,
                      Or(And(Ne(zso, 1), _cond()),
                         And(re(mu + rho + q - p) < 1, Eq(zso, 1))))

        # Since r=k=l=1, in our case there is c14_alt which is the same as calling
        # us with (g1, g2) = (g2, g1). The conditions below enumerate all cases
        # (i.e. we don't have to try arguments reversed by hand), and indeed try
        # all symmetric cases. (i.e. whenever there is a condition involving c14,
        # there is also a dual condition which is exactly what we would get when g1,
        # g2 were interchanged, *but c14 was unaltered*).
        # Hence the following seems correct:
        c14 = Or(c14, c14_alt)

    '''
    When `c15` is NaN (e.g. from `psi` being NaN as happens during
    'test_sympyissue_4992' and/or `theta` is NaN as in 'test_sympyissue_6253',
    both in `test_integrals.py`) the comparison to 0 formerly gave False
    whereas now an error is raised. To keep the old behavior, the value
    of NaN is replaced with False but perhaps a closer look at this condition
    should be made: XXX how should conditions leading to c15=NaN be handled?
    '''
    try:
        lambda_c = (q - p)*abs(omega)**(1/(q - p))*cos(psi) \
            + (v - u)*abs(sigma)**(1/(v - u))*cos(theta)
        # the TypeError might be raised here, e.g. if lambda_c is NaN
        if _eval_cond(lambda_c > 0) != false:
            c15 = (lambda_c > 0)
        else:
            def lambda_s0(c1, c2):
                return c1*(q - p)*abs(omega)**(1/(q - p))*sin(psi) \
                    + c2*(v - u)*abs(sigma)**(1/(v - u))*sin(theta)
            lambda_s = Piecewise(
                ((lambda_s0(+1, +1)*lambda_s0(-1, -1)),
                 And(Eq(arg(sigma), 0), Eq(arg(omega), 0))),
                (lambda_s0(sign(arg(omega)), +1)*lambda_s0(sign(arg(omega)), -1),
                 And(Eq(arg(sigma), 0), Ne(arg(omega), 0))),
                (lambda_s0(+1, sign(arg(sigma)))*lambda_s0(-1, sign(arg(sigma))),
                 And(Ne(arg(sigma), 0), Eq(arg(omega), 0))),
                (lambda_s0(sign(arg(omega)), sign(arg(sigma))), True))
            tmp = [lambda_c > 0,
                   And(Eq(lambda_c, 0), Ne(lambda_s, 0), re(eta) > -1),
                   And(Eq(lambda_c, 0), Eq(lambda_s, 0), re(eta) > 0)]
            c15 = Or(*tmp)
    except TypeError:
        c15 = False
    for cond, i in [(c1, 1), (c2, 2), (c3, 3), (c4, 4), (c5, 5), (c6, 6),
                    (c7, 7), (c8, 8), (c9, 9), (c10, 10), (c11, 11),
                    (c12, 12), (c13, 13), (c14, 14), (c15, 15)]:
        _debug('  c%s:' % i, cond)

    # We will return Or(*conds)
    conds = []

    def pr(count):
        _debug('  case %s:' % count, conds[-1])
    conds += [And(m*n*s*t != 0, bstar.is_positive is True, cstar.is_positive is True, c1, c2, c3, c10,
                  c12)]  # 1
    pr(1)
    conds += [And(Eq(u, v), Eq(bstar, 0), cstar.is_positive is True, sigma.is_positive is True, re(rho) < 1,
                  c1, c2, c3, c12)]  # 2
    pr(2)
    conds += [And(Eq(p, q), Eq(cstar, 0), bstar.is_positive is True, omega.is_positive is True, re(mu) < 1,
                  c1, c2, c3, c10)]  # 3
    pr(3)
    conds += [And(Eq(p, q), Eq(u, v), Eq(bstar, 0), Eq(cstar, 0),
                  sigma.is_positive is True, omega.is_positive is True, re(mu) < 1, re(rho) < 1,
                  Ne(sigma, omega), c1, c2, c3)]  # 4
    pr(4)
    conds += [And(Eq(p, q), Eq(u, v), Eq(bstar, 0), Eq(cstar, 0),
                  sigma.is_positive is True, omega.is_positive is True, re(mu + rho) < 1,
                  Ne(omega, sigma), c1, c2, c3)]  # 5
    pr(5)
    conds += [And(p > q, s.is_positive is True, bstar.is_positive is True, cstar >= 0,
                  c1, c2, c3, c5, c10, c13)]  # 6
    pr(6)
    conds += [And(p < q, t.is_positive is True, bstar.is_positive is True, cstar >= 0,
                  c1, c2, c3, c4, c10, c13)]  # 7
    pr(7)
    conds += [And(u > v, m.is_positive is True, cstar.is_positive is True, bstar >= 0,
                  c1, c2, c3, c7, c11, c12)]  # 8
    pr(8)
    conds += [And(u < v, n.is_positive is True, cstar.is_positive is True, bstar >= 0,
                  c1, c2, c3, c6, c11, c12)]  # 9
    pr(9)
    conds += [And(p > q, Eq(u, v), Eq(bstar, 0), cstar >= 0, sigma.is_positive is True,
                  re(rho) < 1, c1, c2, c3, c5, c13)]  # 10
    pr(10)
    conds += [And(p < q, Eq(u, v), Eq(bstar, 0), cstar >= 0, sigma.is_positive is True,
                  re(rho) < 1, c1, c2, c3, c4, c13)]  # 11
    pr(11)
    conds += [And(Eq(p, q), u > v, bstar >= 0, Eq(cstar, 0), omega.is_positive is True,
                  re(mu) < 1, c1, c2, c3, c7, c11)]  # 12
    pr(12)
    conds += [And(Eq(p, q), u < v, bstar >= 0, Eq(cstar, 0), omega.is_positive is True,
                  re(mu) < 1, c1, c2, c3, c6, c11)]  # 13
    pr(13)
    conds += [And(p < q, u > v, bstar >= 0, cstar >= 0,
                  c1, c2, c3, c4, c7, c11, c13)]  # 14
    pr(14)
    conds += [And(p > q, u < v, bstar >= 0, cstar >= 0,
                  c1, c2, c3, c5, c6, c11, c13)]  # 15
    pr(15)
    conds += [And(p > q, u > v, bstar >= 0, cstar >= 0,
                  c1, c2, c3, c5, c7, c8, c11, c13, c14)]  # 16
    pr(16)
    conds += [And(p < q, u < v, bstar >= 0, cstar >= 0,
                  c1, c2, c3, c4, c6, c9, c11, c13, c14)]  # 17
    pr(17)
    conds += [And(Eq(t, 0), s.is_positive is True, bstar.is_positive is True, phi.is_positive is True, c1, c2, c10)]  # 18
    pr(18)
    conds += [And(Eq(s, 0), t.is_positive is True, bstar.is_positive is True, phi.is_negative is True, c1, c3, c10)]  # 19
    pr(19)
    conds += [And(Eq(n, 0), m.is_positive is True, cstar.is_positive is True, phi.is_negative is True, c1, c2, c12)]  # 20
    pr(20)
    conds += [And(Eq(m, 0), n.is_positive is True, cstar.is_positive is True, phi.is_positive is True, c1, c3, c12)]  # 21
    pr(21)
    conds += [And(Eq(s*t, 0), bstar.is_positive is True, cstar.is_positive is True,
                  c1, c2, c3, c10, c12)]  # 22
    pr(22)
    conds += [And(Eq(m*n, 0), bstar.is_positive is True, cstar.is_positive is True,
                  c1, c2, c3, c10, c12)]  # 23
    pr(23)

    # The following case is from [Luke1969]. As far as I can tell, it is *not*
    # covered by Prudnikov's.
    # Let G1 and G2 be the two G-functions. Suppose the integral exists from
    # 0 to a > 0 (this is easy the easy part), that G1 is exponential decay at
    # infinity, and that the mellin transform of G2 exists.
    # Then the integral exists.
    mt1_exists = _check_antecedents_1(g1, x, helper=True)
    mt2_exists = _check_antecedents_1(g2, x, helper=True)
    conds += [And(mt2_exists, Eq(t, 0), u < s, bstar.is_positive is True, c10, c1, c2, c3)]
    pr('E1')
    conds += [And(mt2_exists, Eq(s, 0), v < t, bstar.is_positive is True, c10, c1, c2, c3)]
    pr('E2')
    conds += [And(mt1_exists, Eq(n, 0), p < m, cstar.is_positive is True, c12, c1, c2, c3)]
    pr('E3')
    conds += [And(mt1_exists, Eq(m, 0), q < n, cstar.is_positive is True, c12, c1, c2, c3)]
    pr('E4')

    # Let's short-circuit if this worked ...
    # the rest is corner-cases and terrible to read.
    r = Or(*conds)
    if _eval_cond(r) != false:
        return r

    conds += [And(m + n > p, Eq(t, 0), Eq(phi, 0), s.is_positive is True, bstar.is_positive is True, cstar.is_negative is True,
                  abs(arg(omega)) < (m + n - p + 1)*pi,
                  c1, c2, c10, c14, c15)]  # 24
    pr(24)
    conds += [And(m + n > q, Eq(s, 0), Eq(phi, 0), t.is_positive is True, bstar.is_positive is True, cstar.is_negative is True,
                  abs(arg(omega)) < (m + n - q + 1)*pi,
                  c1, c3, c10, c14, c15)]  # 25
    pr(25)
    conds += [And(Eq(p, q - 1), Eq(t, 0), Eq(phi, 0), s.is_positive is True, bstar.is_positive is True,
                  cstar >= 0, cstar*pi < abs(arg(omega)),
                  c1, c2, c10, c14, c15)]  # 26
    pr(26)
    conds += [And(Eq(p, q + 1), Eq(s, 0), Eq(phi, 0), t.is_positive is True, bstar.is_positive is True,
                  cstar >= 0, cstar*pi < abs(arg(omega)),
                  c1, c3, c10, c14, c15)]  # 27
    pr(27)
    conds += [And(p < q - 1, Eq(t, 0), Eq(phi, 0), s.is_positive is True, bstar.is_positive is True,
                  cstar >= 0, cstar*pi < abs(arg(omega)),
                  abs(arg(omega)) < (m + n - p + 1)*pi,
                  c1, c2, c10, c14, c15)]  # 28
    pr(28)
    conds += [And(
        p > q + 1, Eq(s, 0), Eq(phi, 0), t.is_positive is True, bstar.is_positive is True, cstar >= 0,
        cstar*pi < abs(arg(omega)),
        abs(arg(omega)) < (m + n - q + 1)*pi,
        c1, c3, c10, c14, c15)]  # 29
    pr(29)
    conds += [And(Eq(n, 0), Eq(phi, 0), s + t > 0, m.is_positive is True, cstar.is_positive is True, bstar.is_negative is True,
                  abs(arg(sigma)) < (s + t - u + 1)*pi,
                  c1, c2, c12, c14, c15)]  # 30
    pr(30)
    conds += [And(Eq(m, 0), Eq(phi, 0), s + t > v, n.is_positive is True, cstar.is_positive is True, bstar.is_negative is True,
                  abs(arg(sigma)) < (s + t - v + 1)*pi,
                  c1, c3, c12, c14, c15)]  # 31
    pr(31)
    conds += [And(Eq(n, 0), Eq(phi, 0), Eq(u, v - 1), m.is_positive is True, cstar.is_positive is True,
                  bstar >= 0, bstar*pi < abs(arg(sigma)),
                  abs(arg(sigma)) < (bstar + 1)*pi,
                  c1, c2, c12, c14, c15)]  # 32
    pr(32)
    conds += [And(Eq(m, 0), Eq(phi, 0), Eq(u, v + 1), n.is_positive is True, cstar.is_positive is True,
                  bstar >= 0, bstar*pi < abs(arg(sigma)),
                  abs(arg(sigma)) < (bstar + 1)*pi,
                  c1, c3, c12, c14, c15)]  # 33
    pr(33)
    conds += [And(
        Eq(n, 0), Eq(phi, 0), u < v - 1, m.is_positive is True, cstar.is_positive is True, bstar >= 0,
        bstar*pi < abs(arg(sigma)),
        abs(arg(sigma)) < (s + t - u + 1)*pi,
        c1, c2, c12, c14, c15)]  # 34
    pr(34)
    conds += [And(
        Eq(m, 0), Eq(phi, 0), u > v + 1, n.is_positive is True, cstar.is_positive is True, bstar >= 0,
        bstar*pi < abs(arg(sigma)),
        abs(arg(sigma)) < (s + t - v + 1)*pi,
        c1, c3, c12, c14, c15)]  # 35
    pr(35)

    return Or(*conds)

    # NOTE An alternative, but as far as I can tell weaker, set of conditions
    #      can be found in [L, section 5.6.2].


def _int0oo(g1, g2, x):
    """
    Express integral from zero to infinity g1*g2 using a G function,
    assuming the necessary conditions are fulfilled.

    >>> from diofant.abc import s
    >>> g1 = meijerg([], [], [-Rational(1, 2), 0], [], s**2*t/4)
    >>> g2 = meijerg([], [], [m/2], [-m/2], t/4)
    >>> _int0oo(g1, g2, t)
    4*meijerg(((1/2, 0), ()), ((m/2,), (-m/2,)), s**(-2))/s**2

    """
    # See: [L, section 5.6.2, equation (1)]
    eta, _ = _get_coeff_exp(g1.argument, x)
    omega, _ = _get_coeff_exp(g2.argument, x)

    def neg(l):
        return [-x for x in l]
    a1 = neg(g1.bm) + list(g2.an)
    a2 = list(g2.aother) + neg(g1.bother)
    b1 = neg(g1.an) + list(g2.bm)
    b2 = list(g2.bother) + neg(g1.aother)
    return meijerg(a1, a2, b1, b2, omega/eta)/eta


def _rewrite_inversion(fac, po, g, x):
    """Absorb ``po`` == x**s into g."""
    _, s = _get_coeff_exp(po, x)
    a, b = _get_coeff_exp(g.argument, x)

    def tr(l):
        return [t + s/b for t in l]
    return (powdenest(fac/a**(s/b), polar=True),
            meijerg(tr(g.an), tr(g.aother), tr(g.bm), tr(g.bother), g.argument))


def _check_antecedents_inversion(g, x):
    """Check antecedents for the laplace inversion integral."""
    from ..functions import re, im, exp
    _debug('Checking antecedents for inversion:')
    z = g.argument
    _, e = _get_coeff_exp(z, x)
    if e < 0:
        _debug('  Flipping G.')
        # We want to assume that argument gets large as |x| -> oo
        return _check_antecedents_inversion(_flip_g(g), x)

    def statement_half(a, b, c, z, plus):
        coeff, exponent = _get_coeff_exp(z, x)
        a *= exponent
        b *= coeff**c
        c *= exponent
        conds = []
        wp = b*exp(I*re(c)*pi/2)
        wm = b*exp(-I*re(c)*pi/2)
        if plus:
            w = wp
        else:
            w = wm
        conds += [And(Or(Eq(b, 0), re(c) <= 0), re(a) <= -1)]
        conds += [And(Ne(b, 0), Eq(im(c), 0), re(c) > 0, re(w) < 0)]
        conds += [And(Ne(b, 0), Eq(im(c), 0), re(c) > 0, re(w) <= 0,
                      re(a) <= -1)]
        return Or(*conds)

    def statement(a, b, c, z):
        """ Provide a convergence statement for z**a * exp(b*z**c),
        c/f sphinx docs.

        """
        return And(statement_half(a, b, c, z, True),
                   statement_half(a, b, c, z, False))

    # Notations from [L], section 5.7-10
    m, n, p, q = sympify([len(g.bm), len(g.an), len(g.ap), len(g.bq)])
    tau = m + n - p
    nu = q - m - n
    rho = (tau - nu)/2
    sigma = q - p
    if sigma == 1:
        epsilon = Rational(1, 2)
    elif sigma > 1:
        epsilon = 1
    else:
        epsilon = nan
    theta = ((1 - sigma)/2 + Add(*g.bq) - Add(*g.ap))/sigma
    delta = g.delta
    _debug('  m=%s, n=%s, p=%s, q=%s, tau=%s, nu=%s, rho=%s, sigma=%s' % (
        m, n, p, q, tau, nu, rho, sigma))
    _debug('  epsilon=%s, theta=%s, delta=%s' % (epsilon, theta, delta))

    # First check if the computation is valid.
    if not (g.delta >= e/2 or (p >= 1 and p >= q)):
        _debug('  Computation not valid for these parameters.')
        return False

    # Now check if the inversion integral exists.

    # Test "condition A"
    for a in g.an:
        for b in g.bm:
            if (a - b).is_integer and a > b:
                _debug('  Not a valid G function.')
                return False

    # There are two cases. If p >= q, we can directly use a slater expansion
    # like [L], 5.2 (11). Note in particular that the asymptotics of such an
    # expansion even hold when some of the parameters differ by integers, i.e.
    # the formula itself would not be valid! (b/c G functions are cts. in their
    # parameters)
    # When p < q, we need to use the theorems of [L], 5.10.

    if p >= q:
        _debug('  Using asymptotic Slater expansion.')
        return And(*[statement(a - 1, 0, 0, z) for a in g.an])

    def E(z):
        return And(*[statement(a - 1, 0, z) for a in g.an])

    def H(z):
        return statement(theta, -sigma, 1/sigma, z)

    def Hp(z):
        return statement_half(theta, -sigma, 1/sigma, z, True)

    def Hm(z):
        return statement_half(theta, -sigma, 1/sigma, z, False)

    # [L], section 5.10
    conds = []
    # Theorem 1
    conds += [And(1 <= n, p < q, 1 <= m, rho*pi - delta >= pi/2, delta > 0,
                  E(z*exp(I*pi*(nu + 1))))]
    # Theorem 2, statements (2) and (3)
    conds += [And(p + 1 <= m, m + 1 <= q, delta > 0, delta < pi/2, n == 0,
                  (m - p + 1)*pi - delta >= pi/2,
                  Hp(z*exp(I*pi*(q - m))), Hm(z*exp(-I*pi*(q - m))))]
    # Theorem 2, statement (5)
    conds += [And(p < q, m == q, n == 0, delta > 0,
                  (sigma + epsilon)*pi - delta >= pi/2, H(z))]
    # Theorem 3, statements (6) and (7)
    conds += [And(Or(And(p <= q - 2, 1 <= tau, tau <= sigma/2),
                     And(p + 1 <= m + n, m + n <= (p + q)/2)),
                  delta > 0, delta < pi/2, (tau + 1)*pi - delta >= pi/2,
                  Hp(z*exp(I*pi*nu)), Hm(z*exp(-I*pi*nu)))]
    # Theorem 4, statements (10) and (11)
    conds += [And(p < q, 1 <= m, rho > 0, delta > 0, delta + rho*pi < pi/2,
                  (tau + epsilon)*pi - delta >= pi/2,
                  Hp(z*exp(I*pi*nu)), Hm(z*exp(-I*pi*nu)))]
    # Trivial case
    conds += [m == 0]

    # TODO
    # Theorem 5 is quite general
    # Theorem 6 contains special cases for q=p+1

    return Or(*conds)


def _int_inversion(g, x, t):
    """Compute the laplace inversion integral, assuming the formula applies."""
    b, a = _get_coeff_exp(g.argument, x)
    C, g = _inflate_fox_h(meijerg(g.an, g.aother, g.bm, g.bother, b/t**a), -a)
    return C/t*g


####################################################################
# Finally, the real meat.
####################################################################

_lookup_table = None


@cacheit
def _rewrite_single(f, x, recursive=True):
    """
    Try to rewrite f as a sum of single G functions of the form
    C*x**s*G(a*x**b), where b is a rational number and C is independent of x.
    We guarantee that result.argument.as_coeff_mul(x) returns (a, (x**b,))
    or (a, ()).
    Returns a list of tuples (C, s, G) and a condition cond.
    Returns None on failure.

    """
    from ..functions import polarify, unpolarify
    global _lookup_table
    if not _lookup_table:
        _lookup_table = defaultdict(list)
        _create_lookup_table(_lookup_table)

    if isinstance(f, meijerg):
        from ..polys import factor
        coeff, m = factor(f.argument, x).as_coeff_mul(x)
        if len(m) > 1:
            return
        m = m[0]
        if m.is_Pow:
            if m.base != x or not m.exp.is_Rational:
                return
        elif m != x:
            return
        return [(1, 0, meijerg(f.an, f.aother, f.bm, f.bother, coeff*m))], True

    f_ = f
    f = f.subs({x: z})
    t = _mytype(f, z)
    if t in _lookup_table:
        l = _lookup_table[t]
        for formula, terms, cond, hint in l:
            subs = f.match(formula)
            if subs:
                subs_ = {}
                for fro, to in subs.items():
                    subs_[fro] = unpolarify(polarify(to, lift=True),
                                            exponents_only=True)
                subs = subs_
                if not isinstance(hint, bool):
                    hint = hint.subs(subs)
                if hint == false:
                    continue
                if not isinstance(cond, (bool, BooleanAtom)):
                    cond = unpolarify(cond.subs(subs))
                if _eval_cond(cond) == false:
                    continue
                if not isinstance(terms, list):
                    terms = terms(subs)
                res = []
                for fac, g in terms:
                    r1 = _get_coeff_exp(unpolarify(fac.subs(subs).subs({z: x}),
                                                   exponents_only=True), x)
                    try:
                        g = g.subs(subs).subs({z: x})
                    except ValueError:
                        continue

                    # NOTE these substitutions can in principle introduce oo,
                    #      zoo and other absurdities. It shouldn't matter,
                    #      but better be safe.
                    if Tuple(*(r1 + (g,))).has(oo, zoo, -oo):
                        continue
                    g = meijerg(g.an, g.aother, g.bm, g.bother,
                                unpolarify(g.argument, exponents_only=True))
                    res.append(r1 + (g,))
                if res:
                    return res, cond

    # try recursive mellin transform
    if not recursive:
        return
    _debug('Trying recursive Mellin transform method.')
    from .transforms import (mellin_transform, inverse_mellin_transform,
                             IntegralTransformError, MellinTransformStripError)
    from ..simplify import simplify
    from ..polys import cancel

    def my_imt(F, s, x, strip):
        """ Calling simplify() all the time is slow and not helpful, since
        most of the time it only factors things in a way that has to be
        un-done anyway. But sometimes it can remove apparent poles.

        """
        # XXX should this be in inverse_mellin_transform?
        try:
            return inverse_mellin_transform(F, s, x, strip,
                                            as_meijerg=True, needeval=True)
        except MellinTransformStripError:
            return inverse_mellin_transform(
                simplify(cancel(expand(F))), s, x, strip,
                as_meijerg=True, needeval=True)
    f = f_
    s = _dummy('s', 'rewrite-single', f)
    # to avoid infinite recursion, we have to force the two g functions case

    def my_integrator(f, x):
        from .integrals import Integral
        from ..simplify import hyperexpand
        r = _meijerint_definite_4(f, x, only_double=True)
        if r is not None:
            res, cond = r
            res = _my_unpolarify(hyperexpand(res, rewrite='nonrepsmall'))
            return Piecewise((res, cond),
                             (Integral(f, (x, 0, oo)), True))
        return Integral(f, (x, 0, oo))
    try:
        F, strip, _ = mellin_transform(f, x, s, integrator=my_integrator,
                                       simplify=False, needeval=True)
        g = my_imt(F, s, x, strip)
    except IntegralTransformError:
        g = None
    if g is None:
        # We try to find an expression by analytic continuation.
        # (also if the dummy is already in the expression, there is no point in
        #  putting in another one)
        a = _dummy_('a', 'rewrite-single')
        if a not in f.free_symbols and _is_analytic(f, x):
            try:
                F, strip, _ = mellin_transform(f.subs({x: a*x}), x, s,
                                               integrator=my_integrator,
                                               needeval=True, simplify=False)
                g = my_imt(F, s, x, strip).subs({a: 1})
            except IntegralTransformError:
                g = None
    if g is None or g.has(oo, nan, zoo):
        _debug('Recursive Mellin transform failed.')
        return
    args = Add.make_args(g)
    res = []
    for f in args:
        c, m = f.as_coeff_mul(x)
        if len(m) > 1:
            raise NotImplementedError('Unexpected form...')
        g = m[0]
        a, b = _get_coeff_exp(g.argument, x)
        res += [(c, 0, meijerg(g.an, g.aother, g.bm, g.bother,
                               unpolarify(polarify(
                                   a, lift=True), exponents_only=True)
                               * x**b))]
    _debug('Recursive Mellin transform worked:', g)
    return res, True


def _rewrite1(f, x, recursive=True):
    """
    Try to rewrite f using a (sum of) single G functions with argument a*x**b.
    Return fac, po, g such that f = fac*po*g, fac is independent of x
    and po = x**s.
    Here g is a result from _rewrite_single.
    Return None on failure.

    """
    fac, po, g = _split_mul(f, x)
    g = _rewrite_single(g, x, recursive)
    if g:
        return fac, po, g[0], g[1]


def _rewrite2(f, x):
    """
    Try to rewrite f as a product of two G functions of arguments a*x**b.
    Return fac, po, g1, g2 such that f = fac*po*g1*g2, where fac is
    independent of x and po is x**s.
    Here g1 and g2 are results of _rewrite_single.
    Returns None on failure.

    """
    fac, po, g = _split_mul(f, x)
    if any(_rewrite_single(expr, x, False) is None for expr in _mul_args(g)):
        return
    l = _mul_as_two_parts(g)
    if not l:
        return
    l = list(ordered(l, [
        lambda p: max(len(_exponents(p[0], x)), len(_exponents(p[1], x))),
        lambda p: max(len(_functions(p[0], x)), len(_functions(p[1], x))),
        lambda p: max(len(_find_splitting_points(p[0], x)),
                      len(_find_splitting_points(p[1], x)))]))

    for recursive in [False, True]:
        for fac1, fac2 in l:
            g1 = _rewrite_single(fac1, x, recursive)
            g2 = _rewrite_single(fac2, x, recursive)
            if g1 and g2:
                cond = And(g1[1], g2[1])
                if cond != false:
                    return fac, po, g1[0], g2[0], cond


def meijerint_indefinite(f, x):
    """
    Compute an indefinite integral of ``f`` by rewriting it as a G function.

    Examples
    ========

    >>> meijerint_indefinite(sin(x), x)
    -cos(x)

    """
    from ..functions import hyper, meijerg

    results = []
    for a in sorted(_find_splitting_points(f, x) | {Integer(0)}, key=default_sort_key):
        res = _meijerint_indefinite_1(f.subs({x: x + a}), x)
        if not res:
            continue
        res = res.subs({x: x - a})
        if _has(res, hyper, meijerg):
            results.append(res)
        else:
            return res
    if f.has(HyperbolicFunction):
        _debug('Try rewriting hyperbolics in terms of exp.')
        rv = meijerint_indefinite(
            _rewrite_hyperbolics_as_exp(f), x)
        if rv:
            if not type(rv) is list:
                return collect(factor_terms(rv),
                               {a for a in rv.atoms(Pow) if a.base is E})
            results.extend(rv)
    if results:
        return next(ordered(results))


def _meijerint_indefinite_1(f, x):
    """Helper that does not attempt any substitution."""
    from .integrals import Integral
    from ..functions import piecewise_fold
    _debug('Trying to compute the indefinite integral of', f, 'wrt', x)

    gs = _rewrite1(f, x)
    if gs is None:
        # Note: the code that calls us will do expand() and try again
        return

    fac, po, gl, cond = gs
    _debug(' could rewrite:', gs)
    res = Integer(0)
    for C, s, g in gl:
        a, b = _get_coeff_exp(g.argument, x)
        _, c = _get_coeff_exp(po, x)
        c += s

        # we do a substitution t=a*x**b, get integrand fac*t**rho*g
        fac_ = fac * C / (b*a**((1 + c)/b))
        rho = (c + 1)/b - 1

        # we now use t**rho*G(params, t) = G(params + rho, t)
        # [L, page 150, equation (4)]
        # and integral G(params, t) dt = G(1, params+1, 0, t)
        #   (or a similar expression with 1 and 0 exchanged ... pick the one
        #    which yields a well-defined function)
        # [R, section 5]
        # (Note that this dummy will immediately go away again, so we
        #  can safely pass Integer(1) for ``expr``.)
        t = _dummy('t', 'meijerint-indefinite', Integer(1))

        def tr(p):
            return [a + rho + 1 for a in p]
        if any(b.is_integer and b.is_nonpositive for b in tr(g.bm)):
            r = -meijerg(
                tr(g.an), tr(g.aother) + [1], tr(g.bm) + [0], tr(g.bother), t)
        else:
            r = meijerg(
                tr(g.an) + [1], tr(g.aother), tr(g.bm), tr(g.bother) + [0], t)

        # The antiderivative is most often expected to be defined
        # in the neighborhood of  x = 0.
        place = 0
        if b < 0 or f.subs({x: 0}).has(nan, zoo):
            place = None
        r = hyperexpand(r.subs({t: a*x**b}), place=place)

        # now substitute back
        # Note: we really do want the powers of x to combine.
        res += powdenest(fac_*r, polar=True)

    def _clean(res):
        """This multiplies out superfluous powers of x we created, and chops off
        constants:

            >> _clean(x*(exp(x)/x - 1/x) + 3)
            exp(x)

        cancel is used before mul_expand since it is possible for an
        expression to have an additive constant that doesn't become isolated
        with simple expansion. Such a situation was identified in issue sympy/sympy#6369:


        >>> a = sqrt(2*x + 1)
        >>> bad = (3*x*a**5 + 2*x - a**5 + 1)/a**2
        >>> bad.expand().as_independent(x)[0]
        0
        >>> cancel(bad).expand().as_independent(x)[0]
        1

        """
        from ..polys import cancel
        res = expand_mul(cancel(res), deep=False)
        return Add._from_args(res.as_coeff_add(x)[1])

    res = piecewise_fold(res)
    if res.is_Piecewise:
        newargs = []
        for expr, cond in res.args:
            expr = _my_unpolarify(_clean(expr))
            newargs += [(expr, cond)]
        res = Piecewise(*newargs)
    else:
        res = _my_unpolarify(_clean(res))
    return Piecewise((res, _my_unpolarify(cond)), (Integral(f, x), True))


def meijerint_definite(f, x, a, b):
    """
    Integrate ``f`` over the interval [``a``, ``b``], by rewriting it as a product
    of two G functions, or as a single G function.

    Return res, cond, where cond are convergence conditions.

    Examples
    ========

    >>> meijerint_definite(exp(-x**2), x, -oo, oo)
    (sqrt(pi), true)

    This function is implemented as a succession of functions
    meijerint_definite, _meijerint_definite_2, _meijerint_definite_3,
    _meijerint_definite_4. Each function in the list calls the next one
    (presumably) several times. This means that calling meijerint_definite
    can be very costly.

    """
    # This consists of three steps:
    # 1) Change the integration limits to 0, oo
    # 2) Rewrite in terms of G functions
    # 3) Evaluate the integral
    #
    # There are usually several ways of doing this, and we want to try all.
    # This function does (1), calls _meijerint_definite_2 for step (2).
    from ..functions import arg, exp, DiracDelta
    _debug('Integrating', f, 'wrt %s from %s to %s.' % (x, a, b))

    if f.has(DiracDelta):
        _debug('Integrand has DiracDelta terms - giving up.')
        return

    f_, x_, a_, b_ = f, x, a, b

    # Let's use a dummy in case any of the boundaries has x.
    d = Dummy('x')
    f = f.subs({x: d})
    x = d

    if a == b:
        return Integer(0), True

    results = []
    if a == -oo and b != oo:
        return meijerint_definite(f.subs({x: -x}), x, -b, -a)

    elif a == -oo:
        # Integrating -oo to oo. We need to find a place to split the integral.
        _debug('  Integrating -oo to +oo.')
        innermost = _find_splitting_points(f, x)
        _debug('  Sensible splitting points:', innermost)
        for c in sorted(innermost, key=default_sort_key, reverse=True) + [Integer(0)]:
            _debug('  Trying to split at', c)
            if not c.is_extended_real:
                _debug('  Non-real splitting point.')
                continue
            res1 = _meijerint_definite_2(f.subs({x: x + c}), x)
            if res1 is None:
                _debug('  But could not compute first integral.')
                continue
            res2 = _meijerint_definite_2(f.subs({x: c - x}), x)
            if res2 is None:
                _debug('  But could not compute second integral.')
                continue
            res1, cond1 = res1
            res2, cond2 = res2
            cond = _condsimp(And(cond1, cond2))
            if cond == false:
                _debug('  But combined condition is always false.')
                continue
            res = res1 + res2
            return res, cond

    elif a == oo:
        res = meijerint_definite(f, x, b, oo)
        return -res[0], res[1]

    elif (a, b) == (0, oo):
        # This is a common case - try it directly first.
        res = _meijerint_definite_2(f, x)
        if res:
            if _has(res[0], meijerg):
                results.append(res)
            else:
                return res

    else:
        if b == oo:
            for split in _find_splitting_points(f, x):
                if (a - split).is_nonnegative:
                    _debug('Trying x -> x + %s' % split)
                    res = _meijerint_definite_2(f.subs({x: x + split})
                                                * Heaviside(x + split - a), x)
                    if res:
                        if _has(res[0], meijerg):
                            results.append(res)
                        else:
                            return res

        f = f.subs({x: x + a})
        b = b - a
        a = 0
        if b != oo:
            phi = exp(I*arg(b))
            b = abs(b)
            f = f.subs({x: phi*x})
            f *= Heaviside(b - x)*phi
            b = oo

        _debug('Changed limits to', a, b)
        _debug('Changed function to', f)
        res = _meijerint_definite_2(f, x)
        if res:
            if _has(res[0], meijerg):
                results.append(res)
            else:
                return res
    if f_.has(HyperbolicFunction):
        _debug('Try rewriting hyperbolics in terms of exp.')
        rv = meijerint_definite(
            _rewrite_hyperbolics_as_exp(f_), x_, a_, b_)
        if rv:
            if not type(rv) is list:
                rv = (collect(factor_terms(rv[0]),
                              {a for a in rv[0].atoms(Pow) if a.base is E}),) + rv[1:]
                return rv
            results.extend(rv)
    if results:
        return next(ordered(results))


def _guess_expansion(f, x):
    """Try to guess sensible rewritings for integrand f(x)."""
    from ..core import expand_trig
    from ..functions.elementary.trigonometric import TrigonometricFunction
    res = [(f, 'original integrand')]

    orig = res[-1][0]
    saw = {orig}
    expanded = expand_mul(orig)
    if expanded not in saw:
        res += [(expanded, 'expand_mul')]
        saw.add(expanded)

    expanded = expand(orig)
    if expanded not in saw:
        res += [(expanded, 'expand')]
        saw.add(expanded)

    if orig.has(TrigonometricFunction, HyperbolicFunction):
        expanded = expand_mul(expand_trig(orig))
        if expanded not in saw:
            res += [(expanded, 'expand_trig, expand_mul')]
            saw.add(expanded)

    if orig.has(cos, sin):
        reduced = sincos_to_sum(orig)
        if reduced not in saw:
            res += [(reduced, 'trig power reduction')]
            saw.add(reduced)

    return res


def _meijerint_definite_2(f, x):
    """
    Try to integrate f dx from zero to infinty.

    The body of this function computes various 'simplifications'
    f1, f2, ... of f (e.g. by calling expand_mul(), trigexpand()
    - see _guess_expansion) and calls _meijerint_definite_3 with each of
    these in succession.
    If _meijerint_definite_3 succeedes with any of the simplified functions,
    returns this result.

    """
    # This function does preparation for (2), calls
    # _meijerint_definite_3 for (2) and (3) combined.

    # use a positive dummy - we integrate from 0 to oo
    # XXX if a nonnegative symbol is used there will be test failures
    dummy = _dummy('x', 'meijerint-definite2', f, positive=True)
    f = f.subs({x: dummy})
    x = dummy

    if f == 0:
        return Integer(0), True

    for g, explanation in _guess_expansion(f, x):
        _debug('Trying', explanation)
        res = _meijerint_definite_3(g, x)
        if res:
            return res


def _meijerint_definite_3(f, x):
    """
    Try to integrate f dx from zero to infinity.

    This function calls _meijerint_definite_4 to try to compute the
    integral. If this fails, it tries using linearity.

    """
    res = _meijerint_definite_4(f, x)
    if res and res[1] != false:
        return res
    if f.is_Add:
        _debug('Expanding and evaluating all terms.')
        ress = [_meijerint_definite_4(g, x) for g in f.args]
        if all(r is not None for r in ress):
            conds = []
            res = Integer(0)
            for r, c in ress:
                res += r
                conds += [c]
            c = And(*conds)
            if c != false:
                return res, c


def _my_unpolarify(f):
    from ..functions import unpolarify
    return _eval_cond(unpolarify(f))


def _meijerint_definite_4(f, x, only_double=False):
    """
    Try to integrate f dx from zero to infinity.

    This function tries to apply the integration theorems found in literature,
    i.e. it tries to rewrite f as either one or a product of two G-functions.

    The parameter ``only_double`` is used internally in the recursive algorithm
    to disable trying to rewrite f as a single G-function.

    """
    # This function does (2) and (3)
    _debug('Integrating', f)
    # Try single G function.
    if not only_double:
        gs = _rewrite1(f, x, recursive=False)
        if gs is not None:
            fac, po, g, cond = gs
            _debug('Could rewrite as single G function:', fac, po, g)
            res = Integer(0)
            for C, s, f in g:
                if C == 0:
                    continue
                C, f = _rewrite_saxena_1(fac*C, po*x**s, f, x)
                res += C*_int0oo_1(f, x)
                cond = And(cond, _check_antecedents_1(f, x))
                if cond == false:
                    break
            cond = _my_unpolarify(cond)
            if cond == false:
                _debug('But cond is always False.')
            else:
                _debug('Result before branch substitutions is:', res)
                return _my_unpolarify(hyperexpand(res)), cond

    # Try two G functions.
    gs = _rewrite2(f, x)
    if gs is not None:
        for full_pb in [False, True]:
            fac, po, g1, g2, cond = gs
            _debug('Could rewrite as two G functions:', fac, po, g1, g2)
            res = Integer(0)
            for C1, s1, f1 in g1:
                for C2, s2, f2 in g2:
                    r = _rewrite_saxena(fac*C1*C2, po*x**(s1 + s2),
                                        f1, f2, x, full_pb)
                    if r is None:
                        _debug('Non-rational exponents.')
                        return
                    C, f1_, f2_ = r
                    _debug('Saxena subst for yielded:', C, f1_, f2_)
                    cond = And(cond, _check_antecedents(f1_, f2_, x))
                    if cond == false:
                        break
                    res += C*_int0oo(f1_, f2_, x)
                else:
                    continue
                break
            cond = _my_unpolarify(cond)
            if cond == false:
                _debug('But cond is always False (full_pb=%s).' % full_pb)
            else:
                _debug('Result before branch substitutions is:', res)
                if only_double:
                    return res, cond
                return _my_unpolarify(hyperexpand(res)), cond


def meijerint_inversion(f, x, t):
    r"""
    Compute the inverse laplace transform
    `\int_{c+i\infty}^{c-i\infty} f(x) e^{tx) dx`,
    for real c larger than the real part of all singularities of f.
    Note that ``t`` is always assumed real and positive.

    Return None if the integral does not exist or could not be evaluated.

    Examples
    ========

    >>> meijerint_inversion(1/x, x, t)
    Heaviside(t)

    """
    from ..core import expand, Add, Mul
    from .integrals import Integral
    from ..functions import exp, log, Heaviside
    f_ = f
    t_ = t
    t = Dummy('t', polar=True)  # We don't want sqrt(t**2) = abs(t) etc
    f = f.subs({t_: t})
    c = Dummy('c')
    _debug('Laplace-inverting', f)
    if not _is_analytic(f, x):
        _debug('But expression is not analytic.')
        return
    # We filter out exponentials here. If we are given an Add this will not
    # work, but the calling code will take care of that.
    shift = 0
    if f.is_Mul:
        args = list(f.args)
        newargs = []
        exponentials = []
        while args:
            arg = args.pop()
            if arg.is_Pow and arg.base is E:
                arg2 = expand(arg)
                if arg2.is_Mul:
                    args += arg2.args
                    continue
                try:
                    a, b = _get_coeff_exp(arg.exp, x)
                except _CoeffExpValueError:
                    b = 0
                if b == 1:
                    exponentials.append(a)
                else:
                    newargs.append(arg)
            elif arg.is_Pow:
                arg2 = expand(arg)
                if arg2.is_Mul:
                    args += arg2.args
                    continue
                if x not in arg.base.free_symbols:
                    try:
                        a, b = _get_coeff_exp(arg.exp, x)
                    except _CoeffExpValueError:
                        b = 0
                    if b == 1:
                        exponentials.append(a*log(arg.base))
                newargs.append(arg)
            else:
                newargs.append(arg)
        shift = Add(*exponentials)
        f = Mul(*newargs)

    gs = _rewrite1(f, x)
    if gs is not None:
        fac, po, g, cond = gs
        _debug('Could rewrite as single G function:', fac, po, g)
        res = Integer(0)
        for C, s, f in g:
            C, f = _rewrite_inversion(fac*C, po*x**s, f, x)
            res += C*_int_inversion(f, x, t)
            cond = And(cond, _check_antecedents_inversion(f, x))
            if cond == false:
                break
        cond = _my_unpolarify(cond)
        if cond == false:
            _debug('But cond is always False.')
        else:
            _debug('Result before branch substitution:', res)
            res = _my_unpolarify(hyperexpand(res))
            if not res.has(Heaviside):
                res *= Heaviside(t)
            res = res.subs({t: t + shift})
            if not isinstance(cond, bool):
                cond = cond.subs({t: t + shift})
            return Piecewise((res.subs({t: t_}), cond),
                             (Integral(f_*exp(x*t),
                                       (x, c - oo*I, c + oo*I)).subs({t: t_}),
                             True))
