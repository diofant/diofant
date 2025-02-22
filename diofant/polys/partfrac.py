"""Algorithms for partial fraction decomposition of rational functions."""

import itertools

from ..core import Add, Dummy, Function, Integer, Lambda, preorder_traversal
from ..core.sympify import sympify
from ..utilities import numbered_symbols
from . import Poly, RootSum, cancel, factor
from .polyerrors import PolynomialError
from .polyoptions import allowed_flags, set_defaults
from .polytools import parallel_poly_from_expr


__all__ = 'apart', 'apart_list', 'assemble_partfrac_list'


def apart(f, x=None, full=False, **options):
    """
    Compute partial fraction decomposition of a rational function.

    Given a rational function ``f``, computes the partial fraction
    decomposition of ``f``. Two algorithms are available: One is based on the
    undertermined coefficients method, the other is Bronstein's full partial
    fraction decomposition algorithm.

    The undetermined coefficients method (selected by ``full=False``) uses
    polynomial factorization (and therefore accepts the same options as
    factor) for the denominator. Per default it works over the rational
    numbers, therefore decomposition of denominators with non-rational roots
    (e.g. irrational, complex roots) is not supported by default (see options
    of factor).

    Bronstein's algorithm can be selected by using ``full=True`` and allows a
    decomposition of denominators with non-rational roots. A human-readable
    result can be obtained via ``doit()`` (see examples below).

    Examples
    ========

    By default, using the undetermined coefficients method:

    >>> apart(y/(x + 2)/(x + 1), x)
    -y/(x + 2) + y/(x + 1)

    The undetermined coefficients method does not provide a result when the
    denominators roots are not rational:

    >>> apart(y/(x**2 + x + 1), x)
    y/(x**2 + x + 1)

    You can choose Bronstein's algorithm by setting ``full=True``:

    >>> apart(y/(x**2 + x + 1), x, full=True)
    RootSum(_w**2 + _w + 1, Lambda(_a, (-2*y*_a/3 - y/3)/(x - _a)))

    Calling ``doit()`` yields a human-readable result:

    >>> apart(y/(x**2 + x + 1), x, full=True).doit()
    (-y/3 - 2*y*(-1/2 - sqrt(3)*I/2)/3)/(x + 1/2 + sqrt(3)*I/2) + (-y/3 -
        2*y*(-1/2 + sqrt(3)*I/2)/3)/(x + 1/2 - sqrt(3)*I/2)

    See Also
    ========

    apart_list, assemble_partfrac_list

    """
    allowed_flags(options, [])

    f = sympify(f)

    if f.is_Atom:
        return f
    P, Q = f.as_numer_denom()

    _options = options.copy()
    options = set_defaults(options, extension=True)
    try:
        (P, Q), _ = parallel_poly_from_expr((P, Q), x, **options)
    except PolynomialError as exc:
        if f.is_commutative:
            raise PolynomialError(exc) from exc
        # non-commutative
        if f.is_Mul:
            c, nc = f.args_cnc(split_1=False)
            nc = f.func(*[apart(i, x=x, full=full, **_options) for i in nc])
            if c:
                c = apart(f.func._from_args(c), x=x, full=full, **_options)
                return c*nc
            return nc
        if f.is_Add:
            c = []
            nc = []
            for i in f.args:
                if i.is_commutative:
                    c.append(i)
                else:
                    nc.append(apart(i, x=x, full=full, **_options))
            return apart(f.func(*c), x=x, full=full, **_options) + f.func(*nc)
        reps = []
        pot = preorder_traversal(f)
        next(pot)
        for e in pot:
            reps.append((e, apart(e, x=x, full=full, **_options)))
            pot.skip()  # this was handled successfully
        return f.xreplace(dict(reps))

    if P.is_multivariate:
        fc = f.cancel()
        if fc != f:
            return apart(fc, x=x, full=full, **_options)

        raise NotImplementedError(
            'multivariate partial fraction decomposition')

    common, P, Q = P.cancel(Q)

    poly, P = P.div(Q, auto=True)
    P, Q = P.rat_clear_denoms(Q)

    if Q.degree() <= 1:
        partial = P/Q
    else:
        if not full:
            partial = apart_undetermined_coeffs(P, Q)
        else:
            partial = apart_full_decomposition(P, Q)

    terms = Integer(0)

    for term in Add.make_args(partial):
        if term.has(RootSum):
            terms += term
        else:
            terms += factor(term)

    return common*(poly.as_expr() + terms)


def apart_undetermined_coeffs(P, Q):
    """Partial fractions via method of undetermined coefficients."""
    X = numbered_symbols(cls=Dummy)
    partial, symbols = [], []

    _, factors = Q.factor_list()

    for f, k in factors:
        n, q = f.degree(), Q

        for i in range(1, k + 1):
            coeffs, q = list(itertools.islice(X, n)), q.quo(f)
            partial.append((coeffs, q, f, i))
            symbols.extend(coeffs)

    dom = Q.domain.inject(*symbols)
    F = Integer(0).as_poly(Q.gen, domain=dom)

    for i, (coeffs, q, f, k) in enumerate(partial):
        h = Poly(coeffs, Q.gen, domain=dom)
        partial[i] = (h, f, k)
        q = q.set_domain(dom)
        F += h*q

    system, result = [], Integer(0)

    for (k,), coeff in F.terms():
        system.append(coeff - P.coeff_monomial((k,)))

    from ..solvers import solve
    solution = solve(system, symbols)[0]

    for h, f, k in partial:
        h = h.as_expr().subs(solution)
        result += h/f.as_expr()**k

    return result


def apart_full_decomposition(P, Q):
    """
    Bronstein's full partial fraction decomposition algorithm.

    Given a univariate rational function ``f``, performing only GCD
    operations over the algebraic closure of the initial ground domain
    of definition, compute full partial fraction decomposition with
    fractions having linear denominators.

    Note that no factorization of the initial denominator of ``f`` is
    performed. The final decomposition is formed in terms of a sum of
    :class:`RootSum` instances.

    References
    ==========

    * :cite:`Bronstein1993partial`

    """
    return assemble_partfrac_list(apart_list(P/Q, P.gens[0]))


def apart_list(f, x=None, dummies=None, **options):
    """
    Compute partial fraction decomposition of a rational function
    and return the result in structured form.

    Given a rational function ``f`` compute the partial fraction decomposition
    of ``f``. Only Bronstein's full partial fraction decomposition algorithm
    is supported by this method. The return value is highly structured and
    perfectly suited for further algorithmic treatment rather than being
    human-readable. The function returns a tuple holding three elements:

    * The first item is the common coefficient, free of the variable `x` used
      for decomposition. (It is an element of the base field `K`.)

    * The second item is the polynomial part of the decomposition. This can be
      the zero polynomial. (It is an element of `K[x]`.)

    * The third part itself is a list of quadruples. Each quadruple
      has the following elements in this order:

      - The (not necessarily irreducible) polynomial `D` whose roots `w_i` appear
        in the linear denominator of a bunch of related fraction terms. (This item
        can also be a list of explicit roots. However, at the moment ``apart_list``
        never returns a result this way, but the related ``assemble_partfrac_list``
        function accepts this format as input.)

      - The numerator of the fraction, written as a function of the root `w`

      - The linear denominator of the fraction *excluding its power exponent*,
        written as a function of the root `w`.

      - The power to which the denominator has to be raised.

    On can always rebuild a plain expression by using the function ``assemble_partfrac_list``.

    Examples
    ========

    A first example:

    >>> f = (2*x**3 - 2*x) / (x**2 - 2*x + 1)
    >>> pfd = apart_list(f)
    >>> pfd
    (1,
    Poly(2*x + 4, x, domain='ZZ'),
    [(Poly(_w - 1, _w, domain='ZZ'), Lambda(_a, 4), Lambda(_a, x - _a), 1)])

    >>> assemble_partfrac_list(pfd)
    2*x + 4 + 4/(x - 1)

    Second example:

    >>> f = (-2*x - 2*x**2) / (3*x**2 - 6*x)
    >>> pfd = apart_list(f)
    >>> pfd
    (-1,
    Poly(2/3, x, domain='QQ'),
    [(Poly(_w - 2, _w, domain='ZZ'), Lambda(_a, 2), Lambda(_a, x - _a), 1)])

    >>> assemble_partfrac_list(pfd)
    -2/3 - 2/(x - 2)

    Another example, showing symbolic parameters:

    >>> pfd = apart_list(t/(x**2 + x + t), x)
    >>> pfd
    (1,
    Poly(0, x, domain='ZZ[t]'),
    [(Poly(_w**2 + _w + t, _w, domain='ZZ[t]'),
    Lambda(_a, -2*t*_a/(4*t - 1) - t/(4*t - 1)),
    Lambda(_a, x - _a), 1)])

    >>> assemble_partfrac_list(pfd)
    RootSum(t + _w**2 + _w, Lambda(_a, (-2*t*_a/(4*t - 1) - t/(4*t - 1))/(x - _a)))

    This example is taken from Bronstein's original paper:

    >>> f = 36 / (x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2)
    >>> pfd = apart_list(f)
    >>> pfd
    (1,
    Poly(0, x, domain='ZZ'),
    [(Poly(_w - 2, _w, domain='ZZ'), Lambda(_a, 4), Lambda(_a, x - _a), 1),
    (Poly(_w**2 - 1, _w, domain='ZZ'), Lambda(_a, -3*_a - 6), Lambda(_a, x - _a), 2),
    (Poly(_w + 1, _w, domain='ZZ'), Lambda(_a, -4), Lambda(_a, x - _a), 1)])

    >>> assemble_partfrac_list(pfd)
    -4/(x + 1) - 3/(x + 1)**2 - 9/(x - 1)**2 + 4/(x - 2)

    See also
    ========

    apart, assemble_partfrac_list

    References
    ==========

    * :cite:`Bronstein1993partial`

    """
    allowed_flags(options, [])

    f = sympify(f)

    if f.is_Atom:
        return f
    P, Q = f.as_numer_denom()

    options = set_defaults(options, extension=True)
    (P, Q), _ = parallel_poly_from_expr((P, Q), x, **options)

    if P.is_multivariate:
        raise NotImplementedError('multivariate partial fraction decomposition')

    common, P, Q = P.cancel(Q)

    poly, P = P.div(Q, auto=True)
    P, Q = P.rat_clear_denoms(Q)

    polypart = poly

    if dummies is None:
        def dummies(name):
            d = Dummy(name)
            while True:
                yield d

        dummies = dummies('w')

    rationalpart = apart_list_full_decomposition(P, Q, dummies)

    return common, polypart, rationalpart


def apart_list_full_decomposition(P, Q, dummygen):
    """
    Bronstein's full partial fraction decomposition algorithm.

    Given a univariate rational function ``f``, performing only GCD
    operations over the algebraic closure of the initial ground domain
    of definition, compute full partial fraction decomposition with
    fractions having linear denominators.

    Note that no factorization of the initial denominator of ``f`` is
    performed. The final decomposition is formed in terms of a sum of
    :class:`RootSum` instances.

    References
    ==========

    * :cite:`Bronstein1993partial`

    """
    f, x, U = P/Q, P.gen, []

    u = Function('u')(x)
    a = Dummy('a')

    Q_c, Q_sqf = Q.sqf_list()
    if Q_c != 1 and Q_sqf:
        if Q_sqf[0][1] == 1:
            Q_sqf[0] = Q_c*Q_sqf[0][0], 1
        else:
            Q_sqf.insert(0, (Q_c.as_poly(x), 1))

    partial = []

    for d, n in Q_sqf:
        b = d.as_expr()
        U += [u.diff((x, n - 1))]

        h = cancel(f*b**n) / u**n

        H, subs = [h], []

        for j in range(1, n):
            H += [H[-1].diff(x) / j]

        for j in range(1, n + 1):
            subs += [(U[j - 1], b.diff((x, j)) / j)]

        for j in range(n):
            P, Q = cancel(H[j]).as_numer_denom()

            for i in range(j + 1):
                P = P.subs([subs[j - i]])

            Q = Q.subs([subs[0]])

            P = P.as_poly(x)
            Q = Q.as_poly(x)

            G = P.gcd(d)
            D = d.quo(G)

            g, B = Q.half_gcdex(D)
            b = (P * B.quo(g)).rem(D)

            Dw = D.subs({x: next(dummygen)})
            numer = Lambda(a, b.as_expr().subs({x: a}))
            denom = Lambda(a, (x - a))
            exponent = n-j

            partial.append((Dw, numer, denom, exponent))

    return partial


def assemble_partfrac_list(partial_list):
    r"""Reassemble a full partial fraction decomposition
    from a structured result obtained by the function ``apart_list``.

    Examples
    ========

    This example is taken from Bronstein's original paper:

    >>> f = 36 / (x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2)
    >>> pfd = apart_list(f)
    >>> pfd
    (1,
    Poly(0, x, domain='ZZ'),
    [(Poly(_w - 2, _w, domain='ZZ'), Lambda(_a, 4), Lambda(_a, x - _a), 1),
    (Poly(_w**2 - 1, _w, domain='ZZ'), Lambda(_a, -3*_a - 6), Lambda(_a, x - _a), 2),
    (Poly(_w + 1, _w, domain='ZZ'), Lambda(_a, -4), Lambda(_a, x - _a), 1)])

    >>> assemble_partfrac_list(pfd)
    -4/(x + 1) - 3/(x + 1)**2 - 9/(x - 1)**2 + 4/(x - 2)

    If we happen to know some roots we can provide them easily inside the structure:

    >>> pfd = apart_list(2/(x**2-2))
    >>> pfd
    (1,
    Poly(0, x, domain='ZZ'),
    [(Poly(_w**2 - 2, _w, domain='ZZ'),
    Lambda(_a, _a/2), Lambda(_a, x - _a),
    1)])

    >>> pfda = assemble_partfrac_list(pfd)
    >>> pfda
    RootSum(_w**2 - 2, Lambda(_a, _a/(x - _a)))/2

    >>> pfda.doit()
    -sqrt(2)/(2*(x + sqrt(2))) + sqrt(2)/(2*(x - sqrt(2)))

    >>> a = Dummy('a')
    >>> pfd = (1, Integer(0).as_poly(x),
    ...        [([sqrt(2), -sqrt(2)],
    ...          Lambda(a, a/2), Lambda(a, -a + x), 1)])

    >>> assemble_partfrac_list(pfd)
    -sqrt(2)/(2*(x + sqrt(2))) + sqrt(2)/(2*(x - sqrt(2)))

    See also
    ========

    apart, apart_list

    """
    # Common factor
    common = partial_list[0]

    # Polynomial part
    polypart = partial_list[1]
    pfd = polypart.as_expr()

    # Rational parts
    for r, nf, df, ex in partial_list[2]:
        if isinstance(r, Poly):
            # Assemble in case the roots are given implicitly by a polynomials
            an, nu = nf.variables, nf.expr
            ad, de = df.variables, df.expr
            # Hack to make dummies equal because Lambda created new Dummies
            de = de.subs({ad[0]: an[0]})
            func = Lambda(an, nu/de**ex)
            pfd += RootSum(r, func, auto=False, quadratic=False)
        else:
            # Assemble in case the roots are given explicitly by a list of algebraic numbers
            for root in r:
                pfd += nf(root)/df(root)**ex

    return common*pfd
