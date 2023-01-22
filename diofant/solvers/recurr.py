"""This module is intended for solving recurrences (difference equations)."""

import collections
import functools

from ..concrete import product
from ..core import (Add, Dummy, Equality, Function, Integer, Lambda, Mul,
                    Rational, Symbol, Wild, oo)
from ..core.compatibility import iterable
from ..core.sympify import sympify
from ..functions import FallingFactorial, RisingFactorial, binomial, factorial
from ..matrices import Matrix, casoratian
from ..polys import Poly, gcd, lcm, quo, resultant, roots
from ..simplify import hypersimilar, hypersimp
from ..utilities import default_sort_key, numbered_symbols
from .ode import constantsimp
from .solvers import solve


def rsolve_poly(coeffs, f, n):
    r"""
    Find polynomial solutions for linear recurrence.

    Given linear recurrence operator `\operatorname{L}` of order
    `k` with polynomial coefficients and inhomogeneous equation
    `\operatorname{L} y = f`, where `f` is a polynomial, we seek for
    all polynomial solutions over field `K` of characteristic zero.

    Notes
    =====

    The algorithm performs two basic steps:

        1. Compute degree `N` of the general polynomial solution.
        2. Find all polynomials of degree `N` or less
           of `\operatorname{L} y = f`.

    There are two methods for computing the polynomial solutions.
    If the degree bound is relatively small, i.e. it's smaller than
    or equal to the order of the recurrence, then naive method of
    undetermined coefficients is being used. This gives system
    of algebraic equations with `N+1` unknowns.

    In the other case, the algorithm performs transformation of the
    initial equation to an equivalent one, for which the system of
    algebraic equations has only `r` indeterminates. This method is
    quite sophisticated (in comparison with the naive one) and was
    invented together by Abramov, Bronstein and Petkovšek.

    It is possible to generalize the algorithm implemented here to
    the case of linear q-difference and differential equations.

    Examples
    ========

    Lets say that we would like to compute `m`-th Bernoulli polynomial
    up to a constant, using `b(n+1) - b(n) = m n^{m-1}` recurrence:

    >>> rsolve_poly([-1, 1], 4*n**3, n)
    (C0 + n**4 - 2*n**3 + n**2, [C0])
    >>> bernoulli(4, n)
    n**4 - 2*n**3 + n**2 - 1/30

    References
    ==========

    * :cite:`Abramov1995polynomial`
    * :cite:`Petkovsek1992hyper`
    * :cite:`Petkovsek1997AeqB`

    """
    f = sympify(f)

    if not f.is_polynomial(n):
        return

    homogeneous = f.is_zero

    r = len(coeffs) - 1

    coeffs = [Poly(coeff, n) for coeff in coeffs]

    g = functools.reduce(lambda x, y: gcd(x, y, n, polys=True), coeffs + [f])
    if not g.is_ground:
        coeffs = [quo(c, g, n, polys=False) for c in coeffs]
        f = quo(f, g, n, polys=False)

    polys = [Poly(0, n)] * (r + 1)
    terms = [(Integer(0), -oo)] * (r + 1)

    for i in range(r + 1):
        for j in range(i, r + 1):
            polys[i] += coeffs[j] * binomial(j, i)

        if not polys[i].is_zero:
            (exp,), coeff = polys[i].LT()
            terms[i] = (coeff, exp)

    d = b = terms[0][1]

    for i in range(1, r + 1):
        if terms[i][1] > d:
            d = terms[i][1]

        if terms[i][1] - i > b:
            b = terms[i][1] - i

    d, b = int(d), int(b)

    x = Dummy('x')

    degree_poly = Integer(0)

    for i in range(r + 1):
        if terms[i][1] - i == b:
            degree_poly += terms[i][0] * FallingFactorial(x, i)

    nni_roots = list(roots(degree_poly, x, filter='Z',
                           predicate=lambda r: r >= 0))

    if nni_roots:
        N = [max(nni_roots)]
    else:
        N = []

    if homogeneous:
        N += [-b - 1]
    else:
        N += [f.as_poly(n).degree() - b, -b - 1]

    N = int(max(N))

    if N < 0:
        if homogeneous:
            return Integer(0), []
        return

    if N <= r:
        C = []
        y = E = Integer(0)

        for i in range(N + 1):
            C.append(Symbol('C' + str(i)))
            y += C[i] * n**i

        for i in range(r + 1):
            E += coeffs[i].as_expr() * y.subs({n: n + i})

        solutions = solve((E - f).as_poly(n).coeffs(), *C)

        if solutions:
            solutions = solutions[0]

        C = [c for c in C if c not in solutions]
        result = y.subs(solutions)
    else:
        A = r
        U = N + A + b + 1

        nni_roots = list(roots(polys[r], filter='Z',
                               predicate=lambda r: r >= 0))

        if nni_roots:
            a = max(nni_roots) + 1
        else:
            a = Integer(0)

        def _zero_vector(k):
            return [Integer(0)] * k

        def _one_vector(k):
            return [Integer(1)] * k

        def _delta(p, k):
            B = Integer(1)
            D = p.subs({n: a + k})

            for i in range(1, k + 1):
                B *= -Rational(k - i + 1, i)
                D += B * p.subs({n: a + k - i})

            return D

        alpha = {}

        for i in range(-A, d + 1):
            E = _one_vector(d + 1)

            for k in range(1, d + 1):
                E[k] = E[k - 1] * (x + i - k + 1) / k

            alpha[i] = Integer(0)

            for j in range(A + 1):
                for k in range(d + 1):
                    B = binomial(k, i + j)
                    D = _delta(polys[j].as_expr(), k)

                    alpha[i] += E[k] * B * D

        V = Matrix(U, A, lambda i, j: int(i == j))

        if homogeneous:
            for i in range(A, U):
                v = _zero_vector(A)

                for k in range(1, A + b + 1):
                    if i - k < 0:
                        break

                    B = alpha[k - A].subs({x: i - k})

                    for j in range(A):
                        v[j] += B * V[i - k, j]

                denom = alpha[-A].subs({x: i})

                for j in range(A):
                    V[i, j] = -v[j] / denom
        else:
            G = _zero_vector(U)

            for i in range(A, U):
                v = _zero_vector(A)
                g = Integer(0)

                for k in range(1, A + b + 1):
                    if i - k < 0:
                        break

                    B = alpha[k - A].subs({x: i - k})

                    for j in range(A):
                        v[j] += B * V[i - k, j]

                    g += B * G[i - k]

                denom = alpha[-A].subs({x: i})

                for j in range(A):
                    V[i, j] = -v[j] / denom

                G[i] = (_delta(f, i - A) - g) / denom

        P, Q = _one_vector(U), _zero_vector(A)

        for i in range(1, U):
            P[i] = (P[i - 1] * (n - a - i + 1) / i).expand()

        for i in range(A):
            Q[i] = Add(*[(v * p).expand() for v, p in zip(V[:, i], P)])

        if not homogeneous:
            h = Add(*[(g * p).expand() for g, p in zip(G, P)])

        C = [Symbol('C' + str(i)) for i in range(A)]

        def g2(i):
            return Add(*[c * _delta(q, i) for c, q in zip(C, Q)])

        if homogeneous:
            E = [g2(i) for i in range(N + 1, U)]
        else:
            E = [g2(i) + _delta(h, i) for i in range(N + 1, U)]

        if E != []:
            solutions = solve(E, *C)
            solutions = solutions[0]
        else:
            solutions = {}

        if homogeneous:
            result = Integer(0)
        else:
            result = h

        for c, q in list(zip(C, Q)):
            if c in solutions:
                s = solutions[c] * q
                C.remove(c)
            else:
                s = c * q

            result += s.expand()

    return result, C


def rsolve_ratio(coeffs, f, n):
    r"""
    Find rational solutions for linear recurrence.

    Given linear recurrence operator `\operatorname{L}` of order `k`
    with polynomial coefficients and inhomogeneous equation
    `\operatorname{L} y = f`, where `f` is a polynomial, we seek
    for all rational solutions over field `K` of characteristic zero.

    Notes
    =====

    The algorithm performs two basic steps:

        1. Compute polynomial `v(n)` which can be used as universal
           denominator of any rational solution of equation
           `\operatorname{L} y = f`.
        2. Construct new linear difference equation by substitution
           `y(n) = u(n)/v(n)` and solve it for `u(n)` finding all its
           polynomial solutions. Return :obj:`None` if none were found.

    Algorithm implemented here is a revised version of the original
    Abramov's algorithm, developed in 1989. The new approach is much
    simpler to implement and has better overall efficiency. This
    method can be easily adapted to q-difference equations case.

    Besides finding rational solutions alone, this functions is
    an important part of the Hyper algorithm were it is used to find
    particular solution of inhomogeneous part of a recurrence.

    Examples
    ========

    >>> rsolve_ratio([-2*n**3 + n**2 + 2*n - 1, 2*n**3 + n**2 - 6*n,
    ...               -2*n**3 - 11*n**2 - 18*n - 9,
    ...               2*n**3 + 13*n**2 + 22*n + 8], 0, n)
    (C2*(2*n - 3)/(2*(n**2 - 1)), [C2])

    References
    ==========

    * :cite:`Abramov1995rational`

    See Also
    ========

    rsolve_hyper

    """
    f = sympify(f)

    if not f.is_polynomial(n):
        return

    coeffs = list(map(sympify, coeffs))

    r = len(coeffs) - 1

    A, B = coeffs[r], coeffs[0]
    A = A.subs({n: n - r}).expand()

    h = Dummy('h')

    res = resultant(A, B.subs({n: n + h}), n)
    assert res.is_polynomial(n)

    nni_roots = list(roots(res, h, filter='Z',
                           predicate=lambda r: r >= 0))

    if not nni_roots:
        return rsolve_poly(coeffs, f, n)

    C, numers = Integer(1), [Integer(0)] * (r + 1)

    for i in range(max(nni_roots), -1, -1):
        d = gcd(A, B.subs({n: n + i}), n)

        A = quo(A, d, n)
        B = quo(B, d.subs({n: n - i}), n)

        C *= Mul(*[d.subs({n: n - j}) for j in range(i + 1)])

    denoms = [C.subs({n: n + i}) for i in range(r + 1)]

    for i in range(r + 1):
        g = gcd(coeffs[i], denoms[i], n)

        numers[i] = quo(coeffs[i], g, n)
        denoms[i] = quo(denoms[i], g, n)

    for i in range(r + 1):
        numers[i] *= Mul(*(denoms[:i] + denoms[i + 1:]))

    result = rsolve_poly(numers, f * Mul(*denoms), n)

    if result is not None:
        return (result[0] / C).simplify(), result[1]


def rsolve_hyper(coeffs, f, n):
    r"""
    Find hypergeometric solutions for linear recurrence.

    Given linear recurrence operator `\operatorname{L}` of order `k`
    with polynomial coefficients and inhomogeneous equation
    `\operatorname{L} y = f` we seek for all hypergeometric solutions
    over field `K` of characteristic zero.

    The inhomogeneous part can be either hypergeometric or a sum
    of a fixed number of pairwise dissimilar hypergeometric terms.

    Notes
    =====

    The algorithm performs three basic steps:

        1. Group together similar hypergeometric terms in the
           inhomogeneous part of `\operatorname{L} y = f`, and find
           particular solution using Abramov's algorithm.
        2. Compute generating set of `\operatorname{L}` and find basis
           in it, so that all solutions are linearly independent.
        3. Form final solution with the number of arbitrary
           constants equal to dimension of basis of `\operatorname{L}`.

    The output of this procedure is a linear combination of fixed
    number of hypergeometric terms.  However the underlying method
    can generate larger class of solutions - D'Alembertian terms.

    This method not only computes the kernel of the
    inhomogeneous equation, but also reduces in to a basis so that
    solutions generated by this procedure are linearly independent.

    Examples
    ========

    >>> rsolve_hyper([-1, 1], 1 + n, n)
    (C0 + n*(n + 1)/2, [C0])

    References
    ==========

    * :cite:`Petkovsek1992hyper`
    * :cite:`Petkovsek1997AeqB`

    """
    coeffs = list(map(sympify, coeffs))

    f = sympify(f)

    r, kernel, symbols = len(coeffs) - 1, [], set()

    if not f.is_zero:
        if f.is_Add:
            similar = {}

            for g in f.expand().args:
                if not g.is_hypergeometric(n):
                    return

                for h in list(similar):
                    if hypersimilar(g, h, n):
                        similar[h] += g
                        break
                else:
                    similar[g] = Integer(0)

            inhomogeneous = []

            for g, h in similar.items():
                inhomogeneous.append(g + h)
        elif f.is_hypergeometric(n):
            inhomogeneous = [f]
        else:
            return

        for i, g in enumerate(inhomogeneous):
            coeff, polys = Integer(1), coeffs[:]
            denoms = [Integer(1)] * (r + 1)

            g = g.simplify()
            s = hypersimp(g, n)

            for j in range(1, r + 1):
                coeff *= s.subs({n: n + j - 1})

                p, q = coeff.as_numer_denom()

                polys[j] *= p
                denoms[j] = q

            for j in range(r + 1):
                polys[j] *= Mul(*(denoms[:j] + denoms[j + 1:]))

            R = rsolve_ratio(polys, Mul(*denoms), n)
            if R is not None:
                R, syms = R
                if syms:
                    R = R.subs(zip(syms, [0] * len(syms)))

            if R:
                inhomogeneous[i] *= R
            else:
                return

            result = Add(*inhomogeneous)
            result = result.simplify()
    else:
        result = Integer(0)

    Z = Dummy('Z')

    p, q = coeffs[0], coeffs[r].subs({n: n - r + 1})

    p_factors = list(roots(p, n))
    q_factors = list(roots(q, n))

    factors = [(Integer(1), Integer(1))]

    for p in p_factors:
        for q in q_factors:
            if p.is_integer and q.is_integer and p <= q:
                continue
            factors += [(n - p, n - q)]

    p = [(n - p, Integer(1)) for p in p_factors]
    q = [(Integer(1), n - q) for q in q_factors]

    factors = p + factors + q

    for A, B in factors:
        polys, degrees = [], []
        D = A * B.subs({n: n + r - 1})

        for i in range(r + 1):
            a = Mul(*[A.subs({n: n + j}) for j in range(0, i)])
            b = Mul(*[B.subs({n: n + j}) for j in range(i, r)])

            poly = quo(coeffs[i] * a * b, D, n)
            polys.append(poly.as_poly(n))

            if not poly.is_zero:
                degrees.append(polys[i].degree())

        d, poly = max(degrees), Integer(0)

        for i in range(r + 1):
            coeff = polys[i].coeff_monomial((d,))

            if coeff != 0:
                poly += coeff * Z**i

        for z in roots(poly, Z):
            if z.is_zero:
                continue

            sol, syms = rsolve_poly([polys[i] * z**i for i in range(r + 1)],
                                    0, n)
            sol = sol.collect(syms)
            sol = [sol.coeff(_) for _ in syms]

            for C in sol:
                ratio = z * A * C.subs({n: n + 1}) / B / C
                ratio = ratio.simplify()

                skip = max([-1] + [v for v in roots(Mul(*ratio.as_numer_denom()), n)
                                   if v.is_Integer]) + 1
                K = product(ratio, (n, skip, n - 1))

                if K.has(factorial, FallingFactorial, RisingFactorial):
                    K = K.simplify()

                if casoratian(kernel + [K], n, zero=False) != 0:
                    kernel.append(K)

    kernel.sort(key=default_sort_key)
    sk = list(zip(numbered_symbols('C'), kernel))

    for C, ker in sk:
        result += C * ker

    symbols |= {s for s, k in sk}
    return result, sorted(symbols, key=default_sort_key)


def rsolve(f, *y, init={}, simplify=True):
    r"""
    Solve recurrence equations.

    The equations can involve objects of the form `y(n + k)`, where
    `k` is a constant.

    Parameters
    ==========

    f : Expr, Equality or iterable of above
        The single recurrence equation or a system of recurrence
        equations.

    \*y : tuple
        Holds function applications `y(n)`, wrt to which the recurrence
        equation(s) will be solved.  If none given (empty tuple), this
        will be guessed from the provided equation(s).

    init : dict, optional
        The initial/boundary conditions for the recurrence equations as
        mapping of the function application `y(n_i)` to its value.
        Default is empty dictionary.

    simplify : bool, optional
        Enable simplification (default) on solutions.

    Examples
    ========

    >>> eq = (n - 1)*f(n + 2) - (n**2 + 3*n - 2)*f(n + 1) + 2*n*(n + 1)*f(n)

    >>> rsolve(eq)
    [{f: Lambda(n, 2**n*C0 + C1*factorial(n))}]
    >>> rsolve(eq, init={f(0): 0, f(1): 3})
    [{f: Lambda(n, 3*2**n - 3*factorial(n))}]

    Notes
    =====

    Currently, the function can handle linear recurrences with polynomial
    coefficients and hypergeometric inhomogeneous part.

    See Also
    ========

    diofant.solvers.ode.dsolve : solving differential equations
    diofant.solvers.solvers.solve : solving algebraic equations

    """
    if not iterable(f):
        f = [f]

    f = [_.lhs - _.rhs if isinstance(_, Equality) else _ for _ in f]
    f = [_.expand() for _ in f]

    if len(f) > 1 or len(y) > 1:
        raise NotImplementedError('Support for systems of recurrence '
                                  'equations is not implemented yet.')
    f = f[0]

    if not y:
        y = sorted(f.atoms(Function), key=default_sort_key)[0]
    else:
        y = y[0]
    n = y.args[0]

    h_part = collections.defaultdict(lambda: Integer(0))
    i_part = Integer(0)

    for h, c in f.collect(y.func(Wild('n')), evaluate=False).items():
        if h.func == y.func:
            k = Wild('k', exclude=(n,))
            r = h.args[0].match(n + k)
            if r:
                c = c.simplify()
                if not c.is_rational_function(n):
                    raise ValueError(f"Rational function of '{n}' expected, got '{c}'")
                h_part[int(r[k])] = c
            else:
                raise ValueError(f"'{y.func}({n} + Integer)' expected, got '{h}'")
        else:
            i_term = h * c
            if i_term.find(y.func(Wild('k'))):
                raise NotImplementedError(f"Linear recurrence for '{y.func}' "
                                          f"expected, got '{f}'")
            i_part -= i_term

    if not i_part.is_zero:
        if not all(p.is_hypergeometric(n) for p in i_part.as_coeff_add(n)[1]):
            raise NotImplementedError('Inhomogeneous part should be a sum of '
                                      f"hypergeometric terms in '{n}', got "
                                      f"'{i_part}'")

    k_min, k_max = min(h_part), max(h_part)

    if k_min < 0:
        return rsolve(f.subs({n: n + abs(k_min)}), y, init=init,
                      simplify=simplify)

    i_numer, i_denom = i_part.as_numer_denom()

    common = functools.reduce(lcm, [x.as_numer_denom()[1]
                                    for x in h_part.values()] + [i_denom])

    if common != 1:
        for k, coeff in h_part.items():
            numer, denom = coeff.as_numer_denom()
            h_part[k] = numer * quo(common, denom, n)

        i_part = i_numer * quo(common, i_denom, n)

    coeffs = [h_part[i] for i in range(k_max + 1)]

    result = rsolve_hyper(coeffs, i_part, n)

    if result is None:
        return

    solution, symbols = result

    if symbols and init != {}:
        equations = []

        for k, v in init.items():
            if k.is_Function and k.func == y.func:
                i = int(k.args[0])
            else:
                raise ValueError(f"'{y.func}(Integer)' expected, got '{k}'")
            eq = solution.limit(n, i) - v
            equations.append(eq)

        result = solve(equations, *symbols)

        if not result:
            return
        solution = solution.subs(result[0])

    if simplify:
        solution = solution.expand(log=True, mul=False)
        solution = constantsimp(solution, symbols)
        solution = solution.simplify()

    return [{y.func: Lambda((n,), solution)}]
