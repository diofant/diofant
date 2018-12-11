r"""
This module is intended for solving recurrences or, in other words,
difference equations. Currently supported are linear, inhomogeneous
equations with polynomial or rational coefficients.

The solutions are obtained among polynomials, rational functions,
hypergeometric terms, or combinations of hypergeometric term which
are pairwise dissimilar.

``rsolve_X`` functions were meant as a low level interface
for ``rsolve`` which would use Mathematica's syntax.

Given a recurrence relation:

    .. math:: a_{k}(n) y(n+k) + a_{k-1}(n) y(n+k-1) +
              ... + a_{0}(n) y(n) = f(n)

where `k > 0` and `a_{i}(n)` are polynomials in `n`. To use
``rsolve_X`` we need to put all coefficients in to a list ``L`` of
`k+1` elements the following way:

    ``L = [ a_{0}(n), ..., a_{k-1}(n), a_{k}(n) ]``

where ``L[i]``, for `i=0, \ldots, k`, maps to
`a_{i}(n) y(n+i)` (`y(n+i)` is implicit).

For example if we would like to compute `m`-th Bernoulli polynomial
up to a constant (example was taken from rsolve_poly docstring),
then we would use `b(n+1) - b(n) = m n^{m-1}` recurrence, which
has solution `b(n) = B_m + C`.

Then ``L = [-1, 1]`` and `f(n) = m n^(m-1)` and finally for `m=4`:

>>> rsolve_poly([-1, 1], 4*n**3, n)
C0 + n**4 - 2*n**3 + n**2

>>> bernoulli(4, n)
n**4 - 2*n**3 + n**2 - 1/30

For the sake of completeness, `f(n)` can be:

    [1] a polynomial              -> rsolve_poly
    [2] a rational function       -> rsolve_ratio
    [3] a hypergeometric function  -> rsolve_hyper
"""

from collections import defaultdict

from ..concrete import product
from ..core import (Add, Dummy, Equality, Mul, Rational, S, Symbol, Wild, oo,
                    sympify)
from ..core.compatibility import default_sort_key
from ..functions import FallingFactorial, RisingFactorial, binomial, factorial
from ..matrices import Matrix, casoratian
from ..polys import Poly, gcd, gcd_list, lcm_list, quo, resultant, roots
from ..simplify import hypersimilar, hypersimp, simplify
from ..utilities import numbered_symbols
from .solvers import solve


def rsolve_poly(coeffs, f, n, **hints):
    r"""
    Given linear recurrence operator `\operatorname{L}` of order
    `k` with polynomial coefficients and inhomogeneous equation
    `\operatorname{L} y = f`, where `f` is a polynomial, we seek for
    all polynomial solutions over field `K` of characteristic zero.

    The algorithm performs two basic steps:

        (1) Compute degree `N` of the general polynomial solution.
        (2) Find all polynomials of degree `N` or less
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

    Lets say that we would like to compute `m`-th Bernoulli polynomial
    up to a constant. For this we can use `b(n+1) - b(n) = m n^{m-1}`
    recurrence, which has solution `b(n) = B_m + C`. For example:

    >>> rsolve_poly([-1, 1], 4*n**3, n)
    C0 + n**4 - 2*n**3 + n**2

    References
    ==========

    * [AbramovBronstein95]_
    * [Petkovšek92]_
    * [Petkovšek97]_
    """
    f = sympify(f)

    if not f.is_polynomial(n):
        return

    homogeneous = f.is_zero

    r = len(coeffs) - 1

    coeffs = [ Poly(coeff, n) for coeff in coeffs ]

    g = gcd_list(coeffs + [f], n, polys=True)
    if not g.is_ground:
        coeffs = [quo(c, g, n, polys=False) for c in coeffs]
        f = quo(f, g, n, polys=False)

    polys = [ Poly(0, n) ] * (r + 1)
    terms = [ (S.Zero, -oo) ] * (r + 1)

    for i in range(r + 1):
        for j in range(i, r + 1):
            polys[i] += coeffs[j]*binomial(j, i)

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

    degree_poly = S.Zero

    for i in range(r + 1):
        if terms[i][1] - i == b:
            degree_poly += terms[i][0]*FallingFactorial(x, i)

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
            if hints.get('symbols', False):
                return S.Zero, []
            else:
                return S.Zero
        else:
            return

    if N <= r:
        C = []
        y = E = S.Zero

        for i in range(N + 1):
            C.append(Symbol('C' + str(i)))
            y += C[i] * n**i

        for i in range(r + 1):
            E += coeffs[i].as_expr()*y.subs({n: n + i})

        solutions = solve((E - f).as_poly(n).coeffs(), *C)

        if solutions:
            solutions = solutions[0]

        C = [c for c in C if (c not in solutions)]
        result = y.subs(solutions)
    else:
        A = r
        U = N + A + b + 1

        nni_roots = list(roots(polys[r], filter='Z',
                               predicate=lambda r: r >= 0))

        if nni_roots != []:
            a = max(nni_roots) + 1
        else:
            a = S.Zero

        def _zero_vector(k):
            return [S.Zero] * k

        def _one_vector(k):
            return [S.One] * k

        def _delta(p, k):
            B = S.One
            D = p.subs({n: a + k})

            for i in range(1, k + 1):
                B *= -Rational(k - i + 1, i)
                D += B * p.subs({n: a + k - i})

            return D

        alpha = {}

        for i in range(-A, d + 1):
            I = _one_vector(d + 1)

            for k in range(1, d + 1):
                I[k] = I[k - 1] * (x + i - k + 1)/k

            alpha[i] = S.Zero

            for j in range(A + 1):
                for k in range(d + 1):
                    B = binomial(k, i + j)
                    D = _delta(polys[j].as_expr(), k)

                    alpha[i] += I[k]*B*D

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
                g = S.Zero

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
            P[i] = (P[i - 1] * (n - a - i + 1)/i).expand()

        for i in range(A):
            Q[i] = Add(*[ (v*p).expand() for v, p in zip(V[:, i], P) ])

        if not homogeneous:
            h = Add(*[ (g*p).expand() for g, p in zip(G, P) ])

        C = [Symbol('C' + str(i)) for i in range(A)]

        def g(i):
            return Add(*[ c*_delta(q, i) for c, q in zip(C, Q) ])

        if homogeneous:
            E = [ g(i) for i in range(N + 1, U) ]
        else:
            E = [ g(i) + _delta(h, i) for i in range(N + 1, U) ]

        if E != []:
            solutions = solve(E, *C)
            solutions = solutions[0]
        else:
            solutions = {}

        if homogeneous:
            result = S.Zero
        else:
            result = h

        for c, q in list(zip(C, Q)):
            if c in solutions:
                s = solutions[c]*q
                C.remove(c)
            else:
                s = c*q

            result += s.expand()

    if hints.get('symbols', False):
        return result, C
    else:
        return result


def rsolve_ratio(coeffs, f, n, **hints):
    r"""
    Given linear recurrence operator `\operatorname{L}` of order `k`
    with polynomial coefficients and inhomogeneous equation
    `\operatorname{L} y = f`, where `f` is a polynomial, we seek
    for all rational solutions over field `K` of characteristic zero.

    This procedure accepts only polynomials, however if you are
    interested in solving recurrence with rational coefficients
    then use ``rsolve`` which will pre-process the given equation
    and run this procedure with polynomial arguments.

    The algorithm performs two basic steps:

        (1) Compute polynomial `v(n)` which can be used as universal
            denominator of any rational solution of equation
            `\operatorname{L} y = f`.

        (2) Construct new linear difference equation by substitution
            `y(n) = u(n)/v(n)` and solve it for `u(n)` finding all its
            polynomial solutions. Return ``None`` if none were found.

    Algorithm implemented here is a revised version of the original
    Abramov's algorithm, developed in 1989. The new approach is much
    simpler to implement and has better overall efficiency. This
    method can be easily adapted to q-difference equations case.

    Besides finding rational solutions alone, this functions is
    an important part of Hyper algorithm were it is used to find
    particular solution of inhomogeneous part of a recurrence.

    Examples
    ========

    >>> rsolve_ratio([-2*x**3 + x**2 + 2*x - 1, 2*x**3 + x**2 - 6*x,
    ... - 2*x**3 - 11*x**2 - 18*x - 9, 2*x**3 + 13*x**2 + 22*x + 8], 0, x)
    C2*(2*x - 3)/(2*(x**2 - 1))

    References
    ==========

    * [Abramov95]_

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
        return rsolve_poly(coeffs, f, n, **hints)
    else:
        C, numers = S.One, [S.Zero]*(r + 1)

        for i in range(int(max(nni_roots)), -1, -1):
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

        result = rsolve_poly(numers, f * Mul(*denoms), n, **hints)

        if result is not None:
            if hints.get('symbols', False):
                return simplify(result[0] / C), result[1]
            else:
                return simplify(result / C)


def rsolve_hyper(coeffs, f, n, **hints):
    r"""
    Given linear recurrence operator `\operatorname{L}` of order `k`
    with polynomial coefficients and inhomogeneous equation
    `\operatorname{L} y = f` we seek for all hypergeometric solutions
    over field `K` of characteristic zero.

    The inhomogeneous part can be either hypergeometric or a sum
    of a fixed number of pairwise dissimilar hypergeometric terms.

    The algorithm performs three basic steps:

        (1) Group together similar hypergeometric terms in the
            inhomogeneous part of `\operatorname{L} y = f`, and find
            particular solution using Abramov's algorithm.

        (2) Compute generating set of `\operatorname{L}` and find basis
            in it, so that all solutions are linearly independent.

        (3) Form final solution with the number of arbitrary
            constants equal to dimension of basis of `\operatorname{L}`.

    Term `a(n)` is hypergeometric if it is annihilated by first order
    linear difference equations with polynomial coefficients or, in
    simpler words, if consecutive term ratio is a rational function.

    The output of this procedure is a linear combination of fixed
    number of hypergeometric terms. However the underlying method
    can generate larger class of solutions - D'Alembertian terms.

    Note also that this method not only computes the kernel of the
    inhomogeneous equation, but also reduces in to a basis so that
    solutions generated by this procedure are linearly independent

    Examples
    ========

    >>> rsolve_hyper([-1, -1, 1], 0, x)
    C0*(1/2 + sqrt(5)/2)**x + C1*(-sqrt(5)/2 + 1/2)**x

    >>> rsolve_hyper([-1, 1], 1 + x, x)
    C0 + x*(x + 1)/2

    References
    ==========

    * [Petkovšek92]_
    * [Petkovšek97]_
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
                    similar[g] = S.Zero

            inhomogeneous = []

            for g, h in similar.items():
                inhomogeneous.append(g + h)
        elif f.is_hypergeometric(n):
            inhomogeneous = [f]
        else:
            return

        for i, g in enumerate(inhomogeneous):
            coeff, polys = S.One, coeffs[:]
            denoms = [ S.One ] * (r + 1)

            g = simplify(g)
            s = hypersimp(g, n)

            for j in range(1, r + 1):
                coeff *= s.subs({n: n + j - 1})

                p, q = coeff.as_numer_denom()

                polys[j] *= p
                denoms[j] = q

            for j in range(r + 1):
                polys[j] *= Mul(*(denoms[:j] + denoms[j + 1:]))

            R = rsolve_ratio(polys, Mul(*denoms), n, symbols=True)
            if R is not None:
                R, syms = R
                if syms:
                    R = R.subs(zip(syms, [0]*len(syms)))

            if R:
                inhomogeneous[i] *= R
            else:
                return

            result = Add(*inhomogeneous)
            result = simplify(result)
    else:
        result = S.Zero

    Z = Dummy('Z')

    p, q = coeffs[0], coeffs[r].subs({n: n - r + 1})

    p_factors = list(roots(p, n))
    q_factors = list(roots(q, n))

    factors = [ (S.One, S.One) ]

    for p in p_factors:
        for q in q_factors:
            if p.is_integer and q.is_integer and p <= q:
                continue
            else:
                factors += [(n - p, n - q)]

    p = [ (n - p, S.One) for p in p_factors ]
    q = [ (S.One, n - q) for q in q_factors ]

    factors = p + factors + q

    for A, B in factors:
        polys, degrees = [], []
        D = A*B.subs({n: n + r - 1})

        for i in range(r + 1):
            a = Mul(*[A.subs({n: n + j}) for j in range(0, i)])
            b = Mul(*[B.subs({n: n + j}) for j in range(i, r)])

            poly = quo(coeffs[i]*a*b, D, n)
            polys.append(poly.as_poly(n))

            if not poly.is_zero:
                degrees.append(polys[i].degree())

        d, poly = max(degrees), S.Zero

        for i in range(r + 1):
            coeff = polys[i].coeff_monomial((d,))

            if coeff is not S.Zero:
                poly += coeff * Z**i

        for z in roots(poly, Z):
            if z.is_zero:
                continue

            sol, syms = rsolve_poly([polys[i]*z**i for i in range(r + 1)],
                                    0, n, symbols=True)
            sol = sol.collect(syms)
            sol = [sol.coeff(_) for _ in syms]

            for C in sol:
                ratio = z * A * C.subs({n: n + 1}) / B / C
                ratio = simplify(ratio)

                skip = max([-1] + [v for v in roots(Mul(*ratio.as_numer_denom()), n)
                                   if v.is_Integer]) + 1
                K = product(ratio, (n, skip, n - 1))

                if K.has(factorial, FallingFactorial, RisingFactorial):
                    K = simplify(K)

                if casoratian(kernel + [K], n, zero=False) != 0:
                    kernel.append(K)

    kernel.sort(key=default_sort_key)
    sk = list(zip(numbered_symbols('C'), kernel))

    for C, ker in sk:
        result += C * ker

    if hints.get('symbols', False):
        symbols |= {s for s, k in sk}
        return result, list(symbols)
    else:
        return result


def rsolve(f, y, init=None):
    r"""
    Solve univariate recurrence with rational coefficients.

    Given `k`-th order linear recurrence `\operatorname{L} y = f`,
    or equivalently:

    .. math:: a_{k}(n) y(n+k) + a_{k-1}(n) y(n+k-1) +
              \ldots + a_{0}(n) y(n) = f(n)

    where `a_{i}(n)`, for `i=0, \ldots, k`, are polynomials or rational
    functions in `n`, and `f` is a hypergeometric function or a sum
    of a fixed number of pairwise dissimilar hypergeometric terms in
    `n`, finds all solutions or returns ``None``, if none were found.

    Initial conditions can be given as a dictionary in two forms:

        (1) ``{   n_0  : v_0,   n_1  : v_1, ...,   n_m  : v_m }``
        (2) ``{ y(n_0) : v_0, y(n_1) : v_1, ..., y(n_m) : v_m }``

    or as a list ``L`` of values:

        ``L = [ v_0, v_1, ..., v_m ]``

    where ``L[i] = v_i``, for `i=0, \ldots, m`, maps to `y(n_i)`.

    Examples
    ========

    Lets consider the following recurrence:

    .. math:: (n - 1) y(n + 2) - (n^2 + 3 n - 2) y(n + 1) +
              2 n (n + 1) y(n) = 0

    >>> y = Function('y')

    >>> f = (n - 1)*y(n + 2) - (n**2 + 3*n - 2)*y(n + 1) + 2*n*(n + 1)*y(n)

    >>> rsolve(f, y(n))
    2**n*C0 + C1*factorial(n)

    >>> rsolve(f, y(n), {y(0): 0, y(1): 3})
    3*2**n - 3*factorial(n)

    See Also
    ========

    rsolve_poly, rsolve_ratio, rsolve_hyper
    """
    if isinstance(f, Equality):
        f = f.lhs - f.rhs
    f = f.expand()

    n = y.args[0]

    h_part = defaultdict(lambda: S.Zero)
    i_part = S.Zero

    for h, c in f.collect(y.func(Wild('n')), evaluate=False).items():
        if h.func == y.func:
            k = Wild('k', exclude=(n,))
            r = h.args[0].match(n + k)
            if r:
                c = simplify(c)
                if not c.is_rational_function(n):
                    raise ValueError("Rational function of '%s' expected, got '%s'" % (n, c))
                h_part[int(r[k])] = c
            else:
                raise ValueError("'%s(%s + Integer)' expected, got '%s'" % (y.func, n, h))
        else:
            i_term = h*c
            if i_term.find(y.func(Wild('k'))):
                raise NotImplementedError("Linear recurrence for '%s' "
                                          "expected, got '%s'" % (y.func, f))
            i_part -= i_term

    if not i_part.is_zero:
        if not all(p.is_hypergeometric(n) for p in i_part.as_coeff_add(n)[1]):
            raise NotImplementedError("Inhomogeneous part should be a sum of "
                                      "hypergeometric terms in '%s', got "
                                      "'%s'" % (n, i_part))

    k_min, k_max = min(h_part), max(h_part)

    if k_min < 0:
        return rsolve(f.subs({n: n + abs(k_min)}), y, init)

    i_numer, i_denom = i_part.as_numer_denom()

    common = lcm_list([x.as_numer_denom()[1] for x in h_part.values()] + [i_denom])

    if common is not S.One:
        for k, coeff in h_part.items():
            numer, denom = coeff.as_numer_denom()
            h_part[k] = numer*quo(common, denom, n)

        i_part = i_numer*quo(common, i_denom, n)

    coeffs = [h_part[i] for i in range(k_max + 1)]

    result = rsolve_hyper(coeffs, i_part, n, symbols=True)

    if result is None:
        return

    solution, symbols = result

    if init == {} or init == []:
        init = None

    if symbols and init is not None:
        if type(init) is list:
            init = {i: init[i] for i in range(len(init))}

        equations = []

        for k, v in init.items():
            try:
                i = int(k)
            except TypeError:
                if k.is_Function and k.func == y.func:
                    i = int(k.args[0])
                else:
                    raise ValueError("Integer or term '%s(Integer)' "
                                     "expected, got '%s'" % (y.func, k))
            eq = solution.limit(n, i) - v
            equations.append(eq)

        result = solve(equations, *symbols)

        if not result:
            return
        else:
            solution = solution.subs(result[0])

    return solution
