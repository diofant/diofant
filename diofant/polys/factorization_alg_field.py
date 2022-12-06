import functools
import operator
import random

from ..config import using
from ..core import Dummy
from ..domains.algebraicfield import AlgebraicElement
from ..integrals.heurisch import _symbols
from ..ntheory import nextprime
from .modulargcd import (_euclidean_algorithm, _gf_gcdex, _minpoly_from_dense,
                         _trunc)
from .polyerrors import NotInvertible, UnluckyLeadingCoefficient
from .polyutils import _sort_factors
from .rings import PolynomialRing
from .solvers import solve_lin_sys


# TODO
# ====

# -) efficiency of _factor can be improved for irreducible polynomials if the
#    univariate factorization is done before the LC is factored


def _alpha_to_z(f, ring):
    r"""
    Change the representation of a polynomial over
    `\mathbb Q(\alpha)` by replacing the algebraic element `\alpha` by a new
    variable `z`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `\mathbb Q(\alpha)[x_0, \ldots, x_n]`
    ring : PolynomialRing
        the polynomial ring `\mathbb Q[x_0, \ldots, x_n, z]`

    Returns
    =======

    f_ : PolyElement
        polynomial in `\mathbb Q[x_0, \ldots, x_n, z]`

    """
    if isinstance(f, AlgebraicElement):
        ring = ring.drop(*ring.gens[:-1])
        f_ = ring(dict(f.rep))

    else:
        f_ = ring.zero

        for monom, coeff in f.items():
            coeff = coeff.rep.all_coeffs()
            n = len(coeff)

            for i in range(n):
                m = monom + (n-i-1,)
                if coeff[n - i - 1]:
                    f_[m] = coeff[n - i - 1]

    return f_


def _z_to_alpha(f, ring):
    r"""
    Change the representation of a polynomial in
    `\mathbb Q[x_0, \ldots, x_n, z]` by replacing the variable `z` by the
    algebraic element `\alpha` of the given ring
    `\mathbb Q(\alpha)[x_0, \ldots, x_n]`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `\mathbb Q[x_0, \ldots, x_n, z]`
    ring : PolynomialRing
        the polynomial ring `\mathbb Q(\alpha)[x_0, \ldots, x_n]`

    Returns
    =======

    f_ : PolyElement
        polynomial in `\mathbb Q(\alpha)[x_0, \ldots, x_n]`

    """
    domain = ring.domain

    f_ = ring.zero
    for monom, coeff in f.items():
        m = monom[:-1]
        c = domain([0]*monom[-1] + [domain.domain(coeff)])

        if m not in f_:
            f_[m] = c
        else:
            f_[m] += c

    return f_


def _distinct_prime_divisors(S, domain):
    r"""
    Try to find pairwise coprime divisors of all elements of a given list
    `S` of integers.

    If this fails, ``None`` is returned.

    References
    ==========

    * :cite:`Javadi2009factor`

    """
    gcd = domain.gcd
    divisors = []

    for i, s in enumerate(S):
        divisors.append(s)

        for j in range(i):
            g = gcd(divisors[i], divisors[j])
            divisors[i] = divisors[i] // g
            divisors[j] = divisors[j] // g
            g1 = gcd(divisors[i], g)
            g2 = gcd(divisors[j], g)

            while g1 != 1:
                g1 = gcd(divisors[i], g1)
                divisors[i] = divisors[i] // g1

            while g2 != 1:
                g2 = gcd(divisors[j], g2)
                divisors[j] = divisors[j] // g2

            if divisors[i] == 1 or divisors[j] == 1:
                return

    return divisors


def _denominator(f):
    r"""
    Compute the denominator `\mathrm{den}(f)` of a polynomial `f` over
    `\mathbb Q`, i.e. the smallest integer such that `\mathrm{den}(f) f` is
    a polynomial over `\mathbb Z`.

    """
    ring = f.ring.domain.ring
    lcm = ring.lcm
    den = ring.one

    for coeff in f.values():
        den = lcm(den, coeff.denominator)

    return den


def _monic_associate(f, ring):
    r"""
    Compute the monic associate of a polynomial `f` over
    `\mathbb Q(\alpha)`, which is defined as

    .. math ::

        \mathrm{den}\left( \frac 1 {\mathrm{lc}(f)} f \right) \cdot \frac 1 {\mathrm{lc}(f)} f.

    The result is a polynomial in `\mathbb Z[x_0, \ldots, x_n, z]`.

    See also
    ========

    _denominator
    _alpha_to_z

    """
    qring = ring.clone(domain=ring.domain.field)
    f = _alpha_to_z(f.monic(), qring)
    f_ = f.clear_denoms()[1].set_ring(ring)

    return f_.primitive()[1]


def _leading_coeffs(f, U, gamma, lcfactors, A, D, denoms, divisors):
    r"""
    Compute the true leading coefficients in `x_0` of the irreducible
    factors of a polynomial `f`.

    If this fails, ``None`` is returned.

    Parameters
    ==========

    f : PolyElement
        squarefree polynomial in `Z[x_0, \ldots, x_n, z]`
    U : list of PolyElement objects
        monic univariate factors of `f(x_0, A)` in `\mathbb Q(\alpha)[x_0]`
    gamma : Integer
        integer content of `\mathrm{lc}_{x_0}(f)`
    lcfactors : list of (PolyElement, Integer) objects
        factorization of `\mathrm{lc}_{x_0}(f)` in
        `\mathbb Z[x_1, \ldots, x_n, z]`
    A : list of Integer objects
        the evaluation point `[a_1, \ldots, a_n]`
    D : Integer
        integral multiple of the defect of `\mathbb Q(\alpha)`
    denoms : list of Integer objects
        denominators of `\frac 1 {l(A)}` for `l` in ``lcfactors``
    divisors : list of Integer objects
        pairwise coprime divisors of all elements of ``denoms``

    Returns
    =======

    f : PolyElement
        possibly updated polynomial `f`
    lcs : list of PolyElement objects
        true leading coefficients of the irreducible factors of `f`
    U_ : list of PolyElement objects
        list of possibly updated monic associates of the univariate factors
        `U`

    References
    ==========

    * :cite:`Javadi2009factor`

    """
    ring = f.ring
    domain = ring.domain
    symbols = f.ring.symbols
    qring = ring.clone(symbols=(symbols[0], symbols[-1]), domain=domain.field)
    gcd = domain.gcd

    U = [_alpha_to_z(u, qring) for u, _ in U]
    denominators = [_denominator(u) for u in U]

    omega = D * gamma

    m = len(denoms)

    for i in range(m):
        pi = gcd(omega, divisors[i])
        divisors[i] //= pi
        if divisors[i] == 1:
            raise NotImplementedError

    e = []

    for dj in denominators:
        ej = []

        for i in range(m):
            eji = 0
            g1 = gcd(dj, divisors[i])

            while g1 != 1:
                eji += 1
                dj = dj // g1
                g1 = gcd(dj, g1)

            ej.append(eji)

        e.append(ej)

    n = len(denominators)
    if any(sum(e[j][i] for j in range(n)) != lcfactors[i][1] for i in range(m)):
        return

    lcring = ring.drop(0)
    lcs = []
    for j in range(n):
        lj = functools.reduce(operator.mul, [lcfactors[i][0]**e[j][i] for i in range(m)], lcring.one)
        lcs.append(lj)

    zring = qring.clone(domain=domain)

    for j in range(n):
        lj = lcs[j]
        dj = denominators[j]
        ljA = lj.eval(list(zip(lcring.gens, A)))

        lcs[j] = lj*dj
        U[j] = (U[j]*dj).set_ring(zring) * ljA.set_ring(zring)

        if omega == 1:
            f *= dj
        else:
            d = gcd(omega, dj)
            f *= (dj // d)

    if omega != 1:
        lcs[0] *= omega
        U[0] *= omega

    return f, lcs, U


def _test_evaluation_points(f, gamma, lcfactors, A, D):
    r"""
    Test if an evaluation point is suitable for _factor.

    If it is not, ``None`` is returned.

    Parameters
    ==========

    f : PolyElement
        squarefree polynomial in `\mathbb Z[x_0, \ldots, x_n, z]`
    gamma : Integer
        leading coefficient of `f` in `\mathbb Z`
    lcfactors : list of (PolyElement, Integer) objects
        factorization of `\mathrm{lc}_{x_0}(f)` in
        `\mathbb Z[x_1, \ldots, x_n, z]`
    A : list of Integer objects
        the evaluation point `[a_1, \ldots, a_n]`
    D : Integer
        integral multiple of the defect of `\mathbb Q(\alpha)`

    Returns
    =======

    fA : PolyElement
        `f` evaluated at `A`, i.e. `f(x_0, A)`
    denoms : list of Integer objects
        the denominators of `\frac 1 {l(A)}` for `l` in ``lcfactors``
    divisors : list of Integer objects
        pairwise coprime divisors of all elements of ``denoms``

    References
    ==========

    * :cite:`Javadi2009factor`

    See also
    ========

    _factor

    """
    ring = f.ring
    qring = ring.clone(domain=ring.domain.field)

    fA = f.eval(list(zip(ring.gens[1:-1], A)))

    if fA.degree() < f.degree():
        return

    if not fA.is_squarefree:
        return

    omega = gamma * D
    denoms = []
    for l, _ in lcfactors:
        lA = l.eval(list(zip(l.ring.gens, A)))  # in Q(alpha)
        denoms.append(_denominator(_alpha_to_z(lA**(-1), qring)))

    if any(denoms.count(denom) > 1 for denom in denoms):
        raise UnluckyLeadingCoefficient

    divisors = _distinct_prime_divisors(denoms, ring.domain)

    if divisors is None:
        return
    elif any(omega % d == 0 for d in divisors):
        return

    return fA, denoms, divisors


def _subs_ground(f, A):
    r"""
    Substitute variables in the coefficients of a polynomial `f` over a
    ``PolynomialRing``.

    """
    f_ = f.ring.zero

    for monom, coeff in f.items():
        if coeff.compose(A):
            f_[monom] = coeff.compose(A)

    return f_


def _padic_lift(f, pfactors, lcs, B, minpoly, p):
    r"""
    Lift the factorization of a polynomial over `\mathbb Z_p[z]/(\mu(z))` to
    a factorization over `\mathbb Z_{p^m}[z]/(\mu(z))`, where `p^m \geq B`.

    If this fails, ``None`` is returned.

    Parameters
    ==========

    f : PolyElement
        squarefree polynomial in `\mathbb Z[x_0, \ldots, x_n, z]`
    pfactors : list of PolyElement objects
        irreducible factors of `f` modulo `p`
    lcs : list of PolyElement objects
        true leading coefficients in `x_0` of the irreducible factors of `f`
    B : Integer
        heuristic numerical bound on the size of the largest integer
        coefficient in the irreducible factors of `f`
    minpoly : PolyElement
        minimal polynomial `\mu` of `\alpha` over `\mathbb Q`
    p : Integer
        prime number

    Returns
    =======

    H : list of PolyElement objects
        factorization of `f` modulo `p^m`, where `p^m \geq B`

    References
    ==========

    * :cite:`Javadi2009factor`

    """
    ring = f.ring
    domain = ring.domain

    x = ring.gens[0]
    tails = [g - g.eject(*ring.gens[1:]).LC.set_ring(ring)*x**g.degree() for g in pfactors]

    coeffs = []
    for i, g in enumerate(tails):
        coeffs += _symbols(f'c{i}', len(g))

    coeffring = PolynomialRing(domain, coeffs)
    ring_ = ring.clone(domain=coeffring)

    S = []
    k = 0
    for t in tails:
        s = ring_.zero
        r = len(t)
        for i, monom in zip(range(k, k + r), t):
            s[monom] = coeffring.gens[i]
        S.append(s)
        k += r

    m = minpoly.set_ring(ring_)
    f = f.set_ring(ring_)
    x = ring_.gens[0]
    H = [t.set_ring(ring_) + li.set_ring(ring_)*x**g.degree() for t, g, li in
         zip(tails, pfactors, lcs)]

    prod = functools.reduce(operator.mul, H)
    e = (f - prod) % m

    P = domain(p)
    while e and P < 2*B:
        poly = e // P

        for s, h in zip(S, H):
            poly -= (prod//h)*s

        poly = _trunc(poly, m, P)

        P_domain = domain.finite_ring(P)

        try:
            solution = solve_lin_sys([_.set_domain(P_domain)
                                      for _ in poly.values()],
                                     coeffring.clone(domain=P_domain))
        except NotInvertible:
            return

        if solution is None:
            return
        else:
            solution = {k.set_domain(domain): v.set_domain(domain).trunc_ground(P)
                        for k, v in solution.items()}
            assert len(solution) == coeffring.ngens

        subs = list(solution.items())

        H = [h + _subs_ground(s, subs)*P for h, s in zip(H, S)]
        P = P**2
        prod = functools.reduce(operator.mul, H)
        e = (f - prod) % m

    if e == 0:
        return [h.set_ring(ring) for h in H]
    else:
        return


def _div(f, g, minpoly, p):
    r"""
    Division with remainder for univariate polynomials over
    `\mathbb Z_p[z]/(\mu(z))`.

    """
    ring = f.ring
    domain = ring.domain

    rem = f
    deg = g.degree(0)
    lcinv, _, gcd = _gf_gcdex(g.eject(*ring.gens[1:]).LC, minpoly, p)

    if not gcd == 1:
        raise NotImplementedError

    quotient = ring.zero

    while True:
        degrem = rem.degree(0)
        if degrem < deg:
            break
        m = ring.from_terms([((degrem - deg, 0), domain.one)])
        quo = (lcinv * rem.eject(*ring.gens[1:]).LC).set_ring(ring)*m
        rem = _trunc(rem - g*quo, minpoly, p)
        quotient += quo

    return _trunc(quotient, minpoly, p), rem


def _extended_euclidean_algorithm(f, g, minpoly, p):
    r"""
    Extended Euclidean Algorithm for univariate polynomials over
    `\mathbb Z_p[z]/(\mu(z))`.

    Returns `s, t, h`, where `h` is the GCD of `f` and `g` and
    `sf + tg = h`.

    """
    ring = f.ring
    zero = ring.zero
    one = ring.one

    f = _trunc(f, minpoly, p)
    g = _trunc(g, minpoly, p)

    s0, s1 = zero, one
    t0, t1 = one, zero

    while g:
        result = _div(f, g, minpoly, p)
        if result is None:
            raise NotImplementedError
        quo, rem = result
        f, g = g, rem
        s0, s1 = s1 - quo*s0, s0
        t0, t1 = t1 - quo*t0, t0

    lcfinv = _gf_gcdex(f.eject(*ring.gens[1:]).LC, minpoly, p)[0].set_ring(ring)

    return (_trunc(s1 * lcfinv, minpoly, p), _trunc(t1 * lcfinv, minpoly, p),
            _trunc(f * lcfinv, minpoly, p))


def _diophantine_univariate(F, m, minpoly, p):
    r"""
    Solve univariate Diophantine equations of the form

    .. math ::

        \sum_{f \in F} \left( h_f(x) \cdot \prod_{g \in F \setminus \lbrace f \rbrace } g(x) \right) = x^m

    over `\mathbb Z_p[z]/(\mu(z))`.

    """
    ring = F[0].ring
    domain = ring.domain
    m = ring.from_terms([((m, 0), domain.one)])

    if len(F) == 2:
        f, g = F
        result = _extended_euclidean_algorithm(g, f, minpoly, p)
        if result is None:
            raise NotImplementedError
        s, t, _ = result

        s *= m
        t *= m

        q, s = _div(s, f, minpoly, p)

        t += q*g

        s = _trunc(s, minpoly, p)
        t = _trunc(t, minpoly, p)

        result = [s, t]
    else:
        G = [F[-1]]

        for f in reversed(F[1:-1]):
            G.insert(0, f * G[0])

        S, T = [], [ring.one]

        for f, g in zip(F, G):
            result = _diophantine([g, f], T[-1], [], 0, minpoly, p)
            if result is None:
                raise NotImplementedError
            t, s = result
            T.append(t)
            S.append(s)

        result, S = [], S + [T[-1]]

        for s, f in zip(S, F):
            r = _div(s*m, f, minpoly, p)[1]
            s = _trunc(r, minpoly, p)

            result.append(s)

    return result


def _diophantine(F, c, A, d, minpoly, p):
    r"""
    Solve multivariate Diophantine equations over `\mathbb Z_p[z]/(\mu(z))`.

    """
    ring = c.ring

    if not A:
        S = [ring.zero for _ in F]
        c = _trunc(c, minpoly, p)

        for (exp,), coeff in c.eject(1).items():
            T = _diophantine_univariate(F, exp, minpoly, p)
            if T is None:
                raise NotImplementedError

            for j, (s, t) in enumerate(zip(S, T)):
                S[j] = _trunc(s + t*coeff.set_ring(ring), minpoly, p)
    else:
        n = len(A)
        e = functools.reduce(operator.mul, F)

        a, A = A[-1], A[:-1]
        B, G = [], []

        for f in F:
            B.append(e//f)
            G.append(f.eval(n, a))

        C = c.eval(n, a)

        S = _diophantine(G, C, A, d, minpoly, p)
        if S is None:
            raise NotImplementedError
        S = [s.set_ring(ring) for s in S]

        for s, b in zip(S, B):
            c = c - s*b

        c = _trunc(c, minpoly, p)

        m = ring.gens[n] - a
        M = ring.one

        for k in range(d):
            if not c:
                break

            M = M * m
            C = c.diff(x=n, m=k + 1).eval(x=n, a=a)

            if C:
                C = C.quo_ground(ring.domain.factorial(k + 1))
                T = _diophantine(G, C, A, d, minpoly, p)
                if T is None:
                    raise NotImplementedError

                for i, t in enumerate(T):
                    T[i] = t.set_ring(ring) * M

                for i, (s, t) in enumerate(zip(S, T)):
                    S[i] = s + t

                for t, b in zip(T, B):
                    c = c - t * b

                c = _trunc(c, minpoly, p)
            else:
                raise NotImplementedError

        S = [_trunc(s, minpoly, p) for s in S]

    return S


def _hensel_lift(f, H, LC, A, minpoly, p):
    r"""
    Parallel Hensel lifting algorithm over `\mathbb  Z_p[z]/(\mu(z))`.

    Parameters
    ==========

    f : PolyElement
        squarefree polynomial in `\mathbb Z[x_0, \ldots, x_n, z]`
    H : list of PolyElement objects
        monic univariate factors of `f(x_0, A)` in
        `\mathbb Z[x_0, z]`
    LC : list of PolyElement objects
        true leading coefficients of the irreducible factors of `f`
    A : list of Integer objects
        the evaluation point `[a_1, \ldots, a_n]`
    p : Integer
        prime number

    Returns
    =======

    pfactors : list of PolyElement objects
        irreducible factors of `f` modulo `p`

    """
    ring = f.ring
    n = len(A)

    S = [f]
    H = list(H)

    for i, a in enumerate(reversed(A[1:])):
        s = S[0].eval(n - i, a)
        S.insert(0, _trunc(s, minpoly, p))

    d = max(f.degree(_) for _ in ring.gens[1:-1])

    for j, s, a in zip(range(1, n + 1), S, A):
        G = list(H)

        I, J = A[:j - 1], A[j:]

        Hring = f.ring
        for _ in range(j, n):
            Hring = Hring.drop(j + 1)

        x = Hring.gens[0]
        evalpoints = list(zip(LC[0].ring.gens[j:-1], J))

        for i, (h, lc) in enumerate(zip(H, LC)):
            if evalpoints:
                lc = lc.eval(evalpoints)
            lc = _trunc(lc, minpoly, p).set_ring(Hring)
            H[i] = h.set_ring(Hring) + (lc - h.eject(*h.ring.gens[1:]).LC.set_ring(Hring))*x**h.degree()

        m = Hring.gens[j] - a
        M = Hring.one

        c = _trunc(s - functools.reduce(operator.mul, H), minpoly, p)

        dj = s.degree(j)

        for k in range(dj):
            if not c:
                break

            M = M * m
            C = c.diff(x=j, m=k + 1).eval(x=j, a=a)

            if C:
                C = C.quo_ground(ring.domain.factorial(k + 1))  # coeff of (x_{j-1} - a_{j-1})^(k + 1) in c
                T = _diophantine(G, C, I, d, minpoly, p)
                if T is None:
                    raise NotImplementedError

                for i, (h, t) in enumerate(zip(H, T)):
                    H[i] = _trunc(h + t.set_ring(Hring)*M, minpoly, p)

                c = _trunc(s - functools.reduce(operator.mul, H), minpoly, p)

    prod = functools.reduce(operator.mul, H)
    if _trunc(prod, minpoly, p) == f.trunc_ground(p):
        return H


def _sqf_p(f, minpoly, p):
    r"""
    Return ``True`` if `f` is square-free in `\mathbb Z_p[z]/(\mu(z))[x]`.

    """
    ring = f.ring
    lcinv, *_ = _gf_gcdex(f.eject(*ring.gens[1:]).LC, minpoly, p)

    f = _trunc(f * lcinv.set_ring(ring), minpoly, p)

    if not f:
        return True
    else:
        return _euclidean_algorithm(f, _trunc(f.diff(0), minpoly, p), minpoly, p) == 1


def _test_prime(fA, D, minpoly, p, domain):
    r"""
    Test if a prime number is suitable for _factor.

    See also
    ========

    _factor

    """
    if fA.LC % p == 0 or minpoly.LC % p == 0:
        return False
    if not _sqf_p(fA, minpoly, p):
        return False
    if D % p == 0:
        return False

    return True


# squarefree f with cont_x0(f) = 1
def _factor(f, save):
    r"""
    Factor a multivariate polynomial `f`, which is squarefree and primitive
    in `x_0`, in `\mathbb Q(\alpha)[x_0, \ldots, x_n]`.

    References
    ==========

    * :cite:`Javadi2009factor`

    """
    ring = f.ring  # Q(alpha)[x_0, ..., x_{n-1}]
    lcring = ring.drop(0)
    uniring = ring.drop(*ring.gens[1:])
    ground = ring.domain.domain
    n = ring.ngens

    z = Dummy('z')
    qring = ring.clone(symbols=ring.symbols + (z,), domain=ground)
    lcqring = qring.drop(0)

    groundring = ground.ring
    zring = qring.clone(domain=groundring)
    lczring = zring.drop(0)

    minpoly = _minpoly_from_dense(ring.domain.mod, zring.drop(*zring.gens[:-1]))
    f_ = _monic_associate(f, zring)

    if save is True:
        D = minpoly.resultant(minpoly.diff(0))
    else:
        D = groundring.one

    # heuristic bound for p-adic lift
    B = (f_.max_norm() + 1)*D

    lc = f_.eject(*zring.gens[1:]).LC
    gamma, lcfactors = efactor(_z_to_alpha(lc, lcring))  # over QQ(alpha)[x_1, ..., x_n]

    gamma = ground.convert(gamma)
    D_ = gamma.denominator
    gamma_ = gamma.numerator
    lcfactors_ = []

    for l, exp in lcfactors:
        den, l_ = _alpha_to_z(l, lcqring).clear_denoms()  # l_ in QQ[x_1, ..., x_n, z], but coeffs in ZZ
        cont, l_ = l_.set_ring(lczring).primitive()
        D_ *= den**exp
        gamma_ *= cont**exp
        lcfactors_.append((l_, exp))

    f_ *= D_
    p = 2

    N = 0
    history = set()
    tries = 5  # how big should this be?

    while True:
        for _ in range(tries):
            A = tuple(random.randint(-N, N) for _ in range(n - 1))

            if A in history:
                continue
            history.add(A)

            try:
                result = _test_evaluation_points(f_, gamma_, lcfactors, A, D)
            except UnluckyLeadingCoefficient:
                # TODO: check interval
                C = [random.randint(1, 3*(N + 1)) for _ in range(n - 1)]
                gens = zring.gens
                x = gens[0]

                for i, ci in zip(range(1, n + 1), C):
                    xi = gens[i]
                    f_ = f_.compose(xi, x + xi*ci)

                lc, factors = _factor(_z_to_alpha(f_, ring), save)
                gens = factors[0].ring.gens
                x = gens[0]

                for i, ci in zip(range(1, n + 1), C):
                    xi = gens[i]
                    factors = [g.compose(xi, (xi - x).quo_ground(ci)) for g in factors]

                return lc, factors

            if result is None:
                continue
            fA, denoms, divisors = result

            with using(aa_factor_method='trager'):
                _, fAfactors = _z_to_alpha(fA, uniring).factor_list()
            if len(fAfactors) == 1:
                g = _z_to_alpha(f_, ring)
                return f.LC, [g.monic()]

            result = _leading_coeffs(f_, fAfactors, gamma_, lcfactors_, A, D, denoms, divisors)
            if result is None:
                continue
            f_, lcs, fAfactors_ = result

            prod = groundring.one
            for lc in lcs:
                prod *= lc.LC
            delta = (ground(prod, f_.LC)).numerator

            f_ *= delta

            while not _test_prime(fA, D, minpoly, p, zring.domain):
                p = nextprime(p)

            pfactors = _hensel_lift(f_, fAfactors_, lcs, A, minpoly, p)
            if pfactors is None:
                p = nextprime(p)
                f_ = f_.primitive()[1]
                continue

            factors = _padic_lift(f_, pfactors, lcs, B, minpoly, p)
            if factors is None:
                p = nextprime(p)
                f_ = f_.primitive()[1]
                B *= B
                continue

            return f.LC, [_z_to_alpha(g.primitive()[1], ring).monic() for g in factors]

        N += 1


# output of the form (lc, [(poly1, exp1), ...])
def efactor(f, save=True):
    r"""Factor a multivariate polynomial `f` in `\mathbb Q(\alpha)[x_0, \ldots, x_n]`.

    By default, an estimate of the defect of the algebraic field is included
    in all computations. If ``save`` is set to ``False``, the defect will be
    treated as one, thus computations are faster. However, if the defect of
    `\alpha` is larger than one, this may lead to wrong results.

    References
    ==========

    * :cite:`Javadi2009factor`

    """
    ring = f.ring

    assert ring.domain.is_AlgebraicField

    if f.is_ground:
        return f[1], []

    n = ring.ngens

    if n == 1:
        with using(aa_factor_method='trager'):
            return f.factor_list()
    else:
        cont, f = f.eject(*ring.gens[1:]).primitive()
        f = f.inject()
        if cont != 1:
            lccont, contfactors = efactor(cont)
            lc, factors = efactor(f)
            contfactors = [(g.set_ring(ring), exp) for g, exp in contfactors]
            return lccont * lc, _sort_factors(contfactors + factors)

        # this is only correct because the content in x_0 is already divided out
        lc, sqflist = f.sqf_list()
        factors = []
        for g, exp in sqflist:
            lcg, gfactors = _factor(g, save)
            lc *= lcg
            factors = factors + [(gi, exp) for gi in gfactors]

        return lc, _sort_factors(factors)
