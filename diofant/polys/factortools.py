"""Polynomial factorization routines in characteristic zero."""

import functools
import math

from ..ntheory import factorint, isprime, nextprime
from ..ntheory.modular import symmetric_residue
from ..utilities import subsets
from .densearith import (dmp_add, dmp_add_mul, dmp_div, dmp_max_norm, dmp_mul,
                         dmp_mul_ground, dmp_neg, dmp_pow, dmp_quo,
                         dmp_quo_ground, dmp_rem, dmp_sub, dmp_sub_mul,
                         dup_add, dup_lshift, dup_mul, dup_sqr, dup_sub)
from .densebasic import (dmp_convert, dmp_degree_in, dmp_degree_list,
                         dmp_ground_LC, dmp_ground_p, dmp_LC, dmp_nest,
                         dmp_normal, dmp_one, dmp_raise, dmp_strip, dmp_TC,
                         dmp_zero_p, dup_inflate)
from .densetools import (dmp_compose, dmp_diff_eval_in, dmp_eval_in,
                         dmp_eval_tail, dmp_ground_content, dmp_ground_monic,
                         dmp_ground_primitive, dmp_ground_trunc)
from .euclidtools import dup_gcdex
from .galoistools import dup_gf_factor_sqf
from .polyconfig import query
from .polyerrors import (CoercionFailed, DomainError, EvaluationFailed,
                         ExtraneousFactors)
from .polyutils import _sort_factors
from .sqfreetools import dmp_sqf_list, dmp_sqf_p, dmp_sqf_part


def dmp_trial_division(f, factors, u, K):
    """Determine multiplicities of factors using trial division."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    factors = list(map(ring.from_list, factors))
    result = ring._trial_division(f, factors)
    return [(f.to_dense(), k) for f, k in result]


def dmp_zz_mignotte_bound(f, u, K):
    """Mignotte bound for multivariate polynomials in `K[X]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    return ring._zz_mignotte_bound(f)


def dup_zz_irreducible_p(f, K):
    """Test irreducibility using Eisenstein's criterion."""
    lc = dmp_LC(f, K)
    tc = dmp_TC(f, K)

    e_fc = dmp_ground_content(f[1:], 0, K)

    if e_fc:
        e_ff = factorint(int(e_fc))

        for p in e_ff:
            if (lc % p) and (tc % p**2):
                return True


def dup_cyclotomic_p(f, K, irreducible=False):
    """
    Efficiently test if ``f`` is a cyclotomic polynomial.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> (x**16 + x**14 - x**10 + x**8 - x**6 + x**2 + 1).is_cyclotomic
    False

    >>> (x**16 + x**14 - x**10 - x**8 - x**6 + x**2 + 1).is_cyclotomic
    True

    """
    if K.is_RationalField:
        try:
            K0, K = K, K.ring
            f = dmp_convert(f, 0, K0, K)
        except CoercionFailed:
            return False
    elif not K.is_IntegerRing:
        return False

    lc = dmp_LC(f, K)
    tc = dmp_TC(f, K)

    if lc != 1 or (tc != -1 and tc != 1):
        return False

    if not irreducible:
        coeff, factors = dmp_factor_list(f, 0, K)

        if coeff != K.one or factors != [(f, 1)]:
            return False

    n = dmp_degree_in(f, 0, 0)
    g, h = [], []

    for i in range(n, -1, -2):
        g.insert(0, f[i])

    for i in range(n - 1, -1, -2):
        h.insert(0, f[i])

    g = dup_sqr(dmp_strip(g, 0), K)
    h = dup_sqr(dmp_strip(h, 0), K)

    F = dup_sub(g, dup_lshift(h, 1, K), K)

    if dmp_LC(F, K) < 0:
        F = dmp_neg(F, 0, K)

    if F == f:
        return True

    g = dmp_compose(f, [-K.one, 0], 0, K)

    if dmp_LC(g, K) < 0:
        g = dmp_neg(g, 0, K)

    if F == g and dup_cyclotomic_p(g, K):
        return True

    G = dmp_sqf_part(F, 0, K)

    if dup_sqr(G, K) == F and dup_cyclotomic_p(G, K):
        return True

    return False


def dup_zz_cyclotomic_poly(n, K):
    """Efficiently generate n-th cyclotomic polynomial."""
    h = [K.one, -K.one]

    for p, k in factorint(n).items():
        h = dmp_quo(dup_inflate(h, p, K), h, 0, K)
        h = dup_inflate(h, p**(k - 1), K)

    return h


def _dup_cyclotomic_decompose(n, K):
    H = [[K.one, -K.one]]

    for p, k in factorint(n).items():
        Q = [dmp_quo(dup_inflate(h, p, K), h, 0, K) for h in H]
        H.extend(Q)

        for i in range(1, k):
            Q = [dup_inflate(q, p, K) for q in Q]
            H.extend(Q)

    return H


def dup_zz_cyclotomic_factor(f, K):
    """
    Efficiently factor polynomials `x**n - 1` and `x**n + 1` in `Z[x]`.

    Given a univariate polynomial `f` in `Z[x]` returns a list of factors
    of `f`, provided that `f` is in the form `x**n - 1` or `x**n + 1` for
    `n >= 1`. Otherwise returns None.

    Factorization is performed using cyclotomic decomposition of `f`,
    which makes this method much faster that any other direct factorization
    approach (e.g. Zassenhaus's).

    References
    ==========

    * :cite:`MathWorld-Cyclotomic-Poly`

    """
    lc_f, tc_f = dmp_LC(f, K), dmp_TC(f, K)

    if dmp_ground_p(f, None, 0):
        return

    if lc_f != 1 or tc_f not in [-1, 1]:
        return

    if any(bool(cf) for cf in f[1:-1]):
        return

    n = dmp_degree_in(f, 0, 0)
    F = _dup_cyclotomic_decompose(n, K)

    if tc_f != K.one:
        return F
    else:
        H = []

        for h in _dup_cyclotomic_decompose(2*n, K):
            if h not in F:
                H.append(h)

        return H


def dup_zz_factor_sqf(f, K):
    """Factor square-free (non-primitive) polynomials in `Z[x]`."""
    ring = K.poly_ring('_0')
    f = ring.from_list(f)
    cont, factors = ring._zz_factor_sqf(f)
    return cont, [_.to_dense() for _ in factors]


def dmp_zz_wang_non_divisors(E, cs, ct, K):
    """Wang/EEZ: Compute a set of valid divisors."""
    result = [cs*ct]

    for q in E:
        q = abs(q)

        for r in reversed(result):
            while r != 1:
                r = K.gcd(r, q)
                q = q // r

            if q == K.one:
                return

        result.append(q)

    return result[1:]


def dmp_zz_wang_test_points(f, T, ct, A, u, K):
    """Wang/EEZ: Test evaluation points for suitability."""
    if not dmp_eval_tail(dmp_LC(f, K), A, u - 1, K):
        raise EvaluationFailed('no luck')

    g = dmp_eval_tail(f, A, u, K)

    if not dmp_sqf_p(g, 0, K):
        raise EvaluationFailed('no luck')

    c, h = dmp_ground_primitive(g, 0, K)

    v = u - 1

    E = [dmp_eval_tail(t, A, v, K) for t, _ in T]
    D = dmp_zz_wang_non_divisors(E, c, ct, K)

    if D is not None:
        return c, h, E
    else:
        raise EvaluationFailed('no luck')


def dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K):
    """Wang/EEZ: Compute correct leading coefficients."""
    C, J, v = [], [0]*len(E), u - 1

    for h in H:
        c = dmp_one(v, K)
        d = dmp_LC(h, K)*cs

        for i in reversed(range(len(E))):
            k, e, t = 0, E[i], T[i][0]

            while not (d % e):
                d, k = d//e, k + 1

            if k != 0:
                c, J[i] = dmp_mul(c, dmp_pow(t, k, v, K), v, K), 1

        C.append(c)

    if any(not j for j in J):  # pragma: no cover
        raise ExtraneousFactors

    CC, HH = [], []

    for c, h in zip(C, H):
        d = dmp_eval_tail(c, A, v, K)
        lc = dmp_LC(h, K)

        if cs == K.one:
            cc = lc//d
        else:
            g = K.gcd(lc, d)
            d, cc = d//g, lc//g
            h, cs = dmp_mul_ground(h, d, 0, K), cs//d

        c = dmp_mul_ground(c, cc, v, K)

        CC.append(c)
        HH.append(h)

    if cs == K.one:
        return f, HH, CC

    CCC, HHH = [], []

    for c, h in zip(CC, HH):
        CCC.append(dmp_mul_ground(c, cs, v, K))
        HHH.append(dmp_mul_ground(h, cs, 0, K))

    f = dmp_mul_ground(f, cs**(len(H) - 1), u, K)

    return f, HHH, CCC


def dup_zz_diophantine(F, m, p, K):
    """Wang/EEZ: Solve univariate Diophantine equations."""
    if len(F) == 2:
        Kp = K.finite_field(p)
        f, g = map(lambda x: dmp_normal(x, 0, Kp), F)

        s, t, _ = dup_gcdex(g, f, Kp)

        s = dup_lshift(s, m, Kp)
        t = dup_lshift(t, m, Kp)

        q, s = dmp_div(s, f, 0, Kp)
        s = dmp_normal(s, 0, K)

        t = dmp_add_mul(t, q, g, 0, Kp)
        t = dmp_normal(t, 0, K)

        result = [s, t]
    else:
        G = [F[-1]]

        for f in reversed(F[1:-1]):
            G.insert(0, dup_mul(f, G[0], K))

        S, T = [], [[1]]

        for f, g in zip(F, G):
            t, s = dmp_zz_diophantine([g, f], T[-1], [], 0, p, 1, K)
            T.append(t)
            S.append(s)

        result, S = [], S + [T[-1]]
        Kp = K.finite_field(p)

        for s, f in zip(S, F):
            s = dmp_rem(dup_lshift(s, m, K), f, 0, Kp)
            s = dmp_normal(s, 0, K)

            result.append(s)

    return result


def dmp_zz_diophantine(F, c, A, d, p, u, K):
    """Wang/EEZ: Solve multivariate Diophantine equations."""
    if not A:
        S = [[] for _ in F]
        n = dmp_degree_in(c, 0, 0)

        for i, coeff in enumerate(c):
            if not coeff:
                continue

            T = dup_zz_diophantine(F, n - i, p, K)

            for j, (s, t) in enumerate(zip(S, T)):
                t = dmp_mul_ground(t, coeff, 0, K)
                S[j] = dmp_ground_trunc(dup_add(s, t, K), p, 0, K)
    else:
        n = len(A)
        e = functools.reduce(lambda x, y: dmp_mul(x, y, u, K), F)

        a, A = A[-1], A[:-1]
        B, G = [], []

        for f in F:
            B.append(dmp_quo(e, f, u, K))
            G.append(dmp_eval_in(f, a, n, u, K))

        C = dmp_eval_in(c, a, n, u, K)

        v = u - 1

        S = dmp_zz_diophantine(G, C, A, d, p, v, K)
        S = [dmp_raise(s, 1, v, K) for s in S]

        for s, b in zip(S, B):
            c = dmp_sub_mul(c, s, b, u, K)

        c = dmp_ground_trunc(c, p, u, K)

        m = dmp_nest([K.one, -a], n, K)
        M = dmp_one(n, K)

        for k in range(d):
            k = K(k)
            if dmp_zero_p(c, u):
                break

            M = dmp_mul(M, m, u, K)
            C = dmp_diff_eval_in(c, k + 1, a, n, u, K)

            if not dmp_zero_p(C, v):
                C = dmp_quo_ground(C, K.factorial(k + 1), v, K)
                T = dmp_zz_diophantine(G, C, A, d, p, v, K)

                for i, t in enumerate(T):
                    T[i] = dmp_mul(dmp_raise(t, 1, v, K), M, u, K)

                for i, (s, t) in enumerate(zip(S, T)):
                    S[i] = dmp_add(s, t, u, K)

                for t, b in zip(T, B):
                    c = dmp_sub_mul(c, t, b, u, K)

                c = dmp_ground_trunc(c, p, u, K)

        S = [dmp_ground_trunc(s, p, u, K) for s in S]

    return S


def dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K):
    """Wang/EEZ: Parallel Hensel lifting algorithm."""
    S, n, v = [f], len(A), u - 1

    H = list(H)

    for i, a in enumerate(reversed(A[1:])):
        s = dmp_eval_in(S[0], a, n - i, u - i, K)
        S.insert(0, dmp_ground_trunc(s, p, v - i, K))

    d = max(dmp_degree_list(f, u)[1:])

    for j, s, a in zip(range(2, n + 2), S, A):
        G, w = list(H), j - 1

        I, J = A[:j - 2], A[j - 1:]

        for i, (h, lc) in enumerate(zip(H, LC)):
            lc = dmp_ground_trunc(dmp_eval_tail(lc, J, v, K), p, w - 1, K)
            H[i] = [lc] + dmp_raise(h[1:], 1, w - 1, K)

        m = dmp_nest([K.one, -a], w, K)
        M = dmp_one(w, K)

        c = functools.reduce(lambda x, y: dmp_mul(x, y, w, K), H)
        c = dmp_sub(s, c, w, K)

        dj = dmp_degree_in(s, w, w)

        for k in range(dj):
            k = K(k)
            if dmp_zero_p(c, w):
                break

            M = dmp_mul(M, m, w, K)
            C = dmp_diff_eval_in(c, k + 1, a, w, w, K)

            if not dmp_zero_p(C, w - 1):
                C = dmp_quo_ground(C, K.factorial(k + 1), w - 1, K)
                T = dmp_zz_diophantine(G, C, I, d, p, w - 1, K)

                for i, (h, t) in enumerate(zip(H, T)):
                    h = dmp_add_mul(h, dmp_raise(t, 1, w - 1, K), M, w, K)
                    H[i] = dmp_ground_trunc(h, p, w, K)

                h = functools.reduce(lambda x, y: dmp_mul(x, y, w, K), H)
                h = dmp_sub(s, h, w, K)
                c = dmp_ground_trunc(h, p, w, K)

    if functools.reduce(lambda x, y: dmp_mul(x, y, u, K), H) != f:
        raise ExtraneousFactors  # pragma: no cover
    else:
        return H


def dmp_zz_wang(f, u, K, mod=None, seed=None):
    """
    Factor primitive square-free polynomials in `Z[X]`.

    Given a multivariate polynomial `f` in `Z[x_1,...,x_n]`, which is
    primitive and square-free in `x_1`, computes factorization of `f` into
    irreducibles over integers.

    The procedure is based on Wang's Enhanced Extended Zassenhaus
    algorithm. The algorithm works by viewing `f` as a univariate polynomial
    in `Z[x_2,...,x_n][x_1]`, for which an evaluation mapping is computed::

                      x_2 -> a_2, ..., x_n -> a_n

    where `a_i`, for `i = 2, ..., n`, are carefully chosen integers.  The
    mapping is used to transform `f` into a univariate polynomial in `Z[x_1]`,
    which can be factored efficiently using Zassenhaus algorithm. The last
    step is to lift univariate factors to obtain true multivariate
    factors. For this purpose a parallel Hensel lifting procedure is used.

    The parameter ``seed`` is passed to _randint and can be used to seed randint
    (when an integer) or (for testing purposes) can be a sequence of numbers.

    References
    ==========

    * :cite:`Wang1978improved`
    * :cite:`Geddes1992algorithms`

    """
    from ..utilities.randtest import _randint

    randint = _randint(seed)

    ct, T = dmp_zz_factor(dmp_LC(f, K), u - 1, K)

    b = dmp_zz_mignotte_bound(f, u, K)
    p = K(nextprime(b))

    if mod is None:
        if u == 1:
            mod = 2
        else:
            mod = 1

    history, configs, A, r = set(), [], [K.zero]*u, None

    try:
        cs, s, E = dmp_zz_wang_test_points(f, T, ct, A, u, K)

        _, H = dup_zz_factor_sqf(s, K)

        r = len(H)

        if r == 1:
            return [f]

        configs = [(s, cs, E, H, A)]
    except EvaluationFailed:
        pass

    eez_num_configs = query('EEZ_NUMBER_OF_CONFIGS')
    eez_num_tries = query('EEZ_NUMBER_OF_TRIES')
    eez_mod_step = query('EEZ_MODULUS_STEP')

    while len(configs) < eez_num_configs:
        for _ in range(eez_num_tries):
            A = [K(randint(-mod, mod)) for _ in range(u)]

            if tuple(A) in history:
                continue
            else:
                history.add(tuple(A))

            try:
                cs, s, E = dmp_zz_wang_test_points(f, T, ct, A, u, K)
            except EvaluationFailed:
                continue

            _, H = dup_zz_factor_sqf(s, K)

            rr = len(H)

            if r is not None:
                if rr != r:  # pragma: no cover
                    if rr < r:
                        configs, r = [], rr
                    else:
                        continue
            else:
                r = rr

            if r == 1:
                return [f]

            configs.append((s, cs, E, H, A))

            if len(configs) == eez_num_configs:
                break
        else:
            mod += eez_mod_step

    s_norm, s_arg, i = None, 0, 0

    for s, _, _, _, _ in configs:
        _s_norm = dmp_max_norm(s, 0, K)

        if s_norm is not None:
            if _s_norm < s_norm:
                s_norm = _s_norm
                s_arg = i
        else:
            s_norm = _s_norm

        i += 1

    _, cs, E, H, A = configs[s_arg]
    orig_f = f

    try:
        f, H, LC = dmp_zz_wang_lead_coeffs(f, T, cs, E, H, A, u, K)
        factors = dmp_zz_wang_hensel_lifting(f, H, LC, A, p, u, K)
    except ExtraneousFactors:  # pragma: no cover
        if query('EEZ_RESTART_IF_NEEDED'):
            return dmp_zz_wang(orig_f, u, K, mod + 1)
        else:
            raise ExtraneousFactors(
                'we need to restart algorithm with better parameters')

    result = []

    for f in factors:
        _, f = dmp_ground_primitive(f, u, K)

        result.append(f)

    return result


def dmp_zz_factor(f, u, K):
    """Factor (non square-free) polynomials in `Z[X]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    lc, factors = ring._zz_factor(f)
    return lc, [(f.to_dense(), k) for f, k in factors]


def dmp_gf_factor(f, u, K):
    """Factor multivariate polynomials over finite fields."""
    if u:
        raise NotImplementedError('multivariate polynomials over finite fields')
    else:
        lc = dmp_ground_LC(f, u, K)
        f = dmp_ground_monic(f, u, K)

        factors = []

        for g, n in dmp_sqf_list(f, 0, K)[1]:
            for h in dup_gf_factor_sqf(g, K):
                factors.append((h, n))

        return lc, factors


def dmp_factor_list(f, u, K):
    """Factor polynomials into irreducibles in `K[X]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    lc, factors = ring.factor_list(f)
    return lc, [(f.to_dense(), k) for f, k in factors]


class _Factor:
    """Mixin class for factorization routines."""

    def _trial_division(self, f, factors):
        result = []

        for factor in factors:
            k = 0

            while f:
                q, r = divmod(f, factor)

                if r.is_zero:
                    f, k = q, k + 1
                else:
                    break

            result.append((factor, k))

        return _sort_factors(result)

    def _aa_factor_trager(self, f):
        """Factor multivariate polynomials over algebraic number fields."""
        domain = self.domain

        lc, f = f.LC, f.monic()

        if f.is_ground:
            return lc, []

        f, F = f.sqf_part(), f
        s, g, r = f.sqf_norm()

        _, factors = r.factor_list()

        if len(factors) == 1:
            factors = [f]
        else:
            for i, (factor, _) in enumerate(factors):
                h = factor.set_domain(domain)
                h, _, g = self.cofactors(h, g)
                h = h.compose({x: x + s*domain.unit for x in self.gens})
                factors[i] = h

        return lc, self._trial_division(F, factors)

    def _zz_factor(self, f):
        """
        Factor (non square-free) polynomials in `Z[X]`.

        Given a multivariate polynomial `f` in `Z[X]` computes its complete
        factorization `f_1, ..., f_n` into irreducibles over integers::

            f = content(f) f_1**k_1 ... f_n**k_n

        The factorization is computed by reducing the input polynomial
        into a primitive square-free polynomial and factoring it using
        Zassenhaus or Enhanced Extended Zassenhaus (EEZ) algorithm. Trial
        division is used to recover the multiplicities of factors.

        The result is returned as a tuple consisting of::

            (content(f), [(f_1, k_1), ..., (f_n, k_n))

        Examples
        ========

        >>> R, x = ring('x', ZZ)

        >>> (2*x**4 - 2).factor_list()
        (2, [(x - 1, 1), (x + 1, 1), (x**2 + 1, 1)])

        >>> R, x, y = ring('x y', ZZ)

        >>> (2*x**2 - 2*y**2).factor_list()
        (2, [(x - y, 1), (x + y, 1)])

        References
        ==========

        * :cite:`Gathen1999modern`

        """
        domain = self.domain

        if self.is_univariate:
            cont, g = f.primitive()

            n = g.degree()

            if n <= 0:
                return cont, []
            elif n == 1:
                return cont, [(g, 1)]

            if query('USE_IRREDUCIBLE_IN_FACTOR'):
                if self.dup_zz_irreducible_p(g):
                    return cont, [(g, 1)]

            g = g.sqf_part()
            H = None

            if query('USE_CYCLOTOMIC_FACTOR'):
                H = self.dup_zz_cyclotomic_factor(g)

            if H is None:
                H = self._zz_zassenhaus(g)

            factors = self._trial_division(f, H)
            return cont, factors

        if f.is_zero:
            return domain.zero, []

        cont, g = f.primitive()

        if g.is_ground:
            return cont, []

        g = g.eject(*self.gens[1:])
        G, g = g.primitive()
        g = g.inject()

        factors = []

        if g.degree() > 0:
            g = g.sqf_part()
            H = self.dmp_zz_wang(g)
            factors = self._trial_division(f, H)

        new_ring = G.ring

        for g, k in new_ring._zz_factor(G)[1]:
            factors.insert(0, (g.set_ring(self), k))

        return cont, _sort_factors(factors)

    def _zz_mignotte_bound(self, f):
        """Mignotte bound for multivariate polynomials in `Z[X]`."""
        domain = self.domain

        a = f.max_norm()
        b = abs(f.LC)
        n = sum(f.degree_list())

        return domain.sqrt(domain(n + 1))*2**n*a*b

    def _zz_zassenhaus(self, f):
        """Factor primitive square-free polynomial in `Z[x]`."""
        domain = self.domain

        assert self.is_univariate and domain.is_IntegerRing

        n = f.degree()

        if n == 1:
            return [f]

        fc = f.coeff(1)
        A = f.max_norm()
        b = f.LC
        B = int(self._zz_mignotte_bound(f))
        C = int((n + 1)**(2*n)*A**(2*n - 1))
        gamma = math.ceil(2*math.log(C, 2))
        bound = int(2*gamma*math.log(gamma))
        a = []

        # choose a prime number `p` such that `f` be square free in Z_p
        # if there are many factors in Z_p, choose among a few different `p`
        # the one with fewer factors
        for p in range(3, bound + 1):  # pragma: no branch
            if not isprime(p) or b % p == 0:
                continue

            p = domain.convert(p)
            p_domain = domain.finite_field(p)

            F = f.set_domain(p_domain)

            if not F.is_squarefree:
                continue

            F = F.monic()
            fsqfx = [_ for _, k in F.factor_list()[1]]

            a.append((p, [_.set_domain(domain) for _ in fsqfx]))
            if len(fsqfx) < 15 or len(a) > 4:
                break
        p, fsqf = min(a, key=lambda x: len(x[1]))

        l = math.ceil(math.log(2*B + 1, p))
        g = self._zz_hensel_lift(p, f, fsqf, l)

        sorted_T = range(len(g))
        T = set(sorted_T)
        factors, s = [], 1
        pl = p**l

        while 2*s <= len(T):
            for S in subsets(sorted_T, s):
                # lift the constant coefficient of the product `G` of the factors
                # in the subset `S`; if it is does not divide `fc`, `G` does
                # not divide the input polynomial

                if b == 1:
                    q = 1
                    for i in S:
                        q = q*g[i].coeff(1)
                    q = q % pl
                    qs = symmetric_residue(q, pl)
                    if qs and fc % qs != 0:
                        continue
                else:
                    G = self.ground_new(b)
                    for i in S:
                        G *= g[i]
                    G = G.trunc_ground(pl)
                    _, G = G.primitive()
                    q = G.coeff(1)
                    if q and fc % q != 0:
                        continue

                H = self.ground_new(b)
                S = set(S)
                T_S = T - S

                if b == 1:
                    G = self.ground_new(b)
                    for i in S:
                        G *= g[i]
                    G = G.trunc_ground(pl)

                for i in T_S:
                    H *= g[i]

                H = H.trunc_ground(pl)

                G_norm = G.l1_norm()
                H_norm = H.l1_norm()

                if G_norm*H_norm <= B:
                    T = T_S
                    sorted_T = [i for i in sorted_T if i not in S]

                    G = G.primitive()[1]
                    f = H.primitive()[1]

                    factors.append(G)
                    b = f.LC

                    break
            else:
                s += 1

        return factors + [f]

    def _zz_factor_sqf(self, f):
        """Factor square-free (non-primitive) polynomials in `Z[x]`."""
        cont, g = f.primitive()

        n = g.degree()

        if n <= 0:
            return cont, []
        elif n == 1:
            return cont, [g]

        if query('USE_IRREDUCIBLE_IN_FACTOR'):
            if self.dup_zz_irreducible_p(g):
                return cont, [g]

        factors = None

        if query('USE_CYCLOTOMIC_FACTOR'):
            factors = self.dup_zz_cyclotomic_factor(g)

        if factors is None:
            factors = self._zz_zassenhaus(g)

        return cont, _sort_factors(factors, multiple=False)

    def _zz_hensel_lift(self, p, f, f_list, l):
        """
        Multifactor Hensel lifting in `Z[x]`.

        Given a prime `p`, polynomial `f` over `Z[x]` such that `lc(f)`
        is a unit modulo `p`, monic pair-wise coprime polynomials `f_i`
        over `Z[x]` satisfying::

            f = lc(f) f_1 ... f_r (mod p)

        and a positive integer `l`, returns a list of monic polynomials
        `F_1`, `F_2`, ..., `F_r` satisfying::

           f = lc(f) F_1 ... F_r (mod p**l)

           F_i = f_i (mod p), i = 1..r

        References
        ==========

        * :cite:`Gathen1999modern`

        """
        domain = self.domain

        r = len(f_list)
        lc = f.LC

        if r == 1:
            F = f * domain.gcdex(lc, p**l)[0]
            return [F.trunc_ground(p**l)]

        m = p
        k = r // 2
        d = math.ceil(math.log(l, 2))
        p_domain = domain.finite_field(p)

        g = self.ground_new(lc)
        g = g.set_domain(p_domain)

        for f_i in f_list[:k]:
            g *= f_i

        h = f_list[k].set_domain(p_domain)

        for f_i in f_list[k + 1:]:
            h *= f_i

        s, t, _ = self.clone(domain=p_domain).gcdex(g, h)

        g, h, s, t = map(lambda x: x.set_domain(domain), (g, h, s, t))

        for _ in range(1, d + 1):
            (g, h, s, t), m = self._zz_hensel_step(m, f, g, h, s, t), m**2

        return (self._zz_hensel_lift(p, g, f_list[:k], l) +
                self._zz_hensel_lift(p, h, f_list[k:], l))

    def _zz_hensel_step(self, m, f, g, h, s, t):
        """
        One step in Hensel lifting in `Z[x]`.

        Given positive integer `m` and `Z[x]` polynomials `f`, `g`, `h`, `s`
        and `t` such that::

            f == g*h (mod m)
            s*g + t*h == 1 (mod m)

            lc(f) is not a zero divisor (mod m)
            lc(h) == 1

            deg(f) == deg(g) + deg(h)
            deg(s) < deg(h)
            deg(t) < deg(g)

        returns polynomials `G`, `H`, `S` and `T`, such that::

            f == G*H (mod m**2)
            S*G + T*H == 1 (mod m**2)

        References
        ==========

        * :cite:`Gathen1999modern`

        """
        M = m**2

        e = f - g*h
        e = e.trunc_ground(M)

        q, r = divmod(s*e, h)

        q = q.trunc_ground(M)
        r = r.trunc_ground(M)

        u = t*e + q*g
        G = (g + u).trunc_ground(M)
        H = (h + r).trunc_ground(M)

        u = s*G + t*H
        b = (u - 1).trunc_ground(M)

        c, d = divmod(s*b, H)

        c = c.trunc_ground(M)
        d = d.trunc_ground(M)

        u = t*b + c*G
        S = (s - d).trunc_ground(M)
        T = (t - u).trunc_ground(M)

        return G, H, S, T

    def factor_list(self, f):
        """Factor polynomials into irreducibles in `K[X]`."""
        domain = self.domain

        cont, f = f.primitive()

        if domain.is_FiniteField:
            coeff, factors = dmp_gf_factor(f.to_dense(), self.ngens-1, domain)
            factors = [(self.from_list(_), k) for _, k in factors]
        elif domain.is_AlgebraicField:
            coeff, factors = self._aa_factor_trager(f)
        else:
            if not domain.is_Exact:
                domain_inexact, domain = domain, domain.get_exact()
                f = f.set_domain(domain)
            else:
                domain_inexact = None

            if domain.is_Field:
                domain1 = domain.ring

                denom, f = f.clear_denoms(convert=True)
            else:
                domain1 = domain

            if domain1.is_IntegerRing:
                ring = self.clone(domain=domain1)
                coeff, factors = ring._zz_factor(f)
            elif domain1.is_PolynomialRing:
                f = f.inject()
                ring = f.ring

                coeff, factors = ring.factor_list(f)

                for i, (f, k) in enumerate(factors):
                    factors[i] = (f.eject(*ring.gens[-domain1.ngens:]), k)

                coeff = domain1.convert(coeff)
            else:  # pragma: no cover
                raise DomainError(f'factorization not supported over {domain}')

            if domain.is_Field:
                for i, (f, k) in enumerate(factors):
                    factors[i] = (f.set_domain(domain), k)

                coeff = domain.convert(coeff)
                coeff = domain.quo(coeff, denom)

                if domain_inexact:
                    for i, (f, k) in enumerate(factors):
                        f = f.quo_ground(denom)
                        f = f.set_domain(domain_inexact)
                        factors[i] = (f, k)
                        coeff *= denom**k

                    coeff = domain_inexact.convert(coeff)

        return coeff*cont, _sort_factors(factors)
