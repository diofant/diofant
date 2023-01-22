"""Polynomial factorization routines in characteristic zero."""

import itertools
import math
import operator

from ..config import query
from ..ntheory import factorint, isprime, nextprime
from ..ntheory.modular import symmetric_residue
from .polyerrors import (CoercionFailedError, DomainError,
                         EvaluationFailedError, ExtraneousFactorsError)
from .polyutils import _sort_factors


class _Factor:
    """Mixin class for factorization routines."""

    def _trial_division(self, f, factors):
        result = []

        for factor in factors:
            k = 0

            while f:
                q, r = divmod(f, factor)

                if not r:
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

        >>> _, x = ring('x', ZZ)

        >>> (2*x**4 - 2).factor_list()
        (2, [(x - 1, 1), (x + 1, 1), (x**2 + 1, 1)])

        >>> _, x, y = ring('x y', ZZ)

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
            if n == 1:
                return cont, [(g, 1)]

            if query('USE_IRREDUCIBLE_IN_FACTOR'):
                if self._zz_irreducible_p(g):
                    return cont, [(g, 1)]

            g = g.sqf_part()
            H = None

            if query('USE_CYCLOTOMIC_FACTOR'):
                H = self._zz_cyclotomic_factor(g)

            if H is None:
                H = self._zz_zassenhaus(g)

            factors = self._trial_division(f, H)
            return cont, factors

        if not f:
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
            H = self._zz_wang(g)
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
        n = sum(f.degree(x) for x in self.gens)

        return domain.sqrt(domain(n + 1))*2**n*a*b

    def _zz_zassenhaus(self, f):
        """Factor primitive square-free polynomial in `Z[x]`."""
        domain = self.domain

        assert self.is_univariate
        assert domain.is_IntegerRing

        n = f.degree()

        if n == 1:
            return [f]

        fc = f[1]
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
            for S in itertools.combinations(sorted_T, s):
                # lift the constant coefficient of the product `G` of the factors
                # in the subset `S`; if it is does not divide `fc`, `G` does
                # not divide the input polynomial

                if b == 1:
                    q = 1
                    for i in S:
                        q = q*g[i][1]
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
                    q = G[1]
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
        domain = self.domain

        assert self.is_univariate
        assert domain.is_IntegerRing

        cont, g = f.primitive()

        n = g.degree()

        if n <= 0:
            return cont, []
        if n == 1:
            return cont, [g]

        if query('USE_IRREDUCIBLE_IN_FACTOR'):
            if self._zz_irreducible_p(g):
                return cont, [g]

        factors = None

        if query('USE_CYCLOTOMIC_FACTOR'):
            factors = self._zz_cyclotomic_factor(g)

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

        g, h, s, t = map(operator.methodcaller('set_domain', domain),
                         (g, h, s, t))

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
            if self.is_multivariate:
                raise NotImplementedError('multivariate polynomials over finite fields')

            coeff = f.LC
            f = f.monic()

            factors = []

            for g, n in f.sqf_list()[1]:
                for h in self._gf_factor_sqf(g):
                    factors.append((h, n))
        elif domain.is_AlgebraicField:
            from .factorization_alg_field import efactor

            _factor_aa_methods = {'trager': self._aa_factor_trager,
                                  'modular': efactor}
            method = _factor_aa_methods[query('AA_FACTOR_METHOD')]
            coeff, factors = method(f)
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
            else:
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

    def _zz_irreducible_p(self, f):
        """Test irreducibility using Eisenstein's criterion."""
        assert self.is_univariate

        lc = f.LC
        tc = f[1]

        f -= f.leading_term()
        e_fc = f.content()

        if e_fc:
            e_ff = factorint(int(e_fc))

            for p in e_ff:
                if lc % p and tc % p**2:
                    return True

    def _gf_irreducible_p_ben_or(self, f):
        """
        Ben-Or's polynomial irreducibility test over finite fields.

        References
        ==========

        * :cite:`Ben-Or1981ff`

        """
        assert self.is_univariate

        domain = self.domain
        n, q = f.degree(), domain.order

        if n <= 1:
            return True

        x = self.gens[0]
        f = f.monic()

        H = h = pow(x, q, f)

        for _ in range(n//2):
            g = h - x

            if self.gcd(f, g) == 1:
                h = h.compose(x, H) % f
            else:
                return False

        return True

    def _gf_irreducible_p_rabin(self, f):
        """
        Rabin's polynomial irreducibility test over finite fields.

        References
        ==========

        * :cite:`Gathen1999modern`, algorithm 14.36

        """
        assert self.is_univariate

        domain = self.domain
        n, q = f.degree(), domain.order

        if n <= 1:
            return True

        x = self.gens[0]
        f = f.monic()

        indices = {n//d for d in factorint(n)}

        H = h = pow(x, q, f)

        for i in range(1, n):
            if i in indices:
                g = h - x

                if self.gcd(f, g) != 1:
                    return False

            h = h.compose(x, H) % f

        return h == x

    def _cyclotomic_decompose(self, n):
        x = self.gens[0]
        H = [x - 1]

        for p, k in factorint(n).items():
            Q = [h.compose(x, x**p) // h for h in H]
            H.extend(Q)

            for _ in range(1, k):
                Q = [q.compose(x, x**p) for q in Q]
                H.extend(Q)

        return H

    def _zz_cyclotomic_factor(self, f):
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
        lc_f, tc_f = f.LC, f[1]

        if f.is_ground:
            return

        if lc_f != 1 or tc_f not in [-1, 1]:
            return

        if any(bool(cf) for cf in f.all_coeffs()[1:-1]):
            return

        n = f.degree()
        F = self._cyclotomic_decompose(n)

        if tc_f != 1:
            return F
        H = []

        for h in self._cyclotomic_decompose(2*n):
            if h not in F:
                H.append(h)

        return H

    def _cyclotomic_p(self, f, irreducible=False):
        """
        Efficiently test if ``f`` is a cyclotomic polynomial.

        Examples
        ========

        >>> _, x = ring('x', ZZ)

        >>> (x**16 + x**14 - x**10 + x**8 - x**6 + x**2 + 1).is_cyclotomic
        False

        >>> (x**16 + x**14 - x**10 - x**8 - x**6 + x**2 + 1).is_cyclotomic
        True

        """
        domain = self.domain

        if domain.is_RationalField:
            try:
                f = f.set_domain(domain.ring)
                return f.is_cyclotomic
            except CoercionFailedError:
                return False
        elif not domain.is_IntegerRing:
            return False

        x = self.gens[0]

        lc = f.LC
        tc = f[1]

        if lc != 1 or tc not in (1, -1):
            return False

        if not irreducible:
            coeff, factors = f.factor_list()

            if coeff != 1 or factors != [(f, 1)]:
                return False

        n = f.degree()
        g, h = self.zero, self.zero

        for j, i in enumerate(range(n, -1, -2)):
            g += f[(i,)]*x**j

        for j, i in enumerate(range(n - 1, -1, -2)):
            h += f[(i,)]*x**j

        g = g**2
        h = h**2

        F = g - h*self.from_terms([((1,), domain.one)])

        if F.LC < 0:
            F = -F

        if F == f:
            return True

        g = f.compose(x, -x)

        if g.LC < 0:
            g = -g

        if F == g and g.is_cyclotomic:
            return True

        G = F.sqf_part()

        if G**2 == F and G.is_cyclotomic:
            return True

        return False

    def _zz_cyclotomic_poly(self, n):
        """Efficiently generate n-th cyclotomic polynomial."""
        x = self.gens[0]

        h = x - 1

        for p, k in factorint(n).items():
            h = h.compose(x, x**p) // h
            h = h.compose(x, x**(p**(k - 1)))

        return h

    def _univar_zz_diophantine(self, F, m, p):
        """Wang/EEZ: Solve univariate Diophantine equations."""
        domain = self.domain
        m = self.from_terms([((m,), domain.one)])

        if len(F) == 2:
            p_domain = domain.finite_field(p)
            p_ring = self.clone(domain=p_domain)
            f, g = map(operator.methodcaller('set_domain', p_domain), F)

            s, t, _ = p_ring.gcdex(g, f)

            s *= m
            t *= m

            q, s = divmod(s, f)
            s = s.set_domain(domain)

            t += q*g
            t = t.set_domain(domain)

            result = [s, t]
        else:
            G = [F[-1]]

            for f in reversed(F[1:-1]):
                G.insert(0, f*G[0])

            S, T = [], [self.one]

            for f, g in zip(F, G):
                t, s = self._zz_diophantine([g, f], T[-1], [], 0, p)
                T.append(t)
                S.append(s)

            result, S = [], S + [T[-1]]
            p_domain = domain.finite_field(p)

            for s, f in zip(S, F):
                s *= m
                s, f = map(operator.methodcaller('set_domain', p_domain),
                           (s, f))
                s = (s % f).set_domain(domain)

                result.append(s)

        return result

    def _zz_diophantine(self, F, c, A, d, p):
        """Wang/EEZ: Solve multivariate Diophantine equations."""
        domain = self.domain

        if not A:
            S = [self.zero for _ in F]
            n = c.degree()

            for i, coeff in enumerate(c.all_coeffs()):
                if not coeff:
                    continue

                T = self._univar_zz_diophantine(F, i, p)

                for j, (s, t) in enumerate(zip(S, T)):
                    t *= coeff
                    S[j] = (s + t).trunc_ground(p)
        else:
            n = len(A)
            e = math.prod(F)

            a, A = A[-1], A[:-1]
            B, G = [], []

            for f in F:
                B.append(e // f)
                G.append(f.eval(x=n, a=a))

            C = c.eval(x=n, a=a)

            S = self.drop(n)._zz_diophantine(G, C, A, d, p)
            S = [s.set_ring(self) for s in S]

            for s, b in zip(S, B):
                c -= s*b

            c = c.trunc_ground(p)

            m = self.gens[n] - a
            M = self.one

            for k in range(d):
                k = domain(k)
                if not c:
                    break

                M *= m
                C = c.diff(x=n, m=int(k + 1)).eval(x=n, a=a)

                if C:
                    C = C.quo_ground(domain.factorial(k + 1))
                    T = C.ring._zz_diophantine(G, C, A, d, p)

                    for i, t in enumerate(T):
                        T[i] = t.set_ring(self)*M

                    for i, (s, t) in enumerate(zip(S, T)):
                        S[i] = s + t

                    for t, b in zip(T, B):
                        c -= t*b

                    c = c.trunc_ground(p)

            S = [s.trunc_ground(p) for s in S]

        return S

    def _zz_wang(self, f, mod=None, seed=None):
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

        ct, T = f.eject(*self.gens[1:]).LC.factor_list()

        domain = self.domain
        uring = self.drop(*self.gens[1:])
        b = self._zz_mignotte_bound(f)
        p = domain(nextprime(b))

        if mod is None:
            if self.ngens == 2:
                mod = 2
            else:
                mod = 1

        history, configs, A, r = set(), [], [domain.zero]*(self.ngens - 1), None

        try:
            cs, s, E = self._zz_wang_test_points(f, T, ct, A)

            _, H = uring._zz_factor_sqf(s)

            r = len(H)

            if r == 1:
                return [f]

            configs = [(s, cs, E, H, A)]
        except EvaluationFailedError:
            pass

        eez_num_configs = query('EEZ_NUMBER_OF_CONFIGS')
        eez_num_tries = query('EEZ_NUMBER_OF_TRIES')
        eez_mod_step = query('EEZ_MODULUS_STEP')

        while len(configs) < eez_num_configs:
            for _ in range(eez_num_tries):
                A = [domain(randint(-mod, mod)) for _ in range(self.ngens - 1)]

                if tuple(A) in history:
                    continue
                history.add(tuple(A))

                try:
                    cs, s, E = self._zz_wang_test_points(f, T, ct, A)
                except EvaluationFailedError:
                    continue

                _, H = uring._zz_factor_sqf(s)

                rr = len(H)

                if r is not None:
                    if rr != r:
                        if rr >= r:
                            continue
                        configs, r = [], rr
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
            _s_norm = s.max_norm()

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
            f, H, LC = self._zz_wang_lead_coeffs(f, T, cs, E, H, A)
            factors = self._zz_wang_hensel_lifting(f, H, LC, A, p)
        except ExtraneousFactorsError as exc:
            if query('EEZ_RESTART_IF_NEEDED'):
                return self._zz_wang(orig_f, mod + 1)
            raise ExtraneousFactorsError('we need to restart algorithm '
                                         'with better parameters') from exc

        result = []

        for f in factors:
            _, f = f.primitive()

            result.append(f)

        return result

    def _zz_wang_test_points(self, f, T, ct, A):
        """Wang/EEZ: Test evaluation points for suitability."""
        if not f.eject(*self.gens[1:]).LC(*A):
            raise EvaluationFailedError('no luck')

        g = f.eject(0)(*A)

        if not g.is_squarefree:
            raise EvaluationFailedError('no luck')

        c, h = g.primitive()

        E = [t(*A) for t, _ in T]
        D = self._zz_wang_non_divisors(E, c, ct)

        if D is not None:
            return c, h, E
        raise EvaluationFailedError('no luck')

    def _zz_wang_non_divisors(self, E, cs, ct):
        """Wang/EEZ: Compute a set of valid divisors."""
        domain = self.domain
        result = [cs*ct]

        for q in E:
            q = abs(q)

            for r in reversed(result):
                while r != 1:
                    r = domain.gcd(r, q)
                    q = q // r

                if q == 1:
                    return

            result.append(q)

        return result[1:]

    def _zz_wang_lead_coeffs(self, f, T, cs, E, H, A):
        """Wang/EEZ: Compute correct leading coefficients."""
        domain = self.domain
        c_ring = self.drop(0)
        C, J = [], [0]*len(E)

        for h in H:
            c = c_ring.one
            d = h.LC*cs

            for i in reversed(range(len(E))):
                k, e, t = 0, E[i], T[i][0]

                while not d % e:
                    d, k = d//e, k + 1

                if k != 0:
                    c *= t**k
                    J[i] = 1

            C.append(c)

        if any(not j for j in J):  # pragma: no cover
            raise ExtraneousFactorsError

        CC, HH = [], []

        for c, h in zip(C, H):
            d = c(*A)
            lc = h.LC

            if cs == 1:
                cc = lc//d
            else:
                g = domain.gcd(lc, d)
                d //= g
                cc = lc//g
                h *= d
                cs //= d

            c *= cc

            CC.append(c)
            HH.append(h)

        if cs == 1:
            return f, HH, CC

        CCC, HHH = [], []

        for c, h in zip(CC, HH):
            CCC.append(c*cs)
            HHH.append(h*cs)

        f *= cs**(len(H) - 1)

        return f, HHH, CCC

    def _zz_wang_hensel_lifting(self, f, H, LC, A, p):
        """Wang/EEZ: Parallel Hensel lifting algorithm."""
        domain = self.domain
        S, n = [f], len(A)

        H = list(H)

        for i, a in enumerate(reversed(A[1:])):
            s = S[0].eval(x=n - i, a=a)
            S.insert(0, s.trunc_ground(p))

        d = max(f.degree(x) for x in self.gens[1:])

        for j, s, a in zip(range(2, n + 2), S, A):
            G, w = list(H), j - 1
            s_ring = s.ring

            I, J = A[:w - 1], A[w:]

            for i, (h, lc) in enumerate(zip(H, LC)):
                if J:
                    lc = lc.eject(*lc.ring.gens[:-len(J)])(*J)
                lc = lc.trunc_ground(p)
                h, lc = map(operator.methodcaller('set_ring', s_ring), (h, lc))
                lt = h.eject(*s_ring.gens[1:]).leading_term().inject()
                H[i] = lc*s_ring.gens[0]**lt.degree() + h - lt

            m = s_ring.gens[-1] - a
            M = s_ring.one

            c = math.prod(H)
            c = s - c

            dj = s.degree(x=w)

            for k in range(dj):
                k = domain(k)
                if not c:
                    break

                M *= m
                C = c.diff(x=w, m=int(k + 1)).eval(x=w, a=a)

                if C:
                    C = C.quo_ground(domain.factorial(k + 1))
                    T = C.ring._zz_diophantine(G, C, I, d, p)

                    for i, (h, t) in enumerate(zip(H, T)):
                        t = t.set_ring(s_ring)
                        h += t*M
                        H[i] = h.trunc_ground(p)

                    h = math.prod(H)
                    h = s - h
                    c = h.trunc_ground(p)

        if math.prod(H) != f:
            raise ExtraneousFactorsError
        return H

    def _gf_Qmatrix(self, f):
        """
        Calculate Berlekamp's ``Q`` matrix.

        Examples
        ========

        >>> R, x = ring('x', FF(5))

        >>> f = 3*x**2 + 2*x + 4
        >>> R._gf_Qmatrix(f)
        [[1, 0], [3, 4]]

        >>> f = x**4 + 1
        >>> R._gf_Qmatrix(f)
        [[1, 0, 0, 0], [0, 4, 0, 0], [0, 0, 1, 0], [0, 0, 0, 4]]

        References
        ==========

        * :cite:`Geddes1992algorithms`, algorithm 8.5

        """
        domain = self.domain
        n, q = f.degree(), domain.order

        r = [domain.one] + [domain.zero]*(n - 1)
        Q = [r.copy()] + [[]]*(n - 1)
        f = f.all_coeffs()

        for i in range(1, (n - 1)*q + 1):
            c, r[1:], r[0] = r[-1], r[:-1], domain.zero
            for j in range(n):
                r[j] -= c*f[j]

            if not i % q:
                Q[i//q] = r.copy()

        return Q

    def _gf_berlekamp(self, f):
        """
        Factor a square-free polynomial over finite fields of small order.

        Examples
        ========

        >>> R, x = ring('x', FF(5))

        >>> R._gf_berlekamp(x**4 + 1)
        [x**2 + 2, x**2 + 3]

        References
        ==========

        * :cite:`Geddes1992algorithms`, algorithm 8.4
        * :cite:`Knuth1985seminumerical`, section 4.6.2

        """
        from .solvers import RawMatrix

        assert self.is_univariate

        domain = self.domain

        Q = self._gf_Qmatrix(f)
        Q = RawMatrix(Q) - RawMatrix.eye(len(Q))
        V = Q.T.nullspace()

        for i, v in enumerate(V):
            V[i] = self.from_list(v)

        factors = [f]

        for v in V[1:]:
            for f in list(factors):
                for s in range(domain.order):
                    h = v - domain(s)
                    g = self.gcd(f, h)

                    if g not in (1, f):
                        factors.remove(f)

                        f //= g
                        factors.extend([f, g])

                    if len(factors) == len(V):
                        return _sort_factors(factors, multiple=False)

        return _sort_factors(factors, multiple=False)

    def _gf_ddf_zassenhaus(self, f):
        """
        Cantor-Zassenhaus: Deterministic Distinct Degree Factorization.

        Given a monic square-free polynomial ``f`` in ``GF(q)[x]``, computes
        partial distinct degree factorization ``f_1 ... f_d`` of ``f`` where
        ``deg(f_i) != deg(f_j)`` for ``i != j``. The result is returned as a
        list of pairs ``(f_i, e_i)`` where ``deg(f_i) > 0`` and ``e_i > 0``
        is an argument to the equal degree factorization routine.

        Examples
        ========

        >>> R, x = ring('x', FF(11))
        >>> R._gf_ddf_zassenhaus(x**15 - 1)
        [(x**5 + 10, 1), (x**10 + x**5 + 1, 2)]

        To obtain factorization into irreducibles, use equal degree factorization
        procedure (EDF) with each of the factors.

        References
        ==========

        * :cite:`Gathen1999modern`, algorithm 14.3
        * :cite:`Geddes1992algorithms`, algorithm 8.8

        See Also
        ========

        _gf_edf_zassenhaus

        """
        domain = self.domain

        factors, q = [], domain.order
        g, x = [self.gens[0]]*2

        for i in range(1, f.degree()//2 + 1):
            g = pow(g, q, f)
            h = self.gcd(f, g - x)

            if h != 1:
                factors.append((h, i))

                f //= h
                g %= f

        if f != 1:
            factors += [(f, f.degree())]

        return factors

    def _gf_edf_zassenhaus(self, f, n):
        """
        Cantor-Zassenhaus: Probabilistic Equal Degree Factorization.

        Given a monic square-free polynomial ``f`` in ``GF(q)[x]`` and
        an integer ``n``, such that ``n`` divides ``deg(f)``, returns all
        irreducible factors ``f_1,...,f_d`` of ``f``, each of degree ``n``.
        EDF procedure gives complete factorization over Galois fields.

        Examples
        ========

        >>> R, x = ring('x', FF(5))
        >>> R._gf_edf_zassenhaus(x**3 + x**2 + x + 1, 1)
        [x + 1, x + 2, x + 3]

        References
        ==========

        * :cite:`Geddes1992algorithms`, algorithm 8.9

        See Also
        ========

        _gf_ddf_zassenhaus

        """
        factors = [f]
        d = f.degree()

        if d <= n:
            return factors

        domain = self.domain
        p, q = domain.characteristic, domain.order
        N = d // n

        while len(factors) < N:
            r = self._gf_random(2*n - 1)

            if p == 2:
                h = r

                for _ in range(1, n):
                    h += pow(r, q, f)
            else:
                h = pow(r, (q**n - 1)//2, f)
                h -= 1

            g = self.gcd(f, h)

            if g not in (1, f):
                factors = (self._gf_edf_zassenhaus(g, n) +
                           self._gf_edf_zassenhaus(f // g, n))

        return _sort_factors(factors, multiple=False)

    def _gf_zassenhaus(self, f):
        """
        Factor a square-free polynomial over finite fields of medium order.

        Examples
        ========

        >>> R, x = ring('x', FF(5))
        >>> R._gf_zassenhaus(x**2 + 4*x + 3)
        [x + 1, x + 3]

        """
        assert self.is_univariate

        factors = []

        for factor, n in self._gf_ddf_zassenhaus(f):
            factors += self._gf_edf_zassenhaus(factor, n)

        return _sort_factors(factors, multiple=False)

    def _gf_ddf_shoup(self, f):
        """
        Kaltofen-Shoup: Deterministic Distinct Degree Factorization.

        Given a monic square-free polynomial ``f`` in ``GF(q)[x]``, computes
        partial distinct degree factorization ``f_1 ... f_d`` of ``f`` where
        ``deg(f_i) != deg(f_j)`` for ``i != j``. The result is returned as a
        list of pairs ``(f_i, e_i)`` where ``deg(f_i) > 0`` and ``e_i > 0``
        is an argument to the equal degree factorization routine.

        Notes
        =====

        This algorithm is an improved version of Zassenhaus algorithm for
        large ``deg(f)`` and order ``q`` (especially for ``deg(f) ~ lg(q)``).

        Examples
        ========

        >>> R, x = ring('x', FF(3))
        >>> R._gf_ddf_shoup(x**6 - x**5 + x**4 + x**3 - x)
        [(x**2 + x, 1), (x**4 + x**3 + x + 2, 2)]

        References
        ==========

        * :cite:`Kaltofen1998subquadratic`, algorithm D
        * :cite:`Shoup1995factor`
        * :cite:`Gathen1992frobenious`

        See Also
        ========

        _gf_edf_shoup

        """
        domain = self.domain

        n, q = f.degree(), domain.order
        k = math.isqrt(n//2 - 1) + 1 if n > 1 else 0
        x = self.gens[0]

        h = pow(x, q, f)

        # U[i] = x**(q**i)
        U = [x, h] + [self.zero]*(k - 1)

        for i in range(2, k + 1):
            U[i] = U[i - 1].compose(x, h) % f

        h, U = U[k], U[:k]
        # V[i] = x**(q**(k*(i+1)))
        V = [h] + [self.zero]*(k - 1)

        for i in range(1, k):
            V[i] = V[i - 1].compose(x, h) % f

        factors = []

        for i, v in enumerate(V):
            h, j = self.one, k - 1

            for u in U:
                g = v - u
                h *= g
                h %= f

            g = self.gcd(f, h)
            f //= g

            for u in reversed(U):
                h = v - u
                F = self.gcd(g, h)

                if F != 1:
                    factors.append((F, k*(i + 1) - j))

                g //= F
                j -= 1

        if f != 1:
            factors.append((f, f.degree()))

        return factors

    def _gf_trace_map(self, a, b, c, n, f):
        """
        Compute polynomial trace map in ``GF(q)[x]/(f)``.

        Given a polynomial ``f`` in ``GF(q)[x]``, polynomials ``a``, ``b``,
        ``c`` in the quotient ring ``GF(q)[x]/(f)`` such that ``b = c**t
        (mod f)`` for some positive power ``t`` of ``q``, and a positive
        integer ``n``, returns a mapping::

           a -> a**t**n, a + a**t + a**t**2 + ... + a**t**n (mod f)

        In factorization context, ``b = x**q mod f`` and ``c = x mod f``.
        This way we can efficiently compute trace polynomials in equal
        degree factorization routine, much faster than with other methods,
        like iterated Frobenius algorithm, for large degrees.

        Examples
        ========

        >>> R, x = ring('x', FF(5))
        >>> a = x + 2
        >>> b = 4*x + 4
        >>> c = x + 1
        >>> f = 3*x**2 + 2*x + 4
        >>> R._gf_trace_map(a, b, c, 4, f)
        (x + 3, x + 3)

        References
        ==========

        * :cite:`Gathen1992ComputingFM`, algorithm 5.2

        """
        u = a.compose(0, b) % f
        v = b

        if n & 1:
            U = a + u
            V = b
        else:
            U = a
            V = c

        n >>= 1

        while n:
            u += u.compose(0, v) % f
            v = v.compose(0, v) % f

            if n & 1:
                U += u.compose(0, V) % f
                V = v.compose(0, V) % f

            n >>= 1

        return a.compose(0, V) % f, U

    def _gf_edf_shoup(self, f, n):
        """
        Gathen-Shoup: Probabilistic Equal Degree Factorization.

        Given a monic square-free polynomial ``f`` in ``GF(q)[x]`` and
        an integer ``n``, such that ``n`` divides ``deg(f)``, returns all
        irreducible factors ``f_1,...,f_d`` of ``f``, each of degree ``n``.
        EDF procedure gives complete factorization over Galois fields.

        Notes
        =====

        This algorithm is an improved version of Zassenhaus algorithm for
        large ``deg(f)`` and order ``q`` (especially for ``deg(f) ~ lg(q)``).

        Examples
        ========

        >>> R, x = ring('x', FF(2917))
        >>> R._gf_edf_shoup(x**2 + 2837*x + 2277, 1)
        [x + 852, x + 1985]

        References
        ==========

        * :cite:`Shoup1991ffactor`
        * :cite:`Gathen1992ComputingFM`, algorithm 3.6

        See Also
        ========

        _gf_ddf_shoup

        """
        domain = self.domain
        q, p = domain.order, domain.characteristic
        N = f.degree()

        if not N:
            return []
        if N <= n:
            return [f]

        factors, x = [f], self.gens[0]

        r = self._gf_random(N - 1)

        h = pow(x, q, f)
        H = self._gf_trace_map(r, h, x, n - 1, f)[1]

        if p == 2:
            h1 = self.gcd(f, H)
            h2 = f // h1

            factors = self._gf_edf_shoup(h1, n) + self._gf_edf_shoup(h2, n)
        else:
            h = pow(H, (q - 1)//2, f)

            h1 = self.gcd(f, h)
            h2 = self.gcd(f, h - 1)
            h3 = f // (h1 * h2)

            factors = (self._gf_edf_shoup(h1, n) + self._gf_edf_shoup(h2, n) +
                       self._gf_edf_shoup(h3, n))

        return _sort_factors(factors, multiple=False)

    def _gf_shoup(self, f):
        """
        Factor a square-free polynomial over finite fields of large order.

        Examples
        ========

        >>> R, x = ring('x', FF(5))
        >>> R._gf_shoup(x**2 + 4*x + 3)
        [x + 1, x + 3]

        """
        assert self.is_univariate

        factors = []

        for factor, n in self._gf_ddf_shoup(f):
            factors += self._gf_edf_shoup(factor, n)

        return _sort_factors(factors, multiple=False)

    def _gf_factor_sqf(self, f):
        """
        Factor a square-free polynomial ``f`` in ``GF(q)[x]``.

        Returns its complete factorization into irreducibles::

                     f_1(x) f_2(x) ... f_d(x)

        where each ``f_i`` is a monic polynomial and ``gcd(f_i, f_j) == 1``,
        for ``i != j``.  The result is given as a list of factors of ``f``.

        Square-free factors of ``f`` can be factored into irreducibles over
        finite fields using three very different methods:

        Berlekamp
            efficient for very small values of order ``q`` (usually ``q < 25``)
        Cantor-Zassenhaus
            efficient on average input and with "typical" ``q``
        Shoup-Kaltofen-Gathen
            efficient with very large inputs and order

        If you want to use a specific factorization method - set
        ``GF_FACTOR_METHOD`` configuration option with one of ``"berlekamp"``,
        ``"zassenhaus"`` or ``"shoup"`` values.

        Examples
        ========

        >>> R, x = ring('x', FF(5))
        >>> f = x**2 + 4*x + 3
        >>> R._gf_factor_sqf(f)
        [x + 1, x + 3]

        References
        ==========

        * :cite:`Gathen1999modern`, chapter 14

        """
        _factor_methods = {
            'berlekamp': self._gf_berlekamp,  # ``p`` : small
            'zassenhaus': self._gf_zassenhaus,  # ``p`` : medium
            'shoup': self._gf_shoup,      # ``p`` : large
        }
        method = query('GF_FACTOR_METHOD')

        return _factor_methods[method](f)
