"""Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences."""

from ..core import cacheit
from ..ntheory import nextprime
from ..ntheory.modular import crt, symmetric_residue
from .polyconfig import query
from .polyerrors import DomainError, HeuristicGCDFailed, HomomorphismFailed


def dup_gcdex(f, g, K):
    """Extended Euclidean algorithm in `F[x]`."""
    ring = K.poly_ring('_0')
    f, g = map(ring.from_list, (f, g))
    return tuple(map(lambda _: _.to_dense(), f.gcdex(g)))


@cacheit
def dmp_resultant(f, g, u, K):
    """Computes resultant of two polynomials in `K[X]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    res = f.resultant(g)
    return res.to_dense()


def dmp_gcd(f, g, u, K):
    """Computes polynomial GCD of `f` and `g` in `K[X]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f, g = map(ring.from_list, (f, g))
    return ring.gcd(f, g).to_dense()


def dmp_primitive(f, u, K):
    """Returns multivariate content and a primitive polynomial."""
    ring = K.poly_ring(*[f'_{i}' for i in range(u + 1)])
    f = ring.from_list(f)
    new_ring, f = map(lambda _: _.eject(*ring.gens[1:]), (ring, f))
    c, f = f.primitive()
    f = f.inject()
    return c.to_dense(), f.to_dense()


class _GCD:
    """Mixin class for computing gcd."""

    def gcd(self, f, g):
        """Returns GCD of ``f`` and ``g``."""
        return self.cofactors(f, g)[0]

    def cofactors(self, f, g):
        """Returns GCD and cofactors of ``f`` and ``g``."""
        if f.is_zero and g.is_zero:
            zero = self.zero
            return zero, zero, zero
        elif f.is_zero:
            h, cff, cfg = self._gcd_zero(g)
            return h, cff, cfg
        elif g.is_zero:
            h, cfg, cff = self._gcd_zero(f)
            return h, cff, cfg

        J, (f, g) = f.deflate(g)
        h, cff, cfg = self._gcd(f, g)

        return h.inflate(J), cff.inflate(J), cfg.inflate(J)

    def _gcd_zero(self, f):
        one, zero = self.one, self.zero
        if self.domain.is_Field:
            return f.monic(), zero, self.ground_new(f.LC)
        else:
            if not self.is_normal(f):
                return -f, zero, -one
            else:
                return f, zero, one

    def _gcd(self, f, g):
        domain = self.domain

        if domain.is_RationalField:
            return self._gcd_QQ(f, g)
        elif domain.is_IntegerRing:
            return self._gcd_ZZ(f, g)
        elif domain.is_AlgebraicField:
            return self._gcd_AA(f, g)
        elif not domain.is_Exact:
            try:
                exact = domain.get_exact()
            except DomainError:
                return self.one, f, g

            f, g = map(lambda x: x.set_domain(exact), (f, g))
            ring = self.clone(domain=exact)

            return tuple(map(lambda x: x.set_domain(domain), ring.cofactors(f, g)))
        elif domain.is_Field:
            return self._ff_prs_gcd(f, g)
        else:
            return self._rr_prs_gcd(f, g)

    def _gcd_ZZ(self, f, g):
        from .modulargcd import modgcd

        if query('USE_HEU_GCD'):
            try:
                return self._zz_heu_gcd(f, g)
            except HeuristicGCDFailed:
                pass

        _gcd_zz_methods = {'modgcd': modgcd,
                           'prs': self._rr_prs_gcd}

        method = _gcd_zz_methods[query('FALLBACK_GCD_ZZ_METHOD')]
        return method(f, g)

    def _gcd_QQ(self, f, g):
        domain = self.domain

        cf, f = f.clear_denoms(convert=True)
        cg, g = g.clear_denoms(convert=True)

        ring = self.clone(domain=domain.ring)

        h, cff, cfg = map(lambda _: _.set_ring(self), ring._gcd_ZZ(f, g))

        c, h = h.LC, h.monic()

        cff = cff.mul_ground(domain.quo(c, cf))
        cfg = cfg.mul_ground(domain.quo(c, cg))

        return h, cff, cfg

    def _gcd_AA(self, f, g):
        from .modulargcd import func_field_modgcd

        _gcd_aa_methods = {'modgcd': func_field_modgcd,
                           'prs': self._ff_prs_gcd}

        method = _gcd_aa_methods[query('GCD_AA_METHOD')]
        return method(f, g)

    def _zz_heu_gcd(self, f, g):
        """
        Heuristic polynomial GCD in ``Z[X]``.

        Given univariate polynomials ``f`` and ``g`` in ``Z[X]``, returns
        their GCD and cofactors, i.e. polynomials ``h``, ``cff`` and ``cfg``
        such that::

              h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

        The algorithm is purely heuristic which means it may fail to compute
        the GCD. This will be signaled by raising an exception. In this case
        you will need to switch to another GCD method.

        The algorithm computes the polynomial GCD by evaluating polynomials
        ``f`` and ``g`` at certain points and computing (fast) integer GCD
        of those evaluations. The polynomial GCD is recovered from the integer
        image by interpolation. The evaluation proces reduces f and g variable
        by variable into a large integer. The final step is to verify if the
        interpolated polynomial is the correct GCD. This gives cofactors of
        the input polynomials as a side effect.

        References
        ==========

        * :cite:`Liao1995heuristic`

        """
        assert self == f.ring == g.ring and self.domain.is_IntegerRing

        ring = self
        x0 = ring.gens[0]
        domain = ring.domain

        gcd, f, g = f.extract_ground(g)

        f_norm = f.max_norm()
        g_norm = g.max_norm()

        B = domain(2*min(f_norm, g_norm) + 29)

        x = max(min(B, 99*domain.sqrt(B)),
                2*min(f_norm // abs(f.LC),
                      g_norm // abs(g.LC)) + 4)

        cofactors = domain.cofactors if ring.is_univariate else ring.drop(0)._zz_heu_gcd

        for i in range(query('HEU_GCD_MAX')):
            ff = f.eval(x0, x)
            gg = g.eval(x0, x)

            if ff and gg:
                h, cff, cfg = cofactors(ff, gg)
                h = ring._gcd_interpolate(h, x)
                h = h.primitive()[1]

                cff_, r = divmod(f, h)

                if not r:
                    cfg_, r = divmod(g, h)

                    if not r:
                        h *= gcd
                        return h, cff_, cfg_

                cff = ring._gcd_interpolate(cff, x)

                h, r = divmod(f, cff)

                if not r:
                    cfg_, r = divmod(g, h)

                    if not r:
                        h *= gcd
                        return h, cff, cfg_

                cfg = ring._gcd_interpolate(cfg, x)

                h, r = divmod(g, cfg)

                if not r:
                    cff_, r = divmod(f, h)

                    if not r:
                        h *= gcd
                        return h, cff_, cfg

            x = 73794*x * domain.sqrt(domain.sqrt(x)) // 27011

        raise HeuristicGCDFailed('no luck')

    def _gcd_interpolate(self, h, x):
        """Interpolate polynomial GCD from integer GCD."""
        f, i = self.zero, 0
        X = self.gens[0]

        while h:
            g = h % x

            if self.is_univariate:
                g = symmetric_residue(g, x)

            h = (h - g) // x

            f += X**i*g
            i += 1

        if not self.is_normal(f):
            f = -f

        return f

    def _rr_prs_gcd(self, f, g):
        """Computes polynomial GCD using subresultants over a ring."""
        ring = self

        if self.is_multivariate:
            ring, f, g = map(lambda _: _.eject(*self.gens[1:]), (ring, f, g))
            return tuple(map(lambda _: _.inject(), ring._rr_prs_gcd(f, g)))

        domain = ring.domain

        fc, ff = f.primitive()
        gc, fg = g.primitive()

        h = ff.subresultants(fg)[-1]
        _, h = h.primitive()

        c = domain.gcd(fc, gc)
        h *= c

        return h, f // h, g // h

    def _ff_prs_gcd(self, f, g):
        """Computes polynomial GCD using subresultants over a field."""
        ring = self

        if ring.is_multivariate:
            ring, F, G = map(lambda _: _.eject(*self.gens[1:]), (ring, f, g))

            fc, F = F.primitive()
            gc, G = G.primitive()

            F, G = map(lambda _: _.inject(), (F, G))

            h = F.subresultants(G)[-1]
            c, _, _ = ring.domain._ff_prs_gcd(fc, gc)
            h = h.eject(*self.gens[1:])
            _, h = h.primitive()
            h = h.inject()
            h *= c
            h = h.monic()

            return h, f // h, g // h

        h = f.subresultants(g)[-1]
        h = h.monic()

        return h, f // h, g // h

    def _primitive_prs(self, f, g):
        """
        Subresultant PRS algorithm in `K[X]`.

        Computes the last non-zero scalar subresultant of `f` and `g`
        and subresultant polynomial remainder sequence (PRS).

        The first subdeterminant is set to 1 by convention to match
        the polynomial and the scalar subdeterminants.
        If 'deg(f) < deg(g)', the subresultants of '(g,f)' are computed.

        Examples
        ========

        >>> R, x = ring('x', ZZ)

        >>> (x**2 + 1).resultant(x**2 - 1, includePRS=True)
        (4, [x**2 + 1, x**2 - 1, -2])

        References
        ==========

        * :cite:`Brown1978prs`
        * :cite:`Geddes1992algorithms`, example 7.6

        """
        ring = self
        domain = ring.domain

        if ring.is_multivariate:
            ring, f, g = map(lambda _: _.eject(*ring.gens[1:]), (ring, f, g))
            res = ring._primitive_prs(f, g)
            return res[0], [_.inject() for _ in res[1]]

        n = f.degree()
        m = g.degree()

        if n < m:
            f, g = g, f
            n, m = m, n

        c, r = domain.zero, []

        if f.is_zero:
            return c, r

        r.append(f)

        if g.is_zero:
            return c, r

        r.append(g)
        d = n - m

        h = f.prem(g)
        h *= (-domain.one)**(d + 1)

        lc = g.LC
        c = -lc**d

        while not h.is_zero:
            k = h.degree()
            r.append(h)

            f, g, m, d = g, h, k, m - k

            h = f.prem(g)
            h = h.quo_ground(-lc*c**d)

            lc = g.LC
            c = domain.quo((-lc)**d, c**(d - 1))

        if r[-1].degree() > 0:
            c = domain.zero
        else:
            c = -c

        return c, r

    def _collins_resultant(self, f, g):
        """
        Collins's modular resultant algorithm in `ZZ[X]` or `QQ[X]`.

        References
        ==========

        * :cite:`Collins1971mod`, algorithm PRES

        """
        ring = self
        domain = ring.domain

        if f.is_zero or g.is_zero:
            return ring.drop(0).zero

        n = f.degree()
        m = g.degree()

        if domain.is_RationalField:
            cf, f = f.clear_denoms(convert=True)
            cg, g = g.clear_denoms(convert=True)

            ring = ring.clone(domain=domain.ring)
            r = ring._collins_resultant(f, g)
            r = r.set_domain(domain)

            c = cf**n * cg**m
            c = domain.convert(c, domain.ring)

            return r.quo_ground(c)

        assert domain.is_IntegerRing

        A = f.max_norm()
        B = g.max_norm()

        a, b = f.LC, g.LC

        B = domain(2)*domain.factorial(domain(n + m))*A**m*B**n
        new_ring = ring.drop(0)
        r, p, P = new_ring.zero, domain.one, domain.one

        while P <= B:
            while True:
                p = domain(nextprime(p))
                if (a % p) and (b % p):
                    break

            p_domain = domain.finite_field(p)
            F, G = map(lambda _: _.set_domain(p_domain), (f, g))

            try:
                R = ring.clone(domain=p_domain)._modular_resultant(F, G)
            except HomomorphismFailed:
                continue

            if P == domain.one:
                r = R
            else:
                def _crt(r, R):
                    return domain(crt([P, p], map(domain.convert, [r, R]),
                                  check=False, symmetric=True)[0])

                if new_ring.is_PolynomialRing:
                    r_new = new_ring.zero

                    for monom in set(r.keys()) | set(R.keys()):
                        r_new[monom] = _crt(r.get(monom, 0), R.get(monom, 0))
                    r = r_new
                else:
                    r = _crt(r, R)

            P *= p

        return r

    def _modular_resultant(self, f, g):
        """
        Compute resultant of `f` and `g` in `GF(p)[X]`.

        References
        ==========

        * :cite:`Collins1971mod`, algorithm CPRES

        """
        ring = self
        domain = ring.domain

        assert domain.is_FiniteField

        if ring.is_univariate:
            return ring._primitive_prs(f, g)[0]

        n = f.degree()
        m = g.degree()

        N = f.degree(1)
        M = g.degree(1)

        B = n*M + m*N

        new_ring = ring.drop(0)
        r = new_ring.zero
        D = ring.eject(1).domain.one
        domain_elts = iter(range(domain.order))

        while D.degree() <= B:
            while True:
                try:
                    a = next(domain_elts)
                except StopIteration:
                    raise HomomorphismFailed('no luck')

                F = f.eval(x=1, a=a)

                if F.degree() == n:
                    G = g.eval(x=1, a=a)

                    if G.degree() == m:
                        break

            R = ring.drop(1)._modular_resultant(F, G)
            e = r.eval(x=0, a=a)

            if new_ring.is_univariate:
                R = new_ring.ground_new(R)
                e = new_ring.ground_new(e)
            else:
                R = R.set_ring(new_ring)
                e = e.set_ring(new_ring)

            d = D * D.eval(x=0, a=a)**-1
            d = d.set_ring(new_ring)
            r += d*(R - e)
            D *= D.ring.gens[0] - a

        return r
