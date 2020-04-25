"""Euclidean algorithms, GCDs, LCMs and polynomial remainder sequences."""

from ..core import cacheit
from ..ntheory import nextprime
from ..ntheory.modular import crt
from .polyerrors import HomomorphismFailed
from .polyutils import dmp_compat_wrapper


def dup_gcdex(f, g, K):
    """Extended Euclidean algorithm in `F[x]`."""
    ring = K.poly_ring(*[f'_{i}' for i in range(1)])
    f, g = map(ring.from_dense, (f, g))
    return tuple(map(ring.to_dense, f.gcdex(g)))


@cacheit
@dmp_compat_wrapper
def dmp_resultant(ring, f, g):
    """Computes resultant of two polynomials in `K[X]`."""
    res = f.resultant(g)
    return ring.drop(0).to_dense(res)


@dmp_compat_wrapper
def dmp_inner_gcd(ring, f, g):
    """Computes polynomial GCD and cofactors of `f` and `g` in `K[X]`."""
    return tuple(map(ring.to_dense, f.cofactors(g)))


def dmp_gcd(f, g, u, K):
    """Computes polynomial GCD of `f` and `g` in `K[X]`."""
    return dmp_inner_gcd(f, g, u, K)[0]


@dmp_compat_wrapper
def dmp_primitive(ring, f):
    """Returns multivariate content and a primitive polynomial."""
    new_ring, f = map(lambda _: _.eject(*ring.gens[1:]), (ring, f))
    c, f = f.primitive()
    f = f.inject()
    return new_ring.domain.to_dense(c), ring.to_dense(f)


class _GCD:
    """Mixin class for computing gcd."""

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
            ring, f, g = map(lambda _: _.eject(*self.gens[1:]), (ring, f, g))
            h, f, g = map(lambda _: _.inject(), ring._rr_prs_gcd(f, g))
            c, h = h.LC, h.monic()
            return h, f.quo_ground(c), g.quo_ground(c)

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
