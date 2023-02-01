"""Square-free decomposition algorithms and related tools."""

from .polyerrors import DomainError


class _SQF:
    """Mixin class to compute square-free decomposition of polynomials."""

    def sqf_list(self, f):
        """
        Return square-free decomposition of a polynomial in ``K[X]``.

        Examples
        ========

        >>> R, x, y = ring('x y', ZZ)

        >>> R.sqf_list(x**5 + 2*x**4*y + x**3*y**2)
        (1, [(x + y, 2), (x, 3)])

        >>> R, x, y = ring('x y', FF(5))

        >>> f = x**5*y**5 + 1

        Note that e.g.:

        >>> f.diff(x)
        0

        This phenomenon doesn't happen in characteristic zero. However we can
        still compute square-free decomposition of ``f``:

        >>> R.sqf_list(f)
        (1, [(x*y + 1, 5)])

        """
        domain = self.domain

        if domain.is_Field:
            coeff, f = f.LC, f.monic()
        else:
            coeff, f = f.primitive()

        if domain.is_FiniteField:
            return coeff, self._gf_sqf_list(f)
        return coeff, self._rr_yun0_sqf_list(f)

    def _gf_sqf_list(self, f):
        """
        Compute square-free decomposition of the monic ``f`` in ``GF(q)[X]``.

        Notes
        =====

        Uses a modified version of Musser's algorithm for square-free
        decomposition of univariate polynomials over finite fields.

        References
        ==========

        * :cite:`Geddes1992algorithms`, algorithm 8.3
        * :cite:`Musser1971factor`, p.54, algorithm V (modified)

        """
        domain = self.domain

        n, factors, p = 1, [], int(domain.characteristic)
        m = int(domain.order // p)

        while not f.is_ground:
            df = [f.diff(x) for x in self.gens]

            if any(_ for _ in df):
                g = f
                for q in df:
                    g = self.gcd(g, q)
                h, f, i = f // g, g, 1

                while h != 1:
                    g = self.gcd(f, h)
                    h //= g

                    if not h.is_ground:
                        factors.append((h, i*n))

                    f //= g
                    h = g
                    i += 1

            n *= p

            g = self.zero
            for monom, coeff in f.items():
                g[(_ // p for _ in monom)] = coeff**m
            f = g

        return factors

    def _rr_yun0_sqf_list(self, f):
        """Compute square-free decomposition of ``f`` in zero-characteristic ring ``K[X]``.

        References
        ==========

        * :cite:`Geddes1992algorithms`, algorithm 8.2
        * :cite:`LeeM2013factor`, algorithm 3.1

        """
        if f.is_ground:
            return []

        result, count = [], 1
        qs = [f.diff(x) for x in self.gens]

        g = f
        for q in qs:
            g = self.gcd(g, q)

        while f != 1:
            qs = [q // g for q in qs]
            f //= g
            qs = [q - f.diff(x) for x, q in zip(self.gens, qs)]

            g = f
            for q in qs:
                g = self.gcd(g, q)
            if g != 1:
                result.append((g, count))

            count += 1

        return result

    def is_squarefree(self, f):
        """
        Return ``True`` if ``f`` is a square-free polynomial in ``K[X]``.

        Examples
        ========

        >>> _, x, y = ring('x y', ZZ)

        >>> ((x + y)**2).is_squarefree
        False
        >>> (x**2 + y**2).is_squarefree
        True

        """
        if f.is_ground:
            return True
        g = f
        for x in self.gens:
            g = self.gcd(g, f.diff(x))
            if g.is_ground:
                return True
        return False

    def sqf_part(self, f):
        """
        Returns square-free part of a polynomial in ``K[X]``.

        Examples
        ========

        >>> R, x, y = ring('x y', ZZ)

        >>> R.sqf_part(x**3 + 2*x**2*y + x*y**2)
        x**2 + x*y

        """
        domain = self.domain

        if domain.is_FiniteField:
            g = self.one
            for f, _ in self.sqf_list(f)[1]:
                g *= f

            return g

        if not f:
            return f

        gcd = f
        for x in self.gens:
            gcd = self.gcd(gcd, f.diff(x))
        sqf = f // gcd

        if domain.is_Field:
            return sqf.monic()
        return sqf.primitive()[1]

    def sqf_norm(self, f):
        """
        Square-free norm of ``f`` in ``K[X]``, useful over algebraic domains.

        Returns ``s``, ``f``, ``r``, such that ``g(x) = f(x-sa)`` and ``r(x) = Norm(g(x))``
        is a square-free polynomial over K, where ``a`` is the algebraic extension of ``K``.

        Examples
        ========

        >>> _, x, y = ring('x y', QQ.algebraic_field(I))

        >>> (x*y + y**2).sqf_norm()
        (1, x*y - I*x + y**2 - 3*I*y - 2,
         x**2*y**2 + x**2 + 2*x*y**3 + 2*x*y + y**4 + 5*y**2 + 4)

        """
        domain = self.domain

        if not domain.is_AlgebraicField:
            raise DomainError(f'ground domain must be algebraic, got {domain}')

        new_ring = self.to_ground().inject(*domain.symbols, front=True)
        g = domain.mod.set_ring(new_ring)
        s = 0

        while True:
            h = f.inject(front=True)
            r = g.resultant(h)

            if r.is_squarefree:
                return s, f, r
            f = f.compose({x: x - domain.unit for x in self.gens})
            s += 1
