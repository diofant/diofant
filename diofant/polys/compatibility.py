"""Compatibility interface between dense and sparse polys."""

from .densearith import dmp_mul, dmp_neg, dmp_sub, dup_lshift, dup_rshift
from .densebasic import dmp_degree_in
from .rootisolation import (dup_count_complex_roots, dup_isolate_all_roots,
                            dup_isolate_all_roots_sqf,
                            dup_isolate_complex_roots_sqf,
                            dup_isolate_real_roots,
                            dup_isolate_real_roots_pair,
                            dup_isolate_real_roots_sqf, dup_real_imag,
                            dup_transform)


__all__ = 'IPolys',


class IPolys:
    """Compatibility class between dense and sparse polynomials."""

    symbols = None
    ngens = None
    domain = None
    order = None
    gens = None

    def dup_lshift(self, f, n):
        return self.from_list(dup_lshift(f.to_dense(), n, self.domain))

    def dup_rshift(self, f, n):
        return self.from_list(dup_rshift(f.to_dense(), n, self.domain))

    def dmp_neg(self, f):
        return self.from_list(dmp_neg(f.to_dense(), self.ngens-1, self.domain))

    def dmp_sub(self, f, g):
        return self.from_list(dmp_sub(f.to_dense(), g.to_dense(), self.ngens-1, self.domain))

    def dmp_mul(self, f, g):
        return self.from_list(dmp_mul(f.to_dense(), g.to_dense(), self.ngens-1, self.domain))

    def dmp_degree_in(self, f, j):
        return dmp_degree_in(f.to_dense(), j, self.ngens-1)

    def dup_real_imag(self, f):
        ring = self
        p, q = dup_real_imag(f.drop(1).to_dense(), ring.domain)
        if ring.domain.is_ComplexAlgebraicField and not ring.domain.is_RealAlgebraicField:
            ring = ring.to_ground()
        return ring.from_list(p), ring.from_list(q)

    def dup_transform(self, f, p, q):
        return self.from_list(dup_transform(f.to_dense(), p.to_dense(), q.to_dense(), self.domain))

    def dup_isolate_real_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_real_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_real_roots(self, f, eps=None, inf=None, sup=None):
        return dup_isolate_real_roots(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup)

    def dup_isolate_real_roots_pair(self, f, g, eps=None, inf=None, sup=None, strict=False, basis=False):
        return dup_isolate_real_roots_pair(*(map(lambda _: _.to_dense(), [f, g])), self.domain, eps=eps, inf=inf, sup=sup, strict=strict, basis=basis)

    def dup_count_complex_roots(self, f, inf=None, sup=None, exclude=None):
        return dup_count_complex_roots(f.to_dense(), self.domain, inf=inf, sup=sup, exclude=exclude)

    def dup_isolate_complex_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_complex_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_all_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_all_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_all_roots(self, f, eps=None, inf=None, sup=None):
        return dup_isolate_all_roots(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup)
