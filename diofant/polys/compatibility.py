"""Compatibility interface between dense and sparse polys."""

from .densebasic import dmp_degree_in
from .rootisolation import (dup_count_complex_roots, dup_isolate_all_roots,
                            dup_isolate_all_roots_sqf,
                            dup_isolate_complex_roots_sqf, dup_real_imag,
                            dup_transform)


__all__ = 'IPolys',


class IPolys:
    """Compatibility class between dense and sparse polynomials."""

    symbols = None
    ngens = None
    domain = None
    order = None
    gens = None

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

    def dup_count_complex_roots(self, f, inf=None, sup=None, exclude=None):
        return dup_count_complex_roots(f.to_dense(), self.domain, inf=inf, sup=sup, exclude=exclude)

    def dup_isolate_complex_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_complex_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_all_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_all_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_all_roots(self, f, eps=None, inf=None, sup=None):
        return dup_isolate_all_roots(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup)
