"""Compatibility interface between dense and sparse polys."""

from .rootisolation import (dup_count_complex_roots,
                            dup_isolate_complex_roots_sqf)


__all__ = 'IPolys',


class IPolys:
    """Compatibility class between dense and sparse polynomials."""

    symbols = None
    ngens = None
    domain = None
    order = None
    gens = None

    def dup_count_complex_roots(self, f, inf=None, sup=None, exclude=None):
        return dup_count_complex_roots(f.to_dense(), self.domain, inf=inf, sup=sup, exclude=exclude)

    def dup_isolate_complex_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_complex_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)
