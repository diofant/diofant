"""Compatibility interface between dense and sparse polys."""

from .densearith import (dmp_abs, dmp_add, dmp_add_mul, dmp_add_term,
                         dmp_exquo_ground, dmp_mul, dmp_mul_ground,
                         dmp_mul_term, dmp_neg, dmp_quo_ground, dmp_sub,
                         dup_lshift, dup_rshift)
from .densebasic import (dmp_degree_in, dmp_degree_list, dmp_ground_TC, dmp_LC,
                         dmp_TC)
from .densetools import (dmp_compose, dup_decompose, dup_real_imag,
                         dup_transform)
from .factortools import (dmp_zz_wang_hensel_lifting, dmp_zz_wang_lead_coeffs,
                          dmp_zz_wang_non_divisors)
from .rootisolation import (dup_count_complex_roots, dup_count_real_roots,
                            dup_isolate_all_roots, dup_isolate_all_roots_sqf,
                            dup_isolate_complex_roots_sqf,
                            dup_isolate_real_roots,
                            dup_isolate_real_roots_pair,
                            dup_isolate_real_roots_sqf, dup_refine_real_root,
                            dup_root_upper_bound, dup_sign_variations,
                            dup_sturm)


__all__ = 'IPolys',


class IPolys:
    """Compatibility class between dense and sparse polynomials."""

    symbols = None
    ngens = None
    domain = None
    order = None
    gens = None

    def dmp_ground_TC(self, f):
        return dmp_ground_TC(f.to_dense(), self.ngens-1, self.domain)

    def dmp_add_term(self, f, c, i):
        c = c.drop(0).to_dense() if self.is_multivariate else c
        return self.from_list(dmp_add_term(f.to_dense(), c, i, self.ngens-1, self.domain))

    def dmp_mul_term(self, f, c, i):
        c = c.drop(0).to_dense() if self.is_multivariate else c
        return self.from_list(dmp_mul_term(f.to_dense(), c, i, self.ngens-1, self.domain))

    def dmp_mul_ground(self, f, c):
        return self.from_list(dmp_mul_ground(f.to_dense(), c, self.ngens-1, self.domain))

    def dmp_quo_ground(self, f, c):
        return self.from_list(dmp_quo_ground(f.to_dense(), c, self.ngens-1, self.domain))

    def dmp_exquo_ground(self, f, c):
        return self.from_list(dmp_exquo_ground(f.to_dense(), c, self.ngens-1, self.domain))

    def dup_lshift(self, f, n):
        return self.from_list(dup_lshift(f.to_dense(), n, self.domain))

    def dup_rshift(self, f, n):
        return self.from_list(dup_rshift(f.to_dense(), n, self.domain))

    def dmp_abs(self, f):
        return self.from_list(dmp_abs(f.to_dense(), self.ngens-1, self.domain))

    def dmp_neg(self, f):
        return self.from_list(dmp_neg(f.to_dense(), self.ngens-1, self.domain))

    def dmp_add(self, f, g):
        return self.from_list(dmp_add(f.to_dense(), g.to_dense(), self.ngens-1, self.domain))

    def dmp_sub(self, f, g):
        return self.from_list(dmp_sub(f.to_dense(), g.to_dense(), self.ngens-1, self.domain))

    def dmp_add_mul(self, f, g, h):
        return self.from_list(dmp_add_mul(f.to_dense(), g.to_dense(), h.to_dense(), self.ngens-1, self.domain))

    def dmp_mul(self, f, g):
        return self.from_list(dmp_mul(f.to_dense(), g.to_dense(), self.ngens-1, self.domain))

    def dmp_LC(self, f):
        LC = dmp_LC(f.to_dense(), self.domain)
        if self.is_multivariate:
            return self.drop(0).from_list(LC)
        else:
            return LC

    def dmp_TC(self, f):
        TC = dmp_TC(f.to_dense(), self.domain)
        if self.is_multivariate:
            return self.drop(0).from_list(TC)
        else:
            return TC

    def dmp_degree_in(self, f, j):
        return dmp_degree_in(f.to_dense(), j, self.ngens-1)

    def dmp_degree_list(self, f):
        return dmp_degree_list(f.to_dense(), self.ngens-1)

    def dup_real_imag(self, f):
        ring = self
        p, q = dup_real_imag(f.drop(1).to_dense(), ring.domain)
        if ring.domain.is_ComplexAlgebraicField and not ring.domain.is_RealAlgebraicField:
            ring = ring.to_ground()
        return ring.from_list(p), ring.from_list(q)

    def dup_transform(self, f, p, q):
        return self.from_list(dup_transform(f.to_dense(), p.to_dense(), q.to_dense(), self.domain))

    def dmp_compose(self, f, g):
        return self.from_list(dmp_compose(f.to_dense(), g.to_dense(), self.ngens-1, self.domain))

    def dup_decompose(self, f):
        components = dup_decompose(f.to_dense(), self.domain)
        return list(map(self.from_list, components))

    def dup_sign_variations(self, f):
        return dup_sign_variations(f.to_dense(), self.domain)

    # E: List[ZZ], cs: ZZ, ct: ZZ
    def dmp_zz_wang_non_divisors(self, E, cs, ct):
        return dmp_zz_wang_non_divisors(E, cs, ct, self.domain)

    # f: Poly, T: List[(Poly, int)], cs: ZZ, E: List[ZZ], H: List[Poly], A: List[ZZ]
    def dmp_zz_wang_lead_coeffs(self, f, T, cs, E, H, A):
        mv = self.drop(0)
        T = [(t.to_dense(), k) for t, k in T]
        uv = self.drop(*range(1, self.ngens))
        H = list(map(lambda _: _.to_dense(), H))
        f, HH, CC = dmp_zz_wang_lead_coeffs(f.to_dense(), T, cs, E, H, A, self.ngens-1, self.domain)
        return self.from_list(f), list(map(uv.from_list, HH)), list(map(mv.from_list, CC))

    # f: Poly, H: List[Poly], LC: List[Poly], A: List[ZZ], p: ZZ
    def dmp_zz_wang_hensel_lifting(self, f, H, LC, A, p):
        H = list(map(lambda _: _.to_dense(), H))
        LC = list(map(lambda _: _.to_dense(), LC))
        result = dmp_zz_wang_hensel_lifting(f.to_dense(), H, LC, A, p, self.ngens-1, self.domain)
        return list(map(self.from_list, result))

    def dup_sturm(self, f):
        seq = dup_sturm(f.to_dense(), self.domain)
        return list(map(self.from_list, seq))

    def dup_root_upper_bound(self, f):
        return dup_root_upper_bound(f.to_dense(), self.domain)

    def dup_refine_real_root(self, f, s, t, eps=None, steps=None, disjoint=None):
        return dup_refine_real_root(f.to_dense(), s, t, self.domain, eps=eps, steps=steps, disjoint=disjoint)

    def dup_isolate_real_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_real_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_real_roots(self, f, eps=None, inf=None, sup=None):
        return dup_isolate_real_roots(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup)

    def dup_isolate_real_roots_pair(self, f, g, eps=None, inf=None, sup=None, strict=False, basis=False):
        return dup_isolate_real_roots_pair(*(map(lambda _: _.to_dense(), [f, g])), self.domain, eps=eps, inf=inf, sup=sup, strict=strict, basis=basis)

    def dup_count_real_roots(self, f, inf=None, sup=None):
        return dup_count_real_roots(f.to_dense(), self.domain, inf=inf, sup=sup)

    def dup_count_complex_roots(self, f, inf=None, sup=None, exclude=None):
        return dup_count_complex_roots(f.to_dense(), self.domain, inf=inf, sup=sup, exclude=exclude)

    def dup_isolate_complex_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_complex_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_all_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_all_roots_sqf(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_all_roots(self, f, eps=None, inf=None, sup=None):
        return dup_isolate_all_roots(f.to_dense(), self.domain, eps=eps, inf=inf, sup=sup)

    def fateman_poly_F_1(self):
        from .specialpolys import dmp_fateman_poly_F_1
        return tuple(map(self.from_list, dmp_fateman_poly_F_1(self.ngens-1, self.domain)))

    def fateman_poly_F_2(self):
        from .specialpolys import dmp_fateman_poly_F_2
        return tuple(map(self.from_list, dmp_fateman_poly_F_2(self.ngens-1, self.domain)))

    def fateman_poly_F_3(self):
        from .specialpolys import dmp_fateman_poly_F_3
        return tuple(map(self.from_list, dmp_fateman_poly_F_3(self.ngens-1, self.domain)))
