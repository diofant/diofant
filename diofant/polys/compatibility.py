"""Compatibility interface between dense and sparse polys. """

from .densearith import (dmp_abs, dmp_add, dmp_add_mul, dmp_add_term, dmp_div,
                         dmp_expand, dmp_exquo, dmp_exquo_ground, dmp_ff_div,
                         dmp_l1_norm, dmp_max_norm, dmp_mul, dmp_mul_ground,
                         dmp_mul_term, dmp_neg, dmp_pdiv, dmp_pexquo, dmp_pow,
                         dmp_pquo, dmp_prem, dmp_quo, dmp_quo_ground, dmp_rem,
                         dmp_rr_div, dmp_sqr, dmp_sub, dmp_sub_mul,
                         dmp_sub_term, dup_add, dup_add_term, dup_lshift,
                         dup_mul, dup_mul_term, dup_pexquo, dup_rshift,
                         dup_sqr, dup_sub, dup_sub_term)
from .densebasic import dmp_degree_in, dmp_LC, dmp_strip, dmp_to_dict
from .densetools import (dmp_clear_denoms, dmp_compose, dmp_diff_eval_in,
                         dmp_diff_in, dmp_eval_in, dmp_eval_tail,
                         dmp_ground_content, dmp_ground_extract,
                         dmp_ground_monic, dmp_ground_primitive,
                         dmp_ground_trunc, dmp_integrate_in, dmp_lift,
                         dmp_trunc, dup_decompose, dup_mirror, dup_real_imag,
                         dup_scale, dup_shift, dup_sign_variations,
                         dup_transform, dup_trunc)
from .euclidtools import (dmp_cancel, dmp_content, dmp_discriminant,
                          dmp_ff_lcm, dmp_ff_prs_gcd, dmp_gcd, dmp_inner_gcd,
                          dmp_inner_subresultants, dmp_lcm, dmp_primitive,
                          dmp_prs_resultant, dmp_qq_collins_resultant,
                          dmp_qq_heu_gcd, dmp_resultant, dmp_rr_lcm,
                          dmp_rr_prs_gcd, dmp_subresultants,
                          dmp_zz_collins_resultant, dmp_zz_heu_gcd,
                          dmp_zz_modular_resultant, dup_euclidean_prs,
                          dup_ff_prs_gcd, dup_gcdex, dup_half_gcdex,
                          dup_inner_subresultants, dup_invert,
                          dup_primitive_prs, dup_prs_resultant, dup_resultant,
                          dup_rr_prs_gcd)
from .factortools import (dmp_ext_factor, dmp_factor_list, dmp_trial_division,
                          dmp_zz_factor, dmp_zz_mignotte_bound, dmp_zz_wang,
                          dmp_zz_wang_hensel_lifting, dmp_zz_wang_lead_coeffs,
                          dmp_zz_wang_non_divisors, dup_cyclotomic_p,
                          dup_zz_cyclotomic_factor, dup_zz_cyclotomic_poly,
                          dup_zz_factor, dup_zz_factor_sqf, dup_zz_hensel_lift,
                          dup_zz_hensel_step, dup_zz_irreducible_p)
from .galoistools import gf_factor_sqf
from .rootisolation import (dup_count_complex_roots, dup_count_real_roots,
                            dup_isolate_all_roots, dup_isolate_all_roots_sqf,
                            dup_isolate_complex_roots_sqf,
                            dup_isolate_imaginary_roots,
                            dup_isolate_real_roots,
                            dup_isolate_real_roots_list,
                            dup_isolate_real_roots_sqf, dup_refine_real_root,
                            dup_root_upper_bound, dup_sturm)
from .sqfreetools import (dmp_sqf_list, dmp_sqf_list_include, dmp_sqf_norm,
                          dmp_sqf_p, dmp_sqf_part)


__all__ = ('IPolys',)


class IPolys:
    symbols = None
    ngens = None
    domain = None
    order = None
    gens = None

    def wrap(self, element):
        from .rings import PolyElement
        if isinstance(element, PolyElement):
            if element.ring == self:
                return element
            else:
                raise NotImplementedError("domain conversions")
        else:
            return self.ground_new(element)

    def to_dense(self, element):
        return self.wrap(element).to_dense()

    def from_dense(self, element):
        return self.from_dict(dmp_to_dict(element, self.ngens-1, self.domain))

    def dup_add_term(self, f, c, i):
        return self.from_dense(dup_add_term(self.to_dense(f), c, i, self.domain))

    def dmp_add_term(self, f, c, i):
        return self.from_dense(dmp_add_term(self.to_dense(f), self.wrap(c).drop(0).to_dense(), i, self.ngens-1, self.domain))

    def dup_sub_term(self, f, c, i):
        return self.from_dense(dup_sub_term(self.to_dense(f), c, i, self.domain))

    def dmp_sub_term(self, f, c, i):
        return self.from_dense(dmp_sub_term(self.to_dense(f), self.wrap(c).drop(0).to_dense(), i, self.ngens-1, self.domain))

    def dup_mul_term(self, f, c, i):
        return self.from_dense(dup_mul_term(self.to_dense(f), c, i, self.domain))

    def dmp_mul_term(self, f, c, i):
        return self.from_dense(dmp_mul_term(self.to_dense(f), self.wrap(c).drop(0).to_dense(), i, self.ngens-1, self.domain))

    def dmp_mul_ground(self, f, c):
        return self.from_dense(dmp_mul_ground(self.to_dense(f), c, self.ngens-1, self.domain))

    def dmp_quo_ground(self, f, c):
        return self.from_dense(dmp_quo_ground(self.to_dense(f), c, self.ngens-1, self.domain))

    def dmp_exquo_ground(self, f, c):
        return self.from_dense(dmp_exquo_ground(self.to_dense(f), c, self.ngens-1, self.domain))

    def dup_lshift(self, f, n):
        return self.from_dense(dup_lshift(self.to_dense(f), n, self.domain))

    def dup_rshift(self, f, n):
        return self.from_dense(dup_rshift(self.to_dense(f), n, self.domain))

    def dmp_abs(self, f):
        return self.from_dense(dmp_abs(self.to_dense(f), self.ngens-1, self.domain))

    def dmp_neg(self, f):
        return self.from_dense(dmp_neg(self.to_dense(f), self.ngens-1, self.domain))

    def dup_add(self, f, g):
        return self.from_dense(dup_add(self.to_dense(f), self.to_dense(g), self.domain))

    def dmp_add(self, f, g):
        return self.from_dense(dmp_add(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_sub(self, f, g):
        return self.from_dense(dup_sub(self.to_dense(f), self.to_dense(g), self.domain))

    def dmp_sub(self, f, g):
        return self.from_dense(dmp_sub(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dmp_add_mul(self, f, g, h):
        return self.from_dense(dmp_add_mul(self.to_dense(f), self.to_dense(g), self.to_dense(h), self.ngens-1, self.domain))

    def dmp_sub_mul(self, f, g, h):
        return self.from_dense(dmp_sub_mul(self.to_dense(f), self.to_dense(g), self.to_dense(h), self.ngens-1, self.domain))

    def dup_mul(self, f, g):
        return self.from_dense(dup_mul(self.to_dense(f), self.to_dense(g), self.domain))

    def dmp_mul(self, f, g):
        return self.from_dense(dmp_mul(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_sqr(self, f):
        return self.from_dense(dup_sqr(self.to_dense(f), self.domain))

    def dmp_sqr(self, f):
        return self.from_dense(dmp_sqr(self.to_dense(f), self.ngens-1, self.domain))

    def dmp_pow(self, f, n):
        return self.from_dense(dmp_pow(self.to_dense(f), n, self.ngens-1, self.domain))

    def dup_pexquo(self, f, g):
        return self.from_dense(dup_pexquo(self.to_dense(f), self.to_dense(g), self.domain))

    def dmp_pdiv(self, f, g):
        q, r = dmp_pdiv(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(q), self.from_dense(r)

    def dmp_prem(self, f, g):
        return self.from_dense(dmp_prem(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dmp_pquo(self, f, g):
        return self.from_dense(dmp_pquo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dmp_pexquo(self, f, g):
        return self.from_dense(dmp_pexquo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dmp_rr_div(self, f, g):
        q, r = dmp_rr_div(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(q), self.from_dense(r)

    def dmp_ff_div(self, f, g):
        q, r = dmp_ff_div(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(q), self.from_dense(r)

    def dmp_div(self, f, g):
        q, r = dmp_div(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(q), self.from_dense(r)

    def dmp_rem(self, f, g):
        return self.from_dense(dmp_rem(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dmp_quo(self, f, g):
        return self.from_dense(dmp_quo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dmp_exquo(self, f, g):
        return self.from_dense(dmp_exquo(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dmp_max_norm(self, f):
        return dmp_max_norm(self.to_dense(f), self.ngens-1, self.domain)

    def dmp_l1_norm(self, f):
        return dmp_l1_norm(self.to_dense(f), self.ngens-1, self.domain)

    def dmp_expand(self, polys):
        return self.from_dense(dmp_expand(list(map(self.to_dense, polys)), self.ngens-1, self.domain))

    def dmp_LC(self, f):
        LC = dmp_LC(self.to_dense(f), self.domain)
        if self.ngens > 1:
            return self.drop(0).from_dense(LC)
        else:
            return LC

    def dmp_degree_in(self, f, j):
        return dmp_degree_in(self.to_dense(f), j, self.ngens-1)

    def dmp_diff_in(self, f, m, j):
        return self.from_dense(dmp_diff_in(self.to_dense(f), m, j, self.ngens-1, self.domain))

    def dmp_integrate_in(self, f, m, j):
        return self.from_dense(dmp_integrate_in(self.to_dense(f), m, j, self.ngens-1, self.domain))

    def dmp_eval_in(self, f, a, j):
        result = dmp_eval_in(self.to_dense(f), a, j, self.ngens-1, self.domain)
        return self.drop(j).from_dense(result)

    def dmp_diff_eval_in(self, f, m, a, j):
        result = dmp_diff_eval_in(self.to_dense(f), m, a, j, self.ngens-1, self.domain)
        return self.drop(j).from_dense(result)

    def dmp_eval_tail(self, f, A):
        result = dmp_eval_tail(self.to_dense(f), A, self.ngens-1, self.domain)
        if isinstance(result, list):
            return self.drop(*range(self.ngens)[-len(A):]).from_dense(result)
        else:
            return result

    def dup_trunc(self, f, p):
        return self.from_dense(dup_trunc(self.to_dense(f), p, self.domain))

    def dmp_trunc(self, f, g):
        return self.from_dense(dmp_trunc(self.to_dense(f), self.drop(0).to_dense(g), self.ngens-1, self.domain))

    def dmp_ground_trunc(self, f, p):
        return self.from_dense(dmp_ground_trunc(self.to_dense(f), p, self.ngens-1, self.domain))

    def dmp_ground_monic(self, f):
        return self.from_dense(dmp_ground_monic(self.to_dense(f), self.ngens-1, self.domain))

    def dmp_ground_extract(self, f, g):
        c, F, G = dmp_ground_extract(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return c, self.from_dense(F), self.from_dense(G)

    def dup_real_imag(self, f):
        ring = self
        p, q = dup_real_imag(ring.wrap(f).drop(1).to_dense(), ring.domain)
        if ring.domain.is_ComplexAlgebraicField:
            ring = ring.to_ground()
        return ring.from_dense(p), ring.from_dense(q)

    def dup_mirror(self, f):
        return self.from_dense(dup_mirror(self.to_dense(f), self.domain))

    def dup_scale(self, f, a):
        return self.from_dense(dup_scale(self.to_dense(f), a, self.domain))

    def dup_shift(self, f, a):
        return self.from_dense(dup_shift(self.to_dense(f), a, self.domain))

    def dup_transform(self, f, p, q):
        return self.from_dense(dup_transform(self.to_dense(f), self.to_dense(p), self.to_dense(q), self.domain))

    def dmp_compose(self, f, g):
        return self.from_dense(dmp_compose(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain))

    def dup_decompose(self, f):
        components = dup_decompose(self.to_dense(f), self.domain)
        return list(map(self.from_dense, components))

    def dmp_lift(self, f):
        result = dmp_lift(self.to_dense(f), self.ngens-1, self.domain)
        return self.to_ground().from_dense(result)

    def dup_sign_variations(self, f):
        return dup_sign_variations(self.to_dense(f), self.domain)

    def dmp_clear_denoms(self, f, convert=False):
        c, F = dmp_clear_denoms(self.to_dense(f), self.ngens-1, self.domain, convert=convert)
        if convert:
            ring = self.clone(domain=self.domain.ring)
        else:
            ring = self
        return c, ring.from_dense(F)

    def dup_half_gcdex(self, f, g):
        s, h = dup_half_gcdex(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(s), self.from_dense(h)

    def dup_gcdex(self, f, g):
        s, t, h = dup_gcdex(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(s), self.from_dense(t), self.from_dense(h)

    def dup_invert(self, f, g):
        return self.from_dense(dup_invert(self.to_dense(f), self.to_dense(g), self.domain))

    def dup_euclidean_prs(self, f, g):
        prs = dup_euclidean_prs(self.to_dense(f), self.to_dense(g), self.domain)
        return list(map(self.from_dense, prs))

    def dup_primitive_prs(self, f, g):
        prs = dup_primitive_prs(self.to_dense(f), self.to_dense(g), self.domain)
        return list(map(self.from_dense, prs))

    def dup_inner_subresultants(self, f, g):
        prs, sres = dup_inner_subresultants(self.to_dense(f), self.to_dense(g), self.domain)
        return list(map(self.from_dense, prs)), sres

    def dmp_inner_subresultants(self, f, g):
        prs, sres = dmp_inner_subresultants(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return list(map(self.from_dense, prs)), sres

    def dmp_subresultants(self, f, g):
        prs = dmp_subresultants(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return list(map(self.from_dense, prs))

    def dup_prs_resultant(self, f, g):
        res, prs = dup_prs_resultant(self.to_dense(f), self.to_dense(g), self.domain)
        return res, list(map(self.from_dense, prs))

    def dmp_prs_resultant(self, f, g):
        res, prs = dmp_prs_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        if isinstance(res, list):
            res = self.drop(0).from_dense(res)
        return res, list(map(self.from_dense, prs))

    def dmp_zz_modular_resultant(self, f, g, p):
        res = dmp_zz_modular_resultant(self.to_dense(f), self.to_dense(g), self.domain_new(p), self.ngens-1, self.domain)
        return self.drop(0).from_dense(res)

    def dmp_zz_collins_resultant(self, f, g):
        res = dmp_zz_collins_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.drop(0).from_dense(res)

    def dmp_qq_collins_resultant(self, f, g):
        res = dmp_qq_collins_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.drop(0).from_dense(res)

    def dup_resultant(self, f, g):
        return dup_resultant(self.to_dense(f), self.to_dense(g), self.domain)

    def dmp_resultant(self, f, g, includePRS=False):
        res = dmp_resultant(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain, includePRS=includePRS)
        if isinstance(res, list):
            return self.drop(0).from_dense(res)
        else:
            return res

    def dmp_discriminant(self, f):
        disc = dmp_discriminant(self.to_dense(f), self.ngens-1, self.domain)
        if isinstance(disc, list):
            return self.drop(0).from_dense(disc)
        else:
            return disc

    def dup_rr_prs_gcd(self, f, g):
        H, F, G = dup_rr_prs_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(H), self.from_dense(F), self.from_dense(G)

    def dup_ff_prs_gcd(self, f, g):
        H, F, G = dup_ff_prs_gcd(self.to_dense(f), self.to_dense(g), self.domain)
        return self.from_dense(H), self.from_dense(F), self.from_dense(G)

    def dmp_rr_prs_gcd(self, f, g):
        H, F, G = dmp_rr_prs_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H), self.from_dense(F), self.from_dense(G)

    def dmp_ff_prs_gcd(self, f, g):
        H, F, G = dmp_ff_prs_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H), self.from_dense(F), self.from_dense(G)

    def dmp_zz_heu_gcd(self, f, g):
        H, F, G = dmp_zz_heu_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H), self.from_dense(F), self.from_dense(G)

    def dmp_qq_heu_gcd(self, f, g):
        H, F, G = dmp_qq_heu_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H), self.from_dense(F), self.from_dense(G)

    def dmp_inner_gcd(self, f, g):
        H, F, G = dmp_inner_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H), self.from_dense(F), self.from_dense(G)

    def dmp_gcd(self, f, g):
        H = dmp_gcd(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)

    def dmp_rr_lcm(self, f, g):
        H = dmp_rr_lcm(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)

    def dmp_ff_lcm(self, f, g):
        H = dmp_ff_lcm(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)

    def dmp_lcm(self, f, g):
        H = dmp_lcm(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain)
        return self.from_dense(H)

    def dmp_content(self, f):
        cont = dmp_content(self.to_dense(f), self.ngens-1, self.domain)
        return self.drop(0).from_dense(cont)

    def dmp_primitive(self, f):
        cont, prim = dmp_primitive(self.to_dense(f), self.ngens-1, self.domain)
        return self.drop(0).from_dense(cont), self.from_dense(prim)

    def dmp_ground_content(self, f):
        cont = dmp_ground_content(self.to_dense(f), self.ngens-1, self.domain)
        return cont

    def dmp_ground_primitive(self, f):
        cont, prim = dmp_ground_primitive(self.to_dense(f), self.ngens-1, self.domain)
        return cont, self.from_dense(prim)

    def dmp_cancel(self, f, g, include=True):
        result = dmp_cancel(self.to_dense(f), self.to_dense(g), self.ngens-1, self.domain, include=include)
        if not include:
            cf, cg, F, G = result
            return cf, cg, self.from_dense(F), self.from_dense(G)
        else:
            F, G = result
            return self.from_dense(F), self.from_dense(G)

    def dmp_trial_division(self, f, factors):
        factors = dmp_trial_division(self.to_dense(f), list(map(self.to_dense, factors)), self.ngens-1, self.domain)
        return [ (self.from_dense(g), k) for g, k in factors ]

    def dmp_zz_mignotte_bound(self, f):
        return dmp_zz_mignotte_bound(self.to_dense(f), self.ngens-1, self.domain)

    def dup_zz_hensel_step(self, m, f, g, h, s, t):
        D = self.to_dense
        G, H, S, T = dup_zz_hensel_step(m, D(f), D(g), D(h), D(s), D(t), self.domain)
        return self.from_dense(G), self.from_dense(H), self.from_dense(S), self.from_dense(T)

    def dup_zz_hensel_lift(self, p, f, f_list, l):
        D = self.to_dense
        polys = dup_zz_hensel_lift(p, D(f), list(map(D, f_list)), l, self.domain)
        return list(map(self.from_dense, polys))

    def dup_zz_irreducible_p(self, f):
        return dup_zz_irreducible_p(self.to_dense(f), self.domain)

    def dup_cyclotomic_p(self, f, irreducible=False):
        return dup_cyclotomic_p(self.to_dense(f), self.domain, irreducible=irreducible)

    def dup_zz_cyclotomic_poly(self, n):
        F = dup_zz_cyclotomic_poly(n, self.domain)
        return self.from_dense(F)

    def dup_zz_cyclotomic_factor(self, f):
        result = dup_zz_cyclotomic_factor(self.to_dense(f), self.domain)
        if result is None:
            return result
        else:
            return list(map(self.from_dense, result))

    # E: List[ZZ], cs: ZZ, ct: ZZ
    def dmp_zz_wang_non_divisors(self, E, cs, ct):
        return dmp_zz_wang_non_divisors(E, cs, ct, self.domain)

    # f: Poly, T: List[(Poly, int)], ct: ZZ, A: List[ZZ]
    # def dmp_zz_wang_test_points(f, T, ct, A):
    #   dmp_zz_wang_test_points(self.to_dense(f), T, ct, A, self.ngens-1, self.domain)

    # f: Poly, T: List[(Poly, int)], cs: ZZ, E: List[ZZ], H: List[Poly], A: List[ZZ]
    def dmp_zz_wang_lead_coeffs(self, f, T, cs, E, H, A):
        mv = self.drop(0)
        T = [ (mv.to_dense(t), k) for t, k in T ]
        uv = self.drop(*range(1, self.ngens))
        H = list(map(uv.to_dense, H))
        f, HH, CC = dmp_zz_wang_lead_coeffs(self.to_dense(f), T, cs, E, H, A, self.ngens-1, self.domain)
        return self.from_dense(f), list(map(uv.from_dense, HH)), list(map(mv.from_dense, CC))

    # f: Poly, H: List[Poly], LC: List[Poly], A: List[ZZ], p: ZZ
    def dmp_zz_wang_hensel_lifting(self, f, H, LC, A, p):
        uv = self.drop(*range(1, self.ngens))
        mv = self.drop(0)
        H = list(map(uv.to_dense, H))
        LC = list(map(mv.to_dense, LC))
        result = dmp_zz_wang_hensel_lifting(self.to_dense(f), H, LC, A, p, self.ngens-1, self.domain)
        return list(map(self.from_dense, result))

    def dmp_zz_wang(self, f, mod=None, seed=None):
        factors = dmp_zz_wang(self.to_dense(f), self.ngens-1, self.domain, mod=mod, seed=seed)
        return [ self.from_dense(g) for g in factors ]

    def dup_zz_factor_sqf(self, f):
        coeff, factors = dup_zz_factor_sqf(self.to_dense(f), self.domain)
        return coeff, [self.from_dense(g) for g in factors]

    def dup_zz_factor(self, f):
        coeff, factors = dup_zz_factor(self.to_dense(f), self.domain)
        return coeff, [(self.from_dense(g), k) for g, k in factors]

    def dmp_zz_factor(self, f):
        coeff, factors = dmp_zz_factor(self.to_dense(f), self.ngens-1, self.domain)
        return coeff, [(self.from_dense(g), k) for g, k in factors]

    def dmp_ext_factor(self, f):
        coeff, factors = dmp_ext_factor(self.to_dense(f), self.ngens-1, self.domain)
        return coeff, [(self.from_dense(g), k) for g, k in factors]

    def dmp_factor_list(self, f):
        coeff, factors = dmp_factor_list(self.to_dense(f), self.ngens-1, self.domain)
        return coeff, [(self.from_dense(g), k) for g, k in factors]

    def dup_sturm(self, f):
        seq = dup_sturm(self.to_dense(f), self.domain)
        return list(map(self.from_dense, seq))

    def dmp_sqf_p(self, f):
        return dmp_sqf_p(self.to_dense(f), self.ngens-1, self.domain)

    def dmp_sqf_norm(self, f):
        s, F, R = dmp_sqf_norm(self.to_dense(f), self.ngens-1, self.domain)
        return s, self.from_dense(F), self.to_ground().from_dense(R)

    def dmp_sqf_part(self, f):
        return self.from_dense(dmp_sqf_part(self.to_dense(f), self.ngens-1, self.domain))

    def dmp_sqf_list(self, f):
        coeff, factors = dmp_sqf_list(self.to_dense(f), self.ngens-1, self.domain)
        return coeff, [(self.from_dense(g), k) for g, k in factors]

    def dmp_sqf_list_include(self, f):
        factors = dmp_sqf_list_include(self.to_dense(f), self.ngens-1, self.domain)
        return [(self.from_dense(g), k) for g, k in factors]

    def dup_root_upper_bound(self, f):
        return dup_root_upper_bound(self.to_dense(f), self.domain)

    def dup_refine_real_root(self, f, s, t, eps=None, steps=None, disjoint=None, fast=False):
        return dup_refine_real_root(self.to_dense(f), s, t, self.domain, eps=eps, steps=steps, disjoint=disjoint, fast=fast)

    def dup_isolate_real_roots_sqf(self, f, eps=None, inf=None, sup=None, fast=False, blackbox=False):
        return dup_isolate_real_roots_sqf(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast, blackbox=blackbox)

    def dup_isolate_real_roots(self, f, eps=None, inf=None, sup=None, fast=False):
        return dup_isolate_real_roots(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast)

    def dup_isolate_imaginary_roots(self, f, eps=None, inf=None, sup=None, fast=False):
        return dup_isolate_imaginary_roots(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast)

    def dup_isolate_real_roots_list(self, polys, eps=None, inf=None, sup=None, strict=False, basis=False, fast=False):
        return dup_isolate_real_roots_list(list(map(self.to_dense, polys)), self.domain, eps=eps, inf=inf, sup=sup, strict=strict, basis=basis, fast=fast)

    def dup_count_real_roots(self, f, inf=None, sup=None):
        return dup_count_real_roots(self.to_dense(f), self.domain, inf=inf, sup=sup)

    def dup_count_complex_roots(self, f, inf=None, sup=None, exclude=None):
        return dup_count_complex_roots(self.to_dense(f), self.domain, inf=inf, sup=sup, exclude=exclude)

    def dup_isolate_complex_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        return dup_isolate_complex_roots_sqf(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, blackbox=blackbox)

    def dup_isolate_all_roots_sqf(self, f, eps=None, inf=None, sup=None, fast=False, blackbox=False):
        return dup_isolate_all_roots_sqf(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast, blackbox=blackbox)

    def dup_isolate_all_roots(self, f, eps=None, inf=None, sup=None, fast=False):
        return dup_isolate_all_roots(self.to_dense(f), self.domain, eps=eps, inf=inf, sup=sup, fast=fast)

    def fateman_poly_F_1(self):
        from .specialpolys import dmp_fateman_poly_F_1
        return tuple(map(self.from_dense, dmp_fateman_poly_F_1(self.ngens-1, self.domain)))

    def fateman_poly_F_2(self):
        from .specialpolys import dmp_fateman_poly_F_2
        return tuple(map(self.from_dense, dmp_fateman_poly_F_2(self.ngens-1, self.domain)))

    def fateman_poly_F_3(self):
        from .specialpolys import dmp_fateman_poly_F_3
        return tuple(map(self.from_dense, dmp_fateman_poly_F_3(self.ngens-1, self.domain)))

    def to_gf_dense(self, element):
        return dmp_strip([self.domain.domain.convert(c, self.domain) for c in self.wrap(element).to_dense()], 0)

    def from_gf_dense(self, element):
        return self.from_dict(dmp_to_dict(element, self.ngens-1, self.domain.domain))

    def gf_factor_sqf(self, f):
        coeff, factors = gf_factor_sqf(self.to_gf_dense(f), self.domain.mod, self.domain.domain)
        return coeff, [ self.from_gf_dense(g) for g in factors ]
