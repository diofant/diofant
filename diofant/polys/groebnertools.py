"""Gröbner bases algorithms."""

from ..config import query
from ..core import Dummy
from .monomials import Monomial
from .orderings import lex


def groebner(seq, ring, method=None):
    """
    Computes Gröbner basis for a set of polynomials in `K[X]`.

    Wrapper around the (default) improved Buchberger and the other algorithms
    for computing Gröbner bases. The choice of algorithm can be changed via
    ``method`` argument or :func:`~diofant.config.setup`,
    where ``method`` can be either ``buchberger`` or ``f5b``.

    """
    if method is None:
        method = query('groebner')

    _groebner_methods = {
        'buchberger': buchberger,
        'f5b': f5b,
    }

    try:
        _groebner = _groebner_methods[method]
    except KeyError as exc:
        raise ValueError(f"'{method}' is not a valid Gröbner "
                         "bases algorithm (valid are 'buchberger'"
                         " and 'f5b')") from exc

    domain, orig = ring.domain, None

    if not domain.is_Field and hasattr(domain, 'field'):
        orig, ring = ring, ring.clone(domain=domain.field)
        seq = [s.set_ring(ring) for s in seq]

    G = _groebner(seq, ring)

    if orig is not None:
        G = [g.clear_denoms()[1].set_ring(orig) for g in G]

    return G


def buchberger(f, ring):
    """
    Computes Gröbner basis for a set of polynomials in `K[X]`.

    Given a set of multivariate polynomials `F`, finds another
    set `G`, such that Ideal `F = Ideal G` and `G` is a reduced
    Gröbner basis.

    The resulting basis is unique and has monic generators if the
    ground domains is a field. Otherwise the result is non-unique
    but Gröbner bases over e.g. integers can be computed (if the
    input polynomials are monic).

    Gröbner bases can be used to choose specific generators for a
    polynomial ideal. Because these bases are unique you can check
    for ideal equality by comparing the Gröbner bases.  To see if
    one polynomial lies in an ideal, divide by the elements in the
    base and see if the remainder vanishes.

    They can also be used to solve systems of polynomial equations
    as,  by choosing lexicographic ordering,  you can eliminate one
    variable at a time, provided that the ideal is zero-dimensional
    (finite number of solutions).

    References
    ==========

    * :cite:`Bose03`
    * :cite:`Giovini1991sugar`
    * :cite:`Ajwa95groebner`
    * :cite:`Cox2015ideals`
    * :cite:`BeckerWeispfenning93`, page 232

    Notes
    =====

    Used an improved version of Buchberger's algorithm
    as presented in :cite:`BeckerWeispfenning93`.

    """
    order = ring.order

    def select(P):
        # normal selection strategy
        # select the pair with minimum LCM(LM(f), LM(g))
        pr = min(P, key=lambda pair: order(f[pair[0]].LM.lcm(f[pair[1]].LM)))
        return pr

    def normal(g, J):
        h = g.div([f[j] for j in J])[1]

        if h:
            h = h.monic()

            if h not in I:
                I[h] = len(f)
                f.append(h)

            return h.LM, I[h]

    def update(G, B, ih):
        # update G using the set of critical pairs B and h
        # [BW] page 230
        h = f[ih]
        mh = h.LM

        # filter new pairs (h, g), g in G
        C = G.copy()
        D = set()

        while C:
            # select a pair (h, g) by popping an element from C
            ig = C.pop()
            g = f[ig]
            mg = g.LM
            LCMhg = mh.lcm(mg)

            def lcm_divides(ip):
                # LCM(LM(h), LM(p)) divides LCM(LM(h), LM(g))
                m = mh.lcm(f[ip].LM)
                return m.divides(LCMhg)

            # HT(h) and HT(g) disjoint: mh*mg == LCMhg
            if mh*mg == LCMhg or (
                not any(lcm_divides(ipx) for ipx in C) and
                    not any(lcm_divides(pr[1]) for pr in D)):
                D.add((ih, ig))

        E = set()

        while D:
            # select h, g from D (h the same as above)
            ih, ig = D.pop()
            mg = f[ig].LM
            LCMhg = mh.lcm(mg)

            if not mh*mg == LCMhg:
                E.add((ih, ig))

        # filter old pairs
        B_new = set()

        while B:
            # select g1, g2 from B (-> CP)
            ig1, ig2 = B.pop()
            mg1 = f[ig1].LM
            mg2 = f[ig2].LM
            LCM12 = mg1.lcm(mg2)

            # if HT(h) does not divide lcm(HT(g1), HT(g2))
            if not mh.divides(LCM12) or mg1.lcm(mh) == LCM12 or mg2.lcm(mh) == LCM12:
                B_new.add((ig1, ig2))

        B_new |= E

        # filter polynomials
        G_new = set()

        while G:
            ig = G.pop()
            mg = f[ig].LM

            if not mh.divides(mg):
                G_new.add(ig)

        G_new.add(ih)

        return G_new, B_new
        # end of update ################################

    if not f:
        return []

    # replace f with a reduced list of initial polynomials; see [BW] page 203
    f1 = f[:]

    while True:
        f = f1[:]
        f1 = []

        for i, p in enumerate(f):
            r = p.div(f[:i])[1]

            if r:
                f1.append(r.monic())

        if f == f1:
            break

    I = {}            # ip = I[p]; p = f[ip]
    F = set()         # set of indices of polynomials
    G = set()         # set of indices of intermediate would-be Gröbner basis
    CP = set()        # set of pairs of indices of critical pairs

    for i, h in enumerate(f):
        I[h] = i
        F.add(i)

    #####################################
    # algorithm GROEBNERNEWS2 in [BW] page 232

    while F:
        # select p with minimum monomial according to the monomial ordering
        h = min((f[x] for x in F), key=lambda f: order(f.LM))
        ih = I[h]
        F.remove(ih)
        G, CP = update(G, CP, ih)

    # count the number of critical pairs which reduce to zero
    reductions_to_zero = 0

    while CP:
        ig1, ig2 = select(CP)
        CP.remove((ig1, ig2))

        h = spoly(f[ig1], f[ig2])
        # ordering divisors is on average more efficient [Cox] page 111
        G1 = sorted(G, key=lambda g: order(f[g].LM))
        ht = normal(h, G1)

        if ht:
            G, CP = update(G, CP, ht[1])
        else:
            reductions_to_zero += 1

    ######################################
    # now G is a Gröbner basis; reduce it
    Gr = set()

    for ig in G:
        ht = normal(f[ig], G - {ig})

        if ht:
            Gr.add(ht[1])

    Gr = [f[ig] for ig in Gr]

    # order according to the monomial ordering
    Gr = sorted(Gr, key=lambda f: order(f.LM), reverse=True)

    return Gr


def spoly(p1, p2):
    """
    Compute LCM(LM(p1), LM(p2))/LM(p1)*p1 - LCM(LM(p1), LM(p2))/LM(p2)*p2.

    This is the S-poly, provided p1 and p2 are monic

    """
    ring = p1.ring
    domain_one = ring.domain.one
    LM1 = p1.LM
    LM2 = p2.LM
    LCM12 = LM1.lcm(LM2)
    m1 = ring.from_terms([(LCM12/LM1, domain_one)])
    m2 = ring.from_terms([(LCM12/LM2, domain_one)])
    s1 = p1*m1
    s2 = p2*m2
    s = s1 - s2
    return s

# F5B

# convenience functions


def Sign(f):
    return f[0]


def Polyn(f):
    return f[1]


def Num(f):
    return f[2]


def sig(monomial, index):
    return monomial, index


def lbp(signature, polynomial, number):
    return signature, polynomial, number

# signature functions


def sig_cmp(u, v, order):
    """
    Compare two signatures by extending the term order to K[X]^n.

    u < v iff
        - the index of v is greater than the index of u
    or
        - the index of v is equal to the index of u and u[0] < v[0] w.r.t. order

    u > v otherwise

    """
    if u[1] > v[1]:
        return -1
    if u[1] == v[1] and order(u[0]) < order(v[0]):
        return -1
    return 1


def sig_key(s, order):
    """
    Key for comparing two signatures.

    s = (m, k), t = (n, l)

    s < t iff [k > l] or [k == l and m < n]
    s > t otherwise

    """
    return -s[1], order(s[0])


def sig_mult(s, m):
    """
    Multiply a signature by a monomial.

    The product of a signature (m, i) and a monomial n is defined as
    (m * t, i).

    """
    return sig(Monomial(s[0])*m, s[1])

# labeled polynomial functions


def lbp_sub(f, g):
    """
    Subtract labeled polynomial g from f.

    The signature and number of the difference of f and g are signature
    and number of the maximum of f and g, w.r.t. lbp_cmp.

    """
    if lbp_cmp(f, g) < 0:
        max_poly = g
    else:
        max_poly = f

    ret = Polyn(f) - Polyn(g)

    return lbp(Sign(max_poly), ret, Num(max_poly))


def lbp_mul_term(f, cx):
    """
    Multiply a labeled polynomial with a term.

    The product of a labeled polynomial (s, p, k) by a monomial is
    defined as (m * s, m * p, k).

    """
    return lbp(sig_mult(Sign(f), cx[0]),
               Polyn(f)*Polyn(f).ring.from_terms([cx]), Num(f))


def lbp_cmp(f, g):
    """
    Compare two labeled polynomials.

    f < g iff
        - Sign(f) < Sign(g)
    or
        - Sign(f) == Sign(g) and Num(f) > Num(g)

    f > g otherwise

    """
    if sig_cmp(Sign(f), Sign(g), Polyn(f).ring.order) == -1:
        return -1
    if Sign(f) == Sign(g) and Num(f) > Num(g):
        return -1
    return 1


def lbp_key(f):
    """Key for comparing two labeled polynomials."""
    return sig_key(Sign(f), Polyn(f).ring.order), -Num(f)

# algorithm and helper functions


def critical_pair(f, g, ring):
    """
    Compute the critical pair corresponding to two labeled polynomials.

    A critical pair is a tuple (um, f, vm, g), where um and vm are
    terms such that um * f - vm * g is the S-polynomial of f and g (so,
    wlog assume um * f > vm * g).
    For performance sake, a critical pair is represented as a tuple
    (Sign(um * f), um, f, Sign(vm * g), vm, g), since um * f creates
    a new, relatively expensive object in memory, whereas Sign(um *
    f) and um are lightweight and f (in the tuple) is a reference to
    an already existing object in memory.

    """
    domain = ring.domain

    ltf = Polyn(f).LT
    ltg = Polyn(g).LT
    lt = ring.from_terms([(Monomial(ltf[0]).lcm(ltg[0]), domain.one)])

    um = lt.quo_term(ltf).LT
    vm = lt.quo_term(ltg).LT

    # The full information is not needed (now), so only the product
    # with the leading term is considered:
    fr = lbp_mul_term(lbp(Sign(f), Polyn(f).leading_term(), Num(f)), um)
    gr = lbp_mul_term(lbp(Sign(g), Polyn(g).leading_term(), Num(g)), vm)

    # return in proper order, such that the S-polynomial is just
    # u_first * f_first - u_second * f_second:
    if lbp_cmp(fr, gr) == -1:
        return Sign(gr), vm, g, Sign(fr), um, f
    return Sign(fr), um, f, Sign(gr), vm, g


def cp_key(c, ring):
    """Key for comparing critical pairs."""
    return lbp_key(lbp(c[0], ring.zero, Num(c[2]))), lbp_key(lbp(c[3], ring.zero, Num(c[5])))


def s_poly(cp):
    """Compute the S-polynomial of a critical pair.

    The S-polynomial of a critical pair cp is cp[1] * cp[2] - cp[4] * cp[5].

    """
    return lbp_sub(lbp_mul_term(cp[2], cp[1]), lbp_mul_term(cp[5], cp[4]))


def is_rewritable_or_comparable(sign, num, B):
    """
    Check if a labeled polynomial is redundant by checking if its
    signature and number imply rewritability or comparability.

    (sign, num) is comparable if there exists a labeled polynomial
    h in B, such that sign[1] (the index) is less than Sign(h)[1]
    and sign[0] is divisible by the leading monomial of h.

    (sign, num) is rewritable if there exists a labeled polynomial
    h in B, such thatsign[1] is equal to Sign(h)[1], num < Num(h)
    and sign[0] is divisible by Sign(h)[0].

    """
    for h in B:
        # comparable
        if sign[1] < Sign(h)[1]:
            if Polyn(h).LM.divides(sign[0]):
                return True

        # rewritable
        if sign[1] == Sign(h)[1]:
            if num < Num(h):
                if Monomial(Sign(h)[0]).divides(sign[0]):
                    return True
    return False


def f5_reduce(f, B):
    """
    F5-reduce a labeled polynomial f by B.

    Continously searches for non-zero labeled polynomial h in B, such
    that the leading term lt_h of h divides the leading term lt_f of
    f and Sign(lt_h * h) < Sign(f). If such a labeled polynomial h is
    found, f gets replaced by f - lt_f / lt_h * h. If no such h can be
    found or f is 0, f is no further F5-reducible and f gets returned.

    A polynomial that is reducible in the usual sense need not be
    F5-reducible, e.g.:

    >>> _, x, y, z = ring('x y z', QQ, lex)

    >>> f = lbp(sig(Monomial((1, 1, 1)), 4), x, 3)
    >>> g = lbp(sig(Monomial((0, 0, 0)), 2), x, 2)

    >>> Polyn(f).div([Polyn(g)])[1]
    0
    >>> f5_reduce(f, [g])
    (((1, 1, 1), 4), x, 3)

    """
    order = Polyn(f).ring.order

    if not Polyn(f):
        return f

    while True:
        g = f

        for h in B:
            if Polyn(h) and Polyn(h).LM.divides(Polyn(f).LM):
                t = Polyn(f).leading_term().quo_term(Polyn(h).LT).LT
                if sig_cmp(sig_mult(Sign(h), t[0]), Sign(f), order) < 0:
                    # The following check need not be done and is in general slower than without.
                    # if not is_rewritable_or_comparable(Sign(gp), Num(gp), B):
                    hp = lbp_mul_term(h, t)
                    f = lbp_sub(f, hp)
                    break

        if g == f or not Polyn(f):
            return f


def f5b(F, ring):
    """
    Computes a reduced Gröbner basis for the ideal generated by F.

    f5b is an implementation of the F5B algorithm by Yao Sun and
    Dingkang Wang. Similarly to Buchberger's algorithm, the algorithm
    proceeds by computing critical pairs, computing the S-polynomial,
    reducing it and adjoining the reduced S-polynomial if it is not 0.

    Unlike Buchberger's algorithm, each polynomial contains additional
    information, namely a signature and a number. The signature
    specifies the path of computation (i.e. from which polynomial in
    the original basis was it derived and how), the number says when
    the polynomial was added to the basis.  With this information it
    is (often) possible to decide if an S-polynomial will reduce to
    0 and can be discarded.

    Optimizations include: Reducing the generators before computing
    a Gröbner basis, removing redundant critical pairs when a new
    polynomial enters the basis and sorting the critical pairs and
    the current basis.

    Once a Gröbner basis has been found, it gets reduced.

    References
    ==========

    * :cite:`SunWang2010f5`
    * :cite:`BeckerWeispfenning93`, pp. 203, 216.

    """
    order = ring.order

    # reduce polynomials (Becker, Weispfenning, p. 203)
    B = F
    while True:
        F = B
        B = []

        for i, p in enumerate(F):
            r = p.div(F[:i])[1]

            if r:
                B.append(r)

        if F == B:
            break

    # basis
    B = [lbp(sig(ring.zero_monom, i + 1), F[i], i + 1) for i in range(len(F))]
    B.sort(key=lambda f: order(Polyn(f).LM), reverse=True)

    # critical pairs
    CP = [critical_pair(B[i], B[j], ring) for i in range(len(B)) for j in range(i + 1, len(B))]
    CP.sort(key=lambda cp: cp_key(cp, ring), reverse=True)

    k = len(B)

    reductions_to_zero = 0

    while len(CP):
        cp = CP.pop()

        # discard redundant critical pairs:
        if any(is_rewritable_or_comparable(x, Num(y), B)
               for x, y in [(cp[0], cp[2]), (cp[3], cp[5])]):
            continue

        s = s_poly(cp)

        p = f5_reduce(s, B)

        p = lbp(Sign(p), Polyn(p).monic(), k + 1)

        if Polyn(p):
            # remove old critical pairs, that become redundant when adding p:
            indices = []
            for i, cp in enumerate(CP):
                if is_rewritable_or_comparable(cp[0], Num(cp[2]), [p]):
                    indices.append(i)
                elif is_rewritable_or_comparable(cp[3], Num(cp[5]), [p]):
                    indices.append(i)

            for i in reversed(indices):
                del CP[i]

            # only add new critical pairs that are not made redundant by p:
            for g in B:
                assert Polyn(g)
                cp = critical_pair(p, g, ring)
                if is_rewritable_or_comparable(cp[0], Num(cp[2]), [p]):
                    continue
                if is_rewritable_or_comparable(cp[3], Num(cp[5]), [p]):
                    continue

                CP.append(cp)

            # sort (other sorting methods/selection strategies were not as successful)
            CP.sort(key=lambda cp: cp_key(cp, ring), reverse=True)

            # insert p into B:
            m = Polyn(p).LM
            if order(m) <= order(Polyn(B[-1]).LM):
                B.append(p)
            else:
                for i, q in enumerate(B):  # pragma: no branch
                    if order(m) > order(Polyn(q).LM):
                        B.insert(i, p)
                        break

            k += 1

        else:
            reductions_to_zero += 1

    # reduce Gröbner basis:
    H = [Polyn(g).monic() for g in B]
    H = red_groebner(H, ring)

    return sorted(H, key=lambda f: order(f.LM), reverse=True)


def red_groebner(G, ring):
    """
    Compute reduced Gröbner basis.

    Selects a subset of generators, that already generate the ideal
    and computes a reduced Gröbner basis for them.

    References
    ==========

    * :cite:`BeckerWeispfenning93`, page 216.

    """
    def reduction(P):
        """The actual reduction algorithm."""
        Q = []
        for i, p in enumerate(P):
            h = p.div(P[:i] + P[i + 1:])[1]
            assert h
            Q.append(h)

        return [p.monic() for p in Q]

    F = G
    H = []

    while F:
        f0 = F.pop()

        if not any(f.LM.divides(f0.LM) for f in F + H):
            H.append(f0)

    # Becker, Weispfenning, p. 217: H is Gröbner basis of the ideal generated by G.
    return reduction(H)


def is_groebner(G):
    """Check if G is a Gröbner basis."""
    for i, Gi in enumerate(G):
        for j in range(i + 1, len(G)):
            s = spoly(Gi, G[j])
            s = s.div(G)[1]
            if s:
                return False

    return True


def is_minimal(G, ring):
    """Checks if G is a minimal Gröbner basis."""
    order = ring.order

    G.sort(key=lambda g: order(g.LM))

    for i, g in enumerate(G):
        if g.LC != 1:
            return False

        for h in G[:i] + G[i + 1:]:
            if h.LM.divides(g.LM):
                return False

    return True


def groebner_lcm(f, g):
    """
    Computes LCM of two polynomials using Gröbner bases.

    The LCM is computed as the unique generator of the intersection
    of the two ideals generated by `f` and `g`. The approach is to
    compute a Gröbner basis with respect to lexicographic ordering
    of `t*f` and `(1 - t)*g`, where `t` is an unrelated variable and
    then filtering out the solution that doesn't contain `t`.

    References
    ==========

    * :cite:`Cox2015ideals`

    """
    if f.ring != g.ring:
        raise ValueError('Values should be equal')

    ring = f.ring
    domain = ring.domain

    if not f or not g:
        return ring.zero

    fc, f = f.primitive()
    gc, g = g.primitive()

    t_ring = ring.clone(order=lex).inject(Dummy('t'), front=True)
    t = t_ring.gens[0]

    basis = groebner([t*f, (1 - t)*g], t_ring)
    H = [h for h in basis if h.degree() <= 0]

    return H[0].drop(0)*domain.lcm(fc, gc)


def groebner_gcd(f, g):
    """Computes GCD of two polynomials using Gröbner bases."""
    if f.ring != g.ring:
        raise ValueError('Values should be equal')
    domain = f.ring.domain

    if not domain.is_Field:
        fc, f = f.primitive()
        gc, g = g.primitive()
        gcd = domain.gcd(fc, gc)

    h = (f*g)//groebner_lcm(f, g)

    if not domain.is_Field:
        return gcd*h
    return h.monic()


def matrix_fglm(F, ring, O_to):
    """
    Converts the reduced Gröbner basis ``F`` of a zero-dimensional
    ideal w.r.t. ``O_from`` to a reduced Gröbner basis
    w.r.t. ``O_to``.

    References
    ==========

    * :cite:`Faugere1993groebner`

    """
    domain = ring.domain
    ngens = ring.ngens

    ring_to = ring.clone(order=O_to)

    old_basis = _basis(F, ring)
    M = _representing_matrices(old_basis, F, ring)

    # V contains the normalforms (wrt O_from) of S
    S = [ring.zero_monom]
    V = [[domain.one] + [domain.zero] * (len(old_basis) - 1)]
    G = []

    L = [(i, 0) for i in range(ngens)]  # (i, j) corresponds to x_i * S[j]
    L.sort(key=lambda k_l: O_to(_incr_k(S[k_l[1]], k_l[0])), reverse=True)
    t = L.pop()

    P = _identity_matrix(len(old_basis), domain)

    while True:
        s = len(S)
        v = _matrix_mul(M[t[0]], V[t[1]])
        _lambda = _matrix_mul(P, v)

        if all(_lambda[i] == 0 for i in range(s, len(old_basis))):
            # there is a linear combination of v by V
            lt = ring.term_new(_incr_k(S[t[1]], t[0]), domain.one)
            rest = ring.from_dict({S[i]: _lambda[i] for i in range(s)})

            g = (lt - rest).set_ring(ring_to)
            if g:
                G.append(g)
        else:
            # v is linearly independant from V
            P = _update(s, _lambda, P)
            S.append(_incr_k(S[t[1]], t[0]))
            V.append(v)

            L.extend([(i, s) for i in range(ngens)])
            L = list(set(L))
            L.sort(key=lambda k_l: O_to(_incr_k(S[k_l[1]], k_l[0])), reverse=True)

        L = [(k, l) for (k, l) in L if all(not g.LM.divides(_incr_k(S[l], k)) for g in G)]

        if not L:
            G = [g.monic() for g in G]
            return sorted(G, key=lambda g: O_to(g.LM), reverse=True)

        t = L.pop()


def _incr_k(m, k):
    return tuple(list(m[:k]) + [m[k] + 1] + list(m[k + 1:]))


def _identity_matrix(n, domain):
    M = [[domain.zero]*n for _ in range(n)]

    for i in range(n):
        M[i][i] = domain.one

    return M


def _matrix_mul(M, v):
    return [sum(row[i] * v[i] for i in range(len(v))) for row in M]


def _update(s, _lambda, P):
    """Update ``P`` such that for the updated `P'` `P' v = e_{s}`."""
    k = min(j for j in range(s, len(_lambda)) if _lambda[j] != 0)

    for r in range(len(_lambda)):
        if r != k:
            P[r] = [P[r][j] - (P[k][j] * _lambda[r]) / _lambda[k] for j in range(len(P[r]))]

    P[k] = [P[k][j] / _lambda[k] for j in range(len(P[k]))]
    P[k], P[s] = P[s], P[k]

    return P


def _representing_matrices(basis, G, ring):
    r"""
    Compute the matrices corresponding to the linear maps `m \mapsto
    x_i m` for all variables `x_i`.

    """
    domain = ring.domain
    u = ring.ngens-1

    def var(i):
        m, = ring.gens[i]
        return m

    def representing_matrix(m):
        M = [[domain.zero] * len(basis) for _ in range(len(basis))]

        for i, v in enumerate(basis):
            r = ring.term_new(m*v, domain.one).div(G)[1]

            for monom, coeff in r.items():
                j = basis.index(monom)
                M[j][i] = coeff

        return M

    return [representing_matrix(var(i)) for i in range(u + 1)]


def _basis(G, ring):
    r"""
    Computes a list of monomials which are not divisible by the leading
    monomials wrt to ``O`` of ``G``. These monomials are a basis of
    `K[X_1, \ldots, X_n]/(G)`.

    """
    order = ring.order

    leading_monomials = [g.LM for g in G]
    candidates = [ring.zero_monom]
    basis = []

    while candidates:
        t = candidates.pop()
        basis.append(t)

        new_candidates = [_incr_k(t, k) for k in range(ring.ngens)
                          if all(not lmg.divides(_incr_k(t, k))
                                 for lmg in leading_monomials)]
        candidates.extend(new_candidates)
        candidates.sort(key=order, reverse=True)

    basis = list(set(basis))

    return sorted(basis, key=order)
