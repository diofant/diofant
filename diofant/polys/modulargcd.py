import operator
import random

from ..core import Dummy
from ..ntheory import nextprime
from ..ntheory.modular import crt, integer_rational_reconstruction
from . import rings
from .polyerrors import ModularGCDFailed


def _swap(f, i):
    """Make the variable `x_i` the leading one in a multivariate polynomial `f`."""
    ring = f.ring
    fswap = ring.zero
    for monom, coeff in f.items():
        monomswap = (monom[i],) + monom[:i] + monom[i+1:]
        fswap[monomswap] = coeff
    return fswap


def _chinese_remainder_reconstruction(hp, hq, p, q):
    r"""
    Construct a polynomial `h_{pq}` in
    `\mathbb{Z}_{p q}[x_0, \ldots, x_{k-1}]` such that

    .. math ::

        h_{pq} = h_p \; \mathrm{mod} \, p

        h_{pq} = h_q \; \mathrm{mod} \, q

    for relatively prime integers `p` and `q` and polynomials
    `h_p` and `h_q` in `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]` and
    `\mathbb{Z}_q[x_0, \ldots, x_{k-1}]` respectively.

    The coefficients of the polynomial `h_{pq}` are computed with the
    Chinese Remainder Theorem. The symmetric representation in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`,
    `\mathbb{Z}_q[x_0, \ldots, x_{k-1}]` and
    `\mathbb{Z}_{p q}[x_0, \ldots, x_{k-1}]` is used.

    Parameters
    ==========

    hp : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    hq : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_q`
    p : Integer
        modulus of `h_p`, relatively prime to `q`
    q : Integer
        modulus of `h_q`, relatively prime to `p`

    Examples
    ========

    >>> _, x, y = ring('x y', ZZ)
    >>> p = 3
    >>> q = 5

    >>> hp = x**3*y - x**2 - 1
    >>> hq = -x**3*y - 2*x*y**2 + 2

    >>> hpq = _chinese_remainder_reconstruction(hp, hq, p, q)
    >>> hpq
    4*x**3*y + 5*x**2 + 3*x*y**2 + 2

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    >>> _, x, y, z = ring('x y z', ZZ)
    >>> p = 6
    >>> q = 5

    >>> hp = 3*x**4 - y**3*z + z
    >>> hq = -2*x**4 + z

    >>> hpq = _chinese_remainder_reconstruction(hp, hq, p, q)
    >>> hpq
    3*x**4 + 5*y**3*z + z

    >>> hpq.trunc_ground(p) == hp
    True
    >>> hpq.trunc_ground(q) == hq
    True

    """
    hpmonoms = set(hp)
    hqmonoms = set(hq)
    monoms = hpmonoms.intersection(hqmonoms)
    hpmonoms.difference_update(monoms)
    hqmonoms.difference_update(monoms)

    zero = hp.ring.domain.zero

    hpq = hp.ring.zero

    if isinstance(hp.ring.domain, rings.PolynomialRing):
        crt_ = _chinese_remainder_reconstruction
    else:
        def crt_(cp, cq, p, q):
            return crt([p, q], [cp, cq], symmetric=True)[0]

    for monom in monoms:
        hpq[monom] = crt_(hp[monom], hq[monom], p, q)
    for monom in hpmonoms:
        hpq[monom] = crt_(hp[monom], zero, p, q)
    for monom in hqmonoms:
        hpq[monom] = crt_(zero, hq[monom], p, q)

    return hpq


def _interpolate(evalpoints, hpeval, ring, i, p, ground=False):
    r"""
    Reconstruct a polynomial `h_p` in `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`
    from a list of evaluation points in `\mathbb{Z}_p` and a list of
    polynomials in
    `\mathbb{Z}_p[x_0, \ldots, x_{i-1}, x_{i+1}, \ldots, x_{k-1}]`, which
    are the images of `h_p` evaluated in the variable `x_i`.

    It is also possible to reconstruct a parameter of the ground domain,
    i.e. if `h_p` is a polynomial over `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`.
    In this case, one has to set ``ground=True``.

    Parameters
    ==========

    evalpoints : list of Integer objects
        list of evaluation points in `\mathbb{Z}_p`
    hpeval : list of PolyElement objects
        list of polynomials in (resp. over)
        `\mathbb{Z}_p[x_0, \ldots, x_{i-1}, x_{i+1}, \ldots, x_{k-1}]`,
        images of `h_p` evaluated in the variable `x_i`
    ring : PolynomialRing
        `h_p` will be an element of this ring
    i : Integer
        index of the variable which has to be reconstructed
    p : Integer
        prime number, modulus of `h_p`
    ground : Boolean
        indicates whether `x_i` is in the ground domain, default is
        ``False``

    Returns
    =======

    hp : PolyElement
        interpolated polynomial in (resp. over)
        `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]`

    """
    hp = ring.zero

    if ground:
        domain = ring.domain.domain
        y = ring.domain.gens[i]
    else:
        domain = ring.domain
        y = ring.gens[i]

    for a, hpa in zip(evalpoints, hpeval):
        numer = ring.one
        denom = domain.one
        for b in evalpoints:
            if b == a:
                continue

            numer *= y - b
            denom *= a - b

        denom = domain.invert(denom, p)
        coeff = numer*denom
        hp += hpa.set_ring(ring) * coeff

    return hp.trunc_ground(p)


def _modgcd_p(f, g, p, degbound, contbound):
    r"""
    Compute the GCD of two polynomials in
    `\mathbb{Z}_p[x0, \ldots, x{k-1}]`.

    The algorithm reduces the problem step by step by evaluating the
    polynomials `f` and `g` at `x_{k-1} = a` for suitable
    `a \in \mathbb{Z}_p` and then calls itself recursively to compute the GCD
    in `\mathbb{Z}_p[x_0, \ldots, x_{k-2}]`. If these recursive calls are
    successful for enough evaluation points, the GCD in `k` variables is
    interpolated, otherwise the algorithm returns ``None``. Every time a GCD
    or a content is computed, their degrees are compared with the bounds. If
    a degree greater then the bound is encountered, then the current call
    returns ``None`` and a new evaluation point has to be chosen. If at some
    point the degree is smaller, the correspondent bound is updated and the
    algorithm fails.

    Parameters
    ==========

    f : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    g : PolyElement
        multivariate integer polynomial with coefficients in `\mathbb{Z}_p`
    p : Integer
        prime number, modulus of `f` and `g`
    degbound : list of Integer objects
        ``degbound[i]`` is an upper bound for the degree of the GCD of `f`
        and `g` in the variable `x_i`
    contbound : list of Integer objects
        ``contbound[i]`` is an upper bound for the degree of the content of
        the GCD in `\mathbb{Z}_p[x_i][x_0, \ldots, x_{i-1}]`,
        ``contbound[0]`` is not used can therefore be chosen
        arbitrarily.

    Returns
    =======

    h : PolyElement
        GCD of the polynomials `f` and `g` or ``None``

    References
    ==========

    * :cite:`Monagan2000Brown`
    * :cite:`Brown1971gcd`

    """
    ring = f.ring
    k = ring.ngens

    domain = ring.domain
    pdomain = domain.finite_field(p)
    pring = ring.clone(domain=pdomain)

    f, g = map(operator.methodcaller('set_domain', pdomain), (f, g))

    if ring.is_univariate:
        h = pring.gcd(f, g)
        degh = h.degree()

        if degh > degbound[0]:
            return
        if degh < degbound[0]:
            degbound[0] = degh
            raise ModularGCDFailed

        return h.set_domain(domain)

    ypring = pring.eject(-1).domain

    degyf = f.degree(-1)
    degyg = g.degree(-1)

    contf, f = f.eject(-1).primitive()
    contg, g = g.eject(-1).primitive()

    f = f.inject()
    g = g.inject()

    conth = ypring.gcd(contf, contg)  # polynomial in Z_p[y]

    degcontf = contf.degree()
    degcontg = contg.degree()
    degconth = conth.degree()

    if degconth > contbound[k-1]:
        return
    if degconth < contbound[k-1]:
        contbound[k-1] = degconth
        raise ModularGCDFailed

    lcf = f.eject(-1).LC
    lcg = g.eject(-1).LC

    delta = ypring.gcd(lcf, lcg)  # polynomial in Z_p[y]

    evaltest = delta

    for i in range(k-1):
        evaltest *= ypring.gcd(_swap(f, i).eject(-1).LC,
                               _swap(g, i).eject(-1).LC)

    degdelta = delta.degree()

    N = min(degyf - degcontf, degyg - degcontg,
            degbound[k-1] - contbound[k-1] + degdelta) + 1

    if p < N:
        return

    n = 0
    d = 0
    evalpoints = []
    heval = []
    points = list(range(p))
    random.shuffle(points)

    while points:
        a = points.pop()

        if not evaltest.eval(0, a):
            continue

        deltaa = delta.eval(0, a)

        fa = f.eval(-1, a)
        ga = g.eval(-1, a)

        # polynomials in Z_p[x_0, ..., x_{k-2}]
        ha = _modgcd_p(fa.set_domain(domain), ga.set_domain(domain),
                       p, degbound, contbound)

        if ha is None:
            if d < n:
                d += 1
                continue
            else:
                return

        ha = ha.set_domain(pdomain)

        if ha.is_ground:
            return conth.set_ring(ring)

        ha *= deltaa

        evalpoints.append(a)
        heval.append(ha)
        n += 1

        if n == N:
            heval = [_.set_domain(domain) for _ in heval]
            h = _interpolate(evalpoints, heval, ring, -1, p)

            h = h.set_domain(pdomain)

            h = h.eject(-1).primitive()[1].inject().set_ring(ring) * conth.set_ring(ring)
            degyh = h.degree(-1)

            assert degyh <= degbound[k-1]
            if degyh < degbound[k-1]:
                degbound[k-1] = degyh
                raise ModularGCDFailed

            return h


def modgcd(f, g):
    r"""
    Compute the GCD of two polynomials in `\mathbb{Z}[x_0, \ldots, x_{k-1}]`
    using a modular algorithm.

    The algorithm computes the GCD of two multivariate integer polynomials
    `f` and `g` by calculating the GCD in
    `\mathbb{Z}_p[x_0, \ldots, x_{k-1}]` for suitable primes `p` and then
    reconstructing the coefficients with the Chinese Remainder Theorem. To
    compute the multivariate GCD over `\mathbb{Z}_p` the recursive
    subroutine ``_modgcd_p`` is used. To verify the result in
    `\mathbb{Z}[x_0, \ldots, x_{k-1}]`, trial division is done, but only for
    candidates which are very likely the desired GCD.

    Parameters
    ==========

    f : PolyElement
        multivariate integer polynomial
    g : PolyElement
        multivariate integer polynomial

    Returns
    =======

    h : PolyElement
        GCD of the polynomials `f` and `g`

    Examples
    ========

    >>> _, x, y = ring('x y', ZZ)

    >>> modgcd((x - y)*(x + y), (x + y)**2)
    x + y

    >>> _, x, y, z = ring('x y z', ZZ)

    >>> modgcd((x - y)*z**2, (x**2 + 1)*z)
    z

    References
    ==========

    * :cite:`Monagan2000Brown`
    * :cite:`Brown1971gcd`

    """
    assert f.ring == g.ring and f.ring.domain.is_IntegerRing

    ring = f.ring
    k = ring.ngens

    # divide out integer content
    cf, f = f.primitive()
    cg, g = g.primitive()
    ch = ring.domain.gcd(cf, cg)

    gamma = ring.domain.gcd(f.LC, g.LC)

    badprimes = ring.domain.one
    for i in range(k):
        badprimes *= ring.domain.gcd(_swap(f, i).LC, _swap(g, i).LC)

    degbound = [min(f.degree(x), g.degree(x)) for x in ring.gens]
    contbound = list(degbound)

    m = 1
    p = 1

    while True:
        p = nextprime(p)
        while badprimes % p == 0:
            p = nextprime(p)

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)

        try:
            # monic GCD of fp, gp in Z_p[x_0, ..., x_{k-2}, y]
            hp = _modgcd_p(fp, gp, p, degbound, contbound)
        except ModularGCDFailed:
            m = 1
            continue

        if hp is None:
            continue

        hp = (hp*gamma).trunc_ground(p)
        if m == 1:
            m = p
            hlastm = hp
            continue

        hm = _chinese_remainder_reconstruction(hp, hlastm, p, m)
        m *= p

        if not hm == hlastm:
            hlastm = hm
            continue

        h = hm.primitive()[1]
        _, frem = divmod(f, h)
        _, grem = divmod(g, h)
        if not frem and not grem:
            return h*ch


def _rational_function_reconstruction(c, p, m):
    r"""
    Reconstruct a rational function `\frac a b` in `\mathbb Z_p(t)` from

    .. math::

        c = \frac a b \; \mathrm{mod} \, m,

    where `c` and `m` are polynomials in `\mathbb Z_p[t]` and `m` has
    positive degree.

    The algorithm is based on the Euclidean Algorithm. In general, `m` is
    not irreducible, so it is possible that `b` is not invertible modulo
    `m`. In that case ``None`` is returned.

    Parameters
    ==========

    c : PolyElement
        univariate polynomial in `\mathbb Z[t]`
    p : Integer
        prime number
    m : PolyElement
        modulus, not necessarily irreducible

    Returns
    =======

    frac : FracElement
        either `\frac a b` in `\mathbb Z(t)` or ``None``

    References
    ==========

    * :cite:`Monagan2004algebraic`

    """
    ring = c.ring
    field = ring.field
    domain = ring.domain
    M = m.degree()
    N = M // 2
    D = M - N - 1

    pdomain = domain.finite_field(p)
    pring = ring.clone(domain=pdomain)
    c, m = map(operator.methodcaller('set_domain', pdomain), (c, m))

    r0, s0 = m, pring.zero
    r1, s1 = c, pring.one

    while r1.degree() > N:
        quo = divmod(r0, r1)[0]
        r0, r1 = r1, r0 - quo*r1
        s0, s1 = s1, s0 - quo*s1

    if s1.degree() > D or pring.gcd(s1, m) != 1:
        return

    lc = s1.LC
    if lc != 1:
        r1 /= lc
        s1 /= lc

    r1, s1 = map(operator.methodcaller('set_domain', domain), (r1, s1))

    return field(r1) / field(s1)


def _rational_reconstruction_func_coeffs(hm, p, m, ring, k):
    r"""
    Reconstruct every coefficient `c_h` of a polynomial `h` in
    `\mathbb Z_p(t_k)[t_1, \ldots, t_{k-1}][x, z]` from the corresponding
    coefficient `c_{h_m}` of a polynomial `h_m` in
    `\mathbb Z_p[t_1, \ldots, t_k][x, z] \cong \mathbb Z_p[t_k][t_1, \ldots, t_{k-1}][x, z]`
    such that

    .. math::

        c_{h_m} = c_h \; \mathrm{mod} \, m,

    where `m \in \mathbb Z_p[t]`.

    The reconstruction is based on the Euclidean Algorithm. In general, `m`
    is not irreducible, so it is possible that this fails for some
    coefficient. In that case ``None`` is returned.

    Parameters
    ==========

    hm : PolyElement
        polynomial in `\mathbb Z[t_1, \ldots, t_k][x, z]`
    p : Integer
        prime number, modulus of `\mathbb Z_p`
    m : PolyElement
        modulus, polynomial in `\mathbb Z[t]`, not necessarily irreducible
    ring : PolynomialRing
        `\mathbb Z(t_k)[t_1, \ldots, t_{k-1}][x, z]`, `h` will be an
        element of this ring
    k : Integer
        index of the parameter `t_k` which will be reconstructed

    Returns
    =======

    h : PolyElement
        reconstructed polynomial in
        `\mathbb Z(t_k)[t_1, \ldots, t_{k-1}][x, z]` or ``None``

    See also
    ========

    _rational_function_reconstruction

    """
    h = ring.zero

    for monom, coeff in hm.items():
        if k == 0:
            coeffh = _rational_function_reconstruction(coeff, p, m)

            if not coeffh:
                return

        else:
            coeffh = ring.domain.zero
            for mon, c in coeff.eject(k).items():
                ch = _rational_function_reconstruction(c, p, m)

                if not ch:
                    return

                coeffh[mon] = ch

        h[monom] = coeffh

    return h


def _gf_gcdex(f, g, p):
    r"""
    Extended Euclidean Algorithm for two univariate polynomials over
    `\mathbb Z_p`.

    Returns polynomials `s, t` and `h`, such that `h` is the GCD of `f` and
    `g` and `sf + tg = h \; \mathrm{mod} \, p`.

    """
    ring = f.ring
    dom = ring.domain
    domp = dom.finite_field(p)
    f, g = map(operator.methodcaller('set_domain', domp), (f, g))
    return tuple(map(operator.methodcaller('set_domain', dom), f.gcdex(g)))


def _trunc(f, minpoly, p):
    r"""
    Compute the reduced representation of a polynomial `f` in
    `\mathbb Z_p[z] / (\check m_{\alpha}(z))[x]`

    Parameters
    ==========

    f : PolyElement
        polynomial in `\mathbb Z[x, z]`
    minpoly : PolyElement
        polynomial `\check m_{\alpha} \in \mathbb Z[z]`, not necessarily
        irreducible
    p : Integer
        prime number, modulus of `\mathbb Z_p`

    Returns
    =======

    ftrunc : PolyElement
        polynomial in `\mathbb Z[x, z]`, reduced modulo
        `\check m_{\alpha}(z)` and `p`

    """
    ring = f.ring
    minpoly = minpoly.set_ring(ring)
    p_ = ring.ground_new(p)

    return f.trunc_ground(p).div([minpoly, p_])[1].trunc_ground(p)


def _euclidean_algorithm(f, g, minpoly, p):
    r"""
    Compute the monic GCD of two univariate polynomials in
    `\mathbb{Z}_p[z]/(\check m_{\alpha}(z))[x]` with the Euclidean
    Algorithm.

    In general, `\check m_{\alpha}(z)` is not irreducible, so it is possible
    that some leading coefficient is not invertible modulo
    `\check m_{\alpha}(z)`. In that case ``None`` is returned.

    Parameters
    ==========

    f, g : PolyElement
        polynomials in `\mathbb Z[x, z]`
    minpoly : PolyElement
        polynomial in `\mathbb Z[z]`, not necessarily irreducible
    p : Integer
        prime number, modulus of `\mathbb Z_p`

    Returns
    =======

    h : PolyElement
        GCD of `f` and `g` in `\mathbb Z[z, x]` or ``None``, coefficients
        are in `\left[ -\frac{p-1} 2, \frac{p-1} 2 \right]`

    """
    ring = f.ring
    domain = ring.domain

    f = _trunc(f, minpoly, p)
    g = _trunc(g, minpoly, p)

    while g:
        rem = f
        deg = g.degree(0)  # degree in x
        lcinv, _, gcd = _gf_gcdex(g.eject(-1).LC, minpoly, p)

        if not gcd == 1:
            return

        while True:
            degrem = rem.degree(0)  # degree in x
            if degrem < deg:
                break
            quo = (lcinv * rem.eject(-1).LC).set_ring(ring)
            m = ring.from_terms([((degrem - deg, 0), domain.one)])
            rem = _trunc(rem - g*m*quo, minpoly, p)

        f = g
        g = rem

    lcfinv = _gf_gcdex(f.eject(-1).LC, minpoly, p)[0].set_ring(ring)

    return _trunc(f * lcfinv, minpoly, p)


def trial_division(f, h, minpoly, p=None):
    r"""
    Check if `h` divides `f` in
    `\mathbb K[t_1, \ldots, t_k][z]/(m_{\alpha}(z))`, where `\mathbb K` is
    either `\mathbb Q` or `\mathbb Z_p`.

    This algorithm is based on pseudo division and does not use any
    fractions. By default `\mathbb K` is `\mathbb Q`, if a prime number `p`
    is given, `\mathbb Z_p` is chosen instead.

    Parameters
    ==========

    f, h : PolyElement
        polynomials in `\mathbb Z[t_1, \ldots, t_k][x, z]`
    minpoly : PolyElement
        polynomial `m_{\alpha}(z)` in `\mathbb Z[t_1, \ldots, t_k][z]`
    p : Integer or None
        if `p` is given, `\mathbb K` is set to `\mathbb Z_p` instead of
        `\mathbb Q`, default is ``None``

    Returns
    =======

    rem : PolyElement
        remainder of `\frac f h`

    References
    ==========

    * :cite:`vanHoeij2002modgcd`

    """
    ring = f.ring
    domain = ring.domain
    zxring = ring.clone(symbols=(ring.symbols[1], ring.symbols[0]))
    minpoly = minpoly.set_ring(ring)
    rem = f

    degrem = rem.degree()
    degh = h.degree()
    degm = minpoly.degree(1)

    lch = h.eject(-1).LC.set_ring(ring)
    lcm = minpoly.LC

    while rem and degrem >= degh:
        # polynomial in Z[t_1, ..., t_k][z]
        lcrem = rem.eject(-1).LC.set_ring(ring)
        rem = rem*lch - h*ring.from_terms([((degrem - degh, 0),
                                            domain.one)])*lcrem
        if p:
            rem = rem.trunc_ground(p)
        degrem = rem.degree(1)

        while rem and degrem >= degm:
            # polynomial in Z[t_1, ..., t_k][x]
            lcrem = rem.set_ring(zxring).eject(-1).LC.set_ring(ring)
            rem = rem*lcm - minpoly*ring.from_terms([((0, degrem - degm),
                                                     domain.one)])*lcrem
            if p:
                rem = rem.trunc_ground(p)
            degrem = rem.degree(1)

        degrem = rem.degree()

    return rem


def _evaluate_ground(f, i, a):
    r"""
    Evaluate a polynomial `f` at `a` in the `i`-th variable of the ground
    domain.

    """
    ring = f.ring.clone(domain=f.ring.domain.ring.drop(i))
    fa = ring.zero

    for monom, coeff in f.items():
        fa[monom] = coeff.eval(i, a)

    return fa


def _func_field_modgcd_p(f, g, minpoly, p):
    r"""
    Compute the GCD of two polynomials `f` and `g` in
    `\mathbb Z_p(t_1, \ldots, t_k)[z]/(\check m_\alpha(z))[x]`.

    The algorithm reduces the problem step by step by evaluating the
    polynomials `f` and `g` at `t_k = a` for suitable `a \in \mathbb Z_p`
    and then calls itself recursively to compute the GCD in
    `\mathbb Z_p(t_1, \ldots, t_{k-1})[z]/(\check m_\alpha(z))[x]`. If these
    recursive calls are successful, the GCD over `k` variables is
    interpolated, otherwise the algorithm returns ``None``. After
    interpolation, Rational Function Reconstruction is used to obtain the
    correct coefficients. If this fails, a new evaluation point has to be
    chosen, otherwise the desired polynomial is obtained by clearing
    denominators. The result is verified with a fraction free trial
    division.

    Parameters
    ==========

    f, g : PolyElement
        polynomials in `\mathbb Z[t_1, \ldots, t_k][x, z]`
    minpoly : PolyElement
        polynomial in `\mathbb Z[t_1, \ldots, t_k][z]`, not necessarily
        irreducible
    p : Integer
        prime number, modulus of `\mathbb Z_p`

    Returns
    =======

    h : PolyElement
        primitive associate in `\mathbb Z[t_1, \ldots, t_k][x, z]` of the
        GCD of the polynomials `f` and `g`  or ``None``, coefficients are
        in `\left[ -\frac{p-1} 2, \frac{p-1} 2 \right]`

    References
    ==========

    * :cite:`Monagan2004algebraic`

    """
    ring = f.ring
    domain = ring.domain  # Z[t_1, ..., t_k]

    if isinstance(domain, rings.PolynomialRing):
        k = domain.ngens
    else:
        return _euclidean_algorithm(f, g, minpoly, p)

    if k == 1:
        qdomain = domain.ring.field
    else:
        qdomain = domain.ring.eject(k - 1)
        qdomain = qdomain.clone(domain=qdomain.domain.ring.field)

    qring = ring.clone(domain=qdomain)  # = Z(t_k)[t_1, ..., t_{k-1}][x, z]

    n = 1
    d = 1

    # polynomial in Z_p[t_1, ..., t_k][z]
    gamma = f.eject(-1).LC * g.eject(-1).LC
    # polynomial in Z_p[t_1, ..., t_k]
    delta = minpoly.LC

    evalpoints = []
    heval = []
    LMlist = []
    points = list(range(p))
    random.shuffle(points)

    while points:
        a = points.pop()

        if k == 1:
            test = delta.eval(k-1, a) % p == 0
        else:
            test = delta.eval(k-1, a).trunc_ground(p) == 0

        if test:
            continue

        gammaa = _evaluate_ground(gamma, k-1, a)
        minpolya = _evaluate_ground(minpoly, k-1, a)

        if gammaa.div([minpolya, gammaa.ring(p)])[1] == 0:
            continue

        fa = _evaluate_ground(f, k-1, a)
        ga = _evaluate_ground(g, k-1, a)

        # polynomial in Z_p[x, t_1, ..., t_{k-1}, z]/(minpoly)
        ha = _func_field_modgcd_p(fa, ga, minpolya, p)

        if ha is None:
            if d < n:
                d += 1
                continue
            else:
                return

        if ha == 1:
            return ha

        LM = [ha.degree()] + [0]*(k-1)
        if k > 1:
            for monom, coeff in ha.items():
                if monom[0] == LM[0] and coeff.LM > tuple(LM[1:]):
                    LM[1:] = coeff.LM

        evalpoints_a = [a]
        heval_a = [ha]
        if k == 1:
            m = qring.domain.ring.one
        else:
            m = qring.domain.domain.ring.one

        t = m.ring.gens[0]

        for b, hb, LMhb in zip(evalpoints, heval, LMlist):
            if LMhb == LM:
                evalpoints_a.append(b)
                heval_a.append(hb)
                m *= (t - b)

        m = m.trunc_ground(p)
        evalpoints.append(a)
        heval.append(ha)
        LMlist.append(LM)
        n += 1

        # polynomial in Z_p[t_1, ..., t_k][x, z]
        h = _interpolate(evalpoints_a, heval_a, ring, k-1, p, ground=True)

        # polynomial in Z_p(t_k)[t_1, ..., t_{k-1}][x, z]
        h = _rational_reconstruction_func_coeffs(h, p, m, qring, k-1)

        if h is None:
            continue

        if k == 1:
            dom = qring.domain.field
            domp = dom.ring.domain.finite_field(p)
            new_dom = dom.clone(domain=domp)
            h = h.set_domain(new_dom)
            _, h = h.clear_denoms()
            h = h.set_domain(domain)
        else:
            dom = qring.domain.domain.field
            domp = dom.ring.domain.finite_field(p)
            new_dom = dom.clone(domain=domp)
            den = new_dom.ring.one

            for c in h.values():
                c = c.set_domain(new_dom)
                c, _ = c.clear_denoms()
                den = den.ring.lcm(den, c)

            den = new_dom((den, new_dom.ring.one))

            for key, coeff in h.items():
                coeff = coeff.set_domain(new_dom)
                h[key] = (coeff*den).set_domain(dom)

            h = ring.from_expr(h.ring.to_expr(h))

        if not trial_division(f, h, minpoly, p) and not trial_division(g, h, minpoly, p):
            return h


def _rational_reconstruction_int_coeffs(hm, m, ring):
    r"""
    Reconstruct every rational coefficient `c_h` of a polynomial `h` in
    `\mathbb Q[t_1, \ldots, t_k][x, z]` from the corresponding integer
    coefficient `c_{h_m}` of a polynomial `h_m` in
    `\mathbb Z[t_1, \ldots, t_k][x, z]` such that

    .. math::

        c_{h_m} = c_h \; \mathrm{mod} \, m,

    where `m \in \mathbb Z`.

    The reconstruction is based on the Euclidean Algorithm. In general,
    `m` is not a prime number, so it is possible that this fails for some
    coefficient. In that case ``None`` is returned.

    Parameters
    ==========

    hm : PolyElement
        polynomial in `\mathbb Z[t_1, \ldots, t_k][x, z]`
    m : Integer
        modulus, not necessarily prime
    ring : PolynomialRing
        `\mathbb Q[t_1, \ldots, t_k][x, z]`, `h` will be an element of this
        ring

    Returns
    =======

    h : PolyElement
        reconstructed polynomial in `\mathbb Q[t_1, \ldots, t_k][x, z]` or
        ``None``

    See also
    ========

    diofant.ntheory.modular.integer_rational_reconstruction

    """
    h = ring.zero

    if isinstance(ring.domain, rings.PolynomialRing):
        reconstruction = _rational_reconstruction_int_coeffs
        domain = ring.domain.ring
    else:
        reconstruction = integer_rational_reconstruction
        domain = hm.ring.domain

    for monom, coeff in hm.items():
        coeffh = reconstruction(coeff, m, domain)

        if not coeffh:
            return

        h[monom] = coeffh

    return h


def _func_field_modgcd_m(f, g, minpoly):
    r"""
    Compute the GCD of two polynomials in
    `\mathbb Q(t_1, \ldots, t_k)[z]/(m_{\alpha}(z))[x]` using a modular
    algorithm.

    The algorithm computes the GCD of two polynomials `f` and `g` by
    calculating the GCD in
    `\mathbb Z_p(t_1, \ldots, t_k)[z] / (\check m_{\alpha}(z))[x]` for
    suitable primes `p` and the primitive associate `\check m_{\alpha}(z)`
    of `m_{\alpha}(z)`. Then the coefficients are reconstructed with the
    Chinese Remainder Theorem and Rational Reconstruction. To compute the
    GCD over `\mathbb Z_p(t_1, \ldots, t_k)[z] / (\check m_{\alpha})[x]`,
    the recursive subroutine ``_func_field_modgcd_p`` is used. To verify the
    result in `\mathbb Q(t_1, \ldots, t_k)[z] / (m_{\alpha}(z))[x]`, a
    fraction free trial division is used.

    Parameters
    ==========

    f, g : PolyElement
        polynomials in `\mathbb Z[t_1, \ldots, t_k][x, z]`
    minpoly : PolyElement
        irreducible polynomial in `\mathbb Z[t_1, \ldots, t_k][z]`

    Returns
    =======

    h : PolyElement
        the primitive associate in `\mathbb Z[t_1, \ldots, t_k][x, z]` of
        the GCD of `f` and `g`

    Examples
    ========

    >>> _, x, z = ring('x z', ZZ)
    >>> minpoly = (z**2 - 2).drop(0)

    >>> f = x**2 + 2*x*z + 2
    >>> g = x + z
    >>> _func_field_modgcd_m(f, g, minpoly)
    x + z

    >>> D, t = ring('t', ZZ)
    >>> _, x, z = ring('x z', D)
    >>> minpoly = (z**2-3).drop(0)

    >>> f = x**2 + (t + 1)*x*z + 3*t
    >>> g = x*z + 3*t
    >>> _func_field_modgcd_m(f, g, minpoly)
    x + t*z

    References
    ==========

    * :cite:`Monagan2004algebraic`

    See also
    ========

    _func_field_modgcd_p

    """
    ring = f.ring
    domain = ring.domain

    if isinstance(domain, rings.PolynomialRing):
        k = domain.ngens
        QQdomain = domain.ring.clone(domain=domain.domain.field)
        QQring = ring.clone(domain=QQdomain)
    else:
        k = 0
        QQring = ring.clone(domain=domain.field)

    cf, f = f.primitive()
    cg, g = g.primitive()

    # polynomial in Z[t_1, ..., t_k][z]
    gamma = f.eject(-1).LC * g.eject(-1).LC
    # polynomial in Z[t_1, ..., t_k]
    delta = minpoly.LC
    assert k > 0 or delta == 1

    p = 1
    primes = []
    hplist = []
    LMlist = []

    while True:
        p = nextprime(p)

        if gamma.trunc_ground(p) == 0:
            continue

        fp = f.trunc_ground(p)
        gp = g.trunc_ground(p)
        minpolyp = minpoly.trunc_ground(p)

        hp = _func_field_modgcd_p(fp, gp, minpolyp, p)

        if hp is None:
            continue

        if hp == 1:
            return ring.one

        LM = [hp.degree()] + [0]*k
        if k > 0:
            for monom, coeff in hp.items():
                if monom[0] == LM[0] and coeff.LM > tuple(LM[1:]):
                    LM[1:] = coeff.LM

        hm = hp
        m = p

        for q, hq, LMhq in zip(primes, hplist, LMlist):
            if LMhq == LM:
                hm = _chinese_remainder_reconstruction(hq, hm, q, m)
                m *= q

        primes.append(p)
        hplist.append(hp)
        LMlist.append(LM)

        hm = _rational_reconstruction_int_coeffs(hm, m, QQring)

        if hm is None:
            continue

        if k == 0:
            h = hm.clear_denoms()[1]
        else:
            den = domain.domain.one
            for coeff in hm.values():
                den = domain.domain.lcm(den, coeff.clear_denoms()[0])
            h = hm*den

        # convert back to Z[t_1, ..., t_k][x, z] from Q[t_1, ..., t_k][x, z]
        h = h.set_ring(ring)
        h = h.primitive()[1]

        if not (trial_division(f*cf, h, minpoly) or
                trial_division(g*cg, h, minpoly)):
            return h


def _to_ZZ_poly(f, ring):
    r"""
    Compute an associate of a polynomial
    `f \in \mathbb Q(\alpha)[x_0, \ldots, x_{n-1}]` in
    `\mathbb Z[x_1, \ldots, x_{n-1}][z] / (\check m_{\alpha}(z))[x_0]`,
    where `\check m_{\alpha}(z) \in \mathbb Z[z]` is the primitive associate
    of the minimal polynomial `m_{\alpha}(z)` of `\alpha` over
    `\mathbb Q`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `\mathbb Q(\alpha)[x_0, \ldots, x_{n-1}]`
    ring : PolynomialRing
        `\mathbb Z[x_1, \ldots, x_{n-1}][x_0, z]`

    Returns
    =======

    f_ : PolyElement
        associate of `f` in
        `\mathbb Z[x_1, \ldots, x_{n-1}][x_0, z]`

    """
    f_ = ring.zero

    if isinstance(ring.domain, rings.PolynomialRing):
        domain = ring.domain.domain
    else:
        domain = ring.domain

    den = domain.one

    for coeff in f.values():
        for c in coeff.rep.all_coeffs():
            if c:
                den = domain.lcm(den, c.denominator)

    for monom, coeff in f.items():
        coeff = coeff.rep.all_coeffs()
        m = ring.domain.one
        if isinstance(ring.domain, rings.PolynomialRing):
            m *= ring.domain.from_terms([(monom[1:], domain.one)])
        n = len(coeff)

        for i in range(n):
            if coeff[n - i - 1]:
                c = domain(coeff[n - i - 1] * den) * m

                if (monom[0], n-i-1) not in f_:
                    f_[(monom[0], n-i-1)] = c
                else:
                    f_[(monom[0], n-i-1)] += c

    return f_


def _to_ANP_poly(f, ring):
    r"""
    Convert a polynomial
    `f \in \mathbb Z[x_1, \ldots, x_{n-1}][z]/(\check m_{\alpha}(z))[x_0]`
    to a polynomial in `\mathbb Q(\alpha)[x_0, \ldots, x_{n-1}]`,
    where `\check m_{\alpha}(z) \in \mathbb Z[z]` is the primitive associate
    of the minimal polynomial `m_{\alpha}(z)` of `\alpha` over
    `\mathbb Q`.

    Parameters
    ==========

    f : PolyElement
        polynomial in `\mathbb Z[x_1, \ldots, x_{n-1}][x_0, z]`
    ring : PolynomialRing
        `\mathbb Q(\alpha)[x_0, \ldots, x_{n-1}]`

    Returns
    =======

    f_ : PolyElement
        polynomial in `\mathbb Q(\alpha)[x_0, \ldots, x_{n-1}]`

    """
    domain = ring.domain
    f_ = ring.zero

    if isinstance(f.ring.domain, rings.PolynomialRing):
        for monom, coeff in f.items():
            for mon, coef in coeff.items():
                m = (monom[0],) + mon
                c = domain([0]*monom[1] + [domain.domain(coef)])

                if m not in f_:
                    f_[m] = c
                else:
                    f_[m] += c

    else:
        for monom, coeff in f.items():
            m = monom[0],
            c = domain([0]*monom[1] + [domain.domain(coeff)])

            if m not in f_:
                f_[m] = c
            else:
                f_[m] += c

    return f_


def _minpoly_from_dense(minpoly, ring):
    r"""
    Change representation of the minimal polynomial from ``Poly`` to
    ``PolyElement`` for a given ring.

    """
    minpoly_ = ring.zero

    for monom, coeff in minpoly.items():
        minpoly_[monom] = ring.domain(coeff)

    return minpoly_


def _primitive_in_x0(f):
    r"""
    Compute the content in `x_0` and the primitive part of a polynomial `f`
    in
    `\mathbb Q(\alpha)[x_0, x_1, \ldots, x_{n-1}] \cong \mathbb Q(\alpha)[x_1, \ldots, x_{n-1}][x_0]`.

    """
    fring = f.ring
    f_ = f.eject(*range(1, fring.ngens))
    dom = f_.ring.domain.ring
    cont = dom.zero

    for coeff in f_.values():
        cont = func_field_modgcd(cont, coeff)
        if cont == 1:
            return cont, f

    return cont, f//cont.set_ring(fring)


def func_field_modgcd(f, g):
    r"""
    Compute the GCD of two polynomials `f` and `g` in
    `\mathbb Q(\alpha)[x_0, \ldots, x_{n-1}]` using a modular algorithm.

    The algorithm first computes the primitive associate
    `\check m_{\alpha}(z)` of the minimal polynomial `m_{\alpha}` in
    `\mathbb{Z}[z]` and the primitive associates of `f` and `g` in
    `\mathbb{Z}[x_1, \ldots, x_{n-1}][z]/(\check m_{\alpha})[x_0]`. Then it
    computes the GCD in
    `\mathbb Q(x_1, \ldots, x_{n-1})[z]/(m_{\alpha}(z))[x_0]`.
    This is done by calculating the GCD in
    `\mathbb{Z}_p(x_1, \ldots, x_{n-1})[z]/(\check m_{\alpha}(z))[x_0]` for
    suitable primes `p` and then reconstructing the coefficients with the
    Chinese Remainder Theorem and Rational Reconstuction. The GCD over
    `\mathbb{Z}_p(x_1, \ldots, x_{n-1})[z]/(\check m_{\alpha}(z))[x_0]` is
    computed with a recursive subroutine, which evaluates the polynomials at
    `x_{n-1} = a` for suitable evaluation points `a \in \mathbb Z_p` and
    then calls itself recursively until the ground domain does no longer
    contain any parameters. For
    `\mathbb{Z}_p[z]/(\check m_{\alpha}(z))[x_0]` the Euclidean Algorithm is
    used. The results of those recursive calls are then interpolated and
    Rational Function Reconstruction is used to obtain the correct
    coefficients. The results, both in
    `\mathbb Q(x_1, \ldots, x_{n-1})[z]/(m_{\alpha}(z))[x_0]` and
    `\mathbb{Z}_p(x_1, \ldots, x_{n-1})[z]/(\check m_{\alpha}(z))[x_0]`, are
    verified by a fraction free trial division.

    Apart from the above GCD computation some GCDs in
    `\mathbb Q(\alpha)[x_1, \ldots, x_{n-1}]` have to be calculated,
    because treating the polynomials as univariate ones can result in
    a spurious content of the GCD. For this ``func_field_modgcd`` is
    called recursively.

    Parameters
    ==========

    f, g : PolyElement
        polynomials in `\mathbb Q(\alpha)[x_0, \ldots, x_{n-1}]`

    Returns
    =======

    h : PolyElement
        monic GCD of the polynomials `f` and `g`

    Examples
    ========

    >>> A = QQ.algebraic_field(sqrt(2))
    >>> _, x = ring('x', A)

    >>> func_field_modgcd(x**2 - 2, x + sqrt(2))
    x + sqrt(2)

    >>> _, x, y = ring('x y', A)

    >>> func_field_modgcd((x + sqrt(2)*y)**2, x + sqrt(2)*y)
    x + sqrt(2)*y

    >>> func_field_modgcd(x + sqrt(2)*y, x + y)
    1

    References
    ==========

    * :cite:`Monagan2004algebraic`

    """
    ring = f.ring
    domain = ring.domain
    n = ring.ngens

    assert ring == g.ring and domain.is_AlgebraicField

    if not f:
        return ring._gcd_zero(g)

    z = Dummy('z')

    ZZring = ring.clone(symbols=ring.symbols + (z,), domain=domain.domain.ring)

    if n == 1:
        f_ = _to_ZZ_poly(f, ZZring)
        g_ = _to_ZZ_poly(g, ZZring)
        minpoly = ZZring.drop(0).from_dict(domain.mod)

        h = _func_field_modgcd_m(f_, g_, minpoly)
        h = _to_ANP_poly(h, ring)

    else:
        # contx0f in Q(a)[x_1, ..., x_{n-1}], f in Q(a)[x_0, ..., x_{n-1}]
        contx0f, f = _primitive_in_x0(f)
        contx0g, g = _primitive_in_x0(g)
        contx0h = func_field_modgcd(contx0f, contx0g)

        ZZring_ = ZZring.eject(*range(1, n))

        f_ = _to_ZZ_poly(f, ZZring_)
        g_ = _to_ZZ_poly(g, ZZring_)
        minpoly = _minpoly_from_dense(domain.mod, ZZring_.drop(0))

        h = _func_field_modgcd_m(f_, g_, minpoly)
        h = _to_ANP_poly(h, ring)

        contx0h_, h = _primitive_in_x0(h)
        h *= contx0h.set_ring(ring)
        f *= contx0f.set_ring(ring)
        g *= contx0g.set_ring(ring)

    return h.quo_ground(h.LC)
