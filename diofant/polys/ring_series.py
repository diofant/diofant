"""Power series manipulating functions acting on polys.ring.PolyElement()"""

import math

from mpmath.libmp.libintmath import giant_steps, ifac

from ..core import Rational
from ..core.compatibility import as_int
from ..domains import QQ
from .monomials import monomial_min, monomial_mul


def _invert_monoms(p1):
    """
    Compute ``x**n * p1(1/x)`` for ``p1`` univariate polynomial.

    Examples
    ========

    >>> R, x = ring('x', ZZ)
    >>> p = x**2 + 2*x + 3
    >>> _invert_monoms(p)
    3*x**2 + 2*x + 1

    See Also
    ========

    diofant.polys.densebasic.dup_reverse
    """
    terms = list(p1.items())
    terms.sort()
    deg = p1.degree()
    ring = p1.ring
    p = ring.zero
    cv = p1.listcoeffs()
    mv = p1.listmonoms()
    for i in range(len(mv)):
        p[(deg - mv[i][0],)] = cv[i]
    return p


def _giant_steps(target):
    """
    list of precision steps for the Newton's method
    """
    res = giant_steps(2, target)
    if res[0] != 2:
        res = [2] + res
    return res


def rs_trunc(p1, x, prec):
    """
    truncate the series in the ``x`` variable with precision ``prec``,
    that is modulo ``O(x**prec)``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = x**10 + x**5 + x + 1
    >>> rs_trunc(p, x, 12)
    x**10 + x**5 + x + 1
    >>> rs_trunc(p, x, 10)
    x**5 + x + 1
    """

    ring = p1.ring
    p = ring.zero
    i = ring.gens.index(x)
    for exp1 in p1:
        if exp1[i] >= prec:
            continue
        p[exp1] = p1[exp1]
    return p


def rs_mul(p1, p2, x, prec):
    """
    product of series modulo ``O(x**prec)``

    ``x`` is the series variable or its position in the generators.

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p1 = x**2 + 2*x + 1
    >>> p2 = x + 1
    >>> rs_mul(p1, p2, x, 3)
    3*x**2 + 3*x + 1
    """

    ring = p1.ring
    p = ring.zero
    if ring.__class__ != p2.ring.__class__ or ring != p2.ring:
        raise ValueError('p1 and p2 must have the same ring')
    iv = ring.gens.index(x)

    get = p.get
    items2 = list(p2.items())
    items2.sort(key=lambda e: e[0][iv])
    if ring.ngens == 1:
        for exp1, v1 in p1.items():
            for exp2, v2 in items2:
                exp = exp1[0] + exp2[0]
                if exp < prec:
                    exp = (exp, )
                    p[exp] = get(exp, 0) + v1*v2
                else:
                    break
    else:
        for exp1, v1 in p1.items():
            for exp2, v2 in items2:
                if exp1[iv] + exp2[iv] < prec:
                    exp = monomial_mul(exp1, exp2)
                    p[exp] = get(exp, 0) + v1*v2
                else:
                    break

    p.strip_zero()
    return p


def rs_square(p1, x, prec):
    """
    square modulo ``O(x**prec)``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = x**2 + 2*x + 1
    >>> rs_square(p, x, 3)
    6*x**2 + 4*x + 1
    """
    ring = p1.ring
    p = ring.zero
    iv = ring.gens.index(x)
    get = p.get
    items = list(p1.items())
    items.sort(key=lambda e: e[0][iv])
    for i in range(len(items)):
        exp1, v1 = items[i]
        for j in range(i):
            exp2, v2 = items[j]
            if exp1[iv] + exp2[iv] < prec:
                exp = monomial_mul(exp1, exp2)
                p[exp] = get(exp, 0) + v1*v2
            else:
                break
    p = p.imul_num(2)
    get = p.get
    for expv, v in p1.items():
        if 2*expv[iv] < prec:
            e2 = monomial_mul(expv, expv)
            p[e2] = get(e2, 0) + v**2
    p.strip_zero()
    return p


def rs_pow(p1, n, x, prec):
    """
    return ``p1**n`` modulo ``O(x**prec)``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = x + 1
    >>> rs_pow(p, 4, x, 3)
    6*x**2 + 4*x + 1
    """
    R = p1.ring
    p = R.zero
    if isinstance(n, Rational):  # pragma: no cover
        raise NotImplementedError

    n = as_int(n)
    if n == 0:
        if p1:
            return R(1)
        else:
            raise ValueError('0**0 is undefined')
    if n < 0:
        p1 = rs_pow(p1, -n, x, prec)
        return rs_series_inversion(p1, x, prec)
    if n == 1:
        return rs_trunc(p1, x, prec)
    if n == 2:
        return rs_square(p1, x, prec)
    if n == 3:
        p2 = rs_square(p1, x, prec)
        return rs_mul(p1, p2, x, prec)
    p = R(1)
    while 1:
        if n & 1:
            p = rs_mul(p1, p, x, prec)
            n -= 1
            if not n:
                break
        p1 = rs_square(p1, x, prec)
        n = n // 2
    return p


def _has_constant_term(p, x):
    """
    test if ``p`` has a constant term in ``x``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = x**2 + x + 1
    >>> _has_constant_term(p, x)
    True
    """
    ring = p.ring
    iv = ring.gens.index(x)
    zm = ring.zero_monom
    a = [0]*ring.ngens
    a[iv] = 1
    miv = tuple(a)
    for expv in p:
        if monomial_min(expv, miv) == zm:
            return True
    return False


def _series_inversion1(p, x, prec):
    """
    univariate series inversion ``1/p`` modulo ``O(x**prec)``

    The Newton method is used.

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = x + 1
    >>> _series_inversion1(p, x, 4)
    -x**3 + x**2 - x + 1
    """
    ring = p.ring
    zm = ring.zero_monom
    assert zm in p and not _has_constant_term(p - p[zm], x)
    if p[zm] != ring(1):
        # TODO add check that it is a unit
        p1 = ring(1)/p[zm]
    else:
        p1 = ring(1)
    for precx in _giant_steps(prec):
        tmp = p1.square()
        tmp = rs_mul(tmp, p, x, precx)
        p1 = 2*p1 - tmp
    return p1


def rs_series_inversion(p, x, prec):
    """
    multivariate series inversion ``1/p`` modulo ``O(x**prec)``

    Examples
    ========

    >>> R, x, y = ring('x, y', QQ)
    >>> rs_series_inversion(1 + x*y**2, x, 4)
    -x**3*y**6 + x**2*y**4 - x*y**2 + 1
    >>> rs_series_inversion(1 + x*y**2, y, 4)
    -x*y**2 + 1
    """
    ring = p.ring
    zm = ring.zero_monom
    ii = ring.gens.index(x)
    m = min(p, key=lambda k: k[ii])[ii]
    if m:  # pragma: no cover
        raise NotImplementedError('no constant term in series')
    if zm not in p:  # pragma: no cover
        raise NotImplementedError('no constant term in series')
    if _has_constant_term(p - p[zm], x):
        raise NotImplementedError('p - p[0] must not have a constant term in the series variables')
    return _series_inversion1(p, x, prec)


def rs_series_from_list(p, c, x, prec, concur=1):
    """
    series ``sum c[n]*p**n`` modulo ``O(x**prec)``

    reduce the number of multiplication summing concurrently
    ``ax = [1, p, p**2, .., p**(J - 1)]``
    ``s = sum(c[i]*ax[i] for i in range(r, (r + 1)*J))*p**((K - 1)*J)``
    with ``K >= (n + 1)/J``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = x**2 + x + 1
    >>> c = [1, 2, 3]
    >>> rs_series_from_list(p, c, x, 4)
    6*x**3 + 11*x**2 + 8*x + 6
    >>> rs_trunc(1 + 2*p + 3*p**2, x, 4)
    6*x**3 + 11*x**2 + 8*x + 6
    >>> pc = R.from_list(list(reversed(c)))
    >>> rs_trunc(pc.compose(x, p), x, 4)
    6*x**3 + 11*x**2 + 8*x + 6

    See Also
    ========

    diofant.polys.rings.PolyElement.compose
    """

    ring = p.ring
    n = len(c)
    if not concur:
        q = ring(1)
        s = c[0]*q
        for i in range(1, n):
            q = rs_mul(q, p, x, prec)
            s += c[i]*q
        return s
    J = int(math.sqrt(n) + 1)
    K, r = divmod(n, J)
    if r:
        K += 1
    ax = [ring(1)]
    b = 1
    q = ring(1)
    if len(p) < 20:
        for i in range(1, J):
            q = rs_mul(q, p, x, prec)
            ax.append(q)
    else:
        for i in range(1, J):
            if i % 2 == 0:
                q = rs_square(ax[i//2], x, prec)
            else:
                q = rs_mul(q, p, x, prec)
            ax.append(q)
    # optimize using rs_square
    pj = rs_mul(ax[-1], p, x, prec)
    b = ring(1)
    s = ring(0)
    for k in range(K - 1):
        r = J*k
        s1 = c[r]
        for j in range(1, J):
            s1 += c[r + j]*ax[j]
        s1 = rs_mul(s1, b, x, prec)
        s += s1
        b = rs_mul(b, pj, x, prec)
        if not b:
            break
    k = K - 1
    r = J*k
    assert r < n
    s1 = c[r]*ring(1)
    for j in range(1, J):
        if r + j >= n:
            break
        s1 += c[r + j]*ax[j]
    s1 = rs_mul(s1, b, x, prec)
    s += s1
    return s


def rs_integrate(self, x):
    """
    integrate ``p`` with respect to ``x``

    Examples
    ========

    >>> R, x, y = ring('x, y', QQ)
    >>> p = x + x**2*y**3
    >>> rs_integrate(p, x)
    1/3*x**3*y**3 + 1/2*x**2
    """
    ring = self.ring
    p1 = ring.zero
    n = ring.gens.index(x)
    mn = [0]*ring.ngens
    mn[n] = 1
    mn = tuple(mn)

    for expv in self:
        e = monomial_mul(expv, mn)
        p1[e] = self[expv]/(expv[n] + 1)
    return p1


def rs_log(p, x, prec):
    """
    logarithm of ``p`` modulo ``O(x**prec)``

    Notes
    =====

    truncation of ``integral dx p**-1*d p/dx`` is used.

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> rs_log(1 + x, x, 8)
    1/7*x**7 - 1/6*x**6 + 1/5*x**5 - 1/4*x**4 + 1/3*x**3 - 1/2*x**2 + x
    """
    ring = p.ring
    if p == 1:
        return ring.zero
    if _has_constant_term(p - 1, x):  # pragma: no cover
        raise NotImplementedError('p - 1 must not have a constant term in the series variables')
    dlog = p.diff(x)
    dlog = rs_mul(dlog, _series_inversion1(p, x, prec), x, prec - 1)
    return rs_integrate(dlog, x)


def _exp1(p, x, prec):
    """
    helper function for ``rs_exp``
    """
    ring = p.ring
    p1 = ring(1)
    for precx in _giant_steps(prec):
        pt = p - rs_log(p1, x, precx)
        tmp = rs_mul(pt, p1, x, precx)
        p1 += tmp
    return p1


def rs_exp(p, x, prec):
    """
    exponentiation of a series modulo ``O(x**prec)``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> rs_exp(x**2, x, 7)
    1/6*x**6 + 1/2*x**4 + x**2 + 1
    """
    ring = p.ring
    if _has_constant_term(p, x):  # pragma: no cover
        raise NotImplementedError
    if len(p) > 20:
        return _exp1(p, x, prec)
    one = ring(1)
    n = 1
    k = 1
    c = []
    for k in range(prec):
        c.append(one/n)
        k += 1
        n *= k

    r = rs_series_from_list(p, c, x, prec)
    return r


def rs_newton(p, x, prec):
    """
    compute the truncated Newton sum of the polynomial ``p``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = x**2 - 2
    >>> rs_newton(p, x, 5)
    8*x**4 + 4*x**2 + 2
    """
    deg = p.degree()
    p1 = _invert_monoms(p)
    p2 = rs_series_inversion(p1, x, prec)
    p3 = rs_mul(p1.diff(x), p2, x, prec)
    res = deg - p3*x
    return res


def rs_hadamard_exp(p1, inverse=False):
    """
    return ``sum f_i/i!*x**i`` from ``sum f_i*x**i``,
    where ``x`` is the first variable.

    If ``invers=True`` return ``sum f_i*i!*x**i``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> p = 1 + x + x**2 + x**3
    >>> rs_hadamard_exp(p)
    1/6*x**3 + 1/2*x**2 + x + 1
    """
    R = p1.ring
    if R.domain != QQ:  # pragma: no cover
        raise NotImplementedError
    p = R.zero
    if not inverse:
        for exp1, v1 in p1.items():
            p[exp1] = v1/int(ifac(exp1[0]))
    else:
        for exp1, v1 in p1.items():
            p[exp1] = v1*int(ifac(exp1[0]))
    return p


def rs_compose_add(p1, p2):
    """
    compute the composed sum ``prod(p2(x - beta) for beta root of p1)``

    Examples
    ========

    >>> R, x = ring('x', QQ)
    >>> f = x**2 - 2
    >>> g = x**2 - 3
    >>> rs_compose_add(f, g)
    x**4 - 10*x**2 + 1

    References
    ==========

    .. [1] A. Bostan, P. Flajolet, B. Salvy and E. Schost "Fast
           Computation with Two Algebraic Numbers", (2002) Research
           Report 4579, Institut National de Recherche en
           Informatique et en Automatique.
    """
    R = p1.ring
    x = R.gens[0]
    prec = p1.degree() * p2.degree() + 1
    np1 = rs_newton(p1, x, prec)
    np1e = rs_hadamard_exp(np1)
    np2 = rs_newton(p2, x, prec)
    np2e = rs_hadamard_exp(np2)
    np3e = rs_mul(np1e, np2e, x, prec)
    np3 = rs_hadamard_exp(np3e, True)
    np3a = (np3[(0,)] - np3)/x
    q = rs_integrate(np3a, x)
    q = rs_exp(q, x, prec)
    q = _invert_monoms(q)
    q = q.primitive()[1]
    dp = p1.degree() * p2.degree() - q.degree()
    # `dp` is the multiplicity of the zeroes of the resultant;
    # these zeroes are missed in this computation so they are put here.
    # if p1 and p2 are monic irreducible polynomials,
    # there are zeroes in the resultant
    # if and only if p1 = p2 ; in fact in that case p1 and p2 have a
    # root in common, so gcd(p1, p2) != 1; being p1 and p2 irreducible
    # this means p1 = p2
    if dp:
        q = q*x**dp
    return q
