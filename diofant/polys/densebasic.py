"""Basic tools for dense recursive polynomials in ``K[x]`` or ``K[X]``. """

import math
import random

from ..core import oo
from .monomials import monomial_div, monomial_min
from .orderings import monomial_key


def dmp_LC(f, K):
    """
    Return leading coefficient of ``f``.

    Examples
    ========

    >>> dmp_LC([], ZZ)
    0
    >>> dmp_LC([ZZ(1), ZZ(2), ZZ(3)], ZZ)
    1
    """
    if not f:
        return K.zero
    else:
        return f[0]


def dmp_TC(f, K):
    """
    Return trailing coefficient of ``f``.

    Examples
    ========

    >>> dmp_TC([], ZZ)
    0
    >>> dmp_TC([ZZ(1), ZZ(2), ZZ(3)], ZZ)
    3
    """
    if not f:
        return K.zero
    else:
        return f[-1]


def dmp_ground_LC(f, u, K):
    """
    Return the ground leading coefficient.

    Examples
    ========

    >>> f = dmp_normal([[[1], [2, 3]]], 2, ZZ)

    >>> dmp_ground_LC(f, 2, ZZ)
    1
    """
    while u:
        f = dmp_LC(f, K)
        u -= 1

    return dmp_LC(f, K)


def dmp_ground_TC(f, u, K):
    """
    Return the ground trailing coefficient.

    Examples
    ========

    >>> f = dmp_normal([[[1], [2, 3]]], 2, ZZ)

    >>> dmp_ground_TC(f, 2, ZZ)
    3
    """
    while u:
        f = dmp_TC(f, K)
        u -= 1

    return dmp_TC(f, K)


def dmp_degree(f, u):
    """
    Return the leading degree of ``f`` in ``x_0`` in ``K[X]``.

    Note that the degree of 0 is negative infinity (the Diofant object -oo).

    Examples
    ========

    >>> dmp_degree([[[]]], 2)
    -oo

    >>> f = dmp_normal([[2], [1, 2, 3]], 1, ZZ)

    >>> dmp_degree(f, 1)
    1
    """
    return -oo if dmp_zero_p(f, u) else len(f) - 1


def dmp_degree_in(f, j, u):
    """
    Return the leading degree of ``f`` in ``x_j`` in ``K[X]``.

    Examples
    ========

    >>> f = dmp_normal([[2], [1, 2, 3]], 1, ZZ)

    >>> dmp_degree_in(f, 0, 1)
    1
    >>> dmp_degree_in(f, 1, 1)
    2
    """
    if not j:
        return dmp_degree(f, u)
    if j < 0 or j > u:
        raise IndexError("0 <= j <= %s expected, got %s" % (u, j))

    def degree_in(g, v, i, j):
        if i == j:
            return dmp_degree(g, v)

        v, i = v - 1, i + 1

        return max(degree_in(c, v, i, j) for c in g)

    return degree_in(f, u, 0, j)


def dmp_degree_list(f, u):
    """
    Return a list of degrees of ``f`` in ``K[X]``.

    Examples
    ========

    >>> f = dmp_normal([[1], [1, 2, 3]], 1, ZZ)
    >>> dmp_degree_list(f, 1)
    (1, 2)
    """
    degs = [-oo]*(u + 1)

    def degree_list(g, v, i, degs):
        degs[i] = max(degs[i], dmp_degree(g, v))

        if v > 0:
            v, i = v - 1, i + 1

            for c in g:
                degree_list(c, v, i, degs)

    degree_list(f, u, 0, degs)
    return tuple(degs)


def dmp_strip(f, u):
    """
    Remove leading zeros from ``f`` in ``K[X]``.

    Examples
    ========

    >>> dmp_strip([[], [0, 1, 2], [1]], 1)
    [[0, 1, 2], [1]]
    """
    if not u:
        for i, c in enumerate(f):
            if c:
                return f[i:]
        else:
            return dmp_zero(u)

    v = u - 1

    for i, c in enumerate(f):
        if not dmp_zero_p(c, v):
            return f[i:]
    else:
        return dmp_zero(u)


def dmp_validate(f, K=None):
    """
    Return the number of levels in ``f`` and recursively strip it.

    Examples
    ========

    >>> dmp_validate([[], [0, 1, 2], [1]])
    ([[1, 2], [1]], 1)

    >>> dmp_validate([[1], 1])
    Traceback (most recent call last):
    ...
    ValueError: invalid data structure for a multivariate polynomial
    """
    def validate(f, g, i, K):
        if type(g) is not list:
            if K is not None and not isinstance(g, K.dtype):
                raise TypeError("%s in %s in not of type %s" % (g, f, K.dtype))

            return {i - 1}
        elif not g:
            return {i}
        else:
            levels = set()

            for c in g:
                levels |= validate(f, c, i + 1, K)

            return levels

    levels = validate(f, f, 0, K)

    u = levels.pop()

    def strip(g, v):
        if not v:
            return dmp_strip(g, 0)

        w = v - 1

        return dmp_strip([strip(c, w) for c in g], v)

    if not levels:
        return strip(f, u), u
    else:
        raise ValueError(
            "invalid data structure for a multivariate polynomial")


def dup_reverse(f):
    """
    Compute ``x**n * f(1/x)``, i.e.: reverse ``f`` in ``K[x]``.

    Examples
    ========

    >>> f = dmp_normal([1, 2, 3, 0], 0, ZZ)
    >>> dup_reverse(f)
    [3, 2, 1]
    """
    return dmp_strip(list(reversed(f)), 0)


def dmp_copy(f, u):
    """
    Create a new copy of a polynomial ``f`` in ``K[X]``.

    Examples
    ========

    >>> f = dmp_normal([[1], [1, 2]], 1, ZZ)
    >>> dmp_copy(f, 1)
    [[1], [1, 2]]
    """
    if not u:
        return list(f)

    v = u - 1
    return [dmp_copy(c, v) for c in f]


def dmp_to_tuple(f, u):
    """
    Convert `f` into a nested tuple of tuples.

    This is needed for hashing.  This is similar to dmp_copy().

    Examples
    ========

    >>> f = dmp_normal([1, 2, 3, 0], 0, ZZ)
    >>> dmp_to_tuple(f, 0)
    (1, 2, 3, 0)

    >>> f = dmp_normal([[1], [1, 2]], 1, ZZ)
    >>> dmp_to_tuple(f, 1)
    ((1,), (1, 2))
    """
    if not u:
        return tuple(f)

    v = u - 1
    return tuple(dmp_to_tuple(c, v) for c in f)


def dmp_normal(f, u, K):
    """
    Normalize a multivariate polynomial in the given domain.

    Examples
    ========

    >>> dmp_normal([[], [0, 1.5, 2]], 1, ZZ)
    [[1, 2]]
    """
    if not u:
        r = [K.normal(c) for c in f]
    else:
        v = u - 1
        r = [dmp_normal(c, v, K) for c in f]

    return dmp_strip(r, u)


def dmp_convert(f, u, K0, K1):
    """
    Convert the ground domain of ``f`` from ``K0`` to ``K1``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> dmp_convert([[R(1)], [R(2)]], 1, R, ZZ)
    [[1], [2]]
    >>> dmp_convert([[ZZ(1)], [ZZ(2)]], 1, ZZ, R)
    [[1], [2]]
    """
    if K0 is not None and K0 == K1:
        return f

    if not u:
        r = [K1.convert(c, K0) for c in f]
    else:
        v = u - 1
        r = [dmp_convert(c, v, K0, K1) for c in f]

    return dmp_strip(r, u)


def dmp_ground_nth(f, N, u, K):
    """
    Return the ground ``n``-th coefficient of ``f`` in ``K[x]``.

    Examples
    ========

    >>> f = dmp_normal([[1], [2, 3]], 1, ZZ)
    >>> dmp_ground_nth(f, (0, 1), 1, ZZ)
    2
    """
    v = u

    for n in N:
        if n < 0:
            raise IndexError("`n` must be non-negative, got %i" % n)
        elif n >= len(f):
            return K.zero
        else:
            d = dmp_degree(f, v)
            if d == -oo:
                d = -1
            f, v = f[d - n], v - 1

    return f


def dmp_zero_p(f, u):
    """
    Return ``True`` if ``f`` is zero in ``K[X]``.

    Examples
    ========

    >>> dmp_zero_p([[[[[]]]]], 4)
    True
    >>> dmp_zero_p([[[[[1]]]]], 4)
    False
    """
    while u:
        if len(f) != 1:
            return False

        f = f[0]
        u -= 1

    return not f


def dmp_zero(u):
    """
    Return a multivariate zero.

    Examples
    ========

    >>> dmp_zero(4)
    [[[[[]]]]]
    """
    r = []

    for i in range(u):
        r = [r]

    return r


def dmp_one_p(f, u, K):
    """
    Return ``True`` if ``f`` is one in ``K[X]``.

    Examples
    ========

    >>> dmp_one_p([[[ZZ(1)]]], 2, ZZ)
    True
    """
    return dmp_ground_p(f, K.one, u)


def dmp_one(u, K):
    """
    Return a multivariate one over ``K``.

    Examples
    ========

    >>> dmp_one(2, ZZ)
    [[[1]]]
    """
    return dmp_ground(K.one, u)


def dmp_ground_p(f, c, u):
    """
    Return True if ``f`` is constant in ``K[X]``.

    Examples
    ========

    >>> dmp_ground_p([[[3]]], 3, 2)
    True
    >>> dmp_ground_p([[[4]]], None, 2)
    True
    """
    if c is not None and not c:
        return dmp_zero_p(f, u)

    while u:
        if len(f) != 1:
            return False
        f = f[0]
        u -= 1

    if c is None:
        return len(f) <= 1
    else:
        return f == [c]


def dmp_ground(c, u):
    """
    Return a multivariate constant.

    Examples
    ========

    >>> dmp_ground(3, 5)
    [[[[[[3]]]]]]
    >>> dmp_ground(1, -1)
    1
    """
    if not c:
        return dmp_zero(u)

    for i in range(u + 1):
        c = [c]

    return c


def dmp_zeros(n, u, K):
    """
    Return a list of multivariate zeros.

    Examples
    ========

    >>> dmp_zeros(3, 2, ZZ)
    [[[[]]], [[[]]], [[[]]]]
    >>> dmp_zeros(3, -1, ZZ)
    [0, 0, 0]
    """
    if not n:
        return []

    if u < 0:
        return [K.zero]*n
    else:
        return [ dmp_zero(u) for i in range(n) ]


def dup_from_dict(f, K):
    """
    Create a ``K[x]`` polynomial from a ``dict``.

    Examples
    ========

    >>> dup_from_dict({(0,): ZZ(7), (2,): ZZ(5), (4,): ZZ(1)}, ZZ)
    [1, 0, 5, 0, 7]
    >>> dup_from_dict({}, ZZ)
    []
    """
    if not f:
        return []

    n, h = max(f), []

    if type(n) is int:
        for k in range(n, -1, -1):
            h.append(f.get(k, K.zero))
    else:
        n, = n

        for k in range(n, -1, -1):
            h.append(f.get((k,), K.zero))

    return dmp_strip(h, 0)


def dmp_from_dict(f, u, K):
    """
    Create a ``K[X]`` polynomial from a ``dict``.

    Examples
    ========

    >>> dmp_from_dict({(0, 0): ZZ(3), (0, 1): ZZ(2), (2, 1): ZZ(1)}, 1, ZZ)
    [[1, 0], [], [2, 3]]
    >>> dmp_from_dict({}, 0, ZZ)
    []
    """
    if not u:
        return dup_from_dict(f, K)
    if not f:
        return dmp_zero(u)

    coeffs = {}

    for monom, coeff in f.items():
        head, tail = monom[0], monom[1:]

        if head in coeffs:
            coeffs[head][tail] = coeff
        else:
            coeffs[head] = { tail: coeff }

    n, v, h = max(coeffs), u - 1, []

    for k in range(n, -1, -1):
        coeff = coeffs.get(k)

        if coeff is not None:
            h.append(dmp_from_dict(coeff, v, K))
        else:
            h.append(dmp_zero(v))

    return dmp_strip(h, u)


def dmp_to_dict(f, u, K=None, zero=False):
    """
    Convert a ``K[X]`` polynomial to a ``dict````.

    Examples
    ========

    >>> dmp_to_dict([[1, 0], [], [2, 3]], 1)
    {(0, 0): 3, (0, 1): 2, (2, 1): 1}
    """
    if dmp_zero_p(f, u) and zero:
        return {(0,)*(u + 1): K.zero}

    n, v, result = dmp_degree(f, u), u - 1, {}

    if n == -oo:
        n = -1

    for k in range(n + 1):
        if f[n - k]:
            if u:
                h = dmp_to_dict(f[n - k], v)
                for exp, coeff in h.items():
                    result[(k,) + exp] = coeff
            else:
                result[(k,)] = f[n - k]

    return result


def dmp_swap(f, i, j, u, K):
    """
    Transform ``K[..x_i..x_j..]`` to ``K[..x_j..x_i..]``.

    Examples
    ========

    >>> f = dmp_normal([[[2], [1, 0]], []], 2, ZZ)

    >>> dmp_swap(f, 0, 1, 2, ZZ)
    [[[2], []], [[1, 0], []]]
    >>> dmp_swap(f, 1, 2, 2, ZZ)
    [[[1], [2, 0]], [[]]]
    >>> dmp_swap(f, 0, 2, 2, ZZ)
    [[[1, 0]], [[2, 0], []]]
    """
    if i < 0 or j < 0 or i > u or j > u:
        raise IndexError("0 <= i < j <= %s expected" % u)
    elif i == j:
        return f

    F, H = dmp_to_dict(f, u), {}

    for exp, coeff in F.items():
        H[exp[:i] + (exp[j],) +
          exp[i + 1:j] +
          (exp[i],) + exp[j + 1:]] = coeff

    return dmp_from_dict(H, u, K)


def dmp_permute(f, P, u, K):
    """
    Return a polynomial in ``K[x_{P(1)},..,x_{P(n)}]``.

    Examples
    ========

    >>> f = dmp_normal([[[2], [1, 0]], []], 2, ZZ)

    >>> dmp_permute(f, [1, 0, 2], 2, ZZ)
    [[[2], []], [[1, 0], []]]
    >>> dmp_permute(f, [1, 2, 0], 2, ZZ)
    [[[1], []], [[2, 0], []]]
    """
    F, H = dmp_to_dict(f, u), {}

    for exp, coeff in F.items():
        new_exp = [0]*len(exp)

        for e, p in zip(exp, P):
            new_exp[p] = e

        H[tuple(new_exp)] = coeff

    return dmp_from_dict(H, u, K)


def dmp_nest(f, l, K):
    """
    Return a multivariate value nested ``l``-levels.

    Examples
    ========

    >>> dmp_nest([[ZZ(1)]], 2, ZZ)
    [[[[1]]]]
    """
    if not isinstance(f, list):
        return dmp_ground(f, l)

    for i in range(l):
        f = [f]

    return f


def dmp_raise(f, l, u, K):
    """
    Return a multivariate polynomial raised ``l``-levels.

    Examples
    ========

    >>> dmp_raise([[], [ZZ(1), ZZ(2)]], 2, 1, ZZ)
    [[[[]]], [[[1]], [[2]]]]
    """
    if not l:
        return f

    if not u:
        if not f:
            return dmp_zero(l)

        k = l - 1

        return [ dmp_ground(c, k) for c in f ]

    v = u - 1

    return [ dmp_raise(c, l, v, K) for c in f ]


def dmp_deflate(f, u, K):
    """
    Map ``x_i**m_i`` to ``y_i`` in a polynomial in ``K[X]``.

    Examples
    ========

    >>> f = dmp_normal([[1, 0, 0, 2], [], [3, 0, 0, 4]], 1, ZZ)

    >>> dmp_deflate(f, 1, ZZ)
    ((2, 3), [[1, 2], [3, 4]])
    """
    if dmp_zero_p(f, u):
        return (1,)*(u + 1), f

    F = dmp_to_dict(f, u)
    B = [0]*(u + 1)

    for M in F:
        for i, m in enumerate(M):
            B[i] = math.gcd(B[i], m)

    for i, b in enumerate(B):
        if not b:
            B[i] = 1

    B = tuple(B)

    if all(b == 1 for b in B):
        return B, f

    H = {}

    for A, coeff in F.items():
        N = [a // b for a, b in zip(A, B)]
        H[tuple(N)] = coeff

    return B, dmp_from_dict(H, u, K)


def dmp_multi_deflate(polys, u, K):
    """
    Map ``x_i**m_i`` to ``y_i`` in a set of polynomials in ``K[X]``.

    Examples
    ========

    >>> f = dmp_normal([[1, 0, 0, 2], [], [3, 0, 0, 4]], 1, ZZ)
    >>> g = dmp_normal([[1, 0, 2], [], [3, 0, 4]], 1, ZZ)

    >>> dmp_multi_deflate((f, g), 1, ZZ)
    ((2, 1), ([[1, 0, 0, 2], [3, 0, 0, 4]], [[1, 0, 2], [3, 0, 4]]))
    """
    F, B = [], [0]*(u + 1)

    for p in polys:
        f = dmp_to_dict(p, u)

        if not dmp_zero_p(p, u):
            for M in f:
                for i, m in enumerate(M):
                    B[i] = math.gcd(B[i], m)

        F.append(f)

    for i, b in enumerate(B):
        if not b:
            B[i] = 1

    B = tuple(B)

    if all(b == 1 for b in B):
        return B, polys

    H = []

    for f in F:
        h = {}

        for A, coeff in f.items():
            N = [ a // b for a, b in zip(A, B) ]
            h[tuple(N)] = coeff

        H.append(dmp_from_dict(h, u, K))

    return B, tuple(H)


def dup_inflate(f, m, K):
    """
    Map ``y`` to ``x**m`` in a polynomial in ``K[x]``.

    Examples
    ========

    >>> f = dmp_normal([1, 1, 1], 0, ZZ)

    >>> dup_inflate(f, 3, ZZ)
    [1, 0, 0, 1, 0, 0, 1]
    """
    if m <= 0:
        raise IndexError("'m' must be positive, got %s" % m)
    if m == 1 or not f:
        return f

    result = [f[0]]

    for coeff in f[1:]:
        result.extend([K.zero]*(m - 1))
        result.append(coeff)

    return result


def dmp_inflate(f, M, u, K):
    """
    Map ``y_i`` to ``x_i**k_i`` in a polynomial in ``K[X]``.

    Examples
    ========

    >>> f = dmp_normal([[1, 2], [3, 4]], 1, ZZ)

    >>> dmp_inflate(f, (2, 3), 1, ZZ)
    [[1, 0, 0, 2], [], [3, 0, 0, 4]]
    """
    if not u:
        return dup_inflate(f, M[0], K)

    def inflate(g, M, v, i, K):
        if not v:
            return dup_inflate(g, M[i], K)
        if M[i] <= 0:
            raise IndexError("all M[i] must be positive, got %s" % M[i])

        w, j = v - 1, i + 1

        g = [inflate(c, M, w, j, K) for c in g]

        result = [g[0]]

        for coeff in g[1:]:
            for _ in range(1, M[i]):
                result.append(dmp_zero(w))

            result.append(coeff)

        return result

    if all(m == 1 for m in M):
        return f
    else:
        return inflate(f, M, u, 0, K)


def dmp_exclude(f, u, K):
    """
    Exclude useless levels from ``f``.

    Return the levels excluded, the new excluded ``f``, and the new ``u``.

    Examples
    ========

    >>> f = dmp_normal([[[1]], [[1], [2]]], 2, ZZ)

    >>> dmp_exclude(f, 2, ZZ)
    ([2], [[1], [1, 2]], 1)
    """
    if not u or dmp_ground_p(f, None, u):
        return [], f, u

    J, F = [], dmp_to_dict(f, u)

    for j in range(u + 1):
        for monom in F:
            if monom[j]:
                break
        else:
            J.append(j)

    if not J:
        return [], f, u

    f = {}

    for monom, coeff in F.items():
        monom = list(monom)

        for j in reversed(J):
            del monom[j]

        f[tuple(monom)] = coeff

    u -= len(J)

    return J, dmp_from_dict(f, u, K), u


def dmp_include(f, J, u, K):
    """
    Include useless levels in ``f``.

    Examples
    ========

    >>> f = dmp_normal([[1], [1, 2]], 1, ZZ)

    >>> dmp_include(f, [2], 1, ZZ)
    [[[1]], [[1], [2]]]
    """
    if not J:
        return f

    F, f = dmp_to_dict(f, u), {}

    for monom, coeff in F.items():
        monom = list(monom)

        for j in J:
            monom.insert(j, 0)

        f[tuple(monom)] = coeff

    u += len(J)

    return dmp_from_dict(f, u, K)


def dmp_inject(f, u, K, front=False):
    """
    Convert ``f`` from ``K[X][Y]`` to ``K[X,Y]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> dmp_inject([R(1), x + 2], 0, R)
    ([[[1]], [[1], [2]]], 2)
    >>> dmp_inject([R(1), x + 2], 0, R, front=True)
    ([[[1]], [[1, 2]]], 2)
    """
    f, h = dmp_to_dict(f, u), {}

    v = K.ngens - 1

    for f_monom, g in f.items():
        g = g.to_dict()

        for g_monom, c in g.items():
            if front:
                h[g_monom + f_monom] = c
            else:
                h[f_monom + g_monom] = c

    w = u + v + 1

    return dmp_from_dict(h, w, K.domain), w


def dmp_eject(f, u, K, front=False):
    """
    Convert ``f`` from ``K[X,Y]`` to ``K[X][Y]``.

    Examples
    ========

    >>> dmp_eject([[[1]], [[1], [2]]], 2, ZZ.poly_ring('x', 'y'))
    [1, x + 2]
    """
    f, h = dmp_to_dict(f, u), {}

    n = K.ngens
    v = u - K.ngens + 1

    for monom, c in f.items():
        if front:
            g_monom, f_monom = monom[:n], monom[n:]
        else:
            g_monom, f_monom = monom[-n:], monom[:-n]

        if f_monom in h:
            h[f_monom][g_monom] = c
        else:
            h[f_monom] = {g_monom: c}

    for monom, c in h.items():
        h[monom] = K(c)

    return dmp_from_dict(h, v - 1, K)


def dmp_terms_gcd(f, u, K):
    """
    Remove GCD of terms from ``f`` in ``K[X]``.

    Examples
    ========

    >>> f = dmp_normal([[1, 0], [1, 0, 0], [], []], 1, ZZ)

    >>> dmp_terms_gcd(f, 1, ZZ)
    ((2, 1), [[1], [1, 0]])
    """
    if dmp_ground_TC(f, u, K) or dmp_zero_p(f, u):
        return (0,)*(u + 1), f

    F = dmp_to_dict(f, u)
    G = monomial_min(*list(F))

    if all(g == 0 for g in G):
        return G, f

    f = {}

    for monom, coeff in F.items():
        f[monomial_div(monom, G)] = coeff

    return G, dmp_from_dict(f, u, K)


def dmp_list_terms(f, u, K, order=None):
    """
    List all non-zero terms from ``f`` in the given order ``order``.

    Examples
    ========

    >>> f = dmp_normal([[1, 1], [2, 3]], 1, ZZ)

    >>> dmp_list_terms(f, 1, ZZ)
    [((1, 1), 1), ((1, 0), 1), ((0, 1), 2), ((0, 0), 3)]
    >>> dmp_list_terms(f, 1, ZZ, order='grevlex')
    [((1, 1), 1), ((1, 0), 1), ((0, 1), 2), ((0, 0), 3)]
    """
    def sort(terms, O):
        return sorted(terms, key=lambda term: O(term[0]), reverse=True)

    def list_terms(g, v, monom):
        d, terms = dmp_degree(g, v), []

        if not v:
            for i, c in enumerate(g):
                if not c:
                    continue

                terms.append((monom + (d - i,), c))
        else:
            w = v - 1

            for i, c in enumerate(g):
                terms.extend(list_terms(c, w, monom + (d - i,)))

        return terms

    terms = list_terms(f, u, ())

    if not terms:
        return [((0,)*(u + 1), K.zero)]

    if order is None:
        return terms
    else:
        return sort(terms, monomial_key(order))


def dmp_apply_pairs(f, g, h, args, u, K):
    """
    Apply ``h`` to pairs of coefficients of ``f`` and ``g``.

    Examples
    ========

    >>> h = lambda x, y, z: 2*x + y - z

    >>> dmp_apply_pairs([[1], [2, 3]], [[3], [2, 1]], h, [1], 1, ZZ)
    [[4], [5, 6]]
    """
    if u < 0:
        return h(f, g, *args)

    n, m, v = len(f), len(g), u - 1

    if n > m:
        g = dmp_zeros(n - m, v, K) + g
    elif n < m:
        f = dmp_zeros(m - n, v, K) + f

    result = []

    for a, b in zip(f, g):
        result.append(dmp_apply_pairs(a, b, h, args, v, K))

    return dmp_strip(result, u)


def dmp_slice(f, m, n, u, K):
    """Take a continuous subsequence of terms of ``f`` in ``K[X]``. """
    return dmp_slice_in(f, m, n, 0, u, K)


def dmp_slice_in(f, m, n, j, u, K):
    """Take a continuous subsequence of terms of ``f`` in ``x_j`` in ``K[X]``. """
    if j < 0 or j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))

    if not u:
        k = len(f)
        M = k - m if k >= m else 0
        N = k - n if k >= n else 0
        f = f[N:M]
        return f + [K.zero]*m if f else []

    f, g = dmp_to_dict(f, u), {}

    for monom, coeff in f.items():
        k = monom[j]

        if k < m or k >= n:
            monom = monom[:j] + (0,) + monom[j + 1:]

        if monom in g:
            g[monom] += coeff
        else:
            g[monom] = coeff

    return dmp_from_dict(g, u, K)


def dup_random(n, a, b, K, percent=None):
    """
    Return a polynomial of degree ``n`` with coefficients in ``[a, b]``.

    If ``percent`` is a natural number less than 100 then only approximately
    the given percentage of elements will be non-zero.

    Examples
    ========

    >>> dup_random(3, -10, 10, ZZ) #doctest: +SKIP
    [-2, -8, 9, -4]
    """
    if percent is None:
        percent = 100//(b - a)
    percent = min(max(0, percent), 100)
    nz = ((n + 1)*percent)//100

    f = []
    while len(f) < n + 1:
        v = K.convert(random.randint(a, b))
        if v:
            f.append(v)

    if nz:
        f[-nz:] = [K.zero]*nz
        lt = f.pop(0)
        random.shuffle(f)
        f.insert(0, lt)

    return f
