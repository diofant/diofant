"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densearith import (dmp_add, dmp_add_term, dmp_div, dmp_exquo_ground,
                         dmp_mul, dmp_mul_ground, dmp_neg, dmp_quo_ground,
                         dmp_sub, dup_add, dup_mul)
from .densebasic import (dmp_convert, dmp_degree_in, dmp_from_dict, dmp_ground,
                         dmp_ground_LC, dmp_LC, dmp_strip, dmp_TC, dmp_to_dict,
                         dmp_zero, dmp_zero_p)
from .polyerrors import DomainError


def dmp_diff_in(f, m, j, u, K):
    """
    ``m``-th order derivative in ``x_j`` of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    >>> R.dmp_diff_in(f, 1, 0)
    y**2 + 2*y + 3
    >>> R.dmp_diff_in(f, 1, 1)
    2*x*y + 2*x + 4*y + 3

    """
    if j < 0 or j > u:
        raise IndexError("0 <= j <= %s expected, got %s" % (u, j))

    if not j:
        if m <= 0:
            return f

        n = dmp_degree_in(f, 0, u)

        if n < m:
            return dmp_zero(u)

        deriv, v = [], u - 1

        if m == 1:
            for coeff in f[:-m]:
                d = dmp_mul_ground(coeff, K(n), v, K) if u else K(n)*coeff
                deriv.append(d)
                n -= 1
        else:
            for coeff in f[:-m]:
                k = n

                for i in range(n - 1, n - m, -1):
                    k *= i

                d = dmp_mul_ground(coeff, K(k), v, K) if u else K(k)*coeff
                deriv.append(d)
                n -= 1

        return dmp_strip(deriv, u)

    def diff_in(f, m, u, i, j, K):
        if i == j:
            return dmp_diff_in(f, m, 0, u, K)

        v, i = u - 1, i + 1

        return dmp_strip([diff_in(c, m, v, i, j, K) for c in f], u)

    return diff_in(f, m, u, 0, j, K)


def dmp_eval_in(f, a, j, u, K):
    """
    Evaluate a polynomial at ``x_j = a`` in ``K[X]`` using the Horner scheme.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = 2*x*y + 3*x + y + 2

    >>> R.dmp_eval_in(f, 2, 0)
    5*y + 8
    >>> R.dmp_eval_in(f, 2, 1)
    7*x + 4

    """
    if j < 0 or j > u:
        raise IndexError("0 <= j <= %s expected, got %s" % (u, j))

    if not j:
        if not a:
            return dmp_TC(f, K)

        result, v = dmp_LC(f, K), u - 1

        if u:
            for coeff in f[1:]:
                result = dmp_mul_ground(result, a, v, K)
                result = dmp_add(result, coeff, v, K)
        else:
            for coeff in f[1:]:
                result *= a
                result += coeff

        return result

    def eval_in(g, a, v, i, j, K):
        if i == j:
            return dmp_eval_in(g, a, 0, v, K)

        v, i = v - 1, i + 1

        return dmp_strip([eval_in(c, a, v, i, j, K) for c in g], v)

    return eval_in(f, a, u, 0, j, K)


def dmp_eval_tail(f, A, u, K):
    """
    Evaluate a polynomial at ``x_j = a_j, ...`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = 2*x*y + 3*x + y + 2

    >>> R.dmp_eval_tail(f, [2])
    7*x + 4
    >>> R.dmp_eval_tail(f, [2, 2])
    18

    """
    if not A:
        return f

    if dmp_zero_p(f, u):
        return dmp_zero(u - len(A))

    def eval_tail(g, i, A, u, K):
        if i == u:
            return dmp_eval_in(g, A[-1], 0, 0, K)
        else:
            h = [eval_tail(c, i + 1, A, u, K) for c in g]

            if i < u - len(A) + 1:
                return h
            else:
                return dmp_eval_in(h, A[-u + i - 1], 0, 0, K)

    e = eval_tail(f, 0, A, u, K)

    if u == len(A) - 1:
        return e
    else:
        return dmp_strip(e, u - len(A))


def dmp_diff_eval_in(f, m, a, j, u, K):
    """
    Differentiate and evaluate a polynomial in ``x_j`` at ``a`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    >>> R.dmp_diff_eval_in(f, 1, 2, 0)
    y**2 + 2*y + 3
    >>> R.dmp_diff_eval_in(f, 1, 2, 1)
    6*x + 11

    """
    if j > u:
        raise IndexError("-%s <= j < %s expected, got %s" % (u, u, j))
    if not j:
        return dmp_eval_in(dmp_diff_in(f, m, 0, u, K), a, 0, u, K)

    def diff_eval(g, m, a, v, i, j, K):
        if i == j:
            return dmp_eval_in(dmp_diff_in(g, m, 0, v, K), a, 0, v, K)

        v, i = v - 1, i + 1

        return dmp_strip([diff_eval(c, m, a, v, i, j, K) for c in g], v)

    return diff_eval(f, m, a, u, 0, j, K)


def dup_trunc(f, p, K):
    """
    Reduce a ``K[x]`` polynomial modulo a constant ``p`` in ``K``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dmp_ground_trunc(2*x**3 + 3*x**2 + 5*x + 7, ZZ(3))
    -x**3 - x + 1

    """
    from ..ntheory.modular import symmetric_residue

    if K.is_IntegerRing:
        g = []

        for c in f:
            c = c % p
            c = symmetric_residue(c, p)
            g.append(c)
    else:
        g = [c % p for c in f]

    return dmp_strip(g, 0)


def dmp_ground_trunc(f, p, u, K):
    """
    Reduce a ``K[X]`` polynomial modulo a constant ``p`` in ``K``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3

    >>> R.dmp_ground_trunc(f, ZZ(3))
    -x**2 - x*y - y

    """
    if not u:
        return dup_trunc(f, p, K)

    v = u - 1

    return dmp_strip([dmp_ground_trunc(c, p, v, K) for c in f], u)


def dmp_ground_monic(f, u, K):
    """
    Divide all coefficients by ``LC(f)`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)
    >>> f = 3*x**2*y + 6*x**2 + 3*x*y + 9*y + 3

    >>> R.dmp_ground_monic(f)
    x**2*y + 2*x**2 + x*y + 3*y + 1

    >>> R, x, y = ring("x y", QQ)
    >>> f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3

    >>> R.dmp_ground_monic(f)
    x**2*y + 8/3*x**2 + 5/3*x*y + 2*x + 2/3*y + 1

    """
    if dmp_zero_p(f, u):
        return f

    lc = dmp_ground_LC(f, u, K)

    if lc == K.one:
        return f
    else:
        return dmp_exquo_ground(f, lc, u, K)


def dmp_ground_content(f, u, K):
    """
    Compute the GCD of coefficients of ``f`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_content(f)
    2

    >>> R, x, y = ring("x y", QQ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_content(f)
    2

    """
    if u < 0:
        return f

    if dmp_zero_p(f, u):
        return K.zero

    cont, v = K.zero, u - 1

    if K.is_RationalField:
        for c in f:
            cont = K.gcd(cont, dmp_ground_content(c, v, K))
    else:
        for c in f:
            cont = K.gcd(cont, dmp_ground_content(c, v, K))

            if cont == K.one:
                break

    if K.is_negative(dmp_ground_LC(f, u, K)):
        cont = -cont

    return cont


def dmp_ground_primitive(f, u, K):
    """
    Compute content and the primitive form of ``f`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_primitive(f)
    (2, x*y + 3*x + 2*y + 6)

    >>> R, x, y = ring("x y", QQ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> R.dmp_ground_primitive(f)
    (2, x*y + 3*x + 2*y + 6)

    """
    if dmp_zero_p(f, u):
        return K.zero, list(f)

    cont = dmp_ground_content(f, u, K)

    if cont != K.one:
        f = dmp_quo_ground(f, cont, u, K)

    return cont, list(f)


def dup_real_imag(f, K):
    """
    Return bivariate polynomials ``f1`` and ``f2``, such that ``f = f1 + f2*I``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dup_real_imag(x**3 + x**2 + x + 1)
    (x**3 + x**2 - 3*x*y**2 + x - y**2 + 1, 3*x**2*y + 2*x*y - y**3 + y)

    >>> R, x, y = ring("x y", QQ.algebraic_field(I))

    >>> R.dup_real_imag(x**2 + I*x - 1)
    (x**2 - y**2 - y - 1, 2*x*y + x)

    """
    if K.is_ComplexAlgebraicField:
        K0 = K.domain
        r1, i1 = dup_real_imag([_.real for _ in f], K0)
        r2, i2 = dup_real_imag([_.imag for _ in f], K0)
        return dmp_add(r1, dmp_neg(i2, 1, K0), 1, K0), dmp_add(r2, i1, 1, K0)
    elif not K.is_IntegerRing and not K.is_RationalField and not K.is_RealAlgebraicField:
        raise DomainError("computing real and imaginary parts is not supported over %s" % K)

    f1 = dmp_zero(1)
    f2 = dmp_zero(1)

    if not f:
        return f1, f2

    g = [[[K.one, K.zero]], [[K.one], []]]
    h = dmp_ground(f[0], 2)

    for c in f[1:]:
        h = dmp_mul(h, g, 2, K)
        h = dmp_add_term(h, dmp_ground(c, 1), 0, 2, K)

    H = dmp_to_dict(h, 0)

    for (k,), h in H.items():
        m = k % 4

        if not m:
            f1 = dmp_add(f1, h, 1, K)
        elif m == 1:
            f2 = dmp_add(f2, h, 1, K)
        elif m == 2:
            f1 = dmp_sub(f1, h, 1, K)
        else:
            f2 = dmp_sub(f2, h, 1, K)

    return f1, f2


def dup_mirror(f, K):
    """
    Evaluate efficiently the composition ``f(-x)`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_mirror(x**3 + 2*x**2 - 4*x + 2)
    -x**3 + 2*x**2 + 4*x + 2

    """
    f = list(f)

    for i in range(len(f) - 2, -1, -2):
        f[i] = -f[i]

    return f


def dup_scale(f, a, K):
    """
    Evaluate efficiently composition ``f(a*x)`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_scale(x**2 - 2*x + 1, ZZ(2))
    4*x**2 - 4*x + 1

    """
    f, n, b = list(f), len(f) - 1, a

    for i in range(n - 1, -1, -1):
        f[i], b = b*f[i], b*a

    return f


def dup_shift(f, a, K):
    """
    Evaluate efficiently Taylor shift ``f(x + a)`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_shift(x**2 - 2*x + 1, ZZ(2))
    x**2 + 2*x + 1

    """
    f, n = list(f), len(f) - 1

    for i in range(n, 0, -1):
        for j in range(i):
            f[j + 1] += a*f[j]

    return f


def dup_transform(f, p, q, K):
    """
    Evaluate functional transformation ``q**n * f(p/q)`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_transform(x**2 - 2*x + 1, x**2 + 1, x - 1)
    x**4 - 2*x**3 + 5*x**2 - 4*x + 4

    """
    if not f:
        return []

    n = len(f) - 1
    h, Q = [f[0]], [[K.one]]

    for i in range(n):
        Q.append(dup_mul(Q[-1], q, K))

    for c, q in zip(f[1:], Q[1:]):
        h = dup_mul(h, p, K)
        q = dmp_mul_ground(q, c, 0, K)
        h = dup_add(h, q, K)

    return h


def dmp_compose(f, g, u, K):
    """
    Evaluate functional composition ``f(g)`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> R.dmp_compose(x*y + 2*x + y, y)
    y**2 + 3*y

    """
    if dmp_zero_p(f, u):
        return f

    h = [f[0]]

    for c in f[1:]:
        h = dmp_mul(h, g, u, K)
        h = dmp_add_term(h, c, 0, u, K)

    return h


def _dup_right_decompose(f, s, K):
    n = len(f) - 1
    lc = dmp_LC(f, K)

    f = dmp_to_dict(f, 0)
    g = {(s,): K.one}

    r = n // s

    for i in range(1, s):
        coeff = K.zero

        for j in range(i):
            if not (n + j - i,) in f:
                continue

            assert (s - j,) in g

            fc, gc = f[(n + j - i,)], g[(s - j,)]
            coeff += (i - r*j)*fc*gc

        g[(s - i,)] = K.quo(coeff, i*r*lc)

    return dmp_from_dict(g, 0, K)


def _dup_left_decompose(f, h, K):
    g, i = {}, 0

    while f:
        q, r = dmp_div(f, h, 0, K)

        if dmp_degree_in(r, 0, 0) > 0:
            return
        else:
            g[(i,)] = dmp_LC(r, K)
            f, i = q, i + 1

    return dmp_from_dict(g, 0, K)


def _dup_decompose(f, K):
    df = len(f) - 1

    for s in range(2, df):
        if df % s != 0:
            continue

        h = _dup_right_decompose(f, s, K)
        g = _dup_left_decompose(f, h, K)

        if g is not None:
            return g, h


def dup_decompose(f, K):
    """
    Compute functional decomposition of ``f`` in ``K[x]``.

    Given a univariate polynomial ``f`` with coefficients in a field of
    characteristic zero, returns list ``[f_1, f_2, ..., f_n]``, where::

              f = f_1 o f_2 o ... f_n = f_1(f_2(... f_n))

    and ``f_2, ..., f_n`` are monic and homogeneous polynomials of at
    least second degree.

    Unlike factorization, complete functional decompositions of
    polynomials are not unique, consider examples:

    1. ``f o g = f(x + b) o (g - b)``
    2. ``x**n o x**m = x**m o x**n``
    3. ``T_n o T_m = T_m o T_n``

    where ``T_n`` and ``T_m`` are Chebyshev polynomials.

    Examples
    ========

    >>> R, x = ring("x", ZZ)

    >>> R.dup_decompose(x**4 - 2*x**3 + x**2)
    [x**2, x**2 - x]

    References
    ==========

    * :cite:`Kozen1989decomposition`

    """
    F = []

    while True:
        result = _dup_decompose(f, K)

        if result is not None:
            f, h = result
            F = [h] + F
        else:
            break

    return [f] + F


def dmp_clear_denoms(f, u, K0, K1=None, convert=False):
    """
    Clear denominators, i.e. transform ``K_0`` to ``K_1``.

    Examples
    ========

    >>> R, x, y = ring("x y", QQ)

    >>> f = x/2 + y/3 + 1

    >>> R.dmp_clear_denoms(f, convert=False)
    (6, 3*x + 2*y + 6)
    >>> R.dmp_clear_denoms(f, convert=True)
    (6, 3*x + 2*y + 6)

    """
    if K1 is None:
        if K0.has_assoc_Ring:
            K1 = K0.ring
        else:
            K1 = K0

    def clear_denoms(g, v, K0, K1):
        common = K1.one

        if not v:
            for c in g:
                common = K1.lcm(common, c.denominator)
        else:
            w = v - 1

            for c in g:
                common = K1.lcm(common, clear_denoms(c, w, K0, K1))

        return common

    common = clear_denoms(f, u, K0, K1)
    f = dmp_mul_ground(f, common, u, K0)

    if not convert:
        return common, f
    else:
        return common, dmp_convert(f, u, K0, K1)
