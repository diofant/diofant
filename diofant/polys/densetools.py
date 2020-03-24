"""Advanced tools for dense recursive polynomials in ``K[x]`` or ``K[X]``."""

from .densearith import (dmp_add, dmp_add_term, dmp_div, dmp_exquo_ground,
                         dmp_mul, dmp_mul_ground, dmp_neg, dmp_sub, dup_add,
                         dup_mul)
from .densebasic import (dmp_degree_in, dmp_from_dict, dmp_ground,
                         dmp_ground_LC, dmp_LC, dmp_strip, dmp_to_dict,
                         dmp_zero, dmp_zero_p)
from .polyerrors import DomainError


def dmp_diff_in(f, m, j, u, K):
    """
    ``m``-th order derivative in ``x_j`` of a polynomial in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    >>> f.diff()
    y**2 + 2*y + 3
    >>> f.diff(y)
    2*x*y + 2*x + 4*y + 3

    """
    ring = K.poly_ring(*['_%d' % i for i in range(u + 1)])
    f = ring.from_dense(f)
    return ring.to_dense(f.diff(x=j, m=m))


def dmp_eval_in(f, a, j, u, K):
    """
    Evaluate a polynomial at ``x_j = a`` in ``K[X]`` using the Horner scheme.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> f = 2*x*y + 3*x + y + 2

    >>> f.eval(x, 2)
    5*y + 8
    >>> f.eval(y, 2)
    7*x + 4

    """
    ring = K.poly_ring(*['_%d' % i for i in range(u + 1)])
    f = ring.from_dense(f)
    r = f.eval(j, a)
    if u:
        return ring.drop(j).to_dense(r)
    else:
        return r


def dmp_eval_tail(f, A, u, K):
    """
    Evaluate a polynomial at ``x_j = a_j, ...`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

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

    >>> R, x, y = ring('x y', ZZ)

    >>> f = x*y**2 + 2*x*y + 3*x + 2*y**2 + 3*y + 1

    >>> f.diff().eval(x, 2)
    y**2 + 2*y + 3
    >>> f.diff(x=y).eval(y, 2)
    6*x + 11

    """
    ring = K.poly_ring(*['_%d' % i for i in range(u + 1)])
    f = ring.from_dense(f)
    r = f.diff(x=j, m=m).eval(j, a)
    return ring.drop(j).to_dense(r) if u else r


def dmp_ground_trunc(f, p, u, K):
    """
    Reduce a ``K[X]`` polynomial modulo a constant ``p`` in ``K``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

    >>> (2*x**3 + 3*x**2 + 5*x + 7).trunc_ground(ZZ(3))
    -x**3 - x + 1

    >>> R, x, y = ring('x y', ZZ)

    >>> f = 3*x**2*y + 8*x**2 + 5*x*y + 6*x + 2*y + 3

    >>> f.trunc_ground(ZZ(3))
    -x**2 - x*y - y

    """
    ring = K.poly_ring(*['_%d' % i for i in range(u + 1)])
    f = ring.from_dense(f)
    return ring.to_dense(f.trunc_ground(p))


def dmp_ground_monic(f, u, K):
    """
    Divide all coefficients by ``LC(f)`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)
    >>> f = 3*x**2*y + 6*x**2 + 3*x*y + 9*y + 3

    >>> R.dmp_ground_monic(f)
    x**2*y + 2*x**2 + x*y + 3*y + 1

    >>> R, x, y = ring('x y', QQ)
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

    >>> R, x, y = ring('x y', ZZ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> f.content()
    2

    >>> R, x, y = ring('x y', QQ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> f.content()
    2

    """
    ring = K.poly_ring(*['_%d' % i for i in range(u + 1)])
    f = ring.from_dense(f)
    return f.content()


def dmp_ground_primitive(f, u, K):
    """
    Compute content and the primitive form of ``f`` in ``K[X]``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> f.primitive()
    (2, x*y + 3*x + 2*y + 6)

    >>> R, x, y = ring('x y', QQ)
    >>> f = 2*x*y + 6*x + 4*y + 12

    >>> f.primitive()
    (2, x*y + 3*x + 2*y + 6)

    """
    ring = K.poly_ring(*['_%d' % i for i in range(u + 1)])
    f = ring.from_dense(f)
    cont, p = f.primitive()
    return cont, ring.to_dense(p)


def dup_real_imag(f, K):
    """
    Return bivariate polynomials ``f1`` and ``f2``, such that ``f = f1 + f2*I``.

    Examples
    ========

    >>> R, x, y = ring('x y', ZZ)

    >>> R.dup_real_imag(x**3 + x**2 + x + 1)
    (x**3 + x**2 - 3*x*y**2 + x - y**2 + 1, 3*x**2*y + 2*x*y - y**3 + y)

    >>> R, x, y = ring('x y', QQ.algebraic_field(I))

    >>> R.dup_real_imag(x**2 + I*x - 1)
    (x**2 - y**2 - y - 1, 2*x*y + x)

    """
    if K.is_ComplexAlgebraicField:
        K0 = K.domain
        r1, i1 = dup_real_imag([_.real for _ in f], K0)
        r2, i2 = dup_real_imag([_.imag for _ in f], K0)
        return dmp_add(r1, dmp_neg(i2, 1, K0), 1, K0), dmp_add(r2, i1, 1, K0)
    elif not K.is_IntegerRing and not K.is_RationalField and not K.is_RealAlgebraicField:
        raise DomainError('computing real and imaginary parts is not supported over %s' % K)

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


def dup_transform(f, p, q, K):
    """
    Evaluate functional transformation ``q**n * f(p/q)`` in ``K[x]``.

    Examples
    ========

    >>> R, x = ring('x', ZZ)

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

    >>> R, x, y = ring('x y', ZZ)

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

    >>> R, x = ring('x', ZZ)

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


def dmp_clear_denoms(f, u, K0, convert=False):
    """
    Clear denominators, i.e. transform ``K_0`` to ``K_1``.

    Examples
    ========

    >>> R, x, y = ring('x y', QQ)

    >>> f = x/2 + y/3 + 1

    >>> f.clear_denoms()
    (6, 3*x + 2*y + 6)

    """
    ring = K0.poly_ring(*['_%d' % i for i in range(u + 1)])
    f = ring.from_dense(f)
    common, f = f.clear_denoms()
    return common, ring.to_dense(f)
