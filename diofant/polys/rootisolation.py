"""Real and complex root isolation and refinement algorithms. """

from ..core import I
from .densearith import dmp_add, dmp_mul_ground, dmp_neg, dmp_rem, dup_rshift
from .densebasic import (dmp_convert, dmp_degree_in, dmp_LC, dmp_permute,
                         dmp_strip, dmp_TC, dmp_terms_gcd, dmp_to_tuple,
                         dup_reverse)
from .densetools import (dmp_clear_denoms, dmp_compose, dmp_diff_in,
                         dmp_eval_in, dup_mirror, dup_real_imag, dup_scale,
                         dup_shift, dup_sign_variations, dup_transform)
from .euclidtools import dmp_gcd, dmp_resultant
from .factortools import dmp_factor_list
from .polyerrors import DomainError, RefinementFailed
from .sqfreetools import dmp_sqf_list, dmp_sqf_part


def dup_sturm(f, K):
    """
    Computes the Sturm sequence of ``f`` in ``F[x]``.

    Given a univariate, square-free polynomial ``f(x)`` returns the
    associated Sturm sequence (see e.g. [Davenport88]_)
    ``f_0(x), ..., f_n(x)`` defined by::

       f_0(x), f_1(x) = f(x), f'(x)
       f_n = -rem(f_{n-2}(x), f_{n-1}(x))

    Examples
    ========

    >>> R, x = ring("x", QQ)

    >>> R.dup_sturm(x**3 - 2*x**2 + x - 3)
    [x**3 - 2*x**2 + x - 3, 3*x**2 - 4*x + 1, 2/9*x + 25/9, -2079/4]
    """
    if not K.is_Field:
        raise DomainError("can't compute Sturm sequence over %s" % K)

    f = dmp_sqf_part(f, 0, K)

    sturm = [f, dmp_diff_in(f, 1, 0, 0, K)]

    while sturm[-1]:
        s = dmp_rem(sturm[-2], sturm[-1], 0, K)
        sturm.append(dmp_neg(s, 0, K))

    return sturm[:-1]


def dup_root_upper_bound(f, K):
    """Compute the LMQ upper bound for the positive roots of `f`.

    LMQ (Local Max Quadratic) bound was developed by
    Akritas-Strzebonski-Vigklas [Alkiviadis09]_.
    """
    n, P = len(f), []
    if K.is_AlgebraicField:
        return
    t = n * [K.one]
    if dmp_LC(f, K) < 0:
        f = dmp_neg(f, 0, K)
    f = list(reversed(f))

    for i in range(n):
        if int(f[i]) >= 0:
            continue

        a, QL = K.log(-f[i], 2), []

        for j in range(i + 1, n):

            if int(f[j]) <= 0:
                continue

            q = t[j] + a - K.log(f[j], 2)
            QL.append([q // (j - i), j])

        assert QL

        q = min(QL)

        t[q[1]] = t[q[1]] + 1

        P.append(q[0])

    if P:
        return K(2)**int(max(P) + 1)


def _mobius_from_interval(I, field):
    """Convert an open interval to a Mobius transform. """
    s, t = I

    a, c = map(field, (s.numerator, s.denominator))
    b, d = map(field, (t.numerator, t.denominator))

    return a, b, c, d


def _mobius_to_interval(M):
    """Convert a Mobius transform to an open interval. """
    a, b, c, d = M

    s, t = a/c, b/d

    return (s, t) if s <= t else (t, s)


def dup_step_refine_real_root(f, M, K, fast=False):
    """One step of positive real root refinement algorithm. """
    a, b, c, d = M

    if a == b and c == d:
        return f, (a, b, c, d)

    A = dup_root_upper_bound(dup_reverse(f), K)

    if A is not None:
        A = 1/K.convert(A)
    else:
        A = K.zero

    if fast and A > 16:
        f = dup_scale(f, A, K)
        a, c, A = A*a, A*c, K.one

    if A >= 1:
        f = dup_shift(f, A, K)
        b, d = A*a + b, A*c + d

        assert dmp_TC(f, K)

    f, g = dup_shift(f, K.one, K), f

    a1, b1, c1, d1 = a, a + b, c, c + d

    if not dmp_eval_in(f, K.zero, 0, 0, K):
        return f, (b1, b1, d1, d1)

    k = dup_sign_variations(f, K)

    if k == 1:
        a, b, c, d = a1, b1, c1, d1
    else:
        f = dup_shift(dup_reverse(g), K.one, K)

        assert dmp_TC(f, K)

        a, b, c, d = b, a + b, d, c + d

    return f, (a, b, c, d)


def dup_inner_refine_real_root(f, M, K, eps=None, steps=None, disjoint=None, fast=False, mobius=False):
    """Refine a positive root of `f` given a Mobius transform or an interval. """
    a, b, c, d = M

    while not c:
        f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c,
                                                        d), K, fast=fast)

    if eps is not None and steps is not None:
        for i in range(steps):
            if abs(a/c - b/d) >= eps:
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)
            else:
                break
    else:
        if eps is not None:
            while abs(a/c - b/d) >= eps:
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)

        if steps is not None:
            for i in range(steps):
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)

    if disjoint is not None:
        while True:
            u, v = _mobius_to_interval((a, b, c, d))

            if u < disjoint < v:
                f, (a, b, c, d) = dup_step_refine_real_root(f, (a, b, c, d), K, fast=fast)
            else:
                break

    return (f, (a, b, c, d)) if mobius else _mobius_to_interval((a, b, c, d))


def dup_outer_refine_real_root(f, s, t, K, eps=None, steps=None, disjoint=None, fast=False):
    """Refine a positive root of `f` given an interval `(s, t)`. """
    a, b, c, d = _mobius_from_interval((s, t), K)

    f = dup_transform(f, dmp_strip([a, b], 0), dmp_strip([c, d], 0), K)

    if dup_sign_variations(f, K) != 1:
        raise RefinementFailed("there should be exactly one root in (%s, %s) interval" % (s, t))

    return dup_inner_refine_real_root(f, (a, b, c, d), K, eps=eps, steps=steps, disjoint=disjoint, fast=fast)


def dup_refine_real_root(f, s, t, K, eps=None, steps=None, disjoint=None, fast=False):
    """Refine real root's approximating interval to the given precision. """
    R, K = K, K.field
    f = dmp_convert(f, 0, R, K)
    f = dmp_clear_denoms(f, 0, K)[1]

    if not (K.is_RationalField or K.is_RealAlgebraicField):
        raise DomainError("real root refinement not supported over %s" % K)

    if s == t:
        return s, t

    if s > t:
        s, t = t, s

    negative = False

    if s < 0:
        if t <= 0:
            f, s, t, negative = dup_mirror(f, K), -t, -s, True
        else:
            raise ValueError("can't refine a real root in (%s, %s)" % (s, t))

    if negative and disjoint is not None:
        if disjoint < 0:
            disjoint = -disjoint
        else:
            disjoint = None

    s, t = dup_outer_refine_real_root(f, s, t, K, eps=eps, steps=steps, disjoint=disjoint, fast=fast)

    return (-t, -s) if negative else (s, t)


def dup_inner_isolate_real_roots(f, K, eps=None, fast=False):
    """Internal function for isolation positive roots up to given precision. """
    a, b, c, d = K.one, K.zero, K.zero, K.one
    k = dup_sign_variations(f, K)

    roots, stack = [], [(a, b, c, d, f, k)]

    while stack:
        a, b, c, d, f, k = stack.pop()

        A = dup_root_upper_bound(dup_reverse(f), K)

        if A is not None:
            A = 1/K.convert(A)
        else:
            A = K.zero

        if fast and A > 16:
            f = dup_scale(f, A, K)
            a, c, A = A*a, A*c, K.one

        if A >= 1:
            f = dup_shift(f, A, K)
            b, d = A*a + b, A*c + d

            assert dmp_TC(f, K)

            k = dup_sign_variations(f, K)

            if k == 0:
                continue
            if k == 1:
                roots.append(dup_inner_refine_real_root(
                    f, (a, b, c, d), K, eps=eps, fast=fast, mobius=True))
                continue

        f1 = dup_shift(f, K.one, K)

        a1, b1, c1, d1, r = a, a + b, c, c + d, 0

        if not dmp_TC(f1, K):
            roots.append((f1, (b1, b1, d1, d1)))
            f1, r = dup_rshift(f1, 1, K), 1

        k1 = dup_sign_variations(f1, K)
        k2 = k - k1 - r

        a2, b2, c2, d2 = b, a + b, d, c + d

        if k2 > 1:
            f2 = dup_shift(dup_reverse(f), K.one, K)

            if not dmp_TC(f2, K):
                f2 = dup_rshift(f2, 1, K)

            k2 = dup_sign_variations(f2, K)
        else:
            f2 = None

        if k1 < k2:
            a1, a2, b1, b2 = a2, a1, b2, b1
            c1, c2, d1, d2 = c2, c1, d2, d1
            f1, f2, k1, k2 = f2, f1, k2, k1

        if not k1:
            continue

        if f1 is None:
            f1 = dup_shift(dup_reverse(f), K.one, K)

            if not dmp_TC(f1, K):
                f1 = dup_rshift(f1, 1, K)

        if k1 == 1:
            roots.append(dup_inner_refine_real_root(
                f1, (a1, b1, c1, d1), K, eps=eps, fast=fast, mobius=True))
        else:
            stack.append((a1, b1, c1, d1, f1, k1))

        if not k2:
            continue

        if f2 is None:
            f2 = dup_shift(dup_reverse(f), K.one, K)

            if not dmp_TC(f2, K):
                f2 = dup_rshift(f2, 1, K)

        if k2 == 1:
            roots.append(dup_inner_refine_real_root(
                f2, (a2, b2, c2, d2), K, eps=eps, fast=fast, mobius=True))
        else:
            stack.append((a2, b2, c2, d2, f2, k2))

    return roots


def _discard_if_outside_interval(f, M, inf, sup, K, negative, fast, mobius):
    """Discard an isolating interval if outside ``(inf, sup)``. """
    while True:
        u, v = _mobius_to_interval(M)

        if negative:
            u, v = -v, -u

        if (inf is None or u >= inf) and (sup is None or v <= sup):
            if not mobius:
                return u, v
            else:
                return f, M
        elif (sup is not None and u > sup) or (inf is not None and v < inf):
            return
        else:
            f, M = dup_step_refine_real_root(f, M, K, fast=fast)


def dup_inner_isolate_positive_roots(f, K, eps=None, inf=None, sup=None, fast=False, mobius=False):
    """Iteratively compute disjoint positive root isolation intervals. """
    if sup is not None and sup < 0:
        return []

    roots = dup_inner_isolate_real_roots(f, K, eps=eps, fast=fast)

    results = []

    if inf is not None or sup is not None:
        for f, M in roots:
            result = _discard_if_outside_interval(f, M, inf, sup, K, False, fast, mobius)

            if result is not None:
                results.append(result)
    elif not mobius:
        for f, M in roots:
            u, v = _mobius_to_interval(M)
            results.append((u, v))
    else:
        results = roots

    return results


def dup_inner_isolate_negative_roots(f, K, inf=None, sup=None, eps=None, fast=False, mobius=False):
    """Iteratively compute disjoint negative root isolation intervals. """
    if inf is not None and inf >= 0:
        return []

    roots = dup_inner_isolate_real_roots(dup_mirror(f, K), K, eps=eps, fast=fast)

    results = []

    if inf is not None or sup is not None:
        for f, M in roots:
            result = _discard_if_outside_interval(f, M, inf, sup, K, True, fast, mobius)

            if result is not None:
                results.append(result)
    elif not mobius:
        for f, M in roots:
            u, v = _mobius_to_interval(M)
            results.append((-v, -u))
    else:
        results = roots

    return results


def _isolate_zero(f, K, inf, sup, sqf=False):
    """Handle special case of CF algorithm when ``f`` is homogeneous. """
    (j,), f = dmp_terms_gcd(f, 0, K)

    if j > 0:
        if (inf is None or inf <= 0) and (sup is None or 0 <= sup):
            if not sqf:
                return [((K.zero, K.zero), j)], f
            else:
                return [(K.zero, K.zero)], f

    return [], f


def dup_isolate_real_roots_sqf(f, K, eps=None, inf=None, sup=None, fast=False, blackbox=False):
    """Isolate real roots of a square-free polynomial. """
    R, K = K, K.field
    f = dmp_convert(f, 0, R, K)
    f = dmp_clear_denoms(f, 0, K)[1]

    if K.is_AlgebraicField and not K.is_RealAlgebraicField:
        A, K = K, K.domain
        polys = [dmp_eval_in(_, K.zero, 1, 1, K) for _ in dup_real_imag(f, A)]
        if not polys[1]:
            f = polys[0]
        else:
            roots = dup_isolate_real_roots_list(polys, K, eps=eps, inf=inf, sup=sup, strict=True)
            roots = [_[0] for _ in roots if _[1].keys() == {0, 1}]
            return [RealInterval((a, b), f, K) for (a, b) in roots] if blackbox else roots

    if not (K.is_RationalField or K.is_RealAlgebraicField):
        raise DomainError("isolation of real roots not supported over %s" % K)

    if dmp_degree_in(f, 0, 0) <= 0:
        return []

    I_zero, f = _isolate_zero(f, K, inf, sup, sqf=True)

    I_neg = dup_inner_isolate_negative_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast)
    I_pos = dup_inner_isolate_positive_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast)

    roots = sorted(I_neg + I_zero + I_pos)
    return [RealInterval((a, b), f, K) for (a, b) in roots] if blackbox else roots


def dup_isolate_real_roots(f, K, eps=None, inf=None, sup=None, fast=False):
    """Isolate real roots.

    Notes
    =====

    Implemented algorithms use Vincent-Akritas-Strzebonski (VAS) continued
    fractions approach [Alkiviadis05]_, [Alkiviadis08]_.
    """
    R, K = K, K.field
    f = dmp_convert(f, 0, R, K)

    _, factors = dmp_sqf_list(f, 0, K)
    factors = [(dmp_clear_denoms(f, 0, K)[1], k) for f, k in factors]

    if len(factors) == 1:
        (f, k), = factors
        return [(r, k) for r in dup_isolate_real_roots_sqf(f, K, eps, inf, sup, fast)]
    else:
        if not (K.is_RationalField or K.is_RealAlgebraicField):
            raise DomainError("isolation of real roots not supported over %s" % K)

        I_zero, f = _isolate_zero(f, K, inf, sup)
        I_neg, I_pos = _real_isolate_and_disjoin(factors, K, eps, inf, sup, fast=fast)
        return sorted(I_neg + I_zero + I_pos)


def dup_isolate_real_roots_list(polys, K, eps=None, inf=None, sup=None, strict=False, basis=False, fast=False):
    """Isolate real roots of a list of polynomials. """
    R, K = K, K.field
    for i, p in enumerate(polys):
        p = dmp_convert(p, 0, R, K)
        polys[i] = dmp_clear_denoms(p, 0, K)[1]

    if not (K.is_RationalField or K.is_RealAlgebraicField):
        raise DomainError("isolation of real roots not supported over %s" % K)

    zeros, factors_dict = False, {}

    if (inf is None or inf <= 0) and (sup is None or 0 <= sup):
        zeros, zero_indices = True, {}

    for i, p in enumerate(polys):
        (j,), p = dmp_terms_gcd(p, 0, K)

        if zeros and j > 0:
            zero_indices[i] = j

        for f, k in dmp_factor_list(p, 0, K)[1]:
            f = tuple(f)

            if f not in factors_dict:
                factors_dict[f] = {i: k}
            else:
                factors_dict[f][i] = k

    factors_list = []

    for f, indices in factors_dict.items():
        factors_list.append((list(f), indices))

    I_neg, I_pos = _real_isolate_and_disjoin(factors_list, K, eps=eps,
                                             inf=inf, sup=sup, strict=strict, basis=basis, fast=fast)

    if not zeros or not zero_indices:
        I_zero = []
    else:
        if not basis:
            I_zero = [((K.zero, K.zero), zero_indices)]
        else:
            I_zero = [((K.zero, K.zero), zero_indices, [K.one, K.zero])]

    return sorted(I_neg + I_zero + I_pos)


def _disjoint_p(M, N, strict=False):
    """Check if Mobius transforms define disjoint intervals. """
    a1, b1, c1, d1 = M
    a2, b2, c2, d2 = N

    a1d1, b1c1 = a1*d1, b1*c1
    a2d2, b2c2 = a2*d2, b2*c2

    if a1d1 == b1c1 and a2d2 == b2c2:
        return True

    if a1d1 > b1c1:
        a1, c1, b1, d1 = b1, d1, a1, c1

    if a2d2 > b2c2:
        a2, c2, b2, d2 = b2, d2, a2, c2

    if not strict:
        return a2*d1 >= c2*b1 or b2*c1 <= d2*a1
    else:
        return a2*d1 > c2*b1 or b2*c1 < d2*a1


def _real_isolate_and_disjoin(factors, K, eps=None, inf=None, sup=None, strict=False, basis=False, fast=False):
    """Isolate real roots of a list of polynomials and disjoin intervals. """
    I_pos, I_neg = [], []

    for i, (f, k) in enumerate(factors):
        for F, M in dup_inner_isolate_positive_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast, mobius=True):
            I_pos.append((F, M, k, f))

        for G, N in dup_inner_isolate_negative_roots(f, K, eps=eps, inf=inf, sup=sup, fast=fast, mobius=True):
            I_neg.append((G, N, k, f))

    for i, (f, M, k, F) in enumerate(I_pos):
        for j, (g, N, m, G) in enumerate(I_pos[i + 1:]):
            while not _disjoint_p(M, N, strict=strict):
                f, M = dup_inner_refine_real_root(f, M, K, steps=1, fast=fast, mobius=True)
                g, N = dup_inner_refine_real_root(g, N, K, steps=1, fast=fast, mobius=True)

            I_pos[i + j + 1] = g, N, m, G

        I_pos[i] = f, M, k, F

    for i, (f, M, k, F) in enumerate(I_neg):
        for j, (g, N, m, G) in enumerate(I_neg[i + 1:]):
            while not _disjoint_p(M, N, strict=strict):
                f, M = dup_inner_refine_real_root(f, M, K, steps=1, fast=fast, mobius=True)
                g, N = dup_inner_refine_real_root(g, N, K, steps=1, fast=fast, mobius=True)

            I_neg[i + j + 1] = g, N, m, G

        I_neg[i] = f, M, k, F

    if strict:
        for i, (f, M, k, F) in enumerate(I_neg):
            if not M[0]:
                while not M[0]:
                    f, M = dup_inner_refine_real_root(f, M, K, steps=1, fast=fast, mobius=True)

                I_neg[i] = f, M, k, F
                break

        for j, (g, N, m, G) in enumerate(I_pos):
            if not N[0]:
                while not N[0]:
                    g, N = dup_inner_refine_real_root(g, N, K, steps=1, fast=fast, mobius=True)

                I_pos[j] = g, N, m, G
                break

    I_neg = [(_mobius_to_interval(M), k, f) for (_, M, k, f) in I_neg]
    I_pos = [(_mobius_to_interval(M), k, f) for (_, M, k, f) in I_pos]

    if not basis:
        I_neg = [((-v, -u), k) for ((u, v), k, _) in I_neg]
        I_pos = [((+u, +v), k) for ((u, v), k, _) in I_pos]
    else:
        I_neg = [((-v, -u), k, f) for ((u, v), k, f) in I_neg]
        I_pos = [((+u, +v), k, f) for ((u, v), k, f) in I_pos]

    return I_neg, I_pos


def dup_count_real_roots(f, K, inf=None, sup=None):
    """Returns the number of distinct real roots of ``f`` in ``[inf, sup]``. """
    if dmp_degree_in(f, 0, 0) <= 0:
        return 0

    if not K.is_Field:
        R, K = K, K.field
        f = dmp_convert(f, 0, R, K)

    if K.is_AlgebraicField:
        return sum(k for *_, k in dup_isolate_real_roots(f, K, inf, sup))

    sturm = dup_sturm(f, K)

    if inf is None:
        signs_inf = dup_sign_variations([dmp_LC(s, K)*(-1)**dmp_degree_in(s, 0, 0) for s in sturm], K)
    else:
        signs_inf = dup_sign_variations([dmp_eval_in(s, inf, 0, 0, K) for s in sturm], K)

    if sup is None:
        signs_sup = dup_sign_variations([dmp_LC(s, K) for s in sturm], K)
    else:
        signs_sup = dup_sign_variations([dmp_eval_in(s, sup, 0, 0, K) for s in sturm], K)

    count = abs(signs_inf - signs_sup)

    if inf is not None and not dmp_eval_in(f, inf, 0, 0, K):
        count += 1

    return count


def dup_isolate_imaginary_roots(f, K, eps=None, inf=None, sup=None, fast=False):
    """Isolate imaginary roots. """
    F = K.algebraic_field(I)
    f = dmp_compose(dmp_convert(f, 0, K, F), [F.unit, 0], 0, F)
    return dup_isolate_real_roots(f, F, eps=eps, inf=inf, sup=sup, fast=fast)


OO = 'OO'  # Origin of (re, im) coordinate system

Q1 = 'Q1'  # Quadrant #1 (++): re > 0 and im > 0
Q2 = 'Q2'  # Quadrant #2 (-+): re < 0 and im > 0
Q3 = 'Q3'  # Quadrant #3 (--): re < 0 and im < 0
Q4 = 'Q4'  # Quadrant #4 (+-): re > 0 and im < 0

A1 = 'A1'  # Axis #1 (+0): re > 0 and im = 0
A2 = 'A2'  # Axis #2 (0+): re = 0 and im > 0
A3 = 'A3'  # Axis #3 (-0): re < 0 and im = 0
A4 = 'A4'  # Axis #4 (0-): re = 0 and im < 0

_rules_simple = {
    # Q --> Q (same) => no change
    (Q1, Q1): 0,
    (Q2, Q2): 0,
    (Q3, Q3): 0,
    (Q4, Q4): 0,

    # A -- CCW --> Q => +1/4 (CCW)
    (A1, Q1): 1,
    (A2, Q2): 1,
    (A3, Q3): 1,
    (A4, Q4): 1,

    # A --  CW --> Q => -1/4 (CCW)
    (A1, Q4): 2,
    (A2, Q1): 2,
    (A3, Q2): 2,
    (A4, Q3): 2,

    # Q -- CCW --> A => +1/4 (CCW)
    (Q1, A2): 3,
    (Q2, A3): 3,
    (Q3, A4): 3,
    (Q4, A1): 3,

    # Q --  CW --> A => -1/4 (CCW)
    (Q1, A1): 4,
    (Q2, A2): 4,
    (Q3, A3): 4,
    (Q4, A4): 4,

    # Q -- CCW --> Q => +1/2 (CCW)
    (Q1, Q2): +5,
    (Q2, Q3): +5,
    (Q3, Q4): +5,
    (Q4, Q1): +5,

    # Q --  CW --> Q => -1/2 (CW)
    (Q1, Q4): -5,
    (Q2, Q1): -5,
    (Q3, Q2): -5,
    (Q4, Q3): -5,
}

_rules_ambiguous = {
    # A -- CCW --> Q => { +1/4 (CCW), -9/4 (CW) }
    (A1, OO, Q1): -1,
    (A2, OO, Q2): -1,
    (A3, OO, Q3): -1,
    (A4, OO, Q4): -1,

    # A --  CW --> Q => { -1/4 (CCW), +7/4 (CW) }
    (A1, OO, Q4): -2,
    (A2, OO, Q1): -2,
    (A3, OO, Q2): -2,
    (A4, OO, Q3): -2,

    # Q -- CCW --> A => { +1/4 (CCW), -9/4 (CW) }
    (Q1, OO, A2): -3,
    (Q2, OO, A3): -3,
    (Q3, OO, A4): -3,
    (Q4, OO, A1): -3,

    # Q --  CW --> A => { -1/4 (CCW), +7/4 (CW) }
    (Q1, OO, A1): -4,
    (Q2, OO, A2): -4,
    (Q3, OO, A3): -4,
    (Q4, OO, A4): -4,

    # A --  OO --> A => { +1 (CCW), -1 (CW) }
    (A1, A3): 7,
    (A2, A4): 7,
    (A3, A1): 7,
    (A4, A2): 7,

    (A1, OO, A3): 7,
    (A2, OO, A4): 7,
    (A3, OO, A1): 7,
    (A4, OO, A2): 7,

    # Q -- DIA --> Q => { +1 (CCW), -1 (CW) }
    (Q1, Q3): 8,
    (Q2, Q4): 8,
    (Q3, Q1): 8,
    (Q4, Q2): 8,

    (Q1, OO, Q3): 8,
    (Q2, OO, Q4): 8,
    (Q3, OO, Q1): 8,
    (Q4, OO, Q2): 8,

    # A --- R ---> A => { +1/2 (CCW), -3/2 (CW) }
    (A1, A2): 9,
    (A2, A3): 9,
    (A3, A4): 9,
    (A4, A1): 9,

    (A1, OO, A2): 9,
    (A2, OO, A3): 9,
    (A3, OO, A4): 9,
    (A4, OO, A1): 9,

    # A --- L ---> A => { +3/2 (CCW), -1/2 (CW) }
    (A1, A4): 10,
    (A2, A1): 10,
    (A3, A2): 10,
    (A4, A3): 10,

    (A1, OO, A4): 10,
    (A2, OO, A1): 10,
    (A3, OO, A2): 10,
    (A4, OO, A3): 10,

    # Q --- 1 ---> A => { +3/4 (CCW), -5/4 (CW) }
    (Q1, A3): 11,
    (Q2, A4): 11,
    (Q3, A1): 11,
    (Q4, A2): 11,

    (Q1, OO, A3): 11,
    (Q2, OO, A4): 11,
    (Q3, OO, A1): 11,
    (Q4, OO, A2): 11,

    # Q --- 2 ---> A => { +5/4 (CCW), -3/4 (CW) }
    (Q1, A4): 12,
    (Q2, A1): 12,
    (Q3, A2): 12,
    (Q4, A3): 12,

    (Q1, OO, A4): 12,
    (Q2, OO, A1): 12,
    (Q3, OO, A2): 12,
    (Q4, OO, A3): 12,

    # A --- 1 ---> Q => { +5/4 (CCW), -3/4 (CW) }
    (A1, Q3): 13,
    (A2, Q4): 13,
    (A3, Q1): 13,
    (A4, Q2): 13,

    (A1, OO, Q3): 13,
    (A2, OO, Q4): 13,
    (A3, OO, Q1): 13,
    (A4, OO, Q2): 13,

    # A --- 2 ---> Q => { +3/4 (CCW), -5/4 (CW) }
    (A1, Q2): 14,
    (A2, Q3): 14,
    (A3, Q4): 14,
    (A4, Q1): 14,

    (A1, OO, Q2): 14,
    (A2, OO, Q3): 14,
    (A3, OO, Q4): 14,
    (A4, OO, Q1): 14,

    # Q --> OO --> Q => { +1/2 (CCW), -3/2 (CW) }
    (Q1, OO, Q2): 15,
    (Q2, OO, Q3): 15,
    (Q3, OO, Q4): 15,
    (Q4, OO, Q1): 15,

    # Q --> OO --> Q => { +3/2 (CCW), -1/2 (CW) }
    (Q1, OO, Q4): 16,
    (Q2, OO, Q1): 16,
    (Q3, OO, Q2): 16,
    (Q4, OO, Q3): 16,

    # A --> OO --> A => { +2 (CCW), 0 (CW) }
    (A1, OO, A1): 17,
    (A2, OO, A2): 17,
    (A3, OO, A3): 17,
    (A4, OO, A4): 17,

    # Q --> OO --> Q => { +2 (CCW), 0 (CW) }
    (Q1, OO, Q1): 18,
    (Q2, OO, Q2): 18,
    (Q3, OO, Q3): 18,
    (Q4, OO, Q4): 18,
}

_values = {
    0: [(+0, 1)],
    1: [(+1, 4)],
    2: [(-1, 4)],
    3: [(+1, 4)],
    4: [(-1, 4)],
    -1: [(+9, 4), (+1, 4)],
    -2: [(+7, 4), (-1, 4)],
    -3: [(+9, 4), (+1, 4)],
    -4: [(+7, 4), (-1, 4)],
    +5: [(+1, 2)],
    -5: [(-1, 2)],
    7: [(+1, 1), (-1, 1)],
    8: [(+1, 1), (-1, 1)],
    9: [(+1, 2), (-3, 2)],
    10: [(+3, 2), (-1, 2)],
    11: [(+3, 4), (-5, 4)],
    12: [(+5, 4), (-3, 4)],
    13: [(+5, 4), (-3, 4)],
    14: [(+3, 4), (-5, 4)],
    15: [(+1, 2), (-3, 2)],
    16: [(+3, 2), (-1, 2)],
    17: [(+2, 1), (+0, 1)],
    18: [(+2, 1), (+0, 1)],
}


def _classify_point(re, im):
    """Return the half-axis (or origin) on which (re, im) point is located. """
    assert not re or not im

    if not re and not im:
        return OO

    if not re:
        if im > 0:
            return A2
        else:
            return A4
    else:
        if re > 0:
            return A1
        else:
            return A3


def _intervals_to_quadrants(intervals, f1, f2, s, t, F):
    """Generate a sequence of extended quadrants from a list of critical points. """
    if not intervals:
        return []

    Q = []

    if not f1:
        (a, b), _, _ = intervals[0]

        if a == b == s:
            if len(intervals) == 1:
                if dmp_eval_in(f2, t, 0, 0, F) > 0:
                    return [OO, A2]
                else:
                    return [OO, A4]
            else:
                (a, _), _, _ = intervals[1]

                if dmp_eval_in(f2, (s + a)/2, 0, 0, F) > 0:
                    Q.extend([OO, A2])
                    f2_sgn = +1
                else:
                    Q.extend([OO, A4])
                    f2_sgn = -1

                intervals = intervals[1:]
        else:
            if dmp_eval_in(f2, s, 0, 0, F) > 0:
                Q.append(A2)
                f2_sgn = +1
            else:
                Q.append(A4)
                f2_sgn = -1

        for (a, _), indices, _ in intervals:
            Q.append(OO)

            assert indices[1] % 2 == 1
            f2_sgn = -f2_sgn

            if a != t:
                if f2_sgn > 0:
                    Q.append(A2)
                else:
                    Q.append(A4)

        return Q

    if not f2:
        (a, b), _, _ = intervals[0]

        if a == b == s:
            if len(intervals) == 1:
                if dmp_eval_in(f1, t, 0, 0, F) > 0:
                    return [OO, A1]
                else:
                    return [OO, A3]
            else:
                (a, _), _, _ = intervals[1]

                if dmp_eval_in(f1, (s + a)/2, 0, 0, F) > 0:
                    Q.extend([OO, A1])
                    f1_sgn = +1
                else:
                    Q.extend([OO, A3])
                    f1_sgn = -1

                intervals = intervals[1:]
        else:
            if dmp_eval_in(f1, s, 0, 0, F) > 0:
                Q.append(A1)
                f1_sgn = +1
            else:
                Q.append(A3)
                f1_sgn = -1

        for (a, _), indices, _ in intervals:
            Q.append(OO)

            assert indices[0] % 2 == 1
            f1_sgn = -f1_sgn

            if a != t:
                if f1_sgn > 0:
                    Q.append(A1)
                else:
                    Q.append(A3)

        return Q

    re = dmp_eval_in(f1, s, 0, 0, F)
    im = dmp_eval_in(f2, s, 0, 0, F)

    if not re or not im:
        Q.append(_classify_point(re, im))

        if len(intervals) == 1:
            re = dmp_eval_in(f1, t, 0, 0, F)
            im = dmp_eval_in(f2, t, 0, 0, F)
        else:
            (a, _), _, _ = intervals[1]

            re = dmp_eval_in(f1, (s + a)/2, 0, 0, F)
            im = dmp_eval_in(f2, (s + a)/2, 0, 0, F)

        intervals = intervals[1:]

    if re > 0:
        f1_sgn = +1
    else:
        f1_sgn = -1

    if im > 0:
        f2_sgn = +1
    else:
        f2_sgn = -1

    sgn = {
        (+1, +1): Q1,
        (-1, +1): Q2,
        (-1, -1): Q3,
        (+1, -1): Q4,
    }

    Q.append(sgn[(f1_sgn, f2_sgn)])

    for (a, b), indices, _ in intervals:
        if a == b:
            re = dmp_eval_in(f1, a, 0, 0, F)
            im = dmp_eval_in(f2, a, 0, 0, F)

            Q.append(_classify_point(re, im))

        if 0 in indices:
            if indices[0] % 2 == 1:
                f1_sgn = -f1_sgn

        if 1 in indices:
            if indices[1] % 2 == 1:
                f2_sgn = -f2_sgn

        if not (a == b and b == t):
            Q.append(sgn[(f1_sgn, f2_sgn)])

    return Q


def _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4, exclude=None):
    """Transform sequences of quadrants to a sequence of rules. """
    if exclude is True:
        edges = [1, 1, 0, 0]

        corners = {
            (0, 1): 1,
            (1, 2): 1,
            (2, 3): 0,
            (3, 0): 1,
        }
    else:
        edges = [0, 0, 0, 0]

        corners = {
            (0, 1): 0,
            (1, 2): 0,
            (2, 3): 0,
            (3, 0): 0,
        }

    if exclude is not None and exclude is not True:
        exclude = set(exclude)

        for i, edge in enumerate(['S', 'E', 'N', 'W']):
            if edge in exclude:
                edges[i] = 1

        for i, corner in enumerate(['SW', 'SE', 'NE', 'NW']):
            if corner in exclude:
                corners[((i - 1) % 4, i)] = 1

    QQ, rules = [Q_L1, Q_L2, Q_L3, Q_L4], []

    for i, Q in enumerate(QQ):
        if not Q:
            continue

        if Q[-1] == OO:
            Q = Q[:-1]

        if Q[0] == OO:
            j, Q = (i - 1) % 4, Q[1:]
            qq = QQ[j][-2], OO, Q[0]
            rules.append((_rules_ambiguous[qq], corners[(j, i)]))

        q1, k = Q[0], 1

        while k < len(Q):
            q2, k = Q[k], k + 1

            if q2 != OO:
                qq = q1, q2
                rules.append((_rules_simple[qq], 0))
            else:
                qq, k = (q1, q2, Q[k]), k + 1
                rules.append((_rules_ambiguous[qq], edges[i]))

            q1 = qq[-1]

    return rules


def _reverse_intervals(intervals):
    """Reverse intervals for traversal from right to left and from top to bottom. """
    return [((b, a), indices, f) for (a, b), indices, f in reversed(intervals)]


def _winding_number(T, field):
    """Compute the winding number of the input polynomial, i.e. the number of roots. """
    return int(sum((field(_values[t][i][0])/field(_values[t][i][1]) for t, i in T), field(0)) / field(2))


def _roots_bound(f, F):
    lc = F.to_expr(dmp_LC(f, F))
    B = 2*max(abs(F.to_expr(c)/lc) for c in f)
    if not F.is_AlgebraicField:
        return F.convert(B)
    else:
        return F.domain(int(100*B) + 1)/F.domain(100)


def _count_roots(f1, f2, F, inf, sup, exclude=None):
    (u, v), (s, t) = inf, sup

    f1L1F = dmp_eval_in(f1, v, 1, 1, F)
    f2L1F = dmp_eval_in(f2, v, 1, 1, F)

    f1L2F = dmp_eval_in(f1, s, 0, 1, F)
    f2L2F = dmp_eval_in(f2, s, 0, 1, F)

    f1L3F = dmp_eval_in(f1, t, 1, 1, F)
    f2L3F = dmp_eval_in(f2, t, 1, 1, F)

    f1L4F = dmp_eval_in(f1, u, 0, 1, F)
    f2L4F = dmp_eval_in(f2, u, 0, 1, F)

    S_L1 = [f1L1F, f2L1F]
    S_L2 = [f1L2F, f2L2F]
    S_L3 = [f1L3F, f2L3F]
    S_L4 = [f1L4F, f2L4F]

    I_L1 = dup_isolate_real_roots_list(S_L1, F, inf=u, sup=s, fast=True, basis=True, strict=True)
    I_L2 = dup_isolate_real_roots_list(S_L2, F, inf=v, sup=t, fast=True, basis=True, strict=True)
    I_L3 = dup_isolate_real_roots_list(S_L3, F, inf=u, sup=s, fast=True, basis=True, strict=True)
    I_L4 = dup_isolate_real_roots_list(S_L4, F, inf=v, sup=t, fast=True, basis=True, strict=True)

    I_L3 = _reverse_intervals(I_L3)
    I_L4 = _reverse_intervals(I_L4)

    Q_L1 = _intervals_to_quadrants(I_L1, f1L1F, f2L1F, u, s, F)
    Q_L2 = _intervals_to_quadrants(I_L2, f1L2F, f2L2F, v, t, F)
    Q_L3 = _intervals_to_quadrants(I_L3, f1L3F, f2L3F, s, u, F)
    Q_L4 = _intervals_to_quadrants(I_L4, f1L4F, f2L4F, t, v, F)

    T = _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4, exclude=exclude)

    return _winding_number(T, F)


def dup_count_complex_roots(f, K, inf=None, sup=None, exclude=None):
    """Count all roots in [u + v*I, s + t*I] rectangle using Collins-Krandick algorithm. """
    F = K.field

    f = dmp_convert(f, 0, K, F)

    if not all(isinstance(_, tuple) for _ in (inf, sup)):
        B = _roots_bound(f, F)

    if isinstance(inf, tuple):
        u, v = inf
    elif inf is not None:
        u, v = inf, -B
    else:
        u, v = -B, -B

    if isinstance(sup, tuple):
        s, t = sup
    elif sup is not None:
        s, t = sup, B
    else:
        s, t = B, B

    f1, f2 = dup_real_imag(f, F)

    if F.is_ComplexAlgebraicField:
        F = F.domain

    return _count_roots(f1, f2, F, (u, v), (s, t), exclude)


def _vertical_bisection(N, a, b, I, Q, F1, F2, f1, f2, F):
    """Vertical bisection step in Collins-Krandick root isolation algorithm. """
    (u, v), (s, t) = a, b

    I_L1, I_L2, I_L3, I_L4 = I
    Q_L1, Q_L2, Q_L3, Q_L4 = Q

    f1L1F, f1L2F, f1L3F, f1L4F = F1
    f2L1F, f2L2F, f2L3F, f2L4F = F2

    x = (u + s) / 2

    f1V = dmp_eval_in(f1, x, 0, 1, F)
    f2V = dmp_eval_in(f2, x, 0, 1, F)

    I_V = dup_isolate_real_roots_list([f1V, f2V], F, inf=v, sup=t, fast=True, strict=True, basis=True)

    I_L1_L, I_L1_R = [], []
    I_L2_L, I_L2_R = I_V, I_L2
    I_L3_L, I_L3_R = [], []
    I_L4_L, I_L4_R = I_L4, _reverse_intervals(I_V)

    for I in I_L1:
        (a, b), indices, h = I

        if a == b:
            if a == x:
                I_L1_L.append(I)
                I_L1_R.append(I)
            elif a < x:
                I_L1_L.append(I)
            else:
                I_L1_R.append(I)
        else:
            if b <= x:
                I_L1_L.append(I)
            elif a >= x:
                I_L1_R.append(I)
            else:
                a, b = dup_refine_real_root(h, a, b, F, disjoint=x, fast=True)

                if b <= x:
                    I_L1_L.append(((a, b), indices, h))
                if a >= x:
                    I_L1_R.append(((a, b), indices, h))

    for I in I_L3:
        (b, a), indices, h = I

        if a == b:
            if a == x:
                I_L3_L.append(I)
                I_L3_R.append(I)
            elif a < x:
                I_L3_L.append(I)
            else:
                I_L3_R.append(I)
        else:
            if b <= x:
                I_L3_L.append(I)
            elif a >= x:
                I_L3_R.append(I)
            else:
                a, b = dup_refine_real_root(h, a, b, F, disjoint=x, fast=True)

                if b <= x:
                    I_L3_L.append(((b, a), indices, h))
                if a >= x:
                    I_L3_R.append(((b, a), indices, h))

    Q_L1_L = _intervals_to_quadrants(I_L1_L, f1L1F, f2L1F, u, x, F)
    Q_L2_L = _intervals_to_quadrants(I_L2_L, f1V, f2V, v, t, F)
    Q_L3_L = _intervals_to_quadrants(I_L3_L, f1L3F, f2L3F, x, u, F)
    Q_L4_L = Q_L4

    Q_L1_R = _intervals_to_quadrants(I_L1_R, f1L1F, f2L1F, x, s, F)
    Q_L2_R = Q_L2
    Q_L3_R = _intervals_to_quadrants(I_L3_R, f1L3F, f2L3F, s, x, F)
    Q_L4_R = _intervals_to_quadrants(I_L4_R, f1V, f2V, t, v, F)

    T_L = _traverse_quadrants(Q_L1_L, Q_L2_L, Q_L3_L, Q_L4_L, exclude=True)
    T_R = _traverse_quadrants(Q_L1_R, Q_L2_R, Q_L3_R, Q_L4_R, exclude=True)

    N_L = _winding_number(T_L, F)
    N_R = _winding_number(T_R, F)

    I_L = I_L1_L, I_L2_L, I_L3_L, I_L4_L
    Q_L = Q_L1_L, Q_L2_L, Q_L3_L, Q_L4_L

    I_R = I_L1_R, I_L2_R, I_L3_R, I_L4_R
    Q_R = Q_L1_R, Q_L2_R, Q_L3_R, Q_L4_R

    F1_L = f1L1F, f1V, f1L3F, f1L4F
    F2_L = f2L1F, f2V, f2L3F, f2L4F

    F1_R = f1L1F, f1L2F, f1L3F, f1V
    F2_R = f2L1F, f2L2F, f2L3F, f2V

    a, b = (u, v), (x, t)
    c, d = (x, v), (s, t)

    D_L = N_L, a, b, I_L, Q_L, F1_L, F2_L
    D_R = N_R, c, d, I_R, Q_R, F1_R, F2_R

    return D_L, D_R


def _horizontal_bisection(N, a, b, I, Q, F1, F2, f1, f2, F):
    """Horizontal bisection step in Collins-Krandick root isolation algorithm. """
    (u, v), (s, t) = a, b

    I_L1, I_L2, I_L3, I_L4 = I
    Q_L1, Q_L2, Q_L3, Q_L4 = Q

    f1L1F, f1L2F, f1L3F, f1L4F = F1
    f2L1F, f2L2F, f2L3F, f2L4F = F2

    y = (v + t) / 2

    f1H = dmp_eval_in(f1, y, 1, 1, F)
    f2H = dmp_eval_in(f2, y, 1, 1, F)

    I_H = dup_isolate_real_roots_list([f1H, f2H], F, inf=u, sup=s, fast=True, strict=True, basis=True)

    I_L1_B, I_L1_U = I_L1, I_H
    I_L2_B, I_L2_U = [], []
    I_L3_B, I_L3_U = _reverse_intervals(I_H), I_L3
    I_L4_B, I_L4_U = [], []

    for I in I_L2:
        (a, b), indices, h = I

        if a == b:
            if a == y:
                I_L2_B.append(I)
                I_L2_U.append(I)
            elif a < y:
                I_L2_B.append(I)
            else:
                I_L2_U.append(I)
        else:
            if b <= y:
                I_L2_B.append(I)
            elif a >= y:
                I_L2_U.append(I)
            else:
                a, b = dup_refine_real_root(h, a, b, F, disjoint=y, fast=True)

                if b <= y:
                    I_L2_B.append(((a, b), indices, h))
                if a >= y:
                    I_L2_U.append(((a, b), indices, h))

    for I in I_L4:
        (b, a), indices, h = I

        if a == b:
            if a == y:
                I_L4_B.append(I)
                I_L4_U.append(I)
            elif a < y:
                I_L4_B.append(I)
            else:
                I_L4_U.append(I)
        else:
            if b <= y:
                I_L4_B.append(I)
            elif a >= y:
                I_L4_U.append(I)
            else:
                a, b = dup_refine_real_root(h, a, b, F, disjoint=y, fast=True)

                if b <= y:
                    I_L4_B.append(((b, a), indices, h))
                if a >= y:
                    I_L4_U.append(((b, a), indices, h))

    Q_L1_B = Q_L1
    Q_L2_B = _intervals_to_quadrants(I_L2_B, f1L2F, f2L2F, v, y, F)
    Q_L3_B = _intervals_to_quadrants(I_L3_B, f1H, f2H, s, u, F)
    Q_L4_B = _intervals_to_quadrants(I_L4_B, f1L4F, f2L4F, y, v, F)

    Q_L1_U = _intervals_to_quadrants(I_L1_U, f1H, f2H, u, s, F)
    Q_L2_U = _intervals_to_quadrants(I_L2_U, f1L2F, f2L2F, y, t, F)
    Q_L3_U = Q_L3
    Q_L4_U = _intervals_to_quadrants(I_L4_U, f1L4F, f2L4F, t, y, F)

    T_B = _traverse_quadrants(Q_L1_B, Q_L2_B, Q_L3_B, Q_L4_B, exclude=True)
    T_U = _traverse_quadrants(Q_L1_U, Q_L2_U, Q_L3_U, Q_L4_U, exclude=True)

    N_B = _winding_number(T_B, F)
    N_U = _winding_number(T_U, F)

    I_B = I_L1_B, I_L2_B, I_L3_B, I_L4_B
    Q_B = Q_L1_B, Q_L2_B, Q_L3_B, Q_L4_B

    I_U = I_L1_U, I_L2_U, I_L3_U, I_L4_U
    Q_U = Q_L1_U, Q_L2_U, Q_L3_U, Q_L4_U

    F1_B = f1L1F, f1L2F, f1H, f1L4F
    F2_B = f2L1F, f2L2F, f2H, f2L4F

    F1_U = f1H, f1L2F, f1L3F, f1L4F
    F2_U = f2H, f2L2F, f2L3F, f2L4F

    a, b = (u, v), (s, y)
    c, d = (u, y), (s, t)

    D_B = N_B, a, b, I_B, Q_B, F1_B, F2_B
    D_U = N_U, c, d, I_U, Q_U, F1_U, F2_U

    return D_B, D_U


def _depth_first_select(rectangles):
    """Find a rectangle of minimum area for bisection. """
    min_area, j = None, None

    for i, (_, (u, v), (s, t), _, _, _, _) in enumerate(rectangles):
        area = (s - u)*(t - v)

        if min_area is None or area < min_area:
            min_area, j = area, i

    return rectangles.pop(j)


def _rectangle_small_p(a, b, eps):
    """Return ``True`` if the given rectangle is small enough. """
    (u, v), (s, t) = a, b

    if eps is not None:
        return s - u < eps and t - v < eps
    else:
        return True


def dup_isolate_complex_roots_sqf(f, K, eps=None, inf=None, sup=None, blackbox=False):
    """Isolate complex roots of a square-free polynomial using Collins-Krandick algorithm. """
    if dmp_degree_in(f, 0, 0) <= 0:
        return []

    F = K.field

    f = dmp_convert(f, 0, K, F)

    if not all(isinstance(_, tuple) for _ in (inf, sup)):
        B = _roots_bound(f, F)

    if isinstance(inf, tuple):
        u, v = inf
    elif inf is not None:
        u, v = inf, -B
    else:
        u, v = -B, -B

    if isinstance(sup, tuple):
        s, t = sup
    elif sup is not None:
        s, t = sup, B
    else:
        s, t = B, B

    if t <= v or s <= u:
        raise ValueError("not a valid complex isolation rectangle")

    if v < 0 < t:
        roots = dup_isolate_complex_roots_sqf(f, F, eps=eps, inf=(u, 0),
                                              sup=(s, t), blackbox=True)
        if F.is_RationalField or F.is_RealAlgebraicField:
            _roots = []
            for root in roots:
                croot = root.conjugate()
                if croot.ay >= v:
                    _roots.append(croot)
                _roots.append(root)
            roots = _roots
        elif F.is_ComplexAlgebraicField:
            # Take conjugated polynomial to get solutions in the
            # bottom half-plane.
            f = [_.conjugate() for _ in f]
            roots += [_.conjugate()
                      for _ in dup_isolate_complex_roots_sqf(f, F, eps=eps,
                                                             inf=(u, 0), sup=(s, -v),
                                                             blackbox=True)]
            roots = sorted(roots, key=lambda r: (r.ax, r.ay))
        else:
            raise NotImplementedError

        return roots if blackbox else [r.as_tuple() for r in roots]

    f1, f2 = dup_real_imag(f, F)

    if F.is_ComplexAlgebraicField:
        F = F.domain

    f1L1 = dmp_eval_in(f1, v, 1, 1, F)
    f2L1 = dmp_eval_in(f2, v, 1, 1, F)

    f1L2 = dmp_eval_in(f1, s, 0, 1, F)
    f2L2 = dmp_eval_in(f2, s, 0, 1, F)

    f1L3 = dmp_eval_in(f1, t, 1, 1, F)
    f2L3 = dmp_eval_in(f2, t, 1, 1, F)

    f1L4 = dmp_eval_in(f1, u, 0, 1, F)
    f2L4 = dmp_eval_in(f2, u, 0, 1, F)

    S_L1 = [f1L1, f2L1]
    S_L2 = [f1L2, f2L2]
    S_L3 = [f1L3, f2L3]
    S_L4 = [f1L4, f2L4]

    I_L1 = dup_isolate_real_roots_list(S_L1, F, inf=u, sup=s, fast=True, strict=True, basis=True)
    I_L2 = dup_isolate_real_roots_list(S_L2, F, inf=v, sup=t, fast=True, strict=True, basis=True)
    I_L3 = dup_isolate_real_roots_list(S_L3, F, inf=u, sup=s, fast=True, strict=True, basis=True)
    I_L4 = dup_isolate_real_roots_list(S_L4, F, inf=v, sup=t, fast=True, strict=True, basis=True)

    I_L3 = _reverse_intervals(I_L3)
    I_L4 = _reverse_intervals(I_L4)

    Q_L1 = _intervals_to_quadrants(I_L1, f1L1, f2L1, u, s, F)
    Q_L2 = _intervals_to_quadrants(I_L2, f1L2, f2L2, v, t, F)
    Q_L3 = _intervals_to_quadrants(I_L3, f1L3, f2L3, s, u, F)
    Q_L4 = _intervals_to_quadrants(I_L4, f1L4, f2L4, t, v, F)

    T = _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4)
    N = _winding_number(T, F)

    if not N:
        return []

    I_L = I_L1, I_L2, I_L3, I_L4
    Q_L = Q_L1, Q_L2, Q_L3, Q_L4

    F1 = f1L1, f1L2, f1L3, f1L4
    F2 = f2L1, f2L2, f2L3, f2L4

    rectangles, roots = [(N, (u, v), (s, t), I_L, Q_L, F1, F2)], []

    while rectangles:
        N, (u, v), (s, t), I_L, Q_L, F1, F2 = _depth_first_select(rectangles)

        if s - u > t - v:
            D_L, D_R = _vertical_bisection(N, (u, v), (s, t), I_L, Q_L, F1, F2, f1, f2, F)

            N_L, a, b, I_L, Q_L, F1_L, F2_L = D_L
            N_R, c, d, I_R, Q_R, F1_R, F2_R = D_R

            if N_L >= 1:
                if N_L == 1 and _rectangle_small_p(a, b, eps):
                    roots.append(ComplexInterval(a, b, I_L, Q_L, F1_L, F2_L, f1, f2, F))
                else:
                    rectangles.append(D_L)

            if N_R >= 1:
                if N_R == 1 and _rectangle_small_p(c, d, eps):
                    roots.append(ComplexInterval(c, d, I_R, Q_R, F1_R, F2_R, f1, f2, F))
                else:
                    rectangles.append(D_R)
        else:
            D_B, D_U = _horizontal_bisection(N, (u, v), (s, t), I_L, Q_L, F1, F2, f1, f2, F)

            N_B, a, b, I_B, Q_B, F1_B, F2_B = D_B
            N_U, c, d, I_U, Q_U, F1_U, F2_U = D_U

            if N_B >= 1:
                if N_B == 1 and _rectangle_small_p(a, b, eps):
                    roots.append(ComplexInterval(a, b, I_B, Q_B, F1_B, F2_B, f1, f2, F))
                else:
                    rectangles.append(D_B)

            if N_U >= 1:
                if N_U == 1 and _rectangle_small_p(c, d, eps):
                    roots.append(ComplexInterval(c, d, I_U, Q_U, F1_U, F2_U, f1, f2, F))
                else:
                    rectangles.append(D_U)

    roots = sorted(roots, key=lambda r: (r.ax, r.ay))
    return roots if blackbox else [r.as_tuple() for r in roots]


def dup_isolate_all_roots_sqf(f, K, eps=None, inf=None, sup=None, fast=False, blackbox=False):
    """Isolate real and complex roots of a square-free polynomial ``f``. """
    return (dup_isolate_real_roots_sqf(f, K, eps=eps, inf=inf, sup=sup, fast=fast, blackbox=blackbox),
            dup_isolate_complex_roots_sqf(f, K, eps=eps, inf=inf, sup=sup, blackbox=blackbox))


def dup_isolate_all_roots(f, K, eps=None, inf=None, sup=None, fast=False):
    """Isolate real and complex roots of a non-square-free polynomial ``f``. """
    if not K.is_IntegerRing and not K.is_RationalField:
        raise DomainError("isolation of real and complex roots is not supported over %s" % K)

    _, factors = dmp_sqf_list(f, 0, K)

    if len(factors) == 1:
        (f, k), = factors

        real_part, complex_part = dup_isolate_all_roots_sqf(f, K, eps=eps, inf=inf, sup=sup, fast=fast)

        real_part = [((a, b), k) for (a, b) in real_part]
        complex_part = [((a, b), k) for (a, b) in complex_part]

        return real_part, complex_part
    else:
        raise NotImplementedError("only trivial square-free polynomials are supported")


class RealInterval:
    """A fully qualified representation of a real isolation interval. """

    def __init__(self, data, f, dom):
        """Initialize new real interval with complete information. """
        if len(data) == 2:
            s, t = data

            self.neg = False

            if s < 0:
                if t <= 0:
                    f, s, t, self.neg = dup_mirror(f, dom), -t, -s, True
                else:
                    raise ValueError("can't refine a real root in (%s, %s)" % (s, t))

            a, b, c, d = _mobius_from_interval((s, t), dom.field)

            f = dup_transform(f, dmp_strip([a, b], 0), dmp_strip([c, d], 0), dom)

            self.mobius = a, b, c, d
        else:
            self.mobius = data[:-1]
            self.neg = data[-1]

        self.f, self.domain = f, dom

    @property
    def a(self):
        """Return the position of the left end. """
        a, b, c, d = self.mobius

        if not self.neg:
            if a*d < b*c:
                return a/c
            return b/d
        else:
            if a*d > b*c:
                return -a/c
            return -b/d

    @property
    def b(self):
        """Return the position of the right end. """
        was = self.neg
        self.neg = not was
        rv = -self.a
        self.neg = was
        return rv

    @property
    def center(self):
        """Return the center of the real isolating interval. """
        return (self.a + self.b)/2

    def as_tuple(self):
        """Return tuple representation of real isolating interval. """
        return self.a, self.b

    def is_disjoint(self, other):
        """Return ``True`` if two isolation intervals are disjoint. """
        return self.b <= other.a or other.b <= self.a

    def refine(self):
        """Perform one step of real root refinement algorithm. """
        assert self.mobius is not None

        f, mobius = dup_inner_refine_real_root(self.f, self.mobius, self.domain, steps=1, mobius=True)

        return RealInterval(mobius + (self.neg,), f, self.domain)


class ComplexInterval:
    """A fully qualified representation of a complex isolation interval.
    The printed form is shown as (x1, y1) x (x2, y2): the southwest x northeast
    coordinates of the interval's rectangle.
    """

    def __init__(self, a, b, I, Q, F1, F2, f1, f2, dom, conj=False):
        """Initialize new complex interval with complete information. """
        self.a, self.b = a, b  # the southwest and northeast corner: (x1, y1), (x2, y2)
        self.I, self.Q = I, Q

        self.f1, self.F1 = f1, F1
        self.f2, self.F2 = f2, F2

        self.domain = dom
        self.conj = conj

    @property
    def ax(self):
        """Return ``x`` coordinate of south-western corner. """
        return self.a[0]

    @property
    def ay(self):
        """Return ``y`` coordinate of south-western corner. """
        if not self.conj:
            return +self.a[1]
        else:
            return -self.b[1]

    @property
    def bx(self):
        """Return ``x`` coordinate of north-eastern corner. """
        return self.b[0]

    @property
    def by(self):
        """Return ``y`` coordinate of north-eastern corner. """
        if not self.conj:
            return +self.b[1]
        else:
            return -self.a[1]

    @property
    def center(self):
        """Return the center of the complex isolating interval. """
        return (self.ax + self.bx)/2, (self.ay + self.by)/2

    def as_tuple(self):
        """Return tuple representation of complex isolating interval. """
        return (self.ax, self.ay), (self.bx, self.by)

    def conjugate(self):
        """Return conjugated isolating interval. """
        return ComplexInterval(self.a, self.b, self.I, self.Q,
                               self.F1, self.F2, self.f1, self.f2,
                               self.domain, not self.conj)

    def is_disjoint(self, other, check_re_refinement=False, re_disjoint=False):
        """Return ``True`` if two isolation intervals are disjoint.

        Parameters
        ==========

        check_re_refinement : bool, optional
            If enabled, test that either real projections of isolation
            intervals are disjoint or roots share common real part.
        re_disjoint : bool, optional
            If enabled, return ``True`` only if real projections of isolation
            intervals are disjoint.
        """
        test_re = self.bx <= other.ax or other.bx <= self.ax
        if test_re or re_disjoint:
            return test_re
        if not check_re_refinement:
            return self.by <= other.ay or other.by <= self.ay

        resultants = []
        for i in (self, other):
            re = dmp_permute(i.f1, [1, 0], 1, i.domain)
            im = dmp_permute(i.f2, [1, 0], 1, i.domain)
            re, im = map(lambda x: dmp_to_tuple(x, 1), [re, im])
            resultants.append(dmp_resultant(re, im, 1, i.domain))
        dom = self.domain.unify(other.domain)
        gcd = dmp_gcd(*resultants, 0, dom)
        gcd_roots = dup_isolate_real_roots(gcd, dom,
                                           inf=max(self.ax, other.ax),
                                           sup=min(self.bx, other.bx))
        if len(gcd_roots) != 1:
            return False

        l, r = gcd_roots[0][0]
        if l == r:
            # Can't use dup_count_complex_roots() here, as the new isolation
            # interval will be degenerate: a vertical line segment.  Make
            # sure we only count roots on the northern/western edges and
            # on the north-western corner of the original isolation rectangle.
            for i in (self, other):
                dom = i.domain.algebraic_field(I)
                f1 = dmp_convert(dmp_eval_in(i.f1, l, 0, 1, i.domain), 0, i.domain, dom)
                f2 = dmp_convert(dmp_eval_in(i.f2, l, 0, 1, i.domain), 0, i.domain, dom)
                f = dmp_add(f1, dmp_mul_ground(f2, dom.unit, 0, dom), 0, dom)
                f = dmp_compose(f, [-dom.unit, 0], 0, dom)
                if i.conj:
                    f = [_.conjugate() for _ in f]
                f = dmp_compose(f, [dom.unit, 0], 0, dom)
                r = dup_isolate_real_roots(f, dom, inf=i.ay, sup=i.by)
                if len([1 for _ in r if ((not i.conj and i.ay < _[0][0]) or
                                         (i.conj and _[0][1] < i.by)) and l < i.bx]) != 1:
                    return False
            else:
                return True
        else:
            return all(_count_roots(i.f1, i.f2, i.domain, (l, i.ay), (r, i.by)) == 1
                       for i in (self, other))

    def refine(self):
        """Perform one step of complex root refinement algorithm. """
        (u, v), (s, t) = self.a, self.b

        I, Q = self.I, self.Q

        f1, F1 = self.f1, self.F1
        f2, F2 = self.f2, self.F2

        dom = self.domain

        if s - u > t - v:
            D_L, D_R = _vertical_bisection(1, (u, v), (s, t), I, Q, F1, F2, f1, f2, dom)

            if D_L[0] == 1:
                _, a, b, I, Q, F1, F2 = D_L
            else:
                _, a, b, I, Q, F1, F2 = D_R
        else:
            D_B, D_U = _horizontal_bisection(1, (u, v), (s, t), I, Q, F1, F2, f1, f2, dom)

            if D_B[0] == 1:
                _, a, b, I, Q, F1, F2 = D_B
            else:
                _, a, b, I, Q, F1, F2 = D_U

        return ComplexInterval(a, b, I, Q, F1, F2, f1, f2, dom, self.conj)
