"""Real and complex root isolation and refinement algorithms."""

import collections
import math
import operator

from ..core import Dummy, I
from .polyerrors import DomainError, RefinementFailedError


def _mobius_from_interval(I, field):
    """Convert an open interval to a Mobius transform."""
    s, t = I

    a, c = map(field, (s.numerator, s.denominator))
    b, d = map(field, (t.numerator, t.denominator))

    return a, b, c, d


def _mobius_to_interval(M):
    """Convert a Mobius transform to an open interval."""
    a, b, c, d = M

    s, t = a/c, b/d

    return (s, t) if s <= t else (t, s)


def _disjoint_p(M, N, strict=False):
    """Check if Mobius transforms define disjoint intervals."""
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
    return a2*d1 > c2*b1 or b2*c1 < d2*a1


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
    """Return the half-axis (or origin) on which (re, im) point is located."""
    assert not re or not im

    if not re and not im:
        return OO

    if not re:
        if im > 0:
            return A2
        return A4
    if re > 0:
        return A1
    return A3


def _intervals_to_quadrants(intervals, f1, f2, s, t):
    """Generate a sequence of extended quadrants from a list of critical points."""
    if not intervals:
        return []

    Q = []

    if not f1:
        (a, b), _, _ = intervals[0]

        if a == b == s:
            if len(intervals) == 1:
                if f2.eval(0, t) > 0:
                    return [OO, A2]
                return [OO, A4]
            (a, _), _, _ = intervals[1]

            if f2.eval(0, (s + a)/2) > 0:
                Q.extend([OO, A2])
                f2_sgn = +1
            else:
                Q.extend([OO, A4])
                f2_sgn = -1

            intervals = intervals[1:]
        else:
            if f2.eval(0, s) > 0:
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
                if f1.eval(0, t) > 0:
                    return [OO, A1]
                return [OO, A3]
            (a, _), _, _ = intervals[1]

            if f1.eval(0, (s + a)/2) > 0:
                Q.extend([OO, A1])
                f1_sgn = +1
            else:
                Q.extend([OO, A3])
                f1_sgn = -1

            intervals = intervals[1:]
        else:
            if f1.eval(0, s) > 0:
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

    re = f1.eval(0, s)
    im = f2.eval(0, s)

    if not re or not im:
        Q.append(_classify_point(re, im))

        if len(intervals) == 1:
            re = f1.eval(0, t)
            im = f2.eval(0, t)
        else:
            (a, _), _, _ = intervals[1]

            re = f1.eval(0, (s + a)/2)
            im = f2.eval(0, (s + a)/2)

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
            re = f1.eval(0, a)
            im = f2.eval(0, a)

            Q.append(_classify_point(re, im))

        if 0 in indices:
            if indices[0] % 2 == 1:
                f1_sgn = -f1_sgn

        if 1 in indices:
            if indices[1] % 2 == 1:
                f2_sgn = -f2_sgn

        if a != b or b != t:
            Q.append(sgn[(f1_sgn, f2_sgn)])

    return Q


def _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4, exclude=None):
    """Transform sequences of quadrants to a sequence of rules."""
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
    """Reverse intervals for traversal from right to left and from top to bottom."""
    return [((b, a), indices, f) for (a, b), indices, f in reversed(intervals)]


def _winding_number(T, field):
    """Compute the winding number of the input polynomial, i.e. the number of roots."""
    return int(sum((field(_values[t][i][0])/field(_values[t][i][1]) for t, i in T), field(0)) / field(2))


def _get_rectangle(f1, f2, inf, sup, exclude=None):
    (u, v), (s, t) = inf, sup

    f1L1 = f1.eval(1, v)
    f2L1 = f2.eval(1, v)

    f1L2 = f1.eval(0, s)
    f2L2 = f2.eval(0, s)

    f1L3 = f1.eval(1, t)
    f2L3 = f2.eval(1, t)

    f1L4 = f1.eval(0, u)
    f2L4 = f2.eval(0, u)

    I_L1 = f1L1.ring._isolate_real_roots_pair(f1L1, f2L1, inf=u, sup=s, basis=True, strict=True)
    I_L2 = f1L2.ring._isolate_real_roots_pair(f1L2, f2L2, inf=v, sup=t, basis=True, strict=True)
    I_L3 = f1L3.ring._isolate_real_roots_pair(f1L3, f2L3, inf=u, sup=s, basis=True, strict=True)
    I_L4 = f1L4.ring._isolate_real_roots_pair(f1L4, f2L4, inf=v, sup=t, basis=True, strict=True)

    I_L3 = _reverse_intervals(I_L3)
    I_L4 = _reverse_intervals(I_L4)

    Q_L1 = _intervals_to_quadrants(I_L1, f1L1, f2L1, u, s)
    Q_L2 = _intervals_to_quadrants(I_L2, f1L2, f2L2, v, t)
    Q_L3 = _intervals_to_quadrants(I_L3, f1L3, f2L3, s, u)
    Q_L4 = _intervals_to_quadrants(I_L4, f1L4, f2L4, t, v)

    T = _traverse_quadrants(Q_L1, Q_L2, Q_L3, Q_L4, exclude=exclude)

    N = _winding_number(T, f1L1.ring.domain)

    I_L = I_L1, I_L2, I_L3, I_L4
    Q_L = Q_L1, Q_L2, Q_L3, Q_L4

    F1 = f1L1, f1L2, f1L3, f1L4
    F2 = f2L1, f2L2, f2L3, f2L4

    return N, inf, sup, I_L, Q_L, F1, F2


def _vertical_bisection(N, a, b, I, Q, F1, F2, f1, f2):
    """Vertical bisection step in Collins-Krandick root isolation algorithm."""
    (u, v), (s, t) = a, b

    I_L1, I_L2, I_L3, I_L4 = I
    _, Q_L2, _, Q_L4 = Q

    f1L1F, f1L2F, f1L3F, f1L4F = F1
    f2L1F, f2L2F, f2L3F, f2L4F = F2
    F = f1L1F.ring.domain

    x = (u + s)/F(2)

    f1V = f1.eval(0, x)
    f2V = f2.eval(0, x)

    I_V = f1V.ring._isolate_real_roots_pair(f1V, f2V, inf=v, sup=t, strict=True, basis=True)

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
                a, b = h.ring._refine_real_root(h, a, b, disjoint=x)

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
                a, b = h.ring._refine_real_root(h, a, b, disjoint=x)

                if b <= x:
                    I_L3_L.append(((b, a), indices, h))
                if a >= x:
                    I_L3_R.append(((b, a), indices, h))

    Q_L1_L = _intervals_to_quadrants(I_L1_L, f1L1F, f2L1F, u, x)
    Q_L2_L = _intervals_to_quadrants(I_L2_L, f1V, f2V, v, t)
    Q_L3_L = _intervals_to_quadrants(I_L3_L, f1L3F, f2L3F, x, u)
    Q_L4_L = Q_L4

    Q_L1_R = _intervals_to_quadrants(I_L1_R, f1L1F, f2L1F, x, s)
    Q_L2_R = Q_L2
    Q_L3_R = _intervals_to_quadrants(I_L3_R, f1L3F, f2L3F, s, x)
    Q_L4_R = _intervals_to_quadrants(I_L4_R, f1V, f2V, t, v)

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


def _horizontal_bisection(N, a, b, I, Q, F1, F2, f1, f2):
    """Horizontal bisection step in Collins-Krandick root isolation algorithm."""
    (u, v), (s, t) = a, b

    I_L1, I_L2, I_L3, I_L4 = I
    Q_L1, _, Q_L3, _ = Q

    f1L1F, f1L2F, f1L3F, f1L4F = F1
    f2L1F, f2L2F, f2L3F, f2L4F = F2
    F = f1L1F.ring.domain

    y = (v + t)/F(2)

    f1H = f1.eval(1, y)
    f2H = f2.eval(1, y)

    I_H = f1H.ring._isolate_real_roots_pair(f1H, f2H, inf=u, sup=s, strict=True, basis=True)

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
                a, b = h.ring._refine_real_root(h, a, b, disjoint=y)

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
                a, b = h.ring._refine_real_root(h, a, b, disjoint=y)

                if b <= y:
                    I_L4_B.append(((b, a), indices, h))
                if a >= y:
                    I_L4_U.append(((b, a), indices, h))

    Q_L1_B = Q_L1
    Q_L2_B = _intervals_to_quadrants(I_L2_B, f1L2F, f2L2F, v, y)
    Q_L3_B = _intervals_to_quadrants(I_L3_B, f1H, f2H, s, u)
    Q_L4_B = _intervals_to_quadrants(I_L4_B, f1L4F, f2L4F, y, v)

    Q_L1_U = _intervals_to_quadrants(I_L1_U, f1H, f2H, u, s)
    Q_L2_U = _intervals_to_quadrants(I_L2_U, f1L2F, f2L2F, y, t)
    Q_L3_U = Q_L3
    Q_L4_U = _intervals_to_quadrants(I_L4_U, f1L4F, f2L4F, t, y)

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
    """Find a rectangle of minimum area for bisection."""
    min_area, j = None, None

    for i, (_, (u, v), (s, t), _, _, _, _) in enumerate(rectangles):
        area = (s - u)*(t - v)

        if min_area is None or area < min_area:
            min_area, j = area, i

    return rectangles.pop(j)


def _rectangle_small_p(a, b, eps):
    """Return ``True`` if the given rectangle is small enough."""
    (u, v), (s, t) = a, b

    if eps is not None:
        return s - u < eps and t - v < eps
    return True


class RealInterval:
    """A fully qualified representation of a real isolation interval."""

    def __init__(self, data, f):
        """Initialize new real interval with complete information."""
        if len(data) == 2:
            s, t = data

            self.neg = False

            ring = f.ring
            domain = ring.domain
            x = ring.gens[0]

            if domain.is_ComplexAlgebraicField and not domain.is_RealAlgebraicField:
                domain = domain.domain

            if s < 0:
                if t <= 0:
                    f, s, t, self.neg = f.compose(x, -x), -t, -s, True
                else:
                    raise ValueError(f"can't refine a real root in ({s}, {t})")

            a, b, c, d = _mobius_from_interval((s, t), domain.field)

            f = ring._transform(f, ring(a)*x + ring(b), ring(c)*x + ring(d))

            self.mobius = a, b, c, d
        else:
            self.mobius = data[:-1]
            self.neg = data[-1]

        self.f = f

    @property
    def a(self):
        """Return the position of the left end."""
        a, b, c, d = self.mobius

        if not self.neg:
            if a*d < b*c:
                return a/c
            return b/d
        if a*d > b*c:
            return -a/c
        return -b/d

    @property
    def b(self):
        """Return the position of the right end."""
        was = self.neg
        self.neg = not was
        rv = -self.a
        self.neg = was
        return rv

    @property
    def center(self):
        """Return the center of the real isolating interval."""
        return (self.a + self.b)/2

    def as_tuple(self):
        """Return tuple representation of real isolating interval."""
        return self.a, self.b

    def is_disjoint(self, other):
        """Return ``True`` if two isolation intervals are disjoint."""
        return self.b < other.a or other.b < self.a

    def refine(self):
        """Perform one step of real root refinement algorithm."""
        assert self.mobius is not None

        f = self.f
        ring = f.ring

        f, mobius = ring._inner_refine_real_root(f, self.mobius, steps=1, mobius=True)

        return RealInterval(mobius + (self.neg,), f)


class ComplexInterval:
    """A fully qualified representation of a complex isolation interval.
    The printed form is shown as (x1, y1) x (x2, y2): the southwest x northeast
    coordinates of the interval's rectangle.

    """

    def __init__(self, a, b, I, Q, F1, F2, f1, f2, conj=False):
        """Initialize new complex interval with complete information."""
        self.a, self.b = a, b  # the southwest and northeast corner: (x1, y1), (x2, y2)
        self.I, self.Q = I, Q

        self.f1, self.F1 = f1, F1
        self.f2, self.F2 = f2, F2

        self.domain = self.f1.ring.domain
        self.conj = conj

    @property
    def ax(self):
        """Return ``x`` coordinate of south-western corner."""
        return self.a[0]

    @property
    def ay(self):
        """Return ``y`` coordinate of south-western corner."""
        if not self.conj:
            return +self.a[1]
        return -self.b[1]

    @property
    def bx(self):
        """Return ``x`` coordinate of north-eastern corner."""
        return self.b[0]

    @property
    def by(self):
        """Return ``y`` coordinate of north-eastern corner."""
        if not self.conj:
            return +self.b[1]
        return -self.a[1]

    @property
    def center(self):
        """Return the center of the complex isolating interval."""
        return (self.ax + self.bx)/2, (self.ay + self.by)/2

    def as_tuple(self):
        """Return tuple representation of complex isolating interval."""
        return (self.ax, self.ay), (self.bx, self.by)

    def conjugate(self):
        """Return conjugated isolating interval."""
        return ComplexInterval(self.a, self.b, self.I, self.Q,
                               self.F1, self.F2, self.f1, self.f2,
                               not self.conj)

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

        dom = self.domain.unify(other.domain)
        ring = self.f1.ring.clone(domain=dom)
        rring = dom.poly_ring(*reversed(ring.symbols))
        resultants = []
        for i in (self, other):
            re, im = map(operator.methodcaller('set_ring', rring), (i.f1, i.f2))
            resultants.append(re.resultant(im))
        gcd = ring.drop(1).gcd(*resultants)
        gcd_roots = ring.drop(1)._isolate_real_roots(gcd,
                                                     inf=max(self.ax, other.ax),
                                                     sup=min(self.bx, other.bx))
        if len(gcd_roots) != 1:
            return False

        l, r = gcd_roots[0][0]
        if l == r:
            # Can't use _count_complex_roots() here, as the new isolation
            # interval will be degenerate: a vertical line segment.  Make
            # sure we only count roots on the northern/western edges and
            # on the north-western corner of the original isolation rectangle.
            for i in (self, other):
                dom = i.domain.algebraic_field(I)
                f1 = i.f1.eval(a=l).set_domain(dom)
                f2 = i.f2.eval(a=l).set_domain(dom)
                f = f1 + f2*dom.unit
                x = f.ring.gens[0]
                f = f.compose(0, -dom.unit*x)
                if i.conj:
                    f = f.ring.from_list([_.conjugate() for _ in f.all_coeffs()])
                f = f.compose(0, dom.unit*x)
                r = f.ring._isolate_real_roots(f, inf=i.ay, sup=i.by)
                if len([1 for _ in r if ((not i.conj and i.ay < _[0][0]) or
                                         (i.conj and _[0][1] < i.by)) and l < i.bx]) != 1:
                    return False
            return True
        return all(_get_rectangle(i.f1, i.f2, (l, i.ay), (r, i.by))[0] == 1
                   for i in (self, other))

    def refine(self, vertical=False):
        """Perform one step of complex root refinement algorithm."""
        a, b = (u, v), (s, t) = self.a, self.b

        I, Q = self.I, self.Q

        f1, F1 = self.f1, self.F1
        f2, F2 = self.f2, self.F2

        if s - u > t - v or vertical:
            D_L, D_R = _vertical_bisection(1, a, b, I, Q, F1, F2, f1, f2)

            if D_L[0] == 1:
                _, a, b, I, Q, F1, F2 = D_L
            else:
                _, a, b, I, Q, F1, F2 = D_R
        else:
            D_B, D_U = _horizontal_bisection(1, a, b, I, Q, F1, F2, f1, f2)

            if D_B[0] == 1:
                _, a, b, I, Q, F1, F2 = D_B
            else:
                _, a, b, I, Q, F1, F2 = D_U

        return ComplexInterval(a, b, I, Q, F1, F2, f1, f2, self.conj)


class _FindRoot:
    """Mixin class for computing polynomial roots."""

    def _sturm(self, f):
        domain = self.domain

        if not domain.is_Field:
            raise DomainError(f"can't compute Sturm sequence over {domain}")

        f = f.sqf_part()

        sturm = [f, f.diff()]

        while sturm[-1]:
            s = sturm[-2] % sturm[-1]
            sturm.append(-s)

        return sturm[:-1]

    def _sign_variations(self, f):
        """
        Compute the number of sign variations of ``f`` in ``K[x]``.

        Examples
        ========

        >>> R, x = ring('x', ZZ)

        >>> R._sign_variations(x**4 - x**2 - x + 1)
        2

        """
        domain = self.domain
        prev, k = domain.zero, 0

        for coeff in f.all_coeffs():
            if coeff*prev < 0:
                k += 1

            if coeff:
                prev = coeff

        return k

    def _reverse(self, f):
        """Compute ``x**n * f(1/x)``, i.e.: reverse ``f`` in ``K[x]``."""
        n = f.degree()
        return self.from_terms([((n - i,), c) for (i,), c in f.items()])

    def _root_upper_bound(self, f):
        """Compute the LMQ upper bound for the positive roots of `f`.

        LMQ (Local Max Quadratic) bound was developed by
        Akritas-Strzebonski-Vigklas :cite:`Alkiviadis2009bounds`.

        """
        domain = self.domain

        n, P = len(f.all_coeffs()), []
        t = n * [1]
        if f.LC < 0:
            f = -f
        f = self._reverse(f)

        def ilog2(a):
            return int(math.log(a, 2))

        for i in range(n):
            b = int(-f[(n - 1 - i,)])
            if b <= 0:
                continue

            a, QL = ilog2(b), []

            for j in range(i + 1, n):
                b = int(f[(n - 1 - j,)])

                if b <= 0:
                    continue

                q = t[j] + a - ilog2(b)
                QL.append([q // (j - i), j])

            if not QL:
                continue

            q = min(QL)

            t[q[1]] = t[q[1]] + 1

            P.append(q[0])

        if P:
            return domain(2)**int(max(P) + 1)

    def _step_refine_real_root(self, f, M):
        """One step of positive real root refinement algorithm."""
        domain = self.domain
        x = self.gens[0]

        a, b, c, d = M

        if a == b and c == d:
            return f, (a, b, c, d)

        A = self._root_upper_bound(self._reverse(f))

        if A is not None:
            A = 1/domain.convert(A)
        else:
            A = domain.zero

        if A > 16:
            f = f.compose(x, A*x)
            a, c, A = A*a, A*c, domain.one

        if A >= 1:
            f = f.compose(x, x + A)
            b, d = A*a + b, A*c + d

            assert f[1]

        f, g = f.compose(x, x + 1), f

        a1, b1, c1, d1 = a, a + b, c, c + d

        if not f.eval(x, 0):
            return f, (b1, b1, d1, d1)

        k = self._sign_variations(f)

        if k == 1:
            a, b, c, d = a1, b1, c1, d1
        else:
            f = self._reverse(g).compose(x, x + 1)

            assert f[1]

            a, b, c, d = b, a + b, d, c + d

        return f, (a, b, c, d)

    def _inner_refine_real_root(self, f, M, eps=None, steps=None, disjoint=None, mobius=False):
        """Refine a positive root of `f` given a Mobius transform or an interval."""
        a, b, c, d = M

        while not c:
            f, (a, b, c, d) = self._step_refine_real_root(f, (a, b, c, d))

        if eps is not None and steps is not None:
            for _ in range(steps):
                if abs(a/c - b/d) >= eps:
                    f, (a, b, c, d) = self._step_refine_real_root(f, (a, b, c, d))
                else:
                    break
        else:
            if eps is not None:
                while abs(a/c - b/d) >= eps:
                    f, (a, b, c, d) = self._step_refine_real_root(f, (a, b, c, d))

            if steps is not None:
                for _ in range(steps):
                    f, (a, b, c, d) = self._step_refine_real_root(f, (a, b, c, d))

        if disjoint is not None:
            while True:
                u, v = _mobius_to_interval((a, b, c, d))

                if u < disjoint < v:
                    f, (a, b, c, d) = self._step_refine_real_root(f, (a, b, c, d))
                else:
                    break

        return (f, (a, b, c, d)) if mobius else _mobius_to_interval((a, b, c, d))

    def _outer_refine_real_root(self, f, s, t, eps=None, steps=None, disjoint=None):
        """Refine a positive root of `f` given an interval `(s, t)`."""
        domain = self.domain
        x = self.gens[0]
        a, b, c, d = _mobius_from_interval((s, t), domain)

        f = self._transform(f, a*x + b, c*x + d)

        if self._sign_variations(f) != 1:
            raise RefinementFailedError(f'there should be exactly one root in ({s}, {t}) interval')

        return self._inner_refine_real_root(f, (a, b, c, d), eps=eps, steps=steps, disjoint=disjoint)

    def _refine_real_root(self, f, s, t, eps=None, steps=None, disjoint=None):
        """Refine real root's approximating interval to the given precision."""
        domain = self.domain.field
        new_ring = self.clone(domain=domain)
        x = new_ring.gens[0]
        f = f.set_domain(domain)
        f = f.clear_denoms()[1]

        if not domain.is_RationalField and not domain.is_RealAlgebraicField:
            raise DomainError(f'real root refinement not supported over {domain}')

        if s == t:
            return s, t

        if s > t:
            s, t = t, s

        negative = False

        if s < 0:
            if t <= 0:
                f, s, t, negative = f.compose(x, -x), -t, -s, True
            else:
                raise ValueError(f"can't refine a real root in ({s}, {t})")

        if negative and disjoint is not None:
            if disjoint < 0:
                disjoint = -disjoint
            else:
                disjoint = None

        s, t = new_ring._outer_refine_real_root(f, s, t, eps=eps, steps=steps, disjoint=disjoint)

        return (-t, -s) if negative else (s, t)

    def _inner_isolate_real_roots(self, f, eps=None):
        """Internal function for isolation positive roots up to given precision."""
        domain = self.domain
        x = self.gens[0]
        a, b, c, d = domain.one, domain.zero, domain.zero, domain.one
        k = self._sign_variations(f)

        roots, stack = [], [(a, b, c, d, f, k)]

        while stack:
            a, b, c, d, f, k = stack.pop()

            A = self._root_upper_bound(self._reverse(f))

            if A is not None:
                A = 1/domain.convert(A)
            else:
                A = domain.zero

            if A > 16:
                f = f.compose(x, A*x)
                a, c, A = A*a, A*c, domain.one

            if A >= 1:
                f = f.compose(x, x + A)
                b, d = A*a + b, A*c + d

                assert f[1]

                k = self._sign_variations(f)

                if k == 0:
                    continue
                if k == 1:
                    roots.append(self._inner_refine_real_root(f, (a, b, c, d),
                                                              eps=eps, mobius=True))
                    continue

            f1 = f.compose(x, x + 1)

            a1, b1, c1, d1, r = a, a + b, c, c + d, 0

            if not f1[1]:
                roots.append((f1, (b1, b1, d1, d1)))
                f1, r = f1 // x, 1

            k1 = self._sign_variations(f1)
            k2 = k - k1 - r

            a2, b2, c2, d2 = b, a + b, d, c + d

            if k2 > 1:
                f2 = self._reverse(f).compose(x, x + 1)

                if not f2[1]:
                    f2 //= x

                k2 = self._sign_variations(f2)
            else:
                f2 = None

            if k1 < k2:
                a1, a2, b1, b2 = a2, a1, b2, b1
                c1, c2, d1, d2 = c2, c1, d2, d1
                f1, f2, k1, k2 = f2, f1, k2, k1

            if not k1:
                continue

            if f1 is None:
                f1 = self._reverse(f).compose(x, x + 1)

                if not f1[1]:
                    f1 //= x

            if k1 == 1:
                roots.append(self._inner_refine_real_root(f1, (a1, b1, c1, d1),
                                                          eps=eps, mobius=True))
            else:
                stack.append((a1, b1, c1, d1, f1, k1))

            if not k2:
                continue

            if f2 is None:
                f2 = self._reverse(f).compose(x, x + 1)

                if not f2[1]:
                    f2 //= x

            if k2 == 1:
                roots.append(self._inner_refine_real_root(f2, (a2, b2, c2, d2),
                                                          eps=eps, mobius=True))
            else:
                stack.append((a2, b2, c2, d2, f2, k2))

        return roots

    def _discard_if_outside_interval(self, f, M, inf, sup, negative, mobius):
        """Discard an isolating interval if outside ``(inf, sup)``."""
        while True:
            u, v = _mobius_to_interval(M)

            if negative:
                u, v = -v, -u

            if (inf is None or u >= inf) and (sup is None or v <= sup):
                if not mobius:
                    return u, v
                return f, M
            if (sup is not None and u > sup) or (inf is not None and v < inf):
                return
            f, M = self._step_refine_real_root(f, M)

    def _inner_isolate_positive_roots(self, f, eps=None, inf=None, sup=None, mobius=False):
        """Iteratively compute disjoint positive root isolation intervals."""
        if sup is not None and sup < 0:
            return []

        roots = self._inner_isolate_real_roots(f, eps=eps)

        results = []

        if inf is not None or sup is not None:
            for f, M in roots:
                result = self._discard_if_outside_interval(f, M, inf, sup, False, mobius)

                if result is not None:
                    results.append(result)
        elif not mobius:
            for f, M in roots:
                u, v = _mobius_to_interval(M)
                results.append((u, v))
        else:
            results = roots

        return results

    def _inner_isolate_negative_roots(self, f, inf=None, sup=None, eps=None, mobius=False):
        """Iteratively compute disjoint negative root isolation intervals."""
        x = self.gens[0]

        if inf is not None and inf >= 0:
            return []

        roots = self._inner_isolate_real_roots(f.compose(x, -x), eps=eps)

        results = []

        if inf is not None or sup is not None:
            for f, M in roots:
                result = self._discard_if_outside_interval(f, M, inf, sup, True, mobius)

                if result is not None:
                    results.append(result)
        elif not mobius:
            for f, M in roots:
                u, v = _mobius_to_interval(M)
                results.append((-v, -u))
        else:
            results = roots

        return results

    def _count_real_roots(self, f, inf=None, sup=None):
        """Returns the number of distinct real roots of ``f`` in ``[inf, sup]``."""
        domain = self.domain.field
        new_ring = self.clone(domain=domain)

        if f.degree() <= 0:
            return 0

        f = f.set_domain(domain)

        if not domain.is_ComplexAlgebraicField and not domain.is_RationalField:
            raise DomainError(f"Can't count real roots in domain {domain}")

        if domain.is_ComplexAlgebraicField and not domain.is_RealAlgebraicField:
            return sum(k for *_, k in new_ring._isolate_real_roots(f, inf, sup))

        sturm = f.sturm()

        if inf is None:
            f_inf = new_ring.from_list([s.LC*(-1)**s.degree() for s in reversed(sturm)])
        else:
            f_inf = new_ring.from_list([s.eval(a=inf) for s in reversed(sturm)])
        signs_inf = new_ring._sign_variations(f_inf)

        if sup is None:
            f_sup = new_ring.from_list([s.LC for s in reversed(sturm)])
        else:
            f_sup = new_ring.from_list([s.eval(a=sup) for s in reversed(sturm)])
        signs_sup = new_ring._sign_variations(f_sup)

        count = abs(signs_inf - signs_sup)

        if inf is not None and not f.eval(a=inf):
            count += 1

        return count

    def _roots_bound(self, f):
        domain = self.domain
        lc = domain.to_expr(f.LC)
        B = 2*max(abs(domain.to_expr(c)/lc) for c in f.values())
        if not domain.is_AlgebraicField:
            return domain.convert(B)
        return domain.domain(int(100*B) + 1)/domain.domain(100)

    def _isolate_zero(self, f, inf, sup, sqf=False):
        """Handle special case of CF algorithm when ``f`` is homogeneous."""
        domain = self.domain
        (j,), f = f.terms_gcd()

        if j > 0:
            if (inf is None or inf <= 0) and (sup is None or 0 <= sup):
                if not sqf:
                    return [((domain.zero, domain.zero), j)], f
                return [(domain.zero, domain.zero)], f

        return [], f

    def _isolate_real_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        """Isolate real roots of a square-free polynomial."""
        domain = self.domain.field
        f = f.set_domain(domain)
        new_ring = self.clone(domain=domain)

        if not domain.is_ComplexAlgebraicField and not domain.is_RationalField:
            raise DomainError(f"Can't isolate real roots in domain {domain}")

        f = f.clear_denoms()[1]

        if domain.is_ComplexAlgebraicField and not domain.is_RealAlgebraicField:
            roots = [r for r, _ in new_ring._isolate_real_roots(f, eps=eps, inf=inf, sup=sup)]
            return [RealInterval((a, b), f) for (a, b) in roots] if blackbox else roots

        if f.degree() <= 0:
            return []

        I_zero, f = new_ring._isolate_zero(f, inf, sup, sqf=True)

        I_neg = new_ring._inner_isolate_negative_roots(f, eps=eps, inf=inf, sup=sup)
        I_pos = new_ring._inner_isolate_positive_roots(f, eps=eps, inf=inf, sup=sup)

        roots = sorted(I_neg + I_zero + I_pos)
        return [RealInterval((a, b), f) for (a, b) in roots] if blackbox else roots

    def _isolate_real_roots(self, f, eps=None, inf=None, sup=None):
        """Isolate real roots.

        Notes
        =====

        Implemented algorithms use Vincent-Akritas-Strzebonski (VAS) continued
        fractions approach :cite:`Alkiviadis2005comp`, :cite:`Alkiviadis2008cf`.

        """
        domain = self.domain.field
        f = f.set_domain(domain)
        new_ring = self.clone(domain=domain)

        if not domain.is_ComplexAlgebraicField and not domain.is_RationalField:
            raise DomainError(f'isolation of real roots not supported over {domain}')

        if domain.is_ComplexAlgebraicField and not domain.is_RealAlgebraicField:
            polys = [_.eval(1, 0) for _ in new_ring._real_imag(f)]
            roots = new_ring.to_ground()._isolate_real_roots_pair(*polys, eps=eps, inf=inf, sup=sup, strict=True)
            return [(_[0], _[1][0]) for _ in roots if _[1].keys() == {0, 1}]

        _, factors = f.sqf_list()
        factors = [(f.clear_denoms()[1], k) for f, k in factors]

        if len(factors) == 1:
            (f, k), = factors
            return [(r, k) for r in new_ring._isolate_real_roots_sqf(f, eps, inf, sup)]
        I_zero, f = new_ring._isolate_zero(f, inf, sup)
        I_neg, I_pos = new_ring._real_isolate_and_disjoin(factors, eps, inf, sup)
        return sorted(I_neg + I_zero + I_pos)

    def _real_isolate_and_disjoin(self, factors, eps=None, inf=None, sup=None, strict=False, basis=False):
        """Isolate real roots of a list of polynomials and disjoin intervals."""
        I_pos, I_neg = [], []

        for i, (f, k) in enumerate(factors):
            f = f.primitive()[1]
            for F, M in self._inner_isolate_positive_roots(f, eps=eps, inf=inf, sup=sup, mobius=True):
                I_pos.append((F, M, k, f))

            for G, N in self._inner_isolate_negative_roots(f, eps=eps, inf=inf, sup=sup, mobius=True):
                I_neg.append((G, N, k, f))

        for i, (f, M, k, F) in enumerate(I_pos):
            for j, (g, N, m, G) in enumerate(I_pos[i + 1:]):
                while not _disjoint_p(M, N, strict=strict):
                    f, M = self._inner_refine_real_root(f, M, steps=1, mobius=True)
                    g, N = self._inner_refine_real_root(g, N, steps=1, mobius=True)

                I_pos[i + j + 1] = g, N, m, G

            I_pos[i] = f, M, k, F

        for i, (f, M, k, F) in enumerate(I_neg):
            for j, (g, N, m, G) in enumerate(I_neg[i + 1:]):
                while not _disjoint_p(M, N, strict=strict):
                    f, M = self._inner_refine_real_root(f, M, steps=1, mobius=True)
                    g, N = self._inner_refine_real_root(g, N, steps=1, mobius=True)

                I_neg[i + j + 1] = g, N, m, G

            I_neg[i] = f, M, k, F

        if strict:
            for i, (f, M, k, F) in enumerate(I_neg):
                if not M[0]:
                    while not M[0]:
                        f, M = self._inner_refine_real_root(f, M, steps=1, mobius=True)

                    I_neg[i] = f, M, k, F
                    break

            for j, (g, N, m, G) in enumerate(I_pos):
                if not N[0]:
                    while not N[0]:
                        g, N = self._inner_refine_real_root(g, N, steps=1, mobius=True)

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

    def _isolate_real_roots_pair(self, f, g, eps=None, inf=None, sup=None, strict=False, basis=False):
        """Isolate real roots of a list of polynomials."""
        domain = self.domain.field
        new_ring = self.clone(domain=domain)

        if not domain.is_RationalField and not domain.is_RealAlgebraicField:
            raise DomainError(f'isolation of real roots not supported over {domain}')

        if (inf is None or inf <= 0) and (sup is None or 0 <= sup):
            zeros, zero_indices = True, {}
        else:
            zeros = False

        polys = [f, g]
        gcd = new_ring.zero

        for i, p in enumerate(polys):
            p = p.set_domain(domain)
            p = p.clear_denoms()[1]
            (j,), p = p.terms_gcd()
            polys[i] = p  # .set_domain(domain)

            if zeros and j > 0:
                zero_indices[i] = j

            gcd = new_ring.gcd(gcd, polys[i])

        polys = [p//gcd for p in polys]

        factors = collections.defaultdict(dict)

        for i, p in enumerate(polys):
            ni = (i + 1) % 2

            if not p and zeros and ni in zero_indices:
                zero_indices[i] = zero_indices[ni]

            for f, _ in new_ring.gcd(p, gcd).sqf_list()[1]:
                k1 = new_ring._trial_division(gcd, [f])[0][1]
                k2 = new_ring._trial_division(p, [f])[0][1]
                factors[f] = {i: k1 + k2, ni: k1}

                gcd //= f**k1
                p //= f**k2

            for f, k in p.sqf_list()[1]:
                factors[f] = {i: k}

        for f, k in gcd.sqf_list()[1]:
            factors[f] = {0: k, 1: k}

        I_neg, I_pos = new_ring._real_isolate_and_disjoin(tuple(factors.items()), eps=eps,
                                                          inf=inf, sup=sup, strict=strict,
                                                          basis=basis)
        I_zero = []

        if zeros and zero_indices:
            if not basis:
                I_zero = [((domain.zero, domain.zero), zero_indices)]
            else:
                I_zero = [((domain.zero, domain.zero), zero_indices, new_ring.gens[0])]

        return sorted(I_neg + I_zero + I_pos)

    def _isolate_all_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        """Isolate real and complex roots of a square-free polynomial ``f``."""
        return (self._isolate_real_roots_sqf(f, eps=eps, inf=inf, sup=sup, blackbox=blackbox),
                self._isolate_complex_roots_sqf(f, eps=eps, inf=inf, sup=sup, blackbox=blackbox))

    def _isolate_all_roots(self, f, eps=None, inf=None, sup=None):
        """Isolate real and complex roots of a non-square-free polynomial ``f``."""
        domain = self.domain

        if not domain.is_IntegerRing and not domain.is_RationalField:
            raise DomainError(f'isolation of real and complex roots is not supported over {domain}')

        _, factors = f.sqf_list()

        if len(factors) == 1:
            (f, k), = factors

            real_part, complex_part = self._isolate_all_roots_sqf(f, eps=eps, inf=inf, sup=sup)

            real_part = [((a, b), k) for (a, b) in real_part]
            complex_part = [((a, b), k) for (a, b) in complex_part]

            return real_part, complex_part
        raise NotImplementedError('only trivial square-free polynomials are supported')

    def _real_imag(self, f, _y=Dummy('y')):
        """
        Return bivariate polynomials ``f1`` and ``f2``, such that ``f = f1 + f2*I``.

        Examples
        ========

        >>> R, x = ring('x', ZZ)

        >>> R._real_imag(x**3 + x**2 + x + 1)
        (x**3 + x**2 - 3*x*_y**2 + x - _y**2 + 1, 3*x**2*_y + 2*x*_y - _y**3 + _y)

        >>> R, x = ring('x', QQ.algebraic_field(I))

        >>> R._real_imag(x**2 + I*x - 1)
        (x**2 - _y**2 - _y - 1, 2*x*_y + x)

        """
        domain = self.domain

        if domain.is_ComplexAlgebraicField and not domain.is_RealAlgebraicField:
            domain = domain.domain
        elif not domain.is_IntegerRing and not domain.is_RationalField and not domain.is_RealAlgebraicField:
            raise DomainError(f'computing real and imaginary parts is not supported over {domain}')

        new_ring = domain.poly_ring(Dummy('z'), self.symbols[0], _y)
        z, x, y = new_ring.gens

        f1 = f2 = new_ring.drop(0).zero

        if not f:
            return f1, f2

        t, h = new_ring.one, new_ring.zero
        g, d = x + y*z, 0

        for (i,), coeff in sorted(f.items(), key=lambda x: x[0]):
            t *= g**(i - d)
            d = i
            h += t*(coeff.real + z*coeff.imag)

        for (k,), h in h.eject(x, y).items():
            m = k % 4

            if not m:
                f1 += h
            elif m == 1:
                f2 += h
            elif m == 2:
                f1 -= h
            else:
                f2 -= h

        return f1, f2

    def _transform(self, f, p, q):
        """
        Evaluate functional transformation ``q**n * f(p/q)`` in ``K[x]``.

        Examples
        ========

        >>> R, x = ring('x', ZZ)

        >>> R._transform(x**2 - 2*x + 1, x**2 + 1, x - 1)
        x**4 - 2*x**3 + 5*x**2 - 4*x + 4

        """
        if not f:
            return self.zero

        n = f.degree()
        h, Q = f.LC*self.one, [self.one]

        for _ in range(n):
            Q.append(Q[-1]*q)

        for c, q in zip(reversed(f.all_coeffs()[:-1]), Q[1:]):
            h *= p
            q *= c
            h += q

        return h

    def _count_complex_roots(self, f, inf=None, sup=None, exclude=None):
        """Count all roots in [u + v*I, s + t*I] rectangle using Collins-Krandick algorithm."""
        domain = self.domain.field
        new_ring = self.clone(domain=domain)
        f = f.set_ring(new_ring)

        if not domain.is_ComplexAlgebraicField and not domain.is_RationalField:
            raise DomainError(f"Can't count complex roots in domain {domain}")

        if not all(isinstance(_, tuple) for _ in (inf, sup)):
            B = new_ring._roots_bound(f)

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

        f1, f2 = new_ring._real_imag(f)

        return _get_rectangle(f1, f2, (u, v), (s, t), exclude)[0]

    def _isolate_complex_roots_sqf(self, f, eps=None, inf=None, sup=None, blackbox=False):
        """Isolate complex roots of a square-free polynomial using Collins-Krandick algorithm."""
        if f.degree() <= 0:
            return []

        domain = self.domain.field
        new_ring = self.clone(domain=domain)
        f = f.set_ring(new_ring)

        if not domain.is_ComplexAlgebraicField and not domain.is_RationalField:
            raise DomainError(f"Can't isolate complex roots in domain {domain}")

        if not all(isinstance(_, tuple) for _ in (inf, sup)):
            B = new_ring._roots_bound(f)

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
            raise ValueError('not a valid complex isolation rectangle')

        if v < 0 < t:
            roots = new_ring._isolate_complex_roots_sqf(f, eps=eps, inf=(u, 0),
                                                        sup=(s, t), blackbox=True)
            if domain.is_RationalField or domain.is_RealAlgebraicField:
                _roots = []
                for root in roots:
                    croot = root.conjugate()
                    if croot.ay >= v:
                        _roots.append(croot)
                    _roots.append(root)
                roots = _roots
            else:
                # Take conjugated polynomial to get solutions in the
                # bottom half-plane.
                f = new_ring.from_list([_.conjugate() for _ in f.all_coeffs()])
                roots += [_.conjugate()
                          for _ in new_ring._isolate_complex_roots_sqf(f, eps=eps,
                                                                       inf=(u, 0), sup=(s, -v),
                                                                       blackbox=True)]
                roots = sorted(roots, key=lambda r: (r.ax, r.ay))

            return roots if blackbox else [r.as_tuple() for r in roots]

        f1, f2 = new_ring._real_imag(f)

        N, a, b, I_L, Q_L, F1, F2 = _get_rectangle(f1, f2, (u, v), (s, t))

        if not N:
            return []

        rectangles, roots = [], []
        if N == 1 and (v > 0 or t < 0):
            roots.append(ComplexInterval(a, b, I_L, Q_L, F1, F2, f1, f2))
        else:
            rectangles.append((N, a, b, I_L, Q_L, F1, F2))

        while rectangles:
            N, (u, v), (s, t), I_L, Q_L, F1, F2 = _depth_first_select(rectangles)

            if s - u > t - v:
                D_L, D_R = _vertical_bisection(N, (u, v), (s, t), I_L, Q_L, F1, F2, f1, f2)

                N_L, a, b, I_L, Q_L, F1_L, F2_L = D_L
                N_R, c, d, I_R, Q_R, F1_R, F2_R = D_R

                if N_L >= 1:
                    if N_L == 1 and _rectangle_small_p(a, b, eps):
                        roots.append(ComplexInterval(a, b, I_L, Q_L, F1_L, F2_L, f1, f2))
                    else:
                        rectangles.append(D_L)

                if N_R >= 1:
                    if N_R == 1 and _rectangle_small_p(c, d, eps):
                        roots.append(ComplexInterval(c, d, I_R, Q_R, F1_R, F2_R, f1, f2))
                    else:
                        rectangles.append(D_R)
            else:
                D_B, D_U = _horizontal_bisection(N, (u, v), (s, t), I_L, Q_L, F1, F2, f1, f2)

                N_B, a, b, I_B, Q_B, F1_B, F2_B = D_B
                N_U, c, d, I_U, Q_U, F1_U, F2_U = D_U

                if N_B >= 1:
                    if N_B == 1 and _rectangle_small_p(a, b, eps):
                        roots.append(ComplexInterval(a, b, I_B, Q_B, F1_B, F2_B, f1, f2))
                    else:
                        rectangles.append(D_B)

                if N_U >= 1:
                    if N_U == 1 and _rectangle_small_p(c, d, eps):
                        roots.append(ComplexInterval(c, d, I_U, Q_U, F1_U, F2_U, f1, f2))
                    else:
                        rectangles.append(D_U)

        roots = sorted(roots, key=lambda r: (r.ax, r.ay))
        return roots if blackbox else [r.as_tuple() for r in roots]
