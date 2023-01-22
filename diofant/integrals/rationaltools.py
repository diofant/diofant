"""This module implements tools for integrating rational functions."""

from ..core import Dummy, I, Integer, Lambda, Symbol, symbols, sympify
from ..domains import ZZ
from ..functions import atan, log
from ..polys import Poly, RootSum, cancel, resultant, roots
from ..simplify import collect
from ..solvers import solve


def ratint(f, x, **flags):
    """Performs indefinite integration of rational functions.

    Given a field `K` and a rational function `f = p/q`,
    where `p` and `q` are polynomials in `K[x]`,
    returns a function `g` such that `f = g'`.

    >>> ratint(36/(x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2), x)
    (12*x + 6)/(x**2 - 1) + 4*log(x - 2) - 4*log(x + 1)

    References
    ==========

    * :cite:`Bronstein2005integration`, pp. 35-70

    See Also
    ========

    diofant.integrals.integrals.Integral.doit
    ratint_logpart
    ratint_ratpart

    """
    if type(f) is not tuple:
        p, q = f.as_numer_denom()
    else:
        p, q = f

    p, q = p.as_poly(x, composite=False, field=True), q.as_poly(x, composite=False, field=True)

    coeff, p, q = p.cancel(q)
    poly, p = p.div(q)

    result = poly.integrate(x).as_expr()

    if p.is_zero:
        return coeff*result

    g, h = ratint_ratpart(p, q, x)

    P, Q = h.as_numer_denom()

    P = P.as_poly(x)
    Q = Q.as_poly(x)

    q, r = P.div(Q)

    result += g + q.integrate(x).as_expr()

    if not r.is_zero:
        symbol = flags.get('symbol', 't')

        if not isinstance(symbol, Symbol):
            t = Dummy(symbol)
        else:
            t = symbol.as_dummy()

        L = ratint_logpart(r, Q, x, t)

        ereal = flags.get('extended_real')

        if ereal is None:
            if type(f) is not tuple:
                atoms = f.atoms()
            else:
                p, q = f

                atoms = p.atoms() | q.atoms()

            for elt in atoms - {x}:
                if not elt.is_extended_real:
                    ereal = False
                    break
            else:
                ereal = True

        eps = Integer(0)

        if not ereal:
            for h, q in L:
                _, h = h.primitive()
                eps += RootSum(
                    q, Lambda(t, t*log(h.as_expr())), quadratic=True)
        else:
            for h, q in L:
                _, h = h.primitive()
                R = log_to_real(h, q, x, t)

                if R is not None:
                    eps += R
                else:
                    eps += RootSum(
                        q, Lambda(t, t*log(h.as_expr())), quadratic=True)

        result += eps

    return coeff*result


def ratint_ratpart(f, g, x):
    """Horowitz-Ostrogradsky algorithm.

    Given a field K and polynomials f and g in K[x], such that f and g
    are coprime and deg(f) < deg(g), returns fractions A and B in K(x),
    such that f/g = A' + B and B has square-free denominator.

    Examples
    ========

    >>> ratint_ratpart(1, x + 1, x)
    (0, 1/(x + 1))
    >>> ratint_ratpart(1, x**2 + y**2, x)
    (0, 1/(x**2 + y**2))
    >>> ratint_ratpart(36, x**5 - 2*x**4 - 2*x**3 + 4*x**2 + x - 2, x)
    ((12*x + 6)/(x**2 - 1), 12/(x**2 - x - 2))

    See Also
    ========

    ratint
    ratint_logpart

    """
    f = sympify(f).as_poly(x)
    g = sympify(g).as_poly(x)

    u, v, _ = g.cofactors(g.diff(x))

    n = u.degree()
    m = v.degree()

    A_coeffs = [Dummy('a' + str(n - i)) for i in range(n)]
    B_coeffs = [Dummy('b' + str(m - i)) for i in range(m)]

    C_coeffs = A_coeffs + B_coeffs

    A = Poly(A_coeffs, x, domain=ZZ.inject(*C_coeffs))
    B = Poly(B_coeffs, x, domain=ZZ.inject(*C_coeffs))

    H = f - A.diff(x)*v + A*(u.diff(x)*v).quo(u) - B*u

    result = solve(H.coeffs(), C_coeffs)[0]

    A = A.as_expr().subs(result)
    B = B.as_expr().subs(result)

    rat_part = cancel(A/u.as_expr(), x)
    log_part = cancel(B/v.as_expr(), x)

    return rat_part, log_part


def ratint_logpart(f, g, x, t=None):
    r"""Lazard-Rioboo-Trager algorithm.

    Given a field K and polynomials f and g in K[x], such that f and g
    are coprime, deg(f) < deg(g) and g is square-free, returns a list
    of tuples (s_i, q_i) of polynomials, for i = 1..n, such that s_i
    in K[t, x] and q_i in K[t], and::

                           ___    ___
                 d  f   d  \  `   \  `
                 -- - = --  )      )   a log(s_i(a, x))
                 dx g   dx /__,   /__,
                          i=1..n a | q_i(a) = 0

    Examples
    ========

    >>> ratint_logpart(1, x**2 + x + 1, x)
    [(Poly(x + 3*_t/2 + 1/2, x, domain='QQ[_t]'),
      Poly(3*_t**2 + 1, _t, domain='ZZ'))]
    >>> ratint_logpart(12, x**2 - x - 2, x)
    [(Poly(x - 3*_t/8 - 1/2, x, domain='QQ[_t]'),
      Poly(_t**2 - 16, _t, domain='ZZ'))]

    See Also
    ========

    ratint
    ratint_ratpart

    """
    f, g = sympify(f).as_poly(x), sympify(g).as_poly(x)

    t = t or Dummy('t')
    a, b = g, f - g.diff(x)*t.as_poly(x)

    res, R = resultant(a, b, includePRS=True)
    res = res.as_poly(t, composite=False)

    assert res, f"BUG: resultant({a}, {b}) can't be zero"

    R_map, H = {}, []

    for r in R:
        R_map[r.degree()] = r

    def _include_sign(c, sqf):
        if c.is_negative:
            h, k = sqf[0]
            sqf[0] = h*c, k

    C, res_sqf = res.sqf_list()
    _include_sign(C, res_sqf)

    for q, i in res_sqf:
        _, q = q.primitive()

        if g.degree() == i:
            H.append((g, q))
        else:
            h = R_map[i]
            h_lc = h.LC().as_poly(t, field=True)

            c, h_lc_sqf = h_lc.sqf_list()
            _include_sign(c, h_lc_sqf)

            for a, j in h_lc_sqf:
                h = h.quo((a.gcd(q)**j).as_poly(x))

            inv, coeffs = h_lc.invert(q), [Integer(1)]

            for coeff in h.coeffs()[1:]:
                T = (inv*coeff).rem(q)
                coeffs.append(T.as_expr())

            h = Poly(dict(zip(h.monoms(), coeffs)), x)

            H.append((h, q))

    return H


def log_to_atan(f, g):
    """Convert complex logarithms to real arctangents.

    Given a real field K and polynomials f and g in K[x], with g != 0,
    returns a sum h of arctangents of polynomials in K[x], such that:

                   dh   d         f + I g
                   -- = -- I log( ------- )
                   dx   dx        f - I g

    Examples
    ========

    >>> log_to_atan(x.as_poly(), Integer(1).as_poly(x))
    2*atan(x)
    >>> log_to_atan((x + Rational(1, 2)).as_poly(x), (sqrt(3)/2).as_poly(x))
    2*atan(2*sqrt(3)*x/3 + sqrt(3)/3)

    See Also
    ========

    log_to_real

    """
    if f.degree() < g.degree():
        f, g = -g, f

    f = f.to_field()
    g = g.to_field()

    p, q = f.div(g)

    if q.is_zero:
        return 2*atan(p.as_expr())
    s, t, h = g.gcdex(-f)
    u = (f*s + g*t).quo(h)
    A = 2*atan(u.as_expr())

    return A + log_to_atan(s, t)


def log_to_real(h, q, x, t):
    r"""Convert complex logarithms to real functions.

    Given real field K and polynomials h in K[t,x] and q in K[t],
    returns real function f such that:
                          ___
                  df   d  \  `
                  -- = --  )  a log(h(a, x))
                  dx   dx /__,
                         a | q(a) = 0

    Examples
    ========

    >>> log_to_real((x + 3*y/2 + Rational(1, 2)).as_poly(x),
    ...             (3*y**2 + 1).as_poly(y), x, y)
    2*sqrt(3)*atan(2*sqrt(3)*x/3 + sqrt(3)/3)/3
    >>> log_to_real((x**2 - 1).as_poly(), (-2*y + 1).as_poly(y), x, y)
    log(x**2 - 1)/2

    See Also
    ========

    log_to_atan

    """
    u, v = symbols('u,v', cls=Dummy)

    H = h.as_expr().subs({t: u + I*v}).expand()
    Q = q.as_expr().subs({t: u + I*v}).expand()

    H_map = collect(H, I, evaluate=False)
    Q_map = collect(Q, I, evaluate=False)

    a, b = H_map.get(Integer(1), Integer(0)), H_map.get(I, Integer(0))
    c, d = Q_map.get(Integer(1), Integer(0)), Q_map.get(I, Integer(0))

    R = resultant(c, d, v).as_poly(u)

    R_u_all = roots(R)
    R_q_all = roots(q)

    if sum(R_u_all.values()) < R.degree() or sum(R_q_all.values()) < q.degree():
        return

    R_u = {k: v for k, v in R_u_all.items() if k.is_extended_real}
    R_q = {k: v for k, v in R_q_all.items() if k.is_extended_real}

    result = Integer(0)

    for r_u in R_u:
        C = c.subs({u: r_u}).as_poly(v, extension=False)

        R_v_all = roots(C)
        if sum(R_v_all.values()) < C.degree():
            return
        R_v = {k: v for k, v in R_v_all.items() if k.is_extended_real is not False}

        R_v_paired = []  # take one from each pair of conjugate roots
        for r_v in R_v:
            if all(_ not in R_v_paired for _ in [+r_v, -r_v]):
                if r_v.could_extract_minus_sign():
                    R_v_paired.append(-r_v)

        for r_v in R_v_paired:
            D = d.subs({u: r_u, v: r_v})

            if D.cancel().evalf(2, chop=True) != 0:
                continue

            A = a.subs({u: r_u, v: r_v}).as_poly(x, extension=False)
            B = b.subs({u: r_u, v: r_v}).as_poly(x, extension=False)

            AB = (A**2 + B**2).as_expr()

            result += r_u*log(AB) + r_v*log_to_atan(A, B)

    for r in R_q:
        result += r*log(h.as_expr().subs({t: r}))

    return result
