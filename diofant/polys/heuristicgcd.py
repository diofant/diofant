"""Heuristic polynomial GCD algorithm (HEUGCD). """

from ..ntheory.modular import symmetric_residue
from .polyconfig import query
from .polyerrors import HeuristicGCDFailed


def heugcd(f, g):
    """
    Heuristic polynomial GCD in ``Z[X]``.

    Given univariate polynomials ``f`` and ``g`` in ``Z[X]``, returns
    their GCD and cofactors, i.e. polynomials ``h``, ``cff`` and ``cfg``
    such that::

          h = gcd(f, g), cff = quo(f, h) and cfg = quo(g, h)

    The algorithm is purely heuristic which means it may fail to compute
    the GCD. This will be signaled by raising an exception. In this case
    you will need to switch to another GCD method.

    The algorithm computes the polynomial GCD by evaluating polynomials
    ``f`` and ``g`` at certain points and computing (fast) integer GCD
    of those evaluations. The polynomial GCD is recovered from the integer
    image by interpolation. The evaluation proces reduces f and g variable
    by variable into a large integer. The final step is to verify if the
    interpolated polynomial is the correct GCD. This gives cofactors of
    the input polynomials as a side effect.

    Examples
    ========

    >>> R, x, y = ring("x y", ZZ)

    >>> heugcd((x + y)**2, x*(x + y))
    (x + y, x + y, x)

    References
    ==========

    * :cite:`Liao1995heuristic`

    """
    assert f.ring == g.ring and f.ring.domain.is_IntegerRing

    ring = f.ring
    x0 = ring.gens[0]
    domain = ring.domain

    gcd, f, g = f.extract_ground(g)

    f_norm = f.max_norm()
    g_norm = g.max_norm()

    B = domain(2*min(f_norm, g_norm) + 29)

    x = max(min(B, 99*domain.sqrt(B)),
            2*min(f_norm // abs(f.LC),
                  g_norm // abs(g.LC)) + 4)

    cofactors = domain.cofactors if ring.is_univariate else heugcd

    for i in range(query('HEU_GCD_MAX')):
        ff = f.eval(x0, x)
        gg = g.eval(x0, x)

        if ff and gg:
            h, cff, cfg = cofactors(ff, gg)
            h = _gcd_interpolate(h, x, ring)
            h = h.primitive()[1]

            cff_, r = divmod(f, h)

            if not r:
                cfg_, r = divmod(g, h)

                if not r:
                    h *= gcd
                    return h, cff_, cfg_

            cff = _gcd_interpolate(cff, x, ring)

            h, r = divmod(f, cff)

            if not r:
                cfg_, r = divmod(g, h)

                if not r:
                    h *= gcd
                    return h, cff, cfg_

            cfg = _gcd_interpolate(cfg, x, ring)

            h, r = divmod(g, cfg)

            if not r:
                cff_, r = divmod(f, h)

                if not r:
                    h *= gcd
                    return h, cff_, cfg

        x = 73794*x * domain.sqrt(domain.sqrt(x)) // 27011

    raise HeuristicGCDFailed('no luck')


def _gcd_interpolate(h, x, ring):
    """Interpolate polynomial GCD from integer GCD."""
    f, i = ring.zero, 0

    # TODO: don't expose poly repr implementation details
    if ring.is_univariate:
        while h:
            g = h % x
            g = symmetric_residue(g, x)
            h = (h - g) // x

            # f += X**i*g
            if g:
                f[(i,)] = g
            i += 1
    else:
        while h:
            g = h.trunc_ground(x)
            h = (h - g).quo_ground(x)

            # f += X**i*g
            if g:
                for monom, coeff in g.items():
                    f[(i,) + monom] = coeff
            i += 1

    if f.LC < 0:
        return -f
    else:
        return f
