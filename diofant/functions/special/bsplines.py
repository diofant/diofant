from ...core import Integer
from ...core.sympify import sympify
from ...sets import Interval
from .. import Piecewise, piecewise_fold


def _add_splines(c, b1, d, b2):
    """Construct c*b1 + d*b2."""
    if b1 == 0 or c == 0:
        rv = piecewise_fold(d*b2)
    elif b2 == 0 or d == 0:
        rv = piecewise_fold(c*b1)
    else:
        new_args = []
        n_intervals = len(b1.args)
        assert n_intervals == len(b2.args)
        new_args.append((c*b1.args[0].expr, b1.args[0].cond))
        for i in range(1, n_intervals - 1):
            new_args.append((
                c*b1.args[i].expr + d*b2.args[i - 1].expr,
                b1.args[i].cond
            ))
        new_args.append((d*b2.args[-2].expr, b2.args[-2].cond))
        new_args.append(b2.args[-1])
        rv = Piecewise(*new_args)

    return rv.expand()


def bspline_basis(d, knots, n, x, close=True):
    """The `n`-th B-spline at `x` of degree `d` with knots.

    B-Splines are piecewise polynomials of degree `d`.  They are defined on
    a set of knots, which is a sequence of integers or floats.

    The 0th degree splines have a value of one on a single interval:

        >>> d = 0
        >>> knots = range(5)
        >>> bspline_basis(d, knots, 0, x)
        Piecewise((1, (x >= 0) & (x <= 1)), (0, true))

    For a given ``(d, knots)`` there are ``len(knots)-d-1`` B-splines defined, that
    are indexed by ``n`` (starting at 0).

    Here is an example of a cubic B-spline:

        >>> bspline_basis(3, range(5), 0, x)
        Piecewise((x**3/6, (x >= 0) & (x < 1)),
                  (-x**3/2 + 2*x**2 - 2*x + 2/3,
                   (x >= 1) & (x < 2)),
                  (x**3/2 - 4*x**2 + 10*x - 22/3,
                   (x >= 2) & (x < 3)),
                  (-x**3/6 + 2*x**2 - 8*x + 32/3,
                   (x >= 3) & (x <= 4)),
                  (0, true))

    By repeating knot points, you can introduce discontinuities in the
    B-splines and their derivatives:

        >>> d = 1
        >>> knots = [0, 0, 2, 3, 4]
        >>> bspline_basis(d, knots, 0, x)
        Piecewise((-x/2 + 1, (x >= 0) & (x <= 2)), (0, true))

    It is quite time consuming to construct and evaluate B-splines. If you
    need to evaluate a B-splines many times, it is best to lambdify them
    first:

        >>> d = 3
        >>> knots = range(10)
        >>> b0 = bspline_basis(d, knots, 0, x)
        >>> f = lambdify(x, b0)
        >>> y = f(0.5)

    See Also
    ========

    diofant.functions.special.bsplines.bspline_basis_set

    References
    ==========

    * https://en.wikipedia.org/wiki/B-spline

    """
    knots = [sympify(k) for k in knots]
    d = int(d)
    n = int(n)
    n_knots = len(knots)
    n_intervals = n_knots - 1
    if n + d + 1 > n_intervals:
        raise ValueError('n + d + 1 must not exceed len(knots) - 1')
    if d == 0:
        result = Piecewise(
            (1, Interval(knots[n], knots[n + 1], False,
                         not close).contains(x)),
            (0, True)
        )
    elif d > 0:
        denom = knots[n + d + 1] - knots[n + 1]
        if denom != 0:
            B = (knots[n + d + 1] - x)/denom
            b2 = bspline_basis(d - 1, knots, n + 1, x, close)
        else:
            b2 = B = Integer(0)

        denom = knots[n + d] - knots[n]
        if denom != 0:
            A = (x - knots[n])/denom
            b1 = bspline_basis(
                d - 1, knots, n, x, close and (B == 0 or b2 == 0))
        else:
            b1 = A = Integer(0)

        result = _add_splines(A, b1, B, b2)
    else:
        raise ValueError(f'degree must be non-negative: {n!r}')
    return result


def bspline_basis_set(d, knots, x):
    """Return the ``len(knots)-d-1`` B-splines at ``x`` of degree ``d`` with ``knots``.

    This function returns a list of Piecewise polynomials that are the
    ``len(knots)-d-1`` B-splines of degree ``d`` for the given knots. This function
    calls ``bspline_basis(d, knots, n, x)`` for different values of ``n``.

    Examples
    ========

    >>> d = 2
    >>> knots = range(5)
    >>> splines = bspline_basis_set(d, knots, x)
    >>> splines
    [Piecewise((x**2/2, (x >= 0) & (x < 1)),
               (-x**2 + 3*x - 3/2, (x >= 1) & (x < 2)),
               (x**2/2 - 3*x + 9/2, (x >= 2) & (x <= 3)),
               (0, true)),
     Piecewise((x**2/2 - x + 1/2, (x >= 1) & (x < 2)),
               (-x**2 + 5*x - 11/2, (x >= 2) & (x < 3)),
               (x**2/2 - 4*x + 8, (x >= 3) & (x <= 4)),
               (0, true))]

    See Also
    ========

    diofant.functions.special.bsplines.bspline_basis

    """
    n_splines = len(knots) - d - 1
    return [bspline_basis(d, knots, i, x) for i in range(n_splines)]
