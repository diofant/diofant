from ..core import Integer, Mul
from ..core.sympify import sympify


def residue(expr, x, x0):
    """
    Finds the residue of ``expr`` at the point ``x=x0``.

    The residue is defined as the coefficient of `1/(x - x_0)`
    in the power series expansion around `x=x_0`.

    This notion is essential for the Residue Theorem.

    Examples
    ========

    >>> residue(1/x, x, 0)
    1
    >>> residue(1/x**2, x, 0)
    0
    >>> residue(2/sin(x), x, 0)
    2

    References
    ==========

    * https://en.wikipedia.org/wiki/Residue_%28complex_analysis%29
    * https://en.wikipedia.org/wiki/Residue_theorem

    """
    # The current implementation uses series expansion to
    # calculate it. A more general implementation is explained in
    # the section 5.6 of the Bronstein's book {M. Bronstein:
    # Symbolic Integration I, Springer Verlag (2005)}. For purely
    # rational functions, the algorithm is much easier. See
    # sections 2.4, 2.5, and 2.7 (this section actually gives an
    # algorithm for computing any Laurent series coefficient for
    # a rational function). The theory in section 2.4 will help to
    # understand why the resultant works in the general algorithm.
    # For the definition of a resultant, see section 1.4 (and any
    # previous sections for more review).

    from ..simplify import collect
    from .order import Order

    expr = sympify(expr)
    if x0 != 0:
        expr = expr.subs({x: x + x0})
    s, n = Order(1, x), 1
    while s.has(Order) and s.getn() <= 0:
        s = expr.nseries(x, n=n)
        n *= 2
    s = collect(s.removeO(), x)
    if s.is_Add:
        args = s.args
    else:
        args = [s]
    res = Integer(0)
    for arg in args:
        c, m = arg.as_coeff_mul(x)
        m = Mul(*m)
        if not (m == 1 or m == x or (m.is_Pow and m.exp.is_Integer)):
            raise NotImplementedError(f'term of unexpected form: {m}')
        if m == 1/x:
            res += c
    return res
