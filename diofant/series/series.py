from ..core.sympify import sympify


def series(expr, x=None, x0=0, n=6, dir='+'):
    """Series expansion of ``expr`` in ``x`` around point ``x0``.

    See Also
    ========

    diofant.core.expr.Expr.series

    """
    expr = sympify(expr)
    return expr.series(x, x0, n, dir)
