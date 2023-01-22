"""Tools for manipulation of rational expressions."""

from ..core import Add, gcd_terms
from ..core.sympify import sympify


def together(expr, deep=False):
    """
    Denest and combine rational expressions using symbolic methods.

    This function takes an expression or a container of expressions
    and puts it (them) together by denesting and combining rational
    subexpressions. No heroic measures are taken to minimize degree
    of the resulting numerator and denominator. To obtain completely
    reduced expression use :func:`~diofant.polys.polytools.cancel`. However, :func:`together`
    can preserve as much as possible of the structure of the input
    expression in the output (no expansion is performed).

    A wide variety of objects can be put together including lists,
    tuples, sets, relational objects, integrals and others. It is
    also possible to transform interior of function applications,
    by setting ``deep`` flag to ``True``.

    By definition, :func:`together` is a complement to :func:`~diofant.polys.partfrac.apart`,
    so ``apart(together(expr))`` should return expr unchanged. Note
    however, that :func:`together` uses only symbolic methods, so
    it might be necessary to use :func:`~diofant.polys.polytools.cancel` to perform algebraic
    simplification and minimise degree of the numerator and denominator.

    Examples
    ========

    >>> together(1/x + 1/y)
    (x + y)/(x*y)
    >>> together(1/x + 1/y + 1/z)
    (x*y + x*z + y*z)/(x*y*z)

    >>> together(1/(x*y) + 1/y**2)
    (x + y)/(x*y**2)

    >>> together(1/(1 + 1/x) + 1/(1 + 1/y))
    (x*(y + 1) + y*(x + 1))/((x + 1)*(y + 1))

    >>> together(exp(1/x + 1/y))
    E**(1/y + 1/x)
    >>> together(exp(1/x + 1/y), deep=True)
    E**((x + y)/(x*y))

    >>> together(1/exp(x) + 1/(x*exp(x)))
    E**(-x)*(x + 1)/x

    >>> together(1/exp(2*x) + 1/(x*exp(3*x)))
    E**(-3*x)*(E**x*x + 1)/x

    """
    def _together(expr):
        if expr.is_Atom or (expr.is_Function and not deep):
            return expr
        if expr.is_Add:
            return gcd_terms(list(map(_together, Add.make_args(expr))))
        if expr.is_Pow:
            base = _together(expr.base)

            if deep:
                exp = _together(expr.exp)
            else:
                exp = expr.exp

            return expr.__class__(base, exp)
        return expr.__class__(*[_together(arg) for arg in expr.args])

    return _together(sympify(expr))
