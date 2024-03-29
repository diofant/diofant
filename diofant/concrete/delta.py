"""This module implements sums and products containing the Kronecker Delta function."""

from ..core import Add, Dummy, Integer, Mul, cacheit
from ..functions import KroneckerDelta, Piecewise, piecewise_fold
from ..logic import true
from ..polys import factor
from ..sets import Interval
from ..utilities import default_sort_key


@cacheit
def _expand_delta(expr, index):
    """Expand the first Add containing a simple KroneckerDelta."""
    if not expr.is_Mul:
        return expr
    delta = None
    func = Add
    terms = [Integer(1)]
    for h in expr.args:
        if delta is None and h.is_Add and _has_simple_delta(h, index):
            delta = True
            func = h.func
            terms = [terms[0]*t for t in h.args]
        else:
            terms = [t*h for t in terms]
    return func(*terms)


@cacheit
def _extract_delta(expr, index):
    """Extract a simple KroneckerDelta from the expression.

    Returns the tuple ``(delta, newexpr)`` where:

      - ``delta`` is a simple KroneckerDelta expression if one was found,
        or ``None`` if no simple KroneckerDelta expression was found.

      - ``newexpr`` is a Mul containing the remaining terms; ``expr`` is
        returned unchanged if no simple KroneckerDelta expression was found.

    Examples
    ========

    >>> from diofant.abc import i, j

    >>> _extract_delta(4*x*y*KroneckerDelta(i, j), i)
    (KroneckerDelta(i, j), 4*x*y)
    >>> _extract_delta(4*x*y*KroneckerDelta(i, j), k)
    (None, 4*x*y*KroneckerDelta(i, j))

    See Also
    ========

    diofant.functions.special.tensor_functions.KroneckerDelta
    deltaproduct
    deltasummation

    """
    if not _has_simple_delta(expr, index):
        return None, expr
    if isinstance(expr, KroneckerDelta):
        return expr, Integer(1)
    assert expr.is_Mul
    delta = None
    terms = []

    for arg in expr.args:
        if delta is None and _is_simple_delta(arg, index):
            delta = arg
        else:
            terms.append(arg)
    return delta, expr.func(*terms)


@cacheit
def _has_simple_delta(expr, index):
    """
    Returns True if ``expr`` is an expression that contains a KroneckerDelta
    that is simple in the index ``index``, meaning that this KroneckerDelta
    is nonzero for a single value of the index ``index``.

    """
    if expr.has(KroneckerDelta):
        if _is_simple_delta(expr, index):
            return True
        if expr.is_Add or expr.is_Mul:
            for arg in expr.args:
                if _has_simple_delta(arg, index):
                    return True
    return False


@cacheit
def _is_simple_delta(delta, index):
    """
    Returns True if ``delta`` is a KroneckerDelta and is nonzero for a single
    value of the index ``index``.

    """
    if isinstance(delta, KroneckerDelta) and delta.has(index):
        p = (delta.args[0] - delta.args[1]).as_poly(index)
        if p:
            return p.degree() == 1
    return False


@cacheit
def _remove_multiple_delta(expr):
    """Evaluate products of KroneckerDelta's."""
    from ..solvers import solve
    if expr.is_Add:
        return expr.func(*list(map(_remove_multiple_delta, expr.args)))
    if not expr.is_Mul:
        return expr
    eqs = []
    newargs = []
    for arg in expr.args:
        if isinstance(arg, KroneckerDelta):
            eqs.append(arg.args[0] - arg.args[1])
        else:
            newargs.append(arg)
    if not eqs:
        return expr
    solns = solve(eqs)
    if len(solns) == 0:
        return Integer(0)
    if len(solns) == 1:
        for key in solns[0]:
            newargs.append(KroneckerDelta(key, solns[0][key]))
        expr2 = expr.func(*newargs)
        if expr != expr2:
            return _remove_multiple_delta(expr2)
    return expr


@cacheit
def _simplify_delta(expr):
    """Rewrite a KroneckerDelta's indices in its simplest form."""
    from ..solvers import solve
    if isinstance(expr, KroneckerDelta):
        slns = solve(expr.args[0] - expr.args[1])
        if slns and len(slns) == 1:
            return Mul(*[KroneckerDelta(*(key, value))
                         for key, value in slns[0].items()])
    return expr


@cacheit
def deltaproduct(f, limit):
    """Handle products containing a KroneckerDelta.

    See Also
    ========

    deltasummation
    diofant.functions.special.tensor_functions.KroneckerDelta
    diofant.concrete.products.product

    """
    from .products import product

    if (limit[2] - limit[1] < 0) == true:
        return Integer(1)

    if not f.has(KroneckerDelta):
        return product(f, limit)

    if f.is_Add:
        # Identify the term in the Add that has a simple KroneckerDelta
        delta = None
        terms = []
        for arg in sorted(f.args, key=default_sort_key):
            if delta is None and _has_simple_delta(arg, limit[0]):
                delta = arg
            else:
                terms.append(arg)
        newexpr = f.func(*terms)
        result = deltaproduct(newexpr, limit)
        if isinstance(limit[1], int) and isinstance(limit[2], int):
            result += sum(deltaproduct(newexpr,
                                       (limit[0], limit[1], ik - 1)) *
                          delta.subs({limit[0]: ik}) *
                          deltaproduct(newexpr,
                                       (limit[0], ik + 1, limit[2]))
                          for ik in range(int(limit[1]), int(limit[2] + 1)))
        else:
            k = Dummy('kprime', integer=True)
            result += deltasummation(deltaproduct(newexpr,
                                                  (limit[0],
                                                   limit[1], k - 1)) *
                                     delta.subs({limit[0]: k}) *
                                     deltaproduct(newexpr, (limit[0],
                                                            k + 1, limit[2])),
                                     (k, limit[1], limit[2]),
                                     no_piecewise=_has_simple_delta(newexpr,
                                                                    limit[0]))
        return _remove_multiple_delta(result)

    delta, _ = _extract_delta(f, limit[0])

    if not delta:
        g = _expand_delta(f, limit[0])
        if f != g:
            return factor(deltaproduct(g, limit))
        return product(f, limit)

    return (_remove_multiple_delta(f.subs({limit[0]: limit[1]}) *
                                   KroneckerDelta(limit[2], limit[1])) +
            _simplify_delta(KroneckerDelta(limit[2], limit[1] - 1)))


@cacheit
def deltasummation(f, limit, no_piecewise=False):
    """Handle summations containing a KroneckerDelta.

    The idea for summation is the following:

    - If we are dealing with a KroneckerDelta expression, i.e.
      KroneckerDelta(g(x), j), we try to simplify it.

      If we could simplify it, then we sum the resulting expression.
      We already know we can sum a simplified expression, because only
      simple KroneckerDelta expressions are involved.

      If we couldn't simplify it, there are two cases:

      1) The expression is a simple expression: we return the summation,
         taking care if we are dealing with a Derivative or with a proper
         KroneckerDelta.

      2) The expression is not simple (i.e. KroneckerDelta(cos(x))): we
         can do nothing at all.

    - If the expr is a multiplication expr having a KroneckerDelta term:

      First we expand it.

      If the expansion did work, then we try to sum the expansion.

      If not, we try to extract a simple KroneckerDelta term, then we
      have two cases:

      1) We have a simple KroneckerDelta term, so we return the summation.

      2) We didn't have a simple term, but we do have an expression with
         simplified KroneckerDelta terms, so we sum this expression.

    Examples
    ========

    >>> i, j = symbols('i, j', integer=True, finite=True)
    >>> deltasummation(KroneckerDelta(i, k), (k, -oo, oo))
    1
    >>> deltasummation(KroneckerDelta(i, k), (k, 0, oo))
    Piecewise((1, 0 <= i), (0, true))
    >>> deltasummation(KroneckerDelta(i, k), (k, 1, 3))
    Piecewise((1, (1 <= i) & (i <= 3)), (0, true))
    >>> deltasummation(k*KroneckerDelta(i, j)*KroneckerDelta(j, k),
    ...                (k, -oo, oo))
    j*KroneckerDelta(i, j)
    >>> deltasummation(j*KroneckerDelta(i, j), (j, -oo, oo))
    i
    >>> deltasummation(i*KroneckerDelta(i, j), (i, -oo, oo))
    j

    See Also
    ========

    deltaproduct
    diofant.functions.special.tensor_functions.KroneckerDelta
    diofant.concrete.sums.summation

    """
    from ..solvers import solve
    from .summations import summation

    if (limit[2] - limit[1] < 0) == true:
        return Integer(0)

    if not f.has(KroneckerDelta):
        return summation(f, limit)

    x = limit[0]

    g = _expand_delta(f, x)
    if g.is_Add:
        return piecewise_fold(g.func(*[deltasummation(h, limit,
                                                      no_piecewise)
                                       for h in g.args]))

    # try to extract a simple KroneckerDelta term
    delta, expr = _extract_delta(g, x)

    if not delta:
        return summation(f, limit)
    solns = solve(delta.args[0] - delta.args[1], x)
    assert len(solns) == 1
    value = solns[0][x]
    if no_piecewise:
        return expr.subs({x: value})
    return Piecewise((expr.subs({x: value}),
                      Interval(*limit[1:3]).as_relational(value)),
                     (0, True))
